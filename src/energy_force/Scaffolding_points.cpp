/**
 * @file Scaffolding_points.cpp
 * @author Y Ying (yying7@jh.edu)
 * @brief This file defines function in Mesh class that deals with
 * harmonic bond energy and force between the membrane and the
 * scaffolding.
 * @date 2023-04-07
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "mesh/Mesh.hpp"
#include "io/io.hpp"
#include <set>
#include <limits>
#include <cctype>
#include <gsl/gsl_eigen.h>

using namespace std;

namespace
{
constexpr double IDEALIZED_COLINEARITY_WARNING_THRESHOLD = 1.0e-2;

struct ReferenceBoundState
{
    int mol1 = -1;
    int mol2 = -1;
    std::string iface1;
    std::string iface2;
};

struct ReferenceMolecule
{
    Matrix com = Matrix(3, 1);
    std::unordered_map<std::string, Matrix> interfaces;
    Matrix referenceVector = Matrix(3, 1);
    bool hasReferenceVector = false;
};

std::string trim_copy(const std::string &value)
{
    const std::string whitespace = " \t\r\n";
    const size_t begin = value.find_first_not_of(whitespace);
    if (begin == std::string::npos)
    {
        return "";
    }
    const size_t end = value.find_last_not_of(whitespace);
    return value.substr(begin, end - begin + 1);
}

double clamp_unit(double value)
{
    if (!std::isfinite(value))
    {
        return 0.0;
    }
    return std::max(-1.0, std::min(1.0, value));
}

double safe_angle_from_vectors(const Matrix &v1, const Matrix &v2)
{
    const double norm1 = v1.calculate_norm();
    const double norm2 = v2.calculate_norm();
    if (norm1 == 0.0 || norm2 == 0.0 || !std::isfinite(norm1) || !std::isfinite(norm2))
    {
        return 0.0;
    }
    const double cosine = clamp_unit(dot_col(v1, v2) / (norm1 * norm2));
    return std::acos(cosine);
}

Matrix safe_cross(const Matrix &v1, const Matrix &v2)
{
    Matrix out(3, 1);
    cross(v1, v2, out);
    return out;
}

double safe_dihedral(const Matrix &p1, const Matrix &p2, const Matrix &p3, const Matrix &p4)
{
    const Matrix b1 = p2 - p1;
    const Matrix b2 = p3 - p2;
    const Matrix b3 = p4 - p3;
    const Matrix n1 = safe_cross(b1, b2);
    const Matrix n2 = safe_cross(b2, b3);
    const double n1norm = n1.calculate_norm();
    const double n2norm = n2.calculate_norm();
    const double b2norm = b2.calculate_norm();
    if (n1norm == 0.0 || n2norm == 0.0 || b2norm == 0.0 ||
        !std::isfinite(n1norm) || !std::isfinite(n2norm) || !std::isfinite(b2norm))
    {
        return 0.0;
    }
    Matrix m1(3, 1);
    cross(n1, b2, m1);
    const double y = dot_col(m1, n2) / b2norm;
    const double x = dot_col(n1, n2);
    if (!std::isfinite(x) || !std::isfinite(y))
    {
        return 0.0;
    }
    return std::atan2(y, x);
}

double wrapped_angle_difference(double value, double reference)
{
    if (!std::isfinite(value) || !std::isfinite(reference))
    {
        return 0.0;
    }
    const double diff = value - reference;
    return std::atan2(std::sin(diff), std::cos(diff));
}

Matrix identity3()
{
    Matrix id(3, 3, true);
    id.set_identity();
    return id;
}

Matrix rotation_matrix_from_axis_angle(const Matrix &axis, double angle)
{
    Matrix axisUnit(axis);
    const double norm = axisUnit.calculate_norm();
    if (norm == 0.0 || !std::isfinite(norm) || !std::isfinite(angle) || std::abs(angle) < 1.0e-16)
    {
        return identity3();
    }
    axisUnit /= norm;
    const double x = axisUnit(0, 0);
    const double y = axisUnit(1, 0);
    const double z = axisUnit(2, 0);
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const double oneMinusC = 1.0 - c;

    Matrix rot(3, 3, true);
    rot.set(0, 0, c + x * x * oneMinusC);
    rot.set(0, 1, x * y * oneMinusC - z * s);
    rot.set(0, 2, x * z * oneMinusC + y * s);
    rot.set(1, 0, y * x * oneMinusC + z * s);
    rot.set(1, 1, c + y * y * oneMinusC);
    rot.set(1, 2, y * z * oneMinusC - x * s);
    rot.set(2, 0, z * x * oneMinusC - y * s);
    rot.set(2, 1, z * y * oneMinusC + x * s);
    rot.set(2, 2, c + z * z * oneMinusC);
    return rot;
}

Matrix rotation_aligning_vectors(const Matrix &fromVec, const Matrix &toVec)
{
    Matrix fromUnit(fromVec);
    Matrix toUnit(toVec);
    const double fromNorm = fromUnit.calculate_norm();
    const double toNorm = toUnit.calculate_norm();
    if (fromNorm == 0.0 || toNorm == 0.0)
    {
        return identity3();
    }
    fromUnit /= fromNorm;
    toUnit /= toNorm;
    const double c = clamp_unit(dot_col(fromUnit, toUnit));
    if (c > 1.0 - 1.0e-12)
    {
        return identity3();
    }
    if (c < -1.0 + 1.0e-12)
    {
        Matrix trialAxis(3, 1);
        trialAxis.set(0, 0, 1.0);
        Matrix axis = safe_cross(fromUnit, trialAxis);
        if (axis.calculate_norm() < 1.0e-8)
        {
            trialAxis.set(0, 0, 0.0);
            trialAxis.set(1, 0, 1.0);
            axis = safe_cross(fromUnit, trialAxis);
        }
        return rotation_matrix_from_axis_angle(axis, M_PI);
    }
    const Matrix axis = safe_cross(fromUnit, toUnit);
    return rotation_matrix_from_axis_angle(axis, std::acos(c));
}

bool parse_reaction_line(const std::string &line, std::string &iface1, std::string &iface2)
{
    const size_t left = line.find("A(");
    if (left == std::string::npos)
    {
        return false;
    }
    const size_t leftEnd = line.find(')', left);
    const size_t plus = line.find('+', leftEnd);
    const size_t right = line.find("A(", plus);
    if (leftEnd == std::string::npos || plus == std::string::npos || right == std::string::npos)
    {
        return false;
    }
    const size_t rightEnd = line.find(')', right);
    if (rightEnd == std::string::npos)
    {
        return false;
    }

    iface1 = line.substr(left + 2, leftEnd - left - 2);
    iface2 = line.substr(right + 2, rightEnd - right - 2);
    const size_t bang1 = iface1.find('!');
    const size_t bang2 = iface2.find('!');
    if (bang1 != std::string::npos)
    {
        iface1 = iface1.substr(0, bang1);
    }
    if (bang2 != std::string::npos)
    {
        iface2 = iface2.substr(0, bang2);
    }
    return true;
}

std::unordered_map<std::string, Mesh::GagReactionTarget> load_gag_reaction_targets(const std::string &filepath)
{
    std::unordered_map<std::string, Mesh::GagReactionTarget> targets;
    if (filepath.empty())
    {
        return targets;
    }

    std::ifstream infile(filepath);
    if (!infile.is_open())
    {
        throw std::runtime_error("Failed to open Gag reaction file: " + filepath);
    }

    Mesh::GagReactionTarget currentTarget;
    bool inReactionBlock = false;
    bool haveReaction = false;
    std::string line;
    while (std::getline(infile, line))
    {
        const std::string trimmed = trim_copy(line);
        if (trimmed.empty() || trimmed[0] == '#')
        {
            continue;
        }
        if (trimmed == "start reactions")
        {
            inReactionBlock = true;
            continue;
        }
        if (trimmed == "end reactions")
        {
            break;
        }
        if (!inReactionBlock)
        {
            continue;
        }

        std::string iface1;
        std::string iface2;
        if (trimmed.find("<->") != std::string::npos && parse_reaction_line(trimmed, iface1, iface2))
        {
            if (haveReaction && currentTarget.sigma > 0.0)
            {
                const std::string key = currentTarget.interface1 + "|" + currentTarget.interface2;
                targets[key] = currentTarget;
                targets[currentTarget.interface2 + "|" + currentTarget.interface1] = currentTarget;
            }
            currentTarget = Mesh::GagReactionTarget();
            currentTarget.interface1 = iface1;
            currentTarget.interface2 = iface2;
            haveReaction = true;
            continue;
        }
        if (!haveReaction)
        {
            continue;
        }
        if (trimmed.find("sigma") == 0 && trimmed.find('=') != std::string::npos)
        {
            currentTarget.sigma = std::stod(trim_copy(trimmed.substr(trimmed.find('=') + 1)));
            continue;
        }
        if (trimmed.find("assocAngles") == 0 && trimmed.find('[') != std::string::npos)
        {
            const size_t leftBracket = trimmed.find('[');
            const size_t rightBracket = trimmed.find(']', leftBracket);
            if (rightBracket == std::string::npos)
            {
                continue;
            }
            std::stringstream ss(trimmed.substr(leftBracket + 1, rightBracket - leftBracket - 1));
            std::string value;
            int idx = 0;
            while (std::getline(ss, value, ',') && idx < 5)
            {
                currentTarget.assocAngles[idx++] = std::stod(trim_copy(value));
            }
        }
    }

    if (haveReaction && currentTarget.sigma > 0.0)
    {
        const std::string key = currentTarget.interface1 + "|" + currentTarget.interface2;
        targets[key] = currentTarget;
        targets[currentTarget.interface2 + "|" + currentTarget.interface1] = currentTarget;
    }

    return targets;
}

Matrix current_interface_coord(const Mesh::GagSubunit &subunit,
                               const std::string &interfaceName,
                               const std::vector<Matrix> &coms)
{
    return coms[subunit.pointIndex] + subunit.rotation * subunit.localInterfaces.at(interfaceName);
}

Matrix current_reference_endpoint(const Mesh::GagSubunit &subunit,
                                  const std::vector<Matrix> &coms)
{
    return coms[subunit.pointIndex] + subunit.rotation * subunit.localReferenceVector;
}

bool parse_reference_state_file(const std::string &filepath,
                                std::unordered_map<int, ReferenceMolecule> &molecules,
                                std::vector<ReferenceBoundState> &boundStates)
{
    std::ifstream infile(filepath);
    if (!infile.is_open())
    {
        return false;
    }

    std::string line;
    bool inMolecules = false;
    bool inBoundStates = false;
    bool inInterfaces = false;
    int currentMol = -1;
    std::string currentInterface;

    while (std::getline(infile, line))
    {
        const std::string trimmed = trim_copy(line);
        if (trimmed == "\"molecules\": {")
        {
            inMolecules = true;
            inBoundStates = false;
            continue;
        }
        if (trimmed == "\"bound_states\": [")
        {
            inMolecules = false;
            inBoundStates = true;
            continue;
        }

        if (inMolecules)
        {
            if (inInterfaces && !currentInterface.empty() && (trimmed == "}," || trimmed == "}"))
            {
                currentInterface.clear();
                continue;
            }
            if (inInterfaces && currentInterface.empty() && (trimmed == "}," || trimmed == "}"))
            {
                inInterfaces = false;
                continue;
            }
            if (!inInterfaces && (trimmed == "}," || trimmed == "}"))
            {
                currentMol = -1;
                continue;
            }
            if (currentMol == -1 && trimmed.size() > 4 && trimmed.front() == '"' && trimmed.find("\": {") != std::string::npos)
            {
                const std::string key = trimmed.substr(1, trimmed.find('"', 1) - 1);
                const bool numericKey = !key.empty() &&
                                        std::all_of(key.begin(), key.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); });
                if (!numericKey)
                {
                    continue;
                }
                currentMol = std::stoi(key);
                continue;
            }
            if (currentMol != -1 && trimmed == "\"com\": [")
            {
                std::string xline;
                std::string yline;
                std::string zline;
                std::getline(infile, xline);
                std::getline(infile, yline);
                std::getline(infile, zline);
                molecules[currentMol].com.set(0, 0, std::stod(trim_copy(xline)));
                molecules[currentMol].com.set(1, 0, std::stod(trim_copy(yline)));
                molecules[currentMol].com.set(2, 0, std::stod(trim_copy(zline)));
                continue;
            }
            if (currentMol != -1 && trimmed == "\"interfaces\": {")
            {
                inInterfaces = true;
                currentInterface.clear();
                continue;
            }
            if (inInterfaces && currentInterface.empty() &&
                trimmed.size() > 4 && trimmed.front() == '"' && trimmed.find("\": {") != std::string::npos)
            {
                currentInterface = trimmed.substr(1, trimmed.find('"', 1) - 1);
                continue;
            }
            if (inInterfaces && !currentInterface.empty() && trimmed == "\"coord\": [")
            {
                std::string xline;
                std::string yline;
                std::string zline;
                std::getline(infile, xline);
                std::getline(infile, yline);
                std::getline(infile, zline);
                Matrix coord(3, 1);
                coord.set(0, 0, std::stod(trim_copy(xline)));
                coord.set(1, 0, std::stod(trim_copy(yline)));
                coord.set(2, 0, std::stod(trim_copy(zline)));
                molecules[currentMol].interfaces[currentInterface] = coord;
                continue;
            }
            if (inInterfaces && !currentInterface.empty() && trimmed == "\"sphere_normal\": [")
            {
                std::string xline;
                std::string yline;
                std::string zline;
                std::getline(infile, xline);
                std::getline(infile, yline);
                std::getline(infile, zline);
                if (!molecules[currentMol].hasReferenceVector)
                {
                    molecules[currentMol].referenceVector.set(0, 0, std::stod(trim_copy(xline)));
                    molecules[currentMol].referenceVector.set(1, 0, std::stod(trim_copy(yline)));
                    molecules[currentMol].referenceVector.set(2, 0, std::stod(trim_copy(zline)));
                    molecules[currentMol].hasReferenceVector = true;
                }
                continue;
            }
            if (inInterfaces && !currentInterface.empty() && trimmed == "\"tangent_direction\": [")
            {
                std::string xline;
                std::string yline;
                std::string zline;
                std::getline(infile, xline);
                std::getline(infile, yline);
                std::getline(infile, zline);
                if (!molecules[currentMol].hasReferenceVector)
                {
                    molecules[currentMol].referenceVector.set(0, 0, std::stod(trim_copy(xline)));
                    molecules[currentMol].referenceVector.set(1, 0, std::stod(trim_copy(yline)));
                    molecules[currentMol].referenceVector.set(2, 0, std::stod(trim_copy(zline)));
                    molecules[currentMol].hasReferenceVector = true;
                }
                continue;
            }
        }

        if (inBoundStates && trimmed == "{")
        {
            ReferenceBoundState boundState;
            std::string mol1Line;
            std::string mol1IdxLine;
            std::string mol1IfaceLine;
            std::string mol2Line;
            std::string mol2IdxLine;
            std::string mol2IfaceLine;
            std::getline(infile, mol1Line);      // "mol1": {
            std::getline(infile, mol1IdxLine);   // molecule_index
            std::getline(infile, mol1IfaceLine); // interface_name
            std::getline(infile, line);          // interface_index
            std::getline(infile, line);          // },
            std::getline(infile, mol2Line);      // "mol2": {
            std::getline(infile, mol2IdxLine);   // molecule_index
            std::getline(infile, mol2IfaceLine); // interface_name
            std::getline(infile, line);          // interface_index
            std::getline(infile, line);          // }

            const size_t mol1Pos = mol1IdxLine.find(':');
            const size_t mol2Pos = mol2IdxLine.find(':');
            const size_t iface1Quote = mol1IfaceLine.find('"', mol1IfaceLine.find(':'));
            const size_t iface1QuoteEnd = mol1IfaceLine.find('"', iface1Quote + 1);
            const size_t iface2Quote = mol2IfaceLine.find('"', mol2IfaceLine.find(':'));
            const size_t iface2QuoteEnd = mol2IfaceLine.find('"', iface2Quote + 1);

            if (mol1Pos == std::string::npos || mol2Pos == std::string::npos ||
                iface1Quote == std::string::npos || iface1QuoteEnd == std::string::npos ||
                iface2Quote == std::string::npos || iface2QuoteEnd == std::string::npos)
            {
                continue;
            }

            boundState.mol1 = std::stoi(trim_copy(mol1IdxLine.substr(mol1Pos + 1)));
            boundState.mol2 = std::stoi(trim_copy(mol2IdxLine.substr(mol2Pos + 1)));
            boundState.iface1 = mol1IfaceLine.substr(iface1Quote + 1, iface1QuoteEnd - iface1Quote - 1);
            boundState.iface2 = mol2IfaceLine.substr(iface2Quote + 1, iface2QuoteEnd - iface2Quote - 1);
            boundStates.push_back(boundState);
        }
    }

    return !molecules.empty() && !boundStates.empty();
}
} // namespace

// --- 1. Find closest rcap from predefined list ---
int Mesh::findClosestRcap(double membraneCapRadius) {
    const std::vector<int> rcap_array = {10, 20, 25, 30, 35, 40};

    auto closest = std::min_element(rcap_array.begin(), rcap_array.end(),
        [membraneCapRadius](int a, int b) {
            return std::abs(a - membraneCapRadius) < std::abs(b - membraneCapRadius);
        });

    return *closest;
}

// --- 2. Interpolation ---
double Mesh::interpolateHeight(double r, const std::vector<double>& r_vals, const std::vector<double>& h_vals) {
    if (r <= r_vals.front()) return h_vals.front();
    if (r >= r_vals.back()) return h_vals.back();

    auto it = std::lower_bound(r_vals.begin(), r_vals.end(), r);
    size_t idx = std::distance(r_vals.begin(), it) - 1;

    double r0 = r_vals[idx], r1 = r_vals[idx + 1];
    double h0 = h_vals[idx], h1 = h_vals[idx + 1];

    return h0 + (r - r0) * (h1 - h0) / (r1 - r0);
}

// --- 3. Load and use in case3 ---
void Mesh::loadMembraneShapeProfile(double membraneCapRadius,
                              std::vector<double>& r_vals,
                              std::vector<double>& h_vals)
{
    int closest_rcap = findClosestRcap(membraneCapRadius);

    std::ostringstream filename;
    filename << "membrane_shape_" << closest_rcap << ".csv";

    std::vector<std::vector<double>> membrane_data = read_data_from_csv<double>(filename.str());

    r_vals.clear();
    h_vals.clear();
    for (const auto& row : membrane_data) {
        if (row.size() >= 2) {
            r_vals.push_back(row[0]);
            h_vals.push_back(row[1]);
        }
    }
}

/**
 * @brief
 * This method takes in the vector of spline point and calculate the
 * average coordinates. Based on the difference between spline points
 * and mesh vertices, a difference vector is calculated and compared to
 * the target bond length. (Supposed only in Z direction). Afterwards,
 * all the mesh points are moved in the direction of the target difference
 * vector.
 * 
 * @note
 * Refer to docs to see a graphic representation of the geometric variables
 *
 * @param fixDir default to true; fix move vector to (0, 0, 50) if set to
 * true; otherwise move vector is determined based on average scaffolding points
 * @return \code{true} if success
 */
bool Mesh::move_vertices_based_on_scaffolding(bool fixDir)
{
    auto safe_sqrt = [](double value) {
        return std::sqrt(std::max(0.0, value));
    };
    auto raise_above_baseline = [](double originalZ, double candidateZ) {
        return std::max(originalZ, candidateZ);
    };

    // If fixed direction, all vertices will be moved upwards by bond length
    // (lbond defined in Param class)
    if (fixDir)
    {
        for (Vertex& vertex : vertices)
        { 
            vertex.coord.set(2, 0, vertex.coord(2, 0) + param.lbond);
        }
        return true;
    }
    else
    {   // If not fixed direction, the vertices will be moved based on the
        // input scaffolding
        std::cout << "[Mesh::move_vertices_based_on_scaffolding] Custom move vertices based on scaffolding."
            << std::endl;
        // Here we move the vertices so that the membrane cap wraps around
        // the gag cap perfectly. 
        bool useDefault = false; // For potential switch
        this->centerScaffoldingSphere = this->find_center_of_scaffolding_sphere(useDefault);

        double gagSphereRadius = param.scaffoldingSphereRaidus;  //R 
        double gagCapRadius = this->approximate_scaffolding_cap_radius(useDefault); //r
        if (!std::isfinite(gagSphereRadius) || gagSphereRadius <= 0.0)
        {
            gagSphereRadius = 50.0;
        }
        if (!std::isfinite(gagCapRadius) || gagCapRadius < 0.0)
        {
            gagCapRadius = 0.0;
        }
        if (gagCapRadius >= gagSphereRadius)
        {
            const double clampedCapRadius = gagSphereRadius * (1.0 - 1.0e-6);
            std::cout << "[Mesh::move_vertices_based_on_scaffolding] Warning: inferred cap radius "
                      << gagCapRadius << " is not smaller than sphere radius "
                      << gagSphereRadius << ". Clamping to " << clampedCapRadius << std::endl;
            gagCapRadius = clampedCapRadius;
        }
        double lbond = param.lbond; //l
        double membraneSphereRadius = gagSphereRadius + lbond; // R+l
        double membraneCapRadiusSquared =  gagCapRadius * gagCapRadius / gagSphereRadius / gagSphereRadius
                                * membraneSphereRadius * membraneSphereRadius;//r^2/R^2*(R+l)^2
        membraneCapRadiusSquared = std::min(membraneCapRadiusSquared,
                                            membraneSphereRadius * membraneSphereRadius);
        double membraneCapRadius = sqrt(membraneCapRadiusSquared);
        // h0 = sqrt(Rs^2 - Rc^2)
        double h0 = safe_sqrt(membraneSphereRadius * membraneSphereRadius - membraneCapRadiusSquared);
        std:: cout << "[Mesh::move_vertices_based_on_scaffolding] Membrane sphere radius approximated : "
            << membraneSphereRadius << std::endl;
        // iterate thru all vertices
        // if vertices is with cap region (x^2 + y^2 <= r^2 / R^2 * (R + l)^2)
        // then move the vertices up
        // z - z0 = sqrt((R + l)^2 - (x - x0)^2 - (y - y0)^2)
        // otherwise use default zShift defined as z - z0 according to the equation above
        // but instead xyQuadrance = r^2/R^2*(R+l)^2 = membraneCapRadiusSquared
        double zShiftRelaxStart = safe_sqrt(membraneSphereRadius * membraneSphereRadius - membraneCapRadiusSquared);

        // Different modes for geometry model
        // Mode 2 tensionless approximation performs the best empirically
        // the rest modes are kept here only for reference / backward
        // compatibility if people want to reproduce old results
        int mode = 2;
        switch (mode)
        {
            // Mode 0: move vertices on cap, ignore the transition region
            case 0:
                for (Vertex& vertex : vertices)
                {
                    const double originalZ = vertex.coord(2, 0);
                    double xShift = vertex.coord(0, 0) - this->centerScaffoldingSphere(0, 0); //x-x0
                    double yShift = vertex.coord(1, 0) - this->centerScaffoldingSphere(1, 0); //y-y0
                    double xyQuadrance = xShift * xShift + yShift * yShift; // (x - x0)^2 + (y - y0)^2
                    //examine if vertex is with in cap region
                    if (xyQuadrance <= membraneCapRadiusSquared)
                    {
                        // z - z0 = sqrt((R + l)^2 - (x - x0)^2 - (y - y0)^2)
                        double zShift = safe_sqrt(membraneSphereRadius * membraneSphereRadius - xyQuadrance);
                        const double candidateZ = this->centerScaffoldingSphere(2, 0) + zShift;
                        vertex.coord.set(2, 0, raise_above_baseline(originalZ, candidateZ));
                    } else {
                        const double candidateZ = this->centerScaffoldingSphere(2, 0) + zShiftRelaxStart;
                        vertex.coord.set(2, 0, raise_above_baseline(originalZ, candidateZ));
                    }
                }
                break;
            // Mode 1: use a naive inversed arc to approximate transition region
            case 1:
                for (Vertex& vertex : vertices)
                {
                    const double originalZ = vertex.coord(2, 0);
                    // end of relaxation : z = 2h0 - Rs
                    double zShiftRelaxEnd = 2 * h0 - membraneSphereRadius;
                    double xShift = vertex.coord(0, 0) - this->centerScaffoldingSphere(0, 0); //x-x0
                    double yShift = vertex.coord(1, 0) - this->centerScaffoldingSphere(1, 0); //y-y0
                    double xyQuadrance = xShift * xShift + yShift * yShift; // (x - x0)^2 + (y - y0)^2
                    //examine if vertex is with in cap region
                    if (xyQuadrance <= membraneCapRadiusSquared)
                    {
                        // z - z0 = sqrt((R + l)^2 - (x - x0)^2 - (y - y0)^2)
                        double zShift = safe_sqrt(membraneSphereRadius * membraneSphereRadius - xyQuadrance);
                        const double candidateZ = this->centerScaffoldingSphere(2, 0) + zShift;
                        vertex.coord.set(2, 0, raise_above_baseline(originalZ, candidateZ));
                    } else if (xyQuadrance <= membraneCapRadiusSquared * 4) {
                        double rMinusTwoRc = sqrt(xyQuadrance) - 2 * membraneCapRadius;
                        // z = 2h0 - sqrt(- (r - 2Rc) ^ 2 + Rs^2)
                        double zShift = 2 * h0 - safe_sqrt( - rMinusTwoRc * rMinusTwoRc
                        + membraneSphereRadius * membraneSphereRadius);
                        const double candidateZ = this->centerScaffoldingSphere(2, 0) + zShift;
                        vertex.coord.set(2, 0, raise_above_baseline(originalZ, candidateZ));
                    } else {
                        const double candidateZ = this->centerScaffoldingSphere(2, 0) + zShiftRelaxEnd;
                        vertex.coord.set(2, 0, raise_above_baseline(originalZ, candidateZ));
                    }
                }
                break;
            // Mode 2: use tensionless approximation with relaxation length = membrane cap radius
            // Eqn: h(r) = c1 * r^2 + c3 ln(r) + c4
            // c1 = r0/(2rR (2r0 + rR)) * tanTheta
            // c3 = - tanTheta * r0 (r0 + rR) ^ 2 / (rR(2r0 + rR))
            case 2:
                for (Vertex& vertex : vertices)
                {
                    const double originalZ = vertex.coord(2, 0);
                    // constants
                    double r0 = membraneCapRadius;
                    double rR = r0 + r0 * std::max(param.relaxLengthRatioApproximation, 1.0e-6);
                    double r0PlusrR = r0 + rR;
                    double twor0PlusrR = 2 * r0 + rR;
                    // tangent theta (sector angle)
                    const double tangentDenominator = safe_sqrt(gagSphereRadius * gagSphereRadius - gagCapRadius * gagCapRadius);
                    double tangentSectorAngle = 0.0;
                    if (tangentDenominator > 1.0e-12)
                    {
                        tangentSectorAngle = gagCapRadius / tangentDenominator;
                    }
                    // c1 = r0/(2rR (2r0 + rR)) * tanTheta
                    double c1 = r0 / (2 * rR * twor0PlusrR) * tangentSectorAngle;
                    // c3 = - t * r0 (r0 + rR) ^ 2 / (rR(2r0 + rR))
                    double c3 = - r0 * r0PlusrR * r0PlusrR / rR / twor0PlusrR * tangentSectorAngle;
                    // c4 = -  c1 r0^2 - c3 ln(r0)
                    double c4 = - c1 * r0 * r0 - c3 * log(std::max(r0, 1.0e-12));

                    // Define geometry shift from origin
                    double xShift = vertex.coord(0, 0) - this->centerScaffoldingSphere(0, 0); //x-x0
                    double yShift = vertex.coord(1, 0) - this->centerScaffoldingSphere(1, 0); //y-y0
                    double xyQuadrance = xShift * xShift + yShift * yShift; // (x - x0)^2 + (y - y0)^2

                    // Use shape function to calcualte end of relaxation
                    double zShiftRelaxEnd = c1 * rR * rR + c3 * log(std::max(rR, 1.0e-12)) + c4;
                    const double capEdgeZ = this->centerScaffoldingSphere(2, 0) + zShiftRelaxStart;
                    const double flatPlaneZ = param.scaffoldingZeroPlaneZ;

                    //examine if vertex is with in cap region
                    if (xyQuadrance <= membraneCapRadiusSquared)
                    {
                        // z - z0 = sqrt((R + l)^2 - (x - x0)^2 - (y - y0)^2)
                        double zShift = safe_sqrt(membraneSphereRadius * membraneSphereRadius - xyQuadrance);
                        const double candidateZ = this->centerScaffoldingSphere(2, 0) + zShift;
                        vertex.coord.set(2, 0, raise_above_baseline(originalZ, candidateZ));
                        // Change C0 of insertion
                        for (int iAdjFaces : vertex.adjacentFaces) {
                            faces[iAdjFaces].spontCurvature = param.insertCurv;
                        }
                    } else if (xyQuadrance <= rR * rR) {
                        // Use shape function
                        
                        double r = sqrt(xyQuadrance);
                        double zShift = c1 * r * r + c3 * log(std::max(r, 1.0e-12)) + c4;
                        double transitionWeight = 0.0;
                        if (std::abs(zShiftRelaxEnd) > 1.0e-12)
                        {
                            transitionWeight = (zShift - zShiftRelaxEnd) / -zShiftRelaxEnd;
                        }
                        else
                        {
                            const double radialWeight = (rR - r) / std::max(rR - r0, 1.0e-12);
                            transitionWeight = radialWeight * radialWeight * (3.0 - 2.0 * radialWeight);
                        }
                        transitionWeight = std::max(0.0, std::min(1.0, transitionWeight));
                        const double transitionZ = flatPlaneZ + transitionWeight * (capEdgeZ - flatPlaneZ);
                        // Rule out boundary point
                        if (vertex.isBoundary)
                        {
                            const double candidateZ = flatPlaneZ;
                            vertex.coord.set(2, 0, raise_above_baseline(originalZ, candidateZ));
                        }
                        else
                        {
                            vertex.coord.set(2, 0, raise_above_baseline(originalZ, transitionZ));
                        }
                    } else {
                        const double candidateZ = flatPlaneZ;
                        vertex.coord.set(2, 0, raise_above_baseline(originalZ, candidateZ));
                    }
                }
                break;
            case 3:
                // Load membrane shape data
                std::vector<double> r_vals, h_vals;
                this->loadMembraneShapeProfile(membraneCapRadius, r_vals, h_vals);

                for (Vertex& vertex : vertices)
                {
                    const double originalZ = vertex.coord(2, 0);
                    // constants
                    double r0 = membraneCapRadius;
                    double rR = 120.0;
                    double r0PlusrR = r0 + rR;
                    double twor0PlusrR = 2 * r0 + rR;
                    // tangent theta (sector angle)
                    const double tangentDenominator = safe_sqrt(gagSphereRadius * gagSphereRadius - gagCapRadius * gagCapRadius);
                    double tangentSectorAngle = 0.0;
                    if (tangentDenominator > 1.0e-12)
                    {
                        tangentSectorAngle = gagCapRadius / tangentDenominator;
                    }

                    // Define geometry shift from origin
                    double xShift = vertex.coord(0, 0) - this->centerScaffoldingSphere(0, 0); //x-x0
                    double yShift = vertex.coord(1, 0) - this->centerScaffoldingSphere(1, 0); //y-y0
                    double xyQuadrance = xShift * xShift + yShift * yShift; // (x - x0)^2 + (y - y0)^2

                    // Use shape function to calcualte end of relaxation
                    double zShiftRelaxEnd = this->interpolateHeight(rR, r_vals, h_vals);;

                    //examine if vertex is with in cap region
                    if (xyQuadrance <= membraneCapRadiusSquared)
                    {
                        // z - z0 = sqrt((R + l)^2 - (x - x0)^2 - (y - y0)^2)
                        double zShift = safe_sqrt(membraneSphereRadius * membraneSphereRadius - xyQuadrance);
                        const double candidateZ = this->centerScaffoldingSphere(2, 0) + zShift;
                        vertex.coord.set(2, 0, raise_above_baseline(originalZ, candidateZ));
                    } else if (xyQuadrance <= membraneCapRadiusSquared * 4) {
                        // Use shape function
                        
                        double r = sqrt(xyQuadrance);
                        double zShift = this->interpolateHeight(r, r_vals, h_vals);
                        // Rule out boundary point
                        if (vertex.isBoundary)
                        {
                            const double candidateZ = this->centerScaffoldingSphere(2, 0) + zShiftRelaxEnd + zShiftRelaxStart;
                            vertex.coord.set(2, 0, raise_above_baseline(originalZ, candidateZ));
                        }
                        else
                        {
                            const double candidateZ = this->centerScaffoldingSphere(2, 0) + zShift + zShiftRelaxStart;
                            vertex.coord.set(2, 0, raise_above_baseline(originalZ, candidateZ));
                        }
                    } else {
                        const double candidateZ = this->centerScaffoldingSphere(2, 0) + zShiftRelaxEnd + zShiftRelaxStart;
                        vertex.coord.set(2, 0, raise_above_baseline(originalZ, candidateZ));
                    }
                }

                // Continue with the rest of your shape model...
                // e.g., use zShift to compute the membrane height at a point (x, y)
                break;
        }
        
        return true;
    }
}

/**
 * @brief A function that finds the center of the scaffolding sphere.
 * Currently the function assumes that the spherical cap is oriented
 * in a way that the point with biggest z-coord corresponds to the
 * x,y of the spherical center.
 * 
 * @note this function does NOT find the radius of the sphere. It takes
 * in the gag sphere radius from the parameter structure in mesh (this->param).
 * 
 * @param use_default if set to true, the function with return the defautl
 * center of scaffolding sphere (0, 0, 0)
 * @return instance of Matrix(3, 1) that represents the center of scaffolding
 * sphere.
 */
Matrix& Mesh::find_center_of_scaffolding_sphere(bool use_default){
    // if enabled, the default center of scaffolding sphere (0, 0, 0) will be used
    if (use_default){
        // return a matrix of (0, 0, 0)
        this->centerScaffoldingSphere = mat_calloc(3, 1);
        return this->centerScaffoldingSphere;
    }

    if (param.scaffoldingPoints.empty())
    {
        this->centerScaffoldingSphere = mat_calloc(3, 1);
        return this->centerScaffoldingSphere;
    }

    // Implement a simple function that determines center x-y based on max(z)
    // Move the center x-y by radius to determine the center coord
    double max_zcoord = param.scaffoldingPoints[0](2, 0);
    int index_with_max_zcoord = 0;
    // Iterate thru all vertices and compare to find maximum
    for (int i = 0; i < param.scaffoldingPoints.size(); i++) {
        if (param.scaffoldingPoints[i](2, 0) > max_zcoord) {
            max_zcoord = param.scaffoldingPoints[i](2, 0);
            index_with_max_zcoord = i;
        }
    }
    // Generate a new copy of matrix representing the center of scaffolding sphere
    this->centerScaffoldingSphere = Matrix(param.scaffoldingPoints[index_with_max_zcoord]);
    this->centerScaffoldingSphere.set(2, 0, this->centerScaffoldingSphere(2, 0) - this->param.scaffoldingSphereRaidus);
    std::cout << "Center of scaffolding sphere" << this->centerScaffoldingSphere << std::endl;
    // Return the matrix
    return this->centerScaffoldingSphere;
}

/**
 * @brief Approximates the radius of the scaffolding cap.
 *
 * This function approximates the radius of the scaffolding cap based on the furthest point from the center of the scaffolding sphere in terms of x and y coordinates.
 * 
 * @param use_default Flag indicating whether to use the default center of scaffolding sphere.
 *                    If set to true, the default center (0,0,0) will be used.
 *                    If set to false, the center will be determined based on the maximum z-coordinate of the scaffolding points.
 * 
 * @return The radius of the scaffolding cap.
 */
double Mesh::approximate_scaffolding_cap_radius(bool use_default){
    // Get center of scaffolding sphere using the member function
    this->find_center_of_scaffolding_sphere(use_default);
    if (param.scaffoldingPoints.empty())
    {
        return 0.0;
    }
    // Iterate thru all gags to find the gag that is furthest away in terms of x,y coord
    // Use quadrance (distance ^ 2)
    double max_quad_2d = 0.0;
    int index_with_max_2d_quad = 0;
    for (int i = 0; i < param.scaffoldingPoints.size(); i++) {
        // Calculate quadrance in x and y direction
        double dx = param.scaffoldingPoints[i](0, 0) - this->centerScaffoldingSphere(0, 0);
        double dy = param.scaffoldingPoints[i](1, 0) - this->centerScaffoldingSphere(1, 0);
        double quad_2d = dx * dx + dy * dy;
        // Get maximum quadrance
        if (quad_2d> max_quad_2d) {
            index_with_max_2d_quad = i;
            max_quad_2d = quad_2d;
        }
    }
    // Cap radius is the assumed to be distance in x,y direction
    double capRadius = sqrt(max_quad_2d);
    return capRadius;
}

/**
 * @brief
 * Used in getClosestVertexIndex
 * Calculates the squared distance between a point represented by a 3x1 matrix
 * and a Vertex object. Returns the squared distance.
 */
double Mesh::get_squared_distance_sp_and_v(const Matrix &scaffoldingPoint, const Vertex &vertex)
{
    // get difference vector
    Matrix diffVec = vertex.coord - scaffoldingPoint;
    return std::pow(diffVec(0, 0), 2) + std::pow(diffVec(1, 0), 2) + std::pow(diffVec(2, 0), 2);
}

/**
 * @brief
 * Gets a vector of indexes of vertices that are closest to
 * the `param.scaffoldingPoints` vector provided. Then sets
 * `param.scaffoldingPoints_correspondingVertexIndex` to represent
 * the vertices bonded with each scaffolding point.
 * 
 * @note
 * This section might be accelerated with ckdtree
 */
void Mesh::set_scaffolding_vertices_correspondence()
{
    // initialize
    vector<int> scaffoldingPoints_correspondingVertexIndex(param.scaffoldingPoints.size(), 0);
    // loop over spline points and find closest vertices for each
    for (int i = 0; i < param.scaffoldingPoints.size(); i++)
    {
        // initialize with the first vertex
        double minDistance2 = get_squared_distance_sp_and_v(param.scaffoldingPoints[i], vertices[0]);
        double newDistance2 = minDistance2;
        int minDistanceVertexIndex = 0;
        // loop over vertices to search for min distance
        for (int j = 1; j < vertices.size(); j++)
        {
            newDistance2 = get_squared_distance_sp_and_v(param.scaffoldingPoints[i], vertices[j]);
            // overide min and index if new d < min
            if (newDistance2 < minDistance2)
            {
                minDistance2 = newDistance2;
                minDistanceVertexIndex = j;
            }
        }
        // set corresponding vertex index to minDistanceVertexIndex
        scaffoldingPoints_correspondingVertexIndex[i] = minDistanceVertexIndex;
    }
    if (param.VERBOSE_MODE)
    {
        for (int i = 0; i < scaffoldingPoints_correspondingVertexIndex.size(); i++)
        {
            for (int j = 0; j < i; j++)
            {
                if (scaffoldingPoints_correspondingVertexIndex[i] == scaffoldingPoints_correspondingVertexIndex[j])
                    std::cout << "Repeated Vertex Index: " << scaffoldingPoints_correspondingVertexIndex[i] << endl;
            }
        }
    }
    param.scaffoldingPoints_correspondingVertexIndex = scaffoldingPoints_correspondingVertexIndex;
}

/*
 * Used in getCloesestVertexIndex
 * Calculate the distance between two points denoted by (3) vector
 * and Vertex respectively.
 * Return the distance
 */
double getDistance(Matrix &scaffoldingPoint, Vertex &vertex)
{
    // get difference vector
    bool singleSide = true; // set r to be negative when z_diffVec is negative
    // (membrane is lower than gag)
    Matrix diffVec = vertex.coord - scaffoldingPoint;
    // std::cout << "Mz: " << vertex.Coord[2] << "; Gz: " << scaffoldingPoint[2] << endl;
    double distance = diffVec.calculate_norm();
    if (singleSide && diffVec(2, 0) < 0.0)
    {
        return (-distance);
    }
    else
        return distance;
}

/**
 *
 * @brief Calculates the energy and force due to the harmonic bond between scaffold points and membrane vertices.
 * @param doLocalSearch Flag indicating whether or not to perform a local search for each vertex (default=false).
 * @return double The total energy of the system due to the scaffold-membrane interactions.
 *
 * The function iterates over each scaffold point, calculates the distance between the point and its corresponding vertex,
 * and calculates the energy and force due to the harmonic bond between the two. If the energy flag is not set to include
 * harmonic bonding, the function returns 0.
 *
 * @param[in] doLocalSearch Flag to indicate whether a local search should be performed for each vertex (default=false).
 *
 *
 * @note The energy is calculated according to the formula: E = 0.5 * k * (r - l)^2, where k is the spring constant, r is the
 * distance between the scaffold point and vertex, and l is the resting length of the bond.
 * @note The force is calculated according to the formula: F = -k * (r - l) * (unit vector pointing from vertex to scaffold point).
 * @note If the distance between the scaffold point and vertex is negative, the force is multiplied by -1 in order to prevent
 * overlapping between the scaffold and membrane.
 * @note If the verbose flag is set, the function prints out information about each vertex's force due to the scaffold-membrane bond.
 */
double Mesh::calculate_scaffolding_energy_force(bool doLocalSearch)
{
    if (!param.isEnergyHarmonicBondIncluded)
    {
        return 0;
    }

    // initialize total energy to zero; reset total force
    double totalEnergy = 0.0;
    forceTotalOnScaffolding = mat_calloc(3, 1);
    forceOnScaffoldingPoints.assign(param.scaffoldingPoints.size(), Matrix(3, 1));

    // iterate over spline points
    for (int i = 0; i < param.scaffoldingPoints.size(); i++)
    {
        int index = param.scaffoldingPoints_correspondingVertexIndex[i];
        double distance = getDistance(param.scaffoldingPoints[i], vertices[index]);
        // @test get absolute distance
        // distance = abs(distance);

        // calculate energy = 0.5 * k * (r - l)^2
        totalEnergy += 0.5 * param.springConst * pow(distance - param.lbond, 2);

        // vertices[index].energy.EnergySpline = 0.5 * springConst * (distance - lbond) * (distance - lbond);
        // vertices[index].energy.EnergyTotal += vertices[index].energy.EnergySpline;
        // calculate force = - k * (r - l) * normalize(vertexvec- splinepointvec)
        Matrix r_unit(3, 1);
        get_unit_vector(vertices[index].coord - param.scaffoldingPoints[i], r_unit);

        // If distance < 0.0, this means that the scaffolding points are above the plane
        // If using standard harmonic bonding equation, it would push the membrane further up
        // However because we do not want any overlapping between the membrane and the scaffolding
        // This is set to negative of original harmonic bond force in this case.
        if (distance > 0.0)
        {
            vertices[index].force.forceHarmonicBond = -param.springConst * (distance - param.lbond) * r_unit;
        }
        else
        {
            //adapt negative distance using:
            vertices[index].force.forceHarmonicBond = param.springConst * (distance - param.lbond) * r_unit;
            //vertices[index].force.forceHarmonicBond = param.springConst * (distance + param.lbond) * r_unit;
        }
        // Add harmonic bond force to total nodal force
        vertices[index].force.forceTotal += vertices[index].force.forceHarmonicBond;
        // Force exerted by membrane on scaffolding = reverse of force exerted by scaffolding on membrane
        forceTotalOnScaffolding -= vertices[index].force.forceHarmonicBond;
        forceOnScaffoldingPoints[i] -= vertices[index].force.forceHarmonicBond;
    }

    if (param.isGagScaffoldingEnergyIncluded || param.isIdealizedProteinLatticeEnergyIncluded)
    {
        const double internalEnergy = calculate_gag_scaffolding_internal_energy();
        totalEnergy += internalEnergy;
        if (param.isGagScaffoldingEnergyIncluded)
        {
            param.energy.energyGagScaffolding = internalEnergy;
        }
        else
        {
            param.energy.energyIdealizedProteinLattice = internalEnergy;
        }
    }
    // Commandline print if enabled
    if (param.VERBOSE_MODE)
    {
    std::cout << "[Mesh::calculate_scaffolding_energy_force] Total force exerted on scaffolding = ("
    << forceTotalOnScaffolding(0, 0) << ", "
    << forceTotalOnScaffolding(1, 0) << ", "
    << forceTotalOnScaffolding(2, 0) << ") " << std::endl;
    }
    return totalEnergy;
}

bool Mesh::propagate_scaffolding()
{
    const bool useInternalLatticeMode = param.isGagScaffoldingEnergyIncluded ||
                                        param.isIdealizedProteinLatticeEnergyIncluded;
    double scaffolding_propagation_stepsize = useInternalLatticeMode
                                            ? param.gagPropagationStepSize
                                            : 1.0e-7;
    int max_scaffolding_propagation_nsteps = param.propagateScaffoldingNstep;
    for (int i = 0; i < max_scaffolding_propagation_nsteps; i++){
        // recalculate scaffolding energy and force
        this->calculate_scaffolding_energy_force(false);
        if (useInternalLatticeMode)
        {
            const std::vector<Matrix> gagTranslationForces = calculate_gag_subunit_translation_forces_fd();
            const std::vector<Matrix> gagRotationTorques = calculate_gag_subunit_rotation_torques_fd();
            Matrix averageDisplacement(3, 1);
            for (int j = 0; j < param.scaffoldingPoints.size(); j++)
            {
                Matrix pointForce = forceOnScaffoldingPoints[j];
                pointForce += gagTranslationForces[j];
                Matrix displacement = scaffolding_propagation_stepsize * pointForce;
                param.scaffoldingPoints[j] += displacement;
                averageDisplacement += displacement;

                Matrix torqueStep = scaffolding_propagation_stepsize * gagRotationTorques[j];
                const double angle = torqueStep.calculate_norm();
                if (j < static_cast<int>(gagSubunits.size()) && angle > 0.0)
                {
                    Matrix rotationDelta = rotation_matrix_from_axis_angle(torqueStep, angle);
                    gagSubunits[j].rotation = rotationDelta * gagSubunits[j].rotation;
                }
            }
            if (!param.scaffoldingPoints.empty())
            {
                averageDisplacement = averageDisplacement / static_cast<double>(param.scaffoldingPoints.size());
            }
            scaffoldingMovementVector += averageDisplacement;
        }
        else
        {
            // move scaffolding based on calculated force
            Matrix displacement_vector = scaffolding_propagation_stepsize * forceTotalOnScaffolding;
            // iterate thru scaffolding and apply displacement vector
            for (int j = 0; j < param.scaffoldingPoints.size(); j++)
            {
                param.scaffoldingPoints[j] += displacement_vector;
            }
            // record movement in scaffoldingMovementVector
            scaffoldingMovementVector += displacement_vector;
        }
    }
    std::cout << "[Mesh::propagate_scaffolding()] Total movement : ("
    << scaffoldingMovementVector(0, 0) << ", "
    << scaffoldingMovementVector(1, 0) << ", "
    << scaffoldingMovementVector(2, 0) << ") " << std::endl;
    return true;
}

bool Mesh::pre_relax_gag_scaffolding()
{
    if ((!param.isGagScaffoldingEnergyIncluded && !param.isIdealizedProteinLatticeEnergyIncluded) ||
        param.gagPreRelaxSteps <= 0)
    {
        return false;
    }

    initialize_gag_scaffolding_topology();

    const double preRelaxStepSize = param.gagPropagationStepSize;
    Matrix totalMovement(3, 1);
    for (int stepIndex = 0; stepIndex < param.gagPreRelaxSteps; ++stepIndex)
    {
        const std::vector<Matrix> gagTranslationForces = calculate_gag_subunit_translation_forces_fd();
        const std::vector<Matrix> gagRotationTorques = calculate_gag_subunit_rotation_torques_fd();
        Matrix averageDisplacement(3, 1);

        for (int subunitIndex = 0; subunitIndex < static_cast<int>(param.scaffoldingPoints.size()); ++subunitIndex)
        {
            const Matrix displacement = preRelaxStepSize * gagTranslationForces[subunitIndex];
            param.scaffoldingPoints[subunitIndex] += displacement;
            averageDisplacement += displacement;

            if (subunitIndex < static_cast<int>(gagSubunits.size()))
            {
                const Matrix torqueStep = preRelaxStepSize * gagRotationTorques[subunitIndex];
                const double angle = torqueStep.calculate_norm();
                if (angle > 0.0)
                {
                    const Matrix rotationDelta = rotation_matrix_from_axis_angle(torqueStep, angle);
                    gagSubunits[subunitIndex].rotation = rotationDelta * gagSubunits[subunitIndex].rotation;
                }
            }
        }

        if (!param.scaffoldingPoints.empty())
        {
            averageDisplacement /= static_cast<double>(param.scaffoldingPoints.size());
        }
        totalMovement += averageDisplacement;
    }

    std::cout << "[Mesh::pre_relax_gag_scaffolding()] Steps: " << param.gagPreRelaxSteps
              << ". eInternal= " << calculate_gag_scaffolding_internal_energy()
              << ". Total average movement: (" << totalMovement(0, 0) << ", "
              << totalMovement(1, 0) << ", " << totalMovement(2, 0) << ")" << std::endl;
    return true;
}

bool Mesh::orient_scaffolding_plane_to_membrane()
{
    if (param.scaffoldingPoints.size() < 3)
    {
        gagInitialAlignmentRotation = identity3();
        return false;
    }

    Matrix centroid(3, 1);
    for (const Matrix &point : param.scaffoldingPoints)
    {
        centroid += point;
    }
    centroid /= static_cast<double>(param.scaffoldingPoints.size());

    gsl_matrix *cov = gsl_matrix_calloc(3, 3);
    for (const Matrix &point : param.scaffoldingPoints)
    {
        Matrix centered = point - centroid;
        for (int r = 0; r < 3; ++r)
        {
            for (int c = 0; c < 3; ++c)
            {
                gsl_matrix_set(cov, r, c, gsl_matrix_get(cov, r, c) + centered(r, 0) * centered(c, 0));
            }
        }
    }

    gsl_vector *eval = gsl_vector_alloc(3);
    gsl_matrix *evec = gsl_matrix_alloc(3, 3);
    gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(3);
    gsl_eigen_symmv(cov, eval, evec, workspace);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

    Matrix normal(3, 1);
    for (int i = 0; i < 3; ++i)
    {
        normal.set(i, 0, gsl_matrix_get(evec, i, 0));
    }
    int highestIndex = 0;
    double highestZ = param.scaffoldingPoints[0](2, 0);
    for (int i = 1; i < static_cast<int>(param.scaffoldingPoints.size()); ++i)
    {
        if (param.scaffoldingPoints[i](2, 0) > highestZ)
        {
            highestZ = param.scaffoldingPoints[i](2, 0);
            highestIndex = i;
        }
    }
    const Matrix apexVec = param.scaffoldingPoints[highestIndex] - centroid;
    if (dot_col(normal, apexVec) < 0.0)
    {
        normal *= -1.0;
    }

    Matrix zAxis(3, 1);
    zAxis.set(2, 0, 1.0);
    gagInitialAlignmentRotation = rotation_aligning_vectors(normal, zAxis);

    for (Matrix &point : param.scaffoldingPoints)
    {
        point = gagInitialAlignmentRotation * (point - centroid) + centroid;
    }

    int capCenterIndex = 0;
    double capCenterZ = param.scaffoldingPoints[0](2, 0);
    for (int i = 1; i < static_cast<int>(param.scaffoldingPoints.size()); ++i)
    {
        if (param.scaffoldingPoints[i](2, 0) > capCenterZ)
        {
            capCenterZ = param.scaffoldingPoints[i](2, 0);
            capCenterIndex = i;
        }
    }
    const double xShiftToMembraneCenter = -param.scaffoldingPoints[capCenterIndex](0, 0);
    const double yShiftToMembraneCenter = -param.scaffoldingPoints[capCenterIndex](1, 0);
    for (Matrix &point : param.scaffoldingPoints)
    {
        point.set(0, 0, point(0, 0) + xShiftToMembraneCenter);
        point.set(1, 0, point(1, 0) + yShiftToMembraneCenter);
    }

    gsl_eigen_symmv_free(workspace);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(cov);
    return true;
}

bool Mesh::initialize_gag_scaffolding_topology()
{
    const bool useGagMode = param.isGagScaffoldingEnergyIncluded;
    const bool useIdealizedMode = param.isIdealizedProteinLatticeEnergyIncluded;
    if (!useGagMode && !useIdealizedMode)
    {
        return false;
    }
    if (useGagMode && useIdealizedMode)
    {
        throw std::runtime_error("Enable only one internal scaffold mode at a time: Gag or idealized protein lattice.");
    }
    if (gagScaffoldingTopologyInitialized)
    {
        return true;
    }

    const std::string &referenceStateFile = useGagMode
                                          ? param.gagReferenceStateFileName
                                          : param.idealizedProteinLatticeFileName;
    if (referenceStateFile.empty())
    {
        throw std::runtime_error(useGagMode
            ? "Gag scaffold energy enabled but gagReferenceStateFileName is empty."
            : "Idealized protein lattice energy enabled but idealizedProteinLatticeFileName is empty.");
    }

    std::unordered_map<int, ReferenceMolecule> referenceMolecules;
    std::vector<ReferenceBoundState> referenceBoundStates;
    if (!parse_reference_state_file(referenceStateFile, referenceMolecules, referenceBoundStates))
    {
        throw std::runtime_error("Failed to parse scaffold reference state file: " + referenceStateFile);
    }

    const std::unordered_map<std::string, GagReactionTarget> reactionTargets =
        useGagMode ? load_gag_reaction_targets(param.gagReactionFileName)
                   : std::unordered_map<std::string, GagReactionTarget>();

    std::vector<int> referenceMoleculeIds;
    referenceMoleculeIds.reserve(referenceMolecules.size());
    for (const auto &entry : referenceMolecules)
    {
        referenceMoleculeIds.push_back(entry.first);
    }
    std::sort(referenceMoleculeIds.begin(), referenceMoleculeIds.end());

    if (referenceMoleculeIds.size() != param.scaffoldingPoints.size())
    {
        throw std::runtime_error("Reference Gag state size does not match imported scaffolding COM count.");
    }

    std::vector<int> referenceToPoint(referenceMoleculeIds.size(), -1);
    std::unordered_map<int, int> pointIndexToReferenceMol;
    std::vector<bool> pointTaken(param.scaffoldingPoints.size(), false);
    for (size_t refIdx = 0; refIdx < referenceMoleculeIds.size(); ++refIdx)
    {
        const Matrix &referencePoint = referenceMolecules[referenceMoleculeIds[refIdx]].com;
        double bestDistance = std::numeric_limits<double>::max();
        int bestPointIndex = -1;
        for (int pointIdx = 0; pointIdx < static_cast<int>(param.scaffoldingPoints.size()); ++pointIdx)
        {
            if (pointTaken[pointIdx])
            {
                continue;
            }
            const double distance = (param.scaffoldingPoints[pointIdx] - referencePoint).calculate_norm();
            if (distance < bestDistance)
            {
                bestDistance = distance;
                bestPointIndex = pointIdx;
            }
        }
        if (bestPointIndex < 0)
        {
            throw std::runtime_error("Failed to map reference Gag COM to current scaffolding COM.");
        }
        pointTaken[bestPointIndex] = true;
        referenceToPoint[refIdx] = bestPointIndex;
        pointIndexToReferenceMol[bestPointIndex] = referenceMoleculeIds[refIdx];
    }

    std::unordered_map<int, int> moleculeIdToPointIndex;
    for (size_t refIdx = 0; refIdx < referenceMoleculeIds.size(); ++refIdx)
    {
        moleculeIdToPointIndex[referenceMoleculeIds[refIdx]] = referenceToPoint[refIdx];
    }

    gagBonds.clear();
    gagAngles.clear();
    gagTorsions.clear();
    gagSubunits.assign(param.scaffoldingPoints.size(), GagSubunit());
    gagInteractions.clear();

    if (gagInitialAlignmentRotation.calculate_quad() == 0.0)
    {
        gagInitialAlignmentRotation = identity3();
    }

    for (const int moleculeId : referenceMoleculeIds)
    {
        const int pointIndex = moleculeIdToPointIndex[moleculeId];
        GagSubunit subunit;
        subunit.pointIndex = pointIndex;
        subunit.moleculeId = moleculeId;
        subunit.rotation = gagInitialAlignmentRotation;
        for (const auto &ifaceEntry : referenceMolecules[moleculeId].interfaces)
        {
            subunit.localInterfaces[ifaceEntry.first] = ifaceEntry.second - referenceMolecules[moleculeId].com;
        }
        if (referenceMolecules[moleculeId].hasReferenceVector)
        {
            subunit.localReferenceVector = referenceMolecules[moleculeId].referenceVector;
            subunit.hasReferenceVector = true;
        }
        gagSubunits[pointIndex] = subunit;
    }

    std::set<std::pair<int, int>> seenPairs;
    int nonColinearInteractionCount = 0;
    double maxColinearityOffset = 0.0;
    for (const ReferenceBoundState &boundState : referenceBoundStates)
    {
        const int point1 = moleculeIdToPointIndex[boundState.mol1];
        const int point2 = moleculeIdToPointIndex[boundState.mol2];
        const std::pair<int, int> pairKey = std::minmax(point1, point2);
        if (seenPairs.count(pairKey) > 0)
        {
            continue;
        }
        seenPairs.insert(pairKey);

        GagBond bond;
        bond.point1 = point1;
        bond.point2 = point2;
        bond.interface1 = boundState.iface1;
        bond.interface2 = boundState.iface2;
        bond.expectedLength = (referenceMolecules[boundState.mol1].com - referenceMolecules[boundState.mol2].com).calculate_norm();
        const std::string key = bond.interface1 + "|" + bond.interface2;
        const auto reactionIt = reactionTargets.find(key);
        if (useGagMode && reactionIt != reactionTargets.end())
        {
            bond.thetaAtPoint1 = reactionIt->second.assocAngles[0];
            bond.thetaAtPoint2 = reactionIt->second.assocAngles[1];
            bond.expectedPhi = 0.5 * (reactionIt->second.assocAngles[2] + reactionIt->second.assocAngles[3]);
            bond.expectedOmega = reactionIt->second.assocAngles[4];
        }
        gagBonds.push_back(bond);

        GagInteraction interaction;
        interaction.subunit1 = point1;
        interaction.subunit2 = point2;
        interaction.interface1 = boundState.iface1;
        interaction.interface2 = boundState.iface2;
        if (useGagMode && reactionIt != reactionTargets.end())
        {
            interaction.expectedSigma = reactionIt->second.sigma;
            interaction.assocAngles = reactionIt->second.assocAngles;
        }
        else
        {
            const GagSubunit &subunit1 = gagSubunits[point1];
            const GagSubunit &subunit2 = gagSubunits[point2];
            const Matrix interface1 = current_interface_coord(subunit1, boundState.iface1, param.scaffoldingPoints);
            const Matrix interface2 = current_interface_coord(subunit2, boundState.iface2, param.scaffoldingPoints);
            const Matrix com1 = param.scaffoldingPoints[subunit1.pointIndex];
            const Matrix com2 = param.scaffoldingPoints[subunit2.pointIndex];
            const Matrix comCom = com2 - com1;
            const double comComNorm = comCom.calculate_norm();
            if (useIdealizedMode && comComNorm > 0.0)
            {
                const double offset1 = safe_cross(interface1 - com1, comCom).calculate_norm() / comComNorm;
                const double offset2 = safe_cross(interface2 - com2, comCom).calculate_norm() / comComNorm;
                const double worstOffset = std::max(offset1, offset2);
                maxColinearityOffset = std::max(maxColinearityOffset, worstOffset);
                if (worstOffset > IDEALIZED_COLINEARITY_WARNING_THRESHOLD)
                {
                    ++nonColinearInteractionCount;
                }
            }
            interaction.expectedSigma = (interface2 - interface1).calculate_norm();
            interaction.assocAngles[0] = safe_angle_from_vectors(com1 - interface1, interface2 - interface1);
            interaction.assocAngles[1] = safe_angle_from_vectors(com2 - interface2, interface1 - interface2);
            interaction.assocAngles[2] = 0.0;
            interaction.assocAngles[3] = 0.0;
            if (subunit1.hasReferenceVector && subunit2.hasReferenceVector)
            {
                const Matrix ref1 = current_reference_endpoint(subunit1, param.scaffoldingPoints);
                const Matrix ref2 = current_reference_endpoint(subunit2, param.scaffoldingPoints);
                interaction.assocAngles[4] = safe_dihedral(ref1, com1, com2, ref2);
            }
            else
            {
                interaction.assocAngles[4] = 0.0;
            }
        }
        gagInteractions.push_back(interaction);

        if (useGagMode && !reactionTargets.empty())
        {
            if (reactionTargets.find(key) == reactionTargets.end() && param.VERBOSE_MODE)
            {
                std::cout << "[Mesh::initialize_gag_scaffolding_topology] Missing reaction target for "
                          << key << " in " << param.gagReactionFileName << std::endl;
            }
        }
    }

    if (useIdealizedMode && nonColinearInteractionCount > 0)
    {
        std::cout << "[Mesh::initialize_gag_scaffolding_topology] Warning: idealized protein lattice mode assumes interfaces are nearly colinear with COM-COM bonds."
                  << " Threshold = " << IDEALIZED_COLINEARITY_WARNING_THRESHOLD
                  << ", offending interactions = " << nonColinearInteractionCount
                  << ", max offset = " << maxColinearityOffset << std::endl;
    }

    gagScaffoldingTopologyInitialized = true;
    return true;
}

double Mesh::calculate_gag_scaffolding_internal_energy() const
{
    if (!param.isGagScaffoldingEnergyIncluded && !param.isIdealizedProteinLatticeEnergyIncluded)
    {
        return 0.0;
    }
    if (!gagScaffoldingTopologyInitialized)
    {
        const_cast<Mesh *>(this)->initialize_gag_scaffolding_topology();
    }
    return calculate_gag_scaffolding_internal_energy(param.scaffoldingPoints, gagSubunits);
}

double Mesh::calculate_gag_scaffolding_internal_energy(const std::vector<Matrix> &scaffoldingPoints) const
{
    return calculate_gag_scaffolding_internal_energy(scaffoldingPoints, gagSubunits);
}

double Mesh::calculate_gag_scaffolding_internal_energy(const std::vector<Matrix> &scaffoldingPoints,
                                                       const std::vector<GagSubunit> &subunits) const
{
    if ((!param.isGagScaffoldingEnergyIncluded && !param.isIdealizedProteinLatticeEnergyIncluded) ||
        !gagScaffoldingTopologyInitialized)
    {
        return 0.0;
    }

    double totalEnergy = 0.0;
    const bool useIdealizedMode = param.isIdealizedProteinLatticeEnergyIncluded;
#ifdef OMP
#pragma omp parallel for reduction(+ : totalEnergy)
#endif
    for (int interactionIndex = 0; interactionIndex < static_cast<int>(gagInteractions.size()); ++interactionIndex)
    {
        const GagInteraction &interaction = gagInteractions[interactionIndex];
        const GagSubunit &subunit1 = subunits[interaction.subunit1];
        const GagSubunit &subunit2 = subunits[interaction.subunit2];
        const Matrix interface1 = current_interface_coord(subunit1, interaction.interface1, scaffoldingPoints);
        const Matrix interface2 = current_interface_coord(subunit2, interaction.interface2, scaffoldingPoints);
        const Matrix com1 = scaffoldingPoints[subunit1.pointIndex];
        const Matrix com2 = scaffoldingPoints[subunit2.pointIndex];

        const Matrix sigmaVec = interface2 - interface1;
        const double sigma = sigmaVec.calculate_norm();
        if (!std::isfinite(sigma))
        {
            continue;
        }
        const double sigmaDelta = sigma - interaction.expectedSigma;
        totalEnergy += 0.5 * param.gagKsigma * sigmaDelta * sigmaDelta;

        const double theta1 = safe_angle_from_vectors(com1 - interface1, interface2 - interface1);
        const double theta2 = safe_angle_from_vectors(com2 - interface2, interface1 - interface2);
        totalEnergy += 0.5 * param.gagKtheta * std::pow(theta1 - interaction.assocAngles[0], 2);
        totalEnergy += 0.5 * param.gagKtheta * std::pow(theta2 - interaction.assocAngles[1], 2);

        if (useIdealizedMode)
        {
            if (subunit1.hasReferenceVector && subunit2.hasReferenceVector)
            {
                const Matrix ref1 = current_reference_endpoint(subunit1, scaffoldingPoints);
                const Matrix ref2 = current_reference_endpoint(subunit2, scaffoldingPoints);
                const double omega = safe_dihedral(ref1, com1, com2, ref2);
                const double omegaDelta = wrapped_angle_difference(omega, interaction.assocAngles[4]);
                totalEnergy += 0.5 * param.gagKomega * omegaDelta * omegaDelta;
            }
        }
        else
        {
            const double torsion = safe_dihedral(com1, interface1, interface2, com2);
            const double phi1Delta = wrapped_angle_difference(torsion, interaction.assocAngles[2]);
            const double phi2Delta = wrapped_angle_difference(torsion, interaction.assocAngles[3]);
            const double omegaDelta = wrapped_angle_difference(torsion, interaction.assocAngles[4]);
            totalEnergy += 0.5 * param.gagKphi * phi1Delta * phi1Delta;
            totalEnergy += 0.5 * param.gagKphi * phi2Delta * phi2Delta;
            totalEnergy += 0.5 * param.gagKomega * omegaDelta * omegaDelta;
        }
    }

    return totalEnergy;
}

std::vector<Matrix> Mesh::calculate_gag_scaffolding_forces_fd() const
{
    std::vector<Matrix> forces(param.scaffoldingPoints.size(), Matrix(3, 1));
    if (!param.isGagScaffoldingEnergyIncluded && !param.isIdealizedProteinLatticeEnergyIncluded)
    {
        return forces;
    }
    if (!gagScaffoldingTopologyInitialized)
    {
        const_cast<Mesh *>(this)->initialize_gag_scaffolding_topology();
    }

    const double step = param.gagFiniteDifferenceStep;
    std::vector<Matrix> plusPoints = param.scaffoldingPoints;
    std::vector<Matrix> minusPoints = param.scaffoldingPoints;
    for (size_t pointIndex = 0; pointIndex < param.scaffoldingPoints.size(); ++pointIndex)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            plusPoints[pointIndex].set(axis, 0, param.scaffoldingPoints[pointIndex](axis, 0) + step);
            minusPoints[pointIndex].set(axis, 0, param.scaffoldingPoints[pointIndex](axis, 0) - step);
            const double ePlus = calculate_gag_scaffolding_internal_energy(plusPoints);
            const double eMinus = calculate_gag_scaffolding_internal_energy(minusPoints);
            forces[pointIndex].set(axis, 0, -(ePlus - eMinus) / (2.0 * step));
            plusPoints[pointIndex].set(axis, 0, param.scaffoldingPoints[pointIndex](axis, 0));
            minusPoints[pointIndex].set(axis, 0, param.scaffoldingPoints[pointIndex](axis, 0));
        }
    }
    return forces;
}

std::vector<Matrix> Mesh::calculate_gag_subunit_translation_forces_fd()
{
    std::vector<Matrix> forces(param.scaffoldingPoints.size(), Matrix(3, 1));
    if (!param.isGagScaffoldingEnergyIncluded && !param.isIdealizedProteinLatticeEnergyIncluded)
    {
        return forces;
    }
    if (!gagScaffoldingTopologyInitialized)
    {
        initialize_gag_scaffolding_topology();
    }

    const double step = param.gagFiniteDifferenceStep;
    const std::vector<Matrix> basePoints = param.scaffoldingPoints;
#ifdef OMP
#pragma omp parallel for
#endif
    for (int pointIndex = 0; pointIndex < static_cast<int>(param.scaffoldingPoints.size()); ++pointIndex)
    {
        std::vector<Matrix> perturbedPoints = basePoints;
        for (int axis = 0; axis < 3; ++axis)
        {
            perturbedPoints[pointIndex].set(axis, 0, basePoints[pointIndex](axis, 0) + step);
            const double ePlus = calculate_gag_scaffolding_internal_energy(perturbedPoints, gagSubunits);
            perturbedPoints[pointIndex].set(axis, 0, basePoints[pointIndex](axis, 0) - step);
            const double eMinus = calculate_gag_scaffolding_internal_energy(perturbedPoints, gagSubunits);
            perturbedPoints[pointIndex].set(axis, 0, basePoints[pointIndex](axis, 0));
            forces[pointIndex].set(axis, 0, -(ePlus - eMinus) / (2.0 * step));
        }
    }
    return forces;
}

std::vector<Matrix> Mesh::calculate_gag_subunit_rotation_torques_fd()
{
    std::vector<Matrix> torques(param.scaffoldingPoints.size(), Matrix(3, 1));
    if (!param.isGagScaffoldingEnergyIncluded && !param.isIdealizedProteinLatticeEnergyIncluded)
    {
        return torques;
    }
    if (!gagScaffoldingTopologyInitialized)
    {
        initialize_gag_scaffolding_topology();
    }

    const double step = param.gagFiniteDifferenceStep;
    const std::vector<Matrix> basePoints = param.scaffoldingPoints;
    const std::vector<GagSubunit> baseSubunits = gagSubunits;
#ifdef OMP
#pragma omp parallel for
#endif
    for (int pointIndex = 0; pointIndex < static_cast<int>(gagSubunits.size()); ++pointIndex)
    {
        std::vector<GagSubunit> perturbedSubunits = baseSubunits;
        const Matrix baseRotation = baseSubunits[pointIndex].rotation;
        for (int axis = 0; axis < 3; ++axis)
        {
            Matrix axisVector(3, 1);
            axisVector.set(axis, 0, 1.0);
            Matrix rotPlus = rotation_matrix_from_axis_angle(axisVector, step);
            Matrix rotMinus = rotation_matrix_from_axis_angle(axisVector, -step);
            perturbedSubunits[pointIndex].rotation = rotPlus * baseRotation;
            const double ePlus = calculate_gag_scaffolding_internal_energy(basePoints, perturbedSubunits);
            perturbedSubunits[pointIndex].rotation = rotMinus * baseRotation;
            const double eMinus = calculate_gag_scaffolding_internal_energy(basePoints, perturbedSubunits);
            perturbedSubunits[pointIndex].rotation = baseRotation;
            torques[pointIndex].set(axis, 0, -(ePlus - eMinus) / (2.0 * step));
        }
    }
    return torques;
}
