#include "io/io.hpp"
#include <cstdio>
#include <iomanip>

using namespace std;

namespace
{
void write_matrix3(std::ostream &out, const Matrix &matrix)
{
    out << matrix(0, 0) << ' ' << matrix(1, 0) << ' ' << matrix(2, 0);
}

bool read_matrix3(std::istream &in, Matrix &matrix)
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    if (!(in >> x >> y >> z))
    {
        return false;
    }
    matrix.set(0, 0, x);
    matrix.set(1, 0, y);
    matrix.set(2, 0, z);
    return true;
}

void write_energy_terms(std::ostream &out, const Energy &energy)
{
    out << energy.energyCurvature << ' '
        << energy.energyArea << ' '
        << energy.energyVolume << ' '
        << energy.energyThickness << ' '
        << energy.energyTilt << ' '
        << energy.energyRegularization << ' '
        << energy.energyHarmonicBond << ' '
        << energy.energyGagScaffolding << ' '
        << energy.energyIdealizedProteinLattice << ' '
        << energy.energyTotal;
}

bool read_energy_terms(std::istream &in, Energy &energy)
{
    return static_cast<bool>(in >> energy.energyCurvature
                                >> energy.energyArea
                                >> energy.energyVolume
                                >> energy.energyThickness
                                >> energy.energyTilt
                                >> energy.energyRegularization
                                >> energy.energyHarmonicBond
                                >> energy.energyGagScaffolding
                                >> energy.energyIdealizedProteinLattice
                                >> energy.energyTotal);
}
}

/**
 * @brief Writes the vertex data for the mesh to a CSV file.
 *
 * This function writes the coordinates of each vertex in the mesh to a CSV file named "vertex[nFile].csv", where nFile is
 * equal to iteration/100. The output is formatted as follows:
 *
 *     x1,y1,z1
 *     x2,y2,z2
 *     ...
 *
 * where xi, yi, and zi are the coordinates of the ith vertex in the mesh.
 *
 * @note Step size and file name are left in definition in case customization is needed.
 *
 * @param mesh The mesh whose vertex data should be written to a CSV file.
 * @param iteration The current iteration number of the optimization algorithm.
 */
void write_vertex_data_to_csv(const Mesh &mesh, const int iteration)
{
    const int STEP_SIZE = 100;
    const std::string VERTEX_FILENAME_PREFIX = "vertex";
    const std::string VERTEX_FILENAME_SUFFIX = ".csv";
    const int nFile = iteration / STEP_SIZE;
    const std::string filename = VERTEX_FILENAME_PREFIX + std::to_string(nFile) + VERTEX_FILENAME_SUFFIX;

    std::ofstream outfile(filename);
    for (const auto &vertex : mesh.vertices)
    {
        outfile << vertex.coord(0, 0) << ',' << vertex.coord(1, 0) << ',' << vertex.coord(2, 0) << '\n';
    }
    outfile.close();
}

/**
 * @brief Writes the energy force data to a CSV file.
 *
 * This function writes the energy and force data for each iteration in the following format:
 *
 * "E_Curvature, E_Area, E_Regularization, E_Total ((pN.nm)), Mean Force (pN)"
 *
 * @note File name is left in definition in case customization is needed.
 *
 * @param model The model that cotains a mesh whose vertex data should be written to a CSV file.
 */
void write_energy_force_data_to_csv(const Model &model)
{
    // Output Energy and meanforce
    std::string ENERGY_FORCE_FILENAME = "EnergyForce.csv";
    ofstream outfileEF(ENERGY_FORCE_FILENAME);
    for (int j = 0; j <= model.iteration; j++)
    {
        if (j == 0)
        {
            outfileEF << "E_Curvature, E_Area, E_Regularization, E_HarmonicBond, E_GagScaffolding, E_IdealizedProteinLattice, E_Total ((pN.nm)), Mean Force (pN)" << '\n';
        }
        outfileEF << model.record.energyVec[j].energyCurvature << ", "
                  << model.record.energyVec[j].energyArea << ", "
                  << model.record.energyVec[j].energyRegularization << ", "
                  << model.record.energyVec[j].energyHarmonicBond << ", "
                  << model.record.energyVec[j].energyGagScaffolding << ", "
                  << model.record.energyVec[j].energyIdealizedProteinLattice << ", "
                  << model.record.energyVec[j].energyTotal << ", "
                  << model.record.meanForce[j] << '\n';
    }
    outfileEF.close();
}

bool write_model_restart_checkpoint(const Model &model,
                                    const std::string &filepath,
                                    const int nextIteration)
{
    if (filepath.empty())
    {
        return false;
    }

    const std::string tempFilepath = filepath + ".tmp";
    std::ofstream outfile(tempFilepath);
    if (!outfile.is_open())
    {
        std::cerr << "[write_model_restart_checkpoint] Could not open "
                  << tempFilepath << " for writing." << std::endl;
        return false;
    }

    outfile << std::setprecision(17);
    outfile << "SLIMED_RESTART_V1\n";
    outfile << "nextIteration " << nextIteration << "\n";
    outfile << "stepSize " << model.stepSize << "\n";
    outfile << "oa "
            << model.oa.usingNCG << ' '
            << model.oa.isNCGstuck << ' '
            << model.oa.usingRpi << ' '
            << model.oa.isCriteriaSatisfied << ' '
            << model.oa.trialStepSize << ' '
            << model.oa.trialIterationInterval << ' '
            << model.oa.c1 << ' '
            << model.oa.c2 << ' '
            << model.oa.stepThreshold << ' '
            << model.oa.nConsecutiveNcgStuck << ' '
            << model.oa.nConsecutiveNcgStuckThreshold << ' '
            << model.oa.energyDiffThreshold << ' '
            << model.oa.forceDiffThreshold << "\n";
    outfile << "thermal "
            << model.isHeating << ' '
            << model.currentCoolingStep << ' '
            << model.coolingStep << ' '
            << model.currentHeatingStep << ' '
            << model.heatingStep << ' '
            << model.highTemperature << ' '
            << model.thermalFluctuationAttemptCount << "\n";
    outfile << "thermalRng " << model.thermalRng << "\n";
    outfile << "springConst " << model.mesh.param.springConst << "\n";
    outfile << "scaffoldingMovement ";
    write_matrix3(outfile, model.mesh.scaffoldingMovementVector);
    outfile << "\n";

    outfile << "vertices " << model.mesh.vertices.size() << "\n";
    for (const Vertex &vertex : model.mesh.vertices)
    {
        outfile << vertex.index << ' ';
        write_matrix3(outfile, vertex.coord);
        outfile << ' ';
        write_matrix3(outfile, vertex.coordPrev);
        outfile << ' ';
        write_matrix3(outfile, vertex.coordRef);
        outfile << ' ';
        write_matrix3(outfile, vertex.force.forceTotal);
        outfile << ' ';
        write_matrix3(outfile, vertex.forcePrev.forceTotal);
        outfile << ' ';
        write_matrix3(outfile, model.ncgDirection0[vertex.index].forceTotal);
        outfile << "\n";
    }

    outfile << "scaffoldingPoints " << model.mesh.param.scaffoldingPoints.size() << "\n";
    for (int i = 0; i < static_cast<int>(model.mesh.param.scaffoldingPoints.size()); ++i)
    {
        outfile << i << ' ';
        write_matrix3(outfile, model.mesh.param.scaffoldingPoints[i]);
        outfile << "\n";
    }

    outfile << "gagRotations " << model.mesh.gagSubunits.size() << "\n";
    for (int i = 0; i < static_cast<int>(model.mesh.gagSubunits.size()); ++i)
    {
        outfile << i;
        for (int row = 0; row < 3; ++row)
        {
            for (int col = 0; col < 3; ++col)
            {
                outfile << ' ' << model.mesh.gagSubunits[i].rotation(row, col);
            }
        }
        outfile << "\n";
    }

    outfile << "records " << model.record.energyVec.size() << "\n";
    for (int i = 0; i < static_cast<int>(model.record.energyVec.size()); ++i)
    {
        outfile << i << ' ' << model.record.areaTotal[i] << ' ';
        write_energy_terms(outfile, model.record.energyVec[i]);
        outfile << ' ' << model.record.meanForce[i] << "\n";
    }
    outfile << "END\n";
    outfile.close();

    if (!outfile)
    {
        std::cerr << "[write_model_restart_checkpoint] Failed while writing "
                  << tempFilepath << "." << std::endl;
        return false;
    }

    if (std::rename(tempFilepath.c_str(), filepath.c_str()) != 0)
    {
        std::cerr << "[write_model_restart_checkpoint] Could not replace "
                  << filepath << "." << std::endl;
        return false;
    }

    std::cout << "[write_model_restart_checkpoint] Wrote checkpoint for next iteration "
              << nextIteration << " to " << filepath << std::endl;
    return true;
}

bool load_model_restart_checkpoint(Model &model, const std::string &filepath)
{
    std::ifstream infile(filepath);
    if (!infile.is_open())
    {
        std::cerr << "[load_model_restart_checkpoint] Could not open "
                  << filepath << " for reading." << std::endl;
        return false;
    }

    std::string tag;
    infile >> tag;
    if (tag != "SLIMED_RESTART_V1")
    {
        std::cerr << "[load_model_restart_checkpoint] Unsupported checkpoint format in "
                  << filepath << "." << std::endl;
        return false;
    }

    infile >> tag >> model.iteration;
    if (tag != "nextIteration")
    {
        return false;
    }
    infile >> tag >> model.stepSize;
    if (tag != "stepSize")
    {
        return false;
    }

    infile >> tag;
    if (tag != "oa")
    {
        return false;
    }
    infile >> model.oa.usingNCG
           >> model.oa.isNCGstuck
           >> model.oa.usingRpi
           >> model.oa.isCriteriaSatisfied
           >> model.oa.trialStepSize
           >> model.oa.trialIterationInterval
           >> model.oa.c1
           >> model.oa.c2
           >> model.oa.stepThreshold
           >> model.oa.nConsecutiveNcgStuck
           >> model.oa.nConsecutiveNcgStuckThreshold
           >> model.oa.energyDiffThreshold
           >> model.oa.forceDiffThreshold;

    infile >> tag;
    if (tag != "thermal")
    {
        return false;
    }
    infile >> model.isHeating
           >> model.currentCoolingStep
           >> model.coolingStep
           >> model.currentHeatingStep
           >> model.heatingStep
           >> model.highTemperature
           >> model.thermalFluctuationAttemptCount;

    infile >> tag;
    if (tag != "thermalRng")
    {
        return false;
    }
    infile >> model.thermalRng;

    infile >> tag >> model.mesh.param.springConst;
    if (tag != "springConst")
    {
        return false;
    }

    infile >> tag;
    if (tag != "scaffoldingMovement" ||
        !read_matrix3(infile, model.mesh.scaffoldingMovementVector))
    {
        return false;
    }

    int count = 0;
    infile >> tag >> count;
    if (tag != "vertices" || count != static_cast<int>(model.mesh.vertices.size()))
    {
        std::cerr << "[load_model_restart_checkpoint] Vertex count mismatch." << std::endl;
        return false;
    }
    for (int row = 0; row < count; ++row)
    {
        int index = -1;
        infile >> index;
        if (index < 0 || index >= static_cast<int>(model.mesh.vertices.size()))
        {
            return false;
        }
        Vertex &vertex = model.mesh.vertices[index];
        if (!read_matrix3(infile, vertex.coord) ||
            !read_matrix3(infile, vertex.coordPrev) ||
            !read_matrix3(infile, vertex.coordRef) ||
            !read_matrix3(infile, vertex.force.forceTotal) ||
            !read_matrix3(infile, vertex.forcePrev.forceTotal) ||
            !read_matrix3(infile, model.ncgDirection0[index].forceTotal))
        {
            return false;
        }
    }

    infile >> tag >> count;
    if (tag != "scaffoldingPoints" ||
        count != static_cast<int>(model.mesh.param.scaffoldingPoints.size()))
    {
        std::cerr << "[load_model_restart_checkpoint] Scaffolding point count mismatch." << std::endl;
        return false;
    }
    for (int row = 0; row < count; ++row)
    {
        int index = -1;
        infile >> index;
        if (index < 0 || index >= static_cast<int>(model.mesh.param.scaffoldingPoints.size()) ||
            !read_matrix3(infile, model.mesh.param.scaffoldingPoints[index]))
        {
            return false;
        }
    }

    infile >> tag >> count;
    if (tag != "gagRotations" || count != static_cast<int>(model.mesh.gagSubunits.size()))
    {
        std::cerr << "[load_model_restart_checkpoint] Gag rotation count mismatch." << std::endl;
        return false;
    }
    for (int row = 0; row < count; ++row)
    {
        int index = -1;
        infile >> index;
        if (index < 0 || index >= static_cast<int>(model.mesh.gagSubunits.size()))
        {
            return false;
        }
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                double value = 0.0;
                infile >> value;
                model.mesh.gagSubunits[index].rotation.set(i, j, value);
            }
        }
    }

    infile >> tag >> count;
    if (tag != "records" || count < 0)
    {
        return false;
    }
    model.record.areaTotal.clear();
    model.record.energyVec.clear();
    model.record.meanForce.clear();
    model.record.areaTotal.reserve(count);
    model.record.energyVec.reserve(count);
    model.record.meanForce.reserve(count);
    for (int row = 0; row < count; ++row)
    {
        int index = -1;
        double area = 0.0;
        Energy energy;
        double meanForce = 0.0;
        infile >> index >> area;
        if (index != row || !read_energy_terms(infile, energy) || !(infile >> meanForce))
        {
            return false;
        }
        model.record.add(area, energy, meanForce);
    }

    infile >> tag;
    if (tag != "END")
    {
        return false;
    }

    std::cout << "[load_model_restart_checkpoint] Restarted from " << filepath
              << " at iteration " << model.iteration << std::endl;
    return true;
}

/**
 * @brief Writes the current iteration element face energy data to a CSV file.
 *
 * This function writes the energy data for each element face of
 *  the current iteration in the following format:
 *
 * "E_Curvature, E_Area, E_Regularization, E_Total"
 *
 * @note @note File name is left in definition in case customization is needed.
 *
 * @param model The model that cotains a mesh whose vertex data should be written to a CSV file.
 *
 * 
 */
void write_element_face_energy_to_csv(const Model &model)
{
    // Output face index and then its corresponding element face energies
    std::string ELEMENT_FACE_ENERGY_FILENAME = "ElementFaceEnergy.csv";
    ofstream outfileEFE(ELEMENT_FACE_ENERGY_FILENAME);
    
    // Add header
    outfileEFE << "Face_index, E_Curvature, E_Area, E_Regularization, E_Total" << "\n";

    // iterate thru face and append to the output file
    for (Face& face: model.mesh.faces)
    {
        outfileEFE << face.index << ","
                << face.energy.energyCurvature << ","
                << face.energy.energyArea << ","
                << face.energy.energyRegularization << ","
                << face.energy.energyTotal << "\n";
    }
    outfileEFE.close();
}

/**
 * @brief Writes the last vertex data for the mesh to a CSV file.
 *
 * This function writes the coordinates of each vertex in the last step
 * in the mesh to a CSV file named "vertexfinal.csv". The output is formatted as follows:
 *
 *     x1,y1,z1
 *     x2,y2,z2
 *     ...
 *
 * where xi, yi, and zi are the coordinates of the ith vertex in the mesh.
 *
 * @note file name are left in definition in case customization is needed.
 *
 * @param mesh The mesh whose vertex data should be written to a CSV file.
 */
void write_final_vertex_data_to_csv(const Mesh &mesh)
{
    const std::string FINAL_VERTEX_FILENANME = "vertexfinal.csv";

    ofstream outfile_final(FINAL_VERTEX_FILENANME);
    for (int j = 0; j < mesh.vertices.size(); j++)
    {
        outfile_final << setprecision(16) << mesh.vertices[j].coordPrev(0, 0) << ',' << mesh.vertices[j].coordPrev(1, 0) << ',' << mesh.vertices[j].coordPrev(2, 0) << '\n';
    }
    outfile_final.close();
}

void dynamics_create_trajectory_files(DynamicMesh &mesh, const std::string &filename)
{
        ofstream surfacepoint_csv("surfacepoint" + filename + ".csv");
        ofstream meshpoint_csv("meshpoint" + filename + ".csv");

        for (int i = 0; i < mesh.vertices.size(); i++)
        {
            meshpoint_csv << setprecision(8) << mesh.matMesh(i, 0) << ',' << mesh.matMesh(i, 1) << ',' << mesh.matMesh(i, 2) << ',';
            surfacepoint_csv << setprecision(8) << mesh.matSurface(i, 0) << ',' << mesh.matSurface(i, 1) << ',' << mesh.matSurface(i, 2) << ',';
        }
        surfacepoint_csv << '\n';
        meshpoint_csv << '\n';
        surfacepoint_csv.close();
        meshpoint_csv.close();
}

void output_trajectory_files(Mesh &mesh, const std::string &input_filename)
{
    ofstream meshpoint_csv_temp;
    meshpoint_csv_temp.open("meshpoint" + input_filename + ".csv", std::ios::app);

    for (int i = 0; i < mesh.vertices.size(); i++) {
        meshpoint_csv_temp << setprecision(8) << mesh.vertices[i].coord(0, 0) << ","
        << mesh.vertices[i].coord(1, 0) << ","
        << mesh.vertices[i].coord(2, 0) << ",";
    }
    meshpoint_csv_temp << '\n';
    meshpoint_csv_temp.close();
}

void dynamics_output_trajectory_files(DynamicMesh &mesh, const std::string &filename)
{
    if (mesh.param.meshpointOutput) {
        ofstream surfacepoint_csv_temp;
        ofstream meshpoint_csv_temp;
        surfacepoint_csv_temp.open("surfacepoint" + filename + ".csv", std::ios::app);
        meshpoint_csv_temp.open("meshpoint" + filename + ".csv", std::ios::app);
        
        for (int i = 0; i < mesh.vertices.size(); i++) {
            meshpoint_csv_temp << setprecision(8) << mesh.matMesh(i, 0) << ',' << mesh.matMesh(i, 1) << ',' << mesh.matMesh(i, 2) << ',';
            surfacepoint_csv_temp << setprecision(8) << mesh.matSurface(i, 0) << ',' << mesh.matSurface(i, 1) << ',' << mesh.matSurface(i, 2) << ',';
        }
        meshpoint_csv_temp << '\n';
        surfacepoint_csv_temp << '\n';
        meshpoint_csv_temp.close();
        surfacepoint_csv_temp.close();
    }
}

/*
std::string format_spaces_decimals(double number) {
    // Limit the values if they exceed the maximum/minimum possible value.
    if (number > 9999) {
        number = 9999.99999;
    }
    if (number < -9999) {
        number = -9999.99999;
    }

    // Format the number to 5 decimal places and 13 characters long.
    std::stringstream ss;
    ss.precision(5);
    ss << std::fixed << std::setw(13) << number;
    return ss.str();
}

// Function to convert a DataFrame to an XYZ string.
std::string to_xyz_string(const std::vector<std::vector<double>>& df) {
    std::stringstream xyz_string;

    // Loop over each row in the DataFrame.
    for (size_t row = 0; row < df.size(); row++) {
        // Calculate the number of vertices per iteration.
        size_t num_vertices = df[row].size() / 3;
        xyz_string << num_vertices << "\n";

        // Write the iteration number.
        xyz_string << "iteration: " << row << "\n";

        // Write the coordinates for each vertex.
        for (size_t i = 0; i < num_vertices; i++) {
            xyz_string << "  do";
            for (size_t j = 0; j < 3; j++) {
                double coord = df[row][i * 3 + j];
                xyz_string << format_spaces_decimals(coord);
            }
            xyz_string << "\n";
        }

        std::cout << "complete:" << row << " out of " << df.size() << std::endl;
    }

    return xyz_string.str();
}

// Main function.
int main() {
    // Read the input DataFrame from a CSV file.
    std::ifstream infile("meshpoint.csv");
    std::vector<std::vector<double>> df;
    std::string line;
    while (std::getline(infile, line)) {
        std::vector<double> row;
        std::istringstream iss(line);
        double val;
        while (iss >> val) {
            row.push_back(val);
            if (iss.peek() == ',')
                iss.ignore();
        }
        df.push_back(row);
    }

    // Convert the DataFrame to an XYZ string.
    std::string xyz_str = to_xyz_string(df);

    // Write the XYZ string to an output file.
    std::ofstream outfile("out_trajectory.xyz");
    outfile << xyz_str;
    outfile.close();

    return 0;
}
*/
