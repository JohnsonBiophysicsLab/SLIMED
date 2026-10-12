/**
 * @file Elastic_box.cpp
 * @brief The elastic periodic box (Param::elasticBox*): the strain of the
 * periodic axes as a Brownian coordinate driven by the frame tension.
 *
 * A periodic box of fixed size fixes the projected area, and whatever
 * tension that projected area implies is what the membrane carries: a sheet
 * that fluctuates, wraps a bead or takes an insertion has to find its excess
 * area by stretching, since there is no reservoir to draw on. The elastic
 * box is that reservoir. The strain lambda of the periodic axes becomes a
 * coordinate with the effective potential
 *
 *     U(lambda) = E(lambda x) - tau A_p(lambda) ,
 *
 * tau the prescribed frame tension and A_p the frame's area, so that its
 * generalized force is -dU/dlambda = W + tau dA_p/dlambda with
 *
 *     W = sum over vertices of x_i F_i          (Mesh::accumulate_virial())
 *
 * the virial of the raw nodal forces along that axis -- exactly
 * -dE/dlambda for an affine stretch of the whole control net, because Loop
 * subdivision is affine invariant and the images and ghosts stretch with
 * their sources. Per unit frame area the force is tau - tau_measured with
 * tau_measured = -W / A_p, the tension the membrane exerts on its frame,
 * so the box takes the same overdamped step the vertices do,
 *
 *     d lambda = mu (-dU/dlambda) dt + sqrt(2 mu kT dt) xi ,
 *
 * and the whole net is rescaled by 1 + d lambda about the origin. The
 * vertices step in absolute coordinates and the box in the strain, which
 * together sample exp(-U/kT) in the vertices' scaled coordinates. There is
 * deliberately no N kT ln(L) term: the in-plane vertex coordinates
 * parameterize a continuum surface, they are not the positions of
 * particles, and the term would push a fictitious vertex gas against the
 * membrane's real tension. The entropy of the height fluctuations, which is
 * real, needs no such term -- z is not rescaled, so the mean virial carries
 * it (see docs/elastic_box.md, section 3).
 */
#include "Dynamics.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>

#include "Counter_rng.hpp"

using slimed::splitmix64;
using slimed::standard_normal;

// ---------------------------------------------------------------------------
// DynamicMesh: what wraps, how long the period is, and the rescale
// ---------------------------------------------------------------------------

void DynamicMesh::detect_periodic_box()
{
    boxPeriodicAxis = {{0, 0, 0}};
    boxLength = {{0.0, 0.0, 0.0}};
    boxScale = {{1.0, 1.0, 1.0}};
    boxWidth = 0.0;

    switch (param.boundaryCondition)
    {
    case BoundaryType::Periodic:
        // The tile the grid lays out: three image rings and a duplicate ring
        // on each side. An imported mesh under this mode has no grid, and
        // then nothing here can say what wraps.
        if (param.nFaceX > 6 && param.nFaceY > 6)
        {
            boxPeriodicAxis[0] = 1;
            boxPeriodicAxis[1] = 1;
            boxLength[0] = (param.nFaceX - 6) * param.dFaceX;
            boxLength[1] = (param.nFaceY - 6) * param.dFaceY;
        }
        break;
    case BoundaryType::Mixed:
        // An axis wraps when some image is offset along it, and the period
        // is the shortest such offset: the generated band and any proper
        // tiling put every image exactly one period from its source.
        for (int v : periodicImageVertices)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                const double offset = std::abs(vertices[v].mirrorOffset[axis]);
                if (offset > 1.0e-9)
                {
                    boxPeriodicAxis[axis] = 1;
                    boxLength[axis] =
                        (boxLength[axis] == 0.0) ? offset : std::min(boxLength[axis], offset);
                }
            }
        }
        break;
    case BoundaryType::Fixed:
    case BoundaryType::Free:
        break;
    }
}

int DynamicMesh::box_periodic_axis_count() const
{
    int count = 0;
    for (int axis = 0; axis < 3; axis++)
    {
        count += boxPeriodicAxis[axis] ? 1 : 0;
    }
    return count;
}

double DynamicMesh::box_projected_area() const
{
    const int nAxes = box_periodic_axis_count();
    if (nAxes == 0)
    {
        return 0.0;
    }
    double area = 1.0;
    for (int axis = 0; axis < 3; axis++)
    {
        if (boxPeriodicAxis[axis])
        {
            area *= boxLength[axis];
        }
    }
    return (nAxes == 1) ? area * boxWidth : area;
}

void DynamicMesh::apply_box_scale(const std::array<double, 3> &factor)
{
    bool identity = true;
    for (int axis = 0; axis < 3; axis++)
    {
        if (!(factor[axis] > 0.0) || !std::isfinite(factor[axis]))
        {
            throw std::invalid_argument(
                "[DynamicMesh::apply_box_scale] the factor along axis " + std::to_string(axis) +
                " is " + std::to_string(factor[axis]) + "; it must be a positive finite number");
        }
        if (factor[axis] != 1.0)
        {
            identity = false;
        }
    }
    if (identity)
    {
        return;
    }

    const int nVertices = static_cast<int>(vertices.size());
    const bool haveMat = !matMesh.empty();
    for (int v = 0; v < nVertices; v++)
    {
        Vertex &vertex = vertices[v];
        for (int axis = 0; axis < 3; axis++)
        {
            const double f = factor[axis];
            if (f == 1.0)
            {
                continue;
            }
            vertex.coord.set(axis, 0, vertex.coord(axis, 0) * f);
            // The previous and reference coordinates stretch with the frame:
            // a uniform stretch leaves the mesh quality where it was, and the
            // solid regularization measures quality against the reference.
            // Neither is allocated until the driver's initialization.
            if (!vertex.coordPrev.empty())
            {
                vertex.coordPrev.set(axis, 0, vertex.coordPrev(axis, 0) * f);
            }
            if (!vertex.coordRef.empty())
            {
                vertex.coordRef.set(axis, 0, vertex.coordRef(axis, 0) * f);
            }
            // An image at source + offset is at f source + f offset afterwards.
            vertex.mirrorOffset[axis] *= f;
            if (haveMat)
            {
                matMesh.set(v, axis, matMesh(v, axis) * f);
            }
        }
    }
    for (int axis = 0; axis < 3; axis++)
    {
        boxLength[axis] *= factor[axis];
        boxScale[axis] *= factor[axis];
    }
    // The solver bakes the offsets into its right-hand side.
    if (param.boundaryCondition == BoundaryType::Mixed)
    {
        surfaceSolver.refresh_image_offsets(*this);
    }
}

// ---------------------------------------------------------------------------
// DynamicModel: the box as a Brownian coordinate
// ---------------------------------------------------------------------------

namespace
{
const char *axis_name(int axis)
{
    return (axis == 0) ? "x" : (axis == 1) ? "y" : "z";
}
} // namespace

void DynamicModel::setup_elastic_box()
{
    boxMobility = 0.0;
    boxVirial = {{0.0, 0.0, 0.0}};
    boxTension = {{0.0, 0.0, 0.0}};
    boxTensionMean = 0.0;
    boxTensionByTerm = {{0.0, 0.0, 0.0, 0.0, 0.0}};
    boxInPlaneVertexCount = 0;
    boxTensionVertexGas = 0.0;
    boxTensionRunningSum = 0.0;
    boxSteps = 0;
    boxClippedSteps = 0;

    const Param &param = mesh.param;
    if (!param.elasticBoxEnabled)
    {
        return;
    }

    const int nAxes = mesh.box_periodic_axis_count();
    if (nAxes == 0)
    {
        throw std::runtime_error(
            "[DynamicModel::setup_elastic_box] elasticBoxEnabled = true, but the sheet wraps "
            "along no axis, so the box has nothing to stretch. The elastic box needs "
            "boundaryType = Periodic, or Mixed with at least one periodic axis (some image "
            "offset along it). See docs/elastic_box.md.");
    }
    if (!(param.elasticBoxMobility > 0.0) || !std::isfinite(param.elasticBoxMobility))
    {
        throw std::runtime_error(
            "[DynamicModel::setup_elastic_box] elasticBoxMobility must be a positive number; "
            "it is " + std::to_string(param.elasticBoxMobility) + ".");
    }
    if (!(param.elasticBoxMaxStrainPerStep > 0.0) || !std::isfinite(param.elasticBoxMaxStrainPerStep))
    {
        throw std::runtime_error(
            "[DynamicModel::setup_elastic_box] elasticBoxMaxStrainPerStep must be a positive "
            "number; it is " + std::to_string(param.elasticBoxMaxStrainPerStep) + ".");
    }
    const bool isotropic = (param.elasticBoxCoupling == "isotropic");
    if (!isotropic && param.elasticBoxCoupling != "anisotropic")
    {
        throw std::runtime_error(
            "[DynamicModel::setup_elastic_box] elasticBoxCoupling must be isotropic or "
            "anisotropic; it is \"" + param.elasticBoxCoupling + "\".");
    }

    if (nAxes == 1)
    {
        // The frame's width across the one periodic axis: the membrane's
        // area per unit period, as built. A flat strip's width; a tube's
        // circumference.
        mesh.calculate_element_area_volume();
        double area = 0.0;
        double volume = 0.0;
        mesh.sum_membrane_area_and_volume(area, volume);
        for (int axis = 0; axis < 3; axis++)
        {
            if (mesh.boxPeriodicAxis[axis])
            {
                mesh.boxWidth = area / mesh.boxLength[axis];
            }
        }
    }
    const double area = mesh.box_projected_area();
    if (!(area > 0.0) || !std::isfinite(area))
    {
        throw std::runtime_error(
            "[DynamicModel::setup_elastic_box] the frame area came out as " +
            std::to_string(area) + " nm^2; the periodic box could not be read off the mesh.");
    }

    // The vertices whose in-plane coordinates the step integrates -- the
    // same test next_step() applies -- are particles in the plane as far as
    // the box is concerned, and their entropy is part of its force.
    if (param.elasticBoxVertexGasTerm && param.inPlaneDynamicsEnabled)
    {
        const bool mixed = (param.boundaryCondition == BoundaryType::Mixed);
        const int nVertices = static_cast<int>(mesh.vertices.size());
        for (int i = 0; i < nVertices; i++)
        {
            const bool isSlaved = (mixed || !param.integratePeriodicDuplicates) &&
                                  !mesh.isSlavedPeriodic.empty() && mesh.isSlavedPeriodic[i] != 0;
            if (mesh.vertices[i].isGhost || isSlaved || mesh.vertices[i].is_fixed())
            {
                continue;
            }
            boxInPlaneVertexCount++;
        }
    }

    // d lambda / dt = m (tau - tau_measured): the strain's mobility is m over
    // the generalized force per unit tension mismatch, n A_p for the one
    // isotropic strain and A_p for each anisotropic one. Fixed now, from
    // the area as built, so that the noise below needs no drift correction
    // -- an Ito mobility that depended on lambda would.
    boxMobility = param.elasticBoxMobility / ((isotropic ? nAxes : 1) * area);
    refresh_box_tension();

    std::cout << "[DynamicModel] Elastic box: " << nAxes << " periodic "
              << (nAxes > 1 ? "axes" : "axis") << " (";
    bool first = true;
    for (int axis = 0; axis < 3; axis++)
    {
        if (!mesh.boxPeriodicAxis[axis])
        {
            continue;
        }
        std::cout << (first ? "" : ", ") << axis_name(axis) << " = " << mesh.boxLength[axis] << " nm";
        first = false;
    }
    if (nAxes == 1)
    {
        std::cout << ", width " << mesh.boxWidth << " nm";
    }
    std::cout << "), frame area " << area << " nm^2; target tension " << param.elasticBoxTension
              << " pN/nm, " << param.elasticBoxCoupling << " coupling, "
              << (param.elasticBoxThermal ? "thermal" : "deterministic") << ", strain mobility "
              << param.elasticBoxMobility << " nm/(pN us), strain cap "
              << param.elasticBoxMaxStrainPerStep << " per step." << std::endl;
    std::cout << "[DynamicModel]   Tension as built: " << boxTensionMean << " pN/nm (bending "
              << boxTensionByTerm[0] << ", area constraint " << boxTensionByTerm[1] << ", mesh term "
              << boxTensionByTerm[2] << ", bead " << boxTensionByTerm[3] << ", vertex gas "
              << boxTensionByTerm[4] << " from " << boxInPlaneVertexCount
              << " in-plane vertices)." << std::endl;
    if (param.uSurf > 0.0)
    {
        const int strainAxes = isotropic ? nAxes : 1;
        const double rate = param.elasticBoxMobility * strainAxes * param.uSurf;
        const double rateDt = rate * param.timeStep;
        std::cout << "[DynamicModel]   With usMembraneStretching = " << param.uSurf
                  << " pN/nm as the area modulus the strain relaxes at " << rate << " per us ("
                  << 1.0 / rate << " us, " << 1.0 / rateDt << " steps); rate x dt = " << rateDt
                  << (rateDt > 0.1 ? " -- WARNING: stiff for an explicit step; lower "
                                     "elasticBoxMobility or timeStep."
                                   : ".")
                  << std::endl;
        if (param.elasticBoxThermal)
        {
            const double strainStd =
                std::sqrt(param.KBT / (static_cast<double>(strainAxes) * strainAxes * area * param.uSurf));
            std::cout << "[DynamicModel]   The thermal strain then has a standard deviation of "
                      << strainStd << " (" << 100.0 * strainAxes * strainStd
                      << "% of the frame area)." << std::endl;
        }
    }
}

void DynamicModel::refresh_box_tension()
{
    const Param &param = mesh.param;
    const bool spring = param.edgeSpringEnabled;
    const bool bead = param.beadEnabled && !beadForce.empty();
    const double area = mesh.box_projected_area();
    const int nAxes = mesh.box_periodic_axis_count();

    double beadVirial[3] = {0.0, 0.0, 0.0};
    if (bead)
    {
        // The bead is not part of the frame: it stays where it is while the
        // membrane stretches past it, so its wall's virial is at the
        // vertices' positions, like every other term.
        const int nVertices = static_cast<int>(mesh.vertices.size());
        for (int v = 0; v < nVertices; v++)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                beadVirial[axis] += mesh.vertices[v].coord(axis, 0) * beadForce(v, axis);
            }
        }
    }

    // -N kT ln(A_p) in the effective potential: -N kT / A_p on every
    // periodic axis, the 2D ideal-gas pressure of the vertices in the plane.
    boxTensionVertexGas = (area > 0.0) ? -boxInPlaneVertexCount * param.KBT / area : 0.0;

    double sumByTerm[4] = {0.0, 0.0, 0.0, 0.0};
    double sumTension = 0.0;
    for (int axis = 0; axis < 3; axis++)
    {
        const double curvature = mesh.virialCurvature[axis];
        const double areaTerm = mesh.virialArea[axis];
        const double meshTerm = spring ? mesh.virialRegularization[axis] : 0.0;
        boxVirial[axis] = curvature + areaTerm + meshTerm + beadVirial[axis];
        boxTension[axis] = (area > 0.0) ? -boxVirial[axis] / area : 0.0;
        if (mesh.boxPeriodicAxis[axis])
        {
            boxTension[axis] += boxTensionVertexGas;
            sumTension += boxTension[axis];
            sumByTerm[0] += curvature;
            sumByTerm[1] += areaTerm;
            sumByTerm[2] += meshTerm;
            sumByTerm[3] += beadVirial[axis];
        }
    }
    boxTensionMean = (nAxes > 0) ? sumTension / nAxes : 0.0;
    for (int k = 0; k < 4; k++)
    {
        boxTensionByTerm[k] = (nAxes > 0 && area > 0.0) ? -sumByTerm[k] / (nAxes * area) : 0.0;
    }
    boxTensionByTerm[4] = boxTensionVertexGas;
}

void DynamicModel::update_elastic_box()
{
    const Param &param = mesh.param;
    if (!param.elasticBoxEnabled)
    {
        return;
    }
    refresh_box_tension();

    const double dt = param.timeStep;
    const double target = param.elasticBoxTension;
    const bool isotropic = (param.elasticBoxCoupling == "isotropic");
    const int nAxes = mesh.box_periodic_axis_count();
    const double area = mesh.box_projected_area();
    const double noiseScale =
        param.elasticBoxThermal ? std::sqrt(2.0 * boxMobility * param.KBT * dt) : 0.0;
    const double cap = param.elasticBoxMaxStrainPerStep;

    // One key per (run, step), as the vertices use, and a lane no vertex
    // has: the box draws as the vertex after the last one.
    const std::uint64_t stepKey = splitmix64(static_cast<std::uint64_t>(randomSeed)) +
                                  static_cast<std::uint64_t>(this->iteration);
    const std::uint64_t boxLane = static_cast<std::uint64_t>(mesh.vertices.size());

    bool clipped = false;
    bool nonFinite = false;
    auto strainIncrement = [&](double force, int lane) {
        double d = boxMobility * force * dt;
        if (noiseScale > 0.0)
        {
            d += noiseScale * standard_normal(stepKey, boxLane, static_cast<std::uint64_t>(lane));
        }
        if (!std::isfinite(d))
        {
            nonFinite = true;
            return 0.0;
        }
        if (d > cap)
        {
            clipped = true;
            return cap;
        }
        if (d < -cap)
        {
            clipped = true;
            return -cap;
        }
        return d;
    };

    std::array<double, 3> factor = {{1.0, 1.0, 1.0}};
    if (isotropic)
    {
        // -dU/dlambda for A_p proportional to lambda^n.
        const double force = nAxes * area * (target - boxTensionMean);
        const double d = strainIncrement(force, 0);
        for (int axis = 0; axis < 3; axis++)
        {
            if (mesh.boxPeriodicAxis[axis])
            {
                factor[axis] = 1.0 + d;
            }
        }
    }
    else
    {
        for (int axis = 0; axis < 3; axis++)
        {
            if (mesh.boxPeriodicAxis[axis])
            {
                const double force = area * (target - boxTension[axis]);
                factor[axis] = 1.0 + strainIncrement(force, axis);
            }
        }
    }

    if (nonFinite)
    {
        // The run will stop on its non-finite energy; leave the box alone
        // rather than fold NaN into every coordinate first.
        std::cerr << "[DynamicModel::update_elastic_box] the box force is not finite at iteration "
                  << this->iteration << "; the box was not moved." << std::endl;
        return;
    }
    if (clipped)
    {
        boxClippedSteps++;
    }
    boxTensionRunningSum += boxTensionMean;
    boxSteps++;
    mesh.apply_box_scale(factor);
}

void DynamicModel::append_box_log(const std::string &filepath, int iteration)
{
    std::ofstream out(filepath, std::ios::app);
    if (!out)
    {
        throw std::runtime_error("[DynamicModel::append_box_log] cannot open " + filepath);
    }
    if (!boxLogHeaderWritten)
    {
        out << "iteration,time_us,scaleX,scaleY,scaleZ,boxLengthX_nm,boxLengthY_nm,boxLengthZ_nm,"
               "frameArea_nm2,membraneArea_nm2,tension_pN_nm,tensionX_pN_nm,tensionY_pN_nm,"
               "tensionZ_pN_nm,tensionBending_pN_nm,tensionAreaConstraint_pN_nm,"
               "tensionMeshTerm_pN_nm,tensionBead_pN_nm,tensionVertexGas_pN_nm,clippedSteps\n";
        boxLogHeaderWritten = true;
    }
    out << iteration << ',' << std::setprecision(12) << iteration * mesh.param.timeStep << ','
        << mesh.boxScale[0] << ',' << mesh.boxScale[1] << ',' << mesh.boxScale[2] << ','
        << mesh.boxLength[0] << ',' << mesh.boxLength[1] << ',' << mesh.boxLength[2] << ','
        << mesh.box_projected_area() << ',' << mesh.param.area << ',' << boxTensionMean << ','
        << boxTension[0] << ',' << boxTension[1] << ',' << boxTension[2] << ','
        << boxTensionByTerm[0] << ',' << boxTensionByTerm[1] << ',' << boxTensionByTerm[2] << ','
        << boxTensionByTerm[3] << ',' << boxTensionByTerm[4] << ',' << boxClippedSteps << '\n';
}
