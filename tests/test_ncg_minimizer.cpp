#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

#include "Parameters.hpp"
#include "mesh/Mesh.hpp"
#include "model/Model.hpp"
#include "model/Record.hpp"

/**
 * @file test_ncg_minimizer.cpp
 * @brief The energy minimizer's NCG switch, its Wolfe test and its stop.
 *
 * Three bugs kept the run_flat minimizer from doing what it says:
 *  - OptimizationAlgorithm::usingNCG defaulted to false and Param::usingNCG
 *    was never copied into it, so NCG never ran;
 *  - the strong-Wolfe curvature term of the NCG line search was taken from
 *    the force at the start point instead of the trial point, so it equalled
 *    the reference term and the test only passed for c2 >= 1;
 *  - the line search returned -1 on failure without setting Model::stepSize,
 *    which is what run_flat tested, so its "could not find an efficient step
 *    size" stop never fired.
 */

namespace
{

/// A 60 nm periodic sheet with no scaffolding.
void configure_sheet(Param &param)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Periodic;
    param.sideX = 60.0;
    param.sideY = 60.0;
    param.lFace = 5.0;
    param.kCurv = 83.4;
    param.uSurf = 250.0;
    param.uVol = 0.0;
    param.isGlobalConstraint = true;
    param.isEnergyHarmonicBondIncluded = false;
}

/// Build the flat sheet as its own reference, lift the centre vertex out of
/// plane so bending and stretching both pull on it, and take that as the
/// start point of a line search: energy, force and the "previous" copies.
void build_lifted_sheet(Mesh &mesh)
{
    mesh.setup_flat();
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();

    const int centre = (mesh.param.nFaceX + 1) * (mesh.param.nFaceY / 2) + mesh.param.nFaceX / 2;
    Vertex &vertex = mesh.vertices[centre];
    vertex.coord.set(2, 0, vertex.coord(2, 0) + 1.0);

    mesh.Compute_Energy_And_Force();
    mesh.update_previous_coord_for_vertex();
    mesh.update_previous_force_for_vertex();
}

/// -sum_i F_i . d_i over the vertices, F taken from force or forcePrev: the
/// derivative of the energy along the search direction d.
double slope_along(const Mesh &mesh, const std::vector<Force> &direction, bool previous)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        const Matrix &force = previous ? mesh.vertices[i].forcePrev.forceTotal
                                       : mesh.vertices[i].force.forceTotal;
        for (int axis = 0; axis < 3; axis++)
        {
            sum -= force.get(axis, 0) * direction[i].forceTotal.get(axis, 0);
        }
    }
    return sum;
}

} // namespace

// ---------------------------------------------------------------------------
// Bug 2: the curvature condition is evaluated at the trial point
// ---------------------------------------------------------------------------

/**
 * With c2 just below 1 the strong-Wolfe curvature test asks only that the
 * slope at the trial point be smaller in magnitude than at the start, which
 * every short enough step along a descent direction satisfies. Taken from
 * the start-point force, the two slopes were equal and the search always ran
 * out of step and fell back to the simple line search.
 */
TEST(NcgMinimizerTest, TheCurvatureConditionUsesTheForceAtTheTrialPoint)
{
    Param param;
    configure_sheet(param);
    Mesh mesh(param);
    build_lifted_sheet(mesh);

    Record record(2);
    record.add(param.area, param.energy, mesh.calculate_mean_force());
    Model model(mesh, record);
    model.oa.usingNCG = true;
    model.oa.c2 = 0.99;
    model.determine_trial_step_size();
    ASSERT_GT(model.oa.trialStepSize, 0.0);

    const double energyStart = param.energy.energyTotal;
    const double slopeStart = slope_along(mesh, model.ncgDirection0, true);
    ASSERT_LT(slopeStart, 0.0) << "the force must be a descent direction";

    const double step = model.linear_search_for_stepsize_to_minimize_energy();

    EXPECT_FALSE(model.oa.isNCGstuck) << "the NCG branch fell back to the simple line search";
    ASSERT_GT(step, 0.0);
    EXPECT_EQ(model.stepSize, step);

    // The search leaves the mesh at the accepted point: both Wolfe
    // conditions hold there, with the slope from the force at that point.
    const double tolerance = 1e-12 * std::abs(energyStart);
    EXPECT_LE(param.energy.energyTotal,
              energyStart + model.oa.c1 * step * slopeStart + tolerance);
    const double slopeAccepted = slope_along(mesh, model.ncgDirection0, false);
    EXPECT_LE(std::abs(slopeAccepted), model.oa.c2 * std::abs(slopeStart));
    EXPECT_GT(std::abs(slopeStart - slopeAccepted), 1e-9 * std::abs(slopeStart))
        << "the step moved the mesh too little to tell the two slopes apart";
}

/// The same search with c2 = 0 cannot meet the curvature condition, so it
/// falls back: the fix must not make the NCG branch accept everything.
TEST(NcgMinimizerTest, AnUnsatisfiableCurvatureConditionStillFallsBack)
{
    Param param;
    configure_sheet(param);
    Mesh mesh(param);
    build_lifted_sheet(mesh);

    Record record(2);
    record.add(param.area, param.energy, mesh.calculate_mean_force());
    Model model(mesh, record);
    model.oa.usingNCG = true;
    model.oa.c2 = 0.0;
    model.determine_trial_step_size();
    ASSERT_GT(model.oa.trialStepSize, 0.0);

    const double step = model.linear_search_for_stepsize_to_minimize_energy();

    EXPECT_TRUE(model.oa.isNCGstuck);
    EXPECT_GT(step, 0.0) << "the simple line search should find a decrease along the force";
}
