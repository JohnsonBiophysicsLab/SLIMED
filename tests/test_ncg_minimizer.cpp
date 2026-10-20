#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

#include "Parameters.hpp"
#include "Run_simulation.hpp"
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
 *
 * Fixing those alone made run_flat stop early, for four more reasons, also
 * pinned here:
 *  - the line search shrank the step by 0.8 before its first trial while the
 *    next trial step was 1.1 times the accepted one, so the step decayed
 *    below the threshold within ~250 iterations and the search failed;
 *  - the Fletcher-Reeves factor was computed from forcePrev after run_flat
 *    had copied the current force into it, so it was 1 at every step;
 *  - nothing restarted the direction from the force after a fallback or
 *    when it pointed uphill;
 *  - pure backtracking with c2 = 0.001 almost never met the strong Wolfe
 *    conditions, so the NCG branch now expands and brackets the step.
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

/// The flat sheet as its own reference with its own area as the relaxed
/// area: a minimum of the energy, with no force to follow.
void build_relaxed_sheet(Mesh &mesh)
{
    mesh.setup_flat();
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.Compute_Energy_And_Force();
    mesh.update_previous_coord_for_vertex();
    mesh.update_previous_force_for_vertex();
}

/// sum_i |F_i|^2 over the vertices.
double force_norm_squared(const Mesh &mesh)
{
    double sum = 0.0;
    for (const Vertex &vertex : mesh.vertices)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            sum += vertex.force.forceTotal.get(axis, 0) * vertex.force.forceTotal.get(axis, 0);
        }
    }
    return sum;
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

/// A fresh temporary directory that is the working directory for the
/// lifetime of the object: run_flat reads ./input.params and writes
/// fixed-name outputs into the current directory.
class ScratchDirectory
{
public:
    ScratchDirectory()
    {
        char cwd[4096];
        if (getcwd(cwd, sizeof cwd) == nullptr)
        {
            throw std::runtime_error("getcwd failed");
        }
        previous_ = cwd;
        const char *tmp = std::getenv("TMPDIR");
        std::string pattern = std::string(tmp ? tmp : "/tmp") + "/slimed_ncg_minimizer_XXXXXX";
        std::vector<char> buffer(pattern.begin(), pattern.end());
        buffer.push_back('\0');
        if (mkdtemp(buffer.data()) == nullptr)
        {
            throw std::runtime_error("mkdtemp failed");
        }
        path_ = buffer.data();
        if (chdir(path_.c_str()) != 0)
        {
            throw std::runtime_error("chdir failed");
        }
    }

    ~ScratchDirectory()
    {
        if (chdir(previous_.c_str()) == 0)
        {
            std::system(("rm -rf '" + path_ + "'").c_str());
        }
    }

private:
    std::string previous_;
    std::string path_;
};

} // namespace

// ---------------------------------------------------------------------------
// Bug 1: the NCG switch comes from the parameters
// ---------------------------------------------------------------------------

TEST(NcgMinimizerTest, TheModelTakesUsingNcgFromTheParameters)
{
    EXPECT_TRUE(Param().usingNCG) << "the documented default is NCG on";

    for (const bool configured : {true, false})
    {
        Param param;
        param.VERBOSE_MODE = false;
        param.usingNCG = configured;
        Mesh mesh(param);
        Record record(1);
        Model model(mesh, record);
        EXPECT_EQ(model.oa.usingNCG, configured);
    }
}

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

/// The same search with c2 = 0 cannot meet the curvature condition: the fix
/// must not make the NCG branch accept everything. It may still take a step
/// with sufficient decrease, but marks the search stuck so that the next
/// direction restarts from the force.
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

// ---------------------------------------------------------------------------
// Bug 3: a failed search is visible through Model::stepSize
// ---------------------------------------------------------------------------

/// At a minimum no step lowers the energy, and the failure must show in the
/// member run_flat reads as well as in the return value.
TEST(NcgMinimizerTest, AFailedLineSearchSetsTheStepSizeToMinusOne)
{
    Param param;
    configure_sheet(param);
    Mesh mesh(param);
    build_relaxed_sheet(mesh);
    ASSERT_EQ(mesh.get_max_force_magnitude(), 0.0) << "the relaxed sheet is not a minimum";

    Record record(2);
    record.add(param.area, param.energy, mesh.calculate_mean_force());
    Model model(mesh, record);
    model.oa.usingNCG = false;
    model.oa.stepThreshold = 1e-8; // fail in a few dozen trials, not ~150
    model.oa.trialStepSize = 0.01;

    testing::internal::CaptureStdout();
    const double step = model.linear_search_for_stepsize_to_minimize_energy();
    const std::string out = testing::internal::GetCapturedStdout();

    EXPECT_EQ(step, -1.0);
    EXPECT_EQ(model.stepSize, -1.0);
    EXPECT_NE(out.find("cannot find an efficient small stepSize"), std::string::npos) << out;
}

// ---------------------------------------------------------------------------
// The step size: tried first, and looked for again before giving up
// ---------------------------------------------------------------------------

/// A trial step that lowers the energy is the step taken. The search used to
/// shrink it by 0.8 before trying it, so with the next trial at 1.1 times the
/// accepted step the step shrank every iteration.
TEST(NcgMinimizerTest, TheSimpleSearchTriesTheTrialStepFirst)
{
    Param param;
    configure_sheet(param);
    Mesh mesh(param);
    build_lifted_sheet(mesh);

    Record record(2);
    record.add(param.area, param.energy, mesh.calculate_mean_force());
    Model model(mesh, record);
    model.oa.usingNCG = false;
    model.determine_trial_step_size();
    model.oa.trialStepSize *= 0.01; // short enough to lower the energy
    const double energyStart = param.energy.energyTotal;

    const double step = model.linear_search_for_stepsize_to_minimize_energy();

    EXPECT_EQ(step, model.oa.trialStepSize);
    EXPECT_LT(param.energy.energyTotal, energyStart);
}

/// When no step along the search direction lowers the energy (here it points
/// uphill), the search tries once more along the force from the rule-of-thumb
/// step, with the regularization reference refreshed, before it reports a
/// failure: a failure means no descent, not a direction or a step that
/// happened to be bad.
TEST(NcgMinimizerTest, AFailedSearchRetriesAlongTheForceBeforeGivingUp)
{
    Param param;
    configure_sheet(param);
    Mesh mesh(param);
    build_lifted_sheet(mesh);

    Record record(2);
    record.add(param.area, param.energy, mesh.calculate_mean_force());
    Model model(mesh, record);
    model.oa.usingNCG = false;
    model.oa.stepThreshold = 1e-8; // fail the first pass in a few dozen trials
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        model.ncgDirection0[i].forceTotal = -1.0 * mesh.vertices[i].forcePrev.forceTotal;
    }
    model.determine_trial_step_size();
    model.oa.trialStepSize *= 1e-3; // the trial step of a collapsed run
    const double energyStart = param.energy.energyTotal;

    testing::internal::CaptureStdout();
    const double step = model.linear_search_for_stepsize_to_minimize_energy();
    const std::string out = testing::internal::GetCapturedStdout();

    EXPECT_EQ(out.find("cannot find an efficient small stepSize"), std::string::npos) << out;
    ASSERT_GT(step, 0.0);
    EXPECT_GT(step, model.oa.trialStepSize) << "the retry should start from the rule of thumb";
    EXPECT_LT(param.energy.energyTotal, energyStart);
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            ASSERT_EQ(model.ncgDirection0[i].forceTotal.get(axis, 0),
                      mesh.vertices[i].forcePrev.forceTotal.get(axis, 0));
            ASSERT_EQ(mesh.vertices[i].coordRef.get(axis, 0),
                      mesh.vertices[i].coordPrev.get(axis, 0));
        }
    }
}

/// A trial step far too short for the curvature condition: pure
/// backtracking could only shrink it further and fall back. The strong-Wolfe
/// search doubles it until a trial overshoots, then accepts a step that meets
/// both conditions.
TEST(NcgMinimizerTest, TheStrongWolfeSearchExpandsAShortTrialStep)
{
    Param param;
    configure_sheet(param);
    // Clamped: on the legacy periodic sheet the ghosts move with their
    // sources but carry none of the force, so -F . d is not the slope of the
    // energy along d and the slope at the line minimum is about half the
    // initial one. No strong-Wolfe step exists there.
    param.boundaryCondition = BoundaryType::Fixed;
    Mesh mesh(param);
    build_lifted_sheet(mesh);

    Record record(2);
    record.add(param.area, param.energy, mesh.calculate_mean_force());
    Model model(mesh, record);
    model.oa.usingNCG = true;
    model.determine_trial_step_size();
    model.oa.trialStepSize *= 1e-3;
    const double energyStart = param.energy.energyTotal;
    const double slopeStart = slope_along(mesh, model.ncgDirection0, true);
    ASSERT_LT(slopeStart, 0.0);

    const double step = model.linear_search_for_stepsize_to_minimize_energy();

    EXPECT_FALSE(model.oa.isNCGstuck);
    EXPECT_GT(step, 2.0 * model.oa.trialStepSize);
    const double tolerance = 1e-12 * std::abs(energyStart);
    EXPECT_LE(param.energy.energyTotal,
              energyStart + model.oa.c1 * step * slopeStart + tolerance);
    EXPECT_LE(std::abs(slope_along(mesh, model.ncgDirection0, false)),
              model.oa.c2 * std::abs(slopeStart));
}

// ---------------------------------------------------------------------------
// The conjugate direction
// ---------------------------------------------------------------------------

/// The Fletcher-Reeves factor is |F_k|^2 / |F_(k-1)|^2 with F_(k-1) the force
/// at the previous update. run_flat copies the current force into forcePrev
/// before it updates the direction (as build_lifted_sheet leaves it here), so
/// taken from forcePrev the factor was 1 at every step.
TEST(NcgMinimizerTest, TheFletcherReevesFactorUsesTheForceOfThePreviousUpdate)
{
    Param param;
    configure_sheet(param);
    Mesh mesh(param);
    build_lifted_sheet(mesh);

    Record record(2);
    Model model(mesh, record);
    model.oa.usingNCG = true;
    model.oa.isNCGstuck = false;
    const double normSquared = force_norm_squared(mesh);
    ASSERT_GT(normSquared, 0.0);
    model.ncgForceNormSquaredPrev = 4.0 * normSquared; // the force has halved
    std::vector<Matrix> previous;
    for (const Force &direction : model.ncgDirection0)
    {
        previous.push_back(direction.forceTotal);
    }

    model.update_ncg_direction();

    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            const double expected = mesh.vertices[i].force.forceTotal.get(axis, 0) +
                                    0.25 * previous[i].get(axis, 0);
            ASSERT_NEAR(model.ncgDirection0[i].forceTotal.get(axis, 0), expected,
                        1e-12 * (std::abs(expected) + 1.0));
        }
    }
    EXPECT_NEAR(model.ncgForceNormSquaredPrev, normSquared, 1e-12 * normSquared);
}

/// After a search that fell back, the previous direction is not a conjugate
/// direction to build on: the next one is the force.
TEST(NcgMinimizerTest, TheDirectionRestartsFromTheForceAfterAFallback)
{
    Param param;
    configure_sheet(param);
    Mesh mesh(param);
    build_lifted_sheet(mesh);

    Record record(2);
    Model model(mesh, record);
    model.oa.usingNCG = true;
    model.oa.isNCGstuck = true;
    model.ncgForceNormSquaredPrev = force_norm_squared(mesh);

    model.update_ncg_direction();

    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            ASSERT_EQ(model.ncgDirection0[i].forceTotal.get(axis, 0),
                      mesh.vertices[i].force.forceTotal.get(axis, 0));
        }
    }
}

/// A Fletcher-Reeves direction that points uphill (sum F . d <= 0) is
/// replaced by the force; searched along, it could only fail.
TEST(NcgMinimizerTest, AnUphillConjugateDirectionRestartsFromTheForce)
{
    Param param;
    configure_sheet(param);
    Mesh mesh(param);
    build_lifted_sheet(mesh);

    Record record(2);
    Model model(mesh, record);
    model.oa.usingNCG = true;
    model.oa.isNCGstuck = false;
    model.ncgForceNormSquaredPrev = force_norm_squared(mesh); // factor 1
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        // F + (-10 F) = -9 F
        model.ncgDirection0[i].forceTotal = -10.0 * mesh.vertices[i].force.forceTotal;
    }

    model.update_ncg_direction();

    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            ASSERT_EQ(model.ncgDirection0[i].forceTotal.get(axis, 0),
                      mesh.vertices[i].force.forceTotal.get(axis, 0));
        }
    }
}

/**
 * run_flat stops at the first failed line search instead of carrying on.
 *
 * The failure here is genuine: on the 40 nm clamped sheet with the smooth
 * regularization (no adaptive cases, whose switching makes the energy jump)
 * NCG converges to the minimum well within 400 iterations, and there no step
 * lowers the energy. Before the fix run_flat went on to maxIterations (the
 * "Stopped" exit never fired); before the step-size and direction fixes it
 * stopped far from the minimum, when the step had decayed or NCG had been
 * turned off.
 */
TEST(NcgMinimizerTest, RunFlatStopsAtTheFirstFailedLineSearch)
{
    ScratchDirectory scratch;
    {
        std::ofstream params("input.params");
        params << "maxIterations = 400\n"
               << "VERBOSE_MODE = false\n"
               << "boundaryType = Fixed\n"
               << "lFace = 5.0\n"
               << "sideX = 40.0\n"
               << "sideY = 40.0\n"
               << "setRelaxAreaToDefault = false\n"
               << "relaxArea = 900.0\n"
               << "usingRpi = false\n"
               << "isEnergyHarmonicBondIncluded = false\n"
               << "thermalFluctuationEnabled = false\n"
               << "meshpointOutput = false\n"
               << "checkpointOutputInterval = 0\n";
    }

    testing::internal::CaptureStdout();
    run_flat("./input");
    const std::string out = testing::internal::GetCapturedStdout();

    const std::size_t failure = out.find("cannot find an efficient small stepSize");
    ASSERT_NE(failure, std::string::npos)
        << "no line search failed in 400 iterations; the sheet did not converge";
    const std::size_t stop = out.find("could not find an efficient step size. Stopped.", failure);
    ASSERT_NE(stop, std::string::npos) << "run_flat carried on after a failed line search";

    // Nothing was recorded between the failure and the stop.
    EXPECT_EQ(out.find("[Record::print_iteration()]", failure), std::string::npos);

    // It stopped early, at the failing step.
    std::istringstream stopLine(out.substr(out.rfind("Step: ", stop)));
    std::string label;
    int stoppedAt = -1;
    stopLine >> label >> stoppedAt;
    EXPECT_GT(stoppedAt, 0);
    EXPECT_LT(stoppedAt, 399);

    // The run stopped at the minimum, and NCG was not turned off on the way
    // there. Once the energy is down to round-off (1e-25 from 1e3 here), the
    // searches stall and NCG may switch itself off just before the stop;
    // which iteration that is depends on the last bits of the force sums.
    auto energy_of_record = [&out](std::size_t at) {
        const std::size_t key = out.find("energy= ", at);
        return std::stod(out.substr(key + 8));
    };
    const double energyFirst = energy_of_record(out.find("[Record::print_iteration()]"));
    const double energyLast = energy_of_record(out.rfind("[Record::print_iteration()]", failure));
    EXPECT_LT(energyLast, 1e-6 * energyFirst);
    const std::size_t disabled = out.find("Disabled NCG");
    if (disabled != std::string::npos)
    {
        EXPECT_LT(energy_of_record(out.rfind("[Record::print_iteration()]", disabled)), 1e-12 * energyFirst)
            << "NCG was turned off before the sheet converged";
    }
}
