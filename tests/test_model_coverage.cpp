#include "test_model_coverage.hpp"

namespace
{

/// A 25 nm x 25 nm free-edged sheet at lFace = 5: a handful of faces, cheap to evaluate.
void configure_free_sheet(Param &param)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Free;
    param.sideX = 25.0;
    param.sideY = 25.0;
    param.lFace = 5.0;
    param.kCurv = 83.4;
    param.uSurf = 250.0;
    param.uVol = 0.0;
    param.usingRpi = false;
    param.isGlobalConstraint = true;
    param.isEnergyHarmonicBondIncluded = false;
    param.randomSeed = 7u;
}

/// Flat sheet whose target area and reference geometry are the flat sheet.
void prepare_sheet(Mesh &mesh)
{
    mesh.setup_flat();
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.Compute_Energy_And_Force();
    mesh.update_previous_force_for_vertex();
}

/// The independent vertex closest to the middle of the sheet.
int central_mobile_vertex(const Mesh &mesh)
{
    const double cx = 0.5 * mesh.param.sideX;
    const double cy = 0.5 * mesh.param.sideY;
    int best = -1;
    double bestDistance = 1e300;
    for (int i = 0; i < static_cast<int>(mesh.vertices.size()); i++)
    {
        if (!mesh.is_independent_vertex(i))
        {
            continue;
        }
        const double dx = mesh.vertices[i].coord.get(0, 0) - cx;
        const double dy = mesh.vertices[i].coord.get(1, 0) - cy;
        const double d = dx * dx + dy * dy;
        if (d < bestDistance)
        {
            bestDistance = d;
            best = i;
        }
    }
    return best;
}

/// Lift one interior vertex out of plane and make that the line search's start point.
void lift_centre_and_rebase(Mesh &mesh, double dz)
{
    const int centre = central_mobile_vertex(mesh);
    ASSERT_GE(centre, 0);
    Vertex &vertex = mesh.vertices[centre];
    vertex.coord.set(2, 0, vertex.coord.get(2, 0) + dz);
    mesh.calculate_element_area_volume();
    mesh.Compute_Energy_And_Force();
    mesh.update_previous_coord_for_vertex();
    mesh.update_previous_force_for_vertex();
}

void set_vector(Matrix &m, double x, double y, double z)
{
    m.set(0, 0, x);
    m.set(1, 0, y);
    m.set(2, 0, z);
}

Energy energy_with_total(double total)
{
    Energy energy;
    energy.energyTotal = total;
    return energy;
}

double sum_force_dot_direction(const Mesh &mesh, const std::vector<Force> &direction)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            sum += mesh.vertices[i].forcePrev.forceTotal.get(axis, 0) *
                   direction[i].forceTotal.get(axis, 0);
        }
    }
    return sum;
}

} // namespace

// ---------------------------------------------------------------------------
// Record
// ---------------------------------------------------------------------------

TEST(ModelCoverageRecordTest, ConstructorsStartEmptyAndPreallocate)
{
    Record empty;
    EXPECT_TRUE(empty.areaTotal.empty());
    EXPECT_TRUE(empty.energyVec.empty());
    EXPECT_TRUE(empty.meanForce.empty());

    Record reserved(64);
    EXPECT_TRUE(reserved.areaTotal.empty());
    EXPECT_GE(reserved.areaTotal.capacity(), 64u);
    EXPECT_GE(reserved.energyVec.capacity(), 64u);
    EXPECT_GE(reserved.meanForce.capacity(), 64u);
}

TEST(ModelCoverageRecordTest, AddAppendsAndSetOverwritesOneIteration)
{
    Record record;
    record.add(10.0, energy_with_total(1.5), 0.25);
    record.add(20.0, energy_with_total(2.5), 0.5);
    ASSERT_EQ(record.areaTotal.size(), 2u);
    ASSERT_EQ(record.energyVec.size(), 2u);
    ASSERT_EQ(record.meanForce.size(), 2u);
    EXPECT_DOUBLE_EQ(record.areaTotal[1], 20.0);
    EXPECT_DOUBLE_EQ(record.energyVec[1].energyTotal, 2.5);
    EXPECT_DOUBLE_EQ(record.meanForce[1], 0.5);

    Energy replacement = energy_with_total(-3.0);
    replacement.energyHarmonicBond = 4.0;
    record.set(0, 11.0, replacement, 0.125);
    EXPECT_EQ(record.areaTotal.size(), 2u);
    EXPECT_DOUBLE_EQ(record.areaTotal[0], 11.0);
    EXPECT_DOUBLE_EQ(record.energyVec[0].energyTotal, -3.0);
    EXPECT_DOUBLE_EQ(record.energyVec[0].energyHarmonicBond, 4.0);
    EXPECT_DOUBLE_EQ(record.meanForce[0], 0.125);
    // The other iteration is untouched.
    EXPECT_DOUBLE_EQ(record.areaTotal[1], 20.0);
    EXPECT_DOUBLE_EQ(record.energyVec[1].energyTotal, 2.5);
}

TEST(ModelCoverageRecordTest, PrintIterationReportsThatIterationsValues)
{
    Record record;
    record.add(1.0, energy_with_total(9.0), 9.0);
    Energy energy = energy_with_total(3.5);
    energy.energyHarmonicBond = 1.25;
    energy.energyGagScaffolding = 0.75;
    energy.energyIdealizedProteinLattice = 0.5;
    record.add(42.0, energy, 0.0625);

    testing::internal::CaptureStdout();
    record.print_iteration(1);
    const std::string out = testing::internal::GetCapturedStdout();

    EXPECT_NE(out.find("[Record::print_iteration()] Step: 1."), std::string::npos) << out;
    EXPECT_NE(out.find("energy= 3.5"), std::string::npos) << out;
    EXPECT_NE(out.find("eHBond= 1.25"), std::string::npos) << out;
    EXPECT_NE(out.find("eGag= 0.75"), std::string::npos) << out;
    EXPECT_NE(out.find("eIdeal= 0.5"), std::string::npos) << out;
    EXPECT_NE(out.find("meanF= 0.0625"), std::string::npos) << out;
    EXPECT_NE(out.find("area= 42"), std::string::npos) << out;
    EXPECT_EQ(out.find("energy= 9"), std::string::npos) << "printed the wrong iteration: " << out;
}

// ---------------------------------------------------------------------------
// OptimizationAlgorithm
// ---------------------------------------------------------------------------

TEST(ModelCoverageOptimizationAlgorithmTest, CustomConstructorSetsTheLineSearchConstants)
{
    OptimizationAlgorithm defaults;
    OptimizationAlgorithm oa(17, 1e-3, 0.5, 1e-9);
    EXPECT_EQ(oa.trialIterationInterval, 17);
    EXPECT_DOUBLE_EQ(oa.c1, 1e-3);
    EXPECT_DOUBLE_EQ(oa.c2, 0.5);
    EXPECT_DOUBLE_EQ(oa.stepThreshold, 1e-9);
    // Everything else keeps its default.
    EXPECT_EQ(oa.usingNCG, defaults.usingNCG);
    EXPECT_EQ(oa.nConsecutiveNcgStuckThreshold, defaults.nConsecutiveNcgStuckThreshold);
    EXPECT_DOUBLE_EQ(oa.energyDiffThreshold, defaults.energyDiffThreshold);
    EXPECT_EQ(oa.nConsecutiveNcgStuck, 0);
}

TEST(ModelCoverageOptimizationAlgorithmTest, ResetClearsOnlyTheStuckFlag)
{
    OptimizationAlgorithm oa;
    oa.usingNCG = true;
    oa.isNCGstuck = true;
    oa.reset_NCG_stuck();
    EXPECT_FALSE(oa.isNCGstuck);
    EXPECT_TRUE(oa.usingNCG);

    oa.usingNCG = false;
    oa.isNCGstuck = true;
    oa.reset_NCG_stuck();
    EXPECT_FALSE(oa.isNCGstuck);
    EXPECT_FALSE(oa.usingNCG);
}

TEST(ModelCoverageOptimizationAlgorithmTest, NcgIsDisabledOnlyAfterMoreThanThresholdConsecutiveStalls)
{
    OptimizationAlgorithm oa;
    oa.usingNCG = true;
    oa.nConsecutiveNcgStuckThreshold = 3;

    testing::internal::CaptureStdout();
    oa.isNCGstuck = true;
    for (int stall = 1; stall <= 3; stall++)
    {
        EXPECT_TRUE(oa.disable_NCG_if_stuck_consecutively()) << "stall " << stall;
        EXPECT_EQ(oa.nConsecutiveNcgStuck, stall);
        EXPECT_TRUE(oa.usingNCG);
    }
    // The fourth consecutive stall exceeds the threshold.
    EXPECT_FALSE(oa.disable_NCG_if_stuck_consecutively());
    EXPECT_EQ(oa.nConsecutiveNcgStuck, 4);
    EXPECT_FALSE(oa.usingNCG);
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_NE(out.find("Disabled NCG"), std::string::npos) << out;
}

TEST(ModelCoverageOptimizationAlgorithmTest, ASuccessfulSearchResetsTheStallCounter)
{
    OptimizationAlgorithm oa;
    oa.usingNCG = true;
    oa.nConsecutiveNcgStuckThreshold = 3;

    oa.isNCGstuck = true;
    oa.disable_NCG_if_stuck_consecutively();
    oa.disable_NCG_if_stuck_consecutively();
    ASSERT_EQ(oa.nConsecutiveNcgStuck, 2);

    oa.isNCGstuck = false;
    EXPECT_TRUE(oa.disable_NCG_if_stuck_consecutively());
    EXPECT_EQ(oa.nConsecutiveNcgStuck, 0);

    // The count starts again: three more stalls are still tolerated.
    oa.isNCGstuck = true;
    for (int stall = 1; stall <= 3; stall++)
    {
        EXPECT_TRUE(oa.disable_NCG_if_stuck_consecutively());
    }
    EXPECT_EQ(oa.nConsecutiveNcgStuck, 3);
}

// ---------------------------------------------------------------------------
// Model: construction, reporting, stopping rule
// ---------------------------------------------------------------------------

TEST(ModelCoverageModelTest, ConstructorCopiesForcesAndConvergenceThresholds)
{
    Param param;
    configure_free_sheet(param);
    param.deltaEnergyConverge = 2.5e-4;
    param.deltaForceScaleConverge = 7.5e-3;
    Mesh mesh(param);
    prepare_sheet(mesh);
    lift_centre_and_rebase(mesh, 0.5);

    param.VERBOSE_MODE = true;
    Record record(2);
    testing::internal::CaptureStdout();
    Model model(mesh, record);
    const std::string out = testing::internal::GetCapturedStdout();

    EXPECT_NE(out.find("[Model::Model()] Model set up with mesh and record"), std::string::npos) << out;
    EXPECT_DOUBLE_EQ(model.oa.energyDiffThreshold, 2.5e-4);
    EXPECT_DOUBLE_EQ(model.oa.forceDiffThreshold, 7.5e-3);
    EXPECT_EQ(model.iteration, 0);
    EXPECT_EQ(model.stepSize, 0.0);
    ASSERT_EQ(model.ncgDirection0.size(), mesh.vertices.size());
    bool anyNonZero = false;
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            const double f = mesh.vertices[i].force.forceTotal.get(axis, 0);
            EXPECT_EQ(model.ncgDirection0[i].forceTotal.get(axis, 0), f);
            anyNonZero = anyNonZero || f != 0.0;
        }
    }
    EXPECT_TRUE(anyNonZero) << "the lifted vertex should leave a force on the sheet";
}

TEST(ModelCoverageModelTest, QuietConstructorPrintsNothing)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    Record record;
    testing::internal::CaptureStdout();
    Model model(mesh, record);
    EXPECT_EQ(testing::internal::GetCapturedStdout(), "");
    EXPECT_TRUE(model.ncgDirection0.empty());
}

TEST(ModelCoverageModelTest, ToStringCurrentStepReportsIterationAndStepSizes)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    Record record;
    Model model(mesh, record);
    model.iteration = 7;
    model.oa.trialStepSize = 0.25;
    model.stepSize = 0.125;
    EXPECT_EQ(model.to_string_current_step(),
              "step: 7, trial StepSize = 0.250000, StepSize = 0.125000");
}

namespace
{
/// A model on an empty mesh whose record holds `n` iterations of the given energy and mean force.
struct StoppingFixture
{
    Param param;
    Mesh mesh;
    Record record;
    Model model;

    explicit StoppingFixture(int maxIterations)
        : param(make_param(maxIterations)), mesh(param), record(maxIterations), model(mesh, record)
    {
    }

    static Param make_param(int maxIterations)
    {
        Param p;
        p.VERBOSE_MODE = false;
        p.maxIterations = maxIterations;
        p.deltaEnergyConverge = 1e-5;
        p.deltaForceScaleConverge = 1e-4;
        return p;
    }

    void fill(int n, double energySlope, double meanForce)
    {
        record.areaTotal.clear();
        record.energyVec.clear();
        record.meanForce.clear();
        for (int i = 0; i < n; i++)
        {
            record.add(100.0, energy_with_total(50.0 + energySlope * i), meanForce);
        }
    }
};
} // namespace

TEST(ModelCoverageStoppingTest, ContinuesDuringTheFirst500IterationsEvenIfConverged)
{
    StoppingFixture f(2000);
    f.fill(600, 0.0, 0.0);
    for (int iteration : {0, 1, 250, 500})
    {
        f.model.iteration = iteration;
        EXPECT_TRUE(f.model.should_continue_optimization()) << "iteration " << iteration;
        EXPECT_FALSE(f.model.oa.isCriteriaSatisfied) << "iteration " << iteration;
    }
}

TEST(ModelCoverageStoppingTest, StopsWhenEnergyPlateausAndForceIsSmall)
{
    StoppingFixture f(2000);
    // Energy drifting by 1e-7 per iteration (below 1e-5) and a mean force below 1e-4.
    f.fill(600, 1e-7, 5e-5);
    f.model.iteration = 501;
    EXPECT_FALSE(f.model.should_continue_optimization());
    EXPECT_TRUE(f.model.oa.isCriteriaSatisfied);

    // A decreasing energy plateaus just the same: the criterion is on |dE|.
    f.fill(600, -1e-7, 5e-5);
    EXPECT_FALSE(f.model.should_continue_optimization());
    EXPECT_TRUE(f.model.oa.isCriteriaSatisfied);
}

TEST(ModelCoverageStoppingTest, ContinuesWhileEnergyStillChangesOrForceIsLarge)
{
    StoppingFixture f(2000);
    f.model.iteration = 501;

    // Average change of 0.01 per iteration over the last 500: not converged.
    f.fill(600, 0.01, 5e-5);
    EXPECT_TRUE(f.model.should_continue_optimization());
    EXPECT_FALSE(f.model.oa.isCriteriaSatisfied);

    // Energy flat, but the mean force is above its threshold.
    f.fill(600, 0.0, 1e-3);
    EXPECT_TRUE(f.model.should_continue_optimization());
    EXPECT_FALSE(f.model.oa.isCriteriaSatisfied);
}

TEST(ModelCoverageStoppingTest, StopsAtTheIterationLimit)
{
    StoppingFixture f(100);
    f.fill(100, 1.0, 10.0); // never converged
    f.model.iteration = 98;
    EXPECT_TRUE(f.model.should_continue_optimization());
    f.model.iteration = 99; // maxIterations - 1
    EXPECT_FALSE(f.model.should_continue_optimization());
    EXPECT_FALSE(f.model.oa.isCriteriaSatisfied);
}

// ---------------------------------------------------------------------------
// Model: trial step size and NCG direction
// ---------------------------------------------------------------------------

namespace
{
void set_all_forces(Mesh &mesh, double x, double y, double z)
{
    for (Vertex &vertex : mesh.vertices)
    {
        set_vector(vertex.force.forceTotal, x, y, z);
    }
}
} // namespace

TEST(ModelCoverageTrialStepTest, RuleOfThumbFromTheLargestForce)
{
    Param param;
    configure_free_sheet(param);
    Mesh mesh(param);
    prepare_sheet(mesh);
    Record record(2);
    Model model(mesh, record);
    const double lFace = mesh.param.lFace;
    model.iteration = 0;

    testing::internal::CaptureStdout();
    // No force at all: no direction to search in.
    set_all_forces(mesh, 0.0, 0.0, 0.0);
    model.determine_trial_step_size();
    EXPECT_EQ(model.oa.trialStepSize, -1.0);

    // Largest force below 5: the step is capped at 0.1 lFace / 5.
    set_all_forces(mesh, 0.0, 0.0, 0.0);
    set_vector(mesh.vertices[3].force.forceTotal, 0.0, 0.0, 2.0);
    model.determine_trial_step_size();
    EXPECT_DOUBLE_EQ(model.oa.trialStepSize, 0.1 * lFace / 5.0);

    // Largest force 10 (a 6-8-0 vector): 0.1 lFace / |F|max.
    set_vector(mesh.vertices[5].force.forceTotal, 6.0, 8.0, 0.0);
    model.determine_trial_step_size();
    EXPECT_DOUBLE_EQ(model.oa.trialStepSize, 0.1 * lFace / 10.0);

    // Between recomputations the trial step grows by 10% per iteration.
    model.iteration = 1;
    model.determine_trial_step_size();
    EXPECT_DOUBLE_EQ(model.oa.trialStepSize, 1.1 * 0.1 * lFace / 10.0);
    model.iteration = 2;
    model.determine_trial_step_size();
    EXPECT_DOUBLE_EQ(model.oa.trialStepSize, 1.1 * 1.1 * 0.1 * lFace / 10.0);

    // ...and is recomputed from the force on every trialIterationInterval-th iteration.
    model.iteration = model.oa.trialIterationInterval;
    model.determine_trial_step_size();
    EXPECT_DOUBLE_EQ(model.oa.trialStepSize, 0.1 * lFace / 10.0);
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_NE(out.find("max_force_scale=10"), std::string::npos) << out;
}

TEST(ModelCoverageNcgDirectionTest, WithoutNcgTheDirectionIsTheCurrentForce)
{
    Param param;
    configure_free_sheet(param);
    Mesh mesh(param);
    prepare_sheet(mesh);
    Record record(2);
    Model model(mesh, record);
    model.oa.usingNCG = false;

    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        set_vector(mesh.vertices[i].forcePrev.forceTotal, 1.0, 0.0, 0.0);
        set_vector(mesh.vertices[i].force.forceTotal, 0.5 * i, -1.0, 2.0);
        set_vector(model.ncgDirection0[i].forceTotal, 9.0, 9.0, 9.0);
    }
    model.update_ncg_direction();
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        EXPECT_EQ(model.ncgDirection0[i].forceTotal.get(0, 0), 0.5 * i);
        EXPECT_EQ(model.ncgDirection0[i].forceTotal.get(1, 0), -1.0);
        EXPECT_EQ(model.ncgDirection0[i].forceTotal.get(2, 0), 2.0);
    }
}

TEST(ModelCoverageNcgDirectionTest, WithNcgTheOldDirectionIsAddedWithTheFletcherReevesFactor)
{
    Param param;
    configure_free_sheet(param);
    Mesh mesh(param);
    prepare_sheet(mesh);
    Record record(2);
    Model model(mesh, record);
    model.oa.usingNCG = true;

    // |F|^2 / |Fprev|^2 = (n * 4) / (n * 1) = 4, exact in floating point.
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        set_vector(mesh.vertices[i].forcePrev.forceTotal, 1.0, 0.0, 0.0);
        set_vector(mesh.vertices[i].force.forceTotal, 0.0, 2.0, 0.0);
        set_vector(model.ncgDirection0[i].forceTotal, 0.0, 0.0, 1.0);
    }
    model.update_ncg_direction();
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        EXPECT_DOUBLE_EQ(model.ncgDirection0[i].forceTotal.get(0, 0), 0.0);
        EXPECT_DOUBLE_EQ(model.ncgDirection0[i].forceTotal.get(1, 0), 2.0);
        EXPECT_DOUBLE_EQ(model.ncgDirection0[i].forceTotal.get(2, 0), 4.0);
    }
}

TEST(ModelCoverageNcgDirectionTest, ResetZeroesTheDirectionAndTheNextUpdateIsTheForce)
{
    Param param;
    configure_free_sheet(param);
    Mesh mesh(param);
    prepare_sheet(mesh);
    lift_centre_and_rebase(mesh, 0.5);
    Record record(2);
    Model model(mesh, record);
    model.oa.usingNCG = true;

    model.reset_ncg_direction();
    for (const Force &direction : model.ncgDirection0)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            EXPECT_EQ(direction.forceTotal.get(axis, 0), 0.0);
        }
    }

    // From a zero direction even the NCG update reduces to the current force.
    model.update_ncg_direction();
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            EXPECT_DOUBLE_EQ(model.ncgDirection0[i].forceTotal.get(axis, 0),
                             mesh.vertices[i].force.forceTotal.get(axis, 0));
        }
    }
}

// ---------------------------------------------------------------------------
// Model: line search
// ---------------------------------------------------------------------------

TEST(ModelCoverageLineSearchTest, NcgBranchAcceptsAStepThatSatisfiesSufficientDecrease)
{
    Param param;
    configure_free_sheet(param);
    Mesh mesh(param);
    prepare_sheet(mesh);
    lift_centre_and_rebase(mesh, 1.0);
    const double energyBefore = mesh.param.energy.energyTotal;

    Record record(2);
    record.add(mesh.param.area, mesh.param.energy, mesh.calculate_mean_force());
    Model model(mesh, record);
    model.oa.usingNCG = true;
    // The curvature test compares the previous force with itself (see the
    // report), so it only passes with c2 >= 1; that leaves the Armijo test.
    model.oa.c2 = 1.0;
    testing::internal::CaptureStdout();
    model.determine_trial_step_size();
    const double trial = model.oa.trialStepSize;
    ASSERT_GT(trial, 0.0);
    const double slope = -sum_force_dot_direction(mesh, model.ncgDirection0);
    ASSERT_LT(slope, 0.0) << "the force direction must be a descent direction";

    const double step = model.linear_search_for_stepsize_to_minimize_energy();
    testing::internal::GetCapturedStdout();

    EXPECT_GT(step, 0.0);
    EXPECT_EQ(step, model.stepSize);
    EXPECT_FALSE(model.oa.isNCGstuck) << "the NCG branch should have accepted a step";
    // The step is the trial step scaled down by 0.8 a whole number of times (at least once).
    const double k = std::log(step / trial) / std::log(0.8);
    EXPECT_GE(k, 1.0 - 1e-9);
    EXPECT_NEAR(k, std::round(k), 1e-9);
    // The accepted energy is lower than the start and satisfies Armijo.
    const double energyAfter = mesh.param.energy.energyTotal;
    EXPECT_LT(energyAfter, energyBefore);
    EXPECT_LE(energyAfter, energyBefore + model.oa.c1 * step * slope + 1e-12 * std::abs(energyBefore));
}

TEST(ModelCoverageLineSearchTest, SimpleSearchAlongTheForceLowersTheEnergy)
{
    Param param;
    configure_free_sheet(param);
    Mesh mesh(param);
    prepare_sheet(mesh);
    lift_centre_and_rebase(mesh, 1.0);
    const double energyBefore = mesh.param.energy.energyTotal;

    Record record(2);
    Model model(mesh, record);
    model.oa.usingNCG = false;
    testing::internal::CaptureStdout();
    model.determine_trial_step_size();
    const double step = model.linear_search_for_stepsize_to_minimize_energy();
    testing::internal::GetCapturedStdout();

    ASSERT_GT(step, 0.0);
    EXPECT_LT(mesh.param.energy.energyTotal, energyBefore);

    // The coordinates the search left are the start plus step * direction.
    const int centre = central_mobile_vertex(mesh);
    const Vertex &vertex = mesh.vertices[centre];
    EXPECT_NEAR(vertex.coord.get(2, 0),
                vertex.coordPrev.get(2, 0) + step * model.ncgDirection0[centre].forceTotal.get(2, 0),
                1e-12);
}

TEST(ModelCoverageLineSearchTest, SimpleSearchUphillGivesUpWithMinusOne)
{
    Param param;
    configure_free_sheet(param);
    Mesh mesh(param);
    prepare_sheet(mesh);
    lift_centre_and_rebase(mesh, 1.0);
    const double energyBefore = mesh.param.energy.energyTotal;

    Record record(2);
    Model model(mesh, record);
    model.oa.usingNCG = false;
    model.oa.stepThreshold = 1e-3; // keep the number of trial steps small
    testing::internal::CaptureStdout();
    model.determine_trial_step_size();
    // Search against the force: every step raises the energy.
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        model.ncgDirection0[i].forceTotal = -1.0 * mesh.vertices[i].forcePrev.forceTotal;
    }
    const double step = model.linear_search_for_stepsize_to_minimize_energy();
    const std::string out = testing::internal::GetCapturedStdout();

    EXPECT_EQ(step, -1.0);
    EXPECT_NE(out.find("cannot find an efficient small stepSize"), std::string::npos) << out;
    EXPECT_GE(mesh.param.energy.energyTotal, energyBefore);
}

// ---------------------------------------------------------------------------
// Model: boundary enforcement after a coordinate update (global Fixed mode)
// ---------------------------------------------------------------------------

TEST(ModelCoverageBoundaryTest, FixedModeRestoresTheVerticesOfBoundaryFaces)
{
    Param param;
    configure_free_sheet(param);
    param.boundaryCondition = BoundaryType::Fixed;
    Mesh mesh(param);
    mesh.setup_flat();
    mesh.update_previous_coord_for_vertex();
    ASSERT_GE(mesh.faces.size(), 4u);

    // Nothing in setup marks a face as boundary; mark two.
    for (Face &face : mesh.faces)
    {
        face.isBoundary = false;
    }
    mesh.faces[0].isBoundary = true;
    mesh.faces[mesh.faces.size() - 1].isBoundary = true;
    std::vector<char> clamped(mesh.vertices.size(), 0);
    for (int f : {0, static_cast<int>(mesh.faces.size()) - 1})
    {
        for (int corner : mesh.faces[f].adjacentVertices)
        {
            clamped[corner] = 1;
        }
    }

    Record record(1);
    Model model(mesh, record);
    model.stepSize = 0.5;
    for (Force &direction : model.ncgDirection0)
    {
        set_vector(direction.forceTotal, 1.0, -2.0, 4.0);
    }
    model.update_vertex_using_NCG();

    int nClamped = 0;
    for (std::size_t v = 0; v < mesh.vertices.size(); v++)
    {
        const Vertex &vertex = mesh.vertices[v];
        if (clamped[v])
        {
            nClamped++;
            for (int axis = 0; axis < 3; axis++)
            {
                EXPECT_EQ(vertex.coord.get(axis, 0), vertex.coordPrev.get(axis, 0)) << "vertex " << v;
            }
        }
        else
        {
            EXPECT_DOUBLE_EQ(vertex.coord.get(0, 0), vertex.coordPrev.get(0, 0) + 0.5);
            EXPECT_DOUBLE_EQ(vertex.coord.get(1, 0), vertex.coordPrev.get(1, 0) - 1.0);
            EXPECT_DOUBLE_EQ(vertex.coord.get(2, 0), vertex.coordPrev.get(2, 0) + 2.0);
        }
    }
    EXPECT_EQ(nClamped, 6);
}

// ---------------------------------------------------------------------------
// Model: thermal (Metropolis) trial move
// ---------------------------------------------------------------------------

namespace
{
void configure_thermal(Param &param)
{
    configure_free_sheet(param);
    param.thermalFluctuationEnabled = true;
    param.thermalFluctuationInterval = 1;
    param.thermalFluctuationTemperatureKelvin = 298.0;
    param.thermalFluctuationMinTemperatureKelvin = 0.0;
    param.thermalFluctuationCoolingRate = 1.0;
    param.thermalFluctuationStepScale = 0.02;
}

std::vector<Matrix> snapshot_coords(const Mesh &mesh)
{
    std::vector<Matrix> coords;
    coords.reserve(mesh.vertices.size());
    for (const Vertex &vertex : mesh.vertices)
    {
        coords.push_back(vertex.coord);
    }
    return coords;
}

bool coords_equal(const Mesh &mesh, const std::vector<Matrix> &coords)
{
    for (std::size_t v = 0; v < mesh.vertices.size(); v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            if (mesh.vertices[v].coord.get(axis, 0) != coords[v].get(axis, 0))
            {
                return false;
            }
        }
    }
    return true;
}
} // namespace

TEST(ModelCoverageThermalTest, OffIntervalIterationsAndNonPositiveIntervalsSkipTheTrial)
{
    Param param;
    configure_thermal(param);
    param.thermalFluctuationInterval = 3;
    Mesh mesh(param);
    Record record(1);
    Model model(mesh, record);

    model.iteration = 4;
    EXPECT_FALSE(model.simulated_annealing_next_step());
    param.thermalFluctuationInterval = 0;
    model.iteration = 3;
    EXPECT_FALSE(model.simulated_annealing_next_step());
    EXPECT_TRUE(model.thermalFluctuationRecords.empty());
    EXPECT_EQ(model.thermalFluctuationAttemptCount, 0);
}

TEST(ModelCoverageThermalTest, ZeroStepOrZeroTemperatureMakesNoAttempt)
{
    Param param;
    configure_thermal(param);
    Mesh mesh(param);
    prepare_sheet(mesh);
    const std::vector<Matrix> before = snapshot_coords(mesh);
    Record record(1);
    Model model(mesh, record);

    param.thermalFluctuationStepScale = 0.0;
    EXPECT_FALSE(model.simulated_annealing_next_step(true));

    param.thermalFluctuationStepScale = 0.02;
    param.thermalFluctuationTemperatureKelvin = 0.0;
    EXPECT_FALSE(model.simulated_annealing_next_step(true));

    EXPECT_TRUE(model.thermalFluctuationRecords.empty());
    EXPECT_EQ(model.thermalFluctuationAttemptCount, 0);
    EXPECT_TRUE(coords_equal(mesh, before));
}

TEST(ModelCoverageThermalTest, NoIndependentVertexMeansNoAttempt)
{
    Param param;
    configure_thermal(param);
    Mesh mesh(param);
    prepare_sheet(mesh);
    for (Vertex &vertex : mesh.vertices)
    {
        vertex.isGhost = true;
    }
    const std::vector<Matrix> before = snapshot_coords(mesh);
    Record record(1);
    Model model(mesh, record);

    EXPECT_FALSE(model.simulated_annealing_next_step(true));
    EXPECT_TRUE(model.thermalFluctuationRecords.empty());
    EXPECT_TRUE(coords_equal(mesh, before));
}

TEST(ModelCoverageThermalTest, CoolingScheduleIsFlooredAtTheMinimumTemperature)
{
    Param param;
    configure_thermal(param);
    param.thermalFluctuationTemperatureKelvin = 400.0;
    param.thermalFluctuationCoolingRate = 0.5;
    param.thermalFluctuationMinTemperatureKelvin = 150.0;
    Mesh mesh(param);
    prepare_sheet(mesh);
    Record record(1);
    Model model(mesh, record);

    testing::internal::CaptureStdout();
    for (int attempt = 0; attempt < 3; attempt++)
    {
        model.iteration = attempt;
        ASSERT_TRUE(model.simulated_annealing_next_step());
    }
    testing::internal::GetCapturedStdout();

    ASSERT_EQ(model.thermalFluctuationRecords.size(), 3u);
    EXPECT_EQ(model.thermalFluctuationAttemptCount, 3);
    EXPECT_DOUBLE_EQ(model.thermalFluctuationRecords[0].temperatureKelvin, 400.0);
    EXPECT_DOUBLE_EQ(model.thermalFluctuationRecords[1].temperatureKelvin, 200.0);
    EXPECT_DOUBLE_EQ(model.thermalFluctuationRecords[2].temperatureKelvin, 150.0); // 100 floored
    EXPECT_DOUBLE_EQ(model.thermalFluctuationRecords[2].kBT, 0.01380649 * 150.0);
    for (int attempt = 0; attempt < 3; attempt++)
    {
        EXPECT_EQ(model.thermalFluctuationRecords[attempt].iteration, attempt);
    }
}

TEST(ModelCoverageThermalTest, AnUphillMoveAtVanishingTemperatureIsRejectedAndUndone)
{
    Param param;
    configure_thermal(param);
    // The flat sheet is the minimum: any move costs energy, and kBT ~ 1e-32 pN nm
    // makes exp(-dE / kBT) exactly 0.
    param.thermalFluctuationTemperatureKelvin = 1e-30;
    Mesh mesh(param);
    prepare_sheet(mesh);
    const std::vector<Matrix> before = snapshot_coords(mesh);
    const double energyBefore = mesh.param.energy.energyTotal;
    Record record(1);
    Model model(mesh, record);

    testing::internal::CaptureStdout();
    ASSERT_TRUE(model.simulated_annealing_next_step(true));
    testing::internal::GetCapturedStdout();

    ASSERT_EQ(model.thermalFluctuationRecords.size(), 1u);
    const Model::ThermalFluctuationRecord &trial = model.thermalFluctuationRecords[0];
    EXPECT_GT(trial.deltaEnergy, 0.0);
    EXPECT_EQ(trial.acceptanceProbability, 0.0);
    EXPECT_FALSE(trial.accepted);
    EXPECT_TRUE(coords_equal(mesh, before));
    EXPECT_NEAR(mesh.param.energy.energyTotal, energyBefore, 1e-12 * std::max(1.0, std::abs(energyBefore)));
}

TEST(ModelCoverageThermalTest, ANonFiniteEnergyChangeIsRejected)
{
    Param param;
    configure_thermal(param);
    Mesh mesh(param);
    prepare_sheet(mesh);
    const std::vector<Matrix> before = snapshot_coords(mesh);
    const double energyBefore = mesh.param.energy.energyTotal;
    Record record(1);
    Model model(mesh, record);

    // A poisoned starting energy makes deltaEnergy NaN whatever the move does.
    mesh.param.energy.energyTotal = std::nan("");
    testing::internal::CaptureStdout();
    ASSERT_TRUE(model.simulated_annealing_next_step(true));
    testing::internal::GetCapturedStdout();

    ASSERT_EQ(model.thermalFluctuationRecords.size(), 1u);
    const Model::ThermalFluctuationRecord &trial = model.thermalFluctuationRecords[0];
    EXPECT_FALSE(std::isfinite(trial.deltaEnergy));
    EXPECT_EQ(trial.acceptanceProbability, 0.0);
    EXPECT_FALSE(trial.accepted);
    EXPECT_TRUE(coords_equal(mesh, before));
    // The energy is recomputed at the restored coordinates.
    EXPECT_NEAR(mesh.param.energy.energyTotal, energyBefore, 1e-12 * std::max(1.0, std::abs(energyBefore)));
}

/**
 * The Gaussian displacement is clamped to +-3 sigma per component. Scan seeds
 * by replaying the exact draws simulated_annealing_next_step() makes (same
 * engine, same distributions, same order), so the test finds a clamping seed
 * on whichever standard library it is built with.
 */
TEST(ModelCoverageThermalTest, DisplacementsAreClampedToThreeSigma)
{
    Param param;
    configure_thermal(param);
    param.thermalFluctuationTemperatureKelvin = 1e30; // accept every move
    Mesh mesh(param);
    prepare_sheet(mesh);

    std::vector<int> mobile;
    for (int i = 0; i < static_cast<int>(mesh.vertices.size()); i++)
    {
        if (mesh.is_independent_vertex(i))
        {
            mobile.push_back(i);
        }
    }
    ASSERT_FALSE(mobile.empty());
    const double sigma = param.thermalFluctuationStepScale * param.lFace;

    struct Draw
    {
        unsigned seed = 0;
        int vertex = -1;
        double d[3] = {0.0, 0.0, 0.0};
    };
    auto replay = [&](unsigned seed) {
        std::mt19937 rng(seed);
        std::uniform_int_distribution<int> vertexDistribution(0, mobile.size() - 1);
        Draw draw;
        draw.seed = seed;
        draw.vertex = mobile[vertexDistribution(rng)];
        std::normal_distribution<double> normal(0.0, sigma);
        for (int j = 0; j < 3; j++)
        {
            draw.d[j] = normal(rng);
        }
        return draw;
    };
    Draw high, low;
    for (unsigned seed = 1; seed < 200000 && (high.vertex < 0 || low.vertex < 0); seed++)
    {
        const Draw draw = replay(seed);
        bool anyHigh = false, anyLow = false;
        for (double d : draw.d)
        {
            anyHigh = anyHigh || d > 3.0 * sigma;
            anyLow = anyLow || d < -3.0 * sigma;
        }
        if (anyHigh && high.vertex < 0)
            high = draw;
        if (anyLow && low.vertex < 0)
            low = draw;
    }
    ASSERT_GE(high.vertex, 0);
    ASSERT_GE(low.vertex, 0);

    Record record(1);
    Model model(mesh, record);
    testing::internal::CaptureStdout();
    for (const Draw &draw : {high, low})
    {
        const std::vector<Matrix> before = snapshot_coords(mesh);
        model.thermalRng.seed(draw.seed);
        ASSERT_TRUE(model.simulated_annealing_next_step(true));
        ASSERT_TRUE(model.thermalFluctuationRecords.back().accepted);
        for (int j = 0; j < 3; j++)
        {
            const double expected = std::max(-3.0 * sigma, std::min(3.0 * sigma, draw.d[j]));
            const double moved = mesh.vertices[draw.vertex].coord.get(j, 0) - before[draw.vertex].get(j, 0);
            EXPECT_NEAR(moved, expected, 1e-12) << "seed " << draw.seed << " axis " << j;
            EXPECT_LE(std::abs(moved), 3.0 * sigma + 1e-12);
        }
    }
    testing::internal::GetCapturedStdout();
}
