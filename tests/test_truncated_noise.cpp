#include "test_truncated_noise.hpp"

/**
 * The truncated noise is a rejection sampler on the counter-based generator,
 * switched on by one key. What can go wrong: the sampler might leak past its
 * cutoff, change draws it should have kept, or get the variance it divides
 * out wrong; the keys might accept a cutoff the sampler cannot serve; and the
 * step might not use the sampler the way the parameters say, or might change
 * a run that leaves the key at 0.
 *
 * @see docs/resolution_limits.md section 2.5
 */

namespace
{

/// A clamped 60 nm sheet under the per-vertex boundary, as in the external
/// force tests: no images, so every vertex is either integrated or held.
void configure_clamped_sheet(Param &param)
{
    param.VERBOSE_MODE = false;
    param.sideX = 60.0;
    param.sideY = 60.0;
    param.lFace = 5.0;
    param.kCurv = 83.4;
    param.KBT = 4.17;
    param.uSurf = 250.0;
    param.uVol = 0.0;
    param.isGlobalConstraint = true;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
    param.randomSeed = 20261021u;
    param.timeStep = 1.0e-3;
    param.diffConst = 1.0;
    param.surfaceSolver = "iterative";
    param.boundaryCondition = BoundaryType::Mixed;
    param.boundaryConditionX = BoundaryType::Fixed;
    param.boundaryConditionY = BoundaryType::Fixed;
    param.fixedBoundaryRings = 2;
    param.inPlaneDynamicsEnabled = true;
}

/// A mesh, its record and its model, built the way run_dynamics_flat() builds them.
struct StepRun
{
    Param param;
    std::unique_ptr<DynamicMesh> mesh;
    Record record{16};
    std::unique_ptr<DynamicModel> model;

    explicit StepRun(const Param &source) : param(source)
    {
        mesh.reset(new DynamicMesh(param));
        mesh->setup_flat();
        for (Vertex &vertex : mesh->vertices)
        {
            vertex.coord.set(2, 0, 10.0);
        }
        mesh->calculate_element_area_volume();
        mesh->sum_membrane_area_and_volume(mesh->param.area0, mesh->param.vol0);
        mesh->param.vol0 = 0.0;
        mesh->update_previous_coord_for_vertex();
        mesh->update_reference_coord_from_previous_coord();
        mesh->Compute_Energy_And_Force();
        mesh->update_previous_force_for_vertex();
        mesh->update_vertices_mat_with_vector();
        mesh->apply_mesh_to_surface();
        record.add(mesh->param.area, mesh->param.energy, mesh->calculate_mean_force());
        model.reset(new DynamicModel(*mesh, record));
    }

    /// One Brownian step from the setup configuration, which the control net
    /// never leaves, so two runs with the same seed see the same forces.
    void step(int iteration)
    {
        mesh->apply_mesh_to_surface();
        model->iteration = iteration;
        model->next_step();
    }
};

/// The key DynamicModel::next_step() draws from at this iteration.
std::uint64_t step_key(unsigned int seed, int iteration)
{
    return slimed::splitmix64(static_cast<std::uint64_t>(seed)) + static_cast<std::uint64_t>(iteration);
}

} // namespace

TEST(TruncatedNoiseTest, ZeroCutoffIsThePlainGaussian)
{
    for (std::uint64_t k = 0; k < 2000; k++)
    {
        const double plain = slimed::standard_normal(k * 977ULL, k % 300, k % 3);
        EXPECT_EQ(slimed::truncated_standard_normal(k * 977ULL, k % 300, k % 3, 0.0), plain);
        EXPECT_EQ(slimed::truncated_standard_normal(k * 977ULL, k % 300, k % 3, -1.0), plain);
    }
    EXPECT_EQ(slimed::truncated_normal_variance(0.0), 1.0);
}

TEST(TruncatedNoiseTest, DrawsStayInsideTheCutoffAndKeepTheOnesAlreadyThere)
{
    for (const double cutoff : {0.1, 1.0, 2.0, 3.0})
    {
        int kept = 0;
        int redrawn = 0;
        for (std::uint64_t k = 0; k < 20000; k++)
        {
            const double plain = slimed::standard_normal(k, k % 500, k % 3);
            const double cut = slimed::truncated_standard_normal(k, k % 500, k % 3, cutoff);
            ASSERT_LT(std::fabs(cut), cutoff);
            if (std::fabs(plain) < cutoff)
            {
                EXPECT_EQ(cut, plain);
                kept++;
            }
            else
            {
                redrawn++;
            }
        }
        // Both branches were exercised at every cutoff.
        EXPECT_GT(kept, 0);
        EXPECT_GT(redrawn, 0);
    }
}

TEST(TruncatedNoiseTest, VarianceMatchesTheClosedForm)
{
    // Known values of 1 - 2 c phi(c) / erf(c / sqrt 2).
    EXPECT_NEAR(slimed::truncated_normal_variance(1.0), 0.29113, 1e-5);
    EXPECT_NEAR(slimed::truncated_normal_variance(2.0), 0.77374, 1e-5);
    EXPECT_NEAR(slimed::truncated_normal_variance(3.0), 0.97334, 1e-5);
    EXPECT_NEAR(slimed::truncated_normal_variance(8.0), 1.0, 1e-12);

    // And the sampler has it: the second moment over 400000 draws, whose
    // standard error is below 0.3 % at every cutoff here.
    for (const double cutoff : {1.0, 2.0, 3.0})
    {
        const int n = 400000;
        double sum = 0.0;
        double sumSquares = 0.0;
        for (int k = 0; k < n; k++)
        {
            const double x = slimed::truncated_standard_normal(static_cast<std::uint64_t>(k / 1000) * 131ULL,
                                                               static_cast<std::uint64_t>(k % 1000), 2, cutoff);
            sum += x;
            sumSquares += x * x;
        }
        const double variance = slimed::truncated_normal_variance(cutoff);
        EXPECT_NEAR(sum / n, 0.0, 0.01) << "cutoff " << cutoff;
        EXPECT_NEAR(sumSquares / n / variance, 1.0, 0.01) << "cutoff " << cutoff;
    }
}

TEST(TruncatedNoiseTest, TheParameterKeysParseAndRejectWhatTheSamplerCannotServe)
{
    Param param;
    EXPECT_EQ(param.noiseTruncation, 0.0);
    EXPECT_TRUE(param.noiseTruncationRenormalize);

    EXPECT_TRUE(import_kv_string("noiseTruncation", "2.5", param));
    EXPECT_DOUBLE_EQ(param.noiseTruncation, 2.5);
    EXPECT_TRUE(import_kv_string("noiseTruncation", "0", param));
    EXPECT_EQ(param.noiseTruncation, 0.0);
    EXPECT_TRUE(import_kv_string("noiseTruncation", "0.1", param));
    EXPECT_DOUBLE_EQ(param.noiseTruncation, 0.1);
    EXPECT_THROW(import_kv_string("noiseTruncation", "-1", param), std::runtime_error);
    EXPECT_THROW(import_kv_string("noiseTruncation", "0.05", param), std::runtime_error);
    EXPECT_THROW(import_kv_string("noiseTruncation", "inf", param), std::runtime_error);

    EXPECT_TRUE(import_kv_string("noiseTruncationRenormalize", "false", param));
    EXPECT_FALSE(param.noiseTruncationRenormalize);
    EXPECT_TRUE(import_kv_string("noiseTruncationRenormalize", "true", param));
    EXPECT_TRUE(param.noiseTruncationRenormalize);
}

/**
 * @brief The step moves each coordinate by exactly sqrt(2 D dt) times the
 * truncated draw more than the Gaussian step does, rescaled or not.
 *
 * Two runs from the same configuration and seed share the force term, so the
 * difference of their surfaces is the noise term alone, and it can be checked
 * draw by draw against the sampler on the key the step uses.
 */
TEST(TruncatedNoiseTest, TheStepUsesTheTruncatedDrawItIsAskedFor)
{
    Param base;
    configure_clamped_sheet(base);
    StepRun gaussian(base);

    for (const bool renormalize : {true, false})
    {
        Param cut = base;
        cut.noiseTruncation = 1.5;
        cut.noiseTruncationRenormalize = renormalize;
        StepRun truncated(cut);
        const double rescale = renormalize ? 1.0 / std::sqrt(slimed::truncated_normal_variance(1.5)) : 1.0;
        EXPECT_DOUBLE_EQ(truncated.model->noiseRescale, rescale);

        const int iteration = 7;
        gaussian.step(iteration);
        truncated.step(iteration);
        const std::uint64_t key = step_key(base.randomSeed, iteration);
        const double scale = gaussian.model->randScaleConst;

        int checked = 0;
        int changed = 0;
        for (int v = 0; v < static_cast<int>(gaussian.mesh->vertices.size()); v++)
        {
            for (int j = 0; j < 3; j++)
            {
                const double difference = truncated.mesh->matSurface(v, j) - gaussian.mesh->matSurface(v, j);
                if (!gaussian.mesh->is_independent_vertex(v))
                {
                    EXPECT_EQ(difference, 0.0);
                    continue;
                }
                const double plain = slimed::standard_normal(key, static_cast<std::uint64_t>(v),
                                                             static_cast<std::uint64_t>(j));
                const double drawn = rescale * slimed::truncated_standard_normal(
                                                   key, static_cast<std::uint64_t>(v),
                                                   static_cast<std::uint64_t>(j), 1.5);
                EXPECT_NEAR(difference, scale * (drawn - plain), 1e-12) << "vertex " << v << " axis " << j;
                checked++;
                changed += (drawn != plain) ? 1 : 0;
            }
        }
        EXPECT_GT(checked, 100);
        EXPECT_GT(changed, 0);
    }
}

/**
 * @brief An explicit cutoff of 0 is the step without the key.
 *
 * Compared to 1e-12 nm rather than bit for bit: under OpenMP the surface
 * solve's reductions differ by a few ulp from call to call even for one
 * configuration. The bit-for-bit check is the byte-identity gate on the
 * shipped workload.
 */
TEST(TruncatedNoiseTest, ACutoffOfZeroIsTheGaussianStep)
{
    Param base;
    configure_clamped_sheet(base);
    Param explicitZero = base;
    explicitZero.noiseTruncation = 0.0;
    explicitZero.noiseTruncationRenormalize = false;
    StepRun a(base);
    StepRun b(explicitZero);
    for (int iteration = 0; iteration < 3; iteration++)
    {
        a.step(iteration);
        b.step(iteration);
        for (int v = 0; v < static_cast<int>(a.mesh->vertices.size()); v++)
        {
            for (int j = 0; j < 3; j++)
            {
                ASSERT_NEAR(a.mesh->matSurface(v, j), b.mesh->matSurface(v, j), 1e-12);
            }
        }
    }
}
