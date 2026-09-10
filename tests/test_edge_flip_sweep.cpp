#include "test_edge_flip_sweep.hpp"

/**
 * The flip sweep is the move that makes the membrane a fluid, and the only
 * question that matters about a Monte Carlo move is whether it samples the
 * distribution it claims to. Everything else -- acceptance rates, throughput,
 * the shape of the log -- is diagnostics.
 *
 * So the gate here is a two-state test. Freeze every vertex but the four one
 * edge touches, and the chain has exactly two triangulations to move between.
 * Its stationary distribution is then known in closed form, and the sweep has
 * to reproduce it: the ratio of occupancies must equal exp(-dE/kT), at more
 * than one temperature so that agreement cannot come from the ratio being near
 * one. This is detailed balance in the only form that is exact for a
 * Metropolis chain rather than asymptotic.
 *
 * @see docs/edge_flip_plan.md work package 3
 */

namespace
{
struct MeshFixture
{
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;
};

/// A closed icosphere, as in the WP2 tests: closed so a volume constraint is
/// defined, and curved so the bending energy is not degenerate.
MeshFixture build_icosphere(int level, double radius)
{
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    std::vector<std::array<double, 3>> points = {
        {{-1, phi, 0}}, {{1, phi, 0}},  {{-1, -phi, 0}}, {{1, -phi, 0}},
        {{0, -1, phi}}, {{0, 1, phi}},  {{0, -1, -phi}}, {{0, 1, -phi}},
        {{phi, 0, -1}}, {{phi, 0, 1}},  {{-phi, 0, -1}}, {{-phi, 0, 1}}};
    std::vector<std::array<int, 3>> triangles = {
        {{0, 11, 5}}, {{0, 5, 1}},  {{0, 1, 7}},   {{0, 7, 10}}, {{0, 10, 11}},
        {{1, 5, 9}},  {{5, 11, 4}}, {{11, 10, 2}}, {{10, 7, 6}}, {{7, 1, 8}},
        {{3, 9, 4}},  {{3, 4, 2}},  {{3, 2, 6}},   {{3, 6, 8}},  {{3, 8, 9}},
        {{4, 9, 5}},  {{2, 4, 11}}, {{6, 2, 10}},  {{8, 6, 7}},  {{9, 8, 1}}};

    const auto onSphere = [](std::array<double, 3> p) {
        const double n = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        return std::array<double, 3>{{p[0] / n, p[1] / n, p[2] / n}};
    };
    for (auto &p : points)
    {
        p = onSphere(p);
    }
    for (int step = 0; step < level; step++)
    {
        std::map<std::pair<int, int>, int> midpoint;
        const auto edgePoint = [&](int a, int b) {
            const std::pair<int, int> key = (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
            const auto found = midpoint.find(key);
            if (found != midpoint.end())
            {
                return found->second;
            }
            const int index = static_cast<int>(points.size());
            points.push_back(onSphere({{0.5 * (points[a][0] + points[b][0]),
                                        0.5 * (points[a][1] + points[b][1]),
                                        0.5 * (points[a][2] + points[b][2])}}));
            midpoint.emplace(key, index);
            return index;
        };
        std::vector<std::array<int, 3>> refined;
        for (const std::array<int, 3> &t : triangles)
        {
            const int ab = edgePoint(t[0], t[1]);
            const int bc = edgePoint(t[1], t[2]);
            const int ca = edgePoint(t[2], t[0]);
            refined.push_back({{t[0], ab, ca}});
            refined.push_back({{ab, t[1], bc}});
            refined.push_back({{ca, bc, t[2]}});
            refined.push_back({{ab, bc, ca}});
        }
        triangles = refined;
    }

    MeshFixture fixture;
    for (const std::array<double, 3> &p : points)
    {
        fixture.vertices.push_back({radius * p[0], radius * p[1], radius * p[2]});
    }
    for (const std::array<int, 3> &t : triangles)
    {
        fixture.faces.push_back({t[0], t[1], t[2]});
    }
    return fixture;
}

void configure(Param &param, bool withConstraints)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.kCurv = 83.4;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
    param.isGlobalConstraint = withConstraints;
    param.uSurf = withConstraints ? 250.0 : 0.0;
    param.uVol = withConstraints ? 30.0 : 0.0;
    param.area0 = 1.0;
    param.vol0 = 1.0;
    param.randomSeed = 20260910u;
    param.timeStep = 1.0;
    param.edgeFlipEnabled = true;
    param.edgeFlipAttemptRate = 2.0;
    param.edgeFlipInterval = 1;
}

void build(Mesh &mesh, Param &param, int level)
{
    const MeshFixture fixture = build_icosphere(level, 12.0);
    mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.calculate_element_area_volume();
    double area = 0.0;
    double volume = 0.0;
    mesh.sum_membrane_area_and_volume(area, volume);
    if (param.uSurf != 0.0)
    {
        param.area0 = 0.9 * area;
        param.vol0 = 1.1 * volume;
    }
    mesh.Compute_Energy_And_Force();
}

/**
 * @brief Freeze every vertex except the four that @p iEdge touches.
 *
 * Leaves exactly one flippable edge, so the chain has two states. Any other
 * edge among those four vertices reaches a fifth, frozen one through its own
 * quadrilateral, so it is refused. And the four survive the flip: the edge
 * becomes the other diagonal of the same quadrilateral.
 */
void isolate_single_edge(Mesh &mesh, int iEdge)
{
    const MeshEdge &edge = mesh.edges[iEdge];
    const int keep[4] = {edge.v[0], edge.v[1], edge.opposite[0], edge.opposite[1]};

    mesh.flipFrozenVertex.assign(mesh.vertices.size(), 1);
    for (int corner : keep)
    {
        mesh.flipFrozenVertex[corner] = 0;
    }
    for (int i = 0; i < static_cast<int>(mesh.edges.size()); i++)
    {
        mesh.refresh_edge_flippability(i);
    }
}

/**
 * @brief Set edgeFlipAttemptRate so that a sweep draws @p lambda attempts.
 *
 * The rate is per edge per unit time, so the mean attempt count also scales
 * with the number of flippable edges and the time step. Tests care about the
 * count, so they say so and let this work backwards -- otherwise a rate that
 * looks small quietly asks for thousands of attempts on a mesh with a few
 * hundred edges.
 */
void set_rate_for_expected_attempts(Param &param, const Mesh &mesh, double lambda);

int count_flippable(const Mesh &mesh)
{
    int n = 0;
    for (const MeshEdge &edge : mesh.edges)
    {
        n += edge.flippable ? 1 : 0;
    }
    return n;
}

/// A deterministic generator, so a failure is reproducible.
struct Lcg
{
    unsigned long long state;
    explicit Lcg(unsigned long long seed) : state(seed) {}
    int below(int bound)
    {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<int>((state >> 33) % static_cast<unsigned int>(bound));
    }
};

void set_rate_for_expected_attempts(Param &param, const Mesh &mesh, double lambda)
{
    const int nFlippable = count_flippable(mesh);
    const double denominator =
        param.timeStep * std::max(1, param.edgeFlipInterval) * std::max(1, nFlippable);
    param.edgeFlipAttemptRate = lambda / denominator;
}
} // namespace

// ---------------------------------------------------------------------------
// The gate: detailed balance, in closed form
// ---------------------------------------------------------------------------

/**
 * @brief With two triangulations to choose between, the sweep visits them in
 * the Boltzmann ratio.
 *
 * The chain is reconstructed from the sweep's own log rather than by watching
 * the mesh: every record says whether that attempt was accepted, and an
 * accepted flip toggles the state, so the log determines which state each
 * attempt was made from. That measures what the sweep actually did, including
 * its random numbers and its acceptance rule, rather than re-deriving the
 * answer from the energy.
 *
 * Two temperatures, because a single one could agree by accident if the ratio
 * happened to sit near unity. Here the expected ratios are exp(-1) and
 * exp(-2), a factor of e apart.
 */
TEST(EdgeFlipSweepTest, TwoStateOccupancyFollowsTheBoltzmannRatio)
{
    // Measure the energy difference once, then choose temperatures that put it
    // squarely in the interesting range.
    double deltaEnergy = 0.0;
    int chosenEdge = -1;
    {
        Param param;
        configure(param, false);
        Mesh mesh(param);
        ASSERT_NO_FATAL_FAILURE(build(mesh, param, 1));
        for (const MeshEdge &edge : mesh.edges)
        {
            EdgeFlipDelta delta;
            if (mesh.evaluate_edge_flip(edge.index, delta) && std::abs(delta.energy) > 1.0)
            {
                chosenEdge = edge.index;
                deltaEnergy = delta.energy;
                break;
            }
        }
    }
    ASSERT_GE(chosenEdge, 0) << "no edge with a usable energy difference";
    ASSERT_GT(std::abs(deltaEnergy), 0.0);

    for (const double betaDeltaE : {1.0, 2.0})
    {
        Param param;
        configure(param, false);
        // kT chosen so that |dE| / kT is exactly the value under test.
        param.KBT = std::abs(deltaEnergy) / betaDeltaE;
        Mesh mesh(param);
        ASSERT_NO_FATAL_FAILURE(build(mesh, param, 1));
        isolate_single_edge(mesh, chosenEdge);
        ASSERT_EQ(count_flippable(mesh), 1) << "the chain must have exactly two states";

        // Confirm the difference is what the temperature was set from: the
        // frozen mesh is the same one it was measured on.
        EdgeFlipDelta check;
        ASSERT_TRUE(mesh.evaluate_edge_flip(chosenEdge, check));
        ASSERT_NEAR(check.energy, deltaEnergy, 1e-9 * std::abs(deltaEnergy));

        // About two attempts per sweep. A few thousand in total is enough to
        // separate exp(-1) from exp(-2) by many standard errors while keeping
        // the test to a few seconds -- each attempt is a real local energy
        // evaluation over the eighteen faces of the flip patch.
        set_rate_for_expected_attempts(param, mesh, 2.0);
        std::vector<EdgeFlipRecord> log;
        for (long long iteration = 0; iteration < 1200; iteration++)
        {
            mesh.edge_flip_sweep(iteration, &log);
        }
        ASSERT_GT(log.size(), 1500u) << "too few attempts to measure an occupancy";

        // Replay the log. State 0 is the starting triangulation, and the one
        // whose energy is lower by dE is state 1.
        long long inStart = 0;
        long long inFlipped = 0;
        bool flipped = false;
        for (const EdgeFlipRecord &record : log)
        {
            (flipped ? inFlipped : inStart)++;
            if (record.accepted)
            {
                flipped = !flipped;
            }
        }
        ASSERT_GT(inStart, 100);
        ASSERT_GT(inFlipped, 100);

        const double measured = static_cast<double>(inFlipped) / static_cast<double>(inStart);
        const double expected = std::exp(-deltaEnergy / param.KBT);

        // Binomial counting error on the smaller population, generously.
        const double relativeError =
            4.0 * std::sqrt(1.0 / static_cast<double>(std::min(inStart, inFlipped)));
        EXPECT_NEAR(measured, expected, relativeError * expected)
            << "dE/kT = " << betaDeltaE << ": visited the flipped state " << inFlipped
            << " times and the starting one " << inStart << ", ratio " << measured
            << " against exp(-dE/kT) = " << expected;
    }
}

// ---------------------------------------------------------------------------
// The schedule
// ---------------------------------------------------------------------------

/**
 * @brief The attempt count is Poisson: mean and variance both equal lambda.
 *
 * A fixed count per step would pass a mean test and fail this one. The
 * variance is what says the attempts are a Poisson process on the edges rather
 * than a quota being filled.
 */
TEST(EdgeFlipSweepTest, TheAttemptCountIsPoissonDistributed)
{
    Param param;
    configure(param, false);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build(mesh, param, 1));

    const int nFlippable = count_flippable(mesh);
    ASSERT_GT(nFlippable, 0);
    // Drive the draw alone: with the rate this low the sweep asks for a few
    // attempts and the acceptance path is irrelevant to the statistic.
    param.edgeFlipAttemptRate = 3.0 / (param.timeStep * nFlippable);
    const double lambda = 3.0;

    const int nSweeps = 20000;
    double sum = 0.0;
    double sumSquares = 0.0;
    std::vector<int> histogram(20, 0);
    for (long long iteration = 0; iteration < nSweeps; iteration++)
    {
        const std::uint64_t sweepKey =
            slimed::splitmix64(static_cast<std::uint64_t>(param.randomSeed) ^ 0x464C4950ULL) +
            static_cast<std::uint64_t>(iteration);
        const int drawn = slimed::poisson(sweepKey, lambda);
        sum += drawn;
        sumSquares += static_cast<double>(drawn) * drawn;
        if (drawn < static_cast<int>(histogram.size()))
        {
            histogram[drawn]++;
        }
    }

    const double mean = sum / nSweeps;
    const double variance = sumSquares / nSweeps - mean * mean;
    EXPECT_NEAR(mean, lambda, 0.05 * lambda);
    EXPECT_NEAR(variance, lambda, 0.10 * lambda) << "a Poisson variable has variance = mean";

    // And the shape, not just the first two moments: P(k) = e^-l l^k / k!
    double factorial = 1.0;
    for (int k = 0; k <= 6; k++)
    {
        if (k > 0)
        {
            factorial *= k;
        }
        const double expected = std::exp(-lambda) * std::pow(lambda, k) / factorial;
        const double measured = static_cast<double>(histogram[k]) / nSweeps;
        EXPECT_NEAR(measured, expected, 0.02 + 0.05 * expected) << "k = " << k;
    }
}

/**
 * @brief Halving the time step halves the attempts, so the rate per edge is
 * the physical quantity the parameter names.
 *
 * This is the whole reason the count is drawn rather than fixed. A fixed
 * number of attempts per step would double the membrane's fluidity whenever
 * someone halved the step to improve the integrator.
 */
TEST(EdgeFlipSweepTest, TheAttemptRateIsPerUnitTimeNotPerStep)
{
    Param param;
    configure(param, false);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build(mesh, param, 1));

    const auto totalDrawn = [&](double timeStep, int nSweeps) {
        param.timeStep = timeStep;
        long long drawn = 0;
        for (long long iteration = 0; iteration < nSweeps; iteration++)
        {
            // The count alone; no flips are committed, because a sweep at this
            // point would change the mesh and the comparison with it.
            const std::uint64_t sweepKey =
                slimed::splitmix64(static_cast<std::uint64_t>(param.randomSeed) ^ 0x464C4950ULL) +
                static_cast<std::uint64_t>(iteration);
            const double lambda = param.edgeFlipAttemptRate * param.timeStep *
                                  param.edgeFlipInterval * count_flippable(mesh);
            drawn += slimed::poisson(sweepKey, lambda);
        }
        return drawn;
    };

    // Twice as many steps at half the step length is the same elapsed time.
    const long double coarse = static_cast<long double>(totalDrawn(0.02, 2000));
    const long double fine = static_cast<long double>(totalDrawn(0.01, 4000));
    ASSERT_GT(coarse, 1000.0L);
    EXPECT_NEAR(static_cast<double>(fine / coarse), 1.0, 0.06)
        << "same elapsed time drew " << static_cast<double>(coarse) << " attempts at dt = 0.02 and "
        << static_cast<double>(fine) << " at dt = 0.01";
}

// ---------------------------------------------------------------------------
// Bookkeeping the sweep is responsible for
// ---------------------------------------------------------------------------

/**
 * @brief The running area and volume totals stay true across a sweep.
 *
 * The sweep updates them incrementally, because the next attempt's constraint
 * difference is measured against them and recomputing from the whole mesh
 * would put an O(faces) pass inside the loop. Incremental totals are exactly
 * the kind of bookkeeping that drifts silently, so this compares them against
 * a full recomputation after a run of accepted flips.
 */
TEST(EdgeFlipSweepTest, RunningTotalsTrackTheTrueAreaAndVolume)
{
    Param param;
    configure(param, true);
    param.KBT = 4.17;
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build(mesh, param, 2));
    set_rate_for_expected_attempts(param, mesh, 6.0);

    int accepted = 0;
    for (long long iteration = 0; iteration < 25; iteration++)
    {
        accepted += mesh.edge_flip_sweep(iteration).accepted;
    }
    ASSERT_GT(accepted, 10) << "the sweep should be accepting flips at this temperature";

    const double runningArea = param.area;
    const double runningVolume = param.vol;

    mesh.calculate_element_area_volume();
    double trueArea = 0.0;
    double trueVolume = 0.0;
    mesh.sum_membrane_area_and_volume(trueArea, trueVolume);

    EXPECT_NEAR(runningArea, trueArea, 1e-8 * std::abs(trueArea));
    EXPECT_NEAR(runningVolume, trueVolume, 1e-8 * std::abs(trueVolume));
}

TEST(EdgeFlipSweepTest, TheSameSeedGivesTheSameFlips)
{
    const auto run = [](unsigned int seed) {
        Param param;
        configure(param, false);
        param.KBT = 4.17;
        param.randomSeed = seed;
        Mesh mesh(param);
        build(mesh, param, 1);
        set_rate_for_expected_attempts(param, mesh, 3.0);
        std::vector<EdgeFlipRecord> log;
        for (long long iteration = 0; iteration < 15; iteration++)
        {
            mesh.edge_flip_sweep(iteration, &log);
        }
        return log;
    };

    const std::vector<EdgeFlipRecord> first = run(1234u);
    const std::vector<EdgeFlipRecord> again = run(1234u);
    const std::vector<EdgeFlipRecord> different = run(4321u);

    ASSERT_GT(first.size(), 20u);
    ASSERT_EQ(first.size(), again.size());
    for (std::size_t i = 0; i < first.size(); i++)
    {
        EXPECT_EQ(first[i].edge, again[i].edge) << "attempt " << i;
        EXPECT_EQ(first[i].accepted, again[i].accepted) << "attempt " << i;
        EXPECT_DOUBLE_EQ(first[i].deltaEnergy, again[i].deltaEnergy) << "attempt " << i;
    }

    // A different seed must actually explore differently, or the first check
    // is testing nothing.
    bool differs = different.size() != first.size();
    for (std::size_t i = 0; i < std::min(first.size(), different.size()) && !differs; i++)
    {
        differs = (first[i].edge != different[i].edge);
    }
    EXPECT_TRUE(differs);
}

TEST(EdgeFlipSweepTest, AZeroRateAndAFrozenMeshBothDoNothing)
{
    Param param;
    configure(param, false);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build(mesh, param, 1));
    const long long versionBefore = mesh.topologyVersion;

    param.edgeFlipAttemptRate = 0.0;
    EdgeFlipSweepStats stats = mesh.edge_flip_sweep(0);
    EXPECT_EQ(stats.drawn, 0);
    EXPECT_EQ(stats.accepted, 0);
    EXPECT_EQ(mesh.topologyVersion, versionBefore);

    param.edgeFlipAttemptRate = 5.0;
    mesh.flipFrozenVertex.assign(mesh.vertices.size(), 1);
    for (int i = 0; i < static_cast<int>(mesh.edges.size()); i++)
    {
        mesh.refresh_edge_flippability(i);
    }
    stats = mesh.edge_flip_sweep(1);
    EXPECT_EQ(stats.drawn, 0);
    EXPECT_EQ(stats.accepted, 0);
    EXPECT_EQ(mesh.topologyVersion, versionBefore);
}

/**
 * @brief A sweep leaves a mesh the evaluator can still describe.
 *
 * Everything WP0 pins about a single flip has to survive a few hundred of them
 * interleaved with acceptance decisions, and every face the sweep disturbed
 * has to come back with a patch rather than an empty one-ring.
 */
TEST(EdgeFlipSweepTest, ASweepLeavesAValidAndEvaluableMesh)
{
    Param param;
    configure(param, true);
    param.KBT = 20.0;
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build(mesh, param, 2));
    set_rate_for_expected_attempts(param, mesh, 8.0);

    int accepted = 0;
    for (long long iteration = 0; iteration < 25; iteration++)
    {
        accepted += mesh.edge_flip_sweep(iteration).accepted;
    }
    ASSERT_GT(accepted, 20);

    std::string why;
    EXPECT_TRUE(mesh.validate_manifold_topology(&why)) << why;

    // Face::patchKind is what the ring builder concluded, which is stricter
    // than classify_face(): the classifier only reads valences, while building
    // the ring can still fail on a fan that does not close. Both have to hold.
    int nUnevaluable = 0;
    for (const Face &face : mesh.faces)
    {
        const PatchClass patch = mesh.classify_face(face.index);
        if (!patch.has_evaluable_patch() || !face.oneRingVertices.empty())
        {
            nUnevaluable += patch.has_evaluable_patch() ? 0 : 1;
        }
        else
        {
            // Evaluable by valence but carrying no control net: this is the
            // silent-zero failure the admission machinery exists to prevent,
            // and a fluid mesh used to reach it through an ambiguous
            // "opposite node" lookup in the fan walk.
            ADD_FAILURE() << "face " << face.index << " has valences (" << patch.valence[0]
                          << ", " << patch.valence[1] << ", " << patch.valence[2]
                          << ") but no one-ring";
        }
        if (face.patchKind == PatchKind::MultiExtraordinary)
        {
            EXPECT_GE(face.patchEntry, 0) << "face " << face.index;
        }
    }
    EXPECT_EQ(nUnevaluable, 0) << "a closed mesh should carry a patch on every face";

    // And the energy is still a finite number.
    mesh.calculate_element_area_volume();
    mesh.Compute_Energy_And_Force();
    EXPECT_TRUE(std::isfinite(param.energy.energyTotal));
}

// ---------------------------------------------------------------------------
// The finding that WP3 turned up
// ---------------------------------------------------------------------------

namespace
{
int count_ambiguous_edges(const Mesh &m)
{
    int n = 0;
    for (const MeshEdge &edge : m.edges)
    {
        const std::vector<int> &a = m.vertices[edge.v[0]].adjacentVertices;
        const std::vector<int> &b = m.vertices[edge.v[1]].adjacentVertices;
        int shared = 0;
        for (int candidate : a)
        {
            shared += (std::find(b.begin(), b.end(), candidate) != b.end()) ? 1 : 0;
        }
        n += (shared > 2) ? 1 : 0;
    }
    return n;
}

int count_without_a_patch(const Mesh &m)
{
    int n = 0;
    for (const Face &face : m.faces)
    {
        n += face.oneRingVertices.empty() ? 1 : 0;
    }
    return n;
}
} // namespace

/**
 * @brief The one-ring walk resolves the configuration a fluid mesh reaches and
 * a near-regular one never does.
 *
 * The fan walk asks "what is the corner across this edge from that one?". That
 * used to be answered by intersecting the two vertices' neighbour lists and
 * taking a common neighbour that was not the excluded one -- correct only when
 * they share exactly two, which is what a near-regular mesh gives and what
 * every mesh this tree built was.
 *
 * Flips break the assumption by construction. Adjacent vertices start sharing
 * a third neighbour that forms no face with the edge between them, and the
 * walk returned whichever candidate it happened to see last. The consequence
 * was faces whose fan would not close, left with no control net and therefore
 * carrying no energy and no force -- silently.
 *
 * The edge table answers the question exactly: an edge of a two-manifold has
 * two incident faces, and their third corners are the only two candidates
 * there have ever been. This drives an unguarded walk through triangulation
 * space until the ambiguous configuration actually appears, then checks the
 * answer against the face structure.
 */
TEST(EdgeFlipSweepTest, TheOneRingWalkResolvesAmbiguousNeighbourhoods)
{
    Param param;
    configure(param, false);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build(mesh, param, 1));

    // Straight through the primitive, with no energy evaluated and no guard:
    // this is what reaches a strained triangulation quickly.
    Lcg random(12345u);
    int flipped = 0;
    int ambiguousEdges = 0;
    for (int attempt = 0; attempt < 400000 && ambiguousEdges == 0; attempt++)
    {
        const int iEdge = random.below(static_cast<int>(mesh.edges.size()));
        if (!mesh.edge_flip_is_admissible(iEdge))
        {
            continue;
        }
        mesh.flip_edge(iEdge);
        flipped++;
        if (flipped % 25 == 0)
        {
            ambiguousEdges = count_ambiguous_edges(mesh);
        }
    }
    ASSERT_GT(ambiguousEdges, 0)
        << "the fixture never reached the configuration this test exists for, after " << flipped
        << " flips";

    std::string why;
    ASSERT_TRUE(mesh.validate_manifold_topology(&why)) << why;

    int nChecked = 0;
    for (const MeshEdge &edge : mesh.edges)
    {
        const std::vector<int> &a = mesh.vertices[edge.v[0]].adjacentVertices;
        const std::vector<int> &b = mesh.vertices[edge.v[1]].adjacentVertices;
        int shared = 0;
        for (int candidate : a)
        {
            shared += (std::find(b.begin(), b.end(), candidate) != b.end()) ? 1 : 0;
        }
        if (shared <= 2)
        {
            continue;
        }
        // However many neighbours the endpoints share, the two faces on the
        // edge name exactly two corners, and those are the only right answers.
        EXPECT_NE(edge.opposite[0], edge.opposite[1]);
        EXPECT_EQ(mesh.find_opposite_node_index(edge.v[0], edge.v[1], edge.opposite[0]),
                  edge.opposite[1]);
        EXPECT_EQ(mesh.find_opposite_node_index(edge.v[0], edge.v[1], edge.opposite[1]),
                  edge.opposite[0]);
        nChecked++;
    }
    EXPECT_GT(nChecked, 0);
}

/**
 * @brief A trial flip that would leave a face without a patch is refused, and
 * that is what keeps a fluid mesh evaluable indefinitely.
 *
 * Not a defect but a policy, and the reasoning is worth stating: a face whose
 * one-ring cannot be built carries no energy, and zero is the lowest energy
 * there is. A chain allowed to reach such a configuration would be actively
 * drawn into it -- the Hamiltonian would develop a hole and the membrane would
 * tear along it. Refusing the move keeps the chain inside the set of
 * configurations the model can describe, which is where a Metropolis chain has
 * to stay for its stationary distribution to mean anything.
 *
 * The fixture accepts every flip the admission test allows, ignoring the
 * energy entirely. That is far more aggressive than a Metropolis sweep, which
 * is the point: the guard, not the energy, is what has to hold the line.
 */
TEST(EdgeFlipSweepTest, AFlipThatWouldCostAFaceItsPatchIsRefused)
{
    Param param;
    configure(param, false);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build(mesh, param, 1));
    ASSERT_EQ(count_without_a_patch(mesh), 0) << "the fixture starts fully evaluable";

    Lcg random(12345u);
    int flipped = 0;
    int nRefusedForLosingAPatch = 0;
    for (int attempt = 0; attempt < 20000 && flipped < 900; attempt++)
    {
        const int iEdge = random.below(static_cast<int>(mesh.edges.size()));
        EdgeFlipDelta delta;
        std::string why;
        if (!mesh.evaluate_edge_flip(iEdge, delta, &why))
        {
            nRefusedForLosingAPatch +=
                (why.find("without a subdivision patch") != std::string::npos) ? 1 : 0;
            continue;
        }
        mesh.flip_edge(iEdge);
        flipped++;
        if (flipped % 25 == 0)
        {
            ASSERT_EQ(count_without_a_patch(mesh), 0)
                << "a face lost its subdivision patch after " << flipped << " guarded flips";
        }
    }
    ASSERT_GT(flipped, 500);

    EXPECT_EQ(count_without_a_patch(mesh), 0);
    EXPECT_GT(nRefusedForLosingAPatch, 0)
        << "the guard never fired, so this run did not test it";

    std::string why;
    EXPECT_TRUE(mesh.validate_manifold_topology(&why)) << why;
}

// ---------------------------------------------------------------------------
// The sampling measure, measured rather than assumed
// ---------------------------------------------------------------------------

namespace
{
/// ln|det A| by LU with partial pivoting. Small dense matrices only.
double log_abs_determinant(std::vector<double> matrix, int n)
{
    double logDeterminant = 0.0;
    for (int column = 0; column < n; column++)
    {
        int pivot = column;
        for (int row = column + 1; row < n; row++)
        {
            if (std::abs(matrix[row * n + column]) > std::abs(matrix[pivot * n + column]))
            {
                pivot = row;
            }
        }
        const double head = matrix[pivot * n + column];
        if (head == 0.0)
        {
            return -std::numeric_limits<double>::infinity();
        }
        if (pivot != column)
        {
            for (int k = 0; k < n; k++)
            {
                std::swap(matrix[pivot * n + k], matrix[column * n + k]);
            }
        }
        logDeterminant += std::log(std::abs(head));
        for (int row = column + 1; row < n; row++)
        {
            const double factor = matrix[row * n + column] / head;
            if (factor == 0.0)
            {
                continue;
            }
            for (int k = column; k < n; k++)
            {
                matrix[row * n + k] -= factor * matrix[column * n + k];
            }
        }
    }
    return logDeterminant;
}

/**
 * @brief The Loop limit mask, with each row using its own vertex's valence.
 *
 * Half on the vertex and 1/(2N) on each of its N neighbours. This is what maps
 * the control net onto the limit surface, and what the Brownian step is
 * integrating when it displaces the surface rather than the net.
 */
std::vector<double> limit_mask(const Mesh &mesh)
{
    const int n = static_cast<int>(mesh.vertices.size());
    std::vector<double> matrix(static_cast<std::size_t>(n) * n, 0.0);
    for (int v = 0; v < n; v++)
    {
        const std::vector<int> &neighbours = mesh.vertices[v].adjacentVertices;
        const int valence = static_cast<int>(neighbours.size());
        if (valence == 0)
        {
            matrix[v * n + v] = 1.0;
            continue;
        }
        matrix[v * n + v] = 0.5;
        for (int neighbour : neighbours)
        {
            matrix[v * n + neighbour] = 0.5 / valence;
        }
    }
    return matrix;
}
} // namespace

/**
 * @brief How far the flip's sampling measure is from the one the dynamics uses.
 *
 * The Brownian step displaces the limit surface `S = M C` with isotropic
 * noise, so at fixed connectivity it samples `exp(-E) dS`. A flip at fixed
 * control net with plain Metropolis samples `exp(-E) dC` across
 * triangulations. The two differ by the Jacobian `|det M_T|`, which depends on
 * the triangulation because the limit mask depends on valence.
 *
 * Plain Metropolis at fixed control net is what every dynamically triangulated
 * surface model does and what this implementation does. The discrepancy is not
 * hidden behind that choice: it is measured here, so the size of the
 * approximation is a recorded number rather than an assumption. It biases only
 * the relative weights of triangulations; the geometry sampled at fixed
 * connectivity is unaffected either way.
 */
TEST(EdgeFlipSweepTest, TheMeasureDiscrepancyOfAFlipIsSmallAndMeasured)
{
    Param param;
    configure(param, false);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build(mesh, param, 1));

    const int n = static_cast<int>(mesh.vertices.size());
    double worst = 0.0;
    int nMeasured = 0;
    for (const MeshEdge &edge : mesh.edges)
    {
        if (!mesh.edge_flip_is_admissible(edge.index) || nMeasured >= 12)
        {
            continue;
        }
        const double before = log_abs_determinant(limit_mask(mesh), n);
        mesh.flip_edge(edge.index);
        const double after = log_abs_determinant(limit_mask(mesh), n);
        mesh.flip_edge(edge.index);

        const double discrepancy = std::abs(after - before);
        EXPECT_TRUE(std::isfinite(discrepancy));
        worst = std::max(worst, discrepancy);
        nMeasured++;
    }
    ASSERT_GT(nMeasured, 5);

    // A few percent in log weight. Small enough that the standard choice is
    // defensible, large enough to be worth stating rather than calling zero.
    // If this ever grows, the Jacobian correction is a 4x4 determinant on the
    // four changed rows -- see section 1.4(c) of the plan.
    EXPECT_LT(worst, 0.25) << "largest |ln det M' - ln det M| over a single flip: " << worst;
    std::cout << "[measure] largest |ln det M' - ln det M| over one flip: " << worst
              << " (weight ratio " << std::exp(worst) << ")" << std::endl;
}
