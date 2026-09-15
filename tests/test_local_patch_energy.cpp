#include "test_local_patch_energy.hpp"

/**
 * A Metropolis sweep asks "what would this move cost?" thousands of times per
 * run. Answering with a whole-mesh energy evaluation would make the sweep cost
 * far more than the dynamics it is interleaved with, and the existing
 * single-vertex Metropolis move in Energy_minimization.cpp does exactly that --
 * it snapshots every coordinate and runs a full Compute_Energy_And_Force() per
 * trial.
 *
 * It does not have to. The energy of a face is a functional of its control
 * net, so a flip can only change the energy of a face incident to one of the
 * four vertices it touches. The claim is that everything else cancels
 * *exactly*, and the test below is that claim: the local difference against
 * the difference of two full evaluations, with both global constraints
 * switched on so the terms that cannot be summed per face are exercised too.
 *
 * @see docs/edge_flip_plan.md work package 2
 */

namespace
{
struct MeshFixture
{
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;
};

/**
 * @brief A closed icosphere: an icosahedron subdivided `level` times with every
 * new vertex pushed out onto the sphere.
 *
 * Closed, because the volume constraint is only defined on a closed surface --
 * validate_volume_constraint_topology() refuses it on anything else, and
 * rightly, since the signed volume of an open sheet is not even independent of
 * where the origin sits. The volume term is one of the two that cannot be
 * differenced face by face, so it has to be in the fixture.
 *
 * At every level the twelve original vertices stay at valence 5 and everything
 * else is 6, and no two valence-5 vertices are adjacent, so the mesh is
 * admissible before any flip.
 */
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

/**
 * @brief A closed icosphere with both global constraints active and a handful
 * of edges flipped.
 *
 * The constraints matter: they are the terms that are quadratic in the
 * mesh-wide totals and so cannot be differenced face by face. With them off,
 * the test would pass on a version of the code that ignored them entirely.
 * The flips matter too, because they put multi-extraordinary faces in the
 * neighbourhoods being measured.
 */
void build_fixture(Mesh &mesh, Param &param, int level, int nFlips)
{
    const MeshFixture fixture = build_icosphere(level, 12.0);
    mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();

    // Reference area and volume deliberately off the current values, so that
    // (A - area0) and (V - vol0) are both large: the cross term in the
    // constraint difference is what a naive dE would drop.
    mesh.calculate_element_area_volume();
    double area = 0.0;
    double volume = 0.0;
    mesh.sum_membrane_area_and_volume(area, volume);
    param.area0 = 0.85 * area;
    param.vol0 = 1.15 * volume;

    Lcg random(4242u);
    int flipped = 0;
    for (int attempt = 0; attempt < 6000 && flipped < nFlips; attempt++)
    {
        const int iEdge = random.below(static_cast<int>(mesh.edges.size()));
        if (!mesh.edge_flip_is_admissible(iEdge))
        {
            continue;
        }
        mesh.flip_edge(iEdge);
        flipped++;
    }
    ASSERT_EQ(flipped, nFlips);
}

/// Everything the fixture needs set before the Mesh is constructed.
void configure(Param &param)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.isGlobalConstraint = true;
    param.uSurf = 250.0;
    param.uVol = 30.0;
    param.kCurv = 83.4;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
    // A closed surface, so the volume constraint is well defined; the
    // references are overwritten from the built mesh in build_fixture().
    param.area0 = 1.0;
    param.vol0 = 1.0;
}

/// The number Metropolis compares against kT: every term of the Hamiltonian.
double total_energy(Mesh &mesh)
{
    mesh.calculate_element_area_volume();
    mesh.Compute_Energy_And_Force();
    return mesh.param.energy.energyTotal;
}
} // namespace

// ---------------------------------------------------------------------------
// The pieces
// ---------------------------------------------------------------------------

/**
 * @brief The standalone regularization term agrees with the one the force pass
 * computes, face by face.
 *
 * The energy is needed without the force, so it is written twice. This is what
 * keeps the two from drifting: a change to one that is not made to the other
 * shows up here rather than as a wrong acceptance rate months later.
 */
TEST(LocalPatchEnergyTest, StandaloneRegularizationMatchesTheForcePass)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_fixture(mesh, param, 2, 4));
    mesh.calculate_element_area_volume();
    mesh.Compute_Energy_And_Force();

    for (int iFace = 0; iFace < static_cast<int>(mesh.faces.size()); iFace++)
    {
        const double standalone = mesh.face_regularization_energy(iFace);
        const double fromForcePass = mesh.faces[iFace].energy.energyRegularization;
        EXPECT_NEAR(standalone, fromForcePass, 1e-12 * std::max(1.0, std::abs(fromForcePass)))
            << "face " << iFace;
    }
}

/**
 * @brief Summing the subset evaluator over every face reproduces the
 * whole-mesh totals.
 *
 * If this drifts, every local difference is wrong by the same amount and the
 * flip test below would still pass, because the error would cancel. So it is
 * worth pinning on its own.
 */
TEST(LocalPatchEnergyTest, TheSubsetEvaluatorAgreesWithTheWholeMeshPass)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_fixture(mesh, param, 2, 4));
    mesh.calculate_element_area_volume();
    mesh.Compute_Energy_And_Force();

    double bending = 0.0;
    double regularization = 0.0;
    for (const Face &face : mesh.faces)
    {
        bending += face.energy.energyCurvature;
        regularization += face.energy.energyRegularization;
    }
    double area = 0.0;
    double volume = 0.0;
    mesh.sum_membrane_area_and_volume(area, volume);

    std::vector<int> everyFace(mesh.faces.size());
    for (int i = 0; i < static_cast<int>(everyFace.size()); i++)
    {
        everyFace[i] = i;
    }
    const FaceSubsetEnergy subset = mesh.evaluate_face_subset(everyFace);

    EXPECT_NEAR(subset.bending, bending, 1e-10 * std::max(1.0, std::abs(bending)));
    EXPECT_NEAR(subset.regularization, regularization,
                1e-10 * std::max(1.0, std::abs(regularization)));
    EXPECT_NEAR(subset.area, area, 1e-10 * std::max(1.0, std::abs(area)));
    EXPECT_NEAR(subset.volume, volume, 1e-10 * std::max(1.0, std::abs(volume)));
}

// ---------------------------------------------------------------------------
// The gate
// ---------------------------------------------------------------------------

/**
 * @brief The local energy difference equals the difference of two full
 * evaluations.
 *
 * This is work package 2's gate, and the reason the sweep can afford to run.
 * Both global constraints are on, with reference values well away from the
 * current area and volume, so the cross term `2 (X - X0) dX` carries real
 * weight -- a version that dropped it would fail here and nowhere else.
 */
TEST(LocalPatchEnergyTest, LocalDeltaMatchesTheDifferenceOfTwoFullEvaluations)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_fixture(mesh, param, 2, 5));

    Lcg random(90210u);
    int nTested = 0;
    for (int attempt = 0; attempt < 400 && nTested < 25; attempt++)
    {
        const int iEdge = random.below(static_cast<int>(mesh.edges.size()));

        // The totals the constraint difference is measured against have to be
        // current, which is what the last full evaluation leaves behind.
        const double before = total_energy(mesh);

        EdgeFlipDelta delta;
        if (!mesh.evaluate_edge_flip(iEdge, delta))
        {
            continue;
        }

        // evaluate_edge_flip() leaves the mesh as it found it, so this is the
        // same configuration `before` was measured on.
        const double restored = total_energy(mesh);
        ASSERT_NEAR(restored, before, 1e-9 * std::max(1.0, std::abs(before)))
            << "evaluating a flip must not change the mesh";

        mesh.flip_edge(iEdge);
        const double after = total_energy(mesh);
        mesh.flip_edge(iEdge);

        const double expected = after - before;
        EXPECT_NEAR(delta.energy, expected, 1e-8 * std::max(1.0, std::abs(expected)))
            << "edge " << iEdge << ": local " << delta.energy << " vs global " << expected
            << " (bending " << delta.bending << ", regularization " << delta.regularization
            << ", area " << delta.areaConstraint << ", volume " << delta.volumeConstraint << ")";
        nTested++;
    }
    EXPECT_GE(nTested, 25) << "not enough admissible flips were found to test";
}

/**
 * @brief The constraint terms are actually being exercised.
 *
 * A guard on the test above rather than on the code: if the fixture happened
 * to make the constraint contributions negligible, the gate would pass without
 * ever testing the part that is hardest to get right.
 */
TEST(LocalPatchEnergyTest, TheGlobalConstraintTermsCarryRealWeight)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_fixture(mesh, param, 2, 5));
    total_energy(mesh);

    Lcg random(555u);
    double largestConstraintShare = 0.0;
    int nTested = 0;
    for (int attempt = 0; attempt < 400 && nTested < 10; attempt++)
    {
        EdgeFlipDelta delta;
        if (!mesh.evaluate_edge_flip(random.below(static_cast<int>(mesh.edges.size())), delta))
        {
            continue;
        }
        const double constraint = std::abs(delta.areaConstraint) + std::abs(delta.volumeConstraint);
        const double total = std::abs(delta.energy);
        if (total > 0.0)
        {
            largestConstraintShare = std::max(largestConstraintShare, constraint / total);
        }
        nTested++;
    }
    ASSERT_GT(nTested, 0);
    EXPECT_GT(largestConstraintShare, 0.01)
        << "the constraint terms contribute nothing measurable, so the gate above is not "
           "testing them";
}

/**
 * @brief Evaluating a flip leaves the mesh evaluable and unchanged.
 *
 * The trial flips and flips back, which exchanges the two incident face labels
 * (see Mesh::flip_edge). Nothing observable may depend on that: the
 * triangulation, every valence, and the total energy must all come back.
 */
TEST(LocalPatchEnergyTest, EvaluatingAFlipRestoresTheMesh)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_fixture(mesh, param, 2, 4));
    const double before = total_energy(mesh);

    std::vector<int> valenceBefore;
    for (const Vertex &vertex : mesh.vertices)
    {
        valenceBefore.push_back(static_cast<int>(vertex.adjacentVertices.size()));
    }

    Lcg random(31337u);
    int nEvaluated = 0;
    for (int attempt = 0; attempt < 2000 && nEvaluated < 60; attempt++)
    {
        EdgeFlipDelta delta;
        if (mesh.evaluate_edge_flip(random.below(static_cast<int>(mesh.edges.size())), delta))
        {
            nEvaluated++;
        }
    }
    EXPECT_GT(nEvaluated, 30);

    std::string why;
    EXPECT_TRUE(mesh.validate_manifold_topology(&why)) << why;
    for (int i = 0; i < static_cast<int>(mesh.vertices.size()); i++)
    {
        EXPECT_EQ(static_cast<int>(mesh.vertices[i].adjacentVertices.size()), valenceBefore[i])
            << "vertex " << i;
    }
    const double after = total_energy(mesh);
    EXPECT_NEAR(after, before, 1e-9 * std::max(1.0, std::abs(before)));
}

/**
 * @brief An inadmissible edge is reported, not evaluated, and the mesh is not
 * touched on the way out.
 */
TEST(LocalPatchEnergyTest, AnInadmissibleFlipIsRefusedWithoutTouchingTheMesh)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_fixture(mesh, param, 1, 0));
    const double before = total_energy(mesh);
    const long long versionBefore = mesh.topologyVersion;

    // A closed surface has no boundary edge, so freeze a vertex instead: a
    // periodic duplicate is not an independent coordinate, and DynamicMesh
    // marks those the same way.
    ASSERT_FALSE(mesh.edges.empty());
    const int frozenEdge = 0;
    mesh.flipFrozenVertex.assign(mesh.vertices.size(), 0);
    mesh.flipFrozenVertex[mesh.edges[frozenEdge].v[0]] = 1;
    for (int iEdge = 0; iEdge < static_cast<int>(mesh.edges.size()); iEdge++)
    {
        mesh.refresh_edge_flippability(iEdge);
    }

    EdgeFlipDelta delta;
    std::string why;
    EXPECT_FALSE(mesh.evaluate_edge_flip(frozenEdge, delta, &why));
    EXPECT_FALSE(why.empty());
    EXPECT_EQ(delta.energy, 0.0);
    EXPECT_EQ(mesh.topologyVersion, versionBefore)
        << "a refused flip must not touch the connectivity";
    EXPECT_NEAR(total_energy(mesh), before, 1e-12 * std::max(1.0, std::abs(before)));
}
