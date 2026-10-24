#include "test_incremental_flip_forces.hpp"

/**
 * The post-flip pass (Mesh::Compute_Energy_And_Force_After_Flips) replaces
 * the whole-mesh re-evaluation that followed every sweep with an acceptance.
 * It re-integrates only the faces whose control nets the flips changed and
 * rebuilds the vertex forces from the per-face store. The claim it has to
 * make good on is exactness: in the mode where it runs, its result and a
 * whole-mesh pass over the same mesh are the same bits. These gates measure
 * that, and the two things the mode changes on its own -- the global
 * constraint factors are applied after the sum instead of inside the
 * integrand, and the sweep records the faces it touched.
 *
 * @see docs/incremental_flip_update.md
 */

namespace
{
struct MeshFixture
{
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;
};

/// A closed icosphere: closed so the volume constraint is defined, curved so
/// the bending energy is not degenerate, and with twelve valence-5 vertices so
/// that irregular and multi-extraordinary faces are exercised from the start.
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

/// A fluid sheet with both global constraints live, so the factors the mode
/// moves out of the integrand are large and far from one.
void configure(Param &param)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.kCurv = 83.4;
    param.KBT = 4.17;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
    param.isGlobalConstraint = true;
    param.uSurf = 250.0;
    param.uVol = 30.0;
    param.area0 = 1.0;
    param.vol0 = 1.0;
    param.randomSeed = 20261024u;
    param.timeStep = 1.0;
    param.edgeFlipEnabled = true;
    param.edgeFlipMinValence = 4;
    param.edgeFlipMaxValence = 8;
    param.edgeFlipInterval = 1;
    // A few attempts per sweep on the level-2 icosphere (480 edges): enough
    // for flips to land, few enough that a sweep's trials leave most of the
    // mesh untouched, which is what the incremental pass is for.
    param.edgeFlipAttemptRate = 0.01;
}

/// Build the icosphere into @p mesh, set the reference area and volume off
/// their relaxed values so both constraints pull, and take the first pass.
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
    param.area0 = 0.9 * area;
    param.vol0 = 1.1 * volume;
    mesh.Compute_Energy_And_Force();
}

/// Everything a force pass leaves on the mesh, flattened for comparison.
struct Snapshot
{
    std::vector<double> vertexForces; ///< 15 per vertex: curvature, area, volume, regularization, total.
    std::vector<double> faceFields;   ///< 6 per face: E_bend, E_reg, E_total, elementArea, elementVolume, meanCurvature.
    std::vector<double> faceNormals;  ///< 3 per face.
    double area = 0.0;
    double vol = 0.0;
    double energyCurvature = 0.0;
    double energyArea = 0.0;
    double energyVolume = 0.0;
    double energyRegularization = 0.0;
    double energyTotal = 0.0;
};

Snapshot snapshot(const Mesh &mesh)
{
    Snapshot s;
    for (const Vertex &vertex : mesh.vertices)
    {
        const Force &f = vertex.force;
        for (int axis = 0; axis < 3; axis++)
        {
            s.vertexForces.push_back(f.forceCurvature.get(axis, 0));
        }
        for (int axis = 0; axis < 3; axis++)
        {
            s.vertexForces.push_back(f.forceArea.get(axis, 0));
        }
        for (int axis = 0; axis < 3; axis++)
        {
            s.vertexForces.push_back(f.forceVolume.get(axis, 0));
        }
        for (int axis = 0; axis < 3; axis++)
        {
            s.vertexForces.push_back(f.forceRegularization.get(axis, 0));
        }
        for (int axis = 0; axis < 3; axis++)
        {
            s.vertexForces.push_back(f.forceTotal.get(axis, 0));
        }
    }
    for (const Face &face : mesh.faces)
    {
        s.faceFields.push_back(face.energy.energyCurvature);
        s.faceFields.push_back(face.energy.energyRegularization);
        s.faceFields.push_back(face.energy.energyTotal);
        s.faceFields.push_back(face.elementArea);
        s.faceFields.push_back(face.elementVolume);
        s.faceFields.push_back(face.meanCurvature);
        for (int axis = 0; axis < 3; axis++)
        {
            s.faceNormals.push_back(face.normVector.get(axis, 0));
        }
    }
    s.area = mesh.param.area;
    s.vol = mesh.param.vol;
    s.energyCurvature = mesh.param.energy.energyCurvature;
    s.energyArea = mesh.param.energy.energyArea;
    s.energyVolume = mesh.param.energy.energyVolume;
    s.energyRegularization = mesh.param.energy.energyRegularization;
    s.energyTotal = mesh.param.energy.energyTotal;
    return s;
}

/**
 * Same bits, or -- only in an OpenMP build -- the same up to the few ulp the
 * runtime's reduction order can move a sum by. sum_membrane_area_and_volume()
 * is an OpenMP reduction, so param.area can differ between two calls on
 * identical faces, and the constraint factor carries that into every area
 * force. The serial build, which CI also runs, demands exact equality.
 */
void expect_same(double a, double b, const std::string &what)
{
#ifdef OMP
    const double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
    EXPECT_NEAR(a, b, 1e-12 * scale) << what;
#else
    EXPECT_EQ(a, b) << what;
#endif
}

void expect_same(const Snapshot &a, const Snapshot &b)
{
    ASSERT_EQ(a.vertexForces.size(), b.vertexForces.size());
    ASSERT_EQ(a.faceFields.size(), b.faceFields.size());
    ASSERT_EQ(a.faceNormals.size(), b.faceNormals.size());
    for (std::size_t i = 0; i < a.vertexForces.size(); i++)
    {
        expect_same(a.vertexForces[i], b.vertexForces[i],
                    "vertex " + std::to_string(i / 15) + " component " + std::to_string(i % 15));
    }
    for (std::size_t i = 0; i < a.faceFields.size(); i++)
    {
        expect_same(a.faceFields[i], b.faceFields[i],
                    "face " + std::to_string(i / 6) + " field " + std::to_string(i % 6));
    }
    for (std::size_t i = 0; i < a.faceNormals.size(); i++)
    {
        expect_same(a.faceNormals[i], b.faceNormals[i], "normal " + std::to_string(i));
    }
    expect_same(a.area, b.area, "param.area");
    expect_same(a.vol, b.vol, "param.vol");
    expect_same(a.energyCurvature, b.energyCurvature, "E_bend");
    expect_same(a.energyArea, b.energyArea, "E_area");
    expect_same(a.energyVolume, b.energyVolume, "E_vol");
    expect_same(a.energyRegularization, b.energyRegularization, "E_reg");
    expect_same(a.energyTotal, b.energyTotal, "E_total");
}

/// Sweep until at least @p wanted flips have landed, or give up.
int sweep_until_flips_land(Mesh &mesh, int wanted, int maxSweeps)
{
    int accepted = 0;
    for (long long iteration = 0; iteration < maxSweeps && accepted < wanted; iteration++)
    {
        accepted += mesh.edge_flip_sweep(iteration).accepted;
    }
    return accepted;
}
} // namespace

// ---------------------------------------------------------------------------
// The store and the touched set
// ---------------------------------------------------------------------------

TEST(IncrementalFlipForcesTest, TheModeIsOnForAFluidCpuRunAndOffOtherwise)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    EXPECT_TRUE(mesh.incremental_flip_forces_active());

    param.edgeFlipIncrementalForces = false;
    EXPECT_FALSE(mesh.incremental_flip_forces_active());
    param.edgeFlipIncrementalForces = true;

    param.edgeFlipEnabled = false;
    EXPECT_FALSE(mesh.incremental_flip_forces_active()) << "a frozen mesh has nothing to patch";
    param.edgeFlipEnabled = true;

    param.lipidFlowMode = "local";
    EXPECT_FALSE(mesh.incremental_flip_forces_active())
        << "a per-face tension is not a global factor";
    param.lipidFlowMode = "global";
    EXPECT_TRUE(mesh.incremental_flip_forces_active());
}

TEST(IncrementalFlipForcesTest, TheParameterRoundTripsThroughTheParser)
{
    Param param;
    EXPECT_TRUE(param.edgeFlipIncrementalForces);
    EXPECT_TRUE(import_kv_string("edgeFlipIncrementalForces", "false", param));
    EXPECT_FALSE(param.edgeFlipIncrementalForces);
    EXPECT_TRUE(import_kv_string("edgeFlipIncrementalForces", "true", param));
    EXPECT_TRUE(param.edgeFlipIncrementalForces);
}

TEST(IncrementalFlipForcesTest, APassStoresOneBlockPerFaceAndLeavesNothingTouched)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build(mesh, param, 1));

    ASSERT_GT(mesh.faces.size(), 0u);
    EXPECT_EQ(mesh.faceForceStore.size() % mesh.faces.size(), 0u);
    // Three kMaxControlPoints x 3 blocks per face.
    EXPECT_EQ(mesh.faceForceStore.size() / mesh.faces.size(), 3u * slimed::kMaxControlPoints * 3u);
    EXPECT_TRUE(mesh.flipTouchedFaces.empty());

    // The store is what the vertex forces are built from: a face's block
    // summed over the faces that reference a vertex is that vertex's force.
    std::vector<double> rebuilt(mesh.vertices.size() * 3, 0.0);
    const std::size_t stride = mesh.faceForceStore.size() / mesh.faces.size();
    for (std::size_t iFace = 0; iFace < mesh.faces.size(); iFace++)
    {
        const Face &face = mesh.faces[iFace];
        if (face.isBoundary || face.patchKind == PatchKind::Boundary)
        {
            continue;
        }
        const double *block = mesh.faceForceStore.data() + iFace * stride;
        for (std::size_t j = 0; j < face.oneRingVertices.size(); j++)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                rebuilt[face.oneRingVertices[j] * 3 + axis] += block[j * 3 + axis];
            }
        }
    }
    int nChecked = 0;
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            const double stored = mesh.vertices[i].force.forceCurvature.get(axis, 0);
            EXPECT_NEAR(rebuilt[i * 3 + axis], stored, 1e-12 * std::max(1.0, std::abs(stored)))
                << "vertex " << i << " axis " << axis;
            nChecked++;
        }
    }
    EXPECT_GT(nChecked, 0);
}

/**
 * @brief The sweep records the faces its trials touched: every face incident
 * to a corner of an attempted flip's quadrilateral, accepted or not.
 *
 * Rejected attempts matter too. A trial flips the edge and flips it back, and
 * that rebuilds the one-rings of the faces around it; the mesh it leaves is
 * the same mesh, but a rebuilt one-ring can list the same control points from
 * a different starting corner, and the kernel's sum over them then rounds
 * differently. Those faces have to be re-integrated for the incremental pass
 * to match a whole-mesh pass to the bit.
 */
TEST(IncrementalFlipForcesTest, TheSweepRecordsTheFacesItsAttemptsTouched)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build(mesh, param, 2));

    std::vector<EdgeFlipRecord> log;
    EdgeFlipSweepStats stats;
    for (long long iteration = 0; iteration < 50 && stats.accepted == 0; iteration++)
    {
        stats = mesh.edge_flip_sweep(iteration, &log);
    }
    ASSERT_GT(stats.accepted, 0) << "no flip landed, so the gate tested nothing";

    // Every face incident to a corner of an attempted flip's quadrilateral is
    // in the set, once, and the marks agree with the list.
    ASSERT_FALSE(mesh.flipTouchedFaces.empty());
    std::vector<char> seen(mesh.faces.size(), 0);
    for (int iFace : mesh.flipTouchedFaces)
    {
        ASSERT_GE(iFace, 0);
        ASSERT_LT(iFace, static_cast<int>(mesh.faces.size()));
        EXPECT_EQ(seen[iFace], 0) << "face " << iFace << " listed twice";
        seen[iFace] = 1;
        EXPECT_EQ(mesh.flipTouchedMark[iFace], 1);
    }
    int nAttempts = 0;
    for (const EdgeFlipRecord &record : log)
    {
        nAttempts++;
        for (int corner : record.vertex)
        {
            for (int iFace : mesh.vertices[corner].adjacentFaces)
            {
                EXPECT_EQ(seen[iFace], 1) << "face " << iFace << " at corner " << corner
                                          << " of an attempted flip";
            }
        }
    }
    EXPECT_GT(nAttempts, stats.accepted) << "every attempt landed; the rejected case went untested";
    // And not the whole mesh, or there would be nothing to save.
    EXPECT_LT(mesh.flipTouchedFaces.size(), mesh.faces.size());

    mesh.Compute_Energy_And_Force_After_Flips();
    EXPECT_TRUE(mesh.flipTouchedFaces.empty());
    for (char mark : mesh.flipTouchedMark)
    {
        EXPECT_EQ(mark, 0);
    }
}

// ---------------------------------------------------------------------------
// Exactness against the whole-mesh pass
// ---------------------------------------------------------------------------

/**
 * @brief The gate. After a sweep that landed flips, the incremental pass and
 * a whole-mesh pass over the same mesh leave the same forces, energies,
 * areas, normals and totals -- bit for bit in a serial build.
 */
TEST(IncrementalFlipForcesTest, ThePostFlipPassMatchesAWholeMeshPass)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build(mesh, param, 2));
    const Snapshot beforeFlips = snapshot(mesh);

    const int accepted = sweep_until_flips_land(mesh, 2, 50);
    ASSERT_GE(accepted, 2) << "too few flips landed for the gate to mean much";
    const std::size_t nTouched = mesh.flipTouchedFaces.size();
    ASSERT_GT(nTouched, 0u);
    ASSERT_LT(nTouched, mesh.faces.size()) << "every face was touched; nothing was kept";

    mesh.Compute_Energy_And_Force_After_Flips();
    const Snapshot incremental = snapshot(mesh);
    EXPECT_TRUE(mesh.flipTouchedFaces.empty());

    // The flips changed the mesh, and the pass noticed: some face energies
    // moved from what the pre-flip pass left.
    int nFaceEnergiesMoved = 0;
    for (std::size_t iFace = 0; iFace < mesh.faces.size(); iFace++)
    {
        if (beforeFlips.faceFields[iFace * 6] != incremental.faceFields[iFace * 6])
        {
            nFaceEnergiesMoved++;
        }
    }
    EXPECT_GT(nFaceEnergiesMoved, 0) << "no face energy moved: the pass recomputed nothing";
    EXPECT_LE(nFaceEnergiesMoved, static_cast<int>(nTouched))
        << "a face outside the touched set changed";

    mesh.Compute_Energy_And_Force();
    const Snapshot whole = snapshot(mesh);
    expect_same(incremental, whole);
}

/**
 * @brief Repeated sweeps and incremental passes do not drift from the
 * whole-mesh pass: the store is rewritten, never incremented, so there is no
 * running sum to accumulate rounding.
 */
TEST(IncrementalFlipForcesTest, ManyIncrementalPassesInARowStayExact)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build(mesh, param, 2));

    int totalAccepted = 0;
    int incrementalPasses = 0;
    for (long long iteration = 0; iteration < 40; iteration++)
    {
        const EdgeFlipSweepStats stats = mesh.edge_flip_sweep(iteration);
        totalAccepted += stats.accepted;
        if (stats.accepted > 0)
        {
            mesh.Compute_Energy_And_Force_After_Flips();
            incrementalPasses++;
        }
    }
    ASSERT_GT(incrementalPasses, 5);
    ASSERT_GT(totalAccepted, incrementalPasses);

    const Snapshot incremental = snapshot(mesh);
    mesh.Compute_Energy_And_Force();
    expect_same(incremental, snapshot(mesh));

    std::string why;
    EXPECT_TRUE(mesh.validate_manifold_topology(&why)) << why;
}

/**
 * @brief Without a sweep there is nothing to patch, and the routine is the
 * whole-mesh pass; with the parameter off it is the whole-mesh pass too, with
 * the kernel scaling the constraint forces itself.
 */
TEST(IncrementalFlipForcesTest, FallsBackToTheWholeMeshPassWhenThereIsNothingToPatch)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build(mesh, param, 1));

    mesh.Compute_Energy_And_Force_After_Flips(); // nothing touched
    const Snapshot a = snapshot(mesh);
    mesh.Compute_Energy_And_Force();
    expect_same(a, snapshot(mesh));

    param.edgeFlipIncrementalForces = false;
    ASSERT_GT(sweep_until_flips_land(mesh, 1, 50), 0);
    EXPECT_FALSE(mesh.flipTouchedFaces.empty()) << "the sweep records faces either way";
    mesh.Compute_Energy_And_Force_After_Flips(); // falls back
    EXPECT_TRUE(mesh.flipTouchedFaces.empty());
    const Snapshot b = snapshot(mesh);
    mesh.Compute_Energy_And_Force();
    expect_same(b, snapshot(mesh));
}

// ---------------------------------------------------------------------------
// What the mode changes on its own
// ---------------------------------------------------------------------------

/**
 * @brief Applying the global constraint factors after the sum instead of
 * inside the integrand is a reassociation: the same terms in a different
 * order. The bending force and every energy are untouched by it, and the
 * area and volume forces move by rounding only.
 */
TEST(IncrementalFlipForcesTest, FactoringTheConstraintScalesIsReassociationOnly)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build(mesh, param, 2));
    ASSERT_TRUE(mesh.incremental_flip_forces_active());
    const Snapshot factored = snapshot(mesh);

    param.edgeFlipIncrementalForces = false;
    ASSERT_FALSE(mesh.incremental_flip_forces_active());
    mesh.Compute_Energy_And_Force();
    const Snapshot inKernel = snapshot(mesh);

    ASSERT_EQ(factored.vertexForces.size(), inKernel.vertexForces.size());
    double maxAreaForce = 0.0;
    double maxVolumeForce = 0.0;
    double maxAreaDifference = 0.0;
    double maxVolumeDifference = 0.0;
    for (std::size_t v = 0; v < mesh.vertices.size(); v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            const std::size_t base = v * 15;
            // Bending: the factor never entered it.
            expect_same(factored.vertexForces[base + axis], inKernel.vertexForces[base + axis],
                        "bending force, vertex " + std::to_string(v));
            maxAreaForce = std::max(maxAreaForce, std::abs(inKernel.vertexForces[base + 3 + axis]));
            maxVolumeForce =
                std::max(maxVolumeForce, std::abs(inKernel.vertexForces[base + 6 + axis]));
            maxAreaDifference =
                std::max(maxAreaDifference, std::abs(factored.vertexForces[base + 3 + axis] -
                                                     inKernel.vertexForces[base + 3 + axis]));
            maxVolumeDifference =
                std::max(maxVolumeDifference, std::abs(factored.vertexForces[base + 6 + axis] -
                                                       inKernel.vertexForces[base + 6 + axis]));
        }
    }
    // Both constraints pull: the forces compared are not trivially zero.
    EXPECT_GT(maxAreaForce, 1.0);
    EXPECT_GT(maxVolumeForce, 1.0);
    EXPECT_LT(maxAreaDifference, 1e-11 * maxAreaForce);
    EXPECT_LT(maxVolumeDifference, 1e-11 * maxVolumeForce);

    // Energies, areas, volumes and normals never saw the factor.
    for (std::size_t i = 0; i < factored.faceFields.size(); i++)
    {
        expect_same(factored.faceFields[i], inKernel.faceFields[i], "face field " + std::to_string(i));
    }
    expect_same(factored.area, inKernel.area, "param.area");
    expect_same(factored.vol, inKernel.vol, "param.vol");
    expect_same(factored.energyTotal, inKernel.energyTotal, "E_total");
}

// ---------------------------------------------------------------------------
// In the dynamics loop
// ---------------------------------------------------------------------------

namespace
{
void configure_flat(Param &param)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Periodic;
    param.sideX = 60.0;
    param.sideY = 60.0;
    param.lFace = 5.0;
    param.uVol = 0.0;
    param.kCurv = 83.4;
    param.KBT = 4.17;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
    param.randomSeed = 20260910u;
    param.timeStep = 1.0e-4;
    param.surfaceSolver = "iterative";
    param.inPlaneDynamicsEnabled = true;
    param.edgeSpringEnabled = true;
    param.edgeFlipEnabled = true;
    param.edgeFlipInterval = 1;
    param.edgeFlipMinValence = 4;
    param.edgeFlipMaxValence = 8;
    param.maxIterations = 12;
}

/// The step loop of run_dynamics_flat() with flips on, as the driver runs it,
/// returning the final control net.
std::vector<double> run_short_fluid_trajectory(bool incremental, int *stepsWithAFlip)
{
    Param param;
    configure_flat(param);
    param.edgeFlipIncrementalForces = incremental;

    DynamicMesh mesh(param);
    mesh.setup_flat();
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.Compute_Energy_And_Force();
    mesh.update_vertices_mat_with_vector();
    mesh.apply_mesh_to_surface();

    int flippable = 0;
    for (const MeshEdge &edge : mesh.edges)
    {
        flippable += edge.flippable ? 1 : 0;
    }
    EXPECT_GT(flippable, 0);
    param.edgeFlipAttemptRate = 6.0 / (param.timeStep * std::max(1, flippable));

    Record record(param.maxIterations);
    DynamicModel model(mesh, record);
    *stepsWithAFlip = 0;
    for (model.iteration = 0; model.iteration < param.maxIterations; model.iteration++)
    {
        mesh.apply_mesh_to_surface();
        model.next_step();
        mesh.apply_surface_to_mesh();
        mesh.postprocess_ghost_periodic();
        mesh.update_vertices_vector_with_mat();
        mesh.Compute_Energy_And_Force();

        const EdgeFlipSweepStats stats = mesh.edge_flip_sweep(model.iteration);
        if (stats.accepted > 0)
        {
            (*stepsWithAFlip)++;
            if (incremental)
            {
                mesh.Compute_Energy_And_Force_After_Flips();
            }
            else
            {
                mesh.Compute_Energy_And_Force();
            }
        }
    }

    std::vector<double> coords;
    for (const Vertex &vertex : mesh.vertices)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            coords.push_back(vertex.coord.get(axis, 0));
        }
    }
    return coords;
}
} // namespace

/**
 * @brief On the periodic flat sheet the dynamics actually runs, a short fluid
 * trajectory with the incremental pass is the trajectory with the whole-mesh
 * pass, up to the rounding the factored constraint scale introduces.
 */
TEST(IncrementalFlipForcesTest, AShortFluidTrajectoryIsUnchangedUpToRounding)
{
    int flipsIncremental = 0;
    int flipsWhole = 0;
    const std::vector<double> incremental = run_short_fluid_trajectory(true, &flipsIncremental);
    const std::vector<double> whole = run_short_fluid_trajectory(false, &flipsWhole);
    ASSERT_FALSE(incremental.empty());
    ASSERT_EQ(incremental.size(), whole.size());
    ASSERT_GT(flipsIncremental, 0) << "no step accepted a flip, so the pass never ran";
    EXPECT_EQ(flipsIncremental, flipsWhole);

    double maxDifference = 0.0;
    for (std::size_t i = 0; i < incremental.size(); i++)
    {
        maxDifference = std::max(maxDifference, std::abs(incremental[i] - whole[i]));
    }
    EXPECT_LT(maxDifference, 1e-9) << "nm, after " << flipsIncremental << " steps with a flip";
}
