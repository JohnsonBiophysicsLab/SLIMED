#include "test_fluid_run_io.hpp"

/**
 * WP5 is the wiring: the parts that turn the flip move into something a run
 * can actually use. Three of them carry a failure mode that no amount of
 * correct physics upstream would catch.
 *
 * **The connectivity has to be written beside the coordinates.** A run with
 * flips has no single face list, and `face.csv` -- written once at setup --
 * describes the mesh the run started with. Pair it with a later coordinate
 * frame and every analysis reads the wrong triangles, silently: the file
 * parses, the triangle count is right, and the spectrum that comes out is
 * wrong by an amount nobody can see.
 *
 * **A restart has to carry the connectivity too.** Same failure, one step
 * worse: the run continues, integrating a surface whose energy is computed
 * from triangles that are not the ones it left off with.
 *
 * **The depth scale has to be a fidelity knob and not a correctness break.**
 * A fluid mesh is mostly irregular, and an irregular face costs `3D` samples
 * per extraordinary corner. The depth is the largest lever on that cost, so
 * it has to be adjustable -- but only if turning it down degrades the answer
 * gracefully rather than producing a different one.
 *
 * @see docs/edge_flip_plan.md work package 5
 */

namespace
{

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

struct MeshFixture
{
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;
};

/// A closed icosphere: closed so the volume constraint is defined, curved so
/// the bending energy is not degenerate, and irregular at the twelve original
/// vertices so the depth scale has something to act on.
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

void configure_sphere(Param &param)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.kCurv = 83.4;
    param.KBT = 4.17;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
    param.isGlobalConstraint = false;
    param.uSurf = 0.0;
    param.uVol = 0.0;
    param.randomSeed = 20260910u;
    param.timeStep = 1.0;
}

void build_sphere(Mesh &mesh, int level)
{
    const MeshFixture fixture = build_icosphere(level, 12.0);
    mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.Compute_Energy_And_Force();
}

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
}

/// Face corner lists as the mesh currently holds them.
std::vector<std::array<int, 3>> corners_of(const Mesh &mesh)
{
    std::vector<std::array<int, 3>> corners(mesh.faces.size());
    for (int i = 0; i < static_cast<int>(mesh.faces.size()); i++)
    {
        for (int k = 0; k < 3; k++)
        {
            corners[i][k] = mesh.faces[i].adjacentVertices[k];
        }
    }
    return corners;
}

/// Read back a face CSV in the layout write_faces_csv() and the face frames use.
std::vector<std::array<int, 3>> read_face_csv(const std::string &path)
{
    std::vector<std::array<int, 3>> corners;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty())
        {
            continue;
        }
        std::array<int, 3> triple = {{-1, -1, -1}};
        std::istringstream fields(line);
        std::string field;
        for (int k = 0; k < 3 && std::getline(fields, field, ','); k++)
        {
            triple[k] = std::stoi(field);
        }
        corners.push_back(triple);
    }
    return corners;
}

/**
 * @brief The mesh's triangles as a set, independent of how they are labelled.
 *
 * Flipping an edge twice returns the same triangulation but exchanges the two
 * faces' slots and rotates their corner lists -- intrinsic to the move, not a
 * defect, and unobservable because the admission test refuses a flip across a
 * spontaneous-curvature boundary, so the two faces that swap carry identical
 * attributes. A comparison that wants to ask whether the *mesh* came back has
 * to ask it this way. Winding is preserved: each triple is rotated, never
 * reordered.
 */
std::vector<std::array<int, 3>> triangle_set_of(const Mesh &mesh)
{
    std::vector<std::array<int, 3>> triangles;
    triangles.reserve(mesh.faces.size());
    for (const Face &face : mesh.faces)
    {
        std::array<int, 3> t = {{face.adjacentVertices[0], face.adjacentVertices[1],
                                 face.adjacentVertices[2]}};
        const int smallest = (t[0] <= t[1] && t[0] <= t[2]) ? 0 : ((t[1] <= t[2]) ? 1 : 2);
        triangles.push_back({{t[smallest], t[(smallest + 1) % 3], t[(smallest + 2) % 3]}});
    }
    std::sort(triangles.begin(), triangles.end());
    return triangles;
}

/// The first flippable edge, so a test can change the connectivity by one move.
int first_flippable_edge(const Mesh &mesh)
{
    for (const MeshEdge &edge : mesh.edges)
    {
        if (edge.flippable)
        {
            return edge.index;
        }
    }
    return -1;
}

/// A scratch path under the current directory, removed by the test that made it.
std::string scratch(const std::string &name)
{
    return "wp5_" + name;
}

int total_children(const IrregularPatchRowTable &table)
{
    int total = 0;
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        total += table.depth_for(valence);
    }
    return total;
}

} // namespace

// ---------------------------------------------------------------------------
// The irregular-patch depth scale
// ---------------------------------------------------------------------------

/**
 * @brief The default and an explicit 1.0 are the same table.
 *
 * The scale exists to be turned down, and a run that does not set it must get
 * exactly what it got before the knob existed -- not something that rounds to
 * the same thing.
 */
TEST(IrregularDepthScaleTest, AScaleOfOneIsTheUnscaledTable)
{
    Param param;
    configure_sphere(param);
    Mesh mesh(param);

    IrregularPatchRowTable scaled;
    scaled.build(param.shapeFunctions, kDefaultIrregularDepth, DepthPolicy::PerValence, 1.0);

    EXPECT_EQ(scaled.depth(), mesh.irregularRows.depth());
    EXPECT_EQ(scaled.memory_bytes(), mesh.irregularRows.memory_bytes());
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        EXPECT_EQ(scaled.depth_for(valence), mesh.irregularRows.depth_for(valence))
            << "valence " << valence;
        EXPECT_EQ(scaled.depth_for(valence), recommended_irregular_depth(valence))
            << "valence " << valence;
    }
}

/**
 * @brief Turning the scale down shrinks what is built as well as what is used.
 *
 * Building twelve levels and consuming six would give the cheaper evaluation
 * and none of the cheaper setup, which on the GPU path is also the upload.
 */
TEST(IrregularDepthScaleTest, AScaleBelowOneShrinksTheBuiltAndConsumedDepth)
{
    Param param;
    configure_sphere(param);
    Mesh reference(param);

    IrregularPatchRowTable half;
    half.build(param.shapeFunctions, kDefaultIrregularDepth, DepthPolicy::PerValence, 0.5);

    // The deepest valence is 4 at 12 levels, so half of it sets the build.
    EXPECT_EQ(half.depth(), 6);
    EXPECT_LT(half.memory_bytes(), reference.irregularRows.memory_bytes());

    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        if (valence == 6)
        {
            continue; // regular; never reads the table
        }
        EXPECT_EQ(half.depth_for(valence),
                  scaled_irregular_depth(recommended_irregular_depth(valence), 0.5))
            << "valence " << valence;
        EXPECT_LT(half.depth_for(valence), reference.irregularRows.depth_for(valence))
            << "valence " << valence;
    }
    EXPECT_LT(total_children(half), total_children(reference.irregularRows));
}

/**
 * @brief A scale small enough to ask for no levels still gets one.
 *
 * Depth 0 is the extraordinary corner itself, with no regular child to
 * evaluate. A face there carries no energy, and the Metropolis sweep would be
 * drawn straight into the configurations that have none.
 */
TEST(IrregularDepthScaleTest, TheDepthNeverFallsBelowOneLevel)
{
    EXPECT_EQ(scaled_irregular_depth(12, 0.001), 1);
    EXPECT_EQ(scaled_irregular_depth(7, 0.01), 1);
    // Rounded, not truncated: 7 * 0.5 = 3.5 is nearer 4 than 3.
    EXPECT_EQ(scaled_irregular_depth(7, 0.5), 4);
    EXPECT_EQ(scaled_irregular_depth(12, 0.5), 6);
    // The identity is exact, not a rounding coincidence.
    EXPECT_EQ(scaled_irregular_depth(12, 1.0), 12);
}

/// A non-positive scale would ask for a table with no levels in it.
TEST(IrregularDepthScaleTest, ANonPositiveScaleIsRejected)
{
    Param param;
    configure_sphere(param);
    Mesh mesh(param);

    IrregularPatchRowTable table;
    EXPECT_THROW(table.build(param.shapeFunctions, kDefaultIrregularDepth,
                             DepthPolicy::PerValence, 0.0),
                 std::invalid_argument);
    EXPECT_THROW(table.build(param.shapeFunctions, kDefaultIrregularDepth,
                             DepthPolicy::PerValence, -1.0),
                 std::invalid_argument);
}

/**
 * @brief The uniform policy ignores the scale.
 *
 * That policy exists so the convergence study can sweep the depth itself. A
 * second knob multiplying the one under study would only confuse it.
 */
TEST(IrregularDepthScaleTest, TheUniformPolicyIgnoresTheScale)
{
    Param param;
    configure_sphere(param);
    Mesh mesh(param);

    IrregularPatchRowTable table;
    table.build(param.shapeFunctions, 5, DepthPolicy::Uniform, 0.25);
    EXPECT_EQ(table.depth(), 5);
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        EXPECT_EQ(table.depth_for(valence), 5) << "valence " << valence;
    }
}

/**
 * @brief The scale degrades the energy gracefully rather than changing it.
 *
 * A cheaper patch is a different patch, so the energies must differ -- if they
 * did not, the depth would not be doing anything and the expensive levels
 * would be waste. But the difference has to be a truncation, small next to the
 * energy itself, or the knob is not a fidelity trade at all.
 */
TEST(IrregularDepthScaleTest, ACheaperDepthPerturbsTheEnergyWithoutChangingIt)
{
    Param full;
    configure_sphere(full);
    Mesh fullMesh(full);
    build_sphere(fullMesh, 1);

    Param cheap;
    configure_sphere(cheap);
    cheap.irregularPatchDepthScale = 0.5;
    Mesh cheapMesh(cheap);
    build_sphere(cheapMesh, 1);

    ASSERT_EQ(cheapMesh.irregularRows.depth(), 6);
    ASSERT_EQ(fullMesh.irregularRows.depth(), kDefaultIrregularDepth);

    const double fullEnergy = fullMesh.param.energy.energyTotal;
    const double cheapEnergy = cheapMesh.param.energy.energyTotal;
    ASSERT_GT(std::abs(fullEnergy), 1.0) << "the fixture has to carry some bending energy";

    EXPECT_NE(fullEnergy, cheapEnergy)
        << "halving the depth changed nothing, so the depth is not being consumed";
    const double relative = std::abs(cheapEnergy - fullEnergy) / std::abs(fullEnergy);
    EXPECT_LT(relative, 0.05)
        << "half depth moved the energy by " << relative
        << " relative; that is a different answer, not a coarser one";
}

/// The parameter reaches Param through the same path a params file takes.
TEST(IrregularDepthScaleTest, TheParameterRoundTripsThroughTheParamsFile)
{
    Param param;
    EXPECT_DOUBLE_EQ(param.irregularPatchDepthScale, 1.0);
    EXPECT_TRUE(import_kv_string("irregularPatchDepthScale", "0.5", param));
    EXPECT_DOUBLE_EQ(param.irregularPatchDepthScale, 0.5);
}

/**
 * @brief Every WP3 and WP4 key is reported as supported.
 *
 * They set their value and then fell out of the else-if chain without
 * returning, so import_kv_string() applied the setting and printed "VARIABLE
 * NOT SUPPORTED" for it. Harmless today only because import_param_file()
 * ignores the result -- but it is the one signal a run has that its parameter
 * file was understood.
 */
TEST(IrregularDepthScaleTest, TheFluidParameterKeysReportThemselvesAsSupported)
{
    Param param;
    EXPECT_TRUE(import_kv_string("inPlaneDynamicsEnabled", "true", param));
    EXPECT_TRUE(param.inPlaneDynamicsEnabled);
    EXPECT_TRUE(import_kv_string("edgeSpringEnabled", "true", param));
    EXPECT_TRUE(param.edgeSpringEnabled);
    EXPECT_TRUE(import_kv_string("edgeSpringConstant", "12.5", param));
    EXPECT_DOUBLE_EQ(param.edgeSpringConstant, 12.5);
    EXPECT_TRUE(import_kv_string("edgeSpringRestLength", "4.5", param));
    EXPECT_DOUBLE_EQ(param.edgeSpringRestLength, 4.5);
    EXPECT_TRUE(import_kv_string("surfaceSolver", "iterative", param));
    EXPECT_EQ(param.surfaceSolver, "iterative");
    EXPECT_TRUE(import_kv_string("edgeFlipEnabled", "true", param));
    EXPECT_TRUE(param.edgeFlipEnabled);
    EXPECT_TRUE(import_kv_string("edgeFlipAttemptRate", "0.25", param));
    EXPECT_DOUBLE_EQ(param.edgeFlipAttemptRate, 0.25);
    EXPECT_TRUE(import_kv_string("edgeFlipInterval", "7", param));
    EXPECT_EQ(param.edgeFlipInterval, 7);
    EXPECT_TRUE(import_kv_string("edgeFlipMinValence", "5", param));
    EXPECT_EQ(param.edgeFlipMinValence, 5);
    EXPECT_TRUE(import_kv_string("edgeFlipMaxValence", "7", param));
    EXPECT_EQ(param.edgeFlipMaxValence, 7);
}

// ---------------------------------------------------------------------------
// Restoring a connectivity
// ---------------------------------------------------------------------------

/**
 * @brief A flipped connectivity restores onto a freshly set up mesh.
 *
 * This is what a restart of a fluid run does, and the check that matters is
 * not that the corner lists match -- they were copied -- but that everything
 * derived from them does, down to the energy. A restore that set the faces and
 * left the one-rings alone would pass a corner comparison and produce a
 * different membrane.
 */
TEST(RestoreConnectivityTest, AFlippedMeshRestoresIntoAFreshOne)
{
    Param flippedParam;
    configure_sphere(flippedParam);
    Mesh flipped(flippedParam);
    build_sphere(flipped, 1);

    const int iEdge = first_flippable_edge(flipped);
    ASSERT_GE(iEdge, 0);
    flipped.flip_edge(iEdge);
    flipped.Compute_Energy_And_Force();

    Param freshParam;
    configure_sphere(freshParam);
    Mesh fresh(freshParam);
    build_sphere(fresh, 1);
    ASSERT_NE(corners_of(fresh), corners_of(flipped)) << "the flip did not change anything";

    std::string why;
    ASSERT_TRUE(fresh.restore_face_connectivity(corners_of(flipped), &why)) << why;
    EXPECT_EQ(corners_of(fresh), corners_of(flipped));
    EXPECT_TRUE(fresh.validate_manifold_topology(&why)) << why;

    // The one-rings, the edge table and the patch tables all came back, which
    // is what the energy is actually a function of.
    fresh.Compute_Energy_And_Force();
    EXPECT_NEAR(fresh.param.energy.energyTotal, flipped.param.energy.energyTotal, 1e-9)
        << "the restored mesh does not carry the same energy as the mesh it was copied from";
    EXPECT_EQ(fresh.edges.size(), flipped.edges.size());
}

/// Restoring the connectivity a mesh already has changes nothing, version
/// included: bumping it would invalidate every cache keyed on it for nothing.
TEST(RestoreConnectivityTest, RestoringTheSameConnectivityIsANoOp)
{
    Param param;
    configure_sphere(param);
    Mesh mesh(param);
    build_sphere(mesh, 1);

    const long long versionBefore = mesh.topologyVersion;
    std::string why;
    EXPECT_TRUE(mesh.restore_face_connectivity(corners_of(mesh), &why)) << why;
    EXPECT_EQ(mesh.topologyVersion, versionBefore);
}

/**
 * @brief A rejected restore leaves the mesh exactly as it was.
 *
 * Half-restoring is worse than refusing: the caller would go on integrating a
 * surface whose one-rings and edge table disagree with its faces, and nothing
 * downstream re-derives them.
 */
TEST(RestoreConnectivityTest, ARejectedRestoreLeavesTheMeshUntouched)
{
    Param param;
    configure_sphere(param);
    Mesh mesh(param);
    build_sphere(mesh, 1);

    const std::vector<std::array<int, 3>> before = corners_of(mesh);
    const long long versionBefore = mesh.topologyVersion;
    mesh.Compute_Energy_And_Force();
    const double energyBefore = mesh.param.energy.energyTotal;

    std::string why;

    // Too few faces.
    std::vector<std::array<int, 3>> truncated = before;
    truncated.pop_back();
    EXPECT_FALSE(mesh.restore_face_connectivity(truncated, &why));
    EXPECT_NE(why.find("faces"), std::string::npos) << why;

    // A vertex that does not exist.
    std::vector<std::array<int, 3>> outOfRange = before;
    outOfRange[0][0] = static_cast<int>(mesh.vertices.size());
    EXPECT_FALSE(mesh.restore_face_connectivity(outOfRange, &why));

    // A repeated corner: a degenerate triangle, not a mesh.
    std::vector<std::array<int, 3>> degenerate = before;
    degenerate[0][1] = degenerate[0][0];
    EXPECT_FALSE(mesh.restore_face_connectivity(degenerate, &why));

    // A well-formed triple list that is not a manifold: give two faces the
    // same corners, so one edge ends up in four triangles.
    std::vector<std::array<int, 3>> nonManifold = before;
    nonManifold[1] = nonManifold[0];
    why.clear();
    EXPECT_FALSE(mesh.restore_face_connectivity(nonManifold, &why));
    EXPECT_FALSE(why.empty()) << "a refusal has to say why";

    EXPECT_EQ(corners_of(mesh), before);
    EXPECT_EQ(mesh.topologyVersion, versionBefore);
    mesh.Compute_Energy_And_Force();
    EXPECT_NEAR(mesh.param.energy.energyTotal, energyBefore, 1e-12)
        << "a refused restore left the mesh in a different state";
}

// ---------------------------------------------------------------------------
// The restart checkpoint's faces block
// ---------------------------------------------------------------------------

/**
 * @brief A checkpoint of a flipped mesh restarts onto the flipped connectivity.
 *
 * Without the faces block the coordinates would land on the triangulation the
 * run was set up with. The mesh would be consistent, the restart would report
 * success, and every energy from there on would be wrong.
 */
TEST(CheckpointFacesTest, ACheckpointCarriesAndRestoresTheConnectivity)
{
    const std::string path = scratch("checkpoint.chk");
    std::remove(path.c_str());

    Param sourceParam;
    configure_sphere(sourceParam);
    Mesh source(sourceParam);
    build_sphere(source, 1);

    const int iEdge = first_flippable_edge(source);
    ASSERT_GE(iEdge, 0);
    source.flip_edge(iEdge);
    source.Compute_Energy_And_Force();
    const long long sourceVersion = source.topologyVersion;

    Record sourceRecord(4);
    sourceRecord.add(source.param.area, source.param.energy, source.calculate_mean_force());
    Model sourceModel(source, sourceRecord);
    sourceModel.stepSize = 0.125;
    ASSERT_TRUE(write_model_restart_checkpoint(sourceModel, path, 7));

    Param targetParam;
    configure_sphere(targetParam);
    Mesh target(targetParam);
    build_sphere(target, 1);
    ASSERT_NE(corners_of(target), corners_of(source));

    Record targetRecord(4);
    Model targetModel(target, targetRecord);
    ASSERT_TRUE(load_model_restart_checkpoint(targetModel, path));

    EXPECT_EQ(targetModel.iteration, 7);
    EXPECT_EQ(corners_of(target), corners_of(source));
    EXPECT_EQ(target.topologyVersion, sourceVersion)
        << "the restarted run did not pick up the connectivity version it left off at";

    std::string why;
    EXPECT_TRUE(target.validate_manifold_topology(&why)) << why;
    target.Compute_Energy_And_Force();
    EXPECT_NEAR(target.param.energy.energyTotal, source.param.energy.energyTotal, 1e-9);

    std::remove(path.c_str());
}

/**
 * @brief A V1 checkpoint still restarts.
 *
 * V1 predates edge flips and carries no connectivity because nothing could
 * change it. Refusing those files would strand every run checkpointed before
 * this branch for no gain.
 */
TEST(CheckpointFacesTest, AVersionOneCheckpointStillLoads)
{
    const std::string v2Path = scratch("v2.chk");
    const std::string v1Path = scratch("v1.chk");
    std::remove(v2Path.c_str());
    std::remove(v1Path.c_str());

    Param sourceParam;
    configure_sphere(sourceParam);
    Mesh source(sourceParam);
    build_sphere(source, 1);
    Record sourceRecord(4);
    sourceRecord.add(source.param.area, source.param.energy, source.calculate_mean_force());
    Model sourceModel(source, sourceRecord);
    sourceModel.stepSize = 0.125;
    ASSERT_TRUE(write_model_restart_checkpoint(sourceModel, v2Path, 3));

    // Strip the block back out, which is exactly the file V1 would have
    // written: the tag, and no faces.
    {
        std::ifstream in(v2Path);
        std::ofstream out(v1Path);
        std::string line;
        bool skipping = false;
        while (std::getline(in, line))
        {
            if (line == "SLIMED_RESTART_V2")
            {
                out << "SLIMED_RESTART_V1\n";
                continue;
            }
            if (line.rfind("faces ", 0) == 0)
            {
                skipping = true;
                continue;
            }
            if (skipping)
            {
                if (line.rfind("scaffoldingPoints ", 0) == 0)
                {
                    skipping = false;
                }
                else
                {
                    continue;
                }
            }
            out << line << '\n';
        }
    }

    Param targetParam;
    configure_sphere(targetParam);
    Mesh target(targetParam);
    build_sphere(target, 1);
    Record targetRecord(4);
    Model targetModel(target, targetRecord);
    EXPECT_TRUE(load_model_restart_checkpoint(targetModel, v1Path));
    EXPECT_EQ(targetModel.iteration, 3);

    std::remove(v2Path.c_str());
    std::remove(v1Path.c_str());
}

/// A checkpoint from a different mesh is refused rather than half-applied.
TEST(CheckpointFacesTest, ACheckpointFromADifferentMeshIsRefused)
{
    const std::string path = scratch("mismatch.chk");
    std::remove(path.c_str());

    Param sourceParam;
    configure_sphere(sourceParam);
    Mesh source(sourceParam);
    build_sphere(source, 1);
    Record sourceRecord(4);
    Model sourceModel(source, sourceRecord);
    ASSERT_TRUE(write_model_restart_checkpoint(sourceModel, path, 1));

    // A coarser sphere: fewer vertices, so the vertex block already refuses.
    // Corrupt the face count instead, on a mesh whose vertices do match, so
    // the faces block is what has to catch it.
    std::string contents;
    {
        std::ifstream in(path);
        std::ostringstream buffer;
        buffer << in.rdbuf();
        contents = buffer.str();
    }
    const std::string marker = "faces " + std::to_string(source.faces.size()) + ' ';
    const std::size_t at = contents.find(marker);
    ASSERT_NE(at, std::string::npos);
    contents.replace(at, marker.size(), "faces " + std::to_string(source.faces.size() + 1) + ' ');
    {
        std::ofstream out(path);
        out << contents;
    }

    Param targetParam;
    configure_sphere(targetParam);
    Mesh target(targetParam);
    build_sphere(target, 1);
    Record targetRecord(4);
    Model targetModel(target, targetRecord);
    EXPECT_FALSE(load_model_restart_checkpoint(targetModel, path));

    std::remove(path.c_str());
}

// ---------------------------------------------------------------------------
// Per-frame connectivity output
// ---------------------------------------------------------------------------

/**
 * @brief A face frame is written once, and again only when the mesh moved.
 *
 * Every frame would be correct and unaffordable -- one face list per frame on
 * a mesh with thousands of triangles, for a run that outputs millions of
 * frames. Never writing again is what the tree did, and it is wrong the
 * moment a flip lands. The rule is: write when it changed.
 */
TEST(FaceFrameTest, AFrameIsWrittenOnceAndAgainOnlyWhenTheConnectivityMoves)
{
    Param param;
    configure_flat(param);
    param.meshpointOutput = true;
    DynamicMesh mesh(param);
    mesh.setup_flat();
    mesh.update_vertices_mat_with_vector();
    mesh.apply_mesh_to_surface();

    const std::string prefix = scratch("frames");
    dynamics_create_trajectory_files(mesh, prefix);

    EXPECT_TRUE(dynamics_output_face_frame(mesh, prefix, 0))
        << "the first frame has no baseline to be compared against and must be written";
    EXPECT_FALSE(dynamics_output_face_frame(mesh, prefix, 1))
        << "nothing changed, so nothing should have been written";
    EXPECT_FALSE(dynamics_output_face_frame(mesh, prefix, 2));

    const int iEdge = first_flippable_edge(mesh);
    ASSERT_GE(iEdge, 0);
    mesh.flip_edge(iEdge);

    EXPECT_TRUE(dynamics_output_face_frame(mesh, prefix, 3));
    EXPECT_FALSE(dynamics_output_face_frame(mesh, prefix, 4));

    // What landed on disk is the connectivity as of that frame, in the layout
    // face.csv uses, so an analysis reads it with what it already has.
    const std::vector<std::array<int, 3>> atZero = read_face_csv(prefix + "face_0.csv");
    const std::vector<std::array<int, 3>> atThree = read_face_csv(prefix + "face_3.csv");
    EXPECT_EQ(atThree, corners_of(mesh));
    EXPECT_EQ(atZero.size(), mesh.faces.size());
    EXPECT_NE(atZero, atThree) << "the two frames recorded the same connectivity";

    for (const std::string &path : {prefix + "face_0.csv", prefix + "face_3.csv",
                                    "meshpoint" + prefix + ".csv",
                                    "surfacepoint" + prefix + ".csv"})
    {
        std::remove(path.c_str());
    }
}

// ---------------------------------------------------------------------------
// The run gate
// ---------------------------------------------------------------------------

/**
 * @brief A refused or rejected trial leaves the connectivity version alone.
 *
 * evaluate_edge_flip() flips the edge, measures, and flips it back. Each of
 * those bumps topologyVersion, so the trial used to leave the mesh exactly as
 * it was and the version two ahead. Nothing became wrong -- over-invalidation
 * is safe -- but the version is what says whether a cache is stale and whether
 * a frame needs its connectivity written beside it, and on a fluid run most
 * attempts are rejected. A 500-step run with flips on and nothing accepted
 * wrote 112 identical face frames and rebuilt the sparse limit mask 121 times.
 */
TEST(FluidRunTest, ARejectedTrialDoesNotMoveTheConnectivityVersion)
{
    Param param;
    configure_sphere(param);
    param.edgeFlipEnabled = true;
    Mesh mesh(param);
    build_sphere(mesh, 1);

    const int iEdge = first_flippable_edge(mesh);
    ASSERT_GE(iEdge, 0);

    const long long versionBefore = mesh.topologyVersion;
    const std::vector<std::array<int, 3>> trianglesBefore = triangle_set_of(mesh);

    EdgeFlipDelta delta;
    ASSERT_TRUE(mesh.evaluate_edge_flip(iEdge, delta));

    EXPECT_EQ(triangle_set_of(mesh), trianglesBefore) << "the trial did not put the mesh back";
    EXPECT_EQ(mesh.topologyVersion, versionBefore)
        << "the trial put the mesh back but left the version " << mesh.topologyVersion
        << " instead of " << versionBefore;

    // And an accepted flip still moves it, by exactly one.
    mesh.flip_edge(iEdge);
    EXPECT_EQ(mesh.topologyVersion, versionBefore + 1);
    EXPECT_NE(triangle_set_of(mesh), trianglesBefore);
}

/**
 * @brief The connectivity-keyed caches rebuild on a version change and not
 * otherwise.
 *
 * The device layout is entirely connectivity-derived, down to the CSR of every
 * one-ring, and an edge flip rewrites the connectivity while leaving the face
 * and vertex counts exactly as they were. A cache that checks counts would go
 * on being used against a mesh it no longer describes -- an out-of-bounds read
 * rather than a failure. topologyVersion is what answers the question, and
 * this is the gate on it.
 *
 * The version is moved by hand rather than by a flip, because a flipped mesh
 * is one the layout cannot build at all: see the GPU test below.
 */
TEST(FluidRunTest, TheDeviceLayoutRebuildsOnAVersionChangeAndNotOtherwise)
{
    Param param;
    configure_flat(param);
    DynamicMesh mesh(param);
    mesh.setup_flat();

    mesh.ensure_device_layout();
    EXPECT_EQ(mesh.deviceLayoutTopologyVersion, mesh.topologyVersion);
    const int nFaces = mesh.deviceLayout.nFaces();
    ASSERT_GT(nFaces, 0);

    // Nothing moved: the same layout is kept.
    const long long versionAfterBuild = mesh.deviceLayoutTopologyVersion;
    mesh.ensure_device_layout();
    EXPECT_EQ(mesh.deviceLayoutTopologyVersion, versionAfterBuild);

    // The connectivity moved. The counts did not, which is the case a
    // count-based check misses entirely.
    mesh.topologyVersion += 1;
    mesh.ensure_device_layout();
    EXPECT_EQ(mesh.deviceLayoutTopologyVersion, mesh.topologyVersion)
        << "the layout did not notice a connectivity change that left the counts alone";
    EXPECT_EQ(mesh.deviceLayout.nFaces(), nFaces);
}

/**
 * @brief Flips and the GPU backend are refused together, before the first step.
 *
 * A flip leaves extraordinary corners at both ends of the new edge, so a fluid
 * mesh is full of faces with more than one of them, and DeviceMeshLayout
 * cannot build those. Without this the run starts, flips, and throws out of
 * the layout builder with a mesh already changed -- a crash partway through
 * rather than a configuration error.
 */
TEST(FluidRunTest, FlipsWithTheGpuBackendAreRefusedAtSetup)
{
    Param param;
    configure_flat(param);
    param.surfaceSolver = "iterative";
    param.edgeFlipEnabled = true;
    param.forceBackend = "gpu";

    DynamicMesh mesh(param);
    EXPECT_THROW(mesh.setup_flat(), std::runtime_error);

    // And the supported combination is not refused.
    Param cpuParam;
    configure_flat(cpuParam);
    cpuParam.surfaceSolver = "iterative";
    cpuParam.edgeFlipEnabled = true;
    cpuParam.forceBackend = "cpu";
    DynamicMesh cpuMesh(cpuParam);
    EXPECT_NO_THROW(cpuMesh.setup_flat());
}

/**
 * @brief A short fluid run writes frames that pair with a consistent mesh.
 *
 * This is the loop of run_dynamics_flat() with flips on, run here so the gate
 * can look inside it. The failure it exists to catch is silent: a coordinate
 * frame whose connectivity was never written is read back against the face
 * list from setup, which parses, has the right triangle count, and describes a
 * different membrane.
 */
TEST(FluidRunTest, AShortRunPairsEveryFrameWithAConsistentConnectivity)
{
    Param param;
    configure_flat(param);
    param.meshpointOutput = true;
    param.surfaceSolver = "iterative";
    param.inPlaneDynamicsEnabled = true;
    param.edgeSpringEnabled = true;
    param.edgeFlipEnabled = true;
    param.edgeFlipInterval = 1;
    param.maxIterations = 12;

    DynamicMesh mesh(param);
    mesh.setup_flat();
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.Compute_Energy_And_Force();
    mesh.update_vertices_mat_with_vector();
    mesh.apply_mesh_to_surface();

    // A rate that gives a handful of attempts per sweep on this mesh, rather
    // than one scaled from a physical nu that would give thousands.
    int flippable = 0;
    for (const MeshEdge &edge : mesh.edges)
    {
        flippable += edge.flippable ? 1 : 0;
    }
    ASSERT_GT(flippable, 0);
    param.edgeFlipAttemptRate = 6.0 / (param.timeStep * flippable);
    mesh.param.edgeFlipAttemptRate = param.edgeFlipAttemptRate;

    Record record(param.maxIterations);
    DynamicModel model(mesh, record);

    const std::string prefix = scratch("run");
    dynamics_create_trajectory_files(mesh, prefix);
    dynamics_output_face_frame(mesh, prefix, 0);

    std::vector<std::string> writtenFaceFrames = {prefix + "face_0.csv"};
    int totalAccepted = 0;
    int stepsWithAFlip = 0;
    int versionMoves = 0;
    long long previousVersion = mesh.topologyVersion;

    for (model.iteration = 0; model.iteration < param.maxIterations; model.iteration++)
    {
        mesh.apply_mesh_to_surface();
        model.next_step();
        mesh.apply_surface_to_mesh();
        mesh.postprocess_ghost_periodic();
        mesh.update_vertices_vector_with_mat();

        const std::string path = prefix + "face_" + std::to_string(model.iteration + 1) + ".csv";
        dynamics_output_trajectory_files(mesh, prefix);
        if (dynamics_output_face_frame(mesh, prefix, model.iteration + 1))
        {
            writtenFaceFrames.push_back(path);
        }

        record.add(mesh.param.area, mesh.param.energy, mesh.calculate_mean_force());
        mesh.Compute_Energy_And_Force();

        const EdgeFlipSweepStats stats = mesh.edge_flip_sweep(model.iteration);
        totalAccepted += stats.accepted;
        if (stats.accepted > 0)
        {
            stepsWithAFlip++;
            mesh.Compute_Energy_And_Force();
        }

        // The signal every connectivity-keyed cache invalidates on. It has to
        // move on exactly the steps a flip landed: not moving leaves a stale
        // cache, and moving without a flip throws away a correct one.
        if (mesh.topologyVersion != previousVersion)
        {
            versionMoves++;
            previousVersion = mesh.topologyVersion;
        }
    }

    ASSERT_GT(totalAccepted, 0) << "no flip was accepted, so the gate tested nothing";
    EXPECT_EQ(versionMoves, stepsWithAFlip)
        << "the connectivity version moved on " << versionMoves << " steps but flips landed on "
        << stepsWithAFlip;

    std::string why;
    EXPECT_TRUE(mesh.validate_manifold_topology(&why)) << why;

    // Every face frame written is a consistent mesh, and the last one is the
    // mesh the run finished with.
    EXPECT_GT(writtenFaceFrames.size(), 1u)
        << "flips landed but no connectivity was ever written beside a frame";
    for (const std::string &path : writtenFaceFrames)
    {
        const std::vector<std::array<int, 3>> corners = read_face_csv(path);
        EXPECT_EQ(corners.size(), mesh.faces.size()) << path;
        for (const std::array<int, 3> &triple : corners)
        {
            for (int k = 0; k < 3; k++)
            {
                EXPECT_GE(triple[k], 0) << path;
                EXPECT_LT(triple[k], static_cast<int>(mesh.vertices.size())) << path;
            }
        }
        // The strong form: a frame reloads into a consistent mesh.
        Param replayParam;
        configure_flat(replayParam);
        DynamicMesh replay(replayParam);
        replay.setup_flat();
        std::string replayWhy;
        EXPECT_TRUE(replay.restore_face_connectivity(corners, &replayWhy)) << path << ": "
                                                                          << replayWhy;
    }

    for (const std::string &path : writtenFaceFrames)
    {
        std::remove(path.c_str());
    }
    std::remove(("meshpoint" + prefix + ".csv").c_str());
    std::remove(("surfacepoint" + prefix + ".csv").c_str());
}
