#include "test_mixed_boundary.hpp"

/**
 * The global boundary modes decide everything from grid position: which rows
 * are ghosts, which vertex duplicates which. That describes a rectangular
 * sheet with the same boundary on every side and nothing else. Under
 * BoundaryType::Mixed each vertex says what it is -- free, fixed, or a
 * periodic image of another vertex -- and the treatment follows from that,
 * so one mesh can wrap along one axis and be clamped along the other, and a
 * mesh read from a file needs no grid at all.
 *
 * Three things have to be right for that to be a boundary condition rather
 * than a labelling, and each has a test below.
 *
 * **An image is not a coordinate.** It sits at its source plus a fixed
 * offset, always: after a minimization step, after a dynamics step, after a
 * surface-to-control solve. If it ever lags, the seam is a crack.
 *
 * **The force on a source is the gradient of the periodic energy.** The
 * energy is a function of the sources alone, so by the chain rule the force
 * on a source is its own accumulated force plus that of every image. The
 * global Periodic mode drops the image forces; this mode folds them, and the
 * finite-difference check is what says the folding is exact -- including at
 * the seam, where the regularization of the copy faces must not be counted
 * twice.
 *
 * **The mesh files round-trip.** A mesh SLIMED wrote must load back as the
 * same mesh, types, mirrors, offsets and copy flags included, or a run from
 * a file is a different run from the one that produced it.
 *
 * @see docs/mixed_boundary_conditions.md
 */

namespace
{

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

void configure_common(Param &param, double side)
{
    param.VERBOSE_MODE = false;
    param.sideX = side;
    param.sideY = side;
    param.lFace = 5.0;
    param.kCurv = 83.4;
    param.KBT = 4.17;
    param.uSurf = 250.0;
    param.uVol = 0.0;
    param.isGlobalConstraint = true;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
    param.randomSeed = 20260912u;
    param.timeStep = 1.0e-3;
    param.diffConst = 1.0;
    param.surfaceSolver = "iterative";
}

/// A flat sheet under the per-vertex mode, with the boundary of each axis.
void configure_sheet(Param &param, BoundaryType x, BoundaryType y, double side = 60.0)
{
    configure_common(param, side);
    param.boundaryCondition = BoundaryType::Mixed;
    param.boundaryConditionX = x;
    param.boundaryConditionY = y;
}

/// The same sheet under the global Periodic mode, for comparison.
void configure_legacy_periodic(Param &param, double side = 60.0)
{
    configure_common(param, side);
    param.boundaryCondition = BoundaryType::Periodic;
}

int grid_index(const Param &param, int i, int j)
{
    return (param.nFaceX + 1) * j + i;
}

/**
 * @brief Lift the sheet into a smooth relief and re-sync the images.
 *
 * z = a sin(2 pi x / lx) cos(2 pi y / ly) + 10, with lx and ly the tile
 * lengths along a periodic axis and the full sheet otherwise, so the relief
 * is continuous across a seam. The images are then put where the mode says
 * they belong -- exactly at their source plus offset -- which the formula
 * only does to rounding.
 */
void undulate(Mesh &mesh, double amplitude)
{
    const Param &param = mesh.param;
    const bool mixed = (param.boundaryCondition == BoundaryType::Mixed);
    const bool periodicX = mixed ? (param.boundaryConditionX == BoundaryType::Periodic)
                                 : (param.boundaryCondition == BoundaryType::Periodic);
    const bool periodicY = mixed ? (param.boundaryConditionY == BoundaryType::Periodic)
                                 : (param.boundaryCondition == BoundaryType::Periodic);
    const double lx = (periodicX ? param.nFaceX - 6 : param.nFaceX) * param.dFaceX;
    const double ly = (periodicY ? param.nFaceY - 6 : param.nFaceY) * param.dFaceY;
    for (Vertex &vertex : mesh.vertices)
    {
        const double x = vertex.coord.get(0, 0);
        const double y = vertex.coord.get(1, 0);
        vertex.coord.set(2, 0,
                         amplitude * std::sin(2.0 * M_PI * x / lx) * std::cos(2.0 * M_PI * y / ly) +
                             10.0);
    }
    mesh.sync_periodic_images();
}

/**
 * @brief Everything the run does before its first force evaluation, with the
 * flat sheet as the reference configuration and a relief on top of it.
 *
 * The reference is taken flat so that the regularization is live in every
 * force below: a term that was zero could not reveal a double count.
 */
void prepare(Mesh &mesh, double amplitude)
{
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.param.vol0 = 0.0;
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    undulate(mesh, amplitude);
    mesh.update_previous_coord_for_vertex();
    mesh.Compute_Energy_And_Force();
    mesh.update_previous_force_for_vertex();
}

double total_energy(Mesh &mesh)
{
    mesh.Compute_Energy_And_Force();
    return mesh.param.energy.energyTotal;
}

double bending_energy(const Mesh &mesh)
{
    double sum = 0.0;
    for (const Face &face : mesh.faces)
    {
        sum += face.energy.energyCurvature;
    }
    return sum;
}

int count_independent(const Mesh &mesh)
{
    int n = 0;
    for (int v = 0; v < static_cast<int>(mesh.vertices.size()); v++)
    {
        n += mesh.is_independent_vertex(v) ? 1 : 0;
    }
    return n;
}

double force_magnitude(const Vertex &vertex)
{
    double sum = 0.0;
    for (int axis = 0; axis < 3; axis++)
    {
        const double f = vertex.force.forceTotal.get(axis, 0);
        sum += f * f;
    }
    return std::sqrt(sum);
}

/// Whether an image sits exactly at its source plus its offset.
::testing::AssertionResult images_are_in_place(const Mesh &mesh)
{
    for (int v : mesh.periodicImageVertices)
    {
        const Vertex &image = mesh.vertices[v];
        const Vertex &source = mesh.vertices[image.reflectiveVertexIndex];
        for (int axis = 0; axis < 3; axis++)
        {
            const double expected = source.coord.get(axis, 0) + image.mirrorOffset[axis];
            if (image.coord.get(axis, 0) != expected)
            {
                return ::testing::AssertionFailure()
                       << "image " << v << " axis " << axis << " is at "
                       << image.coord.get(axis, 0) << ", its source " << source.index
                       << " plus offset is at " << expected;
            }
        }
    }
    return ::testing::AssertionSuccess();
}

/// A small closed solid, for the refusal tests: 6 vertices at valence 4.
void octahedron(std::vector<std::vector<double>> &vertices, std::vector<std::vector<int>> &faces)
{
    vertices = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    faces = {{0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4}, {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}};
}

std::string scratch(const std::string &name)
{
    return "mixed_boundary_" + name;
}

} // namespace

// ---------------------------------------------------------------------------
// Parameters and spellings
// ---------------------------------------------------------------------------

TEST(MixedBoundaryTest, TheParameterKeysRoundTrip)
{
    Param param;
    EXPECT_TRUE(import_kv_string("boundaryType", "Mixed", param));
    EXPECT_EQ(param.boundaryCondition, BoundaryType::Mixed);
    EXPECT_TRUE(import_kv_string("boundaryTypeX", "free", param));
    EXPECT_EQ(param.boundaryConditionX, BoundaryType::Free);
    EXPECT_TRUE(import_kv_string("boundaryTypeY", "Fixed", param));
    EXPECT_EQ(param.boundaryConditionY, BoundaryType::Fixed);
    EXPECT_TRUE(import_kv_string("boundaryTypeY", "Periodic", param));
    EXPECT_EQ(param.boundaryConditionY, BoundaryType::Periodic);
    EXPECT_TRUE(import_kv_string("fixedBoundaryRings", "2", param));
    EXPECT_EQ(param.fixedBoundaryRings, 2);
    EXPECT_TRUE(import_kv_string("meshVerticesFile", "data/example/vertices_flat.csv", param));
    EXPECT_EQ(param.meshVerticesFile, "data/example/vertices_flat.csv");
    EXPECT_TRUE(import_kv_string("meshFacesFile", "data/example/faces_flat.csv", param));
    EXPECT_EQ(param.meshFacesFile, "data/example/faces_flat.csv");

    // A per-axis type is one of the three; a typo stops the run.
    EXPECT_THROW(import_kv_string("boundaryTypeX", "Mixed", param), std::runtime_error);
    EXPECT_THROW(import_kv_string("boundaryTypeY", "sideways", param), std::runtime_error);

    VertexType type = VertexType::Free;
    EXPECT_TRUE(parse_vertex_type("Periodic", type));
    EXPECT_EQ(type, VertexType::Periodic);
    EXPECT_TRUE(parse_vertex_type(" fixed ", type));
    EXPECT_EQ(type, VertexType::Fixed);
    EXPECT_TRUE(parse_vertex_type("3", type)); // the legacy integer spelling
    EXPECT_EQ(type, VertexType::Periodic);
    EXPECT_TRUE(parse_vertex_type("real", type));
    EXPECT_EQ(type, VertexType::Free);
    EXPECT_FALSE(parse_vertex_type("solid", type));
    EXPECT_STREQ(vertex_type_name(VertexType::Fixed), "fixed");
    EXPECT_STREQ(vertex_type_name(VertexType::Periodic), "periodic");

    BoundaryType boundary = BoundaryType::Fixed;
    EXPECT_TRUE(parse_boundary_type("mixed", boundary));
    EXPECT_EQ(boundary, BoundaryType::Mixed);
    EXPECT_FALSE(parse_boundary_type("open", boundary));
    EXPECT_STREQ(boundary_type_name(BoundaryType::Mixed), "Mixed");

    // The enum values the output format relies on have not moved.
    EXPECT_EQ(static_cast<int>(VertexType::Free), 0);
    EXPECT_EQ(static_cast<int>(VertexType::Fixed), 1);
    EXPECT_EQ(static_cast<int>(VertexType::Periodic), 3);
    EXPECT_EQ(static_cast<int>(VertexType::Ghost), 4);
    EXPECT_EQ(VertexType::Real, VertexType::Free);
    EXPECT_EQ(VertexType::PeriodicReflectiveBoundary, VertexType::Periodic);
}

// ---------------------------------------------------------------------------
// The generated sheet
// ---------------------------------------------------------------------------

/**
 * @brief Periodic on both axes, the per-vertex sheet is the global one's
 * layout expressed as images.
 *
 * Same copy faces, an image wherever the global mode has a ghost, sources
 * exactly one period apart from their images, and the same bending energy
 * for the same relief. What differs is what happens to the image forces, and
 * that is the next test.
 */
TEST(MixedBoundaryTest, ThePeriodicSheetLaysOutTheSameBandAsTheGlobalMode)
{
    Param mixedParam;
    configure_sheet(mixedParam, BoundaryType::Periodic, BoundaryType::Periodic);
    Mesh mixed(mixedParam);
    ASSERT_NO_THROW(mixed.setup_flat());

    Param legacyParam;
    configure_legacy_periodic(legacyParam);
    Mesh legacy(legacyParam);
    ASSERT_NO_THROW(legacy.setup_flat());

    ASSERT_EQ(mixed.vertices.size(), legacy.vertices.size());
    ASSERT_EQ(mixed.faces.size(), legacy.faces.size());
    for (std::size_t f = 0; f < mixed.faces.size(); f++)
    {
        EXPECT_EQ(mixed.faces[f].isGhost, legacy.faces[f].isGhost) << "face " << f;
    }

    const int periodX = mixedParam.nFaceX - 6;
    const int periodY = mixedParam.nFaceY - 6;
    ASSERT_GE(periodX, 4);
    ASSERT_GE(periodY, 4);
    const double lx = periodX * mixedParam.dFaceX;
    const double ly = periodY * mixedParam.dFaceY;

    int nSources = 0;
    int nImages = 0;
    int nImagesTouchingRealFaces = 0;
    for (std::size_t v = 0; v < mixed.vertices.size(); v++)
    {
        const Vertex &vertex = mixed.vertices[v];
        if (legacy.vertices[v].isGhost)
        {
            EXPECT_TRUE(vertex.is_periodic_image()) << "vertex " << v;
            EXPECT_TRUE(vertex.isGhost) << "vertex " << v;
        }
        if (vertex.is_periodic_image())
        {
            nImages++;
            nImagesTouchingRealFaces += vertex.isGhost ? 0 : 1;
            const Vertex &source = mixed.vertices[vertex.reflectiveVertexIndex];
            EXPECT_EQ(source.type, VertexType::Free) << "vertex " << v;
            EXPECT_FALSE(source.is_periodic_image()) << "vertex " << v;
            // One lattice period along each wrapped axis, nothing along z.
            const double dx = std::abs(vertex.mirrorOffset[0]);
            const double dy = std::abs(vertex.mirrorOffset[1]);
            EXPECT_TRUE(std::abs(dx) < 1e-9 || std::abs(dx - lx) < 1e-9) << "vertex " << v << " dx " << dx;
            EXPECT_TRUE(std::abs(dy) < 1e-9 || std::abs(dy - ly) < 1e-9) << "vertex " << v << " dy " << dy;
            EXPECT_TRUE(dx > 1e-9 || dy > 1e-9) << "vertex " << v << " mirrors a vertex at the same place";
            EXPECT_EQ(vertex.mirrorOffset[2], 0.0);
        }
        else
        {
            nSources++;
            EXPECT_EQ(vertex.type, VertexType::Free) << "vertex " << v;
            EXPECT_FALSE(vertex.isGhost) << "vertex " << v;
            EXPECT_TRUE(mixed.is_independent_vertex(static_cast<int>(v)));
        }
    }
    EXPECT_EQ(nSources, periodX * periodY);
    EXPECT_EQ(nSources + nImages, static_cast<int>(mixed.vertices.size()));
    EXPECT_EQ(static_cast<int>(mixed.periodicImageVertices.size()), nImages);
    // The duplicate ring: images that corner real faces. On a tile of
    // periodX by periodY it is one row and one column, sharing a corner.
    EXPECT_EQ(nImagesTouchingRealFaces, periodX + periodY + 1);

    // Images are frozen for the flip move, as the global mode's duplicates are
    // once DynamicMesh has marked them; here it holds from setup.
    ASSERT_EQ(mixed.flipFrozenVertex.size(), mixed.vertices.size());
    for (int v : mixed.periodicImageVertices)
    {
        EXPECT_NE(mixed.flipFrozenVertex[v], 0);
    }

    prepare(mixed, 1.0);
    prepare(legacy, 1.0);
    EXPECT_TRUE(images_are_in_place(mixed));
    const double mixedBending = bending_energy(mixed);
    const double legacyBending = bending_energy(legacy);
    ASSERT_GT(std::abs(legacyBending), 1.0) << "the relief has to carry some bending energy";
    EXPECT_NEAR(mixedBending, legacyBending, 1e-9 * std::abs(legacyBending));
}

/**
 * @brief The folded force on a source is minus the gradient of the periodic
 * energy -- at the seam as well as in the interior.
 *
 * The move a periodic membrane can make is "this source, and every image of
 * it, together", so that is what the finite difference moves. The
 * regularization is live (flat reference, undulated sheet) and the copy
 * faces are the ones that would double count it, so the seam probes are the
 * ones that matter.
 */
TEST(MixedBoundaryTest, TheFoldedForceIsMinusTheGradientOfThePeriodicEnergy)
{
    Param param;
    configure_sheet(param, BoundaryType::Periodic, BoundaryType::Periodic);
    Mesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_flat());
    prepare(mesh, 1.5);

    const int nX = param.nFaceX;
    const int nY = param.nFaceY;
    // Sources on the seam -- row 3, column 3, the far row and column -- whose
    // physical fans are split between themselves and their images, and one
    // deep inside for reference.
    const std::vector<int> probes = {grid_index(param, 3, 3),      grid_index(param, 7, 3),
                                     grid_index(param, 3, 6),      grid_index(param, nX - 4, 5),
                                     grid_index(param, 6, nY - 4), grid_index(param, 6, 6)};

    const double h = 1e-5;
    double largestForce = 0.0;
    for (int v : probes)
    {
        ASSERT_TRUE(mesh.is_independent_vertex(v)) << "vertex " << v;
        for (int axis = 0; axis < 3; axis++)
        {
            total_energy(mesh);
            const double analytic = mesh.vertices[v].force.forceTotal.get(axis, 0);
            largestForce = std::max(largestForce, std::abs(analytic));

            const double original = mesh.vertices[v].coord.get(axis, 0);
            mesh.vertices[v].coord.set(axis, 0, original + h);
            mesh.sync_periodic_images();
            const double plus = total_energy(mesh);
            mesh.vertices[v].coord.set(axis, 0, original - h);
            mesh.sync_periodic_images();
            const double minus = total_energy(mesh);
            mesh.vertices[v].coord.set(axis, 0, original);
            mesh.sync_periodic_images();

            const double numeric = -(plus - minus) / (2.0 * h);
            EXPECT_NEAR(analytic, numeric, 1e-5 * std::max(1.0, std::abs(numeric)))
                << "vertex " << v << " axis " << axis;
        }
    }
    EXPECT_GT(largestForce, 1e-3) << "the relief has to put some force on the probes";
}

/**
 * @brief Periodic along x and clamped along y: images on the sides, fixed
 * vertices on the top and bottom, nothing wrapped along y.
 */
TEST(MixedBoundaryTest, APeriodicXFixedYSheetWrapsTheSidesAndClampsTheEdges)
{
    Param param;
    configure_sheet(param, BoundaryType::Periodic, BoundaryType::Fixed);
    param.fixedBoundaryRings = 2;
    Mesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_flat());

    const int nX = param.nFaceX;
    const int nY = param.nFaceY;
    const int periodX = nX - 6;
    int nFixed = 0;
    for (int j = 0; j <= nY; j++)
    {
        for (int i = 0; i <= nX; i++)
        {
            const Vertex &vertex = mesh.vertices[grid_index(param, i, j)];
            const bool imageColumn = (i < 3 || i > nX - 4);
            if (imageColumn)
            {
                ASSERT_TRUE(vertex.is_periodic_image()) << "(" << i << ", " << j << ")";
                const int sourceI = 3 + (((i - 3) % periodX) + periodX) % periodX;
                EXPECT_EQ(vertex.reflectiveVertexIndex, grid_index(param, sourceI, j));
                EXPECT_EQ(vertex.mirrorOffset[1], 0.0) << "nothing wraps along y";
                continue;
            }
            const bool clamped = (j < 2 || j > nY - 2);
            EXPECT_EQ(vertex.type, clamped ? VertexType::Fixed : VertexType::Free)
                << "(" << i << ", " << j << ")";
            EXPECT_EQ(mesh.is_independent_vertex(vertex.index), !clamped);
            nFixed += clamped ? 1 : 0;
        }
    }
    EXPECT_EQ(nFixed, 4 * periodX);
    EXPECT_EQ(count_independent(mesh), periodX * (nY + 1 - 4));

    // Copies are the side bands only; the top and bottom rows are real and,
    // being on an open edge, carry no patch.
    for (int j = 0; j < nY; j++)
    {
        for (int i = 0; i < nX; i++)
        {
            const int f = 2 * nX * j + 2 * i;
            const bool copy = (i < 3 || i > nX - 4);
            EXPECT_EQ(mesh.faces[f].isGhost, copy) << "cell (" << i << ", " << j << ")";
            EXPECT_EQ(mesh.faces[f + 1].isGhost, copy) << "cell (" << i << ", " << j << ")";
            if (copy)
            {
                continue;
            }
            const PatchKind expected =
                (j == 0 || j == nY - 1) ? PatchKind::Boundary : PatchKind::Regular;
            EXPECT_EQ(mesh.faces[f].patchKind, expected) << "cell (" << i << ", " << j << ")";
        }
    }

    // Clamped vertices carry no force; sources next to them do.
    prepare(mesh, 1.0);
    double largestFreeForce = 0.0;
    for (const Vertex &vertex : mesh.vertices)
    {
        if (vertex.is_fixed())
        {
            EXPECT_EQ(force_magnitude(vertex), 0.0) << "vertex " << vertex.index;
        }
        else if (vertex.type == VertexType::Free)
        {
            largestFreeForce = std::max(largestFreeForce, force_magnitude(vertex));
        }
    }
    EXPECT_GT(largestFreeForce, 1e-3);
}

/**
 * @brief An open edge: its vertices are degrees of freedom, its outermost
 * faces carry no patch, and the edge feels the faces one ring in.
 */
TEST(MixedBoundaryTest, AFreeEdgeLeavesItsVerticesAsDegreesOfFreedom)
{
    Param param;
    configure_sheet(param, BoundaryType::Periodic, BoundaryType::Free);
    DynamicMesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_flat());

    const int nX = param.nFaceX;
    const int nY = param.nFaceY;
    const int edgeVertex = grid_index(param, 6, 0);
    EXPECT_EQ(mesh.vertices[edgeVertex].type, VertexType::Free);
    EXPECT_FALSE(mesh.vertices[edgeVertex].isGhost);
    EXPECT_TRUE(mesh.is_independent_vertex(edgeVertex));
    EXPECT_FALSE(mesh.is_interior_vertex(edgeVertex)) << "an edge vertex has an open fan";
    EXPECT_EQ(count_independent(mesh), (nX - 6) * (nY + 1));

    for (int i = 3; i < nX - 3; i++)
    {
        EXPECT_EQ(mesh.faces[2 * nX * 0 + 2 * i].patchKind, PatchKind::Boundary);
        EXPECT_EQ(mesh.faces[2 * nX * 1 + 2 * i].patchKind, PatchKind::Regular);
    }

    prepare(mesh, 1.0);
    EXPECT_GT(force_magnitude(mesh.vertices[edgeVertex]), 1e-6)
        << "the edge vertex is a control point of the patches one ring in";

    // The solver integrates exactly the independent vertices, edge included.
    mesh.update_vertices_mat_with_vector();
    mesh.ensure_surface_solver();
    EXPECT_EQ(mesh.surfaceSolver.nFree(), count_independent(mesh));

    // A sheet with images encloses nothing, so a volume constraint on it is
    // refused at setup like it is for the global modes.
    Param constrained;
    configure_sheet(constrained, BoundaryType::Periodic, BoundaryType::Free);
    constrained.uVol = 1.0;
    Mesh open(constrained);
    EXPECT_THROW(open.setup_flat(), std::runtime_error);
}

// ---------------------------------------------------------------------------
// The images stay in place through every way the mesh moves
// ---------------------------------------------------------------------------

TEST(MixedBoundaryTest, AMinimizationStepHoldsTheClampedVerticesAndCarriesTheImages)
{
    Param param;
    configure_sheet(param, BoundaryType::Periodic, BoundaryType::Fixed);
    Mesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_flat());
    prepare(mesh, 1.5);

    Record record(8);
    record.add(mesh.param.area, mesh.param.energy, mesh.calculate_mean_force());
    Model model(mesh, record);

    std::vector<std::array<double, 3>> before(mesh.vertices.size());
    for (std::size_t v = 0; v < mesh.vertices.size(); v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            before[v][axis] = mesh.vertices[v].coord.get(axis, 0);
        }
    }
    const double energyBefore = mesh.param.energy.energyTotal;

    model.determine_trial_step_size();
    ASSERT_GT(model.oa.trialStepSize, 0.0);
    const double step = model.linear_search_for_stepsize_to_minimize_energy();
    ASSERT_GT(step, 0.0);
    model.update_vertex_using_NCG();
    mesh.Compute_Energy_And_Force();

    EXPECT_LT(mesh.param.energy.energyTotal, energyBefore);
    EXPECT_TRUE(images_are_in_place(mesh));
    bool somethingMoved = false;
    for (std::size_t v = 0; v < mesh.vertices.size(); v++)
    {
        const Vertex &vertex = mesh.vertices[v];
        bool moved = false;
        for (int axis = 0; axis < 3; axis++)
        {
            moved = moved || (vertex.coord.get(axis, 0) != before[v][axis]);
        }
        if (vertex.is_fixed())
        {
            EXPECT_FALSE(moved) << "clamped vertex " << v << " moved";
        }
        somethingMoved = somethingMoved || (moved && vertex.type == VertexType::Free);
    }
    EXPECT_TRUE(somethingMoved);
}

TEST(MixedBoundaryTest, TheDynamicsHoldsTheClampedVerticesAndCarriesTheImages)
{
    Param param;
    configure_sheet(param, BoundaryType::Periodic, BoundaryType::Fixed);
    DynamicMesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_flat());
    prepare(mesh, 1.0);
    mesh.update_vertices_mat_with_vector();

    Record record(16);
    record.add(mesh.param.area, mesh.param.energy, mesh.calculate_mean_force());
    DynamicModel model(mesh, record);

    mesh.ensure_surface_solver();
    EXPECT_EQ(mesh.surfaceSolver.nFree(), count_independent(mesh));
    // The slaved set the step skips is the images, and only the images.
    ASSERT_EQ(mesh.isSlavedPeriodic.size(), mesh.vertices.size());
    for (std::size_t v = 0; v < mesh.vertices.size(); v++)
    {
        EXPECT_EQ(mesh.isSlavedPeriodic[v] != 0, mesh.vertices[v].is_periodic_image()) << v;
    }

    std::vector<std::array<double, 3>> before(mesh.vertices.size());
    for (std::size_t v = 0; v < mesh.vertices.size(); v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            before[v][axis] = mesh.vertices[v].coord.get(axis, 0);
        }
    }

    for (model.iteration = 0; model.iteration < 5; model.iteration++)
    {
        mesh.apply_mesh_to_surface();
        model.next_step();
        mesh.apply_surface_to_mesh();
        mesh.postprocess_boundary();
        mesh.update_vertices_vector_with_mat();
        record.add(mesh.param.area, mesh.param.energy, mesh.calculate_mean_force());
        mesh.Compute_Energy_And_Force();
        EXPECT_TRUE(images_are_in_place(mesh)) << "after step " << model.iteration;
        EXPECT_LT(mesh.surfaceSolver.lastResidual(), 1e-8) << "after step " << model.iteration;
    }
    EXPECT_TRUE(std::isfinite(mesh.param.energy.energyTotal));

    bool somethingMoved = false;
    for (std::size_t v = 0; v < mesh.vertices.size(); v++)
    {
        const Vertex &vertex = mesh.vertices[v];
        bool moved = false;
        for (int axis = 0; axis < 3; axis++)
        {
            moved = moved || (vertex.coord.get(axis, 0) != before[v][axis]);
        }
        if (vertex.is_fixed())
        {
            EXPECT_FALSE(moved) << "clamped vertex " << v << " moved";
        }
        somethingMoved = somethingMoved || (moved && vertex.type == VertexType::Free);
    }
    EXPECT_TRUE(somethingMoved);
}

/**
 * @brief The wrapped solver: M C = S with images substituted, both ways.
 *
 * The mask row of a source next to the seam reaches images, and the solver
 * must read those as the source one period away plus the offset. Applying
 * the mask and solving back has to return the control net it was given,
 * images included, and the images of the surface are their sources' plus
 * the same offset.
 */
TEST(MixedBoundaryTest, TheWrappedSurfaceSolverRoundTripsThroughTheImages)
{
    Param param;
    configure_sheet(param, BoundaryType::Periodic, BoundaryType::Periodic);
    DynamicMesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_flat());
    undulate(mesh, 1.5);
    mesh.update_vertices_mat_with_vector();
    mesh.ensure_surface_solver();

    const int n = static_cast<int>(mesh.vertices.size());
    EXPECT_EQ(mesh.surfaceSolver.nFree(), (param.nFaceX - 6) * (param.nFaceY - 6));

    Matrix surface = mat_calloc(n, 3);
    mesh.surfaceSolver.mesh_to_surface(mesh.matMesh, surface);
    for (int v : mesh.periodicImageVertices)
    {
        const Vertex &image = mesh.vertices[v];
        for (int axis = 0; axis < 3; axis++)
        {
            EXPECT_DOUBLE_EQ(surface(v, axis),
                             surface(image.reflectiveVertexIndex, axis) + image.mirrorOffset[axis])
                << "image " << v << " axis " << axis;
        }
    }

    Matrix control = mat_calloc(n, 3);
    const int iterations = mesh.surfaceSolver.surface_to_mesh(surface, control);
    EXPECT_GT(iterations, 0);
    EXPECT_LT(mesh.surfaceSolver.lastResidual(), 1e-9);
    double worst = 0.0;
    for (int v = 0; v < n; v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            worst = std::max(worst, std::abs(control(v, axis) - mesh.matMesh(v, axis)));
        }
    }
    EXPECT_LT(worst, 1e-7) << "solving M C = S did not return the control net M was applied to";
}

// ---------------------------------------------------------------------------
// The mesh files
// ---------------------------------------------------------------------------

TEST(MixedBoundaryTest, ExportingAndImportingReproducesTheMesh)
{
    const std::string verticesPath = scratch("vertices.csv");
    const std::string facesPath = scratch("faces.csv");

    Param sourceParam;
    configure_sheet(sourceParam, BoundaryType::Periodic, BoundaryType::Fixed);
    Mesh source(sourceParam);
    ASSERT_NO_THROW(source.setup_flat());
    prepare(source, 1.0);
    ASSERT_NO_THROW(export_mesh_to_vertices_faces(source, verticesPath, facesPath));

    // The loaded mesh knows nothing about the grid: no sides, no axis types.
    Param loadedParam;
    loadedParam.VERBOSE_MODE = false;
    loadedParam.boundaryCondition = BoundaryType::Mixed;
    loadedParam.kCurv = sourceParam.kCurv;
    loadedParam.uSurf = sourceParam.uSurf;
    loadedParam.uVol = 0.0;
    loadedParam.isGlobalConstraint = true;
    loadedParam.usingRpi = false;
    loadedParam.isEnergyHarmonicBondIncluded = false;
    Mesh loaded(loadedParam);
    ASSERT_TRUE(import_mesh_from_vertices_faces(loaded, verticesPath, facesPath));

    ASSERT_EQ(loaded.vertices.size(), source.vertices.size());
    ASSERT_EQ(loaded.faces.size(), source.faces.size());
    for (std::size_t v = 0; v < source.vertices.size(); v++)
    {
        const Vertex &a = source.vertices[v];
        const Vertex &b = loaded.vertices[v];
        EXPECT_EQ(a.type, b.type) << "vertex " << v;
        EXPECT_EQ(a.reflectiveVertexIndex, b.reflectiveVertexIndex) << "vertex " << v;
        EXPECT_EQ(a.isGhost, b.isGhost) << "vertex " << v;
        for (int axis = 0; axis < 3; axis++)
        {
            EXPECT_EQ(a.coord.get(axis, 0), b.coord.get(axis, 0)) << "vertex " << v;
            EXPECT_EQ(a.mirrorOffset[axis], b.mirrorOffset[axis]) << "vertex " << v;
        }
    }
    for (std::size_t f = 0; f < source.faces.size(); f++)
    {
        EXPECT_EQ(source.faces[f].adjacentVertices, loaded.faces[f].adjacentVertices) << "face " << f;
        EXPECT_EQ(source.faces[f].isGhost, loaded.faces[f].isGhost) << "face " << f;
        EXPECT_EQ(source.faces[f].patchKind, loaded.faces[f].patchKind) << "face " << f;
    }
    EXPECT_EQ(loaded.periodicImageVertices, source.periodicImageVertices);

    // And it is the same membrane: same bending energy, same area. (The
    // regularization is not compared: the loaded mesh's reference
    // configuration is the relief it was written in, the source's is flat.)
    loaded.param.area0 = source.param.area0;
    loaded.update_previous_coord_for_vertex();
    loaded.update_reference_coord_from_previous_coord();
    loaded.calculate_element_area_volume();
    loaded.Compute_Energy_And_Force();
    source.Compute_Energy_And_Force();
    EXPECT_NEAR(bending_energy(loaded), bending_energy(source), 1e-9 * std::abs(bending_energy(source)));
    EXPECT_NEAR(loaded.param.area, source.param.area, 1e-9 * source.param.area);

    std::remove(verticesPath.c_str());
    std::remove(facesPath.c_str());
}

/**
 * @brief Without a copy-flag column the copies are derived from the mirror
 * map, and the derivation tiles the period exactly once.
 *
 * Which of two copies astride a seam keeps the energy is a convention, so the
 * derived layout need not match the generator's; what must hold is that every
 * physical face is charged once, and then the energy is the same.
 */
TEST(MixedBoundaryTest, CopiesDerivedFromTheMirrorMapTileThePeriodOnce)
{
    const std::string verticesPath = scratch("derived_vertices.csv");
    const std::string flaggedFacesPath = scratch("derived_faces_flagged.csv");
    const std::string bareFacesPath = scratch("derived_faces.csv");

    Param sourceParam;
    configure_sheet(sourceParam, BoundaryType::Periodic, BoundaryType::Periodic);
    Mesh source(sourceParam);
    ASSERT_NO_THROW(source.setup_flat());
    prepare(source, 1.0);
    ASSERT_NO_THROW(export_mesh_to_vertices_faces(source, verticesPath, flaggedFacesPath));
    // Three columns only: the legacy face.csv layout.
    source.write_faces_csv(bareFacesPath);

    Param loadedParam;
    loadedParam.VERBOSE_MODE = false;
    loadedParam.boundaryCondition = BoundaryType::Mixed;
    loadedParam.kCurv = sourceParam.kCurv;
    loadedParam.uSurf = sourceParam.uSurf;
    loadedParam.uVol = 0.0;
    loadedParam.isGlobalConstraint = true;
    loadedParam.usingRpi = false;
    loadedParam.isEnergyHarmonicBondIncluded = false;
    Mesh loaded(loadedParam);
    ASSERT_TRUE(import_mesh_from_vertices_faces(loaded, verticesPath, bareFacesPath));

    const int periodX = sourceParam.nFaceX - 6;
    const int periodY = sourceParam.nFaceY - 6;
    int nReal = 0;
    std::set<std::array<int, 3>> physicalFaces;
    for (const Face &face : loaded.faces)
    {
        if (face.isGhost)
        {
            continue;
        }
        nReal++;
        std::array<int, 3> key;
        for (int k = 0; k < 3; k++)
        {
            const Vertex &corner = loaded.vertices[face.adjacentVertices[k]];
            key[k] = corner.is_periodic_image() ? corner.reflectiveVertexIndex : corner.index;
        }
        std::sort(key.begin(), key.end());
        EXPECT_TRUE(physicalFaces.insert(key).second) << "face " << face.index << " is charged twice";
    }
    EXPECT_EQ(nReal, 2 * periodX * periodY);

    loaded.param.area0 = source.param.area0;
    loaded.update_previous_coord_for_vertex();
    loaded.update_reference_coord_from_previous_coord();
    loaded.calculate_element_area_volume();
    loaded.Compute_Energy_And_Force();
    source.Compute_Energy_And_Force();
    EXPECT_NEAR(bending_energy(loaded), bending_energy(source), 1e-9 * std::abs(bending_energy(source)));
    EXPECT_NEAR(loaded.param.area, source.param.area, 1e-9 * source.param.area);

    std::remove(verticesPath.c_str());
    std::remove(flaggedFacesPath.c_str());
    std::remove(bareFacesPath.c_str());
}

/**
 * @brief The example mesh shipped under data/example loads, and is what the
 * documentation says it is: a 60 nm sheet periodic along x and clamped along y.
 *
 * Written by the dynamics driver itself, so this also pins the exporter's
 * format against a file that is not regenerated by the test.
 */
TEST(MixedBoundaryTest, TheShippedExampleMeshLoads)
{
    Param param;
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Mixed;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
    Mesh mesh(param);
    ASSERT_TRUE(import_mesh_from_vertices_faces(mesh, "./data/example/mixed_sheet_vertices.csv",
                                                "./data/example/mixed_sheet_faces.csv"));
    EXPECT_EQ(mesh.vertices.size(), 195u);
    EXPECT_EQ(mesh.faces.size(), 336u);

    int nFree = 0;
    int nFixed = 0;
    int nImages = 0;
    for (const Vertex &vertex : mesh.vertices)
    {
        nFree += (vertex.type == VertexType::Free) ? 1 : 0;
        nFixed += vertex.is_fixed() ? 1 : 0;
        nImages += vertex.is_periodic_image() ? 1 : 0;
        if (vertex.is_periodic_image())
        {
            // Periodic along x only: every offset is one period along x.
            EXPECT_NEAR(std::abs(vertex.mirrorOffset[0]), 30.0, 1e-9) << "vertex " << vertex.index;
            EXPECT_EQ(vertex.mirrorOffset[1], 0.0) << "vertex " << vertex.index;
            EXPECT_EQ(vertex.mirrorOffset[2], 0.0) << "vertex " << vertex.index;
        }
    }
    EXPECT_EQ(nFree, 78);
    EXPECT_EQ(nFixed, 12);
    EXPECT_EQ(nImages, 105);
    EXPECT_EQ(static_cast<int>(mesh.periodicImageVertices.size()), nImages);

    int nCopies = 0;
    for (const Face &face : mesh.faces)
    {
        nCopies += face.isGhost ? 1 : 0;
    }
    EXPECT_EQ(nCopies, 168);
    EXPECT_TRUE(images_are_in_place(mesh));

    // It carries a limit surface: a flat sheet's area, no bending.
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    EXPECT_GT(mesh.param.area0, 1000.0);
    mesh.param.vol0 = 0.0;
    mesh.Compute_Energy_And_Force();
    EXPECT_NEAR(bending_energy(mesh), 0.0, 1e-9);
}

/// A file with no type column loads under the global modes as it always has.
TEST(MixedBoundaryTest, AnUntypedFileStillLoadsUnderTheGlobalModes)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    EXPECT_EQ(param.boundaryCondition, BoundaryType::Periodic);
    ASSERT_TRUE(import_mesh_from_vertices_faces(mesh, "./data/example/vertices_flat.csv",
                                                "./data/example/faces_flat.csv"));
    EXPECT_EQ(mesh.vertices.size(), 1927u);
    EXPECT_EQ(mesh.faces.size(), 3680u);
    EXPECT_TRUE(mesh.periodicImageVertices.empty());
    for (const Vertex &vertex : mesh.vertices)
    {
        EXPECT_EQ(vertex.type, VertexType::Free);
    }
}

TEST(MixedBoundaryTest, MalformedMeshesAreRefusedWithAReason)
{
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;
    octahedron(vertices, faces);

    // Types under a global mode: refused rather than ignored.
    {
        Param param;
        param.VERBOSE_MODE = false;
        param.boundaryCondition = BoundaryType::Periodic;
        Mesh mesh(param);
        std::vector<VertexType> types(6, VertexType::Free);
        types[2] = VertexType::Fixed;
        EXPECT_THROW(mesh.setup_from_vertices_faces(vertices, faces, types, std::vector<int>(6, -1),
                                                    std::vector<char>()),
                     std::runtime_error);
    }
    // All-free types under a global mode carry no information and are fine.
    {
        Param param;
        param.VERBOSE_MODE = false;
        param.boundaryCondition = BoundaryType::Fixed;
        Mesh mesh(param);
        EXPECT_NO_THROW(mesh.setup_from_vertices_faces(vertices, faces,
                                                       std::vector<VertexType>(6, VertexType::Free),
                                                       std::vector<int>(6, -1), std::vector<char>()));
    }
    // A mirror that names no vertex, a vertex that mirrors itself, a cycle.
    for (int broken = 0; broken < 3; broken++)
    {
        Param param;
        param.VERBOSE_MODE = false;
        param.boundaryCondition = BoundaryType::Mixed;
        Mesh mesh(param);
        std::vector<VertexType> types(6, VertexType::Free);
        std::vector<int> mirrors(6, -1);
        if (broken == 0)
        {
            types[4] = VertexType::Periodic;
            mirrors[4] = 99;
        }
        else if (broken == 1)
        {
            types[4] = VertexType::Periodic;
            mirrors[4] = 4;
        }
        else
        {
            types[4] = VertexType::Periodic;
            mirrors[4] = 5;
            types[5] = VertexType::Periodic;
            mirrors[5] = 4;
        }
        EXPECT_THROW(mesh.setup_from_vertices_faces(vertices, faces, types, mirrors, std::vector<char>()),
                     std::invalid_argument)
            << "case " << broken;
    }
    // A clamped vertex on a closed surface is an ordinary Mixed mesh.
    {
        Param param;
        param.VERBOSE_MODE = false;
        param.boundaryCondition = BoundaryType::Mixed;
        Mesh mesh(param);
        std::vector<VertexType> types(6, VertexType::Free);
        types[0] = VertexType::Fixed;
        ASSERT_NO_THROW(mesh.setup_from_vertices_faces(vertices, faces, types, std::vector<int>(6, -1),
                                                       std::vector<char>()));
        EXPECT_TRUE(mesh.periodicImageVertices.empty());
        EXPECT_EQ(count_independent(mesh), 5);
        mesh.param.uVol = 1.0;
        EXPECT_NO_THROW(mesh.validate_volume_constraint_topology());
    }
    // A periodic axis too short for its band.
    {
        Param param;
        configure_sheet(param, BoundaryType::Periodic, BoundaryType::Free, 25.0);
        Mesh mesh(param);
        EXPECT_THROW(mesh.setup_flat(), std::runtime_error);
    }
    // Only one of the two mesh files named.
    {
        Param param;
        configure_sheet(param, BoundaryType::Periodic, BoundaryType::Periodic);
        param.meshVerticesFile = scratch("only_vertices.csv");
        Mesh mesh(param);
        EXPECT_THROW(setup_mesh_from_parameters(mesh), std::runtime_error);
    }
    // A periodic vertex with no mirror in the file, and a typed file under a
    // global mode, both refused when read.
    {
        const std::string verticesPath = scratch("bad_vertices.csv");
        const std::string facesPath = scratch("bad_faces.csv");
        {
            std::ofstream out(verticesPath);
            out << "# a periodic vertex without a mirror\n";
            out << "0, 0, 0, free\n1, 0, 0, free\n0, 1, 0, periodic\n";
            std::ofstream faceOut(facesPath);
            faceOut << "0, 1, 2\n";
        }
        EXPECT_THROW(read_mesh_vertices_faces_files(verticesPath, facesPath), std::invalid_argument);
        {
            std::ofstream out(verticesPath);
            out << "0, 0, 0, fixed\n1, 0, 0\n0, 1, 0\n";
        }
        MeshFileData data;
        ASSERT_NO_THROW(data = read_mesh_vertices_faces_files(verticesPath, facesPath));
        EXPECT_TRUE(data.has_boundary_types());
        EXPECT_EQ(data.types.size(), 3u);
        EXPECT_EQ(data.types[0], VertexType::Fixed);
        EXPECT_EQ(data.types[1], VertexType::Free);
        EXPECT_TRUE(data.faceIsCopy.empty());
        Param param;
        param.VERBOSE_MODE = false;
        param.boundaryCondition = BoundaryType::Free;
        Mesh mesh(param);
        EXPECT_THROW(import_mesh_from_vertices_faces(mesh, verticesPath, facesPath), std::runtime_error);
        std::remove(verticesPath.c_str());
        std::remove(facesPath.c_str());
    }
}
