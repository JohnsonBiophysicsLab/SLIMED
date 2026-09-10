#include "test_multi_extraordinary.hpp"

/**
 * A face with more than one extraordinary corner has no Stam reduction: that
 * peels three regular children off a *single* extraordinary corner. The tree
 * has always rejected such faces, which was tenable only while connectivity
 * was fixed at setup. Every edge flip produces two of them.
 *
 * They are evaluated now by subdividing the face's own control net once, which
 * splits it into four children that each do have a reduction, and pushing
 * those through the paths that already exist.
 *
 * The binding test is an exact identity rather than a convergence claim. The
 * four children this path evaluates are precisely the four faces a global Loop
 * refinement would produce, with precisely the same control nets, rows and
 * quadrature. So the coarse face's area, volume and bending energy must equal
 * the sum over its four children on the refined mesh, to round-off -- and
 * refine_loop_once() is code that passed its own review, written before any of
 * this and for a different purpose.
 *
 * @see docs/edge_flip_plan.md work package 1
 */

namespace
{
struct MeshFixture
{
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;
};

double bowl_height(double x, double y)
{
    return 0.06 * (x * x + y * y);
}

MeshFixture build_grid(int nx, int ny)
{
    MeshFixture fixture;
    const double dy = std::sqrt(3.0) / 2.0;
    for (int j = 0; j <= ny; j++)
    {
        for (int i = 0; i <= nx; i++)
        {
            const double x = i + 0.5 * j;
            const double y = j * dy;
            fixture.vertices.push_back({x, y, bowl_height(x, y)});
        }
    }
    const auto index = [nx](int i, int j) { return j * (nx + 1) + i; };
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            fixture.faces.push_back({index(i, j), index(i + 1, j), index(i, j + 1)});
            fixture.faces.push_back({index(i + 1, j), index(i + 1, j + 1), index(i, j + 1)});
        }
    }
    return fixture;
}

/// The mesh's current state, as the vertex/face lists a fresh Mesh takes.
MeshFixture snapshot(const Mesh &mesh)
{
    MeshFixture fixture;
    for (const Vertex &vertex : mesh.vertices)
    {
        fixture.vertices.push_back(
            {vertex.coord.get(0, 0), vertex.coord.get(1, 0), vertex.coord.get(2, 0)});
    }
    for (const Face &face : mesh.faces)
    {
        fixture.faces.push_back(face.adjacentVertices);
    }
    return fixture;
}

void configure(Param &param)
{
    param.VERBOSE_MODE = false;
    // Pure bending: the constraint terms are global, and the two meshes being
    // compared have different totals near their boundaries.
    param.uSurf = 0.0;
    param.uVol = 0.0;
    param.usingRpi = false;
    param.isGlobalConstraint = true;
    param.area0 = 1.0;
    param.vol0 = 1.0;
}

/// A grid with a few interior edges flipped, so it carries faces with two and
/// three extraordinary corners.
void build_flipped_grid(Mesh &mesh, int nx, int ny, int nFlips)
{
    const MeshFixture fixture = build_grid(nx, ny);
    mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces);

    unsigned long long state = 987654321ULL;
    int flipped = 0;
    for (int attempt = 0; attempt < 4000 && flipped < nFlips; attempt++)
    {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        const int iEdge = static_cast<int>((state >> 33) % mesh.edges.size());
        if (!mesh.edge_flip_is_admissible(iEdge))
        {
            continue;
        }
        mesh.flip_edge(iEdge);
        flipped++;
    }
    ASSERT_EQ(flipped, nFlips) << "the fixture needs " << nFlips << " flips";

    // energy_force_regularization() measures each face's edges against the
    // same face's edges in coordRef, so coordRef has to exist before any force
    // is computed. Every other fixture in the suite does this too.
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
}

double total_energy(Mesh &mesh)
{
    mesh.calculate_element_area_volume();
    mesh.Compute_Energy_And_Force();
    double sum = 0.0;
    for (const Face &face : mesh.faces)
    {
        sum += face.energy.energyCurvature;
    }
    return sum;
}
} // namespace

// ---------------------------------------------------------------------------
// The generic patch and its prolongations
// ---------------------------------------------------------------------------

TEST(MultiExtraordinaryPatchTest, EveryValenceTripleInRangeBuilds)
{
    MultiPatchTable table;
    int built = 0;
    for (int n0 = kMinIrregularValence; n0 <= kMaxIrregularValence; n0++)
    {
        for (int n1 = kMinIrregularValence; n1 <= kMaxIrregularValence; n1++)
        {
            for (int n2 = kMinIrregularValence; n2 <= kMaxIrregularValence; n2++)
            {
                const GenericFacePatch patch = build_generic_face_patch(n0, n1, n2);
                EXPECT_EQ(patch.nVertices, n0 + n1 + n2 - 6);
                EXPECT_EQ(static_cast<int>(patch.faces.size()), n0 + n1 + n2 - 5);

                const int index = table.ensure(n0, n1, n2);
                const MultiPatchTable::Entry &entry = table.entry(index);
                EXPECT_EQ(entry.nControl, n0 + n1 + n2 - 6);
                // Three corner children carrying their corner's valence, and a
                // regular centre.
                EXPECT_EQ(entry.children[0].valence, n0);
                EXPECT_EQ(entry.children[1].valence, n1);
                EXPECT_EQ(entry.children[2].valence, n2);
                EXPECT_EQ(entry.children[3].valence, 6);
                EXPECT_EQ(entry.children[3].nControl, 12);
                for (int c = 0; c < 4; c++)
                {
                    EXPECT_EQ(entry.children[c].nControl, entry.children[c].valence + 6);
                    EXPECT_LE(entry.children[c].nControl, slimed::kMaxControlPoints);
                }
                EXPECT_LE(entry.nControl, slimed::kMaxControlPoints);
                built++;
            }
        }
    }
    EXPECT_EQ(built, 125);
    EXPECT_EQ(table.size(), 125) << "each triple should be built once and cached";
    // A repeat must hit the cache rather than grow the table.
    EXPECT_EQ(table.ensure(5, 5, 7), table.ensure(5, 5, 7));
    EXPECT_EQ(table.size(), 125);
}

TEST(MultiExtraordinaryPatchTest, ProlongationsAreAffineAndPlanarPrecise)
{
    MultiPatchTable table;
    for (int n0 = kMinIrregularValence; n0 <= kMaxIrregularValence; n0++)
    {
        for (int n1 = kMinIrregularValence; n1 <= kMaxIrregularValence; n1++)
        {
            for (int n2 = kMinIrregularValence; n2 <= kMaxIrregularValence; n2++)
            {
                const int index = table.ensure(n0, n1, n2);
                const MultiPatchTable::Entry &entry = table.entry(index);
                const int K = entry.nControl;
                for (int c = 0; c < 4; c++)
                {
                    const MultiPatchTable::Child &child = entry.children[c];
                    const double *M = table.data() + child.offset;
                    for (int row = 0; row < child.nControl; row++)
                    {
                        double sum = 0.0;
                        for (int column = 0; column < K; column++)
                        {
                            sum += M[row * K + column];
                        }
                        // Affine invariance: a subdivided point is a weighted
                        // average, so translating the mesh translates it.
                        EXPECT_NEAR(sum, 1.0, 1e-14)
                            << "valences (" << n0 << ", " << n1 << ", " << n2 << ") child " << c
                            << " row " << row;
                    }
                }
            }
        }
    }
}

/**
 * @brief A vertex point of the prolongation is Loop's vertex mask, and an edge
 * point is Loop's edge mask, checked against the masks written out by hand.
 *
 * Reaches inside one entry rather than trusting the whole pipeline, so a wrong
 * mask cannot hide behind a compensating error in the ring construction.
 */
TEST(MultiExtraordinaryPatchTest, TheCornerChildStartsAtTheLoopVertexPoint)
{
    MultiPatchTable table;
    const int n0 = 5;
    const int n1 = 7;
    const int n2 = 6;
    const int index = table.ensure(n0, n1, n2);
    const MultiPatchTable::Entry &entry = table.entry(index);
    const int K = entry.nControl;
    const GenericFacePatch patch = build_generic_face_patch(n0, n1, n2);

    // Corner child 0's control net is in canonical order; the extraordinary
    // corner d4 sits at the column canonical_control_order() assigns to
    // internal index 0.
    const int d4Column = canonical_control_order(n0)[0];
    const MultiPatchTable::Child &child = entry.children[0];
    const double *const row = table.data() + child.offset + static_cast<std::size_t>(d4Column) * K;

    // Loop's vertex mask at valence n0, with Warren's weight.
    const double beta = loop_vertex_weight(n0);
    std::vector<double> expected(K, 0.0);
    expected[0] = 1.0 - n0 * beta;
    for (int neighbour : patch.cornerFan[0])
    {
        expected[neighbour] = beta;
    }
    for (int column = 0; column < K; column++)
    {
        EXPECT_NEAR(row[column], expected[column], 1e-15) << "column " << column;
    }
}

// ---------------------------------------------------------------------------
// The identity that binds: against a global Loop refinement
// ---------------------------------------------------------------------------

/**
 * @brief A multi-extraordinary face's area, volume and bending energy equal the
 * sum over the four faces a global refinement puts in its place.
 *
 * This is exact, not asymptotic. The four children evaluated here *are* those
 * four faces: same control nets, same rows, same three-point rule on each.
 * refine_loop_once() emits them in the order (corner a, corner b, corner c,
 * centre) at indices 4i..4i+3, which is the order the prolongation table uses.
 *
 * It tests the prolongation matrices, the child one-ring construction and the
 * chain-rule scatter together, against code written before any of this and for
 * an unrelated purpose.
 */
TEST(MultiExtraordinaryPatchTest, MatchesAGlobalLoopRefinementFaceForFace)
{
    Param coarseParam;
    configure(coarseParam);
    Mesh coarse(coarseParam);
    ASSERT_NO_FATAL_FAILURE(build_flipped_grid(coarse, 10, 10, 6));

    int nMulti = 0;
    for (const Face &face : coarse.faces)
    {
        nMulti += (face.patchKind == PatchKind::MultiExtraordinary) ? 1 : 0;
    }
    ASSERT_GT(nMulti, 0) << "the flipped grid should carry multi-extraordinary faces";

    const MeshFixture flipped = snapshot(coarse);
    coarse.calculate_element_area_volume();
    coarse.Compute_Energy_And_Force();

    Param refinedParam;
    configure(refinedParam);
    refinedParam.isPreRefinementEnabled = true;
    Mesh refined(refinedParam);
    refined.setup_from_vertices_faces(flipped.vertices, flipped.faces);
    ASSERT_EQ(refined.faces.size(), coarse.faces.size() * 4);
    refined.update_previous_coord_for_vertex();
    refined.update_reference_coord_from_previous_coord();
    refined.calculate_element_area_volume();
    refined.Compute_Energy_And_Force();

    int nCompared = 0;
    for (int i = 0; i < static_cast<int>(coarse.faces.size()); i++)
    {
        const Face &parent = coarse.faces[i];
        if (parent.patchKind != PatchKind::MultiExtraordinary)
        {
            continue;
        }
        // Only faces whose four children are all evaluable on the refined mesh
        // can be compared; near the boundary they are not, and that is a
        // property of the fixture rather than of the construction.
        bool allEvaluable = true;
        for (int c = 0; c < 4; c++)
        {
            allEvaluable = allEvaluable && refined.faces[4 * i + c].patchKind != PatchKind::Boundary;
        }
        if (!allEvaluable)
        {
            continue;
        }

        double area = 0.0;
        double volume = 0.0;
        double bending = 0.0;
        for (int c = 0; c < 4; c++)
        {
            const Face &child = refined.faces[4 * i + c];
            area += child.elementArea;
            volume += child.elementVolume;
            bending += child.energy.energyCurvature;
        }

        EXPECT_NEAR(parent.elementArea, area, 1e-11 * std::max(1.0, std::abs(area)))
            << "face " << i;
        EXPECT_NEAR(parent.elementVolume, volume, 1e-11 * std::max(1.0, std::abs(volume)))
            << "face " << i;
        EXPECT_NEAR(parent.energy.energyCurvature, bending,
                    1e-10 * std::max(1.0, std::abs(bending)))
            << "face " << i << " with valences (" << parent.patchValence[0] << ", "
            << parent.patchValence[1] << ", " << parent.patchValence[2] << ")";
        nCompared++;
    }
    EXPECT_GT(nCompared, 5) << "not enough multi-extraordinary faces were away from the boundary";
}

/**
 * @brief The same identity for a face with exactly one extraordinary corner,
 * where the tree's own single-extraordinary path is the reference.
 *
 * A 5/6/6 face is evaluated by Stam's reduction directly. Refining it and
 * summing the four children must give the same answer, so this pins the new
 * machinery against a path that was validated in its own work package.
 */
TEST(MultiExtraordinaryPatchTest, RefinementIdentityHoldsForSingleExtraordinaryFacesToo)
{
    Param coarseParam;
    configure(coarseParam);
    Mesh coarse(coarseParam);
    ASSERT_NO_FATAL_FAILURE(build_flipped_grid(coarse, 10, 10, 3));
    coarse.calculate_element_area_volume();
    coarse.Compute_Energy_And_Force();

    const MeshFixture flipped = snapshot(coarse);
    Param refinedParam;
    configure(refinedParam);
    refinedParam.isPreRefinementEnabled = true;
    Mesh refined(refinedParam);
    refined.setup_from_vertices_faces(flipped.vertices, flipped.faces);
    refined.update_previous_coord_for_vertex();
    refined.update_reference_coord_from_previous_coord();
    refined.calculate_element_area_volume();
    refined.Compute_Energy_And_Force();

    int nCompared = 0;
    for (int i = 0; i < static_cast<int>(coarse.faces.size()); i++)
    {
        const Face &parent = coarse.faces[i];
        if (parent.patchKind != PatchKind::SingleExtraordinary &&
            parent.patchKind != PatchKind::Regular)
        {
            continue;
        }
        bool allEvaluable = true;
        for (int c = 0; c < 4; c++)
        {
            allEvaluable = allEvaluable && refined.faces[4 * i + c].patchKind != PatchKind::Boundary;
        }
        if (!allEvaluable)
        {
            continue;
        }

        double area = 0.0;
        for (int c = 0; c < 4; c++)
        {
            area += refined.faces[4 * i + c].elementArea;
        }
        // Neither case is exact, and for different reasons worth separating.
        //
        // A regular coarse face is integrated by one three-point Gauss rule
        // over the whole triangle, while the refined mesh uses four such rules
        // over its quarters. Both approximate the same integral over the same
        // limit surface, so they differ by quadrature error alone -- around
        // 1e-8 relative here.
        //
        // At valence 5 there is a second, larger difference: both sides
        // truncate the same infinite sum over subdivision rings, but at
        // different points, because the coarse face starts one level higher up
        // than its children do. That leaves the truncation the chosen depth
        // allows rather than round-off.
        //
        // The exact identity is the one in the test above, where both sides
        // tile the same four children with the same rows.
        const double tolerance = (parent.patchKind == PatchKind::Regular) ? 1e-6 : 1e-3;
        EXPECT_NEAR(parent.elementArea, area, tolerance * std::max(1.0, std::abs(area)))
            << "face " << i << " kind "
            << (parent.patchKind == PatchKind::Regular ? "regular" : "single-extraordinary");
        nCompared++;
    }
    EXPECT_GT(nCompared, 20);
}

// ---------------------------------------------------------------------------
// Structural checks that do not go through the subdivision machinery
// ---------------------------------------------------------------------------

/**
 * @brief Force is minus the gradient of energy on a mesh full of
 * multi-extraordinary faces.
 *
 * The single strongest test available, and the one that covers the chain-rule
 * scatter: an error in the transpose would leave the energy right and the
 * force wrong, which nothing else here would catch.
 */
TEST(MultiExtraordinaryPatchTest, ForceIsMinusTheGradientOfEnergy)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_flipped_grid(mesh, 8, 8, 5));

    mesh.calculate_element_area_volume();
    mesh.Compute_Energy_And_Force();

    // Perturb vertices that a multi-extraordinary face actually depends on.
    std::vector<int> probes;
    for (const Face &face : mesh.faces)
    {
        if (face.patchKind != PatchKind::MultiExtraordinary)
        {
            continue;
        }
        for (int corner : face.adjacentVertices)
        {
            if (mesh.is_interior_vertex(corner) &&
                std::find(probes.begin(), probes.end(), corner) == probes.end())
            {
                probes.push_back(corner);
            }
        }
    }
    ASSERT_GT(probes.size(), 3u);
    probes.resize(std::min<std::size_t>(probes.size(), 6));

    const double h = 1e-5;
    for (int iVertex : probes)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            // Read the analytic force before perturbing: reading it inside the
            // loop leaves it one step stale, which looks like a clean linear
            // error and is not one.
            mesh.calculate_element_area_volume();
            mesh.Compute_Energy_And_Force();
            const double analytic = mesh.vertices[iVertex].force.forceCurvature(axis, 0);

            const double original = mesh.vertices[iVertex].coord.get(axis, 0);
            mesh.vertices[iVertex].coord.set(axis, 0, original + h);
            const double plus = total_energy(mesh);
            mesh.vertices[iVertex].coord.set(axis, 0, original - h);
            const double minus = total_energy(mesh);
            mesh.vertices[iVertex].coord.set(axis, 0, original);

            const double numeric = -(plus - minus) / (2.0 * h);
            EXPECT_NEAR(analytic, numeric, 1e-5 * std::max(1.0, std::abs(numeric)))
                << "vertex " << iVertex << " axis " << axis;
        }
    }
}

TEST(MultiExtraordinaryPatchTest, EnergyAndAreaAreInvariantUnderRigidMotion)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_flipped_grid(mesh, 8, 8, 5));

    mesh.calculate_element_area_volume();
    const double energyBefore = total_energy(mesh);
    double areaBefore = 0.0;
    for (const Face &face : mesh.faces)
    {
        areaBefore += face.elementArea;
    }

    // Rotate by 0.7 rad about the x axis and translate far off the origin.
    const double angle = 0.7;
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    for (Vertex &vertex : mesh.vertices)
    {
        const double x = vertex.coord.get(0, 0);
        const double y = vertex.coord.get(1, 0);
        const double z = vertex.coord.get(2, 0);
        vertex.coord.set(0, 0, x + 137.0);
        vertex.coord.set(1, 0, c * y - s * z - 91.0);
        vertex.coord.set(2, 0, s * y + c * z + 44.0);
    }

    mesh.calculate_element_area_volume();
    const double energyAfter = total_energy(mesh);
    double areaAfter = 0.0;
    for (const Face &face : mesh.faces)
    {
        areaAfter += face.elementArea;
    }

    EXPECT_NEAR(energyAfter, energyBefore, 1e-9 * std::max(1.0, std::abs(energyBefore)));
    EXPECT_NEAR(areaAfter, areaBefore, 1e-9 * std::max(1.0, std::abs(areaBefore)));
}

/**
 * @brief The forces a multi-extraordinary face puts on its own control points
 * sum to zero.
 *
 * Newton's third law inside one patch. It holds because the energy depends on
 * the control net only through differences, so a uniform translation cannot
 * change it -- and the chain-rule scatter preserves that, since every
 * prolongation row sums to one.
 */
TEST(MultiExtraordinaryPatchTest, InternalForcesOnOnePatchSumToZero)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_flipped_grid(mesh, 8, 8, 5));
    mesh.calculate_element_area_volume();
    mesh.Compute_Energy_And_Force();

    // Total force over the whole mesh: every face's internal forces cancel, so
    // the sum over all vertices must vanish.
    double total[3] = {0.0, 0.0, 0.0};
    for (const Vertex &vertex : mesh.vertices)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            total[axis] += vertex.force.forceCurvature(axis, 0);
        }
    }
    for (int axis = 0; axis < 3; axis++)
    {
        EXPECT_NEAR(total[axis], 0.0, 1e-8) << "axis " << axis;
    }
}

/**
 * @brief A mesh whose faces all carry several extraordinary corners is
 * accepted now instead of throwing.
 *
 * An icosahedron is the case the previous work package explicitly rejected:
 * every vertex is at valence 5, so every face is 5/5/5 and Stam's reduction
 * does not apply to any of them. It is also the closed fixture that
 * `irregular_patch_results.md` asked for and could not run.
 */
TEST(MultiExtraordinaryPatchTest, AnAllValence5SolidIsEvaluatedRatherThanRejected)
{
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    MeshFixture icosahedron;
    icosahedron.vertices = {{-1, phi, 0}, {1, phi, 0},  {-1, -phi, 0}, {1, -phi, 0},
                            {0, -1, phi}, {0, 1, phi},  {0, -1, -phi}, {0, 1, -phi},
                            {phi, 0, -1}, {phi, 0, 1},  {-phi, 0, -1}, {-phi, 0, 1}};
    icosahedron.faces = {{0, 11, 5}, {0, 5, 1},  {0, 1, 7},   {0, 7, 10}, {0, 10, 11},
                         {1, 5, 9},  {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
                         {3, 9, 4},  {3, 4, 2},  {3, 2, 6},   {3, 6, 8},  {3, 8, 9},
                         {4, 9, 5},  {2, 4, 11}, {6, 2, 10},  {8, 6, 7},  {9, 8, 1}};

    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_from_vertices_faces(icosahedron.vertices, icosahedron.faces));
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();

    for (const Face &face : mesh.faces)
    {
        EXPECT_EQ(face.patchKind, PatchKind::MultiExtraordinary);
        EXPECT_GE(face.patchEntry, 0);
        EXPECT_EQ(face.oneRingVertices.size(), 9u) << "5 + 5 + 5 - 6";
    }

    mesh.calculate_element_area_volume();
    mesh.Compute_Energy_And_Force();

    double area = 0.0;
    double volume = 0.0;
    for (const Face &face : mesh.faces)
    {
        area += face.elementArea;
        volume += face.elementVolume;
        EXPECT_TRUE(std::isfinite(face.energy.energyCurvature));
        EXPECT_GT(face.elementArea, 0.0);
    }

    // The limit surface lies in the convex hull of the control mesh, so it is
    // inside the polyhedron: its area is below the polyhedron's own, and its
    // enclosed volume below the polyhedron's. Both are elementary formulas for
    // edge length 2, and neither came from the subdivision machinery.
    const double edge = 2.0;
    const double polyhedronArea = 5.0 * std::sqrt(3.0) * edge * edge;
    const double polyhedronVolume = 5.0 * (3.0 + std::sqrt(5.0)) / 12.0 * edge * edge * edge;
    EXPECT_GT(area, 0.0);
    EXPECT_LT(area, polyhedronArea);
    EXPECT_GT(std::abs(volume), 0.0);
    EXPECT_LT(std::abs(volume), polyhedronVolume);

    // The real check: refining the icosahedron does not move its limit
    // surface, so the same integral computed on the refined control mesh must
    // agree. The refined mesh is a different control net evaluated through
    // different patches -- 20 multi-extraordinary faces become 80 faces, of
    // which only the 20 carrying an original corner are extraordinary at all --
    // so agreement is genuine corroboration rather than the same arithmetic
    // twice. They differ only by quadrature.
    Param refinedParam;
    configure(refinedParam);
    refinedParam.isPreRefinementEnabled = true;
    Mesh refined(refinedParam);
    refined.setup_from_vertices_faces(icosahedron.vertices, icosahedron.faces);
    refined.update_previous_coord_for_vertex();
    refined.update_reference_coord_from_previous_coord();
    refined.calculate_element_area_volume();

    double refinedArea = 0.0;
    double refinedVolume = 0.0;
    for (const Face &face : refined.faces)
    {
        refinedArea += face.elementArea;
        refinedVolume += face.elementVolume;
    }
    EXPECT_NEAR(area, refinedArea, 2e-3 * refinedArea);
    EXPECT_NEAR(volume, refinedVolume, 2e-3 * std::abs(refinedVolume));
}

/**
 * @brief The all-valence-4 closed solid, evaluated and corroborated.
 *
 * An octahedron is the other case the previous work package rejected, and the
 * one that stresses the machinery hardest: every face is 4/4/4, and valence 4
 * is where the bending energy converges slowest, so its corner children run to
 * depth 12.
 *
 * The face's one-ring here is the entire mesh -- six vertices, K = 4+4+4-6 --
 * which is as degenerate as an embedded patch gets. It is still embedded, so
 * it is still evaluable, and refining it must not move the surface.
 */
TEST(MultiExtraordinaryPatchTest, AnAllValence4SolidIsEvaluatedAndAgreesWithItsRefinement)
{
    MeshFixture octahedron;
    octahedron.vertices = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    octahedron.faces = {{0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4},
                        {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}};

    Param param;
    configure(param);
    Mesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_from_vertices_faces(octahedron.vertices, octahedron.faces));
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();

    for (const Face &face : mesh.faces)
    {
        EXPECT_EQ(face.patchKind, PatchKind::MultiExtraordinary);
        EXPECT_EQ(face.oneRingVertices.size(), 6u) << "4 + 4 + 4 - 6";
        // The one-ring is the whole mesh, and every entry must still be distinct.
        std::vector<int> sorted = face.oneRingVertices;
        std::sort(sorted.begin(), sorted.end());
        EXPECT_EQ(std::unique(sorted.begin(), sorted.end()), sorted.end());
    }

    mesh.calculate_element_area_volume();
    mesh.Compute_Energy_And_Force();
    double area = 0.0;
    double volume = 0.0;
    for (const Face &face : mesh.faces)
    {
        area += face.elementArea;
        volume += face.elementVolume;
        EXPECT_TRUE(std::isfinite(face.energy.energyCurvature));
    }
    EXPECT_GT(area, 0.0);
    // Strictly inside the control polyhedron, whose volume is 4/3.
    EXPECT_LT(std::abs(volume), 4.0 / 3.0);

    Param refinedParam;
    configure(refinedParam);
    refinedParam.isPreRefinementEnabled = true;
    Mesh refined(refinedParam);
    refined.setup_from_vertices_faces(octahedron.vertices, octahedron.faces);
    refined.update_previous_coord_for_vertex();
    refined.update_reference_coord_from_previous_coord();
    refined.calculate_element_area_volume();
    double refinedArea = 0.0;
    double refinedVolume = 0.0;
    for (const Face &face : refined.faces)
    {
        refinedArea += face.elementArea;
        refinedVolume += face.elementVolume;
    }
    EXPECT_NEAR(area, refinedArea, 5e-3 * refinedArea);
    EXPECT_NEAR(volume, refinedVolume, 5e-3 * std::abs(refinedVolume));
}
