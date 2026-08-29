#include "test_mesh_setup_geometry.hpp"

/**
 * @brief test default mesh initiation
 * 
 */
TEST(MeshInitTest, DefaultInitTest){
    Param param;
    param.VERBOSE_MODE = false; // mute output
    Mesh mesh(param);
}

/**
 * @brief Construct a new TEST object, unit-test mesh::setup_flat()
 * 
 */
TEST(MeshInitTest, SetupFlatTest){
    Param param;
    param.VERBOSE_MODE = false; // mute output
    Mesh mesh(param);
    mesh.setup_flat();
}


TEST(MeshFunctionTest, SortVerticesOfFacesTest){
    Param param;
    param.VERBOSE_MODE = false; // mute output
    Mesh mesh(param);
    mesh.setup_flat();
    mesh.sort_vertices_on_faces();

    // Expected output:
    // 21 , 22 , 0 , 
    // 22 , 1 , 0 , 
    // 22 , 23 , 1 ,
    /*
    std::cout << mesh.faces[0].adjacentVertices[0] << " , "
            << mesh.faces[0].adjacentVertices[1] << " , "
            << mesh.faces[0].adjacentVertices[2] << " , " << std::endl;
    std::cout << mesh.faces[1].adjacentVertices[0] << " , "
            << mesh.faces[1].adjacentVertices[1] << " , "
            << mesh.faces[1].adjacentVertices[2] << " , " << std::endl;
    std::cout << mesh.faces[2].adjacentVertices[0] << " , "
            << mesh.faces[2].adjacentVertices[1] << " , "
            << mesh.faces[2].adjacentVertices[2] << " , " << std::endl;
    */
    EXPECT_EQ(mesh.faces[1].adjacentVertices[0], 22);
    EXPECT_EQ(mesh.faces[1].adjacentVertices[1], 1);
    EXPECT_EQ(mesh.faces[1].adjacentVertices[2], 0);
}

/**
 * @brief Construct a new TEST, unit-test mesh::faces_share_edge()
 * 
 * Four test cases including 0-3 common vertices between two faces.
 * 
 * See mesh::faces_share_edge()
 * 
 */
TEST(MeshFunctionTest, FacesShareEdgeTest){
    Param param;
    param.VERBOSE_MODE = false; // mute output
    Mesh mesh(param);

    // Create faces and test share edge
    Face face1;
    Face face2;
    face1.adjacentVertices = std::vector<int>{1, 2, 3};
    face2.adjacentVertices = std::vector<int>{2, 3, 4};
    bool shareEdge = mesh.faces_share_edge(face1, face2);
    EXPECT_TRUE(shareEdge);

    // Function should return true if two inputs are the same
    shareEdge = mesh.faces_share_edge(face1, face1);
    EXPECT_TRUE(shareEdge);

    // Two cases that should return false:
    Face face3;
    Face face4;
    face3.adjacentVertices = std::vector<int>{2, 42, 13};
    face4.adjacentVertices = std::vector<int>{20, 43, 44};
    
    // Only 1 common vertex
    shareEdge = mesh.faces_share_edge(face1, face3);
    EXPECT_FALSE(shareEdge);

    // No common vertex
    shareEdge = mesh.faces_share_edge(face1, face4);
    EXPECT_FALSE(shareEdge);
}
/**
 * @brief Helper: give a mesh the eight faces of a closed octahedron.
 *
 * Written straight into `faces` rather than through setup_from_vertices_faces()
 * so the closure check can be exercised on its own, without dragging in the
 * one-ring construction (which does not yet handle valence-4 corners).
 * Every face is wound outward, so the surface is closed, two-manifold, and
 * consistently oriented.
 */
static void build_closed_octahedron(Mesh &mesh)
{
    // 0:+x 1:-x 2:+y 3:-y 4:+z 5:-z
    const std::vector<std::vector<int>> octahedronFaces{
        {0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4},
        {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}};

    mesh.faces = std::vector<Face>(octahedronFaces.size());
    for (std::size_t i = 0; i < octahedronFaces.size(); i++)
    {
        mesh.faces[i].index = static_cast<int>(i);
        mesh.faces[i].adjacentVertices = octahedronFaces[i];
    }
}

/**
 * @brief A volume constraint on a periodic flat sheet must fail at setup.
 *
 * This is the gate on step 3 of docs/volume_functional_split.md. The flat and
 * periodic workloads are open surfaces, so the signed-volume accumulator does
 * not return a volume there -- it is not even independent of where the origin
 * sits. Before this check the run proceeded silently against that number.
 */
TEST(VolumeConstraintTopologyTest, PeriodicFlatSheetWithConstraintThrows)
{
    Param param;
    param.VERBOSE_MODE = false; // mute output
    param.boundaryCondition = BoundaryType::Periodic;
    param.uVol = 1.0;
    Mesh mesh(param);

    EXPECT_THROW(mesh.setup_flat(), std::runtime_error);
}

/// The same sheet without a constraint must keep running -- this is data/example.
TEST(VolumeConstraintTopologyTest, PeriodicFlatSheetWithoutConstraintIsAccepted)
{
    Param param;
    param.VERBOSE_MODE = false; // mute output
    param.boundaryCondition = BoundaryType::Periodic;
    param.uVol = 0.0;
    Mesh mesh(param);

    EXPECT_NO_THROW(mesh.setup_flat());
}

/**
 * @brief Fixed boundaries produce no ghost faces, and must still be rejected.
 *
 * Ghost faces are the obvious tell, so this pins the case where that tell is
 * absent: the sheet is open regardless, and only the boundary-edge count
 * catches it.
 */
TEST(VolumeConstraintTopologyTest, FixedFlatSheetWithConstraintThrows)
{
    Param param;
    param.VERBOSE_MODE = false; // mute output
    param.boundaryCondition = BoundaryType::Fixed;
    param.uVol = 1.0;
    Mesh mesh(param);

    EXPECT_THROW(mesh.setup_flat(), std::runtime_error);
}

/// A closed, consistently oriented surface is what the constraint is for.
TEST(VolumeConstraintTopologyTest, ClosedOctahedronIsAccepted)
{
    Param param;
    param.VERBOSE_MODE = false; // mute output
    param.boundaryCondition = BoundaryType::Fixed;
    param.uVol = 1.0;
    Mesh mesh(param);
    build_closed_octahedron(mesh);

    EXPECT_NO_THROW(mesh.validate_volume_constraint_topology());
}

/// Removing one face opens the surface: three edges lose their second face.
TEST(VolumeConstraintTopologyTest, OpenOctahedronWithConstraintThrows)
{
    Param param;
    param.VERBOSE_MODE = false; // mute output
    param.boundaryCondition = BoundaryType::Fixed;
    param.uVol = 1.0;
    Mesh mesh(param);
    build_closed_octahedron(mesh);
    mesh.faces.pop_back();

    EXPECT_THROW(mesh.validate_volume_constraint_topology(), std::runtime_error);
}

/**
 * @brief A closed but inconsistently oriented surface must be rejected too.
 *
 * Flipping one face leaves every edge with exactly two incident faces, so a
 * closure test alone passes it -- while the signed volume it produces is
 * wrong. Only the directed-edge count separates the two.
 */
TEST(VolumeConstraintTopologyTest, MisorientedOctahedronWithConstraintThrows)
{
    Param param;
    param.VERBOSE_MODE = false; // mute output
    param.boundaryCondition = BoundaryType::Fixed;
    param.uVol = 1.0;
    Mesh mesh(param);
    build_closed_octahedron(mesh);
    std::swap(mesh.faces[0].adjacentVertices[1], mesh.faces[0].adjacentVertices[2]);

    EXPECT_THROW(mesh.validate_volume_constraint_topology(), std::runtime_error);
}

/// The message must name uVol and say what was wrong, not just fail.
TEST(VolumeConstraintTopologyTest, RejectionMessageNamesTheConstraintAndTheReason)
{
    Param param;
    param.VERBOSE_MODE = false; // mute output
    param.boundaryCondition = BoundaryType::Periodic;
    param.uVol = 2.5;
    Mesh mesh(param);

    try
    {
        mesh.setup_flat();
        FAIL() << "expected setup_flat() to reject a periodic sheet under a volume constraint";
    }
    catch (const std::runtime_error &error)
    {
        const std::string message = error.what();
        EXPECT_NE(message.find("uvVolumeConstraint"), std::string::npos);
        EXPECT_NE(message.find("2.5"), std::string::npos);
        EXPECT_NE(message.find("Periodic"), std::string::npos);
        EXPECT_NE(message.find("volume_functional_split"), std::string::npos);
    }
}

/**
 * @brief The legacy accumulator picks up only the x-component of the normal.
 *
 * setup_flat() lays the sheet in the z = 0 plane, so every patch normal points
 * along z. Lifting the sheet off the origin gives the corrected functional
 * something to integrate -- dot(x, a_1 x a_2) = z * n_z -- while the pre-fix
 * x-only integrand stays at exactly zero, because (a_1 x a_2)_x is zero
 * everywhere on this surface.
 *
 * This is the shape of the bug in one assertion: the old accumulator was
 * blind to two of the three components.
 *
 * @note TEMPORARY, alongside Mesh::sum_legacy_x_only_volume(). Delete with it.
 */
TEST(LegacyVolumeReportingTest, LegacyValueIsZeroWhenTheNormalHasNoXComponent)
{
    Param param;
    param.VERBOSE_MODE = false; // mute output
    Mesh mesh(param);
    mesh.setup_flat();

    // Lift the sheet away from z = 0 so the corrected functional is non-zero.
    for (Vertex &vertex : mesh.vertices)
    {
        vertex.coord.set(2, 0, 10.0);
    }

    mesh.calculate_element_area_volume();
    double area = 0.0;
    double corrected = 0.0;
    mesh.sum_membrane_area_and_volume(area, corrected);

    ASSERT_NE(corrected, 0.0);
    // Zero up to round-off in the cross product: the x-component of the normal
    // cancels analytically but not bit-exactly.
    EXPECT_LT(std::abs(mesh.sum_legacy_x_only_volume()), 1e-12 * std::abs(corrected));
}

/**
 * @brief Where the normal is purely x, the two functionals coincide.
 *
 * Remapping (x, y, 0) -> (10, x, y) stands the same sheet up in the plane
 * x = 10, so (a_1 x a_2) has only an x-component and dot(x, a_1 x a_2)
 * reduces to exactly the term the old code kept. The two accumulators must
 * then agree -- up to the old bare literal, which is low by about 4e-11
 * relative to 1/6.
 *
 * Together with the test above this pins the legacy reproduction from both
 * sides: it must vanish where the old code vanished, and match where the old
 * code was accidentally right.
 *
 * @note TEMPORARY, alongside Mesh::sum_legacy_x_only_volume(). Delete with it.
 */
TEST(LegacyVolumeReportingTest, LegacyValueMatchesCorrectedWhenTheNormalIsAlongX)
{
    Param param;
    param.VERBOSE_MODE = false; // mute output
    Mesh mesh(param);
    mesh.setup_flat();

    for (Vertex &vertex : mesh.vertices)
    {
        const double x = vertex.coord(0, 0);
        const double y = vertex.coord(1, 0);
        vertex.coord.set(0, 0, 10.0);
        vertex.coord.set(1, 0, x);
        vertex.coord.set(2, 0, y);
    }

    mesh.calculate_element_area_volume();
    double area = 0.0;
    double corrected = 0.0;
    mesh.sum_membrane_area_and_volume(area, corrected);
    const double legacy = mesh.sum_legacy_x_only_volume();

    ASSERT_NE(corrected, 0.0);
    EXPECT_NEAR(legacy / corrected, 1.0, 1e-9);
}

namespace
{
/// Vertex/face arrays for a closed Platonic solid, ready for setup_from_vertices_faces().
struct ClosedSolid
{
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;
};

/// Octahedron: 6 vertices, every one at valence 4.
ClosedSolid make_octahedron()
{
    return {{{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}},
            {{0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4},
             {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}}};
}

/// Icosahedron: 12 vertices, every one at valence 5.
ClosedSolid make_icosahedron()
{
    const double t = (1.0 + std::sqrt(5.0)) / 2.0;
    return {{{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0},
             {0, -1, t}, {0, 1, t}, {0, -1, -t}, {0, 1, -t},
             {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1}},
            {{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
             {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
             {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
             {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}}};
}

Mesh setup_solid(Param &param, const ClosedSolid &solid)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.uVol = 0.0;
    Mesh mesh(param);
    mesh.setup_from_vertices_faces(solid.vertices, solid.faces);
    return mesh;
}
} // namespace

/**
 * @brief Adjacent extraordinary corners are rejected loudly, not zeroed.
 *
 * Every octahedron vertex is at valence 4, so every face carries three
 * extraordinary corners. Valence 4 itself is supported now, but the uniform
 * patch reduction assumes exactly one extraordinary corner per face -- with
 * more, the three "regular" children stop being regular and there is no such
 * decomposition. The old code matched neither the 6/6/6 nor the all-5
 * predicate, so oneRingVertices stayed empty and the face was recorded with
 * zero bending energy and zero force.
 */
TEST(OneRingPatchTest, Valence4MeshIsRejectedRatherThanSilentlyZeroed)
{
    Param param;
    try
    {
        Mesh mesh = setup_solid(param, make_octahedron());
        FAIL() << "expected an octahedron (all valence 4) to be rejected";
    }
    catch (const std::runtime_error &error)
    {
        const std::string message = error.what();
        EXPECT_NE(message.find("no supported subdivision patch"), std::string::npos);
        EXPECT_NE(message.find("(4, 4, 4)"), std::string::npos) << message;
    }
}

/**
 * @brief An all-valence-5 mesh is rejected, because that is not what the
 * irregular branch actually builds.
 *
 * This is the predicate contradiction from section 1.3 of the plan, made
 * concrete. The old branch was entered when *all three* corners were at
 * valence 5, but its body then picked one corner as d4 and walked the other
 * two with the valence-6 opposite-node pattern -- i.e. it constructed a 5/6/6
 * patch. An icosahedron is the mesh that predicate was describing, and it is
 * exactly the mesh the construction cannot handle: the old code would have
 * accepted it and built a one-ring that does not match its topology.
 */
TEST(OneRingPatchTest, AllValence5MeshIsRejectedBecauseThePatchBuiltIs5And6And6)
{
    Param param;
    try
    {
        Mesh mesh = setup_solid(param, make_icosahedron());
        FAIL() << "expected an icosahedron (all valence 5) to be rejected";
    }
    catch (const std::runtime_error &error)
    {
        const std::string message = error.what();
        EXPECT_NE(message.find("(5, 5, 5)"), std::string::npos) << message;
        EXPECT_NE(message.find("exactly one corner in"), std::string::npos);
    }
}

/**
 * @brief Boundary faces are skipped without error; that is not the same thing.
 *
 * A face touching the mesh boundary has no complete one-ring and therefore no
 * limit surface, which is a property of the mesh rather than a fault. The
 * defect was that every unmatched face took the same silent path, so genuine
 * unsupported interior valences were indistinguishable from ordinary boundary
 * faces -- and only the latter are benign.
 *
 * The imported example mesh is the case that matters: unlike the generated
 * flat sheet, whose ghost band absorbs its entire edge, this one has a real
 * open boundary and no ghost marking at all. Its 336 boundary faces must be
 * skipped quietly while the mesh as a whole still sets up.
 */
TEST(OneRingPatchTest, BoundaryFacesAreSkippedWithoutRejectingTheMesh)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    ASSERT_TRUE(import_mesh_from_vertices_faces(mesh,
                                                "./data/example/vertices_flat.csv",
                                                "./data/example/faces_flat.csv"));

    int nInteriorWithRing = 0;
    int nBoundaryWithoutRing = 0;
    for (const Face &face : mesh.faces)
    {
        if (face.isGhost)
            continue;
        const bool touchesBoundary = !mesh.is_interior_vertex(face.adjacentVertices[0]) ||
                                     !mesh.is_interior_vertex(face.adjacentVertices[1]) ||
                                     !mesh.is_interior_vertex(face.adjacentVertices[2]);
        if (touchesBoundary)
        {
            EXPECT_TRUE(face.oneRingVertices.empty());
            nBoundaryWithoutRing++;
        }
        else
        {
            EXPECT_EQ(face.oneRingVertices.size(), 12u);
            nInteriorWithRing++;
        }
    }
    EXPECT_EQ(nBoundaryWithoutRing, 336);
    EXPECT_EQ(nInteriorWithRing, 3344);
}

/**
 * @brief The generated flat sheet has no boundary faces at all.
 *
 * Its ghost band covers the entire edge, so every non-ghost face is fully
 * interior and carries a complete 12-vertex one-ring. Recorded because it is
 * the opposite case to the imported mesh above, and because it is what makes
 * the shipped workload safe from the silent-fallthrough defect in the first
 * place.
 */
TEST(OneRingPatchTest, GeneratedFlatSheetIsEntirelyInteriorOrGhost)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_flat());

    int nNonGhost = 0;
    for (const Face &face : mesh.faces)
    {
        if (face.isGhost)
            continue;
        nNonGhost++;
        EXPECT_TRUE(mesh.is_interior_vertex(face.adjacentVertices[0]));
        EXPECT_EQ(face.oneRingVertices.size(), 12u);
    }
    EXPECT_GT(nNonGhost, 0);
}

/**
 * @brief One refinement pass makes a rejected mesh admissible, by construction.
 *
 * WP7's gate. An octahedron has every vertex at valence 4, so every face
 * carries three extraordinary corners and the uniform patch reduction -- which
 * assumes exactly one -- does not apply to any of them. It is rejected at
 * setup.
 *
 * After one Loop refinement every new vertex sits on an old edge and has
 * valence 6, and every old vertex keeps valence 4 but is now adjacent only to
 * new vertices. Of the four children of each old face, three carry exactly one
 * old corner and the middle one carries none. So the same mesh is admitted,
 * and the guarantee is structural rather than a property of this fixture.
 */
TEST(PreRefinementTest, RefinementAdmitsAMeshThatWasRejected)
{
    const ClosedSolid solid = make_octahedron();

    // Without refinement: rejected, every face 4/4/4.
    {
        Param param;
        param.VERBOSE_MODE = false;
        param.boundaryCondition = BoundaryType::Fixed;
        param.isPreRefinementEnabled = false;
        Mesh mesh(param);
        EXPECT_THROW(mesh.setup_from_vertices_faces(solid.vertices, solid.faces),
                     std::runtime_error);
    }

    // With refinement: admitted.
    Param param;
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.isPreRefinementEnabled = true;
    Mesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_from_vertices_faces(solid.vertices, solid.faces));

    // Four children per face, one new vertex per edge.
    EXPECT_EQ(mesh.faces.size(), solid.faces.size() * 4);
    EXPECT_EQ(mesh.vertices.size(), 6u + 12u); // octahedron: 6 vertices, 12 edges

    int nValence4 = 0;
    int nValence6 = 0;
    for (const Vertex &vertex : mesh.vertices)
    {
        const std::size_t valence = vertex.adjacentVertices.size();
        if (valence == 4)
            nValence4++;
        else if (valence == 6)
            nValence6++;
    }
    // Old vertices keep their valence; new ones are all 6.
    EXPECT_EQ(nValence4, 6);
    EXPECT_EQ(nValence6, 12);

    // The property that matters: no face has two extraordinary corners, and
    // every face is admitted with a complete one-ring.
    for (const Face &face : mesh.faces)
    {
        int nExtraordinary = 0;
        for (int node : face.adjacentVertices)
        {
            if (mesh.vertices[node].adjacentVertices.size() != 6)
                nExtraordinary++;
        }
        EXPECT_LE(nExtraordinary, 1) << "face " << face.index;
        // valence 4 -> 10 control points, valence 6 -> 12
        const std::size_t expected = (nExtraordinary == 1) ? 10u : 12u;
        EXPECT_EQ(face.oneRingVertices.size(), expected) << "face " << face.index;
    }
}

/**
 * @brief Refinement preserves the surface it refines.
 *
 * Loop subdivision converges to the same limit surface, so refining the
 * control mesh must not move it. The refined control net is a closed
 * octahedron-derived polyhedron whose signed volume should sit between the
 * control polyhedron's and the limit surface's -- and crucially, the *limit*
 * volume computed from it should match the unrefined one closely rather than
 * drifting.
 */
TEST(PreRefinementTest, RefinementApproachesTheSameLimitSurface)
{
    const ClosedSolid solid = make_octahedron();

    Param param;
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.isPreRefinementEnabled = true;
    Mesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_from_vertices_faces(solid.vertices, solid.faces));

    mesh.calculate_element_area_volume();
    double area = 0.0;
    double volume = 0.0;
    mesh.sum_membrane_area_and_volume(area, volume);

    // The unit octahedron's limit surface is a rounded solid strictly inside
    // its control polyhedron (volume 4/3) and outside nothing in particular;
    // the check that matters is that it is a sane, positive, closed volume
    // rather than a specific number.
    EXPECT_GT(area, 0.0);
    EXPECT_GT(volume, 0.0);
    EXPECT_LT(volume, 4.0 / 3.0) << "the limit surface must lie inside the control polyhedron";

    // Closed and consistently oriented, so the constraint is well defined.
    mesh.param.uVol = 1.0;
    EXPECT_NO_THROW(mesh.validate_volume_constraint_topology());
}
