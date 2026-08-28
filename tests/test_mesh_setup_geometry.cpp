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
