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