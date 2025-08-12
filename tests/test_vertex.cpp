#include "test_vertex.hpp"

/**
 * @brief Construct a new TEST object; unittest Vertex::aproximate_unit_normal_vector()
 * 
 */
/*
TEST(VertexNormalVectorTest, ApproximateUnitNormalVectorTest){
    Param param;
    param.VERBOSE_MODE = false; // mute output
    Mesh mesh(param);

    mesh.vertices[0].approximate_unit_normal_vector(mesh.faces);
    std::cout << mesh.vertices[0].normVector << std::endl;

    // Expect not return nan
    EXPECT_TRUE(!std::isnan(mesh.vertices[0].normVector(0, 0)));
    EXPECT_TRUE(!std::isnan(mesh.vertices[1].normVector(1, 0)));
    EXPECT_TRUE(!std::isnan(mesh.vertices[2].normVector(2, 0)));

    mesh.vertices[0].approximate_unit_normal_vector(mesh.faces);
    std::cout << mesh.vertices[57].normVector << std::endl;

    // Expect not return nan
    EXPECT_TRUE(!std::isnan(mesh.vertices[57].normVector(0, 0)));
    EXPECT_TRUE(!std::isnan(mesh.vertices[57].normVector(1, 0)));
    EXPECT_TRUE(!std::isnan(mesh.vertices[57].normVector(2, 0)));
}
*/