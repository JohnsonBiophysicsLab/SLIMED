#include "test_scaffolding_points.hpp"


// Test fixture for Mesh class
class ScaffoldingPointsTest : public ::testing::Test {
protected:
    void SetUp() override {
        import_param_file(param, "input.params");
        // Set up scaffolding points and radius for testing
        Matrix sp1 = mat_calloc(3, 1);
        Matrix sp2 = mat_calloc(3, 1);
        Matrix sp3 = mat_calloc(3, 1);
        sp2.set(0, 0, 4.0);
        sp2.set(1, 0, 5.0);
        sp2.set(2, 0, 6.0);
        sp3.set(0, 0, 7.0);
        sp3.set(1, 0, 8.0);
        sp3.set(2, 0, 9.0);
        param.scaffoldingPoints = { sp1, sp2, sp3 };
        param.scaffoldingSphereRaidus = 10.0;

        mesh = new Mesh(param);
        mesh->setup_flat();

        //mesh->param.scaffoldingPoints = import_scaffolding_mesh("./tests/COM.csv");

    }

    void TearDown() override {
        delete mesh;
    }

    Param param;
    Mesh* mesh;
};

// Test case for finding center of scaffolding sphere with default option
TEST_F(ScaffoldingPointsTest, FindCenterOfScaffoldingSphereDefault) {
    bool use_default = true;
    Matrix result = mesh->find_center_of_scaffolding_sphere(use_default);
    EXPECT_DOUBLE_EQ(result(0, 0), 0);
    EXPECT_DOUBLE_EQ(result(1, 0), 0);
    EXPECT_DOUBLE_EQ(result(2, 0), 0);
}

// Test case for finding center of scaffolding sphere without default option
TEST_F(ScaffoldingPointsTest, FindCenterOfScaffoldingSphere) {
    bool use_default = false;
    Matrix result = mesh->find_center_of_scaffolding_sphere(use_default);
    // Expected result depends on the test data provided in the fixture
    // You should provide the expected result based on your input data
    // For example, if (7, 8, 9) has the maximum Z coordinate, and radius is 1, the expected result would be (7, 8, 8)
    EXPECT_DOUBLE_EQ(result(0, 0), 7);
    EXPECT_DOUBLE_EQ(result(1, 0), 8);
    EXPECT_DOUBLE_EQ(result(2, 0), -1);
}

// Test case for approximate_scaffolding_cap_radius with default option
TEST_F(ScaffoldingPointsTest, ApproximateScaffoldingCapRadiusDefault) {
    bool use_default = true;
    double result = mesh->approximate_scaffolding_cap_radius(use_default);
    // Provide the expected result based on test data
    // For example, if the default center (0,0,0) is used and all points are on the XY plane, the expected radius would be the distance of the furthest point from the origin in the XY plane
    EXPECT_DOUBLE_EQ(result, sqrt(7 * 7 + 8 * 8)); // Assuming (7,8,9) is the furthest point in the XY plane
}

// Test case for approximate_scaffolding_cap_radius without default option
TEST_F(ScaffoldingPointsTest, ApproximateScaffoldingCapRadius) {
    bool use_default = false;
    double result = mesh->approximate_scaffolding_cap_radius(use_default);
    // Provide the expected result based on test data
    // You should calculate the expected result based on your test data
    // For example, if the center is determined based on the maximum Z coordinate and radius is 1, and (7,8,9) has the maximum Z coordinate, the expected radius would be sqrt(7 * 7 + 8 * 8) - 1
    EXPECT_DOUBLE_EQ(result, sqrt(7 * 7 + 8 * 8));
}
/*
// Test case for move_vertices_based_on_scaffolding with fixDir set to false
TEST_F(ScaffoldingPointsTest, MoveVerticesBasedOnScaffoldingFixDirFalse) {
    bool fixDir = false;
    bool result = mesh->move_vertices_based_on_scaffolding(fixDir);
    // Add assertions to check if vertices are moved correctly based on the implemented logic
    // For example, you can compare the coordinates of vertices after the movement
    EXPECT_FALSE(result); // Placeholder assertion as the implementation is not provided in the test
}
*/