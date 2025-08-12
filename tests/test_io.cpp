#include "test_io.hpp"
/**
 * @brief Test import_kv_string function
 *
 */
TEST(PopSpaceTest, BasicTest)
{
    EXPECT_EQ(pop_space("Hello World"), "HelloWorld");
    EXPECT_EQ(pop_space("Spaces   Removed"), "SpacesRemoved");
    EXPECT_EQ(pop_space("\tTabs\tRemoved"), "TabsRemoved");
    EXPECT_EQ(pop_space("Mixed\t\tSpaces \tRemoved"), "MixedSpacesRemoved");
}

/**
 * @brief Test import_kv_string function
 *
 */
TEST(ImportKVStringTest, BoundaryTypeTest)
{
    Param param;
    import_kv_string("boundaryType", "Periodic", param);
    EXPECT_EQ(param.boundaryCondition, BoundaryType::Periodic);

    import_kv_string("boundaryType", "Free", param);
    EXPECT_EQ(param.boundaryCondition, BoundaryType::Free);

    import_kv_string("boundaryType", "Unknown", param);
    EXPECT_EQ(param.boundaryCondition, BoundaryType::Fixed);
}

/**
 * @brief Test the import_param_file function.
 *
 * Using the following physical constants:
 * c0Insertion = 0.8				# spontaneous curvature of insertion
 *   c0Membrane = 0.0				# spontaneous curvature of membrane
 *   kcMembraneBending = 83.4		# membrane bending constant (pN.nm)
 *   usMembraneStretching = 250.0	# membrane streching modulus (pN/nm)
 *   uvVolumeConstraint = 0.0		# volume constraint coefficient (pN/nm^2)
 *   isGlobalConstraint = false      # true to enable global mode for area and volume
 *   KBT = 4.17                      # 1KbT = 4.17 pN.nm
 *
 *    # dynamics model parameters
 *    timeStep = 0.001                  # in us
 *    diffConst = 1.0                # diffusion constant, nm^2/us
 * This test case checks whether the import_param_file function correctly reads and updates the Param object.
 */
TEST(ParamImportTest, ImportParamFile)
{
    // Create a Param object for testing
    Param testParam;

    // Specify the file path for the parameter file
    std::string filepath = "./data/example/example.params";

    // Call the import_param_file function
    bool result = import_param_file(testParam, filepath);

    // Check if the import was successful
    EXPECT_TRUE(result);

    // Check specific values in the Param object after import
    // Modify these expectations based on the actual structure of your Param object
    EXPECT_EQ(testParam.maxIterations, 10000); // Example expectation
    EXPECT_DOUBLE_EQ(testParam.insertCurv, 0.8);
    EXPECT_DOUBLE_EQ(testParam.spontCurv, 0.0);
    EXPECT_DOUBLE_EQ(testParam.kCurv, 83.4);
    EXPECT_DOUBLE_EQ(testParam.uSurf, 250.0);
    EXPECT_DOUBLE_EQ(testParam.uVol, 0.0);
    EXPECT_FALSE(testParam.isGlobalConstraint);
    EXPECT_DOUBLE_EQ(testParam.KBT, 4.17);
    EXPECT_DOUBLE_EQ(testParam.timeStep, 0.001);
    EXPECT_DOUBLE_EQ(testParam.diffConst, 1.0);

    // Add more expectations based on your Param object - GPT
}

// Test import_mesh_from_vertices_faces function
TEST(ImportMeshTest, BasicTest)
{
    Param param;
    Mesh mesh(param);
    std::string verticesFilepath = "./data/example/vertices_flat.csv";
    std::string facesFilepath = "./data/example/faces_flat.csv";

    EXPECT_TRUE(import_mesh_from_vertices_faces(mesh, verticesFilepath, facesFilepath));

    // Add more assertions or checks based on your specific expectations
}
