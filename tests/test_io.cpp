#include "test_io.hpp"

#include <cstdio>

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
    // A CRLF checkout leaves a trailing carriage return on every value.
    EXPECT_EQ(pop_space("Periodic\r"), "Periodic");
    EXPECT_EQ(pop_space("cpu\r"), "cpu");
    EXPECT_EQ(pop_space(" true \t\r"), "true");
}

TEST(TrimTrailingSemicolonTest, BasicTest)
{
    EXPECT_EQ(trim_trailing_semicolon("slimed_restart.chk;"), "slimed_restart.chk");
    EXPECT_EQ(trim_trailing_semicolon("false;"), "false");
    EXPECT_EQ(trim_trailing_semicolon("slimed_restart.chk"), "slimed_restart.chk");
    EXPECT_EQ(trim_trailing_semicolon(""), "");
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

    // pop_space() strips the carriage return a CRLF checkout leaves behind, so
    // "Periodic\r" must not fall through to the Fixed branch.
    import_kv_string("boundaryType", pop_space("Periodic\r"), param);
    EXPECT_EQ(param.boundaryCondition, BoundaryType::Periodic);
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
    EXPECT_FALSE(testParam.thermalFluctuationEnabled);
    EXPECT_EQ(testParam.thermalFluctuationInterval, 50);
    EXPECT_DOUBLE_EQ(testParam.thermalFluctuationTemperatureKelvin, 298.0);
    EXPECT_DOUBLE_EQ(testParam.thermalFluctuationMinTemperatureKelvin, 298.0);
    EXPECT_DOUBLE_EQ(testParam.thermalFluctuationCoolingRate, 1.0);
    EXPECT_DOUBLE_EQ(testParam.thermalFluctuationStepScale, 0.02);
    EXPECT_DOUBLE_EQ(testParam.timeStep, 0.001);
    EXPECT_DOUBLE_EQ(testParam.diffConst, 1.0);

    // Add more expectations based on your Param object - GPT
}

/**
 * @brief A CRLF-terminated parameter file must parse exactly like an LF one.
 *
 * On a Windows checkout (core.autocrlf=true) every value used to keep a trailing
 * '\r'. Numeric values survived it because std::stod ignores trailing junk, but
 * exact comparisons did not: forceBackend threw "must be cpu, gpu or auto; got
 * 'cpu'", and boundaryType silently fell through to BoundaryType::Fixed. The file
 * is written here rather than committed as a fixture, because .gitattributes now
 * forces every checked-out file to LF.
 */
TEST(ParamImportTest, CrlfParamFileParsesLikeLf)
{
    const std::string filepath = "./test_crlf_line_endings.params";

    // binary mode, so the CRLF endings reach the file verbatim on every platform
    std::ofstream out(filepath, std::ios::out | std::ios::binary);
    ASSERT_TRUE(out.is_open());
    out << "# a parameter file as checked out on Windows\r\n";
    out << "boundaryType = Periodic\r\n";
    out << "forceBackend = cpu\r\n";
    out << "maxIterations = 500000\r\n";
    out << "timeStep = 0.001;\r\n";
    out.close();

    Param param;
    // seed the opposite value so a parse that silently does nothing cannot pass
    param.boundaryCondition = BoundaryType::Fixed;

    EXPECT_TRUE(import_param_file(param, filepath));

    EXPECT_EQ(param.boundaryCondition, BoundaryType::Periodic);
    EXPECT_EQ(param.forceBackend, "cpu");
    EXPECT_EQ(param.maxIterations, 500000);
    EXPECT_DOUBLE_EQ(param.timeStep, 0.001);

    std::remove(filepath.c_str());
}

TEST(ImportKVStringTest, ThermalFluctuationParameters)
{
    Param param;
    EXPECT_TRUE(import_kv_string("thermalFluctuationEnabled", "true", param));
    EXPECT_TRUE(import_kv_string("thermalFluctuationInterval", "7", param));
    EXPECT_TRUE(import_kv_string("thermalFluctuationTemperatureKelvin", "310.0", param));
    EXPECT_TRUE(import_kv_string("thermalFluctuationMinTemperatureKelvin", "298.0", param));
    EXPECT_TRUE(import_kv_string("thermalFluctuationCoolingRate", "0.95", param));
    EXPECT_TRUE(import_kv_string("thermalFluctuationStepScale", "0.01", param));

    EXPECT_TRUE(param.thermalFluctuationEnabled);
    EXPECT_EQ(param.thermalFluctuationInterval, 7);
    EXPECT_DOUBLE_EQ(param.thermalFluctuationTemperatureKelvin, 310.0);
    EXPECT_DOUBLE_EQ(param.thermalFluctuationMinTemperatureKelvin, 298.0);
    EXPECT_DOUBLE_EQ(param.thermalFluctuationCoolingRate, 0.95);
    EXPECT_DOUBLE_EQ(param.thermalFluctuationStepScale, 0.01);
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
