#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Parameters.hpp"
#include "io/io.hpp"
#include "mesh/Mesh.hpp"
#include "model/Model.hpp"
#include "model/Record.hpp"

/**
 * @file test_regularization_settings.cpp
 * @brief The regularization's R-pi settings: who sets them, and who does not.
 *
 * usingRpi, gamaShape and gamaArea are settings of the energy. They come from
 * the parameter file and nothing changes them during a run. The minimizer used
 * to carry its own copy of usingRpi and switch it off when NCG fell back to
 * the simple line search; the regularizer never read that copy, and the tests
 * here pin the behaviour that was kept when the copy was removed: the fallback
 * refreshes the reference and leaves the mode alone.
 */

namespace
{

void configure_sheet(Param &param)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Periodic;
    param.sideX = 60.0;
    param.sideY = 60.0;
    param.lFace = 5.0;
    param.kCurv = 83.4;
    param.uSurf = 250.0;
    param.uVol = 0.0;
    param.isGlobalConstraint = true;
    param.isEnergyHarmonicBondIncluded = false;
}

/// A flat sheet whose reference is the flat sheet.
void build_sheet(Mesh &mesh)
{
    mesh.setup_flat();
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
}

/// A vertex in the middle of the sheet, away from the ghost rings.
int centre_vertex(const Mesh &mesh)
{
    return (mesh.param.nFaceX + 1) * (mesh.param.nFaceY / 2) + mesh.param.nFaceX / 2;
}

/**
 * Slide the centre vertex 2 nm along x. Its six faces go from equilateral to
 * a side-length scatter of roughly 0.1: past a gamaShape of 0.05, short of the
 * default 0.2.
 */
void shear_centre_vertex(Mesh &mesh)
{
    Vertex &vertex = mesh.vertices[centre_vertex(mesh)];
    vertex.coord.set(0, 0, vertex.coord(0, 0) + 2.0);
}

double max_regularization_force(const Mesh &mesh)
{
    double maxComponent = 0.0;
    for (const Vertex &vertex : mesh.vertices)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            maxComponent = std::max(maxComponent,
                                    std::abs(vertex.force.forceRegularization(axis, 0)));
        }
    }
    return maxComponent;
}

std::string scratch_file(const std::string &name)
{
    return "./test_regularization_settings_" + name;
}

void write_param_file(const std::string &filepath, const std::string &body)
{
    std::ofstream out(filepath);
    ASSERT_TRUE(out.is_open());
    out << body;
}

} // namespace

/// The three keys reach Param, with the opposite values seeded first.
TEST(RegularizationSettingsTest, TheKeysAreReadFromTheParameterFile)
{
    const std::string filepath = scratch_file("keys.params");
    write_param_file(filepath,
                     "usingRpi = false      # comment\n"
                     "gamaShape = 0.05\n"
                     "gamaArea = 0.125\n");

    Param param;
    ASSERT_TRUE(param.usingRpi);
    ASSERT_DOUBLE_EQ(param.gamaShape, 0.2);
    ASSERT_DOUBLE_EQ(param.gamaArea, 0.2);

    EXPECT_TRUE(import_param_file(param, filepath));
    EXPECT_FALSE(param.usingRpi);
    EXPECT_DOUBLE_EQ(param.gamaShape, 0.05);
    EXPECT_DOUBLE_EQ(param.gamaArea, 0.125);

    write_param_file(filepath, "usingRpi = true\n");
    EXPECT_TRUE(import_param_file(param, filepath));
    EXPECT_TRUE(param.usingRpi);

    std::remove(filepath.c_str());
}

/// A value that is not a setting stops the run instead of picking one.
TEST(RegularizationSettingsTest, AnInvalidValueIsRefused)
{
    const std::string filepath = scratch_file("invalid.params");
    Param param;

    write_param_file(filepath, "usingRpi = True\n");
    EXPECT_THROW(import_param_file(param, filepath), std::runtime_error);
    EXPECT_TRUE(param.usingRpi);

    write_param_file(filepath, "gamaShape = -0.1\n");
    EXPECT_THROW(import_param_file(param, filepath), std::runtime_error);
    EXPECT_DOUBLE_EQ(param.gamaShape, 0.2);

    write_param_file(filepath, "gamaArea = nan\n");
    EXPECT_THROW(import_param_file(param, filepath), std::runtime_error);
    EXPECT_DOUBLE_EQ(param.gamaArea, 0.2);

    std::remove(filepath.c_str());
}

/**
 * @brief gamaShape decides which faces are adaptive, and usingRpi whether any are.
 *
 * The same sheared sheet under four settings. With the adaptive cases off, or
 * with the threshold out of reach, every face is in the reference-length case
 * and the two energies are the same number.
 */
TEST(RegularizationSettingsTest, TheSettingsReachTheRegularizer)
{
    Param param;
    configure_sheet(param);
    Mesh mesh(param);
    build_sheet(mesh);
    shear_centre_vertex(mesh);

    // Default thresholds: the shear is short of 0.2.
    mesh.Compute_Energy_And_Force();
    EXPECT_EQ(param.deformationCount.shapeDeformCount, 0);
    const double referenceLengthEnergy = param.energy.energyRegularization;
    EXPECT_GT(referenceLengthEnergy, 0.0);

    // A lower shape threshold puts the six sheared faces in the adaptive case.
    param.gamaShape = 0.05;
    mesh.Compute_Energy_And_Force();
    EXPECT_EQ(param.deformationCount.shapeDeformCount, 6);
    const double adaptiveEnergy = param.energy.energyRegularization;
    EXPECT_NE(adaptiveEnergy, referenceLengthEnergy);

    // The same threshold with the adaptive cases off.
    param.usingRpi = false;
    mesh.Compute_Energy_And_Force();
    EXPECT_EQ(param.deformationCount.shapeDeformCount, 0);
    EXPECT_EQ(param.deformationCount.areaDeformCount, 0);
    EXPECT_DOUBLE_EQ(param.energy.energyRegularization, referenceLengthEnergy);

    // The area threshold, alone: shape out of reach, area at zero.
    param.usingRpi = true;
    param.gamaShape = 100.0;
    param.gamaArea = 0.01;
    mesh.Compute_Energy_And_Force();
    EXPECT_EQ(param.deformationCount.shapeDeformCount, 0);
    EXPECT_GT(param.deformationCount.areaDeformCount, 0);
}

/**
 * @brief What usingRpi = false buys: the term vanishes at a refreshed reference.
 *
 * With the adaptive cases on, a shape-deformed face is measured against an
 * equilateral triangle, so refreshing the reference does not release it.
 */
TEST(RegularizationSettingsTest, TheTermVanishesAtARefreshedReferenceOnlyWithoutRpi)
{
    Param param;
    configure_sheet(param);
    param.gamaShape = 0.05;
    Mesh mesh(param);
    build_sheet(mesh);
    shear_centre_vertex(mesh);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();

    param.usingRpi = true;
    mesh.Compute_Energy_And_Force();
    EXPECT_EQ(param.deformationCount.shapeDeformCount, 6);
    EXPECT_GT(param.energy.energyRegularization, 1.0);
    EXPECT_GT(max_regularization_force(mesh), 1.0);

    param.usingRpi = false;
    mesh.Compute_Energy_And_Force();
    EXPECT_EQ(param.energy.energyRegularization, 0.0);
    EXPECT_EQ(max_regularization_force(mesh), 0.0);
}

/**
 * @brief The NCG fallback refreshes the reference and leaves usingRpi alone.
 *
 * c2 = 0 makes the curvature condition unsatisfiable, so the search runs out
 * of step and takes the fallback whatever else the NCG branch does. Run with
 * the mode configured each way: the minimizer must not turn it off (the old
 * dead toggle's intent) and must not turn it on (what resetting it before
 * every search would do to a run configured without it).
 */
TEST(RegularizationSettingsTest, TheLineSearchFallbackKeepsTheConfiguredMode)
{
    for (const bool configured : {true, false})
    {
        Param param;
        configure_sheet(param);
        param.usingRpi = configured;
        param.gamaShape = 0.05;
        Mesh mesh(param);
        build_sheet(mesh);

        // Off the reference, so the fallback's refresh is observable, and out
        // of plane as well, so that bending drives the search too.
        shear_centre_vertex(mesh);
        Vertex &centre = mesh.vertices[centre_vertex(mesh)];
        centre.coord.set(2, 0, centre.coord(2, 0) + 1.0);
        mesh.Compute_Energy_And_Force();
        mesh.update_previous_coord_for_vertex();
        mesh.update_previous_force_for_vertex();
        ASSERT_NE(centre.coordRef(0, 0), centre.coordPrev(0, 0));

        Record record(4);
        record.add(param.area, param.energy, mesh.calculate_mean_force());
        Model model(mesh, record);
        model.oa.usingNCG = true;
        model.oa.c2 = 0.0;
        model.determine_trial_step_size();
        ASSERT_GT(model.oa.trialStepSize, 0.0);

        model.linear_search_for_stepsize_to_minimize_energy();

        EXPECT_TRUE(model.oa.isNCGstuck) << "the fallback was not taken";
        EXPECT_EQ(param.usingRpi, configured);
        for (const Vertex &vertex : mesh.vertices)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                ASSERT_EQ(vertex.coordRef(axis, 0), vertex.coordPrev(axis, 0));
            }
        }
    }
}

/**
 * @brief A restart takes usingRpi from its parameter file, not from the checkpoint.
 *
 * The checkpoint's "oa" line keeps the slot the minimizer's copy of the flag
 * was written to. Checkpoints written before the copy was removed hold a 0
 * there (no line search had run yet) or a 1; both must load, and neither may
 * override the setting.
 */
TEST(RegularizationSettingsTest, ARestartKeepsTheParameterFilesMode)
{
    const std::string path = scratch_file("restart.chk");
    const std::string editedPath = scratch_file("restart_edited.chk");

    Param sourceParam;
    configure_sheet(sourceParam);
    Mesh source(sourceParam);
    build_sheet(source);
    source.Compute_Energy_And_Force();
    Record sourceRecord(4);
    sourceRecord.add(sourceParam.area, sourceParam.energy, source.calculate_mean_force());
    Model sourceModel(source, sourceRecord);
    sourceModel.stepSize = 0.125;
    ASSERT_TRUE(write_model_restart_checkpoint(sourceModel, path, 5));

    for (const char slot : {'0', '1'})
    {
        // Rewrite the third field of the "oa" line.
        {
            std::ifstream in(path);
            std::ofstream out(editedPath);
            std::string line;
            bool edited = false;
            while (std::getline(in, line))
            {
                if (line.rfind("oa ", 0) == 0)
                {
                    std::istringstream fields(line);
                    std::string tag, usingNCG, isNCGstuck, reserved;
                    fields >> tag >> usingNCG >> isNCGstuck >> reserved;
                    std::string rest;
                    std::getline(fields, rest);
                    line = tag + ' ' + usingNCG + ' ' + isNCGstuck + ' ' + slot + rest;
                    edited = true;
                }
                out << line << '\n';
            }
            ASSERT_TRUE(edited);
        }

        for (const bool configured : {true, false})
        {
            Param targetParam;
            configure_sheet(targetParam);
            targetParam.usingRpi = configured;
            Mesh target(targetParam);
            build_sheet(target);
            Record targetRecord(4);
            Model targetModel(target, targetRecord);

            ASSERT_TRUE(load_model_restart_checkpoint(targetModel, editedPath));
            EXPECT_EQ(targetModel.iteration, 5);
            EXPECT_DOUBLE_EQ(targetModel.stepSize, 0.125);
            EXPECT_EQ(targetParam.usingRpi, configured);
        }
    }

    std::remove(path.c_str());
    std::remove(editedPath.c_str());
}
