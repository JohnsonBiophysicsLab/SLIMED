/**
 * @file test_mesh_force_substeps.cpp
 * @brief Gates for the mesh-quality force substeps (Param::meshForceSubsteps).
 */
#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "Dynamics.hpp"
#include "mesh/Mesh.hpp"
#include "Parameters.hpp"

namespace
{
/// A flat periodic sheet with the three wall terms on. With zeroWalls their
/// constants are zero, so the wall force is exactly zero everywhere.
void configure(Param &param, int substeps, const std::string &noise, bool zeroWalls = false)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Periodic;
    param.sideX = 100.0;
    param.sideY = 100.0;
    param.lFace = 5.0;
    param.uVol = 0.0;
    param.kCurv = 83.4;
    param.KBT = 4.17;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
    param.inPlaneDynamicsEnabled = true;
    param.edgeSpringEnabled = true;
    param.triangleShapeEnabled = true;
    param.creaseWallEnabled = true;
    param.timeStep = 0.004;
    param.diffConst = 1.0;
    param.meshForceSubsteps = substeps;
    param.meshForceSubstepNoise = noise;
    if (zeroWalls)
    {
        param.edgeSpringConstant = 0.0;
        param.triangleShapeConstant = 0.0;
        param.creaseWallConstant = 0.0;
    }
}

/// One full step (the driver's steps 1 to 4a) from the flat sheet; returns the
/// surface-point displacement of every integrated vertex, x y z interleaved.
std::vector<double> one_step(int substeps, const std::string &noise, bool zeroWalls = false)
{
    Param param;
    configure(param, substeps, noise, zeroWalls);
    DynamicMesh mesh(param);
    mesh.setup_flat();
    mesh.update_vertices_mat_with_vector();
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(param.area0, param.vol0);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.Compute_Energy_And_Force();

    Record record(4);
    DynamicModel model(mesh, record);
    model.iteration = 0;
    mesh.apply_mesh_to_surface();
    const Matrix before = mesh.matSurface;
    model.next_step();
    mesh.apply_surface_to_mesh();
    mesh.postprocess_boundary();
    model.mesh_force_substeps();
    mesh.apply_mesh_to_surface();

    std::vector<double> moves;
    for (int v = 0; v < static_cast<int>(mesh.vertices.size()); v++)
    {
        const bool slaved = !mesh.isSlavedPeriodic.empty() && mesh.isSlavedPeriodic[v] != 0;
        if (mesh.vertices[v].isGhost || slaved || mesh.vertices[v].is_fixed())
        {
            continue;
        }
        for (int axis = 0; axis < 3; axis++)
        {
            moves.push_back(mesh.matSurface(v, axis) - before(v, axis));
        }
    }
    return moves;
}
} // namespace

TEST(MeshForceSubstepsTest, RefusesAConfigurationItCannotSplit)
{
    Param param;
    configure(param, 4, "split");
    param.fdtConsistentSurfaceUpdate = false; // the substeps move the surface by M^-T F
    DynamicMesh mesh(param);
    mesh.setup_flat();
    mesh.update_vertices_mat_with_vector();
    Record record(4);
    EXPECT_THROW(DynamicModel(mesh, record), std::runtime_error);
}

/**
 * With no wall force the substeps have nothing to integrate: with the noise on
 * the first substep the step must be the shipped step exactly.
 */
TEST(MeshForceSubstepsTest, WithNoWallForceTheOnceStepIsTheShippedStep)
{
    const std::vector<double> shipped = one_step(1, "split", true);
    const std::vector<double> substepped = one_step(4, "once", true);
    ASSERT_EQ(shipped.size(), substepped.size());
    ASSERT_FALSE(shipped.empty());
    for (std::size_t k = 0; k < shipped.size(); k++)
    {
        EXPECT_EQ(shipped[k], substepped[k]) << "component " << k;
    }
}

/**
 * Split, the noise is n draws of variance 2 D dt / n: other numbers than the
 * shipped step, the same variance per step. In-plane, where the mean-height
 * pin does not act, the sample variance of a few hundred components must sit
 * on 2 D dt.
 */
TEST(MeshForceSubstepsTest, SplitNoiseKeepsTheStepVariance)
{
    const std::vector<double> shipped = one_step(1, "split");
    const std::vector<double> split = one_step(4, "split");
    ASSERT_EQ(shipped.size(), split.size());

    const double expected = 2.0 * 1.0 * 0.004;
    double sumShipped = 0.0;
    double sumSplit = 0.0;
    int count = 0;
    bool differs = false;
    for (std::size_t k = 0; k < split.size(); k++)
    {
        if (k % 3 == 2)
        {
            continue; // z carries the pin
        }
        sumShipped += shipped[k] * shipped[k];
        sumSplit += split[k] * split[k];
        differs = differs || (shipped[k] != split[k]);
        count++;
    }
    ASSERT_GT(count, 300);
    EXPECT_TRUE(differs);
    // ~400 components: a 5-sigma band on a chi-square variance is about +-36 %.
    EXPECT_NEAR(sumShipped / count / expected, 1.0, 0.4);
    EXPECT_NEAR(sumSplit / count / expected, 1.0, 0.4);
}

/**
 * With no wall force and split noise the substeps still have the noise to
 * draw: the step differs from the shipped one and keeps its variance.
 */
TEST(MeshForceSubstepsTest, WithNoWallForceTheSplitStepStillDrawsEveryShare)
{
    const std::vector<double> shipped = one_step(1, "split", true);
    const std::vector<double> split = one_step(4, "split", true);
    ASSERT_EQ(shipped.size(), split.size());
    const double expected = 2.0 * 1.0 * 0.004;
    double sum = 0.0;
    int count = 0;
    bool differs = false;
    for (std::size_t k = 0; k < split.size(); k++)
    {
        if (k % 3 == 2)
        {
            continue; // z carries the pin
        }
        sum += split[k] * split[k];
        differs = differs || (shipped[k] != split[k]);
        count++;
    }
    ASSERT_GT(count, 300);
    EXPECT_TRUE(differs);
    EXPECT_NEAR(sum / count / expected, 1.0, 0.4);
}

/// At n = 1 the substep routine is a no-op: the shipped step is the whole step.
TEST(MeshForceSubstepsTest, OneSubstepLeavesTheSurfaceAlone)
{
    Param param;
    configure(param, 1, "split");
    DynamicMesh mesh(param);
    mesh.setup_flat();
    mesh.update_vertices_mat_with_vector();
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(param.area0, param.vol0);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.Compute_Energy_And_Force();
    Record record(4);
    DynamicModel model(mesh, record);
    mesh.apply_mesh_to_surface();
    const Matrix before = mesh.matSurface;
    model.mesh_force_substeps();
    for (int v = 0; v < static_cast<int>(mesh.vertices.size()); v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            EXPECT_EQ(mesh.matSurface(v, axis), before(v, axis));
        }
    }
}
