#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "Parameters.hpp"
#include "io/io.hpp"
#include "mesh/Mesh.hpp"
#include "model/Model.hpp"
#include "model/Record.hpp"

/**
 * @file test_minimizer_gradient.cpp
 * @brief What run_flat's line search needs from the energy model: E is
 * continuous along a search line and -F . d is its slope.
 *
 * With the line search fixed, run_flat stops wherever no step along the force
 * lowers E, and on the shipped bud every remaining stop came from the model.
 * Besides the regularizer's case switch (test_regularization_blend.cpp):
 *
 * - legacy BoundaryType::Periodic: the minimizer moves the ghost rows with
 *   their sources, but their forces were dropped, so -F . d was not the slope
 *   of E. Mesh::foldLegacyPeriodicGhostForces adds them to their sources;
 * - the scaffold's harmonic bond used a signed distance, -|r| below the
 *   point's plane, so E jumped by 2 k l |r| where a vertex beside its point
 *   crossed that plane. Param::harmonicBondSideWidth makes it continuous.
 */

namespace
{

/// A 60 nm legacy periodic sheet, every term on, the regularization blended.
void configure_periodic_sheet(Param &param)
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
    param.usingRpi = true;
    param.regularizationBlendWidth = 0.1;
}

/**
 * Flat sheet as its own reference, then every real vertex jittered, the
 * ghosts moved by the displacement of the vertex they follow (as the
 * minimizer would), and that taken as the start point.
 */
void build_jittered_periodic_sheet(Mesh &mesh, Model &model)
{
    mesh.setup_flat();
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();

    std::mt19937 rng(5u);
    std::uniform_real_distribution<double> jitter(-0.8, 0.8);
    model.ncgDirection0.assign(mesh.vertices.size(), Force());
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        if (!mesh.vertices[i].isGhost)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                model.ncgDirection0[i].forceTotal.set(axis, 0, jitter(rng));
            }
        }
    }
    model.stepSize = 1.0;
    model.update_vertex_using_NCG();
    mesh.update_previous_coord_for_vertex();
    mesh.Compute_Energy_And_Force();
}

/// E at coordPrev + a d, with the boundary condition enforced.
double energy_at(Model &model, double a)
{
    model.stepSize = a;
    model.update_vertex_using_NCG();
    model.mesh.Compute_Energy_And_Force();
    return model.mesh.param.energy.energyTotal;
}

/// Relative error of -F . d against the central difference of E along d.
double slope_error(bool fold)
{
    Param param;
    configure_periodic_sheet(param);
    Mesh mesh(param);
    mesh.foldLegacyPeriodicGhostForces = fold;
    Record record(2);
    Model model(mesh, record);
    build_jittered_periodic_sheet(mesh, model);

    // Another random direction on the real vertices: the ghosts follow.
    std::mt19937 rng(17u);
    std::uniform_real_distribution<double> component(-1.0, 1.0);
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            model.ncgDirection0[i].forceTotal.set(axis, 0, mesh.vertices[i].isGhost ? 0.0 : component(rng));
        }
    }
    mesh.Compute_Energy_And_Force();
    double slope = 0.0;
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            slope -= mesh.vertices[i].force.forceTotal(axis, 0) * model.ncgDirection0[i].forceTotal(axis, 0);
        }
    }
    const double h = 1e-5;
    const double numeric = (energy_at(model, h) - energy_at(model, -h)) / (2.0 * h);
    return std::abs(slope - numeric) / std::abs(numeric);
}

/// A sheet with one scaffold point bonded to its centre vertex.
struct BondedVertex
{
    Param param;
    Mesh mesh;
    int vertex = 0;

    explicit BondedVertex(double sideWidth) : param(make_param(sideWidth)), mesh(param)
    {
        mesh.setup_flat();
        vertex = (param.nFaceX + 1) * (param.nFaceY / 2) + param.nFaceX / 2;
        Matrix point(3, 1);
        point.set(0, 0, mesh.vertices[vertex].coord(0, 0));
        point.set(1, 0, mesh.vertices[vertex].coord(1, 0));
        point.set(2, 0, mesh.vertices[vertex].coord(2, 0));
        mesh.param.scaffoldingPoints.assign(1, point);
        mesh.param.scaffoldingPoints_correspondingVertexIndex.assign(1, vertex);
    }

    static Param make_param(double sideWidth)
    {
        Param param;
        param.VERBOSE_MODE = false;
        param.boundaryCondition = BoundaryType::Fixed;
        param.sideX = 30.0;
        param.sideY = 30.0;
        param.lFace = 5.0;
        param.isEnergyHarmonicBondIncluded = true;
        param.springConst = 50.0;
        param.lbond = 9.0;
        param.harmonicBondSideWidth = sideWidth;
        return param;
    }

    /// Put the vertex at point + offset; the bond energy there.
    double energy_at(const double offset[3])
    {
        for (int axis = 0; axis < 3; axis++)
        {
            mesh.vertices[vertex].coord.set(axis, 0, mesh.param.scaffoldingPoints[0](axis, 0) + offset[axis]);
            mesh.vertices[vertex].force.forceTotal.set(axis, 0, 0.0);
        }
        return mesh.calculate_scaffolding_energy_force(false);
    }
};

} // namespace

/// Folded, -F . d is the slope of E on the legacy periodic sheet; dropped, it is not.
TEST(MinimizerGradientTest, FoldingTheLegacyPeriodicGhostForcesGivesTheSlope)
{
    EXPECT_LT(slope_error(true), 1e-5);
    EXPECT_GT(slope_error(false), 1e-2);
}

/// The ghost map is the one the minimizer's boundary enforcement uses.
TEST(MinimizerGradientTest, TheLegacyPeriodicSourceIsTheVertexAGhostFollows)
{
    Param param;
    configure_periodic_sheet(param);
    Mesh mesh(param);
    Record record(2);
    Model model(mesh, record);
    build_jittered_periodic_sheet(mesh, model);
    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        const int source = mesh.legacy_periodic_source(static_cast<int>(i));
        ASSERT_FALSE(mesh.vertices[source].isGhost) << "vertex " << i;
        if (!mesh.vertices[i].isGhost)
        {
            ASSERT_EQ(source, static_cast<int>(i));
        }
        // The jitter moved every ghost by its source's displacement.
        for (int axis = 0; axis < 3; axis++)
        {
            ASSERT_NEAR(mesh.vertices[i].coord(axis, 0) - mesh.vertices[i].coordRef(axis, 0),
                        mesh.vertices[source].coord(axis, 0) - mesh.vertices[source].coordRef(axis, 0), 1e-12)
                << "vertex " << i;
        }
    }
}

/// Beside its point, a vertex crossing the point's plane: no jump with a side width.
TEST(MinimizerGradientTest, TheBondEnergyIsContinuousAcrossThePointsPlane)
{
    for (const double width : {0.0, 1.0})
    {
        BondedVertex bonded(width);
        const double above[3] = {9.0, 0.0, 1e-9};
        const double below[3] = {9.0, 0.0, -1e-9};
        const double jump = bonded.energy_at(below) - bonded.energy_at(above);
        if (width == 0.0)
        {
            // The legacy signed distance: 2 k l |r|.
            EXPECT_NEAR(jump, 2.0 * 50.0 * 9.0 * 9.0, 1e-3);
        }
        else
        {
            EXPECT_NEAR(jump, 0.0, 1e-6);
        }
    }
}

/// Straight under its point, deeper than the width: the legacy energy less k l w.
TEST(MinimizerGradientTest, DeepUnderItsPointTheBondChargesWhatTheSignedDistanceDid)
{
    BondedVertex legacy(0.0);
    BondedVertex smooth(1.0);
    const double offset[3] = {0.0, 0.0, -4.0};
    EXPECT_NEAR(smooth.energy_at(offset), legacy.energy_at(offset) - 50.0 * 9.0 * 1.0, 1e-9);
}

/// The bond force is minus the gradient of the bond energy, in and out of the wall.
TEST(MinimizerGradientTest, TheBondForceIsMinusTheGradientOfItsEnergy)
{
    BondedVertex bonded(1.0);
    const double offsets[][3] = {{2.0, 1.0, 8.0},   // above, short of l
                                 {9.0, 3.0, 4.0},   // above, past l
                                 {8.0, -2.0, -0.4}, // in the rounded wall
                                 {3.0, 2.0, -6.0}}; // deep below
    for (const auto &offset : offsets)
    {
        bonded.energy_at(offset);
        double force[3];
        for (int axis = 0; axis < 3; axis++)
        {
            force[axis] = bonded.mesh.vertices[bonded.vertex].force.forceHarmonicBond(axis, 0);
            EXPECT_NEAR(bonded.mesh.forceOnScaffoldingPoints[0](axis, 0), -force[axis], 1e-12);
        }
        const double h = 1e-6;
        for (int axis = 0; axis < 3; axis++)
        {
            double plus[3] = {offset[0], offset[1], offset[2]};
            double minus[3] = {offset[0], offset[1], offset[2]};
            plus[axis] += h;
            minus[axis] -= h;
            const double numeric = -(bonded.energy_at(plus) - bonded.energy_at(minus)) / (2.0 * h);
            EXPECT_NEAR(force[axis], numeric, 1e-5 * std::max(1.0, std::abs(numeric)))
                << "offset " << offset[0] << "," << offset[1] << "," << offset[2] << " axis " << axis;
        }
    }
}

/// The three settings reach Param from the parameter file, and nonsense is refused.
TEST(MinimizerGradientTest, TheSettingsAreReadFromTheParameterFile)
{
    const char *tmp = std::getenv("TMPDIR");
    const std::string path = std::string(tmp ? tmp : "/tmp") + "/slimed_minimizer_gradient_keys.params";
    auto write = [&path](const std::string &body) {
        std::ofstream out(path);
        out << body;
    };
    Param param;
    EXPECT_DOUBLE_EQ(param.regularizationBlendWidth, 0.05);
    EXPECT_TRUE(param.periodicGhostForceFolding);
    EXPECT_DOUBLE_EQ(param.harmonicBondSideWidth, 1.0);

    write("regularizationBlendWidth = 0\nperiodicGhostForceFolding = false\nharmonicBondSideWidth = 0.5\n");
    EXPECT_TRUE(import_param_file(param, path));
    EXPECT_DOUBLE_EQ(param.regularizationBlendWidth, 0.0);
    EXPECT_FALSE(param.periodicGhostForceFolding);
    EXPECT_DOUBLE_EQ(param.harmonicBondSideWidth, 0.5);

    for (const char *bad : {"regularizationBlendWidth = -0.1\n", "periodicGhostForceFolding = yes\n",
                            "harmonicBondSideWidth = nan\n"})
    {
        write(bad);
        EXPECT_THROW(import_param_file(param, path), std::runtime_error) << bad;
    }
    std::remove(path.c_str());
}
