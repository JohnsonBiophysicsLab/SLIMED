#include "test_external_force.hpp"

/**
 * The external force pushes a region of the membrane with a constant force
 * per vertex between two iterations. Four things can go wrong, and each has
 * a test: the keys might not parse; the region might select the wrong
 * vertices; the force might act at the wrong time; and it might not be the
 * force it says. A run with the force off must be the run it always was,
 * which the byte-identity gate on the shipped periodic workload covers.
 *
 * @see docs/external_force.md
 */

namespace
{

/// A clamped 60 nm sheet under the per-vertex boundary: no images, so every
/// vertex is either free or held, and the region is easy to reason about.
void configure_clamped_sheet(Param &param)
{
    param.VERBOSE_MODE = false;
    param.sideX = 60.0;
    param.sideY = 60.0;
    param.lFace = 5.0;
    param.kCurv = 83.4;
    param.KBT = 4.17;
    param.uSurf = 250.0;
    param.uVol = 0.0;
    param.isGlobalConstraint = true;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
    param.randomSeed = 20260923u;
    param.timeStep = 1.0e-3;
    param.diffConst = 1.0;
    param.surfaceSolver = "iterative";
    param.boundaryCondition = BoundaryType::Mixed;
    param.boundaryConditionX = BoundaryType::Fixed;
    param.boundaryConditionY = BoundaryType::Fixed;
    param.fixedBoundaryRings = 2;
}

/// What run_dynamics_flat() does between setup and the first step.
void prepare(DynamicMesh &mesh)
{
    for (Vertex &vertex : mesh.vertices)
    {
        vertex.coord.set(2, 0, 10.0);
    }
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.param.vol0 = 0.0;
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.Compute_Energy_And_Force();
    mesh.update_previous_force_for_vertex();
    mesh.update_vertices_mat_with_vector();
    mesh.apply_mesh_to_surface();
}

/// A mesh, its record and its model, built the way the driver builds them.
struct PushRun
{
    Param param;
    std::unique_ptr<DynamicMesh> mesh;
    Record record{16};
    std::unique_ptr<DynamicModel> model;

    explicit PushRun(const Param &source) : param(source)
    {
        mesh.reset(new DynamicMesh(param));
        mesh->setup_flat();
        prepare(*mesh);
        record.add(mesh->param.area, mesh->param.energy, mesh->calculate_mean_force());
        model.reset(new DynamicModel(*mesh, record));
    }

    /**
     * @brief One Brownian step at this iteration, leaving matSurface as the
     * step produced it.
     *
     * The control net is never updated from it, so every call starts from
     * the same configuration and the same forces, and two runs with the same
     * seed draw the same noise at the same iteration: whatever differs
     * between them is the external force alone.
     */
    void step(int iteration)
    {
        mesh->apply_mesh_to_surface();
        model->iteration = iteration;
        model->next_step();
    }
};

int count_independent(const Mesh &mesh)
{
    int n = 0;
    for (int v = 0; v < static_cast<int>(mesh.vertices.size()); v++)
    {
        n += mesh.is_independent_vertex(v) ? 1 : 0;
    }
    return n;
}

} // namespace

TEST(ExternalForceTest, TheParameterKeysParse)
{
    Param param;
    EXPECT_FALSE(param.externalForceEnabled);
    EXPECT_TRUE(import_kv_string("externalForceEnabled", "true", param));
    EXPECT_TRUE(param.externalForceEnabled);
    EXPECT_TRUE(import_kv_string("externalForceX", "0.5", param));
    EXPECT_DOUBLE_EQ(param.externalForceX, 0.5);
    EXPECT_TRUE(import_kv_string("externalForceY", "-1.5", param));
    EXPECT_DOUBLE_EQ(param.externalForceY, -1.5);
    EXPECT_TRUE(import_kv_string("externalForceZ", "4.0", param));
    EXPECT_DOUBLE_EQ(param.externalForceZ, 4.0);
    EXPECT_TRUE(import_kv_string("externalForceShape", "rectangle", param));
    EXPECT_EQ(param.externalForceShape, "rectangle");
    EXPECT_TRUE(import_kv_string("externalForceCenterX", "12.5", param));
    EXPECT_DOUBLE_EQ(param.externalForceCenterX, 12.5);
    EXPECT_TRUE(import_kv_string("externalForceCenterY", "-7.5", param));
    EXPECT_DOUBLE_EQ(param.externalForceCenterY, -7.5);
    EXPECT_TRUE(import_kv_string("externalForceRadius", "15", param));
    EXPECT_DOUBLE_EQ(param.externalForceRadius, 15.0);
    EXPECT_TRUE(import_kv_string("externalForceHalfWidthX", "10", param));
    EXPECT_DOUBLE_EQ(param.externalForceHalfWidthX, 10.0);
    EXPECT_TRUE(import_kv_string("externalForceHalfWidthY", "8", param));
    EXPECT_DOUBLE_EQ(param.externalForceHalfWidthY, 8.0);
    EXPECT_TRUE(import_kv_string("externalForceStartIteration", "100", param));
    EXPECT_EQ(param.externalForceStartIteration, 100);
    EXPECT_TRUE(import_kv_string("externalForceStopIteration", "4000", param));
    EXPECT_EQ(param.externalForceStopIteration, 4000);
    EXPECT_TRUE(import_kv_string("externalForceEnabled", "false", param));
    EXPECT_FALSE(param.externalForceEnabled);
}

/**
 * @brief The region is exactly the free vertices inside the shape.
 *
 * Checked against the definition vertex by vertex, for a disc and for a
 * rectangle, and the count is what a 5 nm lattice puts inside that area.
 */
TEST(ExternalForceTest, TheRegionIsTheFreeVerticesInsideTheShape)
{
    for (const std::string shape : {"disc", "rectangle"})
    {
        Param param;
        configure_clamped_sheet(param);
        param.externalForceEnabled = true;
        param.externalForceZ = 2.0;
        param.externalForceShape = shape;
        param.externalForceCenterX = 2.5;
        param.externalForceCenterY = -2.0;
        param.externalForceRadius = 12.0;
        param.externalForceHalfWidthX = 10.0;
        param.externalForceHalfWidthY = 8.0;
        PushRun run(param);
        const DynamicModel &model = *run.model;
        const DynamicMesh &mesh = *run.mesh;
        ASSERT_EQ(model.externalForceMask.size(), mesh.vertices.size()) << shape;

        int expected = 0;
        int selected = 0;
        for (int v = 0; v < static_cast<int>(mesh.vertices.size()); v++)
        {
            const Vertex &vertex = mesh.vertices[v];
            const double dx = vertex.coord(0, 0) - param.externalForceCenterX;
            const double dy = vertex.coord(1, 0) - param.externalForceCenterY;
            const bool inside = (shape == "disc")
                                    ? (dx * dx + dy * dy <= 12.0 * 12.0)
                                    : (std::fabs(dx) <= 10.0 && std::fabs(dy) <= 8.0);
            const bool independent = mesh.is_independent_vertex(v);
            EXPECT_EQ(model.externalForceMask[v] != 0, inside && independent) << shape << " " << v;
            if (model.externalForceMask[v] != 0)
            {
                selected++;
                EXPECT_FALSE(vertex.is_fixed()) << shape << " " << v;
            }
            if (inside && independent)
            {
                expected++;
            }
        }
        EXPECT_EQ(selected, expected) << shape;
        EXPECT_EQ(model.externalForceVertexCount, expected) << shape;
        // The lattice puts one vertex per (sqrt(3)/2) lFace^2 = 21.65 nm^2:
        // about 21 in the disc (452 nm^2) and 15 in the rectangle (320 nm^2).
        const int lattice = (shape == "disc") ? 21 : 15;
        EXPECT_GE(expected, lattice - 6) << shape;
        EXPECT_LE(expected, lattice + 6) << shape;
    }
}

/**
 * @brief A region that covers the whole sheet selects every free vertex and
 * no clamped one: the clamps are not coordinates and are never pushed.
 */
TEST(ExternalForceTest, TheRegionNeverReachesAClampedVertex)
{
    Param param;
    configure_clamped_sheet(param);
    param.externalForceEnabled = true;
    param.externalForceZ = 1.0;
    param.externalForceRadius = 1000.0;
    PushRun run(param);
    const DynamicModel &model = *run.model;
    const DynamicMesh &mesh = *run.mesh;
    EXPECT_EQ(model.externalForceVertexCount, count_independent(mesh));
    EXPECT_GT(model.externalForceVertexCount, 0);
    for (int v = 0; v < static_cast<int>(mesh.vertices.size()); v++)
    {
        if (mesh.vertices[v].is_fixed())
        {
            EXPECT_EQ(model.externalForceMask[v], 0) << v;
        }
    }
}

TEST(ExternalForceTest, TheForceActsOnlyInsideItsWindow)
{
    Param off;
    configure_clamped_sheet(off);
    Param on = off;
    on.externalForceEnabled = true;
    on.externalForceZ = 2.0;
    on.externalForceRadius = 12.0;
    on.externalForceStartIteration = 3;
    on.externalForceStopIteration = 5;

    PushRun reference(off);
    PushRun pushed(on);
    EXPECT_FALSE(reference.model->external_force_active(3));
    EXPECT_FALSE(pushed.model->external_force_active(2));
    EXPECT_TRUE(pushed.model->external_force_active(3));
    EXPECT_TRUE(pushed.model->external_force_active(4));
    EXPECT_FALSE(pushed.model->external_force_active(5));
    EXPECT_FALSE(pushed.model->external_force_active(6));

    // Outside the window the pushed run takes the reference's step to the bit.
    for (int iteration : {2, 5, 100})
    {
        reference.step(iteration);
        pushed.step(iteration);
        for (int v = 0; v < static_cast<int>(reference.mesh->vertices.size()); v++)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                EXPECT_EQ(pushed.mesh->matSurface(v, axis), reference.mesh->matSurface(v, axis))
                    << "iteration " << iteration << " vertex " << v << " axis " << axis;
            }
        }
    }

    // Inside it, the pushed vertices move and nothing else does.
    reference.step(4);
    pushed.step(4);
    int moved = 0;
    for (int v = 0; v < static_cast<int>(reference.mesh->vertices.size()); v++)
    {
        const double dz = pushed.mesh->matSurface(v, 2) - reference.mesh->matSurface(v, 2);
        if (pushed.model->externalForceMask[v] != 0)
        {
            EXPECT_GT(dz, 0.0) << v;
            moved++;
        }
        else
        {
            EXPECT_EQ(dz, 0.0) << v;
        }
    }
    EXPECT_EQ(moved, pushed.model->externalForceVertexCount);

    // A negative stop never switches off.
    Param forever = on;
    forever.externalForceStopIteration = -1;
    PushRun open(forever);
    EXPECT_TRUE(open.model->external_force_active(3));
    EXPECT_TRUE(open.model->external_force_active(1000000));
    EXPECT_FALSE(open.model->external_force_active(2));
}

/**
 * @brief The force is the force it says: each pushed surface point moves by
 * exactly (D dt / kT) F more than it would have, and nothing else moves.
 */
TEST(ExternalForceTest, AnUpwardForceLiftsThePushedVerticesByTheBrownianMobility)
{
    Param off;
    configure_clamped_sheet(off);
    Param on = off;
    on.externalForceEnabled = true;
    on.externalForceZ = 2.0;
    on.externalForceRadius = 12.0;

    PushRun reference(off);
    PushRun pushed(on);
    reference.step(0);
    pushed.step(0);

    const double expected = pushed.model->forceScaleConst * on.externalForceZ;
    EXPECT_NEAR(expected, on.diffConst * on.timeStep / on.KBT * on.externalForceZ, 1e-15);
    EXPECT_GT(expected, 0.0);
    for (int v = 0; v < static_cast<int>(reference.mesh->vertices.size()); v++)
    {
        const double dz = pushed.mesh->matSurface(v, 2) - reference.mesh->matSurface(v, 2);
        if (pushed.model->externalForceMask[v] != 0)
        {
            EXPECT_NEAR(dz, expected, 1e-12) << v;
        }
        else
        {
            EXPECT_EQ(dz, 0.0) << v;
        }
        // In-plane motion is off, so x and y are untouched either way.
        EXPECT_EQ(pushed.mesh->matSurface(v, 0), reference.mesh->matSurface(v, 0)) << v;
        EXPECT_EQ(pushed.mesh->matSurface(v, 1), reference.mesh->matSurface(v, 1)) << v;
    }
}

/**
 * @brief With in-plane motion off there is no in-plane coordinate to push,
 * so an in-plane force does nothing at all -- like every other in-plane force.
 */
TEST(ExternalForceTest, AnInPlaneForceDoesNothingWhileInPlaneMotionIsOff)
{
    Param off;
    configure_clamped_sheet(off);
    Param on = off;
    on.externalForceEnabled = true;
    on.externalForceX = 3.0;
    on.externalForceY = -3.0;
    on.externalForceRadius = 12.0;

    PushRun reference(off);
    PushRun pushed(on);
    reference.step(0);
    pushed.step(0);
    for (int v = 0; v < static_cast<int>(reference.mesh->vertices.size()); v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            EXPECT_EQ(pushed.mesh->matSurface(v, axis), reference.mesh->matSurface(v, axis))
                << v << " " << axis;
        }
    }
}

TEST(ExternalForceTest, AnEmptyRegionOrABadShapeIsRefused)
{
    Param param;
    configure_clamped_sheet(param);
    param.externalForceEnabled = true;
    param.externalForceZ = 1.0;

    // A centre off the sheet selects nothing: a mis-set centre, not a run.
    param.externalForceCenterX = 500.0;
    param.externalForceCenterY = 500.0;
    EXPECT_THROW(PushRun run(param), std::runtime_error);
    param.externalForceCenterX = 0.0;
    param.externalForceCenterY = 0.0;

    param.externalForceShape = "ring";
    EXPECT_THROW(PushRun run(param), std::runtime_error);
    param.externalForceShape = "disc";

    param.externalForceRadius = 0.0;
    EXPECT_THROW(PushRun run(param), std::runtime_error);
    param.externalForceRadius = 12.0;

    param.externalForceShape = "rectangle";
    param.externalForceHalfWidthY = -1.0;
    EXPECT_THROW(PushRun run(param), std::runtime_error);
    param.externalForceHalfWidthY = 8.0;
    EXPECT_NO_THROW(PushRun run(param));

    // Off, none of that is looked at.
    param.externalForceEnabled = false;
    param.externalForceShape = "ring";
    param.externalForceRadius = 0.0;
    EXPECT_NO_THROW(PushRun run(param));
}

TEST(ExternalForceTest, TheRegionFileListsThePushedVertices)
{
    Param param;
    configure_clamped_sheet(param);
    param.externalForceEnabled = true;
    param.externalForceZ = 1.0;
    param.externalForceRadius = 12.0;
    PushRun run(param);

    const std::string path = "test_external_force_vertices.tmp.csv";
    run.model->write_external_force_vertices_csv(path);
    std::vector<int> listed;
    {
        std::ifstream in(path);
        ASSERT_TRUE(in.is_open());
        int index = 0;
        while (in >> index)
        {
            listed.push_back(index);
        }
    }
    std::remove(path.c_str());

    std::vector<int> expected;
    for (int v = 0; v < static_cast<int>(run.model->externalForceMask.size()); v++)
    {
        if (run.model->externalForceMask[v] != 0)
        {
            expected.push_back(v);
        }
    }
    EXPECT_EQ(listed, expected);
    EXPECT_EQ(static_cast<int>(listed.size()), run.model->externalForceVertexCount);
}
