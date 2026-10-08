#include "test_bead.hpp"

/**
 * The bead is a sphere the membrane cannot enter, moving at constant velocity
 * inside an iteration window. Four things can go wrong: the keys might not
 * parse; the centre might be in the wrong place at the wrong time; the wall
 * might not be the harmonic wall it says, or might reach vertices the step
 * does not move; and the step might not actually keep the membrane out. Each
 * has a test. A run with the bead off must be the run it always was, which
 * the byte-identity gate on the shipped periodic workload covers.
 *
 * @see docs/bead_push_benchmark.md
 */

namespace
{

/// A clamped 60 nm sheet under the per-vertex boundary, lifted to z = 10.
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
struct BeadRun
{
    Param param;
    std::unique_ptr<DynamicMesh> mesh;
    Record record{16};
    std::unique_ptr<DynamicModel> model;

    explicit BeadRun(const Param &source) : param(source)
    {
        mesh.reset(new DynamicMesh(param));
        mesh->setup_flat();
        prepare(*mesh);
        record.add(mesh->param.area, mesh->param.energy, mesh->calculate_mean_force());
        model.reset(new DynamicModel(*mesh, record));
    }

    /// One step as the driver takes it, coordinates mapped back to the net.
    void step(int iteration)
    {
        mesh->apply_mesh_to_surface();
        model->iteration = iteration;
        model->next_step();
        mesh->apply_surface_to_mesh();
        mesh->postprocess_boundary();
        mesh->update_vertices_vector_with_mat();
        mesh->Compute_Energy_And_Force();
    }

    /// One step leaving matSurface as the step produced it, for a comparison
    /// against a reference run with the same seed.
    void step_surface_only(int iteration)
    {
        mesh->apply_mesh_to_surface();
        model->iteration = iteration;
        model->next_step();
    }
};

double distance_to(const Vertex &vertex, const std::array<double, 3> &centre)
{
    const double dx = vertex.coord(0, 0) - centre[0];
    const double dy = vertex.coord(1, 0) - centre[1];
    const double dz = vertex.coord(2, 0) - centre[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

} // namespace

TEST(BeadTest, TheParameterKeysParse)
{
    Param param;
    EXPECT_FALSE(param.beadEnabled);
    EXPECT_TRUE(import_kv_string("beadEnabled", "true", param));
    EXPECT_TRUE(param.beadEnabled);
    EXPECT_TRUE(import_kv_string("beadRadius", "12.5", param));
    EXPECT_DOUBLE_EQ(param.beadRadius, 12.5);
    EXPECT_TRUE(import_kv_string("beadStiffness", "150", param));
    EXPECT_DOUBLE_EQ(param.beadStiffness, 150.0);
    EXPECT_TRUE(import_kv_string("beadCenterX", "1.5", param));
    EXPECT_DOUBLE_EQ(param.beadCenterX, 1.5);
    EXPECT_TRUE(import_kv_string("beadCenterY", "-2.5", param));
    EXPECT_DOUBLE_EQ(param.beadCenterY, -2.5);
    EXPECT_TRUE(import_kv_string("beadCenterZ", "-7", param));
    EXPECT_DOUBLE_EQ(param.beadCenterZ, -7.0);
    EXPECT_TRUE(import_kv_string("beadVelocityX", "0.1", param));
    EXPECT_DOUBLE_EQ(param.beadVelocityX, 0.1);
    EXPECT_TRUE(import_kv_string("beadVelocityY", "-0.2", param));
    EXPECT_DOUBLE_EQ(param.beadVelocityY, -0.2);
    EXPECT_TRUE(import_kv_string("beadVelocityZ", "0.5", param));
    EXPECT_DOUBLE_EQ(param.beadVelocityZ, 0.5);
    EXPECT_TRUE(import_kv_string("beadStartIteration", "10", param));
    EXPECT_EQ(param.beadStartIteration, 10);
    EXPECT_TRUE(import_kv_string("beadStopIteration", "5000", param));
    EXPECT_EQ(param.beadStopIteration, 5000);
    EXPECT_TRUE(import_kv_string("beadEnabled", "false", param));
    EXPECT_FALSE(param.beadEnabled);
}

TEST(BeadTest, TheCentreMovesOnlyInsideItsWindow)
{
    Param param;
    configure_clamped_sheet(param);
    param.beadEnabled = true;
    param.beadRadius = 8.0;
    param.beadCenterX = 1.0;
    param.beadCenterY = 2.0;
    param.beadCenterZ = -30.0; // well clear of the sheet
    param.beadVelocityX = -0.2;
    param.beadVelocityZ = 0.5;
    param.beadStartIteration = 100;
    param.beadStopIteration = 300;
    BeadRun run(param);
    const DynamicModel &model = *run.model;

    const auto expect_centre = [&](int iteration, double x, double y, double z) {
        const std::array<double, 3> c = model.bead_center_at(iteration);
        EXPECT_NEAR(c[0], x, 1e-12) << iteration;
        EXPECT_NEAR(c[1], y, 1e-12) << iteration;
        EXPECT_NEAR(c[2], z, 1e-12) << iteration;
    };
    expect_centre(0, 1.0, 2.0, -30.0);
    expect_centre(100, 1.0, 2.0, -30.0);
    // 50 steps of 0.001 us at (-0.2, 0, 0.5) nm/us
    expect_centre(150, 1.0 - 0.01, 2.0, -30.0 + 0.025);
    expect_centre(300, 1.0 - 0.04, 2.0, -30.0 + 0.1);
    expect_centre(1000, 1.0 - 0.04, 2.0, -30.0 + 0.1);

    Param forever = param;
    forever.beadStopIteration = -1;
    BeadRun open(forever);
    const std::array<double, 3> c = open.model->bead_center_at(1000);
    EXPECT_NEAR(c[0], 1.0 - 0.2 * 0.9, 1e-12);
    EXPECT_NEAR(c[2], -30.0 + 0.5 * 0.9, 1e-12);

    // Nothing touches the sheet from down there.
    run.model->update_bead(150);
    EXPECT_EQ(model.beadContacts, 0);
    EXPECT_DOUBLE_EQ(model.beadEnergy, 0.0);
}

/**
 * @brief The wall is (k/2)(R - d)^2 on every free vertex inside the sphere,
 * radial, and nothing on anything else; the reaction is minus its sum, and
 * the force is minus the gradient of the energy.
 */
TEST(BeadTest, TheWallIsHarmonicRadialAndOnlyOnFreeVertices)
{
    Param param;
    configure_clamped_sheet(param);
    param.beadEnabled = true;
    param.beadRadius = 8.0;
    param.beadStiffness = 150.0;
    param.beadCenterX = 2.5;
    param.beadCenterY = -2.0;
    param.beadCenterZ = 10.0 - 8.0 + 3.0; // the sheet cuts the sphere 3 nm below its top
    BeadRun run(param);
    DynamicModel &model = *run.model;
    const DynamicMesh &mesh = *run.mesh;
    model.update_bead(0);

    const std::array<double, 3> centre = {{2.5, -2.0, 5.0}};
    int contacts = 0;
    double energy = 0.0;
    std::array<double, 3> reaction = {{0.0, 0.0, 0.0}};
    double deepest = 0.0;
    int deepestVertex = -1;
    for (int v = 0; v < static_cast<int>(mesh.vertices.size()); v++)
    {
        const Vertex &vertex = mesh.vertices[v];
        const double d = distance_to(vertex, centre);
        const bool pushed = mesh.is_independent_vertex(v) && d < 8.0;
        if (!pushed)
        {
            EXPECT_EQ(model.beadForce(v, 0), 0.0) << v;
            EXPECT_EQ(model.beadForce(v, 1), 0.0) << v;
            EXPECT_EQ(model.beadForce(v, 2), 0.0) << v;
            continue;
        }
        contacts++;
        const double depth = 8.0 - d;
        energy += 0.5 * 150.0 * depth * depth;
        for (int axis = 0; axis < 3; axis++)
        {
            const double u = (vertex.coord(axis, 0) - centre[axis]) / d;
            EXPECT_NEAR(model.beadForce(v, axis), 150.0 * depth * u, 1e-9) << v << " " << axis;
            reaction[axis] -= 150.0 * depth * u;
        }
        if (depth > deepest)
        {
            deepest = depth;
            deepestVertex = v;
        }
    }
    EXPECT_GE(contacts, 3);
    EXPECT_EQ(model.beadContacts, contacts);
    EXPECT_NEAR(model.beadEnergy, energy, 1e-9);
    EXPECT_NEAR(model.beadMaxPenetration, deepest, 1e-12);
    for (int axis = 0; axis < 3; axis++)
    {
        EXPECT_NEAR(model.beadReaction[axis], reaction[axis], 1e-9) << axis;
    }
    // The sheet is above the centre, so the membrane pushes the bead down.
    EXPECT_LT(model.beadReaction[2], 0.0);

    // -dE/dz at the deepest vertex, by central difference.
    ASSERT_GE(deepestVertex, 0);
    Vertex &probe = run.mesh->vertices[deepestVertex];
    const double z0 = probe.coord(2, 0);
    const double h = 1e-4;
    probe.coord.set(2, 0, z0 + h);
    model.refresh_bead_contact();
    const double up = model.beadEnergy;
    probe.coord.set(2, 0, z0 - h);
    model.refresh_bead_contact();
    const double down = model.beadEnergy;
    probe.coord.set(2, 0, z0);
    model.refresh_bead_contact();
    EXPECT_NEAR(-(up - down) / (2.0 * h), model.beadForce(deepestVertex, 2), 1e-5);
}

/**
 * @brief The step moves the contacts out: against a reference run with the
 * same seed, one step lifts every contact vertex, and a few hundred steps
 * leave nothing more than a fraction of a nanometre inside the wall.
 */
TEST(BeadTest, TheStepKeepsTheMembraneOutOfTheBead)
{
    Param off;
    configure_clamped_sheet(off);
    Param on = off;
    on.beadEnabled = true;
    on.beadRadius = 8.0;
    on.beadStiffness = 200.0;
    on.beadCenterZ = 10.0 - 8.0 + 3.0;

    BeadRun reference(off);
    BeadRun pushed(on);
    reference.step_surface_only(0);
    pushed.step_surface_only(0);
    int lifted = 0;
    for (int v = 0; v < static_cast<int>(reference.mesh->vertices.size()); v++)
    {
        const double dz = pushed.mesh->matSurface(v, 2) - reference.mesh->matSurface(v, 2);
        if (pushed.model->beadForce(v, 2) > 0.0)
        {
            EXPECT_GT(dz, 0.0) << v;
            lifted++;
        }
    }
    EXPECT_EQ(lifted, pushed.model->beadContacts);
    EXPECT_GE(lifted, 3);

    // The wall relaxes at k D dt / kT = 0.048 per step; 400 steps is plenty.
    BeadRun held(on);
    for (int iteration = 0; iteration < 400; iteration++)
    {
        held.step(iteration);
    }
    held.model->refresh_bead_contact();
    EXPECT_GT(held.model->beadContacts, 0);
    EXPECT_LT(held.model->beadMaxPenetration, 0.5);
    EXPECT_LT(held.model->beadReaction[2], 0.0);
    EXPECT_TRUE(std::isfinite(held.mesh->param.energy.energyTotal));
    // Every vertex the wall reached sits at or just inside the sphere; none
    // is left deep inside.
    const std::array<double, 3> centre = {{0.0, 0.0, 5.0}};
    for (int v = 0; v < static_cast<int>(held.mesh->vertices.size()); v++)
    {
        if (held.mesh->is_independent_vertex(v))
        {
            EXPECT_GT(distance_to(held.mesh->vertices[v], centre), 8.0 - 0.5) << v;
        }
    }
}

TEST(BeadTest, ABadBeadIsRefused)
{
    Param param;
    configure_clamped_sheet(param);
    param.beadEnabled = true;
    param.beadCenterZ = -30.0;
    param.beadRadius = 0.0;
    EXPECT_THROW(BeadRun run(param), std::runtime_error);
    param.beadRadius = 8.0;
    param.beadStiffness = -1.0;
    EXPECT_THROW(BeadRun run(param), std::runtime_error);
    param.beadStiffness = 200.0;
    EXPECT_NO_THROW(BeadRun run(param));
    // Off, none of that is looked at.
    param.beadEnabled = false;
    param.beadRadius = 0.0;
    EXPECT_NO_THROW(BeadRun run(param));
}

TEST(BeadTest, TheLogHasAHeaderAndOneRowPerCall)
{
    Param param;
    configure_clamped_sheet(param);
    param.beadEnabled = true;
    param.beadRadius = 8.0;
    param.beadStiffness = 150.0;
    param.beadCenterZ = 10.0 - 8.0 + 3.0;
    param.beadVelocityZ = 0.5;
    BeadRun run(param);
    DynamicModel &model = *run.model;

    const std::string path = "test_bead_log.tmp.csv";
    std::remove(path.c_str());
    model.update_bead(0);
    model.append_bead_log(path, 0);
    model.update_bead(200);
    model.append_bead_log(path, 200);

    std::vector<std::string> lines;
    {
        std::ifstream in(path);
        ASSERT_TRUE(in.is_open());
        std::string line;
        while (std::getline(in, line))
        {
            lines.push_back(line);
        }
    }
    std::remove(path.c_str());
    ASSERT_EQ(lines.size(), 3u);
    EXPECT_EQ(lines[0],
              "iteration,time_us,centerX,centerY,centerZ,energy_pNnm,contacts,maxPenetration_nm,"
              "forceOnBeadX,forceOnBeadY,forceOnBeadZ,adhesion_pNnm,adhered");

    const auto fields = [](const std::string &line) {
        std::vector<double> out;
        std::stringstream stream(line);
        std::string item;
        while (std::getline(stream, item, ','))
        {
            out.push_back(std::stod(item));
        }
        return out;
    };
    const std::vector<double> first = fields(lines[1]);
    const std::vector<double> second = fields(lines[2]);
    ASSERT_EQ(first.size(), 13u);
    ASSERT_EQ(second.size(), 13u);
    EXPECT_DOUBLE_EQ(first[0], 0.0);
    EXPECT_DOUBLE_EQ(first[1], 0.0);
    EXPECT_NEAR(first[4], 5.0, 1e-9);
    EXPECT_GT(first[6], 0.0); // contacts
    EXPECT_DOUBLE_EQ(second[0], 200.0);
    EXPECT_NEAR(second[1], 0.2, 1e-9);
    EXPECT_NEAR(second[4], 5.0 + 0.5 * 0.2, 1e-9);
    EXPECT_NEAR(second[6], model.beadContacts, 1e-12);
    EXPECT_NEAR(second[10], model.beadReaction[2], 1e-6);
    EXPECT_DOUBLE_EQ(second[11], 0.0); // no adhesion in this run
    EXPECT_NEAR(second[12], model.beadContacts, 1e-12);
}

TEST(BeadTest, TheParameterKeysForAdhesionParse)
{
    Param param;
    EXPECT_DOUBLE_EQ(param.beadAdhesionEnergy, 0.0);
    EXPECT_TRUE(import_kv_string("beadAdhesionEnergy", "60", param));
    EXPECT_DOUBLE_EQ(param.beadAdhesionEnergy, 60.0);
    EXPECT_TRUE(import_kv_string("beadAdhesionRange", "2.5", param));
    EXPECT_DOUBLE_EQ(param.beadAdhesionRange, 2.5);
}

/**
 * @brief The well is -w s((d - R) / delta) on every free vertex within the
 * range, its force is minus its gradient, and it pulls toward the centre;
 * nothing beyond the range feels it, and inside the sphere it is the full
 * depth on top of the wall.
 */
TEST(BeadTest, TheAdhesionWellPullsSmoothlyTowardTheBead)
{
    Param param;
    configure_clamped_sheet(param);
    param.beadEnabled = true;
    param.beadRadius = 8.0;
    param.beadStiffness = 150.0;
    param.beadAdhesionEnergy = 60.0;
    param.beadAdhesionRange = 3.0;
    param.beadCenterX = 1.0;
    param.beadCenterY = 0.5;
    param.beadCenterZ = 10.0 - 8.0 - 1.0; // the sheet is 1 nm above the top of the sphere
    BeadRun run(param);
    DynamicModel &model = *run.model;
    const DynamicMesh &mesh = *run.mesh;
    model.update_bead(0);

    const std::array<double, 3> centre = {{1.0, 0.5, 1.0}};
    int adhered = 0;
    double adhesion = 0.0;
    int deepestVertex = -1;
    double nearest = 1e300;
    for (int v = 0; v < static_cast<int>(mesh.vertices.size()); v++)
    {
        const Vertex &vertex = mesh.vertices[v];
        const double d = distance_to(vertex, centre);
        if (!mesh.is_independent_vertex(v) || d >= 8.0 + 3.0)
        {
            EXPECT_EQ(model.beadForce(v, 2), 0.0) << v;
            continue;
        }
        ASSERT_GE(d, 8.0) << "the sheet should be clear of the wall";
        adhered++;
        const double x = (d - 8.0) / 3.0;
        const double s = 1.0 - 3.0 * x * x + 2.0 * x * x * x;
        adhesion -= 60.0 * s;
        const double f = 60.0 * (-6.0 * x + 6.0 * x * x) / 3.0;
        EXPECT_LE(f, 0.0);
        for (int axis = 0; axis < 3; axis++)
        {
            const double u = (vertex.coord(axis, 0) - centre[axis]) / d;
            EXPECT_NEAR(model.beadForce(v, axis), f * u, 1e-9) << v << " " << axis;
        }
        // The sheet is above the bead, so the pull is downward.
        EXPECT_LT(model.beadForce(v, 2), 0.0) << v;
        if (d < nearest)
        {
            nearest = d;
            deepestVertex = v;
        }
    }
    EXPECT_GE(adhered, 3);
    EXPECT_EQ(model.beadAdhered, adhered);
    EXPECT_EQ(model.beadContacts, 0);
    EXPECT_NEAR(model.beadAdhesion, adhesion, 1e-9);
    EXPECT_LT(model.beadAdhesion, 0.0);
    // The membrane pulls the bead up toward itself.
    EXPECT_GT(model.beadReaction[2], 0.0);

    // -dE/dz at the vertex nearest the bead, by central difference on the
    // well's energy.
    ASSERT_GE(deepestVertex, 0);
    Vertex &probe = run.mesh->vertices[deepestVertex];
    const double z0 = probe.coord(2, 0);
    const double h = 1e-4;
    probe.coord.set(2, 0, z0 + h);
    model.refresh_bead_contact();
    const double up = model.beadAdhesion + model.beadEnergy;
    probe.coord.set(2, 0, z0 - h);
    model.refresh_bead_contact();
    const double down = model.beadAdhesion + model.beadEnergy;
    probe.coord.set(2, 0, z0);
    model.refresh_bead_contact();
    EXPECT_NEAR(-(up - down) / (2.0 * h), model.beadForce(deepestVertex, 2), 1e-5);

    // Inside the sphere the well is at full depth, on top of the wall.
    Param inside = param;
    inside.beadCenterZ = 10.0 - 8.0 + 2.0;
    BeadRun pushed(inside);
    pushed.model->update_bead(0);
    EXPECT_GT(pushed.model->beadContacts, 0);
    EXPECT_GE(pushed.model->beadAdhered, pushed.model->beadContacts);
    // Full depth on every contact, and between nothing and full depth on
    // each vertex in the range beyond them.
    EXPECT_LE(pushed.model->beadAdhesion, -60.0 * pushed.model->beadContacts + 1e-9);
    EXPECT_GE(pushed.model->beadAdhesion, -60.0 * pushed.model->beadAdhered - 1e-9);

    // A range that is not positive is refused when there is a well.
    Param bad = param;
    bad.beadAdhesionRange = 0.0;
    EXPECT_THROW(BeadRun run(bad), std::runtime_error);
}
