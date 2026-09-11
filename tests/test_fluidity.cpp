#include "test_fluidity.hpp"

#include <cstdint>

/**
 * Two things are under test here, and the second is the one WP6 exists for.
 *
 * **The tether.** A membrane needs something to stop its triangles
 * degenerating, and the harmonic spring WP4 shipped cannot be that and permit
 * flips at the same time: a flip on a rhombus of two equilateral triangles
 * replaces the short diagonal by the long one, so it has to climb
 * `(k/2)(sqrt(3)-1)^2 l0^2` whatever else the Hamiltonian says. The measured
 * consequence is in docs/edge_flip_plan.md work package 5 -- there is no
 * stiffness that is soft enough to flip and stiff enough to hold. The
 * flat-bottomed tether removes the competition by making the flip free inside
 * an allowed range, and the first tests below are that it does.
 *
 * **That the sweep and the dynamics agree about the Hamiltonian.** The flip
 * trial differences a local energy; the Brownian step integrates a global one.
 * If those are not the same functional then the Metropolis chain samples
 * neither, and nothing about the run means anything -- not the acceptance
 * rate, not the spectrum, not the diffusion. This is not hypothetical: when
 * the tether landed, the trial went on differencing the *reference-length*
 * regularization, and a 3000-step run reported every accepted flip as about
 * -750 pN.nm downhill while the mesh's total energy climbed. The test named
 * TheTrialAndTheDynamicsShareOneRegularization is what that cost.
 *
 * Then the fluidity itself: neighbours that turn over, and in-plane motion
 * that escapes its cage. Those are the properties the flip move was added for,
 * and the only ones that cannot be produced by a move that shuffles without
 * transporting.
 *
 * @see docs/edge_flip_plan.md work package 6
 */

namespace
{

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

void configure_sheet(Param &param, double side = 60.0)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Periodic;
    param.sideX = side;
    param.sideY = side;
    param.lFace = 5.0;
    param.kCurv = 83.4;
    param.KBT = 4.17;
    param.uSurf = 250.0;
    param.uVol = 0.0;
    param.isGlobalConstraint = true;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
    param.randomSeed = 20260910u;
    param.timeStep = 1.0e-3;
    param.diffConst = 1.0;
    param.surfaceSolver = "iterative";
    param.edgeSpringEnabled = true;
    param.edgeSpringConstant = 83.4;
}

/// A flat sheet through the full DynamicMesh setup, ready to be stepped.
void build_sheet(DynamicMesh &mesh)
{
    mesh.setup_flat();
    for (Vertex &vertex : mesh.vertices)
    {
        vertex.coord.set(2, 0, 10.0);
    }
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.Compute_Energy_And_Force();
    mesh.update_previous_coord_for_vertex();
    mesh.update_previous_force_for_vertex();
    mesh.update_vertices_mat_with_vector();
}

double total_regularization(const Mesh &mesh)
{
    double sum = 0.0;
    for (const Face &face : mesh.faces)
    {
        sum += face.energy.energyRegularization;
    }
    return sum;
}

double total_energy(Mesh &mesh)
{
    mesh.Compute_Energy_And_Force();
    return mesh.param.energy.energyTotal;
}

/// Length of the edge between two vertices.
double edge_length(const Mesh &mesh, int a, int b)
{
    double sum = 0.0;
    for (int axis = 0; axis < 3; axis++)
    {
        const double d = mesh.vertices[a].coord.get(axis, 0) - mesh.vertices[b].coord.get(axis, 0);
        sum += d * d;
    }
    return std::sqrt(sum);
}

/// The undirected edges of the mesh, as a set that survives a relabelling.
std::set<std::pair<int, int>> edge_set(const Mesh &mesh)
{
    std::set<std::pair<int, int>> edges;
    for (const Face &face : mesh.faces)
    {
        for (int k = 0; k < 3; k++)
        {
            const int a = face.adjacentVertices[k];
            const int b = face.adjacentVertices[(k + 1) % 3];
            edges.emplace(std::min(a, b), std::max(a, b));
        }
    }
    return edges;
}

int first_flippable_edge(const Mesh &mesh)
{
    for (const MeshEdge &edge : mesh.edges)
    {
        if (edge.flippable)
        {
            return edge.index;
        }
    }
    return -1;
}

/**
 * @brief Set edgeFlipAttemptRate so a sweep draws about @p lambda attempts.
 *
 * The rate is per edge per unit time, so the mean also scales with the number
 * of flippable edges and the time step. A test cares about the count.
 */
void set_rate_for_expected_attempts(Mesh &mesh, double lambda)
{
    int flippable = 0;
    for (const MeshEdge &edge : mesh.edges)
    {
        flippable += edge.flippable ? 1 : 0;
    }
    const double interval = std::max(1, mesh.param.edgeFlipInterval);
    mesh.param.edgeFlipAttemptRate =
        (flippable > 0) ? lambda / (mesh.param.timeStep * interval * flippable) : 0.0;
}

/// Vertices whose fan of faces closes and which the step actually integrates.
std::vector<int> interior_free_vertices(const DynamicMesh &mesh)
{
    std::vector<int> keep;
    for (int v = 0; v < static_cast<int>(mesh.vertices.size()); v++)
    {
        if (mesh.vertices[v].isGhost)
        {
            continue;
        }
        if (!mesh.isSlavedPeriodic.empty() && mesh.isSlavedPeriodic[v])
        {
            continue;
        }
        if (mesh.vertices[v].adjacentFaces.size() != mesh.vertices[v].adjacentVertices.size())
        {
            continue; // a partial fan: on the generated sheet, its outer rim
        }
        keep.push_back(v);
    }
    return keep;
}

/// One Brownian step plus, optionally, one flip sweep. The driver's loop body.
struct Trajectory
{
    std::vector<double> startX, startY;
    std::vector<int> tagged;
    int accepted = 0;
    double meanSquaredDisplacement = 0.0;
};

Trajectory run_steps(DynamicMesh &mesh, DynamicModel &model, Record &record, int nSteps)
{
    Trajectory out;
    out.tagged = interior_free_vertices(mesh);
    for (int v : out.tagged)
    {
        out.startX.push_back(mesh.vertices[v].coord.get(0, 0));
        out.startY.push_back(mesh.vertices[v].coord.get(1, 0));
    }

    for (model.iteration = 0; model.iteration < nSteps; model.iteration++)
    {
        mesh.apply_mesh_to_surface();
        model.next_step();
        mesh.apply_surface_to_mesh();
        mesh.postprocess_ghost_periodic();
        mesh.update_vertices_vector_with_mat();
        record.add(mesh.param.area, mesh.param.energy, mesh.calculate_mean_force());
        mesh.Compute_Energy_And_Force();

        if (mesh.param.edgeFlipEnabled)
        {
            const EdgeFlipSweepStats stats = mesh.edge_flip_sweep(model.iteration);
            out.accepted += stats.accepted;
            if (stats.accepted > 0)
            {
                mesh.Compute_Energy_And_Force();
            }
        }
    }

    // Against the sheet's own drift: the whole membrane translating is a
    // centre-of-mass motion, not a rearrangement, and on a periodic box it is
    // not even physical.
    double driftX = 0.0, driftY = 0.0;
    for (std::size_t i = 0; i < out.tagged.size(); i++)
    {
        driftX += mesh.vertices[out.tagged[i]].coord.get(0, 0) - out.startX[i];
        driftY += mesh.vertices[out.tagged[i]].coord.get(1, 0) - out.startY[i];
    }
    driftX /= out.tagged.size();
    driftY /= out.tagged.size();

    for (std::size_t i = 0; i < out.tagged.size(); i++)
    {
        const double dx = mesh.vertices[out.tagged[i]].coord.get(0, 0) - out.startX[i] - driftX;
        const double dy = mesh.vertices[out.tagged[i]].coord.get(1, 0) - out.startY[i] - driftY;
        out.meanSquaredDisplacement += dx * dx + dy * dy;
    }
    out.meanSquaredDisplacement /= out.tagged.size();
    return out;
}

} // namespace

// ---------------------------------------------------------------------------
// The tether's shape
// ---------------------------------------------------------------------------

/**
 * @brief Inside the allowed range the tether is exactly zero.
 *
 * Not "small": zero. That is the whole difference from the spring, and it is
 * what makes the flip barrier vanish rather than merely shrink.
 */
TEST(FlatTetherTest, TheTetherIsExactlyZeroInsideItsRange)
{
    Param param;
    configure_sheet(param);
    param.edgeTetherMinRatio = 0.6;
    param.edgeTetherMaxRatio = 1.8;

    DynamicMesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_sheet(mesh));

    // The generated sheet is equilateral at lFace, comfortably inside [3, 9].
    for (const MeshEdge &edge : mesh.edges)
    {
        const double length = edge_length(mesh, edge.v[0], edge.v[1]);
        ASSERT_GT(length, 0.6 * param.lFace) << "the fixture is not inside the range";
        ASSERT_LT(length, 1.8 * param.lFace) << "the fixture is not inside the range";
    }

    EXPECT_EQ(total_regularization(mesh), 0.0)
        << "a mesh entirely inside the tether's range carries tether energy";
}

/// Outside it, the tether is the quadratic wall it claims to be.
TEST(FlatTetherTest, TheTetherIsQuadraticOutsideItsRange)
{
    Param param;
    configure_sheet(param);
    param.edgeTetherMinRatio = 0.6;
    param.edgeTetherMaxRatio = 1.8;
    DynamicMesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_sheet(mesh));

    const double k = param.edgeSpringConstant;
    const double lower = 0.6 * param.lFace;
    const double upper = 1.8 * param.lFace;

    EXPECT_DOUBLE_EQ(mesh.edge_tether_energy(lower), 0.0);
    EXPECT_DOUBLE_EQ(mesh.edge_tether_energy(upper), 0.0);
    EXPECT_DOUBLE_EQ(mesh.edge_tether_energy(0.5 * (lower + upper)), 0.0);
    EXPECT_DOUBLE_EQ(mesh.edge_tether_energy(upper + 2.0), 0.5 * k * 4.0);
    EXPECT_DOUBLE_EQ(mesh.edge_tether_energy(lower - 0.5), 0.5 * k * 0.25);

    // And the harmonic shape is still the harmonic shape.
    mesh.param.edgeTetherShape = "harmonic";
    EXPECT_DOUBLE_EQ(mesh.edge_tether_energy(param.lFace + 1.0), 0.5 * k);
    EXPECT_DOUBLE_EQ(mesh.edge_tether_energy(param.lFace), 0.0);
}

/**
 * @brief The tether's force is minus the gradient of its energy, at the wall.
 *
 * The test WP4 ran was on the smooth spring. A piecewise term is where a sign
 * or a bound gets dropped, and inside the flat region the gradient is zero, so
 * a finite difference there would pass against nothing.
 */
TEST(FlatTetherTest, TheTetherForceIsMinusTheGradientOfItsEnergyAtTheWall)
{
    Param param;
    configure_sheet(param);
    param.kCurv = 0.0; // the tether alone, so nothing else can compensate
    param.uSurf = 0.0;
    param.isGlobalConstraint = false;
    // A range narrow enough that the sheet's own edges are already outside it,
    // which is what puts every edge on a wall where the gradient is nonzero.
    param.edgeTetherMinRatio = 1.2;
    param.edgeTetherMaxRatio = 1.4;

    DynamicMesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_sheet(mesh));
    ASSERT_GT(total_regularization(mesh), 0.0) << "the fixture is not against a wall";

    const double h = 1e-6;
    for (int iVertex : {40, 77, 120})
    {
        ASSERT_LT(iVertex, static_cast<int>(mesh.vertices.size()));
        for (int axis = 0; axis < 3; axis++)
        {
            const double original = mesh.vertices[iVertex].coord.get(axis, 0);

            mesh.vertices[iVertex].coord.set(axis, 0, original + h);
            mesh.energy_force_edge_spring();
            const double plus = total_regularization(mesh);

            mesh.vertices[iVertex].coord.set(axis, 0, original - h);
            mesh.energy_force_edge_spring();
            const double minus = total_regularization(mesh);

            mesh.vertices[iVertex].coord.set(axis, 0, original);
            mesh.energy_force_edge_spring();
            const double analytic = mesh.vertices[iVertex].force.forceRegularization.get(axis, 0);

            const double numeric = -(plus - minus) / (2.0 * h);
            EXPECT_NEAR(analytic, numeric, 1e-4 * std::max(1.0, std::abs(numeric)))
                << "vertex " << iVertex << " axis " << axis;
        }
    }
}

/// A range with nothing in it would put every edge against a wall.
TEST(FlatTetherTest, AnEmptyRangeIsRefused)
{
    Param param;
    configure_sheet(param);
    param.edgeTetherMinRatio = 1.5;
    param.edgeTetherMaxRatio = 1.2;
    DynamicMesh mesh(param);
    mesh.setup_flat();
    EXPECT_THROW(mesh.energy_force_edge_spring(), std::runtime_error);
}

/// The parameters reach Param through the path a params file takes.
TEST(FlatTetherTest, TheTetherParametersRoundTrip)
{
    Param param;
    EXPECT_EQ(param.edgeTetherShape, "flat");
    EXPECT_TRUE(import_kv_string("edgeTetherShape", "harmonic", param));
    EXPECT_EQ(param.edgeTetherShape, "harmonic");
    EXPECT_TRUE(import_kv_string("edgeTetherMinRatio", "0.75", param));
    EXPECT_DOUBLE_EQ(param.edgeTetherMinRatio, 0.75);
    EXPECT_TRUE(import_kv_string("edgeTetherMaxRatio", "1.9", param));
    EXPECT_DOUBLE_EQ(param.edgeTetherMaxRatio, 1.9);
}

// ---------------------------------------------------------------------------
// The barrier the tether exists to remove
// ---------------------------------------------------------------------------

/**
 * @brief The claim the whole tether rests on, measured.
 *
 * A flip on the generated sheet replaces an edge of `lFace` with one of
 * `sqrt(3) lFace`. Under the harmonic spring that costs
 * `(k/2)(sqrt(3)-1)^2 lFace^2` -- 559 pN.nm at the shipped constants, 134 kT,
 * which is a membrane that never flips. Under the flat tether with the wall
 * above `sqrt(3) lFace` it costs nothing at all, and what is left to decide
 * the move is the bending energy and the constraints.
 */
TEST(FlatTetherTest, AnEquilateralFlipIsFreeUnderTheTetherAndImpossibleUnderTheSpring)
{
    const auto tether_change_of_one_flip = [](const std::string &shape, double maxRatio) {
        Param param;
        configure_sheet(param);
        param.kCurv = 0.0; // isolate the mesh-quality term
        param.uSurf = 0.0;
        param.isGlobalConstraint = false;
        param.edgeTetherShape = shape;
        param.edgeTetherMinRatio = 0.6;
        param.edgeTetherMaxRatio = maxRatio;

        DynamicMesh mesh(param);
        build_sheet(mesh);
        const int iEdge = first_flippable_edge(mesh);
        EXPECT_GE(iEdge, 0);

        const double before = total_regularization(mesh);
        const double lengthBefore = edge_length(mesh, mesh.edges[iEdge].v[0], mesh.edges[iEdge].v[1]);
        const int t0 = mesh.edges[iEdge].opposite[0];
        const int t1 = mesh.edges[iEdge].opposite[1];
        const double lengthAfter = edge_length(mesh, t0, t1);

        mesh.flip_edge(iEdge);
        mesh.energy_force_edge_spring();
        const double after = total_regularization(mesh);
        return std::array<double, 3>{{after - before, lengthBefore, lengthAfter}};
    };

    const std::array<double, 3> flat = tether_change_of_one_flip("flat", 1.8);
    const std::array<double, 3> spring = tether_change_of_one_flip("harmonic", 1.8);
    const std::array<double, 3> tooNarrow = tether_change_of_one_flip("flat", 1.5);

    // The geometry the argument rests on: the flip really does reach sqrt(3).
    EXPECT_NEAR(flat[1], 5.0, 1e-9);
    EXPECT_NEAR(flat[2], std::sqrt(3.0) * 5.0, 1e-6)
        << "a flip on this sheet does not produce the long diagonal, so the "
           "rest of this test is measuring something else";

    EXPECT_NEAR(flat[0], 0.0, 1e-12)
        << "the flat tether charged " << flat[0] << " pN.nm for a flip inside its range";

    const double expectedSpring =
        0.5 * 83.4 * std::pow(std::sqrt(3.0) - 1.0, 2.0) * 25.0;
    EXPECT_NEAR(spring[0], expectedSpring, 1e-6);
    EXPECT_GT(spring[0] / 4.17, 100.0)
        << "the harmonic barrier should be of order a hundred kT at these constants";

    // Below sqrt(3) the wall is back inside the flip's reach, which is the
    // configuration error the startup diagnostic warns about.
    EXPECT_GT(tooNarrow[0], 0.0)
        << "edgeTetherMaxRatio = 1.5 is below sqrt(3) and must still charge for the flip";
}

// ---------------------------------------------------------------------------
// One Hamiltonian
// ---------------------------------------------------------------------------

/**
 * @brief The flip trial and the force pass difference the same term.
 *
 * The trial evaluates the energy of eighteen faces twice; the dynamics
 * evaluates the whole mesh. Those are two implementations of one functional,
 * and if they drift apart the Metropolis chain samples a distribution that is
 * not the one being integrated -- silently, because both halves keep working.
 *
 * They did drift apart. When the tether landed the trial went on differencing
 * face_regularization_energy(), the reference-length term, which hands a newly
 * created edge the distance between two vertices that were never joined. A
 * 3000-step run then reported every accepted flip as about -750 pN.nm downhill
 * while the mesh's total energy climbed, and the acceptance rate read 17%
 * instead of 40%.
 */
TEST(FluidHamiltonianTest, TheTrialAndTheDynamicsShareOneRegularization)
{
    Param param;
    configure_sheet(param);
    param.edgeTetherMinRatio = 1.2; // a range the sheet is outside, so the
    param.edgeTetherMaxRatio = 1.4; // term is not identically zero
    DynamicMesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_sheet(mesh));

    const double fromForcePass = total_regularization(mesh);
    ASSERT_GT(fromForcePass, 0.0) << "the fixture carries no tether energy to compare";

    // The per-face split itself: an edge's energy shared between the two faces
    // it separates has to add back up to the sum over edges.
    double fromPerFaceTerm = 0.0;
    for (int iFace = 0; iFace < static_cast<int>(mesh.faces.size()); iFace++)
    {
        fromPerFaceTerm += mesh.face_tether_energy(iFace);
    }
    EXPECT_NEAR(fromPerFaceTerm, fromForcePass, 1e-9 * std::max(1.0, fromForcePass))
        << "the per-face tether term does not sum to the term the force pass computes";

    // And that the trial's evaluator actually reaches it. This is the half
    // that was wrong: face_tether_energy() was correct and unused, while
    // evaluate_face_subset() went on calling face_regularization_energy().
    std::vector<int> everyFace(mesh.faces.size());
    for (int iFace = 0; iFace < static_cast<int>(mesh.faces.size()); iFace++)
    {
        everyFace[iFace] = iFace;
    }
    const FaceSubsetEnergy subset = mesh.evaluate_face_subset(everyFace);
    EXPECT_NEAR(subset.regularization, fromForcePass, 1e-9 * std::max(1.0, fromForcePass))
        << "the trial's evaluator returned " << subset.regularization
        << " where the force pass computes " << fromForcePass
        << "; the sweep and the dynamics are differencing different terms";
}

/**
 * @brief With the tether on, a trial's dE is the change in the total energy.
 *
 * The end-to-end form of the test above, and the one that would have caught
 * the defect on its own: it compares the number the Metropolis sweep actually
 * uses against two whole-mesh evaluations.
 */
TEST(FluidHamiltonianTest, ATrialDeltaIsTheChangeInTheWholeMeshEnergy)
{
    Param param;
    configure_sheet(param);
    param.edgeTetherMinRatio = 1.2;
    param.edgeTetherMaxRatio = 1.4;
    DynamicMesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_sheet(mesh));

    const int iEdge = first_flippable_edge(mesh);
    ASSERT_GE(iEdge, 0);

    const double before = total_energy(mesh);

    EdgeFlipDelta delta;
    ASSERT_TRUE(mesh.evaluate_edge_flip(iEdge, delta));

    // evaluate_edge_flip() restores the mesh, so the energy has to be back too.
    EXPECT_NEAR(total_energy(mesh), before, 1e-9 * std::max(1.0, std::abs(before)));

    mesh.flip_edge(iEdge);
    mesh.param.area += delta.area;
    mesh.param.vol += delta.volume;
    const double after = total_energy(mesh);

    EXPECT_NEAR(delta.energy, after - before, 1e-6 * std::max(1.0, std::abs(delta.energy)))
        << "the sweep's dE was " << delta.energy << " but the mesh's energy moved by "
        << (after - before);
}

// ---------------------------------------------------------------------------
// Fluidity
// ---------------------------------------------------------------------------

/**
 * @brief Neighbours are permanent without flips and turn over with them.
 *
 * The microscopic definition of a fluid, and the cleanest of the three
 * measurements because the control is exact rather than statistical: with the
 * connectivity fixed, neighbour survival is 1 by construction.
 */
TEST(FluidityTest, NeighboursArePermanentWithoutFlipsAndTurnOverWithThem)
{
    const auto surviving_fraction = [](bool withFlips) {
        Param param;
        configure_sheet(param);
        param.inPlaneDynamicsEnabled = true;
        param.edgeFlipEnabled = withFlips;
        param.maxIterations = 250;

        DynamicMesh mesh(param);
        build_sheet(mesh);
        if (withFlips)
        {
            set_rate_for_expected_attempts(mesh, 6.0);
        }
        Record record(param.maxIterations);
        DynamicModel model(mesh, record);

        const std::set<std::pair<int, int>> initial = edge_set(mesh);
        run_steps(mesh, model, record, param.maxIterations);
        const std::set<std::pair<int, int>> final = edge_set(mesh);

        int survived = 0;
        for (const std::pair<int, int> &edge : initial)
        {
            survived += final.count(edge);
        }
        return static_cast<double>(survived) / initial.size();
    };

    EXPECT_DOUBLE_EQ(surviving_fraction(false), 1.0)
        << "the connectivity changed with flips off";

    const double fluid = surviving_fraction(true);
    std::cout << "[FluidityTest] neighbour survival after 250 steps: " << fluid << std::endl;
    EXPECT_LT(fluid, 0.99) << "250 steps of flips replaced almost no neighbours";
    EXPECT_GT(fluid, 0.5) << "over half the neighbours turned over in 250 steps, which is "
                             "fast enough to suspect the mesh is coming apart rather than "
                             "rearranging";
}

/**
 * @brief Without flips, in-plane motion has a cage and finds it.
 *
 * Half of the classical fluid signature, and the half a test can reach. A
 * vertex whose neighbours never change is tethered to a fixed set of them, so
 * its in-plane displacement saturates: the MSD stops growing. That saturation
 * is the control the fluid case is read against, and it has to be established
 * before the comparison means anything -- two curves that agree because
 * neither has left free diffusion agree about nothing.
 *
 * The other half -- that flips let a vertex escape the cage -- is a
 * measurement rather than an assertion, because the separation only opens once
 * the caged curve has flattened, which takes of order ten thousand steps. On
 * the 60 nm sheet the two MSDs agree to 5% out to lag 800, cross at about lag
 * 1500, and reach a ratio of 1.23 at lag 4800 with the caged curve flat at 7.6
 * nm^2 and the fluid one still climbing. Those numbers are in
 * docs/fluidity_results.md; asserting them here would mean a ten-minute test.
 */
TEST(FluidityTest, InPlaneMotionIsCagedWithoutFlips)
{
    const auto displacement = [](bool withFlips, int nSteps) {
        Param param;
        configure_sheet(param);
        param.inPlaneDynamicsEnabled = true;
        param.edgeFlipEnabled = withFlips;
        param.maxIterations = nSteps;

        DynamicMesh mesh(param);
        build_sheet(mesh);
        if (withFlips)
        {
            set_rate_for_expected_attempts(mesh, 6.0);
        }
        Record record(nSteps + 1);
        DynamicModel model(mesh, record);
        return run_steps(mesh, model, record, nSteps);
    };

    // At the shipped tether range the cage is small enough to be filled
    // inside a test's step budget. It was not at the range WP6 started from:
    // [0.6, 1.8] gives a cage of about 7.7 nm^2 that takes some five thousand
    // steps to reach, so at any affordable length the vertex is still in free
    // diffusion and a saturation test measures nothing.
    const Trajectory quarter = displacement(false, 150);
    const Trajectory half = displacement(false, 300);
    const Trajectory whole = displacement(false, 600);

    // Printed because this is a measurement as much as an assertion, and a run
    // drifting toward the threshold should be visible before it crosses it.
    std::cout << "[FluidityTest] caged in-plane MSD (nm^2): " << quarter.meanSquaredDisplacement
              << " at 150 steps, " << half.meanSquaredDisplacement << " at 300, "
              << whole.meanSquaredDisplacement << " at 600." << std::endl;

    ASSERT_GT(quarter.meanSquaredDisplacement, 0.0) << "nothing moved in plane at all";

    // Free diffusion doubles the MSD when the time doubles. A vertex that has
    // filled its cage does not: it grows by 1.1 to 1.2 here, and the margin
    // against 2.0 is what says the cage is real rather than that the run was
    // too short to leave it.
    //
    // Only the endpoints are compared, so these ratios carry single-trajectory
    // noise -- enough that the two doublings do not come out reliably ordered.
    // The bound is the claim; the ordering is not.
    const double earlyGrowth = half.meanSquaredDisplacement / quarter.meanSquaredDisplacement;
    const double lateGrowth = whole.meanSquaredDisplacement / half.meanSquaredDisplacement;

    EXPECT_LT(earlyGrowth, 1.5)
        << "the in-plane MSD grew by " << earlyGrowth << "x over the first doubling";
    EXPECT_LT(lateGrowth, 1.5)
        << "the in-plane MSD grew by " << lateGrowth
        << "x over the last doubling, which is close to the 2.0 of free diffusion: this "
           "vertex has not found a cage, so it is not a control for anything";
}

/**
 * @brief At fixed coordinates the sweep leaves the energy stationary.
 *
 * The flip sweep on its own is a Metropolis kernel whose stationary
 * distribution is exp(-E/kT) over triangulations at fixed control net. So the
 * mean energy must not drift: it may wander, but not climb.
 *
 * This is the test that says where a drift lives. A fluid *run* -- flips
 * alternated with the Brownian step -- does drift: on the 60 nm sheet at
 * nu = 2, the tether energy climbs from 600 to 13800 over 20000 steps and the
 * mean control-net edge grows from 5.58 to 6.17 nm, while the same sheet with
 * flips off sits at a tether energy of about 50 indefinitely. If the sweep
 * alone were also to drift here, the fault would be in the sweep. It does not,
 * which puts it in the alternation -- see docs/fluidity_results.md.
 *
 * The coordinates are disordered first. On the pristine equilateral sheet
 * every flip produces the same length and the tether is flat across all of
 * them, so a stationarity test there would be passing against a term that is
 * identically zero.
 */
TEST(FluidityTest, TheSweepAloneLeavesTheEnergyStationary)
{
    Param param;
    configure_sheet(param);
    param.edgeFlipEnabled = true;
    DynamicMesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_sheet(mesh));

    // A reproducible disorder of about a tenth of the edge length, which is
    // enough to put a spread of lengths on the flips without tearing anything.
    std::uint64_t state = 0x9E3779B97F4A7C15ULL;
    const auto next = [&state]() {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>((state >> 11) * (1.0 / 9007199254740992.0)) - 0.5;
    };
    for (Vertex &vertex : mesh.vertices)
    {
        if (vertex.isGhost)
        {
            continue;
        }
        for (int axis = 0; axis < 3; axis++)
        {
            vertex.coord.set(axis, 0, vertex.coord.get(axis, 0) + 0.5 * next());
        }
    }
    mesh.postprocess_ghost_periodic();
    mesh.calculate_element_area_volume();
    mesh.Compute_Energy_And_Force();

    set_rate_for_expected_attempts(mesh, 8.0);

    std::vector<double> trace;
    int accepted = 0;
    for (int sweep = 0; sweep < 400; sweep++)
    {
        accepted += mesh.edge_flip_sweep(sweep).accepted;
        mesh.Compute_Energy_And_Force();
        trace.push_back(mesh.param.energy.energyTotal);
    }
    ASSERT_GT(accepted, 100) << "too few flips landed for a stationarity test to mean anything";

    const std::size_t quarter = trace.size() / 4;
    const auto mean_over = [&](std::size_t from, std::size_t to) {
        double sum = 0.0;
        for (std::size_t i = from; i < to; i++)
        {
            sum += trace[i];
        }
        return sum / (to - from);
    };
    const double first = mean_over(0, quarter);
    const double last = mean_over(trace.size() - quarter, trace.size());

    // Against the series' own fluctuation, not against its magnitude. A
    // stationary chain still wanders, and a fixed percentage would be a
    // statement about how large the energy happens to be rather than about
    // whether it is drifting. The composite run this exists to contrast with
    // moves by about twenty times its own standard deviation, so the
    // distinction is not a close one.
    const double mean = mean_over(0, trace.size());
    double variance = 0.0;
    for (double value : trace)
    {
        variance += (value - mean) * (value - mean);
    }
    const double fluctuation = std::sqrt(variance / trace.size());

    std::cout << "[FluidityTest] sweep-only energy: " << first << " over the first quarter, "
              << last << " over the last, fluctuation " << fluctuation << ", across " << accepted
              << " accepted flips." << std::endl;

    ASSERT_GT(fluctuation, 0.0) << "the energy did not move at all, so nothing was sampled";
    EXPECT_LT(std::abs(last - first), 1.5 * fluctuation)
        << "the energy moved from " << first << " to " << last << " -- "
        << std::abs(last - first) / fluctuation
        << " times its own fluctuation -- under flips alone at fixed coordinates, so the "
           "sweep is not sampling the distribution it claims to";
}

/**
 * @brief More attempts per step means more accepted flips, proportionally.
 *
 * The attempt rate is a physical quantity -- attempts per edge per unit time --
 * and the sweep draws its count from a Poisson distribution with mean
 * `nu * dt * interval * N_flippable`. If the accepted count did not scale with
 * `nu` then `nu` would not be the fluidity knob it is documented as, and WP6's
 * calibration against a neighbour-survival time would have nothing to turn.
 */
TEST(FluidityTest, AcceptedFlipsScaleWithTheAttemptRate)
{
    const auto accepted_at = [](double lambdaPerSweep) {
        Param param;
        configure_sheet(param);
        param.inPlaneDynamicsEnabled = true;
        param.edgeFlipEnabled = true;
        param.maxIterations = 120;

        DynamicMesh mesh(param);
        build_sheet(mesh);
        set_rate_for_expected_attempts(mesh, lambdaPerSweep);
        Record record(param.maxIterations + 1);
        DynamicModel model(mesh, record);
        return run_steps(mesh, model, record, param.maxIterations).accepted;
    };

    const int slow = accepted_at(1.0);
    const int fast = accepted_at(8.0);

    std::cout << "[FluidityTest] accepted flips over 120 steps: " << slow
              << " at one attempt per sweep, " << fast << " at eight." << std::endl;
    ASSERT_GT(slow, 0) << "even the slow rate accepted nothing";
    EXPECT_GT(fast, 3 * slow)
        << "eight times the attempt rate gave " << fast << " accepted flips against " << slow
        << "; the rate is not controlling the fluidity";
}
