#include "test_vertex_move.hpp"

/**
 * A Monte Carlo move is only worth having if it samples the distribution it
 * claims to, and for this one that comes down to three things, each of which
 * fails silently:
 *
 * **The energy difference.** The trial differences two dozen faces instead
 * of the whole mesh. If that local number is not the change of the global
 * energy -- at a periodic seam especially, where a vertex carries copies --
 * the chain samples some other Hamiltonian and every acceptance rate still
 * looks healthy.
 *
 * **The walls.** Hard walls replace three energies by three constraints at
 * the same positions. They have to sit exactly where the soft terms start
 * charging, be checked before a move is taken, and be respected by the flip
 * as well as by the vertex move, or the "admissible set" is not a set.
 *
 * **The stationary distribution.** With one vertex free to move along one
 * axis the Boltzmann distribution is a one-dimensional integral that can be
 * done by quadrature, and the chain has to reproduce it.
 *
 * @see docs/monte_carlo_vertex_moves.md
 */

namespace
{

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

/// A generated flat sheet. 70 nm gives a periodic tile of 8 x 10 vertices:
/// small enough to be quick, wide enough that a vertex and its own copy
/// never share a face.
void configure_sheet(Param &param, BoundaryType boundary)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = boundary;
    param.boundaryConditionX = BoundaryType::Periodic;
    param.boundaryConditionY = BoundaryType::Periodic;
    param.sideX = 70.0;
    param.sideY = 70.0;
    param.lFace = 5.0;
    param.kCurv = 83.4;
    param.KBT = 4.17;
    param.uSurf = 250.0;
    param.uVol = 0.0;
    param.isGlobalConstraint = true;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
    param.randomSeed = 20261001u;
    param.timeStep = 1.0e-3;
    param.diffConst = 1.0;
    param.surfaceSolver = "iterative";
}

/// The fluid mesh-quality term, with the range a fluid run uses.
void configure_fluid_walls(Param &param)
{
    param.inPlaneDynamicsEnabled = true;
    param.edgeSpringEnabled = true;
    param.edgeTetherShape = "flat";
    param.edgeTetherMinRatio = 0.7;
    param.edgeTetherMaxRatio = 1.75;
    param.edgeSpringConstant = 83.4;
    param.triangleShapeEnabled = true;
    param.triangleShapeMinAltitudeRatio = 0.4;
    param.triangleShapeConstant = 83.4;
    param.creaseWallEnabled = true;
    param.creaseWallAngle = 60.0;
    param.creaseWallConstant = 500.0;
}

/// Hard walls, and the pure Monte Carlo run they require.
void configure_hard_walls(Param &param)
{
    param.meshWallType = "hard";
    param.brownianStepEnabled = false;
    param.vertexMoveEnabled = true;
}

/// Everything the run does before its first sweep, on the flat sheet.
void build_sheet(DynamicMesh &mesh)
{
    mesh.setup_flat();
    for (Vertex &vertex : mesh.vertices)
    {
        vertex.coord.set(2, 0, 10.0);
    }
    mesh.sync_periodic_images();
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.param.vol0 = 0.0;
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.build_vertex_move_table();
    mesh.Compute_Energy_And_Force();
}

/// A deterministic generator, so a failure is reproducible.
struct Lcg
{
    unsigned long long state;
    explicit Lcg(unsigned long long seed) : state(seed) {}
    /// Uniform in [-1, 1).
    double symmetric()
    {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>(state >> 11) * (2.0 / 9007199254740992.0) - 1.0;
    }
    int below(int bound)
    {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<int>((state >> 33) % static_cast<unsigned int>(bound));
    }
};

/**
 * @brief Knock every independent vertex off its lattice site, copies with it.
 *
 * Through displace_vertex_with_followers(), so the periodic copies stay
 * exactly one period from their sources and the roughened sheet is still a
 * periodic one.
 */
void roughen(Mesh &mesh, double inPlane, double outOfPlane, unsigned long long seed)
{
    Lcg random(seed);
    for (int v : mesh.vertexMoveCandidates)
    {
        const double d[3] = {inPlane * random.symmetric(), inPlane * random.symmetric(),
                             outOfPlane * random.symmetric()};
        mesh.displace_vertex_with_followers(v, d);
    }
    mesh.Compute_Energy_And_Force();
}

double total_energy(Mesh &mesh)
{
    mesh.Compute_Energy_And_Force();
    return mesh.param.energy.energyTotal;
}

/// The part of the total the Brownian step integrates in solid mode.
double bending_and_area_energy(Mesh &mesh)
{
    mesh.Compute_Energy_And_Force();
    return mesh.param.energy.energyCurvature + mesh.param.energy.energyArea;
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

std::vector<double> all_coordinates(const Mesh &mesh)
{
    std::vector<double> coords;
    coords.reserve(mesh.vertices.size() * 3);
    for (const Vertex &vertex : mesh.vertices)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            coords.push_back(vertex.coord.get(axis, 0));
        }
    }
    return coords;
}

/// A candidate that carries periodic copies, and one that carries none.
void pick_seam_and_interior(const Mesh &mesh, int &seam, int &interior)
{
    seam = -1;
    interior = -1;
    for (int v : mesh.vertexMoveCandidates)
    {
        if (!mesh.vertexMoveFollowers[v].empty() && seam < 0)
        {
            seam = v;
        }
        if (mesh.vertexMoveFollowers[v].empty())
        {
            interior = v; // the last one: well inside the tile
        }
    }
}

const BoundaryType kPeriodicModes[2] = {BoundaryType::Periodic, BoundaryType::Mixed};

const char *mode_name(BoundaryType boundary)
{
    return (boundary == BoundaryType::Mixed) ? "Mixed" : "Periodic";
}

} // namespace

// ---------------------------------------------------------------------------
// The table: who moves, and who follows
// ---------------------------------------------------------------------------

/**
 * @brief Every coordinate of the periodic tile has exactly one owner.
 *
 * Under the global Periodic mode the Brownian step treats both copies of the
 * duplicated ring as slaved to each other; a Monte Carlo move that inherited
 * that would never move the ring at all. The table gives the tile's copy the
 * ownership, so the candidates are the whole tile and nothing else.
 */
TEST(VertexMoveTest, TheCandidatesAreThePeriodicTileAndEveryCopyFollowsOne)
{
    for (BoundaryType boundary : kPeriodicModes)
    {
        SCOPED_TRACE(mode_name(boundary));
        Param param;
        configure_sheet(param, boundary);
        DynamicMesh mesh(param);
        build_sheet(mesh);

        const int tile = (mesh.param.nFaceX - 6) * (mesh.param.nFaceY - 6);
        EXPECT_EQ(static_cast<int>(mesh.vertexMoveCandidates.size()), tile);

        // No vertex follows two sources, and no candidate follows anything.
        std::vector<int> owners(mesh.vertices.size(), 0);
        for (int v : mesh.vertexMoveCandidates)
        {
            for (int follower : mesh.vertexMoveFollowers[v])
            {
                owners[follower]++;
            }
        }
        for (int v : mesh.vertexMoveCandidates)
        {
            EXPECT_EQ(owners[v], 0) << "candidate " << v << " is also a follower";
        }
        for (int v = 0; v < static_cast<int>(mesh.vertices.size()); v++)
        {
            EXPECT_LE(owners[v], 1) << "vertex " << v << " follows more than one source";
        }

        // A copy sits a whole number of periods from its source, and stays
        // there when the source moves.
        int seam = -1;
        int interior = -1;
        pick_seam_and_interior(mesh, seam, interior);
        ASSERT_GE(seam, 0);
        ASSERT_GE(interior, 0);
        const int follower = mesh.vertexMoveFollowers[seam].front();
        double offsetBefore[3];
        for (int axis = 0; axis < 3; axis++)
        {
            offsetBefore[axis] = mesh.vertices[follower].coord.get(axis, 0) -
                                 mesh.vertices[seam].coord.get(axis, 0);
        }
        const double lx = (mesh.param.nFaceX - 6) * mesh.param.dFaceX;
        const double ly = (mesh.param.nFaceY - 6) * mesh.param.dFaceY;
        EXPECT_NEAR(std::remainder(offsetBefore[0], lx), 0.0, 1e-9);
        EXPECT_NEAR(std::remainder(offsetBefore[1], ly), 0.0, 1e-9);
        EXPECT_NEAR(offsetBefore[2], 0.0, 1e-12);

        const double d[3] = {0.31, -0.17, 0.23};
        mesh.displace_vertex_with_followers(seam, d);
        for (int axis = 0; axis < 3; axis++)
        {
            EXPECT_NEAR(mesh.vertices[follower].coord.get(axis, 0) -
                            mesh.vertices[seam].coord.get(axis, 0),
                        offsetBefore[axis], 1e-12);
        }
    }
}

// ---------------------------------------------------------------------------
// The energy difference
// ---------------------------------------------------------------------------

/**
 * @brief The local difference is the change of the whole-mesh energy.
 *
 * On a roughened fluid sheet with every soft wall live, for a vertex in the
 * middle of the tile and for one on the seam whose copies move with it, in
 * both periodic modes. The trial must also leave the mesh exactly as it was.
 */
TEST(VertexMoveTest, TheLocalEnergyDifferenceIsTheChangeOfTheTotal)
{
    for (BoundaryType boundary : kPeriodicModes)
    {
        SCOPED_TRACE(mode_name(boundary));
        Param param;
        configure_sheet(param, boundary);
        configure_fluid_walls(param);
        // Narrow walls, so the roughened sheet is leaning on all three.
        param.edgeTetherMinRatio = 0.95;
        param.edgeTetherMaxRatio = 1.05;
        param.triangleShapeMinAltitudeRatio = 0.85;
        param.creaseWallAngle = 5.0;
        DynamicMesh mesh(param);
        build_sheet(mesh);
        roughen(mesh, 0.5, 0.8, 11u);
        ASSERT_GT(total_regularization(mesh), 1.0);

        int seam = -1;
        int interior = -1;
        pick_seam_and_interior(mesh, seam, interior);
        ASSERT_GE(seam, 0);
        ASSERT_GE(interior, 0);

        const double moves[3][3] = {{0.20, -0.15, 0.25}, {-0.30, 0.10, -0.20}, {0.0, 0.0, 0.4}};
        for (int vertex : {interior, seam})
        {
            for (const double *d : {moves[0], moves[1], moves[2]})
            {
                const std::vector<double> coordsBefore = all_coordinates(mesh);
                const double energyBefore = total_energy(mesh);

                const VertexMoveDelta delta = mesh.evaluate_vertex_move(vertex, d);
                EXPECT_FALSE(delta.wallViolated);
                EXPECT_EQ(all_coordinates(mesh), coordsBefore)
                    << "the trial did not put the mesh back";

                mesh.displace_vertex_with_followers(vertex, d);
                const double energyAfter = total_energy(mesh);

                EXPECT_NEAR(delta.energy, energyAfter - energyBefore, 1e-7)
                    << "vertex " << vertex << (vertex == seam ? " (seam)" : " (interior)");
                EXPECT_NE(delta.regularization, 0.0);
                EXPECT_NE(delta.areaConstraint, 0.0);
            }
        }
    }
}

/**
 * @brief In solid mode the move samples bending and area, as the Brownian
 * step does.
 *
 * The reference-length regularization is evaluated by the force pass and
 * left out of the Brownian drive. A vertex move that differenced it would
 * give the two kernels different Hamiltonians.
 */
TEST(VertexMoveTest, ASolidSheetDifferencesBendingAndAreaOnly)
{
    Param param;
    configure_sheet(param, BoundaryType::Periodic);
    DynamicMesh mesh(param);
    build_sheet(mesh);
    roughen(mesh, 0.0, 0.8, 23u);
    // The reference lengths are the flat sheet's, so the term is live.
    ASSERT_GT(total_regularization(mesh), 1.0);

    int seam = -1;
    int interior = -1;
    pick_seam_and_interior(mesh, seam, interior);
    const double d[3] = {0.0, 0.0, 0.35};
    for (int vertex : {interior, seam})
    {
        const double before = bending_and_area_energy(mesh);
        const VertexMoveDelta delta = mesh.evaluate_vertex_move(vertex, d);
        mesh.displace_vertex_with_followers(vertex, d);
        const double after = bending_and_area_energy(mesh);
        EXPECT_EQ(delta.regularization, 0.0);
        EXPECT_NEAR(delta.energy, after - before, 1e-7);
    }
}

// ---------------------------------------------------------------------------
// The walls
// ---------------------------------------------------------------------------

/**
 * @brief A hard wall sits exactly where the soft one starts charging.
 *
 * Over sheets roughened from barely to badly: the soft energy is zero if and
 * only if no hard wall is broken. Each wall is tested alone, so that one
 * cannot hide another's mistake.
 */
TEST(VertexMoveTest, AHardWallIsBrokenExactlyWhereTheSoftOneCharges)
{
    const char *walls[3] = {"tether", "altitude", "crease"};
    for (int wall = 0; wall < 3; wall++)
    {
        SCOPED_TRACE(walls[wall]);
        int nCharged = 0;
        int nFree = 0;
        for (int level = 0; level < 8; level++)
        {
            Param param;
            configure_sheet(param, BoundaryType::Mixed);
            configure_fluid_walls(param);
            // One wall narrow enough to be reached, the other two out of
            // reach or off.
            param.edgeTetherMinRatio = (wall == 0) ? 0.9 : 0.01;
            param.edgeTetherMaxRatio = (wall == 0) ? 1.1 : 100.0;
            param.triangleShapeEnabled = (wall == 1);
            param.triangleShapeMinAltitudeRatio = 0.8;
            param.creaseWallEnabled = (wall == 2);
            param.creaseWallAngle = 8.0;
            DynamicMesh mesh(param);
            build_sheet(mesh);
            roughen(mesh, 0.08 * level, 0.08 * level, 100u + level);

            const double softEnergy = total_regularization(mesh);
            const WallViolations broken = mesh.count_wall_violations();
            EXPECT_EQ(softEnergy > 0.0, broken.total() > 0)
                << "level " << level << ": soft energy " << softEnergy << ", " << broken.total()
                << " walls broken";

            std::vector<int> everyFace(mesh.faces.size());
            for (int f = 0; f < static_cast<int>(mesh.faces.size()); f++)
            {
                everyFace[f] = f;
            }
            EXPECT_EQ(mesh.faces_satisfy_walls(everyFace), broken.total() == 0);

            (softEnergy > 0.0 ? nCharged : nFree)++;
        }
        // The sweep over levels has to cross the wall, or it tested nothing.
        EXPECT_GT(nCharged, 0);
        EXPECT_GT(nFree, 0);
    }
}

/**
 * @brief Under hard walls the fluid term carries no energy and no force, and
 * a move through a wall is refused before its energy is evaluated.
 */
TEST(VertexMoveTest, AHardWallCarriesNoEnergyAndRefusesTheMoveThroughIt)
{
    Param param;
    configure_sheet(param, BoundaryType::Mixed);
    configure_fluid_walls(param);
    configure_hard_walls(param);
    DynamicMesh mesh(param);
    build_sheet(mesh);
    roughen(mesh, 0.3, 0.3, 5u);
    ASSERT_EQ(mesh.count_wall_violations().total(), 0);

    EXPECT_EQ(total_regularization(mesh), 0.0);
    for (const Vertex &vertex : mesh.vertices)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            EXPECT_EQ(vertex.force.forceRegularization.get(axis, 0), 0.0);
        }
    }

    int seam = -1;
    int interior = -1;
    pick_seam_and_interior(mesh, seam, interior);
    const int neighbour = mesh.vertices[interior].adjacentVertices.front();

    // Straight at a neighbour, far enough to bring the edge under 3.5 nm.
    double toward[3];
    double length = 0.0;
    for (int axis = 0; axis < 3; axis++)
    {
        toward[axis] = mesh.vertices[neighbour].coord.get(axis, 0) -
                       mesh.vertices[interior].coord.get(axis, 0);
        length += toward[axis] * toward[axis];
    }
    length = std::sqrt(length);
    const double reach = length - 0.6 * mesh.param.lFace;
    for (int axis = 0; axis < 3; axis++)
    {
        toward[axis] *= reach / length;
    }

    const std::vector<double> coordsBefore = all_coordinates(mesh);
    const VertexMoveDelta refused = mesh.evaluate_vertex_move(interior, toward);
    EXPECT_TRUE(refused.wallViolated);
    EXPECT_TRUE(std::isinf(refused.energy));
    EXPECT_GT(refused.energy, 0.0);
    EXPECT_EQ(all_coordinates(mesh), coordsBefore);

    // A small move inside the walls is an ordinary trial with a finite cost,
    // and that cost is the bending and area change alone.
    const double small[3] = {0.05, 0.05, 0.05};
    const double before = total_energy(mesh);
    const VertexMoveDelta allowed = mesh.evaluate_vertex_move(interior, small);
    EXPECT_FALSE(allowed.wallViolated);
    EXPECT_EQ(allowed.regularization, 0.0);
    mesh.displace_vertex_with_followers(interior, small);
    EXPECT_NEAR(allowed.energy, total_energy(mesh) - before, 1e-7);
}

/**
 * @brief A flip that would leave an edge outside the tether is rejected, and
 * counted as an attempt.
 *
 * On the flat lattice every flip turns an edge of length l into one of
 * sqrt(3) l. With the upper wall below that no flip can be accepted; with it
 * above, the flip costs nothing and all of them are.
 */
TEST(VertexMoveTest, AFlipThroughAHardWallIsRejected)
{
    for (double maxRatio : {1.5, 1.75})
    {
        SCOPED_TRACE("edgeTetherMaxRatio = " + std::to_string(maxRatio));
        Param param;
        configure_sheet(param, BoundaryType::Mixed);
        configure_fluid_walls(param);
        configure_hard_walls(param);
        param.edgeTetherMaxRatio = maxRatio;
        param.edgeFlipEnabled = true;
        param.edgeFlipMinValence = 4;
        param.edgeFlipMaxValence = 8;
        param.edgeFlipAttemptRate = 50.0; // a few attempts per sweep on this sheet
        param.uSurf = 0.0;                // so that a flip inside the walls is free
        DynamicMesh mesh(param);
        build_sheet(mesh);

        EdgeFlipSweepStats totals;
        for (int sweep = 0; sweep < 20; sweep++)
        {
            const EdgeFlipSweepStats stats = mesh.edge_flip_sweep(sweep);
            totals.attempted += stats.attempted;
            totals.accepted += stats.accepted;
            EXPECT_EQ(mesh.count_wall_violations().total(), 0);
            if (stats.accepted > 0)
            {
                mesh.Compute_Energy_And_Force();
            }
        }
        ASSERT_GT(totals.attempted, 0);
        if (maxRatio < std::sqrt(3.0))
        {
            EXPECT_EQ(totals.accepted, 0);
        }
        else
        {
            EXPECT_GT(totals.accepted, 0);
        }
    }
}

/**
 * @brief Sweeps of vertex moves and flips never leave the admissible set.
 *
 * The invariant the whole scheme rests on: start inside the walls, and no
 * sequence of accepted moves gets out. Run long enough that the walls are
 * actually reached -- wall rejections must occur -- and with flips on, so
 * that the two moves are tested against each other's leftovers.
 */
TEST(VertexMoveTest, SweepsNeverLeaveTheAdmissibleSet)
{
    for (BoundaryType boundary : kPeriodicModes)
    {
        SCOPED_TRACE(mode_name(boundary));
        Param param;
        configure_sheet(param, boundary);
        configure_fluid_walls(param);
        configure_hard_walls(param);
        param.vertexMoveStepRatio = 0.2;
        param.edgeFlipEnabled = true;
        param.edgeFlipMinValence = 4;
        param.edgeFlipMaxValence = 8;
        param.edgeFlipAttemptRate = 100.0;
        DynamicMesh mesh(param);
        build_sheet(mesh);
        ASSERT_EQ(mesh.count_wall_violations().total(), 0);

        int seam = -1;
        int interior = -1;
        pick_seam_and_interior(mesh, seam, interior);
        const int follower = mesh.vertexMoveFollowers[seam].front();
        double offset[3];
        for (int axis = 0; axis < 3; axis++)
        {
            offset[axis] = mesh.vertices[follower].coord.get(axis, 0) -
                           mesh.vertices[seam].coord.get(axis, 0);
        }

        VertexMoveSweepStats totals;
        int flipsAccepted = 0;
        for (int sweep = 0; sweep < 12; sweep++)
        {
            const EdgeFlipSweepStats flips = mesh.edge_flip_sweep(sweep);
            flipsAccepted += flips.accepted;
            totals += mesh.vertex_move_sweep(sweep);
            mesh.Compute_Energy_And_Force();

            const WallViolations broken = mesh.count_wall_violations();
            ASSERT_EQ(broken.total(), 0)
                << "after sweep " << sweep << ": " << broken.tether << " tether, "
                << broken.altitude << " altitude, " << broken.crease << " crease";
            ASSERT_TRUE(std::isfinite(mesh.param.energy.energyTotal));
        }

        EXPECT_GT(totals.accepted, 0);
        EXPECT_GT(totals.wallRejected, 0) << "the walls were never reached";
        EXPECT_LT(totals.accepted, totals.attempted);
        EXPECT_GT(flipsAccepted, 0);

        // The copies are still exactly one period from their source.
        for (int axis = 0; axis < 3; axis++)
        {
            EXPECT_NEAR(mesh.vertices[follower].coord.get(axis, 0) -
                            mesh.vertices[seam].coord.get(axis, 0),
                        offset[axis], 1e-9);
        }

        // The running totals the trials difference against did not drift
        // from what a whole-mesh pass computes.
        const double runningArea = mesh.param.area;
        mesh.Compute_Energy_And_Force();
        EXPECT_NEAR(runningArea, mesh.param.area, 1e-7);
    }
}

// ---------------------------------------------------------------------------
// The stationary distribution
// ---------------------------------------------------------------------------

/**
 * @brief One vertex moving along z samples its Boltzmann distribution.
 *
 * With a single candidate and the move confined to z, the chain's stationary
 * density is exp(-E(z)/kT) / Z, a one-dimensional integral. The energy along
 * that line comes from the trial itself (pinned to the global energy above)
 * and the moments by quadrature; the sweep has to reproduce the variance.
 * Run at two temperatures, so that agreement cannot come from the step size
 * happening to match the width.
 */
TEST(VertexMoveTest, OneVertexSamplesItsBoltzmannDistribution)
{
    for (double kT : {4.17, 16.68})
    {
        SCOPED_TRACE("kT = " + std::to_string(kT));
        Param param;
        configure_sheet(param, BoundaryType::Periodic);
        param.KBT = kT;
        param.inPlaneDynamicsEnabled = false; // z only
        param.vertexMoveEnabled = true;
        param.vertexMoveStepRatio = 0.08;
        param.vertexMoveAttemptsPerVertex = 10.0;
        DynamicMesh mesh(param);
        build_sheet(mesh);
        roughen(mesh, 0.0, 0.3, 77u);

        int seam = -1;
        int interior = -1;
        pick_seam_and_interior(mesh, seam, interior);
        mesh.vertexMoveCandidates.assign(1, interior);
        const double z0 = mesh.vertices[interior].coord.get(2, 0);

        // exp(-E(z)/kT) on a grid wide enough to hold all of its weight.
        const int nGrid = 1601;
        const double halfWidth = 8.0;
        double weightSum = 0.0;
        double meanSum = 0.0;
        double squareSum = 0.0;
        for (int i = 0; i < nGrid; i++)
        {
            const double dz = -halfWidth + 2.0 * halfWidth * i / (nGrid - 1);
            const double d[3] = {0.0, 0.0, dz};
            const double weight = std::exp(-mesh.evaluate_vertex_move(interior, d).energy / kT);
            weightSum += weight;
            meanSum += weight * dz;
            squareSum += weight * dz * dz;
        }
        const double expectedMean = meanSum / weightSum;
        const double expectedVariance = squareSum / weightSum - expectedMean * expectedMean;
        // The grid has to contain the distribution, or the quadrature is not
        // the answer.
        ASSERT_LT(std::sqrt(expectedVariance), halfWidth / 6.0);

        const int nBurn = 100;
        const int nSweeps = 2000;
        double sum = 0.0;
        double sumSquares = 0.0;
        VertexMoveSweepStats totals;
        for (int sweep = 0; sweep < nBurn + nSweeps; sweep++)
        {
            totals += mesh.vertex_move_sweep(sweep);
            if (sweep < nBurn)
            {
                continue;
            }
            const double dz = mesh.vertices[interior].coord.get(2, 0) - z0;
            sum += dz;
            sumSquares += dz * dz;
        }
        const double mean = sum / nSweeps;
        const double variance = sumSquares / nSweeps - mean * mean;

        std::cout << "[VertexMoveTest] kT = " << kT << ": variance sampled " << variance
                  << " nm^2, Boltzmann " << expectedVariance << " nm^2 (ratio "
                  << variance / expectedVariance << "), acceptance " << totals.acceptance()
                  << std::endl;
        EXPECT_GT(totals.acceptance(), 0.2);
        EXPECT_LT(totals.acceptance(), 0.95);
        // 2000 samples ten attempts apart: a few per cent on the variance.
        EXPECT_NEAR(variance / expectedVariance, 1.0, 0.12)
            << "sampled " << variance << ", Boltzmann " << expectedVariance;
        EXPECT_NEAR(mean, expectedMean, 0.12 * std::sqrt(expectedVariance));
    }
}

/// The sweep is a function of the seed and the iteration, and of nothing else.
TEST(VertexMoveTest, TheSweepIsReproducible)
{
    std::vector<double> first;
    std::vector<double> second;
    std::vector<double> otherIteration;
    for (int run = 0; run < 3; run++)
    {
        Param param;
        configure_sheet(param, BoundaryType::Mixed);
        configure_fluid_walls(param);
        param.vertexMoveEnabled = true;
        DynamicMesh mesh(param);
        build_sheet(mesh);
        mesh.vertex_move_sweep(run == 2 ? 8 : 7);
        (run == 0 ? first : run == 1 ? second : otherIteration) = all_coordinates(mesh);
    }
    EXPECT_EQ(first, second);
    EXPECT_NE(first, otherIteration);
}

// ---------------------------------------------------------------------------
// The configurations that are not samplers
// ---------------------------------------------------------------------------

TEST(VertexMoveTest, SetupRefusesTheCombinationsThatDoNotSample)
{
    {
        // A hard wall has no force for the Brownian step to feel.
        Param param;
        configure_sheet(param, BoundaryType::Mixed);
        configure_fluid_walls(param);
        param.meshWallType = "hard";
        param.vertexMoveEnabled = true;
        DynamicMesh mesh(param);
        EXPECT_THROW(mesh.setup_flat(), std::runtime_error);
    }
    {
        // Nothing would move a vertex.
        Param param;
        configure_sheet(param, BoundaryType::Mixed);
        param.brownianStepEnabled = false;
        DynamicMesh mesh(param);
        EXPECT_THROW(mesh.setup_flat(), std::runtime_error);
    }
    {
        // Hard walls are the flat tether's range; without the tether there is none.
        Param param;
        configure_sheet(param, BoundaryType::Mixed);
        param.meshWallType = "hard";
        param.brownianStepEnabled = false;
        param.vertexMoveEnabled = true;
        DynamicMesh mesh(param);
        EXPECT_THROW(mesh.setup_flat(), std::runtime_error);
    }
    {
        // An energy term the trial does not difference.
        Param param;
        configure_sheet(param, BoundaryType::Mixed);
        param.vertexMoveEnabled = true;
        param.beadEnabled = true;
        DynamicMesh mesh(param);
        EXPECT_THROW(mesh.setup_flat(), std::runtime_error);
    }
    {
        Param param;
        configure_sheet(param, BoundaryType::Mixed);
        param.meshWallType = "firm";
        DynamicMesh mesh(param);
        EXPECT_THROW(mesh.setup_flat(), std::runtime_error);
    }
    {
        // The pure Monte Carlo fluid run itself is accepted.
        Param param;
        configure_sheet(param, BoundaryType::Mixed);
        configure_fluid_walls(param);
        param.meshWallType = "hard";
        param.brownianStepEnabled = false;
        param.vertexMoveEnabled = true;
        DynamicMesh mesh(param);
        EXPECT_NO_THROW(mesh.setup_flat());
    }
}
