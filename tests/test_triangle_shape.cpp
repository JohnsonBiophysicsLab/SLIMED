#include "test_triangle_shape.hpp"

/**
 * The tether bounds every edge of a fluid mesh and nothing else, and every
 * long fluid run ended in a folded sliver whose edges were all inside the
 * walls. This term bounds the quantity that actually vanishes in a sliver --
 * the altitude from each corner to its opposite edge -- and these are the
 * checks that it does so without doing anything else:
 *
 * - it is exactly zero on the lattice and charges a sliver what the formula
 *   says, so it is a wall and not a bias;
 * - its force is minus the gradient of its energy, checked by finite
 *   differences with the term active on every face;
 * - the flip trial and the force pass reach one and the same expression --
 *   the WP6 lesson, made a gate before the term is used anywhere;
 * - a flip of an equilateral rhombus, whose new triangles have altitude
 *   exactly half an edge, is free below the default floor and charged above
 *   it, which is what fixes the floor's ceiling;
 * - and in a short fluid run the smallest altitude on the sheet stays above
 *   the floor's shadow, where without the term it does not.
 *
 * @see docs/edge_flip_plan.md work package 7
 */

namespace
{

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
    param.randomSeed = 20260911u;
    param.timeStep = 1.0e-3;
    param.diffConst = 1.0;
    param.surfaceSolver = "iterative";
    param.edgeSpringEnabled = true;
    param.edgeSpringConstant = 83.4;
    param.triangleShapeEnabled = true;
    param.triangleShapeConstant = 83.4;
}

/// Only the shape term in the regularization slot: the tether's walls out of reach.
void tether_out_of_the_way(Param &param)
{
    param.edgeTetherMinRatio = 0.01;
    param.edgeTetherMaxRatio = 100.0;
}

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

/// The three altitudes of a face, computed here independently of the term.
std::array<double, 3> altitudes(const Mesh &mesh, int iFace)
{
    const std::vector<int> &c = mesh.faces[iFace].adjacentVertices;
    double p[3][3];
    for (int k = 0; k < 3; k++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            p[k][axis] = mesh.vertices[c[k]].coord.get(axis, 0);
        }
    }
    double d1[3], d2[3];
    for (int axis = 0; axis < 3; axis++)
    {
        d1[axis] = p[1][axis] - p[0][axis];
        d2[axis] = p[2][axis] - p[0][axis];
    }
    const double nx = d1[1] * d2[2] - d1[2] * d2[1];
    const double ny = d1[2] * d2[0] - d1[0] * d2[2];
    const double nz = d1[0] * d2[1] - d1[1] * d2[0];
    const double twoArea = std::sqrt(nx * nx + ny * ny + nz * nz);
    std::array<double, 3> h;
    for (int i = 0; i < 3; i++)
    {
        double l2 = 0.0;
        for (int axis = 0; axis < 3; axis++)
        {
            const double d = p[(i + 1) % 3][axis] - p[(i + 2) % 3][axis];
            l2 += d * d;
        }
        h[i] = twoArea / std::sqrt(l2);
    }
    return h;
}

double smallest_interior_altitude(const Mesh &mesh)
{
    double smallest = 1e300;
    for (int iFace = 0; iFace < static_cast<int>(mesh.faces.size()); iFace++)
    {
        const Face &face = mesh.faces[iFace];
        if (face.isGhost)
        {
            continue;
        }
        bool touchesCopy = false;
        for (int v : face.adjacentVertices)
        {
            touchesCopy = touchesCopy || mesh.vertices[v].isGhost ||
                          (!mesh.flipFrozenVertex.empty() && mesh.flipFrozenVertex[v]);
        }
        if (touchesCopy)
        {
            continue;
        }
        const std::array<double, 3> h = altitudes(mesh, iFace);
        smallest = std::min(smallest, *std::min_element(h.begin(), h.end()));
    }
    return smallest;
}

int an_interior_face(const Mesh &mesh)
{
    for (int iFace = 0; iFace < static_cast<int>(mesh.faces.size()); iFace++)
    {
        const Face &face = mesh.faces[iFace];
        if (face.isGhost)
        {
            continue;
        }
        bool deep = true;
        for (int v : face.adjacentVertices)
        {
            deep = deep && !mesh.vertices[v].isGhost &&
                   (mesh.flipFrozenVertex.empty() || !mesh.flipFrozenVertex[v]);
        }
        // Near the centre of a sheet that sits about the origin, well away
        // from the duplicate ring.
        if (deep && std::abs(mesh.vertices[face.adjacentVertices[0]].coord.get(0, 0)) < 8.0 &&
            std::abs(mesh.vertices[face.adjacentVertices[0]].coord.get(1, 0)) < 8.0)
        {
            return iFace;
        }
    }
    return -1;
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

/// The driver's loop body, for a short fluid run.
double run_steps(DynamicMesh &mesh, DynamicModel &model, Record &record, int nSteps,
                 double *smallestAltitudeSeen)
{
    double smallest = 1e300;
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
            if (mesh.edge_flip_sweep(model.iteration).accepted > 0)
            {
                mesh.Compute_Energy_And_Force();
            }
        }
        smallest = std::min(smallest, smallest_interior_altitude(mesh));
    }
    if (smallestAltitudeSeen != nullptr)
    {
        *smallestAltitudeSeen = smallest;
    }
    return mesh.param.energy.energyTotal;
}

} // namespace

// ---------------------------------------------------------------------------
// A wall, not a bias
// ---------------------------------------------------------------------------

TEST(TriangleShapeTest, TheLatticeCarriesNoShapeEnergy)
{
    Param param;
    configure_sheet(param);
    tether_out_of_the_way(param);
    DynamicMesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_sheet(mesh));

    // Every altitude on the lattice is 0.866 lFace; the floor is 0.4 lFace.
    EXPECT_EQ(total_regularization(mesh), 0.0)
        << "the lattice, whose triangles are equilateral, carries shape energy";
    for (int iFace = 0; iFace < static_cast<int>(mesh.faces.size()); iFace++)
    {
        ASSERT_EQ(mesh.face_shape_energy(iFace), 0.0) << "face " << iFace;
    }
}

TEST(TriangleShapeTest, ASliverIsChargedWhatTheFormulaSays)
{
    Param param;
    configure_sheet(param);
    tether_out_of_the_way(param);
    DynamicMesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_sheet(mesh));

    const int iFace = an_interior_face(mesh);
    ASSERT_GE(iFace, 0);
    const std::vector<int> c = mesh.faces[iFace].adjacentVertices;

    // Slide corner 0 toward the midpoint of its opposite edge until that
    // altitude is 0.5 nm: a sliver of the kind the fluid runs were full of.
    Matrix &p0 = mesh.vertices[c[0]].coord;
    const Matrix &p1 = mesh.vertices[c[1]].coord;
    const Matrix &p2 = mesh.vertices[c[2]].coord;
    double mid[3], toward[3];
    double norm = 0.0;
    for (int axis = 0; axis < 3; axis++)
    {
        mid[axis] = 0.5 * (p1.get(axis, 0) + p2.get(axis, 0));
        toward[axis] = p0.get(axis, 0) - mid[axis];
        norm += toward[axis] * toward[axis];
    }
    norm = std::sqrt(norm);
    for (int axis = 0; axis < 3; axis++)
    {
        p0.set(axis, 0, mid[axis] + 0.5 * toward[axis] / norm);
    }

    const std::array<double, 3> h = altitudes(mesh, iFace);
    EXPECT_NEAR(h[0], 0.5, 1e-9) << "the sliver was not built as intended";

    const double floor = param.triangleShapeMinAltitudeRatio * param.lFace;
    double expected = 0.0;
    for (double hi : h)
    {
        if (hi < floor)
        {
            expected += 0.5 * param.triangleShapeConstant * (floor - hi) * (floor - hi);
        }
    }
    ASSERT_GT(expected, 0.0);
    EXPECT_NEAR(mesh.face_shape_energy(iFace), expected, 1e-9 * expected);

    // And through the force pass, which sums the term over faces: every
    // other face this corner belongs to was squeezed too, so compare the
    // whole slot against the same hand formula summed over all faces.
    mesh.energy_force_edge_spring();
    mesh.energy_force_triangle_shape();
    double expectedTotal = 0.0;
    for (int f = 0; f < static_cast<int>(mesh.faces.size()); f++)
    {
        if (mesh.faces[f].isGhost)
        {
            continue;
        }
        for (double hi : altitudes(mesh, f))
        {
            if (hi < floor)
            {
                expectedTotal += 0.5 * param.triangleShapeConstant * (floor - hi) * (floor - hi);
            }
        }
    }
    EXPECT_NEAR(total_regularization(mesh), expectedTotal, 1e-9 * expectedTotal);
}

// ---------------------------------------------------------------------------
// The force is the gradient
// ---------------------------------------------------------------------------

TEST(TriangleShapeTest, TheForceIsMinusTheGradientOfTheEnergy)
{
    Param param;
    configure_sheet(param);
    tether_out_of_the_way(param);
    param.kCurv = 0.0; // the shape term alone, so nothing else can compensate
    param.uSurf = 0.0;
    param.isGlobalConstraint = false;
    // A floor above the lattice's own altitude puts every face on a wall, so
    // the gradient is nonzero everywhere and a finite difference tests it
    // rather than passing against zero.
    param.triangleShapeMinAltitudeRatio = 0.95;

    DynamicMesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_sheet(mesh));

    // Some disorder, so the three altitudes of a face differ and the edge
    // and area gradients are both exercised.
    std::uint64_t state = 0x9E3779B97F4A7C15ULL;
    const auto next = [&state]() {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>((state >> 11) * (1.0 / 9007199254740992.0)) - 0.5;
    };
    for (Vertex &vertex : mesh.vertices)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            vertex.coord.set(axis, 0, vertex.coord.get(axis, 0) + 0.8 * next());
        }
    }
    mesh.energy_force_edge_spring();
    mesh.energy_force_triangle_shape();
    ASSERT_GT(total_regularization(mesh), 0.0) << "no face is on the wall";

    const double h = 1e-6;
    int checked = 0;
    for (int iVertex = 30; iVertex < static_cast<int>(mesh.vertices.size()) && checked < 15;
         iVertex += 7)
    {
        if (mesh.vertices[iVertex].isGhost)
        {
            continue;
        }
        for (int axis = 0; axis < 3; axis++)
        {
            const double original = mesh.vertices[iVertex].coord.get(axis, 0);

            mesh.vertices[iVertex].coord.set(axis, 0, original + h);
            mesh.energy_force_edge_spring();
            mesh.energy_force_triangle_shape();
            const double plus = total_regularization(mesh);

            mesh.vertices[iVertex].coord.set(axis, 0, original - h);
            mesh.energy_force_edge_spring();
            mesh.energy_force_triangle_shape();
            const double minus = total_regularization(mesh);

            mesh.vertices[iVertex].coord.set(axis, 0, original);
            mesh.energy_force_edge_spring();
            mesh.energy_force_triangle_shape();
            const double analytic = mesh.vertices[iVertex].force.forceRegularization.get(axis, 0);

            const double numeric = -(plus - minus) / (2.0 * h);
            EXPECT_NEAR(analytic, numeric, 1e-4 * std::max(1.0, std::abs(numeric)))
                << "vertex " << iVertex << " axis " << axis;
            checked++;
        }
    }
    EXPECT_GE(checked, 12);
}

// ---------------------------------------------------------------------------
// One Hamiltonian
// ---------------------------------------------------------------------------

TEST(TriangleShapeTest, TheTrialAndTheDynamicsShareOneShapeTerm)
{
    Param param;
    configure_sheet(param);
    param.triangleShapeMinAltitudeRatio = 0.95; // active on every face
    DynamicMesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_sheet(mesh));

    const double fromForcePass = total_regularization(mesh);
    ASSERT_GT(fromForcePass, 0.0);

    std::vector<int> everyFace(mesh.faces.size());
    for (int iFace = 0; iFace < static_cast<int>(mesh.faces.size()); iFace++)
    {
        everyFace[iFace] = iFace;
    }
    const FaceSubsetEnergy subset = mesh.evaluate_face_subset(everyFace);
    EXPECT_NEAR(subset.regularization, fromForcePass, 1e-9 * fromForcePass)
        << "the trial's evaluator returned " << subset.regularization
        << " where the force pass computes " << fromForcePass;

    // And end to end: a trial's dE is the change of the whole-mesh energy.
    const int iEdge = first_flippable_edge(mesh);
    ASSERT_GE(iEdge, 0);
    const double before = total_energy(mesh);
    EdgeFlipDelta delta;
    ASSERT_TRUE(mesh.evaluate_edge_flip(iEdge, delta));
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
// The floor's ceiling
// ---------------------------------------------------------------------------

TEST(TriangleShapeTest, AFlipOfAnEquilateralRhombusIsFreeBelowHalfAndChargedAbove)
{
    const auto shape_cost_of_one_flip = [](double ratio) {
        Param param;
        configure_sheet(param);
        tether_out_of_the_way(param);
        param.kCurv = 0.0;
        param.uSurf = 0.0;
        param.isGlobalConstraint = false;
        param.triangleShapeMinAltitudeRatio = ratio;
        DynamicMesh mesh(param);
        build_sheet(mesh);
        const int iEdge = first_flippable_edge(mesh);
        EXPECT_GE(iEdge, 0);
        EdgeFlipDelta delta;
        EXPECT_TRUE(mesh.evaluate_edge_flip(iEdge, delta));
        return delta.energy;
    };

    // The two new triangles have edges (l, l, sqrt(3) l) and altitude l/2
    // against the long edge: 2.5 nm here.
    EXPECT_NEAR(shape_cost_of_one_flip(0.4), 0.0, 1e-12)
        << "the default floor charges the flip of an equilateral rhombus";

    // Above half, each new triangle owes one wall of (0.6 l - l/2) = 0.5 nm.
    const double owed = 2.0 * 0.5 * 83.4 * 0.5 * 0.5;
    EXPECT_NEAR(shape_cost_of_one_flip(0.6), owed, 1e-6 * owed);
}

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

TEST(TriangleShapeTest, ShapeWithoutTheTetherIsRefusedAtSetup)
{
    Param param;
    configure_sheet(param);
    param.edgeSpringEnabled = false;
    DynamicMesh mesh(param);
    EXPECT_THROW(mesh.setup_flat(), std::runtime_error);
}

TEST(TriangleShapeTest, TheParametersRoundTrip)
{
    Param param;
    EXPECT_FALSE(param.triangleShapeEnabled);
    EXPECT_DOUBLE_EQ(param.triangleShapeMinAltitudeRatio, 0.4);
    EXPECT_TRUE(import_kv_string("triangleShapeEnabled", "true", param));
    EXPECT_TRUE(param.triangleShapeEnabled);
    EXPECT_TRUE(import_kv_string("triangleShapeMinAltitudeRatio", "0.35", param));
    EXPECT_DOUBLE_EQ(param.triangleShapeMinAltitudeRatio, 0.35);
    EXPECT_TRUE(import_kv_string("triangleShapeConstant", "40", param));
    EXPECT_DOUBLE_EQ(param.triangleShapeConstant, 40.0);
}

// ---------------------------------------------------------------------------
// What it is for
// ---------------------------------------------------------------------------

/**
 * @brief In a fluid run the smallest altitude stays above the floor's shadow.
 *
 * Below the floor the term is a harmonic well of width sqrt(kT/k) = 0.22 nm,
 * so an altitude 1.2 nm under a 2 nm floor is a five-sigma excursion. Without
 * the term the 100 nm sheet reaches 0.1 nm within a few thousand steps.
 */
TEST(TriangleShapeTest, SliversAreHeldAboveTheFloorInAFluidRun)
{
    const auto smallest_altitude_after = [](bool withShape) {
        Param param;
        configure_sheet(param);
        param.inPlaneDynamicsEnabled = true;
        param.edgeFlipEnabled = true;
        param.triangleShapeEnabled = withShape;
        param.maxIterations = 250;
        DynamicMesh mesh(param);
        build_sheet(mesh);
        set_rate_for_expected_attempts(mesh, 6.0);
        Record record(param.maxIterations + 1);
        DynamicModel model(mesh, record);
        double smallest = 0.0;
        run_steps(mesh, model, record, param.maxIterations, &smallest);
        return smallest;
    };

    const double without = smallest_altitude_after(false);
    const double with = smallest_altitude_after(true);
    std::cout << "[TriangleShapeTest] smallest interior altitude over 250 fluid steps: "
              << without << " nm without the term, " << with << " nm with it (floor 2 nm, "
                 "lattice 4.33)." << std::endl;

    EXPECT_GT(with, 0.8)
        << "with the shape term the mesh still reached an altitude of " << with << " nm";
}
