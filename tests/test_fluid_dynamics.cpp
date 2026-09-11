#include "test_fluid_dynamics.hpp"

/**
 * Three changes make the dynamics fluid rather than solid, and each fails in
 * its own way if it is wrong.
 *
 * The **edge spring** replaces a term that measured every face against the
 * configuration it started in. That memory is what a solid has and a fluid does
 * not, and an edge a flip has just created never had a reference length at all.
 * If the spring's force does not match its energy, the Metropolis sweep and the
 * Brownian step disagree about the Hamiltonian and neither samples anything.
 *
 * **In-plane motion** is half of what fluidity means. Flipping the
 * connectivity while the vertices stay pinned laterally is doing half the job.
 *
 * The **surface solver** converts between the control net the energy is a
 * function of and the limit surface the step displaces. The dense path built
 * that map with the valence-6 mask whatever the actual valence was, inverted it
 * once, and could not update it after a flip. Getting the force direction wrong
 * here is not a small bias: a non-reciprocal mobility is not the gradient flow
 * of any potential and can do net work on the membrane.
 *
 * @see docs/edge_flip_plan.md section 3.7
 */

namespace
{
struct MeshFixture
{
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;
};

MeshFixture build_icosphere(int level, double radius)
{
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    std::vector<std::array<double, 3>> points = {
        {{-1, phi, 0}}, {{1, phi, 0}},  {{-1, -phi, 0}}, {{1, -phi, 0}},
        {{0, -1, phi}}, {{0, 1, phi}},  {{0, -1, -phi}}, {{0, 1, -phi}},
        {{phi, 0, -1}}, {{phi, 0, 1}},  {{-phi, 0, -1}}, {{-phi, 0, 1}}};
    std::vector<std::array<int, 3>> triangles = {
        {{0, 11, 5}}, {{0, 5, 1}},  {{0, 1, 7}},   {{0, 7, 10}}, {{0, 10, 11}},
        {{1, 5, 9}},  {{5, 11, 4}}, {{11, 10, 2}}, {{10, 7, 6}}, {{7, 1, 8}},
        {{3, 9, 4}},  {{3, 4, 2}},  {{3, 2, 6}},   {{3, 6, 8}},  {{3, 8, 9}},
        {{4, 9, 5}},  {{2, 4, 11}}, {{6, 2, 10}},  {{8, 6, 7}},  {{9, 8, 1}}};
    const auto onSphere = [](std::array<double, 3> p) {
        const double n = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        return std::array<double, 3>{{p[0] / n, p[1] / n, p[2] / n}};
    };
    for (auto &p : points)
    {
        p = onSphere(p);
    }
    for (int step = 0; step < level; step++)
    {
        std::map<std::pair<int, int>, int> midpoint;
        const auto edgePoint = [&](int a, int b) {
            const std::pair<int, int> key = (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
            const auto found = midpoint.find(key);
            if (found != midpoint.end())
            {
                return found->second;
            }
            const int index = static_cast<int>(points.size());
            points.push_back(onSphere({{0.5 * (points[a][0] + points[b][0]),
                                        0.5 * (points[a][1] + points[b][1]),
                                        0.5 * (points[a][2] + points[b][2])}}));
            midpoint.emplace(key, index);
            return index;
        };
        std::vector<std::array<int, 3>> refined;
        for (const std::array<int, 3> &t : triangles)
        {
            const int ab = edgePoint(t[0], t[1]);
            const int bc = edgePoint(t[1], t[2]);
            const int ca = edgePoint(t[2], t[0]);
            refined.push_back({{t[0], ab, ca}});
            refined.push_back({{ab, t[1], bc}});
            refined.push_back({{ca, bc, t[2]}});
            refined.push_back({{ab, bc, ca}});
        }
        triangles = refined;
    }
    MeshFixture fixture;
    for (const std::array<double, 3> &p : points)
    {
        fixture.vertices.push_back({radius * p[0], radius * p[1], radius * p[2]});
    }
    for (const std::array<int, 3> &t : triangles)
    {
        fixture.faces.push_back({t[0], t[1], t[2]});
    }
    return fixture;
}

void configure(Param &param)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.kCurv = 83.4;
    param.KBT = 4.17;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
    param.isGlobalConstraint = true;
    param.uSurf = 0.0;
    param.uVol = 0.0;
    param.area0 = 1.0;
    param.vol0 = 1.0;
    param.lFace = 5.0;
    param.randomSeed = 5150u;
}

double total_energy(Mesh &mesh)
{
    mesh.calculate_element_area_volume();
    mesh.Compute_Energy_And_Force();
    return mesh.param.energy.energyTotal;
}
} // namespace

// ---------------------------------------------------------------------------
// The edge spring
// ---------------------------------------------------------------------------

/**
 * @brief The spring's force is minus the gradient of its energy.
 *
 * The single test that matters for a new energy term. If it fails, the flip
 * sweep and the Brownian step are working from different Hamiltonians and
 * nothing either of them samples means anything.
 */
TEST(FluidDynamicsTest, TheEdgeSpringForceIsMinusTheGradientOfItsEnergy)
{
    Param param;
    configure(param);
    param.edgeSpringEnabled = true;
    // Explicitly the harmonic shape. The default is now the flat-bottomed
    // tether, whose energy is identically zero for a mesh whose edges are all
    // inside the allowed range -- which this icosphere's are, so these tests
    // would pass against nothing. See WP6.
    param.edgeTetherShape = "harmonic";
    param.edgeSpringConstant = 40.0;
    param.edgeSpringRestLength = 6.0;
    param.kCurv = 0.0; // the spring alone, so nothing else can compensate

    Mesh mesh(param);
    const MeshFixture fixture = build_icosphere(1, 12.0);
    mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();

    const double h = 1e-6;
    for (int iVertex : {0, 7, 19, 33})
    {
        for (int axis = 0; axis < 3; axis++)
        {
            total_energy(mesh);
            const double analytic = mesh.vertices[iVertex].force.forceRegularization(axis, 0);

            const double original = mesh.vertices[iVertex].coord.get(axis, 0);
            mesh.vertices[iVertex].coord.set(axis, 0, original + h);
            const double plus = total_energy(mesh);
            mesh.vertices[iVertex].coord.set(axis, 0, original - h);
            const double minus = total_energy(mesh);
            mesh.vertices[iVertex].coord.set(axis, 0, original);

            const double numeric = -(plus - minus) / (2.0 * h);
            EXPECT_NEAR(analytic, numeric, 1e-5 * std::max(1.0, std::abs(numeric)))
                << "vertex " << iVertex << " axis " << axis;
        }
    }
}

/**
 * @brief The spring energy is a sum over the edges that exist now, and the
 * per-face shares add up to it.
 *
 * Which is the property that lets it survive a flip: the reference-length term
 * it replaces needs an edge to have existed in coordRef, and a flipped edge
 * never did.
 */
TEST(FluidDynamicsTest, TheEdgeSpringIsASumOverTheEdgesThatExistNow)
{
    Param param;
    configure(param);
    param.edgeSpringEnabled = true;
    // Explicitly the harmonic shape. The default is now the flat-bottomed
    // tether, whose energy is identically zero for a mesh whose edges are all
    // inside the allowed range -- which this icosphere's are, so these tests
    // would pass against nothing. See WP6.
    param.edgeTetherShape = "harmonic";
    param.edgeSpringConstant = 40.0;
    param.edgeSpringRestLength = 6.0;

    Mesh mesh(param);
    const MeshFixture fixture = build_icosphere(1, 12.0);
    mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    total_energy(mesh);

    const auto direct_sum = [&]() {
        double sum = 0.0;
        for (const MeshEdge &edge : mesh.edges)
        {
            double along[3];
            for (int axis = 0; axis < 3; axis++)
            {
                along[axis] = mesh.vertices[edge.v[0]].coord.get(axis, 0) -
                              mesh.vertices[edge.v[1]].coord.get(axis, 0);
            }
            const double length =
                std::sqrt(along[0] * along[0] + along[1] * along[1] + along[2] * along[2]);
            const double extension = length - param.edgeSpringRestLength;
            sum += 0.5 * param.edgeSpringConstant * extension * extension;
        }
        return sum;
    };

    double fromFaces = 0.0;
    for (const Face &face : mesh.faces)
    {
        fromFaces += face.energy.energyRegularization;
    }
    EXPECT_NEAR(fromFaces, direct_sum(), 1e-9 * std::max(1.0, direct_sum()));

    // And it still adds up after the connectivity changes, which is the whole
    // point. The reference-length term cannot say that: a flipped edge has no
    // reference length.
    int flipped = 0;
    for (const MeshEdge &edge : mesh.edges)
    {
        if (flipped >= 5)
        {
            break;
        }
        if (mesh.edge_flip_is_admissible(edge.index))
        {
            mesh.flip_edge(edge.index);
            flipped++;
        }
    }
    ASSERT_EQ(flipped, 5);
    total_energy(mesh);
    fromFaces = 0.0;
    for (const Face &face : mesh.faces)
    {
        fromFaces += face.energy.energyRegularization;
    }
    EXPECT_NEAR(fromFaces, direct_sum(), 1e-9 * std::max(1.0, direct_sum()));
}

TEST(FluidDynamicsTest, TheSpringAndTheReferenceLengthTermAreAlternatives)
{
    const auto regularization_energy = [](bool useSpring) {
        Param param;
        configure(param);
        param.edgeSpringEnabled = useSpring;
        param.edgeTetherShape = "harmonic";
        param.edgeSpringConstant = 40.0;
        param.edgeSpringRestLength = 6.0;
        Mesh mesh(param);
        const MeshFixture fixture = build_icosphere(1, 12.0);
        mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces);
        mesh.update_previous_coord_for_vertex();
        mesh.update_reference_coord_from_previous_coord();
        total_energy(mesh);
        double sum = 0.0;
        for (const Face &face : mesh.faces)
        {
            sum += face.energy.energyRegularization;
        }
        return sum;
    };

    // With coordRef equal to the current coordinates, the reference-length term
    // is exactly zero -- every edge is already at its remembered length. The
    // spring is not, because its rest length is a parameter rather than a
    // memory. That difference is the whole distinction between a solid's
    // restoring force and a fluid's.
    EXPECT_NEAR(regularization_energy(false), 0.0, 1e-12);
    EXPECT_GT(regularization_energy(true), 0.0);
}

// ---------------------------------------------------------------------------
// The surface solver
// ---------------------------------------------------------------------------

namespace
{
/// A flat sheet through the full DynamicMesh setup, which is what builds the
/// dense conversion matrices the iterative path is compared against.
void build_flat_dynamic_mesh(DynamicMesh &mesh)
{
    mesh.setup_flat();
    mesh.update_vertices_mat_with_vector();
}

void configure_flat(Param &param)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Periodic;
    param.sideX = 60.0;
    param.sideY = 60.0;
    param.lFace = 5.0;
    param.uVol = 0.0;
    param.kCurv = 83.4;
    param.KBT = 4.17;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
}
} // namespace

/**
 * @brief The sparse solver reproduces the dense one on the mesh the dense one
 * is right for.
 *
 * Everything the flat sheet's mask touches is at valence 6, so the two
 * representations describe the same matrix and must agree to the solver's
 * tolerance in both directions. This is what says the sparse assembly and the
 * conjugate gradient are right before either is used anywhere the dense path
 * cannot follow.
 */
TEST(FluidDynamicsTest, TheSparseSolverAgreesWithTheDenseOneOnARegularMesh)
{
    Param param;
    configure_flat(param);
    DynamicMesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_flat_dynamic_mesh(mesh));

    const int n = static_cast<int>(mesh.vertices.size());
    // Give the sheet some relief, so the comparison is not between two copies
    // of a plane.
    for (int v = 0; v < n; v++)
    {
        const double x = mesh.matMesh(v, 0);
        const double y = mesh.matMesh(v, 1);
        mesh.matMesh.set(v, 2, 0.7 * std::sin(0.21 * x) * std::cos(0.17 * y));
    }

    // Forwards.
    Matrix denseSurface = mat_calloc(n, 3);
    denseSurface = mesh.mesh2surface * mesh.matMesh;

    mesh.ensure_surface_solver();
    Matrix sparseSurface = mat_calloc(n, 3);
    mesh.surfaceSolver.mesh_to_surface(mesh.matMesh, sparseSurface);

    double worstForward = 0.0;
    for (int v = 0; v < n; v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            worstForward =
                std::max(worstForward, std::abs(denseSurface(v, axis) - sparseSurface(v, axis)));
        }
    }
    EXPECT_LT(worstForward, 1e-12) << "M C disagrees between the dense and sparse masks";

    // Backwards. The dense path inverts; the sparse one solves.
    Matrix denseControl = mat_calloc(n, 3);
    denseControl = mesh.surface2mesh * denseSurface;
    Matrix sparseControl = mat_calloc(n, 3);
    const int iterations = mesh.surfaceSolver.surface_to_mesh(sparseSurface, sparseControl);
    EXPECT_GT(iterations, 0);
    EXPECT_LT(mesh.surfaceSolver.lastResidual(), 1e-9)
        << "the solve did not reproduce the surface it was given";

    double worstBackward = 0.0;
    for (int v = 0; v < n; v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            worstBackward = std::max(worstBackward,
                                     std::abs(denseControl(v, axis) - sparseControl(v, axis)));
        }
    }
    // The dense inverse and the iterative solve answer the same question by
    // different routes, so they agree to the tolerance of the looser one.
    EXPECT_LT(worstBackward, 1e-7)
        << "solving M C = S disagrees with multiplying by the stored inverse";
}

/**
 * @brief The round trip is exact: solving after applying returns what went in.
 *
 * Independent of the dense path entirely, and the property the dynamics
 * actually relies on every step.
 */
TEST(FluidDynamicsTest, ApplyingTheMaskAndSolvingItBackIsTheIdentity)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    const MeshFixture fixture = build_icosphere(2, 12.0);
    mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces);

    // A fluid mesh: flips first, so the valences are mixed and the valence-6
    // mask would be the wrong matrix.
    int flipped = 0;
    for (const MeshEdge &edge : mesh.edges)
    {
        if (flipped >= 40)
        {
            break;
        }
        if (mesh.edge_flip_is_admissible(edge.index))
        {
            mesh.flip_edge(edge.index);
            flipped++;
        }
    }
    ASSERT_EQ(flipped, 40);

    int nNotSix = 0;
    for (const Vertex &vertex : mesh.vertices)
    {
        nNotSix += (vertex.adjacentVertices.size() != 6) ? 1 : 0;
    }
    ASSERT_GT(nNotSix, 20) << "the fixture should have mixed valences";

    slimed::SurfaceSolver solver;
    solver.build(mesh);
    ASSERT_EQ(solver.nFree(), solver.nVertices()) << "a closed mesh pins nothing";

    const int n = static_cast<int>(mesh.vertices.size());
    Matrix control = mat_calloc(n, 3);
    for (int v = 0; v < n; v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            control.set(v, axis, mesh.vertices[v].coord.get(axis, 0));
        }
    }

    Matrix surface = mat_calloc(n, 3);
    solver.mesh_to_surface(control, surface);
    Matrix recovered = mat_calloc(n, 3);
    solver.surface_to_mesh(surface, recovered);

    double worst = 0.0;
    for (int v = 0; v < n; v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            worst = std::max(worst, std::abs(control(v, axis) - recovered(v, axis)));
        }
    }
    EXPECT_LT(worst, 1e-8) << "M then M^-1 is not the identity; worst error " << worst;
}

/**
 * @brief The force map is the transpose, and on a mesh with mixed valences that
 * is a different answer from the inverse.
 *
 * The energy is a function of the control points; the step displaces the
 * surface. The chain rule gives F_S = M^-T F_C, and the tree has been using
 * M^-1 -- identical wherever the mask is symmetric, which is every interior
 * valence-6 row, and not identical once the valences are mixed. Getting it
 * wrong is not only a bias on the sampled distribution: a non-reciprocal
 * mobility is not the gradient flow of any potential and can do net work on
 * the membrane.
 */
TEST(FluidDynamicsTest, TheForceMapIsTheTransposeAndDiffersFromTheInverse)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    const MeshFixture fixture = build_icosphere(2, 12.0);
    mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces);
    int flipped = 0;
    for (const MeshEdge &edge : mesh.edges)
    {
        if (flipped >= 40)
        {
            break;
        }
        if (mesh.edge_flip_is_admissible(edge.index))
        {
            mesh.flip_edge(edge.index);
            flipped++;
        }
    }
    ASSERT_EQ(flipped, 40);

    slimed::SurfaceSolver solver;
    solver.build(mesh);
    const int n = solver.nVertices();

    // An arbitrary but reproducible nodal force.
    Matrix nodal = mat_calloc(n, 3);
    for (int v = 0; v < n; v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            nodal.set(v, axis, std::sin(0.37 * v + 1.9 * axis));
        }
    }

    Matrix surfaceForce = mat_calloc(n, 3);
    solver.nodal_force_to_surface(nodal, surfaceForce);

    // The defining property, checked directly: M^T F_S must reproduce F_C.
    // M^T x, assembled from the mask rows rather than from the solver, so a
    // shared error cannot hide.
    Matrix reconstructed = mat_calloc(n, 3);
    for (int v = 0; v < n; v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            reconstructed.set(v, axis, 0.5 * surfaceForce(v, axis));
        }
    }
    for (int v = 0; v < n; v++)
    {
        const int valence = static_cast<int>(mesh.vertices[v].adjacentVertices.size());
        for (int neighbour : mesh.vertices[v].adjacentVertices)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                reconstructed.set(neighbour, axis, reconstructed(neighbour, axis) +
                                                       surfaceForce(v, axis) * 0.5 / valence);
            }
        }
    }

    double worst = 0.0;
    for (int v = 0; v < n; v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            worst = std::max(worst, std::abs(reconstructed(v, axis) - nodal(v, axis)));
        }
    }
    EXPECT_LT(worst, 1e-8) << "M^T F_S does not reproduce F_C; worst error " << worst;

    // And it is genuinely a different map from the inverse. On an all-valence-6
    // mesh the two coincide, which is why nothing noticed until now.
    Matrix throughInverse = mat_calloc(n, 3);
    solver.surface_to_mesh(nodal, throughInverse);
    double largestDifference = 0.0;
    for (int v = 0; v < n; v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            largestDifference =
                std::max(largestDifference, std::abs(throughInverse(v, axis) - surfaceForce(v, axis)));
        }
    }
    EXPECT_GT(largestDifference, 1e-6)
        << "the transpose and the inverse agree, so this mesh does not test the distinction";
}

/**
 * @brief The solver notices a flip and rebuilds; the dense inverse cannot.
 */
TEST(FluidDynamicsTest, TheSolverRebuildsWhenTheConnectivityChanges)
{
    Param param;
    configure(param);
    Mesh mesh(param);
    const MeshFixture fixture = build_icosphere(1, 12.0);
    mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces);

    slimed::SurfaceSolver solver;
    solver.build(mesh);
    const long long versionAtBuild = solver.topologyVersion;
    EXPECT_EQ(versionAtBuild, mesh.topologyVersion);

    int candidate = -1;
    for (const MeshEdge &edge : mesh.edges)
    {
        if (mesh.edge_flip_is_admissible(edge.index))
        {
            candidate = edge.index;
            break;
        }
    }
    ASSERT_GE(candidate, 0);
    const int endpoint = mesh.edges[candidate].v[0];
    const int valenceBefore = static_cast<int>(mesh.vertices[endpoint].adjacentVertices.size());
    mesh.flip_edge(candidate);

    EXPECT_NE(solver.topologyVersion, mesh.topologyVersion)
        << "a flip must be visible to everything derived from the connectivity";

    solver.build(mesh);
    EXPECT_EQ(solver.topologyVersion, mesh.topologyVersion);

    // The rebuilt mask uses the new valence, which the flip changed.
    const int n = solver.nVertices();
    EXPECT_EQ(static_cast<int>(mesh.vertices[endpoint].adjacentVertices.size()),
              valenceBefore - 1);
    Matrix control = mat_calloc(n, 3);
    Matrix surface = mat_calloc(n, 3);
    for (int v = 0; v < n; v++)
    {
        control.set(v, 0, 1.0);
    }
    // The mask is affine: a constant control net has the same constant limit
    // surface, whatever the valences are.
    solver.mesh_to_surface(control, surface);
    for (int v = 0; v < n; v++)
    {
        EXPECT_NEAR(surface(v, 0), 1.0, 1e-14) << "vertex " << v;
    }
}

// ---------------------------------------------------------------------------
// In-plane motion
// ---------------------------------------------------------------------------

/**
 * @brief The in-plane displacement is switched off by default and on by the
 * parameter.
 *
 * Off is what this model has always done, and the fluctuation analysis relies
 * on it: every vertex stays on its ideal lattice site, so the height field can
 * be transformed without resampling. On is what a fluid membrane needs.
 */
TEST(FluidDynamicsTest, InPlaneMotionIsOffByDefaultAndOnWhenAsked)
{
    const auto run_one_step = [](bool inPlane) {
        Param param;
        configure_flat(param);
        param.inPlaneDynamicsEnabled = inPlane;
        param.timeStep = 0.001;
        param.diffConst = 1.0;
        DynamicMesh mesh(param);
        mesh.setup_flat();
        for (Vertex &vertex : mesh.vertices)
        {
            vertex.coord.set(2, 0, 10.0);
        }
        mesh.update_vertices_mat_with_vector();
        mesh.calculate_element_area_volume();
        mesh.sum_membrane_area_and_volume(param.area0, param.vol0);
        mesh.update_previous_coord_for_vertex();
        mesh.update_reference_coord_from_previous_coord();
        mesh.Compute_Energy_And_Force();

        Record record(4);
        DynamicModel model(mesh, record);
        const std::vector<double> beforeX = [&] {
            std::vector<double> x;
            for (const Vertex &vertex : mesh.vertices)
            {
                x.push_back(vertex.coord.get(0, 0));
            }
            return x;
        }();

        model.iteration = 0;
        mesh.apply_mesh_to_surface();
        model.next_step();
        mesh.apply_surface_to_mesh();
        mesh.update_vertices_vector_with_mat();

        double largestInPlaneMove = 0.0;
        for (std::size_t v = 0; v < mesh.vertices.size(); v++)
        {
            largestInPlaneMove = std::max(largestInPlaneMove,
                                          std::abs(mesh.vertices[v].coord.get(0, 0) - beforeX[v]));
        }
        return largestInPlaneMove;
    };

    EXPECT_LT(run_one_step(false), 1e-13) << "x must not move with in-plane dynamics off";
    EXPECT_GT(run_one_step(true), 1e-6) << "x must move with in-plane dynamics on";
}
