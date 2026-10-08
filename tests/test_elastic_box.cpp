#include "test_elastic_box.hpp"

/**
 * The elastic box lets the periodic dimensions stretch until the membrane's
 * frame tension sits at a prescribed value. Four things have to be right for
 * that to be a constant-tension ensemble rather than a knob:
 *
 * **The box force is the derivative of the energy.** The virial of the raw
 * nodal forces has to equal -dE/dlambda for an affine stretch of the whole
 * control net, images and ghosts included. The finite-difference checks
 * below stretch the net through the same routine the box uses and compare.
 * They are run on the per-vertex periodic sheet, on the global Periodic
 * mode's sheet, and on a strip that wraps along one axis only.
 *
 * **The stretch keeps the seams.** After a rescale an image must still be
 * its source plus its offset, the solver must agree, and the global mode's
 * ghost band must still be its partner plus the same lattice translation.
 *
 * **The box finds the target.** A pre-stretched flat sheet under the area
 * constraint has a tension that is known in closed form, so where the
 * deterministic box settles is known too, for a zero and a nonzero target.
 *
 * **The thermal box is canonical.** With its Brownian kick on, the strain of
 * the same flat sheet is a harmonic coordinate whose variance is kT over the
 * curvature of the effective potential.
 *
 * A run with the box off must be the run it always was; the byte-identity
 * gate on the shipped periodic workload covers that.
 *
 * @see docs/elastic_box.md
 */

namespace
{

void configure_common(Param &param, double side)
{
    param.VERBOSE_MODE = false;
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
    param.randomSeed = 20260928u;
    param.timeStep = 1.0e-3;
    param.diffConst = 1.0;
    param.surfaceSolver = "iterative";
    param.elasticBoxEnabled = true;
    param.elasticBoxTension = 0.0;
    param.elasticBoxMobility = 1.0e-3;
    param.elasticBoxThermal = false;
    param.elasticBoxCoupling = "isotropic";
    param.elasticBoxMaxStrainPerStep = 1.0e-3;
}

/// A flat sheet under the per-vertex mode, with the boundary of each axis.
void configure_mixed(Param &param, BoundaryType x, BoundaryType y, double side = 60.0)
{
    configure_common(param, side);
    param.boundaryCondition = BoundaryType::Mixed;
    param.boundaryConditionX = x;
    param.boundaryConditionY = y;
}

/// The same sheet under the global Periodic mode.
void configure_legacy_periodic(Param &param, double side = 60.0)
{
    configure_common(param, side);
    param.boundaryCondition = BoundaryType::Periodic;
}

/**
 * @brief Lift the sheet into a smooth relief that is continuous across the
 * seams, and put the images back on their sources.
 */
void undulate(Mesh &mesh, double amplitude)
{
    const Param &param = mesh.param;
    const bool mixed = (param.boundaryCondition == BoundaryType::Mixed);
    const bool periodicX = mixed ? (param.boundaryConditionX == BoundaryType::Periodic)
                                 : (param.boundaryCondition == BoundaryType::Periodic);
    const bool periodicY = mixed ? (param.boundaryConditionY == BoundaryType::Periodic)
                                 : (param.boundaryCondition == BoundaryType::Periodic);
    const double lx = (periodicX ? param.nFaceX - 6 : param.nFaceX) * param.dFaceX;
    const double ly = (periodicY ? param.nFaceY - 6 : param.nFaceY) * param.dFaceY;
    for (Vertex &vertex : mesh.vertices)
    {
        const double x = vertex.coord.get(0, 0);
        const double y = vertex.coord.get(1, 0);
        vertex.coord.set(2, 0,
                         amplitude * std::sin(2.0 * M_PI * x / lx) * std::cos(2.0 * M_PI * y / ly) +
                             10.0);
    }
    mesh.sync_periodic_images();
}

/**
 * @brief What the driver does before its first step: the flat sheet as the
 * reference and target area, an optional relief on top, and the first force
 * evaluation.
 */
void prepare(DynamicMesh &mesh, double amplitude)
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
    if (amplitude != 0.0)
    {
        undulate(mesh, amplitude);
    }
    mesh.Compute_Energy_And_Force();
    mesh.update_previous_force_for_vertex();
    mesh.update_vertices_mat_with_vector();
    mesh.apply_mesh_to_surface();
}

/// The energy the Brownian step drives with, as the model assembles it.
double drive_energy(const Mesh &mesh)
{
    const Energy &energy = mesh.param.energy;
    return energy.energyCurvature + energy.energyArea +
           (mesh.param.edgeSpringEnabled ? energy.energyRegularization : 0.0);
}

/// The virial of the same terms along one axis.
double drive_virial(const Mesh &mesh, int axis)
{
    return mesh.virialCurvature[axis] + mesh.virialArea[axis] +
           (mesh.param.edgeSpringEnabled ? mesh.virialRegularization[axis] : 0.0);
}

/**
 * @brief Central-difference dE/dlambda along one axis, through the same
 * rescale the box applies, leaving the mesh and its forces as they were.
 */
double finite_difference_dE_dlambda(DynamicMesh &mesh, int axis, double eps)
{
    std::array<double, 3> up = {{1.0, 1.0, 1.0}};
    std::array<double, 3> undoUp = {{1.0, 1.0, 1.0}};
    std::array<double, 3> down = {{1.0, 1.0, 1.0}};
    std::array<double, 3> undoDown = {{1.0, 1.0, 1.0}};
    up[axis] = 1.0 + eps;
    undoUp[axis] = 1.0 / (1.0 + eps);
    down[axis] = 1.0 - eps;
    undoDown[axis] = 1.0 / (1.0 - eps);

    mesh.apply_box_scale(up);
    mesh.Compute_Energy_And_Force();
    const double ePlus = drive_energy(mesh);
    mesh.apply_box_scale(undoUp);

    mesh.apply_box_scale(down);
    mesh.Compute_Energy_And_Force();
    const double eMinus = drive_energy(mesh);
    mesh.apply_box_scale(undoDown);

    mesh.Compute_Energy_And_Force();
    return (ePlus - eMinus) / (2.0 * eps);
}

void expect_virial_is_minus_the_derivative(DynamicMesh &mesh, int axis)
{
    const double eps = 1.0e-4;
    const double derivative = finite_difference_dE_dlambda(mesh, axis, eps);
    const double virial = drive_virial(mesh, axis);
    // The derivative is of order the energy, hundreds of pN.nm; the central
    // difference is good to eps^2 relative, the sum to rounding.
    EXPECT_NEAR(-virial, derivative, 1.0e-4 * std::max(1.0, std::abs(derivative)))
        << "axis " << axis << ": virial " << virial << ", dE/dlambda " << derivative;
    // And a live test: the terms must not all be zero.
    EXPECT_GT(std::abs(derivative), 1.0e-3) << "axis " << axis << " has no energy to differentiate";
}

/// A mesh, its record and its model, built the way the driver builds them.
struct BoxRun
{
    Param param;
    std::unique_ptr<DynamicMesh> mesh;
    Record record{16};
    std::unique_ptr<DynamicModel> model;

    explicit BoxRun(const Param &source, double amplitude = 0.0) : param(source)
    {
        mesh.reset(new DynamicMesh(param));
        mesh->setup_flat();
        prepare(*mesh, amplitude);
        record.add(mesh->param.area, mesh->param.energy, mesh->calculate_mean_force());
        model.reset(new DynamicModel(*mesh, record));
    }

    /// One iteration as the driver takes it.
    void step(int iteration)
    {
        mesh->apply_mesh_to_surface();
        model->iteration = iteration;
        model->next_step();
        mesh->apply_surface_to_mesh();
        mesh->postprocess_boundary();
        model->update_elastic_box();
        mesh->update_vertices_vector_with_mat();
        mesh->Compute_Energy_And_Force();
    }

    /// The box alone: no vertex motion, so a flat sheet stays flat and its
    /// tension is the area constraint's in closed form.
    void box_step(int iteration)
    {
        model->iteration = iteration;
        model->update_elastic_box();
        mesh->update_vertices_vector_with_mat();
        mesh->Compute_Energy_And_Force();
    }
};

} // namespace

TEST(ElasticBoxTest, TheParameterKeysParse)
{
    Param param;
    EXPECT_FALSE(param.elasticBoxEnabled);
    EXPECT_TRUE(import_kv_string("elasticBoxEnabled", "true", param));
    EXPECT_TRUE(param.elasticBoxEnabled);
    EXPECT_TRUE(import_kv_string("elasticBoxTension", "0.25", param));
    EXPECT_DOUBLE_EQ(param.elasticBoxTension, 0.25);
    EXPECT_TRUE(import_kv_string("elasticBoxMobility", "0.002", param));
    EXPECT_DOUBLE_EQ(param.elasticBoxMobility, 0.002);
    EXPECT_TRUE(import_kv_string("elasticBoxThermal", "false", param));
    EXPECT_FALSE(param.elasticBoxThermal);
    EXPECT_TRUE(import_kv_string("elasticBoxCoupling", "anisotropic", param));
    EXPECT_EQ(param.elasticBoxCoupling, "anisotropic");
    EXPECT_TRUE(import_kv_string("elasticBoxMaxStrainPerStep", "0.01", param));
    EXPECT_DOUBLE_EQ(param.elasticBoxMaxStrainPerStep, 0.01);
    EXPECT_TRUE(import_kv_string("elasticBoxEnabled", "false", param));
    EXPECT_FALSE(param.elasticBoxEnabled);
}

TEST(ElasticBoxTest, TheVirialIsMinusTheDerivativeOfTheEnergyOnThePeriodicSheet)
{
    // Bending live through the relief, the area constraint live because the
    // relief has excess area, and the tether live in its harmonic form.
    Param param;
    configure_mixed(param, BoundaryType::Periodic, BoundaryType::Periodic);
    param.edgeSpringEnabled = true;
    param.edgeTetherShape = "harmonic";
    DynamicMesh mesh(param);
    mesh.setup_flat();
    prepare(mesh, 3.0);
    ASSERT_EQ(mesh.box_periodic_axis_count(), 2);
    expect_virial_is_minus_the_derivative(mesh, 0);
    expect_virial_is_minus_the_derivative(mesh, 1);
}

TEST(ElasticBoxTest, TheVirialIsMinusTheDerivativeOfTheEnergyUnderTheGlobalPeriodicMode)
{
    // The global mode zeroes its ghosts' forces after the evaluation; the
    // virial has to be taken before that, or the seam is lost.
    Param param;
    configure_legacy_periodic(param);
    DynamicMesh mesh(param);
    mesh.setup_flat();
    prepare(mesh, 3.0);
    ASSERT_EQ(mesh.box_periodic_axis_count(), 2);
    EXPECT_NEAR(mesh.boxLength[0], (param.nFaceX - 6) * param.dFaceX, 1.0e-12);
    EXPECT_NEAR(mesh.boxLength[1], (param.nFaceY - 6) * param.dFaceY, 1.0e-12);
    expect_virial_is_minus_the_derivative(mesh, 0);
    expect_virial_is_minus_the_derivative(mesh, 1);
}

TEST(ElasticBoxTest, TheVirialIsMinusTheDerivativeOfTheEnergyOnAStrip)
{
    // Wraps along x, clamped along y: one periodic axis, and the clamped
    // rows stretch with the frame like everything else.
    Param param;
    configure_mixed(param, BoundaryType::Periodic, BoundaryType::Fixed);
    param.edgeSpringEnabled = true;
    param.edgeTetherShape = "harmonic";
    DynamicMesh mesh(param);
    mesh.setup_flat();
    prepare(mesh, 3.0);
    ASSERT_EQ(mesh.box_periodic_axis_count(), 1);
    EXPECT_TRUE(mesh.boxPeriodicAxis[0]);
    EXPECT_FALSE(mesh.boxPeriodicAxis[1]);
    expect_virial_is_minus_the_derivative(mesh, 0);
}

TEST(ElasticBoxTest, TheStretchKeepsEveryImageOnItsSourceAndTheSolverAgrees)
{
    Param param;
    configure_mixed(param, BoundaryType::Periodic, BoundaryType::Periodic);
    BoxRun run(param, 2.0);
    DynamicMesh &mesh = *run.mesh;
    const double lengthX = mesh.boxLength[0];
    const double lengthY = mesh.boxLength[1];

    std::vector<std::array<double, 3>> offsetBefore(mesh.vertices.size());
    for (int v : mesh.periodicImageVertices)
    {
        offsetBefore[v] = mesh.vertices[v].mirrorOffset;
    }

    const std::array<double, 3> factor = {{1.01, 0.99, 1.0}};
    mesh.apply_box_scale(factor);
    EXPECT_NEAR(mesh.boxLength[0], 1.01 * lengthX, 1.0e-12);
    EXPECT_NEAR(mesh.boxLength[1], 0.99 * lengthY, 1.0e-12);
    EXPECT_NEAR(mesh.boxScale[0], 1.01, 1.0e-15);
    EXPECT_NEAR(mesh.boxScale[1], 0.99, 1.0e-15);
    for (int v : mesh.periodicImageVertices)
    {
        const Vertex &image = mesh.vertices[v];
        const Vertex &source = mesh.vertices[image.reflectiveVertexIndex];
        for (int axis = 0; axis < 3; axis++)
        {
            EXPECT_NEAR(image.mirrorOffset[axis], factor[axis] * offsetBefore[v][axis], 1.0e-12);
            EXPECT_NEAR(image.coord(axis, 0), source.coord(axis, 0) + image.mirrorOffset[axis],
                        1.0e-9);
            EXPECT_NEAR(mesh.matMesh(v, axis), mesh.matMesh(source.index, axis) + image.mirrorOffset[axis],
                        1.0e-9);
        }
    }

    // A step through the solver: it writes the image rows itself from the
    // offsets it holds, which must be the stretched ones.
    mesh.Compute_Energy_And_Force();
    run.model->iteration = 0;
    mesh.apply_mesh_to_surface();
    run.model->next_step();
    mesh.apply_surface_to_mesh();
    for (int v : mesh.periodicImageVertices)
    {
        const Vertex &image = mesh.vertices[v];
        const int source = image.reflectiveVertexIndex;
        for (int axis = 0; axis < 3; axis++)
        {
            EXPECT_NEAR(mesh.matMesh(v, axis), mesh.matMesh(source, axis) + image.mirrorOffset[axis],
                        1.0e-9)
                << "image " << v << " axis " << axis;
        }
    }
}

TEST(ElasticBoxTest, TheStretchKeepsTheGhostBandOfTheGlobalPeriodicMode)
{
    Param param;
    configure_legacy_periodic(param);
    BoxRun run(param, 2.0);
    DynamicMesh &mesh = *run.mesh;

    // The first step stamps every ghost and duplicate with its partner.
    run.step(0);
    // The box is live in this run and takes a strain step of its own inside
    // step(); the offsets must follow the whole stretch, not only the one
    // applied by hand.
    const std::array<double, 3> scaleBefore = mesh.boxScale;
    std::vector<std::array<double, 3>> offsetBefore(mesh.vertices.size());
    int nPaired = 0;
    for (const Vertex &vertex : mesh.vertices)
    {
        if (vertex.reflectiveVertexIndex < 0)
        {
            continue;
        }
        const Vertex &partner = mesh.vertices[vertex.reflectiveVertexIndex];
        for (int axis = 0; axis < 3; axis++)
        {
            offsetBefore[vertex.index][axis] = vertex.coord(axis, 0) - partner.coord(axis, 0);
        }
        nPaired++;
    }
    ASSERT_GT(nPaired, 0);

    const std::array<double, 3> factor = {{1.02, 0.98, 1.0}};
    mesh.apply_box_scale(factor);
    // Another step: the post-process copies displacements, so the offsets
    // it keeps are whatever the stretch left them.
    run.step(1);
    for (int axis = 0; axis < 2; axis++)
    {
        EXPECT_NEAR(mesh.boxScale[axis] / scaleBefore[axis], factor[axis], 1.0e-4);
    }
    for (const Vertex &vertex : mesh.vertices)
    {
        if (vertex.reflectiveVertexIndex < 0)
        {
            continue;
        }
        const Vertex &partner = mesh.vertices[vertex.reflectiveVertexIndex];
        for (int axis = 0; axis < 2; axis++)
        {
            const double stretch = mesh.boxScale[axis] / scaleBefore[axis];
            EXPECT_NEAR(vertex.coord(axis, 0) - partner.coord(axis, 0),
                        stretch * offsetBefore[vertex.index][axis], 1.0e-9)
                << "vertex " << vertex.index << " axis " << axis;
        }
    }
}

TEST(ElasticBoxTest, TheDeterministicBoxRelaxesAPreStretchedSheetToTheTargetTension)
{
    // A flat sheet stretched by 2% along both axes carries the area
    // constraint's tension K (A - A0) / A0 and nothing else; the box must
    // bring it to A_p = A0 (1 + tau / K).
    Param param;
    configure_mixed(param, BoundaryType::Periodic, BoundaryType::Periodic);
    param.elasticBoxMobility = 0.1;
    param.elasticBoxMaxStrainPerStep = 1.0e-2;
    for (const double target : {0.0, 2.0})
    {
        param.elasticBoxTension = target;
        BoxRun run(param, 0.0);
        DynamicMesh &mesh = *run.mesh;
        const double area0 = mesh.param.area0;
        ASSERT_NEAR(mesh.box_projected_area(), area0, 1.0e-9 * area0);

        mesh.apply_box_scale({{1.02, 1.02, 1.0}});
        mesh.Compute_Energy_And_Force();
        run.model->refresh_box_tension();
        const double expectedStart = param.uSurf * (1.02 * 1.02 - 1.0);
        EXPECT_NEAR(run.model->boxTensionMean, expectedStart, 1.0e-6 * expectedStart);
        EXPECT_NEAR(run.model->boxTensionByTerm[1], expectedStart, 1.0e-6 * expectedStart);
        EXPECT_NEAR(run.model->boxTensionByTerm[0], 0.0, 1.0e-9);

        for (int iteration = 0; iteration < 600; iteration++)
        {
            run.box_step(iteration);
        }
        const double expectedArea = area0 * (1.0 + target / param.uSurf);
        EXPECT_NEAR(mesh.box_projected_area(), expectedArea, 1.0e-6 * expectedArea)
            << "target " << target;
        EXPECT_NEAR(mesh.param.area, expectedArea, 1.0e-6 * expectedArea);
        EXPECT_NEAR(run.model->boxTensionMean, target, 1.0e-4);
        // Isotropic: the box keeps its shape.
        EXPECT_DOUBLE_EQ(mesh.boxScale[0], mesh.boxScale[1]);
        EXPECT_DOUBLE_EQ(mesh.boxScale[2], 1.0);
    }
}

TEST(ElasticBoxTest, TheAnisotropicStripStretchesOnlyItsPeriodicAxis)
{
    Param param;
    configure_mixed(param, BoundaryType::Periodic, BoundaryType::Fixed);
    param.elasticBoxCoupling = "anisotropic";
    param.elasticBoxMobility = 0.1;
    param.elasticBoxMaxStrainPerStep = 1.0e-2;
    BoxRun run(param, 0.0);
    DynamicMesh &mesh = *run.mesh;
    const double area0 = mesh.param.area0;
    ASSERT_EQ(mesh.box_periodic_axis_count(), 1);
    EXPECT_NEAR(mesh.box_projected_area(), area0, 1.0e-9 * area0);
    const double lengthX = mesh.boxLength[0];

    std::vector<double> yBefore;
    for (const Vertex &vertex : mesh.vertices)
    {
        yBefore.push_back(vertex.coord(1, 0));
    }
    mesh.apply_box_scale({{1.02, 1.0, 1.0}});
    mesh.Compute_Energy_And_Force();
    for (int iteration = 0; iteration < 600; iteration++)
    {
        run.box_step(iteration);
    }
    EXPECT_NEAR(mesh.boxLength[0], lengthX, 1.0e-6 * lengthX);
    EXPECT_NEAR(mesh.box_projected_area(), area0, 1.0e-6 * area0);
    EXPECT_DOUBLE_EQ(mesh.boxScale[1], 1.0);
    for (std::size_t v = 0; v < mesh.vertices.size(); v++)
    {
        EXPECT_DOUBLE_EQ(mesh.vertices[v].coord(1, 0), yBefore[v]);
    }
}

TEST(ElasticBoxTest, TheThermalBoxSamplesTheStrainVarianceOfTheEffectivePotential)
{
    // Flat sheet, no vertex motion: U(lambda) = (K A0 / 2) (lambda^2 - 1)^2
    // for the isotropic strain of two axes, so var(lambda) = kT / (4 K A0)
    // -- kT / (n^2 K A_p) in general -- up to the explicit step's factor
    // 2 / (2 - rate dt).
    Param param;
    configure_mixed(param, BoundaryType::Periodic, BoundaryType::Periodic, 50.0);
    param.elasticBoxThermal = true;
    param.elasticBoxMobility = 0.1;
    param.elasticBoxMaxStrainPerStep = 1.0e-1;
    BoxRun run(param, 0.0);
    DynamicMesh &mesh = *run.mesh;
    const double area0 = mesh.param.area0;
    const double rateDt = 2.0 * param.elasticBoxMobility * param.uSurf * param.timeStep;
    const double expectedVariance = param.KBT / (4.0 * param.uSurf * area0) * 2.0 / (2.0 - rateDt);

    const int nSteps = 6000;
    const int nBurn = 500;
    std::vector<double> strain;
    strain.reserve(nSteps - nBurn);
    for (int iteration = 0; iteration < nSteps; iteration++)
    {
        run.box_step(iteration);
        if (iteration >= nBurn)
        {
            strain.push_back(mesh.boxScale[0] - 1.0);
        }
    }
    const double mean = std::accumulate(strain.begin(), strain.end(), 0.0) / strain.size();
    double variance = 0.0;
    for (double s : strain)
    {
        variance += (s - mean) * (s - mean);
    }
    variance /= strain.size();
    // About 275 independent samples at a correlation time of 20 steps.
    EXPECT_NEAR(mean, 0.0, 4.0 * std::sqrt(expectedVariance / 275.0));
    EXPECT_GT(variance / expectedVariance, 0.7);
    EXPECT_LT(variance / expectedVariance, 1.35);
    EXPECT_EQ(run.model->boxClippedSteps, 0);
}

TEST(ElasticBoxTest, TheVertexGasTermShiftsTheBoxByNkTOverTheArea)
{
    // With in-plane motion on, the integrated vertices are particles in the
    // plane and the box carries their -N kT ln(A_p): the flat sheet then
    // settles where K (A - A0) / A0 = tau + N kT / A_p, and where
    // K (A - A0) / A0 = tau with the term switched off.
    Param param;
    configure_mixed(param, BoundaryType::Periodic, BoundaryType::Periodic);
    param.inPlaneDynamicsEnabled = true;
    param.elasticBoxMobility = 0.1;
    param.elasticBoxMaxStrainPerStep = 1.0e-2;
    for (const bool gasTerm : {true, false})
    {
        param.elasticBoxVertexGasTerm = gasTerm;
        BoxRun run(param, 0.0);
        DynamicMesh &mesh = *run.mesh;
        const double area0 = mesh.param.area0;
        // The sources of the periodic sheet, none of them clamped.
        const int expectedCount =
            gasTerm ? (run.param.nFaceX - 6) * (run.param.nFaceY - 6) : 0;
        EXPECT_EQ(run.model->boxInPlaneVertexCount, expectedCount);

        mesh.apply_box_scale({{1.01, 1.01, 1.0}});
        mesh.Compute_Energy_And_Force();
        for (int iteration = 0; iteration < 800; iteration++)
        {
            run.box_step(iteration);
        }
        const double areaP = mesh.box_projected_area();
        const double internalTension = param.uSurf * (mesh.param.area - area0) / area0;
        const double gas = gasTerm ? expectedCount * param.KBT / areaP : 0.0;
        EXPECT_NEAR(internalTension, gas, 1.0e-4) << "gas term " << gasTerm;
        EXPECT_NEAR(run.model->boxTensionMean, 0.0, 1.0e-4);
        EXPECT_NEAR(run.model->boxTensionVertexGas, -gas, 1.0e-9);
        if (gasTerm)
        {
            // About 0.19 pN/nm at K = 250: the box is visibly larger than A0.
            EXPECT_GT(areaP, area0 * (1.0 + 5.0e-4));
        }
        else
        {
            EXPECT_NEAR(areaP, area0, 1.0e-6 * area0);
        }
    }
}

TEST(ElasticBoxTest, TheBoxIsRefusedWhereNothingWraps)
{
    Param param;
    configure_mixed(param, BoundaryType::Fixed, BoundaryType::Fixed);
    EXPECT_THROW(BoxRun run(param, 0.0), std::runtime_error);

    Param bad;
    configure_mixed(bad, BoundaryType::Periodic, BoundaryType::Periodic);
    bad.elasticBoxCoupling = "sideways";
    EXPECT_THROW(BoxRun run(bad, 0.0), std::runtime_error);

    Param zero;
    configure_mixed(zero, BoundaryType::Periodic, BoundaryType::Periodic);
    zero.elasticBoxMobility = 0.0;
    EXPECT_THROW(BoxRun run(zero, 0.0), std::runtime_error);
}

TEST(ElasticBoxTest, WithTheBoxOffNothingIsTouched)
{
    Param param;
    configure_mixed(param, BoundaryType::Periodic, BoundaryType::Periodic);
    param.elasticBoxEnabled = false;
    BoxRun run(param, 2.0);
    DynamicMesh &mesh = *run.mesh;
    for (int axis = 0; axis < 3; axis++)
    {
        EXPECT_DOUBLE_EQ(mesh.virialCurvature[axis], 0.0);
        EXPECT_DOUBLE_EQ(mesh.virialArea[axis], 0.0);
    }
    const double lengthX = mesh.boxLength[0];
    run.step(0);
    run.step(1);
    EXPECT_DOUBLE_EQ(mesh.boxScale[0], 1.0);
    EXPECT_DOUBLE_EQ(mesh.boxScale[1], 1.0);
    EXPECT_DOUBLE_EQ(mesh.boxLength[0], lengthX);
    EXPECT_EQ(run.model->boxSteps, 0);
}
