#include "test_dynamics_totals.hpp"

/**
 * DynamicModel::next_step() drives the membrane with the bending and area
 * forces, and takes the regularization slot in only when the edge spring puts
 * the fluid mesh-quality terms there. With the spring off the slot holds the
 * reference-length regularization: evaluated every step, never felt by the
 * trajectory.
 *
 * It used to be summed into E_Total and the mean force regardless. Measured
 * on the 100 nm periodic sheet, three 3000-step runs with the regularizer
 * scaled by 0, 1 and 25 wrote byte-identical trajectories and mean forces of
 * 3.0, 11.7 and 259 pN; in the 800 000-step run two thirds of the logged
 * total energy was the term the run did not sample.
 *
 * What is pinned here:
 *
 * - a solid Brownian run reports totals without the slot, and those totals
 *   are exactly what the step is driven by;
 * - they are a function of the trajectory alone -- moving the reference
 *   configuration changes the diagnostic column and nothing else;
 * - the slot itself is still written, so the diagnostic is not lost;
 * - the minimizer, which descends the total, keeps the term in it;
 * - a fluid run, whose drive has the slot, keeps it too.
 *
 * @see Mesh::regularizationInTotals
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
    param.KBT = 4.17;
    param.uSurf = 250.0;
    param.uVol = 0.0;
    param.isGlobalConstraint = true;
    param.isEnergyHarmonicBondIncluded = false;
    param.randomSeed = 20261001u;
    param.timeStep = 1.0e-3;
    param.diffConst = 1.0;
}

/// A smooth relief that is continuous across the periodic seam.
void undulate(Mesh &mesh, double amplitude)
{
    const Param &param = mesh.param;
    const double lx = (param.nFaceX - 6) * param.dFaceX;
    const double ly = (param.nFaceY - 6) * param.dFaceY;
    for (Vertex &vertex : mesh.vertices)
    {
        const double x = vertex.coord.get(0, 0);
        const double y = vertex.coord.get(1, 0);
        vertex.coord.set(2, 0,
                         amplitude * std::sin(2.0 * M_PI * x / lx) * std::cos(2.0 * M_PI * y / ly) +
                             10.0);
    }
}

/**
 * @brief What the drivers do before their first force evaluation, with the
 * flat sheet as the reference configuration and a relief on top of it.
 *
 * The reference is taken flat so that the reference-length regularization is
 * live: a term that was zero could be in a total or out of it unnoticed.
 */
void prepare(Mesh &mesh, double amplitude)
{
    for (Vertex &vertex : mesh.vertices)
    {
        vertex.coord.set(2, 0, 10.0);
    }
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    undulate(mesh, amplitude);
    mesh.update_previous_coord_for_vertex();
    mesh.Compute_Energy_And_Force();
    mesh.update_previous_force_for_vertex();
}

/// Every term of the total but the regularization, in the order Energy sums them.
double total_without_regularization(const Energy &energy)
{
    return energy.energyCurvature + energy.energyArea + energy.energyVolume +
           energy.energyThickness + energy.energyTilt + energy.energyHarmonicBond +
           energy.energyGagScaffolding + energy.energyIdealizedProteinLattice;
}

double largest_component(const Mesh &mesh, Matrix Force::*term)
{
    double largest = 0.0;
    for (const Vertex &vertex : mesh.vertices)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            largest = std::max(largest, std::abs((vertex.force.*term).get(axis, 0)));
        }
    }
    return largest;
}

/// The driver's loop body, as run_dynamics_flat() has it.
void step(DynamicMesh &mesh, DynamicModel &model, Record &record)
{
    mesh.apply_mesh_to_surface();
    model.next_step();
    mesh.apply_surface_to_mesh();
    mesh.postprocess_boundary();
    mesh.update_vertices_vector_with_mat();
    record.add(mesh.param.area, mesh.param.energy, mesh.calculate_mean_force());
    mesh.Compute_Energy_And_Force();
}

} // namespace

// ---------------------------------------------------------------------------
// The two sums
// ---------------------------------------------------------------------------

TEST(DynamicsTotalsTest, TheEnergySumTakesTheRegularizationOrLeavesIt)
{
    Energy energy;
    energy.energyCurvature = 505.0;
    energy.energyArea = 11.0;
    energy.energyVolume = 0.25;
    energy.energyRegularization = 978.0;
    energy.energyHarmonicBond = 3.0;

    // The default is the sum every caller had before the switch existed.
    EXPECT_DOUBLE_EQ(energy.calculateTotalEnergy(), 505.0 + 11.0 + 0.25 + 978.0 + 3.0);
    EXPECT_DOUBLE_EQ(energy.calculateTotalEnergy(true), 505.0 + 11.0 + 0.25 + 978.0 + 3.0);
    EXPECT_DOUBLE_EQ(energy.calculateTotalEnergy(false), 505.0 + 11.0 + 0.25 + 3.0);
    EXPECT_DOUBLE_EQ(energy.energyTotal, 505.0 + 11.0 + 0.25 + 3.0);
    // Leaving it out of the sum does not touch the term.
    EXPECT_DOUBLE_EQ(energy.energyRegularization, 978.0);
}

TEST(DynamicsTotalsTest, TheForceSumTakesTheRegularizationOrLeavesIt)
{
    Force force;
    force.forceCurvature.set(2, 0, 1.5);
    force.forceArea.set(2, 0, -0.25);
    force.forceRegularization.set(2, 0, 40.0);
    force.forceRegularization.set(0, 0, -7.0);

    force.calculate_total_force();
    EXPECT_DOUBLE_EQ(force.forceTotal.get(2, 0), 1.5 - 0.25 + 40.0);
    EXPECT_DOUBLE_EQ(force.forceTotal.get(0, 0), -7.0);

    force.calculate_total_force(false);
    EXPECT_DOUBLE_EQ(force.forceTotal.get(2, 0), 1.5 - 0.25);
    EXPECT_DOUBLE_EQ(force.forceTotal.get(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(force.forceRegularization.get(2, 0), 40.0);
}

// ---------------------------------------------------------------------------
// A solid Brownian run: the slot is a diagnostic
// ---------------------------------------------------------------------------

TEST(DynamicsTotalsTest, ASolidRunLeavesTheRegularizationOutOfItsTotals)
{
    Param param;
    configure_sheet(param);
    ASSERT_FALSE(param.edgeSpringEnabled);
    DynamicMesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_flat());
    EXPECT_FALSE(mesh.regularizationInTotals);
    prepare(mesh, 1.5);

    const Energy &energy = mesh.param.energy;
    // The diagnostic is still evaluated, energy and force both.
    ASSERT_GT(energy.energyRegularization, 1.0) << "the fixture's regularization is not live";
    ASSERT_GT(largest_component(mesh, &Force::forceRegularization), 1e-3);
    ASSERT_GT(energy.energyCurvature, 1.0);

    EXPECT_EQ(energy.energyTotal, total_without_regularization(energy));
    for (const Face &face : mesh.faces)
    {
        ASSERT_EQ(face.energy.energyTotal, total_without_regularization(face.energy))
            << "face " << face.index;
    }

    // The total force is the force the Brownian step is driven by: bending
    // plus area. (The sheet encloses no volume, so that term is zero.)
    double sumMagnitude = 0.0;
    for (const Vertex &vertex : mesh.vertices)
    {
        double squared = 0.0;
        for (int axis = 0; axis < 3; axis++)
        {
            const double drive =
                vertex.force.forceCurvature.get(axis, 0) + vertex.force.forceArea.get(axis, 0);
            ASSERT_DOUBLE_EQ(vertex.force.forceTotal.get(axis, 0), drive)
                << "vertex " << vertex.index << ", axis " << axis;
            squared += drive * drive;
        }
        sumMagnitude += std::sqrt(squared);
    }
    const double meanDrive = sumMagnitude / mesh.vertices.size();
    ASSERT_GT(meanDrive, 1e-3);
    EXPECT_NEAR(mesh.calculate_mean_force(), meanDrive, 1e-12 * meanDrive);
}

/**
 * The measurement that found this, as a test: two runs that differ only in
 * the reference configuration the regularization measures against. The
 * trajectories are the same (to the last bit in the serial build, to the
 * order of an OpenMP reduction otherwise), so everything that claims to
 * describe the sampled Hamiltonian has to be too.
 */
TEST(DynamicsTotalsTest, ASolidRunsTotalsDoNotDependOnTheReferenceConfiguration)
{
    Param paramA;
    Param paramB;
    configure_sheet(paramA);
    configure_sheet(paramB);
    DynamicMesh meshA(paramA);
    DynamicMesh meshB(paramB);
    ASSERT_NO_THROW(meshA.setup_flat());
    ASSERT_NO_THROW(meshB.setup_flat());
    prepare(meshA, 1.0);
    prepare(meshB, 1.0);

    // B remembers a different starting configuration: a sheet stretched by a
    // tenth in the plane.
    for (Vertex &vertex : meshB.vertices)
    {
        vertex.coordRef.set(0, 0, 1.1 * vertex.coordRef.get(0, 0));
        vertex.coordRef.set(1, 0, 1.1 * vertex.coordRef.get(1, 0));
    }
    meshB.Compute_Energy_And_Force();
    ASSERT_GT(std::abs(meshB.param.energy.energyRegularization -
                       meshA.param.energy.energyRegularization),
              1.0)
        << "the two references give the same regularization; the test would pass vacuously";

    const int nSteps = 20;
    const double kCoordinateTolerance = 1e-9; // nm
    const double kRelativeTolerance = 1e-9;
    Record recordA(nSteps + 1);
    Record recordB(nSteps + 1);
    recordA.add(meshA.param.area, meshA.param.energy, meshA.calculate_mean_force());
    recordB.add(meshB.param.area, meshB.param.energy, meshB.calculate_mean_force());
    DynamicModel modelA(meshA, recordA);
    DynamicModel modelB(meshB, recordB);
    meshA.update_vertices_mat_with_vector();
    meshB.update_vertices_mat_with_vector();

    for (int iteration = 0; iteration < nSteps; iteration++)
    {
        modelA.iteration = iteration;
        modelB.iteration = iteration;
        step(meshA, modelA, recordA);
        step(meshB, modelB, recordB);

        // In the serial build the two trajectories are the same to the last
        // bit. In the OpenMP build they are not quite: the membrane area that
        // the global area constraint acts on is an OpenMP reduction, whose
        // order varies from one call to the next, and that last-bit noise
        // reaches every force and so every coordinate (about 1e-13 nm a
        // step). The tolerances are set where that noise cannot reach and
        // the regularization force could not hide: with the reference
        // stretched by a tenth it would move a vertex by ~1e-3 nm a step and
        // shift the totals by more than 1 pN.nm (the assertion before the
        // loop), against 1e-9 nm and 1e-9 relative here.
        for (std::size_t v = 0; v < meshA.vertices.size(); v++)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                ASSERT_NEAR(meshA.vertices[v].coord.get(axis, 0), meshB.vertices[v].coord.get(axis, 0),
                            kCoordinateTolerance)
                    << "the trajectories part at step " << iteration << ", vertex " << v;
            }
        }
        EXPECT_NEAR(meshA.param.energy.energyTotal, meshB.param.energy.energyTotal,
                    kRelativeTolerance * std::abs(meshA.param.energy.energyTotal))
            << "step " << iteration;
        const double meanForceA = meshA.calculate_mean_force();
        EXPECT_NEAR(meanForceA, meshB.calculate_mean_force(), kRelativeTolerance * meanForceA)
            << "step " << iteration;
        EXPECT_NE(meshA.param.energy.energyRegularization, meshB.param.energy.energyRegularization)
            << "step " << iteration;
    }
    // What reaches EnergyForce.csv.
    for (std::size_t row = 0; row < recordA.energyVec.size(); row++)
    {
        EXPECT_NEAR(recordA.energyVec[row].energyTotal, recordB.energyVec[row].energyTotal,
                    kRelativeTolerance * std::abs(recordA.energyVec[row].energyTotal))
            << row;
        EXPECT_NEAR(recordA.meanForce[row], recordB.meanForce[row],
                    kRelativeTolerance * recordA.meanForce[row])
            << row;
    }
}

// ---------------------------------------------------------------------------
// Where the slot is in the Hamiltonian, it stays in the totals
// ---------------------------------------------------------------------------

TEST(DynamicsTotalsTest, TheMinimizerKeepsTheRegularizationInItsTotals)
{
    Param param;
    configure_sheet(param);
    Mesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_flat());
    EXPECT_TRUE(mesh.regularizationInTotals);
    prepare(mesh, 1.5);

    const Energy &energy = mesh.param.energy;
    ASSERT_GT(energy.energyRegularization, 1.0);
    // The sum Energy::calculateTotalEnergy() has always made, term for term.
    EXPECT_EQ(energy.energyTotal,
              energy.energyCurvature + energy.energyArea + energy.energyVolume +
                  energy.energyThickness + energy.energyTilt + energy.energyRegularization +
                  energy.energyHarmonicBond + energy.energyGagScaffolding +
                  energy.energyIdealizedProteinLattice);

    bool regularizationReachesTheTotal = false;
    for (const Vertex &vertex : mesh.vertices)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            const Force &force = vertex.force;
            ASSERT_DOUBLE_EQ(force.forceTotal.get(axis, 0),
                             force.forceCurvature.get(axis, 0) + force.forceArea.get(axis, 0) +
                                 force.forceVolume.get(axis, 0) +
                                 force.forceRegularization.get(axis, 0))
                << "vertex " << vertex.index << ", axis " << axis;
            regularizationReachesTheTotal = regularizationReachesTheTotal ||
                                            std::abs(force.forceRegularization.get(axis, 0)) > 1e-3;
        }
    }
    EXPECT_TRUE(regularizationReachesTheTotal);
}

TEST(DynamicsTotalsTest, AFluidRunKeepsTheMeshTermInItsTotals)
{
    Param param;
    configure_sheet(param);
    param.edgeSpringEnabled = true;
    // The lower wall above the lattice spacing, so every tethered edge is
    // against it and the slot is not empty.
    param.edgeTetherMinRatio = 1.05;
    DynamicMesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_flat());
    EXPECT_TRUE(mesh.regularizationInTotals);
    prepare(mesh, 1.5);

    const Energy &energy = mesh.param.energy;
    ASSERT_GT(energy.energyRegularization, 1.0) << "no edge is on the tether's wall";
    EXPECT_EQ(energy.energyTotal,
              energy.energyCurvature + energy.energyArea + energy.energyVolume +
                  energy.energyThickness + energy.energyTilt + energy.energyRegularization +
                  energy.energyHarmonicBond + energy.energyGagScaffolding +
                  energy.energyIdealizedProteinLattice);

    bool meshTermReachesTheTotal = false;
    for (const Vertex &vertex : mesh.vertices)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            const Force &force = vertex.force;
            // What next_step() assembles when the spring is on.
            ASSERT_DOUBLE_EQ(force.forceTotal.get(axis, 0),
                             force.forceCurvature.get(axis, 0) + force.forceArea.get(axis, 0) +
                                 force.forceVolume.get(axis, 0) +
                                 force.forceRegularization.get(axis, 0))
                << "vertex " << vertex.index << ", axis " << axis;
            meshTermReachesTheTotal =
                meshTermReachesTheTotal || std::abs(force.forceRegularization.get(axis, 0)) > 1e-3;
        }
    }
    EXPECT_TRUE(meshTermReachesTheTotal);
}

// ---------------------------------------------------------------------------
// The file says which it is
// ---------------------------------------------------------------------------

TEST(DynamicsTotalsTest, TheCsvHeaderSaysWhenTheRegularizationIsADiagnostic)
{
    const std::string inTotals = energy_force_csv_header(true);
    const std::string diagnostic = energy_force_csv_header(false);

    // The header every existing EnergyForce.csv carries, unchanged.
    EXPECT_EQ(inTotals, "E_Curvature, E_Area, E_Regularization, E_HarmonicBond, E_GagScaffolding, "
                        "E_IdealizedProteinLattice, E_Total ((pN.nm)), Mean Force (pN)");

    EXPECT_NE(diagnostic, inTotals);
    EXPECT_NE(diagnostic.find("E_Regularization (diagnostic; not in E_Total or Mean Force)"),
              std::string::npos);
    // The same eight columns: a reader that splits on commas sees no change.
    EXPECT_EQ(std::count(inTotals.begin(), inTotals.end(), ','), 7);
    EXPECT_EQ(std::count(diagnostic.begin(), diagnostic.end(), ','), 7);
    EXPECT_EQ(diagnostic.find('\n'), std::string::npos);
}
