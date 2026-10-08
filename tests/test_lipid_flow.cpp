#include "test_lipid_flow.hpp"

/**
 * The lipid-flow prototype gives every face its own relaxed area and lets it
 * move. Each piece has a closed-form expectation, and each is tested here:
 *
 * - the force on the control points is minus the gradient of the local
 *   stretching energy, with the relaxed areas held;
 * - uniform relaxed areas on a uniformly stretched sheet are the global
 *   constraint, energy and force;
 * - the lateral flow conserves the total relaxed area to rounding, only ever
 *   lowers the energy, and ends at a uniform exchange tension where the
 *   energy is exactly the global constraint's;
 * - a tension wave on a flat sheet decays at D_sigma q^2;
 * - with the noise on, the relaxed area of a face has the variance kT a / K;
 * - a reservoir brings every face it touches to its tension;
 * - the flux graph has one link per physical edge, across a periodic seam
 *   included.
 *
 * @see docs/lipid_flow_plan.md
 */

namespace
{

void configure(Param &param, BoundaryType x, BoundaryType y, double side)
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
    param.randomSeed = 20261005u;
    param.timeStep = 1.0e-3;
    param.diffConst = 1.0;
    param.surfaceSolver = "iterative";
    param.boundaryCondition = BoundaryType::Mixed;
    param.boundaryConditionX = x;
    param.boundaryConditionY = y;
    param.lipidFlowMode = "local";
    param.lipidFlowDiffusivity = 50.0;
    param.lipidFlowSubsteps = 0;
    param.lipidFlowThermal = false;
    param.lipidReservoirRegion = "none";
}

void undulate(Mesh &mesh, double amplitude)
{
    const Param &param = mesh.param;
    const bool periodicX = (param.boundaryConditionX == BoundaryType::Periodic);
    const bool periodicY = (param.boundaryConditionY == BoundaryType::Periodic);
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

/// What the driver does before its first step, the relaxed areas included.
void prepare(DynamicMesh &mesh, double amplitude)
{
    for (Vertex &vertex : mesh.vertices)
    {
        vertex.coord.set(2, 0, 10.0);
    }
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.param.vol0 = 0.0;
    if (mesh.lipid_flow_on())
    {
        mesh.initialize_lipid_flow();
    }
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

struct LipidRun
{
    Param param;
    std::unique_ptr<DynamicMesh> mesh;
    Record record{16};
    std::unique_ptr<DynamicModel> model;

    explicit LipidRun(const Param &source, double amplitude = 0.0) : param(source)
    {
        mesh.reset(new DynamicMesh(param));
        mesh->setup_flat();
        prepare(*mesh, amplitude);
        record.add(mesh->param.area, mesh->param.energy, mesh->calculate_mean_force());
        model.reset(new DynamicModel(*mesh, record));
    }

    /// The relaxed areas alone, on a net that does not move.
    void flow_step(int iteration)
    {
        model->iteration = iteration;
        model->advance_lipid_flow();
    }
};

/// Sum over faces of (K / 2 a_f) (S_f - a_f)^2 for the areas as they stand.
double stretching_energy(const Mesh &mesh)
{
    double energy = 0.0;
    for (std::size_t f = 0; f < mesh.faces.size(); f++)
    {
        const double a = mesh.faceArea0[f];
        if (a > 0.0)
        {
            const double excess = mesh.faces[f].elementArea - a;
            energy += 0.5 * mesh.param.uSurf / mesh.lipid_modulus_area(static_cast<int>(f)) * excess * excess;
        }
    }
    return energy;
}

/// A deterministic, face-dependent perturbation of the relaxed areas.
void perturb_relaxed_areas(Mesh &mesh, double amplitude)
{
    for (std::size_t f = 0; f < mesh.faces.size(); f++)
    {
        if (mesh.faceArea0[f] > 0.0)
        {
            mesh.faceArea0[f] *= 1.0 + amplitude * std::sin(0.37 * static_cast<double>(f) + 0.2);
        }
    }
}

/// The relaxed areas as a sorted multiset: what a carry must leave unchanged
/// bitwise. (Their sum need not be -- a swap reorders the summation.)
std::vector<double> sorted_relaxed_areas(const Mesh &mesh)
{
    std::vector<double> sorted = mesh.faceArea0;
    std::sort(sorted.begin(), sorted.end());
    return sorted;
}

double centroid_x(const Mesh &mesh, int f)
{
    double x = 0.0;
    for (int corner : mesh.faces[f].adjacentVertices)
    {
        x += mesh.vertices[corner].coord(0, 0) / 3.0;
    }
    return x;
}

} // namespace

TEST(LipidFlowTest, TheParameterKeysParse)
{
    Param param;
    EXPECT_EQ(param.lipidFlowMode, "off");
    EXPECT_TRUE(import_kv_string("lipidFlowMode", "local", param));
    EXPECT_EQ(param.lipidFlowMode, "local");
    EXPECT_TRUE(import_kv_string("lipidFlowDiffusivity", "120.5", param));
    EXPECT_DOUBLE_EQ(param.lipidFlowDiffusivity, 120.5);
    EXPECT_TRUE(import_kv_string("lipidFlowSubsteps", "4", param));
    EXPECT_EQ(param.lipidFlowSubsteps, 4);
    EXPECT_TRUE(import_kv_string("lipidFlowThermal", "true", param));
    EXPECT_TRUE(param.lipidFlowThermal);
    EXPECT_TRUE(import_kv_string("lipidReservoirRegion", "rim", param));
    EXPECT_EQ(param.lipidReservoirRegion, "rim");
    EXPECT_TRUE(import_kv_string("lipidReservoirTension", "0.5", param));
    EXPECT_DOUBLE_EQ(param.lipidReservoirTension, 0.5);
    EXPECT_TRUE(import_kv_string("lipidReservoirRate", "20", param));
    EXPECT_DOUBLE_EQ(param.lipidReservoirRate, 20.0);
    EXPECT_EQ(param.lipidFlipResampling, "auto");
    EXPECT_TRUE(import_kv_string("lipidFlipResampling", "patch", param));
    EXPECT_EQ(param.lipidFlipResampling, "patch");
}

TEST(LipidFlowTest, FlipResamplingDefaultsToPatchExactlyWhenFlipsCanUseIt)
{
    // auto -> patch with flips on under the linear form; none otherwise.
    struct Case
    {
        const char *mode;
        bool flips;
        bool linear;
        const char *expected;
    };
    const Case cases[] = {{"local", true, true, "patch"},
                          {"local", false, true, "none"},
                          {"local", true, false, "none"},
                          {"global", true, true, "none"},
                          {"off", true, true, "none"}};
    for (const Case &c : cases)
    {
        Param param;
        configure(param, BoundaryType::Periodic, BoundaryType::Periodic, 60.0);
        param.lipidFlowMode = c.mode;
        param.inPlaneDynamicsEnabled = c.flips;
        param.edgeSpringEnabled = c.flips;
        param.edgeTetherShape = "flat";
        param.edgeFlipEnabled = c.flips;
        param.lipidFlowLinearStretching = c.linear;
        EXPECT_EQ(param.lipidFlipResampling, "auto");
        LipidRun run(param, 0.0);
        EXPECT_EQ(run.mesh->param.lipidFlipResampling, c.expected)
            << c.mode << (c.flips ? " flips" : " no flips") << (c.linear ? " linear" : " per-lipid");
    }
}

TEST(LipidFlowTest, TheForceIsMinusTheGradientOfTheLocalEnergy)
{
  for (const bool linear : {true, false})
  {
    Param param;
    configure(param, BoundaryType::Fixed, BoundaryType::Fixed, 60.0);
    param.lipidFlowLinearStretching = linear;
    DynamicMesh mesh(param);
    mesh.setup_flat();
    prepare(mesh, 2.0);
    perturb_relaxed_areas(mesh, 0.08);
    mesh.Compute_Energy_And_Force();

    // A few vertices well inside the sheet, every axis.
    const int nx = param.nFaceX + 1;
    const std::vector<int> probes = {nx * (param.nFaceY / 2) + param.nFaceX / 2,
                                     nx * (param.nFaceY / 2 + 1) + param.nFaceX / 2 - 1,
                                     nx * (param.nFaceY / 2 - 2) + param.nFaceX / 2 + 2};
    const double eps = 1.0e-5;
    for (int v : probes)
    {
        ASSERT_TRUE(mesh.is_independent_vertex(v));
        for (int axis = 0; axis < 3; axis++)
        {
            mesh.Compute_Energy_And_Force();
            const double force = mesh.vertices[v].force.forceCurvature(axis, 0) +
                                 mesh.vertices[v].force.forceArea(axis, 0);
            const double areaForce = mesh.vertices[v].force.forceArea(axis, 0);
            const double x0 = mesh.vertices[v].coord(axis, 0);
            mesh.vertices[v].coord.set(axis, 0, x0 + eps);
            mesh.Compute_Energy_And_Force();
            const double ePlus = mesh.param.energy.energyCurvature + mesh.param.energy.energyArea;
            const double aPlus = mesh.param.energy.energyArea;
            mesh.vertices[v].coord.set(axis, 0, x0 - eps);
            mesh.Compute_Energy_And_Force();
            const double eMinus = mesh.param.energy.energyCurvature + mesh.param.energy.energyArea;
            const double aMinus = mesh.param.energy.energyArea;
            mesh.vertices[v].coord.set(axis, 0, x0);
            const double gradient = (ePlus - eMinus) / (2.0 * eps);
            const double areaGradient = (aPlus - aMinus) / (2.0 * eps);
            EXPECT_NEAR(-force, gradient, 1.0e-5 * std::max(1.0, std::abs(gradient)))
                << "vertex " << v << " axis " << axis;
            EXPECT_NEAR(-areaForce, areaGradient, 1.0e-5 * std::max(1.0, std::abs(areaGradient)))
                << "vertex " << v << " axis " << axis << " (stretching term alone)";
        }
    }
    // The stretching term must be live for this to mean anything.
    EXPECT_GT(mesh.param.energy.energyArea, 1.0);
  }
}

TEST(LipidFlowTest, UniformRelaxedAreasAreTheGlobalConstraintOnAStretchedSheet)
{
    Param local;
    configure(local, BoundaryType::Fixed, BoundaryType::Fixed, 60.0);
    Param global = local;
    global.lipidFlowMode = "off";

    DynamicMesh meshLocal(local);
    meshLocal.setup_flat();
    prepare(meshLocal, 0.0);
    DynamicMesh meshGlobal(global);
    meshGlobal.setup_flat();
    prepare(meshGlobal, 0.0);

    meshLocal.apply_box_scale({{1.01, 1.01, 1.0}});
    meshGlobal.apply_box_scale({{1.01, 1.01, 1.0}});
    meshLocal.Compute_Energy_And_Force();
    meshGlobal.Compute_Energy_And_Force();

    const double expectedTension = local.uSurf * (1.01 * 1.01 - 1.0);
    EXPECT_NEAR(meshLocal.param.energy.energyArea, meshGlobal.param.energy.energyArea,
                1.0e-9 * meshGlobal.param.energy.energyArea);
    for (std::size_t f = 0; f < meshLocal.faces.size(); f++)
    {
        if (meshLocal.faceArea0[f] > 0.0)
        {
            EXPECT_NEAR(meshLocal.faceTension[f], expectedTension, 1.0e-9);
        }
    }
    for (std::size_t v = 0; v < meshLocal.vertices.size(); v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            EXPECT_NEAR(meshLocal.vertices[v].force.forceArea(axis, 0),
                        meshGlobal.vertices[v].force.forceArea(axis, 0), 1.0e-8)
                << "vertex " << v << " axis " << axis;
        }
    }
}

TEST(LipidFlowTest, TheFluxGraphHasOneLinkPerPhysicalEdge)
{
    // A sheet that wraps both ways: every face has three neighbours, through
    // the seam included, and there is no rim.
    Param periodic;
    configure(periodic, BoundaryType::Periodic, BoundaryType::Periodic, 60.0);
    LipidRun wrapped(periodic);
    const Mesh &mesh = *wrapped.mesh;
    int nCarrying = 0;
    std::vector<int> degree(mesh.faces.size(), 0);
    for (std::size_t f = 0; f < mesh.faces.size(); f++)
    {
        nCarrying += mesh.face_carries_lipids(static_cast<int>(f)) ? 1 : 0;
        EXPECT_EQ(mesh.lipidRimFace[f], 0);
    }
    for (const Mesh::LipidLink &link : mesh.lipidLinks)
    {
        EXPECT_NE(link.faceA, link.faceB);
        EXPECT_TRUE(mesh.face_carries_lipids(link.faceA));
        EXPECT_TRUE(mesh.face_carries_lipids(link.faceB));
        degree[link.faceA]++;
        degree[link.faceB]++;
    }
    EXPECT_EQ(nCarrying, 2 * (wrapped.param.nFaceX - 6) * (wrapped.param.nFaceY - 6));
    EXPECT_EQ(static_cast<int>(mesh.lipidLinks.size()) * 2, nCarrying * 3);
    for (std::size_t f = 0; f < mesh.faces.size(); f++)
    {
        EXPECT_EQ(degree[f], mesh.face_carries_lipids(static_cast<int>(f)) ? 3 : 0) << "face " << f;
    }

    // A clamped sheet: a rim all the way round, and no face with more than
    // three neighbours.
    Param clamped;
    configure(clamped, BoundaryType::Fixed, BoundaryType::Fixed, 60.0);
    LipidRun held(clamped);
    int nRim = 0;
    for (std::size_t f = 0; f < held.mesh->faces.size(); f++)
    {
        nRim += held.mesh->lipidRimFace[f] ? 1 : 0;
    }
    EXPECT_GT(nRim, 0);
    EXPECT_LT(nRim, static_cast<int>(held.mesh->faces.size()));
}

TEST(LipidFlowTest, TheLateralFlowConservesTheTotalLowersTheEnergyAndEndsAtTheGlobalConstraint)
{
  for (const bool linear : {true, false})
  {
    Param param;
    configure(param, BoundaryType::Periodic, BoundaryType::Periodic, 60.0);
    param.lipidFlowDiffusivity = 400.0;
    param.lipidFlowLinearStretching = linear;
    LipidRun run(param, 2.0);
    Mesh &mesh = *run.mesh;
    const double referenceTotal = mesh.total_face_area0();
    perturb_relaxed_areas(mesh, 0.05);

    const double total = mesh.total_face_area0();
    double totalArea = 0.0;
    for (std::size_t f = 0; f < mesh.faces.size(); f++)
    {
        if (mesh.faceArea0[f] > 0.0)
        {
            totalArea += mesh.faces[f].elementArea;
        }
    }
    double previous = stretching_energy(mesh);
    const double start = previous;
    for (int iteration = 0; iteration < 3000; iteration++)
    {
        run.flow_step(iteration);
        if (iteration % 50 == 49)
        {
            const double now = stretching_energy(mesh);
            EXPECT_LE(now, previous + 1.0e-9 * start) << "iteration " << iteration;
            previous = now;
        }
    }
    EXPECT_NEAR(mesh.total_face_area0(), total, 1.0e-11 * total);
    EXPECT_EQ(run.model->lipidClampedUpdates, 0);

    // Uniform exchange tension, and with it the global constraint's energy.
    double low = 1.0e300;
    double high = -1.0e300;
    for (std::size_t f = 0; f < mesh.faces.size(); f++)
    {
        if (mesh.faceArea0[f] > 0.0)
        {
            const double pi = mesh.lipid_exchange_tension(static_cast<int>(f));
            low = std::min(low, pi);
            high = std::max(high, pi);
        }
    }
    EXPECT_LT(high - low, 1.0e-7);
    // (K / 2 A)(S_total - sum a)^2, with A the sum of the areas the moduli
    // are referred to: the reference areas under the linear form, the
    // relaxed areas themselves otherwise.
    const double modulusTotal = linear ? referenceTotal : total;
    const double globalEnergy = 0.5 * param.uSurf / modulusTotal * (totalArea - total) * (totalArea - total);
    EXPECT_NEAR(stretching_energy(mesh), globalEnergy, 1.0e-8 * std::max(1.0, globalEnergy))
        << (linear ? "linear" : "per-lipid");
    EXPECT_LT(stretching_energy(mesh), start);
  }
}

TEST(LipidFlowTest, ATensionWaveDecaysAtTheDiffusiveRate)
{
    Param param;
    configure(param, BoundaryType::Periodic, BoundaryType::Periodic, 100.0);
    param.lipidFlowDiffusivity = 50.0;
    LipidRun run(param, 0.0);
    Mesh &mesh = *run.mesh;
    const double period = run.mesh->boxLength[0];
    const double q = 2.0 * M_PI / period;
    const double epsilon = 0.01;

    std::vector<double> reference = mesh.faceArea0;
    for (std::size_t f = 0; f < mesh.faces.size(); f++)
    {
        if (mesh.faceArea0[f] > 0.0)
        {
            mesh.faceArea0[f] *= 1.0 + epsilon * std::cos(q * centroid_x(mesh, static_cast<int>(f)));
        }
    }
    const auto amplitude = [&]() {
        double sum = 0.0;
        int n = 0;
        for (std::size_t f = 0; f < mesh.faces.size(); f++)
        {
            if (reference[f] > 0.0)
            {
                sum += (mesh.faceArea0[f] / reference[f] - 1.0) *
                       std::cos(q * centroid_x(mesh, static_cast<int>(f)));
                n++;
            }
        }
        return 2.0 * sum / n;
    };
    const double before = amplitude();
    EXPECT_NEAR(before, epsilon, 1.0e-3 * epsilon);
    const int nSteps = 1000;
    for (int iteration = 0; iteration < nSteps; iteration++)
    {
        run.flow_step(iteration);
    }
    const double after = amplitude();
    const double measuredRate = -std::log(after / before) / (nSteps * param.timeStep);
    const double expectedRate = run.model->lipidDiffusivity * q * q;
    EXPECT_NEAR(measuredRate, expectedRate, 0.03 * expectedRate)
        << "measured " << measuredRate << " per us, D q^2 = " << expectedRate;
}

TEST(LipidFlowTest, TheThermalFlowSamplesTheRelaxedAreaVariance)
{
    Param param;
    configure(param, BoundaryType::Periodic, BoundaryType::Periodic, 60.0);
    param.lipidFlowThermal = true;
    param.lipidFlowDiffusivity = 50.0;
    LipidRun run(param, 0.0);
    Mesh &mesh = *run.mesh;
    const double total = mesh.total_face_area0();
    const int nFaces = static_cast<int>(mesh.faces.size());
    int nCarrying = 0;
    for (int f = 0; f < nFaces; f++)
    {
        nCarrying += mesh.faceArea0[f] > 0.0 ? 1 : 0;
    }

    std::vector<double> sum(nFaces, 0.0);
    std::vector<double> sumSquare(nFaces, 0.0);
    const int nSteps = 40000;
    const int nBurn = 2000;
    for (int iteration = 0; iteration < nSteps; iteration++)
    {
        run.flow_step(iteration);
        if (iteration >= nBurn)
        {
            for (int f = 0; f < nFaces; f++)
            {
                sum[f] += mesh.faceArea0[f];
                sumSquare[f] += mesh.faceArea0[f] * mesh.faceArea0[f];
            }
        }
    }
    EXPECT_NEAR(mesh.total_face_area0(), total, 1.0e-10 * total);
    EXPECT_EQ(run.model->lipidClampedUpdates, 0);

    const int n = nSteps - nBurn;
    double ratio = 0.0;
    for (int f = 0; f < nFaces; f++)
    {
        if (mesh.faceArea0Initial[f] > 0.0)
        {
            const double mean = sum[f] / n;
            const double variance = sumSquare[f] / n - mean * mean;
            ratio += variance / mesh.faceArea0Initial[f];
        }
    }
    ratio /= nCarrying;
    // kT a / K per face, less the one mode the conserved total removes.
    const double expected = param.KBT / param.uSurf * (1.0 - 1.0 / nCarrying);
    EXPECT_NEAR(ratio, expected, 0.05 * expected)
        << "var(a_f)/a_f = " << ratio << " nm^2, kT/K = " << expected;
}

TEST(LipidFlowTest, ThermalRelaxedAreasLeaveTheMeanTensionAloneOnlyUnderTheLinearForm)
{
    // A flat sheet stretched by 1%, net frozen, conserved noise on. The mean
    // relaxed area of a face cannot move, so a tension linear in a_f keeps
    // its mean exactly; a tension that goes as 1/a_f gains K <(da/a)^2> =
    // kT / a_f, which depends on the mesh.
    for (const bool linear : {true, false})
    {
        Param param;
        configure(param, BoundaryType::Periodic, BoundaryType::Periodic, 60.0);
        param.lipidFlowThermal = true;
        param.lipidFlowDiffusivity = 50.0;
        param.lipidFlowLinearStretching = linear;
        LipidRun run(param, 0.0);
        DynamicMesh &mesh = *run.mesh;
        mesh.apply_box_scale({{1.005, 1.005, 1.0}});
        mesh.Compute_Energy_And_Force();
        const double stretch = 1.005 * 1.005;
        const double quietTension = param.uSurf * (stretch - 1.0);

        int nCarrying = 0;
        double faceArea = 0.0;
        for (std::size_t f = 0; f < mesh.faces.size(); f++)
        {
            if (mesh.faceArea0[f] > 0.0)
            {
                nCarrying++;
                faceArea += mesh.faceArea0Initial[f];
            }
        }
        faceArea /= nCarrying;

        double sum = 0.0;
        int samples = 0;
        for (int iteration = 0; iteration < 40000; iteration++)
        {
            run.flow_step(iteration);
            if (iteration >= 2000 && iteration % 10 == 0)
            {
                double mean = 0.0;
                for (std::size_t f = 0; f < mesh.faces.size(); f++)
                {
                    if (mesh.faceArea0[f] > 0.0)
                    {
                        mean += mesh.lipid_face_tension(static_cast<int>(f));
                    }
                }
                sum += mean / nCarrying;
                samples++;
            }
        }
        const double meanTension = sum / samples;
        if (linear)
        {
            EXPECT_NEAR(meanTension, quietTension, 1.0e-9);
        }
        else
        {
            const double shift = stretch * param.KBT / faceArea * (1.0 - 1.0 / nCarrying);
            EXPECT_NEAR(meanTension - quietTension, shift, 0.05 * shift)
                << "mean tension " << meanTension << ", quiet " << quietTension << ", kT/a " << shift;
        }
    }
}

TEST(LipidFlowTest, ABulkReservoirBringsEveryFaceToItsTension)
{
  for (const bool linear : {true, false})
  {
    Param param;
    configure(param, BoundaryType::Fixed, BoundaryType::Fixed, 60.0);
    param.lipidFlowLinearStretching = linear;
    param.lipidFlowDiffusivity = 0.0;
    param.lipidReservoirRegion = "bulk";
    param.lipidReservoirTension = 0.5;
    param.lipidReservoirRate = 20.0;
    LipidRun run(param, 0.0);
    DynamicMesh &mesh = *run.mesh;
    mesh.apply_box_scale({{1.005, 1.005, 1.0}});
    mesh.Compute_Energy_And_Force();
    const double total = mesh.total_face_area0();

    for (int iteration = 0; iteration < 3000; iteration++)
    {
        run.flow_step(iteration);
    }
    mesh.Compute_Energy_And_Force();
    // The reservoir fixes the exchange tension: the tension itself under
    // the linear form, (K/2)(s^2 - 1) otherwise.
    const double expectedTension =
        linear ? 0.5 : param.uSurf * (std::sqrt(1.0 + 2.0 * 0.5 / param.uSurf) - 1.0);
    for (std::size_t f = 0; f < mesh.faces.size(); f++)
    {
        if (mesh.faceArea0[f] > 0.0)
        {
            EXPECT_NEAR(mesh.lipid_exchange_tension(static_cast<int>(f)), 0.5, 1.0e-9);
            EXPECT_NEAR(mesh.faceTension[f], expectedTension, 1.0e-9);
        }
    }
    // The sheet was stretched: the reservoir had to supply area.
    EXPECT_GT(run.model->lipidReservoirExchanged, 0.0);
    EXPECT_NEAR(mesh.total_face_area0() - total, run.model->lipidReservoirExchanged, 1.0e-9 * total);
  }
}

TEST(LipidFlowTest, WhatThePrototypeDoesNotCoverIsRefused)
{
    Param badMode;
    configure(badMode, BoundaryType::Periodic, BoundaryType::Periodic, 60.0);
    badMode.lipidFlowMode = "sideways";
    EXPECT_THROW(
        {
            DynamicMesh mesh(badMode);
            mesh.setup_flat();
        },
        std::runtime_error);

    Param legacy;
    configure(legacy, BoundaryType::Periodic, BoundaryType::Periodic, 60.0);
    legacy.boundaryCondition = BoundaryType::Periodic;
    EXPECT_THROW(
        {
            DynamicMesh mesh(legacy);
            mesh.setup_flat();
        },
        std::runtime_error);

    Param badRegion;
    configure(badRegion, BoundaryType::Fixed, BoundaryType::Fixed, 60.0);
    badRegion.lipidReservoirRegion = "everywhere";
    EXPECT_THROW(LipidRun run(badRegion), std::runtime_error);
}

TEST(LipidFlowTest, TheGlobalModeWithoutAReservoirIsTheGlobalConstraint)
{
    Param off;
    configure(off, BoundaryType::Periodic, BoundaryType::Periodic, 60.0);
    off.lipidFlowMode = "off";
    Param global = off;
    global.lipidFlowMode = "global";
    LipidRun a(off, 2.0);
    LipidRun b(global, 2.0);
    EXPECT_DOUBLE_EQ(a.mesh->param.energy.energyArea, b.mesh->param.energy.energyArea);
    for (std::size_t v = 0; v < a.mesh->vertices.size(); v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            EXPECT_DOUBLE_EQ(a.mesh->vertices[v].force.forceArea(axis, 0),
                             b.mesh->vertices[v].force.forceArea(axis, 0));
        }
    }
    // And the step leaves the relaxed area alone.
    const double area0 = b.mesh->param.area0;
    for (int i = 0; i < 50; i++)
    {
        b.flow_step(i);
    }
    EXPECT_DOUBLE_EQ(b.mesh->param.area0, area0);
}

TEST(LipidFlowTest, TheGlobalModeReservoirBringsTheSheetToItsTension)
{
    Param param;
    configure(param, BoundaryType::Fixed, BoundaryType::Fixed, 60.0);
    param.lipidFlowMode = "global";
    param.lipidReservoirRegion = "bulk";
    param.lipidReservoirTension = 0.5;
    param.lipidReservoirRate = 20.0;
    LipidRun run(param, 0.0);
    DynamicMesh &mesh = *run.mesh;
    mesh.apply_box_scale({{1.005, 1.005, 1.0}});
    mesh.Compute_Energy_And_Force();
    const double normalization = mesh.lipidAreaNormalization;
    EXPECT_GT(normalization, 0.0);
    for (int i = 0; i < 3000; i++)
    {
        run.flow_step(i);
        mesh.Compute_Energy_And_Force();
    }
    const double tension = param.uSurf * (mesh.param.area - mesh.param.area0) / normalization;
    EXPECT_NEAR(tension, 0.5, 1.0e-9);
    EXPECT_GT(run.model->lipidReservoirExchanged, 0.0);
    // The force carries the same tension.
    EXPECT_NEAR(mesh.param.energy.energyArea,
                0.5 * param.uSurf / normalization * std::pow(mesh.param.area - mesh.param.area0, 2),
                1.0e-9);
}

TEST(LipidFlowTest, AnEdgeFlipConservesTheRelaxedAreaAndItsTrialMatchesTheFullEnergy)
{
    // The fluid sheet: in-plane motion, the flat tether, flips admissible.
    Param param;
    configure(param, BoundaryType::Periodic, BoundaryType::Periodic, 60.0);
    param.inPlaneDynamicsEnabled = true;
    param.edgeSpringEnabled = true;
    param.edgeTetherShape = "flat";
    param.edgeTetherMinRatio = 0.7;
    param.edgeTetherMaxRatio = 1.75;
    param.edgeFlipEnabled = true;
    param.edgeFlipMinValence = 4;
    param.edgeFlipMaxValence = 8;
    param.lipidFlipResampling = "none"; // the carry alone; auto would resample
    LipidRun run(param, 2.0);
    Mesh &mesh = *run.mesh;
    perturb_relaxed_areas(mesh, 0.05);
    mesh.Compute_Energy_And_Force();

    int tested = 0;
    for (int iEdge = 0; iEdge < static_cast<int>(mesh.edges.size()) && tested < 6; iEdge++)
    {
        if (!mesh.edges[iEdge].flippable)
        {
            continue;
        }
        EdgeFlipDelta delta;
        if (!mesh.evaluate_edge_flip(iEdge, delta) || !std::isfinite(delta.energy))
        {
            continue;
        }
        const std::vector<double> multiset = sorted_relaxed_areas(mesh);
        const double before = mesh.param.energy.energyCurvature + mesh.param.energy.energyArea +
                              mesh.param.energy.energyRegularization;
        const std::vector<double> saved = mesh.faceArea0;
        const std::vector<double> savedInitial = mesh.faceArea0Initial;
        const long long version = mesh.topologyVersion;

        // Do the flip as the sweep does, and compare against a full evaluation.
        const MeshEdge edgeBefore = mesh.edges[iEdge];
        mesh.flip_edge(iEdge);
        mesh.lipid_carry_after_flip(iEdge, edgeBefore);
        EXPECT_EQ(sorted_relaxed_areas(mesh), multiset) << "carrying is bitwise conservative";
        mesh.Compute_Energy_And_Force();
        const double after = mesh.param.energy.energyCurvature + mesh.param.energy.energyArea +
                             mesh.param.energy.energyRegularization;
        EXPECT_NEAR(after - before, delta.energy, 1.0e-6 * std::max(1.0, std::abs(delta.energy)))
            << "edge " << iEdge;

        // Flip back and restore the field: the mesh is as it was.
        mesh.flip_edge(iEdge);
        mesh.faceArea0 = saved;
        mesh.faceArea0Initial = savedInitial;
        mesh.topologyVersion = version;
        mesh.Compute_Energy_And_Force();
        EXPECT_NEAR(mesh.param.energy.energyCurvature + mesh.param.energy.energyArea +
                        mesh.param.energy.energyRegularization,
                    before, 1.0e-9 * before);
        tested++;
    }
    EXPECT_GE(tested, 3);
}

TEST(LipidFlowTest, AnEdgeFlipCarriesTheFieldAsAnInvolution)
{
    // Flip an edge, then flip the new edge straight back through the same
    // code path the sweep uses. The field must come back bitwise: that is
    // what makes the Metropolis flip reversible with a per-face field. It is
    // not automatic -- flip_edge() leaves the two face indices on each
    // other's triangles after two flips -- so the carry rule's pairing is
    // what this test is really exercising, on every flippable edge.
    Param param;
    configure(param, BoundaryType::Periodic, BoundaryType::Periodic, 60.0);
    param.inPlaneDynamicsEnabled = true;
    param.edgeSpringEnabled = true;
    param.edgeTetherShape = "flat";
    param.edgeTetherMinRatio = 0.7;
    param.edgeTetherMaxRatio = 1.75;
    param.edgeFlipEnabled = true;
    param.edgeFlipMinValence = 4;
    param.edgeFlipMaxValence = 8;
    param.lipidFlipResampling = "none"; // the carry alone; auto would resample
    LipidRun run(param, 2.0);
    Mesh &mesh = *run.mesh;
    // Break the lattice symmetry so the triangle areas are all different,
    // then make the field distinguishable face by face.
    for (int v = 0; v < static_cast<int>(mesh.vertices.size()); v++)
    {
        if (!mesh.is_independent_vertex(v))
        {
            continue;
        }
        const double jitter[3] = {0.21 * std::sin(1.7 * v), 0.19 * std::cos(2.3 * v), 0.0};
        mesh.displace_vertex_with_followers(v, jitter);
    }
    perturb_relaxed_areas(mesh, 0.1);
    mesh.Compute_Energy_And_Force();

    const std::vector<double> field = mesh.faceArea0;
    const std::vector<double> fieldInitial = mesh.faceArea0Initial;
    const std::vector<double> multiset = sorted_relaxed_areas(mesh);
    int tested = 0;
    int swapped = 0;
    for (int iEdge = 0; iEdge < static_cast<int>(mesh.edges.size()); iEdge++)
    {
        if (!mesh.edges[iEdge].flippable)
        {
            continue;
        }
        const long long version = mesh.topologyVersion;

        // Which geometric triangle holds which value, before.
        const MeshEdge before = mesh.edges[iEdge];
        const double onApex0 = mesh.faceArea0[before.face[0]];
        mesh.flip_edge(iEdge);
        mesh.lipid_carry_after_flip(iEdge, before);
        EXPECT_EQ(sorted_relaxed_areas(mesh), multiset) << "edge " << iEdge;
        if (mesh.faceArea0[before.face[0]] != onApex0)
        {
            swapped++;
        }

        const MeshEdge flipped = mesh.edges[iEdge];
        mesh.flip_edge(iEdge);
        mesh.lipid_carry_after_flip(iEdge, flipped);
        mesh.topologyVersion = version;

        // The old triangle of record k had apex before.opposite[k]; find the
        // record that now holds that triangle and check it holds that value.
        for (int k = 0; k < 2; k++)
        {
            const MeshEdge &now = mesh.edges[iEdge];
            int record = -1;
            for (int slot = 0; slot < 2; slot++)
            {
                if (now.opposite[slot] == before.opposite[k])
                {
                    record = now.face[slot];
                }
            }
            ASSERT_GE(record, 0);
            EXPECT_EQ(mesh.faceArea0[record], field[before.face[k]]) << "edge " << iEdge;
            EXPECT_EQ(mesh.faceArea0Initial[record], fieldInitial[before.face[k]])
                << "edge " << iEdge;
        }
        tested++;
    }
    EXPECT_GE(tested, 50);
    // The pairing is decided by the geometry, not by the record indices: on a
    // jittered lattice some flips keep the records on their triangles and
    // some swap them.
    EXPECT_GT(swapped, 0);
    EXPECT_LT(swapped, tested);
}

TEST(LipidFlowTest, FlipResamplingChargesThePooledEnergyAndLeavesTheSetAtOneTension)
{
    // The compound move: flip, then redraw the patch's relaxed areas. The
    // trial must charge the patch's stretching energy minimised over its
    // split, and without noise the redraw is that minimiser, so after the
    // move every face of the patch sits at one tension and the full energy's
    // area term equals what the trial priced.
    for (const char *resampling : {"patch", "pair"})
    {
        Param param;
        configure(param, BoundaryType::Periodic, BoundaryType::Periodic, 60.0);
        param.inPlaneDynamicsEnabled = true;
        param.edgeSpringEnabled = true;
        param.edgeTetherShape = "flat";
        param.edgeTetherMinRatio = 0.7;
        param.edgeTetherMaxRatio = 1.75;
        param.edgeFlipEnabled = true;
        param.edgeFlipMinValence = 4;
        param.edgeFlipMaxValence = 8;
        param.lipidFlipResampling = resampling;
        LipidRun run(param, 2.0);
        Mesh &mesh = *run.mesh;
        perturb_relaxed_areas(mesh, 0.1);
        mesh.Compute_Energy_And_Force();

        int tested = 0;
        for (int iEdge = 0; iEdge < static_cast<int>(mesh.edges.size()) && tested < 4; iEdge++)
        {
            if (!mesh.edges[iEdge].flippable)
            {
                continue;
            }
            EdgeFlipDelta delta;
            if (!mesh.evaluate_edge_flip(iEdge, delta) || !std::isfinite(delta.energy))
            {
                continue;
            }
            const std::vector<int> set = mesh.lipid_resample_set(iEdge);
            const std::vector<int> patch = mesh.flip_patch_faces(iEdge);
            const double total = mesh.total_face_area0();

            // The energy the trial should have priced for the state before:
            // the set pooled at one tension, the rest of the patch as it is.
            double pooledArea = 0.0, pooledRelaxed = 0.0, pooledModulus = 0.0, rest = 0.0;
            for (int f : patch)
            {
                if (std::find(set.begin(), set.end(), f) != set.end())
                {
                    pooledArea += mesh.faces[f].elementArea;
                    pooledRelaxed += mesh.faceArea0[f];
                    pooledModulus += mesh.lipid_modulus_area(f);
                }
                else
                {
                    rest += mesh.lipid_face_stretching_energy(f);
                }
            }
            const double pricedBefore =
                rest + 0.5 * param.uSurf / pooledModulus * std::pow(pooledArea - pooledRelaxed, 2);
            // The area energy of the faces the flip does not touch.
            double outsidePatch = 0.0;
            for (int f = 0; f < static_cast<int>(mesh.faces.size()); f++)
            {
                if (std::find(patch.begin(), patch.end(), f) == patch.end())
                {
                    outsidePatch += mesh.lipid_face_stretching_energy(f);
                }
            }

            // Do the move as the sweep does.
            const std::vector<double> saved = mesh.faceArea0;
            const long long version = mesh.topologyVersion;
            const MeshEdge before = mesh.edges[iEdge];
            mesh.flip_edge(iEdge);
            mesh.lipid_carry_after_flip(iEdge, before);
            std::vector<double> areas;
            mesh.evaluate_face_subset(set, &areas);
            mesh.lipid_resample_after_flip(set, areas, 12345u);
            EXPECT_NEAR(mesh.total_face_area0(), total, 1.0e-10 * total);
            mesh.Compute_Energy_And_Force();

            // One tension across the set.
            const double tension0 = mesh.lipid_face_tension(set[0]);
            for (int f : set)
            {
                EXPECT_NEAR(mesh.lipid_face_tension(f), tension0, 1.0e-8)
                    << resampling << " edge " << iEdge;
            }
            // The full area energy after is outside + rest' + the pooled
            // minimum, and the trial priced (rest' + pooled') - pricedBefore.
            EXPECT_NEAR(mesh.param.energy.energyArea - outsidePatch - pricedBefore,
                        delta.areaConstraint, 1.0e-6 * std::max(1.0, std::abs(delta.energy)))
                << resampling << " edge " << iEdge;

            mesh.flip_edge(iEdge);
            mesh.faceArea0 = saved;
            mesh.topologyVersion = version;
            mesh.Compute_Energy_And_Force();
            tested++;
        }
        EXPECT_GE(tested, 2) << resampling;
    }
}

TEST(LipidFlowTest, ThermalFlipResamplingKeepsTheTotalAndSpreadsTheTension)
{
    Param param;
    configure(param, BoundaryType::Periodic, BoundaryType::Periodic, 60.0);
    param.inPlaneDynamicsEnabled = true;
    param.edgeSpringEnabled = true;
    param.edgeTetherShape = "flat";
    param.edgeFlipEnabled = true;
    param.edgeFlipMinValence = 4;
    param.edgeFlipMaxValence = 8;
    param.lipidFlowThermal = true;
    param.lipidFlipResampling = "patch";
    LipidRun run(param, 2.0);
    Mesh &mesh = *run.mesh;
    mesh.Compute_Energy_And_Force();
    int iEdge = 0;
    while (!mesh.edges[iEdge].flippable)
    {
        iEdge++;
    }
    const std::vector<int> set = mesh.lipid_resample_set(iEdge);
    const double total = mesh.total_face_area0();
    const MeshEdge before = mesh.edges[iEdge];
    mesh.flip_edge(iEdge);
    mesh.lipid_carry_after_flip(iEdge, before);
    std::vector<double> areas;
    mesh.evaluate_face_subset(set, &areas);
    mesh.lipid_resample_after_flip(set, areas, 777u);
    EXPECT_NEAR(mesh.total_face_area0(), total, 1.0e-10 * total);
    mesh.Compute_Energy_And_Force();
    // The draw has the variance kT c_f / K about the equal-tension split, so
    // the tensions scatter by about sqrt(kT K / c): 10 pN/nm here.
    double mean = 0.0, meanSq = 0.0;
    for (int f : set)
    {
        const double t = mesh.lipid_face_tension(f);
        mean += t / set.size();
        meanSq += t * t / set.size();
    }
    const double spread = std::sqrt(std::max(0.0, meanSq - mean * mean));
    EXPECT_GT(spread, 2.0);
    EXPECT_LT(spread, 40.0);
    // Reproducible: the same key gives the same draw (to rounding of the
    // pooled total, which is re-summed from the redrawn values).
    const std::vector<double> drawn = mesh.faceArea0;
    mesh.lipid_resample_after_flip(set, areas, 777u);
    for (int f : set)
    {
        EXPECT_NEAR(mesh.faceArea0[f], drawn[f], 1.0e-9);
    }
}

TEST(LipidFlowTest, AVertexMoveTrialMatchesTheFullEnergyUnderLocalFlow)
{
    Param param;
    configure(param, BoundaryType::Periodic, BoundaryType::Periodic, 60.0);
    param.inPlaneDynamicsEnabled = true;
    param.edgeSpringEnabled = true;
    param.edgeTetherShape = "flat";
    param.vertexMoveEnabled = true;
    param.brownianStepEnabled = false;
    LipidRun run(param, 2.0);
    Mesh &mesh = *run.mesh;
    perturb_relaxed_areas(mesh, 0.05);
    mesh.Compute_Energy_And_Force();
    int tested = 0;
    for (int v = 0; v < static_cast<int>(mesh.vertices.size()) && tested < 5; v += 7)
    {
        if (!mesh.is_independent_vertex(v))
        {
            continue;
        }
        const double displacement[3] = {0.3, -0.2, 0.4};
        const double before = mesh.param.energy.energyCurvature + mesh.param.energy.energyArea +
                              mesh.param.energy.energyRegularization;
        const VertexMoveDelta delta = mesh.evaluate_vertex_move(v, displacement);
        if (delta.wallViolated || !std::isfinite(delta.energy))
        {
            continue;
        }
        const double original[3] = {mesh.vertices[v].coord(0, 0), mesh.vertices[v].coord(1, 0),
                                    mesh.vertices[v].coord(2, 0)};
        mesh.displace_vertex_with_followers(v, displacement);
        mesh.Compute_Energy_And_Force();
        const double after = mesh.param.energy.energyCurvature + mesh.param.energy.energyArea +
                             mesh.param.energy.energyRegularization;
        EXPECT_NEAR(after - before, delta.energy, 1.0e-6 * std::max(1.0, std::abs(delta.energy)))
            << "vertex " << v;
        for (int axis = 0; axis < 3; axis++)
        {
            mesh.vertices[v].coord.set(axis, 0, original[axis]);
        }
        mesh.sync_periodic_images();
        mesh.Compute_Energy_And_Force();
        tested++;
    }
    EXPECT_GE(tested, 3);
}
