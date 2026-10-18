#include "test_divergence_guard.hpp"

/**
 * A diverging run does not always reach NaN: diag_fixedbox_dt0.005 (the frozen
 * flat sheet at timeStep 0.005 with the elastic box off) went to E_total ~
 * 2e206 pN.nm with coordinates ~ 1e48 nm and ran on, finite, to maxIterations.
 * The guard stops it on the coordinates or on the energy's size. These tests
 * check each sign on a mesh made to show it, that a healthy sheet shows none,
 * and that a whole run at a time step far past the stable one stops at the
 * step it goes, with its outputs written.
 *
 * @see docs/resolution_limits.md section 4
 */

namespace
{

/// A small periodic flat sheet, set up and evaluated the way the driver does.
void configure_sheet(Param &param)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Periodic;
    param.sideX = 40.0;
    param.sideY = 40.0;
    param.lFace = 5.0;
    param.kCurv = 83.4;
    param.KBT = 4.17;
    param.uSurf = 250.0;
    param.uVol = 0.0;
    param.isGlobalConstraint = true;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
}

struct Sheet
{
    Param param;
    std::unique_ptr<DynamicMesh> mesh;

    Sheet()
    {
        // The mesh keeps a reference to param, so param is configured
        // first and the mesh built on it afterwards.
        configure_sheet(param);
        mesh.reset(new DynamicMesh(param));
        mesh->setup_flat();
        for (Vertex &vertex : mesh->vertices)
        {
            vertex.coord.set(2, 0, 10.0);
        }
        mesh->calculate_element_area_volume();
        mesh->sum_membrane_area_and_volume(mesh->param.area0, mesh->param.vol0);
        mesh->update_previous_coord_for_vertex();
        mesh->update_reference_coord_from_previous_coord();
        mesh->Compute_Energy_And_Force();
    }
};

/// The run's working directory for the length of a test: run_dynamics_flat()
/// writes fixed-name outputs into the current directory.
class ScratchDirectory
{
public:
    ScratchDirectory()
    {
        char cwd[4096];
        if (getcwd(cwd, sizeof cwd) == nullptr)
        {
            throw std::runtime_error("getcwd failed");
        }
        previous_ = cwd;
        const char *tmp = std::getenv("TMPDIR");
        std::string pattern = std::string(tmp ? tmp : "/tmp") + "/slimed_divergence_XXXXXX";
        std::vector<char> buffer(pattern.begin(), pattern.end());
        buffer.push_back('\0');
        if (mkdtemp(buffer.data()) == nullptr)
        {
            throw std::runtime_error("mkdtemp failed");
        }
        path_ = buffer.data();
        if (chdir(path_.c_str()) != 0)
        {
            throw std::runtime_error("chdir failed");
        }
    }

    ~ScratchDirectory()
    {
        if (chdir(previous_.c_str()) == 0)
        {
            std::system(("rm -rf '" + path_ + "'").c_str());
        }
    }

private:
    std::string previous_;
    std::string path_;
};

/// The rows of EnergyForce.csv after its header, as E_total per row.
std::vector<double> read_total_energies(const std::string &path)
{
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);   // header
    std::vector<double> totals;
    while (std::getline(in, line))
    {
        std::stringstream row(line);
        std::string cell;
        std::vector<std::string> cells;
        while (std::getline(row, cell, ','))
        {
            cells.push_back(cell);
        }
        if (cells.size() >= 7)
        {
            totals.push_back(std::stod(cells[6]));
        }
    }
    return totals;
}

} // namespace

TEST(DivergenceGuardTest, AHealthySheetPasses)
{
    Sheet sheet;
    const DivergenceGuard guard(*sheet.mesh);
    EXPECT_EQ(guard.check(*sheet.mesh), "");

    // Its bounds sit ten mesh sizes out on every axis, and the energy bound
    // is far above a thermal sheet's, which starts at ~ 0 here.
    EXPECT_LT(guard.lower()[0], -300.0);
    EXPECT_GT(guard.upper()[2], 300.0);
    EXPECT_GE(guard.energy_bound(), DivergenceGuard::kEnergyMultiple * 4.17 *
                                        static_cast<double>(sheet.mesh->vertices.size()));

    // A large but physical excursion -- a vertex 50 nm off the sheet -- is
    // not a divergence.
    sheet.mesh->vertices[7].coord.set(2, 0, 60.0);
    EXPECT_EQ(guard.check(*sheet.mesh), "");
}

TEST(DivergenceGuardTest, AFiniteCoordinateFarOutsideTheMeshStopsTheRun)
{
    Sheet sheet;
    const DivergenceGuard guard(*sheet.mesh);
    // The state diag_fixedbox_dt0.005 kept for 80k steps.
    sheet.mesh->vertices[3].coord.set(0, 0, 1.0e48);
    const std::string reason = guard.check(*sheet.mesh);
    EXPECT_NE(reason.find("Vertex 3 has x = 1e+48"), std::string::npos) << reason;
}

TEST(DivergenceGuardTest, ANonFiniteCoordinateStopsTheRun)
{
    Sheet sheet;
    const DivergenceGuard guard(*sheet.mesh);
    sheet.mesh->vertices[5].coord.set(1, 0, std::numeric_limits<double>::quiet_NaN());
    EXPECT_NE(guard.check(*sheet.mesh).find("Vertex 5 has y"), std::string::npos);
}

TEST(DivergenceGuardTest, AHugeFiniteEnergyStopsTheRun)
{
    Sheet sheet;
    const DivergenceGuard guard(*sheet.mesh);
    sheet.mesh->param.energy.energyTotal = 2.0e206;
    const std::string reason = guard.check(*sheet.mesh);
    EXPECT_NE(reason.find("The total energy is 2e+206 pN.nm, beyond"), std::string::npos) << reason;

    // Just inside the bound is still a run.
    sheet.mesh->param.energy.energyTotal = 0.5 * guard.energy_bound();
    EXPECT_EQ(guard.check(*sheet.mesh), "");
}

TEST(DivergenceGuardTest, ANonFiniteEnergyStopsTheRunWithTheOriginalMessage)
{
    Sheet sheet;
    const DivergenceGuard guard(*sheet.mesh);
    sheet.mesh->param.energy.energyTotal = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(guard.check(*sheet.mesh).rfind("The total energy is ", 0), 0u);
    EXPECT_EQ(guard.check(*sheet.mesh).find("beyond"), std::string::npos);
}

TEST(DivergenceGuardTest, TheEnergyBoundScalesWithALargeStartingEnergy)
{
    Sheet sheet;
    sheet.mesh->param.energy.energyTotal = 1.0e9;
    const DivergenceGuard guard(*sheet.mesh);
    EXPECT_DOUBLE_EQ(guard.energy_bound(), DivergenceGuard::kEnergyMultiple * 1.0e9);
}

/**
 * @brief The whole driver on a run that blows up: the frozen periodic sheet of
 * the resolution sweep, shrunk to a 30 nm tile, at timeStep 0.02 -- eight
 * times the stable step. Before the guard this run went 1.7e3 -> 1e7 -> 6e20
 * -> 3e57 -> 1e165 and only stopped at the NaN after; now it stops at the
 * first step whose coordinates leave the bound, with the energy still finite
 * and every output written.
 */
TEST(DivergenceGuardTest, ARunPastTheStableStepStopsAtTheStepItGoes)
{
    ScratchDirectory scratch;
    {
        std::ofstream params("input.params");
        params << "boundaryType = Mixed\n"
                  "boundaryTypeX = Periodic\n"
                  "boundaryTypeY = Periodic\n"
                  "meanHeightPinned = true\n"
                  "surfaceSolver = iterative\n"
                  "lFace = 5.0\n"
                  "isFlat = true\n"
                  "sideX = 60.0\n"
                  "sideY = 51.961524\n"
                  "c0Membrane = 0.0\n"
                  "kcMembraneBending = 83.4\n"
                  "usMembraneStretching = 250.0\n"
                  "setRelaxAreaToDefault = true\n"
                  "uvVolumeConstraint = 0.0\n"
                  "isGlobalConstraint = true\n"
                  "KBT = 4.17\n"
                  "maxIterations = 200\n"
                  "timeStep = 0.02\n"
                  "diffConst = 1.0\n"
                  "randomSeed = 20261012\n"
                  "meshpointOutput = false\n"
                  "inPlaneDynamicsEnabled = true\n"
                  "edgeSpringEnabled = true\n"
                  "edgeTetherShape = flat\n"
                  "edgeTetherMinRatio = 0.7\n"
                  "edgeTetherMaxRatio = 1.75\n"
                  "edgeSpringConstant = 83.4\n"
                  "triangleShapeEnabled = true\n"
                  "triangleShapeMinAltitudeRatio = 0.4\n"
                  "triangleShapeConstant = 83.4\n"
                  "creaseWallEnabled = true\n"
                  "creaseWallAngle = 60.0\n"
                  "creaseWallConstant = 500.0\n"
                  "edgeFlipEnabled = false\n"
                  "elasticBoxEnabled = false\n"
                  "isEnergyHarmonicBondIncluded = false\n"
                  "VERBOSE_MODE = false\n";
    }

    testing::internal::CaptureStderr();
    testing::internal::CaptureStdout();
    run_dynamics_flat("./input");
    testing::internal::GetCapturedStdout();
    const std::string err = testing::internal::GetCapturedStderr();

    // Stopped by the new criteria, not by a NaN.
    EXPECT_NE(err.find("the run has diverged and is stopping"), std::string::npos) << err;
    EXPECT_EQ(err.find("The total energy is nan"), std::string::npos) << err;
    EXPECT_EQ(err.find("The total energy is -nan"), std::string::npos) << err;

    // Everything up to the stop is written, and nothing in it is garbage.
    const std::vector<double> totals = read_total_energies("EnergyForce.csv");
    ASSERT_GE(totals.size(), 2u);
    EXPECT_LT(totals.size(), 12u) << "the run went on after it blew up";
    for (double total : totals)
    {
        EXPECT_TRUE(std::isfinite(total));
        EXPECT_LT(std::fabs(total), 1.0e30);
    }
    EXPECT_TRUE(std::ifstream("inputvertex_final.csv").good());
    EXPECT_TRUE(std::ifstream("inputvertex_type_final.csv").good());
}
