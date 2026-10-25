/**
 * @file test_run_drivers.cpp
 * @brief End-to-end tests of run_flat() and run_dynamics_flat() on tiny sheets.
 *
 * Each test writes a small parameter file into a fresh temporary directory,
 * runs one of the drivers there, and checks the outputs it promises: that the
 * files exist, that their row counts follow from the iteration count and the
 * mesh size, and that the energies are finite (and, for the minimizer, that
 * the energy went down).
 */
#include "test_run_drivers.hpp"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <ftw.h>
#include <limits.h>
#include <unistd.h>

namespace
{

int remove_entry(const char *path, const struct stat *, int, struct FTW *)
{
    return std::remove(path);
}

/// A fresh temporary directory that is the working directory for the
/// lifetime of the object. The destructor restores the previous working
/// directory and deletes the temporary one with everything in it.
class ScratchRunDirectory
{
public:
    ScratchRunDirectory()
    {
        char cwd[PATH_MAX];
        if (getcwd(cwd, sizeof(cwd)) == nullptr)
        {
            throw std::runtime_error("getcwd failed");
        }
        originalDirectory_ = cwd;

        const char *tmpdir = std::getenv("TMPDIR");
        std::string base = (tmpdir != nullptr && tmpdir[0] != '\0') ? tmpdir : "/tmp";
        if (base.back() != '/')
        {
            base += '/';
        }
        std::string pattern = base + "slimed_run_drivers_XXXXXX";
        std::vector<char> buffer(pattern.begin(), pattern.end());
        buffer.push_back('\0');
        if (mkdtemp(buffer.data()) == nullptr)
        {
            throw std::runtime_error(std::string("mkdtemp failed: ") + std::strerror(errno));
        }
        path_ = buffer.data();
        if (chdir(path_.c_str()) != 0)
        {
            throw std::runtime_error("chdir into the scratch directory failed");
        }
    }

    ~ScratchRunDirectory()
    {
        if (chdir(originalDirectory_.c_str()) != 0)
        {
            ADD_FAILURE() << "could not restore the working directory " << originalDirectory_;
        }
        nftw(path_.c_str(), remove_entry, 16, FTW_DEPTH | FTW_PHYS);
    }

    ScratchRunDirectory(const ScratchRunDirectory &) = delete;
    ScratchRunDirectory &operator=(const ScratchRunDirectory &) = delete;

    const std::string &path() const { return path_; }

private:
    std::string originalDirectory_;
    std::string path_;
};

/// An absolute path to a file in the repository: the tests run with the
/// repository root as their working directory, and the drivers run in a
/// scratch directory, so call this before creating the ScratchRunDirectory.
std::string repo_path(const std::string &relative)
{
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) == nullptr)
    {
        throw std::runtime_error("getcwd failed");
    }
    return std::string(cwd) + "/" + relative;
}

void write_params(const std::string &name, const std::vector<std::string> &lines)
{
    std::ofstream out(name + ".params");
    for (const std::string &line : lines)
    {
        out << line << '\n';
    }
}

bool file_exists(const std::string &path)
{
    std::ifstream in(path);
    return in.good();
}

int count_lines(const std::string &path)
{
    std::ifstream in(path);
    int n = 0;
    std::string line;
    while (std::getline(in, line))
    {
        if (!line.empty())
        {
            n++;
        }
    }
    return n;
}

/// Every line of a comma-separated numeric file, after skipping
/// `headerLines` lines. Trailing commas (the trajectory layout) are ignored.
std::vector<std::vector<double>> read_numeric_csv(const std::string &path, int headerLines)
{
    std::ifstream in(path);
    std::vector<std::vector<double>> rows;
    std::string line;
    int lineNumber = 0;
    while (std::getline(in, line))
    {
        if (lineNumber++ < headerLines || line.empty())
        {
            continue;
        }
        std::vector<double> row;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ','))
        {
            if (cell.find_first_not_of(" \t\r") == std::string::npos)
            {
                continue;
            }
            row.push_back(std::strtod(cell.c_str(), nullptr));
        }
        rows.push_back(row);
    }
    return rows;
}

/// E_Total is the seventh column of EnergyForce.csv, Mean Force the eighth.
constexpr int kEnergyTotalColumn = 6;
constexpr int kMeanForceColumn = 7;
constexpr int kEnergyForceColumns = 8;

void expect_all_finite(const std::vector<std::vector<double>> &rows, const std::string &what)
{
    for (std::size_t r = 0; r < rows.size(); r++)
    {
        for (std::size_t c = 0; c < rows[r].size(); c++)
        {
            ASSERT_TRUE(std::isfinite(rows[r][c])) << what << " row " << r << " column " << c;
        }
    }
}

/// The parameters of a small clamped sheet whose relaxed area is below the
/// area it is built with, so the minimizer has a non-zero force to follow.
std::vector<std::string> clamped_sheet_minimization(int maxIterations)
{
    return {
        "maxIterations = " + std::to_string(maxIterations),
        "VERBOSE_MODE = false",
        "boundaryType = Fixed",
        "lFace = 5.0",
        "sideX = 40.0",
        "sideY = 40.0",
        "setRelaxAreaToDefault = false",
        "relaxArea = 900.0",
        "isEnergyHarmonicBondIncluded = false",
        "randomSeed = 11",
    };
}

std::vector<std::string> with(std::vector<std::string> base, const std::vector<std::string> &extra)
{
    base.insert(base.end(), extra.begin(), extra.end());
    return base;
}

} // namespace

// ---------------------------------------------------------------------------
// run_flat
// ---------------------------------------------------------------------------

TEST(RunDriversFlatTest, MinimizationWritesItsOutputsAndLowersTheEnergy)
{
    ScratchRunDirectory dir;
    write_params("flat", with(clamped_sheet_minimization(5),
                              {"meshpointOutput = true",
                               "checkpointOutputInterval = 2",
                               "checkpointOutputFile = flat.chk"}));
    run_flat("flat");

    for (const char *file : {"face.csv", "vertex_begin.csv", "vertex_initial_step_0.csv",
                             "vertex_initial_step_1.csv", "vertex0.csv", "EnergyForce.csv",
                             "ElementFaceEnergy.csv", "vertexfinal.csv", "flat.chk"})
    {
        EXPECT_TRUE(file_exists(file)) << file;
    }
    // Only the Mixed boundary exports the typed mesh; no thermal moves ran.
    EXPECT_FALSE(file_exists("mesh_vertices.csv"));
    EXPECT_FALSE(file_exists("ThermalFluctuation.csv"));

    const int nVertices = count_lines("vertex_begin.csv");
    const int nFaces = count_lines("face.csv");
    ASSERT_GT(nVertices, 0);
    ASSERT_GT(nFaces, 0);
    EXPECT_EQ(count_lines("vertexfinal.csv"), nVertices);
    EXPECT_EQ(count_lines("vertex_initial_step_0.csv"), nVertices);
    EXPECT_EQ(count_lines("ElementFaceEnergy.csv"), nFaces + 1); // plus the header

    // maxIterations = 5: the loop runs while iteration < 4, so four NCG steps
    // and one row per iteration 0..4.
    const auto energy = read_numeric_csv("EnergyForce.csv", 1);
    ASSERT_EQ(energy.size(), 5u);
    for (const auto &row : energy)
    {
        ASSERT_EQ(row.size(), static_cast<std::size_t>(kEnergyForceColumns));
    }
    expect_all_finite(energy, "EnergyForce.csv");
    EXPECT_GT(energy.front()[kEnergyTotalColumn], 0.0);
    EXPECT_LT(energy.back()[kEnergyTotalColumn], energy.front()[kEnergyTotalColumn]);
    EXPECT_GT(energy.front()[kMeanForceColumn], 0.0);

    // The final structure has moved off the initial one.
    const auto begin = read_numeric_csv("vertex_begin.csv", 0);
    const auto final = read_numeric_csv("vertexfinal.csv", 0);
    ASSERT_EQ(begin.size(), final.size());
    double maxDisplacement = 0.0;
    for (std::size_t i = 0; i < begin.size(); i++)
    {
        for (int k = 0; k < 3; k++)
        {
            maxDisplacement = std::max(maxDisplacement, std::abs(final[i][k] - begin[i][k]));
        }
    }
    EXPECT_GT(maxDisplacement, 0.0);
}

TEST(RunDriversFlatTest, RestartContinuesFromTheCheckpoint)
{
    ScratchRunDirectory dir;
    write_params("first", with(clamped_sheet_minimization(3),
                               {"meshpointOutput = false",
                                "checkpointOutputInterval = 2",
                                "checkpointOutputFile = restart.chk"}));
    run_flat("first");
    ASSERT_TRUE(file_exists("restart.chk"));
    const auto firstEnergy = read_numeric_csv("EnergyForce.csv", 1);
    ASSERT_EQ(firstEnergy.size(), 3u); // iterations 0..2

    // A restarted run does not rewrite the initial structure.
    ASSERT_EQ(std::remove("face.csv"), 0);
    ASSERT_EQ(std::remove("vertex_begin.csv"), 0);

    write_params("second", with(clamped_sheet_minimization(6),
                                {"meshpointOutput = false",
                                 "restartInputFile = restart.chk",
                                 "checkpointOutputInterval = 0"}));
    run_flat("second");

    EXPECT_FALSE(file_exists("face.csv"));
    EXPECT_FALSE(file_exists("vertex_begin.csv"));
    const auto secondEnergy = read_numeric_csv("EnergyForce.csv", 1);
    // The record came back from the checkpoint (iterations 0..2) and the run
    // went on to iteration 5.
    ASSERT_EQ(secondEnergy.size(), 6u);
    expect_all_finite(secondEnergy, "EnergyForce.csv after restart");
    for (std::size_t r = 0; r < firstEnergy.size(); r++)
    {
        for (int c = 0; c < kEnergyForceColumns; c++)
        {
            const double expected = firstEnergy[r][c];
            EXPECT_NEAR(secondEnergy[r][c], expected, 1e-5 * std::max(1.0, std::abs(expected)))
                << "row " << r << " column " << c;
        }
    }
    EXPECT_LT(secondEnergy.back()[kEnergyTotalColumn], secondEnergy.front()[kEnergyTotalColumn]);
}

TEST(RunDriversFlatTest, AMissingCheckpointStopsTheRunBeforeAnyOutput)
{
    ScratchRunDirectory dir;
    write_params("broken", with(clamped_sheet_minimization(4),
                                {"restartInputFile = does_not_exist.chk",
                                 "checkpointOutputInterval = 0"}));
    run_flat("broken");
    EXPECT_FALSE(file_exists("EnergyForce.csv"));
    EXPECT_FALSE(file_exists("vertexfinal.csv"));
    EXPECT_FALSE(file_exists("face.csv"));
}

TEST(RunDriversFlatTest, PureMetropolisMonteCarloRecordsEveryTrial)
{
    ScratchRunDirectory dir;
    write_params("mmc", with(clamped_sheet_minimization(4),
                             {"meshpointOutput = true",
                              "checkpointOutputInterval = 1",
                              "checkpointOutputFile = mmc.chk",
                              "thermalFluctuationEnabled = true",
                              "thermalFluctuationPureMMC = true",
                              "thermalFluctuationStepScale = 0.02"}));
    run_flat("mmc");

    // maxIterations = 4: one trial per iteration while iteration < 3.
    const auto thermal = read_numeric_csv("ThermalFluctuation.csv", 1);
    ASSERT_EQ(thermal.size(), 3u);
    for (std::size_t r = 0; r < thermal.size(); r++)
    {
        ASSERT_EQ(thermal[r].size(), 6u);
        EXPECT_EQ(thermal[r][0], static_cast<double>(r));    // iteration
        EXPECT_GT(thermal[r][2], 0.0);                       // kBT
        EXPECT_GE(thermal[r][4], 0.0);                       // acceptance probability
        EXPECT_LE(thermal[r][4], 1.0);
        EXPECT_TRUE(thermal[r][5] == 0.0 || thermal[r][5] == 1.0);
    }
    const auto energy = read_numeric_csv("EnergyForce.csv", 1);
    ASSERT_EQ(energy.size(), 4u);
    expect_all_finite(energy, "EnergyForce.csv");
    for (const char *file : {"vertex_initial_step_0.csv", "vertex_initial_step_1.csv",
                             "vertex_initial_step_3.csv", "mmc.chk", "vertexfinal.csv"})
    {
        EXPECT_TRUE(file_exists(file)) << file;
    }
}

TEST(RunDriversFlatTest, PureMetropolisMonteCarloNeedsThermalFluctuationsEnabled)
{
    ScratchRunDirectory dir;
    write_params("mmc_off", with(clamped_sheet_minimization(4),
                                 {"meshpointOutput = false",
                                  "checkpointOutputInterval = 0",
                                  "thermalFluctuationEnabled = false",
                                  "thermalFluctuationPureMMC = true"}));
    run_flat("mmc_off");
    // The initial structure is written, then the run refuses to start.
    EXPECT_TRUE(file_exists("face.csv"));
    EXPECT_FALSE(file_exists("EnergyForce.csv"));
    EXPECT_FALSE(file_exists("ThermalFluctuation.csv"));
}

TEST(RunDriversFlatTest, PureMetropolisMonteCarloStopsWhenNoTrialCanBeMade)
{
    ScratchRunDirectory dir;
    // A zero step scale makes every trial impossible.
    write_params("mmc_zero", with(clamped_sheet_minimization(4),
                                  {"meshpointOutput = false",
                                   "checkpointOutputInterval = 0",
                                   "thermalFluctuationEnabled = true",
                                   "thermalFluctuationPureMMC = true",
                                   "thermalFluctuationStepScale = 0.0"}));
    run_flat("mmc_zero");
    const auto energy = read_numeric_csv("EnergyForce.csv", 1);
    EXPECT_EQ(energy.size(), 1u); // only the starting configuration
    EXPECT_FALSE(file_exists("ThermalFluctuation.csv"));
    EXPECT_TRUE(file_exists("vertexfinal.csv"));
}

TEST(RunDriversFlatTest, ThermalTrialStandsInWhenThereIsNoForceToFollow)
{
    ScratchRunDirectory dir;
    // A flat periodic sheet at its relaxed area has no force, so the first
    // trial step size is negative and the thermal trial is what moves it.
    write_params("thermal", {"maxIterations = 4",
                             "VERBOSE_MODE = false",
                             "boundaryType = Periodic",
                             "lFace = 5.0",
                             "sideX = 40.0",
                             "sideY = 40.0",
                             "setRelaxAreaToDefault = true",
                             "meshpointOutput = false",
                             "checkpointOutputInterval = 0",
                             "randomSeed = 5",
                             "thermalFluctuationEnabled = true",
                             "thermalFluctuationInterval = 1"});
    run_flat("thermal");
    const auto thermal = read_numeric_csv("ThermalFluctuation.csv", 1);
    ASSERT_GE(thermal.size(), 1u);
    EXPECT_EQ(thermal.front()[0], 0.0);
    const auto energy = read_numeric_csv("EnergyForce.csv", 1);
    ASSERT_GE(energy.size(), 1u);
    expect_all_finite(energy, "EnergyForce.csv");
    EXPECT_EQ(energy.front()[kEnergyTotalColumn], 0.0);
    // At room temperature the first trial is a coin toss (the random streams
    // differ between standard libraries). Rejected, it leaves no force to
    // follow and the run stops with only the starting row.
    if (thermal.front()[5] == 0.0)
    {
        EXPECT_EQ(thermal.size(), 1u);
        EXPECT_EQ(energy.size(), 1u);
    }
}

TEST(RunDriversFlatTest, AcceptedThermalTrialsHandTheRunBackToTheMinimizer)
{
    ScratchRunDirectory dir;
    // The same force-free sheet, at a temperature so high that every trial is
    // accepted: the first one stands in for the missing NCG step, after which
    // there is a force and the minimizer takes over, with a trial after each
    // of its steps.
    write_params("hot", {"maxIterations = 4",
                         "VERBOSE_MODE = false",
                         "boundaryType = Periodic",
                         "lFace = 5.0",
                         "sideX = 40.0",
                         "sideY = 40.0",
                         "setRelaxAreaToDefault = true",
                         "meshpointOutput = false",
                         "checkpointOutputInterval = 0",
                         "randomSeed = 5",
                         "thermalFluctuationEnabled = true",
                         "thermalFluctuationInterval = 1",
                         "thermalFluctuationTemperatureKelvin = 1e12",
                         "thermalFluctuationMinTemperatureKelvin = 1e12"});
    run_flat("hot");
    const auto thermal = read_numeric_csv("ThermalFluctuation.csv", 1);
    ASSERT_EQ(thermal.size(), 3u);
    for (std::size_t r = 0; r < thermal.size(); r++)
    {
        EXPECT_EQ(thermal[r][0], static_cast<double>(r));
        EXPECT_DOUBLE_EQ(thermal[r][1], 1e12);
        EXPECT_EQ(thermal[r][5], 1.0); // accepted
    }
    const auto energy = read_numeric_csv("EnergyForce.csv", 1);
    ASSERT_EQ(energy.size(), 4u);
    expect_all_finite(energy, "EnergyForce.csv");
    EXPECT_EQ(energy.front()[kEnergyTotalColumn], 0.0);
    // The accepted trial left the flat state, so the sheet now has energy.
    EXPECT_GT(energy[1][kEnergyTotalColumn], 0.0);
}

TEST(RunDriversFlatTest, AForceFreeSheetStopsAtTheFirstStep)
{
    ScratchRunDirectory dir;
    write_params("still", {"maxIterations = 4",
                           "VERBOSE_MODE = false",
                           "boundaryType = Periodic",
                           "lFace = 5.0",
                           "sideX = 40.0",
                           "sideY = 40.0",
                           "setRelaxAreaToDefault = true",
                           "meshpointOutput = false",
                           "checkpointOutputInterval = 0",
                           "thermalFluctuationEnabled = false"});
    run_flat("still");
    // No step size can lower a zero energy: only the starting row.
    const auto energy = read_numeric_csv("EnergyForce.csv", 1);
    ASSERT_EQ(energy.size(), 1u);
    EXPECT_EQ(energy.front()[kEnergyTotalColumn], 0.0);
    EXPECT_EQ(energy.front()[kMeanForceColumn], 0.0);
    // The final structure is the starting one.
    const auto begin = read_numeric_csv("vertex_begin.csv", 0);
    const auto final = read_numeric_csv("vertexfinal.csv", 0);
    ASSERT_EQ(begin.size(), final.size());
    for (std::size_t i = 0; i < begin.size(); i++)
    {
        for (int k = 0; k < 3; k++)
        {
            EXPECT_NEAR(final[i][k], begin[i][k], 1e-4) << "vertex " << i;
        }
    }
}

TEST(RunDriversFlatTest, ANearlyRelaxedSheetKeepsLoweringItsEnergy)
{
    ScratchRunDirectory dir;
    // Relaxed area just under the built one: the energy changes by less than
    // 1e-3 per step, which is when the driver refreshes the reference shape.
    // The built sheet is 1039.23 nm^2. At 1039.0 the usable line search
    // lowers the energy by 1.3e-3 a step; 1039.2 keeps every step under the
    // bound and still covers the refresh.
    write_params("near", {"maxIterations = 5",
                          "VERBOSE_MODE = false",
                          "boundaryType = Fixed",
                          "lFace = 5.0",
                          "sideX = 40.0",
                          "sideY = 40.0",
                          "setRelaxAreaToDefault = false",
                          "relaxArea = 1039.2",
                          "meshpointOutput = false",
                          "checkpointOutputInterval = 0"});
    run_flat("near");
    const auto energy = read_numeric_csv("EnergyForce.csv", 1);
    ASSERT_EQ(energy.size(), 5u);
    expect_all_finite(energy, "EnergyForce.csv");
    for (std::size_t r = 1; r < energy.size(); r++)
    {
        EXPECT_LT(std::abs(energy[r][kEnergyTotalColumn] - energy[r - 1][kEnergyTotalColumn]), 1e-3);
    }
    EXPECT_LT(energy.back()[kEnergyTotalColumn], energy.front()[kEnergyTotalColumn]);
}

TEST(RunDriversFlatTest, ScaffoldingRunLiftsTheSheetAndRestartsWithItsBonds)
{
    const std::string scaffold = repo_path("data/csv/GagR50_392.csv");
    ScratchRunDirectory dir;
    const std::vector<std::string> base = {"VERBOSE_MODE = false",
                                           "boundaryType = Fixed",
                                           "lFace = 5.0",
                                           "sideX = 60.0",
                                           "sideY = 60.0",
                                           "setRelaxAreaToDefault = true",
                                           "meshpointOutput = false",
                                           "isEnergyHarmonicBondIncluded = true",
                                           "scaffoldingFileName = " + scaffold,
                                           "scaffoldingZeroPlaneZ = 45.0"};
    write_params("cap", with(base, {"maxIterations = 3",
                                    "checkpointOutputInterval = 2",
                                    "checkpointOutputFile = cap.chk"}));
    run_flat("cap");

    // The sheet starts wrapped around the cap, not flat.
    const auto begin = read_numeric_csv("vertex_begin.csv", 0);
    double maxZ = 0.0;
    for (const auto &row : begin)
    {
        maxZ = std::max(maxZ, row[2]);
    }
    EXPECT_GT(maxZ, 1.0);
    const auto first = read_numeric_csv("EnergyForce.csv", 1);
    ASSERT_EQ(first.size(), 3u);
    expect_all_finite(first, "EnergyForce.csv");
    EXPECT_GT(first.front()[3], 0.0); // E_HarmonicBond
    EXPECT_LT(first.back()[kEnergyTotalColumn], first.front()[kEnergyTotalColumn]);

    // The restart re-establishes the bond correspondence and carries on.
    write_params("cap_restart", with(base, {"maxIterations = 5",
                                            "restartInputFile = cap.chk",
                                            "checkpointOutputInterval = 0"}));
    run_flat("cap_restart");
    const auto second = read_numeric_csv("EnergyForce.csv", 1);
    ASSERT_EQ(second.size(), 5u);
    expect_all_finite(second, "EnergyForce.csv after restart");
    EXPECT_GT(second.back()[3], 0.0);
    EXPECT_NEAR(second.front()[kEnergyTotalColumn], first.front()[kEnergyTotalColumn],
                1e-5 * std::abs(first.front()[kEnergyTotalColumn]));
}

TEST(RunDriversFlatTest, MixedBoundaryExportsTheTypedMesh)
{
    ScratchRunDirectory dir;
    write_params("mixed", {"maxIterations = 4",
                           "VERBOSE_MODE = false",
                           "boundaryType = Mixed",
                           "boundaryTypeX = Periodic",
                           "boundaryTypeY = Fixed",
                           "lFace = 5.0",
                           "sideX = 50.0",
                           "sideY = 40.0",
                           "setRelaxAreaToDefault = false",
                           "relaxArea = 1000.0",
                           "meshpointOutput = false",
                           "checkpointOutputInterval = 0"});
    run_flat("mixed");

    ASSERT_TRUE(file_exists("mesh_vertices.csv"));
    ASSERT_TRUE(file_exists("mesh_faces.csv"));
    // One header line each, then one row per vertex / face.
    EXPECT_EQ(count_lines("mesh_vertices.csv"), count_lines("vertex_begin.csv") + 1);
    EXPECT_EQ(count_lines("mesh_faces.csv"), count_lines("face.csv") + 1);
    // checkpointOutputInterval = 0 disables the checkpoint entirely.
    EXPECT_FALSE(file_exists("slimed_restart.chk"));

    const auto energy = read_numeric_csv("EnergyForce.csv", 1);
    ASSERT_EQ(energy.size(), 4u);
    expect_all_finite(energy, "EnergyForce.csv");
    EXPECT_LT(energy.back()[kEnergyTotalColumn], energy.front()[kEnergyTotalColumn]);
}

// ---------------------------------------------------------------------------
// run_dynamics_flat
// ---------------------------------------------------------------------------

TEST(RunDriversDynamicsTest, BrownianRunWritesTrajectoryAndLocalFrames)
{
    ScratchRunDirectory dir;
    write_params("dyn", {"maxIterations = 4",
                         "VERBOSE_MODE = false",
                         "boundaryType = Periodic",
                         "lFace = 5.0",
                         "sideX = 40.0",
                         "sideY = 40.0",
                         "setRelaxAreaToDefault = true",
                         "meshpointOutput = true",
                         "meshpointOutputInterval = 1",
                         "localOutput = true",
                         "timeStep = 0.001",
                         "diffConst = 1.0",
                         "randomSeed = 7"});
    run_dynamics_flat("dyn");

    for (const char *file : {"dynface.csv", "dynvertex_begin.csv", "dynvertex_type_begin.csv",
                             "dynvertex_final.csv", "dynvertex_type_final.csv",
                             "meshpointdyn.csv", "surfacepointdyn.csv", "EnergyForce.csv"})
    {
        EXPECT_TRUE(file_exists(file)) << file;
    }
    // Not a Mixed boundary, no flips, no bead: none of their outputs.
    for (const char *file : {"dynmesh_vertices.csv", "dynface_0.csv", "dynbead.csv",
                             "dynbox.csv", "dynlipid.csv", "dynEdgeFlips.csv"})
    {
        EXPECT_FALSE(file_exists(file)) << file;
    }

    const int nVertices = count_lines("dynvertex_begin.csv");
    const int nFaces = count_lines("dynface.csv");
    ASSERT_GT(nVertices, 0);
    EXPECT_EQ(count_lines("dynvertex_final.csv"), nVertices);

    // The local files are flushed per frame: the initial frame and one per
    // iteration, four numbers per face and twelve per vertex.
    const auto faceEnergy = read_numeric_csv("faceenergydyn.csv", 0);
    const auto vertexForce = read_numeric_csv("vertexforcedyn.csv", 0);
    ASSERT_EQ(faceEnergy.size(), 5u);
    ASSERT_EQ(vertexForce.size(), 5u);
    for (const auto &frame : faceEnergy)
    {
        EXPECT_EQ(frame.size(), static_cast<std::size_t>(4 * nFaces));
    }
    for (const auto &frame : vertexForce)
    {
        EXPECT_EQ(frame.size(), static_cast<std::size_t>(12 * nVertices));
    }
    expect_all_finite(faceEnergy, "faceenergydyn.csv");
    expect_all_finite(vertexForce, "vertexforcedyn.csv");

    // One record before the first step and one per step.
    const auto energy = read_numeric_csv("EnergyForce.csv", 1);
    ASSERT_EQ(energy.size(), 5u);
    expect_all_finite(energy, "EnergyForce.csv");
    // The noise moved the sheet off z = 10.
    const auto final = read_numeric_csv("dynvertex_final.csv", 0);
    bool moved = false;
    for (const auto &row : final)
    {
        ASSERT_GE(row.size(), 3u);
        moved = moved || row[2] != 10.0;
    }
    EXPECT_TRUE(moved);
}

TEST(RunDriversDynamicsTest, FluidRunWithBeadLipidsAndExternalForceWritesTheirLogs)
{
    ScratchRunDirectory dir;
    write_params("fluid", {"maxIterations = 4",
                           "VERBOSE_MODE = false",
                           "boundaryType = Mixed",
                           "boundaryTypeX = Fixed",
                           "boundaryTypeY = Fixed",
                           "surfaceSolver = iterative",
                           "lFace = 5.0",
                           "sideX = 40.0",
                           "sideY = 40.0",
                           "setRelaxAreaToDefault = true",
                           "meshpointOutput = true",
                           "meshpointOutputInterval = 2",
                           "localOutput = true",
                           "timeStep = 0.001",
                           "diffConst = 1.0",
                           "randomSeed = 7",
                           "inPlaneDynamicsEnabled = true",
                           "edgeSpringEnabled = true",
                           "edgeTetherMinRatio = 0.7",
                           "triangleShapeEnabled = true",
                           "creaseWallEnabled = true",
                           "edgeFlipEnabled = true",
                           "edgeFlipMinValence = 4",
                           "edgeFlipMaxValence = 8",
                           "beadEnabled = true",
                           "beadRadius = 5.0",
                           "beadCenterZ = -10.0",
                           "beadVelocityZ = 1.0",
                           "lipidFlowMode = global",
                           "externalForceEnabled = true",
                           "externalForceZ = 1.0",
                           "externalForceRadius = 8.0"});
    run_dynamics_flat("fluid");

    for (const char *file : {"fluidmesh_vertices.csv", "fluidmesh_faces.csv",
                             "fluidexternal_force_vertices.csv", "fluidface_0.csv",
                             "fluidvertex_final.csv", "EnergyForce.csv"})
    {
        EXPECT_TRUE(file_exists(file)) << file;
    }
    // The initial face frame is the starting triangulation.
    EXPECT_EQ(count_lines("fluidface_0.csv"), count_lines("fluidface.csv"));

    // Frames at iterations 0, 2 and 4: a header and three rows.
    EXPECT_EQ(count_lines("fluidbead.csv"), 4);
    EXPECT_EQ(count_lines("fluidlipid.csv"), 4);
    // Local and lipid face frames: the initial one and two more.
    EXPECT_EQ(count_lines("lipidfacefluid.csv"), 3);
    EXPECT_EQ(count_lines("faceenergyfluid.csv"), 3);
    EXPECT_EQ(count_lines("vertexforcefluid.csv"), 3);
    // The pushed region is a non-empty list of vertices after its header.
    EXPECT_GE(count_lines("fluidexternal_force_vertices.csv"), 2);

    const auto energy = read_numeric_csv("EnergyForce.csv", 1);
    ASSERT_EQ(energy.size(), 5u);
    expect_all_finite(energy, "EnergyForce.csv");
}

TEST(RunDriversDynamicsTest, PureMonteCarloRunLogsFlipsAndVertexMoves)
{
    ScratchRunDirectory dir;
    write_params("mc", {"maxIterations = 4",
                        "VERBOSE_MODE = false",
                        "boundaryType = Mixed",
                        "boundaryTypeX = Periodic",
                        "boundaryTypeY = Periodic",
                        "surfaceSolver = iterative",
                        "lFace = 5.0",
                        "sideX = 50.0",
                        "sideY = 50.0",
                        "setRelaxAreaToDefault = true",
                        "meshpointOutput = true",
                        "meshpointOutputInterval = 1",
                        "timeStep = 1.0",
                        "randomSeed = 7",
                        "brownianStepEnabled = false",
                        "inPlaneDynamicsEnabled = true",
                        "vertexMoveEnabled = true",
                        "vertexMoveStepRatio = 0.1",
                        "meshWallType = hard",
                        "edgeSpringEnabled = true",
                        "edgeTetherMinRatio = 0.7",
                        "triangleShapeEnabled = true",
                        "creaseWallEnabled = true",
                        "edgeFlipEnabled = true",
                        "edgeFlipAttemptRate = 0.1",
                        "edgeFlipMinValence = 4",
                        "edgeFlipMaxValence = 8"});
    run_dynamics_flat("mc");

    // One vertex-move row per trajectory frame after the header.
    EXPECT_EQ(count_lines("mcVertexMoves.csv"), 5);
    // The flip log has a header and one line per attempt.
    ASSERT_TRUE(file_exists("mcEdgeFlips.csv"));
    EXPECT_GE(count_lines("mcEdgeFlips.csv"), 1);
    EXPECT_TRUE(file_exists("mcface_0.csv"));

    // With the Brownian step off the vertex moves are all that moved the net.
    const auto begin = read_numeric_csv("mcvertex_begin.csv", 0);
    const auto final = read_numeric_csv("mcvertex_final.csv", 0);
    ASSERT_EQ(begin.size(), final.size());
    double maxDisplacement = 0.0;
    for (std::size_t i = 0; i < begin.size(); i++)
    {
        for (int k = 0; k < 3; k++)
        {
            maxDisplacement = std::max(maxDisplacement, std::abs(final[i][k] - begin[i][k]));
        }
    }
    EXPECT_GT(maxDisplacement, 0.0);
    // A trial cube of half-width 0.1 lFace, at most four sweeps.
    EXPECT_LT(maxDisplacement, 4 * 0.1 * 5.0 * std::sqrt(3.0) + 1e-9);

    const auto energy = read_numeric_csv("EnergyForce.csv", 1);
    ASSERT_EQ(energy.size(), 5u);
    expect_all_finite(energy, "EnergyForce.csv");
}

TEST(RunDriversDynamicsTest, HardWallsRefuseAStartingMeshOutsideThem)
{
    ScratchRunDirectory dir;
    // Every edge is shorter than 1.2 lFace, so every tether wall is broken.
    write_params("wall", {"maxIterations = 4",
                          "VERBOSE_MODE = false",
                          "boundaryType = Mixed",
                          "boundaryTypeX = Periodic",
                          "boundaryTypeY = Periodic",
                          "surfaceSolver = iterative",
                          "lFace = 5.0",
                          "sideX = 50.0",
                          "sideY = 50.0",
                          "setRelaxAreaToDefault = true",
                          "meshpointOutput = false",
                          "timeStep = 1.0",
                          "randomSeed = 7",
                          "brownianStepEnabled = false",
                          "inPlaneDynamicsEnabled = true",
                          "vertexMoveEnabled = true",
                          "meshWallType = hard",
                          "edgeSpringEnabled = true",
                          "edgeTetherMinRatio = 1.2",
                          "triangleShapeEnabled = true",
                          "creaseWallEnabled = true"});
    EXPECT_THROW(run_dynamics_flat("wall"), std::runtime_error);
    EXPECT_FALSE(file_exists("wallvertex_final.csv"));
}

TEST(RunDriversDynamicsTest, ElasticBoxRunLogsTheBox)
{
    ScratchRunDirectory dir;
    write_params("box", {"maxIterations = 4",
                         "VERBOSE_MODE = false",
                         "boundaryType = Periodic",
                         "lFace = 5.0",
                         "sideX = 40.0",
                         "sideY = 40.0",
                         "setRelaxAreaToDefault = true",
                         "meshpointOutput = true",
                         "meshpointOutputInterval = 2",
                         "timeStep = 0.001",
                         "randomSeed = 7",
                         "elasticBoxEnabled = true",
                         "elasticBoxMobility = 0.001"});
    run_dynamics_flat("box");
    // Frames at iterations 0, 2 and 4.
    EXPECT_EQ(count_lines("boxbox.csv"), 4);
    const auto box = read_numeric_csv("boxbox.csv", 1);
    expect_all_finite(box, "boxbox.csv");
    const auto energy = read_numeric_csv("EnergyForce.csv", 1);
    ASSERT_EQ(energy.size(), 5u);
    expect_all_finite(energy, "EnergyForce.csv");
}

TEST(RunDriversDynamicsTest, LocalLipidFlowWritesPerFaceFrames)
{
    ScratchRunDirectory dir;
    write_params("lip", {"maxIterations = 3",
                         "VERBOSE_MODE = false",
                         "boundaryType = Mixed",
                         "boundaryTypeX = Periodic",
                         "boundaryTypeY = Periodic",
                         "surfaceSolver = iterative",
                         "lFace = 5.0",
                         "sideX = 50.0",
                         "sideY = 50.0",
                         "setRelaxAreaToDefault = true",
                         "meshpointOutput = true",
                         "localOutput = true",
                         "timeStep = 0.001",
                         "randomSeed = 7",
                         "lipidFlowMode = local"});
    run_dynamics_flat("lip");
    // Header plus the initial frame and one per iteration.
    EXPECT_EQ(count_lines("liplipid.csv"), 5);
    EXPECT_EQ(count_lines("lipidfacelip.csv"), 4);
    const auto lipid = read_numeric_csv("liplipid.csv", 1);
    expect_all_finite(lipid, "liplipid.csv");
}

TEST(RunDriversDynamicsTest, ScaffoldingRunStartsWrappedAroundTheCap)
{
    const std::string scaffold = repo_path("data/csv/GagR50_392.csv");
    ScratchRunDirectory dir;
    write_params("gag", {"maxIterations = 3",
                         "VERBOSE_MODE = false",
                         "boundaryType = Fixed",
                         "lFace = 5.0",
                         "sideX = 60.0",
                         "sideY = 60.0",
                         "setRelaxAreaToDefault = true",
                         "meshpointOutput = false",
                         "timeStep = 0.001",
                         "randomSeed = 7",
                         "isEnergyHarmonicBondIncluded = true",
                         "scaffoldingFileName = " + scaffold,
                         "scaffoldingZeroPlaneZ = 45.0"});
    run_dynamics_flat("gag");
    // The generated sheet is put at z = 10 and the cap lifts part of it.
    const auto begin = read_numeric_csv("gagvertex_begin.csv", 0);
    double maxZ = -1e300;
    double minZ = 1e300;
    for (const auto &row : begin)
    {
        maxZ = std::max(maxZ, row[2]);
        minZ = std::min(minZ, row[2]);
    }
    EXPECT_GE(minZ, 10.0 - 1e-9);
    EXPECT_GT(maxZ, 11.0);
    const auto energy = read_numeric_csv("EnergyForce.csv", 1);
    ASSERT_EQ(energy.size(), 4u);
    expect_all_finite(energy, "EnergyForce.csv");
    EXPECT_GT(energy.front()[3], 0.0); // E_HarmonicBond
}

TEST(RunDriversDynamicsTest, ADivergedRunStopsAtTheFirstNonFiniteEnergy)
{
    ScratchRunDirectory dir;
    // An absurd time step blows the sheet up within a few steps.
    write_params("blowup", {"maxIterations = 20",
                            "VERBOSE_MODE = false",
                            "boundaryType = Periodic",
                            "lFace = 5.0",
                            "sideX = 40.0",
                            "sideY = 40.0",
                            "setRelaxAreaToDefault = true",
                            "meshpointOutput = false",
                            "timeStep = 1e60",
                            "randomSeed = 7"});
    run_dynamics_flat("blowup");
    const auto energy = read_numeric_csv("EnergyForce.csv", 1);
    // Stopped well before the 21 rows a full run would write...
    ASSERT_GE(energy.size(), 2u);
    EXPECT_LT(energy.size(), 21u);
    // ...and everything up to the divergence was still written.
    EXPECT_TRUE(file_exists("blowupvertex_final.csv"));
    EXPECT_TRUE(std::isfinite(energy.front()[kEnergyTotalColumn]));
}

namespace
{
/// The flat Mixed periodic tile of docs/resolution_limits.md, shrunk, with the
/// fluid mesh walls in the drive and the walls substepped.
std::vector<std::string> substepped_periodic_sheet(const std::string &noise)
{
    return {"maxIterations = 6",
            "VERBOSE_MODE = false",
            "boundaryType = Mixed",
            "boundaryTypeX = Periodic",
            "boundaryTypeY = Periodic",
            "meanHeightPinned = true",
            "surfaceSolver = iterative",
            "lFace = 5.0",
            "sideX = 80.0",
            "sideY = 77.94228634",
            "setRelaxAreaToDefault = true",
            "meshpointOutput = true",
            "meshpointOutputInterval = 3",
            "timeStep = 0.004",
            "diffConst = 1.0",
            "randomSeed = 7",
            "inPlaneDynamicsEnabled = true",
            "edgeSpringEnabled = true",
            "triangleShapeEnabled = true",
            "creaseWallEnabled = true",
            "meshForceSubsteps = 3",
            "meshForceSubstepNoise = " + noise};
}
} // namespace

TEST(RunDriversDynamicsTest, SubsteppedWallsRunOnAPeriodicMixedSheet)
{
    ScratchRunDirectory dir;
    for (const std::string noise : {"split", "once"})
    {
        const std::string name = "sub" + noise;
        write_params(name, substepped_periodic_sheet(noise));
        run_dynamics_flat(name);
        EXPECT_TRUE(file_exists(name + "vertex_final.csv")) << noise;
        // One record before the first step and one per step, none of them blown up.
        const auto energy = read_numeric_csv("EnergyForce.csv", 1);
        ASSERT_EQ(energy.size(), 7u) << noise;
        expect_all_finite(energy, "EnergyForce.csv");
        const auto final = read_numeric_csv(name + "vertex_final.csv", 0);
        expect_all_finite(final, name + "vertex_final.csv");
    }
}

TEST(RunDriversDynamicsTest, AnUnknownSubstepNoiseModeIsRefused)
{
    ScratchRunDirectory dir;
    write_params("subbad", substepped_periodic_sheet("sometimes"));
    try
    {
        run_dynamics_flat("subbad");
        FAIL() << "an unknown meshForceSubstepNoise was accepted";
    }
    catch (const std::exception &error)
    {
        EXPECT_NE(std::string(error.what()).find("meshForceSubstepNoise"), std::string::npos) << error.what();
    }
}
