#include "test_scaffolding_coverage.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unistd.h>

/**
 * Coverage of src/energy_force/Scaffolding_points.cpp.
 *
 * The file-local parsers and geometry helpers are reached through the public
 * Mesh methods: tiny reference-state (.dat / JSON-like) and NERDSS reaction
 * (.inp) files are written into a fresh temporary directory and the Param is
 * pointed at them. The scaffold is a three-subunit ring with generic (non
 * planar, non colinear) interfaces so every angle term is away from its
 * singular points; the membrane is a small flat periodic sheet.
 */

namespace fs = std::filesystem;

namespace
{

// ---------------------------------------------------------------------------
// Small utilities
// ---------------------------------------------------------------------------

Matrix vec3(double x, double y, double z)
{
    Matrix v(3, 1, true);
    v.set(0, 0, x);
    v.set(1, 0, y);
    v.set(2, 0, z);
    return v;
}

Matrix vec3(const std::array<double, 3> &a) { return vec3(a[0], a[1], a[2]); }

double norm3(const Matrix &v)
{
    return std::sqrt(v(0, 0) * v(0, 0) + v(1, 0) * v(1, 0) + v(2, 0) * v(2, 0));
}

Matrix cross3(const Matrix &a, const Matrix &b)
{
    return vec3(a(1, 0) * b(2, 0) - a(2, 0) * b(1, 0),
                a(2, 0) * b(0, 0) - a(0, 0) * b(2, 0),
                a(0, 0) * b(1, 0) - a(1, 0) * b(0, 0));
}

/// Rodrigues rotation, written independently of the code under test.
Matrix rotation(const std::array<double, 3> &axisIn, double angle)
{
    const double n = std::sqrt(axisIn[0] * axisIn[0] + axisIn[1] * axisIn[1] + axisIn[2] * axisIn[2]);
    const double x = axisIn[0] / n, y = axisIn[1] / n, z = axisIn[2] / n;
    const double c = std::cos(angle), s = std::sin(angle), C = 1.0 - c;
    Matrix r(3, 3, true);
    r.set(0, 0, c + x * x * C);
    r.set(0, 1, x * y * C - z * s);
    r.set(0, 2, x * z * C + y * s);
    r.set(1, 0, y * x * C + z * s);
    r.set(1, 1, c + y * y * C);
    r.set(1, 2, y * z * C - x * s);
    r.set(2, 0, z * x * C - y * s);
    r.set(2, 1, z * y * C + x * s);
    r.set(2, 2, c + z * z * C);
    return r;
}

/// max |R^T R - I| for a 3x3 matrix.
double orthogonality_defect(const Matrix &r)
{
    double worst = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            double dot = 0.0;
            for (int k = 0; k < 3; ++k)
            {
                dot += r(k, i) * r(k, j);
            }
            worst = std::max(worst, std::abs(dot - (i == j ? 1.0 : 0.0)));
        }
    }
    return worst;
}

/// A fresh, uniquely named directory under the system temp dir; removed on destruction.
class TempDir
{
public:
    TempDir()
    {
        static std::atomic<int> counter{0};
        const ::testing::TestInfo *info = ::testing::UnitTest::GetInstance()->current_test_info();
        std::string name = std::string("slimed_scafcov_") + std::to_string(::getpid()) + "_" +
                           std::to_string(counter++) + "_" + (info ? info->name() : "anon");
        // parameterised test names contain '/': keep the directory a single path component
        std::replace_if(name.begin(), name.end(),
                        [](char ch) { return !std::isalnum(static_cast<unsigned char>(ch)) && ch != '_'; }, '_');
        path_ = fs::temp_directory_path() / name;
        fs::remove_all(path_);
        fs::create_directories(path_);
    }
    ~TempDir()
    {
        std::error_code ec;
        fs::remove_all(path_, ec);
    }
    const fs::path &path() const { return path_; }
    std::string file(const std::string &name) const { return (path_ / name).string(); }

private:
    fs::path path_;
};

/// chdir into a directory for the lifetime of the guard; always restores the cwd.
class CwdGuard
{
public:
    explicit CwdGuard(const fs::path &dir) : original_(fs::current_path()) { fs::current_path(dir); }
    ~CwdGuard()
    {
        std::error_code ec;
        fs::current_path(original_, ec);
    }

private:
    fs::path original_;
};

void write_text(const std::string &path, const std::string &text)
{
    std::ofstream out(path);
    out << text;
}

// ---------------------------------------------------------------------------
// Reference-state (.dat) writer in the layout Mesh::write_gag_scaffolding_state_dat emits
// ---------------------------------------------------------------------------

struct RefIface
{
    std::string name;
    std::array<double, 3> coord;
    std::string refKey; // "", "sphere_normal" or "tangent_direction"
    std::array<double, 3> ref{{0.0, 0.0, 0.0}};
};

struct RefMol
{
    int id;
    std::array<double, 3> com;
    std::vector<RefIface> ifaces;
};

struct RefBond
{
    int mol1;
    std::string if1;
    int mol2;
    std::string if2;
    bool malformed = false;
};

void write_triple(std::ostream &o, const std::array<double, 3> &v, const std::string &indent)
{
    o << indent << v[0] << ",\n" << indent << v[1] << ",\n" << indent << v[2] << "\n";
}

void write_reference_state(const std::string &path, const std::vector<RefMol> &mols,
                           const std::vector<RefBond> &bonds, bool withMetaBlock = false)
{
    std::ofstream o(path);
    o << std::setprecision(17);
    o << "{\n  \"molecules\": {\n";
    if (withMetaBlock)
    {
        // a non-numeric key inside "molecules" must be ignored
        o << "    \"meta\": {\n      \"note\": 1\n    },\n";
    }
    for (size_t m = 0; m < mols.size(); ++m)
    {
        const RefMol &mol = mols[m];
        o << "    \"" << mol.id << "\": {\n";
        o << "      \"com\": [\n";
        write_triple(o, mol.com, "        ");
        o << "      ],\n";
        o << "      \"interfaces\": {\n";
        for (size_t i = 0; i < mol.ifaces.size(); ++i)
        {
            const RefIface &iface = mol.ifaces[i];
            o << "        \"" << iface.name << "\": {\n";
            o << "          \"coord\": [\n";
            write_triple(o, iface.coord, "            ");
            o << (iface.refKey.empty() ? "          ]\n" : "          ],\n");
            if (!iface.refKey.empty())
            {
                o << "          \"" << iface.refKey << "\": [\n";
                write_triple(o, iface.ref, "            ");
                o << "          ]\n";
            }
            o << "        }" << (i + 1 < mol.ifaces.size() ? ",\n" : "\n");
        }
        o << "      }\n";
        o << "    }" << (m + 1 < mols.size() ? ",\n" : "\n");
    }
    o << "  },\n  \"bound_states\": [\n";
    for (size_t b = 0; b < bonds.size(); ++b)
    {
        const RefBond &bond = bonds[b];
        o << "    {\n";
        o << "      \"mol1\": {\n";
        if (bond.malformed)
        {
            o << "        \"molecule_index\" " << bond.mol1 << ",\n"; // no ':' -> skipped
        }
        else
        {
            o << "        \"molecule_index\": " << bond.mol1 << ",\n";
        }
        o << "        \"interface_name\": \"" << bond.if1 << "\",\n";
        o << "        \"interface_index\": 0\n";
        o << "      },\n";
        o << "      \"mol2\": {\n";
        o << "        \"molecule_index\": " << bond.mol2 << ",\n";
        o << "        \"interface_name\": \"" << bond.if2 << "\",\n";
        o << "        \"interface_index\": 0\n";
        o << "      }\n";
        o << "    }" << (b + 1 < bonds.size() ? ",\n" : "\n");
    }
    o << "  ]\n}\n";
}

std::array<double, 3> along(const std::array<double, 3> &from, const std::array<double, 3> &to,
                            double frac, const std::array<double, 3> &offset)
{
    return {{from[0] + frac * (to[0] - from[0]) + offset[0],
             from[1] + frac * (to[1] - from[1]) + offset[1],
             from[2] + frac * (to[2] - from[2]) + offset[2]}};
}

/// Three-subunit ring with generic interfaces. Bonds 0.b-1.a, 1.b-2.a, 2.c-0.a,
/// a duplicate of the first (ignored) and a malformed block (ignored).
struct RingLattice
{
    std::vector<RefMol> mols;
    std::vector<RefBond> bonds;
    std::array<double, 3> translation{{0.01, -0.02, 0.03}};

    RingLattice()
    {
        const std::array<double, 3> c0{{0.0, 0.0, 0.0}};
        const std::array<double, 3> c1{{10.0, 0.0, 1.0}};
        const std::array<double, 3> c2{{4.0, 9.0, -1.0}};
        RefMol m0{0, c0, {}};
        // first reference vector wins: sphere_normal on "a", the tangent on "b" is ignored
        m0.ifaces.push_back({"a", along(c0, c2, 0.4, {{0.3, 0.0, 0.5}}), "sphere_normal", {{0.0, 0.0, 1.0}}});
        m0.ifaces.push_back({"b", along(c0, c1, 0.4, {{0.0, 0.4, -0.3}}), "tangent_direction", {{1.0, 0.0, 0.0}}});
        RefMol m1{1, c1, {}};
        // first reference vector wins: tangent on "a", the sphere_normal on "b" is ignored
        m1.ifaces.push_back({"a", along(c1, c0, 0.4, {{0.0, -0.5, 0.2}}), "tangent_direction", {{0.0, 1.0, 0.0}}});
        m1.ifaces.push_back({"b", along(c1, c2, 0.4, {{0.2, 0.1, 0.4}}), "sphere_normal", {{0.0, 0.0, -1.0}}});
        RefMol m2{2, c2, {}};
        // no reference vector at all
        m2.ifaces.push_back({"a", along(c2, c1, 0.4, {{-0.3, 0.2, 0.1}}), "", {{0.0, 0.0, 0.0}}});
        m2.ifaces.push_back({"c", along(c2, c0, 0.4, {{0.1, -0.4, 0.3}}), "", {{0.0, 0.0, 0.0}}});
        mols = {m0, m1, m2};
        bonds = {{0, "b", 1, "a"}, {1, "b", 2, "a"}, {2, "c", 0, "a"}, {1, "a", 0, "b"}, {0, "a", 2, "c", true}};
    }

    /// Scaffold points: the reference COMs translated rigidly, deliberately listed
    /// out of molecule order (point 0 = molecule 2, point 1 = molecule 0, point 2 = molecule 1).
    std::vector<Matrix> points() const
    {
        const int order[3] = {2, 0, 1};
        std::vector<Matrix> out;
        for (int id : order)
        {
            const auto &c = mols[id].com;
            out.push_back(vec3(c[0] + translation[0], c[1] + translation[1], c[2] + translation[2]));
        }
        return out;
    }

    std::array<double, 3> iface(int mol, const std::string &name) const
    {
        for (const RefIface &i : mols[mol].ifaces)
        {
            if (i.name == name)
            {
                return i.coord;
            }
        }
        throw std::runtime_error("no such interface");
    }
};

double distance(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

/// NERDSS-style reaction file: one valid b|a reaction with bangs, a sigma-less
/// c|c reaction, and every kind of malformed line the parser must skip.
const char *kReactionFile =
    "# Gag reactions\n"
    "start molecules\n"
    "  A(a,b,c)\n"
    "end molecules\n"
    "\n"
    "start reactions\n"
    "    sigma = 5.0\n"                              // before any reaction: ignored
    "    B(x) + B(y) <-> B(x!1).B(y!1)\n"            // no A( : not a Gag reaction
    "    A(b!1) + A(a!2) <-> A(b!1).A(a!1)\n"        // the reaction (bangs stripped)
    "        onRate3Dka = 10\n"
    "        sigma = 2.5\n"
    "        assocAngles = [1.1, 2.2, 0.3, 0.5, -0.7]\n"
    "    A(a) + A(c <-> broken\n"                     // right-hand ')' missing
    "    A(x) <-> nothing\n"                          // '+' missing
    "    A(c) + A(c) <-> A(c!1).A(c!1)\n"            // stores b|a; itself has no sigma
    "        assocAngles = [1, 2, 3\n"                // no ']' : skipped
    "end reactions\n"
    "    A(z) + A(z) <-> A(z!1).A(z!1)\n"            // after the block: never read
    "        sigma = 9.0\n";

// ---------------------------------------------------------------------------
// Mesh fixtures
// ---------------------------------------------------------------------------

void configure_sheet(Param &p, double side, double lFace)
{
    p.VERBOSE_MODE = false;
    p.boundaryCondition = BoundaryType::Periodic;
    p.sideX = side;
    p.sideY = side;
    p.lFace = lFace;
    p.isEnergyHarmonicBondIncluded = false;
    p.isGagScaffoldingEnergyIncluded = false;
    p.isIdealizedProteinLatticeEnergyIncluded = false;
    p.lbond = 9.0;
    p.springConst = 4.5;
    p.insertCurv = 0.0;
    p.spontCurv = 0.0;
}

/// Exposes the protected squared-distance helper.
class ScafCovMesh : public Mesh
{
public:
    explicit ScafCovMesh(Param &p) : Mesh(p) {}
    using Mesh::get_squared_distance_sp_and_v;
};

int nearest_interior_vertex(const Mesh &mesh, double x, double y)
{
    int best = -1;
    double bestD = std::numeric_limits<double>::max();
    for (size_t i = 0; i < mesh.vertices.size(); ++i)
    {
        const double dx = mesh.vertices[i].coord(0, 0) - x;
        const double dy = mesh.vertices[i].coord(1, 0) - y;
        const double d = dx * dx + dy * dy;
        if (d < bestD)
        {
            bestD = d;
            best = static_cast<int>(i);
        }
    }
    return best;
}

/// Lattice-mode Param on its own (no mesh needed for the internal energy).
void configure_lattice(Param &p, const std::string &refFile, bool idealized)
{
    p.VERBOSE_MODE = false;
    p.isEnergyHarmonicBondIncluded = false;
    p.isIdealizedProteinLatticeEnergyIncluded = idealized;
    p.isGagScaffoldingEnergyIncluded = !idealized;
    if (idealized)
    {
        p.idealizedProteinLatticeFileName = refFile;
    }
    else
    {
        p.gagReferenceStateFileName = refFile;
    }
    p.gagKsigma = 1.0;
    p.gagKtheta = 0.7;
    p.gagKphi = 0.3;
    p.gagKomega = 0.5;
    p.gagFiniteDifferenceStep = 1.0e-4;
}

/// A deterministic, non-rigid perturbation of the scaffold.
void perturb(std::vector<Matrix> &points)
{
    const double d[3][3] = {{0.31, -0.17, 0.22}, {-0.12, 0.27, -0.19}, {0.08, 0.14, 0.33}};
    for (size_t i = 0; i < points.size(); ++i)
    {
        for (int a = 0; a < 3; ++a)
        {
            points[i].set(a, 0, points[i](a, 0) + d[i % 3][a]);
        }
    }
}

} // namespace

// ===========================================================================
// Membrane shape profile helpers
// ===========================================================================

TEST(ScafCovShapeProfile, FindClosestRcapPicksNearestTabulatedRadius)
{
    Param p;
    Mesh mesh(p);
    EXPECT_EQ(mesh.findClosestRcap(0.0), 10);
    EXPECT_EQ(mesh.findClosestRcap(-50.0), 10);
    EXPECT_EQ(mesh.findClosestRcap(14.9), 10);
    EXPECT_EQ(mesh.findClosestRcap(22.0), 20);
    EXPECT_EQ(mesh.findClosestRcap(23.0), 25);
    EXPECT_EQ(mesh.findClosestRcap(22.5), 20); // a tie keeps the first (smaller) entry
    EXPECT_EQ(mesh.findClosestRcap(36.0), 35);
    EXPECT_EQ(mesh.findClosestRcap(1000.0), 40);
}

TEST(ScafCovShapeProfile, InterpolateHeightClampsAndInterpolatesLinearly)
{
    Param p;
    Mesh mesh(p);
    const std::vector<double> r = {1.0, 2.0, 4.0};
    const std::vector<double> h = {10.0, 6.0, -2.0};
    EXPECT_DOUBLE_EQ(mesh.interpolateHeight(0.0, r, h), 10.0);  // below the table
    EXPECT_DOUBLE_EQ(mesh.interpolateHeight(1.0, r, h), 10.0);  // first node
    EXPECT_DOUBLE_EQ(mesh.interpolateHeight(1.5, r, h), 8.0);
    EXPECT_DOUBLE_EQ(mesh.interpolateHeight(2.0, r, h), 6.0);   // interior node
    EXPECT_DOUBLE_EQ(mesh.interpolateHeight(3.0, r, h), 2.0);
    EXPECT_DOUBLE_EQ(mesh.interpolateHeight(4.0, r, h), -2.0);  // last node
    EXPECT_DOUBLE_EQ(mesh.interpolateHeight(99.0, r, h), -2.0); // above the table
}

TEST(ScafCovShapeProfile, LoadMembraneShapeProfileReadsClosestFileFromCwd)
{
    Param p;
    Mesh mesh(p);
    TempDir tmp;
    // No header; a one-column row and a blank row are dropped.
    write_text(tmp.file("membrane_shape_25.csv"), "25,0.0\n30,-1.5\n7\n\n40,-3.0\n");

    std::vector<double> r = {99.0};
    std::vector<double> h = {99.0, 98.0};
    {
        CwdGuard cwd(tmp.path());
        mesh.loadMembraneShapeProfile(24.0, r, h); // closest tabulated radius is 25
    }
    ASSERT_EQ(r.size(), 3u);
    ASSERT_EQ(h.size(), 3u);
    EXPECT_DOUBLE_EQ(r[0], 25.0);
    EXPECT_DOUBLE_EQ(r[1], 30.0);
    EXPECT_DOUBLE_EQ(r[2], 40.0);
    EXPECT_DOUBLE_EQ(h[0], 0.0);
    EXPECT_DOUBLE_EQ(h[1], -1.5);
    EXPECT_DOUBLE_EQ(h[2], -3.0);
    EXPECT_DOUBLE_EQ(mesh.interpolateHeight(35.0, r, h), -2.25);
}

TEST(ScafCovShapeProfile, LoadMembraneShapeProfileThrowsOnMissingOrHeaderedFile)
{
    Param p;
    Mesh mesh(p);
    const fs::path shipped = fs::absolute("data/shape_profiles/membrane_shape_10.csv");
    TempDir tmp;
    std::vector<double> r, h;
    CwdGuard cwd(tmp.path());
    // membrane_shape_40.csv does not exist in the temp dir
    EXPECT_THROW(mesh.loadMembraneShapeProfile(100.0, r, h), std::invalid_argument);
    // The shipped profiles start with an "r,h" header, which read_data_from_csv
    // cannot convert: document that they cannot be loaded as they stand.
    if (fs::exists(shipped))
    {
        fs::copy_file(shipped, tmp.path() / "membrane_shape_10.csv");
        EXPECT_THROW(mesh.loadMembraneShapeProfile(10.0, r, h), std::invalid_argument);
    }
}

// ===========================================================================
// Membrane initialisation from the scaffold cap
// ===========================================================================

class ScafCovMoveVertices : public ::testing::Test
{
protected:
    static constexpr double kR = 20.0;
    static constexpr double kAlpha = 0.5;
    static constexpr double kCx = 0.7;
    static constexpr double kCy = -0.4;

    void SetUp() override
    {
        configure_sheet(param, 80.0, 8.0);
        param.scaffoldingSphereRaidus = kR;
        param.relaxLengthRatioApproximation = 1.0;
        param.scaffoldingZeroPlaneZ = 0.0;
        param.insertCurv = 0.25;
        param.spontCurv = 0.0;
        mesh.reset(new Mesh(param));
        mesh->setup_flat();
    }

    /// A spherical cap of radius kR, apex at (kCx, kCy, apexZ), half angle kAlpha.
    std::vector<Matrix> cap(double apexZ) const
    {
        std::vector<Matrix> pts;
        pts.push_back(vec3(kCx, kCy, apexZ));
        for (int ring = 1; ring <= 2; ++ring)
        {
            const double a = kAlpha * ring / 2.0;
            for (int k = 0; k < 8; ++k)
            {
                const double phi = 2.0 * M_PI * k / 8.0;
                pts.push_back(vec3(kCx + kR * std::sin(a) * std::cos(phi),
                                   kCy + kR * std::sin(a) * std::sin(phi),
                                   apexZ - kR + kR * std::cos(a)));
            }
        }
        return pts;
    }

    Param param;
    std::unique_ptr<Mesh> mesh;
};

TEST_F(ScafCovMoveVertices, FixedDirectionRaisesEveryVertexByBondLength)
{
    std::vector<double> z0;
    for (const Vertex &v : mesh->vertices)
    {
        z0.push_back(v.coord(2, 0));
    }
    param.lbond = 6.5;
    EXPECT_TRUE(mesh->move_vertices_based_on_scaffolding(true));
    for (size_t i = 0; i < mesh->vertices.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(mesh->vertices[i].coord(2, 0), z0[i] + 6.5);
    }
}

TEST_F(ScafCovMoveVertices, CapVerticesSitOnTheMembraneSphereAndTransitionIsMonotone)
{
    param.scaffoldingPoints = cap(0.0);
    // One far vertex starts high: the initialisation must never lower a vertex.
    const int high = nearest_interior_vertex(*mesh, 36.0, 36.0);
    mesh->vertices[high].coord.set(2, 0, 3.0);

    testing::internal::CaptureStdout();
    ASSERT_TRUE(mesh->move_vertices_based_on_scaffolding(false));
    testing::internal::GetCapturedStdout();

    const double l = param.lbond;
    const double capRadius = kR * std::sin(kAlpha);
    const double r0 = capRadius * (kR + l) / kR;          // membrane cap radius
    const double rR = 2.0 * r0;                            // relax length ratio 1
    const double centerZ = -kR;
    const double capEdgeZ = centerZ + (kR + l) * std::cos(kAlpha);

    EXPECT_NEAR(mesh->centerScaffoldingSphere(0, 0), kCx, 1e-12);
    EXPECT_NEAR(mesh->centerScaffoldingSphere(1, 0), kCy, 1e-12);
    EXPECT_NEAR(mesh->centerScaffoldingSphere(2, 0), centerZ, 1e-12);

    int nCap = 0, nTransition = 0, nOutside = 0;
    std::vector<std::pair<double, double>> transition; // (rho, z)
    for (size_t i = 0; i < mesh->vertices.size(); ++i)
    {
        const Vertex &v = mesh->vertices[i];
        const double dx = v.coord(0, 0) - kCx;
        const double dy = v.coord(1, 0) - kCy;
        const double rho = std::sqrt(dx * dx + dy * dy);
        const double z = v.coord(2, 0);
        ASSERT_TRUE(std::isfinite(z));
        if (rho * rho <= r0 * r0)
        {
            ++nCap;
            // on the sphere of radius R + l about the scaffold sphere centre
            const double dz = z - centerZ;
            EXPECT_NEAR(std::sqrt(rho * rho + dz * dz), kR + l, 1e-9);
            EXPECT_GE(z, capEdgeZ - 1e-9);
            for (int f : v.adjacentFaces)
            {
                EXPECT_DOUBLE_EQ(mesh->faces[f].spontCurvature, param.insertCurv);
            }
        }
        else if (rho * rho <= rR * rR)
        {
            ++nTransition;
            EXPECT_GE(z, 0.0);
            EXPECT_LE(z, capEdgeZ + 1e-9);
            transition.emplace_back(rho, z);
        }
        else
        {
            ++nOutside;
            EXPECT_DOUBLE_EQ(z, static_cast<int>(i) == high ? 3.0 : 0.0);
        }
    }
    EXPECT_GT(nCap, 0);
    EXPECT_GT(nTransition, 2);
    EXPECT_GT(nOutside, 0);

    // The tensionless profile decreases from the cap edge to the flat plane.
    std::sort(transition.begin(), transition.end());
    for (size_t k = 1; k < transition.size(); ++k)
    {
        EXPECT_LE(transition[k].second, transition[k - 1].second + 1e-9);
    }
    EXPECT_GT(transition.front().second, 0.0);
}

TEST_F(ScafCovMoveVertices, BoundaryVertexInTransitionStaysOnTheFlatPlane)
{
    param.scaffoldingPoints = cap(0.0);
    param.scaffoldingZeroPlaneZ = 1.5;
    const double r0 = kR * std::sin(kAlpha) * (kR + param.lbond) / kR;
    // a vertex at about 1.4 r0 from the axis: inside the transition annulus
    const int v = nearest_interior_vertex(*mesh, kCx + 1.4 * r0, kCy);
    const double dx = mesh->vertices[v].coord(0, 0) - kCx;
    const double dy = mesh->vertices[v].coord(1, 0) - kCy;
    const double rho = std::sqrt(dx * dx + dy * dy);
    ASSERT_GT(rho, r0);
    ASSERT_LT(rho, 2.0 * r0);
    mesh->vertices[v].isBoundary = true;

    testing::internal::CaptureStdout();
    ASSERT_TRUE(mesh->move_vertices_based_on_scaffolding(false));
    testing::internal::GetCapturedStdout();

    EXPECT_DOUBLE_EQ(mesh->vertices[v].coord(2, 0), 1.5);
    // and the far field sits on the raised flat plane
    const int far = nearest_interior_vertex(*mesh, 36.0, -36.0);
    EXPECT_DOUBLE_EQ(mesh->vertices[far].coord(2, 0), 1.5);
    for (const Vertex &vertex : mesh->vertices)
    {
        EXPECT_GE(vertex.coord(2, 0), 1.5);
    }
}

TEST_F(ScafCovMoveVertices, DegenerateRadiiAreClampedOrDefaulted)
{
    // Cap wider than the sphere: the cap radius is clamped just below R.
    param.scaffoldingPoints = cap(0.0);
    param.scaffoldingSphereRaidus = 5.0;
    testing::internal::CaptureStdout();
    ASSERT_TRUE(mesh->move_vertices_based_on_scaffolding(false));
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_NE(out.find("Clamping"), std::string::npos);
    for (const Vertex &v : mesh->vertices)
    {
        ASSERT_TRUE(std::isfinite(v.coord(2, 0)));
        EXPECT_GE(v.coord(2, 0), 0.0);
    }

    // Non-positive sphere radius falls back to 50 nm.
    Param p2;
    configure_sheet(p2, 80.0, 8.0);
    p2.scaffoldingSphereRaidus = -1.0;
    p2.scaffoldingPoints = cap(0.0);
    Mesh m2(p2);
    m2.setup_flat();
    testing::internal::CaptureStdout();
    ASSERT_TRUE(m2.move_vertices_based_on_scaffolding(false));
    const std::string out2 = testing::internal::GetCapturedStdout();
    EXPECT_NE(out2.find("Membrane sphere radius approximated : 59"), std::string::npos);
    for (const Vertex &v : m2.vertices)
    {
        ASSERT_TRUE(std::isfinite(v.coord(2, 0)));
        EXPECT_GE(v.coord(2, 0), 0.0);
    }
}

TEST_F(ScafCovMoveVertices, VanishingRelaxProfileUsesSmoothstepBlend)
{
    // A cap so small that the tensionless relax profile is below 1e-12 nm: the
    // transition falls back to a smoothstep in the radial coordinate. Place one
    // vertex half-way across the (very thin) transition annulus.
    param.relaxLengthRatioApproximation = 0.0; // floored to 1e-6
    const double g = 1.0e-4;                    // scaffold cap radius
    const double r0 = g * (kR + param.lbond) / kR;
    const double rho = r0 * (1.0 + 5.0e-7);
    const int v = nearest_interior_vertex(*mesh, 4.0, 4.0);
    const double vx = mesh->vertices[v].coord(0, 0);
    const double vy = mesh->vertices[v].coord(1, 0);
    param.scaffoldingPoints = {vec3(vx + rho, vy, 0.0), vec3(vx + rho + g, vy, -1.0e-3)};

    testing::internal::CaptureStdout();
    ASSERT_TRUE(mesh->move_vertices_based_on_scaffolding(false));
    testing::internal::GetCapturedStdout();

    // capEdge = (apex - R) + sqrt((R+l)^2 - r0^2) ~ l; smoothstep(1/2) = 1/2
    const double capEdge = -kR + std::sqrt((kR + param.lbond) * (kR + param.lbond) - r0 * r0);
    EXPECT_NEAR(mesh->vertices[v].coord(2, 0), 0.5 * capEdge, 1e-3);
    for (size_t i = 0; i < mesh->vertices.size(); ++i)
    {
        if (static_cast<int>(i) != v)
        {
            EXPECT_DOUBLE_EQ(mesh->vertices[i].coord(2, 0), 0.0);
        }
    }
}

// ===========================================================================
// Scaffold-membrane harmonic bond
// ===========================================================================

class ScafCovHarmonicBond : public ::testing::Test
{
protected:
    void SetUp() override
    {
        configure_sheet(param, 64.0, 8.0);
        param.isEnergyHarmonicBondIncluded = true;
        mesh.reset(new ScafCovMesh(param));
        mesh->setup_flat();
        a = nearest_interior_vertex(*mesh, -4.0, 4.0);
        b = nearest_interior_vertex(*mesh, 12.0, -8.0);
        ASSERT_NE(a, b);
        const Matrix &va = mesh->vertices[a].coord;
        const Matrix &vb = mesh->vertices[b].coord;
        // p0 below vertex a (signed distance > 0, stretched beyond lbond),
        // p1 above vertex b (signed distance < 0)
        param.scaffoldingPoints = {vec3(va(0, 0) + 0.3, va(1, 0) - 0.2, va(2, 0) - 12.0),
                                   vec3(vb(0, 0) + 0.1, vb(1, 0) + 0.2, vb(2, 0) + 4.0)};
        mesh->set_scaffolding_vertices_correspondence();
    }

    double energy() { return mesh->calculate_scaffolding_energy_force(false); }

    Param param;
    std::unique_ptr<ScafCovMesh> mesh;
    int a = -1, b = -1;
};

TEST_F(ScafCovHarmonicBond, CorrespondencePicksTheNearestVertex)
{
    ASSERT_EQ(param.scaffoldingPoints_correspondingVertexIndex.size(), 2u);
    EXPECT_EQ(param.scaffoldingPoints_correspondingVertexIndex[0], a);
    EXPECT_EQ(param.scaffoldingPoints_correspondingVertexIndex[1], b);

    // A handful of arbitrary points: compare against a brute-force search.
    param.scaffoldingPoints = {vec3(3.3, -7.1, 2.0), vec3(-25.0, 20.5, -6.0), vec3(0.0, 0.0, 30.0),
                               vec3(17.9, 17.9, 0.1)};
    mesh->set_scaffolding_vertices_correspondence();
    for (size_t i = 0; i < param.scaffoldingPoints.size(); ++i)
    {
        const int chosen = param.scaffoldingPoints_correspondingVertexIndex[i];
        const Matrix &p = param.scaffoldingPoints[i];
        const double dChosen = mesh->get_squared_distance_sp_and_v(p, mesh->vertices[chosen]);
        const Matrix diff = mesh->vertices[chosen].coord - p;
        EXPECT_NEAR(dChosen, norm3(diff) * norm3(diff), 1e-9);
        for (const Vertex &v : mesh->vertices)
        {
            const Matrix d = v.coord - p;
            EXPECT_GE(norm3(d) * norm3(d), dChosen - 1e-9);
        }
    }
}

TEST_F(ScafCovHarmonicBond, VerboseCorrespondenceReportsSharedVertices)
{
    const Matrix &va = mesh->vertices[a].coord;
    param.scaffoldingPoints = {vec3(va(0, 0), va(1, 0), va(2, 0) - 9.0),
                               vec3(va(0, 0) + 0.5, va(1, 0), va(2, 0) - 9.0)};
    param.VERBOSE_MODE = true;
    testing::internal::CaptureStdout();
    mesh->set_scaffolding_vertices_correspondence();
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_NE(out.find("Repeated Vertex Index: " + std::to_string(a)), std::string::npos);
}

TEST_F(ScafCovHarmonicBond, EnergyMatchesTheSignedDistanceFormula)
{
    const double k = param.springConst;
    const double l = param.lbond;
    const double d0 = std::sqrt(0.3 * 0.3 + 0.2 * 0.2 + 12.0 * 12.0);
    const double d1 = -std::sqrt(0.1 * 0.1 + 0.2 * 0.2 + 4.0 * 4.0);
    const double expected = 0.5 * k * ((d0 - l) * (d0 - l) + (d1 - l) * (d1 - l));
    EXPECT_NEAR(energy(), expected, 1e-10 * expected);
    // The local-search flag does not change the result.
    EXPECT_NEAR(mesh->calculate_scaffolding_energy_force(true), expected, 1e-10 * expected);

    // Exactly at the rest length, on the positive side, the bond carries no energy or force.
    const Matrix &va = mesh->vertices[a].coord;
    param.scaffoldingPoints = {vec3(va(0, 0), va(1, 0), va(2, 0) - l)};
    param.scaffoldingPoints_correspondingVertexIndex = {a};
    EXPECT_NEAR(energy(), 0.0, 1e-20);
    EXPECT_NEAR(norm3(mesh->vertices[a].force.forceHarmonicBond), 0.0, 1e-12);
}

TEST_F(ScafCovHarmonicBond, ForceIsMinusTheEnergyGradientAndReactsOnTheScaffold)
{
    energy();
    const Matrix fa = mesh->vertices[a].force.forceHarmonicBond;
    const Matrix fb = mesh->vertices[b].force.forceHarmonicBond;
    EXPECT_GT(norm3(fa), 1.0);
    EXPECT_GT(norm3(fb), 1.0);
    // vertex a is pulled down toward its scaffold point; vertex b, overlapped by its point, is pushed up
    EXPECT_LT(fa(2, 0), 0.0);
    EXPECT_GT(fb(2, 0), 0.0);

    const double h = 1.0e-6;
    for (int vi : {a, b})
    {
        const Matrix f = mesh->vertices[vi].force.forceHarmonicBond;
        for (int axis = 0; axis < 3; ++axis)
        {
            const double x0 = mesh->vertices[vi].coord(axis, 0);
            mesh->vertices[vi].coord.set(axis, 0, x0 + h);
            const double ePlus = energy();
            mesh->vertices[vi].coord.set(axis, 0, x0 - h);
            const double eMinus = energy();
            mesh->vertices[vi].coord.set(axis, 0, x0);
            const double fd = -(ePlus - eMinus) / (2.0 * h);
            EXPECT_NEAR(f(axis, 0), fd, 1e-5 * std::max(1.0, std::abs(fd)))
                << "vertex " << vi << " axis " << axis;
        }
    }

    // Newton's third law on the scaffold side.
    energy();
    ASSERT_EQ(mesh->forceOnScaffoldingPoints.size(), 2u);
    for (int axis = 0; axis < 3; ++axis)
    {
        EXPECT_NEAR(mesh->forceOnScaffoldingPoints[0](axis, 0), -fa(axis, 0), 1e-12);
        EXPECT_NEAR(mesh->forceOnScaffoldingPoints[1](axis, 0), -fb(axis, 0), 1e-12);
        EXPECT_NEAR(mesh->forceTotalOnScaffolding(axis, 0), -(fa(axis, 0) + fb(axis, 0)), 1e-12);
    }
}

TEST_F(ScafCovHarmonicBond, ForceIsAccumulatedIntoTheTotalAndVerboseReportsIt)
{
    for (Vertex &v : mesh->vertices)
    {
        v.force.forceTotal = mat_calloc(3, 1);
    }
    param.VERBOSE_MODE = true;
    testing::internal::CaptureStdout();
    energy();
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_NE(out.find("Total force exerted on scaffolding"), std::string::npos);
    for (int vi : {a, b})
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            EXPECT_DOUBLE_EQ(mesh->vertices[vi].force.forceTotal(axis, 0),
                             mesh->vertices[vi].force.forceHarmonicBond(axis, 0));
        }
    }
}

TEST_F(ScafCovHarmonicBond, DisabledBondReturnsZeroAndTouchesNothing)
{
    param.isEnergyHarmonicBondIncluded = false;
    mesh->forceTotalOnScaffolding = vec3(1.0, 2.0, 3.0);
    EXPECT_EQ(mesh->calculate_scaffolding_energy_force(false), 0.0);
    EXPECT_EQ(mesh->forceTotalOnScaffolding(2, 0), 3.0);
}

TEST_F(ScafCovHarmonicBond, InsertionCurvatureMarksFacesAroundBondedVertices)
{
    param.insertCurv = 0.05;
    param.spontCurv = 0.01;
    // Pick a corner vertex too, whose neighbourhood includes ghost faces.
    std::vector<double> before(mesh->faces.size());
    for (size_t f = 0; f < mesh->faces.size(); ++f)
    {
        before[f] = mesh->faces[f].spontCurvature;
    }
    param.scaffoldingPoints_correspondingVertexIndex = {-1, a, 1000000, 0};
    mesh->set_scaffolding_insertion_curvature();

    std::vector<bool> nearBond(mesh->faces.size(), false);
    for (int vi : {a, 0})
    {
        for (int f : mesh->vertices[vi].adjacentFaces)
        {
            nearBond[f] = true;
        }
    }
    int nGhostNear = 0;
    for (size_t f = 0; f < mesh->faces.size(); ++f)
    {
        const Face &face = mesh->faces[f];
        if (face.isGhost)
        {
            EXPECT_DOUBLE_EQ(face.spontCurvature, before[f]); // ghosts never change
            nGhostNear += nearBond[f] ? 1 : 0;
        }
        else if (nearBond[f] || face.isInsertionPatch)
        {
            EXPECT_DOUBLE_EQ(face.spontCurvature, 0.05);
        }
        else
        {
            EXPECT_DOUBLE_EQ(face.spontCurvature, 0.01);
        }
    }
    EXPECT_GT(nGhostNear, 0);

    // Without the harmonic bond only the insertion-patch rule applies.
    param.isEnergyHarmonicBondIncluded = false;
    mesh->set_scaffolding_insertion_curvature();
    for (const Face &face : mesh->faces)
    {
        if (!face.isGhost && !face.isInsertionPatch)
        {
            EXPECT_DOUBLE_EQ(face.spontCurvature, 0.01);
        }
    }
}

TEST_F(ScafCovHarmonicBond, RigidPropagationTranslatesTheWholeScaffoldDownhill)
{
    param.propagateScaffoldingNstep = 4;
    const std::vector<Matrix> start = param.scaffoldingPoints;
    const double e0 = energy();
    testing::internal::CaptureStdout();
    ASSERT_TRUE(mesh->propagate_scaffolding());
    testing::internal::GetCapturedStdout();
    const double e1 = energy();
    EXPECT_LT(e1, e0);

    const Matrix &move = mesh->scaffoldingMovementVector;
    EXPECT_GT(norm3(move), 0.0);
    for (size_t i = 0; i < start.size(); ++i)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            EXPECT_NEAR(param.scaffoldingPoints[i](axis, 0) - start[i](axis, 0), move(axis, 0), 1e-12);
        }
    }
}

// ===========================================================================
// Reference-state / reaction parsing and topology
// ===========================================================================

class ScafCovGagTopology : public ::testing::Test
{
protected:
    void SetUp() override
    {
        refFile = tmp.file("ring.dat");
        write_reference_state(refFile, ring.mols, ring.bonds, true);
        reactionFile = tmp.file("gag.inp");
        write_text(reactionFile, kReactionFile);
    }

    TempDir tmp;
    RingLattice ring;
    std::string refFile;
    std::string reactionFile;
};

TEST_F(ScafCovGagTopology, ModesOffDoNothing)
{
    Param p;
    p.scaffoldingPoints = ring.points();
    Mesh mesh(p);
    EXPECT_FALSE(mesh.initialize_gag_scaffolding_topology());
    EXPECT_FALSE(mesh.pre_relax_gag_scaffolding());
    EXPECT_EQ(mesh.calculate_gag_scaffolding_internal_energy(), 0.0);
    EXPECT_EQ(mesh.calculate_gag_scaffolding_internal_energy(p.scaffoldingPoints), 0.0);
    for (const auto &forces : {mesh.calculate_gag_scaffolding_forces_fd(),
                               mesh.calculate_gag_subunit_translation_forces_fd(),
                               mesh.calculate_gag_subunit_rotation_torques_fd()})
    {
        ASSERT_EQ(forces.size(), 3u);
        for (const Matrix &f : forces)
        {
            EXPECT_EQ(norm3(f), 0.0);
        }
    }
    EXPECT_FALSE(mesh.gagScaffoldingTopologyInitialized);
}

TEST_F(ScafCovGagTopology, ConfigurationErrorsThrow)
{
    {
        Param p;
        configure_lattice(p, refFile, true);
        p.isGagScaffoldingEnergyIncluded = true; // both modes
        p.scaffoldingPoints = ring.points();
        Mesh mesh(p);
        EXPECT_THROW(mesh.initialize_gag_scaffolding_topology(), std::runtime_error);
    }
    for (bool idealized : {true, false})
    {
        Param p;
        configure_lattice(p, "", idealized); // empty reference file name
        p.scaffoldingPoints = ring.points();
        Mesh mesh(p);
        EXPECT_THROW(mesh.initialize_gag_scaffolding_topology(), std::runtime_error);
    }
    {
        Param p;
        configure_lattice(p, tmp.file("does_not_exist.dat"), true);
        p.scaffoldingPoints = ring.points();
        Mesh mesh(p);
        EXPECT_THROW(mesh.initialize_gag_scaffolding_topology(), std::runtime_error);
    }
    {
        // molecules but no bound states: the parser reports failure
        const std::string noBonds = tmp.file("no_bonds.dat");
        write_reference_state(noBonds, ring.mols, {});
        Param p;
        configure_lattice(p, noBonds, true);
        p.scaffoldingPoints = ring.points();
        Mesh mesh(p);
        EXPECT_THROW(mesh.initialize_gag_scaffolding_topology(), std::runtime_error);
    }
    {
        // reference count does not match the scaffold
        Param p;
        configure_lattice(p, refFile, true);
        p.scaffoldingPoints = ring.points();
        p.scaffoldingPoints.pop_back();
        Mesh mesh(p);
        EXPECT_THROW(mesh.initialize_gag_scaffolding_topology(), std::runtime_error);
    }
    {
        // Gag mode with a reaction file that cannot be opened
        Param p;
        configure_lattice(p, refFile, false);
        p.gagReactionFileName = tmp.file("missing.inp");
        p.scaffoldingPoints = ring.points();
        Mesh mesh(p);
        EXPECT_THROW(mesh.initialize_gag_scaffolding_topology(), std::runtime_error);
    }
}

TEST_F(ScafCovGagTopology, IdealizedTopologyMapsMoleculesAndIsRelaxedAtReference)
{
    Param p;
    configure_lattice(p, refFile, true);
    p.scaffoldingPoints = ring.points();
    Mesh mesh(p);
    // before initialisation the explicit overload reports nothing
    EXPECT_EQ(mesh.calculate_gag_scaffolding_internal_energy(p.scaffoldingPoints, mesh.gagSubunits), 0.0);

    testing::internal::CaptureStdout();
    ASSERT_TRUE(mesh.initialize_gag_scaffolding_topology());
    const std::string out = testing::internal::GetCapturedStdout();
    // the interfaces are offset from the COM-COM lines
    EXPECT_NE(out.find("nearly colinear"), std::string::npos);
    EXPECT_TRUE(mesh.initialize_gag_scaffolding_topology()); // idempotent

    ASSERT_EQ(mesh.gagSubunits.size(), 3u);
    EXPECT_EQ(mesh.gagSubunits[0].moleculeId, 2);
    EXPECT_EQ(mesh.gagSubunits[1].moleculeId, 0);
    EXPECT_EQ(mesh.gagSubunits[2].moleculeId, 1);
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_EQ(mesh.gagSubunits[i].pointIndex, i);
        EXPECT_LT(orthogonality_defect(mesh.gagSubunits[i].rotation), 1e-15); // identity
        EXPECT_DOUBLE_EQ(mesh.gagSubunits[i].rotation(0, 0), 1.0);
    }
    // The first reference vector of each molecule wins; molecule 2 has none.
    EXPECT_TRUE(mesh.gagSubunits[1].hasReferenceVector);
    EXPECT_DOUBLE_EQ(mesh.gagSubunits[1].localReferenceVector(2, 0), 1.0); // sphere_normal
    EXPECT_DOUBLE_EQ(mesh.gagSubunits[1].localReferenceVector(0, 0), 0.0);
    EXPECT_TRUE(mesh.gagSubunits[2].hasReferenceVector);
    EXPECT_DOUBLE_EQ(mesh.gagSubunits[2].localReferenceVector(1, 0), 1.0); // tangent
    EXPECT_DOUBLE_EQ(mesh.gagSubunits[2].localReferenceVector(2, 0), 0.0);
    EXPECT_FALSE(mesh.gagSubunits[0].hasReferenceVector);
    // Local interfaces are relative to the reference COM.
    const auto ifb = ring.iface(0, "b");
    EXPECT_NEAR(mesh.gagSubunits[1].localInterfaces.at("b")(0, 0), ifb[0], 1e-14);
    EXPECT_NEAR(mesh.gagSubunits[2].localInterfaces.at("a")(1, 0), ring.iface(1, "a")[1] - 0.0, 1e-14);

    // Three distinct bonds: the duplicate and the malformed block are dropped.
    ASSERT_EQ(mesh.gagBonds.size(), 3u);
    ASSERT_EQ(mesh.gagInteractions.size(), 3u);
    const double expectedSigma01 = distance(ring.iface(0, "b"), ring.iface(1, "a"));
    EXPECT_EQ(mesh.gagInteractions[0].subunit1, 1);
    EXPECT_EQ(mesh.gagInteractions[0].subunit2, 2);
    EXPECT_NEAR(mesh.gagInteractions[0].expectedSigma, expectedSigma01, 1e-12);
    EXPECT_NEAR(mesh.gagBonds[0].expectedLength, distance(ring.mols[0].com, ring.mols[1].com), 1e-12);
    // Angles are taken from the geometry; only the 0-1 pair has two reference vectors.
    EXPECT_GT(mesh.gagInteractions[0].assocAngles[0], 0.0);
    EXPECT_LT(mesh.gagInteractions[0].assocAngles[0], M_PI);
    EXPECT_NE(mesh.gagInteractions[0].assocAngles[4], 0.0);
    EXPECT_EQ(mesh.gagInteractions[1].assocAngles[4], 0.0);
    EXPECT_EQ(mesh.gagInteractions[2].assocAngles[4], 0.0);

    // At the reference geometry every term vanishes.
    EXPECT_NEAR(mesh.calculate_gag_scaffolding_internal_energy(), 0.0, 1e-20);
    EXPECT_NEAR(mesh.calculate_gag_scaffolding_internal_energy(p.scaffoldingPoints), 0.0, 1e-20);
}

TEST_F(ScafCovGagTopology, GagModeTakesTargetsFromTheReactionFile)
{
    Param p;
    configure_lattice(p, refFile, false);
    p.gagReactionFileName = reactionFile;
    p.VERBOSE_MODE = true;
    p.scaffoldingPoints = ring.points();
    Mesh mesh(p);
    testing::internal::CaptureStdout();
    ASSERT_TRUE(mesh.initialize_gag_scaffolding_topology());
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_NE(out.find("Missing reaction target for c|a"), std::string::npos);
    EXPECT_EQ(out.find("nearly colinear"), std::string::npos); // idealized-only warning

    ASSERT_EQ(mesh.gagInteractions.size(), 3u);
    // b|a bonds: 0.b-1.a and 1.b-2.a
    for (int k : {0, 1})
    {
        const Mesh::GagInteraction &it = mesh.gagInteractions[k];
        EXPECT_DOUBLE_EQ(it.expectedSigma, 2.5);
        EXPECT_DOUBLE_EQ(it.assocAngles[0], 1.1);
        EXPECT_DOUBLE_EQ(it.assocAngles[1], 2.2);
        EXPECT_DOUBLE_EQ(it.assocAngles[2], 0.3);
        EXPECT_DOUBLE_EQ(it.assocAngles[3], 0.5);
        EXPECT_DOUBLE_EQ(it.assocAngles[4], -0.7);
        const Mesh::GagBond &bond = mesh.gagBonds[k];
        EXPECT_DOUBLE_EQ(bond.thetaAtPoint1, 1.1);
        EXPECT_DOUBLE_EQ(bond.thetaAtPoint2, 2.2);
        EXPECT_DOUBLE_EQ(bond.expectedPhi, 0.4);
        EXPECT_DOUBLE_EQ(bond.expectedOmega, -0.7);
    }
    // c|a has no reaction (the c|c block had no sigma): geometric fallback
    const double sigma20 = distance(ring.iface(2, "c"), ring.iface(0, "a"));
    EXPECT_NEAR(mesh.gagInteractions[2].expectedSigma, sigma20, 1e-12);

    // Bond-length term alone, checked by hand at the reference geometry.
    p.gagKtheta = p.gagKphi = p.gagKomega = 0.0;
    const double s01 = distance(ring.iface(0, "b"), ring.iface(1, "a"));
    const double s12 = distance(ring.iface(1, "b"), ring.iface(2, "a"));
    const double expected = 0.5 * ((s01 - 2.5) * (s01 - 2.5) + (s12 - 2.5) * (s12 - 2.5));
    EXPECT_NEAR(mesh.calculate_gag_scaffolding_internal_energy(), expected, 1e-12 * expected);

    // With every constant on, the energy is positive (targets differ from the geometry).
    p.gagKtheta = 0.7;
    p.gagKphi = 0.3;
    p.gagKomega = 0.5;
    EXPECT_GT(mesh.calculate_gag_scaffolding_internal_energy(), expected);
}

TEST_F(ScafCovGagTopology, GagModeWithoutReactionFileUsesGeometry)
{
    Param p;
    configure_lattice(p, refFile, false);
    p.scaffoldingPoints = ring.points();
    Mesh mesh(p);
    // calling the energy initialises the topology on demand
    const double e = mesh.calculate_gag_scaffolding_internal_energy();
    EXPECT_TRUE(mesh.gagScaffoldingTopologyInitialized);
    ASSERT_EQ(mesh.gagInteractions.size(), 3u);
    // sigma and theta targets match the geometry; the torsion targets are zero,
    // so a generic (non-planar) ring still carries torsion energy.
    p.gagKphi = p.gagKomega = 0.0;
    EXPECT_NEAR(mesh.calculate_gag_scaffolding_internal_energy(), 0.0, 1e-20);
    EXPECT_GT(e, 0.0);
}

TEST_F(ScafCovGagTopology, WrittenStateRoundTripsAsAReference)
{
    Param p;
    configure_lattice(p, refFile, true);
    p.scaffoldingPoints = ring.points();
    Mesh mesh(p);
    testing::internal::CaptureStdout();
    mesh.initialize_gag_scaffolding_topology();
    perturb(p.scaffoldingPoints);
    const std::string state = tmp.file("state.dat");
    mesh.write_gag_scaffolding_state_dat(state);
    testing::internal::GetCapturedStdout();
    ASSERT_GT(mesh.calculate_gag_scaffolding_internal_energy(), 0.0);

    Param p2;
    configure_lattice(p2, state, true);
    p2.scaffoldingPoints = p.scaffoldingPoints;
    Mesh mesh2(p2);
    testing::internal::CaptureStdout();
    ASSERT_TRUE(mesh2.initialize_gag_scaffolding_topology());
    testing::internal::GetCapturedStdout();
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_EQ(mesh2.gagSubunits[i].moleculeId, mesh.gagSubunits[i].moleculeId);
    }
    ASSERT_EQ(mesh2.gagInteractions.size(), mesh.gagInteractions.size());
    EXPECT_NEAR(mesh2.calculate_gag_scaffolding_internal_energy(), 0.0, 1e-20);
}

TEST_F(ScafCovGagTopology, InterfaceAtTheComGivesFiniteZeroAngles)
{
    // Two subunits, all on one line, one interface sitting on its own COM:
    // the angle and dihedral helpers must return 0 instead of NaN.
    const std::string file = tmp.file("degenerate.dat");
    std::vector<RefMol> mols = {{0, {{0.0, 0.0, 0.0}}, {{"a", {{0.0, 0.0, 0.0}}, "", {{0, 0, 0}}}}},
                                {1, {{5.0, 0.0, 0.0}}, {{"a", {{3.0, 0.0, 0.0}}, "", {{0, 0, 0}}}}}};
    write_reference_state(file, mols, {{0, "a", 1, "a"}});
    Param p;
    configure_lattice(p, file, false);
    p.scaffoldingPoints = {vec3(0.0, 0.0, 0.0), vec3(5.0, 0.0, 0.0)};
    Mesh mesh(p);
    ASSERT_TRUE(mesh.initialize_gag_scaffolding_topology());
    ASSERT_EQ(mesh.gagInteractions.size(), 1u);
    EXPECT_EQ(mesh.gagInteractions[0].assocAngles[0], 0.0);       // zero-length lever
    EXPECT_NEAR(mesh.gagInteractions[0].assocAngles[1], M_PI, 1e-12); // straight line
    EXPECT_DOUBLE_EQ(mesh.gagInteractions[0].expectedSigma, 3.0);
    EXPECT_NEAR(mesh.calculate_gag_scaffolding_internal_energy(), 0.0, 1e-20);
    p.scaffoldingPoints[1].set(1, 0, 0.5);
    const double e = mesh.calculate_gag_scaffolding_internal_energy();
    EXPECT_TRUE(std::isfinite(e));
    EXPECT_GT(e, 0.0);
}

// ===========================================================================
// Internal energy: invariances and finite-difference forces / torques
// ===========================================================================

class ScafCovGagForces : public ::testing::TestWithParam<bool>
{
protected:
    void SetUp() override
    {
        refFile = tmp.file("ring.dat");
        write_reference_state(refFile, ring.mols, ring.bonds);
        configure_lattice(param, refFile, GetParam());
        param.scaffoldingPoints = ring.points();
        mesh.reset(new Mesh(param));
        testing::internal::CaptureStdout();
        mesh->initialize_gag_scaffolding_topology();
        testing::internal::GetCapturedStdout();
        perturb(param.scaffoldingPoints);
        // rotate one subunit a little so the orientations are not all the same
        mesh->gagSubunits[1].rotation = rotation({{0.3, -0.5, 0.8}}, 0.2);
    }

    TempDir tmp;
    RingLattice ring;
    std::string refFile;
    Param param;
    std::unique_ptr<Mesh> mesh;
};

TEST_P(ScafCovGagForces, EnergyIsInvariantUnderRigidMotion)
{
    const double e0 = mesh->calculate_gag_scaffolding_internal_energy();
    ASSERT_GT(e0, 1e-3);
    const Matrix q = rotation({{1.0, 2.0, -0.5}}, 0.9);
    const Matrix t = vec3(3.0, -7.0, 11.0);
    std::vector<Matrix> moved;
    for (const Matrix &pnt : param.scaffoldingPoints)
    {
        moved.push_back(q * pnt + t);
    }
    std::vector<Mesh::GagSubunit> rotated = mesh->gagSubunits;
    for (Mesh::GagSubunit &s : rotated)
    {
        s.rotation = q * s.rotation;
    }
    EXPECT_NEAR(mesh->calculate_gag_scaffolding_internal_energy(moved, rotated), e0, 1e-10 * e0);
    // translation alone keeps the energy too
    std::vector<Matrix> shifted;
    for (const Matrix &pnt : param.scaffoldingPoints)
    {
        shifted.push_back(pnt + t);
    }
    EXPECT_NEAR(mesh->calculate_gag_scaffolding_internal_energy(shifted), e0, 1e-10 * e0);
}

TEST_P(ScafCovGagForces, FiniteDifferenceForcesAndTorquesBalance)
{
    const std::vector<Matrix> f = mesh->calculate_gag_scaffolding_forces_fd();
    const std::vector<Matrix> ft = mesh->calculate_gag_subunit_translation_forces_fd();
    const std::vector<Matrix> tau = mesh->calculate_gag_subunit_rotation_torques_fd();
    ASSERT_EQ(f.size(), 3u);
    ASSERT_EQ(ft.size(), 3u);
    ASSERT_EQ(tau.size(), 3u);

    double fMax = 0.0, tauMax = 0.0, rMax = 0.0;
    Matrix fSum(3, 1, true), torqueSum(3, 1, true);
    for (int i = 0; i < 3; ++i)
    {
        fMax = std::max(fMax, norm3(f[i]));
        tauMax = std::max(tauMax, norm3(tau[i]));
        rMax = std::max(rMax, norm3(param.scaffoldingPoints[i]));
        for (int a = 0; a < 3; ++a)
        {
            // the two FD force routines perturb the same coordinates
            EXPECT_NEAR(f[i](a, 0), ft[i](a, 0), 1e-9 * std::max(1.0, fMax));
        }
        fSum += f[i];
        torqueSum += cross3(param.scaffoldingPoints[i], f[i]);
        torqueSum += tau[i];
    }
    ASSERT_GT(fMax, 1e-3);
    ASSERT_GT(tauMax, 1e-3);
    // translation invariance: no net force; rotation invariance: no net torque
    EXPECT_LT(norm3(fSum), 1e-6 * fMax);
    EXPECT_LT(norm3(torqueSum), 1e-6 * (fMax * rMax + tauMax));

    // A small step along the force lowers the energy.
    const double e0 = mesh->calculate_gag_scaffolding_internal_energy();
    std::vector<Matrix> stepped = param.scaffoldingPoints;
    for (int i = 0; i < 3; ++i)
    {
        stepped[i] += 1.0e-4 * f[i];
    }
    EXPECT_LT(mesh->calculate_gag_scaffolding_internal_energy(stepped), e0);
}

TEST_P(ScafCovGagForces, NonFiniteSubunitIsSkipped)
{
    // Bond-length term only, so the surviving interaction is easy to check by hand.
    param.gagKtheta = param.gagKphi = param.gagKomega = 0.0;
    std::vector<Matrix> pts = param.scaffoldingPoints;
    pts[0].set(0, 0, std::numeric_limits<double>::infinity()); // subunit of molecule 2
    const double e = mesh->calculate_gag_scaffolding_internal_energy(pts);
    ASSERT_TRUE(std::isfinite(e));
    // Only the molecule 0 - molecule 1 interaction (points 1 and 2) survives.
    const Mesh::GagInteraction &it = mesh->gagInteractions[0];
    ASSERT_EQ(it.subunit1, 1);
    ASSERT_EQ(it.subunit2, 2);
    const Mesh::GagSubunit &s1 = mesh->gagSubunits[1];
    const Mesh::GagSubunit &s2 = mesh->gagSubunits[2];
    const Matrix if1 = pts[1] + s1.rotation * s1.localInterfaces.at(it.interface1);
    const Matrix if2 = pts[2] + s2.rotation * s2.localInterfaces.at(it.interface2);
    const double ds = norm3(if2 - if1) - it.expectedSigma;
    EXPECT_NEAR(e, 0.5 * param.gagKsigma * ds * ds, 1e-12);
}

TEST_P(ScafCovGagForces, BelowResolutionStepGivesZeroTorque)
{
    // An FD step below 1e-16 rad makes the rotation helper return the identity.
    param.gagFiniteDifferenceStep = 1.0e-17;
    const std::vector<Matrix> tau = mesh->calculate_gag_subunit_rotation_torques_fd();
    for (const Matrix &t : tau)
    {
        EXPECT_EQ(norm3(t), 0.0);
    }
}

TEST_P(ScafCovGagForces, PreRelaxationLowersTheInternalEnergy)
{
    param.gagPreRelaxSteps = 0;
    EXPECT_FALSE(mesh->pre_relax_gag_scaffolding());

    param.gagPreRelaxSteps = 15;
    param.gagPropagationStepSize = 0.02;
    const double e0 = mesh->calculate_gag_scaffolding_internal_energy();
    const std::vector<Matrix> before = param.scaffoldingPoints;
    testing::internal::CaptureStdout();
    ASSERT_TRUE(mesh->pre_relax_gag_scaffolding());
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_NE(out.find("Steps: 15"), std::string::npos);
    const double e1 = mesh->calculate_gag_scaffolding_internal_energy();
    EXPECT_LT(e1, 0.9 * e0);
    double moved = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        moved += norm3(param.scaffoldingPoints[i] - before[i]);
        EXPECT_LT(orthogonality_defect(mesh->gagSubunits[i].rotation), 1e-12);
    }
    EXPECT_GT(moved, 0.0);
}

INSTANTIATE_TEST_SUITE_P(ScafCovModes, ScafCovGagForces, ::testing::Values(true, false),
                         [](const ::testing::TestParamInfo<bool> &info) {
                             return info.param ? "Idealized" : "Gag";
                         });

// ===========================================================================
// Scaffold coupled to the membrane in lattice mode
// ===========================================================================

TEST(ScafCovLatticeCoupling, EnergyIncludesInternalTermAndPropagationGoesDownhill)
{
    for (bool idealized : {true, false})
    {
        SCOPED_TRACE(idealized ? "idealized" : "gag");
        TempDir tmp;
        RingLattice ring;
        const std::string refFile = tmp.file("ring.dat");
        write_reference_state(refFile, ring.mols, ring.bonds);

        Param p;
        configure_sheet(p, 64.0, 8.0);
        configure_lattice(p, refFile, idealized);
        p.isEnergyHarmonicBondIncluded = true;
        p.scaffoldingPoints = ring.points();
        Mesh mesh(p);
        mesh.setup_flat();
        for (Vertex &v : mesh.vertices)
        {
            v.coord.set(2, 0, 10.0); // the membrane above the scaffold
        }
        testing::internal::CaptureStdout();
        mesh.initialize_gag_scaffolding_topology();
        testing::internal::GetCapturedStdout();
        perturb(p.scaffoldingPoints);
        mesh.set_scaffolding_vertices_correspondence();

        const double internal = mesh.calculate_gag_scaffolding_internal_energy();
        ASSERT_GT(internal, 0.0);
        const double total = mesh.calculate_scaffolding_energy_force(false);
        if (idealized)
        {
            EXPECT_NEAR(p.energy.energyIdealizedProteinLattice, internal, 1e-12 * internal);
        }
        else
        {
            EXPECT_NEAR(p.energy.energyGagScaffolding, internal, 1e-12 * internal);
        }
        // the harmonic part on its own
        p.isGagScaffoldingEnergyIncluded = false;
        p.isIdealizedProteinLatticeEnergyIncluded = false;
        const double harmonic = mesh.calculate_scaffolding_energy_force(false);
        p.isGagScaffoldingEnergyIncluded = !idealized;
        p.isIdealizedProteinLatticeEnergyIncluded = idealized;
        EXPECT_NEAR(total, harmonic + internal, 1e-12 * total);

        p.propagateScaffoldingNstep = 3;
        p.gagPropagationStepSize = 0.01;
        const std::vector<Matrix> before = p.scaffoldingPoints;
        testing::internal::CaptureStdout();
        ASSERT_TRUE(mesh.propagate_scaffolding());
        testing::internal::GetCapturedStdout();
        const double after = mesh.calculate_scaffolding_energy_force(false);
        EXPECT_LT(after, total);

        // the recorded movement is the average displacement
        Matrix average(3, 1, true);
        for (int i = 0; i < 3; ++i)
        {
            average += p.scaffoldingPoints[i] - before[i];
            EXPECT_LT(orthogonality_defect(mesh.gagSubunits[i].rotation), 1e-12);
        }
        for (int a = 0; a < 3; ++a)
        {
            EXPECT_NEAR(mesh.scaffoldingMovementVector(a, 0), average(a, 0) / 3.0, 1e-12);
        }
    }
}

// ===========================================================================
// Orientation of the imported scaffold
// ===========================================================================

namespace
{

std::vector<Matrix> upright_cap(double R, double alpha, const Matrix &apex)
{
    std::vector<Matrix> pts;
    pts.push_back(Matrix(apex));
    for (int ring = 1; ring <= 2; ++ring)
    {
        const double a = alpha * ring / 2.0;
        for (int k = 0; k < 8; ++k)
        {
            const double phi = 2.0 * M_PI * k / 8.0 + 0.1 * ring;
            pts.push_back(apex + vec3(R * std::sin(a) * std::cos(phi), R * std::sin(a) * std::sin(phi),
                                      -R + R * std::cos(a)));
        }
    }
    return pts;
}

void expect_rigid(const std::vector<Matrix> &a, const std::vector<Matrix> &b)
{
    ASSERT_EQ(a.size(), b.size());
    for (size_t i = 0; i < a.size(); ++i)
    {
        for (size_t j = i + 1; j < a.size(); ++j)
        {
            EXPECT_NEAR(norm3(a[i] - a[j]), norm3(b[i] - b[j]), 1e-9);
        }
    }
}

} // namespace

TEST(ScafCovOrient, TooFewPointsIsRejected)
{
    Param p;
    p.scaffoldingPoints = {vec3(0, 0, 0), vec3(1, 0, 0)};
    Mesh mesh(p);
    EXPECT_FALSE(mesh.orient_scaffolding_plane_to_membrane());
    EXPECT_LT(orthogonality_defect(mesh.gagInitialAlignmentRotation), 1e-15);
    EXPECT_DOUBLE_EQ(mesh.gagInitialAlignmentRotation(1, 1), 1.0);
}

TEST(ScafCovOrient, TiltedCapIsTurnedToFaceTheMembrane)
{
    const double R = 20.0, alpha = 0.5;
    Param p;
    p.scaffoldingSphereRaidus = R;
    p.lbond = 9.0;
    p.scaffoldingZeroPlaneZ = 2.0;
    const Matrix sphereCenter = vec3(5.0, -3.0, 1.0);
    const std::vector<Matrix> upright = upright_cap(R, alpha, sphereCenter + vec3(0, 0, R));
    const std::array<double, 3> tiltAxis{{1.0, 0.4, 0.0}};
    const Matrix tilt = rotation(tiltAxis, 0.3);
    std::vector<Matrix> tilted;
    for (const Matrix &pt : upright)
    {
        tilted.push_back(tilt * (pt - sphereCenter) + sphereCenter);
    }
    p.scaffoldingPoints = tilted;
    Mesh mesh(p);
    testing::internal::CaptureStdout();
    ASSERT_TRUE(mesh.orient_scaffolding_plane_to_membrane());
    testing::internal::GetCapturedStdout();

    const std::vector<Matrix> &out = p.scaffoldingPoints;
    expect_rigid(out, tilted);
    // the cap axis is vertical again: each ring has a single height, the apex on top at x = y = 0
    EXPECT_NEAR(out[0](0, 0), 0.0, 1e-9);
    EXPECT_NEAR(out[0](1, 0), 0.0, 1e-9);
    for (int ring = 0; ring < 2; ++ring)
    {
        for (int k = 1; k < 8; ++k)
        {
            EXPECT_NEAR(out[1 + 8 * ring + k](2, 0), out[1 + 8 * ring](2, 0), 1e-9);
        }
        EXPECT_LT(out[1 + 8 * ring](2, 0), out[0](2, 0));
    }
    // the inferred cap edge of the membrane starts on the zero plane
    const double capRadius = R * std::sin(alpha);
    const double rm2 = capRadius * capRadius * (R + 9.0) * (R + 9.0) / (R * R);
    EXPECT_NEAR(out[0](2, 0) - R + std::sqrt((R + 9.0) * (R + 9.0) - rm2), 2.0, 1e-9);
    // the alignment rotation maps the tilted axis onto +z and is a proper rotation
    const Matrix axis = tilt * vec3(0, 0, 1);
    const Matrix mapped = mesh.gagInitialAlignmentRotation * axis;
    EXPECT_NEAR(mapped(2, 0), 1.0, 1e-9);
    EXPECT_LT(orthogonality_defect(mesh.gagInitialAlignmentRotation), 1e-12);
}

TEST(ScafCovOrient, UprightOrFlatScaffoldStaysHorizontal)
{
    // Already upright: the rotation is the identity, only a shift is applied.
    {
        Param p;
        p.scaffoldingSphereRaidus = 20.0;
        p.scaffoldingPoints = upright_cap(20.0, 0.5, vec3(1.0, 2.0, 30.0));
        const std::vector<Matrix> before = p.scaffoldingPoints;
        Mesh mesh(p);
        testing::internal::CaptureStdout();
        ASSERT_TRUE(mesh.orient_scaffolding_plane_to_membrane());
        testing::internal::GetCapturedStdout();
        expect_rigid(p.scaffoldingPoints, before);
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                EXPECT_NEAR(mesh.gagInitialAlignmentRotation(i, j), i == j ? 1.0 : 0.0, 1e-9);
            }
        }
        const double dz = p.scaffoldingPoints[0](2, 0) - before[0](2, 0);
        for (size_t i = 0; i < before.size(); ++i)
        {
            EXPECT_NEAR(p.scaffoldingPoints[i](0, 0), before[i](0, 0) - 1.0, 1e-9);
            EXPECT_NEAR(p.scaffoldingPoints[i](1, 0), before[i](1, 0) - 2.0, 1e-9);
            EXPECT_NEAR(p.scaffoldingPoints[i](2, 0) - before[i](2, 0), dz, 1e-9);
        }
    }
    // A flat sheet of points (normal sign undetermined): it stays horizontal either way.
    for (double z : {4.0, -3.0})
    {
        Param p;
        p.scaffoldingSphereRaidus = 0.0; // invalid: 50 nm default
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                p.scaffoldingPoints.push_back(vec3(3.0 * i + 0.2 * j, 2.0 * j, z));
            }
        }
        const std::vector<Matrix> before = p.scaffoldingPoints;
        Mesh mesh(p);
        testing::internal::CaptureStdout();
        ASSERT_TRUE(mesh.orient_scaffolding_plane_to_membrane());
        testing::internal::GetCapturedStdout();
        expect_rigid(p.scaffoldingPoints, before);
        for (const Matrix &pt : p.scaffoldingPoints)
        {
            EXPECT_NEAR(pt(2, 0), p.scaffoldingPoints[0](2, 0), 1e-9);
        }
        EXPECT_LT(orthogonality_defect(mesh.gagInitialAlignmentRotation), 1e-12);
    }
}

TEST(ScafCovOrient, CapWiderThanItsSphereIsClamped)
{
    // Sphere radius smaller than the cap radius: the cap radius is clamped to R,
    // so the inferred edge is at the sphere's equator: apex - R + sqrt((R+l)^2 - (R+l)^2) ~ apex - R.
    const double R = 3.0;
    Param p;
    p.scaffoldingSphereRaidus = R;
    p.lbond = 9.0;
    p.scaffoldingZeroPlaneZ = 0.0;
    p.scaffoldingPoints = upright_cap(20.0, 0.5, vec3(0.0, 0.0, 0.0));
    Mesh mesh(p);
    testing::internal::CaptureStdout();
    ASSERT_TRUE(mesh.orient_scaffolding_plane_to_membrane());
    testing::internal::GetCapturedStdout();
    const double rc = R * (1.0 - 1.0e-6);
    const double edge = std::sqrt(std::max(0.0, (R + 9.0) * (R + 9.0) * (1.0 - rc * rc / (R * R))));
    EXPECT_NEAR(p.scaffoldingPoints[0](2, 0) - R + edge, 0.0, 1e-9);
}

TEST(ScafCovOrient, NearlyFlatSheetWithAmbiguousNormalStaysHorizontal)
{
    // Tilts of 1e-8 leave the apex test (normal . (apex - centroid)) at round-off
    // level, so the sign of the best-fit normal is decided by round-off; on macOS
    // arm64 several of these end up anti-parallel to +z and take the half-turn
    // branch. Either way the result must be a rigid motion onto a horizontal plane.
    const double tilts[][3] = {{-2e-8, -1.5e-8, 0}, {-1e-8, -5e-9, 0}, {0.0, 5e-9, 0}, {0.0, -5e-9, 1}};
    for (const auto &tilt : tilts)
    {
        Param p;
        p.scaffoldingSphereRaidus = 20.0;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                double x = 3.0 * i + 0.2 * j, y = 2.0 * j;
                if (tilt[2] != 0.0)
                {
                    std::swap(x, y);
                }
                p.scaffoldingPoints.push_back(vec3(x, y, 4.0 + tilt[0] * x + tilt[1] * y));
            }
        }
        const std::vector<Matrix> before = p.scaffoldingPoints;
        Mesh mesh(p);
        testing::internal::CaptureStdout();
        ASSERT_TRUE(mesh.orient_scaffolding_plane_to_membrane());
        testing::internal::GetCapturedStdout();
        expect_rigid(p.scaffoldingPoints, before);
        for (const Matrix &pt : p.scaffoldingPoints)
        {
            EXPECT_NEAR(pt(2, 0), p.scaffoldingPoints[0](2, 0), 1e-6);
        }
        EXPECT_LT(orthogonality_defect(mesh.gagInitialAlignmentRotation), 1e-12);
        EXPECT_NEAR(std::abs(mesh.gagInitialAlignmentRotation(2, 2)), 1.0, 1e-9);
    }
}

TEST(ScafCovOrient, ApexListedLastIsStillCentred)
{
    Param p;
    p.scaffoldingSphereRaidus = 20.0;
    p.scaffoldingZeroPlaneZ = -1.0;
    std::vector<Matrix> pts = upright_cap(20.0, 0.5, vec3(-6.0, 4.0, 12.0));
    std::rotate(pts.begin(), pts.begin() + 1, pts.end()); // apex now last
    p.scaffoldingPoints = pts;
    Mesh mesh(p);
    testing::internal::CaptureStdout();
    ASSERT_TRUE(mesh.orient_scaffolding_plane_to_membrane());
    testing::internal::GetCapturedStdout();
    const Matrix &apex = p.scaffoldingPoints.back();
    EXPECT_NEAR(apex(0, 0), 0.0, 1e-9);
    EXPECT_NEAR(apex(1, 0), 0.0, 1e-9);
    for (size_t i = 0; i + 1 < p.scaffoldingPoints.size(); ++i)
    {
        EXPECT_LT(p.scaffoldingPoints[i](2, 0), apex(2, 0));
    }
}

// ===========================================================================
// Sphere centre / cap radius corner cases
// ===========================================================================

TEST(ScafCovSphereCenter, DefaultEmptyAndUnorderedScaffolds)
{
    Param p;
    p.scaffoldingSphereRaidus = 15.0;
    Mesh mesh(p);
    // no scaffold: centre at the origin, cap radius zero
    EXPECT_EQ(norm3(mesh.find_center_of_scaffolding_sphere(false)), 0.0);
    EXPECT_EQ(mesh.approximate_scaffolding_cap_radius(false), 0.0);

    // highest point listed last: the centre is R below it
    p.scaffoldingPoints = {vec3(4.0, 0.0, 1.0), vec3(-4.0, 1.0, 2.0), vec3(1.0, 2.0, 7.0)};
    const Matrix c = mesh.find_center_of_scaffolding_sphere(false);
    EXPECT_DOUBLE_EQ(c(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(c(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(c(2, 0), 7.0 - 15.0);
    EXPECT_DOUBLE_EQ(mesh.approximate_scaffolding_cap_radius(false), std::sqrt(26.0));

    // default requested: origin regardless of the scaffold
    EXPECT_EQ(norm3(mesh.find_center_of_scaffolding_sphere(true)), 0.0);
    EXPECT_DOUBLE_EQ(mesh.approximate_scaffolding_cap_radius(true), std::sqrt(17.0));
}

TEST_F(ScafCovGagTopology, LastReactionInTheFileIsKeptWithoutAnEndMarker)
{
    const std::string inp = tmp.file("last.inp");
    write_text(inp,
               "start reactions\n"
               "    A(c) + A(a) <-> A(c!1).A(a!1)\n"
               "        sigma = 1.75\n"
               "        assocAngles = [0.1, 0.2, 0.3, 0.4, 0.5]\n");
    Param p;
    configure_lattice(p, refFile, false);
    p.gagReactionFileName = inp;
    p.scaffoldingPoints = ring.points();
    Mesh mesh(p);
    ASSERT_TRUE(mesh.initialize_gag_scaffolding_topology());
    ASSERT_EQ(mesh.gagInteractions.size(), 3u);
    EXPECT_DOUBLE_EQ(mesh.gagInteractions[2].expectedSigma, 1.75); // c|a from the file
    EXPECT_DOUBLE_EQ(mesh.gagInteractions[2].assocAngles[4], 0.5);
    EXPECT_DOUBLE_EQ(mesh.gagBonds[2].expectedPhi, 0.35);
    // b|a is absent from this file: geometric targets
    EXPECT_NEAR(mesh.gagInteractions[0].expectedSigma, distance(ring.iface(0, "b"), ring.iface(1, "a")), 1e-12);
}

TEST_F(ScafCovGagTopology, FiniteDifferenceRoutinesInitialiseOnDemand)
{
    for (int which = 0; which < 3; ++which)
    {
        Param p;
        configure_lattice(p, refFile, true);
        p.scaffoldingPoints = ring.points();
        Mesh mesh(p);
        ASSERT_FALSE(mesh.gagScaffoldingTopologyInitialized);
        testing::internal::CaptureStdout();
        std::vector<Matrix> out;
        if (which == 0)
        {
            out = mesh.calculate_gag_scaffolding_forces_fd();
        }
        else if (which == 1)
        {
            out = mesh.calculate_gag_subunit_translation_forces_fd();
        }
        else
        {
            out = mesh.calculate_gag_subunit_rotation_torques_fd();
        }
        testing::internal::GetCapturedStdout();
        EXPECT_TRUE(mesh.gagScaffoldingTopologyInitialized);
        ASSERT_EQ(out.size(), 3u);
        // initialised at the current points, so the scaffold is at rest
        for (const Matrix &f : out)
        {
            EXPECT_LT(norm3(f), 1e-6);
        }
    }
}
