// Error-vs-cost study for extraordinary-vertex patches.
// Results: docs/irregular_patch_results.md section 9.
//
// Standalone: links build-omp/libslimed_core.a, uses the library's mesh setup
// (one-rings, patch kinds, multi-EV prolongations) and its kernels, but runs
// its own evaluation loop so every candidate scheme can be expressed as
//   for d < D, c < 3:  kernel(rows_rule(d)[N,d,c], w(d) * gauss_rule(d))
// which is an exact energy with exact (conservative) forces for any fixed w.
//
// Build (from the repository root, after building build-omp in Release):
//   c++ -DOMP -Iinclude -isystem /opt/homebrew/include \
//       -isystem /opt/homebrew/opt/libomp/include -O3 -DNDEBUG -std=c++14 \
//       -Xclang -fopenmp analysis/ev_patch_cost.cpp build-omp/libslimed_core.a \
//       -L/opt/homebrew/lib -L/opt/homebrew/opt/libomp/lib -lgsl -lgslcblas -lomp \
//       -o ev_patch_cost
// Modes (first argument):
//   rates     the per-valence eigenvalues and mode rates
//   fan       two-ring disks around one valence-N vertex, bowl/cyl/saddle/generic
//   reg6      the 3-point rule's own error on regular faces
//   lattice   a flipped hex lattice: single- and multi-EV faces together
//   latsweep  the lattice, one valence at a time, depth sweep
//   real P V F C [perval]
//             a saved run state: input.params, inputvertex_type_final.csv,
//             inputface_<last>.csv, inputmesh_faces.csv (copy flags);
//             "perval" sweeps one valence at a time with the rest at the reference
//   time      per-face wall time against sample count
//   adaptive P V F C
//             per-face depth from a planar-angle curvature proxy, feasibility
#include "mesh/Mesh.hpp"
#include "energy_force/Patch_kernel.hpp"
#include "energy_force/Patch_rows_flat.hpp"
#include "mesh/Gauss_quadrature.hpp"
#include "mesh/Irregular_patch_rows.hpp"
#include "mesh/Multi_extraordinary_patch.hpp"
#include "io/io.hpp"
#include <fstream>
#include <sstream>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace
{
constexpr int kTableDepth = 17;

struct Rule
{
    std::vector<Matrix> sf;
    Matrix weight;
    IrregularPatchRowTable table;
    PatchRowsFlat flat;
    std::vector<double> gauss;
};

Rule *make_rule(int gaussN)
{
    Rule *r = new Rule;
    Matrix vwu;
    get_gauss_quadrature_weight_VWU(gaussN, vwu, r->weight);
    get_shapefunction_vector(vwu, r->sf);
    r->table.build(r->sf, kTableDepth, DepthPolicy::Uniform);
    r->flat.build(r->sf, r->table, r->weight);
    for (int i = 0; i < r->flat.nSamples(); i++)
        r->gauss.push_back(r->flat.gaussCoeff()[i]);
    return r;
}

Rule *R3, *R1, *R12;

/// A candidate evaluation scheme for single-EV patches (and the corner
/// children of multi-EV faces).
struct Scheme
{
    std::string name;
    // depth for valence N (index by N)
    int depth[9] = {0, 0, 0, 0, 12, 8, 6, 7, 7};
    // weight multiplier on level d of a depth-D patch
    std::function<double(int N, int D, int d)> weight = [](int, int, int) { return 1.0; };
    int onePointFrom = 1 << 20; // levels d >= this use the centroid rule
    bool highOrder = false;     // 12-point rule everywhere (quadrature reference)
    int multiCornerReduce = 0;  // corner children of multi-EV faces get depth D - this
    // area/volume weights; empty = the bending weights
    std::function<double(int N, int D, int d)> areaWeight;
};

struct Eval
{
    double eBend = 0.0;
    double area = 0.0;
    long samples = 0;
};

const slimed::PatchParams &bending_params()
{
    static slimed::PatchParams p = [] {
        slimed::PatchParams q;
        q.kCurv = 1.0;
        return q;
    }();
    return p;
}

/// One patch of valence N over control net coords (nCtrl x 3), forces into fOut.
void eval_patch(const Scheme &s, int N, int D, const double *coords, int nCtrl, double *fOut,
                Eval &e)
{
    const slimed::PatchParams &p = bending_params();
    double fArea[slimed::kMaxControlPoints * 3] = {0};
    double fVol[slimed::kMaxControlPoints * 3] = {0};
    double vol = 0.0;
    if (N == 6)
    {
        Rule *r = s.highOrder ? R12 : R3;
        double eb = 0, mc = 0, nv[3] = {0, 0, 0};
        slimed::element_energy_force_patch_pod(r->flat.regular(), r->gauss.data(),
                                               r->flat.nSamples(), coords, nCtrl, p, eb, mc,
                                               nv, fOut, fArea, fVol);
        slimed::element_area_volume_pod(r->flat.regular(), r->gauss.data(), r->flat.nSamples(),
                                        coords, nCtrl, e.area, vol);
        e.eBend += eb;
        e.samples += r->flat.nSamples();
        return;
    }
    for (int d = 0; d < D; d++)
    {
        Rule *r = s.highOrder ? R12 : (d >= s.onePointFrom ? R1 : R3);
        const double w = s.weight(N, D, d);
        const double wa = s.areaWeight ? s.areaWeight(N, D, d) : w;
        double g[16], ga[16];
        for (int q = 0; q < r->flat.nSamples(); q++)
        {
            g[q] = w * r->gauss[q];
            ga[q] = wa * r->gauss[q];
        }
        for (int c = 0; c < kRegularChildrenPerStep; c++)
        {
            double eb = 0, mc = 0, nv[3] = {0, 0, 0};
            slimed::element_energy_force_patch_pod(r->flat.child(N, d, c), g, r->flat.nSamples(),
                                                   coords, nCtrl, p, eb, mc, nv, fOut, fArea,
                                                   fVol);
            slimed::element_area_volume_pod(r->flat.child(N, d, c), ga, r->flat.nSamples(),
                                            coords, nCtrl, e.area, vol);
            e.eBend += eb;
            e.samples += r->flat.nSamples();
        }
    }
}

/// One face of the mesh, any patch kind. Forces scattered into fGlobal (nV x 3).
void eval_face(const Scheme &s, Mesh &mesh, int iFace, std::vector<double> &fGlobal, Eval &e)
{
    Face &face = mesh.faces[iFace];
    const int n = static_cast<int>(face.oneRingVertices.size());
    double coords[slimed::kMaxControlPoints * 3];
    for (int j = 0; j < n; j++)
        for (int a = 0; a < 3; a++)
            coords[j * 3 + a] = mesh.vertices[face.oneRingVertices[j]].coord.get(a, 0);
    double f[slimed::kMaxControlPoints * 3] = {0};
    if (face.patchKind == PatchKind::Regular || face.patchKind == PatchKind::SingleExtraordinary)
    {
        const int N = face.patchKind == PatchKind::Regular ? 6 : n - 6;
        eval_patch(s, N, s.depth[N], coords, n, f, e);
    }
    else if (face.patchKind == PatchKind::MultiExtraordinary && face.patchEntry >= 0)
    {
        const MultiPatchTable::Entry &entry = mesh.multiPatchTable.entry(face.patchEntry);
        for (int c = 0; c < 4; c++)
        {
            const MultiPatchTable::Child &child = entry.children[c];
            const double *M = mesh.multiPatchTable.data() + child.offset;
            double cc[slimed::kMaxControlPoints * 3];
            multi_patch_prolong(M, coords, entry.nControl, child.nControl, cc);
            double cf[slimed::kMaxControlPoints * 3] = {0};
            int D = s.depth[child.valence];
            if (child.valence != 6 && c != 3)
                D = D == 0 ? 0 : std::max(1, D - s.multiCornerReduce);
            eval_patch(s, child.valence, D, cc, child.nControl, cf, e);
            multi_patch_scatter(M, cf, entry.nControl, child.nControl, f);
        }
    }
    else
        return;
    for (int j = 0; j < n; j++)
        for (int a = 0; a < 3; a++)
            fGlobal[face.oneRingVertices[j] * 3 + a] += f[j * 3 + a];
}

struct Total
{
    double eBend = 0, area = 0;
    long samples = 0;
    int faces = 0;
    std::vector<double> force;
};

Total eval_faces(const Scheme &s, Mesh &mesh, const std::vector<int> &faceIds)
{
    Total t;
    t.force.assign(mesh.vertices.size() * 3, 0.0);
    for (int i : faceIds)
    {
        Eval e;
        eval_face(s, mesh, i, t.force, e);
        t.eBend += e.eBend;
        t.area += e.area;
        t.samples += e.samples;
        t.faces++;
    }
    return t;
}

double force_err(const Total &a, const Total &ref)
{
    double num = 0, den = 0;
    for (std::size_t i = 0; i < ref.force.size(); i += 3)
    {
        double dn = 0, rn = 0;
        for (int k = 0; k < 3; k++)
        {
            dn += std::pow(a.force[i + k] - ref.force[i + k], 2);
            rn += std::pow(ref.force[i + k], 2);
        }
        num = std::max(num, std::sqrt(dn));
        den = std::max(den, std::sqrt(rn));
    }
    return den > 0 ? num / den : 0.0;
}

// ---------------------------------------------------------------- geometry
using Shape = std::function<double(double, double)>;

std::map<std::string, Shape> shapes()
{
    const double a = 0.12;
    return {
        {"bowl", [a](double x, double y) { return a * (x * x + y * y); }},
        {"cyl", [a](double x, double y) { return a * x * x; }},
        {"saddle", [a](double x, double y) { return a * (x * x - y * y); }},
        // a rotated cylinder plus a cubic: every angular mode, no symmetry
        {"generic",
         [a](double x, double y) {
             const double c = std::cos(0.37), s = std::sin(0.37);
             const double u = c * x + s * y;
             return a * u * u + 0.03 * (x * x * x - 3 * x * y * y) + 0.02 * x * y * y;
         }},
    };
}

Param make_param()
{
    Param param;
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.uVol = 0.0;
    return param;
}

/// Two-ring disk around one valence-N vertex: spokes u_k, outer ring at 2u_k and u_k + u_{k+1}.
void build_fan(int N, const Shape &z, std::vector<std::vector<double>> &V,
               std::vector<std::vector<int>> &F)
{
    V.clear();
    F.clear();
    V.push_back({0, 0, z(0, 0)});
    std::vector<double> ux(N), uy(N);
    for (int k = 0; k < N; k++)
    {
        ux[k] = std::cos(2 * M_PI * k / N);
        uy[k] = std::sin(2 * M_PI * k / N);
        V.push_back({ux[k], uy[k], z(ux[k], uy[k])});
    }
    for (int k = 0; k < N; k++) // 2u_k at 1+N+k
        V.push_back({2 * ux[k], 2 * uy[k], z(2 * ux[k], 2 * uy[k])});
    for (int k = 0; k < N; k++) // m_k = u_k + u_{k+1} at 1+2N+k
    {
        const int k1 = (k + 1) % N;
        const double x = ux[k] + ux[k1], y = uy[k] + uy[k1];
        V.push_back({x, y, z(x, y)});
    }
    for (int k = 0; k < N; k++)
    {
        const int k1 = (k + 1) % N;
        const int a = 1 + k, b = 1 + k1, ta = 1 + N + k, tb = 1 + N + k1, m = 1 + 2 * N + k;
        F.push_back({0, a, b});
        F.push_back({a, ta, m});
        F.push_back({a, m, b});
        F.push_back({b, m, tb});
    }
}

/// Hex lattice disk with chosen edge flips. Returns vertex positions in-plane.
void build_lattice(int R, const Shape &z,
                   std::vector<std::vector<double>> &V, std::vector<std::vector<int>> &F,
                   int flipSet)
{
    std::map<std::pair<int, int>, int> id;
    V.clear();
    F.clear();
    for (int j = -R; j <= R; j++)
        for (int i = -R; i <= R; i++)
        {
            if (std::abs(i + j) > R)
                continue;
            const double x = i + 0.5 * j, y = j * std::sqrt(3.0) / 2;
            id[{i, j}] = static_cast<int>(V.size());
            V.push_back({x, y, 0.0});
        }
    auto has = [&](int i, int j) { return id.count({i, j}) > 0; };
    for (int j = -R; j <= R; j++)
        for (int i = -R; i <= R; i++)
        {
            if (has(i, j) && has(i + 1, j) && has(i, j + 1))
                F.push_back({id[{i, j}], id[{i + 1, j}], id[{i, j + 1}]});
            if (has(i + 1, j) && has(i + 1, j + 1) && has(i, j + 1))
                F.push_back({id[{i + 1, j}], id[{i + 1, j + 1}], id[{i, j + 1}]});
        }
    // flip edge (A,B), given as lattice-vertex pairs
    auto flip = [&](int A, int B) {
        int f1 = -1, f2 = -1, e1 = -1, e2 = -1;
        for (int f = 0; f < (int)F.size(); f++)
            for (int e = 0; e < 3; e++)
            {
                if (F[f][e] == A && F[f][(e + 1) % 3] == B)
                    f1 = f, e1 = e;
                if (F[f][e] == B && F[f][(e + 1) % 3] == A)
                    f2 = f, e2 = e;
            }
        if (f1 < 0 || f2 < 0)
        {
            std::fprintf(stderr, "flip failed\n");
            std::exit(1);
        }
        const int C = F[f1][(e1 + 2) % 3], D = F[f2][(e2 + 2) % 3];
        F[f1] = {A, D, C};
        F[f2] = {D, B, C};
    };
    auto vid = [&](int i, int j) { return id.at({i, j}); };
    if (flipSet >= 1) // one flip: 5,5,7,7
        flip(vid(0, 0), vid(1, 0));
    if (flipSet >= 2) // a second, elsewhere: 4..8 combinations
    {
        // around (-2,2): lower it twice, raise a common neighbour twice
        flip(vid(-2, 2), vid(-1, 2));
        flip(vid(-2, 2), vid(-2, 3));
        // and a lone 5-7 pair further out
        flip(vid(2, -3), vid(2, -2));
    }
    for (auto &v : V)
        v[2] = z(v[0], v[1]);
}

double now()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch())
        .count();
}

// ---------------------------------------------------------------- schemes
const double kPi = 3.14159265358979323846;
double lam(int N) { return 3.0 / 8 + std::cos(2 * kPi / N) / 4; }
double mu2(int N) { return 3.0 / 8 + std::cos(4 * kPi / N) / 4; }
/// per-level bending ratio of the bowl (k=0) mode: 1/(16 lambda^2) per level
double r_bowl(int N) { return 1.0 / (16 * lam(N) * lam(N)); }
/// saddle (k=2) mode: mu2^2 / lambda^2
double r_sad(int N) { return mu2(N) * mu2(N) / (lam(N) * lam(N)); }
double r_slow(int N) { return std::max(r_bowl(N), r_sad(N)); }

/// Last-level weight 1/(1-r).
std::function<double(int, int, int)> tail1(std::function<double(int)> r)
{
    return [r](int N, int D, int d) { return d == D - 1 ? 1.0 / (1.0 - r(N)) : 1.0; };
}
/// Two-level weights summing both the bowl and saddle tails exactly.
double tail2_weight(int N, int D, int d)
{
    const double r1 = r_bowl(N), r2 = r_sad(N);
    const double g1 = 1 / (1 - r1), g2 = 1 / (1 - r2);
    if (D < 2)
        return d == D - 1 ? 1 / (1 - r_slow(N)) : 1.0;
    const double wb = (g1 - g2) / (r1 - r2);
    const double wa = g1 - wb * r1;
    if (d == D - 1)
        return wb;
    if (d == D - 2)
        return wa;
    return 1.0;
}

Scheme with_depth(Scheme s, int N, int D)
{
    s.depth[N] = D;
    return s;
}

// ---------------------------------------------------------------- adaptive
/// The saved-state loader of the `real` mode: params, vertex_type csv, face
/// csv and (optionally) the copy flags. `param` must outlive the mesh, which
/// keeps a reference to it.
Mesh *load_real_state(Param &param, int argc, char **argv)
{
    import_param_file(param, argv[2]);
    param.VERBOSE_MODE = false;
    std::vector<std::vector<double>> V;
    std::vector<VertexType> types;
    std::vector<int> mirrors;
    {
        std::ifstream in(argv[3]);
        std::string line;
        while (std::getline(in, line))
        {
            std::stringstream ss(line);
            std::string tok;
            std::vector<std::string> t;
            while (std::getline(ss, tok, ','))
                t.push_back(tok);
            if (t.size() < 5)
                continue;
            V.push_back({std::stod(t[0]), std::stod(t[1]), std::stod(t[2])});
            types.push_back(static_cast<VertexType>(std::stoi(t[3])));
            mirrors.push_back(std::stoi(t[4]));
        }
    }
    std::vector<std::vector<int>> F;
    {
        std::ifstream in(argv[4]);
        std::string line;
        while (std::getline(in, line))
        {
            int a, b, c;
            if (std::sscanf(line.c_str(), "%d,%d,%d", &a, &b, &c) == 3)
                F.push_back({a, b, c});
        }
    }
    std::vector<char> isCopy;
    if (argc > 5)
    {
        std::ifstream in(argv[5]);
        std::string line;
        while (std::getline(in, line))
        {
            int a, b, c, cp;
            if (std::sscanf(line.c_str(), "%d,%d,%d,%d", &a, &b, &c, &cp) == 4)
                isCopy.push_back(static_cast<char>(cp));
        }
    }
    Mesh *mesh = new Mesh(param);
    mesh->setup_from_vertices_faces(V, F, types, mirrors, isCopy);
    mesh->calculate_element_area_volume();
    mesh->sum_membrane_area_and_volume(mesh->param.area0, mesh->param.vol0);
    mesh->Compute_Energy_And_Force();
    return mesh;
}

/// eval_face() without the scatter: the face's own force on its control net.
void eval_face_local(const Scheme &s, Mesh &mesh, int iFace, std::vector<double> &fLocal, Eval &e)
{
    Face &face = mesh.faces[iFace];
    const int n = static_cast<int>(face.oneRingVertices.size());
    double coords[slimed::kMaxControlPoints * 3];
    for (int j = 0; j < n; j++)
        for (int a = 0; a < 3; a++)
            coords[j * 3 + a] = mesh.vertices[face.oneRingVertices[j]].coord.get(a, 0);
    double f[slimed::kMaxControlPoints * 3] = {0};
    if (face.patchKind == PatchKind::Regular || face.patchKind == PatchKind::SingleExtraordinary)
    {
        const int N = face.patchKind == PatchKind::Regular ? 6 : n - 6;
        eval_patch(s, N, s.depth[N], coords, n, f, e);
    }
    else if (face.patchKind == PatchKind::MultiExtraordinary && face.patchEntry >= 0)
    {
        const MultiPatchTable::Entry &entry = mesh.multiPatchTable.entry(face.patchEntry);
        for (int c = 0; c < 4; c++)
        {
            const MultiPatchTable::Child &child = entry.children[c];
            const double *M = mesh.multiPatchTable.data() + child.offset;
            double cc[slimed::kMaxControlPoints * 3];
            multi_patch_prolong(M, coords, entry.nControl, child.nControl, cc);
            double cf[slimed::kMaxControlPoints * 3] = {0};
            int D = s.depth[child.valence];
            if (child.valence != 6 && c != 3)
                D = D == 0 ? 0 : std::max(1, D - s.multiCornerReduce);
            eval_patch(s, child.valence, D, cc, child.nControl, cf, e);
            multi_patch_scatter(M, cf, entry.nControl, child.nControl, f);
        }
    }
    fLocal.assign(f, f + n * 3);
}

/// Unit normal of the flat triangle of a face, from its three corners.
void flat_face_normal(const Mesh &mesh, const Face &face, double n[3])
{
    double p[3][3];
    for (int k = 0; k < 3; k++)
        for (int a = 0; a < 3; a++)
            p[k][a] = mesh.vertices[face.adjacentVertices[k]].coord.get(a, 0);
    double u[3], v[3];
    for (int a = 0; a < 3; a++)
    {
        u[a] = p[1][a] - p[0][a];
        v[a] = p[2][a] - p[0][a];
    }
    n[0] = u[1] * v[2] - u[2] * v[1];
    n[1] = u[2] * v[0] - u[0] * v[2];
    n[2] = u[0] * v[1] - u[1] * v[0];
    const double len = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for (int a = 0; a < 3; a++)
        n[a] /= len;
}

double angle_between(const double a[3], const double b[3])
{
    double d = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    d = std::max(-1.0, std::min(1.0, d));
    return std::acos(d);
}

/// Planar-angle curvature proxies of one face, from flat triangle normals only:
///   edge  : largest dihedral angle across the face's three edges
///   ring  : largest angle between an incident face normal and the mean normal
///           at any of the face's three corners (the normal cone of the star)
void face_proxies(const Mesh &mesh, const std::vector<std::vector<double>> &normals, int iFace,
                  double &edgeProxy, double &ringProxy)
{
    const Face &face = mesh.faces[iFace];
    edgeProxy = 0.0;
    for (int j : face.adjacentFaces)
        if (j >= 0 && j < (int)normals.size() && j != iFace)
            edgeProxy = std::max(edgeProxy, angle_between(normals[iFace].data(), normals[j].data()));
    ringProxy = 0.0;
    for (int k = 0; k < 3; k++)
    {
        const Vertex &v = mesh.vertices[face.adjacentVertices[k]];
        double m[3] = {0, 0, 0};
        for (int j : v.adjacentFaces)
            if (j >= 0 && j < (int)normals.size())
                for (int a = 0; a < 3; a++)
                    m[a] += normals[j][a];
        const double len = std::sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
        if (len == 0.0)
            continue;
        for (int a = 0; a < 3; a++)
            m[a] /= len;
        for (int j : v.adjacentFaces)
            if (j >= 0 && j < (int)normals.size())
                ringProxy = std::max(ringProxy, angle_between(normals[j].data(), m));
    }
}

/// Spearman rank correlation.
double spearman(std::vector<double> x, std::vector<double> y)
{
    auto ranks = [](const std::vector<double> &v) {
        std::vector<int> idx(v.size());
        for (std::size_t i = 0; i < v.size(); i++)
            idx[i] = (int)i;
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return v[a] < v[b]; });
        std::vector<double> r(v.size());
        for (std::size_t i = 0; i < idx.size(); i++)
            r[idx[i]] = (double)i;
        return r;
    };
    const std::vector<double> rx = ranks(x), ry = ranks(y);
    const double n = (double)rx.size();
    double mx = 0, my = 0;
    for (std::size_t i = 0; i < rx.size(); i++)
    {
        mx += rx[i];
        my += ry[i];
    }
    mx /= n;
    my /= n;
    double sxy = 0, sxx = 0, syy = 0;
    for (std::size_t i = 0; i < rx.size(); i++)
    {
        sxy += (rx[i] - mx) * (ry[i] - my);
        sxx += (rx[i] - mx) * (rx[i] - mx);
        syy += (ry[i] - my) * (ry[i] - my);
    }
    return sxy / std::sqrt(sxx * syy);
}
} // namespace

int main(int argc, char **argv)
{
    R3 = make_rule(2);
    R1 = make_rule(1);
    R12 = make_rule(6);
    const std::string mode = argc > 1 ? argv[1] : "fan";

    std::vector<Scheme> families;
    {
        Scheme plain;
        plain.name = "plain";
        families.push_back(plain);
        Scheme b1 = plain;
        b1.name = "tail_bowl";
        b1.weight = tail1(r_bowl);
        families.push_back(b1);
        Scheme b2 = plain;
        b2.name = "tail_slow";
        b2.weight = tail1(r_slow);
        families.push_back(b2);
        Scheme b3 = plain;
        b3.name = "tail2";
        b3.weight = tail2_weight;
        families.push_back(b3);
        Scheme b3a = b3;
        b3a.name = "tail2A";
        b3a.areaWeight = [](int N, int D, int d) {
            return d == D - 1 ? 1.0 / (1.0 - lam(N) * lam(N)) : 1.0;
        };
        families.push_back(b3a);
        for (int from : {1, 2, 3})
        {
            Scheme c = plain;
            c.name = "c1pt_from" + std::to_string(from);
            c.onePointFrom = from;
            families.push_back(c);
            Scheme cb = b3;
            cb.name = "tail2_c1pt_from" + std::to_string(from);
            cb.onePointFrom = from;
            families.push_back(cb);
        }
        Scheme hi = plain;
        hi.name = "rule12";
        hi.highOrder = true;
        families.push_back(hi);
    }

    if (mode == "reg6")
    {
        // quadrature floor on regular faces: the same fan at valence 6, 3pt vs 12pt
        for (auto &sh : shapes())
        {
            std::vector<std::vector<double>> V;
            std::vector<std::vector<int>> F;
            build_fan(6, sh.second, V, F);
            Param param = make_param();
            Mesh mesh(param);
            mesh.setup_from_vertices_faces(V, F);
            std::vector<int> fan;
            for (int i = 0; i < (int)mesh.faces.size(); i++)
                if (mesh.faces[i].patchKind == PatchKind::Regular)
                    fan.push_back(i);
            Scheme p3, p12;
            p12.highOrder = true;
            const Total a = eval_faces(p3, mesh, fan), b = eval_faces(p12, mesh, fan);
            std::printf("%s\tN=6 faces %d\t3pt-vs-12pt relE %.3e relA %.3e relF %.3e\n",
                        sh.first.c_str(), a.faces, (a.eBend - b.eBend) / b.eBend,
                        (a.area - b.area) / b.area, force_err(a, b));
        }
        return 0;
    }

    if (mode == "rates")
    {
        for (int N : {4, 5, 7, 8})
            std::printf("N=%d lambda=%.6f mu2=%.6f bowl 1/r=%.4f saddle 1/r=%.4f\n", N, lam(N),
                        mu2(N), 1 / r_bowl(N), 1 / r_sad(N));
        return 0;
    }

    if (mode == "fan")
    {
        std::printf("shape\tN\tscheme\tD\tsamples_per_face\tE\tA\trelE\trelA\trelF\n");
        for (auto &sh : shapes())
            for (int N : {4, 5, 7, 8})
            {
                std::vector<std::vector<double>> V;
                std::vector<std::vector<int>> F;
                build_fan(N, sh.second, V, F);
                Param param = make_param();
                Mesh mesh(param);
                mesh.setup_from_vertices_faces(V, F);
                std::vector<int> fan;
                for (int i = 0; i < (int)mesh.faces.size(); i++)
                    if (mesh.faces[i].patchKind == PatchKind::SingleExtraordinary)
                        fan.push_back(i);
                if ((int)fan.size() != N)
                {
                    std::fprintf(stderr, "N=%d: %zu EV faces\n", N, fan.size());
                    return 1;
                }
                Scheme refS;
                refS.name = "ref";
                refS.weight = tail2_weight;
                refS.areaWeight = [](int N, int D, int d) {
                    return d == D - 1 ? 1.0 / (1.0 - lam(N) * lam(N)) : 1.0;
                };
                const Total ref = eval_faces(with_depth(refS, N, 16), mesh, fan);
                std::printf("%s\t%d\tref16x\t16\t%.1f\t%.12e\t%.12e\t0\t0\t0\n", sh.first.c_str(),
                            N, (double)ref.samples / N, ref.eBend, ref.area);
                for (const Scheme &fam : families)
                    for (int D = 1; D <= 14; D++)
                    {
                        if (fam.highOrder && D != 14)
                            continue;
                        const Total t = eval_faces(with_depth(fam, N, D), mesh, fan);
                        std::printf("%s\t%d\t%s\t%d\t%.1f\t%.12e\t%.12e\t%.3e\t%.3e\t%.3e\n",
                                    sh.first.c_str(), N, fam.name.c_str(), D,
                                    (double)t.samples / N, t.eBend, t.area,
                                    (t.eBend - ref.eBend) / ref.eBend,
                                    (t.area - ref.area) / ref.area, force_err(t, ref));
                    }
            }
        return 0;
    }

    if (mode == "lattice")
    {
        // Mixed-valence lattice: whole-sheet totals, single- and multi-EV
        // faces apart. Depth sets are per valence.
        struct DepthSet
        {
            const char *name;
            int d[9];
        };
        const DepthSet sets[] = {
            {"prod", {0, 0, 0, 0, 12, 8, 6, 7, 7}},
            {"x0.75", {0, 0, 0, 0, 9, 6, 6, 5, 5}},
            {"x0.5", {0, 0, 0, 0, 6, 4, 6, 4, 4}},
            {"6444", {0, 0, 0, 0, 6, 4, 6, 4, 4}},
            {"5333", {0, 0, 0, 0, 5, 3, 6, 3, 3}},
        };
        std::printf("shape\tflips\tscheme\tdepths\tcornerReduce\tkind\tfaces\tsamples\tE\tA\trelE\trelA\trelF\n");
        for (auto &sh : shapes())
            for (int flipSet : {1, 2})
            {
                std::vector<std::vector<double>> V;
                std::vector<std::vector<int>> F;
                build_lattice(7, sh.second, V, F, flipSet);
                Param param = make_param();
                Mesh mesh(param);
                mesh.setup_from_vertices_faces(V, F);
                std::vector<int> single, multi, all;
                for (int i = 0; i < (int)mesh.faces.size(); i++)
                {
                    const PatchKind k = mesh.faces[i].patchKind;
                    if (k == PatchKind::SingleExtraordinary)
                        single.push_back(i);
                    if (k == PatchKind::MultiExtraordinary)
                        multi.push_back(i);
                    if (k == PatchKind::SingleExtraordinary || k == PatchKind::MultiExtraordinary)
                        all.push_back(i);
                }
                Scheme refS;
                for (int N = 4; N <= 8; N++)
                    if (N != 6)
                        refS.depth[N] = 16;
                refS.weight = tail2_weight;
                refS.areaWeight = [](int N, int D, int d) {
                    return d == D - 1 ? 1.0 / (1.0 - lam(N) * lam(N)) : 1.0;
                };
                const Total refAll = eval_faces(refS, mesh, all);
                const Total refMulti = eval_faces(refS, mesh, multi);
                const Total refSingle = eval_faces(refS, mesh, single);
                for (const Scheme &fam0 : families)
                {
                    if (fam0.name != "plain" && fam0.name != "tail2A" && fam0.name != "c1pt_from3")
                        continue;
                    for (const DepthSet &ds : sets)
                        for (int reduce : {0, 1})
                        {
                            Scheme s = fam0;
                            for (int N = 4; N <= 8; N++)
                                s.depth[N] = ds.d[N];
                            s.multiCornerReduce = reduce;
                            const Total ta = eval_faces(s, mesh, all);
                            const Total tm = eval_faces(s, mesh, multi);
                            const Total ts = eval_faces(s, mesh, single);
                            const Total *T[3] = {&ta, &tm, &ts};
                            const Total *Rf[3] = {&refAll, &refMulti, &refSingle};
                            const char *kind[3] = {"all", "multi", "single"};
                            for (int k = 0; k < 3; k++)
                                std::printf("%s\t%d\t%s\t%s\t%d\t%s\t%d\t%ld\t%.12e\t%.12e\t%.3e\t%.3e\t%.3e\n",
                                            sh.first.c_str(), flipSet, fam0.name.c_str(), ds.name,
                                            reduce, kind[k], T[k]->faces, T[k]->samples,
                                            T[k]->eBend, T[k]->area,
                                            (T[k]->eBend - Rf[k]->eBend) / refAll.eBend,
                                            (T[k]->area - Rf[k]->area) / refAll.area,
                                            force_err(*T[k], *Rf[k]));
                        }
                }
            }
        return 0;
    }

    if (mode == "adaptive")
    {
        // Feasibility of a per-face depth chosen from a planar-angle curvature
        // estimate (docs/irregular_patch_results.md section 10.7). Every face is
        // evaluated under each candidate depth set and under the depth-16
        // reference, so any assignment rule can be scored by summing stored
        // per-face results without re-evaluating anything.
        Param param;
        Mesh *meshPtr = load_real_state(param, argc, argv);
        Mesh &mesh = *meshPtr;
        const double kT = param.KBT;
        const double kc = param.kCurv;

        std::vector<int> all;
        for (int i = 0; i < (int)mesh.faces.size(); i++)
        {
            const Face &f = mesh.faces[i];
            if (f.isBoundary)
                continue;
            const bool ok = f.patchKind == PatchKind::Regular ||
                            f.patchKind == PatchKind::SingleExtraordinary ||
                            (f.patchKind == PatchKind::MultiExtraordinary && f.patchEntry >= 0);
            if (ok)
                all.push_back(i);
        }
        const int nF = (int)all.size();

        // Proxies, and their cost.
        std::vector<std::vector<double>> normals(mesh.faces.size(), std::vector<double>(3, 0.0));
        std::vector<double> edgeProxy(nF, 0.0), ringProxy(nF, 0.0);
        const auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < (int)mesh.faces.size(); i++)
            if (!mesh.faces[i].isBoundary)
                flat_face_normal(mesh, mesh.faces[i], normals[i].data());
        for (int k = 0; k < nF; k++)
            face_proxies(mesh, normals, all[k], edgeProxy[k], ringProxy[k]);
        const double proxyMicros =
            std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();

        // Tiers.
        struct Tier
        {
            const char *name;
            int d[9];
        };
        const std::vector<Tier> tiers = {
            {"3333", {0, 0, 0, 0, 3, 3, 6, 3, 3}}, {"4444", {0, 0, 0, 0, 4, 4, 6, 4, 4}},
            {"5444", {0, 0, 0, 0, 5, 4, 6, 4, 4}}, {"6555", {0, 0, 0, 0, 6, 5, 6, 5, 5}},
            {"8666", {0, 0, 0, 0, 8, 6, 6, 6, 6}},
        };
        const int T = (int)tiers.size();
        const int iT3 = 0, iT4 = 1, iT6 = 3, iT8 = 4;
        Scheme tail2A;
        tail2A.name = "tail2A";
        tail2A.weight = tail2_weight;
        tail2A.areaWeight = [](int N, int D, int d) {
            return d == D - 1 ? 1.0 / (1.0 - lam(N) * lam(N)) : 1.0;
        };
        Scheme refS = tail2A;
        for (int N = 4; N <= 8; N++)
            if (N != 6)
                refS.depth[N] = 16;

        // Per-face, per-tier results.
        std::vector<std::vector<double>> E(T, std::vector<double>(nF)), A(T, std::vector<double>(nF));
        std::vector<std::vector<long>> S(T, std::vector<long>(nF));
        std::vector<std::vector<std::vector<double>>> Fl(T, std::vector<std::vector<double>>(nF));
        std::vector<double> Eref(nF), Aref(nF);
        std::vector<std::vector<double>> Fref(nF);
        std::vector<char> isEv(nF, 0);
        for (int k = 0; k < nF; k++)
        {
            isEv[k] = mesh.faces[all[k]].patchKind != PatchKind::Regular;
            Eval e;
            eval_face_local(refS, mesh, all[k], Fref[k], e);
            Eref[k] = e.eBend;
            Aref[k] = e.area;
            for (int t = 0; t < T; t++)
            {
                Scheme s = tail2A;
                for (int N = 4; N <= 8; N++)
                    s.depth[N] = tiers[t].d[N];
                Eval et;
                eval_face_local(s, mesh, all[k], Fl[t][k], et);
                E[t][k] = et.eBend;
                A[t][k] = et.area;
                S[t][k] = et.samples;
            }
        }
        auto faceForceErr = [&](int t, int k) {
            double d = 0;
            for (std::size_t j = 0; j < Fref[k].size(); j++)
                d += std::pow(Fl[t][k][j] - Fref[k][j], 2);
            return std::sqrt(d);
        };
        double sheetE = 0;
        for (int k = 0; k < nF; k++)
            sheetE += Eref[k];
        // The oracle: rank faces by their own depth-3 force error. No cheap
        // estimate can do better than this, so it bounds the adaptive gain.
        std::vector<double> oracleProxy(nF, 0.0);
        for (int k = 0; k < nF; k++)
            oracleProxy[k] = faceForceErr(iT3, k);
        auto proxyByIndex = [&](int which) -> const std::vector<double> & {
            return which == 0 ? edgeProxy : which == 1 ? ringProxy : oracleProxy;
        };
        const char *proxyNames[3] = {"edge", "ring", "oracle_dF3"};

        // A. Does the proxy predict the per-face error?
        std::printf("# proxy cost: %.1f us for %d faces (%.3f us/face)\n", proxyMicros, nF,
                    proxyMicros / nF);
        std::printf("# sheet E_bend (ref) %.6e kc-units = %.1f kT; EV faces %d of %d\n", sheetE,
                    sheetE * kc / kT, (int)std::count(isEv.begin(), isEv.end(), 1), nF);
        for (int which = 0; which < 2; which++)
        {
            const std::vector<double> &px = which == 0 ? edgeProxy : ringProxy;
            std::vector<double> pxEv, dE3, dE4, dF3, dF4, eRef;
            for (int k = 0; k < nF; k++)
                if (isEv[k])
                {
                    pxEv.push_back(px[k]);
                    dE3.push_back(std::fabs(E[iT3][k] - Eref[k]));
                    dE4.push_back(std::fabs(E[iT4][k] - Eref[k]));
                    dF3.push_back(faceForceErr(iT3, k));
                    dF4.push_back(faceForceErr(iT4, k));
                    eRef.push_back(Eref[k]);
                }
            std::printf("# proxy %s: Spearman(proxy, |dE| at 3333) %.3f, (proxy, |dF| at 3333) %.3f, "
                        "(proxy, E_face) %.3f; (E_face, |dE| at 3333) %.3f\n",
                        which == 0 ? "edge" : "ring", spearman(pxEv, dE3), spearman(pxEv, dF3),
                        spearman(pxEv, eRef), spearman(eRef, dE3));
            // quintiles
            std::vector<int> idx(pxEv.size());
            for (std::size_t i = 0; i < idx.size(); i++)
                idx[i] = (int)i;
            std::sort(idx.begin(), idx.end(), [&](int a, int b) { return pxEv[a] < pxEv[b]; });
            std::printf("proxy\tquintile\tproxy_min_deg\tproxy_max_deg\tmean_Eface_kT\tmean_|dE3|_kT\tmax_|dE3|_kT\tmean_|dE4|_kT\tmean_dF3_pN\tmax_dF3_pN\tmean_dF4_pN\n");
            const int nq = 5;
            for (int q = 0; q < nq; q++)
            {
                const int lo = (int)(idx.size() * q / nq), hi = (int)(idx.size() * (q + 1) / nq);
                double mE = 0, m3 = 0, x3 = 0, m4 = 0, f3 = 0, xf3 = 0, f4 = 0;
                for (int i = lo; i < hi; i++)
                {
                    const int k = idx[i];
                    mE += eRef[k];
                    m3 += dE3[k];
                    x3 = std::max(x3, dE3[k]);
                    m4 += dE4[k];
                    f3 += dF3[k];
                    xf3 = std::max(xf3, dF3[k]);
                    f4 += dF4[k];
                }
                const double n = hi - lo;
                const double toKT = kc / kT;
                std::printf("%s\t%d\t%.2f\t%.2f\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",
                            which == 0 ? "edge" : "ring", q + 1, pxEv[idx[lo]] * 180 / kPi,
                            pxEv[idx[hi - 1]] * 180 / kPi, mE / n * toKT, m3 / n * toKT, x3 * toKT,
                            m4 / n * toKT, f3 / n * kc, xf3 * kc, f4 / n * kc);
            }
        }

        // B. Threshold sweep: tier by proxy quantile among EV faces.
        // Step time from section 10.6: t = a + b * samples (6555: 70.6 ms at
        // 79 362 samples, 4333: 49.8 at 49 230, 3333: 49.1 at 48 438).
        const double kMsIntercept = 15.8, kMsPerSample = 6.9e-4;
        auto score = [&](const std::vector<int> &tierOf, const char *label, const char *proxyName,
                         double qLow, double qHigh, double qTop) {
            double e = 0, a = 0;
            long samples = 0;
            int count[8] = {0};
            std::vector<double> force(mesh.vertices.size() * 3, 0.0), fref(mesh.vertices.size() * 3, 0.0);
            double maxJump = 0;
            for (int k = 0; k < nF; k++)
            {
                const int t = tierOf[k];
                count[t]++;
                e += E[t][k];
                a += A[t][k];
                samples += S[t][k];
                const Face &face = mesh.faces[all[k]];
                for (std::size_t j = 0; j < Fref[k].size(); j++)
                {
                    force[face.oneRingVertices[j / 3] * 3 + j % 3] += Fl[t][k][j];
                    fref[face.oneRingVertices[j / 3] * 3 + j % 3] += Fref[k][j];
                }
                if (isEv[k] && t != iT6)
                    maxJump = std::max(maxJump, std::fabs(E[t][k] - E[iT6][k]));
            }
            double eRefAll = 0, aRefAll = 0;
            for (int k = 0; k < nF; k++)
            {
                eRefAll += Eref[k];
                aRefAll += Aref[k];
            }
            double num = 0, den = 0;
            for (std::size_t i = 0; i < fref.size(); i += 3)
            {
                double dn = 0, rn = 0;
                for (int c = 0; c < 3; c++)
                {
                    dn += std::pow(force[i + c] - fref[i + c], 2);
                    rn += std::pow(fref[i + c], 2);
                }
                num = std::max(num, std::sqrt(dn));
                den = std::max(den, std::sqrt(rn));
            }
            std::printf("%s\t%s\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%d\t%ld\t%.3e\t%.3e\t%.3e\t%.3e\t%.1f\n",
                        label, proxyName, qLow, qHigh, qTop, count[0], count[1], count[2], count[3],
                        count[4], samples, (e - eRefAll) / eRefAll, (a - aRefAll) / aRefAll,
                        num / den, maxJump * kc / kT, kMsIntercept + kMsPerSample * samples);
        };
        std::printf("rule\tproxy\tq_low\tq_high\tq_top\tn3333\tn4444\tn5444\tn6555\tn8666\tsamples\trelE\trelA\trelF\tmax_jump_kT\tpred_ms\n");
        for (int t = 0; t < T; t++)
        {
            std::vector<int> tierOf(nF, t);
            score(tierOf, "uniform", tiers[t].name, 0, 0, 0);
        }
        for (int which = 0; which < 3; which++)
        {
            const std::vector<double> &px = proxyByIndex(which);
            std::vector<double> sorted;
            for (int k = 0; k < nF; k++)
                if (isEv[k])
                    sorted.push_back(px[k]);
            std::sort(sorted.begin(), sorted.end());
            auto quant = [&](double q) {
                if (q <= 0)
                    return -1.0;
                if (q >= 1)
                    return 1e9;
                return sorted[std::min(sorted.size() - 1, (std::size_t)(q * sorted.size()))];
            };
            const double qs[] = {0.0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
            for (double qLow : qs)
                for (double qHigh : qs)
                {
                    if (qHigh < qLow)
                        continue;
                    for (double qTop : {1.0, 0.9})
                    {
                        if (qTop < qHigh)
                            continue;
                        std::vector<int> tierOf(nF, iT6);
                        for (int k = 0; k < nF; k++)
                        {
                            if (!isEv[k])
                                continue;
                            if (px[k] <= quant(qLow))
                                tierOf[k] = iT3;
                            else if (px[k] <= quant(qHigh))
                                tierOf[k] = iT4;
                            else if (px[k] > quant(qTop))
                                tierOf[k] = iT8;
                        }
                        score(tierOf, "adaptive", proxyNames[which], qLow, qHigh, qTop);
                    }
                }
        }
        return 0;
    }

    if (mode == "real")
    {
        // A real disordered state: argv[2] params, argv[3] vertex_type csv, argv[4] face csv
        Param param;
        import_param_file(param, argv[2]);
        param.VERBOSE_MODE = false;
        std::vector<std::vector<double>> V;
        std::vector<VertexType> types;
        std::vector<int> mirrors;
        {
            std::ifstream in(argv[3]);
            std::string line;
            while (std::getline(in, line))
            {
                std::stringstream ss(line);
                std::string tok;
                std::vector<std::string> t;
                while (std::getline(ss, tok, ','))
                    t.push_back(tok);
                if (t.size() < 5)
                    continue;
                V.push_back({std::stod(t[0]), std::stod(t[1]), std::stod(t[2])});
                types.push_back(static_cast<VertexType>(std::stoi(t[3])));
                mirrors.push_back(std::stoi(t[4]));
            }
        }
        std::vector<std::vector<int>> F;
        {
            std::ifstream in(argv[4]);
            std::string line;
            while (std::getline(in, line))
            {
                int a, b, c;
                if (std::sscanf(line.c_str(), "%d,%d,%d", &a, &b, &c) == 3)
                    F.push_back({a, b, c});
            }
        }
        std::vector<char> isCopy;
        if (argc > 5)
        {
            std::ifstream in(argv[5]);
            std::string line;
            while (std::getline(in, line))
            {
                int a, b, c, cp;
                if (std::sscanf(line.c_str(), "%d,%d,%d,%d", &a, &b, &c, &cp) == 4)
                    isCopy.push_back(static_cast<char>(cp));
            }
        }
        Mesh mesh(param);
        mesh.setup_from_vertices_faces(V, F, types, mirrors, isCopy);
        mesh.calculate_element_area_volume();
        mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
        mesh.Compute_Energy_And_Force();
        double libE = 0, libA = 0;
        std::vector<int> all, ev;
        std::map<std::string, int> kinds;
        for (int i = 0; i < (int)mesh.faces.size(); i++)
        {
            const Face &f = mesh.faces[i];
            if (f.isBoundary)
                continue;
            const bool ok = f.patchKind == PatchKind::Regular ||
                            f.patchKind == PatchKind::SingleExtraordinary ||
                            (f.patchKind == PatchKind::MultiExtraordinary && f.patchEntry >= 0);
            if (!ok)
                continue;
            libE += f.energy.energyCurvature;
            libA += f.elementArea;
            all.push_back(i);
            if (f.patchKind != PatchKind::Regular)
                ev.push_back(i);
            kinds[f.patchKind == PatchKind::Regular ? "regular"
                  : f.patchKind == PatchKind::SingleExtraordinary ? "single"
                                                                   : "multi"]++;
        }
        std::fprintf(stderr, "faces evaluated %zu:", all.size());
        for (auto &k : kinds)
            std::fprintf(stderr, " %s %d", k.first.c_str(), k.second);
        std::fprintf(stderr, "\n");
        Scheme prod;
        const Total p = eval_faces(prod, mesh, all);
        std::fprintf(stderr, "library E_bend %.12e  fixture prod*kc %.12e  (rel %.1e); area %.12e vs %.12e\n",
                     libE, p.eBend * param.kCurv, p.eBend * param.kCurv / libE - 1, libA, p.area);
        struct DepthSet
        {
            const char *name;
            int d[9];
        };
        const DepthSet sets[] = {
            {"prod", {0, 0, 0, 0, 12, 8, 6, 7, 7}}, {"x0.75", {0, 0, 0, 0, 9, 6, 6, 5, 5}},
            {"x0.5", {0, 0, 0, 0, 6, 4, 6, 4, 4}},  {"5333", {0, 0, 0, 0, 5, 3, 6, 3, 3}},
            {"7555", {0, 0, 0, 0, 7, 5, 6, 5, 5}},  {"8655", {0, 0, 0, 0, 8, 6, 6, 5, 5}},
            {"7444", {0, 0, 0, 0, 7, 4, 6, 4, 4}},  {"8544", {0, 0, 0, 0, 8, 5, 6, 4, 4}},
            {"6555", {0, 0, 0, 0, 6, 5, 6, 5, 5}},  {"8555", {0, 0, 0, 0, 8, 5, 6, 5, 5}},
            {"4333", {0, 0, 0, 0, 4, 3, 6, 3, 3}},  {"4444", {0, 0, 0, 0, 4, 4, 6, 4, 4}},
            {"3333", {0, 0, 0, 0, 3, 3, 6, 3, 3}},  {"5444", {0, 0, 0, 0, 5, 4, 6, 4, 4}},
        };
        Scheme refS;
        for (int N = 4; N <= 8; N++)
            if (N != 6)
                refS.depth[N] = 16;
        refS.weight = tail2_weight;
        refS.areaWeight = [](int N, int D, int d) {
            return d == D - 1 ? 1.0 / (1.0 - lam(N) * lam(N)) : 1.0;
        };
        const Total refAll = eval_faces(refS, mesh, all);
        const Total refEv = eval_faces(refS, mesh, ev);
        if (argc > 6 && std::string(argv[6]) == "perval")
        {
            // error from valence-M patches alone: every other valence at the
            // extrapolated reference
            std::printf("scheme\tN\tD\treduce\tsamples_all\trelE_sheet\trelA_sheet\trelF\n");
            for (const Scheme &fam : families)
            {
                if (fam.name != "plain" && fam.name != "tail2A" && fam.name != "tail2")
                    continue;
                for (int M : {4, 5, 7, 8})
                    for (int D = 2; D <= 14; D++)
                        for (int reduce : {0, 1})
                        {
                            Scheme s = refS;
                            s.depth[M] = D;
                            const auto fw = fam.weight;
                            const auto fa = fam.areaWeight;
                            s.weight = [M, fw](int N, int DD, int d) {
                                return N == M ? fw(N, DD, d) : tail2_weight(N, DD, d);
                            };
                            s.areaWeight = [M, fw, fa](int N, int DD, int d) {
                                if (N != M)
                                    return d == DD - 1 ? 1.0 / (1.0 - lam(N) * lam(N)) : 1.0;
                                return fa ? fa(N, DD, d) : fw(N, DD, d);
                            };
                            // the corner reduction only for the valence under test
                            s.multiCornerReduce = 0;
                            if (reduce)
                            {
                                s.depth[M] = D; // corner children read depth - reduce below
                                s.multiCornerReduce = 1;
                                for (int N = 4; N <= 8; N++)
                                    if (N != M && N != 6)
                                        s.depth[N] = 17; // 17 - 1 = 16, unchanged reference
                            }
                            const Total t = eval_faces(s, mesh, all);
                            // samples spent on valence-M patches
                            Scheme c = s;
                            for (int N = 4; N <= 8; N++)
                                if (N != M && N != 6)
                                    c.depth[N] = 0;
                            const Total tc = eval_faces(c, mesh, all);
                            std::printf("%s\t%d\t%d\t%d\t%ld\t%.3e\t%.3e\t%.3e\n", fam.name.c_str(), M,
                                        D, reduce, tc.samples, (t.eBend - refAll.eBend) / refAll.eBend,
                                        (t.area - refAll.area) / refAll.area, force_err(t, refAll));
                        }
            }
            return 0;
        }
        std::printf("scheme\tdepths\treduce\tsamples_all\tsamples_ev\trelE_sheet\trelE_ev\trelA_sheet\trelF\n");
        for (const Scheme &fam0 : families)
        {
            if (fam0.highOrder)
            {
                const Total t = eval_faces(fam0, mesh, all);
                std::printf("rule12\tprod\t0\t%ld\t-\t%.3e\t-\t%.3e\t%.3e\n", t.samples,
                            (t.eBend - refAll.eBend) / refAll.eBend,
                            (t.area - refAll.area) / refAll.area, force_err(t, refAll));
                continue;
            }
            for (const DepthSet &ds : sets)
                for (int reduce : {0, 1})
                {
                    Scheme s = fam0;
                    for (int N = 4; N <= 8; N++)
                        s.depth[N] = ds.d[N];
                    s.multiCornerReduce = reduce;
                    const Total t = eval_faces(s, mesh, all);
                    const Total te = eval_faces(s, mesh, ev);
                    std::printf("%s\t%s\t%d\t%ld\t%ld\t%.3e\t%.3e\t%.3e\t%.3e\n", fam0.name.c_str(),
                                ds.name, reduce, t.samples, te.samples,
                                (t.eBend - refAll.eBend) / refAll.eBend,
                                (te.eBend - refEv.eBend) / refEv.eBend,
                                (t.area - refAll.area) / refAll.area, force_err(t, refAll));
                }
        }
        return 0;
    }

    if (mode == "latsweep")
    {
        // single-EV faces of one valence on the flipped lattice: depth sweep,
        // relative to that subset's own extrapolated energy
        for (int flipSet : {1, 2})
            for (auto &sh : shapes())
            {
                std::vector<std::vector<double>> V;
                std::vector<std::vector<int>> F;
                build_lattice(7, sh.second, V, F, flipSet);
                Param param = make_param();
                Mesh mesh(param);
                mesh.setup_from_vertices_faces(V, F);
                for (int N : {4, 5, 7, 8})
                {
                    std::vector<int> sel;
                    for (int i = 0; i < (int)mesh.faces.size(); i++)
                        if (mesh.faces[i].patchKind == PatchKind::SingleExtraordinary &&
                            (int)mesh.faces[i].oneRingVertices.size() == N + 6)
                            sel.push_back(i);
                    if (sel.empty())
                        continue;
                    Scheme refS;
                    refS.weight = tail2_weight;
                    const Total ref = eval_faces(with_depth(refS, N, 16), mesh, sel);
                    for (const Scheme &fam : families)
                    {
                        if (fam.name != "plain" && fam.name != "tail2A" && fam.name != "tail_bowl")
                            continue;
                        std::printf("flips%d %-8s N=%d %-9s", flipSet, sh.first.c_str(), N,
                                    fam.name.c_str());
                        for (int D = 2; D <= 12; D++)
                        {
                            const Total t = eval_faces(with_depth(fam, N, D), mesh, sel);
                            std::printf(" %+.1e", (t.eBend - ref.eBend) / ref.eBend);
                        }
                        std::printf("\n");
                    }
                }
            }
        return 0;
    }

    if (mode == "valences")
    {
        for (int flipSet : {1, 2})
        {
            std::vector<std::vector<double>> V;
            std::vector<std::vector<int>> F;
            build_lattice(7, shapes()["bowl"], V, F, flipSet);
            Param param = make_param();
            Mesh mesh(param);
            mesh.setup_from_vertices_faces(V, F);
            std::map<std::string, int> kinds;
            for (auto &f : mesh.faces)
            {
                if (f.patchKind == PatchKind::MultiExtraordinary && f.patchEntry >= 0)
                {
                    auto &e = mesh.multiPatchTable.entry(f.patchEntry);
                    char b[32];
                    std::snprintf(b, 32, "multi %d%d%d", e.valence[0], e.valence[1], e.valence[2]);
                    kinds[b]++;
                }
                else if (f.patchKind == PatchKind::SingleExtraordinary)
                    kinds["single " + std::to_string(f.oneRingVertices.size() - 6)]++;
                else if (f.patchKind == PatchKind::Regular)
                    kinds["regular"]++;
                else
                    kinds["boundary/other"]++;
            }
            std::printf("flipSet %d:\n", flipSet);
            for (auto &k : kinds)
                std::printf("  %-16s %d\n", k.first.c_str(), k.second);
        }
        return 0;
    }

    if (mode == "time")
    {
        // per-face wall time of single-EV evaluation, to check the
        // sample-count cost model (1-point children are cheaper per call?)
        std::printf("N\tscheme\tD\tsamples\tus_per_face\tns_per_sample\n");
        for (int N : {4, 5, 7, 8})
        {
            std::vector<std::vector<double>> V;
            std::vector<std::vector<int>> F;
            build_fan(N, shapes()["generic"], V, F);
            Param param = make_param();
            Mesh mesh(param);
            mesh.setup_from_vertices_faces(V, F);
            std::vector<int> fan;
            for (int i = 0; i < (int)mesh.faces.size(); i++)
                if (mesh.faces[i].patchKind == PatchKind::SingleExtraordinary)
                    fan.push_back(i);
            for (const Scheme &fam : families)
            {
                if (fam.name != "plain" && fam.name != "c1pt_from1" && fam.name != "c1pt_from2" &&
                    fam.name != "c1pt_from3")
                    continue;
                for (int D : {4, 7, 12})
                {
                    const Scheme s = with_depth(fam, N, D);
                    const int reps = 2000;
                    double sink = 0;
                    long samples = 0;
                    const double t0 = now();
                    for (int r = 0; r < reps; r++)
                    {
                        const Total t = eval_faces(s, mesh, fan);
                        sink += t.eBend;
                        samples = t.samples;
                    }
                    const double dt = now() - t0;
                    std::printf("%d\t%s\t%d\t%ld\t%.2f\t%.1f\t%g\n", N, fam.name.c_str(), D,
                                samples / N, dt / reps / N * 1e6,
                                dt / reps / samples * 1e9, sink * 0);
                }
            }
        }
        return 0;
    }
    return 0;
}
