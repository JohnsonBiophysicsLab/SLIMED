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
#include "mesh/Mesh.hpp"
#include "energy_force/Patch_kernel.hpp"
#include "energy_force/Patch_rows_flat.hpp"
#include "mesh/Gauss_quadrature.hpp"
#include "mesh/Irregular_patch_rows.hpp"
#include "mesh/Multi_extraordinary_patch.hpp"
#include "io/io.hpp"
#include <fstream>
#include <sstream>

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
