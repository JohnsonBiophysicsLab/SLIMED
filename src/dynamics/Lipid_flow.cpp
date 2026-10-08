/**
 * @file Lipid_flow.cpp
 * @brief Lipid flow: the relaxed area of the stretching energy as a
 * coordinate that moves (Param::lipidFlowMode). See docs/lipid_flow_plan.md.
 *
 * Two modes. "global" keeps one relaxed area A0 for the whole sheet and
 * moves it toward the tension a reservoir sets: with no reservoir it is the
 * global constraint exactly, since the modulus is normalised by the relaxed
 * area as built. "local" is the per-face field described below.
 *
 * The membrane's stretching energy is, with the global constraint, a function
 * of the total area alone: one tension for the whole sheet, equalized
 * instantly. That is the limit of infinitely fast lateral lipid flow. Here
 * the flow is given a rate. Every face f that carries a patch has its own
 * relaxed area a_f -- the area its lipids would cover at rest -- and
 *
 *     e_f = (K / 2 a_f) (S_f - a_f)^2 ,     sigma_f = de_f/dS_f = K (S_f/a_f - 1) ,
 *
 * so the force on the control points is the per-face tension times the
 * gradient of the face's own area, which is the kernel's existing area term
 * with a different scalar. The relaxed areas are then coordinates of the
 * model in their own right, with the conjugate force
 *
 *     pi_f = -de_f/da_f = (K / 2) ((S_f/a_f)^2 - 1) = sigma_f (1 + sigma_f / 2K) ,
 *
 * and they move by two mechanisms, each the gradient flow of the same energy:
 *
 *  - lateral flow, conserved: across every edge shared by faces f and g a
 *    flux M w_fg (pi_f - pi_g) into f, with w_fg the edge length over the
 *    distance between the two centroids (the two-point flux of a finite
 *    volume scheme) and M = 1/b the mobility of lipids sliding against
 *    whatever the mesh frame stands for. To linear order this is
 *    d sigma/dt = D_sigma Laplacian(sigma), D_sigma = K M: tension diffuses.
 *
 *  - exchange with a reservoir at tension sigma_res, on the rim or
 *    everywhere: da_f/dt = (k a_f / K) (pi_f - sigma_res), which relaxes the
 *    local tension to the reservoir's at the rate k.
 *
 * With the noise on, each flux carries sqrt(2 kT M w dt) xi and each
 * reservoir exchange sqrt(2 kT (k a/K) dt) xi, and the relaxed areas sample
 * exp(-E/kT): the lateral flow with a conserved total, the reservoir with
 * exp(-(E + sigma_res sum a)/kT). The limits are the models the tree already
 * has: M = 0 is purely local elasticity; M -> infinity with no reservoir
 * makes pi uniform and the energy exactly the global constraint's.
 */
#include "Dynamics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_set>

#include "Counter_rng.hpp"

using slimed::splitmix64;
using slimed::standard_normal;

// ---------------------------------------------------------------------------
// Mesh: the relaxed areas, and the graph the lipids flow on
// ---------------------------------------------------------------------------

bool Mesh::face_carries_lipids(int iFace) const
{
    const Face &face = faces[iFace];
    if (face.isGhost || face.adjacentVertices.size() != 3)
    {
        return false;
    }
    return face.patchKind == PatchKind::Regular ||
           face.patchKind == PatchKind::SingleExtraordinary ||
           (face.patchKind == PatchKind::MultiExtraordinary && face.patchEntry >= 0);
}

void Mesh::initialize_lipid_flow()
{
    // The modulus's normalisation under either mode: the sheet's relaxed
    // area as built.
    lipidAreaNormalization = param.area0;
    const int nFaces = static_cast<int>(faces.size());
    if (!lipid_flow_local())
    {
        faceArea0.clear();
        faceArea0Initial.clear();
        faceTension.clear();
        lipidLinks.clear();
        lipidRimFace.clear();
        return;
    }
    faceArea0.assign(nFaces, 0.0);
    faceTension.assign(nFaces, 0.0);

    double total = 0.0;
    for (int f = 0; f < nFaces; f++)
    {
        if (face_carries_lipids(f))
        {
            total += faces[f].elementArea;
        }
    }
    // The uniform strain the global constraint would start from: a_f in
    // proportion to the face's area, summing to the sheet's relaxed area.
    const double scale = (param.area0 > 0.0 && total > 0.0) ? param.area0 / total : 1.0;
    for (int f = 0; f < nFaces; f++)
    {
        if (face_carries_lipids(f))
        {
            faceArea0[f] = faces[f].elementArea * scale;
        }
    }
    faceArea0Initial = faceArea0;
    build_lipid_links();
}

double Mesh::lipid_face_stretching_energy(int iFace) const
{
    const double modulusArea = lipid_modulus_area(iFace);
    if (!(faceArea0[iFace] > 0.0) || !(modulusArea > 0.0))
    {
        return 0.0;
    }
    const double excess = faces[iFace].elementArea - faceArea0[iFace];
    return 0.5 * param.uSurf / modulusArea * excess * excess;
}

double Mesh::flat_triangle_area(int a, int b, int c) const
{
    // Sort the corners so the same three vertices always give the same
    // bits, whichever face record names them and in whatever order. The
    // flip's pairing rule compares these across the flip, and a last-bit
    // difference between "before" and "after" would break the involution.
    int corner[3] = {a, b, c};
    std::sort(corner, corner + 3);
    const Matrix &p0 = vertices[corner[0]].coord;
    const Matrix &p1 = vertices[corner[1]].coord;
    const Matrix &p2 = vertices[corner[2]].coord;
    const double u[3] = {p1(0, 0) - p0(0, 0), p1(1, 0) - p0(1, 0), p1(2, 0) - p0(2, 0)};
    const double v[3] = {p2(0, 0) - p0(0, 0), p2(1, 0) - p0(1, 0), p2(2, 0) - p0(2, 0)};
    const double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                         u[0] * v[1] - u[1] * v[0]};
    return 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
}

void Mesh::lipid_carry_after_flip(int iEdge, const MeshEdge &before)
{
    if (faceArea0.empty())
    {
        return;
    }
    const MeshEdge &after = edges[iEdge];
    // flip_edge() keeps the two face indices in their slots: face[k] names
    // the same record before and after, with the old corner before.opposite[k]
    // and the new corner after.opposite[k]. The old triangle of record k is
    // (before.v, before.opposite[k]); its new triangle (after.v,
    // after.opposite[k]). The field still sits on the records, so what is
    // decided here is only whether the two records' values swap.
    const double oldArea[2] = {flat_triangle_area(before.v[0], before.v[1], before.opposite[0]),
                               flat_triangle_area(before.v[0], before.v[1], before.opposite[1])};
    const double newArea[2] = {flat_triangle_area(after.v[0], after.v[1], after.opposite[0]),
                               flat_triangle_area(after.v[0], after.v[1], after.opposite[1])};

    bool recordKeepsItsTriangle;
    if (oldArea[0] != oldArea[1] && newArea[0] != newArea[1])
    {
        // The larger old triangle takes the larger new one.
        recordKeepsItsTriangle = (oldArea[0] > oldArea[1]) == (newArea[0] > newArea[1]);
    }
    else
    {
        // flip_edge()'s own pairing, which is an involution too.
        recordKeepsItsTriangle = true;
    }
    if (!recordKeepsItsTriangle)
    {
        std::swap(faceArea0[after.face[0]], faceArea0[after.face[1]]);
        std::swap(faceArea0Initial[after.face[0]], faceArea0Initial[after.face[1]]);
    }
}

std::vector<int> Mesh::lipid_resample_set(int iEdge) const
{
    if (param.lipidFlipResampling == "pair")
    {
        return {edges[iEdge].face[0], edges[iEdge].face[1]};
    }
    return flip_patch_faces(iEdge);
}

double Mesh::lipid_stretching_energy_with_pooled(const std::vector<int> &faceList,
                                                 const std::vector<double> &faceAreas,
                                                 const std::vector<int> &pooled) const
{
    double energy = 0.0;
    double pooledArea = 0.0;
    double pooledRelaxed = 0.0;
    double pooledModulus = 0.0;
    for (std::size_t k = 0; k < faceList.size(); k++)
    {
        const int iFace = faceList[k];
        if (faces[iFace].isGhost || !(faceArea0[iFace] > 0.0))
        {
            continue;
        }
        const double modulusArea = lipid_modulus_area(iFace);
        if (!(modulusArea > 0.0))
        {
            continue;
        }
        if (std::find(pooled.begin(), pooled.end(), iFace) != pooled.end())
        {
            pooledArea += faceAreas[k];
            pooledRelaxed += faceArea0[iFace];
            pooledModulus += modulusArea;
        }
        else
        {
            const double excess = faceAreas[k] - faceArea0[iFace];
            energy += 0.5 * param.uSurf / modulusArea * excess * excess;
        }
    }
    if (pooledModulus > 0.0)
    {
        // The minimum over the split: the pooled quadratic at one tension.
        const double excess = pooledArea - pooledRelaxed;
        energy += 0.5 * param.uSurf / pooledModulus * excess * excess;
    }
    return energy;
}

void Mesh::lipid_resample_after_flip(const std::vector<int> &set,
                                     const std::vector<double> &faceAreas, std::uint64_t key)
{
    double sumArea = 0.0;
    double sumRelaxed = 0.0;
    double sumModulus = 0.0;
    for (std::size_t k = 0; k < set.size(); k++)
    {
        sumArea += faceAreas[k];
        sumRelaxed += faceArea0[set[k]];
        sumModulus += lipid_modulus_area(set[k]);
    }
    if (!(sumModulus > 0.0))
    {
        return;
    }
    // The equal-tension split, and the conserved noise around it: independent
    // N(0, kT c_f / K) per face, less their sum shared out in proportion to
    // c_f, which is the conditional of the independent Gaussians on the sum.
    const double tension = param.uSurf * (sumArea - sumRelaxed) / sumModulus;
    std::vector<double> noise(set.size(), 0.0);
    double sumNoise = 0.0;
    if (param.lipidFlowThermal)
    {
        for (std::size_t k = 0; k < set.size(); k++)
        {
            noise[k] = std::sqrt(param.KBT * lipid_modulus_area(set[k]) / param.uSurf) *
                       standard_normal(key, static_cast<std::uint64_t>(set[k]), 0);
            sumNoise += noise[k];
        }
    }
    for (std::size_t k = 0; k < set.size(); k++)
    {
        const double modulusArea = lipid_modulus_area(set[k]);
        faceArea0[set[k]] = faceAreas[k] - tension * modulusArea / param.uSurf + noise[k] -
                            modulusArea / sumModulus * sumNoise;
    }
}

void Mesh::build_lipid_links()
{
    lipidLinksTopologyVersion = topologyVersion;
    lipidLinks.clear();
    const int nFaces = static_cast<int>(faces.size());
    lipidRimFace.assign(nFaces, 0);

    const auto sourceVertex = [this](int v) {
        return vertices[v].is_periodic_image() ? vertices[v].reflectiveVertexIndex : v;
    };
    const auto sourceKey = [&](int f) {
        std::array<int, 3> key = {{-1, -1, -1}};
        for (int c = 0; c < 3; c++)
        {
            key[c] = sourceVertex(faces[f].adjacentVertices[c]);
        }
        std::sort(key.begin(), key.end());
        return key;
    };

    // The lipid-carrying face each face stands for: itself, or for a periodic
    // copy the carrying face with the same three source corners. A face on an
    // open or clamped boundary, which has no patch, stands for nothing.
    std::map<std::array<int, 3>, int> carrying;
    for (int f = 0; f < nFaces; f++)
    {
        if (face_carries_lipids(f))
        {
            carrying.emplace(sourceKey(f), f);
        }
    }
    std::vector<int> physical(nFaces, -1);
    for (int f = 0; f < nFaces; f++)
    {
        if (faces[f].adjacentVertices.size() != 3)
        {
            continue;
        }
        if (face_carries_lipids(f))
        {
            physical[f] = f;
            continue;
        }
        const auto found = carrying.find(sourceKey(f));
        if (found != carrying.end())
        {
            physical[f] = found->second;
        }
    }

    // One link per physical edge. A seam edge is met twice, once from each
    // side, and an edge between two copies duplicates one met elsewhere; the
    // pair of source endpoints names the physical edge.
    std::unordered_set<std::uint64_t> seen;
    for (const MeshEdge &edge : edges)
    {
        const int fa = edge.face[0];
        const int fb = edge.face[1];
        if (fa < 0 || fb < 0)
        {
            const int f = (fa >= 0) ? fa : fb;
            if (f >= 0 && face_carries_lipids(f))
            {
                lipidRimFace[f] = 1;
            }
            continue;
        }
        if (!face_carries_lipids(fa) && !face_carries_lipids(fb))
        {
            continue;
        }
        const int pa = physical[fa];
        const int pb = physical[fb];
        if (pa < 0 || pb < 0)
        {
            // One side has no lipids at all: the membrane ends here.
            if (pa >= 0)
            {
                lipidRimFace[pa] = 1;
            }
            if (pb >= 0)
            {
                lipidRimFace[pb] = 1;
            }
            continue;
        }
        if (pa == pb)
        {
            continue;
        }
        const std::uint64_t key =
            undirected_edge_key(sourceVertex(edge.v[0]), sourceVertex(edge.v[1]));
        if (!seen.insert(key).second)
        {
            continue;
        }
        LipidLink link;
        link.faceA = pa;
        link.faceB = pb;
        link.localFaceA = fa;
        link.localFaceB = fb;
        link.v0 = edge.v[0];
        link.v1 = edge.v[1];
        lipidLinks.push_back(link);
    }
}

double Mesh::lipid_modulus_area(int iFace) const
{
    return param.lipidFlowLinearStretching ? faceArea0Initial[iFace] : faceArea0[iFace];
}

double Mesh::lipid_face_tension(int iFace) const
{
    const double modulusArea = lipid_modulus_area(iFace);
    if (!(faceArea0[iFace] > 0.0) || !(modulusArea > 0.0))
    {
        return 0.0;
    }
    return param.uSurf * (faces[iFace].elementArea - faceArea0[iFace]) / modulusArea;
}

double Mesh::lipid_exchange_tension(int iFace) const
{
    const double relaxedArea = faceArea0[iFace];
    if (!(relaxedArea > 0.0))
    {
        return 0.0;
    }
    if (param.lipidFlowLinearStretching)
    {
        // Quadratic in a_f: the conjugate of a_f is the tension itself.
        return lipid_face_tension(iFace);
    }
    const double stretch = faces[iFace].elementArea / relaxedArea;
    return 0.5 * param.uSurf * (stretch * stretch - 1.0);
}

double Mesh::total_face_area0() const
{
    double total = 0.0;
    for (double a : faceArea0)
    {
        total += a;
    }
    return total;
}

// ---------------------------------------------------------------------------
// DynamicModel: the relaxed areas as coordinates with their own dynamics
// ---------------------------------------------------------------------------

namespace
{
/// Edge length over centroid distance for every link, on the current net.
void lipid_link_weights(const Mesh &mesh, std::vector<double> &weight)
{
    const int nLinks = static_cast<int>(mesh.lipidLinks.size());
    weight.assign(nLinks, 0.0);
    const auto centroid = [&mesh](int f, double out[3]) {
        out[0] = out[1] = out[2] = 0.0;
        for (int corner : mesh.faces[f].adjacentVertices)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                out[axis] += mesh.vertices[corner].coord(axis, 0) / 3.0;
            }
        }
    };
    for (int l = 0; l < nLinks; l++)
    {
        const Mesh::LipidLink &link = mesh.lipidLinks[l];
        double ca[3];
        double cb[3];
        centroid(link.localFaceA, ca);
        centroid(link.localFaceB, cb);
        double length2 = 0.0;
        double distance2 = 0.0;
        for (int axis = 0; axis < 3; axis++)
        {
            const double e = mesh.vertices[link.v0].coord(axis, 0) - mesh.vertices[link.v1].coord(axis, 0);
            const double d = ca[axis] - cb[axis];
            length2 += e * e;
            distance2 += d * d;
        }
        weight[l] = (distance2 > 0.0) ? std::sqrt(length2 / distance2) : 0.0;
    }
}

bool reservoir_acts_on(const Mesh &mesh, const std::string &region, int f)
{
    if (region == "bulk")
    {
        return true;
    }
    if (region == "rim")
    {
        return mesh.lipidRimFace[f] != 0;
    }
    return false;
}
} // namespace

void DynamicModel::setup_lipid_flow()
{
    lipidDiffusivity = 0.0;
    lipidMobility = 0.0;
    lipidSubsteps = 1;
    lipidReservoirExchanged = 0.0;
    lipidClampedUpdates = 0;

    const Param &param = mesh.param;
    if (!mesh.lipid_flow_on())
    {
        return;
    }
    if (!(param.uSurf > 0.0))
    {
        throw std::runtime_error(
            "[DynamicModel::setup_lipid_flow] lipidFlowMode needs usMembraneStretching > 0: "
            "without a stretching modulus there is no tension to drive the lipids.");
    }
    const std::string &region = param.lipidReservoirRegion;
    if (region != "none" && region != "rim" && region != "bulk")
    {
        throw std::runtime_error(
            "[DynamicModel::setup_lipid_flow] lipidReservoirRegion must be none, rim or bulk; it is \"" +
            region + "\".");
    }
    if (param.lipidReservoirRate < 0.0)
    {
        throw std::runtime_error("[DynamicModel::setup_lipid_flow] lipidReservoirRate must not be negative.");
    }
    if (!(mesh.lipidAreaNormalization > 0.0) || (mesh.lipid_flow_local() && mesh.faceArea0.empty()))
    {
        // Built without the driver's initialization: start from the areas as
        // they stand.
        mesh.calculate_element_area_volume();
        if (!(param.area0 > 0.0))
        {
            double volume = 0.0;
            mesh.sum_membrane_area_and_volume(mesh.param.area0, volume);
        }
        mesh.initialize_lipid_flow();
    }

    const bool reservoirWanted = (region != "none") && param.lipidReservoirRate > 0.0;
    if (mesh.lipid_flow_global())
    {
        // One relaxed area for the sheet. Its rate law is the reservoir's;
        // nothing else moves it.
        std::cout << "[DynamicModel] Lipid flow, global mode: relaxed area " << mesh.param.area0
                  << " nm^2 (modulus normalised by " << mesh.lipidAreaNormalization << " nm^2); "
                  << (reservoirWanted ? "reservoir at " + std::to_string(param.lipidReservoirTension) +
                                            " pN/nm, rate " + std::to_string(param.lipidReservoirRate) +
                                            " per us"
                                      : "no reservoir, so exactly the global constraint")
                  << "; " << (param.lipidFlowThermal ? "thermal" : "deterministic") << "." << std::endl;
        if (reservoirWanted && param.lipidReservoirRate * param.timeStep > 1.0)
        {
            std::cout << "[DynamicModel]   WARNING: lipidReservoirRate x dt above 1; the explicit "
                         "update is close to its stability limit of 2."
                      << std::endl;
        }
        return;
    }

    const int nFaces = static_cast<int>(mesh.faces.size());
    int nCarrying = 0;
    int nRim = 0;
    for (int f = 0; f < nFaces; f++)
    {
        if (mesh.face_carries_lipids(f))
        {
            nCarrying++;
            nRim += mesh.lipidRimFace[f] ? 1 : 0;
        }
    }
    const double total = mesh.total_face_area0();
    const double meanFace = (nCarrying > 0) ? total / nCarrying : 0.0;

    lipidDiffusivity = param.lipidFlowDiffusivity;
    if (lipidDiffusivity < 0.0)
    {
        // The friction the vertices themselves feel, per unit area: kT/D per
        // vertex over the two faces' worth of area a vertex stands for.
        lipidDiffusivity = param.uSurf * (2.0 * meanFace) * param.diffConst / param.KBT;
    }
    lipidMobility = lipidDiffusivity / param.uSurf;

    // The fastest mode of the explicit flux update, bounded by twice the
    // largest weighted degree of the face graph (Gershgorin).
    std::vector<double> weight;
    lipid_link_weights(mesh, weight);
    std::vector<double> degree(nFaces, 0.0);
    for (std::size_t l = 0; l < mesh.lipidLinks.size(); l++)
    {
        degree[mesh.lipidLinks[l].faceA] += weight[l];
        degree[mesh.lipidLinks[l].faceB] += weight[l];
    }
    double fastest = 0.0;
    for (int f = 0; f < nFaces; f++)
    {
        if (mesh.faceArea0[f] > 0.0)
        {
            fastest = std::max(fastest, 2.0 * lipidDiffusivity * degree[f] / mesh.faceArea0[f]);
        }
    }
    const bool reservoir = reservoirWanted;
    const double stiffest = fastest + (reservoir ? param.lipidReservoirRate : 0.0);
    lipidSubsteps = (param.lipidFlowSubsteps > 0)
                        ? param.lipidFlowSubsteps
                        : std::max(1, static_cast<int>(std::ceil(stiffest * param.timeStep / 0.5)));

    std::cout << "[DynamicModel] Lipid flow, local mode: " << nCarrying << " faces carry lipids, mean relaxed area "
              << meanFace << " nm^2, total " << total << " nm^2; " << mesh.lipidLinks.size()
              << " edges carry a flux, " << nRim << " faces on the rim." << std::endl;
    std::cout << "[DynamicModel]   Tension diffusivity D_sigma = " << lipidDiffusivity << " nm^2/us"
              << (param.lipidFlowDiffusivity < 0.0 ? " (from the vertices' friction, K a_v D / kT)" : "")
              << ", mobility " << lipidMobility << " nm^3/(pN us); fastest mode " << fastest
              << " per us, " << lipidSubsteps << " substep" << (lipidSubsteps > 1 ? "s" : "")
              << " per step (rate x dt = " << stiffest * param.timeStep / lipidSubsteps << "); "
              << (param.lipidFlowThermal ? "thermal" : "deterministic") << ", "
              << (param.lipidFlowLinearStretching ? "linear" : "per-lipid") << " stretching." << std::endl;
    if (param.edgeFlipEnabled)
    {
        if (param.lipidFlipResampling == "none")
        {
            std::cout << "[DynamicModel]   Edge flips carry each face's relaxed area with it and are "
                         "charged the change at the local modulus: expect a flip acceptance well "
                         "below the global constraint's (measured 0.7% against 29%). "
                         "lipidFlipResampling = patch (the default with the linear form) restores it."
                      << std::endl;
        }
        else
        {
            std::cout << "[DynamicModel]   Edge flips redraw the relaxed areas of the flipped "
                      << param.lipidFlipResampling << " at the new geometry, total fixed." << std::endl;
        }
    }
    if (reservoir)
    {
        std::cout << "[DynamicModel]   Reservoir on the " << region << " at " << param.lipidReservoirTension
                  << " pN/nm, rate " << param.lipidReservoirRate << " per us." << std::endl;
        if (region == "rim" && nRim == 0)
        {
            std::cout << "[DynamicModel]   WARNING: no face is on a rim (a sheet that wraps both ways has "
                         "none); the rim reservoir does nothing."
                      << std::endl;
        }
    }
    if (stiffest * param.timeStep / lipidSubsteps > 1.0)
    {
        std::cout << "[DynamicModel]   WARNING: rate x dt above 1 -- the explicit flux update is close to "
                     "its stability limit of 2. Raise lipidFlowSubsteps or lower lipidFlowDiffusivity."
                  << std::endl;
    }
}

void DynamicModel::advance_lipid_flow()
{
    const Param &param = mesh.param;
    if (!mesh.lipid_flow_on())
    {
        return;
    }
    if (mesh.lipid_flow_global())
    {
        advance_global_lipid_flow();
        return;
    }
    if (mesh.faceArea0.empty())
    {
        return;
    }
    // An accepted flip or a refinement rewired the faces: rebuild the graph.
    if (mesh.lipidLinksTopologyVersion != mesh.topologyVersion)
    {
        mesh.build_lipid_links();
    }
    const int nFaces = static_cast<int>(mesh.faces.size());
    const int nLinks = static_cast<int>(mesh.lipidLinks.size());
    const double dt = param.timeStep / lipidSubsteps;
    const double kT = param.KBT;
    const bool thermal = param.lipidFlowThermal;
    const std::string &region = param.lipidReservoirRegion;
    const bool reservoir = (region != "none") && param.lipidReservoirRate > 0.0;

    std::vector<double> weight;
    lipid_link_weights(mesh, weight);

    // One key per (run, step), as the vertices use, mixed again per substep
    // so the streams differ from the vertices' and from each other.
    const std::uint64_t stepKey = splitmix64(static_cast<std::uint64_t>(randomSeed)) +
                                  static_cast<std::uint64_t>(this->iteration);

    std::vector<double> exchangeTension(nFaces, 0.0);
    std::vector<double> increment(nFaces, 0.0);
    for (int sub = 0; sub < lipidSubsteps; sub++)
    {
        for (int f = 0; f < nFaces; f++)
        {
            exchangeTension[f] = (mesh.faceArea0[f] > 0.0) ? mesh.lipid_exchange_tension(f) : 0.0;
            increment[f] = 0.0;
        }
        const std::uint64_t noiseKey =
            splitmix64(stepKey ^ (0xA5A5A5A55A5A5A5AULL + 0x9E3779B97F4A7C15ULL * static_cast<std::uint64_t>(sub + 1)));

        // Lateral flow: into the face with the higher exchange tension.
        for (int l = 0; l < nLinks; l++)
        {
            const Mesh::LipidLink &link = mesh.lipidLinks[l];
            const double conductance = lipidMobility * weight[l];
            if (!(conductance > 0.0))
            {
                continue;
            }
            double flux = conductance * (exchangeTension[link.faceA] - exchangeTension[link.faceB]) * dt;
            if (thermal)
            {
                flux += std::sqrt(2.0 * kT * conductance * dt) *
                        standard_normal(noiseKey, static_cast<std::uint64_t>(l), 0);
            }
            increment[link.faceA] += flux;
            increment[link.faceB] -= flux;
        }

        // The reservoir. Its mobility is built on the relaxed area as
        // initialized, so the noise carries no drift correction.
        if (reservoir)
        {
            for (int f = 0; f < nFaces; f++)
            {
                if (!(mesh.faceArea0[f] > 0.0) || !reservoir_acts_on(mesh, region, f))
                {
                    continue;
                }
                const double mobility = param.lipidReservoirRate * mesh.faceArea0Initial[f] / param.uSurf;
                double exchange = mobility * (exchangeTension[f] - param.lipidReservoirTension) * dt;
                if (thermal)
                {
                    exchange += std::sqrt(2.0 * kT * mobility * dt) *
                                standard_normal(noiseKey, static_cast<std::uint64_t>(nLinks + f), 1);
                }
                increment[f] += exchange;
                lipidReservoirExchanged += exchange;
            }
        }

        for (int f = 0; f < nFaces; f++)
        {
            if (!(mesh.faceArea0[f] > 0.0))
            {
                continue;
            }
            double updated = mesh.faceArea0[f] + increment[f];
            // A relaxed area cannot go to zero; in a sane run it stays within
            // a few per cent of the face's area and this never acts.
            const double floor = 0.2 * mesh.faceArea0Initial[f];
            if (!(updated > floor))
            {
                updated = floor;
                lipidClampedUpdates++;
            }
            mesh.faceArea0[f] = updated;
        }
    }
    // The sheet's relaxed area, for whatever still reports it.
    mesh.param.area0 = mesh.total_face_area0();
}

void DynamicModel::advance_global_lipid_flow()
{
    const Param &param = mesh.param;
    const std::string &region = param.lipidReservoirRegion;
    if (region == "none" || !(param.lipidReservoirRate > 0.0) || !(mesh.lipidAreaNormalization > 0.0))
    {
        return; // the global constraint, unchanged
    }
    // dA0/dt = (k c / K)(sigma - sigma_res), sigma = K (A - A0) / c: the
    // sheet's tension relaxes to the reservoir's at rate k. With the noise,
    // sqrt(2 kT k c / K dt) per step: the mobility k c / K is a constant.
    const double normalization = mesh.lipidAreaNormalization;
    const double tension = param.uSurf * (param.area - param.area0) / normalization;
    const double mobility = param.lipidReservoirRate * normalization / param.uSurf;
    double exchange = mobility * (tension - param.lipidReservoirTension) * param.timeStep;
    if (param.lipidFlowThermal)
    {
        const std::uint64_t stepKey = splitmix64(static_cast<std::uint64_t>(randomSeed)) +
                                      static_cast<std::uint64_t>(this->iteration);
        const std::uint64_t noiseKey = splitmix64(stepKey ^ 0xA5A5A5A55A5A5A5AULL);
        exchange += std::sqrt(2.0 * param.KBT * mobility * param.timeStep) *
                    standard_normal(noiseKey, 0, 2);
    }
    double updated = param.area0 + exchange;
    const double floor = 0.2 * normalization;
    if (!(updated > floor))
    {
        updated = floor;
        lipidClampedUpdates++;
    }
    lipidReservoirExchanged += updated - param.area0;
    mesh.param.area0 = updated;
}

void DynamicModel::append_lipid_log(const std::string &filepath, int iteration)
{
    std::ofstream out(filepath, std::ios::app);
    if (!out)
    {
        throw std::runtime_error("[DynamicModel::append_lipid_log] cannot open " + filepath);
    }
    if (!lipidLogHeaderWritten)
    {
        out << "iteration,time_us,sumRelaxedArea_nm2,sumArea_nm2,meanTension_pN_nm,minTension_pN_nm,"
               "maxTension_pN_nm,stdTension_pN_nm,rimTension_pN_nm,reservoirExchanged_nm2,clampedUpdates,"
               "aggregateTension_pN_nm\n";
        lipidLogHeaderWritten = true;
    }
    if (mesh.lipid_flow_global())
    {
        const double normalization = mesh.lipidAreaNormalization;
        const double tension =
            (normalization > 0.0) ? mesh.param.uSurf * (mesh.param.area - mesh.param.area0) / normalization : 0.0;
        out << iteration << ',' << std::setprecision(12) << iteration * mesh.param.timeStep << ','
            << mesh.param.area0 << ',' << mesh.param.area << ',' << tension << ',' << tension << ','
            << tension << ',' << 0.0 << ',' << tension << ',' << lipidReservoirExchanged << ','
            << lipidClampedUpdates << ',' << tension << '\n';
        return;
    }
    const int nFaces = static_cast<int>(mesh.faces.size());
    double sumRelaxed = 0.0;
    double sumArea = 0.0;
    double sum = 0.0;
    double sumSquare = 0.0;
    double lowest = std::numeric_limits<double>::infinity();
    double highest = -std::numeric_limits<double>::infinity();
    double rimSum = 0.0;
    int nRim = 0;
    int n = 0;
    for (int f = 0; f < nFaces; f++)
    {
        const double relaxed = mesh.faceArea0.empty() ? 0.0 : mesh.faceArea0[f];
        if (!(relaxed > 0.0))
        {
            continue;
        }
        const double area = mesh.faces[f].elementArea;
        const double tension = mesh.lipid_face_tension(f);
        sumRelaxed += relaxed;
        sumArea += area;
        sum += tension;
        sumSquare += tension * tension;
        lowest = std::min(lowest, tension);
        highest = std::max(highest, tension);
        if (mesh.lipidRimFace[f])
        {
            rimSum += tension;
            nRim++;
        }
        n++;
    }
    // Two means. The face average of sigma_f is what the shape feels: it is
    // the multiplier the force carries, averaged. The aggregate
    // K (sum S / sum a - 1) is the tension of the pooled lipids; the two
    // coincide unless the relaxed areas fluctuate.
    const double aggregate = (sumRelaxed > 0.0) ? mesh.param.uSurf * (sumArea / sumRelaxed - 1.0) : 0.0;
    const double average = (n > 0) ? sum / n : 0.0;
    const double variance = (n > 0) ? std::max(0.0, sumSquare / n - average * average) : 0.0;
    out << iteration << ',' << std::setprecision(12) << iteration * mesh.param.timeStep << ',' << sumRelaxed
        << ',' << sumArea << ',' << average << ',' << (n > 0 ? lowest : 0.0) << ',' << (n > 0 ? highest : 0.0)
        << ',' << std::sqrt(variance) << ',' << (nRim > 0 ? rimSum / nRim : 0.0) << ','
        << lipidReservoirExchanged << ',' << lipidClampedUpdates << ',' << aggregate << '\n';
}

void DynamicModel::append_lipid_face_frame(const std::string &filepath)
{
    std::ofstream out(filepath, std::ios::app);
    if (!out)
    {
        throw std::runtime_error("[DynamicModel::append_lipid_face_frame] cannot open " + filepath);
    }
    out << std::setprecision(8);
    const int nFaces = static_cast<int>(mesh.faces.size());
    for (int f = 0; f < nFaces; f++)
    {
        out << (mesh.faceArea0.empty() ? 0.0 : mesh.faceArea0[f]) << ',' << mesh.faces[f].elementArea << ',';
    }
    out << '\n';
}
