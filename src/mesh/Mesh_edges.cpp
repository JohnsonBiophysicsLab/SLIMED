/**
 * @file Mesh_edges.cpp
 * @brief The edge table and the Monte Carlo edge-flip primitive.
 *
 * A flip swaps the shared edge of two adjacent triangles for the other
 * diagonal of the quadrilateral they form. It is the move that gives a
 * triangulated membrane its in-plane fluidity: without it two vertices that
 * start as neighbours stay neighbours forever, and the sheet carries a shear
 * modulus a lipid bilayer does not have.
 *
 * This file is the topology half of that move -- the half that has to be
 * exactly right before any energy is attached to it. Everything here is
 * connectivity: no coordinates are read, no energy is evaluated, and the
 * acceptance test lives elsewhere.
 *
 * @see docs/edge_flip_plan.md, work package 0
 */

#include "mesh/Mesh.hpp"
#include "mesh/Subdivision_matrices.hpp"

#include <algorithm>
#include <cstdint>
#include <sstream>
#include <stdexcept>

namespace
{
/// Remove the first occurrence of @p value from @p list. Returns whether it was there.
bool erase_first(std::vector<int> &list, int value)
{
    const auto it = std::find(list.begin(), list.end(), value);
    if (it == list.end())
    {
        return false;
    }
    list.erase(it);
    return true;
}

/// Append @p value to @p list unless it is already present.
void push_unique(std::vector<int> &list, int value)
{
    if (std::find(list.begin(), list.end(), value) == list.end())
    {
        list.push_back(value);
    }
}

bool contains(const std::vector<int> &list, int value)
{
    return std::find(list.begin(), list.end(), value) != list.end();
}
} // namespace

void Mesh::build_edge_table()
{
    edges.clear();
    edgeIndex.clear();

    const int nFaces = static_cast<int>(faces.size());
    // A triangle mesh has close to 3F/2 edges; reserving for that keeps the
    // rehashing out of the pass.
    edges.reserve(static_cast<std::size_t>(nFaces) * 2);
    edgeIndex.reserve(static_cast<std::size_t>(nFaces) * 2);

    int nNonManifold = 0;
    for (int iFace = 0; iFace < nFaces; ++iFace)
    {
        const std::vector<int> &corners = faces[iFace].adjacentVertices;
        if (corners.size() != 3)
        {
            continue;
        }
        for (int k = 0; k < 3; ++k)
        {
            const int nodeA = corners[k];
            const int nodeB = corners[(k + 1) % 3];
            const int opposite = corners[(k + 2) % 3];
            const std::uint64_t key = undirected_edge_key(nodeA, nodeB);

            const auto found = edgeIndex.find(key);
            if (found == edgeIndex.end())
            {
                MeshEdge edge;
                edge.index = static_cast<int>(edges.size());
                edge.v[0] = std::min(nodeA, nodeB);
                edge.v[1] = std::max(nodeA, nodeB);
                edge.face[0] = iFace;
                edge.opposite[0] = opposite;
                edgeIndex.emplace(key, edge.index);
                edges.push_back(edge);
            }
            else
            {
                MeshEdge &edge = edges[found->second];
                if (edge.face[1] < 0)
                {
                    edge.face[1] = iFace;
                    edge.opposite[1] = opposite;
                }
                else
                {
                    // A third face on one edge. The record cannot describe it,
                    // and neither can a flip. Say so rather than dropping the
                    // face silently; validate_volume_constraint_topology() is
                    // the designated reporter for a non-manifold mesh and it
                    // only throws when a volume constraint is actually in use.
                    ++nNonManifold;
                }
            }
        }
    }

    if (nNonManifold > 0)
    {
        std::cout << "[Mesh::build_edge_table] WARNING: " << nNonManifold
                  << " edge-face incidence(s) beyond two faces per edge were dropped; "
                     "the mesh is not a two-manifold and those edges cannot be flipped."
                  << std::endl;
    }

    for (int iEdge = 0; iEdge < static_cast<int>(edges.size()); ++iEdge)
    {
        refresh_edge_flippability(iEdge);
    }

    if (param.VERBOSE_MODE)
    {
        int nFlippable = 0;
        for (const MeshEdge &edge : edges)
        {
            nFlippable += edge.flippable ? 1 : 0;
        }
        std::cout << "[Mesh::build_edge_table] " << edges.size() << " edges, " << nFlippable
                  << " structurally flippable." << std::endl;
    }
}

int Mesh::edge_between(int nodeA, int nodeB) const
{
    const auto found = edgeIndex.find(undirected_edge_key(nodeA, nodeB));
    return (found == edgeIndex.end()) ? -1 : found->second;
}

std::vector<int> Mesh::flip_patch_faces(int iEdge) const
{
    std::vector<int> patch;
    if (iEdge < 0 || iEdge >= static_cast<int>(edges.size()))
    {
        return patch;
    }
    const MeshEdge &edge = edges[iEdge];

    // The four vertices the flip touches. Their incident faces are the faces
    // whose control nets change, because a face's energy is a functional of
    // the union of its three corners' one-rings.
    const int quad[4] = {edge.v[0], edge.v[1], edge.opposite[0], edge.opposite[1]};
    for (int corner : quad)
    {
        if (corner < 0)
        {
            continue;
        }
        const std::vector<int> &incident = vertices[corner].adjacentFaces;
        patch.insert(patch.end(), incident.begin(), incident.end());
    }
    std::sort(patch.begin(), patch.end());
    patch.erase(std::unique(patch.begin(), patch.end()), patch.end());
    return patch;
}

void Mesh::refresh_edge_flippability(int iEdge)
{
    MeshEdge &edge = edges[iEdge];
    edge.flippable = false;

    if (edge.is_boundary())
    {
        return;
    }

    // The two vertices being disconnected and the two being joined.
    const int quad[4] = {edge.v[0], edge.v[1], edge.opposite[0], edge.opposite[1]};
    for (int corner : quad)
    {
        if (corner < 0)
        {
            return;
        }
        if (vertices[corner].isGhost)
        {
            return;
        }
        if (!flipFrozenVertex.empty() && flipFrozenVertex[corner] != 0)
        {
            return;
        }
        if (!is_interior_vertex(corner))
        {
            return;
        }
    }

    // Every face whose control net the flip would disturb must itself be a
    // real, evaluable face made of real vertices. This is what keeps the
    // frozen ghost band of a periodic sheet, and any face carrying an
    // insertion (whose spontaneous curvature is indexed by face), out of the
    // move. It costs a walk over roughly eighteen faces per edge, once at
    // build time and once per neighbourhood after an accepted flip.
    for (int iFace : flip_patch_faces(iEdge))
    {
        const Face &face = faces[iFace];
        if (face.isGhost || face.isInsertionPatch)
        {
            return;
        }
        for (int corner : face.adjacentVertices)
        {
            if (vertices[corner].isGhost)
            {
                return;
            }
        }
    }

    edge.flippable = true;
}

bool Mesh::edge_flip_is_admissible(int iEdge, std::string *why) const
{
    auto reject = [&](const std::string &reason) {
        if (why != nullptr)
        {
            *why = reason;
        }
        return false;
    };

    if (iEdge < 0 || iEdge >= static_cast<int>(edges.size()))
    {
        return reject("edge index " + std::to_string(iEdge) + " is out of range");
    }

    const MeshEdge &edge = edges[iEdge];
    const std::string label = "edge " + std::to_string(iEdge);

    if (!edge.flippable)
    {
        return reject(label + " is not structurally flippable: it is on a boundary, or its "
                              "flip patch touches a ghost, frozen or non-interior vertex");
    }

    const int nodeA = edge.v[0];
    const int nodeB = edge.v[1];
    const int target0 = edge.opposite[0];
    const int target1 = edge.opposite[1];

    if (target0 == target1)
    {
        return reject(label + " has the same opposite corner on both sides; flipping it "
                              "would collapse the pair of faces");
    }
    if (edge_between(target0, target1) >= 0)
    {
        return reject(label + " would duplicate the existing edge (" + std::to_string(target0) +
                      ", " + std::to_string(target1) + "), pinching the surface");
    }

    // A flip exchanges which of the two face indices names which triangle --
    // see the note on flip_edge(). That is unobservable as long as the two
    // faces carry the same per-face payload, and spontaneous curvature is the
    // only piece of that payload which is not recomputed from connectivity
    // every step. Enforce the condition rather than assume it: an insertion
    // that made the two sides differ would otherwise migrate by one face on
    // every rejected trial.
    if (faces[edge.face[0]].spontCurvature != faces[edge.face[1]].spontCurvature)
    {
        return reject(label + " separates two faces with different spontaneous curvature; "
                              "flipping it would move the curvature onto the other triangle");
    }

    // The range a flip may leave a vertex in. The parameter can narrow it --
    // a fluid mesh gets cheaper the closer its valences stay to 6 -- but never
    // widen it past what the patch tables can evaluate, because a valence
    // outside that has no patch at all and the face would carry no energy.
    const int minValence = std::max(param.edgeFlipMinValence, kMinIrregularValence);
    const int maxValence = std::min(param.edgeFlipMaxValence, kMaxIrregularValence);

    // Valences after the flip: the two endpoints each lose a neighbour, the two
    // opposite corners each gain one.
    const int valenceAfter[4] = {
        static_cast<int>(vertices[nodeA].adjacentVertices.size()) - 1,
        static_cast<int>(vertices[nodeB].adjacentVertices.size()) - 1,
        static_cast<int>(vertices[target0].adjacentVertices.size()) + 1,
        static_cast<int>(vertices[target1].adjacentVertices.size()) + 1,
    };
    const int quad[4] = {nodeA, nodeB, target0, target1};
    for (int k = 0; k < 4; ++k)
    {
        if (valenceAfter[k] < minValence || valenceAfter[k] > maxValence)
        {
            return reject(label + " would put vertex " + std::to_string(quad[k]) +
                          " at valence " + std::to_string(valenceAfter[k]) +
                          ", outside the supported range [" + std::to_string(minValence) + ", " +
                          std::to_string(maxValence) + "]");
        }
    }

    return true;
}

void Mesh::flip_edge(int iEdge)
{
    if (iEdge < 0 || iEdge >= static_cast<int>(edges.size()))
    {
        throw std::invalid_argument("[Mesh::flip_edge] edge index " + std::to_string(iEdge) +
                                    " is out of range");
    }
    if (edges[iEdge].is_boundary())
    {
        throw std::invalid_argument("[Mesh::flip_edge] edge " + std::to_string(iEdge) +
                                    " is on the mesh boundary and carries only one face");
    }

    const int face0 = edges[iEdge].face[0];
    const int face1 = edges[iEdge].face[1];

    // Name the endpoints by the direction face0 traverses them. The two faces
    // are wound consistently, so exactly one of them carries the directed edge
    // (nodeA -> nodeB); taking face0's direction fixes the orientation of the
    // whole quadrilateral without looking at a single coordinate.
    //
    //        target0                     target0
    //         /   \                       / | \
    //        /     \                     /  |  \
    //   nodeA ----- nodeB    -->    nodeA   |   nodeB
    //        \     /                     \  |  /
    //         \   /                       \ | /
    //        target1                     target1
    //
    // Reading the boundary of the quadrilateral in face order gives the cycle
    // nodeA -> target1 -> nodeB -> target0 -> nodeA, so the two new triangles
    // (nodeA, target1, target0) and (nodeB, target0, target1) are wound along
    // the same cycle the old pair was. Winding is preserved by construction,
    // not by a test.
    const std::vector<int> &face0Corners = faces[face0].adjacentVertices;
    const auto positionOfV0 =
        std::find(face0Corners.begin(), face0Corners.end(), edges[iEdge].v[0]);
    if (positionOfV0 == face0Corners.end() || face0Corners.size() != 3)
    {
        throw std::invalid_argument("[Mesh::flip_edge] edge " + std::to_string(iEdge) +
                                    " is not carried by the face its record names");
    }
    const int slot = static_cast<int>(std::distance(face0Corners.begin(), positionOfV0));

    int nodeA = edges[iEdge].v[0];
    int nodeB = edges[iEdge].v[1];
    if (face0Corners[(slot + 1) % 3] != edges[iEdge].v[1])
    {
        // face0 traverses the edge the other way round.
        std::swap(nodeA, nodeB);
    }

    const int target0 = edges[iEdge].opposite[0]; // the corner of face0
    const int target1 = edges[iEdge].opposite[1]; // the corner of face1

    // The flip patch is invariant under the move: face0 and face1 each keep
    // two of the four quadrilateral corners, so the set of faces incident to
    // {nodeA, nodeB, target0, target1} is the same before and after. Taking it
    // once, before the mutation, means the undo path rebuilds exactly what the
    // do path did.
    const std::vector<int> patch = flip_patch_faces(iEdge);

    // The four edges around the quadrilateral. Their endpoints do not change,
    // so their table entries can be found now and relinked afterwards; only
    // which of the two retriangulated faces holds each of them changes.
    const int edgeAT1 = edge_between(nodeA, target1);
    const int edgeT1B = edge_between(target1, nodeB);
    const int edgeBT0 = edge_between(nodeB, target0);
    const int edgeT0A = edge_between(target0, nodeA);

    // ------------------------------------------------------------------
    // 1. Retriangulate.
    // ------------------------------------------------------------------
    faces[face0].adjacentVertices = {nodeA, target1, target0};
    faces[face1].adjacentVertices = {nodeB, target0, target1};

    // ------------------------------------------------------------------
    // 2. Vertex adjacency.
    // ------------------------------------------------------------------
    erase_first(vertices[nodeA].adjacentVertices, nodeB);
    erase_first(vertices[nodeB].adjacentVertices, nodeA);
    push_unique(vertices[target0].adjacentVertices, target1);
    push_unique(vertices[target1].adjacentVertices, target0);

    // nodeA leaves face1 and nodeB leaves face0; the two targets each pick up
    // the face they were not part of.
    erase_first(vertices[nodeA].adjacentFaces, face1);
    erase_first(vertices[nodeB].adjacentFaces, face0);
    push_unique(vertices[target0].adjacentFaces, face1);
    push_unique(vertices[target1].adjacentFaces, face0);

    // ------------------------------------------------------------------
    // 3. The flipped edge's own record, and its key in the index.
    // ------------------------------------------------------------------
    edgeIndex.erase(undirected_edge_key(nodeA, nodeB));
    {
        MeshEdge &edge = edges[iEdge];
        edge.v[0] = std::min(target0, target1);
        edge.v[1] = std::max(target0, target1);
        // Slot k stays parallel to face[k]: face0 is now (nodeA, target1,
        // target0), whose corner off this edge is nodeA.
        edge.face[0] = face0;
        edge.opposite[0] = nodeA;
        edge.face[1] = face1;
        edge.opposite[1] = nodeB;
    }
    edgeIndex.emplace(undirected_edge_key(target0, target1), iEdge);

    // ------------------------------------------------------------------
    // 4. The four surrounding edges: whichever slot pointed at face0 or face1
    //    is recomputed from the new corner lists.
    // ------------------------------------------------------------------
    auto relink = [&](int iOther) {
        if (iOther < 0)
        {
            return;
        }
        MeshEdge &other = edges[iOther];
        for (int k = 0; k < 2; ++k)
        {
            if (other.face[k] != face0 && other.face[k] != face1)
            {
                continue;
            }
            const int candidates[2] = {face0, face1};
            int chosen = -1;
            for (int candidate : candidates)
            {
                const std::vector<int> &corners = faces[candidate].adjacentVertices;
                if (contains(corners, other.v[0]) && contains(corners, other.v[1]))
                {
                    chosen = candidate;
                    break;
                }
            }
            other.face[k] = chosen;
            other.opposite[k] = -1;
            if (chosen >= 0)
            {
                for (int corner : faces[chosen].adjacentVertices)
                {
                    if (corner != other.v[0] && corner != other.v[1])
                    {
                        other.opposite[k] = corner;
                    }
                }
            }
        }
    };
    relink(edgeAT1);
    relink(edgeT1B);
    relink(edgeBT0);
    relink(edgeT0A);

    // ------------------------------------------------------------------
    // 5. Face::adjacentFaces, for the two retriangulated faces and for the
    //    outer faces across the quadrilateral whose neighbour swapped.
    // ------------------------------------------------------------------
    auto rebuildFaceAdjacency = [&](int iFace) {
        if (iFace < 0)
        {
            return;
        }
        Face &face = faces[iFace];
        // Three slots, value-initialised, packed down from slot 0 in winding
        // order and skipping boundary edges -- the same contract
        // set_adjacent_faces_of_faces() established, which
        // sort_vertices_on_faces() relies on.
        face.adjacentFaces = std::vector<int>(3);
        const std::vector<int> &corners = face.adjacentVertices;
        if (corners.size() != 3)
        {
            return;
        }
        int nFound = 0;
        for (int k = 0; k < 3 && nFound < 3; ++k)
        {
            const int iOther = edge_between(corners[k], corners[(k + 1) % 3]);
            if (iOther < 0)
            {
                continue;
            }
            const MeshEdge &other = edges[iOther];
            for (int s = 0; s < 2; ++s)
            {
                if (other.face[s] >= 0 && other.face[s] != iFace && nFound < 3)
                {
                    face.adjacentFaces[nFound++] = other.face[s];
                }
            }
        }
    };
    rebuildFaceAdjacency(face0);
    rebuildFaceAdjacency(face1);
    for (int iOther : {edgeAT1, edgeT1B, edgeBT0, edgeT0A})
    {
        if (iOther < 0)
        {
            continue;
        }
        for (int s = 0; s < 2; ++s)
        {
            const int neighbour = edges[iOther].face[s];
            if (neighbour >= 0 && neighbour != face0 && neighbour != face1)
            {
                rebuildFaceAdjacency(neighbour);
            }
        }
    }

    // ------------------------------------------------------------------
    // 6. One-rings and flippability over the disturbed neighbourhood.
    // ------------------------------------------------------------------
    for (int iFace : patch)
    {
        build_one_ring_for_face(iFace, nullptr);
    }
    // Only the faces just rebuilt can have acquired a new valence triple.
    ensure_multi_patch_entries(&patch);

    std::vector<int> touchedEdges;
    touchedEdges.reserve(patch.size() * 3);
    for (int iFace : patch)
    {
        const std::vector<int> &corners = faces[iFace].adjacentVertices;
        if (corners.size() != 3)
        {
            continue;
        }
        for (int k = 0; k < 3; ++k)
        {
            const int iOther = edge_between(corners[k], corners[(k + 1) % 3]);
            if (iOther >= 0)
            {
                touchedEdges.push_back(iOther);
            }
        }
    }
    std::sort(touchedEdges.begin(), touchedEdges.end());
    touchedEdges.erase(std::unique(touchedEdges.begin(), touchedEdges.end()), touchedEdges.end());
    for (int iOther : touchedEdges)
    {
        refresh_edge_flippability(iOther);
    }

    ++topologyVersion;
}

bool Mesh::restore_face_connectivity(const std::vector<std::array<int, 3>> &faceCorners,
                                     std::string *why)
{
    auto reject = [&](const std::string &reason) {
        if (why != nullptr)
        {
            *why = reason;
        }
        return false;
    };

    if (faceCorners.size() != faces.size())
    {
        return reject("the connectivity has " + std::to_string(faceCorners.size()) +
                      " faces but the mesh has " + std::to_string(faces.size()) +
                      ". A flip never changes the face count, so this is a different mesh.");
    }

    const int nVertices = static_cast<int>(vertices.size());
    for (int iFace = 0; iFace < static_cast<int>(faceCorners.size()); ++iFace)
    {
        const std::array<int, 3> &corners = faceCorners[iFace];
        for (int k = 0; k < 3; ++k)
        {
            if (corners[k] < 0 || corners[k] >= nVertices)
            {
                return reject("face " + std::to_string(iFace) + " names vertex " +
                              std::to_string(corners[k]) + ", which does not exist");
            }
        }
        if (corners[0] == corners[1] || corners[1] == corners[2] || corners[0] == corners[2])
        {
            return reject("face " + std::to_string(iFace) + " has a repeated corner");
        }
    }

    // Kept so the mesh can be put back exactly as it was if the rebuilt
    // topology turns out not to be a manifold. Half-restoring a mesh is worse
    // than refusing: the caller would go on to integrate a surface whose
    // one-rings and edge table disagree with its faces.
    std::vector<std::vector<int>> previousCorners(faces.size());
    for (int iFace = 0; iFace < static_cast<int>(faces.size()); ++iFace)
    {
        previousCorners[iFace] = faces[iFace].adjacentVertices;
    }

    bool changed = false;
    for (int iFace = 0; iFace < static_cast<int>(faces.size()); ++iFace)
    {
        std::vector<int> &corners = faces[iFace].adjacentVertices;
        corners.assign(faceCorners[iFace].begin(), faceCorners[iFace].end());
        if (corners != previousCorners[iFace])
        {
            changed = true;
        }
    }
    if (!changed)
    {
        // The mesh already had this connectivity. Rebuilding would be correct
        // but would throw away the one-rings and patch tables for nothing, and
        // bumping topologyVersion would invalidate every cache keyed on it.
        return true;
    }

    // The same derivation setup_from_vertices_faces() runs, minus the two
    // steps that do not depend on the adjacency: the ghost flags come from a
    // vertex's grid position, and the pre-refinement has already happened.
    auto rebuild = [this]() {
        // Unsorted, deliberately. The sorted pass reads a vertex's fan off the
        // generated grid's face numbering and drops anything that is not in
        // it, so on a mesh that has flipped -- where a vertex can have seven
        // adjacent faces -- it would silently lose the seventh. flip_edge()
        // maintains these as sets for the same reason, so nothing downstream
        // depends on the fan order once a flip has happened.
        set_adjacent_faces_of_vertices_unsorted();
        set_adjacent_vertices_of_vertices_sorted();
        set_adjacent_faces_of_faces();
        build_edge_table();
        set_one_ring_vertices_sorted();
    };
    auto putBack = [&]() {
        for (int iFace = 0; iFace < static_cast<int>(faces.size()); ++iFace)
        {
            faces[iFace].adjacentVertices = previousCorners[iFace];
        }
        rebuild();
    };

    // The rebuild is not only a source of a wrong answer, it is a source of
    // throws: set_one_ring_vertices_sorted() rejects a mesh whose fans do not
    // close, which is exactly what a corrupt connectivity produces, and it
    // does so before validate_manifold_topology() below ever runs. Catch it
    // and restore, so a bad checkpoint is a refusal rather than a mesh left
    // halfway between two triangulations.
    std::string violation;
    try
    {
        rebuild();
        if (!validate_manifold_topology(&violation))
        {
            putBack();
            return reject("the restored connectivity is not a two-manifold: " + violation);
        }
    }
    catch (const std::exception &error)
    {
        putBack();
        return reject(std::string("the restored connectivity could not be rebuilt: ") +
                      error.what());
    }

    ++topologyVersion;
    return true;
}

bool Mesh::validate_manifold_topology(std::string *why) const
{
    auto reject = [&](const std::string &reason) {
        if (why != nullptr)
        {
            *why = reason;
        }
        return false;
    };

    // 1. Every edge of every face is in the table, carries that face, and is
    //    traversed in opposite directions by its two faces.
    std::unordered_map<std::uint64_t, int> incidence;
    incidence.reserve(faces.size() * 2);
    for (int iFace = 0; iFace < static_cast<int>(faces.size()); ++iFace)
    {
        const std::vector<int> &corners = faces[iFace].adjacentVertices;
        if (corners.size() != 3)
        {
            return reject("face " + std::to_string(iFace) + " does not have three corners");
        }
        if (corners[0] == corners[1] || corners[1] == corners[2] || corners[0] == corners[2])
        {
            return reject("face " + std::to_string(iFace) + " has a repeated corner");
        }
        for (int k = 0; k < 3; ++k)
        {
            ++incidence[undirected_edge_key(corners[k], corners[(k + 1) % 3])];
        }
    }
    for (const auto &entry : incidence)
    {
        if (entry.second > 2)
        {
            return reject("an edge is shared by " + std::to_string(entry.second) +
                          " faces; the mesh is not a two-manifold");
        }
    }

    // 2. Consistent winding: an interior edge must be traversed once in each
    //    direction. Counting directed edges catches a face whose winding was
    //    reversed, which an undirected count cannot see.
    std::unordered_map<std::uint64_t, int> directed;
    directed.reserve(faces.size() * 4);
    for (const Face &face : faces)
    {
        const std::vector<int> &corners = face.adjacentVertices;
        for (int k = 0; k < 3; ++k)
        {
            const std::uint64_t key =
                (static_cast<std::uint64_t>(static_cast<std::uint32_t>(corners[k])) << 32) |
                static_cast<std::uint32_t>(corners[(k + 1) % 3]);
            ++directed[key];
        }
    }
    for (const auto &entry : directed)
    {
        if (entry.second > 1)
        {
            const int from = static_cast<int>(entry.first >> 32);
            const int to = static_cast<int>(entry.first & 0xFFFFFFFFu);
            return reject("directed edge (" + std::to_string(from) + " -> " + std::to_string(to) +
                          ") appears in " + std::to_string(entry.second) +
                          " faces; the winding is inconsistent");
        }
    }

    // 3. The edge table agrees with the face list.
    if (!edges.empty())
    {
        if (edges.size() != incidence.size())
        {
            return reject("the edge table holds " + std::to_string(edges.size()) +
                          " edges but the face list implies " + std::to_string(incidence.size()));
        }
        for (const MeshEdge &edge : edges)
        {
            if (edge.v[0] >= edge.v[1])
            {
                return reject("edge " + std::to_string(edge.index) + " has unsorted endpoints");
            }
            if (edge_between(edge.v[0], edge.v[1]) != edge.index)
            {
                return reject("edge " + std::to_string(edge.index) + " is not reachable through "
                                                                     "edgeIndex under its own key");
            }
            const int expected = incidence.at(undirected_edge_key(edge.v[0], edge.v[1]));
            const int held = (edge.face[0] >= 0 ? 1 : 0) + (edge.face[1] >= 0 ? 1 : 0);
            if (held != expected)
            {
                return reject("edge " + std::to_string(edge.index) + " records " +
                              std::to_string(held) + " face(s) but the face list has " +
                              std::to_string(expected));
            }
            for (int k = 0; k < 2; ++k)
            {
                if (edge.face[k] < 0)
                {
                    continue;
                }
                const std::vector<int> &corners = faces[edge.face[k]].adjacentVertices;
                if (!contains(corners, edge.v[0]) || !contains(corners, edge.v[1]))
                {
                    return reject("edge " + std::to_string(edge.index) + " names face " +
                                  std::to_string(edge.face[k]) + ", which does not carry it");
                }
                if (!contains(corners, edge.opposite[k]) || edge.opposite[k] == edge.v[0] ||
                    edge.opposite[k] == edge.v[1])
                {
                    return reject("edge " + std::to_string(edge.index) +
                                  " has the wrong opposite corner for face " +
                                  std::to_string(edge.face[k]));
                }
            }
        }
    }

    // 4. Vertex adjacency agrees with the face list, and every interior fan is
    //    closed. A vertex whose fan is open is a boundary vertex, which is
    //    legitimate; one whose face and neighbour counts differ by anything
    //    other than 0 or 1 is not.
    std::vector<std::vector<int>> facesOfVertex(vertices.size());
    for (int iFace = 0; iFace < static_cast<int>(faces.size()); ++iFace)
    {
        for (int corner : faces[iFace].adjacentVertices)
        {
            facesOfVertex[corner].push_back(iFace);
        }
    }
    for (int iVertex = 0; iVertex < static_cast<int>(vertices.size()); ++iVertex)
    {
        std::vector<int> recorded = vertices[iVertex].adjacentFaces;
        std::vector<int> actual = facesOfVertex[iVertex];
        std::sort(recorded.begin(), recorded.end());
        std::sort(actual.begin(), actual.end());
        if (recorded != actual)
        {
            return reject("vertex " + std::to_string(iVertex) + " records " +
                          std::to_string(recorded.size()) + " adjacent face(s) but the face list "
                                                            "has " +
                          std::to_string(actual.size()));
        }

        std::vector<int> neighbours;
        for (int iFace : actual)
        {
            for (int corner : faces[iFace].adjacentVertices)
            {
                if (corner != iVertex)
                {
                    push_unique(neighbours, corner);
                }
            }
        }
        std::vector<int> recordedNeighbours = vertices[iVertex].adjacentVertices;
        std::sort(recordedNeighbours.begin(), recordedNeighbours.end());
        std::sort(neighbours.begin(), neighbours.end());
        if (recordedNeighbours != neighbours)
        {
            return reject("vertex " + std::to_string(iVertex) +
                          " records a neighbour set that the face list does not produce");
        }

        const int nFaces = static_cast<int>(actual.size());
        const int nNeighbours = static_cast<int>(neighbours.size());
        if (nFaces > 0 && nNeighbours != nFaces && nNeighbours != nFaces + 1)
        {
            return reject("vertex " + std::to_string(iVertex) + " has " + std::to_string(nFaces) +
                          " face(s) and " + std::to_string(nNeighbours) +
                          " neighbour(s); its fan is neither closed nor a simple boundary fan");
        }
    }

    return true;
}
