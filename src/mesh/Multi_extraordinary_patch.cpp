/**
 * @file Multi_extraordinary_patch.cpp
 * @brief The generic multi-extraordinary patch and its prolongation matrices.
 *
 * @see include/mesh/Multi_extraordinary_patch.hpp for what this computes and why.
 */

#include "mesh/Multi_extraordinary_patch.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>

#include "mesh/Subdivision_matrices.hpp"

namespace
{
using EdgeKey = std::pair<int, int>;

EdgeKey edge_key(int a, int b)
{
    return (a < b) ? EdgeKey{a, b} : EdgeKey{b, a};
}

void require_valence(int valence, int corner)
{
    if (valence < kMinIrregularValence || valence > kMaxIrregularValence)
    {
        throw std::invalid_argument(
            "[build_generic_face_patch] corner " + std::to_string(corner) + " has valence " +
            std::to_string(valence) + ", outside the supported range [" +
            std::to_string(kMinIrregularValence) + ", " + std::to_string(kMaxIrregularValence) +
            "]. Unsupported topology is an error, not a fallback.");
    }
}
} // namespace

GenericFacePatch build_generic_face_patch(int valence0, int valence1, int valence2)
{
    require_valence(valence0, 0);
    require_valence(valence1, 1);
    require_valence(valence2, 2);

    const int N[3] = {valence0, valence1, valence2};
    const int extras[3] = {N[0] - 4, N[1] - 4, N[2] - 4};

    GenericFacePatch patch;
    patch.valence = {{N[0], N[1], N[2]}};
    patch.nVertices = N[0] + N[1] + N[2] - 6;

    // Numbering. The corners come first so that a face's own winding fixes
    // them; then the three vertices opposite its edges, which are shared
    // between two corners' fans; then each corner's private extras.
    const int corner[3] = {0, 1, 2};
    const int opposite01 = 3;
    const int opposite12 = 4;
    const int opposite20 = 5;
    const int edgeOpposite[3] = {opposite01, opposite12, opposite20}; // opposite edge (a, a+1)

    int nextExtra = 6;
    std::array<std::vector<int>, 3> extraOf;
    for (int a = 0; a < 3; a++)
    {
        for (int k = 0; k < extras[a]; k++)
        {
            extraOf[a].push_back(nextExtra++);
        }
    }
    if (nextExtra != patch.nVertices)
    {
        throw std::logic_error("[build_generic_face_patch] numbered " + std::to_string(nextExtra) +
                               " vertices but the patch should have " +
                               std::to_string(patch.nVertices));
    }

    // Corner a's fan, in winding order starting at corner a+1:
    //
    //     [ c_{a+1}, c_{a+2}, e_{a+2,a}, extras..., e_{a,a+1} ]
    //
    // Reading it as consecutive pairs gives corner a's faces:
    // (a, fan[k], fan[k+1]) for every k, cyclically.
    for (int a = 0; a < 3; a++)
    {
        const int next = (a + 1) % 3;
        const int prev = (a + 2) % 3;
        std::vector<int> &fan = patch.cornerFan[a];
        fan.push_back(corner[next]);
        fan.push_back(corner[prev]);
        fan.push_back(edgeOpposite[prev]); // the vertex opposite the edge (prev, a)
        for (int extra : extraOf[a])
        {
            fan.push_back(extra);
        }
        fan.push_back(edgeOpposite[a]); // the vertex opposite the edge (a, next)
        if (static_cast<int>(fan.size()) != N[a])
        {
            throw std::logic_error("[build_generic_face_patch] corner " + std::to_string(a) +
                                   " fan has " + std::to_string(fan.size()) + " entries, expected " +
                                   std::to_string(N[a]));
        }
    }

    // Faces, taken from the fans and deduplicated. Every face is (a, fan[k],
    // fan[k+1]) for some corner a; a face touching two corners is produced by
    // both, so it is filed under its sorted corner set.
    std::map<std::array<int, 3>, std::array<int, 3>> unique;
    for (int a = 0; a < 3; a++)
    {
        const std::vector<int> &fan = patch.cornerFan[a];
        for (int k = 0; k < N[a]; k++)
        {
            const std::array<int, 3> face{{a, fan[k], fan[(k + 1) % N[a]]}};
            std::array<int, 3> key = face;
            std::sort(key.begin(), key.end());
            unique.emplace(key, face);
        }
    }
    for (const auto &entry : unique)
    {
        patch.faces.push_back(entry.second);
    }

    const int expectedFaces = N[0] + N[1] + N[2] - 5;
    if (static_cast<int>(patch.faces.size()) != expectedFaces)
    {
        throw std::logic_error("[build_generic_face_patch] built " +
                               std::to_string(patch.faces.size()) + " faces at valences (" +
                               std::to_string(N[0]) + ", " + std::to_string(N[1]) + ", " +
                               std::to_string(N[2]) + "), expected " +
                               std::to_string(expectedFaces));
    }

    return patch;
}

void MultiPatchTable::clear()
{
    entries_.clear();
    index_.clear();
    buffer_.clear();
}

int MultiPatchTable::ensure(int valence0, int valence1, int valence2)
{
    const std::array<int, 3> key{{valence0, valence1, valence2}};
    const auto found = index_.find(key);
    if (found != index_.end())
    {
        return found->second;
    }

    const GenericFacePatch patch = build_generic_face_patch(valence0, valence1, valence2);
    const int K = patch.nVertices;

    // ---- the local mesh's adjacency, read straight off the face list -------
    std::vector<std::vector<int>> neighbours(K);
    std::map<EdgeKey, std::vector<int>> oppositesOfEdge;
    std::vector<int> faceCount(K, 0);
    for (const std::array<int, 3> &face : patch.faces)
    {
        for (int k = 0; k < 3; k++)
        {
            const int a = face[k];
            const int b = face[(k + 1) % 3];
            const int opposite = face[(k + 2) % 3];
            faceCount[a]++;
            oppositesOfEdge[edge_key(a, b)].push_back(opposite);
            if (std::find(neighbours[a].begin(), neighbours[a].end(), b) == neighbours[a].end())
            {
                neighbours[a].push_back(b);
            }
            if (std::find(neighbours[b].begin(), neighbours[b].end(), a) == neighbours[b].end())
            {
                neighbours[b].push_back(a);
            }
        }
    }

    // ---- child points of one subdivision -----------------------------------
    //
    // A vertex point exists only where the fan closes inside the patch, which
    // is the three corners and nothing else -- everything further out is the
    // boundary of this local mesh. An edge point exists on every interior
    // edge.
    struct ChildPoint
    {
        bool isVertexPoint = false;
        int a = -1;
        int b = -1;
    };
    std::vector<ChildPoint> childPoints;
    std::vector<int> childOfVertex(K, -1);
    std::map<EdgeKey, int> childOfEdge;

    for (int v = 0; v < 3; v++)
    {
        if (faceCount[v] != static_cast<int>(neighbours[v].size()))
        {
            throw std::logic_error("[MultiPatchTable] corner " + std::to_string(v) +
                                   " of the generic patch does not have a closed fan");
        }
        childOfVertex[v] = static_cast<int>(childPoints.size());
        childPoints.push_back({true, v, -1});
    }
    for (const auto &entry : oppositesOfEdge)
    {
        if (entry.second.size() != 2)
        {
            continue; // boundary edge of the local mesh: no child point
        }
        childOfEdge[entry.first] = static_cast<int>(childPoints.size());
        childPoints.push_back({false, entry.first.first, entry.first.second});
    }
    const int nChildPoints = static_cast<int>(childPoints.size());

    auto edgeChild = [&](int a, int b) {
        const auto it = childOfEdge.find(edge_key(a, b));
        return (it == childOfEdge.end()) ? -1 : it->second;
    };

    // ---- A_loc: every child point as a combination of the K parent points --
    std::vector<double> aLoc(static_cast<std::size_t>(nChildPoints) * K, 0.0);
    for (int row = 0; row < nChildPoints; row++)
    {
        double *target = &aLoc[static_cast<std::size_t>(row) * K];
        const ChildPoint &child = childPoints[row];
        if (child.isVertexPoint)
        {
            const int degree = static_cast<int>(neighbours[child.a].size());
            const double beta = loop_vertex_weight(degree);
            target[child.a] = 1.0 - degree * beta;
            for (int neighbour : neighbours[child.a])
            {
                target[neighbour] = beta;
            }
        }
        else
        {
            target[child.a] += 3.0 / 8.0;
            target[child.b] += 3.0 / 8.0;
            for (int opposite : oppositesOfEdge.at(edge_key(child.a, child.b)))
            {
                target[opposite] += 1.0 / 8.0;
            }
        }
    }

    // ---- adjacency of the subdivided local mesh ----------------------------
    std::vector<std::vector<int>> subNeighbours(nChildPoints);
    auto link = [&](int a, int b) {
        if (std::find(subNeighbours[a].begin(), subNeighbours[a].end(), b) ==
            subNeighbours[a].end())
        {
            subNeighbours[a].push_back(b);
        }
    };
    for (const std::array<int, 3> &face : patch.faces)
    {
        const int va = childOfVertex[face[0]];
        const int vb = childOfVertex[face[1]];
        const int vc = childOfVertex[face[2]];
        const int eab = edgeChild(face[0], face[1]);
        const int ebc = edgeChild(face[1], face[2]);
        const int eca = edgeChild(face[2], face[0]);
        const std::array<std::array<int, 3>, 4> subFaces{
            {{va, eab, eca}, {eab, vb, ebc}, {eca, ebc, vc}, {eab, ebc, eca}}};
        for (const std::array<int, 3> &sub : subFaces)
        {
            if (sub[0] < 0 || sub[1] < 0 || sub[2] < 0)
            {
                continue; // hangs off the boundary of the local mesh
            }
            link(sub[0], sub[1]);
            link(sub[1], sub[0]);
            link(sub[1], sub[2]);
            link(sub[2], sub[1]);
            link(sub[2], sub[0]);
            link(sub[0], sub[2]);
        }
    }

    // The common neighbour of node1 and node2 that is not node3 -- the same
    // rule Mesh::find_opposite_node_index() applies on the real mesh.
    auto opposite = [&](int node1, int node2, int node3) {
        if (node1 < 0 || node2 < 0)
        {
            return -1;
        }
        for (int candidate : subNeighbours[node1])
        {
            if (candidate == node3)
            {
                continue;
            }
            if (std::find(subNeighbours[node2].begin(), subNeighbours[node2].end(), candidate) !=
                subNeighbours[node2].end())
            {
                return candidate;
            }
        }
        return -1;
    };

    /**
     * The one-ring of a child triangle whose extraordinary corner is `d4c`,
     * with `fan` its neighbours in winding order starting at `d7c`. Produced
     * in the internal numbering CanonicalPatch documents -- 0 is d4, 1..N the
     * fan, then d6, d9, d10, d11, d12 -- and then permuted into the canonical
     * column order, which is exactly what the mesh does for a real face. The
     * same walk serves a valence-6 child, where it yields the familiar 12.
     */
    auto childOneRing = [&](int d4c, const std::vector<int> &fan) {
        const int valence = static_cast<int>(fan.size());
        std::vector<int> ring(valence + 6, -1);
        ring[0] = d4c;
        for (int k = 0; k < valence; k++)
        {
            ring[1 + k] = fan[k];
        }
        const int d7c = fan[0];
        const int d8c = fan[1];
        const int d5c = fan[2];
        const int d3c = fan[valence - 1];
        const int d11c = opposite(d7c, d8c, d4c);
        ring[valence + 1] = opposite(d3c, d7c, d4c);  // d6
        ring[valence + 2] = opposite(d8c, d5c, d4c);  // d9
        ring[valence + 3] = opposite(d7c, d11c, d8c); // d10
        ring[valence + 4] = d11c;
        ring[valence + 5] = opposite(d8c, d11c, d7c); // d12

        const std::vector<int> columnOf = canonical_control_order(valence);
        std::vector<int> canonical(valence + 6, -1);
        for (int internal = 0; internal < valence + 6; internal++)
        {
            canonical[columnOf[internal]] = ring[internal];
        }
        return canonical;
    };

    Entry entry;
    entry.valence = key;
    entry.nControl = K;

    // ---- the three corner children -----------------------------------------
    //
    // The 1-to-4 split puts corner a's child at (V(ca), E(ca, c_{a+1}),
    // E(ca, c_{a+2})). The parent's own fan for corner a starts at c_{a+1} and
    // runs to c_{a+2}, so mapping it through E(ca, .) gives the child's fan in
    // exactly the order the walk above wants.
    std::array<std::vector<int>, 4> rings;
    std::array<int, 4> ringValence{{0, 0, 0, 0}};
    for (int a = 0; a < 3; a++)
    {
        std::vector<int> fan;
        for (int neighbour : patch.cornerFan[a])
        {
            const int point = edgeChild(a, neighbour);
            if (point < 0)
            {
                throw std::logic_error("[MultiPatchTable] the edge from corner " +
                                       std::to_string(a) + " to " + std::to_string(neighbour) +
                                       " has no child point");
            }
            fan.push_back(point);
        }
        rings[a] = childOneRing(childOfVertex[a], fan);
        ringValence[a] = patch.valence[a];
    }

    // ---- the centre child ---------------------------------------------------
    //
    // (E01, E12, E20), all three at valence 6 in the subdivided mesh, so it is
    // an ordinary regular patch and takes the direct kernel.
    {
        const int e01 = edgeChild(0, 1);
        const int e12 = edgeChild(1, 2);
        const int e20 = edgeChild(2, 0);
        if (e01 < 0 || e12 < 0 || e20 < 0)
        {
            throw std::logic_error("[MultiPatchTable] the centre child is missing a corner");
        }
        // Walk E01's fan starting at E12 -- the winding of the centre child.
        std::vector<int> fan;
        fan.push_back(e12);
        fan.push_back(e20);
        const int degree = static_cast<int>(subNeighbours[e01].size());
        if (degree != 6)
        {
            throw std::logic_error("[MultiPatchTable] the centre child's corner E01 has valence " +
                                   std::to_string(degree) + ", expected 6");
        }
        for (int k = 2; k < degree; k++)
        {
            const int next = opposite(e01, fan[k - 1], fan[k - 2]);
            if (next < 0)
            {
                throw std::logic_error("[MultiPatchTable] the centre child's fan does not close");
            }
            fan.push_back(next);
        }
        rings[3] = childOneRing(e01, fan);
        ringValence[3] = 6;
    }

    // ---- prolongation matrices ---------------------------------------------
    for (int c = 0; c < 4; c++)
    {
        const std::vector<int> &ring = rings[c];
        for (int point : ring)
        {
            if (point < 0)
            {
                throw std::logic_error("[MultiPatchTable] child " + std::to_string(c) +
                                       " of the (" + std::to_string(valence0) + ", " +
                                       std::to_string(valence1) + ", " + std::to_string(valence2) +
                                       ") patch has an incomplete one-ring");
            }
        }

        Child child;
        child.valence = ringValence[c];
        child.nControl = static_cast<int>(ring.size());
        child.offset = buffer_.size();
        buffer_.resize(buffer_.size() + static_cast<std::size_t>(child.nControl) * K, 0.0);

        // Row i of M is row ring[i] of A_loc: the child's i-th control point,
        // written in the parent's control points.
        for (int i = 0; i < child.nControl; i++)
        {
            const double *source = &aLoc[static_cast<std::size_t>(ring[i]) * K];
            double *target = &buffer_[child.offset + static_cast<std::size_t>(i) * K];
            for (int j = 0; j < K; j++)
            {
                target[j] = source[j];
            }
        }
        entry.children[c] = child;
    }

    const int newIndex = static_cast<int>(entries_.size());
    entries_.push_back(entry);
    index_.emplace(key, newIndex);
    return newIndex;
}
