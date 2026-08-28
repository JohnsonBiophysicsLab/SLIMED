#include "mesh/Subdivision_matrices.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>

namespace
{
/// An undirected edge, keyed so the two orientations collide.
using EdgeKey = std::pair<int, int>;

EdgeKey edge_key(int a, int b)
{
    return (a < b) ? EdgeKey{a, b} : EdgeKey{b, a};
}

/// The two faces on an edge, and the corner each of them puts opposite it.
struct EdgeRecord
{
    std::vector<int> opposite;
};

/**
 * @brief A child point of the subdivided patch.
 *
 * Loop subdivision produces one vertex point per original vertex with a
 * complete ring, and one edge point per interior edge. `isVertexPoint`
 * distinguishes them; `a`/`b` name the source vertex or the edge endpoints.
 */
struct ChildPoint
{
    bool isVertexPoint = false;
    int a = -1;
    int b = -1;
    std::array<double, 2> position{{0.0, 0.0}};
};
} // namespace

CanonicalPatch build_canonical_patch(int valence)
{
    const int N = valence;
    CanonicalPatch patch;
    patch.valence = N;
    patch.nVertices = N + 6;

    // Named indices, matching the d1..d12 scheme of the regular patch.
    const int d4 = 0;
    const int d7 = 1;
    const int d8 = 2;
    const int d5 = 3;
    const int d3 = N;      // last entry of the fan
    const int d6 = N + 1;
    const int d9 = N + 2;
    const int d10 = N + 3;
    const int d11 = N + 4;
    const int d12 = N + 5;

    // The fan around d4, cyclic: (d4, fan[i], fan[i+1]).
    for (int i = 0; i < N; i++)
    {
        patch.faces.push_back({d4, 1 + i, 1 + ((i + 1) % N)});
    }
    // d7's four remaining faces (it is valence 6; two are already in the fan).
    patch.faces.push_back({d3, d6, d7});
    patch.faces.push_back({d6, d10, d7});
    patch.faces.push_back({d10, d11, d7});
    patch.faces.push_back({d11, d8, d7});
    // d8's three remaining faces (two in the fan, one shared with d7).
    patch.faces.push_back({d5, d8, d9});
    patch.faces.push_back({d9, d8, d12});
    patch.faces.push_back({d12, d8, d11});

    // Canonical 2D positions, used only to order the children. They reproduce
    // the standard drawing of the regular 12-point patch:
    //
    //          d1    d2                 <- extras live on this row
    //       d3    d4    d5
    //    d6    d7    d8    d9
    //       d10   d11   d12
    //
    // so that ordering children by descending row then ascending column
    // reproduces the child ordering of the literal matrices this replaces.
    patch.layout.assign(patch.nVertices, {{0.0, 0.0}});
    patch.layout[d4] = {{1.5, 2.0}};
    patch.layout[d7] = {{1.0, 1.0}};
    patch.layout[d8] = {{2.0, 1.0}};
    patch.layout[d5] = {{2.5, 2.0}};
    patch.layout[d3] = {{0.5, 2.0}};
    patch.layout[d6] = {{0.0, 1.0}};
    patch.layout[d9] = {{3.0, 1.0}};
    patch.layout[d10] = {{0.5, 0.0}};
    patch.layout[d11] = {{1.5, 0.0}};
    patch.layout[d12] = {{2.5, 0.0}};

    // The N-4 extra fan vertices occupy the top row, right to left, following
    // the fan's counter-clockwise order. At N=6 they land on d2 and d1.
    const int nExtras = N - 4;
    for (int k = 0; k < nExtras; k++)
    {
        const int index = 4 + k; // fan positions 3..N-2 -> vertex indices 4..N-1
        const double x = (nExtras == 1) ? 2.0 : 2.0 - 1.0 * k / (nExtras - 1);
        patch.layout[index] = {{x, 3.0}};
    }

    return patch;
}

SubdivisionMatrices build_subdivision_matrices(int valence)
{
    if (valence < kMinIrregularValence || valence > kMaxIrregularValence)
    {
        throw std::invalid_argument(
            "[build_subdivision_matrices] valence " + std::to_string(valence) +
            " is outside the supported range [" + std::to_string(kMinIrregularValence) + ", " +
            std::to_string(kMaxIrregularValence) +
            "]. Unsupported topology is an error, not a fallback. "
            "See docs/irregular_patch_valence_4_to_8_plan.md.");
    }

    const int N = valence;
    const CanonicalPatch patch = build_canonical_patch(N);
    const int nControl = patch.nVertices; // N + 6

    // Adjacency and edge->opposite-corner tables, both read straight off the
    // face list so nothing is transcribed by hand.
    std::vector<std::vector<int>> neighbours(nControl);
    std::map<EdgeKey, EdgeRecord> edges;
    for (const std::array<int, 3> &face : patch.faces)
    {
        for (int k = 0; k < 3; k++)
        {
            const int a = face[k];
            const int b = face[(k + 1) % 3];
            const int opposite = face[(k + 2) % 3];
            edges[edge_key(a, b)].opposite.push_back(opposite);
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

    // Only a vertex whose fan closes inside the patch has a defined vertex
    // point. That is d4, d7 and d8 and nothing else -- the rest of the patch
    // is the boundary of the local mesh.
    std::vector<int> faceCount(nControl, 0);
    for (const std::array<int, 3> &face : patch.faces)
    {
        for (int k = 0; k < 3; k++)
            faceCount[face[k]]++;
    }

    std::vector<ChildPoint> children;
    for (int v = 0; v < nControl; v++)
    {
        if (faceCount[v] == static_cast<int>(neighbours[v].size()))
        {
            ChildPoint child;
            child.isVertexPoint = true;
            child.a = v;
            child.position = patch.layout[v];
            children.push_back(child);
        }
    }
    for (const auto &entry : edges)
    {
        if (entry.second.opposite.size() != 2)
            continue; // boundary edge of the local patch, no child
        ChildPoint child;
        child.isVertexPoint = false;
        child.a = entry.first.first;
        child.b = entry.first.second;
        child.position = {{0.5 * (patch.layout[child.a][0] + patch.layout[child.b][0]),
                           0.5 * (patch.layout[child.a][1] + patch.layout[child.b][1])}};
        children.push_back(child);
    }

    if (static_cast<int>(children.size()) != N + 12)
    {
        throw std::logic_error("[build_subdivision_matrices] expected " + std::to_string(N + 12) +
                               " child points at valence " + std::to_string(N) + ", built " +
                               std::to_string(children.size()));
    }

    // Order children by row from the top, then left to right. This is the
    // ordering the literal matrices use; deriving it from the layout is what
    // lets the same rule serve every valence.
    std::sort(children.begin(), children.end(),
              [](const ChildPoint &lhs, const ChildPoint &rhs) {
                  if (lhs.position[1] != rhs.position[1])
                      return lhs.position[1] > rhs.position[1];
                  return lhs.position[0] < rhs.position[0];
              });

    // Columns must be in the same canonical order the mesh uses for a one-ring
    // -- {d2..d12} at N=5 -- which is the same "row from the top, then left to
    // right" rule that orders the children. Sorting the layout gives both.
    std::vector<int> controlOrder(nControl);
    for (int v = 0; v < nControl; v++)
        controlOrder[v] = v;
    std::sort(controlOrder.begin(), controlOrder.end(), [&](int lhs, int rhs) {
        if (patch.layout[lhs][1] != patch.layout[rhs][1])
            return patch.layout[lhs][1] > patch.layout[rhs][1];
        return patch.layout[lhs][0] < patch.layout[rhs][0];
    });
    std::vector<int> columnOf(nControl, -1);
    for (int column = 0; column < nControl; column++)
        columnOf[controlOrder[column]] = column;

    // abar: one row per child, one column per control point.
    std::vector<std::vector<double>> abar(N + 12, std::vector<double>(nControl, 0.0));
    for (int row = 0; row < static_cast<int>(children.size()); row++)
    {
        const ChildPoint &child = children[row];
        if (child.isVertexPoint)
        {
            const int degree = static_cast<int>(neighbours[child.a].size());
            const double beta = loop_vertex_weight(degree);
            abar[row][columnOf[child.a]] = 1.0 - degree * beta;
            for (int neighbour : neighbours[child.a])
            {
                abar[row][columnOf[neighbour]] = beta;
            }
        }
        else
        {
            abar[row][columnOf[child.a]] += 3.0 / 8.0;
            abar[row][columnOf[child.b]] += 3.0 / 8.0;
            for (int opposite : edges.at(edge_key(child.a, child.b)).opposite)
            {
                abar[row][columnOf[opposite]] += 1.0 / 8.0;
            }
        }
    }

    // ---- selection matrices -------------------------------------------------
    //
    // Each subdivision step peels three regular children off the patch and
    // hands the extraordinary corner to a smaller copy of itself. Their
    // one-rings are recovered by running the same opposite-node walk the mesh
    // uses in set_one_ring_vertices_sorted(), over the subdivided topology --
    // rather than transcribing index lists, which is what stopped the original
    // port at a single valence.

    // Child index lookups: vertex points exist only for d4, d7 and d8; edge
    // points only for interior edges.
    std::vector<int> childOfVertex(nControl, -1);
    std::map<EdgeKey, int> childOfEdge;
    for (int row = 0; row < static_cast<int>(children.size()); row++)
    {
        if (children[row].isVertexPoint)
            childOfVertex[children[row].a] = row;
        else
            childOfEdge[edge_key(children[row].a, children[row].b)] = row;
    }
    auto edgeChild = [&](int a, int b) {
        const auto found = childOfEdge.find(edge_key(a, b));
        return (found == childOfEdge.end()) ? -1 : found->second;
    };

    // Subdivided faces: each original triangle splits into four. A sub-face is
    // kept only when all three of its corners exist as children, which drops
    // the ones hanging off the boundary of the local patch.
    std::vector<std::vector<int>> subNeighbours(children.size());
    auto link = [&](int a, int b) {
        if (std::find(subNeighbours[a].begin(), subNeighbours[a].end(), b) ==
            subNeighbours[a].end())
            subNeighbours[a].push_back(b);
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
                continue;
            link(sub[0], sub[1]);
            link(sub[1], sub[0]);
            link(sub[1], sub[2]);
            link(sub[2], sub[1]);
            link(sub[2], sub[0]);
            link(sub[0], sub[2]);
        }
    }

    // The common neighbour of node1 and node2 that is not node3 -- the same
    // rule as Mesh::find_opposite_node_index().
    auto opposite = [&](int node1, int node2, int node3) {
        for (int candidate : subNeighbours[node1])
        {
            if (candidate == node3)
                continue;
            if (std::find(subNeighbours[node2].begin(), subNeighbours[node2].end(), candidate) !=
                subNeighbours[node2].end())
                return candidate;
        }
        return -1;
    };

    // The one-ring of a child triangle, in the same d1..d12 order the mesh
    // builds. `regular` drops d1, which a valence-5 ring does not have.
    auto oneRing = [&](int d4c, int d7c, int d8c, bool regular) {
        const int d3c = opposite(d4c, d7c, d8c);
        const int d11c = opposite(d7c, d8c, d4c);
        const int d5c = opposite(d4c, d8c, d7c);
        const int d2c = opposite(d4c, d5c, d8c);
        const int d6c = opposite(d3c, d7c, d4c);
        const int d9c = opposite(d8c, d5c, d4c);
        const int d10c = opposite(d7c, d11c, d8c);
        const int d12c = opposite(d8c, d11c, d7c);
        std::vector<int> ring;
        if (regular)
        {
            ring.push_back(opposite(d3c, d4c, d7c)); // d1
        }
        for (int index : {d2c, d3c, d4c, d5c, d6c, d7c, d8c, d9c, d10c, d11c, d12c})
            ring.push_back(index);
        return ring;
    };

    const int vd4 = childOfVertex[0];
    const int vd7 = childOfVertex[1];
    const int vd8 = childOfVertex[2];
    const int e47 = edgeChild(0, 1);
    const int e48 = edgeChild(0, 2);
    const int e78 = edgeChild(1, 2);

    auto selection = [&](const std::vector<int> &ring) {
        std::vector<std::vector<double>> matrix(ring.size(),
                                                std::vector<double>(children.size(), 0.0));
        for (std::size_t i = 0; i < ring.size(); i++)
        {
            if (ring[i] < 0)
            {
                throw std::logic_error(
                    "[build_subdivision_matrices] the one-ring walk failed on a child patch at "
                    "valence " + std::to_string(N));
            }
            matrix[i][ring[i]] = 1.0;
        }
        return Matrix(matrix);
    };

    SubdivisionMatrices result;
    result.valence = N;
    result.abar = Matrix(abar);
    result.p1 = selection(oneRing(e48, e78, vd8, true));  // corner child at d8
    result.p2 = selection(oneRing(e47, e78, e48, true));  // centre child
    result.p3 = selection(oneRing(e47, vd7, e78, true));  // corner child at d7

    // The irregular child is a canonical patch of the same valence, so its
    // one-ring has to come out in the same order as its parent's -- that
    // self-similarity is what lets (p4 * abar)^d compose. Building it through
    // the fixed d1..d12 walk would hardcode the valence-5 shape: too few
    // points at N=7 and 8, and a repeated point at N=4.
    std::vector<int> childRing(nControl, -1);
    childRing[0] = vd4;                       // the extraordinary corner survives
    for (int i = 0; i < N; i++)
    {
        childRing[1 + i] = edgeChild(0, 1 + i); // its fan, in the parent's fan order
    }
    const int c4 = childRing[0];
    const int c7 = childRing[1];
    const int c8 = childRing[2];
    const int c5 = childRing[3];
    const int c3 = childRing[N];
    childRing[N + 4] = opposite(c7, c8, c4);                  // d11
    childRing[N + 1] = opposite(c3, c7, c4);                  // d6
    childRing[N + 2] = opposite(c8, c5, c4);                  // d9
    childRing[N + 3] = opposite(c7, childRing[N + 4], c8);    // d10
    childRing[N + 5] = opposite(c8, childRing[N + 4], c7);    // d12

    std::vector<int> childRingCanonical(nControl, -1);
    for (int internal = 0; internal < nControl; internal++)
    {
        childRingCanonical[columnOf[internal]] = childRing[internal];
    }
    result.p4 = selection(childRingCanonical);
    return result;
}
