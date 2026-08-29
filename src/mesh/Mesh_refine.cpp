/**
 * @file Mesh_refine.cpp
 * @brief One global Loop refinement pass over the control mesh.
 *
 * For meshes that cannot be generated with isolated extraordinary vertices.
 * Off by default, and deliberately global: Loop's rules are uniform, so
 * refining only part of a control mesh changes the limit surface everywhere
 * else too. Local refinement around adjacent extraordinary vertices is not a
 * cheaper version of this -- it is a different, wrong, thing.
 *
 * The separation guarantee is structural rather than statistical. One pass
 * splits every triangle into four:
 *
 *        v_a                 v_a
 *         /\                  /\
 *        /  \      -->    e_ca--e_ab
 *       /    \             / \  / \
 *   v_c ------ v_b     v_c --e_bc-- v_b
 *
 * Every new vertex sits on an old edge and has valence 6. Every old vertex
 * keeps its valence but is now adjacent only to new vertices. So of the four
 * children, three carry exactly one old corner and the middle one carries
 * none -- no face can have two extraordinary corners, whatever the input was.
 *
 * The cost is real and is why this is opt-in: roughly four times the vertices,
 * and the control mesh IS the dynamical degrees of freedom, so refining
 * rebaselines the run.
 *
 * @note WP7 of docs/irregular_patch_valence_4_to_8_plan.md.
 */

#include "mesh/Mesh.hpp"

#include <algorithm>
#include <cstdint>
#include <map>
#include <utility>

namespace
{
using EdgeKey = std::pair<int, int>;

EdgeKey edge_key(int a, int b)
{
    return (a < b) ? EdgeKey{a, b} : EdgeKey{b, a};
}

/// Incident faces' opposite corners, and the index of this edge's new vertex.
struct EdgeRecord
{
    std::vector<int> opposite;
    int childIndex = -1;
};

Matrix scaled(const Matrix &point, double factor)
{
    Matrix result(3, 1);
    for (int k = 0; k < 3; k++)
        result.set(k, 0, point(k, 0) * factor);
    return result;
}

void add_scaled(Matrix &target, const Matrix &point, double factor)
{
    for (int k = 0; k < 3; k++)
        target.set(k, 0, target(k, 0) + point(k, 0) * factor);
}
} // namespace

void Mesh::refine_loop_once()
{
    const int nOldVertices = static_cast<int>(vertices.size());

    // Topology, read straight off the face list: this runs before the adjacency
    // members are populated, so it cannot use them.
    std::map<EdgeKey, EdgeRecord> edges;
    std::vector<std::vector<int>> neighbours(nOldVertices);
    std::vector<int> incidentFaces(nOldVertices, 0);
    for (const Face &face : faces)
    {
        if (face.adjacentVertices.size() != 3)
            continue;
        for (int k = 0; k < 3; k++)
        {
            const int a = face.adjacentVertices[k];
            const int b = face.adjacentVertices[(k + 1) % 3];
            const int opposite = face.adjacentVertices[(k + 2) % 3];
            edges[edge_key(a, b)].opposite.push_back(opposite);
            incidentFaces[a]++;
            if (std::find(neighbours[a].begin(), neighbours[a].end(), b) == neighbours[a].end())
                neighbours[a].push_back(b);
            if (std::find(neighbours[b].begin(), neighbours[b].end(), a) == neighbours[b].end())
                neighbours[b].push_back(a);
        }
    }
    // Each face counted each of its corners once per edge walk, i.e. twice.
    for (int v = 0; v < nOldVertices; v++)
        incidentFaces[v] /= 2;

    auto is_interior = [&](int v) {
        return incidentFaces[v] == static_cast<int>(neighbours[v].size());
    };

    std::vector<Vertex> refined;
    refined.reserve(nOldVertices + edges.size());

    // Even points: reposition every old vertex. Its valence is unchanged, so
    // an extraordinary vertex stays extraordinary -- the guarantee is about
    // adjacency, not about removing extraordinary vertices.
    for (int v = 0; v < nOldVertices; v++)
    {
        Matrix moved(3, 1);
        moved.set_all(0.0);
        if (is_interior(v))
        {
            const int degree = static_cast<int>(neighbours[v].size());
            const double beta = 3.0 / (8.0 * degree);
            add_scaled(moved, vertices[v].coord, 1.0 - degree * beta);
            for (int neighbour : neighbours[v])
                add_scaled(moved, vertices[neighbour].coord, beta);
        }
        else
        {
            // Loop's boundary rule: the vertex keeps 3/4 and gives 1/8 to each
            // neighbour along the boundary curve, so the boundary refines as a
            // cubic B-spline independent of the interior.
            std::vector<int> boundaryNeighbours;
            for (int neighbour : neighbours[v])
            {
                const EdgeRecord &record = edges.at(edge_key(v, neighbour));
                if (record.opposite.size() == 1)
                    boundaryNeighbours.push_back(neighbour);
            }
            if (boundaryNeighbours.size() == 2)
            {
                add_scaled(moved, vertices[v].coord, 0.75);
                add_scaled(moved, vertices[boundaryNeighbours[0]].coord, 0.125);
                add_scaled(moved, vertices[boundaryNeighbours[1]].coord, 0.125);
            }
            else
            {
                // A corner, or a non-manifold boundary vertex: pin it.
                moved = vertices[v].coord;
            }
        }

        Vertex vertex = vertices[v];
        vertex.index = static_cast<int>(refined.size());
        vertex.coord = moved;
        refined.push_back(vertex);
    }

    // Odd points: one new vertex per edge, each of valence 6 by construction.
    for (auto &entry : edges)
    {
        const int a = entry.first.first;
        const int b = entry.first.second;
        Matrix point(3, 1);
        point.set_all(0.0);
        if (entry.second.opposite.size() == 2)
        {
            add_scaled(point, vertices[a].coord, 3.0 / 8.0);
            add_scaled(point, vertices[b].coord, 3.0 / 8.0);
            add_scaled(point, vertices[entry.second.opposite[0]].coord, 1.0 / 8.0);
            add_scaled(point, vertices[entry.second.opposite[1]].coord, 1.0 / 8.0);
        }
        else
        {
            // Boundary edge: the midpoint, matching the boundary vertex rule.
            point = scaled(vertices[a].coord + vertices[b].coord, 0.5);
        }

        Vertex vertex;
        vertex.index = static_cast<int>(refined.size());
        vertex.coord = point;
        // A new vertex inherits ghost status only when the whole edge is ghost,
        // so refinement cannot quietly widen the physical region.
        vertex.isGhost = vertices[a].isGhost && vertices[b].isGhost;
        entry.second.childIndex = vertex.index;
        refined.push_back(vertex);
    }

    // Four children per face, winding preserved.
    std::vector<Face> refinedFaces;
    refinedFaces.reserve(faces.size() * 4);
    auto emit = [&](int a, int b, int c, const Face &parent) {
        Face child;
        child.index = static_cast<int>(refinedFaces.size());
        child.adjacentVertices = std::vector<int>{a, b, c};
        child.isGhost = parent.isGhost;
        child.spontCurvature = parent.spontCurvature;
        child.isInsertionPatch = parent.isInsertionPatch;
        refinedFaces.push_back(child);
    };
    for (const Face &face : faces)
    {
        if (face.adjacentVertices.size() != 3)
            continue;
        const int a = face.adjacentVertices[0];
        const int b = face.adjacentVertices[1];
        const int c = face.adjacentVertices[2];
        const int ab = edges.at(edge_key(a, b)).childIndex;
        const int bc = edges.at(edge_key(b, c)).childIndex;
        const int ca = edges.at(edge_key(c, a)).childIndex;
        emit(a, ab, ca, face);
        emit(ab, b, bc, face);
        emit(ca, bc, c, face);
        emit(ab, bc, ca, face);
    }

    vertices = std::move(refined);
    faces = std::move(refinedFaces);

    std::cout << "[Mesh::refine_loop_once] Refined to " << vertices.size() << " vertices and "
              << faces.size() << " faces; every face now carries at most one extraordinary corner."
              << std::endl;
}
