/**
 * @file Edge.hpp
 * @brief The undirected edge record the Monte Carlo flip move operates on.
 *
 * Until now this tree had no edge structure at all: HalfedgeMesh.cpp and the
 * old Edge.cpp were commented out wholesale, and the three places that needed
 * edges built a hash map, used it, and threw it away
 * (Mesh::set_adjacent_faces_of_faces(), validate_volume_constraint_topology(),
 * refine_loop_once()). That is fine for a mesh whose connectivity is fixed at
 * setup. An edge flip needs the opposite: a table it can look an edge up in,
 * mutate in place, and hand to the Metropolis sweep as the set to sample from.
 *
 * The record deliberately stores the two *opposite* corners alongside the two
 * incident faces, because those are exactly the flip targets: flipping the
 * edge (v[0], v[1]) replaces it with (opposite[0], opposite[1]).
 *
 * @see docs/edge_flip_plan.md section 3.1
 */

#pragma once

#include <cstdint>

/**
 * @brief One undirected edge of the control mesh.
 *
 * Endpoints are kept sorted so that the record is independent of which
 * incident face is looked at, and so that undirected_edge_key() below is a
 * faithful index into the table.
 *
 * The two slots are parallel: `opposite[k]` is the third corner of `face[k]`,
 * the one not on this edge. On a mesh boundary only slot 0 is filled and
 * `face[1]` stays negative.
 */
struct MeshEdge
{
    int index = -1;             ///< This edge's own position in Mesh::edges
    int v[2] = {-1, -1};        ///< Endpoints, with v[0] < v[1]
    int face[2] = {-1, -1};     ///< Incident faces; face[1] < 0 on a boundary edge
    int opposite[2] = {-1, -1}; ///< Third corner of face[k]; together, the flip targets

    /**
     * @brief Whether a flip of this edge is structurally permitted.
     *
     * The part of the admission test that depends only on connectivity and on
     * which vertices are frozen, so it is settled when the table is built and
     * re-settled when a flip changes the neighbourhood. The part that depends
     * on the valences a flip would produce is checked per attempt, because a
     * neighbouring accepted flip changes it.
     *
     * @see Mesh::edge_flip_is_admissible()
     */
    bool flippable = false;

    /// A boundary edge carries one face and can never be flipped.
    bool is_boundary() const { return face[0] < 0 || face[1] < 0; }

    /// The endpoint that is not @p vertex, or -1 if @p vertex is not an endpoint.
    int other_endpoint(int vertex) const
    {
        if (v[0] == vertex) return v[1];
        if (v[1] == vertex) return v[0];
        return -1;
    }
};

/**
 * @brief Pack an edge into one key with its endpoints in ascending order, so
 * the two faces meeting on an edge land in the same bucket however each of
 * them winds it.
 *
 * Mirrors directed_edge_key() in Mesh_setup_boundary_condition.cpp, which
 * deliberately does not sort because it is counting orientations. This one
 * used to live in an anonymous namespace in Mesh_setup_geometry.cpp; it is
 * here now because the edge table, the face-adjacency pass and the flip all
 * need the same key and must not drift apart.
 */
inline std::uint64_t undirected_edge_key(int nodeA, int nodeB)
{
    const std::uint32_t low = static_cast<std::uint32_t>(nodeA < nodeB ? nodeA : nodeB);
    const std::uint32_t high = static_cast<std::uint32_t>(nodeA < nodeB ? nodeB : nodeA);
    return (static_cast<std::uint64_t>(low) << 32) | high;
}
