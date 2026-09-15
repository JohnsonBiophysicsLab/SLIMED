/**
 * @file Vertex.hpp
 * @author Y Ying (yying7@jh.edu)
 * @brief This file defines the Vertex class which includes information
 * of coordinate and forces exerted at a vertex.
 * @date 2023-01-05
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <array>
#include <string>
#include <vector>

#include "energy_force/Force.hpp"
#include "mesh/Face.hpp"
//#include <Edge.hpp>
#include "linalg/Linear_algebra.hpp"

/**
 * @brief What a vertex is to the boundary condition.
 *
 * Under the global boundary modes (BoundaryType::Fixed, Periodic and Free)
 * the mesh decides everything from grid position and only ever writes Ghost
 * here. Under BoundaryType::Mixed the type is a property of the vertex
 * itself -- set by the flat-sheet generator from the per-axis boundary
 * parameters, or read from the mesh file -- and the whole boundary treatment
 * follows from it: which vertices are degrees of freedom, where the force on
 * an image is folded, which faces carry energy. That is what lets one mesh
 * be periodic along one axis and clamped or open along another, or a closed
 * tube with a free rim. See docs/mixed_boundary_conditions.md.
 *
 * The integer values are part of the vertex-type output format
 * (Mesh::write_vertices_csv_with_type(), read by analysis/fluidity.py) and
 * must not move.
 */
enum class VertexType
{
    /** An ordinary vertex: a degree of freedom wherever it sits, on an open edge included. */
    Free = 0,

    /** Clamped in space: carries no force and never moves. */
    Fixed = 1,

    /** @deprecated Never assigned; kept so the values below do not move. */
    PeriodicBoundary = 2,

    /**
     * A periodic image of the vertex at Vertex::reflectiveVertexIndex. Not a
     * coordinate of its own: it sits at the source's position plus
     * Vertex::mirrorOffset, follows the source wherever it goes, and the force
     * it accumulates is folded onto the source. Under the global Periodic mode
     * the dynamics stamps this on the fourth-ring duplicates as it goes.
     */
    Periodic = 3,

    /** Scaffolding only: pinned, takes no part in the physics. Written by the global modes. */
    Ghost = 4,

    // The names this enum had before the per-vertex boundary mode existed.
    Real = Free,
    FixedBoundary = Fixed,
    PeriodicReflectiveBoundary = Periodic,
};

/// "free", "fixed", "periodic" or "ghost": the spelling the mesh files use.
const char *vertex_type_name(VertexType type);

/**
 * @brief Parse a vertex type as the mesh files spell it.
 *
 * Accepts the names above in any letter case, "real" and "image" as synonyms
 * of free and periodic, and the integer values of the enum.
 *
 * @return false if the text names no type; @p type is then left untouched.
 */
bool parse_vertex_type(const std::string &text, VertexType &type);

/**
 * @brief Defines a vertex in Mesh.
 */
class Vertex
{

public:
    int index;                         ///< index number of the vertex
    // Coordinates
    Matrix coord;                      ///< coordinate of the vertex
    Matrix coordPrev;                  ///< coordinate of previous step
    Matrix coordRef;                   ///< reference coordinate
    // Geometry
    std::vector<int> adjacentVertices; ///< index of adjacent vertices
    std::vector<int> adjacentFaces;    ///< inde of adjacent faces
    Matrix normVector;                 ///< normal vector to the corresponding point on limit surface
    // Layer and Boundary Condition
    int layerIndex = 0;                ///< used in multi-layer (a.k.a. leaflet) model; 0 for middle layer; +1 for upper layer
    VertexType type = VertexType::Free;///< role in the boundary condition; see VertexType
    /**
     * @brief For a Periodic image, the vertex it mirrors: its source. -1 otherwise.
     *
     * Under BoundaryType::Mixed a chain of images (an image of an image) is
     * flattened at setup so this always names a vertex that is not itself an
     * image.
     */
    int reflectiveVertexIndex = -1;
    /**
     * @brief For a Periodic image, coord - coord(source), fixed at setup.
     *
     * The lattice translation that carries the source onto this image. Kept
     * explicitly rather than re-derived from displacements so that syncing an
     * image is one assignment and cannot drift.
     */
    std::array<double, 3> mirrorOffset = {{0.0, 0.0, 0.0}};
    bool isBoundary = false;           ///< @deprecated
                                       ///< point on boundary of the membrane; determined by boundary condition
    // Forces
    Force force;                       ///< force exerted on the vertex
    Force forcePrev;                   ///< force exerted on the vertex from the previous iteration

    //Halfedge* halfedge;                ///< one half edge

    /**
     * @brief “Ghost vertices” are defined as points on the boundary of the triangular mesh
     * that only serve to provide reference when calculating limit surface on the boundary,
     * as calculating position of a point on the limit surface require the coordinates of 12
     * neighboring vertices (if regular). However, the “ghost vertices” themselves do not
     * correspond to real points on the surface.
     *
     * Under BoundaryType::Mixed this is derived rather than declared: a
     * periodic image none of whose faces carries energy, or a vertex typed
     * Ghost outright. See Mesh::apply_vertex_boundary_types().
     */
    bool isGhost = false;

    /// Whether this vertex is a periodic image of another one.
    bool is_periodic_image() const
    {
        return type == VertexType::Periodic && reflectiveVertexIndex >= 0;
    }

    /// Whether this vertex is clamped in space.
    bool is_fixed() const { return type == VertexType::Fixed; }

    /**
     * @brief Construct a new Vertex object with members declared not initialized.
     *
     */
    Vertex();

    /**
     * @brief Constructs a new Vertex object with x,y,z coordinates.
     * @param index The index of the vertex.
     * @param x The x coordinate of the vertex.
     * @param y The y coordinate of the vertex.
     * @param z The z coordinate of the vertex.
     */
    Vertex(const int index, const double x, const double y, const double z);

    /**
     * @brief Constructs a new Vertex object with Matrix of size (3,1)
     * @param index The index of the vertex
     * @param coord Matrix(3,1) coordinates of the vertex
     */
    Vertex(const int index, const Matrix &coord);

    /**
     * @brief Update the current vertex coordinate with the previous coordinate by
     * copying the values of the previous coordinate matrix to the current
     * coordinate.
     *
     */
    void update_coord_with_prev_coord();

    /**
     * @brief Approximate the normal vector of the current vertex by averaging
     * the normal vector of all surrounding faces.
     * The output will be store in member variable normVector
     *
     */
    void approximate_unit_normal_vector(std::vector<Face> &faces);
};
