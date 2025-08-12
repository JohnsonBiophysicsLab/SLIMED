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

#include <vector>

#include "energy_force/Force.hpp"
#include "mesh/Face.hpp"
//#include <Edge.hpp>
#include "linalg/Linear_algebra.hpp"

/**
 * @brief The type of vertex.
 */
enum class VertexType
{
    /** Real vertex that is not on the boundary. */
    Real,

    /** Real vertex that is on the fixed boundary. */
    FixedBoundary,

    /** Real vertex that is on the periodic boundary. */
    PeriodicBoundary,

    /** Real vertex that is reflective of the opposite vertex on the periodic boundary. */
    PeriodicReflectiveBoundary,

    /** Ghost vertex that is reflective of the real or periodic boundary vertex*/
    Ghost
};

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
    VertexType type = VertexType::Real;///< used in boundary condition
    int reflectiveVertexIndex = -1;    ///< reflective vertex used in periodic boundary condition
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
     */
    bool isGhost = false;

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

