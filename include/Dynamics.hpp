/**
 * @file Dynamics.hpp
 * @author Y Ying
 * @brief This file defines dynamic mesh class that
 * inherits mesh with support of time-resolved dynamics.
 *
 * @date 2023-04-17
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <math.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <random>
#include <algorithm>
// model setup
#include "mesh/Mesh.hpp"
#include "mesh/Face.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Gauss_quadrature.hpp"
#include "energy_force/Energy.hpp"
#include "energy_force/Force.hpp"
#include "model/Model.hpp"

// matrix math
#include "linalg/Linear_algebra.hpp"
// parameters
#include "Parameters.hpp"

/**
 * @brief A class representing a dynamic mesh.
 *
 * This class inherits from the Mesh class and adds functionality to simulate the mesh dynamically.
 */
class DynamicMesh : public Mesh
{
public:
    Matrix mesh2surface;
    Matrix surface2mesh;
    Matrix matMesh;
    Matrix matSurface;

    /**
     * @brief Construct a new Dynamic Mesh object. Calls the parent
     * constructor Mesh::Mesh(srcParam).
     *
     * @param srcParam
     */
    DynamicMesh(Param &srcParam);

    /**
     * @brief Overloads the setup_flat in superclass Mesh. Add
     * functionality to assign mesh2surface and surface2mesh.
     *
     */
    void setup_flat();

    /**
     * @brief Assign the mesh2surface member.
     *
     * Assign values to mesh2surface matrix that convertes mesh to surface point
     * matrix.
     *
     */
    void assign_mesh2surface();

    /**
     * @brief Updates the vertices matrix using the values from the vertices vector.
     *
     */
    void update_vertices_mat_with_vector();

    /**
     * @brief Updates the vertices vector using the values from the vertices matrix.
     *
     */
    void update_vertices_vector_with_mat();

    /**
     * @brief Post process ghost vertices in case of periodic boundary condition
     *
     * @note This function does not have check for PBC - it assumes the model
     * is in PBC!
     *
     */
    void postprocess_ghost_periodic();

protected:
    /**
     * @brief calculate real point relative to the given ghost point (index(real) - index(given)) in periodic boundary condition
     * returns 0 if real point (not on 4th ring) given
     * an arbitrary real point is chosen if given 4th
     *
     * @param i
     * @param n
     * @param m
     * @return int
     */
    int get_relative_pt_periodic(int i, int n, int m);
};

class DynamicModel : public Model
{
public:
    DynamicMesh &mesh; /**< Hides Model::mesh; The mesh object. */

    double randScaleConst;
    double forceScaleConst;
    unsigned int randomSeed;
    std::mt19937 gen;
    std::normal_distribution<> normal_dist;

    /**
     * @brief Constructs a new Model object.
     *
     * @param mesh_ The mesh object to be encapsulated.
     * @param record_ The record object to be encapsulated.
     */
    DynamicModel(DynamicMesh &mesh_, Record &record_);

    /**
     * @brief This code iterates over all vertices in a mesh and displaces
     * them randomly based on a force term and a random term. The displacement
     * is added to the original vertex position.
     *
     * The boundary conditions are checked for each vertex before
     * displacement, and if the vertex is a ghost or boundary,
     * and not part of a periodic boundary condition, then the displacement is set to 0.
     *
     * The displacement amount, along with the vertex index,
     * whether it is a boundary or ghost vertex, and the dimension
     * being displaced (x, y, z) is printed to the console using std::cout.
     */
    void next_step();
};

/*
//boundary condition
int get_relatve_pt_periodic(int i, int n, int m);
void postprocess_ghost_periodic(Param& param, gsl_matrix *verticesOnMesh, std::vector<Vertex>& vertex);
*/