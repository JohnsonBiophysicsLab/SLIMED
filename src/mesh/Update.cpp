/**
 * @file Update.cpp
 * @author Y Ying (yying7@jh.edu)
 * @brief Contains member functions for updating mesh properties
 * @date 2023-03-27
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "mesh/Mesh.hpp"

/**
 * @brief Update the `coordPrev` member variable for each vertex in the mesh
 *
 * This function updates the `coordPrev` member variable for each vertex in the mesh. 
 * It does this by copying the current value of `coord` to `coordPrev`.
 */
void Mesh::update_previous_coord_for_vertex()
{
#pragma omp parallel for
    for (Vertex& vertex : vertices)
    {
        vertex.coordPrev = vertex.coord;
    }
}

/**
 * @brief Update the `coordRef` member variable for each vertex in the mesh
 *
 * This function updates the `coordRef` member variable for each vertex in the mesh.
 * It does this by copying the current value of `coordPrev` to `coordRef`.
 */
void Mesh::update_reference_coord_from_previous_coord()
{
#pragma omp parallel for
    for (Vertex& vertex : vertices)
    {
        vertex.coordRef = vertex.coordPrev;
    }
}

/**
 * @brief Update the `forcePrev` member variable for each vertex in the mesh
 *
 * This function updates the `forcePrev` member variable for each vertex in the mesh.
 * It does this by copying the current value of `force` to `forcePrev`.
 */
void Mesh::update_previous_force_for_vertex()
{
#pragma omp parallel for
    for (Vertex& vertex : vertices)
    {   
        vertex.forcePrev = vertex.force;
    }
}

/**
 * @brief Update the `energyPrev` member variable for each face in the mesh
 *
 * This function updates the `energyPrev` member variable for each face in the mesh.
 * It does this by copying the current value of `energy` to `energyPrev`.
 */
void Mesh::update_previous_energy_for_face()
{
#pragma omp parallel for
    for (Face& face : faces)
    {   
        face.energyPrev = face.energy;
    }
}

/**
 * @brief this function sets the force member variable of each vertex, and the energy
 * member variable of each face to their default values. This is useful to clear out
 * any residual forces or energies before computing new ones.
 * 
 */
void Mesh::clear_force_on_vertices_and_energy_on_faces(){
#pragma omp parallel for
    for (int i = 0; i < vertices.size(); i++){
        Force forcetmp;
        vertices[i].force = forcetmp;
    }
#pragma omp parallel for
    for (int i = 0; i < faces.size(); i++){
        Energy energytmp;
        faces[i].energy = energytmp;
    }
}