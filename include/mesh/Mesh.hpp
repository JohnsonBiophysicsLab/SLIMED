/**
 * @file Mesh.hpp
 * @author Y Ying
 * @brief This file defines mesh class containing a
 * vector of vertices and a vector of faces that are
 * required to build up a triangular mesh that represents
 * a limit surface. It provides functions to calculate
 * limit surface from control mesh and control mesh from
 * limit surface.
 *
 * The model class combines mesh and parameter to provide
 * a full environment for simulation.
 * @date 2023-02-01
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
#include <string>
#include <stdexcept>
#include <omp.h>
#include <algorithm>
// model setup
//#include "Edge.hpp"
#include "mesh/Face.hpp"
#include "mesh/Vertex.hpp"
#include "energy_force/Energy.hpp"
#include "energy_force/Force.hpp"
#include "mesh/Gauss_quadrature.hpp"
// matrix math
#include "linalg/Linear_algebra.hpp"
// parameters
#include "Parameters.hpp"


using namespace std;

/**
 * @brief A class representing a triangular mesh that defines a
 * limit surface.
 *
 */
class Mesh
{
public:
    std::vector<Vertex> vertices; ///< Vector to store all vertices in the mesh
    std::vector<Face> faces;      ///< Vector to store all faces in the mesh
    Param& param;                  ///< Object of the Param class containing all necessary parameters for building the Mesh object
    // Scaffolding points, see Scaffolding_points.cpp
    Matrix centerScaffoldingSphere; ///< Center of the scaffolding cap sphere
    Matrix forceTotalOnScaffolding; ///< Total force exerted on the scaffolding lattice
    Matrix scaffoldingMovementVector; ///< Vector representing the movement of scaffolding over the course of simulation

    // New members... for halfedge mesh
    //std::vector<Edge> edges; ///< Vector to store all edges in the mesh
    //std::vector<Halfedge> halfedges; ///< Vector to store all halfedges in the mesh


    /**
     * @brief Construct a new Mesh object with parameters. Initialize
     * vertices and faces with other functions like setVerticesFlat.
     *
     * @param param
     */
    Mesh(Param &srcParam);

    /**
     * 
     * @deprecated This function is deprecated, initialize with Mesh(Param &srcParam)
     * and setup with setup_from_vertices_faces instead.
     * @brief Construct a new mesh object with given vertices, faces,
     * and parameters.
     */
    Mesh(const std::vector<Vertex> &srcVertices, const std::vector<Face> &srcFaces, Param &srcParam)
    __attribute__((deprecated("Initialize with Mesh(Param &srcParam) and setup with setup_from_vertices_faces instead.")));


    /**
     * @brief Initialize halfedges based on the mesh's vertices and faces.
     */
    //void initializeHalfedges();

    /**
     * @brief Helper function to create a new halfedge.
     * 
     * @param vertexIndex Index of the vertex associated with the halfedge.
     * @param faceIndex Index of the face associated with the halfedge.
     * @return Halfedge* Pointer to the created halfedge.
     */
    //Halfedge* createHalfedge(int vertexIndex, int faceIndex);

    // setup

    /**
     * @brief Initialize flat membrane with properities specified in
     * the param} member variable. Equivalent to call:
     * (1) this->set_axes_division_flat
     * (2) this->set_vertices_faces_flat
     * (3) this->set_adjacent_faces_of_vertices_sorted
     * (4) this->determine_ghost_vertices_faces
     *
     */
    void setup_flat();

    /**
     * @brief Initialize membrane with arbitrary vertices and faces imported
     * from files. Sets up the membrane from data and then call:
     * (1) this->set_adjacent_faces_of_vertices_sorted
     * (2) this->determine_ghost_vertices_faces
     * 
     * @param verticesData 2D double vector containing the coordinates of vertices
     * @param facesData 2D int vector containing the indices of vertices on each face
     */
    void setup_from_vertices_faces(const std::vector<std::vector<double>>& verticesData, 
                                   const std::vector<std::vector<int>>& facesData);

    /**
     * @brief Divide x,y axis to nx*dx (number of faces times length of
     * each face) and ny*dy based on X, Y side length of the mesh and
     * side length of faces in parameter. Note:
     *
     * (1) y-axis is in a zig-zag shape and in perpendicular direction,
     * h(y) = sqrt(3)/2 * h(x) assuming equilateral triangles
     * (2) number of vertices along axes = number of faces (edges) + 1
     *
     */
    void set_axes_division_flat();

    /**
     * @brief Set vertices and faces according to mesh
     * and face side lengths in param for flat mesh. This also
     * sets adjacent vertices of faces.
     */
    void set_vertices_faces_flat();

    /**
     * @brief Set adjacentFaces properties based on the current
     * vertices and mesh. Sort the adjacent vertices so that the adjacent
     * vertices property of vertices follow the counterclockwise order
     * and therefore the the adjacent faces with index number difference
     * of one are adjacent to each other. This sorting streamlines the
     * shapefunction calculation.
     *
     */
    void set_adjacent_faces_of_vertices_sorted();

    /**
     * @brief Return true if two faces share edge.
     * 
     * @note Will also return true if face1 and face2 are the same face.
     * 
     * @note This function overloads faces_share_edge(const Face& face1, const Face& face2, )
     * 
     * @param face1 
     * @param face2 
     * @return true 
     * @return false 
     */
    bool faces_share_edge(const Face& face1, const Face& face2);

    /**
     * @brief Return true if two faces share edge. Store the vector of common elements
     * in input commonElements
     * 
     * @note Will also return true if face1 and face2 are the same face.
     * 
     * @note Used in set_adjacent_faces_of_faces()
     * @note see faces_share_edge(const Face& face1, const Face& face2, )
     * 
     * @param face1 
     * @param face2 
     * @param commonElements 
     * @return true 
     * @return false 
     */
    bool faces_share_edge(const Face& face1, const Face& face2, std::vector<int>& commonElements);

    /**
     * @brief Set adjacentFaces properties of faces based on the current
     * geometry of mesh.
     * 
     * This function iterates over the faces of the mesh and populates the
     * adjacentFaces property of each face by finding neighboring faces that
     * share an edge.
     */
    void set_adjacent_faces_of_faces();

    /**
     * @brief Set adjacentVertices of vertices based on current
     * mesh.
     * @todo Sorting to be implemented.
     *
     */
    void set_adjacent_vertices_of_vertices_sorted();

    /**
     * @brief Find the vertex that is adjacent to node1 and node2
     * but not node3. Particularly, when the three nodes are vertices of
     * a triangle, then the function returns the the index of node4 that
     * forms a parallegram with 1->3->2->4.
     *
     * @param node1
     * @param node2
     * @param node3
     * @return int vertex index
     */
    int find_opposite_node_index(const int &node1, const int &node2, const int &node3);

    /**
     * @brief find out the one-ring vertices aound face_i. It should be 12 for the flat surface because we set it up only with regular patch.
     *The boundary faces do not have complete one-ring, neither it will be called in the code, so no need to store their one-ring-vertex
     *
     */
    void set_one_ring_vertices_sorted();

    /**
     * @deprecated Currently the isBoundary property is not used
     * in any part of the model. This is a placeholder in case any future
     * functions need the property.
     * @brief Iterate through vertices and faces and set the isBoundary}
     * property of boundary vertices and faces to true}.
     *
     */
    void determine_boundary_vertices_faces();

    /**
     * @brief Iterate through vertices and faces and set the isGhost
     * property of ghost vertices and face to true. The number of layers
     * of ghost vertices in a flat membrane is dependent upon the boundary conditions:
     * free boundary condition has 1 layer of ghost vertices and faces while
     * periodic bounary condition has 3 layers of ghost vertices.
     *
     */
    void determine_ghost_vertices_faces();

    /**
     * @brief Sort vertices on faces so that the unit normal vector indicates
     * the orientation of the local patch of the membrane.
     * 
     * For example, if a face has vertices A->B->C, then the unit normal vector
     * is calculated as AB x BC. This follows a "half-edge" data structure:
     * if face ABC and face BCD shares edge BC, and ABC has vertices A->B->C, then
     * on BCD, the edge sequence of BC needs to be reverse and therefore BCD has
     * vertices C->B->D.
     * 
     */
    void sort_vertices_on_faces();

    /**
     * @brief Sets the isInsertionPatch flag for each face in the insertion patch.
     * This also sets with spontaneous curvature after setting the flag with
     * set_spontaneous_curvature_for_face function.
     *
     * @param insertionPatch A vector of vectors containing the indices of faces in the insertion patch
     */
    void set_insertion_patch(const vector<vector<int>> &insertionPatch);

    /**
     * @brief Sets the spontaneous curvature for each face in the mesh.
     *
     * This method sets the "spontaneous curvature" of each face in the mesh, which is a
     * parameter used to describe the shape and behavior of lipid membranes. The spontaneous
     * curvature can be different for faces inside and outside of an "insertion patch," which
     * are specified by the `insertCurv` and `spontCurv` parameters, respectively.
     *
     * The method takes two double parameters: `insertCurv` and `spontCurv`. The former specifies
     * the spontaneous curvature for faces within the insertion patch, while the latter specifies
     * the spontaneous curvature for faces outside of the insertion patch.
     *
     * The method loops over all faces in the mesh and sets the `spontCurvature` member variable
     * of each face to either `insertCurv` or `spontCurv`, depending on whether the face is part
     * of the insertion patch or not.
     * @param insertCurv The spontaneous curvature for faces in the insertion patch
     * @param spontCurv The spontaneous curvature for faces outside the insertion patch
     */
    void set_spontaneous_curvature_for_face(const double &insertCurv, const double &spontCurv);

    /**
     * @brief Calculates the area and volume of each element in the mesh.
     *
     * This method computes and sets the area and volume of each element in the mesh,
     * where an "element" refers to a face in the mesh. The method uses Gaussian quadrature
     * to calculate the area and volume of each element, with 3 points used for regular patches
     * and multiple iterations of subdivision used for irregular patches.
     *
     * The method first defines some matrices used for subdivision of the irregular patch,
     * including `M`, `M1`, `M2`, `M3`, and `M4`. It then loops over all faces in the mesh
     * in parallel using OpenMP.
     *
     * For each face, the method checks whether it is a ghost face; if so, it skips to the
     * next face. Otherwise, it initializes variables to accumulate the area and volume of the
     * element and determines the number of one-ring vertices (`nOneRingVertex`) for the face.
     *
     * If `nOneRingVertex` is equal to 12, the face corresponds to a regular patch, and the method
     * calls the `get_one_ring_vertex_matrix` function to compute a matrix representing the
     * coordinates of the one-ring vertices. It then uses Gaussian quadrature with 3 points to
     * compute the area and volume of the element.
     *
     * If `nOneRingVertex` is equal to 11, the face corresponds to an irregular patch. In this case,
     * the method again calls `get_one_ring_vertex_matrix` to compute the matrix of one-ring vertex
     * coordinates, but then performs multiple iterations of subdivision using the `M`, `M1`, `M2`,
     * `M3`, and `M4` matrices to estimate the area and volume of the element.
     *
     * Finally, the method sets the `elementArea` and `elementVolume` member variables of the face
     * to the computed area and volume, respectively.
     */
    void calculate_element_area_volume();

    /**
     * @brief This function calculates the total membrane area and volume of non-ghost faces in a given set of faces.
     * @note This function only sums up to element area and volume of all faces. It does
     * not update the element area and volume based on the current state of the mesh.
     * If you need to update and calculate the total area and volume of the mesh, please
     * call calculate_element_area_volume first before calling this function.
     *
     * @param area  a reference to a double variable to store the computed area
     * @param volume a reference to a double variable to store the computed volume
     */
    void sum_membrane_area_and_volume(double &area, double &volume);

    /**
     *
     * @brief Computes the energy and force on each vertex and face of the mesh.
     *
     * This function performs the following steps:
     *       Calculates the area and volume of each element triangle and sums up the total area and volume of the membrane.
     *       Iterates through faces and calculates the energy and force on each triangular patch.
     *       Regularizes the force and energy.
     *       Sums up the energy and force on each vertex and face.
     *       Calculates the energy due to area constraint.
     *       Calculates the energy due to volume constraint.
     *       Calculates the energy due to scaffolding.
     *
     * @return void
     */
    void Compute_Energy_And_Force();

    /**
     * @brief The purpose of this function is to calculate the energy and forces
     * for a regular element of a given mesh using the provided information about
     * the element and its one-ring neighborhood.
     *
     * @param coordOneRingVertices a constant reference to a vector of matrices
     * representing the coordinates of the vertices in the mesh's one-ring neighborhood.
     * @param spontCurv a constant double representing the spontaneous curvature.
     * @param meanCurv a non-constant double reference representing the mean curvature of the element.
     * @param normVector a non-constant matrix reference representing the normal vector of the element.
     * @param eBend a non-constant double reference representing the bending energy of the element.
     * @param fBend a non-constant matrix reference representing the bending force of the element.
     * @param fArea a non-constant matrix reference representing the area constraint force of the element.
     * @param fVolume a non-constant matrix reference representing the volume constraint force of the element.
     */
    void element_energy_force_regular(const std::vector<Matrix> &coordOneRingVertices,
                                      Face& face,
                                      const double spontCurv,
                                      double &meanCurv,
                                      Matrix &normVector,
                                      double &eBend,
                                      Matrix &fBend,
                                      Matrix &fArea,
                                      Matrix &fVolume);

    /**
     * @brief Calculates the regularization energy and force for each face
     *
     * This function calculates the regularization energy and force
     * for each face in the mesh. The regularization energy is
     * calculated based on the deformation of the face's shape and
     * area relative to an equilateral triangle with the same side length.
     * The regularization force is then calculated based on the energy
     * and applied to the vertices of the face. The function also updates
     * the deformation counts for shapes and areas.
     *
     * @return void
     *
     */
    void energy_force_regularization();

    /**
     * @brief Manages forces depending on different boundary conditions.
     *
     * This method sets the force of nodes that are part of the mesh's boundaries and ghost vertices to zero, based
     * on the type of boundary condition specified in the input parameters.
     *
     * @note See enum BoundaryType in Parameters.hpp
     *
     */
    void manage_force_for_boundary_ghost_vertex();

    /**
     * @brief Get the max force scale of vertices.force.get_total_force_magnitude()
     * 
     * @return double max force magnitude
     */
    double get_max_force_magnitude();

    /**
     * @brief Get the mean force magnitude on all vertices
     * 
     * @return mean force magnitude
     * 
     */
    double calculate_mean_force();

    // defines in update.cpp

    /**
     * @brief This function updates the coordPrev member variable for each vertex
     * in the mesh. It does this by copying the current value of coord to coordPrev
     *
     */
    void update_previous_coord_for_vertex();

    /**
     * @brief Update the reference coordinates for each vertex in the mesh
     *
     * The `coordRef` member variable for each vertex is updated with the
     * current value of `coordPrev`.
     */
    void update_reference_coord_from_previous_coord();

    /**
     * @brief Update the previous force vectors for each vertex in the mesh
     *
     * The `forcePrev` member variable for each vertex is updated with the
     * current value of `force`.
     */
    void update_previous_force_for_vertex();

    /**
     * @brief Update the previous energy values for each face in the mesh
     *
     * The `energyPrev` member variable for each face is updated with the
     * current value of `energy`.
     */
    void update_previous_energy_for_face();

    /**
     * @brief this function sets the force member variable of each vertex, and the energy
     * member variable of each face to their default values. This is useful to clear out
     * any residual forces or energies before computing new ones.
     *
     */
    void clear_force_on_vertices_and_energy_on_faces();

    // implemented in scaffolding_points.cpp

    int findClosestRcap(double membraneCapRadius);

    double interpolateHeight(double r, const std::vector<double>& r_vals, const std::vector<double>& h_vals);

    // --- 3. Load and use in case3 ---
    void loadMembraneShapeProfile(double membraneCapRadius,
                                std::vector<double>& r_vals,
                                std::vector<double>& h_vals);

    /**
     * @brief
     * This method takes in the vector of spline point and calculate the
     * average coordinates. Based on the difference between spline points
     * and mesh vertices, a difference vector is calculated and compared to
     * the target bond length. (Supposed only in Z direction). Afterwards,
     * all the mesh points are moved in the direction of the target difference
     * vector.
     *
     * @param fixDir default to true; fix move vector to (0, 0, 50) if set to
     * true; otherwise move vector is determined based on average scaffolding points
     * @return true if success
     */
    bool move_vertices_based_on_scaffolding(bool fixDir = true);

    /**
     * @brief A function that finds the center of the scaffolding sphere.
     * Currently the function assumes that the spherical cap is oriented
     * in a way that the point with biggest z-coord corresponds to the
     * x,y of the spherical center.
     * 
     * @note this function does NOT find the radius of the sphere. It takes
     * in the gag sphere radius from the parameter structure in mesh (this->param).
     * 
     * @param use_default if set to true, the function with return the defautl
     * center of scaffolding sphere (0, 0, 0)
     * @return instance of Matrix(3, 1) that represents the center of scaffolding
     * sphere.
     */
    Matrix& find_center_of_scaffolding_sphere(bool use_default);


    /**
     * @brief Approximates the radius of the scaffolding cap.
     *
     * This function approximates the radius of the scaffolding cap based on the furthest point from the center of the scaffolding sphere in terms of x and y coordinates.
     * 
     * @param use_default Flag indicating whether to use the default center of scaffolding sphere.
     *                    If set to true, the default center (0,0,0) will be used.
     *                    If set to false, the center will be determined based on the maximum z-coordinate of the scaffolding points.
     * 
     * @return The radius of the scaffolding cap.
     */
    double approximate_scaffolding_cap_radius(bool use_default);

    /**
     * @brief
     * Get a vector of indexes of vertices that are closest to
     * the scaffoldingPoints vector provided. Then set the
     * param.scaffoldingPoints_correspondingVertexIndex}
     * to represent the vertices bonded with each scaffolding point.
     *
     */
    void set_scaffolding_vertices_correspondence();

    /**
     *
     * @brief Calculates the energy and force due to the harmonic bond between scaffold points and membrane vertices.
     * @param doLocalSearch Flag indicating whether or not to perform a local search for each vertex (default=false).
     * @return double The total energy of the system due to the scaffold-membrane interactions.
     *
     * The function iterates over each scaffold point, calculates the distance between the point and its corresponding vertex,
     * and calculates the energy and force due to the harmonic bond between the two. If the energy flag is not set to include
     * harmonic bonding, the function returns 0.
     *
     * @param[in] doLocalSearch Flag to indicate whether a local search should be performed for each vertex (default=false).
     *
     *
     * @note The energy is calculated according to the formula: E = 0.5 * k * (r - l)^2, where k is the spring constant, r is the
     * distance between the scaffold point and vertex, and l is the resting length of the bond.
     * @note The force is calculated according to the formula: F = -k * (r - l) * (unit vector pointing from vertex to scaffold point).
     * @note If the distance between the scaffold point and vertex is negative, the force is multiplied by -1 in order to prevent
     * overlapping between the scaffold and membrane.
     * @note If the verbose flag is set, the function prints out information about each vertex's force due to the scaffold-membrane bond.
     */
    double calculate_scaffolding_energy_force(bool doLocalSearch);

    /**
     * @brief Propagate the scaffolding based on energy and force calculated from
     * calculate_scaffolding_energy_force(bool doLocalSearch).
     * 
     * @return true if success 
     * @return false 
     */
    bool propagate_scaffolding();

    // io

    /**
     * @brief Writes a csv file containing the adjacent vertices for each face in the mesh
     *
     * @param outfile_name The name of the output csv file
     */
    void write_faces_csv(const std::string &outfile_name);

    /**
     * @brief Writes a csv file containing the coordinates of each vertex in the mesh
     *
     * @param outfile_name The name of the output csv file
     */
    void write_vertices_csv(const std::string &outfile_name);

    /**
     * @brief Writes a csv file containing the coordinates and types of each vertex in the mesh
     *
     * @param outfile_name The name of the output csv file
     */
    void write_vertices_csv_with_type(const std::string &outfile_name);

    /**
     * @brief Overrides the operator << in ostream.
     */
    friend std::ostream &operator<<(std::ostream &stream, const Mesh &mesh);

protected:
    /**
     * @brief
     * Private member used in calculating element area volume:
     * Calculates the area and volume at a Gauss quadrature point for a given set of shape functions and dots.
     *
     * @param dots The matrix of dots representing the coordinates of the quadrature point.
     * @param gaussQuadratureCoeff The matrix of Gauss quadrature coefficients.
     * @param area A reference to a double variable storing the accumulated area.
     * @param volume A reference to a double variable storing the accumulated volume.
     */
    void enumerate_gauss_quadrature_point_area_volume(const Matrix &dots,
                                                      double &area,
                                                      double &volume);

    /**
     * @brief
     * Private member used in calculating element area volume:
     * Computes a matrix containing the coordinates of the one-ring vertices for the input face.
     *
     */
    Matrix get_one_ring_vertex_matrix(const Face &face);

    /**
     * @brief
     * Private member used in calculating correspondence between scaffolding points and vertices
     * Calculate the squared distance between two points denoted by (3,1) Matrix
     * and Vertex respectively.
     * @return the squared distance
     *
     */
    double get_squared_distance_sp_and_v(const Matrix &scaffoldingPoint, const Vertex &vertex);
};