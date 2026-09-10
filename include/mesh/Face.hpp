/**
 * @file Face.hpp
 * @author Y Ying
 * @brief This file defines the face class, which is the
 * triangular element defined by vertices.
 * @date 2023-01-31
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <vector>

#include "mesh/Patch_kind.hpp"
#include "energy_force/Force.hpp"
#include "energy_force/Energy.hpp"

/**
 * @brief A triangular mesh face element used in the membrane model.
 *
 * This class represents a triangular face element in a mesh used in a
 * membrane model. It contains information about the vertices and faces
 * adjacent to this face, as well as other geometric properties like
 * curvature and area.
 */
class Face
{
public:
    int index = 0;                              ///< Index of face that is used in adjacentFaces
    int layerIndex = 0;                         ///< Default to 0 for single layer membrane model; +1 for upper and -1 for lower layer
    bool isBoundary = false;                    ///< Face on boundary of the membrane; determined by boundary condition
    bool isGhost = false;                       ///< Ghost face outside boundary of the membrane; determined by boundary condition
    bool isInsertionPatch = false;              ///< True if there is an insertion on this face
    std::vector<int> adjacentVertices{0, 0, 0}; ///< The indices of 3 vertices of the triangle
    std::vector<int> oneRingVertices;           ///< AdjacentVertices + the vertices that are 1 edge away from the face

    /**
     * @brief What kind of subdivision patch this face carries.
     *
     * The evaluator used to infer this from the width of oneRingVertices
     * (`valence = size - 6`). That stops working once a face can have more
     * than one extraordinary corner: a 5/5/7 face and a valence-5 face both
     * have an 11-point control net. Set by Mesh::build_one_ring_for_face().
     */
    PatchKind patchKind = PatchKind::Boundary;

    /// Corner valences, in this face's own winding order. Only meaningful when
    /// patchKind names a patch.
    int patchValence[3] = {0, 0, 0};

    /**
     * @brief Index into Mesh::multiPatchTable, or -1.
     *
     * Resolved once when the face is classified, so the force loop never does
     * a map lookup. Only set when patchKind is MultiExtraordinary.
     */
    int patchEntry = -1;
    /// The faces this one meets along an edge, in the winding order of the
    /// edges they are met on. Three slots; a face on the mesh boundary fills
    /// fewer and the rest stay 0. Set by Mesh::set_adjacent_faces_of_faces().
    std::vector<int> adjacentFaces;
    double spontCurvature = 0.0;                ///< Spontaneous curvature of this face
    double meanCurvature = 0.0;                 ///< Mean curvature of this face
    Matrix normVector;                          ///< Outward pointing normal vector of this face element
    double elementArea = 0.0;                   ///< Local area of this face.
    double elementVolume = 0.0;                 ///< Local volume of this face.
    Energy energyPrev;                          ///< Energy of this face in previous iteration
    Energy energy;                              ///< Current energy of this face
    //Halfedge* halfedge;                         ///< One half edge

    /**
     * @brief Generate a new Face object with all member variable
     * set to default
     */
    Face();
};