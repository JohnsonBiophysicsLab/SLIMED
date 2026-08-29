#include "mesh/Mesh.hpp"

namespace
{
/**
 * Quadrature factor for the signed volume enclosed by the limit surface.
 *
 * The divergence theorem with F = x gives
 *
 *     V = (1/3) * closed_integral( x . n dA )
 *
 * and the Gauss rule below sums over the reference triangle with weights that
 * sum to one, so each sample carries the reference triangle's own area of 1/2.
 * The product (1/3) * (1/2) = 1/6 is this factor.
 *
 * It pairs with the FULL three-component integrand dot(x, a_1 x a_2). A
 * single-axis integrand is also a valid volume form, but it needs 1/2, not
 * 1/6 -- combining the two lands a factor of three low. See
 * docs/volume_functional_split.md.
 */
constexpr double kSignedVolumeQuadratureFactor = 1.0 / 6.0;
} // namespace

Mesh::Mesh(Param &srcParam) : param(srcParam)
{
    // Calculate element triangle area
    param.elementTriangleArea0 = sqrt(3.0) / 4.0 * param.lFace * param.lFace;

    // Compute Gaussian quadrature weights and coordinates
    get_gauss_quadrature_weight_VWU(param.gaussQuadratureN, param.VWU, param.gaussQuadratureCoeff);

    // Compute shape functions for each element and store them in `param`
    //param.shapeFunctions = std::vector<Matrix>(param.VWU.nrow());
    get_shapefunction_vector(param.VWU, param.shapeFunctions);

    // Collapse the irregular-patch recursion into limit-surface rows. Depends
    // only on the valence and the quadrature rule, so it is built once here
    // and shared, immutable, for the life of the mesh.
    irregularRows.build(param.shapeFunctions);

    // initialze scaffolding points matrices
    centerScaffoldingSphere = mat_calloc(3, 1); ///< Center of the scaffolding cap sphere
    forceTotalOnScaffolding = mat_calloc(3, 1); ///< Total force exerted on the scaffolding lattice
    scaffoldingMovementVector = mat_calloc(3, 1); ///< Vector representing the movement of scaffolding over the course of simulation
}

Mesh::Mesh(const std::vector<Vertex> &srcVertices,
           const std::vector<Face> &srcFaces,
           Param &srcParam) : vertices(srcVertices), faces(srcFaces), param(srcParam)
{
    // Calculate element triangle area
    param.elementTriangleArea0 = sqrt(3.0) / 4.0 * param.lFace * param.lFace;

    // Compute Gaussian quadrature weights and coordinates
    get_gauss_quadrature_weight_VWU(param.gaussQuadratureN, param.VWU, param.gaussQuadratureCoeff);

    // Compute shape functions for each element and store them in `param`
    get_shapefunction_vector(param.VWU, param.shapeFunctions);

    // Collapse the irregular-patch recursion into limit-surface rows. Depends
    // only on the valence and the quadrature rule, so it is built once here
    // and shared, immutable, for the life of the mesh.
    irregularRows.build(param.shapeFunctions);

    // initialze scaffolding points matrices
    centerScaffoldingSphere = mat_calloc(3, 1); ///< Center of the scaffolding cap sphere
    forceTotalOnScaffolding = mat_calloc(3, 1); ///< Total force exerted on the scaffolding lattice
    scaffoldingMovementVector = mat_calloc(3, 1); ///< Vector representing the movement of scaffolding over the course of simulation

}

void Mesh::setup_from_vertices_faces(const std::vector<std::vector<double>>& verticesData, 
                                   const std::vector<std::vector<int>>& facesData)
{
    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::setup_from_vertices_faces] Setting up membrane from vertices and faces data." << std::endl;
    }

    // identify if the data contains type and reflective point information
	bool simpleVerticesData = true;
	if (verticesData[0].size() > 3){
		simpleVerticesData = false;
	}

    // step 1. Initialize vertices and faces
    int nVertices = verticesData.size();                       // number of vertices
    vertices = vector<Vertex>(nVertices);                    // declare local vertices list

    int nFaces = facesData.size();                             // number of faces
    faces = vector<Face>(nFaces);

    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::set_vertices_faces_flat] nVertices = " << nVertices
                  << ", nFaces = " << nFaces << std::endl;
    }
    
    // step 2. Copy values from 2D vector to object member variables
#pragma omp parallel for
    for (int iVertex = 0; iVertex < nVertices; iVertex++){
        vertices[iVertex].index = iVertex; // set vertex index
        vertices[iVertex].coord.set(0, 0, verticesData[iVertex][0]); //set vertex coord
        vertices[iVertex].coord.set(1, 0, verticesData[iVertex][1]);
        vertices[iVertex].coord.set(2, 0, verticesData[iVertex][2]);
    }

#pragma omp parallel for
    for (int iFace = 0; iFace < nFaces; iFace++){
        faces[iFace].index = iFace; // assign face index
        // setup adjacent vertex for face
        faces[iFace].adjacentVertices = std::vector<int>(3);
        faces[iFace].adjacentVertices[0] = facesData[iFace][0]; // assign adjacent vertices
        faces[iFace].adjacentVertices[1] = facesData[iFace][1];
        faces[iFace].adjacentVertices[2] = facesData[iFace][2];
    }

    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::setup_from_vertices_faces] Assigned vertices and faces in the member of current mesh object." << std::endl;
    }

    // Optional global refinement, before anything derives from the topology.
    // It replaces both vertices and faces, so every adjacency below is built
    // from the refined mesh.
    if (param.isPreRefinementEnabled)
    {
        refine_loop_once();
    }

    // step 3. Link neighboring geometric components
    set_adjacent_faces_of_vertices_sorted();
    set_adjacent_vertices_of_vertices_sorted();
    determine_ghost_vertices_faces();
    set_one_ring_vertices_sorted();
    validate_volume_constraint_topology();

    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::setup_flat] Finished setting up flat membrane." << std::endl;
    }
}

void Mesh::set_insertion_patch(const vector<vector<int>> &insertionPatch)
{
    for (int i = 0; i < insertionPatch.size(); i++)
    {
        for (int j = 0; j < insertionPatch[0].size(); j++)
        {
            faces[insertionPatch[i][j]].isInsertionPatch = true;
        }
    }
    set_spontaneous_curvature_for_face(param.insertCurv, param.spontCurv);
}

void Mesh::set_spontaneous_curvature_for_face(const double &insertCurv, const double &spontCurv)
{
    // for all faces; set spontaneous curvature to insertion curvature if
    // they are insertion patch; otherwise set to global spontaneous curvature
#pragma omp parallel for
    for (int iFace = 0; iFace < faces.size(); iFace++)
    {
        Face& face = faces[iFace];
        // changed to isGhost - Y Ying
        if (face.isGhost)
        {
            continue;
        }

        if (face.isInsertionPatch)
        {
            face.spontCurvature = insertCurv;
        }
        else
        {
            face.spontCurvature = spontCurv;
        }
    }
}

/**
 *
 * Private member used in calculated element area volume:
 * Calculate and accumulate the area and volume at a Gauss quadrature point
 * for a given set of shape functions and matOneRingVertex.
 *
 * This is referenced constantly so
 */
void Mesh::enumerate_gauss_quadrature_point_area_volume(
    const std::vector<Matrix> &sampleRows,
    const Matrix &matOneRingVertex,
    double &area,
    double &volume)
{
#pragma omp parallel for reduction(+ \
                                   : area, volume)
    for (int j = 0; j < static_cast<int>(sampleRows.size()); j++)
    {
        /*
        Backup alternatives:
        Performance anaylsis shows version 3 is faster (although by <5%)
        Tested and compiled with gcc on SAM
        - Y Ying
        //VERSION 1
        Matrix x = sf.get_row(0) * matOneRingVertex;
        Matrix a_1 = sf.get_row(1) * matOneRingVertex;
        Matrix a_temp = sf.get_row(2) * matOneRingVertex; ///< a_2
        a_temp = cross_row(a_1, a_temp); ///< a_3 = a_1 x a_2
        double sqa = a_temp.calculate_norm();
        //VERSION 2
        Matrix x = sf.get_row(0) * matOneRingVertex;
        Matrix a_1 = sf.get_row(1) * matOneRingVertex;
        Matrix a_2 = sf.get_row(2) * matOneRingVertex;
        Matrix a_3 = cross_row(a_1, a_2);
        double sqa = a_3.calculate_norm();
        */

        // VERSION 3
        const Matrix &sf = sampleRows[j];
        Matrix x = sf.get_row(0) * matOneRingVertex;
        Matrix a_3 = cross_row(sf.get_row(1) * matOneRingVertex, sf.get_row(2) * matOneRingVertex);
        double sqa = a_3.calculate_norm();
        // VERSION 3

        // double s = sqa;
        // vector<double> d = a_3 / sqa;

        double coeff = param.gaussQuadratureCoeff(j, 0);
        area += 0.5 * coeff * sqa; // Update the accumulated area
        // v = 1/3 * s * dot(x, d) <<< tetrahedron volume
        // namely -> double v = dot_row(x, a_3) / 3.0;
        // volume += 0.5 * coeff * v; // Update the accumulated volume
        volume += kSignedVolumeQuadratureFactor * coeff * dot_row(x, a_3); // Update the accumulated volume
    }
}

// Computes a matrix containing the coordinates of the one-ring vertices for the input face.
Matrix Mesh::get_one_ring_vertex_matrix(const Face &face)
{
    const int nOneRingVertices = face.oneRingVertices.size();
    Matrix matOneRingVertex = mat_calloc(nOneRingVertices, 3);

    for (int j = 0; j < nOneRingVertices; j++)
    {
        int oneRingVerticesIndex = face.oneRingVertices[j];
        matOneRingVertex.set_row_from_col(j, vertices[oneRingVerticesIndex].coord, 0);
    }

    return matOneRingVertex;
}

void Mesh::calculate_element_area_volume()
{
#pragma omp parallel for
    for (Face& face : faces)
    {
        // Ghost faces won't contribute to the membrane area
        // I changed this from boundary faces to ghost faces since
        // by def, ghost faces are the ones that do not contribute to physical
        // calculations - Y Ying
        if (face.isGhost)
            continue;

        double area = 0.0;   ///< The accumulated area for a membrane.
        double volume = 0.0; ///< The accumulated volume for a membrane.
        const int nOneRingVertices = static_cast<int>(face.oneRingVertices.size());

        // Geometry reads exactly the rows energy and force read, through the
        // same dispatch. Before WP5 this arm ran its own explicit subdivision
        // recursion over param.subMatrix, bounded by an uninitialized
        // param.subDivideTimes -- so area and volume were integrated by
        // different code, to a different depth, than the energy they are
        // constrained against.
        if (nOneRingVertices == 12)
        {
            const Matrix matOneRingVertex = get_one_ring_vertex_matrix(face);
            enumerate_gauss_quadrature_point_area_volume(param.shapeFunctions, matOneRingVertex,
                                                        area, volume);
        }
        else if (nOneRingVertices >= kMinIrregularValence + 6 &&
                 nOneRingVertices <= kMaxIrregularValence + 6)
        {
            const Matrix matOneRingVertex = get_one_ring_vertex_matrix(face);
            const int valence = nOneRingVertices - 6;
            for (int d = 0; d < irregularRows.depth_for(valence); d++)
            {
                for (int c = 0; c < kRegularChildrenPerStep; c++)
                {
                    enumerate_gauss_quadrature_point_area_volume(
                        irregularRows.rows_for_child(valence, d, c), matOneRingVertex, area,
                        volume);
                }
            }
        }
        // Any other width means no complete one-ring -- a boundary face, with
        // no limit surface and so no area or volume of its own.

        face.elementArea = area;
        face.elementVolume = volume;
    }
}

void Mesh::sum_membrane_area_and_volume(double &area, double &volume)
{
    area = 0.0;
    volume = 0.0;
#pragma omp parallel for reduction(+                   \
                                   : area) reduction(+ \
                                                     : volume)
    for (Face& face : faces)
    {
        if (!face.isGhost)
        {
            area += face.elementArea;
            volume += face.elementVolume;
        }
    }
}
