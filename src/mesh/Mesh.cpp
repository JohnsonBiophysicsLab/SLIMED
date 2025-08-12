#include "mesh/Mesh.hpp"

Mesh::Mesh(Param &srcParam) : param(srcParam)
{
    // Calculate element triangle area
    param.elementTriangleArea0 = sqrt(3.0) / 4.0 * param.lFace * param.lFace;

    // Compute Gaussian quadrature weights and coordinates
    get_gauss_quadrature_weight_VWU(param.gaussQuadratureN, param.VWU, param.gaussQuadratureCoeff);

    // Compute shape functions for each element and store them in `param`
    //param.shapeFunctions = std::vector<Matrix>(param.VWU.nrow());
    get_shapefunction_vector(param.VWU, param.shapeFunctions);

    // Generate matrices used for subdivision of irregular patches in the mesh
    get_subdivision_matrices(param.subMatrix.irregM,
                             param.subMatrix.irregM1,
                             param.subMatrix.irregM2,
                             param.subMatrix.irregM3,
                             param.subMatrix.irregM4);

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

    // Generate matrices used for subdivision of irregular patches in the mesh
    get_subdivision_matrices(param.subMatrix.irregM,
                             param.subMatrix.irregM1,
                             param.subMatrix.irregM2,
                             param.subMatrix.irregM3,
                             param.subMatrix.irregM4);
    
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

    // step 3. Link neighboring geometric components
    set_adjacent_faces_of_vertices_sorted();
    set_adjacent_vertices_of_vertices_sorted();
    determine_ghost_vertices_faces();
    set_one_ring_vertices_sorted();

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
    for (Face face : faces)
    {
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
    const Matrix &matOneRingVertex,
    double &area,
    double &volume)
{
#pragma omp parallel for reduction(+ \
                                   : area, volume)
    for (int j = 0; j < param.gaussQuadratureCoeff.nrow(); j++)
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
        Matrix &sf = param.shapeFunctions[j];
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
        volume += 0.16666666666 * coeff * dot_row(x, a_3); // Update the accumulated volume
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
    // five matrix used for subdivision of the irregular patch
    // M(17,11), M1(12,17), M2(12,17), M3(12,17), M4(11,17);
    // alias for subMatrix - does not cause
    // extra deepcopy for the subMatrix - Y Ying
    const auto &subMat = param.subMatrix; // Get a reference to the subMatrix object
    const Matrix &M = subMat.irregM;
    const Matrix &M1 = subMat.irregM1;
    const Matrix &M2 = subMat.irregM2;
    const Matrix &M3 = subMat.irregM3;
    const Matrix &M4 = subMat.irregM4;
#pragma omp parallel for
    for (Face& face : faces)
    {
        // Ghost faces won't contribute to the membrane area
        // I changed this from boundary faces to ghost faces since
        // by def, ghost faces are the ones that do not contribute to physical
        // calculations - Y Ying
        if (face.isGhost)
            continue;
        

        // Variables initialization
        double area = 0.0;                                        ///< The accumulated area for a membrane.s
        double volume = 0.0;                                      ///< The accumulated volume for a membrane.
        const int nOneRingVertices = face.oneRingVertices.size(); ///< The number of one-ring vertices for a face in a mesh.

        // Compute area and volume for a regular patch
        switch (nOneRingVertices)
        {
        case 12:
            {
                // The matrix representing the coordinates of the one ring vertices.
                Matrix matOneRingVertex = get_one_ring_vertex_matrix(face);

                // Use Gaussian quadrature with 3 points to compute area and volume
                enumerate_gauss_quadrature_point_area_volume(matOneRingVertex,
                                                            area, volume);
            }
            break;

        case 11:
            {
                // irregular patch
                // The matrix representing the coordinates of the one ring vertices.
                Matrix matOrigOneRingVertex = get_one_ring_vertex_matrix(face);
                Matrix matNewOneRingVertex; // store the one ring vertex after subdivision
                Matrix matNewNodes17;       // store the coordinates of 17 new nodes after subdivision

                // using subdivision to estimate irregular patch
                for (int j = 0; j < param.subDivideTimes; j++)
                {
                    //(17 x 3) = (17 x 11) * (11 x 3)
                    matNewNodes17 = M * matOrigOneRingVertex; // 17 new nodes after subdivision

                    matNewOneRingVertex = M1 * matNewNodes17; // element 1
                    enumerate_gauss_quadrature_point_area_volume(matNewOneRingVertex,
                                                                area, volume);

                    matNewOneRingVertex = M2 * matNewNodes17; // element 2
                    enumerate_gauss_quadrature_point_area_volume(matNewOneRingVertex,
                                                                area, volume);

                    matNewOneRingVertex = M3 * matNewNodes17; // element 3
                    enumerate_gauss_quadrature_point_area_volume(matNewOneRingVertex,
                                                                area, volume);

                    matOrigOneRingVertex = M4 * matNewNodes17; // element 4, still irregular patch
                }
            }
            break;
        }
        
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

