#include "Dynamics.hpp"

using namespace std;

// calculate mesh2surface matrix
// surface = mesh2surface * mesh

DynamicMesh::DynamicMesh(Param &srcParam) :
Mesh(srcParam)
{
}

void DynamicMesh::setup_flat() {
    // Call the superclass's setup_flat()
    Mesh::setup_flat();

    // Assign mesh2surface and surface2mesh
    if (param.VERBOSE_MODE)
    {
        std::cout << "[DynamicMesh] Assigning conversion matrices." << std::endl;
    }
    assign_mesh2surface();
    if (param.VERBOSE_MODE)
    {
        std::cout << "[DynamicMesh] mesh2surface matrix : " << std::endl << mesh2surface << std::endl;
    }
    mesh2surface.get_inverted(surface2mesh);
    if (param.VERBOSE_MODE)
    {
        std::cout << "[DynamicMesh] surface2mesh matrix : " << std::endl << surface2mesh << std::endl;
    }
    matMesh = mat_calloc(vertices.size(), 3);
    matSurface = mat_calloc(vertices.size(), 3);
}

void DynamicMesh::assign_mesh2surface()
{
    // Initialize mesh2surface matrix
    mesh2surface = mat_calloc(vertices.size(), vertices.size());
    
    // Iterate through all vertices
    for (int i = 0; i < vertices.size(); i++)
    {
        // find adjacent face that is not ghost
        int indexAdjFace = -1; // -1 means only ghost
        for (int jAdjFace = 0; jAdjFace < vertices[i].adjacentFaces.size(); jAdjFace++)
        {
            const Face& adjF = faces[vertices[i].adjacentFaces[jAdjFace]];
            if (!adjF.isGhost && !adjF.isBoundary)
            {
                // assign non ghost / boundary face index to indexAdjFace
                indexAdjFace = vertices[i].adjacentFaces[jAdjFace];
            }
        }
        if (indexAdjFace < 0)
        {
            // If only ghost or boundary faces are adjacent to vertex, shape function = (1, 0, 0, ...)
            mesh2surface.set(i, i, 1.0);
        }
        else
        {
            // If there exists non-ghost / boundary face
            // find vwu of vertex on face
            //! vwu assumed to be same sequence as defined in adjacent face
            Face &faceAdj = faces[indexAdjFace];
            std::vector<int> &adjVertices = faceAdj.adjacentVertices;
            int vwuInd = -1;
            for (int j = 0; j < adjVertices.size(); j++)
            {
                if ((adjVertices)[j] == i)
                {
                    vwuInd = j;
                }
            }
            Matrix vwuAdj = mat_calloc(1, 3); // get_shapefunction takes row vector
            vwuAdj.set(0, vwuInd, 1.0);
            // std::cout<< i << " ," << vwuAdj << endl;

            // use vwu to get sf - irregular patch -currently not used!
            // @todo implement support for irregular patch
            // Matrix sf(7, 12);
            // get_shapefunction(vwuAdj, sf);
            // vector<double> sfAdj = determine_ShapeFunctions(vwuAdj)[0]; // transposed sf

            // set mat(AB) (mesh2surface) to sf
            //@TODO: currently due to vwu not in order in adjacent vertex of face
            // using 0.5 / 0.0833333 directly for regular patches
            std::vector<int> &iAdjVertices = vertices[i].adjacentVertices;

            mesh2surface.set(i, i, 0.5);
            for (const int &iAdj : iAdjVertices)
            {
                mesh2surface.set(i, iAdj, 0.5 / 6.0);
            }
        }
        // check if ghost / boundary faces are correctly recongnized
        // std::cout << "vertex: " << i << ", indexAdjFace: " << indexAdjFace << endl;
    }
    std::cout << "[DynamicsMesh::asign_mesh2surface] Generate mesh2surface conversion matrix: " << std::endl;
    std::cout << mesh2surface << std::endl;
}

void DynamicMesh::update_vertices_mat_with_vector()
{
    for (int i = 0; i < vertices.size(); i++){
        for (int j = 0; j < 3; j++){
            matMesh.set(i, j, vertices[i].coord(j,0));
        }
    }
}

void DynamicMesh::update_vertices_vector_with_mat()
{
    for (int i = 0; i < vertices.size(); i++){
        for (int j = 0; j < 3; j++){
            vertices[i].coord.set(j, 0, matMesh(i, j));
        }
    }
}

