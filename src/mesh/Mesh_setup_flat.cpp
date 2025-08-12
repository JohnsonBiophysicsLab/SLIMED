#include "mesh/Mesh.hpp"

void Mesh::setup_flat()
{
    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::setup_flat] Setting up flat membrane." << std::endl;
    }
    set_axes_division_flat();
    set_vertices_faces_flat();
    set_adjacent_faces_of_vertices_sorted();
    set_adjacent_vertices_of_vertices_sorted();
    set_adjacent_faces_of_faces();
    sort_vertices_on_faces();
    determine_ghost_vertices_faces();
    set_one_ring_vertices_sorted();
    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::setup_flat] Finished setting up flat membrane." << std::endl;
    }
}

void Mesh::set_axes_division_flat()
{
    param.nFaceX = round(param.sideX / param.lFace); // number of edges along x-axis
    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::set_axes_division] nFaceX = " << param.nFaceX << std::endl;
    }
    param.dFaceX = param.sideX / param.nFaceX; // calculate actual face side length x-axis
    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::set_axes_division] dFaceX = " << param.dFaceX << std::endl;
    }
    param.dFaceY = sqrt(3.0) / 2.0 * param.dFaceX; // sqrt(3.0)/2.0; calculate perpendicular face side length along y-axis
    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::set_axes_division] dFaceY = " << param.dFaceY << std::endl;
    }
    param.nFaceY = round(param.sideY / param.dFaceY); // number of edges along y-axis

    if (param.nFaceY & 1)
    {
        // returns 1 if odd, else 0,
        // @ref https://stackoverflow.com/questions/62030964/
        param.nFaceY += 1; // make nFaceY even for fixed shape in Y axis
    }

    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::set_axes_division] nFaceY = " << param.nFaceY << std::endl;
    }
}

void Mesh::set_vertices_faces_flat()
{
    // step 1. Check axes division
    if (param.nFaceX < 0)
    {
        set_axes_division_flat(); // divide axes if not done yet
    }

    std::cout << "[Mesh::param]" << std::endl << param << std::endl;

    // half lengths used to shift the membrane around (0, 0, 0)
    double lxHalf = param.nFaceX * param.dFaceX * 0.5; // half of actual initial length in x-dir
    double lyHalf = param.nFaceY * param.dFaceY * 0.5; // half of actual initial length in y-dir

    // step 2. Initialize vertices and faces
    int nVertices = (param.nFaceX + 1) * (param.nFaceY + 1); // number of vertices'
    vertices = vector<Vertex>(nVertices);                    // declare local vertices list

    int nFaces = param.nFaceX * param.nFaceY * 2; // number of faces
    faces = vector<Face>(nFaces);

    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::set_vertices_faces_flat] nVertices = " << nVertices
                  << ", nFaces = " << nFaces << std::endl;
    }

#pragma omp parallel for
    for (int j = 0; j < param.nFaceY + 1; j++) // iterate along y-axis
    {
        for (int i = 0; i < param.nFaceX + 1; i++) // iterate along x-axis
        {
            int vertexIndex = (param.nFaceX + 1) * j + i; // index of current vertex
            int faceIndex = 2 * param.nFaceX * j + i * 2; // index of current face
            double xCoord = i * param.dFaceX;             // x coordinate of current vertex
            double yCoord = j * param.dFaceY;             // y coordinate of current vertex
            int node1, node2, node3, node4;               // nodes for faces - sequence of vertices is always defined counter clockwise
                                                          // face1 = triangle1->2->3
                                                          // face2 = triangle2->4->3
            if (j & 1) // true if j is odd
            {
                // odd j: face starts at the lower left (sequence of vertices is always defined counter clockwise)
                node1 = (param.nFaceX + 1) * j + i;
                node2 = (param.nFaceX + 1) * (j + 1) + i;
                node3 = (param.nFaceX + 1) * j + (i + 1);
                node4 = (param.nFaceX + 1) * (j + 1) + (i + 1);
            }
            else // j is even
            {
                xCoord += param.dFaceX / 2.0; // compensate for zig-zag shift along y-axis
                // even j: face starts at the top left
                node3 = (param.nFaceX + 1) * j + i;
                node1 = (param.nFaceX + 1) * (j + 1) + i;
                node4 = (param.nFaceX + 1) * j + (i + 1);
                node2 = (param.nFaceX + 1) * (j + 1) + (i + 1);
            }

            // vertex
            vertices[vertexIndex].index = vertexIndex;
            vertices[vertexIndex].coord.set(0, 0, xCoord - lxHalf); // shifted x-coord
            vertices[vertexIndex].coord.set(1, 0, yCoord - lyHalf); // shifted y-coord
            vertices[vertexIndex].coord.set(2, 0, 0.0);

            if (i != param.nFaceX && j != param.nFaceY)
            {
                // face: triangle 1->2->3
                faces[faceIndex].index = faceIndex; // assign face index
                // setup adjacent vertex for face
                faces[faceIndex].adjacentVertices = std::vector<int>(3);
                faces[faceIndex].adjacentVertices[0] = node1; // assign adjacent vertices
                faces[faceIndex].adjacentVertices[1] = node2;
                faces[faceIndex].adjacentVertices[2] = node3;
                // face: triangle 2->4->3
                faceIndex++;
                faces[faceIndex].index = faceIndex; // assign face index
                // setup adjacent vertex for face
                faces[faceIndex].adjacentVertices = std::vector<int>(3);
                faces[faceIndex].adjacentVertices[0] = node2; // assign adjacent vertices
                faces[faceIndex].adjacentVertices[1] = node4;
                faces[faceIndex].adjacentVertices[2] = node3;
            }
        }
    }
    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::set_vertices_faces_flat] Assigned vertices and faces in the member of current mesh object." << std::endl;
    }
}

