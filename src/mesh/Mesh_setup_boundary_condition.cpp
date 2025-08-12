#include "mesh/Mesh.hpp"

void Mesh::determine_ghost_vertices_faces()
{
    // no ghost vertices for fixed BC
    if (param.boundaryCondition == BoundaryType::Fixed){
        return;
    }

    // 1. set ghost vertex
    int nVertices = (param.nFaceX + 1) * (param.nFaceY + 1);
    vector<int> topBottom;
    vector<int> leftRight;
    switch (param.boundaryCondition)
    {
    case BoundaryType::Periodic:
        topBottom.insert(topBottom.end(), {0, 1, 2, param.nFaceY - 2, param.nFaceY - 1, param.nFaceY});
        leftRight.insert(leftRight.end(), {0, 1, 2, param.nFaceX - 2, param.nFaceX - 1, param.nFaceX});
        break;
    case BoundaryType::Free:
        topBottom.insert(topBottom.end(), {0, param.nFaceY});
        leftRight.insert(leftRight.end(), {0, param.nFaceX});
        break;
    }
    // top and bottom ghost vertex
    for (int k = 0; k < topBottom.size(); k++)
    {
        int j = topBottom[k]; // rows
#pragma omp parallel for
        for (int i = 0; i < param.nFaceX + 1; i++) // iterate columns
        {
            int index = (param.nFaceX + 1) * j + i;
            vertices[index].isGhost = true;
            vertices[index].type = VertexType::Ghost;
        }
    }
    // left and right ghost vertex
    for (int k = 0; k < leftRight.size(); k++)
    {
        int i = leftRight[k]; // columns
#pragma omp parallel for
        for (int j = 0; j < param.nFaceY + 1; j++) // iterate rows
        {
            int index = (param.nFaceX + 1) * j + i;
            vertices[index].isGhost = true;
            vertices[index].type = VertexType::Ghost;
        }
    }
    // 2. set ghost face
    int nFaces = param.nFaceY * param.nFaceX * 2;
    topBottom.clear();
    leftRight.clear();

    switch (param.boundaryCondition)
    {
    case BoundaryType::Periodic:
        topBottom.insert(topBottom.end(), {0, 1, 2, param.nFaceY - 3, param.nFaceY - 2, param.nFaceY - 1});
        leftRight.insert(leftRight.end(), {0, 1, 2, param.nFaceX - 3, param.nFaceX - 2, param.nFaceX - 1});
        break;
    case BoundaryType::Free:
        topBottom.insert(topBottom.end(), {0, param.nFaceY - 1});
        leftRight.insert(leftRight.end(), {0, param.nFaceX - 1});
        break;
    }

    for (int k = 0; k < topBottom.size(); k++)
    {
        int j = topBottom[k]; // rows
#pragma omp parallel for
        for (int i = 0; i < param.nFaceX; i++) // iterate columns
        {
            int index = 2 * param.nFaceX * j + i * 2;
            faces[index].isGhost = true;
            faces[index + 1].isGhost = true;
        }
    }
    for (int k = 0; k < leftRight.size(); k++)
    {
        int i = leftRight[k]; // columns
#pragma omp parallel for
        for (int j = 0; j < param.nFaceY; j++) // iterate rows
        {
            int index = 2 * param.nFaceX * j + i * 2;
            faces[index].isGhost = true;
            faces[index + 1].isGhost = true;
        }
    }
    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::determine_ghost_vertices_faces] Ghost vertices and faces set -- indices:" << std::endl;
        for (const Face& face : faces)
        {
            if (face.isGhost)
            std::cout << face.index << ", ";
        }
        std::cout<< "end of indices" << std::endl;
    }

}