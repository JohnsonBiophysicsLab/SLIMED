#include "Dynamics.hpp"

using namespace std;

// calculate real point relative to the given ghost point (index(real) - index(given)) in periodic boundary condition
// returns 0 if real point (not on 4th ring) given
// an arbitrary real point is chosen if given 4th
// input: i = index of given pt; n = number of segments of columns; m = number of segments of rows
int DynamicMesh::get_relative_pt_periodic(int i, int n, int m)
{
    int indexRelative = 0;
    int colnum = i % (n + 1);
    int rownum = i / (n + 1);
    // special edge case
    if (colnum == n - 3 && rownum == 3){
        return ((m - 6) * (n + 1));
    }
    if (colnum == 3 && rownum == m - 3){
        return (n - 6);
    }
    // get real point based on the position
    if (colnum < 4)
    {
        indexRelative += n - 6;
    }
    else if (colnum > n - 4)
    {
        indexRelative -= n - 6;
    }
    if (rownum < 4)
    {
        indexRelative += (m - 6) * (n + 1);
    }
    else if (rownum > m - 4)
    {
        indexRelative -= (m - 6) * (n + 1);
    }
    return indexRelative;
}

// post move step processing on triangular mesh - syncing ghost vertices and 4th ring real
// vertices with the corresponding mesh points
// G = Ghost
// R = Real
// B = Boundary (which is also real)
// 0 | 1 | 2 | 3 | 4 | 5 .... N - 4 | N - 3 | N - 2 | N - 1 | N
// G | G | G | B | R | R ....  R    |   B   |   G   |   G   | G
void DynamicMesh::postprocess_ghost_periodic()
{
    double original = 0.0;
    // 1. process real point
    for (int i = 0; i < vertices.size(); i++)
    {
        int indexRelative = 0;
        // calculate 4th point pair
        if (!(vertices[i].isGhost))
        {
            indexRelative = get_relative_pt_periodic(i, param.nFaceX, param.nFaceY);
            //std::cout << "index: " << i << "; relative = " << indexRelative << std::endl;
            if (indexRelative != 0)
            {
                if (vertices[i].type != VertexType::PeriodicReflectiveBoundary){
                    vertices[i].type = VertexType::PeriodicReflectiveBoundary;
                    vertices[i].reflectiveVertexIndex = i + indexRelative;
                }
                double displacement = 0.0;
                for (int j = 0; j < 3; j++)
                {
                    // coord before displacement is stored in vertices
                    // coord after displacement is stored in matMesh
                    // get displacement based on relative pt of this pt
                    displacement = matMesh(i + indexRelative, j) - vertices[i + indexRelative].coord(j, 0);
                    original = vertices[i].coord(j, 0); // get original coord of this pt
                    matMesh.set(i, j, original + displacement); // apply same displacment to this pt
                }
            }
        }
    }

    // 2. process ghost point
    for (int i = 0; i < vertices.size(); i++)
    {
        // calcualate real point for periodic boundary condition
        int indexRelative = 0;
        if (vertices[i].isGhost)
        {
            // in order to sync up 4th ring pts
            indexRelative = get_relative_pt_periodic(i, param.nFaceX, param.nFaceY);
            if (indexRelative != 0)
            {

                vertices[i].type = VertexType::Ghost;
                vertices[i].reflectiveVertexIndex = i + indexRelative;

                double displacement = 0.0;
                for (int j = 0; j < 3; j++)
                {
                    // coord before displacement is stored in vertices
                    // coord after displacement is stored in matMesh
                    // get displacement based on relative pt of this pt
                    displacement = matMesh(i + indexRelative, j) - vertices[i + indexRelative].coord(j, 0);
                    original = vertices[i].coord(j, 0); // get original coord of this pt
                    matMesh.set(i, j, original + displacement); // apply same displacment to this pt
                }
            }
        }
    }
}
