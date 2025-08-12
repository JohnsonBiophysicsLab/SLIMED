#include "mesh/Mesh.hpp"


void Mesh::set_adjacent_faces_of_vertices_sorted()
{
    // 1. Get adjacent faces unsorted
    // Initialize vector of empty vectors for adjacent faces
    vector<vector<int>> adjFaces(vertices.size());
    
    // Populate adjacent faces for each vertex
    for (int j = 0; j < faces.size(); j++)
    {
        for (int k = 0; k < 3; k++)
        {
            int i = faces[j].adjacentVertices[k];
            adjFaces[i].push_back(j);
        }
    }
    
    // Transfer adjacent faces from temp vector to each vertex
    for (int i = 0; i < vertices.size(); i++)
    {
        vertices[i].adjacentFaces = std::move(adjFaces[i]);
    }

    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::set_adjacent_faces_of_vertices_sorted] Got adjacent faces unsorted." << std::endl;
    }

    // 2. Sort faces so that adjacentFaces that are adjacent to each
    // other are of +/- 1 index.
    for (int j = 0; j < param.nFaceY + 1; j++) // iterate along y-axis
    {
        for (int i = 0; i < param.nFaceX + 1; i++) // iterate along x-axis
        {
            // initialize adjacent vertex indices
            std::vector<int> adjacentFacesSorted(6, 0);
            if (j & 1) // j is odd
            {
                adjacentFacesSorted[0] = 2 * (param.nFaceX * (j - 1) + i) - 2;
                adjacentFacesSorted[1] = 2 * (param.nFaceX * j + i) - 2;
                adjacentFacesSorted[2] = 2 * (param.nFaceX * j + i) - 1;
                adjacentFacesSorted[3] = 2 * (param.nFaceX * j + i);
                adjacentFacesSorted[4] = 2 * (param.nFaceX * (j - 1) + i);
                adjacentFacesSorted[5] = 2 * (param.nFaceX * (j - 1) + i) - 1;
            }
            else // j is even
            {
                adjacentFacesSorted[0] = 2 * (param.nFaceX * (j - 1) + i) - 1;
                adjacentFacesSorted[1] = 2 * (param.nFaceX * j + i) - 1;
                adjacentFacesSorted[2] = 2 * (param.nFaceX * j + i);
                adjacentFacesSorted[3] = 2 * (param.nFaceX * j + i) + 1;
                adjacentFacesSorted[4] = 2 * (param.nFaceX * (j - 1) + i) + 1;
                adjacentFacesSorted[5] = 2 * (param.nFaceX * (j - 1) + i);
            }

            // pop vertices that do not exist
            for (int k = 5; k >= 0; k--)
            {
                std::vector<int> *adjacentFacesUnsorted = &(vertices[(1 + param.nFaceX) * j + i].adjacentFaces);
                if (std::find(adjacentFacesUnsorted->begin(),
                              adjacentFacesUnsorted->end(),
                              adjacentFacesSorted[k]) == adjacentFacesUnsorted->end()) // if value not in unsorted
                {
                    adjacentFacesSorted.erase(adjacentFacesSorted.begin() + k);
                }
            }
            vertices[(1 + param.nFaceX) * j + i].adjacentFaces = adjacentFacesSorted;
        }
    }

    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::set_adjacent_faces_of_vertices_sorted] Faces sorted." << std::endl;
    }
}

/**
 * @brief Return true if two faces share edge.
 * 
 * @brief Will also return true if face1 and face2 are the same face.
 * 
 * @param face1 
 * @param face2 
 * @return true 
 * @return false 
 */
bool Mesh::faces_share_edge(const Face& face1, const Face& face2){
    // If face1 and face2 have two identical vertices then they 
    // share edge
    // Use std::set_intersection to find common elements
    std::vector<int> commonElements;
    return faces_share_edge(face1, face2, commonElements);
}

bool Mesh::faces_share_edge(const Face& face1, const Face& face2, std::vector<int>& commonElements){
    // If face1 and face2 have two identical vertices then they 
    // share edge
    // Use std::set_intersection to find common elements
    std::set_intersection(face1.adjacentVertices.begin(), face1.adjacentVertices.end(),
                          face2.adjacentVertices.begin(), face2.adjacentVertices.end(),
                          std::back_inserter(commonElements));

    // Check if there are at least two common vertices
    return commonElements.size() >= 2;
}

/**
 * @brief Set adjacentFaces properties of faces based on the current
 * geometry of mesh.
 * 
 * This function iterates over the faces of the mesh and populates the
 * adjacentFaces property of each face by finding neighboring faces that
 * share an edge.
 */
void Mesh::set_adjacent_faces_of_faces(){
    // iterate over faces and add adjacent faces
    for (int i = 0; i < faces.size(); ++i){
        // Initialize the adjacentFaces vector for the current face
        faces[i].adjacentFaces = std::vector<int>(3);
        int adjacentFaceIndex = 0;

        // Iterate over all faces to find adjacent faces
        for (int j = 0; j < faces.size(); ++j){
            // Check if faces i and j share an edge
            if (faces_share_edge(faces[i], faces[j])){
                // Add the index of the adjacent face to the current face's adjacentFaces vector
                faces[i].adjacentFaces[adjacentFaceIndex] = j;
                // Move to the next slot in the adjacentFaces vector
                ++adjacentFaceIndex;
            }
        }
    }
}


void Mesh::set_adjacent_vertices_of_vertices_sorted()
{

#pragma omp parallel for
    // iterate over vertices and add adjacent vertices of adjacent faces
    for (int i = 0; i < vertices.size(); i++)
    {
        vector<int> adjacentVerticesTmp;
        for (int j = 0; j < vertices[i].adjacentFaces.size(); j++)
        {
            int faceIndex = vertices[i].adjacentFaces[j];
            for (int k = 0; k < faces[faceIndex].adjacentVertices.size(); k++)
            {
                int vertexIndex = faces[faceIndex].adjacentVertices[k];
                if (vertexIndex != i)
                {
                    bool isListed = false;
                    // check if both vertices are in adjacent vertices of a face
                    for (int m = 0; m < adjacentVerticesTmp.size(); m++)
                    {
                        if (vertexIndex == adjacentVerticesTmp[m])
                        {
                            isListed = true;
                        }
                    }
                    if (isListed == false)
                    {
                        adjacentVerticesTmp.push_back(vertexIndex);
                    }
                }
            }
        }
        vertices[i].adjacentVertices = adjacentVerticesTmp;
    }
    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::set_adjacent_vertices_of_vertices_sorted] Adjacent vertices set." << std::endl;
    }
}

int Mesh::find_opposite_node_index(const int &node1, const int &node2, const int &node3)
{
    int node = -1;
    for (int i = 0; i < vertices[node1].adjacentVertices.size(); i++)
    {
        int nodetmp1 = vertices[node1].adjacentVertices[i];
        for (int j = 0; j < vertices[node2].adjacentVertices.size(); j++)
        {
            int nodetmp2 = vertices[node2].adjacentVertices[j];
            if (nodetmp1 == nodetmp2 && nodetmp1 != node3)
            {
                node = nodetmp1;
            }
        }
    }
    if (node == -1)
    {
        if (param.VERBOSE_MODE) {

        }
        cout << "No efficent oneRingVerticesIndex is found! Node1 = " << node1 << ", Node2 = " << node2 << ", Node3 = " << node3 << endl;
    }
    return node;
}

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
void Mesh::sort_vertices_on_faces()
{
    bool isAllFacesSorted = false;
    // Initialize a vector of booleans with given length and set all elements to false
    std::vector<bool> isFaceSorted(faces.size(), false);

    // Assume all face sort sequence will be based on faces[0]
    isFaceSorted[0] = true;

    // Loop through all faces in the mesh
    while (!isAllFacesSorted)
    {
        // Loop through all faces in the mesh
        for (int iFace = 0; iFace < faces.size(); iFace++) {
            if (!(isFaceSorted[iFace]))
            {
                for (int j = 0; j < 3; j++){
                    int jAdjFace = faces[iFace].adjacentFaces[j];
                    if (isFaceSorted[jAdjFace] && !isFaceSorted[iFace]){
                        std::vector<int> commonElements;
                        faces_share_edge(faces[iFace], faces[jAdjFace], commonElements);
                        
                        // Concatenate the first vector with itself to check for wrapping-around sequences
                        std::vector<int> extendedAdjFacesj = faces[jAdjFace].adjacentFaces;
                        extendedAdjFacesj.insert(extendedAdjFacesj.end(),
                                faces[jAdjFace].adjacentFaces.begin(),
                                faces[jAdjFace].adjacentFaces.end());
                        std::vector<int> extendedAdjFacesi = faces[iFace].adjacentFaces;
                        extendedAdjFacesi.insert(extendedAdjFacesi.end(),
                                faces[iFace].adjacentFaces.begin(),
                                faces[iFace].adjacentFaces.end());
                        
                        // If extendedAdjFacesj DOES NOT contain commonElements, reverse common Elements
                        if (std::search(extendedAdjFacesj.begin(), extendedAdjFacesj.end(),
                                commonElements.begin(), commonElements.end()) == extendedAdjFacesj.end())
                        {
                            std::reverse(commonElements.begin(), commonElements.end());
                        }

                        // If extendedAdjFacesi CONTAINS commonElemnts reverse it
                        if (std::search(extendedAdjFacesi.begin(), extendedAdjFacesi.end(),
                                commonElements.begin(), commonElements.end()) != extendedAdjFacesi.end())
                        {
                            std::reverse(faces[iFace].adjacentFaces.begin(), faces[iFace].adjacentFaces.end());
                        }

                        // Set processed flag to true
                        isFaceSorted[iFace] = true;
                    }
                }
            }
        }
        // Recalculate isAllFacesSorted
        isAllFacesSorted = true;
        for (bool value : isFaceSorted) {
            if (!value) {
                isAllFacesSorted = false;
                break;
            }
        }
    }
}

// To find out the one-ring vertices aound face_i. It should be 12 for the flat surface because we set it up only with regular patch.
// The boundary faces do not have complete one-ring, neither it will be called in the code, so no need to store their one-ring-vertex
void Mesh::set_one_ring_vertices_sorted()
{
// two types of patch: 1. regular patch with 12 one-ring vertices, each vertex has 6 closest nodes
//                     2. irregular patch with 11 one-ring vertices, one vertex has 5 closest-nodes
#pragma omp parallel for
    for (Face& face : faces)
    {
        // skip for ghost faces -> not enough neighboring mesh triangles!
        if (face.isGhost)
        {
            continue;
        }

        // int d1, d2, d3, d5, d6, d9, d10, d11, d12
        int node0 = face.adjacentVertices[0];
        int node1 = face.adjacentVertices[1];
        int node2 = face.adjacentVertices[2];
        // regular patch, all three nodes have 6 neighbor faces or vertices.
        if (vertices[node0].adjacentVertices.size() == 6 &&
            vertices[node1].adjacentVertices.size() == 6 &&
            vertices[node2].adjacentVertices.size() == 6)
        {
            // make sure d4, d7, d8 are in anti-clock-wise order
            int d4 = node0;
            int d7 = node1;
            int d8 = node2;
            //
            Matrix coord4 = vertices[node0].coord;
            Matrix coord7 = vertices[node1].coord;
            Matrix coord8 = vertices[node2].coord;
            // calculate CoM of the triangle
            Matrix center = 1.0 / 3.0 * (coord4 + coord7 + coord8);
            if (dot_col(center, cross_col(coord7 - coord4, coord8 - coord4)) < 0)
            { // switch d7 and d8 position
                d7 = face.adjacentVertices[2];
                d8 = face.adjacentVertices[1];
                face.adjacentVertices[1] = d7;
                face.adjacentVertices[2] = d8;
            }
            int d3 = find_opposite_node_index(d4, d7, d8);
            int d11 = find_opposite_node_index(d7, d8, d4);
            int d5 = find_opposite_node_index(d4, d8, d7);
            int d1 = find_opposite_node_index(d3, d4, d7);
            int d2 = find_opposite_node_index(d4, d5, d8);
            int d6 = find_opposite_node_index(d3, d7, d4);
            int d9 = find_opposite_node_index(d8, d5, d4);
            int d10 = find_opposite_node_index(d7, d11, d8);
            int d12 = find_opposite_node_index(d8, d11, d7);
            face.oneRingVertices = vector<int>{d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12};
            
            // irregular patch, one node has 5 neighbors and the other two nodes have 6 neighbors.
        }
        else if (vertices[node0].adjacentVertices.size() == 5 &&
                 vertices[node1].adjacentVertices.size() == 5 &&
                 vertices[node2].adjacentVertices.size() == 5)
        {
            int d4, d7, d8;
            // make sure d4 is the one has 5 neighbors, and d4-d7-d8 are in anti-clock-wise order
            if (vertices[node0].adjacentFaces.size() == 5)
            {
                d4 = node0;
                d7 = node1;
                d8 = node2;
            }
            else if (vertices[node1].adjacentFaces.size() == 5)
            {
                d4 = node1;
                d7 = node2;
                d8 = node0;
            }
            else if (vertices[node2].adjacentFaces.size() == 5)
            {
                d4 = node2;
                d7 = node0;
                d8 = node1;
            }
            //
            Matrix coord4 = vertices[node0].coord;
            Matrix coord7 = vertices[node1].coord;
            Matrix coord8 = vertices[node2].coord;
            // calculate CoM of the triangle
            Matrix center = 1.0 / 3.0 * (coord4 + coord7 + coord8);
            if (dot_col(center, cross_col(coord7 - coord4, coord8 - coord4)) < 0)
            { // switch d7 and d8 position
                int nodetmp = d7;
                d7 = d8;
                d8 = nodetmp;
                face.adjacentVertices[0] = d4;
                face.adjacentVertices[1] = d7;
                face.adjacentVertices[2] = d8;
            }
            int d3 = find_opposite_node_index(d4, d7, d8);
            int d11 = find_opposite_node_index(d7, d8, d4);
            int d5 = find_opposite_node_index(d4, d8, d7);
            int d1 = find_opposite_node_index(d3, d4, d7);
            int d2 = find_opposite_node_index(d4, d5, d8);
            int d6 = find_opposite_node_index(d3, d7, d4);
            int d9 = find_opposite_node_index(d8, d5, d4);
            int d10 = find_opposite_node_index(d7, d11, d8);
            int d12 = find_opposite_node_index(d8, d11, d7);
            vector<int> v{d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12};
            face.oneRingVertices = v;
        }
    }
    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::set_one_ring_vertices_sorted] One ring vertices set." << std::endl;
    }
}