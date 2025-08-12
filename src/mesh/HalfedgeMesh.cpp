/*
#include "mesh/Mesh.hpp"

void Mesh::initializeHalfedges() {
    edges.clear(); // Clear existing edges
    halfedges.clear(); // Clear existing halfedges

    // Create halfedges based on the vertices and faces
    for (const auto& face : faces) {
        const int numVertices = static_cast<int>(face.vertices.size());

        for (int i = 0; i < numVertices; ++i) {
            int vertexIndex = face.vertices[i];
            int nextVertexIndex = face.vertices[(i + 1) % numVertices];

            // Create halfedge for the current vertex
            Halfedge* currentHalfedge = createHalfedge(vertexIndex, face.index);
            edges.push_back({currentHalfedge});

            // Create twin halfedge for the next vertex
            Halfedge* twinHalfedge = createHalfedge(nextVertexIndex, -1); // -1 indicates no associated face for the twin
            currentHalfedge->twin = twinHalfedge;
        }
    }
}

Halfedge* Mesh::createHalfedge(int vertexIndex, int faceIndex) {
    Halfedge halfedge;
    halfedge.face = &faces[faceIndex];
    halfedge.faceIndex = faceIndex;
    // Set other halfedge properties based on the vertex and face information
    // ...

    halfedges.push_back(halfedge);
    return &halfedges.back();
}
*/