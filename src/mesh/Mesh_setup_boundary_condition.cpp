#include "mesh/Mesh.hpp"

#include <cstdint>

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
namespace
{
/// Pack an ordered vertex pair into one key, so directed edges can be counted.
inline std::uint64_t directed_edge_key(int from, int to)
{
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(from)) << 32) |
           static_cast<std::uint32_t>(to);
}
} // namespace

void Mesh::validate_volume_constraint_topology() const
{
    // No constraint means the volume is reported but never fed back into the
    // dynamics, so an open surface is merely uninformative rather than wrong.
    if (param.uVol == 0.0)
        return;

    std::vector<std::string> reasons;

    switch (param.boundaryCondition)
    {
    case BoundaryType::Periodic:
        reasons.push_back("the boundary condition is Periodic, which tiles an open sheet");
        break;
    case BoundaryType::Free:
        reasons.push_back("the boundary condition is Free, which leaves the sheet open");
        break;
    case BoundaryType::Fixed:
        break;
    }

    if (faces.empty())
    {
        reasons.push_back("the mesh has no faces");
    }

    int nGhostFaces = 0;
    for (const Face &face : faces)
    {
        if (face.isGhost)
            nGhostFaces++;
    }
    if (nGhostFaces > 0)
    {
        reasons.push_back("the mesh carries " + std::to_string(nGhostFaces) +
                          " ghost faces, so the physical region is bounded");
    }

    // Count each directed edge. On a closed, consistently oriented manifold
    // every edge is walked exactly once in each direction.
    std::unordered_map<std::uint64_t, int> directedEdgeCount;
    for (const Face &face : faces)
    {
        if (face.adjacentVertices.size() != 3)
            continue;
        for (int k = 0; k < 3; k++)
        {
            const int from = face.adjacentVertices[k];
            const int to = face.adjacentVertices[(k + 1) % 3];
            directedEdgeCount[directed_edge_key(from, to)]++;
        }
    }

    int nBoundaryEdges = 0;
    int nNonManifoldEdges = 0;
    int nMisorientedEdges = 0;
    for (const auto &entry : directedEdgeCount)
    {
        const int from = static_cast<int>(entry.first >> 32);
        const int to = static_cast<int>(entry.first & 0xFFFFFFFFu);

        if (entry.second > 1)
            nMisorientedEdges++;

        // Visit each undirected edge once, from its lower-indexed endpoint.
        if (from > to)
            continue;

        const auto opposite = directedEdgeCount.find(directed_edge_key(to, from));
        const int nIncident =
            entry.second + (opposite == directedEdgeCount.end() ? 0 : opposite->second);
        if (nIncident == 1)
            nBoundaryEdges++;
        else if (nIncident > 2)
            nNonManifoldEdges++;
    }

    if (nBoundaryEdges > 0)
    {
        reasons.push_back("the surface is open: " + std::to_string(nBoundaryEdges) +
                          " edges have a single incident face");
    }
    if (nNonManifoldEdges > 0)
    {
        reasons.push_back("the surface is not two-manifold: " + std::to_string(nNonManifoldEdges) +
                          " edges have more than two incident faces");
    }
    if (nMisorientedEdges > 0)
    {
        reasons.push_back("the surface is inconsistently oriented: " +
                          std::to_string(nMisorientedEdges) +
                          " edges are traversed the same way by both of their faces");
    }

    if (reasons.empty())
        return;

    std::ostringstream message;
    message << "[Mesh::validate_volume_constraint_topology] A volume constraint is "
            << "enabled (uvVolumeConstraint = " << param.uVol
            << ") but the mesh does not enclose a volume:";
    for (const std::string &reason : reasons)
    {
        message << "\n  - " << reason;
    }
    message << "\nSigned volume by the divergence theorem is only defined for a closed, "
            << "consistently oriented surface; on anything else the accumulated value is "
            << "not a volume and is not even independent of the coordinate origin. Set "
            << "uvVolumeConstraint = 0.0, or supply a closed mesh. "
            << "See docs/volume_functional_split.md.";
    throw std::runtime_error(message.str());
}
