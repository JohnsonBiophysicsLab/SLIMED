#include "mesh/Mesh.hpp"

#include <cstdint>

void Mesh::determine_ghost_vertices_faces()
{
    // no ghost vertices for fixed BC
    if (param.boundaryCondition == BoundaryType::Fixed){
        return;
    }
    // Under the per-vertex mode the ghost flags come from the vertex types
    // (apply_vertex_boundary_types()), never from grid position.
    if (param.boundaryCondition == BoundaryType::Mixed){
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
    case BoundaryType::Fixed:
    case BoundaryType::Mixed:
        break; // returned above
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
    case BoundaryType::Fixed:
    case BoundaryType::Mixed:
        break; // returned above
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
    case BoundaryType::Mixed:
        // Images make the mesh a tile of an infinite sheet, which encloses
        // nothing. A Mixed mesh without images -- a closed surface with some
        // vertices clamped, say -- is judged on its topology below like any
        // other.
        if (!periodicImageVertices.empty())
        {
            reasons.push_back("the boundary condition is Mixed with " +
                              std::to_string(periodicImageVertices.size()) +
                              " periodic images, which tiles an open sheet");
        }
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

// ---------------------------------------------------------------------------
// BoundaryType::Mixed: the boundary condition as a property of each vertex
//
// The global modes above decide everything from grid position -- which rows
// are ghosts, which vertex duplicates which -- so they can only describe a
// rectangular sheet with the same boundary on every side. Here the mesh is
// read the other way round: each vertex says what it is, and the treatment
// follows. See docs/mixed_boundary_conditions.md.
// ---------------------------------------------------------------------------

void Mesh::apply_vertex_boundary_types(const std::vector<char> *faceIsCopy)
{
    resolve_periodic_images();
    mark_periodic_copy_faces(faceIsCopy);

    // A vertex is ghost -- scaffolding with no limit surface of its own --
    // when it is typed Ghost outright, or when it is a periodic image none of
    // whose faces carries energy. An image that does touch an energy-carrying
    // face (the duplicate ring of a periodic sheet) is not ghost: it is a
    // control point of real patches, slaved to its source rather than pinned.
    // Sources are never ghost; they are the physical vertices.
    for (Vertex &vertex : vertices)
    {
        if (vertex.type == VertexType::Ghost)
        {
            vertex.isGhost = true;
            continue;
        }
        if (!vertex.is_periodic_image())
        {
            vertex.isGhost = false;
            continue;
        }
        bool hasRealFace = false;
        for (int iFace : vertex.adjacentFaces)
        {
            if (!faces[iFace].isGhost)
            {
                hasRealFace = true;
                break;
            }
        }
        vertex.isGhost = !hasRealFace;
    }

    // Images are frozen for the flip move: a flip around one would change its
    // connectivity while its source's stayed put, and the two would stop
    // describing the same neighbourhood. Mirroring a flip onto the images is
    // the extension that would make the seam fluid; it is not done here.
    // Fixed vertices are not frozen -- what is clamped is their position, not
    // their connectivity.
    flipFrozenVertex.assign(vertices.size(), 0);
    for (int v : periodicImageVertices)
    {
        flipFrozenVertex[v] = 1;
    }

    report_vertex_boundary_summary();
}

void Mesh::resolve_periodic_images()
{
    periodicImageVertices.clear();
    const int nVertices = static_cast<int>(vertices.size());
    for (int v = 0; v < nVertices; v++)
    {
        Vertex &vertex = vertices[v];
        if (vertex.type != VertexType::Periodic)
        {
            vertex.reflectiveVertexIndex = -1;
            vertex.mirrorOffset = {{0.0, 0.0, 0.0}};
            continue;
        }

        int root = vertex.reflectiveVertexIndex;
        if (root < 0 || root >= nVertices || root == v)
        {
            throw std::invalid_argument(
                "[Mesh::resolve_periodic_images] vertex " + std::to_string(v) +
                " is a periodic image of vertex " + std::to_string(root) +
                ", which is not another vertex of this mesh (" + std::to_string(nVertices) +
                " vertices)");
        }
        // Follow a chain of images to the source at its end. A chain that
        // comes back on itself has no source, and would otherwise spin here.
        int steps = 0;
        while (vertices[root].type == VertexType::Periodic)
        {
            const int next = vertices[root].reflectiveVertexIndex;
            if (next < 0 || next >= nVertices)
            {
                throw std::invalid_argument(
                    "[Mesh::resolve_periodic_images] vertex " + std::to_string(root) +
                    " is a periodic image of vertex " + std::to_string(next) +
                    ", which is not another vertex of this mesh");
            }
            root = next;
            if (root == v || ++steps > nVertices)
            {
                throw std::invalid_argument(
                    "[Mesh::resolve_periodic_images] the mirror chain from vertex " +
                    std::to_string(v) + " never reaches a source: it is a cycle. Exactly one "
                    "vertex of each periodic family must be free or fixed, and the others "
                    "must mirror it.");
            }
        }

        vertex.reflectiveVertexIndex = root;
        for (int axis = 0; axis < 3; axis++)
        {
            vertex.mirrorOffset[axis] = vertex.coord(axis, 0) - vertices[root].coord(axis, 0);
        }
        periodicImageVertices.push_back(v);
    }
}

void Mesh::mark_periodic_copy_faces(const std::vector<char> *faceIsCopy)
{
    const int nFaces = static_cast<int>(faces.size());

    // The source a corner stands for: itself unless it is an image.
    const auto sourceOf = [this](int v) {
        return vertices[v].is_periodic_image() ? vertices[v].reflectiveVertexIndex : v;
    };

    if (faceIsCopy != nullptr)
    {
        if (static_cast<int>(faceIsCopy->size()) != nFaces)
        {
            throw std::invalid_argument(
                "[Mesh::mark_periodic_copy_faces] " + std::to_string(faceIsCopy->size()) +
                " face flags were given for a mesh of " + std::to_string(nFaces) + " faces");
        }
        for (int f = 0; f < nFaces; f++)
        {
            faces[f].isGhost = ((*faceIsCopy)[f] != 0);
        }
    }
    else
    {
        // Every physical face is a triple of sources. Faces that map to the
        // same triple are copies of one another and exactly one of them keeps
        // the energy: the one with the most source corners, so that the
        // energy sits as close to the sources as the mesh allows, and the
        // lowest index among equals -- the two faces astride a periodic seam
        // are exactly equal, and one of them has to be chosen.
        std::map<std::array<int, 3>, int> representative;
        std::vector<std::array<int, 3>> key(nFaces);
        std::vector<int> nSourceCorners(nFaces, 0);
        for (int f = 0; f < nFaces; f++)
        {
            faces[f].isGhost = false;
            const std::vector<int> &corners = faces[f].adjacentVertices;
            if (corners.size() != 3)
            {
                key[f] = {{-1, -1, -1}};
                continue;
            }
            for (int k = 0; k < 3; k++)
            {
                key[f][k] = sourceOf(corners[k]);
                nSourceCorners[f] += vertices[corners[k]].is_periodic_image() ? 0 : 1;
            }
            std::sort(key[f].begin(), key[f].end());
            const auto found = representative.find(key[f]);
            if (found == representative.end())
            {
                representative.emplace(key[f], f);
            }
            else if (nSourceCorners[f] > nSourceCorners[found->second])
            {
                found->second = f;
            }
        }
        for (int f = 0; f < nFaces; f++)
        {
            if (faces[f].adjacentVertices.size() == 3)
            {
                faces[f].isGhost = (representative.at(key[f]) != f);
            }
        }
    }

    // A vertex typed Ghost is scaffolding, so is every face it corners.
    for (Face &face : faces)
    {
        for (int corner : face.adjacentVertices)
        {
            if (vertices[corner].type == VertexType::Ghost)
            {
                face.isGhost = true;
            }
        }
    }
}

void Mesh::sync_periodic_images()
{
    for (int v : periodicImageVertices)
    {
        Vertex &image = vertices[v];
        const Vertex &source = vertices[image.reflectiveVertexIndex];
        for (int axis = 0; axis < 3; axis++)
        {
            image.coord.set(axis, 0, source.coord(axis, 0) + image.mirrorOffset[axis]);
        }
    }
}

void Mesh::fold_forces_onto_periodic_sources()
{
    // Serial on purpose: several images can share one source.
    for (int v : periodicImageVertices)
    {
        Vertex &image = vertices[v];
        vertices[image.reflectiveVertexIndex].force += image.force;
        image.force.set_all_zero();
    }
}

bool Mesh::is_independent_vertex(int iVertex) const
{
    const Vertex &vertex = vertices[iVertex];
    return !vertex.isGhost && vertex.type != VertexType::Ghost && !vertex.is_fixed() &&
           !vertex.is_periodic_image();
}

void Mesh::report_vertex_boundary_summary() const
{
    int nFree = 0;
    int nFixed = 0;
    int nImages = 0;
    int nGhostImages = 0;
    int nGhostTyped = 0;
    for (const Vertex &vertex : vertices)
    {
        switch (vertex.type)
        {
        case VertexType::Free:
            nFree++;
            break;
        case VertexType::Fixed:
            nFixed++;
            break;
        case VertexType::Periodic:
            nImages++;
            nGhostImages += vertex.isGhost ? 1 : 0;
            break;
        case VertexType::Ghost:
            nGhostTyped++;
            break;
        case VertexType::PeriodicBoundary:
            break;
        }
    }
    int nCopyFaces = 0;
    for (const Face &face : faces)
    {
        nCopyFaces += face.isGhost ? 1 : 0;
    }
    std::cout << "[Mesh::apply_vertex_boundary_types] " << vertices.size() << " vertices: " << nFree
              << " free, " << nFixed << " fixed, " << nImages << " periodic images (" << nGhostImages
              << " of them ghost, touching no energy-carrying face)";
    if (nGhostTyped > 0)
    {
        std::cout << ", " << nGhostTyped << " ghost";
    }
    std::cout << "; " << faces.size() << " faces: " << (faces.size() - nCopyFaces)
              << " carry energy, " << nCopyFaces << " are copies." << std::endl;
}
