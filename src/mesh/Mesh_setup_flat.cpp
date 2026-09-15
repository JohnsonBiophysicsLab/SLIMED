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
    if (param.boundaryCondition == BoundaryType::Mixed)
    {
        // The same band the global mode lays out, expressed as per-vertex
        // types -- and so, from here on, treated exactly like an imported mesh.
        std::vector<char> faceIsCopy;
        set_boundary_types_flat_mixed(faceIsCopy);
        apply_vertex_boundary_types(&faceIsCopy);
    }
    else
    {
        determine_ghost_vertices_faces();
    }
    // After the ghost flags, because an edge's flippability depends on them.
    build_edge_table();
    set_one_ring_vertices_sorted();
    validate_volume_constraint_topology();
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

    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::param]" << std::endl << param << std::endl;
    }

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


void Mesh::set_boundary_types_flat_mixed(std::vector<char> &faceIsCopy)
{
    const int nFaceX = param.nFaceX;
    const int nFaceY = param.nFaceY;
    const BoundaryType axisType[2] = {param.boundaryConditionX, param.boundaryConditionY};
    const int nFaceAlong[2] = {nFaceX, nFaceY};
    const char *axisKey[2] = {"boundaryTypeX", "boundaryTypeY"};

    for (int axis = 0; axis < 2; axis++)
    {
        if (axisType[axis] != BoundaryType::Periodic && axisType[axis] != BoundaryType::Free &&
            axisType[axis] != BoundaryType::Fixed)
        {
            throw std::runtime_error(std::string("[Mesh::set_boundary_types_flat_mixed] ") +
                                     axisKey[axis] +
                                     " must be Periodic, Free or Fixed under boundaryType = Mixed");
        }
        // Three image rings and a duplicate ring on each side, and a period
        // of at least four -- the one-ring of a face spans four rows, and a
        // period shorter than that would put a vertex and its own image in
        // the same control net.
        if (axisType[axis] == BoundaryType::Periodic && nFaceAlong[axis] < 10)
        {
            throw std::runtime_error(
                std::string("[Mesh::set_boundary_types_flat_mixed] a periodic axis needs at "
                            "least 10 faces (three image rings and a duplicate ring on each "
                            "side, and a period of at least four); ") +
                axisKey[axis] + " has " + std::to_string(nFaceAlong[axis]) +
                ". Enlarge the side or reduce lFace.");
        }
    }

    const bool periodic[2] = {axisType[0] == BoundaryType::Periodic,
                              axisType[1] == BoundaryType::Periodic};
    const int period[2] = {nFaceX - 6, nFaceY - 6};
    const int fixedRings = std::max(1, param.fixedBoundaryRings);

    // The sources occupy indices 3 .. nFace - 4 along a periodic axis; every
    // other index wraps onto them. The band this produces -- rings 0, 1, 2
    // and nFace - 2 .. nFace as pure images, rings 3 and nFace - 3 as the
    // duplicate pair -- is exactly the layout the global Periodic mode uses,
    // so its outputs line up index for index with that mode's.
    const auto wrap = [](int index, int periodLength) {
        return 3 + (((index - 3) % periodLength) + periodLength) % periodLength;
    };

    for (int j = 0; j <= nFaceY; j++)
    {
        for (int i = 0; i <= nFaceX; i++)
        {
            Vertex &vertex = vertices[(nFaceX + 1) * j + i];
            vertex.type = VertexType::Free;
            vertex.reflectiveVertexIndex = -1;
            vertex.mirrorOffset = {{0.0, 0.0, 0.0}};
            vertex.isGhost = false;

            const int sourceI = periodic[0] ? wrap(i, period[0]) : i;
            const int sourceJ = periodic[1] ? wrap(j, period[1]) : j;
            if (sourceI != i || sourceJ != j)
            {
                vertex.type = VertexType::Periodic;
                vertex.reflectiveVertexIndex = (nFaceX + 1) * sourceJ + sourceI;
                continue;
            }

            // Only a source can be clamped; an image of a clamped source is
            // held still through it.
            const bool clampedX = (axisType[0] == BoundaryType::Fixed) &&
                                  (i < fixedRings || i > nFaceX - fixedRings);
            const bool clampedY = (axisType[1] == BoundaryType::Fixed) &&
                                  (j < fixedRings || j > nFaceY - fixedRings);
            if (clampedX || clampedY)
            {
                vertex.type = VertexType::Fixed;
            }
        }
    }

    // The faces of the image band duplicate physical faces and carry no
    // energy: cells 0, 1, 2 and nFace - 3, nFace - 2, nFace - 1 along a
    // periodic axis, leaving cells 3 .. nFace - 4 -- one period -- real.
    faceIsCopy.assign(faces.size(), 0);
    for (int j = 0; j < nFaceY; j++)
    {
        for (int i = 0; i < nFaceX; i++)
        {
            const bool copy = (periodic[0] && (i < 3 || i > nFaceX - 4)) ||
                              (periodic[1] && (j < 3 || j > nFaceY - 4));
            const int faceIndex = 2 * nFaceX * j + 2 * i;
            faceIsCopy[faceIndex] = copy ? 1 : 0;
            faceIsCopy[faceIndex + 1] = copy ? 1 : 0;
        }
    }
}
