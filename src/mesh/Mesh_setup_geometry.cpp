#include "mesh/Mesh.hpp"
#include "mesh/Subdivision_matrices.hpp"

#include <cstdint>


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
    std::vector<int> commonElements;
    return faces_share_edge(face1, face2, commonElements);
}

bool Mesh::faces_share_edge(const Face& face1, const Face& face2, std::vector<int>& commonElements){
    // If face1 and face2 have two identical vertices then they share an edge.
    //
    // This was a std::set_intersection over the two adjacentVertices ranges,
    // which that algorithm requires to be sorted -- and they never are.
    // set_vertices_faces_flat() stores each triangle's corners in winding
    // order, e.g. (a, a + nFaceX + 1, a + 1) for a row-odd face, so the merge
    // walk stops at the first descent and under-reports. On the generated flat
    // sheet it missed every genuine neighbour: the only pair it accepted was a
    // face with itself, and set_adjacent_faces_of_faces() gave every face the
    // list {itself, 0, 0}.
    //
    // Three corners against three is small enough to compare directly, which
    // is exact whatever order they are stored in, and keeps face1's ordering
    // in commonElements.
    for (int vertex : face1.adjacentVertices)
    {
        if (std::find(face2.adjacentVertices.begin(), face2.adjacentVertices.end(), vertex) !=
            face2.adjacentVertices.end())
        {
            commonElements.push_back(vertex);
        }
    }

    // Check if there are at least two common vertices
    return commonElements.size() >= 2;
}

namespace
{
/// Pack an edge into one key with its endpoints in ascending order, so the two
/// faces meeting on an edge land in the same bucket however each of them winds
/// it. Mirrors directed_edge_key() in Mesh_setup_boundary_condition.cpp, which
/// deliberately does not sort because it is counting orientations.
inline std::uint64_t undirected_edge_key(int nodeA, int nodeB)
{
    const std::uint32_t low = static_cast<std::uint32_t>(nodeA < nodeB ? nodeA : nodeB);
    const std::uint32_t high = static_cast<std::uint32_t>(nodeA < nodeB ? nodeB : nodeA);
    return (static_cast<std::uint64_t>(low) << 32) | high;
}
} // namespace

/**
 * @brief Set adjacentFaces properties of faces based on the current
 * geometry of mesh.
 *
 * Two faces are adjacent when they meet on an edge, so one pass files each
 * face's three edges into a hash map keyed by that edge, and a second pass
 * reads each face's three edges back out and keeps whichever other face shares
 * the bucket. That is O(F).
 *
 * It replaces a pass that tested every ordered pair of faces with
 * faces_share_edge(). On the 99,360-face sheet (lFace = 5,
 * sideX = sideY = 1035) that was 9.9e9 predicate calls and 61.3 s of a 61.4 s
 * mesh setup -- more than the whole force evaluation the mesh was built for.
 *
 * The pairwise pass also never worked: faces_share_edge() intersected two
 * ranges that are not sorted, and on this mesh accepted only a face with
 * itself, so every face came out holding {itself, 0, 0}. See the comment
 * there. Nothing outside this file reads Face::adjacentFaces, which is why
 * that went unnoticed; its one consumer is sort_vertices_on_faces().
 */
void Mesh::set_adjacent_faces_of_faces(){
    const int nFaces = static_cast<int>(faces.size());

    // A triangle mesh has close to 3F/2 edges, so reserving for that keeps the
    // rehashing out of the pass.
    std::unordered_map<std::uint64_t, std::vector<int>> facesOnEdge;
    facesOnEdge.reserve(static_cast<std::size_t>(nFaces) * 2);

    for (int iFace = 0; iFace < nFaces; ++iFace){
        const std::vector<int>& corners = faces[iFace].adjacentVertices;
        if (corners.size() != 3){
            continue;
        }
        for (int k = 0; k < 3; ++k){
            facesOnEdge[undirected_edge_key(corners[k], corners[(k + 1) % 3])].push_back(iFace);
        }
    }

    for (int iFace = 0; iFace < nFaces; ++iFace){
        // Three slots, value-initialised, exactly as before: a face on the
        // mesh boundary has fewer than three neighbours, and
        // sort_vertices_on_faces() indexes [0], [1] and [2] unconditionally.
        faces[iFace].adjacentFaces = std::vector<int>(3);
        const std::vector<int>& corners = faces[iFace].adjacentVertices;
        if (corners.size() != 3){
            continue;
        }

        // Walk the face's own edges in winding order and pack the neighbours
        // found down from slot 0, skipping any edge on the mesh boundary. The
        // neighbours therefore stay contiguous and in winding order, which is
        // the cyclic ordering sort_vertices_on_faces() assumes when it
        // concatenates the list with itself to look for wrapped-around
        // sequences; leaving a gap mid-cycle for a missing edge would break
        // that, and a gap is indistinguishable from face 0 anyway.
        int nFound = 0;
        for (int k = 0; k < 3 && nFound < 3; ++k){
            // at(): pass one filed every edge of every three-cornered face,
            // so the key is always present.
            const std::vector<int>& onEdge =
                facesOnEdge.at(undirected_edge_key(corners[k], corners[(k + 1) % 3]));
            for (int jFace : onEdge){
                // An edge of a manifold triangle mesh carries at most one other
                // face, so this keeps everything there is to keep. Three slots
                // cannot describe a non-manifold edge in any case;
                // validate_volume_constraint_topology() is what reports those.
                if (jFace != iFace && nFound < 3){
                    faces[iFace].adjacentFaces[nFound++] = jFace;
                }
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

// To find out the one-ring vertices around face_i. A regular patch has 12; an
// irregular patch, with exactly one valence-5 corner, has 11. Anything else is
// rejected rather than silently left empty -- an empty one-ring matches neither
// arm of the dispatch in Compute_Energy_And_Force(), which would store the face
// with zero energy and zero force.
void Mesh::set_one_ring_vertices_sorted()
{
    const int nFaces = static_cast<int>(faces.size());

    // Per-face rejection reason, filled in the parallel pass and reported
    // serially afterwards: writing distinct elements is race-free, throwing out
    // of an OpenMP region is not.
    std::vector<std::string> rejection(nFaces);

#pragma omp parallel for
    for (int iFace = 0; iFace < nFaces; iFace++)
    {
        Face &face = faces[iFace];

        // Ghost faces lack the neighbouring triangles for a complete one-ring,
        // and take no part in any physical calculation.
        if (face.isGhost)
        {
            continue;
        }

        const int node0 = face.adjacentVertices[0];
        const int node1 = face.adjacentVertices[1];
        const int node2 = face.adjacentVertices[2];

        // A face touching the mesh boundary has no complete one-ring and so no
        // limit surface -- that is a property of the mesh, not an error, and it
        // is what the original comment on this function meant by "the boundary
        // faces do not have complete one-ring". Leaving oneRingVertices empty
        // here is deliberate. The defect was that *every* unmatched face took
        // this path, including genuine extraordinary interior vertices.
        if (!is_interior_vertex(node0) || !is_interior_vertex(node1) ||
            !is_interior_vertex(node2))
        {
            continue;
        }

        // Valence is the size of the one-ring. Selection below tests the same
        // quantity the classification does; the two used to disagree
        // (adjacentVertices when classifying, adjacentFaces when selecting), and
        // on any vertex where they differ d4/d7/d8 were left uninitialized and
        // read anyway.
        const int valence0 = static_cast<int>(vertices[node0].adjacentVertices.size());
        const int valence1 = static_cast<int>(vertices[node1].adjacentVertices.size());
        const int valence2 = static_cast<int>(vertices[node2].adjacentVertices.size());

        // d4 anchors the patch -- the extraordinary corner when there is one.
        // Rotating (node0, node1, node2) preserves the face winding, which
        // sort_vertices_on_faces() has already made consistent across the mesh,
        // so no per-face winding decision is taken here.
        //
        // The previous code decided winding with
        //     dot(centre, (c7 - c4) x (c8 - c4)) < 0
        // which asks whether the normal points away from the coordinate ORIGIN.
        // That is only meaningful for a surface star-shaped about the origin. It
        // is identically zero on a flat sheet in the z = 0 plane, arbitrary on a
        // membrane translated off the origin, and actively wrong on a closed
        // surface not enclosing it -- on a torus it flips exactly the inner half
        // of the faces.
        int d4 = -1;
        int d7 = -1;
        int d8 = -1;
        int valence = 0;

        auto isSupportedExtraordinary = [](int v) {
            return v >= kMinIrregularValence && v <= kMaxIrregularValence && v != 6;
        };

        if (valence0 == 6 && valence1 == 6 && valence2 == 6)
        {
            d4 = node0;
            d7 = node1;
            d8 = node2;
            valence = 6;
        }
        else if (isSupportedExtraordinary(valence0) && valence1 == 6 && valence2 == 6)
        {
            d4 = node0;
            d7 = node1;
            d8 = node2;
            valence = valence0;
        }
        else if (isSupportedExtraordinary(valence1) && valence2 == 6 && valence0 == 6)
        {
            d4 = node1;
            d7 = node2;
            d8 = node0;
            valence = valence1;
        }
        else if (isSupportedExtraordinary(valence2) && valence0 == 6 && valence1 == 6)
        {
            d4 = node2;
            d7 = node0;
            d8 = node1;
            valence = valence2;
        }
        else
        {
            // The old predicate required all three corners at valence 5 while
            // the body built a 5/6/6 patch. Those are different topology
            // classes, and only the second matches the subdivision matrices.
            rejection[iFace] = "face " + std::to_string(iFace) + " has valences (" +
                               std::to_string(valence0) + ", " + std::to_string(valence1) +
                               ", " + std::to_string(valence2) +
                               "); a patch needs exactly one corner in [" +
                               std::to_string(kMinIrregularValence) + ", " +
                               std::to_string(kMaxIrregularValence) +
                               "] with the other two at 6";
            continue;
        }

        // Walk d4's fan. fan[0] = d7, fan[1] = d8, and each step takes the
        // neighbour shared with the previous one that is not the one before
        // it -- the same opposite-node rule, applied around the corner. This
        // is what generalises the ring beyond valence 5: the old code
        // enumerated a fixed d1..d12 and could only ever produce 11 or 12
        // points.
        std::vector<int> ringInternal(valence + 6, -1);
        ringInternal[0] = d4;
        ringInternal[1] = d7;
        ringInternal[2] = d8;
        for (int k = 2; k < valence; k++)
        {
            ringInternal[1 + k] = find_opposite_node_index(d4, ringInternal[k], ringInternal[k - 1]);
        }

        const int d3 = ringInternal[valence];      // the fan closes here
        const int d5 = (valence >= 3) ? ringInternal[3] : -1;
        const int d11 = find_opposite_node_index(d7, d8, d4);
        ringInternal[valence + 1] = find_opposite_node_index(d3, d7, d4);   // d6
        ringInternal[valence + 2] = find_opposite_node_index(d8, d5, d4);   // d9
        ringInternal[valence + 3] = find_opposite_node_index(d7, d11, d8);  // d10
        ringInternal[valence + 4] = d11;
        ringInternal[valence + 5] = find_opposite_node_index(d8, d11, d7);  // d12

        // The fan must close back on d7, or the corner is not the simple
        // isolated fan the reduction assumes.
        if (d3 < 0 || find_opposite_node_index(d4, d3, ringInternal[valence - 1]) != d7)
        {
            rejection[iFace] = "face " + std::to_string(iFace) +
                               " has a corner whose fan does not close at valence " +
                               std::to_string(valence);
            continue;
        }

        // Reorder into the canonical one-ring order -- the column order of the
        // rows in IrregularPatchRowTable, shared from one place so the two
        // cannot drift apart.
        const std::vector<int> columnOf = canonical_control_order(valence);
        std::vector<int> oneRing(valence + 6, -1);
        for (int internal = 0; internal < valence + 6; internal++)
        {
            oneRing[columnOf[internal]] = ringInternal[internal];
        }

        // At valence 6 the patch is regular and carries the full d1..d12 ring;
        // the fan walk above produces the same 12 points.
        // find_opposite_node_index() returns -1 when the walk fails, which used
        // to be printed and then used as a vertex index.
        if (std::find(oneRing.begin(), oneRing.end(), -1) != oneRing.end())
        {
            rejection[iFace] = "face " + std::to_string(iFace) +
                               " has an incomplete one-ring: the two-ring walk failed at"
                               " valences (" + std::to_string(valence0) + ", " +
                               std::to_string(valence1) + ", " + std::to_string(valence2) + ")";
            continue;
        }

        face.oneRingVertices = std::move(oneRing);
    }

    report_valence_histogram();

    std::vector<std::string> rejected;
    for (int iFace = 0; iFace < nFaces; iFace++)
    {
        if (!rejection[iFace].empty())
        {
            rejected.push_back(rejection[iFace]);
        }
    }

    if (!rejected.empty())
    {
        std::ostringstream message;
        message << "[Mesh::set_one_ring_vertices_sorted] " << rejected.size()
                << " non-ghost face(s) have no supported subdivision patch:";
        const std::size_t nShown = std::min<std::size_t>(rejected.size(), 10);
        for (std::size_t i = 0; i < nShown; i++)
        {
            message << "\n  - " << rejected[i];
        }
        if (rejected.size() > nShown)
        {
            message << "\n  - ... and " << (rejected.size() - nShown) << " more";
        }
        message << "\nThese faces would otherwise carry zero energy and zero force. "
                   "See docs/irregular_patch_valence_4_to_8_plan.md.";
        throw std::runtime_error(message.str());
    }

    if (param.VERBOSE_MODE)
    {
        std::cout << "[Mesh::set_one_ring_vertices_sorted] One ring vertices set." << std::endl;
    }
}

bool Mesh::is_interior_vertex(int iVertex) const
{
    const Vertex &vertex = vertices[iVertex];
    // adjacentVertices is built as the union of the other two corners of every
    // adjacent face, so a closed fan gives equal counts and an open one -- a
    // boundary vertex -- gives exactly one more vertex than faces.
    return !vertex.adjacentFaces.empty() &&
           vertex.adjacentFaces.size() == vertex.adjacentVertices.size();
}

void Mesh::report_valence_histogram() const
{
    std::map<int, int> histogram;
    int nGhostVertices = 0;
    int nBoundaryVertices = 0;
    for (const Vertex &vertex : vertices)
    {
        if (vertex.isGhost)
        {
            nGhostVertices++;
            continue;
        }
        if (!is_interior_vertex(vertex.index))
        {
            nBoundaryVertices++;
            continue;
        }
        histogram[static_cast<int>(vertex.adjacentVertices.size())]++;
    }

    int nNonGhostFaces = 0;
    int nBoundaryFaces = 0;
    int nMultiExtraordinaryFaces = 0;
    for (const Face &face : faces)
    {
        if (face.isGhost)
            continue;
        nNonGhostFaces++;

        int nExtraordinary = 0;
        bool touchesBoundary = false;
        for (int node : face.adjacentVertices)
        {
            if (!is_interior_vertex(node))
            {
                touchesBoundary = true;
            }
            else if (vertices[node].adjacentVertices.size() != 6)
            {
                nExtraordinary++;
            }
        }
        if (touchesBoundary)
        {
            nBoundaryFaces++;
        }
        // The uniform-patch reduction assumes exactly one extraordinary corner,
        // so this count is what decides whether pre-refinement (WP7) is needed
        // on real inputs at all.
        else if (nExtraordinary > 1)
        {
            nMultiExtraordinaryFaces++;
        }
    }

    // Quiet when there is nothing to report. An all-regular interior -- which
    // every shipped workload has -- says nothing a reader needs, and printing
    // it on every mesh setup buries the cases that do matter. Anything
    // extraordinary, or any face with more than one extraordinary corner, is
    // always worth surfacing.
    const bool hasExtraordinary =
        nMultiExtraordinaryFaces > 0 ||
        std::any_of(histogram.begin(), histogram.end(),
                    [](const std::pair<const int, int> &entry) { return entry.first != 6; });
    if (!hasExtraordinary && !param.VERBOSE_MODE)
    {
        return;
    }

    std::cout << "[Mesh::valence histogram] faces: " << nNonGhostFaces
              << " non-ghost, of which " << nBoundaryFaces
              << " touch the mesh boundary (no limit surface, skipped)" << std::endl;
    std::cout << "[Mesh::valence histogram] vertices: " << nGhostVertices << " ghost, "
              << nBoundaryVertices << " boundary; interior valences:";
    for (const auto &entry : histogram)
    {
        std::cout << " N=" << entry.first << ":" << entry.second;
    }
    std::cout << std::endl;
    std::cout << "[Mesh::valence histogram] interior faces with more than one extraordinary"
                 " corner: " << nMultiExtraordinaryFaces << std::endl;
}
