#include "test_edge_flip.hpp"

/**
 * An edge flip is the move that makes the membrane a fluid rather than a
 * solid, and it is the first thing in this tree that changes connectivity
 * after setup. Everything downstream -- one-rings, the device layout, the
 * limit-surface conversion -- was written on the assumption that connectivity
 * is immutable, so the primitive has to be exactly right before any energy is
 * attached to it.
 *
 * These are the topological gates: no coordinates are read and no energy is
 * evaluated. What they pin is that a flip produces a mesh, that it produces
 * the *same* mesh when applied twice, and that every derived adjacency the
 * rest of the code reads is still true afterwards.
 *
 * @see docs/edge_flip_plan.md work package 0
 */

namespace
{
struct MeshFixture
{
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;
};

/// A curved sheet, so that nothing accidentally depends on planarity.
double bowl_height(double x, double y)
{
    return 0.06 * (x * x + y * y);
}

/**
 * @brief A triangular grid `nx` by `ny` cells wide, lifted into a bowl.
 *
 * Interior vertices are at valence 6, so a flip there produces exactly the
 * 5/5/7/7 pattern that a fluid membrane lives in, and the whole interior is
 * flippable. The same fixture the device-layout tests use.
 */
MeshFixture build_grid(int nx, int ny)
{
    MeshFixture fixture;
    const double dy = std::sqrt(3.0) / 2.0;
    for (int j = 0; j <= ny; j++)
    {
        for (int i = 0; i <= nx; i++)
        {
            const double x = i + 0.5 * j;
            const double y = j * dy;
            fixture.vertices.push_back({x, y, bowl_height(x, y)});
        }
    }
    const auto index = [nx](int i, int j) { return j * (nx + 1) + i; };
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            fixture.faces.push_back({index(i, j), index(i + 1, j), index(i, j + 1)});
            fixture.faces.push_back({index(i + 1, j), index(i + 1, j + 1), index(i, j + 1)});
        }
    }
    return fixture;
}

/**
 * @brief Build a mesh without going through set_one_ring_vertices_sorted().
 *
 * These tests are about topology alone -- no coordinates, no energy -- so they
 * assemble the adjacency directly and stop short of the patch classification
 * and the prolongation tables that a full setup would also build.
 */
void setup_topology_only(Mesh &mesh, const MeshFixture &fixture)
{
    mesh.vertices = std::vector<Vertex>(fixture.vertices.size());
    for (std::size_t i = 0; i < fixture.vertices.size(); i++)
    {
        mesh.vertices[i].index = static_cast<int>(i);
        for (int k = 0; k < 3; k++)
        {
            mesh.vertices[i].coord.set(k, 0, fixture.vertices[i][k]);
        }
    }
    mesh.faces = std::vector<Face>(fixture.faces.size());
    for (std::size_t i = 0; i < fixture.faces.size(); i++)
    {
        mesh.faces[i].index = static_cast<int>(i);
        mesh.faces[i].adjacentVertices = fixture.faces[i];
    }
    mesh.set_adjacent_faces_of_vertices_sorted();
    mesh.set_adjacent_vertices_of_vertices_sorted();
    mesh.set_adjacent_faces_of_faces();
    mesh.build_edge_table();
}

/// The mesh as a labelling-independent object: the set of triangles, each
/// written from its lowest corner so that a rotation does not change it.
std::multiset<std::vector<int>> triangle_set(const Mesh &mesh)
{
    std::multiset<std::vector<int>> triangles;
    for (const Face &face : mesh.faces)
    {
        const std::vector<int> &c = face.adjacentVertices;
        const int lowest = static_cast<int>(
            std::min_element(c.begin(), c.end()) - c.begin());
        triangles.insert({c[lowest], c[(lowest + 1) % 3], c[(lowest + 2) % 3]});
    }
    return triangles;
}

/// Vertex adjacency as sets. Order is not part of the contract: nothing reads
/// adjacentVertices or adjacentFaces in order (the fan walk goes through
/// find_opposite_node_index(), the limit mask is a sum), and a flip appends
/// where it removed.
std::vector<std::pair<std::set<int>, std::set<int>>> adjacency_sets(const Mesh &mesh)
{
    std::vector<std::pair<std::set<int>, std::set<int>>> out;
    for (const Vertex &vertex : mesh.vertices)
    {
        out.push_back({std::set<int>(vertex.adjacentVertices.begin(),
                                     vertex.adjacentVertices.end()),
                       std::set<int>(vertex.adjacentFaces.begin(), vertex.adjacentFaces.end())});
    }
    return out;
}

/// The edge table as a set of (endpoints -> the pair of faces on it).
std::set<std::vector<int>> edge_set(const Mesh &mesh)
{
    std::set<std::vector<int>> out;
    for (const MeshEdge &edge : mesh.edges)
    {
        std::vector<int> incident = {edge.face[0], edge.face[1]};
        std::sort(incident.begin(), incident.end());
        out.insert({edge.v[0], edge.v[1], incident[0], incident[1]});
    }
    return out;
}

int valence(const Mesh &mesh, int iVertex)
{
    return static_cast<int>(mesh.vertices[iVertex].adjacentVertices.size());
}

/**
 * @brief Exchange the labels of two faces throughout the mesh.
 *
 * Used to undo the one thing a double flip does not restore. Applying it to
 * the twice-flipped mesh is what turns "restored up to the two face labels"
 * into an exact comparison, and makes the claim testable rather than a caveat
 * in a comment.
 */
void relabel_two_faces(Mesh &mesh, int faceA, int faceB)
{
    std::swap(mesh.faces[faceA].adjacentVertices, mesh.faces[faceB].adjacentVertices);
    std::swap(mesh.faces[faceA].adjacentFaces, mesh.faces[faceB].adjacentFaces);
    std::swap(mesh.faces[faceA].oneRingVertices, mesh.faces[faceB].oneRingVertices);

    const auto sigma = [faceA, faceB](int f) {
        if (f == faceA) return faceB;
        if (f == faceB) return faceA;
        return f;
    };
    for (Vertex &vertex : mesh.vertices)
    {
        for (int &f : vertex.adjacentFaces)
        {
            f = sigma(f);
        }
    }
    for (MeshEdge &edge : mesh.edges)
    {
        for (int k = 0; k < 2; k++)
        {
            if (edge.face[k] >= 0)
            {
                edge.face[k] = sigma(edge.face[k]);
            }
        }
    }
}

/// A closed tetrahedron. Every edge's two opposite corners are already joined,
/// so it is the smallest mesh on which a flip would create a duplicate edge.
MeshFixture build_tetrahedron()
{
    MeshFixture fixture;
    fixture.vertices = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    // Consistently wound, normals outward.
    fixture.faces = {{0, 2, 1}, {0, 3, 2}, {0, 1, 3}, {1, 2, 3}};
    return fixture;
}

/// A deterministic, dependency-free generator, so a failure is reproducible.
struct Lcg
{
    unsigned long long state;
    explicit Lcg(unsigned long long seed) : state(seed) {}
    unsigned int next()
    {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<unsigned int>(state >> 33);
    }
    int below(int bound) { return static_cast<int>(next() % static_cast<unsigned int>(bound)); }
};
} // namespace

// ---------------------------------------------------------------------------
// The edge table
// ---------------------------------------------------------------------------

TEST(EdgeTableTest, DescribesEveryEdgeOfEveryFace)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    setup_topology_only(mesh, build_grid(5, 5));

    // Euler: V - E + F = 1 for a disk.
    const int nV = static_cast<int>(mesh.vertices.size());
    const int nE = static_cast<int>(mesh.edges.size());
    const int nF = static_cast<int>(mesh.faces.size());
    EXPECT_EQ(nV - nE + nF, 1) << "V=" << nV << " E=" << nE << " F=" << nF;

    for (const Face &face : mesh.faces)
    {
        const std::vector<int> &c = face.adjacentVertices;
        for (int k = 0; k < 3; k++)
        {
            const int iEdge = mesh.edge_between(c[k], c[(k + 1) % 3]);
            ASSERT_GE(iEdge, 0) << "face " << face.index << " has an edge not in the table";
            const MeshEdge &edge = mesh.edges[iEdge];
            EXPECT_TRUE(edge.face[0] == face.index || edge.face[1] == face.index);
            // The opposite corner recorded for this face is the third one.
            const int slot = (edge.face[0] == face.index) ? 0 : 1;
            EXPECT_EQ(edge.opposite[slot], c[(k + 2) % 3]);
        }
    }

    std::string why;
    EXPECT_TRUE(mesh.validate_manifold_topology(&why)) << why;
}

TEST(EdgeTableTest, BoundaryEdgesCarryOneFaceAndAreNotFlippable)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    setup_topology_only(mesh, build_grid(4, 4));

    int nBoundary = 0;
    for (const MeshEdge &edge : mesh.edges)
    {
        if (edge.is_boundary())
        {
            nBoundary++;
            EXPECT_LT(edge.face[1], 0);
            EXPECT_FALSE(edge.flippable) << "a boundary edge has only one face to flip";
        }
    }
    // A 4x4 grid of cells has 16 boundary edges.
    EXPECT_EQ(nBoundary, 16);
}

// ---------------------------------------------------------------------------
// The flip
// ---------------------------------------------------------------------------

TEST(EdgeFlipTest, ChangesTheFourValencesByMinusOneMinusOnePlusOnePlusOne)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    setup_topology_only(mesh, build_grid(5, 5));

    int flipped = -1;
    for (const MeshEdge &edge : mesh.edges)
    {
        if (mesh.edge_flip_is_admissible(edge.index))
        {
            flipped = edge.index;
            break;
        }
    }
    ASSERT_GE(flipped, 0) << "the grid interior should offer a flippable edge";

    const MeshEdge before = mesh.edges[flipped];
    const int a = before.v[0];
    const int b = before.v[1];
    const int t0 = before.opposite[0];
    const int t1 = before.opposite[1];
    const int valenceBefore[4] = {valence(mesh, a), valence(mesh, b), valence(mesh, t0),
                                  valence(mesh, t1)};

    mesh.flip_edge(flipped);

    EXPECT_EQ(valence(mesh, a), valenceBefore[0] - 1);
    EXPECT_EQ(valence(mesh, b), valenceBefore[1] - 1);
    EXPECT_EQ(valence(mesh, t0), valenceBefore[2] + 1);
    EXPECT_EQ(valence(mesh, t1), valenceBefore[3] + 1);

    // The edge now joins the two former opposite corners, and the old pair is
    // no longer joined at all.
    EXPECT_GE(mesh.edge_between(t0, t1), 0);
    EXPECT_LT(mesh.edge_between(a, b), 0);

    std::string why;
    EXPECT_TRUE(mesh.validate_manifold_topology(&why)) << why;
}

/**
 * @brief Every interior flip leaves a consistently wound two-manifold.
 *
 * Winding is the part most easily got wrong, and it fails silently: a reversed
 * face still has three corners and still passes an undirected edge count. The
 * validator counts *directed* edges, so a face wound the wrong way shows up as
 * the same directed edge appearing twice.
 */
TEST(EdgeFlipTest, EveryAdmissibleFlipLeavesAValidManifold)
{
    Param param;
    param.VERBOSE_MODE = false;

    int nTested = 0;
    // A fresh mesh per edge, so each test is a single flip from a known state.
    Mesh probe(param);
    setup_topology_only(probe, build_grid(5, 5));
    const int nEdges = static_cast<int>(probe.edges.size());

    for (int iEdge = 0; iEdge < nEdges; iEdge++)
    {
        Mesh mesh(param);
        setup_topology_only(mesh, build_grid(5, 5));
        if (!mesh.edge_flip_is_admissible(iEdge))
        {
            continue;
        }
        mesh.flip_edge(iEdge);
        std::string why;
        ASSERT_TRUE(mesh.validate_manifold_topology(&why))
            << "flipping edge " << iEdge << ": " << why;
        nTested++;
    }
    EXPECT_GT(nTested, 20) << "the grid should offer many admissible interior flips";
}

/**
 * @brief Flipping twice restores the mesh -- with the two face labels swapped.
 *
 * The swap is intrinsic rather than an artefact of this implementation. The
 * quadrilateral offers no canonical pairing between "the side of target0"
 * before a flip and either side after it, so every consistent assignment rule
 * composes to the exchange. What matters is that nothing else moves, which is
 * what makes a rejected Metropolis trial free to undo.
 */
TEST(EdgeFlipTest, FlippingTwiceRestoresTheMeshUpToTheTwoFaceLabels)
{
    Param param;
    param.VERBOSE_MODE = false;

    Mesh probe(param);
    setup_topology_only(probe, build_grid(5, 5));
    const int nEdges = static_cast<int>(probe.edges.size());

    int nTested = 0;
    for (int iEdge = 0; iEdge < nEdges; iEdge++)
    {
        Mesh mesh(param);
        setup_topology_only(mesh, build_grid(5, 5));
        if (!mesh.edge_flip_is_admissible(iEdge))
        {
            continue;
        }

        const auto trianglesBefore = triangle_set(mesh);
        const auto adjacencyBefore = adjacency_sets(mesh);
        const auto edgesBefore = edge_set(mesh);
        const long long versionBefore = mesh.topologyVersion;
        const int faceA = mesh.edges[iEdge].face[0];
        const int faceB = mesh.edges[iEdge].face[1];

        mesh.flip_edge(iEdge);
        ASSERT_NE(triangle_set(mesh), trianglesBefore)
            << "edge " << iEdge << " reported admissible but the flip changed nothing";

        mesh.flip_edge(iEdge);

        // The triangles come back regardless of labelling.
        EXPECT_EQ(triangle_set(mesh), trianglesBefore) << "edge " << iEdge;

        // Everything else comes back once the two labels are exchanged. Undo
        // that exchange and the match is exact -- which is the precise sense
        // in which a rejected Metropolis trial restores the mesh.
        relabel_two_faces(mesh, faceA, faceB);
        EXPECT_EQ(triangle_set(mesh), trianglesBefore) << "edge " << iEdge;
        EXPECT_EQ(adjacency_sets(mesh), adjacencyBefore) << "edge " << iEdge;
        EXPECT_EQ(edge_set(mesh), edgesBefore) << "edge " << iEdge;
        EXPECT_EQ(mesh.topologyVersion, versionBefore + 2)
            << "every connectivity change must be visible to the caches";
        nTested++;
    }
    EXPECT_GT(nTested, 20);
}

TEST(EdgeFlipTest, FlipPatchCoversEveryFaceWhoseControlNetChanges)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    setup_topology_only(mesh, build_grid(6, 6));

    int flipped = -1;
    for (const MeshEdge &edge : mesh.edges)
    {
        // An edge well inside the grid, so its whole two-ring exists.
        if (mesh.edge_flip_is_admissible(edge.index))
        {
            flipped = edge.index;
        }
    }
    ASSERT_GE(flipped, 0);

    const std::vector<int> patchBefore = mesh.flip_patch_faces(flipped);
    // Record each face's control net -- its three corners' neighbourhoods --
    // which is what its energy is actually a functional of.
    const auto controlNet = [](const Mesh &m, int iFace) {
        std::set<int> net;
        for (int corner : m.faces[iFace].adjacentVertices)
        {
            net.insert(corner);
            for (int neighbour : m.vertices[corner].adjacentVertices)
            {
                net.insert(neighbour);
            }
        }
        return net;
    };
    std::vector<std::set<int>> netsBefore(mesh.faces.size());
    for (std::size_t i = 0; i < mesh.faces.size(); i++)
    {
        netsBefore[i] = controlNet(mesh, static_cast<int>(i));
    }

    mesh.flip_edge(flipped);

    const std::vector<int> patchAfter = mesh.flip_patch_faces(flipped);
    EXPECT_EQ(patchBefore, patchAfter)
        << "the flip patch must be the same set before and after, or the undo path "
           "rebuilds something different from what the do path changed";

    // Around eighteen faces on a near-regular mesh, and certainly more than the
    // two that are retriangulated.
    EXPECT_GT(patchBefore.size(), 10u);

    for (std::size_t i = 0; i < mesh.faces.size(); i++)
    {
        const bool inPatch =
            std::find(patchBefore.begin(), patchBefore.end(), static_cast<int>(i)) !=
            patchBefore.end();
        if (inPatch)
        {
            continue;
        }
        EXPECT_EQ(controlNet(mesh, static_cast<int>(i)), netsBefore[i])
            << "face " << i << " is outside the flip patch but its control net moved";
    }
}

/**
 * @brief A long random walk through triangulation space stays a valid mesh.
 *
 * One flip being right does not make a thousand right: the failure modes that
 * matter are the ones that need an edge to have been flipped before, such as a
 * stale entry in the edge index or an opposite corner that was relinked to the
 * wrong slot.
 */
TEST(EdgeFlipTest, ManyRandomFlipsKeepTheMeshValid)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    setup_topology_only(mesh, build_grid(8, 8));

    const auto trianglesAtStart = triangle_set(mesh);
    Lcg random(20260909u);
    int nAccepted = 0;
    for (int attempt = 0; attempt < 4000; attempt++)
    {
        const int iEdge = random.below(static_cast<int>(mesh.edges.size()));
        if (!mesh.edge_flip_is_admissible(iEdge))
        {
            continue;
        }
        mesh.flip_edge(iEdge);
        nAccepted++;
        if (nAccepted % 50 == 0)
        {
            std::string why;
            ASSERT_TRUE(mesh.validate_manifold_topology(&why))
                << "after " << nAccepted << " flips: " << why;
        }
    }

    EXPECT_GT(nAccepted, 200) << "the walk should keep finding admissible flips";
    std::string why;
    EXPECT_TRUE(mesh.validate_manifold_topology(&why)) << why;
    EXPECT_NE(triangle_set(mesh), trianglesAtStart)
        << "a few hundred accepted flips should have moved the triangulation";

    // Euler characteristic is a flip invariant: a flip conserves the vertex,
    // edge and face counts individually.
    const int nV = static_cast<int>(mesh.vertices.size());
    const int nE = static_cast<int>(mesh.edges.size());
    const int nF = static_cast<int>(mesh.faces.size());
    EXPECT_EQ(nV - nE + nF, 1);

    // Every valence stayed inside the range the admission test promised.
    for (const Vertex &vertex : mesh.vertices)
    {
        if (!mesh.is_interior_vertex(vertex.index))
        {
            continue;
        }
        const int n = static_cast<int>(vertex.adjacentVertices.size());
        EXPECT_GE(n, kMinIrregularValence) << "vertex " << vertex.index;
        EXPECT_LE(n, kMaxIrregularValence) << "vertex " << vertex.index;
    }
}

/**
 * @brief Face::adjacentFaces is maintained incrementally, and must agree with
 * a full rebuild.
 *
 * A flip touches this on six faces: the two it retriangulates and the four
 * across the quadrilateral whose neighbour changed sides. Rebuilding the wrong
 * set would be invisible until sort_vertices_on_faces() -- its only consumer --
 * ran on a mesh nobody re-sorts, which is to say never.
 */
TEST(EdgeFlipTest, IncrementalFaceAdjacencyMatchesAFullRebuild)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    setup_topology_only(mesh, build_grid(7, 7));

    Lcg random(3141u);
    for (int attempt = 0; attempt < 1500; attempt++)
    {
        const int iEdge = random.below(static_cast<int>(mesh.edges.size()));
        if (mesh.edge_flip_is_admissible(iEdge))
        {
            mesh.flip_edge(iEdge);
        }
    }

    std::vector<std::set<int>> incremental;
    for (const Face &face : mesh.faces)
    {
        incremental.push_back(std::set<int>(face.adjacentFaces.begin(), face.adjacentFaces.end()));
    }

    mesh.set_adjacent_faces_of_faces();
    for (std::size_t i = 0; i < mesh.faces.size(); i++)
    {
        const std::set<int> rebuilt(mesh.faces[i].adjacentFaces.begin(),
                                    mesh.faces[i].adjacentFaces.end());
        EXPECT_EQ(incremental[i], rebuilt) << "face " << i;
    }
}

// ---------------------------------------------------------------------------
// Admission
// ---------------------------------------------------------------------------

/**
 * @brief A flip that would join two vertices that are already joined is
 * refused, because it would pinch the surface.
 *
 * The tetrahedron is the smallest mesh where this is the case on every edge:
 * the two corners opposite any edge are themselves an edge, so flipping would
 * produce a second copy of it and two faces with the same three corners. The
 * check has to come before the valence bound, or the mesh would be corrupted
 * on a configuration the valence test happens to let through.
 */
TEST(EdgeFlipAdmissionTest, RefusesAFlipThatWouldDuplicateAnExistingEdge)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    setup_topology_only(mesh, build_tetrahedron());

    std::string why;
    ASSERT_TRUE(mesh.validate_manifold_topology(&why)) << why;
    ASSERT_EQ(mesh.edges.size(), 6u);

    for (const MeshEdge &edge : mesh.edges)
    {
        ASSERT_FALSE(edge.is_boundary()) << "a tetrahedron is closed";
        // The claim the rejection makes must actually be true here.
        ASSERT_GE(mesh.edge_between(edge.opposite[0], edge.opposite[1]), 0);

        std::string reason;
        EXPECT_FALSE(mesh.edge_flip_is_admissible(edge.index, &reason));
        EXPECT_NE(reason.find("would duplicate the existing edge"), std::string::npos) << reason;
    }
}

TEST(EdgeFlipAdmissionTest, RefusesAFlipThatWouldLeaveTheSupportedValenceRange)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    setup_topology_only(mesh, build_grid(6, 6));

    Lcg random(11u);
    bool sawValenceRefusal = false;
    for (int attempt = 0; attempt < 6000 && !sawValenceRefusal; attempt++)
    {
        const int iEdge = random.below(static_cast<int>(mesh.edges.size()));
        std::string why;
        if (mesh.edge_flip_is_admissible(iEdge, &why))
        {
            mesh.flip_edge(iEdge);
            continue;
        }
        if (why.find("outside the supported range") != std::string::npos)
        {
            sawValenceRefusal = true;
        }
    }
    EXPECT_TRUE(sawValenceRefusal)
        << "flips accumulate valence, so the range bound should bind eventually";
}

TEST(EdgeFlipAdmissionTest, RefusesAFlipAcrossAnInsertionPatch)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    setup_topology_only(mesh, build_grid(5, 5));

    int candidate = -1;
    for (const MeshEdge &edge : mesh.edges)
    {
        if (mesh.edge_flip_is_admissible(edge.index))
        {
            candidate = edge.index;
            break;
        }
    }
    ASSERT_GE(candidate, 0);

    // An insertion carries per-face spontaneous curvature, and a flip would
    // hand that curvature to a different triangle.
    mesh.faces[mesh.edges[candidate].face[0]].isInsertionPatch = true;
    mesh.refresh_edge_flippability(candidate);
    std::string why;
    EXPECT_FALSE(mesh.edge_flip_is_admissible(candidate, &why));
    EXPECT_NE(why.find("not structurally flippable"), std::string::npos) << why;
}

TEST(EdgeFlipAdmissionTest, RefusesAFlipThatWouldMoveSpontaneousCurvature)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    setup_topology_only(mesh, build_grid(5, 5));

    int candidate = -1;
    for (const MeshEdge &edge : mesh.edges)
    {
        if (mesh.edge_flip_is_admissible(edge.index))
        {
            candidate = edge.index;
            break;
        }
    }
    ASSERT_GE(candidate, 0);

    mesh.faces[mesh.edges[candidate].face[0]].spontCurvature = 0.02;
    std::string why;
    EXPECT_FALSE(mesh.edge_flip_is_admissible(candidate, &why));
    EXPECT_NE(why.find("spontaneous curvature"), std::string::npos) << why;
}

TEST(EdgeFlipAdmissionTest, FrozenVerticesAreNeverFlippedAround)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    setup_topology_only(mesh, build_grid(5, 5));

    int candidate = -1;
    for (const MeshEdge &edge : mesh.edges)
    {
        if (mesh.edge_flip_is_admissible(edge.index))
        {
            candidate = edge.index;
            break;
        }
    }
    ASSERT_GE(candidate, 0);

    mesh.flipFrozenVertex.assign(mesh.vertices.size(), 0);
    mesh.flipFrozenVertex[mesh.edges[candidate].v[0]] = 1;
    mesh.refresh_edge_flippability(candidate);
    EXPECT_FALSE(mesh.edge_flip_is_admissible(candidate));
}

TEST(EdgeFlipTest, FlippingABoundaryEdgeThrowsRatherThanCorruptingTheMesh)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    setup_topology_only(mesh, build_grid(3, 3));

    int boundaryEdge = -1;
    for (const MeshEdge &edge : mesh.edges)
    {
        if (edge.is_boundary())
        {
            boundaryEdge = edge.index;
            break;
        }
    }
    ASSERT_GE(boundaryEdge, 0);
    EXPECT_THROW(mesh.flip_edge(boundaryEdge), std::invalid_argument);
    EXPECT_THROW(mesh.flip_edge(-1), std::invalid_argument);
    EXPECT_THROW(mesh.flip_edge(static_cast<int>(mesh.edges.size())), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// classify_face(): the same predicate the setup path throws on
// ---------------------------------------------------------------------------

TEST(PatchClassificationTest, AgreesWithTheSetupPathOnARegularGrid)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    setup_topology_only(mesh, build_grid(5, 5));
    mesh.determine_ghost_vertices_faces();

    int nRegular = 0;
    int nBoundary = 0;
    for (int iFace = 0; iFace < static_cast<int>(mesh.faces.size()); iFace++)
    {
        const PatchClass patch = mesh.classify_face(iFace);
        EXPECT_TRUE(patch.kind == PatchKind::Regular || patch.kind == PatchKind::Boundary)
            << "face " << iFace << " on an unflipped grid should be regular or on the boundary";
        EXPECT_TRUE(patch.why.empty()) << patch.why;
        nRegular += (patch.kind == PatchKind::Regular) ? 1 : 0;
        nBoundary += (patch.kind == PatchKind::Boundary) ? 1 : 0;
    }
    EXPECT_GT(nRegular, 0);
    EXPECT_GT(nBoundary, 0);
    // A regular face is evaluable and its one-ring is the full 12 points.
    for (int iFace = 0; iFace < static_cast<int>(mesh.faces.size()); iFace++)
    {
        std::string why;
        const bool built = mesh.build_one_ring_for_face(iFace, &why);
        EXPECT_EQ(built, mesh.classify_face(iFace).kind == PatchKind::Regular);
        if (built)
        {
            EXPECT_EQ(mesh.faces[iFace].oneRingVertices.size(), 12u);
        }
        else
        {
            EXPECT_TRUE(mesh.faces[iFace].oneRingVertices.empty());
        }
    }
}

/**
 * @brief A flip turns two regular faces into faces with three extraordinary
 * corners, and the classifier says so instead of guessing.
 *
 * This is the finding that shapes the whole plan: a flip on a hexagonal
 * lattice takes four valence-6 corners to 5, 5, 7, 7, so the two new triangles
 * each carry three extraordinary corners. The row tables admit exactly one.
 * There is no version of this move that avoids the case, which is why WP1
 * exists.
 */
TEST(PatchClassificationTest, AFlipProducesFacesWithSeveralExtraordinaryCorners)
{
    Param param;
    param.VERBOSE_MODE = false;
    Mesh mesh(param);
    setup_topology_only(mesh, build_grid(6, 6));

    int candidate = -1;
    for (const MeshEdge &edge : mesh.edges)
    {
        if (!mesh.edge_flip_is_admissible(edge.index))
        {
            continue;
        }
        // Insist on the fully regular case: all four corners at valence 6.
        const int quad[4] = {edge.v[0], edge.v[1], edge.opposite[0], edge.opposite[1]};
        bool allSix = true;
        for (int corner : quad)
        {
            allSix = allSix && (valence(mesh, corner) == 6);
        }
        if (allSix)
        {
            candidate = edge.index;
            break;
        }
    }
    ASSERT_GE(candidate, 0);

    const int face0 = mesh.edges[candidate].face[0];
    const int face1 = mesh.edges[candidate].face[1];
    ASSERT_EQ(mesh.classify_face(face0).kind, PatchKind::Regular);
    ASSERT_EQ(mesh.classify_face(face1).kind, PatchKind::Regular);

    mesh.flip_edge(candidate);

    for (int iFace : {face0, face1})
    {
        const PatchClass patch = mesh.classify_face(iFace);
        EXPECT_EQ(patch.kind, PatchKind::MultiExtraordinary)
            << "face " << iFace << " has valences (" << patch.valence[0] << ", "
            << patch.valence[1] << ", " << patch.valence[2] << ")";
        // 5, 5, 7 or 5, 7, 7 depending on which corners it kept -- never 6.
        for (int k = 0; k < 3; k++)
        {
            EXPECT_NE(patch.valence[k], 6);
        }
        // Evaluable, through the local-subdivision path. Note that this
        // fixture does not go through set_one_ring_vertices_sorted(), so the
        // one-ring is only built for the faces the flip itself rebuilt.
        EXPECT_TRUE(patch.has_evaluable_patch());
    }
}

/**
 * @brief The classifier's verdict is what the setup path throws on.
 *
 * setup_from_vertices_faces() rejects an octahedron because every face has
 * three valence-4 corners. classify_face() must reach the same conclusion
 * without the throw, or the flip sweep would admit moves the evaluator cannot
 * describe.
 */
TEST(PatchClassificationTest, AgreesWithTheSetupPathOnSolidsItAcceptsAndRejects)
{
    Param param;
    param.VERBOSE_MODE = false;

    MeshFixture octahedron;
    octahedron.vertices = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    octahedron.faces = {{0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4},
                        {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}};

    Mesh mesh(param);
    setup_topology_only(mesh, octahedron);
    for (int iFace = 0; iFace < static_cast<int>(mesh.faces.size()); iFace++)
    {
        const PatchClass patch = mesh.classify_face(iFace);
        EXPECT_EQ(patch.kind, PatchKind::MultiExtraordinary);
        EXPECT_TRUE(patch.has_evaluable_patch());
        EXPECT_TRUE(patch.why.empty()) << patch.why;
    }

    // A valence outside the supported range is a different matter: there is no
    // patch of any kind, and the classifier says so rather than guessing. A
    // tetrahedron is all valence 3.
    MeshFixture tetrahedron = build_tetrahedron();
    Param tetraParam;
    tetraParam.VERBOSE_MODE = false;
    Mesh tetraMesh(tetraParam);
    setup_topology_only(tetraMesh, tetrahedron);
    for (int iFace = 0; iFace < static_cast<int>(tetraMesh.faces.size()); iFace++)
    {
        const PatchClass patch = tetraMesh.classify_face(iFace);
        EXPECT_EQ(patch.kind, PatchKind::Inadmissible);
        EXPECT_FALSE(patch.has_evaluable_patch());
        EXPECT_NE(patch.why.find("(3, 3, 3)"), std::string::npos) << patch.why;
    }

    // And the whole-mesh path throws on it, exactly as it did before the split.
    Param throwParam;
    throwParam.VERBOSE_MODE = false;
    Mesh throwMesh(throwParam);
    EXPECT_THROW(throwMesh.setup_from_vertices_faces(tetrahedron.vertices, tetrahedron.faces),
                 std::runtime_error);
}
