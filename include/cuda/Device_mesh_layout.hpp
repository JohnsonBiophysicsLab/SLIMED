/**
 * @file Device_mesh_layout.hpp
 * @brief The mesh, repacked into the flat arrays a GPU kernel can read.
 *
 * Mesh, Face and Vertex are pointer-chasing structures -- std::vector members,
 * gsl_matrix coordinates, one heap allocation per 3-vector. None of that can
 * cross onto a device, and none of it can be indexed by a thread block. This
 * class flattens the parts the force evaluation needs into contiguous arrays
 * with explicit index maps.
 *
 * Everything here is ordinary host C++ with no CUDA dependency, and that is
 * deliberate: the index arithmetic is where a port like this actually goes
 * wrong, so it is built and unit-tested on the host, against the production
 * force evaluation, before any of it is handed to nvcc. The .cu file adds
 * memory management and kernel launches over these same arrays and nothing
 * else.
 *
 * Two lifetimes are kept apart:
 *
 *   - Topology (build()). One-ring lists, patch widths, row-block offsets and
 *     the gather maps. Fixed for as long as the mesh connectivity is, so it is
 *     uploaded once and never touched again.
 *   - Coordinates. The only thing a line-search trial changes, and the only
 *     thing that has to be refreshed per force evaluation.
 */
#pragma once

#include <cstddef>
#include <vector>

class Face;
class Mesh;

namespace slimed
{

/// What kind of patch a face carries, and so which rows evaluate it.
enum class PatchKind : int
{
    /// No complete one-ring: the width matches no patch table, so there is no
    /// limit surface to integrate. Says nothing about whether the face is a
    /// boundary or ghost face -- those are separate flags, because the two
    /// passes skip on different ones.
    None = 0,
    /// A 12-point one-ring, evaluated once against the regular shape functions.
    Regular = 1,
    /// A valence+6 point one-ring, tiled by regular children at increasing
    /// depth and evaluated once per child.
    Irregular = 2,
};

/**
 * @brief One face's entry in the flattened mesh.
 *
 * Sixteen bytes, trivially copyable, and laid out so a thread reads its whole
 * descriptor in one go.
 */
struct FacePatchDescriptor
{
    PatchKind kind = PatchKind::None;
    /// Control points in this face's one-ring: 12, or valence + 6.
    int nControlPoints = 0;
    /// Where this face's one-ring indices start in oneRingIndices(), and
    /// equally where its force slots start in a per-slot scratch buffer.
    int oneRingOffset = 0;
    /// Number of shape-function blocks to integrate: 1 for a regular face,
    /// depth_for(valence) * kRegularChildrenPerStep for an irregular one.
    int nChildren = 0;
};

/**
 * @brief The mesh as flat arrays, plus the maps that scatter forces back.
 *
 * Force accumulation is the part of this that a GPU cannot do the way the CPU
 * does. The CPU keeps one full-length buffer per OpenMP thread and sums them;
 * with thousands of threads that is not an option, and atomicAdd would make
 * the result depend on thread scheduling -- unacceptable for a simulation
 * anyone wants to reproduce.
 *
 * Instead each face writes only its own slots, into a scratch buffer laid out
 * exactly like oneRingIndices(): one 3-vector per (face, control point). No
 * two faces ever touch the same slot, so the write needs no synchronisation at
 * all. A second pass has each vertex gather the slots that belong to it,
 * listed in a fixed order by vertexSlotOffsets() and vertexSlots(). The result
 * is deterministic run to run -- something the OpenMP path, with its
 * reduction ordering, is not.
 */
class DeviceMeshLayout
{
public:
    /**
     * @brief Flatten a mesh's topology.
     *
     * Reads connectivity and per-face constants only, never coordinates.
     * Call again if the connectivity changes -- after a refinement, say.
     *
     * @throw std::invalid_argument if a face's one-ring is wider than
     *        slimed::kMaxControlPoints, which no patch table can evaluate.
     */
    void build(const Mesh &mesh);

    bool empty() const { return descriptors_.empty(); }

    int nFaces() const { return nFaces_; }
    int nVertices() const { return nVertices_; }
    /// Total (face, control point) pairs: the length of the force scratch.
    int nSlots() const { return static_cast<int>(oneRingIndices_.size()); }

    const FacePatchDescriptor *descriptors() const { return descriptors_.data(); }
    /// Vertex index per slot, grouped by face.
    const int *oneRingIndices() const { return oneRingIndices_.data(); }
    /// Valence per face, meaningful only where kind is Irregular.
    const int *faceValence() const { return faceValence_.data(); }
    /// Spontaneous curvature per face.
    const double *faceSpontCurvature() const { return faceSpontCurvature_.data(); }
    /**
     * @brief Whether a face is a ghost, which excludes it from area and volume.
     *
     * Deliberately separate from faceIsBoundary(). calculate_element_area_volume()
     * skips ghost faces; Compute_Energy_And_Force() skips boundary faces. They
     * are not the same set, and a boundary face with a complete one-ring does
     * contribute area while contributing no energy -- so folding both into
     * PatchKind would silently drop that area.
     */
    const unsigned char *faceIsGhost() const { return faceIsGhost_.data(); }

    /// Whether a face is on the boundary, which excludes it from energy and
    /// force. See faceIsGhost().
    const unsigned char *faceIsBoundary() const { return faceIsBoundary_.data(); }
    /// The three corner vertices per face, for the regularization term.
    const int *faceCorners() const { return faceCorners_.data(); }

    /// nVertices() + 1 offsets into vertexSlots().
    const int *vertexSlotOffsets() const { return vertexSlotOffsets_.data(); }
    /// Slot indices, so that a vertex sums the scratch entries listed for it.
    const int *vertexSlots() const { return vertexSlots_.data(); }

    /// nVertices() + 1 offsets into vertexCorners().
    const int *vertexCornerOffsets() const { return vertexCornerOffsets_.data(); }
    /// Corner slots (face * 3 + corner) belonging to each vertex, for the
    /// regularization force.
    const int *vertexCorners() const { return vertexCorners_.data(); }

    /**
     * @brief Copy the mesh's current vertex coordinates into a flat buffer.
     *
     * The only per-evaluation transfer. Writes nVertices() * 3 doubles.
     */
    void gather_coordinates(const Mesh &mesh, std::vector<double> &coords) const;

    /// Copy the reference coordinates, which the regularization term needs.
    void gather_reference_coordinates(const Mesh &mesh, std::vector<double> &coords) const;

    /// Total bytes of topology held, for the startup diagnostic.
    std::size_t memory_bytes() const;

private:
    int nFaces_ = 0;
    int nVertices_ = 0;
    std::vector<FacePatchDescriptor> descriptors_;
    std::vector<int> oneRingIndices_;
    std::vector<int> faceValence_;
    std::vector<double> faceSpontCurvature_;
    std::vector<unsigned char> faceIsGhost_;
    std::vector<unsigned char> faceIsBoundary_;
    std::vector<int> faceCorners_;
    std::vector<int> vertexSlotOffsets_;
    std::vector<int> vertexSlots_;
    std::vector<int> vertexCornerOffsets_;
    std::vector<int> vertexCorners_;
};

} // namespace slimed
