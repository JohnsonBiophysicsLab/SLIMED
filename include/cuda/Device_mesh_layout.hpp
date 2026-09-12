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
 *   - Topology (build()). One-ring lists, patch widths, row-block offsets,
 *     the prolongation table for faces with several extraordinary corners,
 *     and the gather maps. Fixed for as long as the mesh connectivity is, so
 *     it is uploaded once per connectivity. An edge flip changes the
 *     connectivity without changing any count; Mesh::topologyVersion is what
 *     tells Mesh::ensure_device_layout() to build this again.
 *   - Coordinates. The only thing a line-search trial changes, and the only
 *     thing that has to be refreshed per force evaluation.
 */
#pragma once

#include <cstddef>
#include <vector>

// GSL-free, so this header stays as portable as it was; it supplies
// kMultiPatchChildren.
#include "energy_force/Patch_kernel.hpp"

class Face;
class Mesh;

namespace slimed
{

/// What kind of patch a face carries, and so which rows evaluate it.
enum class DevicePatchKind : int
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
    /// More than one extraordinary corner: an N0 + N1 + N2 - 6 point one-ring,
    /// evaluated as the four children of one Loop subdivision of its own
    /// control net, each child a Regular or Irregular patch in its own right.
    /// Every edge flip produces two of these; a fluid mesh is mostly them.
    Multi = 3,
};

/**
 * @brief One face's entry in the flattened mesh.
 *
 * Twenty bytes, trivially copyable, and laid out so a thread reads its whole
 * descriptor in one go.
 */
struct FacePatchDescriptor
{
    DevicePatchKind kind = DevicePatchKind::None;
    /// Control points in this face's one-ring: 12, valence + 6, or
    /// N0 + N1 + N2 - 6 for a Multi face.
    int nControlPoints = 0;
    /// Where this face's one-ring indices start in oneRingIndices(), and
    /// equally where its force slots start in a per-slot scratch buffer.
    int oneRingOffset = 0;
    /// Number of shape-function blocks to integrate: 1 for a regular face,
    /// depth_for(valence) * kRegularChildrenPerStep for an irregular one,
    /// kMultiPatchChildren for a Multi face (each of which then integrates
    /// its own blocks, counted in DeviceMultiPatchEntry).
    int nChildren = 0;
    /// Index into multiEntries() for a Multi face; -1 otherwise. The same
    /// number as Face::patchEntry, because the snapshot keeps the table's
    /// order.
    int multiEntry = -1;
};

/**
 * @brief One valence triple's prolongations, as the device reads them.
 *
 * A flattened MultiPatchTable::Entry: the parent width and, per child, the
 * valence the existing kernels evaluate it at, its own width, how many
 * shape-function blocks that evaluation integrates, and where its
 * (child width) x nControl prolongation matrix starts in multiProlongations().
 * Child order is the table's: corner 0, corner 1, corner 2, centre.
 */
struct DeviceMultiPatchEntry
{
    /// K, the parent patch width. Equals the face's nControlPoints.
    int nControl = 0;
    /// 6 for a regular child, else the corner's valence.
    int childValence[kMultiPatchChildren] = {0, 0, 0, 0};
    /// 12, or valence + 6.
    int childNControl[kMultiPatchChildren] = {0, 0, 0, 0};
    /// Blocks the child integrates: 1 when regular, else
    /// depth_for(valence) * kRegularChildrenPerStep -- exactly what a face of
    /// that kind carries in FacePatchDescriptor::nChildren.
    int childNChildren[kMultiPatchChildren] = {0, 0, 0, 0};
    /// First double of the child's prolongation matrix in multiProlongations().
    int childOffset[kMultiPatchChildren] = {0, 0, 0, 0};
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
     *        slimed::kMaxControlPoints, which no patch table can evaluate, or
     *        if a face with several extraordinary corners does not resolve to
     *        an entry of the mesh's MultiPatchTable that matches its width.
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
     * DevicePatchKind would silently drop that area.
     */
    const unsigned char *faceIsGhost() const { return faceIsGhost_.data(); }

    /// Whether a face is on the boundary, which excludes it from energy and
    /// force. See faceIsGhost().
    const unsigned char *faceIsBoundary() const { return faceIsBoundary_.data(); }

    /**
     * @name The prolongation table, for faces with several extraordinary corners
     *
     * A snapshot of the mesh's MultiPatchTable taken at build(), entry for
     * entry, so FacePatchDescriptor::multiEntry and Face::patchEntry are the
     * same number. A snapshot rather than a pointer into the mesh so that the
     * layout is self-contained -- the device gets one copy of exactly what the
     * descriptors index -- and because the table can grow during a rejected
     * flip trial without the topology version moving; a face only ever
     * indexes an entry it was given during a flip that did move it, so a
     * rebuild keyed on the version sees every entry it needs. Both are empty
     * on a mesh with no such faces, and the kernels never index them then.
     * @{
     */
    const DeviceMultiPatchEntry *multiEntries() const
    {
        return multiEntries_.empty() ? nullptr : multiEntries_.data();
    }
    int nMultiEntries() const { return static_cast<int>(multiEntries_.size()); }
    const double *multiProlongations() const
    {
        return multiProlongations_.empty() ? nullptr : multiProlongations_.data();
    }
    std::size_t multiProlongationCount() const { return multiProlongations_.size(); }
    /// How many faces build() classified as Multi, for the startup diagnostic.
    int nMultiFaces() const { return nMultiFaces_; }
    /** @} */
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
    int nMultiFaces_ = 0;
    std::vector<DeviceMultiPatchEntry> multiEntries_;
    std::vector<double> multiProlongations_;
    std::vector<int> vertexSlotOffsets_;
    std::vector<int> vertexSlots_;
    std::vector<int> vertexCornerOffsets_;
    std::vector<int> vertexCorners_;
};

} // namespace slimed
