#include "cuda/Device_mesh_layout.hpp"

#include <stdexcept>
#include <string>

#include "energy_force/Patch_kernel.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Multi_extraordinary_patch.hpp"

namespace slimed
{

void DeviceMeshLayout::build(const Mesh &mesh)
{
    nFaces_ = static_cast<int>(mesh.faces.size());
    nVertices_ = static_cast<int>(mesh.vertices.size());

    descriptors_.assign(nFaces_, FacePatchDescriptor{});
    faceValence_.assign(nFaces_, 0);
    faceSpontCurvature_.assign(nFaces_, 0.0);
    faceIsGhost_.assign(nFaces_, 0);
    faceIsBoundary_.assign(nFaces_, 0);
    faceCorners_.assign(static_cast<std::size_t>(nFaces_) * 3, 0);
    oneRingIndices_.clear();
    nMultiFaces_ = 0;

    // The prolongation table, entry for entry, so a face's patchEntry indexes
    // the snapshot exactly as it indexes the mesh's table. Copied before the
    // face loop so that the loop can check each Multi face against it.
    const MultiPatchTable &table = mesh.multiPatchTable;
    multiEntries_.assign(static_cast<std::size_t>(table.size()), DeviceMultiPatchEntry{});
    multiProlongations_.clear();
    if (table.data() != nullptr)
    {
        multiProlongations_.assign(table.data(), table.data() + table.data_count());
    }
    for (int e = 0; e < table.size(); e++)
    {
        const MultiPatchTable::Entry &entry = table.entry(e);
        DeviceMultiPatchEntry &flat = multiEntries_[e];
        flat.nControl = entry.nControl;
        for (int c = 0; c < kMultiPatchChildren; c++)
        {
            const MultiPatchTable::Child &child = entry.children[c];
            flat.childValence[c] = child.valence;
            flat.childNControl[c] = child.nControl;
            // A child is evaluated exactly as a face of its own kind would be:
            // one block when regular, Stam's tiling otherwise. Same numbers
            // the descriptors below carry for such faces.
            flat.childNChildren[c] =
                (child.valence == 6)
                    ? 1
                    : mesh.irregularRows.depth_for(child.valence) * kRegularChildrenPerStep;
            if (child.offset + static_cast<std::size_t>(child.nControl) * entry.nControl >
                    multiProlongations_.size() ||
                child.offset > static_cast<std::size_t>(0x7fffffff))
            {
                throw std::invalid_argument(
                    "[DeviceMeshLayout] prolongation entry " + std::to_string(e) + " child " +
                    std::to_string(c) + " lies outside the table's buffer");
            }
            flat.childOffset[c] = static_cast<int>(child.offset);
        }
    }

    for (int f = 0; f < nFaces_; f++)
    {
        const Face &face = mesh.faces[f];
        const int width = static_cast<int>(face.oneRingVertices.size());

        faceSpontCurvature_[f] = face.spontCurvature;
        faceIsGhost_[f] = face.isGhost ? 1 : 0;
        faceIsBoundary_[f] = face.isBoundary ? 1 : 0;
        for (int corner = 0; corner < 3; corner++)
        {
            faceCorners_[f * 3 + corner] = face.adjacentVertices[corner];
        }

        FacePatchDescriptor &descriptor = descriptors_[f];
        descriptor.oneRingOffset = static_cast<int>(oneRingIndices_.size());

        // The kind records only what the width can be evaluated as. Whether
        // the face is skipped is a separate question, and the two passes
        // answer it differently: the area pass skips ghost faces, the energy
        // pass skips boundary faces. A boundary face with a complete one-ring
        // contributes area but no energy, so folding the flags in here would
        // quietly drop that area from the constraint the energy is measured
        // against.
        //
        // The kind is read off the face, never inferred from the width: a
        // 6/5/7 face and a regular one are both 12 wide. A face with several
        // extraordinary corners resolves to a prolongation entry when it is
        // classified (Mesh::ensure_multi_patch_entries()); one that has not --
        // patchEntry still -1 -- carries no patch on the CPU either, and the
        // device mirrors that rather than inventing a different answer.
        const bool hasRegularRing = (face.patchKind == ::PatchKind::Regular);
        const bool hasIrregularRing = (face.patchKind == ::PatchKind::SingleExtraordinary);
        const bool hasMultiRing =
            (face.patchKind == ::PatchKind::MultiExtraordinary && face.patchEntry >= 0);
        if (width > slimed::kMaxControlPoints)
        {
            throw std::invalid_argument(
                "[DeviceMeshLayout] face " + std::to_string(f) + " has a one-ring of " +
                std::to_string(width) + " vertices, wider than the kernel's " +
                std::to_string(slimed::kMaxControlPoints) + "-point buffers");
        }
        if (!(hasRegularRing || hasIrregularRing || hasMultiRing))
        {
            descriptor.kind = DevicePatchKind::None;
            descriptor.nControlPoints = 0;
            descriptor.nChildren = 0;
            continue;
        }

        descriptor.nControlPoints = width;
        if (hasRegularRing)
        {
            descriptor.kind = DevicePatchKind::Regular;
            descriptor.nChildren = 1;
        }
        else if (hasIrregularRing)
        {
            const int valence = width - 6;
            descriptor.kind = DevicePatchKind::Irregular;
            faceValence_[f] = valence;
            descriptor.nChildren = mesh.irregularRows.depth_for(valence) * kRegularChildrenPerStep;
        }
        else
        {
            // The entry must exist in the snapshot and describe a patch of
            // this face's width; anything else would index the prolongation
            // buffer somewhere meaningless and produce a plausible wrong
            // force, which is the one outcome worth a throw.
            if (face.patchEntry >= static_cast<int>(multiEntries_.size()))
            {
                throw std::invalid_argument(
                    "[DeviceMeshLayout] face " + std::to_string(f) + " names prolongation entry " +
                    std::to_string(face.patchEntry) + " but the mesh's table holds " +
                    std::to_string(multiEntries_.size()) +
                    "; the face was classified against a table this layout is not seeing");
            }
            if (multiEntries_[face.patchEntry].nControl != width)
            {
                throw std::invalid_argument(
                    "[DeviceMeshLayout] face " + std::to_string(f) + " is " +
                    std::to_string(width) + " wide but its prolongation entry expects " +
                    std::to_string(multiEntries_[face.patchEntry].nControl));
            }
            descriptor.kind = DevicePatchKind::Multi;
            descriptor.nChildren = kMultiPatchChildren;
            descriptor.multiEntry = face.patchEntry;
            nMultiFaces_++;
        }

        for (int j = 0; j < width; j++)
        {
            oneRingIndices_.push_back(face.oneRingVertices[j]);
        }
    }

    // Transpose the one-ring lists: for each vertex, every slot that will hold
    // a force contribution to it. Counting first, then filling, keeps the
    // whole thing two linear passes and one allocation.
    vertexSlotOffsets_.assign(nVertices_ + 1, 0);
    for (int slot = 0; slot < static_cast<int>(oneRingIndices_.size()); slot++)
    {
        vertexSlotOffsets_[oneRingIndices_[slot] + 1]++;
    }
    for (int v = 0; v < nVertices_; v++)
    {
        vertexSlotOffsets_[v + 1] += vertexSlotOffsets_[v];
    }
    vertexSlots_.assign(oneRingIndices_.size(), 0);
    {
        std::vector<int> cursor(vertexSlotOffsets_.begin(), vertexSlotOffsets_.end() - 1);
        for (int slot = 0; slot < static_cast<int>(oneRingIndices_.size()); slot++)
        {
            vertexSlots_[cursor[oneRingIndices_[slot]]++] = slot;
        }
    }

    // The same transpose for the regularization term, whose slots are the
    // three corners of every face rather than a one-ring. Ghost faces are
    // included: energy_force_regularization() runs over every face.
    vertexCornerOffsets_.assign(nVertices_ + 1, 0);
    for (std::size_t corner = 0; corner < faceCorners_.size(); corner++)
    {
        vertexCornerOffsets_[faceCorners_[corner] + 1]++;
    }
    for (int v = 0; v < nVertices_; v++)
    {
        vertexCornerOffsets_[v + 1] += vertexCornerOffsets_[v];
    }
    vertexCorners_.assign(faceCorners_.size(), 0);
    {
        std::vector<int> cursor(vertexCornerOffsets_.begin(), vertexCornerOffsets_.end() - 1);
        for (std::size_t corner = 0; corner < faceCorners_.size(); corner++)
        {
            vertexCorners_[cursor[faceCorners_[corner]]++] = static_cast<int>(corner);
        }
    }
}

void DeviceMeshLayout::gather_coordinates(const Mesh &mesh, std::vector<double> &coords) const
{
    coords.resize(static_cast<std::size_t>(nVertices_) * 3);
    for (int v = 0; v < nVertices_; v++)
    {
        const Matrix &coord = mesh.vertices[v].coord;
        coords[v * 3 + 0] = coord.get(0, 0);
        coords[v * 3 + 1] = coord.get(1, 0);
        coords[v * 3 + 2] = coord.get(2, 0);
    }
}

void DeviceMeshLayout::gather_reference_coordinates(const Mesh &mesh,
                                                    std::vector<double> &coords) const
{
    coords.resize(static_cast<std::size_t>(nVertices_) * 3);
    for (int v = 0; v < nVertices_; v++)
    {
        const Matrix &coord = mesh.vertices[v].coordRef;
        coords[v * 3 + 0] = coord.get(0, 0);
        coords[v * 3 + 1] = coord.get(1, 0);
        coords[v * 3 + 2] = coord.get(2, 0);
    }
}

std::size_t DeviceMeshLayout::memory_bytes() const
{
    return descriptors_.size() * sizeof(FacePatchDescriptor) +
           (oneRingIndices_.size() + faceValence_.size() + faceCorners_.size() +
            vertexSlotOffsets_.size() + vertexSlots_.size() + vertexCornerOffsets_.size() +
            vertexCorners_.size()) *
               sizeof(int) +
           faceSpontCurvature_.size() * sizeof(double) + faceIsGhost_.size() +
           faceIsBoundary_.size() + multiEntries_.size() * sizeof(DeviceMultiPatchEntry) +
           multiProlongations_.size() * sizeof(double);
}

} // namespace slimed
