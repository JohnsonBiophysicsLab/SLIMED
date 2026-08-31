#include "cuda/Device_mesh_layout.hpp"

#include <stdexcept>
#include <string>

#include "energy_force/Patch_kernel.hpp"
#include "mesh/Mesh.hpp"

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
        const bool hasRegularRing = (width == 12);
        const bool hasIrregularRing = (width >= kMinIrregularValence + 6 &&
                                       width <= kMaxIrregularValence + 6);
        if (width > slimed::kMaxControlPoints)
        {
            throw std::invalid_argument(
                "[DeviceMeshLayout] face " + std::to_string(f) + " has a one-ring of " +
                std::to_string(width) + " vertices, wider than the kernel's " +
                std::to_string(slimed::kMaxControlPoints) + "-point buffers");
        }
        if (!(hasRegularRing || hasIrregularRing))
        {
            descriptor.kind = PatchKind::None;
            descriptor.nControlPoints = 0;
            descriptor.nChildren = 0;
            continue;
        }

        descriptor.nControlPoints = width;
        if (hasRegularRing)
        {
            descriptor.kind = PatchKind::Regular;
            descriptor.nChildren = 1;
        }
        else
        {
            const int valence = width - 6;
            descriptor.kind = PatchKind::Irregular;
            faceValence_[f] = valence;
            descriptor.nChildren = mesh.irregularRows.depth_for(valence) * kRegularChildrenPerStep;
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
           faceIsBoundary_.size();
}

} // namespace slimed
