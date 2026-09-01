#include "energy_force/Patch_rows_flat.hpp"

#include <stdexcept>
#include <string>

#include "energy_force/Patch_kernel.hpp"
#include "mesh/Subdivision_matrices.hpp"

namespace
{
constexpr int kValenceCount = kMaxIrregularValence - kMinIrregularValence + 1;
constexpr int kRegularControlPoints = 12;
constexpr std::size_t kNoSlot = static_cast<std::size_t>(-1);

/// Copy one 7 x nCtrl block into dest in row-major order.
void append_block(const Matrix &block, int nCtrl, double *dest)
{
    for (int row = 0; row < slimed::kShapeRows; row++)
    {
        for (int col = 0; col < nCtrl; col++)
        {
            dest[row * nCtrl + col] = block.get(row, col);
        }
    }
}
} // namespace

void PatchRowsFlat::build(const std::vector<Matrix> &regularShapeFunctions,
                          const IrregularPatchRowTable &irregularRows,
                          const Matrix &gaussQuadratureCoeff)
{
    if (regularShapeFunctions.empty())
    {
        throw std::invalid_argument("[PatchRowsFlat] no shape functions supplied");
    }
    for (const Matrix &shapeFunction : regularShapeFunctions)
    {
        // mat is checked before nrow(), which dereferences it unconditionally.
        if (shapeFunction.mat == nullptr || shapeFunction.nrow() != slimed::kShapeRows ||
            shapeFunction.ncol() != kRegularControlPoints)
        {
            throw std::invalid_argument(
                "[PatchRowsFlat] expected 7 x 12 regular shape functions");
        }
    }

    nSamples_ = static_cast<int>(regularShapeFunctions.size());
    if (gaussQuadratureCoeff.mat == nullptr || gaussQuadratureCoeff.nrow() != nSamples_)
    {
        throw std::invalid_argument(
            "[PatchRowsFlat] expected one quadrature weight per shape function, got " +
            std::to_string(gaussQuadratureCoeff.mat == nullptr ? 0 : gaussQuadratureCoeff.nrow()) +
            " for " + std::to_string(nSamples_) + " points");
    }

    gaussCoeff_.resize(nSamples_);
    for (int sample = 0; sample < nSamples_; sample++)
    {
        gaussCoeff_[sample] = gaussQuadratureCoeff.get(sample, 0);
    }

    const std::size_t regularBlock =
        static_cast<std::size_t>(slimed::kShapeRows) * kRegularControlPoints;
    regular_.assign(regularBlock * nSamples_, 0.0);
    for (int sample = 0; sample < nSamples_; sample++)
    {
        append_block(regularShapeFunctions[sample], kRegularControlPoints,
                     regular_.data() + sample * regularBlock);
    }

    // A mesh whose faces are all regular -- every flat sheet away from its
    // boundary -- never builds the irregular table, and asking it for rows
    // would throw. Leaving the offsets empty makes child() report that
    // honestly instead.
    irregular_.clear();
    offsets_.clear();
    depth_ = 0;
    if (irregularRows.empty())
    {
        return;
    }

    depth_ = irregularRows.depth();
    offsets_.assign(static_cast<std::size_t>(kValenceCount) * depth_ * kRegularChildrenPerStep,
                    kNoSlot);

    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        const int nCtrl = valence + 6;
        const std::size_t blockSize = static_cast<std::size_t>(slimed::kShapeRows) * nCtrl;
        // Only the depths this valence actually integrates over are packed;
        // depth_for() truncates deep valences and those blocks are never read.
        for (int depth = 0; depth < irregularRows.depth_for(valence); depth++)
        {
            for (int child = 0; child < kRegularChildrenPerStep; child++)
            {
                const std::vector<Matrix> &rows = irregularRows.rows_for_child(valence, depth, child);
                const std::size_t offset = irregular_.size();
                irregular_.resize(offset + blockSize * nSamples_, 0.0);
                for (int sample = 0; sample < nSamples_; sample++)
                {
                    append_block(rows[sample], nCtrl,
                                 irregular_.data() + offset + sample * blockSize);
                }
                offsets_[slot_index(valence, depth, child)] = offset;
            }
        }
    }
}

int PatchRowsFlat::slot_index(int valence, int depth, int childIndex) const
{
    if (valence < kMinIrregularValence || valence > kMaxIrregularValence || depth < 0 ||
        depth >= depth_ || childIndex < 0 || childIndex >= kRegularChildrenPerStep)
    {
        throw std::out_of_range("[PatchRowsFlat] no rows for valence " + std::to_string(valence) +
                                ", depth " + std::to_string(depth) + ", child " +
                                std::to_string(childIndex));
    }
    const int v = valence - kMinIrregularValence;
    return (v * depth_ + depth) * kRegularChildrenPerStep + childIndex;
}

const double *PatchRowsFlat::child(int valence, int depth, int childIndex) const
{
    const std::size_t offset = offsets_.empty() ? kNoSlot : offsets_[slot_index(valence, depth, childIndex)];
    if (offset == kNoSlot)
    {
        throw std::out_of_range("[PatchRowsFlat] irregular rows were not built for valence " +
                                std::to_string(valence) + ", depth " + std::to_string(depth));
    }
    return irregular_.data() + offset;
}

std::size_t PatchRowsFlat::memory_bytes() const
{
    return (regular_.size() + irregular_.size() + gaussCoeff_.size()) * sizeof(double);
}
