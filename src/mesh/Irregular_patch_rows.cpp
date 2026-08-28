#include "mesh/Irregular_patch_rows.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace
{
/// Valences the table covers, inclusive. 6 is regular but kept so the
/// degeneracy check has something to compare against.
constexpr int kValenceCount = kMaxIrregularValence - kMinIrregularValence + 1;

Matrix identity(int n)
{
    std::vector<std::vector<double>> values(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++)
        values[i][i] = 1.0;
    return Matrix(values);
}
} // namespace

int IrregularPatchRowTable::index(int valence, int depth, int child, int sample) const
{
    if (valence < kMinIrregularValence || valence > kMaxIrregularValence || depth < 0 ||
        depth >= depth_ || child < 0 || child >= kRegularChildrenPerStep || sample < 0 ||
        sample >= nSamples_)
    {
        throw std::out_of_range("[IrregularPatchRowTable] no rows for valence " +
                                std::to_string(valence) + ", depth " + std::to_string(depth) +
                                ", child " + std::to_string(child) + ", sample " +
                                std::to_string(sample));
    }
    const int v = valence - kMinIrregularValence;
    return ((v * depth_ + depth) * kRegularChildrenPerStep + child) * nSamples_ + sample;
}

void IrregularPatchRowTable::build(const std::vector<Matrix> &regularShapeFunctions, int depth)
{
    if (depth <= 0)
    {
        throw std::invalid_argument("[IrregularPatchRowTable] depth must be positive, got " +
                                    std::to_string(depth));
    }
    if (regularShapeFunctions.empty())
    {
        throw std::invalid_argument("[IrregularPatchRowTable] no shape functions supplied");
    }
    for (const Matrix &shapeFunction : regularShapeFunctions)
    {
        // Checked before nrow(), which dereferences unconditionally: without
        // this the validation would crash on exactly the malformed input it
        // exists to reject.
        if (shapeFunction.mat == nullptr)
        {
            throw std::invalid_argument(
                "[IrregularPatchRowTable] an unallocated shape function was supplied");
        }
        if (shapeFunction.nrow() != 7 || shapeFunction.ncol() != 12)
        {
            throw std::invalid_argument(
                "[IrregularPatchRowTable] expected 7 x 12 regular shape functions, got " +
                std::to_string(shapeFunction.nrow()) + " x " +
                std::to_string(shapeFunction.ncol()));
        }
    }

    const int nSamples = static_cast<int>(regularShapeFunctions.size());
    // Built into a local and committed at the end. Anything thrown below --
    // the child-count check inside build_subdivision_matrices(), say -- would
    // otherwise leave the table sized but partly unfilled, and every element
    // still holding a default-constructed Matrix is a null dereference waiting
    // in rows() and memory_bytes().
    //
    // Sized by construction; the elements are filled in below.
    std::vector<Matrix> built(static_cast<std::size_t>(kValenceCount) * depth *
                              kRegularChildrenPerStep * nSamples);

    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        const SubdivisionMatrices matrices = build_subdivision_matrices(valence);
        const Matrix *children[kRegularChildrenPerStep] = {&matrices.p1, &matrices.p2,
                                                           &matrices.p3};

        // S_d, the chain down to the depth-d irregular child. S_0 = I.
        Matrix chain = identity(valence + 6);
        // One subdivision step then into the smaller irregular copy.
        const Matrix descend = matrices.p4 * matrices.abar;

        for (int d = 0; d < depth; d++)
        {
            for (int c = 0; c < kRegularChildrenPerStep; c++)
            {
                // P_c * Abar * S_d, shared across samples: the topology does
                // not depend on the quadrature point.
                const Matrix childRows = (*children[c]) * matrices.abar * chain;
                for (int q = 0; q < nSamples; q++)
                {
                    const int v = valence - kMinIrregularValence;
                    const std::size_t slot =
                        ((static_cast<std::size_t>(v) * depth + d) * kRegularChildrenPerStep + c) *
                            nSamples +
                        q;
                    built[slot] = regularShapeFunctions[q] * childRows;
                }
            }
            chain = descend * chain;
        }
    }

    // Commit only once everything succeeded.
    depth_ = depth;
    nSamples_ = nSamples;
    rows_ = std::move(built);
}

const Matrix &IrregularPatchRowTable::rows(int valence, int depth, int child, int sample) const
{
    return rows_[index(valence, depth, child, sample)];
}

double IrregularPatchRowTable::truncated_fraction() const
{
    return std::pow(0.25, depth_);
}

std::size_t IrregularPatchRowTable::memory_bytes() const
{
    std::size_t total = 0;
    for (const Matrix &matrix : rows_)
    {
        total += static_cast<std::size_t>(matrix.nrow()) * matrix.ncol() * sizeof(double);
    }
    return total;
}
