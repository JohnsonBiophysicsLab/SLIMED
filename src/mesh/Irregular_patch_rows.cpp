#include "mesh/Irregular_patch_rows.hpp"

#include <algorithm>
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

int irregular_depth_index(int valence)
{
    switch (valence)
    {
    case 4:
        return 0;
    case 5:
        return 1;
    case 7:
        return 2;
    case 8:
        return 3;
    default:
        return -1;
    }
}

int recommended_irregular_depth(int valence)
{
    const int index = irregular_depth_index(valence);
    // Valence 6 is regular and never reads the table; anything else is
    // rejected long before it reaches here.
    return index < 0 ? 1 : kDefaultIrregularPatchDepths[index];
}

int legacy_irregular_depth(int valence)
{
    const int index = irregular_depth_index(valence);
    return index < 0 ? 1 : kLegacyIrregularPatchDepths[index];
}

int scaled_irregular_depth(int recommendedDepth, double depthScale)
{
    if (!(depthScale > 0.0) || depthScale == 1.0)
    {
        return recommendedDepth;
    }
    const long scaled = std::lround(static_cast<double>(recommendedDepth) * depthScale);
    return static_cast<int>(std::max<long>(1, scaled));
}

void irregular_tail_weights(int valence, int depth, int level, bool tailCorrection,
                            double &bendWeight, double &areaWeight)
{
    bendWeight = 1.0;
    areaWeight = 1.0;
    if (!tailCorrection || irregular_depth_index(valence) < 0 || depth < 1 || level < depth - 2 ||
        level >= depth)
    {
        return;
    }

    const double pi = std::acos(-1.0);
    const double lambda = 3.0 / 8.0 + std::cos(2.0 * pi / valence) / 4.0;
    const double mu2 = 3.0 / 8.0 + std::cos(4.0 * pi / valence) / 4.0;
    const double lambda2 = lambda * lambda;
    const double rBowl = 1.0 / (16.0 * lambda2);
    const double rSaddle = mu2 * mu2 / lambda2;

    if (level == depth - 1)
    {
        areaWeight = 1.0 / (1.0 - lambda2);
    }

    if (depth == 1)
    {
        // One level cannot carry two corrections; take the slower mode's.
        bendWeight = 1.0 / (1.0 - std::max(rBowl, rSaddle));
        return;
    }
    // w_{D-2} + w_{D-1} r = 1/(1 - r) for r = rBowl and r = rSaddle.
    const double gBowl = 1.0 / (1.0 - rBowl);
    const double gSaddle = 1.0 / (1.0 - rSaddle);
    const double last = (gBowl - gSaddle) / (rBowl - rSaddle);
    bendWeight = (level == depth - 1) ? last : gBowl - last * rBowl;
}

int IrregularPatchRowTable::depth_for(int valence) const
{
    const int index = irregular_depth_index(valence);
    if (policy_ == DepthPolicy::Uniform || index < 0)
    {
        return depth_;
    }
    return std::min(depth_, scaled_irregular_depth(depths_[index], depthScale_));
}

void IrregularPatchRowTable::level_weights(int valence, int level, double &bendWeight,
                                           double &areaWeight) const
{
    irregular_tail_weights(valence, depth_for(valence), level,
                           policy_ == DepthPolicy::PerValence && tailCorrection_, bendWeight,
                           areaWeight);
}

int IrregularPatchRowTable::child_index(int valence, int depth, int child) const
{
    if (valence < kMinIrregularValence || valence > kMaxIrregularValence || depth < 0 ||
        depth >= depth_ || child < 0 || child >= kRegularChildrenPerStep)
    {
        throw std::out_of_range("[IrregularPatchRowTable] no rows for valence " +
                                std::to_string(valence) + ", depth " + std::to_string(depth) +
                                ", child " + std::to_string(child));
    }
    const int v = valence - kMinIrregularValence;
    return (v * depth_ + depth) * kRegularChildrenPerStep + child;
}

void IrregularPatchRowTable::build(const std::vector<Matrix> &regularShapeFunctions, int depth,
                                   DepthPolicy policy, double depthScale)
{
    if (!(depthScale > 0.0))
    {
        throw std::invalid_argument(
            "[IrregularPatchRowTable] irregularPatchDepthScale must be positive, got " +
            std::to_string(depthScale));
    }
    if (policy == DepthPolicy::PerValence && depthScale != 1.0)
    {
        // A scaled table is built to exactly what it hands out.
        IrregularPatchQuadrature quadrature;
        quadrature.depthScale = depthScale;
        build(regularShapeFunctions, quadrature);
        return;
    }
    // Built to the depth asked for; PerValence then consumes the default
    // depths, capped by it, with the tail correction.
    build_rows(regularShapeFunctions, depth, policy, depthScale,
               policy == DepthPolicy::PerValence, kDefaultIrregularPatchDepths);
}

void IrregularPatchRowTable::build(const std::vector<Matrix> &regularShapeFunctions,
                                   const IrregularPatchQuadrature &quadrature)
{
    if (!(quadrature.depthScale > 0.0))
    {
        throw std::invalid_argument(
            "[IrregularPatchRowTable] irregularPatchDepthScale must be positive, got " +
            std::to_string(quadrature.depthScale));
    }
    // The built depth is only ever consumed up to the deepest scaled entry, so
    // build exactly that much: shallower depths cost less to build as well as
    // to evaluate, and deeper ones have the levels to hand out.
    int needed = 1;
    for (int index = 0; index < kExtraordinaryValenceCount; index++)
    {
        const int depth = quadrature.depths[index];
        if (depth < 1 || depth > kMaxIrregularPatchDepth)
        {
            throw std::invalid_argument(
                "[IrregularPatchRowTable] irregularPatchDepths entries must lie in [1, " +
                std::to_string(kMaxIrregularPatchDepth) + "], got " + std::to_string(depth));
        }
        needed = std::max(needed, scaled_irregular_depth(depth, quadrature.depthScale));
    }
    build_rows(regularShapeFunctions, needed, DepthPolicy::PerValence, quadrature.depthScale,
               quadrature.tailCorrection, quadrature.depths);
}

void IrregularPatchRowTable::build_rows(const std::vector<Matrix> &regularShapeFunctions,
                                        int depth, DepthPolicy policy, double depthScale,
                                        bool tailCorrection, const IrregularPatchDepths &depths)
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
    std::vector<std::vector<Matrix>> built(
        static_cast<std::size_t>(kValenceCount) * depth * kRegularChildrenPerStep,
        std::vector<Matrix>(nSamples));

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
                const int v = valence - kMinIrregularValence;
                const std::size_t slot =
                    (static_cast<std::size_t>(v) * depth + d) * kRegularChildrenPerStep + c;
                for (int q = 0; q < nSamples; q++)
                {
                    built[slot][q] = regularShapeFunctions[q] * childRows;
                }
            }
            chain = descend * chain;
        }
    }

    // Commit only once everything succeeded.
    depth_ = depth;
    nSamples_ = nSamples;
    policy_ = policy;
    depthScale_ = depthScale;
    tailCorrection_ = tailCorrection;
    depths_ = depths;
    rows_ = std::move(built);
}

const Matrix &IrregularPatchRowTable::rows(int valence, int depth, int child, int sample) const
{
    const std::vector<Matrix> &samples = rows_for_child(valence, depth, child);
    if (sample < 0 || sample >= nSamples_)
    {
        throw std::out_of_range("[IrregularPatchRowTable] no sample " + std::to_string(sample));
    }
    return samples[sample];
}

const std::vector<Matrix> &IrregularPatchRowTable::rows_for_child(int valence, int depth,
                                                                  int child) const
{
    return rows_[child_index(valence, depth, child)];
}

double IrregularPatchRowTable::truncated_fraction() const
{
    return std::pow(0.25, depth_);
}

std::size_t IrregularPatchRowTable::memory_bytes() const
{
    std::size_t total = 0;
    for (const std::vector<Matrix> &samples : rows_)
    {
        for (const Matrix &matrix : samples)
        {
            total += static_cast<std::size_t>(matrix.nrow()) * matrix.ncol() * sizeof(double);
        }
    }
    return total;
}
