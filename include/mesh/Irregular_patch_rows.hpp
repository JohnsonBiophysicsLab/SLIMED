/**
 * @file Irregular_patch_rows.hpp
 * @brief Collapse the irregular-patch subdivision recursion into a table of
 * limit-surface rows, one per (valence, depth, child, quadrature sample).
 *
 * The entire recursion is pure topology: it depends only on the valence N,
 * never on vertex coordinates. So it can be collapsed once at startup into a
 * set of `7 x (N+6)` matrices -- exactly the object the regular path already
 * consumes as `param.shapeFunctions`, only wider or narrower than 12.
 *
 * Writing `S_d = (P4 * Abar)^d` for the chain down to the depth-`d` irregular
 * child, and `W_q` for the regular `7 x 12` shape function at sample `q`:
 *
 *     R[N,d,c,q]  =  W_q * P_c * Abar * S_d          // shape 7 x (N+6)
 *
 *     geometry    G = R * X                          // X = (N+6) x 3 controls
 *     force       f = R^T * (dE/dG)                  // the regular transpose
 *
 * After this collapse an irregular face costs the same as a regular one and
 * nothing downstream needs to know a face was irregular.
 *
 * @warning There is deliberately no `4^-d` weight rescaling. The derivative
 * rows returned at depth `d` are with respect to the *child's* parameters, so
 * `|a_1 x a_2|` already carries the shrunken metric and the children tile the
 * surface exactly once. Applying an extra `4^-d` would double-count the
 * Jacobian. IrregularPatchRowTableTest.PlanarRegularPatchTilesTheDomainExactly
 * is what holds that line.
 *
 * @note WP3 of docs/irregular_patch_valence_4_to_8_plan.md.
 */
#pragma once

#include <cstddef>
#include <vector>

#include "linalg/Linear_algebra.hpp"
#include "mesh/Subdivision_matrices.hpp"

/// Number of regular children peeled off at each subdivision step.
constexpr int kRegularChildrenPerStep = 3;

/// Default recursion depth. WP6 replaces this with a per-valence default.
constexpr int kDefaultIrregularDepth = 6;

/**
 * @brief Immutable, thread-shared table of limit-surface rows.
 *
 * Built once at startup. It depends only on the valence and the quadrature
 * rule, never on the mesh, so unlike a per-mesh stencil cache it never needs
 * invalidation.
 */
class IrregularPatchRowTable
{
public:
    /**
     * @brief Build the table for every supported valence.
     *
     * @param regularShapeFunctions the `W_q`, one `7 x 12` matrix per
     *        quadrature sample -- the same `param.shapeFunctions` the regular
     *        kernel uses.
     * @param depth how many times to recurse. The children tile all but
     *        `4^-depth` of the parameter domain; the remaining sliver around
     *        the extraordinary corner is truncated.
     *
     * @throw std::invalid_argument on a non-positive depth or a shape function
     *        that is not `7 x 12`.
     */
    void build(const std::vector<Matrix> &regularShapeFunctions,
               int depth = kDefaultIrregularDepth);

    /**
     * @brief The `7 x (N+6)` rows for one child at one sample.
     *
     * @throw std::out_of_range if any key is outside the built range.
     */
    const Matrix &rows(int valence, int depth, int child, int sample) const;

    bool empty() const { return rows_.empty(); }
    int depth() const { return depth_; }
    int nSamples() const { return nSamples_; }

    /// The fraction of the parameter domain the truncated sliver still covers.
    double truncated_fraction() const;

    /// Total size of the stored rows, for the startup diagnostic.
    std::size_t memory_bytes() const;

private:
    int index(int valence, int depth, int child, int sample) const;

    int depth_ = 0;
    int nSamples_ = 0;
    std::vector<Matrix> rows_;
};
