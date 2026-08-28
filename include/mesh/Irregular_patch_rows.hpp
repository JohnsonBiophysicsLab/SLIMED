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

/**
 * @brief Deepest recursion the table is built to.
 *
 * Set by the deepest per-valence recommendation below, so a single build
 * serves every valence. 12 levels across valences 4..8 is about 363 kB, built
 * once at startup.
 */
constexpr int kDefaultIrregularDepth = 12;

/**
 * @brief Subdivision depth chosen for a valence by the WP6 convergence study.
 *
 * Depth is chosen per valence because the convergence rate is per valence. The
 * dropped sliver is 4^-D of the parameter domain at every N, but the surface
 * inside it shrinks like the subdivision matrix's subdominant eigenvalue, and
 * bending energy -- the slowest-converging quantity, and the one most exposed
 * to the truncated corner, since Loop surfaces are C1 but not C2 there --
 * converges at a rate that runs the *opposite* way from area:
 *
 *     N        4       5       6       7       8
 *     area   7.3x    4.9x     --     3.6x    3.4x     per level
 *     bend   2.1x    3.3x     --     4.6x    4.6x     per level
 *
 * So valence 4 is the expensive one, needing roughly four more levels than
 * valence 8 for the same residual, which a single global depth cannot express.
 * Chosen for an estimated remaining bending tail below 1e-4 relative:
 *
 *     N=4  depth 12  residual 7.2e-5
 *     N=5  depth  8  residual 8.3e-5
 *     N=7  depth  7  residual 2.4e-5
 *     N=8  depth  7  residual 2.0e-5
 *
 * Valence 6 is regular and never reads the table, so its depth is irrelevant.
 */
int recommended_irregular_depth(int valence);

/**
 * @brief Immutable, thread-shared table of limit-surface rows.
 *
 * Built once at startup. It depends only on the valence and the quadrature
 * rule, never on the mesh, so unlike a per-mesh stencil cache it never needs
 * invalidation.
 */
/**
 * @brief How much of the built depth a valence should consume.
 *
 * PerValence applies recommended_irregular_depth() -- the production policy.
 * Uniform consumes the full built depth for every valence, which the WP6
 * convergence study needs: it has to be able to sweep past the caps it exists
 * to choose, or it would only ever confirm them.
 */
enum class DepthPolicy
{
    PerValence,
    Uniform
};

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
               int depth = kDefaultIrregularDepth,
               DepthPolicy policy = DepthPolicy::PerValence);

    /**
     * @brief The `7 x (N+6)` rows for one child at one sample.
     *
     * @throw std::out_of_range if any key is outside the built range.
     */
    const Matrix &rows(int valence, int depth, int child, int sample) const;

    /**
     * @brief All samples for one child, in quadrature order.
     *
     * The energy/force kernel consumes one of these per call, in the same
     * shape as `param.shapeFunctions`, so an irregular child is evaluated by
     * the identical code path as a regular face.
     *
     * @throw std::out_of_range if any key is outside the built range.
     */
    const std::vector<Matrix> &rows_for_child(int valence, int depth, int child) const;

    bool empty() const { return rows_.empty(); }
    int depth() const { return depth_; }

    /**
     * @brief Depth to actually consume for this valence.
     *
     * The table is built to one depth for every valence, but callers stop at
     * the depth that valence needs -- see recommended_irregular_depth().
     */
    int depth_for(int valence) const;
    int nSamples() const { return nSamples_; }

    /// The fraction of the parameter domain the truncated sliver still covers.
    double truncated_fraction() const;

    /// Total size of the stored rows, for the startup diagnostic.
    std::size_t memory_bytes() const;

private:
    int child_index(int valence, int depth, int child) const;

    int depth_ = 0;
    int nSamples_ = 0;
    DepthPolicy policy_ = DepthPolicy::PerValence;
    /// Grouped by (valence, depth, child); each entry holds one Matrix per sample.
    std::vector<std::vector<Matrix>> rows_;
};
