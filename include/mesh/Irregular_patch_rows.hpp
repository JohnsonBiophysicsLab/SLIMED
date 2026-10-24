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

#include <array>
#include <cstddef>
#include <vector>

#include "linalg/Linear_algebra.hpp"
#include "mesh/Subdivision_matrices.hpp"

/// Number of regular children peeled off at each subdivision step.
constexpr int kRegularChildrenPerStep = 3;

/**
 * @brief The extraordinary valences a depth list covers, in order: 4, 5, 7, 8.
 *
 * Valence 6 is regular and never reads the table, so it has no entry.
 */
constexpr int kExtraordinaryValenceCount = 4;

/// Position of a valence in an IrregularPatchDepths list, or -1 for 6 and
/// anything outside [kMinIrregularValence, kMaxIrregularValence].
int irregular_depth_index(int valence);

/// Subdivision depth per extraordinary valence, in the order 4, 5, 7, 8.
using IrregularPatchDepths = std::array<int, kExtraordinaryValenceCount>;

/**
 * @brief Default depths, chosen together with the tail correction.
 *
 * docs/irregular_patch_results.md section 9. With the tail correction these
 * leave every valence's bending error at or below 3e-5 of a fluid sheet's
 * bending energy, no worse than the legacy depths below on bending, area or
 * force, at 0.67 of their kernel samples. Without the tail correction they
 * are about 25x less accurate than the legacy depths; do not use them alone.
 */
constexpr IrregularPatchDepths kDefaultIrregularPatchDepths = {{6, 5, 5, 5}};

/**
 * @brief The depths every run used before the tail correction existed.
 *
 * Chosen by the WP6 convergence study for a bending tail below 1e-4 relative,
 * without any tail correction. The rate runs the opposite way from area:
 *
 *     N        4       5       6       7       8
 *     area   7.3x    4.9x     --     3.6x    3.4x     per level
 *     bend   2.1x    3.3x     --     4.6x    4.6x     per level
 *
 * so valence 4 needs the most levels. `irregularPatchDepths = 12, 8, 7, 7`
 * with `irregularPatchTailCorrection = false` reproduces these runs exactly.
 */
constexpr IrregularPatchDepths kLegacyIrregularPatchDepths = {{12, 8, 7, 7}};

/**
 * @brief Deepest level any default table is built to.
 *
 * The deepest entry of kLegacyIrregularPatchDepths. Also the default for a
 * Uniform build, which the convergence study sweeps directly.
 */
constexpr int kDefaultIrregularDepth = 12;

/// The deepest level a depth list may ask for: about 30 kB of rows per level.
constexpr int kMaxIrregularPatchDepth = 24;

/**
 * @brief How an irregular patch is integrated: Param::irregularPatchDepths,
 * Param::irregularPatchTailCorrection and Param::irregularPatchDepthScale.
 */
struct IrregularPatchQuadrature
{
    IrregularPatchDepths depths = kDefaultIrregularPatchDepths;
    bool tailCorrection = true;
    double depthScale = 1.0;
};

/**
 * @brief Default subdivision depth for a valence, kDefaultIrregularPatchDepths.
 *
 * Valence 6 is regular and never reads the table; it, and anything outside
 * [4, 8], returns 1.
 */
int recommended_irregular_depth(int valence);

/// The pre-tail-correction depth for a valence, kLegacyIrregularPatchDepths.
int legacy_irregular_depth(int valence);

/**
 * @brief A recommended depth with Param::irregularPatchDepthScale applied.
 *
 * Rounded to the nearest level and floored at 1: a scale small enough to ask
 * for zero levels would ask for no surface at all, since the depth-0 patch is
 * the extraordinary corner itself and has no regular child to evaluate.
 */
int scaled_irregular_depth(int recommendedDepth, double depthScale);

/**
 * @brief Quadrature weights of one level of a tail-corrected patch.
 *
 * Stam's tiling truncates an infinite sum over levels. Each quantity's terms
 * shrink geometrically, at rates fixed by the subdivision spectrum, so the
 * dropped tail is a known multiple of the last terms kept. These weights put
 * it back:
 *
 * - **Bending.** Two modes converge at different rates: the bowl (k = 0) at
 *   r_b = 1/(16 lambda^2) per level and the saddle (k = 2) at
 *   r_s = mu2^2 / lambda^2, with lambda = 3/8 + cos(2 pi/N)/4 and
 *   mu2 = 3/8 + cos(4 pi/N)/4. At N = 4, 5 the bowl is the slow one; at
 *   N = 7, 8 the saddle is. The last two levels are weighted so that both
 *   tails are summed exactly: w_{D-2} + w_{D-1} r_i = 1/(1 - r_i) for both
 *   rates. One weight on one level corrects only one mode and leaves the
 *   other wrong by 1e-3 at N = 7, 8.
 * - **Area and volume** converge at the single rate lambda^2, so the last
 *   level alone takes 1/(1 - lambda^2).
 *
 * Every weight is positive and depends only on N, D and the level, so the
 * corrected energy is still a positive-weight quadrature of the same
 * integrand: non-negative, with forces that are its exact gradient.
 *
 * With `tailCorrection` false, or for any level before the last two, both
 * weights are exactly 1.
 *
 * @param valence  4, 5, 7 or 8.
 * @param depth    Levels the patch integrates, at least 1.
 * @param level    0 <= level < depth.
 * @param[out] bendWeight  Multiplies that level's bending-energy quadrature.
 * @param[out] areaWeight  Multiplies that level's area and volume quadrature.
 */
void irregular_tail_weights(int valence, int depth, int level, bool tailCorrection,
                            double &bendWeight, double &areaWeight);

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
 * PerValence applies a per-valence depth list -- the production policy.
 * Uniform consumes the full built depth for every valence, with no tail
 * correction, which the WP6 convergence study needs: it has to be able to
 * sweep past the caps it exists to choose, and to see the plain truncated
 * sums, or it would only ever confirm them.
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
     *        the extraordinary corner is truncated. Under `PerValence` a
     *        valence consumes its default depth capped by this, and with a
     *        scale other than 1 the table is built to the deepest scaled
     *        default instead.
     * @param depthScale multiplies every valence's depth --
     *        `Param::irregularPatchDepthScale`. `Uniform` ignores it: that
     *        policy exists so the convergence study can sweep the depth
     *        itself, and a second knob on top of the one under study would
     *        only confuse it.
     *
     * Under `PerValence` this uses kDefaultIrregularPatchDepths with the tail
     * correction, which is what a Mesh built from a default Param gets.
     *
     * @throw std::invalid_argument on a non-positive depth or scale, or a
     *        shape function that is not `7 x 12`.
     */
    void build(const std::vector<Matrix> &regularShapeFunctions,
               int depth = kDefaultIrregularDepth,
               DepthPolicy policy = DepthPolicy::PerValence,
               double depthScale = 1.0);

    /**
     * @brief Build the `PerValence` table a Param asks for.
     *
     * @throw std::invalid_argument on a depth below 1 or above
     *        kMaxIrregularPatchDepth, or a non-positive scale.
     */
    void build(const std::vector<Matrix> &regularShapeFunctions,
               const IrregularPatchQuadrature &quadrature);

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
    /// The depth scale this was built with. See Param::irregularPatchDepthScale.
    double depth_scale() const { return depthScale_; }
    /// Whether the last levels carry tail weights. See irregular_tail_weights().
    bool tail_correction() const { return tailCorrection_; }
    /// The unscaled per-valence depths this was built with.
    const IrregularPatchDepths &depths() const { return depths_; }

    /**
     * @brief Depth to actually consume for this valence.
     *
     * The table is built to one depth for every valence, but callers stop at
     * the depth that valence needs -- its entry in the depth list, as scaled
     * by the depth scale this was built with.
     */
    int depth_for(int valence) const;

    /**
     * @brief The bending and area weights of one level of a valence's patch.
     *
     * irregular_tail_weights() at depth_for(valence). Both are exactly 1
     * without the tail correction, and for every level but the last two.
     */
    void level_weights(int valence, int level, double &bendWeight, double &areaWeight) const;
    int nSamples() const { return nSamples_; }

    /// The fraction of the parameter domain the truncated sliver still covers.
    double truncated_fraction() const;

    /// Total size of the stored rows, for the startup diagnostic.
    std::size_t memory_bytes() const;

private:
    int child_index(int valence, int depth, int child) const;
    void build_rows(const std::vector<Matrix> &regularShapeFunctions, int depth,
                    DepthPolicy policy, double depthScale, bool tailCorrection,
                    const IrregularPatchDepths &depths);

    int depth_ = 0;
    int nSamples_ = 0;
    DepthPolicy policy_ = DepthPolicy::PerValence;
    double depthScale_ = 1.0;
    bool tailCorrection_ = false;
    IrregularPatchDepths depths_ = kDefaultIrregularPatchDepths;
    /// Grouped by (valence, depth, child); each entry holds one Matrix per sample.
    std::vector<std::vector<Matrix>> rows_;
};
