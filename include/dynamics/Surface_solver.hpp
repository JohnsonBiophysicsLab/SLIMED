/**
 * @file Surface_solver.hpp
 * @brief Converting between the control net and the Loop limit surface.
 *
 * The dynamics integrates the limit surface `S`, not the control net `C`, and
 * the two are related by the limit mask `S = M C`. Every step therefore needs
 * `M` forwards, `M⁻¹` backwards, and `M⁻ᵀ` to carry the force from the control
 * points the energy is a function of onto the surface the step displaces.
 *
 * The tree has done that by building `M` as a dense N x N matrix and inverting
 * it once at setup. That has three problems, and edge flips turn all three
 * from untidy into blocking:
 *
 * - **It is not the right matrix.** Every row was written with the valence-6
 *   mask, `(1/2, 1/12, ...)`, whatever the vertex's actual valence. A fluid
 *   membrane is mostly not valence 6.
 * - **It cannot survive a flip.** A flip changes four rows of `M`; the stored
 *   inverse would have to be rebuilt from scratch, at O(N³).
 * - **It is asymmetric where it should not be**, because a vertex whose faces
 *   are all ghost gets an identity row while its neighbours' rows still refer
 *   to it. The existing code reports that asymmetry at startup and then uses
 *   `M⁻¹` in place of `M⁻ᵀ` anyway.
 *
 * ### Both directions are one symmetric system
 *
 * Write `D` for the diagonal of valences and `A` for the adjacency. The mask
 * of a free vertex is `1/2` on itself and `1/(2N)` on each neighbour, so over
 * the free vertices
 *
 * ```text
 *     M = D⁻¹ K,        K = ½ (D + A)
 * ```
 *
 * and `K` is symmetric by construction. It is also positive definite: `D + A`
 * is the signless Laplacian, positive semidefinite for any graph and singular
 * only on a bipartite one, and a triangulation has triangles. So:
 *
 * ```text
 *     M C = S     ⟺   K C = D S            solve, then done
 *     Mᵀ F_S = F_C ⟺  K y = F_C,  F_S = D y
 * ```
 *
 * Both are conjugate gradients on the same matrix. No relaxation parameter to
 * tune, no spectral bound to assume, and the residual says whether it worked.
 *
 * ### Pinned vertices
 *
 * A vertex whose adjacent faces are all ghost or boundary has no limit surface
 * of its own; the mask gives it the identity, `S_v = C_v`. Those are not
 * degrees of freedom -- the Brownian step skips them and the periodic
 * post-process overwrites them -- so they are eliminated rather than solved
 * for. Their known values move to the right-hand side, and what is left is a
 * principal submatrix of `K`, still symmetric and still positive definite.
 *
 * That elimination is also what removes the asymmetry the dense path could
 * only report: there is no asymmetric matrix left to approximate around.
 *
 * @see docs/edge_flip_plan.md section 3.7
 */

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "linalg/Linear_algebra.hpp"

class Mesh;

namespace slimed
{

/**
 * @brief The limit mask, held sparsely, with both solves.
 *
 * Rebuilt when the connectivity changes, which on a fluid membrane is every
 * step that accepts a flip. Building is O(N); the dense path's inverse was
 * O(N³) and could not be updated at all.
 */
class SurfaceSolver
{
public:
    /// Build from the mesh's current adjacency. Safe to call again.
    void build(const Mesh &mesh);

    bool empty() const { return nVertices_ == 0; }
    int nVertices() const { return nVertices_; }
    /// Vertices that are actual degrees of freedom; the rest carry the identity.
    int nFree() const { return static_cast<int>(freeOfCompact_.size()); }

    /// The connectivity this was built for. See Mesh::topologyVersion.
    long long topologyVersion = -1;

    /// Residual the conjugate gradient stops at, relative to the right-hand side.
    double tolerance = 1e-12;
    /// Hard cap, so a pathological mesh fails loudly instead of spinning.
    int maxIterations = 2000;

    /**
     * @brief The limit surface of a control net: S = M C.
     *
     * A plain sparse product; no solve involved.
     *
     * @param control N x 3, row-major in a Matrix.
     * @param surface N x 3. Overwritten.
     */
    void mesh_to_surface(const Matrix &control, Matrix &surface) const;

    /**
     * @brief The control net of a limit surface: solve M C = S.
     *
     * @param initialGuess Optional warm start; the previous step's answer is a
     *        very good one and cuts the iteration count severalfold.
     * @return Iterations taken, summed over the three coordinates.
     * @throw std::runtime_error if the solve does not converge.
     */
    int surface_to_mesh(const Matrix &surface, Matrix &control,
                        const Matrix *initialGuess = nullptr) const;

    /**
     * @brief Carry a nodal force onto the surface: solve Mᵀ F_S = F_C.
     *
     * The energy is a function of the control points, so its gradient lives
     * there; the step displaces the surface. Getting this wrong does not only
     * bias the sampled distribution -- a non-reciprocal mobility is not the
     * gradient flow of any potential and can do net work on the membrane.
     *
     * Pinned vertices are left at zero: they are not integrated.
     *
     * @return Iterations taken.
     * @throw std::runtime_error if the solve does not converge.
     */
    int nodal_force_to_surface(const Matrix &nodalForce, Matrix &surfaceForce) const;

    /// Largest |M x - b| over the last solve, for the startup diagnostic.
    double lastResidual() const { return lastResidual_; }

private:
    /// Solve K x = b over the free vertices, by conjugate gradients.
    int solve_free(const std::vector<double> &rhs, std::vector<double> &solution,
                   const std::vector<double> *initialGuess) const;

    int nVertices_ = 0;
    /// Per vertex: its valence, or 0 when pinned.
    std::vector<int> valence_;
    std::vector<char> pinned_;
    /// Compact index of each free vertex, or -1.
    std::vector<int> compactOfVertex_;
    std::vector<int> freeOfCompact_;
    /// CSR of the free-free adjacency, in compact indices.
    std::vector<int> rowStart_;
    std::vector<int> column_;
    /// Neighbours of each free vertex that are pinned, in vertex indices.
    std::vector<int> pinnedRowStart_;
    std::vector<int> pinnedColumn_;

    mutable double lastResidual_ = 0.0;
};

} // namespace slimed
