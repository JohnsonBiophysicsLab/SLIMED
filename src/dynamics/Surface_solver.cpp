/**
 * @file Surface_solver.cpp
 * @brief Sparse, valence-aware conversion between the control net and the
 * limit surface.
 *
 * @see include/dynamics/Surface_solver.hpp for what this computes and why both
 * directions reduce to one symmetric positive definite system.
 */

#include "dynamics/Surface_solver.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

#include "mesh/Mesh.hpp"

namespace slimed
{

void SurfaceSolver::build(const Mesh &mesh)
{
    nVertices_ = static_cast<int>(mesh.vertices.size());
    valence_.assign(nVertices_, 0);
    pinned_.assign(nVertices_, 0);
    compactOfVertex_.assign(nVertices_, -1);
    freeOfCompact_.clear();
    rowStart_.clear();
    column_.clear();
    pinnedRowStart_.clear();
    pinnedColumn_.clear();
    topologyVersion = mesh.topologyVersion;

    // A vertex with no adjacent face that contributes to the physical surface
    // has no limit point of its own, and the mask gives it the identity. That
    // is the same test the dense path applies, kept identical so the two agree
    // about which vertices are degrees of freedom.
    for (int v = 0; v < nVertices_; v++)
    {
        bool hasRealFace = false;
        for (int iFace : mesh.vertices[v].adjacentFaces)
        {
            const Face &face = mesh.faces[iFace];
            if (!face.isGhost && !face.isBoundary)
            {
                hasRealFace = true;
                break;
            }
        }
        const int valence = static_cast<int>(mesh.vertices[v].adjacentVertices.size());
        pinned_[v] = (!hasRealFace || valence == 0) ? 1 : 0;
        valence_[v] = pinned_[v] ? 0 : valence;
    }

    for (int v = 0; v < nVertices_; v++)
    {
        if (!pinned_[v])
        {
            compactOfVertex_[v] = static_cast<int>(freeOfCompact_.size());
            freeOfCompact_.push_back(v);
        }
    }

    const int nFreeVertices = static_cast<int>(freeOfCompact_.size());
    rowStart_.reserve(nFreeVertices + 1);
    pinnedRowStart_.reserve(nFreeVertices + 1);
    rowStart_.push_back(0);
    pinnedRowStart_.push_back(0);
    for (int compact = 0; compact < nFreeVertices; compact++)
    {
        const int v = freeOfCompact_[compact];
        for (int neighbour : mesh.vertices[v].adjacentVertices)
        {
            if (pinned_[neighbour])
            {
                pinnedColumn_.push_back(neighbour);
            }
            else
            {
                column_.push_back(compactOfVertex_[neighbour]);
            }
        }
        rowStart_.push_back(static_cast<int>(column_.size()));
        pinnedRowStart_.push_back(static_cast<int>(pinnedColumn_.size()));
    }
}

void SurfaceSolver::mesh_to_surface(const Matrix &control, Matrix &surface) const
{
    for (int v = 0; v < nVertices_; v++)
    {
        if (pinned_[v])
        {
            for (int axis = 0; axis < 3; axis++)
            {
                surface.set(v, axis, control(v, axis));
            }
            continue;
        }
        const int compact = compactOfVertex_[v];
        const double neighbourWeight = 0.5 / valence_[v];
        for (int axis = 0; axis < 3; axis++)
        {
            double sum = 0.5 * control(v, axis);
            for (int k = rowStart_[compact]; k < rowStart_[compact + 1]; k++)
            {
                sum += neighbourWeight * control(freeOfCompact_[column_[k]], axis);
            }
            for (int k = pinnedRowStart_[compact]; k < pinnedRowStart_[compact + 1]; k++)
            {
                sum += neighbourWeight * control(pinnedColumn_[k], axis);
            }
            surface.set(v, axis, sum);
        }
    }
}

int SurfaceSolver::solve_free(const std::vector<double> &rhs, std::vector<double> &solution,
                              const std::vector<double> *initialGuess) const
{
    const int n = static_cast<int>(freeOfCompact_.size());
    solution.assign(n, 0.0);
    if (n == 0)
    {
        return 0;
    }
    if (initialGuess != nullptr && static_cast<int>(initialGuess->size()) == n)
    {
        solution = *initialGuess;
    }

    // K x, with K = 1/2 (D + A) restricted to the free vertices. The 1/2 is
    // carried here rather than in the stored structure, which holds only the
    // adjacency.
    const auto applyK = [&](const std::vector<double> &x, std::vector<double> &out) {
        for (int i = 0; i < n; i++)
        {
            double sum = valence_[freeOfCompact_[i]] * x[i];
            for (int k = rowStart_[i]; k < rowStart_[i + 1]; k++)
            {
                sum += x[column_[k]];
            }
            out[i] = 0.5 * sum;
        }
    };

    std::vector<double> residual(n);
    std::vector<double> direction(n);
    std::vector<double> scratch(n);

    applyK(solution, scratch);
    double residualNorm = 0.0;
    double rightHandNorm = 0.0;
    for (int i = 0; i < n; i++)
    {
        residual[i] = rhs[i] - scratch[i];
        direction[i] = residual[i];
        residualNorm += residual[i] * residual[i];
        rightHandNorm += rhs[i] * rhs[i];
    }

    const double target = tolerance * tolerance * std::max(rightHandNorm, 1e-300);
    if (residualNorm <= target)
    {
        return 0;
    }

    int iteration = 0;
    for (; iteration < maxIterations; iteration++)
    {
        applyK(direction, scratch);
        double curvature = 0.0;
        for (int i = 0; i < n; i++)
        {
            curvature += direction[i] * scratch[i];
        }
        if (!(curvature > 0.0))
        {
            // K is positive definite on a non-bipartite graph, and a
            // triangulation has triangles. Reaching here means the structure
            // is not what this assumes, and quietly returning a half-converged
            // answer would put a wrong mobility into the dynamics.
            throw std::runtime_error(
                "[SurfaceSolver] the limit-mask system is not positive definite (direction "
                "curvature " +
                std::to_string(curvature) +
                "). The mesh is probably not a two-manifold triangulation.");
        }
        const double step = residualNorm / curvature;

        double nextResidualNorm = 0.0;
        for (int i = 0; i < n; i++)
        {
            solution[i] += step * direction[i];
            residual[i] -= step * scratch[i];
            nextResidualNorm += residual[i] * residual[i];
        }
        if (nextResidualNorm <= target)
        {
            iteration++;
            break;
        }
        const double beta = nextResidualNorm / residualNorm;
        for (int i = 0; i < n; i++)
        {
            direction[i] = residual[i] + beta * direction[i];
        }
        residualNorm = nextResidualNorm;
    }

    if (iteration >= maxIterations)
    {
        throw std::runtime_error("[SurfaceSolver] the limit-mask solve did not converge in " +
                                 std::to_string(maxIterations) +
                                 " iterations. Raise SurfaceSolver::maxIterations, or check the "
                                 "mesh for a degenerate neighbourhood.");
    }
    return iteration;
}

int SurfaceSolver::surface_to_mesh(const Matrix &surface, Matrix &control,
                                   const Matrix *initialGuess) const
{
    const int n = static_cast<int>(freeOfCompact_.size());

    // A pinned vertex carries the identity, so its control point is its
    // surface point and it is known before anything is solved.
    for (int v = 0; v < nVertices_; v++)
    {
        if (pinned_[v])
        {
            for (int axis = 0; axis < 3; axis++)
            {
                control.set(v, axis, surface(v, axis));
            }
        }
    }

    int totalIterations = 0;
    std::vector<double> rhs(n);
    std::vector<double> solution;
    std::vector<double> guess;
    lastResidual_ = 0.0;

    for (int axis = 0; axis < 3; axis++)
    {
        // Multiply the free rows through by their valence: the mask row
        //     1/2 C_v + sum over neighbours of C_u / (2 N_v) = S_v
        // becomes
        //     N_v C_v + sum over neighbours of C_u = 2 N_v S_v,
        // which is 2 K C = 2 D S with K symmetric. The known pinned
        // neighbours move across to the right.
        for (int i = 0; i < n; i++)
        {
            const int v = freeOfCompact_[i];
            double value = valence_[v] * surface(v, axis);
            for (int k = pinnedRowStart_[i]; k < pinnedRowStart_[i + 1]; k++)
            {
                value -= 0.5 * control(pinnedColumn_[k], axis);
            }
            rhs[i] = value;
        }

        if (initialGuess != nullptr)
        {
            guess.resize(n);
            for (int i = 0; i < n; i++)
            {
                guess[i] = (*initialGuess)(freeOfCompact_[i], axis);
            }
        }
        totalIterations += solve_free(rhs, solution, (initialGuess != nullptr) ? &guess : nullptr);

        for (int i = 0; i < n; i++)
        {
            control.set(freeOfCompact_[i], axis, solution[i]);
        }
    }

    // What the caller actually cares about: how well M C reproduces S.
    for (int v = 0; v < nVertices_; v++)
    {
        if (pinned_[v])
        {
            continue;
        }
        const int compact = compactOfVertex_[v];
        const double neighbourWeight = 0.5 / valence_[v];
        for (int axis = 0; axis < 3; axis++)
        {
            double sum = 0.5 * control(v, axis);
            for (int k = rowStart_[compact]; k < rowStart_[compact + 1]; k++)
            {
                sum += neighbourWeight * control(freeOfCompact_[column_[k]], axis);
            }
            for (int k = pinnedRowStart_[compact]; k < pinnedRowStart_[compact + 1]; k++)
            {
                sum += neighbourWeight * control(pinnedColumn_[k], axis);
            }
            lastResidual_ = std::max(lastResidual_, std::abs(sum - surface(v, axis)));
        }
    }
    return totalIterations;
}

int SurfaceSolver::nodal_force_to_surface(const Matrix &nodalForce, Matrix &surfaceForce) const
{
    const int n = static_cast<int>(freeOfCompact_.size());

    // Mᵀ = (D⁻¹ K)ᵀ = K D⁻¹ over the free vertices, so solving Mᵀ F_S = F_C
    // is solving K y = F_C and then scaling: F_S = D y.
    //
    // Pinned vertices are left at zero. They are not integrated -- the
    // Brownian step skips ghosts and the periodic post-process overwrites the
    // duplicates -- so a force on them would be discarded anyway, and writing
    // one would suggest otherwise.
    for (int v = 0; v < nVertices_; v++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            surfaceForce.set(v, axis, 0.0);
        }
    }

    int totalIterations = 0;
    std::vector<double> rhs(n);
    std::vector<double> solution;
    for (int axis = 0; axis < 3; axis++)
    {
        for (int i = 0; i < n; i++)
        {
            rhs[i] = nodalForce(freeOfCompact_[i], axis);
        }
        totalIterations += solve_free(rhs, solution, nullptr);
        for (int i = 0; i < n; i++)
        {
            const int v = freeOfCompact_[i];
            surfaceForce.set(v, axis, valence_[v] * solution[i]);
        }
    }
    return totalIterations;
}

} // namespace slimed
