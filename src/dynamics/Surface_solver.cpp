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
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "mesh/Mesh.hpp"

namespace slimed
{

namespace
{
/// Pack an ordered pair of compact indices into one key.
inline std::uint64_t ordered_pair_key(int from, int to)
{
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(from)) << 32) |
           static_cast<std::uint32_t>(to);
}
} // namespace

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
    imageRoot_.assign(nVertices_, -1);
    imageOffset_.assign(static_cast<std::size_t>(nVertices_) * 3, 0.0);
    imageVertices_.clear();
    neighbourOffsetSum_.clear();
    topologyVersion = mesh.topologyVersion;

    perVertexBoundary_ = (mesh.param.boundaryCondition == BoundaryType::Mixed);

    // Under the per-vertex boundary an image is an affine function of its
    // source and gets no row: it is written from the source after every
    // solve, and wherever a free row reaches it the source is substituted.
    // The mesh has already flattened every mirror chain, so a root is never
    // itself an image.
    if (perVertexBoundary_)
    {
        for (int v : mesh.periodicImageVertices)
        {
            const Vertex &image = mesh.vertices[v];
            imageRoot_[v] = image.reflectiveVertexIndex;
            for (int axis = 0; axis < 3; axis++)
            {
                imageOffset_[static_cast<std::size_t>(v) * 3 + axis] = image.mirrorOffset[axis];
            }
            imageVertices_.push_back(v);
        }
    }

    for (int v = 0; v < nVertices_; v++)
    {
        if (imageRoot_[v] >= 0)
        {
            // Neither pinned nor free: written from its source.
            valence_[v] = 0;
            continue;
        }

        const Vertex &vertex = mesh.vertices[v];
        const int valence = static_cast<int>(vertex.adjacentVertices.size());
        bool pinned = false;
        if (perVertexBoundary_)
        {
            // A clamped or ghost vertex is known before anything is solved.
            // Everything else with a fan is a degree of freedom -- whether or
            // not the faces around it in this copy of the sheet carry energy,
            // because its physical fan is complete somewhere among its images.
            pinned = vertex.is_fixed() || vertex.type == VertexType::Ghost || vertex.isGhost ||
                     valence == 0;
        }
        else
        {
            // A vertex with no adjacent face that contributes to the physical
            // surface has no limit point of its own, and the mask gives it the
            // identity. That is the same test the dense path applies, kept
            // identical so the two agree about which vertices are degrees of
            // freedom.
            bool hasRealFace = false;
            for (int iFace : vertex.adjacentFaces)
            {
                const Face &face = mesh.faces[iFace];
                if (!face.isGhost && !face.isBoundary)
                {
                    hasRealFace = true;
                    break;
                }
            }
            pinned = (!hasRealFace || valence == 0);
        }
        pinned_[v] = pinned ? 1 : 0;
        valence_[v] = pinned ? 0 : valence;
    }

    for (int v = 0; v < nVertices_; v++)
    {
        if (!pinned_[v] && imageRoot_[v] < 0)
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
    neighbourOffsetSum_.assign(static_cast<std::size_t>(nFreeVertices) * 3, 0.0);
    for (int compact = 0; compact < nFreeVertices; compact++)
    {
        const int v = freeOfCompact_[compact];
        for (int neighbour : mesh.vertices[v].adjacentVertices)
        {
            // An image neighbour stands for its source, and carries its offset
            // into the constant of this row.
            const int root = (imageRoot_[neighbour] >= 0) ? imageRoot_[neighbour] : neighbour;
            if (imageRoot_[neighbour] >= 0)
            {
                for (int axis = 0; axis < 3; axis++)
                {
                    neighbourOffsetSum_[static_cast<std::size_t>(compact) * 3 + axis] +=
                        imageOffset_[static_cast<std::size_t>(neighbour) * 3 + axis];
                }
            }
            if (pinned_[root])
            {
                pinnedColumn_.push_back(root);
            }
            else
            {
                column_.push_back(compactOfVertex_[root]);
            }
        }
        rowStart_.push_back(static_cast<int>(column_.size()));
        pinnedRowStart_.push_back(static_cast<int>(pinnedColumn_.size()));
    }

    // The conjugate gradient below assumes K symmetric, which on the wrapped
    // adjacency holds exactly when the image structure is consistent: if an
    // image of w is adjacent to v, then an image of v is adjacent to w. A
    // mesh whose images say otherwise would be solved silently wrong, so it
    // is refused here, naming the pair.
    if (perVertexBoundary_)
    {
        std::unordered_map<std::uint64_t, int> count;
        count.reserve(column_.size());
        for (int i = 0; i < nFreeVertices; i++)
        {
            for (int k = rowStart_[i]; k < rowStart_[i + 1]; k++)
            {
                count[ordered_pair_key(i, column_[k])]++;
            }
        }
        for (const auto &entry : count)
        {
            const int i = static_cast<int>(entry.first >> 32);
            const int j = static_cast<int>(entry.first & 0xFFFFFFFFu);
            const auto back = count.find(ordered_pair_key(j, i));
            const int reverse = (back == count.end()) ? 0 : back->second;
            if (reverse != entry.second)
            {
                throw std::runtime_error(
                    "[SurfaceSolver] the periodic image structure is not consistent: vertex " +
                    std::to_string(freeOfCompact_[i]) + " is adjacent to vertex " +
                    std::to_string(freeOfCompact_[j]) + " or its images " +
                    std::to_string(entry.second) + " time(s), but vertex " +
                    std::to_string(freeOfCompact_[j]) + " is adjacent to vertex " +
                    std::to_string(freeOfCompact_[i]) + " or its images " +
                    std::to_string(reverse) +
                    " time(s). Every image relation must be mirrored: if an image of one "
                    "vertex neighbours another, an image of the other must neighbour the one.");
            }
        }
    }
}

void SurfaceSolver::mesh_to_surface(const Matrix &control, Matrix &surface) const
{
    for (int v = 0; v < nVertices_; v++)
    {
        if (imageRoot_[v] >= 0)
        {
            continue; // written from its source below
        }
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
            if (perVertexBoundary_)
            {
                sum += neighbourWeight *
                       neighbourOffsetSum_[static_cast<std::size_t>(compact) * 3 + axis];
            }
            surface.set(v, axis, sum);
        }
    }
    for (int v : imageVertices_)
    {
        const int root = imageRoot_[v];
        for (int axis = 0; axis < 3; axis++)
        {
            surface.set(v, axis,
                        surface(root, axis) + imageOffset_[static_cast<std::size_t>(v) * 3 + axis]);
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
        // neighbours move across to the right, and so does the offset an
        // image neighbour carries.
        for (int i = 0; i < n; i++)
        {
            const int v = freeOfCompact_[i];
            double value = valence_[v] * surface(v, axis);
            for (int k = pinnedRowStart_[i]; k < pinnedRowStart_[i + 1]; k++)
            {
                value -= 0.5 * control(pinnedColumn_[k], axis);
            }
            if (perVertexBoundary_)
            {
                value -= 0.5 * neighbourOffsetSum_[static_cast<std::size_t>(i) * 3 + axis];
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

    // An image is its source plus its offset, exactly.
    for (int v : imageVertices_)
    {
        const int root = imageRoot_[v];
        for (int axis = 0; axis < 3; axis++)
        {
            control.set(v, axis,
                        control(root, axis) + imageOffset_[static_cast<std::size_t>(v) * 3 + axis]);
        }
    }

    // What the caller actually cares about: how well M C reproduces S.
    for (int v = 0; v < nVertices_; v++)
    {
        if (pinned_[v] || imageRoot_[v] >= 0)
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
            if (perVertexBoundary_)
            {
                sum += neighbourWeight *
                       neighbourOffsetSum_[static_cast<std::size_t>(compact) * 3 + axis];
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
    // Pinned vertices and images are left at zero. They are not integrated --
    // the Brownian step skips ghosts, clamped vertices and images, and the
    // periodic post-process overwrites the duplicates -- so a force on them
    // would be discarded anyway, and writing one would suggest otherwise. An
    // image's force has already been folded onto its source by the time the
    // nodal force reaches here.
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
