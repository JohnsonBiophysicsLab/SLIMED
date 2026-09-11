/**
 * @file Edge_spring.cpp
 * @brief The fluid-mode mesh-quality term.
 *
 * Every dynamically triangulated surface model needs something to keep its
 * triangles from degenerating, and every one of them uses a term on the edge
 * lengths: hard walls in the Monte Carlo models, a smooth well in the ones
 * that integrate. The Metropolis flip then prefers the shorter diagonal on its
 * own, because the longer one costs energy, so mesh quality comes out of the
 * physics rather than out of a separate repair pass.
 *
 * SLIMED's existing regularization does something different. It measures each
 * face's edges against *the same face's edges in coordRef* -- a solid's memory
 * of the configuration it started in. For a fluid that is the wrong object
 * twice over: the membrane has no reference configuration to remember, and an
 * edge that a flip has just created never had a reference length. coordRef
 * would hand it the distance between two vertices that were not joined, which
 * is not a rest length but an accident of where they happened to be.
 *
 *     E = (k / 2) * sum over edges of (l_e - l0)^2
 *
 * depends only on the edges that exist now, so it survives a flip and means
 * the same thing before and after one. Its physical footprint is a small
 * contribution to the area modulus and, with flips running, none to the shear
 * modulus -- which is the point. The fluctuation spectrum's fitted tension is
 * what says whether the spring is quietly adding one.
 *
 * ### Why that spring does not work, and what replaced it
 *
 * A flip on a rhombus of two equilateral triangles of side `l` replaces the
 * short diagonal by the long one, so it has to climb
 * `(k/2)(sqrt(3) - 1)^2 l^2` no matter what the rest of the Hamiltonian says.
 * Two requirements then pull `k` in opposite directions -- flips need a
 * barrier of a few kT, and a triangulation that does not degenerate needs the
 * bond fluctuation `sqrt(kT/k)` well below `l` -- and at `l = 5 nm` and room
 * temperature there is no `k` that satisfies both. Measured on a 100 nm sheet
 * over 3000 steps: `k = 83.4` accepts 0 flips of 121, `k = 20` accepts 2.1%,
 * `k = 1` accepts 8.9% and the run diverges at step 1535.
 *
 * So the default is the flat-bottomed tether the Monte Carlo models use:
 *
 *     E = (k/2)(l_min - l)^2   below the range
 *         0                    inside it
 *         (k/2)(l - l_max)^2   above it
 *
 * A flip that leaves every edge inside the range costs nothing, so the wall
 * stiffness and the flip barrier stop being the same number. `k` can then be
 * whatever the walls need. The harmonic form is kept behind
 * `edgeTetherShape = harmonic`, because it is the smooth one and a
 * minimization that never flips has no reason to prefer a piecewise term.
 *
 * @see docs/edge_flip_plan.md sections 1.6 and work packages 5 and 6
 */

#include "mesh/Mesh.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

#include "energy_force/Patch_kernel.hpp"

double Mesh::edge_tether_energy(double length) const
{
    const double restLength =
        (param.edgeSpringRestLength >= 0.0) ? param.edgeSpringRestLength : param.lFace;

    // Displacement from the nearest point of the allowed interval: zero inside
    // it, signed outside. The harmonic form is the degenerate case where the
    // interval is the single point l0, which is what lets one expression serve
    // both shapes.
    double extension = 0.0;
    if (param.edgeTetherShape == "harmonic")
    {
        extension = length - restLength;
    }
    else
    {
        const double lowerBound = param.edgeTetherMinRatio * restLength;
        const double upperBound = param.edgeTetherMaxRatio * restLength;
        if (!(lowerBound <= upperBound))
        {
            throw std::runtime_error(
                "[Mesh::edge_tether_energy] edgeTetherMinRatio exceeds edgeTetherMaxRatio; the "
                "tether has no allowed range, so every edge is against a wall.");
        }
        if (length < lowerBound)
        {
            extension = length - lowerBound;
        }
        else if (length > upperBound)
        {
            extension = length - upperBound;
        }
    }
    return 0.5 * param.edgeSpringConstant * extension * extension;
}

double Mesh::face_tether_energy(int iFace) const
{
    const std::vector<int> &corners = faces[iFace].adjacentVertices;
    if (corners.size() != 3)
    {
        return 0.0;
    }

    double sum = 0.0;
    for (int k = 0; k < 3; k++)
    {
        const int a = corners[k];
        const int b = corners[(k + 1) % 3];
        double along[3];
        for (int axis = 0; axis < 3; axis++)
        {
            along[axis] = vertices[a].coord.get(axis, 0) - vertices[b].coord.get(axis, 0);
        }
        const double energy = edge_tether_energy(slimed::v3_norm(along));

        // The same split energy_force_edge_spring() applies: half to each of
        // the two faces an edge separates, all of it to the one face of a
        // boundary edge. Read off the edge table so the two agree even where
        // the table and the corner walk would disagree about incidence.
        const int iEdge = edge_between(a, b);
        int nIncident = 2;
        if (iEdge >= 0)
        {
            const MeshEdge &edge = edges[iEdge];
            nIncident = (edge.face[0] >= 0 ? 1 : 0) + (edge.face[1] >= 0 ? 1 : 0);
        }
        if (nIncident > 0)
        {
            sum += energy / nIncident;
        }
    }
    return sum;
}

void Mesh::energy_force_edge_spring()
{
    if (edges.empty())
    {
        return;
    }

    const double springConstant = param.edgeSpringConstant;
    // A negative rest length means "the edge length the mesh was built for".
    const double restLength =
        (param.edgeSpringRestLength >= 0.0) ? param.edgeSpringRestLength : param.lFace;

    const bool flatBottomed = (param.edgeTetherShape != "harmonic");
    const double lowerBound = param.edgeTetherMinRatio * restLength;
    const double upperBound = param.edgeTetherMaxRatio * restLength;
    if (flatBottomed && !(lowerBound <= upperBound))
    {
        throw std::runtime_error(
            "[Mesh::energy_force_edge_spring] edgeTetherMinRatio exceeds edgeTetherMaxRatio; the "
            "tether has no allowed range, so every edge is against a wall.");
    }

    // Displacement from the nearest point of the allowed interval: zero inside
    // it, signed outside. The harmonic form is the degenerate case where the
    // interval is the single point l0, which is what makes one expression
    // serve both shapes -- and what makes the force below identical in form.
    const auto displacementFromRest = [&](double length) {
        if (!flatBottomed)
        {
            return length - restLength;
        }
        if (length < lowerBound)
        {
            return length - lowerBound;
        }
        if (length > upperBound)
        {
            return length - upperBound;
        }
        return 0.0;
    };

    for (Face &face : faces)
    {
        face.energy.energyRegularization = 0.0;
    }

    const int nVertices = static_cast<int>(vertices.size());
#ifdef OMP
    const int nThreads = omp_get_max_threads();
#else
    const int nThreads = 1;
#endif
    // Per-thread accumulation: an edge writes to both its endpoints, and
    // neighbouring edges share them.
    std::vector<std::vector<double>> forceComponents(
        nThreads, std::vector<double>(static_cast<std::size_t>(nVertices) * 3, 0.0));

    const int nEdges = static_cast<int>(edges.size());
#pragma omp parallel for
    for (int iEdge = 0; iEdge < nEdges; iEdge++)
    {
        const MeshEdge &edge = edges[iEdge];
#ifdef OMP
        const int threadIndex = omp_get_thread_num();
#else
        const int threadIndex = 0;
#endif
        std::vector<double> &local = forceComponents[threadIndex];

        const int a = edge.v[0];
        const int b = edge.v[1];
        double along[3];
        for (int axis = 0; axis < 3; axis++)
        {
            along[axis] = vertices[a].coord.get(axis, 0) - vertices[b].coord.get(axis, 0);
        }
        const double length = slimed::v3_norm(along);
        if (length <= 0.0)
        {
            continue; // coincident vertices carry no direction to push along
        }

        const double extension = displacementFromRest(length);
        const double energy = 0.5 * springConstant * extension * extension;
        // edge_tether_energy() is the same expression, reached by the flip
        // trial. Kept as one definition would be better still, but the force
        // below needs the signed displacement and not just the energy.

        // -dE/dx_a, with the sign flipped on b: a stretched edge pulls its
        // endpoints together.
        const double scale = springConstant * extension / length;
        for (int axis = 0; axis < 3; axis++)
        {
            const double component = scale * along[axis];
            local[static_cast<std::size_t>(a) * 3 + axis] -= component;
            local[static_cast<std::size_t>(b) * 3 + axis] += component;
        }

        // The total energy is summed over faces, so the edge's share has to
        // land on them. Split evenly between the two it separates, or given
        // whole to the one face of a boundary edge.
        const int nIncident = (edge.face[0] >= 0 ? 1 : 0) + (edge.face[1] >= 0 ? 1 : 0);
        if (nIncident == 0)
        {
            continue;
        }
        const double share = energy / nIncident;
        for (int k = 0; k < 2; k++)
        {
            if (edge.face[k] >= 0)
            {
#pragma omp atomic
                faces[edge.face[k]].energy.energyRegularization += share;
            }
        }
    }

#pragma omp parallel for
    for (int i = 0; i < nVertices; i++)
    {
        double sums[3] = {0.0, 0.0, 0.0};
        for (int threadIndex = 0; threadIndex < nThreads; ++threadIndex)
        {
            const std::vector<double> &local = forceComponents[threadIndex];
            for (int axis = 0; axis < 3; axis++)
            {
                sums[axis] += local[static_cast<std::size_t>(i) * 3 + axis];
            }
        }
        for (int axis = 0; axis < 3; axis++)
        {
            vertices[i].force.forceRegularization.set(axis, 0, sums[axis]);
        }
    }

    // The deformation counters describe the reference-length term and mean
    // nothing here; zero them rather than leave the last pass's numbers to be
    // read as this one's.
    param.deformationCount.shapeDeformCount = 0;
    param.deformationCount.areaDeformCount = 0;
    param.deformationCount.noDeformCount = static_cast<int>(faces.size());
}
