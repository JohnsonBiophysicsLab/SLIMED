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
 * @see docs/edge_flip_plan.md section 1.6
 */

#include "mesh/Mesh.hpp"

#include <cmath>

#include "energy_force/Patch_kernel.hpp"

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

        const double extension = length - restLength;
        const double energy = 0.5 * springConstant * extension * extension;

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
