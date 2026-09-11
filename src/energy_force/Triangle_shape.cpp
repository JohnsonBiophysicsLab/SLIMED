/**
 * @file Triangle_shape.cpp
 * @brief The triangle-shape term: the fluid mesh-quality term's second half.
 *
 * The flat tether bounds every edge of a fluid mesh and nothing else, and
 * that is not enough. Under in-plane motion and flips the mesh fills with
 * slivers whose edges all sit inside the walls -- measured on the 100 nm
 * sheet, interior triangles 0.1 nm tall at one degree in every frame. A
 * sliver's normal is ill-conditioned: it turns through tens of degrees under
 * a single 0.05 nm Brownian kick, the limit-surface patch built on it
 * self-intersects, and the bending force on it is not finite. Every long
 * fluid run ended that way. A Monte Carlo model never takes such a step,
 * because the energy rejects it; an explicit Brownian step has no refusal in
 * it, so the Hamiltonian has to carry one.
 *
 * The quantity that vanishes in a sliver is an altitude -- the distance from
 * a corner to its opposite edge, `h_i = 2A / l_i` -- so that is what the term
 * bounds, with the same flat-bottomed wall the tether uses:
 *
 *     E_face = (k / 2) sum over the three corners  max(0, h0 - h_i)^2
 *
 * Zero for any healthy triangle, quadratic below the floor. It adds no
 * tension and does nothing inside the allowed region.
 *
 * The gradient is closed-form. With `n = (p1 - p0) x (p2 - p0)`, `A = |n|/2`
 * and `e_k = p_{k+1} - p_{k+2}` the edge opposite corner k,
 *
 *     dA / dp_k   = (1/2) e_k x n_hat
 *     dl_i / dp_k = +e_i / l_i  for k = i+1,  -e_i / l_i  for k = i+2,  0 for k = i
 *     dh_i / dp_k = (2 / l_i) dA/dp_k  -  (2A / l_i^2) dl_i/dp_k
 *
 * and the force on corner k is `k_s (h0 - h_i)_+ dh_i/dp_k` summed over i.
 * TriangleShapeTest pins it against finite differences.
 *
 * One expression serves the force pass and the flip trial, through
 * triangle_shape() below: the WP6 lesson was that two implementations of a
 * mesh-quality term drift apart, and a Metropolis chain differencing one
 * while the dynamics integrates the other samples neither.
 *
 * @see docs/edge_flip_plan.md work package 7
 */

#include "mesh/Mesh.hpp"

#include <cmath>

#include "energy_force/Patch_kernel.hpp"

namespace
{
/**
 * @brief Energy of one triangle's altitude walls, and optionally the force
 * on its three corners.
 *
 * @param p      The three corners, `p[k][axis]`.
 * @param floor  The altitude `h0` below which a corner is charged.
 * @param k      The wall stiffness.
 * @param force  If not null, filled with minus the gradient, `force[k][axis]`.
 */
double triangle_shape(const double p[3][3], double floor, double k, double force[3][3])
{
    if (force != nullptr)
    {
        for (int corner = 0; corner < 3; corner++)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                force[corner][axis] = 0.0;
            }
        }
    }

    // e_i: the edge opposite corner i, from p_{i+2} to p_{i+1}.
    double e[3][3];
    double l[3];
    for (int i = 0; i < 3; i++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            e[i][axis] = p[(i + 1) % 3][axis] - p[(i + 2) % 3][axis];
        }
        l[i] = slimed::v3_norm(e[i]);
    }

    double d1[3];
    double d2[3];
    for (int axis = 0; axis < 3; axis++)
    {
        d1[axis] = p[1][axis] - p[0][axis];
        d2[axis] = p[2][axis] - p[0][axis];
    }
    double n[3];
    slimed::v3_cross(d1, d2, n);
    const double twoArea = slimed::v3_norm(n);
    const double area = 0.5 * twoArea;
    double nHat[3] = {0.0, 0.0, 0.0};
    if (twoArea > 0.0)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            nHat[axis] = n[axis] / twoArea;
        }
    }

    double energy = 0.0;
    for (int i = 0; i < 3; i++)
    {
        if (l[i] <= 0.0)
        {
            continue; // coincident corners: no edge to measure against
        }
        const double h = 2.0 * area / l[i];
        if (h >= floor)
        {
            continue;
        }
        const double shortfall = floor - h;
        energy += 0.5 * k * shortfall * shortfall;

        if (force == nullptr || twoArea <= 0.0)
        {
            // A collinear triangle has no normal to push along. The energy is
            // charged; the direction is left to the next step, by which time
            // the noise has broken the degeneracy.
            continue;
        }
        for (int corner = 0; corner < 3; corner++)
        {
            double dArea[3];
            slimed::v3_cross(e[corner], nHat, dArea);
            double dLength[3] = {0.0, 0.0, 0.0};
            const double sign = (corner == (i + 1) % 3) ? 1.0 : (corner == (i + 2) % 3) ? -1.0 : 0.0;
            for (int axis = 0; axis < 3; axis++)
            {
                dArea[axis] *= 0.5;
                dLength[axis] = sign * e[i][axis] / l[i];
                const double dh = (2.0 / l[i]) * dArea[axis] - (2.0 * area / (l[i] * l[i])) * dLength[axis];
                force[corner][axis] += k * shortfall * dh; // -dE/dp = k (h0 - h) dh/dp
            }
        }
    }
    return energy;
}
} // namespace

double Mesh::face_shape_energy(int iFace) const
{
    const Face &face = faces[iFace];
    if (face.isGhost || face.adjacentVertices.size() != 3)
    {
        return 0.0;
    }
    double p[3][3];
    for (int corner = 0; corner < 3; corner++)
    {
        const Matrix &coord = vertices[face.adjacentVertices[corner]].coord;
        for (int axis = 0; axis < 3; axis++)
        {
            p[corner][axis] = coord.get(axis, 0);
        }
    }
    return triangle_shape(p, param.triangleShapeMinAltitudeRatio * param.lFace,
                          param.triangleShapeConstant, nullptr);
}

void Mesh::energy_force_triangle_shape()
{
    const double floor = param.triangleShapeMinAltitudeRatio * param.lFace;
    const double k = param.triangleShapeConstant;

    const int nVertices = static_cast<int>(vertices.size());
    const int nFaces = static_cast<int>(faces.size());
#ifdef OMP
    const int nThreads = omp_get_max_threads();
#else
    const int nThreads = 1;
#endif
    // Per-thread accumulation: neighbouring faces share corners.
    std::vector<std::vector<double>> forceComponents(
        nThreads, std::vector<double>(static_cast<std::size_t>(nVertices) * 3, 0.0));

#pragma omp parallel for
    for (int iFace = 0; iFace < nFaces; iFace++)
    {
        Face &face = faces[iFace];
        if (face.isGhost || face.adjacentVertices.size() != 3)
        {
            continue;
        }
#ifdef OMP
        const int threadIndex = omp_get_thread_num();
#else
        const int threadIndex = 0;
#endif
        std::vector<double> &local = forceComponents[threadIndex];

        double p[3][3];
        for (int corner = 0; corner < 3; corner++)
        {
            const Matrix &coord = vertices[face.adjacentVertices[corner]].coord;
            for (int axis = 0; axis < 3; axis++)
            {
                p[corner][axis] = coord.get(axis, 0);
            }
        }
        double force[3][3];
        const double energy = triangle_shape(p, floor, k, force);
        if (energy == 0.0)
        {
            continue;
        }
        // Added to whatever the edge-based term already put there: the two
        // are halves of one mesh-quality term and share the slot.
        face.energy.energyRegularization += energy;
        for (int corner = 0; corner < 3; corner++)
        {
            const std::size_t base = static_cast<std::size_t>(face.adjacentVertices[corner]) * 3;
            for (int axis = 0; axis < 3; axis++)
            {
                local[base + axis] += force[corner][axis];
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
            vertices[i].force.forceRegularization.set(
                axis, 0, vertices[i].force.forceRegularization.get(axis, 0) + sums[axis]);
        }
    }
}
