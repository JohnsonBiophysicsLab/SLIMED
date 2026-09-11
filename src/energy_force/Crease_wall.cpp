/**
 * @file Crease_wall.cpp
 * @brief The crease wall: the fluid mesh-quality term's third half.
 *
 * A dynamically triangulated surface never folds a face onto its neighbour,
 * because its bending energy is a sum over edges of `1 - n1 . n2` on the
 * control net and a fold costs the bending modulus outright. SLIMED's
 * bending energy lives on the limit surface, which is smoother than the
 * control net: a crease in the net is smoothed away until the net is two
 * layers deep, at which point the limit surface pinches, its curvature
 * integrand diverges, and the explicit step blows up. Measured on the 100 nm
 * sheet (WP7): with the altitude floor holding every triangle healthy, the
 * control net carried creases of 177-180 degrees from step 20 000 and folded
 * at 82 000.
 *
 * This term puts the tether's wall on the crease. With `c = n1 . n2` the
 * cosine of the angle between an edge's two face normals,
 *
 *     E_edge = (k / 2) max(0, c0 - c)^2,    c0 = cos(theta_max)
 *
 * Zero while the faces are within theta_max of coplanar; quadratic in the
 * cosine beyond it. In the cosine rather than the angle because the angle's
 * gradient is singular at a full fold, which is exactly where the force has
 * to be well defined.
 *
 * The gradient is closed-form. For a corner p_k of face 1, with e_k the edge
 * opposite it (p_{k+1} - p_{k+2}) and n1 the face's unnormalised normal,
 * `dn1 = dp x e_k`, and with `m1 = (n2 - c n1) / |n1|` (unit normals here)
 *
 *     dc / dp_k = e_k x m1
 *
 * and the same with the faces swapped for a corner of face 2. The two shared
 * corners collect both. The force is `k (c0 - c) dc/dp`. CreaseWallTest pins
 * it against finite differences.
 *
 * Attributed half to each face, so that the sum over faces is the sum over
 * edges and the flip trial differences what the force pass computes. Stale
 * copies -- edges with both endpoints ghost or duplicate -- carry nothing,
 * as for the tether.
 *
 * @see docs/edge_flip_plan.md work package 7
 */

#include "mesh/Mesh.hpp"

#include <cmath>

#include "energy_force/Patch_kernel.hpp"

namespace
{
struct CreaseGeometry
{
    double cosine = 1.0;
    double energy = 0.0;
    bool charged = false;
    /// Minus the gradient with respect to each corner of face 1 and face 2.
    double force1[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double force2[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
};

void unit_normal(const double p[3][3], double nHat[3], double &length)
{
    double d1[3];
    double d2[3];
    for (int axis = 0; axis < 3; axis++)
    {
        d1[axis] = p[1][axis] - p[0][axis];
        d2[axis] = p[2][axis] - p[0][axis];
    }
    double n[3];
    slimed::v3_cross(d1, d2, n);
    length = slimed::v3_norm(n);
    for (int axis = 0; axis < 3; axis++)
    {
        nHat[axis] = (length > 0.0) ? n[axis] / length : 0.0;
    }
}

/// Energy of one edge's crease wall and, if asked, the forces on both faces' corners.
CreaseGeometry crease(const double p1[3][3], const double p2[3][3], double cosineFloor, double k,
                      bool withForce)
{
    CreaseGeometry out;
    double n1[3];
    double n2[3];
    double length1 = 0.0;
    double length2 = 0.0;
    unit_normal(p1, n1, length1);
    unit_normal(p2, n2, length2);
    if (length1 <= 0.0 || length2 <= 0.0)
    {
        return out; // a degenerate face has no normal; the shape term owns that case
    }
    out.cosine = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
    if (out.cosine >= cosineFloor)
    {
        return out;
    }
    const double shortfall = cosineFloor - out.cosine;
    out.energy = 0.5 * k * shortfall * shortfall;
    out.charged = true;
    if (!withForce)
    {
        return out;
    }

    // m1 = (n2 - c n1) / |n1|, m2 = (n1 - c n2) / |n2|.
    double m1[3];
    double m2[3];
    for (int axis = 0; axis < 3; axis++)
    {
        m1[axis] = (n2[axis] - out.cosine * n1[axis]) / length1;
        m2[axis] = (n1[axis] - out.cosine * n2[axis]) / length2;
    }
    for (int corner = 0; corner < 3; corner++)
    {
        double e1[3];
        double e2[3];
        for (int axis = 0; axis < 3; axis++)
        {
            e1[axis] = p1[(corner + 1) % 3][axis] - p1[(corner + 2) % 3][axis];
            e2[axis] = p2[(corner + 1) % 3][axis] - p2[(corner + 2) % 3][axis];
        }
        double g1[3];
        double g2[3];
        slimed::v3_cross(e1, m1, g1); // dc/dp for a corner of face 1
        slimed::v3_cross(e2, m2, g2); // dc/dp for a corner of face 2
        for (int axis = 0; axis < 3; axis++)
        {
            // -dE/dp = k (c0 - c) dc/dp
            out.force1[corner][axis] = k * shortfall * g1[axis];
            out.force2[corner][axis] = k * shortfall * g2[axis];
        }
    }
    return out;
}
} // namespace

double Mesh::face_crease_energy(int iFace) const
{
    const Face &face = faces[iFace];
    if (face.adjacentVertices.size() != 3)
    {
        return 0.0;
    }
    const double cosineFloor = std::cos(param.creaseWallAngle * M_PI / 180.0);
    double sum = 0.0;
    for (int k = 0; k < 3; k++)
    {
        const int iEdge = edge_between(face.adjacentVertices[k], face.adjacentVertices[(k + 1) % 3]);
        if (iEdge < 0)
        {
            continue;
        }
        const MeshEdge &edge = edges[iEdge];
        if (edge.face[0] < 0 || edge.face[1] < 0 || !edge_carries_tether(edge))
        {
            continue;
        }
        double p1[3][3];
        double p2[3][3];
        for (int corner = 0; corner < 3; corner++)
        {
            const Matrix &c1 = vertices[faces[edge.face[0]].adjacentVertices[corner]].coord;
            const Matrix &c2 = vertices[faces[edge.face[1]].adjacentVertices[corner]].coord;
            for (int axis = 0; axis < 3; axis++)
            {
                p1[corner][axis] = c1.get(axis, 0);
                p2[corner][axis] = c2.get(axis, 0);
            }
        }
        sum += 0.5 * crease(p1, p2, cosineFloor, param.creaseWallConstant, false).energy;
    }
    return sum;
}

void Mesh::energy_force_crease_wall()
{
    const double cosineFloor = std::cos(param.creaseWallAngle * M_PI / 180.0);
    const double k = param.creaseWallConstant;

    const int nVertices = static_cast<int>(vertices.size());
    const int nEdges = static_cast<int>(edges.size());
#ifdef OMP
    const int nThreads = omp_get_max_threads();
#else
    const int nThreads = 1;
#endif
    std::vector<std::vector<double>> forceComponents(
        nThreads, std::vector<double>(static_cast<std::size_t>(nVertices) * 3, 0.0));

#pragma omp parallel for
    for (int iEdge = 0; iEdge < nEdges; iEdge++)
    {
        const MeshEdge &edge = edges[iEdge];
        if (edge.face[0] < 0 || edge.face[1] < 0 || !edge_carries_tether(edge))
        {
            continue;
        }
#ifdef OMP
        const int threadIndex = omp_get_thread_num();
#else
        const int threadIndex = 0;
#endif
        std::vector<double> &local = forceComponents[threadIndex];

        const std::vector<int> &corners1 = faces[edge.face[0]].adjacentVertices;
        const std::vector<int> &corners2 = faces[edge.face[1]].adjacentVertices;
        double p1[3][3];
        double p2[3][3];
        for (int corner = 0; corner < 3; corner++)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                p1[corner][axis] = vertices[corners1[corner]].coord.get(axis, 0);
                p2[corner][axis] = vertices[corners2[corner]].coord.get(axis, 0);
            }
        }
        const CreaseGeometry g = crease(p1, p2, cosineFloor, k, true);
        if (!g.charged)
        {
            continue;
        }
        for (int corner = 0; corner < 3; corner++)
        {
            const std::size_t base1 = static_cast<std::size_t>(corners1[corner]) * 3;
            const std::size_t base2 = static_cast<std::size_t>(corners2[corner]) * 3;
            for (int axis = 0; axis < 3; axis++)
            {
                local[base1 + axis] += g.force1[corner][axis];
                local[base2 + axis] += g.force2[corner][axis];
            }
        }
        // Half to each face, as the tether does; atomics because neighbouring
        // edges share faces.
#pragma omp atomic
        faces[edge.face[0]].energy.energyRegularization += 0.5 * g.energy;
#pragma omp atomic
        faces[edge.face[1]].energy.energyRegularization += 0.5 * g.energy;
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
