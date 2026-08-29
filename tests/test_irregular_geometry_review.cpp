#include "test_irregular_geometry_review.hpp"

namespace
{
/// A bowl-shaped disk with one extraordinary corner at its centre.
struct Disk
{
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;
};

Disk build_disk(int valence, double bowl)
{
    const CanonicalPatch patch = build_canonical_patch(valence);
    Disk disk;
    for (int v = 0; v < patch.nVertices; v++)
    {
        const double x = patch.layout[v][0] - 1.5;
        const double y = patch.layout[v][1] - 1.5;
        disk.vertices.push_back({x, y, bowl * (x * x + y * y)});
    }
    for (const std::array<int, 3> &f : patch.faces)
        disk.faces.push_back({f[0], f[1], f[2]});
    return disk;
}

int irregular_face_of(const Mesh &mesh, int valence)
{
    for (std::size_t i = 0; i < mesh.faces.size(); i++)
        if (static_cast<int>(mesh.faces[i].oneRingVertices.size()) == valence + 6)
            return static_cast<int>(i);
    return -1;
}

/**
 * @brief Area of a REGULAR 12-point patch by direct fine quadrature.
 *
 * Independent of the subdivision machinery entirely: it evaluates the regular
 * shape functions at many parameter points and integrates, rather than tiling
 * the domain with subdivided children. This is the only check available that
 * does not go through the thing it is checking.
 */
double area_by_fine_quadrature(const Matrix &control, int nDivisions)
{
    double area = 0.0;
    const double step = 1.0 / nDivisions;
    // Midpoint rule over the barycentric triangle, refined uniformly.
    for (int i = 0; i < nDivisions; i++)
    {
        for (int j = 0; j < nDivisions - i; j++)
        {
            // Two sub-triangles per cell, except the last column.
            const int nSub = (j < nDivisions - i - 1) ? 2 : 1;
            for (int s = 0; s < nSub; s++)
            {
                const double v = (s == 0) ? (i + 1.0 / 3.0) * step : (i + 2.0 / 3.0) * step;
                const double w = (s == 0) ? (j + 1.0 / 3.0) * step : (j + 2.0 / 3.0) * step;
                if (v + w >= 1.0)
                    continue;
                Matrix vwu(1, 3);
                vwu.set(0, 0, v);
                vwu.set(0, 1, w);
                vwu.set(0, 2, 1.0 - v - w);
                Matrix sf(7, 12);
                get_shapefunction(vwu, sf);
                const Matrix a1 = sf.get_row(1) * control;
                const Matrix a2 = sf.get_row(2) * control;
                area += 0.5 * step * step * cross_row(a1, a2).calculate_norm();
            }
        }
    }
    return area;
}

Matrix regular_lattice_control(double bowl)
{
    const int lattice[12][2] = {{0, 2},  {1, 2},  {0, 1}, {1, 1},  {2, 1},  {0, 0},
                               {1, 0},  {2, 0},  {3, 0}, {1, -1}, {2, -1}, {3, -1}};
    Matrix control(12, 3);
    for (int k = 0; k < 12; k++)
    {
        const double x = lattice[k][0] + 0.5 * lattice[k][1];
        const double y = lattice[k][1] * std::sqrt(3.0) / 2.0;
        control.set(k, 0, x);
        control.set(k, 1, y);
        control.set(k, 2, bowl * (x * x + y * y));
    }
    return control;
}
} // namespace

/**
 * @brief The row table's tiling agrees with direct quadrature of the same patch.
 *
 * At valence 6 the patch is regular, so its area can be integrated directly
 * from the shape functions at arbitrarily many parameter points -- no
 * subdivision involved. Tiling the same patch with subdivided children and
 * integrating each with the 3-point rule must converge to the same number.
 *
 * This is the one check in the suite that does not validate the subdivision
 * machinery against itself, and it separates the two error sources: the
 * children tile 1 - 4^-D of the domain (truncation), and each child is
 * integrated with a 3-point rule (quadrature order).
 */
TEST(IrregularGeometryReviewTest, TiledAreaConvergesToDirectQuadrature)
{
    Param param;
    std::vector<Matrix> shapeFunctions;
    get_gauss_quadrature_weight_VWU(param.gaussQuadratureN, param.VWU, param.gaussQuadratureCoeff);
    get_shapefunction_vector(param.VWU, shapeFunctions);

    const Matrix control = regular_lattice_control(0.10);
    const double reference = area_by_fine_quadrature(control, 240);
    ASSERT_GT(reference, 0.0);

    std::printf("\n  direct fine quadrature (240 divisions): %.10f\n", reference);
    std::printf("  D   tiled area      rel. error   truncated fraction\n");

    double previousError = 1.0;
    double worstRatio = 0.0;
    for (int depth = 1; depth <= 8; depth++)
    {
        IrregularPatchRowTable table;
        table.build(shapeFunctions, depth, DepthPolicy::Uniform);
        double tiled = 0.0;
        for (int d = 0; d < depth; d++)
        {
            for (int c = 0; c < kRegularChildrenPerStep; c++)
            {
                for (int q = 0; q < table.nSamples(); q++)
                {
                    const Matrix &rows = table.rows(6, d, c, q);
                    const Matrix a1 = rows.get_row(1) * control;
                    const Matrix a2 = rows.get_row(2) * control;
                    tiled += 0.5 * param.gaussQuadratureCoeff(q, 0) *
                             cross_row(a1, a2).calculate_norm();
                }
            }
        }
        const double error = std::abs(tiled - reference) / reference;
        std::printf("  %d   %.10f  %10.3e  %10.3e\n", depth, tiled, error,
                    table.truncated_fraction());
        // Monotone improvement against an independent reference.
        EXPECT_LT(error, previousError) << "depth " << depth;
        previousError = error;
        if (depth >= 3)
            worstRatio = std::max(worstRatio, std::abs(error / table.truncated_fraction() - 1.0));
    }

    // The strong statement, and the reason this test exists: measured against
    // a reference that never touches subdivision, the relative area error IS
    // the truncated parameter fraction, to about one percent. Two things
    // follow. The tiling is exact -- the children cover the domain once, with
    // no double counting and no gap beyond the sliver. And the 3-point rule
    // contributes nothing measurable at this depth, so the error budget is
    // truncation alone, which is what makes choosing a depth meaningful.
    //
    // An absolute threshold here would be wrong: no depth can do better than
    // its own truncated fraction, which is 1.5e-5 at D=8.
    EXPECT_LT(worstRatio, 0.05) << "area error should track the truncated fraction";
}

/**
 * @brief Rigid motion must not change the geometry or the energy.
 *
 * The plan lists this as one of the tests that bind, and nothing has exercised
 * it on an irregular patch. Area and bending energy are invariant under
 * translation and rotation; a subdivision matrix with a mis-weighted row would
 * still sum to 1 and still be local, but would generally break this.
 */
TEST(IrregularGeometryReviewTest, IrregularGeometryIsRigidMotionInvariant)
{
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        const Disk disk = build_disk(valence, 0.12);

        auto measure = [&](bool moved) {
            Param param;
            param.VERBOSE_MODE = false;
            param.boundaryCondition = BoundaryType::Fixed;
            param.uVol = 0.0;
            Mesh mesh(param);
            mesh.setup_from_vertices_faces(disk.vertices, disk.faces);

            if (moved)
            {
                // Rotate by 0.7 rad about z, then translate a long way.
                const double c = std::cos(0.7);
                const double s = std::sin(0.7);
                for (Vertex &vertex : mesh.vertices)
                {
                    const double x = vertex.coord(0, 0);
                    const double y = vertex.coord(1, 0);
                    const double z = vertex.coord(2, 0);
                    vertex.coord.set(0, 0, c * x - s * y + 137.0);
                    vertex.coord.set(1, 0, s * x + c * y - 59.0);
                    vertex.coord.set(2, 0, z + 23.0);
                }
            }

            mesh.calculate_element_area_volume();
            double area = 0.0;
            double volume = 0.0;
            mesh.sum_membrane_area_and_volume(area, volume);
            mesh.param.area0 = area;
            mesh.param.vol0 = 0.0;
            mesh.update_previous_coord_for_vertex();
            mesh.update_reference_coord_from_previous_coord();
            mesh.Compute_Energy_And_Force();
            const int face = irregular_face_of(mesh, valence);
            return std::make_pair(mesh.faces[face].elementArea,
                                  mesh.faces[face].energy.energyCurvature);
        };

        const auto rest = measure(false);
        const auto moved = measure(true);
        EXPECT_NEAR(moved.first, rest.first, 1e-11 * rest.first)
            << "N=" << valence << ": area changed under rigid motion";
        EXPECT_NEAR(moved.second, rest.second, 1e-10 * rest.second)
            << "N=" << valence << ": bending energy changed under rigid motion";
    }
}

/**
 * @brief Internal forces on an isolated patch must sum to zero.
 *
 * Newton's third law applied to the discretisation: the bending and area
 * forces a face exerts on its own control points are internal, so their vector
 * sum vanishes. It is independent of every other check here -- a row table
 * whose derivative rows failed to sum to zero would produce a net force out of
 * nothing, and no amount of self-consistency would reveal it.
 */
TEST(IrregularGeometryReviewTest, InternalForcesOnAnIrregularPatchSumToZero)
{
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        const Disk disk = build_disk(valence, 0.12);

        Param param;
        param.VERBOSE_MODE = false;
        param.boundaryCondition = BoundaryType::Fixed;
        param.uVol = 0.0;
        Mesh mesh(param);
        mesh.setup_from_vertices_faces(disk.vertices, disk.faces);
        mesh.calculate_element_area_volume();
        double area = 0.0;
        double volume = 0.0;
        mesh.sum_membrane_area_and_volume(area, volume);
        mesh.param.area0 = 0.5 * area; // hold the area constraint away from rest
        mesh.param.vol0 = 0.0;
        mesh.update_previous_coord_for_vertex();
        mesh.update_reference_coord_from_previous_coord();
        mesh.Compute_Energy_And_Force();

        double netBend[3] = {0.0, 0.0, 0.0};
        double netArea[3] = {0.0, 0.0, 0.0};
        double scale = 0.0;
        for (const Vertex &vertex : mesh.vertices)
        {
            for (int k = 0; k < 3; k++)
            {
                const double b = vertex.force.forceCurvature(k, 0);
                const double a = vertex.force.forceArea(k, 0);
                netBend[k] += b;
                netArea[k] += a;
                scale = std::max(scale, std::max(std::abs(b), std::abs(a)));
            }
        }
        ASSERT_GT(scale, 0.0) << "N=" << valence << ": no force to check";
        for (int k = 0; k < 3; k++)
        {
            EXPECT_LT(std::abs(netBend[k]), 1e-9 * scale)
                << "N=" << valence << " axis " << k << ": net bending force is not zero";
            EXPECT_LT(std::abs(netArea[k]), 1e-9 * scale)
                << "N=" << valence << " axis " << k << ": net area force is not zero";
        }
    }
}
