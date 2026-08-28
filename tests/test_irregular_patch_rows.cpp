#include "test_irregular_patch_rows.hpp"

namespace
{
/// The regular 7x12 shape functions at the 3-point Gauss rule, and its weights.
struct RegularQuadrature
{
    std::vector<Matrix> shapeFunctions;
    Matrix coefficients;
};

RegularQuadrature build_regular_quadrature()
{
    Param param;
    RegularQuadrature quadrature;
    get_gauss_quadrature_weight_VWU(param.gaussQuadratureN, param.VWU, param.gaussQuadratureCoeff);
    get_shapefunction_vector(param.VWU, quadrature.shapeFunctions);
    quadrature.coefficients = param.gaussQuadratureCoeff;
    return quadrature;
}

/**
 * @brief The regular 12-point patch on a true triangular lattice.
 *
 * The generator's own layout is only used for ordering and is not metrically a
 * lattice, so the geometry tests need proper coordinates. Lattice basis
 * e1 = (1,0), e2 = (0.5, sqrt(3)/2), with the standard patch numbering:
 *
 *        d1(0,2)  d2(1,2)
 *     d3(0,1)  d4(1,1)  d5(2,1)
 *   d6(0,0)  d7(1,0)  d8(2,0)  d9(3,0)
 *     d10(1,-1) d11(2,-1) d12(3,-1)
 */
Matrix regular_lattice_patch()
{
    const int lattice[12][2] = {{0, 2},  {1, 2},  {0, 1}, {1, 1},  {2, 1},  {0, 0},
                               {1, 0},  {2, 0},  {3, 0}, {1, -1}, {2, -1}, {3, -1}};
    Matrix control(12, 3);
    for (int k = 0; k < 12; k++)
    {
        const double i = lattice[k][0];
        const double j = lattice[k][1];
        control.set(k, 0, i + 0.5 * j);
        control.set(k, 1, j * std::sqrt(3.0) / 2.0);
        control.set(k, 2, 0.0);
    }
    return control;
}

/// Area contributed by one set of rows, using the same rule as the mesh.
double area_from_rows(const Matrix &rows, const Matrix &control, double coefficient)
{
    const Matrix a1 = rows.get_row(1) * control;
    const Matrix a2 = rows.get_row(2) * control;
    const Matrix a3 = cross_row(a1, a2);
    return 0.5 * coefficient * a3.calculate_norm();
}
} // namespace

/**
 * @brief Position rows average, derivative rows cancel.
 *
 * Row 0 of a shape function is a partition of unity, so it must sum to 1; the
 * six derivative rows differentiate a constant and must sum to 0. Holds for
 * every valence, depth, child and sample, and catches a mask that survived
 * composition with the wrong weights.
 */
TEST(IrregularPatchRowTableTest, PositionRowsSumToOneAndDerivativeRowsToZero)
{
    const RegularQuadrature quadrature = build_regular_quadrature();
    IrregularPatchRowTable table;
    table.build(quadrature.shapeFunctions, 4);

    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        for (int d = 0; d < table.depth(); d++)
        {
            for (int c = 0; c < kRegularChildrenPerStep; c++)
            {
                for (int q = 0; q < table.nSamples(); q++)
                {
                    const Matrix &rows = table.rows(valence, d, c, q);
                    ASSERT_EQ(rows.nrow(), 7);
                    ASSERT_EQ(rows.ncol(), valence + 6);

                    for (int r = 0; r < 7; r++)
                    {
                        double sum = 0.0;
                        for (int j = 0; j < rows.ncol(); j++)
                            sum += rows(r, j);
                        const double expected = (r == 0) ? 1.0 : 0.0;
                        EXPECT_NEAR(sum, expected, 1e-12)
                            << "N=" << valence << " depth=" << d << " child=" << c
                            << " sample=" << q << " row=" << r;
                    }
                }
            }
        }
    }
}

/**
 * @brief The children tile the parameter domain exactly once, minus the sliver.
 *
 * This is the gate on the whole work package, and the test that holds the line
 * against a spurious `4^-d` rescaling.
 *
 * On a planar regular control net the Loop limit surface is an affine image of
 * the parameter domain, so |a_1 x a_2| is constant and the 3-point rule is
 * exact rather than approximate. The three regular children peeled off at each
 * step cover three quarters of what remains, so after D steps the table covers
 * exactly `1 - 4^-D` of the domain -- and the area it reports must equal that
 * fraction of the directly evaluated area, to round-off and not merely to
 * convergence.
 *
 * An extra `4^-d` factor anywhere in the composition changes that fraction and
 * fails here immediately; so does a missing one, if the rows were ever assumed
 * to be in parent parameters.
 */
TEST(IrregularPatchRowTableTest, PlanarRegularPatchTilesTheDomainExactly)
{
    const RegularQuadrature quadrature = build_regular_quadrature();
    const Matrix control = regular_lattice_patch();

    // Direct evaluation over the whole patch, the way the regular kernel does.
    double directArea = 0.0;
    for (int q = 0; q < static_cast<int>(quadrature.shapeFunctions.size()); q++)
    {
        directArea += area_from_rows(quadrature.shapeFunctions[q], control,
                                     quadrature.coefficients(q, 0));
    }
    ASSERT_GT(directArea, 0.0);

    for (int depth = 1; depth <= 5; depth++)
    {
        IrregularPatchRowTable table;
        table.build(quadrature.shapeFunctions, depth);

        double tabulatedArea = 0.0;
        for (int d = 0; d < depth; d++)
        {
            for (int c = 0; c < kRegularChildrenPerStep; c++)
            {
                for (int q = 0; q < table.nSamples(); q++)
                {
                    tabulatedArea += area_from_rows(table.rows(6, d, c, q), control,
                                                    quadrature.coefficients(q, 0));
                }
            }
        }

        const double covered = 1.0 - std::pow(0.25, depth);
        EXPECT_NEAR(tabulatedArea, covered * directArea, 1e-12 * directArea)
            << "depth " << depth << ": expected the children to tile " << covered
            << " of the domain";
        EXPECT_NEAR(table.truncated_fraction(), 1.0 - covered, 1e-15);
    }
}

/**
 * @brief A planar control net stays planar through every tabulated row.
 *
 * Independent of the area identity: it checks positions rather than the
 * metric, so a composition that preserves total area while displacing points
 * off the plane still fails.
 */
TEST(IrregularPatchRowTableTest, PlanarControlNetStaysPlanarAtEveryDepth)
{
    const RegularQuadrature quadrature = build_regular_quadrature();
    IrregularPatchRowTable table;
    table.build(quadrature.shapeFunctions, 4);

    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        // Put the control points on the plane z = 2x - 3y + 1.
        const CanonicalPatch patch = build_canonical_patch(valence);
        Matrix control(valence + 6, 3);
        for (int v = 0; v < patch.nVertices; v++)
        {
            const double x = patch.layout[v][0];
            const double y = patch.layout[v][1];
            control.set(v, 0, x);
            control.set(v, 1, y);
            control.set(v, 2, 2.0 * x - 3.0 * y + 1.0);
        }

        for (int d = 0; d < table.depth(); d++)
        {
            for (int c = 0; c < kRegularChildrenPerStep; c++)
            {
                for (int q = 0; q < table.nSamples(); q++)
                {
                    const Matrix point = table.rows(valence, d, c, q).get_row(0) * control;
                    const double expected = 2.0 * point(0, 0) - 3.0 * point(0, 1) + 1.0;
                    EXPECT_NEAR(point(0, 2), expected, 1e-10)
                        << "N=" << valence << " depth=" << d << " child=" << c;
                }
            }
        }
    }
}

/// Bad keys and a bad build are errors, and the table stays small.
TEST(IrregularPatchRowTableTest, RejectsBadInputAndStaysSmall)
{
    const RegularQuadrature quadrature = build_regular_quadrature();

    IrregularPatchRowTable table;
    EXPECT_THROW(table.build(quadrature.shapeFunctions, 0), std::invalid_argument);
    EXPECT_THROW(table.build({}, 4), std::invalid_argument);
    EXPECT_TRUE(table.empty());

    table.build(quadrature.shapeFunctions, kDefaultIrregularDepth);
    EXPECT_FALSE(table.empty());
    EXPECT_THROW(table.rows(3, 0, 0, 0), std::out_of_range);
    EXPECT_THROW(table.rows(9, 0, 0, 0), std::out_of_range);
    EXPECT_THROW(table.rows(6, table.depth(), 0, 0), std::out_of_range);
    EXPECT_THROW(table.rows(6, 0, kRegularChildrenPerStep, 0), std::out_of_range);
    EXPECT_THROW(table.rows(6, 0, 0, table.nSamples()), std::out_of_range);

    // The plan budgets under 200 kB; it is a startup-time constant shared
    // across threads, so this only needs to stay in that neighbourhood.
    EXPECT_LT(table.memory_bytes(), 200u * 1024u);
}

namespace
{
/// Area and signed volume accumulated the way the mesh does, from 7 x K rows.
void accumulate_area_volume(const Matrix &rows, const Matrix &control, double coefficient,
                            double &area, double &volume)
{
    const Matrix x = rows.get_row(0) * control;
    const Matrix a3 = cross_row(rows.get_row(1) * control, rows.get_row(2) * control);
    area += 0.5 * coefficient * a3.calculate_norm();
    volume += (1.0 / 6.0) * coefficient * dot_row(x, a3);
}

/// A deliberately non-planar control net of K = N+6 points, fixed for stability.
Matrix irregular_control_net_for(int valence)
{
    const int K = valence + 6;
    Matrix control(K, 3);
    for (int k = 0; k < K; k++)
    {
        const double t = k;
        control.set(k, 0, std::cos(0.7 * t) * (1.0 + 0.1 * t));
        control.set(k, 1, std::sin(0.7 * t) * (1.0 + 0.1 * t));
        control.set(k, 2, 0.3 * std::sin(1.3 * t) + 0.05 * t * t);
    }
    return control;
}

Matrix irregular_control_net() { return irregular_control_net_for(5); }
} // namespace

/**
 * @brief The collapsed table must reproduce the recursion it replaces.
 *
 * calculate_element_area_volume() already evaluates an 11-control face by
 * explicit recursion -- subdivide with M, evaluate the three regular children
 * through M1..M3, descend into M4, repeat. The row table is that same
 * recursion folded flat, so the two must agree on a real, non-planar control
 * net.
 *
 * The recursion side is driven by get_subdivision_matrices_oracle(), the literal
 * matrices, while the table side is driven by the generated ones. So this
 * checks the generator and the collapse together, against the production code
 * path, on geometry rather than on matrix entries -- which the N=5 parity test
 * cannot do.
 */
TEST(IrregularPatchRowTableTest, CollapsedTableReproducesTheExplicitRecursion)
{
    const RegularQuadrature quadrature = build_regular_quadrature();
    const Matrix control = irregular_control_net();

    Matrix oracleM, oracleM1, oracleM2, oracleM3, oracleM4;
    get_subdivision_matrices_oracle(oracleM, oracleM1, oracleM2, oracleM3, oracleM4);

    for (int depth = 1; depth <= 6; depth++)
    {
        // The explicit recursion, exactly as calculate_element_area_volume()
        // performs it for the 11-control case.
        double recursionArea = 0.0;
        double recursionVolume = 0.0;
        Matrix descending = control;
        for (int d = 0; d < depth; d++)
        {
            const Matrix subdivided = oracleM * descending;
            const Matrix *regularChildren[3] = {&oracleM1, &oracleM2, &oracleM3};
            for (int c = 0; c < 3; c++)
            {
                const Matrix childControl = (*regularChildren[c]) * subdivided;
                for (int q = 0; q < static_cast<int>(quadrature.shapeFunctions.size()); q++)
                {
                    accumulate_area_volume(quadrature.shapeFunctions[q], childControl,
                                           quadrature.coefficients(q, 0), recursionArea,
                                           recursionVolume);
                }
            }
            descending = oracleM4 * subdivided;
        }

        // The same thing, read out of the collapsed table.
        IrregularPatchRowTable table;
        table.build(quadrature.shapeFunctions, depth);
        double tableArea = 0.0;
        double tableVolume = 0.0;
        for (int d = 0; d < depth; d++)
        {
            for (int c = 0; c < kRegularChildrenPerStep; c++)
            {
                for (int q = 0; q < table.nSamples(); q++)
                {
                    accumulate_area_volume(table.rows(5, d, c, q), control,
                                           quadrature.coefficients(q, 0), tableArea, tableVolume);
                }
            }
        }

        ASSERT_GT(recursionArea, 0.0) << "depth " << depth;
        EXPECT_NEAR(tableArea, recursionArea, 1e-12 * std::abs(recursionArea))
            << "area disagrees at depth " << depth;
        EXPECT_NEAR(tableVolume, recursionVolume, 1e-12 * std::abs(recursionVolume))
            << "volume disagrees at depth " << depth;
    }
}

/**
 * @brief Refining the depth converges, and the rate depends on the valence.
 *
 * The dropped sliver is `4^-D` of the *parameter* domain at every valence, but
 * the surface area inside it is not: near an extraordinary corner the metric
 * shrinks like the subdivision matrix's subdominant eigenvalue, which is a
 * function of N. Measured on a fixed non-planar net, successive increments
 * fall by roughly
 *
 *     N=4  7.3x     N=5  5.0x     N=6  4.0x     N=7  3.5x     N=8  3.2x
 *
 * per level -- exactly 4 at N=6, where the patch is regular and the metric is
 * uniform, and progressively slower above it.
 *
 * The consequence is that one depth does not serve every valence equally, so
 * WP6 has to choose D per valence rather than globally. This test pins the
 * qualitative shape of that -- monotone convergence everywhere, exactly 4x at
 * N=6, and slower at N=8 than at N=4 -- so a regression in the recursion shows
 * up as a changed rate rather than only as a changed number.
 */
TEST(IrregularPatchRowTableTest, ConvergenceRateDependsOnValence)
{
    const RegularQuadrature quadrature = build_regular_quadrature();

    auto area_at_depth = [&](int valence, int depth) {
        const Matrix control = irregular_control_net_for(valence);
        IrregularPatchRowTable table;
        table.build(quadrature.shapeFunctions, depth);
        double area = 0.0;
        double volume = 0.0;
        for (int d = 0; d < depth; d++)
            for (int c = 0; c < kRegularChildrenPerStep; c++)
                for (int q = 0; q < table.nSamples(); q++)
                    accumulate_area_volume(table.rows(valence, d, c, q), control,
                                           quadrature.coefficients(q, 0), area, volume);
        return area;
    };

    std::vector<double> ratio(kMaxIrregularValence + 1, 0.0);
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        double previous = area_at_depth(valence, 5);
        const double stepA = area_at_depth(valence, 6) - previous;
        previous += stepA;
        const double stepB = area_at_depth(valence, 7) - previous;
        previous += stepB;
        const double stepC = area_at_depth(valence, 8) - previous;

        // Monotone: each level only adds surface, since the children tile a
        // strictly larger fraction of the domain.
        EXPECT_GT(stepA, 0.0) << "valence " << valence;
        EXPECT_GT(stepB, 0.0) << "valence " << valence;
        EXPECT_GT(stepC, 0.0) << "valence " << valence;
        // And shrinking, so the truncation error actually goes away.
        EXPECT_LT(stepB, stepA) << "valence " << valence;
        EXPECT_LT(stepC, stepB) << "valence " << valence;

        ratio[valence] = stepB / stepC;
    }

    // A regular patch has a uniform metric, so the area in the dropped sliver
    // falls exactly as the parameter fraction does.
    EXPECT_NEAR(ratio[6], 4.0, 0.1);
    // Higher valence converges more slowly; lower, faster.
    EXPECT_LT(ratio[8], ratio[6]);
    EXPECT_GT(ratio[4], ratio[6]);
}

/**
 * @brief An irregular face is now evaluated, not silently zeroed.
 *
 * The canonical patch for valence 5 is itself the fixture: its centre triangle
 * (d4, d7, d8) has one valence-5 corner and two valence-6 corners, so it is
 * admitted with an 11-point one-ring, while every other face touches the
 * boundary and is skipped.
 *
 * Before WP4 this face reached element_energy_force_regular() with a 7x12
 * shape function against an 11-wide control net -- a dimension mismatch inside
 * a gsl_blas_dgemm wrapper that discards its return code. Now it is driven
 * through the same kernel as a regular face, with rows of its own width.
 */
TEST(IrregularPatchRowTableTest, IrregularFaceIsEvaluatedEndToEnd)
{
    const CanonicalPatch patch = build_canonical_patch(5);

    std::vector<std::vector<double>> vertices;
    for (int v = 0; v < patch.nVertices; v++)
    {
        const double x = patch.layout[v][0];
        const double y = patch.layout[v][1];
        // A bowl, so the bending energy has something to see.
        vertices.push_back({x, y, 0.15 * (x * x + y * y)});
    }
    std::vector<std::vector<int>> faces;
    for (const std::array<int, 3> &face : patch.faces)
        faces.push_back({face[0], face[1], face[2]});

    Param param;
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.uVol = 0.0;
    Mesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_from_vertices_faces(vertices, faces));

    // Exactly one face has a complete one-ring, and it is the irregular one.
    int nIrregular = 0;
    int irregularFace = -1;
    for (std::size_t i = 0; i < mesh.faces.size(); i++)
    {
        const std::size_t width = mesh.faces[i].oneRingVertices.size();
        if (width == 0)
            continue;
        EXPECT_EQ(width, 11u) << "face " << i;
        nIrregular++;
        irregularFace = static_cast<int>(i);
    }
    ASSERT_EQ(nIrregular, 1);

    mesh.calculate_element_area_volume();
    double area = 0.0;
    double volume = 0.0;
    mesh.sum_membrane_area_and_volume(area, volume);
    mesh.param.area0 = area;
    mesh.param.vol0 = 0.0;
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    ASSERT_NO_THROW(mesh.Compute_Energy_And_Force());

    // The irregular face now carries real curvature energy rather than the
    // zero initializer it used to be stored with.
    EXPECT_GT(mesh.faces[irregularFace].energy.energyCurvature, 0.0);
    EXPECT_TRUE(std::isfinite(mesh.faces[irregularFace].energy.energyCurvature));
    EXPECT_TRUE(std::isfinite(mesh.faces[irregularFace].meanCurvature));

    // And it exerts force on its own control points.
    double maxForce = 0.0;
    for (int node : mesh.faces[irregularFace].oneRingVertices)
    {
        for (int k = 0; k < 3; k++)
        {
            const double component = mesh.vertices[node].force.forceCurvature(k, 0);
            EXPECT_TRUE(std::isfinite(component)) << "vertex " << node << " axis " << k;
            maxForce = std::max(maxForce, std::abs(component));
        }
    }
    EXPECT_GT(maxForce, 0.0);
}

/**
 * @brief Geometry and energy/force read the same rows for the same face.
 *
 * WP5's gate. Before it, calculate_element_area_volume() ran its own explicit
 * subdivision recursion over param.subMatrix, bounded by an uninitialized
 * param.subDivideTimes, while energy and force went through the row table. So
 * the area a face reported and the area its bending energy was computed from
 * came from different code, integrated to a different depth.
 *
 * Now both go through IrregularPatchRowTable, and the face's stored
 * elementArea must equal the area those same rows give when summed
 * independently.
 */
TEST(IrregularPatchRowTableTest, GeometryAndEnergyReadTheSameRows)
{
    const CanonicalPatch patch = build_canonical_patch(5);

    std::vector<std::vector<double>> vertices;
    for (int v = 0; v < patch.nVertices; v++)
    {
        const double x = patch.layout[v][0];
        const double y = patch.layout[v][1];
        vertices.push_back({x, y, 0.15 * (x * x + y * y)});
    }
    std::vector<std::vector<int>> faces;
    for (const std::array<int, 3> &face : patch.faces)
        faces.push_back({face[0], face[1], face[2]});

    Param param;
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.uVol = 0.0;
    Mesh mesh(param);
    mesh.setup_from_vertices_faces(vertices, faces);
    mesh.calculate_element_area_volume();

    int irregularFace = -1;
    for (std::size_t i = 0; i < mesh.faces.size(); i++)
    {
        if (mesh.faces[i].oneRingVertices.size() == 11u)
            irregularFace = static_cast<int>(i);
    }
    ASSERT_GE(irregularFace, 0);

    // The same rows, summed here rather than inside the mesh.
    Matrix control(11, 3);
    for (int k = 0; k < 11; k++)
    {
        const int node = mesh.faces[irregularFace].oneRingVertices[k];
        for (int axis = 0; axis < 3; axis++)
            control.set(k, axis, mesh.vertices[node].coord(axis, 0));
    }

    double expectedArea = 0.0;
    double expectedVolume = 0.0;
    for (int d = 0; d < mesh.irregularRows.depth(); d++)
    {
        for (int c = 0; c < kRegularChildrenPerStep; c++)
        {
            for (int q = 0; q < mesh.irregularRows.nSamples(); q++)
            {
                accumulate_area_volume(mesh.irregularRows.rows(5, d, c, q), control,
                                       mesh.param.gaussQuadratureCoeff(q, 0), expectedArea,
                                       expectedVolume);
            }
        }
    }

    ASSERT_GT(expectedArea, 0.0);
    EXPECT_NEAR(mesh.faces[irregularFace].elementArea, expectedArea, 1e-12 * expectedArea);
    EXPECT_NEAR(mesh.faces[irregularFace].elementVolume, expectedVolume,
                1e-12 * std::abs(expectedVolume) + 1e-15);

    // And a boundary face, having no limit surface, contributes nothing.
    for (const Face &face : mesh.faces)
    {
        if (face.oneRingVertices.empty())
        {
            EXPECT_DOUBLE_EQ(face.elementArea, 0.0);
            EXPECT_DOUBLE_EQ(face.elementVolume, 0.0);
        }
    }
}
