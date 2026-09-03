#include "test_subdivision_matrices.hpp"

/**
 * @brief The generator must reproduce the literal valence-5 matrix exactly.
 *
 * get_subdivision_matrices_oracle() is a direct port of the reference POC and is the
 * only independently-derived subdivision data in the tree, so it is the oracle
 * for everything the generator produces.
 */
TEST(SubdivisionMatrixTest, AbarMatchesTheLiteralOracleAtValence5)
{
    Matrix oracle, m1, m2, m3, m4;
    get_subdivision_matrices_oracle(oracle, m1, m2, m3, m4);

    const SubdivisionMatrices generated = build_subdivision_matrices(5);

    ASSERT_EQ(generated.abar.nrow(), oracle.nrow());
    ASSERT_EQ(generated.abar.ncol(), oracle.ncol());

    for (int i = 0; i < oracle.nrow(); i++)
    {
        for (int j = 0; j < oracle.ncol(); j++)
        {
            EXPECT_DOUBLE_EQ(generated.abar(i, j), oracle(i, j))
                << "row " << i << " col " << j;
        }
    }
}

/**
 * @brief Every row of abar must sum to 1.
 *
 * Affine invariance: a subdivision step is a weighted average, so translating
 * the control net must translate the subdivided net by the same amount and
 * nothing else. A row that does not sum to 1 is a mask transcribed wrong, and
 * this catches it for the valences that have no literal oracle to compare
 * against.
 */
TEST(SubdivisionMatrixTest, AbarRowsSumToOneForEverySupportedValence)
{
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        const SubdivisionMatrices generated = build_subdivision_matrices(valence);
        ASSERT_EQ(generated.abar.nrow(), valence + 12) << "valence " << valence;
        ASSERT_EQ(generated.abar.ncol(), valence + 6) << "valence " << valence;

        for (int i = 0; i < generated.abar.nrow(); i++)
        {
            double sum = 0.0;
            for (int j = 0; j < generated.abar.ncol(); j++)
            {
                sum += generated.abar(i, j);
            }
            EXPECT_NEAR(sum, 1.0, 1e-14) << "valence " << valence << " row " << i;
        }
    }
}

/**
 * @brief Weights must be non-negative, and each row must be local.
 *
 * A Loop vertex mask touches the vertex and its ring; an edge mask touches
 * four points. Anything wider means the mask picked up a point it should not
 * see -- the failure a hand-written literal makes easy and invisible.
 */
TEST(SubdivisionMatrixTest, AbarRowsAreLocalAndNonNegative)
{
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        const SubdivisionMatrices generated = build_subdivision_matrices(valence);
        int nEdgeRows = 0;
        int nVertexRows = 0;
        for (int i = 0; i < generated.abar.nrow(); i++)
        {
            int nNonZero = 0;
            for (int j = 0; j < generated.abar.ncol(); j++)
            {
                EXPECT_GE(generated.abar(i, j), 0.0) << "valence " << valence << " row " << i;
                if (generated.abar(i, j) != 0.0)
                    nNonZero++;
            }
            if (nNonZero == 4)
                nEdgeRows++;
            else
                nVertexRows++;
        }
        // Three vertex points -- d4, d7, d8 -- and an edge point for every
        // other child.
        EXPECT_EQ(nVertexRows, 3) << "valence " << valence;
        EXPECT_EQ(nEdgeRows, valence + 9) << "valence " << valence;
    }
}

/**
 * @brief A planar control net must subdivide to a planar net.
 *
 * Independent of the row sums: it exercises the masks against actual
 * coordinates rather than against their own weights, so a mask that sums to 1
 * while mixing the wrong points still fails here.
 */
TEST(SubdivisionMatrixTest, PlanarControlNetStaysPlanar)
{
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        const CanonicalPatch patch = build_canonical_patch(valence);
        const SubdivisionMatrices generated = build_subdivision_matrices(valence);

        // Lay the control points on the plane z = 2x - 3y + 1.
        Matrix control(valence + 6, 3);
        for (int v = 0; v < patch.nVertices; v++)
        {
            const double x = patch.layout[v][0];
            const double y = patch.layout[v][1];
            control.set(v, 0, x);
            control.set(v, 1, y);
            control.set(v, 2, 2.0 * x - 3.0 * y + 1.0);
        }

        const Matrix subdivided = generated.abar * control;
        for (int i = 0; i < subdivided.nrow(); i++)
        {
            const double expected =
                2.0 * subdivided(i, 0) - 3.0 * subdivided(i, 1) + 1.0;
            EXPECT_NEAR(subdivided(i, 2), expected, 1e-12)
                << "valence " << valence << " child " << i;
        }
    }
}

/// Out-of-range valence is an error, not a fallback.
TEST(SubdivisionMatrixTest, ValenceOutsideTheSupportedRangeThrows)
{
    EXPECT_THROW(build_subdivision_matrices(3), std::invalid_argument);
    EXPECT_THROW(build_subdivision_matrices(9), std::invalid_argument);
    EXPECT_NO_THROW(build_subdivision_matrices(kMinIrregularValence));
    EXPECT_NO_THROW(build_subdivision_matrices(kMaxIrregularValence));
}

/**
 * @brief All four selection matrices must match the literal oracle at N=5.
 *
 * This is the gate the whole work package rests on: the generator has to
 * reproduce every one of the POC's matrices exactly before it can be trusted
 * at a valence that has no oracle.
 */
TEST(SubdivisionMatrixTest, SelectionMatricesMatchTheLiteralOracleAtValence5)
{
    Matrix oracleM, oracleM1, oracleM2, oracleM3, oracleM4;
    get_subdivision_matrices_oracle(oracleM, oracleM1, oracleM2, oracleM3, oracleM4);

    const SubdivisionMatrices generated = build_subdivision_matrices(5);

    const std::vector<std::pair<const Matrix *, const Matrix *>> pairs{
        {&generated.p1, &oracleM1},
        {&generated.p2, &oracleM2},
        {&generated.p3, &oracleM3},
        {&generated.p4, &oracleM4}};

    const char *names[] = {"p1", "p2", "p3", "p4"};
    for (std::size_t k = 0; k < pairs.size(); k++)
    {
        const Matrix &got = *pairs[k].first;
        const Matrix &want = *pairs[k].second;
        ASSERT_EQ(got.nrow(), want.nrow()) << names[k];
        ASSERT_EQ(got.ncol(), want.ncol()) << names[k];
        for (int i = 0; i < want.nrow(); i++)
        {
            for (int j = 0; j < want.ncol(); j++)
            {
                EXPECT_DOUBLE_EQ(got(i, j), want(i, j))
                    << names[k] << " row " << i << " col " << j;
            }
        }
    }
}

/**
 * @brief The selection matrices must be well-formed at every supported valence.
 *
 * There is no oracle for valences other than 5, so this pins the structural
 * properties instead: the walk completes, each child patch names distinct
 * points, and each selection is a genuine permutation-like pick (one 1 per
 * row, nothing else). Shapes follow Stam's reduction -- three regular
 * 12-point children and one irregular child of N+6.
 */
TEST(SubdivisionMatrixTest, SelectionMatricesAreWellFormedForEverySupportedValence)
{
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        const SubdivisionMatrices generated = build_subdivision_matrices(valence);
        const int nChildren = valence + 12;

        const Matrix *selections[] = {&generated.p1, &generated.p2, &generated.p3,
                                      &generated.p4};
        const int expectedRows[] = {12, 12, 12, valence + 6};
        const char *names[] = {"p1", "p2", "p3", "p4"};

        for (int k = 0; k < 4; k++)
        {
            const Matrix &selection = *selections[k];
            ASSERT_EQ(selection.nrow(), expectedRows[k]) << names[k] << " at N=" << valence;
            ASSERT_EQ(selection.ncol(), nChildren) << names[k] << " at N=" << valence;

            std::vector<int> picked;
            for (int i = 0; i < selection.nrow(); i++)
            {
                int nOnes = 0;
                int column = -1;
                for (int j = 0; j < selection.ncol(); j++)
                {
                    const double value = selection(i, j);
                    if (value == 0.0)
                        continue;
                    EXPECT_DOUBLE_EQ(value, 1.0) << names[k] << " N=" << valence << " row " << i;
                    nOnes++;
                    column = j;
                }
                EXPECT_EQ(nOnes, 1) << names[k] << " N=" << valence << " row " << i;
                picked.push_back(column);
            }

            std::sort(picked.begin(), picked.end());
            EXPECT_EQ(std::adjacent_find(picked.begin(), picked.end()), picked.end())
                << names[k] << " at N=" << valence << " names a child twice";
        }
    }
}

/**
 * @brief Composing a selection with abar preserves affine invariance.
 *
 * P * abar is the object the row table will actually carry, so the row-sum
 * property has to survive the composition, not just hold for abar alone.
 */
TEST(SubdivisionMatrixTest, ComposedRowsSumToOneForEverySupportedValence)
{
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        const SubdivisionMatrices generated = build_subdivision_matrices(valence);
        const Matrix *selections[] = {&generated.p1, &generated.p2, &generated.p3,
                                      &generated.p4};
        const char *names[] = {"p1", "p2", "p3", "p4"};

        for (int k = 0; k < 4; k++)
        {
            const Matrix composed = (*selections[k]) * generated.abar;
            ASSERT_EQ(composed.ncol(), valence + 6) << names[k] << " at N=" << valence;
            for (int i = 0; i < composed.nrow(); i++)
            {
                double sum = 0.0;
                for (int j = 0; j < composed.ncol(); j++)
                    sum += composed(i, j);
                EXPECT_NEAR(sum, 1.0, 1e-14)
                    << names[k] << " at N=" << valence << " row " << i;
            }
        }
    }
}

/**
 * @brief The recursion composes: (p4 * abar)^d stays a valid patch map.
 *
 * Each subdivision step hands the extraordinary corner to a smaller copy of
 * the same patch, so p4 * abar must map the N+6 control points back onto N+6
 * in the same basis. If it did not, the chain the row table is built from
 * would not close after the first step -- and the failure would surface much
 * later, as wrong physics rather than a shape mismatch.
 */
TEST(SubdivisionMatrixTest, IrregularChainComposesToAnyDepth)
{
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        const SubdivisionMatrices generated = build_subdivision_matrices(valence);
        const Matrix step = generated.p4 * generated.abar;
        ASSERT_EQ(step.nrow(), valence + 6) << "valence " << valence;
        ASSERT_EQ(step.ncol(), valence + 6) << "valence " << valence;

        Matrix chain = step;
        for (int depth = 1; depth <= 6; depth++)
        {
            for (int i = 0; i < chain.nrow(); i++)
            {
                double sum = 0.0;
                for (int j = 0; j < chain.ncol(); j++)
                {
                    EXPECT_GE(chain(i, j), 0.0)
                        << "valence " << valence << " depth " << depth << " row " << i;
                    sum += chain(i, j);
                }
                EXPECT_NEAR(sum, 1.0, 1e-12)
                    << "valence " << valence << " depth " << depth << " row " << i;
            }
            // The three regular children must remain extractable at every depth.
            const Matrix regular = generated.p1 * generated.abar * chain;
            ASSERT_EQ(regular.nrow(), 12) << "valence " << valence << " depth " << depth;
            ASSERT_EQ(regular.ncol(), valence + 6) << "valence " << valence << " depth " << depth;
            chain = step * chain;
        }
    }
}
