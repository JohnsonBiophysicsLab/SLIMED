#include "test_core_coverage.hpp"

namespace
{

/**
 * @brief Return the text after "key: " on the first line that starts with
 * that key, or "<missing>" when no line does. Every scalar field of the Param
 * printer is written as one "key: value" line, so this reads one back.
 */
std::string printed_value(const std::string &text, const std::string &key)
{
    std::istringstream in(text);
    std::string line;
    const std::string prefix = key + ": ";
    while (std::getline(in, line))
    {
        if (line.compare(0, prefix.size(), prefix) == 0)
        {
            return line.substr(prefix.size());
        }
    }
    return "<missing>";
}

/** @brief Number of lines in text that start with "key: " (or "key:" at end). */
int count_key_lines(const std::string &text, const std::string &key)
{
    std::istringstream in(text);
    std::string line;
    int count = 0;
    const std::string prefix = key + ":";
    while (std::getline(in, line))
    {
        if (line.compare(0, prefix.size(), prefix) == 0)
        {
            ++count;
        }
    }
    return count;
}

/** @brief A (3,1) column with the given entries. */
Matrix column3(double x, double y, double z)
{
    Matrix m(3, 1);
    m.set(0, 0, x);
    m.set(1, 0, y);
    m.set(2, 0, z);
    return m;
}

/** @brief Expect a (3,1) column to hold exactly (x, y, z). */
void expect_column3(const Matrix &m, double x, double y, double z)
{
    ASSERT_EQ(m.nrow(), 3);
    ASSERT_EQ(m.ncol(), 1);
    EXPECT_DOUBLE_EQ(m(0, 0), x);
    EXPECT_DOUBLE_EQ(m(1, 0), y);
    EXPECT_DOUBLE_EQ(m(2, 0), z);
}

/** @brief A Force whose every component is distinct, for the printers and sums. */
Force numbered_force()
{
    Force f;
    f.forceCurvature = column3(1, 2, 3);
    f.forceArea = column3(10, 20, 30);
    f.forceVolume = column3(100, 200, 300);
    f.forceThickness = column3(0.5, 0.25, 0.125);
    f.forceTilt = column3(-1, -2, -3);
    f.forceRegularization = column3(1000, 2000, 3000);
    f.forceHarmonicBond = column3(4, 5, 6);
    f.forceTotal = column3(7, 8, 9);
    return f;
}

} // namespace

/* ------------------------------------------------------------------------- */
/* Param printer and boundary names (src/parameters/Parameters.cpp)          */
/* ------------------------------------------------------------------------- */

TEST(CoreCovParamPrinter, DefaultParamPrintsItsDefaults)
{
    Param param;
    std::ostringstream os;
    os << param;
    const std::string text = os.str();

    // The printer starts with the verbosity flag and ends with the energies.
    EXPECT_EQ(text.compare(0, 14, "VERBOSE_MODE: "), 0);
    EXPECT_EQ(printed_value(text, "kCurv"), "83.4");
    EXPECT_EQ(printed_value(text, "gaussQuadratureN"), "2");
    EXPECT_EQ(printed_value(text, "numIterations"), "100000"); // maxIterations is printed under this key
    EXPECT_EQ(printed_value(text, "boundaryCondition"), std::to_string(static_cast<int>(BoundaryType::Periodic)));
    EXPECT_EQ(printed_value(text, "boundaryConditionX"), "Periodic");
    EXPECT_EQ(printed_value(text, "boundaryConditionY"), "Periodic");
    EXPECT_EQ(printed_value(text, "edgeTetherShape"), "flat");
    EXPECT_EQ(printed_value(text, "meshWallType"), "soft");
    EXPECT_EQ(printed_value(text, "surfaceSolver"), "dense");
    EXPECT_EQ(printed_value(text, "externalForceShape"), "disc");
    EXPECT_EQ(printed_value(text, "beadPotential"), "wall");
    EXPECT_EQ(printed_value(text, "lipidFlowMode"), "off");
    EXPECT_EQ(printed_value(text, "lipidReservoirRegion"), "none");
    EXPECT_EQ(printed_value(text, "thermalFluctuationTemperatureKelvin"), "298");
    EXPECT_EQ(printed_value(text, "thermalFluctuationStepScale"), "0.02");
    EXPECT_EQ(printed_value(text, "meshVerticesFile"), "");

    // An unallocated VWU / quadrature matrix prints as an empty block, and an
    // empty shape-function list prints nothing between its header and the
    // boundary line.
    EXPECT_NE(text.find("VWU: \n\ngaussQuadratureCoeff: \n\nshapeFunctions: \nboundaryCondition: "),
              std::string::npos);

    // Each key is printed exactly once.
    for (const char *key : {"kCurv", "uSurf", "sideX", "lFace", "timeStep", "KBT", "edgeFlipEnabled",
                            "beadRadius", "elasticBoxTension", "lipidFlipResampling", "localOutput",
                            "energy", "energyPrev"})
    {
        EXPECT_EQ(count_key_lines(text, key), 1) << key;
    }
}

TEST(CoreCovParamPrinter, PrintsTheValuesThatWereSet)
{
    Param param;
    param.kCurv = 12.5;
    param.uSurf = 250;
    param.sideX = 7;
    param.sideY = 9;
    param.lFace = 1.5;
    param.timeStep = 0.002;
    param.maxIterations = 42;
    param.boundaryCondition = BoundaryType::Mixed;
    param.boundaryConditionX = BoundaryType::Free;
    param.boundaryConditionY = BoundaryType::Fixed;
    param.meshVerticesFile = "verts.csv";
    param.meshFacesFile = "faces.csv";
    param.surfaceSolver = "sparse";
    param.lipidFlowMode = "local";
    param.edgeFlipEnabled = true;
    param.beadRadius = 3.25;
    param.energy.energyTotal = 3.5;
    param.energyPrev.energyTotal = -1.25;
    param.VWU = Matrix(std::vector<std::vector<double>>{{1, 2, 3}, {4, 5, 6}});
    param.gaussQuadratureCoeff = Matrix(std::vector<std::vector<double>>{{0.5}, {0.25}});
    param.shapeFunctions.push_back(Matrix(std::vector<std::vector<double>>{{7}, {8}}));
    param.shapeFunctions.push_back(Matrix(std::vector<std::vector<double>>{{9, 10}}));

    std::ostringstream os;
    os << param;
    const std::string text = os.str();

    EXPECT_EQ(printed_value(text, "kCurv"), "12.5");
    EXPECT_EQ(printed_value(text, "uSurf"), "250");
    EXPECT_EQ(printed_value(text, "sideX"), "7");
    EXPECT_EQ(printed_value(text, "sideY"), "9");
    EXPECT_EQ(printed_value(text, "lFace"), "1.5");
    EXPECT_EQ(printed_value(text, "timeStep"), "0.002");
    EXPECT_EQ(printed_value(text, "numIterations"), "42");
    EXPECT_EQ(printed_value(text, "boundaryCondition"), std::to_string(static_cast<int>(BoundaryType::Mixed)));
    EXPECT_EQ(printed_value(text, "boundaryConditionX"), "Free");
    EXPECT_EQ(printed_value(text, "boundaryConditionY"), "Fixed");
    EXPECT_EQ(printed_value(text, "meshVerticesFile"), "verts.csv");
    EXPECT_EQ(printed_value(text, "meshFacesFile"), "faces.csv");
    EXPECT_EQ(printed_value(text, "surfaceSolver"), "sparse");
    EXPECT_EQ(printed_value(text, "lipidFlowMode"), "local");
    EXPECT_EQ(printed_value(text, "edgeFlipEnabled"), "1");
    EXPECT_EQ(printed_value(text, "beadRadius"), "3.25");
    // Energy prints its total through std::to_string.
    EXPECT_EQ(printed_value(text, "energy"), "3.500000");
    EXPECT_EQ(printed_value(text, "energyPrev"), "-1.250000");

    // Matrices: a (2,3) VWU prints row by row, the (2,1) coefficients as one
    // line, and each shape function on its own.
    EXPECT_NE(text.find("VWU: \n1, 2, 3, \n4, 5, 6, \n\n"), std::string::npos);
    EXPECT_NE(text.find("gaussQuadratureCoeff: \n0.5, 0.25, \n"), std::string::npos);
    EXPECT_NE(text.find("shapeFunctions: \n7, 8, \n9, 10, \n\nboundaryCondition: "), std::string::npos);
}

TEST(CoreCovBoundaryNames, NamesRoundTripThroughTheParser)
{
    for (BoundaryType type : {BoundaryType::Fixed, BoundaryType::Periodic, BoundaryType::Free, BoundaryType::Mixed})
    {
        BoundaryType parsed = (type == BoundaryType::Fixed) ? BoundaryType::Free : BoundaryType::Fixed;
        EXPECT_TRUE(parse_boundary_type(boundary_type_name(type), parsed)) << boundary_type_name(type);
        EXPECT_EQ(parsed, type);
    }
    EXPECT_STREQ(boundary_type_name(BoundaryType::Fixed), "Fixed");
    EXPECT_STREQ(boundary_type_name(BoundaryType::Periodic), "Periodic");
    EXPECT_STREQ(boundary_type_name(BoundaryType::Free), "Free");
}

TEST(CoreCovBoundaryNames, OutOfRangeValueIsUnknown)
{
    EXPECT_STREQ(boundary_type_name(static_cast<BoundaryType>(99)), "unknown");
}

TEST(CoreCovBoundaryNames, ParserIgnoresCaseAndWhitespace)
{
    BoundaryType type = BoundaryType::Fixed;
    EXPECT_TRUE(parse_boundary_type(" Periodic\r", type));
    EXPECT_EQ(type, BoundaryType::Periodic);
    EXPECT_TRUE(parse_boundary_type("\tFR EE ", type));
    EXPECT_EQ(type, BoundaryType::Free);
    EXPECT_TRUE(parse_boundary_type("PerVertex", type)); // legacy alias of Mixed
    EXPECT_EQ(type, BoundaryType::Mixed);
    EXPECT_TRUE(parse_boundary_type("FIXED", type));
    EXPECT_EQ(type, BoundaryType::Fixed);

    // A rejected string leaves the output untouched.
    EXPECT_FALSE(parse_boundary_type("", type));
    EXPECT_FALSE(parse_boundary_type("fixed\n", type)); // a newline is not stripped
    EXPECT_FALSE(parse_boundary_type("periodic-x", type));
    EXPECT_EQ(type, BoundaryType::Fixed);
}

/* ------------------------------------------------------------------------- */
/* Force container (src/energy_force/Force.cpp)                              */
/* ------------------------------------------------------------------------- */

TEST(CoreCovForce, DefaultForceIsAllZeroColumns)
{
    Force f;
    for (const Matrix *m : {&f.forceCurvature, &f.forceArea, &f.forceVolume, &f.forceThickness, &f.forceTilt,
                            &f.forceRegularization, &f.forceHarmonicBond, &f.forceTotal})
    {
        expect_column3(*m, 0, 0, 0);
    }
    EXPECT_DOUBLE_EQ(f.get_total_force_magnitude(), 0.0);
}

TEST(CoreCovForce, CopyIsDeep)
{
    Force original = numbered_force();
    Force copy(original);
    expect_column3(copy.forceArea, 10, 20, 30);
    expect_column3(copy.forceTotal, 7, 8, 9);

    copy.forceArea.set(0, 0, -99);
    copy.forceTotal.set_all(0);
    expect_column3(original.forceArea, 10, 20, 30);
    expect_column3(original.forceTotal, 7, 8, 9);
}

TEST(CoreCovForce, TotalSumsComponentsWithOptionalRegularization)
{
    Force f = numbered_force();
    // 1 + 10 + 100 + 0.5 - 1 + 1000 + 4, etc. (all exactly representable)
    Matrix &total = f.calculate_total_force();
    EXPECT_EQ(&total, &f.forceTotal);
    expect_column3(f.forceTotal, 1114.5, 2225.25, 3336.125);

    f.calculate_total_force(false);
    expect_column3(f.forceTotal, 114.5, 225.25, 336.125);

    f.forceTotal = column3(3, 0, -4);
    EXPECT_DOUBLE_EQ(f.get_total_force_magnitude(), 5.0);
}

TEST(CoreCovForce, SetAllZeroClearsEveryComponentInPlace)
{
    Force f = numbered_force();
    f.set_all_zero();
    for (const Matrix *m : {&f.forceCurvature, &f.forceArea, &f.forceVolume, &f.forceThickness, &f.forceTilt,
                            &f.forceRegularization, &f.forceHarmonicBond, &f.forceTotal})
    {
        expect_column3(*m, 0, 0, 0);
    }
}

TEST(CoreCovForce, ToStringFullListsEachComponent)
{
    Force f = numbered_force();
    const std::string expected =
        "force_curve = 1 , 2 , 3 , \n"
        "force_area = 10 , 20 , 30 , \n"
        "force_vol = 100 , 200 , 300 , \n"
        "force_thic = 0.5 , 0.25 , 0.125 , \n"
        "force_tilt = -1 , -2 , -3 , \n"
        "force_reg = 1000 , 2000 , 3000 , \n"
        "force_total = 7 , 8 , 9 , \n";
    EXPECT_EQ(f.to_string_full(), expected);
}

TEST(CoreCovForce, ToStringFullAtVertexTagsTheIndex)
{
    Force f = numbered_force();
    const std::string expected =
        "force_curve @ 7 = 1 , 2 , 3 , \n"
        "force_area @ 7 = 10 , 20 , 30 , \n"
        "force_vol @ 7 = 100 , 200 , 300 , \n"
        "force_thic @ 7 = 0.5 , 0.25 , 0.125 , \n"
        "force_tilt @ 7 = -1 , -2 , -3 , \n"
        "force_reg @ 7 = 1000 , 2000 , 3000 , \n"
        "force_total @ 7 = 7 , 8 , 9 , \n";
    EXPECT_EQ(f.to_string_full_at_vertex(7), expected);
}

TEST(CoreCovForce, StreamOperatorPrintsOnlyTheTotal)
{
    Force f = numbered_force();
    std::ostringstream os;
    os << f;
    EXPECT_EQ(os.str(), "7, 8, 9");
}

TEST(CoreCovForce, PlusAddsComponentwiseAndLeavesOperandsAlone)
{
    Force a = numbered_force();
    Force b;
    b.forceCurvature = column3(1, 1, 1);
    b.forceArea = column3(2, 2, 2);
    b.forceVolume = column3(3, 3, 3);
    b.forceThickness = column3(4, 4, 4);
    b.forceTilt = column3(5, 5, 5);
    b.forceRegularization = column3(6, 6, 6);
    b.forceHarmonicBond = column3(7, 7, 7);
    b.forceTotal = column3(8, 8, 8);

    Force sum = a + b;
    expect_column3(sum.forceCurvature, 2, 3, 4);
    expect_column3(sum.forceArea, 12, 22, 32);
    expect_column3(sum.forceVolume, 103, 203, 303);
    expect_column3(sum.forceThickness, 4.5, 4.25, 4.125);
    expect_column3(sum.forceTilt, 4, 3, 2);
    expect_column3(sum.forceRegularization, 1006, 2006, 3006);
    expect_column3(sum.forceHarmonicBond, 11, 12, 13);
    expect_column3(sum.forceTotal, 15, 16, 17);

    // Neither operand changed.
    expect_column3(a.forceCurvature, 1, 2, 3);
    expect_column3(a.forceTotal, 7, 8, 9);
    expect_column3(b.forceCurvature, 1, 1, 1);

    // += agrees with + and returns its left operand.
    Force &ref = (a += b);
    EXPECT_EQ(&ref, &a);
    for (auto member : {&Force::forceCurvature, &Force::forceArea, &Force::forceVolume, &Force::forceThickness,
                        &Force::forceTilt, &Force::forceRegularization, &Force::forceHarmonicBond,
                        &Force::forceTotal})
    {
        for (int i = 0; i < 3; ++i)
        {
            EXPECT_DOUBLE_EQ((a.*member)(i, 0), (sum.*member)(i, 0));
        }
    }
    expect_column3(b.forceTotal, 8, 8, 8);
}

/* ------------------------------------------------------------------------- */
/* Matrix operators (src/linear_algebra/Linear_algebra.cpp)                  */
/* ------------------------------------------------------------------------- */

TEST(CoreCovMatrix, SelfAssignmentKeepsTheData)
{
    Matrix m(std::vector<std::vector<double>>{{1, 2}, {3, 4}});
    gsl_matrix *storage = m.mat;
    Matrix &alias = m;
    Matrix &ret = (m = alias);
    EXPECT_EQ(&ret, &m);
    EXPECT_EQ(m.mat, storage); // not freed and reallocated
    EXPECT_DOUBLE_EQ(m(0, 1), 2);
    EXPECT_DOUBLE_EQ(m(1, 0), 3);
}

TEST(CoreCovMatrix, TransposeSquareInPlace)
{
    Matrix m(std::vector<std::vector<double>>{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    m.transpose();
    const double expected[3][3] = {{1, 4, 7}, {2, 5, 8}, {3, 6, 9}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(m(i, j), expected[i][j]);
}

TEST(CoreCovMatrix, GetColAndGetTransposed)
{
    Matrix m(std::vector<std::vector<double>>{{1, 2, 3}, {4, 5, 6}});

    Matrix col = m.get_col(1);
    ASSERT_EQ(col.nrow(), 2);
    ASSERT_EQ(col.ncol(), 1);
    EXPECT_DOUBLE_EQ(col(0, 0), 2);
    EXPECT_DOUBLE_EQ(col(1, 0), 5);

    Matrix t = m.get_transposed();
    ASSERT_EQ(t.nrow(), 3);
    ASSERT_EQ(t.ncol(), 2);
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(t(j, i), m(i, j));
    // The source is untouched.
    EXPECT_EQ(m.nrow(), 2);
    EXPECT_DOUBLE_EQ(m(0, 2), 3);
}

TEST(CoreCovMatrix, StreamOperatorFormatsRowsAndColumns)
{
    Matrix grid(std::vector<std::vector<double>>{{1, 2}, {3.5, -4}});
    std::ostringstream os;
    os << grid;
    EXPECT_EQ(os.str(), "1, 2, \n3.5, -4, \n");

    Matrix column(std::vector<std::vector<double>>{{1}, {2}, {3}});
    std::ostringstream oc;
    oc << column;
    EXPECT_EQ(oc.str(), "1, 2, 3, "); // a column prints on one line

    Matrix row(std::vector<std::vector<double>>{{5, 6, 7}});
    std::ostringstream orow;
    orow << row;
    EXPECT_EQ(orow.str(), "5, 6, 7, \n");
}

TEST(CoreCovMatrix, DotIsRowTimesTransposedMatrix)
{
    // dot(m1, m2) = m1 * m2^T into a (1, m1.ncol()) result, which is shaped
    // right when m1 is (1,N) and m2 is (N,N).
    Matrix row(std::vector<std::vector<double>>{{1, 2}});
    Matrix square(std::vector<std::vector<double>>{{3, 4}, {5, 6}});
    Matrix r = dot(row, square);
    ASSERT_EQ(r.nrow(), 1);
    ASSERT_EQ(r.ncol(), 2);
    EXPECT_DOUBLE_EQ(r(0, 0), 1 * 3 + 2 * 4);
    EXPECT_DOUBLE_EQ(r(0, 1), 1 * 5 + 2 * 6);

    Matrix a(std::vector<std::vector<double>>{{2.5}});
    Matrix b(std::vector<std::vector<double>>{{-4}});
    Matrix s = dot(a, b);
    ASSERT_EQ(s.nrow(), 1);
    ASSERT_EQ(s.ncol(), 1);
    EXPECT_DOUBLE_EQ(s(0, 0), -10);
}

TEST(CoreCovMatrix, ScalarDivision)
{
    Matrix m(std::vector<std::vector<double>>{{2, 4}, {-6, 8}});
    Matrix half = m / 2.0;
    EXPECT_DOUBLE_EQ(half(0, 0), 1);
    EXPECT_DOUBLE_EQ(half(0, 1), 2);
    EXPECT_DOUBLE_EQ(half(1, 0), -3);
    EXPECT_DOUBLE_EQ(half(1, 1), 4);
    EXPECT_DOUBLE_EQ(m(1, 0), -6); // operator/ takes its operand by value

    Matrix &ref = (m /= 4.0);
    EXPECT_EQ(&ref, &m);
    EXPECT_DOUBLE_EQ(m(0, 0), 0.5);
    EXPECT_DOUBLE_EQ(m(0, 1), 1);
    EXPECT_DOUBLE_EQ(m(1, 0), -1.5);
    EXPECT_DOUBLE_EQ(m(1, 1), 2);
}

TEST(CoreCovMatrix, CrossColIsTheRightHandedCrossProduct)
{
    expect_column3(cross_col(column3(1, 0, 0), column3(0, 1, 0)), 0, 0, 1);
    expect_column3(cross_col(column3(0, 1, 0), column3(0, 0, 1)), 1, 0, 0);
    expect_column3(cross_col(column3(1, 2, 3), column3(4, 5, 6)), -3, 6, -3);
    expect_column3(cross_col(column3(4, 5, 6), column3(1, 2, 3)), 3, -6, 3);

    // Orthogonal to both inputs, and zero for parallel inputs.
    Matrix a = column3(1.5, -2, 0.25);
    Matrix b = column3(-0.5, 3, 2);
    Matrix c = cross_col(a, b);
    EXPECT_NEAR(dot_col(c, a), 0, 1e-14);
    EXPECT_NEAR(dot_col(c, b), 0, 1e-14);
    expect_column3(cross_col(a, a * 2.0), 0, 0, 0);

    // Agrees with the row-vector version.
    Matrix cr = cross_row(a.get_transposed(), b.get_transposed());
    for (int i = 0; i < 3; ++i)
        EXPECT_DOUBLE_EQ(cr(0, i), c(i, 0));
}

TEST(CoreCovMatrix, NegativeFlipsEverySign)
{
    Matrix m(std::vector<std::vector<double>>{{1, -2, 0}, {3.5, 0, -4}});
    Matrix neg(2, 3);
    negative(m, neg);
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(neg(i, j), -m(i, j));
    EXPECT_DOUBLE_EQ(m(0, 1), -2); // input untouched
}

TEST(CoreCovMatrix, UnitVectorOfZeroIsZero)
{
    Matrix zero = column3(0, 0, 0);
    get_unit_vector(zero);
    expect_column3(zero, 0, 0, 0);
    EXPECT_FALSE(std::isnan(zero(0, 0)));

    Matrix out = column3(9, 9, 9); // stale contents must be overwritten
    get_unit_vector(column3(0, 0, 0), out);
    expect_column3(out, 0, 0, 0);

    get_unit_vector(column3(0, 3, 4), out);
    expect_column3(out, 0, 0.6, 0.8);
}

TEST(CoreCovMatrix, KronOfTwoColumnsIsTheOuterProduct)
{
    Matrix v1 = Matrix(std::vector<std::vector<double>>{{1}, {-2}});
    Matrix v2 = column3(3, 4, 0.5);
    Matrix k = kron(v1, v2);
    ASSERT_EQ(k.nrow(), 2);
    ASSERT_EQ(k.ncol(), 3);
    const double expected[2][3] = {{3, 4, 0.5}, {-6, -8, -1}};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(k(i, j), expected[i][j]);

    // It matches v1 * v2^T.
    Matrix outer = v1 * v2.get_transposed();
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(k(i, j), outer(i, j));
}
