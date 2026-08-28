#include "test_linear_algebra.hpp"

/**
 * @brief Test Matrix Initialization.
 * 
 * Test cases include initializing a column vector by size, initializing a square matrix by size,
 * initializing a matrix by a 2D std::vector, and initializing a matrix by another matrix (copy).
 */
TEST(LinearAlgebraTest, MatrixInitialization){
    // test initialize column vector by size
    Matrix mColVec(3);
    EXPECT_EQ(mColVec.nrow(), 3);
    // test initialize matrix by size
    Matrix mSquare(4, 4);
    // test initialize matrix by 2D std::vector
    std::vector<std::vector<double>> vecValue{
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12}
    };
    Matrix mFromVector(vecValue);
    // test initialize matrix by matrix (copy)
    Matrix mFromMatrix(mFromVector);
    EXPECT_DOUBLE_EQ(mFromVector.get(2,0), 9.0);
}

// Test for the get_unit_vector function with an input matrix
TEST(GetUnitVectorTest, InputMatrixTest)
{
    // Arrange
    Matrix inputMatrix(3, 1);  // Create a 3x1 matrix
    inputMatrix.set(0, 0, 32.0);
    inputMatrix.set(1, 0, 75.0);
    inputMatrix.set(2, 0, -29.0);

    // Act
    get_unit_vector(inputMatrix);

    // Assert
    double magnitude = inputMatrix.calculate_norm();
    EXPECT_DOUBLE_EQ(magnitude, 1.0); // The magnitude of a unit vector should be 1
    EXPECT_DOUBLE_EQ(inputMatrix.get(0, 0), 0.36975075500635607);
    EXPECT_DOUBLE_EQ(inputMatrix.get(1, 0), 0.86660333204614703);
    EXPECT_DOUBLE_EQ(inputMatrix.get(2, 0), -0.3350866217245102);

    // Add additional assertions if needed
}

// Test for the get_unit_vector function with an input and output matrix
TEST(GetUnitVectorTest, InputOutputMatrixTest)
{
    // Arrange
    Matrix inputMatrix(3, 1);  // Create a 3x1 matrix
    Matrix outputMatrix(3, 1); // Create another 3x1 matrix for the output
    // Initialize your matrices with values here if needed
    inputMatrix.set(0, 0, 32.0);
    inputMatrix.set(1, 0, 75.0);
    inputMatrix.set(2, 0, -29.0);

    // Act
    get_unit_vector(inputMatrix, outputMatrix);

    // Assert
    double inputMagnitude = inputMatrix.calculate_norm();
    double outputMagnitude = outputMatrix.calculate_norm();
    EXPECT_DOUBLE_EQ(inputMagnitude, 86.544786093675228); // The magnitude of input should remain unchanged
    EXPECT_DOUBLE_EQ(outputMagnitude, 1.0); // The magnitude of output should be 1
    EXPECT_DOUBLE_EQ(outputMatrix.get(0, 0), 0.36975075500635607);
    EXPECT_DOUBLE_EQ(outputMatrix.get(1, 0), 0.86660333204614703);
    EXPECT_DOUBLE_EQ(outputMatrix.get(2, 0), -0.3350866217245102);
}

/**
 * @brief Test Matrix Addition.
 * 
 * Test cases include adding a 4x4 matrix to itself, checking operator+ override,
 * checking operator+= override, and adding two 2x3 matrices using the addition() function.
 */
TEST(LinearAlgebraTest, MatrixAddition){
    // Sample matrix 1 : 4 x 4 matrix add itself
    std::vector<std::vector<double>> v2d = {{10, -1, 2, 0},
                                            {-1, 11, -1, 3},
                                            {2, -1, 10, -1},
                                            {0, 3, -1, 8}};
    Matrix mat(v2d);
    // check operator+ overide
    Matrix matSum = mat + mat;
    // Reference answer
    std::vector<std::vector<double>> v2dSumRef = {{20, -2, 4, 0},
                                                {-2, 22, -2, 6},
                                                {4, -2, 20, -2},
                                                {0, 6, -2, 16}};
    Matrix matSumRef(v2dSumRef);
    matSum -= matSumRef; // get difference
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            EXPECT_DOUBLE_EQ(matSum.get(i,j), 0.0);
        }
    }
    
    // check operator+= overide
    mat += mat;
    mat -= matSumRef; // get difference
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            EXPECT_DOUBLE_EQ(mat.get(i,j), 0.0);
        }
    }
    
    // Sample matrix 2: 2 x 3 matrix add another 2 x 3 matrix
    // check addition()
    std::vector<std::vector<double>> v2x3a = {{4.35, -1.24, 2.19},
                                            {-91.3, 194.6, 0.0}};
    std::vector<std::vector<double>> v2x3b = {{0.0, -1.24, 1898.53},
                                            {-35.5, -194.6, -0.001}};
    Matrix matA(v2x3a); 
    Matrix matB(v2x3b); 
    Matrix matSum2(2,3); 
    addition(v2x3a, v2x3b, matSum2); 
    // get reference answer
    std::vector<std::vector<double>> v2x3SumRef = {{4.35, -2.48, 1900.72},
                                                {-126.8, 0.0, -0.001}};
    Matrix matSumRef2(v2x3SumRef);
    matSum2 -= matSumRef2;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            EXPECT_DOUBLE_EQ(matSum2.get(i,j), 0.0);
        }
    }
}



/**
 * @brief Test Matrix Inverse.
 * 
 * Test cases include testing the inverse of an identity matrix and a 4x4 matrix.
 */
TEST(LinearAlgebraTest, MatrixInverse){
    // Sample matrix 1 : identity matrix
    Matrix id(15, 15);
    id.set_identity();

    Matrix idInv;
    id.get_inverted(idInv);
    
    // idInv should equal to id
    idInv -= id; // diffVec
    for (int i = 0; i < 15; i++)
    {
        for (int j = 0; j < 15; j++)
        {
            EXPECT_DOUBLE_EQ(idInv.get(i,j), 0.0);
        }
    }

    // Sample matrix 2 : 4 x 4 matrix
    std::vector<std::vector<double>> v2d = {{10, -1, 2, 0},
                                            {-1, 11, -1, 3},
                                            {2, -1, 10, -1},
                                            {0, 3, -1, 8}};
    Matrix mat(v2d);
    // Inverted matrix
    Matrix matInv;
    mat.get_inverted(matInv);
    // Inverse reference result calculated by wolfram alpha
    // https://www.wolframalpha.com/input?i=Inverse+of+%7B%7B10%2C+-1%2C+2%2C+0%7D%2C+%7B-1%2C+11%2C
    // +-1%2C+3%7D%2C+%7B2%2C+-1%2C+10%2C+-1%7D%2C+%7B0%2C+3%2C+-1%2C+8%7D%7D
    std::vector<std::vector<double>> v2dInvRef = {{259.0 / 2465, 23.0 / 2465, -3.0 / 145, -3.0 / 493},
                                                {23.0 / 2465, 758.0 / 7395, 2.0 / 435, -56.0 / 1479},
                                                {-3.0 / 145, 2.0 / 435, 46.0 / 435, 1.0 / 87},
                                                {-3.0 / 493, -56.0 / 1479, 1.0 / 87, 208.0 / 1479}};
    Matrix matInvRef(v2dInvRef);
    matInv -= matInvRef; // diffVec
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            EXPECT_DOUBLE_EQ(idInv.get(i,j), 0.0);
        }
    }
}
/**
 * @brief dot_row must sum over every column of a row vector.
 *
 * The loop bound used to be m1.size1 -- the row count -- while the indices
 * walked the columns, so a 1 x 3 row vector returned m1(0,0)*m2(0,0) alone.
 * That is the volume-functional bug: Mesh::enumerate_gauss_quadrature_point_area_volume()
 * is the sole caller and it passes 1 x 3 vectors, which made param.vol come
 * out as V/3 on a closed surface.
 */
TEST(DotProductTest, DotRowSumsAllColumnsOfRowVector)
{
    Matrix m1(1, 3);
    m1.set(0, 0, 1.0);
    m1.set(0, 1, 2.0);
    m1.set(0, 2, 3.0);

    Matrix m2(1, 3);
    m2.set(0, 0, 4.0);
    m2.set(0, 1, -5.0);
    m2.set(0, 2, 6.0);

    // 1*4 + 2*(-5) + 3*6 == 12, not the x-only product 1*4 == 4.
    EXPECT_DOUBLE_EQ(dot_row(m1, m2), 12.0);
}

/// dot_row must not be limited to three columns either.
TEST(DotProductTest, DotRowHandlesWiderRowVectors)
{
    Matrix m1(1, 5);
    Matrix m2(1, 5);
    for (int j = 0; j < 5; j++)
    {
        m1.set(0, j, static_cast<double>(j + 1)); // 1 2 3 4 5
        m2.set(0, j, 2.0);
    }

    // 2 * (1 + 2 + 3 + 4 + 5) == 30
    EXPECT_DOUBLE_EQ(dot_row(m1, m2), 30.0);
}

/**
 * @brief Orthogonal vectors give zero, and the cancellation needs every term.
 *
 * The x-components alone multiply to +1 here, so this only reaches zero if the
 * y-components are summed in as well.
 */
TEST(DotProductTest, DotProductOfOrthogonalVectorsIsZero)
{
    const double u[3] = {1.0, 1.0, 0.0};
    const double v[3] = {1.0, -1.0, 0.0};

    Matrix r1(1, 3);
    Matrix r2(1, 3);
    Matrix c1(3, 1);
    Matrix c2(3, 1);
    for (int k = 0; k < 3; k++)
    {
        r1.set(0, k, u[k]);
        r2.set(0, k, v[k]);
        c1.set(k, 0, u[k]);
        c2.set(k, 0, v[k]);
    }

    EXPECT_DOUBLE_EQ(dot_row(r1, r2), 0.0);
    EXPECT_DOUBLE_EQ(dot_col(c1, c2), 0.0);
}

/// dot_col sums over every row of a column vector -- the mirror of the above.
TEST(DotProductTest, DotColSumsAllRowsOfColumnVector)
{
    Matrix m1(3, 1);
    m1.set(0, 0, 1.0);
    m1.set(1, 0, 2.0);
    m1.set(2, 0, 3.0);

    Matrix m2(3, 1);
    m2.set(0, 0, 4.0);
    m2.set(1, 0, -5.0);
    m2.set(2, 0, 6.0);

    EXPECT_DOUBLE_EQ(dot_col(m1, m2), 12.0);
}

/// dot_col must not be limited to three rows either.
TEST(DotProductTest, DotColHandlesTallerColumnVectors)
{
    Matrix m1(5, 1);
    Matrix m2(5, 1);
    for (int i = 0; i < 5; i++)
    {
        m1.set(i, 0, static_cast<double>(i + 1)); // 1 2 3 4 5
        m2.set(i, 0, 2.0);
    }

    EXPECT_DOUBLE_EQ(dot_col(m1, m2), 30.0);
}

/**
 * @brief dot_row and dot_col agree on the same vector, transposed.
 *
 * Both helpers are non-square by construction at their call sites, so this
 * pins the two bounds against each other rather than against a literal.
 */
TEST(DotProductTest, DotRowAndDotColAgreeOnTransposedVectors)
{
    const double u[3] = {0.5, -1.25, 3.75};
    const double v[3] = {-2.0, 0.25, 1.5};

    Matrix r1(1, 3);
    Matrix r2(1, 3);
    Matrix c1(3, 1);
    Matrix c2(3, 1);
    for (int k = 0; k < 3; k++)
    {
        r1.set(0, k, u[k]);
        r2.set(0, k, v[k]);
        c1.set(k, 0, u[k]);
        c2.set(k, 0, v[k]);
    }

    EXPECT_DOUBLE_EQ(dot_row(r1, r2), dot_col(c1, c2));
    EXPECT_DOUBLE_EQ(dot_row(r1, r2), u[0] * v[0] + u[1] * v[1] + u[2] * v[2]);
}
