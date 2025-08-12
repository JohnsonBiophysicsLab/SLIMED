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