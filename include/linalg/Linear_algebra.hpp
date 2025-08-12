/**
 * @file Linear_algebra.hpp
 * This file serves as an interface for implementing any
 * linear algebra package for the continuum membrane model.
 *
 * The code here defines a Matrix class that uses the GNU Scientific
 * Library (GSL) to implement various linear algebra operations.
 * The class has four constructors:
 * 
 * @warning this interface between Linalg package (currently GSL linalg) mutes all
 * pre-operation checks to enhance the computational efficiency as much as possible.
 * However, this may also results in "core-dumped" or "double-free" error if the class
 * functions are not handled correctly. Please make sure to read thru the docstring,
 * especially the prerequistes and warnings carefully before you reference to the functions!
 * 
 *
 * @todo benchmark - In this case, if my matrices are small (< 12 * 3), but I need to repeat these matrices
 * calculation for a large number of times, is gsl_blas_dgemm still more eifficent?
 *
 * For small matrices like the ones you've described, it's possible that an element-wise implementation could be faster than using
 * gsl_blas_dgemm() due to the overhead associated with setting up the function call and memory management.
 * However, this depends on the context in which the function is being used and the specific hardware being used.
 *
 * If the calculation of these small matrices is being repeated a large number of times, then the overhead of
 * setting up the function call and memory management might be amortized over multiple calculations, making
 * gsl_blas_dgemm() more efficient.
 *
 * Ultimately, the best approach would be to benchmark both implementations (element-wise and using
 * gsl_blas_dgemm() on your specific system and compare their performance for your use case.
 *
 * -ChatGPT
 *
 * @date 2023-03-20 (Created)
 * @date 2024-01-29 (Edited)
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <vector>
#include <iostream> //for bit push << overriding

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/**
 * @brief Matrix used in continuum membrane model
 * 
 * This class represents a N x M matrix. This header is an interface
 * that can be implemented differently to support the use of different
 * linear algebra pacakages. Currently gsl_blas is used for matrix calculation.
 * Note that it is important to double check the destructor ~Matrix() and assign
 * operator (operator=) if you would like to customize an implementation. Be
 * careful of double free and memory leak that might be caused by misuse of
 * pointers.
 * - Moon Ying
 * 
 * @note To keep all the matrices definition
 * consistent, all the coordinates and forces will be defined as 3 x 1 column vectors.
 * 
 */
class Matrix
{
public:
    gsl_matrix *mat;

    /**
     * @brief Construct a new Matrix object and initialize the members
     * to zeros. The setToZero variable is a place holder to differentiate
     * this function from a malloc or alloc function.
     *
     * @param nrow
     * @param ncol
     * @param setToZero
     */
    Matrix(const int &nrow, const int &ncol, const bool &setToZero);



    /**
     * @brief Default constructor
     * 
     */
    Matrix();


    /**
     * @brief Constructor for creating a Matrix with a specified number of rows and one column
     *
     * This constructor creates a new Matrix object with the specified number of rows and one column.
     * The elements of the matrix are uninitialized and may contain garbage values.
     *
     * @param nrow The number of rows in the matrix
     */
    Matrix(const int &nrow);

    /**
     * @brief Constructor for creating a Matrix with a specified number of rows and columns
     *
     * This constructor creates a new Matrix object with the specified number of rows and columns.
     * The elements of the matrix are uninitialized and may contain garbage values.
     *
     * @param nrow The number of rows in the matrix
     * @param ncol The number of columns in the matrix
     */
    Matrix(const int &nrow, const int &ncol);

    /**
     * @brief Copy constructor for creating a deep copy of a Matrix object
     *
     * This constructor creates a new Matrix object that is a deep copy of the specified Matrix object.
     * All elements of the original matrix are copied to the new matrix, so the two matrices are completely independent.
     *
     * @param mat_src The Matrix object to be copied
     */
    Matrix(const Matrix &mat_src);

    /**
     * @brief Intialize a matrix with a 2d vector
     *
     */
    Matrix(const std::vector<std::vector<double>> &vec);

    /**
     * @brief Assignment operator
     * 
     * @param mat_src 
     * @return Matrix& 
     */
    Matrix& operator=(const Matrix &mat_src);

    /**
     * @brief Destructor for the Matrix class that frees the memory allocated by the 
     * mat member.
     * 
     * A destructor is a special member function of a C++ class that is called automatically
     * when an object of the class goes out of scope or is otherwise destroyed. The purpose
     * of the destructor is to release any resources or memory that the object was holding onto,
     * to ensure that there are no memory leaks or other resource leaks caused by the object's
     * destruction.
     */
    ~Matrix();

    /**
     * @brief Free memory allocation for a Matrix object 
     * 
     */
    void free();

    /**
     * @brief get number of rows
     *
     */
    int nrow() const;

    /**
     * @brief get number of columns
     *
     */
    int ncol() const;

    /**
     * @brief Get element value at i,j
     *
     */
    double get(const int &i, const int &j) const;

    /**
     * @brief Get element value at i,j
     *
     */
    double operator()(const int &i, const int &j) const;

    /**
     * @brief Set element value at i,j to v
     *
     */
    void set(const int &i, const int &j, const double &v);

    /**
     * @brief Set all elements to v
     *
     */
    void set_all(const double &v);

    /**
     * @brief Transpose the matrix. This overwrites the original
     * matrix.
     * @note Use if the source matrix
     * is needed.
     *
     */
    void transpose();

    /**
     * @brief Set the matrix to identity matrix.
     * @note This function checks whether the matrix is a square
     * matrix and will throw error if the matrix is not square
     *
     */
    void set_identity();

    /**
     * @brief Get the row vector in the format of Matrix(1, N)
     *
     * @param irow
     * @return Matrix
     */
    Matrix get_row(const int &irow) const;

    /**
     * @brief Get the col vector in the format of Matrix(N, 1)
     *
     * @param icol
     * @return Matrix
     */
    Matrix get_col(const int &icol) const;

    /**
     * @brief Set the row vector from col vector
     *
     * @param irow
     * @return Matrix
     */
    void set_row_from_col(const int &destRow, const Matrix &srcMat, const int &srcCol);

    /**
     * @brief Get the tranposed matrix
     *
     * @return Matrix
     */
    Matrix get_transposed() const;

    /**
     * @brief This function uses cholesky decomposition to
     * solve for the inverse of the input square matrix.
     * @ref https://www.gnu.org/software/gsl/doc/html/linalg.html#c.gsl_linalg_cholesky_decomp1
     *
     * @param matrix
     * @return Matrix
     */
    void get_inverted(Matrix &matInverted);

    /**
     *
     *  @brief Calculates the Euclidean norm of the member matrix mat.
     */
    double calculate_norm() const;

    /**
     *
     *  @brief Calculates the Euclidean quadrance of the member matrix mat.
     */
    double calculate_quad() const;

};

/**
 * @brief This function allocates a new
 * Matrix
 * object with the specified number of rows and columns, and initializes all elements to zero.
 *
 * @param nrow
 * @param ncol
 * @return Matrix
 */
Matrix mat_calloc(const int &nrow, const int &ncol);

/**
 * @brief Bit push all element of a matrix in the format of:
 *
 * a11, a12, a13, ...,
 * a21, a22, a23, ...,
 * ...,
 *
 * @param stream ostream
 * @param matrix matrix to be output as string
 * @return std::ostream&
 */
std::ostream &operator<<(std::ostream &stream, const Matrix &matrix);

/**
 * @brief This implementation takes two gsl_matrix pointers as
    arguments and adds the elements of the second matrix to the first matrix.
    It first checks that the matrices have the same dimensions, and throws an
    exception if they do not. Then it iterates over the elements of the matrices
    and adds the corresponding elements. Finally, it returns a pointer to the
    modified first matrix.
    Note that this implementation modifies the first matrix in place, rather
    than creating a new matrix. If you want to create a new matrix instead, you
    could modify the implementation accordingly.
 */
Matrix &operator+=(Matrix &m1, const Matrix &m2);

/**

    @brief This implementation checks that the matrices have the
    same dimensions and allocates memory for the result of the
    addition. It then adds the matrices element-wise and stores
    the result in a newly allocated matrix result, which is then
    returned.
*/
Matrix operator+(const Matrix &m1, const Matrix &m2);

/**
 * @brief This implementation takes two gsl_matrix pointers as
    arguments and subtract the elements of the second matrix from the first matrix.
    It first checks that the matrices have the same dimensions, and throws an
    exception if they do not. Then it iterates over the elements of the matrices
    and adds the corresponding elements. Finally, it returns a pointer to the
    modified first matrix.
    Note that this implementation modifies the first matrix in place, rather
    than creating a new matrix. If you want to create a new matrix instead, you
    could modify the implementation accordingly.
 */
Matrix &operator-=(Matrix &m1, const Matrix &m2);

/**

    @brief This implementation checks that the matrices have the
    same dimensions and allocates memory for the result of the
    subtraction. It then adds the matrices element-wise and stores
    the result in a newly allocated matrix result, which is then
    returned.
*/
Matrix operator-(const Matrix &m1, const Matrix &m2);

/**
 * @brief This implementation uses the gsl_matrix_scale function to scale
 * each element in the matrix by the scalar value. The function returns the
 * modified matrix, which allows for method chaining when using the *= operator.
 *
 * @param mat
 * @param scalar
 * @return Matrix
 */
Matrix &operator*=(Matrix &mat, const double &scalar);

/**
 * @brief This operator takes a double scalar and a matrix as its operands. It returns a new
 * gsl_matrix* that is the result of scaling the input matrix by the scalar.
 *
 * @param scalar
 * @param matrix
 * @return Matrix
 */
///@{
Matrix operator*(const Matrix &matrix, const double &scalar);
Matrix operator*(const double &scalar, const Matrix &matrix);
///@}

/**
 * @brief This implementation checks that the matrices have compatible
    dimensions for multiplication and uses the BLAS library to perform
    the matrix multiplication. The result is stored in a temporary matrix,
    which is then copied back into m1. Finally, the temporary memory is freed
    before returning m1.
 *
 * @param m1
 * @param m2
 * @return Matrix
 */
Matrix &operator*=(Matrix m1, const Matrix &m2);

/**
 * @brief This implementation checks that the matrices have compatible
    dimensions for multiplication and uses the BLAS library to perform
    the matrix multiplication. The result is stored in a newly allocated matrix
    result, which is then returned. The temporary memory is not freed here,
    and it is the responsibility of the caller to free the memory once the
    result is no longer needed.
 *
 * @param m1
 * @param m2
 * @return Matrix
 */
///@{
Matrix operator*(const Matrix &m1, const Matrix &m2);
Matrix dot(const Matrix &m1, const Matrix &m2);
///@}

/**
 * @brief This implementation uses the gsl_matrix_scale function to scale
 * each element in the matrix by 1/the scalar value. The function returns the
 * modified matrix, which allows for method chaining when using the *= operator.
 *
 * @param mat
 * @param scalar
 * @return Matrix
 */
Matrix &operator/=(Matrix &mat, const double &scalar);

/**
 * @brief This operator takes a double scalar and a matrix as its operands. It returns a new
 * gsl_matrix* that is the result of scaling the input matrix by 1/scalar.
 *
 * @param scalar
 * @param matrix
 * @return Matrix
 */
Matrix operator/(Matrix matrix, const double &scalar);

/**
 * @brief This implementation takes two column vectors as
 * arguments and returns their dot product as double.
 *
 * @note Note that this implementation assumes that the input matrices are
 * already column vectors, so no additional transposition is necessary.
 */
double dot_col(const Matrix &m1, const Matrix &m2);

/**
 * @brief This implementation takes two row vectors as
 * arguments and returns their dot product as double.
 *
 * @note Note that this implementation converts the row vectors into column
 * vectors before computing the dot product.
 */
double dot_row(const Matrix &m1, const Matrix &m2);

/**
 * @brief This implementation takes two column vectors as
 * arguments and returns their cross product as a new column vector.
 *
 * @note Note that this implementation assumes that the input matrices are
 * already column vectors, so no additional transposition is necessary.
 */
Matrix cross_col(const Matrix &m1, const Matrix &m2);

/**
 * @brief This implementation takes two row vectors  as
 * arguments and returns their cross product as a new row vector.
 *
 *
 * @note Note that this implementation assumes that the input matrices are
 * already row vectors, so no additional transposition is necessary.
 */
Matrix cross_row(const Matrix &m1, const Matrix &m2);

/**
 * @brief Computes the negative of a given matrix.
 *
 * This function computes the element-wise negative of a given matrix and
 * stores the result in another matrix. Element-wise negative is defined as
 * multiplying each element of the matrix with -1. The input matrix 'm1' is
 * not modified by this function.
 *
 * @param m1 The matrix to compute the negative of.
 * @param m_neg The matrix to store the result in.
 *
 * @pre Both matrices 'm1' and 'm_neg' must have been allocated memory and initialized properly.
 * @pre Both matrices 'm1' and 'm_neg' must have the same dimensions.
 *
 * @post The matrix 'm_neg' contains the element-wise negative of 'm1'.
 */
void negative(const Matrix &m1, Matrix &m_neg);

/**
 * @brief Get the unit vector in the direction of m which is assumed to
 * be a (3, 1) Matrix and overwrites the result to m
 *
 * @param m1 The input vector to compute the unit vector of 
 * and the output matrix to store the resulting unit vector.
 *
 * @pre The input matrix 'm' must be a 3x1 matrix.
 *
 * @post The output matrix 'm' contains the unit vector in the direction of 'm'.
 */
void get_unit_vector(Matrix &m);

/**
 * @brief Get the unit vector in the direction of m1 which is a (3, 1) Matrix
 *
 * This function computes the unit vector in the direction of a given input vector,
 * which is assumed to be a 3x1 matrix. The result is stored in another 3x1 matrix.
 *
 * @warning m1 and m_unit cannot be the same Matrix instances!
 * Use get_unit_vector(Matrix &m) instead.
 * 
 * @param m1 The input vector to compute the unit vector of.
 * @param m_unit The output matrix to store the resulting unit vector.
 *
 * @pre The input matrix 'm1' must be a 3x1 matrix.
 *
 * @post The output matrix 'm_unit' contains the unit vector in the direction of 'm1'.
 */
void get_unit_vector(const Matrix &m1, Matrix &m_unit);

/**
 * @brief Compute the cross product of two 3D vectors represented as 3x1 matrices.
 *
 * This function computes the cross product of two input vectors, which are assumed to
 * be 3x1 matrices. The result is stored in another 3x1 matrix passed in as an argument.
 *
 * @param m1 The first input vector as a 3x1 matrix.
 * @param m2 The second input vector as a 3x1 matrix.
 * @param temp The output matrix to store the resulting cross product.
 *
 * @pre All input matrices must have 3 rows and 1 column.
 *
 * @post The output matrix 'temp' contains the cross product of 'm1' and 'm2'.
 */
void cross(const Matrix &m1, const Matrix &m2, Matrix &temp);

/**
 * @brief Compute the element-wise addition of two matrices and store the result in a third matrix.
 *
 * This function computes the element-wise addition of two input matrices of equal dimensions
 * and stores the result in a third matrix. All input matrices are assumed to be of type Matrix,
 * which internally contains a gsl_matrix pointer for efficiency.
 *
 * @param m1 The first input matrix.
 * @param m2 The second input matrix.
 * @param tmp The output matrix to store the resulting sum.
 *
 * @pre All input matrices must have the same dimensions.
 *
 * @post The output matrix 'tmp' contains the element-wise sum of 'm1' and 'm2'.
 */
void addition(const Matrix &m1, const Matrix &m2, Matrix &tmp);

/**
 * @brief Compute the element-wise subtraction of two matrices and store the result in a third matrix.
 *
 * This function computes the element-wise subtraction of two input matrices of equal dimensions
 * and stores the result in a third matrix. All input matrices are assumed to be of type Matrix,
 * which internally contains a gsl_matrix pointer for efficiency.
 *
 * @param m1 The first input matrix.
 * @param m2 The second input matrix.
 * @param tmp The output matrix to store the resulting difference.
 *
 * @pre All input matrices must have the same dimensions.
 *
 * @post The output matrix 'tmp' contains the element-wise difference of 'm1' and 'm2'.
 */
void subtraction(const Matrix &m1, const Matrix &m2, Matrix &tmp);

/**
 * @brief Multiply a matrix by a constant and store the result in a second matrix.
 *
 * This function multiplies the elements of an input matrix with a scalar constant and stores
 * the result in a second output matrix. All input matrices are assumed to be of type Matrix
 *
 * @param m1 The input matrix to be multiplied.
 * @param num The scalar constant to multiply with.
 * @param tmp The output matrix to store the resulting product.
 *
 * @pre The dimensions of 'm1' and 'tmp' must match.
 *
 * @post The output matrix 'tmp' contains the element-wise multiplication of 'm1' and 'num'.
 */
void multiplication(const Matrix &m1, const Matrix &m2, Matrix &tmp);

/**
 * @brief Multiply a matrix by a constant and store the result in a second matrix.
 *
 * This function multiplies the elements of an input matrix with a scalar constant and stores
 * the result in a second output matrix. All input matrices are assumed to be of type Matrix
 *
 * @param m1 The input matrix to be multiplied.
 * @param num The scalar constant to multiply with.
 * @param tmp The output matrix to store the resulting product.
 *
 * @pre The dimensions of 'm1' and 'tmp' must match.
 *
 * @post The output matrix 'tmp' contains the element-wise multiplication of 'm1' and 'num'.
 */
void const_multiplication(const Matrix &m1, const double num, Matrix &tmp);

/**
 * @brief Divide a matrix by a constant and store the result in a second matrix.
 *
 * This function divides the elements of an input matrix by a scalar constant and stores
 * the result in a second output matrix. All input matrices are assumed to be of type Matrix.
 *
 * @param m1 The input matrix to be divided.
 * @param num The scalar constant to divide by.
 * @param tmp The output matrix to store the resulting quotient.
 *
 * @pre The dimensions of 'm1' and 'tmp' must match.
 *
 * @post The output matrix 'tmp' contains the element-wise division of 'm1' and 'num'.
 */
void const_division(const Matrix &m1, const double num, Matrix &tmp);

/**
 * @brief Multiply a row vector by a matrix and store the result in a row vector.
 *
 * This function multiplies a row vector with a matrix and stores the resulting product in
 * another row vector. All input matrices are assumed to be of type Matrix, which internally
 * contains a gsl_matrix pointer for efficiency.
 *
 * @param v1 The input row vector to be multiplied.
 * @param m1  The input matrix.
 * @param tmp The output row vector to store the resulting product.
 *
 * @pre The number of columns of 'v1' must match the number of rows of 'm1', and the dimensions
 *      of 'tmp' must match the number of columns of 'm1'.
 *
 * @post The output row vector 'tmp' contains the product of the input row vector 'v1' and the input
 *       matrix 'm1'.
 *
 * @todo benchmark - In this case, if my matrices are small (< 12 * 3), but I need to repeat these matrices
 * calculation for a large number of times, is gsl_blas_dgemm still more eifficent?
 *
 * For small matrices like the ones you've described, it's possible that an element-wise implementation could be faster than using
 * gsl_blas_dgemm() due to the overhead associated with setting up the function call and memory management.
 * However, this depends on the context in which the function is being used and the specific hardware being used.
 *
 * If the calculation of these small matrices is being repeated a large number of times, then the overhead of
 * setting up the function call and memory management might be amortized over multiple calculations, making
 * gsl_blas_dgemm() more efficient.
 *
 * Ultimately, the best approach would be to benchmark both implementations (element-wise and using
 * gsl_blas_dgemm() on your specific system and compare their performance for your use case.
 */
void colvec_matrix_multiplication(const Matrix &v1, const Matrix &m1, Matrix &tmp);

/**
 * @brief the Kronecker product between two column matrices and returns the
 * result in an output matrix. The Kronecker product is defined as follows:
 *
 *     kron(v1, v2) = [v1(1)*v2(1)  v1(1)*v2(2)  ...  v1(1)*v2(n2)
 *                     v1(2)*v2(1)  v1(2)*v2(2)  ...  v1(2)*v2(n2)
 *                     ...          ...          ...  ...
 *                     v1(n1)*v2(1)  v1(n1)*v2(2)  ...  v1(n1)*v2(n2)]
 *
 * where v1 and v2 are column matrices of size n1 and n2 respectively,
 * and the output matrix has size n1*n2 x n1*n2.
 *
 * @param[in] v1        First input column matrix.
 * @param[in] v2        Second input column matrix.
 * @deprecated Might cause memory leak!
 * @return Output matrix to store the result.
 *
 * @throw std::invalid_argument if the output matrix has incompatible dimensions.
 */
Matrix kron(const Matrix& v1, const Matrix& v2);

/**
 * @brief the Kronecker product between two column matrices and returns the
 * result in an output matrix. The Kronecker product is defined as follows:
 *
 *     kron(v1, v2) = [v1(1)*v2(1)  v1(1)*v2(2)  ...  v1(1)*v2(n2)
 *                     v1(2)*v2(1)  v1(2)*v2(2)  ...  v1(2)*v2(n2)
 *                     ...          ...          ...  ...
 *                     v1(n1)*v2(1)  v1(n1)*v2(2)  ...  v1(n1)*v2(n2)]
 *
 * where v1 and v2 are column matrices of size n1 and n2 respectively,
 * and the output matrix has size n1*n2 x n1*n2.
 *
 * @param[in] v1        First input column matrix.
 * @param[in] v2        Second input column matrix.
 * @param[out] m_result Output matrix to store the result.
 *
 */
Matrix kron(const Matrix &v1, const Matrix &v2, Matrix &m_result);

/**
 * Copies the values in a row vector to a column vector.
 *
 * @param[in] srcRowVec    Row vector containing the values to be copied.
 * @param[in,out] destColVec  Column vector to receive the copied values.
 *
 * @throw std::invalid_argument if the input matrices have incompatible dimensions.
 */
void assign_rowVec_to_colVec(const Matrix &srcRowVec, Matrix &destColVec);

/**
 * Computes the sum of two cross products and stores the result
 * in an output matrix.
 *
 * @param[in] a         First input column matrix.
 * @param[in] b         Second input column matrix.
 * @param[in] c         Third input column matrix.
 * @param[in] d         Fourth input column matrix.
 * @param[in,out] tmp_f Temporary matrix for storing the first cross product.
 * @param[in,out] tmp_l Temporary matrix for storing the second cross product.
 * @param[out] v_result Output matrix to store the final result.
 *
 * @throw std::invalid_argument if the input matrices have incompatible dimensions.
 */
void a_cross_b_plus_c_cross_d(const Matrix &a,
                              const Matrix &b,
                              const Matrix &c,
                              const Matrix &d,
                              Matrix &tmp_f,
                              Matrix &tmp_l,
                              Matrix &v_result);