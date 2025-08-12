#include "linalg/Linear_algebra.hpp"

using namespace std;

Matrix::Matrix() : mat(NULL) {}

Matrix::Matrix(const int &nrow) : mat(gsl_matrix_alloc(nrow, 1))
{
}

Matrix::Matrix(const int &nrow, const int &ncol) : mat(gsl_matrix_alloc(nrow, ncol))
{
}

Matrix::Matrix(const Matrix &mat_src) : mat(gsl_matrix_alloc(mat_src.mat->size1, mat_src.mat->size2))
{
    gsl_matrix_memcpy(mat, mat_src.mat);
}

Matrix& Matrix::operator=(const Matrix &mat_src)
{
    if (this == &mat_src) {
        return *this;  // Avoid self-assignment
    }
    if (mat != NULL) {
        gsl_matrix_free(mat);  // Deallocate existing memory
        mat = NULL;
    }
    mat = gsl_matrix_alloc(mat_src.mat->size1, mat_src.mat->size2);  // Allocate new memory
    gsl_matrix_memcpy(mat, mat_src.mat);  // Copy data
    return *this;
}

Matrix::Matrix(const std::vector<std::vector<double>> &vec) : mat(gsl_matrix_alloc(vec.size(), vec[0].size()))
{
    /*
    It is more efficient to cache
    vec.size() and  vec[0].size() in variables rather than calling
    them repeatedly in the loops. Also, using range-based for loops
    can make the code simpler and more readable. -GPT
    */

    int nrow = vec.size();
    int ncol = vec[0].size();

    for (int i = 0; i < nrow; ++i)
    {
        for (int j = 0; j < ncol; ++j)
        {
            gsl_matrix_set(mat, i, j, vec[i][j]);
        }
    }
}

Matrix::Matrix(const int &nrow, const int &ncol, const bool &setToZero) : mat(gsl_matrix_calloc(nrow, ncol))
{
}

Matrix::~Matrix(){
    if (mat != NULL){
        gsl_matrix_free(mat);
        mat = NULL;
    }
}

int Matrix::nrow() const
{
    return static_cast<int>(mat->size1);
}

int Matrix::ncol() const
{
    return static_cast<int>(mat->size2);
}

double Matrix::get(const int &i, const int &j) const
{
    return gsl_matrix_get(mat, i, j);
}

double Matrix::operator()(const int &i, const int &j) const
{
    return gsl_matrix_get(mat, i, j);
}

void Matrix::set(const int &i, const int &j, const double &v)
{
    gsl_matrix_set(mat, i, j, v);
}

void Matrix::set_all(const double &v)
{
    // fill the matrix with the specified value
    gsl_matrix_set_all(mat, v);
}

void Matrix::transpose()
{
    if (ncol() == nrow()){
        gsl_matrix_transpose(mat);
    }
    else{
        Matrix tMat(mat->size2, mat->size1);
        for (int i = 0; i < mat->size1; i++){
            for (int j =0; j < mat->size2; j++){
                tMat.set(j, i, this->get(i, j));
            }
        }
    }
}

void Matrix::set_identity()
{
    gsl_matrix_set_identity(mat);
}

Matrix Matrix::get_row(const int &irow) const
{
    const size_t num_cols = mat->size2;
    Matrix rowMat(1, num_cols);
    // use submatrix function to get part of the matrix
    gsl_matrix_const_view row_view = gsl_matrix_const_submatrix(mat, irow, 0, 1, num_cols);
    gsl_matrix_memcpy(rowMat.mat, &row_view.matrix);
    return rowMat;
}

Matrix Matrix::get_col(const int &icol) const
{
    const size_t num_rows = mat->size1;
    Matrix colMat(num_rows, 1);
    // use submatrix function to get part of the matrix
    gsl_matrix_const_view col_view = gsl_matrix_const_submatrix(mat, 0, icol, num_rows, 1);
    gsl_matrix_memcpy(colMat.mat, &col_view.matrix);
    return colMat;
}

void Matrix::set_row_from_col(const int &destRow, const Matrix &srcMat, const int &srcCol)
{
    gsl_vector *v = gsl_vector_alloc(srcMat.nrow()); // Allocate memory for vector
    gsl_matrix_get_col(v, srcMat.mat, srcCol); // get row as gsl_vector*
    gsl_matrix_set_row(mat, destRow, v); // set row from vector
    gsl_vector_free(v); // free vector memory
    v = NULL;
}

Matrix Matrix::get_transposed() const
{
    const size_t nrow = mat->size1;
    const size_t ncol = mat->size2;
    
    
        Matrix transMat(ncol, nrow);
    gsl_matrix_transpose_memcpy(transMat.mat, mat);
    return transMat;

}

double Matrix::calculate_quad() const
{
    double squared_sum = 0.0;
    for (int i = 0; i < mat->size1; ++i)
    {
        for (int j = 0; j < mat->size2; ++j)
        {
            double a_ij = gsl_matrix_get(mat, i, j);
            squared_sum += a_ij * a_ij;
        }
    }
    return squared_sum;
}

double Matrix::calculate_norm() const
{
    double norm = sqrt(this->calculate_quad());
    return norm;
}

void Matrix::get_inverted(Matrix& matInverted)
{
    // Deepcopy the input matrix assuming square matrix
    int size = mat->size1; //assume equal rows and cols
    matInverted = Matrix(size, size);

    gsl_matrix* matrix_copy = gsl_matrix_calloc(size, size);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            gsl_matrix_set(matrix_copy, i, j, gsl_matrix_get(mat, i, j));
        }
    }
    gsl_permutation *p = gsl_permutation_alloc(size);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix_copy, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_linalg_LU_invert(matrix_copy, p, matInverted.mat);

    gsl_permutation_free(p);

    
}

void Matrix::free()
{
    if (mat != NULL)
        gsl_matrix_free(mat);
        mat = NULL;
}

Matrix mat_calloc(const int &nrow, const int &ncol)
{
    // create a new Matrix object with dimensions (nrow, ncol)
    Matrix mat(nrow, ncol, true);

    // return the new matrix
    return mat;
}

// Output operator for Matrix objects
std::ostream &operator<<(std::ostream &stream, const Matrix &matrix)
{
    if (matrix.ncol() != 1)
    {
        for (int i = 0; i < matrix.nrow(); ++i)
        {
            for (int j = 0; j < matrix.ncol(); ++j)
            {
                stream << gsl_matrix_get(matrix.mat, i, j) << ", ";
            }
            stream << std::endl;
        }
    } else {
        for (int i = 0; i < matrix.nrow(); ++i)
        {
            stream << gsl_matrix_get(matrix.mat, i, 0) << ", ";
        }
    }
    return stream;
}

// Addition assignment operator for Matrix objects
Matrix &operator+=(Matrix &m1, const Matrix &m2)
{
    gsl_matrix_add(m1.mat, m2.mat);
    return m1;
}

// Addition operator for Matrix objects
Matrix operator+(const Matrix &m1, const Matrix &m2)
{
    Matrix result(m1);
    result += m2;
    return result;
}

// Subtraction assignment operator for Matrix objects
Matrix &operator-=(Matrix &m1, const Matrix &m2)
{
    gsl_matrix_sub(m1.mat, m2.mat);
    return m1;
}

// Subtraction operator for Matrix objects
Matrix operator-(const Matrix &m1, const Matrix &m2)
{
    Matrix result(m1);
    result -= m2;
    return result;
}

// Scalar multiplication assignment operator for Matrix objects
Matrix &operator*=(Matrix &mat, const double &scalar)
{
    gsl_matrix_scale(mat.mat, scalar);
    return mat;
}

// Scalar multiplication operator for Matrix objects
Matrix operator*(const Matrix &matrix, const double &scalar)
{
    Matrix result(matrix);
    result *= scalar;
    return result;
}

// Scalar multiplication operator for Matrix objects (scalar first)
Matrix operator*(const double &scalar, const Matrix &matrix)
{
    return matrix * scalar;
}

// Matrix multiplication assignment operator for Matrix objects
Matrix &operator*=(Matrix &m1, const Matrix &m2)
{
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m1.mat, m2.mat, 0.0, m1.mat);
    return m1;
}

// Matrix multiplication operator for Matrix objects
Matrix operator*(const Matrix &m1, const Matrix &m2)
{
    Matrix result(m1.nrow(), m2.ncol());
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m1.mat, m2.mat, 0.0, result.mat);
    return result;
}

// Dot product function for Matrix objects
Matrix dot(const Matrix &m1, const Matrix &m2)
{
    Matrix result(1, m1.ncol());
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, m1.mat, m2.mat, 0.0, result.mat);
    return result;
}

// Scalar division assignment operator for Matrix objects
Matrix &operator/=(Matrix &mat, const double &scalar)
{
    gsl_matrix_scale(mat.mat, 1.0 / scalar);
    return mat;
}

// Scalar division operator for Matrix objects
Matrix operator/(Matrix matrix, const double &scalar)
{
    matrix /= scalar;
    return matrix;
}

double dot_col(const Matrix &m1, const Matrix &m2)
{
    double result = 0.0;
    for (int i = 0; i < m1.mat->size1; i++)
        result += gsl_matrix_get(m1.mat, i, 0) * gsl_matrix_get(m2.mat, i, 0);
    return result;
}

double dot_row(const Matrix &m1, const Matrix &m2)
{
    double result = 0.0;
    for (int i = 0; i < m1.mat->size1; i++)
        result += gsl_matrix_get(m1.mat, 0, i) * gsl_matrix_get(m2.mat, 0, i);
    return result;
}

Matrix cross_col(const Matrix &m1, const Matrix &m2)
{

    Matrix result(3, 1);

    gsl_matrix_set(result.mat, 0, 0, gsl_matrix_get(m1.mat, 1, 0) * gsl_matrix_get(m2.mat, 2, 0) - gsl_matrix_get(m1.mat, 2, 0) * gsl_matrix_get(m2.mat, 1, 0));
    gsl_matrix_set(result.mat, 1, 0, gsl_matrix_get(m1.mat, 2, 0) * gsl_matrix_get(m2.mat, 0, 0) - gsl_matrix_get(m1.mat, 0, 0) * gsl_matrix_get(m2.mat, 2, 0));
    gsl_matrix_set(result.mat, 2, 0, gsl_matrix_get(m1.mat, 0, 0) * gsl_matrix_get(m2.mat, 1, 0) - gsl_matrix_get(m1.mat, 1, 0) * gsl_matrix_get(m2.mat, 0, 0));

    return result;
}

Matrix cross_row(const Matrix &m1, const Matrix &m2)
{
    Matrix result(1, 3);

    gsl_matrix_set(result.mat, 0, 0, gsl_matrix_get(m1.mat, 0, 1) * gsl_matrix_get(m2.mat, 0, 2) - gsl_matrix_get(m1.mat, 0, 2) * gsl_matrix_get(m2.mat, 0, 1));
    gsl_matrix_set(result.mat, 0, 1, gsl_matrix_get(m1.mat, 0, 2) * gsl_matrix_get(m2.mat, 0, 0) - gsl_matrix_get(m1.mat, 0, 0) * gsl_matrix_get(m2.mat, 0, 2));
    gsl_matrix_set(result.mat, 0, 2, gsl_matrix_get(m1.mat, 0, 0) * gsl_matrix_get(m2.mat, 0, 1) - gsl_matrix_get(m1.mat, 0, 1) * gsl_matrix_get(m2.mat, 0, 0));

    return result;
}

/******************************/
/* High Performance functions */
/******************************/

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
void negative(const Matrix &m1, Matrix &m_neg)
{

    // deep copy m1 to m_neg before scaling
    gsl_matrix_memcpy(m_neg.mat, m1.mat);

    // scale m_neg
    gsl_matrix_scale(m_neg.mat, -1);
}

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
void get_unit_vector(Matrix &m)
{
    // Compute the magntidue of the input vector
    double mag = m.calculate_norm();

    // Compute the unit vector and store it in m
    gsl_matrix_scale(m.mat, 1.0 / mag);
}

/**
 * @brief Get the unit vector in the direction of m1 which is a (3, 1) Matrix
 *
 * This function computes the unit vector in the direction of a given input vector,
 * which is assumed to be a 3x1 matrix. The result is stored in another 3x1 matrix.
 *
 * @param m1 The input vector to compute the unit vector of.
 * @param m_unit The output matrix to store the resulting unit vector.
 *
 * @pre The input matrix 'm1' must be a 3x1 matrix.
 *
 * @post The output matrix 'm_unit' contains the unit vector in the direction of 'm1'.
 */
void get_unit_vector(const Matrix &m1, Matrix &m_unit)
{
    // Compute magnitude of the input vector
    double mag = m1.calculate_norm();

    // Compute the unit vector and store it in m_unit
    gsl_matrix_memcpy(m_unit.mat, m1.mat);
    gsl_matrix_scale(m_unit.mat, 1.0 / mag);
}

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
void cross(const Matrix &m1, const Matrix &m2, Matrix &temp)
{
    double result[3];
    result[0] = gsl_matrix_get(m1.mat, 1, 0) * gsl_matrix_get(m2.mat, 2, 0) - gsl_matrix_get(m1.mat, 2, 0) * gsl_matrix_get(m2.mat, 1, 0);
    result[1] = gsl_matrix_get(m1.mat, 2, 0) * gsl_matrix_get(m2.mat, 0, 0) - gsl_matrix_get(m1.mat, 0, 0) * gsl_matrix_get(m2.mat, 2, 0);
    result[2] = gsl_matrix_get(m1.mat, 0, 0) * gsl_matrix_get(m2.mat, 1, 0) - gsl_matrix_get(m1.mat, 1, 0) * gsl_matrix_get(m2.mat, 0, 0);

    gsl_matrix_set(temp.mat, 0, 0, result[0]);
    gsl_matrix_set(temp.mat, 1, 0, result[1]);
    gsl_matrix_set(temp.mat, 2, 0, result[2]);
}

/**
 * @brief Compute the element-wise addition of two matrices and store the result in a third matrix.
 *
 * This function computes the element-wise addition of two input matrices of equal dimensions
 * and stores the result in a third matrix. All input matrices are assumed to be of type Matrix.
 *
 *
 * @param m1 The first input matrix.
 * @param m2 The second input matrix.
 * @param tmp The output matrix to store the resulting sum.
 *
 * @pre All input matrices must have the same dimensions.
 *
 * @post The output matrix 'tmp' contains the element-wise sum of 'm1' and 'm2'.
 */
void addition(const Matrix &m1, const Matrix &m2, Matrix &tmp)
{
    gsl_matrix_memcpy(tmp.mat, m1.mat);
    gsl_matrix_add(tmp.mat, m2.mat);
}

/**
 * @brief Compute the element-wise subtraction of two matrices and store the result in a third matrix.
 *
 * This function computes the element-wise subtraction of two input matrices of equal dimensions
 * and stores the result in a third matrix. All input matrices are assumed to be of type Matrix.
 *
 *
 * @param m1 The first input matrix.
 * @param m2 The second input matrix.
 * @param tmp The output matrix to store the resulting difference.
 *
 * @pre All input matrices must have the same dimensions.
 *
 * @post The output matrix 'tmp' contains the element-wise difference of 'm1' and 'm2'.
 */
void subtraction(const Matrix &m1, const Matrix &m2, Matrix &tmp)
{
    for (size_t i = 0; i < m1.mat->size1; i++)
    {
        for (size_t j = 0; j < m1.mat->size2; j++)
        {
            double value = gsl_matrix_get(m1.mat, i, j) - gsl_matrix_get(m2.mat, i, j);
            gsl_matrix_set(tmp.mat, i, j, value);
        }
    }
}

/**
 * @brief Compute the matrix multiplication of two matrices and store the result in a third matrix.
 *
 * This function computes the matrix multiplication of two input matrices and stores the result in
 * a third matrix. All input matrices are assumed to be of type Matrix.
 *
 * @param m1 The first input matrix.
 * @param m2 The second input matrix.
 * @param tmp The output matrix to store the resulting product.
 *
 * @pre The number of columns of 'm1' must match the number of rows of 'm2'.
 *
 * @post The output matrix 'tmp' contains the matrix product of 'm1' and 'm2'.
 */
void multiplication(const Matrix &m1, const Matrix &m2, Matrix &tmp)
{
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m1.mat, m2.mat, 0.0, tmp.mat);
}

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
void const_multiplication(const Matrix &m1, const double num, Matrix &tmp)
{
    gsl_matrix_memcpy(tmp.mat, m1.mat);
    gsl_matrix_scale(tmp.mat, num);
}

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
void const_division(const Matrix &m1, const double num, Matrix &tmp)
{
    gsl_matrix_memcpy(tmp.mat, m1.mat);
    double reciprocal = 1.0 / num;
    gsl_matrix_scale(tmp.mat, reciprocal);
}

/**
 * @brief Multiply a row vector by a matrix and store the result in a row vector.
 *
 * This function multiplies a row vector tranposed from input column vector
 * with a matrix and stores the resulting product in
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
void colvec_matrix_multiplication(const Matrix &v1, const Matrix &m1, Matrix &tmp)
{
    double sum = 0.0;
    for (int i = 0; i < m1.mat->size2; i++){
        sum = 0.0;
        for (int j = 0; j < m1.mat->size1; j++){
            sum += v1.get(j, 0) * m1.get(j, i);
        }
        tmp.set(i, 0, sum);
    }
}

/**
 * @brief the Kronecker product between two column vectors and return the
 * result. The Kronecker product is defined as follows:
 *
 *     kron(v1, v2) = [v1(1)*v2(1)  v1(1)*v2(2)  ...  v1(1)*v2(n2)
 *                     v1(2)*v2(1)  v1(2)*v2(2)  ...  v1(2)*v2(n2)
 *                     ...          ...          ...  ...
 *                     v1(n1)*v2(1)  v1(n1)*v2(2)  ...  v1(n1)*v2(n2)]
 *
 * where v1 and v2 are column matrices of size n1 and n2 respectively,
 * and the output matrix has size n1*n2 x n1*n2.
 *
 * @param[in] v1        First input column matrix (M, 1).
 * @param[in] v2        Second input column matrix (N ,1).
 * @return Output matrix to store the result.
 * @note The two inputs are expected to be (M, N) column vectors
 *
 * @throw std::invalid_argument if the output matrix has incompatible dimensions.
 */
Matrix kron(const Matrix &v1, const Matrix &v2)
{
    Matrix m_result(v1.mat->size1, v2.mat->size1); // define a m x n matrix
    for (int i = 0; i < v1.mat->size1; i++)
    {
        for (int j = 0; j < v2.mat->size1; j++)
        {
            m_result.set(i, j, v1(i, 0) * v2(j, 0));
        }
    }
    return m_result;
}

/**
 * @brief the Kronecker product between two column vectors and return the
 * result. The Kronecker product is defined as follows:
 *
 *     kron(v1, v2) = [v1(1)*v2(1)  v1(1)*v2(2)  ...  v1(1)*v2(n2)
 *                     v1(2)*v2(1)  v1(2)*v2(2)  ...  v1(2)*v2(n2)
 *                     ...          ...          ...  ...
 *                     v1(n1)*v2(1)  v1(n1)*v2(2)  ...  v1(n1)*v2(n2)]
 *
 * where v1 and v2 are column matrices of size n1 and n2 respectively,
 * and the output matrix has size n1*n2 x n1*n2.
 *
 * @param[in] v1        First input column matrix (M, 1).
 * @param[in] v2        Second input column matrix (N ,1).
 * @return Output matrix to store the result.
 * @note The two inputs are expected to be (M, N) column vectors
 *
 * @throw std::invalid_argument if the output matrix has incompatible dimensions.
 */
Matrix kron(const Matrix &v1, const Matrix &v2, Matrix &m_result)
{
    for (int i = 0; i < v1.mat->size1; i++)
    {
        for (int j = 0; j < v2.mat->size1; j++)
        {
            m_result.set(i, j, v1(i, 0) * v2(j, 0));
        }
    }
    return m_result;
}

// Complex method
/**
 * Copies the values in a row vector to a column vector.
 *
 * @param[in] srcRowVec    Row vector containing the values to be copied.
 * @param[in,out] destColVec  Column vector to receive the copied values.
 *
 * @throw std::invalid_argument if the input matrices have incompatible dimensions.
 */
void assign_rowVec_to_colVec(const Matrix &srcRowVec, Matrix &destColVec)
{

    // Copy the values from srcRowVec to destColVec.
    for (int i = 0; i < srcRowVec.ncol(); i++)
    {
        destColVec.set(i, 0, srcRowVec.get(0, i));
    }
}

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
                              Matrix &v_result)
{
    // Compute the first cross product between a and b, and store the result in tmp_f.
    cross(a, b, tmp_f);

    // Compute the second cross product between c and d, and store the result in tmp_l.
    cross(c, d, tmp_l);

    // Add the two cross products together and store the result in v_result.
    addition(tmp_f, tmp_l, v_result);
}