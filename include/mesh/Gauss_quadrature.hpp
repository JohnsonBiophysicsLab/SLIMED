/**
 * @file Gauss_quadrature.hpp
 * @author Y Ying
 * @author Y Fu
 * @brief This file defines gauss quadrature and related
 * functions.
 * @version 0.1
 * @date 2023-03-10
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>

// matrix math
#include "linalg/Linear_algebra.hpp"

/**
 * @brief Calculate the VMU (barycentric coordinates) and weight for a Gaussian quadrature point.
 *
 * This function calculates a matrix (N, 3) of VMUs based on the specified number of quadrature points,
 * storing the result in the provided matrix `vmu`. It also calculates a matrix (N, 1) of weights for
 * each point, storing the result in the provided matrix `weight`.
 *
 * Note that the `vmu` and `weight` matrices will be overwritten and should be initialized before calling
 * this function.
 *
 * @param gaussQuadratureN The number of Gauss quadrature points to use.
 * @param vwu A matrix for storing the calculated barycentric coordinates (N, 3).
 * @param weight A matrix for storing the calculated weights (N, 1).
 */
void get_gauss_quadrature_weight_VWU(const int &gaussQuadratureN, Matrix &vwu, Matrix &weight);

/**
 * @brief Calculate the shape function and its derivatives at a given barycentric coordinate.
 *
 * This function calculates the shape function and its derivatives at a given barycentric coordinate `vwu`,
 * storing the result in the provided matrix `shapefunction`. The shape function and its derivatives are
 * stored in rows 0 to 6 of the matrix as follows:
 *
 *   shapefunction(0,:), shape functions;
 *   shapefunction(1,:), differential to v;
 *   shapefunction(2,:), differential to w;
 *   shapefunction(3,:), double differential to v;
 *   shapefunction(4,:), double differential to w;
 *   shapefunction(5,:), differential to v and w;
 *   shapefunction(6,:), differential to w and v;
 *
 * @param vwu The barycentric (1, 3) coordinate to evaluate the shape function at.
 * @param shapefunction A (7, 12) matrix for storing the calculated shape function and its derivatives (7, 12).
 */
void get_shapefunction(const Matrix &vwu, Matrix &shapefunction);

/**
 * @brief Calculate the shape function and its derivatives at a set of barycentric coordinates.
 *
 * This function calculates the shape function and its derivatives at each barycentric coordinate in the
 * provided matrix `vwuMat`, storing the result in the provided vector `sfVec`. Each element of `sfVec`
 * is a matrix representing the shape function and its derivatives at one barycentric coordinate.
 *
 * @param vwuMat A matrix of barycentric coordinates (N,3).
 * @param sfVec A vector for storing the calculated shape functions and their derivatives.
 */
void get_shapefunction_vector(const Matrix &vwuMat, std::vector<Matrix> &sfVec);

/**
 * @brief Calculate the four subdivision matrices for an irregular patch.
 *
 * This function calculates the four subdivision matrices for an irregular patch using the provided matrix `mat`,
 * storing the result in the provided matrices `subMat1`, `subMat2`, `subMat3`, and `subMat4`.
 *
 * @param mat The input matrix from which to compute the subdivision matrices.
 * @param subMat1 A matrix for storing the first subdivision matrix.
 * @param subMat2 A matrix for storing the second subdivision matrix.
 * @param subMat3 A matrix for storing the third subdivision matrix.
 * @param subMat4 A matrix for storing the fourth subdivision matrix.
 */
void get_subdivision_matrices(Matrix &mat,
                                     Matrix &subMat1,
                                     Matrix &subMat2,
                                     Matrix &subMat3,
                                     Matrix &subMat4);