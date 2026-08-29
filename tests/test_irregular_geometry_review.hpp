/**
 * @file test_irregular_geometry_review.hpp
 * @brief Independent validation of the irregular limit-surface geometry.
 *
 * Everything else checks the row table against the recursion it replaces, or
 * against the literal matrices it was generated from. These check the geometry
 * against properties that hold regardless of how it is computed.
 */
#pragma once

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "mesh/Mesh.hpp"
#include "mesh/Subdivision_matrices.hpp"
#include "mesh/Gauss_quadrature.hpp"
#include "Parameters.hpp"
