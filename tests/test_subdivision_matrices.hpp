/**
 * @file test_subdivision_matrices.hpp
 * @brief Tests for the generated Loop subdivision matrices.
 *
 * @note WP2 of docs/irregular_patch_valence_4_to_8_plan.md.
 */
#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <vector>

#include "mesh/Subdivision_matrices.hpp"
#include "mesh/Gauss_quadrature.hpp"
#include "subdivision_oracle.hpp"
