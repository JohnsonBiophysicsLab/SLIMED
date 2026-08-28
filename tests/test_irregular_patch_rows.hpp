/**
 * @file test_irregular_patch_rows.hpp
 * @brief Tests for the collapsed irregular-patch row table.
 *
 * @note WP3 of docs/irregular_patch_valence_4_to_8_plan.md.
 */
#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "mesh/Irregular_patch_rows.hpp"
#include "mesh/Gauss_quadrature.hpp"
#include "Parameters.hpp"
#include "subdivision_oracle.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Subdivision_matrices.hpp"
#include <algorithm>
