/**
 * @file test_irregular_patch_quadrature.hpp
 * @brief The irregular-patch quadrature parameters: per-valence depths and
 * the tail correction.
 *
 * See docs/irregular_patch_results.md section 9 for the measurements the
 * defaults rest on.
 */
#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "energy_force/Patch_kernel.hpp"
#include "io/io.hpp"
#include "mesh/Irregular_patch_rows.hpp"
#include "mesh/Mesh.hpp"
#include "Parameters.hpp"
