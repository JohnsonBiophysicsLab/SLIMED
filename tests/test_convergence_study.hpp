/**
 * @file test_convergence_study.hpp
 * @brief WP6: per-valence convergence of area, volume, bending energy and
 * force-energy conjugacy against subdivision depth.
 *
 * @note WP6 of docs/irregular_patch_valence_4_to_8_plan.md.
 */
#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <vector>

#include "mesh/Mesh.hpp"
#include "mesh/Subdivision_matrices.hpp"
#include "Parameters.hpp"
