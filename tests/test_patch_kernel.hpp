/**
 * @file test_patch_kernel.hpp
 * @brief Pins the flat patch kernel against the Matrix-based original.
 */
#pragma once

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <vector>

#include "energy_force/Patch_kernel.hpp"
#include "energy_force/Patch_rows_flat.hpp"
#include "linalg/Linear_algebra.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Subdivision_matrices.hpp"
#include "Parameters.hpp"
