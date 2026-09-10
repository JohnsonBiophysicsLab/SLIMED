/**
 * @file test_multi_extraordinary.hpp
 * @brief Gates for evaluating faces with several extraordinary corners (WP1).
 */
#pragma once

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "mesh/Mesh.hpp"
#include "mesh/Multi_extraordinary_patch.hpp"
#include "energy_force/Patch_kernel.hpp"
#include "mesh/Subdivision_matrices.hpp"
#include "Parameters.hpp"
