/**
 * @file test_core_coverage.hpp
 * @brief Unit tests for the small core types that the simulation tests only
 * touch indirectly: the Param printer and boundary-name helpers
 * (src/parameters/Parameters.cpp), the Force container
 * (src/energy_force/Force.cpp), and the rarely used Matrix operators in
 * src/linear_algebra/Linear_algebra.cpp.
 */
#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "Parameters.hpp"
#include "energy_force/Force.hpp"
#include "linalg/Linear_algebra.hpp"
