/**
 * @file test_fluid_dynamics.hpp
 * @brief Gates for fluid-mode dynamics: in-plane motion, the edge spring, and
 * the sparse valence-aware surface solver (WP4).
 */
#pragma once

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "Dynamics.hpp"
#include "dynamics/Surface_solver.hpp"
#include "mesh/Mesh.hpp"
#include "Parameters.hpp"
