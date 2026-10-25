/**
 * @file test_incremental_flip_forces.hpp
 * @brief Gates for the post-flip force pass that re-integrates only the faces
 * the accepted flips touched (Param::edgeFlipIncrementalForces).
 */
#pragma once

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "Dynamics.hpp"
#include "Parameters.hpp"
#include "io/io.hpp"
#include "mesh/Mesh.hpp"
