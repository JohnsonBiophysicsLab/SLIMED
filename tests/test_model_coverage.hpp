/**
 * @file test_model_coverage.hpp
 * @brief Unit tests for src/model/: Record, OptimizationAlgorithm, the Model's
 * stopping rule, trial step size, NCG direction bookkeeping, line search,
 * boundary enforcement after a coordinate update, and the thermal trial move.
 */
#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "Parameters.hpp"
#include "energy_force/Energy.hpp"
#include "energy_force/Force.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "model/Model.hpp"
#include "model/Record.hpp"
