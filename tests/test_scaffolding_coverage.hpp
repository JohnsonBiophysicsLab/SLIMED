/**
 * @file test_scaffolding_coverage.hpp
 * @brief Unit tests for src/energy_force/Scaffolding_points.cpp: the membrane
 * cap initialisation, the scaffold-membrane harmonic bond, and the Gag /
 * idealized-protein-lattice internal energy with its reference-state and
 * reaction-file parsers.
 *
 * @copyright Copyright (c) 2026
 */
#pragma once

#include <gtest/gtest.h>

#include <array>
#include <string>
#include <vector>

#include "io/io.hpp"
#include "linalg/Linear_algebra.hpp"
#include "mesh/Face.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "Parameters.hpp"
