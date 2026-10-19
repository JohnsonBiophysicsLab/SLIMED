/**
 * @file test_run_drivers.hpp
 * @brief End-to-end tests of the two simulation drivers, run_flat() and
 * run_dynamics_flat(), on tiny sheets.
 *
 * Both drivers read "<name>.params" from the working directory and write
 * fixed-name outputs into it, so every test here runs inside its own fresh
 * temporary directory and restores the original working directory afterwards.
 *
 * @copyright Copyright (c) 2026
 */
#pragma once

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "Run_simulation.hpp"
