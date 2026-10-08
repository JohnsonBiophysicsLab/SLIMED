/**
 * @file test_lipid_flow.hpp
 * @brief PROTOTYPE -- lipid flow as a relaxed area per face that moves
 * (Param::lipidFlow*, docs/lipid_flow_plan.md).
 */
#pragma once

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "Dynamics.hpp"
#include "Parameters.hpp"
#include "io/io.hpp"
#include "model/Record.hpp"
