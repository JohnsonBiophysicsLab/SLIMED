/**
 * @file test_volume_functional.hpp
 * @brief Permanent tests for the volume functional: origin independence and
 * force-energy conjugacy.
 *
 * @note Step 5 of docs/volume_functional_split.md.
 */
#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "mesh/Mesh.hpp"
#include "model/Model.hpp"
#include "model/Record.hpp"
