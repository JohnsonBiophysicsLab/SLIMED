/**
 * @file test_truncated_noise.hpp
 * @brief The truncated Brownian noise of the vertex step
 * (Param::noiseTruncation, docs/resolution_limits.md section 2.5).
 */
#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>

#include "Counter_rng.hpp"
#include "Dynamics.hpp"
#include "Parameters.hpp"
#include "io/io.hpp"
#include "model/Record.hpp"
