/**
 * @file test_external_force.hpp
 * @brief The external force on a region of the membrane
 * (Param::externalForce*, docs/external_force.md).
 */
#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "Dynamics.hpp"
#include "Parameters.hpp"
#include "io/io.hpp"
#include "model/Record.hpp"
