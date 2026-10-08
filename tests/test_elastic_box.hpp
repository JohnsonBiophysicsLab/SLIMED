/**
 * @file test_elastic_box.hpp
 * @brief The elastic periodic box (Param::elasticBox*, docs/elastic_box.md).
 */
#pragma once

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Dynamics.hpp"
#include "Parameters.hpp"
#include "io/io.hpp"
#include "model/Record.hpp"
