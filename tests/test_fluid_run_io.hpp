/**
 * @file test_fluid_run_io.hpp
 * @brief Gates for wiring a fluid run into the driver: the depth scale, the
 * per-frame connectivity output, and the restart checkpoint's faces block
 * (WP5).
 */
#pragma once

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "Dynamics.hpp"
#include "io/io.hpp"
#include "mesh/Irregular_patch_rows.hpp"
#include "mesh/Mesh.hpp"
#include "model/Model.hpp"
#include "model/Record.hpp"
#include "Parameters.hpp"
