/**
 * @file test_local_output.hpp
 * @brief The per-frame local output of the dynamics (localOutput): the
 * per-face energies and per-vertex forces beside each trajectory frame.
 */
#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "Dynamics.hpp"
#include "Parameters.hpp"
#include "io/io.hpp"
#include "model/Record.hpp"
