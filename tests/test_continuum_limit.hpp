/**
 * @file test_continuum_limit.hpp
 * @brief How far up in q the discrete bending energy still behaves like the
 *        continuum Helfrich one, measured rather than assumed.
 *
 * Used by docs/fluctuation_spectrum.md and by the --qdx-max default in
 * analysis/fluctuation_report.py.
 */
#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <vector>

#include "mesh/Mesh.hpp"
#include "Parameters.hpp"
