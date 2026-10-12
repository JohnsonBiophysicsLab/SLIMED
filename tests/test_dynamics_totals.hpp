/**
 * @file test_dynamics_totals.hpp
 * @brief The logged totals of a Brownian run describe the Hamiltonian the run
 * samples: the regularization slot counts toward E_Total and the mean force
 * only when it is in the drive.
 */
#pragma once

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "Dynamics.hpp"
#include "Parameters.hpp"
#include "energy_force/Energy.hpp"
#include "energy_force/Force.hpp"
#include "io/io.hpp"
#include "mesh/Mesh.hpp"
#include "model/Model.hpp"
#include "model/Record.hpp"
