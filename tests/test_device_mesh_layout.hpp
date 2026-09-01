/**
 * @file test_device_mesh_layout.hpp
 * @brief Pins the flattened, GPU-shaped force pipeline against the production one.
 */
#pragma once

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <vector>

#include "cuda/Cuda_force_backend.hpp"
#include "cuda/Device_mesh_layout.hpp"
#include "cuda/Host_force_backend.hpp"
#include "energy_force/Patch_rows_flat.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Subdivision_matrices.hpp"
#include "Parameters.hpp"
