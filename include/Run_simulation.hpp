/**
 * @file Run_simulation.hpp
 * @author Y Fu (yfu31@jh.edu)
 * @author Y Ying (yying7@jh.edu)
 * @brief This file includes all main simulation function depending
 * on the geometry of the membrane. (Currently only flat)
 * @date 2023-04-04
 *
 * @copyright Copyright (c) 2023
 *
 */
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>   // for accumulate vector
#include <algorithm>
// linear algebra
#include "linalg/Linear_algebra.hpp"
// triangular mesh
#include "mesh/Vertex.hpp"
#include "mesh/Face.hpp"
#include "mesh/Mesh.hpp"
#include "energy_force/Energy.hpp"
#include "energy_force/Force.hpp"
#include "mesh/Gauss_quadrature.hpp"
#include "model/Model.hpp"
#include "model/Record.hpp"
#include "Dynamics.hpp"
// io
#include "io/io.hpp"
// parameters
#include "Parameters.hpp"

using namespace std;

/**
 * @brief Run continuum membrane model on a flat membrane
 * 
 * @param param_filename 
 */
void run_flat(std::string param_filename);

/**
 * @brief Run membrane dynamics model on a flat membrane
 * 
 * @param param_filename 
 */
void run_dynamics_flat(std::string param_filename);