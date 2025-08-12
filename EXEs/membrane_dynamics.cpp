/**
 * @file membrane_dynamics.cpp
 * @author Y Ying (yying7@jh.edu)
 * @brief One instance of membrane dynamics model
 * @date 2023-04-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "Run_simulation.hpp"

//maincode

int main() {
    //ProfilerStart("output.prof");
    run_dynamics_flat("input");
    //ProfilerStop();
}
