/**
 * @file continuum_membrane.cpp
 * @author Y Ying (yying7@jh.edu)
 * @brief One instance of running lowest energy search on a flat membrane model.
 * 
 * @date 2023-04-04
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "Run_simulation.hpp"

using namespace std;

int main()
{
    // ProfilerStart("output.prof");
    run_flat("./input");
    // ProfilerStop();
}