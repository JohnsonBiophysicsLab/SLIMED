/**
 * @file continuum_membrane_multithreading.cpp
 * @author Y Ying (yying7@jh.edu)
 * @brief Embarrassingly parallel implementation of membrane dynamics model.
 * @date 2023-04-14
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include <iostream>
#include <vector>
#include <string>
#include <thread>

#include "Run_simulation.hpp"

int main() {
    // Define a vector of input files
    std::vector<std::string> inputs = {"input1", "input2", "input3", "input4", "input5", "input6", "input7", "input8"};

    // Create an empty vector of threads
    std::vector<std::thread> threads;

    // Loop over each input file in the vector
    for (const auto& input : inputs) {
        // Create a new thread and pass it the run_flat function and the input file
        threads.emplace_back(std::thread(run_dynamics_flat, input));
    }

    // Loop over all the threads and wait for them to finish
    for (auto& thread : threads) {
        thread.join();
    }

    return 0;
}