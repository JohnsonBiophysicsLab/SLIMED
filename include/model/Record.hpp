#pragma once

#include <vector>
#include <iostream>
#include "energy_force/Energy.hpp"

/**
 * @brief Bookkeeps simulation data throughout iterations.
 *
 * This class stores vectors of surface area, energy, and mean force values for each iteration of a simulation.
 */
class Record
{
public:
    std::vector<double> areaTotal; ///< Vector of total surface area values.
    std::vector<Energy> energyVec; ///< Vector of energy values.
    std::vector<double> meanForce; ///< Vector of mean force values.

    /**
     * @brief Constructor for Record class.
     */
    Record();

    /**
     * @brief Constructor for Record class that pre-allocates memory for the data vectors.
     *
     * @param numIterations The number of iterations to allocate memory for.
     */
    Record(const int numIterations);

    /**
     * @brief Adds a new record to the data vectors.
     *
     * This function adds a new record to each of the `areaTotal`, `energyVec`, and `meanForce` vectors. The values added
     * correspond to the current area, energy, and mean force, respectively.
     *
     * @param currentArea The current surface area of the mesh.
     * @param currentEnergy The current energy of the mesh.
     * @param currentMeanForce The current mean force on the mesh vertices.
     */
    void add(const double currentArea, const Energy &currentEnergy, const double currentMeanForce);

    /**
     * @brief Sets the values for a specific record in the data vectors.
     *
     * This function sets the values for a specific record in each of the `areaTotal`, `energyVec`, and `meanForce` vectors.
     * The index of the record to set is given by the `iteration` parameter, and the values to set are given by the
     * `currentArea`, `currentEnergy`, and `currentMeanForce` parameters.
     *
     * @param iteration The index of the record to set.
     * @param currentArea The surface area for the record.
     * @param currentEnergy The energy for the record.
     * @param currentMeanForce The mean force for the record.
     */
    void set(const int iteration, const double currentArea, const Energy &currentEnergy, const double currentMeanForce);

    /**
     * @brief Prints information about a single iteration of the optimization algorithm.
     *
     * This function prints information about a single iteration of the optimization algorithm, including the current step
     * number (`iteration`), the total energy of the mesh at that step (`energyVec[iteration].energyTotal`), the mean force
     * on the mesh vertices at that step (`meanForce[iteration]`), and the surface area of the mesh at that step
     * (`areaTotal[iteration]`).
     *
     * The output is sent to standard output (i.e., to the terminal).
     *
     * @param iteration The step number of the iteration to print information about.
     */
    void print_iteration(const int iteration);
};