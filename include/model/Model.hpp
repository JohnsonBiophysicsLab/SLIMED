/**
 * @file Model.hpp
 * @author Y Ying (yying7@jh.edu)
 * @brief Contains the Model class, which encapsulates a mesh and a record object.
 * @date 2023-04-04
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <ctime>
#include <random>

#include "linalg/Linear_algebra.hpp"
#include "mesh/Mesh.hpp"
#include "model/Model.hpp"
#include "model/Record.hpp"
#include "Parameters.hpp"
#include "energy_force/Force.hpp"

/**
 *
 * @brief A class that defines an optimization algorithm and its settings.
 * This class defines the settings for an optimization algorithm, including whether to use
 * nonlinear conjugate gradient (NCG), whether NCG has gotten stuck and cannot further minimize,
 * and whether to use R-pi adaptive method for regularization energy. These settings can be accessed
 * and modified using the public member variables of the class.
 */
class OptimizationAlgorithm
{
public:
    bool usingNCG = false;                  ///< Flag indicating whether to use the nonlinear conjugate gradient optimizer.
    bool isNCGstuck = false;               ///< Flag indicating whether the NCG optimizer has gotten stuck and cannot further minimize.
    bool usingRpi = false;                  ///< Flag indicating whether to use the R-pi adaptive method for regularization energy.
    bool isCriteriaSatisfied = false;      ///< Flag indicating whether the optimization criteria have been met.
    double trialStepSize = 0.0;            ///< The step size used in the optimization algorithm.
    int trialIterationInterval = 300;       ///< The number of iterations between trial step size adjustments.
    double c1 = 1e-6;                      ///< The c1 parameter in the NCG Wolfe condition.
    double c2 = 0.001;                     ///< The c2 parameter in the NCG Wolfe condition.
    double stepThreshold = 1e-15;          ///< The lower bound for the step size in NCG.
    int nConsecutiveNcgStuck = 0;          ///< The number of consecutive times the NCG optimizer has gotten stuck.
    int nConsecutiveNcgStuckThreshold = 3; ///< The threshold of consecutive NCG stuck iterations before disabling it.
    double energyDiffThreshold = 1e-5;     ///< The energy difference threshold for convergence checking.
    double forceDiffThreshold = 1e-4;      ///< The mean force difference threshold for convergence checking.

    /**
     * @brief Constructor for the OptimizationAlgorithm class.
     *
     * Initializes an instance of the OptimizationAlgorithm class with default settings.
     * Sets usingNCG to true, isNCGstuck to false, and usingRpi to true.
     */
    OptimizationAlgorithm();

    /**
     * @brief Constructor for the OptimizationAlgorithm class.
     *
     * Initializes an instance of the OptimizationAlgorithm class with default settings.
     * Customize iteration interval, c1, c2, and step threshold.
     * Sets usingNCG to true, isNCGstuck to false, and usingRpi to true.
     */
    OptimizationAlgorithm(const int trialIterationInterval,
                          const double c1,
                          const double c2,
                          const double stepThreshold);

    /**
     * @brief Reset NCG and Rpi to default.
     *
     * Sets usingNCG to true, isNCGstuck to false, and usingRpi to true.
     */
    void reset_NCG_Rpi();

    /**
     * @brief Disables the use of NCG if it has been stuck consecutively over a threshold.
     *
     * If the NCG optimizer has been stuck for too many consecutive iterations (as determined by the `isNCGstuck` flag),
     * this function will increment a counter `nConsecutiveNcgStuck`. If that counter exceeds a threshold of nConsecutiveNcgStuckThreshold, the `usingNCG`
     * flag will be set to false, indicating that NCG should not be used. Otherwise, `usingNCG` remains true.
     *
     * If the NCG optimizer is not currently stuck, the `nConsecutiveNcgStuck` counter is reset to 0 and `usingNCG` is set to
     * true.
     *
     * @return True if NCG should be used, false otherwise.
     */
    bool disable_NCG_if_stuck_consecutively();
};

/**
 * @brief The Model class encapsulates a Mesh and a Record object.
 */
class Model
{
public:
    // Geometry, energy, and force
    Mesh& mesh;     /**< The mesh object. */
    Record& record; /**< The record object. */
    // Optimization
    int iteration;                    /**< The current iteration count of the optimization algorithm. */
    double stepSize;                  /**< The current step size */
    OptimizationAlgorithm oa;         /**< The state for the optimization algorithm. */
    std::vector<Force> ncgDirection0; ///< vector to store the initial direction for nonlinear conjugate gradient optimization

    // Simulated Annealing
    bool isHeating = true;
    int currentCoolingStep = 0;
    int coolingStep = 35;
    int currentHeatingStep = 0;
    int heatingStep = 15;
    double highTemperature = 0.0;

    /**
     * @brief Constructs a new Model object.
     *
     * @param mesh_ The mesh object to be encapsulated.
     * @param record_ The record object to be encapsulated.
     */
    Model(Mesh &mesh_, Record &record_);

    /**
     * @brief Determines whether the optimization algorithm should continue running.
     *
     * The function returns true if the model has not yet been optimized and the
     * maximum number of iterations has not been reached. Otherwise, the function
     * returns false and the optimization algorithm stops.
     *
     * @note The implementation of this function can be changed to customize the optimization.
     *
     * @return true if the optimization algorithm should continue running, false otherwise.
     */
    bool should_continue_optimization();

    /**
     * @brief This code determines the step size used in an optimization algorithm.
     * The step size is used to calculate how far to move in each iteration of the
     * optimization algorithm to minimize an energy function.
     *
     * @return double trialStepSize
     */
    void determine_trial_step_size();

    /**
     * @brief Returns a string representation of the current step's iteration number, trial step size, and optimal step size.
     * The function concatenates the values of model.iteration, model.oa.trialStepSize, and
     * model.stepSize to the following format: "step: , trial StepSize = , StepSize = ".
     * @return std::string A string representation of the current step's information.
     */
    std::string to_string_current_step();

    /**
     * @brief Updates the direction for nonlinear conjugate gradient (NCG) optimization.
     *
     * This function calculates the dot product of the previous and current forces to determine the NCG factor. It then updates
     * the NCG direction either as a combination of the current force and the previous NCG direction, or just the current force
     * depending on whether NCG optimization is being used. The NCG direction is updated in the input vector `ncgDirections`.
     *
     */
    void update_ncg_direction();

    /**
     * @brief Reset the direction for nonlinear conjugate gradient (NCG) optimization.
     * 
     */
    void reset_ncg_direction();

    // energy and force
    /**
     *
     * @brief Performs a linear search to find the optimal step size that minimizes energy using
     * either nonlinear conjugate gradient method or simple line search.
     *
     * @return double The optimal step size. Returns -1 if no efficient step size is found with
     * the simple line search method.
     */
    double linear_search_for_stepsize_to_minimize_energy();

    /**
     * @brief Use simulated annealing to sample more point beyond locally gradient descent
     * 
     * @param temperature 
     */
    void simulated_annealing_next_step();

    /**
     * @brief Update vertex coordinates using a non-linear conjugate gradient method.
     *
     */
    void update_vertex_using_NCG();
};