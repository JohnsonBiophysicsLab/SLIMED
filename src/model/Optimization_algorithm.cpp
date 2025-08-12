#include "model/Model.hpp"

/**
 * @brief Constructor for the OptimizationAlgorithm class.
 *
 * Initializes an instance of the OptimizationAlgorithm class with default settings.
 * Sets usingNCG to true, isNCGstuck to false, and usingRpi to true.
 */
OptimizationAlgorithm::OptimizationAlgorithm()
{
}

/**
 * @brief Constructor for the OptimizationAlgorithm class.
 *
 * Initializes an instance of the OptimizationAlgorithm class with default settings.
 * Customize iteration interval, c1, c2, and step threshold.
 * Sets usingNCG to true, isNCGstuck to false, and usingRpi to true.
 */
OptimizationAlgorithm::OptimizationAlgorithm(const int trialIterationInterval,
                                             const double c1,
                                             const double c2,
                                             const double stepThreshold) : trialIterationInterval(trialIterationInterval),
                                                                           c1(c1),
                                                                           c2(c2),
                                                                           stepThreshold(stepThreshold)
{
}

/** @brief Reset NCG and Rpi to default.
 *
 * Sets usingNCG to true, isNCGstuck to false, and usingRpi to true.
 */
void OptimizationAlgorithm::reset_NCG_Rpi()
{
    //usingNCG = true;
    isNCGstuck = false;
    usingRpi = true;
}

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
bool OptimizationAlgorithm::disable_NCG_if_stuck_consecutively()
{
    if (isNCGstuck) {
        nConsecutiveNcgStuck++;
        if (nConsecutiveNcgStuck > nConsecutiveNcgStuckThreshold){
            std::cout << "[OptimizationAlgorithm::disable_NCG_if_stuck_consecutively] Disabled NCG"
            << std::endl;
            usingNCG = false;
        } else {
            usingNCG = true;
        }
    }
    else{
        nConsecutiveNcgStuck = 0;
        //usingNCG = true;
    }
    return usingNCG;
}