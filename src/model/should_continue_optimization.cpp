#include "model/Model.hpp"

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
bool Model::should_continue_optimization()
{

    oa.isCriteriaSatisfied = iteration > 500 &&
                             abs((record.energyVec[iteration].energyTotal - record.energyVec[iteration - 500].energyTotal) / 500) < oa.energyDiffThreshold &&
                             record.meanForce[iteration] < oa.forceDiffThreshold;

    return (!oa.isCriteriaSatisfied) && iteration < mesh.param.maxIterations - 1;
}