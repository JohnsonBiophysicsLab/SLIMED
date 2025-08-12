#include "model/Model.hpp"

/**
 * @brief Constructs a new Model object.
 *
 * @param mesh_ The mesh object to be encapsulated.
 * @param record_ The record object to be encapsulated.
 */
Model::Model(Mesh &mesh_, Record &record_)
    : mesh(mesh_),
      record(record_),
      iteration(0),
      ncgDirection0(std::vector<Force>(mesh.vertices.size()))
{
    std::transform(mesh.vertices.begin(), mesh.vertices.end(),
                   ncgDirection0.begin(), [](Vertex &v)
                   { return v.force; });
    // Output if model setup if verbose mode
    if (mesh.param.VERBOSE_MODE) {
        std::cout << "==================================================" << std::endl;
        std::cout << "[Model::Model()] Model set up with mesh and record" << std::endl;
        std::cout << "==================================================" << std::endl;
    }
    oa.energyDiffThreshold = mesh.param.deltaEnergyConverge;
    oa.forceDiffThreshold = mesh.param.deltaForceScaleConverge;
}

/**
 * @brief This code determines the step size used in an optimization algorithm.
 * The step size is used to calculate how far to move in each iteration of the
 * optimization algorithm to minimize an energy function.
 *
 * @return double trialStepSize
 */
void Model::determine_trial_step_size()
{
    // The step size a0 needs to be determined by rule-of-thumb.
    // Compute trial step size every 'trialIterationInterval' iterations or at the start of optimization (iteration = 0).
    if (iteration % oa.trialIterationInterval == 0)
    {
        double maxForceMag = mesh.get_max_force_magnitude(); ///< Maximum force magnitude acting on vertices in the mesh
        // test
        std::cout << "max_force_scale=" << maxForceMag << ", energy=" << mesh.param.energy.energyTotal << std::endl;
        // exit(0);
        if (maxForceMag == 0)
        {
            oa.trialStepSize = -1;
        }
        else if (maxForceMag < 5.0)
        {
            oa.trialStepSize = 0.1 * mesh.param.lFace / 5.0;
        }
        else
        {
            oa.trialStepSize = 0.1 * mesh.param.lFace / maxForceMag; // a0 = 0.1; try changing a0 to make LGD more stable
        }
        // increase trialIterationInterval
        oa.trialIterationInterval += 0;
    }
    else
    {
        // If not computing trial step size, increase current step size exponentially.
        oa.trialStepSize *= 1.1; // previous 2.0, it was far too big
        // NOTE: value is practical. Maybe ML can help!!
    }
    std::cout << "[Model::determine_trial_step_size] Trial step size = " << oa.trialStepSize << std::endl;
}

/**
 * @brief Returns a string representation of the current step's iteration number, trial step size, and optimal step size.
 * The function concatenates the values of model.iteration, model.oa.trialStepSize, and
 * model.stepSize to the following format: "step: , trial StepSize = , StepSize = ".
 * @return std::string A string representation of the current step's information.
 */
std::string Model::to_string_current_step()
{
    return "step: " + std::to_string(iteration) +
           ", trial StepSize = " + std::to_string(oa.trialStepSize) +
           ", StepSize = " + std::to_string(stepSize);
}

/**
 * @brief Updates the direction for nonlinear conjugate gradient (NCG) optimization.
 *
 * This function calculates the dot product of the previous and current forces to determine the NCG factor. It then updates
 * the NCG direction either as a combination of the current force and the previous NCG direction, or just the current force
 * depending on whether NCG optimization is being used. The NCG direction is updated in the input vector `ncgDirections`.
 *
 * @param[in, out] ncgDirections A vector of NCG directions to be updated.
 */
void Model::update_ncg_direction()
{
    vector<Force> &ncgDirection0 = this->ncgDirection0;
    double ncgFactorPrev = 0.0; // Calculated from previous force
    double ncgFactorCurr = 0.0; // Calculated from current force
    Matrix forceTmp(3, 1);      // Temporary force matrix to calculate dot product

#pragma omp parallel for reduction(+ \
                                   : ncgFactorPrev) private(forceTmp)
    for (Vertex &vertex : mesh.vertices)
    {
        // update ncgFactorPrev with previous force
        forceTmp = vertex.forcePrev.forceTotal;
        ncgFactorPrev += dot_col(forceTmp, forceTmp);
        // update ncgFactorCurr with current force
        forceTmp = vertex.force.forceTotal;
        ncgFactorCurr += dot_col(forceTmp, forceTmp);
    }

    // Calculate new direction for NCG using current force and previous NCG
    // direction scaled by peta1 = ncgFactorCurr / ncgFactorPrev
    if (!oa.usingNCG)
    { // If not using NCG
#pragma omp parallel for
        for (int i = 0; i < mesh.vertices.size(); i++)
        {
            // Update NCG direction to be the current force
            ncgDirection0[i].forceTotal = mesh.vertices[i].force.forceTotal;
        }
    }
    else
    { // If using NCG
        double peta1 = ncgFactorCurr / ncgFactorPrev;

#pragma omp parallel for
        for (int i = 0; i < mesh.vertices.size(); i++)
        {
            // Update NCG direction as a combination of current force and previous NCG direction
            ncgDirection0[i].forceTotal = mesh.vertices[i].force.forceTotal + peta1 * ncgDirection0[i].forceTotal;
        }
    }
}

/**
 * @brief Reset the direction for nonlinear conjugate gradient (NCG) optimization.
 * 
 */
void Model::reset_ncg_direction()
{
#pragma omp parallel for
    for(int i = 0; i < mesh.vertices.size(); i++)
    {
        ncgDirection0[i].forceTotal.set(0, 0, 0.0);
        ncgDirection0[i].forceTotal.set(1, 0, 0.0);
        ncgDirection0[i].forceTotal.set(2, 0, 0.0);
    }
}
