#include <algorithm>
#include <cmath>

#include "model/Model.hpp"

namespace
{
/**
 * @brief The minimizer of the cubic that interpolates phi and phi' at a and b
 * (Nocedal & Wright, eq. 3.59), kept at least a tenth of the interval away
 * from either end; the midpoint when the cubic has no real minimizer.
 */
double interpolate_step(double a, double phiA, double dphiA,
                        double b, double phiB, double dphiB)
{
    const double lower = std::min(a, b);
    const double upper = std::max(a, b);
    const double margin = 0.1 * (upper - lower);
    double step = 0.5 * (a + b);
    const double d1 = dphiA + dphiB - 3.0 * (phiA - phiB) / (a - b);
    const double radicand = d1 * d1 - dphiA * dphiB;
    if (radicand >= 0.0)
    {
        const double d2 = (b > a ? 1.0 : -1.0) * std::sqrt(radicand);
        const double denominator = dphiB - dphiA + 2.0 * d2;
        if (denominator != 0.0)
        {
            const double cubic = b - (b - a) * (dphiB + d2 - d1) / denominator;
            if (std::isfinite(cubic))
            {
                step = cubic;
            }
        }
    }
    return std::min(std::max(step, lower + margin), upper - margin);
}
} // namespace

/**
 * @brief Strong-Wolfe line search along ncgDirection0 (Nocedal & Wright,
 * Algorithms 3.5 and 3.6).
 *
 * phi(a) is the energy at coordPrev + a * ncgDirection0 and phi'(a) is
 * -sum F(a) . ncgDirection0, so every trial costs one energy and force
 * evaluation and gives both. Starting from oa.trialStepSize the step doubles
 * while it still satisfies sufficient decrease and the slope is still steeply
 * negative; once a trial overshoots, the minimizer lies in a bracket that each
 * further trial narrows by cubic interpolation. Pure backtracking cannot do
 * this: when the trial step is too long and the next shorter step already too
 * short, it skips the strong-Wolfe band and runs down to oa.stepThreshold.
 *
 * @param currentEnergy phi(0).
 * @param slope0 phi'(0); negative along a descent direction.
 * @return The accepted step size; the mesh is left at that step with its
 * energy and force. When no trial meets both conditions within
 * oa.nMaxNcgLineSearchTrials, the lowest-energy trial with sufficient decrease is
 * returned, if it gained at least a tenth of the decrease the slope predicts,
 * and oa.isNCGstuck is set so that the next direction restarts from the
 * force; otherwise -1.
 */
double Model::strong_wolfe_line_search(double currentEnergy, double slope0)
{
    const int nVertices = mesh.vertices.size();
    if (!(slope0 < 0.0))
    {
        return -1;
    }

    // lo: the trial with sufficient decrease and the lowest energy so far
    // (step 0 to begin with). hi: once set, the other end of an interval
    // that contains strong-Wolfe steps.
    double stepLo = 0.0, phiLo = currentEnergy, dphiLo = slope0;
    double stepHi = 0.0, phiHi = 0.0, dphiHi = 0.0;
    bool bracketed = false;
    double trial = oa.trialStepSize;

    for (int nTrial = 0; nTrial < oa.nMaxNcgLineSearchTrials; nTrial++)
    {
        this->stepSize = trial;
        update_vertex_using_NCG();
        mesh.Compute_Energy_And_Force();
        const double phi = mesh.param.energy.energyTotal;
        // The directional derivative at the trial point, so from the force
        // just computed there
        double dphi = 0.0;
        #pragma omp parallel for reduction(+ : dphi)
        for (int i = 0; i < nVertices; i++)
        {
            dphi -= dot_col(mesh.vertices[i].force.forceTotal, ncgDirection0[i].forceTotal);
        }

        // Written so that a non-finite energy counts as no decrease
        const bool sufficientDecrease = phi <= currentEnergy + oa.c1 * trial * slope0;
        if (!sufficientDecrease || !(phi < phiLo))
        {
            stepHi = trial, phiHi = phi, dphiHi = dphi;
            bracketed = true;
        }
        else
        {
            if (std::abs(dphi) <= oa.c2 * std::abs(slope0))
            {
                return this->stepSize;
            }
            if (bracketed ? dphi * (stepHi - stepLo) >= 0.0 : dphi >= 0.0)
            {
                stepHi = stepLo, phiHi = phiLo, dphiHi = dphiLo;
                bracketed = true;
            }
            stepLo = trial, phiLo = phi, dphiLo = dphi;
        }

        if (bracketed)
        {
            // A bracket this narrow has closed on a kink of phi (the
            // regularization switching case, a scaffold bond crossing its
            // point), where the slope jumps and no strong-Wolfe step exists.
            if (std::abs(stepHi - stepLo) <= 1e-3 * std::max(stepLo, stepHi))
            {
                break;
            }
            trial = interpolate_step(stepLo, phiLo, dphiLo, stepHi, phiHi, dphiHi);
        }
        else
        {
            trial *= 2.0;
        }
    }

    // Sufficient decrease without the curvature condition: take the step if
    // it gained at least a tenth of what the slope predicts. With c1 = 1e-6 a
    // round-off decrease at a step of 1e-14 passes sufficient decrease, and
    // taking it would make the next trial step just as small. The step does
    // not keep the next Fletcher-Reeves direction downhill, so the direction
    // restarts from the force.
    if (stepLo <= 0.0 || !(phiLo <= currentEnergy + 0.1 * stepLo * slope0))
    {
        return -1;
    }
    oa.isNCGstuck = true;
    if (this->stepSize != stepLo)
    {
        this->stepSize = stepLo;
        update_vertex_using_NCG();
        mesh.Compute_Energy_And_Force();
    }
    return this->stepSize;
}

/**
 *
 * @brief Performs a linear search to find the optimal step size that minimizes energy using
 * either nonlinear conjugate gradient method or simple line search.
 *
 * With NCG on, strong_wolfe_line_search() looks along ncgDirection0. When it
 * finds no step that decreases the energy, the search falls back to the
 * simple line search along the force. The simple line search tries
 * oa.trialStepSize first and shrinks the step by 0.8 after each trial that
 * does not lower the energy. When it reaches oa.stepThreshold it refreshes
 * the regularization reference and searches once more along the force there,
 * from the rule-of-thumb step (the fastest vertex moves 0.1 lFace). So a
 * failure means that no step from the rule of thumb down to
 * oa.stepThreshold lowers the energy -- not that the step decayed over
 * earlier iterations.
 *
 * @param[in,out] ncgDirection0 Reference to a vector of nodal force directions.
 * Will be updated if switching to simple line search.
 * @return double The optimal step size, also stored in stepSize. Returns -1, and sets
 * stepSize to -1, if no efficient step size is found with the simple line search method.
 */
double Model::linear_search_for_stepsize_to_minimize_energy()
{
    // Initialize variables at the beginning of the function
    const int nVertices = mesh.vertices.size();
    double currentEnergy = mesh.param.energy.energyTotal;
    double ncgFactor0 = 0.0;
    // Compute the dot product of the negative force and the element triangle area in parallel
    #pragma omp parallel for reduction(+ : ncgFactor0)
    for (int i = 0; i < nVertices; i++)
    {
        ncgFactor0 -= dot_col(mesh.vertices[i].forcePrev.forceTotal, ncgDirection0[i].forceTotal);
    }
    // Reset flags and variables for NCG optimization
    oa.reset_NCG_stuck();

    if (oa.usingNCG)
    {
        if (strong_wolfe_line_search(currentEnergy, ncgFactor0) > 0.0)
        {
            return this->stepSize;
        }
        cout << "Now change the NCG WolfeConditions to simple line search method!" << endl;
        // Refreshing the reference is the whole fallback: every face
        // in the reference-length case then starts the retry with no
        // regularization energy and no regularization force. The R-pi
        // adaptive cases (Param::usingRpi) stay as configured. Turning
        // them off here and back on for the next search would make
        // that search compare the energy it starts from and the
        // energies it evaluates under two different functions, and
        // where any face is in an adaptive case no step passes that.
        mesh.update_reference_coord_from_previous_coord();
        oa.isNCGstuck = true;
        // Change ncgDirection0 to nodal force in parallel
        #pragma omp parallel for
        for (int i = 0; i < nVertices; i++)
        {
            ncgDirection0[i].forceTotal = mesh.vertices[i].forcePrev.forceTotal;
        }
    }

    // Simple line search when NCG is not used or found no step
    this->stepSize = oa.trialStepSize;
    bool retried = false;
    while (true)
    {
        // Apply displacement and recompute energy and force
        update_vertex_using_NCG();
        mesh.Compute_Energy_And_Force();
        if (mesh.param.energy.energyTotal < currentEnergy)
        {
            break;
        }
        // Scale down step size
        this->stepSize *= 0.8;
        if (this->stepSize < oa.stepThreshold)
        {
            if (retried)
            {
                cout << "Note: cannot find an efficient small stepSize to minimize the Energy even with the simple method!" << endl;
                this->stepSize = -1;
                return -1;
            }
            // Before giving up, search once more from a fresh start. The
            // trial step follows the last accepted step and can have become
            // tiny, so start from the rule of thumb determine_trial_step_size
            // uses: the fastest vertex moves a tenth of lFace. And the force
            // need not be a descent direction of the energy with the
            // regularization reference as it is (this is what the NCG
            // fallback above refreshes it for), so refresh it, re-evaluate
            // the start point under it and search along that force.
            retried = true;
            mesh.update_reference_coord_from_previous_coord();
            this->stepSize = 0.0;
            update_vertex_using_NCG();
            mesh.Compute_Energy_And_Force();
            mesh.update_previous_force_for_vertex();
            currentEnergy = mesh.param.energy.energyTotal;
            double maxForceMag = 0.0;
            #pragma omp parallel for reduction(max : maxForceMag)
            for (int i = 0; i < nVertices; i++)
            {
                ncgDirection0[i].forceTotal = mesh.vertices[i].forcePrev.forceTotal;
                maxForceMag = std::max(maxForceMag, ncgDirection0[i].forceTotal.calculate_norm());
            }
            this->stepSize = 0.1 * mesh.param.lFace / std::max(maxForceMag, 5.0);
        }
    }
    // Update member stepSize
    return this->stepSize;
}

/**
 * @brief Update vertex coordinates using a non-linear conjugate gradient method.
 *
 * @param stepSize The step size for the movement of each vertex
 * @param ncgDirection0 A vector of forces that determine the direction and magnitude of 
 *      movement for each vertex
 */
void Model::update_vertex_using_NCG()
{
    vector<Force>& ncgDirection0 = this->ncgDirection0;
    double stepSize = this->stepSize;

    // Update the position coordinates of vertices in parallel
    // by step size and direction vector (force)
#pragma omp parallel for
    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        // Update coordinate of current vertex using previous coordinate and step size in a specific direction
        mesh.vertices[i].coord = mesh.vertices[i].coordPrev + stepSize * ncgDirection0[i].forceTotal;
    }

    enforce_boundary_conditions_after_coordinate_update();
}

void Model::enforce_boundary_conditions_after_coordinate_update()
{
    int nFaceX = mesh.param.nFaceX; ///< Number of faces in the x-axis division of the mesh.
    int nFaceY = mesh.param.nFaceY; ///< Number of faces in the y-axis division of the mesh.

    // Deal with the boundary/ghost vertex based on boundary condition type defined as
    // enum class BoundaryType
    switch (mesh.param.boundaryCondition)
    {
    // Fixed boundary condition
    case BoundaryType::Fixed:
#pragma omp parallel for
        for (Face &face : mesh.faces)
        {
            // If face is not a boundary, skip to the next iteration
            if (!face.isBoundary)
                continue;

            // Get adjacent vertices of current face
            int node1 = face.adjacentVertices[0];
            int node2 = face.adjacentVertices[1];
            int node3 = face.adjacentVertices[2];

            // Update coordinates of adjacent vertices with previous coordinates
            mesh.vertices[node1].update_coord_with_prev_coord();
            mesh.vertices[node2].update_coord_with_prev_coord();
            mesh.vertices[node3].update_coord_with_prev_coord();
        }
        break;
    // Periodic boundary condition
    case BoundaryType::Periodic:
#pragma omp parallel for
        for (int i = 0; i < nFaceX + 1; i++)
        {
            // sync vertically for periodic boundary condition
            int index0 = (nFaceX + 1) * 0 + i; // ghost points
            int index1 = (nFaceX + 1) * 1 + i;
            int index2 = (nFaceX + 1) * 2 + i;
            int index00 = (nFaceX + 1) * (nFaceY - 6) + i; // corresponding real points
            int index11 = (nFaceX + 1) * (nFaceY - 5) + i;
            int index22 = (nFaceX + 1) * (nFaceY - 4) + i;
            mesh.vertices[index0].coord = mesh.vertices[index0].coordPrev + (mesh.vertices[index00].coord - mesh.vertices[index00].coordPrev);
            mesh.vertices[index1].coord = mesh.vertices[index1].coordPrev + (mesh.vertices[index11].coord - mesh.vertices[index11].coordPrev);
            mesh.vertices[index2].coord = mesh.vertices[index2].coordPrev + (mesh.vertices[index22].coord - mesh.vertices[index22].coordPrev);
            index0 = (nFaceX + 1) * (nFaceY - 2) + i;
            index1 = (nFaceX + 1) * (nFaceY - 1) + i;
            index2 = (nFaceX + 1) * (nFaceY - 0) + i;
            index00 = (nFaceX + 1) * 4 + i;
            index11 = (nFaceX + 1) * 5 + i;
            index22 = (nFaceX + 1) * 6 + i;
            mesh.vertices[index0].coord = mesh.vertices[index0].coordPrev + (mesh.vertices[index00].coord - mesh.vertices[index00].coordPrev);
            mesh.vertices[index1].coord = mesh.vertices[index1].coordPrev + (mesh.vertices[index11].coord - mesh.vertices[index11].coordPrev);
            mesh.vertices[index2].coord = mesh.vertices[index2].coordPrev + (mesh.vertices[index22].coord - mesh.vertices[index22].coordPrev);
        }
#pragma omp parallel for
        for (int j = 0; j < nFaceY + 1; j++)
        {
            // sync horizontally for periodic boundary condition
            int index0 = (nFaceX + 1) * j + 0;
            int index1 = (nFaceX + 1) * j + 1;
            int index2 = (nFaceX + 1) * j + 2;
            int index00 = (nFaceX + 1) * j + nFaceX - 6;
            int index11 = (nFaceX + 1) * j + nFaceX - 5;
            int index22 = (nFaceX + 1) * j + nFaceX - 4;
            mesh.vertices[index0].coord = mesh.vertices[index0].coordPrev + (mesh.vertices[index00].coord - mesh.vertices[index00].coordPrev);
            mesh.vertices[index1].coord = mesh.vertices[index1].coordPrev + (mesh.vertices[index11].coord - mesh.vertices[index11].coordPrev);
            mesh.vertices[index2].coord = mesh.vertices[index2].coordPrev + (mesh.vertices[index22].coord - mesh.vertices[index22].coordPrev);
            index0 = (nFaceX + 1) * j + nFaceX;
            index1 = (nFaceX + 1) * j + nFaceX - 1;
            index2 = (nFaceX + 1) * j + nFaceX - 2;
            index00 = (nFaceX + 1) * j + 6;
            index11 = (nFaceX + 1) * j + 5;
            index22 = (nFaceX + 1) * j + 4;
            mesh.vertices[index0].coord = mesh.vertices[index0].coordPrev + (mesh.vertices[index00].coord - mesh.vertices[index00].coordPrev);
            mesh.vertices[index1].coord = mesh.vertices[index1].coordPrev + (mesh.vertices[index11].coord - mesh.vertices[index11].coordPrev);
            mesh.vertices[index2].coord = mesh.vertices[index2].coordPrev + (mesh.vertices[index22].coord - mesh.vertices[index22].coordPrev);
        }
        break;

    // Per-vertex boundary: clamped vertices go back where they were, and every
    // periodic image is put at its source's new position plus its offset. The
    // images carried no force into the step -- it was folded onto their
    // sources -- so this is the only thing that ever moves them.
    case BoundaryType::Mixed:
        for (Vertex &vertex : mesh.vertices)
        {
            if ((vertex.is_fixed() || vertex.type == VertexType::Ghost) &&
                vertex.coordPrev.nrow() == 3)
            {
                vertex.update_coord_with_prev_coord();
            }
        }
        mesh.sync_periodic_images();
        break;

    // Free boundary condition
    case BoundaryType::Free:

        int index1, index2, index3, index4;
        // left side
        for (int j = 2; j < nFaceY; j++)
        {
            if (j & 1) // if odd
            {
                index1 = (nFaceX + 1) * j + 0;
                index2 = (nFaceX + 1) * j + 1;
                index3 = (nFaceX + 1) * (j - 1) + 0;
                index4 = (nFaceX + 1) * (j - 1) + 1;
                
            }
            else
            {
                index1 = (nFaceX + 1) * j + 0;
                index2 = (nFaceX + 1) * j + 1;
                index3 = (nFaceX + 1) * (j - 1) + 1;
                index4 = (nFaceX + 1) * (j - 1) + 2;
            }
            mesh.vertices[index1].coord = mesh.vertices[index2].coord + (mesh.vertices[index3].coord - mesh.vertices[index4].coord);
        }
        index1 = (nFaceX + 1) * 1 + 0;
        index2 = (nFaceX + 1) * 2 + 0;
        index3 = (nFaceX + 1) * 1 + 1;
        index4 = (nFaceX + 1) * 2 + 1;
        mesh.vertices[index1].coord = mesh.vertices[index2].coord + (mesh.vertices[index3].coord - mesh.vertices[index4].coord);
        // right side
        index1 = (nFaceX + 1) * 1 + nFaceX;
        index2 = (nFaceX + 1) * 2 + nFaceX - 1;
        index3 = (nFaceX + 1) * 1 + nFaceX - 1;
        index4 = (nFaceX + 1) * 2 + nFaceX - 2;
        mesh.vertices[index1].coord = mesh.vertices[index2].coord + (mesh.vertices[index3].coord - mesh.vertices[index4].coord);
        for (int j = 2; j < nFaceY; j++)
        {
            if (j & 1) // if odd
            {
                index1 = (nFaceX + 1) * j + nFaceX;
                index2 = (nFaceX + 1) * j + nFaceX - 1;
                index3 = (nFaceX + 1) * (j - 1) + nFaceX - 1;
                index4 = (nFaceX + 1) * (j - 1) + nFaceX - 2;
                
            }
            else
            {
                index1 = (nFaceX + 1) * j + nFaceX;
                index2 = (nFaceX + 1) * j + nFaceX - 1;
                index3 = (nFaceX + 1) * (j - 1) + nFaceX;
                index4 = (nFaceX + 1) * (j - 1) + nFaceX - 1;
            }
            mesh.vertices[index1].coord = mesh.vertices[index2].coord + (mesh.vertices[index3].coord - mesh.vertices[index4].coord);
        }
        // bottom
        for (int i = 0; i < nFaceX; i++)
        {
            index1 = (nFaceX + 1) * 0 + i;
            index2 = (nFaceX + 1) * 1 + i;
            index3 = (nFaceX + 1) * 1 + i + 1;
            index4 = (nFaceX + 1) * 2 + i;
            mesh.vertices[index1].coord = mesh.vertices[index2].coord + (mesh.vertices[index3].coord - mesh.vertices[index4].coord);
        }
        index1 = (nFaceX + 1) * 0 + nFaceX;
        index2 = (nFaceX + 1) * 0 + nFaceX - 1;
        index3 = (nFaceX + 1) * 1 + nFaceX;
        index4 = (nFaceX + 1) * 1 + nFaceX - 1;
        mesh.vertices[index1].coord = mesh.vertices[index2].coord + (mesh.vertices[index3].coord - mesh.vertices[index4].coord);
        // top
        for (int i = 1; i < nFaceX + 1; i++)
        {
            index1 = (nFaceX + 1) * nFaceY + i;
            index2 = (nFaceX + 1) * (nFaceY - 1) + i;
            index3 = (nFaceX + 1) * (nFaceY - 1) + i + 1;
            index4 = (nFaceX + 1) * (nFaceY - 2) + i;
            mesh.vertices[index1].coord = mesh.vertices[index2].coord + (mesh.vertices[index3].coord - mesh.vertices[index4].coord);
        }
        index1 = (nFaceX + 1) * nFaceY + nFaceX;
        index2 = (nFaceX + 1) * nFaceY + nFaceX - 1;
        index3 = (nFaceX + 1) * (nFaceY - 1) + nFaceX;
        index4 = (nFaceX + 1) * (nFaceY - 1) + nFaceX - 1;
        mesh.vertices[index1].coord = mesh.vertices[index2].coord + (mesh.vertices[index3].coord - mesh.vertices[index4].coord);
        break;
    }
}

bool Model::simulated_annealing_next_step(bool forceAttempt)
{
    Param& param = mesh.param;
    if (!param.thermalFluctuationEnabled)
    {
        return false;
    }
    if (!forceAttempt &&
        (param.thermalFluctuationInterval <= 0 ||
        iteration % param.thermalFluctuationInterval != 0)
    )
    {
        return false;
    }

    const double kBoltzmannPNnmPerK = 0.01380649;
    double temperatureKelvin = param.thermalFluctuationTemperatureKelvin *
                               std::pow(param.thermalFluctuationCoolingRate,
                                        thermalFluctuationAttemptCount);
    if (temperatureKelvin < param.thermalFluctuationMinTemperatureKelvin)
    {
        temperatureKelvin = param.thermalFluctuationMinTemperatureKelvin;
    }

    const double effectiveKBT = kBoltzmannPNnmPerK * temperatureKelvin;
    const double displacementStddev = param.thermalFluctuationStepScale * param.lFace;
    if (!std::isfinite(temperatureKelvin) || !std::isfinite(effectiveKBT) ||
        !std::isfinite(displacementStddev) ||
        effectiveKBT <= 0.0 || displacementStddev <= 0.0)
    {
        return false;
    }

    const double currentEnergy = param.energy.energyTotal;
    std::vector<Matrix> previousCoords(mesh.vertices.size(), Matrix(3, 1));
    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        previousCoords[i] = mesh.vertices[i].coord;
    }

    std::vector<int> mobileVertexIndices;
    mobileVertexIndices.reserve(mesh.vertices.size());
    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        // Ghosts, clamped vertices and periodic images are not coordinates
        // of the model; a trial move on one would be undone by the boundary
        // enforcement below and count as a rejection for nothing.
        if (mesh.is_independent_vertex(i))
        {
            mobileVertexIndices.push_back(i);
        }
    }
    if (mobileVertexIndices.empty())
    {
        return false;
    }

    std::uniform_int_distribution<int> vertexDistribution(0, mobileVertexIndices.size() - 1);
    Vertex& movedVertex = mesh.vertices[mobileVertexIndices[vertexDistribution(thermalRng)]];
    std::normal_distribution<double> normalDistribution(0.0, displacementStddev);
    const double maxDisplacement = 3.0 * displacementStddev;
    for (int j = 0; j < 3; j++)
    {
        double displacement = normalDistribution(thermalRng);
        if (displacement > maxDisplacement)
        {
            displacement = maxDisplacement;
        }
        else if (displacement < -maxDisplacement)
        {
            displacement = -maxDisplacement;
        }
        movedVertex.coord.set(j, 0, movedVertex.coord(j, 0) + displacement);
    }

    enforce_boundary_conditions_after_coordinate_update();
    mesh.Compute_Energy_And_Force();
    const double deltaEnergy = param.energy.energyTotal - currentEnergy;
    double acceptanceProbability = 1.0;
    bool accepted = true;
    if (!std::isfinite(deltaEnergy))
    {
        acceptanceProbability = 0.0;
        accepted = false;
    }
    else if (deltaEnergy > 0.0)
    {
        acceptanceProbability = std::exp(-deltaEnergy / effectiveKBT);
        std::uniform_real_distribution<double> uniformDistribution(0.0, 1.0);
        accepted = uniformDistribution(thermalRng) < acceptanceProbability;
    }

    if (!accepted)
    {
        for (int i = 0; i < mesh.vertices.size(); i++)
        {
            mesh.vertices[i].coord = previousCoords[i];
        }
        mesh.Compute_Energy_And_Force();
    }

    ThermalFluctuationRecord record;
    record.iteration = iteration;
    record.temperatureKelvin = temperatureKelvin;
    record.kBT = effectiveKBT;
    record.deltaEnergy = deltaEnergy;
    record.acceptanceProbability = acceptanceProbability;
    record.accepted = accepted;
    thermalFluctuationRecords.push_back(record);
    thermalFluctuationAttemptCount++;

    std::cout << "[Model::simulated_annealing_next_step] iteration=" << iteration
              << ", temperatureKelvin=" << temperatureKelvin
              << ", kBT=" << effectiveKBT
              << ", deltaEnergy=" << deltaEnergy
              << ", acceptanceProbability=" << acceptanceProbability
              << ", accepted=" << accepted << std::endl;

    return true;
}
