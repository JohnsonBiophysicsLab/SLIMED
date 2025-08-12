#include "model/Model.hpp"

/**
 * 
 * @brief Performs a linear search to find the optimal step size that minimizes energy using
 * either nonlinear conjugate gradient method or simple line search.
 *
 * @param[in,out] ncgDirection0 Reference to a vector of nodal force directions.
 * Will be updated if switching to simple line search.
 * @return double The optimal step size. Returns -1 if no efficient step size is found with
 * the simple line search method.
 */
double Model::linear_search_for_stepsize_to_minimize_energy()
{
    bool ROWOUTPUT = false;
//if (ROWOUTPUT)cout << "R15" << endl;
    // Initialize variables at the beginning of the function
    const int nVertices = mesh.vertices.size();
    this->stepSize = oa.trialStepSize;
    double currentEnergy = mesh.param.energy.energyTotal;
    double ncgFactor0 = 0.0;
//if (ROWOUTPUT)cout << "R21" << endl;
    // Compute the dot product of the negative force and the element triangle area in parallel
    #pragma omp parallel for reduction(+ : ncgFactor0)
    for (int i = 0; i < nVertices; i++)
    {
        ncgFactor0 -= dot_col(mesh.vertices[i].forcePrev.forceTotal, ncgDirection0[i].forceTotal);
    }
//if (ROWOUTPUT)cout << "R28" << endl;
    // Set up the constants for the nonlinear conjugate gradient method
    double newEnergy = 0.0;
    double ncgFactor = 0.0;
//if (ROWOUTPUT)cout << "R32" << endl;
    // Reset flags and variables for NCG optimization
    oa.reset_NCG_Rpi();

    // Start the line search loop until the criteria is satisfied or the step size falls below a threshold
    while (true)
    {
// if (ROWOUTPUT)cout << "R38" << endl;
        // Scale down step size
        this->stepSize *= 0.8;
//if (ROWOUTPUT)cout << "R41" << endl;
        // Apply displacement and recompute energy and force
        update_vertex_using_NCG();//cout << "R43" << endl;
        mesh.Compute_Energy_And_Force();//cout << "R44" << endl;
        newEnergy = mesh.param.energy.energyTotal;
//if (ROWOUTPUT)cout << "R46" << endl;
        // Check if using NCG and not stuck
        if (oa.usingNCG && !oa.isNCGstuck)
        {
            // Calculate NCG factor in parallel
            #pragma omp parallel for reduction(+ : ncgFactor)
            for (int i = 0; i < nVertices; i++)
            {
                ncgFactor -= dot_col(mesh.vertices[i].forcePrev.forceTotal, ncgDirection0[i].forceTotal);
            }
//if (ROWOUTPUT)cout << "R56" << endl;
            // Satisfy NCG Wolfe conditions
            if (newEnergy <= currentEnergy + oa.c1 * this->stepSize * ncgFactor0 && abs(ncgFactor) <= oa.c2 * abs(ncgFactor0))
            {
                break;
            }
//if (ROWOUTPUT)cout << "R62" << endl;
            // Switch to simple linear search if cannot find a step size larger than oa.stepThreshold
            if (this->stepSize < oa.stepThreshold)
            {
                cout << "Now change the NCG WolfeConditions to simple line search method!" << endl;

                // Reinitialize variables for simple line search
                this->stepSize = oa.trialStepSize;
                mesh.update_reference_coord_from_previous_coord();
                oa.isNCGstuck = true;
                oa.usingRpi = false;
//if (ROWOUTPUT)cout << "R73" << endl;
                // Change ncgDirection0 to nodal force in parallel
                #pragma omp parallel for
                for (int i = 0; i < nVertices; i++)
                {
                    ncgDirection0[i].forceTotal = mesh.vertices[i].forcePrev.forceTotal;
                }
            }
        }
        else // Simple line search when NCG is not used or stuck
        {
            if (newEnergy < currentEnergy)
            {
                break;
            }

            if (this->stepSize < oa.stepThreshold)
            {
                cout << "Note: cannot find an efficient small stepSize to minimize the Energy even with the simple method!" << endl;
                return -1;
            }
        }
    }
//if (ROWOUTPUT)cout << "R96" << endl;
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
    int nFaceX = mesh.param.nFaceX; ///< Number of faces in the x-axis division of the mesh.
    int nFaceY = mesh.param.nFaceY; ///< Number of faces in the y-axis division of the mesh.

    // Update the position coordinates of vertices in parallel
    // by step size and direction vector (force)
#pragma omp parallel for
    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        // Update coordinate of current vertex using previous coordinate and step size in a specific direction
        mesh.vertices[i].coord = mesh.vertices[i].coordPrev + stepSize * ncgDirection0[i].forceTotal;
    }

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

void Model::simulated_annealing_next_step()
{   
    std::random_device rd;
    std::mt19937 rng(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);

    if (isHeating)
    {
        for (Vertex& vertex : mesh.vertices)
        {
            for (int j = 0; j < 3; j++)
            {
                //std::cout << vertex.coord(j, 0) << std::endl;
                double random_number = distribution(rng);
                if (random_number > 3.0)
                {
                    random_number = 3.0;
                } else if (random_number < -3.0)
                {
                    random_number = -3.0;
                }
                vertex.coord.set(j, 0, vertex.coord(j, 0) + highTemperature * random_number);
                //std::cout << "anneal to : " << vertex.coord(j, 0) << std::endl;
            }
        }
        // Check if heating is done
        currentHeatingStep++;
        // Change to cooling if heating is done
        if (currentHeatingStep >= heatingStep)
        {
            isHeating = false;
            currentHeatingStep = 0;
            currentCoolingStep = 0;
        }
    } else {
        // Check if cooling is done
        currentCoolingStep++;
        // Change to heating if cooling is done
        if (currentCoolingStep >= coolingStep)
        {
            isHeating = true;
            currentHeatingStep = 0;
            currentCoolingStep = 0;
        }
    }
}

