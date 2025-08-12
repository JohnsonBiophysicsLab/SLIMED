#include "Dynamics.hpp"

/**
 * @brief Constructs a new Model object.
 *
 * @param mesh_ The mesh object to be encapsulated.
 * @param record_ The record object to be encapsulated.
 */
DynamicModel::DynamicModel(DynamicMesh &mesh_, Record &record_)
    : Model(mesh_, record_),
    mesh(mesh_)
{
    // force and rand scale const calculated based on Brownian Dynamics
    forceScaleConst = mesh.param.diffConst * mesh.param.timeStep / mesh.param.KBT;
    randScaleConst = pow(2 * mesh.param.diffConst * mesh.param.timeStep, 0.5); //test override
    
    // commandline print the two scale constants
    if (mesh.param.VERBOSE_MODE) 
    {
        std::cout << "[DynamicModel] forceScaleConst = " << forceScaleConst << std::endl;
        std::cout << "[DynamicModel] randScaleConst = " << randScaleConst << std::endl;
    }
    
    randomSeed = mesh.param.randomSeed;
    std::mt19937 gen{randomSeed};
    std::normal_distribution<> normal_dist{0.0,1.0};
    // Output if model setup if verbose mode
    if (mesh.param.VERBOSE_MODE)
    {
        std::cout << "==================================================" << std::endl;
        std::cout << "[DynamicModel::DynamicModel()] Model set up with mesh and record" << std::endl;
        std::cout << "==================================================" << std::endl;
    }
}

/**
 * @brief This code iterates over all vertices in a mesh and displaces
 * them randomly based on a force term and a random term. The displacement
 * is added to the original vertex position.
 *
 * The boundary conditions are checked for each vertex before
 * displacement, and if the vertex is a ghost or boundary,
 * then the displacement is set to 0.
 *
 * The displacement amount, along with the vertex index,
 * whether it is a boundary or ghost vertex, and the dimension
 * being displaced (x, y, z) is printed to the console using std::cout.
 */
void DynamicModel::next_step()
{
    // Initialize values 

    for (int i = 0; i < mesh.vertices.size(); i++)
    {   
        
        // calculate unit normal vector
        mesh.vertices[i].approximate_unit_normal_vector(mesh.faces);
    }

#pragma omp parallel for shared(mesh, normal_dist, gen) schedule(static)
    for (int i = 0; i < mesh.vertices.size(); i++)
    {   
        Matrix forcetermVector(3, 1);
        Matrix randomtermVector(3, 1);

        if (mesh.vertices[i].isGhost)
        {
            double disp = 0.0; // No displacement for boundary or ghost vertices
        }
        else
        {
            for (int j = 0; j < 3; j++)
            {
                double disp = 0.0; ///< double to store displacement for one coord of one pt
                // boundary condition - periodic is postprocessed differently!
                double forceterm = 0.0;
                double randomterm = 0.0;

                    
                forceterm = mesh.vertices[i].force.forceCurvature(j, 0) +
                                    mesh.vertices[i].force.forceArea(j, 0); // Get force term
                //std::cout << "Force @ " << i << " , "<< j << " = "  << forceterm << std::endl;
                if (std::isnan(forceterm))
                {
                    forceterm = 0.0;
                }
                randomterm = normal_dist(gen); // Get random term
                // x, y direction are trivial for now
                // ! need to calculate the normal vector for the surface in the future
                if (j < 2) {
                    randomterm *= 0.0;
                    forceterm *= 0.0;
                }
                //randomterm *= 0.0;

                randomtermVector.set(j, 0, randomterm);
                forcetermVector.set(j, 0, forceterm);

                //std::cout << "vertex " << i << "@ " << j << ", scal = " << scal << 
                // ", rt =" << randomterm << ",ft = " <<forceterm<< std::endl;

                // forceterm = 0.0;
            
            }

            
            // Only perpendicular displacement is taken into account
            //double scale = dot_col(mesh.vertices[i].normVector, randomtermVector);

            //randomtermVector = scale * mesh.vertices[i].normVector;
            //scale = dot_col(mesh.vertices[i].normVector, forcetermVector);

            //forcetermVector = scale * mesh.vertices[i].normVector;

            //424 for 40 x 40 membrane
            /*if (i == -1){
                std::cout << this->iteration <<
                 " NV: " <<  mesh.vertices[i].normVector(0,0) << ", " <<
                mesh.vertices[i].normVector(1,0) << ", "<<
                mesh.vertices[i].normVector(2,0) << ", "<< "FT: "
                <<  forcetermVector(0,0) << ", "<<
                forcetermVector(1,0) << ", "<<
                forcetermVector(2,0) << ", "<< "RT: "
                <<  randomtermVector(0,0) << ", "<<
               randomtermVector(1,0) << ", "<<
                randomtermVector(2,0) << ", "<< std::endl;
            }*/

            for(int j = 0; j < 3; j++)
            {   
                // for testing only, muting x and y displacement
                double disp = 0.0;
                //if (j == 2)
                //if (i == 425 && j != 2){
                //    disp = 0.1;
                //}
                    disp = randScaleConst * randomtermVector(j,0)
                                + forceScaleConst * forcetermVector(j,0);
                //std::cout << "Disp @ " << i << " , "<< j << " = "  << disp << std::endl;
                //if (j != 2)
                double original = mesh.matSurface(i, j);
                mesh.matSurface.set(i, j, original + disp);

                if (mesh.param.VERBOSE_MODE) 
                {
                    std::cout << "[DynamicModel] Index: " << i <<
                                " , z-Original:" << original <<
                                " , z-Displacement: " << disp <<
                                ", fc: " << mesh.vertices[i].force.forceCurvature(j, 0) <<
                                ", fa: " << mesh.vertices[i].force.forceArea(j, 0) <<
                                ", z-rt: " << randomtermVector(j,0) <<
                                ", G: " << mesh.vertices[i].isGhost << std::endl;
                }
                
            }
        }
    }
    //rsc = 0.04472135955 for ts = 0.001
    //randScaleConst -= 0.0000002;
    //if (randScaleConst <= 0.04462135955){
    //    randScaleConst = 0.0;
    //}
}