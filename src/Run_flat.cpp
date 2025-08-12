#include "Run_simulation.hpp"

void run_flat(std::string param_filename)
{
    Param inputParam;
    import_param_file(inputParam, param_filename + ".params");
    Mesh mesh(inputParam);
    mesh.setup_flat();

    //testing
    /*
    int nFaceX = mesh.param.nFaceX;
    int nFaceY = mesh.param.nFaceY;

    for (int i = 0; i <= nFaceX; i++){
        for (int j = 0; j <= nFaceY; j++){
            if (!(mesh.vertices[(nFaceX+1)*j + i].isGhost)){
                
                double zOrig = mesh.vertices[(nFaceX+1)*j + i].coord.get(2,0);
                mesh.vertices[(nFaceX+1)*j + i].coord.set(2, 0, zOrig + i * j * (nFaceX - i) * (nFaceY - j) / 1000.0);
            }
        }
    }
    */
    // Insertion mode
    if (mesh.param.isInsertionIncluded)
    {   
        vector<vector<int>> insertionPatch; //{{1418, 1419, 1420, 1421, 1422, 1423, 1424, 1425, 1426,
                                        // 1484, 1485, 1486, 1487, 1488, 1489, 1490, 1491, 1492,
                                        //  1550, 1551, 1552, 1553, 1554, 1555, 1556, 1557, 1558}};
        mesh.set_insertion_patch(insertionPatch);
    }

    mesh.calculate_element_area_volume(); // Calculate the elemental area and volume per triangles (faces)
    double init_area = 0.0; // void value, placeholder if mesh.param.area0 should not be changed
    if (mesh.param.setRelaxAreaToDefault) // Set area to default value if flag set to false
    {
        mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0); // Calculate initial area and volume
    }
    else
    {
        mesh.sum_membrane_area_and_volume(init_area, mesh.param.vol0); // Calculate initial area and volume
    }
    
    std::cout << "[Area and Volume]: " << mesh.param.area0 << " ; " << mesh.param.vol0 << std::endl;

    // Scaffolding mode
    if (mesh.param.isEnergyHarmonicBondIncluded)
    { 
        // read scaffolding file
        mesh.param.scaffoldingPoints = import_scaffolding_mesh(mesh.param.scaffoldingFileName);
        // move spline points upwards until the lower boundary is approximately z=0
        mesh.move_vertices_based_on_scaffolding(false);
        // find closest vertex point and save in vector
        mesh.set_scaffolding_vertices_correspondence();
        // move gag up -test
        //for (int i = 0; i < mesh.param.scaffoldingPoints.size(); i++){
        //    mesh.param.scaffoldingPoints[i].set(2, 0, mesh.param.scaffoldingPoints[i](2,0));
        //}
    }

    // Output the vertices and faces matrix
    mesh.write_faces_csv("face.csv");
    mesh.write_vertices_csv("vertex_begin.csv");
    
    // Initialize all value before minimum energy search 
    mesh.calculate_element_area_volume(); // Calculate the elemental area and volume per triangles (faces)
    mesh.sum_membrane_area_and_volume(init_area, mesh.param.vol0); // Calculate initial area and volume
    std::cout << "[Area and Volume]: " << mesh.param.area0 << " ; " << mesh.param.vol0 << std::endl;

    mesh.update_previous_coord_for_vertex(); // Update previous coordinate...
    mesh.update_reference_coord_from_previous_coord(); // ...and reference coord
    mesh.Compute_Energy_And_Force();// Calculate energy on faces and force on vertices
    mesh.update_previous_coord_for_vertex(); // Update previous coordinate
    mesh.update_previous_force_for_vertex(); // Update forces

    // Create a "Record" object to bookkeep the energies and forces throughout the search
    Record record(mesh.param.maxIterations); ///< Bookkeeping for iterations
    record.add(mesh.param.area, mesh.param.energy, mesh.calculate_mean_force()); ///< Record the data before search (N=0)

    // Create a "Model" container object that packs mesh and record
    // This object also contains a member model.oa (OptimizationAlgorithm) which defines
    // the nonlinear conjugate method (NCG) that we use for efficient lowest energy search
    Model model(mesh, record);
 
    while (model.should_continue_optimization())
    {
        // The step size a0 needs a trial value, which is determined by rule-of-thumb.
        model.determine_trial_step_size(); 
        // Search for a step size above model.oa.stepSizeThreshold
        model.linear_search_for_stepsize_to_minimize_energy(); 
        // Simulated annealing
        model.simulated_annealing_next_step();
        if (model.iteration > 1000) model.highTemperature = 0.0;
        // Command line out put for current step
        std::cout << "[main()] " << model.to_string_current_step() << std::endl;
 
        // End program if no efficient step size found; assumed energy minimized (criteria based on model.oa.stepSizeThreshold)
        if (model.stepSize < 0)
        {
            std::cout << "[main()] Step: " << model.iteration << "  -- could not find an efficient step size. Stopped." << std::endl;
            break;
        }
 
        model.update_vertex_using_NCG(); // apply displacement to vertices based on forces

        // Update spring constant if scaffolding is included
        if (mesh.param.isEnergyHarmonicBondIncluded &&
            model.iteration % mesh.param.springConstScalingInterv == 0 &&
            mesh.param.springConst <= mesh.param.springConstUpperBound)
        {
            mesh.param.springConst *= mesh.param.springConstScaler;
            model.reset_ncg_direction();
            std::cout << "[main()] Increase harmonic bond spring constants to : " << mesh.param.springConst << std::endl;
        }

        // Propagate scaffolding
        if (model.iteration % mesh.param.propagateScaffoldingInterv == 0 &&
            model.iteration >= 1000)
        {
            mesh.propagate_scaffolding();
            // reset closest vertex point and save in vector
            mesh.set_scaffolding_vertices_correspondence();
            model.reset_ncg_direction();
            std::cout << "[main()] Reset scaffolding vertices correspondence." << mesh.param.springConst << std::endl;
        }
 
        // calculate the new Force and Energy; sync previous energy and force values with current values
        mesh.Compute_Energy_And_Force();
        mesh.update_previous_coord_for_vertex();
        mesh.update_previous_force_for_vertex();
        mesh.update_previous_energy_for_face();
        mesh.param.energyPrev = mesh.param.energy;

        // Calculate the direction vector for nonlinear conjugate method (NCG)
        model.update_ncg_direction();
        model.oa.trialStepSize = model.stepSize;
 
        // Record the Energy and nodal Force
        record.add(mesh.param.area, mesh.param.energy, mesh.calculate_mean_force());
        
        // Update the reference if energy or force is not changing signficantly
        if (abs(record.energyVec[model.iteration].energyTotal - record.energyVec[model.iteration - 1].energyTotal) < 1e-3 ||
            abs(record.meanForce[model.iteration] - record.meanForce[model.iteration - 1]) < 1e-3)
        {
            mesh.update_reference_coord_from_previous_coord();
            std::cout << "[main()] Small energy or force difference. Updated reference structure." << std::endl;
        }

        // If NCG has been stucked for too many times consecutively, change the flag to not-using-NCG
        model.oa.disable_NCG_if_stuck_consecutively();

        // Output parameters
        record.print_iteration(model.iteration);

        // Output vertices every 100 iterations
        //if (mesh.param.meshpointOutput) {
        //    output_trajectory_files(mesh, "traj");
        //}
        if (model.iteration % 10000 == 0)
        {
            write_vertex_data_to_csv(model.mesh, model.iteration);
        }
        model.iteration++;
    }
    // Output Energy and meanforce
    write_energy_force_data_to_csv(model);

    // Output element face energy
    // write_element_face_energy_to_csv(model);
    std::cout << "====element face energy====" << std::endl;
    std::cout << "face index, eb, ea, et" << std::endl;
    for (Face& face : model.mesh.faces){
        std::cout << face.index << "," << face.energy.energyCurvature << "," << face.energy.energyArea << "," << face.energy.energyTotal << std::endl;
    }

    // Output the final structure
    write_final_vertex_data_to_csv(mesh);
}
