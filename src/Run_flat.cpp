#include "Run_simulation.hpp"

void run_flat(std::string param_filename)
{
    Param inputParam;
    import_param_file(inputParam, param_filename + ".params");
    const bool restartRequested = !inputParam.restartInputFile.empty();
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
    mesh.report_volume_rebaseline(); // TEMPORARY -- see docs/volume_functional_split.md step 4

    // Scaffolding mode
    if (mesh.param.isEnergyHarmonicBondIncluded)
    { 
        // read scaffolding file
        mesh.param.scaffoldingPoints = import_scaffolding_mesh(mesh.param.scaffoldingFileName);
        if (mesh.param.isGagScaffoldingEnergyIncluded ||
            mesh.param.isIdealizedProteinLatticeEnergyIncluded)
        {
            mesh.orient_scaffolding_plane_to_membrane();
            mesh.pre_relax_gag_scaffolding();
        }
        // initialize the membrane shape from the aligned scaffold geometry
        mesh.move_vertices_based_on_scaffolding(false);
        // find closest vertex point and save in vector
        mesh.set_scaffolding_vertices_correspondence();
        // move gag up -test
        //for (int i = 0; i < mesh.param.scaffoldingPoints.size(); i++){
        //    mesh.param.scaffoldingPoints[i].set(2, 0, mesh.param.scaffoldingPoints[i](2,0));
        //}
    }

    // Output the vertices and faces matrix
    if (!restartRequested)
    {
        mesh.write_faces_csv("face.csv");
        mesh.write_vertices_csv("vertex_begin.csv");
    }
    
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

    if (restartRequested)
    {
        if (!load_model_restart_checkpoint(model, mesh.param.restartInputFile))
        {
            std::cerr << "[main()] Restart requested but checkpoint load failed. Stopped." << std::endl;
            return;
        }
        if (mesh.param.isEnergyHarmonicBondIncluded)
        {
            mesh.set_scaffolding_vertices_correspondence();
        }
        mesh.Compute_Energy_And_Force();
        mesh.update_previous_coord_for_vertex();
        mesh.update_previous_force_for_vertex();
        mesh.update_previous_energy_for_face();
        mesh.param.energyPrev = mesh.param.energy;
    }

    if (mesh.param.meshpointOutput)
    {
        model.mesh.write_vertices_csv("vertex_initial_step_0.csv");
        if (mesh.param.isGagScaffoldingEnergyIncluded ||
            mesh.param.isIdealizedProteinLatticeEnergyIncluded)
        {
            model.mesh.write_gag_scaffolding_state_dat("gag_scaffold_initial_step_0.dat");
        }
    }
 
    if (mesh.param.thermalFluctuationPureMMC)
    {
        if (!mesh.param.thermalFluctuationEnabled)
        {
            std::cerr << "[main()] thermalFluctuationPureMMC requires thermalFluctuationEnabled = true. Stopped." << std::endl;
            return;
        }
        std::cout << "[main()] Running pure Metropolis Monte Carlo thermal fluctuation mode. "
                  << "No NCG minimization steps will be performed; one MMC trial is attempted every step."
                  << std::endl;

        while (model.iteration < mesh.param.maxIterations - 1)
        {
            const bool thermalMoveAttempted = model.simulated_annealing_next_step(true);
            if (!thermalMoveAttempted)
            {
                std::cerr << "[main()] Pure MMC thermal move was not attempted. Check thermal parameters. Stopped." << std::endl;
                break;
            }

            mesh.update_previous_coord_for_vertex();
            mesh.update_previous_force_for_vertex();
            mesh.update_previous_energy_for_face();
            mesh.param.energyPrev = mesh.param.energy;

            const int recordIndex = static_cast<int>(record.energyVec.size());
            record.add(mesh.param.area, mesh.param.energy, mesh.calculate_mean_force());
            record.print_iteration(recordIndex);

            if (mesh.param.meshpointOutput && model.iteration < 5)
            {
                const int snapshotIteration = model.iteration + 1;
                model.mesh.write_vertices_csv("vertex_initial_step_" + std::to_string(snapshotIteration) + ".csv");
                if (mesh.param.isGagScaffoldingEnergyIncluded ||
                    mesh.param.isIdealizedProteinLatticeEnergyIncluded)
                {
                    model.mesh.write_gag_scaffolding_state_dat("gag_scaffold_initial_step_" + std::to_string(snapshotIteration) + ".dat");
                }
            }

            const int nextIteration = model.iteration + 1;
            if (nextIteration % 10000 == 0)
            {
                write_vertex_data_to_csv(model.mesh, nextIteration);
                if (mesh.param.isGagScaffoldingEnergyIncluded ||
                    mesh.param.isIdealizedProteinLatticeEnergyIncluded)
                {
                    mesh.write_gag_scaffolding_state_dat("gag_scaffold_" + std::to_string(nextIteration) + ".dat");
                }
            }
            if (mesh.param.checkpointOutputInterval > 0 &&
                nextIteration % mesh.param.checkpointOutputInterval == 0)
            {
                write_model_restart_checkpoint(model,
                                               mesh.param.checkpointOutputFile,
                                               nextIteration);
            }
            model.iteration++;
        }
    }
    else
    {
    while (model.should_continue_optimization())
    {
        // The step size a0 needs a trial value, which is determined by rule-of-thumb.
        model.determine_trial_step_size(); 
        if (model.oa.trialStepSize < 0 && mesh.param.thermalFluctuationEnabled)
        {
            bool thermalMoveAttempted = model.simulated_annealing_next_step();
            if (thermalMoveAttempted)
            {
                if (!model.thermalFluctuationRecords.back().accepted)
                {
                    std::cout << "[main()] Thermal fluctuation trial was rejected and no deterministic force direction was available. Stopped." << std::endl;
                    break;
                }
                mesh.update_previous_coord_for_vertex();
                mesh.update_previous_force_for_vertex();
                mesh.update_previous_energy_for_face();
                mesh.param.energyPrev = mesh.param.energy;
                model.reset_ncg_direction();
                model.update_ncg_direction();
                model.oa.trialStepSize = 0.1 * mesh.param.lFace / 5.0;
                record.add(mesh.param.area, mesh.param.energy, mesh.calculate_mean_force());
                record.print_iteration(model.iteration);
                model.iteration++;
                continue;
            }
        }
        // Search for a step size above model.oa.stepSizeThreshold
        model.linear_search_for_stepsize_to_minimize_energy(); 
        // Command line out put for current step
        std::cout << "[main()] " << model.to_string_current_step() << std::endl;
 
        // End program if no efficient step size found; assumed energy minimized (criteria based on model.oa.stepSizeThreshold)
        if (model.stepSize < 0)
        {
            std::cout << "[main()] Step: " << model.iteration << "  -- could not find an efficient step size. Stopped." << std::endl;
            break;
        }
 
        model.update_vertex_using_NCG(); // apply displacement to vertices based on forces
        bool thermalMoveAttempted = model.simulated_annealing_next_step(); // optional Metropolis thermal trial
        if (thermalMoveAttempted)
        {
            model.reset_ncg_direction();
        }

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
        if (mesh.param.meshpointOutput && model.iteration < 5)
        {
            const int snapshotIteration = model.iteration + 1;
            model.mesh.write_vertices_csv("vertex_initial_step_" + std::to_string(snapshotIteration) + ".csv");
            if (mesh.param.isGagScaffoldingEnergyIncluded ||
                mesh.param.isIdealizedProteinLatticeEnergyIncluded)
            {
                model.mesh.write_gag_scaffolding_state_dat("gag_scaffold_initial_step_" + std::to_string(snapshotIteration) + ".dat");
            }
        }
        
        // Update the reference if energy or force is not changing signficantly
        if (model.iteration > 0 &&
            (abs(record.energyVec[model.iteration].energyTotal - record.energyVec[model.iteration - 1].energyTotal) < 1e-3 ||
             abs(record.meanForce[model.iteration] - record.meanForce[model.iteration - 1]) < 1e-3))
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
            if (mesh.param.isGagScaffoldingEnergyIncluded ||
                mesh.param.isIdealizedProteinLatticeEnergyIncluded)
            {
                mesh.write_gag_scaffolding_state_dat("gag_scaffold_" + std::to_string(model.iteration) + ".dat");
            }
        }
        const int nextIteration = model.iteration + 1;
        if (mesh.param.checkpointOutputInterval > 0 &&
            nextIteration % mesh.param.checkpointOutputInterval == 0)
        {
            write_model_restart_checkpoint(model,
                                           mesh.param.checkpointOutputFile,
                                           nextIteration);
        }
        model.iteration++;
    }
    }
    if (mesh.param.checkpointOutputInterval > 0)
    {
        write_model_restart_checkpoint(model,
                                       mesh.param.checkpointOutputFile,
                                       model.iteration);
    }
    // Output Energy and meanforce
    write_energy_force_data_to_csv(model);
    if (!model.thermalFluctuationRecords.empty())
    {
        std::ofstream thermalOutfile("ThermalFluctuation.csv");
        thermalOutfile << "iteration,temperature_kelvin,kBT_pN_nm,delta_energy,acceptance_probability,accepted\n";
        for (const auto& thermalRecord : model.thermalFluctuationRecords)
        {
            thermalOutfile << thermalRecord.iteration << ","
                           << thermalRecord.temperatureKelvin << ","
                           << thermalRecord.kBT << ","
                           << thermalRecord.deltaEnergy << ","
                           << thermalRecord.acceptanceProbability << ","
                           << thermalRecord.accepted << "\n";
        }
        thermalOutfile.close();
    }

    // Output element face energy
    write_element_face_energy_to_csv(model);

    // Output the final structure
    write_final_vertex_data_to_csv(mesh);
    if (mesh.param.isGagScaffoldingEnergyIncluded ||
        mesh.param.isIdealizedProteinLatticeEnergyIncluded)
    {
        mesh.write_gag_scaffolding_state_dat("gag_scaffold_final.dat");
    }
}
