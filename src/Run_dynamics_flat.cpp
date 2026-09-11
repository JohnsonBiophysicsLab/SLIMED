#include "Run_simulation.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

void run_dynamics_flat(std::string param_filename) {

    Param inputParam;
    import_param_file(inputParam, param_filename + ".params");
    DynamicMesh mesh(inputParam);
    mesh.setup_flat();
    for (Vertex& vertex: mesh.vertices){
        vertex.coord.set(2, 0, 10.0);
    }

    //
    /*
    int n = round(inputParam.sideX/inputParam.lFace); double dx = inputParam.sideX/n; // x axis division
    double a = dx;
    double dy = sqrt(3.0)/2.0 * a; 
    int m = round(inputParam.sideY/dy);      // y axis division

    int nFaceX = n;
    int nFaceY = m;

    for (int i = 0; i <= nFaceX; i++){
        for (int j = 0; j <= nFaceY; j++){
            if (!(mesh.vertices[(nFaceX+1)*j + i].isGhost)){
                
                double zOrig = mesh.vertices[(nFaceX+1)*j + i].coord(2, 0);
                double shift = 80.0 * (i - 8) * (j - 8) * (nFaceX - i - 8) * (nFaceY - j - 8) / 20000.0;
                if (i < 8 || i > nFaceX - 8 || j < 8 || j > nFaceY - 8){
                    shift = 0.0;
                }
                mesh.vertices[(nFaceX+1)*j + i].coord.set(2, 0 , zOrig + shift);
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
        // move spline points upwards until the lower boundary is approximately z=0
        mesh.move_vertices_based_on_scaffolding();
        // find closest vertex point and save in vector
        mesh.set_scaffolding_vertices_correspondence();
    }

    // Output the vertices and faces matrix
    mesh.write_faces_csv(param_filename + "face.csv");
    mesh.write_vertices_csv(param_filename + "vertex_begin.csv");
    mesh.write_vertices_csv_with_type(param_filename + "vertex_type_begin.csv");
    
    // Initialize all value before minimum energy search 
    mesh.calculate_element_area_volume(); // Calculate the elemental area and volume per triangles (faces)
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0); // Calculate initial area and volume
    mesh.report_volume_rebaseline(); // TEMPORARY -- see docs/volume_functional_split.md step 4
    mesh.update_previous_coord_for_vertex(); // Update previous coordinate...
    mesh.update_reference_coord_from_previous_coord(); // ...and reference coord
    mesh.Compute_Energy_And_Force();// Calculate energy on faces and force on vertices
    mesh.update_previous_coord_for_vertex(); // Update previous coordinate
    mesh.update_previous_force_for_vertex(); // Update forces

    // Create a "Record" object to bookkeep the energies and forces
    Record record(mesh.param.maxIterations); ///< Bookkeeping for iterations
    record.add(mesh.param.area, mesh.param.energy, mesh.calculate_mean_force()); ///< Record the data before search (N=0)

    // Create a "Model" container object that packs mesh and record
    DynamicModel model(mesh, record);

    ///////////////////////////////////////
    mesh.update_vertices_mat_with_vector(); //Update matMesh
    mesh.apply_mesh_to_surface(); //Update matSurface -- dense or sparse, per surfaceSolver

    if (mesh.param.meshpointOutput) {
        dynamics_create_trajectory_files(mesh, param_filename);
        if (mesh.param.edgeFlipEnabled) {
            // The baseline a later frame is compared against. Only for a run
            // whose connectivity can move: with flips off this file would be a
            // byte-for-byte copy of the face.csv written above, and a run that
            // was reproducing an earlier one would find an output it did not
            // have before.
            dynamics_output_face_frame(mesh, param_filename, 0);
        }
    }

    // One line per flip attempt, accepted or not. Held for the whole run and
    // appended in blocks rather than per sweep: a sweep writes a handful of
    // lines and the open/close would cost more than the sweep.
    std::vector<EdgeFlipRecord> edgeFlipLog;
    EdgeFlipSweepStats edgeFlipTotals;
    const std::string edgeFlipLogPath = param_filename + "EdgeFlips.csv";
    if (mesh.param.edgeFlipEnabled) {
        // Appending, so a stale log from a previous run of the same name would
        // be read back as part of this one.
        std::remove(edgeFlipLogPath.c_str());
    }


    for (model.iteration = 0; model.iteration < mesh.param.maxIterations; model.iteration ++) {

        //1.control mesh to limit surface
        mesh.apply_mesh_to_surface();

        //2.next time step - calculate displacement on limit surface
        model.next_step(); 

        //3.limit surface to control mesh: verticesOnMesh = surface2mesh * verticesProjSurface
        mesh.apply_surface_to_mesh();

        //4. postprocessing based on boundary condition
        switch (mesh.param.boundaryCondition) {
            case BoundaryType::Periodic:
                mesh.postprocess_ghost_periodic();
                break;
            default:
                // Handle default boundary condition type.
                break;
        }
        
        mesh.update_vertices_vector_with_mat();
        //@todo move this to io.hpp
        //4.update values of vertex (vector of double) with verticesOnMesh (gsl matrix)
        if (mesh.param.meshpointOutput &&
            (model.iteration + 1) % mesh.param.meshpointOutputInterval == 0) {
            // A frame every step is only affordable for short runs: at the
            // 10^6-10^7 iterations an equilibrium spectrum needs, one line per
            // step per vertex runs to tens of gigabytes.
            dynamics_output_trajectory_files(mesh, param_filename);
            if (mesh.param.edgeFlipEnabled) {
                // Beside the frame, and only when the connectivity moved since
                // the last one: a coordinate frame of a fluid run is
                // meaningless without the triangles it belongs to.
                dynamics_output_face_frame(mesh, param_filename, model.iteration + 1);
            }
        }

        // Record the Energy and nodal Force
        record.add(mesh.param.area, mesh.param.energy, mesh.calculate_mean_force());

        mesh.Compute_Energy_And_Force();

        // The flip sweep runs between force evaluations, on a mesh whose
        // energies and forces are current, and at fixed coordinates. It needs
        // the current per-face energies because the trial energy is a local
        // difference against them, and it invalidates them by construction, so
        // an accepted flip is followed by a second evaluation before the next
        // Brownian step reads a force.
        if (mesh.param.edgeFlipEnabled &&
            (model.iteration + 1) % std::max(1, mesh.param.edgeFlipInterval) == 0) {
            const EdgeFlipSweepStats stats = mesh.edge_flip_sweep(model.iteration, &edgeFlipLog);
            edgeFlipTotals.drawn += stats.drawn;
            edgeFlipTotals.attempted += stats.attempted;
            edgeFlipTotals.accepted += stats.accepted;
            edgeFlipTotals.deltaEnergy += stats.deltaEnergy;
            if (stats.accepted > 0) {
                // Every cache keyed on the connectivity -- the one-rings, the
                // device layout, the sparse limit mask -- notices the version
                // bump on its own. What does not is the stored energy and
                // force, which still describe the old triangulation.
                mesh.Compute_Energy_And_Force();
            }
            if (edgeFlipLog.size() >= 4096) {
                write_edge_flip_log_csv(edgeFlipLog, edgeFlipLogPath);
                edgeFlipLog.clear();
            }
        }

        // A run that has gone non-finite has nothing left to say, and it says
        // it for as long as it is given: the 3000-step fluid run that found
        // this wrote NaN into every row after step 1535 and took four minutes
        // doing it. Stop at the first one, keeping what was collected.
        if (!std::isfinite(mesh.param.energy.energyTotal)) {
            std::cerr << "[run_dynamics_flat] The total energy is "
                      << mesh.param.energy.energyTotal << " at iteration " << model.iteration
                      << "; the run has diverged and is stopping. Everything up to here has "
                         "been written. A membrane that blows up this way is usually held "
                         "together too weakly for the moves it is being asked to make -- check "
                         "edgeSpringConstant against the startup barrier report."
                      << std::endl;
            model.iteration++;
            break;
        }

        cout<<"=========ITERATION:" << model.iteration << "==========================" << endl;

    }

    if (mesh.param.edgeFlipEnabled) {
        write_edge_flip_log_csv(edgeFlipLog, edgeFlipLogPath);
        edgeFlipLog.clear();
        std::cout << "[run_dynamics_flat] Edge flips: " << edgeFlipTotals.accepted
                  << " accepted of " << edgeFlipTotals.attempted << " admissible attempts ("
                  << edgeFlipTotals.drawn << " drawn), acceptance "
                  << edgeFlipTotals.acceptance() << ", total dE " << edgeFlipTotals.deltaEnergy
                  << " pN.nm." << std::endl;
    }

    // Output Energy and meanforce
    write_energy_force_data_to_csv(model);

    // write final vertices csv
    mesh.write_vertices_csv(param_filename + "vertex_final.csv");
    mesh.write_vertices_csv_with_type(param_filename + "vertex_type_final.csv");

    //currently running python3 code to convert to xyz file
    /*if (mesh.param.xyzOutput) {
        //@TODO add functionality
    }*/
}
