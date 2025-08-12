#include "Run_simulation.hpp"

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
    mesh.matSurface = mesh.mesh2surface * mesh.matMesh; //Update matSurface

    if (mesh.param.meshpointOutput) {
        dynamics_create_trajectory_files(mesh, param_filename);
    }
    

    for (model.iteration = 0; model.iteration < mesh.param.maxIterations; model.iteration ++) {

        //1.control mesh to limit surface
        mesh.matSurface = mesh.mesh2surface * mesh.matMesh;

        //2.next time step - calculate displacement on limit surface
        model.next_step(); 

        //3.limit surface to control mesh: verticesOnMesh = surface2mesh * verticesProjSurface
        mesh.matMesh = mesh.surface2mesh * mesh.matSurface;

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
        if (mesh.param.meshpointOutput) {
            dynamics_output_trajectory_files(mesh, param_filename);
        }

        // Record the Energy and nodal Force
        record.add(mesh.param.area, mesh.param.energy, mesh.calculate_mean_force());

        mesh.Compute_Energy_And_Force();
        cout<<"=========ITERATION:" << model.iteration << "==========================" << endl;

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