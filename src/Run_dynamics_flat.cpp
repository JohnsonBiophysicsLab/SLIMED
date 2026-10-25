#include "Run_simulation.hpp"
#include "dynamics/Divergence_guard.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

void run_dynamics_flat(std::string param_filename) {

    Param inputParam;
    import_param_file(inputParam, param_filename + ".params");
    DynamicMesh mesh(inputParam);
    // The generated flat sheet, or the mesh files when the parameters name them.
    setup_mesh_from_parameters(mesh);
    if (!mesh_comes_from_files(mesh.param))
    {
        // The generated sheet is lifted off z = 0. An imported mesh keeps the
        // geometry its file gave it.
        for (Vertex& vertex: mesh.vertices){
            vertex.coord.set(2, 0, 10.0);
        }
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

    // Whatever the initialization above moved, the periodic images follow
    // their sources. A no-op on a mesh without images.
    mesh.sync_periodic_images();

    // Output the vertices and faces matrix
    mesh.write_faces_csv(param_filename + "face.csv");
    mesh.write_vertices_csv(param_filename + "vertex_begin.csv");
    mesh.write_vertices_csv_with_type(param_filename + "vertex_type_begin.csv");
    if (mesh.param.boundaryCondition == BoundaryType::Mixed)
    {
        // The mesh with its boundary types, in the form meshVerticesFile and
        // meshFacesFile read back.
        export_mesh_to_vertices_faces(mesh, param_filename + "mesh_vertices.csv",
                                      param_filename + "mesh_faces.csv");
    }
    
    // Initialize all value before minimum energy search 
    mesh.calculate_element_area_volume(); // Calculate the elemental area and volume per triangles (faces)
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0); // Calculate initial area and volume
    if (mesh.lipid_flow_on())
    {
        // Every face's relaxed area, from the areas just computed: the
        // uniform strain the global constraint would start from.
        mesh.initialize_lipid_flow();
    }
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
    if (mesh.param.externalForceEnabled)
    {
        // The pushed region, so an analysis can tell it from the rest.
        model.write_external_force_vertices_csv(param_filename + "external_force_vertices.csv");
    }
    // The bead log: one row per trajectory frame, the first for the
    // configuration before any step.
    const std::string beadLogPath = param_filename + "bead.csv";
    if (mesh.param.beadEnabled)
    {
        std::remove(beadLogPath.c_str());
        model.update_bead(0);
        model.append_bead_log(beadLogPath, 0);
    }
    // The box log: one row per trajectory frame, the first for the
    // configuration before any step, each carrying the tension the step
    // just taken was driven by.
    const std::string boxLogPath = param_filename + "box.csv";
    if (mesh.param.elasticBoxEnabled)
    {
        std::remove(boxLogPath.c_str());
        model.refresh_box_tension();
        model.append_box_log(boxLogPath, 0);
    }
    // The lipid log: totals and tension statistics per trajectory frame, and
    // with localOutput the relaxed and actual area of every face.
    const std::string lipidLogPath = param_filename + "lipid.csv";
    const std::string lipidFacePath = "lipidface" + param_filename + ".csv";
    if (mesh.lipid_flow_on())
    {
        std::remove(lipidLogPath.c_str());
        std::remove(lipidFacePath.c_str());
        model.append_lipid_log(lipidLogPath, 0);
        if (mesh.param.localOutput)
        {
            model.append_lipid_face_frame(lipidFacePath);
        }
    }

    ///////////////////////////////////////
    mesh.update_vertices_mat_with_vector(); //Update matMesh
    mesh.apply_mesh_to_surface(); //Update matSurface -- dense or sparse, per surfaceSolver

    if (mesh.param.meshpointOutput) {
        dynamics_create_trajectory_files(mesh, param_filename);
        if (mesh.param.localOutput) {
            // The local quantities of the same frames: per-face energies and
            // per-vertex forces, as the model stores them.
            dynamics_create_local_files(mesh, model, param_filename);
        }
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


    // The Monte Carlo vertex moves: totals for the run, and since the last
    // trajectory frame for the log beside it.
    VertexMoveSweepStats vertexMoveTotals;
    VertexMoveSweepStats vertexMoveSinceFrame;
    const std::string vertexMoveLogPath = param_filename + "VertexMoves.csv";
    if (mesh.param.vertexMoveEnabled) {
        std::remove(vertexMoveLogPath.c_str());
    }
    // Hard walls are an invariant the moves preserve, not one they establish:
    // a move is rejected when it leaves a wall broken, so a mesh that starts
    // outside the admissible set would have its offending neighbourhood
    // frozen for the whole run and report nothing but a low acceptance.
    if (mesh.walls_are_hard()) {
        const WallViolations broken = mesh.count_wall_violations();
        if (broken.total() > 0) {
            throw std::runtime_error(
                "[run_dynamics_flat] meshWallType = hard, but the starting mesh breaks " +
                std::to_string(broken.tether) + " tether, " + std::to_string(broken.altitude) +
                " altitude and " + std::to_string(broken.crease) +
                " crease walls. Start from a mesh inside the walls, or widen them "
                "(edgeTetherMinRatio/MaxRatio, triangleShapeMinAltitudeRatio, creaseWallAngle).");
        }
    }

    // The divergence bounds, from the mesh as the first step finds it.
    const DivergenceGuard divergenceGuard(mesh);

    for (model.iteration = 0; model.iteration < mesh.param.maxIterations; model.iteration ++) {

        // A pure Monte Carlo run (brownianStepEnabled = false) skips the
        // step: the control net is moved by the sweeps at the bottom of the
        // loop and nothing else.
        if (mesh.param.brownianStepEnabled) {
        //1.control mesh to limit surface
        mesh.apply_mesh_to_surface();

        //2.next time step - calculate displacement on limit surface
        model.next_step();

        //3.limit surface to control mesh: verticesOnMesh = surface2mesh * verticesProjSurface
        mesh.apply_surface_to_mesh();

        //4. postprocessing based on boundary condition
        mesh.postprocess_boundary();

        //4a. the rest of the mesh-quality force's substeps, at
        // meshForceSubsteps > 1. Off, nothing runs.
        if (mesh.param.meshForceSubsteps > 1) {
            model.mesh_force_substeps();
        }

        //4b. the elastic box: one Brownian step of the periodic strain from
        // the same force evaluation the vertex step used, then the rescale
        // of the whole net, images and ghosts included. Off, nothing runs.
        if (mesh.param.elasticBoxEnabled) {
            model.update_elastic_box();
        }
        }

        // The relaxed areas take their step from the same force evaluation:
        // lateral flux down the tension gradient, the reservoir, the noise.
        if (mesh.lipid_flow_on()) {
            model.advance_lipid_flow();
        }

        mesh.update_vertices_vector_with_mat();
        //@todo move this to io.hpp
        //4.update values of vertex (vector of double) with verticesOnMesh (gsl matrix)
        if (mesh.param.meshpointOutput &&
            (model.iteration + 1) % mesh.param.meshpointOutputInterval == 0) {
            if (!mesh.param.brownianStepEnabled) {
                // The frame writes the limit surface beside the control net,
                // and with no Brownian step nothing else has refreshed it.
                mesh.apply_mesh_to_surface();
            }
            if (mesh.param.vertexMoveEnabled) {
                const WallViolations broken =
                    mesh.param.edgeSpringEnabled ? mesh.count_wall_violations() : WallViolations();
                write_vertex_move_log_csv(vertexMoveLogPath, model.iteration + 1,
                                          vertexMoveSinceFrame, broken);
                vertexMoveSinceFrame = VertexMoveSweepStats();
            }
            // A frame every step is only affordable for short runs: at the
            // 10^6-10^7 iterations an equilibrium spectrum needs, one line per
            // step per vertex runs to tens of gigabytes.
            dynamics_output_trajectory_files(mesh, param_filename);
            if (mesh.param.localOutput) {
                dynamics_output_local_frame(mesh, model, param_filename);
            }
            if (mesh.param.beadEnabled) {
                // The wall re-evaluated on the frame just written, at the
                // centre the step used.
                model.refresh_bead_contact();
                model.append_bead_log(beadLogPath, model.iteration + 1);
            }
            if (mesh.param.elasticBoxEnabled) {
                model.append_box_log(boxLogPath, model.iteration + 1);
            }
            if (mesh.lipid_flow_on()) {
                model.append_lipid_log(lipidLogPath, model.iteration + 1);
                if (mesh.param.localOutput) {
                    model.append_lipid_face_frame(lipidFacePath);
                }
            }
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
        // Brownian step reads a force. That second evaluation touches only
        // the faces the flips changed when Param::edgeFlipIncrementalForces
        // allows it; a vertex-move sweep moves coordinates and so needs the
        // whole mesh again.
        bool flipsChangedMesh = false;
        bool movesChangedMesh = false;
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
                flipsChangedMesh = true;
            }
            if (edgeFlipLog.size() >= 4096) {
                write_edge_flip_log_csv(edgeFlipLog, edgeFlipLogPath);
                edgeFlipLog.clear();
            }
        }

        // The vertex-move sweep, after the flips and on the connectivity they
        // left. Each is a kernel that leaves exp(-E/kT) invariant, so the
        // order between them is a matter of bookkeeping: both keep
        // param.area and param.vol current, which is all the next trial's
        // constraint difference reads.
        if (mesh.param.vertexMoveEnabled &&
            (model.iteration + 1) % std::max(1, mesh.param.vertexMoveInterval) == 0) {
            const VertexMoveSweepStats stats = mesh.vertex_move_sweep(model.iteration);
            vertexMoveTotals += stats;
            vertexMoveSinceFrame += stats;
            if (stats.accepted > 0) {
                // The moves act on the vertices; the step and the trajectory
                // read the control net from matMesh.
                mesh.update_vertices_mat_with_vector();
                movesChangedMesh = true;
            }
        }

        if (movesChangedMesh) {
            mesh.Compute_Energy_And_Force();
        } else if (flipsChangedMesh) {
            // Falls back to the whole-mesh pass when the incremental one is
            // not available (parameter off, device backend, local lipid flow).
            mesh.Compute_Energy_And_Force_After_Flips();
        }

        // A run that has gone non-finite has nothing left to say, and it says
        // it for as long as it is given: the 3000-step fluid run that found
        // this wrote NaN into every row after step 1535 and took four minutes
        // doing it. Stop at the first one, keeping what was collected. Nor
        // does a blow-up always reach NaN: diag_fixedbox_dt0.005 stayed finite
        // at E ~ 1e206 for 80k steps, so the guard also stops a run whose
        // coordinates or energy have left any range a healthy run visits.
        const std::string divergence = divergenceGuard.check(mesh);
        if (!divergence.empty()) {
            std::cerr << "[run_dynamics_flat] " << divergence << " at iteration " << model.iteration
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

    if (mesh.param.vertexMoveEnabled) {
        std::cout << "[run_dynamics_flat] Vertex moves: " << vertexMoveTotals.accepted
                  << " accepted of " << vertexMoveTotals.attempted << " attempts, acceptance "
                  << vertexMoveTotals.acceptance() << "; " << vertexMoveTotals.wallRejected
                  << " refused by a hard wall; total dE " << vertexMoveTotals.deltaEnergy
                  << " pN.nm." << std::endl;
        if (mesh.param.edgeSpringEnabled) {
            const WallViolations broken = mesh.count_wall_violations();
            std::cout << "[run_dynamics_flat] Walls at the end of the run: " << broken.tether
                      << " tether, " << broken.altitude << " altitude, " << broken.crease
                      << " crease outside ("
                      << (mesh.walls_are_hard() ? "hard walls: all three must be zero"
                                                : "soft walls: these are leaned on, not broken")
                      << ")." << std::endl;
        }
    }

    if (mesh.lipid_flow_on()) {
        std::cout << "[run_dynamics_flat] Lipid flow: relaxed area "
                  << (mesh.lipid_flow_local() ? mesh.total_face_area0() : mesh.param.area0)
                  << " nm^2 against a membrane area of " << mesh.param.area << " nm^2; "
                  << model.lipidReservoirExchanged << " nm^2 taken from the reservoir; "
                  << model.lipidClampedUpdates << " updates hit the positivity floor." << std::endl;
    }
    if (mesh.param.elasticBoxEnabled) {
        const double meanTension =
            (model.boxSteps > 0) ? model.boxTensionRunningSum / model.boxSteps : 0.0;
        std::cout << "[run_dynamics_flat] Elastic box: strain (" << mesh.boxScale[0] << ", "
                  << mesh.boxScale[1] << "), box " << mesh.boxLength[0] << " x "
                  << mesh.boxLength[1] << " nm, projected area " << mesh.box_projected_area()
                  << " nm^2 against a membrane area of " << mesh.param.area
                  << " nm^2; mean tension over " << model.boxSteps << " box steps "
                  << meanTension << " pN/nm (target " << mesh.param.elasticBoxTension
                  << "); " << model.boxClippedSteps
                  << " steps clipped at elasticBoxMaxStrainPerStep." << std::endl;
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
