#include "Dynamics.hpp"

#include <cstdint>
#include <stdexcept>
#include <vector>

#include "Counter_rng.hpp"

// The counter-based generator these used to define lives in Counter_rng.hpp
// now: the Metropolis flip sweep needs the same primitives, and two copies of
// a generator is two chances for a run to stop being reproducible.
using slimed::splitmix64;
using slimed::standard_normal;


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
    
    // These used to be re-declared here as locals, which shadowed the members
    // of the same name: the seed from the parameter file never reached the
    // generator that next_step() actually drew from. next_step() no longer
    // uses either member -- see standard_normal() above -- but the seed still
    // has to be recorded, because that is what keys the stream.
    randomSeed = mesh.param.randomSeed;
    gen.seed(randomSeed);
    normal_dist = std::normal_distribution<>{0.0, 1.0};

    // The truncated noise, if asked for: see Param::noiseTruncation. The
    // rescale is worked out once here rather than per draw.
    noiseCutoff = mesh.param.noiseTruncation;
    if (noiseCutoff > 0.0)
    {
        const double variance = slimed::truncated_normal_variance(noiseCutoff);
        noiseRescale = mesh.param.noiseTruncationRenormalize ? 1.0 / std::sqrt(variance) : 1.0;
        std::cout << "[DynamicModel] Brownian noise truncated at |xi| < " << noiseCutoff;
        if (mesh.param.noiseTruncationRenormalize)
        {
            std::cout << ", rescaled by " << noiseRescale << " to unit variance" << std::endl;
        }
        else
        {
            std::cout << ", not rescaled: its variance is " << variance
                      << ", so the membrane is sampled at " << variance
                      << " KBT while its relaxation rates are unchanged" << std::endl;
        }
    }

    nodalForce = mat_calloc(mesh.vertices.size(), 3);
    driveForce = mat_calloc(mesh.vertices.size(), 3);
    select_external_force_vertices();
    setup_bead();
    setup_elastic_box();
    setup_lipid_flow();
    if (mesh.param.meshForceSubsteps > 1)
    {
        // The substeps move the limit surface by M^-T F and add up to the full
        // wall force only when next_step() carries the first 1/n of it, which
        // it does on the FDT-consistent path with the spring in the drive.
        if (!mesh.param.fdtConsistentSurfaceUpdate || !mesh.param.edgeSpringEnabled ||
            mesh.walls_are_hard() || !mesh.param.brownianStepEnabled)
        {
            throw std::runtime_error(
                "[DynamicModel] meshForceSubsteps > 1 needs fdtConsistentSurfaceUpdate = true, "
                "edgeSpringEnabled = true, soft walls and brownianStepEnabled = true.");
        }
        std::cout << "[DynamicModel] meshForceSubsteps = " << mesh.param.meshForceSubsteps
                  << ": the mesh-quality force takes substeps of "
                  << mesh.param.timeStep / mesh.param.meshForceSubsteps << " us" << std::endl;
    }
    std::cout << "[DynamicModel] fdtConsistentSurfaceUpdate = "
              << (mesh.param.fdtConsistentSurfaceUpdate ? "true" : "false") << std::endl;
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

    // ---------------------------------------------------------------------
    // Assemble the nodal force, then move it into the coordinates this step
    // advances.
    //
    // The step below displaces `matSurface`, the Loop limit surface S = M C,
    // and Run_dynamics_flat() maps the result back with C = M^-1 S.  The force
    // the mesh reports is -dE/dC, the gradient with respect to the *control*
    // points, so driving S with it makes the mobility and the noise disagree:
    // in C the update reads
    //     dC = (D dt / kT) M^-1 F_C + sqrt(2 D dt) M^-1 xi,
    // whose mobility is (D/kT) M^-1 but whose noise covariance is 2 D M^-2.
    // Fluctuation-dissipation wants 2 kT * mobility = 2 D M^-1, so the sampled
    // distribution is not exp(-E/kT): mode by mode the stationary variance
    // comes out inflated by 1/m(q), where m(q) is the Fourier symbol of the
    // limit mask -- about 1 at long wavelength and about 4 at the zone corner.
    // That is exactly the shape that lifts the tail of <|h_q|^2> off the q^-4
    // line.
    //
    // The gradient with respect to S is F_S = M^-T F_C, so one application of
    // surface2mesh fixes it -- exactly where M is symmetric, which is every
    // interior row, since the valence-6 limit mask is. It is not symmetric
    // everywhere: assign_mesh2surface() gives a vertex whose adjacent faces are
    // all ghost or boundary the identity row, while its neighbours' rows still
    // carry its 1/12, so max|M - M^T| = 1/12 (reported at startup). Those rows
    // all sit in the ghost ring that postprocess_ghost_periodic() overwrites
    // anyway, and the measured S(q)K(q)/kT lands on 1 to within a few percent,
    // but the standing approximation is M^-1 in place of M^-T.
    //
    // Getting this wrong is not only a bias. A non-reciprocal mobility is not
    // the gradient flow of any potential, so it can do net work on the
    // membrane: as shipped, the r.m.s. height grows without settling instead of
    // reaching a stationary spectrum at all.
    // The bead, if there is one: put where this step wants it, its wall
    // evaluated on the control points as a nodal force. Off, nothing here runs.
    const bool beadOn = mesh.param.beadEnabled;
    if (beadOn)
    {
        update_bead(this->iteration);
    }

    if (mesh.param.fdtConsistentSurfaceUpdate)
    {
        for (int i = 0; i < static_cast<int>(mesh.vertices.size()); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                // The mesh-quality term joins the drive only in fluid mode.
                // The reference-length regularization is deliberately left out
                // of it: it is a solid's restoring force toward a
                // configuration a fluid membrane has no reason to return to,
                // and with the in-plane displacement zeroed it had nothing to
                // act on anyway. The edge spring is a real part of the
                // Hamiltonian the flip sweep samples, so leaving it out of the
                // dynamics would make the two disagree about the energy.
                //
                // The logged E_Total and mean force follow this choice: see
                // Mesh::regularizationInTotals, which DynamicMesh sets from
                // the same flag.
                double f = mesh.vertices[i].force.forceCurvature(j, 0) +
                           mesh.vertices[i].force.forceArea(j, 0);
                if (mesh.param.edgeSpringEnabled)
                {
                    // Substepped, this step carries the first of the n
                    // substeps of the mesh-quality force; mesh_force_substeps()
                    // integrates the rest.
                    if (mesh.param.meshForceSubsteps > 1)
                    {
                        f += mesh.vertices[i].force.forceRegularization(j, 0) /
                             static_cast<double>(mesh.param.meshForceSubsteps);
                    }
                    else
                    {
                        f += mesh.vertices[i].force.forceRegularization(j, 0);
                    }
                }
                if (beadOn)
                {
                    f += beadForce(i, j);
                }
                nodalForce.set(i, j, std::isnan(f) ? 0.0 : f);
            }
        }
        mesh.apply_nodal_force_to_surface(nodalForce, driveForce);
    }

    // The external force, if this iteration is inside its window. It acts
    // on the surface point the step displaces, whichever way the membrane's
    // own force reached it, and only on vertices the step integrates: the
    // branch below skips the rest. Off, it changes nothing.
    const bool pushing = external_force_active(this->iteration);
    const double externalForce[3] = {mesh.param.externalForceX, mesh.param.externalForceY,
                                     mesh.param.externalForceZ};

    // One key per (run, step); standard_normal() mixes the vertex and axis in.
    // Split substep noise: this step draws the first of n parts of the
    // variance, mesh_force_substeps() the rest. Otherwise the shipped scale.
    const double noiseScale =
        (mesh.param.meshForceSubsteps > 1 && mesh.param.meshForceSubstepNoise == "split")
            ? randScaleConst / std::sqrt(static_cast<double>(mesh.param.meshForceSubsteps))
            : randScaleConst;
    const std::uint64_t stepKey =
        splitmix64(static_cast<std::uint64_t>(randomSeed)) +
        static_cast<std::uint64_t>(this->iteration);

    // The vertices the loop below integrates, and their mean height before
    // the step, for the pin: see Param::meanHeightPinned.
    const bool pinMean = mesh.param.meanHeightPinned;
    std::vector<int> integrated;
    double meanHeightBefore = 0.0;
    if (pinMean)
    {
        const bool mixedBC = (mesh.param.boundaryCondition == BoundaryType::Mixed);
        for (int i = 0; i < static_cast<int>(mesh.vertices.size()); i++)
        {
            const bool slaved = (mixedBC || !mesh.param.integratePeriodicDuplicates) &&
                                !mesh.isSlavedPeriodic.empty() && mesh.isSlavedPeriodic[i] != 0;
            if (mesh.vertices[i].isGhost || slaved || mesh.vertices[i].is_fixed())
            {
                continue;
            }
            integrated.push_back(i);
            meanHeightBefore += mesh.matSurface(i, 2);
        }
        if (!integrated.empty())
        {
            meanHeightBefore /= static_cast<double>(integrated.size());
        }
    }

#pragma omp parallel for shared(mesh) schedule(static)
    for (int i = 0; i < mesh.vertices.size(); i++)
    {   
        Matrix forcetermVector(3, 1);
        Matrix randomtermVector(3, 1);

        // A ghost, or the fourth-ring duplicate of a vertex on the far side:
        // either way not an independent coordinate.  The duplicates used to be
        // integrated anyway and then overwritten by
        // postprocess_ghost_periodic().  Overwriting them undoes the step at
        // the duplicate itself, but not before surface2mesh has smeared its
        // fresh, independent random kick across the whole mesh -- so the tile
        // was driven by nVertX*nVertY random numbers per step when only
        // (nFaceX-6)*(nFaceY-6) of its coordinates are free, and every mode
        // came out too warm by roughly that ratio.
        //
        // Under the per-vertex boundary the slaved set is the periodic
        // images, and the legacy flag does not reach them: an image is never
        // a coordinate there. A clamped vertex is not one either.
        const bool mixed = (mesh.param.boundaryCondition == BoundaryType::Mixed);
        const bool isSlaved = (mixed || !mesh.param.integratePeriodicDuplicates) &&
                              !mesh.isSlavedPeriodic.empty() &&
                              mesh.isSlavedPeriodic[i] != 0;
        if (mesh.vertices[i].isGhost || isSlaved || mesh.vertices[i].is_fixed())
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

                    
                forceterm = mesh.param.fdtConsistentSurfaceUpdate
                                ? driveForce(i, j)
                                : mesh.vertices[i].force.forceCurvature(j, 0) +
                                      mesh.vertices[i].force.forceArea(j, 0); // Get force term
                if (!mesh.param.fdtConsistentSurfaceUpdate && mesh.param.edgeSpringEnabled)
                {
                    forceterm += mesh.vertices[i].force.forceRegularization(j, 0);
                }
                if (!mesh.param.fdtConsistentSurfaceUpdate && beadOn)
                {
                    forceterm += beadForce(i, j);
                }
                //std::cout << "Force @ " << i << " , "<< j << " = "  << forceterm << std::endl;
                if (std::isnan(forceterm))
                {
                    forceterm = 0.0;
                }
                if (pushing && externalForceMask[i] != 0)
                {
                    forceterm += externalForce[j];
                }
                if (noiseCutoff > 0.0)
                {
                    randomterm = noiseRescale *
                                 slimed::truncated_standard_normal(stepKey, static_cast<std::uint64_t>(i),
                                                                   static_cast<std::uint64_t>(j), noiseCutoff);
                }
                else
                {
                    randomterm = standard_normal(stepKey, static_cast<std::uint64_t>(i),
                                                 static_cast<std::uint64_t>(j)); // Get random term
                }
                // In-plane motion is switched off by default, which is what
                // this model has always done. It is defensible for a sheet
                // whose triangulation cannot rearrange: the in-plane degrees
                // of freedom have nowhere useful to go, and the analysis in
                // docs/fluctuation_spectrum.md relies on every vertex staying
                // on its ideal lattice site so the height field can be
                // transformed without resampling.
                //
                // It is not defensible for a fluid membrane. In-plane motion
                // is half of what fluidity means, and a flip that rearranges
                // the connectivity while the vertices are pinned laterally is
                // only doing half the job.
                if (j < 2 && !mesh.param.inPlaneDynamicsEnabled) {
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
                    disp = noiseScale * randomtermVector(j,0)
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
    // The pin: take the mean vertical displacement of the step back off
    // every integrated vertex. A uniform shift of the surface points is a
    // uniform shift of the control points (the limit mask's rows sum to
    // one), so the shape is untouched; the images and ghosts follow their
    // sources in postprocess_boundary().
    if (pinMean && !integrated.empty())
    {
        double meanHeightAfter = 0.0;
        for (const int i : integrated)
        {
            meanHeightAfter += mesh.matSurface(i, 2);
        }
        meanHeightAfter /= static_cast<double>(integrated.size());
        const double drift = meanHeightAfter - meanHeightBefore;
        for (const int i : integrated)
        {
            mesh.matSurface.set(i, 2, mesh.matSurface(i, 2) - drift);
        }
    }
    //rsc = 0.04472135955 for ts = 0.001
    //randScaleConst -= 0.0000002;
    //if (randScaleConst <= 0.04462135955){
    //    randScaleConst = 0.0;
    //}
}

/**
 * @brief The other n - 1 substeps of the mesh-quality force (Param::meshForceSubsteps).
 *
 * The full step diverges by an explicit overshoot of the walls, not of bending
 * (docs/resolution_limits.md section 2.2): a control-point force sees a
 * mobility up to 16 D/kT at the zone corner, so a few walls leaning on one
 * vertex already put D k dt / kT past 2. The walls are local sums over the
 * edge table, cheap next to the patch kernels, so they can take n Euler
 * substeps of timeStep / n for the price of n - 1 evaluations of themselves
 * and two surface solves each. Bending and area stay on the full step, which
 * next_step() has already taken together with the first wall substep.
 *
 * Each substep is the full step's pipeline with only the wall force: control
 * net to limit surface, the force mapped onto it (M^-T), the move, the mean
 * height taken back off, and the map back with the boundary re-asserted.
 *
 * The noise (Param::meshForceSubstepNoise): "split" adds sqrt(2 D dt / n) of
 * it to every substep, keyed per substep, so the walls take true
 * Euler-Maruyama steps of dt / n while the bending modes still receive one
 * Gaussian of variance 2 D dt per step. "once" leaves all of it on the first
 * substep; the remaining n - 1 then relax every wall violation without
 * refreshing it, and the walls are sampled too cold (measured 2026-10-20 on
 * the frozen 100 nm tile at dt 0.004: wall energy 82 against 262 pN nm at
 * dt 0.0005; "split" gives 297 at n = 2 and 260 at n = 4). Only "once" can
 * stop early, at the first substep whose wall force is exactly zero.
 */
void DynamicModel::mesh_force_substeps()
{
    const int n = mesh.param.meshForceSubsteps;
    if (n <= 1)
    {
        return;
    }
    const double scale = forceScaleConst / static_cast<double>(n);
    const int nVertices = static_cast<int>(mesh.vertices.size());

    // The coordinates next_step() integrates, by the same rule.
    const bool mixed = (mesh.param.boundaryCondition == BoundaryType::Mixed);
    std::vector<int> integrated;
    for (int i = 0; i < nVertices; i++)
    {
        const bool slaved = (mixed || !mesh.param.integratePeriodicDuplicates) &&
                            !mesh.isSlavedPeriodic.empty() && mesh.isSlavedPeriodic[i] != 0;
        if (mesh.vertices[i].isGhost || slaved || mesh.vertices[i].is_fixed())
        {
            continue;
        }
        integrated.push_back(i);
    }
    if (integrated.empty())
    {
        return;
    }
    const int firstAxis = mesh.param.inPlaneDynamicsEnabled ? 0 : 2;
    const bool splitNoise = (mesh.param.meshForceSubstepNoise == "split");
    const double noiseScale = randScaleConst / std::sqrt(static_cast<double>(n));
    const std::uint64_t stepKey =
        splitmix64(static_cast<std::uint64_t>(randomSeed)) +
        static_cast<std::uint64_t>(this->iteration);

    for (int k = 1; k < n; k++)
    {
        // The walls read the vertex copies of the control net.
        mesh.update_vertices_vector_with_mat();
        // Only the regularization slot is rewritten (the edge spring sets it,
        // the floor and the crease wall add to it); the boundary handling is
        // then idempotent on every other slot, which it already folded and
        // zeroed in the full evaluation.
        mesh.energy_force_fluid_terms();
        if (mixed)
        {
            // manage_force_for_boundary_ghost_vertex() on the one slot: the
            // whole-Force version reassigns every slot of every ghost and
            // costs more than the walls themselves. Same order of additions.
            for (const int v : mesh.periodicImageVertices)
            {
                Matrix &imageForce = mesh.vertices[v].force.forceRegularization;
                Matrix &sourceForce =
                    mesh.vertices[mesh.vertices[v].reflectiveVertexIndex].force.forceRegularization;
                for (int axis = 0; axis < 3; axis++)
                {
                    sourceForce.set(axis, 0, sourceForce(axis, 0) + imageForce(axis, 0));
                    imageForce.set(axis, 0, 0.0);
                }
            }
            for (Vertex &vertex : mesh.vertices)
            {
                if (vertex.is_fixed() || vertex.type == VertexType::Ghost || vertex.isGhost)
                {
                    for (int axis = 0; axis < 3; axis++)
                    {
                        vertex.force.forceRegularization.set(axis, 0, 0.0);
                    }
                }
            }
        }
        else
        {
            mesh.manage_force_for_boundary_ghost_vertex();
        }

        bool anyForce = false;
        for (int i = 0; i < nVertices; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                const double f = mesh.vertices[i].force.forceRegularization(j, 0);
                const double clean = std::isnan(f) ? 0.0 : f;
                anyForce = anyForce || (clean != 0.0);
                nodalForce.set(i, j, clean);
            }
        }
        if (!anyForce && !splitNoise)
        {
            break;
        }
        if (anyForce)
        {
            mesh.apply_nodal_force_to_surface(nodalForce, driveForce);
        }
        else
        {
            driveForce.set_all(0.0);
        }
        mesh.apply_mesh_to_surface();
        // A key per substep, hashed so that it collides neither with another
        // step's key (consecutive integers) nor with another substep's.
        const std::uint64_t substepKey =
            splitmix64(stepKey ^ (0x9E3779B97F4A7C15ULL * static_cast<std::uint64_t>(k)));

        double meanHeightBefore = 0.0;
        for (const int i : integrated)
        {
            meanHeightBefore += mesh.matSurface(i, 2);
        }
        meanHeightBefore /= static_cast<double>(integrated.size());
        for (const int i : integrated)
        {
            for (int j = firstAxis; j < 3; j++)
            {
                double disp = scale * driveForce(i, j);
                if (splitNoise)
                {
                    disp += noiseScale * standard_normal(substepKey, static_cast<std::uint64_t>(i),
                                                         static_cast<std::uint64_t>(j));
                }
                mesh.matSurface.set(i, j, mesh.matSurface(i, j) + disp);
            }
        }
        if (mesh.param.meanHeightPinned)
        {
            double meanHeightAfter = 0.0;
            for (const int i : integrated)
            {
                meanHeightAfter += mesh.matSurface(i, 2);
            }
            meanHeightAfter /= static_cast<double>(integrated.size());
            const double drift = meanHeightAfter - meanHeightBefore;
            for (const int i : integrated)
            {
                mesh.matSurface.set(i, 2, mesh.matSurface(i, 2) - drift);
            }
        }
        mesh.apply_surface_to_mesh();
        mesh.postprocess_boundary();
    }
}
