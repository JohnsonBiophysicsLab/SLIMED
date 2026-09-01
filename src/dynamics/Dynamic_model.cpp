#include "Dynamics.hpp"

#include <cstdint>

namespace
{
/// SplitMix64 -- one avalanche round on a 64-bit counter.
inline std::uint64_t splitmix64(std::uint64_t x)
{
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

/// A uniform on (0, 1) -- never exactly 0, so the log() below stays finite.
inline double uniform_open01(std::uint64_t bits)
{
    return (static_cast<double>(bits >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

/// One standard normal, keyed by (run, iteration, vertex, axis).
///
/// Counter-based rather than sequential: what a vertex draws depends only on
/// where it sits in the run, never on the order the OpenMP team happens to
/// reach it. The previous code called a shared std::normal_distribution on a
/// shared std::mt19937 from inside `#pragma omp parallel for`, which is a
/// data race on the generator state: the noise was neither reproducible nor
/// guaranteed to still be Gaussian, and an equilibrium fluctuation spectrum
/// is only ever as good as the noise that drives it.
inline double standard_normal(std::uint64_t stepKey, std::uint64_t vertex, std::uint64_t axis)
{
    const std::uint64_t key = splitmix64(stepKey + splitmix64(vertex * 4ULL + axis));
    const double u1 = uniform_open01(splitmix64(key));
    const double u2 = uniform_open01(splitmix64(key ^ 0xD1B54A32D192ED03ULL));
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
}
} // namespace

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

    nodalForce = mat_calloc(mesh.vertices.size(), 3);
    driveForce = mat_calloc(mesh.vertices.size(), 3);
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
    if (mesh.param.fdtConsistentSurfaceUpdate)
    {
        for (int i = 0; i < static_cast<int>(mesh.vertices.size()); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                double f = mesh.vertices[i].force.forceCurvature(j, 0) +
                           mesh.vertices[i].force.forceArea(j, 0);
                nodalForce.set(i, j, std::isnan(f) ? 0.0 : f);
            }
        }
        driveForce = mesh.surface2mesh * nodalForce;
    }

    // One key per (run, step); standard_normal() mixes the vertex and axis in.
    const std::uint64_t stepKey =
        splitmix64(static_cast<std::uint64_t>(randomSeed)) +
        static_cast<std::uint64_t>(this->iteration);

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
        const bool isSlaved = !mesh.param.integratePeriodicDuplicates &&
                              !mesh.isSlavedPeriodic.empty() &&
                              mesh.isSlavedPeriodic[i] != 0;
        if (mesh.vertices[i].isGhost || isSlaved)
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
                //std::cout << "Force @ " << i << " , "<< j << " = "  << forceterm << std::endl;
                if (std::isnan(forceterm))
                {
                    forceterm = 0.0;
                }
                randomterm = standard_normal(stepKey, static_cast<std::uint64_t>(i),
                                             static_cast<std::uint64_t>(j)); // Get random term
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