/**
 * @file Parameters.hpp
 * @author Y Ying (yying7@jh.edu)
 * @author Y Fu (yfu31@jh.edu)
 * @brief This file defines essential parameters used in continuum membrane
 *        as well as membrane dynamics code. The Param structure defines
 *        all the essential physical constants and simulation parameters for
 *        the model; multiple instances of Param help with parallel
 *        computing of the model.  The Mesh class includes vertices and faces information
 *        for the control mesh of the model.
 * @date 2023-01-05
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <math.h>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "energy_force/Energy.hpp"
#include "energy_force/Force.hpp"

#include "linalg/Linear_algebra.hpp"


/**
 * @brief The type of boundary condition for the simulation.
 */
enum class BoundaryType
{
    /** Fixed boundary condition. */
    Fixed,

    /** Periodic boundary condition. */
    Periodic,

    /** Free boundary condition. */
    Free
};

/**
 * @brief Bookkeeps the shape and area deformation count of Mesh.
 *
 * Shape deformation is calculated by measuring the deviation of normalized side length
 * of mesh triangles from equilateral triangles. Area deformation is calculated by
 * measuring the deviation of triangle areas from their equilibrium values.
 *
 * Note that this struct is only used for development and testing purposes and is not
 * used in the simulation itself.
 * - Yiben
 */
struct DeformationCount
{
    int shapeDeformCount = 0; ///< Number of triangles with shape deformation
    int areaDeformCount = 0;  ///< Number of triangles with area deformation
    int noDeformCount = 0;    ///< Number of triangles with no deformation
};

/**
 * @brief Represents a single diffusing particle in a simulation.
 *
 * A Particle object contains information about the particle's position, velocity,
 * and diffusion constant. The position is represented by a 3x1 matrix of coordinates,
 * while the velocity is represented by another 3x1 matrix of velocities in the x, y,
 * and z directions.
 */
class Particle
{
public:
    int index;                   ///< Index of the particle in the simulation
    int faceIndex;               ///< Index of the face that the particle is located on
    Matrix vwu = Matrix(3, 1);   ///< Velocity vector of the particle in nm/us
    Matrix coord = Matrix(3, 1); ///< Coordinate vector of the particle in nm
    double D = 0.0;              ///< Diffusion constant of the particle in nm^2/us
};

/**
 * @brief Contains simulation parameters and physical constants.
 *
 * The Param struct stores a variety of parameters and physical constants that are used
 * to initialize and run a membrane simulation. This includes properties like the membrane
 * size, shape, and subdivision, as well as boundary conditions, time step, and diffusion
 * constant. Many of these variables have default values but can be customized by changing
 * the corresponding member variables in a Param object.
 *
 * Note that some variables in this struct are deprecated or unused and should not be relied
 * upon for correct behavior of the simulation. These variables may be removed or modified
 * in future versions of the code.
 */
struct Param
{
    // developer options
    bool VERBOSE_MODE = true; ///< Whether to print verbose output during simulation
    int maxIterations = 1E5; ///< Max number of iteration
    std::string restartInputFile = ""; ///< Optional checkpoint file to restart from
    std::string checkpointOutputFile = "slimed_restart.chk"; ///< Checkpoint file written during minimization
    int checkpointOutputInterval = 1000; ///< Iteration interval for checkpoint writes; <=0 disables checkpointing

    // physical constants
    double kCurv = 83.4;  ///< Bending modulus i.e. curvature constant (kc)
    double uSurf = 250.0; ///< Surface area constraint i.e. surface constant (us)
    double uVol = 0.0;    ///< Volume constraint i.e. volume constant (uv)
    double kReg = 83.4;    ///< Coefficient of the regularization constraint (k)
    double kSpring;       ///< Spring constant for insertion zones (K)
    bool setRelaxAreaToDefault = false; ///< true to set area0 equal to area of starting config
    double area0;         ///< Target area for membrane (S0)
    double area;          ///< Total area of the membrane (S)
    double vol0;        ///< Target volume for membrane (V0)
    double vol;         ///< Total volume of the membrane (V)
    double insertCurv;  ///< Spontaneous curvature of insertions (C0)
    double spontCurv;   ///< Spontaneous curvature of membrane (c0)

    // membrane size and axes division
    double sideX = 100.0;                           ///< X-axis length for flat membrane
    double sideY = 100.0;                           ///< Y-axis length for flat membrane
    double radius = 25.0;                          ///< Radius for spherical membrane
    double lFace = 5.0;                           ///< lFace
    int nFaceX = -1;                        ///< Number of faces (edges) along X axis for flat membrane
    int nFaceY = -1;                        ///< Number of faces (edges) along Y axis for flat membrane
    double dFaceX;                          ///< Initial actual face side length along X axis for flat membrane
    double dFaceY;                          ///< Initial actual face side length along Y axis for flat membrane
    double meanL;                           ///< Mean length of edges after subdivision
    double sigma = 0.0;                     ///< Noise level for vertex positions
    bool isInsertionAreaConstraint = false; ///< Whether to apply area constraint to insertions
    bool isAdditiveScheme = false;          ///< Whether to use additive scheme for constraints
    bool isGlobalConstraint = true;                ///< Whether to apply global constraint across entire membrane

    /**
     * @brief Which implementation evaluates the per-face energies and forces.
     *
     * "cpu"  -- the loops in Compute_Energy_And_Force(), the default.
     * "gpu"  -- slimed::CudaForceBackend. Fails loudly if the binary was built
     *           without CUDA or no device is visible, because a run that
     *           silently fell back would be reported as a GPU timing.
     * "auto" -- the GPU when one is usable, the CPU otherwise, saying which it
     *           chose.
     *
     * All three compute the same thing: the device path runs the same kernel
     * bodies over the same flattened mesh, and the tests pin it against the
     * CPU. See docs/gpu_acceleration.md for where the GPU actually wins --
     * below roughly 10^5 faces it does not.
     */
    std::string forceBackend = "cpu";
    double elementTriangleArea0;            ///< Target area for individual triangles

    // gauss quadrature
    int gaussQuadratureN = 2;    ///< Number of Gaussian quadrature points to use
    Matrix VWU;                  ///< (N,3) matrix of vertex coordinates and weights
    Matrix gaussQuadratureCoeff; ///< (N,1) matrix of Gaussian quadrature coefficients

    // shape function
    std::vector<Matrix> shapeFunctions; ///< List of shape functions for each triangle

    // boundary conditions
    BoundaryType boundaryCondition = BoundaryType::Periodic; ///< Type of boundary condition ("Fixed", "Periodic", "Free")

    // optimization methods
    bool usingNCG = true;      ///< Whether to use nonlinear conjugate gradient method
    bool isNCGstuck = false; ///< Whether NCG method has gotten stuck

    // convergence criteria
    double deltaEnergyConverge = 1e-5; ///< Convergence criteria for total energy
    double deltaForceScaleConverge = 1e-5; ///< Convergence criteria for max force scale

    // regularization
    double gamaShape = 0.2;                       ///< Shape deformation coefficient
    double gamaArea = 0.2;                        ///< Area deformation coefficient
    /**
     * @brief Apply one global Loop refinement to the control mesh at setup.
     *
     * Off by default. Only needed for meshes whose extraordinary vertices are
     * adjacent, which the valence histogram reports at setup. It roughly
     * quadruples the vertex count and changes the dynamical degrees of
     * freedom, so enabling it rebaselines the run.
     */
    bool isPreRefinementEnabled = false;

    bool usingRpi = true;             ///< Whether to use R-pi adaptive method for regularization energy;

    // spline points
    bool xyzOutput = true;                                       ///< Whether to output XYZ coordinates of vertices
    bool meshpointOutput = true;                                 ///< Whether to output meshpoints
    int meshpointOutputInterval = 1;                             ///< Write a trajectory frame every N dynamics iterations
    bool isEnergyHarmonicBondIncluded = false;                    ///< Whether to include energy of spline points
    std::vector<Matrix> scaffoldingPoints;                       ///< List of spline points
    std::vector<int> scaffoldingPoints_correspondingVertexIndex; ///< Corresponding vertex indices for spline points
    double scaffoldingSphereRaidus = 50.0;                       ///< spherical radius of the scaffolding lattice cap
    double springConst = 4.5;                                   ///< Spring constant for harmonic bond potential
    double springConstScaler = 1.31607401;                              ///< Spring constant scaler per iteration interval
    int springConstScalingInterv = 200;                          ///< Iteration interval for spring constant scaling
    int springConstUpperBound = 1000.0;                          ///< Max spring constant to stop scaling of spring constant
    int propagateScaffoldingInterv = 7;                         ///< Iteration interval for repositioning scaffolding based on energy and force
    int propagateScaffoldingNstep = 100;                         ///< Number of steps every time propagating scaffolding
    double splinePointsZcoordScaling = 0.0;                     // for z scaling
    double lbond = 9.0;                                          ///< Harminc bond length in nm
    double relaxLengthRatioApproximation = 1.0;                  ///< ratio of relax length to cap length approximation
    double scaffoldingZeroPlaneZ = 0.0;                          ///< Height of the flat membrane plane around the Gag cap
    std::string scaffoldingFileName = "";
    bool isGagScaffoldingEnergyIncluded = false;                 ///< Whether to include Gag-specific internal scaffold energy
    std::string gagReferenceStateFileName = "";                  ///< Extracted selected/json/*.dat file used as Gag reference geometry
    std::string gagReactionFileName = "";                        ///< NERDSS .inp file used to identify Gag reaction blocks
    bool isIdealizedProteinLatticeEnergyIncluded = false;        ///< Whether to include idealized per-instance protein lattice energy
    std::string idealizedProteinLatticeFileName = "";            ///< JSON / .dat file containing the initialized idealized lattice state
    double gagKsigma = 0.0;                                      ///< Spring constant for COM-COM bond length term
    double gagKtheta = 0.0;                                      ///< Spring constant for COM-COM bond angle term
    double gagKphi = 0.0;                                        ///< Spring constant for COM-graph torsion term mapped from phi
    double gagKomega = 0.0;                                      ///< Spring constant for COM-graph torsion term mapped from omega
    double gagFiniteDifferenceStep = 1.0e-4;                     ///< Step size used for scaffold finite-difference gradients
    double gagPropagationStepSize = 1.0e-6;                      ///< Step size used when propagating the scaffold with Gag forces
    int gagPreRelaxSteps = 200;                                  ///< Number of Gag-only pre-relaxation steps before membrane coupling
    

    // dynamics
    double timeStep = 0.1;   // us
    double diffConst = 0.01; // um^2/us
    double KBT = 4.17;       // 1KbT = 4.17 pN.nm
    unsigned int randomSeed = 42; ///< random seed used in normal distribution
    bool fdtConsistentSurfaceUpdate = true; ///< Map the nodal force onto the limit-surface DOFs before the Brownian step
    bool integratePeriodicDuplicates = false; ///< Legacy: give the fourth-ring periodic duplicates their own Brownian kick before overwriting them
    bool surfacepointOutput = true; ///< whether to output surface point file

    // Monte Carlo edge flips -- in-plane fluidity. See docs/edge_flip_plan.md.
    /// Run the Metropolis flip sweep alongside the Brownian step.
    bool edgeFlipEnabled = false;
    /**
     * @brief Attempts per edge per microsecond: the fluidity knob, nu.
     *
     * The number of attempts in a sweep is drawn from a Poisson distribution
     * with mean nu * timeStep * edgeFlipInterval * (flippable edges), so
     * halving the time step halves the attempts and the physical rate per edge
     * is unchanged. That independence is the whole reason the count is drawn
     * rather than fixed.
     *
     * It is a physical quantity, not a numerical one: the accepted-flip rate
     * sets the membrane's in-plane viscosity and the diffusion of a vertex. A
     * vertex at lFace = 5 nm stands for of order a hundred lipids, and with a
     * lipid diffusion constant of 1-10 nm^2/us such a patch exchanges a
     * neighbour every 0.6-6 us, which puts nu in the range 0.1-1 per
     * microsecond. WP6 calibrates it against a measured neighbour-survival
     * time rather than leaving it at a guess.
     */
    double edgeFlipAttemptRate = 0.5;
    /// Steps between sweeps. The attempt count scales with it, so this trades
    /// sweep frequency against sweep size at a fixed physical rate.
    int edgeFlipInterval = 1;
    /// Valences a flip may leave behind, clamped to what the patch tables
    /// support. Narrowing the range cuts the cost of a fluid mesh, at the
    /// price of refusing flips a dynamically triangulated surface would allow.
    int edgeFlipMinValence = 4;
    int edgeFlipMaxValence = 8;

    // Fluid-mode dynamics. See docs/edge_flip_plan.md section 3.7.
    /**
     * @brief Let vertices move in the membrane plane, not only along z.
     *
     * The Brownian step multiplies the x and y displacement by zero. That is
     * defensible for a solid sheet whose triangulation cannot rearrange -- the
     * in-plane degrees of freedom have nowhere useful to go -- and indefensible
     * for a fluid one, where in-plane motion is half of what fluidity means.
     */
    bool inPlaneDynamicsEnabled = false;

    /**
     * @brief Replace the reference-length regularization with an edge spring.
     *
     * The regularization term measures each face's edges against the same
     * face's edges in coordRef. That is a solid's memory of its reference
     * configuration: meaningless for a fluid, and undefined for an edge a flip
     * has just created, which never had a reference length. The spring
     *
     *     E = (k / 2) * sum over edges of (l - l0)^2
     *
     * depends only on the edges that exist now, so it survives a flip. The two
     * are alternatives, not additions: enabling this disables the other.
     */
    bool edgeSpringEnabled = false;
    double edgeSpringConstant = 83.4;   ///< k, in pN/nm.
    /// l0. Negative means "use lFace", the target edge length of the mesh.
    double edgeSpringRestLength = -1.0;

    /**
     * @brief The shape of that term: "flat" (default) or "harmonic".
     *
     * "harmonic" is the spring above, `(k/2)(l - l0)^2`, and it does not work
     * for a fluid membrane. A flip on a rhombus of two equilateral triangles
     * replaces the short diagonal by the long one, so it has to climb
     * `(k/2)(sqrt(3) - 1)^2 l0^2` whatever the rest of the energy says. Two
     * requirements then pull `k` in opposite directions: flips need a barrier
     * of a few kT, and a triangulation that does not degenerate needs a bond
     * fluctuation `sqrt(kT/k)` well below `l0`. At `l0 = 5 nm` and room
     * temperature the first wants `k <= 3.1 pN/nm` and the second wants
     * `k >= 16.7`. The window is empty, and the measurements are in
     * docs/edge_flip_plan.md work package 5.
     *
     * "flat" is the tether every dynamically triangulated surface model uses:
     * zero inside an allowed range, a quadratic wall outside it.
     *
     * ```text
     *     E = (k/2)(l_min - l)^2   l < l_min
     *         0                    l_min <= l <= l_max
     *         (k/2)(l - l_max)^2   l > l_max
     * ```
     *
     * A flip that leaves every edge inside the range costs nothing, so the
     * barrier and the shape constraint stop competing: `k` can be as stiff as
     * the walls need to be, and the flip is decided by the bending energy and
     * the constraints -- which is the physics the sweep exists to sample.
     */
    std::string edgeTetherShape = "flat";
    /**
     * @brief The allowed range, as multiples of the rest length.
     *
     * Squeezed between two requirements that nearly meet.
     *
     * The upper bound must exceed `sqrt(3) = 1.733`: a flip of an equilateral
     * rhombus of side `l` produces an edge of `l sqrt(3)`, so a wall below
     * that forbids exactly the move this exists to permit.
     *
     * The lower bound is a matter of mesh quality, and the range as a whole
     * follows the dynamically triangulated surface literature, whose tether
     * ratio is 1.68-1.73: the narrower the range, the fewer thin triangles a
     * fluid mesh can form, at the cost of acceptance (31% here against 42%
     * at `[0.6, 1.8]`) and a few per cent of bending energy.
     *
     * An earlier version of this comment claimed the range had to be narrow
     * for the mesh to be stationary at all. That was measured with the tether
     * energy summed over every edge, ghost band included, and the ghost
     * band's lattice edges stretch without limit as the interior mixes -- see
     * edge_carries_tether(). Restricted to the membrane, the interior is
     * stationary at either range: docs/fluidity_results.md section 2.
     */
    double edgeTetherMinRatio = 0.95;
    double edgeTetherMaxRatio = 1.75;

    /**
     * @brief How the limit surface and the control net are converted.
     *
     * "dense" builds the whole mask as an N x N matrix and inverts it once,
     * which is what this tree has always done. It costs O(N^2) memory, does
     * not thread, and cannot survive a connectivity change -- an edge flip
     * changes four of its rows and would need the inverse rebuilt from
     * scratch.
     *
     * "iterative" holds the mask sparsely and solves instead. See
     * include/dynamics/Surface_solver.hpp for why both directions reduce to
     * one symmetric positive definite system.
     */
    std::string surfaceSolver = "dense";

    /**
     * @brief Multiplies the per-valence subdivision depth of irregular patches.
     *
     * A face with an extraordinary corner is evaluated by recursing down a
     * chain of regular children, and the depth of that chain is what the patch
     * costs. recommended_irregular_depth() picks a depth per valence for a
     * bending-energy tail below 1e-4 relative, which is the right target for a
     * minimization whose answer is a single converged shape.
     *
     * A Brownian run is not that. It lives in thermal noise several orders
     * above 1e-4, and a fluid mesh is mostly irregular -- measured at WP1 as
     * 36x the cost of an all-regular mesh -- so the depth is the largest lever
     * on the cost of fluidity. This scales every valence's depth together, so
     * a run can trade patch accuracy for throughput and WP6 can measure what
     * the spectrum actually needs.
     *
     * 1.0 is the converged depth and reproduces every existing run exactly.
     * Below 1.0 is cheaper and less accurate; above 1.0 costs build time and
     * memory for accuracy the sampling noise hides.
     */
    double irregularPatchDepthScale = 1.0;

    // thermal fluctuation / annealing for equilibrium searches
    bool thermalFluctuationEnabled = false;           ///< Enable Metropolis thermal trial moves during minimization
    bool thermalFluctuationPureMMC = false;           ///< Run pure Metropolis Monte Carlo trial moves without NCG
    int thermalFluctuationInterval = 50;              ///< Iteration interval between thermal trial moves
    double thermalFluctuationTemperatureKelvin = 298.0;    ///< Effective thermodynamic temperature in Kelvin
    double thermalFluctuationMinTemperatureKelvin = 298.0; ///< Lower bound for annealed temperature in Kelvin
    double thermalFluctuationCoolingRate = 1.0;       ///< Multiplicative cooling rate after each thermal trial
    double thermalFluctuationStepScale = 0.02;        ///< Gaussian displacement std. dev. as a fraction of lFace

    // insertion
    bool isInsertionIncluded = false;
    std::vector<std::vector<int>> insertionPatch;

    Energy energy;
    Energy energyPrev;
    DeformationCount deformationCount;
};

/**
 * @brief Overloaded stream insertion operator to output the Param struct in a readable format
 *
 * @param os The output stream
 * @param param The Param struct to be outputted
 * @return ostream& A reference to the output stream 
 */
std::ostream& operator<<(std::ostream& os, const Param& param);
