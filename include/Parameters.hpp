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
 *
 * The first three are the global modes: the flat-sheet generator lays out
 * ghost rings by grid position and every boundary rule in the code indexes
 * the grid directly. Mixed is the per-vertex mode, where the rule is carried
 * by each vertex (VertexType) and the grid is never consulted -- which is what
 * a mesh with different boundaries on different sides, or an imported mesh
 * with no grid at all, needs. See docs/mixed_boundary_conditions.md.
 */
enum class BoundaryType
{
    /** Fixed boundary condition. */
    Fixed,

    /** Periodic boundary condition. */
    Periodic,

    /** Free boundary condition. */
    Free,

    /**
     * @brief The boundary condition is a property of each vertex.
     *
     * Every vertex is Free, Fixed or a Periodic image of another vertex
     * (VertexType). A periodic image is not a coordinate: it follows its
     * source at a fixed offset and the force it accumulates is folded onto the
     * source, so the periodic energy's gradient is exact at the seam. Faces
     * that duplicate a physical face carry no energy. The generated flat sheet
     * takes its types from boundaryConditionX and boundaryConditionY; an
     * imported mesh takes them from its vertex file.
     */
    Mixed
};

/// "Fixed", "Periodic", "Free" or "Mixed", as the parameter file spells them.
const char *boundary_type_name(BoundaryType type);

/**
 * @brief Parse a boundary type as the parameter file spells it, in any letter case.
 * @return false if the text names no type; @p type is then left untouched.
 */
bool parse_boundary_type(const std::string &text, BoundaryType &type);

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
    double kSpring = 0.0;       ///< Spring constant for insertion zones (K)
    bool setRelaxAreaToDefault = false; ///< true to set area0 equal to area of starting config
    double area0 = 0.0;         ///< Target area for membrane (S0)
    double area = 0.0;          ///< Total area of the membrane (S)
    double vol0 = 0.0;        ///< Target volume for membrane (V0)
    double vol = 0.0;         ///< Total volume of the membrane (V)
    double insertCurv = 0.0;  ///< Spontaneous curvature of insertions (C0)
    double spontCurv = 0.0;   ///< Spontaneous curvature of membrane (c0)

    // membrane size and axes division
    double sideX = 100.0;                           ///< X-axis length for flat membrane
    double sideY = 100.0;                           ///< Y-axis length for flat membrane
    double radius = 25.0;                          ///< Radius for spherical membrane
    double lFace = 5.0;                           ///< lFace
    int nFaceX = -1;                        ///< Number of faces (edges) along X axis for flat membrane
    int nFaceY = -1;                        ///< Number of faces (edges) along Y axis for flat membrane
    double dFaceX = 0.0;                          ///< Initial actual face side length along X axis for flat membrane
    double dFaceY = 0.0;                          ///< Initial actual face side length along Y axis for flat membrane
    double meanL = 0.0;                           ///< Mean length of edges after subdivision
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
     * CPU. That includes a fluid run -- edgeFlipEnabled with faces of several
     * extraordinary corners -- since docs/edge_flip_plan.md WP8; the fluid
     * mesh-quality terms (tether, triangle shape, crease wall) are evaluated
     * on the host after the device returns on either backend. See
     * docs/cuda_implementation.md for what the GPU was measured to gain and
     * for the test to run before trusting a new machine.
     */
    std::string forceBackend = "cpu";
    double elementTriangleArea0 = 0.0;            ///< Target area for individual triangles

    // gauss quadrature
    int gaussQuadratureN = 2;    ///< Number of Gaussian quadrature points to use
    Matrix VWU;                  ///< (N,3) matrix of vertex coordinates and weights
    Matrix gaussQuadratureCoeff; ///< (N,1) matrix of Gaussian quadrature coefficients

    // shape function
    std::vector<Matrix> shapeFunctions; ///< List of shape functions for each triangle

    // boundary conditions
    BoundaryType boundaryCondition = BoundaryType::Periodic; ///< Type of boundary condition ("Fixed", "Periodic", "Free", "Mixed")

    /**
     * @brief The boundary of each axis of the generated flat sheet under
     * BoundaryType::Mixed: Periodic, Free or Fixed.
     *
     * Periodic lays out the same band the global Periodic mode does -- three
     * rings of images and one ring of duplicates on each side, so every real
     * face has a complete one-ring -- but as per-vertex images whose forces
     * fold onto their sources. Fixed clamps the outermost fixedBoundaryRings
     * rings of vertices. Free leaves the edge open: the outermost ring of
     * faces has no complete one-ring and so no limit surface, and the edge
     * vertices move under the forces of the faces one ring in. Ignored by the
     * global modes and by an imported mesh, whose types come from its file.
     */
    BoundaryType boundaryConditionX = BoundaryType::Periodic;
    BoundaryType boundaryConditionY = BoundaryType::Periodic;

    /**
     * @brief How many rings of vertices a Fixed side of the generated sheet
     * clamps.
     *
     * One ring pins the control points on the edge and leaves the slope free
     * (a simply supported edge); two rings pin the tangent as well (a clamped
     * edge). The outermost ring of faces carries no energy either way,
     * because its one-ring is incomplete.
     */
    int fixedBoundaryRings = 1;

    /**
     * @brief Load the mesh from files instead of generating the flat sheet.
     *
     * Both must be set together. The vertex file is "x, y, z" per line with
     * an optional type and mirror column, the faces file "v0, v1, v2" with an
     * optional copy flag; see io.hpp for the format and
     * export_mesh_to_vertices_faces() for the writer.
     */
    std::string meshVerticesFile = "";
    std::string meshFacesFile = "";

    // optimization methods
    bool usingNCG = true;      ///< Whether to use nonlinear conjugate gradient method
    bool isNCGstuck = false; ///< Whether NCG method has gotten stuck

    // convergence criteria
    double deltaEnergyConverge = 1e-5; ///< Convergence criteria for total energy
    double deltaForceScaleConverge = 1e-5; ///< Convergence criteria for max force scale

    // regularization
    /**
     * @brief Shape threshold of the R-pi adaptive regularization (input key gamaShape).
     *
     * A face whose side lengths a, b, c scatter about their mean l by
     * ((a-l)^2 + (b-l)^2 + (c-l)^2) / l^2 > gamaShape counts as shape-deformed
     * and is pulled towards the equilateral triangle of its own area instead
     * of towards its reference side lengths. Read only when usingRpi is true.
     */
    double gamaShape = 0.2;
    /**
     * @brief Area threshold of the R-pi adaptive regularization (input key gamaArea).
     *
     * A face that is not shape-deformed and whose area differs from its
     * reference area by this fraction or more counts as area-deformed and is
     * pulled towards the equilateral triangle of the reference area. Read only
     * when usingRpi is true.
     */
    double gamaArea = 0.2;
    /**
     * @brief Apply one global Loop refinement to the control mesh at setup.
     *
     * Off by default. Only needed for meshes whose extraordinary vertices are
     * adjacent, which the valence histogram reports at setup. It roughly
     * quadruples the vertex count and changes the dynamical degrees of
     * freedom, so enabling it rebaselines the run.
     */
    bool isPreRefinementEnabled = false;

    /**
     * @brief Whether the regularization uses its R-pi adaptive cases (input key usingRpi).
     *
     * true: a face past gamaShape or gamaArea is regularized towards an
     * equilateral triangle, every other face towards its reference side
     * lengths. false: every face is regularized towards its reference side
     * lengths, so the term and its force vanish wherever the reference was
     * just refreshed. A setting of the energy: nothing changes it during a run.
     */
    bool usingRpi = true;

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
    /**
     * @brief Substeps of the mesh-quality force per Brownian step (prototype).
     *
     * 1, the default, is the single Euler step every run has used. n > 1
     * splits the edge tether, altitude floor and crease wall -- the stiff,
     * cheap, local part of the drive -- into n Euler substeps of timeStep / n:
     * the first rides along with the bending, area and noise step from the
     * step's starting force, the other n - 1 re-evaluate only those terms
     * (DynamicModel::mesh_force_substeps()). The noise is drawn once per full
     * step, so the bending modes see exactly the shipped update; only the
     * walls are integrated more finely. Requires fdtConsistentSurfaceUpdate
     * and edgeSpringEnabled, soft walls and the Brownian step.
     */
    int meshForceSubsteps = 1;
    /**
     * @brief Where the substepped step draws its noise: "once" or "split".
     *
     * "split" (the default) draws 2 D timeStep / n on every substep: their
     * sum is still one Gaussian of variance 2 D timeStep, so the bending modes
     * see the shipped step in law, while the walls get a true Euler-Maruyama
     * step of timeStep / n, and the stable step is then timeStep / n <= ~0.003
     * us at lFace 5 up to timeStep 0.02. "once" keeps the shipped single draw
     * on the first substep; the other n - 1 push every wall violation back
     * before the next frame, so the walls are sampled too cold and it needs
     * more substeps to be stable. Read only at meshForceSubsteps > 1.
     */
    std::string meshForceSubstepNoise = "split";
    bool integratePeriodicDuplicates = false; ///< Legacy: give the fourth-ring periodic duplicates their own Brownian kick before overwriting them
    bool surfacepointOutput = true; ///< whether to output surface point file
    /**
     * @brief Write the model's local quantities beside every trajectory frame.
     *
     * Two more trajectory files, one line per frame like meshpoint<prefix>.csv:
     *
     *   faceenergy<prefix>.csv    per face, in face order: bending energy,
     *                             mesh-term energy (tether + altitude floor +
     *                             crease wall, or the reference-length
     *                             regularization), mean curvature, and the
     *                             patch's area (pN.nm, pN.nm, 1/nm, nm^2).
     *   vertexforce<prefix>.csv   per vertex, in vertex order: the bending,
     *                             area-constraint, mesh-term and bead forces,
     *                             three components each (pN).
     *
     * The numbers are the ones the last force evaluation stored, which is
     * the configuration one step before the frame's: at a frame every
     * hundred steps that is invisible, and it keeps the per-face rows on the
     * connectivity the frame's face file describes. The area constraint's
     * energy is global and has no per-face value; the tension it exerts is
     * uniform, and its force per vertex is in the second force block. Off by
     * default. See docs/bead_push_benchmark.md section 8.
     */
    bool localOutput = false;

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
    /**
     * @brief Valences a flip may leave behind, clamped to what the patch
     * tables support (4 to 8).
     *
     * 5 to 7, and not the tables' full 4 to 8, because of geometry. A flat
     * vertex of valence N with legs of length l needs opposite edges of
     * `2 l sin(pi/N)`: with the edges a fluid run actually has, about
     * 1.1 lFace, that is 4.8 nm at valence 7 and 4.2 at valence 8 on a 5 nm
     * mesh, against a tether wall at 4.75. A valence-8 vertex cannot flatten,
     * its surplus angle buckles the neighbourhood, and with the limit-surface
     * bending energy indifferent to a crease in the control net the buckle
     * becomes a flap folded 180 degrees onto its neighbour, whose limit
     * surface pinches and blows the integrator up. Measured (WP7): every
     * fold in the first triangle-shape run sat at a valence-8 or valence-4
     * vertex, and restricting flips to 5-7 gave zero creases where 4-8 gave
     * a dozen. Widening the tether instead makes it worse, because short
     * edges let the bending force crush triangles through the altitude wall.
     *
     * The price is fluidity: fewer flips are admissible, and neighbour
     * survival at a given attempt rate decays more slowly.
     */
    int edgeFlipMinValence = 5;
    int edgeFlipMaxValence = 7;

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
     * @brief Hold the mean height of the integrated vertices where it started.
     *
     * A sheet that is periodic along both axes is held by nothing in z: a
     * bead that pushes on it, or pulls it, moves the whole sheet, at the
     * speed the total force over N kT / D gives, and a push of fixed length
     * is then a push of unknown depth. The surrounding membrane the tile
     * stands for does hold it. With this on, the Brownian step subtracts the
     * mean z displacement of the vertices it integrates, which is the same
     * as spreading the opposite of the net vertical force evenly over them:
     * a uniform counter-force, so the shape and the fluctuations of the
     * tile are untouched and only the drift is removed. Off by default; it
     * does nothing a clamped sheet can notice, and it changes a thermal
     * periodic run by the one Brownian mode it removes.
     */
    bool meanHeightPinned = false;

    /**
     * @brief Truncate the Brownian noise of the vertex step at |xi| < this.
     *
     * The step draws one standard normal xi per coordinate and moves the
     * limit surface by sqrt(2 D dt) xi. A positive value here redraws every
     * xi that falls outside (-noiseTruncation, noiseTruncation), so no single
     * kick exceeds that many standard deviations. 0, the default, is the
     * untruncated Gaussian, and the run is byte-identical to one without the
     * key. The parameter file accepts 0 or a value of at least 0.1.
     *
     * Only the vertex step is affected: the elastic box, the lipid
     * reservoir and the flip sweep draw their own numbers.
     *
     * It does not raise the stable time step: the node cross-over that ends
     * a run at too large a step is the explicit overshoot of the stiff mesh
     * terms, not a large kick, and with the variance held at 1 the divergence
     * rate does not change measurably. See docs/resolution_limits.md
     * section 2.5.
     */
    double noiseTruncation = 0.0;

    /**
     * @brief Divide the truncated noise by its own standard deviation.
     *
     * Truncation alone lowers the variance of xi to
     * truncated_normal_variance(noiseTruncation) -- 0.774 at 2 -- and the
     * step then samples the membrane at that fraction of KBT while its
     * relaxation rates stay where they were, which breaks the
     * fluctuation-dissipation relation between the spectrum and the rates.
     * true (the default) restores unit variance, and with it the spectrum
     * and the rates of the Gaussian step: under Euler-Maruyama both depend on
     * the noise only through its variance. Set false only to study that
     * effect. Ignored when noiseTruncation is 0.
     */
    bool noiseTruncationRenormalize = true;

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
     * @brief Bound triangle shape, not only edge length.
     *
     * The tether bounds every edge and nothing else, and a fluid mesh under
     * in-plane motion and flips fills with slivers inside those bounds: on
     * the 100 nm sheet, interior triangles 0.1 nm tall at an angle of one
     * degree in every frame, with every edge inside its walls. A sliver's
     * normal turns through tens of degrees under a single Brownian kick, its
     * limit-surface patch self-intersects, and the bending force on it is
     * not finite -- every long fluid run ended that way, after 1e4 to 1e5
     * steps. A Monte Carlo model never takes that step; an explicit Brownian
     * one has nothing that refuses it, so the Hamiltonian has to.
     *
     * This term charges a face for each of its three altitudes -- the
     * distance from a corner to the opposite edge, `h_i = 2A / l_i`, which is
     * exactly the quantity that vanishes in a sliver -- that falls below a
     * floor:
     *
     * ```text
     *     E = (k / 2) sum over faces, corners  max(0, h0 - h_i)^2,
     *     h0 = triangleShapeMinAltitudeRatio * lFace
     * ```
     *
     * Zero for any healthy triangle, so it adds no tension and does nothing
     * inside the allowed region, like the tether. It is the second half of
     * the fluid mesh-quality term and requires the tether; setup refuses it
     * without.
     */
    bool triangleShapeEnabled = false;
    /**
     * @brief The altitude floor, as a fraction of lFace.
     *
     * The lattice's altitude is 0.866 lFace. A flip of an equilateral rhombus
     * makes two triangles of altitude exactly 0.5 lFace, so the floor must
     * stay below 0.5 or it forbids the move the model exists to permit; the
     * dynamically triangulated surface literature's tether ratio of 1.7 puts
     * the thinnest allowed triangle at about 0.5 of the minimum edge. 0.4
     * leaves the flip a margin and puts the wall at 2 nm on a 5 nm mesh,
     * where a 0.05 nm Brownian kick turns a normal by under two degrees.
     */
    double triangleShapeMinAltitudeRatio = 0.4;
    double triangleShapeConstant = 83.4;   ///< k, in pN/nm; the tether's stiffness by default.

    /**
     * @brief Bound the crease between adjacent faces.
     *
     * The third mesh-quality term, and the one a dynamically triangulated
     * surface gets for free: its bending energy lives on the control net,
     * so a face folded onto its neighbour costs the bending modulus outright.
     * SLIMED's bending energy lives on the limit surface, which smooths a
     * crease in the control net away until the net is two layers deep and
     * the limit surface pinches. Measured (WP7): with slivers removed by the
     * altitude floor, the sheet still folded at step 82 000 through creases
     * of 177-180 degrees that built up from step 20 000 at valence-8
     * vertices, which the tether's lower wall leaves unable to flatten.
     * Restricting flips to valences 5-7 cut those creases to three in
     * 60 000 steps; this term is what forbids them.
     *
     * For each interior edge, with `c = n1 . n2` the cosine of the angle
     * between its two faces' unit normals,
     *
     * ```text
     *     E = (k / 2) max(0, cos(theta_max) - c)^2
     * ```
     *
     * Zero while the faces are within `theta_max` of coplanar, quadratic in
     * the cosine beyond it, and smooth everywhere -- a wall on the angle
     * itself would have a singular gradient at a full fold, which is where
     * it matters most. Requires the tether, like the shape term.
     */
    bool creaseWallEnabled = false;
    /**
     * @brief The angle between adjacent face normals beyond which the wall
     * starts, in degrees.
     *
     * Thermal undulations on this mesh put adjacent normals a few degrees
     * apart and a fluid run's worst healthy edge under 45; a flap is at 180.
     * 60 leaves the physics alone and puts a full fold at
     * `(k/2)(1 + 1/2)^2 = 1.125 k`.
     */
    double creaseWallAngle = 60.0;
    /**
     * @brief k, in pN.nm (the energy is in the cosine, so the constant carries
     * the units).
     *
     * 500 makes a right-angle crease cost 15 kT and a full fold 135 kT at
     * room temperature: the fold a frustrated vertex drives toward is
     * forbidden, and the crease a thermal kick makes is not noticed.
     */
    double creaseWallConstant = 500.0;

    // Monte Carlo vertex moves and hard walls. See docs/monte_carlo_vertex_moves.md.
    /**
     * @brief How the mesh-quality walls are enforced: "soft" (default) or "hard".
     *
     * "soft" is the three quadratic walls above -- the tether, the altitude
     * floor and the crease wall -- as energies with forces, which is what an
     * explicit Brownian step needs: it has no way to refuse a step, so the
     * Hamiltonian has to push back. The price is a stiffness that caps the
     * time step, and walls that are leaned on rather than respected (a fluid
     * run stores 100-500 pN.nm in them and its shortest edge sits 0.8 nm
     * inside the lower wall).
     *
     * "hard" turns the same three walls, at the same positions, into
     * constraints: they carry no energy and no force, and a Monte Carlo move
     * -- a vertex move or an edge flip -- that would put an edge outside the
     * tether range, an altitude under the floor or a crease past its angle
     * is rejected outright. That is how the dynamically triangulated surface
     * models hold their meshes together. The sampled distribution is then
     * exp(-E/kT) of the membrane's own energy, restricted to the admissible
     * meshes, with nothing stored in the walls and no stiffness to resolve.
     *
     * Which walls exist is still decided by edgeSpringEnabled (required),
     * triangleShapeEnabled and creaseWallEnabled; their stiffness constants
     * are unused. An explicit Brownian step cannot respect a hard wall, so
     * "hard" requires brownianStepEnabled = false and vertexMoveEnabled =
     * true, and the starting mesh must satisfy every wall.
     */
    std::string meshWallType = "soft";

    /**
     * @brief Run a Metropolis sweep of single-vertex moves every
     * vertexMoveInterval iterations.
     *
     * Each attempt picks an independent vertex uniformly at random, displaces
     * its control point by a vector drawn uniformly from the cube
     * [-s, s]^3 with s = vertexMoveStepRatio * lFace (along z only when
     * inPlaneDynamicsEnabled is false), carries its periodic copies along,
     * and accepts with min(1, exp(-dE/kT)), where dE is the exact change of
     * the energy the run samples: bending, the area and volume constraints,
     * and the soft walls if those are in force. Under meshWallType = hard a
     * move that breaks a wall is rejected before its energy is looked at.
     * The proposal is symmetric, so the sweep leaves exp(-E/kT) invariant
     * and composes with the flip sweep and the Brownian step exactly as
     * those compose with each other.
     *
     * Off by default; a run with it off is unchanged to the byte.
     */
    bool vertexMoveEnabled = false;
    /// Half-width of the trial cube, as a fraction of lFace. Tune it for an
    /// acceptance of 30-50%; it is not adapted during the run, because an
    /// adaptive proposal does not satisfy detailed balance.
    double vertexMoveStepRatio = 0.05;
    /// Mean attempts per independent vertex per sweep. 1 is the "Monte Carlo
    /// step" of the triangulated-surface literature.
    double vertexMoveAttemptsPerVertex = 1.0;
    /// Iterations between sweeps.
    int vertexMoveInterval = 1;

    /**
     * @brief Take the Brownian step. On by default.
     *
     * false makes the run a pure Monte Carlo one: an iteration is then a
     * vertex-move sweep and a flip sweep, there is no time in it, and
     * timeStep only converts edgeFlipAttemptRate into attempts per edge per
     * iteration (rate * timeStep * edgeFlipInterval). Requires
     * vertexMoveEnabled = true, since nothing else would move a vertex.
     */
    bool brownianStepEnabled = true;

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

    // External force on a region of the membrane. See docs/external_force.md.
    /**
     * @brief Push a region of the membrane with a constant force per vertex.
     *
     * Every independent vertex whose initial (x, y) lies inside the region
     * receives the force (externalForceX, externalForceY, externalForceZ),
     * in pN, at every Brownian step from externalForceStartIteration up to
     * but not including externalForceStopIteration. The force acts on the
     * limit-surface point the step displaces, so a vertex inside the region
     * moves by (D dt / kT) F per step on top of what the membrane's own
     * forces do to it. With inPlaneDynamicsEnabled = false only the z
     * component acts, as for every other force. Fixed vertices, ghosts and
     * periodic images are never pushed: they are not integrated.
     *
     * A uniform pressure p over the region corresponds to a force per vertex
     * of about p * (sqrt(3)/2) lFace^2 -- 21.7 nm^2 per vertex on a 5 nm mesh.
     *
     * Off by default; a run with it off is unchanged to the byte.
     */
    bool externalForceEnabled = false;
    double externalForceX = 0.0;    ///< pN per vertex
    double externalForceY = 0.0;    ///< pN per vertex
    double externalForceZ = 0.0;    ///< pN per vertex; positive lifts the membrane
    /// "disc" (radius externalForceRadius about the centre) or "rectangle"
    /// (half-widths externalForceHalfWidthX/Y about the centre).
    std::string externalForceShape = "disc";
    /// Centre of the region in the sheet's plane, nm. The generated sheet is
    /// centred on the origin, so (0, 0) is the middle of the membrane.
    double externalForceCenterX = 0.0;
    double externalForceCenterY = 0.0;
    double externalForceRadius = 20.0;       ///< disc radius, nm
    double externalForceHalfWidthX = 20.0;   ///< rectangle half-width along x, nm
    double externalForceHalfWidthY = 20.0;   ///< rectangle half-width along y, nm
    /// The window of iterations the force acts in: [start, stop). A negative
    /// stop means the force never switches off.
    int externalForceStartIteration = 0;
    int externalForceStopIteration = -1;

    // A rigid spherical bead pushed through the membrane. See docs/bead_push_benchmark.md.
    /**
     * @brief A rigid sphere the membrane cannot enter, moving at constant velocity.
     *
     * Every independent vertex closer to the bead's centre than beadRadius
     * is pushed out radially by a harmonic wall,
     *
     * ```text
     *     E = (k / 2) sum over vertices  max(0, R - |x_i - c|)^2,   k = beadStiffness,
     * ```
     *
     * evaluated on the control points. The potential depends on positions
     * only, so an edge flip -- which moves no vertex -- leaves it unchanged
     * and the Metropolis sweep stays exact with the bead present. The force
     * is a nodal force and reaches the Brownian step through the same map as
     * the membrane's own (M^-T under fdtConsistentSurfaceUpdate). Fixed
     * vertices, ghosts and periodic images are not integrated and are not
     * pushed.
     *
     * The centre moves as c = c0 + v (t - t_start) between beadStartIteration
     * and beadStopIteration and stays put outside that window, so a run can
     * approach, push through, and then hold; a negative stop means it never
     * stops moving. The wall acts for the whole run wherever the bead is.
     *
     * Per trajectory frame the dynamics driver writes <prefix>bead.csv: the
     * centre, the wall energy, how many vertices are inside the wall and how
     * deep, and the total force the membrane exerts on the bead.
     *
     * Off by default; a run with it off is unchanged to the byte.
     */
    bool beadEnabled = false;
    double beadRadius = 10.0;       ///< R, nm
    double beadStiffness = 200.0;   ///< k, pN/nm
    double beadCenterX = 0.0;       ///< c0, nm
    double beadCenterY = 0.0;
    double beadCenterZ = 0.0;
    double beadVelocityX = 0.0;     ///< v, nm/us
    double beadVelocityY = 0.0;
    double beadVelocityZ = 0.0;
    int beadStartIteration = 0;
    int beadStopIteration = -1;
    /**
     * @brief Adhesion between the bead and the membrane: a well outside the wall.
     *
     * A vertex within beadAdhesionRange of the sphere's surface gains
     *
     * ```text
     *     E = -w s(x),   x = (d - R) / delta,   s(x) = 1 - 3 x^2 + 2 x^3  on [0, 1],
     * ```
     *
     * with s = 1 inside the sphere and 0 beyond the range: a smooth step of
     * depth w = beadAdhesionEnergy per vertex whose pull peaks at 1.5 w / delta.
     * Zero by default, a bead that only pushes. Wrapping a sphere costs
     * bending of 2 kc / R^2 per area plus the tension, so with the lattice's
     * 21.65 nm^2 per vertex a 10 nm bead at kc = 83.4 needs w above about
     * 40 pN.nm (10 kT) to be wrapped rather than tented over. Like the wall,
     * a function of the control points alone, so the flip sweep stays exact.
     */
    double beadAdhesionEnergy = 0.0;   ///< w, pN.nm per vertex
    double beadAdhesionRange = 2.0;    ///< delta, nm

    /**
     * @brief The shape of the bead's potential: "wall" or "harmonic".
     *
     * "wall" is the harmonic wall inside the sphere plus the adhesion well
     * outside it, above: the bead only pushes, and the well is what holds
     * the membrane on. "harmonic" makes the bead sticky. A vertex that comes
     * within beadCaptureRange of the sphere's surface is bonded to it, and a
     * bonded vertex feels the whole harmonic bond
     *
     *     E = (k / 2) (|x_i - c| - R)^2,   k = beadStiffness,
     *
     * which pushes it out when it is inside the sphere and pulls it back
     * when it is outside -- the bead keeps what it touched. A bond breaks
     * when its vertex is farther than beadBreakRange from the surface;
     * negative, the default, never. The bonded set is state the step
     * carries, changed only by refresh_bead_contact(), which a flip never
     * calls: like the wall the energy is a function of the control points
     * alone, so the flip sweep stays exact. beadAdhesionEnergy is ignored
     * under "harmonic". See docs/bead_push_benchmark.md section 9.
     */
    std::string beadPotential = "wall";
    double beadCaptureRange = 2.0; ///< nm outside the surface within which a vertex is bonded
    double beadBreakRange = -1.0;  ///< nm outside the surface beyond which a bond breaks; negative = never

    /**
     * @name The elastic box (Param::elasticBox*)
     *
     * Let the periodic box relax so that the membrane's frame tension sits
     * at a prescribed value, instead of being whatever the fixed projected
     * area happens to impose. The strain of every periodic axis becomes a
     * Brownian coordinate driven by the difference between the target and
     * the tension the membrane actually exerts on its frame (the in-plane
     * virial of the nodal forces), and the whole control net -- images,
     * ghosts and reference coordinates included -- is rescaled with it.
     * Dynamics only (membrane_dynamics); the minimizer ignores it. Off by
     * default. See docs/elastic_box.md.
     * @{
     */
    bool elasticBoxEnabled = false;
    /// tau, pN/nm: the frame tension the box is driven to. Zero is a
    /// stress-free membrane; positive stretches it.
    double elasticBoxTension = 0.0;
    /// m, nm/(pN us): the strain rate per unit tension mismatch,
    /// d(lambda)/dt = m (tau - tau_measured). Intensive, so the same value
    /// means the same relaxation rate whatever the box size: with an area
    /// modulus K the strain relaxes at m n K per us, n periodic axes.
    double elasticBoxMobility = 1.0e-3;
    /// Give the box its own Brownian kick, sqrt(2 mu kT dt), so that the
    /// box strain samples the constant-tension ensemble. false relaxes the
    /// box deterministically toward the target (a Berendsen-like coupling),
    /// which is smoother but is not an ensemble.
    bool elasticBoxThermal = true;
    /// "isotropic": one strain for every periodic axis, driven by the mean
    /// tension, so the box keeps its shape. "anisotropic": each periodic axis
    /// follows its own tension component; on a fluid sheet the aspect ratio
    /// is then a soft mode that wanders.
    std::string elasticBoxCoupling = "isotropic";
    /// The largest |d lambda| one step may apply. A step that wants more is
    /// clipped and counted; the count is reported at the end, and a run
    /// that clips often has a mobility too large for its time step.
    double elasticBoxMaxStrainPerStep = 1.0e-3;
    /// Count the in-plane entropy of the vertices the step integrates in
    /// the plane (inPlaneDynamicsEnabled): the -N kT ln(A_p) of a box move
    /// on coordinates that are rescaled with it, which every constant-tension
    /// scheme on particles carries (the (A'/A)^N of FreeDTS's box move, the
    /// N kT/V of a Langevin piston). With it the prescribed tension is the
    /// one the fluctuation spectrum reads; without it the box's zero is
    /// offset by N kT / A_p, about 0.2 pN/nm on a 100 nm sheet, and a fluid
    /// sheet run at "zero" tension buckles once the box is large enough.
    /// Moot for a solid sheet, whose in-plane coordinates are not integrated.
    bool elasticBoxVertexGasTerm = true;
    /** @} */

    /**
     * @name Lipid flow: the relaxed area as a coordinate (Param::lipidFlow*)
     *
     * The stretching energy's relaxed area becomes a variable that moves
     * with the lipids. Two modes, see docs/lipid_flow_plan.md:
     *
     * - "global": one relaxed area for the whole sheet, A0(t). The energy is
     *   (K / 2 c)(A - A0)^2 with c the relaxed area as built, so with no
     *   reservoir this is exactly the global constraint the tree always had.
     *   A reservoir moves A0 toward the tension it sets: the limit of
     *   infinitely fast lateral flow, cheap, no tension gradients.
     * - "local": a relaxed area S0_f per face, e_f = (K / 2 c_f)(S_f - S0_f)^2,
     *   the face's own tension in the force, and S0_f moving between
     *   neighbouring faces down the tension gradient (tension diffuses at
     *   D_sigma = K M), with the reservoir per face. Finite-rate tension
     *   propagation; costs 6-10 % of run time.
     *
     * "off" leaves the area term untouched.
     * @{
     */
    std::string lipidFlowMode = "off";
    /// D_sigma, nm^2/us: the tension diffusivity K M of the lateral flow.
    /// 0 freezes the lipids (purely local elasticity). Negative takes the
    /// value the vertices' own friction implies, K a_v D / kT.
    double lipidFlowDiffusivity = -1.0;
    /// Explicit substeps of the transport per dynamics step; 0 chooses the
    /// number that keeps the fastest mode's rate x dt at or under 0.5.
    int lipidFlowSubsteps = 0;
    /// Conserved (edge-flux) and reservoir noise at kT, so that the relaxed
    /// areas sample their Boltzmann distribution.
    bool lipidFlowThermal = false;
    /// The form of the local stretching energy. true: (K / 2 abar_f)(S_f - a_f)^2
    /// with abar_f the face's relaxed area as initialized, so the tension
    /// K (S_f - a_f) / abar_f is linear in a_f, the lipids respond to the
    /// tension itself, and thermal fluctuations of a_f leave its mean alone.
    /// false: (K / 2 a_f)(S_f - a_f)^2, the modulus following the lipids
    /// actually present; the tension is then 1/a_f-like, and with the noise
    /// on its mean is shifted by a mesh-dependent kT / a_f.
    bool lipidFlowLinearStretching = true;
    /// Where the membrane exchanges lipids with a reservoir: none, rim (the
    /// energy-carrying faces along an open or clamped boundary) or bulk
    /// (every face).
    std::string lipidReservoirRegion = "none";
    double lipidReservoirTension = 0.0; ///< sigma_res, pN/nm
    double lipidReservoirRate = 0.0;    ///< k, 1/us: d sigma/dt = -k (sigma - sigma_res) on a reservoir face
    /**
     * Local mode with edge flips: what an accepted flip does to the relaxed
     * areas around it. "none" carries each face's relaxed area with the face
     * (exact, no redistribution), and the flip is charged the change of the
     * two faces' stretching energies at the local modulus -- which makes
     * flips rare under a stiff K. "pair" and "patch" make the flip a compound
     * move: flip, then redraw the relaxed areas of the two faces ("pair") or
     * of every face whose patch the flip touched ("patch") from their
     * Boltzmann distribution at the new geometry with their total fixed. The
     * Metropolis ratio then sees the set's stretching energy minimised over
     * that split, and the acceptance returns toward the global constraint's.
     * Exact for sampling (a partial Gibbs step); for the dynamics it is an
     * instant redistribution of lipids over the set at every accepted flip.
     * Needs lipidFlowLinearStretching. With lipidFlowThermal off the redraw
     * is the equal-tension split itself. "auto", the default, is "patch"
     * when edge flips are enabled under the linear form and "none"
     * otherwise; DynamicMesh::finish_dynamic_setup() resolves it.
     */
    std::string lipidFlipResampling = "auto";
    /** @} */

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
