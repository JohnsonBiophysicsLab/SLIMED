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
