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
