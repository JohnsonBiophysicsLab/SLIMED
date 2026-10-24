/**
 * @file Dynamics.hpp
 * @author Y Ying
 * @brief This file defines dynamic mesh class that
 * inherits mesh with support of time-resolved dynamics.
 *
 * @date 2023-04-17
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <math.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <random>
#include <algorithm>
#include <array>
#include <string>
// model setup
#include "mesh/Mesh.hpp"
#include "mesh/Face.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Gauss_quadrature.hpp"
#include "energy_force/Energy.hpp"
#include "energy_force/Force.hpp"
#include "model/Model.hpp"
#include "dynamics/Surface_solver.hpp"

// matrix math
#include "linalg/Linear_algebra.hpp"
// parameters
#include "Parameters.hpp"

/**
 * @brief A class representing a dynamic mesh.
 *
 * This class inherits from the Mesh class and adds functionality to simulate the mesh dynamically.
 */
class DynamicMesh : public Mesh
{
public:
    Matrix mesh2surface;
    Matrix surface2mesh;
    Matrix matMesh;
    Matrix matSurface;

    /// True where a vertex is a periodic duplicate of another vertex.
    ///
    /// Under BoundaryType::Periodic the fourth ring in from each edge repeats
    /// the ring on the far side, so its coordinates are not independent
    /// degrees of freedom -- postprocess_ghost_periodic() overwrites them from
    /// their partner at the end of every step.  Empty unless the boundary
    /// condition is periodic.
    std::vector<char> isSlavedPeriodic;

    /**
     * @brief Construct a new Dynamic Mesh object. Calls the parent
     * constructor Mesh::Mesh(srcParam).
     *
     * @param srcParam
     */
    DynamicMesh(Param &srcParam);

    /**
     * @brief Overloads the setup_flat in superclass Mesh. Add
     * functionality to assign mesh2surface and surface2mesh.
     *
     */
    void setup_flat();

    /**
     * @brief The same for an imported mesh: Mesh::setup_from_vertices_faces()
     * and then the dynamics' own setup.
     *
     * Both overloads hide Mesh's, so a DynamicMesh loaded through
     * import_mesh_from_vertices_faces() is ready to step, exactly as one built
     * by setup_flat() is.
     */
    void setup_from_vertices_faces(const std::vector<std::vector<double>> &verticesData,
                                   const std::vector<std::vector<int>> &facesData);
    void setup_from_vertices_faces(const std::vector<std::vector<double>> &verticesData,
                                   const std::vector<std::vector<int>> &facesData,
                                   const std::vector<VertexType> &vertexTypes,
                                   const std::vector<int> &mirrorVertices,
                                   const std::vector<char> &faceIsCopy);

    /**
     * @brief What the dynamics adds on top of Mesh's setup, whichever way the
     * mesh was built.
     *
     * The configuration checks a fluid run needs, the dense conversion
     * matrices (unless the iterative solver is in use), the coordinate
     * matrices, and the slaved-vertex map. Called by every setup entry point
     * above.
     */
    void finish_dynamic_setup();

    /**
     * @brief The per-step boundary post-processing for the run's boundary mode.
     *
     * Periodic: postprocess_ghost_periodic(). Mixed:
     * postprocess_periodic_images(). Fixed and Free: nothing.
     */
    void postprocess_boundary();

    /**
     * @brief BoundaryType::Mixed: put every periodic image row of matMesh at
     * its source's row plus its offset.
     *
     * The iterative solver already writes those rows this way, so on the
     * path the driver takes this is a no-op that pins the invariant; it is
     * what keeps a caller that fills matMesh some other way honest.
     */
    void postprocess_periodic_images();

    /**
     * @brief Assign the mesh2surface member.
     *
     * Assign values to mesh2surface matrix that convertes mesh to surface point
     * matrix.
     *
     */
    void assign_mesh2surface();

    /**
     * @brief The sparse, valence-aware form of the same mask.
     *
     * Used when param.surfaceSolver is "iterative". Rebuilt whenever the
     * connectivity changes, which the dense inverse cannot be.
     */
    slimed::SurfaceSolver surfaceSolver;

    /// Build surfaceSolver if it is missing or built for an older topology.
    void ensure_surface_solver();

    /**
     * @name Converting between the control net and the limit surface
     *
     * Three operations the dynamics needs every step, behind one dispatch so
     * that the driver does not care which representation is in use.
     *
     * "dense" reproduces exactly what this tree has always done, including its
     * standing approximation of using M^-1 where M^-T belongs -- the two agree
     * wherever the mask is symmetric, which is every interior valence-6 row.
     * "iterative" is the corrected version: valence-aware, and the force map is
     * the actual transpose.
     * @{
     */
    /// matSurface = M * matMesh.
    void apply_mesh_to_surface();
    /// Solve M * matMesh = matSurface.
    void apply_surface_to_mesh();
    /// Carry a nodal force onto the surface coordinates the step displaces.
    void apply_nodal_force_to_surface(const Matrix &nodalForce, Matrix &surfaceForce);
    /** @} */

    /**
     * @brief Updates the vertices matrix using the values from the vertices vector.
     *
     */
    void update_vertices_mat_with_vector();

    /**
     * @brief Updates the vertices vector using the values from the vertices matrix.
     *
     */
    void update_vertices_vector_with_mat();

    /**
     * @brief Post process ghost vertices in case of periodic boundary condition
     *
     * @note This function does not have check for PBC - it assumes the model
     * is in PBC!
     *
     */
    void postprocess_ghost_periodic();

    /**
     * @brief Fill isSlavedPeriodic from the periodic partner map.
     */
    void mark_slaved_periodic_vertices();

    /**
     * @brief Say at startup whether the flip move can actually move.
     *
     * A run with flips on and a tether too stiff to flip over reports its
     * attempts, accepts none, and looks exactly like a run whose attempt rate
     * is too low. The barrier is known in closed form before the first step,
     * so there is no reason to let the run find out by producing nothing.
     */
    void report_edge_flip_feasibility();

    /**
     * @brief Refuse the Monte Carlo settings that do not make a sampler.
     *
     * Hard walls with the Brownian step on, a run in which nothing moves a
     * vertex, and vertex moves beside an energy term their trial does not
     * difference all run without complaint and sample something other than
     * what was asked for. Builds the vertex-move table when the move is on.
     * See Param::meshWallType and Param::vertexMoveEnabled.
     */
    void check_monte_carlo_configuration();

    /**
     * @name The periodic box (Param::elasticBoxEnabled)
     *
     * What the sheet wraps along, how long the period is, and the affine
     * rescale that stretches it. Detected once by detect_periodic_box() at
     * the end of setup; the lengths and scales are then maintained by
     * apply_box_scale(). See docs/elastic_box.md.
     * @{
     */
    /// Per axis: whether the sheet wraps along it.
    std::array<char, 3> boxPeriodicAxis = {{0, 0, 0}};
    /// Per axis: the period along a periodic axis, nm; zero otherwise.
    std::array<double, 3> boxLength = {{0.0, 0.0, 0.0}};
    /// The frame's width across a single periodic axis, nm: the membrane
    /// area as built divided by the period, which is the strip's width for
    /// a flat strip and the circumference for a tube. Set by
    /// DynamicModel::setup_elastic_box(); unused when two axes wrap.
    double boxWidth = 0.0;
    /// Per axis: the stretch applied since setup, 1 as built.
    std::array<double, 3> boxScale = {{1.0, 1.0, 1.0}};

    /**
     * @brief Read the periodic axes and their periods off the mesh.
     *
     * Under BoundaryType::Periodic both in-plane axes wrap and the period
     * is the tile the grid lays out, (nFace - 6) cells. Under Mixed an axis
     * wraps when some image is offset along it and the period is that
     * offset, which reads an imported mesh as well as the generated sheet.
     * Nothing wraps under Fixed and Free.
     */
    void detect_periodic_box();

    /// How many axes wrap.
    int box_periodic_axis_count() const;

    /// The frame area the tension is referred to: the product of the periods
    /// when two axes wrap, the period times boxWidth when one does. Zero
    /// when none does.
    double box_projected_area() const;

    /**
     * @brief Stretch the whole control net by factor[axis] about the origin.
     *
     * Every vertex's coordinate, previous coordinate and reference
     * coordinate, every row of matMesh, every image offset, and the box
     * lengths and scales -- so that afterwards an image is still its source
     * plus its offset, a ghost of the global Periodic mode is still its
     * partner plus the same lattice translation, and the reference
     * configuration of the solid regularization has stretched with the
     * frame. The surface solver's copy of the offsets is refreshed. Loop
     * subdivision is affine invariant, so this is also the stretch of the
     * limit surface; matSurface is left for the next mesh-to-surface map.
     */
    void apply_box_scale(const std::array<double, 3> &factor);
    /** @} */

protected:
    /**
     * @brief calculate real point relative to the given ghost point (index(real) - index(given)) in periodic boundary condition
     * returns 0 if real point (not on 4th ring) given
     * an arbitrary real point is chosen if given 4th
     *
     * @param i
     * @param n
     * @param m
     * @return int
     */
    int get_relative_pt_periodic(int i, int n, int m);
};

class DynamicModel : public Model
{
public:
    DynamicMesh &mesh; /**< Hides Model::mesh; The mesh object. */

    double randScaleConst;
    double forceScaleConst;
    unsigned int randomSeed;
    /// |xi| cutoff of the vertex noise (Param::noiseTruncation); 0 = Gaussian.
    double noiseCutoff = 0.0;
    /// Factor every truncated xi is multiplied by: 1 / its standard deviation,
    /// or 1 when Param::noiseTruncationRenormalize is false.
    double noiseRescale = 1.0;

    Matrix nodalForce;  /**< -dE/dc_i, assembled per control point, one row per vertex. */
    Matrix driveForce;  /**< The same force expressed in whatever coordinates next_step() advances. */
    std::mt19937 gen;
    std::normal_distribution<> normal_dist;

    /**
     * @brief Constructs a new Model object.
     *
     * @param mesh_ The mesh object to be encapsulated.
     * @param record_ The record object to be encapsulated.
     */
    DynamicModel(DynamicMesh &mesh_, Record &record_);

    /**
     * @brief This code iterates over all vertices in a mesh and displaces
     * them randomly based on a force term and a random term. The displacement
     * is added to the original vertex position.
     *
     * The boundary conditions are checked for each vertex before
     * displacement, and if the vertex is a ghost or boundary,
     * and not part of a periodic boundary condition, then the displacement is set to 0.
     *
     * The displacement amount, along with the vertex index,
     * whether it is a boundary or ghost vertex, and the dimension
     * being displaced (x, y, z) is printed to the console using std::cout.
     */
    /**
     * @name The external force (Param::externalForce*)
     *
     * A constant force per vertex on a region of the membrane, switched on
     * and off by iteration. The region is chosen once, when the model is
     * built, from the vertices' positions then: the sheet's plane does not
     * move, so a vertex is inside or outside for the whole run.
     * @{
     */
    /// One flag per vertex: whether the force acts on it. Empty unless enabled.
    std::vector<char> externalForceMask;
    /// How many vertices the mask selects.
    int externalForceVertexCount = 0;

    /**
     * @brief Fill externalForceMask from the parameters.
     *
     * Only vertices the step integrates are selected: not fixed, not ghost,
     * not a periodic image, not a slaved duplicate. Called by the constructor.
     *
     * @throw std::runtime_error if the shape is not "disc" or "rectangle",
     * if its size is not positive, or if the region selects no vertex at
     * all -- a mis-set centre, not a run worth starting.
     */
    void select_external_force_vertices();

    /// Whether the force acts at this iteration: enabled, and inside
    /// [externalForceStartIteration, externalForceStopIteration).
    bool external_force_active(int iteration) const;

    /// Write the selected vertex indices, one per line, so an analysis can
    /// tell the pushed region from the rest of the membrane.
    void write_external_force_vertices_csv(const std::string &filepath) const;
    /** @} */

    /**
     * @name The rigid bead (Param::bead*)
     *
     * A sphere the membrane cannot enter, moving at constant velocity inside
     * an iteration window. The wall is evaluated on the control points, so it
     * is a nodal force; next_step() adds it to the membrane's own before the
     * map onto the surface.
     * @{
     */
    Matrix beadForce;                                    ///< The wall's force on each control point, one row per vertex.
    std::array<double, 3> beadCenter = {{0.0, 0.0, 0.0}}; ///< Where the last step put the centre.
    double beadEnergy = 0.0;                             ///< The wall energy at the last evaluation, pN.nm.
    int beadContacts = 0;                                ///< Vertices inside the wall at the last evaluation.
    double beadMaxPenetration = 0.0;                     ///< The deepest R - d among them, nm.
    double beadAdhesion = 0.0;                           ///< The adhesion energy at the last evaluation, pN.nm (zero or negative).
    int beadAdhered = 0;                                 ///< Vertices inside the sphere or its adhesion range.
    std::array<double, 3> beadReaction = {{0.0, 0.0, 0.0}}; ///< Force the membrane exerts on the bead: minus the sum of beadForce.
    std::vector<unsigned char> beadBonded;               ///< beadPotential = harmonic: 1 for a vertex bonded to the bead, one entry per vertex.
    int beadBondedCount = 0;                             ///< Bonded vertices at the last evaluation.
    double beadMaxStretch = 0.0;                         ///< The largest d - R among the bonded vertices outside the sphere, nm.

    /**
     * @brief Allocate the force rows and check the bead's parameters.
     *
     * Called by the constructor.
     * @throw std::runtime_error if the bead is enabled with a radius or a
     * stiffness that is not positive.
     */
    void setup_bead();

    /// The centre at this iteration: c0 + v dt (iteration - start), clamped
    /// to the motion window.
    std::array<double, 3> bead_center_at(int iteration) const;

    /// Put the bead where this iteration wants it and evaluate the wall on
    /// the current control points. Called by next_step().
    void update_bead(int iteration);

    /// Evaluate the wall at the current centre on the current control points:
    /// beadForce, beadEnergy, beadContacts, beadMaxPenetration, beadReaction.
    void refresh_bead_contact();

    /// Append one row for this iteration to the bead log, with a header
    /// when the file is new.
    void append_bead_log(const std::string &filepath, int iteration);
    /** @} */

    /**
     * @name The elastic box (Param::elasticBox*)
     *
     * The strain of the periodic axes as a Brownian coordinate of its own,
     * driven by the mismatch between the target tension and the tension the
     * membrane exerts on its frame. See docs/elastic_box.md.
     * @{
     */
    /// mu, 1/(pN nm us): the box coordinate's mobility, fixed at setup from
    /// the initial projected area so that the noise needs no drift correction.
    double boxMobility = 0.0;
    /// Per axis: sum over vertices of coord * F for the terms the step
    /// drives with -- bending, area, the fluid mesh term when the spring is
    /// on, and the bead -- from the last evaluation. pN.nm.
    std::array<double, 3> boxVirial = {{0.0, 0.0, 0.0}};
    /// Per axis: -boxVirial / A_p, the tension the membrane exerts along
    /// that axis, pN/nm. Meaningful on the periodic axes.
    std::array<double, 3> boxTension = {{0.0, 0.0, 0.0}};
    /// The tension the box is driven by: boxTension averaged over the
    /// periodic axes.
    double boxTensionMean = 0.0;
    /// boxTensionMean split by term, for the log: bending, area constraint,
    /// mesh-quality term, bead, and the in-plane vertex gas.
    std::array<double, 5> boxTensionByTerm = {{0.0, 0.0, 0.0, 0.0, 0.0}};
    /// Vertices whose in-plane coordinates the step integrates and the box
    /// rescales: the N of the -N kT ln(A_p) term. Zero on a solid sheet or
    /// with elasticBoxVertexGasTerm off.
    int boxInPlaneVertexCount = 0;
    /// -boxInPlaneVertexCount kT / A_p: the vertex gas's share of every
    /// periodic axis's tension, pN/nm.
    double boxTensionVertexGas = 0.0;
    /// Sum of boxTensionMean over the box steps taken, and their count.
    double boxTensionRunningSum = 0.0;
    long long boxSteps = 0;
    /// Box steps whose strain increment was clipped at elasticBoxMaxStrainPerStep.
    long long boxClippedSteps = 0;

    /**
     * @brief Check the box parameters and fix the mobility. Called by the
     * constructor; does nothing unless the box is enabled.
     * @throw std::runtime_error if the sheet wraps along no axis, the
     *        mobility or the strain cap is not positive, or the coupling is
     *        neither "isotropic" nor "anisotropic".
     */
    void setup_elastic_box();

    /// Read the virials of the last force evaluation into boxVirial,
    /// boxTension, boxTensionMean and boxTensionByTerm.
    void refresh_box_tension();

    /**
     * @brief One Brownian step of the box strain, then the rescale.
     *
     * Called after the vertex step of the same iteration and from the same
     * force evaluation, so that the two together are one Euler step of the
     * joint system. Does nothing unless the box is enabled.
     */
    void update_elastic_box();

    /// Append one row for this iteration to the box log, with a header when
    /// the file is new.
    void append_box_log(const std::string &filepath, int iteration);
    /** @} */

    /**
     * @name Lipid flow (Param::lipidFlow*)
     * @{
     */
    double lipidDiffusivity = 0.0;        ///< D_sigma in force, nm^2/us
    double lipidMobility = 0.0;           ///< M = D_sigma / K, nm^3/(pN us)
    int lipidSubsteps = 1;                ///< explicit substeps per dynamics step
    double lipidReservoirExchanged = 0.0; ///< cumulative relaxed area taken from the reservoir, nm^2
    long long lipidClampedUpdates = 0;    ///< updates that hit the positivity floor
    /// Resolve the diffusivity and the substep count. Called by the constructor.
    void setup_lipid_flow();
    /// One dynamics step of the relaxed areas: lateral flux, reservoir, noise.
    void advance_lipid_flow();
    /// The global mode's step: the sheet's relaxed area toward the reservoir.
    void advance_global_lipid_flow();
    /// One row of scalars for this iteration: totals and tension statistics.
    void append_lipid_log(const std::string &filepath, int iteration);
    /// One line of (a_f, S_f) per face, for a tension map.
    void append_lipid_face_frame(const std::string &filepath);
    /** @} */

    void next_step();

    /**
     * @brief The remaining n - 1 substeps of the mesh-quality force, after
     * next_step() and the map back to the control net (Param::meshForceSubsteps).
     * Each re-evaluates only the fluid mesh terms, moves the surface points by
     * D (timeStep / n) / kT times that force, and maps back. Does nothing at n = 1.
     */
    void mesh_force_substeps();

private:
    bool beadLogHeaderWritten = false;
    bool boxLogHeaderWritten = false;
    bool lipidLogHeaderWritten = false;
};

/*
//boundary condition
int get_relatve_pt_periodic(int i, int n, int m);
void postprocess_ghost_periodic(Param& param, gsl_matrix *verticesOnMesh, std::vector<Vertex>& vertex);
*/