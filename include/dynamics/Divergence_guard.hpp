/**
 * @file Divergence_guard.hpp
 * @brief When run_dynamics_flat() should give up on a run that has blown up.
 *
 * A run that diverges does not always reach NaN. diag_fixedbox_dt0.005 (the
 * frozen flat sheet at timeStep 0.005, docs/resolution_limits.md section 4)
 * went from E_total ~ 2e3 to ~ 2e206 pN.nm in a handful of steps, with
 * coordinates ~ 1e48 nm, and then stayed finite for the remaining 80k steps,
 * writing garbage the whole way. The guard stops such a run at the step it
 * goes, on any of three signs:
 *
 *  - the total energy is not finite (the original check);
 *  - a vertex coordinate is not finite, or lies more than kCoordinateReach
 *    starting-mesh sizes outside the starting mesh's bounding box;
 *  - |E_total| exceeds kEnergyMultiple times the energy scale of the start,
 *    max(|E_total(0)|, nVertices * KBT). The floor matters: a flat sheet
 *    starts at E_total ~ 1e-27.
 *
 * Both bounds are orders of magnitude beyond anything a healthy run reaches
 * (the elastic box would need to stretch tenfold), so the guard only reads
 * the mesh and never changes the numbers of a run it lets through.
 */
#pragma once

#include <array>
#include <string>

#include "mesh/Mesh.hpp"

class DivergenceGuard
{
public:
    /// How far outside the starting bounding box a coordinate may go, in
    /// units of the box's largest side.
    static constexpr double kCoordinateReach = 10.0;
    /// How large |E_total| may grow, in units of the starting energy scale.
    static constexpr double kEnergyMultiple = 1.0e6;

    /**
     * @brief Takes the reference from the mesh as it is about to be stepped:
     * its coordinates and its current param.energy.energyTotal.
     */
    explicit DivergenceGuard(const Mesh &mesh);

    /**
     * @brief Empty while the run is healthy, otherwise what went wrong, as
     * the start of a sentence ("The total energy is nan").
     */
    std::string check(const Mesh &mesh) const;

    /// The box a coordinate must stay inside, per axis.
    std::array<double, 3> lower() const { return lower_; }
    std::array<double, 3> upper() const { return upper_; }
    /// The |E_total| beyond which the run counts as diverged.
    double energy_bound() const { return energyBound_; }

private:
    std::array<double, 3> lower_;
    std::array<double, 3> upper_;
    double energyBound_;
};
