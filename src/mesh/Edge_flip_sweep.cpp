/**
 * @file Edge_flip_sweep.cpp
 * @brief The Metropolis edge-flip sweep: in-plane fluidity.
 *
 * A lipid bilayer is a two-dimensional fluid. A triangulated membrane whose
 * connectivity is fixed is not: two vertices that start as neighbours stay
 * neighbours forever, so the sheet carries an in-plane shear modulus a bilayer
 * does not have, and a large shape change has to drag the triangulation with
 * it. Every dynamically triangulated surface model removes that shear modulus
 * the same way -- by flipping the diagonal of a pair of adjacent triangles as
 * a Monte Carlo move -- and this is that move.
 *
 * Two things about the schedule are worth stating, because the literature
 * leaves them implicit and they are where a fluid membrane model most easily
 * goes wrong.
 *
 * **The attempt count is drawn, not fixed.** A flip attempt rate `nu` per edge
 * per unit time defines a Poisson process on every edge, so over a step of
 * length `dt` on a mesh with `N_E` flippable edges the number of attempts is
 * Poisson with mean `nu * dt * N_E`. Fixing the count per step instead would
 * make the physical rate depend on the time step, and halving the step to
 * improve the integrator would silently halve the membrane's fluidity.
 *
 * **The acceptances are sequential.** Evaluating a batch of flips together and
 * accepting them jointly multiplies their individual acceptance probabilities,
 * which is a different Markov chain with a different stationary distribution.
 * TriMem parallelizes the energy evaluation for exactly this reason and keeps
 * the acceptance serial anyway.
 *
 * Why the combination with Brownian dynamics is sound: the Brownian step at
 * fixed connectivity is a kernel whose stationary density is exp(-E/kT), and
 * this sweep at fixed coordinates is a Metropolis kernel with the same
 * stationary density on the joint space of positions and triangulations. A
 * composition of kernels that each leave a measure invariant leaves it
 * invariant. The composite is not reversible -- alternating two reversible
 * kernels never is -- but stationarity is all equilibrium sampling needs.
 *
 * @see docs/edge_flip_plan.md sections 1.2, 1.3 and 3.5
 */

#include "mesh/Mesh.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>

#include "Counter_rng.hpp"

void write_edge_flip_log_csv(const std::vector<EdgeFlipRecord> &records, const std::string &path)
{
    if (records.empty())
    {
        return;
    }

    // The header goes in only when the file is new, so a run that appends each
    // step does not scatter headers through its own log.
    bool needsHeader = true;
    {
        std::ifstream existing(path);
        needsHeader = !existing.good() || existing.peek() == std::ifstream::traits_type::eof();
    }

    std::ofstream out(path, std::ios::app);
    if (!out.good())
    {
        return;
    }
    if (needsHeader)
    {
        out << "iteration,edge,v0,v1,target0,target1,deltaEnergy,accepted\n";
    }
    for (const EdgeFlipRecord &record : records)
    {
        out << record.iteration << ',' << record.edge << ',' << record.vertex[0] << ','
            << record.vertex[1] << ',' << record.vertex[2] << ',' << record.vertex[3] << ','
            << record.deltaEnergy << ',' << (record.accepted ? 1 : 0) << '\n';
    }
}

EdgeFlipSweepStats Mesh::edge_flip_sweep(long long iteration, std::vector<EdgeFlipRecord> *log)
{
    EdgeFlipSweepStats stats;

    if (edges.empty())
    {
        return stats;
    }

    // Only edges a flip could ever touch enter the mean, so that the rate is
    // per *flippable* edge. Counting the frozen boundary band of a periodic
    // sheet would make the physical rate depend on how much of the mesh is
    // ghost, which is a discretization detail and not a property of the
    // membrane.
    std::vector<int> flippable;
    flippable.reserve(edges.size());
    for (const MeshEdge &edge : edges)
    {
        if (edge.flippable)
        {
            flippable.push_back(edge.index);
        }
    }
    if (flippable.empty())
    {
        return stats;
    }

    const double interval = static_cast<double>(std::max(1, param.edgeFlipInterval));
    const double lambda = param.edgeFlipAttemptRate * param.timeStep * interval *
                          static_cast<double>(flippable.size());
    if (!(lambda > 0.0))
    {
        return stats;
    }
    if (lambda > 700.0)
    {
        // slimed::poisson() caps the mean here, so past this point the sweep
        // would silently stop scaling with the rate. Say so once rather than
        // reporting a fluidity the run did not have.
        static bool warned = false;
        if (!warned)
        {
            warned = true;
            std::cout << "[Mesh::edge_flip_sweep] WARNING: " << lambda
                      << " attempts per sweep exceeds the Poisson draw's supported mean of 700. "
                         "Lower edgeFlipAttemptRate or edgeFlipInterval; the rate this run "
                         "reports will not be the rate it applied."
                      << std::endl;
        }
    }

    // One key per (run, sweep). The tag keeps this stream clear of the
    // Brownian displacement's, which keys on (seed, iteration) alone.
    const std::uint64_t sweepKey =
        slimed::splitmix64(static_cast<std::uint64_t>(param.randomSeed) ^ 0x464C4950ULL) +
        static_cast<std::uint64_t>(iteration);

    stats.drawn = slimed::poisson(sweepKey, lambda);

    const double kT = param.KBT;

    for (int attempt = 0; attempt < stats.drawn; attempt++)
    {
        // A fresh sub-key per attempt, so the draw does not depend on how many
        // attempts happened to be admissible before it.
        const std::uint64_t attemptKey =
            slimed::splitmix64(sweepKey + static_cast<std::uint64_t>(attempt) * 0x2545F491ULL);

        // Uniform over the flippable edges. Symmetric proposal: the reverse
        // flip is offered from the new state with the same probability.
        const int slot = slimed::uniform_below(attemptKey, static_cast<int>(flippable.size()));
        const int iEdge = flippable[slot];

        EdgeFlipDelta delta;
        if (!evaluate_edge_flip(iEdge, delta))
        {
            // Refused by the admission test -- a valence bound, or a flip that
            // would duplicate an edge. It is not an attempt in the Metropolis
            // sense, so it is counted separately rather than as a rejection:
            // reporting it as one would understate the acceptance rate by an
            // amount that depends on the mesh rather than on the physics.
            continue;
        }
        stats.attempted++;

        bool accepted = (delta.energy <= 0.0);
        if (!accepted)
        {
            const double probability = (kT > 0.0) ? std::exp(-delta.energy / kT) : 0.0;
            accepted = slimed::uniform_open01(
                           slimed::splitmix64(attemptKey ^ 0xA5A5A5A5A5A5A5A5ULL)) < probability;
        }

        const MeshEdge &edge = edges[iEdge];
        EdgeFlipRecord record;
        record.iteration = iteration;
        record.edge = iEdge;
        record.vertex[0] = edge.v[0];
        record.vertex[1] = edge.v[1];
        record.vertex[2] = edge.opposite[0];
        record.vertex[3] = edge.opposite[1];
        record.deltaEnergy = delta.energy;
        record.accepted = accepted;
        if (log != nullptr)
        {
            log->push_back(record);
        }

        if (!accepted)
        {
            continue;
        }

        flip_edge(iEdge);
        stats.accepted++;
        stats.deltaEnergy += delta.energy;

        // The constraint difference of the next attempt is measured against
        // these, so they have to move with each accepted flip. Recomputing
        // them from the whole mesh would put an O(F) pass inside the sweep and
        // undo the point of a local trial.
        param.area += delta.area;
        param.vol += delta.volume;

        // An accepted flip can make an edge that was flippable stop being so,
        // or the reverse, anywhere in the neighbourhood. flip_edge() has
        // already re-settled those, but the local list this sweep is drawing
        // from was taken before the first flip. Leaving it is deliberate and
        // harmless: a stale entry is caught by evaluate_edge_flip()'s
        // admission test and counted as inadmissible, and a newly flippable
        // edge simply waits for the next sweep. Rebuilding the list per
        // acceptance would bias the choice toward recently disturbed
        // neighbourhoods, which is a worse trade than a slightly stale
        // uniform.
    }

    return stats;
}
