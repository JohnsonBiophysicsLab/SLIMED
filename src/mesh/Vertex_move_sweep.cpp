/**
 * @file Vertex_move_sweep.cpp
 * @brief Metropolis single-vertex moves, and the mesh-quality walls as hard
 * constraints.
 *
 * The Brownian step holds a fluid mesh together with three quadratic walls --
 * on edge length, on triangle altitude, on the crease between faces -- because
 * an explicit step cannot refuse a move and the Hamiltonian has to push back
 * instead. That costs a stiffness the time step must resolve, and it leaves
 * the walls leaned on: a fluid run keeps 100-500 pN.nm in them.
 *
 * A Monte Carlo move can refuse. The dynamically triangulated surface models
 * (Ho & Baumgaertner 1990; Kroll & Gompper 1992; Gompper & Kroll 1997, 2004)
 * therefore use walls of infinite height: beads that may not overlap, tethers
 * that may not stretch past a maximum, and a vertex displacement or a bond
 * flip that would break either is simply rejected. This file is that move for
 * SLIMED, on the control points of the limit surface:
 *
 *   1. pick an independent vertex uniformly at random;
 *   2. displace it by a vector uniform in the cube [-s, s]^3, its periodic
 *      copies with it;
 *   3. under hard walls, reject at once if a wall is broken;
 *   4. otherwise accept with min(1, exp(-dE / kT)), where dE is the exact
 *      change of the energy the run samples.
 *
 * **Why it samples what it claims.** The vertex is chosen with a probability
 * that does not depend on the configuration, and the displacement d is as
 * likely as -d, so the proposal density q(x -> x') equals q(x' -> x). With the
 * Metropolis acceptance that is detailed balance with respect to
 * exp(-E/kT) on the set of admissible meshes -- a hard wall is the statement
 * that E is infinite outside it. No Hastings factor, and none could be added:
 * a proposal *biased* toward better triangles and accepted on E alone samples
 * exp(-E/kT) times whatever the proposal itself prefers, which is the same as
 * having put that preference into the Hamiltonian.
 *
 * **Why the energy difference is local.** A face's energy is a functional of
 * its control net -- its corners and their neighbours -- so a moved vertex
 * changes exactly the faces incident to its closed one-ring. The two global
 * constraints are quadratic in the mesh-wide area and volume and enter
 * through their change, as for the flip.
 *
 * @see docs/monte_carlo_vertex_moves.md
 */

#include "mesh/Mesh.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

#include "Counter_rng.hpp"
#include "energy_force/Patch_kernel.hpp"

namespace
{
/// The three corners of a face, as p[corner][axis].
void face_corners(const Mesh &mesh, int iFace, double p[3][3])
{
    const std::vector<int> &corners = mesh.faces[iFace].adjacentVertices;
    for (int corner = 0; corner < 3; corner++)
    {
        const Matrix &coord = mesh.vertices[corners[corner]].coord;
        for (int axis = 0; axis < 3; axis++)
        {
            p[corner][axis] = coord.get(axis, 0);
        }
    }
}

/// (p1 - p0) x (p2 - p0): the face's unnormalised normal, twice its area long.
void face_normal(const double p[3][3], double n[3])
{
    double d1[3];
    double d2[3];
    for (int axis = 0; axis < 3; axis++)
    {
        d1[axis] = p[1][axis] - p[0][axis];
        d2[axis] = p[2][axis] - p[0][axis];
    }
    slimed::v3_cross(d1, d2, n);
}

/// The wall positions, read once per call from the same parameters the soft
/// terms read.
struct WallBounds
{
    bool tether = false;
    bool altitude = false;
    bool crease = false;
    double lowerLength = 0.0;
    double upperLength = 0.0;
    double altitudeFloor = 0.0;
    double cosineFloor = -1.0;
};

WallBounds wall_bounds(const Param &param)
{
    WallBounds bounds;
    const double restLength =
        (param.edgeSpringRestLength >= 0.0) ? param.edgeSpringRestLength : param.lFace;
    bounds.tether = param.edgeSpringEnabled;
    bounds.lowerLength = param.edgeTetherMinRatio * restLength;
    bounds.upperLength = param.edgeTetherMaxRatio * restLength;
    bounds.altitude = param.triangleShapeEnabled;
    bounds.altitudeFloor = param.triangleShapeMinAltitudeRatio * param.lFace;
    bounds.crease = param.creaseWallEnabled;
    bounds.cosineFloor = std::cos(param.creaseWallAngle * M_PI / 180.0);
    return bounds;
}

/// Whether an edge's length is inside the tether range. The comparisons are
/// the ones Mesh::edge_tether_energy() makes, so the hard wall sits exactly
/// where the soft one starts charging.
bool edge_length_allowed(const Mesh &mesh, const MeshEdge &edge, const WallBounds &bounds)
{
    double along[3];
    for (int axis = 0; axis < 3; axis++)
    {
        along[axis] = mesh.vertices[edge.v[0]].coord.get(axis, 0) -
                      mesh.vertices[edge.v[1]].coord.get(axis, 0);
    }
    const double length = slimed::v3_norm(along);
    return !(length < bounds.lowerLength) && !(length > bounds.upperLength);
}

/// Whether the crease across an interior edge is within the angle. A face
/// without a normal is the altitude wall's case, as in Crease_wall.cpp.
bool edge_crease_allowed(const Mesh &mesh, const MeshEdge &edge, const WallBounds &bounds)
{
    double p1[3][3];
    double p2[3][3];
    face_corners(mesh, edge.face[0], p1);
    face_corners(mesh, edge.face[1], p2);
    double n1[3];
    double n2[3];
    face_normal(p1, n1);
    face_normal(p2, n2);
    const double length1 = slimed::v3_norm(n1);
    const double length2 = slimed::v3_norm(n2);
    if (length1 <= 0.0 || length2 <= 0.0)
    {
        return true;
    }
    const double cosine = slimed::v3_dot(n1, n2) / (length1 * length2);
    return !(cosine < bounds.cosineFloor);
}

/// Whether all three altitudes of a face clear the floor: h_i = 2A / l_i, as
/// Triangle_shape.cpp defines them.
bool face_altitudes_allowed(const Mesh &mesh, int iFace, const WallBounds &bounds)
{
    double p[3][3];
    face_corners(mesh, iFace, p);
    double n[3];
    face_normal(p, n);
    const double twoArea = slimed::v3_norm(n);
    for (int i = 0; i < 3; i++)
    {
        double e[3];
        for (int axis = 0; axis < 3; axis++)
        {
            e[axis] = p[(i + 1) % 3][axis] - p[(i + 2) % 3][axis];
        }
        const double length = slimed::v3_norm(e);
        if (length <= 0.0)
        {
            continue; // coincident corners: no edge to measure against
        }
        if (twoArea / length < bounds.altitudeFloor)
        {
            return false;
        }
    }
    return true;
}
} // namespace

void write_vertex_move_log_csv(const std::string &path, long long iteration,
                               const VertexMoveSweepStats &sinceLastRow,
                               const WallViolations &walls)
{
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
        out << "iteration,attempted,wallRejected,accepted,acceptance,deltaEnergy,"
               "tetherOutside,altitudeOutside,creaseOutside\n";
    }
    out << iteration << ',' << sinceLastRow.attempted << ',' << sinceLastRow.wallRejected << ','
        << sinceLastRow.accepted << ',' << sinceLastRow.acceptance() << ','
        << sinceLastRow.deltaEnergy << ',' << walls.tether << ',' << walls.altitude << ','
        << walls.crease << '\n';
}

bool Mesh::faces_satisfy_walls(const std::vector<int> &faceList, WallViolations *tally) const
{
    const WallBounds bounds = wall_bounds(param);
    bool satisfied = true;

    for (int iFace : faceList)
    {
        const Face &face = faces[iFace];
        if (face.adjacentVertices.size() != 3)
        {
            continue;
        }

        if (bounds.tether || bounds.crease)
        {
            for (int k = 0; k < 3; k++)
            {
                const int iEdge =
                    edge_between(face.adjacentVertices[k], face.adjacentVertices[(k + 1) % 3]);
                if (iEdge < 0)
                {
                    continue;
                }
                const MeshEdge &edge = edges[iEdge];
                if (!edge_carries_tether(edge))
                {
                    continue; // a stale copy of an interior edge; see the header
                }
                if (bounds.tether && !edge_length_allowed(*this, edge, bounds))
                {
                    satisfied = false;
                    if (tally == nullptr)
                    {
                        return false;
                    }
                    tally->tether++;
                }
                if (bounds.crease && !edge.is_boundary() &&
                    !edge_crease_allowed(*this, edge, bounds))
                {
                    satisfied = false;
                    if (tally == nullptr)
                    {
                        return false;
                    }
                    tally->crease++;
                }
            }
        }

        if (bounds.altitude && !face.isGhost && !face_altitudes_allowed(*this, iFace, bounds))
        {
            satisfied = false;
            if (tally == nullptr)
            {
                return false;
            }
            tally->altitude++;
        }
    }
    return satisfied;
}

WallViolations Mesh::count_wall_violations() const
{
    const WallBounds bounds = wall_bounds(param);
    WallViolations tally;

    for (const MeshEdge &edge : edges)
    {
        if (!edge_carries_tether(edge))
        {
            continue;
        }
        if (bounds.tether && !edge_length_allowed(*this, edge, bounds))
        {
            tally.tether++;
        }
        if (bounds.crease && !edge.is_boundary() && !edge_crease_allowed(*this, edge, bounds))
        {
            tally.crease++;
        }
    }
    if (bounds.altitude)
    {
        for (int iFace = 0; iFace < static_cast<int>(faces.size()); iFace++)
        {
            const Face &face = faces[iFace];
            if (face.isGhost || face.adjacentVertices.size() != 3)
            {
                continue;
            }
            if (!face_altitudes_allowed(*this, iFace, bounds))
            {
                tally.altitude++;
            }
        }
    }
    return tally;
}

void Mesh::build_vertex_move_table()
{
    const int nVertices = static_cast<int>(vertices.size());
    vertexMoveCandidates.clear();
    vertexMoveFollowers.assign(nVertices, std::vector<int>());

    switch (param.boundaryCondition)
    {
    // Per-vertex boundary: the mesh already knows which vertices are
    // coordinates and which are images of which.
    case BoundaryType::Mixed:
    {
        for (int v = 0; v < nVertices; v++)
        {
            if (is_independent_vertex(v))
            {
                vertexMoveCandidates.push_back(v);
            }
        }
        for (int image : periodicImageVertices)
        {
            const int source = vertices[image].reflectiveVertexIndex;
            if (source >= 0 && source < nVertices)
            {
                vertexMoveFollowers[source].push_back(image);
            }
        }
        break;
    }

    // The generated periodic sheet. The periodic tile is the block of
    // (nFaceX - 6) x (nFaceY - 6) vertices starting at index (3, 3); every
    // other vertex -- the three ghost rings, and the ring that duplicates the
    // tile's first row and column -- is the tile vertex a whole number of
    // periods away. The Brownian step treats both copies of the duplicated
    // ring as slaved to each other, which leaves that ring to be moved only
    // by the surface solve; a Monte Carlo move needs every coordinate to have
    // exactly one owner, so the tile's copy is the owner here.
    case BoundaryType::Periodic:
    {
        const int nVertX = param.nFaceX + 1;
        const int nVertY = param.nFaceY + 1;
        const int periodX = param.nFaceX - 6;
        const int periodY = param.nFaceY - 6;
        if (nVertX * nVertY != nVertices || periodX <= 0 || periodY <= 0)
        {
            throw std::runtime_error(
                "[Mesh::build_vertex_move_table] boundaryType = Periodic expects the generated "
                "flat sheet, with at least seven faces along each axis. For any other mesh use "
                "boundaryType = Mixed.");
        }
        const auto wrap = [](int index, int period) {
            const int shifted = (index - 3) % period;
            return 3 + (shifted < 0 ? shifted + period : shifted);
        };
        for (int row = 0; row < nVertY; row++)
        {
            for (int col = 0; col < nVertX; col++)
            {
                const int v = row * nVertX + col;
                const int source = wrap(row, periodY) * nVertX + wrap(col, periodX);
                if (source == v)
                {
                    vertexMoveCandidates.push_back(v);
                }
                else
                {
                    vertexMoveFollowers[source].push_back(v);
                }
            }
        }
        break;
    }

    // Clamped rim: the vertices of the boundary faces are the ones the force
    // pass zeroes, so they are the ones that do not move. A closed surface
    // has no boundary face and every vertex is a candidate.
    case BoundaryType::Fixed:
    {
        std::vector<char> clamped(nVertices, 0);
        for (const Face &face : faces)
        {
            if (!face.isBoundary)
            {
                continue;
            }
            for (int corner : face.adjacentVertices)
            {
                clamped[corner] = 1;
            }
        }
        for (int v = 0; v < nVertices; v++)
        {
            if (clamped[v] == 0 && !vertices[v].isGhost)
            {
                vertexMoveCandidates.push_back(v);
            }
        }
        break;
    }

    case BoundaryType::Free:
        throw std::runtime_error(
            "[Mesh::build_vertex_move_table] vertex moves are not defined under the global "
            "boundaryType = Free, whose rim is an extrapolation of the interior rather than a "
            "set of coordinates. Use boundaryType = Mixed with boundaryTypeX/Y = Free.");
    }
}

std::vector<int> Mesh::vertex_move_patch_faces(int iVertex) const
{
    std::vector<int> patch;
    const auto addClosedOneRing = [&](int moved) {
        const std::vector<int> &own = vertices[moved].adjacentFaces;
        patch.insert(patch.end(), own.begin(), own.end());
        for (int neighbour : vertices[moved].adjacentVertices)
        {
            const std::vector<int> &incident = vertices[neighbour].adjacentFaces;
            patch.insert(patch.end(), incident.begin(), incident.end());
        }
    };
    addClosedOneRing(iVertex);
    if (iVertex < static_cast<int>(vertexMoveFollowers.size()))
    {
        for (int follower : vertexMoveFollowers[iVertex])
        {
            addClosedOneRing(follower);
        }
    }
    std::sort(patch.begin(), patch.end());
    patch.erase(std::unique(patch.begin(), patch.end()), patch.end());
    return patch;
}

void Mesh::displace_vertex_with_followers(int iVertex, const double displacement[3])
{
    const auto displace = [&](int v) {
        Matrix &coord = vertices[v].coord;
        for (int axis = 0; axis < 3; axis++)
        {
            coord.set(axis, 0, coord.get(axis, 0) + displacement[axis]);
        }
    };
    displace(iVertex);
    if (iVertex < static_cast<int>(vertexMoveFollowers.size()))
    {
        for (int follower : vertexMoveFollowers[iVertex])
        {
            displace(follower);
        }
    }
}

VertexMoveDelta Mesh::evaluate_vertex_move(int iVertex, const double displacement[3])
{
    VertexMoveDelta delta;

    if (vertexMoveFollowers.size() != vertices.size())
    {
        build_vertex_move_table();
    }

    // Everything that moves, and where it was: the trial is undone by putting
    // these back rather than by subtracting the displacement, which would not
    // round-trip in floating point and would let rejected moves random-walk
    // the mesh by an ulp at a time.
    std::vector<int> moved;
    moved.push_back(iVertex);
    moved.insert(moved.end(), vertexMoveFollowers[iVertex].begin(),
                 vertexMoveFollowers[iVertex].end());
    std::vector<double> original(moved.size() * 3);
    for (std::size_t i = 0; i < moved.size(); i++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            original[i * 3 + axis] = vertices[moved[i]].coord.get(axis, 0);
        }
    }
    const auto restore = [&]() {
        for (std::size_t i = 0; i < moved.size(); i++)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                vertices[moved[i]].coord.set(axis, 0, original[i * 3 + axis]);
            }
        }
    };

    displace_vertex_with_followers(iVertex, displacement);

    // The walls first: they are a few cross products against two dozen patch
    // integrals, and under hard walls a good share of the rejections are
    // theirs. Only a face with a moved corner can have changed an edge
    // length, an altitude or a crease.
    if (walls_are_hard())
    {
        std::vector<int> incident;
        for (int v : moved)
        {
            incident.insert(incident.end(), vertices[v].adjacentFaces.begin(),
                            vertices[v].adjacentFaces.end());
        }
        std::sort(incident.begin(), incident.end());
        incident.erase(std::unique(incident.begin(), incident.end()), incident.end());
        if (!faces_satisfy_walls(incident))
        {
            restore();
            delta.wallViolated = true;
            delta.energy = std::numeric_limits<double>::infinity();
            return delta;
        }
    }

    const std::vector<int> patch = vertex_move_patch_faces(iVertex);
    const FaceSubsetEnergy after = evaluate_face_subset(patch);
    restore();
    const FaceSubsetEnergy before = evaluate_face_subset(patch);

    delta.bending = after.bending - before.bending;
    // The mesh-quality energy joins only where the dynamics integrates it:
    // the soft fluid walls. Solid mode's reference-length term is evaluated
    // by the force pass but is not in the drive, and a Metropolis move that
    // sampled it would not share a Hamiltonian with the Brownian step.
    if (param.edgeSpringEnabled)
    {
        delta.regularization = after.regularization - before.regularization;
    }
    delta.area = after.area - before.area;
    delta.volume = after.volume - before.volume;

    // The global constraints, exactly as for a flip: quadratic in the
    // mesh-wide total, so the change depends on where the membrane sits.
    if (lipid_flow_local() && !faceArea0.empty())
    {
        // Local lipid flow: a vertex move changes the faces' areas at fixed
        // relaxed areas, and the sum of their own stretching energies is
        // already in the subset.
        delta.areaConstraint = after.stretching - before.stretching;
    }
    else if (lipid_flow_global() && lipidAreaNormalization > 0.0)
    {
        delta.areaConstraint = param.uSurf / (2.0 * lipidAreaNormalization) * delta.area *
                               (delta.area + 2.0 * (param.area - param.area0));
    }
    else if (param.isGlobalConstraint && param.area0 != 0.0)
    {
        delta.areaConstraint = param.uSurf / (2.0 * param.area0) * delta.area *
                               (delta.area + 2.0 * (param.area - param.area0));
    }
    if (param.vol0 != 0.0)
    {
        delta.volumeConstraint = param.uVol / (2.0 * param.vol0) * delta.volume *
                                 (delta.volume + 2.0 * (param.vol - param.vol0));
    }

    delta.energy = delta.bending + delta.regularization + delta.areaConstraint +
                   delta.volumeConstraint;
    return delta;
}

VertexMoveSweepStats Mesh::vertex_move_sweep(long long iteration)
{
    VertexMoveSweepStats stats;

    if (vertexMoveFollowers.size() != vertices.size())
    {
        build_vertex_move_table();
    }
    const int nCandidates = static_cast<int>(vertexMoveCandidates.size());
    if (nCandidates == 0)
    {
        return stats;
    }

    // One key per (run, sweep), tagged so the stream is clear of the Brownian
    // displacement's and the flip sweep's.
    const std::uint64_t sweepKey =
        slimed::splitmix64(static_cast<std::uint64_t>(param.randomSeed) ^ 0x564D4F5645ULL) +
        static_cast<std::uint64_t>(iteration);

    // A whole number of attempts with the requested mean. Rounding the
    // fraction at random rather than always down keeps the mean exact for a
    // rate like 0.5 per vertex on a small mesh.
    const double mean = std::max(0.0, param.vertexMoveAttemptsPerVertex) * nCandidates;
    long long nAttempts = static_cast<long long>(std::floor(mean));
    if (slimed::uniform_open01(slimed::splitmix64(sweepKey ^ 0x5DEECE66DULL)) <
        mean - std::floor(mean))
    {
        nAttempts++;
    }

    const double step = param.vertexMoveStepRatio * param.lFace;
    const double kT = param.KBT;
    const bool inPlane = param.inPlaneDynamicsEnabled;

    for (long long attempt = 0; attempt < nAttempts; attempt++)
    {
        const std::uint64_t attemptKey =
            slimed::splitmix64(sweepKey + static_cast<std::uint64_t>(attempt) * 0x2545F491ULL);

        // Uniform over the independent vertices, whatever the configuration:
        // the reverse move is offered with the same probability.
        const int iVertex = vertexMoveCandidates[slimed::uniform_below(attemptKey, nCandidates)];

        // Uniform in the cube, so d and -d are equally likely. With in-plane
        // motion off the vertex keeps its lattice site, as it does under the
        // Brownian step, and only its height is sampled.
        double displacement[3] = {0.0, 0.0, 0.0};
        for (int axis = (inPlane ? 0 : 2); axis < 3; axis++)
        {
            const double u = slimed::uniform_open01(slimed::splitmix64(
                attemptKey ^ (0xC2B2AE3D27D4EB4FULL * static_cast<std::uint64_t>(axis + 1))));
            displacement[axis] = step * (2.0 * u - 1.0);
        }

        const VertexMoveDelta delta = evaluate_vertex_move(iVertex, displacement);
        stats.attempted++;
        if (delta.wallViolated)
        {
            stats.wallRejected++;
            continue;
        }

        bool accepted = (delta.energy <= 0.0);
        if (!accepted)
        {
            // A non-finite difference compares false both here and above, so
            // a move into a configuration the evaluator cannot describe is
            // rejected rather than accepted by accident.
            const double probability = (kT > 0.0) ? std::exp(-delta.energy / kT) : 0.0;
            accepted = slimed::uniform_open01(
                           slimed::splitmix64(attemptKey ^ 0xA5A5A5A5A5A5A5A5ULL)) < probability;
        }
        if (!accepted)
        {
            continue;
        }

        displace_vertex_with_followers(iVertex, displacement);
        stats.accepted++;
        stats.deltaEnergy += delta.energy;

        // The next attempt's constraint difference is measured against these.
        param.area += delta.area;
        param.vol += delta.volume;
    }

    return stats;
}
