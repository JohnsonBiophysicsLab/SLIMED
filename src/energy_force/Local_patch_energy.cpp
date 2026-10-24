/**
 * @file Local_patch_energy.cpp
 * @brief Energy over a subset of faces, and the energy change of one edge flip.
 *
 * A Metropolis trial has to answer "what would this move cost?" thousands of
 * times per run, and answering it with a whole-mesh energy evaluation would
 * make the sweep cost more than the dynamics it is interleaved with. It does
 * not have to: the energy of a face is a functional of its control net -- the
 * union of its three corners' one-rings -- so a flip can only change the
 * energy of a face incident to one of the four vertices it touches. Everything
 * else cancels exactly, not approximately.
 *
 * The two global constraints are the exception, and they are the reason this
 * file exists rather than a simple sum. Area and volume enter the Hamiltonian
 * quadratically in their mesh-wide totals, so there is no per-face term to
 * difference. Their change is still exact, but it depends on how far the
 * membrane already sits from its reference as well as on how much this flip
 * moves it -- the form OrganL and FreeDTS both use for a local move under a
 * global restraint.
 *
 * @see docs/edge_flip_plan.md work packages 2 and 3
 */

#include "mesh/Mesh.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "energy_force/Patch_kernel.hpp"
#include "mesh/Multi_extraordinary_patch.hpp"

namespace
{
/// out = a - b, for the (3, 1) coordinate columns a Vertex holds.
void coord_difference(const Matrix &a, const Matrix &b, double out[3])
{
    for (int axis = 0; axis < 3; axis++)
    {
        out[axis] = a.get(axis, 0) - b.get(axis, 0);
    }
}
} // namespace

double Mesh::face_regularization_energy(int iFace) const
{
    const Face &face = faces[iFace];
    const double kCurv = param.kCurv;

    // The same rule energy_force_regularization() applies: under the
    // per-vertex boundary a copy face duplicates a physical one and is not
    // charged, so the flip trial and the force pass agree about the total.
    if (param.boundaryCondition == BoundaryType::Mixed && face.isGhost)
    {
        return 0.0;
    }

    const int iVertex0 = face.adjacentVertices[0];
    const int iVertex1 = face.adjacentVertices[1];
    const int iVertex2 = face.adjacentVertices[2];

    double vector10[3];
    double vector21[3];
    double vector02[3];
    coord_difference(vertices[iVertex0].coord, vertices[iVertex1].coord, vector10);
    coord_difference(vertices[iVertex1].coord, vertices[iVertex2].coord, vector21);
    coord_difference(vertices[iVertex2].coord, vertices[iVertex0].coord, vector02);
    const double length10 = slimed::v3_norm(vector10);
    const double length21 = slimed::v3_norm(vector21);
    const double length02 = slimed::v3_norm(vector02);

    double semiperi = (length10 + length21 + length02) / 2.0;
    const double area =
        sqrt(semiperi * (semiperi - length10) * (semiperi - length21) * (semiperi - length02));

    const double meanSideLength = (length10 + length21 + length02) / 3.0;
    const double gama = (pow(length10 - meanSideLength, 2.0) +
                         pow(length21 - meanSideLength, 2.0) +
                         pow(length02 - meanSideLength, 2.0)) /
                        pow(meanSideLength, 2.0);

    double refVector10[3];
    double refVector21[3];
    double refVector02[3];
    coord_difference(vertices[iVertex0].coordRef, vertices[iVertex1].coordRef, refVector10);
    coord_difference(vertices[iVertex1].coordRef, vertices[iVertex2].coordRef, refVector21);
    coord_difference(vertices[iVertex2].coordRef, vertices[iVertex0].coordRef, refVector02);
    const double refLength10 = slimed::v3_norm(refVector10);
    const double refLength21 = slimed::v3_norm(refVector21);
    const double refLength02 = slimed::v3_norm(refVector02);

    semiperi = (refLength10 + refLength21 + refLength02) / 2.0;
    const double refArea = sqrt(semiperi * (semiperi - refLength10) * (semiperi - refLength21) *
                                (semiperi - refLength02));

    const bool isDeformShape = (gama > param.gamaShape && param.usingRpi);
    const bool isDeformArea = (std::abs(area - refArea) / refArea >= param.gamaArea &&
                               param.usingRpi);

    // The same three-way switch energy_force_regularization() takes, and in the
    // same order, so the two cannot disagree about which case a face is in.
    const int deformationCase = (isDeformShape ? 2 : 0) | (isDeformArea ? 1 : 0);
    switch (deformationCase)
    {
    case 0:
        return kCurv / 2.0 *
               (pow(length10 - refLength10, 2.0) + pow(length21 - refLength21, 2.0) +
                pow(length02 - refLength02, 2.0));
    case 1:
    {
        const double meanSideLengthRef = sqrt(4.0 * refArea / sqrt(3.0));
        return kCurv / 2.0 *
               (pow(length10 - meanSideLengthRef, 2.0) + pow(length21 - meanSideLengthRef, 2.0) +
                pow(length02 - meanSideLengthRef, 2.0));
    }
    default:
    {
        const double meanSideLengthOld = sqrt(4.0 * area / sqrt(3.0));
        return kCurv / 2.0 *
               (pow(length10 - meanSideLengthOld, 2.0) + pow(length21 - meanSideLengthOld, 2.0) +
                pow(length02 - meanSideLengthOld, 2.0));
    }
    }
}

FaceSubsetEnergy Mesh::evaluate_face_subset(const std::vector<int> &faceList,
                                            std::vector<double> *faceAreas)
{
    FaceSubsetEnergy result;
    if (faceAreas != nullptr)
    {
        // One entry per face of the list, in its order; a face that carries
        // no area (ghost, boundary, no patch) reads zero.
        faceAreas->assign(faceList.size(), 0.0);
    }

    // Built lazily, and not thread safe, so it happens before anything else.
    ensure_patch_rows_flat();
    const double *const regularRows = patchRowsFlat.regular();
    const double *const gaussCoeff = patchRowsFlat.gaussCoeff();
    const int nSamples = patchRowsFlat.nSamples();

    // The same loop-invariant constants the whole-mesh pass hoists, including
    // the two guards: a surface enclosing nothing has vol0 == 0, and dividing
    // by it would poison the result with NaN.
    slimed::PatchParams patchParams;
    patchParams.kCurv = param.kCurv;
    patchParams.uSurfPerArea = (param.area0 == 0.0) ? 0.0 : param.uSurf / param.area0;
    patchParams.area = param.area;
    patchParams.area0 = param.area0;
    patchParams.uVol = (param.vol0 == 0.0) ? 0.0 : param.uVol / param.vol0;
    patchParams.vol = param.vol;
    patchParams.vol0 = param.vol0;

    for (std::size_t k = 0; k < faceList.size(); k++)
    {
        const int iFace = faceList[k];
        const Face &face = faces[iFace];

        // Whichever mesh-quality term the dynamics is actually integrating.
        // These are alternatives, and a trial that differenced the wrong one
        // would hand the Metropolis sweep a Hamiltonian the Brownian step does
        // not share: the chain would then sample neither. Measured when this
        // read face_regularization_energy() unconditionally, on a run with the
        // tether on -- every accepted flip reported about -750 pN.nm of a
        // reference length that a newly created edge never had, and the mesh's
        // total energy climbed while the flips claimed to be lowering it.
        //
        // Both run over every face, ghosts included, and every face's energy is
        // summed into the total. A flip never touches a ghost, so that only
        // ever matters for consistency.
        //
        // Under hard walls the fluid term carries no energy at all -- the
        // walls are constraints, checked by faces_satisfy_walls() -- which is
        // what energy_force_fluid_terms() does on the whole mesh.
        const bool hardWalls = walls_are_hard();
        if (!param.edgeSpringEnabled)
        {
            result.regularization += face_regularization_energy(iFace);
        }
        else if (!hardWalls)
        {
            result.regularization += face_tether_energy(iFace);
        }
        if (param.triangleShapeEnabled && !hardWalls)
        {
            result.regularization += face_shape_energy(iFace);
        }
        if (param.creaseWallEnabled && !hardWalls)
        {
            result.regularization += face_crease_energy(iFace);
        }

        if (face.isBoundary)
        {
            continue; // no bending energy, matching the force loop
        }

        const int nOneRingVertices = static_cast<int>(face.oneRingVertices.size());
        const bool isRegular = (face.patchKind == PatchKind::Regular);
        const bool isIrregular = (face.patchKind == PatchKind::SingleExtraordinary);
        const bool isMulti =
            (face.patchKind == PatchKind::MultiExtraordinary && face.patchEntry >= 0);
        if (!(isRegular || isIrregular || isMulti) || nOneRingVertices <= 0 ||
            nOneRingVertices > slimed::kMaxControlPoints)
        {
            continue;
        }

        double coordOneRingVertices[slimed::kMaxControlPoints * 3];
        for (int j = 0; j < nOneRingVertices; j++)
        {
            const Matrix &coord = vertices[face.oneRingVertices[j]].coord;
            coordOneRingVertices[j * 3 + 0] = coord.get(0, 0);
            coordOneRingVertices[j * 3 + 1] = coord.get(1, 0);
            coordOneRingVertices[j * 3 + 2] = coord.get(2, 0);
        }

        slimed::PatchParams facePatchParams = patchParams;
        facePatchParams.spontCurv = face.spontCurvature;

        // Forces are computed and discarded. Calling the same kernel rather
        // than a cut-down one is deliberate: an energy-only kernel would be a
        // second implementation of the integrand, and the whole value of this
        // routine is that it agrees with the whole-mesh pass exactly.
        double fBend[slimed::kMaxControlPoints * 3];
        double fArea[slimed::kMaxControlPoints * 3];
        double fVolume[slimed::kMaxControlPoints * 3];

        auto evaluateOnePatch = [&](int valence, const double *coords, int nCtrl, double &bending,
                                    double &area, double &volume) {
            std::fill(fBend, fBend + nCtrl * 3, 0.0);
            std::fill(fArea, fArea + nCtrl * 3, 0.0);
            std::fill(fVolume, fVolume + nCtrl * 3, 0.0);
            double meanCurv = 0.0;
            double normVector[3] = {0.0, 0.0, 0.0};
            bending = 0.0;

            if (valence == 6)
            {
                slimed::element_energy_force_patch_pod(regularRows, gaussCoeff, nSamples, coords,
                                                       nCtrl, facePatchParams, bending, meanCurv,
                                                       normVector, fBend, fArea, fVolume);
                slimed::element_area_volume_pod(regularRows, gaussCoeff, nSamples, coords, nCtrl,
                                                area, volume);
                return;
            }
            for (int d = 0; d < irregularRows.depth_for(valence); d++)
            {
                // The same tail weights the whole-mesh passes use.
                double bendWeight = 1.0;
                double areaWeight = 1.0;
                patchRowsFlat.level_weights(valence, d, bendWeight, areaWeight);
                for (int c = 0; c < kRegularChildrenPerStep; c++)
                {
                    const double *const rows = patchRowsFlat.child(valence, d, c);
                    double childBending = 0.0;
                    double childMeanCurv = 0.0;
                    double childNormVector[3] = {0.0, 0.0, 0.0};
                    slimed::element_energy_force_patch_weighted_pod(
                        rows, gaussCoeff, nSamples, coords, nCtrl, facePatchParams, bendWeight,
                        areaWeight, childBending, childMeanCurv, childNormVector, fBend, fArea,
                        fVolume);
                    bending += childBending;
                    slimed::element_area_volume_weighted_pod(rows, gaussCoeff, nSamples, coords,
                                                             nCtrl, areaWeight, area, volume);
                }
            }
        };

        double bending = 0.0;
        double area = 0.0;
        double volume = 0.0;

        if (isRegular || isIrregular)
        {
            evaluateOnePatch(isRegular ? 6 : nOneRingVertices - 6, coordOneRingVertices,
                             nOneRingVertices, bending, area, volume);
        }
        else
        {
            const MultiPatchTable::Entry &entry = multiPatchTable.entry(face.patchEntry);
            const double *const prolongations = multiPatchTable.data();
            for (int c = 0; c < 4; c++)
            {
                const MultiPatchTable::Child &child = entry.children[c];
                double childCoords[slimed::kMaxControlPoints * 3];
                multi_patch_prolong(prolongations + child.offset, coordOneRingVertices,
                                    entry.nControl, child.nControl, childCoords);
                double childBending = 0.0;
                evaluateOnePatch(child.valence, childCoords, child.nControl, childBending, area,
                                 volume);
                bending += childBending;
            }
        }

        result.bending += bending;
        // Ghost faces contribute no area or volume, exactly as
        // calculate_element_area_volume() and sum_membrane_area_and_volume()
        // agree they do not.
        if (!face.isGhost)
        {
            result.area += area;
            result.volume += volume;
            if (faceAreas != nullptr)
            {
                (*faceAreas)[k] = area;
            }
            // Local lipid flow: the face's own stretching energy, at the
            // area just integrated and the relaxed area it carries now.
            if (lipid_flow_local() && !faceArea0.empty() && faceArea0[iFace] > 0.0)
            {
                const double modulusArea = lipid_modulus_area(iFace);
                if (modulusArea > 0.0)
                {
                    const double excess = area - faceArea0[iFace];
                    result.stretching += 0.5 * param.uSurf / modulusArea * excess * excess;
                }
            }
        }
    }

    return result;
}

bool Mesh::evaluate_edge_flip(int iEdge, EdgeFlipDelta &delta, std::string *why)
{
    delta = EdgeFlipDelta{};

    if (!edge_flip_is_admissible(iEdge, why))
    {
        return false;
    }

    // Invariant under the move: face0 and face1 each keep two of the four
    // quadrilateral corners, so the set of faces incident to those four is the
    // same before and after. Taking it once means the two measurements cover
    // exactly the same faces, which is what makes everything outside cancel.
    const std::vector<int> patch = flip_patch_faces(iEdge);

    // Which faces carry a control net now. A flip must not take that away
    // from any of them -- see the check below.
    std::vector<char> hadPatch(patch.size(), 0);
    for (std::size_t i = 0; i < patch.size(); i++)
    {
        hadPatch[i] = faces[patch[i]].oneRingVertices.empty() ? 0 : 1;
    }

    // Local lipid flow: the two retriangulated faces carry the relaxed areas
    // the accepted flip would hand them, so the trial differences the energy
    // the mesh would actually carry. Put back with the flip. Under flip
    // resampling the per-face areas are needed as well (see below).
    const bool localLipids = lipid_flow_local() && !faceArea0.empty();
    const bool resampleLipids = localLipids && lipid_flip_resampling();
    std::vector<double> areasBefore;
    std::vector<double> areasAfter;
    const FaceSubsetEnergy before =
        evaluate_face_subset(patch, resampleLipids ? &areasBefore : nullptr);

    // The trial below flips and flips back, and each of those bumps
    // topologyVersion. That leaves the mesh exactly as it was but the version
    // two ahead, which is not a harmless overcount: the version is the signal
    // every connectivity-keyed cache invalidates on, and it is what says
    // whether a trajectory frame needs its connectivity written beside it. A
    // rejected attempt would rebuild the device layout and the sparse limit
    // mask, and write a duplicate face frame, for a mesh that never moved --
    // and on a fluid run most attempts are rejected. Put it back with the
    // mesh.
    const long long versionBeforeTrial = topologyVersion;
    const MeshEdge edgeBeforeTrial = edges[iEdge];
    flip_edge(iEdge);

    // A face whose one-ring cannot be built carries no energy at all, and that
    // is not a neutral outcome for a Monte Carlo move: zero is the lowest
    // energy there is, so a chain allowed to reach such a configuration would
    // be actively drawn into it. The Hamiltonian would develop a hole and the
    // membrane would tear along it.
    //
    // Refusing the move is the honest fix. It keeps the chain inside the set
    // of configurations the model can actually describe, which is where a
    // Metropolis chain has to stay for its stationary distribution to mean
    // anything. Rare in practice -- the configurations that trigger it are
    // strained enough that the energy would usually reject them anyway -- but
    // "usually" is not a guarantee when the alternative is a free lunch.
    //
    // Note that this guards evaluate_edge_flip(), not flip_edge(). The
    // primitive stays unguarded on purpose: it is what the trial itself uses
    // to look ahead, and what a test uses to drive the mesh somewhere
    // deliberately.
    bool lostAPatch = false;
    for (std::size_t i = 0; i < patch.size(); i++)
    {
        if (hadPatch[i] != 0 && faces[patch[i]].oneRingVertices.empty())
        {
            lostAPatch = true;
            break;
        }
    }
    if (lostAPatch)
    {
        flip_edge(iEdge); // restore
        topologyVersion = versionBeforeTrial;
        if (why != nullptr)
        {
            *why = "edge " + std::to_string(iEdge) +
                   " would leave a face in its neighbourhood without a subdivision patch, so "
                   "that face would carry no energy and no force";
        }
        return false;
    }

    // Hard walls: a flip that leaves one of its two new triangles outside a
    // wall is an ordinary Metropolis rejection -- the forbidden state has
    // zero weight -- so it is reported as admissible with an infinite energy
    // and counted as an attempt. Only the two retriangulated faces can have
    // changed: their five edges' lengths and creases, and their altitudes.
    if (walls_are_hard())
    {
        const MeshEdge &flipped = edges[iEdge];
        const std::vector<int> retriangulated = {flipped.face[0], flipped.face[1]};
        if (!faces_satisfy_walls(retriangulated))
        {
            flip_edge(iEdge); // restore
            topologyVersion = versionBeforeTrial;
            delta.wallViolated = true;
            delta.energy = std::numeric_limits<double>::infinity();
            return true;
        }
    }

    if (localLipids)
    {
        lipid_carry_after_flip(iEdge, edgeBeforeTrial);
    }
    const FaceSubsetEnergy after =
        evaluate_face_subset(patch, resampleLipids ? &areasAfter : nullptr);
    const MeshEdge edgeAfterTrial = edges[iEdge];
    flip_edge(iEdge); // restore
    if (localLipids)
    {
        // The carry is an involution: carrying back through the reverse flip
        // puts every value on the triangle it came from, bitwise.
        lipid_carry_after_flip(iEdge, edgeAfterTrial);
    }
    topologyVersion = versionBeforeTrial;

    delta.bending = after.bending - before.bending;
    delta.regularization = after.regularization - before.regularization;
    delta.area = after.area - before.area;
    delta.volume = after.volume - before.volume;

    // The global constraints. Quadratic in the mesh-wide total, so the change
    // depends on where the membrane already sits as well as on how far this
    // flip moves it:
    //
    //     E   = (u / 2 X0) (X - X0)^2
    //     dE  = (u / 2 X0) dX (dX + 2 (X - X0))
    //
    // param.area and param.vol are the totals as they stand now, which is what
    // the last Compute_Energy_And_Force() left there. A sweep that accepts
    // several flips has to keep them current between attempts.
    if (resampleLipids)
    {
        // The compound move "flip, then redraw the set's relaxed areas at the
        // new geometry" is accepted on the set's stretching energy minimised
        // over its split -- the Gaussian normalisations cancel -- with the
        // rest of the patch at the relaxed areas it carries.
        const std::vector<int> pooled = lipid_resample_set(iEdge);
        delta.areaConstraint = lipid_stretching_energy_with_pooled(patch, areasAfter, pooled) -
                               lipid_stretching_energy_with_pooled(patch, areasBefore, pooled);
    }
    else if (localLipids)
    {
        // Local lipid flow: the area term is a sum over faces, already in
        // the subset energies.
        delta.areaConstraint = after.stretching - before.stretching;
    }
    else if (lipid_flow_global() && lipidAreaNormalization > 0.0)
    {
        // Global lipid flow: the same quadratic, normalised by the relaxed
        // area as built rather than the one that moves.
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

    // The scaffolding term does not appear: it is a function of vertex
    // positions and a flip moves no vertex. Faces carrying an insertion are
    // refused by the admission test for a different reason -- their per-face
    // spontaneous curvature would migrate to the other triangle.
    delta.energy = delta.bending + delta.regularization + delta.areaConstraint +
                   delta.volumeConstraint;
    return true;
}
