/**
 * @file Patch_kind.hpp
 * @brief What kind of subdivision patch a face carries.
 *
 * Split out of Mesh.hpp so that Face can hold its own classification: the
 * width of oneRingVertices used to be enough to tell the evaluator what it was
 * looking at (`valence = nOneRingVertices - 6`), and once faces can carry more
 * than one extraordinary corner it is not. A 5/5/7 face and a valence-5 face
 * both have an 11-point control net.
 *
 * @see docs/edge_flip_plan.md sections 1.4 and 3.3
 */

#pragma once

#include <string>

/**
 * @brief Which subdivision patch, if any, a face carries.
 *
 * The classification used to be inlined in set_one_ring_vertices_sorted() and
 * its only outputs were "build a one-ring" or "throw". A Monte Carlo flip has
 * to ask the same question speculatively -- would the mesh still be evaluable
 * if this edge were flipped? -- so the question and the answer are named here
 * and the throw is left to the caller that wants one.
 */
enum class PatchKind
{
    /// Ghost face outside the boundary; takes no part in any calculation.
    Ghost,
    /// A corner is not interior, so there is no complete one-ring and no limit
    /// surface. A property of the mesh, not an error: oneRingVertices stays
    /// empty on purpose.
    Boundary,
    /// All three corners at valence 6; the closed-form quartic box spline.
    Regular,
    /// Exactly one corner in [kMinIrregularValence, kMaxIrregularValence]
    /// other than 6, with the other two at exactly 6. Stam's reduction, which
    /// is what IrregularPatchRowTable tabulates.
    SingleExtraordinary,
    /**
     * @brief More than one extraordinary corner, every valence still in range.
     *
     * Unavoidable the moment edges flip: flipping an edge of a hexagonal
     * lattice turns four valence-6 corners into 5, 5, 7, 7, so both new
     * triangles carry three extraordinary corners. Evaluated by subdividing
     * the face's own control net once, which splits it into four children each
     * carrying at most one extraordinary corner, and pushing those through the
     * paths that already exist.
     */
    MultiExtraordinary,
    /// A valence outside the supported range, or a fan that does not close.
    Inadmissible,
};

/**
 * @brief The classification of one face, and the rotation that anchors it.
 */
struct PatchClass
{
    PatchKind kind = PatchKind::Inadmissible;

    /**
     * @brief Which corner (0, 1 or 2 of Face::adjacentVertices) anchors the
     * patch: the extraordinary one when there is exactly one, corner 0 when
     * the face is regular or carries several extraordinary corners, and -1
     * when the face carries no patch.
     *
     * Reading the corners as a rotation starting here preserves the face
     * winding, which sort_vertices_on_faces() has already made consistent, so
     * no per-face winding decision is ever taken.
     */
    int anchor = -1;

    /// Valences of the three corners, in the face's own order (not rotated).
    int valence[3] = {0, 0, 0};

    /// Why the face was rejected; empty unless kind is Inadmissible.
    std::string why;

    /// Whether this face carries a patch the energy kernel can evaluate.
    bool has_evaluable_patch() const
    {
        return kind == PatchKind::Regular || kind == PatchKind::SingleExtraordinary ||
               kind == PatchKind::MultiExtraordinary;
    }
};
