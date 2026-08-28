/**
 * @file Subdivision_matrices.hpp
 * @brief Generate Loop subdivision matrices for an irregular patch of any
 * supported valence, instead of writing them out as literals.
 *
 * The reference implementation this tree inherited writes the valence-5
 * matrices out as literal arrays. Every entry encodes three things at once --
 * the Loop masks, the canonical ordering of the control points, and the
 * ordering of the children -- and none of the three is recoverable from the
 * literal, which is why the port stopped at one valence.
 *
 * Stam's reduction is uniform in the valence N:
 *
 *     control points in the patch      K(N) = N + 6
 *     points after one subdivision     K(N) + 6 = N + 12
 *     subdivision matrix               abar : (N+12) x (N+6)
 *     three regular children           p1, p2, p3 : 12 x (N+12)
 *     one irregular child              p4 : (N+6) x (N+12)
 *
 * At N = 5 those shapes are 17x11, 12x17 and 11x17 -- exactly the matrices in
 * get_subdivision_matrices(). The existing code is one instance of a family,
 * written longhand.
 *
 * Everything here is pure topology: it depends only on N, never on vertex
 * coordinates, so it can be built once at startup and shared.
 *
 * @note WP2 of docs/irregular_patch_valence_4_to_8_plan.md.
 */
#pragma once

#include <array>
#include <vector>

#include "linalg/Linear_algebra.hpp"

/// Valences this generator supports. 6 is regular and handled by the direct kernel.
constexpr int kMinIrregularValence = 4;
constexpr int kMaxIrregularValence = 8;

/**
 * @brief Warren's vertex mask weight for a valence-N vertex.
 *
 * Loop's original rule and Warren's differ; this tree and the reference POC
 * both use Warren's `3 / (8N)`, and get_subdivision_matrices() is built from
 * it, so the N=5 parity gate pins this choice. Exposed as one function so the
 * alternative is a one-line swap rather than a rewrite.
 *
 * At N = 6 it gives 1/16, the regular value, so the all-regular baseline is
 * untouched either way. Note that 3/(8N) is only the standard rule for N > 3,
 * which is one reason valence 3 is out of scope.
 */
inline double loop_vertex_weight(int valence)
{
    return 3.0 / (8.0 * valence);
}

/**
 * @brief The canonical one-ring patch around a single extraordinary vertex.
 *
 * Built as an explicit integer vertex/face list so the Loop masks can be
 * applied generically and the result checked, rather than transcribed.
 *
 * Vertex numbering, using the d1..d12 names the regular patch already uses:
 *
 *     0        d4, the extraordinary corner, valence N
 *     1..N     its fan in cyclic order: d7, d8, d5, <N-4 extras>, d3
 *     N+1      d6
 *     N+2      d9
 *     N+3      d10
 *     N+4      d11
 *     N+5      d12
 *
 * The patch triangle is (d4, d7, d8) = (0, 1, 2). d7 and d8 are both at
 * valence 6; only d4, d7 and d8 have complete one-rings inside the patch, so
 * only those three get vertex points under subdivision.
 */
struct CanonicalPatch
{
    int valence = 0;
    int nVertices = 0;                        ///< N + 6
    std::vector<std::array<int, 3>> faces;    ///< N + 7, consistently wound
    std::vector<std::array<double, 2>> layout; ///< canonical 2D positions, for child ordering
};

/// Build the canonical patch topology for the given valence.
CanonicalPatch build_canonical_patch(int valence);

/**
 * @brief The four matrices describing one subdivision step of an irregular patch.
 *
 * `abar` maps the N+6 control points to the N+12 points of the subdivided
 * patch. `p1`, `p2` and `p3` select the three regular 12-point children;
 * `p4` selects the smaller irregular child, which still carries the
 * extraordinary corner and is what the recursion descends into.
 */
struct SubdivisionMatrices
{
    int valence = 0;
    Matrix abar; ///< (N+12) x (N+6)
    Matrix p1;   ///< 12 x (N+12)
    Matrix p2;   ///< 12 x (N+12)
    Matrix p3;   ///< 12 x (N+12)
    Matrix p4;   ///< (N+6) x (N+12)
};

/**
 * @brief Generate the subdivision matrices for the given valence.
 *
 * @throw std::invalid_argument if the valence is outside [4, 8]. Out-of-range
 * input is an error, not a fallback: silently approximating an unsupported
 * topology is the failure mode this whole work package exists to remove.
 */
SubdivisionMatrices build_subdivision_matrices(int valence);
