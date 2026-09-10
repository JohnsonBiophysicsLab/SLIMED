/**
 * @file Multi_extraordinary_patch.hpp
 * @brief Evaluating a face that has more than one extraordinary corner.
 *
 * Stam's reduction, which IrregularPatchRowTable tabulates, peels three
 * regular children off a *single* extraordinary corner. A face with two or
 * three has no such decomposition -- the three "regular" children stop being
 * regular -- which is why this tree has always rejected them.
 *
 * That was tenable while connectivity was fixed at setup, because a mesh can
 * be generated with its extraordinary vertices isolated. It stops being
 * tenable the moment edges flip: flipping an edge of a hexagonal lattice takes
 * four valence-6 corners to 5, 5, 7, 7, so both new triangles carry three
 * extraordinary corners. Every flip produces the case.
 *
 * The way out is Stam's own observation applied per face rather than per mesh.
 * One Loop subdivision of a face's own control net splits it into four
 * children:
 *
 *            c2                          c2
 *           /  \                        /  \
 *          /    \                     E20--E12
 *         /      \        -->         / \  / \
 *        /        \                  /   \/   \
 *      c0 -------- c1              c0----E01---c1
 *
 * Each corner child keeps one original corner and gains two valence-6 edge
 * points, so it carries at most one extraordinary corner; the centre child
 * carries none. The limit surface does not move, because Loop subdivision is
 * exactly the map that defines it.
 *
 * What makes this local and cheap is that **every control point of every child
 * lies inside the parent's own patch.** A child's control point is either an
 * updated original corner, whose vertex mask reads that corner's one-ring, or
 * the midpoint of an edge incident to a corner, whose edge mask reads the two
 * opposite corners of the faces on that edge -- and those faces are incident
 * to a corner, so their third vertices lie in that corner's one-ring. No
 * two-ring is ever needed.
 *
 * So each child carries a **prolongation matrix** `M` of shape
 * (child width) x K, where K = N0 + N1 + N2 - 6 is the parent's patch width,
 * and the whole evaluation is
 *
 *     Xc = M * X                        the child's control net
 *     Ec, fc = existing_kernel(Xc)      regular, or Stam at the child's valence
 *     E  = sum over children of Ec
 *     f  = sum over children of M^T fc  the chain rule, nothing more
 *
 * Nothing about the kernels changes: a corner child is exactly the kind of
 * patch the tree already evaluates. Only the control net it is handed is a
 * linear image of the parent's rather than the parent's itself. That also
 * keeps the table tiny -- four matrices of at most 14 x 18 per valence
 * triple -- where composing the rows through instead would have stored
 * thousands of 7 x K blocks per triple.
 *
 * There is no Jacobian to apply. The children tile the parent's parameter
 * domain exactly once and the rows returned for a child are with respect to
 * that child's own parameters, so the area element already carries the
 * shrunken metric -- the same argument section 3.1 of
 * `irregular_patch_results.md` makes for why the depth-`d` children need no
 * `4^-d` weight.
 *
 * All of it is pure topology over the valence triple, so it is built once per
 * distinct triple and shared by every face carrying that triple.
 *
 * @see docs/edge_flip_plan.md work package 1
 */

#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <vector>

#include "linalg/Linear_algebra.hpp"

/**
 * @brief The generic one-ring neighbourhood of a triangle, as an explicit
 * vertex and face list.
 *
 * Built for a valence triple rather than read off a mesh, so the Loop masks
 * can be applied generically and the result checked instead of transcribed.
 * The internal numbering is the contract between this and
 * Mesh::build_one_ring_for_face(), which has to list a real face's control net
 * in exactly this order:
 *
 *     0, 1, 2      the three corners, in the face's own winding order
 *     3, 4, 5      the corners opposite the edges (0,1), (1,2), (2,0)
 *     6 ...        corner 0's remaining fan, then corner 1's, then corner 2's
 *
 * Corner `a`'s fan is listed in winding order starting at corner `a+1`, so
 * `cornerFan[a][0]` and `[1]` are the other two corners, `[2]` is the vertex
 * opposite the edge leaving `a` backwards, the middle entries are the extras,
 * and the last is the vertex opposite the edge leaving `a` forwards. At
 * (6, 6, 6) the whole thing is a relabelling of the familiar 12-point patch.
 */
struct GenericFacePatch
{
    std::array<int, 3> valence{{0, 0, 0}};
    int nVertices = 0;                     ///< K = N0 + N1 + N2 - 6
    std::vector<std::array<int, 3>> faces; ///< N0 + N1 + N2 - 5, consistently wound
    /// Cyclic neighbour list of each corner, in winding order from the next corner.
    std::array<std::vector<int>, 3> cornerFan;
};

/**
 * @brief Build the generic patch for a valence triple.
 *
 * @throw std::invalid_argument if any valence is outside [4, 8].
 */
GenericFacePatch build_generic_face_patch(int valence0, int valence1, int valence2);

/**
 * @brief Prolongation matrices for every valence triple a mesh presents.
 *
 * Built lazily -- a mesh with no multi-extraordinary faces never pays for
 * it -- and keyed on the triple alone, so a fluid membrane with thousands of
 * such faces holds a few dozen entries. A face resolves to an entry index once,
 * when it is classified, so the force loop never touches the map.
 */
class MultiPatchTable
{
public:
    /// One child of the subdivided face.
    struct Child
    {
        /// Valence of the child's extraordinary corner; 6 when it is regular.
        int valence = 6;
        /// Control points the child's own evaluation needs: valence + 6, or 12.
        int nControl = 0;
        /// First double of this child's nControl x K prolongation matrix,
        /// row-major, indexing into MultiPatchTable::data().
        std::size_t offset = 0;
    };

    /// Everything a face with one valence triple needs.
    struct Entry
    {
        std::array<int, 3> valence{{0, 0, 0}};
        int nControl = 0;              ///< K, the parent patch width
        std::array<Child, 4> children; ///< corner 0, corner 1, corner 2, centre
    };

    /**
     * @brief Index of the entry for this triple, building it if it is new.
     *
     * Not thread safe: call it while classifying faces, never from inside the
     * force loop.
     *
     * @throw std::invalid_argument if a valence is outside [4, 8].
     * @throw std::logic_error if the generic patch does not subdivide as
     *        expected. That would mean this construction is wrong, not the
     *        mesh, so it is not something a caller can handle.
     */
    int ensure(int valence0, int valence1, int valence2);

    bool empty() const { return entries_.empty(); }
    int size() const { return static_cast<int>(entries_.size()); }
    const Entry &entry(int index) const { return entries_.at(index); }

    /// All prolongation matrices back to back; Child::offset indexes into this.
    const double *data() const { return buffer_.empty() ? nullptr : buffer_.data(); }
    std::size_t memory_bytes() const { return buffer_.size() * sizeof(double); }

    void clear();

private:
    std::vector<Entry> entries_;
    std::map<std::array<int, 3>, int> index_;
    std::vector<double> buffer_;
};

/**
 * @brief The child's control net: Xc = M * X.
 *
 * @param prolongation  nChild x nParent, row-major.
 * @param parentCoords  nParent x 3, row-major.
 * @param childCoords   nChild x 3, row-major. Overwritten.
 */
inline void multi_patch_prolong(const double *prolongation, const double *parentCoords,
                                int nParent, int nChild, double *childCoords)
{
    for (int row = 0; row < nChild; row++)
    {
        const double *weights = prolongation + static_cast<std::size_t>(row) * nParent;
        double accumulated[3] = {0.0, 0.0, 0.0};
        for (int column = 0; column < nParent; column++)
        {
            const double weight = weights[column];
            if (weight == 0.0)
            {
                continue;
            }
            accumulated[0] += weight * parentCoords[column * 3 + 0];
            accumulated[1] += weight * parentCoords[column * 3 + 1];
            accumulated[2] += weight * parentCoords[column * 3 + 2];
        }
        childCoords[row * 3 + 0] = accumulated[0];
        childCoords[row * 3 + 1] = accumulated[1];
        childCoords[row * 3 + 2] = accumulated[2];
    }
}

/**
 * @brief Push a child's force back onto the parent's control points:
 * f += M^T * fc.
 *
 * The chain rule and nothing more. The energy of a child is a function of
 * Xc = M * X, so dE/dX = M^T dE/dXc, and the kernel already returns the
 * negative gradient.
 *
 * @param parentForce  nParent x 3, row-major. Accumulated into.
 */
inline void multi_patch_scatter(const double *prolongation, const double *childForce, int nParent,
                                int nChild, double *parentForce)
{
    for (int row = 0; row < nChild; row++)
    {
        const double *weights = prolongation + static_cast<std::size_t>(row) * nParent;
        for (int column = 0; column < nParent; column++)
        {
            const double weight = weights[column];
            if (weight == 0.0)
            {
                continue;
            }
            parentForce[column * 3 + 0] += weight * childForce[row * 3 + 0];
            parentForce[column * 3 + 1] += weight * childForce[row * 3 + 1];
            parentForce[column * 3 + 2] += weight * childForce[row * 3 + 2];
        }
    }
}
