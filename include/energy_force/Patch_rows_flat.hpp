/**
 * @file Patch_rows_flat.hpp
 * @brief Shape-function rows repacked into flat, contiguous double arrays.
 *
 * The shape functions and the irregular-child rows are built once at mesh
 * setup and never change afterwards, but they are stored as std::vector<Matrix>
 * -- one heap-allocated gsl_matrix per quadrature point per child. Reading
 * them element by element inside the force loop was a large part of the
 * profile, and a pointer-chasing structure cannot be handed to a GPU at all.
 *
 * This repacks the same numbers into contiguous row-major blocks that
 * element_energy_force_patch_pod() consumes directly, and that a future device
 * port can cudaMemcpy in one call. It is a cache: build() copies, it does not
 * take ownership of or alias the Matrix data.
 */
#pragma once

#include <cstddef>
#include <vector>

#include "linalg/Linear_algebra.hpp"
#include "mesh/Irregular_patch_rows.hpp"

/**
 * @brief Flat storage for every shape-function block a patch evaluation reads.
 */
class PatchRowsFlat
{
public:
    /**
     * @brief Repack the regular shape functions, the irregular child rows and
     * the quadrature weights.
     *
     * Safe to call again; it rebuilds from scratch. Throws
     * std::invalid_argument if the regular shape functions are not the
     * expected 7 x 12, or if the quadrature weights do not have one entry per
     * shape function.
     *
     * @param regularShapeFunctions One 7 x 12 block per quadrature point.
     * @param irregularRows         The built irregular-child table. An empty
     *                              table is allowed: a mesh with no irregular
     *                              faces never reads those blocks.
     * @param gaussQuadratureCoeff  Quadrature weights, one row per point.
     */
    void build(const std::vector<Matrix> &regularShapeFunctions,
               const IrregularPatchRowTable &irregularRows,
               const Matrix &gaussQuadratureCoeff);

    bool empty() const { return regular_.empty(); }

    /// Quadrature points per patch evaluation.
    int nSamples() const { return nSamples_; }

    /// nSamples() x 7 x 12 blocks, contiguous and row-major.
    const double *regular() const { return regular_.data(); }

    /// nSamples() weights.
    const double *gaussCoeff() const { return gaussCoeff_.data(); }

    /**
     * @brief The nSamples() x 7 x (valence + 6) blocks for one irregular child.
     *
     * @throw std::out_of_range if the table was built without this slot.
     */
    const double *child(int valence, int depth, int childIndex) const;

    /// Total bytes held, for the startup diagnostic.
    std::size_t memory_bytes() const;

private:
    int slot_index(int valence, int depth, int childIndex) const;

    int nSamples_ = 0;
    int depth_ = 0;
    std::vector<double> regular_;
    /// All irregular blocks back to back; offsets_ indexes into this.
    std::vector<double> irregular_;
    /// One offset per (valence, depth, child) slot, in the same order the
    /// IrregularPatchRowTable uses. kNoSlot marks a slot that was not built.
    std::vector<std::size_t> offsets_;
    std::vector<double> gaussCoeff_;
};
