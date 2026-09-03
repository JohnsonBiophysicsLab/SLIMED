/**
 * @file Host_force_backend.hpp
 * @brief Runs the force kernels over the flattened mesh, on the CPU.
 *
 * This is the driver the CUDA backend is a translation of: same layout, same
 * kernel bodies, same order of stages, with a for loop where the device has a
 * grid. It exists for two reasons.
 *
 * First, it makes the port testable. The hard part of moving this calculation
 * to a GPU is not the arithmetic -- Patch_kernel.inl already carries that, and
 * it is pinned against the original -- but the flattening: one-ring offsets,
 * child-slot arithmetic, and the scatter/gather that replaces per-thread force
 * buffers. All of that is exercised here, on the host, against the production
 * force evaluation, before nvcc sees any of it.
 *
 * Second, it is the fallback. A build without CUDA, or a run on a machine with
 * no device, still has a working implementation of exactly the same pipeline.
 */
#pragma once

#include <vector>

#include "cuda/Device_mesh_layout.hpp"
#include "cuda/Force_kernels.hpp"

class Mesh;
class PatchRowsFlat;

namespace slimed
{

/**
 * @brief Scratch and results for one force evaluation over a flattened mesh.
 *
 * Buffers are sized on the first evaluate() and reused after that, so a line
 * search that calls it a couple of hundred times allocates nothing.
 */
class HostForceBackend
{
public:
    /**
     * @brief Run the full force evaluation and write the results into `mesh`.
     *
     * Reproduces the sequence Compute_Energy_And_Force() runs for the per-face
     * and per-vertex terms: element area and volume, the totals they sum to,
     * the patch energies and forces, and the regularization term. The scalar
     * energy bookkeeping that follows -- the area and volume constraint
     * energies, the scaffolding term, the boundary ghost handling -- stays
     * with the caller.
     *
     * @param layout Must have been built from this mesh's current topology.
     * @param rows   Must have been built from this mesh's shape functions.
     */
    void evaluate(Mesh &mesh, const DeviceMeshLayout &layout, const PatchRowsFlat &rows);

    /// Total scratch held, for the startup diagnostic.
    std::size_t memory_bytes() const;

private:
    /// Point args at buffers sized for this layout, allocating on first use.
    ForceKernelArgs prepare(const DeviceMeshLayout &layout, const PatchRowsFlat &rows);

    std::vector<double> coords_;
    std::vector<double> coordsRef_;
    std::vector<double> faceArea_;
    std::vector<double> faceVolume_;
    std::vector<double> faceEBend_;
    std::vector<double> faceMeanCurv_;
    std::vector<double> faceNormal_;
    std::vector<double> faceERegular_;
    std::vector<unsigned char> faceDeformCase_;
    std::vector<double> slotForceBend_;
    std::vector<double> slotForceArea_;
    std::vector<double> slotForceVolume_;
    std::vector<double> cornerForceRegular_;
    std::vector<double> vertexForceBend_;
    std::vector<double> vertexForceArea_;
    std::vector<double> vertexForceVolume_;
    std::vector<double> vertexForceRegular_;
};

} // namespace slimed
