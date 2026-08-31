/**
 * @file Cuda_force_backend.hpp
 * @brief The force evaluation on an NVIDIA GPU.
 *
 * The device counterpart of HostForceBackend, and deliberately little more
 * than that: same DeviceMeshLayout, same PatchRowsFlat, same kernel bodies
 * from Force_kernels.inl, with a grid where the host has a for loop. Nothing
 * about the physics is restated here, so nothing about the physics can drift.
 *
 * This header is safe to include from a build without CUDA -- it declares the
 * class but nothing CUDA-specific. Whether an implementation is linked in is
 * what cuda_backend_available() answers.
 */
#pragma once

#include <cstddef>
#include <string>

#include "cuda/Device_mesh_layout.hpp"

class Mesh;
class PatchRowsFlat;

namespace slimed
{

/// True if this binary was built with CUDA support compiled in. Says nothing
/// about whether a usable device is present -- CudaForceBackend's constructor
/// answers that.
bool cuda_backend_available();

/**
 * @brief Device-resident mesh state and the kernels that act on it.
 *
 * Topology and shape functions are uploaded once by upload_topology() and
 * reused for the life of the object; only vertex coordinates cross the bus per
 * evaluation. Rebuild after a change in connectivity.
 *
 * @throw std::runtime_error from any method if CUDA is not compiled in, if no
 *        device is visible, or if a CUDA call fails. Callers that want to fall
 *        back to HostForceBackend should construct inside a try block.
 */
class CudaForceBackend
{
public:
    CudaForceBackend();
    ~CudaForceBackend();

    CudaForceBackend(const CudaForceBackend &) = delete;
    CudaForceBackend &operator=(const CudaForceBackend &) = delete;

    /// Name and compute capability of the device in use, for the startup log.
    std::string device_description() const;

    /// Upload the parts that do not change: topology and shape-function rows.
    void upload_topology(const DeviceMeshLayout &layout, const PatchRowsFlat &rows);

    /**
     * @brief Run the force evaluation on the device and write it into `mesh`.
     *
     * Same contract as HostForceBackend::evaluate(): element area and volume
     * and their totals, patch energies and forces, and the regularization
     * term. The scalar energy bookkeeping that follows stays with the caller.
     */
    void evaluate(Mesh &mesh, const DeviceMeshLayout &layout, const PatchRowsFlat &rows);

    /// Device memory currently held, for the startup diagnostic.
    std::size_t device_memory_bytes() const;

private:
    struct Impl;
    Impl *impl_ = nullptr;
};

} // namespace slimed
