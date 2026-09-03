/**
 * @file Cuda_force_backend_stub.cpp
 * @brief Stands in for Cuda_force_backend.cu in a build without CUDA.
 *
 * Keeps CudaForceBackend linkable everywhere, so callers can be written as
 * "try the device, fall back to the host" without an #ifdef at every use.
 * Every entry point reports that CUDA was not compiled in, rather than
 * pretending to work.
 */
#include "cuda/Cuda_force_backend.hpp"

#include <stdexcept>

namespace slimed
{
namespace
{
std::runtime_error not_built()
{
    return std::runtime_error(
        "[CudaForceBackend] this binary was built without CUDA support. Reconfigure with "
        "-DSLIMED_ENABLE_CUDA=ON and a CUDA toolkit on PATH, or use HostForceBackend.");
}
} // namespace

bool cuda_backend_available() { return false; }

CudaForceBackend::CudaForceBackend() { throw not_built(); }
CudaForceBackend::~CudaForceBackend() = default;
std::string CudaForceBackend::device_description() const { return "none (built without CUDA)"; }
void CudaForceBackend::upload_topology(const DeviceMeshLayout &, const PatchRowsFlat &)
{
    throw not_built();
}
void CudaForceBackend::evaluate(Mesh &, const DeviceMeshLayout &, const PatchRowsFlat &)
{
    throw not_built();
}
std::size_t CudaForceBackend::device_memory_bytes() const { return 0; }

} // namespace slimed
