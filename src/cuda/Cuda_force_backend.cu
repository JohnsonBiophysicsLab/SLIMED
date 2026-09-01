/**
 * @file Cuda_force_backend.cu
 * @brief CUDA kernels and device memory management for the force evaluation.
 *
 * Compiled only when SLIMED_ENABLE_CUDA is on; Cuda_force_backend_stub.cpp
 * stands in otherwise. Nothing here restates the physics. Each __global__ is
 * three lines -- work out an index, bounds-check it, call the shared body from
 * Force_kernels.inl -- and the rest of the file is allocation, transfer and
 * launch.
 *
 * That is the whole point of the arrangement: the bodies are the same
 * functions HostForceBackend runs and tests/test_device_mesh_layout.cpp pins
 * against the production path, so the only thing that can be wrong here is the
 * plumbing.
 */
#include "cuda/Cuda_force_backend.hpp"

#include <cuda_runtime.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "cuda/Force_kernels.hpp"
#include "energy_force/Patch_rows_flat.hpp"
#include "mesh/Mesh.hpp"

namespace slimed
{

namespace
{
/// Threads per block. 128 suits a kernel with this register pressure: the
/// patch body holds a dozen or so 3-vectors plus a 14x3 control net in
/// registers and local memory, so a larger block trades occupancy for spills.
constexpr int kBlockSize = 128;

constexpr int grid_for(int n) { return (n + kBlockSize - 1) / kBlockSize; }

/// Throw with the CUDA error text rather than letting a failure go unnoticed.
/// Silently ignoring a cudaMemcpy return is how a GPU port ends up producing
/// plausible garbage.
void check(cudaError_t status, const char *what)
{
    if (status != cudaSuccess)
    {
        throw std::runtime_error(std::string("[CudaForceBackend] ") + what + ": " +
                                 cudaGetErrorString(status));
    }
}

/// Check that the last launch was accepted and, after a sync, that it ran.
void check_launch(const char *what)
{
    check(cudaGetLastError(), what);
    check(cudaDeviceSynchronize(), what);
}

/// A typed device allocation that frees itself.
template <typename T>
class DeviceBuffer
{
public:
    ~DeviceBuffer() { release(); }

    void resize(std::size_t count)
    {
        if (count <= capacity_)
        {
            count_ = count;
            return;
        }
        release();
        if (count > 0)
        {
            check(cudaMalloc(&data_, count * sizeof(T)), "cudaMalloc");
        }
        capacity_ = count;
        count_ = count;
    }

    void upload(const T *host, std::size_t count)
    {
        resize(count);
        if (count > 0)
        {
            check(cudaMemcpy(data_, host, count * sizeof(T), cudaMemcpyHostToDevice),
                  "cudaMemcpy host to device");
        }
    }

    void download(T *host) const
    {
        if (count_ > 0)
        {
            check(cudaMemcpy(host, data_, count_ * sizeof(T), cudaMemcpyDeviceToHost),
                  "cudaMemcpy device to host");
        }
    }

    T *data() { return data_; }
    const T *data() const { return data_; }
    std::size_t count() const { return count_; }
    std::size_t bytes() const { return capacity_ * sizeof(T); }

private:
    void release()
    {
        if (data_ != nullptr)
        {
            cudaFree(data_);
            data_ = nullptr;
        }
        capacity_ = 0;
        count_ = 0;
    }

    T *data_ = nullptr;
    std::size_t capacity_ = 0;
    std::size_t count_ = 0;
};

__global__ void area_volume_kernel(ForceKernelArgs args)
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if (face < args.nFaces)
    {
        area_volume_for_face(args, face);
    }
}

__global__ void patch_force_kernel(ForceKernelArgs args)
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if (face < args.nFaces)
    {
        patch_force_for_face(args, face);
    }
}

__global__ void regularization_kernel(ForceKernelArgs args)
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if (face < args.nFaces)
    {
        regularization_for_face(args, face);
    }
}

__global__ void gather_kernel(ForceKernelArgs args)
{
    const int vertex = blockIdx.x * blockDim.x + threadIdx.x;
    if (vertex < args.nVertices)
    {
        gather_force_for_vertex(args, vertex);
        gather_regularization_for_vertex(args, vertex);
    }
}
} // namespace

/// Every device allocation, plus the host staging buffers the results land in.
struct CudaForceBackend::Impl
{
    int device = 0;
    cudaDeviceProp properties{};
    bool topologyUploaded = false;

    DeviceBuffer<FacePatchDescriptor> descriptors;
    DeviceBuffer<int> oneRingIndices;
    DeviceBuffer<int> faceValence;
    DeviceBuffer<double> faceSpontCurvature;
    DeviceBuffer<unsigned char> faceIsGhost;
    DeviceBuffer<unsigned char> faceIsBoundary;
    DeviceBuffer<int> faceCorners;
    DeviceBuffer<int> vertexSlotOffsets;
    DeviceBuffer<int> vertexSlots;
    DeviceBuffer<int> vertexCornerOffsets;
    DeviceBuffer<int> vertexCorners;

    DeviceBuffer<double> regularRows;
    DeviceBuffer<double> gaussCoeff;
    DeviceBuffer<double> irregularRows;
    DeviceBuffer<std::size_t> irregularOffsets;

    DeviceBuffer<double> coords;
    DeviceBuffer<double> coordsRef;

    DeviceBuffer<double> faceArea;
    DeviceBuffer<double> faceVolume;
    DeviceBuffer<double> faceEBend;
    DeviceBuffer<double> faceMeanCurv;
    DeviceBuffer<double> faceNormal;
    DeviceBuffer<double> faceERegular;
    DeviceBuffer<unsigned char> faceDeformCase;
    DeviceBuffer<double> slotForceBend;
    DeviceBuffer<double> slotForceArea;
    DeviceBuffer<double> slotForceVolume;
    DeviceBuffer<double> cornerForceRegular;
    DeviceBuffer<double> vertexForceBend;
    DeviceBuffer<double> vertexForceArea;
    DeviceBuffer<double> vertexForceVolume;
    DeviceBuffer<double> vertexForceRegular;

    std::vector<double> hostCoords;
    std::vector<double> hostCoordsRef;
    std::vector<double> hostFaceArea;
    std::vector<double> hostFaceVolume;
    std::vector<double> hostFaceEBend;
    std::vector<double> hostFaceMeanCurv;
    std::vector<double> hostFaceNormal;
    std::vector<double> hostFaceERegular;
    std::vector<unsigned char> hostFaceDeformCase;
    std::vector<double> hostVertexForceBend;
    std::vector<double> hostVertexForceArea;
    std::vector<double> hostVertexForceVolume;
    std::vector<double> hostVertexForceRegular;

    /// Point an args struct at the device buffers.
    ForceKernelArgs args(const DeviceMeshLayout &layout, const PatchRowsFlat &rows)
    {
        ForceKernelArgs a;
        a.descriptors = descriptors.data();
        a.oneRingIndices = oneRingIndices.data();
        a.faceValence = faceValence.data();
        a.faceSpontCurvature = faceSpontCurvature.data();
        a.faceIsGhost = faceIsGhost.data();
        a.faceIsBoundary = faceIsBoundary.data();
        a.faceCorners = faceCorners.data();
        a.vertexSlotOffsets = vertexSlotOffsets.data();
        a.vertexSlots = vertexSlots.data();
        a.vertexCornerOffsets = vertexCornerOffsets.data();
        a.vertexCorners = vertexCorners.data();
        a.nFaces = layout.nFaces();
        a.nVertices = layout.nVertices();

        a.regularRows = regularRows.data();
        a.gaussCoeff = gaussCoeff.data();
        a.nSamples = rows.nSamples();
        a.irregularRows = irregularRows.count() > 0 ? irregularRows.data() : nullptr;
        a.irregularOffsets = irregularOffsets.count() > 0 ? irregularOffsets.data() : nullptr;
        a.irregularDepth = rows.depth();

        a.coords = coords.data();
        a.coordsRef = coordsRef.data();

        a.faceArea = faceArea.data();
        a.faceVolume = faceVolume.data();
        a.faceEBend = faceEBend.data();
        a.faceMeanCurv = faceMeanCurv.data();
        a.faceNormal = faceNormal.data();
        a.faceERegular = faceERegular.data();
        a.faceDeformCase = faceDeformCase.data();
        a.slotForceBend = slotForceBend.data();
        a.slotForceArea = slotForceArea.data();
        a.slotForceVolume = slotForceVolume.data();
        a.cornerForceRegular = cornerForceRegular.data();
        a.vertexForceBend = vertexForceBend.data();
        a.vertexForceArea = vertexForceArea.data();
        a.vertexForceVolume = vertexForceVolume.data();
        a.vertexForceRegular = vertexForceRegular.data();
        return a;
    }
};

bool cuda_backend_available() { return true; }

CudaForceBackend::CudaForceBackend() : impl_(new Impl())
{
    int deviceCount = 0;
    const cudaError_t status = cudaGetDeviceCount(&deviceCount);
    if (status != cudaSuccess || deviceCount == 0)
    {
        delete impl_;
        impl_ = nullptr;
        throw std::runtime_error(
            std::string("[CudaForceBackend] no CUDA device is available: ") +
            (status == cudaSuccess ? "the driver reports zero devices"
                                   : cudaGetErrorString(status)));
    }
    check(cudaGetDevice(&impl_->device), "cudaGetDevice");
    check(cudaGetDeviceProperties(&impl_->properties, impl_->device), "cudaGetDeviceProperties");
}

CudaForceBackend::~CudaForceBackend() { delete impl_; }

std::string CudaForceBackend::device_description() const
{
    if (impl_ == nullptr)
    {
        return "none";
    }
    return std::string(impl_->properties.name) + " (compute " +
           std::to_string(impl_->properties.major) + "." +
           std::to_string(impl_->properties.minor) + ", " +
           std::to_string(impl_->properties.multiProcessorCount) + " SMs)";
}

void CudaForceBackend::upload_topology(const DeviceMeshLayout &layout, const PatchRowsFlat &rows)
{
    Impl &impl = *impl_;
    const std::size_t nFaces = static_cast<std::size_t>(layout.nFaces());
    const std::size_t nVertices = static_cast<std::size_t>(layout.nVertices());
    const std::size_t nSlots = static_cast<std::size_t>(layout.nSlots());

    impl.descriptors.upload(layout.descriptors(), nFaces);
    impl.oneRingIndices.upload(layout.oneRingIndices(), nSlots);
    impl.faceValence.upload(layout.faceValence(), nFaces);
    impl.faceSpontCurvature.upload(layout.faceSpontCurvature(), nFaces);
    impl.faceIsGhost.upload(layout.faceIsGhost(), nFaces);
    impl.faceIsBoundary.upload(layout.faceIsBoundary(), nFaces);
    impl.faceCorners.upload(layout.faceCorners(), nFaces * 3);
    impl.vertexSlotOffsets.upload(layout.vertexSlotOffsets(), nVertices + 1);
    impl.vertexSlots.upload(layout.vertexSlots(), nSlots);
    impl.vertexCornerOffsets.upload(layout.vertexCornerOffsets(), nVertices + 1);
    impl.vertexCorners.upload(layout.vertexCorners(), nFaces * 3);

    const std::size_t nSamples = static_cast<std::size_t>(rows.nSamples());
    impl.regularRows.upload(rows.regular(), nSamples * slimed::kShapeRows * 12);
    impl.gaussCoeff.upload(rows.gaussCoeff(), nSamples);
    // A mesh with no irregular faces has no table, and the kernels never index
    // one; leaving both buffers empty makes rows_for_face() return null for a
    // kind it will never be asked about.
    if (rows.irregular() != nullptr)
    {
        impl.irregularRows.upload(rows.irregular(), rows.irregular_count());
        impl.irregularOffsets.upload(rows.irregularOffsets(), rows.irregular_slot_count());
    }

    // Outputs and scratch: allocated, never uploaded. Every entry is written
    // before it is read -- patch_force_for_face() zeroes the slots of a face
    // it skips for exactly that reason.
    impl.coords.resize(nVertices * 3);
    impl.coordsRef.resize(nVertices * 3);
    impl.faceArea.resize(nFaces);
    impl.faceVolume.resize(nFaces);
    impl.faceEBend.resize(nFaces);
    impl.faceMeanCurv.resize(nFaces);
    impl.faceNormal.resize(nFaces * 3);
    impl.faceERegular.resize(nFaces);
    impl.faceDeformCase.resize(nFaces);
    impl.slotForceBend.resize(nSlots * 3);
    impl.slotForceArea.resize(nSlots * 3);
    impl.slotForceVolume.resize(nSlots * 3);
    impl.cornerForceRegular.resize(nFaces * 3 * 3);
    impl.vertexForceBend.resize(nVertices * 3);
    impl.vertexForceArea.resize(nVertices * 3);
    impl.vertexForceVolume.resize(nVertices * 3);
    impl.vertexForceRegular.resize(nVertices * 3);

    impl.hostFaceArea.resize(nFaces);
    impl.hostFaceVolume.resize(nFaces);
    impl.hostFaceEBend.resize(nFaces);
    impl.hostFaceMeanCurv.resize(nFaces);
    impl.hostFaceNormal.resize(nFaces * 3);
    impl.hostFaceERegular.resize(nFaces);
    impl.hostFaceDeformCase.resize(nFaces);
    impl.hostVertexForceBend.resize(nVertices * 3);
    impl.hostVertexForceArea.resize(nVertices * 3);
    impl.hostVertexForceVolume.resize(nVertices * 3);
    impl.hostVertexForceRegular.resize(nVertices * 3);

    impl.topologyUploaded = true;
}

void CudaForceBackend::evaluate(Mesh &mesh, const DeviceMeshLayout &layout,
                                const PatchRowsFlat &rows)
{
    Impl &impl = *impl_;
    if (!impl.topologyUploaded)
    {
        upload_topology(layout, rows);
    }

    layout.gather_coordinates(mesh, impl.hostCoords);
    layout.gather_reference_coordinates(mesh, impl.hostCoordsRef);
    impl.coords.upload(impl.hostCoords.data(), impl.hostCoords.size());
    impl.coordsRef.upload(impl.hostCoordsRef.data(), impl.hostCoordsRef.size());

    ForceKernelArgs args = impl.args(layout, rows);
    Param &param = mesh.param;

    // Stage 1: element area and volume.
    area_volume_kernel<<<grid_for(args.nFaces), kBlockSize>>>(args);
    check_launch("area_volume_kernel");

    // The totals are summed on the host, in face order.
    //
    // A device reduction would save this round trip, and it is the obvious next
    // optimisation -- but its summation order would depend on the block count,
    // so the GPU and CPU would disagree in the last bits and neither would
    // reproduce the other. Until the constraint forces below stop reading these
    // totals mid-evaluation, matching the CPU exactly is worth one transfer of
    // nFaces doubles.
    impl.faceArea.download(impl.hostFaceArea.data());
    impl.faceVolume.download(impl.hostFaceVolume.data());
    double area = 0.0;
    double volume = 0.0;
    for (int face = 0; face < args.nFaces; face++)
    {
        mesh.faces[face].elementArea = impl.hostFaceArea[face];
        mesh.faces[face].elementVolume = impl.hostFaceVolume[face];
        if (!layout.faceIsGhost()[face])
        {
            area += impl.hostFaceArea[face];
            volume += impl.hostFaceVolume[face];
        }
    }
    param.area = area;
    param.vol = volume;

    // Stage 2: energies and forces. No reference quantity means no constraint;
    // the guards match Compute_Energy_And_Force()'s.
    args.patchParams.kCurv = param.kCurv;
    args.patchParams.uSurfPerArea = (param.area0 == 0.0) ? 0.0 : param.uSurf / param.area0;
    args.patchParams.area = param.area;
    args.patchParams.area0 = param.area0;
    args.patchParams.uVol = (param.vol0 == 0.0) ? 0.0 : param.uVol / param.vol0;
    args.patchParams.vol = param.vol;
    args.patchParams.vol0 = param.vol0;
    args.kCurv = param.kCurv;
    args.gamaShape = param.gamaShape;
    args.gamaArea = param.gamaArea;
    args.usingRpi = param.usingRpi;

    patch_force_kernel<<<grid_for(args.nFaces), kBlockSize>>>(args);
    check_launch("patch_force_kernel");
    regularization_kernel<<<grid_for(args.nFaces), kBlockSize>>>(args);
    check_launch("regularization_kernel");
    // The gather must see every scatter, which the sync inside check_launch
    // above guarantees.
    gather_kernel<<<grid_for(args.nVertices), kBlockSize>>>(args);
    check_launch("gather_kernel");

    impl.faceEBend.download(impl.hostFaceEBend.data());
    impl.faceMeanCurv.download(impl.hostFaceMeanCurv.data());
    impl.faceNormal.download(impl.hostFaceNormal.data());
    impl.faceERegular.download(impl.hostFaceERegular.data());
    impl.faceDeformCase.download(impl.hostFaceDeformCase.data());
    impl.vertexForceBend.download(impl.hostVertexForceBend.data());
    impl.vertexForceArea.download(impl.hostVertexForceArea.data());
    impl.vertexForceVolume.download(impl.hostVertexForceVolume.data());
    impl.vertexForceRegular.download(impl.hostVertexForceRegular.data());

    // Stage 3: write the results back into the mesh objects. Identical to
    // HostForceBackend::evaluate()'s final stage.
    int shapeDeformCount = 0;
    int areaDeformCount = 0;
    for (int face = 0; face < args.nFaces; face++)
    {
        Face &meshFace = mesh.faces[face];
        meshFace.energy.energyCurvature = impl.hostFaceEBend[face];
        meshFace.meanCurvature = impl.hostFaceMeanCurv[face];
        meshFace.energy.energyRegularization = impl.hostFaceERegular[face];
        meshFace.normVector.free();
        meshFace.normVector = mat_calloc(3, 1);
        for (int axis = 0; axis < 3; axis++)
        {
            meshFace.normVector.set(axis, 0, impl.hostFaceNormal[face * 3 + axis]);
        }
        if (impl.hostFaceDeformCase[face] == 1)
        {
            areaDeformCount++;
        }
        else if (impl.hostFaceDeformCase[face] == 2)
        {
            shapeDeformCount++;
        }
    }
    param.deformationCount.shapeDeformCount = shapeDeformCount;
    param.deformationCount.areaDeformCount = areaDeformCount;
    param.deformationCount.noDeformCount =
        static_cast<int>(mesh.faces.size()) - shapeDeformCount - areaDeformCount;

    for (int vertex = 0; vertex < args.nVertices; vertex++)
    {
        Force &force = mesh.vertices[vertex].force;
        for (int axis = 0; axis < 3; axis++)
        {
            force.forceCurvature.set(axis, 0, impl.hostVertexForceBend[vertex * 3 + axis]);
            force.forceArea.set(axis, 0, impl.hostVertexForceArea[vertex * 3 + axis]);
            force.forceVolume.set(axis, 0, impl.hostVertexForceVolume[vertex * 3 + axis]);
            force.forceRegularization.set(axis, 0,
                                          impl.hostVertexForceRegular[vertex * 3 + axis]);
        }
    }
}

std::size_t CudaForceBackend::device_memory_bytes() const
{
    if (impl_ == nullptr)
    {
        return 0;
    }
    const Impl &i = *impl_;
    return i.descriptors.bytes() + i.oneRingIndices.bytes() + i.faceValence.bytes() +
           i.faceSpontCurvature.bytes() + i.faceIsGhost.bytes() + i.faceIsBoundary.bytes() +
           i.faceCorners.bytes() + i.vertexSlotOffsets.bytes() + i.vertexSlots.bytes() +
           i.vertexCornerOffsets.bytes() + i.vertexCorners.bytes() + i.regularRows.bytes() +
           i.gaussCoeff.bytes() + i.irregularRows.bytes() + i.irregularOffsets.bytes() +
           i.coords.bytes() + i.coordsRef.bytes() + i.faceArea.bytes() + i.faceVolume.bytes() +
           i.faceEBend.bytes() + i.faceMeanCurv.bytes() + i.faceNormal.bytes() +
           i.faceERegular.bytes() + i.faceDeformCase.bytes() + i.slotForceBend.bytes() +
           i.slotForceArea.bytes() + i.slotForceVolume.bytes() + i.cornerForceRegular.bytes() +
           i.vertexForceBend.bytes() + i.vertexForceArea.bytes() + i.vertexForceVolume.bytes() +
           i.vertexForceRegular.bytes();
}

} // namespace slimed
