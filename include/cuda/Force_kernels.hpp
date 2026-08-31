/**
 * @file Force_kernels.hpp
 * @brief The force evaluation as per-index kernel bodies.
 *
 * Each function here does the work for one face or one vertex and takes its
 * index as an argument. That is the only shape that serves both drivers: the
 * host backend runs them in a loop, and a CUDA translation unit wraps each in
 * a __global__ that computes the index from blockIdx and threadIdx. There is
 * one body per stage, not a host copy and a device copy to keep in sync.
 *
 * Every pointer comes from DeviceMeshLayout and PatchRowsFlat, so the same
 * arguments describe host memory or device memory unchanged; only the driver
 * knows which.
 *
 * Force accumulation is split into a scatter and a gather, for the reason
 * DeviceMeshLayout documents: each face writes only its own slots, so nothing
 * needs atomics and the result does not depend on how threads are scheduled.
 */
#pragma once

#include "cuda/Device_mesh_layout.hpp"
#include "energy_force/Patch_kernel.hpp"

namespace slimed
{

/**
 * @brief Everything one force evaluation reads and writes.
 *
 * Plain pointers and scalars so it can be passed to a kernel by value. The
 * layout and row pointers are fixed for a run; only coords, and the scalars in
 * PatchParams that depend on the current area and volume, change per call.
 */
struct ForceKernelArgs
{
    // --- topology, uploaded once ------------------------------------------
    const FacePatchDescriptor *descriptors = nullptr;
    const int *oneRingIndices = nullptr;
    const int *faceValence = nullptr;
    const double *faceSpontCurvature = nullptr;
    const unsigned char *faceIsGhost = nullptr;
    const unsigned char *faceIsBoundary = nullptr;
    const int *faceCorners = nullptr;
    const int *vertexSlotOffsets = nullptr;
    const int *vertexSlots = nullptr;
    const int *vertexCornerOffsets = nullptr;
    const int *vertexCorners = nullptr;
    int nFaces = 0;
    int nVertices = 0;

    // --- shape functions, uploaded once -----------------------------------
    const double *regularRows = nullptr;
    const double *gaussCoeff = nullptr;
    int nSamples = 0;
    /// Irregular blocks, and the offset of each (valence, depth, child) slot
    /// into them. Indexed exactly as PatchRowsFlat does.
    const double *irregularRows = nullptr;
    const std::size_t *irregularOffsets = nullptr;
    int irregularDepth = 0;

    // --- per-evaluation input ---------------------------------------------
    const double *coords = nullptr;    ///< nVertices * 3
    const double *coordsRef = nullptr; ///< nVertices * 3, for regularization
    PatchParams patchParams;
    double kCurv = 0.0;
    double gamaShape = 0.0;
    double gamaArea = 0.0;
    bool usingRpi = true;

    // --- output ------------------------------------------------------------
    double *faceArea = nullptr;      ///< nFaces
    double *faceVolume = nullptr;    ///< nFaces
    double *faceEBend = nullptr;     ///< nFaces
    double *faceMeanCurv = nullptr;  ///< nFaces
    double *faceNormal = nullptr;    ///< nFaces * 3
    double *faceERegular = nullptr;  ///< nFaces
    /// Which branch the regularization took per face: 0 none, 1 area, 2 shape.
    /// Reduced afterwards into param.deformationCount.
    unsigned char *faceDeformCase = nullptr;

    /// Scratch, one 3-vector per (face, control point), laid out like
    /// oneRingIndices. Three of them, for the bending, area and volume terms.
    double *slotForceBend = nullptr;
    double *slotForceArea = nullptr;
    double *slotForceVolume = nullptr;
    /// Scratch, one 3-vector per (face, corner), for the regularization force.
    double *cornerForceRegular = nullptr;

    double *vertexForceBend = nullptr;    ///< nVertices * 3
    double *vertexForceArea = nullptr;    ///< nVertices * 3
    double *vertexForceVolume = nullptr;  ///< nVertices * 3
    double *vertexForceRegular = nullptr; ///< nVertices * 3
};

/// The shape-function block set for one face, or nullptr if it has none.
SLIMED_HD inline const double *rows_for_face(const ForceKernelArgs &args, int face, int child);

/// Face `face`'s limit-surface area and signed volume. Ghost faces get zero.
SLIMED_HD inline void area_volume_for_face(const ForceKernelArgs &args, int face);

/// Face `face`'s bending energy, mean curvature, normal, and its own slots of
/// the bending, area and volume forces.
SLIMED_HD inline void patch_force_for_face(const ForceKernelArgs &args, int face);

/// Face `face`'s regularization energy and its three corner force slots.
SLIMED_HD inline void regularization_for_face(const ForceKernelArgs &args, int face);

/// Sum the bending, area and volume slots belonging to vertex `vertex`.
SLIMED_HD inline void gather_force_for_vertex(const ForceKernelArgs &args, int vertex);

/// Sum the regularization corners belonging to vertex `vertex`.
SLIMED_HD inline void gather_regularization_for_vertex(const ForceKernelArgs &args, int vertex);

} // namespace slimed

// See the note at the end of Patch_kernel.hpp: the bodies are inline and live
// in a .inl so nvcc can compile the same source for the device.
#include "cuda/Force_kernels.inl"
