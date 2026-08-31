/**
 * @file Force_kernels.inl
 * @brief Definitions of the per-index force kernel bodies.
 *
 * Included from the foot of its own header, so every caller carries the
 * definitions and a CUDA translation unit compiles the same source for the
 * device. See Patch_kernel.inl, which follows the same arrangement.
 */
#pragma once

#include "cuda/Force_kernels.hpp"
#include "energy_force/Patch_kernel.inl"

namespace slimed
{

SLIMED_HD inline const double *rows_for_face(const ForceKernelArgs &args, int face, int child)
{
    const FacePatchDescriptor &descriptor = args.descriptors[face];
    if (descriptor.kind == PatchKind::Regular)
    {
        return args.regularRows;
    }
    if (descriptor.kind != PatchKind::Irregular || args.irregularOffsets == nullptr)
    {
        return nullptr;
    }
    // Unpack the child counter back into the (depth, child) pair that indexes
    // the table, using the same slot arithmetic PatchRowsFlat uses.
    const int depth = child / kChildrenPerSubdivisionStep;
    const int childInStep = child % kChildrenPerSubdivisionStep;
    const int valenceSlot = args.faceValence[face] - kMinPatchValence;
    const int slot =
        (valenceSlot * args.irregularDepth + depth) * kChildrenPerSubdivisionStep + childInStep;
    return args.irregularRows + args.irregularOffsets[slot];
}

/// Copy one face's control points out of the coordinate array.
SLIMED_HD inline void load_control_points(const ForceKernelArgs &args, int face, int nControlPoints,
                                          double ctrlPts[kMaxControlPoints * 3])
{
    const int offset = args.descriptors[face].oneRingOffset;
    for (int j = 0; j < nControlPoints; j++)
    {
        const int vertex = args.oneRingIndices[offset + j];
        ctrlPts[j * 3 + 0] = args.coords[vertex * 3 + 0];
        ctrlPts[j * 3 + 1] = args.coords[vertex * 3 + 1];
        ctrlPts[j * 3 + 2] = args.coords[vertex * 3 + 2];
    }
}

SLIMED_HD inline void area_volume_for_face(const ForceKernelArgs &args, int face)
{
    double area = 0.0;
    double volume = 0.0;

    const FacePatchDescriptor &descriptor = args.descriptors[face];
    // Ghost faces do not contribute to the membrane's area or volume, and a
    // face with no complete one-ring has no limit surface to integrate. The
    // energy pass below skips a different set -- boundary rather than ghost --
    // which is why the two flags are carried separately.
    if (!args.faceIsGhost[face] && descriptor.kind != PatchKind::None)
    {
        double ctrlPts[kMaxControlPoints * 3];
        load_control_points(args, face, descriptor.nControlPoints, ctrlPts);
        for (int child = 0; child < descriptor.nChildren; child++)
        {
            element_area_volume_pod(rows_for_face(args, face, child), args.gaussCoeff,
                                    args.nSamples, ctrlPts, descriptor.nControlPoints, area,
                                    volume);
        }
    }

    args.faceArea[face] = area;
    args.faceVolume[face] = volume;
}

SLIMED_HD inline void patch_force_for_face(const ForceKernelArgs &args, int face)
{
    const FacePatchDescriptor &descriptor = args.descriptors[face];
    const int offset = descriptor.oneRingOffset;
    const int nControlPoints = descriptor.nControlPoints;

    double eBend = 0.0;
    double meanCurv = 0.0;
    double normVector[3] = {0.0, 0.0, 0.0};

    if (descriptor.kind != PatchKind::None && !args.faceIsBoundary[face])
    {
        double ctrlPts[kMaxControlPoints * 3];
        double fBend[kMaxControlPoints * 3] = {0.0};
        double fArea[kMaxControlPoints * 3] = {0.0};
        double fVolume[kMaxControlPoints * 3] = {0.0};
        load_control_points(args, face, nControlPoints, ctrlPts);

        PatchParams patchParams = args.patchParams;
        patchParams.spontCurv = args.faceSpontCurvature[face];

        for (int child = 0; child < descriptor.nChildren; child++)
        {
            double childEBend = 0.0;
            double childMeanCurv = 0.0;
            double childNormVector[3] = {0.0, 0.0, 0.0};
            element_energy_force_patch_pod(rows_for_face(args, face, child), args.gaussCoeff,
                                           args.nSamples, ctrlPts, nControlPoints, patchParams,
                                           childEBend, childMeanCurv, childNormVector, fBend,
                                           fArea, fVolume);
            eBend += childEBend;
            // A regular face has one child and takes its values directly. An
            // irregular face reports the child nearest the patch centre --
            // depth 0, the middle child -- rather than summing a curvature and
            // a normal, which do not add.
            const bool isReportingChild =
                (descriptor.kind == PatchKind::Regular) || (child == 1);
            if (isReportingChild)
            {
                meanCurv = childMeanCurv;
                normVector[0] = childNormVector[0];
                normVector[1] = childNormVector[1];
                normVector[2] = childNormVector[2];
            }
        }

        // Each face owns these slots outright, so no other thread writes here.
        for (int j = 0; j < nControlPoints * 3; j++)
        {
            args.slotForceBend[offset * 3 + j] = fBend[j];
            args.slotForceArea[offset * 3 + j] = fArea[j];
            args.slotForceVolume[offset * 3 + j] = fVolume[j];
        }
    }
    else
    {
        // A skipped face still owns slots, because the one-ring list is built
        // from the width alone. Writing zeros rather than leaving them is what
        // makes the gather below safe to run over every slot: the scratch is
        // never cleared between evaluations.
        for (int j = 0; j < nControlPoints * 3; j++)
        {
            args.slotForceBend[offset * 3 + j] = 0.0;
            args.slotForceArea[offset * 3 + j] = 0.0;
            args.slotForceVolume[offset * 3 + j] = 0.0;
        }
    }

    args.faceEBend[face] = eBend;
    args.faceMeanCurv[face] = meanCurv;
    args.faceNormal[face * 3 + 0] = normVector[0];
    args.faceNormal[face * 3 + 1] = normVector[1];
    args.faceNormal[face * 3 + 2] = normVector[2];
}

SLIMED_HD inline void regularization_for_face(const ForceKernelArgs &args, int face)
{
    const int *corners = args.faceCorners + face * 3;
    const int iVertex0 = corners[0];
    const int iVertex1 = corners[1];
    const int iVertex2 = corners[2];

    const double *c0 = args.coords + iVertex0 * 3;
    const double *c1 = args.coords + iVertex1 * 3;
    const double *c2 = args.coords + iVertex2 * 3;

    // Vectors representing three sides of this face, and their lengths.
    double vector10[3];
    double vector21[3];
    double vector02[3];
    for (int axis = 0; axis < 3; axis++)
    {
        vector10[axis] = c0[axis] - c1[axis];
        vector21[axis] = c1[axis] - c2[axis];
        vector02[axis] = c2[axis] - c0[axis];
    }
    const double length10 = v3_norm(vector10);
    const double length21 = v3_norm(vector21);
    const double length02 = v3_norm(vector02);

    // Heron's formula on the current triangle.
    double semiperi = (length10 + length21 + length02) / 2.0;
    const double area =
        std::sqrt(semiperi * (semiperi - length10) * (semiperi - length21) * (semiperi - length02));

    // Shape deformation: the variance of the three side lengths, which is zero
    // for an equilateral triangle.
    const double meanSideLength = (length10 + length21 + length02) / 3.0;
    const double d10 = length10 - meanSideLength;
    const double d21 = length21 - meanSideLength;
    const double d02 = length02 - meanSideLength;
    const double gama =
        (d10 * d10 + d21 * d21 + d02 * d02) / (meanSideLength * meanSideLength);

    const double *r0 = args.coordsRef + iVertex0 * 3;
    const double *r1 = args.coordsRef + iVertex1 * 3;
    const double *r2 = args.coordsRef + iVertex2 * 3;
    double refVector10[3];
    double refVector21[3];
    double refVector02[3];
    for (int axis = 0; axis < 3; axis++)
    {
        refVector10[axis] = r0[axis] - r1[axis];
        refVector21[axis] = r1[axis] - r2[axis];
        refVector02[axis] = r2[axis] - r0[axis];
    }
    const double refLength10 = v3_norm(refVector10);
    const double refLength21 = v3_norm(refVector21);
    const double refLength02 = v3_norm(refVector02);

    semiperi = (refLength10 + refLength21 + refLength02) / 2.0;
    const double refArea = std::sqrt(semiperi * (semiperi - refLength10) *
                                     (semiperi - refLength21) * (semiperi - refLength02));

    const bool isDeformShape = (gama > args.gamaShape && args.usingRpi);
    const double areaExcess = (area - refArea) / refArea;
    const bool isDeformArea =
        ((areaExcess < 0.0 ? -areaExcess : areaExcess) >= args.gamaArea && args.usingRpi);

    // Unit vectors along the three sides. Reciprocal-then-multiply, as
    // Matrix::operator/= did.
    v3_scale(vector10, 1.0 / length10, vector10);
    v3_scale(vector21, 1.0 / length21, vector21);
    v3_scale(vector02, 1.0 / length02, vector02);

    // Whichever branch is taken, the energy is a spring on the three sides
    // pulling them towards some target length; only the target differs. Shape
    // deformation wins over area deformation.
    double target10 = refLength10;
    double target21 = refLength21;
    double target02 = refLength02;
    unsigned char deformCase = 0;
    if (isDeformShape)
    {
        // Towards the equilateral triangle of the same current area.
        const double meanSideLengthOld = std::sqrt(4.0 * area / std::sqrt(3.0));
        target10 = meanSideLengthOld;
        target21 = meanSideLengthOld;
        target02 = meanSideLengthOld;
        deformCase = 2;
    }
    else if (isDeformArea)
    {
        // Towards the equilateral triangle of the reference area.
        const double meanSideLengthRef = std::sqrt(4.0 * refArea / std::sqrt(3.0));
        target10 = meanSideLengthRef;
        target21 = meanSideLengthRef;
        target02 = meanSideLengthRef;
        deformCase = 1;
    }

    const double e10 = length10 - target10;
    const double e21 = length21 - target21;
    const double e02 = length02 - target02;
    args.faceERegular[face] = args.kCurv / 2.0 * (e10 * e10 + e21 * e21 + e02 * e02);
    args.faceDeformCase[face] = deformCase;

    double *corner0 = args.cornerForceRegular + (face * 3 + 0) * 3;
    double *corner1 = args.cornerForceRegular + (face * 3 + 1) * 3;
    double *corner2 = args.cornerForceRegular + (face * 3 + 2) * 3;
    for (int axis = 0; axis < 3; axis++)
    {
        corner0[axis] = args.kCurv * ((target10 - length10) * vector10[axis] +
                                      (length02 - target02) * vector02[axis]);
        corner1[axis] = args.kCurv * ((target21 - length21) * vector21[axis] +
                                      (length10 - target10) * vector10[axis]);
        corner2[axis] = args.kCurv * ((target02 - length02) * vector02[axis] +
                                      (length21 - target21) * vector21[axis]);
    }
}

SLIMED_HD inline void gather_force_for_vertex(const ForceKernelArgs &args, int vertex)
{
    double bend[3] = {0.0, 0.0, 0.0};
    double area[3] = {0.0, 0.0, 0.0};
    double volume[3] = {0.0, 0.0, 0.0};

    const int begin = args.vertexSlotOffsets[vertex];
    const int end = args.vertexSlotOffsets[vertex + 1];
    for (int i = begin; i < end; i++)
    {
        const int slot = args.vertexSlots[i] * 3;
        for (int axis = 0; axis < 3; axis++)
        {
            bend[axis] += args.slotForceBend[slot + axis];
            area[axis] += args.slotForceArea[slot + axis];
            volume[axis] += args.slotForceVolume[slot + axis];
        }
    }

    for (int axis = 0; axis < 3; axis++)
    {
        args.vertexForceBend[vertex * 3 + axis] = bend[axis];
        args.vertexForceArea[vertex * 3 + axis] = area[axis];
        args.vertexForceVolume[vertex * 3 + axis] = volume[axis];
    }
}

SLIMED_HD inline void gather_regularization_for_vertex(const ForceKernelArgs &args, int vertex)
{
    double regular[3] = {0.0, 0.0, 0.0};

    const int begin = args.vertexCornerOffsets[vertex];
    const int end = args.vertexCornerOffsets[vertex + 1];
    for (int i = begin; i < end; i++)
    {
        const int corner = args.vertexCorners[i] * 3;
        for (int axis = 0; axis < 3; axis++)
        {
            regular[axis] += args.cornerForceRegular[corner + axis];
        }
    }

    for (int axis = 0; axis < 3; axis++)
    {
        args.vertexForceRegular[vertex * 3 + axis] = regular[axis];
    }
}

} // namespace slimed
