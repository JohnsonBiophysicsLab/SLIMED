#include "cuda/Host_force_backend.hpp"

#include "cuda/Force_kernels.hpp"
#include "energy_force/Patch_rows_flat.hpp"
#include "mesh/Mesh.hpp"

namespace slimed
{

ForceKernelArgs HostForceBackend::prepare(const DeviceMeshLayout &layout,
                                          const PatchRowsFlat &rows)
{
    const std::size_t nFaces = static_cast<std::size_t>(layout.nFaces());
    const std::size_t nVertices = static_cast<std::size_t>(layout.nVertices());
    const std::size_t nSlots = static_cast<std::size_t>(layout.nSlots());

    // resize() keeps whatever a previous evaluation left behind, which is what
    // we want: every one of these is fully written before it is read.
    faceArea_.resize(nFaces);
    faceVolume_.resize(nFaces);
    faceEBend_.resize(nFaces);
    faceMeanCurv_.resize(nFaces);
    faceNormal_.resize(nFaces * 3);
    faceERegular_.resize(nFaces);
    faceDeformCase_.resize(nFaces);
    slotForceBend_.resize(nSlots * 3);
    slotForceArea_.resize(nSlots * 3);
    slotForceVolume_.resize(nSlots * 3);
    cornerForceRegular_.resize(nFaces * 3 * 3);
    vertexForceBend_.resize(nVertices * 3);
    vertexForceArea_.resize(nVertices * 3);
    vertexForceVolume_.resize(nVertices * 3);
    vertexForceRegular_.resize(nVertices * 3);

    ForceKernelArgs args;
    args.descriptors = layout.descriptors();
    args.oneRingIndices = layout.oneRingIndices();
    args.faceValence = layout.faceValence();
    args.faceSpontCurvature = layout.faceSpontCurvature();
    args.faceIsGhost = layout.faceIsGhost();
    args.faceIsBoundary = layout.faceIsBoundary();
    args.faceCorners = layout.faceCorners();
    args.vertexSlotOffsets = layout.vertexSlotOffsets();
    args.vertexSlots = layout.vertexSlots();
    args.vertexCornerOffsets = layout.vertexCornerOffsets();
    args.vertexCorners = layout.vertexCorners();
    args.nFaces = layout.nFaces();
    args.nVertices = layout.nVertices();

    args.regularRows = rows.regular();
    args.gaussCoeff = rows.gaussCoeff();
    args.nSamples = rows.nSamples();
    args.irregularRows = rows.irregular();
    args.irregularOffsets = rows.irregularOffsets();
    args.irregularDepth = rows.depth();

    args.faceArea = faceArea_.data();
    args.faceVolume = faceVolume_.data();
    args.faceEBend = faceEBend_.data();
    args.faceMeanCurv = faceMeanCurv_.data();
    args.faceNormal = faceNormal_.data();
    args.faceERegular = faceERegular_.data();
    args.faceDeformCase = faceDeformCase_.data();
    args.slotForceBend = slotForceBend_.data();
    args.slotForceArea = slotForceArea_.data();
    args.slotForceVolume = slotForceVolume_.data();
    args.cornerForceRegular = cornerForceRegular_.data();
    args.vertexForceBend = vertexForceBend_.data();
    args.vertexForceArea = vertexForceArea_.data();
    args.vertexForceVolume = vertexForceVolume_.data();
    args.vertexForceRegular = vertexForceRegular_.data();
    return args;
}

void HostForceBackend::evaluate(Mesh &mesh, const DeviceMeshLayout &layout,
                                const PatchRowsFlat &rows)
{
    ForceKernelArgs args = prepare(layout, rows);
    layout.gather_coordinates(mesh, coords_);
    layout.gather_reference_coordinates(mesh, coordsRef_);
    args.coords = coords_.data();
    args.coordsRef = coordsRef_.data();

    Param &param = mesh.param;

    // Stage 1: element area and volume, then the totals they sum to. The
    // constraint forces below read those totals, so this cannot be merged
    // into the force pass -- on a device it is a separate kernel and a
    // reduction for the same reason.
#pragma omp parallel for
    for (int face = 0; face < args.nFaces; face++)
    {
        area_volume_for_face(args, face);
    }

    double area = 0.0;
    double volume = 0.0;
    for (int face = 0; face < args.nFaces; face++)
    {
        mesh.faces[face].elementArea = faceArea_[face];
        mesh.faces[face].elementVolume = faceVolume_[face];
        // Summed in face order rather than by an OpenMP reduction, so the
        // total does not depend on the thread count.
        if (!layout.faceIsGhost()[face])
        {
            area += faceArea_[face];
            volume += faceVolume_[face];
        }
    }
    param.area = area;
    param.vol = volume;

    // Stage 2: the patch energies and forces.
    //
    // No reference quantity means no constraint: a flat sheet encloses nothing,
    // so vol0 is 0 and uVol / vol0 would be 0.0/0.0, and a parameter file that
    // sets relaxArea to 0 would make the area term infinite. Both guards match
    // the ones Compute_Energy_And_Force() applies.
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

#pragma omp parallel for
    for (int face = 0; face < args.nFaces; face++)
    {
        patch_force_for_face(args, face);
        regularization_for_face(args, face);
    }

#pragma omp parallel for
    for (int vertex = 0; vertex < args.nVertices; vertex++)
    {
        gather_force_for_vertex(args, vertex);
        gather_regularization_for_vertex(args, vertex);
    }

    // Stage 3: write the results back into the mesh objects.
    int shapeDeformCount = 0;
    int areaDeformCount = 0;
    for (int face = 0; face < args.nFaces; face++)
    {
        Face &meshFace = mesh.faces[face];
        meshFace.energy.energyCurvature = faceEBend_[face];
        meshFace.meanCurvature = faceMeanCurv_[face];
        meshFace.energy.energyRegularization = faceERegular_[face];
        meshFace.normVector.free();
        meshFace.normVector = mat_calloc(3, 1);
        for (int axis = 0; axis < 3; axis++)
        {
            meshFace.normVector.set(axis, 0, faceNormal_[face * 3 + axis]);
        }
        if (faceDeformCase_[face] == 1)
        {
            areaDeformCount++;
        }
        else if (faceDeformCase_[face] == 2)
        {
            shapeDeformCount++;
        }
    }
    param.deformationCount.shapeDeformCount = shapeDeformCount;
    param.deformationCount.areaDeformCount = areaDeformCount;
    param.deformationCount.noDeformCount =
        static_cast<int>(mesh.faces.size()) - shapeDeformCount - areaDeformCount;

#pragma omp parallel for
    for (int vertex = 0; vertex < args.nVertices; vertex++)
    {
        Force &force = mesh.vertices[vertex].force;
        for (int axis = 0; axis < 3; axis++)
        {
            force.forceCurvature.set(axis, 0, vertexForceBend_[vertex * 3 + axis]);
            force.forceArea.set(axis, 0, vertexForceArea_[vertex * 3 + axis]);
            force.forceVolume.set(axis, 0, vertexForceVolume_[vertex * 3 + axis]);
            force.forceRegularization.set(axis, 0, vertexForceRegular_[vertex * 3 + axis]);
        }
    }
}

std::size_t HostForceBackend::memory_bytes() const
{
    return (coords_.size() + coordsRef_.size() + faceArea_.size() + faceVolume_.size() +
            faceEBend_.size() + faceMeanCurv_.size() + faceNormal_.size() + faceERegular_.size() +
            slotForceBend_.size() + slotForceArea_.size() + slotForceVolume_.size() +
            cornerForceRegular_.size() + vertexForceBend_.size() + vertexForceArea_.size() +
            vertexForceVolume_.size() + vertexForceRegular_.size()) *
               sizeof(double) +
           faceDeformCase_.size();
}

} // namespace slimed
