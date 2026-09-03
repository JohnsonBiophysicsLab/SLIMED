#include "test_device_mesh_layout.hpp"

/**
 * The CUDA port does not rewrite the physics -- Patch_kernel.inl is the same
 * body on both sides, already pinned against the original. What it does
 * rewrite is the shape of the data: one-ring lists become a CSR array,
 * irregular children become slot arithmetic, and the per-thread force buffers
 * become a scatter into per-face slots followed by a per-vertex gather.
 *
 * That flattening is where a port like this goes wrong, and it is testable
 * without a GPU, because HostForceBackend runs exactly the kernels the device
 * will run over exactly the arrays the device will hold. These tests drive it
 * beside Mesh::Compute_Energy_And_Force() on the same mesh and demand the same
 * forces and energies out.
 */

namespace
{
/// A curved sheet, so bending energy and curvature force are non-trivial.
double bowl_height(double x, double y)
{
    return 0.06 * (x * x + y * y);
}

struct MeshFixture
{
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;
};

/**
 * @brief A triangular grid `nx` by `ny` cells wide, lifted into a bowl.
 *
 * The canonical single-patch fixtures have one interior face each, which
 * exercises the kernels but barely exercises the index maps. A grid gives many
 * interior faces sharing vertices, which is what the scatter and gather are
 * actually for: a vertex here belongs to six faces, and every one of them must
 * land in its gather list exactly once.
 */
MeshFixture build_grid(int nx, int ny)
{
    MeshFixture fixture;
    const double dy = std::sqrt(3.0) / 2.0;
    for (int j = 0; j <= ny; j++)
    {
        for (int i = 0; i <= nx; i++)
        {
            const double x = i + 0.5 * j;
            const double y = j * dy;
            fixture.vertices.push_back({x, y, bowl_height(x, y)});
        }
    }
    const auto index = [nx](int i, int j) { return j * (nx + 1) + i; };
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            // Consistent winding, so the normals do not alternate.
            fixture.faces.push_back({index(i, j), index(i + 1, j), index(i, j + 1)});
            fixture.faces.push_back({index(i + 1, j), index(i + 1, j + 1), index(i, j + 1)});
        }
    }
    return fixture;
}

MeshFixture build_canonical(int valence)
{
    const CanonicalPatch patch = build_canonical_patch(valence);
    MeshFixture fixture;
    for (int v = 0; v < patch.nVertices; v++)
    {
        const double x = patch.layout[v][0];
        const double y = patch.layout[v][1];
        fixture.vertices.push_back({x, y, bowl_height(x, y)});
    }
    for (const std::array<int, 3> &face : patch.faces)
    {
        fixture.faces.push_back({face[0], face[1], face[2]});
    }
    return fixture;
}

/// Everything Compute_Energy_And_Force() writes that the backend also writes.
struct ForceSnapshot
{
    std::vector<double> faceArea, faceVolume, faceEBend, faceMeanCurv, faceERegular, faceNormal;
    std::vector<double> forceBend, forceArea, forceVolume, forceRegular;
    double totalArea = 0.0;
    double totalVolume = 0.0;
    int shapeDeformCount = 0;
    int areaDeformCount = 0;
};

ForceSnapshot snapshot(const Mesh &mesh)
{
    ForceSnapshot out;
    for (const Face &face : mesh.faces)
    {
        out.faceArea.push_back(face.elementArea);
        out.faceVolume.push_back(face.elementVolume);
        out.faceEBend.push_back(face.energy.energyCurvature);
        out.faceMeanCurv.push_back(face.meanCurvature);
        out.faceERegular.push_back(face.energy.energyRegularization);
        for (int axis = 0; axis < 3; axis++)
        {
            out.faceNormal.push_back(face.normVector.get(axis, 0));
        }
    }
    for (const Vertex &vertex : mesh.vertices)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            out.forceBend.push_back(vertex.force.forceCurvature.get(axis, 0));
            out.forceArea.push_back(vertex.force.forceArea.get(axis, 0));
            out.forceVolume.push_back(vertex.force.forceVolume.get(axis, 0));
            out.forceRegular.push_back(vertex.force.forceRegularization.get(axis, 0));
        }
    }
    out.totalArea = mesh.param.area;
    out.totalVolume = mesh.param.vol;
    out.shapeDeformCount = mesh.param.deformationCount.shapeDeformCount;
    out.areaDeformCount = mesh.param.deformationCount.areaDeformCount;
    return out;
}

void expect_vectors_close(const std::vector<double> &reference, const std::vector<double> &actual,
                          const char *what)
{
    ASSERT_EQ(reference.size(), actual.size()) << what;
    double worst = 0.0;
    std::size_t worstAt = 0;
    for (std::size_t i = 0; i < reference.size(); i++)
    {
        const double scale = std::max(std::abs(reference[i]), 1.0);
        const double error = std::abs(reference[i] - actual[i]) / scale;
        if (error > worst)
        {
            worst = error;
            worstAt = i;
        }
    }
    EXPECT_LT(worst, 1e-13) << what << ": worst relative disagreement at index " << worstAt
                            << ", reference " << reference[worstAt] << " vs " << actual[worstAt];
}

/// Set the mesh up the way a real run does before its first force evaluation.
void prime(Mesh &mesh)
{
    mesh.calculate_element_area_volume();
    double area = 0.0;
    double volume = 0.0;
    mesh.sum_membrane_area_and_volume(area, volume);
    // Both constraints held away from equilibrium, so the area and volume
    // force terms are non-zero and actually get compared.
    mesh.param.area0 = 0.9 * area;
    mesh.param.vol0 = (volume == 0.0) ? 1.0 : 0.8 * volume;
    mesh.param.uSurf = 250.0;
    mesh.param.uVol = 30.0;
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
}

/// What runs the flattened pipeline: the host backend, or the device one.
enum class Backend
{
    Host,
    Cuda,
};

/// Run both pipelines over one fixture and compare every result they share.
void compare_pipelines(const MeshFixture &fixture, Backend backend = Backend::Host)
{
    Param param;
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    Mesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces));
    prime(mesh);

    ASSERT_NO_THROW(mesh.Compute_Energy_And_Force());
    const ForceSnapshot reference = snapshot(mesh);

    slimed::DeviceMeshLayout layout;
    ASSERT_NO_THROW(layout.build(mesh));
    PatchRowsFlat rows;
    ASSERT_NO_THROW(rows.build(mesh.param.shapeFunctions, mesh.irregularRows,
                               mesh.param.gaussQuadratureCoeff));

    if (backend == Backend::Host)
    {
        slimed::HostForceBackend hostBackend;
        ASSERT_NO_THROW(hostBackend.evaluate(mesh, layout, rows));
    }
    else
    {
        slimed::CudaForceBackend cudaBackend;
        ASSERT_NO_THROW(cudaBackend.upload_topology(layout, rows));
        ASSERT_NO_THROW(cudaBackend.evaluate(mesh, layout, rows));
    }
    const ForceSnapshot actual = snapshot(mesh);

    expect_vectors_close(reference.faceArea, actual.faceArea, "element area");
    expect_vectors_close(reference.faceVolume, actual.faceVolume, "element volume");
    expect_vectors_close(reference.faceEBend, actual.faceEBend, "bending energy");
    expect_vectors_close(reference.faceMeanCurv, actual.faceMeanCurv, "mean curvature");
    expect_vectors_close(reference.faceERegular, actual.faceERegular, "regularization energy");
    expect_vectors_close(reference.faceNormal, actual.faceNormal, "face normal");
    expect_vectors_close(reference.forceBend, actual.forceBend, "curvature force");
    expect_vectors_close(reference.forceArea, actual.forceArea, "area force");
    expect_vectors_close(reference.forceVolume, actual.forceVolume, "volume force");
    expect_vectors_close(reference.forceRegular, actual.forceRegular, "regularization force");
    EXPECT_NEAR(reference.totalArea, actual.totalArea, 1e-13 * std::abs(reference.totalArea));
    EXPECT_NEAR(reference.totalVolume, actual.totalVolume,
                1e-13 * std::max(std::abs(reference.totalVolume), 1.0));
    EXPECT_EQ(reference.shapeDeformCount, actual.shapeDeformCount);
    EXPECT_EQ(reference.areaDeformCount, actual.areaDeformCount);

    // A mesh where every force came out zero would pass all of the above.
    double largestForce = 0.0;
    for (double component : reference.forceBend)
    {
        largestForce = std::max(largestForce, std::abs(component));
    }
    EXPECT_GT(largestForce, 0.0) << "the fixture produced no curvature force to compare";
}
} // namespace

/**
 * A grid, where vertices are shared by six faces each. This is the case the
 * scatter and gather exist for: the production path sums per-thread buffers,
 * the flattened path sums a per-vertex slot list, and they have to agree.
 */
TEST(DeviceMeshLayoutTest, FlattenedPipelineMatchesProductionOnAGrid)
{
    compare_pipelines(build_grid(6, 6));
}

/**
 * The canonical patches, one per valence. Valence 6 is regular; 4, 5, 7 and 8
 * each exercise the irregular child-slot arithmetic, which is the part of the
 * flattening with real arithmetic in it -- a wrong slot returns plausible rows
 * for the wrong depth rather than crashing.
 */
TEST(DeviceMeshLayoutTest, FlattenedPipelineMatchesProductionOnIrregularPatches)
{
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        SCOPED_TRACE("valence " + std::to_string(valence));
        compare_pipelines(build_canonical(valence));
    }
}

/**
 * The gather map is the transpose of the one-ring lists, and a transpose is
 * easy to get subtly wrong in a way the force comparison would only catch by
 * luck. Check it structurally: every slot appears in exactly one vertex's
 * list, and it appears in the list of the vertex it actually belongs to.
 */
TEST(DeviceMeshLayoutTest, GatherMapIsTheTransposeOfTheOneRings)
{
    const MeshFixture fixture = build_grid(5, 5);
    Param param;
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    Mesh mesh(param);
    ASSERT_NO_THROW(mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces));

    slimed::DeviceMeshLayout layout;
    ASSERT_NO_THROW(layout.build(mesh));

    std::vector<int> timesListed(layout.nSlots(), 0);
    for (int vertex = 0; vertex < layout.nVertices(); vertex++)
    {
        for (int i = layout.vertexSlotOffsets()[vertex];
             i < layout.vertexSlotOffsets()[vertex + 1]; i++)
        {
            const int slot = layout.vertexSlots()[i];
            ASSERT_GE(slot, 0);
            ASSERT_LT(slot, layout.nSlots());
            EXPECT_EQ(layout.oneRingIndices()[slot], vertex)
                << "slot " << slot << " is listed under vertex " << vertex
                << " but belongs to vertex " << layout.oneRingIndices()[slot];
            timesListed[slot]++;
        }
    }
    for (int slot = 0; slot < layout.nSlots(); slot++)
    {
        EXPECT_EQ(timesListed[slot], 1) << "slot " << slot << " is gathered " << timesListed[slot]
                                        << " times; a force contribution is dropped or "
                                           "double-counted";
    }

    // The same for the regularization corners.
    std::vector<int> cornerListed(mesh.faces.size() * 3, 0);
    for (int vertex = 0; vertex < layout.nVertices(); vertex++)
    {
        for (int i = layout.vertexCornerOffsets()[vertex];
             i < layout.vertexCornerOffsets()[vertex + 1]; i++)
        {
            const int corner = layout.vertexCorners()[i];
            EXPECT_EQ(layout.faceCorners()[corner], vertex);
            cornerListed[corner]++;
        }
    }
    for (std::size_t corner = 0; corner < cornerListed.size(); corner++)
    {
        EXPECT_EQ(cornerListed[corner], 1) << "corner " << corner;
    }

    // And the descriptors have to agree with the mesh they came from.
    ASSERT_EQ(layout.nFaces(), static_cast<int>(mesh.faces.size()));
    for (int face = 0; face < layout.nFaces(); face++)
    {
        const slimed::FacePatchDescriptor &descriptor = layout.descriptors()[face];
        const int width = static_cast<int>(mesh.faces[face].oneRingVertices.size());
        if (descriptor.kind == slimed::PatchKind::None)
        {
            EXPECT_EQ(descriptor.nChildren, 0) << "face " << face;
            continue;
        }
        EXPECT_EQ(descriptor.nControlPoints, width) << "face " << face;
        for (int j = 0; j < width; j++)
        {
            EXPECT_EQ(layout.oneRingIndices()[descriptor.oneRingOffset + j],
                      mesh.faces[face].oneRingVertices[j])
                << "face " << face << " slot " << j;
        }
    }
}

/**
 * The GPU against the CPU, on the same meshes and to the same tolerance.
 *
 * This is the test to run first on a machine with a device -- it is the only
 * thing in the suite that executes a CUDA kernel, and so the only thing that
 * can tell you a GPU build is trustworthy:
 *
 *     ctest --test-dir build -R CudaForceBackend
 *
 * It skips rather than passes when CUDA was not compiled in or no device is
 * visible, so a green run on a laptop means nothing was checked. Read the skip
 * message.
 *
 * A failure here is a failure of the plumbing -- transfers, launch bounds,
 * device memory -- not of the physics: the kernel bodies are the same inline
 * functions the host backend runs, and DeviceMeshLayoutTest has already
 * checked those against the production path.
 */
TEST(CudaForceBackendTest, MatchesTheProductionForceEvaluation)
{
    if (!slimed::cuda_backend_available())
    {
        GTEST_SKIP() << "built without CUDA; reconfigure with -DSLIMED_ENABLE_CUDA=ON to "
                        "exercise the GPU backend";
    }
    try
    {
        slimed::CudaForceBackend probe;
        SUCCEED() << "using device: " << probe.device_description();
    }
    catch (const std::exception &error)
    {
        GTEST_SKIP() << "no usable CUDA device: " << error.what();
    }

    compare_pipelines(build_grid(6, 6), Backend::Cuda);
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        SCOPED_TRACE("valence " + std::to_string(valence));
        compare_pipelines(build_canonical(valence), Backend::Cuda);
    }
}
