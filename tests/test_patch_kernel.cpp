#include "test_patch_kernel.hpp"

/**
 * The patch kernel was rewritten from heap-allocated Matrix objects onto flat
 * stack arrays, because malloc/free and cross-library gsl_matrix_get calls
 * were about 70% of the program's runtime. That is a performance change that
 * has every opportunity to become a silent physics change: the two kernels
 * differ in nothing a compiler would catch.
 *
 * So Mesh::element_energy_force_patch_reference() is kept, verbatim, as the
 * oracle, and these tests drive both kernels over the same patches and demand
 * the same numbers out. If the flat kernel is ever edited, this is what says
 * whether the edit changed the model.
 */

namespace
{
/// A curved patch, so bending energy and the curvature force are non-trivial.
double bowl_height(double x, double y)
{
    return 0.15 * (x * x + y * y);
}

/**
 * @brief The canonical patch for `valence`, set up as a Mesh with a real
 * area and volume mismatch.
 *
 * The reference kernel reads the area and volume constraint terms straight off
 * param, so they are given values that differ from their reference quantities:
 * a kernel that dropped the constraint force entirely would still match a
 * fixture where area == area0.
 */
struct PatchFixture
{
    Param param;
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;
};

PatchFixture build_fixture(int valence)
{
    const CanonicalPatch patch = build_canonical_patch(valence);

    PatchFixture fixture;
    fixture.param.VERBOSE_MODE = false;
    fixture.param.boundaryCondition = BoundaryType::Fixed;
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

/// The control-point coordinates of one face, flat and row-major.
std::vector<double> flat_control_points(const Mesh &mesh, const Face &face)
{
    std::vector<double> flat;
    flat.reserve(face.oneRingVertices.size() * 3);
    for (int iVertex : face.oneRingVertices)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            flat.push_back(mesh.vertices[iVertex].coord.get(axis, 0));
        }
    }
    return flat;
}

/// The control points as one (nCtrl, 3) matrix, the shape the area and volume
/// reference wants. Mesh::get_one_ring_vertex_matrix() builds exactly this,
/// but it is protected.
Matrix control_point_matrix(const Mesh &mesh, const Face &face)
{
    const int nCtrl = static_cast<int>(face.oneRingVertices.size());
    Matrix control(nCtrl, 3);
    for (int j = 0; j < nCtrl; j++)
    {
        const Matrix &coord = mesh.vertices[face.oneRingVertices[j]].coord;
        for (int axis = 0; axis < 3; axis++)
        {
            control.set(j, axis, coord.get(axis, 0));
        }
    }
    return control;
}

/// The same control points as the Matrix column vectors the reference wants.
std::vector<Matrix> matrix_control_points(const Mesh &mesh, const Face &face)
{
    std::vector<Matrix> coords;
    coords.reserve(face.oneRingVertices.size());
    for (int iVertex : face.oneRingVertices)
    {
        coords.push_back(mesh.vertices[iVertex].coord);
    }
    return coords;
}

/**
 * @brief Assert that two values agree to rounding.
 *
 * The rewrite preserved the original's operation order, so the two kernels
 * agree far more tightly than this. The tolerance exists because the original
 * contracted its shape functions through gsl_blas_dgemm, and Homebrew's
 * libgslcblas is compiled with fused multiply-add: its rounding cannot be
 * reproduced by an unfused loop, and chasing it would mean forcing std::fma,
 * which is a libm call rather than an instruction on any target without
 * hardware FMA. It is a floor for that, not a licence to drift.
 */
void expect_close(double reference, double actual, const char *what, int face, int index)
{
    const double scale = std::max(std::abs(reference), 1.0);
    EXPECT_NEAR(reference, actual, 1e-12 * scale)
        << what << " on face " << face << ", component " << index;
}
} // namespace

/**
 * Regular faces: a 12-point one-ring evaluated once against the regular shape
 * functions. This is the arm that runs for essentially every face of a flat
 * membrane, and so the arm that the 8.5x speedup was measured on.
 */
TEST(PatchKernelTest, FlatKernelMatchesReferenceOnRegularFaces)
{
    PatchFixture fixture = build_fixture(6);
    Mesh mesh(fixture.param);
    ASSERT_NO_THROW(mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces));

    mesh.calculate_element_area_volume();
    double area = 0.0;
    double volume = 0.0;
    mesh.sum_membrane_area_and_volume(area, volume);
    mesh.param.area = area;
    mesh.param.vol = volume;
    // Both constraints are held away from equilibrium so their force terms are
    // non-zero and actually get compared.
    mesh.param.area0 = 0.9 * area;
    mesh.param.vol0 = (volume == 0.0) ? 1.0 : 0.8 * volume;
    mesh.param.uSurf = 250.0;
    mesh.param.uVol = 30.0;

    PatchRowsFlat flat;
    ASSERT_NO_THROW(flat.build(mesh.param.shapeFunctions, mesh.irregularRows,
                               mesh.param.gaussQuadratureCoeff));

    slimed::PatchParams patchParams;
    patchParams.kCurv = mesh.param.kCurv;
    patchParams.uSurfPerArea = mesh.param.uSurf / mesh.param.area0;
    patchParams.area = mesh.param.area;
    patchParams.area0 = mesh.param.area0;
    patchParams.uVol = mesh.param.uVol / mesh.param.vol0;
    patchParams.vol = mesh.param.vol;
    patchParams.vol0 = mesh.param.vol0;

    int nCompared = 0;
    for (std::size_t i = 0; i < mesh.faces.size(); i++)
    {
        Face &face = mesh.faces[i];
        const int nCtrl = static_cast<int>(face.oneRingVertices.size());
        if (nCtrl != 12)
        {
            continue;
        }
        nCompared++;

        double refEBend = 0.0;
        double refMeanCurv = 0.0;
        Matrix refNorm = mat_calloc(3, 1);
        Matrix refBend = mat_calloc(nCtrl, 3);
        Matrix refArea = mat_calloc(nCtrl, 3);
        Matrix refVol = mat_calloc(nCtrl, 3);
        mesh.element_energy_force_patch_reference(mesh.param.shapeFunctions,
                                                  matrix_control_points(mesh, face), face,
                                                  face.spontCurvature, refMeanCurv, refNorm,
                                                  refEBend, refBend, refArea, refVol);

        const std::vector<double> ctrlPts = flat_control_points(mesh, face);
        double eBend = 0.0;
        double meanCurv = 0.0;
        double norm[3] = {0.0, 0.0, 0.0};
        std::vector<double> fBend(nCtrl * 3, 0.0);
        std::vector<double> fArea(nCtrl * 3, 0.0);
        std::vector<double> fVol(nCtrl * 3, 0.0);
        slimed::PatchParams facePatchParams = patchParams;
        facePatchParams.spontCurv = face.spontCurvature;
        slimed::element_energy_force_patch_pod(flat.regular(), flat.gaussCoeff(), flat.nSamples(),
                                               ctrlPts.data(), nCtrl, facePatchParams, eBend,
                                               meanCurv, norm, fBend.data(), fArea.data(),
                                               fVol.data());

        const int face_index = static_cast<int>(i);
        expect_close(refEBend, eBend, "bending energy", face_index, 0);
        expect_close(refMeanCurv, meanCurv, "mean curvature", face_index, 0);
        for (int axis = 0; axis < 3; axis++)
        {
            expect_close(refNorm.get(axis, 0), norm[axis], "normal", face_index, axis);
        }
        for (int j = 0; j < nCtrl * 3; j++)
        {
            expect_close(refBend.get(j / 3, j % 3), fBend[j], "bending force", face_index, j);
            expect_close(refArea.get(j / 3, j % 3), fArea[j], "area force", face_index, j);
            expect_close(refVol.get(j / 3, j % 3), fVol[j], "volume force", face_index, j);
        }

        refNorm.free();
        refBend.free();
        refArea.free();
        refVol.free();
    }

    // A fixture with no regular face would pass every assertion above without
    // testing anything.
    EXPECT_GT(nCompared, 0) << "the valence-6 fixture produced no 12-point one-ring";
}

/**
 * Irregular faces, valence 4 through 8. Each is tiled by regular children at
 * increasing depth, so this also covers PatchRowsFlat's child indexing: a
 * mis-packed block would show up as the wrong rows rather than as a crash.
 */
TEST(PatchKernelTest, FlatKernelMatchesReferenceOnIrregularFaces)
{
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        if (valence == 6)
        {
            continue; // regular, covered by the test above
        }
        SCOPED_TRACE("valence " + std::to_string(valence));

        PatchFixture fixture = build_fixture(valence);
        Mesh mesh(fixture.param);
        ASSERT_NO_THROW(mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces));

        mesh.calculate_element_area_volume();
        double area = 0.0;
        double volume = 0.0;
        mesh.sum_membrane_area_and_volume(area, volume);
        mesh.param.area = area;
        mesh.param.vol = volume;
        mesh.param.area0 = 0.9 * area;
        mesh.param.vol0 = (volume == 0.0) ? 1.0 : 0.8 * volume;
        mesh.param.uSurf = 250.0;
        mesh.param.uVol = 30.0;

        PatchRowsFlat flat;
        ASSERT_NO_THROW(flat.build(mesh.param.shapeFunctions, mesh.irregularRows,
                                   mesh.param.gaussQuadratureCoeff));

        slimed::PatchParams patchParams;
        patchParams.kCurv = mesh.param.kCurv;
        patchParams.uSurfPerArea = mesh.param.uSurf / mesh.param.area0;
        patchParams.area = mesh.param.area;
        patchParams.area0 = mesh.param.area0;
        patchParams.uVol = mesh.param.uVol / mesh.param.vol0;
        patchParams.vol = mesh.param.vol;
        patchParams.vol0 = mesh.param.vol0;

        int nCompared = 0;
        for (std::size_t i = 0; i < mesh.faces.size(); i++)
        {
            Face &face = mesh.faces[i];
            const int nCtrl = static_cast<int>(face.oneRingVertices.size());
            if (nCtrl != valence + 6)
            {
                continue;
            }
            nCompared++;

            const std::vector<Matrix> matrixCoords = matrix_control_points(mesh, face);
            const std::vector<double> ctrlPts = flat_control_points(mesh, face);

            // Forces accumulate across the children, exactly as
            // Compute_Energy_And_Force() drives them.
            Matrix refBend = mat_calloc(nCtrl, 3);
            Matrix refArea = mat_calloc(nCtrl, 3);
            Matrix refVol = mat_calloc(nCtrl, 3);
            std::vector<double> fBend(nCtrl * 3, 0.0);
            std::vector<double> fArea(nCtrl * 3, 0.0);
            std::vector<double> fVol(nCtrl * 3, 0.0);
            double refEBendTotal = 0.0;
            double eBendTotal = 0.0;

            for (int d = 0; d < mesh.irregularRows.depth_for(valence); d++)
            {
                for (int c = 0; c < kRegularChildrenPerStep; c++)
                {
                    double refEBend = 0.0;
                    double refMeanCurv = 0.0;
                    Matrix refNorm = mat_calloc(3, 1);
                    mesh.element_energy_force_patch_reference(
                        mesh.irregularRows.rows_for_child(valence, d, c), matrixCoords, face,
                        face.spontCurvature, refMeanCurv, refNorm, refEBend, refBend, refArea,
                        refVol);

                    double eBend = 0.0;
                    double meanCurv = 0.0;
                    double norm[3] = {0.0, 0.0, 0.0};
                    slimed::PatchParams facePatchParams = patchParams;
                    facePatchParams.spontCurv = face.spontCurvature;
                    slimed::element_energy_force_patch_pod(
                        flat.child(valence, d, c), flat.gaussCoeff(), flat.nSamples(),
                        ctrlPts.data(), nCtrl, facePatchParams, eBend, meanCurv, norm,
                        fBend.data(), fArea.data(), fVol.data());

                    const int child_index = d * kRegularChildrenPerStep + c;
                    expect_close(refEBend, eBend, "child bending energy", child_index, 0);
                    expect_close(refMeanCurv, meanCurv, "child mean curvature", child_index, 0);
                    for (int axis = 0; axis < 3; axis++)
                    {
                        expect_close(refNorm.get(axis, 0), norm[axis], "child normal", child_index,
                                     axis);
                    }
                    refEBendTotal += refEBend;
                    eBendTotal += eBend;
                    refNorm.free();
                }
            }

            const int face_index = static_cast<int>(i);
            expect_close(refEBendTotal, eBendTotal, "accumulated bending energy", face_index, 0);
            for (int j = 0; j < nCtrl * 3; j++)
            {
                expect_close(refBend.get(j / 3, j % 3), fBend[j], "bending force", face_index, j);
                expect_close(refArea.get(j / 3, j % 3), fArea[j], "area force", face_index, j);
                expect_close(refVol.get(j / 3, j % 3), fVol[j], "volume force", face_index, j);
            }
            // A patch that produced no bending energy would match trivially.
            EXPECT_GT(std::abs(refEBendTotal), 0.0);

            refBend.free();
            refArea.free();
            refVol.free();
        }

        EXPECT_GT(nCompared, 0) << "no face of width " << valence + 6 << " in the fixture";
    }
}

/**
 * The flat row cache has to agree with the Matrix table it repacks, entry for
 * entry. Comparing the kernels would catch a badly mangled block, but a
 * transposed or off-by-one-sample block can still produce plausible forces.
 */
TEST(PatchKernelTest, FlatRowsMatchTheMatrixTable)
{
    PatchFixture fixture = build_fixture(5);
    Mesh mesh(fixture.param);
    ASSERT_NO_THROW(mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces));

    PatchRowsFlat flat;
    ASSERT_NO_THROW(flat.build(mesh.param.shapeFunctions, mesh.irregularRows,
                               mesh.param.gaussQuadratureCoeff));
    ASSERT_EQ(flat.nSamples(), static_cast<int>(mesh.param.shapeFunctions.size()));

    for (int sample = 0; sample < flat.nSamples(); sample++)
    {
        EXPECT_DOUBLE_EQ(flat.gaussCoeff()[sample], mesh.param.gaussQuadratureCoeff.get(sample, 0));

        const Matrix &block = mesh.param.shapeFunctions[sample];
        const double *rows = flat.regular() + sample * slimed::kShapeRows * 12;
        for (int row = 0; row < slimed::kShapeRows; row++)
        {
            for (int col = 0; col < 12; col++)
            {
                EXPECT_DOUBLE_EQ(rows[row * 12 + col], block.get(row, col))
                    << "regular sample " << sample << " (" << row << ", " << col << ")";
            }
        }
    }

    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        const int nCtrl = valence + 6;
        for (int d = 0; d < mesh.irregularRows.depth_for(valence); d++)
        {
            for (int c = 0; c < kRegularChildrenPerStep; c++)
            {
                const std::vector<Matrix> &table = mesh.irregularRows.rows_for_child(valence, d, c);
                const double *rows = flat.child(valence, d, c);
                for (int sample = 0; sample < flat.nSamples(); sample++)
                {
                    const Matrix &block = table[sample];
                    const double *here = rows + sample * slimed::kShapeRows * nCtrl;
                    for (int row = 0; row < slimed::kShapeRows; row++)
                    {
                        for (int col = 0; col < nCtrl; col++)
                        {
                            ASSERT_DOUBLE_EQ(here[row * nCtrl + col], block.get(row, col))
                                << "valence " << valence << " depth " << d << " child " << c
                                << " sample " << sample << " (" << row << ", " << col << ")";
                        }
                    }
                }
            }
        }
    }
}

/**
 * Widths above kMaxControlPoints never reach the kernel, and a mesh without
 * irregular faces never builds the table at all. Both have to say so rather
 * than hand back a pointer into nothing.
 */
TEST(PatchKernelTest, FlatRowsRejectSlotsItDidNotBuild)
{
    PatchFixture fixture = build_fixture(6);
    Mesh mesh(fixture.param);
    ASSERT_NO_THROW(mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces));

    PatchRowsFlat flat;
    ASSERT_NO_THROW(flat.build(mesh.param.shapeFunctions, mesh.irregularRows,
                               mesh.param.gaussQuadratureCoeff));

    EXPECT_THROW(flat.child(kMaxIrregularValence + 1, 0, 0), std::out_of_range);
    EXPECT_THROW(flat.child(kMinIrregularValence, -1, 0), std::out_of_range);
    EXPECT_THROW(flat.child(kMinIrregularValence, 0, kRegularChildrenPerStep), std::out_of_range);

    // The widest patch the mesh can produce still has to fit the kernel's
    // stack buffers; if kMaxIrregularValence ever grows, this is the line that
    // says kMaxControlPoints has to grow with it.
    EXPECT_LE(kMaxIrregularValence + 6, slimed::kMaxControlPoints);
}

/**
 * Area and volume are integrated by a second, cheaper kernel that reads only
 * the first three shape-function rows. They feed the area and volume
 * constraint energies, so a drift here is a drift in the force the patch
 * kernel is checked against above -- worth its own oracle rather than being
 * inferred from the energy.
 *
 * Not exact equality. The code this replaced contracted the shape functions
 * with a gsl_blas_dgemm call, and Homebrew's libgslcblas is compiled with
 * fused multiply-add, so its dgemm rounds each multiply-add differently from
 * any unfused loop -- a property of how the library was built, not of the
 * mathematics. Measured on the 8400-face default mesh the two differ by at
 * most 6.6e-14 on element areas of about 10.8, or 6e-15 relative. The bar
 * below is a hundred times tighter than anything physical and still well
 * clear of that.
 */
TEST(PatchKernelTest, FlatAreaVolumeMatchesReference)
{
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        SCOPED_TRACE("valence " + std::to_string(valence));

        PatchFixture fixture = build_fixture(valence);
        Mesh mesh(fixture.param);
        ASSERT_NO_THROW(mesh.setup_from_vertices_faces(fixture.vertices, fixture.faces));

        PatchRowsFlat flat;
        ASSERT_NO_THROW(flat.build(mesh.param.shapeFunctions, mesh.irregularRows,
                                   mesh.param.gaussQuadratureCoeff));

        int nCompared = 0;
        for (std::size_t i = 0; i < mesh.faces.size(); i++)
        {
            Face &face = mesh.faces[i];
            const int nCtrl = static_cast<int>(face.oneRingVertices.size());
            const bool isRegular = (nCtrl == 12);
            const bool isIrregular = (nCtrl >= kMinIrregularValence + 6 &&
                                      nCtrl <= kMaxIrregularValence + 6);
            if (!isRegular && !isIrregular)
            {
                continue;
            }
            nCompared++;

            const Matrix matrixCoords = control_point_matrix(mesh, face);
            const std::vector<double> ctrlPts = flat_control_points(mesh, face);

            double refArea = 0.0;
            double refVolume = 0.0;
            double area = 0.0;
            double volume = 0.0;

            if (isRegular)
            {
                mesh.enumerate_gauss_quadrature_point_area_volume_reference(
                    mesh.param.shapeFunctions, matrixCoords, refArea, refVolume);
                slimed::element_area_volume_pod(flat.regular(), flat.gaussCoeff(), flat.nSamples(),
                                                ctrlPts.data(), nCtrl, area, volume);
            }
            else
            {
                const int faceValence = nCtrl - 6;
                for (int d = 0; d < mesh.irregularRows.depth_for(faceValence); d++)
                {
                    for (int c = 0; c < kRegularChildrenPerStep; c++)
                    {
                        mesh.enumerate_gauss_quadrature_point_area_volume_reference(
                            mesh.irregularRows.rows_for_child(faceValence, d, c), matrixCoords,
                            refArea, refVolume);
                        slimed::element_area_volume_pod(flat.child(faceValence, d, c),
                                                        flat.gaussCoeff(), flat.nSamples(),
                                                        ctrlPts.data(), nCtrl, area, volume);
                    }
                }
            }

            EXPECT_NEAR(refArea, area, 1e-13 * std::max(std::abs(refArea), 1.0))
                << "area on face " << i;
            EXPECT_NEAR(refVolume, volume, 1e-13 * std::max(std::abs(refVolume), 1.0))
                << "volume on face " << i;
            // A degenerate patch would match at zero without testing anything.
            EXPECT_GT(refArea, 0.0) << "face " << i;
        }

        EXPECT_GT(nCompared, 0) << "no face with a complete one-ring in the fixture";
    }
}
