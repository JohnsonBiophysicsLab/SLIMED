#include "test_volume_functional.hpp"

namespace
{
/**
 * @brief A triangulated torus: closed, orientable, every vertex at valence 6.
 *
 * The fixture has to be closed, because signed volume is only defined there,
 * and it has to be all-regular, because the irregular energy/force path does
 * not exist yet (see docs/irregular_patch_valence_4_to_8_plan.md). Those two
 * requirements together force genus 1: a closed triangulated surface with
 * every vertex at valence 6 has Euler characteristic zero, so no sphere-like
 * fixture can avoid extraordinary vertices. Hence a torus rather than the
 * icosphere one would reach for first.
 */
struct TorusMesh
{
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;

    /// Analytic volume of the smooth torus the control mesh approximates.
    static double analytic_volume(double majorRadius, double minorRadius)
    {
        return 2.0 * M_PI * M_PI * majorRadius * minorRadius * minorRadius;
    }
};

TorusMesh build_torus(int nMajor, int nMinor, double majorRadius, double minorRadius)
{
    TorusMesh torus;
    const double twoPi = 2.0 * M_PI;

    for (int i = 0; i < nMajor; i++)
    {
        const double theta = twoPi * i / nMajor;
        for (int j = 0; j < nMinor; j++)
        {
            const double phi = twoPi * j / nMinor;
            const double ring = majorRadius + minorRadius * std::cos(phi);
            torus.vertices.push_back({ring * std::cos(theta),
                                      ring * std::sin(theta),
                                      minorRadius * std::sin(phi)});
        }
    }

    // Wound so that every face normal points out of the tube.
    auto index = [&](int i, int j) { return ((i % nMajor) * nMinor) + (j % nMinor); };
    for (int i = 0; i < nMajor; i++)
    {
        for (int j = 0; j < nMinor; j++)
        {
            torus.faces.push_back({index(i, j), index(i + 1, j), index(i + 1, j + 1)});
            torus.faces.push_back({index(i, j), index(i + 1, j + 1), index(i, j + 1)});
        }
    }
    return torus;
}

/**
 * @brief Assert that setup left the generator's winding alone.
 *
 * Until WP1, set_one_ring_vertices_sorted() decided face winding with
 * `dot(centre, (c7 - c4) x (c8 - c4)) < 0`, comparing against the coordinate
 * origin. That is only meaningful for a surface star-shaped about the origin,
 * and no torus is: faces on the inner equator have outward normals pointing
 * back toward the axis, so the test flipped them. Exactly half of this fixture
 * -- 128 of 256 faces -- used to come out of setup wound backwards, and these
 * tests carried a repair step to undo it.
 *
 * WP1 removed that test, so the repair is gone and this assertion stands in
 * its place: it fails if per-face winding decisions ever come back.
 */
void expect_winding_preserved(const Mesh &mesh, const TorusMesh &torus)
{
    int nRewound = 0;
    for (std::size_t i = 0; i < mesh.faces.size(); i++)
    {
        if (mesh.faces[i].adjacentVertices != torus.faces[i])
            nRewound++;
    }
    EXPECT_EQ(nRewound, 0) << "setup rewound " << nRewound << " of " << mesh.faces.size()
                           << " faces; the generator already emits consistent outward winding";
}

/// Translate every vertex by the same offset.
void translate(Mesh &mesh, double dx, double dy, double dz)
{
    for (Vertex &vertex : mesh.vertices)
    {
        vertex.coord.set(0, 0, vertex.coord(0, 0) + dx);
        vertex.coord.set(1, 0, vertex.coord(1, 0) + dy);
        vertex.coord.set(2, 0, vertex.coord(2, 0) + dz);
    }
}

double total_volume(Mesh &mesh)
{
    mesh.calculate_element_area_volume();
    double area = 0.0;
    double volume = 0.0;
    mesh.sum_membrane_area_and_volume(area, volume);
    return volume;
}
} // namespace

/**
 * @brief The fixture must be what it claims: closed and consistently oriented.
 *
 * Guards the two tests below. If the repair ever stops working -- or if
 * production starts winding a closed mesh correctly and the repair becomes a
 * no-op that silently does the wrong thing -- this fails first, and with a
 * message that says which.
 */
TEST(VolumeFunctionalTest, TorusFixtureIsClosedAndConsistentlyOriented)
{
    Param param;
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.uVol = 0.0;

    const TorusMesh torus = build_torus(16, 8, 10.0, 3.0);
    Mesh mesh(param);
    mesh.setup_from_vertices_faces(torus.vertices, torus.faces);
    expect_winding_preserved(mesh, torus);

    for (const Face &face : mesh.faces)
    {
        ASSERT_EQ(face.oneRingVertices.size(), 12u)
            << "every patch must be regular; the irregular kernel does not exist yet";

        // A torus this coarse could wrap its own 2-ring and name a vertex
        // twice, which would silently corrupt every patch built on it.
        std::vector<int> sorted = face.oneRingVertices;
        std::sort(sorted.begin(), sorted.end());
        ASSERT_EQ(std::adjacent_find(sorted.begin(), sorted.end()), sorted.end())
            << "face " << face.index << " has a repeated one-ring vertex";
    }

    mesh.param.uVol = 1.0;
    EXPECT_NO_THROW(mesh.validate_volume_constraint_topology());

    // A positive volume confirms the faces are wound outward rather than in.
    EXPECT_GT(total_volume(mesh), 0.0);
}

/**
 * @brief Translating a closed surface must not change its volume.
 *
 * V = (1/3) * closed_integral(x . n dA). Under x -> x + D the extra term is
 * (1/3) D . closed_integral(n dA), and that integral vanishes on any closed
 * surface -- but only if the orientation is consistent. So this one assertion
 * covers an orientation flip and a closure regression at once.
 *
 * It does NOT catch the x-only bug. On a closed surface the single-axis form
 * is translation invariant too, by the same identity. Conjugacy below is what
 * catches that one.
 */
TEST(VolumeFunctionalTest, VolumeIsOriginIndependentOnAClosedSurface)
{
    Param param;
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.uVol = 0.0;

    const TorusMesh torus = build_torus(16, 8, 10.0, 3.0);
    Mesh mesh(param);
    mesh.setup_from_vertices_faces(torus.vertices, torus.faces);
    expect_winding_preserved(mesh, torus);

    const double before = total_volume(mesh);
    ASSERT_GT(before, 0.0);

    translate(mesh, 1000.0, -500.0, 250.0);
    const double after = total_volume(mesh);

    EXPECT_LT(std::abs(after - before) / before, 1e-12);
}

/**
 * @brief A single backwards face breaks origin independence.
 *
 * Without this, the test above could pass on a fixture too degenerate to be
 * sensitive to anything.
 *
 * Flipping *half* the faces -- which is what production's origin-referenced
 * winding test actually does to a torus -- turns out not to work as a
 * counter-example: the flipped set is the inner half, and its symmetry about
 * both the z-axis and the equatorial plane leaves closed_integral(n dA) at
 * zero anyway. So the volume stays translation invariant while being wrong,
 * which is worth knowing on its own. One asymmetric flip is what breaks it.
 */
TEST(VolumeFunctionalTest, ASingleBackwardsFaceBreaksOriginIndependence)
{
    Param param;
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.uVol = 0.0;

    TorusMesh torus = build_torus(16, 8, 10.0, 3.0);
    std::swap(torus.faces[0][1], torus.faces[0][2]); // wind one face backwards

    Mesh mesh(param);
    mesh.setup_from_vertices_faces(torus.vertices, torus.faces);
    expect_winding_preserved(mesh, torus);

    // The closure check from step 3 rejects this outright.
    mesh.param.uVol = 1.0;
    EXPECT_THROW(mesh.validate_volume_constraint_topology(), std::runtime_error);

    const double before = total_volume(mesh);
    translate(mesh, 1000.0, -500.0, 250.0);
    const double after = total_volume(mesh);

    EXPECT_GT(std::abs(after - before) / std::abs(before), 1e-3);
}

/**
 * @brief The volume force is minus the gradient of the volume energy.
 *
 * Central finite difference of `energyVolume` against `forceVolume`, per
 * control point per axis. This is the test whose absence let the split
 * survive: the force differentiates the full-divergence volume, so while
 * `vol` was the x-only one the two were not conjugate and nothing said so.
 *
 * It binds on any mesh, closed or not -- it compares the code's force against
 * the code's own energy -- but it is run on the closed fixture so the
 * quantity being differentiated is a real volume.
 */
TEST(VolumeFunctionalTest, VolumeForceIsMinusGradientOfVolumeEnergy)
{
    Param param;
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.uVol = 0.0;

    const TorusMesh torus = build_torus(8, 6, 10.0, 3.0);
    Mesh mesh(param);
    mesh.setup_from_vertices_faces(torus.vertices, torus.faces);
    expect_winding_preserved(mesh, torus);

    double area = 0.0;
    double volume = 0.0;
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(area, volume);

    // Hold the constraint away from equilibrium, or the force is identically
    // zero and the comparison is vacuous.
    mesh.param.area0 = area;
    mesh.param.vol0 = 0.9 * volume;
    mesh.param.uVol = 1.0;

    // energy_force_regularization() reads coordRef, which Run_flat.cpp seeds
    // before its first force evaluation. Without this the run segfaults.
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.Compute_Energy_And_Force();
    ASSERT_GT(mesh.param.energy.energyVolume, 0.0);

    const int nVertices = static_cast<int>(mesh.vertices.size());
    std::vector<double> analytic(nVertices * 3, 0.0);
    double forceScale = 0.0;
    for (int v = 0; v < nVertices; v++)
    {
        for (int k = 0; k < 3; k++)
        {
            const double component = mesh.vertices[v].force.forceVolume(k, 0);
            analytic[v * 3 + k] = component;
            forceScale = std::max(forceScale, std::abs(component));
        }
    }
    ASSERT_GT(forceScale, 0.0);

    const double h = 1e-5;
    for (int v = 0; v < nVertices; v++)
    {
        for (int k = 0; k < 3; k++)
        {
            const double original = mesh.vertices[v].coord(k, 0);

            mesh.vertices[v].coord.set(k, 0, original + h);
            mesh.Compute_Energy_And_Force();
            const double plus = mesh.param.energy.energyVolume;

            mesh.vertices[v].coord.set(k, 0, original - h);
            mesh.Compute_Energy_And_Force();
            const double minus = mesh.param.energy.energyVolume;

            mesh.vertices[v].coord.set(k, 0, original);

            const double gradient = (plus - minus) / (2.0 * h);
            EXPECT_LT(std::abs(gradient + analytic[v * 3 + k]), 1e-6 * forceScale)
                << "vertex " << v << " axis " << k
                << ": dE/dx = " << gradient
                << ", forceVolume = " << analytic[v * 3 + k];
        }
    }
}

namespace
{
/**
 * @brief A flat sheet: the shipped default workload, and the one whose limit
 * surface lies in the z = 0 plane, so its signed volume is exactly zero.
 */
Mesh build_flat_sheet(Param &param)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Free;
    param.sideX = 25.0;
    param.sideY = 25.0;
    param.lFace = 5.0;

    Mesh mesh(param);
    mesh.setup_flat();
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    return mesh;
}
} // namespace

/**
 * @brief A flat sheet has vol0 == 0, and that must not poison the forces.
 *
 * The volume force scales by uVol / vol0. On a flat sheet vol0 is exactly
 * zero, so with the default uvVolumeConstraint = 0.0 that factor was 0.0/0.0,
 * and the NaN reached every vertex through forceVolume and forceTotal. The
 * energy side has guarded this since before the split; the force side did not.
 *
 * This is the shipped default configuration, so the whole default workload
 * halted at step 0 with "could not find an efficient step size".
 */
TEST(VolumeFunctionalTest, FlatSheetWithNoVolumeConstraintHasFiniteForces)
{
    Param param;
    param.uVol = 0.0;
    Mesh mesh = build_flat_sheet(param);

    ASSERT_EQ(mesh.param.vol0, 0.0) << "a flat sheet must enclose nothing";
    ASSERT_EQ(mesh.param.uVol, 0.0);

    mesh.Compute_Energy_And_Force();

    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        const Force &force = mesh.vertices[i].force;
        for (int k = 0; k < 3; k++)
        {
            EXPECT_EQ(force.forceVolume(k, 0), 0.0)
                << "vertex " << i << " axis " << k
                << ": no reference volume means no volume force";
            EXPECT_TRUE(std::isfinite(force.forceTotal(k, 0)))
                << "vertex " << i << " axis " << k
                << ": forceTotal = " << force.forceTotal(k, 0);
        }
    }

    EXPECT_TRUE(std::isfinite(mesh.calculate_mean_force()));
    EXPECT_TRUE(std::isfinite(mesh.get_max_force_magnitude()));
    EXPECT_DOUBLE_EQ(mesh.param.energy.energyVolume, 0.0);
}

/**
 * @brief A non-finite force field must stop the run, not read as converged.
 *
 * `max_magnitude` is a max reduction, and every comparison against a NaN is
 * false, so a poisoned force field used to reduce to 0. Zero is the value that
 * means "no force left to follow", so the line search reported convergence for
 * what was actually a broken force field -- which is what hid the vol0 bug
 * above behind a plausible-looking halt message.
 */
TEST(VolumeFunctionalTest, NonFiniteForceIsReportedRatherThanReadAsZero)
{
    Param param;
    param.uVol = 0.0;
    Mesh mesh = build_flat_sheet(param);
    mesh.Compute_Energy_And_Force();
    ASSERT_TRUE(std::isfinite(mesh.get_max_force_magnitude()));

    mesh.vertices[0].force.forceTotal.set(0, 0, std::nan(""));
    EXPECT_FALSE(std::isfinite(mesh.get_max_force_magnitude()));

    Record record(1);
    Model model(mesh, record);
    model.iteration = 0;
    EXPECT_THROW(model.determine_trial_step_size(), std::runtime_error);
}

/**
 * @brief A zero reference area must not poison the forces either.
 *
 * The sibling of the vol0 case above, and the same asymmetry: the area energy
 * skips itself when area0 == 0, while the area force divided by it
 * unconditionally. area0 comes from `relaxArea` whenever
 * setRelaxAreaToDefault is false, so this is reachable from a parameter file
 * rather than only from a degenerate mesh -- and because uSurf is non-zero by
 * default it produces inf rather than NaN, which spreads just as far.
 *
 * Kept next to the volume case so the two guards are read together; fixing one
 * of a two-instance pattern is what left this one behind in the first place.
 */
TEST(VolumeFunctionalTest, ZeroReferenceAreaHasFiniteForces)
{
    Param param;
    param.uVol = 0.0;
    Mesh mesh = build_flat_sheet(param);

    // build_flat_sheet() sets area0 from the mesh; override it the way a
    // parameter file with relaxArea = 0.0 would.
    mesh.param.area0 = 0.0;
    ASSERT_NE(mesh.param.uSurf, 0.0) << "uSurf must be non-zero for this to bite";

    mesh.Compute_Energy_And_Force();

    for (std::size_t i = 0; i < mesh.vertices.size(); i++)
    {
        const Force &force = mesh.vertices[i].force;
        for (int k = 0; k < 3; k++)
        {
            EXPECT_EQ(force.forceArea(k, 0), 0.0)
                << "vertex " << i << " axis " << k
                << ": no reference area means no area force";
            EXPECT_TRUE(std::isfinite(force.forceTotal(k, 0)))
                << "vertex " << i << " axis " << k
                << ": forceTotal = " << force.forceTotal(k, 0);
        }
    }

    EXPECT_TRUE(std::isfinite(mesh.get_max_force_magnitude()));
    EXPECT_DOUBLE_EQ(mesh.param.energy.energyArea, 0.0);
}
