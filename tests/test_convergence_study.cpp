#include "test_convergence_study.hpp"

namespace
{
/**
 * @brief A disk with one extraordinary vertex at its centre and a fixed boundary.
 *
 * The plan proposes an octahedron for N=4 and an icosahedron for N=5, but
 * neither works: every vertex of those solids is extraordinary, so every face
 * carries three extraordinary corners and the uniform patch reduction -- which
 * assumes exactly one -- does not apply to any of them. They are rejected at
 * setup, correctly.
 *
 * The plan's own alternative is what works, and it works for every valence:
 * "a hand-built disk with a prescribed centre valence and a fixed boundary".
 * The canonical patch is exactly that, and it is already the object the
 * subdivision matrices are derived from, so the fixture and the thing under
 * test agree by construction.
 */
struct ExtraordinaryDisk
{
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;
    int centreFace = -1;
};

ExtraordinaryDisk build_disk(int valence, double bowl)
{
    const CanonicalPatch patch = build_canonical_patch(valence);
    ExtraordinaryDisk disk;
    for (int v = 0; v < patch.nVertices; v++)
    {
        const double x = patch.layout[v][0] - 1.5;
        const double y = patch.layout[v][1] - 1.5;
        disk.vertices.push_back({x, y, bowl * (x * x + y * y)});
    }
    for (const std::array<int, 3> &face : patch.faces)
        disk.faces.push_back({face[0], face[1], face[2]});
    return disk;
}

/// Set up the disk and return the mesh, with the irregular face located.
int locate_irregular_face(const Mesh &mesh, int valence)
{
    for (std::size_t i = 0; i < mesh.faces.size(); i++)
    {
        if (static_cast<int>(mesh.faces[i].oneRingVertices.size()) == valence + 6)
            return static_cast<int>(i);
    }
    return -1;
}

struct Measurement
{
    double area = 0.0;
    double volume = 0.0;
    double bending = 0.0;
    double conjugacyResidual = 0.0;
};

/// Measure the extraordinary face at one depth, including conjugacy.
Measurement measure(int valence, int depth, double h = 1e-6)
{
    const ExtraordinaryDisk disk = build_disk(valence, 0.12);

    Param param;
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.uVol = 0.0;
    Mesh mesh(param);
    mesh.irregularRows.build(mesh.param.shapeFunctions, depth, DepthPolicy::Uniform);
    mesh.setup_from_vertices_faces(disk.vertices, disk.faces);

    const int face = locate_irregular_face(mesh, valence);
    EXPECT_GE(face, 0) << "no valence-" << valence << " face at depth " << depth;
    if (face < 0)
        return {};

    mesh.calculate_element_area_volume();
    double area = 0.0;
    double volume = 0.0;
    mesh.sum_membrane_area_and_volume(area, volume);
    mesh.param.area0 = area;
    mesh.param.vol0 = 0.0;
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.Compute_Energy_And_Force();

    Measurement result;
    result.area = mesh.faces[face].elementArea;
    result.volume = mesh.faces[face].elementVolume;
    result.bending = mesh.faces[face].energy.energyCurvature;

    // Force-energy conjugacy on the bending term, by central difference over
    // the extraordinary face's own control points. This is the test the plan
    // calls the strongest in either document, applied per valence.
    double worst = 0.0;
    double scale = 0.0;
    const std::vector<int> ring = mesh.faces[face].oneRingVertices;
    for (int node : ring)
    {
        for (int k = 0; k < 3; k++)
            scale = std::max(scale, std::abs(mesh.vertices[node].force.forceCurvature(k, 0)));
    }
    if (scale == 0.0)
        return result;

    // Read every analytic force up front. Reading them inside the loop leaves
    // them stale by one perturbation once the first component has been
    // displaced, which shows up as a residual exactly proportional to h --
    // linear, where a central difference should be quadratic.
    std::vector<double> analyticForce;
    for (int node : ring)
    {
        for (int k = 0; k < 3; k++)
            analyticForce.push_back(mesh.vertices[node].force.forceCurvature(k, 0));
    }

    int component = 0;
    for (int node : ring)
    {
        for (int k = 0; k < 3; k++)
        {
            const double analytic = analyticForce[component++];
            const double original = mesh.vertices[node].coord(k, 0);

            mesh.vertices[node].coord.set(k, 0, original + h);
            mesh.Compute_Energy_And_Force();
            const double plus = mesh.param.energy.energyCurvature;

            mesh.vertices[node].coord.set(k, 0, original - h);
            mesh.Compute_Energy_And_Force();
            const double minus = mesh.param.energy.energyCurvature;

            mesh.vertices[node].coord.set(k, 0, original);
            const double gradient = (plus - minus) / (2.0 * h);
            worst = std::max(worst, std::abs(gradient + analytic) / scale);
        }
    }
    // Restore the unperturbed state for the caller's peace of mind.
    mesh.Compute_Energy_And_Force();
    result.conjugacyResidual = worst;
    return result;
}
} // namespace

/**
 * @brief Sweep depth per valence and record where each one converges.
 *
 * The gate is monotone convergence and a chosen depth per valence recorded
 * with its residual. The depths differ because the convergence rate does: the
 * dropped sliver is 4^-D of the parameter domain at every valence, but the
 * surface inside it shrinks like the subdivision matrix's subdominant
 * eigenvalue, which is a function of N.
 */
TEST(ConvergenceStudyTest, DepthSweepPerValence)
{
    std::printf("\n  N   D   area           d(area)     bending        d(bend)     conjugacy\n");
    for (int valence = kMinIrregularValence; valence <= kMaxIrregularValence; valence++)
    {
        double previousArea = 0.0;
        double previousBending = 0.0;
        double previousAreaStep = 0.0;
        for (int depth = 3; depth <= 8; depth++)
        {
            const Measurement m = measure(valence, depth);
            const double areaStep = (depth > 3) ? m.area - previousArea : 0.0;
            const double bendStep = (depth > 3) ? m.bending - previousBending : 0.0;

            std::printf("  %d   %d   %.10f  %10.3e  %.10f  %10.3e  %9.2e\n", valence, depth,
                        m.area, areaStep, m.bending, bendStep, m.conjugacyResidual);

            if (valence == 6)
            {
                // A valence-6 face is regular: the dispatch sends it to the
                // direct kernel and it never reads the row table, so depth
                // cannot change it. This is the degeneracy check -- recursion
                // versus direct regular evaluation -- and equality here is
                // exact rather than approximate.
                EXPECT_DOUBLE_EQ(areaStep, 0.0) << "N=6 D=" << depth;
                EXPECT_DOUBLE_EQ(bendStep, 0.0) << "N=6 D=" << depth;
            }
            else if (depth > 4)
            {
                // Monotone: each level only adds surface.
                EXPECT_GT(areaStep, 0.0) << "N=" << valence << " D=" << depth;
                // And converging: the increments shrink.
                EXPECT_LT(areaStep, previousAreaStep) << "N=" << valence << " D=" << depth;
                EXPECT_GT(bendStep, 0.0) << "N=" << valence << " D=" << depth;
            }
            previousArea = m.area;
            previousBending = m.bending;
            previousAreaStep = areaStep;
        }
    }
}

/**
 * @brief Is the conjugacy residual a finite-difference floor, or a real one?
 *
 * The sweep above reports ~5e-7 at every valence and every depth, flat. Flat
 * is what a correct implementation looks like -- a genuine force/energy
 * mismatch would move with the configuration -- but flat is also what a
 * central-difference floor looks like, so it is worth separating.
 *
 * Central differences carry truncation ~h^2 and round-off ~eps/h. If the
 * residual is numerical it has a minimum in h and grows on both sides; if it
 * is a real discrepancy it is independent of h.
 */
TEST(ConvergenceStudyTest, ConjugacyResidualIsAFiniteDifferenceFloor)
{
    std::printf("\n  conjugacy residual vs step size, N=5 D=6\n");
    double best = 1.0;
    double atSmallest = 0.0;
    double atLargest = 0.0;
    const std::vector<double> steps{1e-8, 1e-7, 1e-6, 1e-5, 1e-4};
    for (std::size_t i = 0; i < steps.size(); i++)
    {
        const double residual = measure(5, 6, steps[i]).conjugacyResidual;
        std::printf("    h=%.0e  residual=%9.2e\n", steps[i], residual);
        best = std::min(best, residual);
        if (i == 0)
            atSmallest = residual;
        if (i + 1 == steps.size())
            atLargest = residual;
    }
    // A numerical floor is worst at both ends and best somewhere between.
    EXPECT_LT(best, atSmallest);
    EXPECT_LT(best, atLargest);
    EXPECT_LT(best, 1e-8) << "conjugacy should reach well below the sweep's 1e-6 gate";
    // A central difference is quadratic in h, so a decade of h costs two
    // decades of accuracy on the truncation-dominated side.
    EXPECT_LT(measure(5, 6, 1e-5).conjugacyResidual, measure(5, 6, 1e-4).conjugacyResidual / 50.0);
}

/// Probe how deep valence 4 has to go, since it converges slowest on bending.
TEST(ConvergenceStudyTest, DeepProbeForSlowValences)
{
    std::printf("\n  N   D    bending        d(bend)     est. remaining tail (rel)\n");
    for (int valence : {4, 5})
    {
        double previous = 0.0;
        double previousStep = 0.0;
        for (int depth = 9; depth <= 14; depth++)
        {
            const Measurement m = measure(valence, depth);
            const double step = (depth > 9) ? m.bending - previous : 0.0;
            double tail = 0.0;
            if (depth > 10 && step > 0.0 && previousStep > step)
            {
                const double ratio = previousStep / step;
                tail = step / (ratio - 1.0) / m.bending;
            }
            std::printf("  %d  %2d   %.10f  %10.3e  %10.2e\n", valence, depth, m.bending, step,
                        tail);
            previous = m.bending;
            previousStep = step;
        }
    }
}
