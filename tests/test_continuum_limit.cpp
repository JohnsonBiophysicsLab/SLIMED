#include "test_continuum_limit.hpp"

/**
 * Where does the discrete bending energy stop being the continuum one?
 *
 * The Helfrich energy of a nearly flat surface, in the Monge gauge, is
 *
 *     E = (kc/2) integral (laplacian h)^2 dA
 *       = (A/2) sum_q kc q^4 |h_q|^2
 *
 * so a single plane wave h(r) = eps cos(q.r) on a periodic box of projected
 * area A costs exactly
 *
 *     E_continuum(q) = (A/4) kc q^4 eps^2.
 *
 * SLIMED does not evaluate that. It evaluates 2 kc integral H^2 dA by Gauss
 * quadrature over the Loop limit surface, which is a quartic box spline
 * through the control net. The two agree as q -> 0 and part company as the
 * wavelength approaches the mesh spacing, because a quartic box spline cannot
 * carry curvature it has no knots for.
 *
 * This test measures where. It puts an exact plane wave on the control net,
 * one allowed wavevector at a time, reads the bending energy the production
 * code returns, and forms the stiffness K = 4 E / eps^2.
 *
 * One step is needed before that can be compared with anything. The energy is
 * quadratic in the *control* amplitude, but the continuum law is written for
 * the amplitude of the surface, and the two differ by the Fourier symbol of
 * the Loop limit mask,
 *
 *     m(q) = 1/2 + (1/6) [ cos(q.a1) + cos(q.a2) + cos(q.a3) ],
 *
 * which is 1 at long wavelength and about 1/4 at the zone corner. Sampling the
 * limit surface at the vertices returns m(q) times the control amplitude, so
 * the quantity to compare with the continuum law is
 *
 *     K_S(q) / (A kc q^4),    K_S = K / m(q)^2,
 *
 * and it is 1 wherever the discretisation is faithful. That it comes out at
 * 1.000 rather than at some other constant is itself a check on the mask:
 * nothing here was fitted.
 *
 * Nothing in this file samples, integrates in time, or assumes anything about
 * the dynamics. It is a static property of the energy functional and the mesh,
 * which is what makes it the right place to settle how far up in q a measured
 * fluctuation spectrum may be fitted against kT/(A kc q^4).
 *
 * There is no universal constant for that window. The published cutoffs for
 * lipid-bilayer simulations -- Brandt, Braun, Sachs, Nagle and Edholm,
 * Biophys J 100, 2104 (2011), who put it near q = 0.7 /nm for DMPC -- come
 * from molecular structure: the density structure factor, thickness and
 * protrusion modes. None of those exist in a continuum subdivision-surface
 * model, so the number does not transfer. Shiba and Noguchi, Phys Rev E 84,
 * 031926 (2011), show how strongly a fitted bending rigidity depends on where
 * the cutoff is put, which is the reason to derive it rather than pick it.
 * What limits SLIMED is the mesh, and the mesh can simply be asked.
 *
 * See docs/fluctuation_spectrum.md.
 */
namespace
{

struct Wave
{
    double q = 0.0;
    double qdx = 0.0;     ///< q * dFaceX, the natural dimensionless abscissa
    double control = 0.0; ///< K(q) / (A kc q^4), referred to the control net
    double surface = 0.0; ///< the same, referred to the limit surface
};

/// Fourier symbol of the valence-6 Loop limit mask on this lattice.
double limit_mask_symbol(double qx, double qy, double dx, double dy)
{
    // The six neighbours are +/- a1, +/- a2, +/- a3, with a1 = (dx, 0),
    // a2 = (dx/2, dy) and a3 = a2 - a1 = (-dx/2, dy).
    return 0.5 + (std::cos(qx * dx) +
                  std::cos(qx * dx / 2.0 + qy * dy) +
                  std::cos(-qx * dx / 2.0 + qy * dy)) / 6.0;
}

/// Total bending energy of the periodic tile, for h(r) = amplitude*cos(q.r).
double bending_energy_of_plane_wave(Mesh &mesh, double qx, double qy, double amplitude)
{
    for (Vertex &vertex : mesh.vertices)
    {
        const double x = vertex.coord(0, 0);
        const double y = vertex.coord(1, 0);
        vertex.coord.set(2, 0, amplitude * std::cos(qx * x + qy * y));
    }

    mesh.calculate_element_area_volume();
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.Compute_Energy_And_Force();

    // The ghost faces tile the same wave, so summing the physical ones gives
    // the energy of exactly one periodic box.
    double bending = 0.0;
    for (const Face &face : mesh.faces)
    {
        if (!face.isGhost)
        {
            bending += face.energy.energyCurvature;
        }
    }
    return bending;
}

/// A flat periodic box; the default is the shipped resolution, a 14 x 18 tile.
void build_flat_mesh(Param &param, std::unique_ptr<Mesh> &mesh,
                     double side = 100.0, double lFace = 5.0)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Periodic;
    param.sideX = side;
    param.sideY = side;
    param.lFace = lFace;
    param.kCurv = 83.4;
    param.uSurf = 0.0;  // pure bending: no area term in the energy
    param.uVol = 0.0;
    param.spontCurv = 0.0;
    mesh.reset(new Mesh(param));
    mesh->setup_flat();
    mesh->calculate_element_area_volume();
    mesh->sum_membrane_area_and_volume(mesh->param.area0, mesh->param.vol0);
    mesh->param.vol0 = 0.0;
}

} // namespace

/**
 * The plane-wave stiffness has to be quadratic in the amplitude before any of
 * this means anything; otherwise the probe is measuring the Helfrich
 * nonlinearity rather than the discretisation.
 */
TEST(ContinuumLimitTest, PlaneWaveEnergyIsQuadraticInTheAmplitude)
{
    Param param;
    std::unique_ptr<Mesh> mesh;
    build_flat_mesh(param, mesh);

    const double lx = (param.nFaceX - 6) * param.dFaceX;
    const double qx = 2.0 * M_PI / lx;  // the longest wave the box holds

    const double e1 = bending_energy_of_plane_wave(*mesh, qx, 0.0, 0.01);
    const double e2 = bending_energy_of_plane_wave(*mesh, qx, 0.0, 0.02);
    ASSERT_GT(e1, 0.0);
    EXPECT_NEAR(e2 / e1, 4.0, 0.04) << "bending energy is not quadratic in the "
                                       "amplitude; the probe is outside the "
                                       "harmonic regime";
}

/**
 * The measurement, and the assertions that fix the fitting window.
 */
TEST(ContinuumLimitTest, DiscreteBendingStiffnessMatchesTheContinuumAtSmallQdx)
{
    Param param;
    std::unique_ptr<Mesh> mesh;
    build_flat_mesh(param, mesh);

    const int nx = param.nFaceX - 6;
    const int ny = param.nFaceY - 6;
    const double lx = nx * param.dFaceX;
    const double ly = ny * param.dFaceY;
    const double area = lx * ly;
    const double amplitude = 0.01;   // nm, 1/500 of a cell

    // Sweep the allowed wavevectors of the tile along three directions, so any
    // anisotropy of the triangular lattice would show up: x, y, and diagonal.
    const int directions[3][2] = {{1, 0}, {0, 1}, {1, 1}};
    const char *names[3] = {"x", "y", "xy"};
    std::vector<Wave> waves;

    std::printf("\n[ContinuumLimit] %d x %d tile, dFaceX = %.4g nm, "
                "dFaceY = %.4g nm, A = %.5g nm^2, kc = %.4g pN.nm\n",
                nx, ny, param.dFaceX, param.dFaceY, area, param.kCurv);
    std::printf("[ContinuumLimit] %3s %9s %9s %10s %13s %13s\n",
                "dir", "q (1/nm)", "q*dFaceX", "lambda/dx",
                "K_C/(Akcq^4)", "K_S/(Akcq^4)");

    for (int d = 0; d < 3; d++)
    {
        for (int n = 1; n <= 8; n++)
        {
            const double qx = 2.0 * M_PI * (directions[d][0] * n) / lx;
            const double qy = 2.0 * M_PI * (directions[d][1] * n) / ly;
            const double q = std::sqrt(qx * qx + qy * qy);
            const double qdx = q * param.dFaceX;
            if (qdx > 3.3)
            {
                break;  // near the zone boundary the DFT itself is marginal
            }
            const double e = bending_energy_of_plane_wave(*mesh, qx, qy, amplitude);
            const double k = 4.0 * e / (amplitude * amplitude);
            const double continuum = area * param.kCurv * std::pow(q, 4);
            const double m = limit_mask_symbol(qx, qy, param.dFaceX, param.dFaceY);
            const Wave w{q, qdx, k / continuum, k / (continuum * m * m)};
            waves.push_back(w);
            std::printf("[ContinuumLimit] %3s %9.4f %9.3f %10.2f %13.4f %13.4f\n",
                        names[d], w.q, w.qdx, 2.0 * M_PI / w.qdx,
                        w.control, w.surface);
        }
    }
    ASSERT_FALSE(waves.empty());

    // 1. At long wavelength the discrete energy is the continuum one. This is
    //    the convergence statement, and it also says the limit-mask symbol is
    //    the right thing to have divided by: nothing above was fitted.
    double worstLong = 0.0;
    for (const Wave &w : waves)
    {
        if (w.qdx < 0.9)
        {
            worstLong = std::max(worstLong, std::fabs(w.surface - 1.0));
        }
    }
    EXPECT_LT(worstLong, 0.005) << "the discrete bending energy does not "
                                   "converge to (A/4) kc q^4 eps^2 at long "
                                   "wavelength";

    // 2. Out to q*dFaceX = 1.5 -- about four mesh cells per wavelength -- it
    //    is still within 1%. This is where the default fitting window in
    //    analysis/fluctuation_report.py comes from, and the margin says the
    //    window is conservative rather than marginal.
    double worstWindow = 0.0;
    for (const Wave &w : waves)
    {
        if (w.qdx <= 1.5)
        {
            worstWindow = std::max(worstWindow, std::fabs(w.surface - 1.0));
        }
    }
    EXPECT_LT(worstWindow, 0.01)
        << "the q*dFaceX <= 1.5 fitting window is no longer justified on this "
           "mesh; re-derive it from the table above";

    // 3. Still within 3% at q*dFaceX = 2.0, and within 7% at 2.5, so a wider
    //    window is defensible if the extra reach in q is worth a few per cent.
    //    Asserted so those claims stay true.
    double worstAtTwo = 0.0;
    double worstAtTwoAndAHalf = 0.0;
    for (const Wave &w : waves)
    {
        if (w.qdx <= 2.0)
        {
            worstAtTwo = std::max(worstAtTwo, std::fabs(w.surface - 1.0));
        }
        if (w.qdx <= 2.5)
        {
            worstAtTwoAndAHalf =
                std::max(worstAtTwoAndAHalf, std::fabs(w.surface - 1.0));
        }
    }
    EXPECT_LT(worstAtTwo, 0.03);
    EXPECT_LT(worstAtTwoAndAHalf, 0.07);

    // 4. Past that the limit surface is measurably softer, which is the whole
    //    reason a window exists. Assert the sign, so a change that removed the
    //    softening -- or turned it into a stiffening -- would be noticed.
    double atTop = 1.0;
    double qdxAtTop = 0.0;
    for (const Wave &w : waves)
    {
        if (w.qdx > qdxAtTop)
        {
            qdxAtTop = w.qdx;
            atTop = w.surface;
        }
    }
    EXPECT_LT(atTop, 0.95) << "expected the box spline to be softer than the "
                              "continuum near the zone boundary (q*dFaceX = "
                           << qdxAtTop << ")";
}

/**
 * The softening is a function of q times the mesh spacing, not of q alone.
 *
 * That is what makes q*dFaceX the right abscissa for the fitting window. If
 * the deviation depended on q itself, the window would have to be re-derived
 * for every box size rather than once per mesh resolution.
 *
 * The two boxes below are chosen so their allowed wavevectors line up: a
 * 100 nm box at lFace = 5 and an 85 nm box at lFace = 2.5 both tile to
 * Lx = 70 nm, so mode n of the first and mode 2n of the second sit at exactly
 * the same q*dFaceX -- at different physical q, on meshes a factor of two
 * apart.
 */
TEST(ContinuumLimitTest, TheSofteningScalesWithQTimesMeshSpacing)
{
    struct Case
    {
        double side;
        double lFace;
        int mode;
    };
    const Case cases[2] = {{100.0, 5.0, 3}, {85.0, 2.5, 6}};

    double qdx[2] = {0.0, 0.0};
    double surface[2] = {0.0, 0.0};

    for (int c = 0; c < 2; c++)
    {
        Param param;
        std::unique_ptr<Mesh> mesh;
        build_flat_mesh(param, mesh, cases[c].side, cases[c].lFace);

        const int nx = param.nFaceX - 6;
        const int ny = param.nFaceY - 6;
        const double lx = nx * param.dFaceX;
        const double area = lx * (ny * param.dFaceY);
        const double qx = 2.0 * M_PI * cases[c].mode / lx;
        const double amplitude = 0.002 * param.dFaceX;

        const double e = bending_energy_of_plane_wave(*mesh, qx, 0.0, amplitude);
        const double k = 4.0 * e / (amplitude * amplitude);
        const double m = limit_mask_symbol(qx, 0.0, param.dFaceX, param.dFaceY);
        qdx[c] = qx * param.dFaceX;
        surface[c] = k / (area * param.kCurv * std::pow(qx, 4) * m * m);

        std::printf("[ContinuumLimit] lFace = %.3g nm (tile %d x %d, "
                    "Lx = %.4g nm): q = %.4f /nm, q*dFaceX = %.4f, "
                    "K_S/(A kc q^4) = %.4f\n",
                    cases[c].lFace, nx, ny, lx, qx, qdx[c], surface[c]);
    }

    ASSERT_NEAR(qdx[0], qdx[1], 1e-9) << "the two boxes were meant to sample "
                                         "the same q*dFaceX";
    EXPECT_NEAR(surface[0], surface[1], 0.005)
        << "the deviation from the continuum law is not a function of "
           "q*dFaceX alone, so the fitting window cannot be stated in those "
           "terms";
}
