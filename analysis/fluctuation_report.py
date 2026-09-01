#!/usr/bin/env python3
"""Height-spectrum and relaxation report for one SLIMED dynamics run.

    python analysis/fluctuation_report.py <run-dir> [--burn N] [--stride N]
                                          [--npz out.npz] [--png-prefix p]

<run-dir> holds the ``input.params`` the run used and the ``meshpointinput.csv``
/ ``surfacepointinput.csv`` it wrote.

Three things get checked, in increasing order of how much they assume.

1.  The static spectrum against  <|h_q|^2> = kT / (A (kc q^4 + sigma q^2)).
    This is the q^-4 law.  It only holds in the continuum limit, so the fit is
    restricted to q * dFaceX below --qdx-max (1.5 by default); beyond that the
    quartic box spline that Loop subdivision converges to is no longer a
    faithful stand-in for a smooth surface, and the deviation is a property of
    the discretisation, not an error.

    Two columns qualify each bin.  "T/tau_q" is the run length in units of that
    bin's expected relaxation time; below ~10 the bin is flagged UNDER-SAMPLED,
    because a mode the run never relaxed reads low and drags the slope towards
    zero.  "Euler" is 2/(2-lambda), the factor by which an explicit
    Euler-Maruyama step inflates a mode's stationary variance -- ~1 across the
    fit window, growing towards the zone corner.

2.  An internal fluctuation-dissipation check, which assumes nothing about the
    continuum.  Each mode's variance and its own relaxation rate are two
    measurements of the same Boltzmann distribution; their product has to come
    out at 1 for any correct Brownian update, whatever the stiffness really is.

3.  The relaxation rate against the two candidate dynamics:
        free draining (what the code implements)   Gamma = (D a / kT) kc q^4
        Oseen / Zimm hydrodynamics                 Gamma = kc q^3 / (4 eta)
    Distinguishing q^3 from q^4 is the whole content of the h(t)h(t+dt) test.

The legacy reduction (index-space fft2, quarter fold, index-valued |q|) is
recomputed on the same frames so the two can be compared directly.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import membrane_spectrum as ms  # noqa: E402


def read_scalars(params_path):
    text = open(params_path, encoding="utf-8", errors="replace").read()

    def get(name, default=None, cast=float):
        m = re.search(rf"^\s*{re.escape(name)}\s*=\s*(\S+)", text, re.M)
        if m is None:
            if default is None:
                raise KeyError(name)
            return default
        return cast(m.group(1))

    return dict(
        dt=get("timeStep"),
        diff=get("diffConst"),
        kbt=get("KBT"),
        kc=get("kcMembraneBending"),
        us=get("usMembraneStretching", 0.0),
        interval=get("meshpointOutputInterval", 1, int),
        fdt=get("fdtConsistentSurfaceUpdate", "true", str).lower().startswith("t"),
        seed=get("randomSeed", 0, int),
    )


def legacy_reduction(z_tile, lat):
    """The original notebook's reduction, for side-by-side comparison.

    Reproduces, on whatever frames it is given:
      * ``np.fft.fft2`` of the index array, with no zig-zag phase;
      * the ``arr[0:n//2][0:m//2]`` fold, which slices rows twice and so never
        touches the columns;
      * |q| built from the raw row and column index over a "side length" that
        is a count of vertices rather than a distance.
    """
    arr = np.abs(ms.naive_index_fft2(z_tile))
    trimmed = arr[:, : arr.shape[1] // 2][:, : arr.shape[2] // 2]
    s = np.mean(trimmed ** 2, axis=0).ravel()
    n_row, n_col = trimmed.shape[1], trimmed.shape[2]
    l_mean = (lat.n_vert_x + lat.n_vert_y) / 2.0
    q = np.sqrt(np.add.outer((np.arange(n_row) / l_mean) ** 2,
                             (np.arange(n_col) / l_mean) ** 2)).ravel()
    keep = q > 0
    return q[keep], s[keep]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--burn", type=int, default=None,
                    help="frames to discard (default: 20%% of the trajectory)")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--max-frames", type=int, default=None)
    ap.add_argument("--bins", type=int, default=20)
    ap.add_argument("--qdx-max", type=float, default=1.5,
                    help="fit the continuum law only where q*dFaceX is below this. "
                         "1.5 is where the measured discrete stiffness K_S(q) is "
                         "still within ~10%% of A kc q^4 on the default lattice; "
                         "the K_S column in the output is what to check if the "
                         "lattice changes.")
    ap.add_argument("--source", choices=("control", "surface"), default="control",
                    help="control: meshpoint.csv put through the limit mask with "
                         "periodic wrap (exactly periodic).  surface: the "
                         "matSurface the run wrote (not periodic at the seam).")
    ap.add_argument("--eta", type=float, default=1.0e-3,
                    help="solvent viscosity in pN.us/nm^2 (1e-3 = water)")
    ap.add_argument("--npz", default=None)
    ap.add_argument("--png-prefix", default=None)
    args = ap.parse_args(argv)

    lat = ms.lattice_from_params(os.path.join(args.run_dir, "input.params"))
    par = read_scalars(os.path.join(args.run_dir, "input.params"))
    n_mode = lat.nx * lat.ny
    a_cell = lat.area / n_mode

    print("=" * 78)
    print(os.path.abspath(args.run_dir))
    print(lat.describe())
    print(f"dt = {par['dt']} us, D = {par['diff']} nm^2/us, kT = {par['kbt']} pN.nm, "
          f"kc = {par['kc']} pN.nm, us = {par['us']} pN/nm")
    print(f"output every {par['interval']} step(s), seed {par['seed']}, "
          f"fdtConsistentSurfaceUpdate = {par['fdt']}")

    name = "meshpointinput.csv" if args.source == "control" else "surfacepointinput.csv"
    path = os.path.join(args.run_dir, name)
    n_lines = sum(1 for _ in open(path))
    burn = args.burn if args.burn is not None else n_lines // 5
    print(f"reading {name}: {n_lines} frames, discarding the first {burn}")

    # Read from frame 0 so the warm-up is available for the equilibration
    # check below, then slice.  A run that is still heating has no stationary
    # spectrum to measure, and that has to be checked before anything else is
    # believed.
    z_all = ms.read_tile_heights(path, lat, stride=args.stride, start=0,
                                 max_frames=args.max_frames)
    h_all = ms.limit_surface_from_control(z_all, lat) if args.source == "control" else z_all
    keep = max(0, burn // max(args.stride, 1))
    h = h_all[keep:]
    n_frames = h.shape[0]
    dt_frame = par["dt"] * par["interval"] * args.stride
    print(f"analysing {n_frames} frames, {dt_frame} us apart "
          f"({n_frames * dt_frame:.4g} us total)")

    # ------------------------------------------------------- equilibration
    rms_t = h_all.reshape(h_all.shape[0], -1)
    rms_t = np.sqrt(np.mean((rms_t - rms_t.mean(axis=1, keepdims=True)) ** 2, axis=1))
    t_all = np.arange(h_all.shape[0]) * dt_frame
    nblk = 8
    edges = np.linspace(keep, h_all.shape[0], nblk + 1).astype(int)
    blocks = [float(rms_t[a:b].mean()) for a, b in zip(edges[:-1], edges[1:]) if b > a]
    print("\n--- equilibration -----------------------------------------------------")
    print("rms height over 8 equal blocks of the analysed window (nm):")
    print("  " + "  ".join(f"{v:.3f}" for v in blocks))
    if blocks:
        trend = blocks[-1] / blocks[0]
        verdict = ("stationary" if 0.8 <= trend <= 1.25 else
                   "STILL GROWING - no stationary spectrum to measure" if trend > 1.25
                   else "still relaxing downward")
        print(f"  last/first = {trend:.2f}   -> {verdict}")

    qx, qy, q = ms.q_grid(lat)
    m = ms.limit_mask_symbol(lat)
    H = ms.tile_fft2(h, lat)
    S = ms.power_spectrum(H, lat)

    # -------------------------------------------- per-mode rates, needed first
    # The explicit-Euler variance bias below is a function of the measured
    # decay, so the autocorrelation has to be in hand before the static table.
    max_lag = max(4, min(n_frames // 4, 4000))
    acf = ms.mode_autocorrelation(H, max_lag=max_lag)
    lam_frame = 1.0 - acf[1]
    lam_step = 1.0 - np.abs(1.0 - lam_frame) ** (1.0 / (par["interval"] * args.stride))

    # The code's nominal per-step decay for mode q is
    #
    #   FDT-consistent update:   lam = (D dt / kT) * K_S(q) / N
    #   legacy update:           lam = (D dt / kT) * K_C(q) / (N m)
    #                                = (D dt / kT) * K_S(q) * m / N
    #
    # since K_C = K_S m^2.  Inverting for K_S therefore divides by m in the
    # legacy case and not at all in the fixed one, which is what makes the
    # product below an independent test rather than a tautology.
    m_pow = 0.0 if par["fdt"] else -1.0
    K_S = lam_step * par["kbt"] * n_mode / (par["diff"] * par["dt"]) * m ** m_pow
    #  h <- (1-lam) h + noise  overshoots: the stationary variance of the
    #  discrete chain is 2/(2-lam) times the continuous-time one.
    inflation = 2.0 / (2.0 - lam_step)
    fdt_ratio = S * K_S / par["kbt"] / inflation

    # ---------------------------------------------------------------- static
    qc, sm, se, cnt = ms.radial_average(q, S, n_bins=args.bins)
    fit_mask = qc * lat.d_face_x <= args.qdx_max
    kc_fit, sig_fit = ms.fit_kc_sigma(qc[fit_mask], sm[fit_mask], lat.area, par["kbt"])
    slope = ms.log_slope(qc[fit_mask], sm[fit_mask])

    q_leg, s_leg = legacy_reduction(h, lat)
    slope_leg = ms.log_slope(q_leg, s_leg)

    print("\n--- static spectrum -------------------------------------------------")
    print(f"{'q (1/nm)':>9} {'q*dx':>6} {'<|h_q|^2>':>12} {'kT/(A kc q^4)':>14} "
          f"{'ratio':>7} {'modes':>6} {'rel.err':>8} {'T/tau_q':>8} {'Euler':>6}")
    ideal = par["kbt"] / (lat.area * par["kc"] * qc ** 4)
    # How long each bin needs to equilibrate, from the free-draining rate, so a
    # bin the run cannot have sampled gets labelled rather than believed.
    tau_expect = par["kbt"] / (par["diff"] * a_cell * par["kc"] * qc ** 4)
    t_run = n_frames * dt_frame
    _, infl_bin, _, _ = ms.radial_average(q, inflation, n_bins=args.bins)
    for a, b, c, d, e, t, g in zip(qc, sm, ideal, cnt, se, tau_expect, infl_bin):
        flag = "" if a * lat.d_face_x <= args.qdx_max else "  beyond fit range"
        if t_run / t < 10:
            flag += "  UNDER-SAMPLED"
        print(f"{a:>9.4f} {a*lat.d_face_x:>6.2f} {b:>12.4e} {c:>14.4e} "
              f"{b/c:>7.3f} {d:>6d} {e/b if b else 0:>8.3f} {t_run/t:>8.1f} "
              f"{g:>6.3f}{flag}")

    print(f"\nlog-log slope over q*dx <= {args.qdx_max}: {slope:.3f}   (want -4)")
    print(f"fitted kc    = {kc_fit:8.2f} pN.nm   "
          f"({kc_fit/par['kbt']:.1f} kT; input {par['kc']:.2f} pN.nm, "
          f"{par['kc']/par['kbt']:.1f} kT)")
    print(f"fitted sigma = {sig_fit:8.4f} pN/nm")
    if sig_fit > 0:
        print(f"  tension-bending crossover at q = {math.sqrt(sig_fit/par['kc']):.4f} /nm")
    print(f"legacy reduction on the same frames: slope {slope_leg:.3f}")

    # -------------------------------------------------------- dynamics / FDT

    gamma = np.empty_like(q)
    gamma.fill(np.nan)
    it = np.ndindex(q.shape)
    for idx in it:
        if q[idx] > 0:
            gamma[idx] = ms.relaxation_rate(acf[:, idx[0], idx[1]], dt_frame,
                                            n_frames=n_frames)

    print("\n--- internal fluctuation-dissipation check --------------------------")
    print("S(q) * K(q) / kT, corrected for the explicit-Euler bias.")
    print("Must be 1.0 for any correct update; 1/m(q) is the signature of")
    print("driving the limit surface with the control-point gradient.")
    print(f"{'q':>9} {'m(q)':>7} {'lam/step':>10} {'S*K/kT':>9} {'x m(q)':>8} "
          f"{'K_S/(A kc q^4)':>15}")
    qb, rb, _, _ = ms.radial_average(q, fdt_ratio, n_bins=args.bins)
    _, mb, _, _ = ms.radial_average(q, m, n_bins=args.bins)
    _, lb, _, _ = ms.radial_average(q, lam_step, n_bins=args.bins)
    _, kb, _, _ = ms.radial_average(q, K_S / (lat.area * par["kc"] *
                                              np.where(q > 0, q, 1) ** 4),
                                    n_bins=args.bins)
    for a, b, c, d, e in zip(qb, mb, lb, rb, kb):
        print(f"{a:>9.4f} {b:>7.3f} {c:>10.5f} {d:>9.3f} {d*b:>8.3f} {e:>15.3f}")

    ok = fit_mask & np.isfinite(rb)
    print(f"\nmedian S*K/kT over the continuum window: {np.median(rb[ok]):.3f}"
          f"   (x m(q): {np.median((rb*mb)[ok]):.3f})")

    # ------------------------------------------------------------- relaxation
    lam_local = par["diff"] * a_cell / par["kbt"]      # nm^3 / (pN us)
    g_local = lam_local * par["kc"] * q ** 4
    g_oseen = par["kc"] * q ** 3 / (4.0 * args.eta)

    good = np.isfinite(gamma) & (q > 0)
    qg, gg, gge, gcnt = ms.radial_average(q[good], gamma[good], n_bins=args.bins)
    _, gl, _, _ = ms.radial_average(q[good], g_local[good], n_bins=args.bins)
    _, go, _, _ = ms.radial_average(q[good], g_oseen[good], n_bins=args.bins)
    dyn_mask = qg * lat.d_face_x <= args.qdx_max
    dyn_slope = ms.log_slope(qg[dyn_mask], gg[dyn_mask]) if dyn_mask.sum() > 2 else float("nan")

    print("\n--- relaxation rate -------------------------------------------------")
    print(f"{'q':>9} {'Gamma (1/us)':>13} {'free-drain':>12} {'ratio':>7} "
          f"{'Oseen':>12} {'ratio':>9} {'modes':>6}")
    for a, b, c, d, e in zip(qg, gg, gl, go, gcnt):
        print(f"{a:>9.4f} {b:>13.4e} {c:>12.4e} {b/c:>7.3f} {d:>12.4e} "
              f"{b/d:>9.2e} {e:>6d}")
    print(f"\nlog-log slope of Gamma(q) over q*dx <= {args.qdx_max}: {dyn_slope:.3f}")
    print("   free draining predicts 4, Oseen hydrodynamics predicts 3")
    ok = dyn_mask & np.isfinite(gg)
    if ok.any():
        print(f"median Gamma / free-draining prediction: {np.median((gg/gl)[ok]):.3f}")
        print(f"median Gamma / Oseen prediction        : {np.median((gg/go)[ok]):.3e}")

    if args.npz:
        np.savez(args.npz, q=q, qx=qx, qy=qy, S=S, m=m, K_S=K_S, gamma=gamma,
                 lam_step=lam_step, fdt_ratio=fdt_ratio, acf=acf,
                 q_bin=qc, S_bin=sm, S_err=se, S_cnt=cnt,
                 q_gam=qg, gamma_bin=gg, gamma_local=gl, gamma_oseen=go,
                 q_legacy=q_leg, S_legacy=s_leg, rms_t=rms_t, t_all=t_all,
                 area=lat.area, nx=lat.nx, ny=lat.ny, dx=lat.d_face_x,
                 dy=lat.d_face_y, dt_frame=dt_frame, n_frames=n_frames,
                 kc_fit=kc_fit, sigma_fit=sig_fit, slope=slope,
                 slope_legacy=slope_leg, dyn_slope=dyn_slope,
                 **{k: v for k, v in par.items()})
        print(f"\nwrote {args.npz}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
