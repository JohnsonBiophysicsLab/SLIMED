#!/usr/bin/env python3
"""A minimal, exactly-solvable membrane Langevin integrator on the SLIMED lattice.

This is a reference, not a replacement for the C++ dynamics.  It integrates the
Monge-gauge Helfrich model spectrally on the same zig-zag triangular tile, so
its stationary spectrum and its relaxation rates are known in closed form:

    E            = (A/2) sum_q (kc q^4 + sigma q^2) |h_q|^2
    dh_q/dt      = -Lambda(q) (kc q^4 + sigma q^2) h_q + xi_q
    <xi xi*>     = (2 kT Lambda(q) / A) delta(t - t')
    => <|h_q|^2> = kT / (A (kc q^4 + sigma q^2))            (independent of Lambda)
       Gamma(q)  = Lambda(q) (kc q^4 + sigma q^2)

Two mobilities are provided, and they are the whole point of the second half of
this study:

    "local"  Lambda(q) = Lambda0            free draining -- Gamma ~ q^4
    "oseen"  Lambda(q) = 1 / (4 eta q)      a membrane in a solvent -- Gamma ~ q^3

The static spectrum is identical for the two; only the dynamics differ.  That
is why <|h_q|^2> can look perfect while h(t)h(t+dt) says nothing about
hydrodynamics, and why a q^-3 relaxation cannot appear in a simulation whose
mobility is a per-vertex scalar.

Uses
----
* End-to-end validation: run it, analyse the output with the same
  ``membrane_spectrum`` code used on the simulation, and check that q^-4 and
  the expected Gamma(q) come back.
* A quantitative statement of what Oseen hydrodynamics would do to the SLIMED
  relaxation times at the same kc, box and lattice.

    python analysis/reference_langevin.py --mobility oseen --steps 200000
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import membrane_spectrum as ms  # noqa: E402


def mobility(kind: str, q: np.ndarray, lambda0: float, eta: float) -> np.ndarray:
    """Lambda(q), in nm^3 / (pN us)."""
    if kind == "local":
        return np.full_like(q, lambda0)
    if kind == "oseen":
        # 1/(4 eta q): the zz response of a flat sheet in an unbounded solvent,
        # after the in-plane average of the Oseen tensor.  Diverges at q -> 0,
        # so the q = 0 mode (uniform translation) is simply not integrated.
        out = np.zeros_like(q)
        nz = q > 0
        out[nz] = 1.0 / (4.0 * eta * q[nz])
        return out
    raise ValueError(kind)


def run(lat: ms.Lattice, kc: float, sigma: float, kbt: float, dt: float,
        n_steps: int, sample_every: int, kind: str, lambda0: float, eta: float,
        seed: int = 0, burn_relaxation_times: float = 5.0):
    qx, qy, q = ms.q_grid(lat)
    n = lat.nx * lat.ny
    stiff = kc * q**4 + sigma * q**2
    lam = mobility(kind, q, lambda0, eta)
    gamma = lam * stiff                      # 1/us
    decay = 1.0 - gamma * dt

    if np.any(decay < -1.0):
        worst = q.ravel()[np.argmin(decay)]
        raise SystemExit(
            f"dt = {dt} is unstable: Gamma*dt reaches {(gamma*dt).max():.3f} "
            f"at q = {worst:.3f} /nm.  Use dt < {2.0/gamma.max():.3g} us.")

    # Noise increment on h_q, built by filtering real-space white noise so the
    # conjugate structure the zig-zag lattice needs comes out automatically.
    noise_gain = np.sqrt(2.0 * kbt * lam * dt / lat.area) * math.sqrt(n)
    noise_gain[q == 0] = 0.0

    rng = np.random.default_rng(seed)
    h_q = np.zeros((lat.ny, lat.nx), complex)

    with np.errstate(divide="ignore"):
        tau_slow = 1.0 / gamma[q > 0].min()
    n_burn = int(burn_relaxation_times * tau_slow / dt)
    n_burn = min(n_burn, n_steps)

    frames = []
    for step in range(n_steps + n_burn):
        w = rng.standard_normal((lat.ny, lat.nx))
        h_q = decay * h_q + ms.tile_fft2(w, lat) / n * noise_gain
        if step >= n_burn and (step - n_burn) % sample_every == 0:
            frames.append(np.real(ms.tile_ifft2(h_q * n, lat)))
    return np.asarray(frames), q, gamma, tau_slow, n_burn


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-face-x", type=int, default=20)
    ap.add_argument("--n-face-y", type=int, default=24)
    ap.add_argument("--l-face", type=float, default=5.0)
    ap.add_argument("--kc", type=float, default=83.4, help="pN.nm")
    ap.add_argument("--sigma", type=float, default=0.0, help="pN/nm")
    ap.add_argument("--kbt", type=float, default=4.17, help="pN.nm")
    ap.add_argument("--mobility", choices=("local", "oseen"), default="local")
    ap.add_argument("--lambda0", type=float, default=None,
                    help="local mobility in nm^3/(pN us); default D*a/kT with D = 1 "
                         "nm^2/us, matching the SLIMED free-draining update")
    ap.add_argument("--eta", type=float, default=1.0e-3,
                    help="solvent viscosity in pN.us/nm^2 (1e-3 = water)")
    ap.add_argument("--dt", type=float, default=None, help="us; default 0.2/Gamma_max")
    ap.add_argument("--steps", type=int, default=400000)
    ap.add_argument("--sample-every", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--npz", default=None)
    args = ap.parse_args(argv)

    dy = math.sqrt(3) / 2 * args.l_face
    lat = ms.Lattice(args.n_face_x, args.n_face_y, args.l_face, dy)
    a_cell = lat.area / (lat.nx * lat.ny)
    lambda0 = args.lambda0 if args.lambda0 is not None else 1.0 * a_cell / args.kbt

    qx, qy, q = ms.q_grid(lat)
    stiff = args.kc * q**4 + args.sigma * q**2
    gmax = (mobility(args.mobility, q, lambda0, args.eta) * stiff).max()
    dt = args.dt if args.dt is not None else 0.2 / gmax

    print(lat.describe())
    print(f"mobility = {args.mobility}, Lambda0 = {lambda0:.4g} nm^3/(pN.us), "
          f"eta = {args.eta} pN.us/nm^2")
    print(f"kc = {args.kc} pN.nm, sigma = {args.sigma} pN/nm, dt = {dt:.4g} us")

    frames, q, gamma, tau_slow, n_burn = run(
        lat, args.kc, args.sigma, args.kbt, dt, args.steps, args.sample_every,
        args.mobility, lambda0, args.eta, args.seed)
    dt_frame = dt * args.sample_every
    print(f"slowest mode tau = {tau_slow:.4g} us; discarded {n_burn} steps, "
          f"kept {len(frames)} frames {dt_frame:.4g} us apart "
          f"({len(frames)*dt_frame:.4g} us = {len(frames)*dt_frame/tau_slow:.1f} tau)")

    H = ms.tile_fft2(frames, lat)
    S = ms.power_spectrum(H, lat)
    ideal = np.where(q > 0, args.kbt / (lat.area * (args.kc * q**4 + args.sigma * q**2)), 0)

    acf = ms.mode_autocorrelation(H, max_lag=min(len(frames) // 4, 3000))
    meas = np.full_like(q, np.nan)
    for iy in range(lat.ny):
        for ix in range(lat.nx):
            if q[iy, ix] > 0:
                meas[iy, ix] = ms.relaxation_rate(acf[:, iy, ix], dt_frame,
                                                  n_frames=len(frames))

    qc, sm, se, cnt = ms.radial_average(q, S, n_bins=16)
    _, si, _, _ = ms.radial_average(q, ideal, n_bins=16)
    ok = np.isfinite(meas) & (q > 0)
    qg, gm, _, _ = ms.radial_average(q[ok], meas[ok], n_bins=16)
    _, gt, _, _ = ms.radial_average(q[ok], gamma[ok], n_bins=16)

    print(f"\n{'q':>9} {'<|h_q|^2>':>12} {'exact':>12} {'ratio':>7} "
          f"{'Gamma meas':>12} {'Gamma exact':>12} {'ratio':>7}")
    for a, b, c, d, e, f in zip(qc, sm, si, qg, gm, gt):
        print(f"{a:>9.4f} {b:>12.4e} {c:>12.4e} {b/c:>7.3f} "
              f"{e:>12.4e} {f:>12.4e} {e/f:>7.3f}")

    print(f"\nstatic  log-log slope: {ms.log_slope(qc, sm):.3f}   (exact "
          f"{ms.log_slope(qc, si):.3f})")
    print(f"Gamma   log-log slope: {ms.log_slope(qg, gm):.3f}   (exact "
          f"{ms.log_slope(qg, gt):.3f};  4 = free draining, 3 = Oseen)")
    kc_fit, sig_fit = ms.fit_kc_sigma(qc, sm, lat.area, args.kbt)
    print(f"fitted kc = {kc_fit:.2f} pN.nm (input {args.kc}), "
          f"sigma = {sig_fit:.4f} pN/nm (input {args.sigma})")

    if args.npz:
        np.savez(args.npz, q=q, S=S, ideal=ideal, gamma_meas=meas, gamma_exact=gamma,
                 q_bin=qc, S_bin=sm, S_exact=si, q_gam=qg, g_meas=gm, g_exact=gt,
                 mobility=args.mobility, kc=args.kc, sigma=args.sigma,
                 area=lat.area, dt_frame=dt_frame)
        print(f"wrote {args.npz}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
