#!/usr/bin/env python
"""
Choose the vertex diffusion constant D so that the simulation's relaxation
rates match the hydrodynamic (Oseen) prediction at the long wavelengths, and
write the input.params that runs it.

The physics
-----------
SLIMED's vertices are free-draining Brownian particles: a vertex of area a
moves with mobility D/kT and no hydrodynamic coupling to the others. A height
mode h_q then relaxes at

    Gamma_FD(q) = (D a / kT) (kc q^4 + sigma q^2) ,

while a membrane in a solvent of viscosity eta on both sides, coupled through
the Oseen tensor, relaxes at (Brochard & Lennon 1975; Seifert 1997)

    Gamma_Os(q) = (kc q^3 + sigma q) / (4 eta) .

Their ratio is

    Gamma_FD / Gamma_Os = 4 eta D a q / kT ,

independent of kc and of sigma. A single D can therefore match the q^3 law at
exactly one wavevector q*, whatever the tension:

    D* = kT / (4 eta a q*) .

Above q* the simulation relaxes faster than hydrodynamics by q/q*, below it
slower by q*/q. The criteria below choose q*:

    lowest   the lowest mode of the tile, 2 pi / max(Lx, Ly): the slowest mode,
             the one that sets the equilibration time, is then exact.
    q=VALUE  a wavevector of your choice, 1/nm.
    band     the least-squares fit of log Gamma over the modes in a band of
             |q|: q* is the geometric mean of their |q| (every mode one point).

The time step
-------------
The vertex update uses D and dt only through their product, in the drift
(D dt / kT) F and in the noise sqrt(2 D dt). Its stability is a statement
about D dt, so a D larger by a factor c needs a dt smaller by c. Everything
else that is a rate per microsecond -- the flip attempt rate, the box
mobility, the lipid reservoir rate, an explicit lipid diffusivity, the bead's
velocity -- is multiplied by c, so that every step does exactly what it did
before. The number of steps per relaxation time is unchanged: D sets the
clock, not the cost.

Usage
-----
    python analysis/oseen_matched_diffusion.py --params RUN/input.params
    python analysis/oseen_matched_diffusion.py --params RUN/input.params --match band --band 0.06 0.2
    python analysis/oseen_matched_diffusion.py --params RUN/input.params --out NEW_DIR --steps 2400000 --seed 20261016

Without --out it only reports. eta is in mPa s (1.0: water at 20 C; 0.69 at
37 C).

Fixing the time step instead
----------------------------
D dt has a ceiling: on the 100 nm, lFace = 5 frozen-box sheet D dt = 0.002
nm^2 is stable with a 2.5 % excess of bending energy and 0.005 diverges
(measured 2026-10-16). Water at the mesh scale is fast -- a 5 nm bending mode
relaxes in 4 eta / (kc q^3) ~ 0.2 ns -- so matching water forces dt ~ 1e-6 us.
With --dt the script works the other way round: it takes the largest stable D
at that dt and reports the solvent viscosity that D matches at q*,

    eta_eff = kT / (4 D a q*) ,

and the factor eta_eff / eta by which the run is slower than in that solvent:
simulated times divided by it are the solvent's times.
Brownian dynamics without hydrodynamics has no other time scale, so a run at
(D, dt) *is* the matched run with its clock relabelled by that factor.

    python analysis/oseen_matched_diffusion.py --params RUN/input.params --dt 0.001
"""
import argparse
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import membrane_spectrum as ms  # noqa: E402
import membrane_resample as mr  # noqa: E402

MPAS_TO_PN_US_PER_NM2 = 1.0e-3   # 1 mPa s = 1e-3 pN us / nm^2
DDT_MAX = 0.002   # nm^2: the stable ceiling of D dt measured on the 100 nm, lFace = 5 sheet

# Rates per microsecond that must follow D for the per-step dynamics to stay put.
RATE_KEYS = ["edgeFlipAttemptRate", "elasticBoxMobility", "lipidReservoirRate",
             "beadVelocityX", "beadVelocityY", "beadVelocityZ", "externalForceRate"]


def read_params(path):
    out = {}
    for line in open(path):
        if "=" in line and not line.strip().startswith("#"):
            k, v = line.split("=", 1)
            out[k.strip()] = v.split("#")[0].strip()
    return out


def tile(params_path):
    """The tile's lattice, its modes' |q| (one per conjugate pair, q = 0 out) and the area per vertex."""
    lat = ms.lattice_from_params(params_path)
    _, _, q2d = mr.cartesian_q_grid(lat)
    keep = np.zeros((lat.ny, lat.nx), bool)
    seen = set()
    for iy in range(lat.ny):
        for ix in range(lat.nx):
            if (iy, ix) == (0, 0) or ((-iy) % lat.ny, (-ix) % lat.nx) in seen:
                continue
            seen.add((iy, ix))
            keep[iy, ix] = True
    q = np.sort(q2d[keep])
    area_per_vertex = lat.area / (lat.nx * lat.ny)
    return lat, q, area_per_vertex


def matched_diffusion(kbt, eta, a, q_star):
    """D* = kT / (4 eta a q*), nm^2/us, with eta in pN us / nm^2."""
    return kbt / (4.0 * eta * a * q_star)


def choose_q_star(q, match, band=None, q_value=None):
    if match == "lowest":
        return q[0], "the lowest mode of the tile"
    if match == "q":
        return q_value, "a chosen wavevector"
    if match == "band":
        lo, hi = band
        sel = q[(q >= lo) & (q <= hi)]
        if sel.size == 0:
            raise SystemExit(f"no mode in the band [{lo}, {hi}] /nm")
        return float(np.exp(np.mean(np.log(sel)))), f"the geometric mean of {sel.size} modes in [{lo}, {hi}] /nm"
    raise SystemExit(f"unknown criterion {match}")


def rescaled_params(text, c, d_new, dt_new, steps=None, seed=None, note=""):
    """The parameter file with D, dt and every rate per microsecond rescaled; comments kept."""
    lines = []
    for line in text.splitlines():
        m = re.match(r"^(\s*)([A-Za-z_][A-Za-z0-9_]*)(\s*=\s*)([^#\s]+)(\s*)(#.*)?$", line)
        if not m:
            lines.append(line)
            continue
        indent, key, eq, value, gap, comment = m.groups()
        comment = comment or ""
        new = None
        if key == "diffConst":
            new = f"{d_new:.6g}"
            comment = "# nm^2/us: matched to the Oseen rate (analysis/oseen_matched_diffusion.py)"
        elif key == "timeStep":
            new = f"{dt_new:.6g}"
            comment = f"# us: the reference dt / {c:.6g}, so that D dt is unchanged"
        elif key in RATE_KEYS:
            try:
                new = f"{float(value) * c:.6g}"
                comment = f"# x {c:.6g}: a rate per us follows D"
            except ValueError:
                pass
        elif key == "lipidFlowDiffusivity":
            try:
                if float(value) > 0:
                    new = f"{float(value) * c:.6g}"
                    comment = f"# x {c:.6g}: a rate per us follows D"
            except ValueError:
                pass
        elif key == "maxIterations" and steps is not None:
            new = str(steps)
        elif key == "randomSeed" and seed is not None:
            new = str(seed)
        lines.append(f"{indent}{key}{eq}{new if new is not None else value}{gap if gap else ' '}{comment}".rstrip())
    header = ["############################################",
              "## Diffusion constant matched to the Oseen rate at long wavelength",
              "## (analysis/oseen_matched_diffusion.py)" + (f": {note}" if note else ""),
              "############################################", ""]
    return "\n".join(header + lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--params", required=True, help="the reference run's input.params")
    parser.add_argument("--eta", type=float, default=1.0, help="solvent viscosity, mPa s (default water, 1.0)")
    parser.add_argument("--match", choices=["lowest", "q", "band"], default="lowest")
    parser.add_argument("--q", type=float, default=None, help="q* in 1/nm for --match q")
    parser.add_argument("--band", type=float, nargs=2, default=None, metavar=("QMIN", "QMAX"),
                        help="the band of |q| for --match band, 1/nm (default: the lowest mode to pi/(2 dFace))")
    parser.add_argument("--out", default=None, help="write the matched run's input.params into this directory")
    parser.add_argument("--steps", type=int, default=None, help="maxIterations of the matched run")
    parser.add_argument("--seed", type=int, default=None, help="randomSeed of the matched run")
    parser.add_argument("--dt", type=float, default=None,
                        help="fix the time step (us) instead: report the largest stable D and the viscosity it matches")
    args = parser.parse_args()

    p = read_params(args.params)
    kbt = float(p.get("KBT", 4.17))
    kc = float(p["kcMembraneBending"])
    d0 = float(p.get("diffConst", 0.01))
    dt0 = float(p["timeStep"])
    eta = args.eta * MPAS_TO_PN_US_PER_NM2
    lat, q, a = tile(args.params)
    band = args.band or (q[0], np.pi / (2.0 * lat.d_face_x))
    if args.match == "q" and args.q is None:
        raise SystemExit("--match q needs --q")
    q_star, why = choose_q_star(q, args.match, band, args.q)
    d_star = matched_diffusion(kbt, eta, a, q_star)
    c = d_star / d0
    dt_star = dt0 / c

    q_cross_now = kbt / (4.0 * eta * d0 * a)
    tau_slow_os = 4.0 * eta / (kc * q[0] ** 3)          # the slowest mode under Oseen at sigma = 0
    print(f"tile {lat.lx:.1f} x {lat.ly:.1f} nm, {lat.nx * lat.ny} vertices, a = {a:.3f} nm^2 per vertex")
    print(f"kc = {kc} pN nm, kT = {kbt} pN nm, eta = {args.eta} mPa s = {eta:.3g} pN us/nm^2")
    print(f"modes: |q| from {q[0]:.4f} to {q[-1]:.4f} /nm; continuum window up to pi/(2 dFace) = {np.pi / (2 * lat.d_face_x):.4f}")
    print()
    print(f"reference run: D = {d0:g} nm^2/us, dt = {dt0:g} us, D dt = {d0 * dt0:.4g} nm^2")
    print(f"  free-draining meets Oseen at q = kT/(4 eta D a) = {q_cross_now:.3g} /nm;"
          f" at the lowest mode the simulation is {q_cross_now / q[0]:.0f}x slower than hydrodynamics")
    print()
    print(f"matched at q* = {q_star:.4f} /nm ({why}):")
    print(f"  D* = kT / (4 eta a q*) = {d_star:.4g} nm^2/us  (x {c:.4g})")
    print(f"  dt* = {dt_star:.4g} us  (D* dt* = {d_star * dt_star:.4g} nm^2, unchanged)")
    print(f"  slowest mode under Oseen: tau = 4 eta / (kc q_min^3) = {tau_slow_os:.4g} us = {tau_slow_os / dt_star:,.0f} steps")
    print()
    print("  q (1/nm)   Gamma_FD/Gamma_Os   (= q/q*, at any kc and sigma)")
    for qq in sorted({q[0], q_star, *q[np.searchsorted(q, [0.1, 0.2, np.pi / (2 * lat.d_face_x)]).clip(0, q.size - 1)]}):
        print(f"  {qq:8.4f}   {4.0 * eta * d_star * a * qq / kbt:8.3f}")
    rates = {k: float(p[k]) for k in RATE_KEYS if k in p and re.match(r"^[-+0-9.eE]+$", p[k])}
    if rates:
        print("\n  rates per us scaled by the same factor: " + ", ".join(f"{k} {v:g} -> {v * c:.4g}" for k, v in rates.items()))
    if p.get("lipidFlowMode", "off") != "off" and float(p.get("lipidFlowDiffusivity", -1)) < 0:
        print("  lipidFlowDiffusivity = -1 follows D by itself (K a D / kT)")

    if d_star * dt_star > DDT_MAX:
        print(f"\n  WARNING: D* dt* = {d_star * dt_star:.3g} nm^2 exceeds the measured ceiling {DDT_MAX} nm^2")

    if args.dt is not None:
        d_max = DDT_MAX / args.dt
        eta_eff = kbt / (4.0 * d_max * a * q_star)
        print(f"\nat a fixed dt = {args.dt:g} us: the largest stable D is D dt / dt = {d_max:.4g} nm^2/us, "
              f"which matches Oseen at q* = {q_star:.4f} /nm for eta = {eta_eff / MPAS_TO_PN_US_PER_NM2:.4g} mPa s")
        print(f"  the simulation runs {eta_eff / eta:.4g}x slower than the {args.eta:g} mPa s solvent: "
              f"divide simulated times by {eta_eff / eta:.4g} to read them in that solvent's microseconds")
        print(f"  the dt that matches eta = {args.eta:g} mPa s at this D dt is {DDT_MAX / d_star:.3g} us")

    if args.out:
        os.makedirs(args.out, exist_ok=True)
        note = f"eta = {args.eta} mPa s, q* = {q_star:.4f} /nm ({why}), D* = {d_star:.4g}, x {c:.4g} the reference"
        text = rescaled_params(open(args.params).read(), c, d_star, dt_star, args.steps, args.seed, note)
        with open(os.path.join(args.out, "input.params"), "w") as handle:
            handle.write(text)
        print(f"\nwrote {os.path.join(args.out, 'input.params')}")


if __name__ == "__main__":
    main()
