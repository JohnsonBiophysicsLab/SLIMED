#!/usr/bin/env python3
"""Figures from one or more ``fluctuation_report.py --npz`` outputs.

    python analysis/fluctuation_plots.py --out-dir figs \
        --spectrum "fixed=pure_a.npz" --spectrum "as-shipped=legacy_a.npz" \
        --reference "free draining=ref_local.npz" \
        --reference "Oseen=ref_oseen.npz"

Produces, as SVG so they stay sharp when embedded:

  spectrum.svg   <|h_q|^2> against |q|, with kT/(A kc q^4) and the fit
  fdt.svg        S(q) K(q) / kT, the update's own consistency, against 1/m(q)
  dynamics.svg   Gamma(q), against the q^4 and q^3 laws

Needs matplotlib; everything else in this directory does not.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import membrane_spectrum as ms  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

INK = "#1b1b1f"
GRID = "#d8d8de"
SERIES = ["#2b6cb0", "#c53030", "#2f855a", "#b7791f", "#6b46c1", "#0d7490"]


def style(ax, xlabel, ylabel, title=None):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, loc="left", fontsize=11, color=INK)
    ax.grid(True, which="both", lw=0.4, color=GRID)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK, labelsize=9)


def load(spec):
    label, _, path = spec.partition("=")
    if not path:
        label, path = os.path.splitext(os.path.basename(spec))[0], spec
    return label, np.load(path, allow_pickle=False)


def figure_spectrum(runs, out):
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    qmin = qmax = None
    for i, (label, d) in enumerate(runs):
        c = SERIES[i % len(SERIES)]
        q, s, e = d["q_bin"], d["S_bin"], d["S_err"]
        ax.errorbar(q, s, yerr=e, fmt="o", ms=4.5, color=c, lw=0, elinewidth=1,
                    capsize=2, label=f"{label}  (slope {float(d['slope']):.2f})")
        qmin = q.min() if qmin is None else min(qmin, q.min())
        qmax = q.max() if qmax is None else max(qmax, q.max())
        if i == 0:
            ql, sl = d["q_legacy"], d["S_legacy"]
            keep = (ql > 0) & (sl > 0)
            ax.plot(ql[keep], sl[keep], ".", ms=2.5, color="#9aa0a6", alpha=0.55,
                    label=f"same frames, legacy reduction "
                          f"(slope {float(d['slope_legacy']):.2f})", zorder=1)

    d0 = runs[0][1]
    qq = np.logspace(np.log10(qmin), np.log10(qmax), 100)
    th = float(d0["kbt"]) / (float(d0["area"]) * float(d0["kc"]) * qq ** 4)
    ax.plot(qq, th, "-", color=INK, lw=1.6, zorder=0,
            label=r"$k_BT/(A\,\kappa_c q^4)$, $\kappa_c$ = input")

    style(ax, r"$|q|$  (nm$^{-1}$)", r"$\langle |h_q|^2\rangle$  (nm$^4$)",
          "Height fluctuation spectrum")
    ax.legend(fontsize=8.5, frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(out, format="svg")
    plt.close(fig)
    print("wrote", out)


def figure_fdt(runs, out, bins=14):
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for i, (label, d) in enumerate(runs):
        c = SERIES[i % len(SERIES)]
        q2, r2 = d["q"], d["fdt_ratio"]
        good = np.isfinite(r2) & (q2 > 0)
        qb, rb, eb, _ = ms.radial_average(q2[good], r2[good], n_bins=bins)
        ax.errorbar(qb, rb, yerr=eb, fmt="o-", ms=4, lw=1.2, color=c,
                    elinewidth=0.8, capsize=2, label=label)
    d0 = runs[0][1]
    q2, m2 = d0["q"], d0["m"]
    good = q2 > 0
    qb, mb, _, _ = ms.radial_average(q2[good], 1.0 / m2[good], n_bins=bins)
    ax.plot(qb, mb, "--", color=INK, lw=1.4, label=r"$1/m(q)$, the limit-mask symbol")
    ax.axhline(1.0, color=INK, lw=1.0)
    ax.set_yscale("linear")
    style(ax, r"$|q|$  (nm$^{-1}$)", r"$S(q)\,K(q)/k_BT$",
          "Does the Brownian update sample exp(-E/kT)?  (1.0 = yes)")
    ax.set_yscale("linear")
    ax.legend(fontsize=8.5, frameon=False)
    fig.tight_layout()
    fig.savefig(out, format="svg")
    plt.close(fig)
    print("wrote", out)


def figure_dynamics(runs, refs, out):
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    allq = []
    for i, (label, d) in enumerate(runs):
        c = SERIES[i % len(SERIES)]
        q, g = d["q_gam"], d["gamma_bin"]
        ax.plot(q, g, "o", ms=5, color=c,
                label=f"SLIMED {label}  (slope {float(d['dyn_slope']):.2f})")
        allq.append(q)
    for j, (label, d) in enumerate(refs):
        c = SERIES[(len(runs) + j) % len(SERIES)]
        q, g = d["q_gam"], d["g_meas"]
        sl = ms.log_slope(q, g)
        ax.plot(q, g, "s", ms=4, mfc="none", color=c,
                label=f"reference integrator, {label}  (slope {sl:.2f})")
        allq.append(q)

    d0 = runs[0][1] if runs else refs[0][1]
    qq = np.logspace(np.log10(min(a.min() for a in allq)),
                     np.log10(max(a.max() for a in allq)), 100)
    if runs:
        gl = np.interp(qq, d0["q_gam"], d0["gamma_local"])
        go = np.interp(qq, d0["q_gam"], d0["gamma_oseen"])
        ax.plot(qq, gl, "-", color=INK, lw=1.5,
                label=r"free draining  $\Gamma=(Da/k_BT)\kappa_c q^4$")
        ax.plot(qq, go, "-.", color=INK, lw=1.5,
                label=r"Oseen  $\Gamma=\kappa_c q^3/4\eta$  ($\eta$ = water)")
    style(ax, r"$|q|$  (nm$^{-1}$)", r"$\Gamma(q)$  ($\mu$s$^{-1}$)",
          "Mode relaxation rate: q^4 (free draining) vs q^3 (hydrodynamics)")
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(out, format="svg")
    plt.close(fig)
    print("wrote", out)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--spectrum", action="append", default=[],
                    help="label=path.npz from fluctuation_report.py")
    ap.add_argument("--reference", action="append", default=[],
                    help="label=path.npz from reference_langevin.py")
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    runs = [load(s) for s in args.spectrum]
    refs = [load(s) for s in args.reference]

    if runs:
        figure_spectrum(runs, os.path.join(args.out_dir, "spectrum.svg"))
        figure_fdt(runs, os.path.join(args.out_dir, "fdt.svg"))
    if runs or refs:
        figure_dynamics(runs, refs, os.path.join(args.out_dir, "dynamics.svg"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
