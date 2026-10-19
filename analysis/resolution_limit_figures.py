#!/usr/bin/env python
"""
The two summary figures of docs/resolution_limits.md: how large a time step
and how coarse a mesh a run can take.

    python analysis/resolution_limit_figures.py ~/slimed_runs/resolution docs/figures

reads the resolution_summary.npz files that resolution_sweep_analysis.py
reduce writes (plus <root>/_ref_flat/box_* for the dt 0.001 point) and writes
resolution_timestep.png and resolution_meshsize.png to the output directory.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

KBT = 4.17
END, SLACK = 400.0, 0.25
Q_COMPARE = 0.12

# the reference palette (dataviz skill): categorical slots 1-2, the blue ramp from step 250
BLUE, ORANGE = "#2a78d6", "#eb6834"
RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#0d366b"]
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"
DANGER = "#d03b3b"


def style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 200, "font.size": 10, "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE, "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "axes.titlecolor": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
        "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
        "lines.linewidth": 2, "lines.markersize": 7, "legend.labelcolor": INK2,
    })


def plain(axis, ticks):
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    axis.set_minor_locator(NullLocator())


def load(root):
    runs = {}
    for f in list(root.glob("*/resolution_summary.npz")) + list(root.glob("_ref_flat/*/resolution_summary.npz")):
        z = np.load(f)
        runs[f.parent.name if f.parent.parent == root else "ref_" + f.parent.name] = {
            k: (z[k].item() if z[k].ndim == 0 else z[k]) for k in z.files}
    return runs


def default_config(name):
    return not name.startswith(("diag_", "zonly_"))


def survived(r):
    return r["diverged_iteration"] < 0 and r["time_end"] >= END - SLACK


def timestep_figure(runs, out):
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    bounds = {}
    for fluid, color, label in ((False, BLUE, "frozen mesh"), (True, ORANGE, "fluid mesh (ν = 0.5 /edge/µs)")):
        sel = {n: r for n, r in runs.items() if default_config(n) and r["l_face"] == 5.0 and bool(r["flips"]) == fluid}
        ref = sel["dt0.0005_" + ("nu0.5" if fluid else "frozen")]
        by_dt = {}
        for n, r in sel.items():
            by_dt.setdefault(r["dt"], []).append(r)
        xs, ys, es, bad = [], [], [], []
        for dt in sorted(by_dt):
            group = by_dt[dt]
            if any(r["diverged_iteration"] >= 0 for r in group):
                # the mean time to divergence, pooled over the seeds
                T = sum(r["diverged_time"] if r["diverged_iteration"] >= 0 else min(r["time_end"], 600) for r in group)
                bad.append((dt, T / sum(r["diverged_iteration"] >= 0 for r in group)))
                continue
            hi = ref["q"] >= 0.5
            err = [100 * (r["S"][hi] / ref["S"][hi]).mean() - 100 for r in group if survived(r)]
            xs.append(dt); ys.append(np.mean(err)); es.append(np.std(err, ddof=1) / np.sqrt(len(err)) if len(err) > 1 else 0)
        ax.plot(xs, ys, "o-", color=color, mec=SURFACE, mew=1.5, label=f"{label}: no divergence")
        bounds[fluid] = (max(xs), min(d for d, _ in bad))
        # diverged settings: a cross on the row of their mean time to divergence
        for dt, t in bad:
            ax.plot(dt, 14.2 if fluid else 15.4, marker="x", color=color, mew=2, ms=7, ls="none")
    ax.axhline(0, color=MUTED, lw=0.8)
    # the bound: the largest step with no divergence (3 seeds, 1200 us) and the next one up
    for fluid, color in ((False, BLUE), (True, ORANGE)):
        ok, bad_ = bounds[fluid]
        ax.axvspan(ok, 0.06, ymin=0.0, ymax=1.0, color=color, alpha=0.06, lw=0)
        ax.axvline(ok, color=color, lw=1, ls="--")
    ax.axvline(0.001, color=MUTED, lw=1, ls=":")
    ax.text(0.00105, 12.6, "default 0.001", color=INK2, fontsize=8, va="top")
    ax.axvline(0.035, color=INK2, lw=1, ls="-.")
    ax.text(0.0335, 9.0, "height-only motion\ndiverges here (0.03-0.04)", color=INK2, fontsize=8, ha="right")
    ax.text(0.00048, 15.4, "frozen: diverges ✕", color=BLUE, fontsize=8, va="center")
    ax.text(0.00048, 14.2, "fluid: diverges ✕", color=ORANGE, fontsize=8, va="center")
    ax.text(bounds[False][0] * 1.06, 11.6, f"frozen bound\n{bounds[False][0]:g} µs", color=BLUE, fontsize=8, va="top")
    ax.text(bounds[True][0] * 0.94, 11.6, f"fluid bound\n{bounds[True][0]:g} µs", color=ORANGE, fontsize=8, va="top", ha="right")
    ax.set_xscale("log"); ax.set_xlim(0.00045, 0.06); ax.set_ylim(-1.5, 16.5)
    plain(ax.xaxis, (0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05))
    ax.set_xlabel("timeStep (µs), at diffConst = 1 nm²/µs")
    ax.set_ylabel("error of the shortest modes (%)\n⟨|h_q|²⟩ at q ≥ 0.5 nm⁻¹ against dt 0.0005")
    ax.set_title("Time step: the bound is stability, not accuracy (lFace 5 nm)", loc="left", fontsize=11)
    ax.legend(loc="lower left", bbox_to_anchor=(0.4, 0.2), fontsize=8)
    fig.text(0.01, 0.01, "Shaded: runs diverge (1 of 3 seeds one step past the dashed line, all of them beyond).\n"
             "The bound grows about as lFace: 0.005 µs frozen, 0.004 µs fluid at lFace 10.", fontsize=7.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(out / "resolution_timestep.png"); plt.close(fig)


def shells(q, y, edges):
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (q >= a) & (q < b) & np.isfinite(y)
        if m.any():
            out.append((np.exp(np.log(q[m]).mean()), y[m].mean()))
    return np.array(out)


def meshsize_figure(runs, out):
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    mesh = sorted([(r["l_face"], n, r) for n, r in runs.items()
                   if (n.startswith("mesh") or n == "ref_box_frozen") and not r["flips"] and survived(r)], key=lambda t: t[0])
    ref = mesh[0][2]
    edges = np.geomspace(Q_COMPARE, 2.0, 14)
    kref = shells(ref["q"], KBT / (ref["area"] * ref["S"] * ref["q"] ** 4), edges)
    ax.axhspan(-5, 5, color=BLUE, alpha=0.07, lw=0)
    ax.axhline(0, color=MUTED, lw=0.8)
    for y in (5, 10):
        ax.axhline(y, color=MUTED, lw=0.8, ls="--")
        ax.text(5.6, y, f"{y} %", color=INK2, fontsize=8, va="center")
    for (l, n, r), color in zip(mesh[1:], RAMP):
        sh = shells(r["q"], KBT / (r["area"] * r["S"] * r["q"] ** 4), edges)
        kr = np.interp(np.log(sh[:, 0]), np.log(kref[:, 0]), kref[:, 1])
        ax.plot(sh[:, 0] * r["d_face"], 100 * (sh[:, 1] / kr - 1), "o-", color=color, mec=SURFACE, mew=1.2, ms=6,
                lw=1.6, label=f"{l:.3g} nm")
    ax.axvline(2.0, color=INK, lw=1)
    ax.text(1.95, 52, "bound: q·lFace = 2\nλ = π lFace ≈ 3 lFace", color=INK, fontsize=8.5, ha="right", va="top")
    ax.set_xscale("log"); ax.set_xlim(0.55, 5.4); ax.set_ylim(-15, 68)
    plain(ax.xaxis, (0.6, 1, 2, 3, 4, 5))
    ax.set_xlabel("q · lFace  (wavevector × mesh size)")
    ax.set_ylabel(f"error in the bending rigidity a mode reads (%)\nagainst the finest mesh (lFace {ref['l_face']:.3g} nm) at the same q")
    inv = lambda x: 2 * np.pi / np.maximum(x, 1e-9)
    top = ax.secondary_xaxis("top", functions=(inv, inv))
    plain(top.xaxis, (1.5, 2, 3, 5, 10))
    top.set_xlabel("wavelength / lFace", color=INK2)
    ax.set_title("Mesh size: every mesh is right to q·lFace ≈ 2", loc="left", fontsize=11, pad=34)
    ax.legend(title="lFace (frozen mesh)", title_fontsize=8, fontsize=8, loc="upper left", ncol=2)
    fig.text(0.01, 0.01, "Modes at q ≥ 0.12 nm⁻¹ (relaxed within the 300 µs window); the fluid mesh collapses on the same curve.\n"
             "At fixed diffConst a coarser mesh also relaxes faster: Γ ∝ area per vertex ∝ lFace².",
             fontsize=7.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(out / "resolution_meshsize.png"); plt.close(fig)


def main():
    root, out = Path(sys.argv[1]).expanduser(), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    style()
    runs = load(root)
    timestep_figure(runs, out)
    meshsize_figure(runs, out)


if __name__ == "__main__":
    main()
