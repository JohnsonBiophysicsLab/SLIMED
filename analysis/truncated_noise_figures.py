#!/usr/bin/env python
"""
The three figures of docs/resolution_limits.md section 2.5: what truncating the
Brownian noise (noiseTruncation) does to a run's stability and to its spectrum.

    python analysis/truncated_noise_figures.py ~/slimed_runs/trunc_noise/runs docs/figures

reads, under the runs root,

  base_dt0.004_frozen_s1/stepdiag.csv  the per-step force and noise steps of the
                                       diverging reference run (written by the
                                       study's step-diagnostics patch)
  A_dt0.004_*, C_dt0.004_*             the stability runs (run.log, input.params)
  B_dt0.002_*/resolution_summary.npz   the spectrum runs, reduced by
                                       resolution_sweep_analysis.py reduce

and writes truncated_noise_mechanism.png, truncated_noise_stability.png and
truncated_noise_spectrum.png to the output directory.
"""
import math
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

# the reference palette (dataviz skill): categorical slots 1-4 in order, text in ink tokens
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"

SEED0 = 20261012
INVERSION_STEP = 6735   # first inverted face of base_dt0.004_frozen_s1, from its per-step frames
FAILING_VERTEX = 417    # the vertex whose force step runs away first

# one entry per noise law: run-name tag, label, colour, line style
LAWS = {
    "inf": ("Gaussian (default)", BLUE, "-"),
    "c2": ("cut at 2σ, rescaled", ORANGE, "-"),
    "c2raw": ("cut at 2σ, not rescaled", AQUA, "--"),
    "cold": ("Gaussian at 0.774 kT", YELLOW, ":"),
}


def style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 200, "font.size": 10, "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE, "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "axes.titlecolor": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
        "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
        "lines.linewidth": 2, "lines.markersize": 7, "legend.labelcolor": INK2,
        "axes.titlesize": 10.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
    })


def plain(axis, ticks):
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    axis.set_minor_locator(NullLocator())


# --- slimed::standard_normal (include/Counter_rng.hpp), to recompute any kick of a run ---
_M = (1 << 64) - 1


def _splitmix64(x):
    x = (x + 0x9E3779B97F4A7C15) & _M
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & _M
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & _M
    return x ^ (x >> 31)


def standard_normal(seed, iteration, vertex, axis):
    seed, iteration, vertex, axis = int(seed), int(iteration), int(vertex), int(axis)
    key = _splitmix64(((_splitmix64(seed) + iteration) + _splitmix64((vertex * 4 + axis) & _M)) & _M)
    u1 = ((_splitmix64(key) >> 11) + 0.5) / 9007199254740992.0
    u2 = ((_splitmix64(key ^ 0xD1B54A32D192ED03) >> 11) + 0.5) / 9007199254740992.0
    return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)


def read_params(path):
    return dict(re.findall(r"^(\w+)\s*=\s*([^#\s]+)", Path(path).read_text(), re.M))


def mechanism(root, out):
    run = root / "base_dt0.004_frozen_s1"
    d = np.genfromtxt(run / "stepdiag.csv", delimiter=",", names=True)
    p = read_params(run / "input.params")
    seed, dt = int(p["randomSeed"]), float(p["timeStep"])
    scale = math.sqrt(2.0 * float(p["diffConst"]) * dt)
    keep = d["iteration"] >= 6640
    it = d["iteration"][keep].astype(int)
    own_noise = [scale * math.sqrt(sum(standard_normal(seed, i, FAILING_VERTEX, a) ** 2 for a in range(3))) for i in it]
    typical = float(np.median(d["maxForceDisp"][d["iteration"] < 6600]))

    fig, ax = plt.subplots(figsize=(7.2, 3.9), constrained_layout=True)
    ax.plot(it, d["maxForceDisp"][keep], color=BLUE, label="force step, largest over all vertices")
    ax.plot(it, d["maxNoiseDisp"][keep], color=ORANGE, label="noise step, largest over all vertices")
    ax.plot(it, own_noise, color=AQUA, lw=1.5, ls="--", label=f"noise step of vertex {FAILING_VERTEX}")
    ax.axvline(INVERSION_STEP, color=INK2, lw=1, ls=":")
    ax.text(INVERSION_STEP - 2, 3e3, f"step {INVERSION_STEP}:\nfirst face inverts", ha="right", va="top", fontsize=8, color=INK2)
    ax.axhline(typical, color=MUTED, lw=1, ls="--", label=f"median force step before step 6600 ({typical:.2f} nm)")
    ax.axhline(5.0, color=MUTED, lw=0.8)
    ax.text(6642, 5.0 * 1.15, "edge length 5 nm", fontsize=8, color=INK2)
    ax.set_yscale("log")
    ax.set_ylim(0.01, 3e4)
    ax.set_xlim(6640, int(it.max()) + 1)
    ax.set_xlabel("iteration")
    ax.set_ylabel("displacement in one step (nm)")
    ax.set_title(f"The cross-over is the force step's: frozen, dt {dt:g} µs, diverges at step {int(it.max())}")
    ax.legend(loc="upper left", fontsize=8)
    fig.savefig(out / "truncated_noise_mechanism.png")
    plt.close(fig)


def stability_runs(root):
    """{law: [(time of divergence or None, simulated time)]} over the A and C sets."""
    runs = {}
    for d in sorted(list(root.glob("A_dt0.004_*")) + list(root.glob("C_dt0.004_*"))):
        log = d / "run.log"
        if not log.exists() or "exit" not in log.read_text(errors="replace"):
            continue
        p = read_params(d / "input.params")
        dt, cap = float(p["timeStep"]), int(p["maxIterations"]) * float(p["timeStep"])
        m = re.search(r"at iteration (\d+); the run has diverged", log.read_text(errors="replace"))
        law = re.sub(r"_s\d+$", "", d.name).split("_", 2)[2]
        runs.setdefault(law, []).append((int(m.group(1)) * dt if m else None, cap))
    return runs


# where each survival curve's direct label sits: (time, fraction, alignment), clear of the curves
SURVIVAL_LABELS = {"inf": (1.5, 0.05, "left"), "c2": (27.0, 0.36, "left"),
                   "c2raw": (101.5, 0.5, "left"), "cold": (101.5, 0.19, "left")}


def kaplan_meier(events, cap):
    t = sorted(e for e in events if e is not None)
    n, s, xs, ys = len(events), 1.0, [0.0], [1.0]
    for i, ti in enumerate(t):
        xs += [ti, ti]
        ys += [s, s * (n - i - 1) / (n - i)]
        s = ys[-1]
    xs.append(cap)
    ys.append(s)
    return np.array(xs), np.array(ys)


def poisson_interval(k, exposure):
    """Garwood 95 % interval on a Poisson rate, via the Wilson-Hilferty approximation."""
    hi_k = k + 1
    lo = 0.0 if k == 0 else k * (1 - 1 / (9 * k) - 1.96 / (3 * math.sqrt(k))) ** 3
    hi = hi_k * (1 - 1 / (9 * hi_k) + 1.96 / (3 * math.sqrt(hi_k))) ** 3
    return 1000 * lo / exposure, 1000 * hi / exposure


def stability(root, out):
    runs = stability_runs(root)
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 4.0), constrained_layout=True,
                               gridspec_kw=dict(width_ratios=[1.25, 1]))
    for law, (label, color, ls) in LAWS.items():
        r = runs[law]
        x, y = kaplan_meier([e for e, _ in r], max(c for _, c in r))
        a.plot(x, y, color=color, ls=ls, label=f"{label} ({len(r)} seeds)")
        lx, ly, ha = SURVIVAL_LABELS[law]
        a.text(lx, ly, label, color=INK2, fontsize=8, va="center", ha=ha)
    a.set_xlim(0, 135)
    a.set_ylim(-0.02, 1.02)
    plain(a.xaxis, [0, 20, 40, 60, 80, 100])
    a.set_xlabel("time (µs)")
    a.set_ylabel("fraction of runs not yet diverged")
    a.set_title("Survival, frozen sheet, dt 0.004 µs")
    a.legend(loc="upper right", fontsize=8)

    order = [("inf", "Gaussian (default)"), ("c4", "cut at 4σ, rescaled"), ("c3", "cut at 3σ, rescaled"),
             ("c2.5", "cut at 2.5σ, rescaled"), ("c2", "cut at 2σ, rescaled"), ("c1", "cut at 1σ, rescaled"),
             ("c2raw", "cut at 2σ, not rescaled"), ("cold", "Gaussian at 0.774 kT")]
    for row, (law, label) in enumerate(order):
        r = runs[law]
        k = sum(e is not None for e, _ in r)
        exposure = sum(e if e is not None else c for e, c in r)
        rate = 1000 * k / exposure
        lo, hi = poisson_interval(k, exposure)
        color = LAWS[law][1] if law in LAWS else INK2
        b.plot([lo, hi], [row, row], color=color, lw=2, solid_capstyle="round")
        b.plot(rate, row, "o", color=color, mec=SURFACE, mew=1.5)
        b.text(hi * 1.12, row, f"{rate:.0f}  ({len(r)} seeds)", fontsize=8, color=INK2, va="center")
    gauss = runs["inf"]
    g_rate = 1000 * sum(e is not None for e, _ in gauss) / sum(e if e is not None else c for e, c in gauss)
    b.axvline(g_rate, color=BLUE, lw=0.8, ls="--")
    b.set_yticks(range(len(order)))
    b.set_yticklabels([label for _, label in order], fontsize=8.5)
    b.invert_yaxis()
    b.set_xscale("log")
    b.set_xlim(0.8, 600)
    plain(b.xaxis, [1, 3, 10, 30, 100, 300])
    b.grid(axis="y", visible=False)
    b.set_xlabel("divergences per ms (95 % interval)")
    b.set_title("Divergence rate, same runs")
    fig.savefig(out / "truncated_noise_stability.png")
    plt.close(fig)


# q >= 0.12 nm^-1 only: the slower modes are under-sampled in 300 us (section 2.3, Q_COMPARE)
SHELLS = np.array([0.12, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0])


def shell_spectrum(z):
    """Per-shell mean of <|h_q|^2> and its standard error over the 8 time blocks."""
    q = z["q"]
    out = []
    for lo, hi in zip(SHELLS[:-1], SHELLS[1:]):
        blocks = z["S_blocks"][:, (q >= lo) & (q < hi)].mean(axis=1)
        out.append((blocks.mean(), blocks.std(ddof=1) / math.sqrt(len(blocks))))
    return np.array(out)


def shell_rate_ratio(zs, refs):
    """Per-shell mean of Gamma(q) / Gamma_Gaussian(q), taken mode by mode so the
    q^4 spread of Gamma inside a shell cancels; error from the mode-to-mode scatter."""
    q = zs[0]["q"]
    g = np.mean([z["gamma"] for z in zs], axis=0) / np.mean([z["gamma"] for z in refs], axis=0)
    out = []
    for lo, hi in zip(SHELLS[:-1], SHELLS[1:]):
        r = g[(q >= lo) & (q < hi)]
        r = r[np.isfinite(r)]
        out.append((r.mean(), r.std(ddof=1) / math.sqrt(len(r))) if len(r) > 2 else (np.nan, np.nan))
    return np.array(out)


def spectrum(root, out):
    laws = {}
    for f in sorted(root.glob("B_dt0.002_*/resolution_summary.npz")):
        z = np.load(f)
        if "S" not in z.files:
            continue
        laws.setdefault(re.sub(r"_s\d+$", "", f.parent.name).split("_")[-1], []).append(z)

    def pooled(zs):
        v = np.array([shell_spectrum(z) for z in zs])          # seeds x shells x (mean, se)
        return v[:, :, 0].mean(axis=0), np.sqrt((v[:, :, 1] ** 2).sum(axis=0)) / len(zs)

    centres = np.sqrt(SHELLS[:-1] * SHELLS[1:])
    ref_S, ref_Se = pooled(laws["inf"])
    shown = [("c2", "cut at 2σ, rescaled", ORANGE, "-", 0.97), ("c1", "cut at 1σ, rescaled", INK2, "-", 1.03),
             ("c2raw", "cut at 2σ, not rescaled", AQUA, "--", 1.0)]
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 4.0), constrained_layout=True)
    for ax, key in ((a, "S"), (b, "gamma")):
        ax.axhline(1.0, color=BLUE, lw=1.2)
        ax.text(1.02, 1.0, "Gaussian", color=INK2, fontsize=8, va="center", transform=ax.get_yaxis_transform())
        for law, label, color, ls, shift in shown:
            if key == "S":
                m, e = pooled(laws[law])
                r = m / ref_S
                re_ = r * np.sqrt((e / m) ** 2 + (ref_Se / ref_S) ** 2)
            else:
                r, re_ = shell_rate_ratio(laws[law], laws["inf"]).T
            ok = np.isfinite(r)
            ax.errorbar(centres[ok] * shift, r[ok], yerr=re_[ok], color=color, ls=ls, marker="o", mec=SURFACE,
                        mew=1.5, capsize=0, elinewidth=1.5, label=label)
        ax.set_xscale("log")
        plain(ax.xaxis, [0.1, 0.2, 0.3, 0.5, 1.0])
        ax.set_xlabel("q (nm⁻¹)")
    a.axhline(0.7737, color=MUTED, lw=0.8, ls="--")
    a.text(0.125, 0.7737 + 0.008, "variance of the 2σ cut, 0.774", fontsize=8, color=INK2, va="bottom")
    a.set_ylim(0.60, 1.10)
    a.set_xlim(0.11, 1.05)
    plain(a.xaxis, [0.2, 0.3, 0.5, 1.0])
    a.set_ylabel("⟨|h_q|²⟩ / Gaussian")
    a.set_title("Spectrum, frozen, dt 0.002 µs, 2 seeds × 400 µs")
    a.legend(loc="lower right", fontsize=8)
    b.set_ylim(0.90, 1.10)
    b.set_xlim(0.11, 0.33)
    plain(b.xaxis, [0.15, 0.2, 0.3])
    b.set_ylabel("relaxation rate Γ(q) / Gaussian")
    b.set_title("Relaxation rates (the modes the frames can time)")
    b.legend(loc="lower right", fontsize=8)
    fig.savefig(out / "truncated_noise_spectrum.png")
    plt.close(fig)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    root, out = Path(sys.argv[1]).expanduser(), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    style()
    mechanism(root, out)
    stability(root, out)
    spectrum(root, out)
    print("wrote", ", ".join(f"truncated_noise_{n}.png" for n in ("mechanism", "stability", "spectrum")), "to", out)


if __name__ == "__main__":
    main()
