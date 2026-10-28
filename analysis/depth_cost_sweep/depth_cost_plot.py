"""Step cost against the irregular-patch recursion depth (docs/irregular_patch_results.md section 10).

Reads results/results_rep1.txt and results/results_rep2.txt, writes depth_cost_plot.svg and .png
beside this script. Per-step cost = (t_2000 - t_1000) / 1000 steps.
"""
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(path):
    t = {}
    for line in open(path):
        m = re.match(r"(\S+) (?:load=\S+ )?real \S+ user (\S+)", line)
        if m:
            t[m.group(1)] = float(m.group(2))
    return {v: (t[v + "_n2000"] - t[v + "_n1000"]) for v in
            [k[:-6] for k in t if k.endswith("_n1000")]}


reps = [load(os.path.join(HERE, "results", f"results_rep{i}.txt")) for i in (1, 2)]
uniform = [1, 2, 3, 4, 5, 6, 8, 10, 12]
fit_depths = [d for d in uniform if d >= 3]
y_fit = np.array([reps[0][f"d{d}"] for d in fit_depths])
b, a = np.polyfit(fit_depths, y_fit, 1)
frozen = np.mean([r["frozen"] for r in reps])
# effective depth of a per-valence quadrature: (force-pass kernel calls - fixed part) / (3 x corners)
eff = {"default": (25803 - 1092) / (3 * 1629), "legacy": (38666 - 1092) / (3 * 1638)}

BLUE, ORANGE, GRAY, LIGHT = "#2a78d6", "#eb6834", "#5f5e5a", "#b4b2a9"
plt.rcParams.update({"font.size": 11, "font.family": "sans-serif", "svg.fonttype": "none"})
fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=200)

x = np.linspace(0, 12.4, 2)
ax.plot(x, a + b * x, color=GRAY, lw=1.8, zorder=1,
        label=f"fit, depths 3 to 12: {a:.1f} + {b:.1f} x depth ms")
ax.axhline(frozen, color="#1baf7a", lw=2.0, ls=(0, (6, 3)), zorder=2, label=f"frozen reference (regular lattice, no flips), {frozen:.1f} ms per step")
ax.plot(0, a, "o", ms=7, mfc="white", mec=GRAY, mew=1.8, zorder=3)

for i, r in enumerate(reps):
    ax.plot(fit_depths, [r[f"d{d}"] for d in fit_depths], "o", ms=6.5, color=BLUE, mec="white", mew=0.8,
            zorder=4, label="measured, two repetitions" if i == 0 else None)
    ax.plot([1, 2], [r["d1"], r["d2"]], "o", ms=6.5, mfc="white", mec=BLUE, mew=1.8, zorder=4,
            label="depth 1, 2: flips differ, not fitted" if i == 0 else None)

ax.plot(eff["default"], np.mean([r["default"] for r in reps]), "s", ms=8, color=ORANGE, mec="white", mew=0.8,
        zorder=5, label="default 6/5/5/5, tail correction")
ax.plot(eff["legacy"], np.mean([r["legacy"] for r in reps]), "D", ms=7.5, color=ORANGE, mec="white", mew=0.8,
        zorder=5, label="legacy 12/8/7/7")

ax.set_xlim(-0.3, 12.5)
ax.set_ylim(0, 152)
ax.set_xticks(range(0, 13, 2))
ax.set_xlabel("recursion depth (levels per extraordinary corner)")
ax.set_ylabel("ms per step (steps 1000 to 2000, one core)")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.grid(axis="y", color="#e1e0d9", lw=0.8)
ax.set_axisbelow(True)
ax.legend(loc="upper left", frameon=False, fontsize=9.5, handletextpad=0.6)
fig.tight_layout()
for ext in ("svg", "png"):
    fig.savefig(os.path.join(HERE, f"depth_cost_plot.{ext}"))
print(f"fit a={a:.2f} b={b:.3f} frozen={frozen:.2f} eff={eff}")
