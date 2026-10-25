#!/usr/bin/env python
"""
Reduce and compare the runs of analysis/substep_resolution_runs.py.

    python analysis/substep_resolution_report.py reduce <run dir> [...]
    python analysis/substep_resolution_report.py report <root> [--figures DIR] [--fluid-ref DIR]

`reduce` is resolution_sweep_analysis.reduce_run with its averaging window scaled
with the run: 100-400 us x (L/5)^2 for a spec run, 4-40 us x (L/5)^2 for a dyn run,
so every lFace is averaged over the same number of relaxation times.

`report` prints, and writes to <root>/report.md and report.json:
  - the divergence rate of every (lFace, step, substeps) setting;
  - per lFace and step, S(q) and Gamma(q) against the dt = 0.001 (5/L)^-2 reference of the
    same lFace, in bins of q lFace; Gamma from the spec runs where 0.2 us frames resolve
    it and from the dyn runs (0.02 us frames) where they do not;
  - the same errors against Gamma_th(q) dt, Gamma_th = D a kappa_c q^4 / kT with a the
    area per vertex, which is what the Euler-Maruyama bias depends on (Gamma dt / 2);
  - the bending rigidity each mode of the lFace 10 and 20 references reads against the
    lFace 5 reference at the same q (the mesh error);
  - the fluid lFace-5 runs against --fluid-ref (default resolution/dt0.0005_nu0.5) and
    their own dyn reference.
With --figures it draws substep_resolution.png there.
"""
import argparse
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import resolution_sweep_analysis as rsa   # noqa: E402

KBT = 4.17
KC = 83.4
D = 1.0
BINS = [(0.5, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 2.5), (2.5, 3.5), (3.5, 6.0)]
TOL = 0.05
NAME = re.compile(r"^(?P<kind>f?(?:stab|spec|dyn))_l(?P<l>[\d.]+)_f(?P<f>[\d.]+)(?P<sub>_sub[\d.]+)?_s(?P<seed>\d+)$")


def params(d):
    return rsa.read_params(Path(d) / "input.params")


def reduce(dirs):
    for d in map(Path, dirs):
        summary, log = d / "resolution_summary.npz", d / "run.log"
        if summary.exists() and log.exists() and summary.stat().st_mtime > log.stat().st_mtime:
            continue   # already reduced from this run
        p = params(d)
        scale = (float(p["lFace"]) / 5.0) ** 2
        dyn = "dyn_" in d.name
        rsa.BURN, rsa.END = (4.0 * scale, 40.0 * scale) if dyn else (100.0 * scale, 400.0 * scale)
        rsa.FRAME_SLACK = 0.25 * scale * (0.1 if dyn else 1.0)
        rsa.reduce_run(d)


def divergence(d):
    log = (d / "run.log").read_text(errors="replace") if (d / "run.log").exists() else ""
    if "exit" not in log:
        return None, None
    p = params(d)
    dt, nmax = float(p["timeStep"]), int(p["maxIterations"])
    m = re.search(r" at iteration (\d+); the run has diverged", log)
    return (int(m.group(1)) * dt if m else None), nmax * dt


def rate_row(rows):
    div = [t for t, _ in rows if t is not None]
    T = sum(t if t is not None else T0 for t, T0 in rows)
    k = len(div)
    hi = 1000.0 / T * (k + 1) * (1 - 1 / (9 * (k + 1)) + 1.96 / (3 * math.sqrt(k + 1))) ** 3
    return dict(runs=len(rows), diverged=k, simulated_us=T, per_ms=1000.0 * k / T, upper95=hi,
                at_us=sorted(round(t, 2) for t in div))


def load_summary(d):
    f = Path(d) / "resolution_summary.npz"
    if not f.exists():
        return None
    z = np.load(f, allow_pickle=True)
    return z if "S" in z.files else None


def binned_ratio(zs, refs, key, l_face, dt_frame_min=None):
    """Per q lFace bin: mean over seeds of the mean per-mode ratio to the pooled reference."""
    q = refs[0]["q"]
    ref = np.nanmean([r[key] for r in refs], axis=0)
    out = []
    for a, b in BINS:
        m = (q * l_face >= a) & (q * l_face < b) & np.isfinite(ref)
        per_seed = []
        for z in zs:
            mm = m & np.isfinite(z[key])
            if mm.sum() >= 2:
                per_seed.append(float(np.mean(z[key][mm] / ref[mm])))
        if per_seed:
            err = float(np.std(per_seed, ddof=1) / math.sqrt(len(per_seed))) if len(per_seed) > 1 else float("nan")
            out.append(dict(bin=[a, b], ratio=float(np.mean(per_seed)), err=err, per_seed=per_seed, modes=int(m.sum()),
                            q_mid=float(np.mean(q[m])) if m.any() else float("nan")))
    return out


def gamma_theory(q, l_face):
    a = math.sqrt(3.0) / 2.0 * l_face ** 2
    return D * a * KC * q ** 4 / KBT


def report(root, figures, fluid_ref):
    root = Path(root).expanduser()
    groups = defaultdict(list)
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        m = NAME.match(d.name)
        if m:
            groups[(m["kind"] + (m["sub"] or ""), float(m["l"]), float(m["f"]))].append(d)
        elif d.name.startswith("pre_"):
            l_face, dt = re.match(r"pre_l([\d.]+)_dt([\d.]+)_s\d+", d.name).groups()
            groups[("pre", float(l_face), float(dt))].append(d)
    out = dict(stability=[], accuracy=[], mesh=[], fluid=[])
    lines = ["# Substep resolution benchmark", ""]

    # 1. stability
    lines += ["## Divergence", "", "Per ms of simulated time and per million steps (the walls fail per step).", "",
              "| set | lFace | dt (us) | substeps | runs | diverged | simulated (us) | per ms | 95 % upper | per 1e6 steps |",
              "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for (kind, l_face, f), ds in sorted(groups.items()):
        rows = [r for r in map(divergence, ds) if r[1] is not None]
        if not rows:
            continue
        p = params(ds[0])
        row = rate_row(rows)
        row.update(set=kind, l_face=l_face, dt=float(p["timeStep"]), substeps=int(p.get("meshForceSubsteps", 1)))
        row["per_mstep"] = row["diverged"] / (row["simulated_us"] / row["dt"]) * 1e6
        out["stability"].append(row)
        lines.append(f"| {kind} | {l_face:g} | {row['dt']:g} | {row['substeps']} | {row['runs']} | {row['diverged']} | "
                     f"{row['simulated_us']:.0f} | {row['per_ms']:.2f} | {row['upper95']:.2f} | {row['per_mstep']:.1f} |")

    # 2. accuracy against the same lFace's finest step
    def summaries(kind, l_face, f):
        return [z for z in map(load_summary, groups.get((kind, l_face, f), [])) if z is not None]

    lines += ["", "## Accuracy against the dt = 0.001 (L/5)^2 reference of the same lFace", "",
              "Ratio, mean over seeds of the per-mode ratio, by q lFace bin. Gamma: spec runs (0.2 us frames) "
              "below q lFace 1, dyn runs (0.02 us frames x (L/5)^2) above.", ""]
    for kind_spec, kind_dyn, label in (("spec", "dyn", "frozen"), ("fspec", "fdyn", "fluid")):
        for l_face in sorted({k[1] for k in groups if k[0] == kind_spec}):
            fs = sorted(k[2] for k in groups if k[0] == kind_spec and k[1] == l_face)
            refs = summaries(kind_spec, l_face, 0.001)
            refs_dyn = summaries(kind_dyn, l_face, 0.001)
            if kind_spec == "fspec":
                refs = [z for z in [load_summary(fluid_ref)] if z is not None]
            if not refs:
                continue
            for f in fs:
                if kind_spec == "spec" and f == 0.001:
                    continue
                zs, zd = summaries(kind_spec, l_face, f), summaries(kind_dyn, l_face, f)
                if not zs:
                    continue
                dt = float(zs[0]["dt"])
                row = dict(mesh=label, l_face=l_face, f=f, dt=dt, seeds=len(zs),
                           substeps=int(params(groups[(kind_spec, l_face, f)][0]).get("meshForceSubsteps", 1)),
                           S=binned_ratio(zs, refs, "S", l_face),
                           # Gamma: the 400 us runs below q lFace 1, the 36 us dense-frame runs above
                           # (their window is too short for the slower modes, which read as noise)
                           G_slow=[b for b in binned_ratio(zs, refs, "gamma", l_face) if b["bin"][1] <= 1.0],
                           G_fast=[b for b in binned_ratio(zd, refs_dyn, "gamma", l_face) if b["bin"][0] >= 1.0]
                           if (zd and refs_dyn) else [],
                           E_curv=float(np.mean([z["E_curvature"][0] for z in zs])),
                           E_mesh=float(np.mean([z["E_mesh"][0] for z in zs])))
                # worst |error| inside the resolved window q lFace <= 2
                inside = [b for b in row["S"] + row["G_slow"] + row["G_fast"] if b["bin"][1] <= 2.0]
                row["worst_inside"] = max((abs(b["ratio"] - 1) for b in inside), default=float("nan"))
                row["gamma_th_dt_at_qL2"] = gamma_theory(2.0 / l_face, l_face) * dt
                out["accuracy"].append(row)

    lines += ["| mesh | lFace | dt (us) | n | Gamma_th(2/L) dt | S by q lFace bin | Gamma (slow / fast) | worst inside q lFace <= 2 |",
              "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    fmt = lambda bs: " ".join(f"[{b['bin'][0]:g}-{b['bin'][1]:g}] {b['ratio']:.3f}" for b in bs)
    for r in out["accuracy"]:
        lines.append(f"| {r['mesh']} | {r['l_face']:g} | {r['dt']:g} | {r['substeps']} | {r['gamma_th_dt_at_qL2']:.3f} | "
                     f"{fmt(r['S'])} | {fmt(r['G_slow'])} / {fmt(r['G_fast'])} | {100 * r['worst_inside']:.1f} % |")

    # 2b. the frozen lFaces pooled at the same lFace-5 step: one lattice problem as far as bending goes
    lines += ["", "## Frozen, the three lFaces pooled at the same D dt / lFace^2 (step quoted at lFace 5)", "",
              "Mean +- standard error over every seed of lFace 5, 10 and 20. S, then Gamma (slow bins from the "
              "spec runs, fast from the dyn runs).", "",
              "| dt at lFace 5 | Gamma_th(2/L) dt | quantity | " + " | ".join(f"q lFace {a:g}-{b:g}" for a, b in BINS) + " |",
              "| --- | --- | --- | " + " | ".join("---" for _ in BINS) + " |"]
    pooled = defaultdict(lambda: defaultdict(list))
    for r in out["accuracy"]:
        if r["mesh"] != "frozen":
            continue
        for key, bins in (("S", r["S"]), ("G", r["G_slow"] + r["G_fast"])):
            for b in bins:
                pooled[(r["f"], key)][tuple(b["bin"])].extend(b["per_seed"])
    out["pooled"] = []
    for (f, key), by_bin in sorted(pooled.items()):
        cells = []
        for a, b in BINS:
            v = by_bin.get((a, b), [])
            if v:
                mean, sem = float(np.mean(v)), float(np.std(v, ddof=1) / math.sqrt(len(v))) if len(v) > 1 else float("nan")
                out["pooled"].append(dict(f=f, key=key, bin=[a, b], ratio=mean, err=sem, n=len(v)))
                cells.append(f"{mean:.3f} ± {sem:.3f}")
            else:
                cells.append("")
        lines.append(f"| {f:g} | {gamma_theory(2.0 / 5.0, 5.0) * f:.3f} | {'S' if key == 'S' else 'Γ'} | " + " | ".join(cells) + " |")

    # 3. mesh: kappa_app of the coarse references against lFace 5 at the same q
    ref5 = summaries("spec", 5.0, 0.001)
    if ref5:
        q5 = ref5[0]["q"]
        k5 = np.mean([KBT / (z["area"] * z["S"] * z["q"] ** 4) for z in ref5], axis=0)
        lines += ["", "## Mesh: kappa_app against lFace 5 at the same q (finest step of each lFace)", "",
                  "| lFace | q lFace bin | kappa_app / lFace 5 | modes |", "| --- | --- | --- | --- |"]
        for l_face in (10.0, 20.0):
            zs = summaries("spec", l_face, 0.001)
            if not zs:
                continue
            q = zs[0]["q"]
            k = np.mean([KBT / (z["area"] * z["S"] * z["q"] ** 4) for z in zs], axis=0)
            ok = (q >= max(q5.min(), rsa.Q_COMPARE)) & (q <= q5.max())
            # lFace 5 interpolated in log q onto the coarse mesh's modes (shell means first)
            order = np.argsort(q5)
            k5_at = np.interp(np.log(q[ok]), np.log(q5[order]), k5[order])
            for a, b in BINS:
                m = (q[ok] * l_face >= a) & (q[ok] * l_face < b)
                if m.sum():
                    ratio = float(np.mean(k[ok][m] / k5_at[m]))
                    out["mesh"].append(dict(l_face=l_face, bin=[a, b], ratio=ratio, modes=int(m.sum())))
                    lines.append(f"| {l_face:g} | {a:g}-{b:g} | {ratio:.3f} | {int(m.sum())} |")

    (root / "report.md").write_text("\n".join(lines) + "\n")
    (root / "report.json").write_text(json.dumps(out, indent=1))
    print("\n".join(lines))
    if figures:
        figure(out, Path(figures).expanduser())


def figure(out, directory):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from resolution_limit_figures import style, BLUE, ORANGE, INK, INK2, MUTED
    style()
    colors = {5.0: BLUE, 10.0: ORANGE, 20.0: "#1c9e77"}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, key, title in ((axes[0], "S", "static spectrum S(q)"), (axes[1], "G", "relaxation rate Γ(q)")):
        x = np.geomspace(0.003, 3, 50)
        # Euler-Maruyama on a mode of rate Gamma: variance 1/(1 - Gamma dt/2), lag-one decay 1 - Gamma dt
        theory = 1 / (1 - np.minimum(x, 1.9) / 2) - 1 if key == "S" else -np.log(1 - np.minimum(x, 0.99)) / x - 1
        ax.plot(x, 100 * theory, color=INK, lw=1, ls="--", label="Euler–Maruyama, linear mode (≈ Γ dt / 2)")
        for r in out["accuracy"]:
            if r["mesh"] != "frozen":
                continue
            bins = r["S"] if key == "S" else r["G_slow"] + r["G_fast"]
            for b in bins:
                if not math.isfinite(b["q_mid"]):
                    continue
                x0 = gamma_theory(b["q_mid"], r["l_face"]) * r["dt"]
                inside = b["bin"][1] <= 2.0
                ax.plot(x0, 100 * (b["ratio"] - 1), "o", color=colors[r["l_face"]], ms=5 if inside else 3.5,
                        alpha=1.0 if inside else 0.45, mec="none")
        for l_face, c in colors.items():
            ax.plot([], [], "o", color=c, label=f"lFace {l_face:g}")
        ax.axhspan(-5, 5, color=BLUE, alpha=0.07, lw=0)
        ax.axvline(0.1, color=MUTED, lw=1, ls=":")
        ax.set_xscale("log")
        ax.set_xlim(0.003, 3)
        ax.set_ylim(-20, 80 if key == "S" else 40)
        ax.set_xlabel("Γ_th(q) · dt   (Γ_th = D a κ_c q⁴ / kT)")
        ax.set_ylabel("error against the finest step (%)")
        ax.set_title(title, color=INK, fontsize=10)
        ax.legend(fontsize=8, loc="upper left")
    fig.text(0.01, 0.01, "Filled: bins inside q·lFace ≤ 2 (the mesh-resolved window); faint: beyond it. "
             "Shaded: ±5 %. Every run substeps the walls (split noise, dt/n ≤ the plain step's bound).",
             fontsize=7.5, color=INK2)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    directory.mkdir(parents=True, exist_ok=True)
    fig.savefig(directory / "substep_resolution.png")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("reduce")
    r.add_argument("dirs", nargs="+")
    p = sub.add_parser("report")
    p.add_argument("root")
    p.add_argument("--figures", default=None)
    p.add_argument("--fluid-ref", default=str(Path.home() / "slimed_runs/resolution/dt0.0005_nu0.5"))
    a = parser.parse_args()
    if a.cmd == "reduce":
        reduce(a.dirs)
    else:
        report(a.root, a.figures, a.fluid_ref)


if __name__ == "__main__":
    main()
