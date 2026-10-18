#!/usr/bin/env python
"""
Reduce the runs of analysis/resolution_sweep_runs.py and compare them.

    python analysis/resolution_sweep_analysis.py reduce <run dir> [<run dir> ...]
    python analysis/resolution_sweep_analysis.py report <root> [--ref-dt DIR] [--ref-mesh DIR]

`reduce` reads one run directory and writes `resolution_summary.npz` beside it:
whether and when it diverged, the height spectrum of the limit surface on the
tile (the route of analysis/flat_tension_matrix.ipynb: Loop limit surface,
resampled on the cartesian Monge grid scaled with the box, eight time blocks),
the kappa_c / sigma fit over |q| d <= pi/2, each mode's relaxation rate, the
energies, the box's tension and scale, and for a fluid run the flip rate and
the in-plane MSD. It is the expensive step and runs where the runs are.

`report` gathers the summaries under <root> and writes figures and stats.md to
<root>/analysis. The dt sweep is compared with --ref-dt (dt 0.0005) and the mesh
sweep with --ref-mesh (lFace 3.33), by default the runs of those names.

Every run is analysed over the same window, BURN to END us, so the error bars of
runs of different lengths are comparable.
"""
import argparse
import json
import os
import re
import sys
import time
from pathlib import Path

import numpy as np

if not hasattr(np, "trapezoid"):   # numpy < 2, as on the cluster's anaconda module
    np.trapezoid = np.trapz

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fluidity as fl                  # noqa: E402
import membrane_fluid_surface as mfs   # noqa: E402
import membrane_resample as mr         # noqa: E402
import membrane_spectrum as ms         # noqa: E402

KBT = 4.17
NBLOCK = 8
LEVELS = 2
QDX_MAX = np.pi / 2      # the continuum window of the kappa_c fit: four cells per wavelength
BURN = 100.0             # us discarded before every average
END = 400.0              # us, the end of the window (the 600 us references are cut here)
FRAME_SLACK = 0.25        # us: the last frame of a finished run can fall up to one frame short of END
Q_COMPARE = 0.12         # 1/nm: modes at and above this relax within a few us; the slower ones are under-sampled
BLOWUP = 1e7             # pN.nm: |E_total| above this is a divergence (the sheet holds ~1e3 to 1e4)


def read_params(path):
    out = {}
    for line in open(path):
        if "=" in line and not line.strip().startswith("#"):
            k, v = line.split("=", 1)
            out[k.strip()] = v.split("#")[0].strip()
    return out


def divergence(directory):
    """The iteration run_dynamics_flat stopped at for a divergence, or None."""
    log = Path(directory) / "run.log"
    if not log.exists():
        return None
    m = re.search(r" at iteration (\d+); the run has diverged", log.read_text(errors="replace"))
    return int(m.group(1)) if m else None


def unique_modes(lat):
    keep = np.zeros((lat.ny, lat.nx), bool)
    seen = set()
    for iy in range(lat.ny):
        for ix in range(lat.nx):
            if (iy, ix) == (0, 0) or ((-iy) % lat.ny, (-ix) % lat.nx) in seen:
                continue
            seen.add((iy, ix))
            keep[iy, ix] = True
    return keep


def block_edges(n, nblock=NBLOCK):
    return np.linspace(0, n, nblock + 1).astype(int)


def block_mean(y, nblock=NBLOCK):
    e = block_edges(len(y), nblock)
    b = np.array([np.nanmean(y[a:c]) for a, c in zip(e[:-1], e[1:])])
    return float(b.mean()), float(b.std(ddof=1) / np.sqrt(nblock))


def reduce_run(directory):
    t0 = time.perf_counter()
    d = Path(directory)
    p = read_params(d / "input.params")
    dt, interval = float(p["timeStep"]), int(p["meshpointOutputInterval"])
    l_face = float(p["lFace"])
    flips = p.get("edgeFlipEnabled", "false") == "true"
    out = dict(name=d.name, dt=dt, l_face=l_face, flips=flips,
               rate=float(p.get("edgeFlipAttemptRate", 0)) if flips else 0.0,
               kc_in=float(p["kcMembraneBending"]), D=float(p.get("diffConst", 1.0)),
               max_iterations=int(p["maxIterations"]), box=p.get("elasticBoxEnabled", "false") == "true",
               overrides=json.dumps({k: p[k] for k in ("creaseWallEnabled", "triangleShapeEnabled", "edgeSpringEnabled",
                                                       "inPlaneDynamicsEnabled", "elasticBoxEnabled") if k in p}))
    div = divergence(d)

    # energies, one line per step (header, then values)
    energy = d / "EnergyForce.csv"
    if energy.exists():
        e = np.genfromtxt(energy, delimiter=",", skip_header=1, invalid_raise=False)
        e = e if e.ndim == 2 else np.empty((0, 8))
        # The run stops itself only on a non-finite energy; one that blows up to a finite
        # 1e200 carries on to maxIterations. Either is a divergence here.
        blown = np.nonzero(~np.isfinite(e[:, 6]) | (np.abs(e[:, 6]) > BLOWUP))[0]
        if len(blown) and (div is None or blown[0] < div):
            div = int(blown[0])
            out["blowup_without_nan"] = div is not None and divergence(d) is None
        e = e[: blown[0]] if len(blown) else e
        t_e = np.arange(len(e)) * dt
        win = (t_e >= BURN) & (t_e <= END)
        out["steps_logged"] = len(e)
        for col, key in ((0, "E_curvature"), (1, "E_area"), (2, "E_mesh"), (6, "E_total")):
            out[key] = block_mean(e[win, col]) if win.sum() > NBLOCK * 2 else (np.nan, np.nan)
    out["diverged_iteration"] = -1 if div is None else div
    out["diverged_time"] = np.nan if div is None else div * dt

    run = fl.load_run(str(d), "input", interval)
    t = run.frame_iterations * dt
    out["frames"] = len(t)
    out["time_end"] = float(t[-1])
    t_stop = END + 1e-9 if div is None else min(END + 1e-9, div * dt)
    frames = np.where((t >= BURN) & (t <= t_stop))[0]
    out["analysed_frames"] = len(frames)
    if len(frames) < 50:
        np.savez(d / "resolution_summary.npz", **{k: np.asarray(v) for k, v in out.items()})
        print(f"{d.name}: {len(frames)} frames in the window, diverged at {div}; only the survival recorded")
        return

    lat = ms.lattice_from_params(str(d / "input.params"))
    scale = np.ones(len(t))
    box = d / "inputbox.csv"
    if box.exists():
        b = np.genfromtxt(box, delimiter=",", names=True)
        scale = np.interp(run.frame_iterations, b["iteration"], b["scaleX"])
        tension = np.interp(run.frame_iterations, b["iteration"], b["tension_pN_nm"])
        out["box_tension"] = block_mean(tension[frames])
        out["box_scale"] = float(scale[frames].mean())
        for col in ("tensionBending_pN_nm", "tensionAreaConstraint_pN_nm", "tensionMeshTerm_pN_nm", "tensionVertexGas_pN_nm"):
            out["box_" + col.replace("_pN_nm", "")] = float(np.interp(run.frame_iterations, b["iteration"], b[col])[frames].mean())

    # the limit surface on the Monge grid scaled with the box
    net0 = mr.reference_net(str(d / "meshpointinput.csv"), lat)
    X0, Y0 = mr.LimitSurface(lat, net0).monge_grid("cartesian")
    _, _, Q2D = mr.cartesian_q_grid(lat)
    keys = sorted(run.face_frames)
    cache, h = {}, np.empty((len(frames),) + X0.shape)
    for j, k in enumerate(frames):
        i = np.searchsorted(keys, run.frame_iterations[k], side="right") - 1
        key = keys[max(i, 0)]
        if key not in cache:
            cache.clear()   # connectivities are used in order; keep one
            cache[key] = mfs.SubdividedSurface(run.face_frames[key], run.coords.shape[1], LEVELS)
        h[j] = cache[key].interpolate(run.coords[k], X0 * scale[k], Y0 * scale[k])
    bad = ~np.isfinite(h).all(axis=(1, 2))
    if bad.any():
        print(f"{d.name}: {bad.sum()} frames with a grid point outside the mesh, dropped")
        h, frames = h[~bad], frames[~bad]
    h -= h.mean(axis=(1, 2), keepdims=True)
    out["h_rms"] = float(h.std())
    out["h_rms_blocks"] = block_mean(h.std(axis=(1, 2)))
    hf = np.fft.fft2(h, axes=(-2, -1))

    keep = unique_modes(lat)
    s_mean = float(scale[frames].mean())
    q = (Q2D / s_mean)[keep]
    area = lat.area * s_mean ** 2
    order = np.argsort(q)
    q = q[order]
    e = block_edges(len(hf))
    S_blocks = np.array([(np.mean(np.abs(hf[a:c]) ** 2, axis=0) / (lat.nx * lat.ny) ** 2)[keep][order]
                         for a, c in zip(e[:-1], e[1:])])
    win = q * lat.d_face_x * s_mean <= QDX_MAX
    fits = np.array([ms.fit_kc_sigma(q[win], Sb[win], area, KBT) for Sb in S_blocks])
    out.update(q=q, S_blocks=S_blocks, S=S_blocks.mean(axis=0), S_err=S_blocks.std(axis=0, ddof=1) / np.sqrt(NBLOCK),
               area=area, d_face=lat.d_face_x * s_mean, nx=lat.nx, ny=lat.ny,
               kc=(float(fits[:, 0].mean()), float(fits[:, 0].std(ddof=1) / np.sqrt(NBLOCK))),
               sigma=(float(fits[:, 1].mean()), float(fits[:, 1].std(ddof=1) / np.sqrt(NBLOCK))),
               slope=float(np.polyfit(np.log10(q[win]), np.log10(S_blocks.mean(axis=0)[win]), 1)[0]))

    # relaxation rates, mode by mode
    dt_frame = dt * interval
    acf = ms.mode_autocorrelation(hf)
    G = np.full(acf.shape[1:], np.nan)
    for iy, ix in zip(*np.nonzero(keep)):
        G[iy, ix] = ms.relaxation_rate(acf[:, iy, ix], dt_frame, n_frames=len(hf))
    out["gamma"] = G[keep][order]
    out["dt_frame"] = dt_frame

    # the fluid mesh
    if run.flip_accepted is not None and len(run.flip_accepted):
        it = run.flip_iteration * dt
        w = (it >= BURN) & (it <= END)
        out["flip_acceptance"] = float(run.flip_accepted[w].mean()) if w.any() else np.nan
        out["flips_per_us"] = float(run.flip_accepted[w].sum() / (END - BURN)) if w.any() else np.nan
    sub = fl.Run(directory=str(d), coords=run.coords[frames], frame_iterations=run.frame_iterations[frames],
                 free_vertices=run.free_vertices)
    lags = np.unique(np.geomspace(1, len(frames) // 4, 16).astype(int))
    tau, msd = fl.in_plane_msd(sub, lags)
    out["msd_t"], out["msd"] = tau * dt, msd

    np.savez(d / "resolution_summary.npz", **{k: np.asarray(v) for k, v in out.items()})
    print(f"{d.name}: {len(frames)} frames, kc {out['kc'][0]:.1f} +- {out['kc'][1]:.1f}, sigma {out['sigma'][0]:+.3f}, "
          f"h_rms {out['h_rms']:.3f}, {len(cache)} connectivity in cache, {time.perf_counter() - t0:.0f} s")


def load(root):
    runs = {}
    for f in sorted(Path(root).glob("*/resolution_summary.npz")):
        z = np.load(f, allow_pickle=False)
        runs[f.parent.name] = {k: (z[k].item() if z[k].ndim == 0 else z[k]) for k in z.files}
    return runs


def shells(q, y, err, edges):
    """Mean of y over the modes in each q shell, its error, and the shell's mean q."""
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (q >= a) & (q < b) & np.isfinite(y)
        if m.sum() >= 1:
            out.append((np.exp(np.log(q[m]).mean()), y[m].mean(), np.sqrt((err[m] ** 2).sum()) / m.sum(), m.sum()))
    return np.array(out).reshape(-1, 4)


def ref_curve(ref, q):
    """The reference spectrum at q, by log-log interpolation of its shell means."""
    edges = np.geomspace(ref["q"].min() * 0.999, ref["q"].max() * 1.001, 15)
    sh = shells(ref["q"], np.log(ref["S"]), ref["S_err"] / ref["S"], edges)
    return np.exp(np.interp(np.log(q), np.log(sh[:, 0]), sh[:, 1])), np.interp(np.log(q), np.log(sh[:, 0]), sh[:, 2])


def plain_ticks(ax, x=None, y=None):
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    for axis, ticks in ((ax.xaxis, x), (ax.yaxis, y)):
        if ticks:
            axis.set_major_locator(FixedLocator(ticks))
            axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            axis.set_minor_locator(NullLocator())


def survived(r):
    return r["diverged_iteration"] < 0 and r["time_end"] >= END - FRAME_SLACK


def report(root, ref_dt_name, ref_mesh_name, extra_refs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    root = Path(root)
    runs = load(root)
    for path in extra_refs:   # the lFace 5, dt 0.001 point lives in flat_tension
        z = load(Path(path).parent)
        name = Path(path).name
        if name in z:
            runs[f"ref1e-3_{name.split('_')[-1]}"] = z[name]
    out = root / "analysis"
    out.mkdir(exist_ok=True)
    plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 150, "font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "legend.frameon": False})
    fluid_of = lambda r: "nu0.5" if r["flips"] else "frozen"
    lines = []

    # ---- 1. survival -----------------------------------------------------------------
    # Divergence comes at a random time, so the seeds of one setting are pooled into a
    # rate: divergences per us of simulated time. With none, the 95 % upper limit 3/T.
    def variant(name):
        if name.startswith("diag_"):
            return name.split("_")[1]
        return "zonly" if name.startswith("zonly_") else ""

    groups = {}
    for name, r in runs.items():
        if r["diverged_iteration"] < 0 and r["time_end"] < END - FRAME_SLACK:
            continue   # still running
        key = (variant(name), fluid_of(r), float(r["l_face"]), float(r["dt"]))
        g = groups.setdefault(key, dict(n=0, div=0, T=0.0, times=[], names=[]))
        g["n"] += 1
        g["names"].append(name)
        if r["diverged_iteration"] >= 0:
            g["div"] += 1
            g["T"] += r["diverged_time"]
            g["times"].append(r["diverged_time"])
        else:
            g["T"] += min(r["time_end"], 600.0)
    rows = ["| variant | flips | lFace (nm) | dt (us) | runs | diverged | simulated (us) | divergences per ms | diverged at (us) |",
            "|---|---|---|---|---|---|---|---|---|"]
    for key in sorted(groups):
        g = groups[key]
        rate = (f"{1e3 * g['div'] / g['T']:.3g}" if g["div"] else f"< {1e3 * 3.0 / g['T']:.2g}")
        rows.append(f"| {key[0] or 'default'} | {key[1]} | {key[2]:.3g} | {key[3]:g} | {g['n']} | {g['div']} | {g['T']:.0f} | {rate} | "
                    f"{', '.join(f'{t:.3g}' for t in sorted(g['times']))} |")
    lines += ["## Survival, seeds pooled", "", *rows, ""]

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, fluid in zip(axes, ("frozen", "nu0.5")):
        for (var, fl_, l, dt), g in groups.items():
            if fl_ != fluid or var not in ("", "zonly"):
                continue
            x = l * (1.04 if var == "zonly" else 1.0)
            mk = "^" if var == "zonly" else "o"
            if g["div"] == 0:
                ax.scatter(x, dt, marker=mk, color="C2", s=18 + 10 * g["n"], zorder=3)
            else:
                mean_t = g["T"] / g["div"]   # the mean time to divergence
                c = plt.get_cmap("magma")(np.clip(np.log10(mean_t) / 3.0, 0, 0.9))
                ax.scatter(x, dt, marker=mk, color=c, edgecolor="k", lw=0.3, s=18 + 10 * g["n"], zorder=3)
                if g["div"] < g["n"]:
                    ax.scatter(x, dt, marker=mk, facecolor="none", edgecolor="C2", lw=1.2, s=40 + 10 * g["n"], zorder=4)
        ll = np.geomspace(3, 22, 50)
        ax.plot(ll, 0.0025 * (ll / 5) ** 2, "k:", lw=0.8, label=r"$\propto l^2$")
        ax.plot(ll, 0.0025 * (ll / 5), "k--", lw=0.8, label=r"$\propto l$")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("lFace (nm)")
        plain_ticks(ax, x=(3.33, 5, 10, 20), y=(0.001, 0.002, 0.005, 0.01, 0.02, 0.05))
        ax.set_title({"frozen": "frozen mesh", "nu0.5": "fluid mesh, ν = 0.5 /edge/µs"}[fluid], fontsize=9)
    axes[0].set_ylabel("timeStep (µs)")
    axes[0].scatter([], [], marker="^", color="0.5", label="height-only dynamics")
    sm = plt.cm.ScalarMappable(cmap="magma", norm=matplotlib.colors.LogNorm(1, 1000))
    fig.colorbar(sm, ax=axes, shrink=0.8, label="mean time to divergence (µs)")
    axes[0].scatter([], [], color="C2", label="no divergence (size: seeds)")
    axes[0].legend(fontsize=7, loc="upper left")
    fig.savefig(out / "stability_map.png", bbox_inches="tight"); plt.close(fig)

    # ---- 2. the time step ------------------------------------------------------------
    def series(prefix, fluid):
        sel = [(r["dt"], n, r) for n, r in runs.items() if (n.startswith(prefix) or n.startswith("stab_dt") or n.startswith("ref1e-3"))
               and fluid_of(r) == fluid and "S" in r and r["l_face"] == 5.0 and survived(r)]
        return [(n, r) for _, n, r in sorted(sel, key=lambda t: t[0])]

    stat_rows = ["| run | lFace | dt (us) | κ_c (pN·nm) | σ (pN/nm) | slope | h r.m.s. (nm) | E_bend/kT | box tension (pN/nm) | flips/µs | acceptance | MSD(10 µs) nm² |",
                 "|---|---|---|---|---|---|---|---|---|---|---|---|"]

    def stat_row(n, r):
        msd10 = np.interp(10.0, r["msd_t"], r["msd"]) if "msd" in r else np.nan
        bt = r.get("box_tension", (np.nan, np.nan))
        ec = r.get("E_curvature", (np.nan, np.nan))
        return (f"| {n} | {r['l_face']:.3g} | {r['dt']:g} | {r['kc'][0]:.1f} ± {r['kc'][1]:.1f} | {r['sigma'][0]:+.3f} ± {r['sigma'][1]:.3f} | "
                f"{r['slope']:.2f} | {r['h_rms']:.3f} | {ec[0] / KBT:.1f} ± {ec[1] / KBT:.1f} | {bt[0]:+.3f} ± {bt[1]:.3f} | "
                f"{r.get('flips_per_us', np.nan):.0f} | {r.get('flip_acceptance', np.nan):.3f} | {msd10:.1f} |")

    for fluid in ("frozen", "nu0.5"):
        ser = series("dt", fluid)
        ref = runs.get(ref_dt_name.replace("FLUID", fluid))
        if ref is None or "S" not in ref or len(ser) < 2:
            continue
        cmap = plt.get_cmap("viridis")
        fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
        edges = np.geomspace(ref["q"].min() * 0.999, ref["q"].max() * 1.001, 12)
        dts = [r["dt"] for _, r in ser]
        for n, r in ser:
            c = cmap(np.log(r["dt"] / min(dts)) / max(1e-9, np.log(max(dts) / min(dts))))
            ratio = r["S"] / ref["S"]
            err = ratio * np.hypot(r["S_err"] / r["S"], ref["S_err"] / ref["S"])
            sh = shells(r["q"], ratio, err, edges)
            axes[0].errorbar(sh[:, 0], sh[:, 1], sh[:, 2], fmt="o-", ms=3, lw=1, color=c, label=f"dt {r['dt']:g}")
            g = r["gamma"] / ref["gamma"]
            ok = np.isfinite(g)
            if ok.sum() > 2:
                sh = shells(r["q"][ok], g[ok], np.zeros(ok.sum()), edges)
                axes[1].plot(sh[:, 0], sh[:, 1], "o-", ms=3, lw=1, color=c)
            if "msd" in r:
                axes[2].loglog(r["msd_t"], r["msd"], "-", color=c, lw=1)
        axes[0].axhline(1, color="k", lw=0.6)
        axes[0].set_xscale("log"); axes[0].set_xlabel("q (1/nm)"); axes[0].set_ylabel(f"⟨|h_q|²⟩ / reference (dt {ref['dt']:g})")
        axes[0].set_title(f"{fluid}: static spectrum against dt {ref['dt']:g}", fontsize=9)
        plain_ticks(axes[0], x=(0.07, 0.1, 0.2, 0.3, 0.5, 0.7)); plain_ticks(axes[1], x=(0.12, 0.15, 0.2, 0.25))
        axes[0].legend(fontsize=7)
        axes[1].axhline(1, color="k", lw=0.6); axes[1].set_xscale("log"); axes[1].set_xlabel("q (1/nm)")
        axes[1].set_ylabel("Γ(q) / reference"); axes[1].set_title("relaxation rate", fontsize=9)
        axes[2].set_xlabel("lag (µs)"); axes[2].set_ylabel("in-plane MSD (nm²)"); axes[2].set_title("vertex in-plane MSD", fontsize=9)
        fig.tight_layout(); fig.savefig(out / f"dt_{fluid}.png"); plt.close(fig)
        stat_rows += [stat_row(n, r) for n, r in ser]

    # ---- 3. the mesh size ------------------------------------------------------------
    for fluid in ("frozen", "nu0.5"):
        sel = sorted([(r["l_face"], n, r) for n, r in runs.items()
                      if (n.startswith("mesh") or n.startswith("ref1e-3")) and fluid_of(r) == fluid and "S" in r and survived(r)],
                     key=lambda t: t[0])
        ref = runs.get(ref_mesh_name.replace("FLUID", fluid))
        if ref is None or "S" not in ref or len(sel) < 2:
            continue
        # Compared at the same |q|, on the modes that equilibrate inside the window
        # (q >= Q_COMPARE; the slower ones are under-sampled and their block errors
        # too small), as the bending rigidity each mode reads with sigma = 0.
        fig, axes = plt.subplots(1, 3, figsize=(13, 3.7))
        cmap = plt.get_cmap("plasma")
        ls_ = [l for l, _, _ in sel]
        qedges = np.geomspace(Q_COMPARE, 2.0, 14)
        kref = shells(ref["q"], KBT / (ref["area"] * ref["S"] * ref["q"] ** 4), np.zeros(len(ref["q"])), qedges)
        for l, n, r in sel:
            c = cmap(0.85 * np.log(l / min(ls_)) / np.log(max(ls_) / min(ls_)))
            k_app = KBT / (r["area"] * r["S"] * r["q"] ** 4)
            k_err = k_app * r["S_err"] / r["S"]
            sh = shells(r["q"], k_app, k_err, qedges)
            axes[0].errorbar(sh[:, 0], sh[:, 1], sh[:, 2], fmt="o-", ms=3, lw=1, color=c, label=f"l {l:.3g}")
            if r is not ref:
                kr = np.interp(np.log(sh[:, 0]), np.log(kref[:, 0]), kref[:, 1])
                axes[1].plot(sh[:, 0] * r["d_face"], sh[:, 1] / kr, "o-", ms=3, lw=1, color=c)
            a_vertex = r["area"] / (r["nx"] * r["ny"])
            ok = np.isfinite(r["gamma"])
            if ok.sum() > 1:
                sh = shells(r["q"][ok], r["gamma"][ok] / a_vertex, np.zeros(ok.sum()), np.geomspace(0.05, 0.5, 9))
                axes[2].plot(sh[:, 0], sh[:, 1], "o-", ms=3, lw=1, color=c)
        qq = np.geomspace(0.05, 0.5, 20)
        axes[2].plot(qq, ref["kc_in"] * qq ** 4 / KBT * ref["D"], "k:", lw=0.8, label="D κ_in q⁴ / k_BT")
        axes[0].axhline(ref["kc_in"], color="k", lw=0.6)
        axes[0].set_xscale("log"); axes[0].set_ylim(50, 200); axes[0].set_xlabel("q (1/nm)")
        axes[0].set_ylabel("k_BT / (A ⟨|h_q|²⟩ q⁴)  (pN·nm)"); axes[0].legend(fontsize=7, ncol=2)
        axes[0].set_title(f"{fluid}: apparent κ_c per mode (σ = 0)", fontsize=9)
        for y in (1.05, 1.10):
            axes[1].axhline(y, color="0.7", lw=0.6, ls="--")
        axes[1].axhline(1, color="k", lw=0.6); axes[1].set_xscale("log"); axes[1].set_xlabel("q · lFace")
        axes[1].set_ylabel(f"apparent κ_c / finest mesh (l {ref['l_face']:.3g}) at the same q")
        axes[1].set_title("discretization error", fontsize=9)
        axes[2].set_xscale("log"); axes[2].set_yscale("log"); axes[2].set_xlabel("q (1/nm)")
        axes[2].set_ylabel("Γ(q) / a  (1/(µs·nm²)), a = area per vertex"); axes[2].legend(fontsize=7)
        axes[2].set_title("relaxation rate per area per vertex", fontsize=9)
        plain_ticks(axes[0], x=(0.15, 0.2, 0.3, 0.5, 1.0)); plain_ticks(axes[1], x=(0.5, 1, 2, 3, 5))
        plain_ticks(axes[2], x=(0.1, 0.2, 0.3))
        fig.tight_layout(); fig.savefig(out / f"mesh_{fluid}.png"); plt.close(fig)
        stat_rows += [stat_row(n, r) for _, n, r in sel]

    for n, r in sorted(runs.items()):
        if (n.startswith(("big_", "diag_", "cross_", "zonly_", "seed"))) and "S" in r and survived(r):
            stat_rows.append(stat_row(n, r))
    lines += ["## Statistics, window %g-%g us" % (BURN, END), "", *stat_rows, ""]
    (out / "stats.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("reduce")
    r.add_argument("dirs", nargs="+")
    p = sub.add_parser("report")
    p.add_argument("root")
    p.add_argument("--ref-dt", default="dt0.0005_FLUID", help="FLUID is replaced by frozen / nu0.5")
    p.add_argument("--ref-mesh", default="mesh3.33_FLUID")
    p.add_argument("--extra", nargs="*", default=[], help="run dirs of the lFace 5, dt 0.001 point (flat_tension/box_*)")
    args = parser.parse_args()
    if args.cmd == "reduce":
        for directory in args.dirs:
            try:
                reduce_run(directory)
            except Exception as exc:   # one bad run must not stop the rest
                print(f"{directory}: FAILED {type(exc).__name__}: {exc}")
    else:
        report(args.root, args.ref_dt, args.ref_mesh, args.extra)


if __name__ == "__main__":
    main()
