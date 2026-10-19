#!/usr/bin/env python
"""
Animate the last steps of a run that diverges.

    python analysis/blowup_animation.py <run dir> <out.gif> [--quiet-stride 20] [--last 16] [--label NAME]

The run must write a frame every step (meshpointOutputInterval = 1); rerun a
diverging setting that way, it costs seconds. Three panels: the whole sheet's
control net from above with the periodic tile outlined, a zoom on the vertex
whose step is largest just before the end, and the energies per step on a log
scale. The quiet part of the run is shown every --quiet-stride steps, the last
--last steps one by one and slower.
"""
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.collections import LineCollection

import fluidity as fl
import membrane_spectrum as ms

BLUE, ORANGE, RED = "#2a78d6", "#eb6834", "#d03b3b"
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"


def edges_of(faces):
    e = np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    return np.unique(e, axis=0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("run")
    p.add_argument("out")
    p.add_argument("--quiet-stride", type=int, default=20)
    p.add_argument("--last", type=int, default=16)
    p.add_argument("--zoom", type=float, default=14.0, help="half-width of the zoom, nm")
    p.add_argument("--label", default=None, help="the run's name in the title (default: the directory's)")
    args = p.parse_args()
    d = Path(args.run)
    run = fl.load_run(str(d), "input", 1)
    C = run.coords
    n = len(C)
    dt = float([l.split("=")[1].split("#")[0] for l in open(d / "input.params") if l.startswith("timeStep")][0])
    lat = ms.lattice_from_params(str(d / "input.params"))
    energy = np.genfromtxt(d / "EnergyForce.csv", delimiter=",", skip_header=1)[:n]
    types = np.genfromtxt(d / "inputvertex_type_begin.csv", delimiter=",")[:, 3].astype(int)
    keys = sorted(run.face_frames)

    def faces_at(k):
        i = np.searchsorted(keys, k, side="right") - 1
        return run.face_frames[keys[max(i, 0)]]

    # the vertex that moves most in the last finite step, and the seam edge it sits on
    step = np.linalg.norm(C[n - 2] - C[n - 3], axis=1)
    hot = int(np.argmax(step))
    centre = C[n - 4, hot, :2]
    tile_x, tile_y = lat.lx, lat.ly
    x0 = C[0, :, 0].min() + 3 * lat.d_face_x
    y0 = C[0, :, 1].min() + 3 * lat.d_face_y

    frames = list(range(0, n - args.last, args.quiet_stride)) + list(range(n - args.last, n))
    frames += [n - 1] * 6   # hold the end
    plt.rcParams.update({"font.size": 9, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                         "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(11, 5.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], height_ratios=[1.25, 1], wspace=0.22, hspace=0.35)
    ax_all, ax_zoom, ax_e = fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])
    lim = np.abs(C[0, :, :2]).max() + 3
    for ax in (ax_all,):
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
        ax.add_patch(plt.Rectangle((x0, y0), tile_x, tile_y, fill=False, ec=INK, lw=1.0, ls="--"))
        ax.text(x0, y0 + tile_y + 1.5, "periodic tile", color=INK2, fontsize=8)
        ax.set_xlabel("x (nm)"); ax.set_ylabel("y (nm)")
    ax_zoom.set_xlim(centre[0] - args.zoom, centre[0] + args.zoom)
    ax_zoom.set_ylim(centre[1] - args.zoom * 0.62, centre[1] + args.zoom * 0.62)
    ax_zoom.set_aspect("equal")
    ax_zoom.axhline(y0 + tile_y, color=INK, lw=1, ls="--")
    ax_zoom.text(centre[0] - args.zoom + 0.5, y0 + tile_y + 0.4, "seam (top edge of the tile)", color=INK2, fontsize=7.5)
    t = np.arange(len(energy))
    ax_e.semilogy(t, energy[:, 6], color=INK2, lw=1.2, label="total")
    ax_e.semilogy(t, energy[:, 2], color=ORANGE, lw=1.2, label="mesh terms (tether, altitude floor, crease wall)")
    ax_e.semilogy(t, energy[:, 0], color=BLUE, lw=1.2, label="bending")
    ax_e.set_xlabel("step"); ax_e.set_ylabel("energy (pN·nm)")
    ax_e.set_ylim(100, max(1e6, np.nanmax(energy[:, 6]) * 2))
    ax_e.legend(fontsize=7, loc="upper left", frameon=False)
    cursor = ax_e.axvline(0, color=RED, lw=1)

    all_lines = LineCollection([], colors=BLUE, linewidths=0.4)
    zoom_lines = LineCollection([], colors=BLUE, linewidths=1.0)
    hot_lines = LineCollection([], colors=RED, linewidths=2.0)
    ax_all.add_collection(all_lines); ax_zoom.add_collection(zoom_lines); ax_zoom.add_collection(hot_lines)
    pts_img = ax_zoom.scatter([], [], s=14, color=MUTED, zorder=3, label="periodic image")
    pts_free = ax_zoom.scatter([], [], s=14, color=BLUE, zorder=3, label="tile vertex")
    pts_hot = ax_zoom.scatter([], [], s=40, color=RED, zorder=4, label="largest last step")
    ax_zoom.legend(fontsize=7, loc="lower right", frameon=False)
    title = fig.suptitle("", fontsize=11, x=0.02, ha="left", color=INK)

    def draw(k):
        E = edges_of(faces_at(k))
        P = C[k]
        segs = P[E][:, :, :2]
        all_lines.set_segments(segs)
        zoom_lines.set_segments(segs)
        hot_e = E[np.any(E == hot, axis=1)]
        hot_lines.set_segments(P[hot_e][:, :, :2])
        img = types == 3
        pts_img.set_offsets(P[img, :2]); pts_free.set_offsets(P[~img, :2]); pts_hot.set_offsets(P[[hot], :2])
        cursor.set_xdata([k, k])
        last = k >= n - args.last
        # the whole run while it is quiet, then the last steps
        ax_e.set_xlim((n - 3 * args.last, n) if last else (0, n))
        title.set_text(f"{args.label or d.name}: step {k} of {n - 1}  (t = {k * dt:.3f} µs, dt = {dt:g} µs)"
                       + ("   — the last steps, one per frame" if last else ""))
        ax_zoom.set_title(f"vertex {hot} and its edges: max step this frame "
                          f"{np.linalg.norm(C[k] - C[max(k - 1, 0)], axis=1).max():.2g} nm", fontsize=8.5, color=INK2)
        return all_lines, zoom_lines, hot_lines, pts_img, pts_free, pts_hot, cursor, title

    durations = [120 if f < n - args.last else 700 for f in frames]
    anim = animation.FuncAnimation(fig, draw, frames=frames, blit=False)
    writer = animation.PillowWriter(fps=6)
    anim.save(args.out, writer=writer, dpi=75)
    # PillowWriter has one frame rate; re-time the frames so the end is slow
    from PIL import Image
    im = Image.open(args.out)
    out_frames = []
    try:
        while True:
            out_frames.append(im.copy()); im.seek(im.tell() + 1)
    except EOFError:
        pass
    out_frames[0].save(args.out, save_all=True, append_images=out_frames[1:], duration=durations[: len(out_frames)],
                       loop=0, optimize=True)
    print(f"{args.out}: {len(frames)} frames, hot vertex {hot} (type {types[hot]}), blow-up from step {n - args.last}")


if __name__ == "__main__":
    main()
