#!/usr/bin/env python
"""Compare a fluid membrane and a frozen one pushed by a rigid bead.

Two run directories written by ``membrane_dynamics`` with the bead on
(``docs/bead_push_benchmark.md``): one with edge flips on (fluid) and one
with the same tethered mesh and flips off (frozen), so that whatever differs
is the fluidity. Reads from each:

    input.params                   the run's parameters
    <prefix>face.csv               the triangles at setup
    <prefix>face_<iteration>.csv   the triangles after each connectivity change (fluid)
    meshpoint<prefix>.csv          the control net, one frame per line
    surfacepoint<prefix>.csv       the limit surface, one frame per line
    <prefix>vertex_type_begin.csv  which vertices are clamped
    <prefix>bead.csv               the bead per frame
    <prefix>EdgeFlips.csv          every flip attempt (fluid)
    EnergyForce.csv                the energy terms per iteration

and writes into ``--out``:

    bead_push.gif              both membranes with the bead, beside the force on the bead
    bead_push_summary.png      bead and membrane heights, wrapping, waist, force, mesh quality, flips, energy
    bead_push_snapshots.png    both membranes at four bead heights
    bead_push_profiles.png     radius-height cross-sections at the same four heights

and prints a summary per run.

Usage:

    ~/anaconda3/bin/python analysis/bead_push_analysis.py --fluid DIR --frozen DIR \
        [--out DIR] [--prefix input] [--frames 200] [--fps 20] [--no-gif]
"""

import argparse
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import animation  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402
from mpl_toolkits.mplot3d import Axes3D  # noqa: E402,F401

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fluidity  # noqa: E402  the run loader and the mesh-quality helpers

# ---------------------------------------------------------------------------
# Palette: one blue ramp for height, two categorical slots for the two runs,
# a neutral for the bead, neutral inks and surfaces.
# ---------------------------------------------------------------------------
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
             "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
             "#0d366b"]
HEIGHT_CMAP = LinearSegmentedColormap.from_list("slimed_blue", BLUE_RAMP)
RUN_COLOR = {"fluid": "#2a78d6", "frozen": "#eb6834"}
BEAD_COLOR = "#8f8c85"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e5e4e0"
SURFACE = "#fcfcfb"

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK_2, "axes.labelcolor": INK,
    "xtick.color": INK_2, "ytick.color": INK_2, "text.color": INK,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def read_params(path):
    params = {}
    with open(path) as handle:
        for line in handle:
            line = line.split("#", 1)[0].strip()
            if "=" in line:
                key, value = line.split("=", 1)
                params[key.strip()] = value.strip().rstrip(";").strip()
    return params


def read_frames(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            line = line.strip().rstrip(",")
            if line:
                rows.append(np.fromstring(line, sep=","))
    widths = {}
    for row in rows:
        widths[len(row)] = widths.get(len(row), 0) + 1
    width = max(widths, key=widths.get)
    rows = [row for row in rows if len(row) == width]
    return np.array(rows).reshape(len(rows), -1, 3)


class BeadRun:
    """One run: its frames, its bead, and the quantities the figures need."""

    def __init__(self, directory, label, prefix="input"):
        self.directory = directory
        self.label = label
        self.prefix = prefix
        p = read_params(os.path.join(directory, "input.params"))
        self.params = p
        self.dt = float(p.get("timeStep", 0.1))
        self.interval = int(p.get("meshpointOutputInterval", 1))
        self.radius = float(p.get("beadRadius", 10.0))
        self.flips = p.get("edgeFlipEnabled", "false") == "true"
        self.l_face = float(p.get("lFace", 5.0))
        self.tether = (float(p.get("edgeTetherMinRatio", 0.95)) * self.l_face,
                       float(p.get("edgeTetherMaxRatio", 1.75)) * self.l_face)
        self.altitude_floor = float(p.get("triangleShapeMinAltitudeRatio", 0.4)) * self.l_face
        self.vertex_area = (3.0 ** 0.5) / 2.0 * self.l_face ** 2

        self.run = fluidity.load_run(directory, prefix, self.interval)
        self.surface = read_frames(os.path.join(directory, f"surfacepoint{prefix}.csv"))
        n = min(len(self.run.coords), len(self.surface))
        self.n_frames = n
        self.control = self.run.coords[:n]
        self.surface = self.surface[:n]
        self.iterations = self.run.frame_iterations[:n]
        self.time = self.iterations * self.dt
        self.free = np.zeros(self.control.shape[1], dtype=bool)
        self.free[self.run.free_vertices] = True

        bead = np.genfromtxt(os.path.join(directory, f"{prefix}bead.csv"), delimiter=",", names=True)
        bead = np.atleast_1d(bead)
        by_iteration = {int(row["iteration"]): row for row in bead}
        rows = [by_iteration[int(it)] for it in self.iterations if int(it) in by_iteration]
        self.n_frames = n = min(n, len(rows))
        rows = rows[:n]
        self.control, self.surface = self.control[:n], self.surface[:n]
        self.iterations, self.time = self.iterations[:n], self.time[:n]
        self.bead_center = np.array([[r["centerX"], r["centerY"], r["centerZ"]] for r in rows])
        self.bead_energy = np.array([r["energy_pNnm"] for r in rows])
        self.contacts = np.array([r["contacts"] for r in rows])
        # The adhesion columns exist only in logs written with the well; a
        # bead that only pushes counts its contacts as the wrapped vertices.
        names = bead.dtype.names
        self.adhered = np.array([r["adhered"] for r in rows]) if "adhered" in names else self.contacts.copy()
        self.adhesion_energy = (np.array([r["adhesion_pNnm"] for r in rows]) if "adhesion_pNnm" in names
                                else np.zeros(len(rows)))
        self.adhesion = float(p.get("beadAdhesionEnergy", 0.0))
        self.penetration = np.array([r["maxPenetration_nm"] for r in rows])
        self.force_on_bead = np.array([[r["forceOnBeadX"], r["forceOnBeadY"], r["forceOnBeadZ"]] for r in rows])

        # EnergyForce.csv is written when the run ends; a run still going has
        # none, and the figures then leave the energy panel empty.
        energy_path = os.path.join(directory, "EnergyForce.csv")
        if os.path.exists(energy_path):
            energy = np.genfromtxt(energy_path, delimiter=",", skip_header=1)
            index = np.minimum(self.iterations + 1, len(energy) - 1)
            index[0] = 0
            self.energy_bending = energy[index, 0]
            self.energy_area = energy[index, 1]
            self.energy_mesh = energy[index, 2]   # tether + altitude floor + crease wall
            self.energy_total = energy[index, 6]
        else:
            self.energy_bending = self.energy_area = self.energy_mesh = self.energy_total = np.full(n, np.nan)

        self.base_level = float(np.median(self.surface[0, self.free, 2]))
        self.bead_height = self.bead_center[:, 2] - self.base_level
        # Wrapped fraction: the vertices on or at the bead, each standing for
        # one lattice cell of membrane, against the sphere's area.
        self.wrapping = self.adhered * self.vertex_area / (4.0 * np.pi * self.radius ** 2)
        self.apex = self.surface[:, :, 2][:, self.free].max(axis=1) - self.base_level
        self.waist = np.array([self.waist_radius(k) for k in range(n)])
        self.compute_mesh_quality()
        self.compute_flip_rates()

    # The bottleneck: the smallest distance from the bead's axis among the
    # surface points that lie behind the bead -- above the sheet, below the
    # lower half of the sphere.
    def waist_radius(self, k):
        cx, cy, cz = self.bead_center[k]
        pts = self.surface[k][self.free]
        low = self.base_level + 0.2 * self.radius
        high = cz - 0.5 * self.radius
        if high <= low:
            return np.nan
        slab = (pts[:, 2] > low) & (pts[:, 2] < high)
        if not slab.any():
            return np.nan
        rho = np.hypot(pts[slab, 0] - cx, pts[slab, 1] - cy)
        return float(rho.min())

    def compute_mesh_quality(self):
        n = self.n_frames
        self.min_altitude = np.full(n, np.nan)
        self.below_floor = np.full(n, np.nan)
        self.max_edge = np.full(n, np.nan)
        self.min_edge = np.full(n, np.nan)
        self.max_crease = np.full(n, np.nan)
        self.creases_over_90 = np.full(n, np.nan)
        for k in range(n):
            alt = fluidity.altitude_stats(self.run, k, floor=self.altitude_floor)
            edges = fluidity.edge_length_stats(self.run, k)
            crease = fluidity.crease_stats(self.run, k)
            if alt.get("n", 0):
                self.min_altitude[k] = alt["min"]
                self.below_floor[k] = alt["fraction_below"]
            if edges.get("n", 0):
                self.max_edge[k] = edges["max"]
                self.min_edge[k] = edges["min"]
            if crease.get("n", 0):
                self.max_crease[k] = crease["max"]
                self.creases_over_90[k] = crease["over_90"]

    def compute_flip_rates(self):
        """Accepted flips per microsecond, and how far from the bead's axis."""
        n = self.n_frames
        self.flips_per_us = np.zeros(n)
        self.flip_rho = np.array([])
        self.flip_time = np.array([])
        path = os.path.join(self.directory, f"{self.prefix}EdgeFlips.csv")
        if not (self.flips and os.path.exists(path)):
            return
        table = np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True))
        accepted = table[table["accepted"] > 0.5]
        if accepted.size == 0:
            return
        it = accepted["iteration"].astype(int)
        frame = np.clip(np.round(it / self.interval).astype(int), 0, n - 1)
        counts = np.bincount(frame, minlength=n)[:n]
        self.flips_per_us = counts / (self.interval * self.dt)
        v0 = accepted["v0"].astype(int)
        v1 = accepted["v1"].astype(int)
        mid = 0.5 * (self.control[frame, v0] + self.control[frame, v1])
        self.flip_rho = np.hypot(mid[:, 0] - self.bead_center[frame, 0], mid[:, 1] - self.bead_center[frame, 1])
        self.flip_time = it * self.dt

    def faces_at(self, k):
        target = self.iterations[k]
        usable = [it for it in sorted(self.run.face_frames) if it <= target]
        return self.run.face_frames[usable[-1] if usable else min(self.run.face_frames)]

    def frame_at_height(self, height):
        """The first frame with the bead centre at least this high above the sheet."""
        hit = np.nonzero(self.bead_height >= height)[0]
        return int(hit[0]) if len(hit) else self.n_frames - 1


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------
def draw_bead(ax, centre, radius):
    u, v = np.mgrid[0:2 * np.pi:28j, 0:np.pi:14j]
    ax.plot_surface(centre[0] + radius * np.sin(v) * np.cos(u), centre[1] + radius * np.sin(v) * np.sin(u),
                    centre[2] + radius * np.cos(v), color=BEAD_COLOR, alpha=0.55, linewidth=0, shade=True)


def draw_membrane(ax, run, k, zmax, xlim, zoom=1.35):
    ax.cla()
    pts = run.surface[k]
    faces = run.faces_at(k)
    surface = ax.plot_trisurf(pts[:, 0], pts[:, 1], pts[:, 2], triangles=faces, cmap=HEIGHT_CMAP,
                              linewidth=0.2, edgecolor=(0.0, 0.0, 0.0, 0.14), antialiased=True)
    surface.set_clim(run.base_level, run.base_level + zmax)
    draw_bead(ax, run.bead_center[k], run.radius)
    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(-xlim, xlim)
    zlo, zhi = run.base_level - 0.25 * zmax, run.base_level + 1.1 * zmax
    ax.set_zlim(zlo, zhi)
    ax.set_box_aspect((1.0, 1.0, (zhi - zlo) / (2.0 * xlim)), zoom=zoom)
    ax.view_init(elev=16, azim=-58)
    ax.set_facecolor(SURFACE)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        axis._axinfo["grid"]["color"] = GRID
        axis._axinfo["grid"]["linewidth"] = 0.5
    ax.zaxis.set_major_locator(MaxNLocator(4))
    ax.tick_params(labelsize=7, colors=INK_2, pad=-2)
    ax.set_xlabel("x (nm)", fontsize=8, color=INK_2, labelpad=-6)
    ax.set_ylabel("y (nm)", fontsize=8, color=INK_2, labelpad=-6)
    ax.set_zlabel("z (nm)", fontsize=8, color=INK_2, labelpad=-6)


def panel_title(run, k):
    mode = "fluid" if run.flips else "frozen"
    return f"{mode}, bead at {run.bead_height[k]:+.0f} nm"


def style_trace(ax):
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=8.5)


def shade_motion(ax, runs):
    """Grey while the bead is moving (its height changes), on the run that lasts longest."""
    run = max(runs, key=lambda r: r.time[-1])
    moving = np.nonzero(np.diff(run.bead_height) > 1e-9)[0]
    if len(moving):
        ax.axvspan(run.time[moving[0]], run.time[moving[-1] + 1], color="#f0efec", zorder=0, lw=0)


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
def write_summary_figure(runs, path):
    fig, axes = plt.subplots(4, 2, figsize=(11.0, 11.5), dpi=140, sharex=True,
                             gridspec_kw={"hspace": 0.22, "wspace": 0.28})
    a = axes.ravel()
    for ax in a:
        shade_motion(ax, runs)
        style_trace(ax)
    longest = max(runs, key=lambda r: r.time[-1])
    a[0].plot(longest.time, longest.bead_height, color=INK_2, lw=1.2, ls="--", label="bead centre")
    for r in runs:
        c = RUN_COLOR[r.label]
        a[0].plot(r.time, r.apex, color=c, lw=1.5, label=f"membrane apex, {r.label}")
        a[1].plot(r.time, r.wrapping, color=c, lw=1.5, label=r.label)
        a[2].plot(r.time, r.waist, color=c, lw=1.5, label=r.label)
        a[3].plot(r.time, r.force_on_bead[:, 2], color=c, lw=0.6, alpha=0.3)
        a[3].plot(r.time, running_mean(r.force_on_bead[:, 2], r, 2.0), color=c, lw=1.6,
                  label=f"{r.label} (2 \u00b5s mean)")
        a[4].plot(r.time, r.min_altitude, color=c, lw=1.2, label=r.label)
        a[5].plot(r.time, r.max_edge, color=c, lw=1.2, label=f"longest edge, {r.label}")
        a[5].plot(r.time, r.min_edge, color=c, lw=1.2, ls=":", label=f"shortest edge, {r.label}")
        if r.flips:
            a[6].plot(r.time, r.flips_per_us, color=c, lw=1.0, label=f"accepted flips, {r.label}")
        a[7].plot(r.time, r.energy_area, color=c, lw=1.3, label=f"area constraint, {r.label}")
        a[7].plot(r.time, r.energy_mesh, color=c, lw=1.3, ls=":", label=f"tether + shape + crease, {r.label}")
    a[2].axhline(runs[0].radius, color=INK_2, lw=0.9, ls="--")
    a[2].text(0.01, runs[0].radius, " bead radius", va="bottom", fontsize=8, color=INK_2,
              transform=a[2].get_yaxis_transform())
    a[4].axhline(runs[0].altitude_floor, color=INK_2, lw=0.9, ls="--")
    a[4].text(0.01, runs[0].altitude_floor, " altitude floor", va="bottom", fontsize=8, color=INK_2,
              transform=a[4].get_yaxis_transform())
    for wall in runs[0].tether:
        a[5].axhline(wall, color=INK_2, lw=0.9, ls="--")
    a[5].text(0.01, runs[0].tether[1], " tether walls", va="bottom", fontsize=8, color=INK_2,
              transform=a[5].get_yaxis_transform())
    labels = ["height above the sheet (nm)", "wrapped fraction of the bead",
              "waist radius behind the bead (nm)", "force on the bead, z (pN)",
              "smallest triangle altitude (nm)", "edge length (nm)",
              "accepted flips per µs", "energy (pN·nm)"]
    for ax, label in zip(a, labels):
        ax.set_ylabel(label, fontsize=9)
        ax.legend(loc="best", frameon=False, fontsize=7.5)
    a[7].legend(loc="upper left", frameon=False, fontsize=7.5)
    a[2].set_ylim(0, 1.6 * runs[0].radius)
    a[4].set_ylim(0, None)
    a[6].set_ylim(0, None)
    a[7].set_ylim(0, None)
    for ax in a[-2:]:
        ax.set_xlabel("time (µs)")
    fig.subplots_adjust(top=0.95)
    fig.suptitle("A rigid bead pushed through a fluid and a frozen membrane", x=0.5, y=0.985,
                 fontsize=12, color=INK)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def snapshot_frames(run):
    """Four frames: the bead centre at the sheet, one radius up, at its top height, and the end."""
    top = float(run.bead_height.max())
    return [run.frame_at_height(0.0), run.frame_at_height(run.radius), run.frame_at_height(top), run.n_frames - 1]


def write_snapshots(runs, path, zmax, xlim):
    fig = plt.figure(figsize=(15.5, 3.6 * len(runs) + 0.6), dpi=120)
    for i, run in enumerate(runs):
        frames = snapshot_frames(run)
        for j, k in enumerate(frames):
            ax = fig.add_subplot(len(runs), len(frames), i * len(frames) + j + 1, projection="3d")
            draw_membrane(ax, run, k, zmax, xlim, zoom=1.12)
            ax.set_title(f"{panel_title(run, k)}, t = {run.time[k]:.0f} µs", fontsize=9, color=INK, pad=1)
    fig.subplots_adjust(left=0.0, right=1.0, top=0.94, bottom=0.01, wspace=0.0, hspace=0.08)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def write_profiles(runs, path):
    fig, axes = plt.subplots(len(runs), 4, figsize=(14.0, 3.4 * len(runs)), dpi=140,
                             sharex=True, sharey=True)
    axes = np.atleast_2d(axes)
    for i, run in enumerate(runs):
        for j, k in enumerate(snapshot_frames(run)):
            ax = axes[i, j]
            cx, cy, cz = run.bead_center[k]
            pts = run.control[k][run.free]
            rho = np.hypot(pts[:, 0] - cx, pts[:, 1] - cy)
            ax.scatter(rho, pts[:, 2] - run.base_level, s=6, color=RUN_COLOR[run.label], alpha=0.7, lw=0)
            angle = np.linspace(-np.pi / 2, np.pi / 2, 80)
            ax.plot(run.radius * np.cos(angle), cz - run.base_level + run.radius * np.sin(angle),
                    color=BEAD_COLOR, lw=1.6)
            ax.axhline(0.0, color=GRID, lw=0.8)
            ax.set_title(f"{panel_title(run, k)}, t = {run.time[k]:.0f} µs", fontsize=9, color=INK)
            ax.set_aspect("equal")
            style_trace(ax)
            if i == len(runs) - 1:
                ax.set_xlabel("distance from the bead's axis (nm)", fontsize=9)
            if j == 0:
                ax.set_ylabel("height above the sheet (nm)", fontsize=9)
    fig.suptitle("Cross-sections: every free control point by its distance from the axis", fontsize=11,
                 color=INK, y=0.995)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def write_animation(runs, path, n_frames, fps, zmax, xlim):
    n = min(r.n_frames for r in runs)
    picked = np.unique(np.linspace(0, n - 1, min(n_frames, n)).round().astype(int))
    hold = int(round(0.8 * fps))
    sequence = [picked[0]] * hold + list(picked) + [picked[-1]] * (2 * hold)

    fig = plt.figure(figsize=(16.0, 5.6), dpi=100)
    ax3d = [fig.add_axes([0.005 + 0.295 * i, 0.02, 0.29, 0.9], projection="3d") for i in range(len(runs))]
    ax2d = fig.add_axes([0.745, 0.16, 0.245, 0.68])
    shade_motion(ax2d, runs)
    for r in runs:
        ax2d.plot(r.time, r.force_on_bead[:, 2], color=RUN_COLOR[r.label], lw=1.3, label=r.label)
    ax2d.set_xlim(0.0, max(r.time[-1] for r in runs))
    ax2d.set_xlabel("time (µs)")
    ax2d.set_ylabel("force on the bead, z (pN)")
    ax2d.legend(loc="upper left", frameon=False, fontsize=9)
    style_trace(ax2d)
    cursor = ax2d.axvline(runs[0].time[0], color=INK_2, lw=0.9)
    titles = [fig.text(0.02 + 0.295 * i, 0.95, "", fontsize=11, color=INK) for i in range(len(runs))]
    stamp = fig.text(0.745, 0.95, "", fontsize=11, color=INK)

    def update(k):
        for ax, run, title in zip(ax3d, runs, titles):
            draw_membrane(ax, run, k, zmax, xlim)
            title.set_text(panel_title(run, k))
        cursor.set_xdata([runs[0].time[k], runs[0].time[k]])
        stamp.set_text(f"t = {runs[0].time[k]:.1f} µs")
        return tuple(ax3d) + (cursor, stamp)

    anim = animation.FuncAnimation(fig, update, frames=sequence, interval=1000.0 / fps, blit=False)
    anim.save(path, writer=animation.PillowWriter(fps=fps), dpi=100)
    plt.close(fig)
    return len(sequence)


def running_mean(values, run, window_us):
    width = max(1, int(round(window_us / (run.interval * run.dt))))
    kernel = np.ones(width) / width
    padded = np.concatenate([np.full(width // 2, values[0]), values,
                             np.full(width - 1 - width // 2, values[-1])])
    return np.convolve(padded, kernel, mode="valid")


def summarize(run):
    def fmt(v):
        return f"{v:.4g}" if isinstance(v, (float, np.floating)) else str(v)
    waist = run.waist[np.isfinite(run.waist)]
    lines = {
        "label": run.label,
        "frames": run.n_frames,
        "duration_us": run.time[-1],
        "energy_finite_at_end": bool(np.isfinite(run.energy_total[-1])),
        "bead_height_max_nm": run.bead_height.max(),
        "wrapping_max": run.wrapping.max(),
        "wrapping_at_end": run.wrapping[-1],
        "adhered_at_end": int(run.adhered[-1]),
        "contacts_at_end": int(run.contacts[-1]),
        "waist_min_nm": waist.min() if waist.size else float("nan"),
        "waist_min_time_us": run.time[np.nanargmin(run.waist)] if waist.size else float("nan"),
        "waist_at_end_nm": run.waist[-1],
        "force_on_bead_z_max_pN": run.force_on_bead[:, 2].max(),
        "force_on_bead_z_at_end_pN": run.force_on_bead[-1, 2],
        "force_on_bead_z_last_20us_mean_pN": run.force_on_bead[run.time >= run.time[-1] - 20.0, 2].mean(),
        "penetration_max_nm": run.penetration.max(),
        "min_altitude_overall_nm": np.nanmin(run.min_altitude),
        "faces_below_floor_max": np.nanmax(run.below_floor),
        "max_edge_overall_nm": np.nanmax(run.max_edge),
        "min_edge_overall_nm": np.nanmin(run.min_edge),
        "max_crease_deg_overall": np.nanmax(run.max_crease),
        "creases_over_90_at_end": run.creases_over_90[-1],
        "accepted_flips_total": int(run.flips_per_us.sum() * run.interval * run.dt) if run.flips else 0,
        "accepted_flips_within_2R_fraction": (float((run.flip_rho < 2 * run.radius).mean())
                                              if run.flip_rho.size else float("nan")),
        "area_energy_at_end_pNnm": run.energy_area[-1],
        "mesh_energy_at_end_pNnm": run.energy_mesh[-1],
        "bending_energy_at_end_pNnm": run.energy_bending[-1],
    }
    for key, value in lines.items():
        print(f"{key} = {fmt(value)}")
    print()
    return lines


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fluid", required=True)
    parser.add_argument("--frozen", required=True)
    parser.add_argument("--out", default=None, help="output directory (default: the fluid run's)")
    parser.add_argument("--prefix", default="input")
    parser.add_argument("--frames", type=int, default=200)
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument("--no-gif", action="store_true")
    args = parser.parse_args(argv)

    runs = [BeadRun(args.fluid, "fluid", args.prefix), BeadRun(args.frozen, "frozen", args.prefix)]
    out = args.out or args.fluid
    os.makedirs(out, exist_ok=True)
    for r in runs:
        summarize(r)
    # The height range: the bead's rise, and never less than the bead itself.
    zmax = float(np.ceil(max(max(r.apex.max(), r.bead_height.max() + r.radius) for r in runs) * 1.05))
    zmax = max(zmax, 2.5 * runs[0].radius)
    xlim = float(np.abs(runs[0].control[0, :, :2]).max())

    write_summary_figure(runs, os.path.join(out, "bead_push_summary.png"))
    print("wrote", os.path.join(out, "bead_push_summary.png"))
    write_profiles(runs, os.path.join(out, "bead_push_profiles.png"))
    print("wrote", os.path.join(out, "bead_push_profiles.png"))
    write_snapshots(runs, os.path.join(out, "bead_push_snapshots.png"), zmax, xlim)
    print("wrote", os.path.join(out, "bead_push_snapshots.png"))
    if not args.no_gif:
        n = write_animation(runs, os.path.join(out, "bead_push.gif"), args.frames, args.fps, zmax, xlim)
        print("wrote", os.path.join(out, "bead_push.gif"), f"({n} frames at {args.fps} fps)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
