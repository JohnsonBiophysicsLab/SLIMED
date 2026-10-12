#!/usr/bin/env python
"""Animate a SLIMED push-and-release run.

The membrane's limit surface over time, beside the height its apex reaches,
with the window the external force was on marked -- for a run written by
``membrane_dynamics`` with the external force enabled
(``docs/external_force.md``, ``data/example/push_release.params``).

Reads from the run directory:

    input.params                          the run's parameters
    <prefix>face.csv                      the triangles
    surfacepoint<prefix>.csv              one limit-surface frame per line
    <prefix>vertex_type_begin.csv         which vertices are clamped
    <prefix>external_force_vertices.csv   which vertices were pushed
    EnergyForce.csv                       the energy terms per iteration

and writes, into the run directory unless ``--out`` says otherwise:

    push_release.gif              the animation
    push_release_height.png       apex and disc height, and energy, against time
    push_release_snapshots.png    six stills through the push and the recovery

Usage:

    ~/anaconda3/bin/python analysis/push_release_animation.py <run dir> \
        [--prefix input] [--out DIR] [--frames 240] [--fps 20]

``--frames`` is how many trajectory frames the GIF samples (evenly in time,
with a hold on the first frame, the release and the last frame); the two
figures always use every frame.
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
from mpl_toolkits.mplot3d import Axes3D  # noqa: E402,F401

# ---------------------------------------------------------------------------
# Palette: a single blue ramp for height (a magnitude), two categorical slots
# for the two traces of each panel, neutral inks and surfaces.
# ---------------------------------------------------------------------------
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
             "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
             "#0d366b"]
HEIGHT_CMAP = LinearSegmentedColormap.from_list("slimed_blue", BLUE_RAMP)
SERIES_1 = "#2a78d6"   # apex height; bending energy
SERIES_2 = "#eb6834"   # disc mean height; area-constraint energy; the force arrows
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e5e4e0"
SHADE = "#f0efec"
SURFACE = "#fcfcfb"

plt.rcParams.update({
    "font.size": 10,
    "axes.edgecolor": INK_2,
    "axes.labelcolor": INK,
    "xtick.color": INK_2,
    "ytick.color": INK_2,
    "text.color": INK,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})


# ---------------------------------------------------------------------------
# Reading the run
# ---------------------------------------------------------------------------
def read_params(path):
    """The key = value pairs of a parameter file, comments stripped."""
    params = {}
    with open(path) as handle:
        for line in handle:
            line = line.split("#", 1)[0].strip()
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            params[key.strip()] = value.strip().rstrip(";").strip()
    return params


def read_frames(path):
    """A (frames, vertices, 3) array from a trajectory file.

    The modal line width, not the minimum, decides the vertex count, so a
    partial last line from a run still in progress is dropped rather than
    truncating every frame.
    """
    rows = []
    with open(path) as handle:
        for line in handle:
            line = line.strip().rstrip(",")
            if line:
                rows.append(np.fromstring(line, sep=","))
    if not rows:
        raise ValueError(f"{path} has no frames")
    widths = {}
    for row in rows:
        widths[len(row)] = widths.get(len(row), 0) + 1
    width = max(widths, key=widths.get)
    rows = [row for row in rows if len(row) == width]
    return np.array(rows).reshape(len(rows), -1, 3)


class PushRun:
    """A push-and-release run, with the quantities the figures need."""

    def __init__(self, directory, prefix="input"):
        self.directory = directory
        self.prefix = prefix
        params = read_params(os.path.join(directory, "input.params"))
        self.params = params
        self.dt = float(params.get("timeStep", 0.1))
        self.interval = int(params.get("meshpointOutputInterval", 1))
        self.max_iterations = int(params.get("maxIterations", 0))
        self.force_z = float(params.get("externalForceZ", 0.0))
        self.force_start = int(params.get("externalForceStartIteration", 0))
        stop = int(params.get("externalForceStopIteration", -1))
        self.force_stop = stop if stop >= 0 else self.max_iterations
        self.shape = params.get("externalForceShape", "disc")
        self.center = (float(params.get("externalForceCenterX", 0.0)),
                       float(params.get("externalForceCenterY", 0.0)))
        self.radius = float(params.get("externalForceRadius", 20.0))
        self.half_widths = (float(params.get("externalForceHalfWidthX", 20.0)),
                            float(params.get("externalForceHalfWidthY", 20.0)))

        def path(name):
            return os.path.join(directory, name)

        self.faces = np.loadtxt(path(f"{prefix}face.csv"), delimiter=",", dtype=int).reshape(-1, 3)
        self.frames = read_frames(path(f"surfacepoint{prefix}.csv"))
        types = np.loadtxt(path(f"{prefix}vertex_type_begin.csv"), delimiter=",",
                           usecols=[3], dtype=int)
        self.free = (types == 0)
        self.pushed = np.loadtxt(path(f"{prefix}external_force_vertices.csv"), dtype=int).reshape(-1)
        energy = np.genfromtxt(path("EnergyForce.csv"), delimiter=",", skip_header=1)
        if energy.ndim == 1:
            energy = energy.reshape(1, -1)

        n_frames = len(self.frames)
        self.iterations = np.arange(n_frames) * self.interval
        self.time = self.iterations * self.dt
        # Frame 0 is the configuration before the first step; frame k is the
        # configuration after k * interval steps.
        base = self.frames[0, :, 2]
        self.dz = self.frames[:, :, 2] - base[None, :]
        self.base_level = float(np.median(base[self.free]))
        self.apex = self.dz[:, self.free].max(axis=1)
        self.disc_mean = self.dz[:, self.pushed].mean(axis=1)
        # The record lags the configuration by one step: record[j] is the
        # energy of the configuration after j - 1 steps (j >= 1). At a frame
        # every hundred steps that is invisible; it is matched exactly anyway.
        index = np.minimum(self.iterations + 1, len(energy) - 1)
        index[0] = 0
        self.energy_bending = energy[index, 0]
        self.energy_area = energy[index, 1]

    # The force is on during the steps [start, stop). A frame shows the
    # configuration after its last step, so it is "under the force" when that
    # step was: start < iteration <= stop.
    def force_on(self, k):
        it = self.iterations[k]
        return self.force_start < it <= self.force_stop

    @property
    def release_time(self):
        return self.force_stop * self.dt

    @property
    def start_time(self):
        return self.force_start * self.dt

    def frame_at(self, t):
        """The frame nearest a time in microseconds."""
        return int(np.argmin(np.abs(self.time - t)))

    def region_outline(self, n=120):
        cx, cy = self.center
        if self.shape == "rectangle":
            hx, hy = self.half_widths
            x = np.array([cx - hx, cx + hx, cx + hx, cx - hx, cx - hx])
            y = np.array([cy - hy, cy - hy, cy + hy, cy + hy, cy - hy])
            return x, y
        angle = np.linspace(0.0, 2.0 * np.pi, n)
        return cx + self.radius * np.cos(angle), cy + self.radius * np.sin(angle)


# ---------------------------------------------------------------------------
# The surface
# ---------------------------------------------------------------------------
def style_3d_axes(ax):
    ax.set_facecolor(SURFACE)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        axis._axinfo["grid"]["color"] = GRID
        axis._axinfo["grid"]["linewidth"] = 0.5
    ax.tick_params(labelsize=8, colors=INK_2, pad=-1)
    ax.set_xlabel("x (nm)", fontsize=8, color=INK_2, labelpad=-4)
    ax.set_ylabel("y (nm)", fontsize=8, color=INK_2, labelpad=-4)
    ax.set_zlabel("z (nm)", fontsize=8, color=INK_2, labelpad=-4)


def draw_surface(ax, run, k, zmax, show_arrows=True):
    """Draw frame k of the limit surface, coloured by height above the base."""
    ax.cla()
    xyz = run.frames[k]
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    surface = ax.plot_trisurf(x, y, z, triangles=run.faces, cmap=HEIGHT_CMAP,
                              linewidth=0.25, edgecolor=(0.0, 0.0, 0.0, 0.16),
                              antialiased=True)
    surface.set_clim(run.base_level, run.base_level + zmax)

    # Where the force acts, drawn on the base plane; and the force itself,
    # as arrows from below the disc, while it is on.
    ox, oy = run.region_outline()
    ax.plot(ox, oy, np.full_like(ox, run.base_level), color=INK_2, lw=0.9, ls="--", zorder=5)
    if show_arrows and run.force_on(k):
        cx, cy = run.center
        r = 0.5 * run.radius if run.shape == "disc" else 0.5 * min(run.half_widths)
        ax_x = np.array([cx, cx + r, cx - r, cx, cx])
        ax_y = np.array([cy, cy, cy, cy + r, cy - r])
        # Drawn from just above the bump, upward: the camera looks down on
        # the sheet, so anything below it is hidden.
        length = 0.35 * zmax
        z_start = np.full(5, z[run.pushed].max() + 0.06 * zmax)
        ax.quiver(ax_x, ax_y, z_start, np.zeros(5), np.zeros(5), np.ones(5),
                  length=length, arrow_length_ratio=0.3, color=SERIES_2, lw=1.7, zorder=6)

    half_x = np.abs(run.frames[0, :, 0]).max()
    half_y = np.abs(run.frames[0, :, 1]).max()
    ax.set_xlim(-half_x, half_x)
    ax.set_ylim(-half_y, half_y)
    ax.set_zlim(run.base_level - 0.45 * zmax, run.base_level + 1.15 * zmax)
    ax.set_box_aspect((1.0, half_y / half_x, 0.42), zoom=1.3)
    ax.view_init(elev=24, azim=-57)
    style_3d_axes(ax)
    return surface


def frame_title(run, k):
    state = "force ON" if run.force_on(k) else ("force OFF, recovering" if run.time[k] > run.start_time else "at rest")
    return f"t = {run.time[k]:.1f} µs    {state}"


# ---------------------------------------------------------------------------
# The traces
# ---------------------------------------------------------------------------
def shade_force_window(ax, run, label=True):
    ax.axvspan(run.start_time, run.release_time, color=SHADE, zorder=0, lw=0)
    if label:
        ax.text(0.5 * (run.start_time + run.release_time), 0.97, "force on",
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=8.5, color=INK_2)


def style_trace_axes(ax):
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=9)


def draw_height_trace(ax, run, ymax):
    ax.cla()
    shade_force_window(ax, run)
    ax.plot(run.time, run.apex, color=SERIES_1, lw=1.6, label="apex (highest free vertex)")
    ax.plot(run.time, run.disc_mean, color=SERIES_2, lw=1.6, label="pushed disc, mean")
    ax.set_xlim(0.0, run.time[-1])
    ax.set_ylim(min(-1.5, 1.1 * run.dz[:, run.free].min()), ymax)
    ax.set_xlabel("time (µs)")
    ax.set_ylabel("height above the flat sheet (nm)")
    ax.legend(loc="upper right", frameon=False, fontsize=8.5)
    style_trace_axes(ax)


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
def write_animation(run, path, n_frames, fps, zmax):
    n = len(run.frames)
    picked = np.unique(np.linspace(0, n - 1, min(n_frames, n)).round().astype(int))
    release = run.frame_at(run.release_time)
    picked = np.unique(np.concatenate([picked, [release]]))
    # Holds: the flat sheet, the moment of release, and the end.
    hold_first = [picked[0]] * int(round(0.8 * fps))
    hold_release = [release] * int(round(0.8 * fps))
    hold_last = [picked[-1]] * int(round(1.5 * fps))
    sequence = hold_first + [k for k in picked if k < release] + hold_release + \
        [k for k in picked if k > release] + hold_last

    fig = plt.figure(figsize=(12.0, 5.2), dpi=100)
    grid = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], left=0.01, right=0.98,
                            top=0.90, bottom=0.12, wspace=0.08)
    ax3d = fig.add_subplot(grid[0, 0], projection="3d")
    ax2d = fig.add_subplot(grid[0, 1])
    ax2d.set_position([0.60, 0.15, 0.37, 0.70])
    ymax = 1.25 * max(run.apex.max(), 1.0)
    draw_height_trace(ax2d, run, ymax)
    cursor = ax2d.axvline(run.time[0], color=INK_2, lw=0.9)
    marker, = ax2d.plot([run.time[0]], [run.apex[0]], "o", color=SERIES_1, ms=6,
                        mec=SURFACE, mew=1.2)
    title = fig.text(0.02, 0.945, frame_title(run, 0), fontsize=12, color=INK, ha="left")
    n_pushed = len(run.pushed)
    fig.text(0.02, 0.905,
             f"{run.force_z:g} pN per vertex on {n_pushed} vertices of a "
             f"{run.shape} (r = {run.radius:g} nm), then released at "
             f"{run.release_time:g} µs", fontsize=9, color=INK_2, ha="left")

    def update(k):
        surface = draw_surface(ax3d, run, k, zmax)
        cursor.set_xdata([run.time[k], run.time[k]])
        marker.set_data([run.time[k]], [run.apex[k]])
        title.set_text(frame_title(run, k))
        return (surface, cursor, marker, title)

    anim = animation.FuncAnimation(fig, update, frames=sequence, interval=1000.0 / fps,
                                   blit=False)
    anim.save(path, writer=animation.PillowWriter(fps=fps), dpi=100)
    plt.close(fig)
    return len(sequence)


def write_height_figure(run, path):
    fig, (ax_h, ax_e) = plt.subplots(2, 1, figsize=(8.0, 6.0), dpi=150, sharex=True,
                                     gridspec_kw={"hspace": 0.12})
    ymax = 1.25 * max(run.apex.max(), 1.0)
    draw_height_trace(ax_h, run, ymax)
    ax_h.set_xlabel("")
    ax_h.set_title("Deformation and recovery of the pushed membrane", loc="left",
                   fontsize=11, color=INK)

    shade_force_window(ax_e, run, label=False)
    ax_e.plot(run.time, run.energy_bending, color=SERIES_1, lw=1.4, label="bending")
    ax_e.plot(run.time, run.energy_area, color=SERIES_2, lw=1.4, label="area constraint")
    ax_e.set_xlim(0.0, run.time[-1])
    ax_e.set_ylim(0.0, 1.15 * max(run.energy_bending.max(), run.energy_area.max()))
    ax_e.set_xlabel("time (µs)")
    ax_e.set_ylabel("energy (pN·nm)")
    ax_e.legend(loc="upper right", frameon=False, fontsize=8.5)
    style_trace_axes(ax_e)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def write_snapshots(run, path, zmax):
    push_len = run.release_time - run.start_time
    tail = run.time[-1] - run.release_time
    times = [run.start_time, run.start_time + 0.2 * push_len, run.release_time,
             run.release_time + 0.08 * tail, run.release_time + 0.4 * tail, run.time[-1]]
    fig = plt.figure(figsize=(13.0, 6.8), dpi=130)
    for i, t in enumerate(times):
        k = run.frame_at(t)
        ax = fig.add_subplot(2, 3, i + 1, projection="3d")
        draw_surface(ax, run, k, zmax)
        ax.set_title(frame_title(run, k), fontsize=10, color=INK, pad=2)
    fig.subplots_adjust(left=0.0, right=1.0, top=0.91, bottom=0.02, wspace=0.0, hspace=0.1)
    fig.suptitle(f"{run.force_z:g} pN per vertex on {len(run.pushed)} vertices, "
                 f"released at {run.release_time:g} µs", fontsize=11, color=INK_2, y=0.985)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def smooth(values, width):
    if width <= 1:
        return values
    kernel = np.ones(width) / width
    padded = np.concatenate([np.full(width // 2, values[0]), values,
                             np.full(width - 1 - width // 2, values[-1])])
    return np.convolve(padded, kernel, mode="valid")


def summarize(run):
    """The numbers worth quoting, printed as key = value lines."""
    width = max(1, int(round(0.5 / (run.interval * run.dt))))  # a 0.5 us window
    apex = smooth(run.apex, width)
    disc = smooth(run.disc_mean, width)
    release = run.frame_at(run.release_time)
    after = np.arange(len(run.time)) >= release
    lines = {
        "frames": len(run.frames),
        "duration_us": run.time[-1],
        "pushed_vertices": len(run.pushed),
        "free_vertices": int(run.free.sum()),
        "apex_max_nm": run.apex.max(),
        "apex_max_time_us": run.time[int(np.argmax(run.apex))],
        "apex_at_release_nm": apex[release],
        "disc_mean_at_release_nm": disc[release],
        "disc_mean_at_end_nm": disc[-1],
        "apex_at_end_nm": apex[-1],
        "area_energy_at_release_pNnm": run.energy_area[release],
        "area_energy_at_end_pNnm": run.energy_area[-1],
        "bending_energy_mean_pNnm": run.energy_bending[len(run.time) // 10:].mean(),
    }
    for fraction in (0.5, 0.25, 0.1):
        target = fraction * disc[release]
        hit = np.nonzero(after & (disc <= target))[0]
        lines[f"disc_below_{int(fraction * 100)}pct_after_release_us"] = \
            (run.time[hit[0]] - run.release_time) if len(hit) else float("nan")
    for key, value in lines.items():
        print(f"{key} = {value:.4g}" if isinstance(value, float) else f"{key} = {value}")
    return lines


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir")
    parser.add_argument("--prefix", default="input")
    parser.add_argument("--out", default=None, help="output directory (default: the run directory)")
    parser.add_argument("--frames", type=int, default=240)
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument("--no-gif", action="store_true", help="only the two figures")
    args = parser.parse_args(argv)

    run = PushRun(args.run_dir, args.prefix)
    out = args.out or args.run_dir
    os.makedirs(out, exist_ok=True)
    zmax = max(1.0, float(np.ceil(1.05 * run.apex.max())))

    summarize(run)
    write_height_figure(run, os.path.join(out, "push_release_height.png"))
    print("wrote", os.path.join(out, "push_release_height.png"))
    write_snapshots(run, os.path.join(out, "push_release_snapshots.png"), zmax)
    print("wrote", os.path.join(out, "push_release_snapshots.png"))
    if not args.no_gif:
        n = write_animation(run, os.path.join(out, "push_release.gif"), args.frames, args.fps, zmax)
        print("wrote", os.path.join(out, "push_release.gif"), f"({n} frames at {args.fps} fps)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
