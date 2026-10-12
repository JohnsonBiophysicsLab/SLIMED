#!/usr/bin/env python
"""Heat maps of the local energy and force on the membrane.

For runs written with ``localOutput = true`` (docs/bead_push_benchmark.md
section 8), which put beside every trajectory frame each face's bending and
mesh-term energy, mean curvature and area (``faceenergy<prefix>.csv``) and
each vertex's force by term with the bead's force on it
(``vertexforce<prefix>.csv``). Paints three local fields on the membrane:

    bending energy density   the face's bending energy over its area, pN/nm
    in-plane area strain     the face's area against the flat lattice cell,
                             stretched (+) or compressed (-)
    force magnitude          |bending + area-constraint + mesh-term + bead
                             force| on each vertex, pN

and writes into ``--out``:

    local_maps.png        the three fields on both membranes at one time (--at)
    local_profiles.png    ring averages of the fields against distance from the axis, at three times
    local.gif             the first two fields on both membranes, animated

Usage:

    ~/anaconda3/bin/python analysis/local_energy_heatmap.py --fluid DIR --frozen DIR \
        [--out DIR] [--prefix input] [--at END|<time in us>] [--frames 160] [--fps 16] [--no-gif]

One run alone: ``--run DIR``. The bead, if the run has one, is drawn and its
axis is the centre of the radial profiles; without one the sheet's centre is.
"""

import argparse
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import animation  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402
from mpl_toolkits.mplot3d import Axes3D  # noqa: E402,F401

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fluidity  # noqa: E402
from bead_push_analysis import (BLUE_RAMP, HEIGHT_CMAP, RUN_COLOR, BEAD_COLOR, INK, INK_2, GRID,  # noqa: E402
                                SURFACE, read_params, read_frames, draw_bead, style_trace)

# Sequential blue for magnitudes; blue - neutral - red for strain, which has a sign.
MAGNITUDE_CMAP = HEIGHT_CMAP
STRAIN_CMAP = LinearSegmentedColormap.from_list("slimed_diverging", ["#1c5cab", "#2a78d6", "#f0efec", "#e34948", "#a32626"])

FIELDS = {
    "bending": {"label": "bending energy density (pN/nm)", "cmap": MAGNITUDE_CMAP, "diverging": False},
    "strain": {"label": "in-plane area strain", "cmap": STRAIN_CMAP, "diverging": True},
    "force": {"label": "force magnitude (pN)", "cmap": MAGNITUDE_CMAP, "diverging": False},
}


class LocalRun:
    """A run with its local output: frames, faces, per-face and per-vertex fields."""

    def __init__(self, directory, label, prefix="input"):
        self.directory, self.label, self.prefix = directory, label, prefix
        p = read_params(os.path.join(directory, "input.params"))
        self.params = p
        self.dt = float(p.get("timeStep", 0.1))
        self.interval = int(p.get("meshpointOutputInterval", 1))
        self.flips = p.get("edgeFlipEnabled", "false") == "true"
        self.l_face = float(p.get("lFace", 5.0))
        self.run = fluidity.load_run(directory, prefix, self.interval)
        self.surface = read_frames(os.path.join(directory, f"surfacepoint{prefix}.csv"))
        face_rows = read_frames_wide(os.path.join(directory, f"faceenergy{prefix}.csv"), 4)
        vertex_rows = read_frames_wide(os.path.join(directory, f"vertexforce{prefix}.csv"), 12)
        n = min(len(self.run.coords), len(self.surface), len(face_rows), len(vertex_rows))
        self.n_frames = n
        self.control = self.run.coords[:n]
        self.surface = self.surface[:n]
        self.face_data = face_rows[:n]      # (n, faces, 4): E_bend, E_mesh, H, area
        self.vertex_data = vertex_rows[:n]  # (n, vertices, 12)
        self.iterations = self.run.frame_iterations[:n]
        self.time = self.iterations * self.dt
        self.free = np.zeros(self.control.shape[1], dtype=bool)
        self.free[self.run.free_vertices] = True
        self.base_level = float(np.median(self.surface[0, self.free, 2]))

        # The bead, if there is one.
        self.has_bead = p.get("beadEnabled", "false") == "true"
        self.radius = float(p.get("beadRadius", 10.0)) if self.has_bead else 0.0
        bead_path = os.path.join(directory, f"{prefix}bead.csv")
        if self.has_bead and os.path.exists(bead_path):
            bead = np.atleast_1d(np.genfromtxt(bead_path, delimiter=",", names=True))
            by_it = {int(r["iteration"]): r for r in bead}
            rows = [by_it.get(int(it)) for it in self.iterations]
            self.bead_center = np.array([[r["centerX"], r["centerY"], r["centerZ"]] if r is not None else [np.nan] * 3 for r in rows])
        else:
            self.bead_center = np.zeros((n, 3))
        self.bead_height = self.bead_center[:, 2] - self.base_level

        # The reference cell: every face's area on the flat lattice.
        area0 = self.face_data[0, :, 3]
        self.cell_area = float(np.median(area0[area0 > 0]))

    def faces_at(self, k):
        target = self.iterations[k]
        usable = [it for it in sorted(self.run.face_frames) if it <= target]
        return self.run.face_frames[usable[-1] if usable else min(self.run.face_frames)]

    def field(self, name, k):
        """The field at frame k: per face for bending and strain, per vertex for force."""
        if name == "bending":
            e, area = self.face_data[k, :, 0], self.face_data[k, :, 3]
            return np.where(area > 0, e / np.where(area > 0, area, 1.0), 0.0)
        if name == "strain":
            area = self.face_data[k, :, 3]
            return np.where(area > 0, area / self.cell_area - 1.0, 0.0)
        if name == "mesh":
            e, area = self.face_data[k, :, 1], self.face_data[k, :, 3]
            return np.where(area > 0, e / np.where(area > 0, area, 1.0), 0.0)
        if name == "force":
            d = self.vertex_data[k]
            total = d[:, 0:3] + d[:, 3:6] + d[:, 6:9] + d[:, 9:12]
            return np.linalg.norm(total, axis=1)
        raise KeyError(name)

    def per_face(self, name, k):
        """The field as one value per face, for painting.

        A vertex field is averaged over the face's corners. The strain is
        smoothed over the one-ring first -- each vertex takes the mean of its
        faces, each face the mean of its corners -- because face to face it
        is dominated by the thermal spread of edge lengths (a factor of four
        on a 5 nm mesh) and the systematic strain around the bead is a tenth
        of that. The profiles use the raw values: ring means smooth them.
        """
        values = self.field(name, k)
        faces = self.faces_at(k)
        if name == "force":
            return values[faces].mean(axis=1)
        if name == "strain":
            return self.smooth_one_ring(values, faces, k)
        return values

    def smooth_one_ring(self, face_values, faces, k):
        live = self.face_data[k, :, 3] > 0
        n_vertices = self.control.shape[1]
        total = np.zeros(n_vertices)
        count = np.zeros(n_vertices)
        for corner in range(3):
            np.add.at(total, faces[live, corner], face_values[live])
            np.add.at(count, faces[live, corner], 1.0)
        vertex_values = np.where(count > 0, total / np.where(count > 0, count, 1.0), 0.0)
        return np.where(live, vertex_values[faces].mean(axis=1), 0.0)

    def axis_xy(self, k):
        if self.has_bead:
            return self.bead_center[k, 0], self.bead_center[k, 1]
        return 0.0, 0.0

    def frame_at_height(self, height):
        hit = np.nonzero(self.bead_height >= height)[0]
        return int(hit[0]) if len(hit) else self.n_frames - 1

    def frame_at_time(self, t):
        return int(np.argmin(np.abs(self.time - t)))

    def radial_profile(self, name, k, edges):
        """Ring means of the field against distance from the axis, on face centroids or vertices."""
        cx, cy = self.axis_xy(k)
        if name == "force":
            pts = self.control[k][self.free]
            values = self.field(name, k)[self.free]
        else:
            faces = self.faces_at(k)
            pts = self.control[k][faces].mean(axis=1)
            values = self.field(name, k)
            keep = self.face_data[k, :, 3] > 0
            pts, values = pts[keep], values[keep]
        rho = np.hypot(pts[:, 0] - cx, pts[:, 1] - cy)
        which = np.digitize(rho, edges) - 1
        out = np.full(len(edges) - 1, np.nan)
        for b in range(len(edges) - 1):
            sel = which == b
            if sel.any():
                out[b] = values[sel].mean()
        return out


def read_frames_wide(path, width):
    """A (frames, items, width) array from a one-line-per-frame file."""
    rows = []
    with open(path) as handle:
        for line in handle:
            line = line.strip().rstrip(",")
            if line:
                rows.append(np.fromstring(line, sep=","))
    counts = {}
    for row in rows:
        counts[len(row)] = counts.get(len(row), 0) + 1
    length = max(counts, key=counts.get)
    rows = [row for row in rows if len(row) == length]
    return np.array(rows).reshape(len(rows), -1, width)


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------
def field_limits(runs, name, frames=None):
    """One colour scale per field across both runs (and the whole run for the animation)."""
    samples = []
    for r in runs:
        ks = frames if frames is not None else range(0, r.n_frames, max(1, r.n_frames // 40))
        for k in ks:
            if k < r.n_frames:
                samples.append(r.per_face(name, k))
    values = np.concatenate(samples)
    values = values[np.isfinite(values)]
    if FIELDS[name]["diverging"]:
        span = float(max(abs(np.percentile(values, 1)), abs(np.percentile(values, 99)), 1e-9))
        return -span, span
    return 0.0, float(max(np.percentile(values, 99), 1e-9))


def paint_membrane(ax, run, k, name, vmin, vmax, zmax, xlim, zoom=1.3):
    ax.cla()
    pts = run.surface[k]
    faces = run.faces_at(k)
    values = run.per_face(name, k)
    surface = ax.plot_trisurf(pts[:, 0], pts[:, 1], pts[:, 2], triangles=faces, cmap=FIELDS[name]["cmap"],
                              linewidth=0.15, edgecolor=(0.0, 0.0, 0.0, 0.12), antialiased=True)
    surface.set_array(values)
    if FIELDS[name]["diverging"]:
        surface.set_norm(TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax))
    else:
        surface.set_norm(Normalize(vmin=vmin, vmax=vmax))
    if run.has_bead and np.isfinite(run.bead_center[k, 2]):
        draw_bead(ax, run.bead_center[k], run.radius)
    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(-xlim, xlim)
    zlo, zhi = run.base_level - 0.25 * zmax, run.base_level + 1.1 * zmax
    ax.set_zlim(zlo, zhi)
    ax.set_box_aspect((1.0, 1.0, (zhi - zlo) / (2.0 * xlim)), zoom=zoom)
    ax.view_init(elev=30, azim=-58)
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
    return surface


def run_title(run, k):
    mode = ("fluid" if run.flips else "frozen") if len(run.label) == 0 else run.label
    if run.has_bead:
        return f"{mode}, bead at {run.bead_height[k]:+.0f} nm, t = {run.time[k]:.0f} µs"
    return f"{mode}, t = {run.time[k]:.1f} µs"


def geometry(runs):
    zmax = max(float(np.ceil(max(r.surface[:, :, 2][:, r.free].max() - r.base_level,
                                 (r.bead_height.max() + r.radius) if r.has_bead else 0.0) * 1.05)) for r in runs)
    zmax = max(zmax, 2.5 * max(r.radius for r in runs), 5.0)
    xlim = float(np.abs(runs[0].control[0, :, :2]).max())
    return zmax, xlim


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
def write_maps(runs, path, frames, names):
    zmax, xlim = geometry(runs)
    fig = plt.figure(figsize=(5.2 * len(runs) + 1.2, 4.2 * len(names)), dpi=130)
    grid = fig.add_gridspec(len(names), len(runs) + 1, width_ratios=[1.0] * len(runs) + [0.06],
                            left=0.0, right=0.95, top=0.95, bottom=0.02, wspace=0.0, hspace=0.12)
    for i, name in enumerate(names):
        vmin, vmax = field_limits(runs, name, [frames[j] for j in range(len(runs))])
        mappable = None
        for j, run in enumerate(runs):
            ax = fig.add_subplot(grid[i, j], projection="3d")
            mappable = paint_membrane(ax, run, frames[j], name, vmin, vmax, zmax, xlim, zoom=1.25)
            ax.set_title(run_title(run, frames[j]), fontsize=9, color=INK, pad=0)
        cax = fig.add_subplot(grid[i, len(runs)])
        cbar = fig.colorbar(mappable, cax=cax)
        cbar.set_label(FIELDS[name]["label"], fontsize=9, color=INK)
        cbar.ax.tick_params(labelsize=8, colors=INK_2)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def write_profiles(runs, path, times, names):
    radius = max(r.radius for r in runs)
    xlim = float(np.abs(runs[0].control[0, :, :2]).max())
    edges = np.arange(0.0, xlim + 2.5, 2.5)
    centres = 0.5 * (edges[1:] + edges[:-1])
    fig, axes = plt.subplots(len(names), len(times), figsize=(4.6 * len(times), 3.0 * len(names)), dpi=140,
                             sharex=True, gridspec_kw={"hspace": 0.18, "wspace": 0.22})
    axes = np.atleast_2d(axes)
    for j, t in enumerate(times):
        for i, name in enumerate(names):
            ax = axes[i, j]
            for r in runs:
                k = r.frame_at_time(t)
                ax.plot(centres, r.radial_profile(name, k, edges), color=RUN_COLOR.get(r.label, INK), lw=1.6,
                        label=r.label)
            if radius > 0:
                ax.axvline(radius, color=INK_2, lw=0.9, ls="--")
            if FIELDS[name]["diverging"]:
                ax.axhline(0.0, color=GRID, lw=0.8)
            style_trace(ax)
            if i == 0:
                k = runs[0].frame_at_time(t)
                ax.set_title(f"t = {runs[0].time[k]:.0f} µs" +
                             (f", bead at {runs[0].bead_height[k]:+.0f} nm" if runs[0].has_bead else ""),
                             fontsize=10, color=INK)
            if j == 0:
                ax.set_ylabel(FIELDS[name]["label"], fontsize=9)
            if i == len(names) - 1:
                ax.set_xlabel("distance from the axis (nm)", fontsize=9)
            if i == 0 and j == 0:
                ax.legend(loc="upper right", frameon=False, fontsize=8.5)
    if radius > 0:
        axes[0, 0].text(radius, 0.98, " bead radius", transform=axes[0, 0].get_xaxis_transform(),
                        fontsize=8, color=INK_2, va="top")
    fig.suptitle("Ring averages of the local fields", fontsize=11, color=INK, y=1.0)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def write_animation(runs, path, n_frames, fps, names):
    zmax, xlim = geometry(runs)
    n = min(r.n_frames for r in runs)
    picked = np.unique(np.linspace(0, n - 1, min(n_frames, n)).round().astype(int))
    hold = int(round(0.8 * fps))
    sequence = [picked[0]] * hold + list(picked) + [picked[-1]] * (2 * hold)
    limits = {name: field_limits(runs, name, picked[::4]) for name in names}

    fig = plt.figure(figsize=(5.2 * len(runs) + 1.2, 4.2 * len(names)), dpi=100)
    grid = fig.add_gridspec(len(names), len(runs) + 1, width_ratios=[1.0] * len(runs) + [0.06],
                            left=0.0, right=0.88, top=0.93, bottom=0.02, wspace=0.0, hspace=0.12)
    axes = [[fig.add_subplot(grid[i, j], projection="3d") for j in range(len(runs))] for i in range(len(names))]
    caxes = [fig.add_subplot(grid[i, len(runs)]) for i in range(len(names))]
    titles = [[fig.text(0.0, 0.0, "") for _ in runs] for _ in names]
    stamp = fig.text(0.02, 0.965, "", fontsize=12, color=INK)
    first = [True]

    def update(k):
        artists = []
        for i, name in enumerate(names):
            vmin, vmax = limits[name]
            for j, run in enumerate(runs):
                ax = axes[i][j]
                mappable = paint_membrane(ax, run, k, name, vmin, vmax, zmax, xlim, zoom=1.25)
                ax.set_title(run_title(run, k), fontsize=9, color=INK, pad=0)
                artists.append(mappable)
                if first[0] and j == len(runs) - 1:
                    cbar = fig.colorbar(mappable, cax=caxes[i])
                    cbar.set_label(FIELDS[name]["label"], fontsize=9, color=INK)
                    cbar.ax.tick_params(labelsize=8, colors=INK_2)
        first[0] = False
        stamp.set_text(f"t = {runs[0].time[k]:.1f} µs")
        return artists

    anim = animation.FuncAnimation(fig, update, frames=sequence, interval=1000.0 / fps, blit=False)
    anim.save(path, writer=animation.PillowWriter(fps=fps), dpi=100)
    plt.close(fig)
    return len(sequence)


def summarize(run, k):
    cx, cy = run.axis_xy(k)
    faces = run.faces_at(k)
    centroids = run.control[k][faces].mean(axis=1)
    rho_f = np.hypot(centroids[:, 0] - cx, centroids[:, 1] - cy)
    live = run.face_data[k, :, 3] > 0
    near_f = live & (rho_f < 1.5 * max(run.radius, run.l_face))
    far_f = live & (rho_f > 3.0 * max(run.radius, run.l_face))
    pts = run.control[k][run.free]
    rho_v = np.hypot(pts[:, 0] - cx, pts[:, 1] - cy)
    near_v = rho_v < 1.5 * max(run.radius, run.l_face)
    far_v = rho_v > 3.0 * max(run.radius, run.l_face)
    bending, strain, mesh = run.field("bending", k), run.field("strain", k), run.field("mesh", k)
    force = run.field("force", k)[run.free]
    out = {
        "label": run.label, "frame": k, "time_us": run.time[k],
        "bending_density_max_pN_per_nm": float(bending[live].max()),
        "bending_density_near_axis_mean": float(bending[near_f].mean()) if near_f.any() else float("nan"),
        "bending_density_far_mean": float(bending[far_f].mean()) if far_f.any() else float("nan"),
        "area_strain_near_axis_mean": float(strain[near_f].mean()) if near_f.any() else float("nan"),
        "area_strain_near_axis_min_max": (float(strain[near_f].min()), float(strain[near_f].max())) if near_f.any() else (float("nan"),) * 2,
        "area_strain_far_mean": float(strain[far_f].mean()) if far_f.any() else float("nan"),
        "mesh_energy_density_max": float(mesh[live].max()),
        "mesh_energy_density_near_axis_mean": float(mesh[near_f].mean()) if near_f.any() else float("nan"),
        "force_magnitude_near_axis_mean_pN": float(force[near_v].mean()) if near_v.any() else float("nan"),
        "force_magnitude_far_mean_pN": float(force[far_v].mean()) if far_v.any() else float("nan"),
        "force_magnitude_max_pN": float(force.max()),
    }
    for key, value in out.items():
        print(f"{key} = {value}")
    print()
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fluid")
    parser.add_argument("--frozen")
    parser.add_argument("--run", help="a single run, labelled by its flips setting")
    parser.add_argument("--out", default=None)
    parser.add_argument("--prefix", default="input")
    parser.add_argument("--at", default="END", help="time (us) of the static maps, or END")
    parser.add_argument("--frames", type=int, default=160)
    parser.add_argument("--fps", type=int, default=16)
    parser.add_argument("--no-gif", action="store_true")
    args = parser.parse_args(argv)

    runs = []
    if args.run:
        r = LocalRun(args.run, "", args.prefix)
        r.label = "fluid" if r.flips else "frozen"
        runs.append(r)
    else:
        if not (args.fluid and args.frozen):
            parser.error("give --fluid and --frozen, or --run")
        runs = [LocalRun(args.fluid, "fluid", args.prefix), LocalRun(args.frozen, "frozen", args.prefix)]
    out = args.out or runs[0].directory
    os.makedirs(out, exist_ok=True)

    if args.at.upper() == "END":
        frames = [r.n_frames - 1 for r in runs]
    else:
        frames = [r.frame_at_time(float(args.at)) for r in runs]
    for r, k in zip(runs, frames):
        summarize(r, k)

    names = ["bending", "strain", "force"]
    write_maps(runs, os.path.join(out, "local_maps.png"), frames, names)
    print("wrote", os.path.join(out, "local_maps.png"))
    end = min(r.time[-1] for r in runs)
    if runs[0].has_bead:
        times = [runs[0].time[runs[0].frame_at_height(0.0)], runs[0].time[runs[0].frame_at_height(runs[0].radius)], end]
    else:
        times = [0.25 * end, 0.5 * end, end]
    write_profiles(runs, os.path.join(out, "local_profiles.png"), times, names)
    print("wrote", os.path.join(out, "local_profiles.png"))
    if not args.no_gif:
        n = write_animation(runs, os.path.join(out, "local.gif"), args.frames, args.fps, ["bending", "strain"])
        print("wrote", os.path.join(out, "local.gif"), f"({n} frames at {args.fps} fps)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
