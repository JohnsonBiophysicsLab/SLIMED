#!/usr/bin/env python
"""
Analysis of the sticky-bead matrix (analysis/sticky_bead_runs.py): a bead
bonded to a periodic membrane by the full harmonic bond, pushed through it and
held, under two area reservoirs and several fluidities.

    python analysis/sticky_bead_analysis.py --root ~/slimed_runs/sticky_bead/matrix

reads every run directory under --root (or the ones named with --runs) and
writes into --out (default <root>/analysis):

  traces.png            force on the bead, wrapped fraction, tension, reservoir
                        area, neck radius and bending energy against time
  hold_relaxation.png   the hold phase after the bead stops, with exponential fits
  energy_maps_<t>.png   top-view maps of the bending energy density on every run
  energy_profiles.png   ring averages of the bending energy density around the axis
  spectra.png           Fourier analysis of the height field during the hold:
                        the static imprint of the bead and the fluctuation
                        spectrum with its kc, sigma fit
  anim_<run>.gif        the membrane coloured by bending energy density, the
                        bead, and the force trace, per run
  compare.gif           up to four runs side by side (--compare)
  stats.md, stats.csv   one row per run

Run it with ~/anaconda3/bin/python (numpy, scipy, matplotlib).
"""
import argparse
import csv
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fluidity  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import animation  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e5e4e0"
SURFACE = "#fcfcfb"
BEAD_COLOR = "#8f8c85"
ENERGY_CMAP = "magma"
plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK_2, "axes.labelcolor": INK,
    "xtick.color": INK_2, "ytick.color": INK_2, "text.color": INK,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})

# Colour by fluidity, line style by reservoir.
FLUIDITY_COLORS = {"frozen": "#4c5c68", "nu0.05": "#2a7f9e", "nu0.5": "#d9822b", "nu5": "#b5342e"}
RESERVOIR_STYLES = {"box": "-", "s0global": "--", "s0local": ":", "fixed": "-."}
RESERVOIR_LABELS = {"box": "box scaling", "s0global": "adaptive S0 (global)",
                    "s0local": "adaptive S0 (local)", "fixed": "fixed box"}
FLUIDITY_LABELS = {"frozen": "frozen", "nu0.05": "ν = 0.05/µs", "nu0.5": "ν = 0.5/µs", "nu5": "ν = 5/µs"}


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------
def read_params(path):
    params = {}
    with open(path) as handle:
        for line in handle:
            line = line.split("#", 1)[0].strip()
            if "=" in line:
                key, value = line.split("=", 1)
                params[key.strip()] = value.strip()
    return params


def read_frames_wide(path, width):
    """Frames of a per-item file: (n_frames, n_items, width); a partial last line is dropped."""
    rows = []
    with open(path) as handle:
        for line in handle:
            line = line.strip().rstrip(",")
            if line:
                rows.append(np.fromstring(line, sep=","))
    if not rows:
        return np.zeros((0, 0, width))
    counts = {}
    for r in rows:
        counts[len(r)] = counts.get(len(r), 0) + 1
    n = max(counts, key=counts.get)
    rows = [r for r in rows if len(r) == n]
    return np.array(rows).reshape(len(rows), -1, width)


def read_log(path):
    if not os.path.exists(path):
        return None
    table = np.genfromtxt(path, delimiter=",", names=True, invalid_raise=False)
    return np.atleast_1d(table)


def rows_at(table, iterations):
    """The log rows at these iterations (the last row at or before each)."""
    its = table["iteration"].astype(int)
    order = np.argsort(its)
    its = its[order]
    idx = np.searchsorted(its, iterations, side="right") - 1
    idx = np.clip(idx, 0, len(its) - 1)
    return table[order[idx]]


class StickyRun:
    def __init__(self, directory, prefix="input"):
        self.directory = directory
        self.name = os.path.basename(os.path.normpath(directory))
        self.prefix = prefix
        p = read_params(os.path.join(directory, "input.params"))
        self.params = p
        self.dt = float(p.get("timeStep", 0.001))
        self.interval = int(p.get("meshpointOutputInterval", 1))
        self.l_face = float(p.get("lFace", 5.0))
        self.kbt = float(p.get("KBT", 4.17))
        self.kc = float(p.get("kcMembraneBending", 83.4))
        self.K = float(p.get("usMembraneStretching", 250.0))
        self.radius = float(p.get("beadRadius", 10.0))
        self.flips = p.get("edgeFlipEnabled", "false") == "true"
        self.rate = float(p.get("edgeFlipAttemptRate", 0.0)) if self.flips else 0.0
        self.fluidity = "frozen" if not self.flips else {0.05: "nu0.05", 0.5: "nu0.5", 5.0: "nu5"}.get(self.rate, f"nu{self.rate:g}")
        if p.get("elasticBoxEnabled", "false") == "true":
            self.reservoir = "box"
        elif p.get("lipidFlowMode", "off") == "global":
            self.reservoir = "s0global"
        elif p.get("lipidFlowMode", "off") == "local":
            self.reservoir = "s0local"
        else:
            self.reservoir = "fixed"
        self.label = f"{RESERVOIR_LABELS[self.reservoir]}, {FLUIDITY_LABELS.get(self.fluidity, self.fluidity)}"
        self.stop_iteration = int(p.get("beadStopIteration", -1))
        self.t_stop = self.stop_iteration * self.dt
        self.vertex_area = np.sqrt(3.0) / 2.0 * self.l_face ** 2
        self.sphere_cells = 4.0 * np.pi * self.radius ** 2 / self.vertex_area

        self.run = fluidity.load_run(directory, prefix, self.interval)
        # The real tile is the Free vertices of the type file; the images are
        # not (and under the elastic box they scale with it, which defeats the
        # displacement-correlation test fluidity.load_run falls back on).
        self.run.free_vertices = fluidity._free_vertices_from_types(self.run)
        self.surface = read_frames_wide(os.path.join(directory, f"surfacepoint{prefix}.csv"), 3)
        self.face_data = read_frames_wide(os.path.join(directory, f"faceenergy{prefix}.csv"), 4)
        n = min(len(self.run.coords), len(self.surface), len(self.face_data))
        self.n_frames = n
        self.control = self.run.coords[:n]
        self.surface = self.surface[:n]
        self.face_data = self.face_data[:n]
        self.iterations = self.run.frame_iterations[:n]
        self.time = self.iterations * self.dt
        self.free = np.zeros(self.control.shape[1], dtype=bool)
        self.free[self.run.free_vertices] = True
        self.n_free = int(self.free.sum())
        self.base_level = float(np.median(self.control[0, self.free, 2]))

        bead = read_log(os.path.join(directory, f"{prefix}bead.csv"))
        b = rows_at(bead, self.iterations)
        self.bead_center = np.column_stack([b["centerX"], b["centerY"], b["centerZ"]])
        self.force = np.column_stack([b["forceOnBeadX"], b["forceOnBeadY"], b["forceOnBeadZ"]])
        # A sticky bead keeps the vertices it bonded; a bead that only pushes
        # (beadPotential = wall, no well) has contacts instead, the vertices
        # its wall is holding off. The same column serves both as the wrap.
        self.potential = p.get("beadPotential", "wall")
        self.sticky = self.potential == "harmonic"
        self.contact_word = "bonded" if self.sticky else "in contact"
        if self.sticky:
            self.bonded = b["bonded"].astype(int) if "bonded" in b.dtype.names else b["adhered"].astype(int)
        else:
            self.bonded = b["contacts"].astype(int)
        self.bond_energy = b["energy_pNnm"]
        self.max_stretch = b["maxStretch_nm"] if "maxStretch_nm" in b.dtype.names else np.zeros(n)
        self.wrapped = self.bonded / self.sphere_cells

        # Lattice of the real tile, from the free vertices as built.
        x0 = self.control[0, self.free, 0]
        y0 = self.control[0, self.free, 1]
        # Rows alternate by half a cell in x, so the x range is (n_x - 1/2) cells.
        self.n_x = int(round(np.ptp(x0) / self.l_face + 0.5))
        self.n_y = int(round(np.ptp(y0) / (self.l_face * np.sqrt(3.0) / 2.0))) + 1
        self.L0 = np.array([self.n_x * self.l_face, self.n_y * self.l_face * np.sqrt(3.0) / 2.0])
        self.origin = np.array([x0.min(), y0.min()])

        # The reservoir: tension and the area drawn.
        self.scale = np.ones(n)
        box = read_log(os.path.join(directory, f"{prefix}box.csv"))
        lipid = read_log(os.path.join(directory, f"{prefix}lipid.csv"))
        area_all = self.face_data[:, :, 3].sum(axis=1)
        if box is not None:
            r = rows_at(box, self.iterations)
            self.scale = r["scaleX"]
            self.tension = r["tension_pN_nm"]
            self.reservoir_area = r["frameArea_nm2"] - r["frameArea_nm2"][0]
            self.tension_label = "frame tension (box virial)"
        elif lipid is not None:
            r = rows_at(lipid, self.iterations)
            self.tension = r["meanTension_pN_nm"]
            self.reservoir_area = r["sumRelaxedArea_nm2"] - r["sumRelaxedArea_nm2"][0]
            self.tension_label = "mean tension K (S - S0)/S0"
        else:
            self.tension = self.K * (area_all / area_all[0] - 1.0)
            self.reservoir_area = np.zeros(n)
            self.tension_label = "K (A - A0)/A0"
        self.area_all = area_all
        # One measure every reservoir shares: the membrane's excess over its
        # projected area, A - A_p, which is what the bead has pulled out of
        # the plane (the sum of the faces' areas is the tile's membrane area;
        # ghost faces carry none).
        self.projected_area = self.L0.prod() * self.scale ** 2
        self.excess_area = area_all - self.projected_area

        # Energies, one row per step in EnergyForce.csv; sampled at the frames.
        ef = os.path.join(directory, "EnergyForce.csv")
        self.e_bend = np.full(n, np.nan)
        self.e_area = np.full(n, np.nan)
        if os.path.exists(ef):
            table = np.genfromtxt(ef, delimiter=",", skip_header=1, invalid_raise=False)
            table = np.atleast_2d(table)
            idx = np.clip(self.iterations, 0, len(table) - 1)
            self.e_bend = table[idx, 0]
            self.e_area = table[idx, 1]

        # Flips: accepted per microsecond in the windows between frames.
        self.flips_per_us = np.zeros(n)
        if self.run.flip_accepted is not None and n > 1:
            acc = self.run.flip_iteration[self.run.flip_accepted]
            counts, _ = np.histogram(acc, bins=np.append(self.iterations, self.iterations[-1] + self.interval))
            self.flips_per_us = counts / (self.interval * self.dt)

        self.neck, self.apex = self._neck_and_apex()

    # -- geometry ---------------------------------------------------------------
    def faces_at(self, k):
        target = self.iterations[k]
        usable = [it for it in sorted(self.run.face_frames) if it <= target]
        return self.run.face_frames[usable[-1] if usable else min(self.run.face_frames)]

    def bending_density(self, k):
        e, area = self.face_data[k, :, 0], self.face_data[k, :, 3]
        return np.where(area > 0, e / np.where(area > 0, area, 1.0), 0.0)

    def energy_faces(self, k):
        """The energy-carrying faces at frame k and their bending energy density.

        The mesh's periodic image band holds copy faces that carry no energy
        (area 0 in the per-face output). Under edge flips the code never
        re-mirrors them from the flipped real faces, so their connectivity
        goes stale and they stretch across the band, while every force on the
        tile stays exact (docs/mixed_boundary_conditions.md section 9). They
        must not be drawn or binned: this is the set to use."""
        faces = self.faces_at(k)
        carry = self.face_data[k, :, 3] > 0
        return faces[carry], self.bending_density(k)[carry]

    def _neck_and_apex(self):
        """Per frame: the narrowest ring of membrane below the bead's centre
        (mean radial distance of the free vertices in the narrowest 2.5 nm
        slab above the sheet), and the highest free vertex."""
        neck = np.full(self.n_frames, np.nan)
        apex = np.full(self.n_frames, np.nan)
        for k in range(self.n_frames):
            pts = self.control[k, self.free]
            c = self.bead_center[k]
            apex[k] = pts[:, 2].max() - self.base_level
            if c[2] - self.radius < self.base_level + 3.0:
                continue
            r = np.hypot(pts[:, 0] - c[0], pts[:, 1] - c[1])
            z = pts[:, 2]
            sel = (z > self.base_level + 3.0) & (z < c[2]) & (r < 2.5 * self.radius)
            if sel.sum() < 3:
                continue
            edges = np.arange(self.base_level + 3.0, c[2] + 2.5, 2.5)
            best = np.inf
            for lo, hi in zip(edges[:-1], edges[1:]):
                m = sel & (z >= lo) & (z < hi)
                if m.sum() >= 3:
                    best = min(best, r[m].mean())
            if np.isfinite(best):
                neck[k] = best
        return neck, apex

    # -- the height field on a grid, for the Fourier analysis -------------------
    def height_grid(self, k, envelope=True):
        """The membrane's height on the tile's lattice grid (n_y x n_x cells).

        envelope=True takes the highest vertex of each cell (the bead's bump
        included); False the mean of the cell's vertices, which is the sheet
        itself where it is single-valued. Also returns the cell's distance to
        the bead's axis, wrapped, for masking the footprint."""
        pts = self.control[k]
        L = self.L0 * self.scale[k]
        origin = self.origin * self.scale[k] - 0.5 * self.l_face
        u = np.mod(pts[:, 0] - origin[0], L[0])
        v = np.mod(pts[:, 1] - origin[1], L[1])
        ix = np.minimum((u / L[0] * self.n_x).astype(int), self.n_x - 1)
        iy = np.minimum((v / L[1] * self.n_y).astype(int), self.n_y - 1)
        if envelope:
            h = np.full((self.n_y, self.n_x), -np.inf)
            np.maximum.at(h, (iy, ix), pts[:, 2])
            empty = ~np.isfinite(h)
        else:
            # Linear interpolation of the control net at the cell centres: a
            # vertex that drifts across a cell boundary then moves the field
            # continuously, where binning would make it jump. The images
            # around the tile close the Delaunay hull over every cell.
            from scipy.interpolate import LinearNDInterpolator
            cx = (np.arange(self.n_x) + 0.5) * L[0] / self.n_x + origin[0]
            cy = (np.arange(self.n_y) + 0.5) * L[1] / self.n_y + origin[1]
            CX, CY = np.meshgrid(cx, cy)
            interp = LinearNDInterpolator(pts[:, :2], pts[:, 2])
            h = interp(CX, CY)
            empty = ~np.isfinite(h)
            h = np.where(empty, 0.0, h)
        if empty.any():
            from scipy.ndimage import distance_transform_edt
            _, (iy_near, ix_near) = distance_transform_edt(empty, return_indices=True)
            h = h[iy_near, ix_near]
        cx = (np.arange(self.n_x) + 0.5) * L[0] / self.n_x
        cy = (np.arange(self.n_y) + 0.5) * L[1] / self.n_y
        c = self.bead_center[k]
        dx = np.mod(cx - (c[0] - origin[0]) + L[0] / 2.0, L[0]) - L[0] / 2.0
        dy = np.mod(cy - (c[1] - origin[1]) + L[1] / 2.0, L[1]) - L[1] / 2.0
        r = np.hypot(dx[None, :], dy[:, None])
        return h - self.base_level, L, r

    def hold_frames(self, skip_us=20.0):
        """Frames of the hold phase, the first skip_us after the stop left out."""
        if self.stop_iteration < 0:
            return np.arange(self.n_frames)
        return np.nonzero(self.time >= self.t_stop + skip_us)[0]


# ---------------------------------------------------------------------------
# Measurements
# ---------------------------------------------------------------------------
def running_mean(values, window):
    if window <= 1:
        return values
    kernel = np.ones(window) / window
    out = np.convolve(np.nan_to_num(values), kernel, mode="same")
    return out


def fit_relaxation(t, y):
    """y(t) = y_inf + a exp(-(t - t0)/tau) on the hold phase; returns (y_inf, a, tau) or None."""
    good = np.isfinite(y)
    t, y = t[good], y[good]
    if len(t) < 10:
        return None
    t0 = t[0]
    y_inf0 = np.mean(y[-max(5, len(y) // 5):])
    a0 = y[0] - y_inf0
    tau0 = max(1.0, (t[-1] - t0) / 5.0)
    try:
        popt, _ = curve_fit(lambda tt, yi, a, tau: yi + a * np.exp(-(tt - t0) / tau), t, y,
                            p0=(y_inf0, a0, tau0), bounds=([-np.inf, -np.inf, 0.05], [np.inf, np.inf, 1e4]), maxfev=20000)
        return tuple(popt)
    except Exception:
        return None


def ring_profile(run, k, n_bins=12, r_max=None):
    """Mean bending energy density in rings around the bead's axis (face centroids), wrapped into the tile."""
    faces = run.faces_at(k)
    pts = run.control[k]
    centroids = pts[faces].mean(axis=1)
    L = run.L0 * run.scale[k]
    c = run.bead_center[k]
    dx = np.mod(centroids[:, 0] - c[0] + L[0] / 2.0, L[0]) - L[0] / 2.0
    dy = np.mod(centroids[:, 1] - c[1] + L[1] / 2.0, L[1]) - L[1] / 2.0
    r = np.hypot(dx, dy)
    dens = run.bending_density(k)
    area = run.face_data[k, :, 3]
    ok = area > 0
    r_max = r_max or L.min() / 2.0
    edges = np.linspace(0.0, r_max, n_bins + 1)
    prof = np.full(n_bins, np.nan)
    for i in range(n_bins):
        m = ok & (r >= edges[i]) & (r < edges[i + 1])
        if m.any():
            prof[i] = np.average(dens[m], weights=area[m])
    return 0.5 * (edges[:-1] + edges[1:]), prof


def spectrum(run, frames, mask_radii=1.5):
    """Fourier analysis of the hold phase on the tile's lattice grid.

    Static imprint: |<h_q>|^2 of the upper envelope (the bead's bump included).
    Fluctuations: <|h_q - <h_q>|^2> of the sheet's mean height per cell with
    the bead's footprint (mask_radii x R plus one cell around the axis) held
    at its time mean, so that the sphere's and the neck's vertices, which the
    gridding cannot follow, do not put a white floor under the spectrum. The
    fit is kT / (A (kc q^4 + sigma q^2)) for q below 0.35 /nm, clear of the lattice end."""
    env, sheet, radii = [], [], None
    for k in frames:
        h, L, r = run.height_grid(k, envelope=True)
        env.append(h)
        h2, _, _ = run.height_grid(k, envelope=False)
        sheet.append(h2)
        radii = r
    env, sheet = np.array(env), np.array(sheet)
    ny, nx = env.shape[1:]
    L = run.L0 * np.mean(run.scale[frames])
    A = L[0] * L[1]
    hq_env = np.fft.fft2(env, axes=(1, 2)) / (nx * ny)
    static = np.abs(hq_env.mean(axis=0)) ** 2
    footprint = radii < mask_radii * run.radius + run.l_face
    masked = np.where(footprint[None, :, :], sheet.mean(axis=0)[None, :, :], sheet)
    hq = np.fft.fft2(masked, axes=(1, 2)) / (nx * ny)
    fluct = (np.abs(hq - hq.mean(axis=0)) ** 2).mean(axis=0)
    # The masked cells carry no fluctuation: rescale by the sheet's share.
    fluct = fluct / (1.0 - footprint.mean())
    qx = 2.0 * np.pi * np.fft.fftfreq(nx, d=L[0] / nx)
    qy = 2.0 * np.pi * np.fft.fftfreq(ny, d=L[1] / ny)
    QX, QY = np.meshgrid(qx, qy)
    q = np.hypot(QX, QY)
    q_max = min(np.pi * nx / L[0], np.pi * ny / L[1])
    edges = np.linspace(0.0, q_max, 14)[1:]
    qs, st, fl = [], [], []
    for lo, hi in zip(np.append(0.0, edges[:-1]), edges):
        m = (q > max(lo, 1e-9)) & (q <= hi)
        if m.any():
            qs.append(q[m].mean())
            st.append(static[m].mean())
            fl.append(fluct[m].mean())
    qs, st, fl = np.array(qs), np.array(st), np.array(fl)
    fit = None
    m = qs < 0.35
    if m.sum() >= 4:
        try:
            popt, _ = curve_fit(lambda qq, kc, sig: np.log(run.kbt / (A * (kc * qq ** 4 + sig * qq ** 2))),
                                qs[m], np.log(fl[m]), p0=(run.kc, 0.5), bounds=([1.0, -5.0], [2000.0, 50.0]))
            fit = tuple(popt)
        except Exception:
            fit = None
    return qs, st, fl, fit, A


def summarize(run):
    k_stop = int(np.searchsorted(run.iterations, run.stop_iteration)) if run.stop_iteration >= 0 else run.n_frames - 1
    k_stop = min(k_stop, run.n_frames - 1)
    last = slice(max(0, run.n_frames - max(2, run.n_frames // 10)), run.n_frames)
    hold = run.hold_frames(skip_us=0.0)
    out = {
        "run": run.name, "potential": run.potential, "reservoir": run.reservoir, "fluidity": run.fluidity,
        "t_end_us": run.time[-1],
        "bonded_at_stop": int(run.bonded[k_stop]), "bonded_end": int(np.mean(run.bonded[last])),
        "wrapped_end": float(np.mean(run.wrapped[last])),
        "Fz_at_stop_pN": float(running_mean(run.force[:, 2], 5)[k_stop]),
        "Fz_end_pN": float(np.mean(run.force[last, 2])),
        "Flat_end_pN": float(np.mean(np.hypot(run.force[last, 0], run.force[last, 1]))),
        "Fz_peak_pN": float(np.nanmax(running_mean(run.force[:, 2], 5))),
        "Fz_min_pN": float(np.nanmin(running_mean(run.force[:, 2], 5))),
        "tension_end_pN_nm": float(np.nanmean(run.tension[last])),
        "reservoir_area_end_nm2": float(np.mean(run.reservoir_area[last])),
        "excess_area_end_nm2": float(np.mean(run.excess_area[last])),
        "neck_end_nm": float(np.nanmean(run.neck[last])) if np.isfinite(run.neck[last]).any() else np.nan,
        "apex_end_nm": float(np.mean(run.apex[last])),
        "E_bend_end_pNnm": float(np.nanmean(run.e_bend[last])) if np.isfinite(run.e_bend[last]).any() else np.nan,
        "flips_per_us_end": float(np.mean(run.flips_per_us[last])),
        "tau_force_us": np.nan, "tau_bonded_us": np.nan, "tau_tension_us": np.nan,
    }
    if len(hold) > 10:
        t = run.time[hold]
        for key, series in (("tau_force_us", running_mean(run.force[:, 2], 5)[hold]),
                            ("tau_bonded_us", run.bonded[hold].astype(float)),
                            ("tau_tension_us", run.tension[hold])):
            fit = fit_relaxation(t, series)
            if fit is not None and np.ptp(series[np.isfinite(series)]) > 0 and fit[2] < 0.9 * (t[-1] - t[0]) * 10:
                out[key] = float(fit[2])
    try:
        alt = fluidity.altitude_stats(run.run, -1)
        crease = fluidity.crease_stats(run.run, -1)
        out["min_altitude_nm"] = float(alt.get("min", np.nan))
        out["creases_over_90"] = int(crease.get("over_90", crease.get("count_over_90", 0)))
    except Exception:
        out["min_altitude_nm"] = np.nan
        out["creases_over_90"] = -1
    return out


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def style(ax):
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=8.5)


def line_style(run):
    return dict(color=FLUIDITY_COLORS.get(run.fluidity, "#333333"), ls=RESERVOIR_STYLES.get(run.reservoir, "-"), lw=1.3)


def bead_word(runs):
    kinds = {r.sticky for r in runs}
    return "Sticky bead" if kinds == {True} else "Pushing bead" if kinds == {False} else "Sticky and pushing beads"


def write_potential_comparison(sticky, push, path):
    """The sticky bead against the bead that only pushes, run by run: one row
    per reservoir, colour by fluidity, solid = sticky, dashed = pushing."""
    reservoirs = [k for k in RESERVOIR_STYLES if k != "fixed" and any(r.reservoir == k for r in sticky + push)]
    panels = [
        ("force on the bead, z (pN)", lambda r: running_mean(r.force[:, 2], 5)),
        ("vertices on the bead (bonded / in contact)", lambda r: r.bonded.astype(float)),
        ("excess area A − A_p (nm²)", lambda r: r.excess_area),
        ("neck radius (nm)", lambda r: r.neck),
        ("apex height above the sheet (nm)", lambda r: r.apex),
    ]
    fig, axes = plt.subplots(len(reservoirs), len(panels), figsize=(4.0 * len(panels), 3.4 * len(reservoirs)), dpi=110,
                             squeeze=False)
    for i, res in enumerate(reservoirs):
        for j, (label, pick) in enumerate(panels):
            ax = axes[i, j]
            shade_motion(ax, sticky + push)
            for group, ls, lw in ((sticky, "-", 1.4), (push, "--", 1.0)):
                for r in group:
                    if r.reservoir != res:
                        continue
                    ax.plot(r.time, pick(r), color=FLUIDITY_COLORS.get(r.fluidity, "#333333"), ls=ls, lw=lw,
                            label=f"{FLUIDITY_LABELS.get(r.fluidity, r.fluidity)}, {'sticky' if r.sticky else 'pushing'}")
            if j == 0:
                ax.axhline(0.0, color=INK_2, lw=0.6)
                ax.text(0.02, 0.96, RESERVOIR_LABELS[res], transform=ax.transAxes, va="top", fontsize=9, color=INK)
            if i == len(reservoirs) - 1:
                ax.set_xlabel("time (µs)", fontsize=9)
            if i == 0:
                ax.set_title(label, fontsize=9)
            style(ax)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.005))
    fig.suptitle("Sticky bead (solid) against the bead that only pushes (dashed): grey while the bead moves", fontsize=11)
    fig.tight_layout(rect=(0, 0.07, 1, 0.97))
    fig.savefig(path)
    plt.close(fig)


def shade_motion(ax, runs):
    run = max(runs, key=lambda r: r.time[-1])
    if run.stop_iteration >= 0:
        ax.axvspan(0.0, run.t_stop, color="#f0efec", zorder=0, lw=0)


def write_traces(runs, path):
    fig, axes = plt.subplots(3, 3, figsize=(15.0, 10.5), dpi=110)
    panels = [
        ("force on the bead, z (pN)", lambda r: running_mean(r.force[:, 2], 5)),
        ("lateral force on the bead (pN)", lambda r: running_mean(np.hypot(r.force[:, 0], r.force[:, 1]), 5)),
        (f"wrapped fraction ({runs[0].contact_word} vertices / sphere cells)", lambda r: r.wrapped),
        ("tension (pN/nm)", lambda r: running_mean(r.tension, 5)),
        ("excess area A − A_p (nm²)", lambda r: r.excess_area),
        ("neck radius (nm)", lambda r: r.neck),
        ("bending energy (pN·nm)", lambda r: r.e_bend),
        ("apex height above the sheet (nm)", lambda r: r.apex),
        ("accepted flips per µs", lambda r: running_mean(r.flips_per_us, 10)),
    ]
    for ax, (label, pick) in zip(axes.flat, panels):
        shade_motion(ax, runs)
        for r in runs:
            ax.plot(r.time, pick(r), label=r.label, **line_style(r))
        ax.set_ylabel(label, fontsize=9)
        ax.set_xlabel("time (µs)", fontsize=9)
        style(ax)
    axes[0, 0].axhline(0.0, color=INK_2, lw=0.6)
    axes[1, 0].axhline(0.5, color=INK_2, lw=0.6, ls=":")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=min(4, len(runs)), frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, -0.005))
    fig.suptitle(f"{bead_word(runs)} through a periodic membrane: grey while the bead moves, white while it holds",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    fig.savefig(path)
    plt.close(fig)


def write_hold_relaxation(runs, path, stats):
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.4), dpi=110)
    labels = ["force on the bead, z (pN)", f"{runs[0].contact_word} vertices", "tension (pN/nm)"]
    for r in runs:
        hold = r.hold_frames(skip_us=0.0)
        if len(hold) < 10:
            continue
        t = r.time[hold] - r.t_stop
        series = [running_mean(r.force[:, 2], 5)[hold], r.bonded[hold].astype(float), running_mean(r.tension, 5)[hold]]
        for ax, y in zip(axes, series):
            ax.plot(t, y, label=r.label, **line_style(r))
            fit = fit_relaxation(r.time[hold], y)
            if fit is not None:
                yi, a, tau = fit
                ax.plot(t, yi + a * np.exp(-t / tau), color=line_style(r)["color"], lw=0.7, alpha=0.5)
    for ax, label in zip(axes, labels):
        ax.set_xlabel("time after the bead stops (µs)")
        ax.set_ylabel(label)
        style(ax)
    axes[0].legend(frameon=False, fontsize=7.5)
    fig.suptitle("The hold phase: thin lines are the exponential fits y∞ + a·exp(−t/τ); τ in stats.md", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def frame_index(run, t_us):
    return int(np.clip(np.searchsorted(run.time, t_us), 0, run.n_frames - 1))


def energy_map_grid(runs, path, t_us, vmax=None):
    reservoirs = [k for k in RESERVOIR_STYLES if any(r.reservoir == k for r in runs)]
    fluidities = [k for k in FLUIDITY_COLORS if any(r.fluidity == k for r in runs)]
    fig, axes = plt.subplots(len(reservoirs), len(fluidities),
                             figsize=(4.0 * len(fluidities) + 2.0, 3.8 * len(reservoirs) + 0.8), dpi=110, squeeze=False)
    if vmax is None:
        vmax = max(np.percentile(r.bending_density(frame_index(r, t_us)), 99) for r in runs)
    norm = Normalize(0.0, vmax)
    for i, res in enumerate(reservoirs):
        for j, flu in enumerate(fluidities):
            ax = axes[i, j]
            match = [r for r in runs if r.reservoir == res and r.fluidity == flu]
            ax.set_facecolor("#f4f3f0")
            if not match:
                ax.text(0.5, 0.5, "no run", ha="center", va="center", transform=ax.transAxes, color=INK_2)
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            r = match[0]
            k = frame_index(r, t_us)
            pts = r.control[k]
            faces, dens = r.energy_faces(k)
            L = r.L0 * r.scale[k]
            # The margin shows the tile's own periodic continuation: the
            # energy-carrying faces shifted by the period. The mesh's copy
            # faces are not drawn -- see energy_faces().
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    ax.tripcolor(pts[:, 0] + dx * L[0], pts[:, 1] + dy * L[1], faces, facecolors=dens,
                                 cmap=ENERGY_CMAP, norm=norm, edgecolors="none")
            c = r.bead_center[k]
            ax.add_patch(plt.Circle((c[0], c[1]), r.radius, fill=False, ec="white", lw=1.0, ls="--"))
            o = r.origin * r.scale[k] - 0.5 * r.l_face
            ax.add_patch(plt.Rectangle((o[0], o[1]), L[0], L[1], fill=False, ec="#2bb3a5", lw=0.8))
            ax.set_xlim(o[0] - 0.08 * L[0], o[0] + 1.08 * L[0])
            ax.set_ylim(o[1] - 0.08 * L[1], o[1] + 1.08 * L[1])
            ax.set_aspect("equal")
            ax.set_title(f"{r.label}\nbead at {c[2] - r.base_level:+.0f} nm, {r.bonded[k]} {r.contact_word}", fontsize=8)
            ax.tick_params(labelsize=7)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=ENERGY_CMAP)
    fig.subplots_adjust(left=0.05, right=0.86, top=0.92, bottom=0.05, hspace=0.4, wspace=0.3)
    cax = fig.add_axes([0.89, 0.3, 0.015, 0.4])
    fig.colorbar(sm, cax=cax).set_label("bending energy density (pN/nm)")
    fig.suptitle(f"Bending energy density on the control net, top view, t = {t_us:.0f} µs "
                 f"(teal: the periodic tile; white: the bead)", fontsize=10)
    fig.savefig(path)
    plt.close(fig)


def write_profiles(runs, path, times):
    fig, axes = plt.subplots(1, len(times), figsize=(5.0 * len(times), 4.2), dpi=110, squeeze=False)
    for ax, t_us in zip(axes[0], times):
        for r in runs:
            k = frame_index(r, t_us)
            rr, prof = ring_profile(r, k)
            ax.plot(rr, prof, marker="o", ms=3, label=r.label, **line_style(r))
        ax.axvline(runs[0].radius, color=INK_2, lw=0.6, ls=":")
        ax.set_xlabel("distance from the bead's axis (nm)")
        ax.set_ylabel("bending energy density (pN/nm)")
        ax.set_title(f"t = {t_us:.0f} µs", fontsize=10)
        ax.set_yscale("log")
        style(ax)
    axes[0, 0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def write_spectra(runs, path, stats):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0), dpi=110)
    for r in runs:
        frames = r.hold_frames(skip_us=20.0)
        if len(frames) < 10:
            continue
        qs, st, fl, fit, A = spectrum(r, frames)
        axes[0].plot(qs, st, marker="o", ms=3, label=r.label, **line_style(r))
        axes[1].plot(qs, fl, marker="o", ms=3, label=r.label, **line_style(r))
        row = next((s for s in stats if s["run"] == r.name), None)
        if fit is not None:
            kc, sig = fit
            qq = np.linspace(qs[0], qs[-1], 100)
            axes[1].plot(qq, r.kbt / (A * (kc * qq ** 4 + sig * qq ** 2)), color=line_style(r)["color"], lw=0.7, alpha=0.5)
            if row is not None:
                row["spectrum_kc_pNnm"] = float(kc)
                row["spectrum_sigma_pN_nm"] = float(sig)
    q_ref = np.linspace(0.06, 0.6, 50)
    axes[1].plot(q_ref, runs[0].kbt / (runs[0].L0.prod() * runs[0].kc * q_ref ** 4), color=INK_2, lw=0.8, ls="--",
                 label=f"kT/(A kc q⁴), kc = {runs[0].kc:g}, σ = 0")
    axes[0].set_title("static imprint of the bead, |⟨h_q⟩|² over the hold", fontsize=10)
    axes[1].set_title("fluctuations of the sheet outside the bead's footprint (1.5 R), thin lines: kc, σ fit", fontsize=10)
    for ax in axes:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("q (1/nm)")
        ax.set_ylabel("⟨|h_q|²⟩ (nm²), h_q = (1/N) Σ h e^{-iq·r}")
        style(ax)
    axes[1].legend(frameon=False, fontsize=7)
    fig.suptitle("Fourier analysis of the upper envelope of the height field on the tile's lattice grid, hold phase (first 20 µs skipped)", fontsize=9.5)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def draw_bead(ax, centre, radius):
    u, v = np.mgrid[0:2 * np.pi:28j, 0:np.pi:14j]
    ax.plot_surface(centre[0] + radius * np.sin(v) * np.cos(u), centre[1] + radius * np.sin(v) * np.sin(u),
                    centre[2] + radius * np.cos(v), color=BEAD_COLOR, alpha=0.55, linewidth=0, shade=True)


def paint_membrane(ax, run, k, norm, zmax, xlim, zoom=1.3, elev=22, azim=-58):
    ax.cla()
    pts = run.surface[k]
    faces, values = run.energy_faces(k)
    surface = ax.plot_trisurf(pts[:, 0], pts[:, 1], pts[:, 2], triangles=faces, cmap=ENERGY_CMAP,
                              linewidth=0.12, edgecolor=(0.0, 0.0, 0.0, 0.10), antialiased=True)
    surface.set_array(values)
    surface.set_norm(norm)
    draw_bead(ax, run.bead_center[k], run.radius)
    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(-xlim, xlim)
    zlo, zhi = run.base_level - 0.3 * zmax, run.base_level + 1.1 * zmax
    ax.set_zlim(zlo, zhi)
    ax.set_box_aspect((1.0, 1.0, (zhi - zlo) / (2.0 * xlim)), zoom=zoom)
    ax.view_init(elev=elev, azim=azim)
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


def write_animation(runs, path, n_frames, fps, vmax):
    n = min(r.n_frames for r in runs)
    picked = np.unique(np.linspace(0, n - 1, min(n_frames, n)).round().astype(int))
    hold = int(round(0.8 * fps))
    sequence = [picked[0]] * hold + list(picked) + [picked[-1]] * (2 * hold)
    norm = Normalize(0.0, vmax)
    zmax = max(np.nanmax(r.apex) for r in runs) + runs[0].radius
    xlim = 0.5 * max(r.L0.max() for r in runs) * 1.15
    m = len(runs)
    width = 4.6 * m + 4.6
    fig = plt.figure(figsize=(width, 5.4), dpi=100)
    frac = 4.6 / width
    ax3d = [fig.add_axes([0.005 + frac * i, 0.02, frac - 0.01, 0.9], projection="3d") for i in range(m)]
    ax2d = fig.add_axes([frac * m + 0.07, 0.24, 1.0 - frac * m - 0.09, 0.62])
    shade_motion(ax2d, runs)
    for r in runs:
        ax2d.plot(r.time, running_mean(r.force[:, 2], 5), label=r.label, **line_style(r))
    ax2d.axhline(0.0, color=INK_2, lw=0.6)
    ax2d.set_xlim(0.0, max(r.time[-1] for r in runs))
    ax2d.set_xlabel("time (µs)")
    ax2d.set_ylabel("force on the bead, z (pN)")
    ax2d.legend(loc="upper right", frameon=False, fontsize=7.5)
    style(ax2d)
    cursor = ax2d.axvline(runs[0].time[0], color=INK_2, lw=0.9)
    titles = [fig.text(0.02 + frac * i, 0.95, "", fontsize=10, color=INK) for i in range(m)]
    stamp = fig.text(frac * m + 0.05, 0.95, "", fontsize=11, color=INK)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=ENERGY_CMAP)
    cax = fig.add_axes([frac * m + 0.07, 0.07, 1.0 - frac * m - 0.09, 0.025])
    fig.colorbar(sm, cax=cax, orientation="horizontal").set_label("bending energy density (pN/nm)", fontsize=8)

    def update(k):
        for ax, run, title in zip(ax3d, runs, titles):
            paint_membrane(ax, run, k, norm, zmax, xlim)
            title.set_text(f"{run.label}: bead at {run.bead_center[k, 2] - run.base_level:+.0f} nm, "
                           f"{run.bonded[k]} {run.contact_word}")
        cursor.set_xdata([runs[0].time[k], runs[0].time[k]])
        stamp.set_text(f"t = {runs[0].time[k]:.1f} µs")
        return tuple(ax3d) + (cursor, stamp)

    anim = animation.FuncAnimation(fig, update, frames=sequence, interval=1000.0 / fps, blit=False)
    anim.save(path, writer=animation.PillowWriter(fps=fps), dpi=100)
    plt.close(fig)
    return len(sequence)


def write_stats(stats, md_path, csv_path):
    keys = list(stats[0].keys())
    for s in stats:
        for k in s:
            if k not in keys:
                keys.append(k)
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for s in stats:
            writer.writerow({k: s.get(k, "") for k in keys})
    cols = [("run", "run"), ("potential", "bead"), ("bonded_end", "on the bead"), ("wrapped_end", "wrapped"),
            ("Fz_peak_pN", "Fz peak (pN)"),
            ("Fz_at_stop_pN", "Fz at stop"), ("Fz_end_pN", "Fz end"), ("Flat_end_pN", "|F_lat| end"),
            ("tension_end_pN_nm", "σ end"), ("excess_area_end_nm2", "A − A_p (nm²)"), ("reservoir_area_end_nm2", "ΔA reservoir"),
            ("neck_end_nm", "neck (nm)"), ("apex_end_nm", "apex (nm)"), ("E_bend_end_pNnm", "E_bend end"),
            ("flips_per_us_end", "flips/µs"), ("tau_force_us", "τ_F (µs)"), ("tau_bonded_us", "τ_bond (µs)"),
            ("tau_tension_us", "τ_σ (µs)"), ("spectrum_kc_pNnm", "kc fit"), ("spectrum_sigma_pN_nm", "σ fit"),
            ("min_altitude_nm", "min alt (nm)"), ("creases_over_90", "creases>90°")]
    with open(md_path, "w") as handle:
        handle.write("| " + " | ".join(c[1] for c in cols) + " |\n")
        handle.write("|" + "---|" * len(cols) + "\n")
        for s in stats:
            cells = []
            for key, _ in cols:
                v = s.get(key, "")
                if isinstance(v, float):
                    cells.append("nan" if not np.isfinite(v) else (f"{v:.3g}" if abs(v) < 1000 else f"{v:.0f}"))
                else:
                    cells.append(str(v))
            handle.write("| " + " | ".join(cells) + " |\n")


# ---------------------------------------------------------------------------
def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--runs", nargs="*", help="run directory names under --root (default: all with a bead log)")
    parser.add_argument("--out", default=None)
    parser.add_argument("--frames", type=int, default=150)
    parser.add_argument("--fps", type=int, default=15)
    parser.add_argument("--no-gif", action="store_true")
    parser.add_argument("--compare", nargs="*", help="runs for the side-by-side compare.gif (up to 4)")
    parser.add_argument("--map-times", nargs="*", type=float, default=None, help="times (us) of the energy maps")
    parser.add_argument("--compare-root", default=None,
                        help="a second matrix (the sticky bead) to set against --root in potential_comparison.png")
    args = parser.parse_args(argv)
    root = os.path.expanduser(args.root)
    names = args.runs or sorted(os.path.basename(d) for d in glob.glob(os.path.join(root, "*"))
                                if os.path.exists(os.path.join(d, "inputbead.csv")))
    out = args.out or os.path.join(root, "analysis")
    os.makedirs(out, exist_ok=True)
    runs = []
    for name in names:
        try:
            runs.append(StickyRun(os.path.join(root, name)))
            print(f"loaded {name}: {runs[-1].n_frames} frames to {runs[-1].time[-1]:.1f} us, {runs[-1].n_free} free vertices, tile {runs[-1].L0.round(1)}")
        except Exception as exc:
            print(f"skipping {name}: {exc}")
    if not runs:
        sys.exit("no runs")
    runs.sort(key=lambda r: (list(RESERVOIR_STYLES).index(r.reservoir), list(FLUIDITY_COLORS).index(r.fluidity) if r.fluidity in FLUIDITY_COLORS else 9))

    stats = [summarize(r) for r in runs]
    write_traces(runs, os.path.join(out, "traces.png"))
    write_hold_relaxation(runs, os.path.join(out, "hold_relaxation.png"), stats)
    t_end = min(r.time[-1] for r in runs)
    times = args.map_times or sorted({min(runs[0].t_stop, t_end), t_end})
    vmax = max(np.percentile(r.bending_density(frame_index(r, t_end)), 99) for r in runs)
    for t_us in times:
        energy_map_grid(runs, os.path.join(out, f"energy_maps_{t_us:.0f}us.png"), t_us, vmax=vmax)
    write_profiles(runs, os.path.join(out, "energy_profiles.png"), times)
    write_spectra(runs, os.path.join(out, "spectra.png"), stats)
    write_stats(stats, os.path.join(out, "stats.md"), os.path.join(out, "stats.csv"))
    print(open(os.path.join(out, "stats.md")).read())
    if args.compare_root:
        other_root = os.path.expanduser(args.compare_root)
        other = []
        for name in sorted(os.listdir(other_root)):
            if os.path.exists(os.path.join(other_root, name, "inputbead.csv")):
                try:
                    other.append(StickyRun(os.path.join(other_root, name)))
                except Exception as exc:
                    print(f"skipping {name} of {other_root}: {exc}")
        sticky = [r for r in runs + other if r.sticky]
        push = [r for r in runs + other if not r.sticky]
        write_potential_comparison(sticky, push, os.path.join(out, "potential_comparison.png"))
        both = sorted(stats + [summarize(r) for r in other], key=lambda s: (s["reservoir"], s["fluidity"], s["potential"]))
        write_stats(both, os.path.join(out, "stats_comparison.md"), os.path.join(out, "stats_comparison.csv"))
        print("potential_comparison.png and stats_comparison.md written")
    if not args.no_gif:
        for r in runs:
            n = write_animation([r], os.path.join(out, f"anim_{r.name}.gif"), args.frames, args.fps, vmax)
            print(f"anim_{r.name}.gif: {n} frames")
        if args.compare:
            chosen = [r for r in runs if r.name in args.compare][:4]
            if chosen:
                write_animation(chosen, os.path.join(out, "compare.gif"), args.frames, args.fps, vmax)
                print("compare.gif written")


if __name__ == "__main__":
    main()
