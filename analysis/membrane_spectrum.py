"""Height-fluctuation spectroscopy for a periodic flat SLIMED membrane.

What this module is for
-----------------------
A flat fluid membrane in the Monge gauge, parameterised by its height h(x, y)
above the reference plane, has the Helfrich energy

    E = (1/2) * integral dA [ kc (laplacian h)^2 + sigma (grad h)^2 ]

which, for a periodic box of projected area A and the convention

    h(r) = sum_q h_q exp(i q.r),      h_q = (1/A) integral d2r h(r) exp(-i q.r)

diagonalises into E = (A/2) sum_q (kc q^4 + sigma q^2) |h_q|^2.  Equipartition
then gives the relation this module measures:

    <|h_q|^2> = kB T / ( A (kc q^4 + sigma q^2) )                          (1)

With no tension that is the pure q^-4 law: a straight line of slope -4 on a
log-log plot, whose intercept fixes kc.  See Deserno, "Fluid lipid membranes -
a primer", and Brannigan/Brown-style spectral analyses of membrane simulations.

Why the obvious implementation gets the wrong answer
----------------------------------------------------
Four things about a SLIMED trajectory break a naive `np.fft.fft2` of the z
column, and each of them is corrected here.

1. The vertices are not a Cartesian grid.  `Mesh::set_vertices_faces_flat()`
   lays out a *zig-zag triangular* lattice: row j sits at y = j*dFaceY, and its
   x origin is offset by half a cell on even j.  Transforming the (i, j) index
   array as though it were an (x, y) array therefore omits a per-row phase.
   That is not a small error: the omitted factor is a period-2 modulation in j,
   so it mixes every mode n_y with the mode n_y + Ny/2 (see
   `naive_index_fft2`).  Since the spectrum spans four decades, the low-q modes
   leak into the high-q bins and flatten the tail into a plateau.

2. Only part of the mesh is independent.  Under `boundaryType = Periodic` the
   outer three rings are ghosts and the fourth ring is a duplicate of the
   opposite edge (`DynamicMesh::get_relative_pt_periodic`), so the periodic
   tile is the (nFaceX-6) x (nFaceY-6) block starting at index (3, 3).
   Transforming the whole vertex array -- or a block trimmed by an arbitrary
   number of rows -- transforms something that is not periodic, and the
   resulting discontinuity leaks power across the whole spectrum.

3. |q| has to come from the reciprocal lattice, not from the array index.
   q = 2*pi*(n_x/Lx, n_y/Ly) with L a *length*, and with n taken signed
   (`np.fft.fftfreq`): index n_x > Nx/2 is the short wavevector -q, not a long
   one.  Folding "the second half is a mirror image" by slicing the first
   quarter of the array and labelling it with the raw row/column index assigns
   large |q| to modes that physically have small |q|.

4. The DFT needs the 1/N normalisation to be h_q in the sense of (1).
   `np.fft.fft2` returns H_q = sum_j h_j exp(-i q.r_j) = N h_q.

Only the limit surface is the membrane.  The Loop control points are not points
on the surface; `surfacepoint*.csv` already holds the limit positions, and
`limit_surface_from_control` reproduces them from `meshpoint*.csv` for a
cross-check.

Everything in here is numpy-only, so it runs wherever the simulation does.
"""

from __future__ import annotations

import math
import os
import re
from dataclasses import dataclass

import numpy as np

__all__ = [
    "Lattice",
    "lattice_from_params",
    "read_tile_heights",
    "limit_surface_from_control",
    "tile_fft2",
    "naive_index_fft2",
    "q_grid",
    "power_spectrum",
    "radial_average",
    "fit_kc_sigma",
    "mode_autocorrelation",
]


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Lattice:
    """The zig-zag triangular mesh, and the periodic tile inside it.

    Index convention: vertex (i, j) of the full mesh is row-major entry
    ``j * nVertX + i`` of a trajectory line, and sits at

        x = i*dFaceX + (dFaceX/2 if j is even else 0) - lxHalf
        y = j*dFaceY - lyHalf

    The independent (periodic) block is ``i, j in [GHOST_RINGS, GHOST_RINGS+N)``.
    """

    n_face_x: int
    n_face_y: int
    d_face_x: float
    d_face_y: float

    #: Rings of ghost vertices the periodic boundary condition keeps on each
    #: side; the ring just inside them duplicates the opposite edge.
    GHOST_RINGS = 3

    @property
    def n_vert_x(self) -> int:
        return self.n_face_x + 1

    @property
    def n_vert_y(self) -> int:
        return self.n_face_y + 1

    @property
    def n_vertices(self) -> int:
        return self.n_vert_x * self.n_vert_y

    @property
    def nx(self) -> int:
        """Independent columns in the periodic tile."""
        return self.n_face_x - 2 * self.GHOST_RINGS

    @property
    def ny(self) -> int:
        """Independent rows in the periodic tile.  Even, so the zig-zag closes."""
        return self.n_face_y - 2 * self.GHOST_RINGS

    @property
    def lx(self) -> float:
        return self.nx * self.d_face_x

    @property
    def ly(self) -> float:
        return self.ny * self.d_face_y

    @property
    def area(self) -> float:
        """Projected area of the periodic box -- the A in equation (1)."""
        return self.lx * self.ly

    def row_is_offset(self) -> np.ndarray:
        """Which tile rows carry the half-cell x offset, as a bool array (ny,)."""
        j = self.GHOST_RINGS + np.arange(self.ny)
        return (j % 2) == 0

    def tile_slice(self):
        g = self.GHOST_RINGS
        return (slice(g, g + self.ny), slice(g, g + self.nx))

    def describe(self) -> str:
        return (
            f"mesh {self.n_vert_x} x {self.n_vert_y} vertices "
            f"(nFace {self.n_face_x} x {self.n_face_y}), "
            f"cell {self.d_face_x:g} x {self.d_face_y:.6g} nm\n"
            f"periodic tile {self.nx} x {self.ny} = {self.nx * self.ny} vertices, "
            f"Lx = {self.lx:g} nm, Ly = {self.ly:.6g} nm, A = {self.area:.6g} nm^2"
        )


def lattice_from_params(path: str) -> Lattice:
    """Rebuild the lattice from a SLIMED ``.params`` file.

    Mirrors ``Mesh::set_axes_division_flat()`` exactly, including the "make
    nFaceY even" step -- an odd row count would not tile.
    """
    text = open(path, encoding="utf-8", errors="replace").read()

    def get(name: str) -> float:
        m = re.search(rf"^\s*{re.escape(name)}\s*=\s*([-+0-9.eE]+)", text, re.M)
        if m is None:
            raise KeyError(f"{name} not found in {path}")
        return float(m.group(1))

    side_x, side_y, l_face = get("sideX"), get("sideY"), get("lFace")
    n_face_x = round(side_x / l_face)
    d_face_x = side_x / n_face_x
    d_face_y = math.sqrt(3.0) / 2.0 * d_face_x
    n_face_y = round(side_y / d_face_y)
    if n_face_y % 2 == 1:
        n_face_y += 1
    return Lattice(n_face_x, n_face_y, d_face_x, d_face_y)


# --------------------------------------------------------------------------
# reading trajectories
# --------------------------------------------------------------------------
def read_tile_heights(
    csv_path: str,
    lat: Lattice,
    stride: int = 1,
    start: int = 0,
    max_frames: int | None = None,
    dtype=np.float64,
) -> np.ndarray:
    """Read the z column of the periodic tile from a trajectory CSV.

    Each line of ``surfacepoint*.csv`` / ``meshpoint*.csv`` is one frame,
    ``x,y,z,`` repeated over every vertex with a trailing comma.  Only the
    fields the tile needs are converted, which is roughly a sixth of the file.

    Returns an array of shape ``(n_frames, ny, nx)``.
    """
    rows, cols = lat.tile_slice()
    wanted = np.array(
        [(j * lat.n_vert_x + i) * 3 + 2
         for j in range(rows.start, rows.stop)
         for i in range(cols.start, cols.stop)],
        dtype=np.int64,
    )

    n_field = lat.n_vertices * 3
    frames = []
    truncated = 0
    with open(csv_path, "r", encoding="utf-8") as fh:
        for n, line in enumerate(fh):
            if n < start:
                continue
            if (n - start) % stride:
                continue
            if not line.strip():
                continue
            parts = line.rstrip().rstrip(",").split(",")
            if len(parts) != n_field:
                # A trajectory being read while the run is still appending ends
                # in a half-written line.  Skipping it beats crashing, and beats
                # silently reshaping a short row into nonsense.
                truncated += 1
                continue
            frames.append(np.fromiter((float(parts[k]) for k in wanted),
                                      dtype=dtype, count=wanted.size))
            if max_frames is not None and len(frames) >= max_frames:
                break
    if truncated:
        print(f"[read_tile_heights] skipped {truncated} incomplete line(s) "
              f"in {os.path.basename(csv_path)}")
    if not frames:
        raise ValueError(f"no frames read from {csv_path}")
    return np.asarray(frames).reshape(-1, lat.ny, lat.nx)


# --------------------------------------------------------------------------
# control points -> limit surface
# --------------------------------------------------------------------------
def limit_surface_from_control(z_tile: np.ndarray, lat: Lattice) -> np.ndarray:
    """Apply the Loop limit-position mask to control heights, with wrap-around.

    For a valence-6 vertex the Loop limit point is 1/2 of the vertex plus 1/12
    of each of its six neighbours -- the same mask
    ``DynamicMesh::assign_mesh2surface()`` builds, here evaluated with periodic
    indices instead of on the ghost ring, so that it is exactly translation
    invariant across the tile.

    The neighbour columns depend on the row's zig-zag parity: an offset row's
    partners in the rows above and below are columns (i, i+1); a non-offset
    row's are (i-1, i).
    """
    z = np.asarray(z_tile)
    same_row = np.roll(z, 1, axis=-1) + np.roll(z, -1, axis=-1)

    up, down = np.roll(z, -1, axis=-2), np.roll(z, 1, axis=-2)
    offset = lat.row_is_offset()  # (ny,)
    # A row's own parity decides which two columns of the adjacent rows it
    # touches, so build both variants and select per row.
    def cross(a):
        left = np.roll(a, 1, axis=-1)    # column i-1
        right = np.roll(a, -1, axis=-1)  # column i+1
        shape = [1] * a.ndim
        shape[-2] = offset.size
        m = offset.reshape(shape)
        return np.where(m, a + right, left + a)

    return 0.5 * z + (cross(up) + cross(down) + same_row) / 12.0


def limit_mask_symbol(lat: Lattice) -> np.ndarray:
    """Fourier symbol m(q) of the limit mask, on the tile's q grid.

    The mask is a convolution, so on a perfect lattice it multiplies each mode
    by a scalar: m(q) = 1/2 + (1/6)[cos(q.a1) + cos(q.a2) + cos(q.a3)], which
    runs from 1 at q = 0 down to about 1/4 at the zone corner.  Useful for
    separating "the surface really is smoother than the control net" from a
    genuine deviation of the spectrum.
    """
    e = np.zeros((lat.ny, lat.nx))
    e[0, 0] = 1.0
    return np.real(tile_fft2(limit_surface_from_control(e, lat), lat)
                   / tile_fft2(e, lat))


# --------------------------------------------------------------------------
# the transform
# --------------------------------------------------------------------------
def tile_fft2(z_tile: np.ndarray, lat: Lattice) -> np.ndarray:
    """Exact DFT of the tile on its true (zig-zag) positions.

    H(q) = sum_j h_j exp(-i q.r_j).  Factorising r_j = (p*dx + s_offset, s*dy)
    turns this into an FFT along x, a per-row phase, and an FFT along y; the
    phase is exactly what a plain ``fft2`` of the index array drops.

    Accepts ``(..., ny, nx)`` and transforms the last two axes.
    """
    z = np.asarray(z_tile)
    if z.shape[-2:] != (lat.ny, lat.nx):
        raise ValueError(f"expected (..., {lat.ny}, {lat.nx}), got {z.shape}")

    g = np.fft.fft(z, axis=-1)

    # q_x * (dFaceX/2) = pi * n_x / nx for the offset rows.
    #
    # n_x has to be the SIGNED mode number here, the one np.fft.fftfreq
    # reports, not the raw output index.  The row transform itself cannot tell
    # n_x from n_x + nx -- both give the same column of `g` -- but the
    # half-cell phase can, because exp(-i pi (n_x + nx)/nx) = -exp(-i pi n_x/nx).
    # Sampling a half-cell offset is precisely what resolves q_x from
    # q_x + 2*pi/dFaceX, and getting the sign wrong flips that resolution: the
    # whole negative-q_x half of the spectrum comes out shifted by half the
    # zone in n_y, which pairs long-wavelength modes with short-wavelength ones.
    n_x = np.fft.fftfreq(lat.nx) * lat.nx                   # -nx/2 .. nx/2-1
    half_cell = np.exp(-1j * np.pi * n_x / lat.nx)          # (nx,)
    phase = np.where(lat.row_is_offset()[:, None], half_cell[None, :], 1.0)
    return np.fft.fft(g * phase, axis=-2)


def tile_ifft2(h_fft: np.ndarray, lat: Lattice) -> np.ndarray:
    """Inverse of :func:`tile_fft2`."""
    g = np.fft.ifft(np.asarray(h_fft), axis=-2)
    n_x = np.fft.fftfreq(lat.nx) * lat.nx
    half_cell = np.exp(-1j * np.pi * n_x / lat.nx)
    phase = np.where(lat.row_is_offset()[:, None], half_cell[None, :], 1.0)
    return np.fft.ifft(g / phase, axis=-1)


def naive_index_fft2(z_tile: np.ndarray) -> np.ndarray:
    """The uncorrected transform, kept so the error can be shown rather than asserted."""
    return np.fft.fft2(np.asarray(z_tile), axes=(-2, -1))


def synthesize_heights(lat: Lattice, spectrum, n_frames: int, rng) -> np.ndarray:
    """Draw Gaussian height fields with a prescribed <|h_q|^2>.

    Filtering white noise through sqrt(N * S(q)) in the transform of
    :func:`tile_fft2` gives a real field whose spectrum is S by construction,
    and it needs no Hermitian bookkeeping: the transform of a real array
    already has whatever conjugate structure the zig-zag lattice implies, and
    the filter is real and symmetric under q -> -q, so it preserves it.

    ``spectrum`` is an ``(ny, nx)`` array of <|h_q|^2>; its q = 0 entry is
    ignored and returned as zero mean.
    """
    n = lat.nx * lat.ny
    gain = np.sqrt(np.asarray(spectrum, float) * n)
    gain = np.where(np.isfinite(gain), gain, 0.0)
    gain[0, 0] = 0.0

    white = rng.standard_normal((n_frames, lat.ny, lat.nx))
    return np.real(tile_ifft2(tile_fft2(white, lat) * gain, lat))


def q_grid(lat: Lattice, fold_to_brillouin_zone: bool = True):
    """Wavevectors of the tile's modes.

    Returns ``(qx, qy, qmag)``, each of shape ``(ny, nx)`` and indexed the same
    way as the FFT output.  ``np.fft.fftfreq`` supplies the signed mode number,
    so the upper half of each axis is correctly read as a negative wavevector.

    With ``fold_to_brillouin_zone`` each q is replaced by the shortest of its
    reciprocal-lattice images.  A discrete mode *is* the whole family q + G;
    the continuum formula (1) refers to the smallest member, and near the zone
    boundary the rectangular fundamental domain is not the shortest choice.
    """
    qx = 2.0 * np.pi * np.fft.fftfreq(lat.nx, d=lat.d_face_x)
    qy = 2.0 * np.pi * np.fft.fftfreq(lat.ny, d=lat.d_face_y)
    qx, qy = np.meshgrid(qx, qy, indexing="xy")

    if fold_to_brillouin_zone:
        # Reciprocal primitive vectors of the triangular lattice a1 = (dx, 0),
        # a2 = (dx/2, dy).
        dx, dy = lat.d_face_x, lat.d_face_y
        g1 = np.array([2 * np.pi / dx, -np.pi / dy])
        g2 = np.array([0.0, 2 * np.pi / dy])
        best_x, best_y = qx.copy(), qy.copy()
        best = qx**2 + qy**2
        for m1 in (-1, 0, 1):
            for m2 in (-1, 0, 1):
                if m1 == 0 and m2 == 0:
                    continue
                cx = qx + m1 * g1[0] + m2 * g2[0]
                cy = qy + m1 * g1[1] + m2 * g2[1]
                c = cx**2 + cy**2
                take = c < best - 1e-12
                best_x = np.where(take, cx, best_x)
                best_y = np.where(take, cy, best_y)
                best = np.where(take, c, best)
        qx, qy = best_x, best_y

    return qx, qy, np.sqrt(qx**2 + qy**2)


def power_spectrum(h_fft: np.ndarray, lat: Lattice) -> np.ndarray:
    """<|h_q|^2> in nm^4, averaged over frames.

    ``h_fft`` is the output of :func:`tile_fft2`, i.e. H_q = N h_q, so the mean
    square is divided by N^2 to land on the convention of equation (1).
    """
    n = lat.nx * lat.ny
    return np.mean(np.abs(h_fft) ** 2, axis=0) / (n * n)


# --------------------------------------------------------------------------
# reduction and fitting
# --------------------------------------------------------------------------
def radial_average(q: np.ndarray, s: np.ndarray, n_bins: int = 24, q_min=None, q_max=None):
    """Average the 2-D spectrum in logarithmically spaced |q| shells.

    Returns ``(q_centre, s_mean, s_sem, count)``, dropping the q = 0 mode and
    any empty bin.  Log bins because the data are log-spaced in both axes;
    linear bins put nearly every mode in the last few shells.
    """
    q, s = np.asarray(q).ravel(), np.asarray(s).ravel()
    keep = q > 0
    q, s = q[keep], s[keep]
    lo = q.min() * 0.999 if q_min is None else q_min
    hi = q.max() * 1.001 if q_max is None else q_max
    edges = np.logspace(np.log10(lo), np.log10(hi), n_bins + 1)
    idx = np.digitize(q, edges) - 1

    qc, sm, se, cnt = [], [], [], []
    for b in range(n_bins):
        m = idx == b
        k = int(m.sum())
        if k == 0:
            continue
        qc.append(np.exp(np.mean(np.log(q[m]))))
        sm.append(s[m].mean())
        se.append(s[m].std(ddof=1) / math.sqrt(k) if k > 1 else 0.0)
        cnt.append(k)
    return (np.array(qc), np.array(sm), np.array(se), np.array(cnt))


def fit_kc_sigma(q, s, area, kbt, fit_tension=True):
    """Fit equation (1) to a spectrum, returning ``(kc, sigma)`` in pN.nm, pN/nm.

    1/(A S) = (kc/kBT) q^4 + (sigma/kBT) q^2 is linear in the two unknowns, so
    this is one least-squares solve with no starting guess.  Weighted by q^-4
    so that every decade of q counts about equally instead of the fit being
    dominated by the large-q end where 1/(A S) is biggest.
    """
    q, s = np.asarray(q, float), np.asarray(s, float)
    y = 1.0 / (area * s)
    w = q ** -4.0
    cols = [q**4] if not fit_tension else [q**4, q**2]
    a = np.vstack([c * w for c in cols]).T
    coef, *_ = np.linalg.lstsq(a, y * w, rcond=None)
    kc = coef[0] * kbt
    sigma = coef[1] * kbt if fit_tension else 0.0
    return kc, sigma


def log_slope(q, s):
    """Least-squares slope of log10 S against log10 q -- the -4 being looked for."""
    q, s = np.asarray(q, float), np.asarray(s, float)
    m = (q > 0) & (s > 0)
    return float(np.polyfit(np.log10(q[m]), np.log10(s[m]), 1)[0])


# --------------------------------------------------------------------------
# dynamics
# --------------------------------------------------------------------------
def mode_autocorrelation(h_fft: np.ndarray, max_lag: int | None = None) -> np.ndarray:
    """Normalised <h_q(t) h_q*(t+tau)> for every mode, via Wiener-Khinchin.

    ``h_fft`` has shape ``(n_frames, ny, nx)``.  Returns ``(n_lags, ny, nx)``
    real, with lag 0 equal to 1.  The frame mean is removed first, so a mode
    that never decorrelates within the run reads as a slow decay rather than a
    constant offset.
    """
    x = np.asarray(h_fft)
    n = x.shape[0]
    x = x - x.mean(axis=0, keepdims=True)

    n_fft = 1 << int(np.ceil(np.log2(2 * n)))
    f = np.fft.fft(x, n=n_fft, axis=0)
    acf = np.fft.ifft(f * np.conj(f), axis=0)[:n]
    acf = np.real(acf) / (n - np.arange(n))[:, None, None]  # unbiased normalisation

    zero = acf[0].copy()
    zero[zero == 0] = np.nan
    acf = acf / zero
    return acf if max_lag is None else acf[: max_lag + 1]


def relaxation_rate(acf_mode: np.ndarray, dt: float, floor: float = 1.0 / math.e,
                    n_frames: int | None = None, max_tau_fraction: float = 0.05,
                    min_lags: float = 2.0):
    """Decay rate of one mode, from the integral of its autocorrelation.

    Integrating C(tau) down to ``floor`` and inverting is steadier than fitting
    an exponential to noisy tails, and for a genuine exponential the two agree.

    A trajectory can only time a mode that it resolves at both ends, and the
    two guards here are what keep the unresolved ones out of the answer.

    ``max_tau_fraction`` rejects a mode whose correlation time is more than
    that fraction of the record.  It matters more than it looks: subtracting
    the sample mean from a mode whose correlation time approaches the run
    length removes nearly all of its real signal, and what is left decays on
    the sampling timescale.  Reported uncritically that is a spuriously *fast*
    rate at small q -- the wrong way round, and enough to fake a shallower
    power law.  A fraction of 1/50 or so is a reasonable working value.

    ``min_lags`` rejects the opposite failure: a mode that decorrelates inside
    a few sampling intervals is measured by the output interval rather than by
    the physics, and Gamma(q) saturates at 1/dt.  Both show up as a bend
    towards q^3 for entirely artificial reasons, which is exactly the thing
    this measurement is trying to rule in or out.

    Returns ``nan`` when either guard fires, or when the correlation never
    reaches ``floor`` at all.
    """
    c = np.asarray(acf_mode, float)
    below = np.nonzero(c <= floor)[0]
    if below.size == 0:
        return math.nan
    k = int(below[0])
    if k == 0:
        return math.nan
    # trapezoid over [0, k], then the analytic tail of an exponential matched
    # at the crossing point.
    area = np.trapezoid(c[: k + 1], dx=dt)
    tau = area / (1.0 - c[k])
    if not (tau > 0):
        return math.nan
    if n_frames is not None and tau > max_tau_fraction * n_frames * dt:
        return math.nan
    if k < 2 or tau < min_lags * dt:
        return math.nan
    return 1.0 / tau
