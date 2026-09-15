"""Resample the Loop limit surface onto a Monge-gauge (x, y) grid.

Why this module exists
----------------------
`membrane_spectrum.tile_fft2` transforms the height column of the vertex array
and supplies the (x, y) of each sample from the *ideal* lattice.  That is exact
only while the vertices stay where the lattice says they are, which today they
do: `DynamicModel::next_step()` builds a three-component displacement and then
zeroes two of them (`if (j < 2) { randomterm *= 0.0; forceterm *= 0.0; }`), so
x and y never move and the measured in-plane drift over a run is ~1e-10 nm,
pure round-off from the `surface2mesh * mesh2surface` round trip.

The moment anything moves the vertices in plane -- lateral mobility, a
Metropolis move that displaces all three components, a surface that stops being
a graph over its own mesh -- that assumption fails silently, because the height
reader never looks at the x, y columns at all.  A displacement d enters the
transform as a per-site phase exp(i q.d): a Debye-Waller-like loss of coherent
amplitude whose lost power is scattered diffusely over every mode.  The fitting
window barely notices, but the q^-4 tail sits four decades down and is
contaminated by whatever leaks out of the low-q modes -- the same failure mode
as dropping the zig-zag phase.

This module removes the assumption instead of documenting it.  Rather than
reading the height at whatever parameter value a vertex happens to occupy, it
inverts the surface: for a requested (x, y) it finds the patch and the
barycentric coordinate (v, w) whose limit-surface point projects there, and
evaluates z at it.  The field handed to the FFT is then a genuine Monge height
h(x, y) on a grid *we* choose, and the choice of grid is independent of where
the control points have drifted to.

The patch
---------
Every interior vertex of the flat mesh has valence 6, so every face carries the
regular 12-node Loop patch, whose basis is the quartic box spline.  The basis
functions and their v, w derivatives below are transcribed from
`get_shapefunction()` in src/mesh/Gauss_quadrature.cpp, so the surface this
module evaluates is the same surface SLIMED integrates -- not a re-derivation
that could drift from it.

The node layout was derived, not guessed.  Every box spline reproduces affine
functions exactly, so solving

    sum_i N_i(v, w) r_i = affine(v, w)

for the 12 unknown positions r_i has a unique answer, and it puts them on the
integer triangular lattice listed in `PATCH_SHEARED` with the three corners at
nodes 3 (u = 1), 6 (v = 1) and 7 (w = 1).  `test_membrane_resample.py` redoes
that solve and checks the corner masks, so a wrong ordering cannot pass
silently.

Indexing
--------
Two index systems appear here, and keeping them apart is most of the work.

*zig-zag* (i, j) is what the mesh and the trajectory use: row j sits at
y = j*dFaceY and its x origin is offset half a cell on alternating rows.
Periodicity is axis-aligned in these indices -- (i, j+ny) is the same vertex one
box up -- which is what makes the wrap-around trivial.

*sheared* (m, n) is what the patch stencil is naturally written in: position is
m*a1 + n*a2 with a1 = (dx, 0), a2 = (dx/2, dy), so the 12 nodes are a fixed
integer offset pattern with no parity cases.  The two are related by
i = m + shear(n) with shear a floor or ceiling of n/2 depending on which rows
carry the offset -- `LimitSurface._shear` picks the branch from the lattice
itself rather than assuming GHOST_RINGS is odd.

Faces are labelled by their sheared cell (m, n) plus an orientation: UP has
corners (m,n), (m+1,n), (m,n+1); DOWN has (m+1,n), (m,n+1), (m+1,n+1) and is the
same patch seen through the 180-degree rotation that is a symmetry of the
lattice, so it reuses one stencil rather than needing a second.
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass

import numpy as np

__all__ = [
    "shape_functions",
    "PATCH_SHEARED",
    "LimitSurface",
    "Resampler",
    "monge_grid",
    "cartesian_q_grid",
    "cartesian_spectrum",
    "read_tile_control",
    "reference_net",
]

UP, DOWN = 0, 1

#: The 12 regular-patch nodes as (m, n) offsets on the sheared integer lattice,
#: in the column order of `shape_functions`.  Nodes 3, 6, 7 are the corners
#: carrying u = 1, v = 1, w = 1.
PATCH_SHEARED = np.array(
    [(0, -1), (-1, 0), (1, -1), (0, 0), (-1, 1), (2, -1),
     (1, 0), (0, 1), (-1, 2), (2, 0), (1, 1), (0, 2)],
    dtype=np.int64,
)

#: Which basis column is the corner of each barycentric coordinate.
CORNER_U, CORNER_V, CORNER_W = 3, 6, 7


def shape_functions(v, w):
    """The 12 box-spline basis functions and their first derivatives.

    Returns ``(N, N_v, N_w)``, each of shape ``(..., 12)``, for barycentric
    coordinates ``v``, ``w`` with ``u = 1 - v - w``.  Transcribed from rows 0-2
    of `get_shapefunction()` in src/mesh/Gauss_quadrature.cpp.
    """
    v = np.asarray(v, dtype=float)
    w = np.asarray(w, dtype=float)
    u = 1.0 - v - w
    zero = np.zeros(np.broadcast(v, w).shape)

    def _stack(*cols):
        return np.stack([zero + c for c in cols], axis=-1)

    N = _stack(
        (1.0 / 12.0) * (u**4 + 2.0 * u**3 * v),
        (1.0 / 12.0) * (u**4 + 2.0 * u**3 * w),
        (1.0 / 12.0) * (u**4 + 2.0 * u**3 * w + 6.0 * u**3 * v + 6.0 * u**2 * v * w + 12.0 *
            u**2 * v**2 + 6.0 * u * v**2 * w + 6.0 * u * v**3 + 2.0 * v**3 * w + v**4),
        (1.0 / 12.0) * (6.0 * u**4 + 24.0 * u**3 * w + 24.0 * u**2 * w**2 + 8.0 * u * w**3 +
            w**4 + 24.0 * u**3 * v + 60.0 * u**2 * v * w + 36.0 * u * v * w**2 + 6.0 * v * w**3
            + 24.0 * u**2 * v**2 + 36.0 * u * v**2 * w + 12.0 * v**2 * w**2 + 8.0 * u * v**3 +
            6.0 * v**3 * w + v**4),
        (1.0 / 12.0) * (u**4 + 6.0 * u**3 * w + 12.0 * u**2 * w**2 + 6.0 * u * w**3 + w**4 + 2.0
            * u**3 * v + 6.0 * u**2 * v * w + 6.0 * u * v * w**2 + 2.0 * v * w**3),
        (1.0 / 12.0) * (2.0 * u * v**3 + v**4),
        (1.0 / 12.0) * (u**4 + 6.0 * u**3 * w + 12.0 * u**2 * w**2 + 6.0 * u * w**3 + w**4 + 8.0
            * u**3 * v + 36.0 * u**2 * v * w + 36.0 * u * v * w**2 + 8.0 * v * w**3 + 24.0 *
            u**2 * v**2 + 60.0 * u * v**2 * w + 24.0 * v**2 * w**2 + 24.0 * u * v**3 + 24.0 *
            v**3 * w + 6.0 * v**4),
        (1.0 / 12.0) * (u**4 + 8.0 * u**3 * w + 24.0 * u**2 * w**2 + 24.0 * u * w**3 + 6.0 *
            w**4 + 6.0 * u**3 * v + 36.0 * u**2 * v * w + 60.0 * u * v * w**2 + 24.0 * v * w**3
            + 12.0 * u**2 * v**2 + 36.0 * u * v**2 * w + 24.0 * v**2 * w**2 + 6.0 * u * v**3 +
            8.0 * v**3 * w + v**4),
        (1.0 / 12.0) * (2.0 * u * w**3 + w**4),
        (1.0 / 12.0) * (2.0 * v**3 * w + v**4),
        (1.0 / 12.0) * (2.0 * u * w**3 + w**4 + 6.0 * u * v * w**2 + 6.0 * v * w**3 + 6.0 * u *
            v**2 * w + 12.0 * v**2 * w**2 + 2.0 * u * v**3 + 6.0 * v**3 * w + v**4),
        (1.0 / 12.0) * (w**4 + 2.0 * v * w**3),
    )
    N_v = _stack(
        (1.0 / 12.0) * (-2.0 * u**3 - 6.0 * u**2 * v),
        (1.0 / 12.0) * (-4.0 * u**3 - 6.0 * u**2 * w),
        (1.0 / 12.0) * (2.0 * u**3 + 6.0 * u**2 * v - 6.0 * u * v**2 - 2.0 * v**3),
        (1.0 / 12.0) * (-12.0 * u**2 * w - 12.0 * u * w**2 - 2.0 * w**3 - 24.0 * u**2 * v - 48.0
            * u * v * w - 12.0 * v * w**2 - 24.0 * u * v**2 - 18.0 * v**2 * w - 4.0 * v**3),
        (1.0 / 12.0) * (-2.0 * u**3 - 12.0 * u**2 * w - 18.0 * u * w**2 - 4.0 * w**3 - 6.0 *
            u**2 * v - 12.0 * u * v * w - 6.0 * v * w**2),
        (1.0 / 12.0) * (6.0 * u * v**2 + 2.0 * v**3),
        (1.0 / 12.0) * (4.0 * u**3 + 18.0 * u**2 * w + 12.0 * u * w**2 + 2.0 * w**3 + 24.0 *
            u**2 * v + 48.0 * u * v * w + 12.0 * v * w**2 + 24.0 * u * v**2 + 12.0 * v**2 * w),
        (1.0 / 12.0) * (2.0 * u**3 + 12.0 * u**2 * w + 12.0 * u * w**2 + 6.0 * u**2 * v - 12.0 *
            v * w**2 - 6.0 * u * v**2 - 12.0 * v**2 * w - 2.0 * v**3),
        -(1.0 / 6.0) * w**3,
        (1.0 / 12.0) * (6.0 * v**2 * w + 4.0 * v**3),
        (1.0 / 12.0) * (4.0 * w**3 + 18.0 * v * w**2 + 6.0 * u * w**2 + 12.0 * v**2 * w + 12.0 *
            u * v * w + 2.0 * v**3 + 6.0 * u * v**2),
        (1.0 / 6.0) * w**3,
    )
    N_w = _stack(
        (1.0 / 12.0) * (-4.0 * u**3 - 6.0 * u**2 * v),
        (1.0 / 12.0) * (-2.0 * u**3 - 6.0 * u**2 * w),
        (1.0 / 12.0) * (-2.0 * u**3 - 6.0 * u**2 * w - 12.0 * u**2 * v - 12.0 * u * v * w - 18.0
            * u * v**2 - 6.0 * v**2 * w - 4.0 * v**3),
        (1.0 / 12.0) * (-24.0 * u**2 * w - 24.0 * u * w**2 - 4.0 * w**3 - 12.0 * u**2 * v - 48.0
            * u * v * w - 18.0 * v * w**2 - 12.0 * u * v**2 - 12.0 * v**2 * w - 2.0 * v**3),
        (1.0 / 12.0) * (2.0 * u**3 + 6.0 * u**2 * w - 6.0 * u * w**2 - 2.0 * w**3),
        -(1.0 / 6.0) * v**3,
        (1.0 / 12.0) * (2.0 * u**3 + 6.0 * u**2 * w - 6.0 * u * w**2 - 2.0 * w**3 + 12.0 * u**2
            * v - 12.0 * v * w**2 + 12.0 * u * v**2 - 12.0 * v**2 * w),
        (1.0 / 12.0) * (4.0 * u**3 + 24.0 * u**2 * w + 24.0 * u * w**2 + 18.0 * u**2 * v + 48.0
            * u * v * w + 12.0 * v * w**2 + 12.0 * u * v**2 + 12.0 * v**2 * w + 2.0 * v**3),
        (1.0 / 12.0) * (6.0 * u * w**2 + 2.0 * w**3),
        (1.0 / 6.0) * v**3,
        (1.0 / 12.0) * (2.0 * w**3 + 6.0 * u * w**2 + 12.0 * v * w**2 + 12.0 * u * v * w + 18.0
            * v**2 * w + 6.0 * u * v**2 + 4.0 * v**3),
        (1.0 / 12.0) * (4.0 * w**3 + 6.0 * v * w**2),
    )

    return N, N_v, N_w


# --------------------------------------------------------------------------
# the resampling operator
# --------------------------------------------------------------------------
@dataclass
class Resampler:
    """A fixed linear map from tile control heights to grid heights.

    Once the geometry is solved -- which patch and which (v, w) each grid point
    landed on -- the height is just ``sum_k N_k z_k`` over the 12 control
    points of that patch.  The parameter solve depends only on the *in-plane*
    control positions, so while those are frozen it is done once and reused for
    every frame, which turns the whole resampling of a trajectory into one
    gather plus a weighted sum.

    ``node_index`` holds flat indices into the ``(ny*nx)`` tile, already wrapped
    periodically, so a patch that straddles the box edge needs no special case.
    """

    node_index: np.ndarray      #: (npts, 12) int, flat index into the tile
    weights: np.ndarray         #: (npts, 12) the basis values N_k(v, w)
    grid_shape: tuple           #: shape to fold the points back into
    residual: float             #: max |S_xy(v,w) - (x,y)| over the grid, nm
    iterations: int             #: Newton sweeps used
    faces: tuple                #: (m, n, kind) per point, for inspection
    vw: tuple                   #: (v, w) per point, for inspection
    converged: np.ndarray       #: (npts,) bool -- False where no (v, w) was found

    def apply(self, z_tile: np.ndarray) -> np.ndarray:
        """Resample control heights ``(..., ny, nx)`` onto the grid.

        Accumulates node by node rather than gathering all 12 at once: for a
        long trajectory the (frames, points, 12) intermediate is the largest
        array in the calculation and it is not needed.
        """
        z = np.asarray(z_tile, float)
        flat = z.reshape(z.shape[:-2] + (-1,))
        out = np.zeros(flat.shape[:-1] + (self.node_index.shape[0],))
        for k in range(12):
            out += flat[..., self.node_index[:, k]] * self.weights[:, k]
        return out.reshape(out.shape[:-1] + self.grid_shape)


# --------------------------------------------------------------------------
# the surface
# --------------------------------------------------------------------------
class LimitSurface:
    """The Loop limit surface of one periodic tile, as a function of (v, w).

    ``ctrl`` is the tile's control net, ``(ny, nx, 3)``, in the index
    convention of `membrane_spectrum.read_tile_heights` -- entry ``[j, i]`` is
    mesh vertex ``(GHOST_RINGS + i, GHOST_RINGS + j)``.  Only the in-plane
    columns are used for the geometry; heights come later, through
    :meth:`Resampler.apply`.
    """

    def __init__(self, lat, ctrl):
        self.lat = lat
        self.ctrl = np.asarray(ctrl, float)
        if self.ctrl.shape != (lat.ny, lat.nx, 3):
            raise ValueError(f"expected control net {(lat.ny, lat.nx, 3)}, "
                             f"got {self.ctrl.shape}")
        self._ctrl_flat = self.ctrl.reshape(-1, 3)
        # Which parity of tile row carries the half-cell offset.  The tile
        # starts at GHOST_RINGS, so this is not the same parity as the mesh's.
        self._row0_offset = bool(lat.row_is_offset()[0])

        # Origin of the ideal lattice: the position tile vertex (0, 0) would
        # have with its own offset removed, so that x = ox + i*dx + offset.
        dx = lat.d_face_x
        self._ox = float(self.ctrl[0, 0, 0]) - (dx / 2 if self._row0_offset else 0.0)
        self._oy = float(self.ctrl[0, 0, 1])

    # -- index algebra ----------------------------------------------------
    def _shear(self, n):
        """i = m + shear(n): floor(n/2) or ceil(n/2), whichever the tile uses."""
        return (n + (1 if self._row0_offset else 0)) // 2

    def _wrap(self, i, j):
        """Zig-zag index -> (flat tile index, periodic position shift).

        Periodicity is axis-aligned in zig-zag indices: ny is even, so j and
        j+ny have the same parity and therefore the same half-cell offset, and
        i is untouched by a step in y.  That is the whole reason the stencil is
        built in sheared indices but looked up in zig-zag ones.
        """
        lat = self.lat
        ky, jn = np.divmod(j, lat.ny)
        kx, im = np.divmod(i, lat.nx)
        return jn * lat.nx + im, kx * lat.lx, ky * lat.ly

    def face_nodes(self, m, n, kind):
        """The 12 control points of a face: (flat index, xy) each ``(npts, 12)``.

        DOWN is the same stencil read through the 180-degree rotation about the
        cell centre, which is a symmetry of the triangular lattice, so the
        offsets are simply negated about the opposite corner.
        """
        m, n, kind = np.asarray(m), np.asarray(n), np.asarray(kind)
        up = (kind == UP)[..., None]
        base_m = np.where(up, m[..., None], m[..., None] + 1)
        base_n = np.where(up, n[..., None], n[..., None] + 1)
        sgn = np.where(up, 1, -1)

        M = base_m + sgn * PATCH_SHEARED[:, 0]
        N = base_n + sgn * PATCH_SHEARED[:, 1]
        i, j = M + self._shear(N), N

        flat, sx, sy = self._wrap(i, j)
        xy = self._ctrl_flat[flat][..., :2] + np.stack([sx, sy], axis=-1)
        return flat, xy

    def evaluate(self, m, n, kind, v, w):
        """The limit-surface point S(v, w) on a given face, as ``(..., 3)``.

        Mostly for checking the thing works: two faces sharing an edge must
        agree along it, which is what pins the DOWN stencil down independently
        of the UP one it is derived from.
        """
        flat, xy = self.face_nodes(m, n, kind)
        P = np.concatenate([xy, self._ctrl_flat[flat][..., 2:]], axis=-1)
        return np.einsum("...k,...kd->...d", shape_functions(v, w)[0], P)

    # -- locating a target point -----------------------------------------
    def _seed_ideal(self, X, Y):
        """Which face of the *ideal* lattice covers (X, Y), and where in it.

        A seed, not an answer: it is exact while the control net sits on its
        lattice, and merely close once the net deforms, which is what the walk
        in :meth:`solve` is for.
        """
        lat = self.lat
        n_c = (Y - self._oy) / lat.d_face_y
        m_c = (X - self._ox - n_c * (lat.d_face_x / 2)) / lat.d_face_x
        m0, n0 = np.floor(m_c).astype(np.int64), np.floor(n_c).astype(np.int64)
        fm, fn = m_c - m0, n_c - n0

        lower = (fm + fn) <= 1.0
        kind = np.where(lower, UP, DOWN)
        v = np.where(lower, fm, 1.0 - fm)
        w = np.where(lower, fn, 1.0 - fn)
        return m0, n0, kind, v, w

    def _corner_xy(self, m, n, kind, which):
        """In-plane positions of a face's three corners, ordered (u, v, w).

        ``which="mesh"`` takes the control points themselves; ``which="limit"``
        takes the limit points, which is where the patch's parameter corners
        actually are and so is the better affine model of the map being
        inverted.
        """
        _, xy = self.face_nodes(m, n, kind)
        if which == "mesh":
            return xy[:, [CORNER_U, CORNER_V, CORNER_W], :]
        corners = []
        for vv, ww in ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)):
            nc = shape_functions(np.array(vv), np.array(ww))[0]
            corners.append(np.einsum("k,pkd->pd", nc, xy))
        return np.stack(corners, axis=1)

    def _flat_barycentric(self, m, n, kind, X, Y, which):
        """(v, w) of (X, Y) in the affine triangle through a face's corners."""
        c = self._corner_xy(m, n, kind, which)
        e1 = c[:, 1] - c[:, 0]
        e2 = c[:, 2] - c[:, 0]
        d = np.stack([X, Y], axis=-1) - c[:, 0]
        det = e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]
        v = (d[:, 0] * e2[:, 1] - d[:, 1] * e2[:, 0]) / det
        w = (e1[:, 0] * d[:, 1] - e1[:, 1] * d[:, 0]) / det
        return v, w

    def _newton(self, m, n, kind, X, Y, v, w, tol, maxiter):
        """Solve S_xy(v, w) = (X, Y) on a fixed face.

        The map is a quartic in (v, w), so there is no closed-form inverse to
        reach for -- but its Jacobian is the same box spline differentiated,
        which the shape functions already carry, and Newton from the affine
        seed converges quadratically.
        """
        _, P = self.face_nodes(m, n, kind)
        target = np.stack([X, Y], axis=-1)
        used = 0
        for used in range(1, maxiter + 1):
            N, N_v, N_w = shape_functions(v, w)
            F = np.einsum("pk,pkd->pd", N, P) - target
            if np.max(np.abs(F)) < tol:
                break
            Jv = np.einsum("pk,pkd->pd", N_v, P)
            Jw = np.einsum("pk,pkd->pd", N_w, P)
            det = Jv[:, 0] * Jw[:, 1] - Jv[:, 1] * Jw[:, 0]
            v = v - (Jw[:, 1] * F[:, 0] - Jw[:, 0] * F[:, 1]) / det
            w = w - (Jv[:, 0] * F[:, 1] - Jv[:, 1] * F[:, 0]) / det
        return v, w, used

    def _hop(self, m, n, kind, v, w, move):
        """Step to the face across the edge whose barycentric coordinate went negative."""
        k = np.argmin(np.stack([1.0 - v - w, v, w], axis=-1), axis=-1)
        up = kind == UP
        dm = np.where(k == 1, np.where(up, -1, 1), 0)
        dn = np.where(k == 2, np.where(up, -1, 1), 0)
        m = np.where(move, m + dm, m)
        n = np.where(move, n + dn, n)
        kind = np.where(move, 1 - kind, kind)
        return m, n, kind

    # -- the public entry point ------------------------------------------
    def solve(self, X, Y, method="newton", tol=1e-11, max_newton=30, max_walk=8):
        """Build the :class:`Resampler` that reads the surface at (X, Y).

        ``method``:

        ``"newton"``
            Invert S_xy(v, w) = (x, y) exactly, to ``tol``.  This is the answer.
        ``"limit"``
            One affine solve in the triangle through the patch's three *limit*
            points -- no iteration, and the natural first-order version of the
            idea that the barycentric coordinate on the surface is close to the
            one on a flat triangle.
        ``"mesh"``
            The same, but through the three *control* points.  Cheapest, and
            the least accurate: the control triangle is not on the surface.

        The Newton path walks between faces when a point turns out to belong to
        a neighbour, so the seed only has to be close, not right.
        """
        X = np.asarray(X, float)
        grid_shape = X.shape
        Xf, Yf = X.ravel(), np.asarray(Y, float).ravel()

        m, n, kind, v, w = self._seed_ideal(Xf, Yf)
        iterations = 0
        converged = np.ones(Xf.shape, bool)

        if method in ("limit", "mesh"):
            v, w = self._flat_barycentric(m, n, kind, Xf, Yf,
                                          "limit" if method == "limit" else "mesh")
            # No walk here, so a point the affine model puts outside its seed
            # face is being extrapolated.  Say so rather than reporting success;
            # `residual` then shows how far off the answer actually is.
            converged = ~((v < -1e-12) | (w < -1e-12) | (v + w > 1.0 + 1e-12))
        elif method == "newton":
            out = np.zeros(Xf.shape, bool)
            for _ in range(max_walk):
                v, w, used = self._newton(m, n, kind, Xf, Yf, v, w, tol, max_newton)
                iterations += used
                out = (v < -1e-12) | (w < -1e-12) | (v + w > 1.0 + 1e-12)
                if not out.any():
                    break
                m, n, kind = self._hop(m, n, kind, v, w, out)
                v2, w2 = self._flat_barycentric(m, n, kind, Xf, Yf, "limit")
                v, w = np.where(out, v2, v), np.where(out, w2, w)
            if out.any():
                # Usually not a solver failure but a statement about the mesh:
                # where the projected triangles have folded, the surface is not
                # a graph over the plane and no (v, w) maps to that (x, y) at
                # all.  Clamp those points into their last face so the rest of
                # the grid is still usable, and say so rather than returning a
                # silently wrong number -- `inverted_control_triangles` counts
                # the folds if that is the cause.
                converged = ~out
                v = np.clip(v, 0.0, 1.0)
                w = np.clip(w, 0.0, 1.0 - v)
                warnings.warn(
                    f"{int(out.sum())} of {out.size} grid points found no (v, w) "
                    f"in {max_walk} walks; {self.inverted_control_triangles()} "
                    f"control triangles are folded in projection. Those points "
                    f"are clamped to a patch edge and flagged in .converged.",
                    RuntimeWarning, stacklevel=2)
        else:
            raise ValueError(f"unknown method {method!r}")

        flat, P = self.face_nodes(m, n, kind)
        N = shape_functions(v, w)[0]
        residual = float(np.max(np.abs(
            np.einsum("pk,pkd->pd", N, P) - np.stack([Xf, Yf], axis=-1))))
        return Resampler(flat, N, grid_shape, residual, iterations,
                         (m, n, kind), (v, w), converged)

    def face_corners(self, m, n, kind, which="mesh"):
        """The three corner positions of each face, ordered (u, v, w).

        ``which="mesh"`` gives the control points; ``"limit"`` gives the limit
        points, which is where the patch's parameter corners actually sit.
        """
        return self._corner_xy(np.asarray(m), np.asarray(n), np.asarray(kind), which)

    def crosses_seam(self, m, n, kind):
        """Whether each face's 12-node stencil wraps around the periodic box.

        A non-periodic reference field -- a tilted plane, say -- is only a valid
        check on faces that do not wrap, so tests need to know which those are.
        """
        flat, xy = self.face_nodes(m, n, kind)
        raw = self._ctrl_flat[flat][..., :2]
        return np.any(np.abs(xy - raw) > 1e-9, axis=(-2, -1))

    def inverted_control_triangles(self) -> int:
        """How many faces have folded over in projection.

        The Monge gauge presumes the surface is a graph over the (x, y) plane.
        A face whose projected control triangle has flipped orientation is a
        fold, and no resampling can be well defined across it -- so this is the
        first thing to check when :meth:`solve` reports points it could not
        place.  Zero on any mesh that has merely been deformed gently.
        """
        lat = self.lat
        m, n = np.meshgrid(np.arange(lat.nx), np.arange(lat.ny), indexing="xy")
        m, n = m.ravel(), n.ravel()
        folded = 0
        for kind in (UP, DOWN):
            c = self._corner_xy(m, n, np.full(m.shape, kind), "mesh")
            e1, e2 = c[:, 1] - c[:, 0], c[:, 2] - c[:, 0]
            folded += int((e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0] <= 0).sum())
        return folded

    # -- grids ------------------------------------------------------------
    def monge_grid(self, kind="cartesian", factor=1):
        """The (x, y) grid to resample onto, as two ``(ny*f, nx*f)`` arrays.

        ``"zigzag"`` is the vertex set itself -- the control points' own (x, y),
        not an idealisation of them -- so that the resampled field can be
        compared with the vertex-based one point by point.  On a frozen lattice
        the two must agree exactly; once the vertices move in plane they stop
        agreeing, and the difference is the error the vertex-based transform
        was making.

        ``"cartesian"`` is the grid the Monge gauge actually wants: a plain
        rectangular lattice over the same box, which a bare ``np.fft.fft2``
        transforms with no half-cell phase to remember.
        """
        lat = self.lat
        if kind == "zigzag":
            if factor != 1:
                raise ValueError("the zig-zag grid is the vertex set; factor must be 1")
            return self.ctrl[..., 0].copy(), self.ctrl[..., 1].copy()
        if kind != "cartesian":
            raise ValueError(f"unknown grid {kind!r}")
        nx, ny = lat.nx * factor, lat.ny * factor
        dx, dy = lat.d_face_x / factor, lat.d_face_y / factor
        i, j = np.meshgrid(np.arange(nx), np.arange(ny), indexing="xy")
        return self._ox + i * dx, self._oy + j * dy


def cartesian_q_grid(lat, factor=1, fold_to_brillouin_zone=True):
    """Wavevectors for a spectrum taken on the cartesian Monge grid.

    Same box, so the raw bins are the same 2*pi*(n_x/Lx, n_y/Ly) that
    `membrane_spectrum.q_grid` returns -- but the *folding* differs, and that is
    the point of having a separate function.  A sampled mode is the whole family
    q + G, with G the reciprocal lattice of the grid that did the sampling: the
    triangular one for the vertex set, the rectangular one here.  Labelling a
    cartesian-sampled spectrum with triangular folding assigns the wrong |q| to
    the modes near the zone boundary, which is exactly where the two samplings
    differ most.

    Returns ``(qx, qy, qmag)``, each ``(ny*factor, nx*factor)``.
    """
    nx, ny = lat.nx * factor, lat.ny * factor
    dx, dy = lat.d_face_x / factor, lat.d_face_y / factor
    qx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    qy = 2.0 * np.pi * np.fft.fftfreq(ny, d=dy)
    qx, qy = np.meshgrid(qx, qy, indexing="xy")

    if fold_to_brillouin_zone:
        gx, gy = 2.0 * np.pi / dx, 2.0 * np.pi / dy
        bx, by = qx.copy(), qy.copy()
        best = qx**2 + qy**2
        for m1 in (-1, 0, 1):
            for m2 in (-1, 0, 1):
                if m1 == 0 and m2 == 0:
                    continue
                cx, cy = qx + m1 * gx, qy + m2 * gy
                c = cx**2 + cy**2
                take = c < best - 1e-12
                bx, by = np.where(take, cx, bx), np.where(take, cy, by)
                best = np.where(take, c, best)
        qx, qy = bx, by
    return qx, qy, np.sqrt(qx**2 + qy**2)


def cartesian_spectrum(h_grid, lat):
    """<|h_q|^2> in nm^4 from a field on the cartesian grid.

    A plain ``fft2`` is the whole transform here: the grid is rectangular, so
    there is no half-cell phase to put back.  That is the practical dividend of
    resampling -- the correction `tile_fft2` exists to apply stops being needed
    because the sampling positions are ours to choose.
    """
    h = np.asarray(h_grid)
    n = h.shape[-1] * h.shape[-2]
    return np.mean(np.abs(np.fft.fft2(h, axes=(-2, -1))) ** 2, axis=0) / (n * n)


def monge_grid(lat, ctrl, kind="cartesian", factor=1):
    """Convenience wrapper: :meth:`LimitSurface.monge_grid` without the object."""
    return LimitSurface(lat, ctrl).monge_grid(kind, factor)


# --------------------------------------------------------------------------
# reading the in-plane columns
# --------------------------------------------------------------------------
def read_tile_control(csv_path, lat, max_frames=None, start=0, stride=1):
    """Read the full ``(x, y, z)`` of the periodic tile from a trajectory CSV.

    `membrane_spectrum.read_tile_heights` deliberately reads only the z column,
    which is all the vertex-based transform needs.  Resampling needs to know
    where the vertices actually are, so this reads all three.

    Returns ``(n_frames, ny, nx, 3)``.
    """
    rows, cols = lat.tile_slice()
    wanted = np.array(
        [(j * lat.n_vert_x + i) * 3 + c
         for j in range(rows.start, rows.stop)
         for i in range(cols.start, cols.stop)
         for c in (0, 1, 2)],
        dtype=np.int64,
    )
    n_field = lat.n_vertices * 3
    frames = []
    truncated = 0
    with open(csv_path, "r", encoding="utf-8") as fh:
        for k, line in enumerate(fh):
            if k < start or (k - start) % stride or not line.strip():
                continue
            parts = line.rstrip().rstrip(",").split(",")
            if len(parts) != n_field:
                # A trajectory still being appended to ends in a half-written
                # line; skipping beats reshaping it into nonsense.
                truncated += 1
                continue
            frames.append(np.fromiter((float(parts[t]) for t in wanted),
                                      dtype=float, count=wanted.size))
            if max_frames is not None and len(frames) >= max_frames:
                break
    if truncated:
        print(f"[read_tile_control] skipped {truncated} incomplete line(s) "
              f"in {os.path.basename(csv_path)}")
    if not frames:
        raise ValueError(f"no frames read from {csv_path}")
    return np.asarray(frames).reshape(-1, lat.ny, lat.nx, 3)


def reference_net(csv_path, lat):
    """The control net of the first frame, ``(ny, nx, 3)`` -- the geometry to solve on."""
    return read_tile_control(csv_path, lat, max_frames=1)[0]
