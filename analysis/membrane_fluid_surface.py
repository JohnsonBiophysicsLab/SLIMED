"""The Monge height of a Loop surface whose connectivity is not a lattice.

Why this module exists
----------------------
`membrane_resample.LimitSurface` inverts the regular 12-node Loop patch to read
the limit surface at a chosen (x, y).  It is exact, and it is exact *only* for
a mesh in which every interior vertex has valence 6 -- which the flat sheet
has right up to the first accepted edge flip.  A fluid membrane is mostly not
valence 6 (`docs/fluidity_results.md` section 4: about 56% of interior vertices
after 3000 steps), and every face touching an extraordinary vertex carries an
irregular patch the resampler does not implement.

This module reads the surface a different way, one that needs nothing but the
subdivision rules and works for any valence:

1. **Subdivide** the control net `k` times on its actual connectivity, with the
   Loop rules SLIMED itself uses -- `include/mesh/Subdivision_matrices.hpp` and
   `Mesh::refine_loop_once()`: vertex weight `beta = 3/(8n)` (Warren's form),
   edge rule `3/8, 3/8, 1/8, 1/8`, crease rules `3/4, 1/8, 1/8` and the midpoint
   on the sheet's open outer boundary.
2. **Push to the limit.**  With Warren's `beta`, Loop's limit-position mask is
   `1/2` on the vertex and `1/(2n)` on each neighbour at *every* valence -- the
   same mask `slimed::SurfaceSolver` applies -- so the subdivided vertices land
   exactly on the limit surface, not merely near it.
3. **Interpolate linearly** between those limit points onto the Cartesian Monge
   grid.  After `k` levels the triangles have edge `l/2^k`, and the linear
   interpolation error is of order `(q l / 2^k)^2 / 8`: at the top of the fitting
   window, `q l = pi/2`, that is 0.3% for `k = 2` and 0.08% for `k = 3`.

   The interpolation triangulates the *projected* limit points afresh
   (Delaunay) rather than reusing the mesh's own faces.  A fluid control net
   folds in projection -- one to three per cent of its triangles are inverted
   in the xy plane within a few thousand steps, which is in-plane motion doing
   what it does -- and although the limit net inherits almost none of that,
   a single inverted facet is enough to break a point locator built on the
   mesh connectivity.  A Delaunay triangulation of the same points is always
   valid, coincides with the mesh triangulation wherever the mesh is well
   shaped, and regularises across the rare fold, which is the only sensible
   reading of a Monge height there.  `fold_counts()` reports how often it
   happens; ``method="mesh"`` keeps the mesh-connectivity locator for checks
   on the lattice.

Every step is a linear map of the control net, and the maps depend only on the
connectivity, so for each distinct connectivity in a run one sparse matrix is
assembled once -- `L = Limit(F_k) S_k ... S_1` -- and each frame costs one
sparse product plus a point location on the subdivided triangulation.  The
point location does have to be redone per frame, because the vertices move in
plane and the triangles move with them.

What this is checked against
----------------------------
* On a frozen, all-regular mesh the route must reproduce the exact resampler
  to the interpolation error above, and the fitted `kc` to well under a per
  cent -- `test_membrane_fluid_surface.py` and the notebook both do this.
* The limit position of a vertex is invariant under subdivision, so the limit
  points this module computes for the *original* vertices must equal the
  limit positions SLIMED writes to `surfacepoint*.csv`, on a flipped mesh as
  much as on the lattice.  That is an exact cross-check against the C++ on the
  connectivity the resampler cannot handle.
* A planar control net must come back as that plane: every step here
  reproduces affine functions.

One caution about the sheet itself.  Periodicity is realised by duplicated
rings, and a flip is refused wherever it would touch one, so the duplicate
rings keep the lattice connectivity while their partners across the box may
have flipped.  The limit surface is therefore not exactly periodic near the
seam -- a property of the simulation, not of this reader.  `seam_mismatch()`
measures it.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay

try:
    import matplotlib.tri as mtri
except ImportError:  # pragma: no cover
    mtri = None


# --------------------------------------------------------------------------
# connectivity
# --------------------------------------------------------------------------
class Connectivity:
    """Edges, incidence and fans of a triangle mesh, from its face list alone."""

    def __init__(self, faces: np.ndarray, n_vertices: int | None = None):
        F = np.asarray(faces, dtype=np.int64)
        if F.ndim != 2 or F.shape[1] != 3:
            raise ValueError(f"faces must be (M, 3), got {F.shape}")
        self.faces = F
        self.n_vertices = int(F.max()) + 1 if n_vertices is None else int(n_vertices)

        # Directed half-edges (a, b) with the third corner c, one per face side.
        a = F.ravel()
        b = np.roll(F, -1, axis=1).ravel()
        c = np.roll(F, -2, axis=1).ravel()
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        key = lo * self.n_vertices + hi
        uniq, inverse, counts = np.unique(key, return_inverse=True, return_counts=True)
        self.n_edges = uniq.size
        self.edge_lo = uniq // self.n_vertices
        self.edge_hi = uniq % self.n_vertices
        self.edge_face_count = counts
        #: id of the undirected edge for each face side, shape (M, 3): side k is (F[:,k], F[:,k+1]).
        self.face_edge = inverse.reshape(F.shape)
        if counts.max() > 2:
            raise ValueError("an edge is shared by more than two faces; not a two-manifold")

        # The (up to two) opposite corners of each edge, -1 where absent.
        opp = np.full((self.n_edges, 2), -1, dtype=np.int64)
        order = np.argsort(inverse, kind="stable")
        e_sorted, c_sorted = inverse[order], c[order]
        first = np.ones(e_sorted.size, bool)
        first[1:] = e_sorted[1:] != e_sorted[:-1]
        opp[e_sorted[first], 0] = c_sorted[first]
        second = ~first
        opp[e_sorted[second], 1] = c_sorted[second]
        self.edge_opposite = opp

        # Vertex valence (distinct neighbours) and incident-face count.
        self.valence = np.bincount(np.concatenate([self.edge_lo, self.edge_hi]),
                                   minlength=self.n_vertices)
        self.incident_faces = np.bincount(F.ravel(), minlength=self.n_vertices)
        # A vertex on an open boundary touches an edge with one face.
        boundary_edge = counts == 1
        self.on_boundary = np.zeros(self.n_vertices, bool)
        self.on_boundary[self.edge_lo[boundary_edge]] = True
        self.on_boundary[self.edge_hi[boundary_edge]] = True
        #: A closed fan: as many incident faces as neighbours, and not on the boundary.
        self.interior = (~self.on_boundary) & (self.incident_faces == self.valence) \
            & (self.valence > 0)

    # -- neighbour structure as a sparse adjacency ---------------------------
    def adjacency(self) -> sp.csr_matrix:
        rows = np.concatenate([self.edge_lo, self.edge_hi])
        cols = np.concatenate([self.edge_hi, self.edge_lo])
        return sp.csr_matrix((np.ones(rows.size), (rows, cols)),
                             shape=(self.n_vertices, self.n_vertices))

    def boundary_adjacency(self) -> sp.csr_matrix:
        b = self.edge_face_count == 1
        rows = np.concatenate([self.edge_lo[b], self.edge_hi[b]])
        cols = np.concatenate([self.edge_hi[b], self.edge_lo[b]])
        return sp.csr_matrix((np.ones(rows.size), (rows, cols)),
                             shape=(self.n_vertices, self.n_vertices))


# --------------------------------------------------------------------------
# one Loop subdivision step, as a matrix
# --------------------------------------------------------------------------
def loop_subdivision_matrix(conn: Connectivity):
    """The linear map from a control net to its once-subdivided net.

    Returns ``(S, faces_new)``: ``S`` is ``(N + E) x N`` sparse, and the new
    vertex ``N + e`` is the odd (edge) vertex of edge ``e``.  The rules are the
    ones `Mesh::refine_loop_once()` applies, so a mesh refined by SLIMED and
    one refined here are the same mesh.
    """
    N, E = conn.n_vertices, conn.n_edges
    rows, cols, vals = [], [], []

    # --- even vertices ------------------------------------------------------
    A = conn.adjacency().tocoo()
    n = conn.valence.astype(float)

    interior = conn.interior
    beta = np.where(n > 0, 3.0 / (8.0 * np.maximum(n, 1)), 0.0)      # Warren
    # self weight 1 - n beta
    rows.append(np.arange(N)[interior]); cols.append(np.arange(N)[interior])
    vals.append((1.0 - n * beta)[interior])
    # neighbours, beta each
    keep = interior[A.row]
    rows.append(A.row[keep]); cols.append(A.col[keep]); vals.append(beta[A.row[keep]])

    # Boundary: crease rule where exactly two boundary neighbours, else pinned.
    B = conn.boundary_adjacency()
    n_b = np.asarray(B.sum(axis=1)).ravel()
    crease = conn.on_boundary & (n_b == 2)
    pinned = (~interior) & (~crease)
    rows.append(np.arange(N)[crease]); cols.append(np.arange(N)[crease])
    vals.append(np.full(int(crease.sum()), 0.75))
    Bc = B.tocoo()
    keep = crease[Bc.row]
    rows.append(Bc.row[keep]); cols.append(Bc.col[keep])
    vals.append(np.full(int(keep.sum()), 0.125))
    rows.append(np.arange(N)[pinned]); cols.append(np.arange(N)[pinned])
    vals.append(np.ones(int(pinned.sum())))

    # --- odd vertices -------------------------------------------------------
    e = np.arange(E)
    two_faces = conn.edge_face_count == 2
    w_end = np.where(two_faces, 3.0 / 8.0, 0.5)
    rows.append(N + e); cols.append(conn.edge_lo); vals.append(w_end)
    rows.append(N + e); cols.append(conn.edge_hi); vals.append(w_end)
    for k in range(2):
        opp = conn.edge_opposite[:, k]
        ok = two_faces & (opp >= 0)
        rows.append(N + e[ok]); cols.append(opp[ok]); vals.append(np.full(int(ok.sum()), 0.125))

    S = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(N + E, N))

    # --- faces: each (a, b, c) -> (a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)
    F = conn.faces
    ab = N + conn.face_edge[:, 0]
    bc = N + conn.face_edge[:, 1]
    ca = N + conn.face_edge[:, 2]
    a, b, c = F[:, 0], F[:, 1], F[:, 2]
    faces_new = np.concatenate([
        np.stack([a, ab, ca], axis=1),
        np.stack([ab, b, bc], axis=1),
        np.stack([ca, bc, c], axis=1),
        np.stack([ab, bc, ca], axis=1),
    ])
    return S, faces_new


def limit_mask_matrix(conn: Connectivity) -> sp.csr_matrix:
    """Loop's limit-position mask for Warren's weights: 1/2 self, 1/(2n) each neighbour.

    The identity on vertices without a closed fan, which is also what
    `slimed::SurfaceSolver` does with them.
    """
    N = conn.n_vertices
    A = conn.adjacency().tocoo()
    interior = conn.interior
    n = conn.valence.astype(float)
    rows = [np.arange(N)]
    cols = [np.arange(N)]
    vals = [np.where(interior, 0.5, 1.0)]
    keep = interior[A.row]
    rows.append(A.row[keep]); cols.append(A.col[keep])
    vals.append(0.5 / n[A.row[keep]])
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(N, N))


# --------------------------------------------------------------------------
# the surface of one connectivity
# --------------------------------------------------------------------------
class SubdividedSurface:
    """Everything about one connectivity that does not depend on positions.

    ``limit_map`` takes the original control net (N x 3) to the limit points of
    the ``k``-times-subdivided net; ``faces`` is that net's triangulation.
    """

    def __init__(self, faces: np.ndarray, n_vertices: int, levels: int = 2):
        if levels < 0:
            raise ValueError("levels must be >= 0")
        self.levels = levels
        conn = Connectivity(faces, n_vertices)
        self.base = conn
        total = sp.identity(n_vertices, format="csr")
        for _ in range(levels):
            S, faces = loop_subdivision_matrix(conn)
            total = S @ total
            conn = Connectivity(faces, total.shape[0])
        self.fine = conn
        self.faces = conn.faces
        self.limit_map = (limit_mask_matrix(conn) @ total).tocsr()
        self.n_fine = total.shape[0]

    def limit_points(self, control: np.ndarray) -> np.ndarray:
        """Limit points of the fine net, ``(n_fine, 3)``, from the control net ``(N, 3)``."""
        return self.limit_map @ np.asarray(control, float)

    def original_vertex_limits(self, control: np.ndarray) -> np.ndarray:
        """Limit points of the *original* vertices, which persist as the first N of the fine net."""
        return self.limit_points(control)[: self.base.n_vertices]

    def interpolate(self, control: np.ndarray, X: np.ndarray, Y: np.ndarray,
                    method: str = "delaunay") -> np.ndarray:
        """Monge height ``h(X, Y)`` by linear interpolation on the fine limit net.

        ``method="delaunay"`` (default) triangulates the projected limit points
        afresh, which survives a folded mesh; ``"mesh"`` locates points in the
        mesh's own faces and raises on a fold.  ``nan`` where a grid point falls
        outside the mesh; on the periodic sheet the tile is several rings
        inside the outer boundary so that never happens, and the notebook
        asserts it.
        """
        P = self.limit_points(control)
        X, Y = np.asarray(X, float), np.asarray(Y, float)
        if method == "mesh":
            if mtri is None:
                raise RuntimeError("matplotlib is needed for mesh point location")
            tri = mtri.Triangulation(P[:, 0], P[:, 1], self.faces)
            return np.asarray(mtri.LinearTriInterpolator(tri, P[:, 2])(X, Y).filled(np.nan))
        if method != "delaunay":
            raise ValueError(f"unknown method {method!r}")
        interp = LinearNDInterpolator(Delaunay(P[:, :2]), P[:, 2], fill_value=np.nan)
        return interp(np.column_stack([X.ravel(), Y.ravel()])).reshape(X.shape)

    def fold_counts(self, control: np.ndarray) -> dict:
        """How many projected triangles are inverted, in the control net and the fine net.

        A triangle is inverted when its signed xy area has the opposite sign to
        the mesh's majority.  On a graph over the plane the count is zero; on a
        fluid control net it is not, and the limit net inherits a little of it.
        """
        def inverted(P, F):
            a, b, c = P[F[:, 0]], P[F[:, 1]], P[F[:, 2]]
            area = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
            return int((area * np.sign(np.median(area)) <= 0).sum()), len(F)

        P = np.asarray(control, float)
        coarse = inverted(P, self.base.faces)
        fine = inverted(self.limit_points(P), self.faces)
        return {"control_inverted": coarse[0], "control_faces": coarse[1],
                "fine_inverted": fine[0], "fine_faces": fine[1]}


# --------------------------------------------------------------------------
# a whole run
# --------------------------------------------------------------------------
def connectivity_for_frame(face_frames: dict, iteration: int) -> np.ndarray:
    """The face list in force at ``iteration``: the most recent one written at or before it."""
    usable = [it for it in face_frames if it <= iteration]
    return face_frames[max(usable)] if usable else face_frames[min(face_frames)]


def monge_heights(coords: np.ndarray, frame_iterations: np.ndarray, face_frames: dict,
                  X: np.ndarray, Y: np.ndarray, levels: int = 2, verbose: bool = False):
    """Resample every frame of a run onto the grid ``(X, Y)``.

    ``coords`` is ``(n_frames, N, 3)`` -- the *whole* mesh, ghosts included, as
    `fluidity.load_run` returns it.  Returns ``(n_frames, *X.shape)``.  One
    `SubdividedSurface` is built per distinct connectivity and reused.
    """
    n_frames, n_vertices, _ = coords.shape
    cache: dict[int, SubdividedSurface] = {}
    out = np.empty((n_frames,) + X.shape)
    for f in range(n_frames):
        it = int(frame_iterations[f])
        usable = [k for k in face_frames if k <= it]
        key = max(usable) if usable else min(face_frames)
        surf = cache.get(key)
        if surf is None:
            surf = cache[key] = SubdividedSurface(face_frames[key], n_vertices, levels)
            if verbose:
                print(f"  connectivity at iteration {key}: {surf.n_fine} fine vertices, "
                      f"{len(surf.faces)} faces")
        out[f] = surf.interpolate(coords[f], X, Y)
    return out


def seam_mismatch(coords: np.ndarray, frame_iterations, face_frames: dict, lat, ox: float,
                  oy: float, levels: int = 2, stride: int = 20) -> dict:
    """How far from periodic the surface is across the box seam.

    Evaluates ``h`` just inside each edge of the tile and at the same points
    shifted by the box, and returns the r.m.s. difference.  On a lattice this
    is round-off; on a fluid sheet it measures the ghost band's lattice
    connectivity disagreeing with its flipped partner.

    Compare it with the r.m.s. of the *field*, not with ``rms_along_seam``:
    the box's largest modes are constant along a line parallel to one axis, so
    a single edge line under-samples the variance by several times.
    ``rms_along_seam`` is returned only so that a reader can see that.
    """
    n = 16
    eps = 0.25 * lat.d_face_x
    ys = oy + np.linspace(0, lat.ly, n, endpoint=False)
    xs = ox + np.linspace(0, lat.lx, n, endpoint=False)
    X = np.concatenate([np.full(n, ox + eps), xs])
    Y = np.concatenate([ys, np.full(n, oy + eps)])
    Xs = np.concatenate([np.full(n, ox + lat.lx + eps), xs])
    Ys = np.concatenate([ys, np.full(n, oy + lat.ly + eps)])
    diffs, mags = [], []
    frames = range(0, coords.shape[0], max(1, stride))
    h_all = monge_heights(coords[list(frames)], frame_iterations[list(frames)], face_frames,
                          np.concatenate([X, Xs]), np.concatenate([Y, Ys]), levels)
    for h in h_all:
        a, b = h[: 2 * n], h[2 * n:]
        ok = np.isfinite(a) & np.isfinite(b)
        diffs.append(np.sqrt(np.mean((a[ok] - b[ok]) ** 2)))
        # About the frame's mean: the sheet sits at z = 10, and a mismatch
        # measured against that offset would look ten times smaller than it is.
        mags.append(np.sqrt(np.mean((a[ok] - a[ok].mean()) ** 2)))
    return {"rms_mismatch": float(np.mean(diffs)), "rms_along_seam": float(np.mean(mags)),
            "frames": len(diffs)}
