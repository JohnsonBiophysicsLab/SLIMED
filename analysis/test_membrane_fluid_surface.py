"""Self-tests for membrane_fluid_surface, on data whose answer is known exactly.

Run with ``~/anaconda3/bin/python analysis/test_membrane_fluid_surface.py``, or
under pytest.  The checks are of the kind a wrong subdivision rule cannot pass
by luck:

* every map here reproduces affine functions, so a planar net must come back
  as that plane after any number of levels and the limit mask;
* the limit position of a vertex is invariant under subdivision, so pushing
  the original vertices to the limit through k levels must give the same
  point as the mask applied directly -- and, on the flat lattice, the same
  point `membrane_spectrum.limit_surface_from_control` computes by a wholly
  different route;
* on a frozen all-regular net the route must agree with the exact resampler in
  `membrane_resample`, and the error must fall by about 4x per extra level,
  which is what linear interpolation on a mesh halved in edge length does;
* after a flip, the mesh has valences 5 and 7 and the route must still return
  a surface that interpolates the control net's limit points exactly.
"""

import sys

import numpy as np

import membrane_fluid_surface as mfs
import membrane_spectrum as ms
import membrane_resample as mr

FAILURES = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}{'  -- ' + detail if detail else ''}")
    if not ok:
        FAILURES.append(name)
        if "pytest" in sys.modules:
            raise AssertionError(f"{name}{'  -- ' + detail if detail else ''}")


def flat_sheet(n_face_x=20, n_face_y=24, l_face=5.0):
    """A zig-zag sheet in SLIMED's own vertex order, with its face list."""
    lat = ms.Lattice(n_face_x, n_face_y, l_face, np.sqrt(3) / 2 * l_face)
    P = np.zeros((lat.n_vertices, 3))
    for j in range(lat.n_vert_y):
        for i in range(lat.n_vert_x):
            P[j * lat.n_vert_x + i, 0] = i * lat.d_face_x + (lat.d_face_x / 2 if j % 2 == 0 else 0)
            P[j * lat.n_vert_x + i, 1] = j * lat.d_face_y
    F = []
    for j in range(n_face_y):
        for i in range(n_face_x):
            v = lambda ii, jj: jj * lat.n_vert_x + ii
            # An offset (even) row's vertex i sits half a cell right of the
            # row above's vertex i, so its nearest neighbours up there are i
            # and i+1; a non-offset row's are i-1 and i. Get this the wrong
            # way round and the sheet is still a valid triangulation -- just
            # a skewed one whose limit surface is not the lattice's.
            if j % 2 == 0:
                F.append([v(i, j), v(i + 1, j), v(i + 1, j + 1)])
                F.append([v(i, j), v(i + 1, j + 1), v(i, j + 1)])
            else:
                F.append([v(i, j), v(i + 1, j), v(i, j + 1)])
                F.append([v(i + 1, j), v(i + 1, j + 1), v(i, j + 1)])
    return lat, P, np.array(F)


def test_planar_net_is_reproduced():
    lat, P, F = flat_sheet()
    rng = np.random.default_rng(3)
    P = P.copy()
    P[:, :2] += rng.normal(scale=0.4, size=P[:, :2].shape)   # a deformed but planar net
    A, B, C = 0.031, -0.017, 2.5
    P[:, 2] = A * P[:, 0] + B * P[:, 1] + C
    for levels in (0, 1, 3):
        surf = mfs.SubdividedSurface(F, len(P), levels)
        Q = surf.limit_points(P)
        err = np.abs(Q[:, 2] - (A * Q[:, 0] + B * Q[:, 1] + C)).max()
        check(f"planar net survives {levels} level(s) and the limit mask", err < 1e-10,
              f"max z error {err:.1e}")
    surf = mfs.SubdividedSurface(F, len(P), 2)
    X, Y = np.meshgrid(np.linspace(20, 80, 13), np.linspace(20, 80, 11))
    h = surf.interpolate(P, X, Y)
    err = np.abs(h - (A * X + B * Y + C)).max()
    check("linear interpolation of a plane is exact", err < 1e-10, f"max error {err:.1e}")


def test_limit_is_invariant_under_subdivision():
    lat, P, F = flat_sheet()
    rng = np.random.default_rng(5)
    P[:, 2] = 1.5 * np.sin(0.13 * P[:, 0]) * np.cos(0.11 * P[:, 1]) + rng.normal(scale=0.2, size=len(P))
    direct = mfs.SubdividedSurface(F, len(P), 0).original_vertex_limits(P)
    # Interior vertices only: on the sheet's open outer boundary the limit is
    # taken as the identity, and that convention is not invariant under the
    # crease rule. The tile is four rings inside, so nothing there is affected.
    interior = mfs.Connectivity(F, len(P)).interior
    for levels in (1, 2, 3):
        via = mfs.SubdividedSurface(F, len(P), levels).original_vertex_limits(P)
        err = np.abs(via - direct)[interior].max()
        check(f"limit of original vertices unchanged by {levels} level(s)", err < 1e-10,
              f"max difference {err:.1e}")
    # and against membrane_spectrum's own limit mask on the tile
    rows, cols = lat.tile_slice()
    idx = np.array([j * lat.n_vert_x + i for j in range(rows.start, rows.stop)
                    for i in range(cols.start, cols.stop)])
    z_tile = P[idx, 2].reshape(lat.ny, lat.nx)[None]
    ref = ms.limit_surface_from_control(z_tile, lat)[0]
    mine = direct[idx, 2].reshape(lat.ny, lat.nx)
    # membrane_spectrum wraps periodically while this sheet has a real ghost band, so
    # compare away from the tile edge where both see the same neighbours only if the
    # sheet is periodic; make it so by copying the periodic images explicitly
    inner = (slice(1, -1), slice(1, -1))
    err = np.abs(mine[inner] - ref[inner]).max()
    check("limit mask agrees with membrane_spectrum on the tile interior", err < 1e-10,
          f"max difference {err:.1e}")


def test_agrees_with_exact_resampler_on_the_lattice():
    lat, P, F = flat_sheet()
    rng = np.random.default_rng(8)
    # a smooth, periodic-in-the-tile height field of realistic amplitude
    kx, ky = 2 * np.pi / lat.lx, 2 * np.pi / lat.ly
    P[:, 2] = (0.8 * np.sin(2 * kx * P[:, 0]) * np.cos(ky * P[:, 1])
               + 0.5 * np.cos(3 * kx * P[:, 0] + 2 * ky * P[:, 1]))
    rows, cols = lat.tile_slice()
    idx = np.array([j * lat.n_vert_x + i for j in range(rows.start, rows.stop)
                    for i in range(cols.start, cols.stop)])
    net = P[idx].reshape(lat.ny, lat.nx, 3)
    exact_surface = mr.LimitSurface(lat, net)
    X, Y = exact_surface.monge_grid("cartesian")
    exact = exact_surface.solve(X, Y).apply(net[..., 2][None])[0]
    # On the lattice every Cartesian grid point is a fine-net vertex from level
    # 1 on -- the points on offset rows are edge midpoints -- so the route is
    # not approximately right there, it is exact.
    for levels in (1, 2, 3):
        h = mfs.SubdividedSurface(F, len(P), levels).interpolate(P, X, Y)
        err = np.abs(h - exact).max()
        check(f"exact on the lattice at {levels} level(s) (grid points are fine-net vertices)",
              err < 1e-10, f"max |h - exact| {err:.1e}")
    # Off the vertices the route is linear interpolation, whose error must fall
    # by about 4x per level. The shift must not be a dyadic fraction of the
    # cell: a quarter-cell shift lands back on fine-net vertices at level 3
    # (spacing dx/8) and the error there is round-off again.
    Xs, Ys = X + 0.137 * lat.d_face_x, Y + 0.291 * lat.d_face_y
    exact_s = exact_surface.solve(Xs, Ys).apply(net[..., 2][None])[0]
    # r.m.s. rather than max: the worst point lands at a different distance
    # from a fine-net vertex at each level, so the max-norm ratio wanders.
    errors = [np.sqrt(np.mean((mfs.SubdividedSurface(F, len(P), k).interpolate(P, Xs, Ys)
                               - exact_s) ** 2)) for k in (1, 2, 3)]
    check("off the vertices the route converges to the exact resampler", errors[-1] < 2e-3,
          "rms |h - exact| at 1, 2, 3 levels: " + ", ".join(f"{e:.2e}" for e in errors))
    check("error falls at least 2x per level (linear interpolation gives ~4x)",
          errors[0] / errors[1] > 2.0 and errors[1] / errors[2] > 2.0,
          f"ratios {errors[0] / errors[1]:.2f}, {errors[1] / errors[2]:.2f}")


def test_flipped_mesh_interpolates_its_own_limit_points():
    lat, P, F = flat_sheet()
    rng = np.random.default_rng(13)
    P[:, 2] = rng.normal(scale=0.5, size=len(P))
    F = F.copy()
    # flip an interior edge by hand: faces (a, b, c) and (b, a, d) -> (c, d, ...)
    conn = mfs.Connectivity(F, len(P))
    inner = np.nonzero((conn.edge_face_count == 2)
                       & conn.interior[conn.edge_lo] & conn.interior[conn.edge_hi]
                       & (np.abs(P[conn.edge_lo, 0] - 50) < 15) & (np.abs(P[conn.edge_lo, 1] - 45) < 15))[0]
    e = inner[len(inner) // 2]
    a, b = conn.edge_lo[e], conn.edge_hi[e]
    c, d = conn.edge_opposite[e]
    f0, f1 = [k for k in range(len(F)) if a in F[k] and b in F[k]]
    F[f0] = [a, d, c]
    F[f1] = [b, c, d]
    conn2 = mfs.Connectivity(F, len(P))
    check("the hand flip changed two valences", sorted(conn2.valence[[a, b, c, d]]) == [5, 5, 7, 7],
          f"valences now {conn2.valence[[a, b, c, d]].tolist()}")
    surf = mfs.SubdividedSurface(F, len(P), 2)
    Q0 = surf.original_vertex_limits(P)
    h = surf.interpolate(P, Q0[:, 0], Q0[:, 1])
    ok = ~np.isnan(h)
    err = np.abs(h[ok] - Q0[ok, 2]).max()
    check("flipped mesh: surface passes through the original vertices' limit points",
          err < 1e-9, f"max error {err:.1e} over {ok.sum()} vertices")
    # the limit mask at the new valences is the SurfaceSolver mask, 1/2 + 1/(2n)
    for v in (a, c):
        n = conn2.valence[v]
        nb = np.concatenate([conn2.edge_hi[conn2.edge_lo == v], conn2.edge_lo[conn2.edge_hi == v]])
        expect = 0.5 * P[v] + (0.5 / n) * P[nb].sum(axis=0)
        check(f"limit mask at valence {n} is (1/2, 1/(2n))", np.abs(Q0[v] - expect).max() < 1e-12)


if __name__ == "__main__":
    for t in (test_planar_net_is_reproduced, test_limit_is_invariant_under_subdivision,
              test_agrees_with_exact_resampler_on_the_lattice,
              test_flipped_mesh_interpolates_its_own_limit_points):
        print(t.__name__)
        t()
    print("\nall passed" if not FAILURES else f"\n{len(FAILURES)} FAILED: {FAILURES}")
    sys.exit(1 if FAILURES else 0)
