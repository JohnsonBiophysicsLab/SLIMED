"""Self-tests for membrane_resample, on data whose answer is known exactly.

Run with ``python analysis/test_membrane_resample.py``, or under pytest.

The module inverts a quartic map and then indexes a 12-node stencil through two
different coordinate systems and a periodic wrap, which is a lot of places for
an off-by-one to hide and produce a plausible-looking surface.  So the checks
here are all of the kind that a wrong stencil cannot pass by luck:

* the node layout is *solved for* rather than compared against, so a permuted
  ordering does not land on the lattice at all;
* the patch is asked to reproduce affine functions, which it does exactly only
  if every basis function is paired with the right control point;
* two faces sharing an edge are evaluated along it, which pins the DOWN stencil
  independently of the UP one;
* resampling at the vertices is compared with the Loop limit mask computed a
  completely different way, by `membrane_spectrum.limit_surface_from_control`.
"""

import math
import sys

import numpy as np

import membrane_spectrum as ms
import membrane_resample as mr


FAILURES = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}{'  -- ' + detail if detail else ''}")
    if not ok:
        FAILURES.append(name)
        # The standalone runner wants the whole list; pytest wants a failure.
        if "pytest" in sys.modules:
            raise AssertionError(f"{name}{'  -- ' + detail if detail else ''}")


def lattice():
    dx = 5.0
    return ms.Lattice(20, 24, dx, math.sqrt(3) / 2 * dx)


def ideal_net(lat, rng=None, sigma=0.0, z=None):
    """A control net on the ideal lattice, optionally displaced in plane."""
    i, j = np.meshgrid(np.arange(lat.nx), np.arange(lat.ny), indexing="xy")
    x = i * lat.d_face_x + np.where(lat.row_is_offset()[:, None], lat.d_face_x / 2, 0.0)
    y = j * lat.d_face_y
    if sigma:
        x = x + rng.normal(scale=sigma, size=x.shape)
        y = y + rng.normal(scale=sigma, size=y.shape)
    if z is None:
        z = np.zeros_like(x) if rng is None else rng.standard_normal(x.shape)
    return np.stack([x, y, z], axis=-1)


# --------------------------------------------------------------------------
def test_shape_functions():
    print("shape functions")
    rng = np.random.default_rng(0)
    v = rng.uniform(0, 1, 3000)
    w = rng.uniform(0, 1, 3000)
    k = (v + w) < 1
    v, w = v[k], w[k]
    N, Nv, Nw = mr.shape_functions(v, w)

    check("twelve of them", N.shape[-1] == 12)
    check("partition of unity", np.abs(N.sum(-1) - 1).max() < 1e-13,
          f"max err {np.abs(N.sum(-1) - 1).max():.1e}")

    eps = 1e-6
    fd_v = (mr.shape_functions(v + eps, w)[0] - mr.shape_functions(v - eps, w)[0]) / (2 * eps)
    fd_w = (mr.shape_functions(v, w + eps)[0] - mr.shape_functions(v, w - eps)[0]) / (2 * eps)
    check("d/dv matches finite differences", np.abs(fd_v - Nv).max() < 1e-8,
          f"max err {np.abs(fd_v - Nv).max():.1e}")
    check("d/dw matches finite differences", np.abs(fd_w - Nw).max() < 1e-8,
          f"max err {np.abs(fd_w - Nw).max():.1e}")


def test_corner_masks():
    print("corner masks are the Loop limit rule")
    for name, (vv, ww), corner in (("u=1", (0., 0.), mr.CORNER_U),
                                   ("v=1", (1., 0.), mr.CORNER_V),
                                   ("w=1", (0., 1.), mr.CORNER_W)):
        N = mr.shape_functions(np.array(vv), np.array(ww))[0]
        half = [k for k in range(12) if abs(N[k] - 0.5) < 1e-12]
        ring = [k for k in range(12) if abs(N[k] - 1 / 12) < 1e-12]
        check(f"{name}: one node at 1/2, six at 1/12",
              half == [corner] and len(ring) == 6,
              f"half {half}, ring {len(ring)}")


def test_patch_layout_is_recovered():
    print("patch layout, solved for rather than assumed")
    a1, a2 = np.array([1.0, 0.0]), np.array([0.5, math.sqrt(3) / 2])
    rng = np.random.default_rng(1)
    v = rng.uniform(0, 1, 4000)
    w = rng.uniform(0, 1, 4000)
    k = (v + w) < 1
    v, w = v[k], w[k]
    N = mr.shape_functions(v, w)[0]
    target = v[:, None] * a1 + w[:, None] * a2          # corner u=1 at the origin
    r, *_ = np.linalg.lstsq(N, target, rcond=None)

    check("linear precision has an exact solution",
          np.abs(N @ r - target).max() < 1e-12,
          f"residual {np.abs(N @ r - target).max():.1e}")
    coef = np.linalg.solve(np.column_stack([a1, a2]), r.T).T
    check("nodes land on the integer lattice", np.abs(coef - np.round(coef)).max() < 1e-9,
          f"max deviation {np.abs(coef - np.round(coef)).max():.1e}")
    check("recovered layout is PATCH_SHEARED",
          np.array_equal(np.round(coef).astype(int), mr.PATCH_SHEARED))


def test_faces_agree_along_shared_edges():
    print("adjacent patches agree where they meet")
    lat = lattice()
    rng = np.random.default_rng(2)
    surf = mr.LimitSurface(lat, ideal_net(lat, rng))
    t = np.linspace(0.0, 1.0, 21)
    one = np.ones_like(t, dtype=np.int64)

    # UP(m,n) and DOWN(m,n) share the u = 0 edge, traversed in opposite senses.
    m, n = 5 * one, 6 * one
    up = surf.evaluate(m, n, mr.UP * one, t, 1.0 - t)
    down = surf.evaluate(m, n, mr.DOWN * one, 1.0 - t, t)
    check("UP and DOWN agree along the shared edge", np.abs(up - down).max() < 1e-12,
          f"max gap {np.abs(up - down).max():.1e} nm")

    # UP(m,n) and DOWN(m,n-1) share UP's w = 0 edge, which runs between the
    # DOWN patch's u-corner and its v-corner -- so it is DOWN's w = 0 edge too,
    # traversed the other way.
    up = surf.evaluate(m, n, mr.UP * one, t, 0.0 * t)
    down = surf.evaluate(m, n - 1, mr.DOWN * one, 1.0 - t, 0.0 * t)
    check("UP and its lower neighbour agree", np.abs(up - down).max() < 1e-12,
          f"max gap {np.abs(up - down).max():.1e} nm")


def test_resampling_at_vertices_is_the_limit_mask():
    print("resampling at the vertices reproduces the Loop limit points")
    lat = lattice()
    rng = np.random.default_rng(3)
    net = ideal_net(lat, rng)
    surf = mr.LimitSurface(lat, net)
    R = surf.solve(*surf.monge_grid("zigzag"))

    check("every grid point was placed", bool(R.converged.all()))
    check("the inversion is exact", R.residual < 1e-11, f"residual {R.residual:.1e} nm")

    z = rng.standard_normal((7, lat.ny, lat.nx))
    got = R.apply(z)
    want = ms.limit_surface_from_control(z, lat)
    check("matches limit_surface_from_control everywhere, seam included",
          np.abs(got - want).max() < 1e-12,
          f"max diff {np.abs(got - want).max():.1e} nm, field rms {want.std():.3f} nm")


def test_linear_precision_on_the_grid():
    print("an affine surface resamples to itself")
    lat = lattice()
    rng = np.random.default_rng(4)
    A, B, C = 0.037, -0.021, 1.3
    for sigma in (0.0, 0.25, 1.0):
        net = ideal_net(lat, rng, sigma=sigma)
        net[..., 2] = A * net[..., 0] + B * net[..., 1] + C
        surf = mr.LimitSurface(lat, net)
        X, Y = surf.monge_grid("cartesian")
        R = surf.solve(X, Y)
        # A tilted plane is not periodic, so it is only a valid reference on the
        # faces whose stencil does not wrap.
        keep = ~surf.crosses_seam(*R.faces)
        err = np.abs(R.apply(net[..., 2]).ravel()[keep]
                     - (A * X + B * Y + C).ravel()[keep]).max()
        check(f"sigma = {sigma:.2f} nm: xy inverted exactly", R.residual < 1e-10,
              f"residual {R.residual:.1e} nm")
        check(f"sigma = {sigma:.2f} nm: z is the same plane", err < 1e-10,
              f"max err {err:.1e} nm over {int(keep.sum())} points")


def test_periodic_translation_invariance():
    print("the wrap is a real periodic image")
    lat = lattice()
    rng = np.random.default_rng(5)
    net = ideal_net(lat, rng)
    surf = mr.LimitSurface(lat, net)
    X, Y = surf.monge_grid("cartesian")
    z = rng.standard_normal((3, lat.ny, lat.nx))

    base = surf.solve(X, Y).apply(z)
    for dx, dy, tag in ((lat.lx, 0.0, "x"), (0.0, lat.ly, "y"), (-lat.lx, lat.ly, "both")):
        got = surf.solve(X + dx, Y + dy).apply(z)
        check(f"shifting the grid by one box in {tag} changes nothing",
              np.abs(got - base).max() < 1e-11,
              f"max diff {np.abs(got - base).max():.1e} nm")


def test_folds_are_detected_not_hidden():
    print("folded meshes are reported")
    lat = lattice()
    rng = np.random.default_rng(6)
    surf = mr.LimitSurface(lat, ideal_net(lat, rng))
    check("an undeformed net has no folds", surf.inverted_control_triangles() == 0)

    folded = mr.LimitSurface(lat, ideal_net(lat, rng, sigma=3.0))
    n_fold = folded.inverted_control_triangles()
    check("a badly deformed net does have folds", n_fold > 0, f"{n_fold} folds")

    import warnings
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        R = folded.solve(*folded.monge_grid("cartesian"))
    check("solve warns rather than returning silence", len(caught) > 0)
    check("and flags the points it could not place", not bool(R.converged.all()),
          f"{int((~R.converged).sum())} unplaced of {R.converged.size}")


def test_cartesian_spectrum():
    print("the cartesian transform and its wavevectors")
    lat = lattice()
    qx, qy, q = mr.cartesian_q_grid(lat)
    check("bins are 2*pi*n/L",
          abs(qx[0, 1] - 2 * math.pi / lat.lx) < 1e-12
          and abs(qy[1, 0] - 2 * math.pi / lat.ly) < 1e-12)
    q_raw = mr.cartesian_q_grid(lat, fold_to_brillouin_zone=False)[2]
    check("folding never lengthens a wavevector", (q <= q_raw + 1e-12).all())

    # A single cosine of known amplitude has a spectrum we can write down:
    # h = a cos(q.r) puts |h_q|^2 = a^2/4 into each of the two modes +-q.
    i, j = np.meshgrid(np.arange(lat.nx), np.arange(lat.ny), indexing="xy")
    a, nx_mode, ny_mode = 0.37, 2, 3
    h = a * np.cos(2 * math.pi * (nx_mode * i / lat.nx + ny_mode * j / lat.ny))
    S = mr.cartesian_spectrum(h[None], lat)
    check("amplitude lands where it should",
          abs(S[ny_mode, nx_mode] - a ** 2 / 4) < 1e-12,
          f"got {S[ny_mode, nx_mode]:.6f}, want {a ** 2 / 4:.6f}")
    check("and nowhere else", abs(S.sum() - 2 * a ** 2 / 4) < 1e-12)


def test_read_tile_control():
    print("reading all three columns of a trajectory")
    import tempfile, os
    lat = lattice()
    rng = np.random.default_rng(8)
    frames = rng.standard_normal((4, lat.n_vertices, 3))
    fd, path = tempfile.mkstemp(suffix=".csv")
    with os.fdopen(fd, "w") as fh:
        for f in frames:
            fh.write(",".join(f"{v:.12g}" for v in f.ravel()) + ",\n")
        fh.write("1,2,3,\n")           # a half-written trailing line
    try:
        got = mr.read_tile_control(path, lat)
        rows, cols = lat.tile_slice()
        want = frames.reshape(4, lat.n_vert_y, lat.n_vert_x, 3)[:, rows, cols]
        check("shape is (frames, ny, nx, 3)", got.shape == (4, lat.ny, lat.nx, 3),
              str(got.shape))
        check("values are the tile's own", np.abs(got - want).max() < 1e-9)
        check("a truncated last line is skipped, not reshaped", len(got) == 4)
        check("reference_net is the first frame",
              np.abs(mr.reference_net(path, lat) - want[0]).max() < 1e-9)
    finally:
        os.unlink(path)


TESTS = (test_shape_functions, test_corner_masks, test_patch_layout_is_recovered,
         test_faces_agree_along_shared_edges,
         test_resampling_at_vertices_is_the_limit_mask,
         test_linear_precision_on_the_grid, test_periodic_translation_invariance,
         test_folds_are_detected_not_hidden, test_cartesian_spectrum,
         test_read_tile_control)


if __name__ == "__main__":
    for fn in TESTS:
        fn()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILED: {', '.join(FAILURES)}")
        sys.exit(1)
    print("all checks passed")
