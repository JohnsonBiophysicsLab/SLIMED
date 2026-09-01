"""Self-tests for membrane_spectrum, on data whose answer is known exactly.

Run with ``python analysis/test_membrane_spectrum.py``.  No pytest required,
and numpy is the only dependency, so this runs anywhere the simulation does.

The point of the synthetic tests is that they separate two questions that are
easy to confuse when staring at a simulation that does not give q^-4: is the
*analysis* wrong, or is the *simulation* wrong?  Here the field is drawn from a
known spectrum, so any deviation is the analysis's fault.
"""

import math
import sys

import numpy as np

import membrane_spectrum as ms


FAILURES = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}{'  -- ' + detail if detail else ''}")
    if not ok:
        FAILURES.append(name)


def lattice():
    dx = 5.0
    return ms.Lattice(28, 32, dx, math.sqrt(3) / 2 * dx)


# --------------------------------------------------------------------------
def test_geometry():
    print("geometry")
    lat = lattice()
    check("tile is nFace - 6 on each axis", (lat.nx, lat.ny) == (22, 26))
    check("row count is even so the zig-zag closes", lat.ny % 2 == 0)
    check("area is Lx*Ly", abs(lat.area - lat.lx * lat.ly) < 1e-9)

    # The parity pattern has to match Mesh::set_vertices_faces_flat(), which
    # offsets even j.  Tile row s is mesh row 3 + s.
    off = lat.row_is_offset()
    check("even mesh rows are the offset ones", (not off[0]) and off[1])

    # And lattice_from_params has to reproduce Mesh::set_axes_division_flat().
    import tempfile, os
    fd, path = tempfile.mkstemp(suffix=".params", text=True)
    with os.fdopen(fd, "w") as fh:
        fh.write("sideX = 100.0\nsideY = 100.0\nlFace = 5.0\n")
    got = ms.lattice_from_params(path)
    os.unlink(path)
    check("params -> nFaceX 20, nFaceY 24 (rounded up to even)",
          (got.n_face_x, got.n_face_y) == (20, 24),
          f"got {got.n_face_x}, {got.n_face_y}")


# --------------------------------------------------------------------------
def test_transform_is_the_true_dft():
    """tile_fft2 must equal the explicit sum over the real vertex positions."""
    print("transform")
    lat = lattice()
    rng = np.random.default_rng(0)
    z = rng.standard_normal((lat.ny, lat.nx))

    dx, dy = lat.d_face_x, lat.d_face_y
    p = np.arange(lat.nx)[None, :]
    s = np.arange(lat.ny)[:, None]
    x = p * dx + np.where(lat.row_is_offset()[:, None], dx / 2, 0.0)
    y = s * dy

    qx, qy, _ = ms.q_grid(lat, fold_to_brillouin_zone=False)
    brute = np.empty((lat.ny, lat.nx), complex)
    for a in range(lat.ny):
        for b in range(lat.nx):
            brute[a, b] = np.sum(z * np.exp(-1j * (qx[a, b] * x + qy[a, b] * y)))
    fast = ms.tile_fft2(z, lat)
    err = np.abs(brute - fast).max() / np.abs(brute).max()
    check("matches the explicit sum over true (x, y)", err < 1e-12, f"rel err {err:.2e}")

    back = ms.tile_ifft2(fast, lat)
    check("round trips", np.abs(back - z).max() < 1e-12)

    # A plane wave at one lattice wavevector must give exactly one peak.
    iy, ix = 3, 5
    wave = np.cos(qx[iy, ix] * x + qy[iy, ix] * y)
    Hw = np.abs(ms.tile_fft2(wave, lat))
    peak = Hw.max()
    Hw[iy, ix] = 0.0
    Hw[(-iy) % lat.ny, (-ix) % lat.nx] = 0.0
    check("a lattice plane wave is a single pair of peaks",
          Hw.max() < 1e-9 * peak, f"leak {Hw.max()/peak:.2e}")

    # The naive index transform must NOT do this -- that is the bug being
    # guarded against, so assert that it really is different.
    Hn = np.abs(ms.naive_index_fft2(wave))
    Hn[iy, ix] = 0.0
    Hn[(-iy) % lat.ny, (-ix) % lat.nx] = 0.0
    check("naive index fft2 smears that plane wave (the bug)",
          Hn.max() > 1e-3 * peak, f"leak {Hn.max()/peak:.2e}")


def test_mask_symbol():
    print("limit mask")
    lat = lattice()
    dx, dy = lat.d_face_x, lat.d_face_y
    qx, qy, _ = ms.q_grid(lat, fold_to_brillouin_zone=False)
    a1, a2 = np.array([dx, 0.0]), np.array([dx / 2, dy])
    a3 = a2 - a1
    analytic = 0.5 + sum(2 * np.cos(qx * v[0] + qy * v[1]) for v in (a1, a2, a3)) / 12.0
    m = ms.limit_mask_symbol(lat)
    check("symbol equals 1/2 + (1/6) sum cos(q.a)",
          np.abs(m - analytic).max() < 1e-12, f"max dev {np.abs(m-analytic).max():.2e}")
    check("m(0) = 1 (mask preserves a constant)", abs(m[0, 0] - 1.0) < 1e-12)
    check("m > 0 everywhere on this lattice", m.min() > 0, f"min {m.min():.3f}")

    # Applied to a plane in x and y the mask must be a no-op: that is what
    # makes the limit surface sit at the same (x, y) as the control net.
    s = np.arange(lat.ny)[:, None]
    p = np.arange(lat.nx)[None, :]
    plane = 3.0 * (p * dx + np.where(lat.row_is_offset()[:, None], dx / 2, 0.0)) \
        - 2.0 * (s * dy) + 7.0
    # (only the interior is meaningful -- the wrap makes the ramp jump at the
    # seam, so compare away from it)
    got = ms.limit_surface_from_control(plane, lat)[2:-2, 2:-2]
    check("mask reproduces a linear ramp exactly",
          np.abs(got - plane[2:-2, 2:-2]).max() < 1e-9)


# --------------------------------------------------------------------------
def test_recovers_a_known_spectrum():
    """The whole point: feed in kT/(A kc q^4), get kc back."""
    print("spectrum recovery")
    lat = lattice()
    kbt, kc_true, sigma_true = 4.17, 83.4, 0.0
    qx, qy, q = ms.q_grid(lat)
    with np.errstate(divide="ignore"):
        target = kbt / (lat.area * (kc_true * q**4 + sigma_true * q**2))
    target[0, 0] = 0.0

    rng = np.random.default_rng(12345)
    z = ms.synthesize_heights(lat, target, 20000, rng)
    check("synthesised field is real", np.isrealobj(z))

    S = ms.power_spectrum(ms.tile_fft2(z, lat), lat)
    good = q > 0
    dev = np.abs(S[good] / target[good] - 1.0)
    check("per-mode spectrum within sampling error",
          np.median(dev) < 0.03, f"median |S/target - 1| = {np.median(dev):.4f}")

    qc, sm, se, cnt = ms.radial_average(q, S, n_bins=20)
    slope = ms.log_slope(qc, sm)
    check("log-log slope is -4", abs(slope + 4) < 0.03, f"slope {slope:.4f}")

    kc, sig = ms.fit_kc_sigma(qc, sm, lat.area, kbt)
    check("kc recovered", abs(kc / kc_true - 1) < 0.03, f"kc = {kc:.2f} vs {kc_true}")
    check("tension recovered as ~0", abs(sig) < 0.02,
          f"sigma = {sig:.4f} pN/nm")

    # With a tension the two-parameter fit must separate them.
    sigma_true = 0.5
    with np.errstate(divide="ignore"):
        target2 = kbt / (lat.area * (kc_true * q**4 + sigma_true * q**2))
    target2[0, 0] = 0.0
    z2 = ms.synthesize_heights(lat, target2, 20000, np.random.default_rng(999))
    qc2, sm2, _, _ = ms.radial_average(
        q, ms.power_spectrum(ms.tile_fft2(z2, lat), lat), n_bins=20)
    kc2, sig2 = ms.fit_kc_sigma(qc2, sm2, lat.area, kbt)
    check("kc and sigma separated when both are present",
          abs(kc2 / kc_true - 1) < 0.06 and abs(sig2 / sigma_true - 1) < 0.15,
          f"kc = {kc2:.1f}, sigma = {sig2:.3f}")


def test_legacy_pipeline_fails_on_the_same_data():
    """Reproduce the original notebook's reduction; it must miss the answer.

    This is the control that shows the corrections above are load bearing:
    identical input data, only the analysis changed.
    """
    print("legacy reduction on ground-truth data")
    lat = lattice()
    kbt, kc_true = 4.17, 83.4
    qx, qy, q = ms.q_grid(lat)
    with np.errstate(divide="ignore"):
        target = kbt / (lat.area * kc_true * q**4)
    target[0, 0] = 0.0
    z = ms.synthesize_heights(lat, target, 8000, np.random.default_rng(7))

    # (a) index-space fft2, (b) "first quarter" fold written as a repeated row
    #     slice, (c) |q| from the raw row/column index over a length that is a
    #     vertex count rather than a length in nm.
    arr = np.abs(ms.naive_index_fft2(z))
    trimmed = arr[:, : arr.shape[1] // 2][:, : arr.shape[2] // 2]
    s_legacy = np.mean(trimmed**2, axis=0).ravel()
    n_row, n_col = trimmed.shape[1], trimmed.shape[2]
    l_mean = (lat.n_vert_x + lat.n_vert_y) / 2
    q_legacy = np.sqrt(np.add.outer((np.arange(n_row) / l_mean) ** 2,
                                    (np.arange(n_col) / l_mean) ** 2)).ravel()
    keep = q_legacy > 0
    slope = ms.log_slope(q_legacy[keep], s_legacy[keep])
    check("legacy slope is nowhere near -4", abs(slope + 4) > 1.0,
          f"legacy slope {slope:.3f}")

    qc, sm, _, _ = ms.radial_average(q, ms.power_spectrum(ms.tile_fft2(z, lat), lat))
    check("corrected slope on the same frames is -4",
          abs(ms.log_slope(qc, sm) + 4) < 0.05,
          f"corrected slope {ms.log_slope(qc, sm):.3f}")


def test_autocorrelation():
    print("autocorrelation")
    rng = np.random.default_rng(3)
    n, lam = 200000, 0.02
    x = np.zeros(n)
    for i in range(1, n):
        x[i] = (1 - lam) * x[i - 1] + rng.standard_normal()
    a = ms.mode_autocorrelation(x.reshape(-1, 1, 1), max_lag=400)[:, 0, 0]
    check("lag-0 is 1", abs(a[0] - 1) < 1e-12)
    check("one-step decay recovers lambda", abs((1 - a[1]) / lam - 1) < 0.05,
          f"1 - C(1) = {1-a[1]:.5f} vs {lam}")
    rate = ms.relaxation_rate(a, dt=1.0)
    check("integral estimator recovers the rate", abs(rate / lam - 1) < 0.08,
          f"rate {rate:.5f} vs {lam}")


if __name__ == "__main__":
    for fn in (test_geometry, test_transform_is_the_true_dft, test_mask_symbol,
               test_recovers_a_known_spectrum,
               test_legacy_pipeline_fails_on_the_same_data,
               test_autocorrelation):
        fn()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILED: {', '.join(FAILURES)}")
        sys.exit(1)
    print("all checks passed")
