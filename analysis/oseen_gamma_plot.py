#!/usr/bin/env python
"""
Relaxation rate Gamma(q) of a flat sheet's height modes against the
free-draining law (q^4) and the Oseen law (q^3).

For each run: the limit surface on the tile's cartesian Monge grid (Loop
subdivision of the control net, the grid scaled with the box frame by frame,
as analysis/flat_tension_matrix.ipynb does), the height spectrum and its
(kc, sigma) fit over |q| d <= pi/2, and the relaxation rate of every mode from
the integral of its autocorrelation to 1/e (membrane_spectrum.relaxation_rate).
A mode whose correlation time is under 5 % of the analysed record is
"resolved" (filled symbols); between 5 % and --marginal of the record it is
shown open, as a mode the record only just times.

The lines, with D the run's diffConst, a the area per vertex and sigma the
spectrum's tension (floored at zero):

    free-draining  Gamma_FD = (D a / kT)(kc q^4 + sigma q^2)
    Oseen          Gamma_Os = (kc q^3 + sigma q) / (4 eta)

which cross at q* = kT / (4 eta D a) whatever kc and sigma are
(analysis/oseen_matched_diffusion.py). A --reference run made at another D
is drawn rescaled by D_run / D_ref: in overdamped dynamics D only sets the
clock, so the two must collapse.

    python analysis/oseen_gamma_plot.py --run DIR [--reference DIR] [--eta 1.0] --out PREFIX

writes PREFIX.png and PREFIX.svg and prints the per-run numbers. The
resampled heights are cached as <run's parent>/analysis/heights_<run>.npy,
the notebook's convention, so a run the notebook has read is not resampled again.
"""
import argparse
import os
import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fluidity as fl  # noqa: E402
import membrane_fluid_surface as mfs  # noqa: E402
import membrane_resample as mr  # noqa: E402
import membrane_spectrum as ms  # noqa: E402
from oseen_matched_diffusion import MPAS_TO_PN_US_PER_NM2, read_params  # noqa: E402

QDX_MAX = np.pi / 2
NBLOCK = 8
TEAL, BRICK, GREY, INK = "#0f766e", "#b1442a", "#8d9794", "#111827"


def unique_modes(lat):
    keep = np.zeros((lat.ny, lat.nx), bool)
    seen = set()
    for iy in range(lat.ny):
        for ix in range(lat.nx):
            if (iy, ix) == (0, 0) or ((-iy) % lat.ny, (-ix) % lat.nx) in seen:
                continue
            seen.add((iy, ix))
            keep[iy, ix] = True
    return keep


class Run:
    def __init__(self, directory, burn_fraction, levels):
        self.dir = os.path.abspath(os.path.expanduser(directory))
        self.name = os.path.basename(self.dir)
        p = self.params = read_params(os.path.join(self.dir, "input.params"))
        self.kbt = float(p.get("KBT", 4.17))
        self.kc = float(p["kcMembraneBending"])
        self.D = float(p["diffConst"])
        self.dt = float(p["timeStep"])
        self.dt_frame = self.dt * int(p["meshpointOutputInterval"])
        self.lat = ms.lattice_from_params(os.path.join(self.dir, "input.params"))
        self.a = self.lat.area / (self.lat.nx * self.lat.ny)
        self.levels = levels

        t0 = time.perf_counter()
        self.run = fl.load_run(self.dir, "input", int(p["meshpointOutputInterval"]))
        self.n = len(self.run.frame_iterations)
        self.burn = int(burn_fraction * self.n)
        self.scale = np.ones(self.n)
        box = os.path.join(self.dir, "inputbox.csv")
        if os.path.exists(box):
            b = np.genfromtxt(box, delimiter=",", names=True)
            self.scale = np.interp(self.run.frame_iterations, b["iteration"], b["scaleX"])
        print(f"{self.name}: {self.n} frames, D = {self.D:g} nm^2/us, dt = {self.dt:g} us, "
              f"frame every {self.dt_frame:.4g} us, read in {time.perf_counter() - t0:.0f} s")

    def heights(self):
        cache_dir = os.path.join(os.path.dirname(self.dir), "analysis")
        cache = os.path.join(cache_dir, f"heights_{self.name}.npy")
        frames = np.arange(self.burn, self.n)
        if os.path.exists(cache) and np.load(cache, mmap_mode="r").shape[0] == len(frames):
            return np.load(cache)
        t0 = time.perf_counter()
        net0 = mr.reference_net(os.path.join(self.dir, "meshpointinput.csv"), self.lat)
        X0, Y0 = mr.LimitSurface(self.lat, net0).monge_grid("cartesian")
        surfaces = {}
        h = np.empty((len(frames),) + X0.shape)
        for j, k in enumerate(frames):
            usable = [i for i in self.run.face_frames if i <= self.run.frame_iterations[k]]
            key = max(usable) if usable else min(self.run.face_frames)
            surf = surfaces.get(key)
            if surf is None:
                surf = surfaces[key] = mfs.SubdividedSurface(self.run.face_frames[key], self.run.coords.shape[1], self.levels)
            h[j] = surf.interpolate(self.run.coords[k], X0 * self.scale[k], Y0 * self.scale[k])
        if not np.isfinite(h).all():
            raise RuntimeError(f"{self.name}: a grid point fell outside the mesh")
        h -= h.mean(axis=(1, 2), keepdims=True)
        os.makedirs(cache_dir, exist_ok=True)
        np.save(cache, h)
        print(f"{self.name}: {len(frames)} frames resampled in {time.perf_counter() - t0:.0f} s, cached at {cache}")
        return h

    def analyse(self, eta, marginal):
        h = self.heights()
        lat = self.lat
        H = np.fft.fft2(h, axes=(-2, -1))
        keep = unique_modes(lat)
        _, _, q2d = mr.cartesian_q_grid(lat)
        q2d = q2d / self.scale[self.burn:].mean()
        area = lat.area * self.scale[self.burn:].mean() ** 2
        win2d = keep & (q2d * lat.d_face_x * self.scale[self.burn:].mean() <= QDX_MAX)

        # spectrum and (kc, sigma), one fit per time block
        edges = np.linspace(0, len(H), NBLOCK + 1).astype(int)
        fits = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            S = np.mean(np.abs(H[lo:hi]) ** 2, axis=0) / (lat.nx * lat.ny) ** 2
            fits.append(ms.fit_kc_sigma(q2d[win2d], S[win2d], area, self.kbt))
        fits = np.array(fits)
        self.kc_fit, self.sigma_fit = fits.mean(axis=0)
        self.kc_err, self.sigma_err = fits.std(axis=0, ddof=1) / np.sqrt(NBLOCK)
        self.S = np.mean(np.abs(H) ** 2, axis=0)[keep] / (lat.nx * lat.ny) ** 2

        # rates: resolved under the 5 % guard; marginal up to `marginal` of the record
        acf = ms.mode_autocorrelation(H)
        q, G, resolved = [], [], []
        for iy, ix in zip(*np.nonzero(keep)):
            g = ms.relaxation_rate(acf[:, iy, ix], self.dt_frame, n_frames=len(H), max_tau_fraction=marginal)
            if not np.isfinite(g):
                continue
            tau = 1.0 / g
            q.append(q2d[iy, ix]); G.append(g); resolved.append(tau <= 0.05 * len(H) * self.dt_frame)
        o = np.argsort(q)
        self.q, self.G, self.resolved = np.array(q)[o], np.array(G)[o], np.array(resolved)[o]
        self.q_all = np.sort(q2d[keep])
        self.S_q = q2d[keep]
        self.window = self.q * lat.d_face_x * self.scale[self.burn:].mean() <= QDX_MAX
        sig = max(self.sigma_fit, 0.0)
        self.G_fd = self.D * self.a / self.kbt * (self.kc * self.q ** 4 + sig * self.q ** 2)
        self.G_os = (self.kc * self.q ** 3 + sig * self.q) / (4.0 * eta)
        self.q_star = self.kbt / (4.0 * eta * self.D * self.a)
        sel = self.resolved & self.window
        self.ratio_fd = (np.mean(self.G[sel] / self.G_fd[sel]), np.std(self.G[sel] / self.G_fd[sel], ddof=1) / np.sqrt(sel.sum()))
        self.slope = np.polyfit(np.log(self.q[sel]), np.log(self.G[sel]), 1)[0]
        self.record_us = len(H) * self.dt_frame
        return self


def fd_line(run, q, eta):
    return run.D * run.a / run.kbt * (run.kc * q ** 4 + max(run.sigma_fit, 0.0) * q ** 2)


def os_line(run, q, eta):
    return (run.kc * q ** 3 + max(run.sigma_fit, 0.0) * q) / (4.0 * eta)


def plain_log_ticks(ax, ticks=(0.05, 0.07, 0.1, 0.15, 0.2, 0.3)):
    ax.xaxis.set_major_locator(FixedLocator(list(ticks)))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))


def plot(run, ref, eta, eta_mpas, out):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.0))
    qmin = run.q_all[0]
    qq = np.geomspace(qmin * 0.9, QDX_MAX / run.lat.d_face_x * 1.05, 200)
    c = run.D / ref.D if ref is not None else None

    # (a) Gamma(q)
    ax = axes[0]
    if ref is not None:
        w = ref.window
        ax.plot(ref.q[w & ref.resolved], ref.G[w & ref.resolved] * c, "o", ms=4, color=GREY, alpha=0.8,
                label=f"D = {ref.D:g} run, Γ × {c:.0f}")
    w = run.window
    ax.plot(run.q[w & run.resolved], run.G[w & run.resolved], "o", ms=6, color=TEAL, label=f"D = {run.D:.4g} run, resolved")
    if (w & ~run.resolved).any():
        ax.plot(run.q[w & ~run.resolved], run.G[w & ~run.resolved], "o", ms=6, mfc="white", color=TEAL,
                label="marginal (τ 5–15 % of record)")
    ax.plot(qq, fd_line(run, qq, eta), "-", color=TEAL, lw=1.6, label=r"free-draining  $(Da/k_BT)(\kappa q^4+\sigma q^2)$")
    ax.plot(qq, os_line(run, qq, eta), "--", color=BRICK, lw=1.8, label=rf"Oseen  $(\kappa q^3+\sigma q)/4\eta$, η = {eta_mpas:g} mPa·s")
    star_in = qq[0] < run.q_star < qq[-1]
    if star_in:
        ax.axvline(run.q_star, color=INK, lw=0.8, ls=":")
        ax.text(run.q_star * 1.03, 0.97, f"q* = {run.q_star:.4f}", fontsize=8, color=INK, va="top",
                transform=ax.get_xaxis_transform())
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("|q| (1/nm)"); ax.set_ylabel("Γ(q) (1/µs)")
    ax.set_title("relaxation rate of the height modes", fontsize=10)
    ax.legend(fontsize=7.5, loc="upper left")
    plain_log_ticks(ax)

    # (b) Gamma / Gamma_Oseen: q / q* for any free-draining run
    ax = axes[1]
    if ref is not None:
        w = ref.window & ref.resolved
        ax.plot(ref.q[w], ref.G[w] * c / os_line(run, ref.q[w], eta), "o", ms=4, color=GREY, alpha=0.8, label=f"D = {ref.D:g} run, × {c:.0f}")
    w = run.window
    ax.plot(run.q[w & run.resolved], run.G[w & run.resolved] / run.G_os[w & run.resolved], "o", ms=6, color=TEAL, label=f"D = {run.D:.4g} run")
    if (w & ~run.resolved).any():
        ax.plot(run.q[w & ~run.resolved], run.G[w & ~run.resolved] / run.G_os[w & ~run.resolved], "o", ms=6, mfc="white", color=TEAL)
    ax.plot(qq, qq / run.q_star, "-", color=TEAL, lw=1.4, label=r"$4\eta D a q / k_BT = q/q^*$")
    ax.axhline(1.0, color=BRICK, lw=1.2, ls="--", label="Oseen")
    if star_in:
        ax.axvline(run.q_star, color=INK, lw=0.8, ls=":")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("|q| (1/nm)"); ax.set_ylabel("Γ / Γ_Oseen")
    ax.set_title("the simulation against hydrodynamics", fontsize=10)
    ax.legend(fontsize=7.5, loc="upper left")
    plain_log_ticks(ax)

    # (c) Gamma / Gamma_FD: the simulation's own law
    ax = axes[2]
    if ref is not None:
        w = ref.window & ref.resolved
        ax.plot(ref.q[w], ref.G[w] / ref.G_fd[w], "o", ms=4, color=GREY, alpha=0.8, label=f"D = {ref.D:g} run")
    w = run.window
    ax.plot(run.q[w & run.resolved], run.G[w & run.resolved] / run.G_fd[w & run.resolved], "o", ms=6, color=TEAL, label=f"D = {run.D:.4g} run")
    if (w & ~run.resolved).any():
        ax.plot(run.q[w & ~run.resolved], run.G[w & ~run.resolved] / run.G_fd[w & ~run.resolved], "o", ms=6, mfc="white", color=TEAL)
    ax.axhline(1.0, color=TEAL, lw=1.2)
    ax.set_xscale("log"); ax.set_ylim(0, 2)
    ax.set_xlabel("|q| (1/nm)"); ax.set_ylabel("Γ / Γ_free-draining")
    ax.set_title(f"the q⁴ law itself: mean {run.ratio_fd[0]:.2f} ± {run.ratio_fd[1]:.2f}, slope {run.slope:.2f}", fontsize=10)
    ax.legend(fontsize=7.5, loc="upper left")
    plain_log_ticks(ax)

    for ax in axes:
        ax.set_xlim(qq[0], qq[-1])
    fig.suptitle(f"{run.name}: D = {run.D:.4g} nm²/µs; free-draining meets Oseen at q* = kT/(4ηDa) = {run.q_star:.4g} /nm "
                 f"(κ = {run.kc:g} pN·nm, a = {run.a:.2f} nm², spectrum σ = {run.sigma_fit:+.2f} pN/nm)", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for ext in ("png", "svg"):
        fig.savefig(f"{out}.{ext}", dpi=170)
    plt.close(fig)
    print(f"wrote {out}.png and {out}.svg")


def report(r, eta):
    sel = r.resolved & r.window
    lo = np.argmin(r.q[r.window]) if r.window.any() else None
    print(f"\n{r.name}: D = {r.D:g} nm^2/us, record {r.record_us:.4g} us, q* = kT/(4 eta D a) = {r.q_star:.4g} /nm")
    print(f"  spectrum: kc = {r.kc_fit:.1f} +- {r.kc_err:.1f} pN nm (input {r.kc:g}), sigma = {r.sigma_fit:+.3f} +- {r.sigma_err:.3f} pN/nm")
    print(f"  {sel.sum()} resolved and {(r.window & ~r.resolved).sum()} marginal modes in the window, from |q| = {r.q[r.window].min():.4f}")
    print(f"  Gamma / Gamma_FD = {r.ratio_fd[0]:.3f} +- {r.ratio_fd[1]:.3f}, d ln Gamma / d ln q = {r.slope:.2f}")
    if lo is not None:
        qlo = r.q[r.window][lo]
        print(f"  lowest timed mode |q| = {qlo:.4f}: Gamma = {r.G[r.window][lo]:.4g} /us, "
              f"Oseen {((r.kc * qlo ** 3 + max(r.sigma_fit, 0) * qlo) / (4 * eta)):.4g}, free-draining {r.G_fd[r.window][lo]:.4g}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--run", required=True, help="the run to test (normally the matched one)")
    parser.add_argument("--reference", default=None, help="a run at another D, drawn rescaled by D_run / D_ref")
    parser.add_argument("--eta", type=float, default=1.0, help="solvent viscosity, mPa s")
    parser.add_argument("--burn", type=float, default=0.2, help="fraction of each run discarded")
    parser.add_argument("--levels", type=int, default=2, help="subdivision levels of the resampler")
    parser.add_argument("--marginal", type=float, default=0.15, help="largest tau / record still drawn (open symbols)")
    parser.add_argument("--out", required=True, help="output prefix; .png and .svg are written")
    args = parser.parse_args()
    eta = args.eta * MPAS_TO_PN_US_PER_NM2
    run = Run(args.run, args.burn, args.levels).analyse(eta, args.marginal)
    ref = Run(args.reference, args.burn, args.levels).analyse(eta, args.marginal) if args.reference else None
    for r in (run, ref):
        if r is not None:
            report(r, eta)
    if ref is not None:
        sel = ref.resolved & ref.window
        print(f"\nD_run / D_ref = {run.D / ref.D:.4g}; frames per relaxation time are the same when D dt is: "
              f"{run.D * run.dt:.4g} vs {ref.D * ref.dt:.4g} nm^2")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    plot(run, ref, eta, args.eta, args.out)


if __name__ == "__main__":
    main()
