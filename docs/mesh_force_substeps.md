# Raising the time step: substepped mesh walls (`meshForceSubsteps`)

**Status:** implemented and benchmarked, JohnsonBiophysicsLab/SLIMED#19 (branch
`claude/sweet-hawking-46ebb2`, on `main @ 341e6ff`). Off by default; with it off the shipped
legacy periodic and Mixed workloads are byte-identical.
**Companion:** `docs/resolution_limits.md` (sections 2.6 and 3), whose tables this document
derives and explains. The negative result on truncated noise is JohnsonBiophysicsLab/SLIMED#18.
**Runs:** `~/slimed_runs/trunc_noise`, `~/slimed_runs/mts`, `~/slimed_runs/substep_res`
(2026-10-20/21, about 1100 runs in all), generated and reduced by
`analysis/substep_resolution_runs.py` and `analysis/substep_resolution_report.py`.

## 1. Summary

Before this work the Brownian step of the flat sheet was limited by stability: the frozen lFace-5
tile diverged above $\Delta t \approx 0.0025\ \mu s$ and the fluid one above $0.002\ \mu s$, and the
spectrum below those steps was already accurate to a few per cent. The divergence turned out to
have nothing to do with the physics being resolved. It is an explicit-Euler overshoot of the
mesh-quality walls -- the edge tether, the altitude floor and the crease wall -- which act on the
control points, where the limit-surface mapping amplifies their mobility by up to a factor of 16.
Integrating only those walls in $n$ substeps of $\Delta t / n$, while bending, the area constraint
and the noise keep the full step, moves the stability limit out to that of bending itself,
$\Delta t \approx 0.03\,(\ell/5)^2\ \mu s$ with $\ell$ the face length in nm, twelve to fifty times
the old bound.

With stability out of the way the step is limited by accuracy, and the accuracy follows a single
law at every mesh size. Bending is still integrated by Euler-Maruyama, which overstates the
variance and the relaxation rate of a mode of rate $\Gamma(q)$ by about $\Gamma(q)\,\Delta t / 2$.
Keeping both the static spectrum $S(q)$ and the rates $\Gamma(q)$ within 5 % over the range the mesh
resolves, $q\ell \le 2$, requires $\Gamma(2/\ell)\,\Delta t \lesssim 0.1$, which is

$$
\Delta t_{\rm acc} \approx 0.01 \left(\frac{\ell}{5\ \mathrm{nm}}\right)^{2} \mu s
\qquad (D = 1\ \mathrm{nm^2/\mu s},\ \kappa_c = 83.4\ \mathrm{pN\,nm}).
$$

Against the old bound this is a factor of 4 at $\ell = 5$, 8 at $\ell = 10$ and 27 at $\ell = 20$
in the step itself. Because each substep costs about a quarter of a full step on the frozen sheet,
the gain in simulated time per wall-clock second is 2.4, 2.6 and 3.2; on the fluid sheet the
substeps are free and the gain at $\ell = 5$ is the step ratio, 2.5. The mesh-size bound is
unchanged: modes with $q\ell \le 2$ are within 5 %, i.e. $\ell \le \lambda_{\min}/\pi$.

| face length $\ell$ | old bound (stable = accurate) | new stable bound | new accurate bound | gain in $\Delta t$ | gain in throughput |
| --- | --- | --- | --- | --- | --- |
| 5 nm, frozen | 0.0025 µs | 0.03 µs ($n = 12$) | 0.01 µs ($n = 4$) | ×4 | ×2.4 |
| 10 nm, frozen | 0.005 µs | 0.12 µs ($n = 30$) | 0.04 µs ($n = 10$) | ×8 | ×2.6 |
| 20 nm, frozen | 0.006 µs | 0.32 µs ($n = 64$) | 0.16 µs ($n = 32$) | ×27 | ×3.2 |
| 5 nm, fluid ($\nu = 0.5$) | 0.002 µs | ≥ 0.02 µs ($n = 8$) | ≈ 0.005 µs ($n = 2$) | ×2.5 | ×2.5 |

So the accurate step is roughly ten times larger on average across these meshes, but the factor
depends on the mesh: it is smallest on the fine lFace-5 mesh, where bending is stiffest, and grows
as $\ell$ because the old wall bound grew only about as $\ell^{0.5}$ while the bending-limited
accurate step grows as $\ell^2$.

## 2. The update and where its stiffness lives

The Brownian dynamics moves the Loop limit surface $S$, an $N \times 3$ array of surface points,
and recovers the control net from it, $C = M^{-1} S$, where $M$ is the limit mask. The energy is a
function of the control points, so the force the mesh reports is $F_C = -\partial E / \partial C$,
and the force conjugate to $S$ is $F_S = M^{-T} F_C$. One step of length $\Delta t$ is the
Euler-Maruyama update

$$
S^{t+\Delta t} = S^{t} + \frac{D\,\Delta t}{k_BT}\, F_S(S^t) + \sqrt{2 D \Delta t}\;\xi^t,
\qquad \xi^t \sim \mathcal N(0, \mathbb 1),
$$

with the mean height of the integrated vertices pinned and the periodic images re-asserted after
the map back. Expressed on the control points this is a mobility

$$
\mu_C = \frac{D}{k_BT}\, M^{-1} M^{-T}.
$$

On the regular valence-6 lattice $M$ is the convolution with the limit mask (one half on the vertex,
one twelfth on each neighbour), whose Fourier symbol is

$$
m(\mathbf k) = \tfrac12 + \tfrac1{12}\sum_{j=1}^{6} \cos(\mathbf k\cdot \mathbf e_j)
\in \left[\tfrac14,\ 1\right],
$$

equal to 1 at $\mathbf k = 0$ and $\tfrac14$ at the corner of the Brillouin zone. The control-point
mobility therefore reaches $16\,D/k_BT$ for the shortest-wavelength pattern of displacements.
Bending does not feel this: it is computed on the limit surface, and the factors $M^{-1}$ and
$M^{-T}$ cancel the $M$ inside its gradient. The mesh-quality walls do feel it. They are functions
of the control net itself -- a flat-bottomed tether on every edge length,
$E = \tfrac k2 (\ell_{\min} - \ell)^2$ below $0.7\,\ell_0$ and $\tfrac k2(\ell - \ell_{\max})^2$ above
$1.75\,\ell_0$, a floor on every triangle altitude, and a wall on every dihedral crease -- with
stiffnesses of order $k = 83.4$ pN/nm.

An explicit step through a harmonic wall of stiffness $k$ along a direction of mobility $\mu$ is
stable only if $|1 - \mu k \Delta t| < 1$, i.e. $\mu k \Delta t < 2$. With $\mu$ up to $16 D/k_BT$ and
two or three walls pressing on the same vertex, $\Delta t \approx 0.003\ \mu s$ is already enough.
The measured onset agrees: just before a run diverges the mesh-term energy jumps from its usual
300-800 to 3400 pN·nm in one step, an effective stiffness of about 3600 pN/nm, $D k \Delta t/k_BT
\approx 3.5$.

The forensics of one divergence (frozen, $\ell = 5$, $\Delta t = 0.004$, seed 20261012) show the
mechanism directly. At step 6734 the force displacement of one vertex jumps from 0.3 to 1.3 nm while
its noise displacement is 0.19 nm, about $2\sigma$ with $\sigma = \sqrt{2D\Delta t} = 0.09$ nm; at
step 6735 one of its faces inverts, and the run grows by a factor 2-4 per step from there. In the
6700 steps before, the run absorbed a $5.2\sigma$ kick without incident. The trigger is a force
overshoot, not a large random number.

## 3. What does not work: truncating the noise

Since the divergence arrives at a random time, it is tempting to blame the tail of the Gaussian and
to draw $\xi$ from a truncated distribution instead (JohnsonBiophysicsLab/SLIMED#18). For the linear
part of the dynamics this cannot change the physics, provided the variance is kept: the
Euler-Maruyama recursion for a mode of rate $\Gamma$,

$$
h_{t+\Delta t} = (1 - \Gamma \Delta t)\, h_t + \sqrt{2 D_h \Delta t}\;\xi_t ,
$$

has stationary second moments that depend on the noise only through $\operatorname{Var}\xi$. With
the variance renormalised, a cut at $c\sigma$ reproduced $S(q)$ and $\Gamma(q)$ to 1-2 %, as it must.
It also barely helped: at $2\sigma$ the divergence rate fell by a factor 0.58, which, since the rate
rises as $\Delta t^{12.6}$, is a 4 % gain in step. Without the renormalisation the run is simply
colder, $S \to c_v S$ with $c_v = 0.774$ at $2\sigma$ while $\Gamma$ is unchanged, which breaks the
fluctuation-dissipation product $S\,\Gamma$; the same stability gain is obtained by lowering the
temperature by the same factor. Both outcomes point to the force term, not the noise.

## 4. Multiple time stepping of the walls

The force is split into the slow part, bending plus area plus any bead, and the stiff local part
$F_W$, the three walls. With `meshForceSubsteps` $= n > 1$, one Brownian step becomes $n$ substeps.
The first carries the whole slow force, a $1/n$ share of the wall force evaluated at the start of the
step, and the first share of the noise:

$$
S^{(1)} = S^{t} + \frac{D\Delta t}{k_BT}\Big[F_{\rm slow}(S^t) + \tfrac1n F_W(S^t)\Big]
+ \sqrt{\frac{2D\Delta t}{n}}\;\xi^{(0)} .
$$

The remaining $n - 1$ re-evaluate only the walls, on the control net of the current substep:

$$
S^{(k+1)} = S^{(k)} + \frac{D\Delta t}{n\,k_BT}\, F_{W}\!\left(S^{(k)}\right)
+ \sqrt{\frac{2D\Delta t}{n}}\;\xi^{(k)}, \qquad k = 1, \dots, n-1,
\qquad S^{t+\Delta t} = S^{(n)} .
$$

Every substep goes through the same pipeline as a full step: the wall force is evaluated on
$C = M^{-1}S$, carried to the surface as $M^{-T}F_W$ (two conjugate-gradient solves on the sparse
mask), the mean height is pinned, and the periodic images are re-asserted. The walls are closed-form
sums over the edge table, so a substep costs 0.72 ms against 3.08 ms for a full frozen step on the
$20 \times 24$ tile; with flips, where a step costs about 27 ms in the irregular-patch kernels, the
substeps do not show at all.

The $n$ noise draws are keyed per substep from the counter-based generator. Their sum is a single
Gaussian of the full variance, $\sum_k \sqrt{\Delta t / n}\,\xi^{(k)} \sim \mathcal N(0, \Delta t)$,
and the slow force is held at its start-of-step value, so every bending mode receives, in law,
exactly the step it received before; only the walls are now integrated with an Euler-Maruyama step
of length $\Delta t / n$. This is the `split` noise mode and the default.

The alternative, `once`, keeps the original single draw on the first substep. It looks more
conservative but samples the walls wrongly. A wall violation $x$ is refreshed by the noise once per
full step but relaxed $n$ times, by a factor $(1 - \mu k \Delta t / n)^n \approx e^{-\mu k \Delta t}$,
so at the step boundary, where frames and energies are taken, the walls sit at almost no violation.
For a stiff harmonic coordinate the stationary variances of the two schemes are

$$
\sigma^2_{\rm once} = \frac{2 D \Delta t}{1 - a^2}, \quad a = \left(1 - \frac{\mu k\Delta t}{n}\right)^{n},
\qquad
\sigma^2_{\rm split} = \frac{k_BT}{k}\,\frac{1}{1 - \mu k \Delta t / (2n)} ,
$$

the second tending to the exact $k_BT/k$ as $n$ grows and the first not. Measured on the frozen
lFace-5 tile at $\Delta t = 0.004$: the wall energy reads 82 pN·nm with `once` and $n = 2$ (55 at
$\Delta t = 0.005$, 20 at $0.008$ with $n = 8$) against 262 at $\Delta t = 0.0005$, while `split`
gives 297 at $n = 2$ and 260 at $n = 4$, the plain step's value at the substep size. `Once` also moves
the long waves ($S$ +4-5 %, $\Gamma$ −4.5 % at $q = 0.12$-$0.3$ nm⁻¹) and needs about twice the
substeps to be stable. It is kept only for comparison.

## 5. Stability with substeps

The two sources of stiffness now separate. The walls need the substep under the plain step's own
bound $\delta_W(\ell)$, and bending needs the full step under the bound of an explicit step on its
fastest mode, $\Gamma_{\max}\Delta t < 2$:

$$
\frac{\Delta t}{n} \le \delta_W(\ell), \qquad \Delta t \lesssim 0.03 \left(\frac{\ell}{5}\right)^{2} \mu s .
$$

The wall bound was measured with the plain step on a tile 20 faces across at each $\ell$: 0.0025 µs
at $\ell = 5$, 0.005 µs at $\ell = 10$ (one divergence in 1.6 ms), 0.006 µs at $\ell = 20$ (none in
1.9 ms, but about one per $10^5$ steps once used as a substep), so $\delta_W = 0.0025$, 0.004 and
0.005 µs are the values to use. It grows only weakly with $\ell$, because the walls are stiff in
pN/nm whatever the mesh. The bending bound grows as $\ell^2$ because $\Gamma_{\max} \propto D a
\kappa_c q_{\max}^4 / k_BT$ with $a \propto \ell^2$ and $q_{\max} \propto 1/\ell$; it is the
height-only bound of `docs/resolution_limits.md` section 2.2 (0.03-0.04 µs at $\ell = 5$, 0.12-0.16
at $\ell = 10$).

With $n = \lceil \Delta t / \delta_W \rceil$ and split noise, no run diverged on the frozen tile at
$\Delta t = 0.005, 0.01, 0.02, 0.03$ µs ($\ell = 5$, 0.8 ms each), at $0.02$-$0.12$ µs ($\ell = 10$,
3.2 ms each), or at $0.08, 0.16, 0.32$ µs ($\ell = 20$, 25.6 ms each with $\delta_W = 0.005$); at
$\ell = 20$ the step $0.48 = 0.03 \times 4^2$ µs diverges within three steps, as the bending bound
says it should. On the fluid lFace-5 tile, where the plain step's bound is 0.002 µs and the substep
bound about 0.0025 µs, $n = 2$ at 0.004 and 0.005 µs, $n = 4$ at 0.008 and 0.01 µs and $n = 8$ at
0.02 µs ran 320 µs each without a divergence, while the plain step at 0.004 µs diverged in all eight
seeds (409 per ms). Against the old frozen bounds these stable steps are 12 ($\ell = 5$), 24
($\ell = 10$) and 53 ($\ell = 20$) times larger.

## 6. Accuracy: the Euler-Maruyama law

Stable is no longer the same as accurate, because bending is still integrated with the full step.
For a linear mode of rate $\Gamma$ the recursion of section 3 has the stationary variance and lag
correlation

$$
\frac{\langle |h|^2 \rangle_{\Delta t}}{\langle |h|^2 \rangle_{0}} = \frac{1}{1 - \Gamma\Delta t/2},
\qquad
\frac{\langle h_{t+m\Delta t} h_t\rangle}{\langle |h|^2\rangle} = (1 - \Gamma\Delta t)^m
\;\Rightarrow\;
\Gamma_{\Delta t} = -\frac{\ln(1 - \Gamma\Delta t)}{\Delta t},
$$

so both the static spectrum and the apparent relaxation rate are too large by
$\approx \Gamma\Delta t / 2$. For the flat membrane the rate of a mode is the free-draining one,
verified mode by mode in `docs/resolution_limits.md` section 3.2,

$$
\Gamma(q) = \frac{D\, a\, \kappa_c\, q^4}{k_BT}, \qquad a = \frac{\sqrt 3}{2}\,\ell^2 ,
$$

so the error of a mode depends on $q$, $\ell$ and $\Delta t$ only through $\Gamma(q)\Delta t$. At the
edge of the window the mesh resolves, $q = 2/\ell$,

$$
\Gamma(2/\ell) = \frac{8\sqrt3\, D\,\kappa_c}{k_BT\,\ell^2} = 11.1 \left(\frac{5}{\ell}\right)^2 \mu s^{-1},
$$

and a 5 % tolerance, $\Gamma\Delta t \lesssim 0.1$, gives $\Delta t \approx 0.009$-$0.01\,(\ell/5)^2$
µs. More generally, for the largest wave number one cares about,

$$
\Delta t \;\lesssim\; \frac{0.1\, k_BT}{D\, a\, \kappa_c\, q_{\max}^4} .
$$

To test this the benchmark used tiles of 20 faces at $\ell = 5$, 10 and 20 nm, so that a step
$\Delta t$ at face length $\ell$ is the same lattice problem for bending as $\Delta t\,(5/\ell)^2$ at
$\ell = 5$. Long runs of $400\,(\ell/5)^2$ µs with a frame every $0.2\,(\ell/5)^2$ µs gave $S(q)$ and
the slow rates; short runs of $40\,(\ell/5)^2$ µs with a frame every $0.02\,(\ell/5)^2$ µs gave the
fast rates up to $q\ell \approx 2.5$. Each was compared with a reference at $0.001\,(\ell/5)^2$ µs of
the same $\ell$, two seeds each. The three face lengths agree bin by bin, so they are pooled; the
table quotes the step at $\ell = 5$.

| $\Delta t$ at $\ell = 5$ (µs) | $\Gamma(2/\ell)\Delta t$ | | $q\ell$ 1-1.5 | 1.5-2 | 2-2.5 | 2.5-3.5 |
| --- | --- | --- | --- | --- | --- | --- |
| 0.0025 | 0.028 | $S$ | 1.003 ± 0.003 | 1.005 ± 0.003 | 1.019 | 1.044 |
| | | $\Gamma$ | 1.035 ± 0.015 | 1.014 ± 0.007 | 1.023 | |
| 0.005 | 0.055 | $S$ | 1.009 ± 0.004 | 1.016 ± 0.003 | 1.045 | 1.125 |
| | | $\Gamma$ | 1.011 ± 0.021 | 1.032 ± 0.008 | 1.034 | |
| 0.01 | 0.111 | $S$ | 1.021 ± 0.006 | 1.042 ± 0.004 | 1.105 | 1.338 |
| | | $\Gamma$ | 1.019 ± 0.020 | 1.047 ± 0.006 | 1.086 | |
| 0.02 | 0.222 | $S$ | 1.011 ± 0.005 | 1.081 ± 0.004 | 1.284 | 2.48 |
| | | $\Gamma$ | 1.045 ± 0.009 | 1.141 ± 0.005 | 1.255 | |

![Error of S(q) and Γ(q) against Γ(q)Δt at every face length, with the Euler-Maruyama curves](figures/substep_resolution.png)

Plotted against $\Gamma(q)\Delta t$, every bin inside $q\ell \le 2$ at every $\ell$ falls on the two
curves $1/(1 - x/2) - 1$ and $-\ln(1-x)/x - 1$. At $\Delta t = 0.01\,(\ell/5)^2$ the worst bin inside
the window is +4.2 % in $S$ and +4.7 % in $\Gamma$; at twice that step it is +8 % and +14 %. Beyond
$q\ell = 2$ the error grows faster than the curve, because the limit surface stiffens those modes
and $\Gamma$ underestimates their rate, but they lie outside what the mesh resolves anyway. The
lowest modes ($q\ell < 1$, a handful with relaxation times of tens to hundreds of µs) scatter by
about ±3 % at every step and carry no trend.

Two controls separate the scheme from the step. At a fixed step the substeps change nothing:
at $\Delta t = 0.0025$, $\ell = 5$, split with $n = 2$ against the plain step (three seeds each,
400 µs) gives $S$ ratios of $0.996 \pm 0.005$, $0.993 \pm 0.002$ and $0.991 \pm 0.014$ in the bands
$q = 0.12$-$0.3$, $0.3$-$0.6$ and $0.6$-$1.0$ nm⁻¹ and a rate ratio of $1.007 \pm 0.016$. And the
error at larger steps continues the plain step's own trend from below its old bound. The fluid sheet
behaves the same way with a slightly larger error: at $q\ell = 1.5$-$2$ it reads $S = 1.032$,
$\Gamma = 1.059$ at $\Delta t = 0.005$ µs and $S = 1.069$, $\Gamma = 1.078$ at 0.01 µs, so its
accurate step at $\ell = 5$ is about 0.005 µs.

## 7. The mesh size

Substepping does not touch the spatial discretisation. On the 20-face tiles, the bending rigidity
each mode reads, $\kappa_{\rm app}(q) = k_BT / (A\, S(q)\, q^4)$, of the $\ell = 10$ and $20$
references against the $\ell = 5$ reference at the same $q$ is 1.014 and 1.032 for $q\ell = 1$-$1.5$
and $1.5$-$2$, 1.07 at $2$-$2.5$, 1.13 at $2.5$-$3.5$ and 1.4-1.5 beyond: the curve of
`docs/resolution_limits.md` section 3.1. A coarse mesh is a fine mesh truncated, and the bound is

$$
q\,\ell \le 2 \ \ (5\ \%), \qquad q\,\ell \le 2.5 \ \ (10\ \%), \qquad\text{i.e.}\quad
\ell \le \frac{\lambda_{\min}}{\pi}.
$$

The temporal and spatial errors are separate: the first depends on $\Gamma(q)\Delta t$, the second on
$q\ell$, and the step bound of section 6 is stated inside the spatial window.

## 8. Choosing a run

Given the shortest wavelength that matters, $\lambda_{\min}$, set $q_{\max} = 2\pi/\lambda_{\min}$
and take the face length $\ell \le \lambda_{\min}/\pi$. The accurate step is then
$\Delta t \le 0.1\, k_BT / (D a \kappa_c q_{\max}^4)$, which at $q_{\max} = 2/\ell$ is
$0.01\,(\ell/5)^2$ µs for $D = 1$ nm²/µs; for another diffusion constant only the product
$D\,\Delta t$ matters. The substep count is $n = \lceil \Delta t / \delta_W(\ell) \rceil$ with
$\delta_W = 0.0025$, 0.004 and 0.005 µs at $\ell = 5$, 10 and 20 nm on the frozen mesh, and
0.0025 µs on the fluid lFace-5 mesh. For $\ell = 5$:

```
timeStep = 0.01
meshForceSubsteps = 4
meshForceSubstepNoise = split    # the default
surfaceSolver = iterative
fdtConsistentSurfaceUpdate = true
edgeSpringEnabled = true          # the walls must be in the drive
```

The elastic box has its own clock and its own explicit step. Its strain relaxes at
$r = 2 K_A \mu_{\rm box}$, and the benchmark held $r\,\Delta t \le 0.1$ by setting
`elasticBoxMobility` $= \min\!\big(0.02\,(5/\ell)^2,\ 2\times10^{-4}/\Delta t\big)$ nm/(pN µs), with
$K_A = 250$ pN/nm; at large steps that cap
must be applied by hand. The run refuses `meshForceSubsteps > 1` unless the FDT-consistent update,
the edge spring, soft walls and the Brownian step are all on, since the scheme only adds up to the
full wall force under those conditions.

The cost per step is $t_{\rm step} + (n - 1)\, t_{\rm sub}$ with $t_{\rm step} = 3.08$ ms and
$t_{\rm sub} = 0.72$ ms on the frozen tile, so simulated time per wall-clock second improves on the
plain step at its old bound $\delta_0$ by

$$
G = \frac{\Delta t / \left[t_{\rm step} + (n-1)\,t_{\rm sub}\right]}{\delta_0 / t_{\rm step}},
$$

which is 2.4, 2.6 and 3.2 at $\ell = 5$, 10 and 20 for the accurate steps above. On the fluid mesh
$t_{\rm sub} \ll t_{\rm step}$ and $G$ is the step ratio.

## 9. Verification

With `meshForceSubsteps = 1` the code path is the shipped one: the legacy periodic workload (300
steps) and the Mixed frozen tile (2000 steps) produce byte-identical output files against
`origin/main`, 21 files in all. `tests/test_mesh_force_substeps.cpp` checks that an unsupported
configuration is refused, that with zero wall force the `once` scheme at $n = 4$ reproduces the
shipped step bit for bit, and that `split` changes the draws but keeps the per-step variance at
$2D\Delta t$. The full suite passes (408 tests).

## 10. What is left

The remaining limit is the Euler-Maruyama step on bending, so a further gain in accurate step needs
a better integrator for the slow force -- a stochastic Heun or exponential step on the bending
modes, for instance -- not more work on the walls. Three things were not measured: the fluid mesh
at $\ell = 10$ and 20, the substep bound under flips beyond $\ell = 5$, and whether the spectral
tension of the slowest modes, which needs runs of several ms (`docs/resolution_limits.md` section
3.3), depends on the step.
