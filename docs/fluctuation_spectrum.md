# Verifying the membrane against `<|h_q|^2> = kT / (A kc q^4)`

A flat fluid membrane in the Monge gauge has the Helfrich free energy

```
E = (1/2) ∫ dA [ kc (∇²h)² + σ (∇h)² ]
```

which for a periodic box of projected area `A`, with the convention
`h(r) = Σ_q h_q e^{iq·r}` and `h_q = (1/A) ∫ d²r h(r) e^{-iq·r}`, diagonalises to
`E = (A/2) Σ_q (kc q⁴ + σ q²) |h_q|²`. Equipartition then gives the relation this
document is about:

```
<|h_q|²> = kB T / ( A (kc q⁴ + σ q²) )                                     (1)
```

With no tension that is the pure `q^-4` law — a straight line of slope −4 on a
log-log plot whose intercept fixes `kc`. It is the standard end-to-end check on
a membrane model, and it is worth doing precisely because it is sensitive to
almost everything: the geometry of the mesh, the Fourier conventions, the
Brownian integrator, and the boundary condition all have to be right at the
same time before the line comes out straight.

This document records why the check was failing, what was wrong on each side,
and what the corrected pipeline measures. The second half covers the companion
question — what `<h_q(t) h_q*(t+τ)>` can and cannot tell you about hydrodynamics
in this code.

Code:

| file | what it is |
| --- | --- |
| `analysis/membrane_spectrum.py` | the library: geometry, transform, q grid, binning, fits, autocorrelation |
| `analysis/fluctuation_report.py` | per-run report; also recomputes the legacy reduction for comparison |
| `analysis/reference_langevin.py` | an independent spectral integrator with a closed-form answer, for validation and for the Oseen comparison |
| `analysis/fluctuation_plots.py` | figures from the `--npz` output |
| `analysis/test_membrane_spectrum.py` | self-tests on synthetic data whose spectrum is known exactly |

---

## 1. What the mesh actually is

Everything downstream depends on getting this right, and it is not what a
`fft2` of the z column assumes.

`Mesh::set_vertices_faces_flat()` builds a **zig-zag triangular lattice**:

```
x(i, j) = i·dFaceX + (dFaceX/2 if j is even else 0) − lxHalf
y(i, j) = j·dFaceY − lyHalf,          dFaceY = (√3/2)·dFaceX
```

Vertex `(i, j)` is entry `j·nVertX + i` of a trajectory line. Two consequences:

**The sampling positions are not a Cartesian grid.** Alternate rows are offset
by half a cell. The lattice *is* a proper Bravais lattice — `a1 = (dx, 0)`,
`a2 = (dx/2, dy)` — but the map from array index to position is not affine, so
an index-space `fft2` is not the Fourier transform of the sampled field.

**Only part of the array is independent.** Under `boundaryType = Periodic`,
`determine_ghost_vertices_faces()` marks three rings of ghosts on each side, and
`DynamicMesh::get_relative_pt_periodic()` identifies the fourth ring with the
ring on the opposite edge. The periodic tile is therefore the
`(nFaceX−6) × (nFaceY−6)` block starting at index `(3, 3)`, with

```
Lx = (nFaceX−6)·dFaceX      Ly = (nFaceY−6)·dFaceY      A = Lx·Ly
```

`nFaceY` is forced even by `set_axes_division_flat()`, which is exactly what
makes the zig-zag close up under the `Ly` translation.

Two facts that make the analysis easier than expected, both verified against
the trajectories rather than assumed:

* **The in-plane coordinates never move.** `DynamicModel::next_step()` zeroes
  the x and y displacement, so every vertex stays on the ideal lattice for the
  whole run (drift `< 1e-11 nm` over 2000 steps). There is no need to resample
  a moving surface onto a grid.
* **`matSurface` sits at the same `(x, y)` as `matMesh`.** The Loop limit mask
  is affine-invariant and the six neighbour offsets of an interior vertex sum
  to zero, so applying it to a perfect lattice returns the lattice. The limit
  surface differs from the control net only in `z`.

## 2. The five things wrong with the original reduction

`analysis/fft_periodic_result_large(l50).py` gets a slope near −1.4 instead of
−4. Reproducing its reduction inside `fluctuation_report.py` (`legacy_reduction`)
lets the two be run on identical frames, and running it on *synthetic data drawn
from equation (1)* — where the right answer is known exactly — separates
"the analysis is wrong" from "the simulation is wrong". On that synthetic data
the legacy reduction returns **−1.43** and the corrected one returns **−4.00**.

**(a) The zig-zag phase is missing.** The true transform factorises as

```
H(q) = Σ_s e^{-i q_y s·dy} · e^{-i q_x·(dx/2)·[row s offset]} · Σ_p h_{p,s} e^{-i q_x p·dx}
```

— an FFT along x, a per-row phase, an FFT along y. `np.fft.fft2` drops the
middle factor. That is not a small correction. The omitted factor is a period-2
modulation in the row index, so it mixes every mode `n_y` with `n_y + Ny/2`:

```
H_naive[n_y] = cos(θ/2)·H_true[n_y] + i·sin(θ/2)·H_true[n_y + Ny/2],   θ = π n_x / Nx
```

At `n_x ≈ Nx/2` the two contributions are equal. Because the spectrum spans four
decades, this dumps low-q power straight into the high-q bins and flattens the
tail into a plateau. `tile_fft2()` puts the phase back; the self-test checks it
against an explicit sum over the true `(x, y)` and gets agreement to 8e-15.

There is a subtlety worth flagging, because it is easy to get wrong a second
time: the phase must use the **signed** mode number from `np.fft.fftfreq`, not
the raw output index. The row transform cannot distinguish `n_x` from
`n_x + Nx`, but the half-cell phase can, since
`exp(-iπ(n_x+Nx)/Nx) = -exp(-iπ n_x/Nx)`. Resolving `q_x` from `q_x + 2π/dx` is
the entire job of sampling a half-cell offset, and using the unsigned index
flips it — shifting the whole negative-`q_x` half of the spectrum by half the
zone in `n_y`, which pairs long-wavelength modes with short-wavelength ones.

**(b) The wrong block is transformed.** The original trims four rows and
columns from index 0 (`for y in range(0, y_scale-4)`), which is neither the
periodic tile nor aligned with it. Transforming a non-periodic block puts a
step discontinuity at the seam and leaks power across the whole spectrum. The
tile is `[3 : 3+Ny, 3 : 3+Nx]`.

**(c) The fold slices rows twice.** `arr_fft2[0:len(arr_fft2)//2][0:len(arr_fft2[0])//2]`
looks like it keeps the first quadrant. It does not: the second subscript
applies to the *result of the first*, so it is `arr[0:min(N//2, M//2)]` with
every column kept. The array that is then labelled with quarter-fold `|q|`
still contains the full column range.

**(d) `|q|` is built from the array index, not the reciprocal lattice.** Three
errors compound:

* the raw index is used instead of the signed one, so index `n_x > Nx/2` — which
  is the *short* wavevector `−q` — is labelled as a long one;
* the `2π` is missing in one version and present as `4π²` in another;
* `L_mean = (x_scale + y_scale)/2` is a count of vertices, not a length in nm,
  so `q` is not even dimensionally a wavevector.

The right thing is `q = 2π·(fftfreq(Nx, dx), fftfreq(Ny, dy))`, keeping every
mode and binning by `|q|`; conjugate pairs counted twice do a radial mean no
harm. `q_grid()` also folds each `q` to the shortest of its reciprocal-lattice
images, since a discrete mode is the whole family `q + G` and the continuum
formula refers to the smallest member.

**(e) The DFT is not normalised.** `np.fft.fft2` returns
`H_q = Σ_j h_j e^{-iq·r_j} = N·h_q` in the convention of equation (1), so the
measured quantity is `<|H_q|²> / N²`. Without this the fitted `kc` is wrong by
`N²` even when the slope is right.

A sixth item on the original list of suspicions turns out not to matter here:
the control points really are not on the limit surface, but `surfacepoint*.csv`
already holds the limit positions, and the spectra from `surfacepoint*.csv` and
from `meshpoint*.csv` put through the limit mask agree to four significant
figures. (`matSurface` is *not* exactly periodic at the tile seam, because
`assign_mesh2surface()` gives ghost-adjacent vertices a different row; the
difference reaches 0.18 nm against a tile standard deviation of 0.9 nm, but it
is confined to a few points and contributes ~1e-6 nm⁴ to the spectrum. The
report defaults to `--source control`, which reconstructs the limit surface with
periodic wrap and is exactly translation invariant.)

## 3. Three things wrong on the simulation side

Fixing the analysis is necessary but not sufficient. Pointed at an as-shipped
trajectory the corrected reduction returns a slope of −5.8 and a `kc` of 38
pN·nm, which is not an improvement so much as a different kind of wrong: as
section 3(b) shows, that run has no stationary spectrum for any reduction to
measure. The three defects below are what stands between a correct analysis and
a meaningful number.

**(a) The random number generator was raced, and unseeded.**
`DynamicModel::next_step()` called a shared `std::normal_distribution` on a
shared `std::mt19937` from inside `#pragma omp parallel for` — a data race on
the generator state, so the noise was neither reproducible nor guaranteed
Gaussian. Separately, the constructor re-declared `gen` and `normal_dist` as
locals, shadowing the members, so `randomSeed` never reached the generator that
was actually drawn from. Replaced with a counter-based stream keyed by
`(seed, iteration, vertex, axis)`, which is reproducible independently of thread
count and scheduling.

**(b) The Brownian step violated fluctuation-dissipation, and the run does not
converge.** This is the interesting one. `Run_dynamics_flat` runs

```
S = M C ;   S += (D dt/kT)·F + √(2 D dt)·ξ ;   C = M⁻¹ S
```

where `M` = `mesh2surface` is the Loop limit mask, `C` the control net and `S`
the limit surface. But `F` is `-∂E/∂C`, the gradient with respect to the
*control* points. In `C` the update is therefore

```
ΔC = (D dt/kT)·M⁻¹F_C + √(2 D dt)·M⁻¹ξ
```

whose mobility is `(D/kT)M⁻¹` but whose noise covariance is `2D·M⁻²`.
Fluctuation-dissipation wants `2kT·mobility = 2D·M⁻¹`. Working mode by mode
(both `M` and the stiffness are diagonalised by plane waves on this lattice),
the stationary variance comes out

```
<|h_q|²>_measured = [ kT / (A kc q⁴) ] · 1/m(q)
```

where `m(q) = 1/2 + (1/6)[cos(q·a1) + cos(q·a2) + cos(q·a3)]` is the Fourier
symbol of the limit mask — 1 at long wavelength, falling to about 0.25 at the
zone corner. So every mode is too warm, by a factor that grows with `q`: exactly
the shape that lifts the tail off the `q^-4` line.

The gradient with respect to `S` is `M⁻ᵀF_C`, so one application of
`surface2mesh` fixes it — exactly where `M` is symmetric, which every interior
row is, since the valence-6 limit mask is. Controlled by
`fdtConsistentSurfaceUpdate` (default `true`; set `false` to reproduce older
results).

`M` is **not** symmetric everywhere, and that turns the bias into something
worse. `assign_mesh2surface()` gives a vertex whose adjacent faces are all
ghost or boundary the identity row, while its neighbours' rows still carry its
`1/12`; the code now reports `max|M − Mᵀ|` at startup, and on the default
periodic mesh it is exactly `1/12`. A non-symmetric `M` makes the mobility
*non-reciprocal*, and a non-reciprocal mobility is not the gradient flow of any
potential — it can do net work on the membrane instead of relaxing it. That is
what the trajectories show:

Running each defect on its own separates them. The r.m.s. tile height, in the
first and last eighth of each run:

| run | first eighth | last eighth | verdict |
| --- | ---: | ---: | --- |
| as shipped, seed 101 | 2.13 nm | 19.10 nm | runs away |
| as shipped, seed 31337 | 2.27 nm | 25.67 nm | runs away |
| FDT defect only | 0.84 nm | 0.86 nm | large slow excursions, no runaway |
| duplicates defect only | 1.42 nm | 1.11 nm | stationary |
| all three fixed | 0.88 nm | 0.74 nm | stationary |

The runaway is reproducible across seeds, and it needs both defects. On its own
the fluctuation-dissipation defect gives a spectrum that is too warm by `1/m(q)`
and wanders in large slow excursions but comes back; on its own the duplicated
noise is a bounded bias of about 15% in amplitude. Together, the extra heat
feeds the non-reciprocal circulation and the amplitude does not return.

So the as-shipped update does not have a stationary spectrum to measure at all;
anything read off such a run records where it happened to be when it was
stopped. `fluctuation_report.py` now prints this check before anything else.

The remaining approximation is that `surface2mesh` is `M⁻¹`, standing in for
`M⁻ᵀ`. The rows where they differ all lie in the ghost ring that
`postprocess_ghost_periodic()` overwrites anyway, and the measured
`S(q)K(q)/kT` lands on 1 to within a few percent, but it is an approximation
and is flagged as one in the code.

**(c) The periodic duplicates were given their own noise.** The fourth ring in
from each edge is overwritten from its partner by
`postprocess_ghost_periodic()` every step, so it is not an independent
coordinate. It was nevertheless integrated first — and `surface2mesh` is dense,
so its fresh random kick was smeared over the whole mesh *before* being
overwritten. On the 100 nm test box that is 64 extra random numbers per step
against 252 free coordinates: extra heat with no extra friction. Controlled by
`integratePeriodicDuplicates` (default `false`).

The three effects were separated by running the same 60 000-step trajectory
three ways and looking at the internal FDT product described below.

## 4. An internal check that does not assume the continuum limit

The `q^-4` law only holds where the discretisation is a faithful stand-in for a
smooth surface, so a deviation at large `q` is ambiguous: bug, or lattice? The
following test removes the ambiguity, because it never mentions `kc` or `q⁴`.

Each mode of the trajectory yields two independent numbers: its variance `S(q)`
and its own relaxation rate, read off the one-step autocorrelation as
`λ(q) = 1 − C_q(1 step)`. For an overdamped update `h ← (1−λ)h + noise` with
mobility `μ` and stiffness `K`, `λ = μKdt` and the stationary variance is
`(kT/K)·2/(2−λ)`. Reconstructing `K` from `λ` through the scheme's own nominal
mobility and forming

```
S(q) · K(q) / kT / [2/(2−λ)]
```

must give **1** for any correct Brownian update, whatever the true stiffness is.
Measured on the three variants (60 000 steps, output every step, 100 nm box):

| q (nm⁻¹) | m(q) | as shipped | + FDT fix | + duplicates fix |
| ---: | ---: | ---: | ---: | ---: |
| 0.242 | 0.834 | 1.39 | 1.16 | 0.97 |
| 0.278 | 0.785 | 1.47 | 1.15 | 0.96 |
| 0.347 | 0.687 | 1.64 | 1.13 | 0.95 |
| 0.421 | 0.582 | 1.87 | 1.09 | 0.93 |
| 0.501 | 0.475 | 2.22 | 1.05 | 0.91 |
| 0.605 | 0.368 | 2.70 | 0.99 | 0.88 |
| 0.720 | 0.298 | 3.17 | 0.94 | 0.85 |

The as-shipped column tracks `1/m(q)` (1.20, 1.27, 1.46, 1.72, 2.11, 2.72, 3.36)
to within a constant; multiplying it by `m(q)` flattens it to 1.16 … 0.94, which
is the middle column. That is the signature predicted above, measured rather
than argued. The third column is flat at ≈1 — the update samples `exp(-E/kT)`.

The same reconstruction also gives the **discrete stiffness** `K_S(q)`, which is
what says how far up in q the continuum law should be trusted. It is too noisy
to settle that on its own, though — at small q it inherits the same
finite-record bias as the rates it comes from. The next section measures the
same thing statically, and much more precisely.

## 4b. How far up in q the continuum law holds — measured, not assumed

`--qdx-max` needs a justification better than a round number, and the
published cutoffs do not supply one. Brandt, Braun, Sachs, Nagle and Edholm
(*Biophys J* **100**, 2104, 2011) put the limit near `q = 0.7 nm⁻¹` for DMPC,
but their limit is molecular: the density structure factor, thickness
fluctuations, protrusion modes, with an explicit criterion
`q₀ = (K_A k_c ⟨h²⟩)^¼` that has a thickness in it. SLIMED has no thickness, no
tilt and no protrusions — it *is* a Helfrich surface by construction — so
nothing physical limits it and that number does not transfer. Shiba and Noguchi
(*Phys Rev E* **84**, 031926, 2011) make the complementary point: a fitted
bending rigidity depends strongly on where the cutoff is put, and they
recommend fitting as a function of `q_cut` and extrapolating to `q_cut → 0`.

What limits SLIMED is the mesh, and the mesh can be asked directly.
`tests/test_continuum_limit.cpp` puts an exact plane wave `h = ε cos(q·r)` on
the control net, one allowed wavevector at a time, reads the bending energy the
production code returns, and forms `K = 4E/ε²`. The energy is quadratic in the
*control* amplitude while the continuum law is written for the surface, and the
two differ by the limit-mask symbol, so the quantity to compare is
`K_S/(A kc q⁴)` with `K_S = K/m(q)²`. No sampling, no dynamics, no fitting.

| `q·dFaceX` | λ / dFaceX | `K_S/(A kc q⁴)` |
| ---: | ---: | ---: |
| 0.40 | 15.6 | 1.0000 |
| 0.60 | 10.4 | 0.9999 |
| 0.90 | 7.0 | 0.9992 |
| 1.21 | 5.2 | 0.9979 |
| 1.35 | 4.7 | 0.9953 |
| 1.61 | 3.9 | 0.9919 |
| 1.81 | 3.5 | 0.9845 |
| 2.01 | 3.1 | 0.9754 |
| 2.24 | 2.8 | 0.9497 |
| 2.42 | 2.6 | 0.9355 |
| 2.82 | 2.2 | 0.8509 |
| 3.23 | 2.0 | 0.7002 |

So the Loop limit surface reproduces the continuum bending stiffness to **0.5%
at `q·dFaceX = π/2`**, 3% at 2.0 and 7% at 2.5. The default window is
conservative, not marginal, and a wider one is defensible if the extra reach in
q is worth a few per cent. (An earlier estimate of "within 10% out to 2" came
from the dynamic reconstruction above; the static probe is twenty times more
precise and supersedes it.)

Two things fall out that are worth keeping:

* **That the ratio comes out at 1.000 rather than some other constant is a
  check on `m(q)` itself.** Nothing was fitted; the mask symbol was derived
  from the valence-6 limit weights and simply divides the measurement to 1.
* **The softening is a function of `q·dFaceX`, not of `q`.** A 100 nm box at
  `lFace = 5` and an 85 nm box at `lFace = 2.5` both tile to `Lx = 70 nm`, so
  mode 3 of the first and mode 6 of the second sit at exactly the same
  `q·dFaceX` at twice the physical q on half the mesh. They give `0.9953` and
  `0.9953`. The window therefore has to be re-derived per mesh resolution, not
  per box size.

The window is set at `q·dFaceX ≤ π/2` rather than at some decimal near it
because `q·dFaceX = π` is the Nyquist wavevector of the mesh spacing: `π/2` is
exactly half of it, `λ = 4·dFaceX`, four cells per wavelength. That is also
roughly where the usual finite-element points-per-wavelength rule would land
for a C²-continuous quartic element — but neither of those is the argument.
The table is; the round number just happens to sit in a sensible place.

Finally, this closes the loop on the measured spectrum. At `q·dFaceX = 2.5` the
static form factor predicts `⟨|h_q|²⟩` should sit `1/0.93 = 1.07` above the
continuum line; measured, 1.11. At 3.0 it predicts 1.30; measured, 1.34. The
high-q lift in the spectrum is the discretisation, quantitatively.

## 5. Results

A 14x18 periodic tile (Lx = 70 nm, Ly = 77.94 nm, A = 5456 nm^2), kc = 83.4
pN.nm = 20 kT, `usMembraneStretching = 0`, dt = 0.002 us, D = 1 nm^2/us,
818,724 steps, a frame every 100. The reference column is not a fit: it is
`kT/(A kc q^4)` evaluated with the kc that went into the run.

```
 q (1/nm)   q*dx    <|h_q|^2>  kT/(A kc q^4)   ratio  modes  s.e.   T/tau_q
   0.0851   0.43   1.4557e-01     1.7503e-01   0.832      4  5.9%      29.7
   0.1206   0.60   5.1114e-02     4.3257e-02   1.182      4  1.3%     120.2
   0.1612   0.81   1.2012e-02     1.3563e-02   0.886      2    --      383.3
   0.1883   0.94   6.5006e-03     7.2897e-03   0.892     10  5.2%     713.1
   0.2415   1.21   2.3659e-03     2.6954e-03   0.878      6  1.9%    1928.6
   0.2780   1.39   1.4414e-03     1.5336e-03   0.940     14  6.4%    3389.6
   0.3474   1.74   6.1444e-04     6.2940e-04   0.976     24  5.1%    8259.3  (beyond fit)
   0.4207   2.10   2.9750e-04     2.9246e-04   1.017     24  3.2%   17774.5  (beyond fit)
   0.5006   2.50   1.6210e-04     1.4598e-04   1.110     40  3.0%   35610.1  (beyond fit)
   0.6049   3.02   9.1370e-05     6.8429e-05   1.335     62  1.9%   75967.8  (beyond fit)
   0.7197   3.60   6.6525e-05     3.4166e-05   1.947     61  0.9%  152152.3  (beyond fit)

log-log slope over q*dx <= pi/2: -4.022        (want -4)
fitted kc                     :  89.03 pN.nm   (21.3 kT; input 83.40, 20.0 kT)
fitted sigma                  :   0.027 pN/nm  (the run has no area elasticity)
legacy reduction, same frames :  -1.402
```

**Run length is what limits this, not the seed.** The lowest-q bin relaxes as
`q^-4`, so it is the last to converge, and a run that has not sampled it enough
reads low and drags the slope towards zero. Truncating the same trajectory
makes that explicit:

| trajectory | length | T/tau of the lowest bin | slope | fitted kc |
| --- | ---: | ---: | ---: | ---: |
| pure_a, full | 1637 us | 29.7 | **-4.02** | 89.0 pN.nm |
| pure_a, first quarter | 385 us | 7.0 | -3.57 | 76.2 pN.nm |
| pure_b, independent seed | 570 us | 10.3 | -3.25 | 62.3 pN.nm |

The report flags any bin with `T/tau_q < 10` as UNDER-SAMPLED for exactly this
reason. Since `tau ~ L^4` at fixed dt, the cost of the check grows steeply with
box size — see the performance note at the end.

**With area elasticity on, the law is q^-2, not q^-4.** The shipped
`input.params` sets `usMembraneStretching = 250` pN/nm, and
`run_dynamics_flat()` pins the reference area to the flat starting area.
Fluctuations buy excess area that the constraint pulls back, so a real tension
develops. The same box run that way:

```
fitted sigma = 3.74 pN/nm    crossover at q = sqrt(sigma/kc) = 0.212 /nm
measured slope over the fit window = -2.635
```

The crossover sits in the middle of the measurable band (`q` runs from 0.085 to
0.84 /nm here), so most of the spectrum is tension dominated. Nothing in the
analysis can recover `q^-4` from that; the area term has to be off, or the
reference area set to the fluctuating area rather than the flat one.

## 6. The relaxation time, and why Oseen cannot appear

The static spectrum is blind to hydrodynamics: equation (1) has no mobility in
it. All the information about the solvent is in the *rates*.

```
free draining (local mobility Λ)   Γ(q) = Λ (kc q⁴ + σ q²)
Oseen / Zimm (Λ(q) = 1/4ηq)        Γ(q) = kc q³ / 4η + …
```

**The code has no hydrodynamics.** `DynamicModel::next_step()` displaces each
vertex by `(D dt/kT)·F_i + √(2D dt)·ξ_i`: a scalar mobility, per vertex,
independent of every other vertex. There is no Oseen tensor, no
Rotne–Prager–Yamakawa, no mobility matrix anywhere in the repository. So the
implemented dynamics is free draining and `Γ ∝ q⁴` by construction; a `q^-3`
decorrelation cannot come out of it no matter how the trajectory is analysed.

Matching the discrete update to the continuum form gives the local mobility in
terms of the code's own parameters,

```
Λ_local = D·a/kT ,      a = A/N = area per vertex
Γ(q)    = (D·a/kT)·kc·q⁴
```

What the corrected SLIMED runs actually measure, with no fitted prefactor:

```
median Gamma / (D a / kT) kc q^4   =  1.00     (frame every step, q = 0.17 .. 0.75)
median Gamma / (kc q^3 / 4 eta)    =  0.006    (i.e. 170x too slow for water)
log-log slope of Gamma(q)          =  +3.6     over q*dFaceX <= pi/2
```

The prefactor lands on 1 because the free-draining rate is not a free parameter
here — `Gamma = (D a / kT) kc q^4` is fixed by the run's own `diffConst`, area
per vertex and `kcMembraneBending`. Measuring it correctly is a strong check on
the integrator; measuring it as `q^3` was never possible.

Two sampling rates are needed to cover the range, because `tau ~ q^-4` spans
four orders of magnitude across the zone: a trajectory written every 100 steps
cannot time a mode that decorrelates in three, and one written every step is
too short to time the slowest. `relaxation_rate()` returns `nan` rather than
the sampling rate when a mode falls outside what the record resolves — without
that guard `Gamma(q)` saturates at `1/dt_frame` and the measured power law
flattens towards `q^3` for entirely artificial reasons.

`analysis/reference_langevin.py` integrates the same model spectrally with
either mobility, on the same lattice, and confirms both limbs — identical static
spectra, different dynamics:

| | static slope | Γ slope | Γ at q = 0.085 nm⁻¹ |
| --- | ---: | ---: | ---: |
| free draining | −3.98 | **3.84** (exact 4) | 0.095 µs⁻¹ |
| Oseen, η = water | −3.94 | **2.93** (exact 3) | 15.7 µs⁻¹ |
| fitted `kc` | 80.1 / 80.2 pN·nm (input 83.4) | | |

Two further points worth knowing before reading any timescale out of this model:

* **`D` is not connected to a solvent.** It is a free per-vertex parameter in
  nm²/µs. Nothing ties it to a viscosity, so the absolute relaxation times are
  arbitrary until it is calibrated.
* **The free-draining rate is mesh dependent.** `Γ = (D·a/kT)·kc·q⁴` is
  proportional to the area per vertex, so refining the mesh at fixed `D` slows
  the dynamics by the refinement factor squared. To make the dynamics
  mesh-convergent, the quantity to hold fixed is the mobility *per unit area*,
  `Λ = D·a/kT`, i.e. `D ∝ 1/a`.

### What implementing Oseen would take

For a flat periodic sheet the Oseen mobility is diagonal in Fourier space,
`Λ(q) = 1/(4ηq)`, so the natural implementation is not a dense `N×N` mobility at
all:

1. FFT the nodal force to `F_q`;
2. multiply by `Λ(q)·dt`;
3. add noise with variance `2kT·Λ(q)·dt/A` per mode — trivially "square-rooted"
   because the mobility is diagonal, which is the step that is expensive in real
   space;
4. inverse FFT.

That is `O(N log N)` per step and is what `reference_langevin.py` does. A real
space Oseen/RPY implementation would instead need a dense mobility multiply plus
a Cholesky or Chebyshev factorisation for the correlated noise; the mobility
multiply is not a new asymptotic cost here (see below), but the noise
factorisation would be.

## 7. Performance note

`DynamicMesh::setup_flat()` inverts `mesh2surface` densely (`O(N³)`) and the
dynamics loop multiplies by the dense `surface2mesh` every step (`O(N²)`). On
the 100 nm box (525 vertices) that is already the dominant cost, and it does not
thread — the 140 nm box (957 vertices) runs 5.7× slower for 1.8× the vertices.

This is what caps the accessible `q` range: the width of the usable window is
set by the number of vertices per side, the slowest mode relaxes as `L⁴`, and
the cost per step goes as `N²`, so a box twice as wide costs about 250× more.

`mesh2surface` is sparse (7 non-zeros per row) and, on this lattice, very well
conditioned — its eigenvalues `m(q)` lie in `[0.25, 1]`. Applying `M⁻¹` by a
conjugate-gradient solve against the sparse `M` instead of by a dense multiply
would converge in a handful of `O(N)` iterations and remove both the `O(N³)`
setup and the `O(N²)` step. Not done here; recorded because it is the single
change that would most widen the range over which this validation can be run.
