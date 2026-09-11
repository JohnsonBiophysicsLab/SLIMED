# Fluidity results

What the Monte Carlo edge-flip machinery actually does, measured. This is the
results record for work package 6 of [`edge_flip_plan.md`](edge_flip_plan.md);
that document holds the theory and the design, and this one holds the numbers.

Everything below is the CPU backend on one machine (Darwin arm64, clang++ 17,
serial). The analysis is `analysis/fluidity.py`, which reads a run directory
directly:

```bash
~/anaconda3/bin/python analysis/fluidity.py <run directory> <meshpointOutputInterval>
```

Two sheets are used throughout. The **100 nm sheet** is 525 vertices, 960
faces, 221 of them interior and free, and is the one the parameter sweeps run
on. The **60 nm sheet** is 195 vertices with 48 free, small enough to take to
60000 steps, and is the one the long-time measurements run on. Both use
`lFace = 5`, `kCurv = 83.4`, `KBT = 4.17`, `timeStep = 1e-3`, `diffConst = 1`,
`uSurf = 250` with a global area constraint, `surfaceSolver = iterative`,
`inPlaneDynamicsEnabled = true` and `edgeSpringEnabled = true`.

---

## 1. The tether had to be replaced first

WP4 shipped a harmonic mesh-quality term, `E = (k/2) (l - l0)^2`. It cannot
work for a fluid membrane, and the reason is geometric rather than a matter of
tuning. A flip on a rhombus of two equilateral triangles of side `l0` replaces
the short diagonal by the long one, so it must climb

```text
    dE = (k/2) (sqrt(3) - 1)^2 l0^2
```

whatever the rest of the Hamiltonian says. Two requirements then pull `k` in
opposite directions, and at `l0 = 5 nm` and room temperature they do not meet:

| requirement | condition | value |
| --- | --- | --- |
| flips possible (barrier under 5 kT) | `k <= 10 kT / ((sqrt(3)-1)^2 l0^2)` | `k <= 3.1 pN/nm` |
| triangulation survives (bond fluctuation under `0.1 l0`) | `k >= 100 kT / l0^2` | `k >= 16.7 pN/nm` |

Measured on the 100 nm sheet over 3000 steps:

| `edgeSpringConstant` | barrier | accepted / attempted | outcome |
| --- | --- | --- | --- |
| 83.4 (`kCurv`, the WP4 default) | 134 kT | 0 / 121 | frozen solid |
| 20.0 | 32 kT | 17 / 801 (2.1%) | stable, barely fluid |
| 1.0 | 1.6 kT | 37 / 415 (8.9%) | **diverges at step 1535** |

The `k = 1` divergence starts in the dynamics, not the flip move --
`E_curvature` reaches `2.1e7` at step 1535, before the first large flip `dE` at
step 1540 -- and it needs the flips to trigger it: the same configuration with
`edgeFlipEnabled = false` is stable at `E = 2001` after 3000 steps. A tether
too weak to prevent a degenerate triangle leaves the energy unbounded below in
a direction the flip move can reach.

`edgeTetherShape = flat` is the tether every dynamically triangulated surface
model uses: zero inside an allowed range, a quadratic wall outside. A flip that
leaves every edge inside the range costs nothing, so the wall stiffness and the
flip barrier stop being the same number and `k` can be as stiff as the walls
need. The harmonic form stays available for a minimization that never flips.

## 2. The range is narrower than it looks

Two requirements again, and this time they *just* meet.

The upper wall must exceed `sqrt(3) = 1.733`, or it forbids exactly the move it
is there to permit.

The range must also be narrow, which is the part WP6 measured and did not
expect. Inside the flat region there is no restoring force at all, so nothing
sets a length scale for the control net except these walls and the constraint
on the limit surface's area -- and a control net can be wildly non-uniform
while its limit surface stays smooth and the right size. Let the walls stand
far apart and a fluid mesh coarsens into them without ever settling.

60 nm sheet, `nu = 2`, tether energy and control-net edge length against step:

| step | `[0.6, 1.8]` tether E | edge length | `[0.95, 1.75]` tether E | edge length |
| ---: | ---: | --- | ---: | --- |
| 0 | 0 | 5.00 +- 0.00 | 0 | 5.00 +- 0.00 |
| 2000 | 597 | 5.58 +- 1.66 | 1928 | 5.23 +- 0.91 |
| 5000 | 988 | 5.60 +- 1.77 | 5410 | 5.27 +- 1.11 |
| 11000 | 9074 | 5.72 +- 1.93 | — | — |
| 16000 | 9523 | 6.03 +- 2.19 | — | — |
| 20000 | 13806 | 6.17 +- 2.26 | — | — |
| 7000 | — | — | 1681 | 5.34 +- 1.03 |

The wide range never settles: the mean edge is still growing at step 20000 and
the spread with it. The narrow one does: the edge distribution holds at
`5.3 +- 1.0` and the tether energy fluctuates about 3000 pN.nm without trend.
The same sheet with flips off sits at a tether energy of about 50 indefinitely.

So `1.84` is about as wide as the ratio may be and `sqrt(3) = 1.73` is the
floor -- barely 6% apart. **The defaults are now `[0.95, 1.75]`.** That the
window exists at all is what makes the flat tether workable where the harmonic
one is not.

The cost is a little interference with the free energy: at the narrow range the
bending energy runs about 5% higher (614 against 587 pN.nm on the 100 nm sheet
at 3000 steps) and the acceptance drops from 42% to 31%. Both are acceptable;
an unbounded drift is not.

## 3. The sweep and the dynamics were sampling different Hamiltonians

Found by this work package, and the reason the acceptance rate above is 40% and
not 17%.

The flip trial differences a local energy over eighteen faces; the Brownian
step integrates a global one. `evaluate_face_subset()` went on calling
`face_regularization_energy()` -- the *reference-length* term, which hands a
newly created edge the distance between two vertices that were never joined --
while the dynamics integrated the tether. Both halves kept working, and the
Metropolis chain sampled neither distribution.

The symptom was visible in the run and easy to misread as equilibration:

| | before | after |
| --- | --- | --- |
| acceptance | 0.172 | 0.397 |
| mean accepted `dE` | **-787 pN.nm** | **+0.3 pN.nm** (0.08 kT) |
| total `dE` over 224 flips | -95243 | +85 |

Every accepted flip reported a large energy release while the mesh's total
energy climbed. After the fix the accepted energy differences straddle zero, as
a Metropolis chain's must. `FluidHamiltonianTest.ATrialDeltaIsTheChangeInTheWholeMeshEnergy`
is the gate: it compares the number the sweep actually uses against two
whole-mesh evaluations, and against the unfixed code it reports `2e-5` where
the mesh moved by `73.2`.

## 4. Acceptance and valences against `nu`

100 nm sheet, 3000 steps. `nu` is `edgeFlipAttemptRate`, in attempts per edge
per microsecond.

| `nu` | drawn | admissible | accepted | acceptance | survival at 3000 | mean valence |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 (control) | — | — | — | — | 1.000 | 6.000 |
| 0.1 | 158 | 142 | 61 | 0.430 | 0.970 | 5.995 |
| 0.5 | 801 | 564 | 224 | 0.397 | 0.918 | 6.005 |
| 2.0 | 3245 | 2011 | 853 | 0.424 | 0.868 | 6.005 |

Three things to read off it.

**The Poisson schedule is exact.** The drawn counts are `1 : 5.07 : 20.5`
against a nominal `1 : 5 : 20`. That is what makes `nu` a physical rate rather
than a per-step count: halving the time step halves the attempts and leaves the
rate per edge alone.

**Acceptance does not depend on `nu`.** 0.43, 0.40, 0.42 -- as it must, since
acceptance is a property of the Hamiltonian and not of how often a move is
offered. It sits in the 30-50% band the DTS literature reports.

**The admissible fraction falls as the mesh disorders**: 0.90, 0.70, 0.62. More
proposals hit the valence bounds once the mesh is no longer a lattice. This is
why the accepted counts grow sub-linearly in `nu` even though the drawn counts
do not.

Valence histograms over the 221 interior free vertices, at step 3000:

| `nu` | 4 | 5 | 6 | 7 | 8 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | — | — | 1.000 | — | — |
| 0.1 | 0.023 | 0.154 | **0.665** | 0.122 | 0.036 |
| 0.5 | 0.072 | 0.158 | **0.557** | 0.118 | 0.095 |
| 2.0 | 0.068 | 0.149 | **0.561** | 0.154 | 0.068 |

Against the DTS sanity band of roughly 60% valence 6 with 20% each of 5 and 7.
The 5 and 7 fractions are lower and the 4 and 8 fractions higher than that band,
which is what `edgeFlipMinValence = 4` and `edgeFlipMaxValence = 8` permit;
narrowing those to 5 and 7 would tighten the histogram at the cost of refusing
more flips. The mean valence is 6.00 to three figures in every case, which it
must be by Euler's formula and is a check that the counting is right.

## 5. Neighbour survival

The fraction of the initial edges still present. Exactly 1 forever with flips
off, which makes it the cleanest of the three fluidity measures.

| `nu` | survival at 3000 steps |
| ---: | ---: |
| 0 | 1.000 |
| 0.1 | 0.970 |
| 0.5 | 0.918 |
| 2.0 | 0.868 |

Monotone in `nu` and still far from `1/e` at 3000 steps, so the microscopic
fluidity time is much longer than these runs: extrapolating the `nu = 2` decay
puts it near 2e4 steps, or 20 microseconds at `timeStep = 1e-3`. Calibrating
`nu` against a measured lipid neighbour-exchange time -- 0.6 to 6 microseconds
for a patch of a hundred lipids -- would therefore want `nu` of order 10 rather
than the current default of 0.5. That calibration needs runs an order of
magnitude longer than anything here and is not done.

## 6. In-plane diffusion: the gate

The classical signature. Without flips a vertex is tethered to a fixed set of
neighbours and its in-plane mean squared displacement saturates at the cage
size; with flips the cage itself rearranges and the displacement grows without
bound.

60 nm sheet, 60000 steps, `nu = 2`, time-averaged over every pair of frames a
lag apart and over the 48 free vertices, with the sheet's own drift removed:

| lag (steps) | caged (nm^2) | fluid (nm^2) | ratio |
| ---: | ---: | ---: | ---: |
| 100 | 1.129 | 1.083 | 0.96 |
| 300 | 2.485 | 2.337 | 0.94 |
| 800 | 4.215 | 4.041 | 0.96 |
| 1300 | 5.197 | 5.212 | 1.00 |
| 2100 | 5.964 | 6.491 | 1.09 |
| 3300 | 6.439 | 7.990 | 1.24 |
| 5100 | 7.012 | 9.755 | 1.39 |
| 8000 | 7.221 | 12.048 | 1.67 |
| 12400 | 7.589 | 15.121 | 1.99 |
| 19300 | 7.676 | 19.027 | 2.48 |
| 30000 | 7.727 | 23.085 | **2.99** |

Growth exponent over the upper half of the range: **0.097 caged, 0.482 fluid.**

The two curves are indistinguishable out to lag 800 -- both vertices are still
filling their cage and neither has met it -- cross at about lag 1300, and then
separate without limit. The caged curve flattens at 7.7 nm^2; the fluid one
reaches three times that and is still climbing at the end of the run.

**This is the gate, and it is met.** Two cautions on reading it. The fluid
exponent is 0.48 rather than the 1.0 of free diffusion, so the run has not
reached an asymptotic diffusive regime and no diffusion constant should be
quoted from it. And the crossover at lag 1300 is what makes short runs useless
for this: at 600 steps the ratio is 0.94, and a test that asserted a separation
there would be asserting noise. `FluidityTest.InPlaneMotionIsCagedWithoutFlips`
therefore gates only the control's saturation, which is half the signature and
the half a test can reach.

## 7. Throughput

100 nm sheet, 400 steps, `nu = 2`, trajectory output off.

| configuration | steps/s | relative |
| --- | ---: | ---: |
| flips off | 120.8 | 1.0 |
| flips on, `irregularPatchDepthScale = 1.0` | 5.7 | **21x slower** |
| flips on, `irregularPatchDepthScale = 0.5` | 9.7 | 12x slower |

The 21x is the cost of fluidity itself, not of the flip move: a fluid mesh is
mostly irregular, and an irregular face costs `3D` samples per extraordinary
corner instead of 3. WP1 measured 36x for a fully irregular mesh; 400 steps
gets part of the way there.

Halving the depth buys **1.7x** and does not disturb the Monte Carlo: the
acceptance is 0.375 at full depth and 0.369 at half, and the drawn counts are
identical. That makes `irregularPatchDepthScale` the practical lever on fluid
run cost, with the caveat that what a *spectrum* needs from the depth is not
measured -- see below.

## 8. Not done

**The fluctuation spectrum with flips on** -- gate item 4 of the plan. The
existing pipeline needs of order `1e6` steps to fit `kc` from the `q^-4` tail,
and at 5.7 steps/s a fluid run of that length is two days per trajectory on
this machine. The 60000-step runs here give 600 frames, which is short by more
than an order of magnitude. It also needs the resampling path
(`analysis/membrane_resample.py`), because with in-plane motion on, the height
field is no longer sampled on a regular lattice. Nothing here says the spectrum
is wrong; it says it has not been measured.

**Calibration of `nu` against a physical neighbour-exchange time**, for the
reason in section 5.

**The flip move on the GPU.** `DeviceMeshLayout` refuses a face with more than
one extraordinary corner, and a flip creates exactly those, so
`edgeFlipEnabled` with `forceBackend = gpu` is refused at setup. Given the 21x
in section 7, the device path is where a production fluid run should eventually
go.
