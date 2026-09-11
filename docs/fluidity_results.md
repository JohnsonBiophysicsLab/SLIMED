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

## 2. The tether range, and a measurement that was wrong

The upper wall must exceed `sqrt(3) = 1.733`, or it forbids exactly the move
it is there to permit. The lower wall is a matter of mesh quality, and the
range as a whole follows the dynamically triangulated surface literature,
whose tether ratio is 1.68-1.73. **The defaults are `[0.95, 1.75]`**, a ratio
of 1.84.

The first version of this section said something stronger: that the range
*had* to be narrow, because at `[0.6, 1.8]` the tether energy climbed without
settling and the mean control-net edge grew from 5.00 to 6.17 nm over 20 000
steps. Both numbers were computed over every edge of the sheet, ghost band
included. Split by band, they say the opposite:

| 60 nm sheet, `[0.6, 1.8]`, `nu = 2` | ghost-band tether E | interior tether E | interior edge |
| ---: | ---: | ---: | --- |
| step 10 000 | 7 683 | 9 | 5.53 +- 1.72 |
| step 30 000 | 43 421 | 1 | 5.78 +- 1.66 |
| step 40 000 | 85 361 | 9 | 5.66 +- 1.87 |
| step 60 000 | 50 047 | 29 | 5.64 +- 1.76 |

The interior tether energy never exceeds 29 pN.nm -- a few hundredths of kT
per edge -- and the interior edge distribution is stationary from the first
frame. Everything that climbed was in the ghost band, for the reason in
section 8. The wide range is fine; so is the narrow one, whose interior on the
100 nm sheet holds at 5.88 +- 1.3 nm with 0.15 kT per edge over 94 000 steps.

What the narrow range actually costs and buys: acceptance falls from 42% to
31%, the bending energy runs about 5% higher, and the thinnest triangle the
walls permit has a height of 1.85 nm instead of 0.6 nm. It is kept as the
default for the last of those, not for stationarity.

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

100 nm sheet, 300 steps from the flat start, `nu = 0.5`, trajectory output off.
Two binaries: the `Makefile.legacy` build, whose `CXXFLAGS = -std=c++14`
carries **no optimisation flag** and so is a `-O0` build, and the CMake Release
build (`-O3`, OpenMP on). The optimisation matters far more for a fluid mesh
than for a solid one, because a fluid mesh spends its time in the irregular
patch kernels, which are pure arithmetic.

| configuration | `-O0` serial | `-O3`, 1 thread | 2 | 4 | 8 threads |
| --- | ---: | ---: | ---: | ---: | ---: |
| dense solver, no flips (the path every earlier run used) | 139 | 331 | 370 | 396 | 345 |
| iterative solver, no flips | 129 | 539 | 629 | 717 | 523 |
| fluid, `irregularPatchDepthScale = 1.0` | 17.1 | 115 | 184 | 219 | **268** |
| fluid, `irregularPatchDepthScale = 0.5` | 23.0 | 147 | 218 | 259 | **307** |

steps per second. The 27 accepted flips are identical across a row: the thread
count changes nothing but the wall clock.

Read at one thread and `-O3`, **fluidity costs 4.7x** (539 against 115), and
half depth buys 1.3x of that back (147). The `-O0` column exaggerates the
ratio to 7.5x because the irregular kernels are the part `-O0` hurts most.
Threads help the fluid case far more than the solid one -- 2.3x at 8 threads
against 1.3x -- for the same reason: the per-face work is where the time goes,
and it parallelises. The solid case tops out at 4 threads, which is the number
of performance cores on this machine.

These numbers are at the flat start. A fluid mesh disorders over the first few
thousand steps and its per-step cost rises with the fraction of irregular
faces -- in the 3000-step runs of section 4 the average rate was about 2.6x
below the initial one -- so a long fluid run at 8 threads should be planned at
roughly 100 steps/s. That is **400 000 steps in about an hour**, which is what
the fluid spectrum run in `membrane_fluctuation_fluid_cpu.ipynb` uses.

**Correction to the earlier record.** The first version of this section, and
the WP5 and WP6 commit messages, quoted 21x for the cost of fluidity and "two
days per trajectory" for a spectrum run. Both were measured with the `-O0`
Makefile binary after the mesh had disordered, and the ratio was inflated by
the optimiser's absence. The 4.7x here is the number to carry.

## 8. The ghost band: why every long fluid run diverged

Found while setting up the spectrum run, and it took three wrong diagnoses to
get right.

**What happened.** The 100 nm sheet with the fluid flags, `muS = 0`, diverged
at step 8 702 at `dt = 0.002` and at step 94 078 at `dt = 0.001`: a sudden
blow-up over fifty steps, preceded in the log by a face normal reversing in
the last real row of faces. Three 25 000-step probes ruled out an OpenMP race
(the same configuration diverges at step 8 703 on one thread) and showed the
*total* tether energy climbing without settling in every fluid configuration,
area constraint or not:

| configuration | tether E at step 5 000 | at 25 000 |
| --- | ---: | ---: |
| `dt = 0.001`, `muS = 0` | 5 463 | 18 556 |
| `dt = 0.002`, `muS = 250` | 5 583 | 23 974 |

**Where it was.** Split by band over the 94 000-step run:

| step | total (as the C++ reports it) | ghost band, rings 0-2 | seam, ring 3 | interior | longest ghost edge |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10 000 | 7 762 | 7 316 | 65 | 594 | 11.1 nm |
| 50 000 | 17 446 | 16 976 | 176 | 671 | 13.1 nm |
| 90 000 | 29 531 | 29 697 | 49 | 654 | 16.9 nm |

93-98% of it is in the ghost band, and the interior is flat at about 0.15 kT
per edge. The bending energy is flat throughout as well.

**Why.** Periodicity is realised by three rings of ghost vertices and one ring
of duplicates whose *positions* are copied from the far side every step, but
whose *connectivity* is the lattice they were built with -- a flip is refused
wherever it would touch one. Once the interior has mixed, a ghost-band edge
joins two positions that stopped being neighbours long ago, and it stretches
with `1 - survival`: to 17 nm here. The tether charges for that, which is the
bookkeeping half. The other half is that the tether force on a duplicate is
mapped through `M^-T` and spread into the interior *before* the duplicate is
overwritten by its partner, so a force that should not exist is injected at
the seam every step, growing as the ghost edges grow, until a face at the seam
folds and the explicit step blows up.

**The fix.** `Mesh::edge_carries_tether()`: an edge with both endpoints ghost
or duplicate is a stale copy and carries no tether, in the force pass and in
the flip trial alike. On a closed surface, or any mesh without ghosts, every
edge still qualifies. The shipped workload is unaffected (tether off) and
stays byte-identical.

**What it does not fix.** The run with the exclusion in place diverged too,
at step 27 774 -- from a state in which the reported tether energy had been
flat at 500-670 pN.nm and the bending energy flat for the whole run, in about
two hundred steps, and once more with a face normal reversing in the last
real row of faces. The three divergences folded faces 831, 832 and 833: the
same spot every time, row 20 at column 16, where the top duplicate ring meets
the right one. A seam face there is built on three independently copied
vertices, each following a free vertex on the far side of the box, and a
fluid vertex diffuses -- 4 nm r.m.s. over this run. Once a far-side source has
drifted half an edge from its lattice site, the lattice-connected face it is
copied into folds in three dimensions and its patch energy is not finite. The
bottom and left seams are made of ghost faces, whose energy is not counted,
which is why it is always the top-right corner.

So the fluid interior is stationary and correct, and the sheet's periodic
seam is not fluid-safe: it holds for `1e4` to `1e5` steps and then folds. The
remedy is genuinely periodic connectivity -- the duplicate ring and the seam
faces following the far side's flips, or a mesh with no duplicates at all --
and it is the first item of whatever comes after this package. Until then a
fluid run on this sheet is a run of that length, and `analysis/fluidity.py`
and the notebook read it up to the divergence.

## 9. The fluctuation spectrum with flips on

Gate item 4 of the plan, measured in `analysis/membrane_fluctuation_fluid_cpu.ipynb`
against the solid `pure_a` run of `membrane_fluctuation_resample_cpu.ipynb`.
The fluid run is the same 100 nm box with the fluid flags added, `mu_S = 0`,
`dt = 0.001`, `nu = 0.5`; the trajectory analysed is the 93 500 steps written
before the divergence of section 8 -- 94 us, 749 frames after burn-in, 725
distinct connectivities -- over which the bending energy was stationary.

The reader is the subdivision route of `analysis/membrane_fluid_surface.py`,
at level 2, applied to both runs so that the comparison is between membranes
and not between readers. Its own checks: exact on the lattice (4e-8 nm, the
CSV's precision), 8e-7 nm against the C++ `surfacepoint` limit points on 41
flipped frames, and level 3 moves `kc` by 1.3% with an r.m.s. height change of
4e-3 nm.

| | frames | slope | `kc` (pN.nm) | `kc`, block mean +- s.e. | `sigma` (pN/nm) |
| --- | ---: | ---: | ---: | ---: | ---: |
| solid, exact resampler | 6401 | -3.974 | 91.51 | | +0.013 |
| solid, subdivision route | 3201 | -3.968 | 91.47 | 86.70 +- 3.38 | +0.33 +- 0.20 |
| **fluid, subdivision route** | 749 | -4.326 | 91.97 | **86.39 +- 5.00** | +1.05 +- 0.27 |

**`kc` fluid / solid = 0.996 +- 0.070**, 0.1 standard errors from 1. The fluid
membrane returns the same bending modulus as the solid one, to the 7% the run
length allows; both read a few per cent above the input 83.4, which is how
the discrete bending energy relates to the continuum one on this mesh and is
the same in both.

Two things the run is too short to settle. It covers 1.4 relaxation times of
the slowest mode in the box, and the fit from growing prefixes of it is still
moving (83.9, 87.6, 92.0 pN.nm over the last three). And the two-parameter fit
returns `sigma = 1.05 +- 0.27` against the solid's `0.33 +- 0.20`: a small
curvature at the lowest `|q|`, of exactly the kind an under-sampled slowest
mode produces, and also of the kind a real tension would. The tether is zero
inside its range and should add none; the run with the ghost-band fix
(section 8) is going to its full 400 000 steps and is what decides it.

Also measured, and not in the plan: beyond the fitting window the fluid
spectrum carries 0.60 of the solid one's power (median over 105 modes). That is
outside the window because the mesh does not resolve those modes, but it is
not the reader -- levels 2 and 3 agree there -- and it says a fluid mesh's
limit surface is smoother than a lattice's at the mesh scale. The control net
folds in projection 6% of the time and the limit net 0.6%, which is why the
reader triangulates the projected limit points afresh rather than reusing the
mesh's faces.

## 10. Not done

**The fluctuation spectrum at full length**, which waits on the seam. Once a
fluid run reaches 400 000 steps, re-executing
`analysis/membrane_fluctuation_fluid_cpu.ipynb` repeats every number of
section 9 on it; the lowest modes and the fitted `sigma` are what it settles.
Independent replicas of the present length would tighten everything but the
slowest mode.

**Calibration of `nu` against a physical neighbour-exchange time**, for the
reason in section 5.

**A fluid-safe periodic seam**, for the reason at the end of section 8. Until
then every fluid run on the sheet ends in a seam fold after `1e4`-`1e5` steps,
and the spectrum of section 9 is measured on the 94 us the longest of them
gave.

**The flip move on the GPU.** `DeviceMeshLayout` refuses a face with more than
one extraordinary corner, and a flip creates exactly those, so
`edgeFlipEnabled` with `forceBackend = gpu` is refused at setup. Given the 21x
in section 7, the device path is where a production fluid run should eventually
go.
