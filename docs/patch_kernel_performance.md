# Patch kernel performance: the factored force loop

A fluid (edge-flip) run cost about 23x a frozen one on the flat matrix box:
235 s against 10 s for 3000 steps, 960 faces, one thread, CMake Release
`-O3`. The flip sweep is a few per cent of that. About 90% is
`slimed::element_energy_force_patch_pod()` in
`include/energy_force/Patch_kernel.inl`. Once flips disorder the valences
(steady state: about 25% of vertices at valence 5, 20% at 7, 6% at 4 and 7% at
8), only about 10% of faces are still regular. A face touching one
extraordinary vertex is evaluated as 3 x depth regular sub-patches by Stam's
recursion (depth 8, 7, 7 and 12 for valence 5, 7, 8 and 4). The mean is
about 40 kernel calls per face, against 1 on a frozen mesh. Making each call
cheaper is the lever. This change does that without touching the physics.

## 1. What changed

### The per-control-point loop

For each quadrature sample the kernel builds the resultants `n1_be, n2_be`
(bending stress), `m1_be, m2_be` (bending couple), `n1_cons, n2_cons` (area)
and `n1_conv, n2_conv` (volume). It then loops over the patch's control
points. For point `j`, with shape-function values `sf0..sf6` (rows: `x`, `a_1`,
`a_2`, `a_11`, `a_22`, `a_12`, `a_21`), the old loop built two 3x3
outer-product sums:

    da1 = -(sf3 a1 (x) a_3 + sf1 a11 (x) a_3 + sf1 a1 (x) a_31
          + sf6 a2 (x) a_3 + sf2 a21 (x) a_3 + sf2 a2 (x) a_31)
    da2 = -(sf5 a1 (x) a_3 + sf1 a12 (x) a_3 + sf1 a1 (x) a_32
          + sf4 a2 (x) a_3 + sf2 a22 (x) a_3 + sf2 a2 (x) a_32)

It contracted them with `m1_be` and `m2_be`. That is about 400 flops per point
per sample, and the arithmetic is the same for every point except for the
seven scalars `sf_k`.

Contracting `m` into an outer product `u (x) v` gives `(m . u) v`, so the couple
collapses onto eight sample-only dot products:

    m1a1 = m1_be.a1    m1a2 = m1_be.a2    m2a1 = m2_be.a1    m2a2 = m2_be.a2
    m1a11 = m1_be.a11  m1a21 = m1_be.a21  m2a12 = m2_be.a12  m2a22 = m2_be.a22

The force on point `j` is then linear in its shape-function values,
`f_j = sum_k sf_k(j) G_k`. With `c = 0.5 w_s * -sqa`, the sample-only vectors
are:

    gb1 = c (n1_be - (m1a11 + m2a12) a_3 - m1a1 a_31 - m2a1 a_32)
    gb2 = c (n2_be - (m1a21 + m2a22) a_3 - m1a2 a_31 - m2a2 a_32)
    gb3 = -c m1a1 a_3    gb4 = -c m2a2 a_3    gb5 = -c m2a1 a_3    gb6 = -c m1a2 a_3
    ga1 = c n1_cons      ga2 = c n2_cons
    gv0 = c volScale a_3 gv1 = c n1_conv      gv2 = c n2_conv

    fBend_j += sum_{k=1..6} sf_k gb_k
    fArea_j += sf1 ga1 + sf2 ga2
    fVol_j  += sf0 gv0 + sf1 gv1 + sf2 gv2

That is 33 multiply-adds per point (three components of 6 + 2 + 3 terms),
plus a fixed cost per sample to build the eleven `G` vectors. The terms are
exactly those of the expanded form, only reassociated. Energies, mean
curvature, normals, and the area/volume integration are untouched.

### The OpenMP face loops

`Mesh::compute_face_energies_and_forces()` and
`Mesh::calculate_element_area_volume()` now use `schedule(static, 4)` instead of
the default contiguous static blocks. Face cost varies 40-fold between a
regular face and one next to a valence-4 vertex, and on a fluid mesh the
expensive faces cluster, so contiguous blocks leave threads waiting. Small
interleaved chunks balance the load. The schedule is still static, so for a
given thread count the face-to-thread map, and with it the order in which the
per-thread force buffers are summed, is fixed. A dynamic schedule measured
the same speed but is not reproducible, so it is not used.

The area/volume loop writes only per-face fields, so its schedule cannot change
any result. The force loop's schedule changes which thread buffer a face lands
in, so an OpenMP run is not bit-identical to an OpenMP run of the previous tree
at the same thread count. A serial build ignores the pragma entirely.

## 2. Verification

**Rounding, not physics.** The legacy byte-identity gate
(Periodic/Free/Fixed, `sideX = sideY = 100`, harmonic bond off,
`relaxArea = 5456`, serial Release builds of the parent commit and of this
change) can no longer be byte-identical, because the reassociation changes the
last bits of the force. The decision on 2026-10-16 was to **re-baseline**: the
force loop is one body for every boundary mode and for the CUDA backend, and
binaries built from this commit are the reference for future byte-identity
checks. The measured differences against the parent commit:

| workload | identical files | differing files | max abs difference |
| --- | --- | --- | --- |
| legacy dynamics, Periodic, 300 steps | `EnergyForce.csv`, `inputface.csv`, `inputvertex_begin.csv`, `inputvertex_type_begin.csv`, `meshpointinput.csv`, `surfacepointinput.csv` | `inputvertex_final.csv`, `inputvertex_type_final.csv` (485 of 1575 coordinates) | 1.0e-13 nm (relative 1e-14) |
| legacy minimizer, Periodic, 25 iterations | `ElementFaceEnergy.csv`, `EnergyForce.csv`, `face.csv`, `vertex_begin.csv`, `vertex_initial_step_0.csv` | `vertex*.csv` (33 of 1575 numbers), `slimed_restart.chk` | 4e-21 nm on coordinates that are zero to the precision of the sheet (2.6e-18 in the checkpoint) |
| fluid box (`box_nu5`, Mixed periodic, `nu = 5`), 3000 steps | `inputEdgeFlips.csv`, `EnergyForce.csv`, all `inputface*.csv`, `meshpointinput.csv`, `surfacepointinput.csv` | `inputvertex_final.csv`, `inputvertex_type_final.csv`, `inputbox.csv` | 1.3e-12 nm (relative 1.8e-12); box 1e-12 |

The flip history of the fluid run is identical step for step: the drift is far
below anything a Metropolis test resolves.

**Tests.** All 282 unit tests pass in both the serial and the OpenMP CMake
builds (283 run; the CUDA backend test skips without a device).

**CUDA.** The kernel is shared with `src/cuda/Cuda_force_backend.cu`. The new
code uses only `SLIMED_HD` helpers (`v3_dot`) and stack arrays, which is
already true of the rest of the kernel. There is no `nvcc` on the development
Mac, so the `.cu` was compiled as host C++ against a scratch `cuda_runtime.h`
stand-in, with each `kernel<<<g, b>>>(args)` rewritten into a serial loop. It
was linked into the test binary in place of the stub, and
`CudaForceBackendTest.MatchesTheProductionForceEvaluation` passes. A real
device build still needs the cluster. The kernel's per-thread stack is
slightly smaller than before (eleven 3-vectors once per sample in place of
two 3x3 matrices per point).

## 3. Timing

Flat matrix box (`box_nu5` parameters: 960 faces, Mixed periodic, elastic box
on, 3000 steps), CMake Release `-O3` with OpenMP, Apple M5 (4 performance + 6
efficiency cores). Runs were back to back, with no other SLIMED job running.

| run | threads | before | after | speedup |
| --- | ---: | ---: | ---: | ---: |
| frozen (`edgeFlipEnabled = false`) | 1 | 10.5 s | 9.4 s | 1.12x |
| fluid, `nu = 5` | 1 | 235.2 s wall, 234.0 s CPU | 154.9 s CPU | 1.51x |
| fluid, `nu = 5`, prototype measurement | 1 | 235 s | 155 s | 1.52x |
| fluid, `nu = 5`, prototype measurement | 4 | 99 s (default contiguous schedule) | 55 s (`static, 4`) | 1.8x |

Factored kernel alone at 4 threads, prototype: 66 s with the default schedule
and 55 s with `static, 4`. The interleaved schedule is therefore worth about
1.2x on top of the kernel at 4 threads, and nothing at 1 thread.

**How far to trust each row.** The machine was on battery in Low Power Mode
and slept between tool calls during the 2026-10-16 measurement. The frozen
rows and the old fluid run completed while it was awake (wall about equal to
CPU). The new fluid run's wall clock (776 s) includes a sleep, but its CPU time
(154.9 s, single thread, so CPU equals work) agrees with the prototype's
155 s wall to 0.1%. The 4-thread runs of this tree were dominated by sleep and
are not quoted. The 4-thread row is the prototype session's measurement on
the same box. It has the same kernel and schedule as this commit, built
behind `#ifdef SLIMED_FACTORED_PATCH_FORCE` and with `OMP_SCHEDULE=static,4`.
Re-time it on mains power before quoting it elsewhere.

The frozen run gains little because a regular face is one kernel call and the
dense per-step work around it is unchanged. The fluid run is where the
irregular sub-patches multiply the kernel's cost about 40-fold, and the saving
scales with them.

All four completed fluid runs (old and new kernel at 1 and 4 threads) produce
the same flip history and the same `EnergyForce.csv`. Their final coordinates
agree to 1.8e-12 nm, so neither the kernel nor the thread count moved the
trajectory beyond rounding.

## 4. Further levers, not taken here

Both of these change the physics or the sampled ensemble, so they are run
parameters rather than part of this change:

- `edgeFlipMinValence = 5, edgeFlipMaxValence = 7` removes the valence-4
  vertices, which are both the most expensive (depth 12) and the ones with the
  largest valence-dependent bending artefact (measured at about +16%, against
  about 1% at valence 5-7). Measured on the prototype: 155 s to 107 s at one
  thread.
- `irregularPatchDepthScale = 0.5` halves the recursion depth: 155 s to 98 s,
  at an estimated ~1e-2 relative bending tail per extraordinary fan.
