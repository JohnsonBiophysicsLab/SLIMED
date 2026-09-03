# Uniform Irregular Loop Patches — Implementation and Results

**Status:** implemented and merged
**Range:** `cd277a7..f813d8c`, 19 commits, +1969 / −325 lines of source
**Plans:** [`irregular_patch_valence_4_to_8_plan.md`](irregular_patch_valence_4_to_8_plan.md),
[`volume_functional_split.md`](volume_functional_split.md)

This is the record of what was built and what was measured. The two plans say
what was intended; this says what happened, including the places the plans were
wrong.

---

## 1. What exists now

One evaluation path for extraordinary valences 4–8, built as a precomputed row
table on top of the existing regular kernel. No OpenSubdiv, no per-valence
dispatch, no new subdivision scheme.

| WP | Delivered | Where |
| -- | --------- | ----- |
| 1 | One-ring guardrails: 5/6/6 predicate, no uninitialized `d4/d7/d8`, origin-referenced winding test removed, loud rejection, valence histogram | `Mesh_setup_geometry.cpp` |
| 2 | `Ā` generator for any valence, from a synthetic canonical patch | `Subdivision_matrices.{hpp,cpp}` |
| 3 | `IrregularPatchRowTable`, keyed `(valence, depth, child, sample)` | `Irregular_patch_rows.{hpp,cpp}` |
| 4 | One energy/force kernel, width from the rows | `Compute_energy_and_force_on_mesh.cpp` |
| 5 | Area and volume through the same rows | `Mesh.cpp` |
| 6 | Per-valence convergence study and chosen depths | `test_convergence_study.cpp` |
| 7 | Optional global Loop pre-refinement, off by default | `Mesh_refine.cpp` |

The central object is

```text
R[N,d,c,q] = W_q · P_c · Ā · (P4 · Ā)^d          // 7 × (N+6)
```

which has the same type and meaning as `param.shapeFunctions`, only wider or
narrower than 12. After the collapse, nothing downstream knows a face was
irregular.

**Retired:** `get_subdivision_matrices()`, `struct SubdivisionMatrix`,
`Param::subMatrix`, `Param::subDivideTimes`. The valence-5 literals live on in
`tests/subdivision_oracle.cpp` as the permanent oracle — deliberately *not*
regenerated from the generator they check.

---

## 2. What "convergence" means here

The word is used throughout this document and in the plans, and it does not
mean what it usually means in a finite-element context. Three different things
could be called convergence; only one of them is being measured.

### 2.1 The infinite ring decomposition

A regular Loop patch has a closed-form basis — the quartic box spline — so it
can be evaluated at any parameter point directly. That is what
`get_shapefunction()` returns.

An irregular patch has no such basis over the whole domain. Stam's reduction
instead decomposes it. One subdivision splits the parameter triangle into four
children; three are regular and evaluable in closed form, and the fourth still
contains the extraordinary corner. Recursing gives

```text
level 0:   3 regular children covering  3/4      of the domain
level 1:   3 more covering               3/16
level d:   3 more covering               (3/4)·4^-d
after D:   1 - 4^-D covered, 4^-D left
```

The residual is an ever-shrinking triangle at the extraordinary corner, never
covered by any finite number of levels. So an integral over an irregular patch —
area, volume, bending energy — is an **infinite sum over rings**, and evaluating
it at depth `D` truncates that sum.

**"Converging" therefore means: the truncated sum approaching the exact integral
over the limit surface as `D → ∞`.** The limit surface itself never moves. The
control mesh never changes. Depth is a property of how the integral is
evaluated, not of the geometry being integrated.

Two things that are *not* what is converging here:

- **Mesh refinement.** Adding control vertices changes the limit surface, and
  whether that surface approaches a physical membrane shape is a separate
  question this document does not address. WP7's pre-refinement does change the
  surface, which is exactly why it rebaselines a run.
- **Quadrature order.** Each regular child is integrated with a 3-point Gauss
  rule. §3.3 shows that error is negligible against truncation at these depths,
  so the two are cleanly separable.

### 2.2 Why the rate depends on the valence

The per-step map on an irregular patch's control points is `S = P4 · Ā`. Its
spectrum is what sets every rate in this document:

```text
  eigenvalue    1          the affine mode: subdivision preserves the surface
  eigenvalue    λ (×2)     the tangent-plane mode
  the rest      < λ
```

For Loop subdivision with Warren's weights, the subdominant eigenvalue is known
analytically:

```text
  λ(N) = 3/8 + cos(2π/N)/4          λ(6) = 1/2 exactly
```

Under one subdivision the residual patch shrinks by `λ` in linear extent, so its
**area** shrinks by `λ²`, while its **parameter** domain shrinks by `1/4`. The
area still missing after depth `D` therefore falls like `λ^(2D)`, and the error
should shrink by `1/λ²` per level.

Measured against that prediction:

| N | λ analytic | `1/λ²` predicted | area ratio measured |
| - | ---------- | ---------------- | ------------------- |
| 4 | 0.3750000 | 7.111 | 7.3 |
| 5 | 0.4522542 | 4.889 | 4.9 |
| 6 | 0.5000000 | 4.000 | 4.0 |
| 7 | 0.5308725 | 3.548 | 3.58 |
| 8 | 0.5517767 | 3.284 | 3.37 |

The measured ratios drift toward the predicted values with depth, approaching
from above. And the generated matrices carry the spectrum exactly: the leading
eigenvalue of `P4 · Ā` is 1 and the subdominant is a double root equal to
`λ(N)`, to 1e-10, at every valence 4–8.

This closes a loop with §3.3, where the area error at valence 6 was measured to
equal the truncated *parameter* fraction `4^-D`. That equality is specific to
valence 6, and holds precisely because `λ(6) = 1/2`, so `λ^(2D) = 4^-D` there.
At other valences the area error is `λ^(2D)` and departs from `4^-D` in both
directions.

### 2.3 Why bending converges differently, and worse

Area convergence is governed by how fast the residual patch shrinks. Bending
energy also depends on how the **curvature** behaves as the extraordinary point
is approached, and Loop surfaces are C¹ but not C² there — curvature has no
single limiting value at an extraordinary vertex. The result is a rate set by
the ratio of the subdominant eigenvalue to the next one, rather than by `λ`
alone.

Measured, per level:

| | N=4 | N=5 | N=6 | N=7 | N=8 |
| --- | --- | --- | --- | --- | --- |
| area | 7.3× | 4.9× | 4.0× | 3.6× | 3.4× |
| bending | 2.1× | 3.3× | 4.0× | 4.6× | 4.6× |

The valence-6 column comes from §3.3, which drives the table at N=6 deliberately
to have something with a known answer to compare against. Production never does
this: a valence-6 face is regular and takes the direct kernel, which is why
§3.1's sweep shows no depth dependence there at all.

At valence 6 the two rates agree exactly, which is the consistency check: a
regular vertex *is* C², curvature is bounded and continuous, and bending
inherits area's rate. Away from 6 they separate, and in opposite directions — bending is worst
where area is best.

This document does not derive the bending exponent. The depths in §3.2 were
chosen from the measured rates, not from theory, and that is the honest status:
the area law is predicted and confirmed, the bending law is measured only.

---

## 3. Benchmark results

### 3.1 Convergence is per valence, and the two quantities disagree

Measured on a fixed non-planar disk with one extraordinary corner. Factor by
which the increment shrinks per additional level:

| | N=4 | N=5 | N=6 | N=7 | N=8 |
| --- | --- | --- | --- | --- | --- |
| **area** | 7.3× | 4.9× | — | 3.6× | 3.4× |
| **bending** | 2.1× | 3.3× | — | 4.6× | 4.6× |

Area converges *fastest* at low valence; bending converges *slowest* there, by
more than a factor of two per level. Bending is the binding constraint — it is
the integrand most exposed to the truncated corner, Loop surfaces being C¹ but
not C² at an extraordinary vertex.

The plan anticipated slow bending convergence (§8) but expected it to depend on
*distance* from valence 6. It depends on *direction*: valence 4 is the expensive
case and valence 8 is cheap.

### 3.2 Chosen depths

For an estimated remaining bending tail below 1e-4 relative, measured out to
depth 14 and carried by `recommended_irregular_depth()`:

| Valence | Depth | Residual |
| ------- | ----- | -------- |
| 4 | 12 | 7.2e-5 |
| 5 | 8 | 8.3e-5 |
| 7 | 7 | 2.4e-5 |
| 8 | 7 | 2.0e-5 |

A uniform depth 6 — the previous default — would have left valence 4 at 8e-3
relative, nearly one percent.

### 3.3 The error budget is truncation alone

At valence 6 the patch is regular, so its area can be integrated directly from
the shape functions at many parameter points with no subdivision involved.
Against that independent reference:

| D | tiled area error | truncated fraction | ratio |
| - | ---------------- | ------------------ | ----- |
| 5 | 9.864e-04 | 9.766e-04 | 1.010 |
| 6 | 2.466e-04 | 2.441e-04 | 1.010 |
| 7 | 6.163e-05 | 6.104e-05 | 1.010 |
| 8 | 1.537e-05 | 1.526e-05 | 1.007 |

The relative area error **is** the truncated parameter fraction. Two things
follow: the tiling is exact — the children cover the domain once, with no double
counting and no gap beyond the sliver — and the 3-point rule contributes nothing
measurable, so choosing a depth is choosing an error.

No depth can do better than its own truncated fraction. An absolute accuracy
target below it is unsatisfiable.

### 3.4 Force–energy conjugacy

Central finite difference of energy against reported force, per control point
per axis, at every valence 4–8 and every depth 3–8:

```text
residual ≈ 1e-10        (gate: 1e-6)
h-sweep:  8.68e-09 (h=1e-8) → 1.31e-11 (h=1e-5) → 1.22e-09 (h=1e-4)
```

The V shape is the expected central-difference signature: round-off on the left,
`h²` truncation on the right.

It first measured 5e-7, flat across every valence and depth. Flat is what
correct looks like, so it was nearly banked — but it was *exactly* linear in `h`,
which no central difference can be. The cause was the measurement: reading the
analytic force inside the perturbation loop leaves it stale by one step once the
first component has moved.

### 3.5 Cost and size

| | |
| --- | --- |
| Table build, all valences, D=6 | 1.07 ms, once at startup |
| Memory | ~30 kB per level; 177 kB at D=6, 363 kB at D=12 |
| Regular face, 3 samples | 2.9 µs |
| Irregular face, D=1 (3 children × 3) | 7.8 µs |
| Irregular face, D=6 (18 children × 3) | 47.3 µs |
| Explicit recursion, same D=6 | 66.3 µs (table is 1.40× faster) |

**Correction to the plan.** §0 claims that after the collapse "an irregular face
costs the same as a regular one". The *per-sample* cost is the same; the
per-face cost is not, because an irregular face needs `9D` samples against 3 —
18× at depth 6, and 36× at depth 12. This does not threaten the design (a closed
vesicle has ~12 irregular faces, so under a millisecond per energy evaluation),
but the wording sets up a surprise.

### 3.6 Independent structural checks

Not derived from the subdivision machinery, so a shared error would not hide:

- **Rigid-motion invariance** — rotate 0.7 rad, translate far off origin: area
  and bending energy unchanged at every valence 4–8.
- **Net internal force zero** — bending and area forces on a face's own control
  points sum to zero at every valence.
- **Subdivision spectrum** — the leading eigenvalue of `P4 · Ā` is exactly 1 and
  the subdominant is a double root equal to Loop's analytic `3/8 + cos(2π/N)/4`,
  to 1e-10, at every valence. Wrong masks would still give row sums of 1 and
  planar precision, but would not land on this spectrum five times over.
- **N=5 matrix parity** — all four generated matrices bit-exact against the
  literal oracle, entry by entry.
- **Affine and planar precision** — row sums 1.0 to 1e-14 for `Ā` and every
  `P·Ā`; a planar control net stays planar at every depth and valence.

---

## 4. Backward compatibility on regular meshes

Serial build (the OpenMP build is not run-to-run deterministic — see §5). Seven
output files per run, two workloads: the flat/periodic `input.params`, and the
scaffolded one from `02de549` which pulls the membrane out of plane.

| Transition | Flat | Scaffolded |
| ---------- | ---- | ---------- |
| `cd277a7` → `54c6216` (dot_row) | 0/7 | 0/7 |
| `54c6216` → `0bb3392` (volume plan done) | — | 0/7 |
| `0bb3392` → `8f5d1c8` (NaN guards) | **4/7** | 0/7 |
| `8f5d1c8` → `b356642` (Matrix value type) | 0/7 | 0/7 |
| `b356642` → `f813d8c` (WP1–WP7) | 0/7 | 0/7 |
| **end to end** | **4/7** | **0/7** |

**The scaffolded workload is bit-identical end to end** across 19 commits that
rewrote the energy/force kernel, the geometry path, the one-ring builder and a
core linear-algebra type.

**The flat workload changed at exactly one commit**, and the change is the
intended one: 1 recorded row and 45,169 NaN tokens became 30 rows and zero.
Files written before the first step are unchanged even there.

Every comparison was re-derived from an independent second set of runs; all
seven binary/workload pairs reproduced byte-identically in fresh directories.

The reported *volume* did change with the `dot_row` fix (−32,157 → −1,218,840 on
the scaffolded run). With `uVol = 0` it never reaches the trajectory, which is
why the table shows no difference. The claim is "trajectory unchanged", not
"nothing changed".

---

## 5. Defects found that neither plan predicted

| Defect | Consequence |
| ------ | ----------- |
| `uVol / vol0` unguarded when `vol0 == 0` | Every flat run produced NaN forces and halted on iteration 0. The shipped default configuration was not running. |
| `uSurf / area0` unguarded, same asymmetry | Same cascade, reachable from `relaxArea = 0.0`. |
| `get_max_force_magnitude()` reduced NaN to 0 | A broken force field was indistinguishable from a converged one. |
| `subDivideTimes` declared without an initializer, assigned nowhere | The 11-control arm of `calculate_element_area_volume()` had an indeterminate loop bound. Any mesh with an irregular face hung. **Irregular faces were unreachable, not merely wrong.** |
| `Matrix` copy constructor dereferenced its source unconditionally | Copying a default-constructed `Matrix` was a null dereference, so `Matrix` could not live in a standard container. Cost three separate crashes. |
| `Matrix::free()` missing braces | `mat = NULL` was never inside the `if`. Harmless as written; the next edit would have inherited it. |

---

## 6. Environment notes

- **`Makefile.legacy` has no header dependency tracking.** Editing a header
  rebuilds only the `.cpp` files you also touched; the rest keep the old class
  layout. It presents as a hang or nonsense in an unrelated destructor, never as
  a link error. Always `make -f Makefile.legacy clean` after touching a header.
- **The OpenMP build is not run-to-run deterministic.** `reduction(+:)` ordering
  varies, giving ~1e-13 absolute drift in vertex trajectories between identical
  runs. Any bit-identity claim must use the `serial` target.
- `ThermalFluctuationTest.EnabledSwitchAttemptsMetropolisTrial` fails on `main`
  and predates all of this work.

---

## 7. What is not established

Stated plainly, because the test count makes it easy to over-read.

**The independent reference only exists at valence 6.** For N ∈ {4, 5, 7, 8} the
geometry is constrained by structural properties — invariance, force balance,
affine precision, N=5 matrix parity — but is not measured against an
independently computed surface. Closing that would need an external evaluator,
which the plan deliberately rules out.

**Depth 12 at valence 4 rests on an extrapolated tail.** The 7.2e-5 residual
comes from fitting a geometric ratio, not from comparison against a known
answer.

**Nothing has run irregular physics at scale.** Every irregular test is a
single-face fixture. The first closed vesicle will be the first time the row
table, the per-valence depths and the force scatter meet a full mesh under
dynamics.

**The example mesh needs none of this.** `data/example` has 3,680 faces, 336 on
the boundary, interior valences N=6 only, and zero interior faces with more than
one extraordinary corner. The irregular path is capability for closed-vesicle
topologies, not a fix to current runs.

---

## 8. Suggested next steps

1. **Land a closed vesicle workload.** An icosphere has 12 isolated valence-5
   vertices and would exercise the whole path under dynamics. It is the missing
   evidence.
2. **Re-run the depth study against that workload**, since the current numbers
   come from single-face fixtures with a fixed boundary.
3. **Consider a higher-order rule on the regular children** if bending
   convergence at valence 4 proves too slow in practice. The plan's §8 names
   this as the next move, and warns against an adaptive rule, which would make
   energy discontinuous during a timestep.
4. **Fix `ThermalFluctuationTest`**, which has been failing throughout and masks
   any future regression in that area.
