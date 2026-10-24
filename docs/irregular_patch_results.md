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

(Revisited in §9: this residual was measured on a bowl. A cylinder at N=8,
depth 7 leaves 3e-4, and the 1e-4 target is below the 3-point rule's own error
at extraordinary vertices.)

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

---

## 9. How many kernel evaluations an extraordinary face needs

**Measured 2026-10-16.** Fixture: [`analysis/ev_patch_cost.cpp`](../analysis/ev_patch_cost.cpp),
linked against `build-omp/libslimed_core.a` (Release). **Status: implemented as
the default on 2026-10-19 (§9.9).** Sections 9.1–9.8 are the measurements made
before the change, against the 12/8/7/7 default of that time, which they call
"production".

In a fluid run only about 10% of faces stay regular. On the final box_nu5 state
(600k steps, 960 evaluated faces) the split is 644 multi-EV, 236 single-EV and 80
regular faces. **97% of the patch-kernel samples go to the deep children of
extraordinary corners.** The cost unit used throughout this section is the
*sample*: one kernel evaluation at one quadrature point. The fixture measures
270–360 ns per sample (N=4…8), linear in the sample count to within 5%. A
centroid sample costs 3–5% more than a 3-point sample, because of per-call
overhead.

### 9.1 Method

- **Schemes** are written as `for d < D, c < 3: kernel(rows[N,d,c], w_d · gauss)`.
  For any fixed weights `w_d` this is an exact energy, so its forces are its
  exact gradient. Every option below is conservative by construction.
- **Reference.** Depth 16, with the two-rate tail correction of §9.3. Plain
  depth 14 is itself short of the infinite sum by up to 1.1e-5 (N=4 bowl) and
  4e-5 (N=8 saddle), so it is not used as the reference.
- **Fixtures.** (i) The symmetric two-ring fan of the 2026-09-15 study at
  N = 4, 5, 7, 8, under bowl, cylinder, saddle and a generic shape (a rotated
  cylinder plus cubic terms), with amplitude 0.12. The full fan of N faces is
  summed, so cross terms between angular modes cancel. (ii) A hex lattice with
  one or four flipped edges: 16 or 42 EV faces, 6 or 22 of them multi-EV, with
  valences 4–8. (iii) **The final states of the box_nu5 and box_nu0.5 runs**,
  loaded with their copy flags. The fixture reproduces the library's E_bend
  (1208.37 pN·nm) and area to 7e-16.
- **Errors** are relative to the reference. For the real states, "bending" is
  relative to the whole sheet's bending energy. "Force" is
  max |ΔF| / max |F|, taken over vertices.

### 9.2 The 3-point rule's own error is larger than the depth target

The §3.2 depths were chosen for a 1e-4 truncation tail. The quadrature error of
the 3-point rule on the same patches is larger than that. It was measured
against the 12-point degree-6 rule (`gaussQuadratureN = 6`), both at depth 14:

| | bending | area | force |
| --- | --- | --- | --- |
| regular faces (N=6 fan) | 2e-5 – 8e-5 | ≤1e-5 | 4e-4 – 9e-4 |
| EV fans, N=4…8 | 1e-5 – 8e-4 (N=8 worst) | ≤2e-5 | 1e-4 – 1e-3 |
| real box_nu5 sheet | **+3.7e-4** | −6.7e-5 | 1.6e-2 |
| production truncation, same sheet | −1.4e-4 | −7.8e-5 | 3.5e-4 |

Saddle rows are left out above: the saddle's energy is 0.3–3% of the bowl's,
and its 3-point error is ~1% of that, on regular faces as well.

The production truncation error is already below the quadrature error. Both
are far below the valence-dependent artefact of the limit surface itself, ~1%
(+16% at N=4). That figure is the converged fan bending relative to valence 6,
measured on 2026-09-15 with the same fan fixture. A truncation error anywhere
below ~3e-4 is therefore invisible next to the errors the model already
carries.

The §3.2 residual for N=8 (2.0e-5) was measured on a bowl, which never excites
the saddle mode. On a cylinder, depth 7 at N=8 leaves 3e-4.

### 9.3 Option (b), geometric tail reweighting: the single-rate version fails at N=7, 8

Multiplying the last level by `1/(1−r)` corrects one geometric mode. At N=4, 5
the bowl (k=0) mode is the slow one, so `r = r_bowl` works. At N=7, 8 the
saddle (k=2) mode is slower than the bowl mode. The bowl rate leaves the saddle
tail almost uncorrected (fan N=8 cylinder at depth 4: −2e-3). The saddle rate
over-corrects the bowl tail (+5e-3). No single rate serves both.

**Two-rate version, "tail2".** Weight the last two levels so that both
geometric tails are summed exactly: solve `w_{D−2} + w_{D−1}·r_i = 1/(1−r_i)` for
`r_bowl = 1/(16λ²)` and `r_saddle = μ₂²/λ²`, where `μ₂ = 3/8 + cos(4π/N)/4`.
Area and volume converge at the single rate λ², so they take their own weight,
`1/(1−λ²)`, on the last level only ("tail2A"). Without that weight the bending
weights over-count the area tail, by +9e-4 at N=8, depth 4.

| N | r_bowl | r_saddle | bending w_{D−2} | bending w_{D−1} | area/volume w_{D−1} |
| - | ------ | -------- | --------------- | --------------- | ------------------- |
| 4 | 1/2.250 | 1/9.000 | 0.90000 | 2.02500 | 1.16364 |
| 5 | 1/3.273 | 1/6.854 | 0.92483 | 1.68602 | 1.25712 |
| 7 | 1/4.509 | 1/2.763 | 0.83837 | 2.01378 | 1.39242 |
| 8 | 1/4.871 | 1/2.165 | 0.77828 | 2.33838 | 1.43773 |

All weights are positive. The reweighted energy is therefore a positive-weight
quadrature of the same integrand: non-negative, with an exact gradient. The
weights are constants of the valence; nothing is fitted.

On the symmetric fans, tail2A converges about 10× per level. On real geometry a
third component remains, converging at 3.3–3.6× per level. At N=4 this matches
the bowl mode combined with the in-plane k=0 distortion, 2.25 · 4λ = 3.375.
**This residual is what sets the depths below.** A third weight could remove
it, but its rate is close to the other two, so the weights would grow large and
alternate in sign.

**Implementation note.** The kernel takes one weight array for bending and area
together. With tail2A the two differ on the last two levels, so the area and
volume forces of those calls must be rescaled by `w_A/w_B` after the call. The
same weights have to reach `calculate_element_area_volume()`,
`Local_patch_energy.cpp` (the flip ΔE) and the device layout.

### 9.4 Option (c), a centroid rule on deep children: the error does not shrink with depth

Every child at depth ≥2 is a self-similar copy of the one above it. A cheaper
rule there therefore makes a *fixed relative* error on everything past the
switch level. That error plateaus instead of converging:

| switch level | N=4 | N=5 | N=7 | N=8 |
| ------------ | --- | --- | --- | --- |
| d ≥ 2 | −2.7e-3 | −2.7e-4 | −5.7e-5 | −3.2e-4 |
| d ≥ 3 | −1.1e-3 | −6.8e-5 | −3.3e-5 | −1.8e-4 |

These are fan-sum plateaus (cylinder and generic shapes). On the real sheet at
production depths, `d ≥ 3` gives −4.2e-4 bending and a force error of 2.4e-3, at
60% of the samples. tail2A beats it at equal cost by 5–30×, and adding (c) on
top of (b') only adds its plateau. **Rejected.**

### 9.5 Option (d), corner children of multi-EV faces at depth D−1

The geometric claim holds: the truncated region of a corner child at D−1 is the
same parameter region as a single-EV face's at depth D. Today a multi-EV face is
in effect one level deeper than a single-EV face. At fixed depths, (d) cuts
multi-EV samples by 13% and the whole sheet's by 11%, and multiplies the error
by 2–3 (plain production: −1.4e-4 → −3.2e-4). On the error-vs-cost front it
buys about what lowering one depth level buys. Its best showing is
tail2A 6/5/5/5 with (d): 0.56 of the cost at ≤3e-4, against tail2A 6/4/4/4
without it: 0.55 at ≤5e-4. That small gain does not pay for another code path
in the energy loop, `Local_patch_energy.cpp` and the device layout.
**Not recommended on its own.**

### 9.6 Error against cost on the real sheet

The final box_nu5 state, with box_nu0.5 in parentheses. Bending and area errors
add exactly across valences, so per-valence depth sets are assembled from
one-valence-at-a-time sweeps. The direct evaluations agree to three digits.
"Σ|e_N|" sums the per-valence errors by magnitude, a bound that does not rely on
signs cancelling. Cost is the sample count relative to production (126,171 per
force evaluation on this sheet).

| option | depths N=4/5/7/8 | cost | bending (signed) | bending Σ\|e_N\| | area | force |
| ------ | ---------------- | ---- | ---------------- | ------------------- | ---- | ----- |
| production | 12/8/7/7 | 1.00 | −1.4e-4 (−1.3e-4) | 1.4e-4 | −7.8e-5 (−7.6e-5) | 3.5e-4 (2.6e-4) |
| (a) scale 0.75 | 9/6/5/5 | 0.74 | −1.2e-3 (−1.2e-3) | 1.2e-3 | −9.2e-4 | 3.4e-3 |
| (a) scale 0.5 | 6/4/4/4 | 0.55 | −7.9e-3 (−8.0e-3) | 7.9e-3 | −3.5e-3 | 3.2e-2 |
| (b) one rate (bowl) | 6/4/4/4 | 0.55 | −1.1e-3 | 1.1e-3 | −7.8e-4 | 5.3e-3 |
| (b') tail2A | 9/6/5/5 | 0.74 | +1.8e-6 (+3.2e-5) | 5.4e-5 (3.6e-5) | −6e-9 | 2.6e-4 |
| (b') tail2A | 8/5/5/5 | 0.69 | −1.1e-5 (+1.9e-5) | 6.6e-5 (4.9e-5) | −5e-8 | 2.5e-4 (3.8e-4) |
| **(b') tail2A** | **6/5/5/5** | **0.67** | −2.7e-5 (+8.1e-6) | **8.2e-5 (6.0e-5)** | −5e-8 (+6e-7) | 2.7e-4 (3.8e-4) |
| (b') tail2A | 6/4/4/4 | 0.55 | +9.2e-5 (+2.7e-4) | 4.0e-4 (5.1e-4) | −8e-7 | 2.0e-3 |
| (b') tail2A | 5/3/3/3 | 0.42 | +4.8e-4 (+1.8e-3) | 2.6e-3 (3.8e-3) | −7e-6 | 1.2e-2 |
| (c) centroid d≥3 | 12/8/7/7 | 0.60 | −4.2e-4 | — | −8.4e-5 | 2.4e-3 |
| (d) plain | 12/8/7/7 | 0.89 | −3.2e-4 (−2.8e-4) | 3.2e-4 | −2.0e-4 | 6.6e-4 |
| (d) + tail2A | 9/6/5/5 | 0.63 | +8.1e-5 (+1.4e-4) | 1.9e-4 (1.5e-4) | −7e-8 | 9.4e-4 |
| (d) + tail2A | 6/5/5/5 | 0.56 | −2.8e-5 (+4.5e-5) | 3.0e-4 (2.4e-4) | −2e-7 | — |
| 3-point quadrature floor | — | — | +3.7e-4 | — | −6.7e-5 | 1.6e-2 |

Per-valence bending error on the box_nu5 sheet (fraction of sheet bending, other
valences at the reference):

| N | scheme | D=3 | D=4 | D=5 | D=6 | D=7 | D=8 | D=9 | D=10 | D=11 | D=12 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 4 | plain | -1e-02 | -6e-03 | -3e-03 | -1e-03 | -5e-04 | -2e-04 | -1e-04 | -5e-05 | -2e-05 | -9e-06 |
| 4 | tail2A | -1e-03 | -3e-04 | -8e-05 | -2e-05 | -4e-06 | -9e-07 | -2e-07 | -4e-08 | -1e-08 | -2e-09 |
| 5 | plain | -1e-02 | -5e-03 | -2e-03 | -5e-04 | -1e-04 | -5e-05 | -1e-05 | -4e-06 | -1e-06 | -4e-07 |
| 5 | tail2A | -1e-03 | -1e-04 | -1e-05 | -2e-06 | -4e-07 | -7e-08 | -1e-08 | -3e-09 | -6e-10 | -1e-10 |
| 7 | plain | -5e-03 | -1e-03 | -3e-04 | -9e-05 | -3e-05 | -8e-06 | -3e-06 | -9e-07 | -3e-07 | -1e-07 |
| 7 | tail2A | +1e-03 | +2e-04 | +3e-05 | +2e-06 | -3e-07 | -2e-07 | -6e-08 | -2e-08 | -5e-09 | -1e-09 |
| 8 | plain | -2e-03 | -6e-04 | -3e-04 | -1e-04 | -6e-05 | -3e-05 | -1e-05 | -6e-06 | -3e-06 | -1e-06 |
| 8 | tail2A | +9e-05 | -3e-05 | -2e-05 | -9e-06 | -3e-06 | -9e-07 | -3e-07 | -9e-08 | -3e-08 | -8e-09 |

Sample cost per level, per valence, on this sheet: N=4 1,080; N=5 5,535;
N=7 6,489; N=8 2,880. A further 3,348 samples go to regular faces and
multi-EV centre children. The production depths spend 10% / 35% / 36% / 16% of
the samples on valences 4 / 5 / 7 / 8.

The analytic fans, each relative to the fan's own energy (saddle energy is
0.3–3% of the bowl's at the same amplitude):

| option | samples/face N=4/5/7/8 | N=4 bowl / cyl / saddle | N=5 bowl / cyl / saddle | N=7 bowl / cyl / saddle | N=8 bowl / cyl / saddle | area, worst |
|---|---|---|---|---|---|---|
| production: plain 12/8/7/7 | 108/72/63/63 | -6e-05 / -6e-05 / -1e-11 | -7e-05 / -7e-05 / -1e-06 | -3e-05 / -4e-05 / -2e-03 | -2e-05 / -3e-04 / -9e-03 | 2e-04 |
| (a) plain x0.5 = 6/4/4/4 | 54/36/36/36 | -7e-03 / -7e-03 / -6e-06 | -8e-03 / -8e-03 / -2e-03 | -3e-03 / -3e-03 / -4e-02 | -2e-03 / -5e-03 / -9e-02 | 8e-03 |
| (b) tail, bowl rate, 6/4/4/4 | 54/36/36/36 | -1e-04 / -1e-04 / +3e-05 | -2e-04 / -2e-04 / +2e-03 | +3e-05 / -1e-04 / -2e-02 | +5e-05 / -2e-03 / -6e-02 | 3e-03 |
| (b) tail, slowest rate, 6/4/4/4 | 54/36/36/36 | -1e-04 / -1e-04 / +3e-05 | -2e-04 / -2e-04 / +2e-03 | +3e-03 / +3e-03 / -3e-03 | +5e-03 / +5e-03 / -5e-03 | 8e-03 |
| (b') tail2A 6/5/5/5 | 54/45/45/45 | -8e-05 / -8e-05 / +3e-06 | -7e-06 / -6e-06 / +1e-04 | -1e-05 / -1e-05 / +2e-04 | -3e-05 / -2e-05 / +2e-04 | 2e-06 |
| (b') tail2A 6/4/4/4 | 54/36/36/36 | -8e-05 / -8e-05 / +3e-06 | -7e-05 / -7e-05 / +9e-04 | -5e-05 / -4e-05 / +1e-03 | -1e-04 / -1e-04 / +8e-04 | 2e-05 |
| (c) centroid from d=3, 12/8/7/7 | 54/42/39/39 | -1e-03 / -1e-03 / -4e-04 | -1e-04 / -1e-04 / -9e-04 | -3e-05 / -8e-05 / -5e-03 | -3e-05 / -5e-04 / -1e-02 | 2e-04 |
| (c) centroid from d=2, 12/8/7/7 | 48/36/33/33 | -3e-03 / -3e-03 / -2e-03 | -3e-04 / -3e-04 / -4e-03 | -2e-05 / -1e-04 / -9e-03 | -2e-05 / -6e-04 / -2e-02 | 3e-04 |
| 12-point rule minus 3-point, depth 14 | 504/504/504/504 | +1e-06 / -3e-05 / -3e-02 | -4e-05 / +5e-05 / +2e-02 | +1e-04 / +2e-04 / +1e-02 | +2e-04 / +6e-04 / +1e-02 | 2e-05 |

### 9.7 Run time, box_nu5 parameters, 3000 steps

`~/slimed_runs/flat_tension/matrix/box_nu5/input.params` with
`maxIterations = 3000`, the worktree's `build-omp` Release binary,
`OMP_NUM_THREADS=1`, `OMP_SCHEDULE=static`. Each scale was run twice, back to
back, with no other SLIMED job on the machine. The tail weights change only the
weight constants, not the number of kernel calls, so **a tail2A run costs what
a plain run at the same depths costs**. These rows therefore time (a) and (b')
together.

| `irregularPatchDepthScale` | depths | sample cost (600k state) | wall time, rep 1 / rep 2 | speedup (best of 2) |
| -------------------------- | ------ | ------------------------ | ------------------------ | ------------------- |
| 1.0 (production) | 12/8/7/7 | 1.00 | 246.9 / 233.6 s | 1.00× |
| 0.75 | 9/6/5/5 | 0.74 | 182.6 / 180.7 s | 1.29× |
| 0.65 | 8/5/5/5 | 0.69 | 170.7 / 170.3 s | 1.37× |
| 0.5 | 6/4/4/4 | 0.55 | 144.2 / 141.5 s | 1.65× |

The times fit `t ≈ 18 s + 229 s × (sample cost)` to within 3%. The patch
kernel is ~90% of a fluid step, and all of that scales with the sample count. By
the same fit, tail2A at 6/5/5/5 (cost 0.67) would take ~171 s, a speedup of
1.37×. No flag selects that set today, so it was not timed. Option (c) at
production depths (0.60) would come to ~155 s, and (d) at production depths
(0.89) to ~222 s; both are predictions, not measurements.

These numbers stack with the factored-kernel speedup in
the 2026-10-16 fluid-cost measurements (235 → 155 s at 1 thread): that one makes each sample cheaper, this one
needs fewer samples.

### 9.8 Answer and recommendation

**How many evaluations does an EV face need?** Each extraordinary corner needs
about **45–54 samples, i.e. 5–6 levels of 3 children × 3 points, with tail
weights**, to be as accurate as production is today. Production spends 63–108
samples (7–12 levels) without tail weights. At a looser ≤5e-4, which is still
at the 3-point quadrature floor, 36 samples (4 levels) suffice. A multi-EV face
needs 3 samples for its regular centre child plus that per extraordinary
corner (3 per regular corner).

**Recommended default (now the default, §9.9):** tail2A, i.e.
two-rate bending weights plus an area/volume tail weight from §9.3, with
per-valence depths **6/5/5/5** for N = 4/5/7/8.

- **Accuracy.** No worse than production on any metric, on both real states:
  bending ≤8e-5 even when no signs cancel (production 1.4e-4), area 1e-7
  (production 8e-5), force 3e-4 (production 3e-4). On the analytic fans every
  entry is ≤2e-4. The worst entries are saddles, whose energy is 0.3–3% of the
  bowl's. Production reaches 3e-4 on the N=8 cylinder and 9e-3 on the N=8
  saddle.
- **Cost.** 0.67 of production samples. The predicted speedup is ~1.4× on
  box_nu5 fluid runs, from the measured fit above.
- **Physics.** Still an exact energy with exact gradients. All weights are
  positive and fixed per valence, so flips, conservation and the force/energy
  conjugacy test are unaffected in kind.
- **What changes.** It is a new Hamiltonian at the 1e-4 level. Like the
  factored kernel, it breaks byte-identity with existing fluid runs. It does
  not touch meshes without extraordinary vertices. The shipped periodic
  workload is one, so the legacy byte-identity gate is unaffected.

A faster setting for exploratory runs: tail2A at 6/4/4/4
(`irregularPatchDepths = 6, 4, 4, 4`). That is the cost the old
`irregularPatchDepthScale = 0.5` had (1.65× measured) at about 1/20 of its
error (≤5e-4 against 8e-3).

**Not recommended:**
- (a) alone below scale 0.75: errors of 1e-3 to 8e-3.
- (b) with one rate: wrong at N=7, 8.
- (c): its plateau does not shrink with depth.
- (d): about one depth level's worth of savings, for a third code path.

### 9.9 Implemented: the parameters and the new default

Three parameters now set the irregular-patch quadrature:

| key | default | meaning |
| --- | ------- | ------- |
| `irregularPatchDepths` | `6, 5, 5, 5` | levels for valences 4, 5, 7, 8, in that order; each in [1, 24]. `[12, 8, 7, 7]` with brackets also reads |
| `irregularPatchTailCorrection` | `true` | the §9.3 weights on the last two levels (bending) and the last level (area, volume) |
| `irregularPatchDepthScale` | `1.0` | multiplies every entry of `irregularPatchDepths`, rounded, floored at 1 |

**Reproducing a run from before the change:**

```text
irregularPatchDepths = 12, 8, 7, 7
irregularPatchTailCorrection = false
```

This pair is the old evaluation bit for bit: with both weights at 1 the
weighted kernel calls the plain kernel unchanged.
`IrregularPatchQuadratureTest.TheLegacyParametersReproduceThePlainTruncatedSum`
pins that per face. Note that a scale now multiplies the new default:
`irregularPatchDepthScale = 0.5` alone gives 3/3/3/3 with the tail correction,
not the old 6/4/4/4.

**Where the weights live.** `irregular_tail_weights()` in
`Irregular_patch_rows.cpp` holds the formula. `IrregularPatchRowTable` and
`PatchRowsFlat` carry two weights per (valence, level). The force pass,
`calculate_element_area_volume()`, `evaluate_face_subset()` (the flip ΔE) and
the device kernels (`Force_kernels.inl`, uploaded by the CUDA backend) all go
through `element_energy_force_patch_weighted_pod()` and
`element_area_volume_weighted_pod()`. The first of those rescales the area and
volume forces of a weighted level by `w_A/w_B`. The CUDA build was not
compiled or run here (no NVIDIA GPU); the host backend that shares its kernel
bodies is tested against the CPU path on valences 4–8 and on a flipped grid.

**Tests** (`tests/test_irregular_patch_quadrature.cpp`):
- the weights match the §9.3 table;
- the weights sum both geometric bending tails and the area tail exactly at
  every depth 2–24, and are positive;
- defaults and params-file parsing, including the rejected forms;
- the legacy pair is bit-identical per face;
- on fans at N = 4, 5, 7, 8 the default is within 1.5e-4 of a depth-16
  reference, and more than 10× better than the same depths without the
  correction;
- bending, area and volume forces are the gradients of their energies, to
  1e-6, with the correction on.

**Measured after the change.** `box_nu5` parameters, 3000 steps, 1 thread,
`OMP_SCHEDULE=static`, one binary for both rows. This tree includes the
factored kernel of PR #15, so the absolute times are below those in §9.7:

| quadrature | rep 1 | rep 2 | speedup |
| ---------- | ----- | ----- | ------- |
| legacy: `12, 8, 7, 7`, no tail correction | 154.3 s | 153.1 s | 1.00× |
| default: `6, 5, 5, 5`, tail correction | 112.2 s | 112.4 s | **1.37×** |

The §9.7 fit predicted 1.37×. The legacy run's final E_curvature (1068.54) is
the pre-change binary's to the printed digits.

**Byte identity, serial builds, before (5858f10) and after.**
- Shipped periodic workload, new defaults (no extraordinary faces): dynamics,
  8 of 8 output files identical; minimizer, 12 of 12.
- box_nu5 fluid run, 300 steps, with the legacy pair: 14 of 14 identical.

The full test suite passes in the serial and OpenMP builds (290 tests). Run
with the new default, a fluid workload changes at the 1e-4 level and is no
longer byte-identical to earlier runs; the legacy pair restores it.

The depth sweep behind these numbers, the scaling law and the cost of the tail
correction itself are in §10.

## 10. Run time against the recursion depth

**Measured 2026-10-24**, on the tree of §9.9 (commit 4071c4c, PR #21), to
answer two questions for users: how the step time scales with the recursion
depth, and how much of a fluid step the depth controls.

### 10.1 Method

- **Workload.** `~/slimed_runs/sticky_bead/matrix/box_nu5/input.params`
  (960 energy faces, Mixed periodic box, edge flips at rate 5, `forceBackend =
  cpu`, `surfaceSolver = iterative`, `meshpointOutputInterval = 200`,
  `localOutput = true`, `randomSeed = 20261012`), copied into one directory per
  variant with only `maxIterations` and the two lines below changed.
- **Variants.** `irregularPatchDepths = d, d, d, d` with
  `irregularPatchTailCorrection = false` for d = 1, 2, 3, 4, 5, 6, 8, 10, 12;
  the default (`6, 5, 5, 5`, tail correction on); the legacy pair
  (`12, 8, 7, 7`, tail correction off); and a frozen reference,
  `edgeFlipEnabled = false` on the regular lattice (every face regular, no
  recursion at all).
- **Per-step cost.** Flips start from a regular lattice, so the first few
  hundred steps are cheap. Every variant was run to 1000 and to 2000 steps and
  the per-step cost is `(t_2000 − t_1000) / 1000`, the cost of the disordered
  mesh. The §9 valence counts show the mesh is stationary by step 1000.
- **Timing.** The Mac could not be kept quiet (another project's MPI jobs
  held its load average at 15–70; a pass taken under that load ran 1.6–2×
  slower than a clean one and was discarded), so the sweep ran on JHU
  Rockfish: one exclusive `parallel` node (2 × Intel Xeon Gold 6248R,
  Cascade Lake, 24 cores per socket, 3.0 GHz base), Slurm job 32454969. Every
  1-thread run was pinned to its own core with `taskset`, `OMP_NUM_THREADS=1
  OMP_SCHEDULE=static OMP_PROC_BIND=true`, timed with `/usr/bin/time -p`
  (`user` time), eight runs at once on a 48-core node (four streams per socket),
  and the whole 24-run series was run twice (rep 1 on socket 0, rep 2 on
  socket 1). Binary: CMake Release (`-O3 -DNDEBUG`), `SLIMED_ENABLE_OPENMP=ON`,
  GCC 12.3.0, GSL 2.7, sha256 prefix `d67c6d84f74f075e`, built from this tree
  (`rsync` of the worktree, commit 4071c4c). The 4-thread rows are wall time
  on four pinned cores. The job script, stream runner and run tree are in
  `analysis/depth_cost_sweep/`.
- **Kernel calls per step.** The final connectivity of each run was rebuilt by
  replaying the accepted flips of `inputEdgeFlips.csv` on `inputmesh_faces.csv`
  (`analysis/depth_cost_replay.py`), and the force pass's patch-kernel calls
  counted as 1 for a regular face, `3·D(N)` for a face with one extraordinary
  corner of valence N, and for a face with several, the four Loop children:
  1 for the centre, 1 for a regular corner child, `3·D(N)` for an extraordinary
  one (`Mesh::compute_face_energies_and_forces`,
  `include/mesh/Multi_extraordinary_patch.hpp`). The area pass and the flip
  trials call the same rows again; they are not in this count, so the
  per-call cost below is an upper bound on the force kernel alone.

### 10.2 Results: the step is linear in the depth

Per-step cost of steps 1000–2000, `user` seconds / 1000, one thread, one core.
`K/step` is the force pass's patch-kernel calls per step on the final mesh
(mean of the 1000- and 2000-step states), `EV corners` the extraordinary
corners of the 960 energy faces (same mean), flips the accepted flips by step
1000 / 2000.

| depths (valence 4/5/7/8) | tail | ms/step rep 1 | ms/step rep 2 | steps/s | ×frozen | K/step | EV corners | flips |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| frozen reference (regular lattice, no flips) | – | 4.07 | 4.08 | 246 | 1.0 | 960 | 0 | 0 / 0 |
| 1/1/1/1 | off | 14.72 | 14.72 | 67.9 | 3.6 | 3 303 | 777 | 183 / 638 |
| 2/2/2/2 | off | 39.58 | 39.56 | 25.3 | 9.7 | 11 277 | 1 704 | 1 015 / 2 255 |
| 3/3/3/3 | off | 50.09 | 50.14 | 20.0 | 12.3 | 16 455 | 1 710 | 1 324 / 2 538 |
| 4/4/4/4 | off | 61.38 | 61.41 | 16.3 | 15.1 | 20 154 | 1 596 | 1 350 / 2 643 |
| 5/5/5/5 | off | 70.10 | 74.92 | 14.3 | 17.2 | 25 533 | 1 632 | 1 357 / 2 619 |
| 6/6/6/6 | off | 80.03 | 80.13 | 12.5 | 19.7 | 30 738 | 1 650 | 1 336 / 2 565 |
| 8/8/8/8 | off | 101.17 | 101.17 | 9.9 | 24.9 | 40 377 | 1 638 | 1 308 / 2 546 |
| 10/10/10/10 | off | 121.57 | 121.60 | 8.2 | 29.9 | 49 932 | 1 629 | 1 308 / 2 548 |
| 12/12/12/12 | off | 142.45 | 142.12 | 7.0 | 35.0 | 59 706 | 1 629 | 1 308 / 2 548 |
| **6/5/5/5 (default)** | **on** | **79.72** | **79.71** | **12.5** | **19.6** | 25 803 | 1 629 | 1 308 / 2 548 |
| 12/8/7/7 (legacy) | off | 96.63 | 96.71 | 10.3 | 23.7 | 38 666 | 1 638 | 1 308 / 2 546 |

The two repetitions agree to 0.1% except 5/5/5/5 at 2000 steps in rep 2
(+3.5%, one run). The 1000-step totals (not shown; `results/`) agree to 0.1%
in every row.

**The flips are not identical across depths.** The flip ΔE is evaluated with
the same quadrature, so a different depth is a different Metropolis trial
energy: at depth 1 the acceptance collapses to 3–6% (183 / 638 flips) and at
depth 2 it is 21–25% (1 015 / 2 255) against 29–31% from depth 3 up; from
depth 8 up the flip history is identical (1 308 / 2 546–2 548). The rows from
depth 3 up therefore time the same kind of mesh (1 600–1 710 extraordinary
corners, 90% of faces extraordinary); the depth-2 mesh happens to carry the
same corner count by step 1000; the depth-1 mesh has half of it (777 corners,
59% of faces still regular), so its time is not comparable and it is left
out of the fits below. Depth 1 and 2 are also outside the accuracy the model
needs (§9.6), so they are only here to anchor the line.

**Fit.** Over the uniform depths 3–12, `ms/step = a + b · depth`:

| | a (ms/step) | b (ms/step per level) | R² | rms residual | largest residual |
| --- | ---: | ---: | ---: | ---: | ---: |
| rep 1 | 19.45 | 10.22 | 0.99971 | 0.53 ms | 1.7% |
| rep 2 | 21.28 | 10.05 | 0.99762 | 1.49 ms | 4.5% (the 5/5/5/5 outlier) |
| rep 1, depths 2–12 | 19.31 | 10.24 | 0.99977 | 0.51 ms | 1.8% |

**Scaling law: each recursion level adds 10.2 ms to a 960-face step, on top
of a 19–21 ms step at depth → 0.** Linear, no measurable curvature from depth
3 to 12 (the residuals are at the run-to-run noise). The slope is what §9's
cost model predicts: one level is 3 regular children per extraordinary corner,
3 × 1 641 corners = 4 920 extra kernel calls per step, so **one kernel call
costs 2.08 µs (rep 1) / 2.04 µs (rep 2)** on this core. That figure is an
upper bound on the force kernel alone: the area pass
(`calculate_element_area_volume`) walks the same rows with the cheaper
area/volume kernel, and every flip trial re-evaluates its patch, and both
scale with the depth too. Fitting the time against `K/step` instead of the
depth gives 2.1 µs per call and an intercept of 16.6 ms (R² = 0.9987); the
legacy pair lands on that line to 1.3% (predicted 97.9, measured 96.6 ms).
The default landed 11% above it (predicted 70.9, measured 79.7 ms) until the
weighted wrapper was fixed; see §10.4.

For scale, the depth → 0 intercept is itself 4.8× the frozen step: it holds
the first level's bookkeeping that does not scale with depth (the control-net
prolongations of the 540–630 multi-extraordinary faces, the row look-ups, the
flip sweep, the area pass's regular faces) plus everything a frozen step does.

### 10.3 How much of the step the depth controls

Two references: the depth → 0 intercept *a* of the fit (what a fluid step
would cost if the extraordinary faces were free; the honest reference for
"what the depth parameter controls") and the frozen reference (what the same
sheet costs with every face regular; the reference for "what the edge flips
cost in total"). Share = (t − t_ref) / t.

| quadrature | ms/step | share controlled by the depth (vs a = 19.5 ms) | share above frozen (4.07 ms) |
| --- | ---: | ---: | ---: |
| 2/2/2/2 | 39.6 | 51% | 90% |
| 3/3/3/3 | 50.1 | 61% | 92% |
| **6/5/5/5, tail correction (default)** | **79.7** | **76%** | **95%** |
| 12/8/7/7 (legacy) | 96.6 | 80% | 96% |
| 12/12/12/12 | 142.4 | 86% | 97% |

At depth 1 the measured step (14.7 ms) sits below *a* because its mesh is only
half as disordered; on the stationary mesh the fit puts a depth-1 step at
29.7 ms (34% depth-controlled) and a depth-2 step at 39.9 ms (51%).

**Answer: the recursion depth is the dominant cost of a fluid step.** At the
PR default it controls three quarters of the step, at the legacy depths four
fifths, and the whole extraordinary-face machinery (depth plus the fixed part)
is 95–96% of a fluid step, which is 20–24× a frozen step on the same sheet.
Halving the default depths (3/3/3/3, outside the §9 accuracy target) would
buy 1.6×; nothing else in the step is worth more than the remaining 20 ms.
Going the other way, the legacy 12/8/7/7 costs 1.21× the default here for no
accuracy gain (§9.6); §9.9 measured 1.37× on the Mac, see below.

**Four threads (wall time, four pinned cores, same node):**

| quadrature | 1 thread ms/step | 4 threads ms/step | speedup | share above frozen, 4 threads |
| --- | ---: | ---: | ---: | ---: |
| frozen | 4.07 | 2.70 | 1.51× | – |
| 6/5/5/5 default | 79.7 | 28.3 | 2.82× | 90% |
| 12/8/7/7 legacy | 96.6 | 34.6 | 2.80× | 92% |

The legacy/default ratio is 1.21× on one thread and 1.22× on four: the share is
a property of the step, not of the thread count. (The frozen step
parallelises worse because its 4 ms is mostly serial bookkeeping.)

**Against the Mac.** A clean partial pass on the Apple M5 (load < 4.4, same
binary target, Apple clang 21, depths 1–3 only before another project's jobs
took the machine) gave 8.4 / 21.1 / 26.9 ms per step at depths 1 / 2 / 3
against 14.7 / 39.6 / 50.1 here: the M5 core is 1.85× faster. Its two
comparable points (depths 2, 3) give a slope of 5.8 ms per level and an
intercept of 9.5 ms, which put the default at 39 and the legacy pair at 54 ms
per step (ratio 1.39; §9.9 measured 1.37 over 3000 steps) with the same
shares, 76% and 82%. The legacy/default ratio is 1.21 here and 1.37 on the
Mac because the fixed part of the step weighs more on the Xeon; the
depth-controlled share is the same on both.

### 10.4 The tail correction cost about one level under GCC 12; fixed by folding the weights into the kernel

The default's 11% excess over the kernel-call line was isolated on a second
node of the same type (job 32455384, c088), one thread, same protocol:

| depths | tail correction | ms/step | flips (1000 / 2000) |
| --- | --- | ---: | --- |
| 6/5/5/5 | off | 69.74 | 1 349 / 2 622 |
| 6/5/5/5 (default) | **on** | 79.91 | 1 308 / 2 548 |
| 5/5/5/5 | off | 69.87 | 1 357 / 2 619 |
| 5/5/5/5 | **on** | 79.08 | 1 308 / 2 548 |

(The default's 79.91 ms here against 79.72 / 79.71 ms on c453 is the
node-to-node spread, 0.2%.) **Switching the tail correction on costs 9.2–10.2 ms
per step at fixed depths, 13% of the default step, about the price of one more
recursion level**, although it changes no kernel-call count: it only routes
the last two levels of every extraordinary corner (2 × 3 × 1 641 ≈ 9 800 calls
per step) through `element_energy_force_patch_weighted_pod()`, which scales
the Gauss weights, zeroes and rescales two per-child force buffers and divides
the mean curvature (`include/energy_force/Patch_kernel.inl`), and the matching
area wrapper. That is ~0.95 µs per weighted call, nearly half a kernel call,
far more than those loops should cost. **This penalty is specific to the
GCC 12.3 build used here.** On the Mac with AppleClang 21 (same tree, same
workload, 2026-10-24, session notes) an instrumented build counts the same
kernel calls per step with the correction on and off (48 965 against 48 403
force-kernel calls per step, 38% of them through the weighted arm when on), a
micro-benchmark of the two arms gives 0.53 against 0.51–0.56 µs per call, and
a `sample` profile splits the kernel's time between the two groups of call
sites in exactly the plain/weighted call ratio. So under clang the weighted
path costs what the plain path costs and §9.7's statement holds there; under
GCC 12 it is wrong by 13%, most plausibly the inliner treating the two arms
of the wrapper differently or the per-child buffers, not yet confirmed. A
weight-folding kernel (the weights multiplied into the kernel's own `0.5·g`
factors, the plain kernel a forwarder with weights 1, bit-identical when off)
removes the question for every compiler, and is what this branch does.

**Fix and measurement.** `element_energy_force_patch_weighted_pod()` is now
the kernel itself: each sample's `0.5·g` factor is multiplied by `bendWeight`
for the bending terms and by `areaWeight` for the area and volume terms (mean
curvature and the normal stay unweighted), and `element_energy_force_patch_pod()`
forwards to it with both weights 1; the area kernel is folded the same way. A
multiplication by exactly 1 is exact, so the uncorrected path computes bit for
bit what it did: the legacy pair `12/8/7/7`, correction off, writes all 19
output files of a 500-step box_nu5 run byte-identical on the Mac (AppleClang)
and `inputvertex_final.csv` byte-identical after 2000 steps on Rockfish (GCC
12.3). With the correction on, the bending energy is bit-identical to the old
wrapper's and the area/volume forces and mean curvature differ by rounding
only: 500 steps of the default give the same flip history and a maximum
coordinate difference of 6.5e-13 nm. Timed on Rockfish with the §10.1
protocol (job 32456577, node c504, same node type, one pinned thread, `user`
seconds of steps 1000–2000), old binary (commit 4071c4c) against this branch:

| depths | tail correction | before, ms/step | after, ms/step | flips (1000 / 2000) |
| --- | --- | ---: | ---: | --- |
| 6/5/5/5 | off | 68.8 | 69.7 | 1 349 / 2 622 |
| 6/5/5/5 (default) | **on** | **79.1** | **70.5** | 1 308 / 2 548 |

The flip histories are identical between the two binaries in both rows. **The
default step drops from 79.1 to 70.5 ms (1.12×), and the correction now costs
1.2% at fixed depths** (70.5 against 69.7 ms) instead of 15%; the 0.9 ms of
the uncorrected row is at the level of the run-to-run spread seen in §10.2 and
of the four extra multiplications per quadrature sample, and is not resolved
here. The 13% was therefore the wrapper (two zeroed and rescaled per-child
buffers and a second kernel body for the weighted arm), not the weights.
Correcting the §10.3 shares with the new default: depth controls 76% of the
step as before, since both the intercept and the kernel-call line are unchanged;
the default now lands 0.5% below the kernel-call prediction of 70.9 ms instead
of 11% above it. Note also that the instrumented count
of ~49 000 force-kernel calls per step is 1.9× the face-pass count used in
§10.2 (the flip trials and the post-flip re-evaluation call the kernel too),
so the 2.05 µs per call there is nearer 1.1 µs per actual kernel call. The weights also change the flip history (1 308 against
1 349–1 357 accepted flips by step 1000), as any change of the trial energy
does.

### 10.5 Answers

1. **Scaling law.** The step time is linear in the recursion depth,
   `t = (19.5 ± 1) ms + 10.2 ms × depth` per step for this 960-face sheet on one
   Xeon 6248R core (R² > 0.997, depths 3–12; the Apple M5 core is 1.85× faster
   with the same shape). Each level is 3 kernel calls per extraordinary
   corner, 4 900 calls per step on a mesh with 1 640 of them, at 2.05 µs per
   call (force pass plus its share of the area pass and flip trials).
2. **Share.** The depth controls 76% of the step at the default 6/5/5/5 (80% at
   the legacy 12/8/7/7, 51% at a uniform depth 2, 86% at 12), and the
   extraordinary-face evaluation as a whole is 95–96% of a fluid step, which
   is 20–24× the frozen step of the same sheet. It is by far the largest lever
   in the fluid step. The second-largest, the tail correction's own 13% on
   the GCC build (not present under AppleClang), is gone since the weights
   were folded into the kernel (§10.4): the default step is 70.5 ms on the
   Xeon, 1.2% above the same depths uncorrected.

**Commands.** On the worktree of commit 4071c4c:

```text
cmake -S . -B build-omp -DSLIMED_ENABLE_OPENMP=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-omp --target membrane_dynamics -j 4
# per run directory (input.params = box_nu5 + maxIterations + the two depth lines):
OMP_NUM_THREADS=1 OMP_SCHEDULE=static /usr/bin/time -p build-omp/bin/membrane_dynamics
# cluster: sbatch analysis/depth_cost_sweep/depth_sweep.sbatch   (then collect.sh, depth_cost_fit.py)
```
