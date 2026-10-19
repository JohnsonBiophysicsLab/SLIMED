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
