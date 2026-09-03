# The Volume Functional Split

**Status:** implementation plan, not yet implemented
**Base:** `JohnsonBiophysicsLab/SLIMED @ 719025a`
**Companion:** [`irregular_patch_valence_4_to_8_plan.md`](irregular_patch_valence_4_to_8_plan.md)
**Sequencing:** do this first, before any irregular-patch work.

## 0. The finding

The question was whether the volume-functional split documented in `SLIMED_dev`
also exists upstream. It does, with the same root cause, but without the
complications that grew around it downstream: this tree has one call site, no
CUDA duplicate, no per-valence copies, and no compatibility mode to preserve.
The fix is roughly ten lines.

The split is not two competing scientific conventions. It is a copy-paste bug
in a linear algebra helper.

### Root cause

`src/linear_algebra/Linear_algebra.cpp:339`

```cpp
double dot_row(const Matrix &m1, const Matrix &m2)
{
    double result = 0.0;
    for (int i = 0; i < m1.mat->size1; i++)          // size1 is the ROW count
        result += gsl_matrix_get(m1.mat, 0, i) *     // but i indexes COLUMNS
                  gsl_matrix_get(m2.mat, 0, i);
    return result;
}
```

For a `1 x 3` row vector `size1 == 1`, so the loop runs once and the function
returns `m1(0,0) * m2(0,0)` — the x-component product alone. Eight lines above,
`dot_col` uses the same bound while indexing `(i,0)`, which is correct for the
`3 x 1` columns it is called with. The bound was copied; the indices were not.

### Sole call site

`src/mesh/Mesh.cpp:204-206`, inside
`Mesh::enumerate_gauss_quadrature_point_area_volume()`

```cpp
// v = 1/3 * s * dot(x, d)  <<< tetrahedron volume
// namely -> double v = dot_row(x, a_3) / 3.0;
volume += 0.16666666666 * coeff * dot_row(x, a_3);
```

This is the only `dot_row` call in the codebase, and it is the accumulator
behind every `face.elementVolume` and `param.vol`. The comment above it states
the intended functional correctly. The line below it does not compute that
functional.

## 1. What the two functionals are

Both standard divergence-theorem forms are valid for a closed surface, but they
carry different constants. The code has taken the constant from one and the
integrand from the other.

```text
intended  — full divergence, F = x
   V = (1/3) * closed_integral( x . n  dA )
     = (1/6) * sum_q  w_q * dot( x_q , a_1 x a_2 )

also valid — single axis, F = (x, 0, 0)
   V = closed_integral( x_1 * n_1  dA )
     = (1/2) * sum_q  w_q * x_q,1 * (a_1 x a_2)_1

actual — the x-only integrand with the full-divergence constant
     (1/6) * sum_q  w_q * x_q,1 * (a_1 x a_2)_1   ==   V / 3
```

The `1/6` is correct for the full form: it combines the `1/3` from the
divergence theorem with the `1/2` of the reference triangle under the
weights-sum-to-one convention. Paired with the single-axis integrand, which
needs `1/2`, it lands a factor of three low.

### Verified numerically

On a 5,120-face icosphere of unit radius, evaluating both accumulators over
flat triangles:

| Accumulator            | Value         | vs. 4pi/3 |
| ---------------------- | ------------: | --------: |
| full (intended)        |  4.179738948  |   -0.22%  |
| x-only with 1/6 (actual) |  1.393246316  |  -66.74%  |
| **ratio**              | **3.0000000000** |        |

The -0.22% on the full form is the expected polyhedral-inscription deficit, not
error. The ratio is exact.

## 2. Why it matters: energy and force stop being conjugate

The volume *force* is correct. It is the full analytic derivative, built from
three-component `dot_col` calls, and its coefficient is consistent with the
energy.

`src/energy_force/Compute_energy_and_force_on_mesh.cpp:211`, `:266`, `:484-493`

```cpp
// energy
sumEnergy.energyVolume = 0.5 * param.uVol / param.vol0
                       * pow(param.vol - param.vol0, 2.0);
// force, with uVol = param.uVol / param.vol0 rebound at :266
tmp_evol = uVol * (vol - vol0) / 3.0;
const_multiplication(a1, dot_col(x, a_3), tmp_f);   // full 3-vector
```

Nothing is wrong on this side. The coefficients line up:
`dE/dvol = (uVol/vol0)(vol - vol0)` is exactly what the force uses.

So the two halves are individually defensible and jointly inconsistent: **the
force differentiates the full-divergence volume, while the `vol` fed into both
the energy and the force is the x-only one.** Consequences:

- **The reported force is not `-dE/dx`.** A finite-difference check of total
  energy against total force fails on any mesh where the volume constraint is
  active.
- **Energy minimization is internally inconsistent.** The line search evaluates
  one functional; the descent direction comes from another.
- **Brownian dynamics inherits the same defect**, since the displacement
  equation is driven by that force.
- **The target is on the wrong scale.** `vol0` is compared against `V/3`, so a
  target set from a physical volume is off by 3x.

### Triage by run type

| Configuration                        | Effect |
| ------------------------------------ | ------ |
| `uVol == 0` (no volume constraint)   | Trajectories unaffected. The volume written to output is `V/3` on a closed mesh, and not a volume at all on an open one. |
| `uVol != 0`, closed surface          | Trajectories wrong. Constraint targets the wrong scale and the force does not match the energy gradient. |
| `uVol != 0`, flat/periodic workload  | Undefined. Neither functional is a volume on an open surface, and the x-only form is not even origin-independent there. |

## 3. The fix

Deliberately small. The temptation is to build a functional-selection layer;
resist it. Downstream, exactly that impulse produced four routes accumulating
volume four different ways with a pending decision record attached to each.
There is one correct functional here, and the other one is a bug.

### Step 1 — correct `dot_row`

Bound the loop with `size2`. Add unit tests for `dot_row` and `dot_col` on
non-square inputs; the existing bound is right for `dot_col`'s `3 x 1` usage
only by coincidence, and the next `1 x N` caller would hit the same trap.

> Gate: `dot_row` on a `1 x 3` input returns the full dot product.

### Step 2 — name the constant

Replace `0.16666666666` with `1.0 / 6.0` as a named
`kSignedVolumeQuadratureFactor`, commented with the `(1/3) x (1/2)` derivation.
The literal differs from `1/6` by about `4e-11` relative — physically
irrelevant, but it is the kind of stray value that later gets treated as a
compatibility fact worth preserving.

> Gate: no bare volume literal remains in the tree.

### Step 3 — gate the constraint on closure

Signed volume by the divergence theorem requires a closed, consistently
oriented surface. If `uVol != 0` and the mesh is open, has ghost faces, or is
periodic, throw at setup with an explicit message. This is the check whose
absence let the ambiguity persist.

> Gate: a flat/periodic fixture with `uVol != 0` fails at setup, not silently.

### Step 4 — rebaseline `vol0`

Any stored `vol0`, and any published result depending on one, sits on the old
scale. Emit both the old and the corrected volume to output for one release so
existing runs can be mapped across. Then drop the old one.

> Gate: every checked-in `input.params` with `uVol != 0` reviewed and rescaled.

### Step 5 — two permanent tests

**Origin independence.** Translate a closed fixture by a large offset;
`param.vol` unchanged to round-off. One assertion that catches the x-only bug,
an orientation flip, and a closure regression at once.

**Conjugacy.** Central finite difference of `energyVolume` against
`forceVolume`, per control point per axis, to `1e-6` relative.

> Gate: both in CI. Conjugacy is the test that should have failed all along.

## 4. Explicit non-goals

- **No `legacy-x-volume` mode.** A runtime switch between a correct functional
  and a buggy one guarantees that both live forever and that some future route
  picks the wrong default. Reproduce old runs by checking out the old tag.
- **No per-route or per-valence volume selection.** One functional, one
  accumulator, consumed identically by geometry, energy, and force.
- **No change to area.** `area += 0.5 * coeff * |a_1 x a_2|` is correct as
  written.
- **No change to the bending energy or its force.** Out of scope.

## 5. Sequencing

Do this **before** the irregular-patch work, and independently of it.

1. It is about ten lines of source plus tests, with no dependency on anything
   in the companion plan.
2. It moves the numerical baseline once. Landing it afterwards means moving the
   baseline twice, and re-running every convergence study.
3. Every irregular convergence measurement includes volume. Collecting that
   evidence against a functional that is off by a factor of three wastes the
   study.

The conjugacy test from step 5 is also the single most valuable test in either
plan — it is what validates the irregular force rows in the companion document,
and it needs the volume functional correct before it can pass on a closed mesh.

---

*The icosphere figures were computed for this document; they are a check on the
two functionals' ratio, not on SLIMED's Loop quadrature, which uses
limit-surface derivatives rather than flat triangles.*
