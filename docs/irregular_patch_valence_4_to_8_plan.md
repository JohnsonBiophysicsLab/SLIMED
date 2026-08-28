# Uniform Irregular Loop Patches (valence 4-8)

**Status:** implementation plan, not yet implemented
**Base:** `JohnsonBiophysicsLab/SLIMED @ 719025a`
**Companion:** [`volume_functional_split.md`](volume_functional_split.md) — land that first
**Scope:** extraordinary valence `N` in `{4, 5, 7, 8}` (6 is regular)
**Estimated new source:** ~900-1200 lines

One evaluation path for every extraordinary valence from 4 to 8, built as a
precomputed row table on top of the existing regular kernel. No OpenSubdiv, no
per-valence dispatch, no new subdivision scheme.

## 0. The shape of the change

This tree already carries every piece of this except the generalization. It has
the regular 12-control kernel, the 3-point Gauss rule, and — in
`get_subdivision_matrices()` — a working valence-5 recursion used for area and
volume. What it lacks is a way to build those matrices for any other valence,
and a way to feed them to the energy/force kernel.

The plan is to replace the hardcoded matrices with a generator, and to fold the
recursion into a **precomputed row table**. The key observation is that the
entire recursion is *pure topology*: it depends only on the valence `N`, never
on vertex coordinates. So it can be collapsed once at startup into a set of
`7 x (N+6)` matrices — exactly the same object the regular path already
consumes as `param.shapeFunctions`.

After that collapse, an irregular face costs the same as a regular one, the
force scatter is the same transpose, and `element_energy_force_regular` becomes
the single kernel for all valences with one change: replace the hardcoded `12`
with the row width.

## 1. What the tree does today

Five findings, all verified at `719025a`. Read before building.

### 1.1 Irregular faces fall through; unknown valences go silent

`src/energy_force/Compute_energy_and_force_on_mesh.cpp:83-111`

An 11-control face carries a `@todo` and then calls
`element_energy_force_regular()` anyway. Any face that is neither 12- nor
11-control matches no branch at all: `eBend`, `meanCurv`, and all three force
terms stay at their zero initializers and are stored as face state.

### 1.2 The regular kernel is not dimension-safe

`src/energy_force/Compute_energy_and_force_on_mesh.cpp:498`

The kernel hardcodes `for (int j = 0; j < 12; j++)` against a `7 x 12` shape
function. Feeding it an `11 x 3` control matrix is a dimension mismatch inside
a bare `gsl_blas_dgemm` wrapper (`Linear_algebra.cpp`, `multiplication()`) that
discards the return code.

### 1.3 The irregular predicate contradicts its own construction

`src/mesh/Mesh_setup_geometry.cpp:279` — `set_one_ring_vertices_sorted()`

The irregular branch is entered when *all three* corners have valence 5, but
the body then selects one corner as `d4` and walks the other two with the
valence-6 opposite-node pattern. The predicate describes an all-5 patch; the
construction describes a 5/6/6 patch. These are different topology classes and
the matrices only match the second.

### 1.4 Uninitialized reads in the same branch

`int d4, d7, d8;` are assigned inside a three-way `if` with no `else`. The
branch predicate tests `adjacentVertices.size()` while the selection tests
`adjacentFaces.size()`; for any vertex where those differ, all three are read
uninitialized. The winding test also uses `node0/1/2` coordinates rather than
the reassigned `d4/d7/d8`, and compares against the coordinate origin — which
is meaningless for a flat or off-origin membrane.

### 1.5 The valence-5 matrices are an asset

`src/mesh/Gauss_quadrature.cpp:224` — `get_subdivision_matrices()`

A direct port of the reference POC: `M` (17x11), `M1..M3` (12x17), `M4`
(11x17), with Warren's vertex weight `w = 3/(8N)`. Wired into
`calculate_element_area_volume()` for the 11-control case but never into
energy/force.

**This is the oracle the generator must reproduce bit-for-bit at N = 5.**

## 2. Why the POC does not generalize

Reference: `Yibenfu/continuum_membrane_model`,
`functions_for_triangular_mesh.cpp:430` and `myfunctions.cpp:95`.

The reference implementation writes `M` out as a literal 17x11 array. Every
entry encodes three things at once — the Loop masks, the canonical ordering of
the 11 control points, and the ordering of the 17 children — and none of the
three is recoverable from the literal. Extending to `N = 7` means hand-deriving
a 19x13 array with no way to check it short of the physics being wrong.

The structure is fine; only the encoding is the problem. Stam's reduction is
uniform in `N`:

```text
control points in the patch      K(N) = N + 6        // 12 at N=6, 11 at N=5, 14 at N=8
points after one subdivision     K(N) + 6 = N + 12
subdivision matrix               Abar : (N+12) x (N+6)
three regular children           P1,P2,P3 : 12 x (N+12)
one irregular child              P4 : (N+6) x (N+12)
```

At `N = 5` those shapes are 17x11, 12x17, and 11x17 — exactly the POC's
matrices. The existing code is one instance of a family, written out longhand.

## 3. Design: topology-only row tables

Each subdivision step peels off three regular children and hands the
extraordinary corner to a smaller copy of itself. Recursing `D` times tiles all
but a `4^-D` sliver of the parameter domain with regular patches.

```text
        level 0                 level 1                    level 2
      (irregular)        3 regular + 1 irregular      recurse, or truncate

           C                        C                          C
          / \                      / \                        / \
         /   \                    /P2 \                      /P2 \
        /     \                  *-----*                    *-----*
       /       \                / \ P3/ \                  / \ P3/ \
      /  K=N+6  \              /P4 \ /P1 \                 *---*---*  \
     V-----------B            V-----*-----B               V-.-*  P1   B
     ^                        ^                           ^
  valence N                extraordinary corner       residual 4^-d
                           stays in P4                of the domain
```

Every child but one is a regular 12-control patch the existing kernel already
evaluates. The composition `P_c * Abar * (P4 * Abar)^d` is integer-and-rational
topology — no coordinates enter, so it can be computed once per valence at
startup.

### 3.1 The collapse

Write `S_d = (P4 * Abar)^d` for the chain down to the depth-`d` irregular
child, and `W_q` for the existing `7 x 12` regular shape function at quadrature
point `q`. The limit-surface rows for child `c` at depth `d`, sample `q` are:

```text
R[N,d,c,q]  =  W_q * P_c * Abar * S_d          // shape 7 x (N+6)

geometry    G = R * X                          // X = (N+6) x 3 control coords
force       f = R^T * (dE/dG)                  // the transpose the regular path already uses
```

This is the whole design. `R` has the same type and meaning as
`param.shapeFunctions[q]`; it is simply wider or narrower than 12. Nothing
downstream needs to know a face was irregular.

### 3.2 Two things that are easy to get wrong

**No `4^-d` weight rescaling.** The derivative rows returned at depth `d` are
with respect to the *child's* parameters, so `|a_1 x a_2|` already carries the
shrunken metric and the children tile the surface exactly once. Mean curvature
is parameterization-invariant, so it needs no correction either. Applying an
extra `4^-d` would double-count the Jacobian.

**Table size is negligible.** Four valences x depth 6 x 3 children x 3 samples
x `7 x 14` doubles is under 200 kB, computed once at startup, immutable, and
shared across threads. It does not depend on the mesh, only on `N` — so unlike
a per-mesh stencil cache it never needs invalidation.

### 3.3 Building `Abar` without hand-deriving it

Do not write the matrices out. Build a synthetic canonical patch for valence
`N` in memory — the extraordinary vertex, its `N`-ring, and the surrounding
two-ring, as an explicit integer vertex/face list — then apply the generic Loop
masks to produce the `N+12` child points, and recover `P1..P4` by running the
*same* one-ring walk used on the real mesh over the subdivided local topology.
Roughly 300 lines, entirely integer and rational, and every step is checkable:

| Check | Expected |
| ----- | -------- |
| `N = 5` parity against `get_subdivision_matrices()` | exact equality, all four matrices |
| `N = 6` degeneracy — recursion vs. direct regular evaluation | agreement to round-off on a regular fixture |
| Affine invariance — row sums of `Abar`, `P_c`, and every `R` | `1.0` to `1e-14` |
| Planarity — `Abar` applied to a planar control net | stays in the plane |

Vertex weight stays Warren's `beta = 3/(8N)`, matching this tree and the POC,
exposed as one named constant so the alternative Loop rule is a one-line swap
rather than a rewrite. At `N = 6` both give `1/16`, so the all-regular baseline
is untouched either way. Note that `3/(8N)` is only the standard rule for
`N > 3` — another reason to keep `N = 3` out of scope.

## 4. Admission, and the isolation requirement

The reduction assumes **exactly one extraordinary corner**. A patch with two or
three has no such decomposition — the three "regular" children stop being
regular. This is the real constraint on the plan, and it needs a stated policy
rather than a silent fallback.

A face is admitted when it is:

- interior and non-ghost;
- two-manifold with complete one-rings on all three corners;
- carrying exactly one corner with valence in `{4, 5, 7, 8}`;
- with the other two corners at exactly 6.

Everything else throws at setup, naming the face index and its three valences.

### 4.1 When extraordinary vertices are adjacent

1. **Require isolation of the input mesh**, validate at load, and report a
   valence histogram plus a count of faces with more than one extraordinary
   corner. This is the default and it is usually satisfiable at
   mesh-generation time.
2. **One global Loop refinement at setup**, behind an explicit off-by-default
   parameter. Every vertex introduced by a subdivision step has valence 6 and
   every old vertex keeps its valence, so a single pass provably separates all
   extraordinary vertices. The cost is real: 4x the vertices, and the control
   mesh — which *is* the dynamical degrees of freedom — changes, so it
   rebaselines the run.
3. **Reject.** Fail loudly rather than approximate.

Do not attempt local refinement around adjacent extraordinary vertices: Loop
rules are uniform, and refining part of a control mesh changes the limit
surface elsewhere.

Separately, replace the origin-referenced winding test with a consistent
orientation derived once from the half-edge structure in `HalfedgeMesh.cpp`.
The current test is only correct for a closed surface enclosing the coordinate
origin.

## 5. Work packages

### WP1 — Guardrails first

Fix the `set_one_ring_vertices_sorted()` predicate to 5/6/6 (matching what the
code actually builds), initialize `d4/d7/d8` or fail, replace the winding test,
and add explicit rejection plus a valence histogram. No new capability — this
only converts silent wrong answers into errors.

> Gate: `data/example` run bit-identical to the current baseline.

### WP2 — Canonical patch topology and the `Abar` generator

Synthetic patch construction, generic Loop masks, selection matrices for all
four children, for any `N`. Pure integer/rational, no mesh dependency.

> Gate: bit-exact vs. `get_subdivision_matrices()` at `N=5`; row sums `1.0`
> for `N = 4..8`.

### WP3 — Row table

`IrregularPatchRowTable` keyed by `(N, depth, child, sample)`, built at
startup, immutable, thread-shared. Depth is a parameter with a per-valence
default set by WP6.

> Gate: `N=6` rows reproduce `param.shapeFunctions` exactly.

### WP4 — One kernel

Generalize `element_energy_force_regular` to take a row matrix and derive its
width from `R.ncol()`. Drive regular and irregular faces through the identical
call. Delete the 11-control special case and the `@todo` fallthrough.

> Gate: regular path bit-identical; no branch on valence below the dispatch.

### WP5 — Area and volume through the same rows

Replace the case-11 recursion in `calculate_element_area_volume()` with a
row-table lookup. `param.subMatrix` and `get_subdivision_matrices()` retire —
but only after WP2 has copied them into the test suite as the permanent `N=5`
oracle.

> Gate: geometry and energy/force read the same rows for the same face.

### WP6 — Convergence study, per valence

Sweep depth `D` for each `N` and measure area, volume, bending energy, and
force conjugacy. Fixtures: octahedron (`N=4`), icosahedron (`N=5`); `N=7` and
`N=8` need constructed closed fixtures — a subdivided torus, or a hand-built
disk with a prescribed centre valence and a fixed boundary. Report the
truncated parameter fraction and the dropped surface area as run diagnostics.

> Gate: monotone convergence; a chosen `D` per valence, recorded with its
> residual.

### WP7 — Optional pre-refinement

The off-by-default global Loop refinement from section 4, for meshes that
cannot be generated with isolated extraordinary vertices. Only worth building
if WP1's histogram shows real meshes need it.

> Gate: refined mesh passes the isolation validator by construction.

## 6. Tests that actually bind

Four of these would have caught every defect listed in section 1, and they are
cheap enough to run in CI on every commit.

- **Force-energy conjugacy.** Per face, per control point, per axis: central
  finite difference of total energy against the reported force, to `1e-6`
  relative. The single strongest test in the plan, and it applies unchanged to
  every valence.
- **Rigid-motion invariance.** Translate and rotate a fixture; energy
  unchanged, net force and net torque zero to round-off.
- **Regular regression.** The `data/example` workload — all-regular physical
  faces — must stay bit-identical through WP1-WP5.
- **`N=5` matrix parity.** The generator against this tree's literal matrices,
  exact.
- **Convergence.** Energy, area, and volume against depth, with the expected
  `~4^-D` residual.
- **Loud rejection.** A mesh with two adjacent extraordinary vertices throws,
  naming the face; it must never silently produce zeros.

## 7. Where OpenSubdiv fits

Nowhere in production. Everything above is a few hundred lines of rational
arithmetic over local topology, and the strongest validation available —
conjugacy, invariance, and exact `N=5` parity against matrices already in the
tree — needs no external library.

If an independent cross-check is wanted later, it belongs in a standalone test
binary that is not linked into the simulator. Note that it would need its
vertex mask configured to match: stock OpenSubdiv Loop uses
`beta = (1/N)(5/8 - (3/8 + cos(2*pi/N)/4)^2)`, which is `0.0841` at `N = 5`
against SLIMED's `0.075`. An unconfigured comparison disagrees by construction
and proves nothing.

For scale: this plan is roughly 900-1200 lines of new source. The equivalent
program in `SLIMED_dev` is about 8,900 lines of source across three per-valence
routes, three build flags, and three runtime opt-ins, and covers three Platonic
fixtures rather than a valence range.

## 8. Risks and open questions

| Risk | Handling |
| ---- | -------- |
| **Quadrature near the extraordinary point.** Loop surfaces are C1 but not C2 at extraordinary vertices; bending energy is the integrand most sensitive to the truncated corner. A 3-point rule may converge slowly for `N` far from 6. | WP6 measures it per valence rather than assuming. If depth alone is insufficient, the next step is a higher-order symmetric triangle rule on the regular children — not an adaptive rule, which would make energy discontinuous during a timestep. |
| **`N = 7, 8` fixtures.** No Platonic solid gives valence 7 or 8, so closed test fixtures must be constructed. | Scope WP6 to include fixture construction. A closed genus-1 mesh admits both; a fixed-boundary disk with a prescribed centre valence is simpler and enough for row-level validation. |
| **Isolation may not hold on real inputs.** | WP1's histogram answers this before WP7 is scoped. Decide with data. |
| **Volume is measured wrong today**, so any convergence evidence collected before that is fixed has to be re-run. | Do the companion plan first. It is about ten lines. |

### Decisions needed before WP2

- Confirm Warren's `3/(8N)` stays the vertex mask — it matches this tree and
  the POC, so this is a confirmation rather than a change.
- Confirm the valence range stops at 4 and 8, and that out-of-range input is an
  error rather than a fallback.
- Confirm that a mesh failing isolation should throw at setup rather than run
  degraded.
