# Monte Carlo Edge Flips for a Fluid Membrane

**Status:** work packages 0-7 landed
**Base:** `JohnsonBiophysicsLab/SLIMED @ 1fdffbd`
**Builds on:** [`irregular_patch_results.md`](irregular_patch_results.md) (valence 4–8
row tables), [`fluctuation_spectrum.md`](fluctuation_spectrum.md) (the end-to-end
check this work must keep passing)
**Scope:** in-plane fluidity by Metropolis edge flips, scheduled as a Poisson
process in physical time, interleaved with the existing Brownian dynamics.
Fixed surface topology — no fusion or fission (the cut-and-paste moves of
Gompper & Kroll 1998 are out of scope).

The control mesh of a subdivision surface is a solid: two vertices that are
neighbours stay neighbours forever, so the membrane carries an in-plane shear
modulus it should not have, and a large-scale shape change has to drag the
connectivity with it. Every dynamically triangulated surface (DTS) model
removes that shear modulus the same way — by flipping the diagonal of a pair
of adjacent triangles as a Monte Carlo move — and this plan does the same
thing on a Loop limit surface.

Three things make it different from a textbook DTS flip, and they set the
shape of the work:

1. **The energy of a face is not local to the face.** It lives on the union of
   its three corners' one-rings, so a flip changes the energy of every face
   touching any of the four vertices involved, about eighteen faces.
2. **A flip creates faces with two or three extraordinary corners, always.**
   On a hexagonal lattice a flip turns four valence-6 vertices into 5, 5, 7, 7,
   and the two new triangles have three extraordinary corners each. The row
   tables in the tree admit exactly one. This is not a corner case; it is
   every flip. Handling it is a prerequisite, not a follow-up.
3. **The dynamics integrates the limit surface, not the control net.** The
   limit mask depends on valence, so after flips the conversion matrix must
   be rebuilt — and it is currently a dense explicit inverse.

The plan is ordered so that each of these is retired by a work package with a
binding test before the flip move itself is wired in.

---

## 1. Theoretical basis

### 1.1 The flip move in dynamically triangulated surfaces

A DTS represents the membrane as `N` vertices joined into a triangulation `T`.
The statistical weight is the Boltzmann factor of the discrete Hamiltonian,
summed over all triangulations with equal prior weight and integrated over
vertex positions (Gompper & Kroll 2004; Ramakrishnan, Sunil Kumar &
Radhakrishnan 2015):

```text
Z = Σ_T ∫ Π_v dX_v  exp(-E(X, T) / kT)
```

Two Monte Carlo moves sample this: a vertex displacement at fixed `T`, and a
**link (bond, edge) flip** at fixed `X`. In a flip the shared edge `(i, j)` of
two adjacent triangles `(i, j, k)` and `(i, l, j)` is removed and replaced by
`(k, l)`, giving triangles `(i, l, k)` and `(j, k, l)`. The four valences change
by `-1, -1, +1, +1`; the number of vertices, edges and faces does not. The move
is accepted with the Metropolis probability `min[1, exp(-ΔE/kT)]` (Metropolis
et al. 1953). With a
uniformly chosen edge, the proposal is symmetric — the reverse flip is proposed
from the new state with the same probability `1/N_E`, since the flip is an
involution and `N_E` is conserved — so no proposal-ratio correction is needed
(Ramakrishnan, Sunil Kumar & Ipsen 2010, eq. 19; Gompper & Kroll 2004). (there
is an important subtlety: it is the uniform choice among all edges that makes
the proposal symmetric. If you choose uniformly only among currently flippable
edges, the statement is generally no longer true. )

A Monte Carlo sweep in the literature is `N` vertex moves plus one flip attempt
per edge, i.e. `3(N-2)` attempted flips on a closed surface (Ramakrishnan et
al. 2010). TriMem instead sweeps a fraction `γ` of the edges (their target flip
rate, `γ = 10%`) after each hybrid-MC step, and reports the *acceptance*
`ε = accepted / (γ |E|)` separately (Siggel et al. 2022, eq. 19) — with a stiff
tether potential `ε` is well under one percent, and they note that this forces
flips to be evaluated **sequentially**: accepting several flips at once would
multiply the individual acceptance probabilities.

Validity checks before a flip is even proposed are standard and are what keeps
the triangulation a two-manifold:

| Check | Why | Source |
| ----- | --- | ------ |
| Both incident faces exist and are real | boundary / ghost edges are not flippable | all |
| `k` and `l` not already joined | would create a duplicate edge | Gompper & Kroll |
| Valence bounds after the flip | a vertex with fewer than 3 neighbours is not a manifold corner; very high valence degrades the discretization | Gompper & Kroll (3 … 9); FreeDTS |
| Edge length bounds | tether: `l_min ≤ l_e ≤ l_max`, e.g. `l_min = l_dts`, `l_max = 3 l_dts` in FreeDTS, `√3 a_0` in the hard-core models | FreeDTS, Ramakrishnan |
| Dihedral / normal angle | "a mild constraint on the dihedral angle between two faces sharing a tether restores self avoidance" | Ramakrishnan; OrganL (`n_p · n_i > 0.01`) |

In SLIMED the row tables replace the valence bounds with a hard range, `[4, 8]`
today (`Subdivision_matrices.hpp:37-38`), and the edge-length and angle
constraints are replaced by the energy itself — a distorted patch costs
bending and spring energy and is rejected by Metropolis rather than by a
cutoff. Section 3.3 spells out which checks survive.

### 1.2 Flips alongside continuous dynamics

Flips were born in pure Monte Carlo, but every dynamical DTS model combines
them with a continuous integrator, and the combination is what SLIMED needs:

- **Noguchi & Gompper (2004, 2005)** drive the vertices with molecular /
  multiparticle-collision dynamics and flip tethers by Monte Carlo. "These
  bond flips provide also a convenient way to vary the membrane viscosity
  `η_mb`, because it increases with decreasing bond-flip rate." The flip rate
  is a physical knob, not a numerical one.
- **Sadeghi, Weikl & Noé (2018)** integrate a particle membrane with Langevin
  dynamics and propose flips "with a frequency `φ`". The proposal frequency is
  "a control parameter for the model": the in-plane diffusion and the surface
  viscosity are set by it, and they measure `η = η_∞ exp(C_φ / φ)` from a 2D
  Poiseuille flow. They also make the thermodynamic point precisely: an
  accepted flip that lowers the energy extracts work that the thermostat
  replaces as heat — "an entropy production mechanism comparable to viscous
  loss."
- **TriMem** (Siggel et al. 2022) alternates hybrid-MC trajectories with flip
  sweeps and verifies Boltzmann sampling by reproducing the vesicle phase
  diagram.
- **PyMembrane** (2023) offers a fixed-connectivity elastic membrane and a
  bond-flipping liquid one in one framework — the same pair of modes
  §3.10's `edgeFlipEnabled` switches between.
- **FreeDTS** (Pezeshkian & Ipsen 2024) puts `N_T` flip attempts, `N_v` vertex
  updates and the inclusion moves into one MC step, and folds the *global*
  constraints into the local energy change of each move.
- **OrganL (2024)** is the closest analogue to this tree — curved (Nagata
  cubic) elements with dynamic triangulation — and reports that after a flip
  round "neighbor lists [are] rebuilt", that the remeshing move "is executed
  only occasionally and (mostly) serially", and that a global quadratic
  restraint `λ/2 (Φ - Φ_0)²` contributes `ΔF_Φ = λ/2 ΔΦ (ΔΦ - 2(Φ - Φ_0))` to a
  local move — exactly the form SLIMED's area and volume constraints take.

**Why the combination is sound.** The Brownian step at fixed `T` is a Markov
kernel whose stationary density is `exp(-E(·, T)/kT)` up to the usual `O(Δt)`
integrator error, which the tree already accepts (`fdtConsistentSurfaceUpdate`
exists to make that statement true — see `fluctuation_spectrum.md`). The flip
step is a Metropolis kernel with the same stationary density on the joint
space of positions and triangulations. A composition of kernels that each
leave a measure invariant leaves it invariant. The composite chain is not
reversible — alternating two reversible kernels never is — but stationarity is
all that equilibrium sampling requires, and the ordinary detailed-balance
proofs in the DTS literature apply to each kernel on its own.

### 1.3 Scheduling flips as a Poisson process in physical time

What the literature leaves implicit, this plan makes explicit. A flip attempt
rate `ν` per edge per unit time defines a Poisson process on every edge. Over a
Brownian step `Δt` on a mesh with `N_E` flippable edges, the number of attempts
is

```text
n ~ Poisson(λ),      λ = ν · Δt · N_E                                    (1)
```

with the edges drawn uniformly at random (with replacement — an edge may be
attempted twice, harmlessly). Drawing `n` from a Poisson rather than fixing it
is what makes the flip dynamics independent of `Δt`: halving the step halves
`λ`, and the physical attempt rate per edge is unchanged. That is the sense in
which "the number of flips conforms to a distribution scaled with the time step
and the size of the membrane."

`ν` is the fluidity parameter. It is not free: the accepted-flip rate sets the
membrane's in-plane viscosity (Noguchi & Gompper 2004, 2005; quantified as
`η = η_∞ exp(C_φ/φ)` by Sadeghi, Weikl & Noé 2018) and the vertex diffusion,
and should be calibrated (WP6) against the physical neighbour-exchange time of
a mesh vertex. A vertex at `lFace = 5 nm` stands for a patch of order a hundred
lipids; with a lipid diffusion constant of `1–10 nm²/µs` — the fluid-phase
range reviewed by Almeida & Vaz (1995), `10⁻⁸–10⁻⁷ cm²/s` — the time for a
patch to exchange a neighbour is `l² / 4D ≈ 0.6–6 µs`, so `ν` in the range
`0.1–1 /µs` is the physical starting point. For the shipped `input.params`
(`timeStep = 0.001 µs`) on the `data/example` mesh (3,680 faces, `N_E ≈ 5,500`)
that is `λ ≈ 0.5–5` attempts per step. The cost of a flip attempt is a local
energy evaluation (§3.4), so this is cheap; the expensive consequence of
fluidity is elsewhere (§6, risk 1).

### 1.4 What is different on a subdivision surface

**(a) The flip patch.** The energy of face `f` is a functional of its control
net, the union of the one-rings of its three corners
(`Face::oneRingVertices`). A flip of `(i, j) → (k, l)` changes the one-ring of
`i`, `j`, `k` and `l`, so it changes the control net — and the energy — of every
face incident to any of the four. That set, the **flip patch**, is
`|F(i)| + |F(j)| + |F(k)| + |F(l)|` minus overlaps, about 18 faces on a
near-regular mesh. It is the same object TriMem locks for a parallel flip
(Siggel et al. 2022, Fig. 2, which notes it is "significantly larger than the
patch required for a flip subject to the Delaunay criterion"). Everything
outside the flip patch has an unchanged control net and an unchanged energy;
that is what makes `ΔE` local.

**(b) Faces with several extraordinary corners are unavoidable.** Stam's
evaluation (Stam 1998), and the row tables built on it, need exactly one
extraordinary corner per face. That is achievable at mesh generation, and
WP7's global refinement makes any mesh satisfy it. A flip destroys it
immediately:

```text
       k                    k
      / \                  /|\
     /   \                / | \
    i-----j     -->      i  |  j        valences: i,j : 6 -> 5
     \   /                \ | /                   k,l : 6 -> 7
      \ /                  \|/
       l                    l
                          faces (i,l,k), (j,k,l): three extraordinary corners each
```

And it cannot be repaired by refining: the tree's own documentation is right
that a partial refinement changes the limit surface elsewhere
(`Mesh_refine.cpp:5-9`). The way out is Stam's original observation, applied
per face rather than per mesh: **one Loop subdivision of a face's own control
net produces four children, each with at most one extraordinary corner**, and
the limit surface is unchanged because Loop subdivision (Loop 1987) is exactly
the map the limit surface is defined by. The child at corner `v_a` keeps the
valence of `v_a` and gets two valence-6 neighbours; the centre child is
regular.

The key fact that makes this local and cheap: **every control point of every
child lies inside the parent's own patch.** A child's corner is either an old
corner (whose update needs its one-ring — in the patch) or an edge midpoint of
an edge `(v_a, x)` incident to a corner (whose update needs the two opposite
vertices of the faces on that edge — both in the one-ring of `v_a`, hence in
the patch). No two-ring is needed. So for a face with corner valences
`(N_0, N_1, N_2)` and patch width `K = N_0 + N_1 + N_2 - 6` (12 for a regular
face, `N + 6` for the single-extraordinary case the tree already handles) there
is a local subdivision matrix `A_loc` of size `(N_0 + N_1 + N_2) × K`, and the
rows the kernel needs are

```text
corner child a (valence N_a ≠ 6):  R[N_a, d, c, q] · P_a · A_loc     // existing table, then two more selections
corner child a (valence 6):        W_q · P_a · A_loc
centre child:                      W_q · P_ctr · A_loc
```

which have the same type and meaning as `param.shapeFunctions` and the
existing `R[N, d, c, q]` — `7 × K` blocks — and go through the same
`element_energy_force_patch_pod()` with `nCtrl = K`. Pure topology, as before:
`A_loc` depends on the valence triple and the ring overlap pattern, never on
coordinates. WP1 builds it.

**(c) Which coordinates carry the flat measure.** The Brownian step displaces
the limit-surface points `S = M C` with isotropic noise, so at fixed `T` it
samples `exp(-E) dS`. A flip at fixed control net `C` with plain Metropolis
samples `exp(-E) dC` across triangulations. These differ by the Jacobian
`|det M_T|`, which depends on `T`. With Warren's weights (Warren & Weimer 2001)
the limit mask row of a valence-`N` vertex is `1/2` on itself and `1/(2N)` on
each neighbour, so `M = ½ D⁻¹ (D + A)` and

```text
ln det M_T = -Σ_v ln(2 N_v) + ln det(D + A)
```

For a flip `6,6,6,6 → 5,5,7,7` the first term changes by `+0.056` and the
second by a smaller, local amount, so the two measures disagree by roughly
`0.05–0.08` in `ln` weight per flip — a few-percent bias on the relative
weights of triangulations, no bias at all on the geometry sampled at fixed
`T`. Three consistent choices exist:

1. flip at fixed `C`, plain Metropolis — the standard DTS model, every reference
   above, and what this plan implements first;
2. flip at fixed `C`, acceptance multiplied by `|det M_T'| / |det M_T|`, evaluated
   by the matrix determinant lemma on the four changed rows
   (`det(I + Vᵀ M⁻¹ U)`, a 4×4 determinant after four sparse solves);
3. flip at fixed `S` — the material points stay put and the control net is
   re-solved, `C' = M_T'⁻¹ S`; no Jacobian, physically the most natural, but
   `C'` changes in a neighbourhood of the flip and so does the patch that
   `ΔE` must cover.

The plan takes (1), records the size of the discrepancy with a test that
evaluates the log-determinant ratio for a lattice flip (WP3), and leaves (2)
behind a parameter if the measured valence distribution is ever the quantity
of interest. It does not take (3) in this round.

**(d) The limit mask is valence-dependent.** `assign_mesh2surface()` writes the
valence-6 mask `(1/2, 1/12, …)` for every vertex (`Dynamic_mesh.cpp:113-123`).
After flips a vertex of valence `N` needs `(1/2, 1/(2N), …)`. With unequal
valences `M` is no longer symmetric, so the fluctuation–dissipation mapping
must use `M⁻ᵀ` for the force, which the tree already documents as its standing
approximation (`Dynamic_model.cpp:107-138`). Both fixes come with WP4.

### 1.5 Global constraint terms in the local energy change

SLIMED's area and volume constraints are quadratic in the *totals*
(`Compute_energy_and_force_on_mesh.cpp:365-385`):

```text
E_A = (uSurf / 2 area0) (A - area0)²,       E_V = (uVol / 2 vol0) (V - vol0)²
```

A flip changes `A` and `V` only through the flip patch, so with `ΔA` and `ΔV`
summed over the patch the exact change is the OrganL (2024) / FreeDTS
(Pezeshkian & Ipsen 2024) form:

```text
ΔE_A = (uSurf / 2 area0) · ΔA · (ΔA + 2 (A - area0))                      (2)
```

and likewise for volume. This needs the per-face areas and volumes the tree
already stores (`Face::elementArea`, `elementVolume`) and the running totals
`param.area`, `param.vol`, updated after every accepted flip within a sweep so
the next attempt sees the right `A`.

### 1.6 Mesh quality: the tether becomes an edge spring

DTS models keep triangles well-shaped with a tether potential — hard walls in
the Monte Carlo models (Gompper & Kroll 2004), a smooth well in the MD ones
(Noguchi & Gompper 2004 introduced "a smooth bond-interaction potential, which
makes the model amenable for molecular dynamics"; TriMem's eq. 14, Siggel et
al. 2022, is a continuous version for the same reason). The Metropolis flip
then automatically prefers the Delaunay-like diagonal, because the other one is
longer and costs tether energy.

SLIMED's regularization (`energy_force_regularization()`,
`Compute_energy_and_force_on_mesh.cpp:802-975`) is a spring on each face's three
edge lengths against **the same face's edge lengths in `coordRef`**. That is a
solid's memory of its reference configuration: it is meaningless for a fluid,
and a freshly created edge has no reference length at all. Two further
findings: the term uses `param.kCurv` rather than the `kReg` it was given
(`:804`), and its force is not part of the Brownian drive (`Dynamic_model.cpp:145-146`
uses `forceCurvature + forceArea` only), so nothing currently holds the
in-plane mesh in shape during dynamics — which is invisible today only because
the in-plane displacement is zeroed (`Dynamic_model.cpp:203-206`).

Fluid mode therefore uses a **connectivity-independent edge spring**

```text
E_S = (k_S / 2) Σ_{edges e} (l_e - l_0)²,     l_0 = lFace by default             (3)
```

as the Hamiltonian's mesh-quality term, with its force included in the
Brownian drive, and the flip's `ΔE` includes the spring change of the one
removed and one added edge. Its physical footprint is a small contribution to
the area modulus and, with flips, none to the shear modulus; the fluctuation
spectrum fit in `membrane_spectrum.py` reports an effective tension against
the Helfrich `q⁻⁴` law (Helfrich 1973; see
[`fluctuation_spectrum.md`](fluctuation_spectrum.md)), which is the check that
the spring is not adding one (WP6).

---

## 2. What the tree does today

Verified at `1fdffbd`. Items marked **[blocking]** must change before a flip can
be committed at all; the rest are what a correct implementation must respect.

### 2.1 Connectivity is built once and never touched

- **No edge structure. [blocking]** `HalfedgeMesh.cpp` and `Edge.cpp` are
  entirely commented out; `Mesh.hpp:214-216` has `edges`/`halfedges` commented
  out. Edge maps are built ad hoc and discarded three times
  (`Mesh_setup_geometry.cpp:133-166`, `Mesh_setup_boundary_condition.cpp:105-152`,
  `Mesh_refine.cpp:40-52`).
- `Vertex::adjacentVertices` / `adjacentFaces`, `Face::adjacentFaces` and
  `Face::oneRingVertices` are populated by `setup_flat()`
  (`Mesh_setup_flat.cpp:3-22`) and never updated (`Update.cpp` copies
  coordinates, forces and energies only). `setup_from_vertices_faces()`
  (`Mesh.cpp:74-144`) skips `set_adjacent_faces_of_faces()` and
  `sort_vertices_on_faces()` entirely, so imported meshes have empty
  `Face::adjacentFaces`.
- `set_adjacent_faces_of_vertices_sorted()` sorts with hardcoded lattice index
  arithmetic (`Mesh_setup_geometry.cpp:36-74`); `set_adjacent_vertices_of_vertices_sorted()`
  does not sort at all (`:216-254`). Nothing downstream depends on either order
  except `Vertex::approximate_unit_normal_vector`, which is order-independent.
- There is no valence field; valence is `adjacentVertices.size()` everywhere,
  and downstream code infers a face's valence from its patch width,
  `valence = nOneRingVertices - 6` (`Compute_energy_and_force_on_mesh.cpp:175`,
  `Mesh.cpp:308`, `Device_mesh_layout.cpp:74`). A multi-extraordinary face
  cannot be described this way. **[blocking]**

### 2.2 Admission is a throw, not a predicate

`set_one_ring_vertices_sorted()` (`Mesh_setup_geometry.cpp:360-565`) classifies
each face as 6/6/6 or `N`/6/6 with `N ∈ [4, 8]`, rotates the corners so the
extraordinary one anchors the fan walk, permutes the ring into
`canonical_control_order()`, and **throws** at the end of the loop for anything
else (`:542-559`). There is no function that answers "would this face be
evaluable?" without side effects. **[blocking]**

### 2.3 Energy: local per face, global in the constraints

- Per-face energies are stored (`Face::energy.energyCurvature` at `:208`,
  `energyRegularization` at `:950`); `energyArea` and `energyVolume` on faces
  are always zero — those terms are computed once from the totals.
- The bending/area/volume kernel is allocation-free and takes plain arrays:
  `element_energy_force_patch_pod()` and `element_area_volume_pod()`
  (`Patch_kernel.hpp:255-299`). Stack buffers are sized by
  `kMaxControlPoints = 14` (`Patch_kernel.hpp:41`); a three-extraordinary face
  at valence 8 has `K = 18`, and a child of it after local subdivision is
  evaluated on the same `K`-wide net. **[blocking]**
- Nothing evaluates a subset of faces. The force loop is one pass over all
  faces (`compute_face_energies_and_forces()`, `:57-267`).
- `PatchParams` is filled with the global `area`, `area0`, `vol`, `vol0`
  (`:90-97`); the area force on every face scales with `(A - area0)`.

### 2.4 Dynamics

- Loop (`Run_dynamics_flat.cpp:95-133`): `S = M C` → `next_step()` on `S` →
  `C = M⁻¹ S` → periodic post-process → copy to vertices → output →
  **`Compute_Energy_And_Force()` at the end of the iteration** (`:130`).
- `mesh2surface` is a dense `nV × nV` matrix, explicitly inverted once
  (`Dynamic_mesh.cpp:27`); the memory note records that the dense multiply
  already dominates the step time and does not thread. Re-inverting after a
  flip is `O(N³)` and would dominate everything. **[blocking]**
- `x` and `y` displacements are multiplied by zero (`Dynamic_model.cpp:203-206`);
  the drive is `forceCurvature + forceArea` only (`:145-146`, `:192-193`).
- Ghosts and periodic duplicates are skipped (`:173-179`); ghost/slave status
  is by lattice position (`Mesh_setup_boundary_condition.cpp:5-101`,
  `Boundary_processing.cpp:9-39`), not by connectivity, so an interior flip
  does not disturb it.
- RNG is counter-based, keyed by `(seed, iteration, vertex, axis)`
  (`Dynamic_model.cpp:8-37`, `:154-156`). A flip RNG can use the same generator
  with its own tag and be reproducible without state.

### 2.5 Caches that will silently go stale

`ensure_device_layout()` rebuilds only when the face or vertex *count* changes
(`Compute_energy_and_force_on_mesh.cpp:269-280`); `CudaForceBackend::topologyUploaded`
is set once (`Cuda_force_backend.cu:362`); `ensure_patch_rows_flat()` checks
`empty()` only (harmless — the row tables are valence-keyed and do not change).
A flip changes no count. **[blocking for the GPU path]**

### 2.6 Nothing persists connectivity

`face.csv` is written once (`Run_dynamics_flat.cpp:65`); the trajectory is
coordinates only (`output.cpp:505-560`); `slimed_restart.chk` stores no faces and
its reader checks only the vertex count (`output.cpp:308-313`), so a restart of a
flipped run would load coordinates onto the wrong topology.

### 2.7 There is a Metropolis move already

`Model::simulated_annealing_next_step()` (`Energy_minimization.cpp:293-403`) is a
single-vertex Metropolis displacement that snapshots every coordinate, runs a
**full** `Compute_Energy_And_Force()`, and restores everything on rejection. It
is the right acceptance logic and the wrong cost model; the flip must be local.
Its test, `ThermalFluctuationTest.EnabledSwitchAttemptsMetropolisTrial`, has
been failing since before the irregular-patch work.

---

## 3. Design

### 3.1 Edge table

A persistent, undirected edge table on `Mesh`, built once from the faces and
updated in place by a flip:

```cpp
struct Edge {
    int v[2];        // endpoints, v[0] < v[1]
    int face[2];     // incident faces; -1 for a boundary edge
    int opposite[2]; // the third corner of face[k], i.e. the two flip targets
    bool flippable;  // static part of the admission test (§3.3), set at build
};
std::vector<Edge> edges;
std::unordered_map<uint64_t, int> edgeIndex;   // undirected_edge_key -> edge id
```

`undirected_edge_key()` already exists (`Mesh_setup_geometry.cpp:133-138`).
The three ad hoc edge maps are replaced by this one. The table also serves the
edge spring (§3.7) and the per-frame flip log.

### 3.2 The flip primitive

`Mesh::flip_edge(int e)` performs the topology change and nothing else. With
`F1 = (i, j, k)` and `F2 = (i, l, j)` both wound consistently (the boundary of
the quadrilateral runs `i → l → j → k → i`):

```text
F1 := (i, l, k)      F2 := (l, j, k)      edge (i,j) := (k,l)
```

Both new faces are wound along the same boundary, so orientation is preserved
without any geometric test. The bookkeeping, all `O(valence)`:

- `adjacentVertices`: remove `j` from `i`, `i` from `j`; add `l` to `k`, `k` to `l`.
- `adjacentFaces` of `i`, `j`, `k`, `l`: `F2` leaves `i`, `F1` leaves `j`; `F2`
  joins `k`, `F1` joins `l`.
- `Face::adjacentFaces` of `F1`, `F2` and of the two outer neighbours that
  change sides; `Edge::face` / `opposite` of the four boundary edges of the
  quadrilateral.
- Then `rebuild_patches(flipPatch)`: re-run the one-ring construction and the
  classification of §3.3 for every face incident to `i`, `j`, `k`, `l`.
- `topologyVersion++` — the counter every cache checks (§3.6).

`flip_edge(e)` applied twice restores every one of those, with one exception
found while building it: **the two incident face indices come back holding
each other's triangle.** That is intrinsic, not an implementation artefact.
The quadrilateral offers no canonical pairing between "the side of `c0`"
before the flip and either side after it, so every consistent assignment rule
composes to the exchange; OpenMesh's `flip()` does the same (Botsch et al.
2002). It is
unobservable here because the only per-face state that survives a step without
being recomputed from connectivity is the spontaneous curvature, and the
admission test refuses an edge whose two faces disagree about it. The
involution test asserts the exchange rather than papering over it: it undoes
the relabelling explicitly and then demands exact equality.

### 3.3 Admission predicate and multi-extraordinary patches

The classification in `set_one_ring_vertices_sorted()` is split into a pure
function

```cpp
enum class PatchKind { Regular, SingleExtraordinary, MultiExtraordinary, Inadmissible };
PatchKind classify_face(int f, std::string *why = nullptr) const;
```

used both by setup (which keeps throwing on `Inadmissible`) and by the flip
sweep (which rejects instead). A face is admissible when it is non-ghost, all
three corners are interior (`is_interior_vertex()`), and every corner valence
is in the row-table range `[kMinIrregularValence, kMaxIrregularValence]`. The
"other two corners at 6" condition is dropped: that is what WP1 removes.

The static part of the flip test (`Edge::flippable`, computed once): both
incident faces real and non-ghost; none of the four vertices, nor any corner of
any face in the flip patch, is a ghost or a periodic duplicate
(`isSlavedPeriodic`); no face in the flip patch holds a scaffolding point
(`Particle::faceIndex` would need remapping — out of scope). The dynamic part,
evaluated per attempt: valences `N_i - 1`, `N_j - 1 ≥ N_min` and
`N_k + 1`, `N_l + 1 ≤ N_max`; `k` and `l` not already adjacent. With `[4, 8]`
that is more restrictive than the DTS convention `[3, 9]` (Gompper & Kroll
2004); extending the generator to `N = 3` (Loop's `β = 3/16` special case,
Loop 1987) and `N = 9, 10` is a
one-day follow-up once flips are running and the measured valence histogram
says whether it matters.

**Multi-extraordinary rows (WP1).** For a face of kind `MultiExtraordinary`,
`rebuild_patches` constructs, from the live adjacency:

1. the patch `P`: ordered union of the corners and their one-rings, width
   `K = N_0 + N_1 + N_2 - 6` in the generic case (fewer if rings overlap
   beyond the shared corners; the construction uses the actual sets);
2. `A_loc` (`(N_0 + N_1 + N_2) × K`): Loop's vertex rule with `β = 3/(8N)` on the
   three corners, the edge rule `(3/8, 3/8, 1/8, 1/8)` on every edge incident to
   a corner (Loop 1987) — every stencil entry is in `P` by §1.4(b);
3. four child selections `P_c`, produced by running the *same* one-ring walk
   used on the real mesh over the subdivided local topology — the way
   `build_canonical_patch()` already does for the canonical patch;
4. the composed rows per child and quadrature sample, stored flat per face.

Depth for a corner child of valence `N` is `recommended_irregular_depth(N)`,
unchanged. Cost per face: `3 · (Σ_a D_{N_a}) + 3` samples instead of 3 — for a
5/5/7 face at the current depths, about 200 samples, roughly 60× a regular
face. §6 risk 1 addresses what that means at scale.

`A_loc` is pure topology and its composition with the existing table is a
sequence of small dense products, so caching is by a per-face **topology
signature** (the ordered valence triple plus the ring-overlap pattern); faces
with the same signature share one row block. On the flat example mesh with
flips running, a few dozen signatures cover everything.

**Correctness test that binds.** Loop subdivision commutes with itself, so for
any face `f` the sum of the energies of its four children on the globally
refined mesh (`refine_loop_once()`, WP7) must equal the energy of `f` on the
coarse mesh evaluated through `A_loc` — at the same corner depths, to round-off.
This is an exact identity, not a convergence claim, and it tests `A_loc`, the
child selections, and the composition together against code that already
passed its own review.

### 3.4 Local energy and the flip patch

```cpp
struct FlipPatch {
    std::vector<int> faces;    // faces incident to i, j, k or l, deduplicated
    std::vector<int> vertices; // union of their control nets
};
struct PatchEval {
    double eBend, eSpring, area, volume;    // sums over the patch
    // per-face values, to write back on accept
    std::vector<double> faceBend, faceArea, faceVolume;
};
PatchEval evaluate_patch(const FlipPatch &p) const;   // no side effects
```

`evaluate_patch` calls `element_energy_force_patch_pod()` and
`element_area_volume_pod()` on each face's current rows and coordinates, as the
main loop does, but only for the listed faces, into thread-local buffers. It
is `const`: the flip sweep evaluates the patch, flips, rebuilds the affected
rows, evaluates again, and unflips on rejection. Evaluating the *old* state
fresh rather than reading `face.energy` costs one extra local evaluation per
attempt and removes any dependence on the stored values being current — which
after a previous accepted flip in the same sweep they would not be unless
written back carefully. It is also what makes the ΔE test in WP2 clean.

```text
ΔE = (E'_bend - E_bend) + (E'_spring - E_spring)
   + (uSurf / 2 area0) ΔA (ΔA + 2 (A - area0))
   + (uVol  / 2 vol0 ) ΔV (ΔV + 2 (V - vol0))                                (4)
```

with `ΔA = A'_patch - A_patch`, `ΔV` likewise, and `A`, `V` the running totals.
The spring term needs only the two edges `(i,j)` and `(k,l)`.

### 3.5 Acceptance and scheduling

```text
per flip sweep (every edgeFlipInterval steps, after Compute_Energy_And_Force):
    λ  = edgeFlipAttemptRate · timeStep · edgeFlipInterval · N_flippable
    n  ~ Poisson(λ)                                       // counter-based RNG, key (seed, "flip", iteration)
    repeat n times:
        e  ~ Uniform(flippable edges)
        if !dynamic_admissible(e): continue                // counts as an attempt, like a rejected move
        old = evaluate_patch(patch(e)); flip_edge(e); rebuild_patches(patch(e)); new = evaluate_patch(...)
        ΔE by (4);  accept with min(1, exp(-ΔE / kT))
        on accept: param.area += ΔA; param.vol += ΔV; write back per-face values; log
        on reject: flip_edge(e); rebuild_patches(...)      // inverse; see 3.2 on face labels
    if any accepted: Compute_Energy_And_Force()            // forces consistent with the new topology (phase 1)
```

`kT` is `param.KBT`, the same constant the Brownian step uses. Flips are
evaluated sequentially — the TriMem argument in §1.1 — and a sweep is `O(λ)`
local evaluations, independent of mesh size except through `λ`.

The full force recompute after an accepted sweep is the simple, correct choice
for phase 1 and costs at most one extra force evaluation per step. The
incremental alternative — subtract the old patch's force contributions, add the
new, and rescale every face's area force by the changed `(A - area0)` — is a
documented optimization for later, not a correctness requirement.

### 3.6 Invalidation

`Mesh::topologyVersion` (an integer, incremented by `flip_edge`) replaces the
count-based checks: `ensure_device_layout()`, `CudaForceBackend::upload_topology()`
and the dense/sparse conversion matrices each remember the version they were
built for and rebuild when it differs. The GPU path rebuilds and re-uploads
the layout once per step in which a flip was accepted; it never sees a
half-flipped mesh because the sweep runs on the host between force
evaluations.

### 3.7 Fluid-mode dynamics (WP4)

Three changes to the Brownian step, all behind parameters that default to the
current behaviour so the existing spectrum workload is bit-identical with the
flags off:

1. **In-plane motion.** `inPlaneDynamicsEnabled = true` removes the `j < 2`
   zeroing. Ghosts and periodic duplicates stay where they are.
2. **Edge spring in the drive.** `edgeSpringEnabled` adds equation (3) to the
   energy and its force to `nodalForce`. Implemented over the edge table, so
   each edge is counted once.
3. **Valence-aware, sparse conversion.** `assign_mesh2surface()` writes the
   Loop limit mask for the actual valence, `(1/2, 1/(2N))` (Warren & Weimer
   2001), as a sparse
   row structure rebuilt for the four changed rows after a flip. `M` is
   `½(I + W)` with `W` the random-walk matrix of the triangulation, so its
   eigenvalues are `½(1 + μ)` with `μ ∈ (-1, 1]`; `μ = -1` needs a bipartite
   graph and a triangulation never is one, and on the hexagonal lattice
   `μ ≥ -1/2`, giving eigenvalues in `[1/4, 1]`. `C = M⁻¹ S` and
   `F_S = M⁻ᵀ F_C` are then solved by a damped Richardson iteration
   (`ω = 8/5` on the lattice, contraction `0.6` per sweep), warm-started from
   the previous step: about 40 sparse sweeps to `1e-10`, a few milliseconds,
   no `O(N²)` storage and no inverse to invalidate. A Krylov solver
   (BiCGSTAB) is the fallback if a distorted fluid mesh pushes `μ` toward
   `-1` and slows the fixed-point iteration; the residual is checked either
   way. Using `M⁻ᵀ`
   rather than `M⁻¹` for the force is the correction the tree already asked
   for. The dense path stays available (`surfaceSolver = dense | iterative`)
   for the bit-identity gate.

### 3.8 Boundaries

**Periodic flat sheet.** Ghost rings are frozen and the fourth ring is slaved by
lattice position, not by connectivity (§2.4). Phase 1 therefore allows flips
only where the whole flip patch is real, non-slaved and interior, which leaves
a band four vertices wide at each edge of the tile solid. The spectrum
analysis already restricts itself to the interior tile, and the untracked
`analysis/membrane_resample.py` in the working tree is the resampling step a
moving-vertex membrane needs. A properly periodic fluid sheet — ghosts as
wrapped images of real vertices, flips mirrored onto their images — is a
boundary-condition refactor independent of this plan and is listed under §6.

**Closed surfaces** have no boundary and no ghosts; every edge is a candidate.
The icosphere workload that `irregular_patch_results.md` §8 already asks for is
the natural first fluid run.

### 3.9 Output and restart

- `face.csv` becomes `face_<iteration>.csv` written with every trajectory frame
  in which `topologyVersion` changed since the last frame, so a frame can
  always be paired with its connectivity.
- `EdgeFlips.csv`: one line per attempt — iteration, edge, `i j k l`, `ΔE`,
  accepted — the diagnostic every DTS paper reports and the input to the
  fluidity measurements in WP6.
- The restart checkpoint gains a `faces` block and the reader refuses a file
  whose face count or `topologyVersion` does not match (dynamics does not
  checkpoint today; this is for `Run_flat` and for parity).

### 3.10 Parameters

| Name | Default | Meaning |
| ---- | ------- | ------- |
| `edgeFlipEnabled` | `false` | run the flip sweep |
| `edgeFlipAttemptRate` | `0.5` | `ν`, attempts per edge per µs, equation (1) |
| `edgeFlipInterval` | `1` | steps between sweeps; `λ` scales with it |
| `edgeFlipMinValence` / `MaxValence` | `4` / `8` | clamped to the row-table range |
| `inPlaneDynamicsEnabled` | `false` | lift the `x, y` zeroing |
| `edgeSpringEnabled` / `edgeSpringConstant` / `edgeSpringRestLength` | `false` / `kCurv` / `lFace` | the mesh-quality term |
| `edgeTetherShape` | `flat` | `flat` (WP6) or `harmonic` (equation (3), WP4) |
| `edgeTetherMinRatio` / `edgeTetherMaxRatio` | `0.6` / `1.8` | the flat tether's allowed range, in multiples of `l0`; the upper bound must exceed `√3` or a flip of an equilateral rhombus is forbidden |
| `surfaceSolver` | `dense` | `dense` (today) or `iterative` (§3.7) |
| `irregularPatchDepthScale` | `1.0` | multiplies the recommended depths, for the cost study in §6 |

Every new key touches three places — `Param`, `input.cpp`, `Parameters.cpp` —
and the `test_io` round-trip.

---

## 4. Work packages

Ordered so that each one lands with a test that would fail without it, and so
that the flip move is the *last* thing wired in. WP0–WP3 and the host side of
WP5 are the target for this conversation; WP4 and WP6 are what a fluid run
needs on top.

### WP0 — Edge table, flip primitive, admission predicate — **landed**

`Mesh::edges` built from faces; `flip_edge()`; `classify_face()` refactored out
of `set_one_ring_vertices_sorted()` with the setup path calling it;
`Face::adjacentFaces` also populated by `setup_from_vertices_faces()`;
`topologyVersion`, and the caches that were checking counts now checking it.

> Gate: on the bowl grid fixture (`test_device_mesh_layout.cpp:39-63`) every
> interior edge flips and unflips to a bit-identical mesh; Euler
> characteristic, manifoldness (every edge in ≤ 2 faces, every vertex fan
> closed) and consistent winding hold after each flip; `classify_face()` agrees
> with the current throw/no-throw behaviour on every existing fixture.
> `data/example` bit-identical.

**Result: met.** 17 tests in `tests/test_edge_flip.cpp`; the whole suite is 107
passing with the one pre-existing `ThermalFluctuationTest` failure unchanged.
The shipped periodic workload is byte-identical across all eight output files
and stdout, serial build, against a worktree at `1fdffbd`.

Four things the plan did not anticipate, all now pinned by tests:

- **A double flip exchanges the two face labels** (§3.2). Intrinsic; guarded by
  refusing to flip an edge whose faces differ in spontaneous curvature, which
  makes the exchange unobservable rather than merely unlikely to matter.
- **`set_adjacent_faces_of_faces()` was never called on the imported-mesh
  path**, so `Face::adjacentFaces` was empty there. Invisible because its only
  consumer, `sort_vertices_on_faces()`, is not called on that path either.
  Now built, and maintained incrementally by the flip with a test against a
  full rebuild.
- **`find_opposite_node_index()` printed unconditionally** from outside an
  empty `if (param.VERBOSE_MODE) {}` block. Once per boundary face at setup,
  which is why nobody noticed; a flip sweep rebuilds one-rings thousands of
  times per run and would have buried the log.
- **`ensure_device_layout()` and the CUDA topology upload could not see a
  flip.** Both keyed on face and vertex counts, which a flip conserves exactly.
  `topologyVersion` and a new `CudaForceBackend::invalidate_topology()` close
  it.

### WP1 — Multi-extraordinary patches by one local subdivision — **landed**

> Gate: the refinement identity of §3.3 on faces with two and three
> extraordinary corners; a 6/6/6 face pushed through the local-subdivision path
> agrees with the direct kernel to quadrature accuracy; force–energy
> conjugacy; rigid-motion invariance; regular workload bit-identical.

**Result: met.** 10 tests in `tests/test_multi_extraordinary.cpp`; the suite is
118 passing with the same one pre-existing failure. The shipped periodic
workload is byte-identical against the `1fdffbd` baseline across all eight
output files and stdout.

**The design changed, and for the better.** The plan proposed *composing the
rows*: `R[Na,d,c,q] · P_c · A_loc`, cached per valence triple. That works but
stores thousands of `7 × K` blocks per triple. The equivalent formulation is to
stop one step earlier and keep only the **prolongation matrix** of each child,
`M` of shape (child width) × K, then

```text
Xc = M · X                      the child's control net
Ec, fc = existing_kernel(Xc)    regular, or Stam at the child's own valence
E = Σ Ec,   f = Σ Mᵀ fc         the chain rule, nothing more
```

Nothing about the kernels changes — a corner child is exactly the kind of patch
the tree already evaluates, just handed a linear image of the parent's control
net. The table holds four matrices of at most 14 × 18 per triple: all 125
triples in `[4,8]³` together are well under 100 kB, against tens of megabytes
for composed rows. It is also easier to check, because `M` is a Loop mask and
can be compared against one written out by hand.

**What the identity test found.** `refine_loop_once()` — WP7 of the previous
work package, the optional `isPreRefinementEnabled` pass — never moved a single
old vertex. It counted incident faces correctly and then halved the count, so
`is_interior()` read `nFaces / 2 == valence`, false for every interior vertex.
Every vertex took the boundary branch, found no boundary neighbours, fell
through to "pin it", and kept its original position. Loop's even-point rule
never ran.

That is not cosmetic. A refinement that moves the new edge points but leaves
the old vertices alone **describes a different limit surface** — it changes the
geometry rather than only the discretization, which is the one thing
refinement must not do. Every existing pre-refinement test passed throughout,
because they checked face and vertex counts, extraordinary-vertex isolation and
the positivity of the refined volume. None compared the refined surface to the
one it refines. `PreRefinementTest.RefinementApproachesTheSameLimitSurface` now
does, and that is the test whose name always claimed it did.

**Three tests changed meaning**, and they are the ones that recorded the
limitation being removed: an octahedron (all 4/4/4) and an icosahedron (all
5/5/5) were rejected at setup and are evaluated now. Each was rewritten to
assert the new behaviour, and a new test keeps what they were really
protecting — that a face with no usable patch is rejected loudly rather than
silently carrying zero energy — using a tetrahedron, whose valence-3 corners
are genuinely outside the supported range.

The icosahedron is the closed-vesicle fixture `irregular_patch_results.md` §8
asked for and could not run.

**Admission gained one condition.** The construction assumes the face's
one-ring is embedded — that its `K = N0 + N1 + N2 - 6` control points are
distinct. On a closed surface too small to hold it they are not, and such a
face is rejected with that reason rather than mis-evaluated. An octahedron
turns out to sit exactly on the boundary of this: `K = 6` and it has six
vertices, so its one-ring *is* embedded and it evaluates.

**The GPU path refuses rather than guesses.** `DeviceMeshLayout` inferred a
face's patch from the width of its one-ring, which is now ambiguous — a 6/5/7
face and a regular one are both 12 wide. It reads `Face::patchKind` instead and
throws on a multi-extraordinary face, naming the valences and saying to use
`forceBackend = cpu`. A silently wrong force field is the one outcome worth
ruling out, and the device side belongs in WP5. The device's own three-value
`slimed::PatchKind` was renamed `DevicePatchKind` in the same pass, because two
enums of that name in scope together is a trap.

**Cost, measured.** A 30×30 bowl grid, 1800 faces, serial, one force
evaluation:

| edges flipped | regular | single | multi | ms/eval |
| --- | --- | --- | --- | --- |
| 0% | 1568 | 0 | 0 | 2.6 |
| 5% | 430 | 457 | 681 | 50.8 |
| 15% | 49 | 253 | 1266 | 86.5 |
| 30% | 22 | 217 | 1329 | 89.7 |
| 50% | 21 | 192 | 1355 | 93.0 |

**A fluid mesh costs about 36× an all-regular one**, and it reaches that
plateau by 15% of edges flipped — the valence distribution equilibrates fast,
so there is no gentle regime to sit in. This confirms §6 risk 1 with a number
rather than an estimate, and it makes `irregularPatchDepthScale` (WP5) and the
GPU backend load-bearing rather than optional. The depths were chosen for a
1e-4 relative bending tail, which is far below the thermal noise a Brownian run
lives in; WP6 item 5 is now the measurement that matters most.

### WP2 — Local patch evaluation — **landed**

`FaceSubsetEnergy`, `evaluate_face_subset()`, `EdgeFlipDelta` and
`evaluate_edge_flip()` by equation (4), in
`src/energy_force/Local_patch_energy.cpp`.

> Gate: for random interior flips, `ΔE` from the local evaluation equals
> `E_total(after) - E_total(before)` from two full `Compute_Energy_And_Force()`
> calls to `1e-9` relative, with `uSurf` and `uVol` nonzero so the global terms
> are exercised.

**Result: met**, at `1e-8` relative over 25 flips. 6 tests in
`tests/test_local_patch_energy.cpp`; the suite is 124 passing with the same one
pre-existing failure, and the shipped workload stays byte-identical.

**Cost, measured.** A bowl grid with 20% of its edges flipped, serial:

| faces | full evaluation | local trial | speedup |
| --- | --- | --- | --- |
| 800 | 41.8 ms | 1.77 ms | 24× |
| 1800 | 96.0 ms | 1.72 ms | 56× |
| 3872 | 208.5 ms | 1.69 ms | 123× |

The trial cost is **flat in mesh size** — about 1.7 ms, set by the eighteen
faces in the flip patch and nothing else — while a full evaluation is linear.
That is the property the sweep needs: a bigger membrane offers more edges to
flip, not more work per flip. The existing single-vertex Metropolis move in
`Energy_minimization.cpp` takes the other route, a full evaluation per trial,
and this is the measurement of what that costs.

Three notes on what the implementation settled:

- **The fixture is a closed icosphere, not a sheet.** A volume constraint on an
  open surface is refused at setup, and rightly: the signed volume of an open
  sheet is not even independent of where the origin sits. The volume term is
  one of the two that cannot be differenced face by face, so it has to be in
  the fixture, and `area0` and `vol0` are set well away from the mesh's own
  values so the cross term `2 (X − X0) ΔX` carries real weight. A version that
  dropped it would fail this gate and nothing else.
- **The regularization energy is now written twice**, once with its force in
  `energy_force_regularization()` and once alone. A test pins them against each
  other face by face, which is what keeps a change to one from silently missing
  the other.
- **The subset evaluator calls the full kernel and discards the forces.** An
  energy-only kernel would be a second implementation of the integrand, and the
  entire value of this routine is that it agrees with the whole-mesh pass
  exactly.

### WP3 — Metropolis sweep and Poisson schedule — **landed**

`Mesh::edge_flip_sweep()` in `src/mesh/Edge_flip_sweep.cpp`, the counter-based
RNG lifted into `include/Counter_rng.hpp`, `EdgeFlipRecord` and the
`EdgeFlips.csv` writer, and the parameters of §3.10 (`edgeFlipEnabled`,
`edgeFlipAttemptRate`, `edgeFlipInterval`, `edgeFlipMinValence/MaxValence`).

> Gate: **two-state test** — freeze all edges but one, run the sweep at fixed
> coordinates, and check that the occupancy ratio of the two triangulations
> equals `exp(-ΔE/kT)` at two temperatures. Poisson gate: the attempt count has
> mean and variance `λ`. Reproducibility: same seed, same flips.

**Result: met.** The chain is reconstructed from the sweep's own log — each
record says whether that attempt was accepted, and an accepted flip toggles the
state — so what is measured is what the sweep did, random numbers and
acceptance rule included, rather than the energy re-derived. At `ΔE/kT` of 1
and 2 the measured ratios sit within counting error of `exp(-1)` and `exp(-2)`.
10 tests in `tests/test_edge_flip_sweep.cpp`; the suite is 134 passing with the
same one pre-existing failure, and the shipped workload stays byte-identical.

**Two findings, both about the same failure mode: a face with no patch.**

- **The one-ring walk was wrong on an irregular mesh.** It asked "what is the
  corner across this edge from that one?" and answered by intersecting the two
  vertices' neighbour lists, taking a common neighbour that was not the
  excluded one. That is correct only when they share exactly two — the corners
  opposite their shared edge — which is what a near-regular mesh gives and what
  every mesh this tree built was. Flips break it by construction: adjacent
  vertices start sharing a third neighbour that forms no face with the edge
  between them, and the walk returned whichever candidate it saw last. The
  symptom was fans that would not close and faces left with no control net.
  The edge table answers the question exactly, since an edge of a two-manifold
  has two incident faces and their third corners are the only candidates there
  have ever been. Measured: 8 such faces after a sweep before the fix, none
  after.

- **A flip that would cost a face its patch is now refused.** Not a defect but
  a policy, and the reasoning matters: a face whose one-ring cannot be built
  carries *no energy*, and zero is the lowest energy there is. A chain allowed
  to reach such a configuration would be actively drawn into it — the
  Hamiltonian would develop a hole and the membrane would tear along it.
  `evaluate_edge_flip()` checks after the trial flip that no face in the patch
  lost its control net, and refuses the move if one did. The guard is what
  keeps a fluid mesh evaluable indefinitely: driving one with every flip the
  admission test allows and no energy at all, it stays fully evaluable, where
  the unguarded primitive reaches faces with no patch within a few hundred
  flips.

**The measure discrepancy, measured.** §1.4(c) noted that the flip samples
`exp(-E) dC` while the Brownian step samples `exp(-E) dS`, differing by
`|det M_T|`. Estimated there at 0.05–0.08 in log weight from the
`-Σ ln(2 N_v)` term alone. Measured over single flips on an icosphere:

```
largest |ln det M' - ln det M| over one flip: 0.00275   (weight ratio 1.0027)
```

An order of magnitude smaller than the estimate, because the estimate looked at
only one of two terms: `ln det M = -n ln 2 - Σ ln N_v + ln det(D + A)`, and the
two valence-dependent pieces very nearly cancel. A 0.3% bias on the relative
weight of triangulations, and none at all on the geometry sampled at fixed
connectivity. The Jacobian correction stays a documented option and is not
worth taking.

**Acceptance and cost, measured.** A 320-face icosphere with roughly 5 nm
edges, `kc = 83.4 pN·nm` (20 kT, the bilayer value of Rawicz et al. 2000),
`kT = 4.17 pN·nm`, `dt = 1 ns`, serial:

| ν (per edge per µs) | λ per sweep | acceptance | ms per sweep |
| --- | --- | --- | --- |
| 0.1 | 0.048 | 26% | 0.04 |
| 0.5 | 0.24 | 30% | 0.28 |
| 1.0 | 0.48 | 25% | 0.58 |
| 5.0 | 2.4 | 25% | 3.0 |

Acceptance is flat in the rate at 25–30%, which is a healthy Metropolis move —
far above TriMem's 0.17% (Siggel et al. 2022; §1.1), because their tether
potential is a stiff penalty where SLIMED's mesh-quality term is soft. WP4's
edge spring will lower it, and that is the number to watch when it lands. At
the default `ν = 0.5` the sweep costs 0.28 ms against tens of milliseconds for
the force evaluation it sits beside, so fluidity is not what makes a fluid run
expensive — the irregular faces it creates are (WP1).

### WP4 — Fluid-mode dynamics — **landed**

The three items of §3.7: in-plane motion behind `inPlaneDynamicsEnabled`, the
edge spring in `src/energy_force/Edge_spring.cpp`, and the sparse valence-aware
conversion in `src/dynamics/Surface_solver.cpp` behind `surfaceSolver`.

> Gate: with all flags off, bit-identical (serial); with
> `surfaceSolver = iterative` the trajectory agrees with dense to the solver
> tolerance; `M⁻ᵀ` versus `M⁻¹` on a mesh with mixed valences differs and the
> transposed one is right; the edge-spring force passes finite-difference
> conjugacy.

**Result: met.** 8 tests in `tests/test_fluid_dynamics.cpp`; the suite is 142
passing with the same one pre-existing failure. The shipped workload is
byte-identical with the flags off, and a 40-step run with
`surfaceSolver = iterative` agrees with the dense path to `1.7e-13` on
coordinates of order 50.

**The solver is conjugate gradients, not the damped Richardson iteration §3.7
proposed.** Richardson needs a spectral bound assumed in advance, and the one
the plan quoted (`μ ≥ -1/2`) is a property of the hexagonal lattice, not of a
fluid membrane. There is a better observation available. Over the free
vertices the mask factors as

```text
    M = D⁻¹ K,        K = ½ (D + A)
```

with `D` the valences and `A` the adjacency, so `K` is symmetric by
construction. It is also positive definite: `D + A` is the signless Laplacian,
positive semidefinite for any graph and singular only on a bipartite one
(Cvetković, Rowlinson & Simić 2007), and a triangulation has triangles. Both
directions then reduce to the same system:

```text
    M C = S      ⟺   K C = D S
    Mᵀ F_S = F_C  ⟺   K y = F_C,  F_S = D y
```

Conjugate gradients on an SPD matrix, with no parameter to tune, no spectral
bound to assume, and a residual that says whether it worked. Warm-starting from
the previous step's control net is what makes it cheap.

**Pinned vertices are eliminated rather than solved for.** A vertex whose faces
are all ghost has no limit surface of its own and carries the identity row.
Those rows are what made `M` asymmetric, which is why the tree used `M⁻¹` in
place of `M⁻ᵀ` and reported the discrepancy at startup rather than fixing it.
They are not degrees of freedom — the Brownian step skips them — so their known
values move to the right-hand side and what remains is a principal submatrix of
`K`, still symmetric and still positive definite. There is no asymmetry left to
approximate around.

**The dense path was worse than anyone had measured.** Its setup inverts an
N × N matrix and its step multiplies by one, serial in both cases:

| sheet | vertices | dense setup | sparse setup | dense step | sparse step | step speedup |
| --- | --- | --- | --- | --- | --- | --- |
| 60 nm | 195 | 5.9 ms | 0.02 ms | 0.09 ms | 0.016 ms | 5× |
| 120 nm | 725 | 199 ms | 0.06 ms | 1.23 ms | 0.076 ms | 16× |
| 200 nm | 1,927 | 4.3 s | 0.12 ms | 10.5 ms | 0.248 ms | 42× |
| 300 nm | 4,331 | **54.7 s** | 0.24 ms | **139 ms** | 0.590 ms | **236×** |

Fifty-five seconds of setup and 139 ms per step on a 300 nm sheet, against a
quarter of a millisecond and half a millisecond. This is the bottleneck the
build notes recorded as "the dense surface2mesh multiply dominates and does not
thread"; it is `O(N³)` at setup and `O(N²)` per step, and it is gone.

`dense` stays the default so that existing runs reproduce, including its
standing `M⁻¹`-for-`M⁻ᵀ` approximation. It is refused outright when
`edgeFlipEnabled` is set: a stored inverse cannot survive a flip, and a run
configured that way would either be unusably slow or quietly keep using an
inverse that no longer describes its mesh.

**The valence-6 mask is corrected on both paths.** `assign_mesh2surface()` wrote
`1/12` on every neighbour whatever the valence was. Right on a regular mesh,
wrong everywhere else, and a fluid membrane is mostly not valence 6. Every
vertex on the workloads this tree has run is at valence 6 where it matters, so
the correction changes nothing there — which the bit-identity gate confirms.

**The edge spring replaces the reference-length regularization rather than
adding to it.** The old term measures a face's edges against the same face's
edges in `coordRef`, which is a solid's memory of where it started; a flipped
edge has no reference length, and `coordRef` would hand it the distance between
two vertices that were not joined. A test makes the distinction concrete: with
`coordRef` equal to the current coordinates the old term is exactly zero, and
the spring is not, because its rest length is a parameter rather than a memory.
The spring's force joins the Brownian drive, which the old term's never did —
with the in-plane displacement zeroed it had nothing to act on.

### WP5 — Wiring, output, GPU — **landed**

Sweep call in `Run_dynamics_flat.cpp` after the end-of-step force evaluation;
version-based invalidation (§3.6); `face_<iteration>.csv`; checkpoint faces
block; `irregularPatchDepthScale`.

> Gate: a short flat run with flips on writes paired coordinate/face frames
> that reload into a consistent mesh; the layout rebuild fires exactly on
> steps with accepted flips.

**Result: met, with one part of the gate found to be unreachable and replaced.**
21 tests in `tests/test_fluid_run_io.cpp`; the suite is 163, 162 passing and 1
skipped (CUDA, no local device). The shipped periodic workload is byte-identical
to a serial build at `1fdffbd` across all 8 outputs with the new flags off. A
3000-step fluid run on a 100 nm sheet (525 vertices, 960 faces,
`edgeSpringConstant = 20`, `meshpointOutputInterval = 10`) accepted 17 flips and
wrote 15 face frames: none duplicated, all two-manifold, a baseline at iteration
0, so every coordinate frame pairs with a connectivity.

**The device-layout half of the gate cannot be met and was not faked.**
`DeviceMeshLayout` refuses a face with more than one extraordinary corner, and a
flip leaves extraordinary corners at *both* ends of the new edge — so a flipped
mesh is exactly what it cannot build. The GPU backend therefore cannot run a
fluid membrane at all until the device kernel gains WP1's multi-extraordinary
path. `DynamicMesh::setup_flat()` now refuses `edgeFlipEnabled` with
`forceBackend = gpu` before the first step, rather than letting the run flip and
then throw out of the layout builder with a mesh already changed. What is tested
instead is the invalidation predicate itself, on a mesh the layout can build: a
version change rebuilds it and an unchanged version does not.

#### A trial flip was moving the invalidation signal

`evaluate_edge_flip()` flips the edge, measures, and flips it back, and each of
those bumped `topologyVersion`. A *rejected* attempt therefore left the mesh
exactly as it was and the version two ahead. Nothing became wrong —
over-invalidation is safe — but the version is what says whether a cache is
stale and whether a frame needs its connectivity written beside it, and on a
fluid run most attempts are rejected. Measured before the fix: a 500-step run
with 121 attempts and **zero** acceptances wrote **112** identical face frames
and rebuilt the sparse limit mask 121 times. The trial now restores the version
with the mesh; the same run writes one baseline frame and nothing else.

#### Two defects found in the paths this wiring runs through

- **`set_adjacent_faces_of_vertices_sorted()` is only valid for the pristine
  grid.** It reads a vertex's six adjacent faces off the generated flat sheet's
  face numbering and *erases* anything not among them. After a flip a vertex can
  have seven, and the seventh has no grid index, so it was silently dropped —
  the vertex then claimed fewer faces than named it. `flip_edge()` maintains
  these as sets rather than fans for exactly this reason. Split into
  `set_adjacent_faces_of_vertices_unsorted()`, which the connectivity restore
  uses. (`nFaceX` defaults to `-1`, so on an imported mesh the sorting loop
  never ran and the bug was invisible there.)
- **`Model::stepSize` was uninitialized**, and the restart checkpoint writes it.
  A garbage double is very often subnormal, and libc++'s `operator>>` sets
  failbit on a subnormal even though it parses the value correctly — so such a
  checkpoint could not be read back at all. Initialized to `0.0`. The residual
  risk is unfixed: any genuinely subnormal double anywhere in a checkpoint still
  makes it unreadable on this platform.

#### Output and restart

- `<prefix>face_<iteration>.csv`, written beside a trajectory frame whose
  connectivity differs from the last one written, in the same three-column
  layout as `face.csv`. Gated on `edgeFlipEnabled`, so a run with flips off
  produces exactly the files it always did.
- `<prefix>EdgeFlips.csv`, one line per attempt, buffered and flushed in blocks.
- The checkpoint is `SLIMED_RESTART_V2` and carries a `faces` block. The reader
  still accepts V1. `Mesh::restore_face_connectivity()` puts the connectivity
  back and rebuilds everything derived from it; a restore that would not produce
  a two-manifold is refused with the mesh untouched, including when the rebuild
  *throws* rather than returning a verdict.
- A run whose total energy goes non-finite stops at the first such step. The
  3000-step run that found the tether problem below wrote NaN into every row
  after step 1535 and spent four minutes doing it; it now stops in 69 s with no
  NaN written.

#### `irregularPatchDepthScale`

Multiplies every valence's recommended depth, rounded, floored at one level.
Under `PerValence` the *built* depth follows it too, so a scale below 1 is
cheaper to build as well as to evaluate. `1.0` reproduces the unscaled table
exactly. At `0.5` the built depth drops from 12 to 6 and the icosphere's total
energy moves by under 5% — a coarser answer, not a different one. `Uniform`
ignores it, so the convergence study can still sweep the depth itself.

#### The harmonic tether has no usable stiffness — for WP6

The edge spring of §3.7 is harmonic at rest length `l0`. A flip on a rhombus of
two equilateral triangles replaces the short diagonal by the long one, so the
move has to climb

```text
    ΔE = (k / 2) (√3 − 1)² l0²
```

and Metropolis accepts it with probability `exp(−ΔE/kT)`. Two requirements pull
`k` in opposite directions:

| requirement | condition | at `l0 = 5 nm`, `kT = 4.17` |
| --- | --- | --- |
| flips are possible (barrier ≲ 5 kT) | `k ≤ 10 kT / ((√3−1)² l0²)` | `k ≲ 3.1 pN/nm` |
| the triangulation survives (bond fluctuation ≤ 0.1 `l0`) | `k ≥ 100 kT / l0²` | `k ≳ 16.7 pN/nm` |

**The window is empty, by a factor of about five.** Measured on the 100 nm
sheet, 3000 steps:

| `edgeSpringConstant` | barrier | accepted / attempted | outcome |
| --- | --- | --- | --- |
| `83.4` (the default, `kCurv`) | 134 kT | 0 / 121 | frozen solid |
| `20.0` | 32 kT | 17 / 801 (2.1%) | stable, barely fluid |
| `1.0` | 1.6 kT | 37 / 415 (8.9%) | **diverges at step 1535** |

The `k = 1.0` divergence begins in the dynamics, not in the flip move —
`E_curvature` reaches `2.1e7` at step 1535, before the first large flip `ΔE` at
step 1540 — and it needs the flips to trigger it: the same configuration with
`edgeFlipEnabled = false` is stable at `E = 2001` after 3000 steps. A tether too
weak to prevent a degenerate triangle leaves the energy unbounded below in a
direction the flip move can reach, and the move then walks down it.

This is why every dynamically triangulated surface model in §7 uses a
**flat-bottomed** tether — zero energy for `l ∈ [l_min, l_max]`, a wall outside
— rather than a spring (Gompper & Kroll 2004; Ramakrishnan, Sunil Kumar &
Ipsen 2010; Pezeshkian & Ipsen 2024). Inside the range a flip costs nothing, so
the barrier and the constraint stop competing. Implementing that is the first
item of WP6.

Until then, `DynamicMesh::setup_flat()` reports the barrier in kT at startup and
warns above 10 kT, and warns separately when flips are enabled with the spring
off — that configuration leaves the *reference-length* regularization in force,
which measures an edge a flip just created against a length it never had, and
the sweep samples it: measured at about −1200 pN·nm per accepted flip of pure
artifact.

### WP6 — Validation of fluidity and physics — **landed**

**First item, ahead of the measurements below: replace the harmonic edge spring
with a flat-bottomed tether.** WP5 measured the harmonic one to have no usable
stiffness — see the table there. Every measurement in this package is a
measurement of the flip move's behaviour, and with the current tether the move
either never fires or drives the run to divergence.

**Result: the gate is met.** All numbers are in
[`fluidity_results.md`](fluidity_results.md); 12 tests in
`tests/test_fluidity.cpp` and the analysis module `analysis/fluidity.py`. The
suite is 175, 174 passing and 1 skipped (CUDA, no local device), and the
shipped periodic workload is still byte-identical with the new flags off.

The gate, item (3): on the 60 nm sheet over 60000 steps, the in-plane MSD with
flips off saturates at 7.7 nm² (growth exponent 0.097) while with flips it
reaches 23.1 nm² and is still climbing (exponent 0.482). The two are
indistinguishable out to lag 800, cross at about lag 1300, and reach a ratio of
3.0 by lag 30000.

Item (4), the fluctuation spectrum with flips on, is measured separately in
`analysis/membrane_fluctuation_fluid_cpu.ipynb`, once the throughput
correction below made a 400 000-step fluid run an hour rather than the two
days first estimated.

Three findings worth carrying forward:

- **The sweep and the dynamics were sampling different Hamiltonians.**
  `evaluate_face_subset()` went on differencing the reference-length
  regularization after the tether landed, so every accepted flip reported about
  −787 pN·nm while the mesh's energy climbed. Fixed; acceptance went from 17%
  to 40% and the mean accepted ΔE from −787 to +0.3 pN·nm.
- **The ghost band is not the membrane.** Every long fluid run diverged, and
  the tether energy climbed beforehand in every configuration — all of it in
  the ghost band, whose lattice edges join positions that stopped being
  neighbours once the interior mixed, and whose tether force reached the
  interior through `M⁻ᵀ` before the duplicates were overwritten. Fixed by
  excluding edges with two copied endpoints from the tether. The earlier
  claim here that the tether range had to be narrow was the same artifact
  measured a different way; the interior is stationary at either range, and
  `[0.95, 1.75]` is kept for mesh quality.
- **The dynamics is not fluid-safe.** With the exclusion in place the sheet
  still diverges after 10⁴–10⁵ steps, and the full logs put the first fold
  in the interior every time. The flat tether bounds edge lengths and not
  triangle shape; a fluid mesh carries slivers 0.1 nm tall at 1° continuously,
  and an explicit Brownian step — unlike a Metropolis move — eventually walks
  one through a fold. A triangle-shape term, or a rejecting step, is the next
  work package. `docs/fluidity_results.md` section 8 has the numbers.
- **Fluidity costs 4.7× in throughput** at `-O3` on one thread (the 21× first
  recorded here was measured with the `Makefile.legacy` binary, which has no
  optimisation flag and is a `-O0` build). It is the irregular-patch cost of
  WP1 rather than the flip move; `irregularPatchDepthScale = 0.5` buys 1.3× of
  it back without moving the acceptance rate, and OpenMP buys 2.3× at 8
  threads because the per-face work is where the time goes.

On the periodic sheet and on an icosphere:

1. acceptance rate and valence histogram versus `ν` (DTS equilibrium meshes
   are roughly 60% valence 6, 20% each 5 and 7, the coordination distributions
   reported by Gompper & Kroll 2004 and Ramakrishnan, Sunil Kumar &
   Radhakrishnan 2015 — a sanity band, not a target);
2. **neighbour survival**: the fraction of initial edges still present after
   time `t`, whose decay time is the microscopic fluidity time scale and the
   thing `ν` should be calibrated by — the bond-flip time scale Noguchi &
   Gompper (2004, 2005) tie to the membrane viscosity `η_mb`;
3. **in-plane MSD** of tagged vertices: saturating at the cage size without
   flips, linear in `t` with them — the classical signature of a fluid, and the
   in-plane diffusion Sadeghi, Weikl & Noé (2018) tie to the flip frequency;
4. the fluctuation spectrum with flips on, through the resampled pipeline:
   the same `kc` from the Helfrich `q⁻⁴` fit (Helfrich 1973), and an effective
   tension that does not grow with `k_S`;
5. throughput per step against §6 risk 1, at depth scale 1.0 and 0.5.

> Gate: (3) and (4). Numbers recorded in a results document, as
> `irregular_patch_results.md` did for its plan.

### WP7 — The triangle-shape term — **landed**

The tether bounds every edge of a fluid mesh and nothing else, and every long
fluid run ended in a folded sliver whose edges were all inside the walls:
interior triangles 0.1 nm tall at one degree in every frame, whose normals
turn through tens of degrees under a single 0.05 nm Brownian kick, whose
limit-surface patches self-intersect, and whose bending force is then not
finite. A Monte Carlo model never takes that step, because the energy rejects
it; an explicit Brownian step has no refusal in it, so the Hamiltonian has to
carry one. (`fluidity_results.md` §8 has the logs and the altitude tables.)

The term bounds the quantity that vanishes in a sliver — the altitude from
each corner to its opposite edge, `h_i = 2A/l_i` — with the tether's own
flat-bottomed wall:

```text
    E_face = (k/2) Σ_{i=1..3} max(0, h₀ − h_i)²,    h₀ = triangleShapeMinAltitudeRatio · lFace
```

Zero for any healthy triangle, so it adds no tension and does nothing inside
the allowed region. The gradient is closed-form (`∇_{p_k} A = ½ e_k × n̂`,
`∇ l_i` the unit edge), and one expression serves the force pass and the flip
trial, which is the WP6 lesson made structural. The floor's ceiling is fixed
by geometry: a flip of an equilateral rhombus makes two triangles of altitude
exactly `lFace/2`, so the ratio must stay below 0.5 or the term forbids the
move the model exists to permit; 0.4 leaves that flip a margin and puts the
wall at 2 nm on a 5 nm mesh, where a Brownian kick turns a normal by under two
degrees. It requires the tether (setup refuses it without) and shares the
regularization slot with it.

Parameters: `triangleShapeEnabled` (`false`), `triangleShapeMinAltitudeRatio`
(`0.4`), `triangleShapeConstant` (`83.4` pN/nm). `tests/test_triangle_shape.cpp`.

**What the first gate run found.** With the altitude floor holding every
triangle above 1.4 nm, the sheet still folded at step 82 000 — not a sliver
but a *flap*: adjacent faces at 177–180°, a face folded flat onto its
neighbour, growing in number from step 20 000 and every one at a valence-8 or
valence-4 vertex. The reason is geometric. A flat vertex of valence N with
legs of `1.1·lFace` (where a fluid run's edges sit) needs opposite edges of
`2.2·lFace·sin(π/N)`: 4.2 nm at valence 8 against a tether wall at 4.75. It
cannot flatten, its surplus angle buckles the neighbourhood, and with the
limit-surface bending energy indifferent to a crease in the control net the
buckle becomes a flap whose limit surface pinches. Three 60 000-step probes:
widening the tether made it worse (short edges let the bending force crush
triangles through the altitude wall); **flips restricted to valences 5–7**
gave zero creases for three quarters of the run and three at the end. That
is now the default, and setup reports whether the lower wall lets a vertex
of the maximum valence flatten.

**The crease wall.** The third term, and the one a dynamically triangulated
surface gets for free from its `Σ (1 − n̂₁·n̂₂)` discretization of the bending
energy (Gompper & Kroll 2004; Helfrich 1973 for the continuum term it
discretizes): for each interior edge, with `c = n̂₁·n̂₂` the cosine between its
faces' normals,

```text
    E_edge = (k/2) max(0, cos θ_max − c)²,    θ_max = 60°, k = 500 pN·nm
```

zero within 60° of coplanar, smooth everywhere (a wall on the angle itself
has a singular gradient at a full fold, which is where it matters), 15 kT at
a right angle and 135 kT at a full fold. Closed-form gradient
(`∇_p c = e_k × (n̂₂ − c n̂₁)/|n₁|` for a corner of face 1), one expression
for the force pass and the flip trial, half of each edge to each of its
faces. Parameters: `creaseWallEnabled` (`false`), `creaseWallAngle`,
`creaseWallConstant`. `tests/test_crease_wall.cpp`.

**Gate met, and what it then measured.** The 400 000-step run with the
altitude floor and valences 5–7 finished — no crease over 90° in its second
half, altitude floor at 1.4–2.0 nm, survival 0.70, MSD exponent 0.75 — and
with the crease wall on, zero creases at every eighth. Its spectrum, though,
is a crumpled membrane's: κ_c 24 pN·nm at slope −6.8, r.m.s. height 9 nm.
At μ_S = 0 nothing sets the control net's length scale but the tether walls;
the mean edge grows to 6.15 nm and the 50% of excess area buckles out of
plane in a box of fixed projected area. The 94 µs run of WP6 had matched the
solid only because its slivers absorbed that excess as length without area.
The comparison at conserved area — μ_S = 250, the fluid `fluid_c` against
the solid `tension` — is the one the notebook now makes by default:
**κ_c fluid/solid = 1.02** on the whole-trajectory fit (0.83 ± 0.19 by
blocks), tension 3.0 against 3.8 pN/nm, no earlier roll-off, reader
converged. What that run does not yet have is fluidity: at conserved area
with valences 5–7 the neighbour survival is 0.93 after 400 µs. Whether the
crease wall makes the full 4–8 range safe there (`fluid_d`) is the next
measurement; `fluidity_results.md` §10–12.

> Gate: a 400 000-step fluid run on the 100 nm sheet that does not diverge —
> every previous one died at 8 700, 27 800 or 94 000 steps — with the
> smallest interior altitude held above the floor's shadow, the bending
> energy stationary, and the flip acceptance intact; then the fluctuation
> spectrum re-measured on the full trajectory.

**Size.** WP0 ≈ 400 lines, WP1 ≈ 700, WP2–3 ≈ 500, WP4 ≈ 400, WP5 ≈ 300, plus
tests of similar size. The first four are the substance; WP1 is the one that
can go wrong quietly, which is why its gate is an exact identity.

---

## 5. Tests that bind

- **Flip involution** and manifold invariants after every flip (WP0).
- **Refinement identity** for multi-extraordinary faces (WP1). Exact, not
  asymptotic.
- **Local `ΔE` equals global `ΔE`** with the global constraint terms on (WP2).
- **Two-state Boltzmann ratio** (WP3). This is the detailed-balance test in the
  only form that is exact for a Metropolis chain.
- **Force–energy conjugacy** on multi-extraordinary faces and on the edge
  spring (WP1, WP4), reusing `test_convergence_study.cpp:210`.
- **Bit-identity** of the regular workload with every new flag off, through
  every WP.
- **MSD linear with flips, saturating without** (WP6).

---

## 6. Risks and open questions

| Risk | Handling |
| ---- | -------- |
| **A fluid mesh is mostly irregular, and irregular faces are expensive.** **Measured at WP1: 36× an all-regular mesh**, reached by the time 15% of edges have flipped. In a fluid steady state almost every face has at least one extraordinary corner and most have several, each costing `3·D` samples per corner instead of 3. | This is the real cost of fluidity and it is not specific to this design — any subdivision membrane that flips pays it. Levers, in order: `irregularPatchDepthScale` (the `1e-4` bending tail the depths were chosen for is far below the thermal noise a Brownian run lives in; WP6 measures what depth the spectrum actually needs); the GPU backend, which was built for exactly this kind of face-parallel load; a higher-order rule on the children so that a given accuracy needs fewer levels. Now a known quantity rather than a risk, but it moves `irregularPatchDepthScale` from a convenience to a requirement. |
| **The measure question** (§1.4c). | Plain Metropolis at fixed `C` is the literature standard (Gompper & Kroll 2004; Ramakrishnan, Sunil Kumar & Ipsen 2010; Siggel et al. 2022); the log-det diagnostic quantifies the discrepancy; the corrected acceptance is a documented option. Geometry sampled at fixed `T` is unaffected either way. |
| **Periodic band stays solid** (§3.8). | Accepted for phase 1; the interior tile is what the analysis measures. Torus-periodic connectivity is a separate plan. |
| **Valence range `[4, 8]`** rejects flips that DTS models would allow. | Measure the histogram; extend the generator to `3` and `9–10` if the rejection rate at the bounds is material. `N = 3` needs Loop's `β = 3/16`. |
| **Scaffolding and insertion faces.** `Particle::faceIndex` and per-face `spontCurvature` are indexed by face. | Faces carrying either are excluded from flips in phase 1; remapping on flip is straightforward but is its own change. |
| **Regularization term ambiguity.** Fluid mode adds a spring while the old regularization stays in the minimizer. | Fluid mode disables `energy_force_regularization()` and uses the spring; the minimizer is untouched. Documented in the parameter comments. |
| **`ThermalFluctuationTest` already failing.** | Fix or quarantine before adding a second Metropolis path, so the suite is green when the two-state test lands. |

### Decisions needed before WP1

- Confirm phase 1 samples at fixed control net with plain Metropolis (§1.4c,
  option 1), with the Jacobian correction as a later option.
- Confirm the edge spring replaces, rather than supplements, the reference-length
  regularization in fluid mode.
- Confirm `[4, 8]` stays the valence range for the first fluid runs.
- Confirm the periodic sheet keeps its solid band in phase 1, with the icosphere
  as the first fully fluid workload.

---

## 7. References

Cited inline where a result, a model or a parameter value is taken from the
source rather than derived here.

### Dynamically triangulated membranes

- Gompper, G. & Kroll, D. M. *Triangulated-surface models of fluctuating
  membranes*, in *Statistical Mechanics of Membranes and Surfaces* (2004).
  The partition function of §1.1, the flip validity checks, the hard-wall
  tether, the `[3, 9]` valence convention, the equilibrium coordination
  distribution of §4 WP6, and the `Σ (1 − n̂₁·n̂₂)` bending discretization
  of WP7.
  https://www.researchgate.net/publication/285604301_Triangulated-surface_models_of_fluctuating_membranes
- Gompper, G. & Kroll, D. M. *Membranes with fluctuating topology: Monte Carlo
  simulations*, PRL 81, 2284 (1998). The topology-changing moves this plan
  excludes. https://ui.adsabs.harvard.edu/abs/1998PhRvL..81.2284G/abstract
- Ramakrishnan, N., Sunil Kumar, P. B. & Ipsen, J. H. *Monte Carlo simulations
  of fluid vesicles with in-plane orientational ordering*, PRE 81, 041922
  (2010). Eq. 19, the symmetric flip proposal of §1.1; the sweep definition;
  the `√3 a_0` tether. https://arxiv.org/abs/1004.4509
- Ramakrishnan, N., Sunil Kumar, P. B. & Radhakrishnan, R. *Monte Carlo
  simulations of fluid vesicles*, J. Phys.: Condens. Matter 27, 273104 (2015).
  Review; the DTS partition function of §1.1 and the valence statistics of WP6.
  https://iopscience.iop.org/article/10.1088/0953-8984/27/27/273104
- Noguchi, H. & Gompper, G. *Fluid vesicles with viscous membranes in shear
  flow*, PRL 93, 258102 (2004). The bond-flip rate as the membrane-viscosity
  knob (§1.2, §1.3, WP6 item 2) and the smooth bond potential of §1.6.
  https://arxiv.org/abs/cond-mat/0404356
- Noguchi, H. & Gompper, G. *Dynamics of fluid vesicles in shear flow: Effect of
  membrane viscosity and thermal fluctuations*, PRE 72, 011901 (2005).
  https://journals.aps.org/pre/abstract/10.1103/PhysRevE.72.011901
- Sadeghi, M., Weikl, T. R. & Noé, F. *Particle-based membrane model for
  mesoscopic simulation of cellular dynamics*, J. Chem. Phys. 148, 044901
  (2018). Flips at frequency `φ` alongside Langevin dynamics,
  `η = η_∞ exp(C_φ/φ)`, and the entropy-production argument of §1.2.
  https://arxiv.org/abs/1710.06907
- Siggel, M. et al. *TriMem: A parallelized hybrid Monte Carlo software for
  efficient simulations of lipid membranes*, J. Chem. Phys. 157, 174801 (2022).
  Eq. 19 (the acceptance `ε`, and the 0.17% WP3 compares against), eq. 14 (the
  continuous tether), Fig. 2 (the flip patch of §1.4a), and the sequential-flip
  argument of §1.1. https://www.biorxiv.org/content/10.1101/2022.05.25.493239v1
- Pezeshkian, W. & Ipsen, J. H. *Mesoscale simulation of biomembranes with
  FreeDTS*, Nat. Commun. 15, 548 (2024). The `[l_dts, 3 l_dts]` flat tether and
  the global-constraint-in-a-local-move form of §1.5.
  https://www.nature.com/articles/s41467-024-44819-w
- *OrganL: Dynamic triangulation of biomembranes using curved elements* (2024).
  Curved elements with dynamic triangulation; the
  `ΔF_Φ = λ/2 ΔΦ (ΔΦ − 2(Φ − Φ_0))` form equation (2) takes.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11213972/
- *PyMembrane: A flexible framework for efficient simulations of elastic and
  liquid membranes* (2023). Elastic and liquid membranes in one framework, the
  pair `edgeFlipEnabled` switches between. https://arxiv.org/abs/2308.12754

### Membrane physics and parameter values

- Helfrich, W. *Elastic properties of lipid bilayers: theory and possible
  experiments*, Z. Naturforsch. C 28, 693–703 (1973). The bending energy whose
  `q⁻⁴` spectrum §1.6 and WP6 item 4 fit against.
  https://doi.org/10.1515/znc-1973-11-1209
- Rawicz, W., Olbrich, K. C., McIntosh, T., Needham, D. & Evans, E. *Effect of
  chain length and unsaturation on elasticity of lipid bilayers*, Biophys. J.
  79, 328–339 (2000). The `κ ≈ 20 kT` bilayer bending rigidity behind
  `kc = 83.4 pN·nm`.
- Almeida, P. F. F. & Vaz, W. L. C. *Lateral diffusion in membranes*, in
  *Handbook of Biological Physics* vol. 1, ch. 6, 305–357 (1995). The
  `10⁻⁸–10⁻⁷ cm²/s` fluid-phase lipid diffusion constants §1.3 calibrates `ν`
  against.

### Subdivision surfaces, meshes and methods

- Loop, C. T. *Smooth subdivision surfaces based on triangles*, M.S. thesis,
  University of Utah (1987). The subdivision scheme itself: the vertex rule
  `β = 3/(8N)` (and `3/16` at `N = 3`) and the `(3/8, 3/8, 1/8, 1/8)` edge rule
  that build `A_loc` in §3.3.
- Stam, J. *Evaluation of Loop subdivision surfaces*, SIGGRAPH 1998 course
  notes. The exact evaluation the row tables implement, and the
  one-subdivision isolation argument that §1.4(b) applies per face.
- Warren, J. & Weimer, H. *Subdivision methods for geometric design: a
  constructive approach*, Morgan Kaufmann (2001). The `(1/2, 1/(2N))` limit
  mask of §1.4(c) and §3.7.
- Botsch, M., Steinberg, S., Bischoff, S. & Kobbelt, L. *OpenMesh — a generic
  and efficient polygon mesh data structure*, OpenSG Symposium (2002). Its
  `flip()` shows the same face-label exchange as §3.2.
  https://www.graphics.rwth-aachen.de/software/openmesh/
- Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H. & Teller,
  E. *Equation of state calculations by fast computing machines*, J. Chem.
  Phys. 21, 1087–1092 (1953). The acceptance rule of §1.1 and §3.5.
  https://doi.org/10.1063/1.1699114
- Cvetković, D., Rowlinson, P. & Simić, S. K. *Signless Laplacians of finite
  graphs*, Linear Algebra Appl. 423, 155–171 (2007). `D + A` is positive
  semidefinite and singular only on a bipartite graph — what makes the WP4
  solve SPD.
