# Monte Carlo Edge Flips for a Fluid Membrane

**Status:** work package 0 landed; 1-6 planned
**Base:** `JohnsonBiophysicsLab/SLIMED @ 1fdffbd`
**Builds on:** [`irregular_patch_results.md`](irregular_patch_results.md) (valence 4–8
row tables), [`fluctuation_spectrum.md`](fluctuation_spectrum.md) (the end-to-end
check this work must keep passing)
**Scope:** in-plane fluidity by Metropolis edge flips, scheduled as a Poisson
process in physical time, interleaved with the existing Brownian dynamics.
Fixed surface topology — no fusion or fission.

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
vertex positions:

```text
Z = Σ_T ∫ Π_v dX_v  exp(-E(X, T) / kT)
```

Two Monte Carlo moves sample this: a vertex displacement at fixed `T`, and a
**link (bond, edge) flip** at fixed `X`. In a flip the shared edge `(i, j)` of
two adjacent triangles `(i, j, k)` and `(i, l, j)` is removed and replaced by
`(k, l)`, giving triangles `(i, l, k)` and `(j, k, l)`. The four valences change
by `-1, -1, +1, +1`; the number of vertices, edges and faces does not. The move
is accepted with the Metropolis probability `min[1, exp(-ΔE/kT)]`. With a
uniformly chosen edge, the proposal is symmetric — the reverse flip is proposed
from the new state with the same probability `1/N_E`, since the flip is an
involution and `N_E` is conserved — so no proposal-ratio correction is needed
(Ramakrishnan, Sunil Kumar & Ipsen 2010, eq. 19; Gompper & Kroll 2004).

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
- **TriMem (2022)** alternates hybrid-MC trajectories with flip sweeps and
  verifies Boltzmann sampling by reproducing the vesicle phase diagram.
- **FreeDTS (2024)** puts `N_T` flip attempts, `N_v` vertex updates and the
  inclusion moves into one MC step, and folds the *global* constraints into
  the local energy change of each move.
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
membrane's in-plane viscosity (Noguchi & Gompper) and the vertex diffusion,
and should be calibrated (WP6) against the physical neighbour-exchange time of
a mesh vertex. A vertex at `lFace = 5 nm` stands for a patch of order a hundred
lipids; with a lipid diffusion constant of `1–10 nm²/µs` the time for a patch
to exchange a neighbour is `l² / 4D ≈ 0.6–6 µs`, so `ν` in the range
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
near-regular mesh. It is the same object TriMem locks for a parallel flip (their
Fig. 2, which notes it is "significantly larger than the patch required for a
flip subject to the Delaunay criterion"). Everything outside the flip patch
has an unchanged control net and an unchanged energy; that is what makes `ΔE`
local.

**(b) Faces with several extraordinary corners are unavoidable.** Stam's
evaluation, and the row tables built on it, need exactly one extraordinary
corner per face. That is achievable at mesh generation, and WP7's global
refinement makes any mesh satisfy it. A flip destroys it immediately:

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
the limit surface is unchanged because Loop subdivision is exactly the map the
limit surface is defined by. The child at corner `v_a` keeps the valence of
`v_a` and gets two valence-6 neighbours; the centre child is regular.

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
`|det M_T|`, which depends on `T`. With Warren's weights the limit mask row of
a valence-`N` vertex is `1/2` on itself and `1/(2N)` on each neighbour, so
`M = ½ D⁻¹ (D + A)` and

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
summed over the patch the exact change is the OrganL / FreeDTS form:

```text
ΔE_A = (uSurf / 2 area0) · ΔA · (ΔA + 2 (A - area0))                      (2)
```

and likewise for volume. This needs the per-face areas and volumes the tree
already stores (`Face::elementArea`, `elementVolume`) and the running totals
`param.area`, `param.vol`, updated after every accepted flip within a sweep so
the next attempt sees the right `A`.

### 1.6 Mesh quality: the tether becomes an edge spring

DTS models keep triangles well-shaped with a tether potential — hard walls in
the Monte Carlo models, a smooth well in the MD ones (Noguchi & Gompper
introduced "a smooth bond-interaction potential, which makes the model
amenable for molecular dynamics"; TriMem's eq. 14 is a continuous version for
the same reason). The Metropolis flip then automatically prefers the
Delaunay-like diagonal, because the other one is longer and costs tether
energy.

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
spectrum fit in `membrane_spectrum.py` reports an effective tension, which is
the check that the spring is not adding one (WP6).

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
composes to the exchange; OpenMesh's `flip()` does the same. It is
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
that is more restrictive than the DTS convention `[3, 9]`; extending the
generator to `N = 3` (Loop's `β = 3/16` special case) and `N = 9, 10` is a
one-day follow-up once flips are running and the measured valence histogram
says whether it matters.

**Multi-extraordinary rows (WP1).** For a face of kind `MultiExtraordinary`,
`rebuild_patches` constructs, from the live adjacency:

1. the patch `P`: ordered union of the corners and their one-rings, width
   `K = N_0 + N_1 + N_2 - 6` in the generic case (fewer if rings overlap
   beyond the shared corners; the construction uses the actual sets);
2. `A_loc` (`(N_0 + N_1 + N_2) × K`): Loop's vertex rule with `β = 3/(8N)` on the
   three corners, the edge rule `(3/8, 3/8, 1/8, 1/8)` on every edge incident to
   a corner — every stencil entry is in `P` by §1.4(b);
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
   Loop limit mask for the actual valence, `(1/2, 1/(2N))`, as a sparse
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
| `edgeSpringEnabled` / `edgeSpringConstant` / `edgeSpringRestLength` | `false` / `kCurv` / `lFace` | equation (3) |
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

### WP1 — Multi-extraordinary patches by one local subdivision

`A_loc` generator, child selections, per-face composed rows with signature
cache, `PatchKind::MultiExtraordinary` in the CPU kernel loop and in
`DeviceMeshLayout::FacePatchDescriptor` (kind + width + row offset replace the
"valence from width" inference); `kMaxControlPoints` raised to 18 with the
regular path's buffers untouched.

> Gate: the refinement identity of §3.3 to `1e-12` relative on faces with two
> and three extraordinary corners at every valence combination in `[4, 8]`; a
> 6/6/6 face pushed through the local-subdivision path agrees with the direct
> kernel to quadrature accuracy; force–energy conjugacy (the existing test
> pattern) on a 5/5/7 and a 4/8/8 face; rigid-motion invariance; regular
> workload bit-identical.

### WP2 — Local patch evaluation

`FlipPatch`, `evaluate_patch()`, `ΔE` by equation (4).

> Gate: for random interior flips on the bowl grid and on the flat periodic
> fixture, `ΔE` from `evaluate_patch` equals `E_total(after) - E_total(before)`
> from two full `Compute_Energy_And_Force()` calls to `1e-9` relative, with
> `uSurf` and `uVol` nonzero so the global terms are exercised.

### WP3 — Metropolis sweep and Poisson schedule

The sweep of §3.5, the tagged counter-based RNG, `EdgeFlips.csv`, and the
log-determinant diagnostic of §1.4(c) as a test (four sparse solves and a 4×4
determinant on a lattice flip; the number goes into the docs).

> Gate: **two-state test** — freeze all edges but one on a small fixture, run
> the sweep at fixed coordinates for `10⁵` attempts, and check that the
> occupancy ratio of the two triangulations equals `exp(-ΔE/kT)` within
> statistical error, at two temperatures. Poisson gate: the attempt count over
> many sweeps has mean and variance `λ`. Reproducibility: same seed, same flips.

### WP4 — Fluid-mode dynamics

The three items of §3.7.

> Gate: with all flags off, `Run_dynamics_flat` on the spectrum workload is
> bit-identical (serial build) through the sparse-solver change; with
> `surfaceSolver = iterative` the trajectory agrees with dense to the solver
> tolerance; `M⁻ᵀ` versus `M⁻¹` on a mesh with mixed valences differs and the
> transposed one is the one the FDT test in `fluctuation_spectrum.md` accepts;
> the edge-spring force passes finite-difference conjugacy.

### WP5 — Wiring, output, GPU

Sweep call in `Run_dynamics_flat.cpp` after the end-of-step force evaluation;
version-based invalidation (§3.6); `face_<iteration>.csv`; checkpoint faces
block; parameters. Host side only if CUDA cannot be built here (it cannot —
see the GPU memory note); the device layout rebuild is exercised by the
`HostForceBackend` tests that already pin the layout against the CPU loop.

> Gate: a short flat run with flips on writes paired coordinate/face frames
> that reload into a consistent mesh; the layout rebuild fires exactly on
> steps with accepted flips.

### WP6 — Validation of fluidity and physics

On the periodic sheet and on an icosphere:

1. acceptance rate and valence histogram versus `ν` (DTS equilibrium meshes
   are roughly 60% valence 6, 20% each 5 and 7 — a sanity band, not a target);
2. **neighbour survival**: the fraction of initial edges still present after
   time `t`, whose decay time is the microscopic fluidity time scale and the
   thing `ν` should be calibrated by;
3. **in-plane MSD** of tagged vertices: saturating at the cage size without
   flips, linear in `t` with them — the classical signature of a fluid;
4. the fluctuation spectrum with flips on, through the resampled pipeline:
   the same `kc` from the `q⁻⁴` fit, and an effective tension that does not
   grow with `k_S`;
5. throughput per step against §6 risk 1, at depth scale 1.0 and 0.5.

> Gate: (3) and (4). Numbers recorded in a results document, as
> `irregular_patch_results.md` did for its plan.

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
| **A fluid mesh is mostly irregular, and irregular faces are expensive.** Today the row tables are exercised by a dozen faces on an icosphere. In a fluid steady state roughly half the vertices are not valence 6, so most faces have at least one extraordinary corner, each costing `3·D` samples per corner instead of 3 — at the current depths, `10–40×` the regular cost per face, so a global force evaluation on `data/example` goes from ~11 ms to of order 100 ms. | This is the real cost of fluidity and it is not specific to this design — any subdivision membrane that flips pays it. Levers, in order: `irregularPatchDepthScale` (the `1e-4` bending tail the depths were chosen for is far below the thermal noise a Brownian run lives in; WP6 measures what depth the spectrum actually needs); the GPU backend, which was built for exactly this kind of face-parallel load; a higher-order rule on the children so that a given accuracy needs fewer levels. WP6 item 5 gives the number before anyone commits to a workload. |
| **The measure question** (§1.4c). | Plain Metropolis at fixed `C` is the literature standard; the log-det diagnostic quantifies the discrepancy; the corrected acceptance is a documented option. Geometry sampled at fixed `T` is unaffected either way. |
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

- Gompper, G. & Kroll, D. M. *Triangulated-surface models of fluctuating
  membranes*, in *Statistical Mechanics of Membranes and Surfaces* (2004).
  https://www.researchgate.net/publication/285604301_Triangulated-surface_models_of_fluctuating_membranes
- Gompper, G. & Kroll, D. M. *Membranes with fluctuating topology: Monte Carlo
  simulations*, PRL 81, 2284 (1998). https://ui.adsabs.harvard.edu/abs/1998PhRvL..81.2284G/abstract
- Ramakrishnan, N., Sunil Kumar, P. B. & Ipsen, J. H. *Monte Carlo simulations
  of fluid vesicles with in-plane orientational ordering*, PRE 81, 041922
  (2010). https://arxiv.org/abs/1004.4509
- Ramakrishnan, N., Sunil Kumar, P. B. & Radhakrishnan, R. *Monte Carlo
  simulations of fluid vesicles*, J. Phys.: Condens. Matter 27, 273104 (2015).
  https://iopscience.iop.org/article/10.1088/0953-8984/27/27/273104
- Noguchi, H. & Gompper, G. *Fluid vesicles with viscous membranes in shear
  flow*, PRL 93, 258102 (2004). https://arxiv.org/abs/cond-mat/0404356
- Noguchi, H. & Gompper, G. *Dynamics of fluid vesicles in shear flow: Effect of
  membrane viscosity and thermal fluctuations*, PRE 72, 011901 (2005).
  https://journals.aps.org/pre/abstract/10.1103/PhysRevE.72.011901
- Sadeghi, M., Weikl, T. R. & Noé, F. *Particle-based membrane model for
  mesoscopic simulation of cellular dynamics*, J. Chem. Phys. 148, 044901
  (2018). https://arxiv.org/abs/1710.06907
- Siggel, M. et al. *TriMem: A parallelized hybrid Monte Carlo software for
  efficient simulations of lipid membranes*, J. Chem. Phys. 157, 174801 (2022).
  https://www.biorxiv.org/content/10.1101/2022.05.25.493239v1
- Pezeshkian, W. & Ipsen, J. H. *Mesoscale simulation of biomembranes with
  FreeDTS*, Nat. Commun. 15, 548 (2024). https://www.nature.com/articles/s41467-024-44819-w
- *OrganL: Dynamic triangulation of biomembranes using curved elements* (2024).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11213972/
- *PyMembrane: A flexible framework for efficient simulations of elastic and
  liquid membranes* (2023). https://arxiv.org/abs/2308.12754
- Stam, J. *Evaluation of Loop subdivision surfaces* (1998); the one-subdivision
  isolation argument that §1.4(b) applies per face.
