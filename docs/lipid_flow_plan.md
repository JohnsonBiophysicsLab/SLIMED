# Lipid flow: the stretching energy's relaxed area as a coordinate that moves

**Status:** implemented in the tree (2026-10-06), in two modes selected by
`lipidFlowMode` (section 7): `global`, one relaxed area for the whole sheet
with a reservoir rate law, and `local`, a relaxed area per face with
conserved fluxes, optional reservoir and optional noise. Seventeen tests in
`tests/test_lipid_flow.cpp`; the shipped periodic workload is byte-identical
with the feature off (both legacy workloads, 23 files). The local mode works
with edge flips and Monte Carlo vertex moves (section 5.1), and with
`lipidFlipResampling = patch` -- the default when flips are on -- keeps the
global mode's flip acceptance (section 6.5); the GPU backend,
the legacy `Periodic` boundary and the minimizer are still refused. Sections
1-4 and 6 are the plan and the prototype measurements that led here (the
prototype's patch, runs and `tension_propagation.png` are in
`~/slimed_runs/lipid_flow_prototype/`); section 6.5 adds what the in-tree
modes measured.
**Scope:** another way to give the membrane in-plane fluidity than moving
vertices and flipping edges: every face keeps its own relaxed area `S0_f`,
the stretching energy is local, and `S0_f` moves between faces down the
tension gradient and, optionally, in and out of a reservoir. The mesh nodes
stay where they are. Both are available together: the field handles *area*
(density, tension), the fluid-mode moves handle *shape* and connectivity.

## 1. The idea

With the global constraint the stretching energy is `(K/2A0)(A - A0)^2`: one
tension for the whole sheet, equalized instantly. That is the limit of
infinitely fast lateral lipid flow. Giving the flow a rate means tracking
where the lipids are. The natural variable on this mesh is the relaxed area
of each face, the area its lipids would cover unstretched:

```text
    e_f = (K / 2 c_f) (S_f - S0_f)^2 ,     sigma_f = de_f/dS_f = K (S_f - S0_f) / c_f ,
```

with `S_f` the limit-surface area of face `f` (already computed every
step), `S0_f` its relaxed area, and `c_f` a fixed per-face normalisation
(section 3.1). The force on the control points is then the face's own
tension times the gradient of the face's own area, which is the kernel's
existing area term with a different scalar (section 5). `S0_f` is a
coordinate of the model and moves by three mechanisms, each the gradient
flow of the same energy:

- **lateral flow, conserved.** Across every edge shared by faces `f` and `g`
  a flux `M w_fg (sigma_f - sigma_g)` into `f`, with `w_fg` the edge length
  over the distance between the two face centres (the two-point flux of a
  finite-volume scheme) and `M = 1/b` the mobility of lipids sliding
  against whatever the mesh frame stands for. To linear order this is
  `d sigma/dt = D_sigma Laplacian(sigma)`, `D_sigma = K M`: tension
  diffuses (Shi et al. 2018, their eqs. S10-S11).
- **a reservoir**, on the rim or everywhere: `dS0_f/dt = (k c_f / K)(sigma_f -
  sigma_res)`, which relaxes the local tension to `sigma_res` at rate `k`.
- **noise**, so that the relaxed areas sample their Boltzmann distribution:
  `sqrt(2 kT M w dt)` on each flux, `sqrt(2 kT k c_f dt / K)` on each
  reservoir exchange.

The limits are models the tree already has. `M = 0` is purely local
elasticity. `M -> infinity` with no reservoir makes the tension uniform and
the energy exactly `(K/2A0)(A - A0)^2`. A reservoir everywhere at `sigma_res`
is the constant-tension Helfrich term, in a fixed frame and without a box.

This is an Eulerian description of in-plane flow: lipids move through a
mesh whose nodes do not. The existing fluid mode (in-plane Brownian vertices,
tether, flips) is the Lagrangian one. They are not rivals; the Eulerian
field handles *area* (density, tension), the Lagrangian moves handle
*shape* and connectivity.

## 2. What the literature says

A survey of the tension-propagation, density-field and reservoir
literature was made for this plan (2026-10-05; extracted sources under the
session scratchpad, not kept). What matters here:

- **The continuum limit is established.** With friction against an
  immobile frame the equations above reduce to tension diffusion,
  `d sigma/dt = D_sigma Laplacian(sigma)`, `D_sigma = E k/eta` with `k/eta`
  the Darcy mobility of the obstacle array (Shi, Graber, Baumgart, Stone &
  Cohen, *Cell* **175**, 1769, 2018; Cohen & Shi, *BioEssays* **42**, 1900142,
  2020). The per-face energy is the compression energy of Sahu, Sauer &
  Mandadapu (*PRE* **96**, 042409, 2017, eqs. 153-154) with `K = 2 k_c`.
- **The rate is a friction, and a free bilayer has almost none.** In water
  the mean-density mode of a bilayer relaxes at `k q/(eta_2 q + 2 eta)`,
  about 10-100 ns at every sub-micron wavelength (Bitbol & Fournier 2011,
  eq. 26; Bingham, Smye & Olmsted 2015, eq. 15; Rahimi & Arroyo 2012 and
  Torres-Sánchez, Millán & Arroyo, *JFM* **872**, 218, 2019: "t4 ~ 10 ns").
  So for a free bilayer at microsecond steps the existing global constraint
  *is* the right model. Finite-rate propagation is physical where there is
  an external friction: a cortex with immobile obstacles (HeLa: `D_sigma`
  0.024 nm^2/us with the cell's effective modulus, 100-600 nm^2/us with the
  bare one), a substrate (about 50 nm^2/us), a skeleton (red cell, about
  1700 nm^2/us).
- **Nobody stores and transports a per-element reference area.** The
  closest work is the Arroyo group's Loop-subdivision finite elements with
  nodal monolayer densities advected by an ALE continuity equation
  (Torres-Sánchez et al. 2019, implicit with SUPG, mass conserved to 1e-2%);
  Ma & Klug (*J. Comput. Phys.* **227**, 5816, 2008) relax a reference
  configuration as a numerical dashpot; red-cell models (Fedosov, Peng,
  Fai) carry per-triangle area terms with *fixed* rest areas; Kumar,
  Gompper & Lipowsky (*PRL* **86**, 3911, 2001) keep a conserved composition
  on triangles across edge flips, the precedent for a face variable on a
  flipping mesh; Mem3DG transports a protein density on vertices with a
  cotangent Laplacian, non-conserved by default.
- **Reservoir rate laws in cells are slow.** Measured fractional area
  supply rates are 1e-3 to 1e-2 per second per pN/nm (Wolfe et al. 1985;
  Raucher & Sheetz 1999; Sens & Turner 2006; Sinha et al. 2011): seconds,
  frozen on a SLIMED run. A reservoir term in this code is therefore a
  tension-stat or a boundary condition; its physics is in the set-point
  (0.003-0.3 pN/nm in cells) and its capacity, not its rate.
- **Numerics.** Two-point fluxes are conservative, symmetric and monotone
  when the face centres are circumcentres and the mesh is locally Delaunay
  (Droniou 2014; Lenz, Nemadjieu & Rumpf 2011); with a constant mobility,
  one Gaussian per edge satisfies discrete fluctuation-dissipation with no
  drift term (Donev et al. 2010, eqs. 30-32; Delong et al. 2013, eq. 3);
  explicit Euler inflates a mode's variance by `1/(1 - lambda dt/2)`
  (Donev eq. 44), the implicit midpoint rule gets the statics exactly for
  any step (Delong eqs. 22-23).
- **Independent derivation of section 3.1.** The survey derived, without
  seeing the prototype, that integrating out a thermalised conserved `S0`
  field returns exactly the global Boltzmann weight when the energy is
  quadratic in `S0` with a fixed normalisation, and leaves an entropic
  tension offset of about `+kT/a_f` (+0.39 pN/nm at a 5 nm edge) with the
  `S0`-normalised form. Both were then measured (section 6).
- **A conserved face field under edge flips** (second survey, 2026-10-06,
  for section 5.1). No source derives a special acceptance rule; the codes
  split into two camps. *Dynamics* codes pool the two faces' quantity and
  redistribute it to the two new triangles -- Batty, Uribe, Audoly &
  Grinspun (*ACM TOG* **31**, 113, 2012) give each new triangle the same
  thickness, i.e. split by new area; the ALE remapping literature
  (Margolin & Shashkov 2003; Kenamond, Kuzmin & Shashkov 2020) intersects
  old and new cells. Those maps conserve exactly but are averaging
  contractions: not invertible, so a flip and its reverse do not compose to
  the identity, and in a Metropolis chain they would equalise the strain of
  every flipped pair at the flip rate, independent of the energy. *Sampling*
  codes carry the face's value with the face record and let the flip's
  energy difference judge the new areas: Kumar, Gompper & Lipowsky 2001
  (composition on triangles, `A_I` changes, `s_I` stays), the
  Ising-on-triangles dynamical-triangulation literature (David 1993; Janke &
  Weigel 2006, spins on the dual `phi^3` graph, which a flip only rewires),
  and OrganL (Allolio, Fábián & Dostalík 2024: "as the number of faces is
  conserved, face properties can remain associated with a fixed face;
  currently no redistribution occurs"). Baillie & Johnston (1992) are the
  precedent for the one case where arithmetic is unavoidable -- a move that
  creates a cell -- and they draw the new value at random so the reverse
  move has a matching proposal density. Mem3DG (Zhu, Lee & Rangamani 2022)
  calls conservative interpolation across mesh mutations an open problem.
  Involutive MCMC (Neklyudov et al. 2020) is the general statement: a
  deterministic proposal is valid with the plain Metropolis ratio when it
  is its own inverse with unit Jacobian, which a permutation of face values
  is and a pool-and-split is not. Lenz, Nemadjieu & Rumpf (2011) is the
  finite-volume side: cell mass is the state, area change only changes the
  density, mass moves only by edge fluxes, and the vertices are material
  markers -- which is what a vertex move must leave alone.

## 3. Decisions the prototype settled

### 3.1 Normalise the stiffness by a fixed area, not by `S0_f`

Two forms were tried. `(K/2 S0_f)(S_f - S0_f)^2` has the modulus follow the
lipids present, which reads as more physical, but its tension goes as
`1/S0_f`, so thermal fluctuations of `S0_f` shift the *mean* tension by a
mesh-dependent amount: `+kT/a_f` with a conserved total, `-kT/(2 a_f)`
against a reservoir, 0.39 and 0.19 pN/nm at `lFace = 5`. Measured on a
frozen flat sheet: the face-mean tension shifted by exactly `kT/a_f` within
5% (`LipidFlowTest.ThermalRelaxedAreasLeaveTheMeanTensionAloneOnlyUnderTheLinearForm`),
and in a thermal run the face-mean tension read 1.96-2.02 pN/nm against 1.68
for the global constraint on the same sheet. With `(K/2 c_f)(S_f - S0_f)^2`
and `c_f = S0_f(0)` the tension is linear in `S0_f`, its conjugate force is
the tension itself, the fluctuations leave the mean exactly alone, and the
equilibrium is the global constraint's to the digit (section 6.2). The
price is that a face which has gained 35% of relaxed area, as the faces on a
wrapped bead do, is 35% stiffer per unit excess area than the lipids present
would make it; since the constraint only ever holds `S_f` near `S0_f`, that
is a second-order effect. `lipidFlowLinearStretching = true` is the default.

### 3.2 Slow lipids need the noise

Without noise the statics depend on the rate (section 6.2): frozen or slow
lipids over-stiffen the sheet by nearly the factor two a mean-field
estimate gives for a local quartic term, fast flow recovers the global
result to 1-2%, and in between the system is a two-temperature one. With
the noise on, statics are independent of `D_sigma`, as detailed balance
requires. So: noise on for equilibrium work at any rate; noise off is
defensible only for driven problems at a rate fast against the shape
modes, where it saves nothing worth having but the cleaner logs.

### 3.3 The rate need not be a free parameter

If lipids feel the same friction per area as the vertices do, `b = kT/(D a_v)`
and `D_sigma = K a_v D / kT`: 1298 nm^2/us at the defaults (`K = 250`,
`a_v = 21.65 nm^2`, `D = 1`, `kT = 4.17`). That is `lipidFlowDiffusivity = -1`.
Under it the tension mode is faster than the bending mode at every
wavelength on the mesh, by `K/(kappa q^2) >= 4`, which is the right ordering
for a membrane. It is also squarely in the literature's friction-dominated
range (section 2), and it costs three explicit substeps at `dt = 0.001`.

## 4. Ways to implement, and which to take

| way | what moves | what it gives | cost |
| --- | --- | --- | --- |
| **A. one relaxed area for the sheet, with a rate law** `dA0/dt = k A0 (sigma - sigma_res)/K` | nothing local | a constant-tension reservoir for any geometry, the clamped benchmark sheet included; the `D -> infinity` limit of C | **`lipidFlowMode = global`**; the cost of the global constraint; works with flips and vertex moves |
| **B. per-face `S0_f`, local tension, local rate law only** | each face toward `sigma_res` | tension buffering everywhere (a caveolae-like reservoir); no propagation | C with `lipidFlowDiffusivity = 0` and a bulk reservoir |
| **C. per-face `S0_f` with conserved edge fluxes, optional reservoir, optional noise** | lipids between faces | finite-rate tension propagation, exact conservation, the global constraint as a limit | **`lipidFlowMode = local`**; 6-10% of run time |
| D. nodal density with a cotangent or Loop-basis Laplacian | lipids between vertices | second-order accuracy; couples to in-plane node motion | a mass-matrix solve per step; per-quadrature-point tension in the kernel |
| E. two leaflet densities with intermonolayer friction | lipids per leaflet | area-difference dynamics, `kappa~ = kappa + 2 k e^2`, the slip mode (Seifert & Langer 1993; Bitbol & Fournier 2013) | a second field coupled to curvature; the real bilayer physics at these scales |
| F. Lagrangian: `S0_f` attached to faces, nodes moved in-plane by the local tension gradient | the mesh | compressible-fluid kinematics without flips | needs a conservative rule for `S0` under flips |

**Both A and C are in the tree.** A is the `global` mode: it assumes the
relaxed areas equalise the tension -- infinitely fast flow, or any finite
rate with the noise on, which section 6.2 showed gives the same statics --
and keeps only the one coordinate `A0` that the equalised state has left.
It is the cheaper and more approximate dynamics. C is the `local` mode, the
finite-rate one. D is the upgrade if the two-point flux proves too crude on
disordered fluid-mode meshes; E is the physics extension, not a replacement.

## 5. Where it is in the code

| piece | where | change |
| --- | --- | --- |
| mode | `Param::lipidFlowMode` (`off`, `global`, `local`), `Mesh::lipid_flow_global()/local()/on()` | validated in `DynamicMesh::finish_dynamic_setup()`; any other string throws |
| global mode | `Mesh::lipidAreaNormalization`, `DynamicModel::advance_global_lipid_flow()` | the stretching energy is `(K / 2 A0_built)(A - A0)^2`, its modulus normalised by the relaxed area *as built* so that `A0` can move without changing the stiffness; `dA0/dt = (k A0_built / K)(sigma - sigma_res)`, with `sqrt(2 kT k A0_built dt / K)` noise when thermal; without a reservoir it is the global constraint to the byte |
| local state | `Mesh::faceArea0`, `faceArea0Initial`, `faceTension`, `lipidLinks`, `lipidRimFace` | per-face arrays; one link per physical edge, seam included |
| force | `Mesh::compute_face_energies_and_forces()` | the kernel's `areaScale = uSurfPerArea (area - area0)` is made `sigma_f` by setting the three fields per face; `face.energy.energyArea` set per face; step 5's global energy skipped in the local mode, normalised in the global one |
| graph | `Mesh::build_lipid_links()` (`src/dynamics/Lipid_flow.cpp`) | from the edge table and the typed images: a copy face stands for the carrying face with the same three source corners; a seam edge, met from both sides, is deduplicated by its source-vertex pair; rebuilt lazily when `topologyVersion` has moved |
| step | `DynamicModel::advance_lipid_flow()` | explicit flux update with automatic substeps, reservoir, conserved noise through the counter-based generator; after the vertex step and the box, before the vertex copy |
| trials | `Mesh::evaluate_edge_flip()`, `Mesh::evaluate_vertex_move()`, `FaceSubsetEnergy::stretching` | the Metropolis difference uses the per-face stretching energies in the local mode and the normalised global term in the global mode; pinned against the full energy by two tests |
| flips | `Mesh::lipid_carry_after_flip()`, called by the sweep and inside the trial | section 5.1 |
| logs | `<prefix>lipid.csv`, `lipidface<prefix>.csv` | totals, face-mean and pooled tension, min/max/std, rim tension, reservoir exchange; per-face `(S0_f, S_f)` with `localOutput`; the global mode writes its one tension in every tension column |
| parameters | `lipidFlowMode`, `lipidFlowDiffusivity`, `lipidFlowSubsteps`, `lipidFlowThermal`, `lipidFlowLinearStretching`, `lipidReservoirRegion`, `lipidReservoirTension`, `lipidReservoirRate` | parsed, printed, pinned by a test; documented in `input.params` |
| examples | `data/example/lipid_flow_bead_rim_reservoir.params`, `data/example/lipid_flow_global_reservoir.params`, `data/example/lipid_flow_global_flat.params` | the bead of section 6.3 under each mode; the flat sheet with no bead, A0 relaxing to the reservoir tension |

### 5.1 Edge flips with a per-face field

The prototype pooled `S0` and `c` of the two faces and split them by the
new areas, the Batty et al. rule. The literature of section 2 says that is
the wrong rule for a Metropolis move, and the reason is simple: the split
is a projection, so a flip followed by its reverse does not return the
field, there is no reverse proposal for the Hastings ratio to use, and the
chain equalises the strain of every flipped pair at the flip rate whatever
the energy says. What the tree does instead:

- **Carry.** Each old triangle's `(S0_f, c_f)` is re-attached whole to one
  of the two new triangles. No arithmetic: the sum is conserved bitwise, and
  the proposal is a permutation, so the flip is accepted with the plain
  `min(1, e^{-dE/kT})`, `dE` including `(K/2c)[(S_f' - S0_f)^2 - (S_f - S0_f)^2]`
  for both faces -- a flip that changes the two areas is charged its
  stretching cost, which is the Kumar-Gompper-Lipowsky picture. The lipids
  are *declared* to be under the new triangle; the tension imbalance that
  leaves is what the flux step then relaxes. That is the standard compromise
  (OrganL makes it; Kumar et al. implicitly), and it means flips add a
  bookkeeping contribution to the lipid diffusivity at the flip rate.
- **An involutive pairing.** Which old triangle goes to which new one has
  to read the same from both sides of the flip, or the carry is not
  reversible either. The rule is "the larger old triangle takes the larger
  new one", with flat triangle areas computed in a canonical corner order
  so both sides see the same bits; on a tie it defers to `flip_edge()`'s
  own rule. Pairing by areas keeps the lipid mass the flip relabels as
  small as a pairing can (it is the max-overlap matching of the two old
  and two new triangles in the quad's plane). Pairings that minimise the
  *energy* of the proposal are not self-inverse and are not used.
- **`flip_edge()` itself had to change.** Its labelling -- face index 0
  takes the new triangle on the side of its own first endpoint -- composes
  to an exchange: two flips left the two face indices holding each other's
  triangle, which `docs/edge_flip_plan.md` section 3.2 had recorded as
  intrinsic. It is not: pairing the record with the smaller old corner
  index to the new triangle with the smaller corner index is an involution,
  because corners and endpoints trade roles under the flip. Before the
  change, every *rejected* trial silently permuted everything keyed by face
  index without bumping `topologyVersion`: the field, the flux graph built
  on the last accepted flip, and the per-face output frames against a
  connectivity frame that no longer matched. `EdgeFlipTest.FlippingTwiceRestoresTheMeshIncludingTheFaceLabels`
  now demands every face's corner list back (up to a rotation of where the
  list starts), and the trial restores the field by carrying it back
  through the reverse flip rather than by saving it. The legacy workloads
  do not flip and stayed byte-identical.
- **Vertex moves** leave the field alone: faces persist, `S_f` changes, the
  trial charges the stretching difference. Vertices are the material
  markers of the finite-volume picture; all relative lipid motion goes
  through the flux step. The field is never remapped through the swept
  areas of a move.
- **Resampling (`lipidFlipResampling`; `auto`, the default, is `patch`
  whenever flips are on under the linear form).** The carry rule alone
  charges a flip's area change at the local modulus and makes flips rare
  (section 6.5). The exact way to get the fluidity back is a compound
  move: flip, carry, then redraw the relaxed areas of a set `P` -- the two
  faces, or the whole flip patch -- from their conditional Boltzmann
  distribution at the new geometry with `sum_P S0` fixed. Under the linear
  form the set's stretching energy is quadratic with Hessian `diag(K/c_f)`,
  independent of the geometry, so the Gaussian normalisations cancel in
  the Metropolis-Hastings ratio and the acceptance is `exp(-dE/kT)` with
  the set's stretching energy *minimised over its split*,
  `(K / 2 sum c)(sum S - sum S0)^2`, before and after
  (`Mesh::lipid_stretching_energy_with_pooled()`), the rest of the patch at
  the relaxed areas it carries. The redraw
  (`Mesh::lipid_resample_after_flip()`) is the equal-tension split plus
  conserved noise, `S0_f = S_f - sigma* c_f/K + eta_f - (c_f / sum c) sum eta`,
  `eta_f ~ N(0, kT c_f / K)`, keyed on the attempt so it is reproducible;
  with the noise off it is the split alone. This is a partial Gibbs step on
  `(connectivity, S0_P)`: exact for sampling, and for the dynamics an
  instant redistribution of lipids over `P` at every accepted flip -- which
  the flux step at `D = 1298` would do in some 20 (`pair`) to 200 (`patch`)
  steps anyway, but at a cortex-like `D = 50` is a shortcut. Tests:
  `FlipResamplingChargesThePooledEnergyAndLeavesTheSetAtOneTension` (both
  sets; the trial's charge equals the full energy's change given the
  redraw, the set ends at one tension, the total is kept),
  `ThermalFlipResamplingKeepsTheTotalAndSpreadsTheTension`.

Tests: `AnEdgeFlipConservesTheRelaxedAreaAndItsTrialMatchesTheFullEnergy`
(the multiset of relaxed areas is unchanged bitwise, the trial's `dE`
matches a full re-evaluation to 1e-6), `AnEdgeFlipCarriesTheFieldAsAnInvolution`
(every flippable edge of a jittered sheet, flipped and flipped back through
the sweep's code path, puts every value back on its own triangle bitwise;
some flips swap the records and some do not, so the pairing is exercised),
`AVertexMoveTrialMatchesTheFullEnergyUnderLocalFlow`.

### 5.2 Still missing

1. **GPU.** A per-face `areaScale` array in `ForceKernelArgs` instead of one
   scalar. Refused for now (`forceBackend = cpu`), both modes.
2. **Legacy `Periodic` mode.** Its ghost faces are not typed, so the seam
   links cannot be built from the images; the local mode requires `Mixed`
   for periodic sheets. The global mode does not care.
3. **The minimizer.** At a minimum the tension is uniform, so it needs no
   field: the global constraint, plus a constant-tension term `sigma_res A`
   when a reservoir is wanted (the global mode's end state).
4. **Circumcentres and a Delaunay check** for the transmissibilities on
   disordered meshes (Droniou 2014); the links use centroids, which are
   exact on the lattice.
5. **Stiff rates.** Crank-Nicolson on the constant flux matrix for
   `D_sigma dt / a^2` beyond the explicit window, or simply the global mode,
   which is that limit.
6. **Restart files**, the fluid-run analysis (`fluidity.py` reads areas),
   and per-area rather than per-vertex bead adhesion, which becomes natural
   once faces carry areas.
7. **The flip-induced diffusivity** of section 5.1: tag faces with a
   passive scalar carried like `S0` and measure its spread against the flip
   rate at `D = 0`. Not measured yet. With resampling on, the redraw adds
   its own: one `kT c_f / K` kick per face of the set per accepted flip.

## 6. What the prototype measured

All runs: `~/slimed_runs/lipid_flow_prototype/`, serial Release binary of
the patched tree, `K = 250 pN/nm`, pinned vertices (`inPlaneDynamicsEnabled =
false`), `kc = 83.4`, `lFace = 5`, `D = 1`, `kT = 4.17`. Analysis scripts in
its `analysis/`.

### 6.1 Unit tests (`tests/test_lipid_flow.cpp`; the prototype's ten, plus the seven of sections 5 and 5.1)

- the force on the control points is minus the gradient of the local
  energy, both forms, three probe vertices, all axes, 1e-5 relative;
- uniform relaxed areas on a uniformly stretched sheet reproduce the global
  constraint's energy and per-vertex forces to 1e-8;
- the flux graph has one link per physical edge: every face of a sheet
  that wraps both ways has exactly three, through the seam, and a clamped
  sheet has a rim;
- the lateral flow conserves the total to 1e-11, never raises the energy,
  and ends at a uniform exchange tension where the energy equals the
  global constraint's to 1e-8, both forms;
- a tension wave on a 100 nm periodic sheet decays at `D_sigma q^2` within 3%;
- with the noise on, `var(S0_f)/c_f = (kT/K)(1 - 1/N)` within 5%;
- thermal fluctuations leave the face-mean tension alone to 1e-9 under the
  linear form and shift it by `kT/a_f` (within 5%) under the other;
- a bulk reservoir brings every face to its tension exactly;
- a non-CPU backend, legacy `Periodic` under the local mode, and an unknown
  mode string are refused;
- the global mode without a reservoir reproduces the global constraint's
  energy and forces exactly, and with a reservoir brings the sheet to its
  tension to 1e-9.

### 6.2 Statics on a fixed periodic sheet (100 nm, `Mixed`, 200 000 steps at `dt = 0.002`)

Means over 100-400 us; the spectrum from `fluctuation_report.py --burn 500`.

| run | face-mean tension | spectrum tension | `<E_curv>` | rms height (nm) |
| --- | ---: | ---: | ---: | ---: |
| global constraint | 1.68 | 1.42 | 490 | 0.775 |
| local, lipids frozen (`D = 0`) | 1.47 (pooled) | 2.37 | 469 | 0.695 |
| local, `D = 2`, no noise | 1.47 (pooled) | 2.32 | 469 | 0.698 |
| local, `D = 1298`, no noise | 1.65 (pooled) | 1.47 | 485 | 0.766 |
| local, `D = 2`, noise, per-lipid form | 1.96 | 1.80 | 481 | 0.737 |
| local, `D = 1298`, noise, per-lipid form | 2.02 | 1.78 | 480 | 0.739 |
| **local, `D = 2`, noise, linear form** | **1.66** | **1.50** | **489** | **0.766** |
| **local, `D = 1298`, noise, linear form** | **1.67** | **1.44** | **489** | **0.773** |
| bulk reservoir 0.5, noise, linear | 0.50 | 0.25 | 528 | 1.01 |
| bulk reservoir 0.0, noise, linear | 0.00 | (−0.38, slope −4.8) | 557 | 1.37 |

Frozen and slow lipids read a spectrum tension 1.65 times the global one,
the Hartree factor of a local quartic term (`sigma_eff = 2 K <alpha>` against
`K <alpha>`; with the self-consistent excess area the prediction is 1.75).
The per-lipid thermal form shifts the face-mean tension by `+kT/a_f` and
the spectrum by about +0.37, both as predicted. The linear thermal form
reproduces the global constraint to within noise at either rate. A reservoir
at zero *internal* tension in a fixed frame leaves the frame tension
negative (`-N kT/(2 A_p)`, the Schmid residual); the lowest mode is four
times too large and the sheet is near buckling, as Shiba, Noguchi &
Fournier (2016) describe for the `(A_p, sigma)` ensemble. A reservoir in a
fixed frame needs `sigma_res` above about 0.1 pN/nm on this tile.

### 6.3 A bead on the clamped benchmark sheet

The geometry of `docs/bead_push_benchmark.md`, pinned mesh, no adhesion: a
10 nm bead pushed from 5 nm below the sheet at 0.5 nm/us for 70 us (35 nm),
then held 20 us. The clamped sheet has no reservoir, which is why the
benchmark had to soften `K` to 5; the elastic box of `docs/elastic_box.md`
cannot act here because nothing wraps. The lipid flow can: the rim faces
exchange with a reservoir at 0.5 pN/nm (`k = 200 /us`).

| run | tension at the bead / at the rim, 60-70 us | force, 60-70 us | force, held 80-90 us |
| --- | ---: | ---: | ---: |
| global constraint, `K = 250` | 17-19 (uniform) | −772 pN | −779 pN |
| global, `K` softened to 5 | 0.6 (uniform) | −90 | −89 |
| rim reservoir, `D = 10 000` | 0.9 / 0.56 | −98 | −79 |
| rim reservoir, `D = 1298` | 2.1 / 0.62 | −143 | −72 |
| rim reservoir, `D = 100` | 7.6 / 1.1 | −383 | −87 (after 15 us of relaxation) |
| bulk reservoir, `D = 1298` | 1.1 / 0.53 | −98 | −77 |
| rim reservoir, `D = 1298`, noise, linear form | 0.8 / 0.5 | −139 | −73 |

The reservoir supplied about 1050 nm^2, the tent's excess area. During the
push the tension at the bead rises above the rim's by an amount that scales
as `1/D_sigma`, roughly `(Q K / 2 pi D_sigma) ln(R_rim / r_bead)` with `Q`
the area demand (about 25 nm^2/us here): 1.3 pN/nm predicted against 1.5
measured at `D = 1298`. The force on the bead follows the tension at the
bead, so a pushed bead's force depends on how fast it is pushed relative to
the lipid flow, which is the phenomenon Shi et al. measure in cells. In the
hold the tension relaxes to the rim's in about `L^2/D_sigma`. See
`tension_propagation.png`.

### 6.4 Cost

On the clamped sheet (525 vertices, 792 faces), 6000 steps, serial Release:
12.65 s with the global constraint, 13.42 s with the flow at `D = 1298`
(three substeps), 13.97 s at `D = 10 000` (twenty substeps). The flux update
is a sweep over 1148 edges; the force evaluation dominates.

### 6.5 The in-tree modes on a flipping fluid sheet (`~/slimed_runs/lipid_flow_modes/`)

The fluid sheet of `docs/edge_flip_plan.md` (100 nm `Mixed` periodic, 504
carrying faces, in-plane Brownian vertices, flat tether 0.7-1.75, flips at
0.5 attempts per sweep, valences 4-8, `K = 250`), 30 000 steps at
`dt = 0.001`, OpenMP Release build. `fluid_global_ref` is the global
constraint; the two local runs use `D = -1` (1298), noise on, linear form.

| run | flip rule | admissible attempts | accepted | median trial `dE` (pN nm) | `<E_curv>` | `<E_area>` | mean tension | `sum S0` first / last |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| global constraint | -- | 6097 | 1777 (29.1%) | 32 | 536 | 26 | 1.54 | -- |
| local, pool-and-split (prototype rule) | not reversible | 7398 | 101 (1.4%) | 113 | 530 | 1357 | 1.57 | 5455.96004384 / 5455.96004384 |
| **local, carry (section 5.1), `lipidFlipResampling = none`** | involution | 7901 | 57 (0.7%) | 145 | 528 | 1349 | 1.59 | 5455.96004384 / 5455.96004384 |
| local, carry + `pair` resampling | partial Gibbs, 2 faces | 7971 | 93 (1.2%) | 116 | 529 | 1349 | 1.60 | 5455.96004384 / 5455.96004384 |
| **local, carry + `patch` resampling** | partial Gibbs, 18 faces | 6279 | **1868 (29.7%)** | 30 | 535 | 1352 | 1.52 | 5455.96004384 / 5455.96004384 |

The shape statistics and the mean tension agree across all five. `E_area`
of the local runs is the field's own equipartition, 504 faces at `kT/2`
each, not a stretched sheet (per-face tension noise 10.7 pN/nm, section 8).
The conserved total is exact to the printed digit over 30 000 steps in
every run, the resampled ones included; the carry run also exercised
`flip_edge()`'s new labelling through 7900 rejected trials with no drift in
the field.

**`patch` resampling gives the local mode the global mode's fluidity.**
Acceptance 29.7% against 29.1%, the same median trial `dE` (30 against 32),
the same accepted-flip `dE` (0.10 against 0.11), the same bending energy
(535 against 536). That is what the arithmetic says it should do -- with
the eighteen faces pooled the flip is priced at `K / (18 c)` = 1.3 pN/nm
per nm^2 against the pair's 11.6, and most of a flip's area change cancels
inside its own patch -- but it is also a check that the compound move is
what section 5.1 claims: if its acceptance were being bought with a wrong
ratio, the shape statistics would have moved. `pair` resampling lands where
the prototype's pool-and-split did (1.2% against 1.4%), as it must: the two
share the same acceptance formula, and the prototype's rule was the `pair`
redraw with the noise left out and no reverse move.

**Fluidity is the price of a per-face field, whichever flip rule.** A flip
changes the two faces' areas, and the local energy charges that change at
the pair's modulus, `K / (c_f + c_g)` = 11.6 pN/nm per nm^2, against the
global constraint's `K / A0` = 0.046: a 1.5 nm^2 change costs about 3 kT
plus its cross-term with the pair's thermal strain (about 7 pN/nm), where
the global constraint charges `sigma dA`, about 2 pN nm. That is why the
median trial `dE` is 145 against 32 and the acceptance 40 times lower. The
pool-and-split rule halved the penalty (it judged the flip at the pair's
equalised tension) and was still 20 times below the global mode; it was not
a fix, and it was not reversible. The penalty is physical *for the state
space chosen* -- the lipids under a face are declared to belong to it -- but
a flip is mesh bookkeeping, not lipid physics, so the fluidity it removes
is an artefact of pinning lipid patches to triangles.

The fix that keeps detailed balance exact is the resampling of section
5.1, and the table shows it working. Its cost is kinetic, not statistical:
every accepted flip redistributes lipids instantly across its patch --
which at `D = 1298` the flux step does in about 200 steps anyway, but at
the cortex-like `D = 50` is a shortcut of some 5000 steps, applied at about
0.06 accepted flips per step on this sheet, i.e. roughly every face is
redrawn once per 300 steps. So `patch` is the default whenever flips are on
under the linear form (`lipidFlipResampling = auto`): at the natural `D` it
changes nothing measurable and keeps the mesh fluid. Set `none` explicitly
when the tension gradient itself is the object at a slow `D`, accepting the
low flip rate or turning flips off. The global mode remains the cheapest
answer when no gradient is wanted.

### 6.6 The clamped bead with the global mode as a reservoir

`bead_global_reservoir`: the bead of section 6.3 (`K = 250`, pinned mesh,
no adhesion, 35 nm push over 70 us, 20 us hold) with `lipidFlowMode =
global`, bulk reservoir at 0.5 pN/nm, `k = 20 /us`. The tension sits at
0.52 +- 0.11 pN/nm over the second half (the reservoir at 0.5, noise on);
the sheet took 1092 nm^2 from the reservoir, the tent's excess area, as the
rim reservoir of section 6.3 supplied about 1050; the held force is -104 to
-107 pN against -779 for the global constraint at the same `K` and -72 to
-87 for the local rim reservoir. The global mode is the right tool for this
benchmark: it answers "what does a bead cost when the membrane can draw
area at tension `sigma_res`", costs nothing over the global constraint, and
only loses the propagation delay of section 6.3, which is a `1/D_sigma`
effect the local mode exists to show.

## 7. Parameters

```text
lipidFlowMode = off              # off | global | local
lipidFlowDiffusivity = -1        # local: D_sigma, nm^2/us; -1 = K a_v D / kT (1298 at the defaults); 0 = frozen lipids
lipidFlowSubsteps = 0            # local: 0 = automatic, keeping the fastest mode's rate x dt at or under 0.5
lipidFlowThermal = true          # both: noise on the relaxed area(s) at kT (section 3.2)
lipidFlowLinearStretching = true # local: stiffness on the face's initial relaxed area (section 3.1)
lipidReservoirRegion = none      # none | rim | bulk (the global mode treats rim and bulk alike)
lipidReservoirTension = 0.5      # sigma_res, pN/nm
lipidReservoirRate = 20          # k, 1/us; a tension-stat rate, not a biological one (section 2)
lipidFlipResampling = auto       # local + flips: auto | none | pair | patch -- redraw the relaxed areas of the
                                 # flipped pair / the flip patch at each accepted flip (sections 5.1, 6.5); needs
                                 # the linear form. auto = patch when flips are on under the linear form, else none
```

Which mode: `global` when the question is shape, force or a reservoir and
the lipids can be assumed equilibrated across the sheet within a step --
the clamped bead benchmark, a tension-controlled sheet, anything with
flips where fluidity matters (section 6.5). `local` when the tension
*gradient* is the object: a pushed or pulled membrane against a friction,
tension propagation times, a rim reservoir feeding a distant bead.

Ranges: `D_sigma` 100-600 nm^2/us for a cortex-bound membrane with the bare
modulus, about 50 for a supported bilayer, 1298 for the vertices' own
friction; the explicit update is stable below `a^2/(12 D_sigma dt)`, which
the automatic substeps handle up to a few thousand; above that use the
global mode, which is the limit. With a reservoir in a fixed frame keep
`sigma_res >= N kT/(2 A_p)`, about 0.1 pN/nm on the 100 nm tile.

## 8. Limits and open questions

- A single relaxed area is the symmetric bilayer mode. Intermonolayer
  friction, area-difference elasticity and the curvature coupling that
  stiffens `kappa` to `kappa + 2 k e^2` need two leaflet densities (way E);
  at these scales that is the mode that actually relaxes slowly in a free
  bilayer.
- Two-point fluxes with centroids are exact on the lattice and first-order
  on a disordered mesh; a fluid-mode mesh wants circumcentres and a
  Delaunay check.
- The per-face thermal tension noise is `sqrt(K kT / a_f)`, 9.8 pN/nm per
  face: instantaneous tension maps are noise, means converge as
  `1/sqrt(faces x frames)`, and any threshold or exponential reservoir law
  must act on averages.
- The pinned (Monge-like) mesh cannot resolve a wrapped bead; the flow
  removes the area artefacts of that mesh, not its geometric ones.
- With flips and `lipidFlipResampling = none`, the local mode's flip
  acceptance is set by the local modulus, not by the geometry (section
  6.5): 0.7% against 29% on the fluid sheet. `patch` resampling, the
  default with flips, restores it exactly, at the price of an instant
  redistribution over the patch at each accepted flip; at a slow `D` that
  is a kinetic shortcut, and `none` is the setting to study it with.
- Flips add a bookkeeping diffusion to the field (section 5.1); not yet
  measured (section 5.2, item 7).
- Found on the way: with the area term on, the legacy `Periodic` mode's
  spectrum reads a tension of 4.1 pN/nm at an internal tension of 1.37,
  where the `Mixed` periodic mode reads 1.42 at 1.68. The legacy seam, not
  the area constraint, inflates the fixed-box numbers quoted in
  `docs/fluctuation_spectrum.md` section 5 and `docs/elastic_box.md`
  section 6.2; under `Mixed` the elastic box at zero tension gives slope
  −4.15 and a fitted tension of −0.11 pN/nm.
