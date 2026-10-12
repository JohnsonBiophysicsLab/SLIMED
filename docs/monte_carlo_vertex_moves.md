# Monte Carlo vertex moves and hard walls (`vertexMoveEnabled`, `meshWallType`)

**Status:** implemented, tested (`tests/test_vertex_move.cpp`); validation runs
started 2026-10-01 and not yet analysed (section 9)
**Base:** `feature/mixed-boundary-conditions @ 70de8ec`, on top of the edge
flip (`docs/edge_flip_plan.md`) and the fluid mesh-quality terms
(`docs/fluidity_results.md`)
**Scope:** two options for `membrane_dynamics`. `vertexMoveEnabled` adds a
Metropolis sweep of single-vertex moves on the control points, beside the
Brownian step or instead of it (`brownianStepEnabled = false`).
`meshWallType = hard` turns the three fluid mesh-quality walls -- tether
range, altitude floor, crease angle -- from quadratic energies into
constraints that the Monte Carlo moves enforce by rejection. Both are off by
default, and a run with them off is unchanged to the byte (section 8).
Worked example: `data/example/monte_carlo_hard_wall.params`.

## 1. Why

A fluid SLIMED run keeps its triangles healthy with three soft walls
(`Edge_spring.cpp`, `Triangle_shape.cpp`, `Crease_wall.cpp`). They are soft
because the Brownian step is explicit: it cannot refuse a step, so the
Hamiltonian has to push back. That has three costs, all measured:

- **The walls are part of what is sampled.** They are zero inside their
  range, but a fluid run leans on them all the time: 100-470 pN.nm stored per
  frame, the altitude floor live on 5% of faces, the shortest edge 0.8 nm
  inside the lower tether wall (`docs/fluidity_results.md` section 11).
- **Their stiffness caps the time step.** At `timeStep = 0.001` the tether
  and altitude constants can be raised fourfold; at eightfold the run diverges
  after 3 300 steps and at sixteenfold after 180.
- **A regularizer in the drive is a regularizer in the ensemble.** Putting
  the solid mode's reference-length term into the Brownian drive halves the
  bending energy (501 to 275 pN.nm) and reads as a tension of 9 pN/nm.

A Monte Carlo move can refuse. The wall can then be infinitely high, which
means it stores nothing, needs no stiffness and is never leaned on -- the
sampled distribution is `exp(-E/kT)` of the membrane's own energy on the set
of meshes inside the walls. That is how every dynamically triangulated
surface model holds its mesh together (section 2).

What this is *not* is a regularization that leaves the vertex ensemble
untouched. No such thing exists: section 3.3.

## 2. What the triangulated-surface literature does, and what is adopted

The model goes back to the tethered (polymerized) surfaces of Kantor, Kardar
& Nelson [1, 2]: beads joined by tethers of bounded length, hard-core
repulsion between beads, sampled by Metropolis moves of single beads [10].
Making the surface fluid adds the bond flip -- the diagonal of two adjacent
triangles is replaced by the other one -- introduced for dynamically
triangulated random surfaces by Billoire & David [18] and Boulatov et al.
[19], and brought to fluid membranes by Ho & Baumgärtner [3, 4] and Kroll &
Gompper [5]. Reviews: Gompper & Kroll [6, 7]; how the bending energy is
best discretized on such a network is a separate question [11], and not one
SLIMED shares, since its bending energy is evaluated on the limit surface.

| | mesh constraint | vertex move | one Monte Carlo step |
| --- | --- | --- | --- |
| Tether-and-bead model, as reviewed in [6] and stated in [8] | bond potential infinite below the bead diameter `sigma0` and above the tether length `l0`, zero between; `l0 < sqrt(3) sigma0` makes the network self-avoiding | displacement uniform in the cube `[-s, s]^3`, `s` chosen for about 50% acceptance | -- |
| Gompper & Kroll, inflated vesicles [9] | beads of diameter 1, tethers of maximum extension `sqrt(2)`; flip acceptance about 8% | -- | counted in steps per vertex |
| Ramakrishnan, Sunil Kumar & Ipsen [16] | unit beads, maximal tether `sqrt(3)`, plus a mild dihedral-angle constraint | within a cube of side `2 sigma`, `sigma = 0.1`, giving 30-50% acceptance | `N` vertex moves, `3(N - 2)` flip attempts |
| Ramakrishnan, Sunil Kumar & Radhakrishnan [17] | as above | cube of side `2 sigma`, tuned to about 50% | `N` vertex moves and `L` link flips |
| FreeDTS [15] | minimum vertex distance `l_dts`, edges within `l_dts <= l <= sqrt(3) l_dts`, a mild constraint on the dihedral angle between neighbouring triangles | Metropolis vertex update | vertex updates, link flips and inclusion moves in fixed numbers |
| Noguchi & Gompper [12, 13]; TriMem [14] | a *smooth* bond potential, because their vertices are integrated (molecular dynamics; hybrid Monte Carlo) rather than moved one at a time | all vertices together | TriMem: one hybrid Monte Carlo step and a flip sweep over a fraction of the edges |

(The entries for [6] are taken from [8], which cites it for exactly this; the
primary text of [3]-[7] was not available to check, and nothing above is
attributed to those papers beyond what [8], [9] and [16] say of them.)

Two families, then. Models that move one vertex at a time use hard walls and
rejection. Models that integrate all vertices at once need a force, so they
smooth the wall -- which is what SLIMED's Brownian fluid mode does, for the
same reason, and what TriMem says of itself.

**Adopted here**, on SLIMED's own energy (Helfrich bending on the Loop limit
surface, the global area and volume constraints):

- the single-vertex move with a uniform cube proposal and a fixed step,
  exactly as in [6, 8, 16];
- hard walls at the positions the soft walls already have, so the two modes
  describe the same membrane: the tether range `[edgeTetherMinRatio,
  edgeTetherMaxRatio] * lFace` in place of bead diameter and tether length
  (the lower bound plays the bead's role between bonded neighbours), the
  altitude floor, and the crease angle -- the last being the same device as
  the dihedral constraint of [15, 16];
- the flip sweep that already exists, now also rejecting a flip that would
  break a wall.

**Not adopted:** hard-core repulsion between vertices that are *not* joined by
an edge. The bead-and-tether models need it for self-avoidance of a surface
that can fold back on itself; a nearly flat sheet does not get there, and the
bending energy here lives on the limit surface rather than on bead positions.
A closed vesicle or a strongly deformed sheet would need it (section 11).

## 3. The move

### 3.1 Algorithm

One sweep (`Mesh::vertex_move_sweep()`, `src/mesh/Vertex_move_sweep.cpp`)
makes `vertexMoveAttemptsPerVertex * N` attempts, `N` the number of
independent vertices. Each attempt:

1. picks an independent vertex uniformly at random;
2. draws a displacement uniform in the cube `[-s, s]^3`,
   `s = vertexMoveStepRatio * lFace` (along `z` only when
   `inPlaneDynamicsEnabled = false`), and applies it to the vertex's control
   point and to every periodic copy of it;
3. under hard walls, checks the walls around the moved vertices and rejects
   at once if one is broken;
4. otherwise computes the exact energy change `dE` (section 4) and accepts
   with probability `min(1, exp(-dE / kT))`;
5. on rejection, restores the saved coordinates exactly.

The random numbers are counter-based (`Counter_rng.hpp`), keyed on the seed,
the iteration and the attempt, with a tag that keeps the stream clear of the
Brownian noise and the flip sweep. A sweep is reproducible from the seed and
the iteration alone.

### 3.2 Why it samples `exp(-E/kT)`

The vertex is chosen with a probability that does not depend on the
configuration, and the displacement `d` is as likely as `-d`. So the
probability of proposing `x -> x'` equals that of proposing `x' -> x`, and
with the Metropolis acceptance [10] each attempt satisfies detailed balance
with respect to `exp(-E/kT)`. A hard wall is the statement that `E` is
infinite outside it: the acceptance there is zero, which is the rejection in
step 3. The stationary distribution is `exp(-E/kT)` restricted to the
admissible set, provided the run starts inside it -- the run checks that
before the first sweep and refuses to start otherwise.

Composition with the other kernels is the argument `Edge_flip_sweep.cpp`
already makes: the flip sweep, the vertex sweep and the Brownian step each
leave `exp(-E/kT)` invariant, so any sequence of them does. The sequence is
not reversible, which equilibrium sampling does not need [21].

The step is fixed for the whole run. Adapting it to hold an acceptance target
makes the proposal depend on the chain's history and breaks the argument;
tune it in a short run and then leave it.

### 3.3 Why the proposal is not biased toward good triangles

The tempting variant is a move that *proposes* better-shaped triangles and
accepts on the membrane's energy alone, in the hope of regularizing without
touching the ensemble. It does not do that. If the proposal is itself
reversible with respect to some weight `g(x)` -- any stochastic rule that
prefers good triangles is, with `g` larger on good meshes -- then

```text
    g(x) exp(-E(x)/kT) q(x -> x') min(1, exp(-dE/kT))
```

is symmetric in `x` and `x'`, so the chain's stationary distribution is
`g(x) exp(-E/kT)`: the same as adding `-kT ln g` to the Hamiltonian. A
deterministic rule (move the vertex to the centroid of its neighbours) has
no reverse move and satisfies no detailed balance at all. And with the
Hastings correction [20] the chain samples `exp(-E/kT)` exactly, which fixes
the equilibrium distribution of triangle shapes at what `E` alone gives -- the
bias changes the path and nothing else.

So a regularizer either changes the ensemble over vertex positions or does
nothing at equilibrium. The honest options are to change it in a way that can
be stated -- a restriction to an admissible set, which is what a hard wall
is -- or to confine the change to the parametrization (section 11).

## 4. The energy difference

A face's energy is a functional of its control net: its three corners and
their neighbours. A moved vertex therefore changes exactly the faces incident
to its closed one-ring -- 24 on a regular mesh -- and everything else cancels
(`Mesh::vertex_move_patch_faces()`). The trial evaluates that patch at the
trial position and again at the original one, through
`Mesh::evaluate_face_subset()`, the same routine the flip trial uses and the
same kernel as the whole-mesh force pass.

**What is in `dE`** is what the run samples:

| term | fluid mode, soft walls | fluid mode, hard walls | solid mode (`edgeSpringEnabled = false`) |
| --- | --- | --- | --- |
| bending | yes | yes | yes |
| global area and volume constraints | yes | yes | yes |
| tether, altitude floor, crease wall | yes, as energies | no: constraints | -- |
| reference-length regularization | -- | -- | **no**: the Brownian step does not integrate it either |

The two global constraints are quadratic in the mesh-wide totals and enter
through the local change, `dE_A = (uSurf / 2 A0) dA (dA + 2 (A - A0))`,
exactly as for the flip; the sweep keeps `param.area` and `param.vol` current
across accepted moves.

**Periodic copies** move with their source, and the faces around every copy
are part of the patch, so `dE` is the change of the whole-mesh energy as a
function of the independent coordinates. `Mesh::build_vertex_move_table()`
decides who is independent:

- `boundaryType = Mixed`: the mesh's own independent vertices; the followers
  of a source are its periodic images.
- `boundaryType = Periodic`: the `(nFaceX - 6) x (nFaceY - 6)` tile starting
  at index (3, 3); every other vertex follows the tile vertex a whole number
  of periods away. This is *not* the Brownian step's partition: that treats
  both copies of the duplicated ring as slaved to each other, which leaves
  221 of the 252 coordinates integrated directly on the 100 nm sheet and the
  ring moved only through the surface solve. A Monte Carlo move needs every
  coordinate to have one owner, so the tile's copy owns it.
- `boundaryType = Fixed`: every vertex not on a boundary face. `Free` is
  refused; use `Mixed` with free axes.

## 5. Hard walls

`meshWallType = hard` reads the three walls as constraints
(`Mesh::faces_satisfy_walls()`):

| wall | exists when | constraint |
| --- | --- | --- |
| tether | `edgeSpringEnabled` (required) | every membrane edge has `edgeTetherMinRatio * l0 <= l <= edgeTetherMaxRatio * l0` |
| altitude floor | `triangleShapeEnabled` | every altitude `2A / l_i` of every energy-carrying face is at least `triangleShapeMinAltitudeRatio * lFace` |
| crease | `creaseWallEnabled` | across every interior membrane edge, `n1 . n2 >= cos(creaseWallAngle)` |

The positions and the exemptions are the soft terms' own -- an edge joining
two copies, or a ghost face, is not charged in soft mode and not constrained
in hard mode -- and a test pins the two against each other: the soft energy
is zero exactly when no hard wall is broken. The stiffness constants are
unused. The fluid term contributes no energy and no force
(`Mesh::energy_force_fluid_terms()` returns before evaluating it), so
`E_Regularization` in `EnergyForce.csv` is identically zero.

A vertex move can only change the edges, altitudes and creases of faces with
a moved corner, so those are what step 3 checks, before any patch integral is
evaluated. A flip can only change its two new triangles;
`Mesh::evaluate_edge_flip()` checks them and reports a wall-breaking flip as
an attempt with `dE = +inf`, which Metropolis rejects and the flip log
records as `inf`.

Hard walls need `brownianStepEnabled = false`. A hard wall has no force, so
the Brownian step cannot feel it and has no way to refuse a step through it;
once outside, every Monte Carlo move in that neighbourhood is rejected for
good. Setup refuses the combination rather than define it.

## 6. Running it

Pure Monte Carlo with hard walls -- the triangulated-surface sampler -- is
`data/example/monte_carlo_hard_wall.params`:

```text
brownianStepEnabled = false
vertexMoveEnabled = true
vertexMoveStepRatio = 0.1
meshWallType = hard
edgeSpringEnabled = true        # + triangleShapeEnabled, creaseWallEnabled
edgeFlipEnabled = true
timeStep = 1.0                  # attempts per edge per sweep = edgeFlipAttemptRate * timeStep
```

An iteration is then one vertex sweep and one flip sweep. There is no time in
the run: it samples the equilibrium ensemble and says nothing about dynamics.

The move also runs beside the Brownian step with soft walls
(`vertexMoveEnabled = true` alone). The walls are then energies in `dE`, the
two kernels share one Hamiltonian, and the vertex sweep is an extra,
exact relaxation channel on top of the dynamics.

With `inPlaneDynamicsEnabled = false` the move is along `z` only, which makes
it an exact sampler for the solid sheet -- no time-step error -- and the
cleanest check on the Brownian results (section 9).

| parameter | default | meaning |
| --- | --- | --- |
| `vertexMoveEnabled` | `false` | run the vertex sweep |
| `vertexMoveStepRatio` | `0.05` | half-width of the trial cube, in `lFace` |
| `vertexMoveAttemptsPerVertex` | `1.0` | mean attempts per independent vertex per sweep |
| `vertexMoveInterval` | `1` | iterations between sweeps |
| `meshWallType` | `soft` | `soft` or `hard` |
| `brownianStepEnabled` | `true` | take the Brownian step |

Outputs: `<prefix>VertexMoves.csv`, one row per trajectory frame --
`iteration, attempted, wallRejected, accepted, acceptance, deltaEnergy,
tetherOutside, altitudeOutside, creaseOutside` -- the counts since the
previous row and the walls broken on the frame's mesh (always zero under hard
walls; under soft walls, how much the mesh is leaning). The run prints the
totals at the end.

Refused at setup, each with its reason: hard walls with the Brownian step on;
hard walls without the flat tether; no Brownian step and no vertex move;
vertex moves together with the scaffolding terms, the bead, the external
force, the elastic box or the local area constraint (section 11); a starting
mesh that breaks a hard wall.

## 7. Code map

| what | where |
| --- | --- |
| the sweep, the trial, the follower table, the wall checks | `src/mesh/Vertex_move_sweep.cpp` |
| declarations: `VertexMoveDelta`, `VertexMoveSweepStats`, `WallViolations`, `Mesh::vertex_move_sweep()`, `evaluate_vertex_move()`, `faces_satisfy_walls()`, `count_wall_violations()`, `build_vertex_move_table()` | `include/mesh/Mesh.hpp` |
| hard walls carry no energy or force | `Mesh::energy_force_fluid_terms()` (`src/energy_force/Compute_energy_and_force_on_mesh.cpp`), `Mesh::evaluate_face_subset()` (`src/energy_force/Local_patch_energy.cpp`) |
| a flip through a hard wall is rejected | `Mesh::evaluate_edge_flip()` (`src/energy_force/Local_patch_energy.cpp`) |
| configuration checks | `DynamicMesh::check_monte_carlo_configuration()` (`src/dynamics/Dynamic_mesh.cpp`) |
| the loop: optional Brownian step, the sweep, the log | `src/Run_dynamics_flat.cpp` |
| parameters | `include/Parameters.hpp`, `src/io/input.cpp` |

## 8. Tests

`tests/test_vertex_move.cpp`, ten tests:

- **The table.** The candidates are the periodic tile, no vertex follows two
  sources, and a copy stays exactly one period from its source when the
  source moves (both periodic modes).
- **The energy difference.** On a roughened fluid sheet with all three soft
  walls live, the local `dE` equals the change of the whole-mesh energy to
  1e-7 pN.nm, for an interior vertex and for a seam vertex whose copies move
  with it, in both periodic modes; and the trial leaves every coordinate
  bit-identical. In solid mode `dE` is the change of bending plus area, with
  the reference-length term excluded.
- **The walls.** Each wall alone, over sheets roughened from barely to badly:
  the soft energy is non-zero exactly when the hard wall is broken. Under
  hard walls the fluid term's energy and force are zero, a move through a
  wall returns `wallViolated` with the mesh untouched, and a flip that would
  leave an edge above the upper wall is never accepted while the same flip
  inside the walls is.
- **The invariant.** Sweeps of vertex moves and flips from the flat lattice
  never leave the admissible set, with wall rejections actually occurring,
  copies still in place and the running area equal to the recomputed one.
- **The stationary distribution.** One vertex moving along `z`: the sampled
  variance against `exp(-E(z)/kT)` by quadrature, at two temperatures. The
  ratios are 1.026 at `kT = 4.17` and 1.007 at `kT = 16.68` (2 000 samples
  each; the gate is 12%).
- **Reproducibility**, and the refused configurations.

The whole suite: 247 passed, 1 skipped (the CUDA comparison, no device).

**Byte identity.** With both options off, the binary before this change and
the one after produce identical output on the legacy periodic workload
(100 nm sheet, 300 steps: 9 files) and on a fluid workload with flips
(`fluid_e` parameters, 400 steps, 1 064 accepted flips: 51 files).

## 9. Validation

**Not yet analysed.** Three runs on the 100 nm sheet were started and had not
finished when this was written: a solid sheet sampled by `z`-only moves alone
(30 000 sweeps, to compare its bending energy with equipartition, about
523 pN.nm for 251 modes, and its spectrum with the Brownian run's 501 +- 2
pN.nm and `kc` = 83.9), and the fluid sheet under hard walls in both periodic
modes (600 sweeps each). What they showed while running: acceptance 0.79
(solid, trial half-width 0.5 nm) and 0.76 (fluid, 0.3 nm) -- both steps are
on the small side -- about 10% of fluid attempts refused by a wall, and no
wall broken on any frame. The energies and the spectrum are still to be read.

## 10. Cost, and what is traded

**Not yet measured properly.** One attempt evaluates its patch twice, about
48 face integrals on a regular mesh, so a sweep costs on the order of 25
force evaluations on a solid sheet and more on a fluid one, whose irregular
patches are the expensive ones; the trial is serial. Rough rates on a loaded
machine, serial `-O3`, 252 vertices: about 12-30 sweeps/s solid, under
1 sweep/s fluid. A sweep moves a vertex further than a Brownian step does, so
the comparison that matters is per unit of decorrelation, and that has not
been made. Storing the patch's "before" energy instead of re-evaluating it
would halve the cost.

What is traded for the exact walls is the clock: a pure Monte Carlo run has
no dynamics (section 11).

## 11. Not covered

- **Energy terms the trial does not difference**: the scaffolding (harmonic
  bond, Gag, idealized lattice) terms, the bead, the external force, the
  elastic box, the local area constraint. Each is a function of vertex
  positions that would have to be added to `dE`; until it is, setup refuses
  the combination. The bead and the external force are sums over vertices
  and are the easy ones.
- **Self-avoidance between non-bonded vertices.** The walls bound each
  triangle and each pair of adjacent triangles; nothing stops two distant
  parts of a strongly deformed surface from passing through each other. The
  bead-and-tether models get that from the hard-core beads [6, 16].
- **The periodic seam with flips.** Flips are still frozen around the copies,
  in both periodic modes, so the connectivity is periodic only up to one ring
  (`docs/fluidity_results.md` section 12). Under `boundaryType = Periodic`
  the duplicated ring is in addition excluded from the tether, as it is for
  the Brownian step; `Mixed` has no such ring and is the mode to use.
- **Dynamics.** A Monte Carlo sweep is not a time step. Relaxation times,
  diffusion constants and anything else with a clock in it still need the
  Brownian step and its soft walls.
- **Measure factors.** Sampling all `3N` control-point coordinates with a
  flat measure is what the Brownian fluid mode already does, and what the
  triangulated-surface models do; the relation between that measure and the
  one of the continuum Helfrich theory is the subject of [22] and is not
  addressed here.
- **A parametrization-only regularizer.** Tangential moves accepted on the
  membrane's energy would confine the change of ensemble to the in-plane
  parametrization, which for a fluid membrane is nearly a gauge choice. Not
  implemented; it is approximate where the hard wall is exact.
- **Throughput.** The trial evaluates its patch twice and is not threaded
  (section 10).

## 12. References

Bibliographic details checked against Crossref and arXiv on 2026-10-01.

1. Y. Kantor, M. Kardar, D. R. Nelson, "Statistical mechanics of tethered
   surfaces", *Phys. Rev. Lett.* **57**, 791 (1986).
   doi:10.1103/PhysRevLett.57.791
2. Y. Kantor, M. Kardar, D. R. Nelson, "Tethered surfaces: statics and
   dynamics", *Phys. Rev. A* **35**, 3056 (1987). doi:10.1103/PhysRevA.35.3056
3. J.-S. Ho, A. Baumgärtner, "Simulations of fluid self-avoiding membranes",
   *Europhys. Lett.* **12**, 295 (1990). doi:10.1209/0295-5075/12/4/002
4. A. Baumgärtner, J.-S. Ho, "Crumpling of fluid vesicles", *Phys. Rev. A*
   **41**, 5747 (1990). doi:10.1103/PhysRevA.41.5747
5. D. M. Kroll, G. Gompper, "The conformation of fluid membranes: Monte Carlo
   simulations", *Science* **255**, 968 (1992). doi:10.1126/science.1546294
6. G. Gompper, D. M. Kroll, "Network models of fluid, hexatic and polymerized
   membranes", *J. Phys.: Condens. Matter* **9**, 8795 (1997).
   doi:10.1088/0953-8984/9/42/001
7. G. Gompper, D. M. Kroll, "Triangulated-surface models of fluctuating
   membranes", in *Statistical Mechanics of Membranes and Surfaces*, 2nd ed.,
   eds. D. Nelson, T. Piran, S. Weinberg (World Scientific, 2004),
   pp. 359-426. doi:10.1142/9789812565518_0012
8. U. S. Schwarz, G. Gompper, "Bicontinuous surfaces in self-assembling
   amphiphilic systems", in *Morphology of Condensed Matter*, eds. K. R. Mecke,
   D. Stoyan, Lecture Notes in Physics **600** (Springer, 2002), pp. 107-151.
   arXiv:cond-mat/0301325 (section 6.4 for the tether-and-bead model)
9. G. Gompper, D. M. Kroll, "The shape of inflated vesicles", *Phys. Rev. A*
   **46**, 7466 (1992). doi:10.1103/PhysRevA.46.7466, arXiv:hep-lat/9210037
10. N. Metropolis, A. W. Rosenbluth, M. N. Rosenbluth, A. H. Teller,
    E. Teller, "Equation of state calculations by fast computing machines",
    *J. Chem. Phys.* **21**, 1087 (1953). doi:10.1063/1.1699114
11. G. Gompper, D. M. Kroll, "Random surface discretizations and the
    renormalization of the bending rigidity", *J. Phys. I France* **6**, 1305
    (1996). doi:10.1051/jp1:1996246
12. H. Noguchi, G. Gompper, "Fluid vesicles with viscous membranes in shear
    flow", *Phys. Rev. Lett.* **93**, 258102 (2004).
    doi:10.1103/PhysRevLett.93.258102
13. H. Noguchi, G. Gompper, "Dynamics of fluid vesicles in shear flow: effect
    of membrane viscosity and thermal fluctuations", *Phys. Rev. E* **72**,
    011901 (2005). doi:10.1103/PhysRevE.72.011901
14. M. Siggel, S. Kehl, K. Reuter, J. Köfinger, G. Hummer, "TriMem: a
    parallelized hybrid Monte Carlo software for efficient simulations of
    lipid membranes", *J. Chem. Phys.* **157**, 174801 (2022).
    doi:10.1063/5.0101118
15. W. Pezeshkian, J. H. Ipsen, "Mesoscale simulation of biomembranes with
    FreeDTS", *Nat. Commun.* **15**, 548 (2024).
    doi:10.1038/s41467-024-44819-w
16. N. Ramakrishnan, P. B. Sunil Kumar, J. H. Ipsen, "Monte Carlo simulations
    of fluid vesicles with in-plane orientational ordering", *Phys. Rev. E*
    **81**, 041922 (2010). doi:10.1103/PhysRevE.81.041922, arXiv:1004.4509
17. N. Ramakrishnan, P. B. Sunil Kumar, R. Radhakrishnan, "Mesoscale
    computational studies of membrane bilayer remodeling by
    curvature-inducing proteins", *Phys. Rep.* **543**, 1 (2014).
    arXiv:1510.04778
18. A. Billoire, F. David, "Scaling properties of randomly triangulated
    planar random surfaces: a numerical study", *Nucl. Phys. B* **275**, 617
    (1986). doi:10.1016/0550-3213(86)90577-8
19. D. V. Boulatov, V. A. Kazakov, I. K. Kostov, A. A. Migdal, "Analytical
    and numerical study of a model of dynamically triangulated random
    surfaces", *Nucl. Phys. B* **275**, 641 (1986).
    doi:10.1016/0550-3213(86)90578-X
20. W. K. Hastings, "Monte Carlo sampling methods using Markov chains and
    their applications", *Biometrika* **57**, 97 (1970).
    doi:10.1093/biomet/57.1.97
21. V. I. Manousiouthakis, M. W. Deem, "Strict detailed balance is
    unnecessary in Monte Carlo simulation", *J. Chem. Phys.* **110**, 2753
    (1999). doi:10.1063/1.477973
22. W. Cai, T. C. Lubensky, P. Nelson, T. Powers, "Measure factors, tension,
    and correlations of fluid membranes", *J. Phys. II France* **4**, 931
    (1994). doi:10.1051/jp2:1994175
23. P. J. Rossky, J. D. Doll, H. L. Friedman, "Brownian dynamics as smart
    Monte Carlo simulation", *J. Chem. Phys.* **69**, 4628 (1978).
    doi:10.1063/1.436415 (the alternative of a Metropolis test on the
    Brownian step itself; not implemented)
