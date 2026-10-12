# A rigid bead pushed through the membrane (`beadEnabled`), and what fluidity buys

**Status:** implemented, tested (`tests/test_bead.cpp`); benchmark run 2026-09-23
**Base:** `feature/mixed-boundary-conditions @ 70de8ec`, on top of the external
force of `docs/external_force.md`
**Scope:** a rigid sphere the membrane cannot enter, moving at constant
velocity in the Brownian dynamics; and a benchmark that pushes it slowly
through a clamped sheet twice -- once with Metropolis edge flips on (a fluid
membrane) and once with the same tethered mesh and flips off (a frozen one) --
to see whether in-plane fluidity is what lets a triangulated membrane follow a
large deformation with a bottleneck in it. Off by default; the shipped
periodic workload is byte-identical with it off (section 6).

## 1. The bead

A sphere of radius `R` centred at `c` pushes every independent vertex inside
it out radially with a harmonic wall, evaluated on the control points:

```text
    E_bead = (k / 2) sum over vertices  max(0, R - |x_i - c|)^2 ,
    F_i    = k (R - |x_i - c|) (x_i - c) / |x_i - c|          for |x_i - c| < R .
```

Fixed vertices, ghosts, periodic images and the global periodic mode's slaved
duplicates are not integrated and are not pushed. The force is a nodal force
and reaches the Brownian step through the same map as the membrane's own,
`M^-T` under `fdtConsistentSurfaceUpdate`, added in
`DynamicModel::next_step()` before the map (`src/dynamics/Bead.cpp` holds
the evaluation, the centre's motion and the log).

**Why the control points, and why that keeps the flip sweep exact.** The
Metropolis flip trial differences the membrane's Hamiltonian -- bending, area,
tether, altitude floor, crease wall -- between the two triangulations at fixed
coordinates (`Mesh::evaluate_face_subset()`). A potential that depends on
vertex positions alone is the same before and after a flip, so its difference
is exactly zero and leaving it out of the trial changes nothing: the sweep
samples the membrane's Hamiltonian *plus* the wall, which is what the dynamics
integrates. Had the wall been put on the limit-surface points `S = M C`
instead, `M` would change with the connectivity and the sweep would be
differencing a term it never saw. (The external force of
`docs/external_force.md` acts on the surface points, which is why its own
section 3 lists the flip sweep as a limit; a constant force on the *control*
points would have been consistent too, and is the better choice if that
feature is ever used with flips.)

The wall's stiffness sets how far the membrane sinks into the bead under
load: at `k = 200 pN/nm` a 50 pN reaction on ten vertices is a quarter of a
nanometre. Its relaxation rate in the Brownian step is `k D dt / kT`, 0.05
per step at the benchmark's settings, so a vertex that finds itself inside
the wall is out within a few dozen steps.

**Adhesion.** Outside the wall the bead can carry a well, a smooth step of
depth `w = beadAdhesionEnergy` per vertex over `delta = beadAdhesionRange`:

```text
    E_adh = -w s(x),   x = (d - R) / delta,   s(x) = 1 - 3 x^2 + 2 x^3   on [0, 1],
```

`s = 1` inside the sphere and `0` beyond the range, so a vertex on the bead
is held at the wall by a pull that peaks at `1.5 w / delta` and fades to
nothing `delta` away. It is a function of the control points like the wall,
so the flip sweep stays exact with it. Zero by default. Wrapping a sphere
costs bending of `2 kc / R^2` per area plus the tension; with the lattice's
21.65 nm^2 per vertex, a 10 nm bead at `kc = 83.4` and a tension of 0.5
pN/nm needs `w` above about 47 pN.nm (11 kT) to be wrapped rather than
tented over. The benchmark uses 60.

The centre moves as `c = c0 + v (t - t_start)` between `beadStartIteration`
and `beadStopIteration` and stays put outside that window, so a run can
approach, push through, and then hold; a negative stop means it never stops
moving. The wall and the well act for the whole run wherever the bead is.

## 2. Parameters

```text
beadEnabled = true
beadRadius = 10.0                # R, nm
beadStiffness = 200.0            # k, pN/nm
beadCenterX = 0.0                # c0, nm
beadCenterY = 0.0
beadCenterZ = -5.0
beadVelocityX = 0.0              # v, nm/us
beadVelocityY = 0.0
beadVelocityZ = 0.5
beadStartIteration = 0
beadStopIteration = 70000        # exclusive; negative = never stops moving
beadAdhesionEnergy = 60.0        # w, pN.nm per vertex; 0 = a bead that only pushes
beadAdhesionRange = 2.5          # delta, nm
```

`Param` prints them, `import_kv_string()` reads them, and
`BeadTest.TheParameterKeysParse` pins the spelling. A bead with a radius or a
stiffness that is not positive is refused at setup.

## 3. The bead log

Beside every trajectory frame the dynamics driver appends a row to
`<prefix>bead.csv`, the first for the configuration before any step:

```text
iteration,time_us,centerX,centerY,centerZ,energy_pNnm,contacts,maxPenetration_nm,forceOnBeadX,forceOnBeadY,forceOnBeadZ,adhesion_pNnm,adhered
```

`contacts` is the number of vertices inside the wall, `maxPenetration_nm`
the deepest `R - d` among them, `forceOnBead*` the total force the membrane
exerts on the bead, `-sum F_i` -- the quantity a bead-pushing experiment
measures -- `energy_pNnm` the wall's energy, `adhesion_pNnm` the well's, and
`adhered` the vertices inside the sphere or its range: each stands for one
lattice cell of membrane, so `adhered * 21.65 nm^2 / (4 pi R^2)` is the
wrapped fraction the analysis plots. The row is evaluated on the frame just written, at the
centre the step used, so the frame and the row describe the same
configuration.

## 4. The benchmark

`data/example/bead_push_fluid.params` and `bead_push_frozen.params`, which
differ in one line (`edgeFlipEnabled`), run by
`~/anaconda3/bin/python analysis/bead_push_analysis.py --fluid DIR --frozen DIR`.

| | |
| --- | --- |
| membrane | 100 x 100 nm, `lFace = 5`, 525 vertices, `boundaryType = Mixed` with both axes `Fixed`, one clamped ring (88 fixed, 437 free) |
| physics | `kc = 83.4 pN.nm`, `us = 5 pN/nm` (see below), `KBT = 4.17`, `D = 1 nm^2/us`, thermal noise on |
| fluid mesh terms, both runs | in-plane motion; flat tether `[0.7, 1.75] lFace` at 83.4 pN/nm; altitude floor 0.4 lFace; crease wall from 60 degrees at 500 pN.nm; `surfaceSolver = iterative` |
| flips, fluid run only | `nu = 0.5` attempts per edge per us, valences 4-8, `irregularPatchDepthScale = 1` |
| bead | `R = 10 nm` (two edges), wall 200 pN/nm, well 60 pN.nm (14 kT) over 2.5 nm, from 5 nm below the sheet at 0.5 nm/us for 70 us -- 35 nm, to 20 nm above the sheet -- then held for 20 us |
| run | `timeStep = 0.001`, 90 000 steps, a frame every 100; seed 20260923 in both |

The mesh terms are the `fluid_e` configuration of `docs/fluidity_results.md`
section 11, the one measured to be fluid and flat at conserved area: a lower
tether wall of 0.7 (a fluid sheet at 0.95 is a crystal whose flips are bound
dislocation pairs), the altitude floor and the crease wall (without which
every long fluid run folded), and the full valence range the crease wall
makes safe. The frozen run keeps all of that and only stops the sweep, so it
is a tethered network -- vertices free to move in three dimensions, edges held
in the tether's range, triangulation fixed -- and whatever separates the two
runs is the connectivity's freedom to change.

**Why the area constraint is soft.** A clamped sheet has no reservoir of
area, and wrapping a 10 nm bead takes `4 pi R^2 = 1257 nm^2`, 13 % of the
sheet. At the bilayer's stretching modulus of 250 pN/nm that is a tension of
30 pN/nm, several times what tears a bilayer; the membrane would refuse the
wrap and the bead would simply tent it. A real membrane draws the area from
folds, caveolae or the rest of the cell, at a tension of a fraction of a
pN/nm. The soft constraint stands in for that reservoir: at 5 pN/nm the full
wrap costs a tension near 0.7 pN/nm, and the neck a tension sets,
`sqrt(kc / 2 sigma)`, is then 6-8 nm -- a bottleneck the 5 nm mesh can just
resolve, and where fixed connectivity has to stretch or crumple to follow the
shape. Calibrated first without the bead in contact (12 us): the flat fluid
sheet at 5 pN/nm keeps a thermal roughness of 0.6-0.9 nm r.m.s., smallest
altitudes of 1.4-1.6 nm with 3-4 % of faces on the floor, creases under the
60-degree wall, mean edge 5.4 nm -- the same health as `fluid_e` at 250.

**Why the bead adheres.** The first version of this benchmark had no well,
and its frozen control (`~/slimed_runs/bead_frozen_noadhesion`, complete)
showed why one is needed: the membrane never wrapped. It tented -- apex 27 nm
above the sheet with the bead centre at 20, two to seven vertices touching
the sphere, never a waist narrower than the bead -- because at a tension
near 0.5 pN/nm a cone over the bead costs less than the `8 pi kc` of a wrap.
A cone is developable: a flat sheet becomes one without in-plane strain, and
the frozen mesh followed it with its edges inside the tether. Wrapping a
sphere, or forming the neck behind it, is not developable -- it needs the
in-plane rearrangement that only flips (or stretching) can supply -- which
is exactly the deformation this benchmark is meant to test. Adhesion is how
a bead gets wrapped in every nanoparticle-uptake model, and it is what makes
the neck form here.

**Why 0.5 nm/us is slow.** The membrane's own time scales at the bead's size
are the bending relaxation of a 10 nm mode, `1 / ((D/kT) A_vertex kc q^4)`
with `q = 1/R`, about 20 us, and the tension-driven relaxation of the sheet,
about 10 us. The bead takes 20 us to move one radius, and at `nu = 0.5` with
40 % acceptance every edge flips about four times in that span: the mesh has
time to rearrange under the bead. Slower would be better and costs
proportionally more; the fluid run is the expensive one, at 80 steps/s at 8
threads from the flat start and less as the mesh disorders.

## 5. Results

Both runs completed without a divergence: the fluid one in 889 s at 8
threads (101 steps/s), the frozen one in 205 s at 2 threads. Same seed, same
noise stream, same bead. The two membranes are indistinguishable until the
bead's centre reaches the sheet, and then they part.

| | fluid (flips on) | frozen (flips off) |
| --- | --- | --- |
| wrapped fraction, bead centre at the sheet (30 us) | 0.81 | 0.74 |
| wrapped fraction, one radius up (50 us) | 0.90 | 0.81 |
| wrapped fraction, two radii up (70 us) | 0.93 | 0.72 |
| wrapped fraction at the end of the hold (90 us) | 1.10 (64 vertices) | 0.74 (43 vertices) |
| waist behind the bead, one radius up | 12.7 nm | 13.0 nm |
| waist, two radii up; over the hold (min, mean) | 8.6 nm; 8.4, 9.4 nm | 9.4 nm; 9.4, 11.7 nm |
| force on the bead, mean over the hold | -107 +- 6 pN | -136 +- 4 pN |
| excess area and tension at the end | 9.9 %, 0.49 pN/nm | 8.9 %, 0.45 pN/nm |
| bending / adhesion energy at the end | 2219 / -3682 pN.nm | 1753 / -2442 pN.nm |
| smallest altitude over the run; faces under the floor (mean, max) | 1.02 nm; 4.5 %, 7.8 % | 1.13 nm; 2.7 %, 5.1 % |
| edge extremes over the run (walls at 3.5 and 8.75) | 2.31, 9.74 nm | 2.32, 9.67 nm |
| largest crease; creases over 90 degrees | 76 deg; 0 | 81 deg; 0 |
| valences at the end, interior vertices | 4: 8 %, 5: 24 %, 6: 39 %, 7: 21 %, 8: 9 % | all 6 |
| accepted flips | 9903 of 40 010 admissible (25 %), 110 per us | none |

Reading it:

- **The frozen membrane wraps to the equator and stops.** Its wrapped
  fraction rises with the fluid one's to 0.8 while the bead is still below
  the sheet -- the cap is a shallow dome the lattice can follow by stretching
  -- then stalls at 0.7-0.8 and never closes below the equator. In the
  cross-sections it leaves the sphere at the height of the centre and flares
  into a funnel; the waist behind the bead touches the bead's radius once
  (9.4 nm, at 70 us) and opens back to 12-13 nm during the hold. A sphere is
  not developable: past the equator the sheet would have to shrink its
  circumference, and a lattice whose every vertex keeps its six neighbours has
  no vertices to give up.
- **The fluid membrane wraps the whole bead and closes a neck behind it.**
  The wrapped fraction passes 0.9 one radius up and reaches 1.1 by the end of
  the hold -- 64 vertices on a sphere that holds 58 lattice cells, the extra
  being the neck. The waist falls below the bead's radius as soon as the bead
  is one radius clear of the sheet and settles at 8.4-9.4 nm through the hold,
  against the 9.2 nm the tension sets, `sqrt(kc / 2 sigma)` at
  `sigma = 0.49 pN/nm`. That is the bottleneck: a membrane tube narrower than
  the bead it wraps, which this triangulation can only form by rearranging.
- **The rearrangement is not local.** Accepted flips fall inside two radii
  of the bead's axis 15 % of the time, and that disc is 15.5 % of the free
  area: the flip rate under the bead is the background rate, before, during
  and after the wrap. Fluidity here is a property of the whole sheet --
  vertices flow into the neck through a mesh that is remixing everywhere at
  110 flips per us -- not a response localised to the deformation, and the
  valence distribution at the end (61 % of interior vertices off six) is
  `fluid_e`'s equilibrium one.
- **The wrap is paid for in bending and adhesion, not in mesh quality.** The
  fluid membrane's bending energy is 2219 pN.nm against the frozen 1753 -- it
  has more of itself curved around the sphere and the neck -- and its adhesion
  energy is 50 % deeper. The mesh-quality terms sit in their thermal bands for
  both: smallest altitudes near 1 nm at the extreme, 3-8 % of faces on the
  floor (the fluid mesh somewhat more, as in every conserved-area fluid run),
  edge extremes just past both tether walls, which are soft, and no crease
  over 90 degrees anywhere. The frozen mesh is not visibly strained by these
  measures; it simply cannot get further.
- **The force on the bead.** Both membranes pull the bead down -- adhesion
  and tension want it back toward the sheet -- with about 170 pN while it is
  driven through its last radius, and 107 (fluid) against 136 pN (frozen)
  once it holds: the fully wrapped, necked membrane is nearer its rest.

The outputs are in `~/slimed_runs/bead_compare/`: `bead_push.gif` (both
membranes with the bead, beside the force trace), `bead_push_summary.png`,
`bead_push_profiles.png` and `bead_push_snapshots.png`, all from
`analysis/bead_push_analysis.py`.

## 6. What did not change

- With `beadEnabled = false` (the default) `next_step()` evaluates one
  boolean and nothing else. The legacy byte-identity gate passes: the shipped
  periodic workload (100 nm sheet, 300 steps, harmonic bonds off) gives
  byte-identical output files with the binary built before this work and the
  binary that has both the external force and the bead.
- The full test suite passes (221 of 222, the CUDA backend test skipped for
  want of a GPU), the six bead tests included: the keys, the centre's motion
  and window, the wall's force and energy vertex by vertex against the
  definition and against a central difference, the step lifting every contact
  and holding the membrane within half a nanometre of the wall over 400
  steps, refusal of a bad bead, and the log's format.


## 7. Limits

- **One bead, moving on a straight line.** Several beads, a rotating or
  force-driven bead (one that moves at the force the membrane exerts rather
  than at a prescribed velocity) are the extensions a wrapping study would
  want next; the log already records the force such a bead would feel.
- **No topology change.** A neck cannot pinch off, so a bead pushed far enough
  pulls a tube; that is the model, not a failure, and it is where the
  benchmark stops.
- **The wall is on the control points.** The limit surface of a convex
  control net lies inside it, so the membrane proper sits a fraction of a
  nanometre inside the sphere where it wraps -- cosmetic at `R = 2 lFace`,
  and the price of a potential the flip sweep can ignore exactly.

## 8. The local fields: energy and force maps (`localOutput`)

The runs above write only coordinates and the total energies. With
`localOutput = true` the dynamics also writes, beside every trajectory
frame, what the model knows locally (`Param::localOutput`):

```text
faceenergy<prefix>.csv     one line per frame; per face, in face order:
                           bending energy, mesh-term energy (tether +
                           altitude floor + crease wall), mean curvature,
                           patch area            (pN.nm, pN.nm, 1/nm, nm^2)
vertexforce<prefix>.csv    one line per frame; per vertex, in vertex order:
                           the bending, area-constraint, mesh-term and bead
                           forces, three components each               (pN)
```

The rows are the values the last force evaluation stored, one step before
the frame -- invisible at a frame every hundred steps, and it keeps the
per-face rows on the connectivity the frame's face file describes. Both bead
examples set it. The area constraint is global (`isGlobalConstraint`), so
its energy has no per-face share and the tension it exerts is uniform over
the membrane; its force per vertex is in the second block. Tested in
`tests/test_local_output.cpp` against the Face and Vertex objects and the
recorded totals. Off by default, and the byte-identity gate is unaffected.

`analysis/local_energy_heatmap.py --fluid DIR --frozen DIR` paints three
fields on the membrane and averages them in rings around the bead's axis:

- **bending energy density**, the face's bending energy over its patch
  area, in pN/nm (a flat sheet at `kc = 83.4` sits at its thermal value; a
  sphere of radius 10 nm costs `2 kc / R^2 = 1.7 pN/nm`);
- **in-plane area strain**, the face's patch area against the flat
  lattice cell, positive where the membrane is stretched and negative where
  it is compressed -- the map of what fixed connectivity has to absorb;
- **force magnitude**, the net force on each vertex (bending,
  area-constraint, mesh-term and bead forces summed), which in an overdamped
  model is what drives the vertex, so it is large only where the membrane is
  being pushed off equilibrium.

It writes `local_maps.png` (the three fields on both membranes at one time),
`local_profiles.png` (ring averages at three times) and `local.gif` (the
first two fields animated, both membranes side by side).

Measured on the two runs (both re-run with `localOutput = true`; the
re-run reproduced the first run exactly, every flip and every number of
section 5 included), at the end of the hold, in ring means around the axis
-- "near" is within 1.5 radii, "far" beyond 3:

| | fluid (flips on) | frozen (flips off) |
| --- | --- | --- |
| bending energy density, near / far | 1.09 / 0.13 pN/nm | 0.75 / 0.12 pN/nm |
| bending energy density, largest face | 3.8 pN/nm | 3.2 pN/nm |
| area strain, near mean; range on a face | -6.5 %; -49 % to +113 % | +0.1 %; -34 % to +47 % |
| area strain, far mean | +10.3 % | +8.8 % |
| mesh-term energy density, near mean / largest face | 0.076 / 2.2 pN/nm | 0.029 / 1.4 pN/nm |
| force magnitude, near / far mean; largest | 24 / 12 pN; 81 pN | 20 / 12 pN; 70 pN |

Reading the maps and profiles (`local_maps.png`, `local_profiles.png`,
`local.gif`):

- **Bending energy lives on the sphere and nowhere else.** Both membranes
  put 2-2.8 pN/nm on the cap -- a sphere of radius 10 nm costs
  `2 kc / R^2 = 1.7 pN/nm` plus the rim -- against a far field of 0.12-0.13
  pN/nm, which is the thermal value of a 5 nm mesh at `kc = 83.4` (every
  face holds its `kT/2` or so). The fluid membrane's ring inside the bead's
  radius is higher at the end (2.7 against 1.9 pN/nm) because its neck is
  there too: the frozen one has no neck to bend.
- **The strain map means two different things on the two membranes.** On
  the frozen network it is an elastic strain, and it shows what fixed
  connectivity has to absorb: at first contact (30 us) a ring compressed by
  20-25 % at the bead's radius and a ring stretched by 15-20 % just outside
  it, relaxing to -10 % / +20 % by the end as the funnel opens. On the
  fluid membrane vertices flow, so the same quantity is a density: adhesion
  pulls vertices onto the bead until the tether's lower wall says no (64 on
  a sphere that holds 58 lattice cells: the -6 % near the axis) and leaves
  the neck they came through sparse (+50 % at 1.5 radii), which the tether
  allows and costs no elastic energy. The far field of both sits at the
  +9-10 % the soft area constraint lets the sheet expand to.
- **The mesh terms are quiet except in the neck.** Their energy density is
  2-3 times higher near the axis on the fluid membrane than on the frozen
  one -- the neck's sparse faces sit against the tether's upper wall and
  its altitude floor -- and is a tenth of the bending density either way.
- **Force is where the bead is.** The net force on a vertex, which drives
  it, is 20-24 pN in the wrap against 12 pN of thermal noise everywhere
  else, on both membranes; the largest values are at the neck (fluid) and
  at the funnel's rim (frozen).


## 9. The sticky bead on a periodic sheet (`beadPotential = harmonic`)

The wall of section 1 only pushes, and the well of section 5 holds the
membrane on with a force that saturates at `1.5 w / delta`. A *sticky* bead
is the other limit: once the membrane has touched it, it is bonded, and the
bond is the full harmonic potential on both sides of the surface,

    E = (k / 2) sum over bonded vertices (|x_i - c| - R)^2,   k = beadStiffness.

A free vertex is bonded when it comes within `beadCaptureRange` of the
sphere's surface, and a bond breaks when its vertex is farther than
`beadBreakRange` from the surface -- never, when that is negative (the
default), so the bead keeps what it touched. The two ranges give the bond
hysteresis, so a vertex sitting at the capture distance does not flicker.
The bonded set is state the step carries (`DynamicModel::beadBonded`),
changed only in `refresh_bead_contact()`, which a flip never calls: like the
wall the energy is a function of the control points alone and the Metropolis
sweep stays exact. The bead log gains two columns, `bonded` and
`maxStretch_nm` (the largest `d - R` among the bonded vertices outside the
sphere). `beadAdhesionEnergy` is ignored under `harmonic`.
`tests/test_bead.cpp` checks the keys, the force on both sides of the surface
against `-dE/dz`, the capture, the persistence past the capture range, and
the break.

**Choose `k` from the capture energy, not from the wall.** At the wall's
200 pN/nm a vertex captured 2.5 nm out is pulled with 500 pN, and a mesh that
cannot rearrange is folded by it: with flips off the energy went non-finite
at 53-82 us (reversed face normals from 47 us on, as the bead left the plane
and the sheet had to neck), under every reservoir. `k = 20 pN/nm` makes the
capture worth `k delta^2 / 2 = 62 pN.nm = 15 kT` per vertex -- the well of
section 5 -- at a tenth of the force, and with 30-60 vertices sharing the
bead's reaction the wall side still penetrates well under a nanometre.

**A periodic sheet needs a pin.** Nothing holds a sheet periodic along both
axes in z: a bead that pulls on it moves the whole tile at the total force
over `N kT / D` (docs/elastic_box.md section 5, "unheld periodic sheet"), and
a push of fixed length becomes a push of unknown depth. `meanHeightPinned =
true` subtracts the mean vertical displacement of the integrated vertices
from every one of them at each step: a uniform counter-force, which is what
the surrounding membrane the tile stands for exerts, so the shape and the
fluctuations are untouched and only the drift is removed (one Brownian mode
of the tile). Off by default; tested on the periodic sheet with and without a
bead.

### 9.1 The matrix (`analysis/sticky_bead_runs.py`)

A generated 130 nm sheet, which under `Mixed`/`Periodic` is a 100 x 104 nm
tile of 480 free vertices; `kc = 83.4`, `K = 250` (the bilayer's modulus:
the reservoir supplies the area, not stretch), the fluid mesh terms of
section 4, the bead of `bead_push_fluid.params` with `beadPotential =
harmonic`, `k = 20`, capture 2.5 nm, rising 35 nm in 70 us and holding 110 us
(180 000 steps at `dt = 0.001`), frames every 0.2 us with `localOutput`.

- Reservoir: **box** -- the elastic box at `tau = 0.5 pN/nm`, mobility 0.02
  (docs/elastic_box.md); **s0global** -- `lipidFlowMode = global`, bulk
  reservoir at 0.5 pN/nm, rate 20/us, thermal; **s0local** -- the per-face
  relaxed area with tension-driven flow at the natural `D_sigma` and the same
  bulk reservoir; **fixed** -- neither, one run as the control.
- Fluidity: flips off (frozen), and `edgeFlipAttemptRate` 0.05, 0.5 and
  5 per edge per us.

13 runs, made on the cluster (`submit_array.sh`, 4 threads each) with the
`-O3` OpenMP build; `analysis/sticky_bead_analysis.py --root DIR` writes the
traces, the hold-phase relaxation fits, the bending-energy maps and ring
profiles, the Fourier analysis of the hold phase, one animation per run and
`stats.md`.

### 9.2 Results

(filled in below from the runs of 2026-10-12)
