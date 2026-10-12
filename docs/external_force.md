# An external force on a region of the membrane (`externalForceEnabled`)

**Status:** implemented, tested (`tests/test_external_force.cpp`)
**Base:** `feature/mixed-boundary-conditions @ 70de8ec`
**Scope:** a constant force per vertex on a disc or rectangle of the membrane,
switched on and off by iteration, in the Brownian dynamics
(`membrane_dynamics`). Off by default, and a run with it off is the run it
always was: the shipped periodic workload is byte-identical (section 5).
Worked example: `data/example/push_release.params`, animated by
`analysis/push_release_animation.py` (section 4).

## 1. What it is

Every independent vertex whose initial `(x, y)` lies inside the region
receives the force `F = (externalForceX, externalForceY, externalForceZ)`, in
pN, at every Brownian step from `externalForceStartIteration` up to but not
including `externalForceStopIteration`. Fixed vertices, ghosts, periodic
images and the global periodic mode's slaved duplicates are never pushed: the
step does not integrate them, so a force on them would go nowhere.

The force is per vertex, not per area. On the generated lattice each vertex
stands for `(sqrt(3)/2) lFace^2` of membrane -- 21.65 nm^2 at `lFace = 5` --
so a uniform pressure `p` over the region is a force of about `21.65 p` pN per
vertex, and the 5 pN of the example is about 0.23 pN/nm^2.

The region is chosen once, when the model is built, from the vertices'
positions then. A solid sheet's plane does not move (in-plane motion is off),
so a vertex is inside or outside for the whole run. The chosen indices are
written to `<prefix>external_force_vertices.csv`, one per line, and the
startup log says how many were chosen and for which iterations. A region that
selects no vertex is refused at setup, as is a shape other than `disc` or
`rectangle` or a non-positive size: a mis-set centre is not a run worth
starting.

## 2. Parameters

```text
externalForceEnabled = true
externalForceX = 0.0                 # pN per vertex
externalForceY = 0.0                 # pN per vertex
externalForceZ = 5.0                 # pN per vertex; positive lifts the membrane
externalForceShape = disc            # disc | rectangle
externalForceCenterX = 0.0           # nm; the generated sheet is centred on the origin
externalForceCenterY = 0.0
externalForceRadius = 20.0           # nm, for a disc
externalForceHalfWidthX = 20.0       # nm, for a rectangle
externalForceHalfWidthY = 20.0
externalForceStartIteration = 0
externalForceStopIteration = 30000   # exclusive; negative = never switches off
```

`Param` prints them with the rest (`operator<<`), `import_kv_string()` reads
them, and `ExternalForceTest.TheParameterKeysParse` pins the spelling.

## 3. Where the force enters the step

The Brownian step (`DynamicModel::next_step()`) displaces the limit surface
`S = M C`:

```text
    dS = (D dt / kT) F_S + sqrt(2 D dt) xi .
```

The membrane's own force is assembled on the control points, `F_C = -dE/dC`,
and carried onto the surface as `F_S = M^-T F_C` when
`fdtConsistentSurfaceUpdate` is on (see the comment block in `next_step()` for
why). The external force is a force on the membrane surface itself, so it is
added to `F_S` directly, after that map: a vertex inside the region moves by
`(D dt / kT) F` per step on top of whatever the membrane does to it. That is
also what adding its gradient to the nodal force would give -- its energy is
`-F . S_i`, whose control-point gradient is `M^T F`, and `M^-T M^T F = F` --
without a second solve. `ExternalForceTest.AnUpwardForceLiftsThePushedVerticesByTheBrownianMobility`
checks the displacement against `D dt F / kT` to 1e-12 nm, vertex by vertex,
with the noise held fixed by the seed.

Three consequences:

- **With `inPlaneDynamicsEnabled = false` only the z component acts.** The
  step zeroes every in-plane force, and this one is no exception
  (`AnInPlaneForceDoesNothingWhileInPlaneMotionIsOff`).
- **`EnergyForce.csv` holds the membrane's internal energy.** The work the
  force does is not added to it. During a push the bending and
  area-constraint energies rise; that is the membrane storing the work.
- **The Metropolis flip sweep does not see the force.** Its trial energy is
  the membrane's Hamiltonian, which has no `-F . S` term. A fluid run under
  the force is therefore driven consistently but not sampled consistently at
  the pushed vertices; a solid sheet, which is what this was built for, is
  unaffected. See section 6.

Outside its window the pushed run takes the reference run's step to the bit
(`TheForceActsOnlyInsideItsWindow`), and the window is
`[externalForceStartIteration, externalForceStopIteration)` in the driver's
iteration count, so a frame written at iteration `k` is the configuration
after `k` steps and the last pushed step is `stop - 1`.

## 4. The worked example: push a disc, then let go

`data/example/push_release.params`, measured 2026-09-23 with the serial
Release build (`cmake -S . -B build-test -DCMAKE_BUILD_TYPE=Release
-DSLIMED_ENABLE_OPENMP=OFF`), seed 20260923:

| | |
| --- | --- |
| membrane | 100 x 100 nm, `lFace = 5`, 525 vertices, `boundaryType = Mixed` with both axes `Fixed`, two clamped rings (168 fixed, 357 free) |
| physics | `kc = 83.4 pN.nm` (20 kT), `us = 250 pN/nm` global area constraint, `KBT = 4.17`, `D = 1 nm^2/us`, thermal noise on |
| force | 5 pN upward on the 60 vertices of a disc of radius 20 nm at the centre (300 pN, about 0.23 pN/nm^2), on for 0-30 us |
| run | `timeStep = 0.001`, 130 000 steps (30 us pushed, 100 us recovering), a frame every 100 steps; 279 s, 466 steps/s |

What happens (`analysis/push_release_animation.py` prints these):

| quantity | value |
| --- | --- |
| apex height, maximum | 10.8 nm at 23.6 us |
| apex / mean disc height at release (30 us) | 10.4 / 8.1 nm |
| area-constraint energy at release | 615 pN.nm, i.e. 2.2 % excess area, a tension of about 5.4 pN/nm |
| mean disc height below 50 % / 25 % / 10 % of its release value | 12.9 / 43 / 70 us after release |
| apex / mean disc height at the end (100 us after release) | 1.4 / -0.4 nm: the thermal baseline |
| bending energy, whole run | 680 +- 60 pN.nm, flat |

Reading it:

- **The disc rises for about 20 us and then stops rising.** It starts at the
  free-drift rate `(D / kT) F = 1.2 nm/us` and slows as the membrane pushes
  back, reaching 90 % of its final height by about 15 us and a plateau by 20-25
  us. The 4 and 6 pN calibration runs plateau at 9.9 and 11.5 nm (apex), so
  the height is close to linear in the force in this range.
- **The restoring force is mostly tension, not bending.** Lifting the disc 8
  nm adds 2.2 % to the sheet's area, and the global area constraint turns that
  into a tension of about 5 pN/nm, which at the bump's curvature is a much
  larger force than the bending term. The bending energy never leaves its
  thermal level: 357 free height coordinates at `kT/2` each is 744 pN.nm, and
  the run sits at 680 because two clamped rings stiffen the edge modes. A 10
  nm bump on a 100 nm sheet is invisible in `E_bending`; watch `E_area`.
- **Recovery is fast, then slow.** After release the tension pulls the bump
  down with a rate proportional to the excess area, which falls as the height
  squared: the first half is gone in 13 us, the next quarter takes 30 more.
  Once the excess area is thermal the bending term does the rest, and for the
  fundamental mode of a 90 nm clamped patch its rate is
  `(D / kT) A_vertex kc q^4` with `q = pi / 45 nm^-1`, a time constant near
  100 us. That is why the last nanometre takes as long as the first seven, and
  why the disc ends at the thermal baseline only after about 100 us.

The animation (`push_release.gif`, 15 s) shows the limit surface coloured by
height beside the apex trace with the force window shaded;
`push_release_height.png` has the height and energy traces, and
`push_release_snapshots.png` six stills. Run it on any push-release run
directory:

```bash
~/anaconda3/bin/python analysis/push_release_animation.py ~/slimed_runs/push_release
```

## 5. What did not change

- With `externalForceEnabled = false` (the default) `next_step()` evaluates one
  boolean and nothing else, so no run without the key changes. Checked as the
  legacy byte-identity gate prescribes: the shipped periodic workload (100 nm
  sheet, 300 steps, harmonic bonds off) run with the binary built before and
  after this change produces byte-identical `EnergyForce.csv`, `inputface.csv`,
  `inputvertex_begin.csv`, `inputvertex_final.csv`,
  `inputvertex_type_begin.csv`, `inputvertex_type_final.csv`,
  `meshpointinput.csv` and `surfacepointinput.csv`.
- The full test suite passes (215 of 216, the CUDA backend test skipped for
  want of a GPU), the eight new tests included.

## 6. Limits

- **A solid sheet is what this was built and measured for.** Under
  `edgeFlipEnabled` the flip sweep samples the membrane's Hamiltonian without
  the force's potential (section 3), and with `inPlaneDynamicsEnabled` the
  region, chosen from the initial positions, no longer follows the material.
  Both would need the region re-evaluated per step and the potential added to
  the sweep's trial energy.
- **Per-vertex, not per-area.** On a mesh whose vertex density varies the same
  force per vertex is not a uniform pressure.
- **No ramp.** The force switches on and off in one step. The membrane's
  response time is tens of microseconds, so a step is what the example wants;
  a ramp would be a two-line addition to `external_force_active()`.
