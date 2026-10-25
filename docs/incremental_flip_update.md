# The incremental post-flip force pass

A fluid run evaluates the energy and force twice on most steps. The first pass
is the one the Brownian step needs. The second follows the Metropolis flip
sweep: an accepted flip retriangulates two faces and changes the control nets
of the faces around them, so the stored energies and forces no longer describe
the mesh, and `run_dynamics_flat()` called `Compute_Energy_And_Force()` again
over every face. On the `box_nu5` benchmark (960 faces, `nu = 5`,
`docs/patch_kernel_performance.md`) 70% of steps accept at least one flip, and
a 60 s `sample` profile of `main` at 9f32f88 put that second pass at 35-37% of
the run: the same share at step 400 and at step 2000.

| phase (`Run_dynamics_flat.cpp` call site) | share of run time |
| --- | --- |
| `:304` force pass for the Brownian step | 52% |
| `:352` force pass after an accepted flip | 35-37% |
| `:315` the flip sweep itself | 9% |
| `:226` the step (solver, noise, box) | under 2% |

Only the faces incident to the four vertices of a flipped edge's quadrilateral
change -- `Mesh::flip_patch_faces()`, about eighteen faces per flip on a
near-regular mesh, the set the trial itself differences -- and the sweep runs
at fixed coordinates. Everything else the second pass recomputed, it recomputed
to the same bits. This change makes the second pass recompute only what moved.

## 1. What changed

### The per-face force store

`Mesh::compute_face_energies_and_forces()` used to integrate each face into
stack arrays and scatter them onto the vertices inside the same loop. It now
integrates each face into its own block of `Mesh::faceForceStore` -- three
`kMaxControlPoints x 3` blocks per face, bending, area and volume, in the order
of `Face::oneRingVertices` -- and a second loop,
`Mesh::scatter_face_forces_to_vertices()`, sums the blocks onto the vertices.
The scatter uses the face loop's `schedule(static, 4)` and the same per-thread
buffers and reduction, so a face lands in the same thread's buffer in the same
order as before, and the vertex forces are a function of the store alone. The
whole-mesh pass is the same arithmetic it was: the legacy byte-identity gate
holds (section 2). The face loop body is `Mesh::evaluate_face_patch()`, the
area pass's is `Mesh::integrate_face_area_volume()`, and the steps after the
face loop are `Mesh::finalize_energy_and_force()`, so that the incremental
pass can call exactly the pieces it needs.

### The constraint factors, applied after the sum

The area force the kernel integrates is
`uSurfPerArea * (area - area0) * dS_f/dx`, with `area` the whole membrane's.
A flip moves `area`, so even a face the flip never touched has a different
area force afterwards -- which would make every stored block stale. When the
incremental mode is active (`Mesh::incremental_flip_forces_active()`: the
parameter on, flips on, the CPU backend, and not local lipid flow) the kernel
is handed unit factors -- `uSurfPerArea = 1, area = 1, area0 = 0` and
`uVol = 3, vol = 1, vol0 = 0`, both exactly one -- so the blocks hold `dS/dx`
and `dV/dx`, and the scatter multiplies the summed gradient by
`uSurfPerArea * (area - area0)` and `uVol * (vol - vol0) / 3`, the kernel's own
expressions with the kernel's own guards. Same terms, reassociated: the area
and volume forces move from the in-kernel scaling by rounding only (measured
below); the bending force, every energy, every area and normal never saw the
factor and are unchanged. Local lipid flow gives each face its own tension in
place of the global factor, which is not a scalar to pull out, so the mode is
off there and the kernel scales as before.

### The touched set and the incremental pass

`Mesh::evaluate_edge_flip()` records `flip_patch_faces(iEdge)` for every
attempted flip, before its trial (the set is invariant under the flip), as a
union in `Mesh::flipTouchedFaces`. Attempted, not accepted: the trial flips
the edge and flips it back, and each `flip_edge()` rebuilds the one-rings of
the patch faces. The restored mesh is the same mesh, but a rebuilt one-ring can
list the same twelve control points from a different starting corner, and the
kernel's sum over them then rounds differently -- a few ulp on three regular
faces, in the first run of the gate below, which is exactly what a bit-exact
pass cannot leave stale. On `box_nu5` that is about 4.5 attempts per step and
some eighty faces of 960. `Mesh::Compute_Energy_And_Force_After_Flips()` then:

1. zeroes the vertex forces and every face-energy slot except the patch
   energy, which untouched faces keep from the last pass -- the state
   `clear_force_on_vertices_and_energy_on_faces()` leaves, minus what is
   still valid;
2. re-integrates the touched faces' element areas and volumes and re-sums
   `param.area` and `param.vol` over every face, in the order the whole-mesh
   pass sums them;
3. re-integrates the touched faces' patches into their blocks of the store;
4. rebuilds the vertex forces from the store, applying the constraint factors;
5. reruns the closed-form mesh-quality terms (the regularization or the
   tether, the triangle-shape term, the crease wall: sums over the edge table,
   cheap) and `finalize_energy_and_force()`, as the whole-mesh pass does.

It falls back to `Compute_Energy_And_Force()` when the mode is not active or
nothing was touched. `run_dynamics_flat()` calls it when the flip sweep, and
only the flip sweep, changed the mesh; a vertex-move sweep moves coordinates
and still gets the whole-mesh pass. Every whole-mesh pass clears the touched
set.

The store is rewritten, never incremented, so there is no running sum to drift:
after any number of incremental passes the vertex forces are the sum of blocks
each of which is what a whole-mesh pass would have written for that face.

### The parameter

    edgeFlipIncrementalForces = true     # default

`false` restores the whole-mesh pass after every sweep, with the constraint
factors scaled inside the kernel exactly as before: a fluid run with the key
off reproduces the previous tree byte for byte (section 2).

## 2. Verification

**Exactness.** `tests/test_incremental_flip_forces.cpp`, on a level-2 icosphere
with both global constraints pulling (`uSurf = 250`, `uVol = 30`, reference
area and volume 10% off the relaxed ones): after a sweep that landed several
flips, the incremental pass and a whole-mesh pass over the same mesh leave
identical forces (every slot and the total, every vertex), face energies,
element areas and volumes, mean curvatures, normals, totals and constraint
energies -- `EXPECT_EQ` in the serial build, `1e-12` relative in the OpenMP
build, where `sum_membrane_area_and_volume()` is a reduction whose order the
runtime sets. Forty sweeps with an incremental pass after each stay exact. The
gate also checks that the face energies that moved are all inside the touched
set, that the sweep's set contains every face incident to a corner of an
accepted flip, and that the fallback paths (nothing touched; parameter off)
give the whole-mesh pass.

**Reassociation only.** Same mesh, parameter on against off: bending forces
and all energies identical; area and volume forces differ by under `1e-11`
relative.

**The dynamics loop.** Twelve steps of the periodic 60 nm sheet with flips on,
the step loop of `run_dynamics_flat()` run in the test: the trajectory with the
incremental pass is the trajectory with the whole-mesh pass to under `1e-9` nm.

**Legacy byte-identity gate** (`docs/patch_kernel_performance.md` section 2
procedure; serial Release builds of `origin/main` at edb2759 and this branch):

| workload | files compared | result |
| --- | --- | --- |
| shipped periodic sheet, `sideX = sideY = 100`, 300 steps, harmonic bond off, `relaxArea = 5456` | the 8 dynamics outputs | byte-identical |
| `box_nu5` fluid box, 300 steps, `edgeFlipIncrementalForces = false` | all 17 CSV outputs: trajectory, forces, box, bead, `inputEdgeFlips.csv`, every `inputface_*.csv` | byte-identical |

So the store-then-scatter restructuring of the whole-mesh pass changes no
bits, with flips or without, and everything the mode changes is behind the
key.

## 3. Cost

`box_nu5`, 1000 steps from the flat start, `build-omp` Release `-O3`, one
thread, foreground (a job launched from a background shell on this Mac is
niced and reads a third of a core), parameter on against off:

| | `edgeFlipIncrementalForces = false` | `= true` | speedup |
| --- | --- | --- | --- |
| user time, 1000 steps | 68.3 s | 46.6 s | 1.46x |

User time, because a serial test suite was running alongside and stretched the
wall clock of both runs unevenly. The two runs accepted the same 1308 flips
(identical `inputEdgeFlips.csv`), report the same energies to the printed
precision, and their final control nets agree to 1.1e-12 nm after 1000 steps:
the reassociated constraint factor is the only difference, and it stays at the
rounding level.

The profile put the post-flip pass at 35-37% of the run, so removing it
outright would be 1.5-1.6x; the pass still re-integrates the trial faces
(about 80 of 960 per step at `nu = 5`), re-sums the closed-form terms and
redoes the finalisation, which is where the rest went. The saving scales with
the mesh: on a larger sheet the touched set stays about eighteen faces per
attempt while the whole-mesh pass grows with the face count.

What is left is the first force pass, which the whole-mesh profile puts at
over half the run and which is the patch kernel itself: fewer kernel calls per
face (the tail-corrected 6/5/5/5 depths, PR #21) and a vectorised
shape-function contraction are the next levers.
