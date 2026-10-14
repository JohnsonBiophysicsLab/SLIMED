# The elastic periodic box (`elasticBoxEnabled`)

**Status:** implemented, tested (`tests/test_elastic_box.cpp`); validation run
2026-09-28 (section 6)
**Base:** `feature/mixed-boundary-conditions @ 70de8ec`, on top of the bead
(`docs/bead_push_benchmark.md`) and the external force (`docs/external_force.md`)
**Scope:** an option for the Brownian dynamics (`membrane_dynamics`) under
which the periodic box relaxes so that the membrane's *frame tension* sits at
a prescribed value, instead of being whatever the fixed projected area happens
to impose. Off by default; a run with it off is the run it always was: the
shipped periodic workload is byte-identical (section 7). Worked example:
`data/example/elastic_box_zero_tension.params`.

## 1. Why

A periodic box of fixed size fixes the projected area `A_p = Lx Ly`, and the
membrane carries whatever tension that implies. The shipped 100 nm sheet under
`usMembraneStretching = 250` is the standard case: its undulations need about
1% of excess area, the area constraint pins the true area to the flat one, so
the sheet has to stretch to fluctuate at all and reads a fluctuation tension
of 3.7 pN/nm where the input asked for none
(`docs/fluctuation_spectrum.md`, section 5). A bead pushed through a periodic
sheet, or an insertion, has the same problem in a worse form: the deformation
needs area, and with no reservoir the only source is stretch, so the tension
climbs with the deformation instead of staying at the value the experiment
sets. (The bead benchmark worked around this with a clamped sheet and a
deliberately soft area constraint, `docs/bead_push_benchmark.md`.)

The elastic box is the reservoir. The strain of every periodic axis becomes a
Brownian coordinate of its own, driven by the difference between the target
tension and the tension the membrane actually exerts on its frame, and the
whole control net stretches with it. A membrane at `elasticBoxTension = 0` is
stress-free; at a positive value it is held under that tension, whatever it
is asked to do.

## 2. What the literature does, and what is adopted here

Three tensions have to be told apart (Schmid, *EPL* **95**, 28008, 2011;
Shiba, Noguchi & Fournier, *Soft Matter* **12**, 2373, 2016; Farago, *PRE*
**84**, 051914, 2011):

| symbol | name | what it is |
| --- | --- | --- |
| `sigma` | internal tension | the coefficient of the true area in the Hamiltonian: here `usMembraneStretching (A - A0) / A0`, the area constraint's Lagrange multiplier |
| `tau` | frame tension | `dF/dA_p`: the force per unit length the membrane exerts on its frame, the conjugate of the projected area |
| `r` | fluctuation tension | the coefficient of `q^2` in `<|h_q|^2> = kT / [A_p (r q^2 + kc q^4)]` |

They are not equal. `sigma` and `tau` differ by the entropy of the
undulations, at the Gaussian level `tau = sigma - N kT / (2 A_p)` for `N`
height modes at fixed mode count (Fournier & Barbetta, *PRL* **100**, 078103,
2008; Schmid 2011, eq. 4; Shiba et al. 2016, eqs. 9-10, who measure a residual
`sigma = 0.50 kT/a^2` at `tau = 0`). And the spectrum reads the *frame*
tension, `r = tau`, not the internal one: Shiba et al. find the two equal over
five decades of tension on a lattice Monge model with a fluctuating lattice
spacing, Farago 2011 shows it holds for any rotationally invariant Hamiltonian,
and FreeDTS's validation (below) finds the spectrum tension equal to the input
frame tension. So "zero tension" has to mean `tau = 0`, driven by the
derivative of the energy with respect to the projected area, and the test of
success is a `q^-4` spectrum, not a vanishing area-constraint multiplier.

How codes let the box find `tau`:

- **Dynamically triangulated surfaces, Monte Carlo.** FreeDTS (Pezeshkian &
  Ipsen, *Nat. Commun.* **15**, 548, 2024; Supplementary Note 3, eq. 7) adds
  `-tau A_p` to the energy and, with some probability per sweep, proposes
  `Lx -> Lx (1 + d)`, `Ly -> Ly (1 + d)` with every vertex's `(x, y)` rescaled
  by the same factor, accepted with
  `min[1, (A_p' / A_p)^N exp(tau dA_p - dE)]` -- an isotropic in-plane
  dilation, the full energy difference, and the `(A_p'/A_p)^N` Jacobian of
  rescaling `N` vertices' in-plane coordinates. A box move costs a full energy
  evaluation. Shiba et al. 2016 do the same on the lattice Monge model by
  changing the lattice spacing `a` every few sweeps; their heights are not
  rescaled, and no Jacobian appears.
- **Particle membranes, MD.** The constant-surface-tension ensembles of lipid
  simulation are Andersen-type extended systems on the box area (Zhang,
  Feller, Brooks & Pastor, *J. Chem. Phys.* **103**, 10252, 1995), integrated
  as a Langevin piston (Feller, Zhang, Pastor & Brooks, *J. Chem. Phys.*
  **103**, 4613, 1995; Kolb & Dünweg, *J. Chem. Phys.* **111**, 4453, 1999);
  Cooke, Kremer & Deserno (*PRE* **72**, 011506, 2005) run their coarse-grained
  bilayer at zero tension with a modified Andersen barostat on `x` and `y`
  only. The piston's force is the lateral stress from the virial,
  `P_xx = (N kT - sum_i x_i dU/dx_i) / V`, i.e. the affine-scaling derivative
  of the energy plus the ideal-gas term of the rescaled coordinates
  (Grønbech-Jensen & Farago, *J. Chem. Phys.* **141**, 194108, 2014, who write
  the Langevin piston for Brownian particles with exactly that force).
- **Overdamped box.** Bernetti & Bussi (*J. Chem. Phys.* **153**, 114107,
  2020) take the high-friction limit of the Langevin piston: a first-order
  stochastic equation for the box strain with a mobility, the instantaneous
  virial stress as its force, and a fluctuation-dissipation noise, including
  a semi-isotropic surface-tension variant. That is the form the Brownian
  dynamics here needs, since the vertices have no momenta either. The
  deterministic limit is Berendsen's weak coupling (*J. Chem. Phys.* **81**,
  3684, 1984), which relaxes to the target but does not sample an ensemble.

What SLIMED adopts is the overdamped box of Bernetti & Bussi with the force
of the Langevin piston: the affine-scaling derivative of the energy from the
nodal forces, the `-tau dA_p` term, and the ideal-gas term for the vertices
whose in-plane coordinates are actually integrated (section 3.2). The box
steps with the vertices from the same force evaluation, so nothing is
recomputed; with Monte Carlo flips on, the box is a coordinate the flip
sweep sees only through the geometry, which is what the Metropolis trial
already differences.

## 3. The physics as implemented

### 3.1 The box coordinate and its force

Let `lambda` be the strain of the periodic axes (one number for both under
the isotropic coupling, one per axis otherwise). Its effective potential is

```text
    U(lambda) = E(lambda x) - tau A_p(lambda) .
```

Loop subdivision is affine invariant, so an affine stretch of the control net
is the same stretch of the limit surface, and the derivative of the energy
with respect to it is a sum the force evaluation has already made:

```text
    dE/dlambda_a = sum over every vertex of  x_a  dE/dx_a  =  - W_a ,
    W_a = sum over every vertex of  x_a F_a            (the virial along a).
```

The sum runs over *every* vertex the energy-carrying patches reference --
sources, images and ghosts -- on the raw accumulated forces, before the
boundary handling folds the images onto their sources and zeroes the ghosts
(`Mesh::accumulate_virial()`, called from `Compute_Energy_And_Force()` just
before `manage_force_for_boundary_ghost_vertex()`). This is not a detail:
folding first would replace `x_img F_img` by `x_src F_img`, losing
`sum over images of offset * F_img`, which is the force carried across the
seam and most of the answer on a stretched sheet.
`ElasticBoxTest.TheVirialIsMinusTheDerivativeOfTheEnergy...` checks `-W`
against a central difference of the energy under `apply_box_scale()` to
1e-4 relative, on the per-vertex periodic sheet, on the global Periodic
mode's sheet, and on a strip that wraps along one axis only.

With `A_p` proportional to `lambda^n` (`n` periodic axes), the generalized
force is

```text
    -dU/dlambda = sum_a W_a + n tau A_p = n A_p (tau - tau_measured) ,
    tau_measured = -(1/n) sum_a W_a / A_p     [+ the gas term of 3.2] .
```

`tau_measured` is the tension the membrane exerts on its frame, per unit
length, as the lateral virial stress. The box then takes the same overdamped
step the vertices do,

```text
    d lambda = mu (-dU/dlambda) dt + sqrt(2 mu kT dt) xi ,
```

and the whole net is rescaled by `1 + d lambda` about the origin
(`DynamicMesh::apply_box_scale()`: every coordinate, previous and reference
coordinate, every row of `matMesh`, every image offset, the box lengths). The
mobility is set through an intensive parameter, `d lambda / dt = m (tau -
tau_measured)`, so `mu = m / (n A_p)`; it is fixed at setup from the area as
built, which keeps the noise free of the drift correction a
strain-dependent mobility would need.

**The frame area** is `Lx Ly` when both in-plane axes wrap. When one axis
wraps (a strip, or a tube along its axis) the tension needs a width to be a
force per length: `A_p = L W` with `W` the membrane area as built divided by
the period, which is the strip's width for a flat strip and the
circumference for a tube. Under `boundaryType = Periodic` the periods are the
grid's tile, `(nFace - 6) dFace`; under `Mixed` an axis wraps when some image
is offset along it and the period is that offset, which reads an imported
mesh as well as the generated sheet (`DynamicMesh::detect_periodic_box()`).

### 3.2 The vertex-gas term

Which coordinates does the box rescale? The in-plane ones, of every vertex.
On the solid sheet of `docs/fluctuation_spectrum.md` those are not degrees of
freedom -- in-plane motion is off and every vertex stays on its lattice site
-- so the rescale is a change of parameterization and nothing more. On a
fluid sheet (`inPlaneDynamicsEnabled = true`) they are Brownian coordinates
in absolute space, and a coordinate that is both integrated and rescaled
with the box carries the entropy of a particle in the plane: the free energy
of the vertex system is `F = -kT ln integral dx dz exp(-E/kT)`, and writing
`x = lambda s` puts `A_p^N` in front of the integral, so

```text
    tau = dF/dA_p = < -dE/dA_p at fixed s, z > - N kT / A_p .
```

That is the `(A_p'/A_p)^N` of FreeDTS's acceptance rule and the `N kT / V` of
every Langevin piston. It is not optional if the prescribed tension is to be
the one the spectrum reads: the vertex gas's entropy is `-N kT ln A` of the
*true* area too, so it acts on the height field as an internal tension
`-N kT / A_p`, and a box driven by the virial alone would settle where the
membrane's elastic tension cancels it -- a fluctuation tension of
`-N kT / A_p`, about -0.19 pN/nm for the 100 nm sheet's 252 vertices, which
is below the buckling threshold `-kc (2 pi / L)^2` of any sheet wider than
about 120 nm. With the term, the box at `tau = 0` holds the elastic membrane
at `+N kT / A_p` and the two cancel in the spectrum, which is what FreeDTS's
Supplementary Fig. 1 finds.

So `tau_measured` carries `-N kT / A_p` on every periodic axis, with `N` the
number of vertices the step integrates in the plane -- the same test
`next_step()` applies: not ghost, not a slaved duplicate or image, not
clamped -- when `inPlaneDynamicsEnabled` and `elasticBoxVertexGasTerm` are
both on, and zero otherwise. On the solid sheet nothing changes.
`ElasticBoxTest.TheVertexGasTermShiftsTheBoxByNkTOverTheArea` pins the
shift: a flat sheet with in-plane motion on settles where
`K (A - A0)/A0 = N kT / A_p`, and at `A = A0` with the term off.

### 3.3 The height entropy needs no term

The heights are not rescaled, so no Jacobian belongs to them -- and none is
needed, because their entropy reaches the box through the mean of the virial.
For a height mode the bending energy scales as `lambda^-2` under an in-plane
dilation (Monge gauge, small slopes), so its bending virial is
`-W_bend / (n A_p) = -E_bend / A_p` and, by equipartition, averages to
`-N_h kT / (2 A_p)`: the entropic pull of the undulations on the frame, the
`-N kT/(2 A_p)` of Schmid's and Fournier & Barbetta's relation, measured
rather than added. At `tau = 0` the area constraint has to balance it, so the
box settles where `sigma = usMembraneStretching (A - A0)/A0 = +N_h kT / (2 A_p)`
-- about +0.095 pN/nm for the 100 nm sheet -- and the spectrum reads
`r = tau = 0`. The log writes both halves (`tensionBending`,
`tensionAreaConstraint`), so this is checked directly in section 6.

### 3.4 The noise, and what ensemble this is

The vertices step in absolute coordinates and the box in the strain, at
fixed scaled coordinates. Together they are one Euler step of an overdamped
system whose mobility is diagonal and, with `mu` fixed, independent of the
coordinates it multiplies, so the Ito scheme samples `exp(-U/kT)` in the
scaled coordinates without a drift correction, and the vertex-gas term above
turns that into the canonical ensemble in absolute ones. With
`elasticBoxThermal = true` the box strain samples the constant-tension
ensemble: on a flat sheet under the area constraint its variance is
`kT / (n^2 A_p K)` (`n = 2`, `K = usMembraneStretching`), which
`ElasticBoxTest.TheThermalBoxSamplesTheStrainVarianceOfTheEffectivePotential`
recovers to within the explicit step's `2/(2 - rate dt)`; that is a strain
standard deviation of 9e-4 on the 100 nm sheet, 0.17% of the frame area.
With it false the box relaxes deterministically to where the *mean* virial
balances the target -- Berendsen's coupling -- which is smoother and is not
an ensemble; the stationary spectrum is the same in the limit of a slow box.

The time scale: with the area constraint as the area modulus the strain
relaxes at `m n K` per microsecond, 0.5 per us at the defaults (`m = 0.001`,
`K = 250`), i.e. 2 us or 1000 steps at `dt = 0.002` -- fast against the
slowest undulation of the 100 nm box (44 us) and slow against a step. The
startup report prints the rate and warns above `rate dt = 0.1`, where an
explicit step turns stiff; `elasticBoxMaxStrainPerStep` clips any single
increment beyond it and counts the clips, which are reported at the end. A
run that clips is one whose mobility is too large for its step.

### 3.5 Isotropic or anisotropic

`elasticBoxCoupling = isotropic` (the default) gives every periodic axis the
same strain, driven by the mean tension, so the box keeps its shape. The
alternative lets each axis follow its own tension component. On a fluid
sheet the aspect ratio is then a soft mode -- the energy does not care how
the projected area is shared between `x` and `y` -- and it wanders; the
reason to want it is a buckled or otherwise anisotropic state, whose
`tau_x` and `tau_y` differ (Noguchi, *PRE* **83**, 061919, 2011) and which
an isotropic strain cannot relax.

### 3.6 The bead and the external force

The bead is not part of the frame: it stays where it is while the membrane
stretches past it, so its wall's virial enters `tau_measured` at the
vertices' positions like every other term -- a bead pressed into the
membrane raises the tension, and the box feeds it area from the target. It
is logged separately (`tensionBead`). The external force is a constant force
per vertex, not a function of the box, and takes no part in the box's force.
How the box settles, lags and reads with a bead inside, and what the target
tension changes, is `docs/elastic_box_with_bead.md`.

## 4. Parameters

```text
elasticBoxEnabled = true
elasticBoxTension = 0.0             # tau, pN/nm: the frame tension the box is driven to
elasticBoxMobility = 0.001          # m, nm/(pN us): d(strain)/dt = m (tau - tau_measured)
elasticBoxThermal = true            # Brownian kick on the box; false = Berendsen-like relaxation
elasticBoxCoupling = isotropic      # isotropic | anisotropic
elasticBoxMaxStrainPerStep = 0.001  # cap on |d strain| per step; clipped steps are counted
elasticBoxVertexGasTerm = true      # -N kT ln(A_p) of the in-plane vertices (acts only with inPlaneDynamicsEnabled)
```

`Param` prints them with the rest, `import_kv_string()` reads them,
`ElasticBoxTest.TheParameterKeysParse` pins the spelling. Enabling the box on
a mesh that wraps along no axis, with a mobility or a cap that is not
positive, or with a coupling that is neither word, is refused at setup.

## 5. Where it lives

| piece | where |
| --- | --- |
| the virial, per term and axis | `Mesh::virialCurvature/Area/Regularization`, `Mesh::accumulate_virial()` (`src/energy_force/Compute_energy_and_force_on_mesh.cpp`) |
| periodic axes, periods, the rescale | `DynamicMesh::detect_periodic_box()`, `box_projected_area()`, `apply_box_scale()` (`src/dynamics/Elastic_box.cpp`) |
| the solver's copy of the image offsets | `SurfaceSolver::refresh_image_offsets()` |
| the box step and its log | `DynamicModel::setup_elastic_box()`, `refresh_box_tension()`, `update_elastic_box()`, `append_box_log()` |
| the driver | `run_dynamics_flat()`: the box step after `postprocess_boundary()` and before the vertex copy; the log at every trajectory frame; a summary at the end |

The order inside an iteration is: control net to surface, the Brownian step
of the vertices, surface to control net, the boundary post-process, **the box
step and rescale**, the copy into the vertices, the force evaluation. The box
and the vertices therefore move from the same forces, and the rescale is
applied to the post-processed net, so the global mode's ghost band -- which
the post-process rebuilds by copying displacements -- ends up at its partner
plus the stretched lattice translation
(`ElasticBoxTest.TheStretchKeepsTheGhostBandOfTheGlobalPeriodicMode`), and a
per-vertex image at its source plus its stretched offset, in the solver's
right-hand side too (`...KeepsEveryImageOnItsSourceAndTheSolverAgrees`).

**The log**, `<prefix>box.csv`, one row per trajectory frame (the first for
the configuration before any step), each carrying the tension the step just
taken was driven by:

```text
iteration, time_us, scaleX, scaleY, scaleZ, boxLengthX_nm, boxLengthY_nm, boxLengthZ_nm,
frameArea_nm2, membraneArea_nm2, tension_pN_nm, tensionX_pN_nm, tensionY_pN_nm, tensionZ_pN_nm,
tensionBending_pN_nm, tensionAreaConstraint_pN_nm, tensionMeshTerm_pN_nm, tensionBead_pN_nm,
tensionVertexGas_pN_nm, clippedSteps
```

`tension` is `tau_measured`, the mean over the periodic axes of the per-axis
columns; the five `tension*` columns after them are its split by term and
sum to it. `analysis/fluctuation_report.py --box-log <that file>` stretches
the lattice by the log's mean `scaleX, scaleY` before transforming, so `q`
and `A_p` are the box's actual ones.

## 6. Validation

### 6.1 Unit tests (`tests/test_elastic_box.cpp`, 12 tests)

- the virial equals minus the finite-difference derivative of the energy
  under the rescale, on three boundary layouts (1e-4 relative);
- the rescale keeps every image on its source and the solver's offsets in
  step, and keeps the global mode's ghost band on its partners;
- a flat sheet pre-stretched by 2% relaxes to `A_p = A0 (1 + tau / K)` to
  1e-6 for `tau = 0` and `tau = 2 pN/nm`, keeping its aspect ratio; the
  anisotropic strip moves only its periodic axis;
- the thermal box's strain variance is `kT/(n^2 A_p K)`;
- the vertex-gas term shifts the box by `N kT / A_p` when in-plane motion
  is on, and not otherwise;
- the box is refused where nothing wraps, and with the flag off the virials
  stay zero and the box never moves.

### 6.2 The run: the 100 nm sheet at zero frame tension

The configuration of `docs/fluctuation_spectrum.md` section 5 with the area
constraint on -- the "tension" run: `boundaryType = Periodic`, 100 nm,
`lFace = 5`, `kc = 83.4`, `usMembraneStretching = 250`, `dt = 0.002`, the
dense solver -- plus the box at `elasticBoxTension = 0`, thermal, isotropic,
`m = 0.001`; 300 000 steps (600 us), a frame every 100
(`~/slimed_runs/elastic_box_zero`, `data/example/elastic_box_zero_tension.params`,
serial Release binary, 14 min). Nothing was clipped.

**The box.** It shrank to its plateau inside the first 100 us (the box's own
relaxation time is 2 us; the sheet's excess area builds on the slowest
mode's 44 us) and then fluctuated about it with no drift: the mean strain
per 100 us window was 0.99659, 0.99692, 0.99674, 0.99607, 0.99620, 0.99675.
Over 150-600 us, with block standard errors (50 frames a block):

| | value |
| --- | ---: |
| frame area `A_p / A0 - 1` | -0.69% +- 0.02% |
| membrane area `A / A0 - 1` | +0.06% +- 0.01% |
| excess area `A / A_p - 1` | +0.76% |
| strain standard deviation | 1.08e-3 (0.88e-3 from the area constraint alone) |
| frame tension `tau_measured` (target 0) | +0.06 +- 0.03 pN/nm |
| bending virial | -0.0886 +- 0.0003 pN/nm |
| area constraint `sigma = K (A - A0)/A0` | +0.145 +- 0.031 pN/nm |

So the sheet keeps its true area at `A0` and takes the 0.76% of excess area
its undulations need out of the frame instead of out of stretch. The
bending virial is the entropic pull of section 3.3, measured:
`-<E_bend>/A_p = -491/5418 = -0.091 pN/nm` from the energy record, against
the virial's -0.089 (the 2% is the small-slope approximation), and
`-N_h kT/(2 A_p) = -0.097` at equipartition (the record's bending energy sits
at 93% of `252 kT/2`, as the spectrum below also says). The area constraint
balances it: a stress-free frame carries a positive internal tension, here
+0.15 with the box still 0.06 above its target on average -- Shiba, Noguchi
& Fournier's `sigma > 0 at tau = 0`, in a Loop-surface membrane.

**The spectrum.** `fluctuation_report.py --box-log inputbox.csv --burn 750`
on the relaxed lattice, against the same analysis of the fixed-box control
(`~/slimed_runs/tension`) and of a genuinely tensionless fixed-box sheet
(`~/slimed_runs/pure_a`, no area elasticity), all three cut to the same
statistics -- 2251 frames, 450 us, the lowest bin at `T / tau_q = 10`:

| run | box | fitted `sigma` (pN/nm) | log-log slope | fitted `kc` (pN.nm) |
| --- | --- | ---: | ---: | ---: |
| `tension`: area constraint, fixed box | fixed | 4.36 | -2.61 | 36.2 |
| **`elastic_box_zero`: area constraint, box at `tau = 0`** | **elastic** | **0.89** | **-3.45** | **77.3** |
| `pure_a`: no area constraint | fixed | 0.43 | -3.70 | 81.9 |

Ratios `<|h_q|^2> / [kT / (A kc q^4)]` bin by bin, lowest seven `q`:

```text
 q (1/nm)      0.085  0.121  0.171  0.191  0.242  0.270  0.308
 tension       0.125  0.304  0.474  0.562  0.692  0.731  0.794
 elastic box   0.395  0.761  0.751  0.836  0.852  0.888  0.914
 pure_a        0.556  0.880  0.886  0.912  0.906  0.886  0.914
```

**Correction (2026-10-05).** These three runs use the global `Periodic`
mode, whose inexact seam inflates the spectrum tension once the area term
is on: the control's internal tension is only 1.37 pN/nm (`K (A - A0)/A0`
from its energy record) against the 4.36 its spectrum reads, while the same
sheet under `boundaryType = Mixed` reads 1.42 at an internal tension of
1.68. Under `Mixed` the box at `tau = 0` gives slope −4.15 and a fitted
tension of −0.11 pN/nm over 200 000 steps
(`~/slimed_runs/lipid_flow_prototype/A9_static_global_elastic_box_tau0`),
so the residual below is the legacy seam's, not the box's. See
`docs/lipid_flow_plan.md` section 8.

The box takes the sheet from the tension-dominated `q^-2.6` of the control
to within reach of the tensionless reference: the same curve from
`q = 0.24 /nm` up, and lower by 0.1-0.15 in the four longest modes, which
are the ones sampled ten to forty times over in 450 us and carry 10-15%
errors. What is left, 0.9 against the reference's own 0.4 of finite-length
bias, is not resolved by this run; the truncation table of
`docs/fluctuation_spectrum.md` section 5 shows the same reference moving
from slope -3.57 at a quarter of its length to -4.02 at full length, and
telling a residual `r` of 0.1 from 0.5 needs that full length, three to four
times this run. The tension the fixed box imposed is gone.

### 6.3 A fluid sheet with the box on

The same box on the fluid configuration of `docs/fluidity_results.md`
(`fluid_e`: flips on, flat tether with the lower wall at 0.7, altitude floor,
crease wall, in-plane motion on, `usMembraneStretching = 250`), 60 000 steps
at `dt = 0.001` (`~/slimed_runs/elastic_box_fluid`, `build-omp` binary, 8
threads). The run is healthy -- 3566 flips accepted of 12 410 admissible
(29%), no clipped box step, nothing diverged -- and the box does what section
3.2 says it should. Averaged over 20-60 us, block standard errors:

| term | pN/nm |
| --- | ---: |
| bending virial | -0.093 +- 0.001 |
| mesh-quality term (tether, floor, wall) | -0.140 +- 0.003 |
| vertex gas, 221 integrated vertices | -0.170 |
| area constraint, `K (A - A0)/A0` | +0.330 +- 0.123 |
| **frame tension** (target 0) | **-0.07 +- 0.12** |

The area constraint carries the tension that balances the vertex gas, the
mesh terms and the undulations, so the elastic membrane sits at an internal
tension of about +0.3 pN/nm while its frame is stress-free; the box shrank
by 0.47% in area to get there and its strain fluctuates with a standard
deviation of 1.0e-3 (0.9e-3 predicted from the area constraint alone). The
mesh term's -0.14 is the vertex fluid's own pressure at this packing -- the
lower tether wall pushing outward -- which is why a fluid run with the box
should watch `tensionMeshTerm`: a strongly negative value means the tether
range and the vertex density are fighting, and the area constraint is
paying for it. Whether the fluid spectrum then reads `r = 0`, as FreeDTS's
does with the same Jacobian, needs a run ten times longer than this smoke
test and is not claimed here.

## 7. What did not change

- The global boundary modes run the code they always ran when the box is
  off: the shipped periodic workload -- `membrane_dynamics` on the 100 nm
  `Periodic` sheet for 300 steps and `continuum_membrane` for 25 iterations
  -- produces byte-identical output files (all eight of the dynamics', all
  thirteen of the minimizer's) with the binaries built before and after this
  change. The virial is only accumulated with the box on.
- `EnergyForce.csv` and the trajectory files keep their format; the box has a
  file of its own.

## 8. Limits

- **Dynamics only.** The minimizer (`continuum_membrane`) ignores the box; a
  box coordinate in the conjugate-gradient line search is a separate change.
- **The stretch is about the origin.** The generated sheet is centred there;
  an imported mesh that is not will see its centre drift by
  `(lambda - 1)` times its offset, which for a per-cent strain is harmless
  and for a bead at a fixed position is what "the bead stays put" means.
- **The reference configuration stretches with the frame.** The solid
  regularization's `coordRef` is rescaled with the box, so a uniform stretch
  costs no regularization energy; this is the Lagrangian view, and is what
  keeps the mesh-quality terms blind to the box.
- **One box for the whole membrane.** The flip sweep and the box are
  independent moves; the box does not know about connectivity and needs
  not to.
