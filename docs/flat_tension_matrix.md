# The flat periodic sheet at the zero default: three reservoirs, four fluidities

**Runs:** 2026-10-14, `~/slimed_runs/flat_tension/matrix`, cluster copy `slimed/flat_tension`
(job 31990185). `analysis/sticky_bead_runs.py --matrix --potential none --steps 600000`:
the sticky-bead matrix of `docs/bead_push_benchmark.md` section 9 without the bead, at the
code's default target of 0 pN/nm, 600 µs each with frames every 0.2 µs, no per-face output.
**Analysis:** `analysis/flat_tension_matrix.ipynb` (`SLIMED_FLAT_RUNS`); its figures,
`stats.md` and the animations (one per run and `matrix.gif`) are written beside the runs
and are not tracked. The first 120 µs of every run is burn-in.

The question was whether the three reservoirs -- the elastic box, the global relaxed area
and the per-face relaxed area -- put the same sheet on the same footing when each is told
"zero tension", and whether fluidity changes the sheet's spectrum in space or in time.
The answers: the box does what its note promises, the relaxed-area modes do not hold the
same tension at the same setting, and fluidity changes nothing the spectrum can see.

## 1. Each reservoir holds its own variable

Every reservoir sits on its target to within the noise: the box's frame tension at
-0.01 ± 0.01 pN/nm, the face mean of $K(S - S_0)/S_0$ at -0.004 ± 0.008 in all eight
relaxed-area runs, for every fluidity. The fixed box, with no reservoir, carries the
tension its thermal roughness costs at $A = A_0$.

The box's reading decomposes the same way in all four of its runs (pN/nm):

| term | frozen | ν = 0.05 | ν = 0.5 | ν = 5 |
| --- | ---: | ---: | ---: | ---: |
| bending virial (the undulations' pull) | -0.102 | -0.103 | -0.104 | -0.103 |
| mesh terms (tether, altitude floor, crease wall) | -0.164 | -0.142 | -0.144 | -0.143 |
| vertex gas, 480 vertices, $-Nk_BT/A_p$ | -0.194 | -0.194 | -0.194 | -0.194 |
| area constraint, $K(A - A_0)/A_0$ | +0.448 | +0.430 | +0.432 | +0.430 |
| **frame tension** (target 0) | **-0.013** | **-0.010** | **-0.009** | **-0.011** |

So at $\tau = 0$ the elastic membrane under the box sits at an internal stretching
tension of +0.43 pN/nm, balancing the gas, the undulations and the vertex fluid's own
pressure, and the box shrank 0.45 % in area to get there -- the picture of
`docs/elastic_box.md` section 6.3, now over 600 µs and at every fluidity.

## 2. What the spectrum reads

$\langle|h_q|^2\rangle$ on the limit surface, fitted as
$1/(AS) = (\kappa_c/k_BT)q^4 + (\sigma/k_BT)q^2$ over $|q|d \le \pi/2$, eight time blocks:

| run | κ_c (pN·nm) | σ spectrum (pN/nm) | slope | h r.m.s. (nm) |
| --- | ---: | ---: | ---: | ---: |
| box, frozen | 78.9 ± 3.2 | +0.39 ± 0.14 | -4.00 | 1.35 |
| box, ν = 0.05 | 87.6 ± 2.3 | +0.25 ± 0.08 | -4.04 | 1.38 |
| box, ν = 0.5 | 89.6 ± 2.6 | +0.15 ± 0.11 | -4.07 | 1.34 |
| box, ν = 5 | 75.1 ± 6.1 | +0.77 ± 0.21 | -3.82 | 1.10 |
| S0 global, frozen | 75.3 ± 2.0 | -0.11 ± 0.07 | -4.78 | 3.25 |
| S0 global, ν = 0.05 | 83.8 ± 1.0 | -0.18 ± 0.04 | -4.81 | 3.10 |
| S0 global, ν = 0.5 | 84.8 ± 0.9 | -0.23 ± 0.02 | -4.89 | 3.24 |
| S0 global, ν = 5 | 85.9 ± 0.9 | -0.24 ± 0.02 | -4.86 | 3.01 |
| S0 local, frozen | 77.1 ± 1.5 | -0.12 ± 0.05 | -4.74 | 3.00 |
| S0 local, ν = 0.05 | 82.3 ± 3.3 | -0.06 ± 0.10 | -4.64 | 2.60 |
| S0 local, ν = 0.5 | 85.6 ± 1.0 | -0.14 ± 0.06 | -4.62 | 2.34 |
| S0 local, ν = 5 | 84.5 ± 0.6 | -0.11 ± 0.03 | -4.55 | 2.27 |
| fixed box, ν = 0.5 | 78.3 ± 8.0 | +1.88 ± 0.30 | -3.25 | 0.74 |

Input $\kappa_c$ = 83.4 pN·nm.

- **Under the box the sheet is the tension-free Helfrich sheet.** Slope -4.0, $\kappa_c$
  within 5-8 % of the input (the frozen mesh reads low by the lattice's extraordinary-free
  stiffness, the fluid ones high by a few per cent, as in `docs/fluidity_results.md`), and a
  spectral tension of +0.15 to +0.4 pN/nm, small but not zero. The box's own reading is
  zero, so the difference is a term the frame feels and the height modes do not: the mesh
  term's -0.14, the tether's two-dimensional pressure on the frame, which the area
  constraint balances but which does not act on an undulation. The ν = 5 run's +0.77 ± 0.2
  rests on its two lowest modes, which 480 µs covers only three relaxation times of, and
  is not distinguished from the others at this length.
- **Under either relaxed-area mode the sheet is under compression.** The spectral tension
  is -0.06 to -0.24 pN/nm, the slope -4.6 to -4.9, the r.m.s. height 2.3-3.3 nm against
  1.35, and the lowest mode sits five to eight times above the $q^{-4}$ line. This is not a
  failure of the control -- the stretching tension is held at zero exactly -- it is what
  "zero" means to these modes. They hold $K(S - S_0)/S_0$ and know nothing of the vertex
  gas ($Nk_BT/A_p$ = 0.19 pN/nm on this tile) or of the undulations' entropic pull
  ($N_hk_BT/2A_p$ = 0.10); the box's target includes both. A relaxed-area sheet at
  $\sigma_{\rm res} = 0$ therefore has a frame tension of about -0.3 pN/nm, which on a
  100 nm tile is the buckling threshold $-\kappa_c(2\pi/L)^2$ = -0.33 of the lowest mode.
  The sheet does not buckle, but it is as soft as it can be without doing so. **To put the
  relaxed-area modes on the box's footing, set `lipidReservoirTension` to about +0.29
  pN/nm on this tile** (or add the gas term to the reservoir's exchange tension in the
  code, which is not done). The sticky-bead matrix ran all three at +0.5, where the
  offset is a 60 % difference in the frame tension the three sheets actually had; it did
  not change the wrap, the force or the neck there, but it is why the relaxed-area runs'
  spectral fits returned the lower $\sigma$.
- **The fixed box** at $A = A_0$ reads +1.9 pN/nm from the spectrum: the thermal
  roughness stretched against the modulus, with the slope bent to -3.25. The control net's
  own $K(A - A_0)/A_0$ is +4.7 because the net is rougher than the limit surface it
  controls; it is the spectrum's number that the undulations feel.
- **Fluidity does not change the spectrum.** Within each reservoir the frozen and the
  three fluid runs agree in $\sigma$ to the error bars and in $\kappa_c$ to the known few
  per cent, from 0 to 1200 accepted flips per µs.

## 3. How fast the sheet moves

The relaxation rate of every mode, from the integral of its autocorrelation to $1/e$
(`membrane_spectrum.relaxation_rate`), against the free-draining prediction
$\Gamma = (Da/k_BT)(\kappa_c q^4 + \sigma q^2)$ with $a$ = 21.65 nm² per vertex and the
spectrum's $\sigma$ (floored at zero). Over the 25-29 modes the record resolves in the
continuum window:

| | box | S0 global | S0 local | fixed |
| --- | --- | --- | --- | --- |
| Γ / Γ_pred, frozen → ν = 5 | 0.92, 1.00, 1.02, 0.88 | 0.89, 0.93, 0.96, 0.97 | 0.90, 0.95, 0.96, 0.97 | 0.93 |
| d ln Γ / d ln q | 4.0, 4.0, 4.1, 3.8 | 4.2, 4.3, 4.3, 4.3 | 4.2, 4.3, 4.2, 4.3 | 2.9 |

The vertices are overdamped Brownian particles with no hydrodynamic coupling, and that is
what the sheet does: $\Gamma \propto q^4$ at zero tension, at the rate the per-vertex
diffusivity sets, to within 10 % for every reservoir and every fluidity. Flips are
Monte Carlo moves and carry no clock, and the rates say so: the frozen and the fluid
sheets relax identically. The relaxed-area modes' slightly steeper slope is their
negative tension (a $-|\sigma|q^2$ term slows the long modes), and the fixed box's 2.9 is
its positive one. Neither reservoir's own dynamics -- the box strain at $mnK$ = 10 per µs,
the relaxed area at the exchange rate -- shows in the membrane's.

## 4. What to carry forward

- The box at $\tau = 0$ and the relaxed-area modes at $\sigma_{\rm res} = 0$ are not the
  same sheet. The difference is the vertex gas plus the undulation pull, +0.29 pN/nm on
  this tile, and it scales as $N/A_p$: a finer mesh at the same size is further off.
- The spectrum's tension under the box is positive by about the mesh term; a run that
  needs the spectrum at exactly zero should set the box target to minus that term
  (-0.15 pN/nm here) or read `tensionMeshTerm` and correct.
- A fluid sheet's control-net area must be summed over the energy-carrying faces only
  (the copy flags of `*mesh_faces.csv`); the image band's copies go stale under flips and
  make $K(A - A_0)/A_0$ read hundreds of pN/nm.
- 600 µs resolves the spectrum down to the second mode and the rates down to the third;
  the lowest mode of a 100 nm tile relaxes in 150 µs and needs a run several times longer
  for a tension measured to better than 0.1 pN/nm.
