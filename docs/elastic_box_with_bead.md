# The elastic box with a bead inside: how the box moves with the target tension

**Companion to:** `docs/elastic_box.md` (the box itself), `docs/bead_push_benchmark.md`
section 9 (the sticky and pushing bead matrices, run at a target of 0.5 pN/nm).
**Defaults:** `elasticBoxTension = 0` and `lipidReservoirTension = 0`; the matrix generator
`analysis/sticky_bead_runs.py` writes 0 unless `--tension` says otherwise.

This note is about one question: a rigid bead sits inside a periodic tile whose box is
driven to a target frame tension $\tau$. The bead pulls membrane out of the plane, or
pushes it. What does the box do, how fast, and what changes when $\tau$ is set to zero
rather than to the 0.5 pN/nm the benchmarks used? Everything below is in terms of the
quantities the code logs, with the numbers of the 2026-10-12 sticky matrix where they
exist.

## 1. The two balances

The box has one coordinate, the isotropic strain $\lambda$ of the periodic axes. The
projected area of the tile is

$$A_p(\lambda) = A_{p,0}\,\lambda^{2},$$

and the strain moves as an overdamped Brownian coordinate in the effective potential
$U(\lambda) = E(\lambda x) - \tau A_p(\lambda)$:

$$\frac{d\lambda}{dt} = m\,\bigl(\tau - \tau_{\rm box}\bigr) + \sqrt{\frac{2\,m\,k_BT}{n A_p}}\;\xi(t),
\qquad n = 2 .$$

$m$ is `elasticBoxMobility`, in nm per pN per µs, and $\tau_{\rm box}$ is the tension the
membrane exerts on its frame, read off the forces the step has already computed:

$$\tau_{\rm box} = -\frac{1}{n A_p}\sum_{a=x,y}\;\sum_{i} x_{i,a} F_{i,a} \;-\; \frac{N k_BT}{A_p}.$$

The double sum is the in-plane virial of every force on every vertex: bending, the area
constraint, the mesh regularizers, and the bead's wall or bond. The last term is the
vertex gas, present because the in-plane coordinates of a fluid sheet are rescaled with
the box (`docs/elastic_box.md` section 3.2). The box, then, is a controller with one
set-point: it changes the projected area until the *total* in-plane virial of the tile,
bead included, reads $\tau$.

The membrane has its own balance. Its area is the projected area plus whatever is out of
the plane,

$$A = A_p + A_{\rm ex},$$

where $A_{\rm ex}$ is the excess area: the thermal undulations (a few per cent) and, with
a bead, the wrap or the tent. Under the global area constraint the stretching tension is
$\sigma = K\,(A - A_0)/A_0$ with $K$ = `usMembraneStretching` and $A_0$ the relaxed area.

## 2. Where the box settles

Take the bead static: it holds $A_{\rm ex}$ out of the plane and does not let go. The box
stops moving where $\tau_{\rm box} = \tau$. To leading order the frame tension is the
stretching tension (the undulations' entropic pull and the gas term cancel on a sheet of
this size, section 3.3 of the box note; the bead's wall is part of the same virial and
not a separate offset -- see section 5), so

$$K\,\frac{A_p^{*} + A_{\rm ex} - A_0}{A_0} = \tau
\quad\Longrightarrow\quad
A_p^{*} = A_0\Bigl(1 + \frac{\tau}{K}\Bigr) - A_{\rm ex},
\qquad
\lambda^{*} = \sqrt{\frac{A_p^{*}}{A_{p,0}}} .$$

Two things follow. First, **the box's size is set by the bead, not by the target**: the
term $\tau/K$ is 0.5/250 = 0.2 % at the benchmark's target and exactly zero at the
default, while a wrapped 10 nm bead with its neck takes $A_{\rm ex} \approx 1500$ nm$^2$
out of a 10 400 nm$^2$ tile, 14 %. The sticky matrix measured this directly: the frame
shrank by 1485 nm$^2$ ($\lambda^{*} = 0.926$) against a measured excess area of 1516 nm$^2$,
the difference being the thermal roughness the flat sheet already had. **At $\tau = 0$ the
box shrinks by exactly the area the bead holds**, and the sheet around the bead is slack:
the only elastic cost of drawing membrane onto the bead is bending.

Second, what the target *does* set is the shape of the membrane on the bead, because the
tension enters the force balance at the neck and the tent:

- The neck of a wrapped bead has the tether radius $r_{\rm neck} = \sqrt{\kappa_c / 2\sigma}$
  -- 9.1 nm at $\sigma = 0.5$ pN/nm, $\kappa_c = 83.4$ pN nm, which the fluid sticky runs
  reached (7-8.5 nm). At $\sigma \to 0$ this radius has no scale: a catenoidal neck has
  zero mean curvature and costs no bending, so a zero-tension sheet can close any neck it
  likes and the bead is engulfed as far as the mesh lets it. Expect a complete wrap with a
  narrow neck, and a box that keeps shrinking as long as the bonds keep capturing.
- The residual force on a held bead is tension times the neck's circumference plus the
  bending part, $F \simeq 2\pi\sqrt{2\kappa_c\sigma}$ for a tether: 57 pN at 0.5 pN/nm,
  of the same order as the -60 to -100 pN the fluid sticky runs kept pulling with. At
  $\tau = 0$ only the bending part remains, so the held bead feels a much smaller force.
- A bead that only pushes is held by the contact ring, and that ring slips below the
  equator when the tent is tight: at 0.5 pN/nm under the box the neck closed to 9.4 nm,
  below $R = 10$, and the sheet slid over the bead (`bead_push_benchmark.md` section 9.3).
  At $\tau = 0$ the tent is not pulled tight by anything, the contact ring stays on the
  upper cap, and the slip should not happen. This is the one qualitative difference a
  zero target makes to the pushing bead.

The box's own thermal fluctuation is negligible on this scale: in the constant-tension
ensemble $\langle \delta\lambda^2 \rangle = k_BT / (n^2 A_p K)$, a strain of $9\times10^{-4}$,
0.2 % of the area, two orders below the bead's 14 %.

## 3. How fast the box follows

Linearize the controller around its rest point. With the stretching tension as the stiff
term, $d\tau_{\rm box}/d\lambda \simeq n K$, so a disturbance of the frame tension relaxes at
the rate

$$\gamma = m\,n\,K ,$$

10 per µs at the matrix's $m = 0.02$ (0.5 per µs at the code default $m = 0.001$). That is
0.1 µs: against the undulations (the slowest mode of a 100 nm box relaxes in 44 µs) and
against a bead moving at 0.5 nm/µs the box is quasi-static. Its response to a *change of
target* is therefore immediate: switching $\tau$ from 0.5 to 0 would move the frame by
$\tau/K$ = 0.2 % of its area within a microsecond and change nothing else.

What is not immediate is the area the bead demands. A bead rising through the sheet pulls
membrane out of the plane at a rate $\dot A_{\rm ex}$, and to keep the tension at its target
the box must shrink at the same rate. From $\dot A_p = n A_p \dot\lambda$ and the strain
equation, the box can only shrink while it reads a tension above the target, and the
mismatch it needs is

$$\tau_{\rm box} - \tau = \frac{\dot A_{\rm ex}}{n\,m\,A_p} .$$

This is the lag. In the sticky matrix the excess area grew at 25-35 nm$^2$/µs through the
second half of the push, which with $n = 2$, $m = 0.02$ and $A_p \approx 10^4$ nm$^2$ puts the
box 0.06-0.08 pN/nm above its target; the logged frame tension sat at 0.60-0.65 pN/nm
against 0.5 for the 40 µs of the push, and was back at 0.49-0.52 within 20 µs of the
stop. The lag does not depend on the target: a box at $\tau = 0$ pulled on at the same
rate reads +0.06 to +0.1 pN/nm while the bead moves, and 0 once it stops. The lipid-flow
reservoirs do not have it, because their rate is set by the exchange constant and
not by a strain mobility: at the same moment they read 0.50-0.55, and the pushing bead
tented and held under them where it slipped under the box.

If the lag matters, the mobility is the knob. To hold the mismatch below $\Delta\tau$ while
the bead draws area at $\dot A_{\rm ex}$,

$$m \;\ge\; \frac{\dot A_{\rm ex}}{n\,A_p\,\Delta\tau},$$

0.035 nm/(pN µs) for 35 nm$^2$/µs and 0.05 pN/nm on this tile. The ceiling is the explicit
step: `elasticBoxMaxStrainPerStep` clips any increment beyond it, and the startup report
warns when $\gamma\,dt$ passes 0.1, which at $dt$ = 0.001 µs is $m$ = 0.2.

## 4. Pass-through and release

The box only follows the membrane; it has no memory. When the pushing bead slipped
through the sheet at 68 µs, the excess area collapsed from 750 nm$^2$ to 90 within 20 µs
and the frame re-expanded by the same amount, ending 50 nm$^2$ short of its original size
(the thermal roughness of a sheet at 0.46 pN/nm). The same will happen to a sticky bead
whose bonds are allowed to break (`beadBreakRange`) and to a bead that is withdrawn: the
area returns to the plane and the box gives it room at the target tension. Under a zero
target the returned area is absorbed at no tension at all, which is the difference
between a reservoir and a spring: nothing in the box resists the bead's work except the
bending it causes.

## 5. Reading the logged tension with a bead inside

The box controls the total virial, and the bead's wall is part of it. The log splits it:
`tensionBead` is the bead's share, $-\sum_i x_i\cdot F^{\rm bead}_i / (n A_p)$. A sticky
bead holds its bonded vertices just outside the sphere, so its bonds pull inward, the
in-plane virial of those forces is negative, and its share of the tension is positive:
0.12-0.20 pN/nm during the hold of the sticky matrix, 0.04-0.13 during the push; a bead
that only pushes contributes -0.04 while it tents the sheet and nothing once it is
through. None of this is an error in the reading. In mechanical equilibrium the total
virial of the tile is the force transmitted across its periodic seam per unit length,
which is the frame tension by definition, and on a mobile mesh that identity was checked
directly: the box's reading and the stress carried across the seam agree to 0.002 pN/nm.
The bead's share is the part of that transmitted stress the bead's grip accounts for, and
the membrane's elastic terms make up the rest. What it means for the reader is that at a
target of $\tau$ the far field of the sheet is at $\tau$, and the membrane on the bead is
held at whatever the bonds impose locally.

Two cautions carry over from the earlier assessment of the box with a bead. The stretch is
about the origin, so a bead off the centre of the tile sees the membrane slide past it
with every strain step -- a sticky bead 15 nm off-centre felt 17-20 pN toward the origin;
keep the bead at the origin, as the matrices do, or expect that drift. And a sheet that is
periodic along both axes is not held in $z$: with a bead pushing on it the whole tile
would drift, which is why the benchmarks pin the tile's mean height
(`meanHeightPinned`); the pin removes one Brownian mode and leaves the tension untouched.

## 6. What to expect from the matrices at a zero target

Running the sticky matrix at the new default should show, relative to 0.5 pN/nm: the same
frame shrinkage to within a few per cent (the bead sets it), a narrower neck and a more
complete wrap in the fluid runs, a residual force on the held bead that is a fraction of
the -60 to -100 pN measured, the same 0.06-0.1 pN/nm lag during the push, and a frame
tension that reads 0.00 with the same $\pm 0.03$ scatter over an 18 µs window. The
pushing matrix should lose its pass-through under the box and tent like the
lipid-flow runs do. Neither has been run; these are the predictions of the balances above.
