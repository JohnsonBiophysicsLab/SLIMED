# Membrane fluctuations from equipartition: the static spectrum and the free-draining time correlation

This note derives, from first principles, the two relations the fluctuation
analysis in this repository tests:

* the **static spectrum** $\langle h_{\mathbf q} h_{\mathbf q}^{*}\rangle$, from
  the Helfrich energy and equipartition;
* the **temporal correlation** $\langle h_{\mathbf q}(t)\,h_{\mathbf q}^{*}(t+\tau)\rangle$,
  from an overdamped Langevin equation with a *local* mobility and no
  hydrodynamic interaction, which is the dynamics the code actually integrates.

Throughout, every Fourier amplitude is written in **polar form**,

$$
h_{\mathbf q} = r_{\mathbf q}\, e^{i\theta_{\mathbf q}},
\qquad r_{\mathbf q} = |h_{\mathbf q}| \ge 0,\quad \theta_{\mathbf q}\in[0,2\pi),
$$

so that the amplitude and the phase of each mode are the variables, not the real
and imaginary parts. The statics then become a statement about the distribution
of $r_{\mathbf q}$ (Rayleigh) and of $\theta_{\mathbf q}$ (uniform); the dynamics
become a radial equation and a phase-diffusion equation.

The conventions are the ones used in [fluctuation_spectrum.md](fluctuation_spectrum.md)
and [resolution_limits.md](resolution_limits.md), so the formulas can be compared
directly with the measurements there. Symbols:

| symbol | meaning | SLIMED name |
| --- | --- | --- |
| $h(\mathbf r)$ | height of the membrane above the flat reference plane (Monge gauge) | $z$ of a vertex |
| $A = L_x L_y$ | projected area of the periodic box | $(n_{FaceX}-6)\,d_{FaceX}\times(n_{FaceY}-6)\,d_{FaceY}$ |
| $\kappa$ | bending rigidity | `kcMembraneBending` |
| $\sigma$ | lateral tension | no direct key: zero in a fixed box; `elasticBoxTension` / `lipidReservoirTension` set a target (`sigma` in input.params is an unrelated spontaneous-curvature length) |
| $k_BT$ | thermal energy | `KBT` |
| $D$ | per-vertex diffusion constant | `diffConst` |
| $a = A/N$ | projected area per vertex | $\tfrac{\sqrt3}{2}\,l_{Face}^{2}$ on the triangular lattice |
| $\Lambda$ | mobility per unit area of the continuum sheet | $D a / k_BT$ (derived in §4) |
| $K_q$ | mode stiffness $\kappa q^4 + \sigma q^2$ | — |
| $\beta_q$ | inverse variance $A K_q / k_BT$ of a mode | — |
| $\Gamma_q$ | relaxation rate $\Lambda K_q$ of a mode | — |

---

## 1. Setup: Helfrich energy of a nearly flat sheet

For a membrane that is a graph $z = h(x,y)$ over a flat reference plane, with
small slopes $|\nabla h| \ll 1$, the Helfrich energy with a tension term expands
to quadratic order in $h$ as

$$
E[h] \;=\; \frac12 \int_A d^2r \,\Big[\, \kappa\,(\nabla^{2} h)^{2} + \sigma\,(\nabla h)^{2} \,\Big].
$$

The first term is the bending energy ($2H \approx \nabla^2 h$ for small slopes,
$E_{bend} = \tfrac{\kappa}{2}\int (2H)^2 dA$). The second is the tension times the
excess true area over the projected area, $dA - d^2r \approx \tfrac12 (\nabla h)^2 d^2r$.
Spontaneous curvature and Gaussian curvature do not appear: the first is zero for
a symmetric bilayer and the second is a topological constant on a periodic sheet.

### Fourier expansion in polar form

On a periodic box of area $A$ take

$$
h(\mathbf r) = \sum_{\mathbf q} r_{\mathbf q}\, e^{i(\mathbf q\cdot\mathbf r + \theta_{\mathbf q})},
\qquad
r_{\mathbf q} e^{i\theta_{\mathbf q}} = \frac1A \int_A d^2r\; h(\mathbf r)\, e^{-i\mathbf q\cdot\mathbf r},
\qquad
\mathbf q = 2\pi\Big(\frac{n_x}{L_x}, \frac{n_y}{L_y}\Big).
$$

Each mode is a plane wave of amplitude $r_{\mathbf q}$ whose crests sit at
$\mathbf q\cdot\mathbf r = -\theta_{\mathbf q} \pmod{2\pi}$: the phase is the
position of the wave in the box, the amplitude is its height. Because $h$ is
real, $h_{-\mathbf q} = h_{\mathbf q}^{*}$, i.e.

$$
r_{-\mathbf q} = r_{\mathbf q}, \qquad \theta_{-\mathbf q} = -\theta_{\mathbf q}.
$$

The modes at $\mathbf q$ and $-\mathbf q$ are the same standing wave
$2 r_{\mathbf q}\cos(\mathbf q\cdot\mathbf r + \theta_{\mathbf q})$ seen twice, so
the independent variables are one $(r_{\mathbf q},\theta_{\mathbf q})$ pair per
$\pm\mathbf q$ pair. Let $\mathcal H$ denote any half of reciprocal space holding
one member of each pair; $\mathbf q=0$ is a rigid translation and is dropped.

### Diagonalising the energy

Insert the expansion into the energy. Using $\nabla^2 e^{i\mathbf q\cdot\mathbf r} = -q^2 e^{i\mathbf q\cdot\mathbf r}$,
$\nabla e^{i\mathbf q\cdot\mathbf r} = i\mathbf q\, e^{i\mathbf q\cdot\mathbf r}$, and the
orthogonality $\int_A d^2r\, e^{i(\mathbf q+\mathbf q')\cdot\mathbf r} = A\,\delta_{\mathbf q,-\mathbf q'}$:

$$
E = \frac12 \sum_{\mathbf q,\mathbf q'} r_{\mathbf q} r_{\mathbf q'} e^{i(\theta_{\mathbf q}+\theta_{\mathbf q'})}
\big[\kappa\, q^2 q'^2 - \sigma\,\mathbf q\cdot\mathbf q'\big]\, A\,\delta_{\mathbf q,-\mathbf q'}
= \frac{A}{2}\sum_{\mathbf q} \big(\kappa q^4 + \sigma q^2\big)\, r_{\mathbf q} r_{-\mathbf q}\, e^{i(\theta_{\mathbf q}+\theta_{-\mathbf q})}.
$$

With $r_{-\mathbf q} = r_{\mathbf q}$ and $\theta_{-\mathbf q} = -\theta_{\mathbf q}$ the
phase factor is exactly $1$ and

$$
\boxed{\,E = \frac{A}{2} \sum_{\mathbf q} K_q\, r_{\mathbf q}^{2}
          = A \sum_{\mathbf q\in\mathcal H} K_q\, r_{\mathbf q}^{2},
\qquad K_q \equiv \kappa q^{4} + \sigma q^{2}\,}
$$

Two things are now visible that the Cartesian form hides. The energy depends on
the **amplitudes only**; the **phases are cyclic coordinates** (translating a
wave inside a periodic box costs nothing). And the energy of each pair is
$A K_q r_{\mathbf q}^2$ with no $\tfrac12$: the harmonic "spring constant" of the
amplitude is $2AK_q$.

---

## 2. Static spectrum from equipartition

### The Boltzmann measure in polar coordinates

Equipartition is a statement about the canonical measure, so the measure must be
written in the variables being used. The map from the Cartesian pair
$(\mathrm{Re}\,h_{\mathbf q}, \mathrm{Im}\,h_{\mathbf q})$ to $(r_{\mathbf q},\theta_{\mathbf q})$
has Jacobian $r_{\mathbf q}$, so for each $\mathbf q\in\mathcal H$

$$
P(r_{\mathbf q},\theta_{\mathbf q})\, dr_{\mathbf q}\, d\theta_{\mathbf q}
\;\propto\; e^{-A K_q r_{\mathbf q}^{2}/k_BT}\; r_{\mathbf q}\, dr_{\mathbf q}\, d\theta_{\mathbf q}
\;=\; e^{-\beta_q r_{\mathbf q}^{2}}\; r_{\mathbf q}\, dr_{\mathbf q}\, d\theta_{\mathbf q},
\qquad
\beta_q \equiv \frac{A K_q}{k_BT}.
$$

The measure factorises mode by mode, and within a mode it factorises into a
$\theta$ part and an $r$ part:

$$
\int_0^{2\pi} d\theta \int_0^\infty r\, e^{-\beta_q r^{2}}\, dr = 2\pi\cdot\frac{1}{2\beta_q} = \frac{\pi}{\beta_q}.
$$

**Phase.** Nothing depends on $\theta_{\mathbf q}$, so it is uniform on $[0,2\pi)$
and independent of $r_{\mathbf q}$ and of every other mode:

$$
P(\theta_{\mathbf q}) = \frac{1}{2\pi},
\qquad
\langle e^{i n\theta_{\mathbf q}}\rangle = \delta_{n,0}.
$$

**Amplitude.** Normalising the radial factor,

$$
P(r_{\mathbf q}) = 2\beta_q\, r_{\mathbf q}\, e^{-\beta_q r_{\mathbf q}^{2}},
\qquad r_{\mathbf q}\ge 0,
$$

a **Rayleigh distribution** with scale $1/\sqrt{2\beta_q}$. Equivalently
$r_{\mathbf q}^2$ is exponentially distributed with mean $1/\beta_q$.

### Equipartition

The moment that defines the spectrum is

$$
\langle r_{\mathbf q}^{2}\rangle
= \frac{\int_0^\infty r^{3} e^{-\beta_q r^{2}} dr}{\int_0^\infty r\, e^{-\beta_q r^{2}} dr}
= \frac{1/(2\beta_q^{2})}{1/(2\beta_q)}
= \frac{1}{\beta_q}.
$$

Since $h_{\mathbf q}h_{\mathbf q}^{*} = r_{\mathbf q} e^{i\theta_{\mathbf q}}\, r_{\mathbf q} e^{-i\theta_{\mathbf q}} = r_{\mathbf q}^2$
with the phase cancelling identically,

$$
\boxed{\;
\langle h_{\mathbf q} h_{\mathbf q}^{*}\rangle
= \langle r_{\mathbf q}^{2}\rangle
= \frac{k_BT}{A\,(\kappa q^{4} + \sigma q^{2})}
\;}
\tag{1}
$$

The phase uniformity supplies the rest of the second-moment structure for free:

$$
\langle h_{\mathbf q}\rangle = \langle r_{\mathbf q}\rangle\langle e^{i\theta_{\mathbf q}}\rangle = 0,
\qquad
\langle h_{\mathbf q} h_{\mathbf q}\rangle = \langle r_{\mathbf q}^{2}\rangle\langle e^{2i\theta_{\mathbf q}}\rangle = 0,
\qquad
\langle h_{\mathbf q} h_{\mathbf q'}^{*}\rangle = 0 \;\;(\mathbf q'\neq\mathbf q).
$$

The middle identity says $\langle h_{\mathbf q} h_{-\mathbf q}^{*}\rangle = 0$: the
$\pm\mathbf q$ partners are uncorrelated *as complex numbers* even though their
amplitudes are identical. A mode is "circular": no phase is preferred.

### Where the two halves of $k_BT$ come from

The mean energy of a $\pm\mathbf q$ pair is
$\langle A K_q r_{\mathbf q}^{2}\rangle = A K_q/\beta_q = k_BT$. In Cartesian
variables this is the familiar $k_BT/2$ from each of two quadratic coordinates.
In polar variables only $r_{\mathbf q}$ carries energy and $\theta_{\mathbf q}$
carries none, so a naive "$k_BT/2$ per coordinate in the energy" would give
$k_BT/2$, which is wrong. The missing half comes from the Jacobian: the
radial measure $r\,dr$ is equivalent to an effective potential

$$
U_{\rm eff}(r) = A K_q r^{2} - k_BT \ln r,
$$

whose entropic term $-k_BT\ln r$ pushes the amplitude outward. The generalised
equipartition theorem, $\langle r\,\partial_r U_{\rm eff}\rangle = k_BT$, gives
$2AK_q\langle r^2\rangle - k_BT = k_BT$, i.e. $\langle r^2\rangle = k_BT/(AK_q)$
again. The same bookkeeping is why the most probable amplitude,
$r^{*} = 1/\sqrt{2\beta_q}$, is not zero even though the energy is minimal there.

Summed over the $N-1$ non-zero modes this is the check "$E_{bend} = \tfrac12 k_BT$
per free vertex" used in [resolution_limits.md](resolution_limits.md):
$k_BT$ per pair is $k_BT/2$ per $\mathbf q$ in the full sum.

### Higher moments: what the Rayleigh law predicts for an estimator

Everything about the amplitude distribution follows from $\beta_q$ alone:

$$
\langle r_{\mathbf q}\rangle = \frac12\sqrt{\frac{\pi}{\beta_q}},
\qquad
\langle r_{\mathbf q}^{4}\rangle = \frac{2}{\beta_q^{2}},
\qquad
\frac{\langle r_{\mathbf q}^{4}\rangle}{\langle r_{\mathbf q}^{2}\rangle^{2}} = 2,
\qquad
\langle \ln r_{\mathbf q}^{2}\rangle = -\ln\beta_q - \gamma_E .
$$

* The ratio $\langle r^4\rangle/\langle r^2\rangle^2 = 2$ is a parameter-free test
  that a mode is Gaussian-distributed (it is the "$2.00$" reported in
  [resolution_limits.md](resolution_limits.md)). A pinned or driven mode gives
  less than 2; intermittent dynamics give more.
* Because $r^2$ is exponential, its standard deviation equals its mean: an
  estimate of (1) from $M$ independent frames has relative error $1/\sqrt M$,
  per mode, before any averaging over a $|\mathbf q|$ shell.
* Fitting $\ln\langle r^2\rangle$ against $\ln q$ is correct; fitting
  $\langle \ln r^2\rangle$ is biased low by Euler's constant
  $\gamma_E \approx 0.577$, i.e. by a factor $e^{-\gamma_E}\approx 0.56$ in $\kappa$.
  Average the amplitudes squared first, take the log second.
* A histogram of $\theta_{\mathbf q}$ over frames should be flat. A peak means
  the mode is pinned, by a boundary, a defect, or an unrelaxed initial
  condition, and (1) does not apply to it yet.

### What (1) says

* **No tension:** $\langle r_{\mathbf q}^2\rangle = k_BT/(A\kappa q^4)$, a straight
  line of slope $-4$ on a log–log plot whose intercept fixes $\kappa$.
* **With tension** the spectrum crosses over from $q^{-4}$ to $q^{-2}$ at
  $q_\times = \sqrt{\sigma/\kappa}$.
* **The $1/A$ is a convention artefact.** With the normalisation
  $\tilde h_{\mathbf q} = \int d^2r\, h\, e^{-i\mathbf q\cdot\mathbf r} = A\,r_{\mathbf q}e^{i\theta_{\mathbf q}}$
  the result reads $\langle\tilde r_{\mathbf q}^2\rangle = A\,k_BT/K_q$. The phase
  and the Rayleigh shape are convention-free; only the scale changes.
* **Nothing dynamical enters.** Equation (1) holds for any dynamics that samples
  the canonical ensemble, with or without hydrodynamics; it is a test of the
  energy model and the sampler, not of the mobility.

### Real-space roughness (a convention-free corollary)

Since $\langle h(\mathbf r)^2\rangle = \sum_{\mathbf q}\langle r_{\mathbf q}^2\rangle$
(the cross terms vanish by phase uniformity), and the mode density is
$A\,d^2q/(2\pi)^2$,

$$
\langle h^{2}\rangle = \frac{k_BT}{(2\pi)^{2}} \int \frac{d^{2}q}{\kappa q^{4} + \sigma q^{2}}
= \frac{k_BT}{2\pi} \int_{q_{min}}^{q_{max}} \frac{dq}{q\,(\kappa q^{2}+\sigma)} .
$$

With $q_{min} = 2\pi/L$ (box) and $q_{max} \sim \pi/l_{Face}$ (mesh):

$$
\langle h^{2}\rangle =
\begin{cases}
\dfrac{k_BT}{4\pi\kappa}\Big(\dfrac1{q_{min}^{2}} - \dfrac1{q_{max}^{2}}\Big) \;\approx\; \dfrac{k_BT\, L^{2}}{16\pi^{3}\kappa}, & \sigma = 0,\\[3ex]
\dfrac{k_BT}{4\pi\sigma}\,\ln\!\left[\dfrac{q_{max}^{2}\,(\kappa q_{min}^{2}+\sigma)}{q_{min}^{2}\,(\kappa q_{max}^{2}+\sigma)}\right]
\;\xrightarrow{q_{max}\to\infty}\;
\dfrac{k_BT}{4\pi\sigma}\ln\!\Big(1 + \dfrac{\sigma}{\kappa q_{min}^{2}}\Big), & \sigma > 0 .
\end{cases}
$$

The tensionless roughness grows with the box as $L^2$: the longest mode
dominates, which is why the lowest few $q$ bins carry the most statistical
weight and the worst equilibration problems.

---

## 3. Temporal correlation without hydrodynamics

### The dynamics

With no solvent-mediated coupling between distant points, the only thing the
fluid does is provide a local friction. The sheet then obeys an overdamped
Langevin equation,

$$
\frac{\partial h(\mathbf r,t)}{\partial t}
= -\Lambda\, \frac{\delta E}{\delta h(\mathbf r,t)} + \zeta(\mathbf r,t),
\qquad
\langle \zeta(\mathbf r,t)\,\zeta(\mathbf r',t')\rangle = 2\Lambda k_BT\,\delta^{2}(\mathbf r-\mathbf r')\,\delta(t-t'),
$$

where $\Lambda$ is a mobility per unit area (velocity per force density) and
the noise strength is fixed by the fluctuation–dissipation theorem so that the
stationary distribution is $\propto e^{-E/k_BT}$. From §1,
$\delta E/\delta h = \kappa \nabla^{4} h - \sigma \nabla^{2} h$.

### Mode by mode, in polar form

Fourier transforming, with $\nabla^4 \to q^4$ and $-\nabla^2 \to q^2$, every mode
decouples:

$$
\frac{d}{dt}\big(r_{\mathbf q} e^{i\theta_{\mathbf q}}\big)
= -\Gamma_q\, r_{\mathbf q} e^{i\theta_{\mathbf q}} + \zeta_{\mathbf q}(t),
\qquad
\boxed{\,\Gamma_q = \Lambda\,K_q = \Lambda\,(\kappa q^{4} + \sigma q^{2})\,}
$$

The transformed noise, in the $\tfrac1A\int\cdots$ convention, is a complex white
noise with

$$
\langle \zeta_{\mathbf q}(t)\,\zeta^{*}_{\mathbf q'}(t')\rangle
= \frac{2\Lambda k_BT}{A}\,\delta_{\mathbf q,\mathbf q'}\,\delta(t-t'),
\qquad
\langle \zeta_{\mathbf q}(t)\,\zeta_{\mathbf q'}(t')\rangle = 0,
$$

the second identity because the real-space noise is real and $\mathbf q'\neq-\mathbf q$
when both lie in $\mathcal H$. It is convenient to name the per-component noise
intensity

$$
s_q^{2} \equiv \frac{\Lambda k_BT}{A} = \frac{\Gamma_q}{\beta_q},
$$

so that $\zeta_{\mathbf q}\,dt = s_q\,(dB_1 + i\,dB_2)$ with $B_1, B_2$ independent
standard Wiener processes. In the instantaneous frame of the mode, rotate the
noise onto the radial direction ($\parallel e^{i\theta}$) and the tangential
direction ($\parallel i e^{i\theta}$); the two rotated Wiener increments
$dB_r, dB_\theta$ are again independent and standard because the isotropic noise
has no preferred direction.

### The radial and phase equations

Apply Itô's rule to $r = |h|$ and $\theta = \arg h$. Writing $h = x + iy$ for the
moment to take the derivatives, $\partial^2 r/\partial x^2 + \partial^2 r/\partial y^2 = 1/r$
and $\nabla^2\theta = 0$, so with each component carrying variance $s_q^2\,dt$:

$$
\boxed{\;
d r_{\mathbf q} = \Big(-\Gamma_q\, r_{\mathbf q} + \frac{s_q^{2}}{2\, r_{\mathbf q}}\Big) dt + s_q\, dB_r,
\qquad
d\theta_{\mathbf q} = \frac{s_q}{r_{\mathbf q}}\, dB_\theta
\;}
$$

Read in words:

* **The amplitude relaxes toward a finite value.** The drift vanishes at
  $r^{*2} = s_q^2/(2\Gamma_q) = 1/(2\beta_q)$, which is exactly the mode of the
  Rayleigh distribution of §2. The $s_q^2/2r$ term is the dynamical face of the
  $-k_BT\ln r$ entropic potential: a mode is never driven to zero amplitude.
* **The phase performs a free random walk** with no restoring force (the energy
  did not depend on $\theta$) and with a diffusion coefficient
  $D_\theta = s_q^2/(2 r_{\mathbf q}^2)$ that depends on the instantaneous
  amplitude. Phase moves slowly while the wave is tall and slips rapidly when
  the amplitude passes near zero. There is no single "phase relaxation time":
  the mode decorrelates through the joint motion of $r$ and $\theta$.

### Consistency check: the stationary amplitude is equipartition

Itô on $r^2$: $d(r^2) = 2r\,dr + (dr)^2 = \big(-2\Gamma_q r^2 + 2 s_q^2\big)dt + 2 s_q r\, dB_r$.
Taking the average and setting it stationary,

$$
\langle r_{\mathbf q}^{2}\rangle = \frac{s_q^{2}}{\Gamma_q} = \frac{1}{\beta_q} = \frac{k_BT}{A K_q},
$$

which is (1). The mobility cancels: the static spectrum is independent of
$\Lambda$, the rate is linear in it. (The full stationary density of the pair
$(r,\theta)$ obtained from the Fokker–Planck equation of the two SDEs above is
$\propto r\,e^{-\beta_q r^2}$, the Rayleigh–uniform measure of §2, so the dynamics
reproduces not only the second moment but the whole distribution.)

### The time correlation

What is measured is

$$
C_q(\tau) \equiv \langle h_{\mathbf q}(t)\, h_{\mathbf q}^{*}(t+\tau)\rangle
= \Big\langle r_{\mathbf q}(t)\, r_{\mathbf q}(t+\tau)\; e^{\,i[\theta_{\mathbf q}(t)-\theta_{\mathbf q}(t+\tau)]}\Big\rangle .
$$

It is real: under $\theta\to-\theta$ (the $\pm\mathbf q$ swap, a symmetry of the
dynamics) the exponent changes sign, so the imaginary part averages to zero, and

$$
C_q(\tau) = \big\langle r(t)\, r(t+\tau)\, \cos\Delta\theta \big\rangle,
\qquad \Delta\theta \equiv \theta(t+\tau) - \theta(t).
$$

The quantity $r(t+\tau)\cos\Delta\theta$ is the projection of the later mode onto
the phase direction of the earlier one. To find its evolution, fix the reference
phase $\theta_0 = \theta(t)$ and define the projected coordinate
$u(t') \equiv r(t')\cos[\theta(t') - \theta_0]$, which starts at $u(t) = r(t)$.
Applying Itô to $u = \mathrm{Re}\big(r e^{i(\theta-\theta_0)}\big)$ with the two
SDEs above (or more simply noting that $u$ is a fixed linear combination of the
rotated Cartesian components, for which the $s^2/2r$ and $1/r$ terms cancel
exactly) gives

$$
du = -\Gamma_q\, u\, dt + s_q\, dB_u ,
$$

a one-dimensional Ornstein–Uhlenbeck process with no amplitude-dependent terms
left. Its conditional mean is $\langle u(t+\tau)\,|\,u(t)\rangle = u(t)\, e^{-\Gamma_q\tau}$,
since the noise after $t$ has zero mean and is independent of the state at $t$.
Therefore, conditioning on the state at time $t$,

$$
\big\langle r(t+\tau)\cos\Delta\theta \;\big|\; r(t),\theta(t)\big\rangle = r(t)\, e^{-\Gamma_q\tau},
$$

and multiplying by $r(t)$ and averaging over the stationary distribution,

$$
\boxed{\;
\langle h_{\mathbf q}(t)\, h_{\mathbf q}^{*}(t+\tau)\rangle
= \langle r_{\mathbf q}^{2}\rangle\, e^{-\Gamma_q|\tau|}
= \frac{k_BT}{A\,(\kappa q^{4}+\sigma q^{2})}\;
  \exp\!\big[-\Lambda\,(\kappa q^{4}+\sigma q^{2})\,|\tau|\big]
\;}
\tag{2}
$$

The $|\tau|$ follows from stationarity. The exponential is exact even though
neither $r$ nor $\theta$ separately relaxes exponentially: the decay of the
*projection* is clean because the projection obeys a linear equation, while the
amplitude and the phase individually inherit the nonlinear $1/r$ terms.

### Amplitude-only and phase-only correlations (and why not to use them)

Because the decomposition into $r$ and $\theta$ is explicit, it is tempting to
measure them separately. Their correlations are not single exponentials:

* **Intensity (amplitude squared).** For a circular Gaussian mode the fourth
  moment factorises (Isserlis), giving

  $$
  \langle r^{2}(t)\, r^{2}(t+\tau)\rangle = \langle r^{2}\rangle^{2}\,\big(1 + e^{-2\Gamma_q|\tau|}\big).
  $$

  It decays at **twice** the rate and sits on a constant pedestal equal to the
  square of the spectrum. Fitting it as a plain exponential without subtracting
  $\langle r^2\rangle^2$ and halving the rate overestimates $\Gamma_q$.

* **Phase.** $\langle\cos\Delta\theta\rangle$ is an average of $\cos$ of a random
  walk whose step size depends on $r$; it has no closed exponential form, and
  because phase slips are concentrated at small amplitudes, it is dominated by
  rare excursions and converges slowly.

The projected correlation (2) weights every frame by the amplitude, which is
precisely what removes the small-$r$ phase slips from the estimate. It is the
quantity to fit for $\Gamma_q$.

### Two useful rewritings of (2)

**Relaxation time.** $\tau_q = 1/\Gamma_q = 1/[\Lambda(\kappa q^4+\sigma q^2)]$. With
no tension the longest mode in a box of side $L$ relaxes in
$\tau_{max} = L^4/(16\pi^4\Lambda\kappa)$; the shortest mode on a mesh of spacing
$l$ relaxes in $\sim l^4/(\pi^4\Lambda\kappa)$. Their ratio $(L/l)^4$ is what makes
a free-draining membrane simulation stiff.

**Height–height correlation in real space and time.** Summing (2) over modes,

$$
\langle h(\mathbf r,t)\, h(\mathbf r,t+\tau)\rangle
= \frac{k_BT}{2\pi}\int \frac{q\,dq}{\kappa q^{4}+\sigma q^{2}}\, e^{-\Lambda(\kappa q^{4}+\sigma q^{2})|\tau|},
$$

a stretched decay rather than a single exponential, because every mode has its
own rate.

---

## 4. Mapping to the lattice integrator

SLIMED moves every vertex by

$$
\Delta \mathbf x_i = \frac{D\,\Delta t}{k_BT}\,\mathbf F_i + \sqrt{2 D\,\Delta t}\;\boldsymbol\xi_i,
\qquad \boldsymbol\xi_i \sim \mathcal N(0,1),
$$

i.e. an Euler–Maruyama step of $\dot{\mathbf x}_i = (D/k_BT)\mathbf F_i + \sqrt{2D}\,\boldsymbol\eta_i$,
with a scalar per-vertex mobility $\mu = D/k_BT$ and no mobility matrix.

To identify $\Lambda$, note that the discrete energy is a quadrature of the
continuum one, $E \approx \sum_i a\, e(h_i)$ with $a = A/N$ the area per vertex,
so the force on a vertex is $F_i = -\partial E/\partial h_i = -a\,\delta E/\delta h\,|_{\mathbf r_i}$.
The vertex velocity is then

$$
\dot h_i = \mu F_i = -\frac{D a}{k_BT}\,\frac{\delta E}{\delta h}\Big|_{\mathbf r_i},
\qquad\Longrightarrow\qquad
\boxed{\,\Lambda = \frac{D\,a}{k_BT},\qquad
\Gamma_q = \frac{D\,a}{k_BT}\,(\kappa q^{4} + \sigma q^{2})\,}
$$

The noise is consistent with this: $\langle \eta_i\eta_j\rangle = \delta_{ij}\delta(t-t')$
on the lattice is $a\,\delta^2(\mathbf r-\mathbf r')$ in the continuum, so
$2D\,\delta_{ij} \to 2 D a\,\delta^2 = 2\Lambda k_BT\,\delta^2$, as required by
fluctuation–dissipation. Equations (1) and (2) therefore hold on the mesh with
this $\Lambda$, in the range $q\,l_{Face} \lesssim 2$ where the discrete bending
operator still reproduces $\kappa q^4$ (measured in
[fluctuation_spectrum.md §4b](fluctuation_spectrum.md)). On the lattice, the
polar amplitudes are obtained from the sampled heights by the discrete transform
on the actual zig-zag positions (not an index-space `fft2`), and the phase
convention follows the choice of origin; a shift of origin by $\mathbf d$ adds
$-\mathbf q\cdot\mathbf d$ to every $\theta_{\mathbf q}$ and leaves every result
above unchanged.

Two consequences that matter when comparing runs:

* **$\Gamma_q$ is mesh dependent at fixed `diffConst`.** $a \propto l_{Face}^2$, so a
  coarser mesh relaxes faster in physical time at the same $q$. To keep the
  dynamics while coarsening, hold $D a$ fixed, i.e. scale `diffConst` as
  $l_{Face}^{-2}$ ([resolution_limits.md](resolution_limits.md)).
* **The finite step modifies (1) and (2) slightly.** The Euler–Maruyama map on
  the projected coordinate $u$ is $u_{n+1} = (1-\Gamma_q\Delta t)\,u_n + s_q\sqrt{\Delta t}\,\xi_n$,
  so the exact discrete-time results are
  $\langle r_{\mathbf q}^2\rangle_{\Delta t} = \langle r_{\mathbf q}^2\rangle/(1-\Gamma_q\Delta t/2)$
  and $\langle h_{\mathbf q}(t)h_{\mathbf q}^*(t+n\Delta t)\rangle \propto (1-\Gamma_q\Delta t)^n$,
  stable for $\Gamma_q\Delta t < 2$. The amplitude stays Rayleigh (the map is
  linear and Gaussian), only its scale shifts. The spectrum reads high by
  $\approx\Gamma_q\Delta t/2$ at the zone edge, which is the $\approx 3\%$ per
  0.001 µs seen in the dt sweeps.

---

## 5. What changes with hydrodynamics (for contrast only)

The derivation above assumed the solvent exerts only a local drag. If the
membrane instead moves through an incompressible fluid of viscosity $\eta$, the
mobility becomes non-local and, on a flat sheet, diagonal in $\mathbf q$ with the
Oseen kernel $\Lambda(q) = 1/(4\eta q)$. Equation (1), the Rayleigh amplitude law
and the uniform phase are unchanged, and (2) keeps its exponential form with

$$
\Gamma_{hydro}(q) = \frac{\kappa q^{3} + \sigma q}{4\eta}
$$

instead of $\Lambda(\kappa q^4+\sigma q^2)$; the radial and phase SDEs of §3 hold
with $s_q^2 = \Gamma_q/\beta_q$ evaluated at this $\Gamma_q$. The static spectrum
cannot distinguish the two; the exponent of $\Gamma_q$ versus $q$ (3 versus 4,
tensionless) can. The code has no mobility matrix, so what it integrates is the
free-draining $q^4$ law of §3, and the measured slope of $\approx +3.6$ in
[fluctuation_spectrum.md §6](fluctuation_spectrum.md) is the expected one.

---

## Summary

With $h_{\mathbf q} = r_{\mathbf q} e^{i\theta_{\mathbf q}}$ and
$\beta_q = A(\kappa q^4+\sigma q^2)/k_BT$:

$$
P(r_{\mathbf q},\theta_{\mathbf q}) = \frac{\beta_q}{\pi}\, r_{\mathbf q}\, e^{-\beta_q r_{\mathbf q}^{2}}
\quad\Longrightarrow\quad
\langle h_{\mathbf q} h_{\mathbf q}^{*}\rangle = \langle r_{\mathbf q}^{2}\rangle = \frac{k_BT}{A\,(\kappa q^{4}+\sigma q^{2})},
\qquad \theta_{\mathbf q}\ \text{uniform};
$$

$$
dr = \Big(-\Gamma_q r + \frac{s_q^2}{2r}\Big)dt + s_q\,dB_r,\quad
d\theta = \frac{s_q}{r}\,dB_\theta
\quad\Longrightarrow\quad
\langle h_{\mathbf q}(t)\,h_{\mathbf q}^{*}(t+\tau)\rangle
= \langle r(t)\,r(t+\tau)\cos\Delta\theta\rangle
= \frac{k_BT}{A\,(\kappa q^{4}+\sigma q^{2})}\,
  \exp\!\Big[-\frac{D a}{k_BT}\,(\kappa q^{4}+\sigma q^{2})\,|\tau|\Big].
$$

The first follows from the diagonal quadratic energy together with the polar
Jacobian, which supplies the second half of $k_BT$ that the phase does not carry;
the second from the overdamped Langevin equation with local mobility
$\Lambda = Da/k_BT$, whose projected coordinate is a linear Ornstein–Uhlenbeck
process even though the amplitude and phase separately are not.
