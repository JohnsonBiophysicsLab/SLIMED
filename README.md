


# SLIMED

<p align="center">
  <!-- Build / CI -->
  <!--a href="https://github.com/JohnsonBiophysicsLab/SLIMED/actions">
    <img alt="CI" src="https://img.shields.io/github/actions/workflow/status/JohnsonBiophysicsLab/SLIMED/ci.yml?branch=main&label=CI">
  </a>
  <!-- Codecov (or Coveralls) -->
  <!--a href="https://codecov.io/gh/JohnsonBiophysicsLab/SLIMED">
    <img alt="Coverage" src="https://img.shields.io/codecov/c/github/JohnsonBiophysicsLab/SLIMED?logo=codecov">
  </a>
  <!-- Release -->
  <a href="https://github.com/JohnsonBiophysicsLab/SLIMED/releases">
    <img alt="Release" src="https://img.shields.io/github/v/release/JohnsonBiophysicsLab/SLIMED">
  </a>
  <!-- License -->
  <a href="https://github.com/JohnsonBiophysicsLab/SLIMED/blob/main/LICENSE">
    <img alt="License" src="https://img.shields.io/github/license/JohnsonBiophysicsLab/SLIMED">
  </a>
  <!-- C++ standard -->
  <img alt="C++" src="https://img.shields.io/badge/C%2B%2B-%E2%89%A514-blue">
  <!-- CMake minimum -->
  <img alt="CMake" src="https://img.shields.io/badge/CMake-%E2%89%A53.20-informational">
  <!-- Package managers (optional) -->
  <!-- img alt="Conan" src="https://img.shields.io/badge/Conan-ready-0ea5e9">
  <!-- img alt="vcpkg" src="https://img.shields.io/badge/vcpkg-port-22c55e">
  <!-- Platforms -->
  <img alt="Platforms" src="https://img.shields.io/badge/Linux%20|%20macOS%20-supported-success">
  <!-- Code style -->
  <!-- img alt="clang-format" src="https://img.shields.io/badge/clang--format-enforced-brightgreen"-->
</p>


SLIMED (Subdivision-limit membrane dynamics) model is established with triangular mesh and optimized using an energy function. The installation instructions and dependencies are provided in README.rst. The first step requires setting up the triangular mesh model to approximate the membrane's geometry applying Loop's subdivision method. The lowest energy search model minimizes the membrane energy evaluated by the energy function. The Brownian Dynamics model simulates the moving membrane's surface with a displacement equation, performed on the limit surface with the help of a conversion matrix. Three types of boundary conditions are available: Fixed, Periodic, and Free, all defined for "ghost vertices" on the boundary of the triangular mesh.

## Installation
To install the model, follow these steps:

1. Clone this repository to your local machine. You can use ``git`` as follows or download zip from github webpage.

```console
git clone https://github.com/mjohn218/continuum_membrane.git
```

2. Install dependencies (GSL and OpenMP)

This project depends on:

* **GSL** (GNU Scientific Library)
* **OpenMP** (only for the `omp`/parallel build)

### Ubuntu (apt)

```bash
# Update package lists
sudo apt-get update

# GSL (headers + library)
sudo apt-get install -y libgsl-dev

# Compilers
sudo apt-get install -y build-essential clang

# OpenMP
# - If you build with GCC: OpenMP is included (use -fopenmp, links against libgomp).
# - If you build with Clang: install the LLVM OpenMP runtime:
sudo apt-get install -y libomp-dev
```

### macOS (Homebrew)

```bash
# Install Homebrew if needed: https://brew.sh

# GSL
brew install gsl

# Option A (recommended): use Homebrew LLVM (Clang) for OpenMP
brew install llvm libomp
# Option B: Apple Clang + libomp can work, but setup is trickier.
```

**Using Homebrew LLVM for OpenMP (recommended):**

```bash
# Determine prefixes (Apple Silicon uses /opt/homebrew, Intel uses /usr/local)
brew --prefix llvm
brew --prefix libomp
```

## Compile and Run

To compile and run continuum membrane lowest energy conformation search with OpenMP parallelization. See section 5 for model details.

OpenMP parallelization:

```console
make omp
./bin/continuum_membrane
```

Embarrasingly parallel for multiple parameter sets:

```console
make multi
./bin/continuum_membrane_multithreading
```

No parallel:

```console
make serial
./bin/continuum_membrane
```

To compile and run membrane Brownian dynamics simulation. See section 6 for model details.

```console
make dyna
./bin/membrane_dynamics
```

## Parameters

Input parameter file is located in the root directory: ``continuum_membrane/input.params``

## Energy Function and Lowest Energy Search

The goal for the lowest energy search model is to minimize the membrane energy evaluated by the energy function, which is the sum of membrane bending energy, area constraint energy (or elastic area change energy), and volume constraint energy:


$$
E = E_B + E_S + E_V = \int_S \frac{1}{2}\kappa (2H-C_0)^2 dS + \frac{1}{2} \mu_S \frac{(S-S_0)^2}{S_0} + \frac{1}{2} \mu_V \frac{(V-V_0)^2}{V_0}
$$

where:

- $\kappa$ : membrane bending constant ``kcMembraneBending``
- $H$ : mean membrane culvature
- $C_0$ : spontaneous curvature of the membrane ``c0Membrane``
- $\mu_S$ : membrane streching modulus ``usMembraneStretching``
- $S$ : global membrane area
- $S_0$ : target membrane area
- $\mu_V$ : volume constraint coefficient ``uvVolumeConstraint``
- $V$ : global volume
- $V_0$ : target volume
 
## Membrane Brownian Dynamics

Membrane Brownian Dynamics model runs a step-wise simulation of the moving membrane surface with the following equation:



$$
\Delta X = -\frac{D\Delta t}{k_b T} \nabla E + \sqrt{2D\Delta t} (N(0,1))
$$

where:

- $\Delta X$: displacement of point on limit surface
- $D$: diffusion constant of the membrane
- $\Delta t$: time step
- $N(0,1)$: standard normal distribution

Note that the displacement of membrane according to the equation above is performed on the limit surface, not the control mesh.
In this case, a conversion matrix helps to convert between triangular mesh and limit surface, as currently the points on the limit surface
represented by the mesh point are chosen to represent the surface.

$$
M_{s} = C M_{m}
$$

$$
M_{m} = C^{-1} M_{s}
$$

## Boundary Conditions

Three types of boundary conditions are provided currently in both models. Note that "ghost vertices" are defined as points on the boundary of the triangular mesh that only serve to provide reference when calculating limit surface on the boundary, as calculating position of a point on the limit surface require the coordinates of 12 neighboring vertices (if regular). However, the "ghost vertices" themselves do not correspond to real points on the surface.

- Fixed: 2 rings of ghost vertices are fixed in space
- Periodic: 3 rings of ghost vertices that mimics the movement of the vertices on the opposite side of the membrane.
- Free: 2 rings of ghost vertices are generated after movement by forming parallelogram extend from the real points on the control mesh


## For Developers

A detailed doxygen style documentation is available in ``continuum_membrane/html/index.html``. We are working on getting the documentation online, but currently you will need to clone the repository and view it locally.

## Cite Continuum Membrane

If you use or modify continuum membrane model, in addition to citing
NERDSS, please be kind and cite us:

1. For continuum Membrane Implementation: Fu, Y., Yogurtcu, O.N., Kothari,
R., Thorkelsdottir, G., Sodt, A.J. & Johnson, M.E. (2019) An implicit
lipid model for efficient reaction-diffusion simulations of protein binding to surfaces of arbitrary topology. *J Chem Phys.* 151 (12), 124115. 
doi:[10.1063/5.0045867](https://doi.org/10.1063/5.0045867)

3. For membrane energies and insertion: Fu, Y., Zeno, W., Stachowiak, J. &
Johnson, M.E. A continuum membrane model predicts curvature sensing by
helix insertion. Soft Matter. 2021 Dec 8;17(47):10649–10663
doi:[10.1039/d1sm01333e](https://doi.org/10.1039/d1sm01333e)

## License
This project is licensed under the GPL License - see the LICENSE file for details.

## Reference
1. Helfrich, W. (1973). Elastic properties of lipid bilayers: theory and possible experiments. Zeitschrift fur Naturforschung C, 28(11), 693-703.
