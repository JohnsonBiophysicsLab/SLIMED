


# SLIMED

<p align="center">
  <!-- Build / CI -->
  <!--a href="https://github.com/JohnsonBiophysicsLab/SLIMED/actions">
    <img alt="CI" src="https://img.shields.io/github/actions/workflow/status/JohnsonBiophysicsLab/SLIMED/ci.yml?branch=main&label=CI">
  </a>
  <!-- Codecov (or Coveralls) -->
  <a href="https://codecov.io/gh/JohnsonBiophysicsLab/SLIMED">
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
  <img alt="CMake" src="https://img.shields.io/badge/CMake-%E2%89%A53.16-informational">
  <!-- Package managers (optional) -->
  <!-- img alt="Conan" src="https://img.shields.io/badge/Conan-ready-0ea5e9">
  <!-- img alt="vcpkg" src="https://img.shields.io/badge/vcpkg-port-22c55e">
  <!-- Platforms -->
  <img alt="Platforms" src="https://img.shields.io/badge/Linux%20|%20macOS%20-supported-success">
  <!-- Code style -->
  <!-- img alt="clang-format" src="https://img.shields.io/badge/clang--format-enforced-brightgreen"-->
</p>


SLIMED (Subdivision-limit membrane dynamics) model is established with triangular mesh and optimized using an energy function. The installation instructions and dependencies are provided in the [Installation](#installation) section below. The first step requires setting up the triangular mesh model to approximate the membrane's geometry applying Loop's subdivision method. The lowest energy search model minimizes the membrane energy evaluated by the energy function. The Brownian Dynamics model simulates the moving membrane's surface with a displacement equation, performed on the limit surface with the help of a conversion matrix. Three types of boundary conditions are available: Fixed, Periodic, and Free, all defined for "ghost vertices" on the boundary of the triangular mesh.

## Installation
To install the model, follow these steps:

1. Clone this repository to your local machine. You can use ``git`` as follows or download zip from github webpage.

```console
git clone https://github.com/mjohn218/continuum_membrane.git
```

2. Install dependencies (GSL and OpenMP)

This project depends on:

* **CMake** 3.16 or newer
* **GSL** (GNU Scientific Library)
* **OpenMP** -- required for the parallel build; its headers are needed even for
  a serial build, because several headers include `<omp.h>` unconditionally
* **GoogleTest** (optional, for the unit tests)

### Ubuntu (apt)

```bash
# Update package lists
sudo apt-get update

# GSL (headers + library)
sudo apt-get install -y libgsl-dev

# Compilers
sudo apt-get install -y build-essential clang

# Build system
sudo apt-get install -y cmake

# Unit tests (optional)
sudo apt-get install -y libgtest-dev

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

# Build system
brew install cmake

# Unit tests (optional)
brew install googletest

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

### Quickly Compile and Run with CMake Presets

SLIMED builds with CMake. Compile without OpenMP and run membrane energy minimization:

```console
cmake --workflow --preset serial
./build/bin/continuum_membrane
```

Compile with OpenMP and run membrane energy minimization:

```console
cmake --workflow --preset omp
./build-omp/bin/continuum_membrane
```
Other versions of executables are present in the respectively `bin` directory.

| Executable | Purpose | 
| --- | --- | 
| `continuum_membrane` | Membrane energy minimization | 
| `continuum_membrane_multithreading` | Membrane energy minimization with embarassingly parallelization |
| `membrane_dynamics` | Membrane dynamics simulation |
| `membrane_dynamics_multithreading` | Membrane dynamics simulation with embarassingly parallelization |


### Presets (quick way to build and switch between serial and OpenMP)

Because serial and OpenMP builds are both used routinely, `CMakePresets.json`
defines them as named presets, each with its own build directory:

| Preset | Build directory | Builds | Equivalent to |
| --- | --- | --- | --- |
| `serial` | `build/` | simulation programs | `make serial` / `make dyna` |
| `omp` | `build-omp/` | simulation programs | `make omp` / `make dyna_omp` |
| `coverage` | `build-cov/` | programs + unit tests | `make test COVERAGE=1` |

`serial` and `omp` build the four simulation programs only. They set
`SLIMED_BUILD_TESTS=OFF`, so GoogleTest is never looked for and does not need to
be installed to use them.

Configure and build in one command:

```console
cmake --workflow --preset omp
cmake --workflow --preset serial
```

Or drive the steps separately:

```console
cmake --preset omp           # configure
cmake --build --preset omp   # build (all cores)
```

To run the unit tests, either use the `coverage` preset or configure a build
directory by hand -- `SLIMED_BUILD_TESTS` defaults to `ON`, so a plain
configure includes them:

```console
cmake -S . -B build-test
cmake --build build-test -j8
ctest --test-dir build-test --output-on-failure
```

Because each preset owns a separate build directory, switching back and forth
recompiles nothing -- both trees persist side by side. Setting
`-DSLIMED_ENABLE_OPENMP` on a single build directory works too, but changes the
`-DOMP` define for every target and so triggers a full rebuild each time you
flip it.

To list what is available:

```console
cmake --list-presets
```

Presets require CMake 3.25 or newer. With an older CMake, use the `-D` options
below instead -- they are exactly what the presets set. For personal presets
that are not shared, create a `CMakeUserPresets.json`; it is gitignored.

### Configure and build separately

Configure once into a build directory, then build:

```console
cmake -S . -B build
cmake --build build -j8
```

`cmake` locates GSL, sets the C++14 standard, resolves OpenMP, and prints a
summary of what it found. The executables are written to `build/bin/`.

To rebuild after editing code, only the second command is needed:

```console
cmake --build build -j8
```

The configure step re-runs itself when `CMakeLists.txt` changes, and sources are
globbed with `CONFIGURE_DEPENDS`, so adding or removing a `.cpp` under `src/` is
picked up automatically.

### Build options

Options are cached in the build directory, so they are set once rather than on
every invocation.

| Option | Default | Effect |
| --- | --- | --- |
| `SLIMED_ENABLE_OPENMP` | `OFF` | Compile the OpenMP code paths (`-fopenmp` and `-DOMP` together). |
| `SLIMED_ENABLE_COVERAGE` | `OFF` | Instrument for `gcov`/`lcov`. |
| `SLIMED_BUILD_TESTS` | `ON` | Build the GoogleTest suite. |

Standard CMake variables work as usual:

```console
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-march=native"
```

`CMAKE_BUILD_TYPE` defaults to `Release` (`-O3`), matching the hand-written
Makefile. To inspect what a build directory is currently set to:

```console
cmake -LH -N -B build | grep -A1 SLIMED_
```

### The four programs

One build produces all four; each replaces one goal of the hand-written Makefile.

| Program | Old target | What it does |
| --- | --- | --- |
| `continuum_membrane` | `make serial` / `make omp` | Lowest-energy conformation search |
| `continuum_membrane_multithreading` | `make multi` | Same, embarrassingly parallel over parameter sets |
| `membrane_dynamics` | `make dyna` / `make dyna_omp` | Membrane Brownian dynamics |
| `membrane_dynamics_multithreading` | `make dyna_multi` | Same, embarrassingly parallel |

Serial and OpenMP are no longer separate targets. OpenMP is a configure-time
choice, so the same program name covers both:

```console
# Serial
cmake -S . -B build
cmake --build build -j8

# OpenMP
cmake -S . -B build-omp -DSLIMED_ENABLE_OPENMP=ON
cmake --build build-omp -j8
```

These two are what the `serial` and `omp` presets above do, so
`cmake --workflow --preset omp` is the shorter equivalent.

Flipping `SLIMED_ENABLE_OPENMP` on an existing build directory also works, but
it recompiles everything -- keeping two build directories avoids that.

To build a single program instead of all four:

```console
cmake --build build --target membrane_dynamics
```

### Running the programs

The programs read `input.params` from the **current working directory**, so run
them from the top of the source tree rather than from the build directory:

```console
./build/bin/continuum_membrane
./build/bin/membrane_dynamics
```

Running them from inside `build/bin` prints `There was a problem opening the
parameter file!` and silently falls back to built-in defaults.

### OpenMP notes

`-fopenmp` and `-DOMP` are always enabled together and never separately. Most
`#pragma omp parallel for` directives in `src/` are unguarded, but the loop
bodies that need a thread id read it from `omp_get_thread_num()` only under
`#ifdef OMP` and otherwise hard-code `0`, so enabling the flag without the
define would race on the per-thread accumulation buffers.

Serial builds still need `omp.h` on the include path, because `Mesh.hpp`,
`Dynamics.hpp`, `Force.hpp` and `Gauss_quadrature.hpp` include it
unconditionally. CMake attaches the OpenMP include directories without the
`-fopenmp` flag to cover this.

On macOS, Apple Clang has no bundled OpenMP runtime and Homebrew keeps `libomp`
keg-only. CMake retries a failed probe against `brew --prefix libomp`
automatically, so `brew install libomp` is normally all that is required.
Homebrew GCC also works for the simulation programs:

```console
cmake -S . -B build-gcc -DCMAKE_CXX_COMPILER=g++-15
```

Note that the compiler cannot be changed in an existing build directory --
delete it and configure again. Under Homebrew GCC the unit tests will not link,
because Homebrew's `libgtest.a` is built against libc++ while GCC uses
libstdc++; use the default Apple Clang when you need the tests.

### Tests

`SLIMED_BUILD_TESTS` defaults to `ON`, so any build directory configured by hand
includes the unit tests:

```console
cmake -S . -B build-test
cmake --build build-test -j8
ctest --test-dir build-test --output-on-failure
```

Note that the `serial` and `omp` presets deliberately set it to `OFF`, so
`build/` and `build-omp/` contain no tests; use a separate directory as above,
or the `coverage` preset.

Each GoogleTest case is registered with CTest individually. Several tests open
fixture files by paths relative to the working directory, so CTest runs
`test_main` from the top of the source tree automatically. The binary can also
be run directly, from the source root:

```console
./build-test/bin/test_main
```

If GoogleTest is not installed, CMake prints a warning and skips the test
target; the simulation programs build either way. Configure with
`-DSLIMED_BUILD_TESTS=OFF` to silence the warning.

For coverage:

```console
cmake -S . -B build-cov -DSLIMED_ENABLE_COVERAGE=ON
cmake --build build-cov -j8
ctest --test-dir build-cov
```

or simply `cmake --workflow --preset coverage`.

`.gcda`/`.gcno` files are left in the build tree for `lcov` or `gcov` to
collect. There is no packaged `coverage-html` target yet.

### Cleaning

```console
cmake --build build --target clean   # like "make clean"; keeps the configuration
rm -rf build                         # full reset, forgets all cached options
rm -rf build build-omp build-cov build-test   # remove all build directories
```

### The previous Makefile

The hand-written Makefile is still in the tree as `Makefile.legacy` and still
works; CMake neither uses nor interferes with it:

```console
make -f Makefile.legacy serial
./bin/continuum_membrane
```

Its `mpi` target is not carried over to the CMake build: it only added `-DMPI`,
which no source file tests, and its `MPCC` variable was never used by any rule.

There are no `install` or packaging targets in the CMake build yet -- the
programs are run from the build tree, as they were with the Makefile.

## Input Parameters

The input file is `continuum_membrane/input.params`. Parameters are broken down into geometric parameters, physical properties, insertion mode, and advanced parameters.

| Category    | Subcategory / Parameter   | Name / Key              | Description |
|-------------|---------------------------|-------------------------|-------------|
| **Geometric Parameters** | —                         | `lMeshSide`             | Target side length of the triangular mesh (nm). This only serves as a reference scale. The mesh side length set up by the algorithm may vary. |
|             | **Sphere model**          | `isSphere`               | Set `true` to enable sphere mode. |
|             |                           | `rSphere`                | Target radius of sphere (nm). This is the radius of spherical frame to set up the triangular mesh. The radius of the resulting membrane represented by the triangular mesh may vary. |
| **Physical Properties** | —                         | `c0Insertion`            | Curvature of the membrane at the insertion area. |
|             | —                         | `c0Membrane`              | Spontaneous curvature of the membrane. |
|             | —                         | `kcMembraneBending`       | Membrane bending constant in the energy function (pN·nm). |
|             | —                         | `usMembraneStretching`    | Membrane stretching modulus in the energy function (pN/nm). |
|             | —                         | `uvVolumeConstraint`      | Volume constraint coefficient in the energy function (pN/nm²). |
| **Insertion Mode** | —                  | `isInsertionIncluded`     | Set `true` to include insertion. |
|             | —                         | `sigma`                   | 2·sigma (nm) is the length scale of decaying insertion curvature, or in other words expansion of non-spontaneous curvature due to insertion. |
| **Advanced Parameters** | **Optimization** | `numMaxIterations`        | Number of maximum iterations allowed. |
|             |                           | `criterionForce`          | Force criteria to determine if adequate optimization is accomplished (pN). |
|             | **Algorithm**             | `gaussQuadratureN`        | Default Gauss Quadrature used in integral approximation. |


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

A detailed doxygen style documentation is available in [SLIMED developer's documentation](https://johnsonbiophysicslab.github.io/SLIMED/index.html).

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
