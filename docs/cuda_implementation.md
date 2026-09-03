# The CUDA force backend: implementation, performance and consequences

This document describes how the GPU path is built, where each piece lives, why
a speedup is expected, what was actually measured on hardware, and what the
design costs you. It is the reference companion to
[`gpu_acceleration.md`](gpu_acceleration.md), which records the profiling that
motivated the work and the history of the rewrite.

Everything quantitative below was measured on a real device. The hardware and
method are in [section 7](#7-benchmarking).

---

## 1. What is on the GPU, and what is not

Only one thing is on the device: the per-face and per-vertex work inside
`Mesh::Compute_Energy_And_Force()`. That is where the profiler said 95.3% of
wall time went, so it is the only part worth moving.

Everything else — the NCG optimiser, the line search that drives it, mesh
setup, I/O, the area and volume totals — stays on the host. In particular the
optimiser is unchanged, which is why `forceBackend` can be switched without
touching the model.

Selection is at runtime, in `input.params`:

```
forceBackend = cpu      # the CPU loops (default)
forceBackend = gpu      # the CUDA backend; fails the run if no device is usable
forceBackend = auto     # the GPU when one is usable, the CPU otherwise, saying which
```

`Mesh::resolve_force_backend()` in
`src/energy_force/Compute_energy_and_force_on_mesh.cpp` reads it once and
caches the choice. `gpu` **throws** rather than falling back — see
[section 8.3](#83-gpu-fails-loudly).

## 2. File map

| file | role |
| --- | --- |
| `include/energy_force/Patch_kernel.{hpp,inl}` | the patch energy/force and area/volume arithmetic, `__host__ __device__` |
| `include/energy_force/Patch_rows_flat.hpp` | shape-function rows repacked into contiguous blocks |
| `include/cuda/Device_mesh_layout.hpp` | the mesh as flat arrays, plus the gather maps |
| `include/cuda/Force_kernels.{hpp,inl}` | one kernel body per stage, indexed by face or vertex |
| `include/cuda/Host_force_backend.hpp` | runs those bodies in a loop — the fallback and the test oracle |
| `include/cuda/Cuda_force_backend.hpp` | the device driver's interface (no CUDA types leak out) |
| `src/cuda/Cuda_force_backend.cu` | the only file nvcc compiles: kernels, allocation, transfer, launch |
| `src/cuda/Cuda_force_backend_stub.cpp` | stands in when CUDA is not compiled in |

CMake compiles **exactly one** of the last two. `SLIMED_ENABLE_CUDA=ON` swaps
the stub out for the `.cu` and defines `SLIMED_WITH_CUDA`, so call sites need
no `#ifdef`.

## 3. The three-layer design

The central idea is that **the physics exists once in source and is compiled
twice**. There is no host copy of the kernel and a device copy to keep in sync.

### 3.1 Layer one — the kernel bodies

`Patch_kernel.hpp` defines

```cpp
#ifdef __CUDACC__
#define SLIMED_HD __host__ __device__
#else
#define SLIMED_HD
#endif
```

and every arithmetic routine is `SLIMED_HD inline`, written over fixed-size
stack arrays (`double[3]`, `double[kShapeRows][3]`) with no allocation, no STL
and no library call. `kMaxControlPoints = 14` and `kShapeRows = 7` bound every
buffer, so a patch evaluation touches only registers and local memory.

`Force_kernels.hpp` raises this to the stage level. Each stage is a function
taking **one index**:

```cpp
SLIMED_HD inline void area_volume_for_face   (const ForceKernelArgs&, int face);
SLIMED_HD inline void patch_force_for_face   (const ForceKernelArgs&, int face);
SLIMED_HD inline void regularization_for_face(const ForceKernelArgs&, int face);
SLIMED_HD inline void gather_force_for_vertex(const ForceKernelArgs&, int vertex);
SLIMED_HD inline void gather_regularization_for_vertex(const ForceKernelArgs&, int vertex);
```

That signature is the whole trick: a loop can supply the index, and so can
`blockIdx.x * blockDim.x + threadIdx.x`.

Bodies live in `.inl` files included at the end of the headers, so nvcc sees
the definitions in every translation unit and needs no device linking
(`CUDA_SEPARABLE_COMPILATION OFF`).

### 3.2 Layer two — the data layout

`Mesh`, `Face` and `Vertex` are pointer-chasing structures: `std::vector`
members and one heap-allocated `gsl_matrix` per 3-vector. None of that can
cross to a device.

`DeviceMeshLayout::build()` flattens the topology into contiguous arrays with
explicit index maps, and separates two lifetimes:

- **Topology** — one-ring lists as CSR, a per-face `FacePatchDescriptor`
  (`kind`, `nControlPoints`, `oneRingOffset`, `nChildren`), valences, ghost and
  boundary flags, corner lists, and the vertex gather maps. Fixed while
  connectivity is, so it is uploaded **once**.
- **Coordinates** — the only thing a line-search trial changes, and so the only
  thing re-uploaded per evaluation.

`PatchRowsFlat` does the same for the shape functions, which are built once at
setup and then never change: it repacks `std::vector<Matrix>` into row-major
`double` blocks that the kernel indexes directly.

Ghost and boundary flags are kept **separate** rather than folded into
`PatchKind`, because the two passes skip on different ones:
`calculate_element_area_volume()` skips ghost faces, `Compute_Energy_And_Force()`
skips boundary faces, and a boundary face with a complete one-ring does
contribute area while contributing no energy.

`DeviceMeshLayout` deliberately has **no CUDA dependency at all**. It is plain
host C++, so the index arithmetic — where a port like this actually goes wrong
— is unit-tested on the host against the production evaluation before nvcc is
involved.

### 3.3 Layer three — the two drivers

`HostForceBackend` runs the stages in `for` loops. `Cuda_force_backend.cu`
wraps each in a `__global__`. The `.cu` is consequently thin: four three-line
kernels, and otherwise allocation, transfer and launch.

```cpp
__global__ void patch_force_kernel(ForceKernelArgs args)
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if (face < args.nFaces) { patch_force_for_face(args, face); }
}
```

`ForceKernelArgs` is a plain-old-data struct of pointers and scalars, passed to
the kernel **by value**. The same struct describes host memory or device memory
unchanged; only the driver knows which.

## 4. Execution

`kBlockSize = 128`, chosen for register pressure: a patch body holds a dozen
3-vectors plus up to a 14x3 control net, so a larger block trades occupancy for
spills. Grid is `(n + 127) / 128`.

One evaluation is four launches with a host round trip in the middle:

| # | kernel | threads | writes |
| --- | --- | --- | --- |
| 1 | `area_volume_kernel` | one per face | `faceArea`, `faceVolume` |
| — | *host*: download both, sum in face order, set `param.area` / `param.vol` | | |
| 2 | `patch_force_kernel` | one per face | face energy, curvature, normal, and the face's own force **slots** |
| 3 | `regularization_kernel` | one per face | `faceERegular`, `faceDeformCase`, three corner slots |
| 4 | `gather_kernel` | one per vertex | the four per-vertex force vectors |

The host round trip after stage 1 is required, not incidental: the area and
volume **constraint** forces computed in stage 2 need the global totals, which
are a reduction over stage 1's output. See
[section 8.2](#82-the-mid-evaluation-host-reduction).

Every launch is followed by `check_launch()`, which calls `cudaGetLastError()`
and then `cudaDeviceSynchronize()`, throwing with the CUDA error text on
failure. Silently ignoring a `cudaMemcpy` return is how a GPU port ends up
producing plausible garbage. The sync before stage 4 is also what guarantees
the gather sees every scatter.

## 5. Force accumulation: scatter/gather, not atomics

This is the one part a GPU cannot do the way the CPU does.

The CPU keeps a **full-length force buffer per OpenMP thread** and sums them at
the end. With thousands of GPU threads that is not an option — it would be
thousands of full-length buffers.

The obvious replacement, `atomicAdd`, was rejected: floating-point addition is
not associative, so the result would depend on the order the hardware happened
to schedule the atomics, and **the same input would not reproduce the same
output**.

Instead accumulation is split in two:

- **Scatter.** Each face writes only *its own* slots, into a scratch buffer
  laid out exactly like `oneRingIndices()` — one 3-vector per
  (face, control point). No two faces ever touch the same slot, so the write
  needs no synchronisation at all.
- **Gather.** Each vertex sums the slots listed for it by
  `vertexSlotOffsets()` / `vertexSlots()`, in a fixed order.

Both passes are data-race free by construction, and the summation order is a
property of the layout rather than of the schedule. **The GPU result is
therefore bit-identical run to run** — measured, see
[section 8.1](#81-reproducibility-and-bitwise-agreement).

The cost is memory: the scratch is `nSlots * 3` doubles per force term, where
`nSlots = sum of nControlPoints ~= 12 * nFaces`.

## 6. Why a speedup is expected

**The work is embarrassingly parallel.** Each face's patch evaluation is
independent — it reads its own one-ring and shared read-only shape functions,
and writes its own slots. There is no cross-face dependency until the gather.

**The arithmetic intensity is high.** Each face performs a
`kShapeRows x nCtrl` by `nCtrl x 3` product per quadrature sample plus the
curvature algebra, all out of registers, against a one-ring of at most 14
vertices. One force evaluation on the default mesh is roughly 10^8 flops.

**The prerequisite rewrite had already happened.** Before the flat-array work,
about 70% of runtime was `malloc`/`free` and `gsl_matrix` accessor calls on
3-element vectors — none of which is device code. Flattening the kernels was
both the largest single CPU win available (11.3x serial) *and* the thing that
made a port possible at all. The GPU inherits an already-good kernel, which is
also why the measured speedups below are modest rather than the 50x one might
expect against an unoptimised baseline.

### 6.1 ...and why the speedup is bounded

Three structural limits, all visible in the measurements:

1. **Per-call transfers.** The line search calls the force evaluation on the
   order of 250 times per iteration, changing only coordinates. Each call
   uploads coordinates and reference coordinates and downloads the results.
   With `V` vertices and `F` faces that is `48V` bytes up and
   `(6F + 12V) * 8 + F` bytes down, plus `16F` for the mid-evaluation
   reduction — about **1.2 MB per evaluation** on the default 8,400-face mesh.
2. **Four full device syncs per evaluation.** `check_launch()` synchronises
   after every kernel, so launch latency is paid four times and never
   overlapped.
3. **The host-side reduction** forces a device-to-host-to-device turnaround in
   the middle of every evaluation.

None of these grow with mesh size as fast as the kernel work does, which is
exactly why the measured advantage widens from 1.21x to 1.55x between the two
mesh sizes below.

## 7. Benchmarking

### 7.1 Method and hardware

| | |
| --- | --- |
| GPU | NVIDIA GeForce RTX 4050 Laptop, compute capability 8.9, 20 SMs, 6 GB |
| nvcc | 13.3.73 |
| Host | GCC 15.2, GSL 2.8, Intel i5-13450HX (8C/16T), WSL2 |
| Build | `-DSLIMED_ENABLE_CUDA=ON -DSLIMED_CUDA_ARCHITECTURES=89 -DSLIMED_ENABLE_OPENMP=ON`, Release |

Timings are whole-program wall clock for `continuum_membrane`, with
`meshpointOutput` and `xyzOutput` off so file I/O does not enter the
measurement, run outside any synced folder.

**Thread count matters, and 8 is this CPU's optimum.** Benchmarking the GPU
against a badly configured CPU would be meaningless:

| `OMP_NUM_THREADS` | 1 | 2 | 4 | 8 | 16 |
| --- | --- | --- | --- | --- | --- |
| default mesh, 200 iterations | 8.41 s | 4.63 s | 3.23 s | **3.02 s** | 4.10 s |

16 threads is slower than 8 — the cores are saturated and hyperthreading only
adds contention. All comparisons below use 8.

### 7.2 Results

| mesh | iterations | `cpu` | `gpu` | speedup |
| --- | --- | --- | --- | --- |
| 8,399 faces (default input) | 200 | 3.01 s | **2.49 s** | 1.21x |
| 99,359 faces (`sideX = sideY = 1035`) | 500 | 86.57 s | **55.76 s** | 1.55x |

Energy traces are identical between backends at every reported step at both
sizes.

Note that this contradicts the original expectation recorded in
`gpu_acceleration.md`, that the default 8,400-face input would *lose* to the
multi-threaded CPU. On this hardware the crossover is below 8,400 faces. That
prediction was made against a 10-thread Apple M5; this is an 8-core i5 against
a laptop RTX 4050, so the crossover is hardware-dependent and worth
re-measuring rather than assuming.

### 7.3 Reproducing this

```bash
cmake -B build -DSLIMED_ENABLE_CUDA=ON -DSLIMED_CUDA_ARCHITECTURES=<your cc> \
      -DSLIMED_ENABLE_OPENMP=ON
cmake --build build -j

# always verify before trusting a new machine -- see section 9
ctest --test-dir build -R CudaForceBackend

# then, with forceBackend = cpu and = gpu in input.params:
OMP_NUM_THREADS=8 ./build/bin/continuum_membrane
```

For the large mesh set `sideX = sideY = 1035.0` with `lFace = 5.0`. Use enough
iterations that the loop dominates fixed setup; 200 at the default size and 500
at 10^5 faces are sufficient.

### 7.4 A benchmarking trap worth knowing

Before `set_adjacent_faces_of_faces()` was rewritten to use an edge map, mesh
setup was **O(N^2) in face count** and cost ~90 s at 99k faces. Because that
cost is backend-independent, it swamped the comparison:

| 99,359 faces, 500 iterations | cpu | gpu | apparent speedup |
| --- | --- | --- | --- |
| with O(N^2) setup | 167.5 s | 143.6 s | 1.17x |
| per force-evaluation step only | 153.9 ms | 109.4 ms | **1.41x** |
| with O(N) setup (current) | 86.6 s | 55.8 s | **1.55x** |

The lesson generalises: **when a fixed serial cost dominates, whole-program
wall clock understates the backend difference.** If you change mesh size,
re-check that setup is not dominating before drawing conclusions — compare two
iteration counts and take the marginal cost per iteration. At 10^5 faces with
the old setup, a 1-iteration run and a 30-iteration run differed by less than
the run-to-run noise.

## 8. Consequences

### 8.1 Reproducibility and bitwise agreement

Two different questions, with two different answers.

**Run-to-run reproducibility — excellent.** Comparing every output as raw
IEEE-754 bit patterns over 119,174 values on the default mesh:

| comparison | bit-identical |
| --- | --- |
| GPU vs GPU | **100.00%** |
| CPU vs CPU (1/2/4/8/16 threads) | **100.00%** |

The scatter/gather design delivers exactly what it was built for. (The CPU path
was also reproducible in every configuration tested, so the OpenMP reduction is
deterministic in practice here, though nothing in the standard guarantees it.)

**CPU vs GPU — close, but not bit-for-bit.**

| fixture | values | bit-identical | median rel. err | max rel. err |
| --- | --- | --- | --- | --- |
| valence 4-8 canonical patches | 210-290 | 56-70% | ~4e-16 | 3.1e-15 |
| flat mesh, 8,400 faces | 119,174 | 85.9% | 8.2e-16 | **6.2e-13** |

There are exactly two causes, both confirmed by ablation:

1. **FMA contraction.** nvcc defaults to `-fmad=true` and fuses `a*b+c` into a
   single rounding. The host build targets baseline x86-64, where the FMA3
   instruction is unavailable, so GCC emits a multiply and an add with two
   roundings. Same source, different arithmetic.
2. **OpenMP reduction order.** The CPU sums per-thread force buffers in thread
   order; the GPU gathers one-ring slots in fixed index order.

| configuration | patches | flat mesh |
| --- | --- | --- |
| as shipped | 56-70% | 85.9% |
| `-fmad=false` only | **100%** | 93.4% |
| `OMP_NUM_THREADS=1` only | — | 86.0% |
| **both** | **100%** | **100.00%** |

Neither alone is sufficient; together they are exact. Bit-exactness is
therefore *reachable* — but it costs about 14% on the device (`-fmad=false`:
2.40 s to 2.74 s) and requires a single-threaded CPU, which costs far more
(3.02 s to 8.41 s). **It is not recommended.** The FMA result is generally the
*more* accurate one, and matching the CPU would mean deliberately degrading the
GPU.

Two caveats on reading those numbers. About 2,500 of the differing values are
near-zero quantities (|v| < 1e-12) that flip sign; they are zero on both paths,
and ULP is meaningless across zero, so they are counted separately from the
error statistics. And the reason the production mesh (6.2e-13) is so much
looser than an isolated patch (3.1e-15) is **cancellation**: near equilibrium a
vertex force is a sum of large opposing contributions, which amplifies ~1e-16
rounding differences by roughly three orders of magnitude.

None of this is a correctness problem for the science: 500-iteration
trajectories at 99k faces agree to printed precision.

### 8.2 The mid-evaluation host reduction

Summing area and volume on the host costs a device-to-host transfer of
`2 * nFaces` doubles in the middle of every evaluation. It is kept
deliberately, because a device reduction's summation order depends on the block
count — so the totals would change with launch configuration, and the GPU would
no longer reproduce the CPU even to the ~1e-13 the test asserts.

Note the honest consequence: the *totals* agree because they are summed
identically, but as [section 8.1](#81-reproducibility-and-bitwise-agreement)
shows, the per-face inputs to that sum already differ in their last bits. The
host-side reduction removes one source of divergence, not all of them.

### 8.3 `gpu` fails loudly

`forceBackend = gpu` throws if no device is usable rather than falling back. A
job submitted as a GPU run that quietly produced CPU timings would be worse
than a job that failed. `auto` is the option that falls back, and it says so on
stdout. This is intentional and should not be "fixed".

### 8.4 Memory

Device memory is dominated by the per-slot force scratch: three buffers of
`nSlots * 3` doubles, with `nSlots ~= 12 * nFaces`, plus one corner buffer of
`nFaces * 9`. At 10^5 faces that scratch is about 94 MB, and the whole device
layout — topology, shape-function rows, coordinates and outputs — comes to
roughly 120 MB. Comfortable on a 6 GB card, but it scales linearly and is the
first thing to check on a very large mesh.
`CudaForceBackend::device_memory_bytes()` reports the exact total.

### 8.5 Build and toolchain

- **`SLIMED_CUDA_ARCHITECTURES` must match your device.** The historical
  default `70;80;90` **cannot build on CUDA 13.x**, which dropped Volta:
  `compute_70` is no longer a valid target. Pass your own compute capability
  (`89` for Ada, `90` for Hopper, `80` for A100).
- Device code is embedded per translation unit
  (`CUDA_SEPARABLE_COMPILATION OFF`); only `Cuda_force_backend.cu` is compiled
  by nvcc, and everything it includes is header-only.
- A build without `SLIMED_ENABLE_CUDA` links the stub and is fully functional
  on the CPU.

## 9. Verifying a new machine

**Run this before trusting a GPU build:**

```bash
ctest --test-dir build -R CudaForceBackend --output-on-failure
```

`CudaForceBackendTest.MatchesTheProductionForceEvaluation` runs the GPU and the
production CPU path over the same meshes — a grid, and the canonical patch for
every valence 4 through 8 — and compares every element area and volume, every
face energy, curvature and normal, and every vertex force.

**It skips rather than fails when no device is usable, and ctest reports a skip
as `Passed`.** A green line proves nothing on its own. Confirm it actually ran:

```bash
./build/bin/test_main --gtest_filter='CudaForceBackendTest.*'
```

and look for `[ OK ]` rather than `[ SKIPPED ]`. The reliable check is
differential — hiding the device must change the outcome:

```bash
CUDA_VISIBLE_DEVICES=-1 ./build/bin/test_main --gtest_filter='CudaForceBackendTest.*'
# expect: [ SKIPPED ] ... no usable CUDA device
```

**Known gap:** the test asserts `1e-13` relative, but its fixtures are small,
where accumulation is shallow. A production-scale mesh disagrees by up to
6.2e-13 ([section 8.1](#81-reproducibility-and-bitwise-agreement)) and would
*fail* that threshold. The tolerance should be relaxed to ~1e-11 and a large
fixture added.

## 10. Next steps

The optimisation that changes the performance picture is **device residency**:
keep coordinates on the device across the whole line search, move
`Model::update_vertex_using_NCG()` into a kernel, and reduce the totals on the
device. A line-search trial would then cost one kernel sequence and no
transfer, removing all three limits in
[section 6.1](#61-and-why-the-speedup-is-bounded) at once.

The reason it has not been done is
[section 8.2](#82-the-mid-evaluation-host-reduction): a device reduction
changes the totals with block count. Given
[section 8.1](#81-reproducibility-and-bitwise-agreement) has now established
that the two backends are *already* not bit-for-bit, that objection is weaker
than it first appeared — the honest framing is a documented ~1e-12 agreement,
not an exactness that was never actually there. Worth revisiting deliberately.

Smaller, independent wins:

- Overlap the four launches, or drop the per-launch `cudaDeviceSynchronize()`
  in favour of one sync per evaluation plus a single error check.
- Use pinned host staging buffers so the per-call transfers go at full PCIe
  bandwidth.
- Relax the equivalence-test tolerance and add a production-scale fixture
  ([section 9](#9-verifying-a-new-machine)).
