# GPU acceleration

This document records what was measured, what was changed, and — importantly —
what has and has not been run on real hardware.

## Where the time was going

Profiled with `sample` at 1 ms over the default `input.params` (8400 faces,
4331 vertices), single-threaded, on an Apple M5.

| | share of wall time |
| --- | --- |
| `Mesh::Compute_Energy_And_Force()` | 95.3% |
| └ the face loop → `element_energy_force_patch()` | 91.8% |
| everything else | under 5% |

Broken down by what was actually executing inside that loop:

| | ~share of total |
| --- | --- |
| `malloc` / `free` | 34% |
| `gsl_matrix_get` / `gsl_matrix_set` and their cross-library PLT stubs | 30% |
| `bzero` / `memset`, from `gsl_matrix_calloc` | 5% |
| arithmetic (`cblas_dgemm`, `cblas_dcopy`, `kron`) | 6% |

**About 70% of the program's runtime was memory allocation and dynamic-linker
overhead on 3-element vectors.** The cause is `Matrix` in
`include/linalg/Linear_algebra.hpp`: every matrix is a `gsl_matrix*` heap
allocation, and every scalar access is a call into `libgsl`.
`element_energy_force_patch()` declared 43 of them per call and built and threw
away a dozen more per control point, because `kron()` returns its output
matrix by value as well as writing through its out-parameter.

For scale: one force evaluation on that mesh is about 10⁸ flops. It should take
roughly 25 ms on one core. It took about a hundred times that.

None of this could be moved to a GPU as written — `gsl_matrix*`,
`std::vector<Matrix>` and per-3-vector heap allocation are not device code. So
the port's first step was also the largest single win available on the CPU.

## What changed

### Step 1 — the kernels, on flat arrays

`include/energy_force/Patch_kernel.hpp` and its `.inl` hold the patch energy
and force kernel and the area/volume quadrature, written over fixed-size stack
arrays with no allocation, no STL and no library boundary to cross. The
arithmetic mirrors the original operation for operation, including its
association order and `const_division()`'s reciprocal-then-multiply.

`Compute_energy_and_force_on_mesh.cpp`, `Mesh::calculate_element_area_volume()`
and `Mesh::energy_force_regularization()` were moved onto them.

Measured on the default input, 20 iterations:

| | before | after |
| --- | --- | --- |
| 1 thread | 184 s | 17.8 s |
| 10 threads | 35.9 s | 6.9 s |

**11.3× serial, 5.2× at 10 threads**, with no change to the model.

After the rewrite the profile reads as it should: 93% in
`Compute_Energy_And_Force`, of which about 81% is the patch kernel doing real
floating-point work and 3.5% is `malloc`.

### Step 2 — the mesh, flattened

`include/cuda/Device_mesh_layout.hpp` repacks the mesh into contiguous arrays:
one-ring lists as CSR, per-face patch descriptors, and — the part that matters
— the maps that put forces back on vertices.

Force accumulation is the one thing a GPU cannot do the way the CPU does. The
CPU keeps a full-length buffer per OpenMP thread and sums them; with thousands
of threads that is not an option, and `atomicAdd` would make the answer depend
on thread scheduling. Instead each face writes only its own slots, into a
buffer laid out exactly like the one-ring lists, so no two faces ever touch the
same slot and no synchronisation is needed. A second pass has each vertex
gather the slots that belong to it, in a fixed order. **The result is
deterministic run to run, which the OpenMP path is not.**

### Step 3 — one set of kernel bodies, two drivers

`include/cuda/Force_kernels.hpp` expresses each stage as a function taking one
face or vertex index. `HostForceBackend` runs them in a loop;
`Cuda_force_backend.cu` wraps each in a `__global__` that computes the index
from `blockIdx` and `threadIdx`. There is no host copy and device copy of the
physics to keep in sync — there is one body, compiled twice.

The `.cu` is consequently thin: four three-line kernels, and otherwise
allocation, transfer and launch.

## Building with CUDA

```bash
cmake -B build -DSLIMED_ENABLE_CUDA=ON -DSLIMED_CUDA_ARCHITECTURES=80
cmake --build build -j
```

Set `SLIMED_CUDA_ARCHITECTURES` to your device's compute capability without the
dot — `70` for V100, `80` for A100, `90` for H100. The default builds all
three, which is slower to compile. On a cluster you will usually need
`module load cuda` first.

Then select it in `input.params`:

```
forceBackend = gpu      # or: cpu (default), auto
```

`gpu` fails the run if no device is usable, rather than falling back — a job
submitted as a GPU run should not silently report CPU timings. `auto` falls
back and says so.

## Verify before trusting it

**Run this first on any new machine:**

```bash
ctest --test-dir build -R CudaForceBackend
```

`CudaForceBackendTest.MatchesTheProductionForceEvaluation` runs the GPU and the
production CPU path over the same meshes — a grid, and the canonical patch for
every valence 4 through 8 — and compares every element area and volume, every
face energy, curvature and normal, and every vertex force, to 1e-13 relative.

It **skips** rather than passes when CUDA was not compiled in or no device is
visible. A green run on a machine without a GPU means nothing was checked; read
the skip message.

## What has and has not been run

Written and validated on an Apple M5, which has no NVIDIA GPU and no CUDA
toolkit. So, honestly:

**Verified here:**

- The kernel bodies, against the preserved `Matrix`-based originals
  (`tests/test_patch_kernel.cpp`), over regular and valence 4–8 patches.
- The whole flattened pipeline — layout, scatter, gather, child-slot
  arithmetic — against the production force evaluation
  (`tests/test_device_mesh_layout.cpp`).
- End-to-end: 30 iterations of the default input reproduce the pre-rewrite
  `ElementFaceEnergy.csv` bit for bit and every non-zero energy column exactly.
- `Cuda_force_backend.cu` itself: type-checked with `__CUDACC__` defined, then
  compiled and executed against a stand-in CUDA runtime with each launch
  expanded into a serial sweep over the same index range. It reproduced the
  production force evaluation to 1.1e-15 on every fixture. That exercises the
  buffer sizing, the transfers, the argument wiring and the launch bounds.

**Not verified here, because it cannot be:** nvcc actually compiling the file,
real kernel execution, occupancy, and anything about performance on a device.
No GPU speedup is claimed below, because none has been measured.

## What to expect, and the next step

The expectation is not "the GPU is faster". At the default 8400 faces this
backend will most likely **lose** to the 10-thread CPU path, for a structural
reason worth understanding before benchmarking.

The line search in `Model::linear_search_for_stepsize_to_minimize_energy()`
calls the force evaluation on the order of 250 times per iteration, changing
only vertex coordinates. This backend transfers coordinates in and results out
on every one of those calls, and synchronises mid-evaluation to sum the area
and volume totals on the host. At 4331 vertices those transfers are a few
hundred microseconds — comparable to, or larger than, the kernel time.

So: **expect this to pay off at large meshes** — 10⁵ faces and up, where kernel
time grows and the per-call transfer stays roughly proportional but the
arithmetic-to-transfer ratio improves — and expect it to lose on small ones.
Benchmark both before choosing.

The optimisation that changes this picture is device residency: keep vertex
coordinates on the device across the whole line search, move
`Model::update_vertex_using_NCG()` into a kernel, and reduce the area and
volume totals on the device. Then a line-search trial costs one kernel sequence
and no transfer at all. Two things stand in the way, both deliberate rather
than overlooked:

1. The host-side reduction is what makes the GPU and CPU agree bit for bit. A
   device reduction's summation order depends on the block count, so neither
   path would reproduce the other. That is a real cost to accept knowingly, not
   by accident.
2. It reaches up into the optimiser, not just the force evaluation, so it wants
   the equivalence tests above to be passing on real hardware first.

## Files

| | |
| --- | --- |
| `include/energy_force/Patch_kernel.{hpp,inl}` | the patch and area/volume kernels, host and device |
| `include/energy_force/Patch_rows_flat.hpp` | shape-function rows repacked contiguous |
| `include/cuda/Device_mesh_layout.hpp` | the mesh as flat arrays, and the gather maps |
| `include/cuda/Force_kernels.{hpp,inl}` | one body per stage, indexed by face or vertex |
| `include/cuda/Host_force_backend.hpp` | runs them in a loop; the fallback and the test oracle |
| `include/cuda/Cuda_force_backend.hpp` | runs them on a device |
| `src/cuda/Cuda_force_backend.cu` | kernels, allocation, transfer, launch |
| `src/cuda/Cuda_force_backend_stub.cpp` | stands in when CUDA is not compiled in |
