# Run time against the irregular-patch recursion depth

The material behind `docs/irregular_patch_results.md` §10 (measured 2026-10-24
on JHU Rockfish, Slurm job 32455115).

- `depth_sweep.sbatch`, `run_stream.sh`, `streams.txt`: the job. One exclusive
  48-core node; the 24 single-thread runs (12 variants × {1000, 2000} steps)
  are dealt over four streams per socket, each run pinned to one core with
  `taskset`, and the whole series is run twice (rep 1 on socket 0, rep 2 on
  socket 1). Then the 4-thread wall-time rows, two at a time on four cores each.
- `collect.sh`: turns the per-run `time.log` and `inputEdgeFlips.csv` into
  `results/results_<rep>.txt` (one line per run: `user` time, accepted flips).
- `depth_cost_replay.py`: rebuilds the final connectivity of a run by replaying
  its accepted flips on `inputmesh_faces.csv`, and counts the force pass's
  patch-kernel calls per step for a given depth set (1 per regular face,
  `3·D(N)` per extraordinary corner, plus the four Loop children of a
  multi-extraordinary face).
- `depth_cost_fit.py`: run inside a `rep*/` directory with the results file as
  argument; prints the per-step cost of steps 1000–2000, the linear fit against
  depth, the fit against kernel calls, and the depth-controlled share against
  the depth→0 intercept and against the frozen reference.
- `results/`: the raw result tables and the job's header (node, CPU, binary).

The run directories themselves are copies of
`~/slimed_runs/sticky_bead/matrix/box_nu5/input.params` with `maxIterations`
and, per variant, `irregularPatchDepths` / `irregularPatchTailCorrection`
(frozen: `edgeFlipEnabled = false`) appended. Reproduce locally with the same
binary target (`cmake --build build-omp --target membrane_dynamics`) by
running each directory in the foreground with `OMP_NUM_THREADS=1`.
