#!/usr/bin/env python
"""
How coarse can a run be? The flat periodic sheet of docs/flat_tension_matrix.md
(no bead, the elastic box at a frame tension of 0 pN/nm), swept in its time
step and in its mesh size, frozen and fluid (edge flips at 0.5 per edge per us).

    python analysis/resolution_sweep_runs.py --out ~/slimed_runs/resolution

writes one directory per run (input.params, run.sh), runs.txt and a Slurm
array script, submit_array.sh. Three sets:

  dt_*      lFace 5, the 100 nm tile, timeStep from 0.0005 to 0.04 us
  mesh_*    timeStep 0.001 (0.0005 below lFace 5), lFace from 3.33 to 20 nm on
            a tile of about 100 nm (the sheet is the tile plus three image
            rings and a duplicate ring on each side)
  cross_*   lFace 10, frozen, timeStep 0.008 to 0.05: does the stable step grow
            with the mesh size the way free-draining bending says, as lFace^2?

--set stability (runs_stability.txt, bash submit_array.sh runs_stability.txt)
narrows the stable step and finds the term that sets it:

  stab_*    lFace 5, timeStep 0.0025 to 0.0035, frozen and fluid
  diag_*    lFace 5, frozen, timeStep 0.004 to 0.008 with one thing removed: the
            crease wall, the altitude floor, every mesh term (edgeSpringEnabled
            false), every mesh term and the in-plane motion, or the elastic box;
            or with the altitude floor (and the tether) at a quarter of the stiffness
  big_*     lFace 10 and 20, frozen, the step between the stable and diverged ones

--set replicas (runs_replicas.txt): divergence comes at a random time, so the
steps on either side of each boundary again with two more seeds (seed*), and
the height-only dynamics (zonly_*: no in-plane motion, no mesh terms) pushed
until bending itself goes unstable, at lFace 5 and 10.

Every run is the same physical length (--time, us) with a frame every 0.2 us.
The reference point of both sweeps, lFace 5 at dt 0.001, is flat_tension/box_*.
See docs/resolution_limits.md for what the runs showed.
"""
import argparse
import math
import os
import textwrap

from sticky_bead_runs import SEED, write_run

TILE = 100.0          # nm, the periodic tile along x (and about that along y)
FRAME = 0.2           # us between trajectory frames
FLUIDITIES = [("frozen", 0.0), ("nu0.5", 0.5)]
DTS = [0.0005, 0.002, 0.004, 0.005, 0.008, 0.01, 0.0125, 0.02, 0.025, 0.04]
# edges across the tile: lFace = 100 / n
MESH_EDGES = [30, 16, 12, 10, 8, 6, 5]
CROSS_DTS = [0.008, 0.016, 0.025, 0.05]
STAB_DTS = [0.0025, 0.003, 0.0035]
BIG = [(10.0, [0.004, 0.005, 0.006]), (20.0, [0.008, 0.012, 0.016, 0.025])]
# (lFace, flips, dt) near each boundary, run again with these seeds
REPLICA_SEEDS = [20261017, 20261018]
REPLICAS = [(5.0, 0.0, 0.0025), (5.0, 0.0, 0.003), (5.0, 0.5, 0.002), (5.0, 0.5, 0.0025), (5.0, 0.5, 0.003),
            (10.0, 0.0, 0.006), (10.0, 0.0, 0.007), (10.0, 0.5, 0.004), (10.0, 0.5, 0.006),
            (20.0, 0.0, 0.008), (20.0, 0.0, 0.01), (20.0, 0.5, 0.008)]
ZONLY = [(5.0, [0.012, 0.016, 0.02, 0.025, 0.03, 0.04, 0.05]), (10.0, [0.025, 0.04, 0.05, 0.08, 0.12, 0.16, 0.25])]


def sheet(l_face):
    """sideX, sideY of a sheet whose periodic tile is about TILE on a side."""
    n_x = round(TILE / l_face)
    d_y = math.sqrt(3.0) / 2.0 * l_face
    rows = 2 * round(TILE / d_y / 2)          # the tile's row count must be even
    return (n_x + 6) * l_face, (rows + 6) * d_y


DIAGNOSES = {
    "nocrease": {"creaseWallEnabled": "false"},
    "noshape": {"triangleShapeEnabled": "false"},
    "nomesh": {"edgeSpringEnabled": "false", "creaseWallEnabled": "false", "triangleShapeEnabled": "false"},
    "zonly": {"edgeSpringEnabled": "false", "creaseWallEnabled": "false", "triangleShapeEnabled": "false",
              "inPlaneDynamicsEnabled": "false"},
    "fixedbox": {"elasticBoxEnabled": "false"},
    # the walls at a quarter of their stiffness: does the stable step grow with them?
    "softshape": {"triangleShapeConstant": "20.85"},
    "softwalls": {"triangleShapeConstant": "20.85", "edgeSpringConstant": "20.85"},
}


def override(directory, values):
    """Replace `key = value` lines of a written input.params, keeping their comments."""
    path = os.path.join(directory, "input.params")
    lines = open(path).read().splitlines()
    for key, value in values.items():
        hits = [i for i, line in enumerate(lines) if line.split("=")[0].strip() == key]
        assert len(hits) == 1, (key, hits)
        lines[hits[0]] = f"{key} = {value}"
    open(path, "w").write("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--time", type=float, default=400.0, help="physical length of every run, us")
    parser.add_argument("--set", choices=("sweep", "stability", "replicas"), default="sweep")
    args = parser.parse_args()
    root = os.path.expanduser(args.out)
    runs = []   # (name, cpus, hours)

    def add(name, rate, dt, l_face, values=None, seed=SEED):
        side_x, side_y = sheet(l_face)
        steps = round(args.time / dt)
        # a frame every FRAME us, to the nearest step
        directory = write_run(root, name, "box", rate, steps, 0, round(FRAME / dt), side=round(side_x, 6),
                              potential="none", time_step=dt, l_face=round(l_face, 6), side_y=round(side_y, 6), seed=seed)
        if values:
            override(directory, values)
        # cost: ~385 steps/s frozen and ~31 fluid on the 480-vertex tile at 4 cores
        # (flat_tension, job 31990185), scaled with the vertex count
        vertices = (TILE / l_face) ** 2 / 20 ** 2
        hours = steps * vertices / (31.0 if rate > 0 else 385.0) / 3600
        cpus = 8 if hours > 6 else 4
        hours = hours / (1.6 if cpus == 8 else 1.0)
        runs.append((name, cpus, min(36, max(1, math.ceil(1.5 * hours + 0.5)))))

    if args.set == "sweep":
        for fname, rate in FLUIDITIES:
            for dt in DTS:
                add(f"dt{dt:g}_{fname}", rate, dt, 5.0)
            for n in MESH_EDGES:
                l_face = TILE / n
                add(f"mesh{l_face:.3g}_{fname}", rate, 0.0005 if l_face < 5 else 0.001, l_face)
        for dt in CROSS_DTS:
            add(f"cross_l10_dt{dt:g}_frozen", 0.0, dt, 10.0)
    elif args.set == "replicas":
        for l_face, rate, dt in REPLICAS:
            for k, seed in enumerate(REPLICA_SEEDS, start=2):
                add(f"seed{k}_l{l_face:g}_dt{dt:g}_{'nu0.5' if rate else 'frozen'}", rate, dt, l_face, seed=seed)
        for l_face, dts in ZONLY:
            for dt in dts:
                add(f"zonly_l{l_face:g}_dt{dt:g}", 0.0, dt, l_face, DIAGNOSES["zonly"])
    else:
        for fname, rate in FLUIDITIES:
            for dt in STAB_DTS:
                add(f"stab_dt{dt:g}_{fname}", rate, dt, 5.0)
        for diagnosis, values in DIAGNOSES.items():
            for dt in (0.004, 0.005, 0.008):
                add(f"diag_{diagnosis}_dt{dt:g}", 0.0, dt, 5.0, values)
        for l_face, dts in BIG:
            for dt in dts:
                add(f"big_l{l_face:g}_dt{dt:g}_frozen", 0.0, dt, l_face)

    runs_file = "runs.txt" if args.set == "sweep" else f"runs_{args.set}.txt"
    with open(os.path.join(root, runs_file), "w") as handle:
        handle.write("".join(f"{name} {cpus} {hours}\n" for name, cpus, hours in runs))
    with open(os.path.join(root, "submit_array.sh"), "w") as handle:
        handle.write(textwrap.dedent("""\
        #!/bin/bash
        # bash submit_array.sh [runs file]: one sbatch per run, with the cores and hours the file gives it
        cd "$(dirname "$0")"
        BIN=${SLIMED_BIN:-/home/yying7/scr4-mjohn218/yying7/slimed/SLIMED_sticky/build-omp/bin/membrane_dynamics}
        while read -r NAME CPUS HOURS; do
            sbatch --job-name=res_$NAME --partition=shared --account=mjohn218 --nodes=1 --ntasks=1 \\
                --cpus-per-task=$CPUS --time=$HOURS:00:00 --output=$NAME/slurm_%j.out \\
                --export=ALL,SLIMED_BIN=$BIN,OMP_NUM_THREADS=$CPUS --wrap="bash $NAME/run.sh"
        done < "${1:-runs.txt}"
        """))
    for name, cpus, hours in runs:
        print(f"{name:28s} {cpus} cores {hours:3d} h")


if __name__ == "__main__":
    main()
