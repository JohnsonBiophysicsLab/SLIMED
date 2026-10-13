#!/usr/bin/env python
"""
The sticky-bead matrix: a rigid bead bonded to a periodic membrane by the full
harmonic bond (beadPotential = harmonic), pushed 35 nm through it and held,
under two area reservoirs (the elastic box, and the adaptive relaxed area of
the lipid-flow modes) and several fluidities (edge-flip attempt rates).
--potential wall is the same matrix with a bead that only pushes (the
repulsive half of the bond, no adhesion), the control for what stickiness does.

    python analysis/sticky_bead_runs.py --out ~/slimed_runs/sticky_bead --proto
    python analysis/sticky_bead_runs.py --out ~/slimed_runs/sticky_bead --matrix
    python analysis/sticky_bead_runs.py --out ~/slimed_runs/sticky_bead/push_matrix --matrix --potential wall

writes one directory per run with its input.params and a run.sh, and, for the
matrix, a Slurm array script for the cluster. See docs/bead_push_benchmark.md
section 9 for what the runs showed.
"""
import argparse
import os
import stat
import textwrap

PUSH_ITERATIONS = 70000      # 35 nm at 0.5 nm/us
HOLD_ITERATIONS = 110000     # at least 100 us after the bead stops
SEED = 20261012

BASE = """\
############################################
## The sticky bead: docs/bead_push_benchmark.md section 9
############################################

beadEnabled = true
{bead}
beadCenterX = 0.0
beadCenterY = 0.0
beadCenterZ = -5.0               # the sheet is at z = 10: 5 nm of clearance at the start
beadVelocityX = 0.0
beadVelocityY = 0.0
beadVelocityZ = 0.5              # nm/us
beadStartIteration = 0
beadStopIteration = {push}       # 35 nm of travel, then it holds at z = 30

############################################
## Boundary Conditions: periodic along both axes, mean height pinned
############################################

boundaryType = Mixed
boundaryTypeX = Periodic
boundaryTypeY = Periodic
meanHeightPinned = true          # the surrounding membrane holds the tile's mean height
surfaceSolver = iterative

############################################
## Geometric Parameters
############################################

lFace = 5.0
isFlat = true
sideX = {side}
sideY = {side}                    # the generated sheet; the periodic tile is 30 nm narrower (the image bands)

############################################
## Physical Properties
############################################

c0Insertion = 0.0
c0Membrane = 0.0
kcMembraneBending = 83.4         # 20 kT
usMembraneStretching = 250.0     # the bilayer's modulus; the reservoir supplies the area
setRelaxAreaToDefault = true
uvVolumeConstraint = 0.0
isGlobalConstraint = true
forceBackend = cpu
KBT = 4.17

############################################
## Dynamics
############################################

maxIterations = {total}
timeStep = 0.001                 # us
diffConst = 1.0                  # nm^2/us
randomSeed = {seed}
meshpointOutput = true
meshpointOutputInterval = {interval}
xyzOutput = false
fdtConsistentSurfaceUpdate = true
localOutput = true               # per-face energies and per-vertex forces beside every frame

############################################
## The fluid mesh: docs/edge_flip_plan.md, docs/fluidity_results.md
############################################

inPlaneDynamicsEnabled = true
edgeSpringEnabled = true
edgeTetherShape = flat
edgeTetherMinRatio = 0.7
edgeTetherMaxRatio = 1.75
edgeSpringConstant = 83.4
triangleShapeEnabled = true
triangleShapeMinAltitudeRatio = 0.4
triangleShapeConstant = 83.4
creaseWallEnabled = true
creaseWallAngle = 60.0
creaseWallConstant = 500.0
edgeFlipEnabled = {flips}
edgeFlipAttemptRate = {rate}     # attempts per edge per us: the fluidity knob
edgeFlipInterval = 1
edgeFlipMinValence = 4
edgeFlipMaxValence = 8
irregularPatchDepthScale = 1.0

############################################
## The area reservoir: {reservoir_label}
############################################

{reservoir}

############################################
## Off for this run
############################################

isInsertionIncluded = false
isEnergyHarmonicBondIncluded = false
thermalFluctuationEnabled = false
externalForceEnabled = false
VERBOSE_MODE = false
"""

BEADS = {
    "harmonic": """\
beadPotential = harmonic         # the full harmonic bond on every vertex the bead captured
beadRadius = 10.0                # nm, two mesh edges
beadStiffness = {stiffness}      # pN/nm, the bond's stiffness (and the wall's): k/2 (2.5 nm)^2 = 15 kT at capture
beadCaptureRange = 2.5           # nm outside the surface: a vertex this close is bonded
beadBreakRange = -1              # never breaks
beadAdhesionEnergy = 0.0         # ignored under harmonic""",
    "wall": """\
beadPotential = wall             # the repulsive half of the bond only: the bead pushes and holds nothing
beadRadius = 10.0                # nm, two mesh edges
beadStiffness = {stiffness}      # pN/nm, the wall's stiffness
beadAdhesionEnergy = 0.0         # no well: a bead that only pushes""",
}

RESERVOIRS = {
    "fixed": ("no reservoir: a fixed box and the global area constraint",
              "elasticBoxEnabled = false\nlipidFlowMode = off"),
    "box": ("uniform scaling of the box at a frame tension of 0.5 pN/nm (docs/elastic_box.md)",
            textwrap.dedent("""\
            elasticBoxEnabled = true
            elasticBoxTension = 0.5
            elasticBoxMobility = 0.02        # fast enough to follow a bead at 0.5 nm/us
            elasticBoxThermal = true
            elasticBoxCoupling = isotropic
            elasticBoxMaxStrainPerStep = 0.001
            elasticBoxVertexGasTerm = true
            lipidFlowMode = off""")),
    "s0global": ("adaptive relaxed area, one A0 for the sheet, bulk reservoir at 0.5 pN/nm (docs/lipid_flow_plan.md)",
                 textwrap.dedent("""\
                 elasticBoxEnabled = false
                 lipidFlowMode = global
                 lipidFlowThermal = true
                 lipidReservoirRegion = bulk
                 lipidReservoirTension = 0.5
                 lipidReservoirRate = 20""")),
    "s0local": ("adaptive relaxed area per face, tension-driven flow, bulk reservoir at 0.5 pN/nm",
                textwrap.dedent("""\
                elasticBoxEnabled = false
                lipidFlowMode = local
                lipidFlowThermal = true
                lipidFlowLinearStretching = true
                lipidFlowDiffusivity = -1
                lipidFlowSubsteps = 0
                lipidFlipResampling = auto
                lipidReservoirRegion = bulk
                lipidReservoirTension = 0.5
                lipidReservoirRate = 20""")),
}

# (name, attempt rate); 0 = flips off, the frozen mesh
FLUIDITIES = [("frozen", 0.0), ("nu0.05", 0.05), ("nu0.5", 0.5), ("nu5", 5.0)]


def write_run(root, name, reservoir, rate, total, push, interval, side=130.0, stiffness=20.0, seed=SEED, threads=4, binary=None,
              potential="harmonic"):
    label, block = RESERVOIRS[reservoir]
    text = BASE.format(push=push, total=total, seed=seed, interval=interval, side=side,
                       bead=BEADS[potential].format(stiffness=stiffness),
                       flips="true" if rate > 0 else "false", rate=rate if rate > 0 else 0.5,
                       reservoir_label=label, reservoir=block)
    directory = os.path.join(root, name)
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, "input.params"), "w") as handle:
        handle.write(text)
    binary = binary or os.path.expanduser("~/Workspace/SLIMED/build-omp/bin/membrane_dynamics")
    script = os.path.join(directory, "run.sh")
    with open(script, "w") as handle:
        handle.write(f"#!/bin/bash\ncd \"$(dirname \"$0\")\"\n"
                     f"export OMP_NUM_THREADS=${{OMP_NUM_THREADS:-{threads}}}\n"
                     f"BIN=${{SLIMED_BIN:-{binary}}}\n"
                     f"\"$BIN\" > run.log 2>&1\necho \"{name} exit $?\" >> ../status.log\n")
    os.chmod(script, os.stat(script).st_mode | stat.S_IEXEC)
    return directory


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--proto", action="store_true", help="three short runs of --proto-steps for a check")
    parser.add_argument("--proto-steps", type=int, default=10000)
    parser.add_argument("--matrix", action="store_true", help="the whole matrix")
    parser.add_argument("--interval", type=int, default=200)
    parser.add_argument("--side", type=float, default=130.0, help="sideX = sideY of the generated sheet, nm")
    parser.add_argument("--stiffness", type=float, default=20.0, help="beadStiffness, pN/nm (200 folds a frozen mesh)")
    parser.add_argument("--potential", choices=sorted(BEADS), default="harmonic",
                        help="harmonic: the sticky bead; wall: the same bead that only pushes")
    args = parser.parse_args()
    extra = dict(potential=args.potential)
    root = os.path.expanduser(args.out)
    names = []
    if args.proto:
        for reservoir in ("box", "s0global", "s0local"):
            names.append(write_run(root, f"proto_{reservoir}_nu0.5", reservoir, 0.5,
                                   args.proto_steps, min(PUSH_ITERATIONS, args.proto_steps), args.interval, args.side, args.stiffness, **extra))
    if args.matrix:
        total = PUSH_ITERATIONS + HOLD_ITERATIONS
        names.append(write_run(root, "fixed_nu0.5", "fixed", 0.5, total, PUSH_ITERATIONS, args.interval, args.side, args.stiffness, **extra))
        for reservoir in ("box", "s0global", "s0local"):
            for fname, rate in FLUIDITIES:
                names.append(write_run(root, f"{reservoir}_{fname}", reservoir, rate, total,
                                       PUSH_ITERATIONS, args.interval, args.side, args.stiffness, **extra))
        with open(os.path.join(root, "runs.txt"), "w") as handle:
            handle.write("\n".join(os.path.basename(n) for n in names) + "\n")
        with open(os.path.join(root, "submit_array.sh"), "w") as handle:
            handle.write(textwrap.dedent(f"""\
            #!/bin/bash
            #SBATCH --job-name={"sticky" if args.potential == "harmonic" else "push"}_bead
            #SBATCH --partition=shared
            #SBATCH --account=mjohn218
            #SBATCH --nodes=1
            #SBATCH --ntasks=1
            #SBATCH --cpus-per-task=4
            #SBATCH --time=14:00:00
            #SBATCH --array=1-{len(names)}
            #SBATCH --output=slurm_%A_%a.out
            cd "$SLURM_SUBMIT_DIR"
            NAME=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" runs.txt)
            export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
            export SLIMED_BIN=${{SLIMED_BIN:-/home/yying7/scr4-mjohn218/yying7/slimed/SLIMED_sticky/build-omp/bin/membrane_dynamics}}
            bash "$NAME/run.sh"
            """))
    for n in names:
        print(n)


if __name__ == "__main__":
    main()
