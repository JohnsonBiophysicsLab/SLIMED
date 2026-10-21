#!/usr/bin/env python
"""
How far does the time step go once the walls are substepped (Param::meshForceSubsteps)?
The flat periodic sheet of docs/resolution_limits.md, with the tile scaled with the
mesh -- 20 faces across at every lFace -- so that a run at lFace L and step dt is the
same lattice problem as one at lFace 5 and dt (5 / L)^2, as far as bending goes. Only
D dt / lFace^2 enters the bending update; the walls and the area term do not scale,
and that is what the sets below look for.

    python analysis/substep_resolution_runs.py --out ~/slimed_runs/substep_res --set SET

writes one directory per run (input.params, run.sh) and jobs_SET.txt, one run
directory per line, for `xargs -P 8 -I{} bash {}/run.sh < jobs_SET.txt` (local) or
any array scheduler. Every run uses split substep noise with n = ceil(dt / dt_sub(lFace)).

  prescan  the plain step (n = 1) near the wall bound at lFace 10 and 20, 8 seeds x 40000 steps
  stab     the step ladder up to the bending bound with substeps, 8 seeds
  spec     400 us x (L/5)^2, a frame every 0.2 us x (L/5)^2, 2 seeds: S(q), the slow Gamma(q)
  dyn      40 us x (L/5)^2, a frame every 0.02 us x (L/5)^2, 2 seeds: the fast Gamma(q)
  fluid    lFace 5 with flips (nu 0.5): stab, spec and dyn at a few steps
  confirm  lFace 20 at a substep of 0.005 us instead of 0.006, 16 seeds (0.006 diverged
           ~1 in 10^5 steps at dt 0.08-0.32)

The ladder is given as the lFace-5 step f; the run's step is f (L/5)^2. The elastic box
keeps the membrane's clock (mobility 0.02 (5/L)^2) and is capped at rate x dt = 0.1.
Reduce with analysis/resolution_sweep_analysis.py reduce (spec) or
analysis/substep_resolution_report.py reduce (dyn: the short window), and compare with
analysis/substep_resolution_report.py report. See docs/resolution_limits.md section 2.6.
"""
import argparse
import math
import os
import re

from sticky_bead_runs import SEED, write_run

FACES = 20            # faces across the periodic tile, at every lFace
LFACES = [5.0, 10.0, 20.0]
# the largest substep (us) the walls take at each lFace: the plain step's bound on this tile
# (prescan: lFace 10 0.005 had 1 divergence in 1.6 ms, lFace 20 0.006 none in 1.9 ms). With
# substeps 0.006 at lFace 20 diverged ~1 in 1e5 steps; use 0.005 there (the confirm set).
DT_SUB = {5.0: 0.0025, 10.0: 0.004, 20.0: 0.006}
STAB_F = [0.005, 0.01, 0.02, 0.03]
SPEC_F = [0.001, 0.0025, 0.005, 0.01, 0.02]
PRESCAN = [(10.0, [0.005, 0.006, 0.007]), (20.0, [0.005, 0.006, 0.008, 0.01])]
FLUID_STAB_F = [0.01, 0.02]
FLUID_SPEC_F = [0.005, 0.01]
FLUID_DYN_F = [0.001, 0.005, 0.01]
BOX_MOBILITY = 0.02   # nm / (pN us) at lFace 5
BOX_RATE_PER_MOBILITY = 2 * 250.0   # isotropic: 2 axes x usMembraneStretching


def sheet(l_face):
    """sideX, sideY of a sheet whose periodic tile is FACES faces across."""
    d_y = math.sqrt(3.0) / 2.0 * l_face
    rows = 2 * round(FACES / (math.sqrt(3.0) / 2.0) / 2)
    return (FACES + 6) * l_face, (rows + 6) * d_y


def set_keys(directory, values):
    """Replace `key = value` lines of input.params, appending the ones it does not have."""
    path = os.path.join(directory, "input.params")
    text = open(path).read()
    for key, value in values.items():
        text, n = re.subn(rf"^{key}\s*=.*$", f"{key} = {value}", text, flags=re.M)
        if n == 0:
            text += f"{key} = {value}\n"
    open(path, "w").write(text)


def make(root, binary, name, l_face, f, seed, time_lface5, frame_lface5, rate=0.0, substeps=None):
    scale = (l_face / 5.0) ** 2
    dt = round(f * scale, 9)
    steps = round(time_lface5 * scale / dt)
    interval = max(1, round(frame_lface5 * scale / dt))
    side_x, side_y = sheet(l_face)
    directory = write_run(root, name, "box", rate, steps, 0, interval, side=round(side_x, 6), potential="none",
                          time_step=dt, l_face=l_face, side_y=round(side_y, 6), seed=seed)
    n = substeps if substeps is not None else max(1, math.ceil(dt / DT_SUB[l_face] - 1e-9))
    mobility = min(BOX_MOBILITY / scale, 0.1 / (BOX_RATE_PER_MOBILITY * dt))
    set_keys(directory, {"meshForceSubsteps": n, "meshForceSubstepNoise": "split",
                         "elasticBoxMobility": f"{mobility:.6g}"})
    with open(os.path.join(directory, "run.sh"), "w") as handle:
        handle.write(f'#!/bin/bash\ncd "$(dirname "$0")"\nBIN=${{SLIMED_BIN:-{binary}}}\n'
                     f'/usr/bin/time -p -o time.log "$BIN" > run.log 2>&1\necho exit $? >> run.log\n')
    return directory


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--set", required=True, choices=("prescan", "stab", "spec", "dyn", "fluid", "confirm"))
    parser.add_argument("--bin", default=os.path.expanduser("~/Workspace/SLIMED/build/bin/membrane_dynamics"))
    args = parser.parse_args()
    root = os.path.expanduser(args.out)
    os.makedirs(root, exist_ok=True)
    jobs = []

    def add(*a, **k):
        jobs.append(make(root, args.bin, *a, **k))

    if args.set == "prescan":
        for l_face, dts in PRESCAN:
            for dt in dts:
                for k in range(8):
                    # 40000 steps each: the walls fail per step, not per unit of bending time
                    f = dt / (l_face / 5.0) ** 2
                    add(f"pre_l{l_face:g}_dt{dt:g}_s{k}", l_face, f, SEED + k, 40000 * f, 0.2, substeps=1)
    elif args.set == "stab":
        for l_face in LFACES:
            for f in STAB_F:
                for k in range(8):
                    add(f"stab_l{l_face:g}_f{f:g}_s{k}", l_face, f, SEED + k, 100.0, 0.2)
    elif args.set in ("spec", "dyn"):
        long_run = args.set == "spec"
        for l_face in LFACES:
            for f in SPEC_F:
                for k in range(2):
                    add(f"{args.set}_l{l_face:g}_f{f:g}_s{k}", l_face, f, SEED + k,
                        400.0 if long_run else 40.0, 0.2 if long_run else 0.02)
    elif args.set == "confirm":
        for f in (0.005, 0.01, 0.02):
            dt = f * 16.0
            for k in range(16):
                add(f"stab_l20_f{f:g}_sub0.005_s{k}", 20.0, f, SEED + k, 100.0, 0.2,
                    substeps=math.ceil(dt / 0.005 - 1e-9))
    else:
        for f in FLUID_STAB_F:
            for k in range(8):
                add(f"fstab_l5_f{f:g}_s{k}", 5.0, f, SEED + k, 40.0, 0.2, rate=0.5)
        for f in FLUID_SPEC_F:
            for k in range(2):
                add(f"fspec_l5_f{f:g}_s{k}", 5.0, f, SEED + k, 400.0, 0.2, rate=0.5)
        for f in FLUID_DYN_F:
            for k in range(2):
                add(f"fdyn_l5_f{f:g}_s{k}", 5.0, f, SEED + k, 40.0, 0.02, rate=0.5)

    with open(os.path.join(root, f"jobs_{args.set}.txt"), "w") as handle:
        handle.write("\n".join(jobs) + "\n")
    print(f"{len(jobs)} runs -> {root}/jobs_{args.set}.txt")


if __name__ == "__main__":
    main()
