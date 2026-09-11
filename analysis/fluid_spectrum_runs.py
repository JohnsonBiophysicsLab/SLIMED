"""The dynamics runs behind `membrane_fluctuation_fluid_cpu.ipynb`.

Kept in a module rather than a notebook cell so that a run launched by hand --
from a shell, on a cluster, ahead of the notebook -- writes exactly the
parameter file the notebook expects, and is then found "finished" rather than
regenerated.  `membrane_fluctuation_resample_cpu.ipynb` keeps the same logic
inline; the two share `COMMON` deliberately so the fluid runs are the solid
runs with the fluid flags added and nothing else changed.
"""

from __future__ import annotations

import os
import re
import subprocess
import time
from pathlib import Path

KBT = 4.17  # pN.nm, as in input.params

# A 100 nm box on a 5 nm mesh: nFace 20 x 24, so a 14 x 18 periodic tile.
COMMON = {
    "boundaryType": "Periodic",
    "isFlat": "true",
    "lFace": 5.0,
    "sideX": 100.0,
    "sideY": 100.0,
    "kcMembraneBending": 83.4,      # 20 kT
    "c0Membrane": 0.0,
    "uvVolumeConstraint": 0.0,
    "isGlobalConstraint": "true",
    "setRelaxAreaToDefault": "true",
    "KBT": KBT,
    "timeStep": 0.002,              # us
    "diffConst": 1.0,               # nm^2/us
    "forceBackend": "cpu",          # the CPU loops, not the CUDA backend
    "fdtConsistentSurfaceUpdate": "true",
    "integratePeriodicDuplicates": "false",
    "thermalFluctuationEnabled": "false",
    "isInsertionIncluded": "false",
    "isEnergyHarmonicBondIncluded": "false",
    "isGagScaffoldingEnergyIncluded": "false",
    "isIdealizedProteinLatticeEnergyIncluded": "false",
    "meshpointOutput": "true",
    "xyzOutput": "false",
    "VERBOSE_MODE": "false",
}

# Everything a fluid run adds to the solid one. The tether range and stiffness
# are the shipped defaults, written out so the file says what was run.
FLUID = {
    "surfaceSolver": "iterative",       # the dense inverse cannot survive a flip
    "inPlaneDynamicsEnabled": "true",   # half of what fluidity means
    "edgeSpringEnabled": "true",        # the mesh-quality term ...
    "edgeTetherShape": "flat",          # ... as the flat-bottomed tether
    "edgeTetherMinRatio": 0.95,
    "edgeTetherMaxRatio": 1.75,
    "edgeSpringConstant": 83.4,
    "edgeFlipEnabled": "true",
    "edgeFlipAttemptRate": 0.5,         # nu, attempts per edge per us
    "edgeFlipInterval": 1,
    # WP7: valence 8 cannot flatten inside the tether's lower wall and folds;
    # 5-7 gave zero creases where 4-8 gave a dozen.
    "edgeFlipMinValence": 5,
    "edgeFlipMaxValence": 7,
    "irregularPatchDepthScale": 1.0,
    # WP7: the tether bounds edge lengths, this bounds triangle shape. Every
    # long fluid run without it folded a sliver after 1e4-1e5 steps.
    "triangleShapeEnabled": "true",
    "triangleShapeMinAltitudeRatio": 0.4,
    "triangleShapeConstant": 83.4,
}

SPECS = {
    # the solid controls, exactly as membrane_fluctuation_resample_cpu.ipynb made them
    "pure_a":  {"usMembraneStretching": 0.0,   "maxIterations": 800_000,
                "meshpointOutputInterval": 100, "randomSeed": 101},
    "tension": {"usMembraneStretching": 250.0, "maxIterations": 800_000,
                "meshpointOutputInterval": 100, "randomSeed": 202},
    # the fluid membrane: tether, altitude floor, flips leaving valences 5-7
    "fluid_a": dict({"usMembraneStretching": 0.0, "maxIterations": 400_000,
                     "meshpointOutputInterval": 100, "randomSeed": 404,
                     "timeStep": 0.001}, **FLUID),
    # the same with the crease wall, which forbids the flaps the valence
    # restriction only made rare
    "fluid_b": dict({"usMembraneStretching": 0.0, "maxIterations": 400_000,
                     "meshpointOutputInterval": 100, "randomSeed": 404,
                     "timeStep": 0.001}, **FLUID,
                    creaseWallEnabled="true", creaseWallAngle=60.0, creaseWallConstant=500.0),
    # fluid_b at conserved area. With muS = 0 nothing sets the control net's
    # length scale but the tether walls: the mean edge grows to 6.15 nm, the
    # excess area buckles the sheet out of plane (r.m.s. height 9 nm against
    # the solid's 1), and the spectrum is that of a crumpled membrane. The
    # comparison with the solid `tension` run is the one at matched area.
    "fluid_c": dict({"usMembraneStretching": 250.0, "maxIterations": 400_000,
                     "meshpointOutputInterval": 100, "randomSeed": 505,
                     "timeStep": 0.001}, **FLUID,
                    creaseWallEnabled="true", creaseWallAngle=60.0, creaseWallConstant=500.0),
}
SCALABLE = ("pure_a", "tension", "fluid_a", "fluid_b", "fluid_c")
MARKER = ".made-by-this-notebook"


def params_text(name: str, spec: dict, scale: float = 1.0) -> str:
    p = dict(COMMON, **spec)
    if name in SCALABLE:
        p["maxIterations"] = max(2000, int(p["maxIterations"] * scale))
    return "".join(f"{k} = {v}\n" for k, v in p.items())


def state(d: Path, text: str) -> str:
    """"missing", "finished", "stale", or "foreign" -- see the resample notebook."""
    if not ((d / "inputvertex_final.csv").is_file() and (d / "meshpointinput.csv").is_file()):
        return "missing"
    pf = d / "input.params"
    if pf.is_file() and pf.read_text() == text:
        return "finished"
    return "stale" if (d / MARKER).is_file() else "foreign"


def last_iteration(log: Path) -> int:
    try:
        with open(log, "rb") as f:
            f.seek(0, os.SEEK_END)
            f.seek(max(0, f.tell() - 8192))
            hits = re.findall(rb"ITERATION:(\d+)", f.read())
        return int(hits[-1]) + 1 if hits else 0
    except OSError:
        return 0


def generate(runs: Path, exe: Path, specs: dict = SPECS, scale: float = 1.0,
             threads: int | None = None, launch: bool = True, poll: float = 15.0) -> dict:
    """Make every run in `specs` under `runs`, keeping any already finished.

    Returns ``{name: state}`` as found *before* launching.  With
    ``launch=False`` nothing is started, which is how a notebook analyses a run
    still in progress: it reports what it found and goes on with the frames
    that exist.
    """
    found = {}
    jobs = {}
    n_cpu = os.cpu_count() or 4
    for name, spec in specs.items():
        d = runs / name
        d.mkdir(parents=True, exist_ok=True)
        text = params_text(name, spec, scale)
        st = found[name] = state(d, text)
        if st in ("finished", "foreign"):
            note = "" if st == "finished" else (
                "  (not made here, and its parameters differ -- delete the directory to regenerate)")
            print(f"{name:12s} finished already, keeping it{note}")
            continue
        if not launch:
            print(f"{name:12s} {st}: not launched (launch=False); analysing what exists")
            continue
        if st == "stale":
            print(f"{name:12s} parameters changed, running it again")
        (d / "input.params").write_text(text)
        (d / MARKER).touch()
        for stale in ("inputvertex_final.csv", "meshpointinput.csv", "surfacepointinput.csv"):
            (d / stale).unlink(missing_ok=True)
        for old in d.glob("inputface_*.csv"):
            old.unlink()
        # A fluid run threads well (its time is in the per-face patch kernels)
        # and a solid one barely does; give whichever runs the cores.
        n_threads = threads or max(1, min(8, n_cpu // max(1, sum(
            1 for n, s in specs.items() if state(runs / n, params_text(n, s, scale)) not in ("finished", "foreign")))))
        log = open(d / "run.log", "w")
        proc = subprocess.Popen([str(exe)], cwd=d, stdout=log, stderr=subprocess.STDOUT,
                                env=dict(os.environ, OMP_NUM_THREADS=str(n_threads)))
        jobs[name] = (proc, log, int(re.search(r"maxIterations = (\d+)", text).group(1)))
        print(f"{name:12s} started, {jobs[name][2]:,} steps on {n_threads} threads")

    t0 = time.time()
    while any(p.poll() is None for p, _, _ in jobs.values()):
        time.sleep(poll)
        done = " ".join(f"{n} {last_iteration(runs / n / 'run.log') / tot:5.1%}"
                        for n, (_, _, tot) in jobs.items())
        print(f"\r  {time.time() - t0:5.0f} s   {done}", end="", flush=True)
    for name, (proc, log, _) in jobs.items():
        log.close()
        if proc.returncode:
            raise RuntimeError(f"{name} exited {proc.returncode}; see {runs / name / 'run.log'}")
        if "diverged" in (runs / name / "run.log").read_text()[-4000:]:
            raise RuntimeError(f"{name} diverged; see the tail of {runs / name / 'run.log'}")
    if jobs:
        print(f"\r  {time.time() - t0:5.0f} s   all runs finished" + " " * 20)
    return found


if __name__ == "__main__":
    import sys
    # `python fluid_spectrum_runs.py fluid_a` prints the parameter file a run would get.
    name = sys.argv[1] if len(sys.argv) > 1 else "fluid_a"
    print(params_text(name, SPECS[name]), end="")
