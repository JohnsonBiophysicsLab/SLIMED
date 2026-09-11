"""Fluidity measurements for a SLIMED run with Monte Carlo edge flips.

What this module is for
-----------------------
A triangulated membrane whose connectivity is fixed is a solid: two vertices
that start as neighbours stay neighbours forever, so the sheet carries an
in-plane shear modulus a lipid bilayer does not have.  The flip move removes
that shear modulus.  Whether it *has* removed it is not something the
acceptance rate can answer -- a chain can accept moves briskly and still be
sampling the wrong distribution, or accept them and leave the membrane
effectively solid because the moves undo each other.

Three measurements answer it, in increasing order of how hard they are to fake:

1. **Acceptance rate and valence histogram.**  Diagnostics.  An equilibrium
   dynamically triangulated surface sits near 60% valence 6 with about 20%
   each of 5 and 7 -- a sanity band, not a target -- and the accepted energy
   differences should straddle zero.  A chain whose accepted moves are
   systematically downhill is not at equilibrium, or is differencing a
   different Hamiltonian from the one the dynamics integrates.

2. **Neighbour survival.**  The fraction of the initial edges still present at
   time t.  This is the microscopic fluidity time scale, and it is what the
   attempt rate `nu` should be calibrated against: a lipid patch of the size
   one vertex stands for exchanges a neighbour on a known time scale, and the
   model should too.  With flips off it is exactly 1 forever, which makes it a
   clean control.

3. **In-plane mean squared displacement.**  The classical signature.  Without
   flips a vertex is tethered to a fixed set of neighbours, so its in-plane
   MSD saturates at the cage size.  With flips the cage itself rearranges and
   the MSD grows without bound -- linearly in t, in the diffusive regime.
   This is the one that cannot be faked by a move that shuffles without
   transporting.

Inputs are the files a dynamics run writes into its working directory:
`<prefix>EdgeFlips.csv`, `<prefix>face_<iteration>.csv`, `meshpoint<prefix>.csv`
and `<prefix>vertex_type_begin.csv`.  See docs/edge_flip_plan.md work package 6.
"""

from __future__ import annotations

import glob
import os
import re
from dataclasses import dataclass, field

import numpy as np


# ---------------------------------------------------------------------------
# Reading a run
# ---------------------------------------------------------------------------

@dataclass
class Run:
    """Everything the measurements below need from one run directory."""

    directory: str
    prefix: str = "input"
    #: (n_frames, n_vertices, 3) control-net coordinates, one frame per output.
    coords: np.ndarray = field(default=None, repr=False)
    #: Iteration number of each coordinate frame.
    frame_iterations: np.ndarray = field(default=None, repr=False)
    #: iteration -> (n_faces, 3) connectivity, for the iterations one was written.
    face_frames: dict = field(default_factory=dict, repr=False)
    #: (n_attempts, ) arrays from the flip log.
    flip_iteration: np.ndarray = field(default=None, repr=False)
    flip_delta_energy: np.ndarray = field(default=None, repr=False)
    flip_accepted: np.ndarray = field(default=None, repr=False)
    #: Vertices that are neither ghost nor a periodic duplicate.
    free_vertices: np.ndarray = field(default=None, repr=False)

    def path(self, name: str) -> str:
        return os.path.join(self.directory, name)


def _read_face_csv(path: str) -> np.ndarray:
    return np.loadtxt(path, delimiter=",", dtype=int, ndmin=2)


def load_run(directory: str, prefix: str = "input",
             output_interval: int = 1) -> Run:
    """Read a run directory.

    `output_interval` is `meshpointOutputInterval`; the trajectory files carry
    no iteration column, so the frame times have to come from the parameter
    that set them.  Frame 0 is the configuration before the first step.
    """
    run = Run(directory=directory, prefix=prefix)

    meshpoint = run.path(f"meshpoint{prefix}.csv")
    rows = []
    with open(meshpoint) as handle:
        for line in handle:
            line = line.strip().rstrip(",")
            if not line:
                continue
            rows.append(np.fromstring(line, sep=","))
    if not rows:
        raise ValueError(f"{meshpoint} has no frames")
    # The modal width, not the minimum: a run still in progress has a partial
    # last line, and letting it set the width would silently truncate every
    # frame to however far the writer had got.
    widths = {}
    for r in rows:
        widths[len(r)] = widths.get(len(r), 0) + 1
    width = max(widths, key=widths.get)
    rows = [r for r in rows if len(r) == width]
    run.coords = np.array(rows).reshape(len(rows), -1, 3)
    run.frame_iterations = np.arange(len(rows)) * output_interval

    for path in glob.glob(run.path(f"{prefix}face_*.csv")):
        match = re.search(r"face_(\d+)\.csv$", path)
        if match:
            run.face_frames[int(match.group(1))] = _read_face_csv(path)
    if not run.face_frames:
        # A run with flips off writes only the setup-time connectivity, which
        # describes every frame of it.
        run.face_frames[0] = _read_face_csv(run.path(f"{prefix}face.csv"))

    flips = run.path(f"{prefix}EdgeFlips.csv")
    if os.path.exists(flips):
        table = np.genfromtxt(flips, delimiter=",", names=True)
        run.flip_iteration = np.atleast_1d(table["iteration"]).astype(int)
        run.flip_delta_energy = np.atleast_1d(table["deltaEnergy"])
        run.flip_accepted = np.atleast_1d(table["accepted"]).astype(bool)

    run.free_vertices = _free_vertices(run)
    return run


def _free_vertices(run: Run) -> np.ndarray:
    """Vertices whose trajectory is their own.

    Two kinds are not.  A ghost is frozen, so its displacement sequence is
    identically zero.  A periodic duplicate is overwritten from its partner
    every step by `postprocess_ghost_periodic()`, so its displacement sequence
    is its partner's exactly.  Averaging either into an MSD dilutes it toward a
    number about the boundary condition rather than about the membrane: on the
    100 nm sheet used for WP6, 240 of 525 vertices are ghost and another 64 are
    duplicates, so the free fraction is 42%.

    Both are found the same way and without knowing any of that geometry --
    group the vertices by their sequence of frame-to-frame displacements, drop
    the all-zero group, and keep one representative of each remaining group.
    Nothing here depends on which ring of a periodic sheet is which, so it
    works for a run with any boundary condition.
    """
    n = run.coords.shape[1]
    if run.coords.shape[0] < 3:
        return _free_vertices_from_types(run)

    steps = np.diff(run.coords, axis=0)                   # (frames-1, n, 3)
    signature = steps.transpose(1, 0, 2).reshape(n, -1)    # (n, 3*(frames-1))
    magnitude = np.linalg.norm(signature, axis=1)

    # A vertex that never moved is frozen. The scale is set by the mesh's own
    # motion rather than by an absolute number, so this does not need to know
    # the time step or the diffusion constant.
    typical = np.median(magnitude[magnitude > 0]) if np.any(magnitude > 0) else 0.0
    moving = np.nonzero(magnitude > 1e-6 * typical)[0]
    if moving.size == 0:
        return _free_vertices_from_types(run)

    # Among the rest, a duplicate's displacement sequence is its partner's --
    # written out at eight significant figures, so not bit-identical, but the
    # correlation is 1 to many digits where two independent vertices correlate
    # at essentially zero. Anything above 0.999 is a copy; nothing real comes
    # close, so the threshold is not a tuned number.
    unit = signature[moving] / magnitude[moving][:, None]
    keep, taken = [], np.zeros(moving.size, dtype=bool)
    for i in range(moving.size):
        if taken[i]:
            continue
        keep.append(int(moving[i]))
        taken |= (unit @ unit[i]) > 0.999
    return np.array(keep) if keep else _free_vertices_from_types(run)


def _free_vertices_from_types(run: Run) -> np.ndarray:
    """Fallback for a run too short to tell trajectories apart: the type column.

    Column 3 of the vertex-type file is VertexType; 0 is a real interior vertex
    and 4 is a ghost.  This cannot see periodic duplicates, which are real
    vertices, so it is the weaker test and only used when the other cannot run.
    """
    types = run.path(f"{run.prefix}vertex_type_begin.csv")
    n = run.coords.shape[1]
    if not os.path.exists(types):
        return np.arange(n)
    kinds = []
    with open(types) as handle:
        for line in handle:
            fields = [f.strip() for f in line.strip().split(",") if f.strip()]
            if len(fields) >= 4:
                kinds.append(int(fields[3]))
    if len(kinds) != n:
        return np.arange(n)
    keep = [i for i, k in enumerate(kinds) if k == 0]
    return np.array(keep if keep else range(n))


def truncate(run: Run, before_iteration: int) -> Run:
    """The part of a run before `before_iteration`, as a new Run.

    For a run that diverged: the integrator writes a few frames of astronomical
    coordinates on its way to NaN, and one of them in a time average is enough
    to make every number about the run meaningless.  Keep what was written
    while the energy was still finite, which the divergence guard reports.
    """
    keep = run.frame_iterations < before_iteration
    out = Run(directory=run.directory, prefix=run.prefix)
    out.coords = run.coords[keep]
    out.frame_iterations = run.frame_iterations[keep]
    out.face_frames = {it: f for it, f in run.face_frames.items() if it < before_iteration}
    if run.flip_iteration is not None:
        m = run.flip_iteration < before_iteration
        out.flip_iteration = run.flip_iteration[m]
        out.flip_delta_energy = run.flip_delta_energy[m]
        out.flip_accepted = run.flip_accepted[m]
    out.free_vertices = _free_vertices(out)
    return out


# ---------------------------------------------------------------------------
# 1. Acceptance and valences
# ---------------------------------------------------------------------------

def acceptance_summary(run: Run, n_windows: int = 4) -> dict:
    """Acceptance and accepted-energy statistics, split over the run.

    The split matters: a single number averages the relaxation transient
    together with whatever came after it, and the question of interest is
    whether the accepted energy differences settle around zero.
    """
    if run.flip_iteration is None or run.flip_iteration.size == 0:
        return {"attempts": 0, "accepted": 0, "acceptance": float("nan"), "windows": []}

    last = int(run.flip_iteration.max()) + 1
    edges = np.linspace(0, last, n_windows + 1).astype(int)
    windows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        inside = (run.flip_iteration >= lo) & (run.flip_iteration < hi)
        taken = inside & run.flip_accepted
        windows.append({
            "from": int(lo), "to": int(hi),
            "attempts": int(inside.sum()),
            "accepted": int(taken.sum()),
            "acceptance": float(taken.sum() / inside.sum()) if inside.sum() else float("nan"),
            "mean_delta_energy": float(run.flip_delta_energy[taken].mean()) if taken.any() else float("nan"),
            "median_delta_energy": float(np.median(run.flip_delta_energy[taken])) if taken.any() else float("nan"),
        })
    return {
        "attempts": int(run.flip_iteration.size),
        "accepted": int(run.flip_accepted.sum()),
        "acceptance": float(run.flip_accepted.mean()),
        "windows": windows,
    }


def interior_vertices(faces: np.ndarray) -> np.ndarray:
    """Vertices whose fan of faces closes.

    The generated periodic sheet is a rectangle whose periodicity is realised
    by duplicated rings rather than by wrapping the face list, so the vertices
    on its geometric edge have partial fans -- two or three incident faces on
    the pristine grid, before any flip.  Those are not extraordinary vertices
    and counting them into a valence histogram makes an equilibrium fluid mesh
    look pathological: on the 100 nm sheet they drag the mean valence from 6.0
    to 5.6 at step zero.

    A fan closes exactly when the number of incident faces equals the number of
    distinct neighbours, which needs no knowledge of how the sheet was built.
    """
    n = int(faces.max()) + 1
    incident_faces = np.zeros(n, dtype=int)
    neighbours = [set() for _ in range(n)]
    for tri in faces:
        for k in range(3):
            v = int(tri[k])
            incident_faces[v] += 1
            neighbours[v].add(int(tri[(k + 1) % 3]))
            neighbours[v].add(int(tri[(k + 2) % 3]))
    return np.array([v for v in range(n)
                     if incident_faces[v] > 0 and incident_faces[v] == len(neighbours[v])])


def valence_histogram(faces: np.ndarray, vertices: np.ndarray | None = None) -> dict:
    """Valence counts over the interior vertices, intersected with `vertices`.

    The valence of an interior vertex is its number of incident faces, which is
    what this counts.  Restricted to interior vertices because a partial fan
    has no valence in the sense the DTS 60/20/20 band is about.
    """
    n = int(faces.max()) + 1
    valence = np.zeros(n, dtype=int)
    for k in range(3):
        np.add.at(valence, faces[:, k], 1)

    inside = interior_vertices(faces)
    if vertices is not None:
        inside = np.intersect1d(inside, vertices)
    if inside.size == 0:
        return {"counts": {}, "fractions": {}, "mean": float("nan"), "n": 0}
    valence = valence[inside]

    counts = {}
    for v in valence:
        counts[int(v)] = counts.get(int(v), 0) + 1
    total = sum(counts.values())
    return {"counts": dict(sorted(counts.items())),
            "fractions": {k: v / total for k, v in sorted(counts.items())},
            "mean": float(valence.mean()), "n": int(inside.size)}


# ---------------------------------------------------------------------------
# 2. Neighbour survival
# ---------------------------------------------------------------------------

def _edge_set(faces: np.ndarray) -> set:
    edges = set()
    for tri in faces:
        for k in range(3):
            a, b = int(tri[k]), int(tri[(k + 1) % 3])
            edges.add((a, b) if a < b else (b, a))
    return edges


def neighbour_survival(run: Run) -> tuple[np.ndarray, np.ndarray]:
    """Fraction of the initial edges still present, against iteration.

    Returns (iterations, fraction).  Exactly 1 everywhere for a run with flips
    off, which is the control the fluid case is read against.
    """
    iterations = sorted(run.face_frames)
    initial = _edge_set(run.face_frames[iterations[0]])
    if not initial:
        return np.array([]), np.array([])
    times, fractions = [], []
    for it in iterations:
        present = _edge_set(run.face_frames[it])
        times.append(it)
        fractions.append(len(initial & present) / len(initial))
    return np.array(times), np.array(fractions)


def survival_time(times: np.ndarray, fractions: np.ndarray,
                  threshold: float = 1.0 / np.e) -> float:
    """The time at which survival first falls to `threshold`, interpolated.

    NaN when the run never gets there, which is the honest answer: a decay
    constant fitted to a curve that has barely moved is a number about the fit,
    not about the membrane.
    """
    below = np.nonzero(fractions <= threshold)[0]
    if below.size == 0:
        return float("nan")
    i = int(below[0])
    if i == 0:
        return float(times[0])
    x0, x1 = fractions[i - 1], fractions[i]
    if x0 == x1:
        return float(times[i])
    return float(times[i - 1] + (times[i] - times[i - 1]) * (x0 - threshold) / (x0 - x1))


# ---------------------------------------------------------------------------
# 3. In-plane mean squared displacement
# ---------------------------------------------------------------------------

def in_plane_msd(run: Run, lags: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Time-averaged in-plane MSD over the free vertices.

    Averaged over every pair of frames a lag apart as well as over vertices,
    which is what makes a few hundred frames enough to see the shape.  The
    membrane's own drift is removed first: the whole sheet wandering is a
    centre-of-mass motion, not a rearrangement, and on a periodic box it is not
    even physical.
    """
    xy = run.coords[:, run.free_vertices, :2]
    xy = xy - xy.mean(axis=1, keepdims=True)

    n_frames = xy.shape[0]
    if lags is None:
        top = max(1, n_frames // 2)
        lags = np.unique(np.geomspace(1, top, num=min(24, top)).astype(int))

    step = run.frame_iterations[1] - run.frame_iterations[0] if n_frames > 1 else 1
    out = []
    for lag in lags:
        if lag >= n_frames:
            out.append(np.nan)
            continue
        difference = xy[lag:] - xy[:-lag]
        out.append(float((difference ** 2).sum(axis=2).mean()))
    return lags * step, np.array(out)


def msd_growth_exponent(times: np.ndarray, msd: np.ndarray) -> float:
    """Slope of log(MSD) against log(t), over the upper half of the range.

    1 is diffusive, 0 is a saturated cage.  Fitted on the tail because the
    first few lags are the ballistic/cage-filling part in every case and say
    nothing about which of the two this is.
    """
    good = np.isfinite(msd) & (msd > 0) & (times > 0)
    t, m = times[good], msd[good]
    if t.size < 4:
        return float("nan")
    half = t.size // 2
    return float(np.polyfit(np.log(t[half:]), np.log(m[half:]), 1)[0])


# ---------------------------------------------------------------------------
# Edge lengths -- what the tether is actually doing
# ---------------------------------------------------------------------------

def edge_length_stats(run: Run, frame: int = -1, interior_only: bool = True) -> dict:
    """Edge-length distribution at one frame, against the tether's range.

    Restricted by default to edges whose endpoints are both interior and free.
    The ghost band is a copy positioned by the periodic map rather than by the
    dynamics, and its edges are not the ones the tether is acting on -- mixing
    them in makes the distribution look far wider than the membrane's.
    """
    iterations = sorted(run.face_frames)
    target = run.frame_iterations[frame] if frame >= 0 else run.frame_iterations[-1]
    usable = [it for it in iterations if it <= target] or [iterations[0]]
    faces = run.face_frames[usable[-1]]
    points = run.coords[frame]

    keep = None
    if interior_only:
        keep = set(np.intersect1d(interior_vertices(faces), run.free_vertices).tolist())
    lengths = []
    for a, b in _edge_set(faces):
        if keep is not None and (a not in keep or b not in keep):
            continue
        lengths.append(np.linalg.norm(points[a] - points[b]))
    lengths = np.array(lengths)
    if lengths.size == 0:
        return {"n": 0}
    return {"n": lengths.size, "mean": float(lengths.mean()),
            "std": float(lengths.std()), "min": float(lengths.min()),
            "max": float(lengths.max()),
            "percentiles": {p: float(np.percentile(lengths, p)) for p in (1, 25, 50, 75, 99)}}


def report(run: Run) -> str:
    """Everything above, as text for a results document."""
    lines = [f"run: {run.directory}",
             f"  frames {run.coords.shape[0]}, vertices {run.coords.shape[1]} "
             f"({run.free_vertices.size} free), face frames {len(run.face_frames)}"]

    summary = acceptance_summary(run)
    if summary["attempts"]:
        lines.append(f"  flips: {summary['accepted']}/{summary['attempts']} "
                     f"= {summary['acceptance']:.3f} acceptance")
        for w in summary["windows"]:
            lines.append(f"    steps {w['from']:6d}-{w['to']:6d}: "
                         f"{w['accepted']:4d}/{w['attempts']:4d} "
                         f"({w['acceptance']:.3f}), mean dE {w['mean_delta_energy']:+.3f}, "
                         f"median {w['median_delta_energy']:+.3f}")
    else:
        lines.append("  flips: none (control run)")

    last = sorted(run.face_frames)[-1]
    hist = valence_histogram(run.face_frames[last], run.free_vertices)
    lines.append(f"  valences at step {last} ({hist['n']} interior free vertices, "
                 f"mean {hist['mean']:.3f}): " +
                 ", ".join(f"{k}:{v:.3f}" for k, v in hist["fractions"].items()))

    times, fractions = neighbour_survival(run)
    if times.size > 1:
        lines.append(f"  neighbour survival: {fractions[-1]:.3f} at step {times[-1]}, "
                     f"1/e time {survival_time(times, fractions):.1f}")

    t, msd = in_plane_msd(run)
    lines.append(f"  in-plane MSD: {msd[0]:.4f} nm^2 at lag {t[0]}, "
                 f"{msd[-1]:.4f} at lag {t[-1]}, growth exponent "
                 f"{msd_growth_exponent(t, msd):.3f}")

    stats = edge_length_stats(run)
    lines.append(f"  edge lengths: mean {stats['mean']:.3f} +- {stats['std']:.3f} nm, "
                 f"range [{stats['min']:.3f}, {stats['max']:.3f}], "
                 f"1-99% [{stats['percentiles'][1]:.3f}, {stats['percentiles'][99]:.3f}]")
    return "\n".join(lines)


if __name__ == "__main__":
    import sys
    interval = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    print(report(load_run(sys.argv[1], output_interval=interval)))
