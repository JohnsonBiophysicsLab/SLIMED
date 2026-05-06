#!/usr/bin/env python3
"""Plot Gag scaffold COMs and membrane vertices in one equal-scale 3D view."""

import argparse
import json
from pathlib import Path

import numpy as np


def load_vertices(path):
    vertices = np.loadtxt(path, delimiter=",", dtype=float)
    if vertices.ndim == 1:
        vertices = vertices.reshape(1, -1)
    if vertices.shape[1] != 3:
        raise ValueError(f"expected an Nx3 vertex CSV, got shape {vertices.shape}")
    return vertices


def load_gag_coms(path):
    data = json.loads(path.read_text())
    if "molecules" in data:
        molecules = data["molecules"]
    else:
        molecules = data["selected_complex"]["molecules"]

    coms = []
    for molecule_id in sorted(molecules, key=int):
        coms.append(molecules[molecule_id]["com"])

    coms = np.array(coms, dtype=float)
    if coms.ndim != 2 or coms.shape[1] != 3:
        raise ValueError(f"expected Gag COMs to be Nx3, got shape {coms.shape}")
    return coms


def set_equal_3d_axes(ax, points):
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    center = (mins + maxs) / 2.0
    half_range = float(np.max(maxs - mins) / 2.0)
    if half_range == 0:
        half_range = 1.0

    ax.set_xlim(center[0] - half_range, center[0] + half_range)
    ax.set_ylim(center[1] - half_range, center[1] + half_range)
    ax.set_zlim(center[2] - half_range, center[2] + half_range)
    ax.set_box_aspect((1, 1, 1))


def plot(gag_coms, vertices, output=None, show=False):
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")

    ax.scatter(
        vertices[:, 0],
        vertices[:, 1],
        vertices[:, 2],
        s=4,
        c="#4f81bd",
        alpha=0.45,
        label="vertices",
        depthshade=False,
    )
    ax.scatter(
        gag_coms[:, 0],
        gag_coms[:, 1],
        gag_coms[:, 2],
        s=14,
        c="#d62728",
        alpha=0.95,
        label="Gag scaffold COM",
        depthshade=False,
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.legend(loc="best")
    set_equal_3d_axes(ax, np.vstack([vertices, gag_coms]))
    fig.tight_layout()

    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=200)
        print(f"wrote {output}")

    if show:
        plt.show()
    else:
        plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot Gag scaffold COMs and membrane vertices in one 3D plot."
    )
    parser.add_argument("gag_scaffold", type=Path, help="Gag scaffold .dat/.json file.")
    parser.add_argument("vertices", type=Path, help="Three-column vertex CSV file.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Optional image output path, e.g. plot.png.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Open an interactive matplotlib window after plotting.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    gag_coms = load_gag_coms(args.gag_scaffold)
    vertices = load_vertices(args.vertices)
    plot(gag_coms, vertices, output=args.output, show=args.show)
    print(f"Gag COMs: {len(gag_coms)}")
    print(f"vertices: {len(vertices)}")


if __name__ == "__main__":
    main()
