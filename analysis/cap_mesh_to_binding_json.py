#!/usr/bin/env python3
"""Build a coarse-grained binding JSON from a spherical-cap point cloud."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Delaunay


def load_points(path):
    points = np.loadtxt(path, delimiter=",", dtype=float)
    if points.ndim == 1:
        points = points.reshape(1, -1)
    if points.shape[1] != 3:
        raise ValueError(f"expected an Nx3 CSV, got shape {points.shape}")
    if points.shape[0] < 4:
        raise ValueError("at least four points are required for sphere fitting")
    return points


def fit_sphere(points):
    """Least-squares fit of x^2 + y^2 + z^2 + ax + by + cz + d = 0."""
    lhs = np.column_stack([points, np.ones(points.shape[0])])
    rhs = -np.sum(points * points, axis=1)
    coeffs, *_ = np.linalg.lstsq(lhs, rhs, rcond=None)
    center = -0.5 * coeffs[:3]
    radius_sq = float(np.dot(center, center) - coeffs[3])
    if radius_sq <= 0:
        raise ValueError("sphere fit produced a non-positive radius")
    return center, radius_sq**0.5


def orthonormal_basis(axis):
    axis = axis / np.linalg.norm(axis)
    helper = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(axis, helper))) > 0.9:
        helper = np.array([0.0, 1.0, 0.0])
    u = helper - np.dot(helper, axis) * axis
    u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    v /= np.linalg.norm(v)
    return u, v


def project_for_delaunay(points, center):
    cap_axis = points.mean(axis=0) - center
    cap_axis_norm = np.linalg.norm(cap_axis)
    if cap_axis_norm == 0:
        raise ValueError("could not infer cap axis from point cloud")
    cap_axis /= cap_axis_norm
    u, v = orthonormal_basis(cap_axis)
    relative = points - center
    projected = np.column_stack([relative @ u, relative @ v])
    return projected, cap_axis, u, v


def triangulate(projected, points, max_edge_length=None, max_edge_factor=None):
    simplices = Delaunay(projected).simplices

    if max_edge_factor is not None:
        all_lengths = []
        for simplex in simplices:
            for a_pos, b_pos in ((0, 1), (1, 2), (2, 0)):
                all_lengths.append(
                    np.linalg.norm(points[simplex[a_pos]] - points[simplex[b_pos]])
                )
        factor_limit = float(np.median(all_lengths) * max_edge_factor)
        max_edge_length = (
            factor_limit
            if max_edge_length is None
            else min(float(max_edge_length), factor_limit)
        )

    kept = []
    skipped = []
    for simplex in simplices:
        lengths = [
            np.linalg.norm(points[simplex[a_pos]] - points[simplex[b_pos]])
            for a_pos, b_pos in ((0, 1), (1, 2), (2, 0))
        ]
        if max_edge_length is not None and max(lengths) > max_edge_length:
            skipped.append([int(index) for index in simplex])
            continue
        kept.append([int(index) for index in simplex])

    edges = set()
    for a, b, c in kept:
        edges.add(tuple(sorted((a, b))))
        edges.add(tuple(sorted((b, c))))
        edges.add(tuple(sorted((c, a))))

    return kept, sorted(edges), skipped


def connection_line_interface(points, vertex_index, neighbor_index, distance):
    com = points[vertex_index]
    direction = points[neighbor_index] - com
    direction_norm = np.linalg.norm(direction)
    if direction_norm == 0:
        raise ValueError(f"edge {vertex_index}-{neighbor_index} has zero length")
    direction /= direction_norm
    return com + distance * direction, direction


def build_json(points, center, radius, projected, cap_axis, triangles, edges, args):
    molecules = {
        str(index): {
            "mol_index": index,
            "com": points[index].tolist(),
            "interfaces": {},
        }
        for index in range(points.shape[0])
    }

    adjacency = {index: [] for index in range(points.shape[0])}
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)

    for vertex_index, neighbors in adjacency.items():
        for iface_index, neighbor_index in enumerate(sorted(neighbors)):
            coord, connection_direction = connection_line_interface(
                points, vertex_index, neighbor_index, args.interface_distance
            )
            name = f"to_{neighbor_index}"
            molecules[str(vertex_index)]["interfaces"][name] = {
                "name": name,
                "iface_index": iface_index,
                "coord": coord.tolist(),
                "is_bound": True,
                "neighbor_molecule_index": neighbor_index,
                "connection_direction": connection_direction.tolist(),
            }

    bound_states = []
    for left, right in edges:
        left_iface = molecules[str(left)]["interfaces"][f"to_{right}"]
        right_iface = molecules[str(right)]["interfaces"][f"to_{left}"]
        pair = {
            "mol1": {
                "molecule_index": left,
                "interface_name": f"to_{right}",
                "interface_index": left_iface["iface_index"],
            },
            "mol2": {
                "molecule_index": right,
                "interface_name": f"to_{left}",
                "interface_index": right_iface["iface_index"],
            },
        }
        bound_states.append(pair)
        if args.bidirectional:
            bound_states.append(
                {
                    "mol1": pair["mol2"],
                    "mol2": pair["mol1"],
                }
            )

    return {
        "metadata": {
            "source_csv": str(args.input),
            "interface_distance": args.interface_distance,
            "sphere_center": center.tolist(),
            "sphere_radius": radius,
            "cap_axis": cap_axis.tolist(),
            "delaunay_projection": "best-fit cap plane through fitted sphere center",
            "interface_placement": "along COM-COM connection line",
            "vertex_count": int(points.shape[0]),
            "triangle_count": len(triangles),
            "edge_count": len(edges),
            "bidirectional_bound_states": args.bidirectional,
        },
        "mesh": {
            "vertices": points.tolist(),
            "projected_vertices": projected.tolist(),
            "triangles": triangles,
            "edges": [[left, right] for left, right in edges],
        },
        "selected_complex": {
            "complex_index": 0,
            "size": int(points.shape[0]),
            "member_molecule_indices": list(range(points.shape[0])),
            "molecules": molecules,
            "bound_states": bound_states,
        },
    }


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


def save_plot(path, points, molecules, edges):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")

    com_to_interface_segments = []
    interface_pair_segments = []
    interface_points = []

    for molecule_index, molecule in molecules.items():
        com = np.array(molecule["com"], dtype=float)
        for interface in molecule["interfaces"].values():
            coord = np.array(interface["coord"], dtype=float)
            com_to_interface_segments.append([com, coord])
            interface_points.append(coord)

    for left, right in edges:
        left_coord = np.array(
            molecules[str(left)]["interfaces"][f"to_{right}"]["coord"], dtype=float
        )
        right_coord = np.array(
            molecules[str(right)]["interfaces"][f"to_{left}"]["coord"], dtype=float
        )
        interface_pair_segments.append([left_coord, right_coord])

    if interface_pair_segments:
        ax.add_collection3d(
            Line3DCollection(
                interface_pair_segments,
                colors="0.65",
                linewidths=0.8,
                alpha=0.7,
            )
        )

    if com_to_interface_segments:
        ax.add_collection3d(
            Line3DCollection(
                com_to_interface_segments,
                colors="black",
                linewidths=0.7,
                alpha=0.9,
            )
        )

    if interface_points:
        interface_points = np.vstack(interface_points)
        scale_points = np.vstack([points, interface_points])
    else:
        scale_points = points

    ax.scatter(points[:, 0], points[:, 1], points[:, 2], s=20, color="#1f77b4")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    set_equal_3d_axes(ax, scale_points)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Fit a sphere to a cap point cloud, Delaunay-triangulate it, and "
            "emit molecule/interface binding JSON."
        )
    )
    parser.add_argument("input", type=Path, help="Input Nx3 CSV point cloud.")
    parser.add_argument("output", type=Path, help="Output JSON path.")
    parser.add_argument(
        "--interface-distance",
        type=float,
        default=1.5,
        help="Distance from each COM to each interface along the COM-COM line.",
    )
    parser.add_argument(
        "--max-edge-length",
        type=float,
        default=None,
        help="Optional absolute cutoff for triangle edge lengths.",
    )
    parser.add_argument(
        "--max-edge-factor",
        type=float,
        default=None,
        help="Optional cutoff as a multiple of the median Delaunay edge length.",
    )
    parser.add_argument(
        "--one-way",
        action="store_true",
        help="Write one bound_state per mesh edge instead of both directions.",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=None,
        help="Optional path for a 3D equal-scale scatter/mesh plot.",
    )
    args = parser.parse_args()
    if args.interface_distance <= 0:
        parser.error("--interface-distance must be positive")
    if args.max_edge_length is not None and args.max_edge_length <= 0:
        parser.error("--max-edge-length must be positive")
    if args.max_edge_factor is not None and args.max_edge_factor <= 0:
        parser.error("--max-edge-factor must be positive")
    args.bidirectional = not args.one_way
    return args


def main():
    args = parse_args()
    points = load_points(args.input)
    center, radius = fit_sphere(points)
    projected, cap_axis, _, _ = project_for_delaunay(points, center)
    triangles, edges, skipped = triangulate(
        projected,
        points,
        max_edge_length=args.max_edge_length,
        max_edge_factor=args.max_edge_factor,
    )
    data = build_json(points, center, radius, projected, cap_axis, triangles, edges, args)
    if skipped:
        data["mesh"]["skipped_triangles"] = skipped
        data["metadata"]["skipped_triangle_count"] = len(skipped)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2) + "\n")
    if args.plot is not None:
        save_plot(args.plot, points, data["selected_complex"]["molecules"], edges)

    print(f"wrote {args.output}")
    if args.plot is not None:
        print(f"wrote plot {args.plot}")
    print(f"vertices: {points.shape[0]}")
    print(f"triangles: {len(triangles)}")
    print(f"edges: {len(edges)}")
    print(f"sphere_center: {center.tolist()}")
    print(f"sphere_radius: {radius}")


if __name__ == "__main__":
    main()
