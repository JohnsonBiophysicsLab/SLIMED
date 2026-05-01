#!/usr/bin/env python3

import argparse
import json
import math
from pathlib import Path


PAIR_SPECS = {
    "d-p": ("d", "p", 0.7, 1.5, 0.849084),
    "c-q": ("c", "q", 0.6, 1.4, 0.729108),
    "b-b": ("b", "b", 0.3, 1.1, 0.459582),
}


def midpoint(coord_a, coord_b):
    return [(a + b) / 2.0 for a, b in zip(coord_a, coord_b)]


def distance(coord_a, coord_b):
    return math.dist(coord_a, coord_b)


def build_candidates(molecules, left_name, right_name, lo, hi, target):
    left_nodes = []
    right_nodes = []
    same_type = left_name == right_name

    for molecule_id, molecule in molecules.items():
        mol_index = int(molecule_id)
        left_iface = molecule["interfaces"].get(left_name)
        right_iface = molecule["interfaces"].get(right_name)

        if left_iface and not left_iface["is_bound"]:
            left_nodes.append((mol_index, molecule["com"], left_iface["coord"]))
        if right_iface and not right_iface["is_bound"]:
            right_nodes.append((mol_index, molecule["com"], right_iface["coord"]))

    candidates = []
    for left_index, left_com, left_coord in left_nodes:
        for right_index, right_com, right_coord in right_nodes:
            if left_index == right_index:
                continue
            if same_type and left_index > right_index:
                continue

            pair_distance = distance(left_coord, right_coord)
            if not (lo < pair_distance < hi):
                continue

            candidates.append(
                {
                    "score": abs(pair_distance - target),
                    "distance": pair_distance,
                    "left_index": left_index,
                    "right_index": right_index,
                    "left_com": left_com,
                    "right_com": right_com,
                    "left_coord": left_coord,
                    "right_coord": right_coord,
                    "midpoint": midpoint(left_coord, right_coord),
                }
            )

    candidates.sort(key=lambda item: (item["score"], item["distance"]))
    return candidates


def greedy_select(candidates, same_type=False):
    used_left = set()
    used_right = set()
    chosen = []

    for candidate in candidates:
        left_index = candidate["left_index"]
        right_index = candidate["right_index"]

        if left_index in used_left or right_index in used_right:
            continue
        if same_type and (left_index in used_right or right_index in used_left):
            continue

        used_left.add(left_index)
        used_right.add(right_index)
        chosen.append(candidate)

    return chosen


def main():
    parser = argparse.ArgumentParser(
        description="Add local crack-repair interface pairs to a structure JSON file."
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--seed-radius",
        type=float,
        default=5.0,
        help="Keep c-q midpoints within this radius of a d-p seed midpoint.",
    )
    args = parser.parse_args()

    data = json.loads(args.input.read_text())
    selected_complex = data["selected_complex"]
    molecules = selected_complex["molecules"]

    dp_candidates = build_candidates(molecules, *PAIR_SPECS["d-p"])
    chosen_dp = greedy_select(dp_candidates, same_type=False)
    seed_midpoints = [candidate["midpoint"] for candidate in chosen_dp]

    cq_candidates = build_candidates(molecules, *PAIR_SPECS["c-q"])
    local_cq_candidates = []
    for candidate in cq_candidates:
        if not seed_midpoints:
            continue
        midpoint_distance = min(
            distance(candidate["midpoint"], seed_midpoint)
            for seed_midpoint in seed_midpoints
        )
        if midpoint_distance <= args.seed_radius:
            candidate = dict(candidate)
            candidate["seed_distance"] = midpoint_distance
            local_cq_candidates.append(candidate)

    local_cq_candidates.sort(
        key=lambda item: (item["seed_distance"], item["score"], item["distance"])
    )
    chosen_cq = greedy_select(local_cq_candidates, same_type=False)

    # The visible crack in the screenshot sits near the d-p seed cluster, so we
    # intentionally do not add distant b-b bonds from elsewhere in the shell.
    chosen_pairs = [
        ("d", "p", chosen_dp),
        ("c", "q", chosen_cq),
    ]

    summary = []
    added_pairs = []
    for left_name, right_name, chosen in chosen_pairs:
        for candidate in chosen:
            left_molecule = molecules[str(candidate["left_index"])]
            right_molecule = molecules[str(candidate["right_index"])]
            left_iface = left_molecule["interfaces"][left_name]
            right_iface = right_molecule["interfaces"][right_name]

            left_iface["is_bound"] = True
            right_iface["is_bound"] = True

            selected_complex["bound_states"].append(
                {
                    "mol1": {
                        "molecule_index": candidate["left_index"],
                        "interface_name": left_name,
                        "interface_index": left_iface["iface_index"],
                    },
                    "mol2": {
                        "molecule_index": candidate["right_index"],
                        "interface_name": right_name,
                        "interface_index": right_iface["iface_index"],
                    },
                }
            )
            selected_complex["bound_states"].append(
                {
                    "mol1": {
                        "molecule_index": candidate["right_index"],
                        "interface_name": right_name,
                        "interface_index": right_iface["iface_index"],
                    },
                    "mol2": {
                        "molecule_index": candidate["left_index"],
                        "interface_name": left_name,
                        "interface_index": left_iface["iface_index"],
                    },
                }
            )

            added_pairs.append(
                {
                    "pair_type": f"{left_name}-{right_name}",
                    "mol1": candidate["left_index"],
                    "mol2": candidate["right_index"],
                    "distance": candidate["distance"],
                    "distance_delta": candidate["score"],
                    "midpoint": candidate["midpoint"],
                }
            )

        summary.append((f"{left_name}-{right_name}", len(chosen)))

    args.output.write_text(json.dumps(data, indent=2))

    print(f"wrote {args.output}")
    print("summary", summary)
    print("added_pairs", len(added_pairs))
    for pair in added_pairs:
        rounded = {
            "pair_type": pair["pair_type"],
            "mol1": pair["mol1"],
            "mol2": pair["mol2"],
            "distance": round(pair["distance"], 6),
            "distance_delta": round(pair["distance_delta"], 6),
            "midpoint": tuple(round(value, 3) for value in pair["midpoint"]),
        }
        print(rounded)


if __name__ == "__main__":
    main()
