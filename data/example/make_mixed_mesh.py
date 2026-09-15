#!/usr/bin/env python3
"""Write a SLIMED mesh with per-vertex boundary conditions (boundaryType = Mixed).

Two things this produces that the built-in flat-sheet generator cannot:

* a patch whose two axes carry different boundaries -- the generator can do
  that too, through boundaryTypeX / boundaryTypeY, and this script reproduces
  it so the two can be compared;
* a patch with a **hole**, whose rim is a free boundary while the outer border
  stays periodic. That is the bottleneck case: a membrane that wraps in the
  plane and is open at a neck.

The mesh is the same triangular lattice `Mesh::set_vertices_faces_flat()`
builds, with the same index convention (`vertex = (nFaceX + 1) * j + i`), the
same three-rings-plus-a-duplicate periodic band, and the same clockwise
winding the exporter writes -- so a generated sheet and a sheet written by
`export_mesh_to_vertices_faces()` are the same object.

Carving the hole is done on *physical* faces: a face and every periodic copy
of it are removed together, so the seam cannot end up with a hole on one side
and membrane on the other. Vertices left with no face are dropped and the
mirror indices renumbered.

    ./make_mixed_mesh.py --out data/example/mixed_neck \
        --side 120 --pore-radius 12 --neck-height 8 --neck-width 25

See docs/mixed_boundary_conditions.md for the file format and for what the
boundary types mean.
"""

import argparse
import math
import sys

# Vertex types, as the mesh files spell them (include/mesh/Vertex.hpp).
FREE = "free"
FIXED = "fixed"
PERIODIC = "periodic"


def c_round(value):
    """C's round(): halfway cases away from zero, not Python's to-even."""
    return int(math.floor(value + 0.5)) if value >= 0 else -int(math.floor(-value + 0.5))


def axes_division(side_x, side_y, l_face):
    """Mesh::set_axes_division_flat(), including the forced-even row count."""
    n_face_x = c_round(side_x / l_face)
    d_face_x = side_x / n_face_x
    d_face_y = math.sqrt(3.0) / 2.0 * d_face_x
    n_face_y = c_round(side_y / d_face_y)
    if n_face_y % 2:
        # An odd row count leaves the zig-zag half a cell out of step after one
        # period, so a y-periodic image would not be a pure translation.
        n_face_y += 1
    return n_face_x, n_face_y, d_face_x, d_face_y


def build_lattice(n_face_x, n_face_y, d_face_x, d_face_y):
    """Vertex positions and triangles, in Mesh::set_vertices_faces_flat() order."""
    half_x = n_face_x * d_face_x * 0.5
    half_y = n_face_y * d_face_y * 0.5

    coords = []
    for j in range(n_face_y + 1):
        for i in range(n_face_x + 1):
            x = i * d_face_x + (d_face_x / 2.0 if j % 2 == 0 else 0.0) - half_x
            coords.append([x, j * d_face_y - half_y, 0.0])

    def index(i, j):
        return (n_face_x + 1) * j + i

    faces = []
    cell_of_face = []
    for j in range(n_face_y):
        for i in range(n_face_x):
            if j % 2:
                node1, node2 = index(i, j), index(i, j + 1)
                node3, node4 = index(i + 1, j), index(i + 1, j + 1)
            else:
                node3, node1 = index(i, j), index(i, j + 1)
                node4, node2 = index(i + 1, j), index(i + 1, j + 1)
            faces.append([node1, node2, node3])
            faces.append([node2, node4, node3])
            cell_of_face.extend([(i, j), (i, j)])
    return coords, faces, cell_of_face


def orient_clockwise(coords, faces):
    """Wind every triangle clockwise in the xy plane.

    A sheet is a graph over the plane, so one sign for every triangle is a
    consistent orientation -- every shared edge is then traversed in opposite
    directions by its two faces, which is what the one-ring walk needs and
    what the importer does not fix for itself. Clockwise is what
    export_mesh_to_vertices_faces() writes for the generated sheet.
    """
    for face in faces:
        a, b, c = (coords[k] for k in face)
        twice_area = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        if twice_area > 0.0:
            face[1], face[2] = face[2], face[1]


def assign_boundary_types(n_face_x, n_face_y, boundary_x, boundary_y, fixed_rings):
    """Mesh::set_boundary_types_flat_mixed(), in Python.

    Sources occupy indices 3 .. nFace - 4 along a periodic axis; everything
    outside that is an image of the vertex one period in. Three image rings
    and a duplicate ring on each side is what gives every real face a
    complete one-ring.
    """
    n_face = (n_face_x, n_face_y)
    boundary = (boundary_x, boundary_y)
    for axis in (0, 1):
        if boundary[axis] == "periodic" and n_face[axis] < 10:
            raise SystemExit(
                "a periodic axis needs at least 10 faces (three image rings and a "
                "duplicate ring on each side, and a period of at least four); axis "
                "%d has %d. Enlarge --side or reduce --l-face." % (axis, n_face[axis])
            )
    period = (n_face_x - 6, n_face_y - 6)

    def source_of(index, axis):
        if boundary[axis] != "periodic":
            return index
        return 3 + (index - 3) % period[axis]

    types = []
    mirrors = []
    for j in range(n_face_y + 1):
        for i in range(n_face_x + 1):
            source_i = source_of(i, 0)
            source_j = source_of(j, 1)
            if (source_i, source_j) != (i, j):
                types.append(PERIODIC)
                mirrors.append((n_face_x + 1) * source_j + source_i)
                continue
            # Only a source is ever clamped; an image of a clamped source is
            # held still through it.
            clamped_x = boundary_x == "fixed" and (i < fixed_rings or i > n_face_x - fixed_rings)
            clamped_y = boundary_y == "fixed" and (j < fixed_rings or j > n_face_y - fixed_rings)
            types.append(FIXED if (clamped_x or clamped_y) else FREE)
            mirrors.append(-1)
    return types, mirrors


def copy_flags(n_face_x, n_face_y, boundary_x, boundary_y, cell_of_face):
    """The faces of the image band, which duplicate a physical face."""
    flags = []
    for i, j in cell_of_face:
        copy = (boundary_x == "periodic" and (i < 3 or i > n_face_x - 4)) or (
            boundary_y == "periodic" and (j < 3 or j > n_face_y - 4)
        )
        flags.append(1 if copy else 0)
    return flags


def carve_pore(coords, types, mirrors, faces, flags, radius, centre):
    """Remove a disc of faces, leaving a free rim, and renumber what is left.

    The decision is taken on the physical face -- the triple of sources its
    corners stand for -- and applied to every copy of it, so a hole cannot
    appear on one side of a seam and not the other.
    """
    if radius <= 0.0:
        return coords, types, mirrors, faces, flags

    def source_of(vertex):
        return mirrors[vertex] if mirrors[vertex] >= 0 else vertex

    doomed_keys = set()
    keys = []
    for face in faces:
        sources = tuple(sorted(source_of(v) for v in face))
        keys.append(sources)
        centroid_x = sum(coords[v][0] for v in sources) / 3.0
        centroid_y = sum(coords[v][1] for v in sources) / 3.0
        if math.hypot(centroid_x - centre[0], centroid_y - centre[1]) < radius:
            doomed_keys.add(sources)

    kept_faces = [face for face, key in zip(faces, keys) if key not in doomed_keys]
    kept_flags = [flag for flag, key in zip(flags, keys) if key not in doomed_keys]
    if not kept_faces:
        raise SystemExit("--pore-radius removed every face")

    used = sorted({v for face in kept_faces for v in face})
    new_index = {old: new for new, old in enumerate(used)}
    for old in used:
        mirror = mirrors[old]
        if mirror >= 0 and mirror not in new_index:
            # Cannot happen while a face and its copies are removed together:
            # an image outlives its source only if some face still uses the
            # image while none uses the source.
            raise SystemExit(
                "vertex %d survived but its source %d did not" % (old, mirror)
            )
    new_coords = [coords[old] for old in used]
    new_types = [types[old] for old in used]
    new_mirrors = [new_index[mirrors[old]] if mirrors[old] >= 0 else -1 for old in used]
    new_faces = [[new_index[v] for v in face] for face in kept_faces]
    return new_coords, new_types, new_mirrors, new_faces, kept_flags


def lift_neck(coords, types, mirrors, radius, height, width, centre):
    """Raise the membrane into a collar that peaks at the rim of the hole.

    z = height * (1 - s^2)^2 with s = (r - radius) / width, and exactly zero
    beyond r = radius + width, so the far field is flat and a periodic image
    and its source sit at the same height whatever the geometry does in
    between. The images take their height from their sources rather than from
    the formula, so the offset between the two stays a pure in-plane
    translation -- which is what makes it a lattice period.
    """
    if height == 0.0 or width <= 0.0:
        return
    for vertex, coord in enumerate(coords):
        if types[vertex] == PERIODIC:
            continue
        r = math.hypot(coord[0] - centre[0], coord[1] - centre[1])
        s = (r - radius) / width
        coord[2] = height * (1.0 - s * s) ** 2 if 0.0 <= s <= 1.0 else 0.0
    for vertex, mirror in enumerate(mirrors):
        if mirror >= 0:
            coords[vertex][2] = coords[mirror][2]


def write_mesh(prefix, coords, types, mirrors, faces, flags):
    vertices_path = prefix + "_vertices.csv"
    faces_path = prefix + "_faces.csv"
    with open(vertices_path, "w") as out:
        out.write("# SLIMED mesh vertices: x, y, z, type, mirror\n")
        for coord, kind, mirror in zip(coords, types, mirrors):
            out.write(
                "%.17g,%.17g,%.17g,%s,%d\n" % (coord[0], coord[1], coord[2], kind, mirror)
            )
    with open(faces_path, "w") as out:
        out.write("# SLIMED mesh faces: v0, v1, v2, copy\n")
        for face, flag in zip(faces, flags):
            out.write("%d,%d,%d,%d\n" % (face[0], face[1], face[2], flag))
    return vertices_path, faces_path


def report(coords, types, faces, flags):
    """What the mesh turned out to be, in the terms the C++ setup will report."""
    incident = [0] * len(coords)
    neighbours = [set() for _ in coords]
    for face in faces:
        for k in range(3):
            incident[face[k]] += 1
            neighbours[face[k]].update(face[m] for m in range(3) if m != k)
    interior = [incident[v] == len(neighbours[v]) for v in range(len(coords))]
    rim = sum(1 for v in range(len(coords)) if not interior[v] and types[v] == FREE)
    no_patch = sum(
        1
        for face, flag in zip(faces, flags)
        if not flag and not all(interior[v] for v in face)
    )
    counts = {kind: types.count(kind) for kind in (FREE, FIXED, PERIODIC)}
    print(
        "%d vertices: %d free, %d fixed, %d periodic images (%d free vertices on an open rim)"
        % (len(coords), counts[FREE], counts[FIXED], counts[PERIODIC], rim)
    )
    print(
        "%d faces: %d are copies; of the %d that remain, %d carry a patch and %d touch an "
        "open rim and so carry none"
        % (len(faces), sum(flags), flags.count(0), flags.count(0) - no_patch, no_patch)
    )
    valences = sorted({len(neighbours[v]) for v in range(len(coords)) if interior[v]})
    print("interior valences: %s" % (valences,))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="mixed_mesh",
                        help="output prefix; writes <out>_vertices.csv and <out>_faces.csv")
    parser.add_argument("--side", type=float, default=120.0, help="box side in nm (both axes)")
    parser.add_argument("--side-y", type=float, default=None, help="y side, if it differs")
    parser.add_argument("--l-face", type=float, default=5.0, help="target edge length in nm")
    parser.add_argument("--boundary-x", choices=("periodic", "free", "fixed"), default="periodic")
    parser.add_argument("--boundary-y", choices=("periodic", "free", "fixed"), default="periodic")
    parser.add_argument("--fixed-rings", type=int, default=1,
                        help="rings of vertices a fixed side clamps")
    parser.add_argument("--pore-radius", type=float, default=0.0,
                        help="radius of the hole to carve, in nm; 0 for no hole")
    parser.add_argument("--neck-height", type=float, default=0.0,
                        help="height of the collar around the hole, in nm")
    parser.add_argument("--neck-width", type=float, default=25.0,
                        help="how far the collar reaches beyond the rim, in nm")
    args = parser.parse_args(argv)

    side_y = args.side if args.side_y is None else args.side_y
    n_face_x, n_face_y, d_face_x, d_face_y = axes_division(args.side, side_y, args.l_face)
    print("nFaceX = %d, nFaceY = %d, dFaceX = %.6g nm, dFaceY = %.6g nm"
          % (n_face_x, n_face_y, d_face_x, d_face_y))

    coords, faces, cell_of_face = build_lattice(n_face_x, n_face_y, d_face_x, d_face_y)
    orient_clockwise(coords, faces)
    types, mirrors = assign_boundary_types(n_face_x, n_face_y, args.boundary_x,
                                           args.boundary_y, max(1, args.fixed_rings))
    flags = copy_flags(n_face_x, n_face_y, args.boundary_x, args.boundary_y, cell_of_face)

    centre = (0.0, 0.0)
    coords, types, mirrors, faces, flags = carve_pore(
        coords, types, mirrors, faces, flags, args.pore_radius, centre)
    lift_neck(coords, types, mirrors, args.pore_radius, args.neck_height, args.neck_width, centre)

    report(coords, types, faces, flags)
    vertices_path, faces_path = write_mesh(args.out, coords, types, mirrors, faces, flags)
    print("wrote %s and %s" % (vertices_path, faces_path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
