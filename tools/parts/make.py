# -*- coding: utf-8 -*-
"""Write the parts library: six solids a scene can name, in millimetres.

An STL is read as millimetres by `pantometry-world`, and a part is geometry a `block` domain
voxelises. Every solid here is a **polygon extruded along z**, because that is the one primitive
whose area and volume have closed forms an arithmetic nobody has to trust:

    volume = shoelace(profile) * height
    area   = 2 * shoelace(profile) + perimeter(profile) * height

and for the one with a hole, the same with the inner ring subtracted from the first and added to
the second. `a_part_is_the_shape_it_claims_to_be` holds every file against those.

# The triangulation is the part a volume check cannot see

A fan from vertex zero is only valid when every vertex is visible from vertex zero. The L-bracket
happens to be; the heat sink is not, and a fan across its fins would put triangles in the air
between them. **The signed areas still telescope to the shoelace sum**, so the volume comes out
exactly right and the mesh is still closed -- neither of those checks can tell. The surface area
can, because triangles that overlap or leave the polygon do not add up to it. That is why the
test asserts both, and why this writes ears rather than a fan.

Run it from the repository root:

    python tools/parts/make.py
"""
import io
import math
import os

OUT = 'app/pantometry-world/scenes/parts'


# --------------------------------------------------------------------------- geometry


def shoelace(profile):
    """Twice the signed area of a closed polygon, halved. Positive when counter-clockwise."""
    total = 0.0
    n = len(profile)
    for i in range(n):
        ax, ay = profile[i]
        bx, by = profile[(i + 1) % n]
        total += ax * by - bx * ay
    return total / 2.0


def perimeter(profile):
    total = 0.0
    n = len(profile)
    for i in range(n):
        ax, ay = profile[i]
        bx, by = profile[(i + 1) % n]
        total += math.hypot(bx - ax, by - ay)
    return total


def cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def point_in_triangle(p, a, b, c):
    d1 = cross(a, b, p)
    d2 = cross(b, c, p)
    d3 = cross(c, a, p)
    neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    pos = (d1 > 0) or (d2 > 0) or (d3 > 0)
    return not (neg and pos)


def ears(profile):
    """Triangulate a simple counter-clockwise polygon by clipping ears.

    Returns index triples into `profile`. Ear clipping rather than a fan because a fan is wrong
    for a polygon no vertex sees all of, and wrong in a way the volume agrees with.
    """
    n = len(profile)
    idx = list(range(n))
    out = []
    guard = 0
    while len(idx) > 3:
        guard += 1
        if guard > 4 * n * n:
            raise AssertionError('no ear found; is the profile simple and counter-clockwise?')
        for k in range(len(idx)):
            i0 = idx[(k - 1) % len(idx)]
            i1 = idx[k]
            i2 = idx[(k + 1) % len(idx)]
            a, b, c = profile[i0], profile[i1], profile[i2]
            if cross(a, b, c) <= 0.0:
                continue  # reflex, or collinear
            blocked = False
            for j in idx:
                if j in (i0, i1, i2):
                    continue
                if point_in_triangle(profile[j], a, b, c):
                    blocked = True
                    break
            if blocked:
                continue
            out.append((i0, i1, i2))
            idx.pop(k)
            break
    out.append((idx[0], idx[1], idx[2]))
    return out


def extrude(profile, height):
    """A simple counter-clockwise polygon in z = 0, raised to `height`. Outward normals."""
    tris = []
    n = len(profile)
    for i in range(n):
        ax, ay = profile[i]
        bx, by = profile[(i + 1) % n]
        tris.append(((ax, ay, 0.0), (bx, by, 0.0), (bx, by, height)))
        tris.append(((ax, ay, 0.0), (bx, by, height), (ax, ay, height)))
    for i0, i1, i2 in ears(profile):
        a, b, c = profile[i0], profile[i1], profile[i2]
        tris.append(((a[0], a[1], 0.0), (c[0], c[1], 0.0), (b[0], b[1], 0.0)))
        tris.append(((a[0], a[1], height), (b[0], b[1], height), (c[0], c[1], height)))
    return tris


def extrude_ring(outer, inner, height):
    """An annulus: `outer` counter-clockwise, `inner` the hole, same vertex count.

    The caps are bands between corresponding vertices rather than fans, so no triangle crosses the
    hole -- the failure a fan of the outer ring alone would make, and which its volume would hide.
    """
    assert len(outer) == len(inner), 'the two rings are walked together'
    n = len(outer)
    tris = []
    for i in range(n):
        ax, ay = outer[i]
        bx, by = outer[(i + 1) % n]
        tris.append(((ax, ay, 0.0), (bx, by, 0.0), (bx, by, height)))
        tris.append(((ax, ay, 0.0), (bx, by, height), (ax, ay, height)))
    for i in range(n):
        # the hole's wall faces inward, so it is wound the other way
        ax, ay = inner[i]
        bx, by = inner[(i + 1) % n]
        tris.append(((ax, ay, 0.0), (bx, by, height), (bx, by, 0.0)))
        tris.append(((ax, ay, 0.0), (ax, ay, height), (bx, by, height)))
    for i in range(n):
        o0, o1 = outer[i], outer[(i + 1) % n]
        i0, i1 = inner[i], inner[(i + 1) % n]
        tris.append(((o0[0], o0[1], 0.0), (i0[0], i0[1], 0.0), (i1[0], i1[1], 0.0)))
        tris.append(((o0[0], o0[1], 0.0), (i1[0], i1[1], 0.0), (o1[0], o1[1], 0.0)))
        tris.append(((o0[0], o0[1], height), (i1[0], i1[1], height), (i0[0], i0[1], height)))
        tris.append(((o0[0], o0[1], height), (o1[0], o1[1], height), (i1[0], i1[1], height)))
    return tris


def to_origin(tris):
    """Move a solid so its lowest corner sits at the origin.

    **A `block` domain's grid starts at the origin.** `Voxels::onto` reads an STL's coordinates as
    absolute positions -- that is what lets an assembly of several files keep its relative
    placement -- and a mesh reaching outside the grid is *refused* rather than cropped, because a
    part with its corner missing runs and audits and answers about a different shape.

    The hand-made bracket was drawn in the positive octant and nothing said why. `rod` and `pipe`
    are built from `polygon()`, which is centred, and they shipped as two solids no scene could
    use: the volume, the area and the closed-edge count are all translation-invariant and none of
    them can see it. `a_dropped_file_becomes_a_domain` is what found it, by building a scene.
    """
    lo = [min(v[i] for t in tris for v in t) for i in range(3)]
    return [
        tuple(tuple(v[i] - lo[i] for i in range(3)) for v in t)
        for t in tris
    ]


def polygon(sides, radius, turn=0.0):
    """A regular polygon, counter-clockwise, circumradius `radius`."""
    return [
        (
            radius * math.cos(2.0 * math.pi * i / sides + turn),
            radius * math.sin(2.0 * math.pi * i / sides + turn),
        )
        for i in range(sides)
    ]


# --------------------------------------------------------------------------- writing


def number(x):
    """`%g`, which is what the hand-made bracket this library grew around was written with."""
    return '%g' % x


def write_stl(path, name, tris):
    lines = ['solid %s' % name]
    for a, b, c in tris:
        u = [b[i] - a[i] for i in range(3)]
        v = [c[i] - a[i] for i in range(3)]
        n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ]
        length = math.sqrt(sum(t * t for t in n))
        assert length > 0.0, '%s: a triangle with no area' % name
        n = [t / length for t in n]
        lines.append('  facet normal %s %s %s' % tuple(number(t) for t in n))
        lines.append('    outer loop')
        for p in (a, b, c):
            lines.append('      vertex %s %s %s' % tuple(number(t) for t in p))
        lines.append('    endloop')
        lines.append('  endfacet')
    lines.append('endsolid %s' % name)
    io.open(path, 'w', encoding='utf-8', newline='\n').write('\n'.join(lines) + '\n')
    return len(tris)


# --------------------------------------------------------------------------- the library

FINS = 5
FIN_W = 6.0
FIN_H = 16.0
GAP = 5.0
BASE_W = FINS * FIN_W + (FINS + 1) * GAP  # 60
BASE_H = 4.0


def heat_sink_profile():
    """A comb: a base with five fins. No vertex of it sees all the others."""
    p = [(0.0, 0.0), (BASE_W, 0.0), (BASE_W, BASE_H)]
    for k in range(FINS - 1, -1, -1):
        left = GAP + k * (FIN_W + GAP)
        right = left + FIN_W
        p.append((right, BASE_H))
        p.append((right, BASE_H + FIN_H))
        p.append((left, BASE_H + FIN_H))
        p.append((left, BASE_H))
    p.append((0.0, BASE_H))
    return p


ROD_SIDES = 24
ROD_R = 8.0
PIPE_SIDES = 24
PIPE_OUTER = 10.0
PIPE_INNER = 7.0

PARTS = [
    # name, what it is, triangles
    ('plate', 'a flat rectangular slab', extrude([(0.0, 0.0), (60.0, 0.0), (60.0, 40.0), (0.0, 40.0)], 5.0)),
    ('rod', 'a round bar, faceted', extrude(polygon(ROD_SIDES, ROD_R), 60.0)),
    ('wedge', 'a right triangular prism', extrude([(0.0, 0.0), (40.0, 0.0), (0.0, 20.0)], 30.0)),
    (
        'l-bracket',
        'two arms at a right angle, the inner corner chamfered',
        extrude(
            [
                (0.0, 0.0),
                (50.0, 0.0),
                (50.0, 20.0),
                (30.0, 20.0),
                (20.0, 30.0),
                (20.0, 50.0),
                (0.0, 50.0),
            ],
            20.0,
        ),
    ),
    ('heat-sink', 'a finned base', extrude(heat_sink_profile(), 30.0)),
    (
        'pipe',
        'a hollow round tube, faceted',
        extrude_ring(
            polygon(PIPE_SIDES, PIPE_OUTER),
            polygon(PIPE_SIDES, PIPE_INNER),
            50.0,
        ),
    ),
]


def main():
    assert os.path.isdir(OUT), 'run this from the repository root; %s is not there' % OUT
    for name, what, tris in PARTS:
        path = os.path.join(OUT, '%s.stl' % name)
        before = None
        if os.path.exists(path):
            before = io.open(path, 'rb').read()
        n = write_stl(path, name.replace('-', '_'), to_origin(tris))
        after = io.open(path, 'rb').read()
        state = 'new'
        if before is not None:
            state = 'unchanged' if before == after else 'REWRITTEN'
        print('  %-10s %4d triangles  %6d bytes  %s  -- %s' % (name, n, len(after), state, what))


if __name__ == '__main__':
    main()
