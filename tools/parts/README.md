# `tools/parts` — the geometry a scene can name

`make.py` writes six ASCII STL files to `app/pantometry-world/scenes/parts/`, and all six are
**committed**. A `block` domain names one and voxelises it:

```json
{ "kind": "block", "name": "bracket", "cells": [26, 26, 12], "cell_mm": 2.0,
  "parts": [ { "stl": "parts/l-bracket.stl", "material": "aluminium" } ] }
```

```sh
python tools/parts/make.py
```

It is deterministic: running it against an unchanged tree rewrites all six byte for byte.

## The six

Everything here is a **polygon extruded along `z`**, in millimetres, because that is the one
primitive whose two numbers follow from the profile rather than from the code that wrote it:

```
volume = cross_section * height
area   = 2 * cross_section + perimeter * height
```

| part | what it is | mm | triangles | volume mm³ | area mm² |
| --- | --- | --- | --- | --- | --- |
| `plate` | a flat slab | 60 × 40 × 5 | 12 | 12 000 | 5 800 |
| `rod` | a round bar, 24 facets | ⌀16 × 60 | 92 | 11 926.38 | 3 404.87 |
| `wedge` | a right triangular prism | 40 × 20 × 30 | 8 | 12 000 | 3 941.64 |
| `l-bracket` | two arms, inner corner chamfered | 50 × 50 × 20 | 24 | 33 000 | 7 182.84 |
| `heat-sink` | a base with five fins | 60 × 20 × 30 | 92 | 21 600 | 10 080 |
| `pipe` | a hollow tube, 24 facets | ⌀20/⌀14 × 50 | 192 | 7 919.86 | 5 642.26 |

59.9 KiB in total. `a_part_is_the_shape_it_claims_to_be` holds every one of them against the
arithmetic above, plus two things a volume alone cannot say: that every edge is shared by exactly
two triangles, and that the volume is **positive** rather than merely close in magnitude, since an
inside-out export is a real defect and the sign is the only check that sees it.

**The rod and the pipe are polygons and not circles**, so `πr²` is not their closed form and
`(n/2)r² sin(2π/n)` is. A 24-gon is **1.1384 %** under its circle, which is a per-cent effect and
not a rounding one. That shortfall is checked against the expansion of `sin x / x` — two terms give
1.1384 % and the third is `6.4e-8` — so the polygon formula is held against something that is not
the polygon formula.

The tolerance is `1e-5` relative and it comes from `%g`: six significant figures on a coordinate is
`5e-7`, a volume is a product of three of them, so `1.5e-6` is the most the file format can
contribute. The rod is measured `2.6e-7` off. Four of the six have integer coordinates and are
exact.

## A fan gets the volume exactly right and the area wrong

This is the reason `make.py` clips ears instead of fanning from vertex zero, and it is worth
stating because two of the three checks above **cannot tell**.

A fan is a triangulation only when every vertex is visible from vertex zero. The L-bracket happens
to qualify. A comb does not — no vertex of the heat sink sees all the others — so a fan lays
triangles in the air across the gaps between its fins. Measured, by generating exactly that:

```
              triangles   volume mm³   closed   area mm²
eared              92        21 600      yes      10 080
fanned             92        21 600      yes      14 400      <- 43 % high
```

The **signed** areas of a fan telescope to the shoelace sum whatever the polygon looks like, so the
volume is exact to the last digit. And a fan is still a valid *combinatorial* triangulation, so
every edge is shared by exactly two triangles and the mesh is closed. Only the surface area sees
it, because triangles that overlap or leave the polygon do not add up to it.

`a_fan_passes_the_volume_and_the_closed_check_and_fails_the_area_one` builds a non-star-shaped U
in the test and measures all three, so the claim is not a sentence in this file.

Four sabotages, each restored by re-running the generator:

| what was broken | which assertion caught it |
| --- | --- |
| caps fanned rather than eared | the area, and only the area |
| the profile wound the other way | the sign of the volume |
| four fins where the closed form says five | the volume |
| the rod 1 % longer | the volume |

## The L-bracket predates this and had no check at all

It was hand-made, it was the only part in the tree, and nothing in the repository tested it — it is
in the table above now like the rest. Regenerating it changed only the **cap triangulation**: the
walls are identical, both triangulations are valid, and the area comes out at the same
`7 182.8427 mm²`. The check that mattered was the one downstream — `pantometry run` on
`29-a-designed-bracket-becomes-cells.json` writes a **byte-identical** result before and after, so
the voxelisation does not depend on which of the two triangulations it is given.

The solid's internal name changed from `bracket` to `l_bracket`; `Mesh::from_stl` ignores it, as it
ignores the per-facet normals, and takes the winding.

## This is the second exception to "nothing generated is committed"

That rule is stated in `EXAMPLES.md` and `app/pantometry-world/scenes/README.md`, and it is about
what a *run* writes: the SVG an example plots, the glTF a scene exports. A stale one of those is
worse than none.

These are the other thing, and so are `tools/presets`'s: **inputs**. A scene that names
`parts/l-bracket.stl` needs the file to exist before anything runs, and a `block` domain is the
only place in the scene format where geometry comes from a file rather than from numbers — that
and `protein`, whose PDB files are committed for the same reason and are measurements rather than
output.

The guard against the usual price of a committed artefact — that it drifts from its source and
nobody notices — is that the files are checked against **closed forms** and not against the script
that wrote them. `make.py` could be deleted and the test would still say whether the six solids are
the shapes they claim to be.
