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
| `rod` | a round bar, 24 facets | ⌀16 × 60 | 92 | 11 926.40 | 3 404.87 |
| `wedge` | a right triangular prism | 40 × 20 × 30 | 8 | 12 000 | 3 941.64 |
| `l-bracket` | two arms, inner corner chamfered | 50 × 50 × 20 | 24 | 33 000 | 7 182.84 |
| `heat-sink` | a base with five fins | 60 × 20 × 30 | 92 | 21 600 | 10 080 |
| `pipe` | a hollow tube, 24 facets | ⌀20/⌀14 × 50 | 192 | 7 919.92 | 5 642.27 |

58.5 KiB in total. `a_part_is_the_shape_it_claims_to_be` holds every one of them against the
arithmetic above, plus two things a volume alone cannot say: that every edge is shared by exactly
two triangles, and that the volume is **positive** rather than merely close in magnitude, since an
inside-out export is a real defect and the sign is the only check that sees it.

## Every one of them starts at the origin, and two of them did not

A `block`'s grid starts at the origin of its parts' coordinates unless the scene says
`"grid_origin": "parts"`, and `Voxels::onto` reads an STL's coordinates as absolute positions --
which is what lets an assembly of several files keep its relative placement -- so without that key
a solid drawn about its own centre reaches outside its grid and is **refused rather than
cropped**. A dropped file is given the key, and `pantometry fit` writes it. The six here start at
the origin anyway, so a scene written by hand can name one without knowing the key exists.

The hand-made bracket happened to be drawn that way and nothing said why. `rod` and `pipe` are
built from `polygon()`, which is centred, and they **shipped as two solids no scene could use**.
Volume, area and the closed-edge count are all translation-invariant, so none of the three checks
above could see it; it took `a_dropped_file_becomes_a_domain` building a scene around one.
`to_origin` is applied to every part in `main` now, so a new one cannot forget, and the test
asserts the low corner.

**The rod and the pipe are polygons and not circles**, so `πr²` is not their closed form and
`(n/2)r² sin(2π/n)` is. A 24-gon is **1.1384 %** under its circle, which is a per-cent effect and
not a rounding one. That shortfall is checked against the expansion of `sin x / x` — two terms give
1.1384 % and the third is `6.4e-8` — so the polygon formula is held against something that is not
the polygon formula.

The tolerance is `1e-4` relative, and the first derivation of it was wrong. Six significant figures
is a half-ulp of `5e-6` when the mantissa is just above 1 and `5e-7` when it is just below 10, so
the worst per coordinate is **`5e-6`** and not the `5e-7` first written. A volume is a product of
three lengths, taking it to `1.5e-5`, and the **pipe's cross-section is a difference of two nearly
equal polygon areas** -- `310.58 - 152.19` -- which amplifies it by `2.92` to `2.9e-5`. The worst
measured is the pipe's volume at `6.8e-6`; four of the six have integer coordinates and are exact.

That the first bound was too tight showed up by accident: moving the parts to the origin changed
nothing about them but which coordinates `%g` was rounding, and the pipe went from `2.0e-7` to
`6.8e-6` -- most of the way to a bound that was supposed to have sevenfold headroom.

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
