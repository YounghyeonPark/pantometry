//! **A designed part, posed into the world, against the grid it will be rasterised onto.**
//!
//! `check` hands the shell two pictures of the same object: the wireframe of the block's extent,
//! and — now — the triangles of the parts inside it. If those two disagree the viewport is worse
//! than either alone, because a person looking at a part a millimetre outside its own box has no
//! way to know which of the two moved. `mesh`'s header says exactly this about two exporters
//! and it is the same failure from a different direction.
//!
//! So the checks here are about **placement**, not about triangles. That the triangles survive
//! the crossing is `pantometry-view`'s `a_designed_mesh.rs`; that they land in the right place is
//! this.
//!
//! No STL ships in this repository and no scene uses one, so the file is built here — which is
//! also the only way to know the answer exactly.

use editor_core::{check, Checked};
use pantometry_world::Parts;

/// The cube's side in **millimetres**, which is what an STL carries: `Mesh::from_stl` reads the
/// numbers as mm because every mechanical CAD tool writes mm, so 40 here is 0.040 m out.
const SIDE_MM: f32 = 40.0;

/// Where the part sits inside the block, in millimetres — off the block's origin on every axis,
/// so a transform that dropped the offset would be visible rather than cancel.
const AT_MM: f32 = 10.0;

/// A binary STL of an axis-aligned cube spanning `AT_MM .. AT_MM + SIDE_MM`.
///
/// 80-byte header, a `u32` count, then 50 bytes a triangle: a normal and three vertices as
/// `f32`, and a `u16` nobody reads. The normals are written correctly and then ignored —
/// `from_stl` takes the winding as the truth, which is what every reader that has been burned by
/// a disagreeing normal does.
fn cube_stl() -> Vec<u8> {
    let v = |i: usize| {
        [
            AT_MM + if i & 1 == 0 { 0.0 } else { SIDE_MM },
            AT_MM + if i & 2 == 0 { 0.0 } else { SIDE_MM },
            AT_MM + if i & 4 == 0 { 0.0 } else { SIDE_MM },
        ]
    };
    // Six outward quads, each split into two triangles. Corner bits are x, y, z.
    let quads = [
        [0, 2, 3, 1], // z low, outward -z
        [4, 5, 7, 6], // z high, outward +z
        [0, 1, 5, 4], // y low, outward -y
        [2, 6, 7, 3], // y high, outward +y
        [0, 4, 6, 2], // x low, outward -x
        [1, 3, 7, 5], // x high, outward +x
    ];
    let mut tris: Vec<[[f32; 3]; 3]> = Vec::new();
    for q in quads {
        tris.push([v(q[0]), v(q[1]), v(q[2])]);
        tris.push([v(q[0]), v(q[2]), v(q[3])]);
    }

    let mut out = vec![0u8; 80];
    out.extend((tris.len() as u32).to_le_bytes());
    for t in &tris {
        let u = [t[1][0] - t[0][0], t[1][1] - t[0][1], t[1][2] - t[0][2]];
        let w = [t[2][0] - t[0][0], t[2][1] - t[0][1], t[2][2] - t[0][2]];
        let n = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        for c in n {
            out.extend((c / len).to_le_bytes());
        }
        for p in t {
            for c in p {
                out.extend(c.to_le_bytes());
            }
        }
        out.extend(0u16.to_le_bytes());
    }
    out
}

/// One STL under one name, so the scene can reach it without a filesystem.
struct InMemory;

impl Parts for InMemory {
    fn bytes(&self, name: &str) -> Result<Vec<u8>, String> {
        if name == "cube.stl" {
            Ok(cube_stl())
        } else {
            Err(format!("{name}: no such part"))
        }
    }
}

/// A block big enough to hold the cube at its offset: 60 mm on a side, in 2 mm cells.
///
/// `place` goes in verbatim at the **scene** level, which is where a pose lives: `Voxels::onto`
/// rasterises in the domain's own frame and `Pose` moves the domain, so a pose changes where the
/// part is drawn without changing a cell of what is rasterised.
fn scene(place: &str) -> String {
    format!(
        r#"{{
  "title": "one designed part",
  "duration_s": 0.1,
  "frames": 2,
  "conservation_tolerance": 1e-9,{place}
  "domains": [
    {{ "kind": "block", "name": "assembly", "cells": [30, 30, 30], "cell_mm": 2.0,
       "initial_c": 20.0,
       "parts": [ {{ "stl": "cube.stl", "material": "aluminium" }} ] }}
  ]
}}"#
    )
}

fn checked(place: &str) -> Checked {
    let c = check(&scene(place), &InMemory);
    assert!(c.error.is_none(), "the scene did not build: {:?}", c.error);
    c
}

/// The tightest box round every triangle, `[x0, y0, z0, x1, y1, z1]`, in metres.
fn bounds_of(m: &editor_core::PlacedMesh) -> [f64; 6] {
    let mut b = [f64::INFINITY; 3]
        .into_iter()
        .chain([f64::NEG_INFINITY; 3])
        .collect::<Vec<_>>();
    for t in &m.triangles {
        for c in t {
            for a in 0..3 {
                b[a] = b[a].min(c[a]);
                b[a + 3] = b[a + 3].max(c[a]);
            }
        }
    }
    [b[0], b[1], b[2], b[3], b[4], b[5]]
}

/// Enclosed volume by the divergence theorem, from the placed triangles.
fn volume_of(m: &editor_core::PlacedMesh) -> f64 {
    m.triangles
        .iter()
        .map(|t| {
            let cross = [
                t[1][1] * t[2][2] - t[1][2] * t[2][1],
                t[1][2] * t[2][0] - t[1][0] * t[2][2],
                t[1][0] * t[2][1] - t[1][1] * t[2][0],
            ];
            (t[0][0] * cross[0] + t[0][1] * cross[1] + t[0][2] * cross[2]) / 6.0
        })
        .sum()
}

/// The STL's coordinates are `f32` and the arithmetic is `f64`, so the floor is `f32`'s mantissa
/// on a 0.05 m coordinate — about `3e-9` m. A metre-scale tolerance of `1e-9` would be under it;
/// `1e-7` m is a tenth of a micron, far below any placement error worth drawing and far above
/// the representation.
const NEAR_M: f64 = 1e-7;

#[test]
fn a_part_arrives_at_all_and_says_where_it_came_from() {
    // The silent-failure shape first: a scene with a part and no meshes is an empty viewport
    // that looks exactly like a scene with no parts.
    let c = checked("");
    assert_eq!(c.meshes.len(), 1, "the part did not arrive");
    let m = &c.meshes[0];
    assert_eq!(m.name, "assembly", "hiding the domain must hide the part");
    assert_eq!(
        m.site, "assembly/parts[0]",
        "the site the build's notes and verify's findings already use"
    );
    assert_eq!(m.stl, "cube.stl");
    assert_eq!(m.triangles.len(), 12, "twelve triangles of a cube");
}

#[test]
fn an_unposed_part_is_where_the_file_put_it() {
    // No pose, so world *is* the block's local frame, and the answer is the file read as
    // millimetres: a 40 mm cube at 10 mm becomes 0.010 .. 0.050 m on every axis. This is the
    // check that would have caught the thousandfold error the format's own doc comment warns
    // about, and the doubling that once made a 40 mm cube 80 mm across.
    let c = checked("");
    let b = bounds_of(&c.meshes[0]);
    let (lo, hi) = (
        f64::from(AT_MM) / 1000.0,
        f64::from(AT_MM + SIDE_MM) / 1000.0,
    );
    for a in 0..3 {
        assert!((b[a] - lo).abs() < NEAR_M, "axis {a} starts at {}", b[a]);
        assert!(
            (b[a + 3] - hi).abs() < NEAR_M,
            "axis {a} ends at {}",
            b[a + 3]
        );
    }
}

#[test]
fn a_part_is_inside_the_grid_it_will_be_rasterised_onto() {
    // The claim the whole feature rests on. `Voxels::onto` rasterises against a grid whose
    // origin is the block's local zero and whose span is `cells × cell_mm` — 30 × 2 mm = 60 mm —
    // and it *refuses* a mesh that does not fit. So a part drawn outside that span is a part
    // drawn somewhere its voxels can never be, and the two pictures would disagree.
    //
    // Checked against the grid's arithmetic rather than against `Checked::boxes`, deliberately:
    // the box comes from the same `placement` this transform uses, so comparing the two would
    // be a transform agreeing with itself.
    let c = checked("");
    let span = 30.0 * 2.0 / 1000.0;
    let b = bounds_of(&c.meshes[0]);
    for a in 0..3 {
        assert!(b[a] >= -NEAR_M, "axis {a} starts before the grid: {}", b[a]);
        assert!(
            b[a + 3] <= span + NEAR_M,
            "axis {a} ends past the grid: {} > {span}",
            b[a + 3]
        );
    }
}

#[test]
fn a_pose_moves_a_part_and_does_not_resize_it() {
    // A `Pose` is a rigid motion — rotation and translation, no scale and no shear — so the
    // volume it encloses is an invariant. That is the closed form, and it is what separates
    // "placed" from "placed and quietly scaled": a transform applied twice, or applied in
    // millimetres to metres, changes this number and cannot change it back.
    let here = checked("");
    let there = checked(
        r#"
  "poses": { "assembly": { "at_m": [0.100, 0.0, 0.0],
                           "turn": { "axis": [0.0, 0.0, 1.0], "degrees": 90.0 } } },"#,
    );

    let side = f64::from(SIDE_MM) / 1000.0;
    let exact = side * side * side;
    for (what, m) in [("unposed", &here.meshes[0]), ("posed", &there.meshes[0])] {
        let v = volume_of(m);
        assert!(
            (v - exact).abs() / exact < 1e-6,
            "{what} encloses {v}, not {exact}"
        );
    }

    // And it moved *to a stated place*. The volume above is invariant under a pose that was
    // dropped on the floor, which is the likelier bug of the two, so this has to say where.
    //
    // A quarter turn about `+z` is closed form and needs no trigonometry: right-handed,
    // `(x, y) -> (-y, x)`. The cube spans 0.010 .. 0.050 on every axis, so afterwards
    // `x ∈ -0.050 .. -0.010` and `y ∈ 0.010 .. 0.050`, and the 0.100 m translation puts x at
    // 0.050 .. 0.090. z is untouched.
    //
    // Written out rather than compared to a threshold: the first draft of this asserted the
    // part had moved "more than 0.05 m", measured 0.040 because the turn cancels part of the
    // translation, and the fix for *that* is an exact answer rather than a smaller number.
    let want = [0.050, 0.010, 0.010, 0.090, 0.050, 0.050];
    let got = bounds_of(&there.meshes[0]);
    for a in 0..6 {
        assert!(
            (got[a] - want[a]).abs() < NEAR_M,
            "bound {a} is {} and should be {}",
            got[a],
            want[a]
        );
    }
}

#[test]
fn a_part_that_cannot_be_read_leaves_the_error_and_not_a_phantom() {
    // An unreadable STL must not produce a *partial* mesh — half a part on screen is worse than
    // none, because it reads as geometry rather than as a failure. The build says why, and this
    // pins that `check` does not draw anything anyway.
    let missing = scene("").replace("cube.stl", "absent.stl");
    let c = check(&missing, &InMemory);
    assert!(c.meshes.is_empty(), "a phantom part was drawn");
    let e = c.error.expect("an unreadable part must fail the build");
    assert!(e.contains("absent.stl"), "the error does not name it: {e}");
}

#[test]
fn a_scene_with_no_parts_has_no_meshes() {
    // The two-way half: `meshes` must be empty rather than defaulted to something, or the
    // viewport draws a part into every scene that never had one.
    let c = check(
        r#"{ "title": "no parts", "duration_s": 0.1, "frames": 2,
             "conservation_tolerance": 1e-9,
             "domains": [ { "kind": "block", "name": "plain", "cells": [4, 4, 4],
                            "cell_mm": 2.0, "initial_c": 20.0, "material": "aluminium" } ] }"#,
        &InMemory,
    );
    assert!(c.error.is_none(), "{:?}", c.error);
    assert!(c.meshes.is_empty());
    assert!(!c.boxes.is_empty(), "the block still has a box");
}
