//! **A click on a face is the face, the cell under it is the cell, and what is written is what the
//! format reads.**
//!
//! [`editor_core::pick`] is the arithmetic behind pointing at a model in the viewport: a ray against
//! a domain's box in its own frame, the cell under the point, the area a cooled face should state,
//! and the edits. Each is held to a closed form — a ray along an axis enters through the face it
//! points at, at a point and a distance that are subtraction; a turn by a known angle sends a known
//! face to a known side — and the edits to the format's own build.
//!
//! # Not under `wasm32`
//!
//! The scenes are read off a disk, and a `wasm32` target has none.

#![cfg(not(target_family = "wasm"))]

use editor_core::pick::{self, Stated, Target};
use pantometry::scene::Placed;
use pantometry_world::FaceSpec;
use viewer_core::Ray;

/// The shipped scenes.
fn scene(name: &str) -> String {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../pantometry-world/scenes");
    let path = std::fs::read_dir(&dir)
        .expect("the scenes directory reads")
        .map(|e| e.expect("an entry").path())
        .find(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with(name))
        })
        .unwrap_or_else(|| panic!("no scene {name} in {}", dir.display()));
    std::fs::read_to_string(path).expect("the scene reads")
}

fn files() -> editor_core::Beside {
    editor_core::Beside(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../pantometry-world/scenes"),
    )
}

fn checked(text: &str) -> editor_core::Checked {
    editor_core::check(text, &files())
}

fn target(text: &str, name: &str) -> Target {
    checked(text)
        .targets
        .into_iter()
        .find(|t| t.name == name)
        .unwrap_or_else(|| panic!("no target {name}"))
}

/// A unit ray.
fn ray(origin: [f64; 3], toward: [f64; 3]) -> Ray {
    let l = (toward[0] * toward[0] + toward[1] * toward[1] + toward[2] * toward[2]).sqrt();
    Ray {
        origin,
        direction: [toward[0] / l, toward[1] / l, toward[2] / l],
    }
}

/// A 20 x 30 x 40 mm box, its low corner at (1, 2, 3) mm.
const BOX: [f64; 6] = [0.001, 0.002, 0.003, 0.021, 0.032, 0.043];

/// **Along each axis, a ray enters through the face it points at**, at the point and the distance
/// subtraction gives: from 0.1 m out along `-x`, the face is `x-max`, at `(0.021, y, z)`, `0.1 -
/// 0.021` along. Six faces, three points on each, two of them on an edge and a corner where the
/// face is held to the box. Exact to f64 rounding, which is the floor `1e-15` of the box's size.
#[test]
fn along_each_axis_a_ray_enters_through_the_face_it_points_at() {
    let here = Placed::HERE;
    for face in FaceSpec::ALL {
        let (axis, high) = face.normal();
        let [a, b] = face.axes();
        for (fa, fb) in [(0.5, 0.5), (0.25, 0.9), (1.0, 0.0)] {
            let mut origin = [0.0; 3];
            origin[a] = BOX[a] + fa * (BOX[a + 3] - BOX[a]);
            origin[b] = BOX[b] + fb * (BOX[b + 3] - BOX[b]);
            origin[axis] = if high { 0.1 } else { -0.1 };
            let mut toward = [0.0; 3];
            toward[axis] = if high { -1.0 } else { 1.0 };
            let got = pick::hit(BOX, here, &ray(origin, toward))
                .unwrap_or_else(|| panic!("{face:?} at ({fa}, {fb}) was missed"));
            assert_eq!(got.face, face, "entered at ({fa}, {fb})");
            let plane = if high { BOX[axis + 3] } else { BOX[axis] };
            let mut want = origin;
            want[axis] = plane;
            for k in 0..3 {
                assert!(
                    (got.local[k] - want[k]).abs() <= 1e-15,
                    "{face:?}: {:?} against {want:?}",
                    got.local
                );
            }
            let distance = (origin[axis] - plane).abs();
            assert!(
                (got.distance - distance).abs() <= 1e-15,
                "{face:?}: {} m along, and subtraction says {distance}",
                got.distance
            );
        }
    }
}

/// **An oblique ray enters through whichever face its line crosses first.** From `(-0.05, 0.017,
/// 0.023)` toward `+x` and `+z` at one in ten: the line reaches `x = 0.001` at `t ≈ 0.0513`, where
/// `z` is `0.0281` — inside `[0.003, 0.043]` — so it comes in through `x-min`. Turned steeper, one
/// in one, it is at `z = 0.074` by then, past the top: it would have to be inside already, so it
/// misses. Both read off the line's equation.
#[test]
fn an_oblique_ray_enters_where_its_line_says() {
    let origin = [-0.05, 0.017, 0.023];
    let shallow = pick::hit(BOX, Placed::HERE, &ray(origin, [1.0, 0.0, 0.1])).expect("a hit");
    assert_eq!(shallow.face, FaceSpec::XMin);
    let t = (0.001 - origin[0]) / 1.0;
    let z = origin[2] + 0.1 * t;
    assert!(
        (shallow.local[2] - z).abs() < 1e-15,
        "z {} against {z}",
        shallow.local[2]
    );
    assert!(pick::hit(BOX, Placed::HERE, &ray(origin, [1.0, 0.0, 1.0])).is_none());
}

/// **What misses is `None`**: a ray beside the box, one pointing away from it, one parallel to a
/// face outside it, and one that is not a number. And an eye inside the box sees the face in front
/// of it — the one it would leave through.
#[test]
fn what_misses_is_none_and_an_eye_inside_sees_ahead() {
    let here = Placed::HERE;
    for (what, r) in [
        ("beside", ray([0.05, 0.017, 0.1], [0.0, 0.0, -1.0])),
        ("away", ray([0.011, 0.017, 0.1], [0.0, 0.0, 1.0])),
        ("parallel", ray([0.011, 0.017, 0.05], [1.0, 0.0, 0.0])),
        (
            "not a number",
            Ray {
                origin: [f64::NAN, 0.0, 0.0],
                direction: [0.0, 0.0, -1.0],
            },
        ),
    ] {
        assert!(pick::hit(BOX, here, &r).is_none(), "{what} was hit");
    }
    let inside = pick::hit(BOX, here, &ray([0.011, 0.017, 0.023], [0.0, 1.0, 0.0])).expect("a hit");
    assert_eq!(inside.face, FaceSpec::YMax);
    assert!((inside.distance - 0.015).abs() < 1e-15);
}

/// A unit quaternion turning by `degrees` about `axis`, `[x, y, z, w]`.
fn turn(axis: [f64; 3], degrees: f64) -> [f64; 4] {
    let l = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    let h = degrees.to_radians() / 2.0;
    let s = h.sin() / l;
    [axis[0] * s, axis[1] * s, axis[2] * s, h.cos()]
}

/// **A turned domain is hit on its own face.** A quarter turn about `z` sends the box's `-y` to the
/// world's `+x`, so an eye on `+x` looking back sees `y-min`; a third of a turn about `(1, 1, 1)`
/// sends `x` to `y`, `y` to `z` and `z` to `x`, so an eye above sees `y-max`. Both read off what a
/// rotation is. And the point found, placed back in the world, is on the ray.
#[test]
fn a_turned_domain_is_hit_on_its_own_face() {
    let centre = |p: Placed| {
        p.apply([
            0.5 * (BOX[0] + BOX[3]),
            0.5 * (BOX[1] + BOX[4]),
            0.5 * (BOX[2] + BOX[5]),
        ])
    };
    for (place, from, want) in [
        (
            Placed {
                at_m: [0.3, -0.2, 0.05],
                turn: turn([0.0, 0.0, 1.0], 90.0),
            },
            [1.0, 0.0, 0.0],
            FaceSpec::YMin,
        ),
        (
            Placed {
                at_m: [-0.04, 0.11, 0.7],
                turn: turn([1.0, 1.0, 1.0], 120.0),
            },
            [0.0, 0.0, 1.0],
            FaceSpec::YMax,
        ),
    ] {
        let c = centre(place);
        let origin = [c[0] + from[0], c[1] + from[1], c[2] + from[2]];
        let r = ray(origin, [-from[0], -from[1], -from[2]]);
        let got = pick::hit(BOX, place, &r).expect("the turned box is hit");
        assert_eq!(got.face, want, "turned {:?}", place.turn);
        let back = place.apply(got.local);
        for k in 0..3 {
            let on = r.origin[k] + got.distance * r.direction[k];
            assert!(
                (back[k] - on).abs() < 1e-14,
                "the point found is {back:?} and the ray is at {on} along {k}"
            );
        }
    }
}

/// **The cell under a point is the one whose span holds it**, and a point on the far wall is in
/// the last cell rather than one past the grid. Ten cells from 0 to 1: `0.3` is the start of cell
/// 3, a hair under it is cell 2, `1.0` is cell 9, and outside is held to the grid.
#[test]
fn the_cell_under_a_point_is_the_one_that_holds_it() {
    for (x, want) in [
        (0.0, 0),
        (0.05, 0),
        (0.3, 3),
        (0.3 - 1e-12, 2),
        (0.99, 9),
        (1.0, 9),
        (-0.4, 0),
        (7.0, 9),
    ] {
        assert_eq!(pick::cell_along(x, 0.0, 1.0, 10), want, "x = {x}");
    }
    assert_eq!(pick::patch_between([5, 2], [3, 7]), ([3, 2], [6, 8]));
}

/// **A cooled face states the area it covers, and on a plain block that is the format's own new
/// entry to the bit.** Scene 21's wax is 11 x 11 x 11 cells of 2 mm with no void: every face is 121
/// cells of 4 mm², 4.84 cm², which is what `templates::entry_for` writes for it; a 3 x 4 patch is
/// 0.48.
///
/// On a designed part the area is its solid cells. Scene 29's bracket is rasterised from an STL,
/// and `Solid3D::cells_on`'s own account measured 410 of the 676 cells on its `z-min` solid; read
/// here through the void map the build hands back, it is the same 410, and the area is 16.4 cm².
/// The order the map is read in and a high face's layer are held by the housing and the busbars
/// below, whose voids are not symmetric the way the bracket's L is.
#[test]
fn a_cooled_face_states_the_area_it_covers() {
    let wax = target(&scene("21-"), "wax");
    // Nothing in the wax's file can make a cell void, so its map is all solid however it was got.
    assert!(
        wax.void.as_ref().is_some_and(|v| v.iter().all(|x| !*x)),
        "the wax has a void cell, or no map"
    );
    for face in FaceSpec::ALL {
        let (area, solid) = wax.area_cm2(face, None);
        assert_eq!(solid, 121, "{face:?}");
        let fit = pantometry_world::templates::Fit {
            grid: Some((wax.cells, wax.cell_mm)),
            faces: &[],
            material: None,
        };
        let entry = pantometry_world::templates::entry_for("block", "cooling", &fit).expect("one");
        let stated: serde_json::Value = serde_json::from_str(&entry).expect("parses");
        assert_eq!(
            area.to_bits(),
            stated["area_cm2"].as_f64().expect("a number").to_bits(),
            "{face:?}: {area} against the new entry's {}",
            stated["area_cm2"]
        );
    }
    let (patch, solid) = wax.area_cm2(FaceSpec::ZMax, Some(([2, 3], [5, 7])));
    assert_eq!(solid, 12);
    assert!((patch - 0.48).abs() < 1e-15, "{patch}");

    let bracket = target(&scene("29-"), "bracket");
    let void = bracket.void.as_ref().expect("the bracket built");
    assert!(
        void.iter().any(|v| *v),
        "a designed part has void around it"
    );
    let (area, solid) = bracket.area_cm2(FaceSpec::ZMin, None);
    assert_eq!(solid, 410, "solid cells on the bracket's z-min");
    assert!((area - 16.4).abs() < 1e-12, "{area}");
}

/// **A hollow housing counts the walls a face has, on each of its six faces.** Scene 23's housing
/// is 6 x 6 x 7 cells of 8 mm, and its regions make layers `k = 0..3` a ring one cell thick round a
/// hot part and layers `3..5` empty, so — read off the scene's own regions by hand:
///
/// - each side face keeps only the lid's two layers, `k = 5, 6`: 6 x 2 = 12 cells;
/// - `z-min` is the part's 4 x 4 footprint with the ring's 20 cells gone: 16;
/// - `z-max`, the lid, is whole: 36.
///
/// Not a cube and not symmetric in `k`, so a void map read in the wrong order, or a high face read
/// one layer past the grid, gives other numbers.
#[test]
fn a_hollow_housing_counts_the_walls_each_face_has() {
    let housing = target(&scene("23-"), "housing");
    assert!(housing.void.is_some(), "the housing did not build");
    for (face, want) in [
        (FaceSpec::XMin, 12),
        (FaceSpec::XMax, 12),
        (FaceSpec::YMin, 12),
        (FaceSpec::YMax, 12),
        (FaceSpec::ZMin, 16),
        (FaceSpec::ZMax, 36),
    ] {
        let (area, solid) = housing.area_cm2(face, None);
        assert_eq!(solid, want, "{face:?}");
        assert!(
            (area - want as f64 * 0.64).abs() < 1e-12,
            "{face:?}: {area} cm²"
        );
    }
}

/// **The highlight is on the box that is drawn.** The whole of each target, as a box of cells, has
/// the corners the viewport's wireframe is drawn from, and each face's patch corners are four of
/// them — on scene 30, whose busbars are moved and turned thirty degrees.
#[test]
fn the_highlight_is_on_the_box_that_is_drawn() {
    let text = scene("30-");
    let c = checked(&text);
    let busbars = c
        .targets
        .iter()
        .find(|t| t.name == "busbars")
        .expect("busbars");
    assert!(!busbars.place.is_here(), "scene 30 is meant to be posed");
    let drawn = c.boxes.iter().find(|b| b.name == "busbars").expect("a box");
    let corners = busbars.cells_corners([0; 3], busbars.cells);
    let near = |p: [f64; 3], q: [f64; 3]| (0..3).all(|k| (p[k] - q[k]).abs() < 1e-15);
    for (n, p) in corners.iter().enumerate() {
        assert!(
            near(*p, drawn.corners[n]),
            "corner {n}: {p:?} against {:?}",
            drawn.corners[n]
        );
    }
    for face in FaceSpec::ALL {
        for p in busbars.patch_corners(face, busbars.whole(face)) {
            assert!(
                drawn.corners.iter().any(|q| near(p, *q)),
                "{face:?}'s corner {p:?} is no corner of the drawn box"
            );
        }
    }
    // **And inside it, in order round the patch.** Cells `[2, 1]..[5, 3]` on `x-max` are `j` from
    // 8 to 20 mm and `k` from 4 to 12 mm at `x` = 32 mm, the busbars' cells being 4 mm: the
    // corners are those four points, placed, in the order (j0, k0), (j1, k0), (j1, k1), (j0, k1).
    let want = [(8.0, 4.0), (20.0, 4.0), (20.0, 12.0), (8.0, 12.0)]
        .map(|(j, k)| busbars.place.apply([0.032, j * 1e-3, k * 1e-3]));
    let got = busbars.patch_corners(FaceSpec::XMax, ([2, 1], [5, 3]));
    for n in 0..4 {
        assert!(
            near(got[n], want[n]),
            "corner {n}: {:?} against {:?}",
            got[n],
            want[n]
        );
    }
}

/// The text of one entry, `/domains/<d>/<key>/<j>`, as JSON.
fn entry(text: &str, domain: &str, key: &str, j: usize) -> serde_json::Value {
    let root: serde_json::Value = serde_json::from_str(text).expect("parses");
    root["domains"]
        .as_array()
        .expect("domains")
        .iter()
        .find(|d| d["name"] == domain)
        .expect("the domain")[key][j]
        .clone()
}

/// **Pointing a cooled face somewhere else writes the face and its area, and pointing again at what
/// is set writes nothing.** Scene 21's wax has no cooling, so it is given the format's new entry on
/// its top and then pointed at its `x-max`: the face changes, the area is that face's, and the scene
/// still builds. A patch adds `from` and `to` on the entry's own line; the whole face again takes
/// them out, and the file is the one it was, byte for byte.
#[test]
fn a_cooled_face_is_pointed_at_and_back() {
    let start = editor_core::add_entry(&scene("21-"), "wax", "cooling").expect("added");
    let wax = target(&start, "wax");
    let (text, said) = pick::set_face(&start, &wax, "cooling", 0, FaceSpec::XMax, None)
        .expect("written")
        .expect("a change");
    println!("  {said}");
    let e = entry(&text, "wax", "cooling", 0);
    assert_eq!(e["face"], "x-max");
    assert_eq!(
        e["area_cm2"].as_f64(),
        Some(wax.area_cm2(FaceSpec::XMax, None).0)
    );
    assert!(checked(&text).error.is_none(), "{:?}", checked(&text).error);
    assert_eq!(
        pick::set_face(&text, &wax, "cooling", 0, FaceSpec::XMax, None),
        Ok(None),
        "pointing at the face that is set changed the file"
    );

    let patch = ([2, 3], [5, 7]);
    let (patched, _) = pick::set_face(&text, &wax, "cooling", 0, FaceSpec::XMax, Some(patch))
        .expect("written")
        .expect("a change");
    let e = entry(&patched, "wax", "cooling", 0);
    assert_eq!(e["from"], serde_json::json!([2, 3]));
    assert_eq!(e["to"], serde_json::json!([5, 7]));
    assert_eq!(e["area_cm2"].as_f64(), Some(0.48));
    assert!(
        checked(&patched).error.is_none(),
        "{:?}",
        checked(&patched).error
    );
    assert_eq!(
        pick::stated(&patched, "/domains/wax/cooling/0"),
        Some(Stated::Face {
            face: FaceSpec::XMax,
            patch: Some(patch)
        })
    );
    let (whole, _) = pick::set_face(&patched, &wax, "cooling", 0, FaceSpec::XMax, None)
        .expect("written")
        .expect("a change");
    assert_eq!(
        whole, text,
        "taking the patch off did not give the file back"
    );
    // A patch that is the whole face is the whole face.
    assert_eq!(
        pick::set_face(
            &text,
            &wax,
            "cooling",
            0,
            FaceSpec::XMax,
            Some(wax.whole(FaceSpec::XMax))
        ),
        Ok(None)
    );
}

/// **A patch written over a multi-line entry, and taken off again, gives the file back** — scene
/// 29's pad under the bracket, which states its patch on a line of its own. And a patch of only
/// void is refused with the file left alone, because it would cool nothing.
#[test]
fn a_patch_on_a_part_counts_its_solid_cells() {
    let text = scene("29-");
    let bracket = target(&text, "bracket");
    let (moved, said) = pick::set_face(&text, &bracket, "cooling", 0, FaceSpec::ZMin, None)
        .expect("written")
        .expect("a change");
    println!("  {said}");
    let e = entry(&moved, "bracket", "cooling", 0);
    assert!(e.get("from").is_none() && e.get("to").is_none(), "{e}");
    assert!((e["area_cm2"].as_f64().expect("a number") - 16.4).abs() < 1e-12);
    let (back, _) = pick::set_face(
        &moved,
        &bracket,
        "cooling",
        0,
        FaceSpec::ZMin,
        Some(([2, 21], [8, 25])),
    )
    .expect("written")
    .expect("a change");
    assert_eq!(back, text, "the pad put back is not the file that had it");

    // A cell on z-min with no solid under it.
    let empty = (0..bracket.cells[0])
        .flat_map(|i| (0..bracket.cells[1]).map(move |j| (i, j)))
        .find(|&(i, j)| !bracket.is_solid([i, j, 0]))
        .expect("the bracket leaves some of z-min empty");
    let void = ([empty.0, empty.1], [empty.0 + 1, empty.1 + 1]);
    let refused = pick::set_face(&text, &bracket, "cooling", 0, FaceSpec::ZMin, Some(void))
        .expect_err("a patch of void was written");
    assert!(refused.contains("no solid cell"), "{refused}");
}

/// **A held face is a whole face**: a box drawn on one is dropped and said, and what is written is
/// the face alone, which the format reads.
#[test]
fn a_held_face_takes_no_patch() {
    let text = scene("25-");
    let stress = target(&text, "stress");
    let (out, said) = pick::set_face(
        &text,
        &stress,
        "held",
        0,
        FaceSpec::YMax,
        Some(([1, 1], [3, 3])),
    )
    .expect("written")
    .expect("a change");
    assert!(said.contains("whole face"), "{said}");
    let e = entry(&out, "stress", "held", 0);
    assert_eq!(e["face"], "y-max");
    assert!(e.get("from").is_none(), "{e}");
    assert!(serde_json::from_str::<pantometry_world::Scene>(&out).is_ok());
}

/// **A box drawn on a face sets the box in that face's two axes and keeps its depth.** Scene 29's
/// heat source is cells `[20, 2, 0]..[25, 8, 12]`; drawn on the top over `[1, 1]..[4, 4]` it is
/// `[1, 1, 0]..[4, 4, 12]`, and drawn again it is unchanged.
#[test]
fn a_box_drawn_on_a_face_keeps_its_depth() {
    let text = scene("29-");
    let bracket = target(&text, "bracket");
    let (out, said) = pick::set_cells(
        &text,
        &bracket,
        "dissipation",
        0,
        FaceSpec::ZMax,
        ([1, 1], [4, 4]),
    )
    .expect("written")
    .expect("a change");
    println!("  {said}");
    assert_eq!(
        pick::stated(&out, "/domains/bracket/dissipation/0"),
        Some(Stated::Cells {
            from: [1, 1, 0],
            to: [4, 4, 12]
        })
    );
    assert_eq!(
        pick::set_cells(
            &out,
            &bracket,
            "dissipation",
            0,
            FaceSpec::ZMax,
            ([1, 1], [4, 4])
        ),
        Ok(None)
    );
    // On a side, the side's axes: x-max is j and k.
    let (side, _) = pick::set_cells(
        &text,
        &bracket,
        "dissipation",
        0,
        FaceSpec::XMax,
        ([3, 1], [6, 2]),
    )
    .expect("written")
    .expect("a change");
    assert_eq!(
        pick::stated(&side, "/domains/bracket/dissipation/0"),
        Some(Stated::Cells {
            from: [20, 3, 1],
            to: [25, 6, 2]
        })
    );
}

/// **A tree path names an entry, and the format says how it is placed.** Only an entry of a list
/// that says where is anything to point at.
#[test]
fn a_path_names_an_entry_and_the_format_says_how() {
    use pantometry_world::templates::Placing;
    assert_eq!(
        pick::entry_of("/domains/wax/cooling/3"),
        Some(("wax".into(), "cooling".into(), 3))
    );
    for path in [
        "/domains/wax",
        "/domains/wax/cooling",
        "/run/wax",
        "/domains/wax/cooling/x",
    ] {
        assert_eq!(pick::entry_of(path), None, "{path}");
    }
    let text = scene("25-");
    assert_eq!(
        pick::placing_of(&text, "/domains/module/cooling/0"),
        Some(Placing::FacePatch)
    );
    assert_eq!(
        pick::placing_of(&text, "/domains/stress/held/0"),
        Some(Placing::Face)
    );
    assert_eq!(
        pick::placing_of(&text, "/domains/module/dissipation/0"),
        Some(Placing::Cells)
    );
    assert_eq!(pick::placing_of(&text, "/domains/module/contact/0"), None);
}

/// **The cell under a point and its place on a face are what division says.** Scene 25's module is
/// 8 x 8 x 8 cells of 1.5 mm: 4.6 mm is 3.07 cells, 10.6 mm is 7.07, and the far wall at 12 mm is
/// the last cell, 7. A face's two axes are the format's — `x` faces read `[j, k]`, `y` faces
/// `[i, k]`, `z` faces `[i, j]` — so cell `[1, 2, 3]` is `[2, 3]`, `[1, 3]` and `[1, 2]`.
#[test]
fn a_cell_and_its_place_on_a_face_are_what_division_says() {
    let module = target(&scene("25-"), "module");
    assert_eq!(module.cell_at([0.0046, 0.0106, 0.012]), [3, 7, 7]);
    assert_eq!(module.cell_at([0.0, 0.0, 0.0]), [0, 0, 0]);
    for (face, want) in [
        (FaceSpec::XMin, [2, 3]),
        (FaceSpec::XMax, [2, 3]),
        (FaceSpec::YMin, [1, 3]),
        (FaceSpec::YMax, [1, 3]),
        (FaceSpec::ZMin, [1, 2]),
        (FaceSpec::ZMax, [1, 2]),
    ] {
        assert_eq!(module.on(face, [1, 2, 3]), want, "{face:?}");
    }
}

/// `text` with domain `domain`'s `cooling` list replaced by `list`.
fn with_cooling(text: &str, domain: &str, list: serde_json::Value) -> String {
    let mut v: serde_json::Value = serde_json::from_str(text).expect("parses");
    for d in v["domains"].as_array_mut().expect("domains") {
        if d["name"] == domain {
            d["cooling"] = list.clone();
        }
    }
    serde_json::to_string_pretty(&v).expect("serialises")
}

/// **Without a build, a part's solid cells are not known, and no area is guessed.** Scene 29's
/// bracket with its top cooled twice does not build, so its void map is not there; pointing its pad
/// at the top would have written the grid face's 27.04 cm² where the solver shares the area among
/// far fewer cells, and it is refused with the file left alone. Scene 21's wax, broken the same way,
/// has no part and no void region, so its cells are solid whether it builds or not and the pick
/// goes through.
#[test]
fn a_part_that_does_not_build_is_given_no_guessed_area() {
    let twice = serde_json::json!([
        { "face": "z-max", "ambient_c": 20.0, "convection_w_per_m2_k": 12.0, "area_cm2": 16.5 },
        { "face": "z-max", "ambient_c": 20.0, "convection_w_per_m2_k": 12.0, "area_cm2": 16.5 }
    ]);
    let broken = with_cooling(&scene("29-"), "bracket", twice.clone());
    let c = checked(&broken);
    assert!(c.error.is_some(), "the doubled face built");
    let bracket = c
        .targets
        .iter()
        .find(|t| t.name == "bracket")
        .expect("a target");
    assert!(bracket.void.is_none(), "a part's void map with no build");
    let refused = pick::set_face(&broken, bracket, "cooling", 0, FaceSpec::ZMin, None)
        .expect_err("an area was written for cells nobody knows are solid");
    assert!(refused.contains("does not build"), "{refused}");

    let wax = with_cooling(&scene("21-"), "wax", twice);
    let c = checked(&wax);
    assert!(c.error.is_some(), "the doubled face built");
    let t = c
        .targets
        .iter()
        .find(|t| t.name == "wax")
        .expect("a target");
    let (out, _) = pick::set_face(&wax, t, "cooling", 0, FaceSpec::XMax, None)
        .expect("a plain block's cells are known to be solid")
        .expect("a change");
    assert_eq!(
        entry(&out, "wax", "cooling", 0)["area_cm2"].as_f64(),
        Some(4.84)
    );
}

/// **Pointing again at a face whose stated area is not its cells' says so, and a face no cell is on
/// says it cools nothing.** Nothing is written either way — a stated area may be meant.
///
/// Scene 29's bracket is the second case, and was shipped that way. The part fills layers `k =
/// 0..10` of a 12-layer grid with 410 cells each and leaves the top two empty, so the top face,
/// layer 11, has no solid cell, and the 16.5 cm² still-air film stated there acts on nothing. The
/// format refuses a patch of void and not a whole face of it. Measured by counting the void map
/// layer by layer; pinned here as 0, which also holds a high face's layer on a part.
///
/// Its pad under the part is the first case only if its stated area differs, and it does not:
/// 0.96 cm² is its 24 cells. A face whose area is its cells', the format's own new entry on the
/// wax, has nothing to say.
#[test]
fn pointing_again_says_when_an_area_is_not_its_cells() {
    let text = scene("29-");
    let bracket = target(&text, "bracket");
    assert_eq!(
        pick::set_face(&text, &bracket, "cooling", 1, FaceSpec::ZMax, None),
        Ok(None)
    );
    assert_eq!(
        bracket.area_cm2(FaceSpec::ZMax, None).1,
        0,
        "solid cells on the bracket's top"
    );
    let note = pick::area_note(&text, &bracket, "cooling", 1).expect("the top says nothing");
    assert!(note.contains("cools nothing"), "{note}");
    assert_eq!(pick::area_note(&text, &bracket, "cooling", 0), None);

    // A stated area that is not its cells': the pad stated at 1.5 rather than its 0.96.
    let off = text.replacen("\"area_cm2\": 0.96", "\"area_cm2\": 1.5", 1);
    assert_ne!(off, text, "the pad's area was not found to change");
    let note = pick::area_note(&off, &bracket, "cooling", 0).expect("1.5 is not 24 cells");
    assert!(
        note.contains("1.5") && note.contains("24 solid cells") && note.contains("0.96"),
        "{note}"
    );
    let wax = editor_core::add_entry(&scene("21-"), "wax", "cooling").expect("added");
    assert_eq!(
        pick::area_note(&wax, &target(&wax, "wax"), "cooling", 0),
        None
    );
}

/// **Every kind a list can point at has a target.** [`templates::SITES`] names the lists; `check`
/// is what makes a domain of that kind something a click can land on. A site whose kind `check`
/// does not know would select fine and point at nothing.
///
/// [`templates::SITES`]: pantometry_world::templates::SITES
#[test]
fn every_kind_a_list_points_at_has_a_target() {
    let mut kinds: Vec<&str> = pantometry_world::templates::SITES
        .iter()
        .map(|s| s.kind)
        .collect();
    kinds.dedup();
    assert!(kinds.len() >= 2, "{kinds:?}");
    for kind in kinds {
        let text = pantometry_world::templates::scene(&[kind]);
        let c = checked(&text);
        assert!(
            c.targets.iter().any(|t| t.kind == kind),
            "a {kind} has a list a click can set and no target: {:?}",
            c.targets.iter().map(|t| &t.kind).collect::<Vec<_>>()
        );
    }
}

/// **Two crossing bars count their cells by hand, and a patch tells i from j.** Scene 30's busbars
/// are 8 x 8 x 5 cells of 4 mm, void but for one bar along `x` at `j = 3..5, k = 0..2` and one along
/// `y` at `i = 3..5, k = 3..5` — the only shipped void map that is not the same with `i` and `j`
/// swapped. So, from the scene's regions:
///
/// - `x-min` meets the first bar's end, 2 x 2 = 4 cells, and `y-min` the second's, 4;
/// - `z-min` is the first bar's underside, 8 x 2 = 16, and `z-max` the second's top, 16;
/// - cell `[0, 3, 0]` is copper and `[3, 0, 0]` is not;
/// - the patch `[0, 3]..[2, 5]` on `z-min` is 4 cells, 0.64 cm², where a map read with `i` and `j`
///   swapped finds none.
#[test]
fn two_crossing_bars_count_their_cells_by_hand() {
    let busbars = target(&scene("30-"), "busbars");
    for (face, want) in [
        (FaceSpec::XMin, 4),
        (FaceSpec::XMax, 4),
        (FaceSpec::YMin, 4),
        (FaceSpec::YMax, 4),
        (FaceSpec::ZMin, 16),
        (FaceSpec::ZMax, 16),
    ] {
        assert_eq!(busbars.area_cm2(face, None).1, want, "{face:?}");
    }
    assert!(busbars.is_solid([0, 3, 0]) && !busbars.is_solid([3, 0, 0]));
    let (area, solid) = busbars.area_cm2(FaceSpec::ZMin, Some(([0, 3], [2, 5])));
    assert_eq!(solid, 4);
    assert!((area - 0.64).abs() < 1e-15, "{area}");
}
