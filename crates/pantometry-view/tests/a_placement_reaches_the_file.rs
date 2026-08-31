//! **A placed domain has to arrive placed, in every writer.**
//!
//! `Panel::place` landed with only glTF wired to it. USD had no placement code at all — a `Mesh`
//! wrote `xformOpOrder = []` and a `BasisCurves` wrote nothing — so the defect the placement was
//! introduced to fix survived in `.usda`: two domains half a metre apart, both at their own local
//! origin, exported on top of each other. The commit that fixed it in glTF said it was fixed.
//!
//! # Two spellings of one transform
//!
//! glTF and USD compose a node the same way and spell its quaternion differently, so these tests
//! check the two writers against **each other** on the same frame. A copy of the components in the
//! wrong order is the failure with no symptom: the file loads, the prim is there, and the part is
//! turned somewhere nobody asked for.
//!
//! | | glTF | USD |
//! | --- | --- | --- |
//! | quaternion | `[x, y, z, w]` | `(w, x, y, z)` |
//! | composition | node `translation` × `rotation` | `xformOpOrder = [translate, orient]` |
//!
//! Neither row is guessed. `pxr` was asked to write a 60° turn about z and it produced
//! `(0.8660254037844387, 0, 0, 0.49999999999999994)` — every component distinct, so the text
//! cannot be read both ways — and the same file took `(1, 0, 0)` to `(1.0, 0.866025, 0)`, which is
//! the turn applied before the move. That is the one thing here measured outside this workspace;
//! everything below holds it in place without a USD library, which is what the rest of
//! `usd_is_usd.rs` explains.

use pantometry_scene::{Frame, Panel, PanelData, Placed};
use pantometry_view::{gltf, usda};

/// A 60° turn about z: `w = cos 30°`, `z = sin 30°`. Chosen so no two components are equal and no
/// ordering mistake can hide behind a symmetry — the first attempt used a quarter turn, whose `w`
/// and `z` are both `1/√2` and which reads identically in either order.
fn sixty_about_z() -> Placed {
    let h: f64 = std::f64::consts::PI / 6.0;
    Placed {
        at_m: [0.5, 0.0, 0.0],
        turn: [0.0, 0.0, h.sin(), h.cos()],
    }
}

/// One field panel, placed as given.
fn placed(place: Placed) -> Vec<Frame> {
    vec![Frame {
        time_s: 0.0,
        readings: Vec::new(),
        panels: vec![Panel {
            name: "part".into(),
            unit: "K",
            place,
            data: PanelData::Field {
                nx: 2,
                ny: 2,
                nz: 2,
                extent_m: [0.0, 0.0, 0.0, 0.04, 0.04, 0.04],
                values: (0..8).map(|i| 300.0 + i as f64).collect(),
            },
        }],
    }]
}

/// The numbers of the first `attr` in `doc`, in the order the file states them.
fn components(doc: &str, attr: &str) -> Vec<f64> {
    let at = doc
        .find(attr)
        .unwrap_or_else(|| panic!("{attr} is not in the file"));
    let open = doc[at..]
        .find('(')
        .or(doc[at..].find('['))
        .expect("a tuple")
        + at;
    let close = doc[open..].find([')', ']']).expect("a tuple that closes") + open;
    doc[open + 1..close]
        .split(',')
        .map(|s| s.trim().parse::<f64>().expect("a number"))
        .collect()
}

#[test]
fn the_identity_writes_what_it_always_wrote() {
    // Twenty-nine of the thirty shipped scenes are this case, so the bytes must not move: a file
    // that gained a transform saying "no transform" is a diff in every export ever taken. The
    // thirtieth is scene 30, which places two busbars and exists because the twenty-nine could
    // not fail — see `app/pantometry-world/scenes/30-two-phases-crossing-at-a-clearance.json`.
    let doc = usda("a run", &placed(Placed::HERE)).document;
    assert!(
        doc.contains("uniform token[] xformOpOrder = []"),
        "the empty order is what an unplaced prim has always written"
    );
    assert!(
        !doc.contains("xformOp:translate"),
        "an unplaced prim gained a transform"
    );
    assert!(!doc.contains("xformOp:orient"));
}

#[test]
fn usd_writes_the_real_part_first_and_gltf_writes_it_last() {
    // **The reordering, both ways round, out of one `Placed`.** Nothing here reads the writers'
    // own constants: the expected numbers are trigonometry, and each writer is asked for the
    // component the *other* format puts somewhere else.
    let h: f64 = std::f64::consts::PI / 6.0;
    let frames = placed(sixty_about_z());

    // **The tolerance is the file's own precision and nothing looser.** Both writers put numbers
    // through a `{:.6e}`, which is seven significant figures, so half a unit in the last place of
    // a component near 0.87 is 5e-8. Measured error here is 3.8e-9. A first attempt at 1e-9 was
    // tighter than the format can be and failed on a correct writer.
    const ULP: f64 = 5e-8;

    let u = components(&usda("a run", &frames).document, "quatd xformOp:orient");
    assert_eq!(u.len(), 4);
    assert!(
        (u[0] - h.cos()).abs() < ULP,
        "USD's first component is {}, not the real part {}",
        u[0],
        h.cos()
    );
    assert_eq!([u[1], u[2]], [0.0, 0.0]);
    assert!((u[3] - h.sin()).abs() < ULP, "USD's last should be z");

    let g = components(&gltf("a run", &frames[0]).document, "rotation");
    assert_eq!(g.len(), 4);
    assert!(
        (g[3] - h.cos()).abs() < ULP,
        "glTF's last component is {}, not the real part",
        g[3]
    );
    assert!((g[2] - h.sin()).abs() < ULP, "glTF's third should be z");

    // The same rotation, so the two files must be permutations of one another rather than of two
    // different quaternions. This is what a mistyped index breaks while both checks above pass.
    for (a, (x, y)) in u.iter().zip([g[3], g[0], g[1], g[2]]).enumerate() {
        assert!(
            (x - y).abs() < ULP,
            "component {a}: USD says {x}, glTF says {y} — a permutation would differ by ~0.5"
        );
    }
}

#[test]
fn usd_and_gltf_put_the_part_in_the_same_place() {
    let frames = placed(sixty_about_z());
    let u = components(
        &usda("a run", &frames).document,
        "double3 xformOp:translate",
    );
    let g = components(&gltf("a run", &frames[0]).document, "translation");
    assert_eq!(u, vec![0.5, 0.0, 0.0]);
    assert_eq!(u, g, "the two writers disagree about where the part is");
}

#[test]
fn the_order_turns_before_it_moves() {
    // `pxr` composes `[translate, orient]` as turn-then-move, which is a glTF node's order too.
    // Writing the list the other way round puts the part somewhere else entirely, and no amount of
    // checking the four numbers above would notice.
    let doc = usda("a run", &placed(sixty_about_z())).document;
    let order = doc
        .lines()
        .find(|l| l.contains("xformOpOrder"))
        .expect("an order");
    assert!(
        order.contains("xformOp:translate") && order.contains("xformOp:orient"),
        "the order is {order:?}"
    );
    assert!(
        order.find("translate") < order.find("orient"),
        "translate must be listed first, which is outermost: {order:?}"
    );
}

#[test]
fn two_placed_domains_do_not_export_on_top_of_each_other() {
    // **The measured defect, in the writer it survived in.** Before this, both prims wrote
    // `xformOpOrder = []` and a reader opening the file saw one part where the scene has two.
    let panel = |name: &str, place: Placed| Panel {
        name: name.into(),
        unit: "K",
        place,
        data: PanelData::Field {
            nx: 2,
            ny: 2,
            nz: 2,
            extent_m: [0.0, 0.0, 0.0, 0.03, 0.03, 0.03],
            values: (0..8).map(|i| 300.0 + i as f64).collect(),
        },
    };
    let frames = vec![Frame {
        time_s: 0.0,
        readings: Vec::new(),
        panels: vec![
            panel("near", Placed::HERE),
            panel(
                "far",
                Placed {
                    at_m: [0.5, 0.0, 0.0],
                    ..Placed::HERE
                },
            ),
        ],
    }];
    let doc = usda("a run", &frames).document;

    let near_at = doc.find("def Mesh \"near\"").expect("near");
    let far_at = doc.find("def Mesh \"far\"").expect("far");

    // The two are the same solid and their points are identical, so the transform is the only
    // thing that can tell them apart — which is exactly what the old file lacked.
    assert_eq!(
        components(&doc[far_at..], "double3 xformOp:translate"),
        vec![0.5, 0.0, 0.0]
    );
    assert!(
        !doc[near_at..far_at].contains("xformOp:translate"),
        "near was moved: it is at the origin"
    );
}

#[test]
fn a_path_is_placed_too_and_used_not_to_be() {
    // `BasisCurves` had no `xformOpOrder` at all, so a placed domain whose shape is paths was the
    // same defect one prim type further along — and a ray traced in a placed optic is that case.
    let frames = vec![Frame {
        time_s: 0.0,
        readings: Vec::new(),
        panels: vec![Panel {
            name: "rays".into(),
            unit: "nm",
            place: sixty_about_z(),
            data: PanelData::paths(vec![vec![[0.0, 0.0, 0.0], [1.0, 0.2, 0.0]]], vec![486.1]),
        }],
    }];
    let doc = usda("a run", &frames).document;
    let curves = doc.find("def BasisCurves").expect("a curve prim");
    assert_eq!(
        components(&doc[curves..], "double3 xformOp:translate"),
        vec![0.5, 0.0, 0.0],
        "a placed path stayed at the origin"
    );
}
