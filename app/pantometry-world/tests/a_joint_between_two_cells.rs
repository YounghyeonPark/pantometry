//! `contact` and `contact_w_per_m2_k`: the two places a scene can state a joint.
//!
//! A bolted joint, a thermal interface material, a solder layer, an oxide — all are tens of
//! microns in a part discretised at a millimetre, and this format could state none of them. A
//! `regions` entry is a box of *cells*, so the thinnest resistance it can express is `dx/k`: a
//! floor, not a discretisation error, since refining the mesh only lowers it and no grid a real
//! part can afford reaches 100 µm.
//!
//! # The two keys, and why there are two
//!
//! `contact` puts a resistance on a face **between two cells**. A mounting is on the one face that
//! has no cell on the other side of it — the outer one, where a part is bolted to its heatsink —
//! so it belongs to `cooling`, which is what already describes that face. On the power module that
//! ships here the mounting is 33% of the whole junction-to-ambient path and was absent.
//!
//! # What this file holds
//!
//! That the keys reach the domain, and that everything the format refuses elsewhere it refuses
//! here: a selection that picks nothing, a face named twice, a joint across a clearance, an index
//! that is really a boundary. The physics is `pantometry-thermal`'s
//! `a_joint_with_a_contact_resistance`, which checks the series against a closed form.

use pantometry::prelude::{Domain, Exchange};
use pantometry_world::{Scene, World};

/// A one-material block with whatever extra keys the case needs.
fn block(extra: &str) -> String {
    format!(
        r#"{{ "title": "t", "schedule": "multirate", "duration_s": 1e-4, "frames": 2,
          "domains": [{{ "kind": "block", "name": "b", "cells": [4, 4, 8],
            "cell_mm": 1.0, "initial_c": 20.0, "material": "copper"{extra} }}] }}"#
    )
}

/// The block a scene builds, or a panic saying why it did not.
fn built(scene: &str) -> pantometry::thermal::Solid3D {
    let parsed: Scene = serde_json::from_str(scene).expect("it parses");
    let world = World::build(parsed).expect("it builds");
    world
        .simulation()
        .domain_as::<pantometry::thermal::Solid3D>("b")
        .expect("it is a block")
        .clone()
}

/// **A stated joint reaches the face it names, and no other.**
///
/// The whole plane, then a patch of it. Read back through `face_conductance`, which the crate
/// documents as the number the sweep actually uses, so this is the operator and not the parser.
#[test]
fn a_contact_lands_on_the_faces_it_names() {
    let plain = built(&block(""));
    let bare = plain
        .face_conductance((0, 0, 3), (0, 0, 4))
        .expect("neighbours")
        .to_si();

    // A whole plane at z = 4.
    let whole = built(&block(
        r#", "contact": [{ "axis": "z", "at": 4, "w_per_m2_k": 5000.0 }]"#,
    ));
    let dx = 1e-3;
    let area = dx * dx;
    // `1/G = 1/G_bare + 1/(h·A)`, which is the series and nothing else.
    let expected = 1.0 / (1.0 / bare + 1.0 / (5000.0 * area));
    for j in 0..4 {
        for i in 0..4 {
            let g = whole
                .face_conductance((i, j, 3), (i, j, 4))
                .expect("neighbours")
                .to_si();
            assert!(
                (g / expected - 1.0).abs() < 1e-12,
                "({i},{j}): {g:.6e} against {expected:.6e}"
            );
        }
    }
    // And the plane above it is untouched.
    let above = whole
        .face_conductance((0, 0, 4), (0, 0, 5))
        .expect("neighbours")
        .to_si();
    assert!(above - bare == 0.0, "a neighbouring plane moved");

    // A patch: `from`/`to` are the two axes other than `axis`, in x-y-z order, so for z they are
    // `[i, j]`. The corner cell is in it and the far one is not.
    let patch = built(&block(
        r#", "contact": [{ "axis": "z", "at": 4, "w_per_m2_k": 5000.0,
              "from": [0, 0], "to": [2, 2] }]"#,
    ));
    let inside = patch
        .face_conductance((1, 1, 3), (1, 1, 4))
        .expect("neighbours")
        .to_si();
    let outside = patch
        .face_conductance((3, 3, 3), (3, 3, 4))
        .expect("neighbours")
        .to_si();
    println!("  patch: inside {inside:.6e}, outside {outside:.6e}, bare {bare:.6e}");
    assert!((inside / expected - 1.0).abs() < 1e-12, "the patch missed");
    assert!(outside - bare == 0.0, "the patch reached past its bounds");
}

/// **A mounting reaches the cooled face**, and a face with no mounting is unchanged.
#[test]
fn a_mounting_lands_on_the_face_that_is_cooled() {
    let cooling = |extra: &str| {
        block(&format!(
            r#", "cooling": [{{ "face": "z-max", "ambient_c": 20.0,
                 "convection_w_per_m2_k": 3000.0, "area_cm2": 0.16{extra} }}]"#
        ))
    };
    let plain = built(&cooling(""));
    assert_eq!(
        plain.mounting_of(pantometry::thermal::Face::ZMax),
        None,
        "a face with nothing stated carries no mounting"
    );
    let mounted = built(&cooling(r#", "contact_w_per_m2_k": 5000.0"#));
    assert_eq!(
        mounted.mounting_of(pantometry::thermal::Face::ZMax),
        Some(5000.0)
    );
    // The one face, and not the other five.
    assert_eq!(mounted.mounting_of(pantometry::thermal::Face::ZMin), None);
}

/// **Everything that cannot be meant is refused, and the message says which entry.**
///
/// The list is what this format already refuses everywhere else, asked of a new key: a selection
/// that picks nothing, a number that is not a conductance, an index that is really a boundary,
/// and the same face stated twice. A `regions` entry with mistyped bounds is the silent failure
/// this format is most careful about — a scene that runs, audits, renders and answers the wrong
/// question with nothing saying the coating was not applied — and a joint is the same shape.
#[test]
fn a_joint_that_cannot_be_meant_is_refused() {
    for (why, spec, says) in [
        (
            "a boundary rather than a face between two cells",
            r#", "contact": [{ "axis": "z", "at": 0, "w_per_m2_k": 5000.0 }]"#,
            "cooling",
        ),
        (
            "an index the block does not have",
            r#", "contact": [{ "axis": "z", "at": 8, "w_per_m2_k": 5000.0 }]"#,
            "last interior face",
        ),
        (
            "a negative conductance",
            r#", "contact": [{ "axis": "z", "at": 4, "w_per_m2_k": -5000.0 }]"#,
            "zero or more",
        ),
        (
            "a patch that runs backwards",
            r#", "contact": [{ "axis": "z", "at": 4, "w_per_m2_k": 5000.0,
                  "from": [2, 2], "to": [1, 1] }]"#,
            "selects no faces",
        ),
        (
            "a patch outside the plane",
            r#", "contact": [{ "axis": "z", "at": 4, "w_per_m2_k": 5000.0,
                  "from": [0, 0], "to": [9, 9] }]"#,
            "selects no faces",
        ),
        (
            "half a patch",
            r#", "contact": [{ "axis": "z", "at": 4, "w_per_m2_k": 5000.0, "from": [0, 0] }]"#,
            "both from and to",
        ),
        (
            "the same face twice",
            r#", "contact": [{ "axis": "z", "at": 4, "w_per_m2_k": 5000.0 },
                 { "axis": "z", "at": 4, "w_per_m2_k": 9000.0 }]"#,
            "already",
        ),
    ] {
        let scene: Scene =
            serde_json::from_str(&block(spec)).expect("it parses; the check is later");
        let Err(err) = World::build(scene) else {
            panic!("{why} must be refused");
        };
        println!("  {why}: {err}");
        assert!(
            err.contains("b/contact["),
            "{why}: the message should say which entry: {err}"
        );
        assert!(
            err.contains(says),
            "{why}: the message should say {says:?}: {err}"
        );
    }

    // And the mounting's own number.
    let bad = block(
        r#", "cooling": [{ "face": "z-max", "ambient_c": 20.0, "convection_w_per_m2_k": 3000.0,
             "area_cm2": 0.16, "contact_w_per_m2_k": -1.0 }]"#,
    );
    let scene: Scene = serde_json::from_str(&bad).expect("parses");
    let Err(err) = World::build(scene) else {
        panic!("a negative mounting must be refused");
    };
    println!("  a negative mounting: {err}");
    assert!(err.contains("b/cooling[0]"), "which entry: {err}");
}

/// **A joint across a clearance is refused, because it does nothing.**
///
/// The harmonic face mean is already zero where a cell holds nothing, so a contact stated there
/// changes no number at all. That is precisely the shape this format refuses elsewhere: a file
/// that states a joint, runs, audits and answers as though it had none. Refused when *every* face
/// the entry names is dead; noted when only some are, because a patch across a part with voids in
/// it may legitimately cross both.
#[test]
fn a_joint_across_a_clearance_is_refused_and_a_partial_one_is_noted() {
    let dead = block(
        r#", "regions": [{ "material": "void", "from": [0, 0, 4], "to": [4, 4, 5] }],
             "contact": [{ "axis": "z", "at": 4, "w_per_m2_k": 5000.0 }]"#,
    );
    let scene: Scene = serde_json::from_str(&dead).expect("parses");
    let Err(err) = World::build(scene) else {
        panic!("a joint that does nothing must be refused");
    };
    println!("  all dead: {err}");
    assert!(err.contains("b/contact[0]"), "which entry: {err}");
    assert!(
        err.contains("nothing on one side"),
        "the message should say why: {err}"
    );

    // Half of the plane is a clearance and half of it conducts: a note, and it still builds.
    let partial = block(
        r#", "regions": [{ "material": "void", "from": [0, 0, 4], "to": [2, 4, 5] }],
             "contact": [{ "axis": "z", "at": 4, "w_per_m2_k": 5000.0 }]"#,
    );
    let scene: Scene = serde_json::from_str(&partial).expect("parses");
    let world = World::build(scene).expect("half a plane still conducts");
    let notes = world.notes();
    let note = notes
        .iter()
        .find(|n| n.contains("contact[0]"))
        .unwrap_or_else(|| panic!("no note about the half-dead plane: {notes:?}"));
    println!("  half dead: {note}");
    assert!(
        note.contains("16 faces") && note.contains("8 of them"),
        "{note}"
    );
}

/// **A cooled face may name the box on it that is bolted down.**
///
/// The other half of `contact`: a joint between two cells is `contact`, and a joint between a face
/// and whatever cools it is `cooling`'s `contact_w_per_m2_k` — but neither could say that only
/// *part* of a face touches anything. A bracket is bolted at pads, and cooled over its whole
/// footprint it has no route along its own shape: every cell sheds where it stands. Measured on
/// `29-a-designed-bracket-becomes-cells`, which rasterised 4 100 cells from an STL to hold a range
/// of 0.067 K.
///
/// Stating a smaller `area_cm2` is not the same thing, and that is what the second half checks:
/// the area is divided among the cells on the face, so a tenth of the area is a tenth of the
/// conductance **spread over all of it**.
#[test]
fn a_cooled_face_can_name_the_box_on_it_that_is_bolted() {
    // **Hot, or this measures nothing.** The first version left the block at the ambient it was
    // losing to, so both sides lost 0.000000000 J and the equality below was satisfied by two
    // zeroes. A check whose subject is a difference needs the difference to exist.
    let cooled = |extra: &str| {
        block(&format!(
            r#", "cooling": [{{ "face": "z-min", "ambient_c": 20.0,
                 "convection_w_per_m2_k": 3000.0, "area_cm2": 0.04{extra} }}]"#
        ))
        .replace("\"initial_c\": 20.0", "\"initial_c\": 200.0")
    };
    let whole = built(&cooled(""));
    assert_eq!(whole.patch_on(pantometry::thermal::Face::ZMin), None);
    let pad = built(&cooled(r#", "from": [1, 1], "to": [3, 3]"#));
    assert_eq!(
        pad.patch_on(pantometry::thermal::Face::ZMin),
        Some([1, 1, 3, 3]),
        "the patch did not reach the domain"
    );

    // The stated area is the same, so at a uniform temperature the same joules leave; what the
    // patch decides is which cells carry them.
    let (mut whole, mut pad) = (whole, pad);
    let dt = pantometry::units::Time::from_si(
        pad.max_stable_dt(pantometry::units::Time::ZERO).to_si() * 0.5,
    );
    let mut bus = Exchange::new();
    whole
        .step(pantometry::units::Time::ZERO, dt, &mut bus)
        .expect("stable");
    pad.step(pantometry::units::Time::ZERO, dt, &mut bus)
        .expect("stable");
    let (a, b) = (whole.lost_energy().to_si(), pad.lost_energy().to_si());
    println!("  one step: the whole face lost {a:.9} J, the pad lost {b:.9} J");
    assert!(
        a > 0.0,
        "nothing was lost, so nothing was compared: {a:.12} J"
    );
    assert!(
        (a - b).abs() < 1e-12 * a.abs().max(1.0),
        "a patch changed the total: {a:.12} against {b:.12}"
    );
}

/// **A patch that cannot be meant is refused, and one that names no solid cell most of all.**
///
/// A part rasterised from an STL is mostly *not there*, so a pad placed by eye can miss it
/// entirely — and a scene whose cooling names only void runs insulated, warms for as long as it
/// runs and answers a different question in silence. That is the same failure `area_cm2 = 0` is
/// already refused for, arriving by a different door.
#[test]
fn a_cooled_patch_that_cannot_be_meant_is_refused() {
    for (why, spec, says) in [
        (
            "a box that runs backwards",
            r#", "from": [3, 3], "to": [1, 1]"#,
            "selects no cells",
        ),
        (
            "a box outside the face",
            r#", "from": [0, 0], "to": [9, 9]"#,
            "selects no cells",
        ),
        ("half a box", r#", "from": [0, 0]"#, "both from and to"),
    ] {
        let json = block(&format!(
            r#", "cooling": [{{ "face": "z-min", "ambient_c": 20.0,
                 "convection_w_per_m2_k": 3000.0, "area_cm2": 0.04{spec} }}]"#
        ));
        let scene: Scene = serde_json::from_str(&json).expect("it parses; the check is later");
        let Err(err) = World::build(scene) else {
            panic!("{why} must be refused");
        };
        println!("  {why}: {err}");
        assert!(err.contains("b/cooling[0]"), "{why}: which entry: {err}");
        assert!(err.contains(says), "{why}: should say {says:?}: {err}");
    }

    // A pad on a corner of the face that has been voided away: it selects cells, and none of them
    // are there.
    let json = block(
        r#", "regions": [{ "material": "void", "from": [0, 0, 0], "to": [2, 2, 1] }],
             "cooling": [{ "face": "z-min", "ambient_c": 20.0,
               "convection_w_per_m2_k": 3000.0, "area_cm2": 0.04,
               "from": [0, 0], "to": [2, 2] }]"#,
    );
    let scene: Scene = serde_json::from_str(&json).expect("parses");
    let Err(err) = World::build(scene) else {
        panic!("a pad on nothing must be refused");
    };
    println!("  a pad on void: {err}");
    assert!(
        err.contains("names no solid cell"),
        "the message should say why: {err}"
    );
}
