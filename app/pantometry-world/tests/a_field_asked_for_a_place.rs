//! A field could be asked for its mean, its peak and its coldest, and not for a place.
//!
//! A design number is almost always at a place — the junction, the sensor, the point under the
//! baseplate — and `24-a-power-module-junction-to-ambient` read its junction as the block's `peak`,
//! the hottest cell *anywhere*. That is right for a single die by coincidence of that geometry and
//! has no second answer for a module with two.
//!
//! [`Scene::probes`](pantometry_world::Scene::probes) names points. The decision the key turns on
//! is that they are in **millimetres**: a point stated in cells moves when the grid is refined, so
//! a resolution sweep would compare two different places and call the difference discretisation —
//! the same error a notch's blocked cells made until they learned to refine into their eight
//! children.

use pantometry_world::{Scene, World};

fn scene(json: &str) -> Scene {
    serde_json::from_str(json).expect("the test scene parses")
}

/// A 12 mm cube of aluminium at 8 cells, warmed at one end, with whatever probes are asked for.
fn block(probes: &str) -> String {
    format!(
        r#"{{
  "title": "a block with places in it",
  "duration_s": 0.2,
  "frames": 4,
  "domains": [
    {{ "kind": "block", "name": "part", "cells": [8, 8, 8], "cell_mm": 1.5,
       "material": "aluminium", "initial_c": 20.0,
       "regions": [ {{ "material": "aluminium", "from": [0, 0, 6], "to": [8, 8, 8],
                       "initial_c": 120.0 }} ] }}
  ]{probes}
}}"#
    )
}

fn read(frame: &pantometry_world::Frame, domain: &str, label: &str) -> f64 {
    frame
        .readings
        .iter()
        .find(|r| r.domain == domain && r.label == label)
        .unwrap_or_else(|| panic!("no {domain}/{label}"))
        .value
}

/// **A probe at the hottest place reads what `peak` reads**, which is two paths to one number.
///
/// `peak` comes from the cells the domain holds; a probe comes from
/// [`ScalarField::at`](pantometry::core::ScalarField::at), which interpolates. They are different
/// code reaching the same value, and this repository has already found one case where two such
/// paths disagreed — a panel sampled at cell *boundaries* while the readings beside it reported the
/// cell maximum, one frame carrying two answers.
///
/// The same at the cold end against `coldest`.
#[test]
fn a_probe_at_the_hottest_place_agrees_with_the_peak() {
    let text = block(
        r#","probes": {
    "top":    { "in": "part", "at_mm": [6.0, 6.0, 12.0] },
    "bottom": { "in": "part", "at_mm": [6.0, 6.0, 0.0] }
  }"#,
    );
    let mut world = World::build(scene(&text)).expect("it builds");
    let frames = world.run().expect("it conserves");
    let last = frames.last().expect("frames");

    let (top, peak) = (read(last, "part", "top"), read(last, "part", "peak"));
    let (bottom, coldest) = (read(last, "part", "bottom"), read(last, "part", "coldest"));
    println!(
        "  top {top:.6} against peak {peak:.6}; bottom {bottom:.6} against coldest {coldest:.6}"
    );
    assert!(
        (top - peak).abs() < 1e-9,
        "the top face reads {top} and the peak is {peak}"
    );
    assert!(
        (bottom - coldest).abs() < 1e-9,
        "the bottom face reads {bottom} and the coldest is {coldest}"
    );
}

/// **A probe reports celsius, because a reading reports celsius.**
///
/// `Reading::value` states the one exception this workspace makes to SI — "temperatures are
/// celsius, because that is the unit a column of them is read in" — and `peak`, `mean` and
/// `coldest` all follow it. A field's `at` answers in the SI base unit, so a probe reporting what
/// it returns put **500.3026 K beside a peak of 227.1526 C** in one table: the same number twice,
/// in two units, in the column a reader compares down.
#[test]
fn a_probe_reports_the_unit_the_column_is_read_in() {
    let text = block(r#","probes": { "top": { "in": "part", "at_mm": [6.0, 6.0, 12.0] } }"#);
    let mut world = World::build(scene(&text)).expect("it builds");
    let frames = world.run().expect("it conserves");
    let last = frames.last().expect("frames");
    let probe = last
        .readings
        .iter()
        .find(|r| r.domain == "part" && r.label == "top")
        .expect("the probe is reported");
    let peak = last
        .readings
        .iter()
        .find(|r| r.domain == "part" && r.label == "peak")
        .expect("the peak is reported");
    assert_eq!(
        probe.unit, peak.unit,
        "a probe and a peak are the same quantity and must be the same unit"
    );
    assert_eq!(probe.unit, "C");
    // And it is not a kelvin wearing a celsius label.
    assert!(
        probe.value < 200.0,
        "the probe reads {} {}, which is a kelvin",
        probe.value,
        probe.unit
    );
}

/// **A probe is at the same millimetres on a finer grid**, which is why the key takes millimetres.
///
/// `Scene::refined` doubles every grid. A probe stated in cells would follow the index and land
/// somewhere else; stated in millimetres it is the same point, so the two readings are a
/// convergence measurement of one design number rather than a comparison of two places.
#[test]
fn a_probe_is_the_same_place_on_a_finer_grid() {
    let text = block(r#","probes": { "middle": { "in": "part", "at_mm": [6.0, 6.0, 6.0] } }"#);
    let coarse = scene(&text);
    let fine = coarse.refined().expect("a block refines");

    // The scene's own key is untouched by refinement — millimetres do not scale.
    assert_eq!(
        fine.probes, coarse.probes,
        "refinement moved a probe that is stated in millimetres"
    );

    let value = |s: &Scene| {
        let mut world = World::build(s.clone()).expect("it builds");
        let frames = world.run().expect("it conserves");
        read(frames.last().expect("frames"), "part", "middle")
    };
    let (a, b) = (value(&coarse), value(&fine));
    println!("  the middle reads {a:.6} C at 1.5 mm and {b:.6} C at 0.75 mm");
    // The same place, so the two are a convergence pair: close, and not identical, because the
    // field between them is genuinely different.
    assert!(
        (a - b).abs() < 5.0,
        "the same point read {a} and {b}, which is not one point"
    );
    assert!(
        (a - b).abs() > 1e-9,
        "refining changed nothing, so the probe is not reading the grid at all"
    );
}

/// **A probe is at the millimetres it says**, which the faces cannot show.
///
/// `ScalarField::at` **clamps at the faces**, so a point asked for at `z = 12 m` instead of
/// `12 mm` comes back as the top face — the same number. Every probe in the tests above sits on a
/// face or at the centre, so a sabotage that dropped the `1e-3` and read millimetres as metres
/// passed all four. What tells them apart is a point between two cells with a gradient across it.
///
/// The block is 12 mm of aluminium warmed at the top. A quarter of the way up is nowhere near
/// either face and nowhere near the middle, so it has one value and the clamped answers have
/// others.
#[test]
fn a_probe_is_at_the_millimetres_it_says() {
    let text = block(
        r#","probes": {
    "quarter": { "in": "part", "at_mm": [6.0, 6.0, 3.0] },
    "top":     { "in": "part", "at_mm": [6.0, 6.0, 12.0] },
    "middle":  { "in": "part", "at_mm": [6.0, 6.0, 6.0] }
  }"#,
    );
    let mut world = World::build(scene(&text)).expect("it builds");
    let frames = world.run().expect("it conserves");
    let last = frames.last().expect("frames");
    let (quarter, middle, top) = (
        read(last, "part", "quarter"),
        read(last, "part", "middle"),
        read(last, "part", "top"),
    );
    println!("  quarter {quarter:.6}, middle {middle:.6}, top {top:.6} C");

    // Three distinct places in a field with a gradient: three distinct numbers, ordered.
    assert!(
        quarter < middle && middle < top,
        "a warmed-at-the-top block should rise: {quarter} then {middle} then {top}"
    );
    // And the quarter point is neither of the two the clamp would return — which is what a
    // millimetre read as a metre becomes.
    let coldest = read(last, "part", "coldest");
    assert!(
        (quarter - top).abs() > 1.0 && (quarter - coldest).abs() > 0.1,
        "the quarter point reads {quarter}, the top {top} and the coldest {coldest}"
    );
}

/// Every way a probe can be wrong says which, and none of them is a number from nowhere.
#[test]
fn a_probe_that_cannot_be_read_says_why() {
    // Outside the part. This is the one that matters: `at` answers anywhere it is asked, so a
    // point past the die would appear in every frame, in the CSV and in the sweep, and converge.
    let why = World::build(scene(&block(
        r#","probes": { "junction": { "in": "part", "at_mm": [6.0, 6.0, 14.0] } }"#,
    )))
    .err()
    .expect("a probe outside the part must be refused");
    assert!(why.contains("spans"), "{why}");
    assert!(why.contains("[12.0, 12.0, 12.0]"), "{why}");

    // A face is inside, to a nanometre, so a probe written at the surface of a part is not a
    // rounding away from being refused.
    World::build(scene(&block(
        r#","probes": { "junction": { "in": "part", "at_mm": [0.0, 0.0, 12.0] } }"#,
    )))
    .expect("a probe on a face is inside");

    // A name that is not there.
    let why = World::build(scene(&block(
        r#","probes": { "j": { "in": "prt", "at_mm": [1.0, 1.0, 1.0] } }"#,
    )))
    .err()
    .expect("a probe on nothing must be refused");
    assert!(why.contains("does not define"), "{why}");

    // A domain with readings and no places.
    let heated = r#"{
  "title": "a heater and a bar", "duration_s": 4.0, "frames": 4,
  "domains": [
    { "kind": "heater", "name": "element", "watts": 2.0, "reserve_j": 6.0 },
    { "kind": "bar", "name": "bar", "length_mm": 20.0, "cells": 10, "area_mm2": 100.0,
      "initial_c": 20.0 }
  ],
  "probes": { "j": { "in": "element", "at_mm": [0.0, 0.0, 0.0] } }
}"#;
    let why = World::build(scene(heated))
        .err()
        .expect("a heater has no place to read");
    assert!(why.contains("no field to read at a place"), "{why}");
    assert!(why.contains("these do: bar"), "{why}");

    // And a misspelled key is serde's to catch, with what it expected.
    let e = serde_json::from_str::<Scene>(&block(
        r#","probes": { "j": { "in": "part", "at": [1.0, 1.0, 1.0] } }"#,
    ))
    .expect_err("a misspelled probe must not parse");
    assert!(e.to_string().contains("unknown field `at`"), "{e}");
}
