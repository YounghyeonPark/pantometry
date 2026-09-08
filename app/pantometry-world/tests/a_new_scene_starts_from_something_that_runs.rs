//! The editor's "New project" chooser, against the scenes that prove its numbers.
//!
//! [`templates::scene`] turns a set of kinds into a scene to open. Two things in it are numbers
//! rather than text — `duration_s` and `frames` — and a starting point whose duration nothing has
//! ever run is a starting point that opens on a refusal. So they are not chosen here: each is
//! taken from a shipped scene that uses that kind, and this holds them against those scenes.
//!
//! That is the shape `every_domain_has_a_template` established, for the same reason: a table
//! maintained by hand is a table that is silently wrong between the change and the noticing, and
//! the thirty scenes are the set that is maintained by something else.
//!
//! # Not under `wasm32`
//!
//! Every test here reads the scenes off a disk, and a `wasm32` target has none.

#![cfg(not(target_family = "wasm"))]

use pantometry_world::templates::{self, TEMPLATES};
use pantometry_world::{OnDisk, Scene, World};

/// Every `(duration_s, frames)` pair each kind appears with, across the shipped scenes.
fn schedules_per_kind() -> std::collections::BTreeMap<String, Vec<(f64, u64)>> {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("scenes");
    let mut out: std::collections::BTreeMap<String, Vec<(f64, u64)>> = Default::default();
    let mut scenes = 0;
    for entry in std::fs::read_dir(&dir).expect("the scenes directory reads") {
        let path = entry.expect("a directory entry").path();
        if path.extension().is_none_or(|x| x != "json") {
            continue;
        }
        scenes += 1;
        let text = std::fs::read_to_string(&path).expect("a scene reads");
        let value: serde_json::Value =
            serde_json::from_str(&text).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let (Some(secs), Some(frames)) = (value["duration_s"].as_f64(), value["frames"].as_u64())
        else {
            continue;
        };
        for domain in value["domains"].as_array().into_iter().flatten() {
            if let Some(kind) = domain["kind"].as_str() {
                out.entry(kind.to_string())
                    .or_default()
                    .push((secs, frames));
            }
        }
    }
    // A directory that read as empty would make every comparison below vacuously true, which is
    // the shape of a check that turned itself off.
    assert!(
        scenes >= 28,
        "only {scenes} scenes found in {} — this test compares against them",
        dir.display()
    );
    out
}

/// **Every number a new scene starts from is a number some shipped scene runs with.**
#[test]
fn the_schedule_a_template_suggests_is_one_a_scene_uses() {
    let shipped = schedules_per_kind();
    for t in TEMPLATES {
        let used = shipped
            .get(t.kind)
            .unwrap_or_else(|| panic!("{}: no shipped scene uses this kind", t.kind));
        assert!(
            used.contains(&(t.duration_s, t.frames as u64)),
            "{}: the template suggests {} s in {} frames, which no scene that uses it runs. \
             The scenes run it at {used:?}",
            t.kind,
            t.duration_s,
            t.frames
        );
    }
}

/// **A scene made of one kind parses, and builds unless the format will not let it stand alone.**
///
/// `structure` names the block it `follows` and is refused by the build; `beam` names what it
/// shines `onto`, builds, and is stopped by the conservation audit at the first step. Both are
/// left as the templates write them — see [`templates::scene`] — so this pins that the set of
/// kinds a new scene cannot open cleanly is exactly those two.
#[test]
fn a_new_scene_of_one_kind_builds_or_names_what_it_needs() {
    let mut refused = Vec::new();
    for t in TEMPLATES {
        let text = templates::scene(&[t.kind]);
        let scene: Scene = serde_json::from_str(&text)
            .unwrap_or_else(|e| panic!("{}: a new scene does not parse: {e}\n{text}", t.kind));
        assert_eq!(
            scene.domains.len(),
            1,
            "{}: one kind should make one domain",
            t.kind
        );
        match World::build_with(scene, &OnDisk) {
            Ok(_) => {}
            Err(why) if why.contains("follows") => refused.push(t.kind),
            Err(why) => panic!(
                "{}: a new scene of one does not build: {why}\n{text}",
                t.kind
            ),
        }
    }
    assert_eq!(
        refused,
        vec!["structure"],
        "the set of kinds a new scene cannot open has changed"
    );
}

/// **A scene made of several kinds parses and builds**, at the shortest of their durations.
///
/// Three that share a timescale and are already used together in the shipped scenes.
#[test]
fn a_new_scene_of_several_kinds_builds_at_the_shortest_duration() {
    let text = templates::scene(&["bar", "heater", "lump"]);
    let scene: Scene = serde_json::from_str(&text)
        .unwrap_or_else(|e| panic!("three kinds do not parse: {e}\n{text}"));
    assert_eq!(scene.domains.len(), 3, "three kinds, three domains\n{text}");
    // `lump` is the shortest of the three at 1.0 s; `bar` and `heater` are 4.0.
    assert!(
        text.contains("\"duration_s\": 1"),
        "the shortest duration is not the one written:\n{text}"
    );
    World::build_with(scene, &OnDisk).unwrap_or_else(|e| panic!("{e}\n{text}"));
}

/// **Exactly two kinds name a partner**, and they are the two the build and the audit refuse.
///
/// `needs_a_partner` reads the template's own text for `onto` and `follows`. Deriving it means a
/// third such key would be missed, so the set it produces is pinned against the set the test above
/// measures from the format itself.
#[test]
fn exactly_two_kinds_name_a_partner() {
    let named: Vec<&str> = TEMPLATES
        .iter()
        .filter(|t| t.needs_a_partner())
        .map(|t| t.kind)
        .collect();
    assert_eq!(named, vec!["beam", "structure"]);
}

/// **The span is a ratio, and it is 1 for one kind and enormous for the wrong two.**
///
/// The number the chooser says out loud. A quantum `well` settles in 2e-13 s and an `orbit` takes
/// 7200: there is no `duration_s` at which both are a simulation, and the format cannot refuse the
/// scene because it is well formed.
///
/// **This asked about `atoms` against a thermal `network` and wanted `> 1e12`.** That pair was
/// 3e14 apart because `atoms` held `duration_s: 6e-12` against a domain whose own time unit is one
/// second — a units error, not a timescale. Corrected to 6.0 they are **300** apart, under the
/// chooser's own 1e3 threshold, and this test was measuring the bug. The table's real extremes are
/// `well` and `orbit`, and they do not depend on anybody's units being right.
#[test]
fn the_timescale_span_says_when_a_set_cannot_run_together() {
    assert_eq!(templates::timescale_span(&["room"]), 1.0);
    assert!(
        (templates::timescale_span(&["bar", "heater"]) - 1.0).abs() < 1e-12,
        "two kinds a scene already runs together are not one apart"
    );
    let far = templates::timescale_span(&["well", "orbit"]);
    assert!(far > 1e12, "2e-13 s against two hours came out as {far}");
    // And the extremes are the extremes: nothing in the table is further apart than these two, or
    // the sentence above names the wrong pair and the count in `templates`' own doc is stale
    // again. Measured from the table rather than remembered.
    let every: Vec<&str> = templates::TEMPLATES.iter().map(|t| t.kind).collect();
    let widest = templates::timescale_span(&every);
    assert!(
        (widest - far).abs() / far < 1e-12,
        "the widest pair in the table spans {widest}, and `well` against `orbit` spans {far}"
    );
    // Sixteen orders, which is what three doc comments now say.
    assert!(
        (3.6e16..3.7e16).contains(&widest),
        "the table spans {widest}, and the docs say 3.6e16"
    );
}

/// **How many pictures you ask for cannot change the answer.**
///
/// It could, and by a lot. The run advanced by `duration_s / frames` and each domain subdivided
/// that into whole substeps no longer than its own stability limit, so the step was
/// `window / ceil(window / limit)` — a function of the frame count. Measured on
/// `15-a-hot-spot-in-a-block`, whose peak reads 5.177, 5.714, 7.027, 7.696, 7.940, 8.065, 8.129,
/// 8.161, 8.177 and 8.186 K at 2, 5, 11, 26, 51, 101, 201, 401, 801 and 1601 frames: **58%**
/// between the ends, textbook first order, and the shipped `frames: 11` sat 14% below its own
/// grid's converged answer.
///
/// `window_s` is what separates the two. With one stated, the run takes whole steps of it and
/// `frames` only chooses which are photographed.
///
/// The scene here is a block with a hot spot rather than a room mode on purpose: a diffusion
/// problem's error is first order in the step, so a frame count that leaked into the step would
/// show up immediately. A wave problem's would partly cancel by phase and could hide it.
#[test]
fn the_frame_count_is_pictures_and_not_physics() {
    let text = r#"{
  "title": "one hot cell, photographed at several rates",
  "duration_s": 0.006,
  "frames": 11,
  "window_s": 3.75e-6,
  "domains": [
    { "kind": "block", "name": "block", "cells": [9, 9, 9], "cell_mm": 1.0,
      "initial_c": 20.0, "material": "aluminium",
      "hot_spot": { "at": [4, 4, 4], "above_k": 60.0 } }
  ]
}"#;
    let mut answers = Vec::new();
    for frames in [2usize, 11, 88] {
        let mut scene: Scene = serde_json::from_str(text).expect("the scene parses");
        scene.frames = frames;
        let mut world = World::build(scene).expect("it builds");
        let out = world.run().expect("it runs");
        assert_eq!(
            out.len(),
            frames + 1,
            "{frames} frames asked for, {} captured",
            out.len()
        );
        let last = out.last().expect("a run produces frames");
        let peak = last.panels[0]
            .values()
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        answers.push(peak);
    }
    // Bit-for-bit, not to a tolerance: the steps taken are identical, so the arithmetic follows
    // the same path and a tolerance here would hide exactly the coupling this test is about.
    assert_eq!(
        answers[0], answers[1],
        "2 frames and 11 gave {} and {}",
        answers[0], answers[1]
    );
    assert_eq!(
        answers[1], answers[2],
        "11 frames and 88 gave {} and {}",
        answers[1], answers[2]
    );

    // And without the key the old coupling is still there, which is what makes the key
    // load-bearing rather than decorative — and is why every scene written before it is unchanged.
    let mut coupled = Vec::new();
    for frames in [11usize, 88] {
        let mut scene: Scene = serde_json::from_str(text).expect("the scene parses");
        scene.window_s = None;
        scene.frames = frames;
        let mut world = World::build(scene).expect("it builds");
        let out = world.run().expect("it runs");
        let last = out.last().expect("a run produces frames");
        coupled.push(
            last.panels[0]
                .values()
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max),
        );
    }
    assert!(
        (coupled[1] - coupled[0]).abs() > 1.0,
        "without a window the frame count used to move the peak by kelvin; it moved {:.4}",
        (coupled[1] - coupled[0]).abs()
    );
}
