//! **A schedule no domain could survive is refused when the scene is built, not when it is run.**
//!
//! `FRICTION.md`'s finding 8, which stood unactioned since it was written. A scene picks its
//! schedule by name, and `staggered` with a half-second frame is hundreds of times a bar's
//! explicit-diffusion limit. The run was refused — correctly, by name, with the limit and the
//! value, which is the whole argument for this library and it works. But it was refused *when the
//! step was taken*, and the person who wants to know is the one editing the file.
//!
//! `Domain::max_stable_dt` is public, so an application loading scenes from disk can ask. This is
//! the application; the editor checks as you type through the same `World::build_with`.
//!
//! # Necessary and not sufficient, which the tests below say twice
//!
//! `max_stable_dt` takes a `now`. It is the largest step **from here**, and a domain whose
//! conductance rises as it warms has a limit that tightens under it. The build sees only the
//! initial state, so this catches the scene that could never have worked and does not replace the
//! refusal inside `step` — which stays, and is the one that is complete.

use pantometry_world::{Scene, World};

fn scene(json: &str) -> Scene {
    serde_json::from_str(json).expect("the test scene parses")
}

/// A heater and a bar — the pair finding 8 was written about — under a chosen schedule and frame
/// count. 61 cells over 20 mm is a tight limit and deliberately so: the finding's own bar.
fn heated(schedule: &str, seconds: f64, frames: usize) -> String {
    format!(
        r#"{{
  "title": "a heater pays joules onto the bus and a bar takes them",
  "schedule": "{schedule}",
  "duration_s": {seconds},
  "frames": {frames},
  "conservation_tolerance": 1e-9,
  "domains": [
    {{ "kind": "heater", "name": "element", "watts": 2.0, "reserve_j": 6.0 }},
    {{ "kind": "bar", "name": "bar", "length_mm": 20.0, "cells": 61,
      "area_mm2": 100.0, "initial_c": 20.0 }}
  ]
}}"#
    )
}

/// The `N` out of "raise `frames` to at least N".
fn suggested(err: &str) -> usize {
    let at = err
        .find("at least ")
        .unwrap_or_else(|| panic!("the refusal suggests nothing: {err}"));
    let rest = &err[at + "at least ".len()..];
    rest[..rest.find(' ').expect("a number then a space")]
        .parse()
        .expect("a frame count")
}

#[test]
fn a_frame_the_bar_cannot_survive_is_refused_at_build() {
    // Four seconds in eight frames is the half-second frame the finding names.
    let Err(err) = World::build(scene(&heated("staggered", 4.0, 8))) else {
        panic!("a staggered half-second frame was accepted for a 61-cell bar");
    };

    // The domain, the two numbers, and the ratio — everything a person needs to decide what to
    // change. Measured 642x, which is not the finding's "thirty-eight times": that was a different
    // bar, and the number here is this bar's.
    assert!(
        err.contains("bar:"),
        "the refusal does not name the domain: {err}"
    );
    assert!(
        err.contains("5.000e-1 s window"),
        "not the frame window: {err}"
    );
    // The wording `verify` used for the same hazard, kept when the two checks became one.
    assert!(
        err.contains("does not subcycle") && err.contains("silently unstable"),
        "the refusal dropped the battery's words: {err}"
    );
    assert!(err.contains("642"), "not the ratio: {err}");
    // And the caveat, in the message rather than only in a doc comment: somebody reading this in a
    // terminal has to know it is not the whole check.
    assert!(
        err.contains("initial") && err.contains("when the step is taken"),
        "the refusal does not say it is the necessary half: {err}"
    );
}

#[test]
fn the_frame_count_it_suggests_actually_builds() {
    // **The assertion that makes the suggestion a suggestion.** A message that names a number
    // nobody tried is a plausible number, and this repository has shipped one of those before.
    let Err(err) = World::build(scene(&heated("staggered", 4.0, 8))) else {
        panic!("accepted");
    };
    let n = suggested(&err);
    assert!(n > 8, "the suggestion is not an increase: {n}");
    World::build(scene(&heated("staggered", 4.0, n))).unwrap_or_else(|e| {
        panic!("the refusal suggested {n} frames and {n} frames is refused too: {e}")
    });

    // And one fewer than it suggests is still refused, so the number is the threshold rather than
    // a comfortable round-up. `ceil` makes it exact: `duration / limit` frames is the first that
    // fits.
    assert!(
        World::build(scene(&heated("staggered", 4.0, n - 1))).is_err(),
        "{} frames was accepted, so {n} is not the boundary",
        n - 1
    );
}

#[test]
fn one_way_does_not_subcycle_either_and_is_refused_too() {
    // The check `verify` carried exempted **only** `multirate`, and the first draft of this one
    // exempted everything but `staggered` — which would have let a `one-way` scene through, since
    // it takes the whole window in one step just the same. Merging the two is what surfaced it.
    let json = heated("staggered", 4.0, 8).replace("\"staggered\"", "\"one-way\"");
    let Err(err) = World::build(scene(&json)) else {
        panic!("one-way takes the whole window in one step and was accepted");
    };
    assert!(
        err.contains("OneWay"),
        "the refusal names the schedule: {err}"
    );
}

#[test]
fn multirate_is_never_refused_here_because_it_subcycles() {
    // The same scene, the same window, a different schedule. An evolving domain under `multirate`
    // takes as many substeps of the shared window as its own limit requires, so a tight limit is a
    // cost and not a refusal — and this check has nothing to say about cost.
    World::build(scene(&heated("multirate", 4.0, 8)))
        .expect("multirate subcycles; there is no step to refuse");
}

#[test]
fn every_shipped_scene_still_builds() {
    // **The regression this could have been.** `max_stable_dt` is state-dependent, so a domain
    // whose limit is tightest at `t = 0` and loosens as it runs would be refused at build for a
    // run it survives. Measured across all of them: none is.
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("scenes");
    let mut seen = 0;
    for entry in std::fs::read_dir(&dir).expect("the scenes are there") {
        let path = entry.expect("readable").path();
        if path.extension().and_then(|e| e.to_str()) != Some("json") {
            continue;
        }
        let text = std::fs::read_to_string(&path).expect("readable");
        let s: Scene = serde_json::from_str(&text).expect("it parses");
        World::build_with(s, &pantometry_world::Beside::of(&path))
            .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        seen += 1;
    }
    assert!(seen >= 30, "only {seen} scenes were built");
}

#[test]
fn a_scene_with_no_frames_is_left_to_the_other_checks() {
    // Zero frames is a division this must not do, and it is somebody else's refusal: the window is
    // undefined rather than too long. Checked here because a guard nobody exercises is a guard
    // that gets deleted.
    let json = heated("staggered", 4.0, 0);
    // Either it is refused for its own reason or it builds; what it must not do is divide by zero.
    match World::build(scene(&json)) {
        Ok(_) => {}
        Err(e) => assert!(
            !e.contains("staggered frame"),
            "zero frames was refused as a frame-window problem: {e}"
        ),
    }
}
