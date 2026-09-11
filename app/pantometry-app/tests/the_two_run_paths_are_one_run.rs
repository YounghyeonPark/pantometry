//! **The editor's streaming path and the batch call are supposed to be one run.**
//!
//! `run_streaming`'s own documentation says so — "the last payload is byte-identical to what
//! [`run`] returns for the same text, which the tests pin" — and nothing pinned it. The only test
//! that touched the streaming path checked that its JSON *reads back*, which is a different claim
//! and the one `a_streamed_run_reads_back` exists for.
//!
//! What the gap hid: `run_streaming` computed its own step as `duration_s / frames`, ignoring
//! [`Scene::window_s`](pantometry_world::Scene::window_s). That key exists precisely so the step
//! can be shorter than a frame, so any scene using it for what it is for ran one experiment in the
//! CLI and a different one in the editor. No shipped scene did — all thirty have `steps == frames`
//! — which is why nothing said.

#![cfg(not(target_family = "wasm"))]

use pantometry_world::OnDisk;

/// A heater into a bar, with a window four times shorter than a frame.
///
/// `steps` is `ceil(20 / 0.25)` = 80 against 20 frames, so the batch path takes 80 steps of 0.25 s
/// and the streaming path used to take 20 of 1 s. Four times the step is a different answer to the
/// same file.
const WINDOWED: &str = r#"{
  "title": "a window shorter than a frame",
  "duration_s": 20.0,
  "frames": 20,
  "window_s": 0.25,
  "domains": [
    { "kind": "heater", "name": "element", "watts": 50.0, "reserve_j": 100000.0 },
    { "kind": "bar", "name": "bar", "length_mm": 50.0, "cells": 20, "area_mm2": 25.0,
      "initial_c": 20.0 }
  ]
}"#;

fn last_streamed(text: &str) -> String {
    let stop = std::sync::atomic::AtomicBool::new(false);
    let mut last = String::new();
    editor_core::run_streaming(text, &OnDisk, &stop, |json| last = json)
        .expect("the streaming run finishes");
    last
}

/// **The two paths produce the same bytes**, on a scene whose window is shorter than its frame.
#[test]
fn a_windowed_scene_streams_what_it_runs() {
    let batch = editor_core::run(WINDOWED, &OnDisk).expect("the batch run finishes");
    let streamed = last_streamed(WINDOWED);
    assert_eq!(
        batch.len(),
        streamed.len(),
        "the two paths produced runs of different length"
    );
    assert!(
        batch == streamed,
        "the batch run and the streaming run differ on a scene with window_s"
    );
}

/// **And on a scene whose steps are not a whole number of frames**, which is what pins the
/// *schedule* rather than only the step.
///
/// 20 steps photographed at 7 frames: frame `i` lands after `ceil(i * 20 / 7)` steps — 3, 6, 9, 12,
/// 15, 18, 20 — and plain integer division would give 2, 5, 8, 11, 14, 17, 20, every picture one
/// step early and the last one right. The windowed scene above has 80 steps over 20 frames, an
/// exact four, so the two roundings agree there and a sabotage that swapped them passed all three
/// tests. This is the case that tells them apart.
#[test]
fn a_schedule_that_does_not_divide_streams_what_it_runs() {
    let odd = WINDOWED
        .replace("\"frames\": 20,", "\"frames\": 7,")
        .replace("\"window_s\": 0.25,", "\"window_s\": 1.0,");
    assert!(
        odd.contains("\"frames\": 7,"),
        "the frame count was not changed"
    );
    let batch = editor_core::run(&odd, &OnDisk).expect("the batch run finishes");
    let streamed = last_streamed(&odd);
    assert_eq!(
        batch, streamed,
        "20 steps over 7 frames: the two paths photographed different steps"
    );
}

/// And the same on a scene with no window at all, which is the case that always worked — kept so a
/// fix aimed at the windowed case cannot break the ordinary one.
#[test]
fn a_plain_scene_streams_what_it_runs() {
    let plain = WINDOWED.replace("  \"window_s\": 0.25,\n", "");
    assert!(!plain.contains("window_s"), "the window was not removed");
    let batch = editor_core::run(&plain, &OnDisk).expect("the batch run finishes");
    assert_eq!(batch, last_streamed(&plain));
}

/// **The step is the scene's, not the frame count's**, measured rather than inferred.
///
/// The bar's mean after 20 s is the same either way — energy is energy — so the number that
/// separates the two schedules is the *frame* the run is photographed at. With a 0.25 s step the
/// first frame lands after four steps; with a 1 s step it lands after one, and a diffusing bar has
/// a different peak at those two instants.
#[test]
fn the_window_changes_what_the_first_frame_shows() {
    let windowed: viewer_core::Run =
        serde_json::from_str(&last_streamed(WINDOWED)).expect("it reads back");
    let plain_text = WINDOWED.replace("  \"window_s\": 0.25,\n", "");
    let plain: viewer_core::Run =
        serde_json::from_str(&last_streamed(&plain_text)).expect("it reads back");

    let peak = |run: &viewer_core::Run, i: usize| {
        run.frames[i]
            .readings
            .iter()
            .find(|r| r.domain == "bar" && r.label == "peak")
            .expect("the bar reports a peak")
            .value
    };
    // Both schedules photograph the same instants, so the *last* frame agrees to the step's own
    // accuracy; the first is where a four-times-shorter step shows.
    let (a, b) = (peak(&windowed, 1), peak(&plain, 1));
    println!("  first frame: {a:.6} C windowed against {b:.6} C plain");
    assert!(
        (a - b).abs() > 1e-6,
        "a window four times shorter than the frame should move the first frame: {a} against {b}"
    );
}
