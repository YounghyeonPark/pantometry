//! A scene could say what its drive is and not when it changes.
//!
//! Every scene in this repository ran at a constant one, and the questions a design actually asks
//! are not constant: the junction temperature of a module under a duty cycle, a motor at start-up,
//! an espresso pulled with a pre-infusion. [`Scene::stages`] is a load profile, and this holds it
//! against the one number that makes it worth having — **the peak, which averaging the power does
//! not give you**.

use pantometry_world::{Scene, World};

fn scene(json: &str) -> Scene {
    serde_json::from_str(json).expect("the test scene parses")
}

/// A twenty-second run of a heater into a bar, with whatever profile is asked for.
fn heated(stages: &str, extra: &str) -> String {
    format!(
        r#"{{
  "title": "a heater into a bar", "duration_s": 20.0, "frames": 20{extra},
  "domains": [
    {{ "kind": "heater", "name": "element", "watts": 50.0, "reserve_j": 100000.0 }},
    {{ "kind": "bar", "name": "bar", "length_mm": 50.0, "cells": 20, "area_mm2": 25.0,
      "initial_c": 20.0 }}
  ]{stages}
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

/// **The same joules at a different rate are a different peak**, which is the whole reason a
/// profile is worth having.
///
/// A hundred watts for ten seconds and fifty for twenty are the same thousand joules. They leave
/// the bar at the same mean, because energy is energy. They do not leave it at the same peak,
/// because the peak is set by the thermal mass near the source and not by the average:
///
/// ```text
///   pulsed    peak 677.6442 C   final mean 350.6878 C   spent 1000 J
///   averaged  peak 521.4072 C   final mean 350.6878 C   spent 1000 J
/// ```
///
/// **30.0% under**, and a part sized on the average sees the smaller number. The two means agreeing
/// is the check that this is the same experiment in every other respect — if the profile had
/// changed the energy, they would not.
#[test]
fn the_same_joules_at_a_different_rate_are_a_different_peak() {
    let run = |text: &str| {
        let mut world = World::build(scene(text)).expect("it builds");
        let frames = world.run().expect("it conserves");
        let peak = frames
            .iter()
            .map(|f| read(f, "bar", "peak"))
            .fold(f64::MIN, f64::max);
        let last = frames.last().expect("frames");
        (
            peak,
            read(last, "bar", "mean"),
            100_000.0 - read(last, "element", "reserve"),
        )
    };

    let (pulsed, pulsed_mean, pulsed_spent) = run(&heated(
        r#","stages": { "element": [ { "at_s": 0.0, "watts": 100.0 },
                                    { "at_s": 10.0, "watts": 0.0 } ] }"#,
        "",
    ));
    let (flat, flat_mean, flat_spent) = run(&heated("", ""));

    // The same energy, exactly: 100 W for 10 s and 50 W for 20 s are both 1000 J, and the element
    // pays out of a tank that a stage does not touch.
    assert!(
        (pulsed_spent - 1000.0).abs() < 1e-9 && (flat_spent - 1000.0).abs() < 1e-9,
        "{pulsed_spent} J against {flat_spent} J, and both should be 1000"
    );
    // And therefore the same mean, which is the check that nothing else moved.
    assert!(
        (pulsed_mean / flat_mean - 1.0).abs() < 1e-9,
        "the same joules left different means: {pulsed_mean} against {flat_mean}"
    );
    // The peak is where they differ, and by a third.
    let over = pulsed / flat - 1.0;
    println!(
        "  pulsed peaks {:.4} C against {flat:.4} C, {:+.2}%",
        pulsed,
        over * 100.0
    );
    assert!(
        (0.25..0.35).contains(&over),
        "a pulse should peak about 30% over the average: {pulsed:.4} against {flat:.4}"
    );
}

/// **The profile is applied at the instant it names**, checked against `watts * elapsed` rather
/// than against a second run.
///
/// The element pays `watts` joules a second out of a tank, so what it has spent at any moment is
/// the integral of the profile — a closed form with no discretisation in it at all. A stage
/// applied one step late shows up here as exactly one step of the wrong power.
#[test]
fn a_profile_is_integrated_at_the_instants_it_names() {
    let text = heated(
        r#","stages": { "element": [ { "at_s": 0.0, "watts": 100.0 },
                                    { "at_s": 5.0, "watts": 0.0 },
                                    { "at_s": 15.0, "watts": 200.0 } ] }"#,
        "",
    );
    let mut world = World::build(scene(&text)).expect("it builds");
    let frames = world.run().expect("it conserves");

    // 100 W for 5 s, nothing for 10, 200 W for 5: 500, then 500, then 1500 J.
    for (frame, want) in [
        (0usize, 0.0),
        (5, 500.0),
        (10, 500.0),
        (15, 500.0),
        (20, 1500.0),
    ] {
        let spent = 100_000.0 - read(&frames[frame], "element", "reserve");
        assert!(
            (spent - want).abs() < 1e-9,
            "at frame {frame} ({} s) the element had spent {spent} J and the profile says {want}",
            frame as f64
        );
    }
}

/// **Both paths follow the profile**, because it hangs off `advance` and not off `run`.
///
/// The editor streams a scene by driving [`World::advance`] itself, and `verify` sweeps the same
/// way. Hanging a load profile off [`World::run`] would have given the batch path one experiment
/// and the streaming path another — which is the divergence `a_streamed_run_reads_back` already
/// exists for, one level up. This was written that way first and measured here.
#[test]
fn a_streamed_run_follows_the_same_profile() {
    let text = heated(
        r#","stages": { "element": [ { "at_s": 0.0, "watts": 100.0 },
                                    { "at_s": 10.0, "watts": 0.0 } ] }"#,
        "",
    );
    let mut batch = World::build(scene(&text)).expect("it builds");
    let frames = batch.run().expect("it conserves");
    let want = 100_000.0 - read(frames.last().expect("frames"), "element", "reserve");

    let mut streamed = World::build(scene(&text)).expect("it builds");
    let steps = streamed.steps();
    let dt = pantometry::units::Time::from_si(20.0 / steps as f64);
    for _ in 0..steps {
        streamed.advance(dt).expect("it conserves");
    }
    let got = 100_000.0
        - streamed
            .simulation()
            .domains()
            .flat_map(|d| d.readings())
            .find(|r| r.domain == "element" && r.label == "reserve")
            .expect("the element reports its reserve")
            .value;

    assert!(
        (got - want).abs() < 1e-9,
        "streamed spent {got} J and the batch run spent {want}"
    );
    assert!(
        (want - 1000.0).abs() < 1e-9,
        "and both should be 1000: {want}"
    );
}

/// **A step that passes several stages arrives at the last of them.**
///
/// `World::advance` takes whatever `dt` its caller hands it, and a caller is entitled to hand it
/// one coarser than the profile — a preview at ten frames of a scene written for a thousand. The
/// drives here are levels rather than increments, so passing three stages in one step means the
/// third; applying them in turn would put the same value there after two wasted solves, and
/// applying only the *first* unapplied one would leave the drive two stages behind and catch up
/// later, which is a profile the run does not follow.
///
/// Nothing here held that. The six profiles in this file all space their stages wider than the
/// step, so a sabotage that applied `stages[done]` instead of `stages[next - 1]` passed all of
/// them — which is the difference between a claim and a check.
#[test]
fn a_step_that_passes_several_stages_arrives_at_the_last() {
    let text = heated(
        r#","stages": { "element": [ { "at_s": 0.0, "watts": 0.0 },
                                    { "at_s": 5.0, "watts": 1000.0 },
                                    { "at_s": 10.0, "watts": 7.0 } ] }"#,
        "",
    );
    let mut world = World::build(scene(&text)).expect("it builds");

    // One step over the first two stages: 15 s at once, which passes 5.0 and 10.0 together.
    world
        .advance(pantometry::units::Time::from_si(15.0))
        .expect("it conserves");
    let after = 100_000.0
        - world
            .simulation()
            .domains()
            .flat_map(|d| d.readings())
            .find(|r| r.domain == "element" && r.label == "reserve")
            .expect("the element reports its reserve")
            .value;
    // The step ran at the drive in force when it began — 0 W, from the stage at 0 — so nothing
    // was spent, and the drive is now 7 W and not 1000.
    assert!(
        after.abs() < 1e-9,
        "the step should have run at 0 W and spent {after} J"
    );

    world
        .advance(pantometry::units::Time::from_si(5.0))
        .expect("it conserves");
    let spent = 100_000.0
        - world
            .simulation()
            .domains()
            .flat_map(|d| d.readings())
            .find(|r| r.domain == "element" && r.label == "reserve")
            .expect("the element reports its reserve")
            .value;
    assert!(
        (spent - 35.0).abs() < 1e-9,
        "5 s at the last stage reached, 7 W, is 35 J; this spent {spent}"
    );
}

/// **A stage between two steps is refused, and the refusal names a window that works.**
///
/// This is the only rule here that could be wrong quietly. A stage at 10.5 s on a 1 s step lands
/// at 11: the run heats half a second longer than the file says, the conservation audit closes
/// because the joules that were paid were taken, and every reading is a correct run of a
/// different experiment.
///
/// **The first version of this refusal named a window that did not work.** It derived one as
/// `duration / ceil(duration / at_s)` = 10 s, and `World::steps` raises the step count to `frames`
/// — so the step stayed 1 s and 10.5 still landed nowhere. The suggestion is searched now, and
/// this test feeds it back in.
#[test]
fn a_stage_between_steps_is_refused_and_the_window_it_names_works() {
    let off = heated(
        r#","stages": { "element": [ { "at_s": 10.5, "watts": 0.0 } ] }"#,
        "",
    );
    let why = World::build(scene(&off))
        .err()
        .expect("a stage between steps must be refused");
    assert!(why.contains("10.5"), "{why}");
    assert!(why.contains("it would happen at 11 s instead"), "{why}");

    // The window it names, taken out of the message rather than typed here — so a message that
    // changed its number would fail rather than leave this asserting a stale one.
    let quoted = why
        .split("a window_s of ")
        .nth(1)
        .and_then(|rest| rest.split(' ').next())
        .and_then(|n| n.parse::<f64>().ok())
        .unwrap_or_else(|| panic!("the refusal names no window: {why}"));
    let fixed = heated(
        r#","stages": { "element": [ { "at_s": 10.5, "watts": 0.0 } ] }"#,
        &format!(r#", "window_s": {quoted}"#),
    );
    let mut world = World::build(scene(&fixed))
        .unwrap_or_else(|e| panic!("the window the refusal named does not work: {e}"));
    world.run().expect("and it runs");

    // 100 W is the declared 50 until 10.5 s and 0 after: 525 J.
    let frames = World::build(scene(&fixed)).unwrap().run().unwrap();
    let spent = 100_000.0 - read(frames.last().unwrap(), "element", "reserve");
    assert!(
        (spent - 525.0).abs() < 1e-9,
        "50 W for 10.5 s is 525 J, and this spent {spent}"
    );
}

/// **The check and the run share one step count**, rather than agreeing about two.
///
/// `Scene::check_stages` needs a step count before a world exists, and writing the arithmetic
/// twice was the obvious thing to do. A drift between the copies would let a stage pass the check
/// and land off a step at run time — the exact failure the check is for, reintroduced by the
/// check itself. So `World::steps` delegates to `Scene::steps`, and this holds the delegation:
/// every combination below reaches the same number *through the world*, which is the caller that
/// could have kept its own.
#[test]
fn the_check_and_the_run_share_one_step_count() {
    for (frames, window, want) in [
        (20usize, None, 20usize),
        (20, Some(0.5), 40),
        (20, Some(2.0), 20),
        (7, Some(0.3), 67),
        (1, Some(100.0), 1),
        (5, Some(1e-3), 20000),
    ] {
        let extra = match window {
            Some(w) => format!(r#", "window_s": {w}"#),
            None => String::new(),
        };
        let text = format!(
            r#"{{
  "title": "steps", "duration_s": 20.0, "frames": {frames}{extra},
  "domains": [
    {{ "kind": "heater", "name": "element", "watts": 50.0, "reserve_j": 100000.0 }},
    {{ "kind": "bar", "name": "bar", "length_mm": 50.0, "cells": 20, "area_mm2": 25.0,
      "initial_c": 20.0 }}
  ],
  "stages": {{ "element": [ {{ "at_s": 0.0, "watts": 10.0 }} ] }}
}}"#
        );
        let parsed = scene(&text);
        assert_eq!(parsed.steps(), want, "frames {frames}, window {window:?}");
        let world = World::build(parsed).expect("it builds");
        assert_eq!(world.steps(), want, "frames {frames}, window {window:?}");
    }
}

/// Every way a profile can be wrong says which, and none of them is a silent no-op.
///
/// A stage that quietly did nothing is a load profile that reads as applied and is not, which is
/// this format's oldest failure shape. Each of these is refused at build with the domain named.
#[test]
fn a_profile_that_cannot_be_followed_says_why() {
    let cases: [(&str, &str); 7] = [
        (
            r#","stages": { "elment": [ { "at_s": 10.0, "watts": 0.0 } ] }"#,
            "which this scene does not define",
        ),
        (
            r#","stages": { "bar": [ { "at_s": 10.0, "watts": 0.0 } ] }"#,
            "no drive to change",
        ),
        (
            r#","stages": { "element": [ { "at_s": 10.0, "volts": 0.0 } ] }"#,
            "this kind's drive is \"watts\"",
        ),
        (
            r#","stages": { "element": [ { "at_s": 10.0, "watts": 1.0, "volts": 0.0 } ] }"#,
            "Each sets exactly one",
        ),
        (
            r#","stages": { "element": [ { "at_s": 10.0 } ] }"#,
            "Each sets exactly one",
        ),
        (
            r#","stages": { "element": [ { "at_s": 20.0, "watts": 0.0 } ] }"#,
            "would never happen",
        ),
        (
            r#","stages": { "element": [] }"#,
            "a profile that says nothing",
        ),
    ];
    for (stages, says) in cases {
        let why = World::build(scene(&heated(stages, "")))
            .err()
            .unwrap_or_else(|| panic!("this should be refused: {stages}"));
        assert!(why.contains(says), "expected {says:?} in: {why}");
    }

    // Out of order, which needs two stages.
    let why = World::build(scene(&heated(
        r#","stages": { "element": [ { "at_s": 10.0, "watts": 0.0 },
                                    { "at_s": 5.0, "watts": 1.0 } ] }"#,
        "",
    )))
    .err()
    .expect("out of order must be refused");
    assert!(why.contains("in the order written"), "{why}");

    // A misspelled drive is serde's to catch, and it lists what it expected — which is the
    // difference between a typo and a key that silently does nothing.
    let e = serde_json::from_str::<Scene>(&heated(
        r#","stages": { "element": [ { "at_s": 10.0, "wats": 0.0 } ] }"#,
        "",
    ))
    .expect_err("a misspelled drive must not parse");
    assert!(e.to_string().contains("unknown field `wats`"), "{e}");
}
