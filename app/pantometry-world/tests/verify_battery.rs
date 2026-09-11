//! The `verify` battery, held to its own standard.
//!
//! The battery's whole claim is that it measures rather than asserts, so its tests are the
//! places a measurement has a known answer: a refinement that must preserve a total exactly, a
//! window error whose order the kernel documentation derives and measures, a mode whose
//! convergence rate is the scheme's own second order, and a hazard the schedule documentation
//! states in words that a test can arrange in numbers.

use pantometry_world::verify::{verify, Order, SweepOutcome};
use pantometry_world::{DomainSpec, Scene};

fn scene(json: &str) -> Scene {
    serde_json::from_str(json).expect("the test scene parses")
}

/// The sweep behind an outcome, unwrapped where the test knows it must have run.
fn ran(outcome: &SweepOutcome) -> &pantometry_world::verify::Sweep {
    match outcome {
        SweepOutcome::Ran(s) => s,
        SweepOutcome::Skipped(why) => panic!("the sweep was skipped: {why}"),
        SweepOutcome::Failed(why) => panic!("the sweep failed: {why}"),
    }
}

/// The order a sweep measured for one reading, unwrapped where the test knows it must exist.
fn order_of(
    battery: &pantometry_world::verify::Battery,
    which: &str,
    domain: &str,
    label: &str,
) -> f64 {
    let sweep = ran(match which {
        "window" => &battery.window,
        _ => &battery.resolution,
    });
    let (_, _, order) = sweep
        .orders
        .iter()
        .find(|(d, l, _)| d == domain && l == label)
        .expect("the reading was matched across the runs");
    match order {
        Order::Measured(p) => *p,
        Order::BelowFloor => {
            panic!("{domain}/{label}: at the floor, where this test expected a measurable order")
        }
        Order::NotAsymptotic => {
            panic!("{domain}/{label}: not asymptotic, where this test expected a measurable order")
        }
    }
}

/// **Two runs of one scene produce identical bytes, and the battery checks it rather than
/// assuming it.** This is the determinism promise at the scale of a whole run — every panel,
/// body and reading of every frame — not just the pinned generator stream.
#[test]
fn two_runs_are_bit_identical_and_the_battery_says_so() {
    let s = scene(
        r#"{
  "title": "a small room, checked",
  "duration_s": 0.005,
  "frames": 3,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 2.0, "height_m": 1.0,
      "cells_across": 21,
      "release": { "as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0 } }
  ]
}"#,
    );
    let b = verify(&s, false).expect("the battery runs");
    assert!(b.deterministic, "one scene, two different runs");
    assert!(
        b.findings.is_empty(),
        "findings on a clean scene: {:?}",
        b.findings
    );
    // Both sweeps ran and actually compared something — a sweep with zero shifts would be
    // indistinguishable in the report's shape from "swept and nothing moved".
    assert!(!ran(&b.window).shifts.is_empty());
    assert!(!ran(&b.resolution).shifts.is_empty());
    // The margin rows exist: energy was on the ledger and the room reported a limit.
    assert!(b.base.drift.iter().any(|(q, _, _)| q == "energy"));
    assert!(b.base.stability.iter().any(|(n, _, _)| n == "room"));
}

/// **A scene that reports no readings is a finding, not a verification.**
///
/// Neither mechanics domain implements `Domain::readings`, so an orbit scene gives the sweeps
/// nothing to compare. The battery used to run the sweeps anyway, compare zero numbers, print
/// a heading with nothing under it and exit 0 — a report whose shape is identical to "swept
/// and nothing moved". A user would conclude the orbit is insensitive to the coupling window;
/// in fact nothing was measured. Now: both sweeps are skipped with the reason stated, the
/// report says which domain was not swept, and the empty measurement is a finding that
/// carries the exit code.
#[test]
fn a_scene_with_no_readings_is_a_finding_not_a_verification() {
    let s = scene(
        r#"{
  "title": "an orbit nobody reads",
  "duration_s": 4.0,
  "frames": 4,
  "domains": [
    { "kind": "orbit", "name": "orbit", "central_kg": 5e24,
      "radii_m": [7e6], "satellite_kg": 1000.0 }
  ]
}"#,
    );
    let b = verify(&s, false).expect("the run itself is fine");
    assert!(
        b.findings
            .iter()
            .any(|f| f.contains("nothing here was verified")),
        "an empty measurement must be a finding: {:?}",
        b.findings
    );
    assert!(matches!(b.window, SweepOutcome::Skipped(_)));
    assert!(b.unread.contains(&"orbit".to_string()));
    let report = b.render();
    assert!(
        report.contains("SKIPPED: no domain reports a reading"),
        "the report does not say what was not measured:\n{report}"
    );
}

/// **A sweep run the audit refused is a FAILED finding, not a benign skip.**
///
/// A staggered room just under its CFL limit at base resolution is over it at 2×, because
/// refining halves the limit. The refined run is refused on its first advance — and a scene
/// that passes at its own settings and fails when a knob moves is one knob away from wrong,
/// which is about the most alarming fact this battery can learn. It used to render as
/// `SKIPPED:` and exit 0, the same reading as "not applicable".
#[test]
fn a_sweep_run_the_audit_refused_is_a_finding_not_a_skip() {
    let s = scene(
        r#"{
  "title": "stable at its own resolution and not at the next",
  "schedule": "staggered",
  "duration_s": 0.00045,
  "frames": 3,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 4.0, "height_m": 2.0,
      "cells_across": 41,
      "release": { "as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0 } }
  ]
}"#,
    );
    // Base window 1.5e-4 s sits under the 41-sample limit of 2.06e-4 s, so the base run and
    // the hazard check both pass; at 81 samples the limit is 1.03e-4 s and the refined run is
    // refused.
    let b = verify(&s, false).expect("the base scene runs");
    assert!(matches!(b.resolution, SweepOutcome::Failed(_)));
    assert!(
        b.findings
            .iter()
            .any(|f| f.contains("one knob away from wrong")),
        "a refused sweep run must be a finding: {:?}",
        b.findings
    );
    let report = b.render();
    assert!(
        report.contains("FAILED:"),
        "not rendered as a failure:\n{report}"
    );
}

/// **The report states what it did not measure, row by row.**
///
/// A mixed scene: the room is swept and the orbit — which reports no readings — must appear
/// as "no readings — not swept" rather than being simply absent, because a sweep section
/// listing every other domain looks complete. And a flat orbit's `momentum_z` never rises
/// above the audit's scale floor, so it is never audited; the margins table must say that
/// rather than omit the row, because a table of quantities-with-margin that quietly drops one
/// is the shape of a check that turned itself off.
#[test]
fn the_report_states_what_it_did_not_measure() {
    let s = scene(
        r#"{
  "title": "a room and an unread orbit",
  "duration_s": 0.005,
  "frames": 3,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 2.0, "height_m": 1.0,
      "cells_across": 21,
      "release": { "as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0 } },
    { "kind": "orbit", "name": "orbit", "central_kg": 5e24,
      "radii_m": [7e6], "satellite_kg": 1000.0 }
  ]
}"#,
    );
    let b = verify(&s, false).expect("the battery runs");
    let report = b.render();
    assert!(
        report.contains("orbit          no readings — not swept"),
        "the unread domain is invisible:\n{report}"
    );
    assert!(
        report.contains("never above the scale floor — not audited"),
        "the unaudited quantity is invisible:\n{report}"
    );
    // And the room was genuinely swept beside it.
    assert!(!ran(&b.window).shifts.is_empty());
}

/// **A staggered window past a domain's stability limit is refused before anything is built.**
///
/// `ScheduleSpec::Staggered`'s own documentation: a frame interval larger than the limit is
/// silently unstable. The room turns out not to be silent — its step reports the Courant number as
/// a created quantity and the audit refuses the very first advance — but a refusal that only says
/// "Courant number created" leaves the *why* to be worked out.
///
/// This battery used to supply the why. It asked a **built, unrun** world for every domain's
/// `max_stable_dt` and appended the answer to a failed run as "likely why", so its error carried
/// both halves: what the kernel refused, and the schedule choice that made it inevitable.
///
/// **`World::build` refuses that scene now**, so there is no run to explain and the reader gets the
/// cause on its own. `FRICTION.md`'s finding 8 asked for exactly that — a check at build time,
/// where the message can name the file — and half of it turned out to be sitting in this file
/// already, reachable only by somebody who ran `verify` on purpose. Two implementations of one
/// rule became one, and it is the one `--check` and the editor reach while somebody is typing.
///
/// The wording is this battery's, kept deliberately: a reader who knows "does not subcycle" and
/// "silently unstable" meets the same words in the new place.
#[test]
fn a_staggered_window_past_the_limit_is_named_as_the_likely_cause() {
    let unstable = r#"{
  "title": "one step across the whole window",
  "schedule": "staggered",
  "duration_s": 0.002,
  "frames": 3,
  "conservation_tolerance": 1e6,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 4.0, "height_m": 2.0,
      "cells_across": 41,
      "release": { "as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0 } }
  ]
}"#;
    let why = verify(&scene(unstable), false).expect_err("the room cannot be built at that window");
    assert!(
        why.contains("does not subcycle") && why.contains("silently unstable"),
        "the error does not name the schedule hazard: {why}"
    );
    // And it names the domain, both times, and what to change — which "likely why" appended to a
    // crash never did. The battery is passing the build's refusal through, not writing its own.
    assert!(why.starts_with("room:"), "{why}");
    assert!(why.contains("raise `frames`"), "{why}");

    // The same room under multirate subcycles, and the same battery finds nothing.
    let subcycled = unstable.replace("\"staggered\"", "\"multirate\"");
    let b = verify(&scene(&subcycled), false).expect("multirate subcycles the same window");
    assert!(
        !b.findings.iter().any(|f| f.contains("silently unstable")),
        "multirate flagged as unstable: {:?}",
        b.findings
    );
}

/// **A refinement preserves the problem, exactly where exactness is checkable.**
///
/// A bar doubles its cells and halves each face; the beam that lands on it doubles its faces
/// and halves each area. The exposed boundary's total area — `cells × face_area` — is the
/// number both sides must preserve to the bit, because the kernel refuses a face-count
/// mismatch and conserves the flux across whatever count both sides agree on.
#[test]
fn a_bar_and_its_beam_refine_together() {
    let s = scene(
        r#"{
  "title": "a beam on a bar",
  "duration_s": 1.0,
  "frames": 2,
  "domains": [
    { "kind": "beam", "name": "beam", "onto": "face", "faces": 10,
      "face_area_mm2": 2.0, "watts": 5.0, "reserve_j": 100.0, "waist_fraction": 0.2 },
    { "kind": "bar", "name": "bar", "length_mm": 20.0, "cells": 10, "area_mm2": 4.0,
      "initial_c": 20.0, "exposes": { "name": "face", "face_area_mm2": 2.0 } }
  ]
}"#,
    );
    let fine = s.refined().expect("a bar and a beam both refine");
    let mut checked = 0;
    for spec in &fine.domains {
        match spec {
            DomainSpec::Beam {
                faces,
                face_area_mm2,
                watts,
                ..
            } => {
                assert_eq!(*faces, 20);
                assert_eq!(*face_area_mm2, 1.0);
                assert_eq!(*watts, 5.0, "refinement moved a physical parameter");
                checked += 1;
            }
            DomainSpec::Bar { cells, exposes, .. } => {
                assert_eq!(*cells, 20);
                let b = exposes.as_ref().expect("the boundary survives refinement");
                assert_eq!(b.face_area_mm2, 1.0);
                checked += 1;
            }
            other => panic!("an unexpected domain appeared: {other:?}"),
        }
    }
    assert_eq!(checked, 2);
}

/// **A room refines by intervals, not by count.** Its samples sit on the walls, so `n` samples
/// are `n − 1` intervals and halving the spacing is `2n − 1`. `2n` would be a spacing 60/121
/// of the old one — near half, not half — and the residue would be read as physics.
#[test]
fn a_room_refines_by_intervals() {
    let s = scene(
        r#"{
  "title": "a room",
  "duration_s": 0.01,
  "frames": 2,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 4.0, "height_m": 2.0,
      "cells_across": 41 }
  ]
}"#,
    );
    let fine = s.refined().expect("a room refines");
    match &fine.domains[0] {
        DomainSpec::Room { cells_across, .. } => assert_eq!(*cells_across, 81),
        other => panic!("{other:?}"),
    }
}

/// **What cannot be refined honestly refuses, naming the feature.** A hot spot is one cell of
/// excess temperature: at half the cell it holds an eighth of the joules, so a "refined" scene
/// would be a different experiment whose difference the sweep would report as discretisation
/// error. The same for a puck's one-cell channel ring.
///
/// A conductor's blocked cells were on this list for the same stated reason and did not belong:
/// see [`a_notch_keeps_its_millimetres`], which is the measurement that reason was costing.
#[test]
fn one_cell_features_refuse_to_refine() {
    let hot = scene(
        r#"{
  "title": "a warmed cell",
  "duration_s": 1.0,
  "frames": 2,
  "domains": [
    { "kind": "block", "name": "block", "cells": [4, 4, 4], "cell_mm": 2.0,
      "initial_c": 20.0, "hot_spot": { "at": [2, 2, 2], "above_k": 10.0 } }
  ]
}"#,
    );
    let why = hot.refined().expect_err("a hot spot must refuse");
    assert!(
        why.contains("hot spot"),
        "the refusal does not say why: {why}"
    );

    let channelled = scene(
        r#"{
  "title": "a channelled puck",
  "duration_s": 5.0,
  "frames": 2,
  "domains": [
    { "kind": "puck", "name": "puck", "cells": [9, 6, 9], "cell_mm": 2.0,
      "radius_mm": 7.0, "grind_um": 250.0, "porosity": 0.45, "bar": 9.0,
      "brew_c": 93.0, "channel_porosity": 0.6 }
  ]
}"#,
    );
    let why = channelled
        .refined()
        .expect_err("a channel ring must refuse");
    assert!(why.contains("channel"), "{why}");
}

/// **The battery measures the run the CLI performs**, which its loop's documentation has claimed
/// since it was written and nothing held.
///
/// `run_measured`'s doc says "the loop is `World::run`'s". It has drifted from that twice: it
/// computed its own `duration / frames` until `window_s` existed, and it reached past
/// [`World::advance`] to the `Simulation` inside — which steps the domains and does none of the
/// other four things an advance is. A [`stages`](pantometry_world::Scene::stages) profile was
/// never applied, and a structure was never handed its thermal strain or re-solved.
///
/// Measured on a motor with a start-up load: the battery reported its winding at **44.80 °C**
/// where the run ends at **59.64**. The element had run at its start-up 36 W for the whole run,
/// emptied its tank at 800 s and cooled from there — an entirely self-consistent run of a scene
/// nobody wrote, with every audit margin healthy and no finding raised.
///
/// So this compares the two, reading by reading. It is the pin that makes the sentence in that doc
/// comment cost something.
#[test]
fn the_battery_measures_the_run_the_world_performs() {
    let text = r#"{
  "title": "a motor that starts under load",
  "schedule": "multirate",
  "duration_s": 1800.0,
  "frames": 12,
  "window_s": 150.0,
  "conservation_tolerance": 1e-9,
  "stages": { "losses": [ { "at_s": 300.0, "watts": 12.0 } ] },
  "domains": [
    { "kind": "heater", "name": "losses", "watts": 36.0, "reserve_j": 28800.0 },
    { "kind": "network", "name": "motor", "absorbing": "winding",
      "nodes": [
        { "name": "winding", "material": "copper", "volume_cm3": 18.0,
          "thickness_mm": 2.0, "initial_c": 25.0 },
        { "name": "stator", "material": "electrical_steel", "volume_cm3": 140.0,
          "thickness_mm": 8.0, "initial_c": 25.0 },
        { "name": "housing", "material": "aluminium", "volume_cm3": 220.0,
          "thickness_mm": 4.0, "initial_c": 25.0,
          "loses_to": { "ambient_c": 25.0, "area_cm2": 420.0 } }
      ],
      "links": [
        { "from": "winding", "to": "stator", "w_per_k": 0.9 },
        { "from": "stator", "to": "housing", "w_per_k": 2.4 }
      ] }
  ]
}"#;
    let s = scene(text);
    let mut world = pantometry_world::World::build(s.clone()).expect("it builds");
    let frames = world.run().expect("it conserves");
    let ran = &frames.last().expect("frames").readings;

    let b = verify(&s, false).expect("the battery runs");
    assert!(!b.base.readings.is_empty(), "the battery measured nothing");
    assert_eq!(
        b.base.readings.len(),
        ran.len(),
        "the battery and the run report different readings"
    );
    for measured in &b.base.readings {
        let same = ran
            .iter()
            .find(|r| r.domain == measured.domain && r.label == measured.label)
            .unwrap_or_else(|| panic!("the run has no {}/{}", measured.domain, measured.label));
        assert_eq!(
            measured.value, same.value,
            "{}/{}: the battery measured {} and the run gives {}",
            measured.domain, measured.label, measured.value, same.value
        );
    }

    // And the number that made it visible, so a reader sees what was at stake rather than only
    // that two vectors matched.
    let winding = ran
        .iter()
        .find(|r| r.domain == "motor" && r.label == "winding")
        .expect("the motor reports its winding");
    println!("  the winding ends at {:.4} C, in both", winding.value);
    assert!(
        (winding.value - 59.6374).abs() < 1e-3,
        "the profile was not followed: {:.4} C",
        winding.value
    );
}

/// **A blocked cell is a box, and a box refines by keeping its bounds.**
///
/// This refused, and the refusal was reasoned: "a blocked cell is a one-cell notch, so refining
/// shrinks it — a different geometry, not a finer one". That is true of a refinement that keeps
/// the *indices* and of no other. A cell at `(4, 0, 0)` on a 1 mm grid occupies `x ∈ [4, 5] mm`;
/// at 0.5 mm those same millimetres are indices 8 and 9. Each blocked cell becomes its **eight
/// children**, and the notch is the same notch.
///
/// What is asserted here is the millimetres, not the indices: a count of eight is not the claim,
/// and the sabotage that showed why was one that *passed*. Permuting which bit of `n` feeds which
/// axis looks like doubling the wrong corner and is a no-op — `0..8` enumerates all three bits
/// either way, so the same eight children come out in a different order. What does break it is a
/// box that lands half a cell off, an axis left undoubled, or a child count of four; all three
/// fail on the bounds and would pass a count.
#[test]
fn a_notch_keeps_its_millimetres() {
    let notched = scene(
        r#"{
  "title": "a notched conductor",
  "duration_s": 1.0,
  "frames": 2,
  "domains": [
    { "kind": "conductor", "name": "busbar", "cells": [8, 4, 2], "cell_mm": 1.0,
      "resistivity_ohm_m": 1.724e-8, "volts": 0.01, "blocked": [[4, 0, 0], [4, 1, 0]] }
  ]
}"#,
    );

    // The hole the coarse scene cuts, in millimetres: two cells stacked in y, one deep in z.
    let bounds = |spec: &DomainSpec| match spec {
        DomainSpec::Conductor {
            cells,
            cell_mm,
            blocked,
            ..
        } => {
            let lo = |axis: usize| blocked.iter().map(|c| c[axis]).min().unwrap() as f64 * cell_mm;
            let hi =
                |axis: usize| (blocked.iter().map(|c| c[axis]).max().unwrap() + 1) as f64 * cell_mm;
            (
                [lo(0), lo(1), lo(2)],
                [hi(0), hi(1), hi(2)],
                [
                    cells[0] as f64 * cell_mm,
                    cells[1] as f64 * cell_mm,
                    cells[2] as f64 * cell_mm,
                ],
                blocked.len(),
            )
        }
        other => panic!("an unexpected domain appeared: {other:?}"),
    };

    let (lo, hi, bar, count) = bounds(&notched.domains[0]);
    assert_eq!(lo, [4.0, 0.0, 0.0]);
    assert_eq!(hi, [5.0, 2.0, 1.0]);
    assert_eq!(bar, [8.0, 4.0, 2.0]);
    assert_eq!(count, 2);

    let fine = notched.refined().expect("a blocked cell refines");
    let (flo, fhi, fbar, fcount) = bounds(&fine.domains[0]);
    assert_eq!(flo, lo, "the notch moved");
    assert_eq!(fhi, hi, "the notch changed size");
    assert_eq!(fbar, bar, "the bar changed size");
    assert_eq!(fcount, count * 8, "a box has eight children");

    // And the ligament — the metal the current has to squeeze through — is the same metal.
    // This is the number the whole scene is about, so it is asserted from the far side of the
    // hole rather than inferred from the two above.
    let ligament = |spec: &DomainSpec| match spec {
        DomainSpec::Conductor {
            cells,
            cell_mm,
            blocked,
            ..
        } => (cells[1] - (blocked.iter().map(|c| c[1]).max().unwrap() + 1)) as f64 * cell_mm,
        other => panic!("{other:?}"),
    };
    assert_eq!(ligament(&notched.domains[0]), 2.0);
    assert_eq!(ligament(&fine.domains[0]), 2.0);
}

/// **A notch converges at the corner's own exponent, and it is not two.**
///
/// The sweep would not run on `17-a-busbar-with-a-notch` at all, so nobody could ask how much of
/// its resistance was grid. It is **4.65%** at the shipped 1 mm, and the reason is not resolution
/// but the shape.
///
/// Current crowds into a re-entrant corner the way stress does. The conducting wedge at a notch
/// root has an interior angle of `3π/2`, so the potential goes as `r^λ` with `λ = π/ω = 2/3`
/// and the flux as `r^(λ-1)` — unbounded, at every grid. A resistance is a **quadratic**
/// functional of that field, so its error goes as `h^(2λ) = h^(4/3)` and not as the `h²` a smooth
/// field gives.
///
/// ```text
///      h        in-plane     R at t = 5 mm     above the limit    three-point order
///    1 mm        12 x 5      1.2391920e-05        +4.6525%
///    0.5         24 x 10     1.2059472e-05        +1.8449%           1.33534
///    0.25        48 x 20     1.1927724e-05        +0.7323%           1.33296
///    0.125       96 x 40     1.1875426e-05        +0.2906%           1.33318
///    0.0625     192 x 80     1.1854670e-05        +0.1153%           1.33367
///    0.03125    384 x 160    1.1846434e-05        +0.0458%
/// ```
///
/// Read from the ten-digit CSV, not the seven-digit JSON: by `0.0625` the successive differences
/// are 1e-5 of the value and seven digits stop resolving an order.
///
/// **One cell through the thickness**, because `R·t` is exactly independent of it —
/// `6.195960000000e-08` at `nz` = 1, 2, 3, 5 and 8, to twelve digits, the slot spanning the full
/// `z` and the geometry being prismatic. That is what makes the table affordable: `0.0625 mm`
/// costs **1.1 s** this way against the **201 s** the same plane took at `nz = 5`.
///
/// **Which way the error points is the part that was wrong here.** This said "a resistance is an
/// energy-norm quantity", which is the conforming-Galerkin identity `E_h - E = ||u-u_h||²_E` —
/// and by Dirichlet's principle that identity makes a voltage-driven resistor read **low**. Every
/// row above reads high. The crate is not a Galerkin method: it is a cell-centred conductance
/// network with harmonic-mean faces and `Σ g Δφ²`, whose face currents are exactly conservative
/// per cell. So the argument is the **dual** one — a conservative current field is admissible for
/// Thomson's principle, which makes every grid an *upper* bound on the resistance, and the lumped
/// `Σ g Δφ²` sits above the reconstruction's own quadrature by
/// `(a²+b²)/2 - (a²+ab+b²)/3 = (a-b)²/6 ≥ 0`, cell by cell. What is *measured* is the sign, at
/// all six grids; the dual argument is why it has to be that sign.
///
/// **The four orders bracket 4/3 rather than descend to it.** 1.33534, 1.33296, 1.33318, 1.33367
/// against 1.33333: the widest is 0.15% out and the other three are inside 0.03%, but the sequence
/// is **not monotone**, and quoting the first three as a converging one reads tighter than it is.
/// The spread is the subleading term, which at these grids is comparable to what is being
/// measured; it is not uncertainty in the exponent. A free three-parameter fit `R₀ + C·h^p` over
/// the finest five rows puts `p` in **[1.3323, 1.3339]** — 4/3 to within 0.08% — taking the
/// interval where the residual stays inside ten times its minimum. At fixed `p` over all six, 4/3
/// fits **186x** better than 1.25, **335x** better than 1.5, **803x** better than 1.0 and **1131x**
/// better than 2.0. The same fit puts the limit at `1.1841011e-05`, against `1.1841016e-05` from
/// Richardson on the finest pair: seven digits, as the pin says.
///
/// The slot is **one cell wide** at the shipped grid — 1 mm, `x ∈ [6, 7] mm` — so at 1x the two
/// notch roots are exactly `h` apart and the corner is not resolved at all. That is where the
/// 4.65% is.
///
/// The percentages are against the Richardson limit at `p = 4/3`. The finest three pairs give
/// 1.1841016, 1.1841012 and 1.1841016e-05 ohm, so **1.184101e-05** is what is earned and the
/// eighth digit is not. `1.1841016e-05` is the value the pin uses, from the finest pair.
///
/// **And the grid runs out before the answer does, on the audit rather than the clock.** At
/// `0.015625 mm` the scene is *refused*: its own drift budget of `1e-9` reads `1.094e-9` after
/// 95 s, floating-point noise over 246 k cells crossing a tolerance written for a scene a
/// thousandth the size. So the finest statement this scene can make about itself still has
/// **0.046%** in it, and that is the useful thing to know about a notch — the error is not
/// waiting for a finer grid, it is the corner, and `h^(4/3)` is how slowly it goes.
///
/// The tolerance is 1%. The order asserted is the coarsest triplet's, 1.33534, which is 0.15%
/// out — so 1% is six times the deviation and it excludes 1.25 (6.3% away), 1.5 (12.5%), 1.0 (a
/// crack, 25%) and 2.0 (50%). It was 4%, which excluded the same four and said less.
///
/// **The scene is the shipped file, not a copy of it.** This was written as inlined JSON that
/// restated all fifteen blocked cells, and nothing compared the two — so the sentence "17
/// converges at `h^(4/3)`", which four documents now carry, was held by a test measuring a
/// different `Scene` value. `include_str!` cannot drift. It also makes this the only test that
/// puts a *shipped* scene through [`verify`]; the walk in `scene.rs` builds and runs them, and
/// the refinement path this exercises is reached by nothing else.
#[test]
fn a_notch_converges_at_the_corner_exponent() {
    let notched = scene(include_str!("../scenes/17-a-busbar-with-a-notch.json"));
    let b = verify(&notched, true).expect("the notched scene runs at three grids");
    let order = order_of(&b, "resolution", "busbar", "resistance");
    let corner = 4.0 / 3.0;
    assert!(
        (order / corner - 1.0).abs() < 0.01,
        "a re-entrant corner gives 2*lambda = {corner:.4}; the sweep measured {order:.4}"
    );

    // **And the row that started all this says what it is.** The defect was
    // `residual 0.000000 -> 0.000000 (281492.026%)`: two converged residuals printed as zero with
    // a percentage of noise between them. The branch that replaced it was reached by one
    // unrelated test and asserted by nothing, so the text and the format — which are the whole
    // change — were held by nothing either. The negative is the one that would have caught the
    // original row.
    let report = b.render();
    assert!(
        report.contains("a solver diagnostic, not a discretisation"),
        "the residual row does not say what it is:\n{report}"
    );
    assert!(
        !report
            .lines()
            .any(|l| l.contains("residual") && l.contains('%')),
        "a residual is still being given a percentage:\n{report}"
    );
}

/// **A scene of pure systems has no finer statement, and says so instead of "verifying" a
/// no-op.** Refining nothing and comparing a run against itself would print a zero shift and
/// look like the cleanest convergence in the file — the empty-panel failure shape, one level
/// up.
#[test]
fn a_scene_of_pure_systems_has_nothing_to_refine() {
    let s = scene(
        r#"{
  "title": "a heater warming a lump",
  "duration_s": 10.0,
  "frames": 5,
  "domains": [
    { "kind": "heater", "name": "element", "watts": 5.0, "reserve_j": 1000.0 },
    { "kind": "lump", "name": "plate", "volume_cm3": 0.5, "thickness_mm": 2.0,
      "initial_c": 20.0, "ambient_c": 20.0, "area_cm2": 400.0 }
  ]
}"#,
    );
    let why = s.refined().expect_err("nothing here is a discretisation");
    assert!(
        why.contains("nothing in this scene is a discretisation"),
        "{why}"
    );
}

/// **The window error of the once-per-frame feedback is first order, and the battery sees it.**
///
/// The `tracks` loop is closed between frames because it cannot be closed inside the step
/// loop, and its own documentation states the price: the winding's temperature lags by up to
/// one frame interval. A resistance evaluated at a stale temperature dissipates the wrong
/// watts for the whole window, so the accumulated error is first order in the window — the
/// error class the conservation audit structurally cannot see, since every joule the winding
/// spends does arrive. The battery's window sweep must read an order near 1 here; reading 2
/// would mean it cannot see the one error class it most exists to see. (The other candidate,
/// a subcycling consumer taking a deposited interval, stopped being first order when
/// `Exchange::take_share` fixed it — measured by `multirate_timing.rs`, and measured again
/// here in the negative: this scene is the one with a first-order knob left.)
///
/// The band is earned by measurement against the theory's limit, **and it holds only at this
/// pin**. First order is the *limit* as the window shrinks: at 20 s windows against a ~90 s
/// time constant the order read 0.88, at 10 s it read 1.37 — outside this band, one octave
/// from the asymptotic regime — and at the 5 s windows this test uses it reads 1.01, with the
/// instantaneous readings at 1.00 confirming the regime. So `[0.8, 1.2]` around a measured
/// 1.01: the width is against benign code changes moving the arithmetic path, not against
/// platforms — results are bit-identical across platforms by convention 3, so that costs no
/// width at all — and it is narrow enough that second order cannot pass.
#[test]
fn the_window_error_of_the_tracks_feedback_is_first_order() {
    let s = scene(
        r#"{
  "title": "a winding tracking its own heat",
  "duration_s": 60.0,
  "frames": 12,
  "domains": [
    { "kind": "winding", "name": "coil", "length_m": 10.0, "cross_section_mm2": 0.1,
      "amps": 3.0, "at_c": 20.0, "reserve_j": 40000.0, "tracks": "motor/coil" },
    { "kind": "network", "name": "motor",
      "nodes": [ { "name": "coil", "material": "copper", "volume_cm3": 2.0,
                   "thickness_mm": 2.0, "initial_c": 20.0,
                   "loses_to": { "ambient_c": 20.0, "area_cm2": 100.0 } } ],
      "links": [],
      "absorbing": "coil" }
  ]
}"#,
    );
    let b = verify(&s, true).expect("the battery runs");
    let p = order_of(&b, "window", "coil", "spent");
    assert!(
        (0.8..=1.2).contains(&p),
        "the tracks-feedback error should read first order, measured {p:.3}"
    );
}

/// **A room mode with commensurate dimensions converges at the scheme's second order.**
///
/// Width 4 m at 0.1 m spacing and height 2 m — exactly 20 cells — stay commensurate through
/// both refinements, so the quantised geometry does not move and what remains is the scheme:
/// second order in space, second order in time, first-order startup fixed. The measured order
/// of the peak reading should sit near 2. Under the incommensurate default scene it reads
/// near 1, because the *height itself* converges at first order — measured and then confirmed
/// against `f(1,1)` by hand — which is why this test pins the commensurate case: it is the one
/// with a clean known answer.
///
/// Measured: 2.24. Above 2 rather than below because the substep count is a ceiling, so the
/// effective time step shrinks slightly faster than the spacing; 2.4 is that ceiling argument
/// with room, and 1.6 is where a first-order term stops being subdominant. What the band can
/// and cannot refuse is a property of a three-point measurement worth stating exactly: with
/// error `A·h + B·h²`, `p ≥ 1.6` admits a first-order term contributing up to ~40% of the
/// coarse-grid difference — but that term *doubles* its share per refinement, so it is caught
/// one octave later, and a **dominant** first-order defect (both shipped acoustic defects
/// were) reads near 1 and fails here outright.
#[test]
fn a_commensurate_room_mode_converges_at_second_order() {
    let s = scene(
        r#"{
  "title": "a commensurate room",
  "duration_s": 0.02,
  "frames": 5,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 4.0, "height_m": 2.0,
      "cells_across": 41,
      "release": { "as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0 } }
  ]
}"#,
    );
    let b = verify(&s, true).expect("the battery runs");
    let p = order_of(&b, "resolution", "room", "peak");
    assert!(
        (1.6..=2.4).contains(&p),
        "a second-order scheme on unmoving geometry, measured {p:.3}"
    );
}

/// **A scene whose answer depends on how many pictures it was asked for is a finding.**
///
/// The battery has always measured this — the window sweep halves the coupling window and reports
/// what moved — and reporting it was not enough. `15-a-hot-spot-in-a-block` shipped at
/// `frames: 11` and is **14%** below its own grid's converged answer: run at 2, 11, 101 and 1601
/// frames its peak reads 5.177, 7.027, 8.065 and 8.186 K, which is textbook first order in the
/// step. Its window row said **0.212%** and nobody acted, because the shift was being divided by
/// the reading's own magnitude — and for a celsius temperature that is dominated by the 273.15
/// the scale carries. Against what the run actually did it is 0.966%.
///
/// So it carries the exit code now. The threshold is `WINDOW_SHIFT`, half a percent, chosen
/// against the corpus and written down beside itself.
#[test]
fn a_scene_whose_answer_moves_with_the_frame_count_is_a_finding() {
    // One hot cell in aluminium: diffusion, so the error is first order in the step and a coarse
    // window shows immediately. Eleven frames over 6 ms is the shipped scene's own setting.
    let coarse = scene(
        r#"{
  "title": "one hot cell, too few windows",
  "duration_s": 0.006,
  "frames": 11,
  "domains": [
    { "kind": "block", "name": "block", "cells": [9, 9, 9], "cell_mm": 1.0,
      "initial_c": 20.0, "material": "aluminium",
      "hot_spot": { "at": [4, 4, 4], "above_k": 60.0 } }
  ]
}"#,
    );
    let b = verify(&coarse, false).expect("the battery runs");
    let named: Vec<&String> = b
        .findings
        .iter()
        .filter(|f| f.contains("depends on `frames`"))
        .collect();
    assert_eq!(
        named.len(),
        1,
        "a scene 14% from converged reported no finding: {:?}",
        b.findings
    );
    assert!(
        ran(&b.window).worst > 0.005,
        "the sweep measured {:.4} and the finding fired anyway",
        ran(&b.window).worst
    );

    // And the same scene with a window stated is clean, which is what makes the finding a thing a
    // reader can act on rather than one they have to live with.
    let mut settled = coarse.clone();
    settled.window_s = Some(0.006 / 1600.0);
    let b = verify(&settled, false).expect("the battery runs");
    assert!(
        !b.findings.iter().any(|f| f.contains("depends on `frames`")),
        "stating a window did not settle it: {:?}",
        b.findings
    );
}

/// **A block whose cells all hold the same number is a finding.**
///
/// Every other finding in this battery is about arithmetic — determinism, a sweep, a drift. This
/// one is about whether the scene needed the arithmetic. A field solved on a grid and answered as
/// a lump is not wrong; it is a scene that will pass every check it has, converge perfectly, and
/// answer a question one ordinary differential equation answers.
///
/// It shipped twice. `30-two-phases-crossing-at-a-clearance` states two busbars of thirty-two
/// cells each whose temperature spread is **exactly zero** — uniform copper, uniform dissipation,
/// uniform cooling — and `29-a-designed-bracket-becomes-cells` rasterises 4 100 cells from an STL
/// to hold a 0.067 K range, a Biot number of 8.6e-4. Nothing in the report mentioned either.
///
/// # Both directions, in one scene each
///
/// A bar heated at one end and cooled at the other **has** a gradient and must not be named; the
/// same bar heated uniformly and cooled on every face has none and must be. The two differ by
/// where the watts go, so this is measuring the field and not the wording.
#[test]
fn a_field_that_is_flat_is_a_finding_and_a_gradient_is_not() {
    let bar = |dissipation: &str, cooling: &str| {
        scene(&format!(
            r#"{{
  "title": "a bar",
  "duration_s": 20.0,
  "frames": 6,
  "domains": [
    {{ "kind": "block", "name": "bar", "cells": [12, 2, 2], "cell_mm": 3.0,
      "initial_c": 20.0, "material": "copper",
      "dissipation": [{dissipation}],
      "cooling": [{cooling}] }}
  ]
}}"#
        ))
    };

    // Heat in at one end, out at the other: a real gradient, and no finding.
    let gradient = verify(
        &bar(
            r#"{ "watts": 40.0, "from": [0, 0, 0], "to": [1, 2, 2] }"#,
            r#"{ "face": "x-max", "ambient_c": 20.0, "convection_w_per_m2_k": 4000.0,
                 "area_cm2": 0.36 }"#,
        ),
        false,
    )
    .expect("it runs");
    let flat: Vec<&String> = gradient
        .findings
        .iter()
        .filter(|f| f.contains("flat to"))
        .collect();
    assert!(
        flat.is_empty(),
        "a bar with a gradient was called flat: {flat:?}"
    );

    // The bracket's own shape: a copper block starting hot and shedding to air, where the film is
    // so much weaker than the conduction that the whole thing falls together. The Biot number
    // `h * L_c / k` is 3.5e-4 here and 8.6e-4 on the bracket, so this is the shipped defect at
    // test scale rather than a contrivance. Uniform watts with cooling at both ends is *not* it:
    // the middle stays hotter, and the first version of this test measured 0 findings for that
    // reason.
    let cooling_only = scene(
        r#"{
  "title": "a copper block cooling to air",
  "duration_s": 200.0,
  "frames": 6,
  "domains": [
    { "kind": "block", "name": "bar", "cells": [12, 2, 2], "cell_mm": 3.0,
      "initial_c": 200.0, "material": "copper",
      "cooling": [
        { "face": "x-max", "ambient_c": 20.0, "convection_w_per_m2_k": 100.0, "area_cm2": 0.36 },
        { "face": "x-min", "ambient_c": 20.0, "convection_w_per_m2_k": 100.0, "area_cm2": 0.36 }
      ] }
  ]
}"#,
    );
    let lump = verify(&cooling_only, false).expect("it runs");
    let named = lump
        .findings
        .iter()
        .find(|f| f.contains("flat to"))
        .unwrap_or_else(|| {
            panic!(
                "a lump written as a grid was not named: {:?}",
                lump.findings
            )
        });
    println!("  {named}");
    assert!(named.starts_with("bar:"), "which domain: {named}");
    assert!(
        named.contains("lumped node"),
        "the finding should say what it means: {named}"
    );
    // And it reaches the report a person reads, not only the struct.
    assert!(lump.render().contains("flat to"), "{}", lump.render());
}

/// **A run that did nothing cannot report its own flatness.**
///
/// The denominator is the range the readings covered, so a scene that barely moved would divide a
/// rounding difference by a rounding difference and call the answer a structure — or, worse, call
/// a block that never changed "flat" and be technically right about nothing. The guard is the same
/// representation floor the sweeps use.
#[test]
fn a_run_that_did_nothing_is_not_reported_as_flat() {
    let s = scene(
        r#"{
  "title": "a bar at rest",
  "duration_s": 1e-6,
  "frames": 3,
  "domains": [
    { "kind": "block", "name": "still", "cells": [6, 2, 2], "cell_mm": 3.0,
      "initial_c": 20.0, "material": "copper" }
  ]
}"#,
    );
    let b = verify(&s, false).expect("it runs");
    let flat: Vec<&String> = b
        .findings
        .iter()
        .filter(|f| f.contains("flat to"))
        .collect();
    assert!(
        flat.is_empty(),
        "a block that never moved was reported as flat: {flat:?}"
    );
}

/// **How far a run stopped from its own steady state is a measurement, not a judgement.**
///
/// Every reading a design asks for is a steady-state one — a junction temperature, a margin, a rise
/// above ambient — and the only way to know a run had reached one was to look at the last two
/// frames and decide. `Solid3D::steady_state` solves the balance the march converges to, using the
/// same operator, so the distance between them is a number.
///
/// Both directions, in one scene each: a bar run to ten time constants has arrived, and the same
/// bar run to a tenth of one has not. What is asserted is the *ordering* and the scale, not a
/// threshold — the report carries the number and this carries that the number means something.
#[test]
fn the_battery_says_how_far_a_run_stopped_from_its_balance() {
    let bar = |seconds: f64| {
        scene(&format!(
            r#"{{
  "title": "a bar warming",
  "duration_s": {seconds},
  "frames": 6,
  "conservation_tolerance": 1e-6,
  "domains": [
    {{ "kind": "block", "name": "bar", "cells": [10, 2, 2], "cell_mm": 3.0,
      "initial_c": 20.0, "material": "copper",
      "dissipation": [{{ "watts": 8.0, "from": [0, 0, 0], "to": [10, 2, 2] }}],
      "cooling": [{{ "face": "x-max", "ambient_c": 20.0,
        "convection_w_per_m2_k": 2000.0, "area_cm2": 0.36 }}] }}
  ]
}}"#
        ))
    };

    // **Ten time constants, and the time constant is short on purpose.** The bar's capacity
    // is 3.73 J/K against 0.072 W/K of film, so tau is 52 s and 500 s is 9.6 of them, leaving
    // `e^-9.6` of the approach. A gentler film made tau 518 s and this test 4000 s of marching,
    // four times over inside the battery — which is the cost this whole capability exists to
    // avoid, paid by the test for it.
    let arrived = verify(&bar(500.0), false).expect("it runs");
    let (_, moved, travel) = arrived
        .base
        .arrival
        .first()
        .unwrap_or_else(|| panic!("nothing was measured: {:?}", arrived.base.arrival));
    println!("  500 s: {moved:.6} K still to come of {travel:.4} C covered");
    assert!(*travel > 1.0, "the run did nothing: {travel:.6}");
    assert!(
        moved / travel < 1e-3,
        "a run of ten time constants should have arrived: {moved:.6} K of {travel:.4}"
    );

    let cut_short = verify(&bar(5.0), false).expect("it runs");
    let (_, short_moved, short_travel) = cut_short
        .base
        .arrival
        .first()
        .expect("the short run is measured too");
    println!("  5 s:    {short_moved:.6} K still to come of {short_travel:.4} C covered");
    assert!(
        short_moved / short_travel > 0.05,
        "a run of a tenth of a time constant has not arrived: {short_moved:.6} K of \
         {short_travel:.4}"
    );
    // The short run's answer really is the smaller one, which is what "still to come" means.
    assert!(
        short_moved > moved,
        "the shorter run should have further to go: {short_moved:.6} against {moved:.6}"
    );

    // And it reaches the report a person reads.
    let report = cut_short.render();
    assert!(report.contains("still to come"), "{report}");
}

/// **A domain with no steady state to be short of is left out, not reported as arrived.**
///
/// The difference between "it arrived" and "the question does not apply". A block generating heat
/// with every face insulated warms without limit; reporting it as zero kelvin from its balance
/// would be a claim that it had one.
#[test]
fn a_domain_with_no_balance_is_left_out_of_the_arrival() {
    let s = scene(
        r#"{
  "title": "a sealed bar",
  "duration_s": 10.0,
  "frames": 3,
  "conservation_tolerance": 1e-6,
  "domains": [
    { "kind": "block", "name": "sealed", "cells": [4, 2, 2], "cell_mm": 3.0,
      "initial_c": 20.0, "material": "copper",
      "dissipation": [{ "watts": 2.0, "from": [0, 0, 0], "to": [4, 2, 2] }] }
  ]
}"#,
    );
    let b = verify(&s, false).expect("it runs; it just never settles");
    println!("  arrival: {:?}", b.base.arrival);
    assert!(
        b.base.arrival.is_empty(),
        "a block that warms without limit has no balance to be short of: {:?}",
        b.base.arrival
    );
    assert!(
        !b.render().contains("still to come"),
        "the report should not claim an arrival it did not measure:\n{}",
        b.render()
    );
}

/// **A body asked past the strain its material comes back from is a finding.**
///
/// One is a **physical boundary and not a chosen threshold**, which is what makes this different
/// from every other number in this battery. `pantometry-elastic` is linear and has no plasticity,
/// and it documents its own limit in as many words: past yield it "returns a displacement that is
/// arithmetically correct and physically meaningless, and nothing in the answer says which".
///
/// The crate knew. Nothing checked a scene against it, and
/// `25-what-140-kelvin-does-to-the-solder` ships at **5.20×** — SAC305 assembled at its 217 °C
/// reflow and cooled to 40 takes a free strain of 0.3806% against a yield strain of 0.0732%.
///
/// Both directions: the same body with a reference temperature it is already at goes nowhere, and
/// is not named.
#[test]
fn a_body_strained_past_its_yield_is_a_finding() {
    // A copper bar bonded to a block that cools a long way from where it was assembled.
    let assembly = |reference_c: f64| {
        scene(&format!(
            r#"{{
  "title": "a bonded bar",
  "duration_s": 1.0,
  "frames": 3,
  "conservation_tolerance": 1e-6,
  "domains": [
    {{ "kind": "block", "name": "hot", "cells": [4, 2, 2], "cell_mm": 2.0,
      "initial_c": 20.0, "material": "copper" }},
    {{ "kind": "structure", "name": "body", "cells": [4, 2, 2], "cell_mm": 2.0,
      "material": "copper", "follows": "hot", "reference_c": {reference_c},
      "held": [{{ "face": "x-min", "as": "clamp" }}] }}
  ]
}}"#
        ))
    };

    // Copper yields at 0.060% strain and expands at 1.65e-5 per kelvin, so 200 K from its
    // reference is 0.33% — five and a half times over.
    let strained = verify(&assembly(220.0), false).expect("it runs");
    let (name, ratio, at) = strained
        .base
        .past_yield
        .first()
        .unwrap_or_else(|| panic!("nothing was measured: {:?}", strained.base.past_yield));
    println!("  reference 220 C: {name} reaches {ratio:.3}x at {at:?}");
    assert!(*ratio > 1.0, "200 K of copper is past yield: {ratio:.3}x");
    let named = strained
        .findings
        .iter()
        .find(|f| f.contains("comes back from"))
        .unwrap_or_else(|| panic!("it was measured and not raised: {:?}", strained.findings));
    assert!(named.starts_with("body:"), "which domain: {named}");
    assert!(
        strained.render().contains("outside the linear model"),
        "{}",
        strained.render()
    );

    // The same body assembled where it sits takes no free strain at all.
    let easy = verify(&assembly(20.0), false).expect("it runs");
    let (_, easy_ratio, _) = easy.base.past_yield.first().expect("it is measured too");
    println!("  reference 20 C:  {easy_ratio:.6}x");
    assert!(
        *easy_ratio < 1.0,
        "a body at its own reference is inside the model: {easy_ratio:.6}x"
    );
    assert!(
        !easy.findings.iter().any(|f| f.contains("comes back from")),
        "a body inside the model was named: {:?}",
        easy.findings
    );
}

/// **A structure refines with the block it follows**, which took two goes to get right.
///
/// A structure's elements have to be that block's cells — `World::build` refuses the pair
/// otherwise — so refining one without the other is not a finer statement of the same problem. It
/// has been wrong in both directions:
///
/// 1. `DomainSpec::refined` passed a structure through unchanged, with a comment saying it was
///    "reported as unswept". It was not. The block doubled, the structure did not, and the sweep
///    reported the whole refined scene as **refused** — so
///    `25-what-140-kelvin-does-to-the-solder` had been *failing* its own resolution sweep since it
///    shipped, which read as a finding about the scene rather than about the sweep.
/// 2. Then it skipped, with a reason, which was honest and still left the scene unmeasured.
///
/// Refining the pair is a statement about two domains at once, which a function seeing one cannot
/// make. [`Scene::refined`](pantometry_world::Scene::refined) sees all of them, so it doubles the
/// block first, then every structure that follows one, and refuses if the block it names did not
/// double.
///
/// What the skip cost: scene 25's strain energy moves **4.09%** between 512 elements and 4096 —
/// 0.0274238 J against 0.0263010 — and the strain the scene is about moves 1.48% in `z`.
#[test]
fn a_structure_refines_with_the_block_it_follows() {
    let s = scene(
        r#"{
  "title": "a bonded bar",
  "duration_s": 1.0,
  "frames": 3,
  "conservation_tolerance": 1e-6,
  "domains": [
    { "kind": "block", "name": "hot", "cells": [4, 2, 2], "cell_mm": 2.0,
      "initial_c": 20.0, "material": "copper" },
    { "kind": "structure", "name": "body", "cells": [4, 2, 2], "cell_mm": 2.0,
      "material": "copper", "follows": "hot", "reference_c": 40.0,
      "held": [{ "face": "x-min", "as": "clamp" }] }
  ]
}"#,
    );

    // The pair doubles together, and stays the same box.
    let fine = s.refined().expect("a structure refines with its block");
    let mut seen = 0;
    for spec in &fine.domains {
        match spec {
            DomainSpec::Block { cells, cell_mm, .. }
            | DomainSpec::Structure { cells, cell_mm, .. } => {
                assert_eq!(*cells, [8, 4, 4]);
                assert_eq!(*cell_mm, 1.0);
                seen += 1;
            }
            other => panic!("{other:?}"),
        }
    }
    assert_eq!(seen, 2, "both halves of the pair have to refine");
    // Which is the thing that matters: the refined scene has to *build*, and it did not before.
    pantometry_world::World::build(fine).expect("the refined pair is a scene");

    // And the sweep runs rather than skipping or being refused.
    let b = verify(&s, false).expect("it runs");
    let sweep = ran(&b.resolution);
    assert!(
        sweep.shifts.iter().any(|x| x.domain == "body"),
        "the structure was not measured: {:?}",
        sweep.shifts
    );
    assert!(
        !b.findings
            .iter()
            .any(|f| f.contains("sweep's run was refused")),
        "the refined pair was refused: {:?}",
        b.findings
    );
}

/// **A structure whose partner did not refine is refused, naming it.**
///
/// The pair has to double together or `World::build` rejects it with a message about element
/// counts, which says nothing about why. Two ways the partner can fail to double, and only one of
/// them is caught by the block's own refusal:
///
/// - the block **refuses**, for a one-cell feature, and the whole scene refuses with its reason;
/// - the named domain has **no grid at all** — a lump, a heater, a network — so it passes through
///   unchanged and says nothing. That one is this check's, and a sabotage that removed it passed a
///   test using the first case, because the block refused before the structure had its turn.
#[test]
fn a_structure_whose_partner_did_not_refine_is_refused() {
    // A block that refuses for its own reason: the scene refuses with that reason.
    let hot_spot = scene(
        r#"{
  "title": "a body on a spotted block",
  "duration_s": 1.0, "frames": 3, "conservation_tolerance": 1e-6,
  "domains": [
    { "kind": "block", "name": "hot", "cells": [4, 2, 2], "cell_mm": 2.0,
      "initial_c": 20.0, "material": "copper",
      "hot_spot": { "at": [2, 1, 1], "above_k": 10.0 } },
    { "kind": "structure", "name": "body", "cells": [4, 2, 2], "cell_mm": 2.0,
      "material": "copper", "follows": "hot", "reference_c": 40.0,
      "held": [{ "face": "x-min", "as": "clamp" }] }
  ]
}"#,
    );
    let why = hot_spot.refined().expect_err("a hot spot must refuse");
    assert!(why.contains("hot spot"), "{why}");

    // A partner with no grid: it passes through in silence, and the structure has to say so.
    let lump = scene(
        r#"{
  "title": "a body on a lump",
  "duration_s": 1.0, "frames": 3, "conservation_tolerance": 1e-6,
  "domains": [
    { "kind": "lump", "name": "hot", "volume_cm3": 8.0, "thickness_mm": 2.0,
      "initial_c": 20.0, "ambient_c": 20.0, "area_cm2": 24.0 },
    { "kind": "block", "name": "other", "cells": [4, 2, 2], "cell_mm": 2.0,
      "initial_c": 20.0, "material": "copper" },
    { "kind": "structure", "name": "body", "cells": [4, 2, 2], "cell_mm": 2.0,
      "material": "copper", "follows": "hot", "reference_c": 40.0,
      "held": [{ "face": "x-min", "as": "clamp" }] }
  ]
}"#,
    );
    let why = lump
        .refined()
        .expect_err("a structure following something with no grid must be refused");
    assert!(why.contains("did not refine"), "{why}");
    assert!(
        why.contains("\"hot\""),
        "the refusal should name the partner: {why}"
    );
}
