//! **A run states where each panel's frame is, and a reader that ignores it draws a lie.**
//!
//! Run format 2 stopped baking each domain's pose into its numbers and started stating it, so a
//! placed domain's samples and bodies are in its **own** coordinates and `place` is where they go.
//! This crate parsed that key from the day it existed and nothing read it: `Placed::is_here` was
//! defined and never called, `bounds()` said "in world coordinates" and returned the panel's own
//! box, and the editor's viewport drew a placed part's outline where the scene put it — from
//! `editor-core`, which does apply the pose — and its colours at the origin.
//!
//! Two frames in one picture, and nothing on screen saying which a thing was in. That is the same
//! defect the exporters had, in the third reader, and it was found the same way: by placing two
//! identical things and looking at where they landed.
//!
//! # What is checked here, and what cannot be
//!
//! The composition — turn, then move — is a fact about the *format*, and the format's other end is
//! `pantometry_view`, which this crate deliberately does not link. So it is written out by hand
//! here from trigonometry, and `the_wire_format_is_enough.rs` explains why that separation is worth
//! paying for.

use viewer_core::{Panel, Placed};

/// A quarter turn about z, in `[x, y, z, w]`.
fn quarter_about_z() -> [f64; 4] {
    let h = std::f64::consts::FRAC_PI_4;
    [0.0, 0.0, h.sin(), h.cos()]
}

/// An eighth of a turn about z, which is the one that makes a bounding box grow.
///
/// A **quarter** turn does not: it maps axes onto axes, so an axis-aligned box comes back
/// axis-aligned and the same size. The first version of the box test below used one and asserted
/// √2 anyway, and what it measured was 1.0000000000000002 — the test was wrong and the code was
/// right, which is the direction worth saying out loud.
fn eighth_about_z() -> [f64; 4] {
    let h = std::f64::consts::PI / 8.0;
    [0.0, 0.0, h.sin(), h.cos()]
}

/// A unit-cube field panel, placed as given.
fn cube(place: Placed) -> Panel {
    Panel::Field {
        name: "cube".into(),
        unit: "K".into(),
        place,
        nx: 2,
        ny: 2,
        nz: 2,
        extent_m: Some([0.0, 0.0, 0.0, 1.0, 1.0, 1.0]),
        values: vec![300.0; 8],
    }
}

#[test]
fn the_identity_moves_nothing_at_all() {
    // Every run this workspace wrote before format 2 is this case, and so is all but one shipped
    // scene. It has to be **exact**: `world_bounds` returning a rounded copy of `bounds` would put
    // a floating-point wobble into every camera fit in the editor for no reason.
    let here = Placed::default();
    assert!(here.is_here());
    for p in [[0.0, 0.0, 0.0], [1.0, -2.0, 3.5], [1e-9, 1e9, -7.0]] {
        assert_eq!(here.apply(p), p, "the identity moved {p:?}");
    }
    let panel = cube(here);
    assert_eq!(panel.world_bounds(), panel.bounds());
    assert_eq!(panel.place(), Placed::default());
}

#[test]
fn a_quarter_turn_about_z_takes_x_to_y() {
    let place = Placed {
        at_m: [0.0; 3],
        turn: quarter_about_z(),
    };
    let got = place.apply([1.0, 0.0, 0.0]);
    for (a, want) in [0.0, 1.0, 0.0].into_iter().enumerate() {
        assert!(
            (got[a] - want).abs() < 1e-12,
            "x turned to {got:?}, not [0, 1, 0]"
        );
    }
    // And z is the axis, so it does not move.
    assert!((place.apply([0.0, 0.0, 5.0])[2] - 5.0).abs() < 1e-12);
}

#[test]
fn it_turns_and_then_moves_rather_than_the_other_way() {
    // **The assertion that catches the composition being backwards.** A glTF node and a USD
    // `xformOpOrder = [translate, orient]` both turn first, so a viewer that moved first would
    // put a placed part somewhere no exporter agrees with — and for a placement at the origin, or
    // a turn of nothing, the two orders are identical and nothing would say so.
    let place = Placed {
        at_m: [10.0, 0.0, 0.0],
        turn: quarter_about_z(),
    };
    // Turn first: (1,0,0) → (0,1,0), then move → (10, 1, 0).
    // Move first would give (11,0,0) → (0, 11, 0), which is nowhere near.
    let got = place.apply([1.0, 0.0, 0.0]);
    for (a, want) in [10.0, 1.0, 0.0].into_iter().enumerate() {
        assert!(
            (got[a] - want).abs() < 1e-12,
            "got {got:?}; moving before turning would give [0, 11, 0]"
        );
    }
    // The panel's own origin lands on `at_m` under either order, which is why the point above is
    // off the axis.
    assert_eq!(place.apply([0.0, 0.0, 0.0]), place.at_m);
}

#[test]
fn corners_keep_the_shape_and_only_the_bounding_box_grows() {
    // **The information a bounding box loses, and why there are two methods.** An eighth of a
    // turn of a unit cube is still a unit cube; the axis-aligned box *around* it is √2 across on
    // two of three axes. `place` exists because the writer used to store the second and call it
    // the first, so a reader that only had `world_bounds` would be back where it started.
    let place = Placed {
        at_m: [0.0; 3],
        turn: eighth_about_z(),
    };
    let c = place.corners_of([0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);

    // Corner 0 to corner 1 is the box's x edge, whatever it was turned into.
    let edge = |a: usize, b: usize| {
        ((c[a][0] - c[b][0]).powi(2) + (c[a][1] - c[b][1]).powi(2) + (c[a][2] - c[b][2]).powi(2))
            .sqrt()
    };
    for (a, b) in [(0, 1), (0, 2), (0, 4)] {
        assert!(
            (edge(a, b) - 1.0).abs() < 1e-12,
            "edge {a}-{b} is {} long, not 1 — the turn is not rigid",
            edge(a, b)
        );
    }

    let b = cube(place).world_bounds();
    let across = std::f64::consts::SQRT_2;
    assert!(
        (b[3] - b[0] - across).abs() < 1e-12 && (b[4] - b[1] - across).abs() < 1e-12,
        "the bounding box is {} x {} across, not √2",
        b[3] - b[0],
        b[4] - b[1]
    );
    assert!(
        (b[5] - b[2] - 1.0).abs() < 1e-12,
        "z is the axis and must not grow"
    );
}

#[test]
fn two_identical_panels_placed_apart_do_not_share_a_box() {
    // **The defect, in one assertion.** Two of the same solid, half a metre apart. Their own
    // boxes are identical because they are the same solid — that is what made this invisible —
    // and only the placement can tell them apart.
    let near = cube(Placed::default());
    let far = cube(Placed {
        at_m: [0.5, 0.0, 0.0],
        ..Placed::default()
    });

    assert_eq!(
        near.bounds(),
        far.bounds(),
        "the two are the same solid in their own frames"
    );
    assert_ne!(
        near.world_bounds(),
        far.world_bounds(),
        "and a reader using only `bounds` draws them in one place"
    );
    assert_eq!(far.world_bounds()[0], 0.5);
    assert_eq!(far.world_bounds()[3], 1.5);
}

#[test]
fn a_run_written_before_the_key_reads_as_the_identity() {
    // The three recorded runs predate format 2 and state no placement anywhere. Their panels were
    // already in world coordinates, because that is what a run without one meant — so `place`
    // reading as the identity is not a fallback, it is the correct translation of an older file.
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/runs");
    let mut seen = 0;
    for entry in std::fs::read_dir(&dir).expect("the recorded runs are there") {
        let path = entry.expect("readable").path();
        if path.extension().and_then(|e| e.to_str()) != Some("json") {
            continue;
        }
        let text = std::fs::read_to_string(&path).expect("readable");
        assert!(
            !text.contains("\"place\""),
            "{} states a placement; this test is about files that do not",
            path.display()
        );
        let run = viewer_core::Run::from_json(&text).expect("it reads");
        for frame in &run.frames {
            for panel in &frame.panels {
                assert!(panel.place().is_here(), "{}", panel.name());
                assert_eq!(panel.world_bounds(), panel.bounds());
                seen += 1;
            }
        }
    }
    assert!(seen > 0, "no panels were checked");
}

/// A one-run path panel, placed as given, with `vertices` as given.
fn ray(place: Placed, vertices: Vec<f64>) -> Panel {
    Panel::Paths {
        name: "rays".into(),
        unit: "nm".into(),
        place,
        bounds: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        starts: vec![0.0],
        vertices,
        values: vec![486.1],
    }
}

#[test]
fn a_placed_path_projects_where_an_offset_one_does() {
    // **The standalone viewer's half.** `segments` projected the vertices as they came, so a run
    // that stated a placement drew its rays at the origin. That shell draws one panel at a time,
    // so nothing on screen looked wrong -- there was nothing beside it to be wrong against, which
    // is why this is checked by arithmetic rather than by opening it.
    //
    // The equivalence is the assertion: a path placed at `d` must project exactly where a path
    // whose vertices already carry `d` does. That is what "the placement is applied" means, and it
    // cannot be satisfied by applying some other transform.
    let d = [3.0, -2.0, 1.0];
    let local = vec![0.0, 0.0, 0.0, 1.0, 0.5, 0.0];
    let offset: Vec<f64> = local
        .chunks(3)
        .flat_map(|v| [v[0] + d[0], v[1] + d[1], v[2] + d[2]])
        .collect();

    let camera = viewer_core::Camera::default();
    let framing = viewer_core::Framing::of([-5.0, -5.0, -5.0, 5.0, 5.0, 5.0]);
    let span = (0.0, 1000.0);

    let placed = viewer_core::segments(
        &ray(
            Placed {
                at_m: d,
                ..Placed::default()
            },
            local.clone(),
        ),
        &camera,
        &framing,
        16.0 / 9.0,
        span,
    );
    let baked = viewer_core::segments(
        &ray(Placed::default(), offset),
        &camera,
        &framing,
        16.0 / 9.0,
        span,
    );
    let unplaced = viewer_core::segments(
        &ray(Placed::default(), local),
        &camera,
        &framing,
        16.0 / 9.0,
        span,
    );

    assert_eq!(placed.len(), 1);
    assert_eq!(placed, baked, "a placed path does not project where it is");
    assert_ne!(
        placed, unplaced,
        "the placement changed nothing, which is the defect"
    );
}
