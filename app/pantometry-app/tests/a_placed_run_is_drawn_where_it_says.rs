//! **The editor's viewport drew a placed run at the origin.**
//!
//! Run format 2 stopped baking each domain's pose into its numbers and started stating it, so a
//! placed domain's samples, bodies and paths are in its **own** coordinates and `Panel::place` is
//! where they go. `viewer-core` parsed that key and read nothing, and the six places in `batches`
//! that turn a panel into vertices used the numbers as they came.
//!
//! What made it worse than a plain overlay: the wireframe boxes come from `editor-core`, which
//! applies the pose itself, from the scene. So a scene that stated one drew its outline where it
//! belongs and its colours at the origin — two frames in one picture, and nothing on screen saying
//! which anything was in.
//!
//! # Why a subcommand and a number
//!
//! `batches` is a method on the editor's `App`, which is in a `bin` crate, and the alternative is a
//! `lib.rs` invented to let a test in. `--layout-at` set the precedent and `App::new` needs no
//! window. `--drawn-extent` loads a run, builds the batches, and inverts the framing to report what
//! it drew in metres.
//!
//! A number rather than a picture because **five of the six sites would look perfectly fine on any
//! scene that states no pose**, which is twenty-nine of the thirty shipped. A missed site collapses
//! that one panel onto the origin, which is a coordinate and not an appearance.
//!
//! # What this does not reach
//!
//! The **flat painter** — the unshaded fallback that draws with `egui::Painter` — has three sites
//! of its own and no headless hook. They take the same one-line change and are not covered here;
//! said rather than left to be assumed, because an uncovered site is exactly how this defect
//! reached three readers in the first place.

/// The binary, asked what it would draw for a run.
fn drawn(run: &str) -> String {
    let mut p = std::env::current_exe().expect("the test binary knows where it is");
    p.pop();
    if p.ends_with("deps") {
        p.pop();
    }
    let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/runs")
        .join(run);
    let out =
        std::process::Command::new(p.join(format!("pantometry{}", std::env::consts::EXE_SUFFIX)))
            .args(["--drawn-extent", &fixture.to_string_lossy()])
            .output()
            .expect("the binary runs");
    assert!(
        out.status.success(),
        "--drawn-extent failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// One `name=[a b c d e f]` box out of that line.
fn box_of(line: &str, name: &str) -> [f64; 6] {
    let at = line
        .split_once(&format!("{name}=["))
        .unwrap_or_else(|| panic!("no {name} in {line:?}"))
        .1;
    let inner = at.split_once(']').expect("a closed box").0;
    let v: Vec<f64> = inner
        .split_whitespace()
        .map(|n| n.parse().expect("a number"))
        .collect();
    assert_eq!(v.len(), 6, "{name} is {inner:?}");
    [v[0], v[1], v[2], v[3], v[4], v[5]]
}

/// **The batch stores `f32`.** Seven significant figures on a coordinate of 30 is about 2e-6, so
/// this is fifty times the precision the number can carry — and four orders below the ten metres
/// each placement moves things, which is the difference these assertions are about.
const F32_AT_THIRTY: f64 = 1e-4;

#[test]
fn each_shape_is_drawn_where_its_placement_says() {
    // The fixture holds one of each shape, placed on a different axis so that one number tells
    // which site failed: a field at x = 10, bodies at y = 20, a path at z = 30, each a unit box in
    // its own frame. Unplaced, every one of them sits inside `0..1`.
    let line = drawn("three-shapes-placed.json");
    let solid = box_of(&line, "solid");
    let lines = box_of(&line, "lines");

    assert!(
        (solid[3] - 11.0).abs() < F32_AT_THIRTY,
        "the field's far x is {}, not 11 — unplaced it is 1: {line}",
        solid[3]
    );
    // 21.25 and not 21: a body is drawn as a sphere and `mesh::body_radius` sizes it from the
    // run's own bounds, which for four bodies a metre apart is a quarter of a metre. The same
    // radius the exports use, so the viewport and Blender draw one body the same size.
    assert!(
        (solid[4] - 21.25).abs() < F32_AT_THIRTY,
        "the bodies' far y is {}, not 21.25 — unplaced it is 1.25: {line}",
        solid[4]
    );
    assert!(
        (lines[5] - 30.0).abs() < F32_AT_THIRTY,
        "the path's far z is {}, not 30 — unplaced the line batch reaches only the built-in \
         scene's 3.1: {line}",
        lines[5]
    );
}

#[test]
fn the_camera_frames_the_world_and_not_a_union_of_local_boxes() {
    // `Panel::bounds` is the box in whichever frame the panel is in. A union of boxes from
    // different frames is a box around nothing, and that is what the camera fitted to: with the
    // placements dropped this reads `[0 0 0 4.4 3.1 1]` — the built-in scene's room, with all
    // three panels folded into a corner of it.
    let line = drawn("three-shapes-placed.json");
    let world = box_of(&line, "world");
    for (a, want) in [11.0, 21.0, 30.0].into_iter().enumerate() {
        assert!(
            (world[a + 3] - want).abs() < F32_AT_THIRTY,
            "the world's far {a} is {}, not {want}: {line}",
            world[a + 3]
        );
    }
}

#[test]
fn a_run_that_states_no_placement_is_drawn_exactly_as_it_was() {
    // Twenty-nine of the thirty shipped scenes and every run written before format 2. The
    // placement being the identity has to cost nothing: a panel whose numbers went through an
    // identity quaternion and came back rounded would put a wobble into every camera fit for no
    // reason, and would do it everywhere rather than in the one case anybody would look at.
    let line = drawn("three-shapes-here.json");
    let solid = box_of(&line, "solid");
    let lines = box_of(&line, "lines");
    for a in 0..3 {
        assert!(
            solid[a + 3] <= 1.25 + F32_AT_THIRTY,
            "axis {a} reaches {} with nothing placed: {line}",
            solid[a + 3]
        );
    }
    // The line batch still carries the built-in scene's own box, which is not the run's.
    assert!(lines[5] <= 3.1 + F32_AT_THIRTY, "the paths moved: {line}");
}
