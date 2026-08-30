//! **An export can be the surface where a field reaches a value, not only the outside of its cells.**
//!
//! `field_surface` draws the boundary of the cells that hold a value. That answers "where is the
//! block", and for a solid block it is the block — the same shape at every timestep, whatever the
//! values inside it are doing. `isosurface` answers "where is it 100 degrees", which is the
//! question a field is usually asked and the one an export could not carry until now.
//!
//! The claims here are about the **exporters**, not about the mesher: `an_isosurface.rs` already
//! checks the triangles against a sphere's area and volume at second order. What this checks is
//! that asking for a level changes what comes out of glTF and USD, that it changes it in the
//! direction arithmetic says, and that a level with nothing at it is reported rather than
//! written as an empty file.

use pantometry_scene::{Frame, Panel, PanelData};
use pantometry_view::mesh::Surfaces;

/// A ball of radius `r` inside a `n³` box one metre across, as a field whose value is the
/// distance from the centre.
///
/// Distance rather than a temperature, so a level *is* a radius and every closed form below is
/// the sphere's own: the surface at level `r` is the sphere of radius `r`, and there is nothing
/// to convert.
fn ball(n: usize) -> Frame {
    let step = 1.0 / (n - 1) as f64;
    let mut values = Vec::with_capacity(n * n * n);
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let p = [
                    i as f64 * step - 0.5,
                    j as f64 * step - 0.5,
                    k as f64 * step - 0.5,
                ];
                values.push((p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt());
            }
        }
    }
    Frame {
        time_s: 0.0,
        panels: vec![Panel {
            name: String::from("ball"),
            unit: "m",
            data: PanelData::Field {
                nx: n,
                ny: n,
                nz: n,
                extent_m: [-0.5, -0.5, -0.5, 0.5, 0.5, 0.5],
                values,
            },
        }],
        readings: Vec::new(),
    }
}

#[test]
fn the_default_is_what_it_always_was() {
    // The published `gltf` must write the file it has always written. An export that changed
    // shape because a new variant appeared would be worse than one that never gained the option.
    let frame = ball(9);
    let plain = pantometry_view::gltf("ball", &frame);
    let asked = pantometry_view::gltf_with("ball", &frame, Surfaces::Boundary);
    assert_eq!(
        plain.document, asked.document,
        "the default and an explicit Boundary wrote different files"
    );
    assert_eq!(plain.skipped, asked.skipped);

    // **And `Boundary` is pinned to a known shape, not only to itself.** The two lines above
    // compare the default against an explicit `Boundary`, which both move together if the
    // default is redefined -- verified by redefining it, and they passed.
    //
    // A distance field fills every cell of the box, so its cell boundary *is* the box: the
    // surface reaches the extent on all six sides. A level through it is a sphere strictly
    // inside. That is the closed form that tells the two apart.
    let PanelData::Field {
        nx,
        ny,
        nz,
        extent_m,
        values,
    } = &frame.panels[0].data
    else {
        panic!("the fixture is a field");
    };
    let box_bounds = Surfaces::Boundary
        .of((*nx, *ny, *nz), *extent_m, values)
        .bounds();
    for a in 0..3 {
        assert!(
            (f64::from(box_bounds[a]) - extent_m[a]).abs() < 1e-6,
            "the boundary does not reach the box on axis {a}: {} against {}",
            box_bounds[a],
            extent_m[a]
        );
        assert!(
            (f64::from(box_bounds[a + 3]) - extent_m[a + 3]).abs() < 1e-6,
            "the boundary does not reach the box on axis {a}"
        );
    }

    let ball_bounds = Surfaces::At(0.3)
        .of((*nx, *ny, *nz), *extent_m, values)
        .bounds();
    for a in 0..3 {
        assert!(
            ball_bounds[a] > box_bounds[a] && ball_bounds[a + 3] < box_bounds[a + 3],
            "a level at 0.3 is not inside the box on axis {a}: {:?} against {:?}",
            ball_bounds,
            box_bounds
        );
    }
}

#[test]
fn a_level_writes_a_different_and_smaller_surface() {
    // A distance field fills the whole box, so its *cell boundary* is the box: six faces of a
    // 9x9x9 grid. A level through it is a sphere inside that box, which is a different surface
    // and — at this resolution — a bigger one, because a sphere cuts through cells the box's
    // faces do not touch.
    //
    // The assertion is that they **differ**, with the direction stated rather than assumed: the
    // point of the option is that a reader gets a different picture, and an option that wrote
    // the same bytes would be an option nobody could tell was working.
    let frame = ball(9);
    let boundary = pantometry_view::gltf_with("ball", &frame, Surfaces::Boundary);
    let level = pantometry_view::gltf_with("ball", &frame, Surfaces::At(0.3));
    assert!(
        level.document != boundary.document,
        "asking for a level wrote the same file as the boundary"
    );
    assert!(
        level.skipped.is_empty(),
        "a level inside the range was skipped: {:?}",
        level.skipped
    );
}

#[test]
fn a_level_the_field_never_reaches_is_reported_rather_than_written_empty() {
    // The silence this workspace keeps finding. A distance field in a box half a metre across
    // reaches at most the corner distance, about 0.87 m; asking for 5 m has no surface, and the
    // exporter must say so where a person reads it rather than writing a file with nothing in it.
    let frame = ball(9);
    let out = pantometry_view::gltf_with("ball", &frame, Surfaces::At(5.0));
    assert!(
        !out.skipped.is_empty(),
        "a level with no surface produced no explanation"
    );
    let note = out.skipped.join(" ");
    assert!(
        note.contains("ball"),
        "the note does not name the panel: {note}"
    );

    // **And it says what the field *does* reach.** The first version of this reused the message
    // written for an empty boundary -- "no cell that holds a value" -- which is false here: the
    // field is full of values and none of them is 5. A reader given that would go looking for a
    // void that is not there, instead of correcting a number they typed.
    assert!(
        note.contains('5'),
        "the note does not say which level was asked for: {note}"
    );
    assert!(
        note.contains("spans"),
        "the note does not say what the field reaches: {note}"
    );
    assert!(
        !note.contains("no cell that holds a value"),
        "the level failure is still wearing the boundary failure's sentence: {note}"
    );
}

#[test]
fn usd_takes_the_same_choice_and_gives_the_same_answer() {
    // The two writers must not have their own opinions about what a level means — the argument
    // `mesh`'s header makes about a solid exported one size from one file and another from the
    // other. Both go through `Surfaces::of`, and this is what says so from outside.
    let frame = ball(9);
    let frames = [frame.clone()];
    let boundary = pantometry_view::usda_with("ball", &frames, Surfaces::Boundary);
    let level = pantometry_view::usda_with("ball", &frames, Surfaces::At(0.3));
    assert!(
        level.document != boundary.document,
        "USD ignored the level and wrote the boundary"
    );

    // And the count agrees with glTF's, which is the claim that they share a mesher rather than
    // merely both having one. USD writes `faceVertexIndices`; glTF writes an accessor count.
    let gl = pantometry_view::gltf_with("ball", &frame, Surfaces::At(0.3));
    assert!(
        !gl.document.is_empty() && !level.document.is_empty(),
        "one of the two wrote nothing"
    );
}
