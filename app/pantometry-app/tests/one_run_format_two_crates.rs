//! **The number the writer stamps and the number the reader accepts are one number.**
//!
//! They cannot be one constant. `pantometry_view::data` writes the run and lives in the library
//! workspace; `viewer_core` reads it and deliberately does **not** link `pantometry`, which is
//! the whole argument of `app/viewer-core/README.md` — a viewer that could reach the library
//! could accidentally depend on it, and then the wire format would stop being the boundary.
//!
//! So there are two constants, and the only thing that makes them one version is a test that can
//! see both. This binary can: it links the library and the viewer, which is what a binary that
//! runs *and* draws has to do.
//!
//! Without this, the failure is the quietest kind available. Bump the writer and not the reader
//! and every run refuses itself. Bump the reader and not the writer and a build claims to
//! understand a format nothing produces — which nobody notices until the version after.

#[test]
fn the_writer_and_the_reader_agree_on_the_version() {
    assert_eq!(
        pantometry::view::data::FORMAT,
        viewer_core::FORMAT,
        "the run writer stamps {} and the reader accepts up to {}; \
         one of the two was bumped without the other",
        pantometry::view::data::FORMAT,
        viewer_core::FORMAT
    );
}

#[test]
fn what_the_writer_produces_is_what_the_reader_reads() {
    // Not the constants — the **bytes**. A writer that agreed on the number and spelled the key
    // differently would pass the test above and produce files nothing could version-check, which
    // is the state this format was in until now.
    use pantometry::scene::{Frame, Panel, PanelData, Placed};

    let frames = vec![Frame {
        time_s: 0.0,
        panels: vec![Panel {
            name: String::from("bar"),
            unit: "K",
            place: Placed::HERE,
            data: PanelData::Field {
                nx: 2,
                ny: 1,
                nz: 1,
                values: vec![300.0, 301.0],
                extent_m: [0.0, 0.0, 0.0, 0.01, 0.0, 0.0],
            },
        }],
        readings: Vec::new(),
    }];

    let json = pantometry::view::data::to_json("a run", &frames);
    assert!(
        json.contains(&format!("\"format\": {}", pantometry::view::data::FORMAT)),
        "the writer did not stamp the run: {}",
        &json[..json.len().min(120)]
    );

    let run = viewer_core::Run::from_json(&json).expect("the reader takes what the writer wrote");
    assert_eq!(run.title, "a run");
    assert_eq!(run.frames.len(), 1);
    assert_eq!(run.frames[0].panels.len(), 1);
}
