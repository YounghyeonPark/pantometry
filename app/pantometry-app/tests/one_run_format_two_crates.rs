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
                lattice: pantometry::scene::Lattice::Nodal,
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

/// Every key the writer emits, per place it can appear: the file, a frame, a reading, and a panel of
/// each kind. A panel's own `kind` is its row's name.
fn shape_of(json: &str) -> std::collections::BTreeMap<String, Vec<String>> {
    let v: serde_json::Value = serde_json::from_str(json).expect("the writer writes JSON");
    let keys = |o: &serde_json::Value| -> Vec<String> {
        let mut k: Vec<String> = o.as_object().expect("an object").keys().cloned().collect();
        k.sort();
        k
    };
    let mut out = std::collections::BTreeMap::new();
    out.insert("run".to_string(), keys(&v));
    let frame = &v["frames"][0];
    out.insert("frame".to_string(), keys(frame));
    out.insert("reading".to_string(), keys(&frame["readings"][0]));
    for panel in frame["panels"].as_array().expect("panels") {
        let kind = panel["kind"].as_str().expect("a panel names its kind");
        out.insert(format!("panel {kind}"), keys(panel));
        // The one nested object, and the same for every kind; the last one seen stands for all.
        if panel.get("place").is_some() {
            out.insert("place".to_string(), keys(&panel["place"]));
        }
    }
    out
}

/// What format 3 is, written down.
///
/// There is no row for 1 or 2. Format 1 had no `format` key at all and format 2 was not one shape:
/// files stamped 2 were written both before and after `surface`, `labels` and `bonds` existed,
/// which is the defect this table is here to stop. A reader of 3 reads all of them.
fn format_3() -> std::collections::BTreeMap<String, Vec<String>> {
    let row = |keys: &str| keys.split(' ').map(str::to_string).collect::<Vec<_>>();
    [
        ("run", "format frames title"),
        ("frame", "panels readings t"),
        ("reading", "domain label unit value"),
        ("place", "at_m turn"),
        (
            "panel field",
            "extent_m kind lattice name nx ny nz place unit values",
        ),
        (
            "panel paths",
            "bounds kind name place starts unit values vertices",
        ),
        (
            "panel points",
            "bonds bounds boxed kind labels name place positions unit values",
        ),
        (
            "panel surface",
            "bounds kind name place positions triangles unit values",
        ),
    ]
    .into_iter()
    .map(|(place, keys)| (place.to_string(), row(keys)))
    .collect()
}

/// One frame holding one panel of every shape the format has, each with every optional key present.
///
/// **The optional keys are written only when they differ from their default**, so a panel at the
/// origin has no `place` and a nodal field no `lattice`. The first version of this used both
/// defaults and the table below came out without either key — two keys the guard would never have
/// watched. So every panel is placed off the origin and the field is cell-centred.
fn every_shape() -> Vec<pantometry::scene::Frame> {
    use pantometry::scene::{Frame, Lattice, Panel, PanelData, Placed};
    let panel = |name: &str, data: PanelData| Panel {
        name: name.to_string(),
        unit: "K",
        place: Placed {
            at_m: [0.01, 0.0, 0.0],
            turn: [0.0, 0.0, 0.0, 1.0],
        },
        data,
    };
    vec![Frame {
        time_s: 0.0,
        panels: vec![
            panel(
                "field",
                PanelData::Field {
                    nx: 2,
                    ny: 1,
                    nz: 1,
                    lattice: Lattice::Centred,
                    values: vec![300.0, 301.0],
                    extent_m: [0.0, 0.0, 0.0, 0.01, 0.0, 0.0],
                },
            ),
            panel(
                "paths",
                PanelData::paths(vec![vec![[0.0; 3], [1.0, 0.0, 0.0]]], vec![1.0]),
            ),
            panel(
                "surface",
                PanelData::surface(
                    vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                    vec![[0, 1, 2]],
                    vec![1.0, 2.0, 3.0],
                ),
            ),
            panel(
                "points",
                PanelData::Points {
                    positions: vec![[0.0; 3], [1.0, 0.0, 0.0]],
                    values: vec![1.0, 2.0],
                    bounds: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                    boxed: false,
                    labels: vec!["A:1:GLY".to_string(), "A:2:ALA".to_string()],
                    bonds: vec![[0, 1]],
                },
            ),
        ],
        readings: vec![pantometry::core::Reading {
            domain: "field".to_string(),
            label: "peak".to_string(),
            value: 301.0,
            unit: "K",
        }],
    }]
}

/// **What each format was, key by key, so a shape that grows without its number is caught.**
///
/// `FORMAT` exists so that a reader meeting a newer file can say "upgrade" rather than "broken",
/// and it only does that job if the number moves whenever the shape does. It did not, twice: the
/// `surface` panel kind and the `labels` and `bonds` keys on `points` both shipped under format 2,
/// and a reader of format 2 — which is `deny_unknown_fields` on purpose — met them as an unknown
/// variant and an unknown field. Correct to refuse, and wrong about why.
///
/// This is a **pin** and not a closed form: the table is what the format *is*, written down.
/// Adding a key or a panel kind fails it, and the fix is two lines — bump `FORMAT` on both sides
/// and add the new row — rather than editing the old one, because a file stamped with an old
/// number has the old number's shape.
#[test]
fn the_format_number_moves_whenever_the_shape_does() {
    let json = pantometry::view::data::to_json("every shape", &every_shape());
    let seen = shape_of(&json);
    for (place, keys) in &seen {
        println!("  {place:<16} {}", keys.join(" "));
    }
    let format = pantometry::view::data::FORMAT;
    let want = match format {
        3 => format_3(),
        n => panic!(
            "FORMAT is {n} and there is no table of what format {n} is. Write it down beside \
             format_3, from what this test prints"
        ),
    };
    assert_eq!(
        seen, want,
        "the writer's shape is not format {format}'s. If a key or a panel kind was added, bump \
         FORMAT in pantometry-view and viewer-core and add a row for the new number — a file \
         stamped {format} has to keep format {format}'s shape"
    );
    let run =
        viewer_core::Run::from_json(&json).expect("the reader takes every shape the writer has");
    assert_eq!(run.frames[0].panels.len(), 4, "one panel of each kind");
}
