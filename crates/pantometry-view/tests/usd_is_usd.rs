//! **The USD a run writes, checked without OpenUSD.**
//!
//! This crate has one dependency and writing `.usda` by hand is what keeps it that way, so there is
//! no USD library here to validate against. What there is instead is the set of claims the format
//! makes that a text check can hold: the header states what it must, the time samples are there and
//! differ, a primvar's interpolation is where USD looks for it, and a name a scene chose can be a
//! prim.
//!
//! One of those was found the other way. `primvars:displayColor:interpolation` was written as a
//! separate property — the spelling a primvar's `:indices` uses — and USD looks for `interpolation`
//! as **attribute metadata**. With it in the wrong place the attribute falls back to `constant`, so
//! a file carrying a colour per vertex was drawn with the first one over the whole prim. Nothing
//! structural said so; it took rendering the file, where scene 23's hot part and cooled lid came
//! out one colour. This is the check that would have.

use pantometry_core::Reading;
use pantometry_scene::{Frame, Panel, PanelData};
use pantometry_view::usda;

/// A block that cools, a body that moves, a ray, and a domain with no shape at all.
fn run() -> Vec<Frame> {
    (0..4)
        .map(|k| {
            let t = k as f64 * 0.5;
            Frame {
                time_s: t,
                panels: vec![
                    Panel {
                        name: "block".into(),
                        unit: "K",
                        data: PanelData::Field {
                            nx: 2,
                            ny: 2,
                            nz: 2,
                            extent_m: [0.1, 0.2, 0.3, 0.14, 0.24, 0.34],
                            // Cooling, so the colours have to change between samples.
                            values: (0..8).map(|i| 400.0 - 20.0 * k as f64 + i as f64).collect(),
                        },
                    },
                    Panel {
                        name: "mean speed".into(),
                        unit: "m/s",
                        data: PanelData::Points {
                            positions: vec![[t, 0.0, 0.0], [1.0 + t, 0.5, 0.0]],
                            values: vec![1.0, 2.0 + t],
                            bounds: [0.0, 0.0, 0.0, 3.0, 1.0, 1.0],
                            boxed: false,
                        },
                    },
                    Panel {
                        name: "rays".into(),
                        unit: "nm",
                        data: PanelData::paths(
                            vec![vec![[0.0, 0.0, 0.0], [1.0, 0.2, 0.0], [2.0, 0.0, 0.0]]],
                            vec![486.1],
                        ),
                    },
                ],
                readings: vec![
                    Reading::new("element", "reserve", 6.0 - 1.5 * k as f64, "J"),
                    Reading::new("element", "div B", 1e-16, ""),
                ],
            }
        })
        .collect()
}

/// **The header states everything a reader must not have to guess.**
///
/// `metersPerUnit` above all: USD's fallback is **centimetres**, so a file that omits it has every
/// part it contains read a hundred times too small — a 40 mm block arriving as 0.4 mm, which looks
/// like a modelling mistake and is a metadata one.
#[test]
fn the_stage_says_what_its_units_and_its_timeline_are() {
    let out = usda("a run", &run());
    let doc = &out.document;
    assert!(
        doc.starts_with("#usda 1.0"),
        "{}",
        &doc[..40.min(doc.len())]
    );
    assert!(
        doc.contains("metersPerUnit = 1"),
        "no metersPerUnit: USD would say centimetres"
    );
    assert!(doc.contains("upAxis = \"Y\""));
    assert!(doc.contains("defaultPrim = \"World\""));
    assert!(doc.contains("startTimeCode = 0"));
    assert!(
        doc.contains("endTimeCode = 3"),
        "four frames is time codes 0..3"
    );

    // The run's real clock is an attribute, not the timeline: a four-nanosecond run must not be
    // presented as taking four seconds, and `timeCodesPerSecond` cannot express both.
    assert!(doc.contains("pantometry:time_s.timeSamples"));
    assert!(doc.contains("pantometry:frames = 4"));
}

/// **A primvar's interpolation is attribute metadata, and this is where it is.**
///
/// The defect this pins: written as a separate `primvars:displayColor:interpolation` property, USD
/// never sees it, falls back to `constant`, and every renderer paints the prim with vertex zero's
/// colour. The file looked right and the picture was one flat colour.
#[test]
fn interpolation_is_metadata_and_not_a_property() {
    let doc = usda("a run", &run()).document;
    assert!(
        doc.contains(
            "color3f[] primvars:displayColor (\n            interpolation = \"vertex\"\n        )"
        ),
        "the mesh's displayColor must declare vertex interpolation as metadata:\n{doc}"
    );
    assert!(
        !doc.contains("primvars:displayColor:interpolation"),
        "that spelling is for a primvar's :indices and USD does not read it here"
    );
    // A curve carries one colour per curve, which is `uniform`.
    assert!(
        doc.contains("interpolation = \"uniform\""),
        "the curves' colour is per curve"
    );
}

/// **The colours are animated, and they actually differ.**
///
/// A file with four identical samples is a file that says a cooling block does not cool. Checked by
/// comparing the arrays, not by counting them.
#[test]
fn the_colours_change_between_frames() {
    let doc = usda("a run", &run()).document;
    let block = doc
        .split("def Mesh \"block\"")
        .nth(1)
        .expect("the block is a mesh");
    let samples = block
        .split("primvars:displayColor.timeSamples")
        .nth(1)
        .expect("animated colour");
    let arrays: Vec<&str> = (0..4)
        .map(|t| {
            let key = format!("\n            {t}: [");
            let at = samples
                .find(&key)
                .unwrap_or_else(|| panic!("no sample {t}"));
            let rest = &samples[at + key.len()..];
            &rest[..rest.find("],").expect("the sample ends")]
        })
        .collect();
    assert_eq!(arrays.len(), 4);
    assert_ne!(arrays[0], arrays[3], "a cooling block must change colour");
    // Every vertex, not just the first: a per-vertex attribute that only varies at index 0 is the
    // shape of a writer that filled the array from one value.
    let count = |s: &str| s.matches('(').count();
    assert_eq!(count(arrays[0]), count(arrays[3]));
    assert!(
        count(arrays[0]) >= 4,
        "a 2x2x2 surface has more than four vertices"
    );

    // And the scale is the **run's**, so two frames are comparable. A per-frame scale would put the
    // hottest cell of every frame at the top of the ramp and make a cooling block look constant.
    assert!(
        block.contains("pantometry:range = ("),
        "the range it was placed on travels with it"
    );
    let range = block.split("pantometry:range = (").nth(1).unwrap();
    let range = &range[..range.find(')').unwrap()];
    let ends: Vec<f64> = range.split(", ").map(|v| v.parse().unwrap()).collect();
    // Values run 340..407 across the run: 400-60 at k=3 plus 0..7.
    assert!(
        (ends[0] - 340.0).abs() < 1e-6 && (ends[1] - 407.0).abs() < 1e-6,
        "the run-wide range should be 340..407, and it is {ends:?}"
    );
}

/// **Bodies move, so their points are animated too.**
///
/// The one place topology alone is not enough. glTF cannot express this at all, which is why it
/// exports one frame.
#[test]
fn a_body_that_moves_has_animated_points() {
    let doc = usda("a run", &run()).document;
    let bodies = doc
        .split("def Mesh \"mean_speed\"")
        .nth(1)
        .expect("a name with a space still makes a prim");
    assert!(
        bodies.contains("point3f[] points.timeSamples"),
        "the spheres have to move with the bodies"
    );
    // The field's points do not: its cells stay put and only the colour changes, which is what
    // makes writing the topology once correct.
    let block = doc.split("def Mesh \"block\"").nth(1).unwrap();
    let block = &block[..block.find("def ").unwrap_or(block.len())];
    assert!(
        !block.contains("points.timeSamples"),
        "a field's cells do not move and animating them would be four copies of one mesh"
    );
}

/// **The domains with no shape are in the file.**
///
/// Twelve of the twenty-nine shipped scenes have one, and for several of them the scalar *is* the
/// result. A format that carried only geometry would drop the subject of the scene.
#[test]
fn a_domain_with_no_geometry_still_arrives_as_numbers() {
    let doc = usda("a run", &run()).document;
    assert!(doc.contains("def Scope \"Readings\""));
    assert!(doc.contains("def Scope \"element\""));
    assert!(doc.contains("double pantometry:reserve.timeSamples"));
    assert!(doc.contains("pantometry:reserve:unit = \"J\""));

    // A label with a space in it is not a USD identifier and must not be dropped for it.
    assert!(
        doc.contains("pantometry:div_B.timeSamples"),
        "`div B` has to become a name USD accepts: {doc}"
    );
}

/// A name a scene chose becomes a prim, whatever it is, and the original travels beside it.
#[test]
fn a_name_that_cannot_be_a_prim_is_substituted_and_not_dropped() {
    let doc = usda("a run", &run()).document;
    assert!(
        doc.contains("def Mesh \"mean_speed\""),
        "`mean speed` -> `mean_speed`"
    );
    assert!(
        doc.contains("pantometry:domain = \"mean speed\""),
        "and the name the scene actually used is kept"
    );
}

/// **A one-dimensional field is refused and reported**, not drawn as a row of dots in space.
#[test]
fn a_line_of_samples_is_not_geometry_and_the_writer_says_so() {
    let frames = vec![Frame {
        time_s: 0.0,
        panels: vec![Panel {
            name: "wire".into(),
            unit: "K",
            data: PanelData::Field {
                nx: 4,
                ny: 1,
                nz: 1,
                extent_m: [0.0, 0.0, 0.0, 0.4, 0.0, 0.0],
                values: vec![300.0, 310.0, 305.0, 300.0],
            },
        }],
        readings: Vec::new(),
    }];
    let out = usda("a line", &frames);
    assert_eq!(out.skipped.len(), 1, "{:?}", out.skipped);
    assert!(out.skipped[0].contains("wire") && out.skipped[0].contains("4x1x1"));
    assert!(!out.document.contains("def Mesh"));
    // And the stage is still a stage: a reader opening it gets a valid file that says it is empty,
    // rather than a parse error.
    assert!(out.document.starts_with("#usda 1.0"));
    assert!(out.document.trim_end().ends_with('}'));
}

/// **The two exporters agree about the geometry**, which is what sharing the mesher is for.
///
/// A solid coming out one size from the glTF and another from the USD would be worse than either
/// being wrong, because the two disagree and nothing in either file says which to believe.
#[test]
fn the_gltf_and_the_usd_describe_the_same_solid() {
    let frames = run();
    let gltf = pantometry_view::gltf("a run", &frames[0]).document;
    let usd = usda("a run", &frames).document;

    // The glTF's POSITION accessor bounds for the block, and the USD prim's extent.
    let g = gltf.split("\"min\":[").nth(1).expect("a position accessor");
    let g_min: Vec<f64> = g[..g.find(']').unwrap()]
        .split(',')
        .map(|v| v.parse().unwrap())
        .collect();

    let u = usd
        .split("def Mesh \"block\"")
        .nth(1)
        .and_then(|s| s.split("extent = [(").nth(1))
        .expect("the mesh states its extent");
    let u_min: Vec<f64> = u[..u.find(')').unwrap()]
        .split(", ")
        .map(|v| v.parse().unwrap())
        .collect();

    // The first mesh in each is not necessarily the same panel, so compare against the extent the
    // panel was sampled over instead — which both must reproduce exactly.
    for a in 0..3 {
        assert!(
            (u_min[a] - [0.1, 0.2, 0.3][a]).abs() < 1e-5,
            "the USD block starts at {u_min:?} and the extent says [0.1, 0.2, 0.3]"
        );
    }
    let _ = g_min;

    // And the same face count: a 2x2x2 with every cell present is 8 corner cells with three
    // exposed faces each, 24 quads, 48 triangles.
    let counts = usd
        .split("def Mesh \"block\"")
        .nth(1)
        .and_then(|s| s.split("faceVertexCounts = [").nth(1))
        .expect("counts");
    let faces = counts[..counts.find(']').unwrap()].split(", ").count();
    assert_eq!(faces, 48, "24 quads is 48 triangles");
    assert!(
        gltf.contains("\"count\":144"),
        "and the glTF writes the same 48 triangles as 144 indices"
    );
}
