//! What each writer does with a sample that is not a number.
//!
//! A field gained cells with nothing in them, and every writer here had an answer for them that
//! looked like an answer for something. Measured on a real scene before the fix: the glTF carried
//! 252 points for a grid with 120 solid cells, so two thirds of what went to Blender was material
//! that is not there; the SVG drew the empty ones in the **middle** colour of the ramp, because
//! `NaN as i32` saturates to zero and zero is the middle bucket, so an empty cell was the same
//! colour as a cell at the average; and the JSON writers would have emitted the literal `NaN`,
//! which is not a JSON token at all.
//!
//! One rule, four writers: **absence is not a value.** These are on the writers. The sampler that
//! decides what is absent has its own, in `pantometry-thermal`.

use pantometry_scene::{Frame, Panel, PanelData};

/// A 3x1x3 field with the middle *column* empty, so the hole survives any projection.
///
/// The values are `300 + 10k`, distinct per slice, so an assertion can say where in the array a
/// hole is and not merely how many there are.
fn holed() -> Frame {
    let mut values = Vec::new();
    for k in 0..3 {
        for i in 0..3 {
            values.push(if i == 1 {
                f64::NAN
            } else {
                300.0 + 10.0 * k as f64
            });
        }
    }
    Frame {
        time_s: 0.0,
        panels: vec![Panel {
            name: "block".into(),
            unit: "K",
            data: PanelData::Field {
                nx: 3,
                ny: 1,
                nz: 3,
                extent_m: [0.0, 0.0, 0.0, 0.06, 0.0, 0.06],
                values,
            },
        }],
        readings: Vec::new(),
    }
}

/// **The glTF has surface where the material is and a hole where it is not.**
///
/// This is the one output that leaves the workspace — Blender, three.js, a USD tool — so a
/// clearance exported as geometry becomes a wrong model in somebody else's software, where nothing
/// about it says where it came from.
///
/// Checked **geometrically** rather than by counting vertices, because the count is a property of
/// the meshing and the claim is a property of the object: no vertex may lie inside the empty
/// column, and there must be faces on **both** of its walls, because a hole through a solid has
/// two sides and a mesher that only wrote the outside would leave the block hollow-looking and
/// closed.
///
/// The fixture is 3x1x3 over 60 mm x 60 mm with the middle column empty, so the samples are at
/// x = 0, 30 and 60 mm, the cells are 30 mm wide, and the gap is the open interval
/// (15 mm, 45 mm).
#[test]
fn the_gltf_leaves_the_holes_out() {
    let out = pantometry_view::gltf("holed", &holed());
    let doc = &out.document;

    // Positions come back out of the base64 rather than being trusted from the source: this is a
    // test about what a reader will load, and the encoder is between here and there.
    let xs = exported_x(doc);
    assert!(!xs.is_empty(), "the export has no geometry at all: {doc}");
    let inside: Vec<f32> = xs
        .iter()
        .copied()
        .filter(|x| *x > 0.0151 && *x < 0.0449)
        .collect();
    assert!(
        inside.is_empty(),
        "{} vertices sit inside the empty column: {inside:?}",
        inside.len()
    );

    // Both walls of the gap. Within a tenth of a millimetre, which is float32 at this magnitude
    // with room, and far below the 30 mm cell the assertion is about.
    let near = |t: f32| xs.iter().any(|x| (x - t).abs() < 1e-4);
    assert!(near(0.015), "no face on the near wall of the gap: {xs:?}");
    assert!(near(0.045), "no face on the far wall of the gap: {xs:?}");

    // And it is a surface, with normals, rather than the point cloud this used to write. A point
    // cloud takes no light and casts no shadow, so a clearance in one is invisible either way.
    assert!(doc.contains("\"NORMAL\""), "no normals: {doc}");
    assert!(doc.contains("\"mode\":4"), "not triangles: {doc}");

    // Six present cells: four with one exposed z face and two with none, all six exposed on both
    // x and both y, which is 28 quads and 168 indices. Written out because a mesher that stopped
    // culling interior faces would still pass every assertion above.
    assert!(
        doc.contains("\"count\":168"),
        "28 quads is 168 indices, and this is not that: {doc}"
    );
}

/// The x coordinate of every exported vertex, decoded from the document's own buffer.
fn exported_x(doc: &str) -> Vec<f32> {
    let key = "base64,";
    let from = doc.find(key).expect("a buffer") + key.len();
    let to = from + doc[from..].find('"').expect("the buffer ends");
    let bytes = un_base64(&doc[from..to]);

    // The POSITION accessor's view: the first VEC3 accessor, which the writer emits first.
    let vi = doc
        .find("\"bufferView\":")
        .and_then(|i| doc[i + 13..].split(',').next()?.parse::<usize>().ok())
        .expect("a view index");
    let view = doc
        .match_indices("{\"buffer\":0,")
        .nth(vi)
        .map(|(i, _)| &doc[i..])
        .expect("the view");
    let field = |name: &str| -> usize {
        let k = view.find(name).expect(name) + name.len() + 1;
        view[k..]
            .split(|c: char| !c.is_ascii_digit())
            .next()
            .unwrap()
            .parse()
            .unwrap()
    };
    let (offset, length) = (field("\"byteOffset\""), field("\"byteLength\""));

    bytes[offset..offset + length]
        .chunks_exact(4)
        .enumerate()
        .filter(|(i, _)| i % 3 == 0)
        .map(|(_, c)| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect()
}

/// Standard base64 back to bytes. The writer has its own encoder and no decoder; this is the
/// decoder, and it lives in the test because only a test needs to read what was written.
fn un_base64(s: &str) -> Vec<u8> {
    let val = |c: u8| -> Option<u32> {
        Some(match c {
            b'A'..=b'Z' => (c - b'A') as u32,
            b'a'..=b'z' => (c - b'a') as u32 + 26,
            b'0'..=b'9' => (c - b'0') as u32 + 52,
            b'+' => 62,
            b'/' => 63,
            _ => return None,
        })
    };
    let mut out = Vec::with_capacity(s.len() / 4 * 3);
    for chunk in s.as_bytes().chunks(4) {
        let mut n = 0u32;
        let mut got = 0;
        for &c in chunk {
            n <<= 6;
            if let Some(v) = val(c) {
                n |= v;
                got += 1;
            }
        }
        out.push((n >> 16) as u8);
        if got > 2 {
            out.push((n >> 8) as u8);
        }
        if got > 3 {
            out.push(n as u8);
        }
    }
    out
}

/// **The JSON writers say `null`, in the right places, and never say `NaN`.**
///
/// `NaN` is not a token this format has. A document containing one is refused outright by a strict
/// parser — which turns a hole in a field into a file nobody can open — and accepted by a lenient
/// one, which is worse, because it comes back as a float and gets plotted.
///
/// Positional rather than a count: the array has to be `300, nothing, 300` and not three nulls
/// swept to the end, because a writer that dropped the empties would shift every later cell.
#[test]
fn the_json_writers_spell_absence_null_and_leave_it_where_it_was() {
    for (what, text) in [
        ("to_json", pantometry_view::to_json("holed", &[holed()])),
        ("html", pantometry_view::html("holed", &[holed()])),
    ] {
        // Spelled as it would appear *in an array*, because the report's own script names the
        // JavaScript constant legitimately — a bare search for `NaN` matches that and says the
        // writer is broken when it is not.
        assert!(
            !text.contains(",NaN") && !text.contains("[NaN") && !text.contains("NaN]"),
            "{what}: `NaN` is not a JSON token and must not reach a values array"
        );
        // `3e2` rather than `3.000000e2`: both writers share one encoder now, and it drops
        // trailing zeros that say nothing — the same number, 8.3% smaller on the largest report
        // this workspace produces. What matters here is unchanged: the hole is between the two
        // values it was between, and not swept to the end of the array.
        assert!(
            text.contains("3e2,null,3e2")
                && text.contains("3.1e2,null,3.1e2")
                && text.contains("3.2e2,null,3.2e2"),
            "{what}: the hole stays in the slot it was in"
        );
    }
}

/// **The SVG draws no square where there is nothing**, and the background shows through.
///
/// The failure was silent in the exact way that matters: nothing about a mid-ramp square says it is
/// not a measurement. Counted by the squares the raster emits — the filmstrip draws the middle
/// z-slice of a volume, which here is three cells with one of them empty.
#[test]
fn the_svg_draws_no_square_for_an_empty_cell() {
    let svg = pantometry_view::svg("holed", &[holed()], 1);
    assert_eq!(
        svg.matches("h1v1h-1z").count(),
        2,
        "one square per cell that has a value, and none for the one that does not"
    );
}
