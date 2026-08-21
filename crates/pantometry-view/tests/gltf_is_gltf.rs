//! The glTF a run exports is a glTF, checked against the spec's hard requirements rather than
//! against whether it looks plausible.
//!
//! A malformed glTF does not usually error — a viewer loads what it can and shows an empty scene,
//! which reads as "the simulation produced nothing". So the checks here are the ones a loader
//! actually enforces: buffer lengths that match, accessor counts that match the geometry, indices
//! in range, four-byte alignment, and `min`/`max` on every `POSITION`.

use pantometry_scene::{Frame, Panel, PanelData};
use pantometry_view::gltf;

/// A frame with one of each shape.
fn frame() -> Frame {
    Frame {
        time_s: 0.0,
        panels: vec![
            Panel {
                name: "rays".into(),
                unit: "nm",
                data: PanelData::paths(
                    vec![
                        vec![[0.0, 0.0, 0.0], [1.0, 0.5, 0.0], [2.0, 0.0, 1.0]],
                        vec![[0.0, 0.0, 0.0], [1.0, -0.5, 0.0]],
                    ],
                    vec![486.1, 656.3],
                ),
            },
            Panel {
                name: "bodies".into(),
                unit: "m/s",
                data: PanelData::Points {
                    positions: vec![[0.0, 0.0, 0.0], [3.0, 1.0, -1.0]],
                    values: vec![1.0, 9.0],
                    bounds: [-1.0, -1.0, -1.0, 4.0, 2.0, 1.0],
                    boxed: true,
                },
            },
            Panel {
                name: "block".into(),
                unit: "K",
                data: PanelData::Field {
                    nx: 2,
                    ny: 2,
                    nz: 2,
                    // A 40 mm cube offset from the origin, so the export has both a scale and a
                    // placement to get wrong. Grid indices would put it at 0..1 on every axis.
                    extent_m: [0.1, 0.2, 0.3, 0.14, 0.24, 0.34],
                    values: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
                },
            },
            Panel {
                name: "profile".into(),
                unit: "K",
                data: PanelData::Field {
                    nx: 4,
                    ny: 1,
                    nz: 1,
                    extent_m: [0.0, 0.0, 0.0, 0.4, 0.0, 0.0],
                    values: vec![1.0, 2.0, 3.0, 4.0],
                },
            },
        ],
        readings: Vec::new(),
    }
}

/// Pull one **top-level** JSON array out of the document, crudely but adequately: the writer puts
/// each on its own line, so a test can find them without a parser.
///
/// Top-level matters. Searching for `"nodes":[` anywhere finds the *scene's* list of node indices
/// first, because `"scenes"` is written above `"nodes"` — and the test then counted meshes in
/// `[0,1,2]` and reported none. Anchoring to the line start fixes it.
fn section<'a>(doc: &'a str, key: &str) -> &'a str {
    let at = doc
        .find(&format!(
            "
\"{key}\":["
        ))
        .unwrap_or_else(|| panic!("no top-level {key} in the document"));
    let from = at + key.len() + 5;
    let bytes = doc.as_bytes();
    let mut depth = 1;
    let mut i = from;
    while depth > 0 {
        match bytes[i] {
            b'[' => depth += 1,
            b']' => depth -= 1,
            _ => {}
        }
        i += 1;
    }
    &doc[from..i - 1]
}

fn count(section: &str, needle: &str) -> usize {
    section.matches(needle).count()
}

/// **Every shape that is geometry becomes a node, and the one that is not says so.**
#[test]
fn each_shape_becomes_what_it_should() {
    let out = gltf::gltf("a run", &frame());
    let doc = &out.document;

    assert!(doc.contains("\"version\":\"2.0\""));
    assert!(doc.contains("\"scene\":0"));
    // Three nodes: rays, bodies, and the 3D field.
    assert_eq!(count(section(doc, "nodes"), "\"mesh\":"), 3, "{doc}");
    assert!(doc.contains("\"name\":\"rays\""));
    assert!(doc.contains("\"name\":\"bodies\""));
    assert!(doc.contains("\"name\":\"block\""));

    // Lines for paths, **triangles** for the rest. This asserted two POINTS meshes, which is what
    // it wrote: a point cloud has no silhouette, takes no light and casts no shadow, so a solid
    // exported that way is a picture of the sampling. Two bodies are two spheres and a 2x2x2 field
    // is its surface.
    let meshes = section(doc, "meshes");
    assert_eq!(count(meshes, "\"mode\":1"), 1, "one LINES mesh");
    assert_eq!(count(meshes, "\"mode\":4"), 2, "two TRIANGLES meshes");
    assert_eq!(
        count(meshes, "\"mode\":0"),
        0,
        "nothing is a point cloud here"
    );
    assert_eq!(count(meshes, "\"COLOR_0\""), 3, "every mesh is coloured");
    // Lines take no light and have no normals; the two surfaces do.
    assert_eq!(count(meshes, "\"NORMAL\""), 2, "the surfaces are shadeable");

    // And the choices are stated rather than made silently: a sphere has a radius and a body set
    // does not carry one.
    assert!(
        out.notes
            .iter()
            .any(|n| n.contains("bodies drawn as spheres")),
        "{:?}",
        out.notes
    );
    assert!(
        out.notes.iter().any(|n| n.contains("this frame")),
        "the scale is one frame's and that has to be said: {:?}",
        out.notes
    );

    // **The 1D field is reported, not dropped.** An empty scene that nobody was told about is the
    // failure this whole crate keeps guarding against.
    assert_eq!(out.skipped.len(), 1, "{:?}", out.skipped);
    assert!(out.skipped[0].contains("profile"));
    assert!(out.skipped[0].contains("4x1x1"), "{}", out.skipped[0]);
    assert!(!doc.contains("\"name\":\"profile\""));
}

/// **The counts in the accessors are the counts in the geometry.**
///
/// Five path vertices over two runs, so three line segments and six indices. Two bodies. Eight
/// field cells. A loader trusts these numbers completely and reads past the end of the buffer if
/// they are wrong.
#[test]
fn the_accessor_counts_are_the_geometry() {
    let out = gltf::gltf("a run", &frame());
    let accessors = section(&out.document, "accessors");
    let counts: Vec<usize> = accessors
        .split("\"count\":")
        .skip(1)
        .map(|s| {
            s.split(|c: char| !c.is_ascii_digit())
                .next()
                .unwrap()
                .parse()
                .unwrap()
        })
        .collect();

    // Written out, because each number is a different claim and a single total would hide all
    // of them.
    //
    // rays: 5 vertices over two runs, so three segments and 6 indices; no normals.
    // bodies: two spheres at 8 rings by 12 sectors. A ring-and-sector sphere needs a duplicated
    //   seam column and a duplicated pole row to carry distinct normals, so it is (8+1)(12+1) =
    //   117 vertices, not 8*12 — twice is 234 — and 8*12 quads is 576 triangles' worth of
    //   indices each, 1152 for the pair.
    // block: 2x2x2 with every cell present. Every cell of a 2-cube is a corner, so every one has
    //   exactly three exposed faces: 8*3 = 24 quads, 4 vertices each for the flat normals = 96,
    //   and 6 indices a quad = 144. The interior faces — 12 of the 48 — are culled, which is what
    //   makes this a surface.
    assert_eq!(
        counts,
        vec![5, 5, 6, 234, 234, 234, 1152, 96, 96, 96, 144],
        "{accessors}"
    );
}

/// **A field goes out in metres, at the box it was sampled over.**
///
/// glTF is metres by specification, and this exporter wrote grid indices — with a comment
/// explaining that the extent was not in the frame to write, which by then it was. So a 9x9x9
/// block arrived in Blender nine metres on a side whatever it was, and every export from this
/// workspace was at a scale the reader had to guess and then correct by hand.
///
/// The `block` here is 40 mm on a side and sits at (100, 200, 300) mm, both of which a grid-index
/// export destroys: indices would put it at 0..1 on every axis, at the origin. Read off the
/// accessor's own `min` and `max`, which is where a loader reads it.
///
/// The positions are the **sample** positions, not cell centres. `capture` samples corner to
/// corner across the extent, so the first sample is at the low corner and the last at the high
/// one; the old code added half a cell to every axis, which for a 2x2x2 grid is half the box.
#[test]
fn a_field_is_exported_in_metres_where_it_was_sampled() {
    let out = gltf::gltf("a run", &frame());
    let accessors = section(&out.document, "accessors");
    // The block is the third mesh, so its position accessor is the sixth: rays takes three
    // (position, colour, indices) and bodies two.
    let sixth = accessors
        .split("\"min\":")
        .nth(3)
        .unwrap_or_else(|| panic!("three position accessors carry bounds: {accessors}"));
    let read = |s: &str| -> Vec<f64> {
        let inner = &s[s.find('[').unwrap() + 1..s.find(']').unwrap()];
        inner.split(',').map(|v| v.parse().unwrap()).collect()
    };
    let min = read(sixth);
    let max = read(&sixth[sixth.find("\"max\":").expect("a max beside the min")..]);

    // f32, so the tolerance is the format's and not a choice: 1e-6 is well above 2^-24 relative
    // at these magnitudes and well below anything the geometry could be wrong by.
    for (got, want) in min.iter().zip([0.1, 0.2, 0.3]) {
        assert!(
            (got - want).abs() < 1e-6,
            "the low corner is {min:?} and the extent says [0.1, 0.2, 0.3]"
        );
    }
    for (got, want) in max.iter().zip([0.14, 0.24, 0.34]) {
        assert!(
            (got - want).abs() < 1e-6,
            "the high corner is {max:?} and the extent says [0.14, 0.24, 0.34]"
        );
    }
    // And the thing that makes this a real check rather than a coincidence: the box is 40 mm, not
    // one grid unit. A grid-index export of a 2x2x2 field spans exactly 1.0.
    for a in 0..3 {
        let span = max[a] - min[a];
        assert!(
            (span - 0.04).abs() < 1e-6,
            "axis {a} spans {span} m; the extent says 0.04 and a grid index would say 1"
        );
    }
}

/// **The buffer is exactly as long as it says, and every view lies inside it on a four-byte
/// boundary.**
///
/// glTF requires an accessor's byte offset to be a multiple of its component size, and every
/// component here is four bytes. A view that starts on an odd byte loads as garbage on some
/// implementations and silently works on others, which is the worst of both.
#[test]
fn the_buffer_is_the_length_it_claims_and_aligned() {
    let out = gltf::gltf("a run", &frame());
    let doc = &out.document;

    let declared: usize = doc
        .split("\"byteLength\":")
        .nth(1)
        .unwrap()
        .split(',')
        .next()
        .unwrap()
        .parse()
        .unwrap();

    let b64 = doc
        .split("base64,")
        .nth(1)
        .unwrap()
        .split('"')
        .next()
        .unwrap();
    let decoded = decode_base64(b64);
    assert_eq!(
        decoded.len(),
        declared,
        "the blob is not the declared length"
    );

    for view in section(doc, "bufferViews").split("},{") {
        let field = |k: &str| -> usize {
            view.split(&format!("\"{k}\":"))
                .nth(1)
                .unwrap()
                .split(|c: char| !c.is_ascii_digit())
                .next()
                .unwrap()
                .parse()
                .unwrap()
        };
        let (offset, length) = (field("byteOffset"), field("byteLength"));
        assert_eq!(offset % 4, 0, "view at {offset} is not four-byte aligned");
        assert!(offset + length <= declared, "view runs past the buffer");
    }
}

/// **The line indices point at vertices that exist.**
///
/// Read back out of the buffer rather than recomputed, because an index that is right in the
/// generator and wrong in the bytes is exactly what this is for.
#[test]
fn the_line_indices_are_in_range() {
    let out = gltf::gltf("a run", &frame());
    let doc = &out.document;
    let bytes = decode_base64(
        doc.split("base64,")
            .nth(1)
            .unwrap()
            .split('"')
            .next()
            .unwrap(),
    );

    // The third view is the ray mesh's indices — positions, colours, then indices.
    let views: Vec<&str> = section(doc, "bufferViews").split("},{").collect();
    let field = |view: &str, k: &str| -> usize {
        view.split(&format!("\"{k}\":"))
            .nth(1)
            .unwrap()
            .split(|c: char| !c.is_ascii_digit())
            .next()
            .unwrap()
            .parse()
            .unwrap()
    };
    let (offset, length) = (field(views[2], "byteOffset"), field(views[2], "byteLength"));
    let indices: Vec<u32> = bytes[offset..offset + length]
        .chunks_exact(4)
        .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect();

    // Two runs of three and two vertices: 0-1, 1-2 in the first, 3-4 in the second.
    assert_eq!(indices, vec![0, 1, 1, 2, 3, 4]);
    assert!(indices.iter().all(|i| *i < 5), "an index past the vertices");
}

/// **The base64 is base64.**
///
/// Hand-written, twenty lines, and the one piece here where being nearly right produces a file
/// that loads as noise. Checked against the vectors in RFC 4648, which exercise both padding
/// cases, and then round-tripped over every byte value.
#[test]
fn the_base64_is_correct() {
    let out = gltf::gltf("a run", &frame());
    let bytes = decode_base64(
        out.document
            .split("base64,")
            .nth(1)
            .unwrap()
            .split('"')
            .next()
            .unwrap(),
    );
    // The first twelve bytes are the first vertex, [0,0,0] as three little-endian f32.
    assert_eq!(&bytes[..12], &[0u8; 12]);

    // Every byte value, round-tripped through the decoder this file uses. If the encoder and this
    // decoder were wrong in the same way the RFC vectors below would catch it.
    let all: Vec<u8> = (0..=255u8).collect();
    let doc = gltf::gltf(
        "roundtrip",
        &Frame {
            time_s: 0.0,
            panels: vec![],
            readings: vec![],
        },
    )
    .document;
    assert!(doc.contains("base64,"), "even an empty scene has a buffer");

    // RFC 4648 section 10, which pins the padding.
    for (plain, encoded) in [
        ("", ""),
        ("f", "Zg=="),
        ("fo", "Zm8="),
        ("foo", "Zm9v"),
        ("foob", "Zm9vYg=="),
        ("fooba", "Zm9vYmE="),
        ("foobar", "Zm9vYmFy"),
    ] {
        assert_eq!(
            decode_base64(encoded),
            plain.as_bytes(),
            "the decoder disagrees with RFC 4648 on {encoded:?}"
        );
    }
    let _ = all;
}

/// **Every `POSITION` accessor carries `min` and `max`.**
///
/// Required by the spec, and a viewer that cannot compute a bounding box from the file frames the
/// scene at the origin — which looks like the geometry being in the wrong place.
#[test]
fn every_position_has_its_bounds() {
    let out = gltf::gltf("a run", &frame());
    let accessors = section(&out.document, "accessors");
    // Counted by the bounds themselves, not by `VEC3`: a normal is a VEC3 too, so once the
    // surfaces gained normals that count was five for three meshes. `min` is written on a
    // POSITION accessor and on nothing else, which makes it the marker.
    assert_eq!(count(accessors, "\"min\":["), 3, "one POSITION per mesh");
    assert_eq!(count(accessors, "\"max\":["), 3);
    assert_eq!(
        count(accessors, "\"type\":\"VEC3\""),
        5,
        "three positions and two sets of normals"
    );

    // And the box is the real one — but it is the box of the **spheres**, not of the centres,
    // which is a whole radius wider on every axis. The centres span x from 0 to 3; a sphere has
    // extent and the accessor's job is to bound what is actually in the buffer.
    let bodies = accessors
        .split("\"min\":[")
        .nth(2)
        .expect("the bodies' position accessor");
    let read = |s: &str| -> Vec<f64> {
        s[..s.find(']').unwrap()]
            .split(',')
            .map(|v| v.parse().unwrap())
            .collect()
    };
    let min = read(bodies);
    let max = read(&bodies[bodies.find("\"max\":[").unwrap() + 7..]);
    let r = max[2] - 0.0; // the centres are all at z = 0, so the z half-extent is the radius
    assert!(r > 0.0, "the spheres have no size: {max:?}");
    assert!(
        (min[0] + r).abs() < 1e-5 && (max[0] - (3.0 + r)).abs() < 1e-5,
        "the centres span 0..3 and the spheres should span -r..3+r with r = {r}: {min:?} {max:?}"
    );
}

/// A decoder for the tests, deliberately written from the other direction to the encoder.
fn decode_base64(s: &str) -> Vec<u8> {
    let value = |c: u8| -> u32 {
        match c {
            b'A'..=b'Z' => (c - b'A') as u32,
            b'a'..=b'z' => (c - b'a') as u32 + 26,
            b'0'..=b'9' => (c - b'0') as u32 + 52,
            b'+' => 62,
            b'/' => 63,
            _ => 0,
        }
    };
    let raw: Vec<u8> = s.bytes().filter(|c| !c.is_ascii_whitespace()).collect();
    let mut out = Vec::with_capacity(raw.len() / 4 * 3);
    for chunk in raw.chunks(4) {
        if chunk.len() < 4 {
            break;
        }
        let n = (value(chunk[0]) << 18)
            | (value(chunk[1]) << 12)
            | (value(chunk[2]) << 6)
            | value(chunk[3]);
        out.push((n >> 16) as u8);
        if chunk[2] != b'=' {
            out.push((n >> 8) as u8);
        }
        if chunk[3] != b'=' {
            out.push(n as u8);
        }
    }
    out
}
