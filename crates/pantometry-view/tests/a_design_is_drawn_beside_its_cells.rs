//! **The exporters draw the shape somebody designed, beside the cells it became.**
//!
//! Every export of a designed part showed its *rasterisation* — a staircase on the solver's grid —
//! and the thing it was rasterised from sat on disk beside the scene, reachable by nothing that
//! writes a file. `pantometry_world::Rasterised` has reported the volume error as a number since
//! designed parts existed, and a number in a terminal is not what a designer looks at.
//!
//! # It is not in the run, and that is the decision
//!
//! A run is the simulation's **output**; an STL is its input. So [`Drawing::designed`] is an
//! argument to the writer rather than a fourth `PanelData`, and `mesh::Designed`'s own doc lists
//! the three ways of putting it in the run instead and why each is worse. What is checked here is
//! that the two arrive together and stay distinguishable.
//!
//! # Uncoloured, and that is a claim about honesty
//!
//! The solver ran on cells. Tinting a smooth surface from the field would look like a
//! higher-resolution answer than the one computed — the same mistake as renormalising a colour
//! scale per frame, wearing different clothes. So the design gets one flat grey, and the test
//! below asserts it carries **no** variation rather than asserting its exact shade.

use pantometry_scene::{Frame, Panel, PanelData, Placed};
use pantometry_view::mesh::{self, Designed, Drawing, Surfaces};
use pantometry_view::{gltf_with, usda_with};

/// A 2×2×2 field over a 40 mm box, every cell full, so its boundary is the whole box.
fn cells() -> Frame {
    Frame {
        time_s: 0.0,
        readings: Vec::new(),
        panels: vec![Panel {
            name: "part".into(),
            unit: "K",
            place: Placed::HERE,
            data: PanelData::Field {
                nx: 2,
                ny: 2,
                nz: 2,
                extent_m: [0.0, 0.0, 0.0, 0.04, 0.04, 0.04],
                // Distinct, so the field's own colours vary and the design's must not.
                values: (0..8).map(|i| 300.0 + 10.0 * i as f64).collect(),
            },
        }],
    }
}

/// A unit tetrahedron of side `s`, as four triangles — the smallest closed shape there is.
fn tetra(s: f64) -> Vec<[[f64; 3]; 3]> {
    let (a, b, c, d) = ([0.0, 0.0, 0.0], [s, 0.0, 0.0], [0.0, s, 0.0], [0.0, 0.0, s]);
    vec![[a, c, b], [a, b, d], [a, d, c], [b, c, d]]
}

fn designed(place: Placed, s: f64) -> Designed {
    Designed {
        name: "part/parts[0]".into(),
        domain: "part".into(),
        place,
        surface: mesh::mesh_surface(&tetra(s), 0),
    }
}

/// The `n`th `"name":"..."` in the document, in order.
fn names(doc: &str) -> Vec<String> {
    doc.match_indices("\"name\":\"")
        .map(|(i, m)| {
            let rest = &doc[i + m.len()..];
            rest[..rest.find('"').expect("a closed name")].to_string()
        })
        .collect()
}

#[test]
fn nothing_is_drawn_that_was_not_asked_for() {
    // The default carries no designed parts, so an export takes the shape it always had. Checked
    // against the **panels**, not against another call with the default — comparing a default with
    // itself is a test that passes when the default is redefined, which this file's sibling
    // `an_export_can_be_a_level.rs` found out the hard way.
    let out = gltf_with("a run", &cells(), &Drawing::default());
    assert!(out.skipped.is_empty(), "{:?}", out.skipped);
    let named = names(&out.document);
    assert!(
        named.iter().filter(|n| n.as_str() == "part").count() >= 1,
        "the panel is missing: {named:?}"
    );
    assert!(
        !named.iter().any(|n| n.contains("parts[")),
        "a designed part appeared with none given: {named:?}"
    );
    assert!(
        !out.notes.iter().any(|n| n.contains("designed part")),
        "{:?}",
        out.notes
    );
}

#[test]
fn a_design_arrives_beside_its_cells_in_both_writers() {
    let drawing = Drawing::of(Surfaces::Boundary).and(vec![designed(Placed::HERE, 0.04)]);

    let gl = gltf_with("a run", &cells(), &drawing).document;
    let named = names(&gl);
    assert!(named.iter().any(|n| n == "part"), "{named:?}");
    assert!(named.iter().any(|n| n == "part/parts[0]"), "{named:?}");

    let us = usda_with("a run", &[cells()], &drawing).document;
    assert!(us.contains("def Mesh \"part\""), "the cells are missing");
    // USD identifiers are `[A-Za-z_][A-Za-z0-9_]*`, so the site is substituted — and the original
    // is kept beside it, which is what `pantometry:designed` is for.
    assert!(
        us.contains("def Mesh \"part_parts_0_\""),
        "the design is missing"
    );
    assert!(
        us.contains("uniform string pantometry:designed = \"part/parts[0]\""),
        "the name a finding uses is not in the file"
    );
}

#[test]
fn the_design_carries_no_value_where_the_cells_carry_eight() {
    // **The honesty claim, as an assertion.** The field's vertices differ from one another because
    // they hold values; the design's do not, because it holds none. Written as "no variation"
    // rather than "this exact grey" on purpose: what matters is that nothing about the picture
    // invites reading a temperature off a shape that has none.
    let drawing = Drawing::of(Surfaces::Boundary).and(vec![designed(Placed::HERE, 0.04)]);
    let us = usda_with("a run", &[cells()], &drawing).document;

    let design = &us[us.find("def Mesh \"part_parts_0_\"").expect("the design")..];
    let block = &design[..design.find("    }").expect("a closed prim")];
    assert!(
        block.contains("interpolation = \"constant\""),
        "the design's colour is per-vertex: {block}"
    );
    assert!(
        !block.contains("timeSamples"),
        "the design's colour animates, and a design does not change during a run"
    );
    // One colour, and it is the grey the editor's viewport uses, so the two agree about what a
    // designed part looks like.
    let want = format!(
        "[({}, {}, {})]",
        mesh::DESIGNED_GREY[0],
        mesh::DESIGNED_GREY[1],
        mesh::DESIGNED_GREY[2]
    );
    assert!(block.contains(&want), "expected {want} in: {block}");

    // The cells, by contrast, animate a colour per vertex.
    let cells_prim = &us[us.find("def Mesh \"part\"").expect("the cells")..];
    assert!(cells_prim.contains("interpolation = \"vertex\""));
}

#[test]
fn a_placed_design_says_where_it_is_in_both_spellings() {
    // The same `Placed` a panel carries, so a scene that poses a domain poses its parts with it.
    // glTF writes `[x, y, z, w]` and USD writes `(w, x, y, z)`; the reordering is checked in
    // `a_placement_reaches_the_file.rs` and what matters here is that a *design* goes through it
    // at all rather than being written at the origin.
    let h: f64 = std::f64::consts::PI / 6.0;
    let place = Placed {
        at_m: [0.5, 0.0, 0.0],
        turn: [0.0, 0.0, h.sin(), h.cos()],
    };
    let drawing = Drawing::of(Surfaces::Boundary).and(vec![designed(place, 0.04)]);

    let gl = gltf_with("a run", &cells(), &drawing).document;
    // The numbers, not their spelling: `compact` writes `5e-1` for a half, and a test that
    // matched `0.5` would be checking the formatter rather than the placement. It did, and failed
    // on a writer that was doing the right thing.
    let at = gl
        .find("\"translation\":[")
        .unwrap_or_else(|| panic!("the design is at the origin in the glTF: {gl:.400}"));
    let rest = &gl[at + "\"translation\":[".len()..];
    let nums: Vec<f64> = rest[..rest.find(']').expect("a closed array")]
        .split(',')
        .map(|n| n.parse().expect("a number"))
        .collect();
    assert_eq!(nums, vec![0.5, 0.0, 0.0]);

    let us = usda_with("a run", &[cells()], &drawing).document;
    let design = &us[us.find("def Mesh \"part_parts_0_\"").expect("the design")..];
    assert!(
        design.contains("double3 xformOp:translate = (5.000000e-1, 0, 0)"),
        "the design is at the origin in the USD"
    );
}

#[test]
fn the_cells_are_bigger_than_the_design_and_that_is_the_loss() {
    // **What the pair is for.** A shape rasterised onto a grid is drawn on cell boundaries, so it
    // reaches past the surface it came from — that overshoot *is* the rasterisation, and until
    // both were in one file the only way to see it was a printed volume error.
    //
    // The tetrahedron is 40 mm on a side and the field's cells fill a 40 mm box, so the cells are
    // a strictly worse answer than the design by construction: same corner, same reach, but the
    // design is a tetrahedron and the cells are the whole box.
    let drawing = Drawing::of(Surfaces::Boundary).and(vec![designed(Placed::HERE, 0.04)]);
    let out = gltf_with("a run", &cells(), &drawing);
    assert!(
        out.notes
            .iter()
            .any(|n| n.contains("1 designed part") && n.contains("rasterisation loss")),
        "the file does not say what the pair is for: {:?}",
        out.notes
    );

    // The design's own volume is a sixth of the box it sits in; the cells are all of it. A reader
    // opening the file sees both and the difference is five sixths of the box.
    let s: f64 = 0.04;
    let surface = mesh::mesh_surface(&tetra(s), 0);
    let mut signed = 0.0;
    for f in 0..surface.faces() {
        let p = |k: usize| {
            let i = surface.indices[3 * f + k] as usize;
            [
                f64::from(surface.positions[i][0]),
                f64::from(surface.positions[i][1]),
                f64::from(surface.positions[i][2]),
            ]
        };
        let (a, b, c) = (p(0), p(1), p(2));
        signed += a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]);
    }
    let volume = (signed / 6.0).abs();
    let exact = s * s * s / 6.0;
    assert!(
        (volume - exact).abs() / exact < 1e-6,
        "the designed tetrahedron is {volume} m³ against {exact}"
    );
}

#[test]
fn an_outline_is_the_meshs_edges_and_eulers_formula_says_how_many() {
    // **The closed form for a wireframe.** For any closed polyhedron `V - E + F = 2`, so a
    // tetrahedron's four faces and four corners force exactly six edges. Nothing about the
    // implementation is consulted: if it returned three per face it would say twelve, and if it
    // deduplicated by index rather than by position it would say twelve too — `mesh_surface` gives
    // every face its own vertices so each can carry a flat normal, and that is the trap.
    let d = designed(Placed::HERE, 0.04);
    let pts = d.outline_points();
    let edges = d.outline().expect("a tetrahedron is under any budget");
    assert_eq!(d.surface.faces(), 4, "F");
    assert_eq!(
        pts.len(),
        4,
        "V — four corners, from twelve stored vertices"
    );
    assert_eq!(edges.len(), 6, "E — Euler: 4 - E + 4 = 2");
    assert_eq!(
        pts.len() as i64 - edges.len() as i64 + d.surface.faces() as i64,
        2
    );

    // Every pair indexes a point that exists, and no edge joins a point to itself.
    for [a, b] in &edges {
        assert!((*a as usize) < pts.len() && (*b as usize) < pts.len());
        assert_ne!(a, b);
    }
}

#[test]
fn the_outline_comes_out_placed() {
    // A wireframe is a vertex list, and a rigid motion multiplied into one loses nothing — unlike a
    // box, which is why a *field* states its placement instead. So these are in the world, and the
    // report draws them without a transform of its own.
    let here = designed(Placed::HERE, 0.04).outline_points();
    let there = designed(
        Placed {
            at_m: [1.0, 2.0, 3.0],
            ..Placed::HERE
        },
        0.04,
    )
    .outline_points();
    assert_eq!(here.len(), there.len());
    for (a, b) in here.iter().zip(&there) {
        assert_eq!([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [1.0, 2.0, 3.0]);
    }
}

#[test]
fn an_outline_past_the_budget_is_refused_rather_than_thinned() {
    // Every second edge of a wireframe is not a coarser picture of the part, it is a picture of a
    // different one — which is why this refuses where `MAX_FACES` subsamples. A field's boundary at
    // a stride is still that boundary at a stride; half an outline is nothing.
    //
    // One triangle contributes three edges, so the budget is crossed a little over a third of the
    // way through that many triangles. Built at the smallest size that does it rather than at a
    // round number, because the assertion is about the threshold.
    let n = mesh::MAX_OUTLINE_EDGES / 3 + 1;
    let far = |i: usize| i as f64 * 1e-3;
    let tris: Vec<[[f64; 3]; 3]> = (0..n)
        .map(|i| {
            [
                [far(i), 0.0, 0.0],
                [far(i) + 1e-4, 0.0, 0.0],
                [far(i), 1e-4, 0.0],
            ]
        })
        .collect();
    let big = Designed {
        name: "huge".into(),
        domain: "huge".into(),
        place: Placed::HERE,
        surface: mesh::mesh_surface(&tris, 0),
    };
    assert!(big.surface.faces() > mesh::MAX_OUTLINE_EDGES / 3);
    assert!(
        big.outline().is_none(),
        "{} faces should be past the budget",
        big.surface.faces()
    );

    // And just under it comes back whole.
    let small = Designed {
        surface: mesh::mesh_surface(&tris[..n - 2], 0),
        ..big
    };
    assert!(small.outline().is_some());
}

#[test]
fn the_report_carries_the_outline_and_says_which_domain_it_is_for() {
    // **Where the expectation lives.** `tools/report-check` runs the report's JavaScript and
    // asserts it drew one line per edge — but only when the run *has* a design, so dropping the
    // outline from the JSON made it assert nothing at all and pass. Measured: that sabotage
    // survived until this test existed.
    //
    // A JS harness can check "the viewer drew what the data said". It cannot check "the data
    // should have been there", because the page is its only input. That belongs here.
    let drawing = Drawing::of(Surfaces::Boundary).and(vec![designed(Placed::HERE, 0.04)]);
    let page = pantometry_view::html_with("a run", &[cells()], &drawing);

    let at = page
        .find("\"design\":")
        .expect("the run JSON has no design key at all");
    let head = &page[at..at + 200.min(page.len() - at)];
    assert!(
        head.contains("\"part\":"),
        "the outline is not keyed by its domain: {head}"
    );
    assert!(
        head.contains("part/parts[0]"),
        "the outline does not carry the site a finding names: {head}"
    );
    // Four points and six edges, flattened: twelve numbers and twelve.
    assert!(head.contains("\"p\":["), "no points: {head}");
    assert!(head.contains("\"e\":["), "no edges: {head}");
}

#[test]
fn a_report_of_a_run_with_no_design_says_so_rather_than_omitting_the_key() {
    // An empty object and not a missing key. A reader — the viewer, or the next person writing a
    // check — can tell "this run has no design" from "this build does not write designs" only if
    // the key is always there.
    let page = pantometry_view::html_with("a run", &[cells()], &Drawing::default());
    assert!(
        page.contains("\"design\":{}"),
        "the key is missing entirely"
    );
}
