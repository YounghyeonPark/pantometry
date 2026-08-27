//! **A field's surface fills the box it was placed in, and faces outwards.**
//!
//! Both claims have a history. The size one has been wrong three times in this workspace — a tube,
//! a room, and a glTF export that made a 40 mm cube 80 mm across — always the same arithmetic: a
//! field is sampled corner to corner, so an end node owns half a cell and a surface drawn on cell
//! *centres* is one whole cell too big. The outward one is new here, because the editor's boxes can
//! be **placed**: a rotated box's face normals are the cross products of its own edge vectors, and
//! getting a sign wrong there lights every face from inside and the solid reads as a hole.

use editor_core::{corners_of, field_shell};

/// A cube of one value, all cells present.
fn solid(n: usize) -> Vec<f64> {
    vec![20.0; n * n * n]
}

fn bounds(positions: &[[f64; 3]]) -> [f64; 6] {
    let mut b = [f64::MAX, f64::MAX, f64::MAX, f64::MIN, f64::MIN, f64::MIN];
    for p in positions {
        for a in 0..3 {
            b[a] = b[a].min(p[a]);
            b[a + 3] = b[a + 3].max(p[a]);
        }
    }
    b
}

/// The surface's box is the placed box, to the last significant digit.
#[test]
fn it_is_exactly_as_big_as_the_box_it_was_placed_in() {
    // 40 mm across, offset from the origin, so a bug that dropped the offset and a bug that
    // doubled the size are different failures.
    let extent = [0.1, -0.02, 0.5, 0.14, 0.02, 0.54];
    let corners = corners_of(extent);

    for n in [2usize, 3, 9, 20] {
        let out = field_shell(
            &corners,
            (n, n, n),
            &solid(n),
            "C",
            Some((0.0, 100.0)),
            None,
        );
        assert!(
            !out.indices.is_empty(),
            "a solid {n}^3 has a surface: {}",
            out.note
        );
        let got = bounds(&out.positions);
        for a in 0..6 {
            let e = (got[a] - extent[a]).abs();
            // The mesh carries positions as `f32`, so the budget is that type's, on numbers of
            // order 0.5 m: about `6e-8` relative, so `1e-7` m absolute with room. Not a physics
            // tolerance — the width of the type the geometry is stored in.
            assert!(
                e < 1e-7,
                "{n}^3 axis {a}: surface at {} against a box at {} — off by {e:.3e} m",
                got[a],
                extent[a]
            );
        }
    }
}

/// Every normal points away from the centre of a solid convex block.
///
/// Which is the whole of what "outward" means for one, and it is checked on a **rotated** box
/// because that is where the cross products earn their keep.
#[test]
fn the_faces_look_outwards_even_when_the_box_is_turned() {
    let n = 6;
    let values = solid(n);

    // Corner 0 is the low one and bits 0, 1, 2 step one axis each. Turned 40° about z and tilted,
    // with the three edges kept perpendicular so "outward" is unambiguous.
    let (c, s) = (0.4f64.cos(), 0.4f64.sin());
    let o = [0.3, 0.15, -0.4];
    let e = [
        [0.05 * c, 0.05 * s, 0.0],
        [-0.03 * s, 0.03 * c, 0.0],
        [0.0, 0.0, 0.08],
    ];
    let corners: [[f64; 3]; 8] = std::array::from_fn(|i| {
        [0, 1, 2].map(|k| {
            let mut v = o[k];
            for (bit, edge) in e.iter().enumerate() {
                if i & (1 << bit) != 0 {
                    v += edge[k];
                }
            }
            v
        })
    });

    let out = field_shell(&corners, (n, n, n), &values, "C", Some((0.0, 100.0)), None);
    assert!(
        !out.indices.is_empty(),
        "a turned solid still has a surface"
    );

    let centre = [
        o[0] + 0.5 * (e[0][0] + e[1][0] + e[2][0]),
        o[1] + 0.5 * (e[0][1] + e[1][1] + e[2][1]),
        o[2] + 0.5 * (e[0][2] + e[1][2] + e[2][2]),
    ];

    let mut worst = f64::MAX;
    let mut inward = 0usize;
    for (p, nrm) in out.positions.iter().zip(&out.normals) {
        let away = [p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]];
        let dot = away
            .iter()
            .zip(nrm)
            .map(|(a, b)| a * *b as f64)
            .sum::<f64>();
        // Normalise by the box's own size so the number is a direction and not a length.
        let reach = (away[0] * away[0] + away[1] * away[1] + away[2] * away[2]).sqrt();
        let cosine = dot / reach.max(f64::MIN_POSITIVE);
        if cosine <= 0.0 {
            inward += 1;
        }
        worst = worst.min(cosine);
    }

    println!(
        "  {} vertices, worst outward cosine {worst:.4}",
        out.normals.len()
    );
    assert_eq!(
        inward,
        0,
        "{inward} of {} vertices face into the block",
        out.normals.len()
    );
    // **The worst case has a closed form, so it is computed rather than guessed.** A corner vertex
    // sits at `(±a, ±b, ±c)` from the centre and carries one of the three face normals, so its
    // cosine is that half-extent over the half-diagonal — and the worst over the whole surface is
    // the *smallest* half-extent over it. A threshold picked by eye was `0.3` against a measured
    // `0.3030`, which is a tolerance that passes on this box and fails on the next one.
    let half = [
        0.5 * (e[0][0] * e[0][0] + e[0][1] * e[0][1] + e[0][2] * e[0][2]).sqrt(),
        0.5 * (e[1][0] * e[1][0] + e[1][1] * e[1][1] + e[1][2] * e[1][2]).sqrt(),
        0.5 * (e[2][0] * e[2][0] + e[2][1] * e[2][1] + e[2][2] * e[2][2]).sqrt(),
    ];
    let diagonal = (half[0] * half[0] + half[1] * half[1] + half[2] * half[2]).sqrt();
    let expected = half.iter().copied().fold(f64::MAX, f64::min) / diagonal;
    println!("  the corner cosine this box predicts is {expected:.4}");
    assert!(
        (worst - expected).abs() < 1e-3,
        "worst outward cosine {worst:.4} against {expected:.4} predicted"
    );
}

/// An absent cell in the middle adds the six faces its neighbours grew.
#[test]
fn a_void_inside_becomes_a_cavity_with_a_surface() {
    let n = 5;
    let corners = corners_of([0.0, 0.0, 0.0, 0.05, 0.05, 0.05]);
    let scale = Some((0.0, 100.0));

    let whole = field_shell(&corners, (n, n, n), &solid(n), "C", scale, None);

    // A void is *not a value*, which is how `Solid3D` reports one: a void has no temperature and a
    // zero there is a number somebody would plot.
    let mut with_hole = solid(n);
    with_hole[2 + n * (2 + n * 2)] = f64::NAN;
    let holed = field_shell(&corners, (n, n, n), &with_hole, "C", scale, None);

    let (before, after) = (whole.indices.len() / 3, holed.indices.len() / 3);
    println!("  {before} triangles solid, {after} with one cell absent");
    assert_eq!(
        after,
        before + 12,
        "an interior cell going absent uncovers its six neighbours' faces — twelve triangles"
    );
}

/// A field one cell thick in two directions is a graph, and it says so rather than drawing dots.
#[test]
fn a_line_of_samples_is_not_geometry() {
    let corners = corners_of([0.0, 0.0, 0.0, 0.1, 0.001, 0.001]);
    let out = field_shell(
        &corners,
        (32, 1, 1),
        &vec![20.0; 32],
        "C",
        Some((0.0, 100.0)),
        None,
    );
    assert!(out.indices.is_empty(), "a line has no surface");
    assert!(
        out.note.contains("graph"),
        "and the canvas is told why: {}",
        out.note
    );
}

/// The surface's colours follow its values, and follow the same ramp everything else uses.
///
/// Worth its own test because a shaded viewport can look completely convincing while painting one
/// colour: geometry, lighting and depth are all visible in a screenshot and a constant colour is
/// not. A block diffusing a hot spot really *is* almost uniform after a while — the run-wide scale
/// says so — so the picture cannot tell you whether the colour path works.
#[test]
fn the_colours_follow_the_values() {
    let n = 6;
    let corners = corners_of([0.0, 0.0, 0.0, 0.06, 0.06, 0.06]);

    // A ramp along x, so opposite faces of the block are at the two ends of the scale.
    let mut values = vec![0.0; n * n * n];
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                values[i + n * (j + n * k)] = 100.0 * i as f64 / (n - 1) as f64;
            }
        }
    }
    let out = field_shell(&corners, (n, n, n), &values, "C", Some((0.0, 100.0)), None);
    assert!(!out.colours.is_empty(), "there are colours: {}", out.note);

    // The two ends of the scale, as the ramp gives them.
    let cold = out
        .source
        .iter()
        .position(|c| values[*c as usize] == 0.0)
        .expect("a cell at the bottom of the scale is on the surface");
    let hot = out
        .source
        .iter()
        .position(|c| values[*c as usize] == 100.0)
        .expect("a cell at the top of the scale is on the surface");

    let (a, b) = (out.colours[cold], out.colours[hot]);
    let apart = (0..3).map(|i| (a[i] - b[i]).abs()).fold(0.0f32, f32::max);
    println!("  cold {a:?}  hot {b:?}  furthest channel apart by {apart:.3}");
    // The library's scale runs 17.9 to 92.0 in L*, so its ends are nowhere near each other in any
    // colour space. A tenth of the unit cube is a floor far below that and far above the `1e-7` a
    // constant-colour bug would give.
    assert!(
        apart > 0.1,
        "the two ends of the scale are the same colour: {a:?} against {b:?}"
    );

    // And it is the same ramp the rest of the workspace draws with, converted to linear. Checked
    // against `value_colour` rather than against a second copy of the formula.
    let srgb = editor_core::value_colour(100.0, Some((0.0, 100.0)));
    let expect: [f32; 3] = std::array::from_fn(|i| {
        let c = srgb[i] as f32 / 255.0;
        if c <= 0.04045 {
            c / 12.92
        } else {
            ((c + 0.055) / 1.055).powf(2.4)
        }
    });
    for i in 0..3 {
        assert!(
            (b[i] - expect[i]).abs() < 1e-6,
            "channel {i}: surface {} against the ramp's {}",
            b[i],
            expect[i]
        );
    }
}
