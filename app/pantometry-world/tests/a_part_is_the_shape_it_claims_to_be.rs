//! **Every shipped part, against a closed form for its volume and its surface area.**
//!
//! A part is an STL a scene names, and `tools/parts/make.py` writes the six of them. They are all
//! polygons extruded along `z`, which is the one primitive whose two numbers follow from the
//! profile without trusting the generator:
//!
//! ```text
//! volume = cross_section * height
//! area   = 2 * cross_section + perimeter * height
//! ```
//!
//! The cross-sections are rectangles, triangles, a chamfered L, a comb and two regular polygons,
//! and each of those has an area and a perimeter of its own. Nothing here compares the mesh to a
//! second mesh.
//!
//! The L-bracket predates the generator and was the only part in the tree for a long time, with
//! **no check at all**. It is in the table below like the rest.
//!
//! # Why this reads files rather than building them
//!
//! `an_assembly_from_files` writes its bricks in the test, and says why: geometry under test
//! should be beside the assertion. These are the opposite case — the files are the artefact, they
//! are committed, and what has to be true is a property of *them*. A fixture written here would
//! check the fixture.
#![cfg(not(target_family = "wasm"))]

use pantometry::shape::Mesh;

/// A millimetre in metres.
///
/// An STL is unitless and every CAD tool writes millimetres, so `Mesh::from_stl` scales by `1e-3`
/// and everything it hands back is SI. The closed forms below are in millimetres because that is
/// what the profiles are drawn in.
const MM: f64 = 1e-3;

/// Six significant figures, which is what `%g` writes and what the committed files therefore hold.
///
/// **Derived, and the first derivation of it was wrong.** Six significant figures is a half-ulp of
/// `5e-6` relative when the mantissa is just above 1 and `5e-7` when it is just below 10, so the
/// worst case per coordinate is `5e-6` and not the `5e-7` this said. A volume is a product of
/// three lengths, which takes it to `1.5e-5`; and the **pipe's cross-section is a difference of
/// two nearly equal polygon areas**, `310.58 - 152.19`, which amplifies it by
/// `(A_out + A_in) / (A_out - A_in) = 2.92` to `2.9e-5`.
///
/// `1e-4` is about three times the worst the file format can produce and fifteen times the worst
/// measured, which is the pipe's volume at `6.8e-6`. It is still a hundred times tighter than the
/// smallest defect any sabotage of this library has produced — the rod one per cent long, at
/// `1e-2`. Four of the six parts have integer coordinates and are exact.
///
/// The old value, `1e-5`, sat **below** what the format can do: the parts happened to land inside
/// it, and moving them to the origin moved the pipe from `2.0e-7` to `6.8e-6` for no reason but
/// the rounding of a different coordinate.
const ROUNDING: f64 = 1e-4;

fn part(name: &str) -> Mesh {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("scenes")
        .join("parts")
        .join(format!("{name}.stl"));
    let bytes = std::fs::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    Mesh::from_stl(&bytes).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// The volume and surface area of a prism, from its cross-section and perimeter. In millimetres.
fn prism(cross: f64, perimeter: f64, height: f64) -> (f64, f64) {
    (cross * height, 2.0 * cross + perimeter * height)
}

/// A regular `n`-gon's area and perimeter at circumradius `r`: `(n/2)r² sin(2π/n)`, `2nr sin(π/n)`.
///
/// The faceted rod and pipe are polygons and not circles, so this is the closed form and `πr²` is
/// not — `the_rods_facets_fall_short_of_a_circle_by_what_the_series_says` is where the two meet.
fn regular(n: f64, r: f64) -> (f64, f64) {
    let pi = std::f64::consts::PI;
    (
        0.5 * n * r * r * (2.0 * pi / n).sin(),
        2.0 * n * r * (pi / n).sin(),
    )
}

/// Name, triangles, volume and area in millimetres — each from the profile it is drawn from.
fn shipped() -> Vec<(&'static str, usize, f64, f64)> {
    let mut out: Vec<(&'static str, usize, f64, f64)> = Vec::new();

    // A rectangular slab, 60 by 40 by 5.
    let (v, a) = prism(60.0 * 40.0, 2.0 * (60.0 + 40.0), 5.0);
    out.push(("plate", 12, v, a));

    // A round bar: a 24-gon of circumradius 8, 60 long.
    let (c, p) = regular(24.0, 8.0);
    let (v, a) = prism(c, p, 60.0);
    out.push(("rod", 92, v, a));

    // A right triangular prism on a 40 by 20 triangle, 30 long. The hypotenuse is the third side.
    let (v, a) = prism(
        0.5 * 40.0 * 20.0,
        40.0 + 20.0 + (40.0f64 * 40.0 + 20.0 * 20.0).sqrt(),
        30.0,
    );
    out.push(("wedge", 8, v, a));

    // Two arms at a right angle: a 50 square less a 30 notch, the inner corner chamfered by a
    // triangle of legs 10. Six sides of the seven are axis-aligned and sum to 180; the chamfer is
    // the hypotenuse of that triangle.
    let (v, a) = prism(
        50.0 * 50.0 - 30.0 * 30.0 + 0.5 * 10.0 * 10.0,
        180.0 + 10.0 * 2.0f64.sqrt(),
        20.0,
    );
    out.push(("l-bracket", 24, v, a));

    // A finned base: 60 by 4, carrying five fins 6 wide and 16 tall with 5 between them and at
    // each end. The perimeter is the bottom, the two ends, the six gaps and the five fins.
    let (v, a) = prism(
        60.0 * 4.0 + 5.0 * 6.0 * 16.0,
        60.0 + 2.0 * 4.0 + 6.0 * 5.0 + 5.0 * (2.0 * 16.0 + 6.0),
        30.0,
    );
    out.push(("heat-sink", 92, v, a));

    // A hollow tube: a 24-gon of 10 with a 24-gon hole of 7, 50 long. The hole subtracts from the
    // cross-section and *adds* to the perimeter, which is the whole difference between a tube and
    // a bar.
    let (co, po) = regular(24.0, 10.0);
    let (ci, pi) = regular(24.0, 7.0);
    let (v, a) = prism(co - ci, po + pi, 50.0);
    out.push(("pipe", 192, v, a));

    out
}

/// **Each part's volume and area are the ones its profile says.**
#[test]
fn every_part_is_the_shape_it_claims_to_be() {
    for (name, triangles, want_v, want_a) in shipped() {
        let mesh = part(name);
        let v = mesh.volume().to_si() / (MM * MM * MM);
        let a = mesh.area().to_si() / (MM * MM);
        println!(
            "  {name:<10} {:>3} triangles   V {v:>10.4} mm3 (closed form {want_v:>10.4})   A {a:>10.4} mm2 ({want_a:>10.4})",
            mesh.triangles().len()
        );

        // Closed first: for an open mesh the volume is meaningless rather than approximate, so a
        // volume that agrees would be agreeing by accident.
        assert!(
            mesh.is_closed(),
            "{name} has an edge not shared by exactly two triangles, so it encloses nothing"
        );
        // Positive, and not merely close in magnitude. An inside-out export is a real defect and
        // the only check that sees it is the sign.
        assert!(
            v > 0.0,
            "{name}'s volume is {v:.4} mm3, so its winding is inside out"
        );
        assert!(
            (v - want_v).abs() < ROUNDING * want_v,
            "{name}'s volume is {v:.6} mm3 and its profile says {want_v:.6}"
        );
        assert!(
            (a - want_a).abs() < ROUNDING * want_a,
            "{name}'s area is {a:.6} mm2 and its profile says {want_a:.6}"
        );
        // **Its lowest corner is the origin**, which is the check that was missing and the one
        // that matters for using a part at all. `Voxels::onto` reads an STL's coordinates as
        // absolute positions and a `block` domain's grid starts at the origin, so a solid drawn
        // around the origin reaches outside its grid and the build refuses it. The rod and the
        // pipe shipped that way: volume, area and the closed-edge count are all
        // translation-invariant, so none of the three assertions above could see it, and it took
        // `a_dropped_file_becomes_a_domain` building a scene to find it.
        let (low, _) = mesh.bounds().expect("a closed part has bounds");
        let low = low.to_si();
        for (axis, at) in [("x", low.x), ("y", low.y), ("z", low.z)] {
            assert!(
                at.abs() < 1e-12,
                "{name} starts at {at:e} m on {axis} rather than at the origin, so a block domain                  would refuse it for reaching outside its own grid"
            );
        }
        // A pin and not a closed form: how many triangles a profile needs is the generator's
        // choice, and this is here so that changing it is deliberate.
        assert_eq!(
            mesh.triangles().len(),
            triangles,
            "{name} has {} triangles and was written with {triangles}",
            mesh.triangles().len()
        );
    }
}

/// Write an ASCII STL of a polygon extruded along `z`, with the caps **fanned** from vertex zero.
///
/// The normals are written as zero because `Mesh::from_stl` ignores them and takes the winding —
/// exporters disagree about normals often enough that the vertices are the more reliable of the
/// two, which that type's documentation says and this relies on.
fn fanned(profile: &[[f64; 2]], height: f64) -> String {
    fn facet(s: &mut String, a: [f64; 3], b: [f64; 3], c: [f64; 3]) {
        s.push_str("  facet normal 0 0 0\n    outer loop\n");
        for p in [a, b, c] {
            s.push_str(&format!("      vertex {} {} {}\n", p[0], p[1], p[2]));
        }
        s.push_str("    endloop\n  endfacet\n");
    }
    let n = profile.len();
    let mut s = String::from("solid fanned\n");
    for i in 0..n {
        let (ax, ay) = (profile[i][0], profile[i][1]);
        let (bx, by) = (profile[(i + 1) % n][0], profile[(i + 1) % n][1]);
        facet(&mut s, [ax, ay, 0.0], [bx, by, 0.0], [bx, by, height]);
        facet(&mut s, [ax, ay, 0.0], [bx, by, height], [ax, ay, height]);
    }
    for i in 1..n - 1 {
        let (p0, p1, p2) = (profile[0], profile[i], profile[i + 1]);
        facet(
            &mut s,
            [p0[0], p0[1], 0.0],
            [p2[0], p2[1], 0.0],
            [p1[0], p1[1], 0.0],
        );
        facet(
            &mut s,
            [p0[0], p0[1], height],
            [p1[0], p1[1], height],
            [p2[0], p2[1], height],
        );
    }
    s.push_str("endsolid fanned\n");
    s
}

/// **A fan gets the volume exactly right and the area wrong**, which is why the generator clips
/// ears instead.
///
/// A fan from vertex zero is only a triangulation when every vertex is visible from vertex zero.
/// The U below is not: the straight line from its corner to the far top corner leaves the solid
/// through the slot. Fanning it puts triangles in the air across that slot.
///
/// **Neither of the other two checks can tell.** The signed areas of a fan telescope to the
/// shoelace sum whatever the polygon looks like, so the volume is exact; and a fan is still a
/// valid combinatorial triangulation, so every edge is shared by exactly two triangles and the
/// mesh is closed. Only the surface area sees it, because triangles that overlap or leave the
/// polygon do not add up to it.
///
/// The heat sink is the shipped part this applies to — no vertex of a comb sees all the others.
#[test]
fn a_fan_passes_the_volume_and_the_closed_check_and_fails_the_area_one() {
    // A U: a 30 by 20 block with a slot 10 wide cut 15 deep into its top.
    let u = [
        [0.0, 0.0],
        [30.0, 0.0],
        [30.0, 20.0],
        [20.0, 20.0],
        [20.0, 5.0],
        [10.0, 5.0],
        [10.0, 20.0],
        [0.0, 20.0],
    ];
    let height = 10.0;
    let cross = 30.0 * 20.0 - 10.0 * 15.0;
    let perimeter = 30.0 + 20.0 + 10.0 + 15.0 + 10.0 + 15.0 + 10.0 + 20.0;
    let (want_v, want_a) = prism(cross, perimeter, height);

    let mesh = Mesh::from_stl(fanned(&u, height).as_bytes()).expect("the fanned U parses");
    let v = mesh.volume().to_si() / (MM * MM * MM);
    let a = mesh.area().to_si() / (MM * MM);
    println!("  fanned U   V {v:.4} mm3 (closed form {want_v:.4})   A {a:.4} mm2 ({want_a:.4})");

    assert!(
        mesh.is_closed(),
        "a fan is a valid combinatorial triangulation, and the point here is that being closed \
         does not make it a correct one"
    );
    assert!(
        (v - want_v).abs() < 1e-9 * want_v,
        "a fan's signed areas telescope to the shoelace sum, so its volume should be exact: \
         {v:.9} against {want_v:.9}"
    );
    // The slot is 10 by 15, and the fan lays triangles across it twice over — once per cap. The
    // assertion is that the area is wrong, which is the finding; how wrong is printed.
    assert!(
        a > want_a * 1.05,
        "a fan across the slot should overstate the area, and this one gives {a:.4} mm2 against \
         {want_a:.4}. If they now agree, either the polygon became star-shaped or `Mesh::area` \
         stopped summing what it sums"
    );

    // And the shipped comb is the same shape of polygon, so the generator must not fan it.
    let sink = part("heat-sink");
    let (_, _, want_v, want_a) = shipped()
        .into_iter()
        .find(|p| p.0 == "heat-sink")
        .expect("the heat sink is shipped");
    let sink_a = sink.area().to_si() / (MM * MM);
    println!(
        "  heat-sink area {sink_a:.4} mm2 against {want_a:.4}, volume closed form {want_v:.1}"
    );
    assert!(
        (sink_a - want_a).abs() < ROUNDING * want_a,
        "the heat sink's area is {sink_a:.4} mm2 and its profile says {want_a:.4}; a fanned comb \
         would read high here and nowhere else"
    );
}

/// **The rod's facets fall short of a circle by what the series says.**
///
/// `(n/2)r² sin(2π/n)` is the polygon's area and `πr²` is the circle's, and the first approaches
/// the second. Checking the shortfall against the expansion of `sin x / x` is the polygon formula
/// held against something that is not the polygon formula — the closed form is checked, rather
/// than being the thing that checks.
#[test]
fn the_rods_facets_fall_short_of_a_circle_by_what_the_series_says() {
    let n = 24.0;
    let r = 8.0;
    let (polygon, _) = regular(n, r);
    let circle = std::f64::consts::PI * r * r;
    let short = 1.0 - polygon / circle;

    // sin x / x = 1 - x²/6 + x⁴/120 - x⁶/5040, with x = 2π/n.
    let x = 2.0 * std::f64::consts::PI / n;
    let predicted = x * x / 6.0 - x.powi(4) / 120.0;
    let next_term = x.powi(6) / 5040.0;
    println!(
        "  a {n}-gon is {:.4}% under its circle; two terms of the series say {:.4}%, and the third is {:.2e}",
        100.0 * short,
        100.0 * predicted,
        next_term
    );

    assert!(
        (short - predicted).abs() < next_term,
        "the shortfall is {short:.9} and two terms of the series give {predicted:.9}; they should \
         differ by less than the term that was dropped, {next_term:.3e}"
    );
    // And it is a per-cent effect and not a rounding one, so the closed forms above had to use the
    // polygon and not the circle.
    assert!(
        short > 0.01,
        "a 24-gon is about 1.14% under its circle, and this says {:.4}% — if that has become \
         small, the rod is not the polygon these closed forms are written for",
        100.0 * short
    );
}
