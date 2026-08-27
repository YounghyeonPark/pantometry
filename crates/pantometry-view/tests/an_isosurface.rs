//! The surface where a field reaches a value, against shapes whose area and volume are known.
//!
//! `field_surface` draws the outside of the cells that hold a value, which is a picture of the
//! grid as much as of the object. An isosurface draws where the *values* reach a number, and the
//! difference is that it has a right answer: sample a sphere and the surface should have the area
//! of a sphere.
//!
//! Every check here is against a closed form rather than against a second implementation, and the
//! two that matter are **rates** rather than sizes — a piecewise-linear surface is never exactly
//! a curved one, so the question is whether refining the grid closes the gap at the right speed.

use pantometry_view::mesh::isosurface;

/// Sample `f` on an `n³` grid spanning `extent` and return the values in the layout the mesher
/// expects.
fn sample(n: usize, extent: [f64; 6], f: impl Fn(f64, f64, f64) -> f64) -> Vec<f64> {
    let step = |lo: f64, hi: f64| (hi - lo) / (n - 1) as f64;
    let (dx, dy, dz) = (
        step(extent[0], extent[3]),
        step(extent[1], extent[4]),
        step(extent[2], extent[5]),
    );
    let mut out = Vec::with_capacity(n * n * n);
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                out.push(f(
                    extent[0] + dx * i as f64,
                    extent[1] + dy * j as f64,
                    extent[2] + dz * k as f64,
                ));
            }
        }
    }
    out
}

/// Total area of the triangles.
fn area(s: &pantometry_view::mesh::Surface) -> f64 {
    s.indices
        .chunks_exact(3)
        .map(|f| {
            let p = |n: u32| s.positions[n as usize];
            let (a, b, c) = (p(f[0]), p(f[1]), p(f[2]));
            let u = [
                (b[0] - a[0]) as f64,
                (b[1] - a[1]) as f64,
                (b[2] - a[2]) as f64,
            ];
            let v = [
                (c[0] - a[0]) as f64,
                (c[1] - a[1]) as f64,
                (c[2] - a[2]) as f64,
            ];
            let n = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
        })
        .sum()
}

/// The volume the mesh encloses, by the divergence theorem: a sixth of the sum of the signed
/// triple products. **Only meaningful if the winding is consistent** — a mesh with faces pointing
/// both ways sums to something near zero, which is why this is also the winding test.
fn volume(s: &pantometry_view::mesh::Surface) -> f64 {
    s.indices
        .chunks_exact(3)
        .map(|f| {
            let p = |n: u32| {
                let q = s.positions[n as usize];
                [q[0] as f64, q[1] as f64, q[2] as f64]
            };
            let (a, b, c) = (p(f[0]), p(f[1]), p(f[2]));
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum()
}

/// **A plane comes out flat, exactly.**
///
/// The one case with no discretisation error at all: the field is linear, so the linear
/// interpolation along each edge is not an approximation. Every vertex must sit on `x = 0.25` to
/// floating point, and a surface that is nearly flat would mean the interpolation is wrong rather
/// than coarse.
#[test]
fn a_linear_field_gives_an_exactly_flat_surface() {
    let extent = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let values = sample(9, extent, |x, _, _| x);
    let s = isosurface((9, 9, 9), extent, &values, 0.25);

    assert!(s.faces() > 0, "the level is inside the range");
    for p in &s.positions {
        assert!(
            (p[0] as f64 - 0.25).abs() < 1e-6,
            "a vertex sits at x = {}, and the plane is at 0.25",
            p[0]
        );
    }
    // And it is the whole cross-section: a unit square.
    assert!(
        (area(&s) - 1.0).abs() < 1e-6,
        "the plane through a unit cube has unit area, got {}",
        area(&s)
    );
}

/// **A sphere's area and volume converge on the closed forms, at second order.**
///
/// `4πR²` and `4/3πR³`. A piecewise-linear surface never reaches either exactly, so the assertion
/// is on the **rate**: halving the spacing should quarter the error. Measured across three
/// refinements rather than asserted at one size, because a first-order error is small at one size
/// too — and first order here would mean the edge interpolation is not doing what it claims.
#[test]
fn a_sphere_converges_on_its_area_and_volume() {
    let r = 0.7;
    let extent = [-1.0, -1.0, -1.0, 1.0, 1.0, 1.0];
    let exact_area = 4.0 * std::f64::consts::PI * r * r;
    let exact_volume = 4.0 / 3.0 * std::f64::consts::PI * r * r * r;

    let mut errors: Vec<(f64, f64)> = Vec::new();
    for n in [21, 41, 81] {
        let values = sample(n, extent, |x, y, z| (x * x + y * y + z * z).sqrt());
        let s = isosurface((n, n, n), extent, &values, r);
        errors.push((
            (area(&s) - exact_area).abs() / exact_area,
            (volume(&s) - exact_volume).abs() / exact_volume,
        ));
    }

    for (label, pick) in [("area", 0usize), ("volume", 1usize)] {
        let e: Vec<f64> = errors
            .iter()
            .map(|p| if pick == 0 { p.0 } else { p.1 })
            .collect();
        assert!(
            e[2] < 0.01,
            "{label}: the finest grid should be within a percent, was {:.3e}",
            e[2]
        );
        for step in 0..2 {
            let order = (e[step] / e[step + 1]).log2();
            assert!(
                (1.5..=2.6).contains(&order),
                "{label}: halving the spacing should quarter the error; \
                 step {step} measured order {order:.2} ({:.3e} -> {:.3e})",
                e[step],
                e[step + 1]
            );
        }
    }
}

/// **The volume is positive, which is the winding test.**
///
/// The divergence-theorem sum is signed. A mesh whose faces point both ways cancels toward zero,
/// and one wound entirely inside-out comes out negative — so a volume that is both positive and
/// the right size says every triangle agrees which way is out. A renderer with back-face culling
/// depends on that, and so does anything that lights the surface.
#[test]
fn every_face_agrees_which_way_is_out() {
    let extent = [-1.0, -1.0, -1.0, 1.0, 1.0, 1.0];
    let values = sample(31, extent, |x, y, z| (x * x + y * y + z * z).sqrt());
    let s = isosurface((31, 31, 31), extent, &values, 0.6);
    let v = volume(&s);
    let exact = 4.0 / 3.0 * std::f64::consts::PI * 0.6f64.powi(3);
    assert!(
        v > 0.0,
        "an inside-out or mixed mesh gives a volume of {v:.3e}"
    );
    assert!(
        (v - exact).abs() / exact < 0.02,
        "{v:.4} against {exact:.4}"
    );
}

/// **A void takes its whole cell out, and does not become a value.**
///
/// A `NaN` sample is what `Solid3D` returns for an emptied cell, deliberately. Interpolating an
/// edge that ends in one would invent the number the `NaN` exists to refuse, and the invented
/// number would be plotted as a surface somebody believes.
#[test]
fn a_void_is_not_interpolated_through() {
    let extent = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let n = 11;
    let whole = sample(n, extent, |x, _, _| x);
    let full = isosurface((n, n, n), extent, &whole, 0.5);

    // Empty a slab on one side of the surface.
    let mut holed = whole.clone();
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                if j >= n / 2 {
                    holed[i + n * (j + n * k)] = f64::NAN;
                }
            }
        }
    }
    let cut = isosurface((n, n, n), extent, &holed, 0.5);

    assert!(cut.faces() > 0, "half the field is still there");
    assert!(
        cut.faces() < full.faces(),
        "the voided half should be missing: {} against {}",
        cut.faces(),
        full.faces()
    );
    assert!(
        cut.positions
            .iter()
            .all(|p| p.iter().all(|c| c.is_finite())),
        "a vertex was interpolated into the void"
    );
    // Nothing may sit inside the voided region.
    let highest = cut.positions.iter().map(|p| p[1]).fold(f32::MIN, f32::max);
    assert!(
        (highest as f64) <= 0.5 + 1e-6,
        "a vertex at y = {highest} is inside the void"
    );
}

/// **A level outside the field's range is an empty surface, not a guess.**
#[test]
fn a_level_nothing_reaches_draws_nothing() {
    let extent = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let values = sample(9, extent, |x, _, _| x);
    for level in [-1.0, 2.0, f64::NAN] {
        let s = isosurface((9, 9, 9), extent, &values, level);
        assert_eq!(s.faces(), 0, "level {level} is not in the field");
    }
}
