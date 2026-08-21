//! **The two things both shells were doing themselves, and one of them was doing wrong.**
//!
//! A viewport that composites translucent splats has to paint them far to near, and a viewport
//! that colours a value has to place it on a scale. Both were written twice — once in the native
//! shell and once in the browser's — which is precisely the shape `editor-core` exists to prevent,
//! and precisely the shape that let one of the two be wrong for as long as it was.
//!
//! The camera here is `viewer-core`'s real one, not a stand-in, because the claim is about what
//! the projection returns and a fake projection cannot fail the way the real one did.

use viewer_core::{Camera, Framing};

/// A 6x6x6 lattice of splat centres in a unit box, which is the thing that actually gets sorted.
fn lattice(n: usize) -> Vec<[f64; 3]> {
    let f = |i: usize| i as f64 / (n - 1) as f64 - 0.5;
    (0..n)
        .flat_map(|k| (0..n).flat_map(move |j| (0..n).map(move |i| [f(i), f(j), f(k)])))
        .collect()
}

/// The stand-in the native shell used: project the point again a millimetre along **world z** and
/// take the reciprocal of how far the two land apart on screen.
///
/// Reproduced here so the test measures the thing that was wrong rather than describing it.
fn old_proxy(cam: &Camera, frame: &Framing, aspect: f64, p: [f64; 3]) -> f64 {
    let a = cam.project(p, frame, aspect);
    let b = cam.project([p[0], p[1], p[2] + 1e-3], frame, aspect);
    let sep = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2)).sqrt();
    if sep > 0.0 {
        1.0 / sep
    } else {
        f64::MAX
    }
}

/// How many pairs one ordering puts the other way round from the true depth.
fn disagreements(truth: &[f64], claim: &[f64]) -> (usize, usize) {
    let n = truth.len();
    let mut wrong = 0;
    for i in 0..n {
        for j in i + 1..n {
            if (truth[i] - truth[j]) * (claim[i] - claim[j]) < 0.0 {
                wrong += 1;
            }
        }
    }
    (wrong, n * (n - 1) / 2)
}

/// **`far_to_near` orders by the depth the camera computed, exactly.**
///
/// And the proxy it replaces does not: on this lattice it disagrees with the truth on 26% of
/// pairs. The bound below is 1% rather than 0 only for the ties a lattice produces — points at
/// genuinely equal depth, which no ordering can get wrong — and the measurement is 0.
#[test]
fn the_splat_order_is_the_camera_s_depth_and_the_old_proxy_was_not() {
    let pts = lattice(6);
    let frame = Framing::of([-0.5, -0.5, -0.5, 0.5, 0.5, 0.5]);
    // **The camera the app opens at**, untouched. An earlier draft called `turn(0.7, 0.4)` on top
    // of it, and `turn` *adds* — so it measured a camera at 1.4 and 0.8 radians and reported 4.5%
    // where the default reports 26%. A number about a view nobody looks at is not a measurement,
    // and 26% is the figure `far_to_near`'s own documentation quotes.
    let cam = Camera::default();
    let aspect = 1.6;

    let truth: Vec<f64> = pts
        .iter()
        .map(|p| cam.project(*p, &frame, aspect).depth)
        .collect();
    let order = editor_core::far_to_near(&truth);

    // Furthest first, and no pair out of place.
    for w in order.windows(2) {
        assert!(
            truth[w[0]] >= truth[w[1]],
            "{} at depth {} was painted before {} at depth {}",
            w[0],
            truth[w[0]],
            w[1],
            truth[w[1]]
        );
    }

    // And the proxy, measured on the same points and the same camera.
    let proxy: Vec<f64> = pts
        .iter()
        .map(|p| old_proxy(&cam, &frame, aspect, *p))
        .collect();
    let (wrong, total) = disagreements(&truth, &proxy);
    // 6,026 of 23,220 — 25.9% — at `Camera::default()`. A *lower* bound rather than the exact
    // count, so a change to the default view does not fail this: the point is that the proxy was
    // badly wrong, and if it ever measures as correct this test has stopped measuring anything.
    assert!(
        wrong * 100 / total >= 20,
        "the proxy disagreed on only {wrong} of {total} pairs; if it has become correct this \
         test has stopped measuring anything"
    );

    // The failure is worst where the object is. A point near the view axis has almost no screen
    // motion to measure, so its recovered depth is enormous and it paints first whatever its
    // real depth: two points on the axis at different depths are ordered by the proxy's noise.
    let near_axis = [0.0, 0.0, -0.45];
    let far_axis = [0.0, 0.0, 0.45];
    let (dn, df) = (
        cam.project(near_axis, &frame, aspect).depth,
        cam.project(far_axis, &frame, aspect).depth,
    );
    let (pn, pf) = (
        old_proxy(&cam, &frame, aspect, near_axis),
        old_proxy(&cam, &frame, aspect, far_axis),
    );
    assert!(dn < df, "the second point is further: {dn} then {df}");
    assert!(
        pn > 1e3 && pf > 1e3,
        "both axis points should have collapsed under the proxy: {pn}, {pf}"
    );
}

/// A `NaN` depth does not take the paint loop down with it.
///
/// `sort_by` over a partial order containing a `NaN` panics, and a panic in a paint loop is a
/// closed window. There is no depth for a point the projection could not place, and drawing it
/// first is the answer that loses least.
#[test]
fn a_depth_that_is_not_a_number_sorts_rather_than_panics() {
    let order = editor_core::far_to_near(&[2.0, f64::NAN, 5.0, 1.0]);
    // Furthest first, and the one with no depth first of all: in painter's order that means drawn
    // first, so anything with a real depth paints over it rather than the other way round.
    assert_eq!(
        order,
        vec![1, 2, 0, 3],
        "NaN to the far end, the rest in order"
    );
}

/// **A larger value is never darker**, in the app as well as in the report.
///
/// The scale the shells used covered 16.6 L\* — 44.1 to 60.7 — and carried its signal on the
/// blue-red axis, which is the pair a colour-vision deficiency collapses and a greyscale print
/// loses entirely. This checks the property on the colours `value_colour` actually returns, so it
/// fails if the library's scale is swapped for something that does not keep it.
#[test]
fn the_fallback_scale_climbs_in_lightness() {
    let scale = Some((293.0, 353.0));
    let mut worst_fall = 0.0f64;
    let mut prev = pantometry::view::ramp::lightness(editor_core::value_colour(293.0, scale));
    let first = prev;
    let mut last = prev;
    for i in 1..=200 {
        let v = 293.0 + 60.0 * i as f64 / 200.0;
        let l = pantometry::view::ramp::lightness(editor_core::value_colour(v, scale));
        worst_fall = worst_fall.max(prev - l);
        prev = l;
        last = l;
    }
    assert!(
        worst_fall < 0.2,
        "lightness fell by {worst_fall} L* somewhere"
    );
    assert!(
        last - first > 60.0,
        "the scale covers only {} L*, and the one it replaced covered 16.6",
        last - first
    );
}

/// **A signed field's neutral colour is at the value zero, not at the middle of its range.**
///
/// For -100 to +300 those are a quarter of the scale apart, and colouring the midpoint neutral
/// draws +100 as though it were the datum.
#[test]
fn zero_is_the_neutral_of_a_signed_scale_wherever_it_falls() {
    let lopsided = Some((-100.0, 300.0));
    assert!(editor_core::scale_is_signed(lopsided));

    let chroma = |v: f64| {
        let lab = pantometry::view::ramp::to_lab(editor_core::value_colour(v, lopsided));
        lab[1].hypot(lab[2])
    };
    assert!(chroma(0.0) < 2.0, "zero carries chroma {}", chroma(0.0));
    assert!(
        chroma(100.0) > 8.0,
        "the range midpoint should be well off neutral, not the datum: {}",
        chroma(100.0)
    );

    // And equal deflections either way are equally far from neutral.
    let l = |v: f64| pantometry::view::ramp::lightness(editor_core::value_colour(v, lopsided));
    assert!(
        (l(-80.0) - l(80.0)).abs() < 0.5,
        "the two arms are L* {} and {}",
        l(-80.0),
        l(80.0)
    );

    // A one-sided range is not signed, however wide it is.
    assert!(!editor_core::scale_is_signed(Some((0.0, 4e6))));
    assert!(!editor_core::scale_is_signed(Some((293.0, 353.0))));
}

/// **The colour bar is readable back to a number**, which is the only thing that makes it a bar
/// and not a decoration.
#[test]
fn the_bar_inverts_the_scale_it_shows() {
    let one_sided = Some((293.0, 353.0));
    assert!((editor_core::bar_value(0.0, one_sided) - 293.0).abs() < 1e-9);
    assert!((editor_core::bar_value(1.0, one_sided) - 353.0).abs() < 1e-9);
    assert!((editor_core::bar_value(0.5, one_sided) - 323.0).abs() < 1e-9);

    // A signed bar spans the deflection, symmetric about zero, because the scale is.
    let lopsided = Some((-100.0, 300.0));
    assert!(
        (editor_core::bar_value(0.5, lopsided)).abs() < 1e-9,
        "zero at the middle"
    );
    assert!((editor_core::bar_value(1.0, lopsided) - 300.0).abs() < 1e-9);
    assert!((editor_core::bar_value(0.0, lopsided) + 300.0).abs() < 1e-9);
}

/// **A reading keeps its magnitude.** `{:.4}` printed a cavity's 3.19e-10 J as `0.0000`.
#[test]
fn a_reading_that_is_very_small_is_not_printed_as_zero() {
    for v in [3.19e-10, -7.5e-9, 4.2e7, 1.6e-19] {
        let s = editor_core::magnitude(v);
        let back: f64 = s
            .parse()
            .unwrap_or_else(|_| panic!("{s} does not parse back"));
        assert!(
            back != 0.0 && (back / v - 1.0).abs() < 1e-3,
            "{v:e} printed as {s}"
        );
    }
    // And the ordinary sizes stay readable rather than all becoming exponents.
    assert_eq!(editor_core::magnitude(0.0), "0");
    assert_eq!(editor_core::magnitude(20.25), "20.2500");
    assert_eq!(editor_core::magnitude(1013.25), "1013.25");
    assert_eq!(editor_core::magnitude(f64::NAN), "-");
}
