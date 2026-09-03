//! **Where sample `i` sits in the world**, which two pictures of one run have to agree about.
//!
//! `viewer_core::field_points` is the only answer now: the editor's splats and the viewer's
//! crosses both walk through it. Before that the editor had the loop and the viewer was about to
//! get a second copy, and two copies of this arithmetic are two pictures that can disagree about
//! where a hot spot is.
//!
//! Checked against the closed form rather than against a rendering. `every_scene_can_be_drawn`
//! counts lit pixels and says in its own header that it cannot see whether the picture is *right*
//! — a renderer that put every sample at the origin would pass it, which was measured by
//! sabotaging exactly that. This is what notices.

use viewer_core::field_points;

/// A unit box, corner 0 at the origin, one axis per bit of the corner index.
fn unit_box() -> [[f64; 3]; 8] {
    let mut corners = [[0.0; 3]; 8];
    for (i, c) in corners.iter_mut().enumerate() {
        *c = [(i & 1) as f64, ((i >> 1) & 1) as f64, ((i >> 2) & 1) as f64];
    }
    corners
}

/// **A sample sits at `i / (n - 1)` along its axis, and the ends land on the box.**
#[test]
fn the_first_sample_is_the_low_corner_and_the_last_is_the_high_one() {
    let corners = unit_box();
    let (nx, ny, nz) = (3, 2, 2);
    let values: Vec<f64> = (0..nx * ny * nz).map(|i| i as f64).collect();
    let out = field_points(&corners, (nx, ny, nz), &values, 1);

    assert_eq!(out.len(), nx * ny * nz, "one point per sample");
    assert_eq!(out[0].0, [0.0, 0.0, 0.0], "sample 0 is corner 0");
    assert_eq!(out[0].1, 0.0, "and carries its own value");
    assert_eq!(
        out[out.len() - 1].0,
        [1.0, 1.0, 1.0],
        "the last sample is the opposite corner"
    );

    // x is fastest, so index 1 is the middle of three along x and nothing else has moved.
    assert_eq!(out[1].0, [0.5, 0.0, 0.0], "x is the fast axis");
    // Then y, then z.
    assert_eq!(out[nx].0, [0.0, 1.0, 0.0], "y is next");
    assert_eq!(out[nx * ny].0, [0.0, 0.0, 1.0], "z is slowest");
}

/// **One sample along an axis sits in the middle of it**, which is the only place it can be.
#[test]
fn a_single_sample_along_an_axis_is_at_its_middle() {
    let out = field_points(&unit_box(), (1, 1, 1), &[7.0], 1);
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].0, [0.5, 0.5, 0.5]);
    assert_eq!(out[0].1, 7.0);
}

/// **A cell the domain does not occupy is skipped**, and that is how a hole reaches a picture.
///
/// JSON has no NaN, so `pantometry-view` writes a non-finite value as `null` and this crate reads
/// it back as `f64::NAN`. A sample with no value has no position worth drawing.
#[test]
fn a_hole_is_not_a_point() {
    let corners = unit_box();
    let values = [1.0, f64::NAN, 3.0, f64::INFINITY, 5.0, 6.0, 7.0, 8.0];
    let out = field_points(&corners, (2, 2, 2), &values, 1);
    assert_eq!(out.len(), 6, "two of the eight are not finite");
    assert!(
        out.iter().all(|(_, v)| v.is_finite()),
        "a non-finite value reached the output"
    );
}

/// **A stride takes every `stride`-th sample along each axis**, and 0 means 1.
#[test]
fn a_stride_thins_each_axis_and_never_takes_nothing() {
    let corners = unit_box();
    let (nx, ny, nz) = (5, 5, 5);
    let values = vec![1.0; nx * ny * nz];
    assert_eq!(field_points(&corners, (nx, ny, nz), &values, 1).len(), 125);
    // 0, 2, 4 along each axis.
    assert_eq!(field_points(&corners, (nx, ny, nz), &values, 2).len(), 27);
    // A stride of zero would step by nothing and never finish; it is read as one.
    assert_eq!(field_points(&corners, (nx, ny, nz), &values, 0).len(), 125);
}

/// **A grid that does not match its values is refused rather than read past its end.**
#[test]
fn a_short_value_array_draws_nothing() {
    let corners = unit_box();
    assert!(field_points(&corners, (2, 2, 2), &[1.0, 2.0], 1).is_empty());
    assert!(field_points(&corners, (0, 2, 2), &[1.0; 4], 1).is_empty());
}

/// **The box is the placed one, so a turned domain's samples turn with it.**
///
/// Corner 0 is the origin and the axes come from corners 1, 2 and 4, so a box handed in rotated
/// hands back points in the rotated frame — the placement is applied before this is called and
/// there is nothing here that assumes the axes are the world's.
#[test]
fn the_samples_follow_the_box_they_are_given() {
    // The unit box turned a quarter turn about z: x goes to y, y goes to -x.
    let turned = |p: [f64; 3]| [-p[1], p[0], p[2]];
    let mut corners = [[0.0; 3]; 8];
    for (i, c) in corners.iter_mut().enumerate() {
        *c = turned([(i & 1) as f64, ((i >> 1) & 1) as f64, ((i >> 2) & 1) as f64]);
    }
    let out = field_points(&corners, (2, 1, 1), &[1.0, 2.0], 1);
    assert_eq!(out.len(), 2);

    // Derived rather than guessed, because the first draft of this guessed and was wrong. The
    // axes are `corner1 - corner0`, `corner2 - corner0` and `corner4 - corner0`, which under this
    // turn are `[0, 1, 0]`, `[-1, 0, 0]` and `[0, 0, 1]`. The single sample along y and along z
    // sits at the middle of *those* axes, so both samples carry `ay * 0.5 = [-0.5, 0, 0]` and
    // `az * 0.5 = [0, 0, 0.5]` — the displacement the first draft forgot — and the two differ
    // only by one whole step along the turned x.
    assert_eq!(out[0].0, [-0.5, 0.0, 0.5]);
    assert_eq!(out[1].0, [-0.5, 1.0, 0.5]);
}
