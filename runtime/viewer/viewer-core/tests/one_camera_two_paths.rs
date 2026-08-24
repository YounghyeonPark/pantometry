//! **The matrix and the projection are the same camera.**
//!
//! A 2D painter walks points through [`Camera::project`]; a GPU multiplies every vertex by
//! [`Camera::matrix`]. The editor does both in one viewport — shaded surfaces underneath, labels
//! and a colour bar and a probe on top — so a disagreement between them is not a rounding
//! curiosity, it is a caption pointing at empty space.
//!
//! This is a consistency check and says so. It is not a claim about whether the projection is
//! *right*; it is the claim that there is only one of it.

use viewer_core::{Camera, Framing};

/// `m` is column-major, so `m[4 * column + row]`. `p` is already framing-local.
fn apply(m: &[f32; 16], p: [f32; 3]) -> [f64; 4] {
    let mut out = [0.0f64; 4];
    for row in 0..4 {
        out[row] = (0..3)
            .map(|c| m[4 * c + row] as f64 * p[c] as f64)
            .sum::<f64>()
            + m[4 * 3 + row] as f64;
    }
    out
}

/// A spread of cameras and framings, so agreement is not an accident of one arrangement.
fn cases() -> Vec<(Camera, Framing, f64)> {
    let mut out = Vec::new();
    for (az, el) in [(0.0, 0.0), (0.7, 0.4), (-2.1, -1.3), (3.0, 1.49)] {
        for distance in [1.2, 2.5, 9.0] {
            for scale in [0.5, 1.0, 4.0] {
                for (centre, span) in [
                    ([0.0, 0.0, 0.0], 1.0),
                    ([0.04, -0.011, 0.2], 0.009),
                    ([-3.0, 12.0, 0.5], 4.4),
                ] {
                    for aspect in [0.5, 1.0, 2.37] {
                        out.push((
                            Camera {
                                azimuth: az,
                                elevation: el,
                                distance,
                                scale,
                            },
                            Framing { centre, span },
                            aspect,
                        ));
                    }
                }
            }
        }
    }
    out
}

/// Points spread over and outside the framed box, in each framing's own units.
fn points(frame: &Framing) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    for fx in [-0.5, -0.13, 0.0, 0.37, 0.5, 1.4] {
        for fy in [-0.5, 0.0, 0.29, 0.5] {
            for fz in [-0.5, 0.0, 0.5] {
                out.push([
                    frame.centre[0] + fx * frame.span,
                    frame.centre[1] + fy * frame.span,
                    frame.centre[2] + fz * frame.span,
                ]);
            }
        }
    }
    out
}

/// `x/w` and `y/w` from the matrix are `project`'s `x` and `y`; `w` is its depth.
#[test]
fn the_matrix_and_the_projection_agree() {
    let mut worst = 0.0f64;
    let mut worst_at = String::new();
    let mut checked = 0usize;

    for (camera, frame, aspect) in cases() {
        let m = camera.matrix(aspect);
        for p in points(&frame) {
            let q = camera.project(p, &frame, aspect);
            let c = apply(&m, frame.local(p));

            // `project` clamps its divisor at the near plane where the matrix does not, because a
            // painter has nothing to clip against. Those points are the near-plane cases and the
            // GPU discards them; comparing them would be comparing a clamp to a clip.
            if c[3] <= Camera::NEAR {
                continue;
            }
            checked += 1;

            // Relative to the scale of the picture, which is order one — an absolute tolerance
            // here would be a tolerance about `scale`, and `scale` runs over a factor of eight.
            for (got, want) in [(c[0] / c[3], q.x), (c[1] / c[3], q.y), (c[3], q.depth)] {
                let e = (got - want).abs() / want.abs().max(1.0);
                if e > worst {
                    worst = e;
                    worst_at = format!("{camera:?} {frame:?} aspect {aspect} at {p:?}");
                }
            }
        }
    }

    println!("  {checked} points, worst relative disagreement {worst:.3e}");
    assert!(
        checked > 5_000,
        "the sweep has to actually sweep: {checked}"
    );
    // The two expressions differ only in the order they multiply the same terms, and the matrix
    // rounds to `f32` on the way out. `f32` carries about `6e-8`, and a few operations of it is
    // the whole budget — this is not a physics tolerance, it is the width of the type.
    assert!(
        worst < 1e-6,
        "the matrix and the projection are different cameras: {worst:.3e} at {worst_at}"
    );
}

/// Depth into the buffer: `-1` at the near plane, `1` at the far one, and monotone between.
#[test]
fn the_depth_buffer_runs_the_right_way() {
    let camera = Camera::default();
    let m = camera.matrix(1.0);

    // Straight down the view direction, so `w` is the only thing changing.
    let ndc_at = |depth: f64| {
        // A point `t` along the view direction has `w = t + distance`, so solve for `t`. The view
        // direction is the matrix's own `w` row, which is the third rotation row.
        let t = depth - camera.distance;
        let dir = [m[3], m[7], m[11]];
        let c = apply(
            &m,
            [dir[0] * t as f32, dir[1] * t as f32, dir[2] * t as f32],
        );
        c[2] / c[3]
    };

    let near = ndc_at(Camera::NEAR);
    let far = ndc_at(camera.far());
    assert!(
        (near + 1.0).abs() < 1e-5,
        "the near plane is not at -1: {near}"
    );
    assert!((far - 1.0).abs() < 1e-5, "the far plane is not at 1: {far}");

    let mut previous = f64::NEG_INFINITY;
    for i in 0..64 {
        let depth = Camera::NEAR + (camera.far() - Camera::NEAR) * i as f64 / 63.0;
        let z = ndc_at(depth);
        assert!(
            z > previous,
            "depth {depth} went backwards: {z} after {previous}"
        );
        previous = z;
    }
}
