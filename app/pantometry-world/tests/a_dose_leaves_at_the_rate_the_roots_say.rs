//! The bi-exponential's slow root, and the order at which the scene's step reaches it.
//!
//! A two-compartment model's plasma curve is `A e^{−αt} + B e^{−βt}`, and `α` and `β` are the roots
//! of
//!
//! ```text
//!   λ² − (k₁₀ + k₁₂ + k₂₁) λ + k₁₀ k₂₁ = 0
//! ```
//!
//! computed from the scene's own numbers and nothing else. For
//! `32-a-dose-distributing-and-leaving` — 5 L of plasma clearing at 3 L/h, 20 L of tissue, 8 L/h
//! between them — that is `k₁₀ = 0.6`, `k₁₂ = 1.6`, `k₂₁ = 0.4` per hour, giving `α = 2.504 /h`
//! and `β = 0.0958 /h`. After four hours the fast term is `4.5e-5` of where it started, so the
//! late frames are a single exponential and their log-linear slope is `β`.
//!
//! # Why the check is a rate and not a tolerance
//!
//! `CompartmentModel` steps with explicit Euler, so its error is first order in the step and a
//! single number would be a number somebody picked. Measured here instead: refine the step and the
//! disagreement with `β` has to fall **in proportion**, which is a statement about the scheme that
//! no single run can make.
//!
//! # `max_stable_dt` is not an accuracy limit
//!
//! The model reports `1/turnover_rate`, which for this scene's plasma compartment is
//! `1 / ((3 + 8)/5 per hour)` = **1636 s**. The scene as first written stepped at `43200/13` =
//! 3323 s — over twice that — and read `β` 1.5% high. But stepping *at* the limit would still read
//! it badly: `α · 1636 s` is `1.14`, and a first-order scheme wants that number small rather than
//! merely under one. The shipped `window_s` is 10 s.
//!
//! # Not under `wasm32`
//!
//! It reads the shipped scene off a disk.

#![cfg(not(target_family = "wasm"))]

use pantometry_world::{DomainSpec, OnDisk, Scene, World};

/// The scene, with its window replaced.
fn scene_with(window_s: Option<f64>) -> Scene {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("scenes/32-a-dose-distributing-and-leaving.json");
    let text = std::fs::read_to_string(&path).expect("the shipped scene reads");
    let mut scene: Scene = serde_json::from_str(&text).expect("the shipped scene parses");
    scene.window_s = window_s;
    scene
}

/// `α` and `β`, from the scene's own compartments rather than from anything measured.
fn roots(scene: &Scene) -> (f64, f64) {
    let DomainSpec::Compartments { volumes, links, .. } = &scene.domains[0] else {
        panic!("the scene's first domain is not the compartment model it was written to be");
    };
    let central = &volumes[0];
    let peripheral = &volumes[1];
    let q = links[0].clearance_l_per_h;
    let (k10, k12, k21) = (
        central.clearance_l_per_h / central.volume_l,
        q / central.volume_l,
        q / peripheral.volume_l,
    );
    let b = k10 + k12 + k21;
    let root = (b * b - 4.0 * k10 * k21).sqrt();
    ((b + root) / 2.0, (b - root) / 2.0)
}

/// The log-linear slope of the plasma concentration over the last `n` frames, per hour.
fn terminal_rate(scene: Scene, n: usize) -> f64 {
    let mut world = World::build_with(scene, &OnDisk).expect("the scene builds");
    let mut hours = Vec::new();
    let mut logs = Vec::new();
    let mut take = |w: &World, t: f64| {
        let c = w
            .capture()
            .readings
            .iter()
            .find(|r| r.label == "plasma")
            .expect("a plasma reading")
            .value;
        hours.push(t / 3600.0);
        logs.push(c.ln());
    };
    let steps = world.steps();
    let dt = world.scene().duration_s / steps as f64;
    for i in 1..=steps {
        world
            .advance(pantometry::units::Qty::from_si(dt))
            .expect("the model steps");
        // Only the tail, where the fast term is gone.
        if i * 4 >= steps * 3 && i % (steps / 40).max(1) == 0 {
            take(&world, i as f64 * dt);
        }
    }
    assert!(logs.len() >= 4, "only {} points in the tail", logs.len());
    let n = logs.len().min(n.max(4));
    let (h, l) = (&hours[hours.len() - n..], &logs[logs.len() - n..]);
    let m = n as f64;
    let (sx, sy): (f64, f64) = (h.iter().sum(), l.iter().sum());
    let sxx: f64 = h.iter().map(|x| x * x).sum();
    let sxy: f64 = h.iter().zip(l).map(|(x, y)| x * y).sum();
    -(m * sxy - sx * sy) / (m * sxx - sx * sx)
}

/// **The terminal rate converges on the quadratic's slow root at first order in the step.**
///
/// Explicit Euler is first order, so halving the step has to halve the disagreement. Measured,
/// against `β = 0.09584 /h`:
///
/// | steps | terminal rate | relative error |
/// | --- | --- | --- |
/// | 144 | `0.09623` | `4.0e-3` |
/// | 720 | `0.09592` | `8.0e-4` |
/// | 4320 | `0.09585` | `1.3e-4` |
///
/// Five times the steps divided the error by 5.0, and six times by 6.0. That is the claim: not a
/// tolerance, a rate.
#[test]
fn the_terminal_rate_converges_on_the_slow_root_at_first_order() {
    let (alpha, beta) = roots(&scene_with(None));
    println!("  alpha {alpha:.5} /h, beta {beta:.5} /h, from the scene's own numbers");

    let mut errors = Vec::new();
    for window in [300.0, 60.0, 10.0] {
        let scene = scene_with(Some(window));
        let steps = scene.steps();
        let rate = terminal_rate(scene, 8);
        let error = (rate / beta - 1.0).abs();
        println!("  window {window:5.0} s, {steps:5} steps: {rate:.5} /h, off by {error:.2e}");
        errors.push((steps as f64, error));
    }

    // First order: the error falls in proportion to the step. Asked as a *ratio* of ratios, so
    // nothing here is a number somebody picked — a second-order scheme would give 25 and 36, and
    // a scheme that had stopped converging would give 1.
    for pair in errors.windows(2) {
        let (coarse, fine) = (pair[0], pair[1]);
        let refined = fine.0 / coarse.0;
        let improved = coarse.1 / fine.1;
        let order = improved.ln() / refined.ln();
        println!("  {refined:.1}x the steps improved it {improved:.2}x — order {order:.2}");
        assert!(
            (0.85..1.15).contains(&order),
            "explicit Euler is first order and this measured {order:.2}; a scheme that had \
             stopped converging would read 0 and a second-order one 2"
        );
    }

    // And the shipped window is the finest of the three, so the scene ships at that error.
    let shipped = errors.last().expect("three windows").1;
    assert!(
        shipped < 1e-3,
        "the shipped scene reads the terminal rate {shipped:.2e} from beta, and it was 1.3e-4"
    );
}

/// **The step the scene ships with is well inside the model's own stability limit, and that limit
/// would not have been enough.**
///
/// `max_stable_dt` is `1/turnover_rate` — for this scene's plasma compartment, `1636 s`. The first
/// version of the scene had no `window_s` and stepped at `43200/13 = 3323 s`, over twice the
/// limit, and read the terminal rate 1.5% high. Stepping exactly at the limit is still wrong for a
/// first-order scheme: `α · 1636 s = 1.14`, and the error is proportional to that.
///
/// This pins the distinction, because a scene that merely respected `max_stable_dt` would look
/// correct and read a rate that is 1% out.
#[test]
fn the_stability_limit_is_not_the_accuracy_limit() {
    let scene = scene_with(None);
    let (alpha, _) = roots(&scene);
    let world = World::build_with(scene_with(None), &OnDisk).expect("the scene builds");
    let limit = world
        .simulation()
        .domains()
        .map(|d| {
            d.max_stable_dt(pantometry::units::Qty::from_si(0.0))
                .to_si()
        })
        .fold(f64::INFINITY, f64::min);
    let shipped = scene_with(Some(10.0));
    let dt = shipped.duration_s / shipped.steps() as f64;
    println!(
        "  max_stable_dt {limit:.0} s, shipped step {dt:.0} s, alpha*dt {:.4}",
        alpha / 3600.0 * dt
    );

    assert!(
        (1500.0..1800.0).contains(&limit),
        "the limit is {limit:.0} s and 1/((3+8)/5 per hour) is 1636"
    );
    assert!(
        dt < limit / 100.0,
        "the shipped step is {dt:.0} s against a stability limit of {limit:.0} s, which is not \
         the margin a first-order scheme needs"
    );
    // At the limit itself, the dimensionless step is of order one — which is the point.
    let at_limit = alpha / 3600.0 * limit;
    assert!(
        at_limit > 1.0,
        "alpha times the stability limit is {at_limit:.2}; if it were small, respecting the limit \
         would have been enough and this test would be saying nothing"
    );
}
