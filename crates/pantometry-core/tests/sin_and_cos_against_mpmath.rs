//! `pantometry_core::math::sin_cos_turns` against mpmath at 50 significant digits.
//!
//! The references are `tests/data/{sin,cos}-reference.bin`, written by
//! `tools/math-reference/generate.py` as mpmath's `sinpi(2t)` and `cospi(2t)` over one list of
//! `t`, which is drawn here by the script's recipe and checked against the digest each table
//! carries, as `exp_and_ln_against_mpmath.rs` does for its own; the comparison is the same code.
//!
//! **The bounds are derived from the method, not fitted to what it measured**: 0.501 ulp where
//! the result comes from sin's series and 0.503 where it comes from cos's — see [`sin_bound`] and
//! [`cos_bound`] — and nothing at all at and above `2⁵²`, where the result is exact.

use pantometry_core::math::sin_cos_turns;

mod reference;
use reference::*;

/// `generate.py`'s `turn_inputs`, line for line.
fn turn_inputs(scale: u64) -> Vec<u64> {
    let mut xs = Vec::new();
    // What the samplers pass: `Rng::unit`'s values.
    xs.extend((0..16384 * scale).map(|i| b((splitmix(21, i) >> 11) as f64 * pow2(-53))));
    xs.extend(uniform(b(pow2(-60)), b(pow2(53)), 22, 16384 * scale, 2));
    xs.extend(
        uniform(0, MASK, 23, 4096 * scale, 0)
            .into_iter()
            .filter(|&x| f64::from_bits(x).is_finite()),
    );
    // Every reduction boundary and every zero, |m| ≤ 256.
    for m in -256i32..=256 {
        xs.extend(around(b(f64::from(2 * m + 1) / 8.0), 4, &[0]));
        xs.extend(around(b(f64::from(m) / 4.0), 4, &[0]));
    }
    for e in 0..52 {
        xs.extend(around(b(pow2(e)), 4, &[0, 1]));
        if e < 50 {
            xs.extend(around(b(pow2(e) + 0.125), 4, &[0, 1]));
        }
    }
    xs.extend(uniform(1, 0x000F_FFFF_FFFF_FFFF, 24, 2048 * scale, 2));
    xs.extend(around(
        b(pow2(-1022) / (2.0 * std::f64::consts::PI)),
        64,
        &[0, 1],
    ));
    xs.extend(around(b(pow2(-502)), 64, &[0, 1]));
    xs.extend(around(b(pow2(52)), 64, &[0, 1]));
    // The worst of the run at twenty times the random points, one from each series.
    xs.extend([0x3FD7_7AD6_59A0_A538, 0x3FBB_8EED_AA51_E530]);
    xs.into_iter()
        .filter(|&x| f64::from_bits(x).is_finite())
        .collect()
}

/// Which series gives sin of `t`: `true` for sin's, in an even quadrant. The reduction is
/// written out again here, `4|t| = n + u` with `u` in `(−½, ½]`, rather than read from `math`.
fn sin_series_for_sin(t: f64) -> bool {
    if t.abs() >= pow2(52) {
        // Exact there, by either name; and 4|t| as an integer would not fit.
        return true;
    }
    let x = 4.0 * t.abs();
    let n = x as u64;
    let n = if x - n as f64 > 0.5 { n + 1 } else { n };
    n % 2 == 0
}

/// **sin's series, 0.501 ulp, derived.** `math`'s `sin_quarter` computes `sin((π/2) u)`, `|u| ≤
/// ½`, as `S₁u + S₃u³ + S₅u⁵ + u⁷ P(u²)`, rounded once at the end. Every other error is bounded
/// below as a fraction of the result, in units of `2⁻⁵³`, and is therefore at most that many ulps,
/// since an ulp is more than `2⁻⁵³` of the value. Each fraction is largest at `|u| = ½`, where
/// the angle is `π/4` and the result `sin(π/4)`:
///
/// - **the first three terms are exact**, less what is `2⁻¹⁰⁴` or smaller: `u²`, `u³` and `u⁵`
///   are formed as doubles and their errors by Dekker's product, each product with a
///   coefficient's high part likewise, and the three summed by two exact additions — `|S₁u| >
///   |S₃u³| > |S₅u⁵|`. What is left is the coefficients' low parts (`2⁻¹⁰⁶`), the products that
///   carry the low parts (`2⁻⁵³` of a term `2⁻⁵²` of the result), and the products of two low
///   parts, which are dropped (`2⁻¹⁰⁴`);
/// - **the tail** `u⁷ P(u²)` is at most `(π/4)⁷/5040/sin(π/4) = 5.2 · 10⁻⁵` of the result.
///   It is computed from `u⁵` and `u²` rounded to doubles (three roundings in `u⁵`, and `u²`'s
///   and the product's in `u⁷`), by Horner's rule (one rounding of the leading sum, one of the
///   coefficient `S₇`, and `u²`'s error carried in at `(π/4)²/72`), and multiplied once: eight
///   roundings, `8.1 · 2⁻⁵³` of the tail, **0.00042 ulp**;
/// - **the sum the tail is added into** rounds once: `2⁻⁵³` of the tail and of the exact
///   additions' errors, **0.00005 ulp**;
/// - **the truncation** after `u¹⁹`, the next term `(π/4)²¹/21!` against `sin(π/4)`: `1.7 ·
///   10⁻²²`, **0.000002 ulp**.
///
/// So 0.5005, held at 0.501. Below `|u| = 2⁻⁵⁰⁰` the result is `(π/2) u` — `u³` is `2⁻¹⁰⁰⁰` of
/// it — as a double-double rounded once, onto the subnormal grid where the result is subnormal:
/// within `0.5 + 2⁻⁵²`.
fn sin_bound(_: f64) -> (&'static str, f64) {
    ("sin's series", 0.501)
}

/// **cos's series, 0.503 ulp, derived.** `cos_quarter` computes `1 + C₂u² + C₄u⁴ + u⁶ R(u²)` with
/// the first three terms exact as above, and the result is in `[cos(π/4), 1]`, where an ulp is
/// `2⁻⁵³` absolutely (it is exactly 1 only at `u = 0`, exactly). At `|u| = ½`:
///
/// - **the tail** is at most `(π/4)⁶/720 = 3.3 · 10⁻⁴`, and is computed from `u⁴` and `u²` as
///   doubles (three roundings in `u⁴`, two more in `u⁶`), by Horner's rule (two, and `u²`'s error
///   at `(π/4)²/56`), and one product: `8.1 · 2⁻⁵³` of the tail with the second-order terms,
///   **0.0026 ulp**;
/// - **the sum it is added into**, `2⁻⁵³` of the tail: **0.0003 ulp**;
/// - **the truncation** after `u²⁰`, `(π/4)²²/22! = 4.4 · 10⁻²⁴`: nothing to speak of.
///
/// So 0.50297, held at 0.503. The tail is where the bound comes from, and it is ten times sin's
/// because cos's series starts its tail one power lower.
fn cos_bound(_: f64) -> (&'static str, f64) {
    ("cos's series", 0.503)
}

fn by_series(t: f64, sin_series: bool) -> (&'static str, f64) {
    if t.abs() >= pow2(52) {
        ("exact, |t| ≥ 2^52", 0.0)
    } else if sin_series {
        sin_bound(t)
    } else {
        cos_bound(t)
    }
}

fn sin_regions(t: f64) -> (&'static str, f64) {
    by_series(t, sin_series_for_sin(t))
}

fn cos_regions(t: f64) -> (&'static str, f64) {
    by_series(t, !sin_series_for_sin(t))
}

/// sin and cos of `2πt` are exact only where `4t` is an integer — every `|t| ≥ 2⁵²` among them.
/// Elsewhere they are irrational at a rational `t` except at sixths and twelfths of a turn, where
/// sin or cos is `±½` (Niven), and no double is one of those.
fn sin_is_exact(t: f64) -> bool {
    let x = 4.0 * t.abs();
    t.abs() >= pow2(52) || x == x.trunc()
}

fn cos_is_exact(t: f64) -> bool {
    sin_is_exact(t)
}

fn sin(t: f64) -> f64 {
    sin_cos_turns(t).0
}

fn cos(t: f64) -> f64 {
    sin_cos_turns(t).1
}

/// **sin 2πt is within its bounds of mpmath at every point**: what `Rng::unit` returns, the whole
/// finite domain, every reduction boundary `(2m + 1)/8` and every zero `m/4` for `|m| ≤ 256`, the
/// same about every power of two that leaves a fraction, the subnormals, where the result turns
/// subnormal, where the small path begins, `2⁵²`, and the two worst points of a run at twenty
/// times the random points. And the comparison can see one ulp, and the table marks exactly the
/// exact values, as for `exp` and `ln`.
#[test]
fn sin_is_within_its_bounds_of_mpmath() {
    let (inputs, table) = (turn_inputs(1), include_bytes!("data/sin-reference.bin"));
    compare("sin", &sin, &inputs, table, Some(sin_regions));
    let shifted = compare(
        "sin moved up one ulp",
        &one_ulp_up(sin),
        &inputs,
        table,
        None,
    );
    assert!(shifted.worst() > 1.0, "one ulp up is still within one ulp");
    table_marks_exact_values("sin", &inputs, table, sin_is_exact);
}

/// **cos 2πt is within its bounds of mpmath** at the same points, with the same two checks of the
/// comparison.
///
/// Neither comparison can see the budget beyond the last rounding's ½ except near a tie, so the
/// 0.503 of cos's series is never needed by a committed point: every one of them, and every one
/// of the run at twenty times, is within 0.501. `math`'s
/// `each_series_is_within_its_budget_before_rounding` holds the budgets themselves.
#[test]
fn cos_is_within_its_bounds_of_mpmath() {
    let (inputs, table) = (turn_inputs(1), include_bytes!("data/cos-reference.bin"));
    compare("cos", &cos, &inputs, table, Some(cos_regions));
    let shifted = compare(
        "cos moved up one ulp",
        &one_ulp_up(cos),
        &inputs,
        table,
        None,
    );
    assert!(shifted.worst() > 1.0, "one ulp up is still within one ulp");
    table_marks_exact_values("cos", &inputs, table, cos_is_exact);
}

/// The same against a larger table, written outside the repository by `generate.py --scale`, as
/// in `exp_and_ln_against_mpmath.rs`. Ignored, and not run in CI, because the table is not
/// committed.
#[cfg(not(target_family = "wasm"))]
#[test]
#[ignore = "needs a table written by generate.py --scale: see exp_and_ln_against_mpmath.rs"]
fn against_a_larger_table() {
    let dir = std::env::var("PANTOMETRY_MATH_REFERENCE").expect("PANTOMETRY_MATH_REFERENCE");
    let scale: u64 = std::env::var("PANTOMETRY_MATH_REFERENCE_SCALE")
        .expect("PANTOMETRY_MATH_REFERENCE_SCALE")
        .parse()
        .unwrap();
    let read = |name: &str| std::fs::read(format!("{dir}/{name}-reference.bin")).unwrap();
    let inputs = turn_inputs(scale);
    compare("sin", &sin, &inputs, &read("sin"), Some(sin_regions));
    compare("cos", &cos, &inputs, &read("cos"), Some(cos_regions));
}

/// **The cost**, in release, over `Rng::unit`'s values, beside the platform's `sin` and `cos` of
/// `2πt`. Asserted: nothing, since it is a property of the machine; printed.
#[cfg(not(target_family = "wasm"))]
#[test]
#[ignore = "a timing: run with --release -- --ignored --nocapture"]
fn the_cost() {
    use std::hint::black_box;
    use std::time::Instant;
    let args: Vec<f64> = (0..4096u64)
        .map(|i| (splitmix(51, i) >> 11) as f64 * pow2(-53))
        .collect();
    let time = |f: &dyn Fn(f64) -> f64| -> f64 {
        let rounds = 2000;
        let mut best = f64::INFINITY;
        for _ in 0..5 {
            let t = Instant::now();
            let mut sum = 0.0;
            for _ in 0..rounds {
                for &x in black_box(&args) {
                    sum += f(x);
                }
            }
            black_box(sum);
            best = best.min(t.elapsed().as_secs_f64() * 1e9 / (rounds as f64 * args.len() as f64));
        }
        best
    };
    let ours = time(&|t| {
        let (s, c) = sin_cos_turns(t);
        s + c
    });
    let platform = time(&|t: f64| {
        let (s, c) = (std::f64::consts::TAU * t).sin_cos();
        s + c
    });
    println!("ns per call, sin and cos together: ours {ours:.2}, platform {platform:.2}");
}
