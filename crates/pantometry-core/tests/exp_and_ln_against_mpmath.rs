//! `pantometry_core::math::{exp, ln}` against mpmath at 50 significant digits.
//!
//! The reference is `tests/data/{exp,ln}-reference.bin`, written by
//! `tools/math-reference/generate.py`: for each point, the correctly rounded value and where the
//! exact value lies from it, in ulps. **The points themselves are not stored**: they are drawn
//! here by the script's own recipe of integer operations on bit patterns, and each table carries a
//! digest of the inputs it was made for, which is checked first. A recipe that drifted on one side
//! fails by name rather than comparing one point's reference with another point's value.
//!
//! The error is `(computed − exact) / ulp(exact)`, `ulp` the spacing of doubles at the exact
//! value (`2⁻¹⁰⁷⁴` below the smallest normal), and **the bounds asserted are derived from the
//! methods, not fitted to what they measured**: 0.520 ulp for `exp`; for `ln`, 0.508 where the
//! binade is that of 1 (`k = 0`) and 0.501 elsewhere — see [`exp_bound`] and [`ln_bound`]. The
//! measured worst case of each region is printed, and is the figure in `math`'s documentation.
//!
//! An ignored test does the same against a larger table, written elsewhere by `--scale`; see
//! [`against_a_larger_table`].

use pantometry_core::math::{exp, ln};

const LN_2: f64 = std::f64::consts::LN_2;
const MASK: u64 = u64::MAX;

fn splitmix(seed: u64, i: u64) -> u64 {
    let mut z = seed.wrapping_add((i + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// `n` bit patterns uniform in `[lo, hi]`; `sign` 0 or 1 is the sign bit, 2 alternates.
fn uniform(lo: u64, hi: u64, seed: u64, n: u64, sign: u64) -> Vec<u64> {
    let span = (hi - lo).wrapping_add(1);
    (0..n)
        .map(|i| {
            let r = splitmix(seed, i);
            let b = lo + if span == 0 { r } else { r % span };
            let s = if sign == 2 { i & 1 } else { sign };
            b | (s << 63)
        })
        .collect()
}

fn around(base: u64, reach: i64, signs: &[u64]) -> Vec<u64> {
    let mut out = Vec::new();
    for &s in signs {
        for j in -reach..=reach {
            out.push((base.wrapping_add(j as u64) & MASK) | (s << 63));
        }
    }
    out
}

fn b(x: f64) -> u64 {
    x.to_bits()
}

/// `generate.py`'s `exp_inputs`, line for line.
fn exp_inputs(scale: u64) -> Vec<u64> {
    let mut xs = Vec::new();
    xs.extend(uniform(b(pow2(-30)), b(746.0), 1, 16384 * scale, 2));
    xs.extend(uniform(b(1.0), b(709.8), 2, 8192 * scale, 0));
    xs.extend(uniform(b(1.0), b(745.2), 3, 8192 * scale, 1));
    xs.extend(
        uniform(0, MASK, 4, 4096 * scale, 0)
            .into_iter()
            .filter(|&x| !f64::from_bits(x).is_nan()),
    );
    for k in -1075..1025 {
        xs.extend(around(b((f64::from(k) + 0.5) * LN_2), 4, &[0]));
    }
    for base in [
        0x3FD6_2E43_0000_0000,
        0x3FF0_A2B2_0000_0000,
        0x3E30_0001_0000_0000,
        0x4086_232B_0000_0000,
        b(708.396_418_532_264_1),
        0x4086_2E42_FEFA_39EF,
        0x4087_4910_D52D_3051,
    ] {
        xs.extend(around(base, 64, &[0, 1]));
    }
    xs.extend(uniform(
        b(708.396_418_532_264_1),
        b(745.133_219_101_941_2),
        5,
        2048 * scale,
        1,
    ));
    // The seams of `math`'s method the lines above do not reach: every boundary of the table
    // index, `(n + ½) ln 2/128` for `|n| ≤ 512`, and `|x| = 708`, where the main path ends.
    for n in -512..512 {
        xs.extend(around(b((f64::from(n) + 0.5) * LN_2 / 128.0), 2, &[0]));
    }
    xs.extend(around(b(708.0), 64, &[0, 1]));
    xs
}

/// `generate.py`'s `ln_inputs`, line for line.
fn ln_inputs(scale: u64) -> Vec<u64> {
    let mut xs = Vec::new();
    xs.extend(uniform(1, 0x7FEF_FFFF_FFFF_FFFF, 11, 16384 * scale, 0));
    xs.extend(uniform(1, 0x000F_FFFF_FFFF_FFFF, 12, 2048 * scale, 0));
    xs.extend(around(b(1.0), 4096, &[0]));
    xs.extend(uniform(b(0.9375), b(1.0625), 13, 8192 * scale, 0));
    for e in 1u64..2047 {
        xs.extend(around((e << 52) | (0x6A09E << 32), 3, &[0]));
    }
    xs.extend(uniform(b(1.0), b(64.0), 14, 8192 * scale, 0));
    xs.extend(around(0x0010_0000_0000_0000, 64, &[0]));
    xs.extend((0..129).map(|j| 0x7FEF_FFFF_FFFF_FFFF - j));
    xs.extend(1..65);
    // Every power of two, 2^-1074 to 2^1023, by its bits.
    xs.extend((0..52).map(|j| 1u64 << j));
    xs.extend((1u64..2047).map(|e| e << 52));
    // Every edge of the first implementation's 128 cells, k = -1, 0 and 1: kept as points.
    for k in [-1i64, 0, 1] {
        for i in 0u64..129 {
            let base = (0x3FE6_0000_0000_0000u64 + (i << 45)).wrapping_add((k << 52) as u64);
            xs.extend(around(base, 8, &[0]));
        }
    }
    // The seams of `math`'s method: every edge of its 256 cells at k = −1, 0 and 1, and where k
    // steps, at 0.703125 · 2^e, in every binade.
    for k in [-1i64, 0, 1] {
        for i in 0u64..256 {
            let base =
                (0x3FE6_8000_0000_0000u64 + (i << 44) + (1 << 43)).wrapping_add((k << 52) as u64);
            xs.extend(around(base, 2, &[0]));
        }
    }
    for e in 1u64..2047 {
        xs.extend(around((e << 52) | 0x0006_8000_0000_0000, 1, &[0]));
    }
    xs
}

fn fnv1a(words: &[u64]) -> u64 {
    let mut h: u64 = 0xCBF2_9CE4_8422_2325;
    for w in words {
        for byte in w.to_le_bytes() {
            h ^= u64::from(byte);
            h = h.wrapping_mul(0x0000_0100_0000_01B3);
        }
    }
    h
}

/// The table: its input digest, and per point the correctly rounded bits and the exact value's
/// offset from them in ulps.
fn parse(bytes: &[u8]) -> (u64, Vec<(u64, f64)>) {
    assert_eq!(&bytes[..8], b"PMREF1\0\0", "not a reference table");
    let digest = u64::from_le_bytes(bytes[8..16].try_into().unwrap());
    let n = u32::from_le_bytes(bytes[16..20].try_into().unwrap()) as usize;
    let body = &bytes[20..];
    assert_eq!(body.len(), n * 10, "the table's length is not its count");
    // By index rather than `as_chunks`, which is newer than the 1.78 floor.
    let rows = (0..n)
        .map(|k| &body[10 * k..10 * k + 10])
        .map(|c| {
            let hi = u64::from_le_bytes(c[..8].try_into().unwrap());
            let frac = f64::from(i16::from_le_bytes([c[8], c[9]])) / 32768.0;
            (hi, frac)
        })
        .collect();
    (digest, rows)
}

/// The spacing of doubles at the exact value, rebuilt from the rounded value `hi` and the sign of
/// `frac = (exact − hi)/ulp`: half `hi`'s own spacing when `hi` is a power of two and the exact
/// value lies below it in magnitude, where the doubles are twice as dense.
fn ulp_of_exact(hi: f64, frac: f64) -> f64 {
    let bits = hi.to_bits() & 0x7FFF_FFFF_FFFF_FFFF;
    let e = (bits >> 52) as i32;
    if e <= 1 {
        // Subnormal, zero, or the smallest normal binade: the spacing is 2^-1074 throughout.
        return f64::from_bits(1);
    }
    // 2^(e − 1075): a normal double's bits above the 53rd binade, a subnormal's below it.
    let u = if e > 52 {
        f64::from_bits(((e - 52) as u64) << 52)
    } else {
        f64::from_bits(1u64 << (e - 1))
    };
    let below = frac != 0.0 && (frac < 0.0) != (hi < 0.0);
    if bits & 0x000F_FFFF_FFFF_FFFF == 0 && below {
        0.5 * u
    } else {
        u
    }
}

/// Where a function's error is bounded, and by how much: a name and a bound in ulps, per input.
type Bound = fn(f64) -> (&'static str, f64);

/// What a comparison found: per region, its worst error and where; and how many were not
/// correctly rounded.
struct Found {
    worst: Vec<(&'static str, f64, f64)>,
    incorrect: usize,
    points: usize,
}

impl Found {
    fn worst(&self) -> f64 {
        self.worst.iter().map(|w| w.1).fold(0.0, f64::max)
    }
}

/// `f` at every input against the table. With a `bound`, each point is held to its region's bound;
/// without one the comparison only measures, for the check that it can see one ulp.
fn compare(
    name: &str,
    f: &dyn Fn(f64) -> f64,
    inputs: &[u64],
    table: &[u8],
    bound: Option<Bound>,
) -> Found {
    let (digest, rows) = parse(table);
    assert_eq!(
        fnv1a(inputs),
        digest,
        "{name}: the inputs drawn here are not the ones the table was made for: \
         the recipe here and in generate.py have drifted apart"
    );
    assert_eq!(rows.len(), inputs.len(), "{name}: one row per input");
    let strict = bound.is_some();
    let mut found = Found {
        worst: Vec::new(),
        incorrect: 0,
        points: rows.len(),
    };
    for (&x, &(hi_bits, frac)) in inputs.iter().zip(&rows) {
        let x = f64::from_bits(x);
        let hi = f64::from_bits(hi_bits);
        let y = f(x);
        if hi.is_infinite() {
            if strict {
                assert_eq!(y, hi, "{name}({x:e}) overflows to {hi}, and gave {y:e}");
            }
            continue;
        }
        if !y.is_finite() {
            assert!(!strict, "{name}({x:e}) = {y}, and exactly it is {hi:e}");
            continue;
        }
        if y.to_bits() != hi_bits {
            found.incorrect += 1;
        }
        // y − hi is exact: the two are within a factor of two of each other, or both subnormal.
        let error = ((y - hi) / ulp_of_exact(hi, frac) - frac).abs();
        let (region, limit) = bound.map_or(("all", f64::INFINITY), |b| b(x));
        match found.worst.iter_mut().find(|w| w.0 == region) {
            Some(w) if error > w.1 => (w.1, w.2) = (error, x),
            Some(_) => {}
            None => found.worst.push((region, error, x)),
        }
        assert!(
            error <= limit,
            "{name}({x:e} = {:#018x}) = {y:e}, {error:.4} ulp from the exact {hi:e}{frac:+.4}, \
             over the {limit} ulp bound for {region}",
            x.to_bits()
        );
    }
    for (region, worst, at) in &found.worst {
        println!(
            "{name}, {region}: worst {worst:.4} ulp at x = {at:e} ({:#018x})",
            at.to_bits()
        );
    }
    println!(
        "{name}: {} points, {} ({:.3}%) not correctly rounded",
        found.points,
        found.incorrect,
        100.0 * found.incorrect as f64 / found.points as f64,
    );
    found
}

/// **`exp`'s bound, 0.520 ulp everywhere, derived** from the method rather than read off the
/// measurement. `math::exp` computes `2ᵏ (T_hi + c)`, `c = T_hi p + T_lo`, `p ≈ eʳ − 1`, and
/// rounds the sum once. Before that rounding, absolutely and with `T_hi` in `[1, 2)`:
///
/// - `r` is off by at most `2⁻⁶²`: the head of the Cody–Waite product and its difference are
///   exact, the tail's product rounds by `2⁻⁸⁰`, `ln 2`'s split is short by `2⁻⁸⁰` over
///   `|n| < 2¹⁸`, and the last subtraction rounds by half an ulp of `|r| < 2⁻⁸`, `2⁻⁶²`;
/// - `p` is off from `eʳ − 1` by the truncation `r⁶/720 · (1 + r) < 2⁻⁶⁰·⁶⁶` at
///   `|r| ≤ ln 2/256 · (1 + 2⁻³⁴)`, the last addition's `2⁻⁶²`, and `2⁻⁶⁹·⁶` inside;
/// - `T_hi p` and `+ T_lo` round by at most `2⁻⁶¹` each, being under `2⁻⁷·⁵`; `T_lo p`, left out,
///   is under `2⁻⁶¹·⁵`; `T_hi + T_lo` is `2^(j/128)` to `2⁻¹⁰⁶`.
///
/// Relative to the result, at least `T_hi (1 − 2⁻⁸·⁵)`, that is at most `2.16 × 10⁻¹⁸`, largest at
/// `T_hi = 1`; an ulp is at least `2⁻⁵³` of the value, so it is **0.0195 ulp**, and the last
/// rounding's ½ makes 0.5195. Scaling by `2ᵏ` is exact. Below `2⁻¹⁰²²` the sum is rounded once
/// onto the subnormal grid instead, by adding it to 1 at `2¹⁰²²` times its size: the same error
/// carried in, now at most half as large against the grid, and `2⁻⁹` of the grid for the rounding
/// of the small part — under 0.512.
fn exp_bound(_: f64) -> (&'static str, f64) {
    ("every x", 0.520)
}

/// **`ln`'s bounds, by region of its method, derived.** `math::ln` writes `x = 2ᵏ m`, `m` in
/// `[0.703125, 1.40625)`, rounds `m` to `F` on a grid of `2⁻⁹` below 1 and `2⁻⁸` above, takes
/// `c ≈ 1/F` of nine bits from its table, and computes `r = m c − 1` **exactly**, `|r| ≤ 1.71 ·
/// 2⁻⁹`. It returns `s + (e + (lo + p))`: `s = h + r` with `h = k ln2_hi + (−ln c)_hi` exact
/// and `e` the sum's error, exact by Dekker; `lo = k ln2_lo + (−ln c)_lo`; `p ≈ ln(1 + r) − r`.
/// Every error is then `p`'s and three small roundings, against the last rounding's ½:
///
/// - `p = r² Q(r)` is evaluated to `4 · 2⁻⁵³` of itself — `r²`, the two sums of `Q`, whose
///   magnitude is about ½, and the product, `2⁻⁵³` each — and truncated after `r⁷`, `|r|⁸/8`.
/// - `lo + p` and `e + (lo + p)` are under `2⁻¹⁷` and round by `2⁻⁷¹` each; `(−ln c)_lo` is
///   represented to `2⁻⁸⁹`, and for `k ≠ 0`, `k ln2_lo`'s rounding and `ln2_lo`'s
///   representation are `2⁻⁸⁰` each at `|k| ≤ 1075`.
///
/// By region:
///
/// - **the cells about 1, `c = 1`**, `x` in `[1 − 3 · 2⁻¹⁰, 1 + 2⁻⁹)`: `h = e = lo = 0`, `r = m − 1`
///   and the result is `r + p` rounded once. Against `|ln x| ≥ |r| (1 − |r|/2)`, `p`'s errors are
///   at most `4 · 2⁻⁵³ · |r|/2 · 1.003 + |r|⁷/8` relatively, at `|r| ≤ 3 · 2⁻¹⁰`: **0.0080 ulp**,
///   so 0.508.
/// - **elsewhere in `k = 0`**: the absolute errors above, taken cell by cell at the cell's
///   largest `|r|` and smallest `|ln x|`. The worst is `F = 1 + 2⁻⁸`, `|ln x| ≥ ln(1 + 2⁻⁹)`:
///   **0.0080 ulp**, so 0.508.
/// - **`k ≠ 0`**: `|ln x| ≥ ln 1.40625 = 0.34`, and the same absolute errors are **0.0002 ulp**
///   of it: 0.501. A subnormal argument is scaled by `2⁵²`, exactly, first.
fn ln_bound(x: f64) -> (&'static str, f64) {
    if (1.0 - 3.0 * pow2(-10)..1.0 + pow2(-9)).contains(&x) {
        ("the cells about 1, c = 1", 0.508)
    } else if (0.703_125..1.406_25).contains(&x) {
        ("k = 0, c ≠ 1", 0.508)
    } else {
        ("k ≠ 0", 0.501)
    }
}

/// One step up from `f(x)`, in place of it: the error the comparison must see.
fn one_ulp_up(f: fn(f64) -> f64) -> impl Fn(f64) -> f64 {
    move |x| {
        let y = f(x);
        if y.is_finite() {
            f64::from_bits(if y >= 0.0 {
                y.to_bits() + 1
            } else {
                y.to_bits() - 1
            })
        } else {
            y
        }
    }
}

/// **`exp` is within its bound of mpmath at every point**: the whole domain, every reduction
/// boundary `(k + ½) ln 2` for `k` from −1075 to 1024, every boundary of the table index
/// `(n + ½) ln 2/128` for `|n| ≤ 512`, the branch thresholds each side, and the subnormal results.
///
/// **And the comparison can see one ulp**: the same values each moved up one double are more than
/// one ulp from the exact value somewhere, which a comparison that read the wrong row, or rounded
/// both sides alike, would not report.
#[test]
fn exp_is_within_its_bound_of_mpmath() {
    let (inputs, table) = (exp_inputs(1), include_bytes!("data/exp-reference.bin"));
    compare("exp", &exp, &inputs, table, Some(exp_bound));
    let shifted = compare(
        "exp moved up one ulp",
        &one_ulp_up(exp),
        &inputs,
        table,
        None,
    );
    assert!(shifted.worst() > 1.0, "one ulp up is still within one ulp");
}

/// **`ln` is within its bounds of mpmath at every point**: every positive double, every subnormal,
/// densely about 1, every binade at √2 and at 0.703125 (where `math`'s `k` steps), every edge of
/// its 256 cells at `k = −1, 0, 1`, the descreening integral's range, and every power of two. And
/// the comparison can see one ulp, as for `exp`.
#[test]
fn ln_is_within_its_bounds_of_mpmath() {
    let (inputs, table) = (ln_inputs(1), include_bytes!("data/ln-reference.bin"));
    compare("ln", &ln, &inputs, table, Some(ln_bound));
    let shifted = compare("ln moved up one ulp", &one_ulp_up(ln), &inputs, table, None);
    assert!(shifted.worst() > 1.0, "one ulp up is still within one ulp");
}

/// The same against a larger table, written outside the repository:
///
/// ```text
/// python tools/math-reference/generate.py --scale 10 --out <dir>
/// PANTOMETRY_MATH_REFERENCE=<dir> PANTOMETRY_MATH_REFERENCE_SCALE=10 \
///     cargo test --release -p pantometry-core --test exp_and_ln_against_mpmath -- --ignored --nocapture
/// ```
///
/// Ignored, and not run in CI, because the table is not committed.
#[cfg(not(target_family = "wasm"))]
#[test]
#[ignore = "needs a table written by generate.py --scale: see its documentation"]
fn against_a_larger_table() {
    let dir = std::env::var("PANTOMETRY_MATH_REFERENCE").expect("PANTOMETRY_MATH_REFERENCE");
    let scale: u64 = std::env::var("PANTOMETRY_MATH_REFERENCE_SCALE")
        .expect("PANTOMETRY_MATH_REFERENCE_SCALE")
        .parse()
        .unwrap();
    let read = |name: &str| std::fs::read(format!("{dir}/{name}-reference.bin")).unwrap();
    compare(
        "exp",
        &exp,
        &exp_inputs(scale),
        &read("exp"),
        Some(exp_bound),
    );
    compare("ln", &ln, &ln_inputs(scale), &read("ln"), Some(ln_bound));
}

/// **The cost**, in release: nanoseconds per call over arguments where generalized Born uses
/// them -- `exp` on `[-50, 0]`, Still's `-r^2/4R_iR_j`, and `ln` on `[1, 64]`, the descreening
/// integral's `U/L` -- beside the platform's `f64::exp` and `f64::ln`. Asserted about the cost:
/// nothing, since it is a property of the machine; printed.
#[cfg(not(target_family = "wasm"))]
#[test]
#[ignore = "a timing: run with --release -- --ignored --nocapture"]
fn the_cost() {
    use std::hint::black_box;
    use std::time::Instant;
    let n = 4096u64;
    let u = |seed: u64| -> Vec<f64> {
        (0..n)
            .map(|i| (splitmix(seed, i) >> 11) as f64 / (1u64 << 53) as f64)
            .collect()
    };
    let exp_args: Vec<f64> = u(41).iter().map(|v| -50.0 * v).collect();
    let ln_args: Vec<f64> = u(42).iter().map(|v| 1.0 + 63.0 * v).collect();
    let time = |f: &dyn Fn(f64) -> f64, args: &[f64]| -> f64 {
        let rounds = 2000;
        let mut best = f64::INFINITY;
        for _ in 0..5 {
            let t = Instant::now();
            let mut sum = 0.0;
            for _ in 0..rounds {
                for &x in black_box(args) {
                    sum += f(x);
                }
            }
            black_box(sum);
            best = best.min(t.elapsed().as_secs_f64() * 1e9 / (rounds as f64 * args.len() as f64));
        }
        best
    };
    let ours_exp = time(&|x| exp(x), &exp_args);
    let platform_exp = time(&|x: f64| x.exp(), &exp_args);
    let ours_ln = time(&|x| ln(x), &ln_args);
    let platform_ln = time(&|x: f64| x.ln(), &ln_args);
    println!(
        "ns per call: exp {ours_exp:.2} (platform {platform_exp:.2}), ln {ours_ln:.2} (platform {platform_ln:.2})"
    );
}

/// `2^e` for a normal `e`, exactly, from its bits. `f64::powi`'s precision is unspecified by Rust,
/// and on two CI runners it gave 0 for `2^-1074`, so nothing here relies on it for an input.
fn pow2(e: i32) -> f64 {
    assert!((-1022..=1023).contains(&e), "2^{e} is not a normal double");
    f64::from_bits(((e + 1023) as u64) << 52)
}
