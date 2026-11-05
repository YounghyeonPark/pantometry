//! What the tests against mpmath share: `generate.py`'s recipe for drawing points, the table
//! format, and the comparison. Each test file uses part of it.
#![allow(dead_code)]

pub const MASK: u64 = u64::MAX;

pub fn splitmix(seed: u64, i: u64) -> u64 {
    let mut z = seed.wrapping_add((i + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// `n` bit patterns uniform in `[lo, hi]`; `sign` 0 or 1 is the sign bit, 2 alternates.
pub fn uniform(lo: u64, hi: u64, seed: u64, n: u64, sign: u64) -> Vec<u64> {
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

pub fn around(base: u64, reach: i64, signs: &[u64]) -> Vec<u64> {
    let mut out = Vec::new();
    for &s in signs {
        for j in -reach..=reach {
            out.push((base.wrapping_add(j as u64) & MASK) | (s << 63));
        }
    }
    out
}

pub fn b(x: f64) -> u64 {
    x.to_bits()
}

pub fn fnv1a(words: &[u64]) -> u64 {
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
pub fn parse(bytes: &[u8]) -> (u64, Vec<(u64, f64)>) {
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
pub fn ulp_of_exact(hi: f64, frac: f64) -> f64 {
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
pub type Bound = fn(f64) -> (&'static str, f64);

/// What a comparison found: per region, its worst error and where; and how many were not
/// correctly rounded.
pub struct Found {
    pub worst: Vec<(&'static str, f64, f64)>,
    pub incorrect: usize,
    pub points: usize,
}

impl Found {
    pub fn worst(&self) -> f64 {
        self.worst.iter().map(|w| w.1).fold(0.0, f64::max)
    }
}

/// `f` at every input against the table. With a `bound`, each point is held to its region's bound;
/// without one the comparison only measures, for the check that it can see one ulp.
pub fn compare(
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
        // By value rather than by bits: mpmath has no −0, so the table's zeros are all +0, and
        // `sin_cos_turns` signs its exact zeros as IEEE 754's sinPi does.
        if y != hi {
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

/// One step up from `f(x)`, in place of it: the error the comparison must see.
pub fn one_ulp_up(f: impl Fn(f64) -> f64) -> impl Fn(f64) -> f64 {
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

/// **A row says its value is exact exactly when it is.** `frac` is stored in steps of `2⁻¹⁵`
/// ulp, and the spacing [`ulp_of_exact`] rebuilds is read from its sign: just below a power of
/// two it is half the rounded value's. A `frac` rounded to 0 there would make a whole ulp of
/// error measure as half — a value of `1 − 2⁻⁵³` where 1 was due passed sin's comparison at 2048
/// committed points, when `generate.py` still rounded such a `frac` to 0. So a row with `frac`
/// 0 must be a point where the function is exact, and a point where it is exact must have a row
/// with `frac` 0; `is_exact` says which points those are, from the mathematics rather than the
/// table. Rows that overflow are left out: their `frac` is 0 by construction.
pub fn table_marks_exact_values(
    name: &str,
    inputs: &[u64],
    table: &[u8],
    is_exact: fn(f64) -> bool,
) {
    let (digest, rows) = parse(table);
    assert_eq!(fnv1a(inputs), digest, "{name}: not the table's inputs");
    let mut exact = 0usize;
    for (&x, &(hi_bits, frac)) in inputs.iter().zip(&rows) {
        let x = f64::from_bits(x);
        if f64::from_bits(hi_bits).is_infinite() {
            continue;
        }
        assert_eq!(
            frac == 0.0,
            is_exact(x),
            "{name}({x:e} = {:#018x}): the row's frac is {frac}",
            x.to_bits()
        );
        exact += usize::from(frac == 0.0);
    }
    println!("{name}: {exact} rows exact, every other row on a side");
}

/// `2^e` for a normal `e`, exactly, from its bits. `f64::powi`'s precision is unspecified by Rust,
/// and on two CI runners it gave 0 for `2^-1074`, so nothing here relies on it for an input.
pub fn pow2(e: i32) -> f64 {
    assert!((-1022..=1023).contains(&e), "2^{e} is not a normal double");
    f64::from_bits(((e + 1023) as u64) << 52)
}
