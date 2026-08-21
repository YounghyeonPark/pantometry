//! The colour scales, constructed in CIE LCh rather than picked.
//!
//! # Why this is not a gradient someone liked
//!
//! A colour scale over a scalar field makes two promises, and both are measurable:
//!
//! 1. **A larger value is never darker.** Otherwise the scale folds back on itself and two
//!    different values are the same colour to anyone reading it in greyscale, in a print, or
//!    with a colour-vision deficiency.
//! 2. **Equal steps in the value look like equal steps.** Otherwise the same physical change
//!    looks large in one part of the range and small in another, which is the colour axis of
//!    exactly the sin [`crate::report`] refuses in the time axis when it fixes the scale across
//!    a run.
//!
//! Both are statements about CIE L\* and about the distance between colours in CIELAB, so both
//! are things a test can fail. The scale this crate shipped until now was a four-stop linear
//! interpolation in sRGB — dark blue, blue, amber, orange-red — and it kept neither. Measured
//! over 256 steps:
//!
//! | | the four-stop ramp | [`sequential`] |
//! | --- | --- | --- |
//! | steps where a larger value is **darker** | **89 of 255** | 0 of 255 |
//! | where the scale is brightest | t = 0.67, then it falls back | t = 1 |
//! | L\* covered | 22.9 → 82.0, folding | 18 → 92, monotone |
//! | largest single step, ΔE₇₆ | 2.34 | 1.26 |
//! | how much a 10% change of value varies in apparent size | **3.07×** | 1.85× |
//!
//! The 89 backwards steps are the amber-to-red arm: the top third of that scale is darker than
//! its middle, so a value at 0.55 and a value at 1.0 print as the same grey.
//!
//! # How these are built
//!
//! Lightness is **linear in the value**, by construction — that is promise 1, made structurally
//! rather than checked. Hue sweeps once, monotonically. Chroma follows a smooth hump, and is
//! then reduced until the colour exists in sRGB at all.
//!
//! That last reduction is the only approximation here, and it is where the residual 1.85× comes
//! from: along this hue path sRGB holds a chroma of about 27 in the blues and about 100 in the
//! greens, so mid-scale steps carry less colour difference than the ends however they are asked
//! for. It is a property of the sRGB solid and not of the construction, and no map that stays in
//! sRGB and keeps L\* linear avoids it.
//!
//! # Which map, and chosen from the data
//!
//! [`sequential`] for a field with a floor — a temperature, a speed, an intensity. [`diverging`]
//! for a field that swings about zero — a pressure, a component of **E**, a displacement — where
//! the neutral point is a *number* and not a place in the range, and where +x and −x must look
//! equally far from it. `report` picks between them by whether the run's range straddles zero,
//! which is the same test [`crate::report`] already used to choose a volume's opacity curve. The
//! opacity knew about the sign before the colour did.

/// How many entries a rendered table has. Eight bits per channel is what a canvas holds, so a
/// finer table would quantise to the same colours.
pub const STEPS: usize = 256;

// --- sRGB <-> CIELAB, IEC 61966-2-1 and CIE 15, D65 2-degree ------------------------------------

/// The D65 white point, as CIE XYZ with Y = 1 — **taken from the matrix below rather than from
/// the standard's table**.
///
/// The two disagree in the seventh place, because the sRGB primaries matrix is published rounded
/// to seven figures and its rows therefore sum to D65 only to that. Using the tabulated white
/// against the rounded matrix put sRGB white at L\* 100.000004, which is not wrong by anything
/// that matters and is wrong by more than zero — and the one property this file is for is that
/// its claims are exact where they can be. Row sums make the round trip exact by construction.
const WHITE: [f64; 3] = [
    RGB_TO_XYZ[0][0] + RGB_TO_XYZ[0][1] + RGB_TO_XYZ[0][2],
    RGB_TO_XYZ[1][0] + RGB_TO_XYZ[1][1] + RGB_TO_XYZ[1][2],
    RGB_TO_XYZ[2][0] + RGB_TO_XYZ[2][1] + RGB_TO_XYZ[2][2],
];

/// Linear-light RGB from CIE XYZ. The inverse of the sRGB primaries matrix.
const XYZ_TO_RGB: [[f64; 3]; 3] = [
    [3.240_454_2, -1.537_138_5, -0.498_531_4],
    [-0.969_266_0, 1.876_010_8, 0.041_556_0],
    [0.055_643_4, -0.204_025_9, 1.057_225_2],
];

/// CIE XYZ from linear-light RGB. The sRGB primaries under D65.
const RGB_TO_XYZ: [[f64; 3]; 3] = [
    [0.412_456_4, 0.357_576_1, 0.180_437_5],
    [0.212_672_9, 0.715_152_2, 0.072_175_0],
    [0.019_333_9, 0.119_192_0, 0.950_304_1],
];

fn lab_to_linear_rgb(l: f64, a: f64, b: f64) -> [f64; 3] {
    let fy = (l + 16.0) / 116.0;
    let fx = fy + a / 500.0;
    let fz = fy - b / 200.0;
    // The CIE piecewise inverse, with the linear segment below the cube-root's knee.
    let g = |t: f64| {
        let c = t * t * t;
        if c > 216.0 / 24389.0 {
            c
        } else {
            (116.0 * t - 16.0) * 108.0 / 841.0
        }
    };
    let (x, y, z) = (WHITE[0] * g(fx), WHITE[1] * g(fy), WHITE[2] * g(fz));
    [
        XYZ_TO_RGB[0][0] * x + XYZ_TO_RGB[0][1] * y + XYZ_TO_RGB[0][2] * z,
        XYZ_TO_RGB[1][0] * x + XYZ_TO_RGB[1][1] * y + XYZ_TO_RGB[1][2] * z,
        XYZ_TO_RGB[2][0] * x + XYZ_TO_RGB[2][1] * y + XYZ_TO_RGB[2][2] * z,
    ]
}

fn in_gamut(rgb: [f64; 3]) -> bool {
    rgb.iter().all(|c| (-1e-9..=1.0 + 1e-9).contains(c))
}

fn encode(c: f64) -> u8 {
    let c = c.clamp(0.0, 1.0);
    let s = if c <= 0.003_130_8 {
        12.92 * c
    } else {
        1.055 * c.powf(1.0 / 2.4) - 0.055
    };
    (s * 255.0).round().clamp(0.0, 255.0) as u8
}

fn decode(c: u8) -> f64 {
    let c = c as f64 / 255.0;
    if c <= 0.040_45 {
        c / 12.92
    } else {
        ((c + 0.055) / 1.055).powf(2.4)
    }
}

/// The CIELAB coordinates of an 8-bit sRGB triple.
///
/// Public because the properties this module claims are claims about *the colours that come out*,
/// not about the parameters that went in, and a caller checking them should not have to
/// reimplement the transform. [`lightness`] is the common case.
pub fn to_lab(rgb: [u8; 3]) -> [f64; 3] {
    let (r, g, b) = (decode(rgb[0]), decode(rgb[1]), decode(rgb[2]));
    let x = (RGB_TO_XYZ[0][0] * r + RGB_TO_XYZ[0][1] * g + RGB_TO_XYZ[0][2] * b) / WHITE[0];
    let y = (RGB_TO_XYZ[1][0] * r + RGB_TO_XYZ[1][1] * g + RGB_TO_XYZ[1][2] * b) / WHITE[1];
    let z = (RGB_TO_XYZ[2][0] * r + RGB_TO_XYZ[2][1] * g + RGB_TO_XYZ[2][2] * b) / WHITE[2];
    let f = |t: f64| {
        if t > 216.0 / 24389.0 {
            t.cbrt()
        } else {
            (841.0 / 108.0) * t + 4.0 / 29.0
        }
    };
    let (fx, fy, fz) = (f(x), f(y), f(z));
    [116.0 * fy - 16.0, 500.0 * (fx - fy), 200.0 * (fy - fz)]
}

/// CIE L\* of an 8-bit sRGB triple — how light it is, on the scale a person perceives.
///
/// 0 is black and 100 is the white point. This is the number the "a larger value is never
/// darker" promise is about, and the number a greyscale print keeps.
pub fn lightness(rgb: [u8; 3]) -> f64 {
    to_lab(rgb)[0]
}

/// CIE76 colour difference — the Euclidean distance in CIELAB.
///
/// Used here rather than ΔE₀₀ because the claim being checked is about the *spread* of step
/// sizes across one scale, and ΔE₀₀'s corrections change the individual numbers without
/// reordering them. A single ΔE₇₆ of about 1 is roughly the smallest difference a person can see
/// side by side.
pub fn difference(a: [u8; 3], b: [u8; 3]) -> f64 {
    let (p, q) = (to_lab(a), to_lab(b));
    ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
}

/// The sRGB triple for a lightness, a chroma and a hue, with the chroma reduced until it fits.
fn lch(l: f64, chroma: f64, hue_deg: f64) -> [u8; 3] {
    let hr = hue_deg.to_radians();
    let at = |c: f64| lab_to_linear_rgb(l, c * hr.cos(), c * hr.sin());
    let rgb = at(chroma.min(max_chroma(l, hue_deg)));
    [encode(rgb[0]), encode(rgb[1]), encode(rgb[2])]
}

/// The largest chroma sRGB holds at this lightness and hue.
///
/// Bisected rather than solved: the sRGB solid in LCh has no closed boundary, and 40 halvings
/// place it to better than a thousandth of a chroma unit, which is far below what eight bits per
/// channel can express anyway.
fn max_chroma(l: f64, hue_deg: f64) -> f64 {
    let hr = hue_deg.to_radians();
    let at = |c: f64| lab_to_linear_rgb(l, c * hr.cos(), c * hr.sin());
    if in_gamut(at(MAX_CHROMA_SEARCH)) {
        return MAX_CHROMA_SEARCH;
    }
    let (mut lo, mut hi) = (0.0, MAX_CHROMA_SEARCH);
    for _ in 0..40 {
        let mid = 0.5 * (lo + hi);
        if in_gamut(at(mid)) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

/// Where the bisection starts. No sRGB colour exceeds a chroma of about 132 (blue, at L\* 32),
/// so 150 is above the solid everywhere and the search always brackets.
const MAX_CHROMA_SEARCH: f64 = 150.0;

// --- the two scales -----------------------------------------------------------------------------

/// L\* at the bottom and top of [`sequential`]. 18 rather than 0 so the low end still separates
/// from the report's panel background, which sits at L\* ≈ 4; 92 rather than 100 because the
/// sRGB solid closes up above it and the chroma would collapse over the last two per cent.
const SEQ_L: (f64, f64) = (18.0, 92.0);
/// Hue in degrees, swept once: violet-blue, through blue, teal and green, to yellow.
const SEQ_H: (f64, f64) = (300.0, 95.0);
/// The chroma asked for, before the gamut has its say.
const SEQ_C: f64 = 55.0;

/// A scale for a field with a floor — a temperature, a speed, an intensity.
///
/// `t` outside `[0, 1]` is clamped, so a caller does not have to. Lightness rises linearly from
/// L\* 18 to L\* 92 across the range: **the top is always the lightest colour on the scale**, and
/// a greyscale print of it is still a scale.
pub fn sequential(t: f64) -> [u8; 3] {
    let t = t.clamp(0.0, 1.0);
    let l = SEQ_L.0 + (SEQ_L.1 - SEQ_L.0) * t;
    let h = SEQ_H.0 + (SEQ_H.1 - SEQ_H.0) * t;
    // A hump: least chroma at the two ends, where the gamut is narrowest and a saturated colour
    // would be clipped to whatever the boundary happens to be, and most in the middle.
    let c = SEQ_C * (0.55 + 0.45 * (core::f64::consts::PI * (0.12 + 0.76 * t)).sin());
    lch(l, c, h)
}

/// L\* at the neutral point and at full deflection. **Dark in the middle**, which inverts the
/// usual convention and does so for two reasons.
///
/// The report is dark, and the neutral point should be the quietest thing on it: a cavity at
/// rest, or a room before the pulse arrives, is a panel of zeros, and a light neutral makes that
/// panel glow. The convention is light-in-the-middle because diverging scales were designed for
/// paper. The second reason is that the volume renderer already decided this — its opacity curve
/// is `|2t − 1|`, transparent at the neutral and opaque at both extremes — and a colour scale
/// that brightened where the opacity vanishes would be arguing with it.
///
/// The cost is the one every diverging scale pays: lightness carries *how far from zero* and
/// hue carries *which way*, so a greyscale print of a signed field shows magnitude and loses the
/// sign. There is no scale that keeps both.
const DIV_L: (f64, f64) = (20.0, 68.0);
/// Hue for the low arm and for the high arm. They do not blend: they meet at a neutral.
///
/// Blue against orange rather than the conventional blue against red. Red and blue at equal
/// lightness are the pair that deuteranopia and protanopia collapse — which is most of the eight
/// per cent of men with a colour-vision deficiency — and this scale puts its two arms at equal
/// lightness *by construction*, so it would collapse them completely. The blue-yellow axis
/// survives both.
const DIV_H: (f64, f64) = (258.0, 68.0);
/// Chroma at full deflection.
const DIV_C: f64 = 47.0;

/// A scale for a field that swings about zero — a pressure, a component of **E**, a displacement.
///
/// `t = 0.5` is the neutral point and must be **where zero falls in the range**, which is the
/// caller's job: a range of −100 to +300 puts zero at 0.25, and passing the value's position in
/// the range without re-centring would colour zero as if it were negative.
///
/// Lightness is symmetric about the middle by construction, so a value and its negative are
/// equally far from neutral to look at. Measured over the scale, the largest gap between the two
/// arms is 0.32 ΔE₇₆, which is 8-bit quantisation and not a bias.
pub fn diverging(t: f64) -> [u8; 3] {
    let t = t.clamp(0.0, 1.0);
    let s = (2.0 * t - 1.0).abs();
    let l = DIV_L.0 + (DIV_L.1 - DIV_L.0) * s;
    let h = if t < 0.5 { DIV_H.0 } else { DIV_H.1 };
    // **The chroma both arms can hold, not the chroma each one can.** sRGB is not symmetric in
    // hue, and clipping each arm to its own ceiling is how a diverging scale acquires a bias with
    // no counterpart in the data: an earlier set of parameters here put a value and its negative
    // 65.0 and 75.6 from neutral, a 16% asymmetry that came entirely from the gamut.
    //
    // With the constants above nothing is clipped at all -- 47 chroma at L* 68 fits inside both
    // arms, and the two spellings produce byte-identical tables, measured. This is the guard that
    // keeps that true if the constants move, and it is cheap: the bisection runs 256 times to
    // build a table that is then used for every pixel of every frame.
    let ceiling = max_chroma(l, DIV_H.0).min(max_chroma(l, DIV_H.1));
    lch(l, (DIV_C * s).min(ceiling), h)
}

/// Which scale a run wants, from the run's own range.
///
/// The same test the volume renderer already used to choose an opacity curve: a range that
/// straddles zero is a signed quantity and gets [`diverging`]; anything else gets
/// [`sequential`]. Nothing here reads the domain's name or its unit — a view chosen by the name
/// of the physics is a view that a new domain does not get.
pub fn is_signed(lo: f64, hi: f64) -> bool {
    lo < 0.0 && hi > 0.0
}

/// [`STEPS`] entries of a scale, as lower-case hex with no separators — `rrggbb` repeated.
///
/// This is how the scale reaches the HTML report: the table is built here, so the picture in the
/// browser and the properties this module's tests pin are the same numbers rather than two
/// implementations that agree until one is edited. 1,536 characters a scale.
pub fn hex_table(map: fn(f64) -> [u8; 3]) -> String {
    let mut out = String::with_capacity(STEPS * 6);
    for i in 0..STEPS {
        let c = map(i as f64 / (STEPS - 1) as f64);
        out.push_str(&format!("{:02x}{:02x}{:02x}", c[0], c[1], c[2]));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Round-tripping a known colour through both transforms, against a **published** value:
    /// sRGB white is L\* = 100 exactly, and the mid grey `#808080` is L\* ≈ 53.585.
    ///
    /// This is the check that the matrices and the two piecewise curves are right. Everything
    /// else in this module is measured with [`to_lab`], so a transform that was wrong would make
    /// every other test here agree with itself and with nothing outside.
    #[test]
    fn the_lab_transform_agrees_with_published_values() {
        let white = lightness([255, 255, 255]);
        assert!((white - 100.0).abs() < 1e-6, "white is L* {white}");
        assert!(lightness([0, 0, 0]).abs() < 1e-9);
        // #808080 -> Y = 0.2158605, L* = 116*(Y)^(1/3) - 16 = 53.5850.
        let grey = lightness([128, 128, 128]);
        assert!((grey - 53.585).abs() < 1e-3, "mid grey is L* {grey}");
        // A saturated primary, chosen because it exercises the off-diagonal matrix terms:
        // sRGB red is L* 53.2408, a* 80.0925, b* 67.2032.
        let red = to_lab([255, 0, 0]);
        assert!(
            (red[0] - 53.2408).abs() < 1e-3
                && (red[1] - 80.0925).abs() < 1e-2
                && (red[2] - 67.2032).abs() < 1e-2,
            "sRGB red is {red:?}"
        );
    }

    /// **Promise 1**, and the one the old scale broke 89 times: over the whole table, lightness
    /// never falls as the value rises.
    ///
    /// The tolerance is one 8-bit quantisation step, not zero. L\* is computed from the rounded
    /// triple, and two adjacent entries 0.3 L\* apart can round to the same byte and then differ
    /// by a few thousandths in the wrong direction. Measured across the table the largest such
    /// wobble is under 0.05 L\*; 0.1 is that with room, and is still forty times smaller than the
    /// 2.3 L\* a single visible step carries.
    #[test]
    fn the_sequential_scale_never_gets_darker() {
        let mut worst = 0.0f64;
        for i in 0..STEPS - 1 {
            let a = lightness(sequential(i as f64 / (STEPS - 1) as f64));
            let b = lightness(sequential((i + 1) as f64 / (STEPS - 1) as f64));
            worst = worst.max(a - b);
        }
        assert!(worst < 0.1, "lightness falls by {worst} L* somewhere");
        // And it covers the range it claims, which is what makes it a scale and not a tint.
        let (lo, hi) = (lightness(sequential(0.0)), lightness(sequential(1.0)));
        assert!(
            (lo - 18.0).abs() < 0.6 && (hi - 92.0).abs() < 0.6,
            "L* runs {lo} to {hi}"
        );
    }

    /// **Promise 2**: the same change of value looks about the same size wherever it happens.
    ///
    /// Measured as the apparent distance across a window of 10% of the range, swept over the
    /// scale. The bound is 2.0× and the measurement is 1.85×, and the gap is deliberately small:
    /// the residual is the sRGB solid — about 27 chroma available in the blues against about 100
    /// in the greens — so a bound of 5× would sit still through a change that broke the
    /// construction. The four-stop scale this replaced measured 3.07×.
    #[test]
    fn a_fixed_change_of_value_looks_the_same_size_wherever_it_is() {
        let window = 0.10;
        let mut lo = f64::MAX;
        let mut hi: f64 = 0.0;
        for i in 0..=36 {
            let t = i as f64 / 40.0;
            let d = difference(sequential(t), sequential(t + window));
            lo = lo.min(d);
            hi = hi.max(d);
        }
        let spread = hi / lo;
        assert!(
            spread < 2.0,
            "a 10% window varies {spread:.2}x ({lo:.1}..{hi:.1})"
        );
    }

    /// No step is an edge that is not in the data.
    ///
    /// A scale assembled from stops has a corner at every stop, and a corner reads as a contour
    /// line in a field that has none. The bound is 1.5 ΔE₇₆ against a measured 1.26 and a mean of
    /// about 0.74 — so the largest step is under 1.7 times the average one.
    #[test]
    fn no_single_step_stands_out() {
        let table: Vec<[u8; 3]> = (0..STEPS)
            .map(|i| sequential(i as f64 / (STEPS - 1) as f64))
            .collect();
        let steps: Vec<f64> = (0..STEPS - 1)
            .map(|i| difference(table[i], table[i + 1]))
            .collect();
        let worst = steps.iter().cloned().fold(0.0f64, f64::max);
        let mean = steps.iter().sum::<f64>() / steps.len() as f64;
        assert!(worst < 1.5, "largest step {worst:.2} dE, mean {mean:.2}");
    }

    /// A value and its negative must be equally far from neutral, or a picture of a standing wave
    /// shows one half of it more strongly than the other and the asymmetry is the renderer's.
    ///
    /// Checked on lightness and on total colour difference from the neutral point. The tolerance
    /// is 8-bit quantisation: 0.3 L\* is one step of the low byte at this lightness.
    #[test]
    fn the_diverging_scale_is_symmetric_about_its_middle() {
        let mid = diverging(0.5);
        for k in 1..=10 {
            let d = k as f64 / 20.0;
            let (below, above) = (diverging(0.5 - d), diverging(0.5 + d));
            let (lb, la) = (lightness(below), lightness(above));
            assert!(
                (lb - la).abs() < 0.3,
                "at +/-{d} the arms are L* {lb:.2} and {la:.2}"
            );
            // Equal lightness and equal chroma, so the distance from neutral is equal too --
            // exactly, up to the rounding to eight bits. The bound is 1 dE, which is about the
            // smallest difference a person can see side by side and about three quantisation
            // steps; the measurement is under 0.4.
            let (db, da) = (difference(mid, below), difference(mid, above));
            assert!(
                (db - da).abs() < 1.0,
                "at +/-{d} the arms are {db:.2} and {da:.2} from neutral"
            );
        }
        // The neutral point is neutral: nearly no chroma, so zero does not read as either sign.
        let lab = to_lab(mid);
        assert!(
            lab[1].hypot(lab[2]) < 2.0,
            "the middle carries chroma {:.1}",
            lab[1].hypot(lab[2])
        );
    }

    /// The table that reaches the browser is the table this module's other tests measured.
    ///
    /// A separate implementation in JavaScript is the failure this avoids: it would agree on the
    /// day it was written and drift on the first edit, and nothing here would notice because
    /// every test above would still be measuring the Rust.
    #[test]
    fn the_hex_table_is_the_scale() {
        let table = hex_table(sequential);
        assert_eq!(table.len(), STEPS * 6);
        assert!(table
            .chars()
            .all(|c| c.is_ascii_hexdigit() && !c.is_uppercase()));
        for i in [0usize, 1, 128, STEPS - 1] {
            let c = sequential(i as f64 / (STEPS - 1) as f64);
            let want = format!("{:02x}{:02x}{:02x}", c[0], c[1], c[2]);
            assert_eq!(&table[i * 6..i * 6 + 6], want, "entry {i}");
        }
    }

    /// The choice of scale is made by the numbers and not by the name of the domain.
    #[test]
    fn a_range_that_straddles_zero_is_signed_and_nothing_else_is() {
        assert!(is_signed(-1.0, 1.0));
        assert!(is_signed(-1e-30, 4.0));
        assert!(
            !is_signed(0.0, 400.0),
            "a floor at zero is not a swing about it"
        );
        assert!(!is_signed(293.0, 353.0));
        assert!(
            !is_signed(-9.0, -1.0),
            "all-negative has a ceiling, not a centre"
        );
    }
}
