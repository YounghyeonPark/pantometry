//! Numbers and units drawn as strokes, so a picture can be read without its caption.
//!
//! **A colour bar with no numbers on it is a gradient.** The viewer's snapshot is the only static
//! figure this workspace produces — the editor has its own painter and the HTML report has a
//! browser — and it had no legend, no scale and no units at all: a reader could see where the hot
//! end was and not whether it was 119 °C or 1190. That is the difference between a picture and a
//! reading, and it is the same sentence `colour_bar` in the editor was written under.
//!
//! # Why strokes rather than a font
//!
//! The viewer is a wgpu surface with a line pipeline and a triangle pipeline, and nothing else. A
//! real font means an atlas, a shaper, a licence and a file to ship; what is needed here is
//! sixteen digits and a handful of letters at one size, and a table of line segments is the whole
//! of it. It is also the only form that survives being drawn by the pipeline that is already
//! there.
//!
//! # The grid
//!
//! Each glyph is strokes on a 4-wide by 6-tall lattice, origin bottom left, and is drawn into a
//! caller-chosen box. Only what the labels need is here — the digits, a sign, a point, an
//! exponent, and the letters that spell the units this workspace reports. A character with no
//! entry draws nothing rather than a box, because a missing glyph in a legend should cost a
//! reader a character and not their confidence in the number beside it.

/// A stroke, as two points on the glyph lattice.
type Stroke = (f32, f32, f32, f32);

/// The strokes for one character, or an empty slice.
///
/// Lower case is folded to upper: the units are `K`, `C`, `Pa`, `V/m`, `m/s`, `mm`, `s`, `J`, and a
/// legend in one case reads as a legend rather than as prose.
fn strokes(c: char) -> &'static [Stroke] {
    match c.to_ascii_uppercase() {
        '0' => &[
            (0.0, 0.0, 4.0, 0.0),
            (4.0, 0.0, 4.0, 6.0),
            (4.0, 6.0, 0.0, 6.0),
            (0.0, 6.0, 0.0, 0.0),
            (0.0, 0.0, 4.0, 6.0),
        ],
        '1' => &[
            (1.0, 4.0, 2.0, 6.0),
            (2.0, 6.0, 2.0, 0.0),
            (0.0, 0.0, 4.0, 0.0),
        ],
        '2' => &[
            (0.0, 6.0, 4.0, 6.0),
            (4.0, 6.0, 4.0, 3.0),
            (4.0, 3.0, 0.0, 3.0),
            (0.0, 3.0, 0.0, 0.0),
            (0.0, 0.0, 4.0, 0.0),
        ],
        '3' => &[
            (0.0, 6.0, 4.0, 6.0),
            (4.0, 6.0, 4.0, 0.0),
            (4.0, 0.0, 0.0, 0.0),
            (0.0, 3.0, 4.0, 3.0),
        ],
        '4' => &[
            (0.0, 6.0, 0.0, 3.0),
            (0.0, 3.0, 4.0, 3.0),
            (4.0, 6.0, 4.0, 0.0),
        ],
        '5' => &[
            (4.0, 6.0, 0.0, 6.0),
            (0.0, 6.0, 0.0, 3.0),
            (0.0, 3.0, 4.0, 3.0),
            (4.0, 3.0, 4.0, 0.0),
            (4.0, 0.0, 0.0, 0.0),
        ],
        '6' => &[
            (4.0, 6.0, 0.0, 6.0),
            (0.0, 6.0, 0.0, 0.0),
            (0.0, 0.0, 4.0, 0.0),
            (4.0, 0.0, 4.0, 3.0),
            (4.0, 3.0, 0.0, 3.0),
        ],
        '7' => &[(0.0, 6.0, 4.0, 6.0), (4.0, 6.0, 1.0, 0.0)],
        '8' => &[
            (0.0, 0.0, 4.0, 0.0),
            (4.0, 0.0, 4.0, 6.0),
            (4.0, 6.0, 0.0, 6.0),
            (0.0, 6.0, 0.0, 0.0),
            (0.0, 3.0, 4.0, 3.0),
        ],
        '9' => &[
            (4.0, 0.0, 4.0, 6.0),
            (4.0, 6.0, 0.0, 6.0),
            (0.0, 6.0, 0.0, 3.0),
            (0.0, 3.0, 4.0, 3.0),
        ],
        '.' => &[(1.5, 0.0, 2.5, 0.0)],
        '-' => &[(0.5, 3.0, 3.5, 3.0)],
        '+' => &[(0.5, 3.0, 3.5, 3.0), (2.0, 1.5, 2.0, 4.5)],
        '/' => &[(0.0, 0.0, 4.0, 6.0)],
        // The degree ring, small and high, for a celsius legend.
        '\u{b0}' => &[
            (1.0, 5.0, 3.0, 5.0),
            (3.0, 5.0, 3.0, 6.0),
            (3.0, 6.0, 1.0, 6.0),
            (1.0, 6.0, 1.0, 5.0),
        ],
        'A' => &[
            (0.0, 0.0, 2.0, 6.0),
            (2.0, 6.0, 4.0, 0.0),
            (0.7, 2.0, 3.3, 2.0),
        ],
        'C' => &[
            (4.0, 6.0, 0.0, 6.0),
            (0.0, 6.0, 0.0, 0.0),
            (0.0, 0.0, 4.0, 0.0),
        ],
        'E' => &[
            (4.0, 6.0, 0.0, 6.0),
            (0.0, 6.0, 0.0, 0.0),
            (0.0, 0.0, 4.0, 0.0),
            (0.0, 3.0, 3.0, 3.0),
        ],
        'J' => &[
            (4.0, 6.0, 4.0, 1.0),
            (4.0, 1.0, 2.0, 0.0),
            (2.0, 0.0, 0.0, 1.0),
        ],
        'K' => &[
            (0.0, 6.0, 0.0, 0.0),
            (4.0, 6.0, 0.0, 3.0),
            (0.0, 3.0, 4.0, 0.0),
        ],
        'M' => &[
            (0.0, 0.0, 0.0, 6.0),
            (0.0, 6.0, 2.0, 3.0),
            (2.0, 3.0, 4.0, 6.0),
            (4.0, 6.0, 4.0, 0.0),
        ],
        'P' => &[
            (0.0, 0.0, 0.0, 6.0),
            (0.0, 6.0, 4.0, 6.0),
            (4.0, 6.0, 4.0, 3.0),
            (4.0, 3.0, 0.0, 3.0),
        ],
        'S' => &[
            (4.0, 6.0, 0.0, 6.0),
            (0.0, 6.0, 0.0, 3.0),
            (0.0, 3.0, 4.0, 3.0),
            (4.0, 3.0, 4.0, 0.0),
            (4.0, 0.0, 0.0, 0.0),
        ],
        'T' => &[(0.0, 6.0, 4.0, 6.0), (2.0, 6.0, 2.0, 0.0)],
        'V' => &[(0.0, 6.0, 2.0, 0.0), (2.0, 0.0, 4.0, 6.0)],
        _ => &[],
    }
}

/// How wide one character is on the lattice, including the gap after it.
const ADVANCE: f32 = 5.5;

/// The strokes for `text`, in units of the lattice, with the baseline at `y = 0` and the first
/// glyph starting at `x = 0`.
///
/// The caller scales and places them; this only knows shapes. Returned as a flat list of segments
/// because that is what the line pipeline takes.
pub fn text(s: &str) -> Vec<Stroke> {
    let mut out = Vec::new();
    let mut x = 0.0;
    for c in s.chars() {
        if c == ' ' {
            x += ADVANCE;
            continue;
        }
        for (ax, ay, bx, by) in strokes(c) {
            out.push((ax + x, *ay, bx + x, *by));
        }
        x += ADVANCE;
    }
    out
}

/// How wide `text` is on the lattice.
pub fn width(s: &str) -> f32 {
    if s.is_empty() {
        0.0
    } else {
        s.chars().count() as f32 * ADVANCE - (ADVANCE - 4.0)
    }
}

/// The lattice's height, so a caller can size a box for one line.
pub const HEIGHT: f32 = 6.0;

#[cfg(test)]
mod tests {
    use super::*;

    /// **Every character a label can carry draws something.**
    ///
    /// A glyph with no entry draws nothing, which is the right failure for an unexpected character
    /// and the wrong one for a digit. This is the list the legends actually use.
    #[test]
    fn every_character_a_legend_uses_has_strokes() {
        for c in "0123456789.-+/°ACEJKMPSTV".chars() {
            assert!(
                !strokes(c).is_empty(),
                "{c:?} is in a legend and draws nothing"
            );
        }
        // And lower case folds, so `Pa` and `m/s` read.
        for c in "acejkmpstv".chars() {
            assert!(!strokes(c).is_empty(), "{c:?} did not fold to upper case");
        }
    }

    /// **A string's strokes are its characters', shifted.**
    #[test]
    fn text_lays_its_glyphs_out_in_a_row() {
        let one = text("8");
        let two = text("88");
        assert_eq!(two.len(), one.len() * 2);
        // The second glyph is exactly one advance to the right of the first.
        let (a, b) = (one[0], two[one.len()]);
        assert!((b.0 - a.0 - ADVANCE).abs() < 1e-6, "{a:?} then {b:?}");
        // A space moves the pen and draws nothing.
        assert_eq!(text("8 8").len(), one.len() * 2);
        assert!(width("8 8") > width("88"));
    }

    /// **An unknown character costs a space and not a box.**
    #[test]
    fn an_unknown_character_draws_nothing() {
        assert!(strokes('\u{4e2d}').is_empty());
        assert_eq!(text("8\u{4e2d}8").len(), text("8").len() * 2);
    }
}
