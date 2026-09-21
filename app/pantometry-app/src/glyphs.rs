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
//! caller-chosen box. The digits, a sign, a point, an exponent, and **the whole alphabet**.
//!
//! # It held ten letters, and the argument for that was measured wrong
//!
//! It carried "the letters that spell the units this workspace reports", which was `A C E J K M P
//! S T V`, and a character with no entry drew nothing "because a missing glyph in a legend should
//! cost a reader a character and not their confidence in the number beside it". Both halves were
//! wrong together. A unit here is any string a panel carries, and they had grown past that
//! alphabet without anything saying so: `refractive index` drew as `_E__ACT_VE ___E_` on the
//! front page, `A moved` as `A M_VE_`, and `deg field, 4 = glass` as `_E_ __E___ 4 _ __ASS`.
//! Costing a reader a character is what it was supposed to do; costing them most of the word is
//! what it did, silently, in every committed figure.
//!
//! So: sixteen more letters, and a character this does not know draws a box with both diagonals
//! through it. A hole that can be seen is a hole somebody fixes.

/// A stroke, as two points on the glyph lattice.
type Stroke = (f32, f32, f32, f32);

/// What a character with no entry draws: a box with both diagonals, which is nothing this
/// alphabet spells.
///
/// **Not an empty slice, which is what it was.** A legend that silently drops the characters it
/// cannot draw reads as a legend, and the reader has no way to know a letter is missing rather
/// than absent. Three committed figures shipped that way — the front page's bench render spelled
/// its own unit `_E__ACT_VE ___E_` — and none of the three was noticed by anybody looking at it,
/// including in the session that added the third.
const UNKNOWN: &[Stroke] = &[
    (0.0, 0.0, 4.0, 0.0),
    (4.0, 0.0, 4.0, 6.0),
    (4.0, 6.0, 0.0, 6.0),
    (0.0, 6.0, 0.0, 0.0),
    (0.0, 0.0, 4.0, 6.0),
    (0.0, 6.0, 4.0, 0.0),
];

/// The strokes for one character, or the box for one this does not know.
///
/// Lower case is folded to upper: a legend in one case reads as a legend rather than as prose.
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
        'B' => &[
            (0.0, 0.0, 0.0, 6.0),
            (0.0, 6.0, 3.5, 6.0),
            (3.5, 6.0, 3.5, 3.0),
            (0.0, 3.0, 4.0, 3.0),
            (4.0, 3.0, 4.0, 0.0),
            (4.0, 0.0, 0.0, 0.0),
        ],
        'C' => &[
            (4.0, 6.0, 0.0, 6.0),
            (0.0, 6.0, 0.0, 0.0),
            (0.0, 0.0, 4.0, 0.0),
        ],
        'D' => &[
            (0.0, 0.0, 0.0, 6.0),
            (0.0, 6.0, 3.0, 6.0),
            (3.0, 6.0, 4.0, 4.0),
            (4.0, 4.0, 4.0, 2.0),
            (4.0, 2.0, 3.0, 0.0),
            (3.0, 0.0, 0.0, 0.0),
        ],
        'E' => &[
            (4.0, 6.0, 0.0, 6.0),
            (0.0, 6.0, 0.0, 0.0),
            (0.0, 0.0, 4.0, 0.0),
            (0.0, 3.0, 3.0, 3.0),
        ],
        'F' => &[
            (0.0, 0.0, 0.0, 6.0),
            (0.0, 6.0, 4.0, 6.0),
            (0.0, 3.0, 3.0, 3.0),
        ],
        'G' => &[
            (4.0, 6.0, 0.0, 6.0),
            (0.0, 6.0, 0.0, 0.0),
            (0.0, 0.0, 4.0, 0.0),
            (4.0, 0.0, 4.0, 3.0),
            (4.0, 3.0, 2.0, 3.0),
        ],
        'H' => &[
            (0.0, 0.0, 0.0, 6.0),
            (4.0, 0.0, 4.0, 6.0),
            (0.0, 3.0, 4.0, 3.0),
        ],
        'I' => &[
            (0.0, 6.0, 4.0, 6.0),
            (2.0, 6.0, 2.0, 0.0),
            (0.0, 0.0, 4.0, 0.0),
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
        'L' => &[(0.0, 6.0, 0.0, 0.0), (0.0, 0.0, 4.0, 0.0)],
        'M' => &[
            (0.0, 0.0, 0.0, 6.0),
            (0.0, 6.0, 2.0, 3.0),
            (2.0, 3.0, 4.0, 6.0),
            (4.0, 6.0, 4.0, 0.0),
        ],
        'N' => &[
            (0.0, 0.0, 0.0, 6.0),
            (0.0, 6.0, 4.0, 0.0),
            (4.0, 0.0, 4.0, 6.0),
        ],
        'O' => &[
            (0.0, 0.0, 4.0, 0.0),
            (4.0, 0.0, 4.0, 6.0),
            (4.0, 6.0, 0.0, 6.0),
            (0.0, 6.0, 0.0, 0.0),
        ],
        'P' => &[
            (0.0, 0.0, 0.0, 6.0),
            (0.0, 6.0, 4.0, 6.0),
            (4.0, 6.0, 4.0, 3.0),
            (4.0, 3.0, 0.0, 3.0),
        ],
        'Q' => &[
            (0.0, 0.0, 4.0, 0.0),
            (4.0, 0.0, 4.0, 6.0),
            (4.0, 6.0, 0.0, 6.0),
            (0.0, 6.0, 0.0, 0.0),
            (2.0, 2.0, 4.0, 0.0),
        ],
        'R' => &[
            (0.0, 0.0, 0.0, 6.0),
            (0.0, 6.0, 4.0, 6.0),
            (4.0, 6.0, 4.0, 3.0),
            (4.0, 3.0, 0.0, 3.0),
            (0.0, 3.0, 4.0, 0.0),
        ],
        // **Chamfered, because square it was `5` exactly.** The two drew the same five strokes,
        // so no legend here could tell `5 M` from `S M` and nothing said so until a test asked
        // whether any two characters render alike. The digits keep the square idiom the rest of
        // them share; the letter gives up its corners.
        'S' => &[
            (4.0, 6.0, 1.0, 6.0),
            (1.0, 6.0, 0.0, 5.0),
            (0.0, 5.0, 0.0, 3.0),
            (0.0, 3.0, 4.0, 3.0),
            (4.0, 3.0, 4.0, 1.0),
            (4.0, 1.0, 3.0, 0.0),
            (3.0, 0.0, 0.0, 0.0),
        ],
        'T' => &[(0.0, 6.0, 4.0, 6.0), (2.0, 6.0, 2.0, 0.0)],
        'U' => &[
            (0.0, 6.0, 0.0, 1.0),
            (0.0, 1.0, 2.0, 0.0),
            (2.0, 0.0, 4.0, 1.0),
            (4.0, 1.0, 4.0, 6.0),
        ],
        'V' => &[(0.0, 6.0, 2.0, 0.0), (2.0, 0.0, 4.0, 6.0)],
        'W' => &[
            (0.0, 6.0, 1.0, 0.0),
            (1.0, 0.0, 2.0, 3.0),
            (2.0, 3.0, 3.0, 0.0),
            (3.0, 0.0, 4.0, 6.0),
        ],
        'X' => &[(0.0, 6.0, 4.0, 0.0), (0.0, 0.0, 4.0, 6.0)],
        'Y' => &[
            (0.0, 6.0, 2.0, 3.0),
            (4.0, 6.0, 2.0, 3.0),
            (2.0, 3.0, 2.0, 0.0),
        ],
        'Z' => &[
            (0.0, 6.0, 4.0, 6.0),
            (4.0, 6.0, 0.0, 0.0),
            (0.0, 0.0, 4.0, 0.0),
        ],
        // **A space is the one character that draws nothing.** Everything else that is not in
        // the table draws the box below, because drawing nothing is what let three figures ship
        // with holes in their legends and nobody read them as holes.
        ' ' => &[],
        _ => UNKNOWN,
    }
}

/// Whether every character of `s` is one this can draw, spaces aside.
///
/// **A test's question, not a renderer's**, which is why it is not compiled into the binary. A
/// unit that arrives from somebody's scene file and contains a character this has no glyph for is
/// *supposed* to draw boxes — that is the whole of what the boxes are for, and refusing to draw
/// it would be worse. What must never contain one is a label this workspace *composes*: the scale
/// bar picks its own unit, `UM` and `NM` are in the table because it needed them, and a unit whose
/// letters did not draw would read as `M` — a bar wrong by a factor of a million, in a legend a
/// reader has no way to doubt.
#[cfg(test)]
pub fn can_draw(s: &str) -> bool {
    s.chars().all(|c| c == ' ' || strokes(c) != UNKNOWN)
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

    /// **The alphabet and the digits all draw something, and none of them draws the box.**
    ///
    /// This used to walk the string `"0123456789.-+/°ACEJKMPSTV"`, under a doc calling it "the
    /// list the legends actually use". It was the *table's own contents* transcribed: adding a
    /// letter to the table and to the string kept it green, and a legend using a letter in
    /// neither was invisible to it. It was green while the figure on the front page spelled its
    /// own unit `_E__ACT_VE ___E_` and the viewer's colour bar read `A M_VE_` for `A moved`.
    ///
    /// A unit string is a string. The specification is the alphabet, not the table.
    #[test]
    fn the_alphabet_and_the_digits_all_draw_something() {
        let missing: Vec<char> = ('A'..='Z')
            .chain('0'..='9')
            .chain(".-+/°".chars())
            .filter(|c| strokes(*c).is_empty() || strokes(*c) == UNKNOWN)
            .collect();
        assert!(
            missing.is_empty(),
            "these draw nothing a reader can tell from a gap: {missing:?}"
        );
        // And lower case folds, so `Pa` and `m/s` read.
        for c in 'a'..='z' {
            assert_eq!(
                strokes(c),
                strokes(c.to_ascii_uppercase()),
                "{c:?} did not fold to upper case"
            );
        }
    }

    /// **No glyph leaves its own lattice**, or it runs into the character beside it and a word
    /// reads as a smear. Sixteen glyphs written in one sitting is sixteen chances to mistype a
    /// coordinate, and the lattice is the only thing that says which ones are wrong.
    #[test]
    fn every_stroke_stays_on_the_lattice() {
        let inside = |v: f32, hi: f32| (0.0..=hi).contains(&v);
        for c in ('A'..='Z').chain('0'..='9').chain(".-+/°".chars()) {
            for &(x0, y0, x1, y1) in strokes(c) {
                assert!(
                    inside(x0, 4.0) && inside(x1, 4.0) && inside(y0, 6.0) && inside(y1, 6.0),
                    "{c:?} has a stroke ({x0}, {y0}) to ({x1}, {y1}), off a 4 by 6 lattice"
                );
            }
        }
    }

    /// **No two characters draw the same strokes.** Two letters that render identically are a
    /// legend that lies rather than one that is hard to read, and a paste from the line above is
    /// how that happens. `0` and `O` are one stroke apart on purpose; everything else is further.
    #[test]
    fn no_two_characters_draw_the_same_strokes() {
        let all: Vec<char> = ('A'..='Z')
            .chain('0'..='9')
            .chain(".-+/°".chars())
            .collect();
        for (i, &a) in all.iter().enumerate() {
            for &b in &all[i + 1..] {
                assert!(
                    strokes(a) != strokes(b),
                    "{a:?} and {b:?} draw the same strokes"
                );
            }
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

    /// **An unknown character costs a box, and it used to cost nothing.**
    ///
    /// The argument for nothing was that a missing glyph should cost a reader one character
    /// rather than their confidence in the number beside it. What it cost instead was three
    /// committed figures with most of their units missing, and nobody reading a gap as a gap. A
    /// space is still the one thing that is blank.
    #[test]
    fn an_unknown_character_leaves_a_mark() {
        for c in ['\u{4e2d}', '#', '%', '@', '\u{3bc}'] {
            assert_eq!(
                strokes(c),
                UNKNOWN,
                "{c:?} draws nothing, which reads as absent"
            );
        }
        assert_eq!(
            text("8\u{4e2d}8").len(),
            text("8").len() * 2 + UNKNOWN.len(),
            "the box is not in the line"
        );
        assert!(strokes(' ').is_empty(), "a space is a space");
        assert_eq!(text("8 8").len(), text("8").len() * 2);
    }
}
