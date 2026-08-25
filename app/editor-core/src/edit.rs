//! Changing one number in the scene, without touching any other byte of it.
//!
//! The editor's inspector could show a value and not change it, so the only editable widget in
//! the whole shell was a text box holding the raw JSON. This is the machinery that fixes that,
//! and it lives here rather than in the shell for the reason everything else does: a splice that
//! silently corrupts a file is a thing to have a test for, not an event handler.
//!
//! # Why a byte splice and not a round trip
//!
//! The obvious implementation is `from_str::<Value>`, `pointer_mut`, `to_string_pretty`. It is
//! wrong here and the reason is measurable: `serde_json` without its `preserve_order` feature —
//! which this workspace does not enable — backs an object with a `BTreeMap`, so a round trip
//! **alphabetises every key**. The shipped scenes are hand-formatted, several keys to a line,
//! `kind` and `name` first because that is the order a person reads them in. Dragging one slider
//! would reflow the file and reorder every object in it, and the diff of a scene would stop being
//! readable the first time anybody touched the editor.
//!
//! So a change here replaces exactly the bytes of the value and copies the rest verbatim. The
//! tests assert that literally: every byte outside the span is identical.
//!
//! # A literal's kind is the file's to state, not this module's to guess
//!
//! `frames` is a `usize` in the scene format and `cell_mm` is an `f64`, and nothing in a pointer
//! says which. Writing `11.0` where the format wants a count makes the scene stop parsing, and
//! writing `3` where it had `2.5` quietly changes a float literal into an integer one.
//!
//! Neither is decided from a list of key names, which would be one more place to forget a key.
//! It is read off **the literal that is already there**: a value written without `.`, `e` or `E`
//! is a whole number and stays one — a fractional drag on it is refused rather than rounded —
//! and one written with them keeps a decimal point even when it lands on an integer.

use std::ops::Range;

/// One number in the scene text that the inspector can change.
#[derive(Clone, Debug, PartialEq)]
pub struct Editable {
    /// The key, as the scene spells it: `cell_mm`, or `cells[0]` for an element of an array.
    pub label: String,
    /// Where it is, as an RFC 6901 JSON pointer: `/domains/2/cell_mm`.
    ///
    /// A pointer rather than the outliner's path, because the outliner is keyed by *name* and
    /// the file is keyed by position — and the two disagree the moment a domain has no extent,
    /// since [`crate::check`] skips those when it builds its boxes.
    pub pointer: String,
    /// What it is now.
    pub value: f64,
    /// Whether the literal in the file is a whole number, so a shell knows to step by one and
    /// not to offer a fractional drag. See this module's note on why it is read and not guessed.
    pub integral: bool,
    /// The unit the key names, or empty when the key does not name one. **Empty rather than
    /// guessed**: a unit invented for a key this table has not met is a wrong label on a number,
    /// which is worse than no label on a number.
    pub unit: &'static str,
}

/// The unit a key names, by the suffix convention the scene format already uses.
///
/// Longest suffix first, because `_m` matches the tail of `_mm` and would win otherwise. Keys
/// this does not recognise get no unit at all rather than a plausible one.
fn unit_of(key: &str) -> &'static str {
    const SUFFIXES: [(&str, &str); 14] = [
        ("_per_s2", "m/s²"),
        ("_per_m_k", "W/m·K"),
        ("_um", "µm"),
        ("_mm", "mm"),
        ("_hz", "Hz"),
        ("_pa", "Pa"),
        ("_kg", "kg"),
        ("_m", "m"),
        ("_s", "s"),
        ("_c", "°C"),
        ("_k", "K"),
        ("_j", "J"),
        ("_w", "W"),
        ("_v", "V"),
    ];
    if key == "watts" {
        return "W";
    }
    for (suffix, unit) in SUFFIXES {
        if key.ends_with(suffix) {
            return unit;
        }
    }
    ""
}

/// Every number the inspector can change for the selected outliner row.
///
/// `path` is a [`crate::Node::path`]: `/` for the scene itself, `/extents/<name>` for a domain.
/// Anything else — a run's field, a reading — yields nothing, because those are output and the
/// only honest way to change one is to change the scene that produced it.
///
/// **Nothing here enumerates domain kinds.** Every number at the top level of a domain object is
/// offered, whatever the object is, which is ARCHITECTURE.md's rule for the inspection half: a
/// panel written per domain shows an out-of-tree domain nothing at all. The cost is that a key
/// this crate has never heard of is draggable, which is the right failure — the scene's own
/// check is the authority on whether the result is legal, and it runs on every keystroke already.
///
/// Empty when the text does not parse. A file that is not JSON has no values to point at, and
/// the text box is where that gets fixed; a build failure is different and still yields the full
/// set, which is the case a person most wants to edit their way out of.
pub fn editable(text: &str, path: &str) -> Vec<Editable> {
    let Ok(root) = serde_json::from_str::<serde_json::Value>(text) else {
        return Vec::new();
    };
    let mut out = Vec::new();

    if path == "/" {
        collect(&root, "", &mut out);
        return out;
    }

    let Some(name) = domain_of(path) else {
        return out;
    };
    let Some(domains) = root.get("domains").and_then(|d| d.as_array()) else {
        return out;
    };
    // By name, not by position: `Checked::boxes` skips a domain with no extent, so the outliner's
    // nth row is not the file's nth domain.
    let Some(i) = domains
        .iter()
        .position(|d| d.get("name").and_then(|n| n.as_str()) == Some(name))
    else {
        return out;
    };
    collect(&domains[i], &format!("/domains/{i}"), &mut out);
    out
}

/// The domain an outliner path is about, if it is about one.
///
/// A domain occupies **three** rows: `/extents/<name>` before a run, and `/run/<name>` and
/// `/readings/<name>` after one. All three are the same object and all three offer the same
/// inputs, because the row a person is looking at after a run is one of the last two and making
/// them hunt for the first would be friction with nothing on the other side of it.
///
/// Deeper paths are not: `/readings/<name>/<label>` is one scalar a domain reported, which is
/// output, and the group rows `/extents`, `/run` and `/readings` are about a set rather than a
/// domain.
fn domain_of(path: &str) -> Option<&str> {
    let rest = ["/extents/", "/run/", "/readings/"]
        .iter()
        .find_map(|prefix| path.strip_prefix(prefix))?;
    (!rest.is_empty() && !rest.contains('/')).then_some(rest)
}

/// Every number directly inside `value`, and inside any array of numbers directly inside it.
///
/// One level of array, not a general walk: `cells` and a colour are lists of numbers a person
/// edits, and a nested object is a different row of the outliner's job rather than this one's.
fn collect(value: &serde_json::Value, prefix: &str, out: &mut Vec<Editable>) {
    let Some(object) = value.as_object() else {
        return;
    };
    for (key, member) in object {
        match member {
            serde_json::Value::Number(n) => out.push(Editable {
                label: key.clone(),
                pointer: format!("{prefix}/{}", escape(key)),
                value: n.as_f64().unwrap_or(f64::NAN),
                integral: n.is_i64() || n.is_u64(),
                unit: unit_of(key),
            }),
            serde_json::Value::Array(items)
                if !items.is_empty() && items.iter().all(|i| i.is_number()) =>
            {
                for (j, item) in items.iter().enumerate() {
                    let n = item.as_number().expect("checked above");
                    out.push(Editable {
                        label: format!("{key}[{j}]"),
                        pointer: format!("{prefix}/{}/{j}", escape(key)),
                        value: n.as_f64().unwrap_or(f64::NAN),
                        integral: n.is_i64() || n.is_u64(),
                        unit: unit_of(key),
                    });
                }
            }
            _ => {}
        }
    }
}

/// A key as an RFC 6901 pointer segment: `~` becomes `~0` and `/` becomes `~1`.
///
/// No scene key has either today. It is here because the *decoding* half has to exist for the
/// pointer to be well defined, and an encoder that did not match it would be a bug waiting for
/// the first key with a slash in it.
fn escape(key: &str) -> String {
    key.replace('~', "~0").replace('/', "~1")
}

/// The segments of an RFC 6901 pointer, unescaped.
fn segments(pointer: &str) -> Option<Vec<String>> {
    if pointer.is_empty() {
        return Some(Vec::new());
    }
    let rest = pointer.strip_prefix('/')?;
    Some(
        rest.split('/')
            .map(|s| s.replace("~1", "/").replace("~0", "~"))
            .collect(),
    )
}

/// Replace the number at `pointer` with `value`, leaving every other byte of `text` alone.
///
/// Fails rather than writes when the pointer names nothing, when the value there is not a
/// number, when `value` is not finite — `NaN` and the infinities have no JSON spelling, and a
/// scene containing one of them would be a file no reader can load — or when a whole number in
/// the file is handed a fraction.
pub fn set_number(text: &str, pointer: &str, value: f64) -> Result<String, String> {
    let span = value_span(text, pointer)
        .ok_or_else(|| format!("{pointer}: no such value in this scene"))?;
    let old = &text[span.clone()];
    if old.parse::<f64>().is_err() {
        return Err(format!("{pointer}: {old} is not a number"));
    }
    if !value.is_finite() {
        return Err(format!(
            "{pointer}: {value} has no JSON spelling, and a scene holding one would not load"
        ));
    }
    let integral = !old.contains(['.', 'e', 'E']);
    let rendered = if integral {
        if value.fract() != 0.0 {
            return Err(format!(
                "{pointer}: {old} is a whole number in this file and {value} is not one"
            ));
        }
        format!("{}", value as i64)
    } else {
        let s = format!("{value}");
        if s.contains(['.', 'e', 'E']) {
            s
        } else {
            format!("{s}.0")
        }
    };

    let mut out = String::with_capacity(text.len() + rendered.len());
    out.push_str(&text[..span.start]);
    out.push_str(&rendered);
    out.push_str(&text[span.end..]);
    Ok(out)
}

/// The byte range of the value at `pointer`, or `None` if the pointer names nothing.
///
/// A scanner rather than a parse, because the whole point is to know *where* a value is and a
/// parsed `Value` has thrown that away. It walks the text once per segment and never allocates
/// except to decode a key.
fn value_span(text: &str, pointer: &str) -> Option<Range<usize>> {
    let segments = segments(pointer)?;
    span_of(text.as_bytes(), 0, &segments)
}

fn skip_ws(b: &[u8], mut i: usize) -> usize {
    while i < b.len() && matches!(b[i], b' ' | b'\t' | b'\n' | b'\r') {
        i += 1;
    }
    i
}

/// Index just past the string whose opening quote is at `i`.
fn end_of_string(b: &[u8], i: usize) -> Option<usize> {
    let mut j = i + 1;
    while j < b.len() {
        match b[j] {
            // A backslash escapes whatever follows, including a quote and including another
            // backslash — skipping two is what keeps `"a\\"` from looking unterminated.
            b'\\' => j += 2,
            b'"' => return Some(j + 1),
            _ => j += 1,
        }
    }
    None
}

/// Index just past the value starting at `i`.
fn end_of_value(b: &[u8], i: usize) -> Option<usize> {
    match *b.get(i)? {
        b'"' => end_of_string(b, i),
        b'{' | b'[' => {
            let mut depth = 0usize;
            let mut j = i;
            while j < b.len() {
                match b[j] {
                    // Strings are skipped whole: a brace inside one is a character, not a nesting
                    // level, and counting it would end the object at the wrong byte.
                    b'"' => {
                        j = end_of_string(b, j)?;
                        continue;
                    }
                    b'{' | b'[' => depth += 1,
                    b'}' | b']' => {
                        depth -= 1;
                        if depth == 0 {
                            return Some(j + 1);
                        }
                    }
                    _ => {}
                }
                j += 1;
            }
            None
        }
        // A number, `true`, `false` or `null`: runs to the first delimiter.
        _ => {
            let mut j = i;
            while j < b.len() && !matches!(b[j], b',' | b'}' | b']' | b' ' | b'\t' | b'\n' | b'\r')
            {
                j += 1;
            }
            (j > i).then_some(j)
        }
    }
}

fn span_of(b: &[u8], start: usize, segments: &[String]) -> Option<Range<usize>> {
    let i = skip_ws(b, start);
    let Some((want, rest)) = segments.split_first() else {
        return Some(i..end_of_value(b, i)?);
    };
    match *b.get(i)? {
        b'{' => {
            let mut j = skip_ws(b, i + 1);
            loop {
                if *b.get(j)? != b'"' {
                    return None;
                }
                let key_end = end_of_string(b, j)?;
                let key: String =
                    serde_json::from_str(std::str::from_utf8(b.get(j..key_end)?).ok()?).ok()?;
                let colon = skip_ws(b, key_end);
                if *b.get(colon)? != b':' {
                    return None;
                }
                let value_at = skip_ws(b, colon + 1);
                if key == *want {
                    return span_of(b, value_at, rest);
                }
                j = skip_ws(b, end_of_value(b, value_at)?);
                if *b.get(j)? != b',' {
                    return None;
                }
                j = skip_ws(b, j + 1);
            }
        }
        b'[' => {
            let want: usize = want.parse().ok()?;
            let mut j = skip_ws(b, i + 1);
            let mut n = 0usize;
            loop {
                if *b.get(j)? == b']' {
                    return None;
                }
                if n == want {
                    return span_of(b, j, rest);
                }
                j = skip_ws(b, end_of_value(b, j)?);
                if *b.get(j)? != b',' {
                    return None;
                }
                j = skip_ws(b, j + 1);
                n += 1;
            }
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Hand-formatted the way the shipped scenes are: several keys to a line, `kind` and `name`
    /// first. Every test that asserts nothing else moved is asserting about *this* shape.
    const SCENE: &str = r#"{
  "title": "a block and a lamp",
  "duration_s": 0.006,
  "frames": 11,
  "domains": [
    { "kind": "lump", "name": "lamp", "volume_cm3": 12.0, "thickness_mm": 3.0,
      "initial_c": 80.0, "ambient_c": 20.0, "area_cm2": 30.0 },
    { "kind": "block", "name": "buffer", "cells": [11, 11, 11], "cell_mm": 2.0,
      "initial_c": 20.0, "material": "aluminium" }
  ]
}"#;

    /// **A change replaces the number and copies every other byte.**
    ///
    /// Asserted literally rather than by re-parsing and comparing values: a round trip that
    /// happened to produce an equivalent scene would pass a value comparison while having
    /// reflowed the file, which is the failure this module exists to avoid.
    #[test]
    fn a_change_moves_exactly_the_bytes_of_the_value() {
        let out = set_number(SCENE, "/domains/1/cell_mm", 2.5).expect("the pointer resolves");
        assert_eq!(out, SCENE.replace("\"cell_mm\": 2.0", "\"cell_mm\": 2.5"));
        assert!(
            out.contains("\"kind\": \"block\", \"name\": \"buffer\""),
            "{out}"
        );
    }

    /// **A whole number stays whole, and a fraction on one is refused rather than rounded.**
    ///
    /// `frames` is a `usize` in the format, so `11.0` is a scene that stops parsing. Nothing here
    /// knows that `frames` is a count — it reads the literal, which has no decimal point.
    #[test]
    fn a_count_keeps_its_kind_and_refuses_a_fraction() {
        let out = set_number(SCENE, "/frames", 24.0).expect("a whole number is fine");
        assert!(out.contains("\"frames\": 24,"), "{out}");
        assert!(
            !out.contains("24.0"),
            "a count must not grow a decimal point: {out}"
        );
        serde_json::from_str::<serde_json::Value>(&out).expect("and it still parses");

        let why = set_number(SCENE, "/frames", 11.5).expect_err("a fractional count is refused");
        assert!(why.contains("whole number"), "{why}");
    }

    /// **A float that lands on an integer keeps its decimal point.**
    ///
    /// `2.0` dragged to `3` must not become the literal `3`. Both parse, but the file's character
    /// is the thing being preserved, and a value that changes kind under a drag is the kind of
    /// churn that makes a diff unreadable.
    #[test]
    fn a_float_that_lands_on_a_whole_number_stays_a_float() {
        let out = set_number(SCENE, "/domains/1/cell_mm", 3.0).expect("the pointer resolves");
        assert!(out.contains("\"cell_mm\": 3.0"), "{out}");
    }

    /// **An element of an array is addressable on its own.**
    #[test]
    fn one_cell_count_can_change_without_the_others() {
        let out = set_number(SCENE, "/domains/1/cells/1", 21.0).expect("the pointer resolves");
        assert!(out.contains("\"cells\": [11, 21, 11]"), "{out}");
    }

    /// **A pointer that names nothing fails and says so**, rather than writing somewhere else.
    #[test]
    fn a_pointer_that_names_nothing_is_refused() {
        for bad in [
            "/domains/9/cell_mm",
            "/domains/1/no_such_key",
            "/domains/1/cells/7",
            "/title",
        ] {
            let r = set_number(SCENE, bad, 1.0);
            assert!(r.is_err(), "{bad} should not resolve to a number");
        }
    }

    /// **`NaN` and the infinities are refused**, because there is no JSON for them and a scene
    /// carrying one is a file nothing can load.
    #[test]
    fn a_value_json_cannot_spell_is_refused() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let why = set_number(SCENE, "/domains/1/cell_mm", bad).expect_err("refused");
            assert!(why.contains("no JSON spelling"), "{why}");
        }
    }

    /// **The inspector's list is resolved by name, not by position.**
    ///
    /// `lamp` is a `lump`, which has no extent, so `Checked::boxes` does not contain it and the
    /// outliner's first *extent* row is the second domain in the file. A list built by counting
    /// rows would edit the lamp's watts while the reader was looking at the block.
    #[test]
    fn a_domain_is_found_by_name_even_when_the_outliner_skipped_one() {
        let fields = editable(SCENE, "/extents/buffer");
        let cell = fields
            .iter()
            .find(|e| e.label == "cell_mm")
            .expect("the block's cell size is editable");
        assert_eq!(cell.pointer, "/domains/1/cell_mm");
        assert_eq!(cell.value, 2.0);
        assert!(!cell.integral);
        assert_eq!(cell.unit, "mm");

        // And the array came through element by element, in order.
        let cells: Vec<&Editable> = fields
            .iter()
            .filter(|e| e.label.starts_with("cells"))
            .collect();
        assert_eq!(cells.len(), 3);
        assert_eq!(cells[0].pointer, "/domains/1/cells/0");
        assert!(cells[0].integral, "a cell count is a whole number");
    }

    /// **The scene's own numbers are editable from the root row**, and its strings are not.
    #[test]
    fn the_root_row_offers_the_scenes_own_numbers() {
        let fields = editable(SCENE, "/");
        let labels: Vec<&str> = fields.iter().map(|e| e.label.as_str()).collect();
        assert!(labels.contains(&"duration_s"), "{labels:?}");
        assert!(labels.contains(&"frames"), "{labels:?}");
        assert!(
            !labels.contains(&"title"),
            "a string is not a number: {labels:?}"
        );
        assert!(
            !labels.contains(&"domains"),
            "an array of objects is not: {labels:?}"
        );
    }

    /// **All three of a domain's rows offer the same inputs.**
    ///
    /// After a run the reader is looking at `/run/buffer` or `/readings/buffer`, not at the
    /// extent row, and the scene values are the same object's either way.
    #[test]
    fn every_row_that_names_a_domain_offers_its_inputs() {
        let from_extent = editable(SCENE, "/extents/buffer");
        assert!(!from_extent.is_empty());
        for same in ["/run/buffer", "/readings/buffer"] {
            assert_eq!(editable(SCENE, same), from_extent, "{same}");
        }
    }

    /// **A row that is output, or a group, offers nothing.** A reading is a result; the way to
    /// change one is to change the scene that produced it.
    #[test]
    fn an_output_or_group_row_is_not_editable() {
        for not_a_domain in [
            "/readings/buffer/peak",
            "/extents",
            "/run",
            "/readings",
            "/extents/nothing-called-this",
            "/buffer",
        ] {
            assert!(
                editable(SCENE, not_a_domain).is_empty(),
                "{not_a_domain} should offer nothing"
            );
        }
    }

    /// **What comes out is still a scene the builder accepts**, not merely still JSON.
    ///
    /// The splice is a text operation and every other test here judges it as one. This is the
    /// only one that hands the result to the thing that has to load it, which is the difference
    /// between "the bytes look right" and "the editor did not just break the file".
    #[test]
    fn a_changed_scene_still_builds() {
        let fixture = crate::check(SCENE, &crate::OnDisk);
        assert!(
            fixture.error.is_none(),
            "the fixture has to be a real scene for this test to mean anything: {:?}",
            fixture.error
        );
        for (pointer, value) in [
            ("/domains/1/cell_mm", 2.5),
            ("/domains/1/cells/0", 9.0),
            ("/frames", 3.0),
            ("/duration_s", 0.01),
        ] {
            let out = set_number(SCENE, pointer, value).expect(pointer);
            let checked = crate::check(&out, &crate::OnDisk);
            assert!(checked.error.is_none(), "{pointer}: {:?}", checked.error);
        }
    }

    /// **Text that does not parse offers nothing**, rather than pointing into a file whose shape
    /// is unknown.
    #[test]
    fn unparseable_text_offers_nothing() {
        assert!(editable("{ \"title\": ", "/").is_empty());
    }

    /// **A brace inside a string does not end an object.**
    ///
    /// The scanner counts nesting, and a title containing `}` would close the scene early and put
    /// every span after it in the wrong place. Worth its own case because the failure is silent:
    /// the splice would still produce a file, just not the right one.
    #[test]
    fn a_brace_inside_a_string_is_not_a_nesting_level() {
        let awkward = r#"{ "title": "a } and a { and a \" quote", "frames": 3 }"#;
        let out = set_number(awkward, "/frames", 4.0).expect("the pointer resolves");
        assert_eq!(
            out,
            r#"{ "title": "a } and a { and a \" quote", "frames": 4 }"#
        );
    }
}
