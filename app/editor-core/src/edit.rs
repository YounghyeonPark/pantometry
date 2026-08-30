//! Changing one value in the scene, without touching any other byte of it.
//!
//! The editor's inspector could show a value and not change it, so the only editable widget in
//! the whole shell was a text box holding the raw JSON. This is the machinery that fixes that,
//! and it lives here rather than in the shell for the reason everything else does: a splice that
//! silently corrupts a file is a thing to have a test for, not an event handler.
//!
//! # What it covers, counted rather than assumed
//!
//! A census of every value in the twenty-eight scenes that shipped when it was taken: **358 numbers, 237 strings, 59
//! arrays of numbers, 46 arrays of objects, 44 nested objects, one array of arrays — and no
//! booleans at all.** Three things follow, and none of them was obvious before the count.
//!
//! Strings are not a minority case; they are two fifths of the file, and an inspector without
//! them is an inspector that cannot change a material. Nesting is not exotic either: a
//! `hot_spot`'s `above_k` and a `region`'s `material` are exactly the numbers a person reaches
//! for, and a one-level walk — which is what this module did first — shows a domain while hiding
//! half of what the domain says. So the walk is a full one.
//!
//! And [`Value`] has no flag variant, because there is nothing for it to point at.
//! It goes in the day the format grows a boolean.
//!
//! The walk is uncapped, which is also measured: the widest domain in any shipped scene yields
//! **55** fields and the median is **8**. A cap would have to drop rows, and a panel that
//! silently omits part of the scene it claims to describe is this workspace's oldest failure.
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

/// One value in the scene text that the inspector can change.
#[derive(Clone, Debug, PartialEq)]
pub struct Editable {
    /// The key as the scene spells it, with the route to it when it is not at the top level:
    /// `cell_mm`, `cells[0]`, `hot_spot.above_k`, `regions[1].material`.
    pub label: String,
    /// Where it is, as an RFC 6901 JSON pointer: `/domains/2/cell_mm`.
    ///
    /// A pointer rather than the outliner's path, because the outliner is keyed by *name* and
    /// the file is keyed by position — and the two disagree the moment a domain has no extent,
    /// since [`crate::check`] skips those when it builds its boxes.
    pub pointer: String,
    /// What it is now, and what kind of widget it wants.
    pub value: Value,
    /// The unit the key names, or empty when the key does not name one. **Empty rather than
    /// guessed**: a unit invented for a key this table has not met is a wrong label on a number,
    /// which is worse than no label on a number.
    pub unit: &'static str,
}

/// What kind of value an [`Editable`] is, and what it would take to change it.
///
/// Two variants and not three: a census of the twenty-eight scenes of the time counts **358 numbers,
/// 237 strings and no booleans at all**, so a flag variant would be a widget with nothing to
/// point at. It goes in the day the format grows one.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    /// A number, and whether the file writes it as a whole one — see this module's note on why
    /// that is read off the literal rather than decided from the key.
    Number {
        /// What it is now.
        now: f64,
        /// Whether the literal has no `.`, `e` or `E`, so a shell steps it by one and this
        /// module refuses a fraction for it.
        integral: bool,
    },
    /// A string.
    Text {
        /// What it is now.
        now: String,
        /// The values this key is known to accept, or empty when it is free text.
        ///
        /// Non-empty for exactly one key today, `material`, and the list is **the catalogue plus
        /// whatever the scene declared** rather than the catalogue alone — a menu that omitted a
        /// scene's own `materials` and `composites` would be a menu that cannot express the file
        /// it is editing.
        choices: Vec<String>,
    },
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

/// Every value the inspector can change for the selected outliner row.
///
/// `path` is a [`crate::Node::path`]: `/` for the scene itself, `/extents/<name>` for a domain —
/// and `/run/<name>` or `/readings/<name>`, which are the same domain seen after a run. Anything
/// else yields nothing, because a reading is output and the only honest way to change one is to
/// change the scene that produced it.
///
/// **Nothing here enumerates domain kinds.** Every value in a domain's object is offered, at
/// whatever depth, whatever the object is — which is ARCHITECTURE.md's rule for the inspection
/// half: a panel written per domain shows an out-of-tree domain nothing at all. The two keys held
/// back, `kind` and `name`, are about the format rather than about any domain. The cost is that a
/// key this crate has never heard of is editable, which is the right failure — the scene's own
/// check is the authority on whether the result is legal, and it runs on every keystroke already.
///
/// Empty when the text does not parse. A file that is not JSON has no values to point at, and
/// the text box is where that gets fixed; a build failure is different and still yields the full
/// set, which is the case a person most wants to edit their way out of.
pub fn editable(text: &str, path: &str) -> Vec<Editable> {
    let Ok(root) = serde_json::from_str::<serde_json::Value>(text) else {
        return Vec::new();
    };
    let choices = Choices::of(&root);
    let mut out = Vec::new();

    if path == "/" {
        // Everything the scene says about itself, and **not** `domains` — every domain has rows
        // of its own, and folding them in here would put the whole file under one selection.
        // `materials` and `composites` are collected, because a declared substance is a thing the
        // scene states and nothing else in the outliner speaks for it.
        if let Some(object) = root.as_object() {
            for (key, member) in object {
                if key == "domains" || STRUCTURAL.contains(&key.as_str()) {
                    continue;
                }
                collect(
                    member,
                    &format!("/{}", escape(key)),
                    key,
                    &choices,
                    &mut out,
                );
            }
        }
        return out;
    }

    let Some(name) = domain_named(path) else {
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
    collect(
        &domains[i],
        &format!("/domains/{i}"),
        "",
        &choices,
        &mut out,
    );
    out
}

/// The domain an outliner path is about, if it is about one.
///
/// Public because a shell needs the same answer for a different question — whether the selected
/// row is a thing that can be deleted — and two crates deciding separately what counts as a
/// domain row is two chances to disagree about it.
///
/// A domain occupies **three** rows: `/extents/<name>` before a run, and `/run/<name>` and
/// `/readings/<name>` after one. All three are the same object and all three offer the same
/// inputs, because the row a person is looking at after a run is one of the last two and making
/// them hunt for the first would be friction with nothing on the other side of it.
///
/// Deeper paths are not: `/readings/<name>/<label>` is one scalar a domain reported, which is
/// output, and the group rows `/extents`, `/run` and `/readings` are about a set rather than a
/// domain.
pub fn domain_named(path: &str) -> Option<&str> {
    let rest = ["/extents/", "/run/", "/readings/"]
        .iter()
        .find_map(|prefix| path.strip_prefix(prefix))?;
    (!rest.is_empty() && !rest.contains('/')).then_some(rest)
}

/// Keys that are **structure** rather than value, and are not offered.
///
/// `kind`, because every other key in the object is only meaningful for the kind it is and the
/// format is `deny_unknown_fields` — changing `"room"` to `"block"` makes the eight keys beside
/// it unknown at once, so every edit of it would be an error. A widget whose every use fails is
/// a trap rather than a feature.
///
/// `name`, for two reasons that are both about it being an **identifier and not a value**.
/// Counted across the shipped scenes it appears 65 times and is referred to by five other keys —
/// `tracks`, `follows`, `onto`, `from`, `to` — plus every key of `poses`, so a rename is a
/// multi-site edit this cannot make atomically. And the outliner's own selection is keyed by the
/// name, so typing a new one a character at a time would delete the row being typed into after
/// the first keystroke: the field would vanish under the cursor.
///
/// Both stay jobs for the text box, which is where a change that rewrites more than one place
/// belongs. Naming two keys of the format is not enumerating domains: nothing here knows what
/// kinds exist, only that whatever they are, these two keys are how the file refers to them.
const STRUCTURAL: [&str; 2] = ["kind", "name"];

/// Every value inside `value` a person could change, however deeply it sits.
///
/// A full walk, not one level. Measured on the shipped scenes, the levels matter: a `hot_spot`'s
/// `above_k` and a `region`'s `material` are nested, and there are **44 nested objects and 46
/// arrays of objects** across the twenty-eight it covered. A one-level inspector shows a domain and
/// hides half of what the domain says.
///
/// No cap on how many come out, and that is measured rather than hoped: the widest domain in any
/// shipped scene yields **55** fields and the median is **8**, which is a scroll area's business
/// and not this function's. A cap would have to drop rows, and a panel that silently omits part
/// of the scene it claims to describe is the failure shape this workspace keeps finding.
fn collect(
    value: &serde_json::Value,
    pointer: &str,
    label: &str,
    choices: &Choices,
    out: &mut Vec<Editable>,
) {
    let join = |key: &str| {
        if label.is_empty() {
            key.to_string()
        } else {
            format!("{label}.{key}")
        }
    };
    match value {
        serde_json::Value::Object(object) => {
            for (key, member) in object {
                if STRUCTURAL.contains(&key.as_str()) {
                    continue;
                }
                collect(
                    member,
                    &format!("{pointer}/{}", escape(key)),
                    &join(key),
                    choices,
                    out,
                );
            }
        }
        serde_json::Value::Array(items) => {
            for (j, item) in items.iter().enumerate() {
                collect(
                    item,
                    &format!("{pointer}/{j}"),
                    &format!("{label}[{j}]"),
                    choices,
                    out,
                );
            }
        }
        serde_json::Value::Number(n) => out.push(Editable {
            label: label.to_string(),
            pointer: pointer.to_string(),
            value: Value::Number {
                now: n.as_f64().unwrap_or(f64::NAN),
                integral: n.is_i64() || n.is_u64(),
            },
            unit: unit_of(last_key(label)),
        }),
        serde_json::Value::String(s) => out.push(Editable {
            label: label.to_string(),
            pointer: pointer.to_string(),
            value: Value::Text {
                now: s.clone(),
                choices: choices.for_key(last_key(label)),
            },
            unit: "",
        }),
        // `true`, `false` and `null`. No shipped scene contains one; a widget for a shape the
        // format does not use would be a widget nothing has ever exercised.
        _ => {}
    }
}

/// The key at the end of a dotted label, with any `[i]` stripped.
///
/// `hot_spot.above_k` is a `_k`, and `radii_m[2]` is an `_m`. Taking the unit from the whole
/// label would find neither, and taking it from the outermost key would call a nested
/// temperature by its parent's name.
fn last_key(label: &str) -> &str {
    let tail = label.rsplit('.').next().unwrap_or(label);
    tail.split('[').next().unwrap_or(tail)
}

/// The values a string key is known to accept.
///
/// Built from the scene rather than from a constant, because a scene declares its own materials
/// and composites and a menu that offered only the catalogue could not express the file it was
/// editing. Empty for every key but `material`, which is the one place the format has a set and
/// [`crate::MATERIALS`] is already re-exported for a shell's menus.
struct Choices {
    materials: Vec<String>,
}

impl Choices {
    fn of(root: &serde_json::Value) -> Choices {
        let mut materials: Vec<String> = crate::MATERIALS.iter().map(|m| m.to_string()).collect();
        for declared in ["materials", "composites"] {
            if let Some(object) = root.get(declared).and_then(|d| d.as_object()) {
                materials.extend(object.keys().cloned());
            }
        }
        materials.sort();
        materials.dedup();
        Choices { materials }
    }

    fn for_key(&self, key: &str) -> Vec<String> {
        if key == "material" {
            self.materials.clone()
        } else {
            Vec::new()
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

    Ok(splice(text, span, &rendered))
}

/// Replace the string at `pointer` with `value`, leaving every other byte of `text` alone.
///
/// The counterpart to [`set_number`], and the same promise. `value` is written through
/// `serde_json`, so a name with a quote, a backslash or a newline in it comes out as valid JSON
/// rather than as a file that stops parsing — which is the one way a text field can break a scene
/// that a number field cannot.
///
/// Fails when the pointer names nothing or names something that is not a string. Refusing on the
/// second is the useful half: a pointer aimed at a number would otherwise turn `2.0` into `"2.0"`,
/// and the scene would fail to load with an error about a type rather than about an edit.
pub fn set_text(text: &str, pointer: &str, value: &str) -> Result<String, String> {
    let span = value_span(text, pointer)
        .ok_or_else(|| format!("{pointer}: no such value in this scene"))?;
    if !text[span.clone()].starts_with('"') {
        return Err(format!(
            "{pointer}: {} is not a string",
            &text[span.clone()]
        ));
    }
    let rendered = serde_json::to_string(value).map_err(|e| format!("{pointer}: {e}"))?;
    Ok(splice(text, span, &rendered))
}

/// `text` with `span` replaced by `rendered`, and every other byte copied.
///
/// One place, so the two setters cannot differ about what "leaving the rest alone" means.
fn splice(text: &str, span: Range<usize>, rendered: &str) -> String {
    let mut out = String::with_capacity(text.len() + rendered.len());
    out.push_str(&text[..span.start]);
    out.push_str(rendered);
    out.push_str(&text[span.end..]);
    out
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

/// Add a domain of `kind` to the scene, with a name nothing else is using.
///
/// The template comes from [`pantometry_world::templates::TEMPLATES`], which lives beside the
/// format rather than here — see that module for how a list of nineteen is kept in step with an
/// enum nobody can enumerate.
///
/// # What "a name nothing else is using" costs
///
/// A template's name is its kind, so adding two blocks would produce two domains called `block`,
/// and the format's own error for that arrives at build time as a name collision. The second gets
/// `block 2`, the third `block 3`. The number is a suffix rather than a rename of the first,
/// because renaming a domain the scene already has would break every key pointing at it — the
/// same reason `name` is not draggable in the inspector.
///
/// # Two kinds arrive incomplete, and that is the format's doing
///
/// `beam` states what it shines `onto` and `structure` states the block it `follows`. Both are
/// required and neither can be guessed, so the template carries a placeholder and the scene has
/// to be told where to point it. `structure`'s is refused by the build, which the inspector shows
/// immediately; `beam`'s builds and is stopped by the conservation audit at the first step, which
/// only appears when the scene is run. Both are pinned by tests in `pantometry-world`.
pub fn add_domain(text: &str, kind: &str) -> Result<String, String> {
    let template = pantometry_world::templates::TEMPLATES
        .iter()
        .find(|(k, _)| *k == kind)
        .map(|(_, t)| *t)
        .ok_or_else(|| format!("no template for a {kind}"))?;

    let root: serde_json::Value =
        serde_json::from_str(text).map_err(|e| format!("this scene does not parse: {e}"))?;
    let domains = root
        .get("domains")
        .and_then(|d| d.as_array())
        .ok_or_else(|| "this scene has no `domains` array to add to".to_string())?;

    let taken: Vec<&str> = domains
        .iter()
        .filter_map(|d| d.get("name").and_then(|n| n.as_str()))
        .collect();
    let mut name = kind.to_string();
    let mut n = 1;
    while taken.contains(&name.as_str()) {
        n += 1;
        name = format!("{kind} {n}");
    }
    let template = set_text(template, "/name", &name)?;

    // Where the last element starts on its line, so the new one lines up with it rather than
    // with the margin. A scene written by hand keeps its shape when the editor adds to it.
    let array = value_span(text, "/domains").ok_or_else(|| "no `domains` array".to_string())?;
    match domains.len() {
        0 => {
            // `[]` or `[ ]`: put the element between the brackets and let the one-line form be.
            let inner = array.start + 1..array.end - 1;
            Ok(splice(text, inner, &format!(" {template} ")))
        }
        last => {
            let end = value_span(text, &format!("/domains/{}", last - 1))
                .ok_or_else(|| "the last domain has no span".to_string())?;
            let indent = indent_of(text, end.start);
            Ok(splice(
                text,
                end.end..end.end,
                &format!(",\n{indent}{template}"),
            ))
        }
    }
}

/// Remove the domain called `name`, and the comma that held it in place.
///
/// **Nothing here follows the references.** A `tracks`, a `follows`, an `onto` or a `poses` key
/// naming the removed domain is left as it is, because the alternative is this function deciding
/// what a dangling coupling should become — and there is no answer to that which is right more
/// often than the person's. The check runs on the same frame and names what is now dangling,
/// which is the report a person can act on.
pub fn remove_domain(text: &str, name: &str) -> Result<String, String> {
    let root: serde_json::Value =
        serde_json::from_str(text).map_err(|e| format!("this scene does not parse: {e}"))?;
    let domains = root
        .get("domains")
        .and_then(|d| d.as_array())
        .ok_or_else(|| "this scene has no `domains` array".to_string())?;
    let i = domains
        .iter()
        .position(|d| d.get("name").and_then(|n| n.as_str()) == Some(name))
        .ok_or_else(|| format!("this scene defines no domain called {name:?}"))?;

    let span = value_span(text, &format!("/domains/{i}"))
        .ok_or_else(|| format!("{name}: no span for domains[{i}]"))?;

    // The element alone is not a valid removal: the comma beside it would be left behind. Take
    // the gap on whichever side has one, which is the previous element's end when there is a
    // previous element and the next element's start otherwise.
    let cut = if i > 0 {
        let prev = value_span(text, &format!("/domains/{}", i - 1))
            .ok_or_else(|| "no span for the previous domain".to_string())?;
        prev.end..span.end
    } else if i + 1 < domains.len() {
        let next = value_span(text, &format!("/domains/{}", i + 1))
            .ok_or_else(|| "no span for the next domain".to_string())?;
        span.start..next.start
    } else {
        span
    };
    Ok(splice(text, cut, ""))
}

/// Where a domain's origin sits in the world, in metres.
///
/// `[0, 0, 0]` when the scene says nothing, which is what an absent `poses` entry means — the
/// guarantee `a_placed_part` pins, and the reason every scene written before the key existed
/// still places its domains where it always did.
pub fn pose_of(text: &str, name: &str) -> [f64; 3] {
    let mut at = [0.0; 3];
    let Ok(root) = serde_json::from_str::<serde_json::Value>(text) else {
        return at;
    };
    let Some(stated) = root
        .get("poses")
        .and_then(|p| p.get(name))
        .and_then(|p| p.get("at_m"))
    else {
        return at;
    };
    for (i, slot) in at.iter_mut().enumerate() {
        if let Some(v) = stated.get(i).and_then(|v| v.as_f64()) {
            *slot = v;
        }
    }
    at
}

/// Move a domain, writing whatever part of `poses` is not there yet.
///
/// # The first write that creates rather than replaces
///
/// [`set_number`] and [`set_text`] replace the bytes of a value that already exists, which is
/// every edit the inspector makes to a domain's own fields. A position is different: **no shipped
/// scene states `poses` at all** — zero of the twenty-nine — so moving anything means writing a
/// key that is not in the file. Three levels of it can be missing, and each is handled where it
/// is found: the `at_m` array, the domain's entry, or the whole `poses` object.
///
/// What does not change is the promise. A new member is appended after the last one at the
/// indentation of the object it goes into, and every other byte is copied — a scene that gains a
/// pose keeps its formatting, and one that already had the key has three numbers replaced.
///
/// `poses` is a map beside `materials` rather than a field on each domain, so this writes at the
/// top level and not inside the domain's object. That is the format's decision and the reason is
/// in [`pantometry_world::Scene::poses`]: `serde(flatten)` would have silently disabled
/// `deny_unknown_fields`.
pub fn set_pose(text: &str, name: &str, at_m: [f64; 3]) -> Result<String, String> {
    if let Some(bad) = at_m.iter().find(|v| !v.is_finite()) {
        return Err(format!(
            "{name}: a position of {bad} has no JSON spelling, and a scene holding one would              not load"
        ));
    }
    let root: serde_json::Value =
        serde_json::from_str(text).map_err(|e| format!("this scene does not parse: {e}"))?;
    if !root
        .get("domains")
        .and_then(|d| d.as_array())
        .is_some_and(|d| {
            d.iter()
                .any(|x| x.get("name").and_then(|n| n.as_str()) == Some(name))
        })
    {
        // A pose for a domain that is not there builds nothing and is refused by the check with
        // the same sentence; catching it here means the file never holds the mistake at all.
        return Err(format!("this scene defines no domain called {name:?}"));
    }

    let numbers = format!(
        "[{}, {}, {}]",
        as_float(at_m[0]),
        as_float(at_m[1]),
        as_float(at_m[2])
    );
    let escaped = escape(name);

    if let Some(span) = value_span(text, &format!("/poses/{escaped}/at_m")) {
        return Ok(splice(text, span, &numbers));
    }
    if let Some(entry) = value_span(text, &format!("/poses/{escaped}")) {
        return append_member(text, entry, &format!(r#""at_m": {numbers}"#));
    }
    if let Some(poses) = value_span(text, "/poses") {
        return append_member(
            text,
            poses,
            &format!(r#""{name}": {{ "at_m": {numbers} }}"#),
        );
    }
    let root_span = value_span(text, "").ok_or_else(|| "this scene has no root".to_string())?;
    append_member(
        text,
        root_span,
        &format!(r#""poses": {{ "{name}": {{ "at_m": {numbers} }} }}"#),
    )
}

/// How a domain is turned: an axis, and degrees about it.
///
/// `([0, 0, 1], 0)` when the scene says nothing — an identity rotation with a real axis, because
/// the format refuses a zero one and a default that cannot be written back is not a default.
pub fn turn_of(text: &str, name: &str) -> ([f64; 3], f64) {
    let mut axis = [0.0, 0.0, 1.0];
    let mut degrees = 0.0;
    let Ok(root) = serde_json::from_str::<serde_json::Value>(text) else {
        return (axis, degrees);
    };
    let Some(turn) = root
        .get("poses")
        .and_then(|p| p.get(name))
        .and_then(|p| p.get("turn"))
    else {
        return (axis, degrees);
    };
    if let Some(stated) = turn.get("axis") {
        for (i, slot) in axis.iter_mut().enumerate() {
            if let Some(v) = stated.get(i).and_then(|v| v.as_f64()) {
                *slot = v;
            }
        }
    }
    if let Some(d) = turn.get("degrees").and_then(|d| d.as_f64()) {
        degrees = d;
    }
    (axis, degrees)
}

/// Turn a domain, writing whatever part of `poses` is not there yet.
///
/// The counterpart to [`set_pose`], one level deeper: a `turn` sits inside the entry that
/// `set_pose` writes, so there are four places the path can stop rather than three.
///
/// **Nothing else has to change for the result to be editable.** Once the key exists the
/// inspector's generic walk finds `poses.<name>.turn.axis[0]` and `poses.<name>.turn.degrees`
/// like any other value — measured rather than assumed. Creating the key is the whole job.
///
/// Refuses what the format refuses, before the file is touched: a value that is not finite, and
/// an axis with no direction. `PoseSpec::to_pose` gives the same two answers a frame later, and
/// its note on the second is the reason — normalising a zero vector gives a `NaN` rather than an
/// error, so it has to be caught by somebody.
pub fn set_turn(text: &str, name: &str, axis: [f64; 3], degrees: f64) -> Result<String, String> {
    if !degrees.is_finite() || !axis.iter().all(|v| v.is_finite()) {
        return Err(format!(
            "{name}: a turn of {degrees} degrees about {axis:?} has no JSON spelling"
        ));
    }
    if axis.iter().all(|v| *v == 0.0) {
        return Err(format!(
            "{name}: an axis of {axis:?} has no direction, and a rotation needs one"
        ));
    }
    let root: serde_json::Value =
        serde_json::from_str(text).map_err(|e| format!("this scene does not parse: {e}"))?;
    if !root
        .get("domains")
        .and_then(|d| d.as_array())
        .is_some_and(|d| {
            d.iter()
                .any(|x| x.get("name").and_then(|n| n.as_str()) == Some(name))
        })
    {
        return Err(format!("this scene defines no domain called {name:?}"));
    }

    let body = format!(
        r#"{{ "axis": [{}, {}, {}], "degrees": {} }}"#,
        as_float(axis[0]),
        as_float(axis[1]),
        as_float(axis[2]),
        as_float(degrees)
    );
    let escaped = escape(name);

    if let Some(span) = value_span(text, &format!("/poses/{escaped}/turn")) {
        return Ok(splice(text, span, &body));
    }
    if let Some(entry) = value_span(text, &format!("/poses/{escaped}")) {
        return append_member(text, entry, &format!(r#""turn": {body}"#));
    }
    if let Some(poses) = value_span(text, "/poses") {
        return append_member(text, poses, &format!(r#""{name}": {{ "turn": {body} }}"#));
    }
    let root_span = value_span(text, "").ok_or_else(|| "this scene has no root".to_string())?;
    append_member(
        text,
        root_span,
        &format!(r#""poses": {{ "{name}": {{ "turn": {body} }} }}"#),
    )
}

/// A float that stays a float: `3` would be a different literal, and `at_m` is three `f64`.
fn as_float(v: f64) -> String {
    let s = format!("{v}");
    if s.contains(['.', 'e', 'E']) {
        s
    } else {
        format!("{s}.0")
    }
}

/// Put `member` into the object at `object`, after whatever is already in it.
///
/// Appended rather than inserted in sorted position, because the file's key order is the author's
/// and this is not the code that gets to reorder it. The indentation is the object's own plus two
/// spaces, which is what the scenes use.
fn append_member(text: &str, object: Range<usize>, member: &str) -> Result<String, String> {
    if !text[object.clone()].starts_with('{') {
        return Err(format!(
            "cannot add a key to {}, which is not an object",
            &text[object.clone()]
        ));
    }
    let inner = object.start + 1..object.end - 1;
    if text[inner.clone()].trim().is_empty() {
        return Ok(splice(text, inner, &format!(" {member} ")));
    }
    // Just past the last member's final byte, which is the last non-space before the brace.
    let last = object.start + 1 + text[inner].trim_end().len();
    let indent = format!("{}  ", indent_of(text, object.start));
    Ok(splice(
        text,
        last..last,
        &format!(
            ",
{indent}{member}"
        ),
    ))
}

/// How far along `axis` a drag moves the thing at `origin`, in world units.
///
/// The arithmetic behind a translate handle, kept here rather than in the shell because it is a
/// function of a camera and two vectors and nothing about an event loop.
///
/// # The formula, and the case it refuses
///
/// A handle is an axis through the object. Project the origin and a point one unit along the
/// axis; the difference is where that unit of world motion *goes on screen*. The drag is then
/// projected onto that direction — `t = (d · s) / (s · s)` — which is the least-squares answer to
/// "how far along the axis did the pointer mean to go", and the only answer that behaves when the
/// pointer wanders off the handle.
///
/// **An axis pointing at the camera returns zero.** Its screen direction is a point, `s · s` is
/// nearly nothing, and the division would turn a one-pixel twitch into a leap across the scene.
/// A handle you cannot see is a handle you cannot drag, and refusing is the honest answer; the
/// alternative is a gizmo that occasionally throws the object out of the world.
///
/// # Exact for a small drag, second order for a large one
///
/// The projection divides by depth, so moving along the axis changes how far the object is and
/// the screen-to-world map is not linear. `s` is the map at the *starting* position, so the error
/// grows with the square of the drag — measured in `a_drag_lands_where_the_pointer_went`, which
/// halves the drag and watches the error quarter. That is the right shape for a gizmo: a pointer
/// is moved continuously, and each frame's drag is a few pixels.
pub fn drag_along_axis(
    camera: &viewer_core::Camera,
    frame: &viewer_core::Framing,
    aspect: f64,
    origin: [f64; 3],
    axis: [f64; 3],
    screen: [f64; 2],
) -> f64 {
    // **A short step, not a whole unit.** `s` has to be the screen displacement *per unit of
    // axis at this position*, and the projection divides by depth, so the secant across a whole
    // unit is not that — it is the average over a stretch where the object gets nearer or further.
    // Measured with the whole-unit version: the pointer asked for a move and the object went
    // **77%** of the way, at every drag size, converging to a fixed 23.5% relative error rather
    // than to zero. A step of a ten-thousandth of the framed span makes this the derivative.
    let length = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    if length < 1e-30 {
        return 0.0;
    }
    let eps = 1e-4 * frame.span.max(1e-12) / length;
    let along = [
        origin[0] + axis[0] * eps,
        origin[1] + axis[1] * eps,
        origin[2] + axis[2] * eps,
    ];
    let here = camera.project(origin, frame, aspect);
    let there = camera.project(along, frame, aspect);
    let s = [(there.x - here.x) / eps, (there.y - here.y) / eps];
    let len2 = s[0] * s[0] + s[1] * s[1];
    // Edge-on. `s` is now normalised device units *per unit of axis*, and the window is two
    // across, so this is a handle one unit of which draws shorter than a thousandth of the screen.
    if len2 < 1e-6 {
        return 0.0;
    }
    (screen[0] * s[0] + screen[1] * s[1]) / len2
}

/// The whitespace between the start of `at`'s line and `at` itself.
fn indent_of(text: &str, at: usize) -> String {
    let line_start = text[..at].rfind('\n').map_or(0, |n| n + 1);
    text[line_start..at]
        .chars()
        .take_while(|c| c.is_whitespace())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn labels(fields: &[Editable]) -> Vec<String> {
        fields.iter().map(|e| e.label.clone()).collect()
    }

    /// Hand-formatted the way the shipped scenes are: several keys to a line, `kind` and `name`
    /// first. Every test that asserts nothing else moved is asserting about *this* shape.
    ///
    /// It carries one of each shape that census found — a
    /// nested object (`hot_spot`), an array of objects (`regions`), an array of numbers
    /// (`cells`), strings and numbers — because a fixture that only holds the easy shapes is a
    /// fixture that agrees with a one-level walk.
    const SCENE: &str = r#"{
  "title": "a block and a lamp",
  "duration_s": 0.006,
  "frames": 11,
  "domains": [
    { "kind": "lump", "name": "lamp", "volume_cm3": 12.0, "thickness_mm": 3.0,
      "initial_c": 80.0, "ambient_c": 20.0, "area_cm2": 30.0 },
    { "kind": "block", "name": "buffer", "cells": [11, 11, 11], "cell_mm": 2.0,
      "initial_c": 20.0, "material": "aluminium",
      "hot_spot": { "at": [5, 5, 5], "above_k": 60.0 },
      "regions": [
        { "material": "copper", "from": [0, 0, 0], "to": [2, 2, 2] }
      ] }
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
        assert_eq!(
            cell.value,
            Value::Number {
                now: 2.0,
                integral: false
            }
        );
        assert_eq!(cell.unit, "mm");

        // And the array came through element by element, in order.
        let cells: Vec<&Editable> = fields
            .iter()
            .filter(|e| e.label.starts_with("cells"))
            .collect();
        assert_eq!(cells.len(), 3);
        assert_eq!(cells[0].pointer, "/domains/1/cells/0");
        assert_eq!(
            cells[0].value,
            Value::Number {
                now: 11.0,
                integral: true
            },
            "a cell count is a whole number"
        );
    }

    /// **The root row offers what the scene says about itself, and not its domains.**
    ///
    /// The title is there — a string is editable now — and nothing under `domains` is, because
    /// every domain has rows of its own and folding them in here would put the whole file under
    /// one selection.
    #[test]
    fn the_root_row_offers_the_scene_and_not_its_domains() {
        let fields = editable(SCENE, "/");
        let labels: Vec<&str> = fields.iter().map(|e| e.label.as_str()).collect();
        assert!(labels.contains(&"duration_s"), "{labels:?}");
        assert!(labels.contains(&"frames"), "{labels:?}");
        assert!(labels.contains(&"title"), "{labels:?}");
        assert!(
            !labels.iter().any(|l| l.starts_with("domains")),
            "a domain has its own row: {labels:?}"
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

    /// **A value nested inside an object is reached, and it is labelled by its route.**
    ///
    /// `hot_spot.above_k` is the case the one-level walk missed: a number a person plainly wants
    /// to change, sitting one object down, invisible to an inspector that only read the top.
    #[test]
    fn a_nested_object_is_walked_into() {
        let fields = editable(SCENE, "/extents/buffer");
        let above = fields
            .iter()
            .find(|e| e.label == "hot_spot.above_k")
            .unwrap_or_else(|| panic!("{:?}", labels(&fields)));
        assert_eq!(above.pointer, "/domains/1/hot_spot/above_k");
        assert_eq!(
            above.value,
            Value::Number {
                now: 60.0,
                integral: false
            }
        );
        // The unit comes from the *last* key. Taking it from the whole label would find nothing,
        // and taking it from the outermost key would call this a `hot_spot`.
        assert_eq!(above.unit, "K");

        let at = fields
            .iter()
            .find(|e| e.label == "hot_spot.at[2]")
            .unwrap_or_else(|| panic!("{:?}", labels(&fields)));
        assert_eq!(at.pointer, "/domains/1/hot_spot/at/2");
    }

    /// **An array of objects is walked into too**, which is how `regions`, `parts` and `held` get
    /// an inspector at all.
    #[test]
    fn an_array_of_objects_is_walked_into() {
        let fields = editable(SCENE, "/extents/buffer");
        let material = fields
            .iter()
            .find(|e| e.label == "regions[0].material")
            .unwrap_or_else(|| panic!("{:?}", labels(&fields)));
        assert_eq!(material.pointer, "/domains/1/regions/0/material");
        let Value::Text { now, choices } = &material.value else {
            panic!("a material is text: {:?}", material.value);
        };
        assert_eq!(now, "copper");
        assert!(choices.contains(&"aluminium".to_string()), "{choices:?}");
    }

    /// **The two identifying keys are not offered**, and the second one is a defect caught before
    /// it shipped rather than a design.
    ///
    /// `kind` because every other key in the object depends on it. `name` because the outliner's
    /// selection is keyed by it: a text field committing per keystroke would rename the domain to
    /// `buffe` on the first backspace, the row `/extents/buffer` would stop existing, and the
    /// field being typed into would vanish. Five other keys point at names as well.
    #[test]
    fn the_identifying_keys_are_not_editable() {
        let fields = labels(&editable(SCENE, "/extents/buffer"));
        for identifier in ["kind", "name"] {
            assert!(
                !fields.iter().any(|l| l == identifier),
                "{identifier} is structure, not a value: {fields:?}"
            );
        }
        // The values beside them still are — this must not have turned the panel off.
        assert!(fields.iter().any(|l| l == "cell_mm"), "{fields:?}");
        assert!(fields.iter().any(|l| l == "material"), "{fields:?}");
    }

    /// **A material menu offers what the scene declared, not only the catalogue.**
    ///
    /// A menu built from `MATERIALS` alone could not express a file that names its own substance,
    /// which is a menu that silently cannot say what the scene already says.
    #[test]
    fn the_material_menu_includes_the_scenes_own_substances() {
        let declared = SCENE.replace(
            r#""domains": ["#,
            r#""materials": { "mystery_alloy": { "density_kg_per_m3": 3000.0,
      "specific_heat_j_per_kg_k": 800.0, "conductivity_w_per_m_k": 50.0 } },
  "domains": ["#,
        );
        let fields = editable(&declared, "/extents/buffer");
        let material = fields
            .iter()
            .find(|e| e.label == "material")
            .unwrap_or_else(|| panic!("{:?}", labels(&fields)));
        let Value::Text { choices, .. } = &material.value else {
            panic!("a material is text");
        };
        assert!(
            choices.contains(&"mystery_alloy".to_string()),
            "{choices:?}"
        );
        assert!(choices.contains(&"copper".to_string()), "{choices:?}");
    }

    /// **A string is spliced like a number, and JSON-escaped on the way in.**
    ///
    /// The one way a text field can break a scene that a number field cannot: a quote or a
    /// backslash written raw would end the string early and the file would stop parsing.
    #[test]
    fn a_string_is_replaced_and_escaped() {
        let out = set_text(SCENE, "/title", r#"a "quoted" title\"#).expect("the pointer resolves");
        let round: serde_json::Value = serde_json::from_str(&out).expect("it still parses");
        assert_eq!(round["title"], r#"a "quoted" title\"#);
        assert!(
            out.contains(r#""duration_s": 0.006"#),
            "nothing else moved: {out}"
        );
    }

    /// **A pointer aimed at a number refuses a string**, rather than quoting it and producing a
    /// scene that fails to load with an error about a type.
    #[test]
    fn a_number_cannot_be_turned_into_a_string() {
        let why = set_text(SCENE, "/frames", "eleven").expect_err("refused");
        assert!(why.contains("not a string"), "{why}");
        let why = set_number(SCENE, "/title", 3.0).expect_err("refused");
        assert!(why.contains("not a number"), "{why}");
    }

    /// **A domain comes out, and the file is still a scene.**
    ///
    /// The comma is the whole difficulty: removing the element alone leaves one behind and the
    /// file stops parsing. Checked from all three positions, because each takes a different side.
    #[test]
    fn a_domain_can_be_removed_from_any_position() {
        let three = r#"{
  "title": "three",
  "duration_s": 1.0,
  "frames": 1,
  "domains": [
    { "kind": "heater", "name": "a", "watts": 20.0, "reserve_j": 1200.0 },
    { "kind": "heater", "name": "b", "watts": 20.0, "reserve_j": 1200.0 },
    { "kind": "heater", "name": "c", "watts": 20.0, "reserve_j": 1200.0 }
  ]
}"#;
        for (gone, left) in [("a", ["b", "c"]), ("b", ["a", "c"]), ("c", ["a", "b"])] {
            let out = remove_domain(three, gone).unwrap_or_else(|e| panic!("{gone}: {e}"));
            let root: serde_json::Value = serde_json::from_str(&out).unwrap_or_else(|e| {
                panic!(
                    "{gone}: {e}
{out}"
                )
            });
            let names: Vec<&str> = root["domains"]
                .as_array()
                .expect("still an array")
                .iter()
                .map(|d| d["name"].as_str().expect("a name"))
                .collect();
            assert_eq!(
                names, left,
                "removing {gone} left {names:?}
{out}"
            );
        }
    }

    /// **The last one out leaves an empty array, not a broken one.**
    #[test]
    fn removing_the_only_domain_leaves_a_scene_that_parses() {
        let one = r#"{ "title": "one", "duration_s": 1.0, "frames": 1,
  "domains": [ { "kind": "heater", "name": "a", "watts": 20.0, "reserve_j": 1200.0 } ] }"#;
        let out = remove_domain(one, "a").expect("the only domain comes out");
        let root: serde_json::Value = serde_json::from_str(&out).expect("and it still parses");
        assert!(
            root["domains"].as_array().expect("an array").is_empty(),
            "{out}"
        );
    }

    /// **Removing something that is not there is refused by name.**
    #[test]
    fn removing_a_domain_that_is_not_there_says_so() {
        let why = remove_domain(SCENE, "ghost").expect_err("refused");
        assert!(why.contains("ghost"), "{why}");
    }

    /// **Every template can be added to a scene, and the result still parses.**
    ///
    /// Nineteen kinds through the splice, which is the check that the insertion point and the
    /// indentation are right for all of them rather than for the one that was tried by hand.
    #[test]
    fn every_kind_can_be_added() {
        for (kind, _) in pantometry_world::templates::TEMPLATES {
            let out = add_domain(SCENE, kind).unwrap_or_else(|e| panic!("{kind}: {e}"));
            let root: serde_json::Value = serde_json::from_str(&out).unwrap_or_else(|e| {
                panic!(
                    "{kind}: {e}
{out}"
                )
            });
            let domains = root["domains"].as_array().expect("an array");
            assert_eq!(domains.len(), 3, "{kind}: the scene had two");
            assert_eq!(domains[2]["kind"].as_str(), Some(kind), "{kind}");
            // And nothing that was there moved.
            assert!(out.contains(r#""cell_mm": 2.0"#), "{kind}: {out}");
        }
    }

    /// **A second domain of the same kind gets a name of its own.**
    ///
    /// The template's name is its kind, so adding two would otherwise make two domains called
    /// `heater` — which the format refuses at build time, from a menu click that looks harmless.
    #[test]
    fn a_repeated_kind_is_numbered_rather_than_colliding() {
        let mut text = SCENE.to_string();
        for expected in ["heater", "heater 2", "heater 3"] {
            text = add_domain(&text, "heater").expect("a heater goes in");
            let root: serde_json::Value = serde_json::from_str(&text).expect("it parses");
            let names: Vec<&str> = root["domains"]
                .as_array()
                .expect("an array")
                .iter()
                .filter_map(|d| d["name"].as_str())
                .collect();
            assert!(names.contains(&expected), "{expected} not in {names:?}");
        }
    }

    /// **Added and removed is the file it started as**, byte for byte.
    ///
    /// The strongest thing either operation can promise, and it holds only if the insertion takes
    /// exactly the bytes the removal gives back — the comma, the newline and the indent included.
    #[test]
    fn adding_then_removing_returns_the_original_text() {
        for (kind, _) in pantometry_world::templates::TEMPLATES {
            let added = add_domain(SCENE, kind).unwrap_or_else(|e| panic!("{kind}: {e}"));
            let back = remove_domain(&added, kind).unwrap_or_else(|e| panic!("{kind}: {e}"));
            assert_eq!(
                back, SCENE,
                "{kind}: the round trip did not land where it started"
            );
        }
    }

    /// **A scene with no `poses` gains one, and keeps everything else.**
    ///
    /// The path no other write takes: three levels missing at once. No shipped scene states
    /// `poses`, so this is what moving anything in a real file actually does.
    #[test]
    fn moving_a_domain_creates_the_key_the_scene_never_had() {
        assert!(
            !SCENE.contains("poses"),
            "the fixture has to start without one"
        );
        let out = set_pose(SCENE, "buffer", [0.01, 0.0, -0.002]).expect("the domain is there");
        assert!(
            out.contains(r#""poses": { "buffer": { "at_m": [0.01, 0.0, -0.002] } }"#),
            "{out}"
        );
        // Everything that was there is still there, byte for byte, and it still builds.
        assert!(out.contains(r#""cell_mm": 2.0"#), "{out}");
        assert!(
            out.contains(r#""kind": "block", "name": "buffer""#),
            "{out}"
        );
        let checked = crate::check(&out, &crate::OnDisk);
        assert!(checked.error.is_none(), "{:?}", checked.error);
        assert_eq!(pose_of(&out, "buffer"), [0.01, 0.0, -0.002]);
    }

    /// **The three levels are each filled in where they are found.**
    ///
    /// A `poses` object with somebody else in it, an entry with no `at_m`, and an `at_m` already
    /// written — the second and third are the branches a scene reaches once it has been moved
    /// once, so they are the common case rather than the exotic one.
    #[test]
    fn each_missing_level_of_a_pose_is_written_where_it_is_missing() {
        let with_other = SCENE.replace(
            r#"  "domains": ["#,
            "  \"poses\": { \"lamp\": { \"at_m\": [1.0, 0.0, 0.0] } },
  \"domains\": [",
        );
        let out = set_pose(&with_other, "buffer", [2.0, 0.0, 0.0]).expect("added beside the other");
        assert_eq!(
            pose_of(&out, "lamp"),
            [1.0, 0.0, 0.0],
            "the other one did not move"
        );
        assert_eq!(pose_of(&out, "buffer"), [2.0, 0.0, 0.0]);

        let empty_entry = SCENE.replace(
            r#"  "domains": ["#,
            "  \"poses\": { \"buffer\": {} },
  \"domains\": [",
        );
        let out = set_pose(&empty_entry, "buffer", [3.0, 0.0, 0.0]).expect("filled the entry");
        assert_eq!(pose_of(&out, "buffer"), [3.0, 0.0, 0.0], "{out}");

        // And moving twice replaces rather than appending a second `at_m`.
        let again = set_pose(&out, "buffer", [4.0, 0.0, 0.0]).expect("moved again");
        assert_eq!(pose_of(&again, "buffer"), [4.0, 0.0, 0.0]);
        assert_eq!(again.matches("at_m").count(), 1, "{again}");
        assert!(crate::check(&again, &crate::OnDisk).error.is_none());
    }

    /// **The box actually moves**, which nothing else here would have noticed.
    ///
    /// Every other pose test asserts the file is valid JSON, that it still checks, and that
    /// `pose_of` reads back what was written. **All of them would pass if `set_pose` wrote to a
    /// key the builder ignores** — the scene would parse, build, draw, and sit exactly where it
    /// started. This is the one that reads the geometry back out of `check` and compares it to
    /// what was asked for, which is the only claim the feature actually makes.
    #[test]
    fn a_posed_domain_is_drawn_where_the_pose_says() {
        let before = crate::check(SCENE, &crate::OnDisk)
            .bounds
            .expect("the block has geometry");
        let moved = set_pose(SCENE, "buffer", [0.5, -0.25, 0.125]).expect("moved");
        let after = crate::check(&moved, &crate::OnDisk)
            .bounds
            .expect("and still has it");

        // The box is the same size and sits exactly the stated offset away. Both halves matter:
        // a pose that scaled or reshaped anything would be a different defect from one that did
        // nothing at all.
        for axis in 0..3 {
            let shift = [0.5, -0.25, 0.125][axis];
            assert!(
                (after[axis] - (before[axis] + shift)).abs() < 1e-12,
                "axis {axis}: low corner went {} -> {}, expected a shift of {shift}",
                before[axis],
                after[axis]
            );
            let span_before = before[axis + 3] - before[axis];
            let span_after = after[axis + 3] - after[axis];
            assert!(
                (span_after - span_before).abs() < 1e-12,
                "axis {axis}: the box changed size, {span_before} -> {span_after}"
            );
        }
    }

    /// **A turn is created where it is missing, and is editable the moment it exists.**
    ///
    /// The second half is the point: nothing was added to the inspector for rotation, because the
    /// generic walk finds `poses.<name>.turn.*` like any other value once the key is there. This
    /// asserts that rather than assuming it — if the walk ever stopped reaching into `poses`, a
    /// rotation control would silently become read-only.
    #[test]
    fn a_turn_is_created_and_then_walks_like_anything_else() {
        let turned = set_turn(SCENE, "buffer", [0.0, 0.0, 1.0], 30.0).expect("the domain is there");
        assert_eq!(turn_of(&turned, "buffer"), ([0.0, 0.0, 1.0], 30.0));
        assert!(crate::check(&turned, &crate::OnDisk).error.is_none());

        let labels: Vec<String> = editable(&turned, "/")
            .into_iter()
            .map(|e| e.label)
            .collect();
        for wanted in [
            "poses.buffer.turn.axis[0]",
            "poses.buffer.turn.axis[2]",
            "poses.buffer.turn.degrees",
        ] {
            assert!(
                labels.iter().any(|l| l == wanted),
                "{wanted} not in {labels:?}"
            );
        }
    }

    /// **A turn goes in beside a position without disturbing it**, at whichever level is missing.
    #[test]
    fn a_turn_and_a_position_coexist() {
        let moved = set_pose(SCENE, "buffer", [0.01, 0.0, 0.0]).expect("moved");
        let both = set_turn(&moved, "buffer", [0.0, 1.0, 0.0], 45.0).expect("turned");
        assert_eq!(
            pose_of(&both, "buffer"),
            [0.01, 0.0, 0.0],
            "the position moved: {both}"
        );
        assert_eq!(turn_of(&both, "buffer"), ([0.0, 1.0, 0.0], 45.0));
        assert_eq!(both.matches("turn").count(), 1, "{both}");

        // And the other order, which takes a different branch at every level.
        let turned = set_turn(SCENE, "buffer", [0.0, 1.0, 0.0], 45.0).expect("turned");
        let both = set_pose(&turned, "buffer", [0.01, 0.0, 0.0]).expect("moved");
        assert_eq!(pose_of(&both, "buffer"), [0.01, 0.0, 0.0]);
        assert_eq!(turn_of(&both, "buffer"), ([0.0, 1.0, 0.0], 45.0));
        assert!(crate::check(&both, &crate::OnDisk).error.is_none());
    }

    /// **The box actually turns**, and a full circle puts it back.
    ///
    /// The same trap the pose tests fell into: every assertion above passes if `set_turn` writes
    /// a key the builder never reads. A full turn returning the original bounds is the closed
    /// form — 360 degrees is the identity for any axis — and a quarter turn moving them is what
    /// says the rotation is applied at all rather than the test being blind to both.
    #[test]
    fn a_turned_domain_is_drawn_turned() {
        let plain = crate::check(SCENE, &crate::OnDisk)
            .bounds
            .expect("geometry");
        let bounds_of = |degrees: f64| {
            let text = set_turn(SCENE, "buffer", [0.0, 0.0, 1.0], degrees).expect("turned");
            crate::check(&text, &crate::OnDisk)
                .bounds
                .expect("still geometry")
        };

        let quarter = bounds_of(90.0);
        assert!(
            plain
                .iter()
                .zip(quarter.iter())
                .any(|(a, b)| (a - b).abs() > 1e-9),
            "a quarter turn about z moved nothing: {plain:?}"
        );

        let full = bounds_of(360.0);
        for (a, b) in plain.iter().zip(full.iter()) {
            assert!(
                (a - b).abs() < 1e-9,
                "a full turn should be the identity: {plain:?} against {full:?}"
            );
        }
    }

    /// **An axis with no direction is refused before the file is touched.**
    ///
    /// `PoseSpec::to_pose` refuses it too, and says why: normalising a zero vector gives a `NaN`
    /// rather than an error, so somebody has to catch it.
    #[test]
    fn a_turn_about_nothing_is_refused() {
        let why = set_turn(SCENE, "buffer", [0.0, 0.0, 0.0], 30.0).expect_err("refused");
        assert!(why.contains("no direction"), "{why}");
        let why = set_turn(SCENE, "buffer", [0.0, 0.0, 1.0], f64::NAN).expect_err("refused");
        assert!(why.contains("no JSON spelling"), "{why}");
        let why = set_turn(SCENE, "ghost", [0.0, 0.0, 1.0], 30.0).expect_err("refused");
        assert!(why.contains("ghost"), "{why}");
    }

    /// **An absent pose reads as the origin**, which is what the format means by silence.
    #[test]
    fn a_domain_with_no_pose_is_at_the_origin() {
        assert_eq!(pose_of(SCENE, "buffer"), [0.0, 0.0, 0.0]);
        assert_eq!(pose_of("not json at all", "buffer"), [0.0, 0.0, 0.0]);
    }

    /// **Posing something that is not in the scene is refused before the file is touched.**
    #[test]
    fn a_pose_for_a_domain_that_is_not_there_is_refused() {
        let why = set_pose(SCENE, "ghost", [1.0, 0.0, 0.0]).expect_err("refused");
        assert!(why.contains("ghost"), "{why}");
        let why = set_pose(SCENE, "buffer", [f64::NAN, 0.0, 0.0]).expect_err("refused");
        assert!(why.contains("no JSON spelling"), "{why}");
    }

    /// **A position keeps its decimal point**, so a whole-number move does not change the
    /// literal's kind under a format that reads three `f64`.
    #[test]
    fn a_whole_number_position_is_still_written_as_a_float() {
        let out = set_pose(SCENE, "buffer", [1.0, 2.0, 3.0]).expect("moved");
        assert!(out.contains("[1.0, 2.0, 3.0]"), "{out}");
    }

    /// A camera and a framing to drag against: the default three-quarter view, on a box one
    /// metre across at the origin.
    fn a_view() -> (viewer_core::Camera, viewer_core::Framing) {
        (
            viewer_core::Camera::default(),
            viewer_core::Framing::of([-0.5, -0.5, -0.5, 0.5, 0.5, 0.5]),
        )
    }

    /// The screen displacement per unit of `axis`, at the origin — the same derivative
    /// `drag_along_axis` takes, by the same short step.
    ///
    /// Deliberately not the secant across a whole unit. That is what the function used to do, and
    /// a helper that repeated the mistake would have agreed with it: the parallel-drag test below
    /// passed the whole time the gizmo was moving things 77% of the way.
    fn on_screen(axis: [f64; 3]) -> [f64; 2] {
        let (camera, frame) = a_view();
        let eps = 1e-4 * frame.span;
        let here = camera.project([0.0; 3], &frame, 1.6);
        let there = camera.project([axis[0] * eps, axis[1] * eps, axis[2] * eps], &frame, 1.6);
        [(there.x - here.x) / eps, (there.y - here.y) / eps]
    }

    /// **A drag exactly along the handle moves exactly one unit**, and one across it moves none.
    ///
    /// Both are exact rather than approximate, because they are what the projection `(d·s)/(s·s)`
    /// *is* — the parallel case is `s·s` over itself and the perpendicular case is a zero dot
    /// product. A gizmo that failed either would be one that drifts off its own axis.
    #[test]
    fn a_drag_along_the_handle_moves_along_the_handle() {
        let (camera, frame) = a_view();
        for axis in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
            let s = on_screen(axis);
            let along = drag_along_axis(&camera, &frame, 1.6, [0.0; 3], axis, s);
            assert!(
                (along - 1.0).abs() < 1e-9,
                "{axis:?}: a drag of one unit's worth moved {along}"
            );
            let across = drag_along_axis(&camera, &frame, 1.6, [0.0; 3], axis, [-s[1], s[0]]);
            assert!(
                across.abs() < 1e-12,
                "{axis:?}: a drag across the handle moved {across}"
            );
        }
    }

    /// **A handle pointing at the camera refuses rather than exploding.**
    ///
    /// Its screen direction is a point, so `s·s` is nearly nothing and the division would turn a
    /// pixel into a leap across the scene. Set up by aiming the camera straight down the axis:
    /// azimuth zero and elevation zero looks along one of them.
    #[test]
    fn an_edge_on_handle_does_not_move() {
        let frame = viewer_core::Framing::of([-0.5, -0.5, -0.5, 0.5, 0.5, 0.5]);
        let camera = viewer_core::Camera {
            azimuth: 0.0,
            elevation: 0.0,
            ..viewer_core::Camera::default()
        };
        // Whichever axis is along the view direction has no screen length. Find it rather than
        // asserting which one it is, because that depends on the rotation's convention.
        let flattest = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            .into_iter()
            .min_by(|a, b| {
                let l = |v: [f64; 3]| {
                    let here = camera.project([0.0; 3], &frame, 1.6);
                    let there = camera.project(v, &frame, 1.6);
                    (there.x - here.x).hypot(there.y - here.y)
                };
                l(*a).total_cmp(&l(*b))
            })
            .expect("three axes");
        let moved = drag_along_axis(&camera, &frame, 1.6, [0.0; 3], flattest, [0.5, 0.5]);
        assert_eq!(moved, 0.0, "{flattest:?} is edge-on and must not move");
    }

    /// **A drag lands where the pointer went, and the error is second order.**
    ///
    /// The projection divides by depth, so the map from screen to world changes as the object
    /// moves and a finite drag cannot be inverted exactly. What can be checked is the *rate*: the
    /// map is the derivative at the starting position, so halving the drag quarters the error.
    ///
    /// **This is the test that found the defect, and no tolerance would have.** The first version
    /// of `drag_along_axis` took its screen direction as the secant across a whole unit of the
    /// axis rather than the derivative, and the object went **77%** of the way the pointer asked
    /// — at every drag size, converging to a fixed 23.5% relative error instead of to zero. An
    /// assertion of "close enough at one drag size" would have been written around it. The order
    /// measured 0.52, 0.80, 0.91, 0.95 — creeping toward 1 where a correct inverse gives 2 — and
    /// that shape is what said the map itself was wrong rather than merely coarse.
    #[test]
    fn a_drag_lands_where_the_pointer_went() {
        let (camera, frame) = a_view();
        let axis = [1.0, 0.0, 0.0];
        let s = on_screen(axis);

        let error_at = |fraction: f64| {
            let wanted = [s[0] * fraction, s[1] * fraction];
            let t = drag_along_axis(&camera, &frame, 1.6, [0.0; 3], axis, wanted);
            let landed = camera.project([axis[0] * t, axis[1] * t, axis[2] * t], &frame, 1.6);
            let here = camera.project([0.0; 3], &frame, 1.6);
            // How far the object actually moved on screen, against how far it was asked to.
            ((landed.x - here.x) - wanted[0]).hypot((landed.y - here.y) - wanted[1])
        };

        let mut previous = error_at(0.4);
        for step in 1..=3 {
            let fraction = 0.4 / 2f64.powi(step);
            let now = error_at(fraction);
            let ratio = previous / now.max(1e-300);
            assert!(
                (3.0..=5.0).contains(&ratio),
                "halving the drag should quarter the error; step {step} gave {ratio:.2}                  ({previous:.3e} -> {now:.3e})"
            );
            previous = now;
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

/// How many states a [`History`] keeps.
///
/// A scene is a few kilobytes and sixty-four of them is a few hundred, which is nothing beside
/// the run this editor holds in memory. The bound exists so a long session cannot grow without
/// one, not because the memory matters.
pub const HISTORY: usize = 64;

/// The sequence of texts the scene has been, so an edit can be taken back.
///
/// # Why the shell cannot own this
///
/// Every edit in this module is a **whole new string**: a splice returns the text it produced
/// rather than mutating one in place. So the history of a scene is a list of strings and undo is
/// an index into it, and the only interesting part is the bookkeeping — which is exactly the
/// kind of thing that is wrong in one direction until somebody writes a test that walks it both
/// ways.
///
/// # The buffer changes two ways and only one of them comes through here
///
/// The shell binds the text pane straight to its buffer, so **typing mutates the scene without
/// passing any function in this crate**. A history that recorded only the splices would take a
/// snapshot from before the typing, and undoing a deleted domain would silently discard whatever
/// had been typed since.
///
/// So [`History::commit`] is called with the current text *before* an edit as well as after it.
/// The first call folds any typing into a state of its own and is a no-op when there was none;
/// the second records the edit. Undo then walks back through the typing rather than over it.
///
/// egui's text widget keeps its own fine-grained undo for what it has focus over, and that is
/// left alone: a focused text field owning its own undo is what every editor does, and fighting
/// it would mean reimplementing keystroke coalescing to no benefit.
#[derive(Clone, Debug)]
pub struct History {
    states: Vec<String>,
    at: usize,
}

impl History {
    /// A history holding one state and nothing to undo to.
    pub fn new(initial: impl Into<String>) -> History {
        History {
            states: vec![initial.into()],
            at: 0,
        }
    }

    /// The text now.
    pub fn current(&self) -> &str {
        &self.states[self.at]
    }

    /// Record that the text is now `next`.
    ///
    /// **Committing the same text again is not a step.** The shell calls this before every edit
    /// to fold in typing, and on a scene nobody typed into that call must not leave an undo that
    /// appears to do nothing — an undo that does nothing is indistinguishable from an undo that
    /// is broken.
    ///
    /// A commit after an undo discards the redo tail, which is what a linear history means: the
    /// future you did not take stops existing the moment you do something else.
    pub fn commit(&mut self, next: impl Into<String>) {
        let next = next.into();
        if next == self.states[self.at] {
            return;
        }
        self.states.truncate(self.at + 1);
        self.states.push(next);
        self.at += 1;
        // Drop from the front, oldest first, and move the cursor with the window. Trimming
        // without moving `at` would leave it pointing at somebody else's state.
        if self.states.len() > HISTORY {
            let over = self.states.len() - HISTORY;
            self.states.drain(..over);
            self.at -= over;
        }
    }

    /// Step back, returning the text to show.
    pub fn undo(&mut self) -> Option<&str> {
        if self.at == 0 {
            return None;
        }
        self.at -= 1;
        Some(&self.states[self.at])
    }

    /// Step forward again.
    pub fn redo(&mut self) -> Option<&str> {
        if self.at + 1 >= self.states.len() {
            return None;
        }
        self.at += 1;
        Some(&self.states[self.at])
    }

    /// Whether there is anything to step back to.
    pub fn can_undo(&self) -> bool {
        self.at > 0
    }

    /// Whether there is anything to step forward to.
    pub fn can_redo(&self) -> bool {
        self.at + 1 < self.states.len()
    }

    /// Start again from `text`, forgetting everything.
    ///
    /// For a **load**, which is not an edit: the file on disk is a different document, and
    /// letting undo walk from one scene back into another would offer a text that belongs to no
    /// file and no history that explains it.
    pub fn reset(&mut self, text: impl Into<String>) {
        self.states = vec![text.into()];
        self.at = 0;
    }

    /// How many states are held, for a test to reason about the bound.
    pub fn len(&self) -> usize {
        self.states.len()
    }

    /// How many steps back are available.
    pub fn steps_back(&self) -> usize {
        self.at
    }

    /// How many steps forward are available.
    pub fn steps_forward(&self) -> usize {
        self.states.len() - self.at - 1
    }

    /// Whether the history holds nothing but its initial state.
    pub fn is_empty(&self) -> bool {
        self.states.len() <= 1
    }
}
