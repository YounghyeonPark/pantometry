//! A run as data: a table for a plot, and a document for a viewer this crate does not contain.
//!
//! The two assets that are not pictures. Both exist because a picture is not always the answer —
//! a researcher with a plotting stack they already trust wants the numbers, and one building
//! their own viewer wants the frames.
//!
//! Written by hand rather than derived from `serde`, and that is deliberate for the JSON: the
//! moment anything reads it, it is a wire format, and a wire format should look like a decision
//! somebody made rather than like whatever the field names happened to be.
//!
//! It is also the **only** asset here that carries a three-dimensional field whole. The two
//! pictures show one slice of one, because a flat canvas cannot do otherwise; this carries every
//! sample, and is what a caller with a volume renderer should take.

use pantometry_scene::{Frame, PanelData};

/// Every domain's scalars, one row per frame.
///
/// **The asset for the domains that have no picture at all** — a heater, a lamp, a winding, a
/// thermal network — and for several of those the scalar *is* the result. It is also the shape
/// a plot wants and the shape a spreadsheet wants.
///
/// Columns are `domain.label` so two networks with a node called `winding` do not collide, and
/// the unit is in the header rather than in a separate legend nobody reads.
pub fn readings_csv(frames: &[Frame]) -> String {
    let Some(first) = frames.first() else {
        return String::from("t_s\n");
    };
    let mut out = String::from("t_s");
    for r in &first.readings {
        out.push_str(&format!(",{}.{} [{}]", r.domain, r.label, r.unit));
    }
    out.push('\n');
    for frame in frames {
        // The time gets the same treatment, and leaving it out was a mistake this file caught an
        // hour later: a 4 ns run printed **every** timestamp as `0.000000000` or `0.000000004`,
        // so the axis a reader plots against had two distinct values in two hundred rows. Fixed
        // notation is only readable at the scale it was chosen for, and a scene format that spans
        // nanoseconds to hours has no such scale.
        out.push_str(&format!("{:.9e}", frame.time_s));
        for r in &frame.readings {
            // **Significant figures, not decimal places.** A fixed `{:.9}` prints every value
            // under a nanounit as `0.000000000`, and a column of zeros that is not zero is the
            // shape this workspace calls a silent failure: measured on a cavity holding 3.2e-10 J,
            // whose entire energy history came out as a column of zeros beside a field the run had
            // just reported at 921 V/m. `{:e}` keeps the magnitude wherever it is, and a
            // spreadsheet reads it.
            out.push_str(&format!(",{:.9e}", r.value));
        }
        out.push('\n');
    }
    out
}

/// What [`to_json`] writes in its `format` field, and the highest a reader of this version
/// understands.
///
/// **One, and it has always been one** - this is the first version *written*, not the first
/// version of the shape. Every run file produced before the key existed is a format 1 file, and
/// a reader treats an absent `format` as 1 for exactly that reason.
///
/// The scene format has carried a version since it had consumers; this one did not, and the gap
/// only mattered once something wanted to add a shape to it. A reader that meets a panel kind it
/// does not know cannot tell "this file is newer than me" from "this file is broken", and those
/// two want different words in front of a person.
pub const FORMAT: u32 = 1;

/// The frames as JSON, for a viewer this crate does not contain.
///
/// Fields as grids, bodies as positions in space, and the readings beside them. Written by hand
/// rather than derived, so the shape is chosen here and stays where a reader can see it — this
/// is a wire format the moment anything consumes it, and it should look deliberate.
pub fn to_json(title: &str, frames: &[Frame]) -> String {
    let mut out = format!(
        "{{\n  \"format\": {FORMAT},\n  \"title\": {},\n  \"frames\": [\n",
        quote(title)
    );
    for (fi, frame) in frames.iter().enumerate() {
        // **Significant figures, not decimal places** — the same defect `readings_csv` above was
        // fixed for and this writer was not. `{:.6}` prints every timestamp of a four-nanosecond
        // run as `0.000000`, so a consumer reading this file back gets two hundred frames all at
        // t = 0 and a frequency read off it is meaningless. Found while sharing the encoder; it
        // had been here since the format was written.
        out.push_str(&format!(
            "    {{ \"t\": {}, \"panels\": [",
            compact(frame.time_s)
        ));
        for (pi, panel) in frame.panels.iter().enumerate() {
            out.push_str(&format!(
                "\n      {{ \"name\": {}, \"unit\": {}, ",
                quote(&panel.name),
                quote(panel.unit)
            ));
            match &panel.data {
                PanelData::Field {
                    nx,
                    ny,
                    nz,
                    extent_m,
                    values,
                } => out.push_str(&format!(
                    "\"kind\": \"field\", \"nx\": {nx}, \"ny\": {ny}, \"nz\": {nz}, \
                     \"extent_m\": {}, \"values\": {}",
                    numbers(extent_m),
                    numbers(values)
                )),
                PanelData::Paths {
                    vertices,
                    starts,
                    values,
                    bounds,
                } => {
                    let flat: Vec<f64> = vertices.iter().flatten().copied().collect();
                    let heads: Vec<f64> = starts.iter().map(|k| *k as f64).collect();
                    out.push_str(&format!(
                        "\"kind\": \"paths\", \"bounds\": {}, \"starts\": {}, \
                         \"vertices\": {}, \"values\": {}",
                        numbers(bounds),
                        numbers(&heads),
                        numbers(&flat),
                        numbers(values)
                    ));
                }
                PanelData::Points {
                    positions,
                    values,
                    bounds,
                    boxed,
                } => {
                    let flat: Vec<f64> = positions.iter().flatten().copied().collect();
                    out.push_str(&format!(
                        "\"kind\": \"points\", \"boxed\": {boxed}, \"bounds\": {}, \
                         \"positions\": {}, \"values\": {}",
                        numbers(bounds),
                        numbers(&flat),
                        numbers(values)
                    ));
                }
            }
            out.push_str(if pi + 1 == frame.panels.len() {
                " }"
            } else {
                " },"
            });
        }
        out.push_str("\n    ], \"readings\": [");
        for (ri, r) in frame.readings.iter().enumerate() {
            out.push_str(&format!(
                "\n      {{ \"domain\": {}, \"label\": {}, \"unit\": {}, \"value\": {} }}{}",
                quote(&r.domain),
                quote(&r.label),
                quote(r.unit),
                compact(r.value),
                if ri + 1 == frame.readings.len() {
                    ""
                } else {
                    ","
                }
            ));
        }
        out.push_str(if fi + 1 == frames.len() {
            "\n    ] }\n"
        } else {
            "\n    ] },\n"
        });
    }
    out.push_str("  ]\n}\n");
    out
}

/// One number, at six significant figures, with nothing spent on saying so.
///
/// Six figures because a picture needs no more. The trailing zeros go because `{:.6e}` writes
/// them whether or not they carry anything: a field of zeros came out as `0.000000e0`, ten
/// characters at a time, and a report is a file somebody has to open.
///
/// **Measured, and smaller than it looks like it should be.** Over the run block of three real
/// reports, against the same block written `{:.6e}`:
///
/// | scene | numbers | saved |
/// | --- | --- | --- |
/// | a cavity ringing, 201 frames of 24x4x24 | 466,320 | 432 kB, **8.3%** |
/// | a room with a ceiling | 61,128 | 7.5 kB, 1.0% |
/// | a world | 27,310 | 5.1 kB, 1.5% |
///
/// The first estimate written here was "5.21 MB to 1.48 MB" and it was a guess, not a
/// measurement. A field of physically varying values *needs* most of its six figures, so the trim
/// takes the odd trailing zero and no more; the cavity gains because its first frames are exactly
/// zero everywhere and `0.000000e0` becomes `0`. Where the data is dense the win is one per cent,
/// and the honest summary is that this makes the largest report noticeably smaller and the rest
/// barely.
///
/// Nothing here changes a value. `3.000000e2` and `3e2` are the same number and parse to the same
/// f64; what is dropped is characters that say nothing about it.
///
/// `null` for anything not finite. `NaN` is not a JSON literal and would make the whole document
/// unreadable to a strict parser; `null` is the spelling for a sample that is not there, and
/// numpy, pandas and `JSON.parse` all take it.
pub(crate) fn compact(x: f64) -> String {
    if !x.is_finite() {
        return "null".to_string();
    }
    if x == 0.0 {
        return "0".to_string();
    }
    let s = format!("{x:.6e}");
    let (mantissa, exponent) = s.split_once('e').unwrap_or((s.as_str(), "0"));
    let mantissa = if mantissa.contains('.') {
        mantissa.trim_end_matches('0').trim_end_matches('.')
    } else {
        mantissa
    };
    if exponent == "0" {
        mantissa.to_string()
    } else {
        format!("{mantissa}e{exponent}")
    }
}

pub(crate) fn numbers(v: &[f64]) -> String {
    let mut s = String::from("[");
    for (i, x) in v.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&compact(*x));
    }
    s.push(']');
    s
}

fn quote(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if c.is_control() => out.push(' '),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}
