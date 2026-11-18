//! **A parametric sweep: one value of a scene, several runs, and the readings against it.**
//!
//! A design question is rarely "what does this scene do". It is "what does it do as this number
//! changes" — the film coefficient a heatsink needs, the width at which a well's level crosses an
//! energy, the current at which a busbar runs too hot. A simulation tool answers that by running the
//! model once per value and tabulating what came out, and that is what this is:
//!
//! - **which values** — every number the scene states, by the same listing the settings form uses
//!   ([`parameters`]), so anything a person can type a value into they can sweep;
//! - **the variants** — the text with that one literal replaced, by the byte splice every other
//!   edit uses ([`variant`]), so each run is of the file as written with one thing changed, and a
//!   whole-number literal refuses a fraction rather than being rounded into a different question;
//! - **what came out** — every reading of each run's last frame ([`Point`]), the state the run
//!   reached, which is the answer a sweep is usually asking for;
//! - **the picture** — one reading against the value, through the same chart the readings use.
//!
//! No thread and no window here, so the editor, the CLI and the browser make the same calls. The
//! runs are the run button's, so a sweep point is the run a person would get by typing that value
//! in and pressing run, which `the_scenes_own_value_is_the_scenes_own_run` holds — read off its
//! frames at full precision rather than back out of the run file's seven figures.

use crate::chart::{self, Across, Line};
use crate::edit;

/// One number the scene states, which a sweep can vary.
#[derive(Clone, Debug, PartialEq)]
pub struct Parameter {
    /// Where it is: an RFC 6901 pointer, `/domains/0/width_nm`.
    pub pointer: String,
    /// What it is called: `particle.width_nm`, or `duration_s` for the scene's own.
    pub label: String,
    /// The unit its key names, or empty.
    pub unit: &'static str,
    /// What it is now.
    pub now: f64,
    /// Whether the file writes it as a whole number, which a sweep must keep it.
    pub integral: bool,
}

/// Every number the scene states, the scene's own first and then each domain's in file order.
///
/// The settings form's own listing, so a value a person can type into is a value they can sweep,
/// and one the form holds back — a domain's `kind`, its `name` — is not. Empty for text that does not
/// parse.
pub fn parameters(text: &str) -> Vec<Parameter> {
    let mut out = Vec::new();
    let as_parameter = |e: edit::Editable, prefix: &str| match e.value {
        edit::Value::Number { now, integral } => Some(Parameter {
            pointer: e.pointer,
            label: if prefix.is_empty() {
                e.label
            } else {
                format!("{prefix}.{}", e.label)
            },
            unit: e.unit,
            now,
            integral,
        }),
        edit::Value::Text { .. } => None,
    };
    out.extend(
        edit::editable(text, "/")
            .into_iter()
            .filter_map(|e| as_parameter(e, "")),
    );
    let Ok(root) = serde_json::from_str::<serde_json::Value>(text) else {
        return out;
    };
    let names: Vec<String> = root
        .get("domains")
        .and_then(|d| d.as_array())
        .into_iter()
        .flatten()
        .filter_map(|d| d.get("name").and_then(|n| n.as_str()).map(str::to_string))
        .collect();
    for name in names {
        out.extend(
            edit::editable(text, &format!("/domains/{name}"))
                .into_iter()
                .filter_map(|e| as_parameter(e, &name)),
        );
    }
    out
}

/// `count` values from `from` to `to`, both included, evenly spaced — or evenly in their logarithm
/// when `log`, which is how a film coefficient or a cell size is usually swept.
///
/// Refused for no values, for a range that is not two finite numbers, and for a logarithmic one
/// that reaches zero or below. A whole-number parameter takes whole numbers: its ends must be whole,
/// as a fraction written into it is refused, and the values between are rounded to them and a value
/// that two steps round to kept once, so asking for ten values between 2 and 5 gives four and does
/// not run three of them twice.
///
/// **The kernel's `ln` and `exp`**, not the platform's. A value here becomes a literal in the scene,
/// so an ulp's difference between Windows, Linux and WebAssembly would be a different run of the
/// same sweep; the kernel's are the same bits everywhere, which is what they are for.
pub fn spaced(
    from: f64,
    to: f64,
    count: usize,
    log: bool,
    integral: bool,
) -> Result<Vec<f64>, String> {
    if count == 0 {
        return Err("a sweep of no values runs nothing".to_string());
    }
    if !from.is_finite() || !to.is_finite() {
        return Err(format!("{from} to {to} is not a range"));
    }
    if log && (from <= 0.0 || to <= 0.0) {
        return Err(format!(
            "a logarithmic sweep needs both ends above zero, and this one is {from} to {to}"
        ));
    }
    if integral && (from.fract() != 0.0 || to.fract() != 0.0) {
        return Err(format!(
            "a whole number is swept between whole numbers, and this is {from} to {to}"
        ));
    }
    use pantometry::core::math::{exp, ln};
    let mut out: Vec<f64> = (0..count)
        .map(|i| {
            if count == 1 {
                return from;
            }
            let f = i as f64 / (count - 1) as f64;
            // The ends exactly, so a sweep from 5 to 20 runs 5 and 20 and not 19.999999999999996.
            if i == 0 {
                from
            } else if i + 1 == count {
                to
            } else if log {
                exp(ln(from) + f * (ln(to) - ln(from)))
            } else {
                from + f * (to - from)
            }
        })
        .collect();
    if integral {
        for v in &mut out {
            *v = v.round();
        }
        out.dedup();
    }
    Ok(out)
}

/// Values written out, `5,10,20`, as a person types them after `=`.
pub fn listed(said: &str) -> Result<Vec<f64>, String> {
    let values: Vec<f64> = said
        .split(',')
        .map(|v| {
            v.trim()
                .parse::<f64>()
                .map_err(|_| format!("{v:?} in {said:?} is not a number"))
        })
        .collect::<Result<_, _>>()?;
    if values.iter().any(|v| !v.is_finite()) {
        return Err(format!("{said:?} holds a value that is not finite"));
    }
    if values.is_empty() {
        return Err("a sweep of no values runs nothing".to_string());
    }
    Ok(values)
}

/// A sweep as a person writes one: `particle.width_nm=5,10,20`, or `…=5:20:4`, or
/// `…=1:1000:4:log` — the number by the label [`parameters`] gives it or by its pointer, and the
/// values written out or spaced. The CLI, the editor's window and `--ui-dump` all read this one
/// spelling.
///
/// Refused, with what to write instead, for a number the scene does not state, values that do not
/// parse, and a fraction for a number the file writes whole — said once here rather than once per
/// point.
pub fn asked(text: &str, spec: &str) -> Result<(Parameter, Vec<f64>), String> {
    let Some((name, said)) = spec.split_once('=') else {
        return Err(format!(
            "{spec} names no values: write NAME=5,10,20 or NAME=from:to:count[:log]"
        ));
    };
    let Some(p) = parameters(text)
        .into_iter()
        .find(|p| p.label == name || p.pointer == name)
    else {
        return Err(format!(
            "the scene has no number called {name}; `pantometry sweep scene.json --list` names them"
        ));
    };
    let parts: Vec<&str> = said.split(':').collect();
    let values = match parts[..] {
        [from, to, count] | [from, to, count, _] => {
            let log = match parts.get(3) {
                None => false,
                Some(&"log") => true,
                Some(other) => {
                    return Err(format!(
                        "{other:?} in {said:?}: the fourth part is `log` or nothing"
                    ))
                }
            };
            let number = |s: &str| {
                s.parse::<f64>()
                    .map_err(|_| format!("{s:?} in {said:?} is not a number"))
            };
            let count = count
                .parse::<usize>()
                .map_err(|_| format!("{count:?} in {said:?} is not a count"))?;
            spaced(number(from)?, number(to)?, count, log, p.integral)?
        }
        [_] => listed(said)?,
        _ => {
            return Err(format!(
                "{said:?}: values are 5,10,20 or from:to:count[:log]"
            ))
        }
    };
    if let Some(v) = values.iter().find(|v| p.integral && v.fract() != 0.0) {
        return Err(format!(
            "{} is a whole number in this file and {v} is not one",
            p.label
        ));
    }
    Ok((p, values))
}

/// The scene with the number at `pointer` set to `value`, every other byte as it was.
pub fn variant(text: &str, pointer: &str, value: f64) -> Result<String, String> {
    edit::set_number(text, pointer, value)
}

/// One reading of a run's last frame.
#[derive(Clone, Debug, PartialEq)]
pub struct Reading {
    /// The domain that reported it.
    pub domain: String,
    /// What it is.
    pub label: String,
    /// Its unit.
    pub unit: String,
    /// Its value.
    pub value: f64,
}

/// What one value of a sweep came to: the readings of its run's last frame, or why there were
/// none.
#[derive(Clone, Debug, PartialEq)]
pub struct Point {
    /// The value the parameter was set to.
    pub value: f64,
    /// The last frame's readings, in the order the frame lists them. Empty when the run refused.
    pub readings: Vec<Reading>,
    /// Why there are no readings: the edit refused, the scene did not build, or the audit
    /// stopped the run. `None` when it ran.
    pub error: Option<String>,
}

/// Run the scene with `pointer` at `value` and keep its last frame's readings.
///
/// The run is the one [`crate::run`] writes out, so a point is the run a person would get by typing
/// the value in; its readings are taken before they are written, at full precision. A refusal is
/// kept beside the value rather than ending the sweep: a value where the scene stops building is
/// an answer about the design.
pub fn run_point(text: &str, pointer: &str, value: f64, files: &dyn crate::Parts) -> Point {
    let fail = |error: String| Point {
        value,
        readings: Vec::new(),
        error: Some(error),
    };
    let scene = match variant(text, pointer, value) {
        Ok(s) => s,
        Err(e) => return fail(e),
    };
    // **A panic is a refusal too.** A sweep hands the scene values nobody wrote, and a run that
    // panicked on one took its thread down with nothing sent: the editor sat at "sweeping: k of n"
    // with Stop doing nothing. Caught here, it is that value's answer, and the sweep goes on.
    let ran = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        crate::run_frames(&scene, files)
    }));
    let frames = match ran {
        Ok(Ok((_, f))) => f,
        Ok(Err(e)) => return fail(e),
        Err(panic) => {
            let why = panic
                .downcast_ref::<String>()
                .map(String::as_str)
                .or_else(|| panic.downcast_ref::<&str>().copied())
                .unwrap_or("no message");
            return fail(format!("the run panicked: {why}"));
        }
    };
    let Some(last) = frames.last() else {
        return fail("the run wrote no frames".to_string());
    };
    Point {
        value,
        readings: last
            .readings
            .iter()
            .map(|r| Reading {
                domain: r.domain.clone(),
                label: r.label.clone(),
                unit: r.unit.to_string(),
                value: r.value,
            })
            .collect(),
        error: None,
    }
}

/// Every value in turn, handing each point to `each` as it arrives, and stopping between runs when
/// `stop` is raised. The points run, in order: a stopped sweep is a prefix, and says nothing it did
/// not run.
pub fn run_all(
    text: &str,
    pointer: &str,
    values: &[f64],
    files: &dyn crate::Parts,
    stop: &std::sync::atomic::AtomicBool,
    mut each: impl FnMut(usize, &Point),
) -> Vec<Point> {
    let mut out = Vec::new();
    for (i, v) in values.iter().enumerate() {
        if stop.load(std::sync::atomic::Ordering::Relaxed) {
            break;
        }
        let p = run_point(text, pointer, *v, files);
        each(i, &p);
        out.push(p);
    }
    out
}

/// One reading a sweep tabulates: which domain, what it is called, its unit, and which of the
/// readings of that name it is when a frame reports more than one.
#[derive(Clone, Debug, PartialEq)]
pub struct Column {
    /// The domain that reports it.
    pub domain: String,
    /// What it is called.
    pub label: String,
    /// Its unit.
    pub unit: String,
    /// Which of a frame's readings with this domain and label it is, counting from zero: almost
    /// always 0, and more only where a domain reports one name twice, which `verify` flags.
    pub nth: usize,
    /// How many of the sweep's points report it.
    pub reported: usize,
}

impl Column {
    /// What a table heads it and a menu offers it as: `particle <E> / J`.
    pub fn name(&self) -> String {
        let nth = if self.nth > 0 {
            format!(" ({})", self.nth + 1)
        } else {
            String::new()
        };
        if self.unit.is_empty() {
            format!("{} {}{nth}", self.domain, self.label)
        } else {
            format!("{} {}{nth} / {}", self.domain, self.label, self.unit)
        }
    }
}

/// The readings any point reports, in the order they first appear, each by its domain, its label
/// and which of that name it is.
///
/// **From every point, by name.** A frame need not report the same readings at every value: a
/// block reports `generated` only while its source is on, so on scene 25 a sweep of the source's
/// watts from 0 has one more reading at 45 W than at 0, and every reading after it sits one place
/// further down. Taken from the first point and matched by place, every `stress` reading at 45 W was
/// blank. By name, each lands in its own column, and [`Column::reported`] says how many points have
/// one.
pub fn columns(points: &[Point]) -> Vec<Column> {
    let mut out: Vec<Column> = Vec::new();
    for p in points.iter().filter(|p| p.error.is_none()) {
        let mut seen: Vec<(&str, &str)> = Vec::new();
        for r in &p.readings {
            let nth = seen
                .iter()
                .filter(|(d, l)| *d == r.domain && *l == r.label)
                .count();
            seen.push((&r.domain, &r.label));
            match out
                .iter_mut()
                .find(|c| c.domain == r.domain && c.label == r.label && c.nth == nth)
            {
                Some(c) => c.reported += 1,
                None => out.push(Column {
                    domain: r.domain.clone(),
                    label: r.label.clone(),
                    unit: r.unit.clone(),
                    nth,
                    reported: 1,
                }),
            }
        }
    }
    out
}

/// Point `p`'s value in column `k`: its reading of that domain and label, the column's `nth` of
/// them, wherever in the frame it is.
pub fn at(p: &Point, columns: &[Column], k: usize) -> Option<f64> {
    let c = columns.get(k)?;
    p.readings
        .iter()
        .filter(|r| r.domain == c.domain && r.label == c.label)
        .nth(c.nth)
        .map(|r| r.value)
}

/// How many points ran.
pub fn ran(points: &[Point]) -> usize {
    points.iter().filter(|p| p.error.is_none()).count()
}

/// The first column whose values differ across the points, to open a table and a chart on: a
/// column that holds one number at every value — a norm, a reserve — is the least interesting one.
///
/// A column moves when its spread is more than a billionth of its largest magnitude, which a
/// column of rounding around one number does not reach. `None` when nothing moves, or before two
/// points have run, when nothing can have.
pub fn moving(points: &[Point], columns: &[Column]) -> Option<usize> {
    (0..columns.len()).find(|&k| {
        let v: Vec<f64> = points.iter().filter_map(|p| at(p, columns, k)).collect();
        let (lo, hi) = v
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
                (lo.min(*x), hi.max(*x))
            });
        let size = v.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        v.len() >= 2 && size > 0.0 && hi - lo > 1e-9 * size
    })
}

/// The sweep as CSV: the parameter, then every reading, one row a value; a row that did not run
/// says why in the last column, `error`, and a cell a run did not report is empty.
pub fn csv(parameter: &str, points: &[Point]) -> String {
    let cols = columns(points);
    let quote = |s: &str| {
        if s.contains([',', '"', '\n', '\r']) {
            format!("\"{}\"", s.replace('"', "\"\""))
        } else {
            s.to_string()
        }
    };
    let mut out = quote(parameter);
    for c in &cols {
        out.push(',');
        let nth = if c.nth > 0 {
            format!(" ({})", c.nth + 1)
        } else {
            String::new()
        };
        out.push_str(&quote(&if c.unit.is_empty() {
            format!("{}.{}{nth}", c.domain, c.label)
        } else {
            format!("{}.{}{nth} [{}]", c.domain, c.label, c.unit)
        }));
    }
    out.push_str(",error\n");
    for p in points {
        out.push_str(&format!("{:e}", p.value));
        for k in 0..cols.len() {
            out.push(',');
            if let Some(v) = at(p, &cols, k) {
                out.push_str(&format!("{v:e}"));
            }
        }
        out.push(',');
        out.push_str(&quote(p.error.as_deref().unwrap_or("")));
        out.push('\n');
    }
    out
}

/// Column `k` against the parameter, as a chart, or `None` when no point has a value there.
///
/// **In the value's order, and broken where a value has no reading.** The points are in the order
/// they ran, which for values written out is the order they were written — `20, 5, 10` drew a line
/// from right to left and back. And a value the scene refused would have been joined across: the
/// line now stops at it, so nothing is drawn where nothing ran.
pub fn chart(parameter: &Parameter, points: &[Point], k: usize) -> Option<chart::Chart> {
    let cols = columns(points);
    let c = cols.get(k)?.clone();
    let mut drawn: Vec<[f64; 2]> = points
        .iter()
        .map(|p| [p.value, at(p, &cols, k).unwrap_or(f64::NAN)])
        .collect();
    drawn.sort_by(|a, b| a[0].total_cmp(&b[0]));
    if drawn.iter().all(|p| !p[1].is_finite()) {
        return None;
    }
    chart::graph(
        Across::Stated {
            name: parameter.label.clone(),
            unit: parameter.unit.to_string(),
        },
        &c.unit,
        vec![Line {
            domain: c.domain,
            label: c.label,
            points: drawn,
        }],
    )
}
