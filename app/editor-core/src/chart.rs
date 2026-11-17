//! **Readings over time, for a run that has no picture.**
//!
//! Five shipped scenes draw nothing — three thermal networks and two compartment models — because
//! their domains have no positions to put anything at, and `as_field` and `as_bodies` decline
//! rather than invent one. What they *do* have is readings: every frame carries each
//! compartment's concentration and each node's temperature. The viewport used to say "this scene
//! reports readings, not places" and stop. This is what it draws instead: every reading against
//! time, one plot per unit, a cursor at the frame on screen.
//!
//! It is here rather than in either shell for the reason [`crate::Colouring`] is: two shells
//! drawing one run would otherwise be two chances to put a point somewhere else. A shell asks
//! [`chart`] for the data and [`Chart::draw`] for the strokes and the text, in pixels, and paints
//! them. Neither shell computes a coordinate.
//!
//! # What it promises
//!
//! - **The points are the readings, bit for bit.** A line holds `(t, value)` straight off the
//!   frames, and the drawing maps them; nothing is resampled or smoothed.
//! - **A log axis drops what it cannot place.** A value at or below zero has no logarithm, and
//!   clamping it to the bottom of the axis would draw a cliff the run does not have — a
//!   compartment's `cleared` is exactly zero before the first step. Such a point is left out and
//!   the line breaks there.
//! - **It dispatches on the data**: units group the plots and labels name the lines. Nothing here
//!   names a domain kind.

use viewer_core::Run;

/// One reading's history.
#[derive(Clone, Debug, PartialEq)]
pub struct Line {
    /// The domain that reports it.
    pub domain: String,
    /// What it is called there.
    pub label: String,
    /// `(t in seconds, value)`, one per frame, as the frames hold them.
    pub points: Vec<[f64; 2]>,
}

/// Every line in one unit, on one vertical axis.
#[derive(Clone, Debug, PartialEq)]
pub struct Axis {
    /// The unit, as the readings state it.
    pub unit: String,
    /// The bottom of the axis, in the unit.
    pub lo: f64,
    /// The top of the axis, in the unit.
    pub hi: f64,
    /// The lines drawn against it, in the order their readings first appear.
    pub lines: Vec<Line>,
}

/// What a chart's horizontal axis measures.
///
/// The readings are against time and a profile along a line is against distance. One chart and
/// two axes, rather than a second chart that would be a second chance to place a point wrong.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Across {
    /// Seconds, labelled in whatever unit the span is said in: `time / h`.
    Time,
    /// Metres, labelled the same way: `distance / mm`.
    Distance,
}

/// A run's readings, ready to draw.
#[derive(Clone, Debug, PartialEq)]
pub struct Chart {
    /// Whether the vertical axes are logarithmic.
    pub log: bool,
    /// What the horizontal axis measures.
    pub across: Across,
    /// Where the horizontal axis starts: seconds for time, metres for distance.
    pub x0: f64,
    /// Where it ends, in the same unit.
    pub x1: f64,
    /// One per unit, in the order the units first appear.
    pub axes: Vec<Axis>,
}

/// A run's readings as a chart, or `None` when there is nothing to chart: no frames, or no
/// readings in any of them, or — on a log axis — no positive value in any unit.
pub fn chart(run: &Run, log: bool) -> Option<Chart> {
    let mut axes: Vec<Axis> = Vec::new();
    for frame in &run.frames {
        for r in &frame.readings {
            let a = match axes.iter().position(|a| a.unit == r.unit) {
                Some(a) => a,
                None => {
                    axes.push(Axis {
                        unit: r.unit.clone(),
                        lo: 0.0,
                        hi: 0.0,
                        lines: Vec::new(),
                    });
                    axes.len() - 1
                }
            };
            let lines = &mut axes[a].lines;
            let l = match lines
                .iter()
                .position(|l| l.domain == r.domain && l.label == r.label)
            {
                Some(l) => l,
                None => {
                    lines.push(Line {
                        domain: r.domain.clone(),
                        label: r.label.clone(),
                        points: Vec::new(),
                    });
                    lines.len() - 1
                }
            };
            lines[l].points.push([frame.t, r.value]);
        }
    }
    for axis in &mut axes {
        let values = axis
            .lines
            .iter()
            .flat_map(|l| l.points.iter().map(|p| p[1]))
            .filter(|v| v.is_finite() && (!log || *v > 0.0));
        let (min, max) = values.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
            (lo.min(v), hi.max(v))
        });
        (axis.lo, axis.hi) = if min > max {
            (f64::NAN, f64::NAN)
        } else if log {
            decades(min, max)
        } else {
            let (lo, hi) = (min.min(0.0), max.max(0.0));
            if lo == hi {
                (lo, lo + 1.0)
            } else {
                (lo, hi)
            }
        };
    }
    axes.retain(|a| a.lo.is_finite());
    let x0 = run.frames.first()?.t;
    let x1 = run.frames.last()?.t;
    if axes.is_empty() {
        return None;
    }
    Some(Chart {
        log,
        across: Across::Time,
        x0,
        x1,
        axes,
    })
}

/// Lines in one unit as a chart against `across`, or `None` when no line has a finite point.
///
/// **The vertical axis fits the data**, where [`chart`]'s takes in zero. A reading is usually an
/// amount, and an amount's distance from nothing is part of what it says; a profile through a
/// part at 293 K is not, and an axis from zero would draw it as a flat line along the top of the
/// plot. A twentieth of the range is left above and below, and a line with no range at all gets
/// a thousandth of its value either side, so it is drawn as a line in the middle and not as an
/// axis of zero height.
///
/// **No range means less than a billionth of the value.** A run writes seven figures, so a field
/// that is flat in the file and read between its samples varies only by the interpolation's
/// rounding: scene 24's profile through its module is 488.5735 K to within 1.7e-13. An axis fitted
/// to that is labelled to twelve decimals — scene 25's read `435.018200000000`, 106 points wide,
/// out of the plot and over the outliner beside it. A range that small is not the field's and is
/// drawn as flat.
pub fn graph(across: Across, unit: &str, lines: Vec<Line>) -> Option<Chart> {
    let finite = || {
        lines
            .iter()
            .flat_map(|l| l.points.iter())
            .filter(|p| p[0].is_finite() && p[1].is_finite())
    };
    let (mut x0, mut x1) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for p in finite() {
        x0 = x0.min(p[0]);
        x1 = x1.max(p[0]);
        lo = lo.min(p[1]);
        hi = hi.max(p[1]);
    }
    if lo > hi {
        return None;
    }
    let size = lo.abs().max(hi.abs());
    let pad = if hi - lo > 1e-9 * size {
        0.05 * (hi - lo)
    } else if size > 0.0 {
        // Flat: centred on its middle, so a line rounding-wide is drawn as one line.
        let mid = 0.5 * (lo + hi);
        (lo, hi) = (mid, mid);
        1e-3 * size
    } else {
        1.0
    };
    Some(Chart {
        log: false,
        across,
        x0,
        x1,
        axes: vec![Axis {
            unit: unit.to_string(),
            lo: lo - pad,
            hi: hi + pad,
            lines,
        }],
    })
}

/// The whole decades that hold `min..=max`, both positive: `10^⌊log min⌋` to `10^⌈log max⌉`, and
/// at least one decade apart.
pub(crate) fn decades(min: f64, max: f64) -> (f64, f64) {
    let mut lo = 1.0;
    while lo > min {
        lo /= 10.0;
    }
    while lo * 10.0 <= min {
        lo *= 10.0;
    }
    let mut hi = lo * 10.0;
    while hi < max {
        hi *= 10.0;
    }
    (lo, hi)
}

impl Axis {
    /// Where `v` sits on this axis, 0 at the bottom and 1 at the top, or `None` when it cannot be
    /// placed: not finite, or not positive on a log axis.
    pub fn place(&self, v: f64, log: bool) -> Option<f64> {
        if !v.is_finite() || (log && v <= 0.0) {
            return None;
        }
        Some(if log {
            (v.log10() - self.lo.log10()) / (self.hi.log10() - self.lo.log10())
        } else {
            (v - self.lo) / (self.hi - self.lo)
        })
    }

    /// Where the axis is labelled: each decade on a log axis, and steps of 1, 2 or 5 times a power
    /// of ten on a linear one, about five of them.
    pub fn ticks(&self, log: bool) -> Vec<f64> {
        let mut out = Vec::new();
        // **An axis no wider than its numbers' rounding has one tick, and an axis of no number
        // has none.** A run of one frame is a time axis from `t` to `t`; the step below was a
        // fifth of nothing, and `nice_step` looked for a power of ten under zero until it reached
        // zero and then multiplied zero by ten for ever. And an axis a rounding wide is worse:
        // the results panel's line tab on scene 24 drew a profile flat at 488.5735 K to within
        // 1.7e-13, the step was 5e-14 and the first tick's count 9.77e15 — past 2^53, 9.01e15,
        // where adding one to it changes nothing — and the loop pushed the same tick until it
        // asked for 34 GB and the editor aborted.
        if !(self.lo.is_finite() && self.hi.is_finite()) {
            return out;
        }
        // **And the ends of the doubles.** A fifth of the span is the step's first guess, and
        // `numerics-reviewer` found two axes where it is not a normal number: `[0, 1e-323]`,
        // whose fifth rounds to zero, and `±9.35e307`, whose span overflows to infinity. Both
        // looped for ever in `nice_step` without allocating. They take one tick too.
        let size = self.lo.abs().max(self.hi.abs());
        let fifth = (self.hi - self.lo) / 5.0;
        if self.hi - self.lo <= 64.0 * f64::EPSILON * size || !fifth.is_normal() {
            out.push(self.lo);
            return out;
        }
        if log {
            // A decade is a positive number times ten, and the loop ends at the first one past
            // the top or at infinity, whichever comes first: an axis built by hand with a top
            // near the largest double made `hi * (1 + 1e-9)` infinite and the loop endless.
            // `lo` is finite here, so this is every bottom with no logarithm.
            if self.lo <= 0.0 {
                return out;
            }
            let mut v = self.lo;
            while v.is_finite() && v <= self.hi * (1.0 + 1e-9) {
                out.push(v);
                v *= 10.0;
            }
            return out;
        }
        let step = nice_step(fifth);
        // Past the guard above the span is over 64 ulps of the value, so the step is over a
        // fifth of that and `k` starts below 5 / (64 ε) ≈ 3.5e14 — well under 2^53, where
        // `k += 1.0` still counts.
        let mut k = (self.lo / step).ceil();
        while k * step <= self.hi + step * 1e-9 {
            // `k * step` rather than a running sum, so the tenth tick is not ten roundings out.
            out.push(k * step);
            k += 1.0;
        }
        out
    }
}

/// The step a tick's label is written to: the distance between two ticks, or for an axis with
/// one, the tick's own size.
///
/// **One, it was**, which labels a run of one frame at 0.5 s as `0` — `unearned-pass-hunter` read
/// it off the code, and a check that could not fail had hidden it.
fn step_of(ticks: &[f64]) -> f64 {
    match ticks {
        [a, b, ..] => b - a,
        [v] if *v != 0.0 => v.abs(),
        _ => 1.0,
    }
}

/// One connected run of a line's points, drawn: a stroke for two or more, and a dot for one.
///
/// **A dot, because a stroke of one point is not painted.** The first frame of every run streams
/// in alone, and a line that touches a domain at one point has one value; each drew a legend and
/// nothing beside it.
fn trace(d: &mut Drawing, run: Vec<[f64; 2]>, colour: [u8; 3]) {
    match run.as_slice() {
        [] => {}
        [only] => d.discs.push(Disc {
            centre: *only,
            radius: 2.5,
            fill: Some(colour),
            outline: colour,
        }),
        _ => d.strokes.push(Stroke {
            points: run,
            colour,
            width: 2.0,
        }),
    }
}

/// The smallest of 1, 2 and 5 times a power of ten that is at least `raw`, which must be a
/// positive number: [`Axis::ticks`] does not call it for an axis of no length.
fn nice_step(raw: f64) -> f64 {
    let mut p = 1.0;
    while p > raw {
        p /= 10.0;
    }
    while p * 10.0 <= raw {
        p *= 10.0;
    }
    [1.0, 2.0, 5.0, 10.0]
        .into_iter()
        .map(|m| m * p)
        .find(|s| *s >= raw)
        .unwrap_or(10.0 * p)
}

/// The unit time is labelled in, and seconds in one of it: hours for a run of hours, days for a
/// run of days. A run of four days labelled in seconds is a row of six-digit numbers.
pub fn time_unit(span_s: f64) -> (&'static str, f64) {
    if span_s >= 4.0 * 86_400.0 {
        ("d", 86_400.0)
    } else if span_s >= 2.0 * 3_600.0 {
        ("h", 3_600.0)
    } else if span_s >= 120.0 {
        ("min", 60.0)
    } else {
        ("s", 1.0)
    }
}

/// The unit a length is labelled in, and metres in one of it: millimetres for a part, microns
/// for a channel. Spelled `um`, as [`crate::metres`] spells it.
pub fn length_unit(span_m: f64) -> (&'static str, f64) {
    let a = span_m.abs();
    if a >= 1.0 {
        ("m", 1.0)
    } else if a >= 1e-3 {
        ("mm", 1e-3)
    } else if a >= 1e-6 {
        ("um", 1e-6)
    } else {
        ("nm", 1e-9)
    }
}

/// A line's colour, by its place on its axis: the Okabe–Ito palette, which stays distinct under
/// the three common colour-vision deficiencies, in the order that puts the most distinct first.
pub fn line_colour(i: usize) -> [u8; 3] {
    const PALETTE: [[u8; 3]; 7] = [
        [230, 159, 0],
        [86, 180, 233],
        [0, 158, 115],
        [240, 228, 66],
        [0, 114, 178],
        [213, 94, 0],
        [204, 121, 167],
    ];
    PALETTE[i % PALETTE.len()]
}

/// Where a piece of text is anchored: which of its corners or edges sits on the point.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Anchor {
    /// Its left edge, vertically centred.
    Left,
    /// Its right edge, vertically centred.
    Right,
    /// Its top edge, horizontally centred.
    Top,
    /// Its top-left corner.
    TopLeft,
    /// Its top-right corner.
    TopRight,
    /// Its bottom edge, horizontally centred.
    Bottom,
    /// Its bottom-left corner.
    BottomLeft,
    /// Its bottom-right corner.
    BottomRight,
}

/// A connected run of points, in pixels from the top-left of the drawing.
#[derive(Clone, Debug, PartialEq)]
pub struct Stroke {
    /// The points, in order.
    pub points: Vec<[f64; 2]>,
    /// Its colour.
    pub colour: [u8; 3],
    /// Its width, in pixels.
    pub width: f64,
}

/// A piece of text at a point, in pixels from the top-left of the drawing.
#[derive(Clone, Debug, PartialEq)]
pub struct Text {
    /// Where it goes.
    pub at: [f64; 2],
    /// What it says.
    pub text: String,
    /// Which part of it sits on `at`.
    pub anchor: Anchor,
    /// Its colour.
    pub colour: [u8; 3],
}

/// A filled circle, in pixels from the top-left of the drawing.
#[derive(Clone, Debug, PartialEq)]
pub struct Disc {
    /// Its centre.
    pub centre: [f64; 2],
    /// Its radius, in pixels.
    pub radius: f64,
    /// What it is filled with, or `None` for an outline alone: a value the scale cannot place.
    pub fill: Option<[u8; 3]>,
    /// The colour of its rim.
    pub outline: [u8; 3],
}

/// Everything a chart or a diagram puts on the screen, for a shell to paint as it is.
///
/// **Painted in the order of the fields**: every stroke, then every disc over them, then every
/// text over both. A diagram's vessels run centre to centre and the discs cover their ends.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Drawing {
    /// Lines: the readings, the axes and gridlines, and the cursor.
    pub strokes: Vec<Stroke>,
    /// Circles: a diagram's compartments, and a chart's points that stand alone — a line's only
    /// point, or one between two its axis cannot place.
    pub discs: Vec<Disc>,
    /// Tick labels, units, the time axis and each line's name.
    pub texts: Vec<Text>,
}

impl Drawing {
    /// The same drawing moved down by `dy` pixels, for stacking two in one canvas.
    pub fn shifted(mut self, dy: f64) -> Drawing {
        for s in &mut self.strokes {
            for p in &mut s.points {
                p[1] += dy;
            }
        }
        for c in &mut self.discs {
            c.centre[1] += dy;
        }
        for t in &mut self.texts {
            t.at[1] += dy;
        }
        self
    }

    /// Everything in `other` added after everything here, each kind keeping its painting order.
    pub fn append(&mut self, other: Drawing) {
        self.strokes.extend(other.strokes);
        self.discs.extend(other.discs);
        self.texts.extend(other.texts);
    }
}

/// The grey the frame and tick labels are drawn in, and the fainter one of the gridlines.
pub(crate) const INK: [u8; 3] = [150, 158, 170];
pub(crate) const GRID: [u8; 3] = [58, 64, 74];
/// Room for the tick labels at the left and the time axis below, in pixels.
const LEFT: f64 = 70.0;
const RIGHT: f64 = 14.0;
const TOP: f64 = 28.0;
const BOTTOM: f64 = 34.0;
/// Between two plots, room for the lower one's unit; the first one's is inside `TOP`.
const GAP: f64 = 28.0;

/// A tick's label: a decade as `1e-3`, anything else with as many decimals as the step needs.
fn tick_label(v: f64, step: f64, log: bool) -> String {
    if log {
        let e = v.log10().round() as i32;
        return if (-2..=3).contains(&e) {
            let decimals = (-e).max(0) as usize;
            format!("{v:.decimals$}")
        } else {
            format!("1e{e}")
        };
    }
    // **Zero is `0`**, not twelve decimals of it: the protein's time axis starts at 0 with a step
    // of 5e-12 s, and the decimals the step asks for made its first label `0.000000000000`.
    if v == 0.0 {
        return "0".to_string();
    }
    let a = v.abs();
    if !(1e-3..1e6).contains(&a) {
        // **As many figures as tell two ticks apart.** Two decimals always made the protein's
        // softest mode, flat at 2.1837e11 Hz on an axis stepping 2e8, four labels of `2.18e11`.
        // A tick is a whole number of steps, so the figures between the value's exponent and the
        // step's are exactly what distinguishes it from its neighbours.
        let e = a.log10().floor();
        let s = step.abs().log10().floor();
        let figures = (e - s).clamp(1.0, 12.0) as usize;
        return format!("{v:.figures$e}");
    }
    let mut decimals = 0usize;
    let mut s = step;
    while s < 1.0 - 1e-9 && decimals < 12 {
        s *= 10.0;
        decimals += 1;
    }
    format!("{v:.decimals$}")
}

impl Chart {
    /// Where `x` — seconds or metres, as [`Chart::across`] says — sits across the plots, 0 at
    /// the start of the axis and 1 at its end.
    pub fn place_x(&self, x: f64) -> f64 {
        if self.x1 > self.x0 {
            (x - self.x0) / (self.x1 - self.x0)
        } else {
            0.0
        }
    }

    /// The plot rectangle each axis gets in a `width × height` drawing, as `[left, top, right,
    /// bottom]` in pixels: stacked top to bottom, the same width, sharing one time axis.
    pub fn plots(&self, width: f64, height: f64) -> Vec<[f64; 4]> {
        let n = self.axes.len() as f64;
        let room = (height - TOP - BOTTOM - GAP * (n - 1.0)).max(n);
        let each = room / n;
        (0..self.axes.len())
            .map(|i| {
                let top = TOP + i as f64 * (each + GAP);
                [LEFT, top, (width - RIGHT).max(LEFT + 1.0), top + each]
            })
            .collect()
    }

    /// The chart as strokes and text in a `width × height` drawing, with a cursor at `cursor` —
    /// in the horizontal axis's own unit, seconds or metres — when there is one: the frame on
    /// screen.
    pub fn draw(&self, width: f64, height: f64, cursor: Option<f64>) -> Drawing {
        let mut d = Drawing::default();
        let plots = self.plots(width, height);
        let (quantity, (unit, per)) = match self.across {
            Across::Time => ("time", time_unit(self.x1 - self.x0)),
            Across::Distance => ("distance", length_unit(self.x1 - self.x0)),
        };
        let x = |r: &[f64; 4], t: f64| r[0] + self.place_x(t) * (r[2] - r[0]);
        for (axis, r) in self.axes.iter().zip(&plots) {
            let y = |u: f64| r[3] - u * (r[3] - r[1]);
            let ticks = axis.ticks(self.log);
            let step = step_of(&ticks);
            for v in &ticks {
                let Some(u) = axis.place(*v, self.log) else {
                    continue;
                };
                d.strokes.push(Stroke {
                    points: vec![[r[0], y(u)], [r[2], y(u)]],
                    colour: GRID,
                    width: 1.0,
                });
                d.texts.push(Text {
                    at: [r[0] - 6.0, y(u)],
                    text: tick_label(*v, step, self.log),
                    anchor: Anchor::Right,
                    colour: INK,
                });
            }
            d.strokes.push(Stroke {
                points: vec![[r[0], r[1]], [r[0], r[3]], [r[2], r[3]]],
                colour: INK,
                width: 1.0,
            });
            d.texts.push(Text {
                at: [r[0], r[1] - 20.0],
                text: if self.log {
                    format!("{} (log)", axis.unit)
                } else {
                    axis.unit.clone()
                },
                anchor: Anchor::TopLeft,
                colour: INK,
            });
            for (i, line) in axis.lines.iter().enumerate() {
                let colour = line_colour(i);
                let mut run: Vec<[f64; 2]> = Vec::new();
                for p in &line.points {
                    match axis.place(p[1], self.log) {
                        Some(u) => run.push([x(r, p[0]), y(u)]),
                        None => trace(&mut d, std::mem::take(&mut run), colour),
                    }
                }
                trace(&mut d, run, colour);
                d.texts.push(Text {
                    at: [r[2] - 6.0, r[1] + 4.0 + 15.0 * i as f64],
                    text: format!("{} {}", line.domain, line.label),
                    anchor: Anchor::TopRight,
                    colour,
                });
            }
        }
        if let (Some(first), Some(last)) = (plots.first(), plots.last()) {
            let span = self.x1 - self.x0;
            let axis = Axis {
                unit: String::new(),
                lo: self.x0 / per,
                hi: self.x1 / per,
                lines: Vec::new(),
            };
            let ticks = axis.ticks(false);
            let step = step_of(&ticks);
            for v in ticks {
                let px = x(last, v * per);
                // **A tick at the end of the axis is labelled inside it.** Centred on the plot's
                // right edge, half of `0.006` hung past the drawing's 14-pixel margin and was
                // cut to `0.00` by whatever held the drawing — the results panel's first chart
                // showed it. Right-aligned there, it ends on its tick.
                let near_end = last[2] - px < 24.0;
                d.texts.push(Text {
                    at: [px, last[3] + 6.0],
                    text: tick_label(v, step, false),
                    anchor: if near_end {
                        Anchor::TopRight
                    } else {
                        Anchor::Top
                    },
                    colour: INK,
                });
            }
            d.texts.push(Text {
                at: [last[2], last[3] + 20.0],
                text: format!("{quantity} / {unit}"),
                anchor: Anchor::TopRight,
                colour: INK,
            });
            if let Some(t) = cursor.filter(|_| span > 0.0) {
                let px = x(first, t);
                d.strokes.push(Stroke {
                    points: vec![[px, first[1]], [px, last[3]]],
                    colour: [220, 224, 230],
                    width: 1.0,
                });
            }
        }
        d
    }
}
