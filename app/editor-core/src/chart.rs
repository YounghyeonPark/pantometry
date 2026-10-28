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

/// A run's readings, ready to draw.
#[derive(Clone, Debug, PartialEq)]
pub struct Chart {
    /// Whether the vertical axes are logarithmic.
    pub log: bool,
    /// The time the run starts at, in seconds.
    pub t0: f64,
    /// The time it ends at, in seconds.
    pub t1: f64,
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
    let t0 = run.frames.first()?.t;
    let t1 = run.frames.last()?.t;
    if axes.is_empty() {
        return None;
    }
    Some(Chart { log, t0, t1, axes })
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
        if log {
            let mut v = self.lo;
            while v <= self.hi * (1.0 + 1e-9) {
                out.push(v);
                v *= 10.0;
            }
            return out;
        }
        let step = nice_step((self.hi - self.lo) / 5.0);
        let mut k = (self.lo / step).ceil();
        while k * step <= self.hi + step * 1e-9 {
            // `k * step` rather than a running sum, so the tenth tick is not ten roundings out.
            out.push(k * step);
            k += 1.0;
        }
        out
    }
}

/// The smallest of 1, 2 and 5 times a power of ten that is at least `raw`.
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
    /// Circles: a diagram's compartments. A chart has none.
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
    let mut decimals = 0usize;
    let mut s = step;
    while s < 1.0 - 1e-9 && decimals < 12 {
        s *= 10.0;
        decimals += 1;
    }
    let a = v.abs();
    if a != 0.0 && !(1e-3..1e6).contains(&a) {
        format!("{v:.2e}")
    } else {
        format!("{v:.decimals$}")
    }
}

impl Chart {
    /// Where time `t` (seconds) sits across the plots, 0 at the run's start and 1 at its end.
    pub fn place_time(&self, t: f64) -> f64 {
        if self.t1 > self.t0 {
            (t - self.t0) / (self.t1 - self.t0)
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

    /// The chart as strokes and text in a `width × height` drawing, with a cursor at `cursor`
    /// seconds when there is one — the frame on screen.
    pub fn draw(&self, width: f64, height: f64, cursor: Option<f64>) -> Drawing {
        let mut d = Drawing::default();
        let plots = self.plots(width, height);
        let (unit, per) = time_unit(self.t1 - self.t0);
        let x = |r: &[f64; 4], t: f64| r[0] + self.place_time(t) * (r[2] - r[0]);
        for (axis, r) in self.axes.iter().zip(&plots) {
            let y = |u: f64| r[3] - u * (r[3] - r[1]);
            let ticks = axis.ticks(self.log);
            let step = if ticks.len() > 1 {
                ticks[1] - ticks[0]
            } else {
                1.0
            };
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
                        None => {
                            if !run.is_empty() {
                                d.strokes.push(Stroke {
                                    points: std::mem::take(&mut run),
                                    colour,
                                    width: 2.0,
                                });
                            }
                        }
                    }
                }
                if !run.is_empty() {
                    d.strokes.push(Stroke {
                        points: run,
                        colour,
                        width: 2.0,
                    });
                }
                d.texts.push(Text {
                    at: [r[2] - 6.0, r[1] + 4.0 + 15.0 * i as f64],
                    text: format!("{} {}", line.domain, line.label),
                    anchor: Anchor::TopRight,
                    colour,
                });
            }
        }
        if let (Some(first), Some(last)) = (plots.first(), plots.last()) {
            let span = self.t1 - self.t0;
            let axis = Axis {
                unit: String::new(),
                lo: self.t0 / per,
                hi: self.t1 / per,
                lines: Vec::new(),
            };
            let ticks = axis.ticks(false);
            let step = if ticks.len() > 1 {
                ticks[1] - ticks[0]
            } else {
                1.0
            };
            for v in ticks {
                let px = x(last, v * per);
                d.texts.push(Text {
                    at: [px, last[3] + 6.0],
                    text: tick_label(v, step, false),
                    anchor: Anchor::Top,
                    colour: INK,
                });
            }
            d.texts.push(Text {
                at: [last[2], last[3] + 20.0],
                text: format!("time / {unit}"),
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
