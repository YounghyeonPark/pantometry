//! **A compartment model drawn as the body it describes.**
//!
//! A `compartments` domain has no geometry, and the run says so: it offers no field and no bodies,
//! only readings. [`crate::chart`] draws those against time, which is where a rate is read. What a
//! chart does not show is *where the drug is*: which pool feeds which, and how the dose spreads
//! out of the blood into the liver, the lean tissue and the fat and comes back. This draws that —
//! every compartment as a disc, every perfusion and link as a vessel, every depot as the route a
//! dose arrives by — coloured by the concentration in the frame on screen.
//!
//! # Where its shape comes from, and why not from the run
//!
//! The run carries concentrations and nothing about what is connected to what; the topology is in
//! the scene. So this reads `DomainSpec::Compartments` from the parsed scene, which the editor is
//! allowed to do — ARCHITECTURE.md's authoring half is the composition root and may name domains,
//! as [`crate::check`] already does to lay out boxes — and takes every value from the run's
//! readings, matched by domain and label.
//!
//! **No position here reaches any physics.** The library's `Placement::marker` was written for a
//! diagram like this, and does not fit it: it is one position per *domain*, and the four
//! compartments of scene 34 are one domain. The layout is computed instead: the first compartment —
//! the central one, by the format's own rule — at the middle, the rest on a ring around it, each
//! depot outside the compartment it empties into.
//!
//! # What it promises
//!
//! - **Every compartment is a disc and every vessel a line**, one each, in the scene's order. A
//!   compartment with no reading in the frame is drawn as an outline labelled `no reading`, not
//!   left out.
//! - **A disc's area is its volume.** The radius goes as `√V`, so a reader comparing two discs is
//!   comparing two volumes; a vessel's width goes with its flow.
//! - **A disc's colour is its concentration on one scale for the whole run**: whole decades from
//!   the lowest positive concentration any compartment reaches to the highest, logarithmic,
//!   through the library's sequential ramp. A scale renormalised per frame would make a body that
//!   is emptying look the same at every frame. Zero has no logarithm and is drawn as an outline.

use crate::chart::{self, Anchor, Disc, Drawing, Stroke, Text, GRID, INK};
use pantometry_world::{DomainSpec, Scene};
use viewer_core::{Frame, Run};

/// The unit a compartment's concentration is reported in.
const CONCENTRATION: &str = "mg/L";
/// The radius of the largest compartment, in the layout's units, where the ring is radius one.
const R_MAX: f64 = 0.42;
/// A depot's radius, and how far outside its compartment's rim it sits, in the same units.
const DEPOT_R: f64 = 0.07;
const DEPOT_REACH: f64 = 0.55;
/// How far a clearance arrow reaches out of its compartment's rim.
const ARROW: f64 = 0.32;
/// The colour a compartment's name is written in: darker than [`INK`], and readable on the light
/// viewport and the dark one alike.
const NAME: [u8; 3] = [112, 120, 134];
/// The grey of a vessel.
const VESSEL: [u8; 3] = [96, 106, 120];
/// Room around the layout for the labels, and at the right for the colour bar, in pixels.
const PAD_X: f64 = 76.0;
const PAD_TOP: f64 = 44.0;
const PAD_BOTTOM: f64 = 26.0;
const BAR: f64 = 84.0;
/// How many turns of the layout are tried, to find the one the rectangle holds largest.
const TURNS: usize = 24;

/// One compartment, placed.
#[derive(Clone, Debug, PartialEq)]
pub struct Organ {
    /// Its label, which is the label its concentration is reported under.
    pub label: String,
    /// Its volume, in litres.
    pub volume_l: f64,
    /// Its clearance out of the body, in litres per hour; zero when it does not eliminate.
    pub clearance_l_per_h: f64,
    /// Its centre, in the layout's units: the central compartment at the origin, the rest on the
    /// unit circle, `y` down as on a screen.
    pub at: [f64; 2],
    /// Its radius in the same units, `R_MAX·√(V/V_max)`.
    pub radius: f64,
    /// Which way its clearance arrow points, a unit vector.
    pub clears_towards: [f64; 2],
}

/// What a vessel carries.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Carries {
    /// Blood flow from a pool through a tissue: `from` is the blood, `to` the tissue.
    Perfusion,
    /// An intercompartmental clearance, symmetric.
    Link,
}

/// A perfusion or a link between two compartments.
#[derive(Clone, Debug, PartialEq)]
pub struct Vessel {
    /// Index into [`Body::organs`]: the blood pool of a perfusion, either end of a link.
    pub from: usize,
    /// The other end: the tissue of a perfusion.
    pub to: usize,
    /// The flow of a perfusion or the clearance of a link, in litres per hour.
    pub flow_l_per_h: f64,
    /// Which of the two it is.
    pub carries: Carries,
}

/// An absorption route, placed.
#[derive(Clone, Debug, PartialEq)]
pub struct Depot {
    /// Its label.
    pub label: String,
    /// Index into [`Body::organs`] of the compartment it empties into.
    pub into: usize,
    /// Its centre, in the layout's units.
    pub at: [f64; 2],
}

/// One `compartments` domain of a scene, laid out to draw.
#[derive(Clone, Debug, PartialEq)]
pub struct Body {
    /// The domain's name, which is the `domain` of its readings.
    pub domain: String,
    /// Every compartment, in the scene's order.
    pub organs: Vec<Organ>,
    /// Every perfusion, then every link, in the scene's order. One whose labels name no
    /// compartment is not drawn; the build refuses such a scene and says why.
    pub vessels: Vec<Vessel>,
    /// Every depot, in the scene's order.
    pub depots: Vec<Depot>,
}

/// Every `compartments` domain in `scene`, laid out. Empty for a scene with none.
pub fn bodies(scene: &Scene) -> Vec<Body> {
    let mut out = Vec::new();
    for spec in &scene.domains {
        let DomainSpec::Compartments {
            name,
            volumes,
            links,
            perfusions,
            depots,
            ..
        } = spec
        else {
            continue;
        };
        if volumes.is_empty() {
            continue;
        }
        let index = |label: &str| volumes.iter().position(|c| c.label == label);
        let ring = volumes.len() - 1;
        let angles: Vec<f64> = (0..ring)
            .map(|k| std::f64::consts::PI * (1.0 + 2.0 * k as f64 / ring as f64))
            .collect();
        let largest = volumes.iter().map(|c| c.volume_l).fold(0.0_f64, f64::max);
        let mut organs: Vec<Organ> = volumes
            .iter()
            .enumerate()
            .map(|(i, c)| {
                let at = if i == 0 {
                    [0.0, 0.0]
                } else {
                    let a = angles[i - 1];
                    [a.cos(), a.sin()]
                };
                let radius = if largest > 0.0 && c.volume_l > 0.0 {
                    R_MAX * (c.volume_l / largest).sqrt()
                } else {
                    0.0
                };
                Organ {
                    label: c.label.clone(),
                    volume_l: c.volume_l,
                    clearance_l_per_h: c.clearance_l_per_h,
                    at,
                    radius,
                    clears_towards: [1.0, 0.0],
                }
            })
            .collect();
        // Outward along its own ray for a compartment on the ring; into the widest gap between
        // the ring's compartments for the central one, which has no ray of its own.
        let gap = widest_gap(&angles);
        let fed: Vec<usize> = depots.iter().filter_map(|d| index(&d.into)).collect();
        for (i, organ) in organs.iter_mut().enumerate() {
            let base = if i == 0 { gap } else { angles[i - 1] };
            // A depot takes the ray, so the clearance arrow turns off it rather than through it.
            let turn = if fed.contains(&i) { 0.96 } else { 0.0 };
            organ.clears_towards = [(base + turn).cos(), (base + turn).sin()];
        }
        let mut vessels = Vec::new();
        for p in perfusions {
            if let (Some(from), Some(to)) = (index(&p.blood), index(&p.tissue)) {
                vessels.push(Vessel {
                    from,
                    to,
                    flow_l_per_h: p.flow_l_per_h,
                    carries: Carries::Perfusion,
                });
            }
        }
        for l in links {
            if let (Some(from), Some(to)) = (index(&l.between[0]), index(&l.between[1])) {
                vessels.push(Vessel {
                    from,
                    to,
                    flow_l_per_h: l.clearance_l_per_h,
                    carries: Carries::Link,
                });
            }
        }
        let depots = depots
            .iter()
            .filter_map(|d| {
                let into = index(&d.into)?;
                let o = &organs[into];
                let a = if into == 0 { gap } else { angles[into - 1] };
                let reach = o.radius + DEPOT_REACH;
                Some(Depot {
                    label: d.label.clone(),
                    into,
                    at: [o.at[0] + reach * a.cos(), o.at[1] + reach * a.sin()],
                })
            })
            .collect();
        out.push(Body {
            domain: name.clone(),
            organs,
            vessels,
            depots,
        });
    }
    out
}

/// The middle of the widest angular gap between `angles`, which are ascending in `[π, 3π)`; to
/// the right when there are none.
fn widest_gap(angles: &[f64]) -> f64 {
    if angles.is_empty() {
        return 0.0;
    }
    let tau = 2.0 * std::f64::consts::PI;
    let mut best = (f64::NEG_INFINITY, 0.0);
    for (k, a) in angles.iter().enumerate() {
        let next = angles.get(k + 1).copied().unwrap_or(angles[0] + tau);
        if next - a > best.0 {
            best = (next - a, a + 0.5 * (next - a));
        }
    }
    best.1
}

impl Body {
    /// Each compartment's concentration in `frame`, in mg/L, in the order of [`Body::organs`]:
    /// the reading under this domain and the compartment's label, or `None` when the frame has
    /// none.
    pub fn concentrations(&self, frame: &Frame) -> Vec<Option<f64>> {
        self.organs
            .iter()
            .map(|o| {
                frame
                    .readings
                    .iter()
                    .find(|r| {
                        r.domain == self.domain && r.label == o.label && r.unit == CONCENTRATION
                    })
                    .map(|r| r.value)
            })
            .collect()
    }

    /// The colour scale for the whole run, `(lo, hi)` in mg/L: the whole decades holding every
    /// positive concentration any compartment reaches in any frame. `None` when none is positive.
    pub fn scale(&self, run: &Run) -> Option<(f64, f64)> {
        let (min, max) = run
            .frames
            .iter()
            .flat_map(|f| self.concentrations(f))
            .flatten()
            .filter(|v| v.is_finite() && *v > 0.0)
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
                (lo.min(v), hi.max(v))
            });
        (min <= max).then(|| chart::decades(min, max))
    }

    /// The colour of concentration `v` on `scale`, or `None` when it cannot be placed: not
    /// positive, not finite, or no scale.
    pub fn colour(v: f64, scale: Option<(f64, f64)>) -> Option<[u8; 3]> {
        let (lo, hi) = scale?;
        if !v.is_finite() || v <= 0.0 {
            return None;
        }
        let u = (v.log10() - lo.log10()) / (hi.log10() - lo.log10());
        Some(pantometry::view::ramp::sequential(u.clamp(0.0, 1.0)))
    }

    /// The diagram of frame `frame` of `run` in the rectangle `[left, top, right, bottom]`, in
    /// pixels.
    ///
    /// The discs come in the order of [`Body::organs`] and then of [`Body::depots`], which is what
    /// a test matches them by.
    pub fn draw(&self, run: &Run, frame: usize, rect: [f64; 4]) -> Drawing {
        let mut d = Drawing::default();
        let shown = run.frames.get(frame);
        let values =
            shown.map_or_else(|| vec![None; self.organs.len()], |f| self.concentrations(f));
        let scale = self.scale(run);

        // The layout's extent in its own units, arrows and depots included, fitted to the room
        // left once the labels and the colour bar have theirs — turned to whichever of `TURNS`
        // angles that room holds largest, so a tall viewport gets a tall body. A turn moves every
        // point together and changes no distance, so no area and no width means anything else.
        let room = [
            (rect[2] - rect[0] - 2.0 * PAD_X - BAR).max(1.0),
            (rect[3] - rect[1] - PAD_TOP - PAD_BOTTOM).max(1.0),
        ];
        let fit = |turn: f64| {
            let (sin, cos) = turn.sin_cos();
            let mut lo = [f64::INFINITY; 2];
            let mut hi = [f64::NEG_INFINITY; 2];
            let mut grow = |p: [f64; 2], r: f64| {
                let q = [cos * p[0] - sin * p[1], sin * p[0] + cos * p[1]];
                for a in 0..2 {
                    lo[a] = lo[a].min(q[a] - r);
                    hi[a] = hi[a].max(q[a] + r);
                }
            };
            for o in &self.organs {
                grow(o.at, o.radius);
                if o.clearance_l_per_h > 0.0 {
                    grow(tip(o), 0.0);
                }
            }
            for p in &self.depots {
                grow(p.at, DEPOT_R);
            }
            let span = [(hi[0] - lo[0]).max(1e-9), (hi[1] - lo[1]).max(1e-9)];
            let s = (room[0] / span[0]).min(room[1] / span[1]);
            (s, lo, span, [sin, cos])
        };
        let mut best = fit(0.0);
        for k in 1..TURNS {
            let tried = fit(2.0 * std::f64::consts::PI * k as f64 / TURNS as f64);
            // Strictly larger by more than rounding, so a tie keeps the earlier and unturned one.
            if tried.0 > best.0 * (1.0 + 1e-9) {
                best = tried;
            }
        }
        let (s, lo, span, [sin, cos]) = best;
        let origin = [
            rect[0] + PAD_X + 0.5 * (room[0] - s * span[0]) - s * lo[0],
            rect[1] + PAD_TOP + 0.5 * (room[1] - s * span[1]) - s * lo[1],
        ];
        let turned = |p: [f64; 2]| [cos * p[0] - sin * p[1], sin * p[0] + cos * p[1]];
        let px = |p: [f64; 2]| {
            let q = turned(p);
            [origin[0] + s * q[0], origin[1] + s * q[1]]
        };

        let widest = self
            .vessels
            .iter()
            .map(|v| v.flow_l_per_h)
            .fold(0.0_f64, f64::max);
        for v in &self.vessels {
            let width = if widest > 0.0 {
                1.5 + 6.0 * v.flow_l_per_h / widest
            } else {
                1.5
            };
            d.strokes.push(Stroke {
                points: vec![px(self.organs[v.from].at), px(self.organs[v.to].at)],
                colour: VESSEL,
                width,
            });
        }
        for p in &self.depots {
            let o = &self.organs[p.into];
            let dir = unit([o.at[0] - p.at[0], o.at[1] - p.at[1]]);
            let from = [p.at[0] + DEPOT_R * dir[0], p.at[1] + DEPOT_R * dir[1]];
            let to = [o.at[0] - o.radius * dir[0], o.at[1] - o.radius * dir[1]];
            arrow(&mut d, px(from), px(to), INK);
        }
        for o in self.organs.iter().filter(|o| o.clearance_l_per_h > 0.0) {
            let t = turned(o.clears_towards);
            let from = [
                o.at[0] + o.radius * o.clears_towards[0],
                o.at[1] + o.radius * o.clears_towards[1],
            ];
            let end = px(tip(o));
            arrow(&mut d, px(from), end, INK);
            d.texts.push(Text {
                at: [end[0] + 6.0 * t[0].signum(), end[1] + 6.0 * t[1]],
                text: format!("CL {} L/h", sig(o.clearance_l_per_h)),
                anchor: if t[0] >= 0.0 {
                    Anchor::Left
                } else {
                    Anchor::Right
                },
                colour: INK,
            });
        }

        for (i, (o, v)) in self.organs.iter().zip(&values).enumerate() {
            let c = px(o.at);
            let r = s * o.radius;
            d.discs.push(Disc {
                centre: c,
                radius: r,
                fill: v.and_then(|v| Body::colour(v, scale)),
                outline: INK,
            });
            let value = match v {
                Some(v) => format!("{} {CONCENTRATION}", sig(*v)),
                None => "no reading".to_string(),
            };
            let away: Vec<[f64; 2]> = self.spokes(i).into_iter().map(turned).collect();
            label(&mut d, c, r, &away, &o.label, Some(&value));
        }
        // What is in the depot is reported for all of them together, so it can be put against a
        // depot only when there is one.
        let depot_mass = shown.and_then(|f| {
            f.readings
                .iter()
                .find(|r| r.domain == self.domain && r.label == "in the depot")
                .map(|r| r.value)
        });
        for p in &self.depots {
            let c = px(p.at);
            let r = s * DEPOT_R;
            d.discs.push(Disc {
                centre: c,
                radius: r,
                fill: None,
                outline: INK,
            });
            let amount = match (depot_mass, self.depots.len()) {
                (Some(m), 1) => Some(format!("{} mg", sig(m))),
                _ => None,
            };
            let o = &self.organs[p.into];
            let away = [turned([o.at[0] - p.at[0], o.at[1] - p.at[1]])];
            label(&mut d, c, r, &away, &p.label, amount.as_deref());
        }

        // The domain's name, and nothing more: its books in milligrams are the chart's lower plot,
        // and a line of five of them ran off a viewport the width of the editor's.
        d.texts.push(Text {
            at: [rect[0] + 8.0, rect[1] + 6.0],
            text: self.domain.clone(),
            anchor: Anchor::TopLeft,
            colour: INK,
        });

        if let Some((lo, hi)) = scale {
            colour_bar(&mut d, rect, (lo, hi));
        }
        d
    }
}

impl Body {
    /// Every direction something leaves compartment `i` by, in the layout's units: a vessel to
    /// each neighbour, a depot's arrow, the clearance arrow. Its label goes where these are not.
    fn spokes(&self, i: usize) -> Vec<[f64; 2]> {
        let o = &self.organs[i];
        let towards = |p: [f64; 2]| [p[0] - o.at[0], p[1] - o.at[1]];
        let mut out: Vec<[f64; 2]> = self
            .vessels
            .iter()
            .filter_map(|v| match (v.from == i, v.to == i) {
                (true, false) => Some(towards(self.organs[v.to].at)),
                (false, true) => Some(towards(self.organs[v.from].at)),
                _ => None,
            })
            .collect();
        out.extend(
            self.depots
                .iter()
                .filter(|p| p.into == i)
                .map(|p| towards(p.at)),
        );
        if o.clearance_l_per_h > 0.0 {
            out.push(o.clears_towards);
        }
        out
    }
}

/// A disc's name and value as a block of one or two lines just outside its rim, on whichever of
/// the eight compass points is furthest in angle from every direction in `away` — straight below
/// when nothing is in the way, since that is the first one tried.
fn label(d: &mut Drawing, c: [f64; 2], r: f64, away: &[[f64; 2]], name: &str, value: Option<&str>) {
    const LINE: f64 = 14.0;
    const H: f64 = std::f64::consts::FRAC_1_SQRT_2;
    // Screen compass points, y down: S, N, E, W, SE, SW, NE, NW.
    const POINTS: [[f64; 2]; 8] = [
        [0.0, 1.0],
        [0.0, -1.0],
        [1.0, 0.0],
        [-1.0, 0.0],
        [H, H],
        [-H, H],
        [H, -H],
        [-H, -H],
    ];
    let away: Vec<[f64; 2]> = away.iter().map(|a| unit(*a)).collect();
    let mut best = (f64::NEG_INFINITY, 0usize);
    for (k, p) in POINTS.iter().enumerate() {
        // The largest cosine to anything in the way; the point with the smallest wins.
        let near = away
            .iter()
            .map(|a| a[0] * p[0] + a[1] * p[1])
            .fold(-1.0_f64, f64::max);
        if -near > best.0 + 1e-9 {
            best = (-near, k);
        }
    }
    let p = POINTS[best.1];
    let at = [c[0] + (r + 3.0) * p[0], c[1] + (r + 3.0) * p[1]];
    let lines: Vec<(&str, [u8; 3])> = std::iter::once((name, NAME))
        .chain(value.map(|v| (v, INK)))
        .collect();
    let n = lines.len() as f64;
    for (j, (text, colour)) in lines.into_iter().enumerate() {
        let j = j as f64;
        let (dy, anchor) = match best.1 {
            0 => (j * LINE, Anchor::Top),
            1 => ((j - n + 1.0) * LINE, Anchor::Bottom),
            2 => ((j - 0.5 * (n - 1.0)) * LINE, Anchor::Left),
            3 => ((j - 0.5 * (n - 1.0)) * LINE, Anchor::Right),
            4 => (j * LINE, Anchor::TopLeft),
            5 => (j * LINE, Anchor::TopRight),
            6 => ((j - n + 1.0) * LINE, Anchor::BottomLeft),
            _ => ((j - n + 1.0) * LINE, Anchor::BottomRight),
        };
        d.texts.push(Text {
            at: [at[0], at[1] + dy],
            text: text.to_string(),
            anchor,
            colour,
        });
    }
}

/// Where a compartment's clearance arrow ends, in the layout's units.
fn tip(o: &Organ) -> [f64; 2] {
    let t = o.clears_towards;
    let reach = o.radius + ARROW;
    [o.at[0] + reach * t[0], o.at[1] + reach * t[1]]
}

fn unit(v: [f64; 2]) -> [f64; 2] {
    let n = (v[0] * v[0] + v[1] * v[1]).sqrt();
    if n > 0.0 {
        [v[0] / n, v[1] / n]
    } else {
        [1.0, 0.0]
    }
}

/// A line from `from` to `to` with a head at `to`, in pixels.
fn arrow(d: &mut Drawing, from: [f64; 2], to: [f64; 2], colour: [u8; 3]) {
    d.strokes.push(Stroke {
        points: vec![from, to],
        colour,
        width: 1.5,
    });
    let u = unit([to[0] - from[0], to[1] - from[1]]);
    let (len, half) = (8.0, 4.0);
    let back = [to[0] - len * u[0], to[1] - len * u[1]];
    d.strokes.push(Stroke {
        points: vec![
            [back[0] - half * u[1], back[1] + half * u[0]],
            to,
            [back[0] + half * u[1], back[1] - half * u[0]],
        ],
        colour,
        width: 1.5,
    });
}

/// The scale, as a column of short bands down the right of `rect` with each decade labelled.
fn colour_bar(d: &mut Drawing, rect: [f64; 4], (lo, hi): (f64, f64)) {
    const BANDS: usize = 48;
    let x = rect[2] - BAR + 12.0;
    let (top, bottom) = (rect[1] + PAD_TOP, rect[3] - PAD_BOTTOM);
    if bottom <= top {
        return;
    }
    let each = (bottom - top) / BANDS as f64;
    for k in 0..BANDS {
        let u = (k as f64 + 0.5) / BANDS as f64;
        let y1 = bottom - k as f64 * each;
        d.strokes.push(Stroke {
            points: vec![[x, y1], [x, y1 - each]],
            colour: pantometry::view::ramp::sequential(u),
            width: 12.0,
        });
    }
    let decades = (hi.log10() - lo.log10()).round().max(1.0);
    let mut v = lo;
    let mut k = 0.0;
    while k <= decades {
        let y = bottom - (k / decades) * (bottom - top);
        d.strokes.push(Stroke {
            points: vec![[x + 6.0, y], [x + 10.0, y]],
            colour: GRID,
            width: 1.0,
        });
        d.texts.push(Text {
            at: [x + 12.0, y],
            text: sig(v),
            anchor: Anchor::Left,
            colour: INK,
        });
        v *= 10.0;
        k += 1.0;
    }
    d.texts.push(Text {
        at: [x, top - 6.0],
        text: CONCENTRATION.to_string(),
        anchor: Anchor::Bottom,
        colour: INK,
    });
}

/// Three significant figures, which is what a label beside a disc has room for.
fn sig(v: f64) -> String {
    if !v.is_finite() {
        return "-".to_string();
    }
    if v == 0.0 {
        return "0".to_string();
    }
    let a = v.abs();
    if !(1e-3..1e5).contains(&a) {
        return format!("{v:.2e}");
    }
    let mut decimals = 0usize;
    let mut x = a;
    while x < 99.95 && decimals < 6 {
        x *= 10.0;
        decimals += 1;
    }
    format!("{v:.decimals$}")
}

/// What the viewport draws for a run with nothing to place: each compartment model as a body
/// above, when the scene has one, and the readings against time below — or the chart alone, in a
/// `width × height` canvas at frame `frame`. `None` when there is neither.
pub fn view(
    bodies: &[Body],
    run: &Run,
    log: bool,
    width: f64,
    height: f64,
    frame: usize,
) -> Option<Drawing> {
    let cursor = run.frames.get(frame).map(|f| f.t);
    let readings = chart::chart(run, log);
    // A body is drawn only when the run reports for it, so a scene edited since the run does not
    // draw a diagram of compartments the run never had.
    let drawn: Vec<&Body> = bodies
        .iter()
        .filter(|b| {
            run.frames
                .iter()
                .any(|f| f.readings.iter().any(|r| r.domain == b.domain))
        })
        .collect();
    if drawn.is_empty() {
        return readings.map(|c| c.draw(width, height, cursor));
    }
    let split = if readings.is_some() {
        (0.45 * height).round()
    } else {
        height
    };
    let mut d = Drawing::default();
    let each = width / drawn.len() as f64;
    for (i, b) in drawn.iter().enumerate() {
        let left = i as f64 * each;
        d.append(b.draw(run, frame, [left, 0.0, left + each, split]));
    }
    if let Some(c) = readings {
        d.strokes.push(Stroke {
            points: vec![[0.0, split], [width, split]],
            colour: GRID,
            width: 1.0,
        });
        d.append(c.draw(width, height - split, cursor).shifted(split));
    }
    Some(d)
}
