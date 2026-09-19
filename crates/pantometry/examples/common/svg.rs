//! Enough SVG to show a physics result, and no more.
//!
//! There is no plotting dependency here on purpose. This workspace has twelve external
//! crates, three of which reach a built artifact, and a chart library would be the largest
//! thing in the tree by an order of magnitude — for output that is a few hundred lines of
//! text. SVG *is* text: a `format!` and a file write, no encoder, no font handling, and it
//! opens by double-click on every platform.
//!
//! What it does: axes with ticks, polylines, filled rectangles for a raster, labels, and a
//! camera that puts points in space onto the page.
//! What it does not: legends laid out automatically, log axes with minor ticks, anything
//! interactive. When an example needs one of those, add it here rather than reaching for a
//! crate — and when *that* stops being reasonable, the answer is `pantometry-world`, not a
//! bigger version of this file.

#![allow(dead_code)]

use std::fmt::Write as _;

/// A drawing, in user coordinates that map linearly onto the canvas.
///
/// The mapping is set once by [`Plot::new`] and every draw call goes through it, so an
/// example never converts a temperature into a pixel by hand — which is the one place a
/// plot silently lies.
pub struct Plot {
    body: String,
    width: f64,
    height: f64,
    /// Data-space bounds: x from `x0` to `x1`, y from `y0` to `y1`.
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    /// The drawing area in canvas pixels: left, top, width, height. Set from sensible
    /// padding by [`Plot::new`] and overridable by [`Plot::viewport`], which is how several
    /// plots share one canvas.
    area: (f64, f64, f64, f64),
}

/// A camera for drawing points in space on a flat page.
///
/// Rotate about the vertical axis, tilt, then divide by depth: the same three steps
/// `pantometry-view`'s HTML viewer takes, so a still and the page a reader can spin agree about
/// what the scene looks like. No lighting and no hidden-surface removal, because what this is for
/// is a **layout** -- the lines are the subject and a solid would hide them.
pub struct View3 {
    azimuth: f64,
    elevation: f64,
    distance: f64,
    centre: [f64; 3],
    span: f64,
    /// The box the projected scene occupies on the page, measured at construction.
    page: (f64, f64, f64, f64),
}

impl View3 {
    /// A camera framed on everything it will be asked to draw.
    ///
    /// `azimuth` and `elevation` are radians; `distance` is in units of the scene's own size, so
    /// the framing does not move when the scene grows.
    ///
    /// The points are walked twice -- once for the box they occupy in space, once for the box
    /// they occupy on the page after projection. The second is not the first: perspective is not
    /// a scaling, and a bounding box projected corner by corner is not the projection's bounding
    /// box.
    pub fn framing(
        points: impl IntoIterator<Item = [f64; 3]>,
        azimuth: f64,
        elevation: f64,
        distance: f64,
    ) -> View3 {
        let pts: Vec<[f64; 3]> = points.into_iter().collect();
        let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
        for p in &pts {
            for a in 0..3 {
                lo[a] = lo[a].min(p[a]);
                hi[a] = hi[a].max(p[a]);
            }
        }
        let centre = if pts.is_empty() {
            [0.0; 3]
        } else {
            [0, 1, 2].map(|a| (lo[a] + hi[a]) / 2.0)
        };
        let span = if pts.is_empty() {
            1.0
        } else {
            (0..3).fold(0.0f64, |s, a| s.max(hi[a] - lo[a])).max(1e-30)
        };
        let mut view = View3 {
            azimuth,
            elevation,
            distance,
            centre,
            span,
            page: (-0.5, 0.5, -0.5, 0.5),
        };
        let (mut x0, mut x1, mut y0, mut y1) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        for p in &pts {
            let (x, y, _) = view.project(*p);
            x0 = x0.min(x);
            x1 = x1.max(x);
            y0 = y0.min(y);
            y1 = y1.max(y);
        }
        if x0 <= x1 {
            view.page = (x0, x1, y0, y1);
        }
        view
    }

    /// A point on the page, and its depth. Nearer is smaller.
    pub fn project(&self, p: [f64; 3]) -> (f64, f64, f64) {
        let n = [0, 1, 2].map(|a| (p[a] - self.centre[a]) / self.span);
        let (ca, sa) = (self.azimuth.cos(), self.azimuth.sin());
        let (x1, z1) = (n[0] * ca - n[2] * sa, n[0] * sa + n[2] * ca);
        let (ce, se) = (self.elevation.cos(), self.elevation.sin());
        let (y1, z2) = (n[1] * ce - z1 * se, n[1] * se + z1 * ce);
        // A point level with the eye has no place on the page at all. Clamping is what the
        // report's viewer does, and it keeps a line that reaches the eye plane from flying off
        // to infinity and taking the framing with it.
        let d = (z2 + self.distance).max(0.05);
        (x1 / d, y1 / d, d)
    }

    /// The projected box, **grown** until it is `aspect` wide for its height.
    ///
    /// [`Plot`] maps x and y onto the viewport independently, so a projection handed a viewport
    /// of a different shape comes out stretched -- a lens drawn as an ellipse, and nothing in the
    /// picture to say it happened. Growing rather than cropping means the whole scene is still
    /// there; the cost is margin, which is cheap.
    pub fn page_bounds(&self, aspect: f64) -> ((f64, f64), (f64, f64)) {
        let (x0, x1, y0, y1) = self.page;
        let (mut w, mut h) = ((x1 - x0).max(1e-30), (y1 - y0).max(1e-30));
        if w / h < aspect {
            w = h * aspect;
        } else {
            h = w / aspect;
        }
        let (cx, cy) = ((x0 + x1) / 2.0, (y0 + y1) / 2.0);
        ((cx - w / 2.0, cx + w / 2.0), (cy - h / 2.0, cy + h / 2.0))
    }
}

/// A run resolved for drawing: how far away it is, where its points land, and how to ink it.
///
/// Named rather than left a tuple because it is four fields deep, and clippy is right that
/// nobody can read `Vec<(f64, Vec<(f64, f64)>, String, f64)>` and say which `f64` is which.
struct Drawn {
    depth: f64,
    points: Vec<(f64, f64)>,
    colour: String,
    width: f64,
}

/// A colour, written the way SVG wants it.
pub fn rgb(r: u8, g: u8, b: u8) -> String {
    format!("#{r:02x}{g:02x}{b:02x}")
}

/// Black-body-ish ramp from cold to hot, for temperature and intensity rasters.
///
/// Perceptually rough — a real visualiser would use something with uniform lightness — but
/// it has the property that matters for a heat map: it is monotonic in brightness, so a
/// greyscale print still reads the right way round.
pub fn heat(v: f64) -> String {
    // Quantised to 48 steps before colouring. More than the eye resolves in a small panel,
    // and it is what lets `Plot::raster` merge runs -- an unquantised ramp changes colour
    // every single cell and run-length encoding then saves nothing at all.
    let v = (v.clamp(0.0, 1.0) * 48.0).round() / 48.0;
    let r = (255.0 * (v * 2.2).clamp(0.0, 1.0)).round() as u8;
    let g = (255.0 * ((v - 0.28) * 1.7).clamp(0.0, 1.0)).round() as u8;
    let b = (255.0 * ((v - 0.68) * 3.0).clamp(0.0, 1.0)).round() as u8;
    rgb(r, g, b)
}

/// A blue–white–red ramp for a signed field, with `v = 0` mapping to the middle.
///
/// [`heat`] is wrong for a pressure and the reason is not aesthetic. A temperature is
/// positive and climbs, so a ramp anchored at the minimum tells the truth about it. A
/// pressure swings either side of zero and averages out, so the same ramp puts the midpoint
/// wherever the extremes happen to fall this frame — a room at rest would be drawn as though
/// half of it were cold, and the colours would shift as the wave passed rather than staying
/// pinned to the physics.
///
/// So zero is fixed at white and `v` is in −1 to 1 against a symmetric scale the caller
/// chooses, usually the peak magnitude. Quantised for the same reason [`heat`] is.
pub fn diverging(v: f64) -> String {
    let v = (v.clamp(-1.0, 1.0) * 32.0).round() / 32.0;
    let m = v.abs();
    // Toward red for positive, toward blue for negative, white in between.
    let (hot, cold) = ((178.0, 34.0, 34.0), (28.0, 76.0, 168.0));
    let (r, g, b) = if v >= 0.0 { hot } else { cold };
    let mix = |c: f64| (255.0 + (c - 255.0) * m).round() as u8;
    rgb(mix(r), mix(g), mix(b))
}

impl Plot {
    pub fn new(width: f64, height: f64, x: (f64, f64), y: (f64, f64)) -> Plot {
        Plot {
            body: String::new(),
            width,
            height,
            x0: x.0,
            x1: x.1,
            y0: y.0,
            y1: y.1,
            area: (64.0, 34.0, width - 82.0, height - 78.0),
        }
    }

    /// Put this plot's drawing area somewhere specific on the canvas, so several can share
    /// one document. Combine the results with [`document`].
    pub fn viewport(mut self, x: f64, y: f64, width: f64, height: f64) -> Plot {
        self.area = (x, y, width, height);
        self
    }

    /// Data x to canvas x.
    fn px(&self, x: f64) -> f64 {
        let span = if self.x1 == self.x0 {
            1.0
        } else {
            self.x1 - self.x0
        };
        self.area.0 + (x - self.x0) / span * self.area.2
    }

    /// Data y to canvas y, flipped — SVG counts downwards and physics does not.
    fn py(&self, y: f64) -> f64 {
        let span = if self.y1 == self.y0 {
            1.0
        } else {
            self.y1 - self.y0
        };
        self.area.1 + self.area.3 - (y - self.y0) / span * self.area.3
    }

    pub fn polyline(&mut self, points: impl IntoIterator<Item = (f64, f64)>, colour: &str, w: f64) {
        let pts: Vec<String> = points
            .into_iter()
            .map(|(x, y)| format!("{:.2},{:.2}", self.px(x), self.py(y)))
            .collect();
        if pts.is_empty() {
            return;
        }
        let _ = write!(
            self.body,
            r#"<polyline fill="none" stroke="{colour}" stroke-width="{w}" stroke-linejoin="round" points="{}"/>"#,
            pts.join(" ")
        );
    }

    /// Runs of points in space, drawn back to front through a [`View3`].
    ///
    /// Sorted by each run's mean depth, which is the only hidden-line rule a set of polylines can
    /// have: a line in front covers one behind, rather than whichever happened to be last in the
    /// array. It is what `pantometry-view`'s viewer does with a `Paths` panel, for the same
    /// reason and with the same limitation -- two lines that cross are ordered as wholes.
    ///
    /// A run of fewer than two points is drawn as nothing by [`Plot::polyline`] and is **not**
    /// dropped here, so a caller counting what it handed over gets the count back.
    pub fn paths3(
        &mut self,
        view: &View3,
        runs: impl IntoIterator<Item = (Vec<[f64; 3]>, String, f64)>,
    ) {
        let mut ordered: Vec<Drawn> = runs
            .into_iter()
            .map(|(pts, colour, width)| {
                let projected: Vec<(f64, f64, f64)> =
                    pts.iter().map(|&p| view.project(p)).collect();
                // An empty run has no depth to sort by, and furthest is the harmless answer: it
                // is drawn first and draws nothing.
                let depth = if projected.is_empty() {
                    f64::MAX
                } else {
                    projected.iter().map(|p| p.2).sum::<f64>() / projected.len() as f64
                };
                Drawn {
                    depth,
                    points: projected.into_iter().map(|(x, y, _)| (x, y)).collect(),
                    colour,
                    width,
                }
            })
            .collect();
        // Furthest first. `total_cmp` rather than `partial_cmp`, because a NaN depth would make
        // the comparison inconsistent and the sort's order arbitrary.
        ordered.sort_by(|a, b| b.depth.total_cmp(&a.depth));
        for run in ordered {
            self.polyline(run.points, &run.colour, run.width);
        }
    }

    /// A filled cell of a raster, given in data coordinates.
    ///
    /// Half a pixel of overlap is added on purpose: adjacent rectangles that share an exact
    /// edge show hairline seams in most renderers, and a heat map full of white lines looks
    /// like structure that is not in the data.
    pub fn cell(&mut self, x: (f64, f64), y: (f64, f64), colour: &str) {
        let (a, b) = (self.px(x.0), self.px(x.1));
        let (c, d) = (self.py(y.1), self.py(y.0));
        let _ = write!(
            self.body,
            r#"<rect x="{:.2}" y="{:.2}" width="{:.2}" height="{:.2}" fill="{colour}" shape-rendering="crispEdges"/>"#,
            a.min(b) - 0.25,
            c.min(d) - 0.25,
            (b - a).abs() + 0.5,
            (d - c).abs() + 0.5
        );
    }

    /// A whole raster over a rectangular data region, run-length encoded by row.
    ///
    /// One `<rect>` per cell is the obvious implementation and it produced a 1.4 MB file for
    /// a 120×120 Airy pattern. Merging horizontally adjacent cells of the same colour is the
    /// cheap half of the fix.
    ///
    /// Only the cheap half, and the measurements are worth writing down because the first
    /// version of this was disappointing. Run-length encoding alone saved **17%** — a smooth
    /// radial gradient changes colour at every single cell, so there are no runs to merge.
    /// Quantising the ramp first (see [`heat`]) took the saving to **32%**, and the rest of
    /// the way was fewer cells: 96 across rather than 120. 1.4 MB became 620 KB, most of it
    /// from the resolution rather than from the clever part.
    ///
    /// If a raster ever needs to be genuinely small, the answer is an `<image>` with an
    /// encoded bitmap, not a better encoding of rectangles.
    ///
    /// `shade` is called once per cell with its column and row.
    pub fn raster<F>(
        &mut self,
        cols: usize,
        rows: usize,
        x: (f64, f64),
        y: (f64, f64),
        mut shade: F,
    ) where
        F: FnMut(usize, usize) -> String,
    {
        let edge = |lo: f64, hi: f64, k: usize, n: usize| lo + (hi - lo) * k as f64 / n as f64;
        for row in 0..rows {
            let (y0, y1) = (edge(y.0, y.1, row, rows), edge(y.0, y.1, row + 1, rows));
            let mut start = 0usize;
            let mut current = shade(0, row);
            for col in 1..=cols {
                let next = if col < cols {
                    Some(shade(col, row))
                } else {
                    None
                };
                if next.as_deref() != Some(current.as_str()) {
                    self.cell(
                        (edge(x.0, x.1, start, cols), edge(x.0, x.1, col, cols)),
                        (y0, y1),
                        &current,
                    );
                    start = col;
                    if let Some(n) = next {
                        current = n;
                    }
                }
            }
        }
    }

    pub fn text(&mut self, x: f64, y: f64, s: &str, size: f64, colour: &str, anchor: &str) {
        let _ = write!(
            self.body,
            r#"<text x="{:.2}" y="{:.2}" font-family="ui-sans-serif,system-ui,sans-serif" font-size="{size}" fill="{colour}" text-anchor="{anchor}">{}</text>"#,
            self.px(x),
            self.py(y),
            escape(s)
        );
    }

    /// A label placed in canvas coordinates, for titles and captions that should not move
    /// when the data range does.
    pub fn label(&mut self, cx: f64, cy: f64, s: &str, size: f64, colour: &str, anchor: &str) {
        let _ = write!(
            self.body,
            r#"<text x="{cx:.2}" y="{cy:.2}" font-family="ui-sans-serif,system-ui,sans-serif" font-size="{size}" fill="{colour}" text-anchor="{anchor}">{}</text>"#,
            escape(s)
        );
    }

    /// Axes with ticks and a light grid. `fmt` turns a tick value into its label, so the
    /// caller decides on units and decimals rather than this file guessing.
    pub fn axes(
        &mut self,
        xticks: &[f64],
        yticks: &[f64],
        fmt_x: impl Fn(f64) -> String,
        fmt_y: impl Fn(f64) -> String,
    ) {
        let (grid, axis, ink) = ("#00000018", "#7a7a7a", "#3a3a3a");
        for &x in xticks {
            let (a, b) = (self.px(x), self.py(self.y0));
            let _ = write!(
                self.body,
                r#"<line x1="{a:.2}" y1="{b:.2}" x2="{a:.2}" y2="{:.2}" stroke="{grid}" stroke-width="1"/>"#,
                self.py(self.y1)
            );
            let _ = write!(
                self.body,
                r#"<text x="{a:.2}" y="{:.2}" font-family="ui-sans-serif,system-ui,sans-serif" font-size="12" fill="{ink}" text-anchor="middle">{}</text>"#,
                b + 18.0,
                escape(&fmt_x(x))
            );
        }
        for &y in yticks {
            let (a, b) = (self.px(self.x0), self.py(y));
            let _ = write!(
                self.body,
                r#"<line x1="{a:.2}" y1="{b:.2}" x2="{:.2}" y2="{b:.2}" stroke="{grid}" stroke-width="1"/>"#,
                self.px(self.x1)
            );
            let _ = write!(
                self.body,
                r#"<text x="{:.2}" y="{:.2}" font-family="ui-sans-serif,system-ui,sans-serif" font-size="12" fill="{ink}" text-anchor="end">{}</text>"#,
                a - 8.0,
                b + 4.0,
                escape(&fmt_y(y))
            );
        }
        // The two axis lines last, so they sit over the grid.
        let _ = write!(
            self.body,
            r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="{axis}" stroke-width="1.2"/><line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="{axis}" stroke-width="1.2"/>"#,
            self.px(self.x0),
            self.py(self.y0),
            self.px(self.x1),
            self.py(self.y0),
            self.px(self.x0),
            self.py(self.y0),
            self.px(self.x0),
            self.py(self.y1)
        );
    }

    pub fn title(&mut self, s: &str) {
        self.label(self.area.0, 20.0, s, 15.0, "#1b1b1b", "start");
    }

    pub fn caption(&mut self, s: &str) {
        self.label(self.width - 18.0, 20.0, s, 12.0, "#6a6a6a", "end");
    }

    /// A caption under this plot's own area, for a panel in a multi-plot figure.
    pub fn footnote(&mut self, s: &str) {
        let (x, y, w, h) = self.area;
        self.label(x + w / 2.0, y + h + 34.0, s, 11.5, "#6a6a6a", "middle");
    }

    /// The drawing without the document around it, for combining with [`document`].
    pub fn into_body(self) -> String {
        self.body
    }

    /// The finished document.
    ///
    /// Both light and dark viewers get a readable page: the background is painted, so a
    /// dark-mode browser does not leave black text on black.
    pub fn finish(self) -> String {
        let (w, h) = (self.width, self.height);
        document(w, h, [self.into_body()])
    }
}

/// Wrap one or more plot bodies in a document.
///
/// Doubled hashes on the literal: a `"#` inside would close an `r#"..."#` early, and the
/// background colour has one.
pub fn document(width: f64, height: f64, parts: impl IntoIterator<Item = String>) -> String {
    let body: String = parts.into_iter().collect();
    format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}"><rect width="{width}" height="{height}" fill="#fbfbfa"/>{body}</svg>"##
    )
}

fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// Round numbers covering a range, for axis ticks.
///
/// Steps of 1, 2 or 5 times a power of ten, which is what makes an axis readable — 0, 25,
/// 50, 75 rather than 0, 23.7, 47.4.
pub fn ticks(lo: f64, hi: f64, target: usize) -> Vec<f64> {
    if hi.is_nan() || lo.is_nan() || hi <= lo || target == 0 {
        return vec![lo];
    }
    let raw = (hi - lo) / target as f64;
    let magnitude = 10f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 5.0, 10.0]
        .iter()
        .map(|m| m * magnitude)
        .find(|s| *s >= raw)
        .unwrap_or(10.0 * magnitude);
    let first = (lo / step).ceil() * step;
    let mut out = Vec::new();
    let mut v = first;
    while v <= hi * (1.0 + 1e-12) {
        // Snap values that are a rounding error away from zero, so an axis shows "0" and
        // not "-0.0000000000000001".
        out.push(if v.abs() < step * 1e-9 { 0.0 } else { v });
        v += step;
    }
    out
}
