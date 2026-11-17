//! **Results, read off a run's fields.**
//!
//! A run is frames, and a field in a frame is samples on a box. What a person asks of it once the
//! picture is on screen is a number: the value at a point, a profile along a line, a cut through
//! a plane, the hottest place in a part and its mean, and how each of those went over time. Every
//! commercial package has these under a results node, and until this module the editor had a
//! readout under the pointer and nothing else.
//!
//! It is here rather than in either shell for the reason [`crate::chart`] is: two shells reading
//! one run would otherwise be two chances to put a number somewhere else.
//!
//! # What it promises
//!
//! - **Between samples, trilinear.** Exact for a field linear in x, y and z, and at a sample the
//!   sample itself, bit for bit: the weight of every other neighbour is zero and is not
//!   multiplied in.
//! - **A hole is a hole.** A cell the domain does not occupy is `NaN` in the run. A point whose
//!   neighbours include one has no value, rather than an average of the others: a probe on the
//!   edge of a bracket would otherwise report a temperature made up of the half of its neighbours
//!   that exist, which is a number no cell has.
//! - **Between the last sample and the wall, the last sample.** A cell-centred field has no
//!   sample on the box's faces. Inside the box and past the outermost sample, the value is that
//!   sample's — the cell's own average — and not a line extended past the data, which can
//!   overshoot what any cell holds.
//! - **An integral is the lattice's own rule.** Cell averages are summed by the midpoint rule and
//!   samples on the faces by the trapezoid rule. Each is exact for a linear field and second
//!   order otherwise, and each gives a sample the share of the box its spacing says it stands
//!   for. Applied the other way round they differ by a quarter of a cell at every face: on scene
//!   15's block, trapezoid weights on cell averages put the mean 0.035 K above the 293.2323 K the
//!   initial condition fixes.
//! - **In the panel's own frame.** Coordinates are the domain's own, in metres — the frame its
//!   grid is stated in. A placed domain is drawn where its pose puts it, through
//!   [`viewer_core::Placed::apply`], and asked about in its own coordinates, so a probe stays on
//!   the same material when the part is moved.
//! - **A frame without the domain is a gap, not a zero.** A history has one point per frame, and
//!   a frame that holds no such field gives `NaN`, which a chart draws as a break.

use viewer_core::{Frame, Lattice, Panel, Placed, Run};

/// One field of one frame, ready to be asked about.
#[derive(Clone, Copy, Debug)]
pub struct Grid<'a> {
    /// The domain it came from.
    pub name: &'a str,
    /// What its values are in.
    pub unit: &'a str,
    /// Samples along x, y and z.
    pub counts: [usize; 3],
    /// The box it was sampled over, `[x0, y0, z0, x1, y1, z1]` in metres, in its own frame.
    pub extent: [f64; 6],
    /// Where the samples sit inside the box.
    pub lattice: Lattice,
    /// `counts` multiplied together, x fastest; `NaN` where the domain is not.
    pub values: &'a [f64],
    /// Where its own frame sits in the world.
    pub place: Placed,
}

/// The statistics of a field over a region of it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Stats {
    /// The smallest value.
    pub min: f64,
    /// Where it is: the first sample, in grid order, that holds it.
    pub min_at: [f64; 3],
    /// The largest value.
    pub max: f64,
    /// Where it is, the same way.
    pub max_at: [f64; 3],
    /// The integral divided by the measure: the mean **by volume**, not by sample.
    pub mean: f64,
    /// The value integrated over the region, in the unit times metres to the power of `dims`.
    pub integral: f64,
    /// The region's length, area or volume in metres to the power of `dims`: what the samples
    /// counted stand for.
    pub measure: f64,
    /// How many of the box's axes have a thickness: 3 for a solid, 2 for a sheet, 1 for a line.
    pub dims: usize,
    /// How many samples were counted.
    pub samples: usize,
    /// How many in the region were holes and were not.
    pub holes: usize,
}

/// A field cut by a plane at right angles to one axis.
#[derive(Clone, Debug, PartialEq)]
pub struct Slice {
    /// The axis the plane is at right angles to: 0 for x, 1 for y, 2 for z.
    pub axis: usize,
    /// Where along it, in metres.
    pub at: f64,
    /// The two axes that lie in the plane, in x-y-z order: the image's across and up.
    pub axes: [usize; 2],
    /// Samples across and up.
    pub counts: [usize; 2],
    /// The plane's extent, `[u0, v0, u1, v1]` in metres.
    pub extent: [f64; 4],
    /// Where the samples sit inside it.
    pub lattice: Lattice,
    /// `counts[0] * counts[1]` values, across fastest; `NaN` where the domain is not.
    pub values: Vec<f64>,
}

/// The name of an axis, for a label.
pub const AXIS_NAMES: [&str; 3] = ["x", "y", "z"];

impl<'a> Grid<'a> {
    /// The field in `panel`, or why there is none to ask about.
    ///
    /// Every refusal names the domain, so a shell can put it where the result would have been.
    pub fn of(panel: &'a Panel) -> Result<Grid<'a>, String> {
        let Panel::Field {
            name,
            unit,
            place,
            nx,
            ny,
            nz,
            extent_m,
            lattice,
            values,
        } = panel
        else {
            return Err(format!(
                "{} is bodies, paths or a surface, not a field — a slice, a line and a probe read a field",
                panel.name()
            ));
        };
        let Some(extent) = *extent_m else {
            return Err(format!(
                "{name}: this run does not say where its samples are — it was written before the \
                 format carried a field's extent, and a position in cell units is not one in metres"
            ));
        };
        let counts = [*nx, *ny, *nz];
        if counts.contains(&0) || values.len() != nx * ny * nz {
            return Err(format!(
                "{name}: {} values for a {nx} x {ny} x {nz} grid",
                values.len()
            ));
        }
        Ok(Grid {
            name,
            unit,
            counts,
            extent,
            lattice: *lattice,
            values,
            place: *place,
        })
    }

    /// Where sample `i` sits along `axis`, in metres.
    pub fn sample(&self, axis: usize, i: usize) -> f64 {
        let (lo, hi) = (self.extent[axis], self.extent[axis + 3]);
        lo + self.lattice.along(i, self.counts[axis]) * (hi - lo)
    }

    /// The sample at grid index `[i, j, k]`.
    pub fn value(&self, at: [usize; 3]) -> f64 {
        self.values[at[0] + self.counts[0] * (at[1] + self.counts[1] * at[2])]
    }

    /// The middle of the box.
    pub fn centre(&self) -> [f64; 3] {
        std::array::from_fn(|a| 0.5 * (self.extent[a] + self.extent[a + 3]))
    }

    /// The distance between neighbouring samples along `axis`, or infinity along an axis with
    /// one sample, which no line crosses a cell of.
    pub fn spacing(&self, axis: usize) -> f64 {
        let n = self.counts[axis];
        let len = self.extent[axis + 3] - self.extent[axis];
        match (self.lattice, n) {
            (_, 1) => f64::INFINITY,
            (Lattice::Nodal, _) => len / (n - 1) as f64,
            (Lattice::Centred, _) => len / n as f64,
        }
    }

    /// How many of the box's axes have a thickness.
    pub fn dims(&self) -> usize {
        (0..3)
            .filter(|&a| self.extent[a + 3] > self.extent[a])
            .count()
    }

    /// The two samples either side of `x` along `axis` and the weight of the second, or `None`
    /// when `x` is outside the box.
    ///
    /// An axis with no thickness takes any `x`: a sheet has one place across it and every point
    /// asked about is there. An axis with one sample and a thickness holds `x` to the box, like
    /// any other.
    fn bracket(&self, axis: usize, x: f64) -> Option<(usize, usize, f64)> {
        let n = self.counts[axis];
        let (lo, hi) = (self.extent[axis], self.extent[axis + 3]);
        if hi == lo {
            return Some((0, 0, 0.0));
        }
        if !(x >= lo && x <= hi) {
            return None;
        }
        if n == 1 {
            return Some((0, 0, 0.0));
        }
        if x <= self.sample(axis, 0) {
            return Some((0, 0, 0.0));
        }
        if x >= self.sample(axis, n - 1) {
            return Some((n - 1, n - 1, 0.0));
        }
        // **The pair from the lattice, and the weight from the positions.** The index can come
        // out one low of the sample `x` sits on, and then the weight is `(b - a) / (b - a)`,
        // which is exactly one; or `x` can be a rounding short of a sample and land on it, where
        // the weight clamps to zero. Either way a point at a sample reads that sample to the bit,
        // and a point between two reads them in proportion. A correction loop against the
        // positions was here and the bit-for-bit test passed without it: `numerics-reviewer`
        // found no input it changed.
        let f = (x - lo) / (hi - lo);
        let u = match self.lattice {
            Lattice::Nodal => f * (n - 1) as f64,
            Lattice::Centred => f * n as f64 - 0.5,
        };
        let i = (u.floor().max(0.0) as usize).min(n - 2);
        let (a, b) = (self.sample(axis, i), self.sample(axis, i + 1));
        Some((i, i + 1, ((x - a) / (b - a)).clamp(0.0, 1.0)))
    }

    /// The value at `p`, in the panel's own frame, or `None` outside the box or beside a hole.
    pub fn at(&self, p: [f64; 3]) -> Option<f64> {
        let b = [
            self.bracket(0, p[0])?,
            self.bracket(1, p[1])?,
            self.bracket(2, p[2])?,
        ];
        let pairs = |(lo, hi, w): (usize, usize, f64)| [(lo, 1.0 - w), (hi, w)];
        let mut sum = 0.0;
        for (i, wi) in pairs(b[0]) {
            for (j, wj) in pairs(b[1]) {
                for (k, wk) in pairs(b[2]) {
                    let w = wi * wj * wk;
                    // Not multiplied in: zero times a hole is a hole, and a weight of zero is a
                    // neighbour that does not count.
                    if w == 0.0 {
                        continue;
                    }
                    let v = self.value([i, j, k]);
                    if !v.is_finite() {
                        return None;
                    }
                    sum += w * v;
                }
            }
        }
        Some(sum)
    }

    /// How many points a line from `a` to `b` is read at: four in every cell it crosses, along
    /// whichever axis it crosses most of, and never fewer than two or more than 4001.
    ///
    /// Trilinear values are straight between samples and bend at them, so a line read once a
    /// cell would cut every bend off; four times shows each.
    pub fn points_along(&self, a: [f64; 3], b: [f64; 3]) -> usize {
        let cells = (0..3)
            .map(|axis| (b[axis] - a[axis]).abs() / self.spacing(axis))
            .fold(0.0, f64::max);
        if !cells.is_finite() {
            return 2;
        }
        ((4.0 * cells).ceil() as usize + 1).clamp(2, 4001)
    }

    /// The field along the straight line from `a` to `b`, at `n` points evenly spaced and at
    /// least two: `(distance from a in metres, value)`, `NaN` where the line is outside the box
    /// or beside a hole, which a chart draws as a break.
    pub fn line(&self, a: [f64; 3], b: [f64; 3], n: usize) -> Vec<[f64; 2]> {
        let n = n.max(2);
        let length = distance(a, b);
        (0..n)
            .map(|k| {
                let f = k as f64 / (n - 1) as f64;
                // `a (1 - f) + b f` and not `a + (b - a) f`: the last point is then `b` itself.
                // The difference misses it by an ulp on eleven pairs of the test grid's samples,
                // which moves no value read off them — a ulp of position is 1e-16 of a cell — so
                // this is tidiness, and a sabotage of it passes every test here.
                let p = std::array::from_fn(|i| a[i] * (1.0 - f) + b[i] * f);
                [length * f, self.at(p).unwrap_or(f64::NAN)]
            })
            .collect()
    }

    /// The field cut by the plane at right angles to `axis` at `at` metres, or `None` when the
    /// plane is outside the box.
    ///
    /// Read at the samples of the other two axes and interpolated along `axis` only, so a plane
    /// through a row of samples is that row, bit for bit.
    pub fn slice(&self, axis: usize, at: f64) -> Option<Slice> {
        let (k0, k1, w) = self.bracket(axis, at)?;
        let axes = match axis {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        let counts = [self.counts[axes[0]], self.counts[axes[1]]];
        let mut values = Vec::with_capacity(counts[0] * counts[1]);
        for j in 0..counts[1] {
            for i in 0..counts[0] {
                let mut idx = [0usize; 3];
                idx[axes[0]] = i;
                idx[axes[1]] = j;
                let mut sum = 0.0;
                let mut hole = false;
                for (k, wk) in [(k0, 1.0 - w), (k1, w)] {
                    if wk == 0.0 {
                        continue;
                    }
                    idx[axis] = k;
                    let v = self.value(idx);
                    if !v.is_finite() {
                        hole = true;
                        break;
                    }
                    sum += wk * v;
                }
                values.push(if hole { f64::NAN } else { sum });
            }
        }
        Some(Slice {
            axis,
            at,
            axes,
            counts,
            extent: [
                self.extent[axes[0]],
                self.extent[axes[1]],
                self.extent[axes[0] + 3],
                self.extent[axes[1] + 3],
            ],
            lattice: self.lattice,
            values,
        })
    }

    /// What sample `i` along `axis` stands for: its share of the box's length along that axis.
    ///
    /// One along an axis with no thickness, which is then not part of the measure at all: a
    /// sheet's integral is over its area.
    pub fn weight(&self, axis: usize, i: usize) -> f64 {
        let n = self.counts[axis];
        let len = self.extent[axis + 3] - self.extent[axis];
        if len == 0.0 {
            return 1.0;
        }
        match (self.lattice, n) {
            (_, 1) => len,
            (Lattice::Centred, _) => len / n as f64,
            (Lattice::Nodal, _) => {
                let h = len / (n - 1) as f64;
                if i == 0 || i == n - 1 {
                    0.5 * h
                } else {
                    h
                }
            }
        }
    }

    /// The statistics of the field over `within` — a box `[x0, y0, z0, x1, y1, z1]` in the
    /// panel's frame, either corner first — or over the whole domain for `None`. `None` when the
    /// region holds no sample that is not a hole: there is no mean of nothing, and zero is a
    /// temperature.
    ///
    /// A box counts the samples inside it, each with the share of the domain it has in the whole
    /// grid: the region a box selects is the cells whose sample it holds.
    pub fn region(&self, within: Option<[f64; 6]>) -> Option<Stats> {
        let inside = |axis: usize, x: f64| match within {
            None => true,
            Some(b) => {
                if self.extent[axis + 3] == self.extent[axis] {
                    return true;
                }
                let (lo, hi) = (b[axis].min(b[axis + 3]), b[axis].max(b[axis + 3]));
                x >= lo && x <= hi
            }
        };
        let mut stats = Stats {
            min: f64::INFINITY,
            min_at: [f64::NAN; 3],
            max: f64::NEG_INFINITY,
            max_at: [f64::NAN; 3],
            mean: f64::NAN,
            integral: 0.0,
            measure: 0.0,
            dims: self.dims(),
            samples: 0,
            holes: 0,
        };
        let [nx, ny, nz] = self.counts;
        for k in 0..nz {
            let z = self.sample(2, k);
            if !inside(2, z) {
                continue;
            }
            for j in 0..ny {
                let y = self.sample(1, j);
                if !inside(1, y) {
                    continue;
                }
                let wyz = self.weight(1, j) * self.weight(2, k);
                for i in 0..nx {
                    let x = self.sample(0, i);
                    if !inside(0, x) {
                        continue;
                    }
                    let v = self.value([i, j, k]);
                    if !v.is_finite() {
                        stats.holes += 1;
                        continue;
                    }
                    let w = self.weight(0, i) * wyz;
                    stats.samples += 1;
                    stats.integral += w * v;
                    stats.measure += w;
                    if v < stats.min {
                        stats.min = v;
                        stats.min_at = [x, y, z];
                    }
                    if v > stats.max {
                        stats.max = v;
                        stats.max_at = [x, y, z];
                    }
                }
            }
        }
        if stats.samples == 0 {
            return None;
        }
        stats.mean = stats.integral / stats.measure;
        Some(stats)
    }
}

impl Slice {
    /// The value at sample `i` across and `j` up.
    pub fn value(&self, i: usize, j: usize) -> f64 {
        self.values[i + self.counts[0] * j]
    }

    /// The rectangle an image of this slice covers, `[u0, v0, u1, v1]` in metres, with one pixel
    /// per sample: each pixel's **centre** on its sample.
    ///
    /// The plane's own extent for cell averages, whose samples are cell centres already. Half a
    /// spacing wider on every side for samples on the faces, whose outermost samples sit on the
    /// edges — drawn into the plane's own extent, a nodal field's picture would put every sample
    /// up to half a pixel out, the error `Lattice` was introduced to stop.
    pub fn image_extent(&self) -> [f64; 4] {
        let mut e = self.extent;
        if self.lattice == Lattice::Nodal {
            for a in 0..2 {
                let n = self.counts[a];
                if n > 1 {
                    let half = 0.5 * (e[a + 2] - e[a]) / (n - 1) as f64;
                    e[a] -= half;
                    e[a + 2] += half;
                }
            }
        }
        e
    }

    /// The slice as an image, `counts[0]` wide and `counts[1]` high, four bytes a pixel, the
    /// first row at the top — the largest `v` — and a hole transparent.
    ///
    /// **In the view's colours, dimming included.** A glowing field is drawn darker where it is
    /// cooler — [`crate::Colouring::brightness`] — and the view applies that to linear light. So
    /// does this, before encoding back to sRGB; a ramp's brightness is one and passes through.
    pub fn rgba(&self, colouring: &crate::Colouring) -> Vec<u8> {
        let [w, h] = self.counts;
        let mut out = Vec::with_capacity(w * h * 4);
        for row in 0..h {
            let j = h - 1 - row;
            for i in 0..w {
                let v = self.value(i, j);
                if v.is_finite() {
                    let k = colouring.brightness(v);
                    let [r, g, b] = colouring.srgb(v).map(|c| dimmed(c, k));
                    out.extend_from_slice(&[r, g, b, 255]);
                } else {
                    out.extend_from_slice(&[0, 0, 0, 0]);
                }
            }
        }
        out
    }

    /// The point in the panel's own frame at `u`, `v` metres in the plane.
    pub fn point(&self, u: f64, v: f64) -> [f64; 3] {
        let mut p = [0.0; 3];
        p[self.axes[0]] = u;
        p[self.axes[1]] = v;
        p[self.axis] = self.at;
        p
    }
}

/// An sRGB channel at `k` of its brightness: decoded to linear light, scaled, encoded again. One
/// is the channel itself, which is every ramp.
fn dimmed(c: u8, k: f64) -> u8 {
    if k >= 1.0 {
        return c;
    }
    let s = f64::from(c) / 255.0;
    let linear = if s <= 0.04045 {
        s / 12.92
    } else {
        ((s + 0.055) / 1.055).powf(2.4)
    } * k.max(0.0);
    let encoded = if linear <= 0.003_130_8 {
        12.92 * linear
    } else {
        1.055 * linear.powf(1.0 / 2.4) - 0.055
    };
    (encoded * 255.0).round().clamp(0.0, 255.0) as u8
}

/// The straight-line distance from `a` to `b`.
pub fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d: [f64; 3] = std::array::from_fn(|i| b[i] - a[i]);
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

/// The field called `name` in `frame`, or why there is none.
pub fn grid<'a>(frame: &'a Frame, name: &str) -> Result<Grid<'a>, String> {
    let panel = frame
        .panels
        .iter()
        .find(|p| p.name() == name)
        .ok_or_else(|| format!("{name} is not in this frame"))?;
    Grid::of(panel)
}

/// Every domain the run has a field for, in the order the run first names them: the domains a
/// slice, a line, a probe or a region can be asked of.
///
/// **Every field, including one [`Grid::of`] refuses.** A field whose values do not fill its grid
/// or whose run does not say where its samples are was left out, so the refusal written for it
/// could never be shown and the panel said "nothing in this run is a field", which was false.
/// Listed, the tab says why it cannot read it.
pub fn fields(run: &Run) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for frame in &run.frames {
        for panel in &frame.panels {
            if matches!(panel, Panel::Field { .. }) && !out.iter().any(|n| n == panel.name()) {
                out.push(panel.name().to_string());
            }
        }
    }
    out
}

/// The value at `p` in the field called `name`, in every frame: `(t in seconds, value)`, with
/// `NaN` for a frame that has no such field or no value there.
pub fn history(run: &Run, name: &str, p: [f64; 3]) -> Vec<[f64; 2]> {
    run.frames
        .iter()
        .map(|frame| {
            let v = grid(frame, name)
                .ok()
                .and_then(|g| g.at(p))
                .unwrap_or(f64::NAN);
            [frame.t, v]
        })
        .collect()
}

/// [`Grid::region`] in every frame: `(t in seconds, the statistics or None)`.
pub fn region_history(
    run: &Run,
    name: &str,
    within: Option<[f64; 6]>,
) -> Vec<(f64, Option<Stats>)> {
    run.frames
        .iter()
        .map(|frame| {
            (
                frame.t,
                grid(frame, name).ok().and_then(|g| g.region(within)),
            )
        })
        .collect()
}

/// The unit an integral over `dims` axes is in: the value's, times metres that many times.
pub fn integral_unit(unit: &str, dims: usize) -> String {
    match dims {
        0 => unit.to_string(),
        1 => format!("{unit} m"),
        n => format!("{unit} m^{n}"),
    }
}

/// What a measure over `dims` axes is called and the unit it is in.
pub fn measure_name(dims: usize) -> (&'static str, &'static str) {
    match dims {
        0 => ("samples", ""),
        1 => ("length", "m"),
        2 => ("area", "m^2"),
        _ => ("volume", "m^3"),
    }
}
