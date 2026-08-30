//! Reads a pantometry run and turns it into something a renderer can draw.
//!
//! No GPU, no window, no `pantometry` dependency. That last one is deliberate and is the point of the
//! crate: if a viewer can be written against the **file** a run produces and nothing else, then
//! the wire format is complete. If it needed to reach back into the library for something the
//! file did not carry, the format would be the thing to fix — and finding that out is worth more
//! than the convenience of linking.
//!
//! # What a run is
//!
//! What `pantometry_view::to_json` writes: a title and a list of frames, each with a time, some
//! panels and some readings. A panel is one of three shapes — a field on a grid, a set of points,
//! or runs of connected points — and the reader below accepts all three by name.
//!
//! # What this crate is responsible for
//!
//! Everything a renderer would otherwise get wrong the same way twice:
//!
//! - **One colour scale across the whole run.** A frame normalised to itself makes a decay look
//!   like a steady state, which is the one thing a picture of a simulation must never do. The
//!   scale is computed once, over every frame, per panel.
//! - **One framing across the whole run**, from the geometry rather than from the current frame,
//!   so a body moving is a body moving and not the camera chasing it.
//! - **The projection**, so the window and any test agree about where a point lands.

#![deny(missing_docs)]
#![forbid(unsafe_code)]

use serde::Deserialize;

/// The highest run format this build understands.
///
/// Written by `pantometry_view::data::FORMAT` and read here, and the two are held together by
/// `the_writer_and_the_reader_agree_on_the_version` rather than by anyone remembering. They are
/// in different workspaces on purpose -- this crate does not link the library -- so a constant
/// each is the only way, and a test comparing them is the only thing that makes it one number.
pub const FORMAT: u32 = 1;

/// A whole run, as read from a file.
#[derive(Clone, Debug, Deserialize)]
pub struct Run {
    /// What the run was called.
    pub title: String,
    /// Every captured instant, in order.
    pub frames: Vec<Frame>,
}

/// One instant.
#[derive(Clone, Debug, Deserialize)]
pub struct Frame {
    /// Simulation time, in seconds.
    pub t: f64,
    /// One per drawable domain.
    pub panels: Vec<Panel>,
    /// Named scalars, drawable or not.
    #[serde(default)]
    pub readings: Vec<Reading>,
}

/// One named scalar.
#[derive(Clone, Debug, Deserialize)]
pub struct Reading {
    /// Which domain it came from.
    pub domain: String,
    /// What it is.
    pub label: String,
    /// Its value.
    pub value: f64,
    /// Its unit.
    pub unit: String,
}

/// One domain, captured, in whichever shape it has.
///
/// Tagged by `kind`, which is what the writer emits. An unknown kind is an **error** rather than a
/// panel that is silently skipped: a viewer that quietly drops a shape it does not know shows a
/// page with something missing and no way to tell.
#[derive(Clone, Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "lowercase", deny_unknown_fields)]
pub enum Panel {
    /// A field sampled on a grid of one, two or three dimensions.
    Field {
        /// Which domain this came from.
        name: String,
        /// What the values are in.
        unit: String,
        /// Samples along x.
        nx: usize,
        /// Along y.
        ny: usize,
        /// Along z.
        nz: usize,
        /// The box the field was sampled over, `[x0, y0, z0, x1, y1, z1]` in metres — or `None`
        /// for a run written before the format carried it.
        ///
        /// **Optional, and that is the compatibility story.** This enum is
        /// `deny_unknown_fields`, deliberately: a wire format that silently discards a key it
        /// does not know is how a renamed field once went unnoticed here. The price is that
        /// adding a key breaks every existing reader, and it did — the day `pantometry-view`
        /// began writing `extent_m`, this crate refused every run file the library produced,
        /// with the library's own twenty-step gate unable to see it because the viewer is a
        /// separate workspace with a separate job. `Option` with a `default` is what makes the
        /// direction work both ways: a new reader takes an old file, and an old reader is the
        /// thing that cannot be helped.
        #[serde(default)]
        extent_m: Option<[f64; 6]>,
        /// `nx*ny*nz` values, x fastest then y then z.
        values: Vec<f64>,
    },
    /// A countable set of bodies.
    Points {
        /// Which domain this came from.
        name: String,
        /// What the values are in.
        unit: String,
        /// Whether the bounding box is a real wall.
        boxed: bool,
        /// `[x0,y0,z0,x1,y1,z1]`.
        bounds: [f64; 6],
        /// Flattened `xyz` per body.
        positions: Vec<f64>,
        /// One per body.
        values: Vec<f64>,
    },
    /// Runs of connected points — rays, trajectories, field lines.
    Paths {
        /// Which domain this came from.
        name: String,
        /// What the values are in.
        unit: String,
        /// `[x0,y0,z0,x1,y1,z1]`.
        bounds: [f64; 6],
        /// Where each run begins, as an index into `vertices` divided by three.
        starts: Vec<f64>,
        /// Flattened `xyz` per vertex.
        vertices: Vec<f64>,
        /// One per run.
        values: Vec<f64>,
    },
}

impl Panel {
    /// Which domain this panel came from.
    pub fn name(&self) -> &str {
        match self {
            Panel::Field { name, .. } | Panel::Points { name, .. } | Panel::Paths { name, .. } => {
                name
            }
        }
    }

    /// What its values are in.
    pub fn unit(&self) -> &str {
        match self {
            Panel::Field { unit, .. } | Panel::Points { unit, .. } | Panel::Paths { unit, .. } => {
                unit
            }
        }
    }

    /// The values, whichever shape this is.
    pub fn values(&self) -> &[f64] {
        match self {
            Panel::Field { values, .. }
            | Panel::Points { values, .. }
            | Panel::Paths { values, .. } => values,
        }
    }

    /// The box this panel occupies, in world coordinates and in metres.
    ///
    /// A field says so directly now. It did not: the run file recorded a grid and not the extent
    /// it was sampled over, so this returned **the grid in cell units** — a 9x9x9 block framed as
    /// a nine-metre cube, at the origin however far from it the part really was, and a slab that
    /// was asked for at more slices drawn taller for it. That was stated in a doc comment rather
    /// than fixed because the number was genuinely not in the file. It is now.
    ///
    /// The cell-unit box is still the answer for a file written before the format carried the
    /// extent, and [`Panel::extent_m`] is how a caller tells the two apart rather than being
    /// handed a plausible number with no provenance.
    pub fn bounds(&self) -> [f64; 6] {
        match self {
            Panel::Field {
                nx,
                ny,
                nz,
                extent_m,
                ..
            } => extent_m.unwrap_or([0.0, 0.0, 0.0, *nx as f64, *ny as f64, *nz as f64]),
            Panel::Points { bounds, .. } | Panel::Paths { bounds, .. } => *bounds,
        }
    }

    /// The box a field was sampled over, in metres, or `None` if this is not a field or the run
    /// predates the format carrying it.
    ///
    /// The distinction matters to anything that would put a number on screen: [`Panel::bounds`]
    /// falls back to cell units, and a scale bar reading "9 m" under a 40 mm part is worse than
    /// no scale bar.
    pub fn extent_m(&self) -> Option<[f64; 6]> {
        match self {
            Panel::Field { extent_m, .. } => *extent_m,
            _ => None,
        }
    }
}

impl Run {
    /// Read a run from the JSON a `pantometry` run wrote.
    ///
    /// # The version is read before the panels, and that ordering is the point
    ///
    /// A run states a `format`. This reads **that key alone** first, and refuses a number it
    /// does not understand before trying to parse anything else.
    ///
    /// The order matters because of what fails otherwise. `PanelData` is `deny_unknown_fields`
    /// and tagged by `kind`, so a run carrying a shape this build has never heard of fails
    /// inside serde with a message about an unknown variant — which is what a *corrupt* file
    /// says too, and the two want different words in front of a person. The scene format checks
    /// its version after parsing and has exactly this gap; this does not.
    ///
    /// An absent `format` is **1**, because every run written before the key existed is one.
    pub fn from_json(text: &str) -> Result<Run, String> {
        /// Just the version, from a reader that ignores everything else in the file.
        #[derive(Deserialize)]
        struct Stamp {
            #[serde(default = "one")]
            format: u32,
        }
        fn one() -> u32 {
            1
        }

        // A file that is not JSON at all fails here, and says so as "not a pantometry run"
        // rather than as a complaint about a version.
        let stamp: Stamp =
            serde_json::from_str(text).map_err(|e| format!("not a pantometry run: {e}"))?;
        if stamp.format == 0 || stamp.format > FORMAT {
            return Err(format!(
                "this run is format {}, and this build reads {FORMAT}. Upgrade pantometry to open it",
                stamp.format
            ));
        }
        serde_json::from_str(text).map_err(|e| format!("not a pantometry run: {e}"))
    }

    /// The value range of one panel **across every frame**.
    ///
    /// The whole run, not the current frame. A per-frame scale makes a decaying mode look
    /// constant and a constant one look like noise, and it is what a renderer does if nobody
    /// stops it. Returns `None` for a panel name the run does not have.
    pub fn scale_of(&self, panel: &str) -> Option<(f64, f64)> {
        let mut lo = f64::MAX;
        let mut hi = f64::MIN;
        let mut seen = false;
        for frame in &self.frames {
            for p in &frame.panels {
                if p.name() != panel {
                    continue;
                }
                seen = true;
                for v in p.values() {
                    lo = lo.min(*v);
                    hi = hi.max(*v);
                }
            }
        }
        seen.then_some((lo, hi))
    }

    /// The box one panel occupies **across every frame**, widened to hold all of them.
    ///
    /// Same reason as the scale: a camera framed to the current frame follows a moving body and
    /// makes it look still.
    pub fn framing_of(&self, panel: &str) -> Option<[f64; 6]> {
        let mut out = [f64::MAX, f64::MAX, f64::MAX, f64::MIN, f64::MIN, f64::MIN];
        let mut seen = false;
        for frame in &self.frames {
            for p in &frame.panels {
                if p.name() != panel {
                    continue;
                }
                seen = true;
                let b = p.bounds();
                for a in 0..3 {
                    out[a] = out[a].min(b[a]);
                    out[a + 3] = out[a + 3].max(b[a + 3]);
                }
            }
        }
        seen.then_some(out)
    }

    /// Every panel name in the run, in the order they first appear.
    pub fn panels(&self) -> Vec<String> {
        let mut names: Vec<String> = Vec::new();
        for frame in &self.frames {
            for p in &frame.panels {
                if !names.iter().any(|n| n == p.name()) {
                    names.push(p.name().to_string());
                }
            }
        }
        names
    }
}

/// Where the eye is, in orbit around what it is looking at.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Camera {
    /// Azimuth, radians.
    pub azimuth: f64,
    /// Elevation, radians, clamped to just under a right angle so the up vector never degenerates.
    pub elevation: f64,
    /// Distance, in units of the framed box's longest side. This is how strong the perspective
    /// is; it is **not** the knob for how big the subject looks — see [`Camera::fit`].
    pub distance: f64,
    /// Focal length, as a multiple of the default. One is the field of view the HTML report uses,
    /// so the two open on the same picture; [`Camera::fit`] moves it to frame the subject.
    pub scale: f64,
}

impl Default for Camera {
    /// Three-quarter view from slightly above — the same angles the HTML report opens at, so a
    /// reader moving between the two is not re-orienting.
    fn default() -> Camera {
        Camera {
            azimuth: 0.7,
            elevation: 0.4,
            distance: 2.5,
            scale: 1.0,
        }
    }
}

impl Camera {
    /// Turn by a drag, in radians. Elevation is clamped; azimuth wraps.
    pub fn turn(&mut self, d_az: f64, d_el: f64) {
        self.azimuth += d_az;
        self.elevation = (self.elevation + d_el).clamp(-1.5, 1.5);
    }

    /// Zoom by a scroll. Bounded, so a viewer cannot be scrolled inside the subject or into
    /// the next county and left with a blank window and no way back.
    pub fn zoom(&mut self, factor: f64) {
        self.distance = (self.distance * factor).clamp(1.2, 9.0);
    }

    /// Choose the focal length that makes the run's bounding box fill the frame.
    ///
    /// # Why this is the focal length and not the distance
    ///
    /// [`Framing`] normalises by the box's **longest** side, so a run that is tall and thin fills
    /// the unit along one axis and a fraction of it along the other two. A fixed camera is
    /// therefore a camera set up for a cube: a portafilter, 67 mm across and 118 mm tall, came out
    /// at 15% of the frame height and 0.29% of its pixels.
    ///
    /// Backing *in* to fix that does not work. At this field of view a subject one unit across
    /// fills the frame at a distance of about 0.35, which is inside its own bounding box — the
    /// eye ends up between the near and far faces, [`Camera::project`]'s depth clamp starts firing
    /// and the near geometry stretches across the window. Distance is the wrong knob; it sets how
    /// strong the perspective is, and here that was never the problem.
    ///
    /// So this sets the focal length, and the projection is **linear** in it. One pass gives the
    /// exact answer, with no iteration and no convergence test to get wrong.
    ///
    /// `fill` is where the furthest corner should land, with 1.0 the edge of the shorter side.
    pub fn fit(&mut self, bounds: [f64; 6], frame: &Framing, aspect: f64, fill: f64) {
        self.scale = 1.0;
        let mut worst: f64 = 0.0;
        for i in 0..8 {
            let corner = [
                if i & 1 == 0 { bounds[0] } else { bounds[3] },
                if i & 2 == 0 { bounds[1] } else { bounds[4] },
                if i & 4 == 0 { bounds[2] } else { bounds[5] },
            ];
            let q = self.project(corner, frame, aspect);
            worst = worst.max(q.x.abs()).max(q.y.abs());
        }
        if worst > 0.0 {
            self.scale = fill.clamp(0.05, 0.98) / worst;
        }
    }

    /// The rotation this camera applies, as three rows.
    ///
    /// **Written once** because two things need it: [`Camera::project`], which a 2D painter uses a
    /// point at a time, and [`Camera::matrix`], which a GPU uses for every vertex at once. A camera
    /// with two ideas of which way is up puts the labels somewhere the geometry is not.
    fn rotation(&self) -> [[f64; 3]; 3] {
        let (ca, sa) = (self.azimuth.cos(), self.azimuth.sin());
        let (ce, se) = (self.elevation.cos(), self.elevation.sin());
        [
            [ca, 0.0, -sa],
            [-sa * se, ce, -ca * se],
            [sa * ce, se, ca * ce],
        ]
    }

    /// The near plane. The same number [`Camera::project`] clamps its divisor to, so the flat path
    /// and the GPU path lose a point at the same place.
    pub const NEAR: f64 = 0.05;

    /// The far plane. The framed box is a unit cube, so its half-diagonal is 0.87 and two units
    /// past the pivot clears it with room for a run that reached outside the box it was framed on.
    pub fn far(&self) -> f64 {
        self.distance + 2.0
    }

    /// The projection as a 4×4, column-major, ready for `glUniformMatrix4fv`.
    ///
    /// Clip space, so `w` is exactly the depth [`Camera::project`] reports and `x/w`, `y/w` are
    /// exactly its `x` and `y` — `the_matrix_and_the_projection_agree` holds that. What the matrix
    /// adds is a `z` for a depth buffer and honest clipping at the near plane, where `project` can
    /// only clamp because a 2D painter has nothing to clip against.
    ///
    /// # It takes framing-local coordinates, not world ones
    ///
    /// Feed it [`Framing::local`]. There is no centre or span in here, and that is deliberate: a
    /// matrix that folded them in would carry `centre / span` as an `f32`, and a pantometry scene is
    /// routinely a 9 mm block sitting 200 mm from the origin. Measured, that arrangement disagreed
    /// with [`Camera::project`] by `4.4e-6` — two numbers near 22 subtracting to 0.1, which is
    /// `f32` keeping two digits of seven. `runtime/gpu`'s README has the same finding about storing
    /// absolute kelvin, and it is the same fix: subtract in `f64`, hand the GPU the small number.
    pub fn matrix(&self, aspect: f64) -> [f32; 16] {
        let r = self.rotation();
        let fy = 0.6 * self.scale;
        let fx = fy / aspect.max(1e-6);

        // Depth: `z/w` runs -1 at the near plane to 1 at the far one, the OpenGL convention, and
        // `w` is the eye distance so the mapping is the usual `1/z` one.
        let (n, f) = (Camera::NEAR, self.far());
        let a = (f + n) / (f - n);
        let b = -2.0 * f * n / (f - n);

        let rows = [
            [fx * r[0][0], fx * r[0][1], fx * r[0][2], 0.0],
            [fy * r[1][0], fy * r[1][1], fy * r[1][2], 0.0],
            [a * r[2][0], a * r[2][1], a * r[2][2], a * self.distance + b],
            [r[2][0], r[2][1], r[2][2], self.distance],
        ];

        // Column-major: `m[4 * column + row]`.
        let mut m = [0.0f32; 16];
        for (i, r) in rows.iter().enumerate() {
            for (j, v) in r.iter().enumerate() {
                m[4 * j + i] = *v as f32;
            }
        }
        m
    }

    /// Project a world point into normalised device coordinates, with its depth.
    ///
    /// `x` and `y` run `-1..1` across the shorter side; `depth` grows away from the eye and is
    /// what a renderer sorts or tests by. `aspect` is width over height.
    pub fn project(&self, p: [f64; 3], frame: &Framing, aspect: f64) -> Projected {
        let d = frame.span.max(1e-12);
        let v = [
            (p[0] - frame.centre[0]) / d,
            (p[1] - frame.centre[1]) / d,
            (p[2] - frame.centre[2]) / d,
        ];
        let r = self.rotation();
        let at = |row: [f64; 3]| row[0] * v[0] + row[1] * v[1] + row[2] * v[2];
        let (x1, y1, z2) = (at(r[0]), at(r[1]), at(r[2]));

        // Never divide by a depth at or behind the eye: a point there has no screen position, and
        // returning a huge coordinate instead of clamping is how a renderer draws a streak across
        // the window and calls it geometry.
        let depth = (z2 + self.distance).max(Camera::NEAR);
        let f = 0.6 * self.scale / depth;
        Projected {
            x: x1 * f / aspect.max(1e-6),
            y: y1 * f,
            depth,
        }
    }
}

/// A projected point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Projected {
    /// Normalised device x.
    pub x: f64,
    /// Normalised device y.
    pub y: f64,
    /// Distance from the eye. Larger is further.
    pub depth: f64,
}

impl Framing {
    /// A world point in this framing's own units: centred on the subject, one unit across.
    ///
    /// The subtraction happens in `f64` and only the result narrows, which is the whole point —
    /// see [`Camera::matrix`]. Vertices are then order one however far from the origin the scene
    /// sits and however small it is, which is the range `f32` is good at.
    pub fn local(&self, p: [f64; 3]) -> [f32; 3] {
        let d = self.span.max(1e-12);
        [
            ((p[0] - self.centre[0]) / d) as f32,
            ((p[1] - self.centre[1]) / d) as f32,
            ((p[2] - self.centre[2]) / d) as f32,
        ]
    }
}

/// The box a panel is drawn in, reduced to a centre and one span.
///
/// **One span for all three axes**, so a cube is drawn as a cube. Scaling each axis to fill the
/// window independently is the easiest way to make a picture that is not of the thing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Framing {
    /// The centre of the box.
    pub centre: [f64; 3],
    /// The longest side.
    pub span: f64,
}

impl Framing {
    /// From a bounding box.
    pub fn of(bounds: [f64; 6]) -> Framing {
        let centre = [
            0.5 * (bounds[0] + bounds[3]),
            0.5 * (bounds[1] + bounds[4]),
            0.5 * (bounds[2] + bounds[5]),
        ];
        let span = (bounds[3] - bounds[0])
            .max(bounds[4] - bounds[1])
            .max(bounds[5] - bounds[2]);
        Framing {
            centre,
            span: if span > 0.0 { span } else { 1.0 },
        }
    }
}

/// One line segment, ready to hand to a renderer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Segment {
    /// Both ends, projected.
    pub from: Projected,
    /// The other end.
    pub to: Projected,
    /// Where this segment's value sits in the run's range, `0..1`.
    pub shade: f64,
}

/// Build the segments for one panel at one frame, shaded against `span`.
///
/// Sorted **back to front**, so a renderer without a depth buffer still draws them in an order
/// that reads. A renderer with one can ignore the order and lose nothing.
///
/// Returns an empty list for a field, which is not a set of lines — a field view is a different
/// pipeline and pretending otherwise here would produce a plausible picture of nothing.
///
/// # `span` is the run's, not the frame's
///
/// It was the frame's, computed here from the values in hand, and that is the mistake
/// [`Run::scale_of`] exists to prevent — a method this crate had, tested, and never called from
/// the renderer. A scale that re-fits every frame makes a quantity look **constant while it
/// changes by orders of magnitude**: water that leaves a shower screen clean and arrives at the
/// spout at 83 kg/m³ renders mid-ramp the whole way down, because at every instant it is halfway
/// between that instant's lightest and darkest.
///
/// Pass `run.scale_of(name)`. Passing a degenerate span is allowed and shades everything at the
/// bottom of the ramp, which is the honest picture of a quantity that does not vary.
pub fn segments(
    panel: &Panel,
    camera: &Camera,
    framing: &Framing,
    aspect: f64,
    span: (f64, f64),
) -> Vec<Segment> {
    let mut out = Vec::new();
    if let Panel::Paths {
        starts,
        vertices,
        values,
        ..
    } = panel
    {
        let (lo, hi) = span;
        let width = if hi > lo { hi - lo } else { 1.0 };
        for (k, start) in starts.iter().enumerate() {
            let from = *start as usize;
            let to = starts
                .get(k + 1)
                .map(|s| *s as usize)
                .unwrap_or(vertices.len() / 3);
            let shade = ((values.get(k).copied().unwrap_or(0.0) - lo) / width).clamp(0.0, 1.0);
            for i in from..to.saturating_sub(1) {
                let a = camera.project(vertex(vertices, i), framing, aspect);
                let b = camera.project(vertex(vertices, i + 1), framing, aspect);
                out.push(Segment {
                    from: a,
                    to: b,
                    shade,
                });
            }
        }
    }
    out.sort_by(|a, b| {
        let (da, db) = (a.from.depth + a.to.depth, b.from.depth + b.to.depth);
        db.total_cmp(&da)
    });
    out
}

fn vertex(flat: &[f64], i: usize) -> [f64; 3] {
    [flat[3 * i], flat[3 * i + 1], flat[3 * i + 2]]
}
