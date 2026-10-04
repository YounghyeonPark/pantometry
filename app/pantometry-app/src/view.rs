//! A native window for a pantometry run.
//!
//! ```text
//! cargo run --release -- run.json
//! ```
//!
//! Drag to rotate, scroll to zoom, space to play, left and right to scrub, F to frame each panel in
//! turn. `--frame-panel NAME` opens framed on one panel, in the window or a snapshot.
//!
//! ```text
//! cargo run --release -- run.json --snapshot out.ppm
//! ```
//!
//! renders one frame **without opening a window** and writes it out. That mode is not a
//! convenience: it is the only way anything can check that this file draws. A window nobody can
//! photograph proves it did not panic, which is a much weaker claim than it looks — every wrong
//! projection, every empty vertex buffer and every silently-failed pipeline also does not panic.
//!
//! # What is here and what is deliberately not
//!
//! This is the shell. Everything that could be got wrong twice — the colour scale across a run,
//! the framing, the projection — is in `viewer-core`, which has no GPU dependency and is tested
//! against real run files. What is left here is a surface, two pipelines and an event loop.
//!
//! That split is not tidiness. A renderer is the one place where a wrong answer looks like a
//! picture, so the arithmetic lives where a test can reach it and this file only draws what it is
//! handed.
//!
//! # The geometry is the exporters', and that is why this links `pantometry`
//!
//! **This said "it does not depend on `pantometry`, deliberately", and the claim has moved.** What
//! it bought was that a viewer written against the run file alone demonstrates the wire format
//! carries enough to draw a run — and that property lives in `viewer-core`, which still links
//! nothing, and is held by `one_run_format_two_crates` rather than by this file.
//!
//! What this file gained by giving it up is that a field is drawn as a **solid**. The triangles
//! come from `editor_core::field_shell`, which is `pantometry_view::mesh` — the same function that
//! writes the glTF and the USD and that the editor's viewport shades. The alternative was a second
//! implementation of where a field's boundary is, and that arithmetic has already been wrong once
//! here: a 40 mm cube exported 80 mm across. One picture, not three.

use std::sync::Arc;

use image::ImageEncoder;
use viewer_core::{segments, Camera, Framing, Run};
use winit::application::ApplicationHandler;
use winit::event::{ElementState, MouseScrollDelta, WindowEvent};
use winit::event_loop::{ActiveEventLoop, ControlFlow, EventLoop};
use winit::keyboard::{KeyCode, PhysicalKey};
use winit::window::{Window, WindowId};

/// One vertex, as the shader wants it: screen `x`, screen `y`, and depth for the buffer.
///
/// **Three components, because a solid needs a depth test.** The camera projects on the CPU — the
/// same `viewer_core::Camera` a test can reach — so what reaches the GPU is already in clip space
/// and the third component is `Projected::depth`, which the pass compares rather than sorts.
#[repr(C)]
#[derive(Clone, Copy)]
struct Vertex {
    position: [f32; 3],
    colour: [f32; 3],
}

/// A stand-in for `bytemuck`, which is one dependency for one trait.
mod bytemuck_lite {
    /// Marker for a type that can be reinterpreted as bytes.
    ///
    /// # Safety
    ///
    /// Implementors must be `#[repr(C)]`, contain no padding and no references. Both types here
    /// are arrays of `f32`, which satisfies all three.
    pub unsafe trait Pod: Copy {}

    /// The bytes of a slice.
    pub fn cast_slice<T: Pod>(v: &[T]) -> &[u8] {
        // SAFETY: `T: Pod` promises a plain-data layout, and the length is derived from the
        // slice's own, so the range is exactly the memory the slice owns.
        unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, std::mem::size_of_val(v)) }
    }
}

// SAFETY: `#[repr(C)]`, five `f32` with no padding and no references.
unsafe impl bytemuck_lite::Pod for Vertex {}

/// `pantometry view <run.json> [--snapshot out.ppm]`.
pub fn run(args: &[String]) -> i32 {
    let path = match args.first() {
        Some(p) => p.clone(),
        None => {
            eprintln!(
                "usage: pantometry view <run.json> [--frame-panel NAME] [--snapshot out.png [--frame N | --all-frames]]"
            );
            eprintln!(
                "  --frame-panel NAME frames one panel rather than the whole run; F cycles in the window"
            );
            eprintln!("  produced by `pantometry run <scene> out.json`, or by any run that calls");
            eprintln!("  pantometry_view::to_json");
            return 2;
        }
    };
    let text = match std::fs::read_to_string(&path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("cannot read {path}: {e}");
            std::process::exit(1);
        }
    };
    let run = match Run::from_json(&text) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    };

    // Which panel to show: the first one there is.
    //
    // **It used to be the first one with `Paths`**, and it refused everything else — "fields and
    // point clouds are different pipelines and are not built yet". They are not different
    // pipelines: a body and a field sample are points, and a point is two short segments in
    // screen space, which is the pipeline that was already here. What the old rule cost is that
    // `PanelData::paths` is built by two examples and by tests and by no `Domain` at all, so
    // **none of the shipped scenes could be drawn by this shell** — the refusal was the
    // only thing it ever said about a scene.
    let Some(panel) = run.panels().into_iter().next() else {
        eprintln!("{}: this run has no panels at all", run.title);
        std::process::exit(1);
    };

    println!(
        "{} — {} frames, showing {panel}",
        run.title,
        run.frames.len()
    );
    println!(
        "  drag to rotate, scroll to zoom, space to play, left/right to scrub, F to frame each panel in turn"
    );

    // **Which box the camera frames: the whole run, or one panel.** A run of several scales is
    // unreadable framed whole -- scene 14 puts a 20 mm bar beside orbits spanning 2.7e11 m, and the
    // bar is a point. Read here, before either mode, so the window and the snapshot frame alike.
    let framed: Option<String> = {
        let rest: Vec<&str> = args.iter().skip(1).map(String::as_str).collect();
        match rest.iter().position(|a| *a == "--frame-panel") {
            None => None,
            Some(i) => match rest.get(i + 1) {
                Some(name) if run.panels().iter().any(|p| p == name) => Some((*name).to_string()),
                said => {
                    eprintln!(
                        "--frame-panel takes a panel this run has ({}); got {:?}",
                        run.panels().join(", "),
                        said
                    );
                    return 2;
                }
            },
        }
    };

    // Headless: render one frame to a texture, read it back, write a PPM. No window, no display,
    // so a machine with neither can still check that the renderer puts something on the canvas.
    // From the dispatcher, past the run file. `std::env::args().skip(2)` was right when this was
    // its own binary and is off by one as a subcommand — it would have read the run file's path as
    // the flag.
    // **Found by scanning, like its siblings.** `--snapshot` had to be the argument straight after
    // the run file while `--frame` and `--thumbnail` are looked for anywhere, and the asymmetry
    // was silent both ways: `view run.json --thumbnail --snapshot out.png` fell past this branch
    // and opened a window, and `view run.json --snapshot --frame 5` took `--frame` as the output
    // path and wrote a 1100x720 PPM into a file called `--frame`.
    let rest: Vec<&str> = args.iter().skip(1).map(String::as_str).collect();
    if let Some(i) = rest.iter().position(|a| *a == "--snapshot") {
        let out = match rest.get(i + 1) {
            Some(p) if !p.starts_with("--") => (*p).to_string(),
            _ => {
                eprintln!(
                    "--snapshot takes a file to write; got {:?}",
                    rest.get(i + 1)
                );
                return 2;
            }
        };
        // Which frame, because the first one is often the least interesting: a run that fills up
        // over its length — a shot, a spreading spot — has nothing in it at `t = 0`, and a
        // snapshot of that is a picture of an empty box that still counts as "the renderer works".
        let at: usize = rest
            .iter()
            .position(|a| *a == "--frame")
            .and_then(|i| rest.get(i + 1))
            .and_then(|n| n.parse().ok())
            .unwrap_or(0)
            .min(run.frames.len().saturating_sub(1));
        // **A tile for the New-project chooser**: cropped to what was drawn and shrunk. The
        // cropping is the point — the camera fits the panel's bounds, and a one-dimensional bar
        // is a hairline in a frame that is otherwise background.
        let tile = rest.contains(&"--thumbnail");
        // **Every frame, from one load.** A forty-eight frame run written by forty-eight
        // invocations reparses the file and brings up a GPU each time: measured at 20 s a frame
        // against an 18.2 MB run, which is sixteen minutes. `docs/README.md` gives every figure a
        // command that refreshes it, and a figure whose command is a sixteen-minute shell loop is
        // a figure that ages quietly — which is the thing that page exists to argue against.
        let every = rest.contains(&"--all-frames");
        // **A `.gif` path collects instead of writing.** One command for a figure that moves:
        // the alternative is forty-eight numbered PNGs and an assembly step that lives somewhere
        // else, which is a figure refreshed by a recipe rather than by a command — and those are
        // the ones that go stale.
        let animate = every && out.to_ascii_lowercase().ends_with(".gif");
        let mut collected: Vec<Vec<u8>> = Vec::new();
        let mut size = (0u32, 0u32);
        let mut app = App::new(run, panel);
        app.frame_on(framed.clone());
        app.frame = at;
        app.legend = !tile;
        let count = app.run.frames.len();
        let wanted: Vec<usize> = if every {
            (0..count).collect()
        } else {
            vec![at]
        };
        if every {
            println!("  every one of {count} frames, from one load");
        } else {
            println!("  snapshot of frame {at} of {count}");
        }
        // Rendered three times the tile and averaged down, because a marker one pixel wide
        // vanishes under nearest-neighbour and a tile is not where to discover that.
        let (rw, rh) = if tile {
            (THUMB.0 * 3, THUMB.1 * 3)
        } else {
            (1100, 720)
        };
        for at in wanted {
            app.frame = at;
            let out = if every {
                numbered(&out, at)
            } else {
                out.clone()
            };
            match app.snapshot(rw, rh) {
                Ok(pixels) => {
                    let (pw, ph, pixels) = if tile {
                        let (cw, ch, cropped) = crop_to_content(rw, rh, &pixels);
                        (THUMB.0, THUMB.1, shrink(cw, ch, &cropped, THUMB.0, THUMB.1))
                    } else {
                        (rw, rh, pixels)
                    };
                    if animate {
                        collected.push(pixels.clone());
                    } else {
                        write_image(&out, pw, ph, &pixels);
                    }
                    // **Compared against the corner pixel, not against a constant.** The target is
                    // sRGB, so the clear colour is stored far brighter than the linear number the
                    // pass was given — 56,66,77 rather than 10,14,19. A fixed threshold called every
                    // pixel in the image a line and reported 100%, which is exactly the kind of
                    // "measurement" a renderer check exists to avoid.
                    let background = [pixels[0], pixels[1], pixels[2]];
                    // `as_chunks` rather than `chunks_exact`: each pixel is a `[u8; 4]` by type, which
                    // is what it is. This workspace has no MSRV to keep it off — see the note beside
                    // the same call in `pantometry-gpu`.
                    let (rgba, _) = pixels.as_chunks::<4>();
                    let lit = rgba.iter().filter(|p| p[..3] != background).count();
                    let mut shades: Vec<[u8; 3]> = Vec::new();
                    for p in rgba {
                        let c = [p[0], p[1], p[2]];
                        if c != background && !shades.contains(&c) {
                            shades.push(c);
                        }
                    }
                    // **The image's own size, not the default one.** `lit` is counted over
                    // `pixels`, which under `--thumbnail` is 240x156 and not 1100x720 — so a tile
                    // lighting 173 of its 37 440 was reported as "173 of 792000 (0.02%)", a fifth of
                    // a percent read as a fiftieth. The count was right and the denominator was a
                    // constant from before there was a second size.
                    let all = pw as u64 * ph as u64;
                    if animate {
                        println!(
                        "  frame {at}: {lit} of {all} pixels carry a line ({:.2}%), in {} shades",
                        100.0 * lit as f64 / all as f64,
                        shades.len()
                    );
                    } else {
                        println!(
                        "  wrote {out} — {lit} of {all} pixels carry a line ({:.2}%), in {} shades",
                        100.0 * lit as f64 / all as f64,
                        shades.len()
                    );
                    }
                    size = (pw, ph);
                    // **What is in the picture, in numbers.** A shaded solid says where the hot end is
                    // and cannot say whether it is 119 °C or 1190; the size of a thing on screen is
                    // whatever the camera chose; and a frame index is not a time. All three are known
                    // here and none of them was printed, so a person with the file still had to open
                    // the run to read its own picture.
                    // **Per frame, because an empty one in the middle is the failure that hides.**
                    // A run whose renderer stopped working on frame 31 writes thirty good pictures
                    // and seventeen of nothing, and a check at the end sees only the last.
                    if lit == 0 {
                        eprintln!("  nothing was drawn");
                        std::process::exit(1);
                    }
                }
                Err(e) => {
                    eprintln!("  no GPU available for a snapshot: {e}");
                    std::process::exit(3);
                }
            }
        }
        if animate {
            write_gif(&out, size.0, size.1, &collected, GIF_CENTISECONDS);
        }
        app.describe();
        return 0;
    }

    let event_loop = EventLoop::new().expect("an event loop");
    event_loop.set_control_flow(ControlFlow::Poll);
    let mut app = App::new(run, panel);
    app.frame_on(framed);
    event_loop.run_app(&mut app).expect("the window runs");
    0
}

/// Everything the window holds between frames.
struct App {
    run: Run,
    panel: String,
    /// The panel the window opened on, which the colour bar returns to when the framing does.
    home: String,
    /// The panel the camera is framed on, or `None` for the whole run.
    ///
    /// **One framing for everything drawn**, whichever this is: every panel is projected through
    /// `framing`, so framing one panel moves nothing relative to anything else, and what is far
    /// from it leaves the screen. The scale bar measures the framed box, because that is the box
    /// the camera was fitted to.
    framed: Option<String>,
    framing: Framing,
    camera: Camera,
    /// The run-wide range the shading is measured against.
    span: (f64, f64),
    frame: usize,
    /// Whether a legend is drawn over the picture.
    ///
    /// **Off for a tile.** At 240x156 the labels are a fifth of the size they are laid out for and
    /// read as speckle, and worse: the bar spans the full width, so `crop_to_content` can no
    /// longer find the object and every tile is cropped to the legend instead. A tile's job is to
    /// tell one scene from another at a glance, and the chooser prints the title and the kinds
    /// beside it.
    legend: bool,
    playing: bool,
    dragging: Option<(f64, f64)>,
    gpu: Option<Gpu>,
}

/// `out.png` and frame 7 give `out-007.png`, so a sequence sorts and globs in frame order.
///
/// **The dot has to be after the last separator.** `../docs/out.png` holds dots in its `..` as
/// well, and splitting on the last one anywhere would have written `.` beside `./docs/out-007`
/// — a file in the wrong directory under a name nothing would look for. Three digits because a
/// run with a thousand frames is a different problem than this.
fn numbered(path: &str, frame: usize) -> String {
    let after = path.rfind(['/', '\\']).map_or(0, |i| i + 1);
    match path[after..].rfind('.') {
        Some(dot) => {
            let at = after + dot;
            format!("{}-{frame:03}{}", &path[..at], &path[at..])
        }
        None => format!("{path}-{frame:03}"),
    }
}

/// The size to draw a colour bar's two end labels at, so they do not run through each other.
///
/// `ends` is their combined width in glyph-lattice units, `bar` the bar's width in normalised
/// device coordinates, `em` the legend's own size and `aspect` the canvas's.
///
/// **Shrinks rather than drops.** The editor's colour bar carries three labels and can lose its
/// middle one; a bar has only two ends and a bar showing one of them is not a scale. Never grows:
/// a short pair is drawn at the legend's size like everything else, because a bar whose numbers
/// changed size with their own length would read as two different kinds of label.
///
/// One character of clear air between them, so the two read as two.
fn ends_that_fit(ends: f32, bar: f32, em: f32, aspect: f32) -> f32 {
    let want = (ends + crate::glyphs::width("0")) * em / aspect;
    if want <= bar || want <= 0.0 {
        em
    } else {
        em * bar / want
    }
}

/// The shortest the scale bar is drawn, in normalised device `x`: a twelfth of the width.
///
/// Rounding *up* a 1-2-5 ladder from here gives at most 2.5 times this, so the bar is always
/// between 0.167 and 0.417 of the half-width -- the range it had when it was a fraction of the
/// subject, now earned on the screen rather than assumed.
const SHORTEST_BAR: f64 = 1.0 / 6.0;

/// A round number of metres, and how much of the screen's half-width it covers.
///
/// **It takes the camera's metres-to-screen and nothing else**, which is
/// [`viewer_core::Camera::across_per_metre`]: two points a known distance apart at the framing's
/// centre, square to the view, projected through the same `Camera::project` every vertex goes
/// through. So the bar is the length of what it names *in the picture*, and follows a zoom.
///
/// **It assumed half the screen was half the subject**, and it was not. That ignored the focal
/// length `Camera::fit` chooses — the whole point of which is that the subject does *not* land at
/// a fixed fraction — and the perspective at the distance the camera stands. Measured framed on
/// scene 14's 20 mm bar: the bar 1018 px long on screen, the `10 MM` ruler under it 137 px, 3.7
/// times short. The test beside this encoded the same identity, so it agreed with it.
///
/// **And it still cannot be handed a frame's box.** The camera is fitted once over every frame of
/// the run, which is what fixed a ruler that ran 108 px at the ends of the protein's swing and 119
/// in the middle under one `2 NM`; it takes the camera now, which is fitted once, rather than the
/// box, which was.
///
/// `None` when there is nothing to measure.
fn scale_bar(per_metre: f64) -> Option<(f64, f32)> {
    if !per_metre.is_finite() || per_metre <= 0.0 {
        return None;
    }
    let raw = SHORTEST_BAR / per_metre;
    if !raw.is_finite() || raw <= 0.0 {
        return None;
    }
    let decade = 10f64.powf(raw.log10().floor());
    let nice = [1.0, 2.0, 5.0, 10.0]
        .into_iter()
        .map(|m| m * decade)
        .find(|n| *n >= raw)
        .unwrap_or(decade);
    Some((nice, (nice * per_metre) as f32))
}

/// The scale bar's label: a round number of metres with the SI prefix that number is near.
///
/// **Plain formatting, not `magnitude`.** The length is a round number by construction — one, two
/// or five times a decade — and the general formatter printed `5.0000 MM` for 5.
///
/// **And more than two prefixes, which is what it had.** Metres and millimetres, formatted to
/// three decimals: everything below a micrometre came out as `0 MM`. A protein 4.5 nm across is
/// the first thing this workspace has drawn small enough to show that, and it showed it on the
/// figure going onto the front page — a legend reporting zero length for a real object, beside a
/// library whose front page is about dimensions living in the type system.
fn scale_label(metres: f64) -> String {
    let round = |v: f64| {
        let s = format!("{v:.3}");
        let s = s.trim_end_matches('0').trim_end_matches('.').to_string();
        if s.is_empty() {
            "0".to_string()
        } else {
            s
        }
    };
    // Micro is `U`, because the glyph lattice has no Greek and a letter every reader of a plot has
    // seen stand in for it beats the box an unknown character now draws.
    for (floor, per_metre, unit) in [
        (1e3, 1e-3, "KM"),
        (1.0, 1.0, "M"),
        (1e-3, 1e3, "MM"),
        (1e-6, 1e6, "UM"),
        (1e-9, 1e9, "NM"),
    ] {
        if metres >= floor {
            return format!("{} {unit}", round(metres * per_metre));
        }
    }
    format!("{} PM", round(metres * 1e12))
}

struct Gpu {
    window: Arc<Window>,
    surface: wgpu::Surface<'static>,
    device: wgpu::Device,
    queue: wgpu::Queue,
    config: wgpu::SurfaceConfiguration,
    /// Triangles, for the solid a field is.
    solid: wgpu::RenderPipeline,
    /// Lines, for a `Paths` panel, a body's cross and the box.
    pipeline: wgpu::RenderPipeline,
    /// Rebuilt whenever the surface is, because it has to match its size.
    depth: wgpu::TextureView,
}

impl App {
    fn new(run: Run, panel: String) -> App {
        let mut app = App {
            run,
            home: panel.clone(),
            panel,
            framed: None,
            framing: Framing::of([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0]),
            span: (0.0, 1.0),
            camera: Camera::default(),
            legend: true,
            frame: 0,
            playing: true,
            dragging: None,
            gpu: None,
        };
        app.frame_on(None);
        app
    }

    /// The box the camera frames: every panel's, or one panel's, over **every frame**.
    ///
    /// **The box every panel occupies, not the first one's**, by default. A run can hold several
    /// and this shell draws all of them; framing on one unasked would put the others off the
    /// screen, which is how a picture loses something without saying it has. Framing on one is
    /// something a person asks for, with `--frame-panel` or `F`.
    ///
    /// Over every frame for the reason `Run::framing_of` gives: a camera framed to the current
    /// frame follows a moving body and makes it look still. This is that function's loop with a
    /// filter on it, kept on each panel's own box exactly as that one is.
    fn box_of(run: &Run, which: Option<&str>) -> Option<[f64; 6]> {
        let mut out = [f64::MAX, f64::MAX, f64::MAX, f64::MIN, f64::MIN, f64::MIN];
        for name in run.panels() {
            if which.is_some_and(|w| w != name) {
                continue;
            }
            if let Some(b) = run.framing_of(&name) {
                for a in 0..3 {
                    out[a] = out[a].min(b[a]);
                    out[a + 3] = out[a + 3].max(b[a + 3]);
                }
            }
        }
        (out[0] <= out[3]).then_some(out)
    }

    /// Frame the camera on one panel, or on the whole run with `None`.
    ///
    /// Sets the framing and the fit together, from one box, so the two cannot describe different
    /// things; the scale bar is read off the camera this fits, so it cannot either. Framing a
    /// panel also hands it the colour bar, because a reader who asked to look at the bar is
    /// asking about the bar's values; the whole run hands it back to the panel the window opened
    /// on.
    fn frame_on(&mut self, which: Option<String>) {
        let whole = App::box_of(&self.run, which.as_deref())
            .or_else(|| App::box_of(&self.run, None))
            .unwrap_or([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0]);
        self.framing = Framing::of(whole);
        self.panel = which.clone().unwrap_or_else(|| self.home.clone());
        self.framed = which;
        // Once, from the whole run. Re-fitting it per frame is what `Run::scale_of` exists to
        // stop, and for a while nothing called it.
        self.span = self.run.scale_of(&self.panel).unwrap_or((0.0, 1.0));
        // Framed to what is actually there. The window can still be zoomed; `--snapshot` cannot,
        // and a fixed distance is a distance chosen for a cube. The angles are kept: changing
        // what is framed is not turning to look at it from somewhere else.
        self.camera.fit(whole, &self.framing, 16.0 / 9.0, 0.85);
    }

    /// The scale bar: a round number of metres, and how much of the half-width it covers.
    ///
    /// Read off the camera and the framing as they are now, at the aspect being drawn, so a zoom
    /// or a framed panel moves it with the picture.
    fn ruler(&self, aspect: f64) -> Option<(f64, f32)> {
        scale_bar(self.camera.across_per_metre(&self.framing, aspect)?)
    }

    /// The framing `F` steps to next: the whole run, then each panel in the order they appear,
    /// then the whole run again.
    fn next_framing(&self) -> Option<String> {
        let names = self.run.panels();
        match &self.framed {
            None => names.first().cloned(),
            Some(now) => names
                .iter()
                .position(|n| n == now)
                .and_then(|i| names.get(i + 1))
                .cloned(),
        }
    }

    /// The legend, drawn into the picture: a colour bar with its ends and unit, a scale bar in
    /// metres, and the time.
    ///
    /// **A picture nobody can read a value off is a gradient.** The snapshot is the only static
    /// figure this workspace makes, and it had no legend, no scale and no units — a reader could
    /// see where the hot end was and not whether it was 119 °C or 1190, nor whether the thing was
    /// 12 mm or 4 m across.
    ///
    /// In screen space at the near plane, so it sits over the solid whatever the camera does. The
    /// colours are `editor_core::value_colour` over `bar_value`, which is the editor's own bar, so
    /// the two legends cannot disagree about which colour a number is.
    fn legend(&self, aspect: f64) -> (Vec<Vertex>, Vec<Vertex>) {
        let Some(panel) = self
            .run
            .frames
            .get(self.frame)
            .and_then(|f| f.panels.iter().find(|p| p.name() == self.panel))
        else {
            return (Vec::new(), Vec::new());
        };
        let (lo, hi) = self.span;
        // **A flat colour range skips the colour bar and nothing else.** This returned an *empty
        // legend* on `hi <= lo` — no scale bar, no clock, nothing — and the reasoning only ever
        // held for the colour bar: a run whose values never change has no scale to draw, and has
        // a size and a time like any other run. What it looked like was a picture with no ruler,
        // which is the thing the legend exists to stop.
        let shading = lo.is_finite() && hi.is_finite() && hi > lo;
        let (mut tris, mut lines) = (Vec::new(), Vec::new());

        // Everything here is normalised device coordinates: x and y run -1..1, and the depth is
        // the near end so nothing in the scene can cover it.
        const NEAR: f32 = 0.001;
        let ink = [0.82f32, 0.85, 0.88];
        // One lattice unit, sized so a label is legible at 1100 wide and still is at 240.
        let em = 0.011f32;
        // How tall one line of it is, from the font rather than from a number that clears it
        // today. Every label below something is placed at `-= row + GAP`, so changing `HEIGHT`
        // moves the layout instead of overlapping what the label is for.
        let row = crate::glyphs::HEIGHT * em;
        const GAP: f32 = 0.024;
        let quad = |v: &mut Vec<Vertex>, x0: f32, y0: f32, x1: f32, y1: f32, c: [f32; 3]| {
            for (x, y) in [(x0, y0), (x1, y0), (x1, y1), (x0, y0), (x1, y1), (x0, y1)] {
                v.push(Vertex {
                    position: [x, y, NEAR],
                    colour: c,
                });
            }
        };
        let write = |v: &mut Vec<Vertex>, s: &str, x: f32, y: f32| {
            for (ax, ay, bx, by) in crate::glyphs::text(s) {
                v.push(Vertex {
                    position: [x + ax * em / aspect as f32, y + ay * em, NEAR],
                    colour: ink,
                });
                v.push(Vertex {
                    position: [x + bx * em / aspect as f32, y + by * em, NEAR],
                    colour: ink,
                });
            }
        };

        // **The colour bar**, bottom right, in the same steps the editor uses.
        let scale = Some((lo, hi));
        if shading {
            // Clear of the bottom edge: a label whose descender is off the canvas is a label
            // a reader distrusts, and the first attempt put both rows there.
            let (bx0, bx1, by0, by1) = (0.34f32, 0.94, -0.80, -0.755);
            let steps = 96;
            for i in 0..steps {
                let u = i as f64 / (steps - 1) as f64;
                let [r, g, b] = editor_core::value_colour(editor_core::bar_value(u, scale), scale);
                let x0 = bx0 + (bx1 - bx0) * i as f32 / steps as f32;
                let x1 = bx0 + (bx1 - bx0) * (i + 1) as f32 / steps as f32;
                quad(
                    &mut tris,
                    x0,
                    by0,
                    x1,
                    by1,
                    [
                        (r as f32 / 255.0).powf(2.2),
                        (g as f32 / 255.0).powf(2.2),
                        (b as f32 / 255.0).powf(2.2),
                    ],
                );
            }
            let unit = panel.unit();
            // **Both ends, at a size that fits between them.** They went in at the legend's own
            // size with no regard for how wide they are: `3.9674E-11` and `9.5054E-11` want 0.79
            // of the canvas at the two ends of a bar 0.60 wide, and came out as
            // `3.9674E-5054E-11` — two numbers through each other, which a reader cannot tell
            // from one. The editor's colour bar had the same defect and drops its *middle* label;
            // here there is no middle to drop, because a bar showing one end is not a scale.
            let (los, his) = (editor_core::magnitude(lo), editor_core::magnitude(hi));
            let ends = crate::glyphs::width(&los) + crate::glyphs::width(&his);
            let fits = ends_that_fit(ends, bx1 - bx0, em, aspect as f32);
            let small = |v: &mut Vec<Vertex>, s: &str, x: f32, y: f32| {
                for (ax, ay, bx, by) in crate::glyphs::text(s) {
                    v.push(Vertex {
                        position: [x + ax * fits / aspect as f32, y + ay * fits, NEAR],
                        colour: ink,
                    });
                    v.push(Vertex {
                        position: [x + bx * fits / aspect as f32, y + by * fits, NEAR],
                        colour: ink,
                    });
                }
            };
            small(&mut lines, &los, bx0, by0 - row - GAP);
            small(
                &mut lines,
                &his,
                bx1 - crate::glyphs::width(&his) * fits / aspect as f32,
                by0 - row - GAP,
            );
            // **Centred on the bar, and kept inside the frame.** The unit went in at the bar's
            // midpoint used as a *left* edge, so a long one started halfway along the bar and ran off
            // the canvas: the front page's bench render read `REFRACTIVE` and then the edge of the
            // picture. That was true before the glyph table could spell the word, which is why
            // nobody read it as clipping — it looked like the holes everything else had.
            let span = |s: &str| crate::glyphs::width(s) * em / aspect as f32;
            let left = ((bx0 + bx1) * 0.5 - span(unit) * 0.5)
                .min(0.98 - span(unit))
                .max(-0.98);
            write(&mut lines, unit, left, by1 + 0.02);
        }

        // **The scale bar**, bottom left: a round number of metres across the object, so a reader
        // knows whether they are looking at a die or a room.
        //
        // **From the camera, and not from a box.** It was the run's longest side, on the
        // assumption that half the screen is half the subject, which the fit's focal length makes
        // false — 3.7 times short framed on scene 14's bar. Before that it was this frame's box,
        // which the camera is not fitted to either. `scale_bar` takes how far a metre at the
        // framing's centre reaches on this screen, which is the only number a ruler can be.
        if let Some((nice, across)) = self.ruler(aspect) {
            let (sx, sy) = (-0.94f32, -0.79);
            quad(&mut tris, sx, sy, sx + across, sy + 0.008, ink);
            for x in [sx, sx + across] {
                quad(&mut tris, x, sy - 0.012, x + 0.004, sy + 0.02, ink);
            }
            write(&mut lines, &scale_label(nice), sx, sy - row - GAP);
        }

        // The time, top left, because a frame index is not one.
        let t = self.run.frames.get(self.frame).map_or(0.0, |f| f.t);
        write(
            &mut lines,
            &format!("T {} S", editor_core::magnitude(t)),
            -0.94,
            1.0 - row - 0.02,
        );

        // **What the run measured, under the clock.** A run carries `readings` — the numbers the
        // simulation took, frame by frame — and this file did not contain the word. So a window
        // showing a protein closing over a molecule could not say how far it had closed, with the
        // number sitting in the file it had open. The editor drew them and the viewer dropped
        // them, which is the shape of a picture that is a gradient rather than a reading.
        //
        // Smaller than the rest of the legend, because a label like `the furthest residue has
        // moved` is 38 characters and at the legend's own size that is three quarters of the
        // canvas. Measured, and the value is right-aligned to a column so the digits line up
        // rather than wandering with the label.
        let small = em * 0.72;
        let readings = self
            .run
            .frames
            .get(self.frame)
            .map_or(&[][..], |f| &f.readings);
        // **The column is the widest label in the frame, not a constant.** It was -0.20, chosen
        // by eye, and `the furthest residue has moved` runs past it — so its value landed on top
        // of its own label and the picture read `M(+)VE(+)77(+)` where a number should be. The
        // legend knows how wide a string is; it had not been asked.
        let wide = |s: &str| crate::glyphs::width(s) * small / aspect as f32;
        let widest_value = readings
            .iter()
            .map(|r| wide(&format!("{} {}", editor_core::magnitude(r.value), r.unit)))
            .fold(0.0f32, f32::max);
        // **Clamped so the value stays on the canvas.** A column placed purely by the label runs
        // the number off the right-hand edge as soon as a label is long enough, and a number that
        // is not there looks exactly like a number that is short. If the two cannot both fit they
        // overlap, which is ugly and visible; going off the edge is neither.
        let column = (-0.94
            + readings
                .iter()
                .map(|r| wide(&r.label.to_ascii_uppercase()))
                .fold(0.0f32, f32::max)
            + wide("  "))
        .min(0.94 - widest_value);
        let mut line = 1.0 - row - 0.02 - crate::glyphs::HEIGHT * small - GAP;
        // **A label longer than its column is cut, and says it was.** Sixty-nine characters runs
        // past the right-hand edge of the canvas, where the part that is missing looks like a
        // label that was simply shorter. Two dots is the difference between a legend that lost
        // something and a legend that says it did.
        let room = column - (-0.94) - wide(" ");
        let fitted = |s: &str| {
            if wide(s) <= room {
                return s.to_string();
            }
            let mut cut = String::new();
            for c in s.chars() {
                if wide(&format!("{cut}{c}..")) > room {
                    break;
                }
                cut.push(c);
            }
            format!("{cut}..")
        };
        for r in readings {
            let label = fitted(&r.label.to_ascii_uppercase());
            let value = format!("{} {}", editor_core::magnitude(r.value), r.unit);
            for (s, x) in [(&label, -0.94f32), (&value, column)] {
                for (ax, ay, bx, by) in crate::glyphs::text(s) {
                    lines.push(Vertex {
                        position: [x + ax * small / aspect as f32, line + ay * small, NEAR],
                        colour: ink,
                    });
                    lines.push(Vertex {
                        position: [x + bx * small / aspect as f32, line + by * small, NEAR],
                        colour: ink,
                    });
                }
            }
            line -= crate::glyphs::HEIGHT * small + GAP * 0.5;
        }

        // **How to drive it, in the window.** These went to stdout, where a person looking at the
        // picture is not. It is most of what there is to know about using this, and a screenshot
        // that does not carry it teaches nothing about the application it is a screenshot of.
        //
        // Top right, in two lines, because the corners are taken: the scale bar is bottom left and
        // the colour bar bottom right. One line of all four is 52 characters, which at this size
        // is three quarters of the canvas and ran off the edge — the first attempt drew it across
        // both bars and lost `ARROWS SCRUB` past the right-hand side.
        //
        // **A third line only while one panel is framed**, saying which: a picture of one panel
        // with the rest gone off-screen is otherwise indistinguishable from a run that only had
        // one. Not drawn for the whole run, so a whole-run picture is the one it always was and
        // the committed figures `docs/README.md` refreshes by command are still true of it.
        let framed = self
            .framed
            .as_ref()
            .map(|name| format!("FRAMED ON {}   F NEXT", name.to_ascii_uppercase()));
        let hints = ["DRAG ROTATE   SCROLL ZOOM", "SPACE PLAY   ARROWS SCRUB"];
        for (k, hint) in hints.into_iter().chain(framed.as_deref()).enumerate() {
            let at = 0.94 - wide(hint);
            let y = 1.0 - row - 0.02 - (crate::glyphs::HEIGHT * small + GAP * 0.5) * k as f32;
            for (ax, ay, bx, by) in crate::glyphs::text(hint) {
                lines.push(Vertex {
                    position: [at + ax * small / aspect as f32, y + ay * small, NEAR],
                    colour: ink,
                });
                lines.push(Vertex {
                    position: [at + bx * small / aspect as f32, y + by * small, NEAR],
                    colour: ink,
                });
            }
        }
        (tris, lines)
    }

    /// Print what is in the picture, in numbers.
    ///
    /// **A shaded solid says where the hot end is and cannot say how hot.** The snapshot reported
    /// how many pixels it had lit and nothing else, so a reader with the file could not tell 119 °C
    /// from 1190, could not tell a 12 mm module from a 4 m room, and had a frame index where a time
    /// belongs. All three are known here.
    ///
    /// Printed rather than drawn, which is the cheap half of the answer: the colour bar and the
    /// scale go **into** the picture beside this, and a line of text is what a person pastes into a
    /// note.
    fn describe(&self) {
        let Some(panel) = self
            .run
            .frames
            .get(self.frame)
            .and_then(|f| f.panels.iter().find(|p| p.name() == self.panel))
        else {
            return;
        };
        let (lo, hi) = self.span;
        println!(
            "  {} spans {} to {} {} over the run",
            panel.name(),
            editor_core::magnitude(lo),
            editor_core::magnitude(hi),
            panel.unit()
        );
        // The box in metres, from the panel's own extent rather than from the camera — what is on
        // screen is whatever the framing chose, and this is the thing itself.
        let b = panel.world_bounds();
        let (dx, dy, dz) = (b[3] - b[0], b[4] - b[1], b[5] - b[2]);
        if dx.max(dy).max(dz) > 0.0 {
            println!(
                "  {} x {} x {} m, at t = {} s",
                editor_core::magnitude(dx),
                editor_core::magnitude(dy),
                editor_core::magnitude(dz),
                editor_core::magnitude(self.run.frames.get(self.frame).map_or(0.0, |f| f.t))
            );
        }
    }

    /// The triangles and the lines for this frame, already projected.
    ///
    /// **A field is a solid and is drawn as one.** This shell drew every panel as line segments —
    /// a body and a field sample each became a small cross — and a block of cells came out as a
    /// cloud of `+` glyphs with no surface, no shading and no way to tell which way is up. That is
    /// the same criticism `render.rs` was written to answer for the editor's viewport, and it
    /// says so in its own first paragraph; the shell a person reaches by typing `pantometry view`
    /// kept the old picture. Measured on `24-a-power-module-junction-to-ambient`, a stack of four
    /// materials: 5 045 lit pixels **in eight shades**, and nothing in it said where the die was.
    ///
    /// The geometry is `editor_core::field_shell`, which is `pantometry_view::mesh` — the same
    /// function that writes glTF and USD and that the editor shades. One picture, not three, and
    /// a size that cannot disagree with an export.
    fn vertices(&self, aspect: f64) -> (Vec<Vertex>, Vec<Vertex>) {
        let Some(frame) = self.run.frames.get(self.frame) else {
            return (Vec::new(), Vec::new());
        };
        let mut tris = Vec::new();
        let mut out: Vec<Vertex> = Vec::new();
        // **Every panel.** This drew `self.panel` and nothing else -- "the first one there is" --
        // which was fair while every shape was lines and a run held one domain worth drawing. A
        // bench whose glass is a solid and whose rays are paths is two panels, and showing one of
        // them is a picture with the light missing and no note to say so. The legend and the
        // readout still name the selected panel, because a colour bar belongs to one quantity.
        for panel in &frame.panels {
            let (t, l) = self.panel_vertices(panel, aspect);
            tris.extend(t);
            out.extend(l);
        }
        self.assemble(tris, out, aspect)
    }

    /// One panel's triangles and lines.
    fn panel_vertices(
        &self,
        panel: &viewer_core::Panel,
        aspect: f64,
    ) -> (Vec<Vertex>, Vec<Vertex>) {
        let mut tris = Vec::new();
        // **A solid comes with its triangles.** The only producer of a lit surface here was
        // `field_shell`, an isosurface through somebody's grid, so an instrument drawn as its own
        // surfaces arrived as lines and left as lines however carefully it had been meshed.
        // **Each panel on its own scale.** `self.span` is the selected panel's, which was the
        // only panel there was; with two in a frame it painted a field angle as a refractive
        // index and the rays came out the colours of glass. The legend still belongs to the
        // selected one, because a colour bar names one quantity and there is room for one bar.
        let span = self.run.scale_of(panel.name()).unwrap_or((0.0, 1.0));
        // **Whether a point is between the planes a GPU would clip at.** Everything here is
        // projected on the CPU and reaches the pass already in screen space, so nothing clips it:
        // `Camera::project` clamps a point behind the eye to the near plane and hands back an
        // enormous coordinate, and a triangle with such a corner is drawn across the window. That
        // never happened while the framing was the whole run, because a fitted run is inside both
        // planes. Framed on one panel of a run with several scales, nearly everything else is
        // outside them. A primitive with any corner outside is dropped rather than cut, which
        // loses the part of it that straddles a plane; nothing that does is the thing framed.
        let seen = |p: [f64; 3]| {
            self.camera
                .sees(&self.camera.project(p, &self.framing, aspect))
        };
        if let viewer_core::Panel::Surface { values, .. } = panel {
            let placed = panel.placed_surface_points();
            let faces = panel.surface_faces();
            let (lo, hi) = span;
            let mut normals = vec![[0.0f64; 3]; placed.len()];
            for t in &faces {
                let (a, b, c) = (placed[t[0]], placed[t[1]], placed[t[2]]);
                let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
                let n = [
                    u[1] * v[2] - u[2] * v[1],
                    u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0],
                ];
                for &i in t {
                    for (slot, add) in normals[i].iter_mut().zip(n) {
                        *slot += add;
                    }
                }
            }
            // The same key light the field's solid takes, so two solids in one frame are lit
            // alike rather than by two conventions.
            let key = [0.35f32, 0.45, 0.82];
            for t in &faces {
                // Dropped whole when any corner is outside the clip planes -- see `seen` below.
                if !t.iter().all(|&i| seen(placed[i])) {
                    continue;
                }
                for &i in t {
                    let c = self.camera.project(placed[i], &self.framing, aspect);
                    let n = normals[i];
                    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
                    let unit_n = if len > 0.0 {
                        [
                            (n[0] / len) as f32,
                            (n[1] / len) as f32,
                            (n[2] / len) as f32,
                        ]
                    } else {
                        [0.0, 1.0, 0.0]
                    };
                    let lit = (unit_n[0] * key[0] + unit_n[1] * key[1] + unit_n[2] * key[2])
                        .abs()
                        .mul_add(0.65, 0.35);
                    let v = values.get(i).copied().unwrap_or(lo);
                    let base = ramp(((v - lo) / (hi - lo).max(1e-30)).clamp(0.0, 1.0));
                    tris.push(Vertex {
                        position: [c.x as f32, c.y as f32, c.depth as f32],
                        colour: [base[0] * lit, base[1] * lit, base[2] * lit],
                    });
                }
            }
        }
        if let viewer_core::Panel::Field {
            nx,
            ny,
            nz,
            unit,
            lattice,
            values,
            ..
        } = panel
        {
            // Two axes at least, or there is no surface: a row of samples along a line is a graph,
            // and the cross pipeline below is the honest picture of it.
            if [*nx, *ny, *nz].iter().filter(|&&n| n > 1).count() >= 2 {
                if let Some(corners) = panel.placed_corners() {
                    let shell = editor_core::field_shell(
                        &corners,
                        (*nx, *ny, *nz),
                        *lattice,
                        values,
                        unit,
                        Some(span),
                        None,
                    );
                    // **Lambert on the CPU**, because the camera already projects here and the
                    // shader takes a colour rather than a light. One key light down the eye and a
                    // floor under it, which is what keeps a face turned away from being black and
                    // a solid from reading as a silhouette.
                    let key = [0.35f32, 0.45, 0.82];
                    let kept: Vec<u32> = shell
                        .indices
                        .as_chunks::<3>()
                        .0
                        .iter()
                        .filter(|t| t.iter().all(|&i| seen(shell.positions[i as usize])))
                        .flatten()
                        .copied()
                        .collect();
                    for i in &kept {
                        let i = *i as usize;
                        let p = shell.positions[i];
                        let c = self.camera.project(p, &self.framing, aspect);
                        let n = shell.normals.get(i).copied().unwrap_or([0.0, 0.0, 1.0]);
                        let lit = (n[0] * key[0] + n[1] * key[1] + n[2] * key[2])
                            .abs()
                            .mul_add(0.65, 0.35);
                        let base = shell.colours.get(i).copied().unwrap_or([1.0, 1.0, 1.0]);
                        tris.push(Vertex {
                            position: [c.x as f32, c.y as f32, c.depth as f32],
                            colour: [base[0] * lit, base[1] * lit, base[2] * lit],
                        });
                    }
                }
            }
        }

        let mut out: Vec<Vertex> = Vec::new();
        // The lines are what is left over: a `Paths` panel, and a point set, which has no surface
        // to build. A field that got a solid above skips them, or every cell would carry a cross
        // inside the block that is drawn over it.
        if tris.is_empty() {
            for s in segments(panel, &self.camera, &self.framing, aspect, span) {
                if !self.camera.sees(&s.from) || !self.camera.sees(&s.to) {
                    continue;
                }
                let colour = ramp(s.shade);
                out.push(Vertex {
                    position: [s.from.x as f32, s.from.y as f32, s.from.depth as f32],
                    colour,
                });
                out.push(Vertex {
                    position: [s.to.x as f32, s.to.y as f32, s.to.depth as f32],
                    colour,
                });
            }
        }
        (tris, out)
    }

    /// Both sets, mapped into the depth range wgpu keeps.
    ///
    /// Split out when this shell began drawing every panel rather than one: the mapping has to
    /// see the whole frame, so doing it per panel would have put each panel's own near face at
    /// `z = 0` and stacked them in the order they were listed rather than where they are.
    fn assemble(
        &self,
        tris: Vec<Vertex>,
        out: Vec<Vertex>,
        aspect: f64,
    ) -> (Vec<Vertex>, Vec<Vertex>) {
        let (mut tris, mut out) = (tris, out);
        // **`Projected::depth` is a distance from the eye, not a clip `z`.** wgpu keeps
        // `0 <= z <= 1` and discards the rest, so handing it metres drew nothing at all — the
        // first run of this rendered an empty frame and said so. Mapped over the range this frame
        // actually spans, and over **both** sets together: two mappings would put every line
        // either in front of or behind every triangle regardless of where it is.
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for v in tris.iter().chain(out.iter()) {
            let z = v.position[2] as f64;
            if z.is_finite() {
                lo = lo.min(z);
                hi = hi.max(z);
            }
        }
        let span = if hi > lo { hi - lo } else { 1.0 };
        // A hair inside the range at each end, so the nearest surface is not on the clip plane
        // and the furthest is not lost to the `LessEqual` compare against the cleared 1.0.
        for v in tris.iter_mut().chain(out.iter_mut()) {
            let z = v.position[2] as f64;
            v.position[2] = if z.is_finite() {
                (0.01 + 0.98 * (z - lo) / span) as f32
            } else {
                0.99
            };
        }
        // The legend last, at the near plane, so it is over whatever the camera is looking at.
        if self.legend {
            let (mut bar, mut labels) = self.legend(aspect);
            tris.append(&mut bar);
            out.append(&mut labels);
        }
        (tris, out)
    }
}

impl ApplicationHandler for App {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        if self.gpu.is_some() {
            return;
        }
        let attrs = Window::default_attributes()
            .with_title(format!("pantometry — {}", self.run.title))
            .with_inner_size(winit::dpi::LogicalSize::new(1100.0, 720.0));
        let window = Arc::new(event_loop.create_window(attrs).expect("a window"));

        let instance = wgpu::Instance::default();
        let surface = instance
            .create_surface(window.clone())
            .expect("a surface on that window");
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: Some(&surface),
            force_fallback_adapter: false,
        }))
        .expect("a GPU adapter");
        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("pantometry viewer"),
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::downlevel_defaults(),
                memory_hints: wgpu::MemoryHints::default(),
            },
            None,
        ))
        .expect("a device");

        let size = window.inner_size();
        let caps = surface.get_capabilities(&adapter);
        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format: caps.formats[0],
            width: size.width.max(1),
            height: size.height.max(1),
            present_mode: wgpu::PresentMode::Fifo,
            alpha_mode: caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);

        let solid = pipeline_for(
            &device,
            config.format,
            wgpu::PrimitiveTopology::TriangleList,
        );
        let pipeline = pipeline_for(&device, config.format, wgpu::PrimitiveTopology::LineList);
        let depth = depth_texture(&device, config.width, config.height);

        self.gpu = Some(Gpu {
            window,
            surface,
            device,
            queue,
            config,
            solid,
            pipeline,
            depth,
        });
    }

    fn window_event(&mut self, event_loop: &ActiveEventLoop, _id: WindowId, event: WindowEvent) {
        let Some(gpu) = self.gpu.as_mut() else {
            return;
        };
        match event {
            WindowEvent::CloseRequested => event_loop.exit(),
            WindowEvent::Resized(size) => {
                gpu.config.width = size.width.max(1);
                gpu.config.height = size.height.max(1);
                gpu.surface.configure(&gpu.device, &gpu.config);
                // The depth buffer is an attachment and has to match the colour one's size, or
                // the pass is refused. A resize that rebuilt one and not the other would take the
                // window down on the next frame.
                gpu.depth = depth_texture(&gpu.device, gpu.config.width, gpu.config.height);
            }
            WindowEvent::MouseInput { state, .. } => {
                self.dragging = matches!(state, ElementState::Pressed).then_some((0.0, 0.0));
            }
            WindowEvent::CursorMoved { position, .. } => {
                if let Some((lx, ly)) = self.dragging {
                    let (x, y) = (position.x, position.y);
                    if lx != 0.0 || ly != 0.0 {
                        self.camera.turn((x - lx) * 0.008, (y - ly) * 0.006);
                    }
                    self.dragging = Some((x, y));
                }
            }
            WindowEvent::MouseWheel { delta, .. } => {
                let step = match delta {
                    MouseScrollDelta::LineDelta(_, y) => y as f64,
                    MouseScrollDelta::PixelDelta(p) => p.y / 60.0,
                };
                self.camera.zoom(1.0 - step * 0.1);
            }
            WindowEvent::KeyboardInput { event, .. } if event.state == ElementState::Pressed => {
                let n = self.run.frames.len().max(1);
                match event.physical_key {
                    PhysicalKey::Code(KeyCode::Space) => self.playing = !self.playing,
                    PhysicalKey::Code(KeyCode::ArrowRight) => {
                        self.playing = false;
                        self.frame = (self.frame + 1) % n;
                    }
                    PhysicalKey::Code(KeyCode::ArrowLeft) => {
                        self.playing = false;
                        self.frame = (self.frame + n - 1) % n;
                    }
                    PhysicalKey::Code(KeyCode::KeyF) => {
                        let next = self.next_framing();
                        self.frame_on(next);
                        match &self.framed {
                            Some(name) => println!("  framed on {name}"),
                            None => println!("  framed on the whole run"),
                        }
                    }
                    PhysicalKey::Code(KeyCode::Escape) => event_loop.exit(),
                    _ => {}
                }
            }
            WindowEvent::RedrawRequested => self.draw(),
            _ => {}
        }
    }

    fn about_to_wait(&mut self, _event_loop: &ActiveEventLoop) {
        if let Some(gpu) = &self.gpu {
            gpu.window.request_redraw();
        }
    }
}

impl App {
    fn draw(&mut self) {
        let Some(gpu) = self.gpu.as_ref() else {
            return;
        };
        let aspect = gpu.config.width as f64 / gpu.config.height.max(1) as f64;
        let (tris, verts) = self.vertices(aspect);

        let Ok(surface_texture) = gpu.surface.get_current_texture() else {
            return;
        };
        let view = surface_texture
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());
        let solid_buffer = (!tris.is_empty()).then(|| {
            gpu.device
                .create_buffer_init_lite(bytemuck_lite::cast_slice(&tris))
        });
        let buffer = if verts.is_empty() {
            None
        } else {
            Some(
                gpu.device
                    .create_buffer_init_lite(bytemuck_lite::cast_slice(&verts)),
            )
        };

        let mut encoder = gpu
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("lines"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(BACKGROUND),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &gpu.depth,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            // The solid first and the lines over it, both depth-tested, so a marker inside a
            // block is hidden by it rather than drawn through it.
            if let Some(b) = &solid_buffer {
                pass.set_pipeline(&gpu.solid);
                pass.set_vertex_buffer(0, b.slice(..));
                pass.draw(0..tris.len() as u32, 0..1);
            }
            if let Some(b) = &buffer {
                pass.set_pipeline(&gpu.pipeline);
                pass.set_vertex_buffer(0, b.slice(..));
                pass.draw(0..verts.len() as u32, 0..1);
            }
        }
        gpu.queue.submit(Some(encoder.finish()));
        surface_texture.present();

        if self.playing && !self.run.frames.is_empty() {
            self.frame = (self.frame + 1) % self.run.frames.len();
        }
    }
}

impl App {
    /// Render one frame with no window at all, and hand back the pixels as RGBA.
    ///
    /// The same camera, the same `viewer-core` segments and the same pipeline the window uses —
    /// only the target differs. A snapshot path that built its own vertices would be checking
    /// itself.
    fn snapshot(&self, width: u32, height: u32) -> Result<Vec<u8>, String> {
        let instance = wgpu::Instance::default();
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        }))
        .ok_or("no adapter")?;
        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("snapshot"),
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::downlevel_defaults(),
                memory_hints: wgpu::MemoryHints::default(),
            },
            None,
        ))
        .map_err(|e| e.to_string())?;

        let format = wgpu::TextureFormat::Rgba8UnormSrgb;
        let target = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("snapshot"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let view = target.create_view(&wgpu::TextureViewDescriptor::default());
        let solid_pipeline = pipeline_for(&device, format, wgpu::PrimitiveTopology::TriangleList);
        let pipeline = pipeline_for(&device, format, wgpu::PrimitiveTopology::LineList);
        let depth = depth_texture(&device, width, height);

        let (tris, verts) = self.vertices(width as f64 / height as f64);
        let solid_buffer = (!tris.is_empty())
            .then(|| device.create_buffer_init_lite(bytemuck_lite::cast_slice(&tris)));
        let buffer = (!verts.is_empty())
            .then(|| device.create_buffer_init_lite(bytemuck_lite::cast_slice(&verts)));

        // Copies out of a texture want rows padded to 256 bytes, so the readback is wider than the
        // image and the padding is dropped below rather than left in the file as a smear.
        let unpadded = width * 4;
        let padded = unpadded.div_ceil(256) * 256;
        let readback = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback"),
            size: (padded * height) as u64,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        let mut encoder =
            device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("snapshot"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(BACKGROUND),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &depth,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            // The same order the window draws in, so the picture a test reads is the picture a
            // person sees.
            if let Some(b) = &solid_buffer {
                pass.set_pipeline(&solid_pipeline);
                pass.set_vertex_buffer(0, b.slice(..));
                pass.draw(0..tris.len() as u32, 0..1);
            }
            if let Some(b) = &buffer {
                pass.set_pipeline(&pipeline);
                pass.set_vertex_buffer(0, b.slice(..));
                pass.draw(0..verts.len() as u32, 0..1);
            }
        }
        encoder.copy_texture_to_buffer(
            wgpu::ImageCopyTexture {
                texture: &target,
                mip_level: 0,
                origin: wgpu::Origin3d::ZERO,
                aspect: wgpu::TextureAspect::All,
            },
            wgpu::ImageCopyBuffer {
                buffer: &readback,
                layout: wgpu::ImageDataLayout {
                    offset: 0,
                    bytes_per_row: Some(padded),
                    rows_per_image: Some(height),
                },
            },
            wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
        );
        queue.submit(Some(encoder.finish()));

        let slice = readback.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        device.poll(wgpu::Maintain::Wait);
        rx.recv()
            .map_err(|e| e.to_string())?
            .map_err(|e| e.to_string())?;

        let mapped = slice.get_mapped_range();
        let mut out = Vec::with_capacity((unpadded * height) as usize);
        for row in 0..height {
            let from = (row * padded) as usize;
            out.extend_from_slice(&mapped[from..from + unpadded as usize]);
        }
        drop(mapped);
        readback.unmap();
        Ok(out)
    }
}

/// A plain PPM, because a picture nothing can open is not evidence and an encoder is a dependency.
/// A thumbnail's size. Wide enough to tell a room's standing wave from a block's hot core in a
/// grid of thirty, small enough that thirty of them are under two hundred kilobytes.
pub const THUMB: (u32, u32) = (240, 156);

/// Trim the background around what was drawn, then pad back to `THUMB`'s aspect.
///
/// **The camera fits the panel's bounds, and a bar is a line across an otherwise empty frame.**
/// At thumbnail size that is a picture of nothing: measured on `05-beam-on-bar`, where the drawn
/// content was a diagonal hairline in a frame that was 95% background. Cropping to the content is
/// what makes a grid of these worth looking at.
///
/// The background is the corner pixel, which is what the snapshot's own "not background" count
/// already assumes: the shell follows no theme, but the clear colour is one place and this reads
/// it rather than naming it twice.
///
/// # The magnification is clamped, and the number is chosen rather than derived
///
/// A scene with one body has a content box one marker wide, and cropping to that fills the tile
/// with a single cross — a picture of the *marker* rather than of the run. Measured on
/// `07-bouncing-ball`, which is exactly one body. **Three** is a rendering taste; it is written
/// down as one so the next reader does not go looking for where it came from.
fn crop_to_content(w: u32, h: u32, rgba: &[u8]) -> (u32, u32, Vec<u8>) {
    const MAGNIFY: f64 = 3.0;
    const MARGIN: f64 = 0.06;
    let at = |x: u32, y: u32| ((y * w + x) * 4) as usize;
    let bg = [rgba[0], rgba[1], rgba[2]];
    let (mut x0, mut y0, mut x1, mut y1) = (w as i64, h as i64, -1i64, -1i64);
    for y in 0..h {
        for x in 0..w {
            let k = at(x, y);
            if (0..3).any(|c| rgba[k + c].abs_diff(bg[c]) > 6) {
                x0 = x0.min(x as i64);
                y0 = y0.min(y as i64);
                x1 = x1.max(x as i64);
                y1 = y1.max(y as i64);
            }
        }
    }
    if x1 < x0 || y1 < y0 {
        return (w, h, rgba.to_vec());
    }

    let (mut cw, mut ch) = ((x1 - x0 + 1) as f64, (y1 - y0 + 1) as f64);
    let (cx, cy) = (x0 as f64 + cw / 2.0, y0 as f64 + ch / 2.0);
    let pad = cw.max(ch) * MARGIN;
    cw += 2.0 * pad;
    ch += 2.0 * pad;
    cw = cw.max(w as f64 / MAGNIFY);
    ch = ch.max(h as f64 / MAGNIFY);

    // To the tile's aspect, so nothing is stretched.
    let want = THUMB.0 as f64 / THUMB.1 as f64;
    if cw / ch < want {
        cw = ch * want;
    } else {
        ch = cw / want;
    }
    let (x0, y0) = (cx - cw / 2.0, cy - ch / 2.0);

    let (ow, oh) = (cw.round() as u32, ch.round() as u32);
    let mut out = vec![0u8; (ow * oh * 4) as usize];
    for y in 0..oh {
        for x in 0..ow {
            let (sx, sy) = (x0 as i64 + x as i64, y0 as i64 + y as i64);
            let k = ((y * ow + x) * 4) as usize;
            if sx >= 0 && sy >= 0 && (sx as u32) < w && (sy as u32) < h {
                let j = at(sx as u32, sy as u32);
                out[k..k + 4].copy_from_slice(&rgba[j..j + 4]);
            } else {
                out[k..k + 3].copy_from_slice(&bg);
                out[k + 3] = 255;
            }
        }
    }
    (ow, oh, out)
}

/// Box-average down to `THUMB`.
///
/// Averaged rather than sampled, because a marker one pixel wide vanishes under
/// nearest-neighbour: five bodies at 1100 points across become nothing at 240.
fn shrink(w: u32, h: u32, rgba: &[u8], tw: u32, th: u32) -> Vec<u8> {
    let mut out = vec![255u8; (tw * th * 4) as usize];
    for y in 0..th {
        let y0 = y * h / th;
        let y1 = (((y + 1) * h / th).max(y0 + 1)).min(h);
        for x in 0..tw {
            let x0 = x * w / tw;
            let x1 = (((x + 1) * w / tw).max(x0 + 1)).min(w);
            let (mut r, mut g, mut b, mut n) = (0u64, 0u64, 0u64, 0u64);
            for yy in y0..y1 {
                for xx in x0..x1 {
                    let k = ((yy * w + xx) * 4) as usize;
                    r += rgba[k] as u64;
                    g += rgba[k + 1] as u64;
                    b += rgba[k + 2] as u64;
                    n += 1;
                }
            }
            let k = ((y * tw + x) * 4) as usize;
            out[k] = (r / n) as u8;
            out[k + 1] = (g / n) as u8;
            out[k + 2] = (b / n) as u8;
        }
    }
    out
}

/// Write the pixels, as whatever the path's extension asks for.
///
/// PPM is what this wrote and is a format you can produce in six lines and read in a hex editor,
/// which is why the snapshot started there. PNG is what a picture that goes anywhere else has to
/// be, and `image` is already in this binary's tree through `eframe` — declaring it added no crate
/// to the lockfile, which is the only reason it is here.
/// How long each frame of a written GIF is shown, in hundredths of a second.
///
/// Six, so a forty-eight frame run loops in 2.9 s. The GIF format counts in centiseconds and
/// nothing rounder is close: five is 2.4 s and reads as a twitch, ten is 4.8 s and reads as a
/// slideshow.
const GIF_CENTISECONDS: u32 = 6;

/// Every frame as one looping GIF.
///
/// **The encoder quantises each frame on its own**, which is the format's doing and not a choice
/// here: a GIF frame carries at most 256 colours. On a run like this one — a flat background and
/// a smoothly shaded solid — that is invisible at the sizes these are looked at, and the
/// alternative is a palette built across the whole run, which is a colour quantiser this
/// workspace would then own and test.
fn write_gif(path: &str, width: u32, height: u32, frames: &[Vec<u8>], centiseconds: u32) {
    let file = match std::fs::File::create(path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("  cannot write {path}: {e}");
            std::process::exit(1);
        }
    };
    let mut encoder = image::codecs::gif::GifEncoder::new(std::io::BufWriter::new(file));
    if let Err(e) = encoder.set_repeat(image::codecs::gif::Repeat::Infinite) {
        eprintln!("  cannot write {path}: {e}");
        std::process::exit(1);
    }
    for (i, rgba) in frames.iter().enumerate() {
        let Some(buffer) = image::RgbaImage::from_raw(width, height, rgba.clone()) else {
            eprintln!("  frame {i} is not {width}x{height}");
            std::process::exit(1);
        };
        let delay = image::Delay::from_numer_denom_ms(centiseconds * 10, 1);
        if let Err(e) = encoder.encode_frame(image::Frame::from_parts(buffer, 0, 0, delay)) {
            eprintln!("  cannot write frame {i} of {path}: {e}");
            std::process::exit(1);
        }
    }
    drop(encoder);
    let bytes = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    println!(
        "  wrote {path} — {} frames at {centiseconds} centiseconds, {:.2} MB",
        frames.len(),
        bytes as f64 / 1_048_576.0
    );
}

fn write_image(path: &str, width: u32, height: u32, rgba: &[u8]) {
    let png = std::path::Path::new(path)
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("png"));
    let wrote = if png {
        // **Three channels, not four.** The pass writes an opaque frame and the alpha is 255
        // everywhere; carrying it cost 536 kB across thirty tiles against 188 kB without, for a
        // channel that says nothing.
        let (pixels, _) = rgba.as_chunks::<4>();
        let rgb: Vec<u8> = pixels.iter().flat_map(|p| p[..3].to_vec()).collect();
        // **Best compression, chosen rather than defaulted.** These are written once and read
        // from a binary that embeds them, so encoder time is free and every kilobyte is carried
        // in the repository for as long as the scene exists: the default setting cost 520 kB
        // across twenty-seven tiles against 188 kB at this one.
        std::fs::File::create(path)
            .map_err(|e| e.to_string())
            .and_then(|f| {
                image::codecs::png::PngEncoder::new_with_quality(
                    std::io::BufWriter::new(f),
                    image::codecs::png::CompressionType::Best,
                    image::codecs::png::FilterType::Adaptive,
                )
                .write_image(&rgb, width, height, image::ExtendedColorType::Rgb8)
                .map_err(|e| e.to_string())
            })
    } else {
        let mut out = format!("P6\n{width} {height}\n255\n").into_bytes();
        let (pixels, _) = rgba.as_chunks::<4>();
        for p in pixels {
            out.extend_from_slice(&p[..3]);
        }
        std::fs::write(path, out).map_err(|e| e.to_string())
    };
    if let Err(e) = wrote {
        eprintln!("cannot write {path}: {e}");
        std::process::exit(1);
    }
}

/// The clear colour, shared by the window and the snapshot so the two agree.
const BACKGROUND: wgpu::Color = wgpu::Color {
    r: 0.039,
    g: 0.055,
    b: 0.075,
    a: 1.0,
};

/// What the depth buffer is, in both paths.
const DEPTH: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;

/// A depth texture the size of what is being drawn into.
fn depth_texture(device: &wgpu::Device, width: u32, height: u32) -> wgpu::TextureView {
    device
        .create_texture(&wgpu::TextureDescriptor {
            label: Some("depth"),
            size: wgpu::Extent3d {
                width: width.max(1),
                height: height.max(1),
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: DEPTH,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        })
        .create_view(&wgpu::TextureViewDescriptor::default())
}

/// The pipeline, built once here so the window and the snapshot cannot drift apart.
///
/// **Two topologies, one shader.** A field is a solid and is drawn as one — triangles from
/// `pantometry_view::mesh`, the same geometry the glTF and USD exporters write and the editor's
/// viewport shades — and the lines are what is left: a `Paths` panel, a body's cross, the box.
/// Both are depth-tested, which is what a solid needs and what sorting on the CPU could only
/// approximate.
fn pipeline_for(
    device: &wgpu::Device,
    format: wgpu::TextureFormat,
    topology: wgpu::PrimitiveTopology,
) -> wgpu::RenderPipeline {
    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("lines"),
        source: wgpu::ShaderSource::Wgsl(SHADER.into()),
    });
    let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: None,
        bind_group_layouts: &[],
        push_constant_ranges: &[],
    });
    device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("lines"),
        layout: Some(&layout),
        vertex: wgpu::VertexState {
            module: &shader,
            entry_point: "vs",
            buffers: &[wgpu::VertexBufferLayout {
                array_stride: std::mem::size_of::<Vertex>() as u64,
                step_mode: wgpu::VertexStepMode::Vertex,
                attributes: &[
                    wgpu::VertexAttribute {
                        offset: 0,
                        shader_location: 0,
                        format: wgpu::VertexFormat::Float32x3,
                    },
                    wgpu::VertexAttribute {
                        offset: 12,
                        shader_location: 1,
                        format: wgpu::VertexFormat::Float32x3,
                    },
                ],
            }],
            compilation_options: Default::default(),
        },
        fragment: Some(wgpu::FragmentState {
            module: &shader,
            entry_point: "fs",
            targets: &[Some(format.into())],
            compilation_options: Default::default(),
        }),
        primitive: wgpu::PrimitiveState {
            topology,
            // **No culling.** A field's boundary is closed, so back faces are hidden by the depth
            // test anyway; an isosurface's is not, and culling it would put holes in a level set
            // depending on which way the camera happened to be.
            cull_mode: None,
            ..Default::default()
        },
        depth_stencil: Some(wgpu::DepthStencilState {
            format: DEPTH,
            depth_write_enabled: true,
            depth_compare: wgpu::CompareFunction::LessEqual,
            stencil: wgpu::StencilState::default(),
            bias: wgpu::DepthBiasState::default(),
        }),
        multisample: wgpu::MultisampleState::default(),
        multiview: None,
        cache: None,
    })
}

/// A vertex buffer from bytes, without pulling in `wgpu::util`.
trait BufferInit {
    fn create_buffer_init_lite(&self, bytes: &[u8]) -> wgpu::Buffer;
}

impl BufferInit for wgpu::Device {
    fn create_buffer_init_lite(&self, bytes: &[u8]) -> wgpu::Buffer {
        let buffer = self.create_buffer(&wgpu::BufferDescriptor {
            label: Some("lines"),
            size: bytes.len() as u64,
            usage: wgpu::BufferUsages::VERTEX,
            mapped_at_creation: true,
        });
        buffer
            .slice(..)
            .get_mapped_range_mut()
            .copy_from_slice(bytes);
        buffer.unmap();
        buffer
    }
}

/// The same ramp the HTML report uses, so the two views of one run agree about colour.
fn ramp(t: f64) -> [f32; 3] {
    let t = t.clamp(0.0, 1.0);
    let stops = [
        [24.0, 52.0, 110.0],
        [61.0, 139.0, 255.0],
        [255.0, 194.0, 71.0],
        [255.0, 82.0, 33.0],
    ];
    let (a, b, u) = if t < 0.34 {
        (stops[0], stops[1], t / 0.34)
    } else if t < 0.67 {
        (stops[1], stops[2], (t - 0.34) / 0.33)
    } else {
        (stops[2], stops[3], (t - 0.67) / 0.33)
    };
    [
        ((a[0] + (b[0] - a[0]) * u) / 255.0) as f32,
        ((a[1] + (b[1] - a[1]) * u) / 255.0) as f32,
        ((a[2] + (b[2] - a[2]) * u) / 255.0) as f32,
    ]
}

const SHADER: &str = r#"
struct In { @location(0) pos: vec3<f32>, @location(1) colour: vec3<f32> };
struct Out { @builtin(position) clip: vec4<f32>, @location(0) colour: vec3<f32> };

@vertex
fn vs(v: In) -> Out {
    var out: Out;
    out.clip = vec4<f32>(v.pos, 1.0);
    out.colour = v.colour;
    return out;
}

@fragment
fn fs(v: Out) -> @location(0) vec4<f32> {
    return vec4<f32>(v.colour, 1.0);
}
"#;

#[cfg(test)]
mod tests {
    use super::{ends_that_fit, numbered, scale_bar, scale_label};

    /// **A field is coloured on its own scale whichever panel is selected.**
    ///
    /// `panel_vertices` worked out each panel's own span — its comment said so — and then handed
    /// the field branch `self.span`, the selected panel's. Scene 25 has a temperature in kelvin and
    /// a displacement in metres on one grid: the displacement, some 1e-6 m, was painted on a scale
    /// running from 300 K, every value fell off the bottom of it, and the tile showed a solid
    /// block of the coldest colour where it had shown the temperature. A per-panel colour is the
    /// same colour whichever panel happens to own the legend, so the two renderings must agree.
    #[test]
    fn a_field_takes_its_own_scale_and_not_the_selected_panels() {
        let field = |name: &str, unit: &str, values: [f64; 8]| {
            format!(
                r#"{{"name": "{name}", "unit": "{unit}", "kind": "field", "nx": 2, "ny": 2,
                    "nz": 2, "extent_m": [0, 0, 0, 0.01, 0.01, 0.01], "values": {values:?}}}"#
            )
        };
        let hot = field("hot", "K", [300., 310., 320., 330., 340., 350., 360., 400.]);
        let small = field(
            "small",
            "m",
            [1e-6, 1.1e-6, 1.2e-6, 1.3e-6, 1.4e-6, 1.5e-6, 1.6e-6, 2e-6],
        );
        let json = format!(
            r#"{{"format": 3, "title": "two scales", "frames": [{{"t": 0, "panels": [{hot}, {small}],
                "readings": []}}]}}"#
        );
        let colours = |selected: &str| {
            let run = viewer_core::Run::from_json(&json).expect("the run parses");
            let app = super::App::new(run, selected.to_string());
            let panel = &app.run.frames[0].panels[1];
            let (tris, _) = app.panel_vertices(panel, 16.0 / 9.0);
            assert!(!tris.is_empty(), "a 2x2x2 field draws a shell");
            tris.iter().map(|v| v.colour).collect::<Vec<_>>()
        };
        let (owned, alone) = (colours("hot"), colours("small"));
        let differ = owned.iter().zip(&alone).filter(|(a, b)| a != b).count();
        assert!(
            owned.len() == alone.len() && differ == 0,
            "the metres field changed colour when the kelvin field owned the legend: {differ} of \
             {} vertices differ",
            alone.len()
        );
    }

    /// **A colour bar's two ends do not run through each other, at any width they need.**
    ///
    /// They were drawn at the legend'''s own size whatever they said, so `3.9674E-11` and
    /// `9.5054E-11` at the ends of a 0.60-wide bar -- which want 0.79 between them -- came out as
    /// `3.9674E-5054E-11`. Two numbers through each other read as one, and nothing about the
    /// picture says otherwise. The editor'''s colour bar had the same defect and was fixed first;
    /// this is the second implementation, and finding it only after fixing the first is the
    /// argument for this test rather than for that one.
    #[test]
    fn a_colour_bars_two_ends_always_fit_between_them() {
        let (em, aspect, bar) = (0.011f32, 1100.0 / 720.0, 0.60f32);
        // Every pair the formatter can produce, from a short integer to a long exponent.
        for ends in [4.0f32, 12.0, 40.0, 110.0, 220.0, 400.0] {
            let fits = ends_that_fit(ends, bar, em, aspect);
            let taken = (ends + crate::glyphs::width("0")) * fits / aspect;
            assert!(
                taken <= bar + 8.0 * f32::EPSILON * bar,
                "labels {ends} lattice units wide take {taken} of a {bar} bar at {fits}"
            );
            assert!(fits > 0.0, "a label drawn at {fits} is not drawn");
            assert!(
                fits <= em,
                "{fits} is larger than the legend'''s own {em}: a short pair should not grow"
            );
        }
        // A pair that already fits is left alone, so the numbers on most bars are the size every
        // other label on the picture is.
        assert_eq!(ends_that_fit(12.0, bar, em, aspect), em);
    }

    /// The scale bar's metres per half-width, against the drawn geometry's, on two points placed
    /// symmetrically about the framing's centre a known distance apart.
    ///
    /// **The geometry is what reached the GPU**: the two crosses' centres and depths as `f32`
    /// vertices out of `panel_vertices`, not a call to the function the bar uses. A segment of
    /// length `l` through the centre has its ends at depths `d_a` and `d_b`, and the part of it
    /// square to the view is `h = sqrt(l^2 - ((d_b - d_a) * span)^2)`. Each end lands
    /// `(h / 2) * k * D / d` from the centre's image, where `k` is the bar's metres-to-screen at
    /// the centre's depth `D` -- which is the camera's distance, by the projection's definition --
    /// so the drawn length is `(h / 2) * k * D * (1 / d_a + 1 / d_b)`, exactly, and solving for
    /// `k` gives what the ruler must say. No linearisation: perspective is in it, and so is the
    /// foreshortening of a segment that is not square to the view.
    ///
    /// Returns `(bar, geometry, tolerance)`, each a relative quantity.
    fn bar_against_geometry(app: &super::App, name: &str, length_m: f64) -> (f64, f64, f64) {
        let lines = lines_of(app, name);
        assert!(
            lines.len() >= 8,
            "{name} drew {} line vertices",
            lines.len()
        );
        let (first, last) = (&lines[..4], &lines[lines.len() - 4..]);
        let centre = |c: &[[f32; 3]]| ((c[0][0] + c[1][0]) / 2.0, (c[0][1] + c[1][1]) / 2.0);
        let ((ax, ay), (bx, by)) = (centre(first), centre(last));
        let (da, db) = (first[0][2] as f64, last[0][2] as f64);
        // In half-widths both ways: device `y` is a half-height, and the projection divides `x`
        // by the aspect so that a pixel is a pixel in both directions.
        let dx = (bx - ax) as f64;
        let dy = (by - ay) as f64 / ASPECT;
        let drawn = (dx * dx + dy * dy).sqrt();
        let along = (db - da) * app.framing.span;
        let square = (length_m * length_m - along * along).sqrt();
        let depth = app.camera.distance;
        let geometry = drawn / (square / 2.0 * depth * (1.0 / da + 1.0 / db));

        let (nice, across) = app.ruler(ASPECT).expect("a framed run has a ruler");
        let bar = across as f64 / nice;
        // **Earned from `f32`.** Every vertex coordinate and depth, and the bar's own width, is an
        // `f32` of order one, so each carries at most half an `f32` epsilon relative to unity; the
        // two crosses' centres are an average of two such numbers. The drawn length is a
        // difference of coordinates no larger than 1, so its relative error is that epsilon times
        // `1 / drawn`. Eight epsilons against unity in all, divided by the shorter of the two
        // lengths measured, bounds it -- a few parts in a million here, where the old identity
        // was out by a factor of 3.7.
        let tolerance = 8.0 * f32::EPSILON as f64 / drawn.min(across as f64);
        (bar, geometry, tolerance)
    }

    /// A box of known size, 30 x 10 x 5 mm well away from the origin, marked by a body at each of
    /// two opposite corners -- its space diagonal, which runs through the framing's centre.
    fn a_box() -> super::App {
        let json = r#"{"format": 3, "title": "a box", "frames": [{"t": 0, "panels": [
            {"name": "corners", "unit": "K", "kind": "points", "boxed": false,
             "bounds": [0.2, -0.1, 0.05, 0.23, -0.09, 0.055],
             "positions": [0.2, -0.1, 0.05, 0.23, -0.09, 0.055],
             "values": [1, 2]}],
            "readings": []}]}"#;
        let run = viewer_core::Run::from_json(json).expect("the run parses");
        let mut app = super::App::new(run, "corners".to_string());
        app.legend = false;
        app
    }

    /// **The bar covers, on the screen, the metres it names: measured against the geometry.**
    ///
    /// The test that was here checked `across * widest * 2 == nice` -- "half the screen is half
    /// the subject" -- which is the identity the bar was computed from, so it could not disagree
    /// with it. The identity was false: `Camera::fit` sets the focal length so the subject lands
    /// wherever its furthest corner says, and framed on scene 14's 20 mm bar the `10 MM` ruler
    /// was 137 px under a bar 1018 px long.
    ///
    /// So this draws things of a known size and measures them: a box's diagonal from the opening
    /// camera, turned and zoomed in and out -- a ruler that does not follow a zoom is the wrong
    /// one after the first scroll -- and scene 14's bar framed on itself.
    #[test]
    fn the_bar_measures_what_the_geometry_measures() {
        let diagonal = (0.03f64 * 0.03 + 0.01 * 0.01 + 0.005 * 0.005).sqrt();
        let mut cases: Vec<(String, super::App, &str, f64)> = Vec::new();
        cases.push(("a box, as it opens".into(), a_box(), "corners", diagonal));
        let mut turned = a_box();
        turned.camera.turn(1.1, -0.5);
        cases.push(("a box, turned".into(), turned, "corners", diagonal));
        let mut near = a_box();
        near.camera.zoom(0.6);
        cases.push(("a box, zoomed in".into(), near, "corners", diagonal));
        let mut far = a_box();
        far.camera.zoom(2.9);
        cases.push(("a box, zoomed out".into(), far, "corners", diagonal));
        // Scene 14's problem: the run's bar, framed on itself. The 61 samples are cell-centred,
        // so the first and last sit half a cell in from each end and symmetric about the centre.
        let mut bar = two_scales();
        bar.frame_on(Some("bar".to_string()));
        cases.push((
            "scene 14's bar, framed".into(),
            bar,
            "bar",
            0.02 * 60.0 / 61.0,
        ));

        for (what, app, name, length) in &cases {
            let (bar, geometry, tolerance) = bar_against_geometry(app, name, *length);
            let off = (bar - geometry).abs() / geometry;
            let (nice, across) = app.ruler(ASPECT).expect("a ruler");
            println!(
                "  {what}: {} over {across:.4} of the half-width; the bar says {bar:.6e} per \
                 metre, the geometry {geometry:.6e} ({off:.1e} apart, {tolerance:.1e} allowed)",
                scale_label(nice)
            );
            assert!(
                off <= tolerance,
                "{what}: the scale bar draws {bar:e} of the half-width per metre and the geometry \
                 is drawn at {geometry:e} -- {:.3}x -- so the bar does not cover the {} it names",
                bar / geometry,
                scale_label(nice)
            );
        }
    }

    /// **The number on the bar is round, and the bar is a readable length**, at every scale from
    /// a molecule to a solar system.
    ///
    /// Of `scale_bar` alone, so given how far a metre reaches rather than measuring it: one, two
    /// or five times a decade, and between a twelfth and five twelfths of the width.
    #[test]
    fn the_bar_is_round_and_readable_at_every_scale() {
        for exponent in -3..=13 {
            for mantissa in [1.0, 1.7, 3.0, 4.9, 7.3, 9.9] {
                let per_metre = mantissa * 10f64.powi(exponent);
                let (nice, across) = scale_bar(per_metre).expect("a positive scale has a bar");
                let mantissa_of = nice / 10f64.powf(nice.log10().floor());
                assert!(
                    [1.0, 2.0, 5.0, 10.0]
                        .iter()
                        .any(|m| (m - mantissa_of).abs() < 1e-9),
                    "{nice:e} is not a round number"
                );
                // The length is the label times the scale, and nothing else.
                assert!(
                    ((across as f64) - nice * per_metre).abs() <= f32::EPSILON as f64,
                    "{across} is not {nice:e} m at {per_metre:e} per metre"
                );
                // A bar the width of the screen measures nothing and one a pixel long measures
                // nothing either.
                assert!(
                    (0.05..=0.60).contains(&across),
                    "{per_metre:e} per metre gives a bar {across} of the half-width"
                );
            }
        }
        assert!(
            scale_bar(0.0).is_none(),
            "nothing to measure is not a bar of zero"
        );
        assert!(scale_bar(f64::NAN).is_none());
        assert!(scale_bar(f64::INFINITY).is_none());
        assert!(scale_bar(-1.0).is_none());
    }

    /// **A sequence numbers before the extension, and the dot it finds is the extension's.**
    ///
    /// `../docs/out.png` carries dots in its `..` too. Splitting on the last dot anywhere would
    /// have written `.` beside `./docs/out-007` -- a file in the wrong directory under a name
    /// nothing would go looking for, and nothing about the run would have said so.
    #[test]
    fn a_frame_is_numbered_before_its_extension() {
        assert_eq!(numbered("out.png", 7), "out-007.png");
        assert_eq!(numbered("../docs/out.png", 7), "../docs/out-007.png");
        assert_eq!(numbered(r"..\docs\out.png", 7), r"..\docs\out-007.png");
        // No extension at all: the number goes on the end rather than inventing one.
        assert_eq!(numbered("../out", 7), "../out-007");
        assert_eq!(numbered("frames/f", 0), "frames/f-000");
        // And the order is the frame order, which is the whole reason for three digits.
        let mut names: Vec<String> = [9usize, 10, 100, 1]
            .iter()
            .map(|f| numbered("a.png", *f))
            .collect();
        let sorted = {
            let mut c = names.clone();
            c.sort();
            c
        };
        names.sort_by_key(|n| n.clone());
        assert_eq!(names, sorted);
        assert_eq!(sorted[0], "a-001.png");
    }

    /// **The scale bar names a length at every size this workspace models.**
    ///
    /// It had metres and millimetres and three decimal places, so everything below a micrometre
    /// printed `0 MM` — which is what the protein figure showed, 4.5 nm across, while it was
    /// being put on the front page. A legend reporting zero length for a real object is worse
    /// than no legend.
    #[test]
    fn the_scale_bar_names_a_length_at_every_size() {
        // One, two and five times every decade from a picometre to a kilometre. The bar is
        // snapped to that ladder before it is labelled, so these are the only values it is ever
        // built from — an atom at one end and a room at the other.
        for exponent in -12..=3 {
            for mantissa in [1.0, 2.0, 5.0] {
                let metres = mantissa * 10f64.powi(exponent);
                let label = scale_label(metres);
                let (number, unit) = label.split_once(' ').expect("a number and a unit");
                let value: f64 = number.parse().expect("a number");
                assert!(
                    value > 0.0,
                    "{metres:e} m is labelled {label:?}, which is no length at all"
                );
                // **And the number means what its unit says.** A prefix off by a thousand still
                // prints a positive number, and that is the failure a reader cannot catch: the
                // bar drawn one size and labelled another. Read back through its own unit, it has
                // to give the metres it was made from.
                let per_metre = match unit {
                    "KM" => 1e-3,
                    "M" => 1.0,
                    "MM" => 1e3,
                    "UM" => 1e6,
                    "NM" => 1e9,
                    "PM" => 1e12,
                    other => panic!("{other:?} is not a unit this knows"),
                };
                let back = value / per_metre;
                // The label carries three decimals of its own unit, so that is its resolution and
                // the whole of what the round trip can lose.
                assert!(
                    (back - metres).abs() <= 5e-4 / per_metre,
                    "{metres:e} m is labelled {label:?}, which reads back as {back:e}"
                );
            }
        }
    }

    /// **Every unit the bar can print is one the glyph table can draw.** `UM` and `NM` were added
    /// for this, and a unit whose letters render as gaps would read as `M` — metres — which is a
    /// scale bar that lies by a factor of a million.
    #[test]
    fn every_unit_the_bar_prints_can_be_drawn() {
        for exponent in -12..=3 {
            let label = scale_label(10f64.powi(exponent));
            assert!(
                crate::glyphs::can_draw(&label),
                "{label:?} has a character the glyph table draws as a box"
            );
        }
    }

    /// The width over height `App::new` fits at.
    const ASPECT: f64 = 16.0 / 9.0;

    /// The smallest run with scene 14's problem: a 20 mm bar of 61 samples, and two planets
    /// 1.5e11 m away. Framed whole, the bar is `1e-13` of the framing.
    fn two_scales() -> super::App {
        let values: Vec<String> = (0..61)
            .map(|i| format!("{}", 293.0 + i as f64 / 60.0))
            .collect();
        let json = format!(
            r#"{{"format": 3, "title": "two scales", "frames": [{{"t": 0, "panels": [
                {{"name": "bar", "unit": "K", "kind": "field", "nx": 61, "ny": 1, "nz": 1,
                  "extent_m": [0, 0, 0, 0.02, 0, 0], "lattice": "centred",
                  "values": [{}]}},
                {{"name": "sky", "unit": "m/s", "kind": "points", "boxed": false,
                  "bounds": [-1.0e11, -1.0e11, -1.0e10, 1.5e11, 1.0e11, 1.0e10],
                  "positions": [1.5e11, 0, 0, -1.0e11, 1.0e11, 1.0e10],
                  "values": [30000, 25000]}}],
                "readings": []}}]}}"#,
            values.join(", ")
        );
        let run = viewer_core::Run::from_json(&json).expect("the run parses");
        let mut app = super::App::new(run, "bar".to_string());
        app.legend = false;
        app
    }

    /// One panel's line vertices, in normalised device coordinates, as they go to the GPU.
    fn lines_of(app: &super::App, name: &str) -> Vec<[f32; 3]> {
        let panel = app.run.frames[0]
            .panels
            .iter()
            .find(|p| p.name() == name)
            .expect("the panel is in the run");
        let (_, lines) = app.panel_vertices(panel, ASPECT);
        lines.iter().map(|v| v.position).collect()
    }

    /// The centre of each cross, as the GPU receives it: a horizontal arm's two ends, averaged.
    fn crosses(lines: &[[f32; 3]]) -> Vec<(f32, f32)> {
        lines
            .chunks(4)
            .map(|c| ((c[0][0] + c[1][0]) / 2.0, (c[0][1] + c[1][1]) / 2.0))
            .collect()
    }

    /// **Framed on the bar, its 61 samples are 61 places on the screen.**
    ///
    /// They reach the GPU as `f32` screen positions, so this is after the narrowing. Framed whole,
    /// the run's longest side is 2.5e11 m and the bar is `8e-14` of it: every one of its 61
    /// crosses was drawn at **one** position, and the bar was a single `+`. Framed on itself it
    /// spans most of the window: `Camera::fit` puts its outer end at 0.85 of the half-frame, and
    /// the claim is that the 61 centres reach across at least 0.7 of the viewport along one axis.
    #[test]
    fn framing_one_panel_spreads_it_across_the_window() {
        let mut app = two_scales();
        let distinct = |c: &[(f32, f32)]| {
            let mut seen: Vec<(f32, f32)> = Vec::new();
            for p in c {
                if !seen.contains(p) {
                    seen.push(*p);
                }
            }
            seen.len()
        };
        let whole = crosses(&lines_of(&app, "bar"));
        println!(
            "  framed whole, the bar's 61 samples are at {} places",
            distinct(&whole)
        );

        app.frame_on(Some("bar".to_string()));
        let framed = crosses(&lines_of(&app, "bar"));
        assert_eq!(framed.len(), 61, "a sample was dropped");
        assert_eq!(
            distinct(&framed),
            61,
            "framed on the bar, its 61 samples reach the GPU at {} distinct positions",
            distinct(&framed)
        );
        let reach = |axis: usize| {
            let v: Vec<f32> = framed
                .iter()
                .map(|c| if axis == 0 { c.0 } else { c.1 })
                .collect();
            let lo = v.iter().copied().fold(f32::MAX, f32::min);
            let hi = v.iter().copied().fold(f32::MIN, f32::max);
            (hi - lo) / 2.0
        };
        let fraction = reach(0).max(reach(1));
        println!("  framed on the bar, it spans {fraction:.4} of the viewport");
        assert!(
            fraction >= 0.7,
            "framed on the bar, its samples span {fraction:.3e} of the viewport"
        );
        // The scale bar measures what is framed: a ruler for metres-across-the-run under a
        // picture of a 20 mm bar is the wrong ruler by thirteen orders.
        let (nice, _) = app.ruler(ASPECT).expect("framed on the bar, a ruler");
        assert!(
            nice < 0.02,
            "framed on a 20 mm bar the ruler names {nice:e} m: it is not the framed box's"
        );
    }

    /// **What is far from the framed panel is not drawn**, and the whole run brings it back.
    ///
    /// Both planets are 7.5e12 bar-lengths away. A GPU would clip them; this shell projects on
    /// the CPU and hands the pass screen positions, so nothing clipped them and a point behind
    /// the eye — its depth clamped to the near plane by `Camera::project` — arrived at a coordinate
    /// in the 1e14s. They are dropped by `Camera::sees` now. Framing the whole run again puts both
    /// back inside the viewport.
    #[test]
    fn framing_one_panel_drops_what_is_far_and_the_whole_run_restores_it() {
        let mut app = two_scales();
        app.frame_on(Some("bar".to_string()));
        let far = lines_of(&app, "sky");
        assert!(
            far.is_empty(),
            "framed on the bar, the planets 1.5e11 m away were drawn: {far:?}"
        );

        app.frame_on(None);
        assert_eq!(
            app.panel, "bar",
            "the colour bar did not return to the panel it opened on"
        );
        let back = crosses(&lines_of(&app, "sky"));
        assert_eq!(
            back.len(),
            2,
            "framed whole, the planets are not both drawn"
        );
        for (x, y) in back {
            assert!(
                x.abs() <= 1.0 && y.abs() <= 1.0,
                "framed whole, a planet is off-screen at ({x}, {y})"
            );
        }
        let (nice, _) = app.ruler(ASPECT).expect("framed whole, a ruler");
        assert!(
            nice > 1e10,
            "framed on a 2.5e11 m run the ruler names {nice:e} m: it is not the whole run's"
        );
    }

    /// **F walks the whole run, then each panel in order, then the whole run again.**
    #[test]
    fn f_cycles_through_the_whole_run_and_each_panel() {
        let mut app = two_scales();
        let mut seen = vec![app.framed.clone()];
        for _ in 0..3 {
            let next = app.next_framing();
            app.frame_on(next);
            seen.push(app.framed.clone());
        }
        assert_eq!(
            seen,
            [None, Some("bar".into()), Some("sky".into()), None],
            "F did not cycle whole, bar, sky, whole"
        );
    }
}
