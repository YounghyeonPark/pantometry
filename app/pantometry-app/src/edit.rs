//! Edit a pantometry scene beside a 3D view of it — and leave it open while a script does the
//! editing.
//!
//! ```text
//! cargo run --release                 # opens on the built-in room
//! cargo run --release -- scene.json   # opens on a file
//! ```
//!
//! The left pane is the scene's JSON, checked as you type with the same two steps
//! `pantometry check` runs; a parse error is shown with its `line:column`. The viewport
//! draws every placed extent as a wireframe — live, from the text, before anything runs — and
//! **Run** streams the run in as it computes: each frame appears when it is captured, the
//! slider grows, and **stop** ends a long run between frames. **Verify** runs the battery
//! from `pantometry-world verify` and shows the same report the CLI prints.
//!
//! # The live loop
//!
//! With **watch file** on (it is on by default), the editor polls the file's modified time
//! and reloads when something else writes it — a script, an agent, another editor. With
//! **run on change** on as well, the reload runs the scene, so the loop
//! `script writes → editor rechecks → runs → draws` closes with no hand on the window.
//! An in-flight run is stopped and superseded when the file changes again, so the picture
//! converges on the latest text rather than queueing history.
//!
//! One rule keeps that honest: **unsaved edits in the pane are never clobbered.** If the pane
//! is dirty and the disk changes, the status line says so, loudly, and the disk's version
//! waits for an explicit `load`. The alternative silently discards somebody's typing, and
//! which of the two writers meant it is not the editor's call to make.
//!
//! Drag to rotate, scroll to zoom, and the camera is `viewer-core`'s — the same fit, the same
//! projection, the same clamps, because that arithmetic has been wrong here before and is not
//! being written a third time.
//!
//! # The two halves, kept apart
//!
//! ARCHITECTURE.md's platform rules: the authoring half is the composition root and may name
//! domains; the inspection half dispatches on the shape of the data. This file is the shell
//! over both and holds neither — authoring machinery lives in `editor-core`, and everything
//! painted below is painted from a *shape* (a box, points, paths, a reading) with the domain
//! name used only as a label. A domain added next year is drawn by the code below unchanged.

use crate::render;
use editor_core::Beside;
use eframe::egui;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{mpsc, Arc, Mutex};
use std::time::{Duration, Instant, SystemTime};

/// `pantometry edit [scene.json] [--run]`.
pub fn run(args: &[String]) -> i32 {
    match open(args) {
        Ok(()) => 0,
        Err(e) => {
            eprintln!("editor: {e}");
            1
        }
    }
}

fn open(passed: &[String]) -> eframe::Result {
    // `--run` starts the run on load, so `pantometry edit scene.json --run` opens on the picture
    // rather than on three clicks. Written because there was no way to *see* the shaded viewport
    // without a human at the keyboard, which made a rendering change unverifiable.
    //
    // The arguments come from the dispatcher, not from `std::env` — this is a subcommand now, and
    // reading the process's own arguments here would see `edit` as a filename.
    let mut path = None;
    let mut run_at_once = false;
    for arg in passed {
        match arg.as_str() {
            "--run" => run_at_once = true,
            other if other.starts_with('-') => {
                eprintln!("pantometry edit: unknown option {other}");
                eprintln!("usage: pantometry edit [scene.json] [--run]");
                std::process::exit(2);
            }
            other => path = Some(other.to_string()),
        }
    }
    eframe::run_native(
        "pantometry editor",
        eframe::NativeOptions {
            // **Big enough that the viewport exists.** Three side panels have minimum widths, and a
            // `SidePanel` is served before the `CentralPanel` is — so on the default 533-point
            // window the outliner, the text and the inspector took all of it and the paint callback
            // was handed a rect of **zero width**. It drew nothing, correctly, and looked exactly
            // like a renderer that did not work.
            viewport: egui::ViewportBuilder::default().with_inner_size([1500.0, 950.0]),
            ..eframe::NativeOptions::default()
        },
        Box::new(move |_cc| {
            // **No file named means no file open.** It used to mean the built-in scene, which is
            // an answer to a question nobody asked: somebody who typed `pantometry edit` either
            // has a scene somewhere or wants to make one, and both of those are on the start
            // screen. `App::new(None)` still builds the built-in scene underneath, so `New
            // project` has something to hand back and the fields are all populated.
            let mut app = match path {
                Some(p) => App::new(Some(p)),
                None => App::empty(),
            };
            app.run_at_once = run_at_once;
            Ok(Box::new(app))
        }),
    )
}

/// What a background job sent back.
enum Job {
    /// The run so far, as JSON — one of these per captured frame, plus the settled final.
    Frames(String),
    /// How the streaming run ended: finished, stopped, or the violation that refused it.
    RunEnded(Result<editor_core::RunEnd, String>),
    /// The battery's rendered report and its findings count, or why it could not start.
    Verified(Result<(String, usize), String>),
}

/// Where the shaded pass actually puts a run's geometry, in world metres.
///
/// The `--drawn-extent` subcommand's whole body, and the only headless way to ask the editor's
/// viewport what it drew. `batches` is a method on `App`, `App::new` needs no window -- which is
/// what `--layout-at` already relies on -- and `Batch::positions` is framing-local, so multiplying
/// by the span and adding the centre puts it back in metres.
///
/// **Why a number and not a picture.** Six places in `batches` turn a run panel into geometry and
/// each has to apply `Panel::place`; five of six would look completely fine on any scene that
/// states no pose, which is all but one of them. A missed site collapses that panel onto the
/// origin, and the box below is what says so.
///
/// The two batches are reported apart on purpose. The lines carry the *scene's* boxes, which
/// `editor-core` has always placed correctly, so they would span the assembly even with every run
/// panel drawn at the origin -- and a single combined box would hide exactly the defect this
/// exists to catch.
pub fn drawn_extent(scene: Option<String>, run_path: &str) -> Result<String, String> {
    let mut app = App::new(scene);
    let text = std::fs::read_to_string(run_path).map_err(|e| format!("{run_path}: {e}"))?;
    let run = viewer_core::Run::from_json(&text)?;
    app.run = Some(RunView {
        run,
        frame: 0,
        playing: false,
        last_step: None,
        partial: false,
    });

    let world = app.world().unwrap_or([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0]);
    let framing = viewer_core::Framing::of(world);
    let (solid, lines, _, _) = app.batches(&framing);

    let box_of = |b: &render::Batch| {
        if b.positions.is_empty() {
            return None;
        }
        let mut lo = [f64::MAX; 3];
        let mut hi = [f64::MIN; 3];
        for v in b.positions.as_chunks::<3>().0 {
            for a in 0..3 {
                // `Framing::local` subtracts the centre and divides by the span; this is that,
                // backwards.
                let w = v[a] as f64 * framing.span + framing.centre[a];
                lo[a] = lo[a].min(w);
                hi[a] = hi[a].max(w);
            }
        }
        Some([lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]])
    };
    let say = |b: Option<[f64; 6]>| match b {
        Some(v) => format!(
            "[{:.6} {:.6} {:.6} {:.6} {:.6} {:.6}]",
            v[0], v[1], v[2], v[3], v[4], v[5]
        ),
        None => "empty".to_string(),
    };
    Ok(format!(
        "world={} solid={} lines={}",
        say(Some(world)),
        say(box_of(&solid)),
        say(box_of(&lines))
    ))
}

/// One frame of the editor's interface, as text: every string it drew and where.
///
/// The `--ui-dump` subcommand's whole body, in the pattern `--layout-at` set. `App::new` needs no
/// window and `App::ui` needs no `eframe::Frame`, so a frame can be built, laid out and read
/// without a display — and the shaded viewport, which does need a GPU, arrives as a
/// `Shape::Callback` that is recorded and never run.
///
/// # Two passes, and the second is the one reported
///
/// egui settles: a panel's width, a scroll area's extent and a menu's size can depend on what the
/// previous frame measured. One pass reports a layout mid-settle. Two is what the interaction
/// model needs and what a reader sees.
///
/// # Why text and rects rather than pixels
///
/// A screenshot answers "does it look right" only to somebody looking. This answers the questions
/// a test can hold: **is the control there at this width**, **was the label elided**, **is the
/// error where a reader will find it**. Elision is exact rather than inferred — egui ends a
/// truncated galley with `…`, so a label that did not fit says so in its own text.
pub fn ui_dump(
    path: Option<String>,
    width: f32,
    height: f32,
    choosing: Option<Vec<String>>,
    click: Option<(f32, f32)>,
) -> String {
    // The same two doors the window opens by, so what this reports is what a person would see.
    let mut app = match path {
        Some(p) => App::new(Some(p)),
        None => App::empty(),
    };
    // And the third, which a window reaches by clicking New project. Without it the chooser is
    // the one screen this hook cannot see, which is the position the whole editor was in.
    if let Some(ticked) = choosing {
        app.start = true;
        app.choosing = Some(ticked.into_iter().collect());
    }
    let ctx = egui::Context::default();
    let input = || egui::RawInput {
        screen_rect: Some(egui::Rect::from_min_size(
            egui::Pos2::ZERO,
            egui::vec2(width, height),
        )),
        ..Default::default()
    };
    let _ = ctx.run(input(), |c| app.ui(c));
    // **A click, when one is asked for.** A menu's items live inside a closed menu, so a frame
    // built without input can see six menu titles and nothing under them — which made "this
    // command still has a home in a menu" a claim no test could reach. A press and a release at
    // a point is the smallest input that opens one, and egui needs the frame after to lay the
    // opened menu out.
    if let Some((x, y)) = click {
        let at = egui::pos2(x, y);
        let press = |down| egui::Event::PointerButton {
            pos: at,
            button: egui::PointerButton::Primary,
            pressed: down,
            modifiers: egui::Modifiers::default(),
        };
        let mut with_click = input();
        with_click.events = vec![egui::Event::PointerMoved(at), press(true), press(false)];
        let _ = ctx.run(with_click, |c| app.ui(c));
        let _ = ctx.run(input(), |c| app.ui(c));
    }
    let out = ctx.run(input(), |c| app.ui(c));

    let mut texts: Vec<Drawn> = Vec::new();
    let mut callbacks: Vec<egui::Rect> = Vec::new();
    for clipped in &out.shapes {
        collect(
            &clipped.shape,
            clipped.clip_rect,
            &mut texts,
            &mut callbacks,
        );
    }
    // Reading order: down the screen, then across. A dump a person can scan and a test can index.
    texts.sort_by(|a: &Drawn, b| {
        (a.at.top() as i32)
            .cmp(&(b.at.top() as i32))
            .then((a.at.left() as i32).cmp(&(b.at.left() as i32)))
    });

    // **egui's own flag, not the ellipsis at the end of the string.** Read off the suffix, this
    // counted `New project…` and `Open a scene…` as truncated: a trailing `…` is the convention
    // for a control that opens a dialog, and two of them made a start screen report two elided
    // labels with nothing wrong. `Galley::elided` is set by the layout that did the eliding.
    let elided = texts.iter().filter(|d| d.elided).count();
    // **Cut off by the window rather than by a widget that scrolls.** The first form of this asked
    // whether a string reached past the right edge of the window, and the toolbar's path box failed
    // it honestly: a `TextEdit` lays its whole line out and clips it to its own box, so an absolute
    // path in a temp directory is 431 points of galley inside 220 points of field and nothing is
    // wrong. What is wrong is a string clipped by something that goes to the window's own edge --
    // a panel -- which is what `fit view` at x=704 of a 700-point window was.
    // Both edges. Written for the right one only, it stopped seeing anything the moment the
    // toolbar's two longest labels moved to a menu: what runs off a narrow window now is the
    // status bar's right-aligned indicator, pushed off the *left* at 40 points. A string clipped
    // by the window is lost whichever side it went out of.
    let cut = |d: &Drawn| {
        (d.at.right() > d.clip.right() + 0.5 && d.clip.right() >= width - 0.5)
            || (d.at.left() < d.clip.left() - 0.5 && d.clip.left() <= 0.5)
    };
    let cut_off = texts.iter().filter(|d| cut(d)).count();
    // **The viewport's rect, because that rect is the defect this whole hook exists for.** It
    // was handed zero width once and drew nothing, correctly, for an afternoon in the renderer.
    // Reported rather than inferred from what happens to be drawn inside it.
    let mut s = format!(
        "size={width}x{height}\ntexts={}\nelided={elided}\ncut={cut_off}\ncallbacks={}\n",
        texts.len(),
        callbacks.len()
    );
    for r in &callbacks {
        s.push_str(&format!(
            "viewport={:.0},{:.0} {:.0}x{:.0}\n",
            r.left(),
            r.top(),
            r.width(),
            r.height()
        ));
    }
    for d in &texts {
        s.push_str(&format!(
            "{} {:>5.0} {:>5.0} {:>5.0}  {}\n",
            if cut(d) { "!" } else { " " },
            d.at.left(),
            d.at.top(),
            d.at.width(),
            d.text.replace('\n', "\\n")
        ));
    }
    s
}

/// One string the frame drew, where it went, and what was allowed to show of it.
struct Drawn {
    /// Where the string was laid out, in window points.
    at: egui::Rect,
    /// What it was clipped to. A widget that scrolls clips to itself; a panel clips to the window.
    clip: egui::Rect,
    /// The string, with newlines escaped.
    text: String,
    /// egui's own flag, set by the layout that shortened the line.
    elided: bool,
}

/// Walk a shape tree, keeping the text and counting the paint callbacks.
///
/// `Shape::Vec` nests, so this recurses; a flat scan over `FullOutput::shapes` finds the panels'
/// own shapes and none of the widgets inside them.
fn collect(
    shape: &egui::Shape,
    clip: egui::Rect,
    texts: &mut Vec<Drawn>,
    callbacks: &mut Vec<egui::Rect>,
) {
    match shape {
        egui::Shape::Text(t) => {
            let text = t.galley.text().trim().to_string();
            if !text.is_empty() {
                texts.push(Drawn {
                    at: t.galley.rect.translate(t.pos.to_vec2()),
                    clip,
                    text,
                    elided: t.galley.elided,
                });
            }
        }
        egui::Shape::Vec(v) => {
            for s in v {
                collect(s, clip, texts, callbacks);
            }
        }
        egui::Shape::Callback(c) => callbacks.push(c.rect),
        _ => {}
    }
}

/// What the viewport must keep, whatever the panels want.
///
/// 280 points is about the smallest a shaded scene reads at — narrower and the colour bar, the
/// readings and the geometry compete for the same forty pixels. Named here because two places
/// have to agree about it: `panels_that_fit`, which decides how many panels there is room for,
/// and `panel_budget`, which decides how wide each of them may be. They disagreed, and the
/// viewport was zero points wide at four of the eleven widths `--ui-dump` was asked about.
pub const VIEWPORT_FLOOR: f32 = 280.0;

/// The narrowest useful outliner: a tree row with a disclosure arrow, a dot and a truncated name.
pub const OUTLINER_MIN: f32 = 170.0;

/// The narrowest useful scene text: JSON at this indent wraps unreadably below it.
pub const TEXT_MIN: f32 = 240.0;

/// The narrowest useful inspector: a label and a number-drag side by side.
pub const INSPECTOR_MIN: f32 = 200.0;

/// What a side panel costs beyond its own width: the gap egui leaves beside it.
///
/// Read from the style rather than written down — it is egui's number, not this program's, and
/// this program sets no style, so `Style::default()` is the one in force. Left out of both sums
/// at first, and with three panels up the viewport settled at **264** points against a floor of
/// 280: three of these, which is exactly the shortfall. A floor that is missed by the width of
/// the separators is still a floor that does not hold.
fn panel_gap() -> f32 {
    egui::Style::default().spacing.item_spacing.x
}

/// What the panels do at a given window width, as one line.
///
/// The `--layout-at` subcommand's whole body. `panels_that_fit` is a method on `App` because it
/// reads the `show_*` flags; this asks it with all three on, which is the state a reader starts in
/// and the only one where the fitting has anything to decide.
pub fn layout_at(width: f32) -> String {
    let app = App::new(None);
    let (fits, dropped) = app.panels_that_fit(width);
    let taken = (if fits.outliner { 170.0 } else { 0.0 })
        + (if fits.text { 240.0 } else { 0.0 })
        + (if fits.inspector { 200.0 } else { 0.0 });
    format!(
        "width={width} view={} outliner={} text={} inspector={} dropped={}",
        width - taken,
        if fits.outliner { "on" } else { "off" },
        if fits.text { "on" } else { "off" },
        if fits.inspector { "on" } else { "off" },
        dropped.unwrap_or("none"),
    )
}

/// Which side panels are on screen this frame.
///
/// Not the same as the `show_*` flags: those are what the reader asked for, this is what there is
/// room for. Keeping them apart is what makes a widened window restore the panel by itself.
#[derive(Clone, Copy, Debug)]
struct Panels {
    outliner: bool,
    text: bool,
    inspector: bool,
}

/// One thing the cursor can be over in the shaded view, with its readout already written.
///
/// The text is built when the geometry is, because that is when the value is in hand and because a
/// readout assembled per paint is a `format!` per candidate per frame.
struct Probed {
    at: [f64; 3],
    path: String,
    label: String,
}

/// One change the inspector collected this frame, held until the widgets let go of the text.
///
/// A pointer and a value rather than the spliced result: the splice reads `self.text`, and the
/// widgets are still borrowing it while the grid is being drawn. Applying it inside the loop
/// would also invalidate every pointer after the edited one, since a pointer is a position in a
/// string that just changed length.
enum Edited {
    /// A number a `DragValue` moved.
    Number(String, f64),
    /// A string a menu picked or a text field typed.
    Text(String, String),
}

/// The three directions a translate handle points, and the colour each is drawn in.
///
/// x red, y green, z blue: the convention every 3D tool shares, and the one place in this shell
/// where a colour means something other than a value.
const AXES: [([f64; 3], egui::Color32); 3] = [
    ([1.0, 0.0, 0.0], egui::Color32::from_rgb(226, 96, 88)),
    ([0.0, 1.0, 0.0], egui::Color32::from_rgb(120, 196, 116)),
    ([0.0, 0.0, 1.0], egui::Color32::from_rgb(108, 152, 232)),
];

/// How near the pointer has to be to a handle to take hold of it, in points.
const GRAB_RADIUS: f32 = 9.0;

/// What the three axis handles do.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Gizmo {
    /// Arrows along the axes: drag to move the domain.
    Move,
    /// Rings about the axes: drag to turn it.
    Turn,
}

/// How far `at` is from the segment `a`..`b`, in points.
///
/// A handle is a line and a pointer is a point, so the distance to the *segment* is the question
/// — not the distance to either end, which would let a click far past the tip take hold of it.
fn near_segment(at: egui::Pos2, a: egui::Pos2, b: egui::Pos2) -> f32 {
    let ab = b - a;
    let len2 = ab.length_sq();
    if len2 <= f32::EPSILON {
        return (at - a).length();
    }
    let t = (((at - a).dot(ab)) / len2).clamp(0.0, 1.0);
    (at - (a + ab * t)).length()
}

/// A parsed run being scrubbed through.
struct RunView {
    run: viewer_core::Run,
    frame: usize,
    /// Whether the frame advances on its own.
    ///
    /// A run is a sequence, and a slider alone makes a reader drag to see it move. What a reader
    /// is usually here to do is compare two frames, which is a key pressed twice — so there is a
    /// clock, a step either way, and the space bar.
    playing: bool,
    /// When the frame last advanced, so playback runs on a clock and not on the paint rate. At
    /// sixty frames a second a twelve-frame run is over in a fifth of a second.
    last_step: Option<Instant>,
    /// Whether what is on screen is a prefix of a stopped run rather than the run. Set from
    /// [`editor_core::RunEnd`], drawn on the canvas, because a prefix that looks complete is
    /// a picture of something that did not happen.
    partial: bool,
}

struct App {
    /// The scene text being edited, and where it loads from and saves to.
    text: String,
    path: String,
    /// Every text this scene has been, so an edit can be taken back.
    ///
    /// **Committed on both sides of an edit.** The text pane is bound straight to `text`, so
    /// typing never passes through `editor_core::edit` at all; committing the current text
    /// *before* a splice folds that typing into a state of its own, and committing after
    /// records the splice. Undoing a removed domain therefore walks back through what was typed
    /// rather than discarding it. See [`editor_core::edit::History`].
    history: editor_core::edit::History,
    /// The result of checking `text`, refreshed whenever the text changes.
    checked: editor_core::Checked,
    /// The last run — possibly still growing, streamed frame by frame.
    run: Option<RunView>,
    /// The verify report, its findings count, and whether its window is open.
    verify: Option<(String, usize)>,
    verify_open: bool,
    deep: bool,
    /// The in-flight background job, if any. One at a time: a second run pressed mid-run
    /// would race two worlds for one pane.
    busy: Option<(&'static str, mpsc::Receiver<Job>)>,
    /// Raised to end a streaming run between frames; replaced on every spawn.
    stop: Arc<AtomicBool>,
    status: String,
    /// The viewport camera — `viewer-core`'s, shared arithmetic with the viewer and the HTML
    /// report so all three open on the same picture.
    camera: viewer_core::Camera,
    /// Fit the camera on the next paint. Set when geometry appears or is replaced, not every
    /// frame — a camera that re-fits while you drag is a camera you cannot aim.
    needs_fit: bool,
    /// Frame the selection on the next paint. Separate from `needs_fit` because it fits a
    /// different box, and because the aspect ratio it needs only exists inside the viewport.
    pending_frame: bool,

    // The live loop.
    /// Poll the file for outside writes.
    watch: bool,
    /// Run automatically after an outside write is loaded.
    auto_run: bool,
    /// The pane has edits the disk does not. While true, outside writes are announced and
    /// never applied.
    dirty: bool,
    /// The modified time last loaded or saved, so an outside write is a *different* mtime
    /// rather than any mtime.
    known_mtime: Option<SystemTime>,
    /// When the file was last polled; polling is cheap but not free, and sixty times a
    /// second buys nothing over twice a second.
    last_poll: Option<Instant>,
    /// An auto-run is owed as soon as the current job ends — set when the file changed while
    /// something was in flight, so the picture converges on the latest text.
    rerun_owed: bool,

    // The outliner and the inspector.
    /// The selected row's path, or `None`. **A path, not an index**: the tree is rebuilt on every
    /// check and every streamed frame, and a selection keyed by row number moves under the hand
    /// that made it.
    selected: Option<String>,
    /// Paths collapsed in the outliner. Collapsed rather than expanded, so a tree that grows a
    /// branch shows it instead of hiding it.
    collapsed: std::collections::BTreeSet<String>,
    /// Paths the reader has hidden. Everything under one is hidden with it.
    hidden: std::collections::BTreeSet<String>,
    /// Draw only the selection, for reading one thing out of a crowded scene.
    solo: bool,
    // The shaded viewport.
    /// The meshes on the GPU, shared with the paint callback.
    ///
    /// Behind a lock because `egui::PaintCallback` holds an `Arc<dyn Any + Send + Sync>`: the
    /// closure cannot borrow from `App`, so the two take turns on this instead.
    gpu: Arc<Mutex<render::Shared>>,
    /// What the meshes on the GPU were built from. A drag repaints sixty times a second and the
    /// geometry has not changed; rebuilding a 200 000-face boundary each time is the difference
    /// between a viewport you can aim and one you fight.
    built: Option<u64>,
    /// Draw surfaces at all. Off falls back to the flat splat painter, which is a *different
    /// picture* and not a style — see `viewport`.
    shaded: bool,
    /// Start a run on the first update. Set by `--run`; cleared once it has fired.
    run_at_once: bool,
    /// What the cursor can be over in the shaded view, and the note each field's canvas carries.
    ///
    /// **Built with the meshes, not per paint.** The first version rebuilt every field's surface
    /// inside the readout, so a 128 cubed boundary — 200 000 faces — was meshed sixty times a
    /// second to answer a question about one pixel.
    shaded_probes: Vec<Probed>,
    shaded_notes: Vec<(String, &'static str)>,

    /// What the three handles do: move the domain along an axis, or turn it about one.
    ///
    /// One gizmo at a time, which is what every DCC tool does and for the same reason: arrows
    /// and rings occupying the same space would be six controls in the room a person aims at
    /// with one pointer.
    gizmo: Gizmo,
    /// Where on a rotation ring the drag started, in world metres.
    ///
    /// The angle a drag asks for depends on **where** the ring was taken hold of — the far side
    /// of a ring turns the other way — so the grip is remembered rather than recomputed from the
    /// pointer, which moves.
    grip: Option<[f64; 3]>,
    /// Which translate handle the pointer has hold of, as an index into [`AXES`].
    ///
    /// Held across frames because a drag is not one event: the press picks the handle, every
    /// frame after it moves the domain, and the release lets go. While it is `Some` the camera
    /// does **not** turn — a gizmo that spun the view while you dragged along it would be
    /// unusable, and the two gestures are the same gesture as far as egui is concerned.
    grabbed: Option<usize>,

    /// The framing this viewport is using, held still rather than recomputed from the geometry
    /// every frame.
    ///
    /// **The camera moves when it is asked to, and not when the scene changes shape.** Every
    /// projection here is relative to [`viewer_core::Framing::centre`], the centre of everything
    /// the scene contains, so any edit that moves or resizes geometry slides the whole picture
    /// under the reader. Dragging a handle was the case where that is fatal: moving the only
    /// domain in a scene moves the centre by the same amount, and the object does not appear to
    /// move **at all** — measured on scene 15, the box goes 3 mm, the centre goes 3 mm, apparent
    /// motion `0.000000000 m`. It reads as the box resisting the pointer, and it is the camera
    /// following it.
    ///
    /// It is the general rule rather than a special case for dragging, because the smaller
    /// version of it is the same surprise: a `cells` drag should make the block bigger, not make
    /// the room around it shrink.
    ///
    /// Cleared exactly where the camera is *asked* to move — **Fit view**, **Frame selection**,
    /// loading a file, and the first run arriving — all of which already set
    /// [`App::needs_fit`] or [`App::pending_frame`]. So there is one rule and no second list of
    /// places to keep in step.
    framing_hold: Option<viewer_core::Framing>,

    /// The value whose surface to draw instead of the object's outside, if any.
    ///
    /// `None` draws the boundary of the cells that hold a value — where the block *is*. `Some(v)`
    /// draws where the field reaches `v`, which is a different question and the one a field is
    /// usually asked: not where the block is, but where it is 100 °C. The viewport has had a depth
    /// buffer since it was written and nothing here had a surface that needed one; this is that
    /// thing, and `ARCHITECTURE.md` names it as the condition.
    iso: Option<f64>,

    /// Whether the outliner and the inspector are on screen at all.
    /// Whether the start screen is up instead of a scene.
    ///
    /// Set when the editor is opened with no file named. Before this the editor opened onto the
    /// built-in scene whatever you meant by starting it, and a scene on disk was reachable only by
    /// typing its path into the toolbar's text box. See [`crate::start`].
    start: bool,
    /// What the start screen lists, read once at startup and after each open.
    recent: Vec<crate::start::Recent>,
    /// The kinds ticked in the New-project chooser, or `None` when it is not up.
    choosing: Option<std::collections::BTreeSet<String>>,
    show_outliner: bool,
    show_inspector: bool,
    show_text: bool,
}

impl App {
    fn new(path: Option<String>) -> App {
        let (text, path, status) = match path {
            Some(p) => match std::fs::read_to_string(&p) {
                Ok(t) => (t, p.clone(), format!("loaded {p}")),
                Err(e) => (
                    default_scene(),
                    p.clone(),
                    format!("{p}: {e}; opened the built-in scene"),
                ),
            },
            None => (
                default_scene(),
                String::from("scene.json"),
                String::from("the built-in scene"),
            ),
        };
        let checked = editor_core::check(&text, &Beside::of(&path));
        let known_mtime = mtime_of(&path);
        App {
            history: editor_core::edit::History::new(text.clone()),
            text,
            path,
            checked,
            run: None,
            verify: None,
            verify_open: false,
            deep: false,
            busy: None,
            stop: Arc::new(AtomicBool::new(false)),
            status,
            camera: viewer_core::Camera::default(),
            needs_fit: true,
            pending_frame: false,
            watch: true,
            auto_run: false,
            dirty: false,
            known_mtime,
            last_poll: None,
            rerun_owed: false,
            selected: None,
            collapsed: std::collections::BTreeSet::new(),
            hidden: std::collections::BTreeSet::new(),
            solo: false,
            gpu: Arc::new(Mutex::new(render::Shared::default())),
            built: None,
            shaded: true,
            run_at_once: false,
            shaded_probes: Vec::new(),
            shaded_notes: Vec::new(),
            gizmo: Gizmo::Move,
            grip: None,
            grabbed: None,
            iso: None,
            framing_hold: None,
            start: false,
            recent: Vec::new(),
            choosing: None,
            show_outliner: true,
            show_inspector: true,
            show_text: true,
        }
    }

    /// The editor with nothing open: the start screen, and the list of what was open before.
    fn empty() -> App {
        let mut app = App::new(None);
        app.start = true;
        app.recent = crate::start::remembered();
        app.status = String::new();
        app
    }

    /// Open `path`, and stay where we are if it will not open.
    ///
    /// **The start screen is left only on a successful read.** Leaving it on the attempt would put
    /// an empty editor and an error message where the two ways forward used to be.
    fn open(&mut self, path: String) {
        let was = std::mem::replace(&mut self.path, path);
        if !self.load_from_disk() {
            self.path = was;
            return;
        }
        self.start = false;
        if let Err(e) = crate::start::remember(&self.path) {
            self.status = format!("{}; the recent list was not written: {e}", self.status);
        }
        self.recent = crate::start::remembered();
    }

    fn recheck(&mut self) {
        self.checked = editor_core::check(&self.text, &Beside::of(&self.path));
        self.run = None;
    }

    /// Step the scene back through the history, or forward again.
    ///
    /// **The current text is folded in first.** Otherwise undo pressed after typing would step
    /// back from a state the history has never seen, discarding the typing with no way to reach
    /// it again. `commit` is a no-op when nothing was typed, so the ordinary case costs a string
    /// comparison.
    ///
    /// The selection is dropped. A path into the scene is only meaningful against the text it
    /// was taken from, and undoing a removed domain can move every path after it -- keeping the
    /// old one would leave the inspector pointed at a different value than the one highlighted.
    fn take_back(&mut self, backwards: bool) {
        // **Both directions**, and the redo case is the one that needed thinking about. If the
        // user typed and then pressed redo, the buffer holds a state the history has never
        // seen; stepping forward would overwrite it. Committing first makes the typing the
        // newest state, which truncates the redo tail -- so redo then correctly reports that
        // there is nothing ahead, instead of silently throwing the typing away to get there.
        //
        // In the ordinary case, where nothing was typed, this is a string comparison and no
        // step, so redo is unaffected.
        self.history.commit(self.text.clone());
        let stepped = if backwards {
            self.history.undo()
        } else {
            self.history.redo()
        };
        let Some(text) = stepped.map(str::to_string) else {
            self.status = String::from(if backwards {
                "nothing to undo"
            } else {
                "nothing to redo"
            });
            return;
        };
        self.text = text;
        self.selected = None;
        self.dirty = true;
        self.recheck();
        self.status = format!(
            "{} — {} back, {} forward",
            if backwards { "undid" } else { "redid" },
            self.history.steps_back(),
            self.history.steps_forward()
        );
    }

    /// Start a background job, refusing a second while one is in flight.
    fn spawn(&mut self, label: &'static str, job: impl FnOnce(mpsc::Sender<Job>) + Send + 'static) {
        if self.busy.is_some() {
            self.status = format!("still busy with the last job; {label} not started");
            return;
        }
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || job(tx));
        self.busy = Some((label, rx));
        self.status = format!("{label}…");
    }

    fn start_run(&mut self) {
        let text = self.text.clone();
        // The path travels with the text, because a `parts` entry is resolved beside the scene
        // and the thread has no idea where the scene came from.
        let beside = Beside::of(&self.path);
        self.stop = Arc::new(AtomicBool::new(false));
        let stop = self.stop.clone();
        self.spawn("running", move |tx| {
            let end = editor_core::run_streaming(&text, &beside, &stop, |json| {
                let _ = tx.send(Job::Frames(json));
            });
            let _ = tx.send(Job::RunEnded(end));
        });
    }

    fn start_verify(&mut self) {
        let text = self.text.clone();
        let deep = self.deep;
        let beside = Beside::of(&self.path);
        self.spawn("verifying", move |tx| {
            let _ = tx.send(Job::Verified(editor_core::verify(&text, deep, &beside)));
        });
    }

    /// Drain everything the background job has sent. Drained, not sampled: a fast run can
    /// produce several frames per paint, and showing only the first of them would make the
    /// stream look slower than the simulation.
    fn poll_jobs(&mut self) {
        let Some((_, rx)) = &self.busy else { return };
        let mut ended = None;
        let mut latest_frames = None;
        while let Ok(job) = rx.try_recv() {
            match job {
                Job::Frames(json) => latest_frames = Some(json),
                other => {
                    ended = Some(other);
                    break;
                }
            }
        }
        if let Some(json) = latest_frames {
            match viewer_core::Run::from_json(&json) {
                Ok(run) => {
                    let last = run.frames.len().saturating_sub(1);
                    // Follow the tail while the run grows, unless the person has scrubbed
                    // back — a slider that snatches itself out of a hand is worse than one
                    // that lags.
                    let follow = self
                        .run
                        .as_ref()
                        .is_none_or(|v| v.frame + 1 >= v.run.frames.len());
                    let frame = if follow {
                        last
                    } else {
                        self.run.as_ref().map_or(last, |v| v.frame.min(last))
                    };
                    if self.run.is_none() {
                        self.needs_fit = true;
                    }
                    self.run = Some(RunView {
                        run,
                        frame,
                        playing: false,
                        last_step: None,
                        partial: true,
                    });
                    self.status = format!("running: {} frame(s) so far", last + 1);
                }
                // The editor wrote this JSON one call ago, so the viewer failing to read it
                // back is a wire-format defect, not a user mistake — say so.
                Err(e) => self.status = format!("the run's own JSON did not read back: {e}"),
            }
        }
        match ended {
            None => {}
            Some(Job::RunEnded(end)) => {
                self.busy = None;
                match end {
                    Ok(editor_core::RunEnd::Finished) => {
                        if let Some(v) = &mut self.run {
                            v.partial = false;
                            // **And show the last frame.** Pressing run means "compute this and
                            // show me", and what was on screen when a run finished was frame
                            // zero — the initial condition, which is the one frame the reader
                            // already knew. Only on a *finished* run: a stopped one is a prefix,
                            // and jumping to the end of a prefix says it ended there.
                            v.frame = v.run.frames.len().saturating_sub(1);
                            self.status = format!("ran: {} frames", v.run.frames.len());
                        }
                    }
                    Ok(editor_core::RunEnd::Stopped) => {
                        self.status =
                            String::from("stopped — what is on screen is a prefix, not the run");
                    }
                    Err(e) => {
                        // The frames already streamed stay on screen beside the reason the
                        // run ended, which is the view somebody debugging a violation wants.
                        self.status = e;
                    }
                }
                if std::mem::take(&mut self.rerun_owed) {
                    self.start_run();
                }
            }
            Some(Job::Verified(result)) => {
                self.busy = None;
                match result {
                    Ok((report, findings)) => {
                        self.verify = Some((report, findings));
                        self.verify_open = true;
                        self.status = match findings {
                            0 => String::from("verified: no structural findings"),
                            n => format!("verify: {n} finding(s) — see the report"),
                        };
                    }
                    Err(e) => self.status = e,
                }
                if std::mem::take(&mut self.rerun_owed) {
                    self.start_run();
                }
            }
            Some(Job::Frames(_)) => unreachable!("frames are drained above"),
        }
    }

    /// Notice an outside write, and apply it only when the pane has nothing to lose.
    fn poll_disk(&mut self) {
        if !self.watch {
            return;
        }
        let due = self
            .last_poll
            .is_none_or(|t| t.elapsed() > Duration::from_millis(400));
        if !due {
            return;
        }
        self.last_poll = Some(Instant::now());
        let Some(disk) = mtime_of(&self.path) else {
            return;
        };
        if Some(disk) == self.known_mtime {
            return;
        }
        if self.dirty {
            // Announced, sticky, and not applied. Which of two writers meant it is not the
            // editor's call; the disk's version waits for an explicit `load`.
            self.status = format!(
                "{} changed on disk while the pane has unsaved edits — press load to take \
                 the disk's version",
                self.path
            );
            self.known_mtime = Some(disk);
            return;
        }
        match std::fs::read_to_string(&self.path) {
            Ok(t) => {
                self.known_mtime = Some(disk);
                self.text = t;
                // A reload from disk is not an edit -- see `History::reset`.
                self.history.reset(self.text.clone());
                self.recheck();
                self.status = format!("reloaded {} (changed outside)", self.path);
                if self.auto_run && self.checked.error.is_none() {
                    if self.busy.is_some() {
                        // Converge on the latest text: end the in-flight run at the next
                        // frame boundary and owe a fresh one.
                        self.stop.store(true, Ordering::Relaxed);
                        self.rerun_owed = true;
                    } else {
                        self.start_run();
                    }
                }
            }
            Err(e) => self.status = format!("{}: {e}", self.path),
        }
    }
}

impl eframe::App for App {
    /// Give the GPU objects back while there is still a context to give them to.
    ///
    /// The only place it can be done: a `Drop` on the meshes would run after the window is gone and
    /// would have nothing to call, so the buffers and the two programs would stay allocated for the
    /// process's life. Harmless at exit and not harmless as a habit — the same reason
    /// `runtime/gpu`'s context destroys its buffers rather than leaving them to a lazy reclaim.
    fn on_exit(&mut self, gl: Option<&eframe::glow::Context>) {
        if let (Some(gl), Ok(shared)) = (gl, self.gpu.lock()) {
            shared.destroy(gl);
        }
    }

    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.ui(ctx);
    }
}

impl App {
    /// One frame of interface, from a context alone.
    ///
    /// Split out of [`eframe::App::update`] so a frame can be built **without a window**. The
    /// `eframe::Frame` argument was already unused — it is `_frame` — so nothing here ever needed
    /// eframe, and the shaded viewport is a `PaintCallback` that a headless run records and never
    /// executes, so no GL context is needed either.
    ///
    /// That is what `--ui-dump` runs. Before it existed the editor could not be looked at except
    /// by opening it, which is the same position the HTML report was in before `tools/report-check`
    /// — and this repository has already paid for that once: the viewport was handed a rect of zero
    /// width, drew nothing correctly, and cost an afternoon in the renderer before anybody printed
    /// the rect.
    fn ui(&mut self, ctx: &egui::Context) {
        // **Before anything else, and instead of it.** No menu bar, no toolbar, no panels: there is
        // nothing yet for any of them to act on, and a window of disabled controls is a worse
        // answer to "what now" than two buttons.
        // The chooser, which is the second half of the start screen rather than a window over
        // it: there is still nothing open, and the way back is a button on it.
        if let Some(ticked) = self.choosing.as_mut() {
            let mut made = crate::start::Made::Nothing;
            let mut ticked = std::mem::take(ticked);
            egui::CentralPanel::default().show(ctx, |ui| {
                made = crate::start::chooser(ui, &mut ticked);
            });
            match made {
                crate::start::Made::Back => self.choosing = None,
                crate::start::Made::Create => {
                    let kinds: Vec<&str> = ticked.iter().map(String::as_str).collect();
                    self.text = pantometry_world::templates::scene(&kinds);
                    self.history.reset(self.text.clone());
                    self.path = String::from("scene.json");
                    self.dirty = true;
                    self.known_mtime = None;
                    self.recheck();
                    self.needs_fit = true;
                    self.start = false;
                    self.choosing = None;
                    self.status = format!(
                        "a new scene of {} — save it to give it a name",
                        kinds.join(", ")
                    );
                }
                crate::start::Made::Nothing => self.choosing = Some(ticked),
            }
            return;
        }
        if self.start {
            let mut chose = crate::start::Chose::Nothing;
            egui::CentralPanel::default().show(ctx, |ui| {
                chose = crate::start::screen(ui, &self.recent);
            });
            match chose {
                crate::start::Chose::New => {
                    self.choosing = Some(std::collections::BTreeSet::new());
                }
                crate::start::Chose::Open(p) => self.open(p),
                crate::start::Chose::Nothing => {}
            }
            return;
        }

        if self.run_at_once {
            self.run_at_once = false;
            if self.checked.error.is_none() {
                self.start_run();
            } else {
                self.status = format!("--run: the scene does not check out: {}", self.status);
            }
        }
        self.poll_jobs();
        self.poll_disk();
        if self.busy.is_some() {
            // A stream draws itself; without this, frames wait for a mouse wiggle.
            ctx.request_repaint_after(Duration::from_millis(60));
        } else if self.watch {
            ctx.request_repaint_after(Duration::from_millis(400));
        }

        // A menu bar. What belongs in one is everything a reader needs occasionally and should
        // not have to find on a crowded strip: which panels are on screen, what the viewport
        // draws, and the two file operations.
        egui::TopBottomPanel::top("menu").show(ctx, |ui| {
            // Wrapped for the reason the toolbar is: `menu::bar` is a plain horizontal row, and
            // below about 160 points `Domain` was laid out at x=114..157 of the window and clipped
            // by its edge. Measured through `--ui-dump`'s `cut`, which is what that column counts.
            ui.horizontal_wrapped(|ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Open…").clicked() {
                        if let Some(p) = crate::start::pick() {
                            self.open(p);
                        }
                        ui.close_menu();
                    }
                    // **Named for what it does.** It re-reads the path already in the toolbar's
                    // box, which is a revert; called `Load`, it was the only thing on this menu
                    // that looked like a way to open a different file, and it was not one.
                    if ui.button("Revert").clicked() {
                        self.load_from_disk();
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui.button("Save").clicked() {
                        self.save_to_disk();
                        ui.close_menu();
                    }
                    if ui.button("Save as…").clicked() {
                        if let Some(p) = crate::start::pick_save(&self.path) {
                            self.path = p;
                            self.save_to_disk();
                            if let Err(e) = crate::start::remember(&self.path) {
                                self.status = format!(
                                    "{}; the recent list was not written: {e}",
                                    self.status
                                );
                            }
                            self.recent = crate::start::remembered();
                        }
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui.button("Quit").clicked() {
                        ui.ctx().send_viewport_cmd(egui::ViewportCommand::Close);
                    }
                });
                ui.menu_button("Edit", |ui| {
                    // **Undo is over the editor's own edits, not over typing.** The text pane
                    // is bound straight to the buffer and egui keeps a fine-grained undo for
                    // what it has focus over; taking that over would mean reimplementing
                    // keystroke coalescing to no benefit. What this covers is everything the
                    // editor does to the scene on the user's behalf -- a value, a string, a
                    // domain added or removed, a placement, a drag, a nudge -- and those are
                    // the ones with nothing on screen to type back.
                    let can_undo = self.history.can_undo();
                    let can_redo = self.history.can_redo();
                    if ui
                        .add_enabled(can_undo, egui::Button::new("Undo	Ctrl+Z"))
                        .clicked()
                    {
                        self.take_back(true);
                        ui.close_menu();
                    }
                    if ui
                        .add_enabled(can_redo, egui::Button::new("Redo	Ctrl+Shift+Z"))
                        .clicked()
                    {
                        self.take_back(false);
                        ui.close_menu();
                    }
                    ui.separator();
                    // **The gizmo, in Edit rather than View**, because it changes what a drag
                    // does to the scene rather than what the scene looks like. `View` is the
                    // menu of pictures.
                    ui.label("Handles");
                    if ui.radio(self.gizmo == Gizmo::Move, "Move	W").clicked() {
                        self.gizmo = Gizmo::Move;
                        self.grabbed = None;
                        self.grip = None;
                    }
                    if ui.radio(self.gizmo == Gizmo::Turn, "Turn	E").clicked() {
                        self.gizmo = Gizmo::Turn;
                        self.grabbed = None;
                        self.grip = None;
                    }
                    ui.separator();
                    ui.label(format!(
                        "{} back, {} forward",
                        self.history.steps_back(),
                        self.history.steps_forward()
                    ));
                });
                ui.menu_button("View", |ui| {
                    // **Two pictures, not two styles.** Shaded draws the boundary of what is
                    // present, lit and depth-tested, which is what a solid looks like. Cells draws
                    // every sample as a translucent splat composited far to near, which is what is
                    // *inside* it. Neither is a rendering of the other and the labels say which
                    // question each answers.
                    ui.label(egui::RichText::new("Viewport").weak().size(11.0));
                    ui.radio_value(&mut self.shaded, true, "Shaded surfaces");
                    ui.radio_value(&mut self.shaded, false, "Cells (see inside)");
                    ui.separator();
                    ui.checkbox(&mut self.show_outliner, "Outliner");
                    ui.checkbox(&mut self.show_inspector, "Inspector");
                    ui.checkbox(&mut self.show_text, "Scene text");
                    ui.separator();
                    ui.checkbox(&mut self.solo, "Solo the selection");
                    // A third picture, and the menu says what each is for the same reason the
                    // other two do: a surface says where the object is, a splat cloud says what is
                    // inside it, and an isosurface says where a *value* is. None is a rendering of
                    // either other.
                    let mut on = self.iso.is_some();
                    if ui.checkbox(&mut on, "Isosurface at a value").changed() {
                        // Opens at the middle of whatever the run's scale is, because a level
                        // outside the range draws nothing and an empty viewport reads as broken.
                        self.iso = on.then(|| self.iso_default()).flatten();
                    }
                    let (lo, hi) = self.iso_range();
                    let frame = self.frame_range();
                    if let Some(level) = &mut self.iso {
                        ui.add(egui::Slider::new(level, lo..=hi).text("level"));
                        // **What this frame actually reaches**, because the slider spans the
                        // *run* and a frame is usually narrower. Scene 15's hot spot starts 60 K
                        // above ambient and has diffused to 7 K by the last frame, so seven
                        // eighths of the slider's travel draw the same tiny blob — which reads as
                        // a broken control rather than as a level nothing has reached.
                        //
                        // The range stays the run's: rescaling per frame would make "the surface
                        // at 300 K" a different temperature at every frame, which is the mistake
                        // every other picture here avoids. So the silence is labelled instead.
                        match frame {
                            Some((flo, fhi)) => {
                                let inside = (flo..=fhi).contains(level);
                                let note = format!("this frame: {flo:.6} to {fhi:.6}");
                                if inside {
                                    ui.label(egui::RichText::new(note).weak().size(11.0));
                                } else {
                                    ui.label(
                                        egui::RichText::new(format!(
                                            "{note} — nothing reaches this level here"
                                        ))
                                        .weak()
                                        .size(11.0)
                                        .color(egui::Color32::from_rgb(230, 180, 60)),
                                    );
                                }
                            }
                            None => {
                                ui.label(
                                    egui::RichText::new("no field in this frame")
                                        .weak()
                                        .size(11.0),
                                );
                            }
                        }
                    }
                    ui.separator();
                    if ui.button("Show everything").clicked() {
                        self.hidden.clear();
                        self.solo = false;
                    }
                    ui.separator();
                    if ui.button("Fit view").clicked() {
                        self.needs_fit = true;
                        ui.close_menu();
                    }
                    if ui
                        .add_enabled(
                            self.selected.is_some(),
                            egui::Button::new("Frame selection"),
                        )
                        .clicked()
                    {
                        self.frame_selection();
                        ui.close_menu();
                    }
                });
                // Adding and removing a domain, which used to mean typing an object into the
                // text pane from memory. The list is `pantometry-world`'s, so a twentieth domain
                // appears here without this file learning about it.
                ui.menu_button("Domain", |ui| {
                    let mut wanted: Option<&str> = None;
                    ui.menu_button("Add", |ui| {
                        for t in editor_core::TEMPLATES {
                            if ui.button(t.kind).on_hover_text(t.about).clicked() {
                                wanted = Some(t.kind);
                                ui.close_menu();
                            }
                        }
                    });
                    if let Some(kind) = wanted {
                        match editor_core::add_domain(&self.text, kind) {
                            Ok(text) => {
                                self.history.commit(self.text.clone());
                                self.text = text;
                                self.history.commit(self.text.clone());
                                self.dirty = true;
                                self.recheck();
                                self.status = format!("added a {kind}");
                            }
                            Err(why) => self.status = why,
                        }
                        ui.close_menu();
                    }
                    ui.separator();
                    // Only a row that *is* a domain. `selected` can be a reading or a group, and
                    // the core answers `None` for those rather than this file deciding.
                    let selected = self
                        .selected
                        .as_deref()
                        .and_then(editor_core::domain_named)
                        .map(str::to_string);
                    let label = match &selected {
                        Some(name) => format!("Delete {name}"),
                        None => "Delete (select a domain)".to_string(),
                    };
                    if ui
                        .add_enabled(selected.is_some(), egui::Button::new(label))
                        .clicked()
                    {
                        let name = selected.expect("enabled only when there is one");
                        match editor_core::remove_domain(&self.text, &name) {
                            Ok(text) => {
                                self.history.commit(self.text.clone());
                                self.text = text;
                                self.history.commit(self.text.clone());
                                self.dirty = true;
                                self.selected = None;
                                self.recheck();
                                self.status = format!("removed {name}");
                            }
                            Err(why) => self.status = why,
                        }
                        ui.close_menu();
                    }
                });
                ui.menu_button("Run", |ui| {
                    let runnable = self.checked.error.is_none() && self.busy.is_none();
                    if ui.add_enabled(runnable, egui::Button::new("Run")).clicked() {
                        self.start_run();
                        ui.close_menu();
                    }
                    if ui
                        .add_enabled(self.busy.is_some(), egui::Button::new("Stop"))
                        .clicked()
                    {
                        self.stop.store(true, Ordering::Relaxed);
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui
                        .add_enabled(runnable, egui::Button::new("Verify"))
                        .clicked()
                    {
                        self.start_verify();
                        ui.close_menu();
                    }
                    ui.checkbox(&mut self.deep, "Deep");
                });
                ui.menu_button("Watch", |ui| {
                    ui.checkbox(&mut self.watch, "Watch file");
                    ui.checkbox(&mut self.auto_run, "Run on change");
                });
            });
        });

        egui::TopBottomPanel::top("bar").show(ctx, |ui| {
            // **Wrapped, because a row of controls that does not fit loses its end silently.**
            // Measured through `--ui-dump` at the moment it gained one: in a 700-point window
            // `fit view` was laid out at x=704..746 — outside the window, drawn, clipped by the
            // panel and so invisible and unclickable; at 500 points the last three controls
            // (`watch file`, `run on change`, `fit view`) were not laid out at all. Neither
            // showed up as an elided label, because egui does not elide a button's text: the row
            // overflows instead. Wrapping is the one-word fix that keeps every control reachable
            // at every width, and `nothing_is_drawn_past_the_right_edge` is the check that failed
            // before it and holds it now.
            ui.horizontal_wrapped(|ui| {
                // **The name, not a field.** This was an editable path 220 points wide, and it
                // was the only way into a file on disk — which is why it was a control at all.
                // `Open…` and `Save as…` are that now, so what is left is the question it also
                // answered: which file am I in. A name answers that; a path a person has to read
                // sideways does not, and the full one is a hover away.
                let name = std::path::Path::new(&self.path)
                    .file_name()
                    .map_or_else(|| self.path.clone(), |n| n.to_string_lossy().into_owned());
                ui.label(egui::RichText::new(name).strong())
                    .on_hover_text(&self.path);
                if self.dirty {
                    ui.label(egui::RichText::new("edited").weak());
                }
                ui.separator();
                let runnable = self.checked.error.is_none() && self.busy.is_none();
                if ui.add_enabled(runnable, egui::Button::new("run")).clicked() {
                    self.start_run();
                }
                if let Some(("running", _)) = self.busy {
                    if ui.button("stop").clicked() {
                        self.stop.store(true, Ordering::Relaxed);
                    }
                }
                if ui
                    .add_enabled(runnable, egui::Button::new("verify"))
                    .clicked()
                {
                    self.start_verify();
                }
                ui.checkbox(&mut self.deep, "deep");
                ui.separator();
                if ui.button("fit view").clicked() {
                    self.needs_fit = true;
                }
                // The transport, beside the run controls rather than buried under the canvas: a
                // reader stepping through a run reaches for it constantly.
                if let Some(view) = &mut self.run {
                    let frames = view.run.frames.len();
                    if frames > 1 {
                        ui.separator();
                        if ui
                            .button("|<")
                            .on_hover_text("first frame (Home)")
                            .clicked()
                        {
                            view.frame = 0;
                            view.playing = false;
                        }
                        if ui
                            .button("<")
                            .on_hover_text("previous frame (left)")
                            .clicked()
                        {
                            view.frame = (view.frame + frames - 1) % frames;
                            view.playing = false;
                        }
                        let label = if view.playing { "pause" } else { "play" };
                        if ui.button(label).on_hover_text("space").clicked() {
                            view.playing = !view.playing;
                            view.last_step = None;
                        }
                        if ui.button(">").on_hover_text("next frame (right)").clicked() {
                            view.frame = (view.frame + 1) % frames;
                            view.playing = false;
                        }
                        if ui.button(">|").on_hover_text("last frame (End)").clicked() {
                            view.frame = frames - 1;
                            view.playing = false;
                        }
                    }
                }
            });
        });

        // Keys, because comparing two frames means pressing the same thing twice and a small
        // button is the wrong target for that. Skipped while the text pane has focus, or the
        // space bar would step the run instead of typing a space into the scene.
        let typing = ctx.memory(|m| m.focused().is_some());

        // **Undo, and the focus guard is the function's own.** `shortcut` refuses while a widget
        // has focus, because a focused text field owns its undo; putting this inside the
        // `if !typing` block below would pass a constant and leave that refusal untested in
        // everything but a test.
        //
        // Outside the `frames > 1` guard too: undo is about the scene, and the whole point is
        // that it works before anything has been run.
        let step = ctx.input(|i| {
            editor_core::edit::shortcut(
                i.modifiers.command,
                i.modifiers.shift,
                i.key_pressed(egui::Key::Z),
                i.key_pressed(egui::Key::Y),
                typing,
            )
        });
        if let Some(step) = step {
            self.take_back(step == editor_core::edit::Step::Back);
        }

        // W and E, which is what a person arrives from Blender or Unreal expecting. Guarded on
        // focus like everything else here, or they would be two letters somebody cannot type.
        if !typing {
            let wanted = ctx.input(|i| {
                if i.key_pressed(egui::Key::W) {
                    Some(Gizmo::Move)
                } else if i.key_pressed(egui::Key::E) {
                    Some(Gizmo::Turn)
                } else {
                    None
                }
            });
            if let Some(g) = wanted {
                if self.gizmo != g {
                    self.gizmo = g;
                    // Dropped, not carried: a grip is a point on a ring and means nothing to an
                    // arrow, and a handle index held across a mode change would be a drag the
                    // user started on a control that is no longer there.
                    self.grabbed = None;
                    self.grip = None;
                    self.status = match g {
                        Gizmo::Move => String::from("handles: move"),
                        Gizmo::Turn => String::from("handles: turn"),
                    };
                }
            }
        }

        if !typing {
            let frames = self.run.as_ref().map_or(0, |v| v.run.frames.len());
            if frames > 1 {
                let (mut step, mut jump, mut toggle) = (0i64, None, false);
                ctx.input(|i| {
                    if i.key_pressed(egui::Key::ArrowLeft) {
                        step -= 1;
                    }
                    if i.key_pressed(egui::Key::ArrowRight) {
                        step += 1;
                    }
                    if i.key_pressed(egui::Key::Home) {
                        jump = Some(0);
                    }
                    if i.key_pressed(egui::Key::End) {
                        jump = Some(frames - 1);
                    }
                    if i.key_pressed(egui::Key::Space) {
                        toggle = true;
                    }
                });
                if let Some(view) = &mut self.run {
                    if toggle {
                        view.playing = !view.playing;
                        view.last_step = None;
                    }
                    if let Some(f) = jump {
                        view.frame = f;
                        view.playing = false;
                    }
                    if step != 0 {
                        let n = frames as i64;
                        view.frame = (((view.frame as i64 + step) % n + n) % n) as usize;
                        view.playing = false;
                    }
                }
            }
        }

        // Playback, on a clock. Sixteen frames a second: a run is not a film and a reader is
        // comparing frames, not watching one.
        if let Some(view) = &mut self.run {
            let frames = view.run.frames.len();
            if view.playing && frames > 1 {
                let due = view
                    .last_step
                    .is_none_or(|t| t.elapsed() >= Duration::from_millis(62));
                if due {
                    view.frame = (view.frame + 1) % frames;
                    view.last_step = Some(Instant::now());
                }
                ctx.request_repaint_after(Duration::from_millis(20));
            }
        }

        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            // **Right to left, and the message last.** Wrapped in a plain `horizontal` this row
            // stopped wrapping — a `horizontal` lays out on one line by definition — and the
            // narrow-window notice, which is 380 points long, was drawn past the edge of a
            // 420-point window and clipped. `--ui-dump`'s `cut` column found it in the same commit
            // that caused it. So the indicator is placed first, against the right edge, and the
            // message takes what is left and wraps into it.
            ui.with_layout(egui::Layout::right_to_left(egui::Align::TOP), |ui| {
                // **What the editor does when nobody clicks anything, said where a person looks
                // to find out what it is doing.** These were two checkboxes on the toolbar, which
                // is the wrong place for the same reason `deep` is the right one: `deep` changes
                // what a button does when you press it, so it belongs beside that button, and
                // these change what happens without anyone pressing anything. The toggles are in
                // the Watch menu; this is their state.
                let watching = match (self.watch, self.auto_run) {
                    (true, true) => Some("watching the file, running on change"),
                    (true, false) => Some("watching the file"),
                    (false, true) => Some("running on change, but not watching the file"),
                    (false, false) => None,
                };
                if let Some(what) = watching {
                    ui.label(egui::RichText::new(what).weak());
                }
                match self.checked.error.as_deref() {
                    Some(e) => {
                        ui.add(
                            egui::Label::new(
                                egui::RichText::new(e).color(egui::Color32::from_rgb(220, 80, 80)),
                            )
                            .wrap(),
                        );
                    }
                    None => {
                        ui.add(egui::Label::new(&self.status).wrap());
                    }
                }
            });
        });

        // Rebuilt every paint, from the checked scene and whatever the run has produced so far.
        // Cheap — a few hundred rows at most — and the alternative is a cache that goes stale
        // exactly when a streaming run is adding to it.
        let frame_at = self.run.as_ref().map_or(0, |v| v.frame);
        let tree = editor_core::tree(&self.checked, self.run.as_ref().map(|v| &v.run), frame_at);

        // **A width range on every panel, because `default_width` is only where they start.**
        // A `SidePanel` grows to whatever its widest child asks for, and one outliner row is a
        // scene's title: "a hot spot in aluminium, meeting a wall of borosilicate halfway: 1
        // domain(s), 0.500 s in 11 frames". Measured on a 1706-pixel window, that row took the
        // outliner to **940 px** and left the viewport nothing at all — a 3D editor showing no 3D,
        // because of a caption. The rows truncate now and the panels have a ceiling.
        //
        // **And the viewport gets a floor**, which the ceiling alone does not give it. Every
        // `SidePanel` is laid out before the `CentralPanel`, so the panels take what they ask for
        // and the viewport gets the remainder — including when the remainder is **nothing**. On the
        // default window before this it was exactly that: the paint callback was handed a rect of
        // zero width and drew nothing, correctly, looking precisely like a renderer that did not
        // work. Opening the window larger hid that rather than fixing it.
        //
        // So the panels yield, in the order they are worth least to somebody looking at a picture,
        // and the status bar says which one went. Their `show_*` flags are untouched, so widening
        // the window brings the panel back without anybody having to find the menu item again.
        let (fits, squeezed) = self.panels_that_fit(ctx.screen_rect().width());
        if let Some(dropped) = squeezed {
            self.status =
                format!("{dropped} hidden — the window is too narrow for it and the view");
        }
        // Each panel's ceiling is its minimum plus its share of what is over the floor, so the
        // three of them together can never take the viewport below `VIEWPORT_FLOOR`.
        //
        // **The ceiling alone, because `width_range` clamps the default too.** Both were written
        // — the ceiling and a `default_width` clamped to it — and sabotage found each of them
        // sufficient on its own: removing either left the panel exactly as wide, at all fourteen
        // widths, because the other did the clamping. Two mechanisms neither of which any test
        // could show to be necessary. The ceiling is the one kept, since it also holds a width a
        // reader has dragged, which `default_width` is not consulted about after the first frame.
        let share = self.panel_budget(ctx.screen_rect().width(), &fits);
        let cap = |min: f32, most: f32| (min + share).min(most);
        if fits.outliner {
            let most = cap(OUTLINER_MIN, 440.0);
            egui::SidePanel::left("outliner")
                .resizable(true)
                .default_width(260.0)
                .width_range(OUTLINER_MIN..=most)
                .show(ctx, |ui| self.outliner(ui, &tree));
        }
        if fits.inspector {
            let most = cap(INSPECTOR_MIN, 460.0);
            egui::SidePanel::right("inspector")
                .resizable(true)
                .default_width(320.0)
                .width_range(INSPECTOR_MIN..=most)
                .show(ctx, |ui| self.inspector(ui, &tree));
        }

        if fits.text {
            let most = cap(TEXT_MIN, 560.0);
            egui::SidePanel::left("text")
                .resizable(true)
                .default_width(430.0)
                .width_range(TEXT_MIN..=most)
                .show(ctx, |ui| {
                    if let Some(summary) = &self.checked.summary {
                        ui.label(summary.clone());
                        for note in &self.checked.notes {
                            ui.colored_label(
                                egui::Color32::from_rgb(230, 180, 60),
                                format!("note: {note}"),
                            );
                        }
                        ui.separator();
                    }
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        let edit = ui.add(
                            egui::TextEdit::multiline(&mut self.text)
                                .code_editor()
                                .desired_width(f32::INFINITY)
                                .desired_rows(34),
                        );
                        if edit.changed() {
                            self.dirty = true;
                            self.recheck();
                        }
                    });
                });
        }

        if self.verify_open {
            let mut open = true;
            if let Some((report, findings)) = &self.verify {
                let title = match findings {
                    0 => "verify — no structural findings".to_string(),
                    n => format!("verify — {n} FINDING(S)"),
                };
                egui::Window::new(title).open(&mut open).show(ctx, |ui| {
                    egui::ScrollArea::both().show(ui, |ui| {
                        ui.label(egui::RichText::new(report.clone()).monospace());
                    });
                });
            }
            self.verify_open = open;
        }

        egui::CentralPanel::default().show(ctx, |ui| {
            self.viewport(ui);
        });
    }
}

impl App {
    /// Load the file at `path` into the pane, replacing whatever is there.
    fn load_from_disk(&mut self) -> bool {
        match std::fs::read_to_string(&self.path) {
            Ok(t) => {
                self.text = t;
                self.history.reset(self.text.clone());
                self.dirty = false;
                self.known_mtime = mtime_of(&self.path);
                self.recheck();
                self.needs_fit = true;
                self.status = format!("loaded {}", self.path);
                true
            }
            Err(e) => {
                self.status = format!("{}: {e}", self.path);
                false
            }
        }
    }

    /// Write the pane to `path`.
    fn save_to_disk(&mut self) {
        match std::fs::write(&self.path, &self.text) {
            Ok(()) => {
                self.dirty = false;
                self.known_mtime = mtime_of(&self.path);
                self.status = format!("saved {}", self.path);
            }
            Err(e) => self.status = format!("{}: {e}", self.path),
        }
    }

    /// Put the camera on the selection's box, if it has one.
    fn frame_selection(&mut self) {
        self.pending_frame = true;
    }

    /// Whether a row is drawn, given what is hidden and whether the selection is soloed.
    ///
    /// Hiding is by **path prefix**, so hiding a domain hides its field and its readings with it,
    /// which is what a reader means by hiding a domain.
    fn visible(&self, path: &str) -> bool {
        if self
            .hidden
            .iter()
            .any(|h| editor_core::Tree::contains(h, path))
        {
            return false;
        }
        match (&self.solo, &self.selected) {
            (true, Some(sel)) => {
                editor_core::Tree::contains(sel, path) || editor_core::Tree::contains(path, sel)
            }
            _ => true,
        }
    }

    /// Whether a *panel or extent named `name`* is drawn. The viewport asks this.
    ///
    /// A name can appear in the tree twice — as a placed extent and as the panel a run filled it
    /// with — and hiding either should hide the drawing, which is one thing.
    fn draws(&self, name: &str) -> bool {
        self.visible(&format!("/run/{name}")) && self.visible(&format!("/extents/{name}"))
    }

    /// The outliner: what is in this scene, as a tree you can select, hide and collapse.
    fn outliner(&mut self, ui: &mut egui::Ui, tree: &editor_core::Tree) {
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new("Outliner").strong());
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.checkbox(&mut self.solo, "solo")
                    .on_hover_text("draw only the selection");
            });
        });
        ui.separator();
        egui::ScrollArea::vertical()
            .auto_shrink([false, false])
            .show(ui, |ui| {
                // Rows whose parent is collapsed are skipped, and so are their children: a
                // collapsed branch has to hide the whole branch, not one level of it.
                let mut hide_below: Option<usize> = None;
                for (i, node) in tree.nodes.iter().enumerate() {
                    if let Some(d) = hide_below {
                        if node.depth > d {
                            continue;
                        }
                        hide_below = None;
                    }
                    let has_children = tree.nodes.get(i + 1).is_some_and(|n| n.depth > node.depth);
                    let is_collapsed = self.collapsed.contains(&node.path);
                    if has_children && is_collapsed {
                        hide_below = Some(node.depth);
                    }
                    let selected = self.selected.as_deref() == Some(node.path.as_str());
                    let shown = self.visible(&node.path);

                    ui.horizontal(|ui| {
                        ui.add_space(node.depth as f32 * 12.0);
                        if has_children {
                            let glyph = if is_collapsed { "\u{25b8}" } else { "\u{25be}" };
                            if ui.small_button(glyph).clicked() {
                                if is_collapsed {
                                    self.collapsed.remove(&node.path);
                                } else {
                                    self.collapsed.insert(node.path.clone());
                                }
                            }
                        } else {
                            ui.add_space(18.0);
                        }
                        // The eye. Greyed rather than removed for a row hidden by an ancestor,
                        // so a reader can see *that* it is hidden and where from.
                        let own = self.hidden.contains(&node.path);
                        let eye = if own {
                            "\u{2013}"
                        } else if shown {
                            "\u{25cf}"
                        } else {
                            "\u{25cb}"
                        };
                        if ui
                            .small_button(eye)
                            .on_hover_text(if own { "show" } else { "hide" })
                            .clicked()
                        {
                            if own {
                                self.hidden.remove(&node.path);
                            } else {
                                self.hidden.insert(node.path.clone());
                            }
                        }
                        let mut text = egui::RichText::new(&node.name);
                        if !shown {
                            text = text.weak();
                        }
                        // The kind tag first, from the right, so the name gets what is left and
                        // truncates into it rather than widening the panel. The full name is on
                        // hover, which is where a name too long to show belongs.
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            ui.label(
                                egui::RichText::new(node.kind.label())
                                    .weak()
                                    .monospace()
                                    .size(10.0),
                            );
                            ui.with_layout(
                                egui::Layout::left_to_right(egui::Align::Center),
                                |ui| {
                                    // `Button::selected` with the frame only when selected,
                                    // rather than `SelectableLabel`: the two draw the same thing
                                    // and only one of them can truncate.
                                    let row = ui
                                        .add(
                                            egui::Button::new(text)
                                                .selected(selected)
                                                .frame(selected)
                                                .truncate(),
                                        )
                                        .on_hover_text(&node.name);
                                    if row.clicked() {
                                        self.selected = Some(node.path.clone());
                                    }
                                },
                            );
                        });
                    });
                }
            });
    }

    /// The range an isosurface level may take, from the run's own scale.
    ///
    /// The scale rather than this frame's spread, for the reason every picture here uses it: a
    /// slider that renormalised per frame would move the surface while the field stood still.
    fn iso_range(&self) -> (f64, f64) {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        if let Some(view) = &self.run {
            for panel in view.run.panels() {
                if let Some((a, b)) = view.run.scale_of(&panel) {
                    lo = lo.min(a);
                    hi = hi.max(b);
                }
            }
        }
        if lo.is_finite() && hi.is_finite() && hi > lo {
            (lo, hi)
        } else {
            (0.0, 1.0)
        }
    }

    /// What the *current frame's* fields actually span, as against the run's scale.
    ///
    /// The two differ by a lot whenever a field settles: a hot spot released 60 K above ambient is
    /// 7 K above it by the end, so most of a run-wide slider is outside every late frame. The
    /// slider keeps the run's range on purpose — a level has to mean one temperature across the
    /// whole run — and this is what lets the control say so rather than appear stuck.
    fn frame_range(&self) -> Option<(f64, f64)> {
        let view = self.run.as_ref()?;
        let frame = view.run.frames.get(view.frame)?;
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for panel in &frame.panels {
            if matches!(panel, viewer_core::Panel::Field { .. }) {
                for v in panel.values().iter().filter(|v| v.is_finite()) {
                    lo = lo.min(*v);
                    hi = hi.max(*v);
                }
            }
        }
        (lo.is_finite() && hi.is_finite()).then_some((lo, hi))
    }

    /// Where the slider opens: the middle of the range, or nothing when there is no run to have a
    /// range. Turning it on before running would draw an empty viewport, which reads as broken.
    fn iso_default(&self) -> Option<f64> {
        self.run.as_ref()?;
        // The middle of **this frame**, not of the run. The run's range is what the slider spans,
        // but opening at the middle of it lands outside every late frame of a settling field —
        // scene 15's run is 60 K wide and its last frame is 7 K, so the midpoint reaches nothing
        // and the control would turn on showing an empty viewport. Which is the exact failure
        // this function was written to avoid, committed by the function itself.
        let (lo, hi) = self.frame_range().unwrap_or_else(|| self.iso_range());
        Some(0.5 * (lo + hi))
    }

    /// The selected row's domain, if the selection is one.
    fn selected_domain(&self) -> Option<String> {
        self.selected
            .as_deref()
            .and_then(editor_core::domain_named)
            .map(str::to_string)
    }

    /// Where this domain's handles sit: the centre of its placed box.
    ///
    /// The centre rather than the origin the pose states, because a pose of `[0, 0, 0]` puts the
    /// origin at the world origin, which for most scenes is a corner of the geometry and for some
    /// is outside it entirely. The handles are a thing to grab, and what they move is a delta, so
    /// where they are drawn is a matter of being reachable rather than of being correct.
    fn handle_origin(&self, name: &str) -> Option<[f64; 3]> {
        let placed = self.checked.boxes.iter().find(|b| b.name == name)?;
        let mut centre = [0.0; 3];
        for corner in &placed.corners {
            for (c, v) in centre.iter_mut().zip(corner) {
                *c += v / 8.0;
            }
        }
        Some(centre)
    }

    /// The inspector: everything the core knows about the selected row.
    ///
    /// The rows come from `editor_core::Node::detail` rather than being assembled here, because
    /// "what is this thing" is a question two shells should not answer differently.
    fn inspector(&mut self, ui: &mut egui::Ui, tree: &editor_core::Tree) {
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new("Inspector").strong());
        });
        ui.separator();
        let Some(path) = self.selected.clone() else {
            ui.weak(
                "nothing selected — pick a row in the outliner, or click something in the viewport",
            );
            return;
        };
        let Some(node) = tree.find(&path) else {
            // A selection that no longer exists is worth saying rather than blanking: it means
            // the scene was edited under it, which is a thing the reader did.
            ui.weak(format!("{path} is no longer in this scene"));
            return;
        };
        ui.label(egui::RichText::new(&node.name).heading());
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new(node.kind.label()).monospace().weak());
            ui.label(
                egui::RichText::new(&node.path)
                    .monospace()
                    .weak()
                    .size(10.0),
            );
        });
        if node.bounds.is_some() && ui.button("frame this").clicked() {
            self.frame_selection();
        }
        ui.separator();

        // What the scene *states* about this row, and can be changed. Above the read-only detail
        // because it is the half a person came here to act on, and separated from it because the
        // difference between "this is what you asked for" and "this is what came out" is the
        // difference between an input and a result — an inspector that mixed them would invite
        // dragging a peak temperature.
        let fields = editor_core::editable(&self.text, &path);
        // One change per frame, applied after the grid. Collected rather than applied in place
        // because the closure cannot hold `&mut self` while `self.text` is being read for the
        // widgets, and because a splice mid-iteration would invalidate every pointer after it.
        let mut change: Option<Edited> = None;
        if !fields.is_empty() {
            ui.label(egui::RichText::new("scene").weak().size(11.0));
            egui::Grid::new("scene values")
                .num_columns(2)
                .striped(true)
                .spacing([10.0, 3.0])
                .show(ui, |ui| {
                    for field in &fields {
                        let label = if field.unit.is_empty() {
                            field.label.clone()
                        } else {
                            format!("{} ({})", field.label, field.unit)
                        };
                        ui.label(egui::RichText::new(label).weak().monospace().size(11.0));
                        match &field.value {
                            editor_core::Value::Number { now, integral } => {
                                let mut value = *now;
                                // A rate, not a range. Nothing here clamps: the scene's own check
                                // is the authority on what is legal and it runs on every change
                                // already, so a limit invented in the shell would be a second
                                // opinion that can only be wrong — and wrong in the direction of
                                // refusing a value the format takes.
                                let drag = if *integral {
                                    egui::DragValue::new(&mut value)
                                        .speed(1.0)
                                        .fixed_decimals(0)
                                } else {
                                    egui::DragValue::new(&mut value)
                                        .speed(now.abs().max(1.0) * 0.01)
                                };
                                if ui.add(drag).changed() {
                                    change = Some(Edited::Number(field.pointer.clone(), value));
                                }
                            }
                            // A menu where the format has a set and a text box where it does not.
                            // The menu is not a restriction: `choices` is the catalogue plus the
                            // scene's own declarations, so anything the file could already say is
                            // in it, and a key with no known set stays free text rather than
                            // being offered a list this shell invented.
                            editor_core::Value::Text { now, choices } if !choices.is_empty() => {
                                let mut picked = now.clone();
                                egui::ComboBox::from_id_salt(&field.pointer)
                                    .selected_text(now.as_str())
                                    .show_ui(ui, |ui| {
                                        for choice in choices {
                                            ui.selectable_value(
                                                &mut picked,
                                                choice.clone(),
                                                choice.as_str(),
                                            );
                                        }
                                    });
                                if picked != *now {
                                    change = Some(Edited::Text(field.pointer.clone(), picked));
                                }
                            }
                            editor_core::Value::Text { now, .. } => {
                                let mut typed = now.clone();
                                if ui
                                    .add(
                                        egui::TextEdit::singleline(&mut typed)
                                            .desired_width(f32::INFINITY),
                                    )
                                    .changed()
                                {
                                    change = Some(Edited::Text(field.pointer.clone(), typed));
                                }
                            }
                        }
                        ui.end_row();
                    }
                });
            ui.separator();
        }

        // Where the domain sits, which is the scene's business rather than the domain's: `poses`
        // is a map beside `materials`, not a field inside the object above. So it is its own
        // control with its own write, and it appears for a domain row whether or not the file has
        // ever said anything about placement — no shipped scene does.
        let mut moved: Option<(String, [f64; 3])> = None;
        let mut turned: Option<(String, [f64; 3], f64)> = None;
        if let Some(domain) = editor_core::domain_named(&path).map(str::to_string) {
            let at = editor_core::pose_of(&self.text, &domain);
            let mut to = at;
            // A millimetre a pixel would be wrong for a room and right for a die, so the rate
            // comes from what this scene is the size of. Still a rate and not a limit: nothing
            // here bounds where a domain may go.
            let speed = self
                .checked
                .bounds
                .map(|b| ((b[3] - b[0]).max(b[4] - b[1]).max(b[5] - b[2]) / 400.0).max(1e-6))
                .unwrap_or(1e-3);
            ui.label(egui::RichText::new("placement").weak().size(11.0));
            egui::Grid::new("placement")
                .num_columns(2)
                .striped(true)
                .spacing([10.0, 3.0])
                .show(ui, |ui| {
                    for (i, axis) in ["x", "y", "z"].iter().enumerate() {
                        ui.label(
                            egui::RichText::new(format!("at_m.{axis} (m)"))
                                .weak()
                                .monospace()
                                .size(11.0),
                        );
                        if ui
                            .add(egui::DragValue::new(&mut to[i]).speed(speed))
                            .changed()
                        {
                            moved = Some((domain.clone(), to));
                        }
                        ui.end_row();
                    }
                    // The turn, beside the position because they are one key in the file. Once it
                    // exists the inspector's generic walk offers it from the scene row as well —
                    // this is the control that brings it into being, which is the only part that
                    // needed writing.
                    let (axis, degrees) = editor_core::turn_of(&self.text, &domain);
                    let (mut spin, mut about) = (degrees, axis);
                    ui.label(
                        egui::RichText::new("turn (deg)")
                            .weak()
                            .monospace()
                            .size(11.0),
                    );
                    if ui.add(egui::DragValue::new(&mut spin).speed(1.0)).changed() {
                        turned = Some((domain.clone(), about, spin));
                    }
                    ui.end_row();
                    for (i, name) in ["axis.x", "axis.y", "axis.z"].iter().enumerate() {
                        ui.label(egui::RichText::new(*name).weak().monospace().size(11.0));
                        if ui
                            .add(egui::DragValue::new(&mut about[i]).speed(0.05))
                            .changed()
                        {
                            turned = Some((domain.clone(), about, spin));
                        }
                        ui.end_row();
                    }
                });
            ui.separator();
        }
        if let Some((domain, to)) = moved {
            match editor_core::set_pose(&self.text, &domain, to) {
                Ok(text) => {
                    // Fold in anything typed since the last edit, then record this one. The
                    // first commit is a no-op when nobody typed, by construction.
                    self.history.commit(self.text.clone());
                    self.text = text;
                    self.history.commit(self.text.clone());
                    self.dirty = true;
                    self.recheck();
                }
                Err(why) => self.status = why,
            }
        }
        if let Some((domain, axis, degrees)) = turned {
            // An axis dragged through zero on every component is refused by the format rather
            // than normalised into a `NaN`, and the status line says so instead of the file
            // taking it. Dragging back out of zero recovers.
            match editor_core::set_turn(&self.text, &domain, axis, degrees) {
                Ok(text) => {
                    // Fold in anything typed since the last edit, then record this one. The
                    // first commit is a no-op when nobody typed, by construction.
                    self.history.commit(self.text.clone());
                    self.text = text;
                    self.history.commit(self.text.clone());
                    self.dirty = true;
                    self.recheck();
                }
                Err(why) => self.status = why,
            }
        }
        {}
        if let Some(edit) = change {
            let spliced = match &edit {
                Edited::Number(pointer, value) => {
                    editor_core::set_number(&self.text, pointer, *value)
                }
                Edited::Text(pointer, value) => editor_core::set_text(&self.text, pointer, value),
            };
            match spliced {
                Ok(text) => {
                    // Fold in anything typed since the last edit, then record this one. The
                    // first commit is a no-op when nobody typed, by construction.
                    self.history.commit(self.text.clone());
                    self.text = text;
                    self.history.commit(self.text.clone());
                    self.dirty = true;
                    self.recheck();
                }
                // A refused change says why and leaves the file alone. The reachable case is a
                // fraction dragged onto a count, which the widget's own step of one already makes
                // hard — but "hard to reach" is not "cannot happen", and a splice that half
                // applied would be the worst outcome available here.
                Err(why) => self.status = why,
            }
        }

        egui::ScrollArea::vertical()
            .auto_shrink([false, false])
            .show(ui, |ui| {
                egui::Grid::new("props")
                    .num_columns(2)
                    .striped(true)
                    .spacing([10.0, 3.0])
                    .show(ui, |ui| {
                        for (k, v) in &node.detail {
                            ui.label(egui::RichText::new(k).weak().monospace().size(11.0));
                            ui.label(egui::RichText::new(v).monospace().size(11.0));
                            ui.end_row();
                        }
                    });
            });
    }

    /// What the meshes on the GPU depend on, as one number.
    ///
    /// Everything that changes the *geometry* and nothing that does not: the camera is a uniform and
    /// the frame counter is not, so a drag rebuilds nothing and a step rebuilds everything.
    ///
    /// A hash rather than a comparison of the parts, because the parts are a scene's text and a set
    /// of names and the alternative is keeping a copy of both. `DefaultHasher` is deterministic
    /// within a process, which is all a cache key needs — this is not a result, and nothing is
    /// pinned to it.
    fn geometry_key(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut h = std::collections::hash_map::DefaultHasher::new();
        self.text.hash(&mut h);
        self.hidden.hash(&mut h);
        self.solo.hash(&mut h);
        self.selected.hash(&mut h);
        // **The isosurface level changes the mesh, so it belongs in the key.** Left out, the
        // slider moved and the viewport kept the geometry it had already built — the control
        // worked, the mesher worked, and nothing happened. Hashed as bits because `f64` is not
        // `Hash`, and `None` as a distinct byte so "off" is not the same key as some level.
        match self.iso {
            Some(level) => {
                1u8.hash(&mut h);
                level.to_bits().hash(&mut h);
            }
            None => 0u8.hash(&mut h),
        }
        if let Some(view) = &self.run {
            view.frame.hash(&mut h);
            view.run.frames.len().hash(&mut h);
            // A streaming run replaces its frames as it goes, so the count alone would leave the
            // last frame's geometry stale by one. The final frame's time settles it.
            if let Some(f) = view.run.frames.last() {
                f.t.to_bits().hash(&mut h);
            }
        } else {
            0u8.hash(&mut h);
        }
        h.finish()
    }

    /// The box everything on screen sits in: the scene's placed extents and the run's panels.
    ///
    /// **World, not each panel's own.** `Panel::bounds` is the box in whichever frame the panel
    /// is in, and a union of boxes from different frames is a box around nothing -- which is what
    /// the camera fitted to for as long as a run could state a placement and this could not read
    /// one.
    ///
    /// A method rather than four lines inside `viewport` so that `drawn_extent` frames what the
    /// viewport frames. A hook that computed its own world would be checking a second
    /// implementation and would agree with the first exactly until one of them changed.
    fn world(&self) -> Option<[f64; 6]> {
        let mut bounds = self.checked.bounds;
        if let Some(view) = &self.run {
            if let Some(frame) = view.run.frames.first() {
                for panel in &frame.panels {
                    bounds = Some(union(bounds, panel.world_bounds()));
                }
            }
        }
        bounds
    }

    /// The two batches the shaded pass draws: surfaces, and lines.
    ///
    /// Framing-local, in `f64` until the last step — see `render`'s module docs on why the centre
    /// must not be folded into the matrix instead.
    #[allow(clippy::type_complexity)]
    fn batches(
        &self,
        framing: &viewer_core::Framing,
    ) -> (
        render::Batch,
        render::Batch,
        Vec<Probed>,
        Vec<(String, &'static str)>,
    ) {
        let mut solid = render::Batch::new();
        let mut lines = render::Batch::new();
        let mut probes: Vec<Probed> = Vec::new();
        let mut notes: Vec<(String, &'static str)> = Vec::new();

        // The scene's placed extents, as wire boxes. Depth-tested with everything else, so a box
        // behind a solid is behind it.
        let wire = [0.42f32, 0.45, 0.5];
        let picked_colour = [1.0f32, 0.62, 0.12];
        let selected_name = self.selected.as_deref().and_then(|p| {
            p.strip_prefix("/run/")
                .or_else(|| p.strip_prefix("/extents/"))
                .or_else(|| p.strip_prefix("/readings/"))
                .map(|rest| rest.split('/').next().unwrap_or(rest).to_string())
        });
        for placed in &self.checked.boxes {
            if !self.draws(&placed.name) {
                continue;
            }
            let c = if selected_name.as_deref() == Some(placed.name.as_str()) {
                picked_colour
            } else {
                wire
            };
            let base = lines.vertices();
            for corner in placed.corners {
                lines.push(framing.local(corner), [0.0, 0.0, 1.0], c);
            }
            for (a, b) in editor_core::EDGES {
                lines.indices.push(base + a as u32);
                lines.indices.push(base + b as u32);
            }
        }

        // **The designed parts**, in the same place their voxels will be.
        //
        // Everything else in this viewport is a picture of a *run*: a field's cells, a level
        // through its values, a body's sphere. A part is a picture of the **scene**, so it draws
        // before there is a run and is the only three-dimensional thing on screen while somebody
        // is still authoring — which is most of the time an editor is open.
        //
        // Drawn only for a domain the current frame has no panel for. With a run there is
        // already a surface in that space and two surfaces in one place z-fight into a picture
        // of neither; the run wins, because a reader looking at results is asking about the
        // values and not about the shape they were sampled on.
        let running: std::collections::BTreeSet<&str> = self
            .run
            .as_ref()
            .and_then(|v| {
                v.run
                    .frames
                    .get(v.frame.min(v.run.frames.len().saturating_sub(1)))
            })
            .map(|f| f.panels.iter().map(|p| p.name()).collect())
            .unwrap_or_default();
        for m in &self.checked.meshes {
            if !self.draws(&m.name) || running.contains(m.name.as_str()) {
                continue;
            }
            let c = if selected_name.as_deref() == Some(m.name.as_str()) {
                picked_colour
            } else {
                // Deliberately unlike any colour scale in this workspace. A designed part
                // carries no value, and giving it one from the same palette a field uses would
                // invite reading a temperature off a thing that has none.
                [0.55f32, 0.57, 0.60]
            };
            let surface = pantometry::view::mesh::mesh_surface(&m.triangles, 0);
            let base = solid.vertices();
            for (p, n) in surface.positions.iter().zip(&surface.normals) {
                // The triangles are the STL's own; `m.place` is where the scene puts them. They
                // used to arrive here already in world metres, which drew correctly and was the
                // second convention in a file that has one.
                let at = editor_core::place(
                    m.place,
                    [f64::from(p[0]), f64::from(p[1]), f64::from(p[2])],
                );
                solid.push(framing.local(at), *n, c);
            }
            solid
                .indices
                .extend(surface.indices.iter().map(|i| base + i));
            notes.push((m.name.clone(), "part: as designed, not as rasterised"));
        }

        let Some(view) = &self.run else {
            return (solid, lines, probes, notes);
        };
        let Some(frame) = view
            .run
            .frames
            .get(view.frame.min(view.run.frames.len().saturating_sub(1)))
        else {
            return (solid, lines, probes, notes);
        };

        for panel in &frame.panels {
            if !self.draws(panel.name()) {
                continue;
            }
            let scale = view.run.scale_of(panel.name());
            match panel {
                viewer_core::Panel::Field {
                    name,
                    unit,
                    nx,
                    ny,
                    nz,
                    values,
                    ..
                } => {
                    // The run's own extent first, the scene's placed box second — the same order
                    // the flat path uses, and for the same reason: a run opened without the file
                    // that produced it still knows where it was.
                    // `Panel::placed_corners` and not the raw extent: the extent is the
                    // domain's own box and the placement is where that box is. Eight corners rather
                    // than a box, so a turned part stays the shape it is — `field_shell` maps a unit
                    // cube through them and never needed it to be axis-aligned.
                    let from_run = panel
                        .placed_corners()
                        .map(|corners| editor_core::PlacedBox {
                            name: name.clone(),
                            corners,
                        });
                    let Some(b) = from_run
                        .as_ref()
                        .or_else(|| self.checked.boxes.iter().find(|b| &b.name == name))
                    else {
                        continue;
                    };
                    let shell = editor_core::field_shell(
                        &b.corners,
                        (*nx, *ny, *nz),
                        values,
                        unit,
                        scale,
                        self.iso,
                    );
                    let base = solid.vertices();
                    for i in 0..shell.positions.len() {
                        solid.push(
                            framing.local(shell.positions[i]),
                            shell.normals[i],
                            shell.colours[i],
                        );
                    }
                    solid.indices.extend(shell.indices.iter().map(|i| base + i));

                    notes.push((name.clone(), shell.note));
                    // The readout offers the **surface's** vertices, because that is what is on
                    // screen: naming an interior cell a reader cannot see would answer about a
                    // different object.
                    let path = format!("/run/{name}");
                    for (at, cell) in shell.probes() {
                        let c = cell as usize;
                        let (i, j, k) = (
                            c % (*nx).max(1),
                            c / (*nx).max(1) % (*ny).max(1),
                            c / (nx * ny).max(1),
                        );
                        let v = values.get(c).copied().unwrap_or(f64::NAN);
                        probes.push(Probed {
                            at,
                            path: path.clone(),
                            label: format!(
                                "{name} [{i},{j},{k}]   {} {unit}",
                                editor_core::magnitude(v)
                            ),
                        });
                    }
                }
                viewer_core::Panel::Points {
                    positions, values, ..
                } => {
                    let pts = panel.placed_positions();
                    if pts.is_empty() {
                        continue;
                    }
                    let _ = positions;
                    // The radius the exports use, so a body is the same size in the viewport, in
                    // Blender and in usdview. `mesh::body_radius` sizes it from the run's own
                    // bounds rather than from a constant, which is why an orbit and a block do not
                    // both get a dot.
                    // World positions, so the box paired with them is the world one.
                    // `body_radius` reads the box only as a fallback when there are fewer than two
                    // bodies -- the median nearest-neighbour distance it normally uses is invariant
                    // under a rigid motion -- but pairing a world point set with a local box is the
                    // kind of mismatch this whole change is about.
                    let bounds = panel.world_bounds();
                    let radius = pantometry::view::mesh::body_radius(&pts, &bounds);
                    let colouring = editor_core::Colouring::of(panel.unit(), values, scale);
                    for (i, centre) in pts.iter().enumerate() {
                        let sphere = pantometry::view::mesh::body_spheres(&[*centre], radius);
                        let colour = colouring.linear(values[i]);
                        let base = solid.vertices();
                        for (p, n) in sphere.positions.iter().zip(&sphere.normals) {
                            solid.push(
                                framing.local([p[0] as f64, p[1] as f64, p[2] as f64]),
                                *n,
                                colour,
                            );
                        }
                        solid
                            .indices
                            .extend(sphere.indices.iter().map(|k| base + k));
                    }
                }
                viewer_core::Panel::Paths {
                    starts,
                    vertices,
                    values,
                    ..
                } => {
                    let colouring = editor_core::Colouring::of(panel.unit(), values, scale);
                    let placed = panel.placed_vertices();
                    let _ = vertices;
                    for (r, value) in values.iter().enumerate() {
                        let lo = starts[r] as usize;
                        let hi = starts.get(r + 1).map_or(placed.len(), |s| *s as usize);
                        let colour = colouring.linear(*value);
                        let base = lines.vertices();
                        for at in placed.iter().take(hi).skip(lo) {
                            lines.push(framing.local(*at), [0.0, 0.0, 1.0], colour);
                        }
                        for step in 0..hi.saturating_sub(lo).saturating_sub(1) {
                            lines.indices.push(base + step as u32);
                            lines.indices.push(base + step as u32 + 1);
                        }
                    }
                }
            }
        }
        (solid, lines, probes, notes)
    }

    /// Which side panels there is room for, and the first one that had to go.
    ///
    /// Widths are the minimums the panels are held to below; [`VIEWPORT_FLOOR`] is what the
    /// viewport must keep.
    fn panels_that_fit(&self, width: f32) -> (Panels, Option<&'static str>) {
        const FLOOR: f32 = VIEWPORT_FLOOR;
        let mut fits = Panels {
            outliner: self.show_outliner,
            text: self.show_text,
            inspector: self.show_inspector,
        };
        // The viewport has a gap of its own, on the side no panel is against: with one panel up
        // the chrome measured 16 points, not 8, and the viewport settled at 272 against a floor
        // of 280. `+ gap` unconditionally is that side.
        let gap = panel_gap();
        let taken = |p: &Panels| {
            gap + (if p.outliner { OUTLINER_MIN + gap } else { 0.0 })
                + (if p.text { TEXT_MIN + gap } else { 0.0 })
                + (if p.inspector {
                    INSPECTOR_MIN + gap
                } else {
                    0.0
                })
        };
        let mut dropped = None;

        // **Three steps, written out.** Least useful first for somebody looking at a picture: the
        // inspector describes the selection, which the viewport is already showing; the text is the
        // scene, which is not what a view is for; the outliner is how you *choose* what to look at,
        // so it goes last. A loop over closures here was a type clippy refused and three constants
        // nobody can get wrong quietly are what the twelve edges of a box are written out for.
        if width - taken(&fits) < FLOOR && fits.inspector {
            fits.inspector = false;
            dropped = Some("the inspector");
        }
        if width - taken(&fits) < FLOOR && fits.text {
            fits.text = false;
            dropped = dropped.or(Some("the scene text"));
        }
        if width - taken(&fits) < FLOOR && fits.outliner {
            fits.outliner = false;
            dropped = dropped.or(Some("the outliner"));
        }
        (fits, dropped)
    }

    /// What each shown side panel may grow to, so that the viewport keeps its floor.
    ///
    /// **`panels_that_fit` reserved the minimums and the panels spent their defaults.** Measured
    /// through `--ui-dump`, which reports the callback's own rect: at a 900-point window the
    /// arithmetic said `view=290` and the rect was **0 x 870** — the panels had taken 260, 430 and
    /// 320 against the 170, 240 and 200 that were budgeted, and the viewport got the remainder,
    /// which was nothing. Zero at 1000, 950, 900 and 700 points. That is the same zero-width rect
    /// the panel comment above was written about; dropping a panel had made it rarer and the sum
    /// still had the hole in it.
    ///
    /// So the budget is handed to the spenders. Every shown panel keeps its minimum and takes an
    /// equal share of what is left over the floor, which makes the sum of the maxima exactly
    /// `width - VIEWPORT_FLOOR`: the floor now holds however the panels are dragged, rather than
    /// only while nobody drags them. `width_range` clamps a remembered width into the range each
    /// frame, so narrowing the window pulls a panel in rather than eating the view.
    fn panel_budget(&self, width: f32, fits: &Panels) -> f32 {
        let mins = (if fits.outliner { OUTLINER_MIN } else { 0.0 })
            + (if fits.text { TEXT_MIN } else { 0.0 })
            + (if fits.inspector { INSPECTOR_MIN } else { 0.0 });
        let shown =
            usize::from(fits.outliner) + usize::from(fits.text) + usize::from(fits.inspector);
        if shown == 0 {
            return 0.0;
        }
        let gaps = (shown as f32 + 1.0) * panel_gap();
        (width - VIEWPORT_FLOOR - mins - gaps).max(0.0) / shown as f32
    }

    /// The 3D viewport: wireframe extents from the text, run panels over them by shape.
    fn viewport(&mut self, ui: &mut egui::Ui) {
        let (response, painter) =
            ui.allocate_painter(ui.available_size(), egui::Sense::click_and_drag());
        let rect = response.rect;
        let aspect = (rect.width() / rect.height().max(1.0)) as f64;

        // **The world this paint frames**, hoisted above the camera because the handles need it
        // and the handles have to answer before the camera does: one drag is either turning the
        // view or moving a domain, and whichever claims it, the other must not also act on it.
        let bounds = self.world();

        // Asking the camera to fit is asking the framing to follow the geometry again. Done here,
        // above both readers of it, so the handles and the paint cannot disagree within a frame.
        if self.needs_fit || self.pending_frame {
            self.framing_hold = None;
        }

        // The translate handles, if a domain is selected and the scene has somewhere to put them.
        let mut moved: Option<(String, [f64; 3])> = None;
        // A turn is a third thing the drag can produce: a name, the world axis it was about, and
        // how far. Applied after the hit test for the same reason `moved` is -- the splice needs
        // `self.text` while the projection above still borrows `self.camera`.
        let mut turned: Option<(String, [f64; 3], f64)> = None;
        let mut holding = false;
        if let (Some(b), Some(name)) = (bounds, self.selected_domain()) {
            let framing = self
                .framing_hold
                .unwrap_or_else(|| viewer_core::Framing::of(b));
            if let Some(origin) = self.handle_origin(&name) {
                let to_screen = |p: [f64; 3]| {
                    let q = self.camera.project(p, &framing, aspect);
                    egui::pos2(
                        rect.center().x + (q.x as f32) * rect.width() * 0.5,
                        rect.center().y - (q.y as f32) * rect.height() * 0.5,
                    )
                };
                let reach = framing.span * 0.35;
                let tip = |axis: [f64; 3]| {
                    [
                        origin[0] + axis[0] * reach,
                        origin[1] + axis[1] * reach,
                        origin[2] + axis[2] * reach,
                    ]
                };

                if response.drag_started() {
                    if let Some(at) = response.interact_pointer_pos() {
                        match self.gizmo {
                            Gizmo::Move => {
                                let base = to_screen(origin);
                                self.grabbed = AXES
                                    .iter()
                                    .enumerate()
                                    .map(|(i, (axis, _))| {
                                        (i, near_segment(at, base, to_screen(tip(*axis))))
                                    })
                                    .filter(|(_, d)| *d <= GRAB_RADIUS)
                                    .min_by(|a, b| a.1.total_cmp(&b.1))
                                    .map(|(i, _)| i);
                                self.grip = None;
                            }
                            Gizmo::Turn => {
                                // Nearest segment of any ring, and **which point of it** —
                                // the far side of a ring turns the other way, so the grip is
                                // part of the answer rather than a detail of finding it.
                                let mut best: Option<(f32, usize, [f64; 3])> = None;
                                for (i, (axis, _)) in AXES.iter().enumerate() {
                                    let pts = editor_core::edit::ring_points(origin, *axis, reach);
                                    for w in 0..pts.len() {
                                        let a = pts[w];
                                        let b = pts[(w + 1) % pts.len()];
                                        let d = near_segment(at, to_screen(a), to_screen(b));
                                        if d <= GRAB_RADIUS && best.is_none_or(|(bd, _, _)| d < bd)
                                        {
                                            best = Some((d, i, a));
                                        }
                                    }
                                }
                                self.grabbed = best.map(|(_, i, _)| i);
                                self.grip = best.map(|(_, _, p)| p);
                            }
                        }
                    }
                }
                if let Some(i) = self.grabbed {
                    holding = true;
                    if response.dragged() {
                        let d = response.drag_delta();
                        // Points to normalised device: the window is two across in both, and y
                        // runs the other way. Same mapping `project`'s screen conversion uses,
                        // read backwards.
                        let ndc = [
                            d.x as f64 / (rect.width() as f64 * 0.5).max(1.0),
                            -d.y as f64 / (rect.height() as f64 * 0.5).max(1.0),
                        ];
                        let axis = AXES[i].0;
                        match (self.gizmo, self.grip) {
                            (Gizmo::Move, _) => {
                                let t = editor_core::drag_along_axis(
                                    &self.camera,
                                    &framing,
                                    aspect,
                                    origin,
                                    axis,
                                    ndc,
                                );
                                if t != 0.0 {
                                    moved = Some((
                                        name.clone(),
                                        [axis[0] * t, axis[1] * t, axis[2] * t],
                                    ));
                                }
                            }
                            (Gizmo::Turn, Some(grip)) => {
                                let radians = editor_core::edit::turn_about_axis(
                                    &self.camera,
                                    &framing,
                                    aspect,
                                    origin,
                                    axis,
                                    grip,
                                    ndc,
                                );
                                if radians != 0.0 {
                                    turned = Some((name.clone(), axis, radians.to_degrees()));
                                }
                            }
                            // Held a ring and lost the grip: refuse rather than guess one.
                            (Gizmo::Turn, None) => {}
                        }
                    }
                }
            }
        }
        if response.drag_stopped() {
            self.grabbed = None;
            self.grip = None;
        }
        // The pose is absolute in the file and the drag is a delta, so it is read, added to and
        // written back each frame. Applied after the hit test rather than inside it, because the
        // splice needs `self.text` and the projection above is still borrowing `self.camera`.
        if let Some((name, delta)) = moved {
            let at = editor_core::pose_of(&self.text, &name);
            let to = [at[0] + delta[0], at[1] + delta[1], at[2] + delta[2]];
            match editor_core::set_pose(&self.text, &name, to) {
                Ok(text) => {
                    // Fold in anything typed since the last edit, then record this one. The
                    // first commit is a no-op when nobody typed, by construction.
                    self.history.commit(self.text.clone());
                    self.text = text;
                    self.history.commit(self.text.clone());
                    self.dirty = true;
                    self.recheck();
                }
                Err(why) => self.status = why,
            }
        }

        // **A turn composes rather than replaces.** The file states one axis and one angle, and
        // the ring asks for a rotation about a *world* axis applied on top of whatever is there;
        // rotations do not commute, so the two are composed as rotations rather than added.
        // `editor_core::edit::compose_turns` is where that lives and where its tests are.
        if let Some((name, axis, degrees)) = turned {
            let now = editor_core::turn_of(&self.text, &name);
            let (next_axis, next_deg) = editor_core::edit::compose_turns(now, (axis, degrees));
            match editor_core::set_turn(&self.text, &name, next_axis, next_deg) {
                Ok(text) => {
                    self.history.commit(self.text.clone());
                    self.text = text;
                    self.history.commit(self.text.clone());
                    self.dirty = true;
                    self.recheck();
                }
                Err(why) => self.status = why,
            }
        }

        if response.dragged() && !holding {
            let d = response.drag_delta();
            self.camera.turn(d.x as f64 * 0.01, d.y as f64 * 0.01);
        }
        if response.hovered() {
            let scroll = ui.input(|i| i.raw_scroll_delta.y) as f64;
            if scroll != 0.0 {
                self.camera.zoom(0.999_f64.powf(scroll * 3.0));
            }
        }
        // Where the cursor is, for the probe. `None` while dragging: a reader turning the model
        // is not asking what is under their finger, and a readout that flickers through every
        // body on the way is noise.
        let pointer = if response.dragged() {
            None
        } else {
            response.hover_pos()
        };
        let mut probe = Probe::default();
        // What the selection is, resolved to a name the panels can be matched against: the
        // outliner's paths are `/run/<name>` and `/extents/<name>`, and the viewport draws by
        // name.
        let selected_name = self.selected.as_deref().and_then(|p| {
            p.strip_prefix("/run/")
                .or_else(|| p.strip_prefix("/extents/"))
                .or_else(|| p.strip_prefix("/readings/"))
                .map(|rest| rest.split('/').next().unwrap_or(rest).to_string())
        });

        let Some(bounds) = bounds else {
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "nothing in this scene has geometry — sources, lumps and networks are \
                 readings, not places; run to see what they report",
                egui::FontId::proportional(14.0),
                ui.visuals().weak_text_color(),
            );
            return;
        };
        let framing = *self
            .framing_hold
            .get_or_insert_with(|| viewer_core::Framing::of(bounds));
        if self.needs_fit {
            self.camera.fit(bounds, &framing, aspect, 0.85);
            self.needs_fit = false;
        }
        // **Framing the selection keeps the whole scene's framing.** `Framing` is the run-wide
        // centre and span every projection here shares, and re-centring it on one object would
        // move everything else relative to it — the camera would appear to fit while the picture
        // silently changed what it was of. So only the focal length moves: the selection's box is
        // fitted *within* the scene's framing.
        if self.pending_frame {
            self.pending_frame = false;
            let want = self
                .selected
                .as_deref()
                .and_then(|p| tree_bounds(&self.checked, self.run.as_ref(), p));
            if let Some(b) = want {
                self.camera.fit(b, &framing, aspect, 0.85);
            }
        }

        // **The shaded pass.** Everything with a surface goes to the GPU with a depth buffer;
        // everything that is a *number* — labels, the colour bar, the readout — stays on egui on
        // top. The split is not stylistic: a caption has to be legible over whatever is behind it,
        // and geometry has to be occluded by whatever is in front of it, and one painter cannot do
        // both.
        if self.shaded {
            let key = self.geometry_key();
            if self.built != Some(key) {
                let (solid, lines, probes, notes) = self.batches(&framing);
                if let Ok(mut gpu) = self.gpu.lock() {
                    gpu.pending_solid = Some(solid);
                    gpu.pending_lines = Some(lines);
                }
                self.shaded_probes = probes;
                self.shaded_notes = notes;
                self.built = Some(key);
            }
            let clip = self.camera.matrix(aspect);
            let gpu = self.gpu.clone();
            if std::env::var_os("PANTOMETRY_VIEWPORT").is_some() {
                eprintln!("viewport: adding a callback for rect {rect:?}");
            }
            painter.add(egui::PaintCallback {
                rect,
                callback: std::sync::Arc::new(eframe::egui_glow::CallbackFn::new(
                    move |info, glow_painter| {
                        if std::env::var_os("PANTOMETRY_VIEWPORT").is_some() {
                            eprintln!("viewport: the callback ran");
                        }
                        let px = info.viewport_in_pixels();
                        if let Ok(mut shared) = gpu.lock() {
                            shared.paint(
                                glow_painter.gl(),
                                &clip,
                                [px.left_px, px.from_bottom_px, px.width_px, px.height_px],
                            );
                        }
                    },
                )),
            });
            // A driver that would not give a 3.3 core context loses the surfaces and keeps the
            // editor — and is told which of those happened, here, rather than shown an empty
            // rectangle that looks exactly like a scene with nothing in it.
            let failed = self.gpu.lock().ok().and_then(|g| g.error.clone());
            if let Some(why) = failed {
                self.shaded = false;
                self.status = format!("{why} — flat view instead");
            }
        }

        // **The projection returns a depth and this keeps it.** `Camera::project` computes
        // "distance from the eye, larger is further"; the closure here returned only a position,
        // so paint order was recovered afterwards by projecting each point a second time a
        // millimetre along world z and taking the reciprocal of the screen separation. That is
        // proportional to how far *off the view axis* a point is, not to how far away, and it
        // ordered 26% of a 6x6x6 splat lattice backwards — worst at the centre of the screen,
        // which is where the object is. The browser shell had always sorted by the real depth;
        // both call `editor_core::far_to_near` now, so there is one copy to be wrong.
        let project = |p: [f64; 3]| -> (egui::Pos2, f64) {
            let q = self.camera.project(p, &framing, aspect);
            (
                egui::pos2(
                    rect.center().x + (q.x as f32) * rect.width() * 0.5,
                    rect.center().y - (q.y as f32) * rect.height() * 0.5,
                ),
                q.depth,
            )
        };
        let to_screen = |p: [f64; 3]| -> egui::Pos2 { project(p).0 };

        // **The translate handles**, on top of everything because they are a control rather than
        // a picture — the same reason every caption and the colour bar are egui over the shaded
        // pass. Painted last, below, so nothing occludes what the pointer has to hit; the hit
        // test above uses the same origin and the same reach, so what is drawn is what is grabbed.
        //
        // In `Turn` the same list holds the **rings**: each is a closed polyline, and it comes
        // from `ring_points` — the one the hit test above walks — so what is drawn is what is
        // grabbed. That equality is the reason both go through the same function rather than
        // each computing a circle.
        let reach = framing.span * 0.35;
        let handles: Vec<(usize, Vec<egui::Pos2>)> = match self.selected_domain() {
            Some(name) => match self.handle_origin(&name) {
                Some(origin) => AXES
                    .iter()
                    .enumerate()
                    .map(|(i, (axis, _))| {
                        let world: Vec<[f64; 3]> = match self.gizmo {
                            Gizmo::Move => vec![
                                origin,
                                [
                                    origin[0] + axis[0] * reach,
                                    origin[1] + axis[1] * reach,
                                    origin[2] + axis[2] * reach,
                                ],
                            ],
                            Gizmo::Turn => {
                                let mut p = editor_core::edit::ring_points(origin, *axis, reach);
                                // Closed: the last segment back to the first is a segment a
                                // pointer can be over, and a ring with a gap in it is a ring
                                // that refuses in one place for no reason a person can see.
                                if let Some(first) = p.first().copied() {
                                    p.push(first);
                                }
                                p
                            }
                        };
                        (i, world.into_iter().map(to_screen).collect())
                    })
                    .collect(),
                None => Vec::new(),
            },
            None => Vec::new(),
        };

        // The scene's own geometry: every placed extent, wireframed. Drawn from the text, not
        // the run, so layout is visible while editing and before anything is computed.
        //
        // **Painted flat only when the shaded pass is off.** With surfaces on screen a wireframe
        // has to be behind them where it is behind them, and egui has no depth buffer to do that
        // with — so it moves into the line pass, which shares one. A box outline drawn over a solid
        // it is inside says the solid is not there.
        let wire = ui.visuals().weak_text_color();
        let highlight = egui::Color32::from_rgb(255, 190, 60);
        for placed in &self.checked.boxes {
            if !self.draws(&placed.name) {
                continue;
            }
            let picked = selected_name.as_deref() == Some(placed.name.as_str());
            let stroke = if picked {
                egui::Stroke::new(2.0_f32, highlight)
            } else {
                egui::Stroke::new(1.0_f32, wire)
            };
            if !self.shaded {
                for (a, b) in editor_core::EDGES {
                    painter.line_segment(
                        [to_screen(placed.corners[a]), to_screen(placed.corners[b])],
                        stroke,
                    );
                }
            }
            // The name, held inside the viewport. It hangs off a *projected* corner, so the
            // camera decides where it lands: in a 700-point window the only domain in scene 07
            // put its label at x=698 of 700 and the word was clipped to its first letters. Laid
            // out first and then clamped, rather than drawn and left where it fell.
            let colour = if picked {
                highlight
            } else {
                ui.visuals().text_color()
            };
            let galley = painter.layout_no_wrap(
                placed.name.clone(),
                egui::FontId::proportional(12.0),
                colour,
            );
            let at = to_screen(placed.corners[0]);
            let x =
                at.x.min(rect.right() - galley.rect.width())
                    .max(rect.left());
            painter.galley(egui::pos2(x, at.y - galley.rect.height()), galley, colour);
        }

        // How big any of this is. A viewport with no scale is a picture of a shape and not of a
        // part: the same wireframe serves a 40 mm die and a 4 m room, and until now nothing on
        // screen said which. Drawn from the framing's own span, so it is the model's metres and
        // not the window's pixels.
        scale_bar(
            &painter,
            rect,
            &project,
            &framing,
            ui.visuals().weak_text_color(),
        );

        // The run, by shape. Points are circles, paths are polylines, fields are splats, all
        // coloured by value on the run-wide scale — one scale across the run, never per frame,
        // for the reason viewer-core states; while a run streams, "the run" is the run so far,
        // and the colours settle when it does.
        // Answered before `self.run` is borrowed mutably below: `draws` reads `self`, and the
        // borrow checker is right that it cannot do so while `view` is out.
        let visible_names: Vec<String> = self
            .run
            .as_ref()
            .and_then(|v| v.run.frames.first())
            .map(|f| {
                f.panels
                    .iter()
                    .filter(|p| self.draws(p.name()))
                    .map(|p| p.name().to_string())
                    .collect()
            })
            .unwrap_or_default();
        // **Everything below is about a run, and there may not be one.** This used to be an
        // early `return`, which quietly took the rest of the function with it: the translate
        // handles are painted after this point, so on a scene that had not been run they were
        // computed, hit-tested and never drawn. Compiles, passes, does nothing.
        if let Some(view) = &mut self.run {
            if view.partial {
                painter.text(
                    egui::pos2(rect.left() + 8.0, rect.top() + 8.0),
                    egui::Align2::LEFT_TOP,
                    "streaming — a prefix of the run",
                    egui::FontId::proportional(12.0),
                    egui::Color32::from_rgb(230, 180, 60),
                );
            }
            let frames = view.run.frames.len();
            if frames > 1 {
                ui.scope_builder(egui::UiBuilder::new().max_rect(rect.shrink(8.0)), |ui| {
                    ui.with_layout(egui::Layout::bottom_up(egui::Align::LEFT), |ui| {
                        ui.horizontal(|ui| {
                            ui.add(
                                egui::Slider::new(&mut view.frame, 0..=frames - 1).text("frame"),
                            );
                            ui.label(format!(
                                "t = {} s",
                                editor_core::magnitude(view.run.frames[view.frame].t)
                            ));
                        });
                    });
                });
            }
            let frame = &view.run.frames[view.frame.min(frames - 1)];

            // The scale the colour bar will show. The first panel that has one: a viewport holding
            // two quantities cannot draw one bar for both, and drawing the first with its name on it
            // beats drawing none.
            let mut legend: Option<(f64, f64, String, bool)> = None;

            for panel in &frame.panels {
                if !visible_names.iter().any(|n| n == panel.name()) {
                    continue;
                }
                let scale = view.run.scale_of(panel.name());
                if legend.is_none() {
                    if let Some((lo, hi)) = scale {
                        legend = Some((
                            lo,
                            hi,
                            format!("{} / {}", panel.name(), panel.unit()),
                            editor_core::scale_is_signed(scale),
                        ));
                    }
                }
                match panel {
                    viewer_core::Panel::Points {
                        name,
                        unit,
                        positions,
                        values,
                        ..
                    } => {
                        // Far to near, so a body in front covers one behind rather than whichever
                        // happened to be last in the array.
                        let pts = panel.placed_positions();
                        let _ = positions;
                        let depths: Vec<f64> = pts.iter().map(|p| project(*p).1).collect();
                        for i in editor_core::far_to_near(&depths) {
                            let at = to_screen(pts[i]);
                            if !self.shaded {
                                painter.circle_filled(at, 3.5, shade(values[i], scale));
                            }
                            probe.offer(
                                pointer,
                                at,
                                6.0,
                                depths[i],
                                &format!("/run/{name}"),
                                || {
                                    format!(
                                        "{name} body {i}   {} {unit}",
                                        editor_core::magnitude(values[i])
                                    )
                                },
                            );
                        }
                    }
                    viewer_core::Panel::Paths {
                        starts,
                        vertices,
                        values,
                        ..
                    } => {
                        for (r, value) in values.iter().enumerate() {
                            let lo = starts[r] as usize;
                            let hi = starts
                                .get(r + 1)
                                .map_or(vertices.len() / 3, |s| *s as usize);
                            if self.shaded {
                                continue;
                            }
                            let placed = panel.placed_vertices();
                            for w in lo..hi.saturating_sub(1) {
                                let (a, b) = (placed[w], placed[w + 1]);
                                painter.line_segment(
                                    [to_screen(a), to_screen(b)],
                                    egui::Stroke::new(1.5_f32, shade(*value, scale)),
                                );
                            }
                        }
                    }
                    viewer_core::Panel::Field {
                        name,
                        unit,
                        nx,
                        ny,
                        nz,
                        values,
                        ..
                    } => {
                        // **The run's own extent first, the scene's placed box second.** A field
                        // carries the box it was sampled over now, so a run opened without the file
                        // that produced it still draws in the right place and at the right size. The
                        // placed box is the fallback for a run written before the format carried it,
                        // and a field with neither has nowhere to be drawn — which is said rather
                        // than shown as an absence.
                        let from_run =
                            panel
                                .placed_corners()
                                .map(|corners| editor_core::PlacedBox {
                                    name: name.clone(),
                                    corners,
                                });
                        let placed = from_run
                            .as_ref()
                            .or_else(|| self.checked.boxes.iter().find(|b| &b.name == name));
                        let Some(b) = placed else {
                            continue;
                        };
                        let note = if self.shaded {
                            self.shaded_notes
                                .iter()
                                .find(|(n, _)| n == name)
                                .map(|(_, note)| *note)
                        } else {
                            draw_field(
                                &painter,
                                &project,
                                b,
                                (*nx, *ny, *nz),
                                values,
                                unit,
                                scale,
                                pointer,
                                &mut probe,
                                name,
                            )
                        };
                        if let Some(note) = note {
                            painter.text(
                                to_screen(b.corners[0]),
                                egui::Align2::LEFT_TOP,
                                note,
                                egui::FontId::proportional(10.0),
                                ui.visuals().weak_text_color(),
                            );
                        }
                    }
                }
            }

            // The colour bar. Without it the viewport says *more* and *less* and never *how much* —
            // which is the whole difference between a picture and a reading, and the report grew one
            // for every view for the same reason.
            if let Some((lo, hi, label, signed)) = legend {
                colour_bar(&painter, rect, lo, hi, &label, signed, ui.visuals());
            }

            // The frame's readings, top-right: the numbers for everything that has no picture.
            //
            // `editor_core::magnitude`, not `{:.4}` — a cavity holding 3.19e-10 J printed `0.0000`
            // here beside a field the same run reported at 921 V/m.
            let readings = &frame.readings;
            if !readings.is_empty() {
                let mut y = rect.top() + 8.0;
                for r in readings {
                    painter.text(
                        egui::pos2(rect.right() - 8.0, y),
                        egui::Align2::RIGHT_TOP,
                        format!(
                            "{} {} {} {}",
                            r.domain,
                            r.label,
                            editor_core::magnitude(r.value),
                            r.unit
                        ),
                        egui::FontId::monospace(11.0),
                        ui.visuals().text_color(),
                    );
                    y += 14.0;
                }
            }

            // The shaded view's readout, from the cache rather than from a fresh mesh.
            if self.shaded {
                for t in &self.shaded_probes {
                    let (at, depth) = project(t.at);
                    probe.offer(pointer, at, 6.0, depth, &t.path, || t.label.clone());
                }
            }

            // What the shaded pass actually drew. Every DCC viewport has this readout, and it is also
            // the only thing that distinguishes a pass drawing nothing from a pass never asked to run.
            if self.shaded {
                if let Ok(gpu) = self.gpu.lock() {
                    let (tris, lines, paints) = gpu.drawn;
                    painter.text(
                        egui::pos2(rect.left() + 8.0, rect.bottom() - 8.0),
                        egui::Align2::LEFT_BOTTOM,
                        format!("{tris} triangles, {lines} lines, {paints} paints"),
                        egui::FontId::monospace(10.0),
                        ui.visuals().weak_text_color(),
                    );
                }
            }
        }

        // The handles, over everything. Drawn from the positions computed before the geometry so
        // the hit test and the picture cannot drift apart.
        for (i, points) in &handles {
            let held = self.grabbed == Some(*i);
            let colour = AXES[*i].1;
            let stroke = egui::Stroke::new(if held { 3.5_f32 } else { 2.0_f32 }, colour);
            for pair in points.windows(2) {
                painter.line_segment([pair[0], pair[1]], stroke);
            }
            // An arrow gets a head; a ring does not, because a ring has no end and a blob on one
            // would read as the handle rather than as the whole circle being it.
            if self.gizmo == Gizmo::Move {
                if let Some(tip) = points.last() {
                    painter.circle_filled(*tip, if held { 6.0_f32 } else { 4.5_f32 }, colour);
                }
            }
        }

        // And what the cursor is over, if anything. A click on it selects.
        //
        // **Not while a handle is held.** A drag that ends on empty space reads as a click to
        // egui, so letting go of a handle over nothing would clear the selection and take the
        // handles away from under the pointer that was just using them.
        //
        // **And not at all without a run.** Everything that offers itself to the probe is a run
        // panel, so with no run `probe.path` is always `None` and this would clear the selection
        // on every click in the viewport — including the click that just selected something. It
        // was unreachable before the block above stopped being an early `return`, which is the
        // sort of thing removing a `return` quietly switches on.
        if response.clicked() && self.grabbed.is_none() && self.run.is_some() {
            // A click on nothing clears the selection, which is what every outliner does and is
            // the only way to get back to "no selection" without a keyboard.
            self.selected = probe.path.clone();
        }
        probe.draw(&painter, ui.visuals());
    }
}

/// The box a tree path refers to, for framing it.
///
/// Rebuilt from the same two sources the tree is, rather than carried on the selection: a
/// selection is a path and paths outlive the objects they name, so resolving late is what makes
/// a stale selection a missing box instead of a wrong one.
fn tree_bounds(
    checked: &editor_core::Checked,
    run: Option<&RunView>,
    path: &str,
) -> Option<[f64; 6]> {
    let frame = run.map_or(0, |v| v.frame);
    let tree = editor_core::tree(checked, run.map(|v| &v.run), frame);
    tree.find(path).and_then(|n| n.bounds)
}

/// The nearest thing under the cursor, and what to say about it.
///
/// A viewport you cannot point at makes a reader export a CSV to answer "what is that". Nearest
/// **to the eye** rather than nearest on screen among things within reach, so the thing named is
/// the thing they can see: two bodies overlapping in projection are two different answers and
/// only one of them is visible.
#[derive(Default)]
struct Probe {
    best: Option<(f64, egui::Pos2, String)>,
    /// The tree path of whatever the readout named, so a click selects the thing under the
    /// cursor. The outliner and the viewport are two views of one selection; a viewport you can
    /// read but not select from makes the reader find the row by hand.
    path: Option<String>,
}

impl Probe {
    fn offer(
        &mut self,
        pointer: Option<egui::Pos2>,
        at: egui::Pos2,
        reach: f32,
        depth: f64,
        path: &str,
        text: impl FnOnce() -> String,
    ) {
        let Some(p) = pointer else { return };
        if (p - at).length() > reach {
            return;
        }
        if self.best.as_ref().is_none_or(|(d, _, _)| depth < *d) {
            self.best = Some((depth, at, text()));
            self.path = Some(path.to_string());
        }
    }

    fn draw(&self, painter: &egui::Painter, visuals: &egui::Visuals) {
        let Some((_, at, text)) = &self.best else {
            return;
        };
        let font = egui::FontId::monospace(11.0);
        let galley = painter.layout_no_wrap(text.clone(), font, visuals.strong_text_color());
        let pad = egui::vec2(6.0, 3.0);
        let origin = *at + egui::vec2(10.0, -galley.size().y - 10.0);
        let box_rect = egui::Rect::from_min_size(origin, galley.size() + pad * 2.0);
        painter.rect_filled(box_rect, 3.0, visuals.extreme_bg_color);
        painter.rect_stroke(
            box_rect,
            3.0,
            egui::Stroke::new(1.0_f32, visuals.weak_text_color()),
            egui::StrokeKind::Inside,
        );
        painter.galley(origin + pad, galley, visuals.strong_text_color());
        painter.circle_stroke(
            *at,
            5.0,
            egui::Stroke::new(1.5_f32, visuals.strong_text_color()),
        );
    }
}

/// A bar of the scale in use, with the numbers a reader needs to invert it.
///
/// The ends are `editor_core::bar_value`'s, so a diverging scale is labelled with the deflection
/// it actually covers rather than with the data's range — the two differ whenever the range is
/// not symmetric about zero, and labelling a symmetric scale with an asymmetric range is a lie
/// about where the neutral colour sits.
fn colour_bar(
    painter: &egui::Painter,
    rect: egui::Rect,
    lo: f64,
    hi: f64,
    label: &str,
    signed: bool,
    visuals: &egui::Visuals,
) {
    let scale = Some((lo, hi));
    let w = (rect.width() * 0.28).clamp(120.0, 280.0);
    let x0 = rect.right() - 16.0 - w;
    let y = rect.bottom() - 34.0;
    let steps = w as usize;
    for i in 0..steps {
        let u = i as f64 / (steps - 1).max(1) as f64;
        let [r, g, b] = editor_core::value_colour(editor_core::bar_value(u, scale), scale);
        painter.rect_filled(
            egui::Rect::from_min_size(egui::pos2(x0 + i as f32, y), egui::vec2(1.0, 9.0)),
            0.0,
            egui::Color32::from_rgb(r, g, b),
        );
    }
    let font = egui::FontId::monospace(10.0);
    for (u, align) in [
        (0.0, egui::Align2::LEFT_TOP),
        (0.5, egui::Align2::CENTER_TOP),
        (1.0, egui::Align2::RIGHT_TOP),
    ] {
        painter.text(
            egui::pos2(x0 + (u as f32) * w, y + 11.0),
            align,
            editor_core::magnitude(editor_core::bar_value(u, scale)),
            font.clone(),
            visuals.weak_text_color(),
        );
    }
    painter.text(
        egui::pos2(x0 - 8.0, y + 4.0),
        egui::Align2::RIGHT_CENTER,
        if signed {
            format!("{label}  (0 at the middle)")
        } else {
            label.to_string()
        },
        egui::FontId::monospace(10.0),
        visuals.weak_text_color(),
    );
}

/// A bar of known length in model metres, so the viewport says how big the thing is.
///
/// Measured by projecting two points a known distance apart at the centre of the framing and
/// reading how far they land apart on screen, then rounding that length down a 1-2-5 ladder —
/// the same reason the report's ticks use one, which is that a bar labelled 0.3714 m is a bar
/// nobody reads.
fn scale_bar(
    painter: &egui::Painter,
    rect: egui::Rect,
    project: &impl Fn([f64; 3]) -> (egui::Pos2, f64),
    framing: &viewer_core::Framing,
    colour: egui::Color32,
) {
    let c = framing.centre;
    let probe = framing.span.max(f64::MIN_POSITIVE);
    let a = project(c).0;
    let b = project([c[0] + probe, c[1], c[2]]).0;
    let px_per_metre = ((b - a).length() as f64) / probe;
    if !(px_per_metre.is_finite() && px_per_metre > 0.0) {
        return;
    }
    // Aim for about a fifth of the window, then round down to 1, 2 or 5 times a power of ten.
    let want = (rect.width() as f64 * 0.2) / px_per_metre;
    let mag = 10f64.powf(want.log10().floor());
    let n = want / mag;
    let length = mag
        * if n >= 5.0 {
            5.0
        } else if n >= 2.0 {
            2.0
        } else {
            1.0
        };
    let px = (length * px_per_metre) as f32;
    if !(px.is_finite() && (12.0..rect.width() * 0.8).contains(&px)) {
        return;
    }
    let y = rect.bottom() - 14.0;
    let x0 = rect.left() + 16.0;
    let stroke = egui::Stroke::new(1.5_f32, colour);
    painter.line_segment([egui::pos2(x0, y), egui::pos2(x0 + px, y)], stroke);
    for x in [x0, x0 + px] {
        painter.line_segment([egui::pos2(x, y - 4.0), egui::pos2(x, y + 4.0)], stroke);
    }
    painter.text(
        egui::pos2(x0 + px * 0.5, y - 6.0),
        egui::Align2::CENTER_BOTTOM,
        metres(length),
        egui::FontId::monospace(10.0),
        colour,
    );
}

/// A length in metres, on the unit an engineer would say it in.
fn metres(m: f64) -> String {
    let a = m.abs();
    if a >= 1e4 {
        format!("{m:.2e} m")
    } else if a >= 1.0 {
        format!("{} m", trim(format!("{m:.3}")))
    } else if a >= 1e-3 {
        format!("{} mm", trim(format!("{:.3}", m * 1e3)))
    } else {
        format!("{} um", trim(format!("{:.3}", m * 1e6)))
    }
}

/// Trailing zeros only count as trailing **after a decimal point**. Stripping them without that
/// check turned 500 into 5 in the report, and captioned a 500 um channel as a 5 um one.
fn trim(s: String) -> String {
    if !s.contains('.') {
        return s;
    }
    s.trim_end_matches('0').trim_end_matches('.').to_string()
}

#[allow(clippy::too_many_arguments)]
fn draw_field(
    painter: &egui::Painter,
    project: &impl Fn([f64; 3]) -> (egui::Pos2, f64),
    placed: &editor_core::PlacedBox,
    counts: (usize, usize, usize),
    values: &[f64],
    unit: &str,
    scale: Option<(f64, f64)>,
    pointer: Option<egui::Pos2>,
    probe: &mut Probe,
    name: &str,
) -> Option<&'static str> {
    let out = editor_core::field_splats(&placed.corners, counts, values, unit, scale);
    if out.splats.is_empty() {
        return Some(out.note);
    }

    // One splat's screen size: the box's own projected span divided by its grid, so a coarse
    // field draws fat cells and a fine one draws small ones instead of both drawing dots.
    let span_px = {
        let a = project(placed.corners[0]).0;
        let b = project(placed.corners[7]).0;
        (a - b).length()
    };
    let (nx, ny, nz) = counts;
    let radius =
        (span_px / (nx.max(ny).max(nz) as f32) * 0.7 / out.stride.max(1) as f32).clamp(1.0, 40.0);

    // Painter's algorithm, far to near, which is what makes a translucent volume composite in
    // the right order — on the depth the projection returned, not on a stand-in for it.
    let placed_at: Vec<(egui::Pos2, f64)> = out.splats.iter().map(|s| project(s.at)).collect();
    let depths: Vec<f64> = placed_at.iter().map(|(_, d)| *d).collect();
    let stride = out.stride.max(1);
    for i in editor_core::far_to_near(&depths) {
        let s = &out.splats[i];
        let (at, depth) = placed_at[i];
        let c = egui::Color32::from_rgba_unmultiplied(s.rgba[0], s.rgba[1], s.rgba[2], s.rgba[3]);
        painter.circle_filled(at, radius, c);
        // The grid index this splat came from, recovered from its place in the strided walk, so
        // the readout names a cell a reader can find in the CSV rather than a serial number.
        let per_row = nx.div_ceil(stride);
        let per_plane = per_row * ny.div_ceil(stride);
        let (gi, gj, gk) = (
            (i % per_row) * stride,
            (i / per_row % ny.div_ceil(stride)) * stride,
            (i / per_plane) * stride,
        );
        let v = values
            .get(gi + nx * (gj + ny * gk))
            .copied()
            .unwrap_or(f64::NAN);
        probe.offer(
            pointer,
            at,
            radius.max(4.0),
            depth,
            &format!("/run/{name}"),
            || {
                format!(
                    "{name} [{gi},{gj},{gk}]   {} {unit}",
                    editor_core::magnitude(v)
                )
            },
        );
    }
    Some(out.note)
}

/// A value on the run-wide scale, as a colour.
///
/// `editor_core::value_colour`, which is `pantometry::view::ramp` — the same scale the HTML report
/// draws with, built in CIE LCh with lightness linear in the value. What was here was a straight
/// line in sRGB from blue to red, duplicated in the browser shell, covering 16.6 L\* against the
/// library scale's 74: a scale that a greyscale print, a projector and a colour-vision deficiency
/// all reduce to one flat block.
fn shade(value: f64, scale: Option<(f64, f64)>) -> egui::Color32 {
    let [r, g, b] = editor_core::value_colour(value, scale);
    egui::Color32::from_rgb(r, g, b)
}

/// The file's modified time, or `None` for a path that does not resolve to one.
fn mtime_of(path: &str) -> Option<SystemTime> {
    std::fs::metadata(path).and_then(|m| m.modified()).ok()
}

/// Union of an optional box with another box.
fn union(a: Option<[f64; 6]>, b: [f64; 6]) -> [f64; 6] {
    match a {
        None => b,
        Some(a) => [
            a[0].min(b[0]),
            a[1].min(b[1]),
            a[2].min(b[2]),
            a[3].max(b[3]),
            a[4].max(b[4]),
            a[5].max(b[5]),
        ],
    }
}

/// The scene the editor opens on with no file: the same built-in room `pantometry-world` runs
/// with no arguments, so the two front ends agree about where "hello" is.
fn default_scene() -> String {
    String::from(
        r#"{
  "title": "a small room ringing in its (1,1) mode",
  "schedule": "multirate",
  "duration_s": 0.02,
  "frames": 11,
  "conservation_tolerance": 1e-6,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 4.4, "height_m": 3.1,
      "cells_across": 61,
      "release": { "as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0 } }
  ]
}
"#,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **A pointer is measured against the handle's line, not against its ends.**
    ///
    /// The distinction is the whole function. Measuring to the nearer endpoint would make the
    /// middle of a long handle unreachable, and measuring to the infinite line would let a click
    /// far past the tip take hold of it — a handle that grabs from off-screen is worse than one
    /// that does not grab at all, because the thing it moves is somewhere the pointer is not.
    #[test]
    fn a_handle_is_grabbed_by_its_length() {
        let a = egui::pos2(100.0, 100.0);
        let b = egui::pos2(200.0, 100.0);

        // Beside the middle: the perpendicular distance, exactly.
        assert!((near_segment(egui::pos2(150.0, 107.0), a, b) - 7.0).abs() < 1e-4);
        // On it: nothing.
        assert!(near_segment(egui::pos2(150.0, 100.0), a, b) < 1e-6);
        // Past the tip: the distance to the tip, not to the line it lies on.
        assert!((near_segment(egui::pos2(260.0, 100.0), a, b) - 60.0).abs() < 1e-4);
        assert!((near_segment(egui::pos2(40.0, 100.0), a, b) - 60.0).abs() < 1e-4);
        // Diagonally past the tip: the hypotenuse, which the infinite-line answer would miss.
        let corner = near_segment(egui::pos2(203.0, 104.0), a, b);
        assert!((corner - 5.0).abs() < 1e-4, "got {corner}");
    }

    /// **A handle seen exactly end-on is a point, and does not divide by zero.**
    ///
    /// It happens whenever an axis points at the camera, which is a rotation away at any moment.
    #[test]
    fn a_handle_with_no_length_is_still_measurable() {
        let p = egui::pos2(50.0, 50.0);
        assert!(near_segment(egui::pos2(53.0, 54.0), p, p) - 5.0 < 1e-4);
        assert!(near_segment(p, p, p) < 1e-6);
    }
}
