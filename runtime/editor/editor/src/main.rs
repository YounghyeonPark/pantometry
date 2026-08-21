//! Edit a pantometry scene beside a 3D view of it — and leave it open while a script does the
//! editing.
//!
//! ```text
//! cargo run --release                 # opens on the built-in room
//! cargo run --release -- scene.json   # opens on a file
//! ```
//!
//! The left pane is the scene's JSON, checked as you type with the same two steps
//! `pantometry-world --check` runs; a parse error is shown with its `line:column`. The viewport
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

use editor_core::OnDisk;
use eframe::egui;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{mpsc, Arc};
use std::time::{Duration, Instant, SystemTime};

fn main() -> eframe::Result {
    let path = std::env::args().nth(1);
    eframe::run_native(
        "pantometry editor",
        eframe::NativeOptions::default(),
        Box::new(|_cc| Ok(Box::new(App::new(path)))),
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
    /// Whether the outliner and the inspector are on screen at all.
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
        let checked = editor_core::check(&text, &OnDisk);
        let known_mtime = mtime_of(&path);
        App {
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
            show_outliner: true,
            show_inspector: true,
            show_text: true,
        }
    }

    fn recheck(&mut self) {
        self.checked = editor_core::check(&self.text, &OnDisk);
        self.run = None;
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
        self.stop = Arc::new(AtomicBool::new(false));
        let stop = self.stop.clone();
        self.spawn("running", move |tx| {
            let end = editor_core::run_streaming(&text, &OnDisk, &stop, |json| {
                let _ = tx.send(Job::Frames(json));
            });
            let _ = tx.send(Job::RunEnded(end));
        });
    }

    fn start_verify(&mut self) {
        let text = self.text.clone();
        let deep = self.deep;
        self.spawn("verifying", move |tx| {
            let _ = tx.send(Job::Verified(editor_core::verify(&text, deep, &OnDisk)));
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
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
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
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Load").clicked() {
                        self.load_from_disk();
                        ui.close_menu();
                    }
                    if ui.button("Save").clicked() {
                        self.save_to_disk();
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui.button("Quit").clicked() {
                        ui.ctx().send_viewport_cmd(egui::ViewportCommand::Close);
                    }
                });
                ui.menu_button("View", |ui| {
                    ui.checkbox(&mut self.show_outliner, "Outliner");
                    ui.checkbox(&mut self.show_inspector, "Inspector");
                    ui.checkbox(&mut self.show_text, "Scene text");
                    ui.separator();
                    ui.checkbox(&mut self.solo, "Solo the selection");
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
            ui.horizontal(|ui| {
                ui.label("file");
                ui.add(egui::TextEdit::singleline(&mut self.path).desired_width(220.0));
                if ui.button("load").clicked() {
                    self.load_from_disk();
                }
                if ui.button("save").clicked() {
                    self.save_to_disk();
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
                ui.checkbox(&mut self.watch, "watch file");
                ui.checkbox(&mut self.auto_run, "run on change");
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
            match self.checked.error.as_deref() {
                Some(e) => {
                    ui.colored_label(egui::Color32::from_rgb(220, 80, 80), e);
                }
                None => {
                    ui.label(&self.status);
                }
            }
        });

        // Rebuilt every paint, from the checked scene and whatever the run has produced so far.
        // Cheap — a few hundred rows at most — and the alternative is a cache that goes stale
        // exactly when a streaming run is adding to it.
        let frame_at = self.run.as_ref().map_or(0, |v| v.frame);
        let tree = editor_core::tree(&self.checked, self.run.as_ref().map(|v| &v.run), frame_at);

        if self.show_outliner {
            egui::SidePanel::left("outliner")
                .resizable(true)
                .default_width(260.0)
                .show(ctx, |ui| self.outliner(ui, &tree));
        }
        if self.show_inspector {
            egui::SidePanel::right("inspector")
                .resizable(true)
                .default_width(320.0)
                .show(ctx, |ui| self.inspector(ui, &tree));
        }

        if self.show_text {
            egui::SidePanel::left("text")
                .resizable(true)
                .default_width(430.0)
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
    fn load_from_disk(&mut self) {
        match std::fs::read_to_string(&self.path) {
            Ok(t) => {
                self.text = t;
                self.dirty = false;
                self.known_mtime = mtime_of(&self.path);
                self.recheck();
                self.needs_fit = true;
                self.status = format!("loaded {}", self.path);
            }
            Err(e) => self.status = format!("{}: {e}", self.path),
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
                        if ui.selectable_label(selected, text).clicked() {
                            self.selected = Some(node.path.clone());
                        }
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            ui.label(
                                egui::RichText::new(node.kind.label())
                                    .weak()
                                    .monospace()
                                    .size(10.0),
                            );
                        });
                    });
                }
            });
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

    /// The 3D viewport: wireframe extents from the text, run panels over them by shape.
    fn viewport(&mut self, ui: &mut egui::Ui) {
        let (response, painter) =
            ui.allocate_painter(ui.available_size(), egui::Sense::click_and_drag());
        let rect = response.rect;
        let aspect = (rect.width() / rect.height().max(1.0)) as f64;

        if response.dragged() {
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

        // The world this paint frames: scene geometry, widened by whatever the run reached —
        // an orbit's bodies live far outside any placed extent.
        let mut bounds = self.checked.bounds;
        if let Some(view) = &self.run {
            if let Some(frame) = view.run.frames.first() {
                for panel in &frame.panels {
                    bounds = Some(union(bounds, panel.bounds()));
                }
            }
        }
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
        let framing = viewer_core::Framing::of(bounds);
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

        // The scene's own geometry: every placed extent, wireframed. Drawn from the text, not
        // the run, so layout is visible while editing and before anything is computed.
        let wire = ui.visuals().weak_text_color();
        let highlight = egui::Color32::from_rgb(255, 190, 60);
        for placed in &self.checked.boxes {
            if !self.draws(&placed.name) {
                continue;
            }
            let picked = selected_name.as_deref() == Some(placed.name.as_str());
            let stroke = if picked {
                egui::Stroke::new(2.0, highlight)
            } else {
                egui::Stroke::new(1.0, wire)
            };
            for (a, b) in editor_core::EDGES {
                painter.line_segment(
                    [to_screen(placed.corners[a]), to_screen(placed.corners[b])],
                    stroke,
                );
            }
            painter.text(
                to_screen(placed.corners[0]),
                egui::Align2::LEFT_BOTTOM,
                &placed.name,
                egui::FontId::proportional(12.0),
                if picked {
                    highlight
                } else {
                    ui.visuals().text_color()
                },
            );
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
        let Some(view) = &mut self.run else {
            return;
        };
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
                        ui.add(egui::Slider::new(&mut view.frame, 0..=frames - 1).text("frame"));
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
                    let pts: Vec<[f64; 3]> = (0..values.len())
                        .map(|i| [positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]])
                        .collect();
                    let depths: Vec<f64> = pts.iter().map(|p| project(*p).1).collect();
                    for i in editor_core::far_to_near(&depths) {
                        let at = to_screen(pts[i]);
                        painter.circle_filled(at, 3.5, shade(values[i], scale));
                        probe.offer(pointer, at, 6.0, depths[i], &format!("/run/{name}"), || {
                            format!(
                                "{name} body {i}   {} {unit}",
                                editor_core::magnitude(values[i])
                            )
                        });
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
                        for w in lo..hi.saturating_sub(1) {
                            let a = [vertices[3 * w], vertices[3 * w + 1], vertices[3 * w + 2]];
                            let b = [
                                vertices[3 * w + 3],
                                vertices[3 * w + 4],
                                vertices[3 * w + 5],
                            ];
                            painter.line_segment(
                                [to_screen(a), to_screen(b)],
                                egui::Stroke::new(1.5, shade(*value, scale)),
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
                    let from_run = panel.extent_m().map(|e| editor_core::PlacedBox {
                        name: name.clone(),
                        corners: editor_core::corners_of(e),
                    });
                    let placed = from_run
                        .as_ref()
                        .or_else(|| self.checked.boxes.iter().find(|b| &b.name == name));
                    let Some(b) = placed else {
                        continue;
                    };
                    if let Some(note) = draw_field(
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
                    ) {
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

        // And what the cursor is over, if anything. A click on it selects.
        if response.clicked() {
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
            egui::Stroke::new(1.0, visuals.weak_text_color()),
            egui::StrokeKind::Inside,
        );
        painter.galley(origin + pad, galley, visuals.strong_text_color());
        painter.circle_stroke(
            *at,
            5.0,
            egui::Stroke::new(1.5, visuals.strong_text_color()),
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
    let stroke = egui::Stroke::new(1.5, colour);
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
