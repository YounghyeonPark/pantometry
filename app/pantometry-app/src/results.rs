//! **Results, under the view.**
//!
//! A run on screen is a picture and a slider. What a person running a simulation tool does next
//! is ask it for numbers: the hottest place and how hot, the mean over a part, a profile along a
//! line, a cut through the middle, and how a point went over the run. This panel asks
//! [`editor_core::post`] for each of them and draws what comes back; it computes nothing itself,
//! so the browser editor can draw the same numbers from the same calls.
//!
//! Five tabs, one domain at a time:
//!
//! - **Region** — the largest and smallest value and where they are, the mean by volume, the
//!   integral and the volume, over the domain or a box of it, and the three of them over the run.
//! - **Slice** — the field on a plane across one axis, in the view's own colours.
//! - **Line** — the field along a straight line, against distance.
//! - **Point** — the value at a point now, and at every frame.
//! - **Readings** — every reading against time: the chart a run with nothing to place gets in the
//!   viewport, here for every run.
//!
//! Positions are in the domain's own frame, the one its grid is stated in, and in the length unit
//! its size is said in. What the tab is reading is drawn in the viewport as well — the plane, the
//! line, the point, the box — so a number is never far from the place it is about.

use editor_core::chart::{self, Across, Line};
use editor_core::post::{self, Grid, Stats, AXIS_NAMES};
use eframe::egui;
use std::hash::{Hash, Hasher};

/// Which question the panel is asking.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tab {
    /// Statistics over the domain or a box of it.
    Region,
    /// A plane across one axis.
    Slice,
    /// A straight line.
    Line,
    /// One point, over time.
    Point,
    /// Every reading, over time.
    Readings,
}

impl Tab {
    /// In the order the panel offers them.
    pub const ALL: [Tab; 5] = [
        Tab::Region,
        Tab::Slice,
        Tab::Line,
        Tab::Point,
        Tab::Readings,
    ];

    /// What the tab is called.
    pub fn label(self) -> &'static str {
        match self {
            Tab::Region => "Region",
            Tab::Slice => "Slice",
            Tab::Line => "Line",
            Tab::Point => "Point",
            Tab::Readings => "Readings",
        }
    }

    /// The tab a word names, ignoring case, for `--ui-dump --results`.
    pub fn named(word: &str) -> Option<Tab> {
        Tab::ALL
            .into_iter()
            .find(|t| t.label().eq_ignore_ascii_case(word))
    }
}

/// What the panel is asking, and where.
#[derive(Default)]
pub struct Results {
    /// The tab on screen, or `None` until a run arrives: then the region when the run has a field,
    /// and the readings when it has none.
    pub tab: Option<Tab>,
    /// The domain asked about.
    pub domain: Option<String>,
    /// The domain and the extent the positions below were set for, to the bit. A different
    /// domain or a grid of a different size puts them back in its middle; a rerun of the same
    /// grid leaves a probe where it was put.
    set_for: Option<(String, [u64; 6])>,
    /// The slice's axis, and where along it, in metres.
    pub axis: usize,
    /// Where along the axis the plane is.
    pub at: f64,
    /// The line's two ends, in metres.
    pub a: [f64; 3],
    /// The far end.
    pub b: [f64; 3],
    /// The probe, in metres.
    pub point: [f64; 3],
    /// Whether the region is the whole domain rather than the box below.
    pub whole: bool,
    /// The box's near corner, in metres.
    pub lo: [f64; 3],
    /// Its far corner.
    pub hi: [f64; 3],
    /// The reading the readings tab charts, by its place in the frame's list.
    ///
    /// **By place, not by name.** Two readings can share a domain and a label — `verify` flags
    /// exactly that — and keyed by the pair both rows lit and only the first could be charted.
    pub reading: Option<usize>,
    /// The outliner's selection when it was last applied, so it is applied when it changes and
    /// not on every frame.
    prefer_was: Option<String>,
    /// The slice's image, and what it was made from.
    texture: Option<(u64, egui::TextureHandle)>,
    /// The region's statistics in every frame, and what they were computed from. A pass over
    /// every sample of every frame is not a thing to do sixty times a second.
    over_time: Option<(u64, OverTime)>,
}

/// A region's statistics in every frame: `(t in seconds, the statistics or None)`, as
/// [`post::region_history`] gives them.
type OverTime = Vec<(f64, Option<Stats>)>;

impl Results {
    /// Put the positions in `g`'s middle, unless they were set for this domain and this extent.
    fn fit_to(&mut self, g: &Grid) {
        let key = (g.name.to_string(), g.extent.map(f64::to_bits));
        if self.set_for.as_ref() == Some(&key) {
            return;
        }
        self.set_for = Some(key);
        let c = g.centre();
        let span = |a: usize| g.extent[a + 3] - g.extent[a];
        // The longest axis, the first of equals: a line through the middle along it shows the
        // most of the part.
        let long = (1..3).fold(0, |best, a| if span(a) > span(best) { a } else { best });
        self.axis = 2;
        self.at = c[2];
        self.a = c;
        self.a[long] = g.extent[long];
        self.b = c;
        self.b[long] = g.extent[long + 3];
        self.point = c;
        self.whole = true;
        self.lo = [g.extent[0], g.extent[1], g.extent[2]];
        self.hi = [g.extent[3], g.extent[4], g.extent[5]];
    }

    /// The box the region tab reads, or `None` for the whole domain.
    fn within(&self) -> Option<[f64; 6]> {
        (!self.whole).then_some([
            self.lo[0], self.lo[1], self.lo[2], self.hi[0], self.hi[1], self.hi[2],
        ])
    }
}

/// The length unit a domain's positions are shown in: the one its largest side is said in.
fn length_of(g: &Grid) -> (&'static str, f64) {
    let span = (0..3)
        .map(|a| g.extent[a + 3] - g.extent[a])
        .fold(0.0, f64::max);
    chart::length_unit(span)
}

/// A position along `axis` of `g`, dragged in the display unit and held to the box.
fn coordinate(ui: &mut egui::Ui, value: &mut f64, g: &Grid, axis: usize, per: f64) -> bool {
    let (lo, hi) = (g.extent[axis] / per, g.extent[axis + 3] / per);
    let mut shown = *value / per;
    let changed = ui
        .add(
            egui::DragValue::new(&mut shown)
                .speed(((hi - lo) / 200.0).max(1e-9))
                .range(lo..=hi)
                .custom_formatter(|v, _| short(v))
                .custom_parser(|s| s.trim().parse::<f64>().ok()),
        )
        .changed();
    // Written back only when moved: `shown * per` is not `value` to the bit, and a probe that
    // drifted a rounding a frame would stop being on the sample it was put on.
    if changed {
        *value = shown * per;
    }
    changed
}

/// A row of controls that wraps when the panel is narrow, with every label laid out whole.
///
/// **Whole, because a wrapping row lays a wrapping label out from the row's left edge** and starts
/// its text at the cursor, so each label's box covers every widget before it on the line. Nothing
/// is drawn over anything, and `--ui-dump` cannot tell: the first frames of this panel reported
/// sixty overlaps on the line tab, every one of them a label's box over the controls to its left.
fn row(ui: &mut egui::Ui, add: impl FnOnce(&mut egui::Ui)) {
    ui.horizontal_wrapped(|ui| {
        ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
        add(ui);
    });
}

/// A table and a chart beside it, in all the height the panel has left.
///
/// **Not `ui.horizontal`, which is one line tall until something makes it taller.** The region's
/// eight rows of statistics made it tall enough for the chart beside them; the readings tab's
/// table of three made it 64 points, under the chart's floor, and the chart was not drawn.
fn side_by_side(ui: &mut egui::Ui, add: impl FnOnce(&mut egui::Ui)) {
    ui.allocate_ui_with_layout(
        ui.available_size(),
        egui::Layout::left_to_right(egui::Align::Min),
        add,
    );
}

/// Three positions in a row, each labelled with its axis.
fn triple(ui: &mut egui::Ui, label: &str, p: &mut [f64; 3], g: &Grid, per: f64) {
    ui.label(label);
    for (a, name) in AXIS_NAMES.iter().enumerate() {
        ui.weak(*name);
        coordinate(ui, &mut p[a], g, a, per);
    }
}

/// A position, to four significant figures and no trailing zeros: `4.5`, `0.5`, `12.25`.
///
/// Not [`editor_core::magnitude`], which is for a value and gives every one the same width:
/// a cell centre half a millimetre in read `0.500000`, three times over on one line.
fn short(v: f64) -> String {
    if !v.is_finite() {
        return "-".to_string();
    }
    let a = v.abs();
    if a == 0.0 {
        return "0".to_string();
    }
    if !(1e-3..1e5).contains(&a) {
        return format!("{v:.3e}");
    }
    let mut decimals = 0usize;
    let mut scaled = a;
    while scaled < 1000.0 && decimals < 6 {
        scaled *= 10.0;
        decimals += 1;
    }
    let t = format!("{v:.decimals$}");
    if t.contains('.') {
        t.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        t
    }
}

/// A point, in the display unit, for a caption.
fn spelled(p: [f64; 3], unit: &str, per: f64) -> String {
    format!(
        "({}, {}, {}) {unit}",
        short(p[0] / per),
        short(p[1] / per),
        short(p[2] / per)
    )
}

/// The panel. `run` and `frame` are what is on screen, `prefer` the domain the outliner has
/// selected, and `generation` a number that changes whenever the run does.
pub fn panel(
    ui: &mut egui::Ui,
    r: &mut Results,
    run: &viewer_core::Run,
    frame: usize,
    prefer: Option<&str>,
    generation: u64,
) {
    // **The panel is as tall as it was made, whatever the tab draws.** A panel keeps the size of
    // what was drawn in it, and the readings tab's table drew less than its chart needed: the
    // panel shrank to its 120-point floor, the chart then had no room and drew nothing, and the
    // tab was a table of three rows over an empty strip.
    ui.set_min_height(ui.available_height());
    let names = post::fields(run);
    // **The outliner's selection, when it changes to a field.** Applied every frame it undid a
    // pick in the menu below the frame after it was made: with `module` selected, choosing
    // `stress` showed it once and went back. Otherwise whatever was being read, and failing
    // that the first.
    let prefer = prefer.map(str::to_string);
    if prefer != r.prefer_was {
        if let Some(p) = prefer.as_ref().filter(|p| names.contains(p)) {
            r.domain = Some(p.clone());
        }
        r.prefer_was = prefer;
    }
    if !r.domain.as_ref().is_some_and(|d| names.contains(d)) {
        r.domain = names.first().cloned();
    }
    let tab = *r.tab.get_or_insert(if names.is_empty() {
        Tab::Readings
    } else {
        Tab::Region
    });
    // Wrapping, like every row here: at a 500-point window the tabs and the domain's menu are
    // wider than the panel, and an unwrapped row widens everything under it with it.
    row(ui, |ui| {
        ui.strong("Results");
        ui.separator();
        for t in Tab::ALL {
            if ui.selectable_label(tab == t, t.label()).clicked() {
                r.tab = Some(t);
            }
        }
        if tab != Tab::Readings && !names.is_empty() {
            ui.separator();
            let mut picked = r.domain.clone().unwrap_or_default();
            egui::ComboBox::from_id_salt("results-domain")
                .selected_text(picked.as_str())
                .show_ui(ui, |ui| {
                    for n in &names {
                        ui.selectable_value(&mut picked, n.clone(), n.as_str());
                    }
                });
            r.domain = Some(picked);
        }
    });
    ui.separator();

    let tab = r.tab.unwrap_or(Tab::Readings);
    let frame = frame.min(run.frames.len().saturating_sub(1));
    let Some(now) = run.frames.get(frame) else {
        ui.weak("the run has no frames yet");
        return;
    };
    if tab == Tab::Readings {
        readings(ui, r, run, now);
        return;
    }
    let Some(name) = r.domain.clone() else {
        ui.weak(
            "nothing in this run is a field — a region, a slice, a line and a probe read one; \
             what it reports is under Readings",
        );
        return;
    };
    let g = match post::grid(now, &name) {
        Ok(g) => g,
        Err(why) => {
            ui.weak(why);
            return;
        }
    };
    r.fit_to(&g);
    // A placed domain's positions are its own, before its pose: say so, or a probe's `(50, 20,
    // 10) mm` reads as a place in the world it is not.
    if !g.place.is_here() {
        ui.weak(format!(
            "positions are in the frame of {}, before its pose moves and turns it",
            g.name
        ));
    }
    match tab {
        Tab::Region => region(ui, r, run, &g, frame, now.t, generation),
        Tab::Slice => slice(ui, r, run, &g),
        Tab::Line => line(ui, r, &g),
        Tab::Point => point(ui, r, run, &g, now.t),
        Tab::Readings => {}
    }
}

fn region(
    ui: &mut egui::Ui,
    r: &mut Results,
    run: &viewer_core::Run,
    g: &Grid,
    frame: usize,
    t: f64,
    generation: u64,
) {
    let (unit, per) = length_of(g);
    row(ui, |ui| {
        ui.checkbox(&mut r.whole, "the whole domain");
        if !r.whole {
            triple(ui, "from", &mut r.lo, g, per);
            triple(ui, "to", &mut r.hi, g, per);
            ui.label(unit);
        }
    });
    let within = r.within();
    let key = {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        generation.hash(&mut h);
        run.frames.len().hash(&mut h);
        g.name.hash(&mut h);
        within.map(|b| b.map(f64::to_bits)).hash(&mut h);
        h.finish()
    };
    if r.over_time.as_ref().map(|c| c.0) != Some(key) {
        r.over_time = Some((key, post::region_history(run, g.name, within)));
    }
    let Some((_, history)) = &r.over_time else {
        return;
    };
    let Some(stats) = history.get(frame).and_then(|h| h.1) else {
        ui.weak(if r.whole {
            "every cell of this domain is a hole in this frame"
        } else {
            "the box holds no cell of this domain"
        });
        return;
    };
    let (measure, measure_unit) = post::measure_name(stats.dims);
    // Three hundred points, or all there is: a column wider than its panel is cut off by the
    // window rather than by anything that scrolls.
    let column = ui.available_width().min(300.0);
    let fits = chart_fits(egui::vec2(
        ui.available_width() - column - ui.spacing().item_spacing.x,
        ui.available_height(),
    ));
    side_by_side(ui, |ui| {
        ui.vertical(|ui| {
            ui.set_width(column);
            if !fits {
                ui.weak(NO_ROOM);
            }
            // **In a scroll area of its own.** A short panel holds fewer rows than the table has,
            // and the table ran on under the panel's edge: at a 600-point window its last rows
            // were clipped, and the chart beside it laid its time axis under the status bar.
            egui::ScrollArea::vertical()
                .id_salt("results-region-table")
                .auto_shrink([false, true])
                .show(ui, |ui| {
                    egui::Grid::new("results-region")
                        .num_columns(2)
                        .spacing([12.0, 3.0])
                        .show(ui, |ui| {
                            let row = |ui: &mut egui::Ui, label: &str, value: String| {
                                ui.weak(label);
                                ui.monospace(value);
                                ui.end_row();
                            };
                            row(
                                ui,
                                "max",
                                format!("{} {}", editor_core::magnitude(stats.max), g.unit),
                            );
                            row(ui, "  at", spelled(stats.max_at, unit, per));
                            row(
                                ui,
                                "min",
                                format!("{} {}", editor_core::magnitude(stats.min), g.unit),
                            );
                            row(ui, "  at", spelled(stats.min_at, unit, per));
                            row(
                                ui,
                                "mean",
                                format!("{} {}", editor_core::magnitude(stats.mean), g.unit),
                            );
                            row(
                                ui,
                                "integral",
                                format!(
                                    "{} {}",
                                    editor_core::magnitude(stats.integral),
                                    post::integral_unit(g.unit, stats.dims)
                                ),
                            );
                            row(
                                ui,
                                measure,
                                format!("{} {measure_unit}", editor_core::magnitude(stats.measure)),
                            );
                            row(
                                ui,
                                "cells",
                                if stats.holes == 0 {
                                    stats.samples.to_string()
                                } else {
                                    format!(
                                        "{} and {} not in the domain",
                                        stats.samples, stats.holes
                                    )
                                },
                            );
                        });
                });
        });
        if !fits {
            return;
        }
        // The three over the run, against the cursor at the frame on screen.
        let line = |label: &str, pick: fn(&Stats) -> f64| Line {
            domain: g.name.to_string(),
            label: label.to_string(),
            points: history
                .iter()
                .map(|(t, s)| [*t, s.as_ref().map_or(f64::NAN, pick)])
                .collect(),
        };
        let lines = vec![
            line("max", |s| s.max),
            line("mean", |s| s.mean),
            line("min", |s| s.min),
        ];
        match chart::graph(Across::Time, g.unit, lines) {
            Some(c) => chart_in(ui, &c, Some(t)),
            None => {
                ui.weak("no frame has a value here to chart");
            }
        }
    });
}

fn slice(ui: &mut egui::Ui, r: &mut Results, run: &viewer_core::Run, g: &Grid) {
    let (unit, per) = length_of(g);
    let scale = run.scale_of(g.name);
    row(ui, |ui| {
        ui.label("plane across");
        let was = r.axis;
        egui::ComboBox::from_id_salt("results-axis")
            .selected_text(AXIS_NAMES[r.axis])
            .width(40.0)
            .show_ui(ui, |ui| {
                for (a, name) in AXIS_NAMES.iter().enumerate() {
                    ui.selectable_value(&mut r.axis, a, *name);
                }
            });
        // A new axis has a new range, and the old position on it means nothing.
        if r.axis != was {
            r.at = g.centre()[r.axis];
        }
        ui.label("at");
        coordinate(ui, &mut r.at, g, r.axis, per);
        ui.label(unit);
        ui.separator();
        ui.weak("in the view's colours, over the whole run");
    });
    let Some(s) = g.slice(r.axis, r.at) else {
        ui.weak("the plane is outside the domain");
        return;
    };
    // **The whole field decides the colouring, as it does in the view.** Whether a field is
    // drawn in Planck's colours depends on its hottest cell, and a plane that misses that cell
    // would otherwise be drawn on the ramp beside a view drawn in glow.
    let colouring = editor_core::Colouring::of(g.unit, g.values, scale);
    let key = {
        let mut h = std::collections::hash_map::DefaultHasher::new();
        colouring.physical().hash(&mut h);
        // And the hottest cell, which a glowing field's brightness is relative to: a plane whose
        // own values did not change is still drawn differently when the hottest cell did.
        g.values
            .iter()
            .copied()
            .filter(|v| v.is_finite())
            .fold(f64::NEG_INFINITY, f64::max)
            .to_bits()
            .hash(&mut h);
        s.counts.hash(&mut h);
        for v in &s.values {
            v.to_bits().hash(&mut h);
        }
        scale.map(|(a, b)| (a.to_bits(), b.to_bits())).hash(&mut h);
        h.finish()
    };
    if r.texture.as_ref().map(|t| t.0) != Some(key) {
        let image = egui::ColorImage::from_rgba_unmultiplied(s.counts, &s.rgba(&colouring));
        let texture = ui
            .ctx()
            .load_texture("results-slice", image, egui::TextureOptions::NEAREST);
        r.texture = Some((key, texture));
    }
    let Some((_, texture)) = &r.texture else {
        return;
    };

    let [ua, va] = s.axes;
    let caption = format!(
        "{} {} to {} across, {} {} to {} up, {unit}",
        AXIS_NAMES[ua],
        short(s.extent[0] / per),
        short(s.extent[2] / per),
        AXIS_NAMES[va],
        short(s.extent[1] / per),
        short(s.extent[3] / per),
    );
    let e = s.image_extent();
    let (mut w, mut h) = (e[2] - e[0], e[3] - e[1]);
    // A plane across a sheet is one sample tall and no metres: drawn a sample's width tall, so it
    // is a strip and not nothing.
    if h <= 0.0 {
        h = w / s.counts[0].max(1) as f64;
    }
    if w <= 0.0 {
        w = h / s.counts[1].max(1) as f64;
    }
    // **A plane along a bar meets it at one sample**, a point with no size to draw: a bar is
    // flat across, so a plane across x has nothing in it but where the bar is cut. Said, with
    // the value there, rather than drawn as nothing.
    if !(w > 0.0 && h > 0.0) {
        let v = s.values.first().copied().unwrap_or(f64::NAN);
        ui.label(format!(
            "the plane meets {} at one sample: {} {}",
            g.name,
            editor_core::magnitude(v),
            g.unit
        ));
        return;
    }
    // Room for the picture, its caption below and its colour bar to the right.
    const BAR: f32 = 110.0;
    let space = ui.available_size();
    if space.x < 8.0 + BAR + 24.0 || space.y < 4.0 + 22.0 + 24.0 {
        ui.weak("the panel is too small for the picture");
        return;
    }
    let (rect, response) = ui.allocate_exact_size(space, egui::Sense::hover());
    let room = egui::Rect::from_min_max(
        rect.min + egui::vec2(8.0, 4.0),
        rect.max - egui::vec2(8.0 + BAR, 22.0),
    );
    let px = (room.width() as f64 / w).min(room.height() as f64 / h);
    let size = egui::vec2((w * px) as f32, (h * px) as f32);
    let image = egui::Rect::from_min_size(room.min, size);
    let painter = ui.painter_at(rect);
    painter.image(
        texture.id(),
        image,
        egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
        egui::Color32::WHITE,
    );
    painter.rect_stroke(
        image,
        0.0,
        egui::Stroke::new(1.0, ui.visuals().weak_text_color()),
        egui::StrokeKind::Outside,
    );
    painter.text(
        egui::pos2(image.left(), image.bottom() + 4.0),
        egui::Align2::LEFT_TOP,
        caption,
        egui::FontId::monospace(10.0),
        ui.visuals().weak_text_color(),
    );
    // **The colour bar**, beside the picture and as tall as it: the run's range, bottom to top,
    // through the same `Colouring` the pixels went through. A picture in false colour without
    // one is a picture of more and less and not of how much.
    if let Some((lo, hi)) = scale.filter(|(lo, hi)| hi > lo) {
        // As tall as the picture, and never less than room for its two labels: a bar's slice is a
        // strip a few points high, and its ends' labels were drawn over each other.
        let tall = image.height().max(48.0).min(room.height());
        let bar = egui::Rect::from_min_size(
            egui::pos2(image.right() + 14.0, image.top()),
            egui::vec2(12.0, tall),
        );
        const STEPS: usize = 48;
        for i in 0..STEPS {
            let f = (i as f64 + 0.5) / STEPS as f64;
            let [cr, cg, cb] = colouring.srgb(lo + f * (hi - lo));
            let top = bar.bottom() - bar.height() * (i + 1) as f32 / STEPS as f32;
            let bottom = bar.bottom() - bar.height() * i as f32 / STEPS as f32;
            painter.rect_filled(
                egui::Rect::from_min_max(
                    egui::pos2(bar.left(), top),
                    egui::pos2(bar.right(), bottom),
                ),
                0.0,
                egui::Color32::from_rgb(cr, cg, cb),
            );
        }
        let ink = ui.visuals().weak_text_color();
        painter.text(
            egui::pos2(bar.right() + 6.0, bar.top()),
            egui::Align2::LEFT_TOP,
            format!("{} {}", editor_core::magnitude(hi), g.unit),
            egui::FontId::monospace(10.0),
            ink,
        );
        painter.text(
            egui::pos2(bar.right() + 6.0, bar.bottom()),
            egui::Align2::LEFT_BOTTOM,
            format!("{} {}", editor_core::magnitude(lo), g.unit),
            egui::FontId::monospace(10.0),
            ink,
        );
    }
    // Under the pointer: the value there, read off the field and not off the pixel.
    if let Some(at) = response.hover_pos().filter(|p| image.contains(*p)) {
        let u = e[0] + f64::from((at.x - image.left()) / image.width()) * (e[2] - e[0]);
        let v = e[3] - f64::from((at.y - image.top()) / image.height()) * (e[3] - e[1]);
        let u = u.clamp(s.extent[0], s.extent[2]);
        let v = v.clamp(s.extent[1], s.extent[3]);
        let p = s.point(u, v);
        let said = match g.at(p) {
            Some(value) => format!(
                "{}: {} {}",
                spelled(p, unit, per),
                editor_core::magnitude(value),
                g.unit
            ),
            None => format!("{}: not in the domain", spelled(p, unit, per)),
        };
        painter.text(
            at + egui::vec2(12.0, -12.0),
            egui::Align2::LEFT_BOTTOM,
            said,
            egui::FontId::monospace(11.0),
            ui.visuals().strong_text_color(),
        );
    }
}

fn line(ui: &mut egui::Ui, r: &mut Results, g: &Grid) {
    let (unit, per) = length_of(g);
    row(ui, |ui| {
        triple(ui, "from", &mut r.a, g, per);
        triple(ui, "to", &mut r.b, g, per);
        ui.label(unit);
        ui.separator();
        for (a, name) in AXIS_NAMES.iter().enumerate() {
            if g.extent[a + 3] > g.extent[a] && ui.small_button(format!("along {name}")).clicked() {
                let c = g.centre();
                r.a = c;
                r.a[a] = g.extent[a];
                r.b = c;
                r.b[a] = g.extent[a + 3];
            }
        }
    });
    let points = g.line(r.a, r.b, g.points_along(r.a, r.b));
    let lines = vec![Line {
        domain: g.name.to_string(),
        label: "along the line".to_string(),
        points,
    }];
    match chart::graph(Across::Distance, g.unit, lines) {
        Some(c) if chart_fits(ui.available_size()) => chart_in(ui, &c, None),
        Some(_) => {
            ui.weak(NO_ROOM);
        }
        None => {
            ui.weak("the line does not cross the domain anywhere it has a value");
        }
    }
}

fn point(ui: &mut egui::Ui, r: &mut Results, run: &viewer_core::Run, g: &Grid, t: f64) {
    let (unit, per) = length_of(g);
    row(ui, |ui| {
        triple(ui, "at", &mut r.point, g, per);
        ui.label(unit);
        ui.separator();
        match g.at(r.point) {
            Some(v) => ui.strong(format!("{} {}", editor_core::magnitude(v), g.unit)),
            None => ui.weak("no value here: outside the domain or beside a hole"),
        };
    });
    let lines = vec![Line {
        domain: g.name.to_string(),
        label: "at the point".to_string(),
        points: post::history(run, g.name, r.point),
    }];
    match chart::graph(Across::Time, g.unit, lines) {
        Some(c) if chart_fits(ui.available_size()) => chart_in(ui, &c, Some(t)),
        Some(_) => {
            ui.weak(NO_ROOM);
        }
        None => {
            ui.weak("the point has no value in any frame");
        }
    }
}

/// **A table of every reading, and the one picked against time.** The viewport's chart of a run
/// with nothing to place stacks a plot per unit, which is right for the full height of a view and
/// wrong for a strip under it: scene 17 reports a dozen readings in four units, and in two hundred
/// points the plots were drawn over each other fifty times. A table reads at any height, and a
/// chart of one reading has room to be read.
fn readings(ui: &mut egui::Ui, r: &mut Results, run: &viewer_core::Run, now: &viewer_core::Frame) {
    if now.readings.is_empty() {
        ui.weak("this run reports no readings");
        return;
    }
    if !r.reading.is_some_and(|i| i < now.readings.len()) {
        r.reading = Some(0);
    }
    let column = ui.available_width().min(380.0);
    let fits = chart_fits(egui::vec2(
        ui.available_width() - column - ui.spacing().item_spacing.x,
        ui.available_height(),
    ));
    side_by_side(ui, |ui| {
        ui.vertical(|ui| {
            ui.set_width(column);
            if !fits {
                ui.weak(NO_ROOM);
            }
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    egui::Grid::new("results-readings")
                        .num_columns(2)
                        .spacing([12.0, 2.0])
                        .striped(true)
                        .show(ui, |ui| {
                            for (i, x) in now.readings.iter().enumerate() {
                                let on = r.reading == Some(i);
                                if ui
                                    .selectable_label(on, format!("{} {}", x.domain, x.label))
                                    .clicked()
                                {
                                    r.reading = Some(i);
                                }
                                ui.monospace(format!(
                                    "{} {}",
                                    editor_core::magnitude(x.value),
                                    x.unit
                                ));
                                ui.end_row();
                            }
                        });
                });
        });
        let Some(picked) = r.reading.and_then(|i| now.readings.get(i)) else {
            return;
        };
        if !fits {
            return;
        }
        // The k-th reading of this domain and label in every frame, where it is the k-th here:
        // two of one name stay two lines, and a frame whose list is ordered differently still
        // finds the right one. A frame without it is a gap.
        let same = |x: &viewer_core::Reading| x.domain == picked.domain && x.label == picked.label;
        let i = r.reading.unwrap_or(0);
        let k = now.readings[..i].iter().filter(|x| same(x)).count();
        let points = run
            .frames
            .iter()
            .map(|f| {
                let v = f
                    .readings
                    .iter()
                    .filter(|x| same(x))
                    .nth(k)
                    .map_or(f64::NAN, |x| x.value);
                [f.t, v]
            })
            .collect();
        let line = Line {
            domain: picked.domain.clone(),
            label: picked.label.clone(),
            points,
        };
        match chart::graph(Across::Time, &picked.unit, vec![line]) {
            Some(c) => chart_in(ui, &c, Some(now.t)),
            None => {
                ui.weak("this reading has no value in any frame");
            }
        }
    });
}

/// The least room a chart is drawn in. Its value labels take 70 points on the left and its time
/// axis five labels of up to 40 points, and in less than this they are drawn over each other: at
/// a 700-point window the region's chart beside its statistics laid `0.000` to `0.020` across 64
/// points.
const CHART_MIN: egui::Vec2 = egui::vec2(320.0, 100.0);

/// Whether a chart fits in `room`.
fn chart_fits(room: egui::Vec2) -> bool {
    room.x >= CHART_MIN.x && room.y >= CHART_MIN.y
}

/// What a tab says in place of a chart it has no room for.
///
/// **Said, not left blank.** Under 921 points the inspector goes and the outliner narrows, and
/// the region's chart and the readings' had less than [`CHART_MIN`] beside their tables: each
/// was a table and an empty half, which is what a chart with nothing to draw looks like too.
const NO_ROOM: &str = "the chart needs a wider or taller panel";

/// A chart drawn into whatever room is left. The caller has asked [`chart_fits`].
fn chart_in(ui: &mut egui::Ui, chart: &chart::Chart, cursor: Option<f64>) {
    let (rect, _) = ui.allocate_exact_size(ui.available_size(), egui::Sense::hover());
    let drawing = chart.draw(f64::from(rect.width()), f64::from(rect.height()), cursor);
    paint(&ui.painter_at(rect), rect.min, &drawing);
}

/// A chart or a diagram, painted as `editor_core` drew it, from `origin`.
///
/// One painter for both places a drawing goes — the viewport, for a run with nothing to place,
/// and this panel — so they cannot come to disagree about what an anchor means.
pub fn paint(painter: &egui::Painter, origin: egui::Pos2, drawing: &chart::Drawing) {
    let at = |p: [f64; 2]| egui::pos2(origin.x + p[0] as f32, origin.y + p[1] as f32);
    let rgb = |c: [u8; 3]| egui::Color32::from_rgb(c[0], c[1], c[2]);
    for s in &drawing.strokes {
        painter.add(egui::Shape::line(
            s.points.iter().map(|p| at(*p)).collect(),
            egui::Stroke::new(s.width as f32, rgb(s.colour)),
        ));
    }
    for c in &drawing.discs {
        let outline = egui::Stroke::new(1.0, rgb(c.outline));
        match c.fill {
            Some(fill) => painter.circle(at(c.centre), c.radius as f32, rgb(fill), outline),
            None => painter.circle_stroke(at(c.centre), c.radius as f32, outline),
        };
    }
    for t in &drawing.texts {
        let align = match t.anchor {
            chart::Anchor::Left => egui::Align2::LEFT_CENTER,
            chart::Anchor::Right => egui::Align2::RIGHT_CENTER,
            chart::Anchor::Top => egui::Align2::CENTER_TOP,
            chart::Anchor::TopLeft => egui::Align2::LEFT_TOP,
            chart::Anchor::TopRight => egui::Align2::RIGHT_TOP,
            chart::Anchor::Bottom => egui::Align2::CENTER_BOTTOM,
            chart::Anchor::BottomLeft => egui::Align2::LEFT_BOTTOM,
            chart::Anchor::BottomRight => egui::Align2::RIGHT_BOTTOM,
        };
        painter.text(
            at(t.at),
            align,
            &t.text,
            egui::FontId::monospace(11.0),
            rgb(t.colour),
        );
    }
}

/// What the tab is reading, in world coordinates, for the viewport to draw over the picture:
/// polylines, and points to mark.
#[derive(Default)]
pub struct Overlay {
    /// Each a polyline.
    pub lines: Vec<Vec<[f64; 3]>>,
    /// Each a point to circle.
    pub marks: Vec<[f64; 3]>,
}

/// The overlay for what `r` is reading in `frame` of `run`: nothing until the panel has set its
/// positions for the domain, and nothing for the readings.
pub fn overlay(r: &Results, run: &viewer_core::Run, frame: usize) -> Overlay {
    let mut out = Overlay::default();
    let (Some(tab), Some(name)) = (r.tab, r.domain.as_deref()) else {
        return out;
    };
    let Some(now) = run
        .frames
        .get(frame.min(run.frames.len().saturating_sub(1)))
    else {
        return out;
    };
    let Ok(g) = post::grid(now, name) else {
        return out;
    };
    if r.set_for.as_ref().map(|s| (s.0.as_str(), s.1)) != Some((name, g.extent.map(f64::to_bits))) {
        return out;
    }
    let world = |p: [f64; 3]| g.place.apply(p);
    match tab {
        Tab::Slice => {
            if let Some(s) = g.slice(r.axis, r.at) {
                let [u0, v0, u1, v1] = s.extent;
                out.lines.push(
                    [(u0, v0), (u1, v0), (u1, v1), (u0, v1), (u0, v0)]
                        .into_iter()
                        .map(|(u, v)| world(s.point(u, v)))
                        .collect(),
                );
            }
        }
        Tab::Line => {
            out.lines.push(vec![world(r.a), world(r.b)]);
            out.marks.push(world(r.a));
            out.marks.push(world(r.b));
        }
        Tab::Point => out.marks.push(world(r.point)),
        Tab::Region => {
            // Where the largest value is, read from the statistics the panel last drew.
            if let Some(s) = r
                .over_time
                .as_ref()
                .and_then(|(_, h)| h.get(frame))
                .and_then(|h| h.1)
            {
                out.marks.push(world(s.max_at));
            }
            if let Some(b) = r.within() {
                let corners = g.place.corners_of(b);
                for (i, j) in editor_core::EDGES {
                    out.lines.push(vec![corners[i], corners[j]]);
                }
            }
        }
        Tab::Readings => {}
    }
    out
}
