//! **A parametric sweep, in a window.**
//!
//! Run → Sweep a value… opens it. One number of the scene, chosen from every number the scene states
//! — the settings form's own listing, so its names are the form's — a range of values for it, and
//! a run per value in the background, with the table and the chart filling in as each one arrives.
//! [`editor_core::sweep`] does every part of it that is not a widget, which is the part the CLI's
//! `pantometry sweep` calls too: one sweep, two doors.
//!
//! **The table is of the sweep that ran, not the one being set up.** Changing the value or the range
//! after a sweep leaves its table and chart as they were, labelled with what they swept, until the
//! next one runs. A table relabelled by an edit it was not run with would be a table of numbers for
//! a question nobody asked.

use crate::results;
use editor_core::sweep::{self, Parameter, Point};
use eframe::egui;

/// What the window has been asked to do, for the editor to do it: it owns the job and the stop flag.
pub enum Asked {
    /// Run this value at these values.
    Run(Parameter, Vec<f64>),
    /// Stop between runs.
    Stop,
}

/// The window's state.
#[derive(Default)]
pub struct Sweeping {
    /// Whether it is on screen.
    pub open: bool,
    /// The value to sweep, by pointer, or `None` until one is chosen — then the first a domain
    /// states.
    pub pointer: Option<String>,
    /// The first value.
    pub from: f64,
    /// The last.
    pub to: f64,
    /// How many, both ends included.
    pub count: usize,
    /// Spaced evenly in the logarithm rather than the value.
    pub log: bool,
    /// The values written out, `5,10,20`, which take the place of the range when there are any:
    /// a sweep is as often a few values somebody has in mind as an even spread.
    pub listed: String,
    /// The pointer the range was set for, so choosing another value starts from its own.
    set_for: Option<String>,
    /// The value the table is of, as it was when the sweep started.
    pub ran: Option<Parameter>,
    /// How many values the sweep that ran asked for.
    pub asked: usize,
    /// What came back, in order.
    pub points: Vec<Point>,
    /// The reading the table and the chart show, when one was picked. `None` is the first that
    /// moves, worked out afresh as points arrive: worked out once, on the first point, nothing could
    /// have moved yet and every sweep opened on its first reading.
    pub column: Option<usize>,
    /// Stop was pressed and the value already running has not finished.
    pub stopping: bool,
    /// How the sweep that ran ended, said above its table: finished, stopped, refused.
    pub ended: Option<String>,
}

impl Sweeping {
    /// Start a range for `p`: half to twice what it is, five values, or the next four whole
    /// numbers above a count. A value of zero has no half, so it runs from zero to one.
    fn range_for(&mut self, p: &Parameter) {
        self.set_for = Some(p.pointer.clone());
        self.log = false;
        self.listed.clear();
        if p.integral {
            self.from = p.now;
            self.to = p.now + 4.0;
            self.count = 5;
        } else if p.now == 0.0 {
            (self.from, self.to, self.count) = (0.0, 1.0, 5);
        } else {
            (self.from, self.to, self.count) = (p.now * 0.5, p.now * 2.0, 5);
        }
    }

    /// The sweep a `--ui-dump --sweep` asked for, as though it had been set up in the window.
    pub fn set(&mut self, p: &Parameter, values: &[f64]) {
        self.open = true;
        self.pointer = Some(p.pointer.clone());
        self.set_for = Some(p.pointer.clone());
        self.from = values.first().copied().unwrap_or(p.now);
        self.to = values.last().copied().unwrap_or(p.now);
        self.count = values.len();
        self.log = false;
        self.listed = values
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");
    }

    /// A sweep starting: the table empties and is about this value now.
    pub fn started(&mut self, p: Parameter, asked: usize) {
        self.ran = Some(p);
        self.asked = asked;
        self.points.clear();
        self.column = None;
        self.stopping = false;
        self.ended = None;
    }

    /// A sweep over, and the sentence that says how: every value ran, some refused, or a stop cut
    /// it short. Kept above the table, because the status bar's copy of it is overwritten by the
    /// next thing the editor does — a rerun a watched file owes, for one.
    pub fn finished(&mut self) -> String {
        let (ran, refused) = (
            self.points.len(),
            self.points.iter().filter(|p| p.error.is_some()).count(),
        );
        let said = match (ran == self.asked, refused) {
            (true, 0) => format!("swept {ran} value(s)"),
            (true, n) => format!("swept {ran} value(s), {n} refused — the table says why"),
            (false, _) => format!(
                "stopped after {ran} of {} — the table is what ran",
                self.asked
            ),
        };
        self.stopping = false;
        self.ended = Some(said.clone());
        said
    }

    /// The window, when it is open. `running` is whether a sweep is in flight, and `blocked` why
    /// one cannot start — a run or a verify the editor is already doing.
    pub fn show(
        &mut self,
        ctx: &egui::Context,
        text: &str,
        running: bool,
        blocked: Option<&str>,
    ) -> Option<Asked> {
        if !self.open {
            return None;
        }
        let mut open = true;
        let mut asked = None;
        // Over the middle of the window, where the view is, rather than over the menus and the
        // tree at the top left, which is where egui puts a window it is told nothing about.
        egui::Window::new("Sweep")
            .open(&mut open)
            .pivot(egui::Align2::CENTER_CENTER)
            .default_pos(ctx.content_rect().center())
            .default_size([560.0, 460.0])
            .show(ctx, |ui| asked = self.body(ui, text, running, blocked));
        self.open = open;
        asked
    }

    fn body(
        &mut self,
        ui: &mut egui::Ui,
        text: &str,
        running: bool,
        blocked: Option<&str>,
    ) -> Option<Asked> {
        let mut asked = None;
        let params = sweep::parameters(text);
        // **Setting up needs a scene that parses; a sweep in flight does not.** Returning here hid
        // the Stop button, the count and the table behind a stray keystroke mid-sweep.
        let values = if params.is_empty() {
            ui.weak("this scene states no number to sweep — it has to parse first");
            None
        } else {
            Some(self.setup(ui, text, &params))
        };
        row(ui, |ui| {
            if running {
                if ui.button("Stop").clicked() {
                    self.stopping = true;
                    asked = Some(Asked::Stop);
                }
                ui.label(format!("{} of {}", self.points.len(), self.asked));
                if self.stopping {
                    ui.weak("stopping after this value…");
                }
            } else if let Some((chosen, values)) = values {
                let ready = values.is_ok() && blocked.is_none();
                let button = ui.add_enabled(ready, egui::Button::new("Run sweep"));
                if let Some(why) = blocked {
                    ui.weak(format!("not while {why}"));
                }
                if button.clicked() {
                    if let Ok(v) = values {
                        asked = Some(Asked::Run(chosen, v));
                    }
                }
            }
        });
        ui.separator();
        self.results(ui);
        asked
    }

    /// The value, the range and the preview of what a run would do: the chosen value and the
    /// values, or why there are none.
    fn setup(
        &mut self,
        ui: &mut egui::Ui,
        text: &str,
        params: &[Parameter],
    ) -> (Parameter, Result<Vec<f64>, String>) {
        let chosen = self
            .pointer
            .as_ref()
            .and_then(|ptr| params.iter().find(|p| &p.pointer == ptr))
            .or_else(|| params.iter().find(|p| p.pointer.starts_with("/domains/")))
            .unwrap_or(&params[0])
            .clone();
        if self.set_for.as_ref() != Some(&chosen.pointer) {
            self.range_for(&chosen);
        }
        let mut picked = chosen.pointer.clone();
        row(ui, |ui| {
            ui.label("sweep");
            egui::ComboBox::from_id_salt("sweep-value")
                .selected_text(&chosen.label)
                .width(240.0)
                .show_ui(ui, |ui| {
                    for p in params {
                        ui.selectable_value(&mut picked, p.pointer.clone(), &p.label);
                    }
                });
            ui.weak(format!(
                "now {} {}",
                results::short(chosen.now),
                chosen.unit
            ));
        });
        self.pointer = Some(picked);

        let step = if chosen.integral {
            1.0
        } else {
            (self.to - self.from).abs() / 100.0
        };
        row(ui, |ui| {
            ui.label("from");
            ui.add(egui::DragValue::new(&mut self.from).speed(step.max(1e-12)));
            ui.label("to");
            ui.add(egui::DragValue::new(&mut self.to).speed(step.max(1e-12)));
            ui.label(chosen.unit);
            ui.label("values");
            ui.add(egui::DragValue::new(&mut self.count).range(1..=50));
            ui.checkbox(&mut self.log, "log");
        });
        row(ui, |ui| {
            ui.label("or the values");
            ui.add(
                egui::TextEdit::singleline(&mut self.listed)
                    .hint_text("5, 10, 20")
                    .desired_width(160.0),
            );
        });
        // The same reading `pantometry sweep` gives `NAME=…`, so a whole number refuses a fraction
        // here as it does there.
        let values = if self.listed.trim().is_empty() {
            sweep::spaced(self.from, self.to, self.count, self.log, chosen.integral)
        } else {
            sweep::asked(text, &format!("{}={}", chosen.pointer, self.listed.trim()))
                .map(|(_, v)| v)
        };
        match &values {
            Ok(v) => {
                let shown: Vec<String> = v.iter().take(6).map(|x| results::short(*x)).collect();
                let more = if v.len() > 6 { ", …" } else { "" };
                ui.weak(format!(
                    "{} run(s), at {}{more} {}",
                    v.len(),
                    shown.join(", "),
                    chosen.unit
                ));
            }
            Err(why) => {
                ui.colored_label(egui::Color32::from_rgb(230, 180, 60), why);
            }
        }
        (chosen, values)
    }

    /// The table and the chart of the sweep that ran.
    fn results(&mut self, ui: &mut egui::Ui) {
        let Some(ran) = self.ran.clone() else {
            ui.weak("no sweep yet: choose a value and a range, and run it");
            return;
        };
        if let Some(said) = &self.ended {
            ui.label(said);
        }
        if self.points.is_empty() {
            ui.weak(format!("sweeping {}…", ran.label));
            return;
        }
        let cols = sweep::columns(&self.points);
        let ran_ok = sweep::ran(&self.points);
        if cols.is_empty() {
            // Two different emptinesses, said apart: nothing ran, or it all ran and the scene
            // reports nothing a table can hold.
            ui.weak(if ran_ok == 0 {
                format!("no value of {} ran — each says why below", ran.label)
            } else {
                format!(
                    "{ran_ok} value(s) ran, and this scene reports no readings to tabulate — \
                     its answer is a picture: run it to see"
                )
            });
        }
        let pick = self
            .column
            .filter(|k| *k < cols.len())
            .or_else(|| sweep::moving(&self.points, &cols))
            .unwrap_or(0);
        if !cols.is_empty() {
            let mut chosen = pick;
            row(ui, |ui| {
                ui.label("reading");
                egui::ComboBox::from_id_salt("sweep-reading")
                    .selected_text(cols[pick].name())
                    .width(240.0)
                    .show_ui(ui, |ui| {
                        for (k, c) in cols.iter().enumerate() {
                            ui.selectable_value(&mut chosen, k, c.name());
                        }
                    });
                if cols[pick].reported < ran_ok {
                    ui.weak(format!("reported by {} of {ran_ok}", cols[pick].reported));
                }
            });
            // Kept only when somebody chose it; otherwise worked out again as points arrive.
            if chosen != pick {
                self.column = Some(chosen);
            }
        }
        egui::ScrollArea::vertical()
            .id_salt("sweep-table")
            .max_height(120.0)
            .auto_shrink([false, true])
            .show(ui, |ui| {
                egui::Grid::new("sweep-rows")
                    .num_columns(2)
                    .striped(true)
                    .show(ui, |ui| {
                        ui.strong(if ran.unit.is_empty() {
                            ran.label.clone()
                        } else {
                            format!("{} / {}", ran.label, ran.unit)
                        });
                        ui.strong(cols.get(pick).map(|c| c.name()).unwrap_or_default());
                        ui.end_row();
                        for p in &self.points {
                            ui.monospace(results::short(p.value));
                            match (&p.error, sweep::at(p, &cols, pick)) {
                                (Some(e), _) => {
                                    ui.add(
                                        egui::Label::new(
                                            egui::RichText::new(format!("refused: {e}"))
                                                .color(egui::Color32::from_rgb(230, 180, 60)),
                                        )
                                        .truncate(),
                                    )
                                    .on_hover_text(e);
                                }
                                (None, Some(v)) => {
                                    ui.monospace(editor_core::magnitude(v));
                                }
                                (None, None) => {
                                    ui.weak("not reported at this value");
                                }
                            }
                            ui.end_row();
                        }
                    });
            });
        match sweep::chart(&ran, &self.points, pick) {
            Some(c) if results::chart_fits(ui.available_size()) => results::chart_in(ui, &c, None),
            Some(_) => {
                ui.weak(results::NO_ROOM);
            }
            None => {
                ui.weak("nothing to chart: no value reported this reading");
            }
        }
    }
}

/// A wrapping row, its labels laid out whole — the results panel's [`results::row`].
fn row(ui: &mut egui::Ui, add: impl FnOnce(&mut egui::Ui)) {
    results::row(ui, add);
}
