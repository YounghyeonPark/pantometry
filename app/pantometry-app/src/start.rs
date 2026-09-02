//! **The window before a scene**: make one, open one, or return to one you had open.
//!
//! The editor used to open onto the built-in scene whatever you meant, and the only way to reach a
//! file was to type its path into a text box on the toolbar — `File > Load` re-reads the path
//! already in that box, which is a revert and is now named one. So a person arriving with a scene
//! on disk had to know that the box was the way in.
//!
//! # What is on it
//!
//! Two verbs and a list. Nothing is on this screen that is not a way to start work: no banner, no
//! version, no sample gallery, no tips. The list is what this program can offer that the file
//! dialog cannot — the scenes *you* were working on — and it is short enough to read at a glance.
//!
//! # Where the list lives
//!
//! A JSON array of paths beside the platform's other per-user configuration:
//! `%APPDATA%\pantometry\recent.json` on Windows, `$XDG_CONFIG_HOME/pantometry/recent.json` or
//! `~/.config/pantometry/recent.json` elsewhere. Written with `serde_json`, which is already here,
//! rather than by turning on eframe's `persistence` feature, which would add `ron` and
//! `directories` to carry eight strings.
//!
//! `PANTOMETRY_RECENT` names the file instead when it is set. That is what the tests use, so a run
//! of the suite neither reads nor writes the list belonging to whoever is running it.

use eframe::egui;

/// One line of the list: a path, and whether it is still there.
///
/// **A remembered file that has been deleted is shown and refused, not quietly dropped.** A list
/// that silently shortens itself looks like a list that forgot, and the reader is left wondering
/// whether the editor lost their work or they imagined opening it.
pub struct Recent {
    /// The path as it was opened.
    pub path: String,
    /// Whether something is at that path now.
    pub there: bool,
}

/// What the reader asked the start screen for.
pub enum Chose {
    /// Nothing yet.
    Nothing,
    /// A new scene, which the caller decides the contents of.
    New,
    /// This file.
    Open(String),
}

/// How many paths are kept. Eight is about what fits without the list needing its own scrollbar.
const KEEP: usize = 8;

/// The file the list is kept in, or `None` if the platform will not say where configuration goes.
///
/// `PANTOMETRY_RECENT` wins when set, which is how the tests keep out of a real one.
pub fn list_file() -> Option<std::path::PathBuf> {
    if let Ok(p) = std::env::var("PANTOMETRY_RECENT") {
        return Some(std::path::PathBuf::from(p));
    }
    let dir = if cfg!(windows) {
        std::env::var_os("APPDATA").map(std::path::PathBuf::from)
    } else {
        std::env::var_os("XDG_CONFIG_HOME")
            .map(std::path::PathBuf::from)
            .or_else(|| {
                std::env::var_os("HOME").map(|h| std::path::PathBuf::from(h).join(".config"))
            })
    }?;
    Some(dir.join("pantometry").join("recent.json"))
}

/// The remembered paths, most recent first, each with whether it is still on disk.
///
/// A file that will not parse is the same as no file: the list is a convenience, and refusing to
/// start because of it would be the tail wagging the dog.
pub fn remembered() -> Vec<Recent> {
    let Some(file) = list_file() else {
        return Vec::new();
    };
    let Ok(text) = std::fs::read_to_string(&file) else {
        return Vec::new();
    };
    let paths: Vec<String> = serde_json::from_str(&text).unwrap_or_default();
    paths
        .into_iter()
        .take(KEEP)
        .map(|p| {
            let there = std::path::Path::new(&p).is_file();
            Recent { path: p, there }
        })
        .collect()
}

/// Put `path` at the front of the list, and say what went wrong if the write did not happen.
///
/// The error is returned rather than swallowed because the status bar can carry it: a list that
/// quietly stops remembering is indistinguishable from a list nobody used.
pub fn remember(path: &str) -> Result<(), String> {
    let Some(file) = list_file() else {
        return Err(String::from(
            "nowhere to keep the recent list: neither APPDATA nor XDG_CONFIG_HOME nor HOME is set",
        ));
    };
    let mut paths: Vec<String> = remembered().into_iter().map(|r| r.path).collect();
    paths.retain(|p| p != path);
    paths.insert(0, path.to_string());
    paths.truncate(KEEP);
    let text = serde_json::to_string_pretty(&paths).map_err(|e| e.to_string())?;
    if let Some(dir) = file.parent() {
        std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    }
    std::fs::write(&file, text).map_err(|e| format!("{}: {e}", file.display()))
}

/// Draw the screen and report what was clicked.
///
/// Held to a column of a fixed width in the middle of the window rather than stretched across it:
/// a button the width of a 1500-point window is a button whose label is a long way from its edge.
pub fn screen(ui: &mut egui::Ui, recent: &[Recent]) -> Chose {
    let mut chose = Chose::Nothing;
    ui.vertical_centered(|ui| {
        ui.add_space(ui.available_height() * 0.18);
        ui.allocate_ui_with_layout(
            egui::vec2(320.0, ui.available_height()),
            egui::Layout::top_down(egui::Align::Min),
            |ui| {
                ui.label(egui::RichText::new("pantometry").size(28.0).strong());
                ui.label(
                    egui::RichText::new("a scene is a physics, a geometry and a schedule").weak(),
                );
                ui.add_space(24.0);

                if ui
                    .add_sized([320.0, 32.0], egui::Button::new("New project…"))
                    .clicked()
                {
                    chose = Chose::New;
                }
                ui.add_space(6.0);
                if ui
                    .add_sized([320.0, 32.0], egui::Button::new("Open a scene…"))
                    .clicked()
                {
                    if let Some(p) = pick() {
                        chose = Chose::Open(p);
                    }
                }

                if recent.is_empty() {
                    return;
                }
                ui.add_space(24.0);
                ui.label(egui::RichText::new("Recent").weak());
                ui.add_space(4.0);
                for r in recent {
                    let name = std::path::Path::new(&r.path)
                        .file_name()
                        .map_or_else(|| r.path.clone(), |n| n.to_string_lossy().into_owned());
                    if r.there {
                        if ui
                            .add(egui::Button::new(&name).frame(false))
                            .on_hover_text(&r.path)
                            .clicked()
                        {
                            chose = Chose::Open(r.path.clone());
                        }
                    } else {
                        // Shown, greyed, and not a button. See `Recent`.
                        ui.add_enabled(
                            false,
                            egui::Button::new(format!("{name} — missing")).frame(false),
                        )
                        .on_disabled_hover_text(&r.path);
                    }
                }
            },
        );
    });
    chose
}

/// The platform's own open dialog, or `None` if it was dismissed.
///
/// Native and modal. `rfd` is two crates here — everything it needs on this platform was already
/// in the tree for the window — and on Linux it talks to the desktop portal rather than linking
/// GTK, so it adds nothing to what CI installs.
pub fn pick() -> Option<String> {
    rfd::FileDialog::new()
        .add_filter("scene", &["json"])
        .pick_file()
        .map(|p| p.to_string_lossy().into_owned())
}

/// The platform's own save dialog, or `None` if it was dismissed.
pub fn pick_save(suggest: &str) -> Option<String> {
    let path = std::path::Path::new(suggest);
    let mut d = rfd::FileDialog::new().add_filter("scene", &["json"]);
    if let Some(name) = path.file_name() {
        d = d.set_file_name(name.to_string_lossy().into_owned());
    }
    if let Some(dir) = path.parent().filter(|p| p.is_dir()) {
        d = d.set_directory(dir);
    }
    d.save_file().map(|p| p.to_string_lossy().into_owned())
}
