//! **One binary above the library.** Run a scene, check it, verify it, watch it, edit it.
//!
//! ```text
//! pantometry                              the built-in scene, checked, nothing written
//! pantometry check   scene.json           does it parse and build, without running it
//! pantometry run     scene.json out.html  run it and write an asset, chosen by the extension
//! pantometry verify  scene.json [--deep]  the battery: margins, determinism, sweeps, geometry
//! pantometry emit    scene.json           write the built-in scene out to start from
//! pantometry fit     part.stl ...         choose a cell size by measuring the geometry
//! pantometry view    run.json             a window: rotate, zoom, scrub
//! pantometry edit    scene.json [--run]   the editor: check as you type, run, verify
//! ```
//!
//! # Why one binary and not four
//!
//! There were four, in four workspaces: the scene format's CLI, a wgpu viewer, an egui editor and a
//! compute accelerator. They were split for a reason that is still true — **a GUI stack, a GPU stack
//! and a libpython link are dependency trees the published crates must not carry** — but that is an
//! argument about the boundary between the library and everything above it, and it was never an
//! argument for three boundaries *above* it. All four are things you do to a run; they share the
//! camera, the colour scale and the scene format; and they were held apart by nothing but the order
//! they were written in.
//!
//! What it bought, concretely: `"device": "gpu"` in a scene works from the command line now. The
//! scene format refuses that key on its own — it has no device and cannot acquire one — so honouring
//! it needed an application linking both the format and an accelerator, and there was no such
//! application. See `cli::build`.
//!
//! # The old spellings still work
//!
//! `pantometry --check s.json`, `verify s.json` and `--emit-default s.json` are what
//! `pantometry-world` took, and a scene file is not the only thing with users. Anything this
//! dispatcher does not recognise goes to the CLI unchanged.

mod cli;
mod edit;
mod render;
mod view;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let rest: Vec<String> = args.iter().skip(1).cloned().collect();

    let code = match args.first().map(String::as_str) {
        Some("view") => view::run(&rest),
        Some("edit") => edit::run(&rest),

        // The subcommands the CLI already understood, under the spellings a person would guess.
        // Translated rather than reimplemented: `cli` is the same code that shipped, and a second
        // parser for the same arguments is a second thing to get out of step.
        Some("check") => cli::run(&with("--check", &rest)),
        Some("emit") => cli::run(&with("--emit-default", &rest)),
        Some("run") => cli::run(&rest),

        Some("--help" | "-h" | "help") => {
            print!("{}", usage());
            0
        }

        // **What the editor's panels do at a given window width, without a window.** The viewport
        // once got a rect of zero width because three side panels are laid out before the central
        // one, and that decision is a pure function of one number — so it is checked as one, by
        // `the_viewport_always_has_room`, rather than by resizing a window and looking.
        Some("--layout-at") => match rest.first().and_then(|w| w.parse::<f32>().ok()) {
            Some(width) => {
                println!("{}", edit::layout_at(width));
                0
            }
            None => {
                eprintln!("usage: pantometry --layout-at <width in points>");
                2
            }
        },

        // **Where the editor's shaded pass puts a run's geometry, without a window.** Six places
        // turn a panel into vertices and each has to apply the panel's placement; a missed one
        // draws that panel at the origin and looks right on every scene that states no pose. See
        // `a_placed_run_is_drawn_where_it_says`.
        Some("--drawn-extent") => {
            let mut it = rest.iter();
            match (it.next(), it.next()) {
                (Some(run), scene) => match edit::drawn_extent(scene.cloned(), run) {
                    Ok(line) => {
                        println!("{line}");
                        0
                    }
                    Err(e) => {
                        eprintln!("{e}");
                        1
                    }
                },
                _ => {
                    eprintln!("usage: pantometry --drawn-extent <run.json> [scene.json]");
                    2
                }
            }
        }

        // **One frame of the editor's interface, as text, without a window.** The editor could
        // not be looked at except by opening it, which is where the HTML report was before
        // `tools/report-check` — and the viewport-with-no-width bug is what that costs. See
        // `a_frame_of_the_editor_can_be_read`.
        Some("--ui-dump") => {
            let mut it = rest.iter();
            let path = it.next().filter(|p| !p.starts_with('-')).cloned();
            let width = rest
                .iter()
                .position(|a| a == "--width")
                .and_then(|i| rest.get(i + 1))
                .and_then(|w| w.parse().ok())
                .unwrap_or(1500.0);
            let height = rest
                .iter()
                .position(|a| a == "--height")
                .and_then(|i| rest.get(i + 1))
                .and_then(|h| h.parse().ok())
                .unwrap_or(950.0);
            print!("{}", edit::ui_dump(path, width, height));
            0
        }

        // `verify`, `fit`, `--check`, `--emit-default`, a bare scene path, or nothing.
        _ => cli::run(&args),
    };
    std::process::exit(code);
}

/// `flag` in front of `rest`, which is how the CLI already spells these.
fn with(flag: &str, rest: &[String]) -> Vec<String> {
    let mut out = vec![flag.to_string()];
    out.extend_from_slice(rest);
    out
}

fn usage() -> String {
    String::from(
        "pantometry — run, check, verify, view and edit scenes\n\
         \n\
         \x20 pantometry                              the built-in scene, checked, nothing written\n\
         \x20 pantometry check   scene.json           does it parse and build, without running it\n\
         \x20 pantometry run     scene.json out.html  run it and write an asset, by extension:\n\
         \x20                                         .html .svg .csv .json .gltf .usda\n\
         \x20 pantometry verify  scene.json [--deep]  the battery: margins, determinism, sweeps, geometry\n\
         \x20 pantometry emit    scene.json           write the built-in scene out to start from\n\
         \x20 pantometry fit     part.stl ...         choose a cell size by measuring the geometry\n\
         \x20 pantometry view    run.json             a window: rotate, zoom, scrub\n\
         \x20 pantometry edit    scene.json [--run]   check as you type, run, verify\n\
         \n\
         A scene may say \"device\": \"gpu\" on a block; this binary honours it. The CPU is the\n\
         reference either way and the difference is measured, not asserted — see pantometry-gpu.\n",
    )
}
