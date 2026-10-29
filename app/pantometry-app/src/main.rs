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
mod glyphs;
mod render;
mod start;
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
        Some("--ui-dump") => match dump_args(&rest) {
            Err(said) => {
                eprintln!("{said}");
                eprintln!("{DUMP_USAGE}");
                2
            }
            // **A refused run exits 1.** It printed its refusal as the frame and exited 0, so
            // `take.ps1`, reading the exit code, wrote it over a caption.
            Ok(asked) => match edit::ui_dump(asked) {
                Ok(frame) => {
                    print!("{frame}");
                    0
                }
                Err(refused) => {
                    eprintln!("{refused}");
                    1
                }
            },
        },

        // `verify`, `fit`, `--check`, `--emit-default`, a bare scene path, or nothing.
        _ => cli::run(&args),
    };
    std::process::exit(code);
}

const DUMP_USAGE: &str = concat!(
    "usage: pantometry --ui-dump [scene.json] [--width W] [--height H] [--new [kind,kind]]",
    " [--open AREA] [--custom] [--click X,Y]... [--drag X0,Y0,X1,Y1]... [--ran] [--iso]",
    " [--solo] [--run-after-clicks]"
);

/// **`--ui-dump`'s arguments, every one accounted for.**
///
/// Each flag used to be looked up on its own with `position`, so whatever no lookup asked for was
/// never seen. `--rann`, a scene path written after a flag, `--width 1500px` and `--click 10;20`
/// all produced a frame and exit 0: the first two the shut start screen, the third a 1500-point
/// window whatever was meant, the fourth no click. And a scene path that did not exist opened the
/// built-in scene, which is right for a window — it says so on the status line — and wrong here,
/// where `--ran` then runs the built-in scene and a test written against a mistyped path passes.
/// This is a door built for tests, and the one thing it must not do is answer a question it was
/// not asked. So every argument is consumed by exactly one rule, and anything left over is refused.
///
/// So is a value that is consumed and could not mean anything: a kind the chooser does not have,
/// `--open` or `--custom` without the chooser they act on, and a point that is not finite or not
/// in the window.
fn dump_args(rest: &[String]) -> Result<edit::Dump, String> {
    let last = pantometry_world::presets::AREA_COUNT - 1;
    let mut asked = edit::Dump {
        path: None,
        width: 1500.0,
        height: 950.0,
        choosing: None,
        click: Vec::new(),
        drag: Vec::new(),
        ran: false,
        iso: false,
        solo: false,
        open: None,
        custom: false,
        run_after: false,
    };
    // A flag that sets one value is refused the second time, rather than one of the two winning.
    let mut once: Vec<&str> = Vec::new();
    let mut i = 0;
    while i < rest.len() {
        let a = rest[i].as_str();
        let said = rest.get(i + 1).map(String::as_str);
        if matches!(a, "--width" | "--height" | "--open" | "--new") {
            if once.contains(&a) {
                return Err(format!("{a} was given twice"));
            }
            once.push(a);
        }
        let mut took = 1;
        match a {
            "--width" | "--height" => {
                let n = said
                    .and_then(|v| v.parse::<f32>().ok())
                    .filter(|n| n.is_finite() && *n > 0.0)
                    .ok_or_else(|| {
                        format!(
                            "{a} takes a positive number of points; got {}",
                            said.unwrap_or("nothing")
                        )
                    })?;
                if a == "--width" {
                    asked.width = n;
                } else {
                    asked.height = n;
                }
                took = 2;
            }
            // `--click x,y` may be given more than once, and they land in the order written.
            "--click" => {
                match said.and_then(numbers).as_deref() {
                    Some(&[x, y]) => asked.click.push((x, y)),
                    _ => {
                        return Err(format!(
                            "--click takes x,y; got {}",
                            said.unwrap_or("nothing")
                        ))
                    }
                }
                took = 2;
            }
            // `--drag x0,y0,x1,y1` may be given more than once too, after the clicks.
            "--drag" => {
                match said.and_then(numbers).as_deref() {
                    Some(&[x0, y0, x1, y1]) => asked.drag.push(((x0, y0), (x1, y1))),
                    _ => {
                        return Err(format!(
                            "--drag takes x0,y0,x1,y1; got {}",
                            said.unwrap_or("nothing")
                        ))
                    }
                }
                took = 2;
            }
            // **A bad index is refused, not ignored.** `--open 99`, `--open banana`, `--open -1`
            // and a bare `--open` all used to produce a dump byte-identical to the shut screen and
            // exit 0 — and the shut screen satisfies every geometric assertion the tests make
            // about an open one, so a mistyped index bought a green.
            "--open" => {
                match said
                    .and_then(|n| n.parse::<usize>().ok())
                    .filter(|n| *n <= last)
                {
                    Some(n) => asked.open = Some(n),
                    None => {
                        return Err(format!(
                            "--open takes an area from 0 to {last}; got {}",
                            said.unwrap_or("nothing")
                        ))
                    }
                }
                took = 2;
            }
            // `--new` opens on the chooser, and takes the kinds to tick as one comma-separated
            // word so a test can ask for a pair whose timescales are far apart.
            "--new" => match said.filter(|w| !w.starts_with('-')) {
                Some(w) => {
                    // **A kind is one the chooser has.** It ticks a box by asking whether the set
                    // contains that box's kind, so a name that is no kind ticks nothing and says
                    // nothing: `--new welll,orbit` was the screen with `orbit` alone.
                    let kinds = pantometry_world::templates::TEMPLATES.map(|t| t.kind);
                    for k in w.split(',') {
                        if !kinds.contains(&k) {
                            return Err(format!(
                                "--new takes kinds from {}; got {k:?}",
                                kinds.join(",")
                            ));
                        }
                    }
                    asked.choosing = Some(w.split(',').map(str::to_string).collect());
                    took = 2;
                }
                None => asked.choosing = Some(Vec::new()),
            },
            // The kinds screen with *nothing* ticked had no door: `--new` alone is the preset
            // screen and `--new bar` is the kinds screen with `bar` ticked, so the empty chooser —
            // which is what `Custom…` actually opens on — could not be looked at.
            "--custom" => asked.custom = true,
            "--ran" => asked.ran = true,
            "--iso" => {
                asked.ran = true;
                asked.iso = true;
            }
            "--solo" => asked.solo = true,
            // `--run-after-clicks` runs whatever the clicks opened. `--ran` runs first, which
            // cannot reach a preset: a preset exists only once its tile has been clicked.
            "--run-after-clicks" => asked.run_after = true,
            _ if i == 0 && !a.starts_with('-') => asked.path = Some(a.to_string()),
            _ if a.starts_with('-') => return Err(format!("--ui-dump does not take {a}")),
            _ => {
                return Err(format!(
                    "--ui-dump takes a scene only first; got {a} after a flag"
                ))
            }
        }
        i += took;
    }
    // **Two flags that only mean something on the chooser.** The dump reads them only inside
    // `--new`, so without it `--open 3` and `--custom` were accepted and drew the editor.
    if asked.choosing.is_none() {
        if asked.open.is_some() {
            return Err(String::from("--open takes effect only with --new"));
        }
        if asked.custom {
            return Err(String::from("--custom takes effect only with --new"));
        }
    }
    // **A point outside the window lands on nothing**, which is what a mistyped click used to do.
    let points = asked
        .click
        .iter()
        .chain(asked.drag.iter().flat_map(|(a, b)| [a, b]));
    for &(x, y) in points {
        if !(0.0..=asked.width).contains(&x) || !(0.0..=asked.height).contains(&y) {
            return Err(format!(
                "{x},{y} is outside the {}x{} window",
                asked.width, asked.height
            ));
        }
    }
    // `read_to_string`, because that is the read the editor opens a scene by: a file that exists
    // and is not UTF-8 — what PowerShell 5.1's `>` writes — failed it and opened the built-in
    // scene, and `fs::read` would have let that through.
    if let Some(p) = &asked.path {
        if let Err(e) = std::fs::read_to_string(p) {
            return Err(format!(
                "{p}: {e}; a dump does not stand the built-in scene in for it"
            ));
        }
    }
    Ok(asked)
}

/// Comma-separated finite numbers, or `None`. `f32` parses `NaN` and `inf`, and a click at either
/// lands nowhere.
fn numbers(v: &str) -> Option<Vec<f32>> {
    v.split(',')
        .map(|n| n.parse::<f32>().ok().filter(|n| n.is_finite()))
        .collect()
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
