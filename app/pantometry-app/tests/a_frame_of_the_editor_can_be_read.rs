//! **The editor could not be looked at except by opening it.**
//!
//! `--layout-at` checks the arithmetic that decides which side panels fit; `--drawn-extent` checks
//! where the shaded pass puts a run's geometry. Neither says what is *on the screen*. So every
//! change to the interface was unverifiable, and this repository has already paid for that once:
//! three side panels took the whole of a 533-point window, the paint callback was handed a rect of
//! **zero width**, it drew nothing — correctly — and an afternoon went into the renderer before
//! anybody printed the rect.
//!
//! `--ui-dump` builds one frame from an `egui::Context` alone and reports every string it drew and
//! where. `eframe::App::update`'s `Frame` argument was already unused, so nothing here ever needed
//! eframe; the shaded viewport is a `PaintCallback`, which a headless run records and never
//! executes, so no GPU is needed either.
//!
//! # Text and rects, not pixels
//!
//! A screenshot answers "does it look right" only to somebody looking. These are the questions a
//! test can hold: **is the control inside the window at this width**, **is the viewport wider than
//! nothing**, **is the error where a reader will find it**.
//!
//! Not "was the label elided", which is what was written here first. egui does not truncate a
//! button's text — the row overflows instead — so `elided` stayed at 0 with a hundred-character
//! label in a 500-point window, and an assertion that it was 0 could not have failed. It is still
//! in the dump as a number to read. What replaced it measures each string's own right edge, and
//! that found `fit view` laid out at x=704 of a 700-point window on the first run.

/// The binary, asked for one frame.
fn dump(args: &[&str]) -> String {
    dump_listing("", args)
}

/// The binary, asked for one frame, with `recent` named as the file the start screen reads.
///
/// **Every dump names one**, including the ones that do not care: without it the start screen would
/// read whichever list belongs to whoever is running the suite, and a test whose result depends on
/// the machine it runs on is not a test. The default is a path that does not exist, which the
/// editor reads as an empty list.
fn dump_listing(recent: &str, args: &[&str]) -> String {
    let mut p = std::env::current_exe().expect("the test binary knows where it is");
    p.pop();
    if p.ends_with("deps") {
        p.pop();
    }
    let mut argv = vec!["--ui-dump".to_string()];
    argv.extend(args.iter().map(|a| (*a).to_string()));
    let none = std::env::temp_dir().join("pantometry-no-such-recent-list.json");
    let out =
        std::process::Command::new(p.join(format!("pantometry{}", std::env::consts::EXE_SUFFIX)))
            .args(&argv)
            .env(
                "PANTOMETRY_RECENT",
                if recent.is_empty() {
                    none.to_string_lossy().into_owned()
                } else {
                    recent.to_string()
                },
            )
            .output()
            .expect("the binary runs");
    assert!(
        out.status.success(),
        "--ui-dump {argv:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// A scene on disk to open, written by the binary's own `emit`.
///
/// The built-in scene, which is what the dump used to show when no file was named. Written out
/// rather than reached at a path relative to this file, because a test binary that knows where the
/// source tree is, is a test binary that keeps working after the tree has moved — this workspace
/// has already had six tests pass against a checkout that was no longer there.
fn scene() -> String {
    static PATH: std::sync::OnceLock<String> = std::sync::OnceLock::new();
    PATH.get_or_init(|| {
        let mut p = std::env::current_exe().expect("the test binary knows where it is");
        p.pop();
        if p.ends_with("deps") {
            p.pop();
        }
        let dir = std::env::temp_dir().join("pantometry-ui-dump-scene");
        std::fs::create_dir_all(&dir).expect("a temp dir");
        let out = dir.join("built-in.json");
        let done = std::process::Command::new(
            p.join(format!("pantometry{}", std::env::consts::EXE_SUFFIX)),
        )
        .args(["emit", &out.to_string_lossy()])
        .output()
        .expect("the binary runs");
        assert!(
            done.status.success() && out.is_file(),
            "emit did not write a scene: {}",
            String::from_utf8_lossy(&done.stderr)
        );
        out.to_string_lossy().into_owned()
    })
    .clone()
}

/// The viewport's rect — `left, top, width, height` — out of the header.
///
/// The rect the paint callback was handed, which is the number the original defect was: it was
/// zero points wide and the pass drew nothing, correctly, for an afternoon.
fn viewport(dump: &str) -> (f32, f32, f32, f32) {
    let line = dump
        .lines()
        .find(|l| l.starts_with("viewport="))
        .unwrap_or_else(|| panic!("no viewport in:\n{dump}"));
    let (pos, size) = line["viewport=".len()..]
        .split_once(' ')
        .expect("left,top WxH");
    let (l, t) = pos.split_once(',').expect("left,top");
    let (w, h) = size.split_once('x').expect("WxH");
    let n = |s: &str| s.parse::<f32>().expect("a number");
    (n(l), n(t), n(w), n(h))
}

/// The `n` out of a `key=n` line in the header.
fn count(dump: &str, key: &str) -> usize {
    let line = dump
        .lines()
        .find(|l| l.starts_with(&format!("{key}=")))
        .unwrap_or_else(|| panic!("no {key} in:\n{dump}"));
    line[key.len() + 1..].parse().expect("a number")
}

#[test]
fn a_frame_has_something_on_it_and_a_viewport_to_draw_in() {
    let d = dump(&[&scene()]);
    assert!(d.starts_with("size=1500x950"), "{d:.80}");
    // Thirty-six strings on the default window. A floor rather than the number, because the count
    // moves whenever a label is added and the assertion worth having is "the frame is not empty".
    assert!(count(&d, "texts") > 20, "the frame is nearly empty:\n{d}");
    // **The viewport is queued**, which is one of the two things the zero-width defect needed and
    // the weaker one: a rect of zero width is still a callback, and this assertion passed at every
    // width while the rect was empty. How wide it is, is `the_viewport_keeps_its_floor_at_every_width`.
    assert_eq!(count(&d, "callbacks"), 1, "no shaded viewport:\n{d}");
}

#[test]
fn the_viewport_keeps_its_floor_at_every_width() {
    // **The assertion the arithmetic could not make.** `--layout-at 900` said `view=290` while the
    // callback's rect was `0 x 870`: `panels_that_fit` reserved the panels' minimums — 170, 240,
    // 200 — and the panels then took their defaults — 260, 430, 320. The sum had a hole in it at
    // 1000, 950, 900 and 700 points: the same zero-width viewport that once cost an afternoon in
    // the renderer, still there, under a guard that was checking its own arithmetic rather than
    // the program. `the_viewport_always_has_room` passes against every one of those widths.
    //
    // So this reads the rect the paint callback was actually handed. A window narrower than the
    // floor plus its two gaps cannot honour the floor and is not asked to — at 260 points every
    // panel is gone and the viewport gets the 244 that exist.
    let floor = 280.0_f32;
    let chrome = 16.0_f32;
    for width in [
        1500, 1200, 1000, 950, 900, 800, 700, 600, 533, 500, 420, 360, 300, 260,
    ] {
        let d = dump(&[&scene(), "--width", &width.to_string()]);
        let (_, _, w, h) = viewport(&d);
        let want = floor.min(width as f32 - chrome);
        assert!(
            w >= want,
            "the viewport is {w} points wide at {width}, under {want}:\n{d}"
        );
        assert!(
            h > 100.0,
            "the viewport is {h} points tall at {width}:\n{d}"
        );
    }
}

#[test]
fn nothing_the_reader_should_see_is_cut_off_by_the_window() {
    // **A control that does not fit is not elided, it is lost.** The toolbar was a plain
    // `horizontal`, and at 700 points `fit view` was laid out at x=704..746 of a 700-point window
    // — drawn, clipped by the panel, invisible and unclickable. `elided` saw none of it: egui does
    // not truncate a button's text, the row overflows instead, so the assertion that stood here
    // counted a thing that cannot happen and read zero with a hundred-character label at 500.
    //
    // Nor is it "reaches past the right edge of the window", which was the next try and which the
    // path box fails honestly: a `TextEdit` lays its whole line out and clips it to its own field,
    // so an absolute path is 431 points of galley inside 220 points of box and nothing is wrong.
    // What `cut` counts is a string clipped by something that reaches the window's own edge.
    for width in [1500, 1100, 900, 700, 620, 500, 420, 300, 260] {
        let d = dump(&[&scene(), "--width", &width.to_string()]);
        assert_eq!(
            count(&d, "cut"),
            0,
            "something is cut off at {width}; the lines marked `!`:\n{d}"
        );
    }

    // **And the count is one that can move.** Disabling the condition it is computed from changed
    // nothing at any width above, because nothing is cut at any of them — which is the state to be
    // in and also the shape of an assertion that cannot fail. So a width where something *must* be
    // cut: at 40 points `watch file` is a wider string than the whole window, and wrapping a row
    // cannot help when one item is wider than the row. If this stops being 2, `cut` has stopped
    // measuring and the zeros above mean nothing.
    let tiny = dump(&[&scene(), "--width", "40"]);
    assert!(
        count(&tiny, "cut") > 0,
        "nothing is cut in a 40-point window, so `cut` is not measuring:\n{tiny}"
    );
}

#[test]
fn every_control_is_still_there_at_every_width() {
    // **The other half, and the one `cut` cannot see.** At 500 points the unwrapped toolbar did not
    // draw its last three controls at all — `watch file`, `run on change` and `fit view` were not
    // laid out, not clipped, and the status bar said only that the inspector had gone. A control
    // that is absent reads as a feature the program does not have.
    let want = [
        "open",
        "revert",
        "save",
        "run",
        "verify",
        "deep",
        "watch file",
        "run on change",
        "fit view",
    ];
    for width in [1500, 1100, 900, 700, 620, 500, 420] {
        let d = dump(&[&scene(), "--width", &width.to_string()]);
        // The text begins at column 21: a one-character marker, three five-wide numbers and two
        // spaces. Sliced rather than split, because a label has spaces in it.
        let on_screen: Vec<&str> = d.lines().filter_map(|l| l.get(21..)).collect();
        let missing: Vec<&&str> = want
            .iter()
            .filter(|w| !on_screen.iter().any(|t| t == *w))
            .collect();
        assert!(
            missing.is_empty(),
            "gone from the toolbar at {width}: {missing:?}\n{d}"
        );
    }
}

#[test]
fn a_dropped_panel_says_so_rather_than_vanishing() {
    // At 500 points the inspector does not fit. The failure this guards is the silent one — a
    // panel that is simply not there, which reads as a missing feature rather than as a window
    // that is too narrow.
    let wide = dump(&[&scene(), "--width", "1500"]);
    let narrow = dump(&[&scene(), "--width", "500"]);
    assert!(wide.contains("Inspector"), "no inspector at 1500:\n{wide}");
    assert!(!narrow.contains("Inspector"), "the inspector fits at 500?");
    assert!(
        narrow.contains("the inspector hidden"),
        "the inspector left without a word:\n{narrow}"
    );
}

#[test]
fn opening_with_no_file_offers_the_two_ways_in_and_no_editor() {
    // **What `pantometry edit` with no argument used to do was open the built-in scene**, which is
    // an answer to a question nobody asked: whoever typed it either has a scene somewhere or wants
    // to make one. Neither of those was on the screen, and the way to a file on disk was to know
    // that the toolbar's text box was a way in.
    let d = dump(&[]);
    assert!(d.contains("New project"), "no way to make one:\n{d}");
    assert!(d.contains("Open a scene"), "no way to open one:\n{d}");
    // And none of the editor, because there is nothing yet for it to act on. A window of disabled
    // controls is a worse answer to "what now" than two buttons.
    assert_eq!(count(&d, "callbacks"), 0, "a viewport with no scene:\n{d}");
    assert!(!d.contains("Outliner"), "the outliner over no scene:\n{d}");
    assert!(
        !d.contains("a small room ringing"),
        "still the built-in scene:\n{d}"
    );
}

#[test]
fn the_chooser_offers_every_kind_the_format_defines() {
    // **Derived from the table, not listed here.** A twentieth domain reaches the chooser without
    // this file learning about it, and a nineteenth that stops being offered fails rather than
    // quietly disappearing from a menu — which is how the editor's Add list is checked too.
    let d = dump(&["--new"]);
    assert_eq!(count(&d, "callbacks"), 0, "a viewport with no scene:\n{d}");
    for t in pantometry_world::templates::TEMPLATES {
        assert!(d.contains(t.kind), "{} is not offered:\n{d}", t.kind);
        assert!(
            d.contains(t.about),
            "{} is offered with nothing said about it:\n{d}",
            t.kind
        );
    }
}

#[test]
fn the_kinds_that_need_a_partner_say_so_and_only_they_do() {
    // `beam` names what it shines `onto` and `structure` names the block it `follows`. Neither is
    // wired up — three fields have to agree, not one name — so the chooser says it and the scene
    // opens on the format's own complaint. Exactly two, because a third that stopped saying it
    // would be a row that looks complete and is not.
    let d = dump(&["--new"]);
    let said = d
        .lines()
        .filter(|l| l.contains("names another domain"))
        .count();
    let want = pantometry_world::templates::TEMPLATES
        .iter()
        .filter(|t| t.needs_a_partner())
        .count();
    assert_eq!(said, want, "{said} rows say it, {want} kinds need it:\n{d}");
    assert_eq!(want, 2, "the set of kinds needing a partner has changed");
}

#[test]
fn a_set_that_cannot_share_one_duration_says_so_before_it_is_made() {
    // **The thing a chooser can say and a checker cannot.** Every combination of these kinds is a
    // *well-formed* scene, so nothing downstream refuses it: what happens is that `atoms` finishes
    // in picoseconds while a thermal `network` advances by a ten-trillionth of what it needs, and
    // the picture looks like a bug in the physics. Measured from the same numbers the scene is
    // built with, and said where the choice is being made.
    let far = dump(&["--new", "atoms,network"]);
    assert!(
        far.contains("3e14 times apart"),
        "picoseconds against half an hour passed without comment:\n{far}"
    );

    // And not said about two kinds a shipped scene already runs together, or the warning is
    // decoration rather than a measurement.
    let near = dump(&["--new", "bar,heater"]);
    assert!(
        !near.contains("times apart"),
        "two kinds that run together were warned about:\n{near}"
    );
}

#[test]
fn a_remembered_file_that_is_gone_is_shown_and_refused() {
    // A list that silently shortens itself looks like a list that forgot. The reader is then left
    // wondering whether the editor lost their scene or they imagined opening it — so a path that
    // has gone is on the screen, named, and not clickable.
    let dir = std::env::temp_dir().join("pantometry-recent-test");
    std::fs::create_dir_all(&dir).expect("a temp dir");
    let list = dir.join("recent.json");
    let here = scene();
    let gone = dir.join("deleted-since.json");
    let _ = std::fs::remove_file(&gone);
    std::fs::write(
        &list,
        serde_json::to_string(&[here.clone(), gone.to_string_lossy().into_owned()])
            .expect("two strings"),
    )
    .expect("writable");

    let d = dump_listing(&list.to_string_lossy(), &[]);
    assert!(d.contains("Recent"), "no list at all:\n{d}");
    assert!(d.contains("built-in.json"), "the scene that is there:\n{d}");
    assert!(
        d.contains("deleted-since.json — missing"),
        "the scene that is gone left without a word:\n{d}"
    );
}

#[test]
fn a_scene_that_does_not_check_out_says_why_on_the_screen() {
    // The parser knows `line:column` and the missing key by name. What this asserts is that the
    // sentence reaches the *screen* — a diagnostic that only ever appears on stderr is a
    // diagnostic an editor's user never sees.
    let dir = std::env::temp_dir().join("pantometry-ui-dump-test");
    std::fs::create_dir_all(&dir).expect("a temp dir");
    let path = dir.join("broken.json");
    std::fs::write(
        &path,
        r#"{ "title": "broken", "duration_s": 1.0, "frames": 2, "domains": [ { "kind": "block", "name": "b" } ] }"#,
    )
    .expect("writable");

    let d = dump(&[&path.to_string_lossy()]);
    assert!(
        d.contains("missing field `cells`"),
        "the parser's reason is not on screen:\n{d}"
    );
    assert!(d.contains("1:100"), "the place is not on screen:\n{d}");
}

#[test]
fn the_two_command_rows_are_both_there_and_are_not_the_same_row() {
    // **What the dump can and cannot see.** A menu's *items* live inside a closed menu and are not
    // on the screen; what is on the screen is six menu titles and the strip of buttons under them.
    // So this asserts the two rows, not their contents — an earlier version of this test was named
    // for a claim about duplication that its assertion did not make, and could not have: it counted
    // whether five strings appeared *anywhere*, which is true of five strings that appear once.
    //
    // The duplication is real and is read from the source rather than from here: `File` holds Load
    // and Save while the toolbar holds `load` and `save`, and `View` holds Fit view while the
    // toolbar holds `fit view`. Whether that is redundancy worth removing is a design question this
    // test does not answer.
    let d = dump(&[&scene()]);
    let row = |top: i32| -> Vec<String> {
        d.lines()
            .filter_map(|l| {
                let mut it = l.split_whitespace();
                let (_, y) = (it.next()?, it.next()?.parse::<i32>().ok()?);
                (y == top).then(|| l.split_whitespace().skip(3).collect::<Vec<_>>().join(" "))
            })
            .collect()
    };
    let menu = row(4);
    let bar = row(26);
    assert_eq!(
        menu,
        ["File", "Edit", "View", "Domain", "Run", "Watch"],
        "the menu row is not what it was"
    );
    assert!(
        bar.contains(&"run".to_string()) && bar.contains(&"fit view".to_string()),
        "the toolbar row is not below the menu row: {bar:?}"
    );
    assert_ne!(menu, bar, "the two rows collapsed into one");
}
