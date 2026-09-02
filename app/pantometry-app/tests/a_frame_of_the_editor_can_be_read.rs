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
    let mut p = std::env::current_exe().expect("the test binary knows where it is");
    p.pop();
    if p.ends_with("deps") {
        p.pop();
    }
    let mut argv = vec!["--ui-dump".to_string()];
    argv.extend(args.iter().map(|a| (*a).to_string()));
    let out =
        std::process::Command::new(p.join(format!("pantometry{}", std::env::consts::EXE_SUFFIX)))
            .args(&argv)
            .output()
            .expect("the binary runs");
    assert!(
        out.status.success(),
        "--ui-dump {argv:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
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
    let d = dump(&[]);
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
        let d = dump(&["--width", &width.to_string()]);
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
fn nothing_is_drawn_past_the_right_edge() {
    // **A control that does not fit is not elided, it is lost.** The toolbar was a plain
    // `horizontal`, and at 700 points `fit view` was laid out at x=704..746 of a 700-point window
    // — drawn, clipped by the panel, invisible and unclickable; at 500 the last three controls
    // were not laid out at all. `elided` saw none of that, because egui does not truncate a
    // button's text: the row overflows instead. The check that stood here counted a thing that
    // cannot happen, and a 100-character label at 500 points left it reading zero.
    //
    // Measured with that label instead: the widest right edge at 900 points went 885 -> 1057 in a
    // plain row, and stays inside the window in a wrapped one. `elided` is still in the dump, as
    // a number to read rather than one to assert.
    for width in [1500, 1100, 900, 700, 620, 500, 420] {
        let d = dump(&["--width", &width.to_string()]);
        let over: Vec<&str> = d
            .lines()
            .filter(|l| l.starts_with("  "))
            .filter(|l| {
                let mut it = l.split_whitespace();
                match (it.next(), it.next(), it.next()) {
                    (Some(left), Some(_), Some(w)) => {
                        match (left.parse::<f32>(), w.parse::<f32>()) {
                            (Ok(left), Ok(w)) => left + w > width as f32,
                            _ => false,
                        }
                    }
                    _ => false,
                }
            })
            .collect();
        assert!(
            over.is_empty(),
            "drawn past the right edge at {width}: {over:?}"
        );
    }
}

#[test]
fn a_dropped_panel_says_so_rather_than_vanishing() {
    // At 500 points the inspector does not fit. The failure this guards is the silent one — a
    // panel that is simply not there, which reads as a missing feature rather than as a window
    // that is too narrow.
    let wide = dump(&["--width", "1500"]);
    let narrow = dump(&["--width", "500"]);
    assert!(wide.contains("Inspector"), "no inspector at 1500:\n{wide}");
    assert!(!narrow.contains("Inspector"), "the inspector fits at 500?");
    assert!(
        narrow.contains("the inspector hidden"),
        "the inspector left without a word:\n{narrow}"
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
    let d = dump(&[]);
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
