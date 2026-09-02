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

/// Where a dump line's text begins: two marker columns, three five-wide numbers, two spaces.
const TEXT_AT: usize = 22;

/// What a dump line says, or `None` if it is not a text line.
///
/// **Written once.** This column was written down in four places and then the dump grew a second
/// marker — `!` for a string the window cut off, `=` for a string drawn on another — and four
/// tests began failing with *"no menu titled File"*. A constant nobody can get wrong in one place
/// is worth more than four correct copies.
fn said(line: &str) -> Option<&str> {
    line.get(TEXT_AT..).filter(|t| !t.is_empty())
}

/// A dump line's `left`, `top` and `width`.
///
/// Parsed out of the columns before the text, so a marker in either of them is skipped rather
/// than counted as a field.
fn box_of(line: &str) -> Option<(f32, f32, f32)> {
    let mut n = line
        .get(..TEXT_AT)?
        .split_whitespace()
        .filter_map(|w| w.parse::<f32>().ok());
    Some((n.next()?, n.next()?, n.next()?))
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
fn nothing_is_drawn_on_top_of_anything_else() {
    // **The viewport's bottom was claimed four times.** The scale bar, the colour bar, the count
    // of what the shaded pass drew and the frame transport each anchored to `rect.bottom()` on
    // their own, and with a run open the slider was laid out through the middle of the count.
    // Legible in a screenshot of the window and invisible to every other column of the dump,
    // which is why `overlap` exists: `13`, `frame` and `0 triangles, 0 lines, 0 paints` were
    // marked `=` on three lines four points apart.
    //
    // A band each now — 26 points, which is what `colour_bar` and `scale_bar` occupy — and the
    // transport is two rows rather than one so the reservation above it can be counted rather
    // than measured. That count is the fragile part, and this is what holds it.
    //
    // **Two of the five placements are slack, measured.** Sabotage moved the draw count back to
    // the viewport's bottom edge and halved the transport's reserved rows, and neither made this
    // fail: without the ladder those two only *abut* their neighbours, by a point or two, and
    // `overlap` shrinks each rect by 2 before intersecting so that two adjacent rows of text are
    // not a collision. So this holds the reservation and the colour bar's band, and the other two
    // lines buy room to read rather than correctness. Saying which is which is cheaper than a
    // reader assuming all four are load-bearing.
    for width in [1500, 1100, 900, 700, 620, 500, 420, 300, 260] {
        for flag in ["--ran", "--iso"] {
            let d = dump(&[&scene(), flag, "--width", &width.to_string()]);
            assert_eq!(
                count(&d, "overlap"),
                0,
                "strings are drawn on each other at {width} with {flag}; \
                 the lines marked `=`:\n{d}"
            );
            assert_eq!(
                count(&d, "cut"),
                0,
                "something is cut off at {width} with {flag}; the lines marked `!`:\n{d}"
            );
        }
    }

    // **And the count is one that can move.** Every width above reads zero, which is the state to
    // be in and also the shape of an assertion that cannot fail. A menu drawn over the panel
    // behind it is an overlap by construction — that is what a menu is — so opening one is the
    // case where the number must not be zero.
    let plain = dump(&[&scene()]);
    let (x, y) = menu_at(&plain, "View");
    let opened = dump(&[&scene(), "--click", &format!("{x},{y}")]);
    assert!(
        count(&opened, "overlap") > 0,
        "an open menu overlaps nothing, so `overlap` is not measuring:\n{opened}"
    );
}

#[test]
fn every_control_is_still_there_at_every_width() {
    // **The other half, and the one `cut` cannot see.** At 500 points the unwrapped toolbar did not
    // draw its last three controls at all — `watch file`, `run on change` and `fit view` were not
    // laid out, not clipped, and the status bar said only that the inspector had gone. A control
    // that is absent reads as a feature the program does not have.
    // What is left after seven of the nine went back to the menus that already held them. `deep`
    // stays because it changes what `verify` does when you press it; `watch file` and `run on
    // change` went because they change what happens when nobody presses anything, and that is a
    // thing to read in the status bar rather than a control to keep room for.
    let want = ["run", "verify", "deep", "fit view"];
    for width in [1500, 1100, 900, 700, 620, 500, 420] {
        let d = dump(&[&scene(), "--width", &width.to_string()]);
        // The text begins at column 21: a one-character marker, three five-wide numbers and two
        // spaces. Sliced rather than split, because a label has spaces in it.
        let on_screen: Vec<&str> = d.lines().filter_map(said).collect();
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

/// Where a menu title sits, so a click can be aimed at it.
fn menu_at(dump: &str, title: &str) -> (f32, f32) {
    let line = dump
        .lines()
        .find(|l| said(l) == Some(title) && box_of(l).is_some_and(|(_, top, _)| top == 4.0))
        .unwrap_or_else(|| panic!("no menu titled {title}:\n{dump}"));
    let (left, _, w) = box_of(line).expect("a menu title's box");
    (left + w / 2.0, 10.0)
}

#[test]
fn what_left_the_toolbar_is_in_the_menu_it_was_taken_from() {
    // **The claim that justified the removal, checked rather than asserted.** Seven controls came
    // off the toolbar because a menu already held them, and until `--click` there was no way to
    // see inside a menu: a frame built with no input shows six titles and nothing under them, so
    // "it still has a home" was exactly the kind of sentence this repository has learned to
    // distrust. A press and a release at the title opens one.
    let plain = dump(&[&scene()]);
    for (title, items) in [
        ("File", &["Open…", "Revert", "Save", "Save as…"][..]),
        ("Watch", &["Watch file", "Run on change"][..]),
        ("View", &["Fit view"][..]),
    ] {
        let (x, y) = menu_at(&plain, title);
        let opened = dump(&[&scene(), "--click", &format!("{x},{y}")]);
        for item in items {
            assert!(
                opened.contains(item),
                "{title} does not hold {item} — it left the toolbar for nowhere:\n{opened}"
            );
        }
        // And the menu really opened, rather than the item having been on screen all along.
        assert!(
            items.iter().any(|i| !plain.contains(*i)),
            "{title}'s items were already on screen without opening it, so this proves nothing"
        );
    }
}

#[test]
fn a_run_can_be_looked_at_without_an_event_loop() {
    // **Half the editor does not exist until a run does**: the colour bar, the readings, the
    // transport, the probe and the isosurface strip. The window reaches them by streaming on a
    // thread and reporting through a channel, which a frame built with no event loop has nowhere
    // to receive — so until `--ran` the dump could see the empty half and nothing else, which is
    // the position the whole editor was in before this hook existed.
    let before = dump(&[&scene()]);
    let after = dump(&[&scene(), "--ran"]);
    assert!(
        !before.contains("frame") || !before.contains("ran: "),
        "the editor had already run without being asked:\n{before}"
    );
    assert!(
        after.contains("ran: 12 frames"),
        "no run happened:\n{after}"
    );
    assert!(
        after.contains("t = 0.020000 s"),
        "no time on the transport:\n{after}"
    );
}

#[test]
fn the_isosurface_level_is_on_the_picture_it_changes() {
    // **A menu is where a thing is turned on; it is not where a continuous value is set.** The
    // level's whole point is watching the surface move as you drag it, and it lived inside the
    // View menu — which is drawn over the viewport it moves in. The toggle stayed; the slider and
    // the two notes that go with it are on the strip at the bottom of the viewport, under the
    // colour bar whose scale they share.
    let d = dump(&[&scene(), "--iso"]);
    assert_eq!(
        d.matches("isosurface level").count(),
        1,
        "the level is in two places, or none:\n{d}"
    );

    // On the picture: inside the viewport's own rect, not in a panel beside it.
    let (vx, vy, vw, vh) = viewport(&d);
    let line = d
        .lines()
        .find(|l| said(l) == Some("isosurface level"))
        .expect("the slider's line");
    let (x, y, _) = box_of(line).expect("the slider's box");
    assert!(
        x >= vx && x <= vx + vw && y >= vy && y <= vy + vh,
        "the level is at {x},{y}, outside the viewport at {vx},{vy} {vw}x{vh}:\n{d}"
    );

    // And the toggle is still in the menu, where turning a thing on belongs.
    let (mx, my) = menu_at(&d, "View");
    let opened = dump(&[&scene(), "--iso", "--click", &format!("{mx},{my}")]);
    assert!(
        opened.contains("Isosurface at a value"),
        "the toggle left with the slider:\n{opened}"
    );
    assert_eq!(
        opened.matches("isosurface level").count(),
        1,
        "opening the menu produced a second level control:\n{opened}"
    );
}

#[test]
fn the_panels_are_in_a_menu_about_the_window_and_not_about_the_view() {
    // **`View` was doing two jobs.** Eleven items across three subjects — how the viewport draws,
    // which panels are on screen, and where the camera goes — while every other menu did one.
    // The three that are not about the picture at all are the window's furniture, and they have a
    // menu named for that now: a seventh name to read against three fewer items in the largest
    // menu, which is the trade this makes on purpose.
    //
    // Counted rather than searched, because `Outliner` and `Inspector` are also panel headings on
    // the screen behind the menu: opening the right menu makes each appear a second time.
    let plain = dump(&[&scene()]);
    let count = |d: &str, needle: &str| d.matches(needle).count();

    let (wx, wy) = menu_at(&plain, "Window");
    let window = dump(&[&scene(), "--click", &format!("{wx},{wy}")]);
    for item in ["Outliner", "Inspector", "Scene text"] {
        assert!(
            count(&window, item) > count(&plain, item),
            "Window does not hold {item}:\n{window}"
        );
    }

    let (vx, vy) = menu_at(&plain, "View");
    let view = dump(&[&scene(), "--click", &format!("{vx},{vy}")]);
    for item in ["Outliner", "Inspector", "Scene text"] {
        assert_eq!(
            count(&view, item),
            count(&plain, item),
            "View still offers {item}:\n{view}"
        );
    }
    // And View opened, or the three assertions above are about a menu that never appeared.
    assert!(
        count(&view, "Fit view") > count(&plain, "Fit view"),
        "the View menu did not open, so this proves nothing:\n{view}"
    );
}

#[test]
fn solo_is_offered_once_and_where_the_selection_is_made() {
    // **The eighth duplicate.** Seven came off the toolbar because a menu already held them; this
    // one went the other way. `solo` acts on the *selection*, the outliner is where a selection is
    // made, and its header carries the checkbox with the state visible — so the View menu's copy
    // was the one to go.
    let plain = dump(&[&scene()]);
    assert!(
        plain.contains("solo"),
        "the outliner's own checkbox is gone:\n{plain}"
    );
    let (x, y) = menu_at(&plain, "View");
    let opened = dump(&[&scene(), "--click", &format!("{x},{y}")]);
    assert!(
        !opened.contains("Solo the selection"),
        "the menu still offers it:\n{opened}"
    );
    // The menu did open — otherwise this would pass with the item still there.
    assert!(
        opened.contains("Isosurface at a value"),
        "the View menu did not open, so this proves nothing:\n{opened}"
    );
}

#[test]
fn a_picture_filtered_by_a_hidden_control_says_so() {
    // **What removing the menu copy would otherwise have left.** The checkbox is in the outliner's
    // header, so hiding the outliner — narrowing the window until it is squeezed out, or turning
    // it off in the View menu — leaves the viewport drawing one domain out of five with nothing on
    // screen saying why. A picture filtered by a control nobody can see is the silence this
    // repository hunts.
    //
    // 460 points is the first width the outliner does not survive; 1500 is one where it does.
    let hidden = dump(&[&scene(), "--solo", "--width", "420"]);
    assert!(
        hidden.contains("solo — drawing only the selection"),
        "the outliner is gone and nothing says the view is filtered:\n{hidden}"
    );

    // And it is not said when the checkbox is on screen to say it, or when solo is off — a notice
    // that is always there is decoration rather than a measurement.
    let shown = dump(&[&scene(), "--solo", "--width", "1500"]);
    assert!(
        !shown.contains("solo — drawing only the selection"),
        "said while the checkbox is on screen:\n{shown}"
    );
    let off = dump(&[&scene(), "--width", "420"]);
    assert!(
        !off.contains("solo — drawing only the selection"),
        "said while solo is off:\n{off}"
    );
}

#[test]
fn the_editor_says_what_it_does_when_nobody_is_looking() {
    // Watching a file and running on change happen with no click, so the status bar is where a
    // person finds out they are on. They were two checkboxes taking permanent room on the
    // toolbar, which showed the same state and cost four times the width.
    let d = dump(&[&scene()]);
    assert!(
        d.contains("watching the file"),
        "nothing says the editor is watching:\n{d}"
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

    // **And what the editor did is still on screen while the document is wrong.** The status bar
    // was `match error { Some => error, None => status }`, so a parse error hid every report the
    // editor made about its own actions — saving a scene with a missing field wrote the file and
    // said nothing. The two are facts about different things, and the second is wanted most
    // exactly when the first is true.
    assert!(
        d.contains("loaded "),
        "the error hid what the editor just did:\n{d}"
    );

    // Reading order: the document's problem first, then what happened. The error was on the
    // right-hand side of the window for one commit, because a right-to-left row puts the first
    // widget added at the right.
    let x = |needle: &str| -> f32 {
        let line = d
            .lines()
            .find(|l| said(l).is_some_and(|t| t.starts_with(needle)))
            .unwrap_or_else(|| panic!("no {needle} in:\n{d}"));
        box_of(line).expect("a left edge").0
    };
    assert!(
        x("1:100") < x("loaded "),
        "the error is not the first thing on the line:\n{d}"
    );
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
                (y == top).then(|| said(l).unwrap_or("").to_string())
            })
            .collect()
    };
    let menu = row(4);
    let bar = row(26);
    assert_eq!(
        menu,
        ["File", "Edit", "View", "Window", "Domain", "Run", "Watch"],
        "the menu row is not what it was"
    );
    assert!(
        bar.contains(&"run".to_string()) && bar.contains(&"fit view".to_string()),
        "the toolbar row is not below the menu row: {bar:?}"
    );
    assert_ne!(menu, bar, "the two rows collapsed into one");
}
