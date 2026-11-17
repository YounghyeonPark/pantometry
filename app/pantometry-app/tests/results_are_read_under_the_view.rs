//! **The results panel, read off the frame.**
//!
//! `editor_core::post` is held to closed forms in its own tests. None of that matters if the
//! panel draws a number `post` did not compute, opens on a tab with nothing in it, or lays its
//! numbers over each other — so this reads the panel through `--ui-dump`, the way a person reads
//! it, and checks what it says against things that do not come from the panel: the block's own
//! reported peak, the size its scene states, and where its scene put the heat.

#![cfg(not(target_family = "wasm"))]

/// **In a view too short for its bands, a domain's name is still on the view.** The corner labels
/// are held above the bands the readouts own, and in a view 214 points tall with an isosurface
/// those bands reach above its top: every name was held there and clipped away. Overlaps are
/// allowed here — a view this short is crowded — and absence is not.
#[test]
fn a_short_view_keeps_its_names_on_it() {
    let d = dump(&[&scene(BLOCK), "--ran", "--iso", "--height", "300"]);
    let view = d
        .lines()
        .find_map(|l| l.strip_prefix("viewport="))
        .unwrap_or_else(|| panic!("no viewport:\n{d}"));
    let (at, size) = view.split_once(' ').expect("x,y wxh");
    let (x, y) = at.split_once(',').expect("x,y");
    let (w, _) = size.split_once('x').expect("wxh");
    let (x, y, w): (f32, f32, f32) = (
        x.parse().expect("x"),
        y.parse().expect("y"),
        w.parse().expect("w"),
    );
    let name = drawn(&d)
        .into_iter()
        .find(|r| r.text == "block" && r.x >= x && r.x <= x + w)
        .unwrap_or_else(|| panic!("the view has no name:\n{d}"));
    assert!(
        name.y >= y,
        "the name is at y={}, above the view's top at {y}:\n{d}",
        name.y
    );
}

fn dump(args: &[&str]) -> String {
    let out = run(args);
    assert!(
        out.status.success(),
        "--ui-dump {args:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

fn run(args: &[&str]) -> std::process::Output {
    let mut argv = vec!["--ui-dump".to_string()];
    argv.extend(args.iter().map(|a| (*a).to_string()));
    let none = std::env::temp_dir().join("pantometry-no-such-recent-list.json");
    std::process::Command::new(env!("CARGO_BIN_EXE_pantometry"))
        .args(&argv)
        .env("PANTOMETRY_RECENT", none.to_string_lossy().into_owned())
        .output()
        .expect("the binary runs")
}

/// One drawn string: where it is and what it says, from a dump line's two flag columns, `x`, `y`,
/// `width` and the text — the reading `a_scene_is_set_up_from_its_tree` uses.
struct Drawn {
    x: f32,
    y: f32,
    text: String,
}

fn drawn(dump: &str) -> Vec<Drawn> {
    let mut out = Vec::new();
    for line in dump.lines() {
        if line.len() < 3 || !line.is_char_boundary(2) {
            continue;
        }
        let rest = &line[2..];
        let mut it = rest.split_whitespace();
        let (Some(x), Some(y), Some(w)) = (it.next(), it.next(), it.next()) else {
            continue;
        };
        let (Ok(x), Ok(y), Ok(_)) = (x.parse::<f32>(), y.parse::<f32>(), w.parse::<f32>()) else {
            continue;
        };
        let after = rest
            .split_whitespace()
            .take(3)
            .fold(0usize, |at, f| rest[at..].find(f).unwrap() + at + f.len());
        out.push(Drawn {
            x,
            y,
            text: rest[after..].trim_start().to_string(),
        });
    }
    out
}

/// The panel's own strings: everything drawn at or below its `Results` heading and not left of
/// it. The outliner runs the window's full height beside the panel, and on a short window its
/// `mean` or `peak` rows would sit below the heading and be found first.
fn in_panel(dump: &str) -> Vec<Drawn> {
    let rows = drawn(dump);
    let Some((top, left)) = rows
        .iter()
        .find(|d| d.text == "Results")
        .map(|d| (d.y, d.x))
    else {
        panic!("no results panel:\n{dump}");
    };
    rows.into_iter()
        .filter(|d| d.y >= top - 1.0 && d.x >= left - 1.0)
        .collect()
}

/// What is beside the first row labelled `label` below the row labelled `after`: the `at` under
/// `max` and not the `at` under `min`, which says the same kind of thing about another place.
fn beside_after(rows: &[Drawn], after: &str, label: &str) -> Vec<String> {
    let Some(start) = rows.iter().find(|d| d.text == after).map(|d| d.y) else {
        panic!("nothing in the panel says {after:?}");
    };
    let Some(key) = rows
        .iter()
        .filter(|d| d.text == label && d.y > start + 1.0)
        .min_by(|a, b| a.y.total_cmp(&b.y))
    else {
        panic!("no {label:?} under {after:?}");
    };
    rows.iter()
        .filter(|d| (d.y - key.y).abs() < 6.0 && d.x > key.x && d.x < key.x + 320.0)
        .map(|d| d.text.clone())
        .collect()
}

/// What a tab says in place of a chart it has no room for.
const NO_ROOM: &str = "the chart needs a wider or taller panel";

fn says(rows: &[Drawn], text: &str) -> bool {
    rows.iter().any(|d| d.text == text)
}

/// What is drawn to the right of `label` on its line, in the statistics' column.
fn beside(rows: &[Drawn], label: &str) -> Vec<String> {
    let Some(key) = rows.iter().find(|d| d.text == label) else {
        panic!("nothing in the panel says {label:?}");
    };
    rows.iter()
        .filter(|d| (d.y - key.y).abs() < 6.0 && d.x > key.x && d.x < key.x + 320.0)
        .map(|d| d.text.clone())
        .collect()
}

fn header(dump: &str, key: &str) -> usize {
    dump.lines()
        .find_map(|l| l.strip_prefix(&format!("{key}=")))
        .and_then(|v| v.parse().ok())
        .unwrap_or_else(|| panic!("no {key}= line:\n{dump}"))
}

fn scene(name: &str) -> String {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../pantometry-world/scenes")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

const BLOCK: &str = "15-a-hot-spot-in-a-block.json";

/// **The region's numbers are the block's.** Scene 15 is a 9 x 9 x 9 block of 1 mm cells, heated
/// at the middle cell and cooled nowhere. Its volume is (9 mm)^3, its mean is fixed by its initial
/// condition at 293.15 + 60/729 K for the whole run, the hottest place is the cell the heat was
/// put in, and the domain reports its own peak in Celsius every frame — which the field's
/// maximum has to be, read off a different path.
#[test]
fn a_field_opens_on_its_region_and_the_numbers_are_the_blocks() {
    let d = dump(&[&scene(BLOCK), "--ran"]);
    assert_eq!(header(&d, "overlap"), 0, "{d}");
    let rows = in_panel(&d);
    assert!(
        says(&rows, "the whole domain"),
        "it opens on the region:\n{d}"
    );

    assert_eq!(beside(&rows, "volume"), ["7.2900e-7 m^3"]);
    assert_eq!(beside(&rows, "cells"), ["729"]);
    let mean = 293.15 + 60.0 / 729.0;
    assert_eq!(beside(&rows, "mean"), [format!("{mean:.2} K")]);

    // The block's own reading, on the same frame, from the readout the viewport draws.
    let peak: f64 = drawn(&d)
        .iter()
        .find_map(|r| {
            r.text
                .strip_prefix("block peak ")?
                .strip_suffix(" C")?
                .parse()
                .ok()
        })
        .unwrap_or_else(|| panic!("no peak reading:\n{d}"));
    assert_eq!(beside(&rows, "max"), [format!("{:.2} K", peak + 273.15)]);
    // Where the heat was put, and where it stays by symmetry: on the row under `max`, since
    // the row under `min` names a place in the same words.
    assert_eq!(beside_after(&rows, "max", "at"), ["(4.5, 4.5, 4.5) mm"]);
    assert_eq!(beside_after(&rows, "min", "at"), ["(0.5, 0.5, 0.5) mm"]);
    assert!(says(&rows, "block max") && says(&rows, "block mean") && says(&rows, "block min"));
    assert!(says(&rows, "time / s"), "and the three over the run:\n{d}");
}

/// **A part with holes is counted by its cells, and agrees with its own readings.** Scene 29's
/// bracket is an STL filled into 2 mm cells, and the check's note says how many it filled; the
/// rest of its grid is holes. The block reports its own peak, coldest and mean in Celsius every
/// frame, computed by the domain over the cells it has — a different path to the same three
/// numbers the panel reads off the field.
#[test]
fn a_part_with_holes_agrees_with_its_own_readings() {
    let d = dump(&[&scene("29-a-designed-bracket-becomes-cells.json"), "--ran"]);
    let rows = in_panel(&d);
    let filled: usize = d
        .split("filled ")
        .nth(1)
        .and_then(|rest| rest.split_whitespace().next())
        .and_then(|n| n.parse().ok())
        .unwrap_or_else(|| panic!("no note of the cells filled:\n{d}"));
    let cells = beside(&rows, "cells");
    assert_eq!(cells.len(), 1);
    assert!(
        cells[0].starts_with(&format!("{filled} and ")) && cells[0].ends_with(" not in the domain"),
        "{cells:?} against {filled} filled"
    );
    let cell = 0.002f64 * 0.002 * 0.002;
    let volume = format!("{} m^3", editor_core::magnitude(filled as f64 * cell));
    assert_eq!(beside(&rows, "volume"), [volume]);

    let reading = |label: &str| -> f64 {
        drawn(&d)
            .iter()
            .find_map(|r| {
                r.text
                    .strip_prefix(&format!("bracket {label} "))?
                    .strip_suffix(" C")?
                    .parse()
                    .ok()
            })
            .unwrap_or_else(|| panic!("no {label} reading:\n{d}"))
    };
    for (row, label) in [("max", "peak"), ("min", "coldest"), ("mean", "mean")] {
        let kelvin = reading(label) + 273.15;
        assert_eq!(
            beside(&rows, row),
            [format!("{} K", editor_core::magnitude(kelvin))],
            "{row} against the bracket's own {label}"
        );
    }
}

/// **Each tab draws what it reads, and nothing on top of anything.**
#[test]
fn every_tab_draws_its_own_answer() {
    let tab = |t: &str| {
        let d = dump(&[&scene(BLOCK), "--ran", "--results", t]);
        assert_eq!(
            (header(&d, "overlap"), header(&d, "cut")),
            (0, 0),
            "{t}:\n{d}"
        );
        (in_panel(&d), d)
    };

    let (rows, d) = tab("slice");
    assert!(says(&rows, "x 0 to 9 across, y 0 to 9 up, mm"), "{d}");
    // The colour bar's ends are the run's range: where the block began and the hot cell began.
    assert!(says(&rows, "353.15 K") && says(&rows, "293.15 K"), "{d}");

    let (rows, d) = tab("line");
    assert!(
        says(&rows, "distance / mm") && says(&rows, "block along the line"),
        "{d}"
    );
    // Through the middle along x, wall to wall, read off the plot's own axis: its labels are
    // fourteen points above its title. `0` and `9` alone were the line's two ends in the
    // controls, which say where it was asked to go and not what was drawn.
    let title = rows
        .iter()
        .find(|r| r.text == "distance / mm")
        .expect("the distance axis");
    let mut ticks: Vec<(f32, String)> = rows
        .iter()
        .filter(|r| (r.y - (title.y - 14.0)).abs() < 3.0)
        .map(|r| (r.x, r.text.clone()))
        .collect();
    ticks.sort_by(|a, b| a.0.total_cmp(&b.0));
    let ticks: Vec<String> = ticks.into_iter().map(|t| t.1).collect();
    assert_eq!(ticks, ["0", "2", "4", "6", "8"], "{d}");

    let (rows, d) = tab("point");
    // The probe opens in the middle, which is the hottest cell, so it reads the region's maximum.
    let (region, _) = tab("region");
    let max = beside(&region, "max");
    assert!(says(&rows, &max[0]), "the middle reads {max:?}:\n{d}");
    assert!(
        says(&rows, "block at the point") && says(&rows, "time / s"),
        "{d}"
    );

    let (rows, d) = tab("readings");
    assert!(
        says(&rows, "block peak") && says(&rows, "block absorbed"),
        "{d}"
    );
    assert!(!says(&rows, "the whole domain"), "{d}");
    // And the reading picked, charted: four rows of a table left a row too short for a chart
    // beside it, and the tab was a table over an empty strip.
    assert!(says(&rows, "time / s"), "the readings chart:\n{d}");
}

/// **A run with nothing to place opens on what it reports.** Scene 11 is a thermal network:
/// nodes and conductances, no grid. The region, the slice, the line and the probe each read a
/// field, and asked of this one they say there is none rather than drawing an empty plot.
#[test]
fn a_run_with_no_field_opens_on_its_readings() {
    let network = scene("11-motor-thermal-network.json");
    let d = dump(&[&network, "--ran"]);
    let rows = in_panel(&d);
    assert!(!says(&rows, "the whole domain"), "{d}");
    assert!(
        rows.iter().any(|r| r.text.starts_with("time / ")),
        "the readings:\n{d}"
    );

    let empty = dump(&[&network, "--ran", "--results", "region"]);
    assert!(
        empty.contains("nothing in this run is a field"),
        "the region says why it is empty:\n{empty}"
    );
    // And says it in a panel as tall as the others: a panel keeps the size of what is drawn in
    // it, and one sentence is a strip.
    let heading = |d: &str| {
        drawn(d)
            .into_iter()
            .find(|r| r.text == "Results")
            .map(|r| r.y)
            .unwrap_or_else(|| panic!("no panel:\n{d}"))
    };
    assert_eq!(heading(&empty), heading(&d), "{empty}");
}

#[test]
fn there_are_no_results_before_a_run() {
    let d = dump(&[&scene(BLOCK)]);
    assert!(!drawn(&d).iter().any(|r| r.text == "Results"), "{d}");
}

/// **A tab of a panel that is not there is refused.** Without a run there is no panel, and a frame
/// with none drawn and exit 0 is a test that forgot `--ran` passing over nothing.
#[test]
fn a_tab_without_a_run_is_refused() {
    let out = run(&[&scene(BLOCK), "--results", "line"]);
    assert!(!out.status.success());
    let said = String::from_utf8_lossy(&out.stderr);
    assert!(
        said.contains("--results takes effect only with --ran or --run-after-clicks"),
        "{said}"
    );
}

#[test]
fn a_tab_that_is_not_there_is_refused() {
    let out = run(&[&scene(BLOCK), "--ran", "--results", "contour"]);
    assert!(!out.status.success());
    let said = String::from_utf8_lossy(&out.stderr);
    assert!(
        said.contains("--results takes one of Region, Slice, Line, Point, Readings; got contour"),
        "{said}"
    );
}

/// **The domain asked about follows the outliner.** Scene 25 has two fields, a block's
/// temperature and a structure's stress; selecting the structure's row turns the panel to it.
#[test]
fn the_panel_reads_the_domain_the_outliner_selects() {
    let two = scene("25-what-140-kelvin-does-to-the-solder.json");
    let first = in_panel(&dump(&[&two, "--ran"]));
    let picked = in_panel(&dump(&[&two, "--ran", "--select", "/domains/stress"]));
    let menu = |rows: &[Drawn]| {
        let tabs = rows
            .iter()
            .find(|r| r.text == "Readings")
            .expect("the tabs");
        rows.iter()
            .filter(|r| (r.y - tabs.y).abs() < 6.0 && r.x > tabs.x)
            .map(|r| r.text.clone())
            .collect::<Vec<_>>()
    };
    assert_eq!(menu(&first), ["module"]);
    assert_eq!(menu(&picked), ["stress"]);

    // **And a pick in the menu stands.** The selection was applied every frame, so with `module`
    // selected, choosing `stress` showed it for one frame and went back.
    let selected = ["--ran", "--select", "/domains/module"];
    let base = |extra: &[&str]| {
        let mut args = vec![two.as_str()];
        args.extend_from_slice(&selected);
        args.extend_from_slice(extra);
        dump(&args)
    };
    let shut = base(&[]);
    let combo = in_panel(&shut)
        .into_iter()
        .find(|r| r.text == "module")
        .map(|r| format!("{},{}", r.x + 3.0, r.y + 3.0))
        .expect("the menu");
    let open = base(&["--click", &combo]);
    // The menu's item, which is in the panel. The lowest `stress` on screen was the outliner's
    // readings row, and clicking that selected the structure through the outliner: the test
    // passed with the menu's pick undone every frame, which is what it was written to see.
    let pick = in_panel(&open)
        .into_iter()
        .find(|r| r.text == "stress")
        .map(|r| format!("{},{}", r.x + 3.0, r.y + 3.0))
        .unwrap_or_else(|| panic!("the open menu has no stress:\n{open}"));
    let chose = base(&["--click", &combo, "--click", &pick]);
    assert_eq!(menu(&in_panel(&chose)), ["stress"], "{chose}");

    // **Two domains on one corner have two names.** The structure follows its block on the
    // block's own grid, so both corners are one point, and their names were drawn on it.
    //
    // **And a selected domain's inspector fits its panel.** Its bounds line did not wrap and ran
    // past the window, and the viewport beside it ran under the inspector: every selected domain
    // read `cut=1`, and this frame, with a run, `overlap=15`. It was the first frame here with a
    // domain selected and a run on screen, which is how it was found.
    for d in [dump(&[&two, "--ran"]), chose] {
        assert_eq!((header(&d, "overlap"), header(&d, "cut")), (0, 0), "{d}");
    }
}

/// **A placed domain's positions say whose frame they are in.** Scene 30's busbars are turned and
/// moved, so a probe at `(50, 20, 10) mm` is a place in the busbars and not in the world.
#[test]
fn a_placed_domain_says_its_positions_are_its_own() {
    let placed = dump(&[
        &scene("30-two-phases-crossing-at-a-clearance.json"),
        "--ran",
    ]);
    assert!(
        placed
            .contains("positions are in the frame of busbars, before its pose moves and turns it"),
        "{placed}"
    );
    let here = dump(&[&scene(BLOCK), "--ran"]);
    assert!(!here.contains("positions are in the frame of"), "{here}");
}

/// **A plane along a bar meets it at one sample, and says so.** Scene 04's bar is flat across, so
/// a plane across x has nothing in it but the point where the bar is cut, and it drew nothing.
#[test]
fn a_plane_along_a_bar_says_its_one_value() {
    let bar = scene("04-heater-and-bar.json");
    let shut = dump(&[&bar, "--ran", "--results", "slice"]);
    let axis = in_panel(&shut)
        .into_iter()
        .find(|r| r.text == "z")
        .map(|r| format!("{},{}", r.x + 3.0, r.y + 3.0))
        .unwrap_or_else(|| panic!("no axis menu:\n{shut}"));
    let open = dump(&[&bar, "--ran", "--results", "slice", "--click", &axis]);
    let x = drawn(&open)
        .into_iter()
        .filter(|r| r.text == "x")
        .max_by(|a, b| a.y.total_cmp(&b.y))
        .map(|r| format!("{},{}", r.x + 3.0, r.y + 3.0))
        .unwrap_or_else(|| panic!("the open menu has no x:\n{open}"));
    let chose = dump(&[
        &bar,
        "--ran",
        "--results",
        "slice",
        "--click",
        &axis,
        "--click",
        &x,
    ]);
    assert!(
        chose.contains("the plane meets bar at one sample: "),
        "{chose}"
    );
}

/// **Window, Results puts it away**, the way the other panels go.
#[test]
fn the_window_menu_puts_the_results_away() {
    let plain = dump(&[&scene(BLOCK), "--ran"]);
    let at = |d: &str, text: &str, below: f32| {
        drawn(d)
            .into_iter()
            .filter(|r| r.text == text && r.y > below)
            .min_by(|a, b| a.y.total_cmp(&b.y))
            .map(|r| format!("{},{}", r.x + 3.0, r.y + 3.0))
            .unwrap_or_else(|| panic!("nothing says {text:?}:\n{d}"))
    };
    let window = at(&plain, "Window", -1.0);
    let open = dump(&[&scene(BLOCK), "--ran", "--click", &window]);
    // The menu's item, which is above the panel's own heading.
    let item = at(&open, "Results", 20.0);
    let shut = dump(&[&scene(BLOCK), "--ran", "--click", &window, "--click", &item]);
    assert!(!says(&drawn(&shut), "the whole domain"), "{shut}");
    assert!(says(&drawn(&plain), "the whole domain"));
}

/// **At every size the panel draws its answer, or says it has no room for it, and nothing is on
/// top of anything.** The line tab has the widest row, two corners of three coordinates and three
/// buttons; under 921 points the inspector goes and the charts beside the tables lose their room.
/// Asserting only that nothing overlapped passed over a panel that drew nothing at all.
#[test]
fn the_panel_holds_at_every_size() {
    let answer = |t: &str| -> &[&str] {
        match t {
            "region" | "point" => &["time / s", NO_ROOM],
            "line" => &["distance / mm", NO_ROOM],
            "slice" => &["x 0 to 9 across, y 0 to 9 up, mm"],
            _ => &["time / s", NO_ROOM],
        }
    };
    for (width, height) in [
        ("700", "950"),
        // About 233 points beside the region's table, under the chart's floor of 320. With the
        // floor sabotaged to 120 the chart was drawn in them and this size failed, where 700,
        // with 116 points, leaves the chart out either way.
        ("740", "950"),
        ("900", "950"),
        ("1000", "950"),
        ("1200", "950"),
        ("1500", "950"),
        ("1800", "950"),
        ("1500", "700"),
        ("1500", "600"),
    ] {
        for t in ["region", "slice", "line", "point", "readings"] {
            let d = dump(&[
                &scene(BLOCK),
                "--ran",
                "--results",
                t,
                "--width",
                width,
                "--height",
                height,
            ]);
            assert_eq!(
                (header(&d, "overlap"), header(&d, "cut")),
                (0, 0),
                "{t} at {width} x {height}:\n{d}"
            );
            let rows = in_panel(&d);
            assert!(
                answer(t).iter().any(|a| says(&rows, a)),
                "{t} at {width} x {height} drew neither {:?}:\n{d}",
                answer(t)
            );
        }
    }
    // Too short for the panel and the view both: it goes, and the status bar says why.
    let short = dump(&[&scene(BLOCK), "--ran", "--height", "400"]);
    assert!(
        !drawn(&short).iter().any(|r| r.text == "Results"),
        "{short}"
    );
    assert!(
        short.contains("the results hidden — the window is too short for them and the view"),
        "{short}"
    );
}
