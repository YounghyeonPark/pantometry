//! **The asset panel, and the one thing it is for: putting a file into the scene.**
//!
//! An asset browser that lists files is a list. What makes it a browser is that a row becomes a
//! domain, and this drives that through the interface rather than through `editor_core` —
//! `a_dropped_file_becomes_a_domain` already holds the arithmetic, and none of it would matter if
//! the panel drew nothing or the button called nobody.
//!
//! # The drag and the button, both
//!
//! The row is draggable, which is what an asset browser is for, and it is also a button; both call
//! the same `App::add_asset`. A drag is a press, a movement and a release, and `--ui-dump` could
//! only drive one point at a time — so for the first version of this file the drag was code
//! nothing here could reach, and the header said so. `--drag` walks a press across frames to a
//! release now, and the drop is checked the way the button is: the scene text, the outliner and
//! the build's own note about the cells the part became. And against a control, because a drop
//! that fired wherever the pointer was let go would pass the first test too.
#![cfg(not(target_family = "wasm"))]

/// One frame of the editor, as text. The pattern `a_frame_of_the_editor_can_be_read` set, and the
/// `PANTOMETRY_RECENT` for the reason it gives: without it the start screen reads whichever list
/// belongs to whoever is running the suite.
fn dump(args: &[&str]) -> String {
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
            .env("PANTOMETRY_RECENT", none.to_string_lossy().into_owned())
            .output()
            .expect("the binary runs");
    assert!(
        out.status.success(),
        "--ui-dump {argv:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// `--layout-at`, which is a subcommand of its own and not a flag of `--ui-dump`.
fn layout_at(width: u32) -> String {
    let mut p = std::env::current_exe().expect("the test binary knows where it is");
    p.pop();
    if p.ends_with("deps") {
        p.pop();
    }
    let out =
        std::process::Command::new(p.join(format!("pantometry{}", std::env::consts::EXE_SUFFIX)))
            .args(["--layout-at", &width.to_string()])
            .output()
            .expect("the binary runs");
    assert!(out.status.success(), "--layout-at {width} failed");
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// One drawn string: where it is and what it says.
///
/// The dump writes two flag columns, then `x`, `y`, `width`, then the text — which may itself
/// contain spaces, so the text is what is left after the three numbers rather than the last field.
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
        let (Some(x), Some(y), Some(_w)) = (it.next(), it.next(), it.next()) else {
            continue;
        };
        let (Ok(x), Ok(y)) = (x.parse::<f32>(), y.parse::<f32>()) else {
            continue;
        };
        // Past the three numbers and the two spaces the format leaves before the text.
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

/// The shipped scene with a part in it, which is also the folder the assets live beside.
fn scene() -> String {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../pantometry-world/scenes/29-a-designed-bracket-becomes-cells.json")
        .to_string_lossy()
        .into_owned()
}

/// Every asset that ships beside the scenes, which is what the panel has to show.
const SHIPPED: [&str; 7] = [
    "parts/heat-sink.stl",
    "parts/l-bracket.stl",
    "parts/pipe.stl",
    "parts/plate.stl",
    "parts/rod.stl",
    "parts/wedge.stl",
    "structures/1CRN.pdb",
];

/// **The panel lists every file beside the scene that a domain could read, and no others.**
///
/// Both halves matter. A missing row is a part nobody can reach from the interface; an extra row
/// is worse — a file the drop refuses, drawn to look exactly like the six that work.
#[test]
fn the_panel_lists_every_asset_beside_the_scene() {
    let text = dump(&[&scene(), "--width", "1400", "--height", "800"]);
    let rows = drawn(&text);
    assert!(
        rows.iter().any(|d| d.text == "Assets"),
        "the panel's own heading is missing, so what follows is not evidence of anything"
    );
    for want in SHIPPED {
        assert!(
            rows.iter().any(|d| d.text == want),
            "{want} ships beside the scenes and the panel does not list it"
        );
    }
    // And nothing else in the panel's column looks like a file. The column is found from the
    // heading rather than written down, so a panel that moves does not quietly stop being checked.
    let left = rows
        .iter()
        .find(|d| d.text == "Assets")
        .map(|d| d.x)
        .expect("the heading was just asserted");
    let listed: Vec<&str> = rows
        .iter()
        .filter(|d| d.x >= left && (d.text.contains(".stl") || d.text.contains(".pdb")))
        .filter(|d| !d.text.contains(' '))
        .map(|d| d.text.as_str())
        .collect();
    assert_eq!(
        listed.len(),
        SHIPPED.len(),
        "the panel lists {listed:?} and the folder holds {SHIPPED:?}"
    );
}

/// **Every listed asset has a button beside it**, so none of them is a row that only looks like
/// the others.
#[test]
fn every_listed_asset_has_a_way_to_add_it() {
    let text = dump(&[&scene(), "--width", "1400", "--height", "800"]);
    let rows = drawn(&text);
    for want in SHIPPED {
        let row = rows
            .iter()
            .find(|d| d.text == want)
            .unwrap_or_else(|| panic!("{want} is listed"));
        let button = rows
            .iter()
            .find(|d| d.text == "+" && (d.y - row.y).abs() < 4.0 && d.x < row.x);
        assert!(
            button.is_some(),
            "{want} is drawn at y={} with no + to its left, so a reader has only the drag",
            row.y
        );
    }
}

/// **Pressing the button puts the domain in the scene**, on a grid fitted to the part.
///
/// Three independent things are checked because each could be true without the others: the scene
/// *text* gained a domain, the *outliner* shows it, and the *build* voxelised the file. A splice
/// that produced text nothing could build would pass the first alone.
#[test]
fn pressing_add_puts_the_domain_in_the_scene() {
    let before = dump(&[&scene(), "--width", "1400", "--height", "800"]);
    let rows = drawn(&before);
    let row = rows
        .iter()
        .find(|d| d.text == "parts/rod.stl")
        .expect("the rod is listed");
    let button = rows
        .iter()
        .find(|d| d.text == "+" && (d.y - row.y).abs() < 4.0 && d.x < row.x)
        .expect("the rod has a button");
    // A few points into it: the dump reports a string's own corner and the button's rect is
    // larger than the glyph it holds.
    let at = format!("{},{}", button.x + 3.0, button.y + 3.0);

    let after = dump(&[
        &scene(),
        "--width",
        "1400",
        "--height",
        "800",
        "--click",
        &at,
    ]);
    assert!(
        !before.contains("\"name\": \"rod\""),
        "the scene already had a rod, so this proves nothing"
    );
    assert!(
        after.contains("\\\"name\\\": \\\"rod\\\"") || after.contains("\"name\": \"rod\""),
        "the scene text should have gained a rod domain:\n{after}"
    );
    assert!(
        after.contains("added parts/rod.stl"),
        "the status bar should say what it did"
    );
    // The build ran on it: this line is the rasterisation report, which only exists if the file
    // was found, read, voxelised and kept.
    assert!(
        after.contains("rod/parts[0]: parts/rod.stl filled"),
        "the dropped part should have become cells:\n{after}"
    );
    let after_rows = drawn(&after);
    assert!(
        after_rows.iter().any(|d| d.text == "rod"),
        "the outliner should show the new domain"
    );
}

/// **The assets go first when the window narrows**, and the viewport figure counts them.
///
/// `--layout-at` wrote the three panel widths out as literals — a second copy of the constants —
/// so a fourth panel would have left it reporting a viewport 190 points wider than the one the
/// program lays out, with `the_viewport_always_has_room` still passing. It reads the constants now
/// and this is what holds that.
#[test]
fn the_assets_go_first_when_the_window_narrows() {
    let wide = layout_at(1400);
    println!("  {}", wide.trim());
    assert!(
        wide.contains("assets=on") && wide.contains("inspector=on"),
        "at 1400 points there is room for all four: {wide}"
    );

    // **The narrowest window each panel survives.** The assets go first as the window shrinks, so
    // they come back *last* as it grows: the width at which they appear must be the higher of the
    // two. Scanned rather than written down, because the thresholds are the sum of four constants
    // and a gap, and a number copied here would be a fifth copy of that arithmetic.
    let mut assets_from = None;
    let mut inspector_from = None;
    for width in (500..=1500).step_by(20) {
        let line = layout_at(width);
        if assets_from.is_none() && line.contains("assets=on") {
            assets_from = Some(width);
        }
        if inspector_from.is_none() && line.contains("inspector=on") {
            inspector_from = Some(width);
        }
    }
    let (a, i) = (
        assets_from.expect("the assets fit somewhere below 1500"),
        inspector_from.expect("the inspector fits somewhere below 1500"),
    );
    println!("  assets appear at {a} points, the inspector at {i}");
    assert!(
        a > i,
        "the assets are supposed to be the first panel to go and so the last to come back, and \
         they appear at {a} points against the inspector's {i}"
    );

    // And the viewport keeps its floor with four panels up, which is the whole reason the order
    // exists. Read from the line rather than recomputed, because recomputing it here would be the
    // third copy of the arithmetic.
    let view: f32 = wide
        .split_once("view=")
        .and_then(|(_, rest)| rest.split_whitespace().next())
        .and_then(|n| n.parse().ok())
        .expect("the line reports a viewport width");
    assert!(
        view >= 280.0,
        "the viewport is {view} points wide with four panels up, and its floor is 280"
    );
}

/// **A folder with nothing in it says where it looked**, rather than drawing an empty panel.
///
/// An empty list and a broken scan look identical, and this workspace has shipped the second
/// wearing the first more than once.
#[test]
fn a_folder_with_no_assets_says_where_it_looked() {
    let dir = std::env::temp_dir().join(format!("pantometry-bare-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("somewhere to put a scene with nothing beside it");
    let path = dir.join("bare.json");
    std::fs::write(
        &path,
        r#"{ "title": "nothing beside me", "duration_s": 1.0, "frames": 2, "domains": [] }"#,
    )
    .expect("the scene is written");

    let text = dump(&[
        &path.to_string_lossy(),
        "--width",
        "1400",
        "--height",
        "800",
    ]);
    let rows = drawn(&text);
    assert!(
        rows.iter().any(|d| d.text == "Assets"),
        "the panel is still there when it has nothing to show"
    );
    let note = rows
        .iter()
        .find(|d| d.text.starts_with("no ."))
        .unwrap_or_else(|| panic!("the empty panel says nothing at all:\n{text}"));
    println!("  {}", note.text);
    assert!(
        note.text.contains(".stl") && note.text.contains(".pdb"),
        "it should say what it was looking for: {}",
        note.text
    );
    assert!(
        note.text.contains("bare") || note.text.contains("pantometry-bare"),
        "it should say where it looked: {}",
        note.text
    );

    let _ = std::fs::remove_dir_all(&dir);
}

/// The viewport's rect from the dump's own `viewport=x,y WxH` line, as its centre.
fn viewport_centre(dump: &str) -> (f32, f32) {
    let line = dump
        .lines()
        .find_map(|l| l.strip_prefix("viewport="))
        .unwrap_or_else(|| panic!("the dump reports no viewport:\n{dump}"));
    let (at, size) = line.split_once(' ').expect("x,y WxH");
    let (x, y) = at.split_once(',').expect("x,y");
    let (w, h) = size.split_once('x').expect("WxH");
    let n = |v: &str| v.trim().parse::<f32>().expect("a number");
    (n(x) + n(w) / 2.0, n(y) + n(h) / 2.0)
}

/// Where the rod's row is drawn, a little inside its label: on the label and not on its `+`, so a
/// release that landed as a click would press nothing.
fn rod_row(dump: &str) -> (f32, f32) {
    let rows = drawn(dump);
    let row = rows
        .iter()
        .find(|d| d.text == "parts/rod.stl")
        .expect("the rod is listed");
    (row.x + 20.0, row.y + 6.0)
}

/// **Dragging a row onto the view puts the domain in the scene**, the same three ways the button does.
#[test]
fn dragging_a_row_onto_the_view_puts_the_domain_in_the_scene() {
    let before = dump(&[&scene(), "--width", "1400", "--height", "800"]);
    let (from, to) = (rod_row(&before), viewport_centre(&before));
    let drag = format!("{},{},{},{}", from.0, from.1, to.0, to.1);
    println!("  dragged from {from:?} to {to:?}");

    let after = dump(&[
        &scene(),
        "--width",
        "1400",
        "--height",
        "800",
        "--drag",
        &drag,
    ]);
    assert!(
        !before.contains("added parts/rod.stl"),
        "the status said it before anything was dragged, so it proves nothing"
    );
    assert!(
        after.contains("added parts/rod.stl"),
        "the status bar should say what the drop did:\n{after}"
    );
    assert!(
        after.contains("rod/parts[0]: parts/rod.stl filled"),
        "the dropped part should have become cells:\n{after}"
    );
    assert!(
        drawn(&after).iter().any(|d| d.text == "rod"),
        "the outliner should show the new domain"
    );
}

/// **A drag let go anywhere but the view adds nothing.**
///
/// The control for the test above. The drop zone is the viewport, and a drag released over the
/// outliner — the same press, the same walk, a different place to let go — has to leave the scene
/// as it was. Without this, a drop that fired on any release would pass above just as well.
#[test]
fn a_drag_let_go_anywhere_but_the_view_adds_nothing() {
    let before = dump(&[&scene(), "--width", "1400", "--height", "800"]);
    let from = rod_row(&before);
    let outliner = drawn(&before)
        .iter()
        .find(|d| d.text == "Outliner")
        .map(|d| (d.x + 20.0, d.y + 120.0))
        .expect("the outliner is on screen at 1400");
    let drag = format!("{},{},{},{}", from.0, from.1, outliner.0, outliner.1);

    let after = dump(&[
        &scene(),
        "--width",
        "1400",
        "--height",
        "800",
        "--drag",
        &drag,
    ]);
    assert!(
        !after.contains("added parts/rod.stl") && !after.contains("rod/parts[0]"),
        "a drag let go over the outliner added the rod anyway:\n{after}"
    );
}
