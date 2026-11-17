//! **The settings pane, driven through the interface.**
//!
//! `editor_core`'s own tests hold the forms and the splices against every shipped scene. None of
//! that matters if the pane draws nothing for a row, or a button calls nobody, or a menu's pick
//! never reaches the file — so this presses the buttons through `--ui-dump`, the way
//! `an_asset_is_added_from_the_panel` presses the asset panel's, and reads what came of it from
//! three places: the tree, the status line and the scene's text.

#![cfg(not(target_family = "wasm"))]

fn dump(args: &[&str]) -> String {
    let mut argv = vec!["--ui-dump".to_string()];
    argv.extend(args.iter().map(|a| (*a).to_string()));
    let none = std::env::temp_dir().join("pantometry-no-such-recent-list.json");
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_pantometry"))
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

/// One drawn string: where it is and what it says. The same reading of a dump line as
/// `an_asset_is_added_from_the_panel`'s: two flag columns, `x`, `y`, `width`, then the text.
struct Drawn {
    x: f32,
    y: f32,
    w: f32,
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
        let (Ok(x), Ok(y), Ok(w)) = (x.parse::<f32>(), y.parse::<f32>(), w.parse::<f32>()) else {
            continue;
        };
        let after = rest
            .split_whitespace()
            .take(3)
            .fold(0usize, |at, f| rest[at..].find(f).unwrap() + at + f.len());
        out.push(Drawn {
            x,
            y,
            w,
            text: rest[after..].trim_start().to_string(),
        });
    }
    out
}

/// Where to click to press the string `text`: a few points into it, the way the asset test does.
fn at(dump: &str, text: &str) -> String {
    let rows = drawn(dump);
    let d = rows
        .iter()
        .find(|d| d.text == text)
        .unwrap_or_else(|| panic!("nothing drawn says {text:?}:\n{dump}"));
    format!("{},{}", d.x + 3.0, d.y + 3.0)
}

/// The same, for the value drawn on `label`'s line of a form. A face is drawn three times: as the
/// entry's name in the tree, as the inspector's heading and as the menu's text, and only the
/// last is on the line that says `face`.
fn beside(dump: &str, label: &str, text: &str) -> String {
    let rows = drawn(dump);
    let key = rows
        .iter()
        .find(|d| d.text == label)
        .unwrap_or_else(|| panic!("nothing drawn says {label:?}:\n{dump}"));
    let d = rows
        .iter()
        .find(|d| d.text == text && (d.y - key.y).abs() < 6.0 && d.x > key.x)
        .unwrap_or_else(|| panic!("nothing says {text:?} on the {label:?} line:\n{dump}"));
    format!("{},{}", d.x + 3.0, d.y + 3.0)
}

fn says(dump: &str, text: &str) -> bool {
    drawn(dump).iter().any(|d| d.text == text)
}

fn scene() -> String {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../pantometry-world/scenes/29-a-designed-bracket-becomes-cells.json")
        .to_string_lossy()
        .into_owned()
}

/// **The pane opens on a form, and the text is one tab away.**
///
/// With nothing selected it is the scene's own form. The scene's text is not drawn there — it is
/// on the Source tab — so the form is not a caption over a text box.
#[test]
fn the_pane_is_a_form_and_the_source_is_a_tab() {
    let plain = dump(&[&scene()]);
    assert!(
        says(&plain, "Settings") && says(&plain, "Source"),
        "two tabs:\n{plain}"
    );
    assert!(
        says(&plain, "duration_s (s)"),
        "the scene's own settings:\n{plain}"
    );
    assert!(
        says(&plain, "Cooling (2)"),
        "the tree shows the bracket's lists:\n{plain}"
    );
    assert!(
        !plain.contains("\"convection_w_per_m2_k\""),
        "the text is not on the form tab"
    );

    let source = dump(&[&scene(), "--source"]);
    assert!(
        source.contains("\"convection_w_per_m2_k\""),
        "and is on the other:\n{source}"
    );
}

/// **An entry's row opens that entry's form**, with the format's set for its face and the two
/// things an entry can have done to it.
#[test]
fn an_entry_opens_its_own_form() {
    let d = dump(&[&scene(), "--select", "/domains/bracket/cooling/0"]);
    // `z-min` is also the entry's name in the tree and the inspector's heading, so the face is
    // read off its own line of the form.
    assert_eq!(on_line(&d, "face"), ["z-min"], "{d}");
    for want in ["ambient_c (°C)", "duplicate", "delete"] {
        assert!(says(&d, want), "{want} is not drawn:\n{d}");
    }
    // The domain's own values are on the domain's form, not this one.
    assert!(
        !says(&d, "cell_mm (mm)"),
        "the entry's form holds the domain's values:\n{d}"
    );
}

/// A row that is not in the scene is refused rather than drawn as an empty pane.
#[test]
fn a_row_that_is_not_there_is_refused() {
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_pantometry"))
        .args([
            "--ui-dump",
            &scene(),
            "--select",
            "/domains/bracket/cooling/9",
        ])
        .output()
        .expect("the binary runs");
    assert!(!out.status.success(), "a selection of nothing was drawn");
    // Its own words: a wrong scene path or a parser refusing every path also exits non-zero.
    let said = String::from_utf8_lossy(&out.stderr);
    assert!(
        said.contains("--select: no row /domains/bracket/cooling/9"),
        "refused, but not for the row: {said}"
    );
}

/// **Pressing add puts an entry in the scene, and the tree, the status line and the text agree.**
///
/// Then the Source tab, clicked after the add, shows the third face in the file. Three readings
/// of one change, because each could be true without the others.
#[test]
fn pressing_add_puts_an_entry_in_the_scene() {
    let list = "/domains/bracket/cooling";
    let before = dump(&[&scene(), "--select", list]);
    let add = at(&before, "+ cooling");
    let source = at(&before, "Source");
    let after = dump(&[&scene(), "--select", list, "--click", &add]);
    assert!(
        says(&after, "Cooling (3)"),
        "the tree has three faces:\n{after}"
    );
    assert!(
        after.contains("added a cooling entry to bracket"),
        "the status line says so"
    );
    let text = dump(&[
        &scene(),
        "--select",
        list,
        "--click",
        &add,
        "--click",
        &source,
    ]);
    assert_eq!(
        text.matches("convection_w_per_m2_k").count(),
        3,
        "the file holds three faces:\n{text}"
    );
}

/// **Delete takes an entry out**, and the tree says one face is left.
#[test]
fn delete_takes_an_entry_out() {
    let entry = "/domains/bracket/cooling/1";
    let before = dump(&[&scene(), "--select", entry]);
    assert!(
        says(&before, "z-max"),
        "z-max is drawn before the delete:\n{before}"
    );
    let after = dump(&[
        &scene(),
        "--select",
        entry,
        "--click",
        &at(&before, "delete"),
    ]);
    assert!(says(&after, "Cooling (1)"), "one face is left:\n{after}");
    assert!(
        !says(&after, "z-max"),
        "and it is the one that was not deleted:\n{after}"
    );
}

/// **A menu's pick reaches the file.** The face menu is opened, `x-min` is picked from it, and
/// the Source tab shows the word written where `z-min` was.
#[test]
fn a_face_picked_from_its_menu_is_written() {
    let entry = "/domains/bracket/cooling/0";
    let before = dump(&[&scene(), "--select", entry]);
    let menu = beside(&before, "face", "z-min");
    let open = dump(&[&scene(), "--select", entry, "--click", &menu]);
    // The menu's items are drawn once it is open, `x-min` among them; it was drawn nowhere before.
    assert!(
        !says(&before, "x-min") && says(&open, "x-min"),
        "the menu did not open:\n{open}"
    );
    let pick = at(&open, "x-min");
    let after = dump(&[
        &scene(),
        "--select",
        entry,
        "--click",
        &menu,
        "--click",
        &pick,
    ]);
    let source = at(&after, "Source");
    let text = dump(&[
        &scene(),
        "--select",
        entry,
        "--click",
        &menu,
        "--click",
        &pick,
        "--click",
        &source,
    ]);
    assert!(
        text.contains("\"face\": \"x-min\""),
        "the file says x-min:\n{text}"
    );
    assert!(
        !text.contains("\"face\": \"z-min\""),
        "and z-min is gone:\n{text}"
    );
}

/// What is drawn on the line of a form that says `label`, to its right and inside the pane.
///
/// Inside, because the asset panel's rows sit at the same heights: the face line of a cooled
/// face's form shares its y with an asset's `+` and its file name, and this read all three.
/// The pane ends where the asset panel's heading begins.
fn on_line(dump: &str, label: &str) -> Vec<String> {
    let rows = drawn(dump);
    let Some(key) = rows.iter().find(|d| d.text == label) else {
        panic!("nothing drawn says {label:?}:\n{dump}");
    };
    let edge = rows
        .iter()
        .find(|d| d.text == "Assets")
        .map_or(f32::INFINITY, |d| d.x);
    rows.iter()
        .filter(|d| (d.y - key.y).abs() < 6.0 && d.x > key.x && d.x < edge)
        .map(|d| d.text.clone())
        .collect()
}

/// **A small number keeps its magnitude.** The form's drag printed two decimals, so the protein's
/// 2e-11 s duration and the bracket's 1e-8 conservation tolerance were both drawn `0.00` on the
/// front page of the editor — the figure `docs/editor-protein.png` was retaken showing it.
#[test]
fn a_small_number_is_drawn_with_its_magnitude() {
    let protein = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../pantometry-world/scenes/31-a-protein-shaking-at-body-temperature.json")
        .to_string_lossy()
        .into_owned();
    let d = dump(&[&protein]);
    assert_eq!(on_line(&d, "duration_s (s)"), ["2.0000e-11"], "{d}");

    let d = dump(&[&scene()]);
    assert_eq!(on_line(&d, "conservation_tolerance"), ["1.0000e-8"], "{d}");
    assert_eq!(on_line(&d, "duration_s (s)"), ["200.00"], "{d}");
}

/// **Placement is offered for a domain that is at somewhere, and only there.** A thermal network
/// is not, and `build` takes a pose for one and does nothing with it: dragging its position wrote
/// the file, recorded an undo step and moved nothing.
#[test]
fn placement_is_offered_only_where_it_moves_something() {
    let bracket = dump(&[&scene(), "--select", "/domains/bracket"]);
    assert!(
        says(&bracket, "at_m (m)"),
        "the bracket can be placed:\n{bracket}"
    );

    let motor = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../pantometry-world/scenes/11-motor-thermal-network.json")
        .to_string_lossy()
        .into_owned();
    let network = dump(&[&motor, "--select", "/domains/motor"]);
    assert!(
        !says(&network, "at_m (m)"),
        "a network offered a position:\n{network}"
    );
    assert!(
        network.contains("no placement: this domain is not at anywhere"),
        "and the form says why:\n{network}"
    );
}

/// **A long value is read whole, inside its own pane.** A single-line field lays its string out
/// in full and clips it to the box, so the protein's title showed its first third, and the rest
/// lay under the asset panel's rows at the same height -- `docs/editor-protein.png` was retaken
/// with `overlap=2` from exactly that. A value on a line of its own wraps to its field now.
#[test]
fn a_long_value_wraps_inside_its_pane() {
    let protein = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../pantometry-world/scenes/31-a-protein-shaking-at-body-temperature.json")
        .to_string_lossy()
        .into_owned();
    let d = dump(&[&protein]);
    let rows = drawn(&d);
    let label = rows
        .iter()
        .find(|r| r.text == "title")
        .unwrap_or_else(|| panic!("the scene's form has no title line:\n{d}"));
    let edge = rows
        .iter()
        .find(|r| r.text == "Assets")
        .unwrap_or_else(|| panic!("no asset panel to stay left of:\n{d}"))
        .x;
    let title = "crambin at body temperature, moving along the collective modes its own fold has";
    let value = rows
        .iter()
        .find(|r| r.text == title && r.x > label.x && r.x < edge)
        .unwrap_or_else(|| panic!("the title is not drawn in its field:\n{d}"));
    assert!(
        value.x + value.w <= edge,
        "the title runs from x={} to x={}, past the asset panel at x={edge}:\n{d}",
        value.x,
        value.x + value.w
    );
}
