//! **Pointing at the model, through the window.**
//!
//! `editor_core::pick` is held to closed forms in its own tests, and `viewer_core::Camera::ray` to
//! the projection. None of that matters if a click in the viewport does not reach them, writes a
//! face other than the one under the pointer, or takes the selection away — so this clicks and
//! drags through `--ui-dump`, the way a person does, and reads back what the file then says.
//!
//! The scene is 25: a module of 8 x 8 x 8 cells of 1.5 mm with no void, so every face is 64
//! cells of 2.25 mm², 1.44 cm², and a patch of `n` cells is `n x 2.25 / 100` — closed forms that
//! do not come from the editor. Its structure is the same grid and holds three faces.
//!
//! What a pixel is over is read off the frame itself, from the face names the view writes on the
//! faces toward the eye and the readout it draws under a hovered pointer — the things a person
//! reads to decide where to click.

#![cfg(not(target_family = "wasm"))]

const SCENE: &str = "25-what-140-kelvin-does-to-the-solder.json";

fn scene() -> String {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../pantometry-world/scenes")
        .join(SCENE)
        .to_string_lossy()
        .into_owned()
}

fn dump(args: &[&str]) -> String {
    let mut argv = vec!["--ui-dump".to_string(), scene()];
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

/// One drawn string: its left, top and width, and what it says.
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

/// The viewport's rect, `(x, y, w, h)`.
fn viewport(d: &str) -> (f32, f32, f32, f32) {
    let v = d
        .lines()
        .find_map(|l| l.strip_prefix("viewport="))
        .unwrap_or_else(|| panic!("no viewport:\n{d}"));
    let (at, size) = v.split_once(' ').expect("x,y wxh");
    let (x, y) = at.split_once(',').expect("x,y");
    let (w, h) = size.split_once('x').expect("wxh");
    (
        x.parse().expect("x"),
        y.parse().expect("y"),
        w.parse().expect("w"),
        h.parse().expect("h"),
    )
}

const FACES: [&str; 6] = ["x-min", "x-max", "y-min", "y-max", "z-min", "z-max"];

/// The face names written in the view, with the point at their middle.
fn face_names(d: &str) -> Vec<(String, (f32, f32))> {
    let (vx, vy, vw, vh) = viewport(d);
    drawn(d)
        .into_iter()
        .filter(|r| FACES.contains(&r.text.as_str()))
        .filter(|r| r.x >= vx && r.x + r.w <= vx + vw && r.y >= vy && r.y <= vy + vh)
        .map(|r| (r.text, (r.x + r.w / 2.0, r.y + 7.0)))
        .collect()
}

fn at(p: (f32, f32)) -> String {
    format!("{:.0},{:.0}", p.0, p.1)
}

/// The scene's text as the Source tab shows it, parsed. Without its carriage returns: a checkout
/// may hold the scene with CRLF endings, and the frame reports the text as it is.
fn source(d: &str) -> serde_json::Value {
    let text = drawn(d)
        .into_iter()
        .map(|r| r.text.replace('\r', ""))
        .find(|t| t.starts_with("{\\n") && t.contains("\"domains\""))
        .unwrap_or_else(|| panic!("no source in the frame:\n{d}"));
    serde_json::from_str(&text.replace("\\n", "\n")).unwrap_or_else(|e| panic!("{e}: {text}"))
}

fn entry<'a>(v: &'a serde_json::Value, domain: &str, key: &str, j: usize) -> &'a serde_json::Value {
    &v["domains"]
        .as_array()
        .expect("domains")
        .iter()
        .find(|d| d["name"] == domain)
        .unwrap_or_else(|| panic!("no domain {domain}"))[key][j]
}

/// The status line's sentence: the string starting with `start`.
fn said(d: &str, start: &str) -> String {
    drawn(d)
        .into_iter()
        .find(|r| r.text.starts_with(start))
        .unwrap_or_else(|| panic!("nothing says {start:?}:\n{d}"))
        .text
}

/// What the readout under a hovered pointer says: its two lines.
fn readout(d: &str) -> (String, String) {
    let line = drawn(d)
        .into_iter()
        .find(|r| r.text.contains("\\nclick:") || r.text.contains("\\nshift-drag:"))
        .unwrap_or_else(|| panic!("no readout under the pointer:\n{d}"))
        .text;
    let (what, hint) = line.split_once("\\n").expect("two lines");
    (what.to_string(), hint.to_string())
}

/// `F · cell u, v` read back as the face and the two numbers, with a `· void` after them dropped.
fn cell_of(what: &str) -> (String, [usize; 2]) {
    let what = what.trim_end_matches(" · void");
    let (face, cell) = what
        .split_once(" · cell ")
        .unwrap_or_else(|| panic!("{what}"));
    let n: Vec<usize> = cell
        .split(", ")
        .map(|v| v.parse().unwrap_or_else(|_| panic!("{what}")))
        .collect();
    (face.to_string(), [n[0], n[1]])
}

/// **The faces toward the eye are named on them, and a click on a name writes that face.** Each
/// name is clicked in turn with the module's cooled face selected; the file then says that face
/// and the area 1.44 cm², and the scene still builds. Clicking the face it already has says so and
/// writes nothing.
#[test]
fn a_click_on_a_named_face_writes_that_face() {
    let picking = ["--select", "/domains/module/cooling/0"];
    let d = dump(&picking);
    let names = face_names(&d);
    assert!(
        (1..=3).contains(&names.len()),
        "a box shows one to three faces to an eye outside it, and this view names {}:\n{d}",
        names.len()
    );
    assert!(
        d.contains("overlap=0\n"),
        "the names sit on something:\n{d}"
    );
    // The face the entry has is among them, so the branch below that clicks it runs.
    assert!(
        names.iter().any(|(f, _)| f == "z-min"),
        "the cooled face is not one the eye sees: {:?}",
        names.iter().map(|(f, _)| f).collect::<Vec<_>>()
    );
    for (face, middle) in &names {
        let mut args = picking.to_vec();
        let click = at(*middle);
        args.extend(["--click", &click, "--source"]);
        let d = dump(&args);
        assert!(!d.contains("\nerror="), "the scene stopped building:\n{d}");
        let e = entry(&source(&d), "module", "cooling", 0).clone();
        assert_eq!(e["face"], face.as_str(), "clicked {face} at {click}");
        assert_eq!(e["area_cm2"].as_f64(), Some(1.44), "{e}");
        if face == "z-min" {
            said(&d, "module cooling 1 already says that");
        } else {
            let s = said(&d, "module cooling 1:");
            assert!(s.contains("over 64 solid cells"), "{s}");
        }
    }
    // **An edited frame is as legible as the one before it.** The toolbar's `edited` was laid out
    // over the file's name on every frame after an edit, which a frame with nothing edited cannot
    // show. Read without the Source tab, whose text runs past its pane in the dump's count.
    let (face, middle) = names
        .iter()
        .find(|(f, _)| f != "z-min")
        .expect("a face other than the cooled one is named");
    let mut args = picking.to_vec();
    let click = at(*middle);
    args.extend(["--click", &click]);
    let d = dump(&args);
    said(&d, &format!("module cooling 1: {face}"));
    assert!(
        d.contains("overlap=0\n"),
        "the edited frame has strings over each other:\n{d}"
    );
    // And with nothing selected, the faces are not named: the names are for pointing.
    assert!(face_names(&dump(&[])).is_empty());
}

/// **A shift-drag draws the patch under its two ends.** The readout under the pointer is read at
/// each end first; the drag between them then writes `from` and `to` spanning those two cells, and
/// an area of the cells it covers times 2.25 mm².
#[test]
fn a_shift_drag_draws_the_patch_under_its_two_ends() {
    let picking = ["--select", "/domains/module/cooling/0"];
    let names = face_names(&dump(&picking));
    let (face, (cx, cy)) = names.first().expect("a face is named").clone();
    let (a, b) = ((cx - 20.0, cy - 20.0), (cx + 20.0, cy + 20.0));
    let mut ends = Vec::new();
    for p in [a, b] {
        let mut args = picking.to_vec();
        let hover = at(p);
        args.extend(["--hover", &hover]);
        let (what, hint) = readout(&dump(&args));
        assert!(hint.contains("shift-drag: a patch"), "{hint}");
        let (on, cell) = cell_of(&what);
        assert_eq!(on, face, "{hover} is not on {face}");
        ends.push(cell);
    }
    let from = [ends[0][0].min(ends[1][0]), ends[0][1].min(ends[1][1])];
    let to = [
        ends[0][0].max(ends[1][0]) + 1,
        ends[0][1].max(ends[1][1]) + 1,
    ];
    let cells = (to[0] - from[0]) * (to[1] - from[1]);
    assert!(
        cells > 1 && cells < 64,
        "the two ends span {cells} cells, which tests nothing"
    );

    let mut args = picking.to_vec();
    let drag = format!("{},{}", at(a), at(b));
    args.extend(["--shift-drag", &drag, "--source"]);
    let d = dump(&args);
    let e = entry(&source(&d), "module", "cooling", 0).clone();
    assert_eq!(e["face"], face.as_str());
    assert_eq!(e["from"], serde_json::json!(from), "{e}");
    assert_eq!(e["to"], serde_json::json!(to), "{e}");
    let area = cells as f64 * 1.5 * 1.5 / 100.0;
    assert!(
        (e["area_cm2"].as_f64().expect("a number") - area).abs() < 1e-12,
        "{e} against {cells} cells, {area} cm²"
    );
    assert!(!d.contains("\nerror="), "the patch does not build:\n{d}");

    // **Without Shift the same drag turns the view** and writes nothing: a person turning the
    // model over a selected face is not drawing on it.
    let mut args = picking.to_vec();
    args.extend(["--drag", &drag, "--source"]);
    let e = entry(&source(&dump(&args)), "module", "cooling", 0).clone();
    assert_eq!(e["face"], "z-min", "a plain drag moved the face: {e}");
    assert!(e.get("from").is_none(), "a plain drag drew a patch: {e}");
}

/// **A held face is set by a click, and a heat source's box by a drag that keeps its depth.** The
/// structure's first held face, clicked on a name, is that face. The module's heat source is cells
/// `[0, 0, 7]..[8, 8, 8]`: drawn across a face, its two axes are the cells under the drag's ends,
/// and the third is what it was.
#[test]
fn a_held_face_and_a_heat_source_are_pointed_at() {
    let held = ["--select", "/domains/stress/held/0"];
    let names = face_names(&dump(&held));
    // Not the face it holds already, or a click that did nothing would pass.
    let (face, middle) = names
        .iter()
        .find(|(f, _)| f != "x-min")
        .expect("a face other than the held one is named")
        .clone();
    let mut args = held.to_vec();
    let click = at(middle);
    args.extend(["--click", &click, "--source"]);
    let e = entry(&source(&dump(&args)), "stress", "held", 0).clone();
    assert_eq!(e["face"], face.as_str());
    assert_eq!(e["as"], "roller", "the hold changed with the face: {e}");

    let source_box = ["--select", "/domains/module/dissipation/0"];
    let d = dump(&source_box);
    let names = face_names(&d);
    let (face, (cx, cy)) = names.first().expect("a face is named").clone();
    let (a, b) = ((cx - 20.0, cy - 15.0), (cx + 20.0, cy + 15.0));
    let mut ends = Vec::new();
    for p in [a, b] {
        let mut args = source_box.to_vec();
        let hover = at(p);
        args.extend(["--hover", &hover]);
        let (what, hint) = readout(&dump(&args));
        assert!(hint.contains("footprint"), "{hint}");
        // `F · cell i, j, k` for a box: all three, of which the face's two are what a drag sets.
        let (on, rest) = what.split_once(" · cell ").expect("a cell");
        assert_eq!(on, face);
        let c: Vec<usize> = rest
            .split(", ")
            .map(|v| v.parse().expect("a number"))
            .collect();
        ends.push(c);
    }
    let normal = ["x", "y", "z"]
        .iter()
        .position(|a| face.starts_with(a))
        .expect("an axis");
    let mut args = source_box.to_vec();
    let drag = format!("{},{}", at(a), at(b));
    args.extend(["--shift-drag", &drag, "--source"]);
    let e = entry(&source(&dump(&args)), "module", "dissipation", 0).clone();
    let (was_from, was_to) = ([0, 0, 7], [8, 8, 8]);
    // A drag that spans the whole face would leave the box as it was and pass doing nothing.
    assert!(
        (0..3)
            .filter(|k| *k != normal)
            .any(|k| ends[0][k].min(ends[1][k]) != 0 || ends[0][k].max(ends[1][k]) + 1 != 8),
        "the drag's ends {ends:?} span the whole face"
    );
    for k in 0..3 {
        let (f, t) = if k == normal {
            (was_from[k], was_to[k])
        } else {
            (ends[0][k].min(ends[1][k]), ends[0][k].max(ends[1][k]) + 1)
        };
        assert_eq!(e["from"][k], f, "from[{k}] of {e}");
        assert_eq!(e["to"][k], t, "to[{k}] of {e}");
    }
}

/// **A click beside the face keeps the entry and says it missed.** With a run on screen a click on
/// nothing clears the selection; with an entry being pointed, a click that misses its domain was
/// meant for it.
///
/// Three frames make that a claim about this point and not about any click: hovered, nothing is
/// under it, so it is in the view and off the module; with a domain row selected, which aims at
/// nothing, the same click clears the selection, so it reached the view; and with the cooled face
/// selected it keeps the entry, says it missed, and the file still says `z-min` with no patch.
#[test]
fn a_click_beside_the_face_keeps_the_entry() {
    let pointing = ["--ran", "--select", "/domains/module/cooling/0"];
    let d = dump(&pointing);
    let (vx, vy, _, _) = viewport(&d);
    let beside = at((vx + 12.0, vy + 40.0));

    let mut args = pointing.to_vec();
    args.extend(["--hover", &beside]);
    let hovered = dump(&args);
    assert!(
        !hovered.contains("\\nclick:") && !hovered.contains("\\nshift-drag:"),
        "{beside} is over the module, so a click there is not beside it:\n{hovered}"
    );

    let control = dump(&["--ran", "--select", "/domains/module", "--click", &beside]);
    said(&control, "the scene itself");

    let mut args = pointing.to_vec();
    args.extend(["--click", &beside, "--source"]);
    let after = dump(&args);
    assert!(
        !face_names(&after).is_empty(),
        "a click at {beside}, beside the module, took the cooled face's selection away:\n{after}"
    );
    said(&after, "that click missed module");
    let e = entry(&source(&after), "module", "cooling", 0).clone();
    assert_eq!(e["face"], "z-min", "{e}");
    assert!(e.get("from").is_none(), "{e}");
}

/// **A turned domain is pointed at in its own frame, by both paths.** Scene 30's busbars are moved
/// and turned thirty degrees. The results panel reads them through the pose the run file carries,
/// and a selected cooled face through the pose the scene states: two copies of one placement, by
/// two routes. Under one pixel they must name one place — the face the probe's readout puts at a
/// wall, and the cell its other two coordinates fall in, at 4 mm a cell.
#[test]
fn a_turned_domain_is_one_place_to_both_paths() {
    let posed = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../pantometry-world/scenes/30-two-phases-crossing-at-a-clearance.json")
        .to_string_lossy()
        .into_owned();
    let run = |extra: &[&str]| {
        let mut args = vec![posed.as_str(), "--ran"];
        args.extend(extra);
        // `dump` puts scene 25 first; this one names its own.
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
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        String::from_utf8_lossy(&out.stdout).into_owned()
    };
    let d = run(&["--results", "point"]);
    let (vx, vy, vw, vh) = viewport(&d);
    let middle = at((vx + vw / 2.0, vy + vh / 2.0));

    let probe = run(&["--results", "point", "--hover", &middle]);
    let (what, _) = readout(&probe);
    let p = triple(&what);
    let entry = run(&["--select", "/domains/busbars/cooling/0", "--hover", &middle]);
    assert_eq!(
        viewport(&entry),
        viewport(&probe),
        "two layouts, two pixels"
    );
    let (face, cell) = cell_of(&readout(&entry).0);
    // The busbars are 8 x 8 x 5 cells of 4 mm: 32 x 32 x 20 mm.
    let size = [32.0, 32.0, 20.0];
    let axis = ["x", "y", "z"]
        .iter()
        .position(|a| face.starts_with(a))
        .expect("an axis");
    let wall = if face.ends_with("max") {
        size[axis]
    } else {
        0.0
    };
    assert_eq!(p[axis], wall, "the probe's {what} is not on {face}");
    let others: Vec<usize> = (0..3).filter(|k| *k != axis).collect();
    for (n, k) in others.iter().enumerate() {
        let want = ((p[*k] / 4.0).floor() as usize).min(size[*k] as usize / 4 - 1);
        assert_eq!(
            cell[n], want,
            "axis {k}: the probe at {what}, the face's cell {cell:?}"
        );
    }
}

/// The numbers inside `(a, b, c) mm`.
fn triple(s: &str) -> [f64; 3] {
    let inner = s
        .split_once('(')
        .and_then(|(_, r)| r.split_once(')'))
        .unwrap_or_else(|| panic!("no (x, y, z) in {s}"))
        .0;
    let v: Vec<f64> = inner
        .split(", ")
        .map(|x| x.parse().expect("a number"))
        .collect();
    [v[0], v[1], v[2]]
}

/// **The results panel is pointed at the module.** On the probe tab the readout under the
/// pointer is a point on the module's surface — one coordinate at a wall, 0 or 12 mm — and a click
/// there puts the probe at those coordinates. A drag on the line tab runs from one readout to the
/// other; on the region tab it draws a box. On the slice tab a click takes the plane through the
/// point.
#[test]
fn the_results_are_pointed_at_where_the_pointer_says() {
    let d = dump(&["--ran", "--results", "point"]);
    let (vx, vy, vw, vh) = viewport(&d);
    let middle = (vx + vw / 2.0, vy + vh / 2.0);
    let near = (middle.0 + 25.0, middle.1 + 20.0);

    let mut reads = Vec::new();
    for p in [middle, near] {
        let hover = at(p);
        let (what, hint) = readout(&dump(&["--ran", "--results", "point", "--hover", &hover]));
        assert_eq!(hint, "click: the probe here");
        let q = triple(&what);
        assert!(
            q.iter().any(|v| *v == 0.0 || *v == 12.0) && q.iter().all(|v| (0.0..=12.0).contains(v)),
            "{what} is not on the module's surface"
        );
        reads.push(what);
    }
    let click = at(middle);
    let probe = said(
        &dump(&["--ran", "--results", "point", "--click", &click]),
        "probe at ",
    );
    assert_eq!(probe, format!("probe at {}", reads[0]));

    let drag = format!("{},{}", at(middle), at(near));
    let line = said(
        &dump(&["--ran", "--results", "line", "--shift-drag", &drag]),
        "line from ",
    );
    assert_eq!(line, format!("line from {} to {}", reads[0], reads[1]));
    // Two clicks are the line's two ends, start and then end.
    let second = at(near);
    let two = dump(&[
        "--ran",
        "--results",
        "line",
        "--click",
        &click,
        "--click",
        &second,
    ]);
    assert_eq!(
        said(&two, "line end at "),
        format!("line end at {}", reads[1])
    );

    let region = said(
        &dump(&["--ran", "--results", "region", "--shift-drag", &drag]),
        "region from ",
    );
    let (lo, hi) = region
        .strip_prefix("region from ")
        .and_then(|r| r.split_once(" to "))
        .expect("from and to");
    let (lo, hi, a, b) = (triple(lo), triple(hi), triple(&reads[0]), triple(&reads[1]));
    // Both ends on one face, so the box is drawn on it: the axis across that face is the one both
    // readouts put at a wall.
    let across = (0..3)
        .find(|&k| a[k] == b[k] && (a[k] == 0.0 || a[k] == 12.0))
        .unwrap_or_else(|| panic!("{a:?} and {b:?} are not on one face, which this needs"));
    for k in 0..3 {
        let want = if k == across {
            // The box keeps the depth it had, the whole domain's.
            (0.0, 12.0)
        } else {
            (a[k].min(b[k]), a[k].max(b[k]))
        };
        assert_eq!((lo[k], hi[k]), want, "axis {k} of {region}");
    }

    let plane = said(
        &dump(&["--ran", "--results", "slice", "--click", &click]),
        "plane at ",
    );
    let z = triple(&reads[0])[2];
    assert!(
        plane.starts_with(&format!("plane at z = {z} "))
            || plane.starts_with(&format!("plane at z = {} ", z as i64)),
        "{plane} against the readout's z of {z}"
    );
}

/// **The form says how to point, and what is wrong with what is set.** A selected cooled face's
/// form tells a person they can click a face or draw a patch, a heat source's that a drag draws its
/// box. And scene 29's bracket, whose part fills ten of its grid's twelve layers, states a film on
/// its top layer, which no solid cell is on: the form says it cools nothing, before any click.
#[test]
fn the_form_says_how_to_point_and_what_is_wrong() {
    let d = dump(&["--select", "/domains/module/cooling/0"]);
    said(&d, "in the view: click a face of module to cool all of it");
    assert!(
        !d.contains("cools nothing"),
        "scene 25's film is on solid cells:\n{d}"
    );
    let d = dump(&["--select", "/domains/module/dissipation/0"]);
    said(&d, "in the view: shift-drag across a face of module");

    let bracket = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../pantometry-world/scenes/29-a-designed-bracket-becomes-cells.json")
        .to_string_lossy()
        .into_owned();
    let none = std::env::temp_dir().join("pantometry-no-such-recent-list.json");
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_pantometry"))
        .args([
            "--ui-dump",
            &bracket,
            "--select",
            "/domains/bracket/cooling/1",
        ])
        .env("PANTOMETRY_RECENT", none.to_string_lossy().into_owned())
        .output()
        .expect("the binary runs");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let d = String::from_utf8_lossy(&out.stdout);
    said(
        &d,
        "bracket cooling 2; no solid cell of bracket is on z-max",
    );
}
