//! **The viewer refused every scene that ships.**
//!
//! `pantometry view` drew `Panel::Paths` and nothing else, and said so rather than opening a blank
//! window — which was the right thing to say and hid how far it went: `PanelData::paths` is built
//! by two examples and by tests and by **no `Domain` at all**. Not one of the thirty shipped
//! scenes produced a panel this shell could draw, so the only thing it ever said about a scene
//! was that it could not draw it.
//!
//! A body and a field sample are points, and a point is two short segments in screen space, which
//! is the pipeline that was already there. This runs every scene, renders frame 0 offscreen, and
//! asks how much of the canvas came back with something on it.
//!
//! # What this cannot see
//!
//! Whether the picture is *right*. It counts lit pixels: a renderer that drew the correct number
//! of crosses in the wrong places would pass. What it holds is the thing that was wrong — a shell
//! that drew nothing at all for every input it ships with.
//!
//! # Not under `wasm32`, and not without an adapter
//!
//! The scenes are read off a disk and the render wants a GPU. A machine with no adapter skips
//! loudly, which is a result; a skip that says nothing is the shape of a suite that has stopped
//! testing anything.

#![cfg(not(target_family = "wasm"))]

/// The binary, and a place to put what it writes.
fn tools() -> (std::path::PathBuf, std::path::PathBuf) {
    let mut exe = std::env::current_exe().expect("the test binary knows where it is");
    exe.pop();
    if exe.ends_with("deps") {
        exe.pop();
    }
    let bin = exe.join(format!("pantometry{}", std::env::consts::EXE_SUFFIX));
    let dir = std::env::temp_dir().join("pantometry-draws-every-scene");
    std::fs::create_dir_all(&dir).expect("a temp dir");
    (bin, dir)
}

/// **Every shipped scene draws something.**
#[test]
fn every_shipped_scene_puts_something_on_the_canvas() {
    let (bin, out) = tools();
    let scenes = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("pantometry-app has a parent")
        .join("pantometry-world/scenes");
    let mut paths: Vec<_> = std::fs::read_dir(&scenes)
        .expect("the scenes directory reads")
        .filter_map(Result::ok)
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "json"))
        .collect();
    paths.sort();
    assert!(
        paths.len() >= 28,
        "only {} scenes in {}",
        paths.len(),
        scenes.display()
    );

    let mut drawn = 0;
    let mut blank = Vec::new();
    let mut empty = Vec::new();
    let mut skipped = 0;
    for scene in paths {
        let name = scene
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .into_owned();
        let run = out.join(format!("{name}.json"));
        let ran = std::process::Command::new(&bin)
            .args(["run", &scene.to_string_lossy(), &run.to_string_lossy()])
            .output()
            .expect("the binary runs");
        assert!(
            ran.status.success(),
            "{name}: the run refused: {}",
            String::from_utf8_lossy(&ran.stderr)
        );

        let shot = out.join(format!("{name}.ppm"));
        let drew = std::process::Command::new(&bin)
            .args([
                "view",
                &run.to_string_lossy(),
                "--snapshot",
                &shot.to_string_lossy(),
            ])
            .output()
            .expect("the binary runs");
        let said = String::from_utf8_lossy(&drew.stdout).into_owned()
            + &String::from_utf8_lossy(&drew.stderr);
        // A machine with no adapter is not a failing renderer. It says so, and this says which.
        if said.contains("no adapter") || said.contains("no GPU") {
            skipped += 1;
            continue;
        }
        // A scene with nothing to draw is a fact about the scene, not a failing renderer: a
        // network, a lump and a source are readings rather than places. Collected and pinned
        // below rather than treated as either a pass or a failure.
        if said.contains("no panels at all") {
            empty.push(name);
            continue;
        }
        assert!(drew.status.success(), "{name}: the viewer refused: {said}");

        // `--snapshot` reports how much of the canvas is not background, which is the number this
        // test is about. Parsed from its own line rather than by reading the PPM again.
        let lit = said
            .split_once(" — ")
            .and_then(|(_, rest)| rest.split_whitespace().next())
            .and_then(|n| n.parse::<u64>().ok())
            .unwrap_or_else(|| panic!("{name}: no pixel count in: {said}"));
        if lit == 0 {
            blank.push(name);
        } else {
            drawn += 1;
        }
    }

    if skipped > 0 {
        println!("  {skipped} scenes skipped: this machine has no GPU adapter");
    }
    println!("  {drawn} scenes drew something");
    println!("  {} with no panel at all: {empty:?}", empty.len());
    assert!(
        blank.is_empty(),
        "these scenes rendered an empty canvas: {blank:?}"
    );
    // **Pinned, in both directions.** Twenty-seven of the thirty draw; three carry no panel at
    // all, and they are the three whose domains are a `network` and a `winding` — readings, not
    // places, which the editor's viewport says in as many words. A fourth arriving means either
    // a scene lost its geometry or a domain stopped reporting one, and both are worth a failure
    // rather than a quieter list.
    assert_eq!(
        empty,
        [
            "11-motor-thermal-network",
            "12-winding-heats-a-motor",
            "13-winding-that-heats-itself"
        ],
        "the set of scenes with nothing to draw has changed"
    );
    // A run in which every scene skipped would report nothing and pass, which is the shape of a
    // check that turned itself off.
    assert!(
        drawn > 0 || skipped > 0,
        "no scene was either drawn or skipped — the walk broke"
    );
}
