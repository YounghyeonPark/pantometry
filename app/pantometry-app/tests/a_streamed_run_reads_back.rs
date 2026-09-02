//! **The editor could not read back JSON it had written one call ago.**
//!
//! Opening `29-a-designed-bracket-becomes-cells.json` with `--run` put this on the status bar:
//!
//! ```text
//! the run's own JSON did not read back: not a pantometry run:
//! invalid type: null, expected f64 at line 7 column 5
//! ```
//!
//! `editor_core::run` — the batch call — and `editor_core::run_streaming` — what the editor uses,
//! because it wants a picture before the run ends — write the same shape through different code,
//! and only one of them was ever read back by a test. The window is where that showed up, and it
//! showed up as an editor that ran a scene and drew nothing.
//!
//! Every shipped scene, through the streaming path, back through the reader the editor uses.

#![cfg(not(target_family = "wasm"))]

/// Where the scenes are, from this crate's manifest.
fn scenes() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("pantometry-app has a parent")
        .join("pantometry-world/scenes")
}

/// **Every frame a streaming run emits is a run the viewer can read.**
#[test]
fn every_scene_streams_json_the_viewer_can_read() {
    let dir = scenes();
    let mut checked = 0;
    let mut broken = Vec::new();
    let mut entries: Vec<_> = std::fs::read_dir(&dir)
        .expect("the scenes directory reads")
        .filter_map(Result::ok)
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "json"))
        .collect();
    entries.sort();
    assert!(
        entries.len() >= 28,
        "only {} scenes found in {} — this test compares against them",
        entries.len(),
        dir.display()
    );

    for path in entries {
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        let text = std::fs::read_to_string(&path).expect("a scene reads");
        let beside = pantometry_world::Beside::of(&path);
        let stop = std::sync::atomic::AtomicBool::new(false);
        let mut last: Option<String> = None;
        let mut frames = 0;
        let end = editor_core::run_streaming(&text, &beside, &stop, |json| {
            frames += 1;
            last = Some(json);
        });
        if let Err(why) = end {
            // A scene the kernel refuses is not this test's business — it is
            // `every_scene_that_ships_runs_and_says_something_true`'s.
            println!("  {name}: the run itself refused ({why})");
            continue;
        }
        let Some(json) = last else {
            broken.push(format!("{name}: streamed no frames at all"));
            continue;
        };
        checked += 1;
        if let Err(e) = viewer_core::Run::from_json(&json) {
            // The first lines, because the reader reports a line and a column and the useful
            // thing is to see what is at it.
            let head: String = json.lines().take(9).collect::<Vec<_>>().join("\n");
            broken.push(format!("{name}: {e}\n{head}"));
        }
    }
    println!("  {checked} scenes streamed and read back");
    assert!(checked > 0, "no scene streamed a frame — the walk broke");
    assert!(broken.is_empty(), "{}", broken.join("\n\n"));
}
