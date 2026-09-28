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

/// A reader's error, with the first lines of what it was reading: the reader reports a line and a
/// column, and the useful thing is to see what is at it.
fn read_error(json: &str, e: &str) -> String {
    let head: String = json.lines().take(9).collect::<Vec<_>>().join("\n");
    format!("{e}\n{head}")
}

/// **Every frame a streaming run emits is a run the viewer can read**, and appended they are the
/// run.
///
/// **This took twenty-six minutes, and the test was not what was slow.** The stream sent the whole
/// run so far after every frame, so the cavity scene's 1600 frames serialised about 31 GB to
/// reach a 38.5 MB run, while every shipped scene *runs* in under two minutes of a debug build.
/// The same cost fell on the editor whenever somebody pressed run. Each payload is one frame now.
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

    let scenes = entries.len();
    for path in entries {
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        let text = std::fs::read_to_string(&path).expect("a scene reads");
        let beside = pantometry_world::Beside::of(&path);
        let stop = std::sync::atomic::AtomicBool::new(false);
        let mut appended: Vec<f64> = Vec::new();
        let mut whole: Option<String> = None;
        let mut unread = Vec::new();
        let end = editor_core::run_streaming(&text, &beside, &stop, |streamed| match streamed {
            editor_core::Streamed::Frame(json) => match viewer_core::Run::from_json(&json) {
                Ok(run) => appended.extend(run.frames.iter().map(|f| f.t)),
                Err(e) => unread.push(read_error(&json, &e.to_string())),
            },
            editor_core::Streamed::Whole(json) => whole = Some(json),
        });
        // **A refused run is a failure here.** It was a `continue` with a line printed, on the
        // argument that a scene the kernel refuses is another test's business — which made this
        // one pass over any scene that stopped running, having checked nothing about it.
        if let Err(why) = end {
            broken.push(format!("{name}: the run refused ({why})"));
            continue;
        }
        broken.extend(unread.into_iter().map(|e| format!("{name}: a frame: {e}")));
        let Some(json) = whole else {
            broken.push(format!("{name}: the stream ended without its whole run"));
            continue;
        };
        match viewer_core::Run::from_json(&json) {
            Err(e) => broken.push(format!(
                "{name}: the whole run: {}",
                read_error(&json, &e.to_string())
            )),
            // The frames appended one at a time are the run's frames, in order: each carries
            // its own instant, so a frame streamed twice or skipped shows up here.
            Ok(run) => {
                let want: Vec<f64> = run.frames.iter().map(|f| f.t).collect();
                if appended != want {
                    broken.push(format!(
                        "{name}: {} frames appended, {} in the run",
                        appended.len(),
                        want.len()
                    ));
                }
            }
        }
        checked += 1;
    }
    println!("  {checked} scenes streamed and read back");
    assert!(broken.is_empty(), "{}", broken.join("\n\n"));
    assert_eq!(checked, scenes, "not every scene was streamed");
}
