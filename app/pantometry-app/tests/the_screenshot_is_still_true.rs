//! **`docs/editor.png` is the one figure no command can refresh, so nothing noticed it ageing.**
//!
//! The other two are an example's output — `cargo run --example lens_spots -- docs/lens-achromat.svg`
//! — and CI runs those examples on every commit, so a figure that stopped being true would take a
//! failing example with it. The editor's is a photograph of a window. It needs a display and a GPU,
//! CI has neither, and every change to the interface since it was taken would have left it quietly
//! wrong: the toolbar it shows, the seven menus, the bands along the bottom of the viewport.
//!
//! So the picture has a caption that is machine-readable. `tools/screenshot/take.ps1` writes both:
//! the PNG, and `docs/editor.txt` — the same frame through `--ui-dump`, which is the same egui
//! layout with no window and no GPU. This regenerates that text and compares it.
//!
//! # What this can and cannot say
//!
//! It says **the frame the picture shows is still the frame the editor draws**: the same strings in
//! the same places, the same viewport rect, the same counts. It cannot say the *pixels* are right —
//! a change to a colour, a font or the shaded pass moves nothing in the dump. What it holds is the
//! failure that was actually coming: a screenshot of an interface that has since been rearranged.
//!
//! When it fails, the fix is to retake the picture, not to update the text: the text is written by
//! the same script that takes the picture, so editing it alone would restore the green and leave
//! the PNG exactly as stale as it was.
//!
//! # Not under `wasm32`
//!
//! It reads a document off a disk and runs a scene.

#![cfg(not(target_family = "wasm"))]

/// The scene the picture is of, from this crate's manifest.
fn scene() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("pantometry-app has a parent")
        .join("pantometry-world/scenes/29-a-designed-bracket-becomes-cells.json")
}

/// **The frame the screenshot shows is the frame the editor draws.**
#[test]
fn the_screenshot_shows_the_editor_as_it_is() {
    let repo = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("app/pantometry-app has two ancestors")
        .to_path_buf();
    let Ok(stored) = std::fs::read_to_string(repo.join("docs/editor.txt")) else {
        // A checkout without the figure is not this test's business — `docs/` is documentation and
        // a packaging arrangement that dropped it is somebody else's problem.
        return;
    };
    assert!(
        repo.join("docs/editor.png").is_file(),
        "docs/editor.txt is here and docs/editor.png is not — the caption outlived its picture"
    );

    let mut bin = std::env::current_exe().expect("the test binary knows where it is");
    bin.pop();
    if bin.ends_with("deps") {
        bin.pop();
    }
    let out =
        std::process::Command::new(bin.join(format!("pantometry{}", std::env::consts::EXE_SUFFIX)))
            .args(["--ui-dump", &scene().to_string_lossy(), "--ran"])
            .output()
            .expect("the binary runs");
    assert!(
        out.status.success(),
        "--ui-dump refused: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let fresh = String::from_utf8_lossy(&out.stdout).into_owned();

    // Carriage returns are the checkout's, not the frame's: `core.autocrlf` rewrites the stored
    // file and the binary's own output has none.
    let flat = |s: &str| s.replace('\r', "");
    let (stored, fresh) = (flat(&stored), flat(&fresh));
    if stored == fresh {
        return;
    }

    // The first line that differs, because eighty-six lines of diff is not what a reader needs.
    let at = stored
        .lines()
        .zip(fresh.lines())
        .position(|(a, b)| a != b)
        .unwrap_or_else(|| stored.lines().count().min(fresh.lines().count()));
    let show = |s: &str| {
        s.lines()
            .skip(at.saturating_sub(2))
            .take(5)
            .collect::<Vec<_>>()
            .join("\n")
    };
    panic!(
        "the editor's frame has changed since `docs/editor.png` was taken, from line {}:\n\
         \n  docs/editor.txt says\n{}\n\n  the editor now draws\n{}\n\n\
         Retake the picture — `cd app && powershell -File ../tools/screenshot/take.ps1` — which \
         writes both files. Editing `editor.txt` alone would restore this green and leave the PNG \
         as stale as it is.",
        at + 1,
        show(&stored),
        show(&fresh)
    );
}
