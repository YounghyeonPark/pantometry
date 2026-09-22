//! **The editor's figures are photographs of a window, and no command in CI can retake them.**
//!
//! The three SVGs on the front page are an example's output — `cargo run --example lens_spots --
//! docs/lens-achromat.svg` — and CI runs those examples on every commit, so a figure that stopped
//! being true would take a failing example with it. The editor's need a display and a GPU, CI has
//! neither, and every change to the interface since they were taken would have left them quietly
//! wrong: the toolbar they show, the seven menus, the bands along the bottom of the viewport.
//!
//! So each picture has a caption that is machine-readable. `tools/screenshot/take.ps1` writes
//! both: the PNG, and the `.txt` beside it — the same frame through `--ui-dump`, which is the same
//! egui layout with no window and no GPU. This regenerates that text and compares it, for **both**
//! figures: the bracket, which is a field and a mesher, and the protein, which is bodies and a
//! different set of readings.
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
fn scene(file: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("pantometry-app has a parent")
        .join("pantometry-world/scenes")
        .join(file)
}

/// **Both editor figures, because there are two of them now.**
///
/// `docs/editor.png` is the bracket, which is about the mesher and a field; `docs/editor-protein.png`
/// is a protein, which is about bodies and a different set of readings. A guard written for one of
/// them by name is a guard the second figure quietly does not have.
#[test]
fn the_screenshots_show_the_editor_as_it_is() {
    for (file, stem) in [
        ("29-a-designed-bracket-becomes-cells.json", "editor"),
        (
            "31-a-protein-shaking-at-body-temperature.json",
            "editor-protein",
        ),
    ] {
        holds(file, stem);
    }
}

/// One figure against a fresh `--ui-dump` of the scene it is of.
fn holds(file: &str, stem: &str) {
    let repo = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("app/pantometry-app has two ancestors")
        .to_path_buf();
    let Ok(stored) = std::fs::read_to_string(repo.join(format!("docs/{stem}.txt"))) else {
        // A checkout without the figure is not this test's business — `docs/` is documentation and
        // a packaging arrangement that dropped it is somebody else's problem.
        return;
    };
    assert!(
        repo.join(format!("docs/{stem}.png")).is_file(),
        "docs/{stem}.txt is here and docs/{stem}.png is not — the caption outlived its picture"
    );

    let mut bin = std::env::current_exe().expect("the test binary knows where it is");
    bin.pop();
    if bin.ends_with("deps") {
        bin.pop();
    }
    let out =
        std::process::Command::new(bin.join(format!("pantometry{}", std::env::consts::EXE_SUFFIX)))
            .args(["--ui-dump", &scene(file).to_string_lossy(), "--ran"])
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
        "the editor's frame has changed since docs/{stem}.png was taken, from line {}:\n\
         \n  docs/{stem}.txt says\n{}\n\n  the editor now draws\n{}\n\n\
         Retake it from app/ with `powershell -File ../tools/screenshot/take.ps1 -Scene \
         pantometry-world/scenes/{file} -Out ..\\docs\\{stem}.png`, which writes both files. \
         Editing the text alone would restore this green and leave the PNG as stale as it is.",
        at + 1,
        show(&stored),
        show(&fresh)
    );
}
