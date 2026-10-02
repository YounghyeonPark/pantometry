//! **Every preset, opened the way a person opens one, checks without an error.**
//!
//! The New-project screen offers the thirty-three shipped scenes, and picking one opens it as an
//! unsaved `scene.json`. Three of them name a file beside themselves — the protein's
//! `structures/1CRN.pdb`, the bracket's `parts/l-bracket.stl` and, since scene 33, aspirin's
//! `structures/AIN.cif` — and an unsaved scene's "beside"
//! is whatever directory the editor was started from. So the protein opened as
//! `crambin: structures/1CRN.pdb (at structures/1CRN.pdb): The system cannot find the path
//! specified` and a viewport saying the scene had no geometry.
//!
//! Nothing had opened a preset through the editor. The test that built every preset — then
//! `every_preset_builds_and_says_when_it_needs_a_file` — built each one *beside the scenes
//! directory*, which is exactly the place an opened preset is not, and it passed.
//!
//! These click each tile and read the frame that results, from a directory of their own so what
//! the editor finds on disk is what the test put there. The first found the protein and the
//! bracket, and only those two, before the files travelled with the presets.
//!
//! # What `unearned-pass-hunter` found in the first version
//!
//! It looked for ` (at ` — the wording of a file *not found* — and passed with each carried file
//! served half-length, because a truncated PDB is a different error; the dump prints the check's
//! own error now, and the control below proves that line appears when there is one. Nothing ran
//! a preset, since a dump's run came before its clicks. And Revert kept a preset's files attached
//! to whatever `scene.json` it loaded. Each has a test here.

use pantometry_world::presets::{Preset, AREAS, PRESETS};
use std::path::{Path, PathBuf};

/// The binary, asked for one frame from `dir`, with an empty recent list.
fn dump_in(dir: &Path, args: &[&str]) -> String {
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
            .current_dir(dir)
            .env("PANTOMETRY_RECENT", none)
            .output()
            .expect("the binary runs");
    assert!(
        out.status.success(),
        "--ui-dump {argv:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// A fresh, empty directory for one test to start the editor in.
fn empty_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("pantometry-preset-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("a temporary directory");
    dir
}

/// Where a dump line's text begins: two marker columns, three five-wide numbers, two spaces.
const TEXT_AT: usize = 22;

/// The `left` and `top` of the line that says exactly `text`.
fn at(dump: &str, text: &str) -> Option<(f32, f32)> {
    let line = dump.lines().find(|l| l.get(TEXT_AT..) == Some(text))?;
    let mut n = line
        .get(..TEXT_AT)?
        .split_whitespace()
        .filter_map(|w| w.parse::<f32>().ok());
    Some((n.next()?, n.next()?))
}

/// The middle of the line that says exactly `text`, which is where a person clicks it.
fn centre(dump: &str, text: &str) -> String {
    let (left, top) = at(dump, text).unwrap_or_else(|| panic!("nothing says `{text}`:\n{dump}"));
    format!("{},{}", left + 8.0, top + 7.0)
}

/// Tall enough that every tile of the largest area is laid out without scrolling.
const HEIGHT: &str = "3000";

/// The arguments that open `p` by clicking its tile, from a frame of its area.
fn opening(dir: &Path, p: &Preset) -> Vec<String> {
    let a = AREAS
        .iter()
        .position(|(key, _, _)| *key == p.area)
        .expect("a preset's area is one of the areas")
        .to_string();
    let area = dump_in(dir, &["--new", "--open", &a, "--height", HEIGHT]);
    let (left, top) = at(&area, p.title)
        .unwrap_or_else(|| panic!("area {a} does not show `{}`:\n{area}", p.title));
    // The tile is a 240 x 156 button above its title; its middle is what a person clicks.
    let tile = format!("{},{}", left + 120.0, top - 80.0);
    ["--new", "--open", &a, "--height", HEIGHT, "--click", &tile]
        .map(String::from)
        .to_vec()
}

fn args(v: &[String]) -> Vec<&str> {
    v.iter().map(String::as_str).collect()
}

/// The check's own error, which the dump prints in its header when there is one.
fn error(frame: &str) -> Option<&str> {
    frame.lines().find_map(|l| l.strip_prefix("error="))
}

/// **The control: a scene whose file is missing does print `error=`.** Without this, "no `error=`
/// line" below would be as true of a dump that had stopped printing the line as of a scene that
/// checked clean.
#[test]
fn a_missing_file_is_reported_as_the_checks_error() {
    let dir = empty_dir("control");
    let protein = PRESETS
        .iter()
        .find(|p| p.files.iter().any(|(n, _)| n.ends_with(".pdb")))
        .expect("a preset carries a PDB");
    std::fs::write(dir.join("scene.json"), protein.json).expect("the scene, without its file");
    let frame = dump_in(&dir, &["scene.json"]);
    let said = error(&frame).unwrap_or_else(|| panic!("no error= line:\n{frame}"));
    assert!(
        said.contains(protein.files[0].0),
        "the error names another thing: {said}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn every_preset_opens_without_an_error() {
    let dir = empty_dir("open");
    let mut opened = 0;
    let mut refused = Vec::new();
    for p in &PRESETS {
        let frame = dump_in(&dir, &args(&opening(&dir, p)));
        // **It opened this preset**, or every assertion below is about the chooser. The titles
        // are unique and none is a prefix of another, so this cannot be another tile's.
        let said = format!("a new scene from {}", p.title);
        assert!(
            frame.contains(&said),
            "{}: clicking its tile did not open it:\n{frame}",
            p.file
        );
        opened += 1;
        if let Some(e) = error(&frame) {
            refused.push(format!("{}: {e}", p.file));
        }
        // And the asset panel says what the scene carries, since the directory it lists
        // cannot: "no .pdb beside ." under a protein that has one reads as the file missing.
        for (name, _) in p.files {
            assert!(
                frame.contains(&format!("{name} — carried until saved")),
                "{}: the asset panel does not say it carries {name}:\n{frame}",
                p.file
            );
        }
    }
    println!("  {opened} presets opened");
    assert_eq!(opened, PRESETS.len(), "not every preset was reached");
    assert!(
        refused.is_empty(),
        "these presets open as an error:\n{}",
        refused.join("\n")
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// **A preset that carries a file runs from it.** Checking and running are separate calls that
/// each take the scene's files, and only the first was reached by a dump until a run could come
/// after a click.
#[test]
fn a_preset_that_carries_a_file_runs() {
    let dir = empty_dir("run");
    let carrying: Vec<&Preset> = PRESETS.iter().filter(|p| !p.files.is_empty()).collect();
    assert!(!carrying.is_empty(), "no preset carries a file");
    for p in carrying {
        let mut a = opening(&dir, p);
        a.push("--run-after-clicks".into());
        let frame = dump_in(&dir, &args(&a));
        assert!(!frame.starts_with("the run refused"), "{}: {frame}", p.file);
        assert!(
            frame.contains("ran: "),
            "{}: nothing says it ran:\n{frame}",
            p.file
        );
        // **And once run, it draws.** Every empty-viewport sentence is absent, which is the
        // protein's forty-six residues and the bracket's field having somewhere to be.
        for checks in [true, false] {
            for ran in [true, false] {
                let empty = editor_core::nothing_to_draw(checks, ran);
                assert!(
                    !frame.contains(empty),
                    "{}: ran, and the viewport says {empty:?}",
                    p.file
                );
            }
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// **An empty viewport says why it is empty.**
///
/// It said one sentence in every case — "nothing in this scene has geometry — sources, lumps and
/// networks are readings, not places" — and the protein opened from New project got it twice
/// wrong: the scene had not checked, and a protein is not a reading. Each of the three cases,
/// shown where it happens:
///
/// - a scene that does not check: the protein, written out without its PDB;
/// - one that checks and has not run: the protein opened from its tile;
/// - one that has run and has nothing to draw: a preset with no tile, which is exactly the set
///   of scenes whose run carries no panel — a network and two windings, reporting readings.
#[test]
fn an_empty_viewport_says_why_it_is_empty() {
    // **Three sentences, not one.** Each check below compares a frame against the function's own
    // answer, so a function answering one sentence for every case — which is the defect — passed
    // all three: measured, with `(true, false) | (true, true)` sharing an arm.
    let said: Vec<&str> = [(false, false), (true, false), (true, true)]
        .into_iter()
        .map(|(checks, ran)| editor_core::nothing_to_draw(checks, ran))
        .collect();
    assert!(
        said[0] != said[1] && said[1] != said[2] && said[0] != said[2],
        "the three cases share a sentence: {said:?}"
    );
    assert!(
        !said[1].contains("readings"),
        "a scene that has not run is told it reports readings: {:?}",
        said[1]
    );

    let dir = empty_dir("empty");
    let protein = PRESETS
        .iter()
        .find(|p| p.files.iter().any(|(n, _)| n.ends_with(".pdb")))
        .expect("a preset carries a PDB");

    std::fs::write(dir.join("broken.json"), protein.json).expect("the scene, without its file");
    let broken = dump_in(&dir, &["broken.json"]);
    let said = editor_core::nothing_to_draw(false, false);
    assert!(
        broken.contains(said),
        "a scene that does not check:\n{broken}"
    );

    let unrun = dump_in(&dir, &args(&opening(&dir, protein)));
    let said = editor_core::nothing_to_draw(true, false);
    assert!(
        unrun.contains(said),
        "the protein, before it runs:\n{unrun}"
    );

    let readings_only = PRESETS
        .iter()
        .find(|p| p.thumb.is_none())
        .expect("a preset with nothing to draw");
    std::fs::write(dir.join("readings.json"), readings_only.json).expect("the scene");
    let readings = dump_in(&dir, &["readings.json", "--ran"]);
    let said = editor_core::nothing_to_draw(true, true);
    assert!(
        readings.contains(said),
        "a scene that reports readings, after its run:\n{readings}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// **Saving writes what the preset carried beside the scene, and lets go of it.** From then on
/// the scene is a file with a directory, and reads its files from there.
#[test]
fn saving_a_preset_writes_its_files_and_lets_go_of_them() {
    let dir = empty_dir("save");
    let p = PRESETS
        .iter()
        .find(|p| !p.files.is_empty())
        .expect("a preset carries a file");
    let mut a = opening(&dir, p);
    let file_menu = centre(&dump_in(&dir, &args(&a)), "File");
    a.extend(["--click".into(), file_menu]);
    let save = centre(&dump_in(&dir, &args(&a)), "Save");
    a.extend(["--click".into(), save]);
    let frame = dump_in(&dir, &args(&a));

    assert_eq!(
        std::fs::read_to_string(dir.join("scene.json")).expect("the scene was saved"),
        p.json
    );
    for (name, bytes) in p.files {
        assert_eq!(
            std::fs::read(dir.join(name)).unwrap_or_else(|e| panic!("{name}: {e}")),
            *bytes,
            "{name} was not written beside the scene"
        );
        assert!(
            !frame.contains("carried until saved"),
            "saved, and still carrying:\n{frame}"
        );
    }
    assert_eq!(error(&frame), None, "saved and reread, it does not check");
    let _ = std::fs::remove_dir_all(&dir);
}

/// **Revert lets go of a preset's files.** After a preset the path is `scene.json` in the
/// working directory, and Revert loads whatever is there — a different scene, which was being
/// served the preset's files ahead of its own directory's and would have had them saved beside it.
#[test]
fn revert_after_a_preset_lets_go_of_its_files() {
    let dir = empty_dir("revert");
    let other = PRESETS
        .iter()
        .find(|p| p.files.is_empty())
        .expect("a preset that carries nothing");
    std::fs::write(dir.join("scene.json"), other.json).expect("somebody's own scene.json");
    let p = PRESETS
        .iter()
        .find(|p| !p.files.is_empty())
        .expect("a preset carries a file");
    let mut a = opening(&dir, p);
    let file_menu = centre(&dump_in(&dir, &args(&a)), "File");
    a.extend(["--click".into(), file_menu]);
    let revert = centre(&dump_in(&dir, &args(&a)), "Revert");
    a.extend(["--click".into(), revert]);
    let frame = dump_in(&dir, &args(&a));

    assert!(
        frame.contains(other.title),
        "Revert did not load the scene.json that was there:\n{frame}"
    );
    assert!(
        !frame.contains("carried until saved"),
        "the reverted scene is still carrying the preset's files:\n{frame}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}
