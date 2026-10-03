//! **Scene 14 opened on three planets' markers, and its bar could not be looked at.**
//!
//! The scene puts a 20 mm bar, a 4.4 m room and orbits spanning 2.7e11 m in one run. The viewer
//! framed every panel at once, which is the right default and the wrong only choice: the bar is
//! `7.5e-14` of that framing, every one of its 61 samples reached the GPU at one `f32` position,
//! and zoom could not help — `Camera::zoom` clamps the distance to 1.2..9, and a focal length
//! cannot spread a point.
//!
//! `--frame-panel NAME` frames one panel, and `F` cycles through them in the window. This checks
//! the snapshot end to end: the picture of the bar lights far more of the frame than the picture
//! of the whole run, and a name the run does not have is refused rather than ignored.
//!
//! # What the count is of
//!
//! Mostly not the bar. Framed on it, the bar's own crosses are a diagonal line of 61 small `+`
//! marks, and the 4.4 m room the bar lies along the edge of fills the half of the frame beside
//! it — which is the honest picture: the room is not far from the bar, it touches it. The bar
//! itself is checked on the CPU path in `view.rs`, where its 61 samples are 61 positions; this is
//! the check that the flag reaches the renderer and the renderer draws what it was framed on.
//!
//! # Not under `wasm32`, and not without an adapter
//!
//! A machine with no GPU skips loudly, as its neighbours do. The refusal needs no GPU and runs
//! everywhere.

#![cfg(not(target_family = "wasm"))]

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// The `pantometry` binary beside this test.
fn binary() -> PathBuf {
    let mut bin = std::env::current_exe().expect("the test binary knows where it is");
    bin.pop();
    if bin.ends_with("deps") {
        bin.pop();
    }
    bin.join(format!("pantometry{}", std::env::consts::EXE_SUFFIX))
}

/// Scene 14, run into a run file named for the test asking, so two tests running at once do not
/// write one file under each other.
fn scene_14_run(name: &str) -> (PathBuf, PathBuf) {
    let dir = std::env::temp_dir().join("pantometry-one-panel-framed");
    std::fs::create_dir_all(&dir).expect("a temporary directory");
    let scene =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../pantometry-world/scenes/14-a-world.json");
    let run = dir.join(format!("{name}.json"));
    let out = Command::new(binary())
        .args(["run", &scene.to_string_lossy(), &run.to_string_lossy()])
        .output()
        .expect("the binary runs");
    assert!(
        out.status.success(),
        "scene 14 did not run: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    (dir, run)
}

fn said(out: &Output) -> String {
    String::from_utf8_lossy(&out.stdout).into_owned() + &String::from_utf8_lossy(&out.stderr)
}

/// One headless frame as `(width, height, rgb)`, or `None` when there is no adapter.
fn snapshot(run: &Path, shot: &Path, extra: &[&str]) -> Option<(usize, usize, Vec<u8>)> {
    let mut args = vec![
        "view".to_string(),
        run.to_string_lossy().into_owned(),
        "--snapshot".to_string(),
        shot.to_string_lossy().into_owned(),
    ];
    args.extend(extra.iter().map(|s| s.to_string()));
    let out = Command::new(binary())
        .args(&args)
        .output()
        .expect("the binary runs");
    let text = said(&out);
    // Exit 3 is the viewer's own "no GPU available for a snapshot".
    if out.status.code() == Some(3) || text.contains("no adapter") {
        return None;
    }
    assert!(out.status.success(), "the snapshot refused:\n{text}");
    Some(ppm(&std::fs::read(shot).expect("the snapshot is written")))
}

/// A binary PPM, read token by token: `P6`, width, height, maximum, one whitespace byte, pixels.
fn ppm(bytes: &[u8]) -> (usize, usize, Vec<u8>) {
    let mut at = 0;
    let mut token = || {
        while bytes[at].is_ascii_whitespace() {
            at += 1;
        }
        let start = at;
        while !bytes[at].is_ascii_whitespace() {
            at += 1;
        }
        String::from_utf8_lossy(&bytes[start..at]).into_owned()
    };
    assert_eq!(token(), "P6", "not a binary PPM");
    let w: usize = token().parse().expect("a width");
    let h: usize = token().parse().expect("a height");
    assert_eq!(token(), "255", "not eight bits a channel");
    let pixels = bytes[at + 1..].to_vec();
    assert_eq!(pixels.len(), w * h * 3, "the PPM is not the size it says");
    (w, h, pixels)
}

/// Pixels that are not `background`.
fn lit(rgb: &[u8], background: [u8; 3]) -> usize {
    rgb.as_chunks::<3>()
        .0
        .iter()
        .filter(|p| **p != background)
        .count()
}

#[test]
fn a_panel_the_run_does_not_have_is_refused() {
    let (dir, run) = scene_14_run("refused");
    let out = Command::new(binary())
        .args([
            "view",
            &run.to_string_lossy(),
            "--frame-panel",
            "nothing-by-that-name",
            "--snapshot",
            &dir.join("refused.ppm").to_string_lossy(),
        ])
        .output()
        .expect("the binary runs");
    let text = said(&out);
    // **Refused, and with the names it would have taken.** A flag that ignored a name it did not
    // know would write the whole run's picture and exit 0, which is this test's other half
    // passing for the wrong reason.
    assert_eq!(
        out.status.code(),
        Some(2),
        "a wrong name was not refused:\n{text}"
    );
    assert!(
        text.contains("bar, room, sky"),
        "the refusal did not say which panels there are:\n{text}"
    );
}

#[test]
fn framing_the_bar_lights_far_more_of_the_frame_than_the_whole_run() {
    let (dir, run) = scene_14_run("lit");
    let Some((w, h, whole)) = snapshot(&run, &dir.join("whole.ppm"), &[]) else {
        println!(
            "no GPU adapter on this machine — the framed snapshot was not rendered or measured"
        );
        return;
    };
    let (fw, fh, framed) = snapshot(&run, &dir.join("bar.ppm"), &["--frame-panel", "bar"])
        .expect("the adapter is still here");
    assert_eq!((w, h), (fw, fh), "the two snapshots are different sizes");

    // **The clear colour, from the picture of the whole run.** Its corner is empty because the fit
    // puts the subject inside 0.85 of the half-frame and the legend starts inside the margin; the
    // framed picture is read against the same colour rather than against its own corner, which a
    // framed subject can cover.
    let background = [whole[0], whole[1], whole[2]];
    let (a, b) = (lit(&whole, background), lit(&framed, background));
    let all = w * h;
    println!(
        "  whole run: {a} of {all} pixels lit ({:.2}%); framed on the bar: {b} ({:.2}%)",
        100.0 * a as f64 / all as f64,
        100.0 * b as f64 / all as f64
    );
    // Measured at 13 318 against 400 616, a factor of 30.1. Ten is the claim: framing one panel
    // fills the frame with it and its neighbourhood, where the whole run is a few markers and a
    // legend. An unrecognised flag reproduces the whole run's picture exactly, a factor of one.
    assert!(
        b >= 10 * a,
        "framed on the bar the snapshot lit {b} pixels against the whole run's {a}; the flag did \
         not reach the renderer, or the renderer did not draw what it was framed on"
    );
}
