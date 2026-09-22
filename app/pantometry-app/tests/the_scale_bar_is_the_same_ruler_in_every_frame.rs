//! **A scale bar is a measurement, and it was a different one in every frame.**
//!
//! The camera is fitted once, over the union of every frame's bounds, so metres-per-pixel is
//! fixed for a whole run. The bar divided by *this frame's* box instead — so on a run whose
//! contents move it came out the wrong length by the ratio between the two, with the number
//! beside it unchanged. Measured on the committed protein animation: a bar labelled `2 NM` ran
//! 108 px at the ends of the swing and 119 in the middle. 119/108 is 1.102; the frames' own boxes
//! are 5.1351e-9 and 4.6684e-9, a ratio of 1.100.
//!
//! A ruler that changes length is worse than no ruler, and nothing in the picture says so: it
//! looks exactly like a bar either way.
//!
//! # Why this is not the unit test beside it
//!
//! `the_bar_covers_the_metres_it_names` checks the arithmetic, and `scale_bar` now takes the
//! run's longest side and nothing else, so it *cannot* be handed a frame's box. Both of those are
//! about the function. This is about the wiring — that what reaches it is the run's box — and the
//! only way to see the wiring is to render two frames and measure the bar in each.
//!
//! The test that was here before this one asked whether the bar's *label* named a positive
//! length. It did, in all forty-eight frames.
//!
//! # Not under `wasm32`, and not without an adapter
//!
//! It writes a file and renders it twice. A machine with no GPU skips loudly, as its neighbours
//! do: a skip that says why is a result, and one that says nothing is a suite that has stopped
//! testing anything.

#![cfg(not(target_family = "wasm"))]

/// Two frames whose panel declares a **different box in each**, which is what a run of anything
/// moving looks like — `PanelData::surface` computes a panel's bounds from its own positions, so
/// a solid that swings has a box that swings with it.
///
/// The boxes are 1.0 and 1.49 along their longest side. Both fall in the same rung of the bar's
/// ladder, so the *label* is `0.5 M` in both frames and only the bar's length could differ: with
/// the frame's box it is `0.5/1.0 * 0.5 = 0.250` of the half-width against `0.5/1.49 * 0.5 =
/// 0.168`, a bar half again as long in the first frame for the same stated metres. With the run's
/// box, which is the union, both are 0.168.
///
/// **Two runs with different values**, because a run whose colour range was flat used to draw no
/// legend at all — the viewer returned early on `hi <= lo` and the scale bar went with the colour
/// bar. A first version of this had one run valued `1.0` and rendered a line and nothing else,
/// which the "0 px is not a scale bar" assertion below caught rather than reporting `0 == 0` as
/// agreement. That is fixed, and [`a_run_with_one_value_still_gets_a_ruler`] holds it; this run
/// keeps two values so it is testing one thing.
fn run_with_a_moving_box() -> String {
    let frame = |t: f64, far: f64| {
        format!(
            "{{\"t\":{t},\"panels\":[{{\"name\":\"mover\",\"unit\":\"m\",\"kind\":\"paths\",\
             \"bounds\":[0,0,0,{far},0.4,0.4],\"starts\":[0,2],\
             \"vertices\":[0,0.2,0.2, {far},0.2,0.2, 0,0.3,0.2, {far},0.3,0.2],\
             \"values\":[0.0,1.0]}}],\
             \"readings\":[]}}"
        )
    };
    format!(
        "{{\"format\":2,\"title\":\"a moving box\",\"frames\":[{},{}]}}",
        frame(0.0, 1.0),
        frame(1.0, 1.49)
    )
}

#[test]
fn the_bar_is_the_same_length_in_both_frames() {
    let dir = std::env::temp_dir().join("pantometry-scale-bar-test");
    std::fs::create_dir_all(&dir).expect("a temporary directory");
    let path = dir.join("moving.json");
    std::fs::write(&path, run_with_a_moving_box()).expect("the run writes");

    let mut lengths = Vec::new();
    for at in [0usize, 1] {
        let shot = dir.join(format!("frame-{at}.ppm"));
        let Some(pixels) = render(&path, &shot, at) else {
            println!("no GPU adapter on this machine — the scale bar was not rendered or measured");
            return;
        };
        let (width, height, rgb) = image_of(&pixels);
        lengths.push(longest_horizontal_run(width, height, &rgb));
    }
    println!("  the bar is {} px and {} px", lengths[0], lengths[1]);

    // **Exactly equal**, not nearly. The bar is a run of whole pixels from one quad, and the two
    // renders differ only in a panel's stated bounds — nothing about the camera, the aspect or
    // the rasteriser moves between them, so any difference at all is the frame's box reaching the
    // bar. Before this was wired to the run's box the two measured 137 and 92.
    assert_eq!(
        lengths[0], lengths[1],
        "the ruler changed length between two frames of one run, which is the defect this exists \
         to catch: the camera is fitted once over the whole run and the bar was sized by each \
         frame's own contents"
    );
    // And it is a bar rather than a stray glyph stroke: the label's characters are four units
    // wide on a lattice drawn at a fraction of this, and nothing else in the corner is a solid
    // horizontal run.
    assert!(
        lengths[0] > 40,
        "{} px is not a scale bar, so this measured something else",
        lengths[0]
    );
}

/// **A run whose values never change still has a size.**
///
/// The legend returned empty on `hi <= lo` — no scale bar, no clock, nothing — and that reasoning
/// only ever held for the colour bar, which genuinely has no scale to draw. What it produced was
/// a picture with no ruler on it, and a run with one value in it is not exotic: a panel drawn as
/// geometry rather than as a field has exactly that shape.
///
/// This is the run the first draft of the test above wrote by accident, which rendered a line and
/// nothing else.
#[test]
fn a_run_with_one_value_still_gets_a_ruler() {
    let dir = std::env::temp_dir().join("pantometry-flat-scale-test");
    std::fs::create_dir_all(&dir).expect("a temporary directory");
    let path = dir.join("flat.json");
    std::fs::write(
        &path,
        "{\"format\":2,\"title\":\"one value\",\"frames\":[{\"t\":0,\"panels\":[{\
         \"name\":\"flat\",\"unit\":\"m\",\"kind\":\"paths\",\"bounds\":[0,0,0,1,0.4,0.4],\
         \"starts\":[0],\"vertices\":[0,0.2,0.2, 1,0.2,0.2],\"values\":[1.0]}],\
         \"readings\":[]}]}",
    )
    .expect("the run writes");

    let shot = dir.join("flat.ppm");
    let Some(pixels) = render(&path, &shot, 0) else {
        println!("no GPU adapter on this machine — the flat run was not rendered or measured");
        return;
    };
    let (width, height, rgb) = image_of(&pixels);
    let bar = longest_horizontal_run(width, height, &rgb);
    println!("  the bar is {bar} px on a run with one value");
    assert!(
        bar > 40,
        "a run whose colour range is flat drew no scale bar at all ({bar} px): the guard for the \
         colour bar took the ruler with it"
    );
}

/// One frame rendered headless, as the bytes of a binary PPM. `None` when there is no adapter.
fn render(run: &std::path::Path, shot: &std::path::Path, at: usize) -> Option<Vec<u8>> {
    let mut bin = std::env::current_exe().expect("the test binary knows where it is");
    bin.pop();
    if bin.ends_with("deps") {
        bin.pop();
    }
    let out =
        std::process::Command::new(bin.join(format!("pantometry{}", std::env::consts::EXE_SUFFIX)))
            .args([
                "view",
                &run.to_string_lossy(),
                "--snapshot",
                &shot.to_string_lossy(),
                "--frame",
                &at.to_string(),
            ])
            .output()
            .expect("the binary runs");
    let said =
        String::from_utf8_lossy(&out.stdout).into_owned() + &String::from_utf8_lossy(&out.stderr);
    if said.contains("no adapter") || said.contains("no GPU") {
        return None;
    }
    assert!(out.status.success(), "the render refused:\n{said}");
    Some(std::fs::read(shot).expect("the snapshot is written"))
}

/// A binary PPM as `(width, height, rgb)`.
fn image_of(bytes: &[u8]) -> (usize, usize, Vec<u8>) {
    let text = String::from_utf8_lossy(&bytes[..bytes.len().min(64)]).into_owned();
    let mut fields = text.split_ascii_whitespace();
    assert_eq!(fields.next(), Some("P6"), "not a binary PPM");
    let w: usize = fields.next().and_then(|s| s.parse().ok()).expect("a width");
    let h: usize = fields
        .next()
        .and_then(|s| s.parse().ok())
        .expect("a height");
    let max = fields.next().expect("a maximum value");
    let start = text
        .find(max)
        .map(|i| i + max.len() + 1)
        .expect("the header ends");
    (w, h, bytes[start..].to_vec())
}

/// The longest unbroken horizontal run of non-background pixels in the bottom-left corner.
///
/// **Found rather than measured at a fixed place.** The bar's position is a pair of constants in
/// the renderer and a test that hardcoded them would pass while the bar moved off the canvas. In
/// that corner the only solid horizontal run is the bar: the label's glyphs are strokes on a four
/// unit lattice and the colour bar is on the other side.
///
/// The background is the corner pixel rather than a constant, which is the rule `--snapshot` uses
/// and for the reason it gives: the target is sRGB, so the clear colour is stored far brighter
/// than the linear number the pass was handed, and a fixed threshold calls the whole image a line.
fn longest_horizontal_run(width: usize, height: usize, rgb: &[u8]) -> usize {
    let at = |x: usize, y: usize| {
        let i = (y * width + x) * 3;
        [rgb[i], rgb[i + 1], rgb[i + 2]]
    };
    let background = at(0, 0);
    let mut best = 0;
    for y in height / 2..height {
        let mut run = 0;
        for x in 0..width / 2 {
            if at(x, y) == background {
                run = 0;
            } else {
                run += 1;
                best = best.max(run);
            }
        }
    }
    best
}
