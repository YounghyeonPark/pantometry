//! **A run carries the numbers the simulation took, and the viewer drew none of them.**
//!
//! `Frame::readings` is in the wire format, the editor draws them, and `view.rs` did not contain
//! the word. So a window showing a protein closing over a molecule could not say how far it had
//! closed, with the number sitting in the file it had open — and a snapshot of any run was a
//! picture of its geometry with its measurements dropped.
//!
//! The controls had the same shape: `drag to rotate, scroll to zoom, space to play, left/right to
//! scrub` went to **stdout**, where a person looking at the picture is not. That is most of what
//! there is to know about driving this, and a screenshot without it teaches nothing about the
//! application it is a screenshot of.
//!
//! # What these check, and what they cannot
//!
//! That the text is *drawn*, and that it stays on the canvas. Both failures here were of the kind
//! that looks like success: nothing drawn reads as a run with nothing to say, and a value pushed
//! past the right edge reads as a value that is simply short. Neither is visible in a picture
//! without the file beside it.
//!
//! They do not check that the words are the right words — the glyph table's own tests cover what
//! a character draws, and no pixel comparison here would say more than a person reading it.
//!
//! # Not under `wasm32`, and not without an adapter
//!
//! A machine with no GPU skips loudly, as its neighbours do.

#![cfg(not(target_family = "wasm"))]

/// One frame, with whatever readings are handed in.
///
/// The geometry is the same in every run here — two short lines with two values, which is the
/// least that draws a legend at all — so the only thing that varies between renders is the text.
fn run_with(readings: &str) -> String {
    format!(
        "{{\"format\":2,\"title\":\"a measured run\",\"frames\":[{{\"t\":0,\"panels\":[{{\
         \"name\":\"thing\",\"unit\":\"m\",\"kind\":\"paths\",\"bounds\":[0,0,0,1,0.4,0.4],\
         \"starts\":[0,2],\
         \"vertices\":[0,0.2,0.2, 1,0.2,0.2, 0,0.3,0.2, 1,0.3,0.2],\"values\":[0.0,1.0]}}],\
         \"readings\":[{readings}]}}]}}"
    )
}

/// `{"domain":..,"label":..,"value":..,"unit":..}`, as the wire format spells one.
fn reading(label: &str, value: f64) -> String {
    format!("{{\"domain\":\"d\",\"label\":\"{label}\",\"value\":{value},\"unit\":\"A\"}}")
}

#[test]
fn the_readings_are_drawn_and_an_empty_list_draws_none_of_them() {
    let dir = std::env::temp_dir().join("pantometry-readings-test");
    std::fs::create_dir_all(&dir).expect("a temporary directory");

    let with = run_with(&format!(
        "{},{}",
        reading("from the closed form", 4.2926),
        reading("excursion", 84.632)
    ));
    let Some(lit_with) = lit_in_the_text_band(&dir, "with", &with) else {
        println!("no GPU adapter on this machine — the readings were not rendered or measured");
        return;
    };
    let lit_without =
        lit_in_the_text_band(&dir, "without", &run_with("")).expect("the adapter is still here");

    println!("  the band carries {lit_with} lit pixels with readings, {lit_without} without");
    // **Strictly more**, and by a lot: two lines of text against none. The clock and the control
    // hints are in this band in both runs, which is why the comparison is against the same render
    // with the readings removed rather than against zero.
    assert!(
        lit_with > lit_without + 200,
        "two readings added {} lit pixels to the band, which is not two lines of text — the \
         viewer drew none of them for years and that looked exactly like a run with nothing to say",
        lit_with - lit_without
    );
}

#[test]
fn a_long_label_does_not_push_its_value_off_the_canvas() {
    let dir = std::env::temp_dir().join("pantometry-readings-test");
    std::fs::create_dir_all(&dir).expect("a temporary directory");

    // Longer than anything this workspace produces — the longest real one is `the furthest
    // residue has moved`, at thirty. The column is placed from the widest label, so without a
    // clamp a label like this puts its number past the right-hand edge, where a value that is not
    // drawn looks exactly like a value that is short.
    let run = run_with(&reading(
        "a label considerably longer than any this workspace has ever produced",
        18.5778,
    ));
    let shot = dir.join("long.ppm");
    let path = dir.join("long.json");
    std::fs::write(&path, run).expect("the run writes");
    let Some(pixels) = render(&path, &shot) else {
        println!("no GPU adapter on this machine — the long label was not rendered or measured");
        return;
    };
    let (width, height, rgb) = image_of(&pixels);

    // Nothing lit in the last few columns, anywhere: a glyph that reached the edge was clipped by
    // it, and one that went past is simply gone.
    let background = pixel(&rgb, width, 0, 0);
    let mut at_the_edge = 0;
    for y in 0..height {
        for x in width - 3..width {
            if pixel(&rgb, width, x, y) != background {
                at_the_edge += 1;
            }
        }
    }
    println!("  {at_the_edge} lit pixels in the last three columns");
    assert_eq!(
        at_the_edge, 0,
        "a label of 69 characters pushed something to the canvas edge, so a value was clipped or \
         lost there"
    );
}

/// Lit pixels in the top band, where the clock, the readings and the control hints are drawn.
fn lit_in_the_text_band(dir: &std::path::Path, name: &str, body: &str) -> Option<usize> {
    let path = dir.join(format!("{name}.json"));
    let shot = dir.join(format!("{name}.ppm"));
    std::fs::write(&path, body).expect("the run writes");
    let pixels = render(&path, &shot)?;
    let (width, height, rgb) = image_of(&pixels);
    let background = pixel(&rgb, width, 0, 0);
    // The top sixth, which holds the legend's text and none of the scene: the camera fits the
    // subject to the middle of the canvas.
    let mut lit = 0;
    for y in 0..height / 6 {
        for x in 0..width {
            if pixel(&rgb, width, x, y) != background {
                lit += 1;
            }
        }
    }
    Some(lit)
}

fn pixel(rgb: &[u8], width: usize, x: usize, y: usize) -> [u8; 3] {
    let i = (y * width + x) * 3;
    [rgb[i], rgb[i + 1], rgb[i + 2]]
}

/// One frame rendered headless, as the bytes of a binary PPM. `None` when there is no adapter.
fn render(run: &std::path::Path, shot: &std::path::Path) -> Option<Vec<u8>> {
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
