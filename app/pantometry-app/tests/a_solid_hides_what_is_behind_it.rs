//! **A line behind a solid is hidden by it, and a line in front is not.**
//!
//! `pantometry view` gained a shape that can be lit — a `surface` panel, triangles rather than
//! lines — and with it the first frame in this workspace where two pipelines have to agree about
//! depth. They share a depth buffer and a `LessEqual` compare, so they should; "should" is not a
//! measurement, and the first picture of an optical bench with a solid lens in it *looked* as
//! though the rays were painted over the glass.
//!
//! # The comparison, and the one that was not a comparison
//!
//! The first version of this rendered one frame and asserted that the far line covered less than
//! half the pixels of the near one. That is not the same claim: the two lines are the same length
//! in space and **not** on the screen, because the near one is nearer and perspective is what a
//! camera does. Measured, 99 against 124 — which the assertion read as a depth test that was not
//! running, and which is what a correct render of that scene looks like.
//!
//! So the frame is rendered **twice**, differing only in whether the wall is there, and the
//! question is what happens to one line's own pixels. The near line is the control: hiding
//! something behind a wall should not change what is in front of it.
//!
//! # Why a hand-written run rather than a scene
//!
//! No domain produces a `surface` panel — the shape arrived for an *example*, which lives outside
//! this workspace — so there is no scene to point at. The run is four vertices and two lines,
//! written here, which is also what makes this a test of the renderer rather than of a physics.
//!
//! # Not under `wasm32`, and not without an adapter
//!
//! It writes a file and renders it. A machine with no GPU skips loudly, as the scene walk beside
//! it does: a skip that says why is a result, and one that says nothing is a suite that has
//! stopped testing anything.

#![cfg(not(target_family = "wasm"))]

/// The two lines, and the wall that is only in one of the two runs.
///
/// Different values so the renderer paints them different colours, which is how the counting below
/// tells one from the other without knowing where either landed.
///
/// **The lines' stated bounds already contain the wall**, so the two runs frame identically. They
/// did not at first, and the camera pulled back to fit the wall: *both* lines shrank, the one in
/// front from 189 pixels to 124, and a control that should have been untouched moved by a third.
/// A comparison whose two halves are photographed from different distances measures the distance.
///
/// **A fifth of a unit each side of the wall, not two.** At two the far line stood far enough
/// behind that perspective carried most of it out past the wall's silhouette: 109 pixels alone and
/// 99 with the wall, a nine per cent difference that is real occlusion and a poor experiment. Just
/// behind, it is covered end to end, and the claim is one a count cannot argue with.
const LINES: &str = r#"{"name":"lines","unit":"m","kind":"paths","bounds":[-1,-1,-0.2,1,1,0.2],
   "starts":[0,2],"vertices":[-0.5,0,0.2, 0.5,0,0.2, -0.5,0,-0.2, 0.5,0,-0.2],
   "values":[0.0,1.0]}"#;
const WALL: &str = r#"{"name":"wall","unit":"n","kind":"surface","bounds":[-1,-1,0,1,1,0],
   "positions":[-1,-1,0, 1,-1,0, 1,1,0, -1,1,0],"triangles":[0,1,2, 0,2,3],
   "values":[1,1,1,1]},"#;

#[test]
fn a_solid_hides_the_line_behind_it_and_not_the_one_in_front() {
    let dir = std::env::temp_dir().join("pantometry-depth-test");
    std::fs::create_dir_all(&dir).expect("a temporary directory");

    let Some(bare) = render(&dir, "bare", &run(false)) else {
        println!("  skipped: this machine has no GPU adapter");
        return;
    };
    let walled = render(&dir, "walled", &run(true)).expect("the second render, on the same device");

    // The two lines are the largest things in the bare frame: the colour bar is two hundred pixels
    // of gradient and no one of its shades covers more than a dozen. Taking the top two by area
    // names them without this test knowing the ramp.
    let mut seen: Vec<([u8; 3], usize)> = tally(&bare).into_iter().collect();
    // The legend's text is the largest non-background thing in a frame this empty -- 1334 pixels
    // of glyph against a line's 167 -- which is why `tally` counts a window in the middle and not
    // the whole image. Nothing but the subject is in it.
    seen.sort_by_key(|(_, n)| std::cmp::Reverse(*n));
    assert!(
        seen.len() >= 2,
        "a frame with two lines in it has at least two colours and this had {}",
        seen.len()
    );
    let (a, b) = (seen[0].0, seen[1].0);
    let after = tally(&walled);
    let (a_before, a_after) = (seen[0].1, after.get(&a).copied().unwrap_or(0));
    let (b_before, b_after) = (seen[1].1, after.get(&b).copied().unwrap_or(0));

    // Which of the two is behind the wall: the one that lost pixels to it.
    let (hidden_before, hidden_after, shown_before, shown_after) =
        if a_before.saturating_sub(a_after) > b_before.saturating_sub(b_after) {
            (a_before, a_after, b_before, b_after)
        } else {
            (b_before, b_after, a_before, a_after)
        };
    println!("  behind the wall: {hidden_before} px alone, {hidden_after} px with the wall");
    println!("  in front of it:  {shown_before} px alone, {shown_after} px with the wall");

    // **The measurement.** The far line sits a fifth of a unit behind a wall two units across
    // and is inside its silhouette from end to end, so with the wall there it is not visible.
    assert!(
        hidden_before > 50,
        "the line behind the wall covered {hidden_before} pixels with nothing in front of \
         it, which is not a line -- this measured nothing before it measured occlusion"
    );
    assert_eq!(
        hidden_after, 0,
        "the line behind the wall still covers {hidden_after} pixels of {hidden_before}, so \
         the two pipelines are not sharing a depth buffer and a solid does not hide what is \
         behind it"
    );
    // The control. Nothing was put in front of this one, so nothing should have happened to it.
    assert_eq!(
        shown_before, shown_after,
        "the line in *front* of the wall changed from {shown_before} pixels to {shown_after} when \
         the wall was added behind it, so what the depth test is doing is not occlusion"
    );
}

/// The run file, with or without the wall.
fn run(walled: bool) -> String {
    format!(
        "{{\"title\":\"depth\",\"format\":2,\"frames\":[{{\"t\":0.0,\"panels\":[{}{}],\
         \"readings\":[]}}]}}",
        if walled { WALL } else { "" },
        LINES
    )
}

/// One render, as `(width, pixels)`. `None` when this machine has no GPU.
fn render(dir: &std::path::Path, name: &str, body: &str) -> Option<Vec<u8>> {
    let run = dir.join(format!("{name}.json"));
    let shot = dir.join(format!("{name}.ppm"));
    std::fs::write(&run, body).expect("the run writes");

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
    Some(std::fs::read(&shot).expect("the snapshot is written"))
}

/// How many pixels each non-background colour covers.
///
/// The background is the corner pixel rather than a constant, which is the rule `--snapshot` uses
/// and for the reason it gives: the target is sRGB, so the clear colour is stored far brighter
/// than the linear number the pass was handed, and a fixed threshold calls the whole image a line.
fn tally(bytes: &[u8]) -> std::collections::BTreeMap<[u8; 3], usize> {
    let text = String::from_utf8_lossy(&bytes[..bytes.len().min(64)]).into_owned();
    let mut fields = text.split_ascii_whitespace();
    assert_eq!(fields.next(), Some("P6"), "not a binary PPM");
    let _w: usize = fields.next().and_then(|s| s.parse().ok()).expect("a width");
    let _h: usize = fields
        .next()
        .and_then(|s| s.parse().ok())
        .expect("a height");
    assert!(_w > 0 && _h > 0, "an empty image");
    let max = fields.next().expect("a maximum value");
    let start = text
        .find(max)
        .map(|i| i + max.len() + 1)
        .expect("the header ends");
    let pixels = &bytes[start..];
    let background = [pixels[0], pixels[1], pixels[2]];
    let mut counts = std::collections::BTreeMap::new();
    // **The middle of the frame.** The scale bar is bottom left, the colour bar bottom right and
    // the clock top left; all three are text, and text is the largest non-background thing in a
    // picture of two lines. Counted over the whole image the glyphs outvoted the subject four to
    // one and this test named them as one of its lines.
    let (x0, x1) = (_w / 5, _w * 4 / 5);
    let (y0, y1) = (_h * 3 / 20, _h * 3 / 4);
    for y in y0..y1 {
        for x in x0..x1 {
            let i = (y * _w + x) * 3;
            let c = [pixels[i], pixels[i + 1], pixels[i + 2]];
            if c != background {
                *counts.entry(c).or_default() += 1;
            }
        }
    }
    counts
}
