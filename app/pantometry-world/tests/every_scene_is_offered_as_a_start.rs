//! The chooser's starting points, against the scenes they are.
//!
//! `presets::PRESETS` embeds every shipped scene so the editor can offer it without knowing where
//! the repository is. That list is written by hand, which is a list that goes wrong quietly
//! between a scene being added and somebody remembering — the same failure `every_domain_has_a_
//! template` exists for, and the same answer: compare it with the directory **in both
//! directions**.
//!
//! What is *not* checked here is which area a scene was filed under. That is a judgement — a power
//! module is heat in a solid and is offered under power electronics because that is what somebody
//! opening it is doing — and a test that agreed with the judgement would only be repeating it.
//! What is checked is that every area is named, non-empty, and reached by something.
//!
//! # Not under `wasm32`
//!
//! It reads the scenes directory. The presets themselves are embedded and work anywhere.

#![cfg(not(target_family = "wasm"))]

use pantometry_world::presets::{AREAS, PRESETS};
use pantometry_world::{OnDisk, Scene, World};
use std::collections::BTreeSet;

/// The scene files on disk.
fn shipped() -> BTreeSet<String> {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("scenes");
    let names: BTreeSet<String> = std::fs::read_dir(&dir)
        .expect("the scenes directory reads")
        .filter_map(Result::ok)
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".json"))
        .collect();
    assert!(
        names.len() >= 28,
        "only {} scenes in {}",
        names.len(),
        dir.display()
    );
    names
}

/// **Every shipped scene is offered, and every offer is a shipped scene.**
#[test]
fn every_shipped_scene_is_offered() {
    let on_disk = shipped();
    let offered: BTreeSet<String> = PRESETS.iter().map(|p| p.file.to_string()).collect();
    assert_eq!(
        offered.len(),
        PRESETS.len(),
        "a scene is offered twice: {:?}",
        PRESETS.iter().map(|p| p.file).collect::<Vec<_>>()
    );
    let missing: Vec<_> = on_disk.difference(&offered).collect();
    let ghosts: Vec<_> = offered.difference(&on_disk).collect();
    assert!(
        missing.is_empty(),
        "these scenes ship and are not offered: {missing:?}"
    );
    assert!(
        ghosts.is_empty(),
        "these are offered and do not ship: {ghosts:?}"
    );
}

/// **Every preset's title and kinds are the scene's own**, not a second copy of them.
///
/// They are generated from the file, and this is what says they still match it: a scene renamed in
/// its own JSON and not here would offer a title nothing else in the workspace uses.
#[test]
fn a_presets_title_and_kinds_are_the_scenes_own() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("scenes");
    for p in PRESETS {
        let text = std::fs::read_to_string(dir.join(p.file)).expect("the scene reads");
        let value: serde_json::Value = serde_json::from_str(&text).expect("the scene parses");
        assert_eq!(
            value["title"].as_str(),
            Some(p.title),
            "{}: the offered title is not the scene's",
            p.file
        );
        let kinds: BTreeSet<&str> = value["domains"]
            .as_array()
            .into_iter()
            .flatten()
            .filter_map(|d| d["kind"].as_str())
            .collect();
        let offered: BTreeSet<&str> = p.kinds.iter().copied().collect();
        assert_eq!(
            kinds, offered,
            "{}: the offered kinds are not the scene's",
            p.file
        );
        // And the text is the file, which `include_str!` gives by construction — checked anyway,
        // because a path typed one character wrong embeds a different scene and still compiles.
        assert_eq!(
            p.json, text,
            "{}: the embedded text is another file",
            p.file
        );
    }
}

/// **Every preset builds**, which is what makes it a starting point rather than a sample.
///
/// One names a file beside itself — an STL for a designed part — and says so, because a preset
/// opened as an unsaved project has nowhere to resolve that from. It builds here because the test
/// runs beside the scenes.
#[test]
fn every_preset_builds_and_says_when_it_needs_a_file() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("scenes");
    let mut needs = Vec::new();
    for p in PRESETS {
        let scene: Scene = serde_json::from_str(p.json)
            .unwrap_or_else(|e| panic!("{}: the preset does not parse: {e}", p.file));
        let beside = pantometry_world::Beside::of(dir.join(p.file));
        World::build_with(scene, &beside)
            .unwrap_or_else(|e| panic!("{}: the preset does not build: {e}", p.file));
        if p.needs_a_part {
            needs.push(p.file);
        }
        // A scene that reaches for a file says so, and one that does not, does not.
        let reaches = p.json.contains("\"stl\"");
        assert_eq!(
            reaches,
            p.needs_a_part,
            "{}: `needs_a_part` is {} and the scene {} a part",
            p.file,
            p.needs_a_part,
            if reaches { "names" } else { "does not name" }
        );
    }
    assert_eq!(
        needs,
        ["29-a-designed-bracket-becomes-cells.json"],
        "the set of presets that need a file beside them has changed"
    );
}

/// **Every area is named, and every area has something in it.**
///
/// An area with nothing under it is a heading a reader can open and find empty, which is the
/// shape of a control that looks like a feature. `OnDisk` is unused here and named so the import
/// stays honest about what this file needs.
#[test]
fn every_area_is_named_and_none_is_empty() {
    let _ = OnDisk;
    let named: BTreeSet<&str> = AREAS.iter().map(|(k, _, _)| *k).collect();
    assert_eq!(named.len(), AREAS.len(), "an area key is used twice");
    for (key, name, about) in AREAS {
        assert!(!name.is_empty() && !about.is_empty(), "{key} has no words");
        let n = PRESETS.iter().filter(|p| p.area == key).count();
        assert!(n > 0, "{key} ({name}) is offered and holds nothing");
    }
    let used: BTreeSet<&str> = PRESETS.iter().map(|p| p.area).collect();
    let unnamed: Vec<_> = used.difference(&named).collect();
    assert!(
        unnamed.is_empty(),
        "presets filed under unnamed areas: {unnamed:?}"
    );
}

/// **A preset has a picture exactly when its scene has something to draw.**
///
/// `thumb: None` is what the chooser turns into "no picture — this one reports readings, not
/// places", which is a *claim about the scene*. Nothing checked it. A generator that failed on
/// one scene — the render refused, the file was not written, the path was typed wrong — would
/// have made that sentence appear under a scene that draws perfectly well, and the screen would
/// have looked deliberate.
///
/// # This is a pin on a generated file, and not corroboration
///
/// It was written as one — "two lists from two different sources, one from the tiles on disk and
/// one from a GPU render" — and that was wrong twice over, which is worth keeping because the
/// sentence was persuasive.
///
/// Both lists trace to one `eprintln!`. `pantometry view` says `this run has no panels at all`
/// (`view.rs`), `make.py` greps for that string to decide `thumb: None`, and
/// `every_shipped_scene_puts_something_on_the_canvas` greps for the same string at test time. The
/// tiles on disk *are* that check, frozen at generation. And the render is not where it comes
/// from: the no-panel branch runs before an adapter is ever requested, so it reports the same
/// three on a machine with no GPU — which is CI's case, stated in `ci.yml`.
///
/// What this test is, then, is a pin: a compile-time constant against three literals. It cannot
/// fail on its own, and it does not have to. It fails when somebody regenerates after a scene
/// lost its geometry, which is exactly when the screen would otherwise start claiming that scene
/// reports readings. The live check is the other file's, and it is live on any runner because of
/// where that branch sits.
#[test]
fn a_preset_has_a_picture_unless_its_scene_has_nothing_to_draw() {
    let without: Vec<&str> = PRESETS
        .iter()
        .filter(|p| p.thumb.is_none())
        .map(|p| p.file)
        .collect();
    assert_eq!(
        without,
        [
            "11-motor-thermal-network.json",
            "12-winding-heats-a-motor.json",
            "13-winding-that-heats-itself.json"
        ],
        "the set of presets offered without a picture has changed"
    );
}

/// **Every tile decodes, is the size the chooser draws it at, and is not empty.**
///
/// The chooser hands the bytes to `image` and asks for a texture at 240x156; a PNG of another size
/// is silently *stretched*, and a PNG of nothing but the background is a grey rectangle that reads
/// as a scene which failed rather than as a scene with little in it.
///
/// # The two sparse ones are pinned by rank, not by a threshold
///
/// "Not background" is counted the way the renderer's own snapshot count is: pixels differing from
/// the corner, out of 37 440. The bottom of the spread is very low — `07-bouncing-ball` lights
/// **5** and `06-orbits` **29**, against 10 704 for a room mode. Both are a handful of point
/// bodies, and a point is a cross about 1.2% of the frame across whatever the camera is fitted to,
/// so magnifying the crop does not rescue them: measured at a twentyfold cap, the bouncing ball is
/// a small cross in the middle of a grey field, which is a picture of the *marker*.
/// `crop_to_content` caps at three for that reason and this is the cost of it.
///
/// This was a `SPARSE = 90` threshold, and the third-lowest tile lights **94** — four pixels of
/// 37 440 of headroom, so a hair's movement in one scene's render would have fired an assertion
/// saying "the set of nearly-empty tiles has changed", which would have been a lie about which
/// thing moved. Rank has no such edge: the two lowest are named with their counts, and the third
/// is required to be clear of them.
#[test]
fn every_tile_decodes_to_something() {
    /// What the chooser draws a tile at, and what `pantometry view --thumbnail` writes.
    const TILE: (u32, u32) = (240, 156);
    let mut lit_by: Vec<(u64, &str)> = Vec::new();
    for p in PRESETS {
        let Some(bytes) = p.thumb else { continue };
        let img = image::load_from_memory(bytes)
            .unwrap_or_else(|e| panic!("{}: the tile does not decode: {e}", p.file))
            .to_rgb8();
        assert_eq!(
            (img.width(), img.height()),
            TILE,
            "{}: the tile is {}x{} and the chooser draws it {}x{}",
            p.file,
            img.width(),
            img.height(),
            TILE.0,
            TILE.1
        );
        let bg = *img.get_pixel(0, 0);
        let lit = img.pixels().filter(|px| **px != bg).count() as u64;
        assert!(
            lit > 0,
            "{}: every pixel is the background — the tile is blank",
            p.file
        );
        lit_by.push((lit, p.file));
    }
    // A run in which every preset had `thumb: None` would light nothing and pass the loop above.
    assert_eq!(
        lit_by.len(),
        27,
        "27 presets carry a tile, and this run found {}",
        lit_by.len()
    );
    lit_by.sort_unstable();
    assert_eq!(
        &lit_by[..2],
        &[(5, "07-bouncing-ball.json"), (29, "06-orbits.json")],
        "the two sparsest tiles are not the two that were measured"
    );
    assert!(
        lit_by[2].0 >= 90,
        "the third-sparsest tile is {} at {} pixels, down among the two that are a single body",
        lit_by[2].1,
        lit_by[2].0
    );
}

/// **No two tiles are the same picture** — and four pairs of them are.
///
/// A chooser whose pictures do not tell two scenes apart is a chooser with no pictures, and the
/// per-tile checks above cannot see it: each of these decodes, is the right size, and is far from
/// blank. Measured over the 27 committed tiles: **23 distinct images**.
///
/// The collisions are explicable rather than a rendering fault, which is why they are pinned
/// rather than fixed here. A tile's shade is the sample's value normalised over the run's own
/// range, and the frame is the last one — so three heat scenes in the same block, driven to the
/// same steady state by the same heater, produce the same *normalised* field whatever the material
/// between them is. `--frame` is not the culprit: rendering `20-melting-a-block-of-ice` at frames
/// 0, 6 and 11 gives three different pictures, and 21 and 22 reach the shared one by frame 6 while
/// 20 does not until 11.
///
/// So this is honest and it is still bad on the screen. What is asserted is the *set* of
/// collisions: a new one is a scene that stopped being distinguishable from another, and that is
/// worth a failure rather than something to notice while scrolling.
#[test]
fn the_tiles_tell_the_scenes_apart_or_say_which_they_do_not() {
    let mut by_bytes: std::collections::BTreeMap<&[u8], Vec<&str>> =
        std::collections::BTreeMap::new();
    for p in PRESETS {
        if let Some(bytes) = p.thumb {
            by_bytes.entry(bytes).or_default().push(p.file);
        }
    }
    let mut same: Vec<Vec<&str>> = by_bytes.values().filter(|v| v.len() > 1).cloned().collect();
    same.sort_unstable();
    assert_eq!(
        same,
        [
            vec!["08-atoms-crystal.json", "09-atoms-liquid.json"],
            vec![
                "20-melting-a-block-of-ice.json",
                "21-a-wax-thermal-buffer.json",
                "22-wax-in-an-aluminium-matrix.json"
            ],
            vec![
                "24-a-power-module-junction-to-ambient.json",
                "25-what-140-kelvin-does-to-the-solder.json"
            ],
        ],
        "the set of scenes whose tiles are the same picture has changed"
    );
    assert_eq!(
        by_bytes.len(),
        23,
        "27 tiles, and {} of them are distinct pictures",
        by_bytes.len()
    );
}
