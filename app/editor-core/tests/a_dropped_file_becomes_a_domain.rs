//! **A file dropped on the editor becomes a domain, and the scene it lands in still builds.**
//!
//! Only two kinds in the whole scene format take a file: a `block` reads STL geometry through its
//! `parts`, and a `protein` reads a PDB. The other nineteen are numbers. So an asset browser has
//! two entries in its table and a refusal for everything else, and this is that table held against
//! the files that actually ship.
//!
//! # The interesting half is the grid
//!
//! A `block` needs `cells` and `cell_mm`, and a wrong pair is not a cosmetic problem: too coarse
//! and the part voxelises to no cells, which the builder refuses by name; too fine and the scene
//! is a grid nobody can run. `asset_domain` fits the grid to the geometry through the same
//! `fit::propose` that `pantometry fit` uses, so the assertion below is that a dropped part lands
//! on a grid that **holds** it — measured against the STL's own extent, not against a number this
//! file agrees with.
#![cfg(not(target_family = "wasm"))]

use editor_core::{add_domain_json, asset_domain, check};
use pantometry_world::Beside;

/// A scene with an empty `domains`, which is also the branch `add_domain_json` splices into
/// differently.
const EMPTY: &str = r#"{
  "title": "somewhere to drop a thing",
  "duration_s": 1.0,
  "frames": 2,
  "domains": []
}"#;

/// The shipped scenes directory, so a relative `parts/rod.stl` resolves the way it does in a scene.
fn beside() -> Beside {
    Beside::of(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../pantometry-world/scenes/any-scene.json"),
    )
}

fn asset(relative: &str) -> Vec<u8> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../pantometry-world/scenes")
        .join(relative);
    std::fs::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// Drop `relative` into `EMPTY` and give back the scene text.
fn drop_into(text: &str, relative: &str) -> String {
    let domain = asset_domain(relative, &asset(relative))
        .unwrap_or_else(|e| panic!("{relative} should become a domain: {e}"));
    add_domain_json(text, &domain)
        .unwrap_or_else(|e| panic!("{relative}'s domain should splice in: {e}"))
}

/// **A dropped part lands on a grid that holds it, and the scene builds.**
#[test]
fn a_dropped_part_lands_on_a_grid_that_holds_it() {
    let text = drop_into(EMPTY, "parts/rod.stl");
    println!("  {text}");

    let checked = check(&text, &beside());
    assert!(
        checked.error.is_none(),
        "a scene with a dropped part should build: {:?}",
        checked.error
    );
    // The mesh arriving is the part being *found*. A scene that named a file nobody could read
    // would be refused above, but a scene whose part silently contributed nothing is the shape
    // this repository keeps meeting, so the triangles are counted.
    assert_eq!(
        checked.meshes.len(),
        1,
        "the dropped part should be in the scene's geometry"
    );
    assert_eq!(
        checked.meshes[0].triangles.len(),
        92,
        "the rod is 92 triangles"
    );

    // And the grid covers the geometry. The rod is a 24-gon of circumradius 8 mm, 60 long, so its
    // box is 16 x 16 x 60 mm — read here from the closed form the parts library is built on, not
    // from the file and not from the fit.
    let parsed: serde_json::Value = serde_json::from_str(&text).expect("the scene parses");
    let d = &parsed["domains"][0];
    let cell_mm = d["cell_mm"].as_f64().expect("cell_mm is a number");
    let cells: Vec<f64> = d["cells"]
        .as_array()
        .expect("cells is an array")
        .iter()
        .map(|c| c.as_f64().expect("a count"))
        .collect();
    let want = [16.0, 16.0, 60.0];
    for a in 0..3 {
        let span = cells[a] * cell_mm;
        println!(
            "  axis {a}: {} cells x {cell_mm} mm = {span:.3} mm, part {}",
            cells[a], want[a]
        );
        assert!(
            span >= want[a],
            "the grid is {span:.3} mm across axis {a} and the rod is {} mm, so part of it is \
             outside the box it was given",
            want[a]
        );
    }
    assert_eq!(d["kind"], "block", "an STL is geometry, so it is a block");
    assert_eq!(d["name"], "rod", "the name is the file's stem");
    assert_eq!(d["parts"][0]["stl"], "parts/rod.stl");
}

/// **Every shipped part can be dropped**, which is the thing an asset panel will offer.
///
/// Not one example: the panel lists all six, and a part that produced an unbuildable scene would
/// be a tile that looks like the others and fails when used.
#[test]
fn every_shipped_part_can_be_dropped() {
    for name in ["plate", "rod", "wedge", "l-bracket", "heat-sink", "pipe"] {
        let relative = format!("parts/{name}.stl");
        let text = drop_into(EMPTY, &relative);
        let checked = check(&text, &beside());
        let parsed: serde_json::Value = serde_json::from_str(&text).expect("the scene parses");
        let d = &parsed["domains"][0];
        println!(
            "  {name:<10} cells {} cell_mm {}   {}",
            d["cells"],
            d["cell_mm"],
            checked
                .error
                .clone()
                .unwrap_or_else(|| "builds".to_string())
        );
        assert!(
            checked.error.is_none(),
            "{relative} should drop into a scene that builds: {:?}",
            checked.error
        );
        assert_eq!(
            checked.meshes.len(),
            1,
            "{relative} should be in the geometry"
        );
    }
}

/// **A dropped structure names its file and takes its name from it.**
#[test]
fn a_dropped_structure_names_its_file() {
    let text = drop_into(EMPTY, "structures/1CRN.pdb");
    println!("  {text}");
    let parsed: serde_json::Value = serde_json::from_str(&text).expect("the scene parses");
    let d = &parsed["domains"][0];
    assert_eq!(d["kind"], "protein", "a PDB is a structure");
    assert_eq!(d["pdb"], "structures/1CRN.pdb");
    assert_eq!(d["name"], "1CRN", "the name is the file's stem");

    let checked = check(&text, &beside());
    assert!(
        checked.error.is_none(),
        "a scene with a dropped structure should build: {:?}",
        checked.error
    );
}

/// **The `block` template cannot do this**, which is why `asset_domain` builds the domain instead
/// of patching a template.
///
/// A `block` template is a solid of one `material` with no `parts` array, so a pointer aimed at
/// `/parts/0/stl` names nothing and `set_text` refuses it. That refusal is the check working: a
/// version of it that created the missing keys would turn a typo into a silently different scene.
#[test]
fn the_block_template_has_nowhere_to_put_a_file() {
    let with_block = editor_core::add_domain(EMPTY, "block").expect("a block can be added");
    let parsed: serde_json::Value = serde_json::from_str(&with_block).expect("it parses");
    assert!(
        parsed["domains"][0].get("parts").is_none(),
        "the block template is a solid, and if it has grown a `parts` array then this whole \
         function is obsolete and `asset_domain` may be able to use it"
    );

    let why = editor_core::set_text(&with_block, "/domains/0/parts/0/stl", "parts/rod.stl")
        .expect_err("there is nothing at that pointer");
    println!("  {why}");
}

/// **Two drops of one file get two names**, by the same rule two blocks do.
#[test]
fn two_drops_of_one_file_get_two_names() {
    let once = drop_into(EMPTY, "parts/rod.stl");
    let twice = drop_into(&once, "parts/rod.stl");
    let parsed: serde_json::Value = serde_json::from_str(&twice).expect("the scene parses");
    println!(
        "  {} then {}",
        parsed["domains"][0]["name"], parsed["domains"][1]["name"]
    );
    assert_eq!(parsed["domains"][0]["name"], "rod");
    assert_eq!(parsed["domains"][1]["name"], "rod 2");

    // And it still builds, which is the half a rename exists for: two domains of one name is a
    // collision the format refuses.
    let checked = check(&twice, &beside());
    assert!(
        checked.error.is_none(),
        "two of the same part should be a scene: {:?}",
        checked.error
    );
}

/// **A file the format reads for nothing is refused by name.**
///
/// Not defaulted to a block, and not ignored. A drop that quietly did nothing would be the same
/// silence as a panel that renders empty.
#[test]
fn a_file_the_format_reads_for_nothing_is_refused_by_name() {
    for name in ["model.obj", "part.step", "notes.txt", "noextension"] {
        let why = asset_domain(name, b"whatever").expect_err("this is not a kind");
        println!("  {why}");
        assert!(
            why.contains(name) && why.contains(".stl") && why.contains(".pdb"),
            "the refusal should name the file and what it could have been: {why}"
        );
    }
    // And an STL that is not one fails on the geometry rather than on the extension, which is a
    // different sentence and the one a person can act on.
    let why = asset_domain("broken.stl", b"not a solid at all").expect_err("that is not an STL");
    println!("  {why}");
    assert!(why.contains("broken.stl"), "{why}");
}

/// **A path with a quote in it stays valid JSON**, which is the one way a file name can break a
/// scene that a fitted number cannot.
#[test]
fn a_path_with_a_quote_in_it_stays_valid_json() {
    let rod = asset("parts/rod.stl");

    // A quote in the file's own name. The stem keeps it, and the JSON survives it.
    let quoted = "parts/a \"quoted\" part.stl";
    let domain = asset_domain(quoted, &rod).expect("the extension is still .stl");
    println!("  {domain}");
    let parsed: serde_json::Value =
        serde_json::from_str(&domain).expect("a name with a quote must not break the JSON");
    assert_eq!(parsed["parts"][0]["stl"], quoted, "the path round-trips");
    assert_eq!(parsed["name"], "a \"quoted\" part");

    // And a backslash, which is a **path separator** here and not part of a name: on Windows
    // `a\\b.stl` is `b.stl` in a directory called `a`, so the stem is `b`. Stated because the
    // first version of this test expected the whole thing and was wrong about the platform, not
    // about the escaping — which held either way.
    let windows = "parts\\sub\\rod.stl";
    let domain = asset_domain(windows, &rod).expect("the extension is still .stl");
    let parsed: serde_json::Value = serde_json::from_str(&domain).expect("still valid JSON");
    assert_eq!(
        parsed["parts"][0]["stl"], windows,
        "the path round-trips whole"
    );
    assert_eq!(
        parsed["name"], "rod",
        "the directory is not part of the name"
    );
}

/// An ASCII STL of an axis-aligned brick with its low corner at `low`, in millimetres.
fn brick(low: [f64; 3], size: [f64; 3]) -> Vec<u8> {
    let c = |i: usize| {
        [
            if i & 1 == 0 { low[0] } else { low[0] + size[0] },
            if i & 2 == 0 { low[1] } else { low[1] + size[1] },
            if i & 4 == 0 { low[2] } else { low[2] + size[2] },
        ]
    };
    // Each face as two triangles, wound so the normals point out.
    const FACES: [[usize; 4]; 6] = [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
    ];
    let mut s = String::from("solid brick\n");
    for f in FACES {
        for (a, b, cc) in [(f[0], f[1], f[2]), (f[0], f[2], f[3])] {
            s.push_str("  facet normal 0 0 0\n    outer loop\n");
            for p in [c(a), c(b), c(cc)] {
                s.push_str(&format!("      vertex {} {} {}\n", p[0], p[1], p[2]));
            }
            s.push_str("    endloop\n  endfacet\n");
        }
    }
    s.push_str("endsolid brick\n");
    s.into_bytes()
}

/// Two bricks of one size held in memory, the way the web shell holds a dropped file: one drawn from
/// the origin, one drawn about it.
fn bricks() -> pantometry_world::Uploaded {
    pantometry_world::Uploaded::new()
        .with("at-origin.stl", brick([0.0, 0.0, 0.0], [20.0, 20.0, 20.0]))
        .with(
            "centred.stl",
            brick([-10.0, -10.0, -10.0], [20.0, 20.0, 20.0]),
        )
}

/// How many cells the scene's one part filled, from the build's own note.
fn filled(checked: &editor_core::Checked) -> usize {
    checked
        .notes
        .iter()
        .find_map(|n| n.split(" filled ").nth(1)?.split(' ').next()?.parse().ok())
        .unwrap_or_else(|| panic!("no rasterisation note in {:?}", checked.notes))
}

/// The smallest and largest coordinate of every triangle drawn for the scene's one part.
fn drawn_box(checked: &editor_core::Checked) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for t in &checked.meshes[0].triangles {
        for p in t {
            for a in 0..3 {
                lo[a] = lo[a].min(p[a]);
                hi[a] = hi[a].max(p[a]);
            }
        }
    }
    (lo, hi)
}

/// **A part drawn around the origin lands on a grid around it**, and fills exactly what the same
/// part drawn from the origin fills.
///
/// Until `grid_origin` existed a `block`'s grid started at the origin of the parts' coordinates
/// and this brick, spanning `-10..10` mm, was refused for reaching outside it — this test pinned
/// that refusal. The drop now writes `"grid_origin": "parts"`, which the builder resolves from the
/// meshes themselves.
///
/// # What is checked, and against what
///
/// **The count is a closed form.** A cube that fills its grid fills every cell, so the cells filled
/// must be the product of the three counts the scene states — whatever grid the fit chose, and for
/// both bricks. A grid that started anywhere but the brick's corner would cut a slab off one face
/// and the product would not be reached.
///
/// **The surface is drawn on the cells.** The editor draws a part in its grid's frame, which runs
/// from `0` to `cells × cell_mm`. Drawn where the file put it, the centred brick would sit 10 mm
/// down each axis from the cells it became — the tolerance below is a picometre, ten orders under
/// that, and far above the rounding of `±0.01 m`.
#[test]
fn a_part_drawn_around_the_origin_is_dropped_on_a_grid_around_it() {
    use pantometry_world::Parts;
    let files = bricks();
    let drop = |name: &str| {
        let bytes = files.bytes(name).expect("held in memory");
        let domain = asset_domain(name, &bytes).expect("a brick is an STL");
        add_domain_json(EMPTY, &domain).expect("it splices")
    };

    for name in ["at-origin.stl", "centred.stl"] {
        let text = drop(name);
        let checked = check(&text, &files);
        assert!(
            checked.error.is_none(),
            "{name} should build on a grid around it: {:?}",
            checked.error
        );

        let parsed: serde_json::Value = serde_json::from_str(&text).expect("the scene parses");
        let d = &parsed["domains"][0];
        assert_eq!(d["grid_origin"], "parts", "the drop says where its grid is");
        let cells: Vec<usize> = d["cells"]
            .as_array()
            .expect("cells")
            .iter()
            .map(|c| c.as_u64().expect("a count") as usize)
            .collect();
        let cell_m = d["cell_mm"].as_f64().expect("cell_mm") * 1e-3;
        let every = cells[0] * cells[1] * cells[2];
        let got = filled(&checked);
        println!(
            "  {name:<14} cells {cells:?} at {} mm, filled {got} of {every}",
            cell_m * 1e3
        );
        assert_eq!(
            got, every,
            "{name}: a cube filling its grid fills every cell, and {got} of {every} were"
        );

        let (lo, hi) = drawn_box(&checked);
        for a in 0..3 {
            let edge = cells[a] as f64 * cell_m;
            assert!(
                lo[a].abs() < 1e-12 && (hi[a] - edge).abs() < 1e-12,
                "{name}: drawn from {} to {} m on axis {a}, and its cells run from 0 to {edge}",
                lo[a],
                hi[a]
            );
        }
    }
}

/// **Without `grid_origin` the refusal still stands, and still names both boxes.**
///
/// The key is opt-in, so every block written before it keeps its grid where it was. This is the
/// same scene as above with the one key taken out, and it has to be refused the way it always was:
/// a part with its corner missing runs and audits and answers about a different shape.
#[test]
fn without_grid_origin_a_part_around_the_origin_is_still_refused() {
    use pantometry_world::Parts;
    let files = bricks();
    let bytes = files.bytes("centred.stl").expect("held in memory");
    let domain = asset_domain("centred.stl", &bytes).expect("a brick is an STL");
    let text = add_domain_json(EMPTY, &domain)
        .expect("it splices")
        .replace("\"grid_origin\": \"parts\", ", "");
    assert!(
        !text.contains("grid_origin"),
        "the key was meant to be removed, and the test proves nothing if it was not"
    );
    let why = check(&text, &files)
        .error
        .expect("a part drawn about the origin reaches outside a grid that starts there");
    println!("  {why}");
    assert!(
        why.contains("cut off") && why.contains("refused rather than cropped"),
        "the refusal should say both boxes and that it did not crop: {why}"
    );
}

/// **What `pantometry fit` recommends for a part off the origin builds, pasted as it comes.**
///
/// It did not. `fit::propose` measured every candidate grid from the parts' lowest corner and
/// `scene_fragment` left the corner out, on the stated belief that the rasteriser "places every
/// part against the grid's own corner" — it does not, the builder puts that corner at the origin.
/// The CLI printed `from (-10.0, -10.0, -10.0) mm` above its own recommendation, and the
/// recommendation was refused. The web shell's `fit` is the same function, so this is both.
#[test]
fn what_fit_recommends_for_a_part_off_the_origin_builds() {
    let files = bricks();
    let only = pantometry_world::Uploaded::new().with(
        "centred.stl",
        brick([-10.0, -10.0, -10.0], [20.0, 20.0, 20.0]),
    );
    let json = editor_core::fit(&only, 2_000_000, "aluminium").expect("a brick has a grid");
    let v: serde_json::Value = serde_json::from_str(&json).expect("fit answers in JSON");
    let fragment = v["fragment"].as_str().expect("a recommended fragment");
    println!("{fragment}");
    assert!(
        fragment.contains("\"grid_origin\": \"parts\""),
        "the fragment should say where the grid it measured starts"
    );
    let text = format!(
        "{{ \"title\": \"as fit said\", \"duration_s\": 1.0, \"frames\": 2, \"domains\": [ \
         {{ \"kind\": \"block\", \"name\": \"cube\", \"initial_c\": 20.0,\n{fragment} }} ] }}"
    );
    let checked = check(&text, &files);
    assert!(
        checked.error.is_none(),
        "the grid fit recommended should be a grid that builds: {:?}",
        checked.error
    );
}

/// **A `grid_origin` on a block with no parts is refused**, because it would change nothing.
///
/// A key the format reads and then ignores is a key somebody believes is doing something.
#[test]
fn a_grid_origin_with_no_parts_is_refused() {
    let text = r#"{ "title": "a grid among no parts", "duration_s": 1.0, "frames": 2,
      "domains": [ { "kind": "block", "name": "solid", "cells": [2, 2, 2], "cell_mm": 1.0,
                     "initial_c": 20.0, "grid_origin": "parts" } ] }"#;
    let why = check(text, &bricks())
        .error
        .expect("grid_origin without parts is refused");
    println!("  {why}");
    assert!(
        why.contains("grid_origin") && why.contains("solid"),
        "the refusal should name the key and the block: {why}"
    );
}
