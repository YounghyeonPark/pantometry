//! **The model tree and its forms, held against the thirty-four shipped scenes.**
//!
//! The settings pane replaces typing JSON for most of what a person does to a scene, so the ways
//! it could be quietly wrong are the ways a file could be quietly misrepresented: a value no form
//! shows, a value two forms both claim, an add that leaves the file different from before after
//! it is taken back, a menu that offers a word the format refuses. Each of those is checked here
//! over every shipped scene, and none of them needs a window.
//!
//! # Not under `wasm32`
//!
//! It reads the scenes off a disk.

#![cfg(not(target_family = "wasm"))]

use editor_core::{form, model, Row, Value};
use pantometry_world::Scene;

/// Every shipped scene, as `(file name, text)`.
fn shipped() -> Vec<(String, String)> {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("editor-core has a parent")
        .join("pantometry-world/scenes");
    let mut out = Vec::new();
    for entry in std::fs::read_dir(&dir).expect("the scenes directory reads") {
        let path = entry.expect("a directory entry").path();
        if path.extension().is_some_and(|x| x == "json") {
            let name = path.file_name().unwrap().to_string_lossy().to_string();
            out.push((name, std::fs::read_to_string(&path).expect("a scene reads")));
        }
    }
    out.sort();
    // An empty directory would make every loop below vacuously true.
    assert!(
        out.len() >= 34,
        "only {} scenes in {}",
        out.len(),
        dir.display()
    );
    out
}

fn one(file: &str) -> String {
    shipped()
        .into_iter()
        .find(|(n, _)| n == file)
        .unwrap_or_else(|| panic!("no shipped scene {file}"))
        .1
}

/// Every row of `text`'s model tree, as the outliner gets them: from `editor_core::tree`, the
/// function the shell calls, and not rebuilt here from the model it is built from.
fn rows_of(text: &str) -> Vec<String> {
    let checked = editor_core::Checked {
        model: model(text),
        ..Default::default()
    };
    editor_core::tree(&checked, None, 0)
        .nodes
        .into_iter()
        .map(|n| n.path)
        .collect()
}

/// Every leaf of `value` as an RFC 6901 pointer, walked here rather than taken from the crate.
///
/// Numbers, strings, booleans and nulls, at any depth, except under `kind` and `name` — the two
/// keys the forms hold back, by a rule `editor_core::edit` states. A boolean is in this list and in
/// no form, because no form has a widget for one, so the day the format grows a flag this fails
/// instead of the flag quietly having no home.
fn leaves(value: &serde_json::Value, at: &str, out: &mut Vec<String>) {
    match value {
        serde_json::Value::Object(o) => {
            for (k, m) in o {
                if k == "kind" || k == "name" {
                    continue;
                }
                let key = k.replace('~', "~0").replace('/', "~1");
                leaves(m, &format!("{at}/{key}"), out);
            }
        }
        serde_json::Value::Array(a) => {
            for (j, m) in a.iter().enumerate() {
                leaves(m, &format!("{at}/{j}"), out);
            }
        }
        _ => out.push(at.to_string()),
    }
}

/// How many values `text` states, after checking each is in exactly one row's form.
fn homes(file: &str, text: &str) -> usize {
    let root: serde_json::Value = serde_json::from_str(text).expect("a scene is JSON");
    let mut expected = Vec::new();
    leaves(&root, "", &mut expected);
    let mut found: std::collections::BTreeMap<String, String> = Default::default();
    for row in rows_of(text) {
        for v in form(text, &row).values() {
            if let Some(first) = found.insert(v.pointer.clone(), row.clone()) {
                panic!("{file}: {} is in two forms, {first} and {row}", v.pointer);
            }
        }
    }
    expected.sort();
    let found: Vec<String> = found.into_keys().collect();
    assert_eq!(
        found, expected,
        "{file}: the forms and the scene's values differ"
    );
    expected.len()
}

/// **Every value the scene states is in exactly one form.**
///
/// The tree splits a scene across rows: a domain keeps what it states outside its lists and each
/// entry takes its own. So the union of every row's form has to be every value in the file. A
/// value missing would be part of the file no form can change; a value in two forms would be two
/// widgets for one number, which drift when one of them is stale.
///
/// Both halves are independent of what they check. The rows are `tree`'s, so a row the outliner
/// lost is a value with no home here too; and the values are a walk of the JSON written in this
/// file, so a value the crate's own walk dropped is not dropped from the expectation as well.
#[test]
fn every_value_has_one_home() {
    let mut checked = 0;
    for (file, text) in shipped() {
        checked += homes(&file, &text);
    }
    // 865 on 2026-11-17. A floor rather than the count, so a new scene does not fail this; a
    // number far below it would mean the walk stopped walking.
    println!("  {checked} values, each in one form");
    assert!(
        checked >= 800,
        "only {checked} values checked across the scenes"
    );
}

/// **A name with a slash in it keeps its values.** The format refuses only an empty material
/// name, so `FR-4/Cu` is a material — and its row's path has one more `/` than a group entry's,
/// which split as a path left that row's form empty and its values with no home anywhere.
#[test]
fn a_name_with_a_slash_keeps_its_values() {
    let text = r#"{ "title": "t", "duration_s": 1.0, "frames": 2,
  "materials": { "FR-4/Cu": { "density": 1850.0, "specific_heat": 1100.0 } },
  "domains": [] }"#;
    assert!(rows_of(text).contains(&"/scene/materials/FR-4/Cu".to_string()));
    assert_eq!(
        homes("a slash", text),
        5,
        "title, duration, frames and the two properties"
    );
}

/// **An entry added and then taken out leaves the file as it was, byte for byte**, for every
/// list in every shipped scene — and so does a copy of the first and of the last entry.
///
/// And in between, the scene still parses: an entry the format refused a key of would be an add
/// button that breaks the file. Whether it *builds* is the scene's business — a dose into a
/// compartment this scene does not have is refused by name, which is the check doing its job.
#[test]
fn an_entry_added_or_copied_and_removed_leaves_the_file_as_it_was() {
    let (mut lists, mut adds) = (0, 0);
    for (file, text) in shipped() {
        // Every shipped scene parses, so the parse check below is never skipped.
        serde_json::from_str::<Scene>(&text).unwrap_or_else(|e| panic!("{file}: {e}"));
        for d in model(&text).domains {
            for l in &d.lists {
                let n = l.entries.len();
                lists += 1;
                let mut changed = Vec::new();
                if pantometry_world::templates::entry(&d.kind, &l.key).is_some() {
                    let added = editor_core::add_entry(&text, &d.name, &l.key)
                        .unwrap_or_else(|e| panic!("{file} {} {}: {e}", d.name, l.key));
                    changed.push(("added", added, n));
                    adds += 1;
                }
                if n > 0 {
                    for j in [0, n - 1] {
                        let copied = editor_core::duplicate_entry(&text, &d.name, &l.key, j)
                            .unwrap_or_else(|e| panic!("{file} {} {}[{j}]: {e}", d.name, l.key));
                        changed.push(("copied", copied, j + 1));
                    }
                }
                for (how, after, at) in changed {
                    let grown = model(&after);
                    let list = grown
                        .domains
                        .iter()
                        .find(|x| x.name == d.name)
                        .and_then(|x| x.lists.iter().find(|y| y.key == l.key))
                        .unwrap_or_else(|| panic!("{file} {}: the list went missing", d.name));
                    assert_eq!(
                        list.entries.len(),
                        n + 1,
                        "{file} {} {} {how}",
                        d.name,
                        l.key
                    );
                    serde_json::from_str::<Scene>(&after).unwrap_or_else(|e| {
                        panic!(
                            "{file} {} {}: {how}, the scene stops parsing: {e}",
                            d.name, l.key
                        )
                    });
                    let back = editor_core::remove_entry(&after, &d.name, &l.key, at)
                        .unwrap_or_else(|e| panic!("{file} {} {}: {e}", d.name, l.key));
                    assert_eq!(back, text, "{file} {} {}: {how} and removed", d.name, l.key);
                }
            }
        }
    }
    println!("  {lists} lists added to, copied in and restored");
    // Every list in a shipped scene has a starting entry, so every one went through the add.
    assert_eq!(adds, lists, "a list was copied in but never added to");
    assert!(lists >= 30, "only {lists} lists across the scenes");
}

/// **A list a domain does not hold yet is started, and taking its one entry out leaves `[]`**,
/// which the format reads as the list being absent. Every list the format keeps an entry for,
/// started in every shipped domain of its kind that lacks it.
///
/// Held to bytes: the emptied file is the original with one member inserted, `"key": []` after
/// a comma, and nothing else anywhere changed. Both the started and the emptied scene parse.
#[test]
fn a_list_a_domain_lacks_is_started_and_emptied() {
    let mut started = 0;
    for (file, text) in shipped() {
        for d in model(&text).domains.iter() {
            let held: Vec<&str> = d.lists.iter().map(|l| l.key.as_str()).collect();
            let lacking: Vec<String> = pantometry_world::templates::ENTRIES
                .iter()
                .filter(|e| e.kind == d.kind && !held.contains(&e.key))
                .map(|e| e.key.to_string())
                .collect();
            // The domain row offers exactly the lists it lacks: no fewer, and none it holds.
            let offered = form(&text, &format!("/domains/{}", d.name)).adds;
            assert_eq!(
                offered, lacking,
                "{file} {}: what the domain row offers",
                d.name
            );
            for key in &lacking {
                let added = editor_core::add_entry(&text, &d.name, key)
                    .unwrap_or_else(|why| panic!("{file} {} {key}: {why}", d.name));
                serde_json::from_str::<Scene>(&added).unwrap_or_else(|e| {
                    panic!("{file} {} {key}: started, it stops parsing: {e}", d.name)
                });
                let emptied = editor_core::remove_entry(&added, &d.name, key, 0)
                    .unwrap_or_else(|why| panic!("{file} {} {key}: {why}", d.name));
                serde_json::from_str::<Scene>(&emptied).unwrap_or_else(|e| {
                    panic!("{file} {} {key}: emptied, it stops parsing: {e}", d.name)
                });
                let same = text
                    .bytes()
                    .zip(emptied.bytes())
                    .take_while(|(a, b)| a == b)
                    .count();
                let grew = emptied.len() - text.len();
                assert_eq!(
                    &emptied[same + grew..],
                    &text[same..],
                    "{file} {} {key}: more than one place changed",
                    d.name
                );
                let inserted = &emptied[same..same + grew];
                assert_eq!(
                    inserted.trim_start_matches(',').trim(),
                    format!("\"{key}\": []"),
                    "{file} {} {key}: inserted {inserted:?}",
                    d.name
                );
                started += 1;
            }
        }
    }
    // 42 on 2026-11-17.
    println!("  {started} lists started and emptied");
    assert!(
        started >= 40,
        "only {started} lists started across the scenes"
    );
}

/// The choices a form offers for the value at the end of `label`.
fn choices(text: &str, row: &str, label: &str) -> Vec<String> {
    let f = form(text, row);
    let v = f
        .values()
        .into_iter()
        .find(|v| v.label == label)
        .unwrap_or_else(|| panic!("{row} has no {label}"))
        .clone();
    match v.value {
        Value::Text { choices, .. } => choices,
        Value::Number { .. } => panic!("{row} {label} is a number"),
    }
}

/// **A menu offers what the format accepts there, and nothing else.**
///
/// The sets below are the format's own, written out from its doc comments rather than from this
/// code: six faces, two holds, three schedules. A room's `"as": "mode"` is a tag whose other value
/// needs other keys, so it is no choice and stays free text. A material keeps the catalogue.
#[test]
fn a_menu_offers_the_formats_own_words() {
    let bracket = one("29-a-designed-bracket-becomes-cells.json");
    assert_eq!(
        choices(&bracket, "/domains/bracket/cooling/0", "face"),
        ["x-min", "x-max", "y-min", "y-max", "z-min", "z-max"]
    );
    assert_eq!(
        choices(&bracket, "/", "schedule"),
        ["one-way", "staggered", "multirate"]
    );
    let material = choices(&bracket, "/domains/bracket/parts/0", "material");
    assert!(
        material.contains(&"aluminium".to_string()) && material.len() > 5,
        "{material:?}"
    );
    assert!(
        choices(&bracket, "/domains/bracket/parts/0", "stl").is_empty(),
        "a path is free text"
    );

    let solder = one("25-what-140-kelvin-does-to-the-solder.json");
    assert_eq!(
        choices(&solder, "/domains/stress/held/0", "as"),
        ["clamp", "roller"]
    );

    let room = one("01-room-mode.json");
    let tag = form(&room, "/domains/room")
        .values()
        .into_iter()
        .find(|v| v.label == "release.as")
        .cloned()
        .expect("the room states how it is released");
    assert_eq!(
        tag.value,
        Value::Text {
            now: "mode".into(),
            choices: Vec::new()
        },
        "pulse is a variant, but not one these keys can be"
    );
}

/// **A form is in the file's order**, with an array's elements on one line.
#[test]
fn a_form_reads_in_the_files_order() {
    let bracket = one("29-a-designed-bracket-becomes-cells.json");
    let f = form(&bracket, "/domains/bracket");
    let labels: Vec<&str> = f
        .rows
        .iter()
        .map(|r| match r {
            Row::One(e) => e.label.as_str(),
            Row::Many { label, .. } => label.as_str(),
        })
        .collect();
    assert_eq!(
        labels,
        ["cells", "cell_mm", "initial_c"],
        "the domain, outside its lists"
    );
    assert!(matches!(&f.rows[0], Row::Many { items, .. } if items.len() == 3));
    assert!(f.adds.contains(&"regions".to_string()) && !f.adds.contains(&"cooling".to_string()));

    let face = form(&bracket, "/domains/bracket/cooling/0");
    let first = face.values().first().map(|v| v.label.clone());
    assert_eq!(
        first.as_deref(),
        Some("face"),
        "a face before the numbers that describe it"
    );
    assert_eq!(face.entry, Some(("cooling".to_string(), 0)));
}

/// **An entry is called what it says**: its name or label, else its strings, else its place.
#[test]
fn an_entry_is_named_by_what_it_says() {
    let bracket = model(&one("29-a-designed-bracket-becomes-cells.json"));
    let lists = &bracket.domains[0].lists;
    let names = |key: &str| lists.iter().find(|l| l.key == key).unwrap().entries.clone();
    assert_eq!(names("cooling"), ["z-min", "z-max"]);
    assert_eq!(names("dissipation"), ["dissipation 1"]);
    assert_eq!(names("parts"), ["parts/l-bracket.stl, aluminium"]);

    let motor = model(&one("11-motor-thermal-network.json"));
    let nodes = motor
        .domains
        .iter()
        .find(|d| d.kind == "network")
        .unwrap()
        .lists
        .iter()
        .find(|l| l.key == "nodes")
        .unwrap();
    assert!(
        nodes.entries.contains(&"winding".to_string()),
        "{:?}",
        nodes.entries
    );
}

/// **A copy is renamed away from its original**, the way a second domain of a kind is: two
/// network nodes called `winding` are refused by the build, and the copy is `winding 2`.
#[test]
fn a_copy_of_a_named_entry_gets_a_name_of_its_own() {
    let motor = one("11-motor-thermal-network.json");
    let d = &model(&motor)
        .domains
        .into_iter()
        .find(|d| d.kind == "network")
        .unwrap();
    let j = d
        .lists
        .iter()
        .find(|l| l.key == "nodes")
        .unwrap()
        .entries
        .iter()
        .position(|e| e == "winding")
        .unwrap();
    let copied = editor_core::duplicate_entry(&motor, &d.name, "nodes", j).expect("copies");
    let after = model(&copied);
    let nodes = &after
        .domains
        .iter()
        .find(|d| d.kind == "network")
        .unwrap()
        .lists
        .iter()
        .find(|l| l.key == "nodes")
        .unwrap()
        .entries;
    assert_eq!(nodes[j + 1], "winding 2", "{nodes:?}");
}

/// **A domain written out of tree gets a tree and a form**, because nothing here names a kind.
/// It gets no add, because the format has no starting entry for a list it has never heard of —
/// and its entries can still be copied.
#[test]
fn a_domain_written_out_of_tree_has_a_tree_and_a_form() {
    let text = r#"{ "title": "t", "duration_s": 1.0, "frames": 2, "domains": [
  { "kind": "made-up", "name": "x", "things": [ { "label": "a", "weight_kg": 1.0 } ] }
] }"#;
    let m = model(text);
    assert_eq!(m.domains[0].lists[0].key, "things");
    assert_eq!(m.domains[0].lists[0].entries, ["a"]);
    let f = form(text, "/domains/x/things/0");
    assert_eq!(
        f.values()
            .iter()
            .map(|v| v.label.as_str())
            .collect::<Vec<_>>(),
        ["label", "weight_kg"]
    );
    assert!(
        form(text, "/domains/x/things").adds.is_empty(),
        "no starting entry for `things`"
    );
    assert!(editor_core::add_entry(text, "x", "things").is_err());
    let copied = editor_core::duplicate_entry(text, "x", "things", 0).expect("copies");
    assert_eq!(model(&copied).domains[0].lists[0].entries, ["a", "a 2"]);
}

/// **A new face avoids the faces every list names, not only its own.** A structure that holds its
/// top is not pressed there: a load on a clamped face moves nothing, and nothing says so.
#[test]
fn a_pressed_face_is_not_a_held_one() {
    let text = r#"{ "title": "t", "duration_s": 1.0, "frames": 2, "domains": [
  { "kind": "structure", "name": "s", "cells": [4, 4, 4], "cell_mm": 1.0, "material": "copper",
    "follows": "b", "reference_c": 20.0,
    "held": [ { "face": "z-max", "as": "clamp" }, { "face": "z-min", "as": "clamp" } ] }
] }"#;
    let pressed = editor_core::add_entry(text, "s", "pressed").expect("adds");
    let value: serde_json::Value = serde_json::from_str(&pressed).expect("JSON");
    let face = value["domains"][0]["pressed"][0]["face"]
        .as_str()
        .expect("a face");
    assert!(
        !["z-max", "z-min"].contains(&face),
        "pressed on a held face: {face}"
    );
}
