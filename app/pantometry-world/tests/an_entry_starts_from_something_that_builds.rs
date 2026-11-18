//! The editor's "add an entry" — a cooled face, a heat source, a dose — against the format.
//!
//! [`templates::ENTRIES`] is the starting point for one entry of each list a domain can hold, and
//! like [`templates::TEMPLATES`] it is a table kept by hand. So it is held against two things that
//! are not: the lists the thirty-four shipped scenes use, and the format's own build.
//!
//! # Not under `wasm32`
//!
//! The first test reads the scenes off a disk, and a `wasm32` target has none.

#![cfg(not(target_family = "wasm"))]

use pantometry_world::templates::{self, Entry, ENTRIES};
use pantometry_world::{Beside, Scene, World};

/// Where the shipped scenes and their files are.
fn scenes() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("scenes")
}

/// Every `(kind, key)` whose value is a non-empty list of objects, across the shipped scenes.
fn lists_the_scenes_use() -> std::collections::BTreeSet<(String, String)> {
    let mut out = std::collections::BTreeSet::new();
    let mut read = 0;
    for entry in std::fs::read_dir(scenes()).expect("the scenes directory reads") {
        let path = entry.expect("a directory entry").path();
        if path.extension().is_none_or(|x| x != "json") {
            continue;
        }
        read += 1;
        let text = std::fs::read_to_string(&path).expect("a scene reads");
        let value: serde_json::Value =
            serde_json::from_str(&text).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        for domain in value["domains"].as_array().into_iter().flatten() {
            let (Some(kind), Some(object)) = (domain["kind"].as_str(), domain.as_object()) else {
                continue;
            };
            for (key, member) in object {
                if member
                    .as_array()
                    .is_some_and(|l| !l.is_empty() && l.iter().all(|x| x.is_object()))
                {
                    out.insert((kind.to_string(), key.clone()));
                }
            }
        }
    }
    // An empty directory would make the comparison below vacuously true.
    assert!(
        read >= 28,
        "only {read} scenes found in {}",
        scenes().display()
    );
    out
}

/// **Every list a shipped scene uses has a starting entry.**
///
/// Fifteen, measured when this was written. The other direction — an entry for a list no scene
/// uses — is allowed and is two of them, `infusions` and `pressed`; what keeps those honest is
/// the build below, since a key the format does not take is refused there.
#[test]
fn every_list_a_scene_uses_has_an_entry() {
    let used = lists_the_scenes_use();
    assert!(
        used.len() >= 15,
        "only {} lists found across the scenes — this compares against them",
        used.len()
    );
    for (kind, key) in &used {
        assert!(
            templates::entry(kind, key).is_some(),
            "a shipped scene's {kind} has a `{key}` list and the editor has no entry to start one \
             from. Add it to `templates::ENTRIES`"
        );
    }
}

/// A scene holding `kind`'s template, in which that template builds on its own.
///
/// `structure` cannot stand alone — it names the block it `follows`, and that block's grid has to
/// be its own — so it gets one, which is what [`templates::scene`] leaves for a person to do.
fn scene_for(kind: &str) -> serde_json::Value {
    if kind != "structure" {
        return serde_json::from_str(&templates::scene(&[kind])).expect("a template scene parses");
    }
    let mut value: serde_json::Value =
        serde_json::from_str(&templates::scene(&["block", "structure"])).expect("parses");
    for d in value["domains"].as_array_mut().expect("domains") {
        if d["kind"] == "structure" {
            d["follows"] = "block".into();
        }
        if d["kind"] == "block" {
            d["cells"] = serde_json::json!([8, 8, 8]);
            d["cell_mm"] = serde_json::json!(1.5);
        }
    }
    value
}

/// `value` with `entry` appended to its `kind` domain's list, which is created when absent.
fn with(mut value: serde_json::Value, entry: &Entry) -> serde_json::Value {
    let d = value["domains"]
        .as_array_mut()
        .expect("domains")
        .iter_mut()
        .find(|d| d["kind"] == entry.kind)
        .expect("the kind is in its own scene");
    let item: serde_json::Value = serde_json::from_str(entry.json).unwrap_or_else(|e| {
        panic!(
            "{} {}: the entry does not parse: {e}",
            entry.kind, entry.key
        )
    });
    match d.get_mut(entry.key).and_then(|l| l.as_array_mut()) {
        Some(list) => list.push(item),
        None => d[entry.key] = serde_json::Value::Array(vec![item]),
    }
    value
}

fn build(value: &serde_json::Value) -> Result<(), String> {
    let text = serde_json::to_string(value).expect("serialises");
    let scene: Scene = serde_json::from_str(&text).map_err(|e| format!("parse: {e}"))?;
    World::build_with(scene, &Beside(scenes())).map(|_| ())
}

/// **Each entry builds in its kind's template, and three name what they need instead.**
///
/// Measured on 2026-11-17: fourteen build as they are. `parts` names `parts/rod.stl`, 16 x 16 x
/// 60 mm, which an 11 mm block cannot hold, and is refused as cut off rather than cropped. A
/// `depots` entry with no oral dose through it, and an `oral` dose with no depot, are each refused
/// for lacking the other. The set of three is pinned, so a fourth is a failure here rather than a
/// button in the editor that only ever produces an error.
///
/// Each of the three then builds once its need is met: the block grown to the rod's size, and
/// the depot and the dose added together.
#[test]
fn each_entry_builds_in_its_template_or_names_what_it_needs() {
    let mut named = Vec::new();
    for e in &ENTRIES {
        match build(&with(scene_for(e.kind), e)) {
            Ok(()) => {}
            Err(why) => {
                let expected = match e.key {
                    "parts" => "cut off",
                    // The refusal's own words, not the key's name: a serde error that listed
                    // the fields it expected would name `oral` too.
                    "depots" => "has no dose in `oral`",
                    "oral" => "no depot",
                    _ => panic!(
                        "{} {}: the entry does not build in its template: {why}",
                        e.kind, e.key
                    ),
                };
                assert!(
                    why.contains(expected),
                    "{} {}: refused for another reason than the one pinned: {why}",
                    e.kind,
                    e.key
                );
                named.push(e.key);
            }
        }
    }
    assert_eq!(
        named,
        vec!["parts", "depots", "oral"],
        "the set that needs something changed"
    );

    let mut big = scene_for("block");
    big["domains"][0]["cells"] = serde_json::json!([16, 16, 60]);
    let parts = templates::entry("block", "parts").expect("a parts entry");
    build(&with(big, parts)).unwrap_or_else(|e| panic!("a block that holds the rod: {e}"));

    let depot = templates::entry("compartments", "depots").expect("a depots entry");
    let oral = templates::entry("compartments", "oral").expect("an oral entry");
    build(&with(with(scene_for("compartments"), depot), oral))
        .unwrap_or_else(|e| panic!("a depot and a dose through it: {e}"));
}

/// **Sorted, unique, and one line each**, because a menu is built from this and an entry is
/// spliced into a list at whatever indentation the file has.
#[test]
fn the_entries_are_sorted_unique_and_one_line() {
    let pairs: Vec<(&str, &str)> = ENTRIES.iter().map(|e| (e.kind, e.key)).collect();
    let mut sorted = pairs.clone();
    sorted.sort();
    sorted.dedup();
    assert_eq!(
        pairs, sorted,
        "ENTRIES should be sorted by kind and key, with no repeats"
    );
    for e in &ENTRIES {
        assert!(
            !e.json.contains('\n'),
            "{} {}: an entry is one line",
            e.kind,
            e.key
        );
        let value: serde_json::Value = serde_json::from_str(e.json).expect("parses");
        assert!(
            value.is_object(),
            "{} {}: an entry is an object",
            e.kind,
            e.key
        );
    }
}

/// A domain's grid, `cells` and `cell_mm`, when it states both.
fn grid_of(d: &serde_json::Value) -> Option<([usize; 3], f64)> {
    let cells: Vec<usize> = d
        .get("cells")?
        .as_array()?
        .iter()
        .map(|v| v.as_u64().map(|v| v as usize))
        .collect::<Option<_>>()?;
    Some((
        <[usize; 3]>::try_from(cells).ok()?,
        d.get("cell_mm")?.as_f64()?,
    ))
}

/// The faces domain `d`'s lists already name, every one of them.
fn faces_in(d: &serde_json::Value) -> Vec<&str> {
    d.as_object()
        .into_iter()
        .flat_map(|o| o.values())
        .filter_map(|m| m.as_array())
        .flatten()
        .filter_map(|e| e.get("face").and_then(|f| f.as_str()))
        .collect()
}

/// What `entry_for` is told about domain `d`.
fn fit_of<'a>(d: &'a serde_json::Value, faces: &'a [&'a str]) -> templates::Fit<'a> {
    templates::Fit {
        grid: grid_of(d),
        faces,
        material: d.get("material").and_then(|m| m.as_str()),
    }
}

/// Every shipped scene, parsed, with its file name.
fn shipped() -> Vec<(String, serde_json::Value)> {
    let mut out = Vec::new();
    for entry in std::fs::read_dir(scenes()).expect("the scenes directory reads") {
        let path = entry.expect("a directory entry").path();
        if path.extension().is_some_and(|x| x == "json") {
            let text = std::fs::read_to_string(&path).expect("a scene reads");
            let value = serde_json::from_str(&text).expect("a scene is JSON");
            out.push((
                path.file_name().unwrap().to_string_lossy().to_string(),
                value,
            ));
        }
    }
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

/// **What an entry says about where in a grid is that grid's own.**
///
/// For each kind's template, [`templates::entry_for`] gives exactly what [`ENTRIES`] says, so the
/// table and the fitting cannot drift apart. For every block in the shipped scenes, a cooled face
/// added to it states the area of the whole face it covers — the template's `1.21` would have been
/// a silent wrong answer on all but an 11 mm block — and a box or a plane added to it lies inside
/// its grid. And the cooled face builds, in every one of them: it is the entry a person adds first
/// and the one whose mistake the build could not have caught.
#[test]
fn an_entry_fits_the_grid_it_is_added_to() {
    for e in &ENTRIES {
        let scene = scene_for(e.kind);
        let d = scene["domains"]
            .as_array()
            .expect("domains")
            .iter()
            .find(|d| d["kind"] == e.kind)
            .expect("the kind is in its own scene");
        let faces = faces_in(d);
        assert_eq!(
            templates::entry_for(e.kind, e.key, &fit_of(d, &faces)).as_deref(),
            Some(e.json),
            "{} {}: fitted to its own template, an entry should be the table's",
            e.kind,
            e.key
        );
    }

    let mut blocks = 0;
    for (file, value) in shipped() {
        let domains = value["domains"].as_array().cloned().unwrap_or_default();
        for (i, d) in domains.iter().enumerate() {
            if d["kind"] != "block" {
                continue;
            }
            let n: Vec<usize> = d["cells"]
                .as_array()
                .expect("a block states its cells")
                .iter()
                .map(|v| v.as_u64().expect("a count") as usize)
                .collect();
            let mm = d["cell_mm"].as_f64().expect("a block states its cell");
            let faces = faces_in(d);
            let fitted = |key: &str| -> serde_json::Value {
                let text =
                    templates::entry_for("block", key, &fit_of(d, &faces)).expect("an entry");
                serde_json::from_str(&text).unwrap_or_else(|e| panic!("{file} {key}: {e}"))
            };

            let cooled = fitted("cooling");
            let side = cooled["face"].as_str().expect("a face").to_string();
            let taken = &faces;
            assert!(
                !taken.contains(&side.as_str()),
                "{file} {}: a new cooled face on {side}, which {taken:?} already holds",
                d["name"]
            );
            // The two counts that lie in the face, worked out here rather than by the code under
            // test: an x face spans y and z, a y face x and z, a z face x and y.
            let (p, q) = match &side[..1] {
                "x" => (n[1], n[2]),
                "y" => (n[0], n[2]),
                _ => (n[0], n[1]),
            };
            let area = cooled["area_cm2"].as_f64().expect("an area");
            let face = p as f64 * q as f64 * mm * mm / 100.0;
            assert!(
                (area - face).abs() <= 1e-12 * face,
                "{file} {}: a cooled {side} face states {area} cm² of a {face} cm² face",
                d["name"]
            );
            // A region of the block's own material would change nothing, and say nothing.
            let region = fitted("regions");
            assert_ne!(
                region["material"], d["material"],
                "{file} {}: a region of what the block already is",
                d["name"]
            );
            for key in ["dissipation", "regions"] {
                let entry = fitted(key);
                for a in 0..3 {
                    let (from, to) = (
                        entry["from"][a].as_u64().unwrap(),
                        entry["to"][a].as_u64().unwrap(),
                    );
                    assert!(
                        from < to && to as usize <= n[a],
                        "{file} {key}: cells {from}..{to} on an axis of {}",
                        n[a]
                    );
                }
            }
            let at = fitted("contact")["at"].as_u64().unwrap() as usize;
            assert!(
                at >= 1 && at < n[2].max(2),
                "{file}: a contact at {at} of {}",
                n[2]
            );

            let mut cooled = value.clone();
            let list = &mut cooled["domains"][i]["cooling"];
            if !list.is_array() {
                *list = serde_json::json!([]);
            }
            list.as_array_mut().unwrap().push(fitted("cooling"));
            build(&cooled).unwrap_or_else(|why| panic!("{file}: a cooled face added: {why}"));
            blocks += 1;
        }
    }
    // Every block that ships, which is more than ten; a count near zero would mean the walk found
    // no blocks and checked nothing.
    assert!(blocks >= 10, "only {blocks} blocks found across the scenes");
}

/// **Every face's area, on a grid with three different counts.** The shipped blocks reach only a
/// z face and an x face, and both that reach an x face are square in x and y — so an x face given
/// a y face's area, or any area at all for a y face, passed every test above. Here each face is
/// taken in turn, so the next entry lands on the next one, and each area is the two counts that
/// lie in that face, worked out here: 3 x 5 x 7 cells of 1 mm.
#[test]
fn every_face_states_its_own_area() {
    let n = [3usize, 5, 7];
    let mut taken: Vec<String> = Vec::new();
    for _ in 0..6 {
        let held: Vec<&str> = taken.iter().map(String::as_str).collect();
        let fit = templates::Fit {
            grid: Some((n, 1.0)),
            faces: &held,
            material: None,
        };
        let text = templates::entry_for("block", "cooling", &fit).expect("an entry");
        let entry: serde_json::Value = serde_json::from_str(&text).expect("parses");
        let face = entry["face"].as_str().expect("a face").to_string();
        assert!(!taken.contains(&face), "{face} again, after {taken:?}");
        let (p, q) = match &face[..1] {
            "x" => (n[1], n[2]),
            "y" => (n[0], n[2]),
            _ => (n[0], n[1]),
        };
        let area = entry["area_cm2"].as_f64().expect("an area");
        let want = (p * q) as f64 / 100.0;
        assert!(
            (area - want).abs() <= 1e-12,
            "{face}: {area} cm², not {want}"
        );
        taken.push(face);
    }
    taken.sort();
    assert_eq!(
        taken,
        ["x-max", "x-min", "y-max", "y-min", "z-max", "z-min"]
    );
}

/// **A pressed face is not a held one.** The template holds three faces with rollers; a
/// structure that also holds the top is not pressed there, since a load on a clamped face moves
/// nothing and nothing says so.
#[test]
fn a_new_face_avoids_every_list_not_only_its_own() {
    let held = ["x-min", "y-min", "z-min", "z-max"];
    let fit = templates::Fit {
        grid: None,
        faces: &held,
        material: None,
    };
    let pressed = templates::entry_for("structure", "pressed", &fit).expect("an entry");
    assert!(pressed.contains("\"face\": \"x-max\""), "{pressed}");
}

/// **Every list that says where on a grid is one a click can set, and says it the way its type
/// does.**
///
/// [`templates::SITES`] is a table kept by hand, like [`ENTRIES`], so it is held against the
/// format three ways. Each entry naming a `face` is a face site and each naming three-cell `from`
/// and `to` is a cells site, so a list that says where cannot be missing from the table; every
/// site has an entry, so the table names no list the format lacks; and a face site takes a patch
/// exactly when its type does — `from` and `to` added to a cooled face builds, and added to a held
/// or pressed one is refused as an unknown field, which is what a click would otherwise write.
#[test]
fn every_list_that_says_where_is_one_a_click_can_set() {
    use templates::Placing;
    let mut faces = 0;
    let mut boxes = 0;
    for entry in &ENTRIES {
        let item: serde_json::Value = serde_json::from_str(entry.json).expect("parses");
        let three = |k: &str| item[k].as_array().is_some_and(|a| a.len() == 3);
        let says = if item.get("face").is_some() {
            faces += 1;
            Some("face")
        } else if three("from") && three("to") {
            boxes += 1;
            Some("cells")
        } else {
            None
        };
        let placing = templates::placing(entry.kind, entry.key);
        match (says, placing) {
            (Some("face"), Some(Placing::Face | Placing::FacePatch))
            | (Some("cells"), Some(Placing::Cells))
            | (None, None) => {}
            _ => panic!(
                "{} {}: the entry says {says:?} and the table says {placing:?}",
                entry.kind, entry.key
            ),
        }
        if !matches!(placing, Some(Placing::Face | Placing::FacePatch)) {
            continue;
        }
        let mut patched = item.clone();
        patched["from"] = serde_json::json!([0, 0]);
        patched["to"] = serde_json::json!([1, 1]);
        let json = serde_json::to_string(&patched).expect("serialises");
        let with_patch = Entry {
            kind: entry.kind,
            key: entry.key,
            json: Box::leak(json.into_boxed_str()),
        };
        let result = build(&with(scene_for(entry.kind), &with_patch));
        match placing {
            Some(Placing::FacePatch) => {
                result.unwrap_or_else(|e| panic!("{} {} with a patch: {e}", entry.kind, entry.key));
            }
            _ => {
                let e = result.expect_err("a whole-face entry took a patch");
                assert!(
                    e.contains("unknown field"),
                    "{} {} with a patch was refused, but not for the patch: {e}",
                    entry.kind,
                    entry.key
                );
            }
        }
    }
    for site in &templates::SITES {
        assert!(
            templates::entry(site.kind, site.key).is_some(),
            "{} {} is a site with no entry in the format",
            site.kind,
            site.key
        );
    }
    // Three face lists and three box lists when this was written; none would make it vacuous.
    assert_eq!((faces, boxes), (3, 3), "face and cells entries");
    // The words a click writes are the format's, and each face is across the axis it says.
    for (n, face) in pantometry_world::FaceSpec::ALL.into_iter().enumerate() {
        assert_eq!(
            serde_json::to_value(face).expect("serialises"),
            serde_json::Value::from(face.spelled()),
            "{face:?}"
        );
        let (axis, high) = face.normal();
        assert_eq!((axis, high), (n / 2, n % 2 == 1), "{face:?}");
        assert_eq!(pantometry_world::FaceSpec::across(axis, high), face);
        assert!(
            !face.axes().contains(&axis),
            "{face:?} lies across its own axis"
        );
    }
}
