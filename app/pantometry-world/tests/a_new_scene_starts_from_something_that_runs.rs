//! The editor's "New project" chooser, against the scenes that prove its numbers.
//!
//! [`templates::scene`] turns a set of kinds into a scene to open. Two things in it are numbers
//! rather than text — `duration_s` and `frames` — and a starting point whose duration nothing has
//! ever run is a starting point that opens on a refusal. So they are not chosen here: each is
//! taken from a shipped scene that uses that kind, and this holds them against those scenes.
//!
//! That is the shape `every_domain_has_a_template` established, for the same reason: a table
//! maintained by hand is a table that is silently wrong between the change and the noticing, and
//! the thirty scenes are the set that is maintained by something else.
//!
//! # Not under `wasm32`
//!
//! Every test here reads the scenes off a disk, and a `wasm32` target has none.

#![cfg(not(target_family = "wasm"))]

use pantometry_world::templates::{self, TEMPLATES};
use pantometry_world::{OnDisk, Scene, World};

/// Every `(duration_s, frames)` pair each kind appears with, across the shipped scenes.
fn schedules_per_kind() -> std::collections::BTreeMap<String, Vec<(f64, u64)>> {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("scenes");
    let mut out: std::collections::BTreeMap<String, Vec<(f64, u64)>> = Default::default();
    let mut scenes = 0;
    for entry in std::fs::read_dir(&dir).expect("the scenes directory reads") {
        let path = entry.expect("a directory entry").path();
        if path.extension().is_none_or(|x| x != "json") {
            continue;
        }
        scenes += 1;
        let text = std::fs::read_to_string(&path).expect("a scene reads");
        let value: serde_json::Value =
            serde_json::from_str(&text).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let (Some(secs), Some(frames)) = (value["duration_s"].as_f64(), value["frames"].as_u64())
        else {
            continue;
        };
        for domain in value["domains"].as_array().into_iter().flatten() {
            if let Some(kind) = domain["kind"].as_str() {
                out.entry(kind.to_string())
                    .or_default()
                    .push((secs, frames));
            }
        }
    }
    // A directory that read as empty would make every comparison below vacuously true, which is
    // the shape of a check that turned itself off.
    assert!(
        scenes >= 28,
        "only {scenes} scenes found in {} — this test compares against them",
        dir.display()
    );
    out
}

/// **Every number a new scene starts from is a number some shipped scene runs with.**
#[test]
fn the_schedule_a_template_suggests_is_one_a_scene_uses() {
    let shipped = schedules_per_kind();
    for t in TEMPLATES {
        let used = shipped
            .get(t.kind)
            .unwrap_or_else(|| panic!("{}: no shipped scene uses this kind", t.kind));
        assert!(
            used.contains(&(t.duration_s, t.frames as u64)),
            "{}: the template suggests {} s in {} frames, which no scene that uses it runs. \
             The scenes run it at {used:?}",
            t.kind,
            t.duration_s,
            t.frames
        );
    }
}

/// **A scene made of one kind parses, and builds unless the format will not let it stand alone.**
///
/// `structure` names the block it `follows` and is refused by the build; `beam` names what it
/// shines `onto`, builds, and is stopped by the conservation audit at the first step. Both are
/// left as the templates write them — see [`templates::scene`] — so this pins that the set of
/// kinds a new scene cannot open cleanly is exactly those two.
#[test]
fn a_new_scene_of_one_kind_builds_or_names_what_it_needs() {
    let mut refused = Vec::new();
    for t in TEMPLATES {
        let text = templates::scene(&[t.kind]);
        let scene: Scene = serde_json::from_str(&text)
            .unwrap_or_else(|e| panic!("{}: a new scene does not parse: {e}\n{text}", t.kind));
        assert_eq!(
            scene.domains.len(),
            1,
            "{}: one kind should make one domain",
            t.kind
        );
        match World::build_with(scene, &OnDisk) {
            Ok(_) => {}
            Err(why) if why.contains("follows") => refused.push(t.kind),
            Err(why) => panic!(
                "{}: a new scene of one does not build: {why}\n{text}",
                t.kind
            ),
        }
    }
    assert_eq!(
        refused,
        vec!["structure"],
        "the set of kinds a new scene cannot open has changed"
    );
}

/// **A scene made of several kinds parses and builds**, at the shortest of their durations.
///
/// Three that share a timescale and are already used together in the shipped scenes.
#[test]
fn a_new_scene_of_several_kinds_builds_at_the_shortest_duration() {
    let text = templates::scene(&["bar", "heater", "lump"]);
    let scene: Scene = serde_json::from_str(&text)
        .unwrap_or_else(|e| panic!("three kinds do not parse: {e}\n{text}"));
    assert_eq!(scene.domains.len(), 3, "three kinds, three domains\n{text}");
    // `lump` is the shortest of the three at 1.0 s; `bar` and `heater` are 4.0.
    assert!(
        text.contains("\"duration_s\": 1"),
        "the shortest duration is not the one written:\n{text}"
    );
    World::build_with(scene, &OnDisk).unwrap_or_else(|e| panic!("{e}\n{text}"));
}

/// **Exactly two kinds name a partner**, and they are the two the build and the audit refuse.
///
/// `needs_a_partner` reads the template's own text for `onto` and `follows`. Deriving it means a
/// third such key would be missed, so the set it produces is pinned against the set the test above
/// measures from the format itself.
#[test]
fn exactly_two_kinds_name_a_partner() {
    let named: Vec<&str> = TEMPLATES
        .iter()
        .filter(|t| t.needs_a_partner())
        .map(|t| t.kind)
        .collect();
    assert_eq!(named, vec!["beam", "structure"]);
}

/// **The span is a ratio, and it is 1 for one kind and enormous for the wrong two.**
///
/// The number the chooser says out loud. `atoms` settles in picoseconds and a thermal `network` in
/// half an hour: there is no `duration_s` at which both are a simulation, and the format cannot
/// refuse the scene because it is well formed.
#[test]
fn the_timescale_span_says_when_a_set_cannot_run_together() {
    assert_eq!(templates::timescale_span(&["room"]), 1.0);
    assert!(
        (templates::timescale_span(&["bar", "heater"]) - 1.0).abs() < 1e-12,
        "two kinds a scene already runs together are not one apart"
    );
    let far = templates::timescale_span(&["atoms", "network"]);
    assert!(
        far > 1e12,
        "picoseconds against half an hour came out as {far}"
    );
}
