//! The editor's "add a domain" list, against the format it claims to cover.
//!
//! Rust cannot enumerate an enum's variants, so nothing can prove [`TEMPLATES`] is complete by
//! construction. What can be checked is that it agrees with a set maintained somewhere else, and
//! the shipped scenes are that set: every kind the format defines appears in at least one of
//! them. Comparing the two **in both directions** is what makes the pair load-bearing — one half
//! catches a domain added without a template, the other catches a template for something that is
//! not there.
//!
//! That is `counts_in_prose`'s shape, and it is here for the reason that test exists: a list that
//! only ever gets longer by hand is a list that is silently wrong between the addition and the
//! noticing.
//!
//! # Not under `wasm32`
//!
//! Every test here reads the scenes off a disk, and a `wasm32` target has none. `editor-core`
//! compiles to wasm, which is exactly why the templates live in this crate rather than beside the
//! code that uses them.

#![cfg(not(target_family = "wasm"))]

use pantometry_world::templates::TEMPLATES;
use pantometry_world::{OnDisk, Scene, World};
use std::collections::BTreeSet;

/// The `kind` of every domain in every shipped scene.
fn kinds_in_scenes() -> BTreeSet<String> {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("scenes");
    let mut kinds = BTreeSet::new();
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
        for domain in value["domains"].as_array().into_iter().flatten() {
            if let Some(kind) = domain["kind"].as_str() {
                kinds.insert(kind.to_string());
            }
        }
    }
    // A directory that read as empty would make every comparison below vacuously true, which is
    // the shape of a check that turned itself off.
    assert!(
        scenes >= 28,
        "only {scenes} scenes found in {} — this test compares against them and cannot do it \
         from an empty directory",
        dir.display()
    );
    kinds
}

/// **Every kind the scenes use has a template, and every template is a kind the scenes use.**
///
/// Both directions, because each catches a different mistake and neither catches the other's.
#[test]
fn the_templates_and_the_scenes_name_the_same_domains() {
    let in_scenes = kinds_in_scenes();
    let in_templates: BTreeSet<String> = TEMPLATES.iter().map(|(k, _)| k.to_string()).collect();

    let missing: Vec<&String> = in_scenes.difference(&in_templates).collect();
    assert!(
        missing.is_empty(),
        "a scene uses {missing:?} and the editor has no template for it — add one to \
         pantometry_world::templates"
    );

    let extra: Vec<&String> = in_templates.difference(&in_scenes).collect();
    assert!(
        extra.is_empty(),
        "there is a template for {extra:?} and no shipped scene uses it — either the kind is \
         gone, or it needs a scene"
    );

    assert_eq!(
        TEMPLATES.len(),
        in_templates.len(),
        "two templates share a kind, so one of them can never be reached"
    );
}

/// **Every template is a domain the format accepts**, checked by parsing it as one.
///
/// The templates are text, so nothing but this says they are still valid after a field is renamed
/// or a default removed. Parsing is the cheap half; the next test does the expensive half.
#[test]
fn every_template_parses_as_the_domain_it_claims_to_be() {
    for (kind, text) in TEMPLATES {
        let value: serde_json::Value = serde_json::from_str(text)
            .unwrap_or_else(|e| panic!("{kind}: the template is not JSON: {e}\n{text}"));
        assert_eq!(
            value["kind"].as_str(),
            Some(kind),
            "{kind}: the template's own `kind` disagrees with the key it is filed under"
        );
        assert_eq!(
            value["name"].as_str(),
            Some(kind),
            "{kind}: a template's name is its kind, which the inserter then makes unique"
        );
    }
}

/// **Every template builds alone, except the one the format will not let stand alone.**
///
/// Two templates name another domain — `beam` states what it shines `onto`, `structure` states
/// the block it `follows` — and both keys are required, so no template for them can be complete
/// on its own. Whatever name is written here is one the receiving scene probably does not have.
///
/// **The two are caught at different times, which is measured here rather than assumed.**
/// `structure`'s dangling `follows` is refused by `World::build`. `beam`'s dangling `onto` is
/// **not**: it builds, and the kernel's conservation audit stops the run at `t = 0` with
/// *"published but not consumed"*, naming the absent domain. Nothing is silent either way, which
/// is what matters — but the second message is about a face and a relative change of `1.0`, where
/// the first is about a name, so the two are not equally easy to act on.
///
/// This pins the set, so a third reference arriving unnoticed is a test failure rather than a
/// surprise in the editor's menu.
#[test]
fn every_template_builds_except_the_one_that_cannot() {
    let mut refused_at_build = Vec::new();
    for (kind, text) in TEMPLATES {
        let json = format!(
            r#"{{ "title": "one {kind}", "duration_s": 0.001, "frames": 1,
  "domains": [ {text} ] }}"#
        );
        let scene: Scene = serde_json::from_str(&json)
            .unwrap_or_else(|e| panic!("{kind}: the template does not fit a scene: {e}\n{json}"));
        match World::build_with(scene, &OnDisk) {
            Ok(_) => {}
            Err(why) if why.contains("follows") => refused_at_build.push(kind),
            Err(why) => panic!("{kind}: the template does not build: {why}"),
        }
    }
    assert_eq!(
        refused_at_build,
        vec!["structure"],
        "the set of templates the build refuses has changed"
    );
}

/// **A beam aimed at nothing is refused, by the audit rather than by the build.**
///
/// The other half of the sentence above, checked rather than described. It matters for the
/// editor: adding a `beam` from the menu produces a scene that *checks clean* and then fails the
/// moment it is run, so the inspector shows no error and the Run button does.
#[test]
fn a_beam_aimed_at_a_domain_that_is_not_there_is_stopped_at_the_first_step() {
    let beam = TEMPLATES
        .iter()
        .find(|(k, _)| *k == "beam")
        .expect("there is a beam template")
        .1;
    let json = format!(
        r#"{{ "title": "a beam pointed at nothing", "duration_s": 0.05, "frames": 3,
  "domains": [
    {{ "kind": "bar", "name": "bar", "length_mm": 20.0, "cells": 61, "area_mm2": 100.0,
      "initial_c": 20.0 }},
    {beam} ] }}"#
    );
    let scene: Scene = serde_json::from_str(&json).expect("the scene parses");
    let mut world = World::build_with(scene, &OnDisk).expect("and it builds, which is the point");
    let why = world
        .run()
        .expect_err("a beam publishing into nowhere cannot conserve");
    assert!(
        why.to_string().contains("not consumed"),
        "the audit should say the energy went nowhere: {why}"
    );
}
