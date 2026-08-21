//! **What the outliner and the inspector say, checked without a window.**
//!
//! The tree is the part of a scene browser that can be wrong quietly: a row missing is a domain
//! the reader never learns about, a stale bound is a "frame this" that flies off, and a path that
//! moves between rebuilds is a selection that jumps under the hand that made it. None of those is
//! visible in a screenshot, and all of them are checkable here.
//!
//! Built from a real scene and a real run, through `check` and `run`, because the claim is about
//! what this crate does with what `pantometry-world` produces.

use editor_core::{NodeKind, Tree};

/// A field with an extent, and a domain that has **no geometry at all** — which is the case an
/// outliner built from placed boxes alone would lose entirely.
///
/// This is scene 04 as shipped, rather than a scene written here. Two earlier drafts were
/// physically invalid — one used a field name the format does not have, the other paid joules
/// onto the bus that nothing consumed and the conservation audit stopped the run at t = 0. A test
/// about an outliner should not also be a test of whether its author can write a scene.
const SCENE: &str = r#"{
  "title": "a heater pays joules onto the bus and a bar takes them",
  "schedule": "multirate",
  "duration_s": 4.0,
  "frames": 5,
  "conservation_tolerance": 1e-9,
  "domains": [
    { "kind": "heater", "name": "element", "watts": 2.0, "reserve_j": 6.0 },
    { "kind": "bar", "name": "bar", "length_mm": 20.0, "cells": 61,
      "area_mm2": 100.0, "initial_c": 20.0 }
  ]
}"#;

fn built() -> (editor_core::Checked, viewer_core::Run) {
    let files = editor_core::OnDisk;
    let checked = editor_core::check(SCENE, &files);
    assert!(checked.error.is_none(), "{:?}", checked.error);
    let json = editor_core::run(SCENE, &files).expect("the scene runs");
    let run = viewer_core::Run::from_json(&json).expect("and the viewer reads it");
    (checked, run)
}

/// **Every shape in the scene and the run has a row, including the one with no geometry.**
#[test]
fn nothing_in_the_scene_is_missing_from_the_tree() {
    let (checked, run) = built();
    let tree = editor_core::tree(&checked, Some(&run), 0);
    let paths: Vec<&str> = tree.nodes.iter().map(|n| n.path.as_str()).collect();

    assert_eq!(tree.nodes[0].kind, NodeKind::Scene, "the root is the scene");
    assert!(
        paths.contains(&"/extents/bar"),
        "the placed block: {paths:?}"
    );
    assert!(
        paths.contains(&"/run/bar"),
        "and what the run put in it: {paths:?}"
    );

    // The heater has no extent and no panel. It reports scalars, and those are the only trace of
    // it in the run — so a tree that walks geometry loses a domain, which is the whole reason
    // readings are in here.
    assert!(
        paths.iter().any(|p| p.starts_with("/readings/element")),
        "the heater has no geometry and must still appear: {paths:?}"
    );

    // Parents come before children, which is what lets a shell draw this with one pass and an
    // indent counter rather than a recursive walk.
    for (i, n) in tree.nodes.iter().enumerate() {
        if let Some(parent) = n.parent {
            assert!(parent < i, "{} comes before its parent", n.path);
            assert_eq!(
                tree.nodes[parent].depth + 1,
                n.depth,
                "{} is indented {} under a parent at {}",
                n.path,
                n.depth,
                tree.nodes[parent].depth
            );
        }
    }
}

/// **A path is the same path after a rebuild.**
///
/// The tree is rebuilt on every check and every streamed frame. A selection is stored as a path,
/// so if a path moved between rebuilds the selection would follow it to a different object — the
/// exact failure that makes a row index the wrong handle.
#[test]
fn a_path_survives_a_rebuild_and_a_different_frame() {
    let (checked, run) = built();
    let a = editor_core::tree(&checked, Some(&run), 0);
    let b = editor_core::tree(&checked, Some(&run), run.frames.len() - 1);
    let pa: Vec<&String> = a.nodes.iter().map(|n| &n.path).collect();
    let pb: Vec<&String> = b.nodes.iter().map(|n| &n.path).collect();
    assert_eq!(
        pa, pb,
        "the same scene at a different frame is the same tree"
    );

    // And a tree built before the run exists is a prefix of the paths, not a different naming.
    let bare = editor_core::tree(&checked, None, 0);
    for n in &bare.nodes {
        assert!(
            pa.contains(&&n.path),
            "{} exists before the run and vanishes after it",
            n.path
        );
    }
}

/// **Hiding a domain hides what is under it, and nothing that merely starts with its name.**
///
/// `starts_with` alone says `/room` contains `/roomier`. The separator is the whole test.
#[test]
fn containment_is_by_path_and_not_by_prefix() {
    assert!(Tree::contains("/run", "/run/bar"));
    assert!(Tree::contains("/run/bar", "/run/bar"));
    assert!(Tree::contains(
        "/readings/element",
        "/readings/element/reserve"
    ));
    assert!(!Tree::contains("/run/bar", "/run/barometer"));
    assert!(!Tree::contains("/run/bar", "/run"));
    assert!(!Tree::contains("/extents/bar", "/run/bar"));
}

/// **The inspector's numbers are the run's, not the frame's, where it says so.**
///
/// A range taken over one frame makes a decaying quantity look constant — the same reason the
/// colour scale spans the run. Both are shown, and both are labelled, because they answer
/// different questions and a picture drawn on one is often read as the other.
#[test]
fn the_inspector_separates_the_runs_range_from_this_frames() {
    let (checked, run) = built();
    let last = run.frames.len() - 1;
    let first = editor_core::tree(&checked, Some(&run), 0);
    let later = editor_core::tree(&checked, Some(&run), last);

    let get = |t: &Tree, key: &str| -> String {
        t.find("/run/bar")
            .expect("the block has a row")
            .detail
            .iter()
            .find(|(k, _)| k == key)
            .unwrap_or_else(|| panic!("no {key} row"))
            .1
            .clone()
    };

    assert_eq!(
        get(&first, "range (run)"),
        get(&later, "range (run)"),
        "the run's range must not depend on which frame is open"
    );
    // The heater is paying joules in, so the bar is not the same at the end as at the start.
    assert_ne!(
        get(&first, "this frame"),
        get(&later, "this frame"),
        "2 W into a bar for four seconds has to move something"
    );

    // And the grid is the grid, stated in cells and in metres.
    assert_eq!(get(&first, "grid"), "61 x 1 x 1");
    let size = get(&first, "size");
    assert!(
        size.contains("20 mm"),
        "the bar is 20 mm long and the inspector says {size}"
    );
}

/// **A row that can be framed has a box, and a row that cannot does not claim one.**
///
/// "Frame this" on a scalar would fly the camera to the origin.
#[test]
fn only_the_rows_with_geometry_offer_a_box() {
    let (checked, run) = built();
    let tree = editor_core::tree(&checked, Some(&run), 0);
    for n in &tree.nodes {
        match n.kind {
            NodeKind::Reading => assert!(n.bounds.is_none(), "{} claims a box", n.path),
            NodeKind::Placed | NodeKind::Field => {
                let b = n.bounds.unwrap_or_else(|| panic!("{} has no box", n.path));
                // At least one axis, not all three: a bar is a line and a plane cut through a
                // solid is a plane, and both are legitimately flat in the axes they collapse.
                assert!(
                    (0..3).any(|a| b[a + 3] > b[a]),
                    "{} has a box with no extent on any axis: {b:?}",
                    n.path
                );
            }
            _ => {}
        }
    }
}

/// A length is written the way an engineer says it, and 500 does not become 5.
///
/// The trailing-zero trim that produced that shipped in the HTML report and captioned a 500 um
/// channel as a 5 um one, which is the kind of wrong that looks like a plausible number.
#[test]
fn a_length_keeps_its_zeros_when_they_are_not_after_a_point() {
    assert_eq!(editor_core::metres(0.5), "500 mm");
    assert_eq!(editor_core::metres(0.0005), "500 um");
    assert_eq!(editor_core::metres(0.02), "20 mm");
    assert_eq!(editor_core::metres(1.2), "1.2 m");
    assert_eq!(editor_core::metres(0.0), "0 m");
    assert!(editor_core::metres(3.11e11).contains('e'));
}
