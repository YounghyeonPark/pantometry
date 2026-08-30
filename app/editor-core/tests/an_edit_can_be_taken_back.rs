//! **Undo, walked in both directions and past both ends.**
//!
//! The editor gained six ways to change a scene before it gained any way to change one back:
//! a value, a string, a domain added, a domain removed, a placement, and a drag on an axis. Two
//! of those destroy something — a removed domain takes its whole block of JSON with it, and a
//! drag rewrites a number nobody is necessarily looking at.
//!
//! The logic is a list and an index, which is the kind of thing that is right in one direction
//! and wrong in the other until a test walks it both. So every case here goes forward and back,
//! and three of them go past an end on purpose, because the interesting bugs in a history are at
//! its edges: undoing at the start, redoing at the end, and editing after an undo.

use editor_core::edit::{History, HISTORY};

#[test]
fn a_fresh_history_has_nothing_to_undo() {
    let mut h = History::new("{}");
    assert_eq!(h.current(), "{}");
    assert!(!h.can_undo() && !h.can_redo());
    assert_eq!(h.undo(), None, "undoing at the start invented a state");
    assert_eq!(h.redo(), None);
    assert_eq!(h.current(), "{}", "a refused step still moved the cursor");
}

#[test]
fn one_edit_goes_back_and_forward() {
    let mut h = History::new("a");
    h.commit("b");
    assert_eq!(h.current(), "b");
    assert!(h.can_undo() && !h.can_redo());

    assert_eq!(h.undo(), Some("a"));
    assert_eq!(h.current(), "a");
    assert!(!h.can_undo() && h.can_redo());

    assert_eq!(h.redo(), Some("b"));
    assert_eq!(h.current(), "b");
    assert!(h.can_undo() && !h.can_redo());
}

#[test]
fn committing_the_same_text_is_not_a_step() {
    // The shell commits the current text before every edit, to fold in typing. On a scene
    // nobody typed into, that call must leave nothing behind -- an undo that appears to do
    // nothing cannot be told apart from an undo that is broken, and the user presses it again.
    let mut h = History::new("a");
    h.commit("a");
    h.commit("a");
    assert!(!h.can_undo(), "a no-op commit became an undo step");
    assert_eq!(h.len(), 1);

    h.commit("b");
    h.commit("b");
    assert_eq!(h.len(), 2);
    assert_eq!(h.undo(), Some("a"));
}

#[test]
fn typing_before_an_edit_is_folded_in_and_undone_separately() {
    // The case the two-commit protocol exists for. Start at `a`; the user types, making the
    // buffer `a typed`; then removes a domain, which splices `a typed` into `spliced`.
    //
    // Undo must walk back to `a typed` and *then* to `a`. A history that recorded only the
    // splice would jump straight to `a` and the typing would be gone with no way back to it.
    let mut h = History::new("a");

    h.commit("a typed"); // the fold-in, called before the edit
    h.commit("spliced"); // the edit itself

    assert_eq!(h.undo(), Some("a typed"), "undo skipped over the typing");
    assert_eq!(h.undo(), Some("a"));
    assert_eq!(h.undo(), None);
    assert_eq!(h.redo(), Some("a typed"));
    assert_eq!(h.redo(), Some("spliced"));
}

#[test]
fn an_edit_after_an_undo_discards_the_future() {
    // What a linear history means. Anything else is a tree, and a tree needs a way to show the
    // branches or it is a history that silently keeps states nobody can reach.
    let mut h = History::new("a");
    h.commit("b");
    h.commit("c");
    assert_eq!(h.undo(), Some("b"));
    assert!(h.can_redo());

    h.commit("d");
    assert!(!h.can_redo(), "c survived an edit made from before it");
    assert_eq!(h.current(), "d");
    assert_eq!(h.undo(), Some("b"));
    assert_eq!(h.undo(), Some("a"));
}

#[test]
fn redoing_past_the_end_does_nothing_and_stays_put() {
    let mut h = History::new("a");
    h.commit("b");
    assert_eq!(h.redo(), None);
    assert_eq!(h.current(), "b", "a refused redo moved the cursor");
    h.undo();
    assert_eq!(h.redo(), Some("b"));
    assert_eq!(h.redo(), None);
    assert_eq!(h.current(), "b");
}

#[test]
fn the_bound_drops_the_oldest_and_keeps_the_cursor_on_the_newest() {
    // Trimming the front without moving the index is the bug this is shaped to catch: the
    // cursor would keep its number and that number would name somebody else's state, so the
    // editor would show a text the user never made.
    let mut h = History::new("state 0");
    for n in 1..=HISTORY + 10 {
        h.commit(format!("state {n}"));
    }
    assert_eq!(h.len(), HISTORY, "the bound did not hold");
    assert_eq!(
        h.current(),
        format!("state {}", HISTORY + 10),
        "trimming moved the cursor off the newest state"
    );

    // And the whole window is walkable: HISTORY states means HISTORY - 1 steps back, and the
    // oldest survivor is the one arithmetic says it is.
    let mut steps = 0;
    while h.undo().is_some() {
        steps += 1;
    }
    assert_eq!(steps, HISTORY - 1);
    assert_eq!(h.current(), format!("state {}", 11));
}

#[test]
fn a_load_forgets_everything() {
    // A load is not an edit. Undoing from a freshly opened file back into the previous one
    // would offer a text belonging to no file, and the editor would then be holding a scene
    // with a path that does not describe it.
    let mut h = History::new("first");
    h.commit("first edited");
    h.reset("second");
    assert_eq!(h.current(), "second");
    assert!(!h.can_undo() && !h.can_redo());
    assert_eq!(h.len(), 1);
}

#[test]
fn a_real_splice_round_trips_through_the_history() {
    // The others are about the bookkeeping. This one runs an actual edit through it, so the
    // claim is about the thing the editor does rather than about strings named `a` and `b`.
    let scene = r#"{ "title": "t", "duration_s": 1.0, "frames": 2,
  "conservation_tolerance": 1e-9,
  "domains": [ { "kind": "bar", "name": "bar", "length_mm": 20.0, "cells": 61,
                 "area_mm2": 100.0, "initial_c": 20.0 } ] }"#;
    let mut h = History::new(scene);

    let hotter = editor_core::edit::set_number(scene, "/domains/0/initial_c", 55.0)
        .expect("initial_c is a number in this scene");
    assert!(hotter.contains("55"), "the splice did nothing");
    h.commit(&hotter);

    let gone = editor_core::edit::remove_domain(&hotter, "bar").expect("the bar is there");
    h.commit(&gone);
    assert!(!gone.contains("length_mm"), "the domain is still there");

    assert_eq!(h.undo().map(str::to_string), Some(hotter.clone()));
    assert!(
        h.current().contains("length_mm"),
        "the domain did not come back"
    );
    assert!(h.current().contains("55"), "the value edit was undone too");

    assert_eq!(h.undo().map(str::to_string), Some(scene.to_string()));
    assert!(
        !h.current().contains("55"),
        "the value edit survived its undo"
    );
}

#[test]
fn redo_after_typing_cannot_silently_discard_it() {
    // The shell folds the current text in before stepping **either** way, and this is why the
    // redo direction needs it too. Typing leaves the buffer holding a state the history has
    // never seen; stepping forward from there would overwrite it.
    //
    // Folding in first makes the typing the newest state, which truncates the tail -- so redo
    // then reports there is nothing ahead, which is true, rather than reaching a future by
    // discarding the present. Typing *ends* a redo run, and it ends it without losing anything.
    //
    // The first draft committed only before undo, and redo after typing lost it.
    let mut h = History::new("a");
    h.commit("b");
    assert_eq!(h.undo(), Some("a"));
    assert!(h.can_redo(), "b should be ahead");

    // Typing lands as a commit from the shell's next edit, and that is what closes the future.
    h.commit("a typed");
    assert!(!h.can_redo(), "b survived a commit made from before it");
    assert_eq!(h.undo(), Some("a"), "the typing is still one step back");
    assert_eq!(h.redo(), Some("a typed"), "and reachable again");
}
