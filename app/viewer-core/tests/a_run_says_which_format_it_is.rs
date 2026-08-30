//! **A run states its format, and a reader that cannot read it says so in words.**
//!
//! The run format went unversioned for as long as nothing wanted to add a shape to it. The
//! moment something did, the failure it would have produced became visible: `PanelData` is
//! tagged and `deny_unknown_fields`, so a run carrying a panel kind this build has never heard
//! of fails inside serde with a complaint about an unknown variant — which is also what a
//! truncated file, a half-written file and a file of the wrong kind entirely say.
//!
//! Those want different words. "This run is newer than your pantometry" is a thing a person can
//! act on; "unknown variant `surface`, expected one of ..." is a thing they file a bug about.
//!
//! So the version is read **first**, from a reader that ignores every other key, and the checks
//! here are about that ordering as much as about the number.

use viewer_core::{Run, FORMAT};

/// A minimal run at whatever `format` line is given, with one panel this build does understand.
fn run_with(format_line: &str) -> String {
    format!(
        r#"{{ {format_line}"title": "a run",
  "frames": [
    {{ "t": 0.0, "panels": [
      {{ "name": "bar", "unit": "K", "kind": "field",
         "nx": 2, "ny": 1, "nz": 1, "values": [300.0, 301.0] }}
    ], "readings": [] }}
  ] }}"#
    )
}

#[test]
fn a_run_without_a_format_is_format_one() {
    // Every run written before the key existed. There are recorded ones in `tests/runs/` and
    // they are the reason this is not a formality: a reader that demanded the key would refuse
    // every file this workspace has ever produced.
    let run = Run::from_json(&run_with("")).expect("an unstamped run still reads");
    assert_eq!(run.title, "a run");
    assert_eq!(run.frames.len(), 1);
}

#[test]
fn the_format_this_build_writes_is_read() {
    let run = Run::from_json(&run_with(&format!("\"format\": {FORMAT}, ")))
        .expect("the current format reads");
    assert_eq!(run.frames.len(), 1);
}

#[test]
fn a_run_from_the_future_is_refused_by_number() {
    let err = Run::from_json(&run_with("\"format\": 99, ")).expect_err("99 must be refused");
    assert!(err.contains("99"), "the error does not say which: {err}");
    assert!(
        err.to_lowercase().contains("upgrade"),
        "the error does not say what to do: {err}"
    );
}

#[test]
fn zero_is_not_a_version_this_format_ever_had() {
    // What an uninitialised field looks like, and what a writer that meant to stamp a run and
    // did not would leave. Accepting it would make "no version" and "a broken version" the same
    // answer, which is the distinction this whole file is about.
    let err = Run::from_json(&run_with("\"format\": 0, ")).expect_err("0 must be refused");
    assert!(err.contains('0'), "{err}");
}

#[test]
fn the_version_is_read_before_the_panels() {
    // **The ordering, which is the point.** A future run will carry shapes this build does not
    // know, and serde would fail on the unknown variant before anything looked at the version.
    // Stamping it 99 must produce the version message and not the variant one.
    let future = r#"{ "format": 99, "title": "a run",
  "frames": [
    { "t": 0.0, "panels": [
      { "name": "part", "unit": "m", "kind": "a-shape-from-later",
        "whatever": [1, 2, 3] }
    ], "readings": [] }
  ] }"#;
    let err = Run::from_json(future).expect_err("a future run must be refused");
    assert!(
        err.contains("99"),
        "the version was not what refused it: {err}"
    );
    assert!(
        !err.contains("unknown variant"),
        "serde got there first, which is the bug this ordering exists to prevent: {err}"
    );
}

#[test]
fn a_file_that_is_not_a_run_at_all_still_says_so() {
    // The version check must not turn every failure into a version failure. A truncated file has
    // no `format` to read either, and the message it deserves is the old one.
    for bad in ["{", "[]", "not json", r#"{"title": 3}"#] {
        let err = Run::from_json(bad).expect_err("{bad} is not a run");
        assert!(
            err.contains("not a pantometry run"),
            "{bad:?} produced {err}"
        );
    }
}
