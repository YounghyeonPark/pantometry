//! **`viewer-core` does not link the library, and that is a claim about the run file.**
//!
//! If a viewer can be written against the **file** and nothing else, the wire format is complete.
//! If it needed the library for something the file did not carry, the format would be the thing to
//! fix — and finding that out is worth more than the convenience. It has been worth it once
//! already: a field carried a grid and not the extent it was sampled over, so a 9×9×9 block framed
//! as a nine-metre cube at the origin. The fix went into the *format*.
//!
//! # Why this is a test now
//!
//! It used to be a **workspace boundary**: `runtime/viewer` was its own cargo workspace with its own
//! lockfile, and nothing in it could reach `crates/` without somebody adding a path dependency and
//! noticing. The viewer, the editor, the CLI and the GPU accelerator are one workspace now — see
//! `app/Cargo.toml` for why — and in one workspace `pantometry = { workspace = true }` is one line
//! away in any manifest. A property that used to be structural is a discipline, and a discipline
//! that nothing checks is a comment.
//!
//! The editor is the crate that makes this live rather than theoretical: `editor-core` links
//! `pantometry` deliberately, sits beside this one, and shares a `[workspace.dependencies]` table
//! with it.

use std::path::Path;

/// Every dependency `viewer-core` declares, from its own manifest.
fn declared() -> Vec<String> {
    let text = std::fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml"))
        .expect("this crate's manifest is beside it");
    let mut out = Vec::new();
    let mut inside = false;
    for line in text.lines() {
        let line = line.trim();
        if line.starts_with('[') {
            // Every dependency table, including `dev-` and `build-`: a dev-dependency on the
            // library would let a *test* here be written against something the file does not
            // carry, which is the same claim going quietly wrong.
            inside = line.ends_with("dependencies]");
            continue;
        }
        if !inside || line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some((name, _)) = line.split_once('=') {
            out.push(name.trim().trim_matches('"').to_string());
        }
    }
    out
}

/// Nothing from `crates/`, by any spelling.
#[test]
fn it_names_nothing_from_the_library() {
    let deps = declared();
    assert!(
        !deps.is_empty(),
        "the manifest parser found no dependencies at all, which means it stopped working \
         rather than that the crate has none"
    );
    println!("  viewer-core declares: {}", deps.join(", "));

    let reaching: Vec<&String> = deps
        .iter()
        .filter(|d| d.starts_with("pantometry") || d.as_str() == "editor-core")
        .collect();
    assert!(
        reaching.is_empty(),
        "viewer-core reaches into the library through {reaching:?} — if the run file is missing \
         something, the file is what to fix"
    );
}

/// And it stays small enough that the claim means something.
///
/// Two serde crates. A viewer that had grown a dozen dependencies could satisfy the test above and
/// still have stopped being an independent reading of the format.
#[test]
fn it_is_still_two_serde_crates() {
    let deps = declared();
    assert!(
        deps.len() <= 3,
        "viewer-core has grown to {} dependencies: {deps:?}",
        deps.len()
    );
    for d in &deps {
        assert!(
            d.starts_with("serde"),
            "{d} is not serde — see this test's docs before adding it"
        );
    }
}
