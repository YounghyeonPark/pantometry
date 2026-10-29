//! **A crate went to crates.io twice without its licence, and nothing said so.**
//!
//! Every crate here is published as `MIT OR Apache-2.0`, and nineteen of the twenty carry the two
//! texts beside their manifest, which is how a `.crate` comes to hold them. `pantometry-protein`
//! never did: the packages for 0.21.0 and 0.22.0 hold `Cargo.toml`, `README.md`, `src`,
//! `structures` and `tests`, and no licence. The manifest's `license` field is an identifier, and
//! Apache-2.0 asks for the text itself to go with the work. A published version cannot be
//! changed, so those two stay as they are; this holds the next one.
//!
//! So: every directory under `crates/` with a manifest has both files, byte for byte the
//! repository's own. A copy that has drifted from the root is a crate licensed under a text the
//! project does not state.
//!
//! # Not under `wasm32`
//!
//! The files are read off a disk, and a `wasm32` target has none. `counts_in_prose` and
//! `citation_is_valid` carry the same line for the same reason.

#![cfg(not(target_family = "wasm"))]

#[test]
fn every_crate_carries_both_licence_texts_as_the_root_states_them() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("crates/pantometry has two ancestors")
        .to_path_buf();
    let names = ["LICENSE-MIT", "LICENSE-APACHE"];
    let texts: Vec<Vec<u8>> = names
        .iter()
        .map(|n| std::fs::read(root.join(n)).unwrap_or_else(|e| panic!("the root's {n}: {e}")))
        .collect();

    let mut crates: Vec<std::path::PathBuf> = std::fs::read_dir(root.join("crates"))
        .expect("crates/ is readable")
        .map(|e| e.expect("a directory entry").path())
        .filter(|p| p.join("Cargo.toml").is_file())
        .collect();
    crates.sort();
    // **A walk that found nothing would pass forever.** Twenty is the count CLAUDE.md states and
    // the workspace's members; fewer means the walk is looking in the wrong place.
    assert!(
        crates.len() >= 20,
        "only {} crates found under {}",
        crates.len(),
        root.join("crates").display()
    );

    let mut wrong = Vec::new();
    for c in &crates {
        for (n, text) in names.iter().zip(&texts) {
            match std::fs::read(c.join(n)) {
                Err(_) => wrong.push(format!("{} has no {n}", c.display())),
                Ok(b) if b != *text => {
                    wrong.push(format!("{}'s {n} is not the root's", c.display()))
                }
                Ok(_) => {}
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
