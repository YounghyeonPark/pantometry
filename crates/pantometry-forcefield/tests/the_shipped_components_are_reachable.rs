//! `pantometry_forcefield::components`: the entries this crate ships, reachable by identifier.
//!
//! The list is written out in `src/components.rs`, so what holds it to the directory is this
//! test, which reads `components/` and compares. Each listed text is also parsed, and must be the
//! entry its key names.

use pantometry_forcefield::{components, Component};

/// Read off the disk, so not on a target with none.
#[cfg(not(target_family = "wasm"))]
#[test]
fn every_shipped_component_is_listed() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("components");
    let mut on_disk: Vec<String> = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()))
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter_map(|n| n.strip_suffix(".cif").map(str::to_owned))
        .collect();
    on_disk.sort();
    assert!(on_disk.len() >= 52, "only {} entries found", on_disk.len());
    let listed: Vec<String> = components::ALL
        .iter()
        .map(|(id, _)| id.to_string())
        .collect();
    assert_eq!(listed, on_disk, "components::ALL is not the directory");
    // And each listed text is that file's bytes.
    for (id, text) in components::ALL {
        let file = std::fs::read_to_string(dir.join(format!("{id}.cif"))).unwrap();
        assert_eq!(text, file, "{id}");
    }
}

#[test]
fn each_entry_is_the_component_its_key_names() {
    let mut previous = "";
    for (id, text) in components::ALL {
        assert!(id > previous, "{id} is out of order or repeated");
        previous = id;
        let c = Component::from_ccd(text).unwrap_or_else(|e| panic!("{id}: {e}"));
        assert_eq!(c.id(), id);
        assert_eq!(components::text(id), Some(text));
    }
}

#[test]
fn a_lookup_is_exact() {
    assert_eq!(
        components::text("BNZ"),
        Some(include_str!("../components/BNZ.cif"))
    );
    assert_eq!(components::text("bnz"), None);
    assert_eq!(components::text("BNZ "), None);
    assert_eq!(components::text(""), None);
    assert_eq!(
        components::text("OCT"),
        None,
        "n-octane is not shipped; OC9 is octanol"
    );
}
