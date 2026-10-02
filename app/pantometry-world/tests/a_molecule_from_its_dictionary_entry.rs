//! A `molecule` reads one Chemical Component Dictionary entry, and refuses what it cannot run.
//!
//! The scene that uses it, `33-aspirin-relaxing-out-of-a-clash`, is checked in `scene.rs` against
//! what minimisation guarantees. What is here is the format's half: the entry beside the scenes is
//! the one the library is checked against, and an entry the build cannot honour is an error naming
//! the domain and the file rather than a run that does nothing.
//!
//! # Not under `wasm32`
//!
//! It reads the shipped entry off a disk, and a `wasm32` target has none.

#![cfg(not(target_family = "wasm"))]

use pantometry_world::{Scene, Uploaded, World};

/// A scene of one molecule called `drug`, reading `name`.
fn scene(name: &str) -> Scene {
    serde_json::from_str(&format!(
        r#"{{ "title": "a molecule", "duration_s": 3.0, "frames": 3,
             "domains": [ {{ "kind": "molecule", "name": "drug", "ccd": "{name}" }} ] }}"#
    ))
    .expect("the scene parses")
}

/// Dimethyl sulfone, written out: C1–S1–C2 with two S=O. Its sulfur's bond orders sum to six,
/// which is UFF's `S_3+6`, and this force field has no parameters for it.
const SULFONE: &str = r#"data_TST
_chem_comp.id TST
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
TST C1 C 0 N -1.5 0.0 0.0
TST C2 C 0 N 1.5 0.0 0.0
TST S1 S 0 N 0.0 0.5 0.0
TST O1 O 0 N 0.0 1.2 1.2
TST O2 O 0 N 0.0 1.2 -1.2
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
C1 S1 SING N
S1 C2 SING N
S1 O1 DOUB N
S1 O2 DOUB N
"#;

/// **The entry beside the scenes is the one the library crate is checked against.**
///
/// `crates/pantometry-forcefield/components/AIN.cif` is a test fixture and this is a scene's input,
/// and they are the same dictionary entry as RCSB serves it. Two copies of a file are two files
/// that can drift, and the one that would drift is the one nothing compares — the reason
/// `the_two_copies_of_crambin_are_one_file` exists, and this is the same check.
#[test]
fn the_two_copies_of_aspirin_are_one_file() {
    let here = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let scene_copy = here.join("scenes/structures/AIN.cif");
    let crate_copy = here.join("../../crates/pantometry-forcefield/components/AIN.cif");
    let (a, b) = (
        std::fs::read(&scene_copy).expect("the scene's copy"),
        std::fs::read(&crate_copy).expect("the crate's copy"),
    );
    assert_eq!(
        a.len(),
        b.len(),
        "{} is {} bytes and {} is {}",
        scene_copy.display(),
        a.len(),
        crate_copy.display(),
        b.len()
    );
    assert!(a == b, "the two copies of AIN.cif have diverged");
    // And neither is empty or a stub, which two missing or truncated files would also satisfy.
    let text = String::from_utf8(a).expect("text");
    assert!(
        text.contains("_chem_comp.formula") && text.contains("C9 H8 O4"),
        "AIN.cif does not state aspirin's formula"
    );
    println!("  both copies of AIN.cif are {} bytes", text.len());
}

/// **A molecule the force field refuses is refused at build, by name.** `Molecule::new` never
/// fails — a library caller who asked for the molecule gets it, drawn, with NaN energies — but a
/// scene asked for a relaxation, and a run in which nothing moves and the only sign is a NaN is
/// the failure this format refuses everywhere else.
#[test]
fn a_molecule_the_force_field_cannot_describe_is_refused() {
    let files = Uploaded::new().with("sulfone.cif", SULFONE);
    let refused = match World::build_with(scene("sulfone.cif"), &files) {
        Ok(_) => panic!("a sulfone built, and would have run without moving"),
        Err(e) => e,
    };
    for said in ["drug", "sulfone.cif", "S1", "S_3+6"] {
        assert!(
            refused.contains(said),
            "the refusal does not say {said:?}: {refused}"
        );
    }

    // The control: the same skeleton with no oxygens is dimethyl sulfide, which UFF describes, so
    // the refusal is the sulfur's valence and not the file's shape.
    let sulfide: String = SULFONE
        .lines()
        .filter(|l| !l.contains("O1") && !l.contains("O2"))
        .map(|l| format!("{l}\n"))
        .collect();
    let files = Uploaded::new().with("sulfide.cif", sulfide);
    World::build_with(scene("sulfide.cif"), &files).expect("dimethyl sulfide builds");
}

/// **A missing entry and one that is not an entry are both refused naming the domain and the
/// file**, the way a missing `pdb` or `stl` is.
#[test]
fn a_missing_or_malformed_entry_is_refused_by_name() {
    let missing = match World::build_with(scene("absent.cif"), &Uploaded::new()) {
        Ok(_) => panic!("a molecule built from a file nobody supplied"),
        Err(e) => e,
    };
    assert!(
        missing.contains("drug") && missing.contains("absent.cif"),
        "{missing}"
    );

    let files = Uploaded::new().with("notes.cif", "this is not a dictionary entry\n");
    let malformed = match World::build_with(scene("notes.cif"), &files) {
        Ok(_) => panic!("a molecule built from prose"),
        Err(e) => e,
    };
    assert!(
        malformed.contains("drug") && malformed.contains("notes.cif"),
        "{malformed}"
    );
}
