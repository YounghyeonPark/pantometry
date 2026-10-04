//! **Every refusal of [`System::from_pdb`], fed the input it refuses.**
//!
//! Each malformed file is 181L with one edit made in memory, and every edit asserts that its anchor
//! picks exactly one line, so an edit that missed cannot leave the file intact and pass. Each test
//! asserts the whole error, names included: a refusal that named the wrong residue would be one a
//! reader could not act on.

mod protein;

use pantometry_forcefield::pdb::Dropped;
use pantometry_forcefield::{Component, Element, Histidine, PdbError, Selection, System};
use protein::*;

fn refused(text: &str) -> PdbError {
    refused_with(text, &templates(), &selection())
}

fn refused_with(text: &str, templates: &[Component], selection: &Selection) -> PdbError {
    match System::from_pdb(text, templates, selection) {
        Ok(_) => panic!("the edited file was accepted"),
        Err(e) => e,
    }
}

/// The atom line for `name` of residue `number`, as it is in 181L.
fn line(number: i32, name: &str) -> &'static str {
    let hits: Vec<&str> = atom_lines(PDB_181L)
        .filter(|l| is_atom(l, number, name))
        .collect();
    assert_eq!(hits.len(), 1, "{number} {name}");
    hits[0]
}

/// `line` with columns `from` to `to` (from one) replaced by `with`, which must be as wide.
fn set_columns(line: &str, from: usize, to: usize, with: &str) -> String {
    assert_eq!(with.len(), to - from + 1);
    format!("{}{with}{}", &line[..from - 1], &line[to..])
}

/// `text` with every atom line of residue `number` deleted; at least one must go.
fn without_residue(text: &str, number: i32) -> String {
    let is = |l: &str| {
        (l.starts_with("ATOM  ") || l.starts_with("HETATM"))
            && columns(l, 23, 26).trim() == number.to_string()
    };
    assert!(text.lines().any(is), "no residue {number}");
    text.lines()
        .filter(|l| !is(l))
        .map(|l| format!("{l}\n"))
        .collect()
}

/// A dropped type of several atoms, twice: 181L's repeated types (water, chloride) are one atom
/// each, so only a second HED tells atoms counted per record from atoms counted per residue.
#[test]
fn a_dropped_type_counts_every_atom_of_every_residue() {
    let copies: String = atom_lines(PDB_181L)
        .filter(|l| columns(l, 18, 20) == "HED")
        .map(|l| format!("{}\n", set_columns(l, 23, 26, " 501")))
        .collect();
    assert_eq!(copies.lines().count(), 8);
    let c1 = line(400, "C1");
    let text = edit_one(PDB_181L, |l| l == c1, Some(&format!("{copies}{c1}")));
    let s = System::from_pdb(&text, &templates(), &selection()).unwrap();
    let hed = s.dropped().iter().find(|d| d.name == "HED").unwrap();
    assert_eq!((hed.residues, hed.atoms), (2, 16));
}

/// The unedited file is accepted, so every refusal below is the edit's.
#[test]
fn the_unedited_file_is_accepted() {
    let s = system();
    assert_eq!(
        s.dropped()
            .iter()
            .map(|d| d.name.as_str())
            .collect::<Vec<_>>(),
        ["CL", "HED", "HOH"]
    );
    assert_eq!(
        s.dropped()[1],
        Dropped {
            name: "HED".into(),
            residues: 1,
            atoms: 8
        }
    );
}

#[test]
fn a_residue_missing_a_heavy_atom_is_refused_with_every_missing_name() {
    let text = edit_one(PDB_181L, |l| is_atom(l, 70, "CB"), None);
    assert_eq!(
        refused(&text),
        PdbError::MissingAtoms {
            residue: "A:ASP70".into(),
            atoms: vec!["CB".into()]
        }
    );
    let text = edit_one(&text, |l| is_atom(l, 70, "OD1"), None);
    assert_eq!(
        refused(&text),
        PdbError::MissingAtoms {
            residue: "A:ASP70".into(),
            atoms: vec!["CB".into(), "OD1".into()]
        }
    );
    // A missing backbone C is found where the chain is joined, and named the same way.
    let text = edit_one(PDB_181L, |l| is_atom(l, 49, "C"), None);
    assert_eq!(
        refused(&text),
        PdbError::MissingAtoms {
            residue: "A:ALA49".into(),
            atoms: vec!["C".into()]
        }
    );
}

#[test]
fn a_gap_in_the_chain_is_refused() {
    let text = without_residue(PDB_181L, 50);
    match refused(&text) {
        PdbError::ChainBreak { after, before, .. } => {
            assert_eq!((after.as_str(), before.as_str()), ("A:ALA49", "A:GLY51"));
        }
        other => panic!("{other:?}"),
    }
    // Numbered consecutively but 3 Å too far apart.
    let n = line(51, "N");
    let x: f64 = columns(n, 31, 38).trim().parse().unwrap();
    let moved = set_columns(n, 31, 38, &format!("{:>8.3}", x + 3.0));
    let text = edit_one(PDB_181L, |l| is_atom(l, 51, "N"), Some(&moved));
    match refused(&text) {
        PdbError::ChainBreak {
            after,
            before,
            distance,
        } => {
            assert_eq!((after.as_str(), before.as_str()), ("A:ILE50", "A:GLY51"));
            assert!(distance > 2.0, "{distance}");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_residue_that_is_not_seqres_is_refused() {
    let first = PDB_181L
        .lines()
        .find(|l| l.starts_with("SEQRES   1 A"))
        .unwrap();
    let edited = first.replacen("  MET ASN", "  ALA ASN", 1);
    assert_ne!(edited, first);
    let text = edit_one(PDB_181L, |l| l == first, Some(&edited));
    assert_eq!(
        refused(&text),
        PdbError::SequenceMismatch {
            residue: "A:MET1".into(),
            seqres: Some("ALA".into())
        }
    );
    let text: String = PDB_181L
        .lines()
        .filter(|l| !l.starts_with("SEQRES"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(refused(&text), PdbError::NoSequence { chain: 'A' });
}

#[test]
fn an_atom_the_template_does_not_have_is_refused() {
    let nz = line(162, "NZ");
    let text = edit_one(
        PDB_181L,
        |l| l == nz,
        Some(&set_columns(nz, 13, 16, " NQ ")),
    );
    assert_eq!(
        refused(&text),
        PdbError::UnknownAtom {
            residue: "A:LYS162".into(),
            atom: "NQ".into()
        }
    );
    // An OXT on a residue that is not the C-terminus has no place in a chain.
    let o = line(100, "O");
    let oxt = set_columns(o, 13, 16, " OXT");
    let text = edit_one(PDB_181L, |l| l == o, Some(&format!("{o}\n{oxt}")));
    assert_eq!(
        refused(&text),
        PdbError::UnknownAtom {
            residue: "A:ILE100".into(),
            atom: "OXT".into()
        }
    );
}

#[test]
fn an_atom_of_the_wrong_element_is_refused() {
    let sd = line(1, "SD");
    let text = edit_one(PDB_181L, |l| l == sd, Some(&set_columns(sd, 77, 78, " O")));
    assert_eq!(
        refused(&text),
        PdbError::ElementMismatch {
            residue: "A:MET1".into(),
            atom: "SD".into(),
            file: "O".into(),
            template: Element::S
        }
    );
}

#[test]
fn a_hydrogen_in_the_file_is_refused() {
    let n = line(2, "N");
    let h = set_columns(&set_columns(n, 13, 16, " H  "), 77, 78, " H");
    let text = edit_one(PDB_181L, |l| l == n, Some(&format!("{n}\n{h}")));
    assert_eq!(
        refused(&text),
        PdbError::HydrogenInFile {
            residue: "A:ASN2".into(),
            atom: "H".into()
        }
    );
}

#[test]
fn a_hetero_residue_not_named_is_refused() {
    let s = Selection::new("BNZ").dropping(&["HOH", "CL"]);
    assert_eq!(
        refused_with(PDB_181L, &templates(), &s),
        PdbError::Unexpected {
            residue: "A:HED170".into()
        }
    );
}

#[test]
fn a_ligand_present_other_than_once_is_refused() {
    let s = Selection::new("XYZ").dropping(&["HOH", "CL", "HED", "BNZ"]);
    assert_eq!(
        refused_with(PDB_181L, &templates(), &s),
        PdbError::LigandCount {
            ligand: "XYZ".into(),
            count: 0
        }
    );
    let copies: String = atom_lines(PDB_181L)
        .filter(|l| columns(l, 18, 20) == "BNZ")
        .map(|l| format!("{}\n", set_columns(l, 23, 26, " 401")))
        .collect();
    let c1 = line(400, "C1");
    let text = edit_one(PDB_181L, |l| l == c1, Some(&format!("{copies}{c1}")));
    assert_eq!(
        refused(&text),
        PdbError::LigandCount {
            ligand: "BNZ".into(),
            count: 2
        }
    );
}

#[test]
fn a_residue_number_with_two_names_is_refused() {
    let cb = line(10, "CB");
    let text = edit_one(PDB_181L, |l| l == cb, Some(&set_columns(cb, 18, 20, "GLU")));
    assert_eq!(
        refused(&text),
        PdbError::MixedResidue {
            residue: "A:10".into(),
            names: ["ASP".into(), "GLU".into()]
        }
    );
}

#[test]
fn an_atom_given_twice_is_refused() {
    let cb = line(10, "CB");
    let text = edit_one(PDB_181L, |l| l == cb, Some(&format!("{cb}\n{cb}")));
    assert_eq!(
        refused(&text),
        PdbError::DuplicateAtom {
            residue: "A:ASP10".into(),
            atom: "CB".into()
        }
    );
}

/// The first alternate location met is kept, whichever letter it is, and the other counted.
#[test]
fn the_first_alternate_location_is_kept_and_the_other_counted() {
    let cb = line(70, "CB");
    let x: f64 = columns(cb, 31, 38).trim().parse().unwrap();
    let a = set_columns(cb, 17, 17, "A");
    let b = set_columns(
        &set_columns(cb, 17, 17, "B"),
        31,
        38,
        &format!("{:>8.3}", x + 0.5),
    );
    for (first, second, kept) in [(&a, &b, x), (&b, &a, x + 0.5)] {
        let text = edit_one(PDB_181L, |l| l == cb, Some(&format!("{first}\n{second}")));
        let s = System::from_pdb(&text, &templates(), &selection()).unwrap();
        assert_eq!(s.alternates_dropped(), 1);
        let i = (0..s.component().atoms().len())
            .find(|&i| s.residue_of(i).number == 70 && s.atom_name(i) == "CB")
            .unwrap();
        let got = s.component().atoms()[i].at[0] / 1e-10;
        assert!((got - kept).abs() < 1e-9, "{got} against {kept}");
    }
}

#[test]
fn an_insertion_code_is_refused() {
    let cb = line(10, "CB");
    let text = edit_one(PDB_181L, |l| l == cb, Some(&set_columns(cb, 27, 27, "A")));
    assert_eq!(
        refused(&text),
        PdbError::InsertionCode {
            residue: "A:ASP10A".into()
        }
    );
}

#[test]
fn more_than_one_model_is_refused() {
    let n = line(1, "N");
    let cb = line(10, "CB");
    let text = edit_one(PDB_181L, |l| l == n, Some(&format!("MODEL        1\n{n}")));
    let text = edit_one(&text, |l| l == cb, Some(&format!("MODEL        2\n{cb}")));
    assert_eq!(refused(&text), PdbError::Models { count: 2 });
}

#[test]
fn an_unreadable_coordinate_is_refused() {
    let n = line(1, "N");
    let text = edit_one(
        PDB_181L,
        |l| l == n,
        Some(&set_columns(n, 31, 38, "  4x.982")),
    );
    match refused(&text) {
        PdbError::Record { field, text, .. } => {
            assert_eq!((field, text.as_str()), ("x", "4x.982"));
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_residue_without_a_template_is_refused() {
    let without_trp: Vec<Component> = templates()
        .into_iter()
        .filter(|t| t.id() != "TRP")
        .collect();
    assert_eq!(
        refused_with(PDB_181L, &without_trp, &selection()),
        PdbError::NoTemplate {
            residue: "A:TRP126".into()
        }
    );
    let text: String = PDB_181L
        .lines()
        .map(|l| {
            if is_atom(l, 1, columns(l, 13, 16).trim()) && columns(l, 18, 20) == "MET" {
                format!("{}\n", set_columns(l, 18, 20, "MSE"))
            } else {
                format!("{l}\n")
            }
        })
        .collect();
    assert_eq!(
        refused(&text),
        PdbError::NoTemplate {
            residue: "A:MSE1".into()
        }
    );
    let mut twice = templates();
    twice.push(template("ALA"));
    assert_eq!(
        refused_with(PDB_181L, &twice, &selection()),
        PdbError::DuplicateTemplate {
            template: "ALA".into()
        }
    );
}

/// A template without the hydrogen the protonation rule removes is not the dictionary's entry.
#[test]
fn a_template_without_what_the_rules_name_is_refused() {
    let asp = ccd("ASP");
    assert!(asp.matches("HD2").count() >= 2);
    let edited = Component::from_ccd(&asp.replace("HD2", "HDX")).unwrap();
    let mut ts: Vec<Component> = templates()
        .into_iter()
        .filter(|t| t.id() != "ASP")
        .collect();
    ts.push(edited);
    assert_eq!(
        refused_with(PDB_181L, &ts, &selection()),
        PdbError::TemplateMismatch {
            template: "ASP".into(),
            what: "atom HD2".into()
        }
    );
}

/// The template a δ-histidine's double bond moves across must have both bonds it moves.
#[test]
fn a_histidine_template_without_its_ring_bond_is_refused() {
    let his = ccd("HIS");
    let row = "HIS CE1 NE2 SING Y N 18 \n";
    assert_eq!(his.matches(row).count(), 1);
    let edited = Component::from_ccd(&his.replace(row, "")).unwrap();
    let mut ts: Vec<Component> = templates()
        .into_iter()
        .filter(|t| t.id() != "HIS")
        .collect();
    ts.push(edited);
    assert_eq!(
        refused_with(
            PDB_181L,
            &ts,
            &selection().histidine('A', 31, Histidine::Delta)
        ),
        PdbError::TemplateMismatch {
            template: "HIS".into(),
            what: "the ND1-CE1 and CE1-NE2 bonds".into()
        }
    );
}

/// An N-terminal template that already has an atom named H3 has no name free for the third
/// hydrogen: Met with its HE3 renamed H3.
#[test]
fn a_terminal_template_with_an_h3_already_is_refused() {
    let met = ccd("MET");
    assert_eq!(met.matches("MET HE3 3HE").count(), 1);
    assert_eq!(met.matches("MET CE  HE3").count(), 1);
    let edited = met
        .replace("MET HE3 3HE", "MET H3  3HE")
        .replace("MET CE  HE3", "MET CE  H3 ");
    let edited = Component::from_ccd(&edited).unwrap();
    let mut ts: Vec<Component> = templates()
        .into_iter()
        .filter(|t| t.id() != "MET")
        .collect();
    ts.push(edited);
    assert_eq!(
        refused_with(PDB_181L, &ts, &selection()),
        PdbError::TemplateMismatch {
            template: "MET".into(),
            what: "free name H3 for the third N-terminal hydrogen".into()
        }
    );
}

#[test]
fn an_n_terminal_proline_is_refused() {
    let mut text = PDB_181L.to_string();
    for n in 1..=36 {
        text = without_residue(&text, n);
    }
    assert_eq!(
        refused(&text),
        PdbError::UnsupportedTerminus {
            residue: "A:PRO37".into()
        }
    );
}

#[test]
fn a_histidine_that_is_not_one_is_refused() {
    let s = selection().histidine('A', 32, Histidine::Both);
    assert_eq!(
        refused_with(PDB_181L, &templates(), &s),
        PdbError::NoSuchHistidine {
            residue: "A:32".into()
        }
    );
}

/// Methanol as the ligand: the hydrogens on its carbon have a parent with one heavy neighbour
/// and nothing beyond it, two points, which fix no orientation.
#[test]
fn a_hydrogen_no_superposition_can_place_is_refused() {
    let methanol = [
        atom_line(true, 9001, "C", "MOH", 'A', 500, [60.0, 60.0, 60.0], "C"),
        atom_line(true, 9002, "O", "MOH", 'A', 500, [61.43, 60.0, 60.0], "O"),
    ]
    .join("\n");
    let c1 = line(400, "C1");
    let text = edit_one(PDB_181L, |l| l == c1, Some(&format!("{methanol}\n{c1}")));
    let mut ts = templates();
    ts.push(Component::from_ccd(include_str!("../components/MOH.cif")).unwrap());
    let s = Selection::new("MOH").dropping(&["HOH", "CL", "HED", "BNZ"]);
    match refused_with(&text, &ts, &s) {
        PdbError::Unplaceable { residue, atom } => {
            assert_eq!(residue, "A:MOH500");
            assert!(atom.starts_with('H'), "{atom}");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_file_with_no_protein_is_refused() {
    let text: String = atom_lines(PDB_181L)
        .filter(|l| columns(l, 18, 20) == "BNZ")
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(refused(&text), PdbError::NoProtein);
}

/// Every refusal says what it refused, in words that name it.
#[test]
fn every_refusal_names_what_it_refused() {
    let text = edit_one(PDB_181L, |l| is_atom(l, 70, "CB"), None);
    let message = refused(&text).to_string();
    assert!(
        message.contains("A:ASP70") && message.contains("CB"),
        "{message}"
    );
}
