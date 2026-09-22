//! **A run holding a protein used to be a point cloud.**
//!
//! `Bodies` gave a count, a position and one number per body, so a file holding forty-six points
//! and forty-six numbers could not be matched against the entry it was read from, against a
//! crystallographer's temperature factors, or against a second structure. `Residue` has the chain,
//! the number and the name; `Network` keeps only geometry, and the identity was dropped at that
//! boundary. That is the difference between a picture and a measurement.
//!
//! `Protein::with_residues` carries it across, and the backbone comes out as `Bodies::bonds`.
//! These check both against the file rather than against a second implementation: crambin's own
//! sequence, and the numbering of a structure with a hole in it.
//!
//! # Why a bond needs all three conditions
//!
//! Same chain, numbered one apart, and within `PEPTIDE_REACH` of each other. Each covers what the
//! others miss — a disordered loop leaves a gap the numbering shows and the distance may not, two
//! fragments that end near each other are close without being joined, and a file numbered with
//! insertion codes repeats a number where the chain does continue. The arrangement prefers a
//! **missing** bond, which is visible in a picture, over a drawn one nobody measured.

use pantometry_protein::{Network, Protein, Structure};
use pantometry_units::{Mass, Qty, Temperature};

/// Ångström, in metres.
const A: f64 = 1e-10;
/// `1CRN`, which is 46 residues of one chain with no gaps in it.
const CRAMBIN: &str = include_str!("../structures/1CRN.pdb");

/// A protein at body temperature, built the way a scene builds one.
fn shaking(structure: &Structure) -> Protein {
    let network = Network::new(structure, Qty::from_si(15.0 * A), Qty::from_si(1.0));
    Protein::new(
        "p",
        &network,
        Temperature::from_si(310.15),
        Mass::from_si(1.8e-25),
        1,
    )
    .with_residues(structure)
}

/// **Crambin comes out as crambin**, in the file's own words.
#[test]
fn every_body_says_which_residue_it_is() {
    use pantometry_core::Bodies;
    let structure = Structure::from_pdb(CRAMBIN).expect("crambin parses");
    let p = shaking(&structure);

    assert_eq!(p.count(), 46, "crambin is 46 residues");
    // The first six of `1CRN`'s sequence, which is T T C C P S — a fact about the molecule and
    // not about this code. A reader with the entry open can check it by eye.
    let first: Vec<String> = (0..6).map(|i| p.label(i).expect("labelled")).collect();
    assert_eq!(
        first,
        ["A:1:THR", "A:2:THR", "A:3:CYS", "A:4:CYS", "A:5:PRO", "A:6:SER"],
        "these are the residues the file names, in the order it names them"
    );
    assert_eq!(p.label(45).as_deref(), Some("A:46:ASN"), "and the last one");

    // Every body labelled, because a half-labelled set is worse than none: a reader matching on
    // the labels silently loses the bodies that have none.
    assert!(
        (0..p.count()).all(|i| p.label(i).is_some()),
        "a body with no label is a body a reader cannot match"
    );

    // **One unbroken chain**, which is what 46 residues numbered 1 to 46 with no gaps means.
    let bonds = p.bonds();
    assert_eq!(bonds.len(), 45, "45 bonds join 46 residues in one chain");
    for (i, b) in bonds.iter().enumerate() {
        assert_eq!(
            *b,
            [i as u32, i as u32 + 1],
            "the chain is in sequence order"
        );
    }
}

/// **A gap in the numbering is not a bond.**
///
/// A crystal structure whose density ran out over a loop has residues on both sides of it and
/// nothing in between. Joining them draws a bond nobody measured, across a distance nobody saw —
/// and in a picture that reads as a chain rather than as two fragments.
#[test]
fn a_hole_in_the_numbering_is_not_bonded_across() {
    // Four alpha carbons in a line, 3.8 Å apart, numbered 1, 2, **5**, 6. The geometry is
    // continuous on purpose: only the numbering says the chain is not.
    let text = "\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C
ATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00 10.00           C
ATOM      3  CA  ALA A   5       7.600   0.000   0.000  1.00 10.00           C
ATOM      4  CA  ALA A   6      11.400   0.000   0.000  1.00 10.00           C
";
    let structure = Structure::from_pdb(text).expect("four alpha carbons");
    let p = shaking(&structure);
    use pantometry_core::Bodies;

    assert_eq!(
        p.bonds(),
        vec![[0, 1], [2, 3]],
        "residues 2 and 5 are 3.8 Å apart and three residues apart; the distance alone would \
         have joined them"
    );
}

/// **Two chains do not join where they meet, and only the chain letter says so.**
///
/// An asymmetric unit holding two copies puts the end of one near the start of the next: as close
/// as any two bonded residues, and not bonded.
///
/// **Chain B is numbered 3 and 4 on purpose.** The first version of this test numbered it 1 and 2,
/// so `B:1` did not follow `A:2` in *numbering* either — and the numbering condition blocked the
/// bond while the chain condition never ran. Measured: deleting `a.0 == b.0` from `bonds` left all
/// six tests in this file green. Numbering B from 3 leaves the chain letter as the only thing
/// between these two residues and a bond across them.
#[test]
fn two_chains_do_not_join_where_they_touch() {
    let text = "\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C
ATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00 10.00           C
ATOM      3  CA  ALA B   3       7.600   0.000   0.000  1.00 10.00           C
ATOM      4  CA  ALA B   4      11.400   0.000   0.000  1.00 10.00           C
";
    let structure = Structure::from_pdb(text).expect("two chains of two");
    let p = shaking(&structure);
    use pantometry_core::Bodies;

    assert_eq!(
        p.bonds(),
        vec![[0, 1], [2, 3]],
        "A:2 and B:3 are 3.8 Å apart and numbered one apart; they are in different chains"
    );
    assert_eq!(
        p.label(2).as_deref(),
        Some("B:3:ALA"),
        "and the chain is in the label"
    );
}

/// **Consecutive numbering across a distance nothing spans is not a bond either.**
///
/// The numbering is a convention and a file can be wrong about it; the peptide bond is 3.8 Å and
/// cannot be 12. This is the backstop for a structure numbered as one chain across a piece that
/// was never modelled.
#[test]
fn residues_too_far_apart_are_not_bonded_however_they_are_numbered() {
    let text = "\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C
ATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00 10.00           C
ATOM      3  CA  ALA A   3      15.800   0.000   0.000  1.00 10.00           C
";
    let structure = Structure::from_pdb(text).expect("three alpha carbons");
    let p = shaking(&structure);
    use pantometry_core::Bodies;

    assert_eq!(
        p.bonds(),
        vec![[0, 1]],
        "residues 2 and 3 are numbered as neighbours and are 12 Å apart, which no peptide bond is"
    );
}

/// **A protein with no residues supplied is honest about it** rather than inventing names.
#[test]
fn an_unlabelled_protein_says_nothing_rather_than_guessing() {
    use pantometry_core::Bodies;
    let structure = Structure::from_pdb(CRAMBIN).expect("crambin parses");
    let network = Network::new(&structure, Qty::from_si(15.0 * A), Qty::from_si(1.0));
    let bare = Protein::new(
        "p",
        &network,
        Temperature::from_si(310.15),
        Mass::from_si(1.8e-25),
        1,
    );
    assert_eq!(bare.label(0), None, "nothing said is nothing invented");
    assert!(
        bare.bonds().is_empty(),
        "and a chain nobody described is not a chain"
    );
}

/// **Labelling from a structure of the wrong length is refused.**
///
/// It would name every body wrongly and say nothing about it, which is the failure this whole
/// change exists to remove.
#[test]
#[should_panic(expected = "would name every one of them wrongly")]
fn a_structure_of_the_wrong_length_is_refused() {
    let structure = Structure::from_pdb(CRAMBIN).expect("crambin parses");
    let network = Network::new(&structure, Qty::from_si(15.0 * A), Qty::from_si(1.0));
    let shorter = structure
        .chain('A')
        .expect("crambin has a chain A")
        .residues()[..20]
        .iter()
        .map(|r| r.at)
        .collect::<Vec<_>>();
    Protein::new(
        "p",
        &network,
        Temperature::from_si(310.15),
        Mass::from_si(1.8e-25),
        1,
    )
    .with_residues(&Structure::from_positions(shorter));
}
