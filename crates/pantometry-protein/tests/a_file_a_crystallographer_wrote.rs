//! The PDB reader against the things the format does that a reader can miss quietly.
//!
//! Every failure here produces a structure rather than an error, which is why each has a test of
//! its own: a calcium ion taken as a residue is one extra node, an alternate location taken twice
//! is one duplicated node, and both leave a network that diagonalises, reports six rigid modes,
//! and predicts fluctuations for a protein that is not the one in the file.

use pantometry_protein::{ParseError, Structure};

/// Ångström, in metres.
const A: f64 = 1e-10;

/// Two residues of a real backbone, in the columns the format puts them in.
const TWO: &str = "\
ATOM      1  N   MET A   1      20.154  29.699   5.276  1.00 49.05           N
ATOM      2  CA  MET A   1      21.618  29.291   5.205  1.00 48.22           C
ATOM      3  C   MET A   1      22.219  29.591   6.569  1.00 47.51           C
ATOM      9  CA  GLY B   7      25.112  30.010   7.881  1.00 12.34           C
";

/// **The columns are read where the format puts them.**
#[test]
fn a_record_is_read_out_of_its_columns() {
    let s = Structure::from_pdb(TWO).expect("two alpha carbons");
    assert_eq!(
        s.len(),
        2,
        "the N and C backbone atoms were taken as residues"
    );

    let met = &s.residues()[0];
    assert_eq!(met.name, "MET");
    assert_eq!(met.chain, 'A');
    assert_eq!(met.number, 1);
    for (got, want) in met.at.iter().zip(&[21.618 * A, 29.291 * A, 5.205 * A]) {
        assert!((got - want).abs() < 1e-16, "{got:e} against {want:e}");
    }
    // 48.22 Å² of B-factor, in m².
    assert!((met.b_factor / (48.22 * A * A) - 1.0).abs() < 1e-12);

    let gly = &s.residues()[1];
    assert_eq!((gly.name.as_str(), gly.chain, gly.number), ("GLY", 'B', 7));
    println!("  two chains, two residues, MET A1 at {:?}", met.at);
}

/// **A calcium ion writes `CA` in the atom-name column and is not an alpha carbon.**
///
/// Told apart two ways, because a file may use either: the record is `HETATM` rather than `ATOM`,
/// and the element column says `CA` rather than `C`. A reader that keyed on the name alone gains
/// one node per bound calcium, at a place no residue is, and nothing downstream notices.
#[test]
fn a_calcium_ion_is_not_a_residue() {
    let with_metal = format!(
        "{TWO}\
HETATM  200  CA   CA A 301      30.000  30.000  30.000  1.00 20.00          CA
ATOM    201  CA   CA A 302      40.000  40.000  40.000  1.00 20.00          CA
"
    );
    let s = Structure::from_pdb(&with_metal).expect("the two alpha carbons");
    assert_eq!(
        s.len(),
        2,
        "a calcium was taken as a residue: {:?}",
        s.residues()
            .iter()
            .map(|r| (&r.name, r.number))
            .collect::<Vec<_>>()
    );
    println!("  a HETATM calcium and an ATOM calcium, both refused");
}

/// **A residue modelled in two places is one residue.**
///
/// `altLoc` `A` and `B` are two conformations of the same side chain, and both carry an alpha
/// carbon a fraction of an ångström apart. Taking both puts two nodes almost on top of each
/// other, and a pair of nodes that close has an antisymmetric motion that costs nearly nothing —
/// so the spectrum gains a near-zero eigenvalue that reads exactly like a soft collective mode.
#[test]
fn a_residue_modelled_twice_is_one_node() {
    let both = "\
ATOM      2  CA AMET A   1      21.618  29.291   5.205  0.60 48.22           C
ATOM      3  CA BMET A   1      21.701  29.180   5.310  0.40 50.10           C
ATOM      9  CA  GLY A   2      25.112  30.010   7.881  1.00 12.34           C
";
    let s = Structure::from_pdb(both).expect("two residues");
    assert_eq!(s.len(), 2, "both alternate locations were taken");
    // The first one, which is the one with the higher occupancy in every file that orders them.
    assert!((s.residues()[0].at[0] / (21.618 * A) - 1.0).abs() < 1e-12);
    println!("  altLoc A and B collapse to one node");
}

/// **An NMR ensemble is twenty copies of the same protein, and the first one is the structure.**
#[test]
fn only_the_first_model_is_read() {
    let ensemble = format!("MODEL        1\n{TWO}ENDMDL\nMODEL        2\n{TWO}ENDMDL\n");
    let s = Structure::from_pdb(&ensemble).expect("the first model");
    assert_eq!(s.len(), 2, "every model's atoms were piled together");
    println!("  two models, one structure");
}

/// **A file with no alpha carbons is an error and not an empty structure.**
///
/// The quiet failure: a file of the wrong kind, or one that is all `HETATM`, parses without
/// complaint into a list of nothing. An empty network has an empty spectrum, no rigid modes, no
/// fluctuations and no error anywhere in the chain.
#[test]
fn a_file_with_nothing_in_it_says_so() {
    assert_eq!(Structure::from_pdb(""), Err(ParseError::NoAlphaCarbons));
    assert_eq!(
        Structure::from_pdb("HEADER    HYDROLASE\nREMARK   2 RESOLUTION.    1.50 ANGSTROMS.\n"),
        Err(ParseError::NoAlphaCarbons)
    );
    assert_eq!(
        Structure::from_pdb(
            "HETATM  200  CA   CA A 301      30.000  30.000  30.000  1.00 20.00          CA  \n"
        ),
        Err(ParseError::NoAlphaCarbons)
    );

    // A line that begins ATOM and stops before the coordinates is a truncated file, which is a
    // different thing from a file about something else.
    let short = Structure::from_pdb("ATOM      2  CA  MET A   1      21.618\n");
    assert_eq!(short, Err(ParseError::ShortRecord { line: 1 }));
    println!("  empty, wrong-kind, all-HETATM and truncated all refuse");
}

/// **Two residues at the same place make the network singular in a way that reads as physics.**
#[test]
fn two_residues_in_one_place_are_refused() {
    let stacked = "\
ATOM      2  CA  MET A   1      21.618  29.291   5.205  1.00 48.22           C
ATOM      9  CA  GLY A   2      21.618  29.291   5.205  1.00 12.34           C
";
    assert_eq!(
        Structure::from_pdb(stacked),
        Err(ParseError::Coincident {
            first: 0,
            second: 1
        })
    );
    println!("  coincident residues refused");
}

/// **The B-factors come back as mean-square displacements, and their absence comes back as
/// absence.**
///
/// `B = 8π²⟨Δr²⟩/3`, inverted. A file without the column would otherwise give a vector of zeros,
/// which correlates with nothing and whose correlation coefficient is `NaN` — and `NaN` compares
/// false against every threshold, so a test written with `<` passes on it.
#[test]
fn the_temperature_factors_become_displacements_or_nothing() {
    let s = Structure::from_pdb(TWO).expect("two residues");
    let f = s
        .experimental_fluctuations()
        .expect("a file with B-factors");
    assert_eq!(f.len(), 2);
    // 48.22 Å² -> 3/(8π²) × 48.22 Å² = 1.832 Å².
    let expected = 3.0 / (8.0 * std::f64::consts::PI * std::f64::consts::PI) * 48.22 * A * A;
    assert!(
        (f[0] / expected - 1.0).abs() < 1e-12,
        "{:e} against {expected:e}",
        f[0]
    );

    // A file whose B-factor column is all zeros, which is what a model deposited without
    // refinement looks like.
    let flat = "\
ATOM      2  CA  MET A   1      21.618  29.291   5.205  1.00  0.00           C
ATOM      9  CA  GLY A   2      25.112  30.010   7.881  1.00  0.00           C
";
    assert_eq!(
        Structure::from_pdb(flat)
            .expect("two residues")
            .experimental_fluctuations(),
        None
    );
    println!(
        "  B-factors become {:.3} A^2, and a flat column becomes None",
        f[0] / (A * A)
    );
}

/// **A structure superposed onto a rotated copy of itself comes back exactly.**
///
/// The closed form a fitting routine has: move a structure by a rotation and a translation
/// somebody chose, fit it back, and the residual is not "small" but **zero to rounding**, because
/// an exact fit exists. A fitting routine that is subtly wrong — a transposed correlation matrix,
/// the wrong eigenvector of the four, a quaternion convention reversed — still produces a
/// rotation, still reduces the distance, and fails this.
#[test]
fn a_rotated_copy_fits_back_onto_the_original() {
    let s = Structure::from_pdb(TWO).expect("two residues");
    let long = Structure::from_positions(
        (0..40)
            .map(|i| {
                let turn = i as f64 * 100.0f64.to_radians();
                [
                    2.3 * A * turn.cos(),
                    2.3 * A * turn.sin(),
                    1.5 * A * i as f64,
                ]
            })
            .collect(),
    );

    // A rotation nothing here can simplify, about an axis that is not a coordinate axis, and a
    // translation far larger than the structure.
    let (c, sn) = (0.8f64.cos(), 0.8f64.sin());
    let moved = Structure::from_positions(
        long.residues()
            .iter()
            .map(|r| {
                let (x, y, z) = (r.at[0], r.at[1], r.at[2]);
                let (y1, z1) = (c * y - sn * z, sn * y + c * z);
                let (x2, y2) = (
                    0.5f64.cos() * x - 0.5f64.sin() * y1,
                    0.5f64.sin() * x + 0.5f64.cos() * y1,
                );
                [x2 + 137.0 * A, y2 - 61.0 * A, z1 + 909.0 * A]
            })
            .collect(),
    );

    let apart = long.rmsd_to(&moved).expect("the same length");
    let back = moved.superposed_onto(&long).expect("the same length");
    let left = long.rmsd_to(&back).expect("the same length");
    assert!(apart > 100.0 * A, "the copy was not moved: {:e}", apart);
    assert!(
        left < 1e-15 * apart,
        "an exact fit exists and {:e} m of it was left, against a starting {:e}",
        left,
        apart
    );

    // Onto itself, and onto something of another size.
    assert_eq!(long.rmsd_to(&long), Some(0.0));
    assert!(long.superposed_onto(&s).is_none());
    assert!(
        long.direction_to(&long).is_none(),
        "a structure has no direction to itself"
    );
    println!("  a rotated, translated copy fits back with {left:.3e} m left of {apart:.3e}");
}

/// **A chain is selected, and a chain that is not there is `None`.**
#[test]
fn one_chain_of_two_comes_out_alone() {
    let s = Structure::from_pdb(TWO).expect("two residues");
    assert_eq!(s.len(), 2);
    assert_eq!(s.chain('A').map(|c| c.len()), Some(1));
    assert_eq!(s.chain('B').map(|c| c.len()), Some(1));
    assert_eq!(
        s.chain('C').map(|c| c.len()),
        None,
        "an absent chain gave an empty structure"
    );
    assert_eq!(s.chain('A').unwrap().residues()[0].name, "MET");
    println!("  chains A and B separate, chain C is None");
}
