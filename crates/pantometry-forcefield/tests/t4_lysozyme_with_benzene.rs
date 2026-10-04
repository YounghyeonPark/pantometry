//! **T4 lysozyme L99A with benzene, PDB 181L, as one typed system.**
//!
//! Each test holds the system to a fact the crate did not produce: the entry's own records
//! (`ATOM`, `HETATM`, `SEQRES`, `REMARK 465`, `SEQADV`), read here with string operations; the
//! textbook formula of each amino-acid residue; the valence of each element; and the geometry of
//! the dictionary templates the hydrogens came from. Measured numbers that nothing bounds in
//! advance — the whole-residue fits, the closest contacts — are printed and not asserted.

mod protein;

use pantometry_forcefield::pdb::{Histidine, Part, Placement, PEPTIDE_BOND_LIMIT};
use pantometry_forcefield::{qeq, uff, BondOrder, Element, ForceField, UffType};
use protein::*;
use std::collections::BTreeMap;

const ANGSTROM: f64 = 1e-10;
const DEG: f64 = std::f64::consts::PI / 180.0;

/// The residues the file models: 162 of `SEQRES`'s 164.
const MODELLED: usize = 162;

/// The binding pocket's radius for QEq, Å: a residue is in it when any heavy atom is this near a
/// benzene atom. 6 Å is 19 residues and 322 atoms, 8.5 s unoptimised; 8 Å is 44 and 719 atoms
/// and 66 s, too long for a test the gate runs unoptimised, so that comparison is ignored by
/// default (`benzenes_charges_hardly_move_with_the_pockets_edge`).
const POCKET: f64 = 6.0;

/// What was dropped is exactly the file's `HETATM` residues other than benzene, counted here from
/// the records; and 181L has no alternate locations, so none was dropped by that rule.
#[test]
fn the_waters_ions_and_additive_are_dropped_and_counted() {
    let s = system();
    let mut expected: BTreeMap<String, (std::collections::BTreeSet<String>, usize)> =
        BTreeMap::new();
    for l in atom_lines(PDB_181L).filter(|l| l.starts_with("HETATM")) {
        let name = columns(l, 18, 20).trim().to_string();
        if name == "BNZ" {
            continue;
        }
        let e = expected.entry(name).or_default();
        e.0.insert(columns(l, 22, 26).to_string());
        e.1 += 1;
    }
    let got: BTreeMap<String, (usize, usize)> = s
        .dropped()
        .iter()
        .map(|d| (d.name.clone(), (d.residues, d.atoms)))
        .collect();
    let want: BTreeMap<String, (usize, usize)> = expected
        .into_iter()
        .map(|(k, (r, a))| (k, (r.len(), a)))
        .collect();
    assert_eq!(got, want);
    eprintln!("dropped: {got:?}");
    assert_eq!(got["HOH"], (136, 136));
    assert_eq!(got["CL"], (2, 2));
    assert_eq!(got["HED"], (1, 8));
    let alternates = atom_lines(PDB_181L)
        .filter(|l| columns(l, 17, 17) != " ")
        .count();
    assert_eq!(alternates, 0, "181L has no alternate locations");
    assert_eq!(s.alternates_dropped(), 0);
}

/// Every residue has every heavy atom its template has — less `OXT` except at the C-terminus —
/// and every heavy atom in the file is in the system at the file's coordinates, bit for bit.
#[test]
fn every_residue_has_its_heavy_atoms_at_the_files_coordinates() {
    let s = system();
    let c = s.component();
    for r in s.residues() {
        let t = template(&r.name);
        let mut want: Vec<&str> = t
            .atoms()
            .iter()
            .filter(|a| a.element != Element::H)
            .map(|a| a.name.as_str())
            .filter(|n| r.c_terminal || r.part == Part::Ligand || *n != "OXT")
            .collect();
        let mut got: Vec<&str> = r
            .atoms
            .clone()
            .filter(|&i| c.atoms()[i].element != Element::H)
            .map(|i| s.atom_name(i))
            .collect();
        want.sort_unstable();
        got.sort_unstable();
        assert_eq!(got, want, "{}", r.label());
    }
    let mut from_file = 0;
    for l in atom_lines(PDB_181L) {
        let res = columns(l, 18, 20).trim();
        if l.starts_with("HETATM") && res != "BNZ" {
            continue;
        }
        from_file += 1;
        let number: i32 = columns(l, 23, 26).trim().parse().unwrap();
        let name = columns(l, 13, 16).trim();
        let at: Vec<f64> = [31, 39, 47]
            .iter()
            .map(|&k| columns(l, k, k + 7).trim().parse::<f64>().unwrap() * ANGSTROM)
            .collect();
        let i = (0..c.atoms().len())
            .find(|&i| s.residue_of(i).number == number && s.atom_name(i) == name)
            .unwrap_or_else(|| panic!("{res}{number} {name} is not in the system"));
        assert_eq!(c.atoms()[i].at.to_vec(), at, "{res}{number} {name}");
        assert_eq!(s.placements()[i], Placement::Crystal);
    }
    let heavy = c.atoms().iter().filter(|a| a.element != Element::H).count();
    // 1289 ATOM records and benzene's six, and the one OXT placed on Lys162.
    assert_eq!(from_file, 1289 + 6);
    assert_eq!(heavy, from_file + 1);
    let placed: Vec<String> = (0..c.atoms().len())
        .filter(|&i| c.atoms()[i].element != Element::H && s.placements()[i] != Placement::Crystal)
        .map(|i| c.atoms()[i].name.clone())
        .collect();
    assert_eq!(placed, ["A:LYS162:OXT"]);
}

/// **The dictionary's own protonation, which the rules start from.** Every amino-acid template is
/// a free amino acid — a neutral NH₂ (`H`, `H2`; proline's one `H`) and a neutral COOH (`OXT`,
/// `HXT`) — and its side chain is drawn neutral for Asp (`HD2`) and Glu (`HE2`) and charged for
/// Lys, Arg and His, whose `_chem_comp.pdbx_formal_charge` is +1. Benzene is C₆H₆ and neutral.
#[test]
fn the_dictionarys_templates_are_free_amino_acids() {
    for code in CODES {
        let t = template(code);
        let has = |n: &str| t.atoms().iter().any(|a| a.name == n);
        let charge: i32 = t.atoms().iter().map(|a| a.charge).sum();
        assert_eq!(t.formal_charge(), Some(charge), "{code}");
        if code == "BNZ" {
            assert_eq!(charge, 0);
            assert_eq!(t.formula(), Some("C6 H6"));
            continue;
        }
        assert!(has("N") && has("H") && has("OXT") && has("HXT"), "{code}");
        assert_eq!(has("H2"), code != "PRO", "{code}");
        let charged = ["LYS", "ARG", "HIS"].contains(&code);
        assert_eq!(charge, i32::from(charged), "{code}");
    }
    assert!(template("ASP").atoms().iter().any(|a| a.name == "HD2"));
    assert!(template("GLU").atoms().iter().any(|a| a.name == "HE2"));
    let his = template("HIS");
    assert!(["HD1", "HE2"]
        .iter()
        .all(|n| his.atoms().iter().any(|a| a.name == *n)));
}

/// Protein first, benzene last: the twelve atoms of one residue, C₆H₆, and every other atom
/// protein — the partition a binding energy subtracts by.
#[test]
fn the_partition_is_the_protein_and_then_benzene() {
    let s = system();
    let n = s.component().atoms().len();
    let ligand = s.atoms_in(Part::Ligand);
    let protein = s.atoms_in(Part::Protein);
    assert_eq!(ligand, (n - 12..n).collect::<Vec<_>>());
    assert_eq!(protein, (0..n - 12).collect::<Vec<_>>());
    let elements: Vec<Element> = ligand
        .iter()
        .map(|&i| s.component().atoms()[i].element)
        .collect();
    assert_eq!(elements.iter().filter(|e| **e == Element::C).count(), 6);
    assert_eq!(elements.iter().filter(|e| **e == Element::H).count(), 6);
    for &i in &ligand {
        assert_eq!(s.residue_of(i).label(), "A:BNZ400");
    }
    // No bond crosses the partition.
    for b in s.component().bonds() {
        assert_eq!(s.part(b.atoms[0]), s.part(b.atoms[1]));
    }
    eprintln!(
        "{n} atoms: {} protein, {} ligand",
        protein.len(),
        ligand.len()
    );
}

/// The modelled residues are `SEQRES` from the first, in order and numbered from one; the two
/// left over are exactly what `REMARK 465` says was not located.
#[test]
fn the_sequence_is_the_entrys_seqres() {
    let s = system();
    let seq = seqres(PDB_181L, 'A');
    assert_eq!(seq.len(), 164);
    let protein: Vec<_> = s
        .residues()
        .iter()
        .filter(|r| r.part == Part::Protein)
        .collect();
    assert_eq!(protein.len(), MODELLED);
    assert_eq!(s.component().id(), "181L");
    let termini = |f: fn(&&&pantometry_forcefield::Residue) -> bool| -> Vec<i32> {
        protein.iter().filter(f).map(|r| r.number).collect()
    };
    assert_eq!(termini(|r| r.n_terminal), [1]);
    assert_eq!(termini(|r| r.c_terminal), [MODELLED as i32]);
    for (k, r) in protein.iter().enumerate() {
        assert_eq!(r.name, seq[k], "residue {}", k + 1);
        assert_eq!(r.number, k as i32 + 1);
        assert_eq!(r.chain, 'A');
    }
    let missing: Vec<(i32, String)> = PDB_181L
        .lines()
        .filter(|l| l.starts_with("REMARK 465     "))
        .filter(|l| columns(l, 20, 20) == "A")
        .map(|l| {
            (
                columns(l, 22, 26).trim().parse().unwrap(),
                columns(l, 16, 18).to_string(),
            )
        })
        .collect();
    assert_eq!(
        missing,
        [(163, "ASN".to_string()), (164, "LEU".to_string())]
    );
    let unmodelled: Vec<(i32, String)> = s
        .unmodelled()
        .iter()
        .map(|u| (u.number, u.name.clone()))
        .collect();
    assert_eq!(unmodelled, missing);
    assert_eq!(seq[162..], ["ASN", "LEU"]);
}

/// One bond per hydrogen, and every heavy atom's bond orders — Kekulé, as the dictionary writes an
/// aromatic ring — summing to its element's valence shifted by its formal charge: C 4, N 3 + q,
/// O 2 + q, S 2.
#[test]
fn every_atom_has_its_valence_with_its_formal_charge() {
    let s = system();
    let c = s.component();
    for (i, a) in c.atoms().iter().enumerate() {
        let bonds: Vec<_> = c.neighbours(i).collect();
        let sum: i32 = bonds.iter().map(|(_, b)| b.order.count() as i32).sum();
        let want = match a.element {
            Element::H => {
                assert_eq!(bonds.len(), 1, "{} has {} bonds", a.name, bonds.len());
                assert_eq!(bonds[0].1.order, BondOrder::Single, "{}", a.name);
                assert_eq!(a.charge, 0, "{}", a.name);
                continue;
            }
            Element::C => {
                assert_eq!(a.charge, 0, "{}", a.name);
                4
            }
            Element::N => 3 + a.charge,
            Element::O => 2 + a.charge,
            Element::S => {
                assert_eq!(a.charge, 0, "{}", a.name);
                2
            }
            other => panic!(
                "{} is {other}, which T4 lysozyme and benzene do not have",
                a.name
            ),
        };
        assert_eq!(sum, want, "{}: bond orders {sum}, valence {want}", a.name);
    }
}

/// The textbook formula of each residue in a chain — the amino acid less one water — as
/// `(C, H, N, O, S)`, typed from the neutral amino acids' formulas and nothing in this crate.
fn residue_formula(name: &str) -> [usize; 5] {
    match name {
        "GLY" => [2, 3, 1, 1, 0],
        "ALA" => [3, 5, 1, 1, 0],
        "SER" => [3, 5, 1, 2, 0],
        "PRO" => [5, 7, 1, 1, 0],
        "VAL" => [5, 9, 1, 1, 0],
        "THR" => [4, 7, 1, 2, 0],
        "CYS" => [3, 5, 1, 1, 1],
        "LEU" | "ILE" => [6, 11, 1, 1, 0],
        "ASN" => [4, 6, 2, 2, 0],
        "ASP" => [4, 5, 1, 3, 0],
        "GLN" => [5, 8, 2, 2, 0],
        "LYS" => [6, 12, 2, 1, 0],
        "GLU" => [5, 7, 1, 3, 0],
        "MET" => [5, 9, 1, 1, 1],
        "HIS" => [6, 7, 3, 1, 0],
        "PHE" => [9, 9, 1, 1, 0],
        "ARG" => [6, 12, 4, 1, 0],
        "TYR" => [9, 9, 1, 2, 0],
        "TRP" => [11, 10, 2, 1, 0],
        other => panic!("{other} is not an amino acid"),
    }
}

/// Each residue's formula and charge after protonation are what the rule says, worked here from
/// the textbook formula: Asp and Glu lose a proton (−1), Lys and Arg gain one (+1), a salt-bridged
/// His gains one (+1, by [`salt_bridged_histidines`]) and any other is neutral, the N-terminus
/// gains two hydrogens (NH₃⁺, +1) and the C-terminus an oxygen (COO⁻, −1). Benzene is C₆H₆ and
/// neutral.
#[test]
fn each_residue_is_the_formula_its_protonation_says() {
    let s = system();
    let c = s.component();
    for r in s.residues() {
        let mut count = [0usize; 5];
        for i in r.atoms.clone() {
            let k = match c.atoms()[i].element {
                Element::C => 0,
                Element::H => 1,
                Element::N => 2,
                Element::O => 3,
                Element::S => 4,
                other => panic!("{other}"),
            };
            count[k] += 1;
        }
        if r.part == Part::Ligand {
            assert_eq!(r.name, "BNZ");
            assert_eq!(count, [6, 6, 0, 0, 0]);
            assert_eq!(r.charge, 0);
            continue;
        }
        let mut want = residue_formula(&r.name);
        let mut charge = 0;
        match r.name.as_str() {
            "ASP" | "GLU" => {
                want[1] -= 1;
                charge -= 1;
            }
            "LYS" | "ARG" => {
                want[1] += 1;
                charge += 1;
            }
            "HIS" if salt_bridged_histidines().contains(&r.number) => {
                want[1] += 1;
                charge += 1;
            }
            _ => {}
        }
        if r.n_terminal {
            want[1] += 2;
            charge += 1;
        }
        if r.c_terminal {
            want[3] += 1;
            charge -= 1;
        }
        assert_eq!(count, want, "{} (C, H, N, O, S)", r.label());
        assert_eq!(r.charge, charge, "{}", r.label());
    }
}

/// The histidines the stated rule calls salt bridges, worked here from the file's records alone:
/// those with ND1 or NE2 within 3.2 Å of an Asp OD1/OD2 or Glu OE1/OE2.
fn salt_bridged_histidines() -> Vec<i32> {
    let atom = |l: &str| -> (String, i32, String, [f64; 3]) {
        (
            columns(l, 18, 20).to_string(),
            columns(l, 23, 26).trim().parse().unwrap(),
            columns(l, 13, 16).trim().to_string(),
            [31, 39, 47].map(|k| columns(l, k, k + 7).trim().parse::<f64>().unwrap()),
        )
    };
    let records: Vec<_> = atom_lines(PDB_181L)
        .filter(|l| l.starts_with("ATOM  "))
        .map(atom)
        .collect();
    let oxygens: Vec<[f64; 3]> = records
        .iter()
        .filter(|(res, _, name, _)| {
            (res == "ASP" && (name == "OD1" || name == "OD2"))
                || (res == "GLU" && (name == "OE1" || name == "OE2"))
        })
        .map(|r| r.3)
        .collect();
    let mut out: Vec<i32> = records
        .iter()
        .filter(|(res, _, name, at)| {
            res == "HIS"
                && (name == "ND1" || name == "NE2")
                && oxygens.iter().any(|o| len(sub(*at, *o)) <= 3.2)
        })
        .map(|r| r.1)
        .collect();
    out.dedup();
    out
}

/// The total formal charge is the sequence's, counted from `SEQRES` text — Arg + Lys − Asp − Glu,
/// the termini cancelling — plus one for each histidine the rule calls a salt bridge, counted from
/// the file's coordinates. The two unmodelled residues, Asn and Leu, carry no charge, so the
/// modelled chain's total is the whole sequence's.
#[test]
fn the_formal_charge_is_the_sequences() {
    let s = system();
    let seq = seqres(PDB_181L, 'A');
    let n = |name: &str| seq.iter().filter(|x| *x == name).count() as i32;
    let sequence_charge = n("ARG") + n("LYS") - n("ASP") - n("GLU");
    eprintln!(
        "SEQRES: Arg {} Lys {} Asp {} Glu {} His {} Cys {}; charge {sequence_charge:+}",
        n("ARG"),
        n("LYS"),
        n("ASP"),
        n("GLU"),
        n("HIS"),
        n("CYS")
    );
    let bridged = salt_bridged_histidines();
    assert_eq!(bridged, [31]);
    assert_eq!(sequence_charge, 8);
    assert_eq!(s.formal_charge(), sequence_charge + bridged.len() as i32);
    assert_eq!(s.formal_charge(), 9);
    for u in s.unmodelled() {
        assert!(!["ARG", "LYS", "ASP", "GLU", "HIS"].contains(&u.name.as_str()));
    }
    let ligand: i32 = s
        .atoms_in(Part::Ligand)
        .iter()
        .map(|&i| s.component().atoms()[i].charge)
        .sum();
    assert_eq!(ligand, 0);
}

/// T4 lysozyme's one histidine, His31, is a salt bridge by the rule — its Nδ1 is 2.66 Å from
/// Asp70's Oδ2 in the file — so it is doubly protonated, +1. Overridden to either neutral state,
/// the total falls by one, and each keeps the one hydrogen its state names.
#[test]
fn histidine_31_is_a_salt_bridge_by_the_rule() {
    let s = system();
    let his: Vec<_> = s.residues().iter().filter(|r| r.name == "HIS").collect();
    assert_eq!(his.len(), 1);
    assert_eq!(his[0].number, 31);
    assert_eq!(his[0].histidine, Some(Histidine::Both));
    assert_eq!(his[0].charge, 1);
    let names: Vec<&str> = his[0].atoms.clone().map(|i| s.atom_name(i)).collect();
    assert!(
        names.contains(&"HD1") && names.contains(&"HE2"),
        "{names:?}"
    );
    let at = |number: i32, name: &str| -> [f64; 3] {
        let l = atom_lines(PDB_181L)
            .find(|l| is_atom(l, number, name))
            .unwrap();
        [31, 39, 47].map(|k| columns(l, k, k + 7).trim().parse::<f64>().unwrap())
    };
    let d = len(sub(at(31, "ND1"), at(70, "OD2")));
    assert!((d - 2.66).abs() < 0.005, "{d}");
    for (state, kept, gone) in [
        (Histidine::Delta, "HD1", "HE2"),
        (Histidine::Epsilon, "HE2", "HD1"),
    ] {
        let neutral = pantometry_forcefield::System::from_pdb(
            PDB_181L,
            &templates(),
            &selection().histidine('A', 31, state),
        )
        .unwrap();
        assert_eq!(neutral.formal_charge(), s.formal_charge() - 1);
        let r = neutral.residues().iter().find(|r| r.name == "HIS").unwrap();
        let names: Vec<&str> = r.atoms.clone().map(|i| neutral.atom_name(i)).collect();
        assert!(names.contains(&kept) && !names.contains(&gone), "{names:?}");
        // Each ring atom's Kekulé bond orders meet its valence: the δ state moves the double bond.
        let c = neutral.component();
        for i in r.atoms.clone() {
            let a = &c.atoms()[i];
            let sum: i32 = c.neighbours(i).map(|(_, b)| b.order.count() as i32).sum();
            let want = match a.element {
                Element::H => 1,
                Element::C => 4,
                Element::N => 3 + a.charge,
                Element::O => 2 + a.charge,
                other => panic!("{other}"),
            };
            assert_eq!(sum, want, "{state:?}: {}", a.name);
        }
    }
}

/// One peptide bond per consecutive pair, `C(i)–N(i+1)`, residues − 1 of them, and no disulfide:
/// 181L is the cysteine-free pseudo-wild type, C54T and C97A by its `SEQADV`.
#[test]
fn peptide_bonds_join_the_chain_and_there_is_no_disulfide() {
    let s = system();
    let c = s.component();
    assert_eq!(s.peptide_bonds().len(), MODELLED - 1);
    let mut longest: f64 = 0.0;
    for (k, &b) in s.peptide_bonds().iter().enumerate() {
        let [i, j] = c.bonds()[b].atoms;
        assert_eq!(s.atom_name(i), "C");
        assert_eq!(s.atom_name(j), "N");
        assert_eq!(s.residue_of(i).number, k as i32 + 1);
        assert_eq!(s.residue_of(j).number, k as i32 + 2);
        let d = len(sub(c.atoms()[i].at, c.atoms()[j].at)) / ANGSTROM;
        longest = longest.max(d);
    }
    eprintln!("longest peptide C-N: {longest:.3} Å (limit {PEPTIDE_BOND_LIMIT})");
    assert!(s.disulfides().is_empty());
    assert_eq!(
        seqres(PDB_181L, 'A').iter().filter(|x| *x == "CYS").count(),
        0
    );
    let seqadv: Vec<&str> = PDB_181L
        .lines()
        .filter(|l| l.starts_with("SEQADV"))
        .map(|l| columns(l, 13, 49))
        .collect();
    assert!(seqadv
        .iter()
        .any(|l| l.contains("THR A   54") && l.contains("CYS")));
    assert!(seqadv
        .iter()
        .any(|l| l.contains("ALA A   97") && l.contains("CYS")));
}

/// The template atom of the same residue type and name.
fn template_at(code: &str, name: &str) -> [f64; 3] {
    template(code)
        .atoms()
        .iter()
        .find(|a| a.name == name)
        .unwrap_or_else(|| panic!("{code} has no {name}"))
        .at
}

/// **Every placed hydrogen keeps its template's geometry, to a bound the fit earns.**
///
/// A hydrogen carried by a superposition is `P + R (H_t − P_t)`, so its bond length is the
/// template's exactly (to rounding), and so is its angle to any other hydrogen the same fit placed
/// on the same parent. Its angle to a heavy neighbour X differs from the template's by at most the
/// angle between `R (X_t − P_t)` and `X − P` (the triangle inequality on the sphere); both ends of
/// that vector are fitted points, each at most `√k · rmsd` from its crystal position for a k-point
/// fit, so the angle is at most `asin(2 √k rmsd / |X − P|)`. That bound is asserted, atom by atom.
/// The backbone H lies in the `C(i−1)–N–CA` plane at equal angles to both, outside them (the three
/// angles at N make a full turn), at the template's N–H length; the N-terminal H3 likewise.
#[test]
fn every_placed_hydrogen_keeps_its_templates_geometry() {
    let s = system();
    let c = s.component();
    let at = |i: usize| c.atoms()[i].at;
    let (mut worst_angle, mut worst_fraction, mut worst_length) = (0.0f64, 0.0f64, 0.0f64);
    let mut worst_name = String::new();
    let (mut checked_backbone, mut checked_completion) = (0usize, 0usize);
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for (h, a) in c.atoms().iter().enumerate() {
        if a.element != Element::H {
            continue;
        }
        let r = s.residue_of(h);
        let (p, _) = c.neighbours(h).next().unwrap();
        let (hn, pn) = (s.atom_name(h), s.atom_name(p));
        match s.placements()[h] {
            Placement::Template { rmsd, points } => {
                *counts.entry("template").or_default() += 1;
                let ht = sub(template_at(&r.name, hn), template_at(&r.name, pn));
                let hc = sub(at(h), at(p));
                let rel = (len(hc) - len(ht)).abs() / len(ht);
                worst_length = worst_length.max(rel);
                assert!(rel < 1e-12, "{}: length off by {rel:e}", a.name);
                for (x, _) in c.neighbours(p) {
                    if x == h || s.residue_of(x).number != r.number {
                        continue;
                    }
                    let xn = s.atom_name(x);
                    if xn == "H3" {
                        continue;
                    }
                    let xt = sub(template_at(&r.name, xn), template_at(&r.name, pn));
                    let xc = sub(at(x), at(p));
                    let off = (angle_between(hc, xc) - angle_between(ht, xt)).abs();
                    if c.atoms()[x].element == Element::H {
                        assert!(off < 1e-9, "{} to {xn}: {} °", a.name, off / DEG);
                    } else {
                        let bound = (2.0 * (points as f64).sqrt() * rmsd / len(xc)).asin();
                        assert!(
                            off <= bound + 1e-12,
                            "{} to {xn}: {:.3}° against a bound of {:.3}°",
                            a.name,
                            off / DEG,
                            bound / DEG
                        );
                        if off > worst_angle {
                            worst_angle = off;
                            worst_name = format!("{}-{pn}-{xn}", a.name);
                        }
                        worst_fraction = worst_fraction.max(off / bound);
                    }
                }
            }
            Placement::Backbone => {
                *counts.entry("backbone").or_default() += 1;
                assert_eq!(pn, "N");
                let ca = c
                    .neighbours(p)
                    .map(|(x, _)| x)
                    .find(|&x| s.atom_name(x) == "CA")
                    .unwrap();
                let cp = c
                    .neighbours(p)
                    .map(|(x, _)| x)
                    .find(|&x| s.atom_name(x) == "C")
                    .unwrap();
                let (u, v, w) = (sub(at(cp), at(p)), sub(at(ca), at(p)), sub(at(h), at(p)));
                let normal = cross(u, v);
                let out_of_plane = dot(normal, w) / (len(normal) * len(w));
                assert!(out_of_plane.abs() < 1e-12, "{}: {out_of_plane:e}", a.name);
                let (a1, a2) = (angle_between(u, w), angle_between(v, w));
                assert!((a1 - a2).abs() < 1e-12, "{}", a.name);
                // The external bisector: the three angles at a planar N make a full turn.
                let full = a1 + a2 + angle_between(u, v);
                assert!(
                    (full - 2.0 * std::f64::consts::PI).abs() < 1e-12,
                    "{}",
                    a.name
                );
                let nh = len(sub(template_at(&r.name, "H"), template_at(&r.name, "N")));
                assert!((len(w) - nh).abs() < 1e-12 * nh, "{}", a.name);
                checked_backbone += 1;
            }
            Placement::Completion => {
                *counts.entry("completion").or_default() += 1;
                assert_eq!((hn, pn), ("H3", "N"));
                let nh = len(sub(template_at(&r.name, "H"), template_at(&r.name, "N")));
                let w = sub(at(h), at(p));
                assert!((len(w) - nh).abs() < 1e-12 * nh);
                let angles: Vec<f64> = c
                    .neighbours(p)
                    .filter(|&(x, _)| x != h)
                    .map(|(x, _)| angle_between(w, sub(at(x), at(p))) / DEG)
                    .collect();
                eprintln!("N-terminal H3: angles to the other three {angles:.2?}°");
                // Opposite the other three bonds, not among them: at a tetrahedral N every angle
                // is 109.5°, and an H3 placed along their sum instead would be within about 70°
                // of each.
                assert_eq!(angles.len(), 3);
                for angle in &angles {
                    assert!(*angle > 100.0, "{}: {angles:.2?}", a.name);
                }
                checked_completion += 1;
            }
            Placement::Crystal => panic!("{} is a hydrogen read from the file", a.name),
        }
    }
    // The placement rules ran on the atoms they name, counted from `SEQRES`: one backbone H on
    // every residue but the first and the prolines, and one H3 on the one N-terminus.
    let seq = seqres(PDB_181L, 'A');
    let amides = seq[1..MODELLED].iter().filter(|r| *r != "PRO").count();
    assert_eq!(amides, 158);
    assert_eq!(counts.get("backbone"), Some(&amides));
    assert_eq!(checked_backbone, amides);
    assert_eq!(counts.get("completion"), Some(&1));
    assert_eq!(checked_completion, 1);
    eprintln!(
        "hydrogens: {counts:?}; worst length {worst_length:e}, worst angle to a heavy \
         neighbour {:.3}° ({worst_name}), the worst {:.3} of its bound",
        worst_angle / DEG,
        worst_fraction
    );
}

/// `v` turned by `phi` about the unit vector `u` (Rodrigues), written here and not taken from the
/// crate.
fn turned(v: [f64; 3], u: [f64; 3], phi: f64) -> [f64; 3] {
    let (s, c) = phi.sin_cos();
    let along = dot(u, v) * (1.0 - c);
    let x = cross(u, v);
    [
        v[0] * c + x[0] * s + u[0] * along,
        v[1] * c + x[1] * s + u[1] * along,
        v[2] * c + x[2] * s + u[2] * along,
    ]
}

/// **Every rotor was turned to a step at least as clear as the dictionary's.** Each hydrogen
/// [`System::turns`] reports turned is put back, here, by its turn about its parent's one heavy
/// bond; its group's nearest heavy atom (other than the two on the bond) must then be no farther
/// than it is where the rule left it. Every turn is a whole number of 10° steps, inside a third of
/// a revolution for three hydrogens and a whole one otherwise, and one turn per group. **The case
/// the rule was written for**: Lys162's NH₃⁺ at the dictionary's torsion puts a hydrogen within
/// 1.5 Å of Asp159's OD1; turned, every one of them is farther than 1.8 Å, a hydrogen bond's
/// length and not an overlap (measured: 1.923 Å).
///
/// [`System::turns`]: pantometry_forcefield::System::turns
#[test]
fn every_rotor_was_turned_to_a_clearer_step() {
    let s = system();
    let c = s.component();
    let at = |i: usize| c.atoms()[i].at.map(|x| x / ANGSTROM);
    let heavy: Vec<usize> = (0..c.atoms().len())
        .filter(|&i| c.atoms()[i].element != Element::H)
        .collect();
    let mut parents: Vec<usize> = Vec::new();
    for (h, &t) in s.turns().iter().enumerate() {
        if t != 0.0 {
            assert_eq!(c.atoms()[h].element, Element::H, "{}", c.atoms()[h].name);
            let (p, _) = c.neighbours(h).next().unwrap();
            if !parents.contains(&p) {
                parents.push(p);
            }
        }
    }
    let step = 10.0 * DEG;
    let mut turned_hydrogens = 0;
    for &p in &parents {
        let hs: Vec<usize> = c
            .neighbours(p)
            .map(|(x, _)| x)
            .filter(|&x| c.atoms()[x].element == Element::H)
            .collect();
        let heavies: Vec<usize> = c
            .neighbours(p)
            .map(|(x, _)| x)
            .filter(|&x| c.atoms()[x].element != Element::H)
            .collect();
        assert_eq!(heavies.len(), 1, "{}", c.atoms()[p].name);
        let a = heavies[0];
        let turn = s.turns()[hs[0]];
        for &h in &hs {
            assert_eq!(
                s.turns()[h],
                turn,
                "{}: one turn per rotor",
                c.atoms()[h].name
            );
        }
        let k = turn / step;
        assert!(
            (k - k.round()).abs() < 1e-9,
            "{}: {k} steps",
            c.atoms()[p].name
        );
        let span = if hs.len() == 3 { 120.0 } else { 360.0 };
        assert!(
            turn > 0.0 && turn < span * DEG - 1e-9,
            "{}",
            c.atoms()[p].name
        );
        let d = sub(at(p), at(a));
        let axis = [d[0] / len(d), d[1] / len(d), d[2] / len(d)];
        let back: Vec<[f64; 3]> = hs
            .iter()
            .map(|&h| {
                let v = turned(sub(at(h), at(p)), axis, -turn);
                [at(p)[0] + v[0], at(p)[1] + v[1], at(p)[2] + v[2]]
            })
            .collect();
        let now: Vec<[f64; 3]> = hs.iter().map(|&h| at(h)).collect();
        let clearance = |pts: &[[f64; 3]]| {
            pts.iter()
                .flat_map(|q| {
                    heavy
                        .iter()
                        .filter(|&&x| x != p && x != a)
                        .map(move |&x| len(sub(*q, at(x))))
                })
                .fold(f64::INFINITY, f64::min)
        };
        let (chosen, start) = (clearance(&now), clearance(&back));
        assert!(
            chosen >= start - 1e-9,
            "{}: {chosen:.3} Å turned, {start:.3} Å at the dictionary's torsion",
            c.atoms()[p].name
        );
        turned_hydrogens += hs.len();
    }
    assert_eq!(
        turned_hydrogens,
        s.turns().iter().filter(|t| **t != 0.0).count()
    );
    eprintln!(
        "{} rotors turned, {turned_hydrogens} hydrogens",
        parents.len()
    );

    let named = |n: &str| {
        (0..c.atoms().len())
            .find(|&i| c.atoms()[i].name == n)
            .unwrap()
    };
    let (nz, od1) = (named("A:LYS162:NZ"), named("A:ASP159:OD1"));
    assert!(parents.contains(&nz), "Lys162's NH3+ was not turned");
    let ce = named("A:LYS162:CE");
    let d = sub(at(nz), at(ce));
    let axis = [d[0] / len(d), d[1] / len(d), d[2] / len(d)];
    for h in ["A:LYS162:HZ1", "A:LYS162:HZ2", "A:LYS162:HZ3"].map(named) {
        assert!(s.turns()[h] != 0.0);
        let now = len(sub(at(h), at(od1)));
        assert!(now > 1.8, "{}: {now:.3} Å", c.atoms()[h].name);
        eprintln!("{} to Asp159 OD1 turned: {now:.3} Å", c.atoms()[h].name);
    }
    let closest_before = ["A:LYS162:HZ1", "A:LYS162:HZ2", "A:LYS162:HZ3"]
        .map(named)
        .iter()
        .map(|&h| {
            let v = turned(sub(at(h), at(nz)), axis, -s.turns()[h]);
            len(sub(
                [at(nz)[0] + v[0], at(nz)[1] + v[1], at(nz)[2] + v[2]],
                at(od1),
            ))
        })
        .fold(f64::INFINITY, f64::min);
    assert!(closest_before < 1.5, "{closest_before}");
    eprintln!("Lys162 NH3+ to Asp159 OD1: {closest_before:.3} Å at the dictionary's torsion");
}

/// **No hydrogen sits on a heavy atom.** Every hydrogen–heavy-atom pair that is neither bonded nor
/// 1-3 is counted by distance and the closest are printed. None may be within 1.5 Å: the shortest
/// hydrogen-bond H···O contacts in proteins are about 1.6 Å, so anything shorter is two atoms on
/// top of each other.
#[test]
fn no_hydrogen_sits_on_a_heavy_atom() {
    let s = system();
    let c = s.component();
    let n = c.atoms().len();
    let mut near: Vec<std::collections::BTreeSet<usize>> = vec![Default::default(); n];
    for (i, set) in near.iter_mut().enumerate() {
        for (j, _) in c.neighbours(i) {
            set.insert(j);
            for (k, _) in c.neighbours(j) {
                set.insert(k);
            }
        }
    }
    let mut contacts: Vec<(f64, usize, usize)> = Vec::new();
    for h in (0..n).filter(|&i| c.atoms()[i].element == Element::H) {
        for x in (0..n).filter(|&i| c.atoms()[i].element != Element::H) {
            if near[h].contains(&x) {
                continue;
            }
            let d = len(sub(c.atoms()[h].at, c.atoms()[x].at)) / ANGSTROM;
            if d < 2.5 {
                contacts.push((d, h, x));
            }
        }
    }
    contacts.sort_by(|a, b| a.0.total_cmp(&b.0));
    for limit in [1.5, 1.6, 1.8, 2.0, 2.2, 2.5] {
        let k = contacts.iter().filter(|c| c.0 < limit).count();
        eprintln!("H-heavy contacts below {limit} Å: {k}");
    }
    for (d, h, x) in contacts.iter().take(12) {
        eprintln!(
            "  {d:.3} Å  {} ({})  {}",
            c.atoms()[*h].name,
            s.placements()[*h],
            c.atoms()[*x].name
        );
    }
    let clashes: Vec<String> = contacts
        .iter()
        .filter(|c| c.0 < 1.5)
        .map(|&(d, h, x)| format!("{} {} {d:.3}", c.atoms()[h].name, c.atoms()[x].name))
        .collect();
    assert!(clashes.is_empty(), "{clashes:?}");
}

/// The RMS distance after the rigid motion that puts `t[0]` on `c[0]`, the direction to `t[1]`
/// on the direction to `c[1]`, and the plane of `t[0..3]` on that of `c[0..3]` — one particular
/// superposition, built from three atoms, so the least-squares one can only be better.
fn frame_rmsd(t: &[[f64; 3]], c: &[[f64; 3]]) -> f64 {
    let frame = |p: &[[f64; 3]]| {
        let e1 = sub(p[1], p[0]);
        let e1 = [e1[0] / len(e1), e1[1] / len(e1), e1[2] / len(e1)];
        let v = sub(p[2], p[0]);
        let along = dot(v, e1);
        let e2 = sub(v, [e1[0] * along, e1[1] * along, e1[2] * along]);
        let e2 = [e2[0] / len(e2), e2[1] / len(e2), e2[2] / len(e2)];
        [e1, e2, cross(e1, e2)]
    };
    let (ft, fc) = (frame(t), frame(c));
    let sq: f64 = t
        .iter()
        .zip(c)
        .map(|(tp, cp)| {
            let d = sub(*tp, t[0]);
            let local = [dot(d, ft[0]), dot(d, ft[1]), dot(d, ft[2])];
            let mut moved = c[0];
            for k in 0..3 {
                for (m, f) in moved.iter_mut().zip(fc[k]) {
                    *m += local[k] * f;
                }
            }
            dot(sub(moved, *cp), sub(moved, *cp))
        })
        .sum();
    (sq / t.len() as f64).sqrt()
}

/// A lower bound on the least-squares RMSD of `t` onto `c`, from distances alone. A rigid motion
/// keeps `|t_i − t_j|`, so with `e_i` the fitted `t_i`'s miss of `c_i`,
/// `| |t_i − t_j| − |c_i − c_j| | ≤ |e_i − e_j| ≤ |e_i| + |e_j| ≤ √2 √(|e_i|² + |e_j|²) ≤ √(2k) rmsd`,
/// the last because `Σ |e|² = k rmsd²`. So `rmsd ≥ max |Δd_ij| / √(2k)` for every pair — no fit
/// can do better, however it is computed, and a fit that reported less would be wrong.
fn distance_floor(t: &[[f64; 3]], c: &[[f64; 3]]) -> f64 {
    let mut worst: f64 = 0.0;
    for i in 0..t.len() {
        for j in i + 1..t.len() {
            worst = worst.max((len(sub(t[i], t[j])) - len(sub(c[i], c[j]))).abs());
        }
    }
    worst / (2.0 * t.len() as f64).sqrt()
}

/// **Every superposition is at least as good as one built by hand, and the hydrogens' never span
/// a torsion.** The least-squares fit cannot be worse than any other rigid motion, so each reported
/// residual is held below the residual of a three-atom frame superposition computed here: for the
/// whole residue on N, CA and C (benzene on C1, C2 and C3), and for each hydrogen's fragment —
/// rebuilt here by the stated rule, the parent and its heavy neighbours, then theirs when fewer
/// than three, and required to have the size the fit reports — on its first three atoms. The whole
/// residue's residual says how far the crystal side chain is from the dictionary conformer, which
/// nothing bounds above, so its size is printed; the fragments' are bond lengths and angles alone.
/// **And no better than distances allow**: every reported residual, whole or fragment, is held
/// above [`distance_floor`], so a residual reported too small is caught as well as one too large.
#[test]
fn every_superposition_is_no_worse_than_a_frame_fit() {
    let s = system();
    let c = s.component();
    let angstrom = |p: [f64; 3]| p.map(|x| x / ANGSTROM);
    let mut tightest = f64::INFINITY;
    let mut whole: Vec<(f64, f64, String)> = Vec::new();
    for r in s.residues() {
        let t = template(&r.name);
        let mut names: Vec<&str> = vec![];
        let first: &[&str] = if r.part == Part::Ligand {
            &["C1", "C2", "C3"]
        } else {
            &["N", "CA", "C"]
        };
        names.extend(first);
        for a in t.atoms() {
            if a.element != Element::H && !names.contains(&a.name.as_str()) {
                names.push(a.name.as_str());
            }
        }
        let crystal = |name: &str| {
            r.atoms
                .clone()
                .find(|&i| s.atom_name(i) == name && s.placements()[i] == Placement::Crystal)
        };
        names.retain(|n| crystal(n).is_some());
        let tp: Vec<[f64; 3]> = names
            .iter()
            .map(|n| angstrom(template_at(&r.name, n)))
            .collect();
        let cp: Vec<[f64; 3]> = names
            .iter()
            .map(|n| angstrom(c.atoms()[crystal(n).unwrap()].at))
            .collect();
        let by_hand = frame_rmsd(&tp, &cp);
        let fit = r.fit_rmsd / ANGSTROM;
        let floor = distance_floor(&tp, &cp);
        assert!(fit >= floor - 1e-12, "{}: {fit} below {floor}", r.label());
        tightest = tightest.min(fit / floor);
        assert!(
            fit <= by_hand + 1e-12,
            "{}: {fit} against {by_hand}",
            r.label()
        );
        whole.push((fit, by_hand, r.label()));
    }
    whole.sort_by(|a, b| b.0.total_cmp(&a.0));
    let mean = whole.iter().map(|w| w.0).sum::<f64>() / whole.len() as f64;
    eprintln!(
        "whole-residue RMSD: mean {mean:.3} Å; largest (fit, frame) {:.3?}",
        &whole[..6]
    );
    eprintln!("whole-residue fit / its distance floor: at least {tightest:.2}");
    for code in CODES {
        let v: Vec<f64> = whole
            .iter()
            .filter(|w| w.2.contains(code))
            .map(|w| w.0)
            .collect();
        if !v.is_empty() {
            eprintln!(
                "  {code} x{}: max {:.3} Å",
                v.len(),
                v.iter().copied().fold(0.0, f64::max)
            );
        }
    }

    let mut local: Vec<(f64, f64, String)> = Vec::new();
    for (h, p) in s.placements().iter().enumerate() {
        let Placement::Template { rmsd, points } = *p else {
            continue;
        };
        let r = s.residue_of(h);
        let t = template(&r.name);
        let kept = |name: &str| r.atoms.clone().any(|i| s.atom_name(i) == name);
        let name_of = |k: usize| t.atoms()[k].name.as_str();
        let heavy_kept = |k: usize| t.atoms()[k].element != Element::H && kept(name_of(k));
        let me = t
            .atoms()
            .iter()
            .position(|a| a.name == s.atom_name(h))
            .unwrap();
        let parent = t
            .neighbours(me)
            .map(|(k, _)| k)
            .find(|&k| heavy_kept(k))
            .unwrap();
        let mut fragment = vec![parent];
        fragment.extend(
            t.neighbours(parent)
                .map(|(k, _)| k)
                .filter(|&k| k != me && heavy_kept(k)),
        );
        if fragment.len() < 3 {
            let shell: Vec<usize> = fragment[1..].to_vec();
            for j in shell {
                for k in t.neighbours(j).map(|(k, _)| k) {
                    if k != me && heavy_kept(k) && !fragment.contains(&k) {
                        fragment.push(k);
                    }
                }
            }
        }
        assert_eq!(fragment.len(), points, "{}", c.atoms()[h].name);
        let tp: Vec<[f64; 3]> = fragment
            .iter()
            .map(|&k| angstrom(t.atoms()[k].at))
            .collect();
        let cp: Vec<[f64; 3]> = fragment
            .iter()
            .map(|&k| {
                let i = r
                    .atoms
                    .clone()
                    .find(|&i| s.atom_name(i) == name_of(k))
                    .unwrap();
                angstrom(c.atoms()[i].at)
            })
            .collect();
        let by_hand = frame_rmsd(&tp, &cp);
        let fit = rmsd / ANGSTROM;
        let floor = distance_floor(&tp, &cp);
        assert!(
            fit >= floor - 1e-12,
            "{}: {fit} below {floor}",
            c.atoms()[h].name
        );
        assert!(
            fit <= by_hand + 1e-12,
            "{}: {fit} against {by_hand}",
            c.atoms()[h].name
        );
        local.push((fit, by_hand, c.atoms()[h].name.clone()));
    }
    local.sort_by(|a, b| b.0.total_cmp(&a.0));
    let mean = local.iter().map(|w| w.0).sum::<f64>() / local.len() as f64;
    eprintln!(
        "local fits: {} hydrogens, mean {mean:.4} Å; largest (fit, frame) {:.4?}",
        local.len(),
        &local[..4]
    );
    let turned = s.turns().iter().filter(|t| **t != 0.0).count();
    eprintln!("rotor hydrogens turned from the dictionary's torsion: {turned}");
}

/// Every atom typed, and the force field built for all of them: no hypervalent sulfur (the five
/// methionines' SD are divalent), and every type one UFF has.
#[test]
fn every_atom_is_typed_and_the_force_field_accepts_the_system() {
    let s = system();
    let c = s.component();
    let types = uff::assign(c);
    assert_eq!(types.len(), c.atoms().len());
    let mut histogram: BTreeMap<String, usize> = BTreeMap::new();
    for t in &types {
        *histogram.entry(t.label().to_string()).or_default() += 1;
    }
    eprintln!("{} atoms: {histogram:?}", c.atoms().len());
    for (i, t) in types.iter().enumerate() {
        if c.atoms()[i].element == Element::S {
            assert_eq!(*t, UffType::S3Divalent, "{}", c.atoms()[i].name);
        }
    }
    // The peptide N and C are the amide's N_R and C_R; benzene is C_R throughout.
    for (i, t) in types.iter().enumerate() {
        let r = s.residue_of(i);
        match (r.part, s.atom_name(i)) {
            (Part::Protein, "N") if !r.n_terminal => assert_eq!(*t, UffType::NR, "{}", r.label()),
            (Part::Protein, "C") if !r.c_terminal => assert_eq!(*t, UffType::CR, "{}", r.label()),
            (Part::Ligand, _) if c.atoms()[i].element == Element::C => {
                assert_eq!(*t, UffType::CR)
            }
            _ => {}
        }
    }
    let ff = ForceField::new(c, &types);
    assert!(ff.is_ok(), "{:?}", ff.err());
    let ff = ff.unwrap();
    eprintln!(
        "force field: {} stretches, {} bends, {} torsions, {} inversions, {} pairs",
        ff.stretches().len(),
        ff.bends().len(),
        ff.torsions().len(),
        ff.inversions().len(),
        ff.pairs().len()
    );
}

/// Benzene and every residue with a heavy atom within `radius` Å of it, whole residues so the
/// total is an integer, equilibrated by QEq: the residues, the atom count, the total, the
/// charges, and benzene's charges in its atom order.
fn pocket(
    s: &pantometry_forcefield::System,
    radius: f64,
) -> (Vec<String>, usize, i32, qeq::Charges, Vec<f64>) {
    let c = s.component();
    let ligand = s.atoms_in(Part::Ligand);
    let atoms: Vec<usize> = s
        .residues()
        .iter()
        .filter(|r| {
            r.part == Part::Ligand
                || r.atoms.clone().any(|i| {
                    c.atoms()[i].element != Element::H
                        && ligand.iter().any(|&l| {
                            len(sub(c.atoms()[i].at, c.atoms()[l].at)) < radius * ANGSTROM
                        })
                })
        })
        .flat_map(|r| r.atoms.clone())
        .collect();
    let residues: Vec<String> = s
        .residues()
        .iter()
        .filter(|r| atoms.contains(&r.atoms.start))
        .map(|r| r.label())
        .collect();
    let elements: Vec<Element> = atoms.iter().map(|&i| c.atoms()[i].element).collect();
    let at: Vec<[f64; 3]> = atoms.iter().map(|&i| c.atoms()[i].at).collect();
    let total: i32 = atoms.iter().map(|&i| c.atoms()[i].charge).sum();
    let q = qeq::Qeq::default()
        .equilibrate(&elements, &at, f64::from(total))
        .unwrap_or_else(|e| panic!("{radius} Å: {e}"));
    let sum: f64 = q.charges.iter().sum();
    assert!(
        (sum - f64::from(total)).abs() < 1e-9,
        "{sum} against {total}"
    );
    let on_ligand: Vec<f64> = atoms
        .iter()
        .zip(&q.charges)
        .filter(|(i, _)| s.part(**i) == Part::Ligand)
        .map(|(_, q)| *q)
        .collect();
    eprintln!(
        "{radius} Å pocket: {} residues, {} atoms, total {total:+}; QEq {} solves, {} at a bound; \
         benzene's charges {on_ligand:.3?}; residues {residues:?}",
        residues.len(),
        atoms.len(),
        q.iterations,
        q.at_bound.len()
    );
    (residues, atoms.len(), total, q, on_ligand)
}

/// **QEq on the binding pocket**: benzene and every residue with a heavy atom within [`POCKET`] of
/// it. Each solve is dense, O(N³): scaling the pocket's 8.5 s by (2616 / 322)³ puts the whole
/// protein past an hour unoptimised, so the pocket is what is solved here. Asserted: its size and
/// total, counted for this entry — 18 residues and benzene, 322 atoms, +1 — that His31, the one
/// charged histidine, is outside it, that no charge reached an end of its range, and that the
/// iteration converged. **Open**: it takes 31 solves, against the six to ten the QEq paper reports
/// (p. 3361) for its molecules; it converges, and why it is slower here is not yet known.
#[test]
fn qeq_converges_on_the_binding_pocket() {
    let s = system();
    let (residues, atoms, total, q, on_ligand) = pocket(&s, POCKET);
    assert_eq!((residues.len(), atoms, total), (19, 322, 1));
    assert!(!residues.iter().any(|r| r == "A:HIS31"), "{residues:?}");
    assert!(q.at_bound.is_empty(), "{:?}", q.at_bound);
    assert!(q.iterations < qeq::Qeq::default().max_iterations);
    assert_eq!(on_ligand.len(), 12);
}

/// **The pocket's edge barely reaches benzene**: 8 Å instead of 6 — 43 residues and benzene, 719
/// atoms — moves no benzene charge by more than 0.004 e (measured 0.0023 e). Ignored: it takes
/// 66 s unoptimised and the gate runs every test that way; run it with
/// `cargo test -p pantometry-forcefield --release --test t4_lysozyme_with_benzene -- --ignored`.
#[test]
#[ignore = "66 s unoptimised; run with --release -- --ignored"]
fn benzenes_charges_hardly_move_with_the_pockets_edge() {
    let s = system();
    let (_, _, _, _, near) = pocket(&s, POCKET);
    let (residues, atoms, total, q, far) = pocket(&s, 8.0);
    assert_eq!((residues.len(), atoms, total), (44, 719, 3));
    assert!(q.at_bound.is_empty(), "{:?}", q.at_bound);
    let moved = near
        .iter()
        .zip(&far)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    eprintln!("benzene's largest charge change, 6 to 8 Å: {moved:.5} e");
    assert!(moved <= 0.004, "{moved}");
}
