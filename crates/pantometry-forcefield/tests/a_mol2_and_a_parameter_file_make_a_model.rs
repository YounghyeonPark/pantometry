//! **A `.mol2` and an AMBER parameter file make a non-bonded model, and every wrong file is refused
//! by name.** The library's GAFF/AM1-BCC path ([`pantometry_forcefield::mol2`],
//! [`pantometry_forcefield::amber`]) with no GAFF value and no FreeSolv file: every `.mol2` and
//! every parameter line here is written for this file, so it ships with the crate.
//!
//! - **The reader** gives each column to the bit, and refuses each malformed file with the
//!   variant that names it — not merely an error, because a reader that refused everything would
//!   pass a test that only asked for one.
//! - **The map is by bonding, and its counts are closed forms**: a graph's automorphisms are 24
//!   for methane (4!), 72 for ethane (2 · 3! · 3!) and 12 for benzene (the hexagon's dihedral
//!   group). A benzene written out of order puts every hydrogen's charge on the hydrogen of the
//!   right carbon. A benzene with one bond moved, the same atoms and bond count, has none.
//! - **The parameter file's `MOD4 RE` block** is read and nothing else, `R*/2` and `ε` become UFF's
//!   `[x, D]` at `x = 2 R*/2`, and each malformed block is refused by name.
//! - **A model applied** is its values on the solute and nothing else changed: the waters' van der
//!   Waals and charges are the field's bits, the rules are Lorentz–Berthelot and `[½, 1/1.2]`, and
//!   benzene's 1-4 pairs are its 21.
//!
//! Every edit checks its anchor occurs exactly once first, as `a_wrong_file_is_refused_by_name.rs`
//! does: an edit whose anchor did not match leaves the file as it was, and the test then fails
//! for the wrong reason — or passes.

use pantometry_forcefield::mol2::{element_of, isomorphisms};
use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::{
    Component, Element, EwaldParameters, Mixing, ModelError, Mol2, Mol2Error, NonbondedModel,
    ParmError, PeriodicBox, PeriodicForceField, VanDerWaals, AMBER_ONE_FOUR,
};

const BNZ: &str = include_str!("../components/BNZ.cif");
const ANGSTROM: f64 = 1e-10;

/// A water, written for this file: three atoms, two bonds, the second bond written backwards.
const WATER: &str =
    "@<TRIPOS>MOLECULE\nwater\n    3     2     1     0     0\nSMALL\nUSER_CHARGES\n\n\n\
     @<TRIPOS>ATOM\n\
     \x20     1 O1          0.0000    0.0000    0.0000 ow        1 WAT     -0.834000\n\
     \x20     2 H1          0.9572    0.0000    0.0000 hw        1 WAT      0.417000\n\
     \x20     3 H2         -0.2400    0.9266    0.0000 hw        1 WAT      0.417000\n\
     @<TRIPOS>BOND\n\
     \x20    1    1    2 1\n\
     \x20    2    3    1 1\n\
     @<TRIPOS>SUBSTRUCTURE\n\
     \x20    1 WAT         1 TEMP              0 ****  ****    0 ROOT\n";

/// A benzene written out of order, for this file: carbons and hydrogens interleaved, the bonds
/// shuffled, each carbon's charge its own (−0.1301 to −0.1306 by its name) and each hydrogen's
/// minus its carbon's. Mapped by file order, the dictionary's C1 would take a hydrogen's charge.
const SHUFFLED_BENZENE: &str =
    "@<TRIPOS>MOLECULE\nshuffled\n   12    12     1     0     0\nSMALL\n\
     USER_CHARGES\n\n\n@<TRIPOS>ATOM\n\
     \x20 1 H4  0 0 0 ha 1 MOL  0.1304\n\
     \x20 2 C2  0 0 0 ca 1 MOL -0.1302\n\
     \x20 3 H1  0 0 0 ha 1 MOL  0.1301\n\
     \x20 4 C5  0 0 0 ca 1 MOL -0.1305\n\
     \x20 5 C1  0 0 0 ca 1 MOL -0.1301\n\
     \x20 6 H6  0 0 0 ha 1 MOL  0.1306\n\
     \x20 7 C4  0 0 0 ca 1 MOL -0.1304\n\
     \x20 8 H3  0 0 0 ha 1 MOL  0.1303\n\
     \x20 9 C6  0 0 0 ca 1 MOL -0.1306\n\
     \x2010 H2  0 0 0 ha 1 MOL  0.1302\n\
     \x2011 C3  0 0 0 ca 1 MOL -0.1303\n\
     \x2012 H5  0 0 0 ha 1 MOL  0.1305\n\
     @<TRIPOS>BOND\n\
     \x20 1  7  1 1\n\x20 2  2 11 ar\n\x20 3  9  5 ar\n\x20 4  5  3 1\n\
     \x20 5 11  7 ar\n\x20 6  4 12 1\n\x20 7  5  2 ar\n\x20 8  9  6 1\n\
     \x20 9  2 10 1\n\x2010  4  9 ar\n\x2011 11  8 1\n\x2012  7  4 ar\n";

/// A parameter file written for this file, in `gaff.dat`'s layout: a title, a mass section whose
/// `ca` line must not be read as van der Waals, a blank line, then the `MOD4 RE` block and `END`.
/// Every value is invented, the mass lines' too: `gaff.dat`'s own mass lines are not copied here.
const PARM: &str = "INVENTED PARAMETERS FOR A TEST\n\
     ca 99.99         9.999               a mass line, not a van der Waals one\n\
     ha 9.999         9.999\n\
     \n\
     MOD4      RE\n\
     \x20 ca          1.9000  0.0900             invented, with a comment of words\n\
     \x20 ha          1.4500  0.0100             invented\n\
     \x20 hz          0.0000  0.0000             no van der Waals at all\n\
     \n\
     END\n";

fn edit(text: &str, from: &str, to: &str) -> String {
    assert_eq!(
        text.matches(from).count(),
        1,
        "anchor {from:?} must occur once"
    );
    text.replacen(from, to, 1)
}

fn graph(m: &Mol2) -> (Vec<Element>, Vec<[usize; 2]>) {
    (m.atoms.iter().map(|a| a.element).collect(), m.bonds.clone())
}

/// **Each column of a `.mol2`, to the bit**, and each malformed file refused by the variant that
/// names it.
#[test]
fn a_mol2_is_read_with_its_columns_and_refused_by_name() {
    let m = Mol2::parse(WATER).expect("a water");
    assert_eq!(m.name, "water");
    let names: Vec<&str> = m.atoms.iter().map(|a| a.name.as_str()).collect();
    assert_eq!(names, ["O1", "H1", "H2"]);
    let elements: Vec<Element> = m.atoms.iter().map(|a| a.element).collect();
    assert_eq!(elements, [Element::O, Element::H, Element::H]);
    let kinds: Vec<&str> = m.atoms.iter().map(|a| a.kind.as_str()).collect();
    assert_eq!(kinds, ["ow", "hw", "hw"]);
    let charges: Vec<u64> = m.atoms.iter().map(|a| a.charge.to_bits()).collect();
    assert_eq!(charges, [-0.834f64, 0.417, 0.417].map(f64::to_bits));
    assert_eq!(
        m.atoms[2].at.map(f64::to_bits),
        [-0.24f64, 0.9266, 0.0].map(f64::to_bits)
    );
    assert_eq!(m.bonds, [[0, 1], [2, 0]]);
    // The element by each of its three readings: the whole letters, their first two, their first.
    for (name, want) in [
        ("Cl12", Element::Cl),
        ("CL1", Element::Cl),
        ("c1", Element::C),
        ("Br", Element::Br),
        ("CLA1", Element::Cl),
        ("HO", Element::H),
        ("CA", Element::C),
        ("OAB", Element::O),
        ("1H2", Element::H),
        ("NA1", Element::N),
    ] {
        assert_eq!(element_of(name), Ok(want), "{name}");
    }

    let refused = |text: String| Mol2::parse(&text).expect_err("this file should be refused");
    let e = refused(edit(WATER, "WAT      0.417000\n      3", "\n      3"));
    assert!(
        matches!(
            e,
            Mol2Error::Columns {
                row: 2,
                columns: 7,
                ..
            }
        ),
        "{e:?}"
    );
    let e = refused(edit(WATER, "-0.834000", "-0.83x"));
    assert!(
        matches!(&e, Mol2Error::NotFinite { what, .. } if what == "atom O1's charge"),
        "{e:?}"
    );
    let e = refused(edit(WATER, "0.9572", "inf"));
    assert!(
        matches!(&e, Mol2Error::NotFinite { what, .. } if what == "a coordinate"),
        "{e:?}"
    );
    let e = refused(edit(WATER, "    3     2  ", "    4     2  "));
    assert_eq!(
        e,
        Mol2Error::Records {
            what: "atoms",
            found: 3,
            counted: 4
        }
    );
    let e = refused(edit(WATER, "    3     2  ", "    3     3  "));
    assert_eq!(
        e,
        Mol2Error::Records {
            what: "bonds",
            found: 2,
            counted: 3
        }
    );
    let e = refused(edit(WATER, "    2    3    1 1", "    2    3"));
    assert!(matches!(e, Mol2Error::ShortBond { row: 2, .. }), "{e:?}");
    let e = refused(edit(WATER, "    3     2  ", "    3x    2  "));
    assert_eq!(e, Mol2Error::Count { text: "3x".into() });
    let e = refused(edit(
        WATER,
        "water\n    3     2     1     0     0\n",
        "water\n    3\n",
    ));
    assert_eq!(e, Mol2Error::NoCounts);
    let e = refused(edit(WATER, "    2    3    1 1", "    2    4    1 1"));
    assert!(matches!(e, Mol2Error::NoSuchAtom { bond: 2, .. }), "{e:?}");
    let e = refused(edit(WATER, "    2    3    1 1", "    2    2    1 1"));
    assert_eq!(
        e,
        Mol2Error::RepeatedBond {
            bond: 2,
            atoms: [2, 1]
        }
    );
    let e = refused(edit(WATER, "    2    3    1 1", "    2    1    1 1"));
    assert_eq!(e, Mol2Error::SelfBond { bond: 2, atom: 1 });
    let e = refused(edit(WATER, "      2 H1", "      7 H1"));
    assert!(
        matches!(&e, Mol2Error::Numbered { what: "atom", row: 2, number } if number == "7"),
        "{e:?}"
    );
    let e = refused(edit(WATER, "    2    3    1 1", "    5    3    1 1"));
    assert!(
        matches!(&e, Mol2Error::Numbered { what: "bond", row: 2, number } if number == "5"),
        "{e:?}"
    );
    let e = refused(edit(WATER, "O1  ", "Xx1 "));
    assert!(
        matches!(&e, Mol2Error::Element { name } if name == "Xx1"),
        "{e:?}"
    );
    assert!(e.to_string().contains("Xx1"), "{e}");
    assert_eq!(
        refused(format!("{WATER}{WATER}")),
        Mol2Error::SecondMolecule
    );
    assert_eq!(
        refused(edit(WATER, "@<TRIPOS>MOLECULE\n", "")),
        Mol2Error::NoMolecule
    );
}

/// **A graph's automorphisms are closed forms**: methane 4! = 24, ethane 2 · 3! · 3! = 72,
/// benzene's hexagon 12. Each graph is a `.mol2` matched onto itself.
#[test]
fn the_automorphisms_are_the_symmetry_groups_orders() {
    let molecule = |name: &str, atoms: &[&str], bonds: &[[usize; 2]]| -> String {
        let mut s = format!(
            "@<TRIPOS>MOLECULE\n{name}\n {} {} 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n",
            atoms.len(),
            bonds.len()
        );
        for (k, a) in atoms.iter().enumerate() {
            s += &format!("{} {a} 0 0 0 x 1 MOL 0.0\n", k + 1);
        }
        s += "@<TRIPOS>BOND\n";
        for (k, b) in bonds.iter().enumerate() {
            s += &format!("{} {} {} 1\n", k + 1, b[0] + 1, b[1] + 1);
        }
        s
    };
    let methane = molecule(
        "methane",
        &["C1", "H1", "H2", "H3", "H4"],
        &[[0, 1], [0, 2], [0, 3], [0, 4]],
    );
    let ethane = molecule(
        "ethane",
        &["C1", "C2", "H1", "H2", "H3", "H4", "H5", "H6"],
        &[[0, 1], [0, 2], [0, 3], [0, 4], [1, 5], [1, 6], [1, 7]],
    );
    let benzene = Mol2::parse(SHUFFLED_BENZENE).expect("the shuffled benzene");
    for (name, m, want) in [
        ("methane", Mol2::parse(&methane).expect("methane"), 24),
        ("ethane", Mol2::parse(&ethane).expect("ethane"), 72),
        ("benzene", benzene, 12),
    ] {
        let g = graph(&m);
        let all = isomorphisms((&g.0, &g.1), (&g.0, &g.1));
        assert_eq!(all.len(), want, "{name}");
        let mut distinct = all.clone();
        distinct.sort();
        distinct.dedup();
        assert_eq!(distinct.len(), want, "{name}: every match once");
        for map in &all {
            for &[a, b] in &g.1 {
                let [x, y] = [map[a], map[b]];
                assert!(
                    g.1.contains(&[x, y]) || g.1.contains(&[y, x]),
                    "{name}: a match that is not a bond's"
                );
            }
        }
    }
}

/// **The map is by bonding, not by file order**: [`SHUFFLED_BENZENE`] onto the dictionary's `BNZ`
/// puts on every hydrogen minus its own carbon's charge, every carbon a carbon's, every type by
/// element, and finds the 12 matches with a spread of the carbons' range, 0.0005 to its rounding.
/// A benzene with one ring bond moved across the ring, the same atoms and bond count, is refused.
#[test]
fn the_map_is_by_bonding_not_by_file_order() {
    let bnz = Component::from_ccd(BNZ).expect("BNZ");
    let file = Mol2::parse(SHUFFLED_BENZENE).expect("the shuffled benzene");
    let m = file.onto(&bnz).expect("benzene onto benzene");
    assert_eq!(m.automorphisms, 12);
    let carbons = [-0.1301, -0.1302, -0.1303, -0.1304, -0.1305, -0.1306];
    for (i, a) in bnz.atoms().iter().enumerate() {
        let q = m.charges[i];
        match a.element {
            Element::C => {
                assert!(carbons.contains(&q), "{}: {q}", a.name);
                assert_eq!(m.kinds[i], "ca", "{}", a.name);
            }
            Element::H => {
                let (c, _) = bnz.neighbours(i).next().expect("a hydrogen's carbon");
                assert_eq!(q, -m.charges[c], "{}", a.name);
                assert_eq!(m.kinds[i], "ha", "{}", a.name);
            }
            other => panic!("{other:?}"),
        }
        assert_eq!(file.atoms[m.map[i]].charge, q, "{}", a.name);
    }
    for b in bnz.bonds() {
        let [x, y] = b.atoms.map(|i| m.map[i]);
        assert!(
            file.bonds.contains(&[x, y]) || file.bonds.contains(&[y, x]),
            "{:?}",
            b.atoms
        );
    }
    // The first match, by the search's order. The dictionary's C1 (atom 0) tries the file's atoms
    // in turn and takes the first carbon, C2 (index 1, −0.1302). Its neighbour C2 (atom 1) then
    // takes the first carbon bonded to that one, C1 (index 4, −0.1301). The last match would put
    // C3 on atom 0, so every assertion above, true of all twelve, cannot tell them apart.
    assert_eq!(&m.map[..2], &[1, 4], "the first match the search finds");
    assert_eq!([m.charges[0], m.charges[1]], [-0.1302, -0.1301]);
    let g = graph(&file);
    let bnz_graph = pantometry_forcefield::mol2::graph(&bnz);
    let all = isomorphisms((&bnz_graph.0, &bnz_graph.1), (&g.0, &g.1));
    assert_eq!(all.len(), 12);
    assert_eq!(m.map, all[0]);
    let range = -0.1301f64 - -0.1306;
    assert!(
        (m.spread - range).abs() <= 4.0 * f64::EPSILON,
        "{} against {range}",
        m.spread
    );

    // C1–C6 (5–9) becomes C1–C4 (5–7): still twelve bonds, and now a ring of four.
    let moved = Mol2::parse(&edit(
        SHUFFLED_BENZENE,
        "\x20 3  9  5 ar",
        "\x20 3  7  5 ar",
    ))
    .expect("still a well-formed file");
    let e = moved.onto(&bnz).expect_err("not benzene");
    assert_eq!(
        e,
        Mol2Error::NotThisMolecule {
            name: "shuffled".into(),
            file: "C6H6".into(),
            target: "C6H6".into(),
        }
    );
    assert!(e.to_string().contains("bonded differently"), "{e}");
    // A hydrogen made a chlorine: the compositions say what differs.
    let e = Mol2::parse(&edit(SHUFFLED_BENZENE, " 6 H6 ", " 6 Cl6"))
        .expect("a chlorobenzene")
        .onto(&bnz)
        .expect_err("not benzene");
    assert!(
        e.to_string()
            .contains("the file is C6H5Cl and the molecule C6H6"),
        "{e}"
    );
}

/// **The `MOD4 RE` block is read and nothing else**, `[x, D]` is `[2 R*/2, ε]` in metres and
/// joules, and each malformed block is refused by the variant that names it.
#[test]
fn a_parameter_files_van_der_waals_are_read_and_refused_by_name() {
    let v = VanDerWaals::from_parm(PARM).expect("the invented file");
    let kinds: Vec<&str> = v.entries().iter().map(|e| e.0.as_str()).collect();
    assert_eq!(kinds, ["ca", "ha", "hz"], "the mass lines are not read");
    assert_eq!(v.get("ca"), Some((1.9, 0.09)));
    assert_eq!(v.get("ha"), Some((1.45, 0.01)));
    assert_eq!(v.get("hz"), Some((0.0, 0.0)));
    assert_eq!(v.get("c3"), None);
    let [x, d] = v.uff_form("ca").expect("ca");
    assert!((x - 3.8 * ANGSTROM).abs() <= 2.0 * f64::EPSILON * x, "{x}");
    assert!(
        (d - 0.09 * KCAL_PER_MOL).abs() <= 2.0 * f64::EPSILON * d,
        "{d}"
    );
    assert_eq!(v.uff_form("hz"), Some([0.0, 0.0]));

    let refused = |text: String| VanDerWaals::from_parm(&text).expect_err("refused");
    assert_eq!(
        refused(edit(PARM, "MOD4      RE\n", "")),
        ParmError::NoBlock
    );
    assert_eq!(
        refused(edit(PARM, "MOD4      RE", "MOD4      SK")),
        ParmError::Kind {
            line: 5,
            kind: "SK".into()
        }
    );
    let e = refused(edit(PARM, "1.4500  0.0100             invented", "1.4500"));
    assert!(matches!(e, ParmError::Short { line: 7, .. }), "{e:?}");
    let e = refused(edit(PARM, "1.4500  0.0100", "1.4500  NaN"));
    assert!(
        matches!(&e, ParmError::NotFinite { line: 7, text } if text == "NaN"),
        "{e:?}"
    );
    for (from, to) in [
        ("1.4500  0.0100", "-1.4500  0.0100"),
        ("1.4500  0.0100", "1.4500  -0.0100"),
        ("1.4500  0.0100", "0.0000  0.0100"),
    ] {
        let e = refused(edit(PARM, from, to));
        assert!(
            matches!(&e, ParmError::Value { kind, .. } if kind == "ha"),
            "{to}: {e:?}"
        );
    }
    assert_eq!(
        refused(edit(PARM, "  hz          0.0000", "  ca          0.0000")),
        ParmError::Repeated { kind: "ca".into() }
    );
    assert!(matches!(
        VanDerWaals::from_entries(vec![("ca".into(), f64::INFINITY, 0.1)]),
        Err(ParmError::Value { .. })
    ));
    assert_eq!(
        refused(edit(PARM, "MOD4      RE", "MOD4")),
        ParmError::Kind {
            line: 5,
            kind: String::new()
        }
    );

    // A second block is read too, and a type in both is refused.
    let two = format!("{PARM}\nMOD4      RE\n  c3   1.7000  0.1200   invented\n\nEND\n");
    let v2 = VanDerWaals::from_parm(&two).expect("two blocks");
    let kinds: Vec<&str> = v2.entries().iter().map(|e| e.0.as_str()).collect();
    assert_eq!(kinds, ["ca", "ha", "hz", "c3"]);
    assert_eq!(
        VanDerWaals::from_parm(&edit(&two, "  c3   1.7000", "  ha   1.7000")),
        Err(ParmError::Repeated { kind: "ha".into() })
    );

    // A frcmod's NONBON block, read alone and laid over the file it modifies.
    let frcmod = "remark goes here\nMASS\nca 99.99 9.999\n\nNONBON\n  ca   2.0000  0.0800\n  n4   1.8000  0.1700\n\n";
    let over = VanDerWaals::from_parm(frcmod).expect("a frcmod");
    let kinds: Vec<&str> = over.entries().iter().map(|e| e.0.as_str()).collect();
    assert_eq!(kinds, ["ca", "n4"], "the MASS line is not read");
    let laid = v.with_overrides(&over);
    let kinds: Vec<&str> = laid.entries().iter().map(|e| e.0.as_str()).collect();
    assert_eq!(kinds, ["ca", "ha", "hz", "n4"]);
    assert_eq!(laid.get("ca"), Some((2.0, 0.08)));
    assert_eq!(laid.get("ha"), v.get("ha"));
}

/// The dictionary's benzene in a box with two waters, as [`NonbondedModel::apply`] takes it.
fn benzene_with_waters(bnz: &Component) -> PeriodicForceField {
    let cell = PeriodicBox::cubic(20.0 * ANGSTROM);
    let p = EwaldParameters::for_accuracy(&cell, 9.0 * ANGSTROM, 1e-5);
    PeriodicForceField::solvated(bnz, &uff::assign(bnz), 2, cell, p).expect("supported")
}

/// **A model is the file's columns and the parameter file's values on the component**, and
/// **applied, it changes the solute and nothing else**: the waters' van der Waals and charges keep
/// the field's bits, the rules are Lorentz–Berthelot and [`AMBER_ONE_FOUR`], and benzene's 21
/// 1-4 pairs are scaled. A type the parameter file lacks, and a model of the wrong size, are
/// refused by name.
#[test]
fn a_model_applied_changes_the_solute_and_nothing_else() {
    let bnz = Component::from_ccd(BNZ).expect("BNZ");
    let file = Mol2::parse(SHUFFLED_BENZENE).expect("the shuffled benzene");
    let parm = VanDerWaals::from_parm(PARM).expect("the invented file");
    let model = NonbondedModel::new(&file, &bnz, &parm).expect("a model");
    let mapped = file.onto(&bnz).expect("benzene onto benzene");
    assert_eq!(model.atoms(), 12);
    assert_eq!(model.charges, mapped.charges);
    assert_eq!(model.types, mapped.kinds);
    for (i, a) in bnz.atoms().iter().enumerate() {
        let want = parm
            .uff_form(if a.element == Element::C { "ca" } else { "ha" })
            .expect("both types");
        assert_eq!(model.van_der_waals[i], want, "{}", a.name);
    }

    let before = benzene_with_waters(&bnz);
    let after = model.apply(before.clone()).expect("applied");
    let n = 12;
    assert_eq!(after.solute_atoms(), n);
    assert_eq!(&after.van_der_waals()[..n], &model.van_der_waals[..]);
    assert_eq!(&after.charges()[..n], &model.charges[..]);
    assert_eq!(&after.van_der_waals()[n..], &before.van_der_waals()[n..]);
    assert_eq!(&after.charges()[n..], &before.charges()[n..]);
    assert_eq!(after.van_der_waals().len(), n + 6);
    assert_eq!(after.mixing(), Mixing::LorentzBerthelot);
    assert_eq!(after.one_four(), AMBER_ONE_FOUR);
    assert_eq!(AMBER_ONE_FOUR, [0.5, 1.0 / 1.2]);
    assert_eq!(after.one_four_pairs().len(), 21);
    assert_eq!(before.mixing(), Mixing::Geometric);
    assert_eq!(before.one_four_pairs().len(), 0);
    // Undone by hand, the field is the one it was applied to, every field of it: the cell, the
    // bonded terms, the exclusions, the waters, Ewald.
    let undone = after
        .with_nonbonded(0..n, before.van_der_waals()[..n].to_vec())
        .expect("UFF's values")
        .with_solute_charges(before.charges()[..n].to_vec())
        .with_mixing(Mixing::Geometric)
        .with_one_four(1.0, 1.0);
    assert!(undone == before, "apply changed more than the model");

    // Two copies: the model repeated, on a field of two benzenes.
    let twice = model.repeated(2);
    assert_eq!(twice.atoms(), 2 * n);
    assert_eq!(&twice.charges[n..], &model.charges[..]);
    assert_eq!(&twice.van_der_waals[..n], &model.van_der_waals[..]);
    assert_eq!(&twice.types[n..], &model.types[..]);

    let no_ha = VanDerWaals::from_parm(&edit(
        PARM,
        "  ha          1.4500  0.0100             invented\n",
        "",
    ))
    .expect("ca and hz");
    let first_h = bnz
        .atoms()
        .iter()
        .position(|a| a.element == Element::H)
        .expect("a hydrogen");
    assert_eq!(
        NonbondedModel::new(&file, &bnz, &no_ha),
        Err(ModelError::NoSuchType {
            atom: first_h,
            name: bnz.atoms()[first_h].name.clone(),
            kind: "ha".into()
        })
    );
    let sybyl = Mol2::parse(
        &SHUFFLED_BENZENE
            .replace(" ha ", " H ")
            .replace(" ca ", " C.ar "),
    )
    .expect("SYBYL's types");
    let e = NonbondedModel::new(&sybyl, &bnz, &parm).expect_err("no C.ar in the file");
    assert!(e.to_string().contains("SYBYL type"), "{e}");

    // Charges that do not add to the molecule's net charge: one hydrogen's +0.1304 made +0.5304.
    let charged =
        Mol2::parse(&edit(SHUFFLED_BENZENE, "MOL  0.1304", "MOL  0.5304")).expect("well-formed");
    let e = NonbondedModel::new(&charged, &bnz, &parm).expect_err("a net charge of 0.4");
    assert!(
        matches!(e, ModelError::NetCharge { sum, formal: 0 } if (sum - 0.4).abs() < 1e-12),
        "{e:?}"
    );

    let mut ragged = model.clone();
    ragged.charges.pop();
    assert_eq!(
        ragged.apply(before.clone()),
        Err(ModelError::Ragged {
            types: 12,
            van_der_waals: 12,
            charges: 11
        })
    );

    let mut short = model.clone();
    short.types.pop();
    short.van_der_waals.pop();
    short.charges.pop();
    assert_eq!(
        short.apply(before),
        Err(ModelError::Atoms {
            model: 11,
            solute: 12
        })
    );
}
