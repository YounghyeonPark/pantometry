//! **Aspirin, read from its own dictionary entry and checked against the entry's other statements.**
//!
//! `AIN` is acetylsalicylic acid, C₉H₈O₄. Every check here is against something this crate did not
//! compute: the formula the same file states, the valence every organic chemist knows, and the
//! types a chemist reads off the structure. None of them is a second implementation of the parser.
//!
//! The fixture is embedded with `include_str!`, so these run under `wasm32` as well.

use pantometry_forcefield::{Component, Coordinates, Element, UffType};
use std::collections::BTreeMap;

const AIN: &str = include_str!("../components/AIN.cif");

fn aspirin() -> Component {
    Component::from_ccd(AIN).expect("AIN parses")
}

fn index(c: &Component, name: &str) -> usize {
    c.atoms()
        .iter()
        .position(|a| a.name == name)
        .unwrap_or_else(|| panic!("no atom {name}"))
}

/// **Twenty-one atoms and twenty-one bonds.** C₉H₈O₄ is 9 + 8 + 4 atoms, and a molecule with one
/// ring has exactly one more bond than a tree on the same atoms would: 20 + 1.
#[test]
fn aspirin_has_the_atoms_and_bonds_its_formula_and_one_ring_require() {
    let c = aspirin();
    assert_eq!(c.id(), "AIN");
    assert_eq!(c.atoms().len(), 9 + 8 + 4);
    // A connected molecule of n atoms with r independent rings has n - 1 + r bonds.
    assert_eq!(c.bonds().len(), c.atoms().len() - 1 + 1);
    assert_eq!(c.coordinates(), Coordinates::Ideal);
}

/// **The composition counted from the atoms is the formula the file states.** Two statements in one
/// file, made separately by the dictionary: the formula line and the atom list. The formula is
/// parsed here, independently of the crate, from its `C9 H8 O4` form.
#[test]
fn the_atoms_add_up_to_the_stated_formula() {
    let c = aspirin();
    let formula = c.formula().expect("AIN states a formula");
    assert_eq!(formula, "C9 H8 O4");
    let mut stated = BTreeMap::new();
    for term in formula.split_whitespace() {
        let split = term
            .find(|ch: char| ch.is_ascii_digit())
            .unwrap_or(term.len());
        let (symbol, count) = term.split_at(split);
        let count: usize = if count.is_empty() {
            1
        } else {
            count.parse().expect("a count")
        };
        let element = Element::from_symbol(symbol).expect("a known element");
        stated.insert(element, count);
    }
    assert_eq!(c.composition(), stated);

    // And the formal charges sum to the net charge the file states.
    let net: i32 = c.atoms().iter().map(|a| a.charge).sum();
    assert_eq!(Some(net), c.formal_charge());
    assert_eq!(net, 0, "aspirin as deposited is the neutral acid");
}

/// **Every hydrogen has exactly one bond**, which is what hydrogen is.
#[test]
fn every_hydrogen_has_one_bond() {
    let c = aspirin();
    for (i, atom) in c.atoms().iter().enumerate() {
        if atom.element == Element::H {
            assert_eq!(
                c.neighbours(i).count(),
                1,
                "{} has the wrong bonds",
                atom.name
            );
        }
    }
}

/// **Every heavy atom's bond orders sum to its valence** — carbon four, oxygen two — with an
/// aromatic bond counted as one and a half, the bond order of benzene.
///
/// The dictionary writes the ring in Kekulé form, so summing the *stated* orders would also come
/// to four; counting aromatic bonds as 1.5 is what makes this check depend on the aromatic flags
/// being right as well. A ring carbon with one aromatic flag missing sums to 3.5 or 4.5.
#[test]
fn every_heavy_atom_has_its_valence() {
    let c = aspirin();
    for (i, atom) in c.atoms().iter().enumerate() {
        let valence = match atom.element {
            Element::H => continue,
            Element::C => 4.0,
            Element::N => 3.0,
            Element::O => 2.0,
            other => panic!("aspirin has no {other}"),
        };
        let sum: f64 = c
            .neighbours(i)
            .map(|(_, b)| {
                if b.aromatic {
                    1.5
                } else {
                    f64::from(b.order.count())
                }
            })
            .sum();
        assert_eq!(sum, valence, "{} sums to {sum}", atom.name);
    }
}

/// **Aspirin's atoms have the types a chemist reads off its structure.**
///
/// Six ring carbons resonant; the acid's and the ester's carbonyl carbons sp²; the acetyl methyl
/// sp³; both carbonyl oxygens sp²; every hydrogen `H_`. **The acid's hydroxyl oxygen and the
/// ester's bridging oxygen are resonant, `O_R`** — each is divalent with single bonds and bonded
/// to an sp² atom (the acid's O1 to the carboxyl carbon C7; the ester's O3 to ring carbon C2 and
/// carbonyl carbon C8), so its lone pair is conjugated, and `uff::assign`'s resonant-heteroatom
/// rule types it as the paper types methyl vinyl ether's oxygen. Until that rule this test said
/// `O_3` for both, the reading of a structure drawn with sp³ ether oxygens; the paper's own
/// numbers (anisole's barrier, methyl vinyl ether's angle and bond) are what changed it. The atom
/// names are the dictionary's — `C7` is the acid carbon, `C8` the ester's, `C9` the methyl.
#[test]
fn aspirins_atoms_are_typed_as_its_structure_says() {
    let c = aspirin();
    let types = pantometry_forcefield::uff::assign(&c);
    let want = [
        ("C1", UffType::CR),
        ("C2", UffType::CR),
        ("C3", UffType::CR),
        ("C4", UffType::CR),
        ("C5", UffType::CR),
        ("C6", UffType::CR),
        ("C7", UffType::C2),
        ("C8", UffType::C2),
        ("C9", UffType::C3),
        ("O2", UffType::O2),
        ("O4", UffType::O2),
        ("O1", UffType::OR),
        ("O3", UffType::OR),
    ];
    for (name, t) in want {
        assert_eq!(types[index(&c, name)], t, "{name}");
    }
    for (i, atom) in c.atoms().iter().enumerate() {
        if atom.element == Element::H {
            assert_eq!(types[i], UffType::H, "{}", atom.name);
        }
    }
    // Thirteen heavy atoms and all thirteen are named above: nothing typed by default.
    assert_eq!(c.atoms().len() - 8, want.len());
}

/// **Positions are the ideal coordinates, in metres.** `O1` is written at
/// (1.731, 0.062, −2.912) Å.
#[test]
fn positions_are_the_ideal_coordinates_in_metres() {
    let c = aspirin();
    let o1 = c.atoms()[index(&c, "O1")].at;
    for (got, want) in o1.iter().zip([1.731e-10, 0.062e-10, -2.912e-10]) {
        assert!((got - want).abs() < 1e-22, "{got} against {want}");
    }
}

/// **One missing ideal coordinate moves the whole molecule to the model frame, and says so.**
/// Mixing the two would tear the molecule apart: in this entry they are twenty ångström apart.
#[test]
fn a_missing_ideal_coordinate_falls_back_to_the_model_for_every_atom() {
    let anchor = "2.659  0.080  -3.183";
    assert_eq!(
        AIN.matches(anchor).count(),
        1,
        "the anchor is HO1's ideal x y z"
    );
    let text = AIN.replace(anchor, "?      ?      ?");
    let c = Component::from_ccd(&text).expect("still parses, from the model coordinates");
    assert_eq!(c.coordinates(), Coordinates::Model);
    // O1's *model* coordinate, not its ideal one.
    let o1 = c.atoms()[index(&c, "O1")].at;
    for (got, want) in o1.iter().zip([13.907e-10, 16.130e-10, 0.624e-10]) {
        assert!((got - want).abs() < 1e-22, "{got} against {want}");
    }
}
