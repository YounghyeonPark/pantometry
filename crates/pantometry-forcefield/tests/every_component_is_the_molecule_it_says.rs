//! **Every chemical component the tests use is the molecule it is used as**: its parsed
//! composition is its own stated formula, every hydrogen has exactly one bond, and every atom is
//! typed — a fixture that is a different molecule would make every comparison built on it
//! meaningless without failing one. The dictionary entries state their formula; the hand-built
//! ones are checked against the formula of the molecule their first line names, typed here.

mod protein;

use pantometry_forcefield::{uff, Component, Element};
use std::collections::BTreeMap;

/// `"C2 H6 O"` as element counts.
fn formula(text: &str) -> BTreeMap<Element, usize> {
    text.trim_matches('"')
        .split_whitespace()
        .map(|part| {
            let digits = part
                .find(|c: char| c.is_ascii_digit())
                .unwrap_or(part.len());
            let (symbol, count) = part.split_at(digits);
            (
                Element::from_symbol(symbol).expect("an element this crate reads"),
                if count.is_empty() {
                    1
                } else {
                    count.parse().expect("a count")
                },
            )
        })
        .collect()
}

fn check(id: &str, text: &str, stated: &str) {
    let c = Component::from_ccd(text).unwrap_or_else(|e| panic!("{id}: {e}"));
    assert_eq!(c.composition(), formula(stated), "{id}");
    for (i, a) in c.atoms().iter().enumerate() {
        if a.element == Element::H {
            assert_eq!(c.neighbours(i).count(), 1, "{id}: {} has one bond", a.name);
        }
    }
    assert_eq!(uff::assign(&c).len(), c.atoms().len(), "{id}");
}

#[test]
fn every_component_is_the_molecule_it_says() {
    let mut seen: Vec<&str> = Vec::new();
    for (id, text) in [
        ("NME", include_str!("../components/NME.cif")),
        ("MOH", include_str!("../components/MOH.cif")),
        ("MEE", include_str!("../components/MEE.cif")),
        ("PEO", include_str!("../components/PEO.cif")),
        ("S2H", include_str!("../components/S2H.cif")),
        ("A1JFW", include_str!("../components/A1JFW.cif")),
        ("16R", include_str!("../components/16R.cif")),
        ("ACE", include_str!("../components/ACE.cif")),
        ("61G", include_str!("../components/61G.cif")),
        ("PYJ", include_str!("../components/PYJ.cif")),
        ("2F2", include_str!("../components/2F2.cif")),
        ("2ME", include_str!("../components/2ME.cif")),
        ("ACN", include_str!("../components/ACN.cif")),
        ("DMN", include_str!("../components/DMN.cif")),
        ("KEN", include_str!("../components/KEN.cif")),
        ("CCN", include_str!("../components/CCN.cif")),
        ("ACM", include_str!("../components/ACM.cif")),
        ("TME", include_str!("../components/TME.cif")),
        ("6AC", include_str!("../components/6AC.cif")),
        ("DMF", include_str!("../components/DMF.cif")),
        ("ALA", include_str!("../components/ALA.cif")),
        ("ARG", include_str!("../components/ARG.cif")),
        ("ASN", include_str!("../components/ASN.cif")),
        ("ASP", include_str!("../components/ASP.cif")),
        ("CYS", include_str!("../components/CYS.cif")),
        ("GLN", include_str!("../components/GLN.cif")),
        ("GLU", include_str!("../components/GLU.cif")),
        ("GLY", include_str!("../components/GLY.cif")),
        ("HIS", include_str!("../components/HIS.cif")),
        ("ILE", include_str!("../components/ILE.cif")),
        ("LEU", include_str!("../components/LEU.cif")),
        ("LYS", include_str!("../components/LYS.cif")),
        ("MET", include_str!("../components/MET.cif")),
        ("PHE", include_str!("../components/PHE.cif")),
        ("PRO", include_str!("../components/PRO.cif")),
        ("SER", include_str!("../components/SER.cif")),
        ("THR", include_str!("../components/THR.cif")),
        ("TRP", include_str!("../components/TRP.cif")),
        ("TYR", include_str!("../components/TYR.cif")),
        ("VAL", include_str!("../components/VAL.cif")),
        ("BNZ", include_str!("../components/BNZ.cif")),
        ("BZF", include_str!("../components/BZF.cif")),
        ("DEN", include_str!("../components/DEN.cif")),
        ("I4B", include_str!("../components/I4B.cif")),
        ("IND", include_str!("../components/IND.cif")),
        ("N4B", include_str!("../components/N4B.cif")),
        ("PXY", include_str!("../components/PXY.cif")),
        ("OXE", include_str!("../components/OXE.cif")),
    ] {
        let c = Component::from_ccd(text).expect("parses");
        let stated = c
            .formula()
            .expect("the entry states its formula")
            .to_string();
        assert_eq!(c.id(), id);
        check(id, text, &stated);
        seen.push(id);
    }
    // Every template the protein tests use is held here too, so the two lists cannot drift.
    for code in protein::CODES {
        assert!(
            seen.contains(&code),
            "{code} is a protein template and is not checked here"
        );
    }
    for (id, text, stated) in [
        ("ethane", include_str!("hand_built/eta.cif"), "C2 H6"),
        (
            "methylphosphine",
            include_str!("hand_built/mph.cif"),
            "C H5 P",
        ),
        (
            "N-methylformamide",
            include_str!("hand_built/nmf.cif"),
            "C2 H5 N O",
        ),
        (
            "methyl formate",
            include_str!("hand_built/mfo.cif"),
            "C2 H4 O2",
        ),
        (
            "methyl vinyl ether",
            include_str!("hand_built/mve.cif"),
            "C3 H6 O",
        ),
        (
            "dimethyldiazene",
            include_str!("hand_built/dmd.cif"),
            "C2 H6 N2",
        ),
        ("propene", include_str!("hand_built/ppe.cif"), "C3 H6"),
        ("butadiene", include_str!("hand_built/bde.cif"), "C4 H6"),
        ("propyne", include_str!("hand_built/pyn.cif"), "C3 H4"),
    ] {
        assert!(text.starts_with("# HAND-BUILT"), "{id} says what it is");
        check(id, text, stated);
    }
}
