//! **A disulfide, which 181L does not have.** T4 lysozyme's pseudo-wild type has no cysteine, so
//! the disulfide path is held here on a structure built for it: two one-residue chains, each the
//! dictionary's cysteine at its ideal coordinates, the second turned half a revolution about an
//! axis through the point 1.02 Å beyond the first's SG along CB→SG — which puts the two SG
//! 2.04 Å apart, a disulfide's length, and keeps the second cysteine L (a turn, not a mirror).
//! Benzene is far away, because a system has a ligand.

mod protein;

use pantometry_forcefield::{
    uff, Element, ForceField, PdbError, Placement, Selection, System, UffType,
};
use protein::*;

fn unit(a: [f64; 3]) -> [f64; 3] {
    let l = len(a);
    [a[0] / l, a[1] / l, a[2] / l]
}

/// `x` turned half a revolution about the unit axis `u` through `p`.
fn half_turn(x: [f64; 3], p: [f64; 3], u: [f64; 3]) -> [f64; 3] {
    let d = sub(x, p);
    let along = dot(u, d);
    [
        p[0] + 2.0 * along * u[0] - d[0],
        p[1] + 2.0 * along * u[1] - d[1],
        p[2] + 2.0 * along * u[2] - d[2],
    ]
}

/// Where a chain's cysteine goes: a map from the template's ideal coordinates, Å.
type Place<'a> = &'a dyn Fn([f64; 3]) -> [f64; 3];

/// The text of a PDB file with one cysteine per chain in `chains`, each turned by its function,
/// and benzene 30 Å away.
fn structure(chains: &[(char, Place)]) -> String {
    let cys = template("CYS");
    let bnz = template("BNZ");
    let mut text = String::new();
    for (chain, _) in chains {
        text += &format!("SEQRES   1 {chain}    1  CYS\n");
    }
    let mut serial = 1;
    for (chain, place) in chains {
        for a in cys.atoms() {
            if a.element == Element::H || a.name == "OXT" {
                continue;
            }
            let at = place(a.at.map(|x| x / 1e-10));
            text += &atom_line(
                false,
                serial,
                &a.name,
                "CYS",
                *chain,
                1,
                at,
                a.element.symbol(),
            );
            text.push('\n');
            serial += 1;
        }
    }
    for a in bnz.atoms().iter().filter(|a| a.element != Element::H) {
        let at = a.at.map(|x| x / 1e-10 + 30.0);
        text += &atom_line(true, serial, &a.name, "BNZ", 'A', 900, at, "C");
        text.push('\n');
        serial += 1;
    }
    text
}

/// The point 1.02 Å past SG along CB→SG, and two axes through it perpendicular to that bond.
fn geometry() -> ([f64; 3], [f64; 3], [f64; 3]) {
    let at = |n: &str| template_at_angstrom(n);
    let d = unit(sub(at("SG"), at("CB")));
    let p = [
        at("SG")[0] + 1.02 * d[0],
        at("SG")[1] + 1.02 * d[1],
        at("SG")[2] + 1.02 * d[2],
    ];
    let u = unit(cross(d, [0.0, 0.0, 1.0]));
    let w = unit(cross(d, u));
    (p, u, w)
}

fn template_at_angstrom(name: &str) -> [f64; 3] {
    template("CYS")
        .atoms()
        .iter()
        .find(|a| a.name == name)
        .unwrap()
        .at
        .map(|x| x / 1e-10)
}

#[test]
fn a_disulfide_bonds_two_chains_and_takes_their_thiol_hydrogens() {
    let (p, u, _) = geometry();
    let same = |x: [f64; 3]| x;
    let turned = move |x: [f64; 3]| half_turn(x, p, u);
    let text = structure(&[('A', &same), ('B', &turned)]);
    let s = System::from_pdb(&text, &templates(), &Selection::new("BNZ")).unwrap();
    let c = s.component();
    assert_eq!(s.disulfides().len(), 1);
    let [i, j] = c.bonds()[s.disulfides()[0]].atoms;
    assert_eq!(
        [c.atoms()[i].name.as_str(), c.atoms()[j].name.as_str()],
        ["A:CYS1:SG", "B:CYS1:SG"]
    );
    let d = len(sub(c.atoms()[i].at, c.atoms()[j].at)) / 1e-10;
    // The format writes three decimals: each coordinate is within 0.0005 Å of the one built, so
    // the distance is within 2 √3 × 0.0005 = 0.0017 Å of 2.04.
    assert!((d - 2.04).abs() <= 2.0 * 3f64.sqrt() * 0.0005, "{d}");
    assert!(s.peptide_bonds().is_empty());
    let types = uff::assign(c);
    let mut checked = 0;
    for r in s.residues().iter().filter(|r| r.name == "CYS") {
        assert!(r.disulfide && r.n_terminal && r.c_terminal, "{}", r.label());
        let names: Vec<&str> = r.atoms.clone().map(|i| s.atom_name(i)).collect();
        assert!(!names.contains(&"HG"), "{names:?}");
        // C3H5NOS in a chain, NH3+ (+2 H) and COO- (+1 O), less the thiol H: C3 H6 N O2 S.
        let count = |e: Element| {
            r.atoms
                .clone()
                .filter(|&i| c.atoms()[i].element == e)
                .count()
        };
        assert_eq!(
            [Element::C, Element::H, Element::N, Element::O, Element::S].map(count),
            [3, 6, 1, 2, 1],
            "{}",
            r.label()
        );
        assert_eq!(r.charge, 0);
        let sg = r.atoms.clone().find(|&i| s.atom_name(i) == "SG").unwrap();
        let valence: u32 = c.neighbours(sg).map(|(_, b)| b.order.count()).sum();
        assert_eq!(valence, 2);
        assert_eq!(types[sg], UffType::S3Divalent);
        let oxt = r.atoms.clone().find(|&i| s.atom_name(i) == "OXT").unwrap();
        assert!(matches!(s.placements()[oxt], Placement::Template { .. }));
        checked += 1;
    }
    assert_eq!(checked, 2, "both cysteines were checked");
    assert!(ForceField::new(c, &types).is_ok());
    assert_eq!(s.formal_charge(), 0);
}

/// A third cysteine whose SG lands on the second's: the first SG then has two within reach.
#[test]
fn a_sulfur_with_two_partners_is_refused() {
    let (p, u, w) = geometry();
    let same = |x: [f64; 3]| x;
    let one = move |x: [f64; 3]| half_turn(x, p, u);
    let other = move |x: [f64; 3]| half_turn(x, p, w);
    let text = structure(&[('A', &same), ('B', &one), ('C', &other)]);
    assert_eq!(
        System::from_pdb(&text, &templates(), &Selection::new("BNZ")).unwrap_err(),
        PdbError::SulfurPartners {
            residue: "A:CYS1".into()
        }
    );
}
