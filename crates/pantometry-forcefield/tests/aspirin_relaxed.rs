//! **Aspirin relaxed from the dictionary's ideal coordinates: what the minimum looks like, term by
//! term, and how far it is from the crystal.**
//!
//! Reported, not asserted, beyond what minimisation guarantees: the energy falls and the force
//! ends below the tolerance. The numbers are facts about this force field's minimum nearest the
//! dictionary's geometry, and the comparison is with a crystal structure (PDB `1OXR`, which the
//! entry's `model_Cartn_*` columns come from), where neighbours and the protein pack against it;
//! no closed form says how close a gas-phase minimum should be to that.

mod common;

use common::{distance, index, positions, relaxed, rmsd, ANGSTROM};
use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::{Component, Element, Energy, ForceField, UffType};

const AIN: &str = include_str!("../components/AIN.cif");

/// The entry's `model_Cartn_x/y/z` for each atom, Å, read from its atom loop: the tenth to
/// twelfth fields of each row, which this entry's column order puts there (checked against
/// O1's 13.907, 16.130, 0.624).
fn model_coordinates(c: &Component) -> Vec<[f64; 3]> {
    let mut out = vec![[f64::NAN; 3]; c.atoms().len()];
    for line in AIN.lines() {
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() >= 18 && f[0] == "AIN" && f[17].parse::<usize>().is_ok() {
            if let Some(i) = c.atoms().iter().position(|a| a.name == f[1]) {
                out[i] = [
                    f[9].parse().expect("x"),
                    f[10].parse().expect("y"),
                    f[11].parse().expect("z"),
                ];
            }
        }
    }
    let o1 = index(c, "O1");
    assert_eq!(out[o1], [13.907, 16.130, 0.624]);
    assert!(out.iter().all(|p| p[0].is_finite()));
    out
}

fn heavy(c: &Component, at: &[[f64; 3]]) -> Vec<[f64; 3]> {
    c.atoms()
        .iter()
        .zip(at)
        .filter(|(a, _)| a.element != Element::H)
        .map(|(_, p)| *p)
        .collect()
}

fn kcal(e: Energy) -> [f64; 7] {
    let k = |x: f64| x / KCAL_PER_MOL;
    [
        k(e.bond),
        k(e.angle),
        k(e.torsion),
        k(e.inversion),
        k(e.van_der_waals),
        k(e.electrostatic),
        k(e.total),
    ]
}

#[test]
fn aspirin_relaxed_from_the_ideal_coordinates() {
    let c = Component::from_ccd(AIN).expect("AIN parses");
    let model = model_coordinates(&c);
    let start = positions(&c);
    let in_angstrom = |at: &[[f64; 3]]| -> Vec<[f64; 3]> {
        at.iter()
            .map(|p| [p[0] / ANGSTROM, p[1] / ANGSTROM, p[2] / ANGSTROM])
            .collect()
    };
    let (o4, h1) = (index(&c, "O4"), index(&c, "H1"));
    let rmsd_start = rmsd(&heavy(&c, &in_angstrom(&start)), &heavy(&c, &model));
    println!(
        "start (ideal): O4-H1 {:.3} Å, heavy-atom RMSD to the model coordinates {rmsd_start:.3} Å",
        distance(&start, o4, h1)
    );
    let mut types = uff::assign(&c);
    for (label, retype) in [
        ("this crate's typing", false),
        (
            "O1 and O3 typed O_R, as the paper appears to type an O on an sp2 carbon",
            true,
        ),
    ] {
        if retype {
            types[index(&c, "O1")] = UffType::OR;
            types[index(&c, "O3")] = UffType::OR;
        }
        let ff = ForceField::with_variant(&c, &types, Default::default()).expect("supported");
        let before = kcal(ff.energy(&start));
        let (e_after, at) = relaxed(&ff, &start);
        let after = kcal(ff.energy(&at));
        assert!(after[6] < before[6]);
        assert!((after[6] - e_after).abs() < 1e-9);
        let r = rmsd(&heavy(&c, &in_angstrom(&at)), &heavy(&c, &model));
        println!("{label}:");
        println!("| term | ideal start | relaxed |");
        println!("| --- | --- | --- |");
        for (k, name) in [
            "bond",
            "angle",
            "torsion",
            "inversion",
            "van der Waals",
            "electrostatic",
            "total",
        ]
        .iter()
        .enumerate()
        {
            println!("| {name} | {:.4} | {:.4} |", before[k], after[k]);
        }
        println!(
            "relaxed: O4-H1 {:.3} Å; heavy-atom RMSD to the model coordinates {r:.3} Å (from the \
             ideal start, {rmsd_start:.3} Å); C2-O3-C8-O4 {:.1}°, C3-C7-O1-HO1 {:.1}°",
            distance(&at, o4, h1),
            pantometry_forcefield::minimise::dihedral(
                &at,
                ["C2", "O3", "C8", "O4"].map(|n| index(&c, n))
            )
            .to_degrees(),
            pantometry_forcefield::minimise::dihedral(
                &at,
                ["C3", "C7", "O1", "HO1"].map(|n| index(&c, n))
            )
            .to_degrees()
        );
    }
}
