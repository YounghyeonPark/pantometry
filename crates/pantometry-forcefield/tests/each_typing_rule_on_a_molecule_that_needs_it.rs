//! **Each typing rule, on the smallest molecule whose chemistry decides it.**
//!
//! Aspirin exercises three carbon types and two oxygen ones. The rest — the amide nitrogen, a
//! nitrile, an imine, a pyridine, a phosphate, carbon monoxide, a thiol, a thiophene, the halogens —
//! each have a textbook molecule whose hybridisation is not in question, and that molecule is the
//! check. The entries are written here in the dictionary's format, with explicit hydrogens and the
//! dictionary's conventions (Kekulé rings with aromatic flags), and coordinates that are irrelevant
//! to typing.

use pantometry_forcefield::{uff, Component, UffType};

/// A one-component CCD file: atoms as `(name, element, aromatic)`, bonds as
/// `(first, second, order, aromatic)`.
fn entry(atoms: &[(&str, &str, bool)], bonds: &[(&str, &str, &str, bool)]) -> String {
    let yn = |b: bool| if b { "Y" } else { "N" };
    let mut s = String::from(
        "data_TST\n_chem_comp.id TST\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n\
         _chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.pdbx_aromatic_flag\n\
         _chem_comp_atom.pdbx_model_Cartn_x_ideal\n_chem_comp_atom.pdbx_model_Cartn_y_ideal\n\
         _chem_comp_atom.pdbx_model_Cartn_z_ideal\n",
    );
    for (k, (name, element, aromatic)) in atoms.iter().enumerate() {
        s += &format!("TST {name} {element} 0 {} {k}.0 0.0 0.0\n", yn(*aromatic));
    }
    s += "loop_\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n\
          _chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n";
    for (a, b, order, aromatic) in bonds {
        s += &format!("{a} {b} {order} {}\n", yn(*aromatic));
    }
    s
}

fn types(atoms: &[(&str, &str, bool)], bonds: &[(&str, &str, &str, bool)]) -> Vec<UffType> {
    let c = Component::from_ccd(&entry(atoms, bonds)).expect("the fixture parses");
    uff::assign(&c)
}

/// **Acetamide's nitrogen is resonant, and methylamine's is not.** CH₃C(=O)NH₂ is planar at the
/// nitrogen because the C–N bond is partly double; CH₃NH₂ is pyramidal. Their bond orders, as
/// written, are identical at the nitrogen — only the neighbouring C=O tells them apart.
#[test]
fn an_amide_nitrogen_is_resonant_and_an_amine_is_not() {
    let acetamide = types(
        &[
            ("C1", "C", false),
            ("C2", "C", false),
            ("O", "O", false),
            ("N", "N", false),
            ("H1", "H", false),
            ("H2", "H", false),
            ("H3", "H", false),
            ("HN1", "H", false),
            ("HN2", "H", false),
        ],
        &[
            ("C1", "C2", "SING", false),
            ("C2", "O", "DOUB", false),
            ("C2", "N", "SING", false),
            ("C1", "H1", "SING", false),
            ("C1", "H2", "SING", false),
            ("C1", "H3", "SING", false),
            ("N", "HN1", "SING", false),
            ("N", "HN2", "SING", false),
        ],
    );
    assert_eq!(
        acetamide[..4],
        [UffType::C3, UffType::C2, UffType::O2, UffType::NR]
    );

    let methylamine = types(
        &[
            ("C", "C", false),
            ("N", "N", false),
            ("H1", "H", false),
            ("H2", "H", false),
            ("H3", "H", false),
            ("HN1", "H", false),
            ("HN2", "H", false),
        ],
        &[
            ("C", "N", "SING", false),
            ("C", "H1", "SING", false),
            ("C", "H2", "SING", false),
            ("C", "H3", "SING", false),
            ("N", "HN1", "SING", false),
            ("N", "HN2", "SING", false),
        ],
    );
    assert_eq!(methylamine[..2], [UffType::C3, UffType::N3]);
}

/// **Hydrogen cyanide is linear at both heavy atoms**; **methanimine is trigonal**; **carbon
/// monoxide's oxygen is sp**; **carbon dioxide's carbon is sp** by its two double bonds.
#[test]
fn multiple_bonds_set_the_hybridisation() {
    let hcn = types(
        &[("H", "H", false), ("C", "C", false), ("N", "N", false)],
        &[("H", "C", "SING", false), ("C", "N", "TRIP", false)],
    );
    assert_eq!(hcn, [UffType::H, UffType::C1, UffType::N1]);

    let methanimine = types(
        &[
            ("C", "C", false),
            ("N", "N", false),
            ("H1", "H", false),
            ("H2", "H", false),
            ("HN", "H", false),
        ],
        &[
            ("C", "N", "DOUB", false),
            ("C", "H1", "SING", false),
            ("C", "H2", "SING", false),
            ("N", "HN", "SING", false),
        ],
    );
    assert_eq!(methanimine[..2], [UffType::C2, UffType::N2]);

    let co = types(
        &[("C", "C", false), ("O", "O", false)],
        &[("C", "O", "TRIP", false)],
    );
    assert_eq!(co, [UffType::C1, UffType::O1]);

    let co2 = types(
        &[("O1", "O", false), ("C", "C", false), ("O2", "O", false)],
        &[("O1", "C", "DOUB", false), ("C", "O2", "DOUB", false)],
    );
    assert_eq!(co2, [UffType::O2, UffType::C1, UffType::O2]);
}

/// **Pyridine's nitrogen is aromatic although the file writes it with a double bond**, and
/// **thiophene's and furan's heteroatoms are aromatic although the file writes them single.** The
/// aromatic flag is read before the order, which is the whole reason for its rule's position.
#[test]
fn an_aromatic_heteroatom_is_resonant_whatever_its_kekule_order() {
    // Pyridine, Kekulé: N1=C2, C3=C4, C5=C6. Hydrogens left off — typing only reads the ring.
    let ring = |hetero: (&str, &str)| {
        let atoms = [
            (hetero.0, hetero.1, true),
            ("C2", "C", true),
            ("C3", "C", true),
            ("C4", "C", true),
            ("C5", "C", true),
            ("C6", "C", true),
        ];
        let first = if hetero.1 == "N" { "DOUB" } else { "SING" };
        let bonds = [
            (hetero.0, "C2", first, true),
            ("C2", "C3", "SING", true),
            ("C3", "C4", "DOUB", true),
            ("C4", "C5", "SING", true),
            ("C5", "C6", "DOUB", true),
            ("C6", hetero.0, "SING", true),
        ];
        types(&atoms, &bonds)
    };
    assert_eq!(ring(("N1", "N"))[0], UffType::NR);
    assert!(ring(("N1", "N"))[1..].iter().all(|t| *t == UffType::CR));

    // Five-membered: S1-C2=C3-C4=C5-S1.
    let five = |hetero: (&str, &str)| {
        let atoms = [
            (hetero.0, hetero.1, true),
            ("C2", "C", true),
            ("C3", "C", true),
            ("C4", "C", true),
            ("C5", "C", true),
        ];
        let bonds = [
            (hetero.0, "C2", "SING", true),
            ("C2", "C3", "DOUB", true),
            ("C3", "C4", "SING", true),
            ("C4", "C5", "DOUB", true),
            ("C5", hetero.0, "SING", true),
        ];
        types(&atoms, &bonds)[0]
    };
    assert_eq!(five(("S1", "S")), UffType::SR);
    assert_eq!(five(("O1", "O")), UffType::OR);
}

/// **Methanethiol's sulfur is divalent sp³; thioformaldehyde's is sp²; phosphoric acid's
/// phosphorus is +5 and phosphine's +3; and each halogen has its one type.**
#[test]
fn sulfur_phosphorus_and_the_halogens() {
    let thiol = types(
        &[("C", "C", false), ("S", "S", false), ("HS", "H", false)],
        &[("C", "S", "SING", false), ("S", "HS", "SING", false)],
    );
    assert_eq!(thiol[1], UffType::S3Divalent);

    let thioformaldehyde = types(
        &[("C", "C", false), ("S", "S", false)],
        &[("C", "S", "DOUB", false)],
    );
    assert_eq!(thioformaldehyde, [UffType::C2, UffType::S2]);

    // H3PO4: P=O and three P-OH.
    let phosphoric = types(
        &[
            ("P", "P", false),
            ("O1", "O", false),
            ("O2", "O", false),
            ("O3", "O", false),
            ("O4", "O", false),
        ],
        &[
            ("P", "O1", "DOUB", false),
            ("P", "O2", "SING", false),
            ("P", "O3", "SING", false),
            ("P", "O4", "SING", false),
        ],
    );
    assert_eq!(
        phosphoric,
        [
            UffType::P3Pentavalent,
            UffType::O2,
            UffType::O3,
            UffType::O3,
            UffType::O3
        ]
    );

    let phosphine = types(
        &[
            ("P", "P", false),
            ("H1", "H", false),
            ("H2", "H", false),
            ("H3", "H", false),
        ],
        &[
            ("P", "H1", "SING", false),
            ("P", "H2", "SING", false),
            ("P", "H3", "SING", false),
        ],
    );
    assert_eq!(phosphine[0], UffType::P3Trivalent);

    // CF Cl Br I on one carbon, the dictionary's upper-case spelling of the symbols.
    let halo = types(
        &[
            ("C", "C", false),
            ("F", "F", false),
            ("CL", "CL", false),
            ("BR", "BR", false),
            ("I", "I", false),
        ],
        &[
            ("C", "F", "SING", false),
            ("C", "CL", "SING", false),
            ("C", "BR", "SING", false),
            ("C", "I", "SING", false),
        ],
    );
    assert_eq!(
        halo,
        [
            UffType::C3,
            UffType::F,
            UffType::Cl,
            UffType::Br,
            UffType::I
        ]
    );
}

/// **A quoted atom name with a prime in it**, as nucleotides are written (`"C1'"`), is one name and
/// bonds to it resolve.
#[test]
fn a_primed_atom_name_is_one_name() {
    let text = entry(
        &[("\"C1'\"", "C", false), ("O", "O", false)],
        &[("\"C1'\"", "O", "SING", false)],
    );
    let c = Component::from_ccd(&text).expect("parses");
    assert_eq!(c.atoms()[0].name, "C1'");
    assert_eq!(c.bonds()[0].atoms, [0, 1]);
}
