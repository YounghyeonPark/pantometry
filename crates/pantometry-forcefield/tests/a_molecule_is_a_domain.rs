//! **A molecule runs in a simulation, and answers as bodies with names and bonds.**
//!
//! What `pantometry-scene` reads to draw it — the bodies, their labels, their bonds — and what the
//! kernel reads to run it. For this step it does not move, and these hold it to that as well: a
//! domain claiming "static" that drifted would be lying about the one thing it says.

use pantometry_core::{Domain, Schedule, Simulation};
use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::{Component, ForceField, Molecule, Unsupported};
use pantometry_units::Time;

const AIN: &str = include_str!("../components/AIN.cif");

fn aspirin() -> Molecule {
    Molecule::new("aspirin", Component::from_ccd(AIN).expect("AIN parses"))
}

#[test]
fn the_bodies_are_the_atoms_with_their_names_and_bonds() {
    let m = aspirin();
    let bodies = m.as_bodies().expect("a molecule is bodies");
    assert_eq!(bodies.count(), 21);
    assert_eq!(bodies.label(0).as_deref(), Some("O1"));
    assert_eq!(bodies.label(20).as_deref(), Some("H93"));
    let bonds = bodies.bonds();
    assert_eq!(bonds.len(), 21);
    assert!(bonds.iter().flatten().all(|&i| (i as usize) < 21));
    // O1 is bonded to C7 and HO1, the first two bonds in the file.
    assert_eq!(bonds[0], [0, 1]);
    assert_eq!(bonds[1], [0, 13]);
    // Coloured by atomic number: O1 is oxygen, C7 carbon, HO1 hydrogen.
    assert_eq!(
        [bodies.value(0), bodies.value(1), bodies.value(13)],
        [8.0, 6.0, 1.0]
    );
    let o1 = bodies.position(0);
    assert!((o1.x().to_si() - 1.731e-10).abs() < 1e-22);
}

#[test]
fn it_runs_balances_and_stays_where_it_is() {
    let before: Vec<[f64; 3]> = {
        let m = aspirin();
        (0..21).map(|i| m.at(i)).collect()
    };
    let mut sim = Simulation::new(Schedule::OneWay).with(aspirin());
    for _ in 0..3 {
        sim.advance(Time::s(1e-15)).expect("a static molecule runs");
    }
    let m = sim
        .domain_as::<Molecule>("aspirin")
        .expect("reachable as itself");
    let after: Vec<[f64; 3]> = (0..21).map(|i| m.at(i)).collect();
    assert_eq!(before, after, "step 1a moves nothing");

    let readings = m.readings();
    let value = |label: &str| {
        readings
            .iter()
            .find(|r| r.label == label)
            .unwrap_or_else(|| panic!("no {label} reading"))
            .value
    };
    assert_eq!(value("atoms"), 21.0);
    assert_eq!(value("heavy atoms"), 13.0);
    assert_eq!(value("bonds"), 21.0);
    assert_eq!(value("formal charge"), 0.0);
    assert!(m.books_balance());
    assert!(m.supports_restore());
}

/// **The energy readings are the force field's energy, in kcal/mol**, and the total is the sum
/// of the six terms read beside it. The conversion is one division by an exact constant, so the
/// readings equal the evaluation divided by it, bit for bit.
#[test]
fn the_energy_is_read_out_in_kcal_per_mol() {
    let m = aspirin();
    let e = m.evaluate().expect("aspirin is supported").energy;
    let readings = m.readings();
    let reading = |label: &str| {
        readings
            .iter()
            .find(|r| r.label == label)
            .unwrap_or_else(|| panic!("no {label} reading"))
    };
    for (label, joules) in [
        ("energy", e.total),
        ("bond stretch", e.bond),
        ("angle bend", e.angle),
        ("torsion", e.torsion),
        ("inversion", e.inversion),
        ("van der Waals", e.van_der_waals),
        ("electrostatic", e.electrostatic),
    ] {
        let r = reading(label);
        assert_eq!(r.unit, "kcal/mol", "{label}");
        assert_eq!(r.value, joules / KCAL_PER_MOL, "{label}");
        assert!(r.value.is_finite(), "{label}");
    }
    assert_eq!(
        reading("electrostatic").value,
        0.0,
        "charges default to zero"
    );
}

/// A dimethyl sulfide skeleton, C1–S1–C2, with `oxygens` S=O bonds added (0, 1 or 2) and, if
/// `third_methyl`, a third carbon C3 single-bonded to the sulfur with formal charge +1 on it.
fn sulfur_compound(oxygens: usize, third_methyl: bool) -> Component {
    let mut s = String::from(
        "data_TST\n_chem_comp.id TST\nloop_\n_chem_comp_atom.comp_id\n\
         _chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n\
         _chem_comp_atom.pdbx_aromatic_flag\n_chem_comp_atom.pdbx_model_Cartn_x_ideal\n\
         _chem_comp_atom.pdbx_model_Cartn_y_ideal\n_chem_comp_atom.pdbx_model_Cartn_z_ideal\n\
         TST C1 C 0 N -1.5 0.0 0.0\nTST C2 C 0 N 1.5 0.0 0.0\n",
    );
    s += if third_methyl {
        "TST S1 S 1 N 0.0 0.5 0.0\nTST C3 C 0 N 0.0 1.0 -1.6\n"
    } else {
        "TST S1 S 0 N 0.0 0.5 0.0\n"
    };
    for (k, z) in [1.2, -1.2].iter().take(oxygens).enumerate() {
        s += &format!("TST O{} O 0 N 0.0 1.2 {z}\n", k + 1);
    }
    s += "loop_\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n\
          _chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n\
          C1 S1 SING N\nS1 C2 SING N\n";
    if third_methyl {
        s += "S1 C3 SING N\n";
    }
    for k in 0..oxygens {
        s += &format!("S1 O{} DOUB N\n", k + 1);
    }
    Component::from_ccd(&s).expect("parses")
}

/// **A sulfoxide and a sulfone are refused by name, and still drawn.** Dimethyl sulfoxide's
/// sulfur has one S=O, a bond-order sum of four; dimethyl sulfone's has two, six. UFF's types for
/// them (`S_3+4`, `S_3+6`) are not in this crate, and the divalent type they would otherwise get
/// would give them the wrong radius. So [`ForceField::new`] refuses each, naming the type it would
/// need; the molecule is still a molecule, and its energy readings are `NaN` rather than a number
/// computed wrong. Dimethyl sulfide, the same skeleton without the oxygens, is accepted — so the
/// threshold sits between two and four, not merely below six.
#[test]
fn a_hypervalent_sulfur_is_refused_by_name() {
    for (oxygens, valence, names) in [(1, 4, "S_3+4"), (2, 6, "S_3+6")] {
        let c = sulfur_compound(oxygens, false);
        let refused = ForceField::new(&c, &uff::assign(&c));
        let expected = Unsupported::HypervalentSulfur {
            atom: "S1".into(),
            valence,
        };
        assert_eq!(refused, Err(expected.clone()));
        let message = expected.to_string();
        assert!(message.contains(names), "{message}");
        assert!(message.contains("S1"), "{message}");

        let m = Molecule::new("hypervalent", c);
        assert!(m.force_field().is_err());
        assert_eq!(
            m.as_bodies().expect("still bodies").count(),
            3 + oxygens,
            "{valence}"
        );
        let energy = m
            .readings()
            .into_iter()
            .find(|r| r.label == "energy")
            .expect("an energy reading");
        assert!(energy.value.is_nan());
    }

    let sulfide = sulfur_compound(0, false);
    assert!(ForceField::new(&sulfide, &uff::assign(&sulfide)).is_ok());
}

/// **A sulfonium — three single bonds, sum three — is refused too, and the message does not call
/// it a sulfoxide or sulfone.** It is neither `S_3+2` nor either hypervalent type this crate lacks,
/// and a message naming `S_3+4` for it would send a reader to the wrong parameters.
#[test]
fn a_sulfonium_is_refused_and_named_for_what_it_is() {
    let c = sulfur_compound(0, true);
    let expected = Unsupported::HypervalentSulfur {
        atom: "S1".into(),
        valence: 3,
    };
    assert_eq!(ForceField::new(&c, &uff::assign(&c)), Err(expected.clone()));
    let message = expected.to_string();
    assert!(message.contains("sum 3"), "{message}");
    assert!(message.contains("neither divalent"), "{message}");
    assert!(
        !message.contains("S_3+4") && !message.contains("S_3+6"),
        "{message}"
    );
}
