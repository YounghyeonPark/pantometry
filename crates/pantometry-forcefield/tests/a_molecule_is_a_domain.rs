//! **A molecule runs in a simulation, and answers as bodies with names and bonds.**
//!
//! What `pantometry-scene` reads to draw it — the bodies, their labels, their bonds — and what the
//! kernel reads to run it. For this step it does not move, and these hold it to that as well: a
//! domain claiming "static" that drifted would be lying about the one thing it says.

use pantometry_core::{Domain, Schedule, Simulation};
use pantometry_forcefield::{Component, Molecule};
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
