//! **A molecule reaches a picture without the scene layer knowing what a molecule is.**
//!
//! `pantometry-scene` names no domain: it asks each one what it offers and draws the shape that
//! comes back. So the claim that `pantometry-forcefield` is drawable is a claim about what
//! `capture` produces from it, and this is the test of that claim — the panel the viewer reads,
//! with every atom, its name, and the dictionary's bonds.

use pantometry::prelude::*;
use std::collections::BTreeMap;

const AIN: &str = include_str!("../../pantometry-forcefield/components/AIN.cif");

#[test]
fn aspirin_is_a_points_panel_with_its_names_and_bonds() {
    let component = Component::from_ccd(AIN).expect("AIN parses");
    let sim = Simulation::new(Schedule::OneWay).with(Molecule::new("aspirin", component));
    let mut placed = BTreeMap::new();
    placed.insert("aspirin".to_string(), Placement::default());
    let frame = capture(&sim, &placed);

    let panel = frame
        .panels
        .iter()
        .find(|p| p.name == "aspirin")
        .expect("the molecule was drawn");
    let PanelData::Points {
        positions,
        values,
        labels,
        bonds,
        ..
    } = &panel.data
    else {
        panic!("a molecule is bodies, not a field");
    };
    assert_eq!(positions.len(), 21);
    assert_eq!(values.len(), 21);
    assert_eq!(labels.len(), 21);
    assert_eq!(labels[0], "O1");
    assert_eq!(bonds.len(), 21);
    // C7=O2 is the third bond in the file: atoms 1 and 2.
    assert_eq!(bonds[2], [1, 2]);

    // And its numbers are in the frame, drawable or not.
    let atoms = frame
        .readings
        .iter()
        .find(|r| r.domain == "aspirin" && r.label == "atoms")
        .expect("an atom count");
    assert_eq!(atoms.value, 21.0);
}
