//! **A minimising molecule is a picture that changes: each captured frame is one minimiser step
//! further downhill.**
//!
//! `pantometry-forcefield` makes a kernel step one L-BFGS iteration, so the viewer's frames are
//! the molecule relaxing. This is the test of that claim from where the viewer stands — what
//! `capture` produces — rather than from the domain's own accessors: the positions in the panel
//! move from frame to frame, and the `energy` reading beside them never rises, which is the one
//! thing the minimiser guarantees (an Armijo line search accepts a step only if it lowers the
//! energy). Aspirin from the dictionary's ideal coordinates, which start with a 134 kcal/mol
//! clash, so the first frames move visibly.

use pantometry::prelude::*;
use std::collections::BTreeMap;

const AIN: &str = include_str!("../../pantometry-forcefield/components/AIN.cif");

#[test]
fn aspirin_relaxes_frame_by_frame() {
    let component = Component::from_ccd(AIN).expect("AIN parses");
    let mut sim = Simulation::new(Schedule::OneWay).with(Molecule::new("aspirin", component));
    let mut placed = BTreeMap::new();
    placed.insert("aspirin".to_string(), Placement::default());

    let frame = |sim: &Simulation| {
        let frame = capture(sim, &placed);
        let panel = frame
            .panels
            .iter()
            .find(|p| p.name == "aspirin")
            .expect("the molecule was drawn");
        let PanelData::Points { positions, .. } = &panel.data else {
            panic!("a molecule is bodies");
        };
        let energy = frame
            .readings
            .iter()
            .find(|r| r.domain == "aspirin" && r.label == "energy")
            .expect("an energy reading")
            .value;
        (positions.clone(), energy)
    };

    let (mut last_positions, mut last_energy) = frame(&sim);
    let first_energy = last_energy;
    for k in 0..10 {
        sim.advance(Time::s(1.0))
            .expect("a minimising molecule runs");
        let (positions, energy) = frame(&sim);
        assert_ne!(positions, last_positions, "frame {k}: nothing moved");
        assert!(
            energy < last_energy,
            "frame {k}: {energy} after {last_energy}"
        );
        last_positions = positions;
        last_energy = energy;
    }
    println!("aspirin over ten frames: {first_energy:.3} → {last_energy:.3} kcal/mol");
}
