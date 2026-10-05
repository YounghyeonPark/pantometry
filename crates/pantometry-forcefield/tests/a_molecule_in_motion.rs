//! **A molecule in dynamics is a domain whose step is a time step**: the kernel subdivides a
//! frame into steps no longer than the molecule's, the energy books are audited — and refused at
//! the kernel's default tolerance, because no integrator conserves energy to 1e-9 — and the
//! readings say how hot it is and what the bath has done.
//!
//! The minimiser stays the default: a molecule that is not given dynamics steps as it always has,
//! which `a_molecule_is_a_domain.rs` and `the_minimiser.rs` hold.

use pantometry_core::conserved::quantity;
use pantometry_core::{Domain, Schedule, Simulation};
use pantometry_forcefield::{Bath, Component, Molecule};
use pantometry_units::Time;

const AIN: &str = include_str!("../components/AIN.cif");
const FS: f64 = 1e-15;

/// Aspirin relaxed, in a 300 K Langevin bath at 1 ps⁻¹.
fn moving(seed: u64) -> Molecule {
    let mut m = Molecule::new("aspirin", Component::from_ccd(AIN).expect("AIN parses"));
    m.minimise(50_000);
    m.thermalised(
        300.0,
        Bath::Langevin {
            temperature: 300.0,
            friction: 1e12,
            seed,
        },
        seed,
    )
}

fn value(m: &Molecule, label: &str) -> f64 {
    m.readings()
        .into_iter()
        .find(|r| r.label == label)
        .unwrap_or_else(|| panic!("no {label} reading"))
        .value
}

/// **Without dynamics nothing has changed**: the step is a minimiser iteration, the step length is
/// unbounded, the ledger is empty and the diagnostics are the minimiser's four.
#[test]
fn the_minimiser_is_still_the_default() {
    let m = Molecule::new("aspirin", Component::from_ccd(AIN).expect("parses"));
    assert!(m.dynamics().is_none());
    assert!(m.max_stable_dt(Time::ZERO).to_si().is_infinite());
    assert!(m.ledger().is_empty());
    assert_eq!(
        m.diagnostics(),
        ["max force", "rms force", "converged", "minimiser steps"]
    );
    assert!(m.readings().iter().all(|r| r.label != "temperature"));
}

/// **The default audit refuses a molecule in motion, by name**, and one set to the integrator's
/// band runs — in frames of one step and of 50 fs. Velocity Verlet keeps energy to `O(h²)`, not to
/// the 1e-9 a `Simulation` checks by default, so a run that did not say what it was willing to
/// accept would be checking nothing it could pass.
///
/// And the tolerance it is set to can see a bath whose work is half counted: on these runs that
/// moves the books by at least 1.06e-2 of their scale across a one-step frame and 6.7e-2 across a
/// 50 fs frame, against [`Molecule::ENERGY_TOLERANCE`]'s 7e-3 (`the_energy_tolerance_measured`).
#[test]
fn the_books_are_audited_and_the_tolerance_has_to_be_said() {
    let mut strict = Simulation::new(Schedule::OneWay).with(moving(1));
    let refused = (0..20)
        .find_map(|_| strict.advance(Time::s(0.5 * FS)).err())
        .expect("1e-9 is tighter than any integrator");
    assert_eq!(refused.quantity, quantity::ENERGY);

    let mut sim = Simulation::new(Schedule::OneWay)
        .with(moving(1))
        .conservation_tolerance_for(quantity::ENERGY, Molecule::ENERGY_TOLERANCE);
    for _ in 0..200 {
        sim.advance(Time::s(0.5 * FS))
            .expect("inside the integrator's band");
    }
    for _ in 0..20 {
        sim.advance(Time::s(50.0 * FS))
            .expect("inside the integrator's band");
    }
    let m = sim.domain_as::<Molecule>("aspirin").expect("itself");
    assert_eq!(m.dynamics().expect("in dynamics").steps(), 2200);
    assert!(value(m, "temperature").is_finite());
    for label in m.diagnostics() {
        assert!(m.readings().iter().any(|r| r.label == *label), "{label}");
    }
}

/// **A frame longer than the molecule's step is cut into steps of at most that length**, and a
/// run cut into frames of one length or two is the same run, to the bit.
#[test]
fn a_frame_is_cut_into_steps_and_the_cut_does_not_matter() {
    let run = |frame_fs: f64, frames: usize| {
        let mut sim = Simulation::new(Schedule::OneWay)
            .with(moving(2))
            .conservation_tolerance_for(quantity::ENERGY, Molecule::ENERGY_TOLERANCE);
        for _ in 0..frames {
            sim.advance(Time::s(frame_fs * FS)).expect("runs");
        }
        let m = sim.domain_as::<Molecule>("aspirin").expect("itself");
        let steps = m.dynamics().expect("in dynamics").steps();
        let bits: Vec<u64> = (0..21).flat_map(|i| m.at(i)).map(f64::to_bits).collect();
        (steps, bits)
    };
    let (s1, b1) = run(0.5, 100);
    let (s2, b2) = run(1.0, 50);
    let (s3, b3) = run(5.0, 10);
    assert_eq!((s1, s2, s3), (100, 100, 100));
    assert_eq!(b1, b2);
    assert_eq!(b1, b3);

    // The molecule says how long a step it can take, and cuts a longer one itself when it is
    // stepped directly rather than through a simulation that has already cut it.
    let mut m = moving(2);
    assert_eq!(
        m.max_stable_dt(Time::ZERO).to_si(),
        Molecule::DEFAULT_TIME_STEP
    );
    let mut bus = pantometry_core::Exchange::new();
    for _ in 0..10 {
        m.step(Time::ZERO, Time::s(5.0 * FS), &mut bus).unwrap();
    }
    assert_eq!(m.dynamics().unwrap().steps(), 100);
    let direct: Vec<u64> = (0..21).flat_map(|i| m.at(i)).map(f64::to_bits).collect();
    assert_eq!(direct, b3);
}

/// **A restored molecule repeats itself to the bit**, noise included: the step count the kicks are
/// keyed by is part of what is saved.
#[test]
fn a_restored_run_repeats_itself() {
    let mut m = moving(3);
    let mut bus = pantometry_core::Exchange::new();
    let h = Time::s(0.5 * FS);
    for _ in 0..50 {
        m.step(Time::ZERO, h, &mut bus).expect("steps");
    }
    m.checkpoint();
    let mut first = Vec::new();
    for _ in 0..50 {
        m.step(Time::ZERO, h, &mut bus).expect("steps");
        first.push(m.at(5));
    }
    m.restore();
    let mut second = Vec::new();
    for _ in 0..50 {
        m.step(Time::ZERO, h, &mut bus).expect("steps");
        second.push(m.at(5));
    }
    assert_eq!(first, second);
}

/// **Changing the charges of a molecule in motion changes the forces it moves under.** The
/// dynamics keeps its forces between steps, keyed on the positions alone, so a molecule that
/// changed its potential and kept them would take its next step under the old one. Charging
/// before and after thermalising has to give the same run.
#[test]
fn new_charges_are_new_forces() {
    let q: Vec<f64> = (0..21).map(|i| 0.05 * ((i % 5) as f64 - 2.0)).collect();
    let thermostat = Bath::Langevin {
        temperature: 300.0,
        friction: 1e12,
        seed: 4,
    };
    let c = Component::from_ccd(AIN).expect("parses");
    let before = Molecule::new("aspirin", c.clone())
        .with_charges(q.clone())
        .thermalised(300.0, thermostat, 4);
    let after = Molecule::new("aspirin", c)
        .thermalised(300.0, thermostat, 4)
        .with_charges(q);
    let mut bus = pantometry_core::Exchange::new();
    let (mut a, mut b) = (before, after);
    for _ in 0..20 {
        a.step(Time::ZERO, Time::s(0.5 * FS), &mut bus).unwrap();
        b.step(Time::ZERO, Time::s(0.5 * FS), &mut bus).unwrap();
    }
    let bits =
        |m: &Molecule| -> Vec<u64> { (0..21).flat_map(|i| m.at(i)).map(f64::to_bits).collect() };
    assert_eq!(bits(&a), bits(&b));
}

/// The books as `(total, scale)`, the scale as the kernel's audit takes it from the entries:
/// kinetic, the potential counted from `zero`, and minus `work_counted` of the bath's work.
fn books(m: &Molecule, zero: f64, work_counted: f64) -> (f64, f64) {
    let d = m.dynamics().expect("in dynamics");
    let entries = [
        d.kinetic_energy(),
        d.potential_energy().expect("prepared") - zero,
        -work_counted * d.thermostat_work(),
    ];
    let total = entries[0] + entries[1] + entries[2];
    let scale = entries.iter().fold(total.abs(), |s, e| s.max(e.abs()));
    (total, scale)
}

/// The largest relative change of the books across one frame, judged as the kernel's audit judges
/// it — against the largest of the totals and of the entries before and after — for the books as
/// they are kept, and for two that are wrong: the bath's work half counted, with the potential
/// counted from its start as the books count it (`[1]`), and with the potential booked whole, from
/// UFF's own zero (`[2]`), as it was before the books counted from the start.
fn worst_frame_change(frame_fs: f64, picoseconds: f64, seed: u64) -> [f64; 3] {
    let mut m = moving(seed);
    let start = m.dynamics().unwrap().potential_reference().unwrap();
    let mut bus = pantometry_core::Exchange::new();
    let frames = (picoseconds * 1000.0 / frame_fs).round() as usize;
    let mut worst = [0.0f64; 3];
    let ways = [(start, 1.0), (start, 0.5), (0.0, 0.5)];
    for _ in 0..frames {
        let before = ways.map(|(z, w)| books(&m, z, w));
        let mut t = 0.0;
        while t < frame_fs * FS * (1.0 - 1e-9) {
            m.step(Time::ZERO, Time::s(Molecule::DEFAULT_TIME_STEP), &mut bus)
                .unwrap();
            t += Molecule::DEFAULT_TIME_STEP;
        }
        let after = ways.map(|(z, w)| books(&m, z, w));
        for k in 0..3 {
            let scale = before[k].1.max(after[k].1);
            worst[k] = worst[k].max((after[k].0 - before[k].0).abs() / scale);
        }
    }
    // The first is the ledger itself, which `books` must reproduce bit for bit.
    let ledger = m.ledger();
    assert_eq!(
        ledger.get(quantity::ENERGY).unwrap(),
        books(&m, start, 1.0).0
    );
    worst
}

/// What [`Molecule::ENERGY_TOLERANCE`] is earned from: over 20 ps at 300 K and four seeds, for
/// frames of one step, 50 fs and 1 ps, the largest relative change of the books across a frame —
/// and, for the same runs, what the books would show with the bath's work half counted, against
/// the scale the books use now and against the scale with UFF's zero in it.
#[test]
#[ignore = "a measurement that prints a table and asserts nothing: run with --release -- --ignored --nocapture"]
fn the_energy_tolerance_measured() {
    println!("| frame (fs) | books, worst | half work, least worst over seeds | the same, UFF zero in the scale |");
    for frame in [0.5, 50.0, 1000.0] {
        let runs: Vec<[f64; 3]> = (1..=4)
            .map(|seed| worst_frame_change(frame, 20.0, seed))
            .collect();
        let worst = runs.iter().map(|r| r[0]).fold(0.0f64, f64::max);
        let half = runs.iter().map(|r| r[1]).fold(f64::INFINITY, f64::min);
        let whole = runs.iter().map(|r| r[2]).fold(f64::INFINITY, f64::min);
        println!("| {frame} | {worst:.3e} | {half:.3e} | {whole:.3e} |");
    }
}

/// **The dynamics readings are the dynamics' own numbers**, in the units they say: the
/// temperature in °C from the kinetic energy after the bath's `O`, the kinetic energy and the
/// thermostat work in kcal/mol, and the conserved energy their books.
#[test]
fn the_readings_say_how_hot_it_is_and_what_the_bath_did() {
    use pantometry_forcefield::uff::KCAL_PER_MOL;
    let mut m = moving(5);
    let mut bus = pantometry_core::Exchange::new();
    for _ in 0..40 {
        m.step(Time::ZERO, Time::s(0.5 * FS), &mut bus).unwrap();
    }
    let d = m.dynamics().expect("in dynamics");
    let unit = |label: &str| {
        m.readings()
            .into_iter()
            .find(|r| r.label == label)
            .map(|r| r.unit)
            .unwrap()
    };
    assert_eq!(value(&m, "temperature"), d.half_step_temperature() - 273.15);
    assert_eq!(unit("temperature"), "C");
    assert_eq!(
        value(&m, "kinetic energy"),
        d.kinetic_energy() / KCAL_PER_MOL
    );
    assert_eq!(
        value(&m, "thermostat work"),
        d.thermostat_work() / KCAL_PER_MOL
    );
    let books = m.ledger().get(quantity::ENERGY).unwrap();
    // The same three numbers summed in the same order, so the same bits.
    assert_eq!(value(&m, "conserved energy"), books / KCAL_PER_MOL);
    assert_eq!(value(&m, "dynamics steps"), 40.0);
    assert!(m.readings().iter().all(|r| r.label != "converged"));
    // At 300 K, after 20 fs, a molecule of 21 atoms is somewhere between frozen and boiling.
    let t = value(&m, "temperature") + 273.15;
    assert!(t > 100.0 && t < 600.0, "{t} K");
}
