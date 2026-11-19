//! The audit passed a ledger that had stopped being a number.
//!
//! `audit` compared `|after − before| / scale` against the tolerance, and every comparison with
//! `NaN` is false, so a total that became `NaN` was never over the line. The scale was no help:
//! `f64::max` discards a `NaN` operand, so a ledger that was `NaN` at both ends had a scale of
//! zero and was skipped as "two denormals". A run whose books were poisoned at its first step
//! audited clean to its last.
//!
//! `pantometry-world`'s verification battery knew, said so in its own comments, and recorded it
//! from a copy of the audit's arithmetic. A run that is not verified never ran that copy.

use pantometry_core::conserved::quantity;
use pantometry_core::units::{Area, Time};
use pantometry_core::{
    audit, Domain, Exchange, Flux, Interface, Kind, Ledger, Schedule, Simulation, Violation,
};

fn energy(v: f64) -> Ledger {
    Ledger::new().with(quantity::ENERGY, v)
}

/// **Every way a total can stop being finite is refused**, and the message says which.
///
/// The five rows are the ways the old arithmetic passed one: a finite total that became `NaN`
/// (relative change `NaN`, compared false), one that was `NaN` at both ends (scale zero,
/// skipped), one that came back from `NaN` (the same), one that overflowed to an infinity
/// (relative change `NaN` again, since `∞ / ∞` is), and one infinite at both ends (`∞ − ∞` is
/// `NaN`) — the shape a source with an unlimited reserve gave its ledger.
///
/// At a tolerance above zero. At exactly zero a refusal of `NaN` to `NaN` and a
/// [`Violation::at`] detail of `NaN` are the same value — both store the number twice beside a
/// zero tolerance — and the message reads as the detail. The run is still refused; telling the
/// two apart needs a field `Violation` does not have, and is recorded rather than done here.
#[test]
fn a_total_that_is_not_finite_is_refused() {
    let cases = [
        (1.0, f64::NAN),
        (f64::NAN, f64::NAN),
        (f64::NAN, 1.0),
        (1.0, f64::INFINITY),
        (f64::NEG_INFINITY, f64::NEG_INFINITY),
    ];
    for (before, after) in cases {
        let err = audit("poisoned", &energy(before), &energy(after), 1e-9)
            .err()
            .unwrap_or_else(|| panic!("{before} to {after} passed the audit"));
        assert_eq!(err.quantity, quantity::ENERGY);
        assert_eq!(err.site, "poisoned");
        let text = err.to_string();
        assert!(
            text.contains("not finite"),
            "{before} to {after} should say the total is not finite: {text}"
        );
        // "Destroyed" and "created" are claims about which way an amount moved, and neither
        // means anything for an amount that is not one.
        assert!(
            !text.contains("destroyed") && !text.contains("created"),
            "{before} to {after}: {text}"
        );
    }
}

/// **The finite cases are unchanged**: the line is where it was, both sides of it.
#[test]
fn a_finite_total_is_judged_as_before() {
    assert!(audit("fine", &energy(3.7), &energy(3.7 + 4e-16), 1e-12).is_ok());
    let err = audit("leak", &energy(1.0), &energy(0.6), 1e-9).expect_err("40% is a leak");
    assert!(err.to_string().contains("destroyed"), "{err}");
    // Two denormals are still equal for every purpose a simulation has.
    assert!(audit("tiny", &energy(1e-310), &energy(2e-310), 1e-9).is_ok());
}

/// **An amount left on the bus that is not a number is refused too.** `audit_transfers` compared
/// `|left| > tolerance`, which a `NaN` never is — on a plain channel and, in a loop of its own, on
/// every face of a spatial one, so both are put to it. The spatial one is poisoned on one face
/// of three, because it is audited face by face and a face is what the refusal must name.
#[test]
fn a_channel_left_holding_nan_is_refused() {
    let mut bus = Exchange::new();
    bus.publish(quantity::ENERGY, f64::NAN);
    let err = bus
        .audit_transfers("bus", 1e-9)
        .expect_err("NaN published and never consumed is not nothing");
    assert_eq!(err.quantity, quantity::ENERGY);
    assert!(err.to_string().contains("not finite"), "{err}");

    let skin = Interface::uniform("skin", 3, Area::from_si(1e-4));
    let mut spatial = Exchange::new();
    spatial
        .publish_on(
            &skin,
            quantity::ENERGY,
            &Flux::from_faces(vec![0.0, f64::NAN, 0.0]),
        )
        .expect("three faces onto three");
    let err = spatial
        .audit_transfers("bus", 1e-9)
        .expect_err("NaN on a face published and never consumed is not nothing");
    assert_eq!(err.quantity, "skin/energy face 1");
    assert!(err.to_string().contains("not finite"), "{err}");

    let mut fine = Exchange::new();
    fine.publish(quantity::ENERGY, 0.0);
    fine.publish_on(&skin, quantity::ENERGY, &Flux::zeros(3))
        .expect("three faces onto three");
    assert!(fine.audit_transfers("bus", 1e-9).is_ok());
}

/// A domain whose stored energy becomes `NaN` at a chosen step.
struct Poisoned {
    energy: f64,
    steps: u32,
    at: u32,
}

impl Domain for Poisoned {
    fn name(&self) -> &str {
        "poisoned"
    }
    fn kind(&self) -> Kind {
        Kind::Evolving
    }
    fn step(&mut self, _t: Time, _dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        self.steps += 1;
        if self.steps == self.at {
            self.energy = f64::NAN;
        }
        Ok(())
    }
    fn ledger(&self) -> Ledger {
        energy(self.energy)
    }
}

/// **The run stops at the step that poisoned it**, not at the end and not never.
#[test]
fn a_run_stops_at_the_step_its_books_stop_being_numbers() {
    let mut sim = Simulation::new(Schedule::Staggered).with(Poisoned {
        energy: 1000.0,
        steps: 0,
        at: 3,
    });
    let dt = Time::from_si(1.0);
    sim.advance(dt).expect("step 1 is finite");
    sim.advance(dt).expect("step 2 is finite");
    let err = sim
        .advance(dt)
        .expect_err("step 3 made the books NaN, and the audit must say so there");
    assert_eq!(err.quantity, quantity::ENERGY);
    assert!(err.after.is_nan(), "{err:?}");
}

/// **A refusal that carries a `NaN` detail reads as that refusal**, not as a conservation law.
///
/// `Violation::at` stores its detail as both `before` and `after` and its `Display` recognised
/// that by `before == after` — which is false for `NaN`. A pressure solve that never converged
/// was reported as `pressure residual destroyed at puck: NaN became NaN`, which names a law
/// nobody stated and a direction nothing moved in.
#[test]
fn a_detail_that_is_not_a_number_is_still_a_detail() {
    let v = Violation::at("puck", "pressure residual", f64::NAN);
    let text = v.to_string();
    assert_eq!(text, "at puck: pressure residual (NaN)");
}
