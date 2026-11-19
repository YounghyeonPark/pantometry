//! Electricity publishing heat and thermal consuming it, over a bus neither crate names.
//!
//! Lives in the facade rather than in either domain, because a domain that dev-depended on
//! another would break the property the crate split exists to hold — and `invariant-guard`
//! greps the manifest, so a dev-dependency is caught the same as a real one. The README says
//! this is where the cross-domain integration tests live; this is the sixth pair.

use pantometry::prelude::*;
use pantometry_electrical::Winding;

/// **A winding warms a plate, and the joules that arrive are `I²R·t` to the bit.**
///
/// The coupling that motivated the crate: electrical publishes, thermal consumes, over the
/// channel that already existed, with neither crate naming the other. Run under a live
/// conservation audit at 1e-9.
///
/// The rise is checked against `Q/C` with the capacity computed here from the substance, so
/// this is not the simulation compared with itself.
#[test]
fn a_winding_warms_a_plate_by_the_joules_it_dissipated() {
    let volume = Volume::from_si(60e-3 * 60e-3 * 3e-3);
    let start = Temperature::celsius(20.0);
    let seconds = 5.0;

    let coil = Winding::of_resistance("coil", Resistance::ohm(2.0), 0.0, start)
        .driven_at(Current::a(1.5))
        .with_reserve(1_000.0);
    let watts = coil.dissipation().to_si();
    assert!(
        (watts - 4.5).abs() < 1e-12,
        "I^2 R = 1.5^2 * 2 = 4.5, got {watts}"
    );

    let plate = LumpedMass::new(
        "plate",
        Substance::aluminium_6061(),
        volume,
        Length::mm(1.5),
        start,
        // Ambient *at* the starting temperature and no radiation would still leave convection
        // shedding heat as it warms, so this instead gives it a vanishing area: the plate keeps
        // every joule, and the closed form below is exact rather than approximate.
        Environment::still_air(start, Area::from_si(0.0)),
    );

    let mut sim = Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-9)
        .with(coil)
        .with(plate);
    let dt = 0.01;
    for _ in 0..(seconds / dt) as usize {
        sim.advance(Time::s(dt)).expect("the books close");
    }

    let capacity = Substance::aluminium_6061()
        .heat_capacity(volume)
        .expect("aluminium has a specific heat")
        .to_si();
    let want = watts * seconds / capacity;
    let got = sim
        .domain_as::<LumpedMass>("plate")
        .expect("the plate is there")
        .temperature()
        .to_si()
        - start.to_si();
    assert!(
        (got / want - 1.0).abs() < 1e-9,
        "rise {got:.9} K against {want:.9} K"
    );

    // The winding's books say the same thing from the other side.
    let coil = sim.domain_as::<Winding>("coil").expect("the coil is there");
    assert!((coil.dissipated_energy().to_si() - watts * seconds).abs() < 1e-9);
    assert!((coil.reserve().to_si() - (1_000.0 - watts * seconds)).abs() < 1e-9);
}

/// **A source with an unlimited supply is refused by the domain, and by name.**
///
/// This test was written expecting the *audit* to catch it and it did not: the ledger reported
/// `inf` before and `inf` after, `inf` equals itself, and a winding pouring joules into a plate
/// ran green. An infinite reserve is not a large number; it is a number that cannot be
/// subtracted, so it switched the check off rather than failing it. The audit now refuses a total
/// that is not finite, but a winding with no reserve reports no entry rather than an infinite
/// one — see the next test for why — so the refusal is still the domain's.
///
/// The same shape as `pantometry-world`'s lamp, where a skipped `with_reserve` left a scene
/// auditing clean at tolerance zero with the lamp doing nothing. So `Winding::step` refuses
/// while the reserve is infinite, and this pins that it is the domain saying so and by name.
#[test]
fn a_winding_with_no_reserve_is_caught_by_the_audit() {
    let coil = Winding::of_resistance(
        "coil",
        Resistance::ohm(2.0),
        0.0,
        Temperature::celsius(20.0),
    )
    .driven_at(Current::a(1.5));
    let plate = LumpedMass::new(
        "plate",
        Substance::aluminium_6061(),
        Volume::from_si(1e-4),
        Length::mm(1.5),
        Temperature::celsius(20.0),
        Environment::still_air(Temperature::celsius(20.0), Area::from_si(0.0)),
    );
    let mut sim = Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-9)
        .with(coil)
        .with(plate);

    let violation = sim
        .advance(Time::s(1.0))
        .expect_err("joules from an infinite tank are joules from nowhere");
    assert_eq!(violation.site, "coil", "the winding names itself");
    assert!(
        violation.quantity.contains("with_reserve"),
        "{}",
        violation.quantity
    );
    assert_eq!(sim.time().to_si(), 0.0, "a refused step keeps the clock");
}

/// A store of joules that loses one a step to nowhere: the leak an audit exists to catch.
struct Leaky {
    held: f64,
}

impl Domain for Leaky {
    fn name(&self) -> &str {
        "leaky"
    }
    fn kind(&self) -> Kind {
        Kind::Evolving
    }
    fn step(&mut self, _t: Time, _dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        self.held -= 1.0;
        Ok(())
    }
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.held)
    }
}

/// **A winding nobody drives, with no reserve set, neither stops a run nor hides a leak in it.**
///
/// `Winding::step` deliberately demands a reserve only of a winding that would publish, so one at
/// zero current may leave it infinite. Its ledger reported that `inf`, and two things followed
/// from one entry. While the audit compared `inf` with itself, the energy total of the *whole*
/// simulation was `inf`, so a domain beside it destroying a joule a step went unseen. Once the
/// audit refused a total that is not finite, the same winding stopped every run it was in, at its
/// first step, for books that were not a number. It now reports no entry: it publishes nothing,
/// so there is nothing to account for.
///
/// So beside a plate the run is clean, and beside a leak the audit names the leak.
#[test]
fn an_undriven_winding_with_no_reserve_leaves_the_audit_to_the_rest() {
    let idle = || {
        Winding::of_resistance(
            "coil",
            Resistance::ohm(2.0),
            0.0,
            Temperature::celsius(20.0),
        )
    };
    let start = Temperature::celsius(20.0);
    let plate = LumpedMass::new(
        "plate",
        Substance::aluminium_6061(),
        Volume::from_si(1e-4),
        Length::mm(1.5),
        start,
        Environment::still_air(start, Area::from_si(0.0)),
    );
    let mut sim = Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-9)
        .with(idle())
        .with(plate);
    for n in 0..10 {
        sim.advance(Time::s(0.1))
            .unwrap_or_else(|e| panic!("step {n} beside an idle winding was refused: {e}"));
    }

    let mut sim = Simulation::new(Schedule::Staggered)
        .conservation_tolerance(1e-9)
        .with(idle())
        .with(Leaky { held: 1000.0 });
    let violation = sim
        .advance(Time::s(0.1))
        .expect_err("a joule a step from nowhere is a leak");
    let said = violation.to_string();
    assert_eq!(violation.quantity, quantity::ENERGY, "{said}");
    assert!(
        violation.before.is_finite() && violation.after.is_finite(),
        "the leak is judged as numbers, not refused as not being any: {said}"
    );
    assert!(violation.after < violation.before, "{said}");
}
