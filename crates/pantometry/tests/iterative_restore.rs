//! The restore path of `Schedule::Iterative`, which nothing in this workspace had ever taken.
//!
//! `supports_restore` is a promise every restorable domain makes, and `iterate` is the only
//! caller that collects on it — by rewinding every domain before each sweep after the first.
//! But no coupling in this workspace is strong enough to need iterating, so `Domain::residual`
//! is the default `0.0` everywhere, `iterate` converges on its first sweep, and the restore
//! branch has never executed. A design review noticed that and it turned out to be the more
//! useful half of what it was looking for.
//!
//! So this file supplies the missing ingredient: a domain with a genuine residual, which forces
//! real iterations, which forces real restores.

use pantometry::core::integrator::substeps_for;
use pantometry::prelude::*;

/// A domain that converges over sweeps, so `iterate` has to iterate.
///
/// The point is the split between what is rewound and what is not. `x` is state and is
/// checkpointed; `estimate` is what the sweep *learned* and is deliberately not, because that is
/// exactly what makes an iterative coupling converge — each sweep restarts from the same
/// physical state with a better guess about its neighbours. A domain that rewound its estimate
/// too would repeat the same sweep forever.
struct Converging {
    /// Rewound by `restore`.
    x: f64,
    /// Carried across restores, on purpose.
    estimate: f64,
    last_change: f64,
    saved: Option<f64>,
}

impl Converging {
    fn new() -> Converging {
        Converging {
            x: 0.0,
            estimate: 0.0,
            last_change: f64::INFINITY,
            saved: None,
        }
    }
}

impl Domain for Converging {
    fn name(&self) -> &str {
        "converging"
    }

    /// `e ← (e + 1)/2`, whose fixed point is 1 and whose error halves every sweep. Chosen
    /// because the number of iterations to a given tolerance is then exactly predictable, so
    /// the test can assert it rather than observe it.
    fn step(&mut self, _t: Time, _dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        let next = 0.5 * (self.estimate + 1.0);
        self.last_change = (next - self.estimate).abs();
        self.estimate = next;
        self.x = next;
        Ok(())
    }

    fn residual(&self) -> f64 {
        self.last_change
    }

    /// Also where the estimate is reset, because `iterate` calls this exactly once per
    /// `advance` — so it is the only hook a domain has for "a new outer step is starting".
    ///
    /// Without it the estimate stays converged from the previous step, every later `advance`
    /// finishes in one sweep, and the restore path goes unexercised again. Which is a small
    /// version of how it came to be unexercised in the first place.
    fn checkpoint(&mut self) {
        self.saved = Some(self.x);
        self.estimate = 0.0;
        self.last_change = f64::INFINITY;
    }

    fn restore(&mut self) {
        if let Some(x) = self.saved {
            self.x = x;
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }
}

fn plate() -> LumpedMass {
    LumpedMass::new(
        "plate",
        Substance::aluminium_6061(),
        Volume::from_si(60e-3 * 60e-3 * 3e-3),
        Length::mm(1.5),
        Temperature::celsius(80.0),
        Environment::still_air(
            Temperature::celsius(20.0),
            Area::from_si(2.0 * 60e-3 * 60e-3),
        ),
    )
}

/// **The restore branch executes, and it is asserted that it does.**
///
/// Starting from zero, `e ← (e+1)/2` changes by 0.5, then 0.25, then 0.125, so a tolerance of
/// 0.2 is met on the third sweep and `iterate` restores twice on the way. Before this file the
/// same call converged on sweep one and restored zero times, which is why `supports_restore`
/// had never been collected on.
#[test]
fn an_iterative_coupling_actually_iterates_and_actually_restores() {
    let mut sim = Simulation::new(Schedule::Iterative {
        max_iter: 8,
        tol: 0.2,
    })
    .with(Converging::new());

    let report = sim.advance(Time::s(1.0)).expect("it converges");
    assert_eq!(
        report.iterations, 3,
        "0.5, 0.25, 0.125 against a tolerance of 0.2 is three sweeps"
    );

    // And a coupling that cannot reach the tolerance is refused rather than accepted, which is
    // the other half of the promise.
    let mut stubborn = Simulation::new(Schedule::Iterative {
        max_iter: 2,
        tol: 1e-9,
    })
    .with(Converging::new());
    let violation = stubborn
        .advance(Time::s(1.0))
        .expect_err("two sweeps cannot reach 1e-9");
    assert_eq!(violation.quantity, "coupling residual");
    assert_eq!(
        stubborn.time().to_si(),
        0.0,
        "a refused step keeps the clock"
    );
}

/// **A `LumpedMass` rewound mid-coupling does not report energy it did not lose.**
///
/// The hypothesis this file was written to test. `LumpedMass::checkpoint` saves the temperature
/// and nothing else, while `ledger()` reports `stored + lost` — so if `lost` kept accumulating
/// across sweeps whose temperatures were rewound, the books would grow by one sweep's losses
/// every iteration and the audit would see energy created out of nothing.
///
/// It is not hypothetical any more: with `Converging` in the simulation this runs three sweeps
/// and two restores per `advance`, and the conservation audit is live at 1e-9 throughout.
#[test]
fn a_rewound_lumped_mass_does_not_invent_the_heat_it_shed() {
    let mut sim = Simulation::new(Schedule::Iterative {
        max_iter: 8,
        tol: 0.2,
    })
    .conservation_tolerance(1e-9)
    .with(Converging::new())
    .with(plate());

    // The plate starts at 80 C against 20 C ambient, so it sheds heat every substep of every
    // sweep and `lost` is large. That is load-bearing: the first version of this test started
    // it *at* ambient, where `lost` stays identically zero — and a double count of zero is
    // zero, so it passed while measuring nothing. The configuration-with-no-physics trap, in
    // the test written to catch a different one.
    let start = sim
        .domain_as::<LumpedMass>("plate")
        .expect("the plate is there")
        .lost_energy()
        .to_si();
    assert_eq!(start, 0.0);
    for _ in 0..40 {
        let report = sim
            .advance(Time::s(5.0))
            .expect("a rewound sweep must not create energy");
        assert_eq!(report.iterations, 3);
    }

    // And it really did shed something, so the audit had something to be wrong about: three
    // sweeps an advance for forty advances, every one of them rewinding the temperature while
    // this total kept climbing.
    let shed = sim
        .domain_as::<LumpedMass>("plate")
        .expect("the plate is there")
        .lost_energy()
        .to_si();
    assert!(
        shed > 100.0,
        "the plate should have shed real joules, got {shed:.3} J"
    );
}

/// **A room and a hall that audit their own books can be iterated.**
///
/// Both claim `books_balance`, so every sweep compares the energy each holds at the start of its
/// turn with what it holds after and what crossed the bus. Their restores put back neither of two
/// things only the energy reads. The startup offset, which the first step sets, while an advance's
/// first checkpoint is taken before it: the second sweep of the first window was refused, at the
/// default 1e-9, a room reporting 1.1987e-5 J against the 1.1959e-5 it held. And the pressure a step
/// before, which the energy is a product with: with the offset put back, the second window was
/// refused instead. Each claimed `supports_restore` and neither could be iterated past one sweep.
///
/// So each runs beside `Converging` for three windows of three sweeps, and ends every window on the
/// readings and the books one sweep a window gives under `Schedule::Multirate`, to the bit. More
/// than one, because the first is checkpointed before any step, where the pressure a step before is
/// the pressure now, and a window of half a limit is one substep, which leaves it there: leaving it
/// out of the restore costs nothing until the second window.
#[test]
fn a_room_that_audits_its_own_books_can_be_iterated() {
    let released = |name: &str| -> Box<dyn Domain> {
        if name == "room" {
            Box::new(
                Room::of_air("room", Length::m(4.4), Length::m(3.1), 21).released_in_mode(
                    1,
                    1,
                    Pressure::from_si(1.0),
                ),
            )
        } else {
            Box::new(
                Hall::of_air("hall", Length::m(4.4), Length::m(3.1), Length::m(2.4), 9)
                    .released_in_mode((1, 1, 1), Pressure::from_si(1.0)),
            )
        }
    };
    let view = |sim: &Simulation, name: &str| -> Vec<u64> {
        let d = sim.domain(name).expect("it is there");
        let mut bits: Vec<u64> = d.readings().iter().map(|r| r.value.to_bits()).collect();
        bits.extend(d.ledger().quantities().map(|(_, v)| v.to_bits()));
        bits
    };
    for name in ["room", "hall"] {
        let (alone, coupled) = (released(name), released(name));
        assert!(
            alone.books_balance(),
            "{name}, precondition: it audits its own books"
        );
        let window = alone.max_stable_dt(Time::ZERO) * 0.5;
        let mut once = Simulation::new(Schedule::Multirate).with_boxed(alone);
        let mut sim = Simulation::new(Schedule::Iterative {
            max_iter: 8,
            tol: 0.2,
        })
        .with(Converging::new())
        .with_boxed(coupled);
        for advance in 1..=3 {
            once.advance(window).expect("one window");
            let report = sim
                .advance(window)
                .unwrap_or_else(|v| panic!("{name}, advance {advance}: refused: {v:?}"));
            assert_eq!(
                report.iterations, 3,
                "{name}, advance {advance}: three sweeps, so two restores"
            );
            assert!(
                view(&sim, name) == view(&once, name),
                "{name}, advance {advance}: three sweeps and two restores did not end where one \
                 sweep does"
            );
        }
    }
}

/// A warming espresso bed: 8 x 6 x 8 cells of 1 mm at 20 C under a 93 C inlet, so a step moves its
/// viscosity, its flow and its stability limit, which falls as it warms.
fn warming_puck() -> Puck {
    let dx = 1e-3;
    let mut p = Puck::new(
        "puck",
        Basket {
            counts: (8, 6, 8),
            cell: Length::from_si(dx),
            radius: Length::from_si(4.0 * dx),
            porosity: 0.45,
            grind: Grind::sieved(Length::from_si(800e-6)),
            ..Basket::espresso()
        },
    );
    p.set_temperature(Temperature::celsius(20.0));
    p.set_inlet_temperature(Temperature::celsius(93.0));
    p
}

/// What a bed reads out, as bits: each cell's temperature, extraction, concentration and pressure,
/// the flow, the delivered totals and both books. Not its `Debug`, which the porous tests also
/// compare: a bed in a coupling holds the checkpoint the coupling took, and a bed alone holds none.
fn bed_bits(sim: &Simulation) -> Vec<u64> {
    let p = sim.domain_as::<Puck>("puck").expect("the puck is there");
    let (nx, ny, nz) = p.counts();
    let mut bits = Vec::new();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                bits.push(p.temperature_at(i, j, k).to_si().to_bits());
                bits.push(p.extraction_at(i, j, k).to_bits());
                bits.push(p.concentration_at(i, j, k).to_si().to_bits());
                bits.push(p.pressure_at(i, j, k).to_si().to_bits());
            }
        }
    }
    let books = p.ledger();
    for x in [
        p.flow_rate().to_si(),
        p.delivered().to_si(),
        p.delivered_solute().to_si(),
        p.delivered_enthalpy().to_si(),
        books.get(quantity::MASS).expect("a mass book"),
        books.get(quantity::ENERGY).expect("an energy book"),
    ] {
        bits.push(x.to_bits());
    }
    bits
}

/// **A puck rewound by an iterative coupling moves one window an advance, however many sweeps it
/// takes.**
///
/// `Converging` makes `iterate` sweep three times and restore twice an advance. The puck is not
/// coupled to it, so every sweep is the same window from the same bed, and the bed it ends on must
/// be the one that window gives under `Schedule::Multirate`, which sweeps once and restores
/// nothing — to the bit, and divided into the same substeps. The window is two and a half of the
/// bed's first limit, and the limit falls as the bed warms, so each window is divided again inside
/// itself, into more substeps than its first read plans, and how depends on the bed the sweep
/// began from. Two advances, so the second iterates from a checkpoint taken over the first's: one
/// that kept the first would rewind the second advance to where the first began.
///
/// The puck's restore took its checkpoint as it went, so the second restore did nothing and the
/// third sweep began where the second ended: two windows in one. It did not claim
/// `supports_restore` then, and `iterate` refused it.
#[test]
fn a_puck_rewound_by_an_iterative_coupling_moves_one_window() {
    let t = Time::from_si(0.0);
    let puck = warming_puck();
    let window = Time::from_si(2.5 * puck.max_stable_dt(t).to_si());

    let mut once = Simulation::new(Schedule::Multirate).with(puck.clone());
    let mut sim = Simulation::new(Schedule::Iterative {
        max_iter: 8,
        tol: 0.2,
    })
    .with(Converging::new())
    .with(puck);
    let substeps = |r: &Report| {
        r.substeps
            .iter()
            .find(|(name, _)| name == "puck")
            .map(|(_, n)| *n)
    };
    for advance in 1..=2 {
        let limit = once
            .domain_as::<Puck>("puck")
            .expect("the puck is there")
            .max_stable_dt(t);
        let planned = substeps_for(window, limit);
        let alone = once.advance(window).expect("one window");
        let coupled = sim
            .advance(window)
            .expect("a restorable puck in a coupling that converges");
        println!(
            "  advance {advance}: {planned} substeps planned, {:?} alone, {:?} in the coupling",
            substeps(&alone),
            substeps(&coupled)
        );
        assert_eq!(
            coupled.iterations, 3,
            "advance {advance}: three sweeps, so two restores"
        );
        assert!(
            substeps(&alone) > Some(planned),
            "advance {advance}, precondition: the window must be divided again inside itself"
        );
        assert_eq!(
            substeps(&coupled),
            substeps(&alone),
            "advance {advance}: the third sweep divided the window differently from one sweep \
             from the same bed"
        );
        assert!(
            bed_bits(&sim) == bed_bits(&once),
            "advance {advance}: the bed after three sweeps and two restores is not the bed one \
             window gives"
        );
    }
}
