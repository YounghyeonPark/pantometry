//! The two ways a cooled face could be handed a step it could not take.
//!
//! Both were found by review, both were reproduced before they were fixed, and both had the
//! same signature: the block does not diverge loudly, it **overshoots the air it is losing to**
//! while the conservation audit stays perfectly happy — because the heat really did leave.

use pantometry_core::conserved::quantity;
use pantometry_core::units::{Area, Length, Temperature, Time};
use pantometry_core::{Domain, Exchange, Substance};
use pantometry_thermal::{Environment, Face, Solid3D, HEAT};

/// **A cell colder than its surroundings is not handed a step sized for the wrong curve.**
///
/// The radiative derivative `4εσT³` dominates the secant only for `T ≥ T∞`; below ambient it
/// reverses without bound. Measured before the fix: a 20 mm cell at 20 °C inside a 700 °C
/// radiant enclosure had a true step ratio of **13**, and one accepted step took it to 8847 °C.
#[test]
fn a_cell_colder_than_its_surroundings_is_not_overdriven() {
    let glass = Substance::from_name("borosilicate").expect("the catalogue has it");
    let mut block = Solid3D::new(
        "part",
        glass,
        (1, 1, 1),
        Length::mm(20.0),
        Temperature::celsius(20.0),
    )
    .losing_from(
        Face::ZMax,
        Environment {
            // A radiant enclosure far hotter than the part: the loss runs the other way.
            ambient: Temperature::celsius(700.0),
            convection_w_per_m2_k: 0.0,
            area: Area::from_si(20.0e-3 * 20.0e-3),
        },
    );

    let limit = block.max_stable_dt(Time::ZERO);
    let mut bus = Exchange::new();
    block
        .step(Time::ZERO, limit, &mut bus)
        .expect("the reported limit must be a step this domain can take");

    let after = block.temperature_at(0, 0, 0).to_si();
    let ambient = Temperature::celsius(700.0).to_si();
    assert!(
        after > Temperature::celsius(20.0).to_si() - 1e-9 && after <= ambient,
        "a cell warming toward its enclosure must not pass it: {after} K against {ambient} K"
    );
}

/// **Heat off the bus arrives after the sweep, so a slug is added to the step and none of it is
/// shed by the step that brings it.**
///
/// `stability_ratio` reads the temperatures the step *started* with. Bus heat was spread over
/// every cell before the film was applied, so a slug could make a cell far stiffer than the limit
/// was sized for. Measured on this cube when this was written, at ambient and given its rise in one
/// accepted step: 300 K took a corner to −137 °C and 900 K to −7007 °C, with the ledger balancing
/// perfectly throughout.
///
/// It was refused for that, by a re-check in the film's own pass after the sweep — and the
/// re-check went when the film moved into the sweep. This test allowed a refusal *or* a step that
/// left no cell below the air, and by then the cube no longer reached the air, so it passed on its
/// second branch with no refusal anywhere — its coldest cell 506 K above the air — while one black
/// cell of 1 mm and `k = 200` at 1000 K, stepped at its own limit under a slug of 1000 K, landed at
/// −1982 K. The heat goes in after the sweep now — see
/// `tests/heat_off_the_bus_arrives_after_the_sweep.rs` — and nothing is refused: every cell lands
/// where the step without the slug puts it, plus the rise, to the rounding of adding 900 K.
#[test]
fn a_bus_slug_is_added_to_every_cell_of_its_step() {
    let pla = Substance::from_name("pla").expect("the catalogue has it");
    let side = 4.0 * 5.0e-3;
    // Hot, so the step without the slug does something: its films shed, its corners most. From
    // the air's temperature it was the identity, and the check came down to the rise alone.
    let build = || {
        let mut block = Solid3D::new(
            "part",
            pla.clone(),
            (4, 4, 4),
            Length::mm(5.0),
            Temperature::celsius(200.0),
        );
        for face in Face::ALL {
            block = block.losing_from(
                face,
                Environment {
                    ambient: Temperature::celsius(20.0),
                    convection_w_per_m2_k: 0.0,
                    area: Area::from_si(side * side),
                },
            );
        }
        block
    };

    // A slug big enough that the film, met at the post-slug temperature, would overshoot.
    let block = build();
    let limit = block.max_stable_dt(Time::ZERO);
    // From the substance rather than from the block, so a block that spread the heat over the
    // wrong total would not also be the one saying how much heat that is.
    let capacity = pla.density.to_si()
        * pla
            .thermal
            .expect("the catalogue gives it thermal properties")
            .specific_heat
            .to_si()
        * side
        * side
        * side;
    let rise = 900.0;

    let mut bare = block.clone();
    bare.step(Time::ZERO, limit, &mut Exchange::new())
        .expect("its own limit is a step it can take");
    let mut slugged = block.clone();
    let mut bus = Exchange::new();
    bus.publish(HEAT, rise * capacity);
    slugged
        .step(Time::ZERO, limit, &mut bus)
        .expect("a slug is added to the step, and the step is the one the limit was taken for");

    let ambient = Temperature::celsius(20.0).to_si();
    for k in 0..4 {
        for j in 0..4 {
            for i in 0..4 {
                let (without, with) = (
                    bare.temperature_at(i, j, k).to_si(),
                    slugged.temperature_at(i, j, k).to_si(),
                );
                // 900 K added to a cell near 473 K and read back: an ulp or two of 1400 K.
                assert!(
                    (with - (without + rise)).abs() < 1e-11,
                    "cell ({i},{j},{k}) lands at {with} K, where the step without the slug puts \
                     it at {without} K and the slug raises it {rise} K"
                );
                assert!(
                    with >= ambient,
                    "cell ({i},{j},{k}) was shed to {with} K, below the {ambient} K air"
                );
            }
        }
    }
}

/// **The books still balance through steps at the limit.** The audit could never see either defect
/// above — the heat genuinely left — so this is not the check that catches them; it is the check
/// that the fixes did not break the thing that was right. It was named for a refusal that the bus
/// slug once met, and that nothing meets now.
#[test]
fn the_books_balance_through_steps_at_the_limit() {
    let glass = Substance::from_name("borosilicate").expect("the catalogue has it");
    let mut block = Solid3D::new(
        "part",
        glass,
        (2, 2, 2),
        Length::mm(10.0),
        Temperature::celsius(400.0),
    )
    .losing_from(
        Face::ZMax,
        Environment::still_air(Temperature::celsius(20.0), Area::from_si(4.0e-4)),
    );

    let books = |b: &Solid3D| {
        b.ledger()
            .get(quantity::ENERGY)
            .expect("a block keeps its books in energy")
    };
    let opening = books(&block);
    let dt = block.max_stable_dt(Time::ZERO);
    // A slug off the bus in the first step, which goes in after its sweep: what crosses the bus is
    // all the books may move by, whatever the film shed meanwhile.
    let slug = 100.0;
    for step in 0..200 {
        let mut bus = Exchange::new();
        if step == 0 {
            bus.publish(HEAT, slug);
        }
        block.step(Time::ZERO, dt, &mut bus).expect("a stable step");
    }
    let closing = books(&block);
    let moved = block.lost_energy().to_si();
    assert!(moved > 0.0);
    assert!(
        (closing - opening - slug).abs() / moved.max(slug) < 1e-12,
        "the ledger moved by {:e} J, where {slug} J came over the bus and {moved:e} J was shed",
        closing - opening
    );
}
