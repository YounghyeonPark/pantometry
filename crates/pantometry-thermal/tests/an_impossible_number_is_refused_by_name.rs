//! A number no body can have is refused by the step that would have used it, by name and with its
//! value, before anything moves.
//!
//! **Most of them ran, and said nothing.** Measured on every thermal domain before the refusal
//! existed, stepped by hand on a block of four cells with one at 100 °C and the rest at 20 °C:
//!
//! - A conductivity that is `NaN`, infinite or negative made a `Solid3D` an **insulator**. The
//!   series mean of a face is guarded at zero, and every one of those reads as zero there, so the
//!   limit was infinite and nothing moved. Through `Simulation` the audit agreed, because nothing
//!   moving conserves.
//! - A negative specific heat or density moved heat up its gradient, to 1e37 K in a block, and gave
//!   a `Bar1D` a negative limit that it stepped under.
//! - An emissivity of 1.5 radiated half as much again as a black body, and one of −0.5 took in what
//!   it sent. A `NaN` one made the field `NaN`, which `Simulation` refused at its second advance as
//!   "energy at simulation", naming neither the block nor the number.
//! - A liquid whose specific heat was negative took a block with no way to lose heat from 100 °C
//!   to −8.6 °C, through `Simulation` and past `Substance::check`. A melting point that was not a
//!   number made every cell `NaN`, and a latent heat of zero or less held a cell in the phase it
//!   began in: ice at −1 °C given eleven kelvin of its own heat sat solid at +10 °C.
//! - A film coefficient or area that was negative lost nothing, and a `NaN` one made the field
//!   `NaN`. A mounting of `NaN` was read as no mounting, and a negative one as a joint carrying
//!   nothing. Air of `NaN` or −10 W·m⁻²·K⁻¹ lost nothing, through `Simulation` too.
//! - A face's or the air's ambient of `NaN` made `Solid3D::steady_state` **panic**, in the clamp
//!   that damps its uniform shift, and an ambient below absolute zero was taken as one.
//! - A `LumpedMass` or a `ThermalNetwork` whose film or emissivity was `NaN` stepped to `NaN`
//!   temperatures, and one whose film, area or emissivity was negative stepped as if it were
//!   smaller. A lumped film without limit gave a limit of zero and a step of `NaN`.
//!
//! **And a limit that is allowed did not work.** A film without limit on a block's face holds the
//! surface at the room, and on a face that radiates it stepped its cell to `NaN`: the surface was
//! solved for as `∞/∞`, while the limit read the half cell and said the step was fine.
//!
//! # What is checked against what
//!
//! **A refusal names the number and carries it.** Each case is asked of `step`, and of
//! `steady_state` where the domain has one, and the message has to name what was wrong and its
//! `detail` has to be the value, bit for bit — a refusal that says "conductivity" for a bad film
//! would send a caller to the wrong place, which is the reason the block's joint refusal names
//! itself. **Wherever the number is**: in an insert rather than the block's own substance, on a
//! second face, on a second node. **And nothing moved**: every temperature is the one it started
//! at, bit for bit, the heat waiting on the bus is still there, and nothing has been lost to a
//! film — on films that act, the exposed cells sitting away from the room — because a refused step
//! that has half happened is a state nobody asked for.
//!
//! **The limits a caller can mean still run, against their closed forms.** A refusal is a bound,
//! and the bound in each is drawn where a body stops being possible rather than where it stops
//! being ordinary: a conductivity of zero, a film, area, mounting or air of zero, an emissivity of
//! zero or one, a film, mounting or air without limit, an ambient of absolute zero. Each is
//! stepped and checked against what it has to do — nothing moving but what was given, nothing
//! lost, one step of `G·(T − T∞)·dt/C` behind a half cell or of `σA(T⁴ − T∞⁴)·dt/C` from a black
//! body, to the ulp, or the steady line `T∞ + P·x/(kA)` of a bar whose surface is held at the
//! ambient.

use pantometry_core::substance::ThermalProps;
use pantometry_core::units::{
    Area, Conductance, Density, LatentHeat, Length, Power, SpecificHeat, Temperature,
    ThermalConductivity, ThermalExpansion, Time, Volume, STEFAN_BOLTZMANN,
};
use pantometry_core::{Domain, Exchange, Schedule, Simulation, Substance, Violation};
use pantometry_thermal::{Bar1D, Environment, Face, LumpedMass, Solid3D, ThermalNetwork, HEAT};

const NAN: f64 = f64::NAN;
const INF: f64 = f64::INFINITY;
/// Joules waiting on the bus when a step is refused, which it must leave there.
const WAITING: f64 = 1e-3;

fn props(k: f64, c: f64, emissivity: f64) -> ThermalProps {
    ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(k),
        specific_heat: SpecificHeat::j_per_kg_k(c),
        expansion: ThermalExpansion::ppm_per_k(50.0),
        emissivity,
    }
}

fn stuff(rho: f64, k: f64, c: f64, emissivity: f64) -> Substance {
    Substance::bulk("stuff", Density::from_si(rho)).with_thermal(props(k, c, emissivity))
}

/// A substance every domain can step: water-like, and radiating, so the emissivity is read.
fn fine() -> Substance {
    stuff(1000.0, 50.0, 4000.0, 0.5)
}

/// Ice whose liquid is `liquid`. The blocks start at 20 °C, so every cell is the liquid.
fn melting(liquid: ThermalProps) -> Substance {
    let ice = Substance::ice();
    let fusion = ice.fusion.expect("ice melts");
    ice.with_fusion(fusion.with_liquid(liquid))
}

/// Ice that melts at `kelvin` and takes `latent` J/kg to.
fn fusing(kelvin: f64, latent: f64) -> Substance {
    let ice = Substance::ice();
    let mut fusion = ice.fusion.expect("ice melts");
    fusion.melting_point = Temperature::from_si(kelvin);
    fusion.latent_heat = LatentHeat::from_si(latent);
    ice.with_fusion(fusion)
}

fn film(h: f64, area: f64, ambient: Temperature) -> Environment {
    Environment {
        ambient,
        convection_w_per_m2_k: h,
        area: Area::from_si(area),
    }
}

/// The room every film and the air sit at: below the blocks' coldest cell, so a film acts on every
/// cell it touches from the first step.
fn cold() -> Temperature {
    Temperature::celsius(0.0)
}

fn still() -> Environment {
    film(10.0, 1e-6, cold())
}

/// Four cells of 1 mm, the first at 100 °C and the rest at 20 °C, so a refusal that let a step
/// through would have something to move.
fn block(s: Substance) -> Solid3D {
    let mut b = Solid3D::new(
        "block",
        s,
        (4, 1, 1),
        Length::mm(1.0),
        Temperature::celsius(20.0),
    );
    b.set_temperature(0, 0, 0, Temperature::celsius(100.0));
    b
}

/// A good block with `s` in its third cell: the mistake in an insert, not in the block.
fn inserted(s: Substance) -> Solid3D {
    let mut b = block(fine());
    b.fill(s, |i, _, _| i == 2);
    b
}

/// The same with a fifth cell, void, for the air to fill.
fn voided(s: Substance) -> Solid3D {
    let mut b = Solid3D::new(
        "block",
        s,
        (5, 1, 1),
        Length::mm(1.0),
        Temperature::celsius(20.0),
    );
    b.set_temperature(0, 0, 0, Temperature::celsius(100.0));
    b.empty(|i, _, _| i == 4)
}

/// A bus with heat waiting on it.
fn waiting() -> Exchange {
    let mut bus = Exchange::new();
    bus.publish(HEAT, WAITING);
    bus
}

/// Every cell's temperature, as bits — a void's `NaN` included, which compares as itself.
fn field(b: &Solid3D) -> Vec<u64> {
    let (nx, ny, nz) = b.counts();
    let mut out = Vec::with_capacity(nx * ny * nz);
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                out.push(b.temperature_at(i, j, k).to_si().to_bits());
            }
        }
    }
    out
}

/// The same value: equal bits, or both not a number.
fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan())
}

/// `refusal` names `names`, carries `value` as its detail, and was raised at `site`.
fn names(why: &str, refusal: &Violation, site: &str, names: &str, value: f64) {
    println!("  {why:34} -> {refusal}");
    assert!(
        refusal.quantity.contains(names),
        "{why}: the refusal must name {names:?}, says {refusal}"
    );
    assert!(
        same(refusal.before, value),
        "{why}: the refusal must carry {value}, carries {}",
        refusal.before
    );
    assert_eq!(refusal.site, site, "{why}: refused at the wrong place");
}

/// What a refused step must have left: the heat on the bus, and nothing lost or taken.
fn untouched(why: &str, bus: &Exchange, lost: f64, absorbed: f64) {
    assert_eq!(
        bus.peek(HEAT),
        WAITING,
        "{why}: a refused step took heat off the bus"
    );
    assert_eq!(lost, 0.0, "{why}: a refused step lost heat");
    assert_eq!(absorbed, 0.0, "{why}: a refused step counted heat taken");
}

/// A step that is long enough to move anything and short enough to be stable where there is a
/// limit at all.
fn a_step(d: &dyn Domain) -> Time {
    let limit = d.max_stable_dt(Time::ZERO).to_si();
    Time::from_si(if limit.is_finite() && limit > 0.0 {
        0.5 * limit
    } else {
        1e-3
    })
}

#[test]
fn a_block_refuses_every_number_no_block_can_have() {
    let below_zero = Temperature::celsius(-300.0);
    let on = |e: Environment| block(fine()).losing_from(Face::XMax, e);
    // A good face first, so the bad one is not the first a loop over faces meets.
    let second = |e: Environment| {
        block(fine())
            .losing_from(Face::XMin, still())
            .losing_from(Face::XMax, e)
    };
    let cases: Vec<(&str, Solid3D, &str, f64)> = vec![
        // The solid's own numbers.
        (
            "conductivity NaN",
            block(stuff(1000.0, NAN, 4000.0, 0.0)),
            "stuff's conductivity",
            NAN,
        ),
        (
            "conductivity +inf",
            block(stuff(1000.0, INF, 4000.0, 0.0)),
            "stuff's conductivity",
            INF,
        ),
        (
            "conductivity -50",
            block(stuff(1000.0, -50.0, 4000.0, 0.0)),
            "stuff's conductivity",
            -50.0,
        ),
        (
            "specific heat NaN",
            block(stuff(1000.0, 50.0, NAN, 0.0)),
            "stuff's specific heat",
            NAN,
        ),
        (
            "specific heat +inf",
            block(stuff(1000.0, 50.0, INF, 0.0)),
            "stuff's specific heat",
            INF,
        ),
        (
            "specific heat 0",
            block(stuff(1000.0, 50.0, 0.0, 0.0)),
            "stuff's specific heat",
            0.0,
        ),
        (
            "specific heat -4000",
            block(stuff(1000.0, 50.0, -4000.0, 0.0)),
            "stuff's specific heat",
            -4000.0,
        ),
        (
            "density NaN",
            block(stuff(NAN, 50.0, 4000.0, 0.0)),
            "stuff's density",
            NAN,
        ),
        (
            "density 0",
            block(stuff(0.0, 50.0, 4000.0, 0.0)),
            "stuff's density",
            0.0,
        ),
        (
            "density -1000",
            block(stuff(-1000.0, 50.0, 4000.0, 0.0)),
            "stuff's density",
            -1000.0,
        ),
        (
            "emissivity NaN",
            block(stuff(1000.0, 50.0, 4000.0, NAN)),
            "stuff's emissivity",
            NAN,
        ),
        (
            "emissivity 1.5",
            block(stuff(1000.0, 50.0, 4000.0, 1.5)),
            "stuff's emissivity",
            1.5,
        ),
        (
            "emissivity -0.5",
            block(stuff(1000.0, 50.0, 4000.0, -0.5)),
            "stuff's emissivity",
            -0.5,
        ),
        // In an insert, which is where a mistyped substance usually is: the block is fine.
        (
            "an insert's conductivity -50",
            inserted(
                Substance::bulk("insert", Density::from_si(1000.0))
                    .with_thermal(props(-50.0, 4000.0, 0.0)),
            ),
            "insert's conductivity",
            -50.0,
        ),
        // The liquid's, which is what every one of these cells is at 20 °C.
        (
            "liquid conductivity NaN",
            block(melting(props(NAN, 4182.0, 0.0))),
            "ice's liquid conductivity",
            NAN,
        ),
        (
            "liquid conductivity +inf",
            block(melting(props(INF, 4182.0, 0.0))),
            "ice's liquid conductivity",
            INF,
        ),
        (
            "liquid conductivity -0.6",
            block(melting(props(-0.6, 4182.0, 0.0))),
            "ice's liquid conductivity",
            -0.6,
        ),
        (
            "liquid specific heat NaN",
            block(melting(props(0.6, NAN, 0.0))),
            "ice's liquid specific heat",
            NAN,
        ),
        (
            "liquid specific heat -4182",
            block(melting(props(0.6, -4182.0, 0.0))),
            "ice's liquid specific heat",
            -4182.0,
        ),
        (
            "liquid emissivity NaN",
            block(melting(props(0.6, 4182.0, NAN))),
            "ice's liquid emissivity",
            NAN,
        ),
        (
            "liquid emissivity 1.5",
            block(melting(props(0.6, 4182.0, 1.5))),
            "ice's liquid emissivity",
            1.5,
        ),
        // Where it melts and what that takes.
        (
            "melting point NaN",
            block(fusing(NAN, 3.34e5)),
            "ice's melting point",
            NAN,
        ),
        (
            "latent heat 0",
            block(fusing(273.15, 0.0)),
            "ice's latent heat",
            0.0,
        ),
        (
            "latent heat -3e5",
            block(fusing(273.15, -3.0e5)),
            "ice's latent heat",
            -3.0e5,
        ),
        // A face's film.
        (
            "film NaN",
            on(film(NAN, 1e-6, cold())),
            "XMax: the film's convection",
            NAN,
        ),
        (
            "film -10",
            on(film(-10.0, 1e-6, cold())),
            "XMax: the film's convection",
            -10.0,
        ),
        (
            "film -10 on the second face",
            second(film(-10.0, 1e-6, cold())),
            "XMax: the film's convection",
            -10.0,
        ),
        (
            "film area NaN",
            on(film(10.0, NAN, cold())),
            "XMax: the film's area",
            NAN,
        ),
        (
            "film area -1e-6",
            on(film(10.0, -1e-6, cold())),
            "XMax: the film's area",
            -1e-6,
        ),
        (
            "film area +inf",
            on(film(10.0, INF, cold())),
            "XMax: the film's area",
            INF,
        ),
        (
            "film ambient NaN",
            on(film(10.0, 1e-6, Temperature::from_si(NAN))),
            "XMax: the film's ambient",
            NAN,
        ),
        (
            "film ambient -300 C",
            on(film(10.0, 1e-6, below_zero)),
            "XMax: the film's ambient",
            below_zero.to_si(),
        ),
        // A mounting on it.
        (
            "mounting NaN",
            on(still()).mounted_on(Face::XMax, NAN),
            "XMax: the mounting",
            NAN,
        ),
        (
            "mounting -100",
            on(still()).mounted_on(Face::XMax, -100.0),
            "XMax: the mounting",
            -100.0,
        ),
        (
            "mounting NaN on the second face",
            second(still())
                .mounted_on(Face::XMin, 100.0)
                .mounted_on(Face::XMax, NAN),
            "XMax: the mounting",
            NAN,
        ),
        // The air in a void.
        (
            "air NaN",
            voided(fine()).air_in(cold(), NAN),
            "the air in the void must be zero",
            NAN,
        ),
        (
            "air -10",
            voided(fine()).air_in(cold(), -10.0),
            "the air in the void must be zero",
            -10.0,
        ),
        (
            "air ambient NaN",
            voided(fine()).air_in(Temperature::from_si(NAN), 10.0),
            "the air in the void must be finite",
            NAN,
        ),
        (
            "air ambient -300 C",
            voided(fine()).air_in(below_zero, 10.0),
            "the air in the void must be finite",
            below_zero.to_si(),
        ),
    ];
    for (why, mut b, named, value) in cases {
        let before = field(&b);
        let dt = a_step(&b);
        let mut bus = waiting();
        let refusal = b
            .step(Time::ZERO, dt, &mut bus)
            .expect_err(&format!("{why} must be refused"));
        names(why, &refusal, "block", named, value);
        assert_eq!(field(&b), before, "{why}: a refused step moved the block");
        untouched(
            why,
            &bus,
            b.lost_energy().to_si(),
            b.absorbed_energy().to_si(),
        );
        // And the balance, which panicked on two of these rather than answering.
        let solve = b
            .steady_state()
            .expect_err(&format!("{why}: the steady state must be refused too"));
        assert_eq!(
            solve.to_string(),
            refusal.to_string(),
            "{why}: the solve must refuse it as the step does"
        );
    }
}

/// **Every substance, including one no cell holds now.** Two substances, the second filled over
/// the whole block: the first was stated, and a mistake in it is refused rather than forgotten.
/// And a substance with no thermal properties at all, filled over in the same way, is the exception
/// the refusal says it is: only the cells that hold one refuse it, and there are none.
#[test]
fn a_block_refuses_a_substance_it_no_longer_holds_and_runs_one_with_nothing_to_say() {
    let mut b = block(stuff(1000.0, -50.0, 4000.0, 0.0));
    b.fill(fine(), |_, _, _| true);
    let mut bus = Exchange::new();
    let refusal = b
        .step(Time::ZERO, a_step(&b), &mut bus)
        .expect_err("a substance that conducts -50 W/m/K was stated");
    names(
        "filled over",
        &refusal,
        "block",
        "stuff's conductivity",
        -50.0,
    );

    let mut b = block(Substance::bulk("unknown", Density::from_si(1000.0)));
    b.fill(fine(), |_, _, _| true);
    let dt = a_step(&b);
    for n in 0..10 {
        b.step(Time::from_si(n as f64 * dt.to_si()), dt, &mut bus)
            .expect("every cell holds a substance that says what it is");
    }
}

/// Through `Simulation`, which let three of these through for twenty advances and refused the
/// fourth at its second, as "energy at simulation": the block's own refusal arrives at the first
/// advance and names the block and the number. The films and the air sit below every cell, so each
/// acts in the first advance whether or not the refusal is there.
#[test]
fn through_the_simulation_the_refusal_is_the_blocks_and_comes_first() {
    let cases: Vec<(&str, Solid3D, &str, f64)> = vec![
        (
            "conductivity NaN",
            block(stuff(1000.0, NAN, 4000.0, 0.0)),
            "stuff's conductivity",
            NAN,
        ),
        (
            "air -10",
            voided(fine()).air_in(cold(), -10.0),
            "the air in the void",
            -10.0,
        ),
        (
            "liquid specific heat -4182",
            block(melting(props(0.6, -4182.0, 0.0))),
            "ice's liquid specific heat",
            -4182.0,
        ),
        (
            "emissivity NaN",
            block(stuff(1000.0, 50.0, 4000.0, NAN)).losing_from(Face::XMax, still()),
            "stuff's emissivity",
            NAN,
        ),
    ];
    for (why, b, named, value) in cases {
        let mut sim = Simulation::new(Schedule::Multirate).with(b);
        let refusal = sim
            .advance(Time::from_si(0.1))
            .expect_err(&format!("{why} must be refused at the first advance"));
        names(why, &refusal, "block", named, value);
    }
}

#[test]
fn a_bar_refuses_every_number_no_bar_can_have() {
    let bar = |s: Substance| {
        Bar1D::new(
            "bar",
            s,
            8,
            Length::mm(1.0),
            Area::from_si(1e-6),
            Temperature::celsius(20.0),
        )
    };
    let cases: Vec<(&str, Bar1D, &str, f64)> = vec![
        (
            "conductivity NaN",
            bar(stuff(1000.0, NAN, 4000.0, 0.0)),
            "stuff's conductivity",
            NAN,
        ),
        (
            "conductivity -50",
            bar(stuff(1000.0, -50.0, 4000.0, 0.0)),
            "stuff's conductivity",
            -50.0,
        ),
        (
            "specific heat +inf",
            bar(stuff(1000.0, 50.0, INF, 0.0)),
            "stuff's specific heat",
            INF,
        ),
        (
            "specific heat -4000",
            bar(stuff(1000.0, 50.0, -4000.0, 0.0)),
            "stuff's specific heat",
            -4000.0,
        ),
        (
            "density -1000",
            bar(stuff(-1000.0, 50.0, 4000.0, 0.0)),
            "stuff's density",
            -1000.0,
        ),
        // A bar does not radiate, and the substance is wrong all the same.
        (
            "emissivity 1.5",
            bar(stuff(1000.0, 50.0, 4000.0, 1.5)),
            "stuff's emissivity",
            1.5,
        ),
    ];
    for (why, mut b, named, value) in cases {
        let before: Vec<u64> = (0..b.cell_count())
            .map(|i| b.temperature_at(i).to_si().to_bits())
            .collect();
        // Heat waiting on the bus, which a step that had begun would have taken.
        let mut bus = waiting();
        let refusal = b
            .step(Time::ZERO, Time::from_si(1e-3), &mut bus)
            .expect_err(&format!("{why} must be refused"));
        names(why, &refusal, "bar", named, value);
        let after: Vec<u64> = (0..b.cell_count())
            .map(|i| b.temperature_at(i).to_si().to_bits())
            .collect();
        assert_eq!(after, before, "{why}: a refused step moved the bar");
        untouched(why, &bus, 0.0, b.absorbed_energy().to_si());
    }
}

fn lump(s: Substance, e: Environment) -> LumpedMass {
    LumpedMass::new(
        "lump",
        s,
        Volume::from_si(1e-6),
        Length::mm(1.0),
        Temperature::celsius(100.0),
        e,
    )
}

#[test]
fn a_lumped_mass_refuses_every_number_no_body_can_have() {
    let below_zero = Temperature::celsius(-300.0);
    let at = |h: f64, area: f64| film(h, area, cold());
    let cases: Vec<(&str, LumpedMass, &str, f64)> = vec![
        (
            "specific heat -4000",
            lump(stuff(1000.0, 50.0, -4000.0, 0.5), still()),
            "stuff's specific heat",
            -4000.0,
        ),
        (
            "emissivity NaN",
            lump(stuff(1000.0, 50.0, 4000.0, NAN), still()),
            "stuff's emissivity",
            NAN,
        ),
        (
            "emissivity 1.5",
            lump(stuff(1000.0, 50.0, 4000.0, 1.5), still()),
            "stuff's emissivity",
            1.5,
        ),
        (
            "emissivity -0.5",
            lump(stuff(1000.0, 50.0, 4000.0, -0.5), still()),
            "stuff's emissivity",
            -0.5,
        ),
        (
            "film NaN",
            lump(fine(), at(NAN, 1e-6)),
            "the film's convection",
            NAN,
        ),
        (
            "film -10",
            lump(fine(), at(-10.0, 1e-6)),
            "the film's convection",
            -10.0,
        ),
        // Nothing in series with a lumped film, so a film without limit is a limit of zero.
        (
            "film +inf",
            lump(fine(), at(INF, 1e-6)),
            "the film's convection",
            INF,
        ),
        (
            "area NaN",
            lump(fine(), at(10.0, NAN)),
            "the film's area",
            NAN,
        ),
        (
            "area -1e-6",
            lump(fine(), at(10.0, -1e-6)),
            "the film's area",
            -1e-6,
        ),
        (
            "ambient NaN",
            lump(fine(), film(10.0, 1e-6, Temperature::from_si(NAN))),
            "the film's ambient",
            NAN,
        ),
        (
            "ambient -300 C",
            lump(fine(), film(10.0, 1e-6, below_zero)),
            "the film's ambient",
            below_zero.to_si(),
        ),
    ];
    for (why, mut m, named, value) in cases {
        let before = m.temperature().to_si().to_bits();
        let mut bus = waiting();
        let refusal = m
            .step(Time::ZERO, Time::from_si(1e-3), &mut bus)
            .expect_err(&format!("{why} must be refused"));
        names(why, &refusal, "lump", named, value);
        assert_eq!(
            m.temperature().to_si().to_bits(),
            before,
            "{why}: a refused step moved the body"
        );
        untouched(
            why,
            &bus,
            m.lost_energy().to_si(),
            m.absorbed_energy().to_si(),
        );
    }
}

/// A network of two nodes joined by a fine link, the second absorbing: the first made of `a` and
/// losing to `ea`, the second made of `b` and losing to `eb`.
fn two_nodes(a: Substance, ea: Environment, b: Substance, eb: Environment) -> ThermalNetwork {
    let mut net = ThermalNetwork::new("net");
    let first = net.node_losing_to(
        "a",
        a,
        Volume::from_si(1e-6),
        Length::mm(1.0),
        Temperature::celsius(100.0),
        ea,
    );
    let second = net.node_losing_to(
        "b",
        b,
        Volume::from_si(1e-6),
        Length::mm(1.0),
        Temperature::celsius(20.0),
        eb,
    );
    net.link(first, second, Conductance::from_si(0.01))
        .expect("a fine link");
    net.absorbing(second).expect("its own node");
    net
}

/// The mistake on the first node.
fn network(s: Substance, e: Environment) -> ThermalNetwork {
    two_nodes(s, e, fine(), still())
}

#[test]
fn a_network_refuses_every_number_no_node_can_have() {
    let at = |h: f64, area: f64| film(h, area, cold());
    let cases: Vec<(&str, ThermalNetwork, &str, &str, f64)> = vec![
        (
            "specific heat -4000",
            network(stuff(1000.0, 50.0, -4000.0, 0.5), still()),
            "net/a",
            "stuff's specific heat",
            -4000.0,
        ),
        (
            "emissivity NaN",
            network(stuff(1000.0, 50.0, 4000.0, NAN), still()),
            "net/a",
            "stuff's emissivity",
            NAN,
        ),
        (
            "emissivity 1.5",
            network(stuff(1000.0, 50.0, 4000.0, 1.5), still()),
            "net/a",
            "stuff's emissivity",
            1.5,
        ),
        (
            "film NaN",
            network(fine(), at(NAN, 1e-6)),
            "net/a",
            "the film's convection",
            NAN,
        ),
        (
            "film -10",
            network(fine(), at(-10.0, 1e-6)),
            "net/a",
            "the film's convection",
            -10.0,
        ),
        (
            "film +inf",
            network(fine(), at(INF, 1e-6)),
            "net/a",
            "the film's convection",
            INF,
        ),
        (
            "area -1e-6",
            network(fine(), at(10.0, -1e-6)),
            "net/a",
            "the film's area",
            -1e-6,
        ),
        (
            "ambient NaN",
            network(fine(), film(10.0, 1e-6, Temperature::from_si(NAN))),
            "net/a",
            "the film's ambient",
            NAN,
        ),
        // The second node, so a loop that stops at the first does not pass for one over all.
        (
            "the second node's film -10",
            two_nodes(fine(), still(), fine(), at(-10.0, 1e-6)),
            "net/b",
            "the film's convection",
            -10.0,
        ),
    ];
    for (why, mut net, site, named, value) in cases {
        let nodes: Vec<_> = net.handles().map(|(n, _)| n).collect();
        let before: Vec<u64> = nodes
            .iter()
            .map(|n| net.temperature(*n).to_si().to_bits())
            .collect();
        let mut bus = waiting();
        let refusal = net
            .step(Time::ZERO, Time::from_si(1e-3), &mut bus)
            .expect_err(&format!("{why} must be refused"));
        names(why, &refusal, site, named, value);
        let after: Vec<u64> = nodes
            .iter()
            .map(|n| net.temperature(*n).to_si().to_bits())
            .collect();
        assert_eq!(after, before, "{why}: a refused step moved the network");
        untouched(
            why,
            &bus,
            net.lost_energy().to_si(),
            net.absorbed_energy().to_si(),
        );
        let solve = match net.steady_state(Power::w(0.01)) {
            Ok(_) => panic!("{why}: the steady state must be refused too"),
            Err(v) => v,
        };
        assert_eq!(
            solve.to_string(),
            refusal.to_string(),
            "{why}: the solve must refuse it as the step does"
        );
    }
}

/// **`steady_state` answers rather than panicking, whatever the field.** A cell given a temperature
/// of `NaN` reached the clamp that damps the uniform shift with every number the block was built
/// from possible, and `clamp` panics on a bound that is not a number; then it was refused as a
/// solve that diverged. It is refused as the step refuses it now, by the cell, before any sweep.
///
/// **And a field the solve itself takes past a number still diverges**, rather than panicking: a
/// source of 1e300 W, finite and a number, overflows the field it drives.
#[test]
fn a_field_that_is_not_a_number_is_refused_rather_than_panicking() {
    let mut b = block(fine()).losing_from(Face::XMax, still());
    b.set_temperature(1, 0, 0, Temperature::from_si(NAN));
    let refusal = b
        .steady_state()
        .expect_err("a field with a NaN in it has no balance");
    names(
        "a NaN cell",
        &refusal,
        "block",
        "cell (1, 0, 0)'s temperature must be finite",
        NAN,
    );
    let b = block(fine())
        .losing_from(Face::XMax, still())
        .dissipating(1e300, |i, _, _| i == 0);
    let refusal = b
        .steady_state()
        .expect_err("a field driven past a number has no balance");
    println!("  a source of 1e300 W -> {refusal}");
    assert!(
        refusal.quantity.contains("diverged"),
        "it is a divergence, says {refusal}"
    );
}

/// The steady state of a bar of `n` cells along x whose surface at `x = 0` is held at `ambient`
/// and which generates `p` watts in its last cell: the whole of `p` crosses every face between
/// that cell and the surface, so the line `T∞ + p·x/(kA)`, with `x` from the surface to a cell's
/// centre, is exact for the discrete operator as for the continuum — the half cell and every face
/// are each `kA/dx` per cell length.
fn held_line(n: usize, p: f64, k: f64, dx: f64, ambient: f64) -> Vec<f64> {
    (0..n)
        .map(|i| ambient + p * (i as f64 + 0.5) * dx / (k * dx * dx))
        .collect()
}

/// Steps `b` `steps` times at half its limit, refusing nothing: a limit a caller can mean.
fn march(why: &str, b: &mut dyn Domain, steps: usize, bus: &mut Exchange) {
    let dt = a_step(&*b);
    for n in 0..steps {
        b.step(Time::from_si(n as f64 * dt.to_si()), dt, bus)
            .unwrap_or_else(|v| panic!("{why}: a limit a caller can mean: {v}"));
    }
}

#[test]
fn the_limits_a_caller_can_mean_still_run() {
    let (k, dx, ambient, p) = (50.0, 1e-3, 300.0, 0.01);
    let mut bus = Exchange::new();

    // **A conductivity of zero** is a block that conducts nothing: every cell where it started.
    let mut b = block(stuff(1000.0, 0.0, 4000.0, 0.0));
    let before = field(&b);
    march("an insulating block", &mut b, 50, &mut bus);
    assert_eq!(field(&b), before, "an insulator moved");
    // And a bar, given heat: it lands in the first cell and goes nowhere.
    let area = 1e-6;
    let mut bar = Bar1D::new(
        "bar",
        stuff(1000.0, 0.0, 4000.0, 0.0),
        8,
        Length::from_si(dx),
        Area::from_si(area),
        Temperature::celsius(20.0),
    );
    let mut fed = waiting();
    march("an insulating bar", &mut bar, 50, &mut fed);
    let start = Temperature::celsius(20.0).to_si();
    let first = start + WAITING / ((1000.0 * (area * dx)) * 4000.0);
    let got = bar.temperature_at(0).to_si();
    assert!(
        (got - first).abs() <= 2.0 * f64::EPSILON * first,
        "the insulating bar's first cell is at {got} K, and the heat it was given puts it at {first} K"
    );
    for i in 1..8 {
        assert_eq!(
            bar.temperature_at(i).to_si(),
            start,
            "an insulating bar moved heat into cell {i}"
        );
    }

    // **A film, an area, a mounting or air of zero** loses nothing, exactly — the substance does
    // not radiate, so the film is the only path and it is closed. Every exposed cell sits 20 K
    // above the room, so a film that were open would lose from the first step.
    let dark = || stuff(1000.0, 50.0, 4000.0, 0.0);
    let closed: Vec<(&str, Solid3D)> = vec![
        (
            "film 0",
            block(dark()).losing_from(Face::XMax, film(0.0, 1e-6, cold())),
        ),
        (
            "area 0",
            block(dark()).losing_from(Face::XMax, film(10.0, 0.0, cold())),
        ),
        (
            "mounting 0",
            block(dark())
                .losing_from(Face::XMax, still())
                .mounted_on(Face::XMax, 0.0),
        ),
        ("air 0", voided(dark()).air_in(cold(), 0.0)),
    ];
    for (why, mut b) in closed {
        march(why, &mut b, 50, &mut bus);
        assert_eq!(b.lost_energy().to_si(), 0.0, "{why} lost heat");
    }

    // **A mounting without limit** is the face bolted straight to its film, bit for bit — on a
    // film that loses, so the two cannot agree by both losing nothing or both being `NaN`.
    let mut bolted = block(fine()).losing_from(Face::XMax, still());
    let mut mounted = block(fine())
        .losing_from(Face::XMax, still())
        .mounted_on(Face::XMax, INF);
    march("bolted", &mut bolted, 50, &mut bus);
    march("a mounting without limit", &mut mounted, 50, &mut bus);
    let lost = bolted.lost_energy().to_si();
    assert!(
        lost.is_finite() && lost > 0.0,
        "the bolted face lost {lost} J"
    );
    assert_eq!(
        field(&mounted),
        field(&bolted),
        "a mounting without limit is none"
    );

    // **Air at absolute zero** takes heat, and says nothing is wrong.
    let mut frozen = voided(dark()).air_in(Temperature::from_si(0.0), 10.0);
    march("air at absolute zero", &mut frozen, 50, &mut bus);
    let lost = frozen.lost_energy().to_si();
    assert!(
        lost.is_finite() && lost > 0.0,
        "air at absolute zero took {lost} J"
    );

    // **A film without limit on a face that radiates**: one cell, whose surface the film holds at
    // the room, so what leaves it is the half cell's `G = 2kA/dx` times its excess, and its limit
    // is `C/G`. One step at half of it is `T − G·(T − T∞)·dt/C`, to the ulp. It stepped to `NaN`.
    let (t0, room, eps) = (400.0, 300.0, 0.5);
    let mut cell = Solid3D::new(
        "a radiating face held at the room",
        stuff(1000.0, k, 4000.0, eps),
        (1, 1, 1),
        Length::from_si(dx),
        Temperature::from_si(t0),
    )
    .losing_from(Face::XMin, film(INF, dx * dx, Temperature::from_si(room)));
    let g = 2.0 * k * dx * dx / dx;
    let c = 1000.0 * (dx * dx * dx) * 4000.0;
    let limit = cell.max_stable_dt(Time::ZERO).to_si();
    assert!(
        (limit / (c / g) - 1.0).abs() < 1e-14,
        "the limit is C/G = {} s, and the block reports {limit} s",
        c / g
    );
    let dt = 0.5 * limit;
    cell.step(Time::ZERO, Time::from_si(dt), &mut bus)
        .unwrap_or_else(|v| panic!("a radiating face held at the room: {v}"));
    let expected = t0 - g * (t0 - room) * dt / c;
    let got = cell.temperature_at(0, 0, 0).to_si();
    println!("  a radiating face held at the room: {got:.12} K against {expected:.12} K");
    assert!(
        (got - expected).abs() <= 4.0 * f64::EPSILON * t0,
        "one step of a radiating face held at the room is at {got} K, and should be at \
         {expected} K"
    );

    // **A film without limit, and air without limit**, hold a surface at their ambient, so the
    // steady state is the held line. Eight cells; for the air a ninth, void, at the cold end. The
    // film's block radiates, which a surface at the room makes no difference to. Each is marched
    // first, which is the path a scene takes.
    let n = 8;
    let line = held_line(n, p, k, dx, ambient);
    let held = Solid3D::new(
        "held by a film",
        stuff(1000.0, k, 4000.0, 0.5),
        (n, 1, 1),
        Length::from_si(dx),
        Temperature::from_si(ambient),
    )
    .losing_from(
        Face::XMin,
        film(INF, dx * dx, Temperature::from_si(ambient)),
    )
    .dissipating(p, |i, _, _| i == n - 1);
    let aired = Solid3D::new(
        "held by the air",
        stuff(1000.0, k, 4000.0, 0.0),
        (n + 1, 1, 1),
        Length::from_si(dx),
        Temperature::from_si(ambient),
    )
    .empty(|i, _, _| i == 0)
    .air_in(Temperature::from_si(ambient), INF)
    .dissipating(p, |i, _, _| i == n);
    for (b, first) in [(&held, 0), (&aired, 1)] {
        let name = b.name().to_string();
        // A copy is marched, so the solve below starts at the room as the measurements did.
        let mut marched = b.clone();
        march(&name, &mut marched, 20, &mut bus);
        assert!(
            (first..first + n).all(|i| marched.temperature_at(i, 0, 0).to_si().is_finite()),
            "{name}: a cell is not a number after the march"
        );
        let t = b.steady_state().unwrap_or_else(|v| panic!("{name}: {v}"));
        assert!(
            t[first..].iter().all(|v| v.is_finite()),
            "{name}: the steady state has a cell that is not a number"
        );
        let worst = (0..n)
            .map(|i| (t[first + i] - line[i]).abs())
            .fold(0.0, f64::max);
        println!(
            "  {name:16} the held line {:.6} .. {:.6} K above, worst {worst:.3e} K off",
            line[0] - ambient,
            line[n - 1] - ambient
        );
        // **The solver's alone.** The line is exact for the operator, so what is left is the
        // solve's stopping rule: it stops when no round moves a cell by 1e-6 K, and what that
        // leaves is the slow mode's remaining decay. A round here is an over-relaxed sweep — at
        // `ω = 2/(1 + sin(π/m))` for the longest axis, 1.4465 for the film's eight cells and
        // 1.4903 for the air's nine — and then the uniform shift, and that iteration contracts at
        // `ρ = 0.73` and 0.70, from its own matrix, by an emulation in `f64` that reproduces both
        // measured errors to four figures. Measured, 3.6e-6 K and 2.8e-6 K: 4.4 and 3.8 times the
        // last round's move, rather than the `ρ/(1 − ρ)` of 2.7 and 2.3, because `ρ` is a rate and
        // not a norm. So 1e-5 K: ten times the stopping rule and 2.8 times the larger error. A film
        // read at the centre rather than the surface would put every cell 0.1 K off, and a film
        // without limit read as none would have no balance at all.
        assert!(worst < 1e-5, "{name}: {worst:e} K off the held line");
    }

    // **An emissivity of one and an ambient of absolute zero**, on a lumped body with no
    // convection: one step is `T − σA(T⁴ − T∞⁴)·dt/C`. To two ulps of the start, which a step
    // that moves 1.7e-7 K resolves to nine figures; measured, to the bit.
    let (area, volume, t0, dt) = (1e-6, 1e-6, Temperature::celsius(100.0).to_si(), 1e-3);
    for (why, emissivity, sky) in [
        ("black, to a room", 1.0, Temperature::celsius(20.0).to_si()),
        ("black, to absolute zero", 1.0, 0.0),
    ] {
        let mut m = lump(
            stuff(1000.0, 50.0, 4000.0, emissivity),
            film(0.0, area, Temperature::from_si(sky)),
        );
        m.step(Time::ZERO, Time::from_si(dt), &mut bus)
            .unwrap_or_else(|v| panic!("{why}: {v}"));
        let fourth = |t: f64| (t * t) * (t * t);
        let capacity = 1000.0 * 4000.0 * volume;
        let expected = t0
            - emissivity * STEFAN_BOLTZMANN.to_si() * area * (fourth(t0) - fourth(sky)) * dt
                / capacity;
        let got = m.temperature().to_si();
        println!("  {why:24} {got:.12} K against {expected:.12} K");
        assert!(
            (got - expected).abs() <= 2.0 * f64::EPSILON * t0,
            "{why}: {got} K against {expected} K"
        );
    }
    // And one that does not radiate, which does not move at all.
    let mut m = lump(
        stuff(1000.0, 50.0, 4000.0, 0.0),
        film(0.0, area, Temperature::celsius(20.0)),
    );
    m.step(Time::ZERO, Time::from_si(dt), &mut bus)
        .unwrap_or_else(|v| panic!("not radiating: {v}"));
    assert_eq!(m.temperature().to_si(), t0, "a body that cannot lose moved");

    // **A network** with a film and an area of zero loses nothing either, its nodes 100 K and
    // 20 K above the room.
    for (why, e) in [
        ("film 0", film(0.0, 1e-6, cold())),
        ("area 0", film(10.0, 0.0, cold())),
    ] {
        let mut net = two_nodes(dark(), e, dark(), e);
        march(why, &mut net, 50, &mut bus);
        assert_eq!(net.lost_energy().to_si(), 0.0, "{why} lost heat");
    }
}
