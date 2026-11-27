//! A temperature below absolute zero is refused: one a caller set, before a step reads it, and one a
//! step would make, with the body, its books and the bus as the step found them.
//!
//! **Each of these answered `Ok`.** Measured before this change:
//!
//! - **A state set below zero.** A lumped body, a bar, a block and a network started at −10 K each
//!   stepped, the lumped body to 5.43 K and the bar holding −10 K, and started at `NaN` each stepped
//!   to `NaN`. A block given −50 K by `set_temperature` stepped on from it, and one given a
//!   kilojoule taken out by `deposit` stepped with a cell at −206 386 K.
//! - **A step that goes below zero.** A megajoule off 100 cm³ of water at 20 °C put a lumped body
//!   and a network node at −2206.85 K; a kilojoule off a bar's first cell put it at −3840.45 K; a
//!   block whose sink of 50 mW outran the 2 mW its film could bring in stepped twenty thousand
//!   times to −354.67 K. A lumped body stepped one and a half time constants from 1000 K overshot
//!   its room to −90.3 K. Water at 1 °C between two cells of ice at 1 K, stepped at its block's
//!   limit, froze to −108.4 K.
//! - **A balance below zero.** 100 W taken out of a node a still-air room can bring 20.5 W to was
//!   answered at −1135.42 K, and its path conductance there as 0.0700 W/K; with a surface that
//!   radiates, the solve went off and was refused as not finite, and with one that only radiates,
//!   as singular. A block's balance with that sink was refused as not converging.
//!
//! # What is checked against what
//!
//! **A refusal names the cell or the node and carries the temperature**, bit for bit: the one set,
//! or the one the step would have made, which is the closed form `T + (Q − L)/C` of a body nothing
//! else moves, `L` its film's loss over the step. **Nothing moved**: every temperature, every
//! melted fraction, the heat waiting on the bus and on the faces and what the bus has counted as
//! taken, and every running total — from totals that are not zero — are as the step found them, so
//! the refused step can be asked again, and is refused again. **A part in a billion short of the
//! bound runs**, to the closed form, and a part in a billion past it is refused. Every question
//! asked of a body set below zero refuses with its step's refusal.
//!
//! **Only what can take a cell below zero is asked**, and the rest is held to the reason: a block
//! with no sink and no phase change, in a room at absolute zero, cooling through films and a
//! clearance at its limit — near zero, and from 1000 K where radiation carries it — never goes
//! below it. Nor does a liquid freezing at its block's limit, which the limit bounds and a pass over
//! every cell of a two-phase block used to refuse. A balance below zero is refused against its
//! closed form, `Tₐ + (P + εσA·Tₐ⁴)/hA` for a node alone, to the solve's own error.

use pantometry_core::substance::ThermalProps;
use pantometry_core::units::{
    Area, Conductance, Density, Energy, Length, Power, SpecificHeat, Temperature,
    ThermalConductivity, ThermalExpansion, Time, Volume, STEFAN_BOLTZMANN,
};
use pantometry_core::{Domain, Exchange, Flux, Substance, Violation};
use pantometry_thermal::{Bar1D, Environment, Face, LumpedMass, Solid3D, ThermalNetwork, HEAT};

const NAN: f64 = f64::NAN;
const INF: f64 = f64::INFINITY;

fn stuff(k: f64, emissivity: f64) -> Substance {
    Substance::bulk("stuff", Density::from_si(1000.0)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(k),
        specific_heat: SpecificHeat::j_per_kg_k(4000.0),
        expansion: ThermalExpansion::ppm_per_k(50.0),
        emissivity,
    })
}

/// Water-like, conducting, and radiating, so the emissivity is read.
fn fine() -> Substance {
    stuff(50.0, 0.5)
}

/// Ice with water named as its liquid, which makes a block of it two-phase; and the liquid holding
/// `c_l` J/kg/K, which for water is twice what ice holds.
fn two_phase_ice(c_l: f64) -> Substance {
    let base = Substance::ice();
    let fusion = base.fusion.expect("ice melts");
    base.with_fusion(fusion.with_liquid(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(0.6),
        specific_heat: SpecificHeat::j_per_kg_k(c_l),
        expansion: ThermalExpansion::ppm_per_k(207.0),
        emissivity: 0.96,
    }))
}

fn film(h: f64, area: f64, ambient: Temperature) -> Environment {
    Environment {
        ambient,
        convection_w_per_m2_k: h,
        area: Area::from_si(area),
    }
}

/// A surface that loses nothing, so a step moves a body by exactly the heat it is given.
fn sealed() -> Environment {
    film(0.0, 0.0, Temperature::celsius(20.0))
}

fn still() -> Environment {
    film(10.0, 1e-6, Temperature::celsius(0.0))
}

fn waiting(joules: f64) -> Exchange {
    let mut bus = Exchange::new();
    bus.publish(HEAT, joules);
    bus
}

fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan())
}

/// `refusal` names `names`, carries `value`, and was raised at `site`.
fn names(why: &str, refusal: &Violation, site: &str, names: &str, value: f64) {
    println!("  {why:56} -> {refusal}");
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

fn identical(why: &str, asked: &str, got: &Violation, step: &Violation) {
    assert!(
        got.quantity == step.quantity && got.site == step.site && same(got.before, step.before),
        "{why}: {asked} refused with {got}, where the step refused with {step}"
    );
}

const BELOW: &str = "below absolute zero, to this many kelvin";
const SET: &str = "must be finite and at or above absolute zero, in kelvin";

/// What a bus holds and has counted as taken, as bits: a refused step that took heat and put it
/// back would leave the first and not the second.
fn bus_of(bus: &Exchange) -> [u64; 2] {
    [bus.peek(HEAT).to_bits(), bus.total_consumed(HEAT).to_bits()]
}

fn lump(initial: Temperature, environment: Environment) -> LumpedMass {
    LumpedMass::new(
        "lump",
        fine(),
        Volume::from_si(1e-6),
        Length::mm(1.0),
        initial,
        environment,
    )
}

/// Eight cells of 1 mm with a section of 1 mm², exposed one face per cell.
fn bar(initial: Temperature) -> Bar1D {
    Bar1D::new(
        "bar",
        fine(),
        8,
        Length::mm(1.0),
        Area::from_si(1e-6),
        initial,
    )
    .exposing("side", Area::from_si(1e-6))
}

/// Four cells of 1 mm at 20 °C, the first at 100 °C, under a film on the far face.
fn block() -> Solid3D {
    let mut b = Solid3D::new(
        "block",
        fine(),
        (4, 1, 1),
        Length::mm(1.0),
        Temperature::celsius(20.0),
    )
    .losing_from(Face::XMax, still());
    b.set_temperature(0, 0, 0, Temperature::celsius(100.0));
    b
}

/// `b` stepped once with nothing on the bus, so what is set after it is the only thing set since
/// a step last read the cells, and its books are not zero.
fn stepped(mut b: Solid3D) -> Solid3D {
    let dt = a_step(&b);
    b.step(Time::ZERO, dt, &mut Exchange::new())
        .expect("a fine block steps");
    b
}

/// Two nodes joined by a link: the first interior and absorbing, the second under a film and
/// starting at `second`.
fn network(second: Temperature) -> ThermalNetwork {
    let mut net = ThermalNetwork::new("net");
    let a = net.node(
        "a",
        fine(),
        Volume::from_si(1e-6),
        Length::mm(1.0),
        Temperature::celsius(100.0),
    );
    let b = net.node_losing_to(
        "b",
        fine(),
        Volume::from_si(1e-6),
        Length::mm(1.0),
        second,
        still(),
    );
    net.link(a, b, Conductance::from_si(0.01))
        .expect("a fine link");
    net.absorbing(a).expect("its own node");
    net
}

fn cells_of_bar(b: &Bar1D) -> Vec<u64> {
    (0..b.cell_count())
        .map(|i| b.temperature_at(i).to_si().to_bits())
        .collect()
}

fn cells_of_block(b: &Solid3D) -> Vec<u64> {
    let (nx, ny, nz) = b.counts();
    let mut out = Vec::new();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                out.push(b.temperature_at(i, j, k).to_si().to_bits());
                out.push(b.melted_fraction_at(i, j, k).to_bits());
            }
        }
    }
    out
}

fn books_of_block(b: &Solid3D) -> [u64; 3] {
    [
        b.absorbed_energy().to_si().to_bits(),
        b.lost_energy().to_si().to_bits(),
        b.generated_energy().to_si().to_bits(),
    ]
}

fn a_step(d: &dyn Domain) -> Time {
    let limit = d.max_stable_dt(Time::ZERO).to_si();
    Time::from_si(if limit.is_finite() && limit > 0.0 {
        0.5 * limit
    } else {
        1e-3
    })
}

/// Step a block that must be refused, check that nothing moved — cells, phases, books and bus —
/// and that it is refused again, as it was.
fn refused_block(why: &str, b: &mut Solid3D, bus: &mut Exchange, dt: Time) -> Violation {
    let (cells, books, heat) = (cells_of_block(b), books_of_block(b), bus_of(bus));
    let refusal = b
        .step(Time::ZERO, dt, bus)
        .expect_err(&format!("{why} must be refused"));
    assert_eq!(
        cells_of_block(b),
        cells,
        "{why}: a refused step moved the block"
    );
    assert_eq!(
        books_of_block(b),
        books,
        "{why}: a refused step moved the books"
    );
    assert_eq!(bus_of(bus), heat, "{why}: a refused step took heat");
    let again = b
        .step(Time::ZERO, dt, bus)
        .expect_err(&format!("{why} must be refused again"));
    identical(why, "a second step", &again, &refusal);
    refusal
}

#[test]
fn a_temperature_set_below_absolute_zero_is_refused_before_a_step_reads_it() {
    for kelvin in [-10.0, NAN, INF, -1e-4] {
        let at = Temperature::from_si(kelvin);
        let why = format!("a lumped body started at {kelvin} K");
        let mut m = lump(at, still());
        let mut bus = waiting(1e-3);
        let refusal = m
            .step(Time::ZERO, Time::from_si(1e-3), &mut bus)
            .expect_err(&format!("{why} must be refused"));
        names(&why, &refusal, "lump", "the body's temperature", kelvin);
        assert!(refusal.quantity.contains(SET), "{why}: {refusal}");
        assert!(same(m.temperature().to_si(), kelvin), "{why}: moved");
        assert_eq!(bus_of(&bus), bus_of(&waiting(1e-3)), "{why}: took heat");
        identical(&why, "biot_number", &m.biot_number().unwrap_err(), &refusal);
        let tau = m.time_constant().unwrap_err();
        identical(&why, "time_constant", &tau, &refusal);
        let rise = m.equilibrium_rise(Power::w(1.0)).unwrap_err();
        identical(&why, "equilibrium_rise", &rise, &refusal);

        let why = format!("a bar started at {kelvin} K");
        let mut b = bar(at);
        let mut bus = waiting(1e-3);
        let refusal = b
            .step(Time::ZERO, Time::from_si(1e-3), &mut bus)
            .expect_err(&format!("{why} must be refused"));
        names(&why, &refusal, "bar", "cell 0's temperature", kelvin);
        assert_eq!(bus_of(&bus), bus_of(&waiting(1e-3)), "{why}: took heat");

        let why = format!("a block started at {kelvin} K");
        let mut b = Solid3D::new("block", fine(), (4, 1, 1), Length::mm(1.0), at)
            .losing_from(Face::XMax, still());
        let refusal = refused_block(&why, &mut b, &mut waiting(1e-3), Time::from_si(1e-3));
        names(
            &why,
            &refusal,
            "block",
            "cell (0, 0, 0)'s temperature",
            kelvin,
        );
        let balance = b.steady_state().unwrap_err();
        identical(&why, "steady_state", &balance, &refusal);

        // On the second node, which a check of the first would miss.
        let why = format!("a network's second node started at {kelvin} K");
        let mut net = network(at);
        let refusal = net
            .step(Time::ZERO, Time::from_si(1e-3), &mut waiting(1e-3))
            .expect_err(&format!("{why} must be refused"));
        names(&why, &refusal, "net/b", "the node's temperature", kelvin);
        let a = net.node_named("a").expect("the first node");
        let got = net.steady_state(Power::w(1.0)).unwrap_err();
        identical(&why, "steady_state", &got, &refusal);
        let tau = net.time_constant(a).unwrap_err();
        identical(&why, "time_constant", &tau, &refusal);
        identical(
            &why,
            "biot_number",
            &net.biot_number(a).unwrap_err(),
            &refusal,
        );
    }

    // Within a microkelvin below zero is the rounding a step at its limit can leave, and runs.
    let hair = Temperature::from_si(-5e-7);
    lump(hair, still())
        .step(Time::ZERO, Time::from_si(1e-3), &mut Exchange::new())
        .expect("a lumped body a hair below zero runs");
    bar(hair)
        .step(Time::ZERO, Time::from_si(1e-3), &mut Exchange::new())
        .expect("a bar a hair below zero runs");
    network(hair)
        .step(Time::ZERO, Time::from_si(1e-3), &mut Exchange::new())
        .expect("a network a hair below zero runs");

    // Set on a block between steps, in a cell other than the first, every way a caller can, each
    // after a step that read the cells.
    let dt = a_step(&block());
    let mut set = stepped(block());
    set.set_temperature(2, 0, 0, Temperature::from_si(-50.0));
    let refusal = refused_block("set to -50 K", &mut set, &mut waiting(1e-3), dt);
    names(
        "set to -50 K",
        &refusal,
        "block",
        "cell (2, 0, 0)'s temperature",
        -50.0,
    );

    let mut given = stepped(block());
    given.deposit(1, 0, 0, Energy::from_si(-2.0));
    let would = given.temperature_at(1, 0, 0).to_si();
    assert!(
        would < 0.0,
        "two joules did not take the cell below zero: {would} K"
    );
    let why = "two joules taken out of a cell by deposit";
    let refusal = refused_block(why, &mut given, &mut waiting(1e-3), dt);
    names(
        why,
        &refusal,
        "block",
        "cell (1, 0, 0)'s temperature",
        would,
    );

    let mut moded = stepped(block());
    moded.release_mode((1, 0, 0), Temperature::from_si(10.0), 20.0);
    let why = "a mode whose amplitude is more than its mean";
    let would = moded.temperature_at(3, 0, 0).to_si();
    let refusal = refused_block(why, &mut moded, &mut waiting(1e-3), dt);
    names(
        why,
        &refusal,
        "block",
        "cell (3, 0, 0)'s temperature",
        would,
    );

    // A void cell holds a placeholder: below zero it is not asked, and the block runs with its
    // neighbour as it would; filled after a step, it is a temperature again, and is asked.
    let mut hollow = block().empty(|i, _, _| i == 3);
    hollow.set_temperature(3, 0, 0, Temperature::from_si(-50.0));
    let mut hollow = stepped(hollow);
    let beside = hollow.temperature_at(2, 0, 0).to_si();
    assert!(
        beside.is_finite() && beside > 0.0,
        "{beside} K beside a void"
    );
    assert!(hollow.lost_energy().to_si().is_finite());
    hollow.fill(fine(), |i, _, _| i == 3);
    let why = "a void cell at -50 K filled after a step";
    let refusal = refused_block(why, &mut hollow, &mut waiting(1e-3), dt);
    names(
        why,
        &refusal,
        "block",
        "cell (3, 0, 0)'s temperature",
        -50.0,
    );

    // A void cell's placeholder has only to be a number, which `NaN` and an infinity are not: a
    // face that carries nothing still multiplies it. Set to `NaN`, and deposited into, whose
    // capacity of nothing makes a millijoule `+∞`. (A void cell reads as `NaN` whatever it holds, so
    // what it holds is said here rather than read.)
    for (why, set_it, held) in [
        ("a void cell set to NaN", true, NAN),
        ("a void cell deposited into", false, INF),
    ] {
        let mut b = stepped(block().empty(|i, _, _| i == 3));
        if set_it {
            b.set_temperature(3, 0, 0, Temperature::from_si(NAN));
        } else {
            b.deposit(3, 0, 0, Energy::from_si(1e-3));
        }
        let refusal = refused_block(why, &mut b, &mut waiting(1e-3), dt);
        names(why, &refusal, "block", "cell (3, 0, 0) is void", held);
    }

    // Put right, it runs; and a checkpoint taken before it was put right brings it back, refused.
    let mut fixed = stepped(block());
    fixed.set_temperature(2, 0, 0, Temperature::from_si(-50.0));
    fixed.checkpoint();
    fixed.set_temperature(2, 0, 0, Temperature::celsius(20.0));
    fixed
        .step(Time::ZERO, dt, &mut Exchange::new())
        .expect("put right, it runs");
    fixed.restore();
    let refusal = refused_block("restored to -50 K", &mut fixed, &mut Exchange::new(), dt);
    names(
        "restored to -50 K",
        &refusal,
        "block",
        "cell (2, 0, 0)'s temperature",
        -50.0,
    );
}

#[test]
fn heat_taken_out_below_absolute_zero_is_refused_with_nothing_taken() {
    let start = Temperature::celsius(20.0);
    let t = start.to_si();
    let c = fine().heat_capacity(Volume::from_si(1e-6)).unwrap().to_si();

    // A lumped body under a film, which has lost heat already: one step of it is
    // `T + (Q − L)/C`, `L` what the film takes over it.
    let dt = Time::from_si(1.0);
    let mut m = lump(start, still());
    m.step(Time::ZERO, dt, &mut Exchange::new()).unwrap();
    let now = m.temperature().to_si();
    let loss = still().loss_from(m.temperature(), 0.5).to_si() * dt.to_si();
    let (past, short) = (
        (loss - c * now) * (1.0 + 1e-9),
        (loss - c * now) * (1.0 - 1e-9),
    );
    let why = "a lumped body given a part in a billion past its heat";
    let (lost, mut bus) = (m.lost_energy().to_si(), waiting(past));
    assert!(lost > 0.0, "the body has lost nothing to put back");
    let refusal = m
        .step(Time::ZERO, dt, &mut bus)
        .expect_err(&format!("{why} must be refused"));
    names(why, &refusal, "lump", BELOW, now + (past - loss) / c);
    assert_eq!(m.temperature().to_si(), now, "{why}: moved");
    assert_eq!(bus_of(&bus), bus_of(&waiting(past)), "{why}: took heat");
    assert_eq!(m.absorbed_energy().to_si(), 0.0, "{why}: counted heat");
    assert_eq!(m.lost_energy().to_si(), lost, "{why}: counted heat lost");
    m.step(Time::ZERO, dt, &mut waiting(short))
        .expect("a part in a billion short of its heat runs");
    assert_eq!(m.temperature().to_si(), now + (short - loss) / c);
    assert!(m.temperature().to_si() >= 0.0);

    // Nothing taken out and a step longer than the body's time constant, past which its film
    // overshoots its room: from 1000 K, one and a half of them came to −90.3 K.
    let mut hot = LumpedMass::new(
        "lump",
        stuff(50.0, 0.0),
        Volume::from_si(1e-6),
        Length::mm(1.0),
        Temperature::from_si(1000.0),
        still(),
    );
    let long = Time::from_si(1.5 * hot.time_constant().unwrap().to_si());
    let loss = still().loss_from(hot.temperature(), 0.0).to_si() * long.to_si();
    let would = 1000.0 + (0.0 - loss) / c;
    let why = "a lumped body stepped one and a half time constants from 1000 K";
    let refusal = hot
        .step(Time::ZERO, long, &mut Exchange::new())
        .expect_err(&format!("{why} must be refused"));
    names(
        why,
        &refusal,
        "lump",
        "longer than the body's time constant",
        would,
    );
    assert!(
        (would + 90.275).abs() < 1e-9,
        "{why}: the closed form is {would} K"
    );
    assert_eq!(hot.temperature().to_si(), 1000.0, "{why}: moved");
    assert_eq!(hot.lost_energy().to_si(), 0.0, "{why}: counted heat lost");

    // A share of an interval: the first half runs, the second is refused, and its share stays.
    let mut m = lump(start, sealed());
    let mut bus = Exchange::new();
    bus.covering(Time::from_si(2.0));
    bus.publish(HEAT, -1.5 * c * t);
    m.step(Time::ZERO, Time::from_si(1.0), &mut bus)
        .expect("the first half takes out three quarters of its heat");
    let left = bus_of(&bus);
    let refusal = m
        .step(Time::from_si(1.0), Time::from_si(1.0), &mut bus)
        .expect_err("the second half takes out more than is left");
    assert!(refusal.quantity.contains(BELOW), "{refusal}");
    assert_eq!(bus_of(&bus), left, "a refused substep took its share");

    // A network node under a film, the second of two, which has lost heat already.
    let mut net = ThermalNetwork::new("net");
    net.node(
        "idle",
        fine(),
        Volume::from_si(1e-6),
        Length::mm(1.0),
        start,
    );
    let n = net.node_losing_to(
        "n",
        fine(),
        Volume::from_si(1e-6),
        Length::mm(1.0),
        start,
        still(),
    );
    net.absorbing(n).unwrap();
    net.step(Time::ZERO, dt, &mut Exchange::new()).unwrap();
    let now = net.temperature(n).to_si();
    let loss = still().loss_from(net.temperature(n), 0.5).to_si() * dt.to_si();
    let (past, short) = (
        (loss - c * now) * (1.0 + 1e-9),
        (loss - c * now) * (1.0 - 1e-9),
    );
    let why = "a network node given a part in a billion past its heat";
    let (lost, mut bus) = (net.lost_energy().to_si(), waiting(past));
    let refusal = net
        .step(Time::ZERO, dt, &mut bus)
        .expect_err(&format!("{why} must be refused"));
    names(why, &refusal, "net/n", BELOW, now + (past - loss) / c);
    assert_eq!(net.temperature(n).to_si(), now, "{why}: moved");
    assert_eq!(bus_of(&bus), bus_of(&waiting(past)), "{why}: took heat");
    assert_eq!(net.absorbed_energy().to_si(), 0.0, "{why}: counted heat");
    assert_eq!(net.lost_energy().to_si(), lost, "{why}: counted heat lost");
    net.step(Time::ZERO, dt, &mut waiting(short))
        .expect("a part in a billion short runs");
    assert_eq!(net.temperature(n).to_si(), now + (short - loss) / c);

    // A bar: uniform, so its sweep moves nothing and a cell's step is `T + Q/C`, on its first
    // cell from the channel and on a face other than the first from its boundary.
    let cell = fine()
        .heat_capacity(Volume::from_si(1e-6 * 1e-3))
        .unwrap()
        .to_si();
    let (past, short) = (-cell * t * (1.0 + 1e-9), -cell * t * (1.0 - 1e-9));
    let why = "a bar's first cell given a part in a billion past its heat";
    let mut b = bar(start);
    let dt = a_step(&b);
    let mut bus = waiting(past);
    let refusal = b
        .step(Time::ZERO, dt, &mut bus)
        .expect_err(&format!("{why} must be refused"));
    names(
        why,
        &refusal,
        "bar",
        &format!("cell 0 {BELOW}"),
        t + past / cell,
    );
    assert_eq!(cells_of_bar(&b), cells_of_bar(&bar(start)), "{why}: moved");
    assert_eq!(bus_of(&bus), bus_of(&waiting(past)), "{why}: took heat");
    assert_eq!(b.absorbed_energy().to_si(), 0.0, "{why}: counted heat");
    b.step(Time::ZERO, dt, &mut waiting(short))
        .expect("a part in a billion short runs");
    assert_eq!(b.temperature_at(0).to_si(), t + short / cell);

    let why = "a bar's fourth face given a part in a billion past its heat";
    let mut b = bar(start);
    let side = b.boundary().unwrap().clone();
    let mut faces = vec![0.0; 8];
    faces[1] = 1e-6;
    faces[3] = past;
    let mut bus = waiting(1e-6);
    bus.publish_on(&side, HEAT, &Flux::from_faces(faces.clone()))
        .unwrap();
    let refusal = b
        .step(Time::ZERO, dt, &mut bus)
        .expect_err(&format!("{why} must be refused"));
    names(
        why,
        &refusal,
        "bar",
        &format!("cell 3 {BELOW}"),
        t + past / cell,
    );
    assert_eq!(cells_of_bar(&b), cells_of_bar(&bar(start)), "{why}: moved");
    assert_eq!(bus.peek(HEAT), 1e-6, "{why}: took heat off the channel");
    let offered = bus.peek_on(&side, HEAT).expect("still on the faces");
    assert_eq!(
        offered.per_face().to_vec(),
        faces,
        "{why}: took heat off the faces"
    );
    assert_eq!(b.absorbed_energy().to_si(), 0.0, "{why}: counted heat");
    faces[3] = short;
    let mut bus = Exchange::new();
    bus.publish_on(&side, HEAT, &Flux::from_faces(faces))
        .unwrap();
    b.step(Time::ZERO, dt, &mut bus)
        .expect("a part in a billion short runs");
    assert_eq!(b.temperature_at(3).to_si(), t + short / cell);

    // A bar with a gradient, so a refused step that swept would show: the whole bar as it was.
    let mut b = bar(start);
    b.step(Time::ZERO, dt, &mut waiting(1e-3)).unwrap();
    let before = cells_of_bar(&b);
    let why = "a bar with a gradient given a kilojoule out";
    let refusal = b
        .step(Time::ZERO, dt, &mut waiting(-1e3))
        .expect_err(&format!("{why} must be refused"));
    assert!(
        refusal.quantity.contains(&format!("cell 0 {BELOW}")),
        "{refusal}"
    );
    assert_eq!(cells_of_bar(&b), before, "{why}: moved");

    // A block of one cell with nothing to lose to, by its source and by the bus.
    let one = || Solid3D::new("block", fine(), (1, 1, 1), Length::mm(1.0), start);
    let c1 = one().heat_capacity().to_si();
    let dt = Time::from_si(1e-3);
    let watts = |share: f64| -c1 * t * share / dt.to_si();
    let why = "a block's sink a part in a billion past its heat";
    let mut b = one().dissipating(watts(1.0 + 1e-9), |_, _, _| true);
    let refusal = refused_block(why, &mut b, &mut Exchange::new(), dt);
    let would = t + watts(1.0 + 1e-9) * dt.to_si() / c1;
    // Said as the heat taken out, which it is, and not as a cell that froze.
    names(
        why,
        &refusal,
        "block",
        &format!("the heat taken out in this step would take cell (0, 0, 0) {BELOW}"),
        would,
    );
    let mut b = one().dissipating(watts(1.0 - 1e-9), |_, _, _| true);
    b.step(Time::ZERO, dt, &mut Exchange::new())
        .expect("a part in a billion short runs");
    assert_eq!(
        b.temperature_at(0, 0, 0).to_si(),
        t + watts(1.0 - 1e-9) * dt.to_si() / c1
    );

    let (past, short) = (-c1 * t * (1.0 + 1e-9), -c1 * t * (1.0 - 1e-9));
    let why = "a block given a part in a billion past its heat off the bus";
    let refusal = refused_block(why, &mut one(), &mut waiting(past), dt);
    names(
        why,
        &refusal,
        "block",
        &format!("cell (0, 0, 0) {BELOW}"),
        t + past / c1,
    );
    let mut b = one();
    b.step(Time::ZERO, dt, &mut waiting(short))
        .expect("a part in a billion short runs");
    assert_eq!(b.temperature_at(0, 0, 0).to_si(), t + short / c1);

    // In a block of four, by the cell the sink is in, and from books that are not zero.
    let why = "a sink in the third of four cells";
    let dt = a_step(&block());
    let mut b = stepped(block()).dissipating(-1e3, |i, _, _| i == 2);
    let refusal = refused_block(why, &mut b, &mut waiting(1e-3), dt);
    assert!(
        refusal
            .quantity
            .contains(&format!("cell (2, 0, 0) {BELOW}")),
        "{refusal}"
    );

    // The bus takes from every cell: the one it takes below zero is named, whichever it is, and a
    // void cell's placeholder below zero is not asked. A step too short to conduct anything.
    let tiny = Time::from_si(1e-12);
    let mut b = stepped(block());
    b.set_temperature(2, 0, 0, Temperature::from_si(0.5));
    b.step(Time::ZERO, tiny, &mut Exchange::new()).unwrap();
    let rise = -b.heat_capacity().to_si();
    let why = "a kelvin off every cell, one of them at half a kelvin";
    let refusal = refused_block(why, &mut b, &mut waiting(rise), tiny);
    assert!(
        refusal
            .quantity
            .contains(&format!("cell (2, 0, 0) {BELOW}")),
        "{refusal}"
    );
    let mut b = block().empty(|i, _, _| i == 3);
    b.set_temperature(3, 0, 0, Temperature::from_si(-50.0));
    let mut b = stepped(b);
    let rise = -0.1 * b.heat_capacity().to_si();
    b.step(Time::ZERO, tiny, &mut waiting(rise))
        .expect("a void cell's placeholder below zero is not asked");
}

/// **A refused step puts the phases back too.** Ice with a cell half melted beside one a sink
/// takes below zero: the sweep froze some of the mushy cell, by conduction to its colder
/// neighbour, before the sink's cell was found below zero.
#[test]
fn a_refused_step_puts_a_melting_cell_back() {
    let mut b = Solid3D::new(
        "ice",
        Substance::ice(),
        (3, 1, 1),
        Length::mm(1.0),
        Temperature::celsius(-1.0),
    );
    b.set_melted_fraction(1, 0, 0, 0.5);
    let mut b = b.dissipating(-1e3, |i, _, _| i == 0);
    let dt = a_step(&b);
    let refusal = refused_block("ice with a sink", &mut b, &mut Exchange::new(), dt);
    assert!(
        refusal
            .quantity
            .contains(&format!("cell (0, 0, 0) {BELOW}")),
        "{refusal}"
    );
    assert_eq!(
        b.melted_fraction_at(1, 0, 0),
        0.5,
        "the mushy cell was not put back"
    );

    // And the same step without the sink freezes the mushy cell, so the fraction put back is one
    // that moved.
    let mut free = Solid3D::new(
        "ice",
        Substance::ice(),
        (3, 1, 1),
        Length::mm(1.0),
        Temperature::celsius(-1.0),
    );
    free.set_melted_fraction(1, 0, 0, 0.5);
    free.step(Time::ZERO, dt, &mut Exchange::new()).unwrap();
    assert!(
        free.melted_fraction_at(1, 0, 0) < 0.5,
        "the sweep froze nothing"
    );

    // And a block with water as its liquid, whose operator moves with its phases: the only
    // refusal a two-phase block has is heat taken out, since its limit bounds a freeze, and it
    // puts back what its sweep moved as a one-phase block's does.
    let two = || {
        let mut b = Solid3D::new(
            "two-phase",
            two_phase_ice(4182.0),
            (3, 1, 1),
            Length::mm(1.0),
            Temperature::celsius(-1.0),
        );
        b.set_melted_fraction(1, 0, 0, 0.5);
        b
    };
    let mut b = two().dissipating(-1e3, |i, _, _| i == 0);
    let dt = a_step(&b);
    let refusal = refused_block(
        "two-phase ice with a sink",
        &mut b,
        &mut Exchange::new(),
        dt,
    );
    assert!(
        refusal
            .quantity
            .contains(&format!("cell (0, 0, 0) {BELOW}")),
        "{refusal}"
    );
    assert_eq!(
        b.melted_fraction_at(1, 0, 0),
        0.5,
        "the two-phase block's mushy cell was not put back"
    );
    let mut free = two();
    free.step(Time::ZERO, dt, &mut Exchange::new()).unwrap();
    assert!(
        free.melted_fraction_at(1, 0, 0) < 0.5,
        "the two-phase sweep froze nothing"
    );
}

/// **A liquid freezing at its block's limit stays above absolute zero**, with nothing taken out. A
/// step moves a cell's heat in enthalpy at its liquid's capacity and gives it back as temperature
/// at its solid's, and water holds twice what ice does: water at 1 °C between two cells of ice at
/// 1 K, stepped at its block's limit, froze to −108.4 K, and one cell of it under films to a room
/// at absolute zero to −120.1 K. A pass over every cell of a two-phase block refused those. The
/// limit holds such a cell to what it meets now, and bounds it — see
/// `a_phase_change_is_bounded_by_what_the_cell_meets.rs` — so each stepped at its limit stays at or
/// above zero, and the step the limit was, the ice's `ρdx²c_S/k_f` in closed form, is refused as
/// 1.2506 times it. A liquid that holds what its solid does stays above zero, and so does water at
/// half the limit.
#[test]
fn a_liquid_frozen_at_its_limit_stays_above_absolute_zero() {
    let why = "water at 1 C between ice at 1 K";
    let row = || {
        let mut row = Solid3D::new(
            "ice",
            two_phase_ice(4182.0),
            (3, 1, 1),
            Length::mm(1.0),
            Temperature::from_si(1.0),
        );
        row.set_temperature(1, 0, 0, Temperature::from_si(274.0));
        row
    };
    let mut stepped = row();
    let dt = stepped.max_stable_dt(Time::ZERO);
    stepped
        .step(Time::ZERO, dt, &mut Exchange::new())
        .expect("water between ice at its limit");
    for i in 0..3 {
        let t = stepped.temperature_at(i, 0, 0).to_si();
        assert!(t >= 0.0, "{why}: cell {i} at {t} K");
    }
    // The ice's row, one face to the water: `k_f = 2·k_ice·k_water/(k_ice + k_water)`.
    let ice = Substance::ice();
    let thermal = ice.thermal.unwrap();
    let k_ice = thermal.conductivity.to_si();
    let k_f = 2.0 * k_ice * 0.6 / (k_ice + 0.6);
    let was = ice.density.to_si() * 1e-6 * thermal.specific_heat.to_si() / k_f;
    let refusal = refused_block(why, &mut row(), &mut Exchange::new(), Time::from_si(was));
    let past = refusal.after / refusal.before;
    println!("  {why}, at the ice's limit of {was} s -> {refusal}");
    assert!(
        refusal.quantity == "Fourier number" && (past - 1.2506).abs() < 1e-4,
        "{why}: at the ice's limit, {refusal}"
    );

    let one = |c_l: f64, room: f64| {
        let mut b = Solid3D::new(
            "water",
            two_phase_ice(c_l),
            (1, 1, 1),
            Length::mm(1.0),
            Temperature::from_si(274.15),
        );
        for face in [
            Face::XMin,
            Face::XMax,
            Face::YMin,
            Face::YMax,
            Face::ZMin,
            Face::ZMax,
        ] {
            b = b.losing_from(face, film(1e3, 1e-6, Temperature::from_si(room)));
        }
        b
    };
    let why = "a cell of water under films to a room at 0 K";
    let mut b = one(4182.0, 0.0);
    let dt = b.max_stable_dt(Time::ZERO);
    b.step(Time::ZERO, dt, &mut Exchange::new()).expect(why);
    // On the room, which the limit lands the cell on, to the rounding of an enthalpy near the
    // melting point: measured, one unit in the last place of 273.15 K below it.
    let t = b.temperature_at(0, 0, 0).to_si();
    assert!(t >= -16.0 * f64::EPSILON * 273.15, "{why}: at {t} K");
    let ice = Substance::ice().thermal.unwrap().specific_heat.to_si();
    let mut b = one(ice, 0.0);
    let dt = b.max_stable_dt(Time::ZERO);
    b.step(Time::ZERO, dt, &mut Exchange::new())
        .expect("a liquid that holds what its solid does stays above zero");
    assert!(b.temperature_at(0, 0, 0).to_si() > 0.0);
    let mut b = one(4182.0, 0.0);
    let dt = Time::from_si(0.5 * b.max_stable_dt(Time::ZERO).to_si());
    b.step(Time::ZERO, dt, &mut Exchange::new())
        .expect("water at half the limit stays above zero");
    assert!(b.temperature_at(0, 0, 0).to_si() > 0.0);
}

/// **A sink that outruns its film is refused at the step that would cross zero**, and every step
/// before it ran: a block whose 50 mW sink outran the 2 mW its film could bring in stepped twenty
/// thousand times to −354.67 K. The refused step leaves the block, its books — which by then are
/// not zero — and the bus as they were. **Its balance is refused at its one root below absolute
/// zero**, where it was refused as not converging, and then as having run out of sweeps while the
/// steady solve read its film's slope at the room.
#[test]
fn a_sink_that_outruns_its_film_is_refused_where_it_crosses_zero() {
    let mut b = Solid3D::new(
        "block",
        fine(),
        (4, 1, 1),
        Length::mm(1.0),
        Temperature::celsius(20.0),
    )
    .losing_from(Face::XMax, film(7.0, 1e-6, Temperature::celsius(0.0)))
    .dissipating(-0.05, |_, _, _| true);
    let refusal = b.steady_state().unwrap_err();
    // The root: the film's surface `s` where `hA(T∞ − s) + εσA(T∞⁴ − s·|s|³)` brings the 50 mW
    // in, found by halving, and the first cell 2 K below it — half a cell of `k = 50` carrying all
    // of it, 0.5 K, then three faces carrying three quarters, a half and a quarter, 0.75, 0.5 and
    // 0.25 K.
    let (room, area) = (Temperature::celsius(0.0).to_si(), 1e-6);
    let brings = |s: f64| {
        7.0 * area * (room - s)
            + 0.5
                * STEFAN_BOLTZMANN.to_si()
                * area
                * ((room * room) * (room * room) - s * s * s * s.abs())
    };
    let (mut lo, mut hi) = (-1e4, room);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if brings(mid) > 0.05 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let root = 0.5 * (lo + hi) - 2.0;
    assert!(
        refusal.quantity.starts_with(
            "no balance above absolute zero was found: cell (0, 0, 0) came to below it"
        ) && refusal
            .quantity
            .ends_with("which take more than can reach them")
            && (refusal.before - root).abs() <= 1e-5,
        "refused at {} K where its root is {root} K: {refusal}",
        refusal.before
    );
    let dt = a_step(&b);
    let mut bus = Exchange::new();
    let mut steps = 0;
    let refusal = loop {
        let (cells, books, heat) = (cells_of_block(&b), books_of_block(&b), bus_of(&bus));
        match b.step(Time::ZERO, dt, &mut bus) {
            Ok(()) => steps += 1,
            Err(v) => {
                assert_eq!(
                    cells_of_block(&b),
                    cells,
                    "the refused step moved the block"
                );
                assert_eq!(
                    books_of_block(&b),
                    books,
                    "the refused step moved the books"
                );
                assert_eq!(bus_of(&bus), heat, "the refused step took heat");
                break v;
            }
        }
        assert!(steps < 100_000, "never refused");
    };
    println!("  refused after {steps} steps: {refusal}");
    assert!(refusal.quantity.contains(BELOW), "{refusal}");
    assert!(refusal.before < 0.0);
    assert!(b.lost_energy().to_si() < 0.0 && b.generated_energy().to_si() < 0.0);
    for i in 0..4 {
        let t = b.temperature_at(i, 0, 0).to_si();
        assert!(t.is_finite() && t >= 0.0, "a step before it went to {t} K");
    }
}

/// Every cell of `b` is a number at or above zero.
fn every_cell_above_zero(why: &str, b: &Solid3D) {
    let (nx, ny, nz) = b.counts();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                if b.is_void(i, j, k) {
                    continue;
                }
                let t = b.temperature_at(i, j, k).to_si();
                assert!(
                    t.is_finite() && t >= 0.0,
                    "{why}: cell ({i}, {j}, {k}) at {t} K"
                );
            }
        }
    }
}

/// **Only what can take a cell below zero is asked, and nothing else goes below it.** A block with
/// no sink and no phase change, in a room at absolute zero, cooling through films on every face
/// and across a clearance, stepped at its own limit: near zero, where conduction and convection
/// carry it, for twenty thousand steps until it has all but reached its room; and from 1000 K,
/// black and barely conducting, under films that only radiate, where radiation carries every row.
/// Every cell is a number at or above zero after every step, which is the reason a step that takes
/// nothing out of a one-phase block is not asked.
#[test]
fn without_a_sink_nothing_goes_below_absolute_zero() {
    let zero = Temperature::from_si(0.0);
    let faces = [
        Face::XMin,
        Face::XMax,
        Face::YMin,
        Face::YMax,
        Face::ZMin,
        Face::ZMax,
    ];
    for (why, substance, start, h, steps) in [
        ("near zero", stuff(50.0, 1.0), 1.0, 1e3, 20_000),
        (
            "radiating from 1000 K",
            stuff(0.05, 1.0),
            1000.0,
            0.0,
            4_000,
        ),
    ] {
        let mut b = Solid3D::new(
            "block",
            substance,
            (5, 2, 2),
            Length::mm(1.0),
            Temperature::from_si(start),
        );
        for face in faces {
            b = b.losing_from(face, film(h, 1e-6, zero));
        }
        let mut b = b.empty(|i, _, _| i == 2);
        for _ in 0..steps {
            let dt = b.max_stable_dt(Time::ZERO);
            b.step(Time::ZERO, dt, &mut Exchange::new())
                .expect("a block with no sink is not refused");
            every_cell_above_zero(why, &b);
        }
        let warmest = b.peak_temperature().to_si();
        println!("  {why}: warmest at the end {warmest:e} K");
        assert!(
            warmest.is_finite() && warmest < start,
            "{why}: it did not cool"
        );
    }
}

/// **A network's balance below absolute zero is refused**, against its closed form. A node alone
/// under a still-air film, its radiation held at absolute zero below it, balances at
/// `Tₐ + (P + εσA·Tₐ⁴)/hA`. Below zero it is refused with that, to the solve's own error — its stop
/// leaves `1e-12·|P|` W of the balance, which over `hA` is 1.4e-9 K here, so 1e-11 of the answer —
/// for a surface that convects and one that radiates too, and with its path conductance; a part in
/// a billion past the bound is refused, and a part in a billion short of it is answered at its
/// closed form, `1e-9·Tₐ` on convection alone. In a room at 77 K a part in a million past is
/// refused: a bound of 1e-4 K for every node answered it at zero. A node that only radiates is
/// refused as having no balance above zero, where it was refused as singular; and in a chain the
/// node named is the coldest.
#[test]
fn a_network_balance_below_absolute_zero_is_refused() {
    let ta = Temperature::celsius(20.0).to_si();
    let (h, area) = (7.0, 0.01);
    let alone = |emissivity: f64, h: f64, room: f64| {
        let mut net = ThermalNetwork::new("net");
        let n = net.node_losing_to(
            "n",
            stuff(50.0, emissivity),
            Volume::from_si(1e-4),
            Length::mm(5.0),
            Temperature::from_si(room),
            film(h, area, Temperature::from_si(room)),
        );
        net.absorbing(n).unwrap();
        (net, n)
    };
    for emissivity in [0.0, 0.5] {
        let (net, n) = alone(emissivity, h, ta);
        let er = emissivity * STEFAN_BOLTZMANN.to_si() * area;
        let s = ta * ta;
        let most = h * area * ta + er * s * s;
        let why = format!("100 W out of a node, emissivity {emissivity}");
        let refusal = net.steady_state(Power::w(-100.0)).unwrap_err();
        let want = ta + (-100.0 + er * s * s) / (h * area);
        assert!(
            refusal.site == "net/n"
                && refusal.quantity.contains("no balance above absolute zero")
                && (refusal.before - want).abs() <= 1e-11 * want.abs(),
            "{why}: refused with {refusal}, where the balance is at {want} K"
        );
        let path = net.path_conductance(n, Power::w(-100.0)).unwrap_err();
        identical(&why, "path_conductance", &path, &refusal);

        let past = net
            .steady_state(Power::from_si(-most * (1.0 + 1e-9)))
            .unwrap_err();
        assert!(
            past.quantity.contains("no balance above absolute zero") && past.before < 0.0,
            "{why}: a part in a billion past the bound: {past}"
        );
        let short = -most * (1.0 - 1e-9);
        let got = net
            .steady_state(Power::from_si(short))
            .unwrap()
            .temperature(n)
            .to_si();
        let want = ta + (short + er * s * s) / (h * area);
        assert!(got >= 0.0, "{why}: {got} K");
        if emissivity == 0.0 {
            assert!(
                (got - want).abs() <= 10.0 * 1e-12 * most / (h * area),
                "{why}: a part in a billion short came to {got} K, where it is {want} K"
            );
        }
    }

    // At the bound to the bit the root is absolute zero, and the solve's comes out a rounding below
    // it — `−hA·Tₐ/hA` is a hair past `−Tₐ` — under a film of 1e-7 over a square millimetre: it is
    // answered as zero, and not as a temperature below it.
    let mut net = ThermalNetwork::new("net");
    let n = net.node_losing_to(
        "n",
        stuff(50.0, 0.0),
        Volume::from_si(1e-6),
        Length::mm(1.0),
        Temperature::from_si(ta),
        film(1e-7, 1e-6, Temperature::from_si(ta)),
    );
    net.absorbing(n).unwrap();
    let at = net
        .steady_state(Power::from_si(-(1e-7 * 1e-6 * ta)))
        .unwrap()
        .temperature(n)
        .to_si();
    assert_eq!(at, 0.0, "at the bound to the bit: {at:e} K");

    let (net, _) = alone(0.0, h, 77.0);
    let most = h * area * 77.0;
    let past = net
        .steady_state(Power::from_si(-most * (1.0 + 1e-6)))
        .unwrap_err();
    assert!(
        past.quantity.contains("no balance above absolute zero"),
        "a part in a million past the bound in a room at 77 K: {past}"
    );

    let (net, _) = alone(0.5, 0.0, ta);
    let why = "100 W out of a node that only radiates";
    let refusal = net.steady_state(Power::w(-100.0)).unwrap_err();
    println!("  {why} -> {refusal}");
    assert!(
        refusal.site == "net/n"
            && refusal
                .quantity
                .contains("no balance above absolute zero was found")
            && refusal.before < 0.0,
        "{why}: {refusal}"
    );

    // A chain: the heat taken out of an interior node behind a joint from the one that loses, so
    // much that both come to below zero, the first −135.4 K and the second −435.4 K: the one named
    // is the coldest, not the first.
    let mut net = ThermalNetwork::new("net");
    let room = net.node_losing_to(
        "room-side",
        stuff(50.0, 0.0),
        Volume::from_si(1e-4),
        Length::mm(5.0),
        Temperature::from_si(ta),
        film(h, area, Temperature::from_si(ta)),
    );
    let deep = net.node(
        "deep",
        stuff(50.0, 0.0),
        Volume::from_si(1e-4),
        Length::mm(5.0),
        Temperature::from_si(ta),
    );
    net.link(room, deep, Conductance::w_per_k(0.1)).unwrap();
    net.absorbing(deep).unwrap();
    let refusal = net.steady_state(Power::w(-30.0)).unwrap_err();
    assert!(
        refusal.site == "net/deep" && refusal.quantity.contains("no balance above absolute zero"),
        "the coldest node is the deep one: {refusal}"
    );
}
