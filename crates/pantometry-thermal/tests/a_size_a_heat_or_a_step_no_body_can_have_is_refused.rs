//! A size, a source, a heat or a step that no body can have is refused by name before anything
//! moves, and every question asked of a body refuses what its step refuses, with the same refusal.
//!
//! **Each of these answered, and none said anything.** Measured before this change, on the four
//! thermal domains stepped by hand:
//!
//! - **A size.** A bar whose cell length or cross-section was negative was given a millijoule and
//!   its heated cell cooled, 293.1459 K against 293.1541 K, with its books saying the millijoule
//!   had arrived. A cross-section of zero put the cell at `+∞`, one of `NaN` at `NaN`, and one of
//!   `+∞` moved nothing. A block of cells a millimetre less than nothing cooled under a milliwatt,
//!   to 292.90 K where it warms to 295.34 K. An exposed face's area of zero, −1 or `NaN` was read
//!   as zero and stepped. A lumped body's or a node's thickness was never asked at all, and a
//!   node's volume of zero, `NaN` or less was solved for by `ThermalNetwork::steady_state` as if
//!   it were fine. A cell size or a cell length of zero or `NaN` was refused as a substance with no
//!   diffusivity, which sent its caller to a substance that was fine.
//! - **A source.** `+∞` or `−∞` watts in a block's cell stepped it there and answered `Ok`. `NaN`
//!   made the cell `NaN`, which the block's warmest and coldest readings step over.
//! - **A heat.** `NaN` or `+∞` joules on the channel, or on one of a bar's faces, stepped every
//!   domain to `NaN` or `+∞`. A network with no node to absorb heat refused it after taking it off
//!   the bus and counting it as absorbed.
//! - **A step.** `NaN` seconds stepped a lumped body and a network to `NaN`. −1 s ran a block's
//!   conduction backwards, from 350 K beside 293 K to 4274.5 K beside −3631.2 K, and a network
//!   stepped −1 s did nothing and called it success.
//! - **A question.** A lumped body's `equilibrium_rise` answered 71.43 K for an emissivity of `NaN`
//!   where the body settles 37.81 K up, −293.15 K for a film area of `NaN`, `+∞` for a power of
//!   `NaN` or `−∞`, and −14 231.87 K for a kilowatt taken out of a body its room could bring 24 W
//!   to. Its Biot number was −0.058 for a conductivity of −0.6, which passes the `< 0.1` a caller
//!   asks of it; its time constant −3294 s for a volume below zero. A network's time constant was
//!   0 s, `+∞` and −800 s for a node's volume of zero, `NaN` and less.
//!
//! # What is checked against what
//!
//! **A refusal names the number and carries it**, bit for bit, at the domain or the node it is in,
//! and a source by its cell, a face's heat by its face. **Nothing moved**: every temperature is the
//! one it started at, bit for bit, the heat waiting on the bus is still there, on the channel and on
//! the faces, and nothing is counted as taken, lost or generated. **Wherever the number is**: on
//! the second node, in a cell other than the first, in a source written after another, on a face
//! other than the first.
//!
//! **Every question refuses what the step refuses, with the same refusal**: the same words, the same
//! site and the same value. The step is the reference because it is the one that was already right
//! for a substance and a film, and a question that answers for a body the step refuses is an answer
//! about a body that cannot be one.
//!
//! **The limits a caller can mean still run**, against their closed forms: a step of no length, the
//! balance down to the edge of absolute zero against the root of `hA·ΔT + εσA((Tₐ+ΔT)⁴ − Tₐ⁴) = P`
//! and against marching there, each loss path alone, and a Biot number where nothing convects.

use pantometry_core::substance::ThermalProps;
use pantometry_core::units::{
    Area, Conductance, Density, Length, Power, SpecificHeat, Temperature, ThermalConductivity,
    ThermalExpansion, Time, Volume, STEFAN_BOLTZMANN,
};
use pantometry_core::{Domain, Exchange, Flux, Schedule, Simulation, Substance, Violation};
use pantometry_thermal::{
    Axis, Bar1D, Environment, Face, LumpedMass, Node, Solid3D, ThermalNetwork, HEAT,
};

const NAN: f64 = f64::NAN;
const INF: f64 = f64::INFINITY;
/// Joules waiting on the bus when a step is refused, which it must leave there.
const WAITING: f64 = 1e-3;

fn stuff(rho: f64, k: f64, c: f64, emissivity: f64) -> Substance {
    Substance::bulk("stuff", Density::from_si(rho)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(k),
        specific_heat: SpecificHeat::j_per_kg_k(c),
        expansion: ThermalExpansion::ppm_per_k(50.0),
        emissivity,
    })
}

/// A substance every domain can step: water-like, and radiating, so the emissivity is read.
fn fine() -> Substance {
    stuff(1000.0, 50.0, 4000.0, 0.5)
}

/// A substance that says nothing of its heat.
fn silent() -> Substance {
    Substance::bulk("silent", Density::from_si(1000.0))
}

fn film(h: f64, area: f64, ambient: Temperature) -> Environment {
    Environment {
        ambient,
        convection_w_per_m2_k: h,
        area: Area::from_si(area),
    }
}

fn room() -> Temperature {
    Temperature::celsius(0.0)
}

/// A film colder than every body here, so it acts from the first step.
fn still() -> Environment {
    film(10.0, 1e-6, room())
}

fn waiting(joules: f64) -> Exchange {
    let mut bus = Exchange::new();
    bus.publish(HEAT, joules);
    bus
}

/// The same value: equal bits, or both not a number.
fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan())
}

/// `refusal` names `names`, carries `value`, and was raised at `site`.
fn names(why: &str, refusal: &Violation, site: &str, names: &str, value: f64) {
    println!("  {why:44} -> {refusal}");
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

/// Two refusals that are the same refusal: the same words, the same site, the same value.
fn identical(why: &str, asked: &str, got: &Violation, step: &Violation) {
    assert!(
        got.quantity == step.quantity && got.site == step.site && same(got.before, step.before),
        "{why}: {asked} refused with {got}, where the step refused with {step}"
    );
}

/// What a refused step must have left on the bus: the heat that was waiting, all of it.
fn left_waiting(why: &str, bus: &Exchange, joules: f64) {
    assert!(
        same(bus.peek(HEAT), joules),
        "{why}: a refused step took heat off the bus, {} of {joules} left",
        bus.peek(HEAT)
    );
}

fn lump(substance: Substance, environment: Environment, volume: f64, thickness: f64) -> LumpedMass {
    LumpedMass::new(
        "lump",
        substance,
        Volume::from_si(volume),
        Length::from_si(thickness),
        Temperature::celsius(100.0),
        environment,
    )
}

fn a_lump() -> LumpedMass {
    lump(fine(), still(), 1e-6, 1e-3)
}

/// Eight cells of `dx` at 20 °C.
fn bar(dx: f64, area: f64) -> Bar1D {
    Bar1D::new(
        "bar",
        fine(),
        8,
        Length::from_si(dx),
        Area::from_si(area),
        Temperature::celsius(20.0),
    )
}

/// `b` given a millijoule in its first cell and stepped once, so it has a gradient: on a uniform
/// bar a sweep moves nothing, and a refused step that had swept would look like one that had not.
fn warmed(mut b: Bar1D) -> Bar1D {
    let dt = a_step(&b);
    b.step(Time::ZERO, dt, &mut waiting(1e-3))
        .expect("a fine bar steps");
    assert!(
        b.temperature_at(0).to_si() > b.temperature_at(1).to_si(),
        "the bar has no gradient to move"
    );
    b
}

/// Four cells, the first at 100 °C and the rest at 20 °C, under a film on the far face.
fn block(dx: f64) -> Solid3D {
    let mut b = Solid3D::new(
        "block",
        fine(),
        (4, 1, 1),
        Length::from_si(dx),
        Temperature::celsius(20.0),
    )
    .losing_from(Face::XMax, still());
    b.set_temperature(0, 0, 0, Temperature::celsius(100.0));
    b
}

/// Two nodes joined by a link, the first interior and absorbing at 100 °C, the second losing to
/// `environment` at 20 °C and made of `substance` with this volume and thickness: the mistake on
/// the second node, where a check that read only the first would miss it.
fn network(
    substance: Substance,
    environment: Environment,
    volume: f64,
    thickness: f64,
) -> ThermalNetwork {
    network_of(1e-3, substance, environment, volume, thickness)
}

/// The same, with the interior node's thickness given: a thickness only the Biot number reads, of
/// a node that has no Biot number.
fn network_of(
    first_thickness: f64,
    substance: Substance,
    environment: Environment,
    volume: f64,
    thickness: f64,
) -> ThermalNetwork {
    let mut net = ThermalNetwork::new("net");
    let first = net.node(
        "a",
        fine(),
        Volume::from_si(1e-6),
        Length::from_si(first_thickness),
        Temperature::celsius(100.0),
    );
    let second = net.node_losing_to(
        "b",
        substance,
        Volume::from_si(volume),
        Length::from_si(thickness),
        Temperature::celsius(20.0),
        environment,
    );
    net.link(first, second, Conductance::from_si(0.01))
        .expect("a fine link");
    net.absorbing(first).expect("its own node");
    net
}

fn a_network() -> ThermalNetwork {
    network(fine(), still(), 1e-6, 1e-3)
}

fn node(net: &ThermalNetwork, label: &str) -> Node {
    net.node_named(label).expect("a node of this network")
}

/// Every node's temperature, as bits.
fn nodes(net: &ThermalNetwork) -> Vec<u64> {
    net.handles()
        .map(|(n, _)| net.temperature(n).to_si().to_bits())
        .collect()
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
            }
        }
    }
    out
}

/// A step long enough to move anything and short enough to be stable where the limit is a number.
fn a_step(d: &dyn Domain) -> Time {
    let limit = d.max_stable_dt(Time::ZERO).to_si();
    Time::from_si(if limit.is_finite() && limit > 0.0 {
        0.5 * limit
    } else {
        1e-3
    })
}

/// Step a lumped body that must be refused and check that nothing moved. Returns the refusal.
fn refused_lump(why: &str, m: &mut LumpedMass, bus: &mut Exchange, dt: Time) -> Violation {
    let (before, heat) = (m.temperature().to_si().to_bits(), bus.peek(HEAT));
    let refusal = m
        .step(Time::ZERO, dt, bus)
        .expect_err(&format!("{why} must be refused"));
    assert_eq!(
        m.temperature().to_si().to_bits(),
        before,
        "{why}: a refused step moved the body"
    );
    left_waiting(why, bus, heat);
    assert_eq!(
        m.absorbed_energy().to_si(),
        0.0,
        "{why}: counted heat taken"
    );
    assert_eq!(m.lost_energy().to_si(), 0.0, "{why}: counted heat lost");
    refusal
}

/// Every question asked of a lumped body its step refuses refuses as the step did, and its limit is
/// no limit, so a scheduler steps it once and the refusal is the body's.
fn every_lump_question(why: &str, m: &LumpedMass, step: &Violation) {
    identical(why, "biot_number", &m.biot_number().unwrap_err(), step);
    identical(why, "time_constant", &m.time_constant().unwrap_err(), step);
    for watts in [0.0, 1.0, -1.0] {
        let asked = format!("equilibrium_rise({watts} W)");
        let got = m.equilibrium_rise(Power::w(watts)).unwrap_err();
        identical(why, &asked, &got, step);
    }
    assert_eq!(
        m.max_stable_dt(Time::ZERO).to_si(),
        INF,
        "{why}: a body its step refuses has no limit to plan by"
    );
}

fn refused_network(why: &str, net: &mut ThermalNetwork, bus: &mut Exchange, dt: Time) -> Violation {
    let (before, heat) = (nodes(net), bus.peek(HEAT));
    let refusal = net
        .step(Time::ZERO, dt, bus)
        .expect_err(&format!("{why} must be refused"));
    assert_eq!(nodes(net), before, "{why}: a refused step moved a node");
    left_waiting(why, bus, heat);
    assert_eq!(
        net.absorbed_energy().to_si(),
        0.0,
        "{why}: counted heat taken"
    );
    assert_eq!(net.lost_energy().to_si(), 0.0, "{why}: counted heat lost");
    refusal
}

/// The same for a network, asked of the first node when the mistake is on the second.
fn every_network_question(why: &str, net: &ThermalNetwork, step: &Violation) {
    let a = node(net, "a");
    let got = net.steady_state(Power::w(1.0)).unwrap_err();
    identical(why, "steady_state", &got, step);
    // Whatever it is asked: a power that is not a number was refused first, as that.
    let got = net.steady_state(Power::from_si(NAN)).unwrap_err();
    identical(why, "steady_state(NaN W)", &got, step);
    let got = net.path_conductance(a, Power::w(1.0)).unwrap_err();
    identical(why, "path_conductance", &got, step);
    identical(
        why,
        "time_constant",
        &net.time_constant(a).unwrap_err(),
        step,
    );
    identical(why, "biot_number", &net.biot_number(a).unwrap_err(), step);
}

/// The bar may have been [`warmed`], so what it had absorbed before is what it must have after.
fn refused_bar(why: &str, b: &mut Bar1D, bus: &mut Exchange, dt: Time) -> Violation {
    let (before, heat) = (cells_of_bar(b), bus.peek(HEAT));
    let absorbed = b.absorbed_energy().to_si();
    let refusal = b
        .step(Time::ZERO, dt, bus)
        .expect_err(&format!("{why} must be refused"));
    assert_eq!(
        cells_of_bar(b),
        before,
        "{why}: a refused step moved the bar"
    );
    left_waiting(why, bus, heat);
    assert_eq!(
        b.absorbed_energy().to_si(),
        absorbed,
        "{why}: counted heat taken"
    );
    refusal
}

/// Step a block that must be refused, check that nothing moved, and that `steady_state` refuses as
/// the step did.
fn refused_block(why: &str, b: &mut Solid3D, bus: &mut Exchange, dt: Time) -> Violation {
    let (before, heat) = (cells_of_block(b), bus.peek(HEAT));
    let refusal = b
        .step(Time::ZERO, dt, bus)
        .expect_err(&format!("{why} must be refused"));
    assert_eq!(
        cells_of_block(b),
        before,
        "{why}: a refused step moved the block"
    );
    left_waiting(why, bus, heat);
    assert_eq!(
        b.absorbed_energy().to_si(),
        0.0,
        "{why}: counted heat taken"
    );
    assert_eq!(b.lost_energy().to_si(), 0.0, "{why}: counted heat lost");
    assert_eq!(
        b.generated_energy().to_si(),
        0.0,
        "{why}: counted heat made"
    );
    refusal
}

/// Zero, less than nothing, not a number, and without limit: every way a size can be none.
fn sizes_no_body_has(normal: f64) -> [(&'static str, f64); 4] {
    [
        ("zero", 0.0),
        ("negative", -normal),
        ("NaN", NAN),
        ("+inf", INF),
    ]
}

const STEP: Time = Time::from_si(1e-3);

#[test]
fn a_size_no_body_can_have_is_refused_by_name_before_anything_moves() {
    for (how, v) in sizes_no_body_has(1e-6) {
        let why = format!("lump volume {how}");
        let mut m = lump(fine(), still(), v, 1e-3);
        let refusal = refused_lump(&why, &mut m, &mut waiting(WAITING), STEP);
        names(&why, &refusal, "lump", "the body's volume", v);
        every_lump_question(&why, &m, &refusal);
    }
    for (how, v) in sizes_no_body_has(1e-3) {
        let why = format!("lump thickness {how}");
        let mut m = lump(fine(), still(), 1e-6, v);
        let refusal = refused_lump(&why, &mut m, &mut waiting(WAITING), STEP);
        names(&why, &refusal, "lump", "the body's thickness", v);
        every_lump_question(&why, &m, &refusal);
    }
    for (how, v) in sizes_no_body_has(1e-3) {
        let why = format!("bar cell length {how}");
        let refusal = refused_bar(&why, &mut bar(v, 1e-6), &mut waiting(WAITING), STEP);
        names(&why, &refusal, "bar", "the bar's cell length", v);
    }
    for (how, v) in sizes_no_body_has(1e-6) {
        let why = format!("bar cross-section {how}");
        let refusal = refused_bar(&why, &mut bar(1e-3, v), &mut waiting(WAITING), STEP);
        names(&why, &refusal, "bar", "the bar's cross-section", v);
    }
    for (how, v) in sizes_no_body_has(1e-6) {
        let why = format!("bar face area {how}");
        let mut b = bar(1e-3, 1e-6).exposing("side", Area::from_si(v));
        let refusal = refused_bar(&why, &mut b, &mut waiting(WAITING), STEP);
        names(&why, &refusal, "bar", "the exposed face's area", v);
    }
    for (how, v) in sizes_no_body_has(1e-3) {
        let why = format!("block cell size {how}");
        let mut b = block(v);
        let refusal = refused_block(&why, &mut b, &mut waiting(WAITING), STEP);
        names(&why, &refusal, "block", "the block's cell size", v);
        identical(
            &why,
            "steady_state",
            &b.steady_state().unwrap_err(),
            &refusal,
        );
    }
    for (how, v) in sizes_no_body_has(1e-6) {
        let why = format!("second node's volume {how}");
        let mut net = network(fine(), still(), v, 1e-3);
        let refusal = refused_network(&why, &mut net, &mut waiting(WAITING), STEP);
        names(&why, &refusal, "net/b", "the node's volume", v);
        every_network_question(&why, &net, &refusal);
    }
    for (how, v) in sizes_no_body_has(1e-3) {
        let why = format!("second node's thickness {how}");
        let mut net = network(fine(), still(), 1e-6, v);
        let refusal = refused_network(&why, &mut net, &mut waiting(WAITING), STEP);
        names(&why, &refusal, "net/b", "the node's thickness", v);
        every_network_question(&why, &net, &refusal);
    }
    // And on the interior node, which has no film for a Biot number to read it against: a number
    // nobody can have is refused wherever it was given.
    for (how, v) in sizes_no_body_has(1e-3) {
        let why = format!("interior node's thickness {how}");
        let mut net = network_of(v, fine(), still(), 1e-6, 1e-3);
        let refusal = refused_network(&why, &mut net, &mut waiting(WAITING), STEP);
        names(&why, &refusal, "net/a", "the node's thickness", v);
        every_network_question(&why, &net, &refusal);
    }
}

/// **Every question refuses what the step refuses**, for the numbers the step already refused —
/// a substance, a film — and for a substance that says nothing of its heat, which the step refuses
/// as having no heat capacity.
#[test]
fn every_question_refuses_what_the_step_refuses() {
    let lumps: Vec<(&str, LumpedMass, &str, f64)> = vec![
        (
            "emissivity NaN",
            lump(stuff(1000.0, 50.0, 4000.0, NAN), still(), 1e-6, 1e-3),
            "stuff's emissivity",
            NAN,
        ),
        (
            "conductivity -0.6",
            lump(stuff(1000.0, -0.6, 4000.0, 0.5), still(), 1e-6, 1e-3),
            "stuff's conductivity",
            -0.6,
        ),
        (
            "specific heat -4000",
            lump(stuff(1000.0, 50.0, -4000.0, 0.5), still(), 1e-6, 1e-3),
            "stuff's specific heat",
            -4000.0,
        ),
        (
            "film -7",
            lump(fine(), film(-7.0, 1e-6, room()), 1e-6, 1e-3),
            "the film's convection",
            -7.0,
        ),
        (
            "film area NaN",
            lump(fine(), film(10.0, NAN, room()), 1e-6, 1e-3),
            "the film's area",
            NAN,
        ),
        (
            "ambient -10 K",
            lump(
                fine(),
                film(10.0, 1e-6, Temperature::from_si(-10.0)),
                1e-6,
                1e-3,
            ),
            "the film's ambient",
            -10.0,
        ),
        (
            "no thermal properties",
            lump(silent(), still(), 1e-6, 1e-3),
            "substance has no heat capacity",
            INF,
        ),
    ];
    for (why, mut m, named, value) in lumps {
        let refusal = refused_lump(why, &mut m, &mut waiting(WAITING), STEP);
        names(why, &refusal, "lump", named, value);
        every_lump_question(why, &m, &refusal);
    }

    let networks: Vec<(&str, ThermalNetwork, &str, f64)> = vec![
        (
            "second node's emissivity NaN",
            network(stuff(1000.0, 50.0, 4000.0, NAN), still(), 1e-6, 1e-3),
            "stuff's emissivity",
            NAN,
        ),
        (
            "second node's film -7",
            network(fine(), film(-7.0, 1e-6, room()), 1e-6, 1e-3),
            "the film's convection",
            -7.0,
        ),
        // Interior nodes ask nothing of an emissivity, and the second node loses to a film: a
        // network solved its balance as if a node with no thermal properties did not radiate.
        (
            "second node with no thermal properties",
            network(silent(), still(), 1e-6, 1e-3),
            "substance has no heat capacity",
            f64::NAN,
        ),
    ];
    for (why, mut net, named, value) in networks {
        let refusal = refused_network(why, &mut net, &mut waiting(WAITING), STEP);
        names(why, &refusal, "net/b", named, value);
        every_network_question(why, &net, &refusal);
    }

    // A block its step cannot step, for what makes it so: `steady_state` said "this block cannot
    // be stepped" of each, and the step's refusal of a joint carried `+∞`, a value nobody stated.
    let joint =
        |h: f64| block(1e-3).joined(move |axis, i, _, _| (axis == Axis::X && i == 2).then_some(h));
    let blocks: Vec<(&str, Solid3D, &str, f64)> = vec![
        (
            "a block with no thermal properties",
            Solid3D::new("block", silent(), (4, 1, 1), Length::mm(1.0), room())
                .losing_from(Face::XMax, still()),
            "substance has no diffusivity",
            INF,
        ),
        (
            "a joint of -5",
            joint(-5.0),
            "a joint's contact conductance is negative",
            -5.0,
        ),
        (
            "a joint of NaN",
            joint(NAN),
            "a joint's contact conductance is not a number",
            NAN,
        ),
    ];
    for (why, mut b, named, value) in blocks {
        let refusal = refused_block(why, &mut b, &mut waiting(WAITING), STEP);
        names(why, &refusal, "block", named, value);
        identical(
            why,
            "steady_state",
            &b.steady_state().unwrap_err(),
            &refusal,
        );
    }
}

/// **Heat taken out of a body nothing brings heat to cools it without limit**, and is refused
/// saying so, as heat given to one nothing takes heat from warms it: a network with no film said
/// both warm.
#[test]
fn a_balance_nothing_can_strike_is_refused_the_way_it_goes() {
    let mut net = ThermalNetwork::new("net");
    let a = net.node("a", fine(), Volume::from_si(1e-6), Length::mm(1.0), room());
    net.absorbing(a).expect("its own node");
    let sealed = lump(fine(), film(0.0, 0.0, room()), 1e-6, 1e-3);
    for (watts, how) in [(1.0, "warms without limit"), (-1.0, "cools without limit")] {
        let why = format!("{watts} W, with no film");
        let refusal = net.steady_state(Power::w(watts)).unwrap_err();
        names(&why, &refusal, "net", how, watts);
        let refusal = sealed.equilibrium_rise(Power::w(watts)).unwrap_err();
        names(&why, &refusal, "lump", how, watts);
    }
}

#[test]
fn a_source_that_is_not_a_number_is_refused_by_the_cell_it_is_in() {
    for watts in [NAN, INF, -INF] {
        // A good source first, so a check that read only the first cell with one would pass.
        let why = format!("source {watts} in the third cell");
        let mut b = block(1e-3)
            .dissipating(1e-3, |i, _, _| i == 0)
            .dissipating(watts, |i, _, _| i == 2);
        let refusal = refused_block(&why, &mut b, &mut waiting(WAITING), STEP);
        names(
            &why,
            &refusal,
            "block",
            "the source in cell (2, 0, 0)",
            watts,
        );
        identical(
            &why,
            "steady_state",
            &b.steady_state().unwrap_err(),
            &refusal,
        );
    }

    // Named by all three of its indices, in a block whose counts and whose cell's indices all
    // differ, so no two of them can be swapped unseen.
    let why = "source NaN at (1, 2, 3)";
    let mut b = Solid3D::new("block", fine(), (2, 3, 4), Length::mm(1.0), room())
        .losing_from(Face::XMax, still())
        .dissipating(NAN, |i, j, k| (i, j, k) == (1, 2, 3));
    let refusal = refused_block(why, &mut b, &mut waiting(WAITING), STEP);
    names(why, &refusal, "block", "the source in cell (1, 2, 3)", NAN);

    // The first in the cells' order, not the last one written.
    let why = "+inf at the second cell, then NaN at the fourth";
    let mut b = block(1e-3)
        .dissipating(INF, |i, _, _| i == 1)
        .dissipating(NAN, |i, _, _| i == 3);
    let refusal = refused_block(why, &mut b, &mut waiting(WAITING), STEP);
    names(why, &refusal, "block", "the source in cell (1, 0, 0)", INF);

    // A cell emptied after its source was written is still where the mistake was made, and a
    // source written after that, elsewhere, does not forget it.
    let why = "NaN in a cell emptied after";
    let mut b = block(1e-3)
        .dissipating(NAN, |i, _, _| i == 3)
        .empty(|i, _, _| i == 3)
        .dissipating(1e-3, |i, _, _| i == 0);
    let refusal = refused_block(why, &mut b, &mut waiting(WAITING), STEP);
    names(why, &refusal, "block", "the source in cell (3, 0, 0)", NAN);

    // Written over with a number, it is a number: the block runs, and generates what it was told.
    let mut b = block(1e-3)
        .dissipating(NAN, |i, _, _| i == 2)
        .dissipating(1e-3, |i, _, _| i == 2);
    let dt = a_step(&b);
    b.step(Time::ZERO, dt, &mut Exchange::new())
        .expect("a source written over with a number runs");
    assert_eq!(b.generated_energy().to_si(), 1e-3 * dt.to_si());
    b.steady_state()
        .expect("and has a balance, its film carrying the milliwatt away");
}

#[test]
fn heat_that_is_not_a_number_is_refused_before_any_of_it_is_taken() {
    let what = "the heat waiting on the bus";
    for joules in [NAN, INF, -INF] {
        let why = format!("{joules} J on the channel");
        let refusal = refused_lump(&why, &mut a_lump(), &mut waiting(joules), STEP);
        names(&why, &refusal, "lump", what, joules);
        let refusal = refused_bar(
            &why,
            &mut warmed(bar(1e-3, 1e-6)),
            &mut waiting(joules),
            STEP,
        );
        names(&why, &refusal, "bar", what, joules);
        let refusal = refused_block(&why, &mut block(1e-3), &mut waiting(joules), STEP);
        names(&why, &refusal, "block", what, joules);
        let refusal = refused_network(&why, &mut a_network(), &mut waiting(joules), STEP);
        names(&why, &refusal, "net", what, joules);
    }

    // On a bar's faces, by the face, and left where it was offered, with heat waiting on the
    // channel beside it. A good face first, and a second bad one of another kind after, so the
    // refusal can be told from the faces' sum and the first bad face from the last.
    for joules in [NAN, INF, -INF] {
        let why = format!("{joules} J on the third face");
        let mut b = warmed(bar(1e-3, 1e-6).exposing("side", Area::from_si(1e-6)));
        let side = b.boundary().expect("exposed").clone();
        let mut faces = vec![0.0; 8];
        faces[0] = 1e-4;
        faces[2] = joules;
        faces[5] = if joules.is_nan() { INF } else { NAN };
        let mut bus = waiting(WAITING);
        bus.publish_on(&side, HEAT, &Flux::from_faces(faces.clone()))
            .expect("one face per cell");
        let refusal = refused_bar(&why, &mut b, &mut bus, STEP);
        names(&why, &refusal, "bar", "the heat arriving on face 2", joules);
        let offered = bus.peek_on(&side, HEAT).expect("still on the faces");
        let still_there = offered
            .per_face()
            .iter()
            .zip(&faces)
            .all(|(a, b)| same(*a, *b));
        assert!(still_there, "{why}: a refused step took heat off the faces");
    }

    // On another grid it is refused as another grid, by the faces' count: looked at face by face
    // it named a face the bar does not have.
    let why = "NaN J on face 40 of a 64-face grid";
    let mut b = warmed(bar(1e-3, 1e-6).exposing("side", Area::from_si(1e-6)));
    let camera = pantometry_core::Interface::uniform("side", 64, Area::from_si(1e-7));
    let mut faces = vec![1e-6; 64];
    faces[40] = NAN;
    let mut bus = waiting(WAITING);
    bus.publish_on(&camera, HEAT, &Flux::from_faces(faces))
        .expect("the camera's own grid");
    let refusal = refused_bar(why, &mut b, &mut bus, STEP);
    assert!(
        refusal.quantity.contains("expected 8"),
        "{why}: refused for {refusal}, not for its grid"
    );
    assert!(
        bus.peek_on(&camera, HEAT).is_some(),
        "{why}: a refused step took the flux"
    );

    // Heat with nowhere to land, refused before the take: it was taken and counted first.
    let mut net = ThermalNetwork::new("net");
    net.node_losing_to(
        "a",
        fine(),
        Volume::from_si(1e-6),
        Length::mm(1.0),
        Temperature::celsius(100.0),
        still(),
    );
    let why = "heat for a network with no node to absorb it";
    let refusal = refused_network(why, &mut net, &mut waiting(WAITING), STEP);
    names(
        why,
        &refusal,
        "net",
        "no node was named to absorb it",
        WAITING,
    );
}

#[test]
fn a_step_no_body_can_take_is_refused() {
    let what = "the step must be finite and zero or more seconds";
    for seconds in [NAN, -1e-3, -INF, INF] {
        let why = format!("a step of {seconds} s");
        let dt = Time::from_si(seconds);
        let refusal = refused_lump(&why, &mut a_lump(), &mut waiting(WAITING), dt);
        names(&why, &refusal, "lump", what, seconds);
        let refusal = refused_bar(
            &why,
            &mut warmed(bar(1e-3, 1e-6)),
            &mut waiting(WAITING),
            dt,
        );
        names(&why, &refusal, "bar", what, seconds);
        let refusal = refused_block(&why, &mut block(1e-3), &mut waiting(WAITING), dt);
        names(&why, &refusal, "block", what, seconds);
        let refusal = refused_network(&why, &mut a_network(), &mut waiting(WAITING), dt);
        names(&why, &refusal, "net", what, seconds);
    }
}

/// **A step of no length is a step**, and each domain takes it as it always did: a lumped body, a
/// bar and a block take the whole channel, which is what [`Exchange::take_share`] gives a step it
/// is not told an interval for, and add it with nothing lost or conducted; a network takes nothing
/// and moves nothing. Minus zero is zero. And it still refuses a body no body can be.
#[test]
fn a_step_of_no_length_runs() {
    for seconds in [0.0, -0.0] {
        let dt = Time::from_si(seconds);
        let mut m = a_lump();
        let capacity = m.heat_capacity().to_si();
        let start = m.temperature().to_si();
        m.step(Time::ZERO, dt, &mut waiting(WAITING))
            .expect("a lumped body takes a step of no length");
        assert_eq!(m.temperature().to_si(), start + WAITING / capacity);

        let mut b = bar(1e-3, 1e-6);
        let start = b.temperature_at(0).to_si();
        let mut bus = waiting(WAITING);
        b.step(Time::ZERO, dt, &mut bus)
            .expect("a bar takes a step of no length");
        assert_eq!(bus.peek(HEAT), 0.0, "and the whole channel with it");
        assert!(b.temperature_at(0).to_si() > start);
        assert_eq!(b.temperature_at(7).to_si(), start, "with nothing conducted");

        let mut k = block(1e-3);
        let start = cells_of_block(&k);
        let mut bus = waiting(WAITING);
        k.step(Time::ZERO, dt, &mut bus)
            .expect("a block takes a step of no length");
        assert_eq!(bus.peek(HEAT), 0.0, "and the whole channel with it");
        assert_ne!(cells_of_block(&k), start);
        assert_eq!(k.lost_energy().to_si(), 0.0, "with nothing lost");

        let mut net = a_network();
        let start = nodes(&net);
        let mut bus = waiting(WAITING);
        net.step(Time::ZERO, dt, &mut bus)
            .expect("a network takes a step of no length");
        assert_eq!(bus.peek(HEAT), WAITING, "and leaves the channel as it was");
        assert_eq!(nodes(&net), start);

        // And a step of no length still refuses a body that cannot be one. A network returned
        // from a step of no length before it asked anything, as success.
        let why = format!("a step of {seconds} s of a node of no volume");
        let mut net = network(fine(), still(), 0.0, 1e-3);
        let refusal = refused_network(&why, &mut net, &mut waiting(WAITING), dt);
        names(&why, &refusal, "net/b", "the node's volume", 0.0);
    }
}

/// The body the balance is asked of: a cubic centimetre of `fine()` under a film of `h` over
/// `area`, in a room at 20 °C, with emissivity `emissivity`.
fn settling(h: f64, area: f64, emissivity: f64) -> LumpedMass {
    LumpedMass::new(
        "lump",
        stuff(1000.0, 50.0, 4000.0, emissivity),
        Volume::from_si(1e-6),
        Length::mm(1.0),
        Temperature::celsius(20.0),
        film(h, area, Temperature::celsius(20.0)),
    )
}

/// `(x·x)·(x·x)`, the arithmetic the balance is solved and bounded in, so a power at the bound to
/// the bit is the bound the body was asked against.
fn fourth(x: f64) -> f64 {
    let s = x * x;
    s * s
}

/// What `hA·ΔT + εσA((Tₐ + ΔT)⁴ − Tₐ⁴) − P` comes to, in watts, and the sum of its four terms'
/// sizes, which is what its rounding is a fraction of.
fn residual(h: f64, area: f64, emissivity: f64, rise: f64, p: f64) -> (f64, f64) {
    let ta = Temperature::celsius(20.0).to_si();
    let er = emissivity * STEFAN_BOLTZMANN.to_si() * area;
    let terms = [
        h * area * rise,
        er * fourth(ta + rise),
        -er * fourth(ta),
        -p,
    ];
    (terms.iter().sum(), terms.iter().map(|x| x.abs()).sum())
}

/// The rounding of a sum of four terms, each a few operations deep, as a fraction of their sizes:
/// sixteen units in the last place. The most measured on these bodies was 0.37 of one.
const BALANCE_ROUNDING: f64 = 16.0 * f64::EPSILON;

/// `rise` is a root of the balance: what is left of it is within the rounding of its terms.
fn is_the_root(why: &str, h: f64, area: f64, emissivity: f64, rise: f64, p: f64) {
    let (off, size) = residual(h, area, emissivity, rise, p);
    assert!(
        off.abs() <= BALANCE_ROUNDING * size,
        "{why}: answered {rise} K, where the balance is off by {off} W of terms of {size} W"
    );
}

/// How close a march comes to the root it settles on: each step closes a tenth of the gap, so it
/// stops moving once a tenth of the gap rounds away — within five units in the last place of the
/// temperature, and eight is that with the root's own rounding. The most measured was 4.12.
const MARCH_ROUNDING: f64 = 8.0 * f64::EPSILON;

/// **A balance below absolute zero has none**, and one above it is the balance: the root of the
/// loss, and where the body marches to.
///
/// The most the room brings in is to a body at absolute zero, `hA·Tₐ + εσA·Tₐ⁴`, computed here
/// rather than read. A part in a billion past it is refused, and so is one unit in the last place;
/// at it to the bit, a unit short of it, a part in a million short of it, at half of it and given a
/// watt, the answer is at or above absolute zero and is the root of the balance to the rounding of
/// the balance's own terms; and at half of it and given a watt, it is where `step` settles with the
/// heat taken out or put in over the bus.
///
/// On each loss path alone, so each start the solve takes is met, and under films that carry a
/// sliver beside radiation, where a start from convection ran out of steps: under 1e-9
/// W·m⁻²·K⁻¹ a body that settles 46.64 K below its room was answered 1115.6 K above it, and given a
/// watt, 9797.6 K for 481.5 K. At the bound the answer came out a rounding below absolute zero,
/// −5.7e-14 K under a film of 1 and −1.1e-8 K under one of 1e-6, and on convection alone over a
/// square millimetre under a film of 1e-7, −5.7e-14 K.
#[test]
fn a_balance_below_absolute_zero_is_refused_and_one_above_it_is_answered() {
    let ta = Temperature::celsius(20.0).to_si();
    let sigma = STEFAN_BOLTZMANN.to_si();
    let bodies = [
        (10.0, 0.5, 1e-4),
        (10.0, 0.0, 1e-4),
        (1e-7, 0.0, 1e-6),
        (0.0, 0.5, 1e-4),
        (1.0, 0.5, 1e-4),
        (1e-6, 0.5, 1e-4),
        (1e-7, 0.5, 1e-4),
        (1e-9, 0.5, 1e-4),
    ];
    for (h, emissivity, area) in bodies {
        let body = settling(h, area, emissivity);
        let most = h * area * ta + emissivity * sigma * area * fourth(ta);
        let why = format!("h {h:e}, emissivity {emissivity}, area {area:e}");

        // Negative numbers grow away from zero as their bits do.
        let (at, ulp_past, ulp_short) = (
            -most,
            f64::from_bits((-most).to_bits() + 1),
            f64::from_bits((-most).to_bits() - 1),
        );
        for past in [-most * (1.0 + 1e-9), ulp_past] {
            let refusal = body.equilibrium_rise(Power::from_si(past)).unwrap_err();
            names(
                &why,
                &refusal,
                "lump",
                "no balance above absolute zero",
                past,
            );
        }

        for p in [at, ulp_short, -most * (1.0 - 1e-6), -0.5 * most, 1.0] {
            let rise = body.equilibrium_rise(Power::from_si(p)).unwrap().to_si();
            assert!(ta + rise >= 0.0, "{why}, {p} W: answered {} K", ta + rise);
            is_the_root(&format!("{why}, {p} W"), h, area, emissivity, rise, p);
        }

        // Marched there with the heat over the bus: explicit Euler's fixed point is the balance
        // itself, so the march comes to the root whatever its step.
        for p in [-0.5 * most, 1.0] {
            let want = body.equilibrium_rise(Power::from_si(p)).unwrap().to_si();
            let mut marched = settling(h, area, emissivity);
            let mut bus = Exchange::new();
            for _ in 0..4000 {
                let dt = marched.max_stable_dt(Time::ZERO);
                bus.publish(HEAT, p * dt.to_si());
                marched.step(Time::ZERO, dt, &mut bus).unwrap();
            }
            let got = marched.rise().to_si();
            assert!(
                (got - want).abs() <= MARCH_ROUNDING * (ta + want),
                "{why}, {p} W: marched to {got} K, the balance is {want} K"
            );
        }
    }

    // A power that is not a number of watts.
    for watts in [NAN, INF, -INF] {
        let why = format!("{watts} W");
        let refusal = settling(10.0, 1e-4, 0.5)
            .equilibrium_rise(Power::from_si(watts))
            .unwrap_err();
        names(
            &why,
            &refusal,
            "lump",
            "the absorbed power must be finite",
            watts,
        );
    }
}

/// **A Biot number where nothing convects is zero**, whatever conducts: it was `0/0` for a perfect
/// insulator in still air. Under a film the insulator's is `+∞`, a body no lumped model describes.
#[test]
fn a_biot_number_where_nothing_convects_is_zero() {
    let insulator = stuff(1000.0, 0.0, 4000.0, 0.5);
    let calm = lump(insulator.clone(), film(0.0, 1e-6, room()), 1e-6, 1e-3);
    assert_eq!(calm.biot_number().unwrap(), 0.0);
    let filmed = lump(insulator, still(), 1e-6, 1e-3);
    assert_eq!(filmed.biot_number().unwrap(), INF);
    // And the closed form where it is one: hL/k.
    assert_eq!(a_lump().biot_number().unwrap(), 10.0 * 1e-3 / 50.0);

    // A network says `None` for a node with no film, and the Biot number of one with a film: the
    // same number the lumped body gives, an insulator's included, which was `None`.
    let net = a_network();
    assert_eq!(
        net.biot_number(node(&net, "a")).unwrap(),
        None,
        "an interior node"
    );
    assert_eq!(
        net.biot_number(node(&net, "b")).unwrap(),
        Some(10.0 * 1e-3 / 50.0)
    );
    let insulator = stuff(1000.0, 0.0, 4000.0, 0.5);
    let filmed = network(insulator.clone(), still(), 1e-6, 1e-3);
    assert_eq!(
        filmed.biot_number(node(&filmed, "b")).unwrap(),
        Some(INF),
        "an insulator under a film"
    );
    let calm = network(insulator, film(0.0, 1e-6, room()), 1e-6, 1e-3);
    assert_eq!(
        calm.biot_number(node(&calm, "b")).unwrap(),
        Some(0.0),
        "an insulator where nothing convects"
    );
}

/// **Through a simulation the refusal is the body's**, by name, and comes out of
/// `Simulation::advance` as the step gave it. Whatever limit a refused body reports, the kernel
/// steps it once — a limit of zero, `NaN` or less plans one substep — so what this holds is the
/// step's own refusal on the way out; that a refused lumped body's limit is `+∞` is held where
/// every question is asked of it.
#[test]
fn through_the_simulation_the_refusal_is_the_bodys() {
    let run = |why: &str, d: Box<dyn Domain>, site: &str, named: &str, value: f64| {
        let mut sim = Simulation::new(Schedule::Multirate).with_boxed(d);
        let refusal = sim
            .advance(Time::from_si(0.1))
            .expect_err(&format!("{why} must be refused at the first advance"));
        names(why, &refusal, site, named, value);
    };
    run(
        "lump volume NaN",
        Box::new(lump(fine(), still(), NAN, 1e-3)),
        "lump",
        "the body's volume",
        NAN,
    );
    run(
        "lump volume -1e-6",
        Box::new(lump(fine(), still(), -1e-6, 1e-3)),
        "lump",
        "the body's volume",
        -1e-6,
    );
    run(
        "bar cell length NaN",
        Box::new(bar(NAN, 1e-6)),
        "bar",
        "the bar's cell length",
        NAN,
    );
    run(
        "block cell size 0",
        Box::new(block(0.0)),
        "block",
        "the block's cell size",
        0.0,
    );
    run(
        "block source NaN",
        Box::new(block(1e-3).dissipating(NAN, |i, _, _| i == 1)),
        "block",
        "the source in cell (1, 0, 0)",
        NAN,
    );
    run(
        "second node's volume NaN",
        Box::new(network(fine(), still(), NAN, 1e-3)),
        "net/b",
        "the node's volume",
        NAN,
    );
}
