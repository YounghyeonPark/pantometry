//! Heat off the bus is added after the sweep, so the step that brings it does not also shed it.
//!
//! A forward step is `T′ = T + (dt/C)(P + F(T))`: what arrives and what the operator does are both
//! the step's, and the operator is evaluated at the state the step began with. `Solid3D` put the
//! bus's heat in first, before the sweep read its state, so the film met the heated state and shed
//! a share `θ = dt·G/C` of the heat in the step that brought it. The books balanced throughout,
//! because the heat really did leave. Two things followed:
//!
//! - **The steady state depended on the step.** One cell convecting at `h = 100` under constant
//!   bus power settled at 0.900, 0.500 and 0.100 of `P/G` above the air at a tenth, a half and
//!   nine tenths of its limit — and at its limit it stayed at the air.
//! - **The step outran its limit.** The limit is taken from the state the step begins in, and a
//!   radiating film met at a temperature the heat had raised is stiffer than that: a black cell of
//!   1 mm and `k = 200` at 1000 K, stepped at its own limit while the bus raised it by 600 K, landed
//!   at −33.5 K, and by 1000 K at −1982 K, each step accepted.
//!
//! `dissipating`'s heat had been through the first of these already and went in after the sweep
//! for it; the bus's had not.
//!
//! # What is checked against what
//!
//! **Constant bus power, steady.** At every step up to the limit, the march settles where the
//! scheme's own balance puts it, which no step size enters: one cell convecting, at `P/G` above the
//! air, `G` the film behind its half cell; one cell radiating, where the surface sends `P` and the
//! half cell carries it, found by bisection; and a column under a film, the heat spread evenly over
//! it as the bus spreads it, at the parabola plus `q'''dx²/8k` at every cell.
//!
//! **A slug.** One step that brings heat is the same step without it, with the rise added: to
//! rounding, at the limit, however large the slug — and so never below the room. In a black cell
//! of two conductivities, and in a cell of water, whose substance melts and so takes the block's
//! other sweep, in place through the enthalpy.
//!
//! **The phase it leaves.** Heat that melts a cell in the step it arrives in must do so before the
//! step re-derives the block, so that the next step reads the phase it left: after a step that
//! melts half-melted ice to nine tenths, the block's limit and a second step's loss are a block
//! set to that phase from the start's, to rounding.

use pantometry_core::substance::ThermalProps;
use pantometry_core::units::{
    Area, Density, Length, SpecificHeat, Temperature, ThermalConductivity, ThermalExpansion, Time,
};
use pantometry_core::{Domain, Exchange, Substance};
use pantometry_thermal::{Environment, Face, Solid3D, HEAT};

/// The Stefan–Boltzmann constant, W/m²·K⁴, as CODATA 2018 states it exactly.
const SIGMA: f64 = 5.670_374_419e-8;
const RHO: f64 = 1000.0;
const C: f64 = 4000.0;

fn stuff(k: f64, emissivity: f64) -> Substance {
    Substance::bulk("stuff", Density::from_si(RHO)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(k),
        specific_heat: SpecificHeat::j_per_kg_k(C),
        expansion: ThermalExpansion::ppm_per_k(50.0),
        emissivity,
    })
}

fn fourth(t: f64) -> f64 {
    let s = t * t;
    s * s
}

/// Whether a residual is out: at or past `tolerance`, or not a number.
fn off(residual: f64, tolerance: f64) -> bool {
    residual.is_nan() || residual >= tolerance
}

/// `block` with a film on `face`, convecting at `h` to `room`, over `area`.
fn filmed(block: Solid3D, face: Face, h: f64, room: f64, area: f64) -> Solid3D {
    block.losing_from(
        face,
        Environment {
            ambient: Temperature::from_si(room),
            convection_w_per_m2_k: h,
            area: Area::from_si(area),
        },
    )
}

/// March `block` under `watts` off the bus, each step `fraction` of its own limit read before it,
/// for `steps` steps, and return its temperatures.
fn settle(mut block: Solid3D, watts: f64, fraction: f64, steps: usize) -> Result<Solid3D, String> {
    for step in 0..steps {
        let dt = fraction * block.max_stable_dt(Time::ZERO).to_si();
        let mut bus = Exchange::new();
        bus.publish(HEAT, watts * dt);
        block
            .step(Time::ZERO, Time::from_si(dt), &mut bus)
            .map_err(|v| format!("step {step} at {fraction} of its limit was refused: {v:?}"))?;
    }
    Ok(block)
}

/// **Constant power off the bus settles a block where its balance puts it, at any step.**
#[test]
fn bus_power_settles_where_the_balance_does_at_every_step() {
    let (dx, room) = (1e-3, 300.0);
    let area = dx * dx;
    let mut wrong = Vec::new();
    // At a tenth of its limit a single cell falls toward its balance by 0.9 a step, and twenty
    // thousand steps leave nothing of its start. What is left is where the march stalls: once a
    // step's update rounds to nothing it moves no more, about an ulp of the temperature over twice
    // the fraction of the slowest mode a step removes. For a single cell that is a few ulps,
    // measured 6.8e-13 K at most; for the column at a tenth of its limit, a thousand steps an
    // e-fold, 2.7e-11 K predicted and 3.7e-11 K measured. 1e-9 K is 27 times the worst of them, and
    // grows smaller against a slower march: at a two-hundred-and-fiftieth of the limit the column
    // would reach it. Against 0.21 K and more where the step shed what it brought.
    let held = 1e-9;
    let fractions = [0.1, 0.5, 0.9, 1.0];

    // One cell convecting: `P/G` above the air, `G` the film behind the half cell.
    let (k, h, watts) = (200.0, 100.0, 1e-3);
    let g = 1.0 / (1.0 / (h * area) + dx / (2.0 * k * area));
    for fraction in fractions {
        let block = filmed(
            Solid3D::new(
                "cell",
                stuff(k, 0.0),
                (1, 1, 1),
                Length::from_si(dx),
                Temperature::from_si(room),
            ),
            Face::XMin,
            h,
            room,
            area,
        );
        match settle(block, watts, fraction, 20_000) {
            Err(e) => wrong.push(format!("a convecting cell: {e}")),
            Ok(b) => {
                let rise = b.temperature_at(0, 0, 0).to_si() - room;
                println!(
                    "  a convecting cell at {fraction} of its limit: {rise:.12} K over {:.12}",
                    watts / g
                );
                if off((rise - watts / g).abs(), held) {
                    wrong.push(format!(
                        "a convecting cell at {fraction} of its limit settles {rise:.9} K above \
                         the air, where P/G is {:.9} K",
                        watts / g
                    ));
                }
            }
        }
    }

    // One black cell radiating: its surface where it sends `P`, its centre `P` over the half cell
    // above that.
    let (k, watts) = (0.6, 0.02);
    let half = 2.0 * k * area / dx;
    let (mut lo, mut hi) = (room, 5000.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if SIGMA * area * (fourth(mid) - fourth(room)) > watts {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let centre = 0.5 * (lo + hi) + watts / half;
    for fraction in fractions {
        let block = filmed(
            Solid3D::new(
                "cell",
                stuff(k, 1.0),
                (1, 1, 1),
                Length::from_si(dx),
                Temperature::from_si(room),
            ),
            Face::XMin,
            0.0,
            room,
            area,
        );
        match settle(block, watts, fraction, 20_000) {
            Err(e) => wrong.push(format!("a radiating cell: {e}")),
            Ok(b) => {
                let t = b.temperature_at(0, 0, 0).to_si();
                println!("  a radiating cell at {fraction} of its limit: {t:.12} K against {centre:.12} K");
                if off((t - centre).abs(), held) {
                    wrong.push(format!(
                        "a radiating cell at {fraction} of its limit settles at {t:.9} K, where \
                         its surface sends what arrives at {centre:.9} K"
                    ));
                }
            }
        }
    }

    // A column of eight cells under a film at its far end, the heat spread over it as the bus
    // spreads it: the parabola from the film's surface, plus `q'''dx²/8k` at every cell.
    let (n, k, h, watts) = (8usize, 0.6, 200.0, 2e-2);
    let depth = n as f64 * dx;
    let heating = watts / (depth * area);
    let surface = room + watts / (h * area);
    for fraction in fractions {
        let block = filmed(
            Solid3D::new(
                "column",
                stuff(k, 0.0),
                (n, 1, 1),
                Length::from_si(dx),
                Temperature::from_si(room),
            ),
            Face::XMax,
            h,
            room,
            area,
        );
        // The column's slowest mode decays over 318 s, a thousand steps at a tenth of its limit:
        // twenty thousand left 1.8e-7 K of its start, and sixty thousand leave nothing but the
        // stall above.
        match settle(block, watts, fraction, 60_000) {
            Err(e) => wrong.push(format!("a column: {e}")),
            Ok(b) => {
                let mut worst: f64 = 0.0;
                for i in 0..n {
                    let x = (i as f64 + 0.5) * dx;
                    let want = surface
                        + heating / (2.0 * k) * (depth * depth - x * x)
                        + heating * dx * dx / (8.0 * k);
                    let miss = b.temperature_at(i, 0, 0).to_si() - want;
                    if miss.is_nan() || miss.abs() > worst.abs() {
                        worst = miss;
                    }
                }
                println!(
                    "  a column at {fraction} of its limit: its worst cell {worst:+.3e} K off"
                );
                if off(worst.abs(), held) {
                    wrong.push(format!(
                        "a column at {fraction} of its limit has a cell {worst:+.3e} K from the \
                         parabola plus q'''dx²/8k"
                    ));
                }
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// **A slug of heat off the bus is the step without it, with the rise added.**
#[test]
fn a_slug_is_added_to_the_step_it_arrives_in() {
    let (dx, room, from) = (1e-3, 300.0, 1000.0);
    let area = dx * dx;
    let capacity = RHO * C * dx * dx * dx;
    let mut wrong = Vec::new();
    let cell = |k: f64| {
        filmed(
            Solid3D::new(
                "cell",
                stuff(k, 1.0),
                (1, 1, 1),
                Length::from_si(dx),
                Temperature::from_si(from),
            ),
            Face::XMin,
            0.0,
            room,
            area,
        )
    };
    for k in [200.0, 0.6] {
        let block = cell(k);
        let limit = Time::from_si(block.max_stable_dt(Time::ZERO).to_si());
        let mut bare = block.clone();
        bare.step(Time::ZERO, limit, &mut Exchange::new())
            .expect("its own limit");
        let without = bare.temperature_at(0, 0, 0).to_si();
        for rise in [300.0, 600.0, 1000.0, 2000.0] {
            let mut slugged = block.clone();
            let mut bus = Exchange::new();
            bus.publish(HEAT, rise * capacity);
            if let Err(v) = slugged.step(Time::ZERO, limit, &mut bus) {
                wrong.push(format!(
                    "k = {k}, a slug of {rise} K: refused at its limit: {v:?}"
                ));
                continue;
            }
            let with = slugged.temperature_at(0, 0, 0).to_si();
            // The rise, added to a temperature near 1000 K and read back: an ulp or two of 3000 K.
            let apart = (with - (without + rise)).abs();
            println!("  k = {k}, a slug of {rise} K: {with:.9} K, the step without it plus the rise {:.9} K", without + rise);
            if off(apart, 1e-11) {
                wrong.push(format!(
                    "k = {k}, a slug of {rise} K: the step lands at {with:.9} K, where the step \
                     without it lands at {without:.9} K"
                ));
            }
            if with < room {
                wrong.push(format!(
                    "k = {k}, a slug of {rise} K: shed to {with} K, below the room"
                ));
            }
        }
    }

    // A block that melts takes its own sweep, in place and through each cell's enthalpy, and the
    // heat has to go in after that one as well: a cell of the water ice melts to, under a film of
    // `h = 10⁴` to a room 73 K below it, every bit of it liquid before and after.
    let water = filmed(
        Solid3D::new(
            "cell",
            ice_and_water(),
            (1, 1, 1),
            Length::mm(1.0),
            Temperature::from_si(373.0),
        ),
        Face::XMin,
        1e4,
        room,
        area,
    );
    let limit = Time::from_si(water.max_stable_dt(Time::ZERO).to_si());
    let capacity = water.heat_capacity().to_si();
    let mut bare = water.clone();
    bare.step(Time::ZERO, limit, &mut Exchange::new())
        .expect("its own limit");
    let without = bare.temperature_at(0, 0, 0).to_si();
    for rise in [100.0, 300.0] {
        let mut slugged = water.clone();
        let mut bus = Exchange::new();
        bus.publish(HEAT, rise * capacity);
        if let Err(v) = slugged.step(Time::ZERO, limit, &mut bus) {
            wrong.push(format!(
                "water, a slug of {rise} K: refused at its limit: {v:?}"
            ));
            continue;
        }
        let with = slugged.temperature_at(0, 0, 0).to_si();
        println!(
            "  water, a slug of {rise} K: {with:.9} K, the step without it plus the rise {:.9} K",
            without + rise
        );
        if off((with - (without + rise)).abs(), 1e-11) {
            wrong.push(format!(
                "water, a slug of {rise} K: the step lands at {with:.9} K, where the step without \
                 it lands at {without:.9} K"
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Ice with water named as its liquid, neither radiating: water conducts 0.598 W/m·K against the
/// ice's 2.22, so a cell's half cell and capacity move a long way with its phase.
fn ice_and_water() -> Substance {
    let base = Substance::ice().with_emissivity(0.0);
    let fusion = base.fusion.expect("ice melts");
    base.with_fusion(fusion.with_liquid(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(0.598),
        specific_heat: SpecificHeat::j_per_kg_k(4182.0),
        expansion: ThermalExpansion::ppm_per_k(207.0),
        emissivity: 0.0,
    }))
}

/// **The phase the bus heat leaves is the one the next step reads.**
///
/// The heat goes in after the sweep and before the step re-derives its phases, and that second
/// half of the order matters as much as the first: added after the re-derivation, a cell the heat
/// melted would be read by the next step as the phase it was, its limit and its film a step
/// behind. One cell of half-melted ice under a film at `h = 10⁴`, fed four tenths of its latent
/// heat in a step: after it, the block's limit and what a second step loses are a block's that was
/// set to the phase the first step left, to the rounding of two derivations of the same numbers.
#[test]
fn the_phase_the_bus_heat_leaves_is_the_one_the_next_step_reads() {
    let room = Temperature::celsius(-10.0).to_si();
    let build = |fraction: f64| {
        let mut block = filmed(
            Solid3D::new(
                "cell",
                ice_and_water(),
                (1, 1, 1),
                Length::mm(1.0),
                Temperature::celsius(0.0),
            ),
            Face::XMin,
            1e4,
            room,
            1e-6,
        );
        block.set_melted_fraction(0, 0, 0, fraction);
        block
    };
    let latent = Substance::ice().density.to_si()
        * 1e-9
        * Substance::ice()
            .fusion
            .expect("ice melts")
            .latent_heat
            .to_si();

    let mut fed = build(0.5);
    let dt = Time::from_si(0.25 * fed.max_stable_dt(Time::ZERO).to_si());
    let mut bus = Exchange::new();
    bus.publish(HEAT, 0.4 * latent);
    fed.step(Time::ZERO, dt, &mut bus)
        .expect("a quarter of its limit");
    let left = fed.melted_fraction_at(0, 0, 0);
    let fresh = build(left);

    // What a block a step behind would read: the phase it began in.
    let behind = build(0.5).max_stable_dt(Time::ZERO).to_si();
    let (now, set) = (
        fed.max_stable_dt(Time::ZERO).to_si(),
        fresh.max_stable_dt(Time::ZERO).to_si(),
    );
    println!(
        "  melted {left:.12} by the step; the limit after it {now:.12e} s, set to that phase \
         {set:.12e} s, a step behind {behind:.12e} s"
    );
    assert!(
        left > 0.8 && left < 0.95,
        "the step was to melt the cell to about nine tenths, and left it {left}"
    );
    // Two derivations of the same mixed numbers, a handful of roundings apart; and the phase a
    // step behind is the half-melted one, which conducts 1.8 times as much.
    assert!(
        (now / set - 1.0).abs() < 1e-14,
        "after the step the limit is {now:e} s, and a block set to the phase it left has {set:e} s"
    );
    assert!(
        (behind / set - 1.0).abs() > 1e-3,
        "the phase a step behind gives {behind:e} s, too close to {set:e} s to tell them apart"
    );

    // And the next step sheds what that phase's film sheds.
    let second = Time::from_si(0.25 * set);
    let (mut a, mut b) = (fed.clone(), fresh.clone());
    let (lost_a, lost_b) = (a.lost_energy().to_si(), b.lost_energy().to_si());
    a.step(Time::ZERO, second, &mut Exchange::new())
        .expect("a quarter of its limit");
    b.step(Time::ZERO, second, &mut Exchange::new())
        .expect("a quarter of its limit");
    let (shed_a, shed_b) = (
        a.lost_energy().to_si() - lost_a,
        b.lost_energy().to_si() - lost_b,
    );
    println!(
        "  the next step sheds {shed_a:.12e} J, and the block set to the phase {shed_b:.12e} J"
    );
    assert!(
        (shed_a / shed_b - 1.0).abs() < 1e-12,
        "the next step shed {shed_a:e} J, where the phase the bus heat left sheds {shed_b:e} J"
    );
}
