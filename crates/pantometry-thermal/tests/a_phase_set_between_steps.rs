//! A phase named between steps is the phase the next step conducts through.
//!
//! For a substance with a liquid phase, a cell's conductivity and capacity follow how much of it
//! has melted, so anything that moves a fraction moves the operator: the face conductances, the
//! capacities and the stability limit. A step re-derives them from the fractions it leaves. The
//! three mutators that move a fraction **between** steps — `set_melted_fraction`,
//! `set_temperature` and `deposit` — did not, and the operator the block had before them stayed
//! until the end of the next step.
//!
//! The worst of it was the limit. Water at 5 °C declared ice at its melting point reported 3.206 s
//! for a column whose own limit was 0.423 s, and one step at the limit it reported took the cell
//! set to −20 °C up to 0 °C and its neighbour down to −20.4 °C: past each other, out of the range
//! the column started in, with energy conserved to the bit, so no audit saw it.
//!
//! `fill` had a neighbouring defect: a cell given a substance that melts kept the fraction of what
//! was there before. A wax melting at 28.15 °C filled into aluminium at 40 °C was a solid at 40 °C,
//! its books opened 0.199 J short — its latent heat — and with nothing flowing it came out of one
//! step at 28.15 °C and 9.4% melted.
//!
//! # What is checked against what
//!
//! A mutator now re-derives the cell it moved and nothing else, so the first test holds that to the
//! bit against a full re-derivation of the same state, reached two ways: every face, every capacity
//! and mobility, the limit, and one step on from each. That is an exactness claim about an
//! incremental update, which is the one place comparing two computations is the check rather than
//! a substitute for one. The second holds the deferred limit to a refusal that was already there. The
//! third and fifth are closed forms — the limit of a uniform column, the state one step at that
//! limit leaves, the phase a temperature names — and the fourth is the restore contract.

use pantometry_core::conserved::quantity;
use pantometry_core::substance::{FusionProps, ThermalProps};
use pantometry_core::units::{
    Area, Density, Energy, LatentHeat, Length, SpecificHeat, Temperature, ThermalConductivity,
    ThermalExpansion, Time,
};
use pantometry_core::{Domain, Exchange, Substance};
use pantometry_thermal::{Axis, Environment, Face, Solid3D};

/// Ice with water named as its liquid, which is what makes a block two-phase.
fn two_phase_ice() -> Substance {
    let base = Substance::ice();
    let fusion = base.fusion.expect("ice melts");
    base.with_fusion(fusion.with_liquid(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(0.598),
        specific_heat: SpecificHeat::j_per_kg_k(4182.0),
        expansion: ThermalExpansion::ppm_per_k(207.0),
        emissivity: 0.96,
    }))
}

/// A paraffin melting at 28.15 °C, with `21-a-wax-thermal-buffer`'s solid and none of its liquid:
/// one-phase, so what is checked of it is the phase and not an operator.
fn wax() -> Substance {
    Substance::bulk("wax", Density::from_si(814.0))
        .with_thermal(ThermalProps {
            conductivity: ThermalConductivity::w_per_m_k(0.358),
            specific_heat: SpecificHeat::j_per_kg_k(1934.0),
            expansion: ThermalExpansion::ppm_per_k(800.0),
            emissivity: 0.9,
        })
        .with_fusion(FusionProps::new(
            Temperature::from_si(301.3),
            LatentHeat::from_si(244_000.0),
        ))
}

fn bits(x: f64) -> u64 {
    x.to_bits()
}

/// The same block, re-derived from scratch, by two roads: a `fill` with nothing, which re-runs the
/// whole derivation without changing a cell, and a checkpoint and a restore, which re-derives a
/// two-phase block's operator from the fractions it puts back. Two, so a shortcut taken by either —
/// a `fill` that returns early when it matched nothing — does not quietly make the reference a copy
/// of the block under test.
fn rederived(block: &Solid3D) -> [(&'static str, Solid3D); 2] {
    let mut filled = block.clone();
    let substance = filled.substance_at(0, 0, 0).clone();
    filled.fill(substance, |_, _, _| false);
    let mut restored = block.clone();
    restored.checkpoint();
    restored.restore();
    [
        ("a fill of nothing", filled),
        ("a checkpoint and restore", restored),
    ]
}

/// Where two blocks' operators differ: every face, capacity and mobility, and the limit — as
/// `max_stable_dt` gives it, and as the stability ratio of a second, which keeps a `NaN` the limit
/// reports as an infinity.
fn operator_differences(got: &Solid3D, want: &Solid3D, against: &str) -> Vec<String> {
    let (g, w) = (got.coefficients(), want.coefficients());
    let mut out = Vec::new();
    for (name, a, b) in [
        ("kx", g.kx, w.kx),
        ("ky", g.ky, w.ky),
        ("kz", g.kz, w.kz),
        ("mobility", g.mobility, w.mobility),
        ("capacity", got.cell_capacities(), want.cell_capacities()),
    ] {
        if let Some(at) = (0..a.len()).find(|&f| bits(a[f]) != bits(b[f])) {
            out.push(format!(
                "{name}[{at}] {:e} where {against} has {:e}",
                a[at], b[at]
            ));
        }
    }
    let (lg, lw) = (
        got.max_stable_dt(Time::from_si(0.0)).to_si(),
        want.max_stable_dt(Time::from_si(0.0)).to_si(),
    );
    if bits(lg) != bits(lw) {
        out.push(format!("limit {lg} s where {against} has {lw} s"));
    }
    let (rg, rw) = (
        got.stability_ratio(Time::from_si(1.0)),
        want.stability_ratio(Time::from_si(1.0)),
    );
    if bits(rg) != bits(rw) {
        out.push(format!(
            "stability ratio {rg} at 1 s where {against} has {rw}"
        ));
    }
    out
}

/// Where two blocks' states differ, cell by cell, and in their books.
fn state_differences(got: &Solid3D, want: &Solid3D) -> Vec<String> {
    let (nx, ny, nz) = got.counts();
    let mut out = Vec::new();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let (tg, tw) = (got.temperature_at(i, j, k), want.temperature_at(i, j, k));
                let (mg, mw) = (
                    got.melted_fraction_at(i, j, k),
                    want.melted_fraction_at(i, j, k),
                );
                if bits(tg.to_si()) != bits(tw.to_si()) || bits(mg) != bits(mw) {
                    out.push(format!(
                        "cell ({i},{j},{k}) {:.6} K, {mg} melted where {:.6} K, {mw}",
                        tg.to_si(),
                        tw.to_si()
                    ));
                }
            }
        }
    }
    let (eg, ew) = (
        got.ledger().get(quantity::ENERGY),
        want.ledger().get(quantity::ENERGY),
    );
    if eg.map(bits) != ew.map(bits) {
        out.push(format!("books {eg:?} J where {ew:?} J"));
    }
    out
}

/// A block of two-phase ice with everything the operator is built from in it: three axes, a cell
/// of void, contacts of a different conductance on each axis, and a film on one face.
fn awkward(initial_c: f64) -> Solid3D {
    Solid3D::new(
        "awkward",
        two_phase_ice(),
        (4, 3, 2),
        Length::mm(1.0),
        Temperature::celsius(initial_c),
    )
    .empty(|i, j, k| (i, j, k) == (3, 0, 1))
    .joined(|axis, i, j, _| match axis {
        Axis::X if i == 2 => Some(4.0e4),
        Axis::Y if j == 1 => Some(9.0e4),
        Axis::Z if i == 1 => Some(2.5e4),
        _ => None,
    })
    .losing_from(
        Face::YMax,
        Environment::still_air(Temperature::celsius(-20.0), Area::cm2(0.12)),
    )
}

/// One block before its mutators and after them, and whether the mutators have to move the limit
/// for the scenario to mean anything.
struct Scenario {
    name: &'static str,
    before: Solid3D,
    after: Solid3D,
    moves_the_limit: bool,
}

fn scenario(
    name: &'static str,
    before: Solid3D,
    moves_the_limit: bool,
    mutate: impl FnOnce(&mut Solid3D),
) -> Scenario {
    let mut after = before.clone();
    mutate(&mut after);
    Scenario {
        name,
        before,
        after,
        moves_the_limit,
    }
}

/// Every mutator that moves a phase between steps leaves the operator a full re-derivation of the
/// same state would, and the next step goes where that one's does.
///
/// Compared after each scenario's mutators, not after each call: in the one that declares every
/// cell, a later neighbour's re-derivation rewrites most of what an earlier one wrote, which is why
/// the isolated scenarios are there as well. Each scenario has to move the operator, or it would
/// pass with nothing to compare.
#[test]
fn every_mutator_leaves_the_operator_a_full_resolve_would() {
    let latent_per_cell = 917.0 * 1e-9 * 333_550.0;
    // Enough to take a cell of ice at −5 °C to water at 5 °C.
    let melt = 917.0 * 1e-9 * (2050.0 * 5.0 + 4182.0 * 5.0) + latent_per_cell;
    let mut scenarios = vec![
        scenario("set_melted_fraction", awkward(5.0), true, |b| {
            b.set_melted_fraction(0, 0, 0, 0.0);
            b.set_melted_fraction(1, 1, 0, 0.5);
            b.set_melted_fraction(3, 2, 1, 0.0);
            b.set_melted_fraction(2, 0, 1, 1.0); // already liquid: moves nothing
            b.set_melted_fraction(1, 1, 0, 0.25); // the same cell again
        }),
        scenario("set_temperature", awkward(0.0), false, |b| {
            b.set_temperature(0, 0, 0, Temperature::celsius(5.0));
            b.set_temperature(2, 1, 1, Temperature::celsius(3.0));
            b.set_temperature(3, 2, 0, Temperature::celsius(-5.0)); // stays solid
            b.set_temperature(2, 1, 1, Temperature::celsius(-1.0)); // and back
            b.set_temperature(1, 2, 1, Temperature::celsius(2.0));
        }),
        // Half as much to another, and some taken back out of the first.
        scenario("deposit", awkward(-5.0), false, |b| {
            b.deposit(1, 1, 1, Energy::from_si(melt));
            b.deposit(2, 2, 0, Energy::from_si(0.5 * melt));
            b.deposit(1, 1, 1, Energy::from_si(-0.6 * latent_per_cell));
        }),
        // Every cell, one at a time and in index order, so each cell's neighbours move after it.
        scenario("every cell declared solid", awkward(5.0), true, |b| {
            let (nx, ny, nz) = b.counts();
            for k in 0..nz {
                for j in 0..ny {
                    for i in 0..nx {
                        b.set_melted_fraction(i, j, k, 0.0);
                    }
                }
            }
            b.set_temperature(0, 0, 0, Temperature::celsius(-20.0));
        }),
    ];

    // Two substances, so a cell's pair has to be its own substance's and a face can be between
    // two: aluminium, and ice in every cell past the first plane, joined to it on that face.
    let mut two = Solid3D::new(
        "two substances",
        Substance::aluminium_6061(),
        (4, 3, 2),
        Length::mm(1.0),
        Temperature::celsius(-5.0),
    )
    .joined(|axis, i, _, _| (axis == Axis::X && i == 1).then_some(3.0e4));
    two.fill(two_phase_ice(), |i, _, _| i >= 1);
    scenarios.push(scenario("two substances", two, false, |b| {
        b.set_temperature(1, 1, 0, Temperature::celsius(5.0)); // on the face with aluminium
        b.set_melted_fraction(2, 2, 1, 0.5);
        b.deposit(1, 0, 1, Energy::from_si(melt));
    }));

    // A cell's neighbours' conduction sums are seen only through the limit, and only where one of
    // them sets it. In a 3×3×3 block of ice the one cell with six neighbours does, and melting any
    // cell beside it lowers its sum without moving it off the maximum — so each of the six puts
    // one direction of the update in the limit, and is required to move it.
    for (name, at) in [
        ("melted below the centre in x", (0, 1, 1)),
        ("melted above the centre in x", (2, 1, 1)),
        ("melted below the centre in y", (1, 0, 1)),
        ("melted above the centre in y", (1, 2, 1)),
        ("melted below the centre in z", (1, 1, 0)),
        ("melted above the centre in z", (1, 1, 2)),
    ] {
        let cube = Solid3D::new(
            "cube",
            two_phase_ice(),
            (3, 3, 3),
            Length::mm(1.0),
            Temperature::celsius(-5.0),
        );
        scenarios.push(scenario(name, cube, true, |b| {
            b.set_temperature(at.0, at.1, at.2, Temperature::celsius(5.0));
        }));
    }

    let mut failures = Vec::new();
    for s in scenarios {
        let moved = operator_differences(&s.after, &s.before, "the block before");
        if moved.is_empty() {
            failures.push(format!("{}: moved nothing this test can see", s.name));
        }
        if s.moves_the_limit && !moved.iter().any(|d| d.starts_with("limit ")) {
            failures.push(format!(
                "{}: did not move the limit, so a stale sum would not show",
                s.name
            ));
        }
        let [by_fill, by_restore] = rederived(&s.after);
        for (against, full) in [&by_fill, &by_restore] {
            for d in operator_differences(&s.after, full, against) {
                failures.push(format!("{}: {d}", s.name));
            }
        }
        let dt = Time::from_si(0.9 * by_fill.1.max_stable_dt(Time::from_si(0.0)).to_si());
        let (mut a, mut b) = (s.after, by_fill.1);
        b.step(Time::from_si(0.0), dt, &mut Exchange::new())
            .expect("the full re-derivation is stable at 0.9 of its own limit");
        // A refusal is a difference too, and the one a wrong operator tighter than the right one
        // gives: recorded with the rest rather than ending the test before they are said.
        match a.step(Time::from_si(0.0), dt, &mut Exchange::new()) {
            Err(v) => failures.push(format!(
                "{}: refused a step its full re-derivation takes, {} {} against {}",
                s.name, v.quantity, v.after, v.before
            )),
            Ok(()) => {
                for d in state_differences(&a, &b) {
                    failures.push(format!("{}, one step on: {d}", s.name));
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "a mutator left an operator a full resolve would not have:\n  {}",
        failures.join("\n  ")
    );
}

/// **A block that cannot be stepped still cannot after a phase moves.**
///
/// A substance that does not say what it conducts, and a contact that is not a number, both leave
/// a cell's rate `NaN`, and the block is refused rather than stepped. Moving a phase defers the
/// limit to a search over every cell, and that search has to keep the `NaN` the derivation kept:
/// `f64::max` passes over one, so a search made of it would find the finite rates of the other
/// cells, report a limit, and let the step write `NaN` into the cells that have none.
#[test]
fn a_block_that_cannot_be_stepped_still_cannot_after_a_phase_moves() {
    let mystery = Substance::bulk("mystery", Density::from_si(1000.0));
    let mut unknown = Solid3D::new(
        "unknown",
        two_phase_ice(),
        (3, 1, 1),
        Length::mm(1.0),
        Temperature::celsius(-5.0),
    );
    unknown.fill(mystery, |i, _, _| i == 0);
    let bad_joint = Solid3D::new(
        "bad joint",
        two_phase_ice(),
        (3, 1, 1),
        Length::mm(1.0),
        Temperature::celsius(-5.0),
    )
    .joined(|axis, i, _, _| (axis == Axis::X && i == 1).then_some(f64::NAN));

    let mut wrong = Vec::new();
    for (block, why) in [
        (unknown, "substance has no diffusivity"),
        (bad_joint, "a joint's contact conductance is not a number"),
    ] {
        for (when, mut b) in [("before", block.clone()), ("after a phase moved", block)] {
            if when != "before" {
                b.set_temperature(2, 0, 0, Temperature::celsius(5.0));
            }
            let name = b.name().to_string();
            let limit = b.max_stable_dt(Time::from_si(0.0)).to_si();
            if limit != f64::INFINITY {
                wrong.push(format!(
                    "{name}, {when}: a limit of {limit} s for a block with a cell whose rate is not a number"
                ));
            }
            match b.step(
                Time::from_si(0.0),
                Time::from_si(0.01),
                &mut Exchange::new(),
            ) {
                Err(v) if v.quantity == why => {}
                Err(v) => wrong.push(format!(
                    "{name}, {when}: refused for {:?}, not {why:?}",
                    v.quantity
                )),
                Ok(()) => wrong.push(format!(
                    "{name}, {when}: stepped, and cell 0 is now {:?}",
                    b.temperature_at(0, 0, 0)
                )),
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The column the limit was wrong for: water declared ice, one cell set to −20 °C.
fn declared_column() -> Solid3D {
    let mut c = Solid3D::new(
        "column",
        two_phase_ice(),
        (1, 1, 8),
        Length::mm(1.0),
        Temperature::celsius(5.0),
    );
    for k in 0..8 {
        c.set_melted_fraction(0, 0, k, 0.0);
    }
    c.set_temperature(0, 0, 0, Temperature::celsius(-20.0));
    c
}

/// **A liquid declared solid is stepped as the solid it is.**
///
/// Against two closed forms. The limit of a uniform column of ice is `dx²/(2α)` — Gershgorin's
/// bound on an interior row, which is what the block reports as its limit, with ice's own
/// `α = k/(ρc)`. One step at exactly that limit moves each face's two cells half way to each other,
/// so the cell at −20 °C and its neighbour at the melting point both land on −10 °C and nothing
/// melts. Only those two cells, the face between them and the limit carry information here; the
/// first test is the one that would see a single cell left as it was.
#[test]
fn a_liquid_declared_solid_is_stepped_as_the_solid_it_is() {
    let dx = 1e-3;
    let alpha = 2.22 / (917.0 * 2050.0);
    let want = dx * dx / (2.0 * alpha);
    let c = declared_column();
    let limit = c.max_stable_dt(Time::from_si(0.0)).to_si();
    // Collected, so a column wrong in both ways says both.
    let mut wrong = Vec::new();
    // Two expressions for one number, rounded in different orders: measured one ulp apart, and
    // `+ − × ÷` round the same everywhere, so four is room rather than hope.
    if (limit / want - 1.0).abs() >= 1e-15 {
        wrong.push(format!(
            "the limit is ice's dx²/(2α) = {want:.6} s, and the block reports {limit:.6} s"
        ));
    }

    let mut stepped = c.clone();
    match stepped.step(
        Time::from_si(0.0),
        Time::from_si(limit),
        &mut Exchange::new(),
    ) {
        Err(v) => wrong.push(format!("refused at its own limit: {v:?}")),
        Ok(()) => {
            for k in 0..8 {
                let t = stepped.temperature_at(0, 0, k).to_si() - 273.15;
                let expected = if k < 2 { -10.0 } else { 0.0 };
                // A step an ulp off the limit moves ten kelvin by 2e-15 K, and the enthalpy round
                // trip rounds at 273 K, whose ulp is 5.7e-14 K. Measured: exactly zero.
                if (t - expected).abs() >= 1e-12 {
                    wrong.push(format!(
                        "cell {k} should be at {expected} °C one step on, and is at {t:.6} °C"
                    ));
                }
                let melted = stepped.melted_fraction_at(0, 0, k);
                if melted != 0.0 {
                    wrong.push(format!(
                        "cell {k} is {melted} melted, in a column nothing warmer than its melting \
                         point touches"
                    ));
                }
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// **A restore after a declaration gives back the operator the declared fractions name.**
///
/// Ice declared water at its melting point, with one cell set to −20 °C: a step freezes part of
/// the cell beside the cold one, which moves its operator, so a restore to a checkpoint taken after
/// the declaration has to give back the one the declared fractions name — and the step after it
/// has to be the step after the checkpoint. As shipped the first step ran on the ice's operator,
/// which the constructor had derived at 0 °C, and the step after a restore on the water's: the
/// cold cell came out at −15.000 °C after one and −17.878 °C after the other.
#[test]
fn a_restore_after_a_declaration_gives_back_the_declared_operator() {
    let mut restored = Solid3D::new(
        "column",
        two_phase_ice(),
        (1, 1, 8),
        Length::mm(1.0),
        Temperature::celsius(0.0),
    );
    for k in 0..8 {
        restored.set_melted_fraction(0, 0, k, 1.0);
    }
    restored.set_temperature(0, 0, 0, Temperature::celsius(-20.0));
    restored.checkpoint();
    let dt = Time::from_si(0.5 * restored.max_stable_dt(Time::from_si(0.0)).to_si());
    restored
        .step(Time::from_si(0.0), dt, &mut Exchange::new())
        .expect("stable");
    let first = restored.clone();
    assert!(
        first.melted_fraction_at(0, 0, 1) < 1.0,
        "the step should freeze some of the cell beside the cold one, or nothing here moves the \
         operator and the restore below checks nothing"
    );
    restored.restore();
    restored
        .step(Time::from_si(0.0), dt, &mut Exchange::new())
        .expect("stable");
    let differences = state_differences(&restored, &first);
    assert!(
        differences.is_empty(),
        "a step after restoring to a checkpoint taken after the declaration is not the step after \
         the checkpoint:\n  {}",
        differences.join("\n  ")
    );
}

/// **A substance filled into a block takes the phase its temperature names**, and a cell whose
/// substance did not change keeps the phase it had.
#[test]
fn a_substance_filled_in_takes_the_phase_its_temperature_names() {
    // Into aluminium at 40 °C, above the wax's melting point: liquid, and nothing has arrived, so
    // the books are exactly zero and a step with nothing flowing moves nothing.
    let mut warm = Solid3D::new(
        "warm",
        Substance::aluminium_6061(),
        (3, 1, 1),
        Length::mm(1.0),
        Temperature::celsius(40.0),
    );
    warm.fill(wax(), |i, _, _| i == 1);
    // Collected, so a block wrong in all three ways says all three.
    let mut wrong = Vec::new();
    let phase = warm.melted_fraction_at(1, 0, 0);
    if phase != 1.0 {
        wrong.push(format!(
            "wax at 40 °C, twelve kelvin over its melting point, is liquid, and is {phase} melted"
        ));
    }
    let books = warm.ledger().get(quantity::ENERGY);
    if books != Some(0.0) {
        wrong.push(format!(
            "nothing has arrived, so the books open at zero, and they open at {books:?} J"
        ));
    }
    let dt = Time::from_si(0.5 * warm.max_stable_dt(Time::from_si(0.0)).to_si());
    warm.step(Time::from_si(0.0), dt, &mut Exchange::new())
        .expect("stable");
    let t = warm.temperature_at(1, 0, 0).to_si() - 273.15;
    let after = warm.melted_fraction_at(1, 0, 0);
    // Through the enthalpy and back, a multiply and a divide by one capacity and the cell's latent
    // heat added and taken away: an ulp or two of 313 K, and 2e-14 K more.
    if (t - 40.0).abs() >= 1e-12 || after != 1.0 {
        wrong.push(format!(
            "a uniform block moves nothing, and the wax went from 40 °C to {t:.6} °C, {after:.4} \
             melted"
        ));
    }
    assert!(wrong.is_empty(), "{}", wrong.join("; "));

    // The same composition built the other way round is the same block.
    let mut other = Solid3D::new(
        "warm",
        wax(),
        (3, 1, 1),
        Length::mm(1.0),
        Temperature::celsius(40.0),
    );
    other.fill(Substance::aluminium_6061(), |i, _, _| i != 1);
    let mut once = Solid3D::new(
        "warm",
        Substance::aluminium_6061(),
        (3, 1, 1),
        Length::mm(1.0),
        Temperature::celsius(40.0),
    );
    once.fill(wax(), |i, _, _| i == 1);
    let differences = state_differences(&once, &other);
    assert!(
        differences.is_empty(),
        "wax filled into aluminium is not aluminium filled into wax:\n  {}",
        differences.join("\n  ")
    );

    // Into a block below the melting point it is solid, which is the same rule and not a default.
    let mut cold = Solid3D::new(
        "cold",
        Substance::aluminium_6061(),
        (3, 1, 1),
        Length::mm(1.0),
        Temperature::celsius(-10.0),
    );
    cold.fill(wax(), |i, _, _| i == 1);
    assert_eq!(
        cold.melted_fraction_at(1, 0, 0),
        0.0,
        "wax at −10 °C is solid"
    );

    // A half-melted cell keeps its fraction through a fill somewhere else, and through a fill of
    // its own cell with the substance it already holds: nothing new is in it.
    let mut slush = Solid3D::new(
        "slush",
        two_phase_ice(),
        (4, 1, 1),
        Length::mm(1.0),
        Temperature::celsius(0.0),
    );
    slush.set_melted_fraction(1, 0, 0, 0.5);
    slush.fill(Substance::aluminium_6061(), |i, _, _| i == 3);
    assert_eq!(
        slush.melted_fraction_at(1, 0, 0),
        0.5,
        "a fill that did not touch a half-melted cell re-derived its phase anyway"
    );
    slush.fill(two_phase_ice(), |i, _, _| i == 1);
    assert_eq!(
        slush.melted_fraction_at(1, 0, 0),
        0.5,
        "a fill of a half-melted cell with the ice it already holds re-derived its phase"
    );

    // And a cell refilled from void is new matter: the fraction it held before is about nothing.
    let mut refilled = Solid3D::new(
        "refilled",
        two_phase_ice(),
        (4, 1, 1),
        Length::mm(1.0),
        Temperature::celsius(0.0),
    );
    refilled.set_melted_fraction(1, 0, 0, 0.5);
    let mut refilled = refilled.empty(|i, _, _| i == 1);
    refilled.fill(two_phase_ice(), |i, _, _| i == 1);
    assert_eq!(
        refilled.melted_fraction_at(1, 0, 0),
        0.0,
        "ice put into a void at its melting point is the unmelted phase, as a temperature names it"
    );
}
