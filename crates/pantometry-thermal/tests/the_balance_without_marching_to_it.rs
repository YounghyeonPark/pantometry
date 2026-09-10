//! `Solid3D::steady_state` against the closed forms, and against the march it replaces.
//!
//! A design answer is a steady-state answer — a junction temperature, a margin, a rise above
//! ambient — and the only way to get one was to run for several time constants and check that it
//! had stopped moving. `29-a-designed-bracket-becomes-cells` needed **900 s** from a cold start,
//! nine test binaries walk every shipped scene, and the app gate stopped fitting inside ten
//! minutes for it.
//!
//! # What agreement with the march is worth, and what it is not
//!
//! Marching is the **reference**: it is the thing already checked against series resistances,
//! exponential decay and measured convergence rates. So the tests below compare against a formula
//! computed here, or against a march that has been run long enough to have arrived — never against
//! the solve alone. That is the same rule `steady_state.rs` states for `ThermalNetwork`.
//!
//! The two are not independent implementations either: the residual is `flux_at`, the function the
//! sweep marches with. They agreeing says the fixed point was found, not that two people wrote the
//! same physics twice — which is why every test here also has a closed form.

use pantometry_core::{
    units::{Area, Length, Temperature, Time},
    Domain, Exchange, Substance,
};
use pantometry_thermal::{Environment, Face, Solid3D};

/// Aluminium with radiation switched off, so a balance is linear and its closed form is the
/// problem's rather than a linearisation's.
fn grey() -> Substance {
    let mut s = Substance::aluminium_6061();
    if let Some(t) = s.thermal.as_mut() {
        t.emissivity = 0.0;
    }
    s
}

fn film(h: f64, area: f64) -> Environment {
    Environment {
        ambient: Temperature::celsius(20.0),
        convection_w_per_m2_k: h,
        area: Area::from_si(area),
    }
}

/// March until the mean stops moving, checking that it arrived rather than marching a while and
/// calling it settled. Returns the seconds it took.
///
/// **The bound is generous because the thing being replaced is slow.** A body whose only path out
/// is a clearance has a radiative conductance of a third of a milliwatt per kelvin, so its time
/// constant is half an hour of simulated time and the march needs eight of those to stop moving.
/// That is the case `steady_state` exists for, and a helper that gave up at twenty thousand
/// seconds reported "the march never settled" for a scene that settles perfectly well.
fn settle_by_marching(block: &mut Solid3D) -> f64 {
    let mut bus = Exchange::new();
    let mut previous = block.mean_temperature().to_si();
    let mut elapsed = 0.0;
    for round in 0..60_000 {
        let mut window = 0.0;
        while window < 1.0 {
            let dt = block.max_stable_dt(Time::ZERO);
            block.step(Time::ZERO, dt, &mut bus).expect("a stable step");
            window += dt.to_si();
        }
        elapsed += window;
        let now = block.mean_temperature().to_si();
        if round > 0 && (now - previous).abs() / window < 1e-8 {
            return elapsed;
        }
        previous = now;
    }
    panic!("the march never settled in {elapsed:.0} s");
}

/// **A bar with a source at one end and a film at the other, solved rather than marched.**
///
/// The closed form is the one `a_block_that_generates` uses: every watt crosses every face and
/// leaves through the film, so the rise from the source cell's centre is
/// `P·((n − ½)·dx/(k·A) + 1/(h·A))`, the half cell being the solid the film sits behind.
///
/// Both routes are checked against it, and then against each other — which is the claim that
/// matters, because the solve is meant to replace the march and not to sit beside it.
#[test]
fn a_driven_bar_is_solved_to_the_same_place_it_marches_to() {
    let (watts, n, dx, h) = (5.0, 16usize, 3e-3, 4000.0);
    let area = dx * dx;
    let k = grey()
        .thermal
        .expect("aluminium conducts")
        .conductivity
        .to_si();
    let build = || {
        Solid3D::new(
            "bar",
            grey(),
            (n, 1, 1),
            Length::from_si(dx),
            Temperature::celsius(20.0),
        )
        .dissipating(watts, |i, _, _| i == 0)
        .losing_from(Face::XMax, film(h, area))
    };

    let closed = watts * ((n as f64 - 0.5) * dx / (k * area) + 1.0 / (h * area));

    let mut solved = build();
    let moved = solved.settle().expect("it has a steady state");
    let by_solve = solved.temperature_at(0, 0, 0).to_si() - Temperature::celsius(20.0).to_si();

    let mut marched = build();
    let marched_for = settle_by_marching(&mut marched);
    let by_march = marched.temperature_at(0, 0, 0).to_si() - Temperature::celsius(20.0).to_si();

    println!("  solved {by_solve:.6} K, marched {by_march:.6} K, closed form {closed:.6} K");
    println!(
        "  the solve moved the field by {moved:.4} K from cold; the march took {marched_for:.0} s"
    );
    assert!(
        (by_solve / closed - 1.0).abs() < 1e-6,
        "the solve says {by_solve:.6} K where the series says {closed:.6}"
    );
    assert!(
        (by_march / closed - 1.0).abs() < 2e-3,
        "the march says {by_march:.6} K where the series says {closed:.6}"
    );
    // And the two agree far more tightly than either's own tolerance, because they are the same
    // operator at the same fixed point.
    assert!(
        (by_solve - by_march).abs() < 1e-3,
        "solved {by_solve:.6} against marched {by_march:.6}"
    );

    // Settling a block that is already settled moves nothing.
    let again = solved.settle().expect("still has one");
    println!("  settling it again moved {again:.3e} K");
    assert!(again < 1e-5, "it was not at its own balance: {again:.3e} K");
}

/// **A block with radiation is solved too**, and the non-linearity is where a linear solver would
/// quietly stop being right.
///
/// The film's radiative half is `εσA(T⁴ − T∞⁴)`, so the balance is the root of a quartic rather
/// than a division. Solved here by bisection from the constants — the same `balance` function
/// `a_block_that_generates` uses — and the block is one cell so the lumped form is the exact
/// answer rather than an approximation to it.
#[test]
fn a_radiating_block_is_solved_to_the_root_of_its_quartic() {
    const SIGMA: f64 = 5.670_374_419e-8;
    let (watts, dx, h) = (12.0, 20e-3, 15.0);
    let area = 6.0 * dx * dx;
    let emissivity = Substance::aluminium_6061()
        .thermal
        .expect("aluminium radiates")
        .emissivity;
    let ambient = Temperature::celsius(20.0).to_si();

    let mut block = Solid3D::new(
        "lump",
        Substance::aluminium_6061(),
        (1, 1, 1),
        Length::from_si(dx),
        Temperature::celsius(20.0),
    )
    .dissipating(watts, |_, _, _| true)
    .losing_from(Face::XMin, film(h, area));

    // The root of `hA(T − T∞) + εσA(T⁴ − T∞⁴) = P`, **in series with the half cell of solid
    // the film sits behind**, by bisection from the constants. That half cell is `2kA/dx`, four
    // orders above the film here, so it moves the answer by 0.04% — which is nothing to the
    // physics and everything to a check written at one part in a million. Leaving it out measured
    // 588.520 against 588.255.
    let k = Substance::aluminium_6061()
        .thermal
        .expect("aluminium conducts")
        .conductivity
        .to_si();
    let half = 2.0 * k * area / dx;
    let (mut lo, mut hi) = (ambient, ambient + 10_000.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        let gap = mid - ambient;
        let bare = h * area * gap + emissivity * SIGMA * area * (mid.powi(4) - ambient.powi(4));
        let shed = if gap > 0.0 {
            gap / (gap / bare + 1.0 / half)
        } else {
            0.0
        };
        if shed < watts {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let closed = 0.5 * (lo + hi);

    block.settle().expect("it has a steady state");
    let got = block.temperature_at(0, 0, 0).to_si();
    println!("  solved {got:.6} K against the quartic's root {closed:.6} K");
    assert!(
        (got - closed).abs() < 1e-6 * closed,
        "{got:.6} against {closed:.6}"
    );
    // A solver that had linearised the radiation would land somewhere else entirely: at this
    // temperature the radiative term is a third of the loss.
    let linear = ambient + watts / (h * area);
    println!("  a purely convective block would sit at {linear:.4} K");
    assert!(
        (got - linear).abs() > 20.0,
        "the radiation is doing nothing, so this proves nothing: {got:.4} against {linear:.4}"
    );
}

/// **Heat with nowhere to go is refused, and so is a block that melts.**
///
/// Returning the last iterate would be a plausible temperature for a balance that was never
/// struck. A block generating heat with every face insulated warms without limit and has no steady
/// state at all; a block that melts has one, and it is not a statement about temperature, because
/// a mushy cell's state is its enthalpy.
#[test]
fn a_block_with_no_balance_to_strike_is_refused() {
    let mut sealed = Solid3D::new(
        "sealed",
        grey(),
        (4, 4, 4),
        Length::mm(4.0),
        Temperature::celsius(20.0),
    )
    .dissipating(3.0, |_, _, _| true);
    let why = sealed.settle().expect_err("it warms without limit");
    println!("  sealed: {why:?}");
    assert!(
        format!("{why:?}").contains("warms without limit"),
        "the refusal should say why: {why:?}"
    );

    // The same block with a face exposed does have one.
    let mut open = Solid3D::new(
        "open",
        grey(),
        (4, 4, 4),
        Length::mm(4.0),
        Temperature::celsius(20.0),
    )
    .dissipating(3.0, |_, _, _| true)
    .losing_from(Face::ZMax, film(50.0, 2.56e-4));
    open.settle().expect("a face to lose through is enough");

    let mut icy = Solid3D::new(
        "ice",
        Substance::ice(),
        (4, 4, 4),
        Length::mm(4.0),
        Temperature::celsius(-10.0),
    )
    .losing_from(Face::ZMax, film(50.0, 2.56e-4));
    let why = icy
        .settle()
        .expect_err("a melting block has no steady field");
    println!("  melting: {why:?}");
    assert!(
        format!("{why:?}").contains("enthalpy"),
        "the refusal should say why: {why:?}"
    );
}

/// **A clearance's pair exchange is solved with everything else**, which is the term a
/// cell-by-cell relaxation is least likely to carry.
///
/// A pair is not a face: it couples two cells that are not neighbours, through `T⁴`. Checked
/// against the march, which is the reference here, and against the physics — the hot side must sit
/// above the cold one and both above ambient, and a solve that dropped the pairing would leave the
/// fed half with nowhere to send its heat and refuse.
#[test]
fn a_clearance_is_solved_along_with_the_faces() {
    let dx = 4e-3;
    // **Blackened, and a fifth of a watt.** A pair carries `σA(T₁⁴−T₂⁴)/(1/ε+1/ε−1)`, and bare
    // aluminium at ε = 0.09 gives a coefficient of 4.3e-14 W/K⁴ across these two cells — so 0.4 W
    // would need the fed half at 1750 K, and the first version of this test asked exactly that.
    // The solve reported diverging, correctly: it is a scene whose answer is a melted bar.
    let mut black = Substance::aluminium_6061();
    if let Some(th) = black.thermal.as_mut() {
        th.emissivity = 0.9;
    }
    let build = || {
        Solid3D::new(
            "pair",
            black.clone(),
            (2, 1, 5),
            Length::from_si(dx),
            Temperature::celsius(20.0),
        )
        .empty(|_, _, k| k == 2)
        .dissipating(0.02, |_, _, k| k < 2)
        .losing_from(Face::ZMax, film(30.0, 2.0 * dx * dx))
    };
    let mut solved = build();
    solved.settle().expect("the pair carries it");
    let (hot, cold) = (
        solved.temperature_at(0, 0, 1).to_si(),
        solved.temperature_at(0, 0, 3).to_si(),
    );

    let mut marched = build();
    let marched_for = settle_by_marching(&mut marched);
    let (mhot, mcold) = (
        marched.temperature_at(0, 0, 1).to_si(),
        marched.temperature_at(0, 0, 3).to_si(),
    );
    println!(
        "  solved {hot:.4}/{cold:.4} K, marched {mhot:.4}/{mcold:.4} K after {marched_for:.0} s"
    );
    assert!(
        (hot - mhot).abs() < 0.05 && (cold - mcold).abs() < 0.05,
        "the solve and the march disagree: {hot:.4}/{cold:.4} against {mhot:.4}/{mcold:.4}"
    );
    // The fed half is hotter and both are above the air, which is the only arrangement the pair
    // and the film together allow.
    assert!(
        hot > cold && cold > Temperature::celsius(20.0).to_si() + 1.0,
        "heat flowed the wrong way: {hot:.4} and {cold:.4}"
    );
}

/// **A long bar needs the over-relaxation, and a short one does not say so.**
///
/// Relaxation's convergence is set by its spectral radius, and for plain Gauss-Seidel on an
/// `n`-cell chain that is `cos²(π/2n)` — so the sweeps needed grow as `n²`. At sixteen cells that
/// is a few thousand and fits inside the bound; at sixty-four it is **thirty-four thousand** and
/// does not. Over-relaxing by `ω = 2/(1 + sin(π/n))` makes it grow as `n` instead.
///
/// Every other test here is sixteen cells or fewer, and all four pass with the over-relaxation
/// removed — which is what this exists to fix. The scenes that matter are bigger: the bracket is
/// twenty-six cells along its longest axis and eight thousand in total.
#[test]
fn a_long_bar_needs_the_over_relaxation() {
    let (watts, n, dx, h) = (5.0, 64usize, 1e-3, 4000.0);
    let area = dx * dx;
    let k = grey()
        .thermal
        .expect("aluminium conducts")
        .conductivity
        .to_si();
    let mut bar = Solid3D::new(
        "long",
        grey(),
        (n, 1, 1),
        Length::from_si(dx),
        Temperature::celsius(20.0),
    )
    .dissipating(watts, |i, _, _| i == 0)
    .losing_from(Face::XMax, film(h, area));

    bar.settle()
        .expect("sixty-four cells is not a hard problem");
    let rise = bar.temperature_at(0, 0, 0).to_si() - Temperature::celsius(20.0).to_si();
    let closed = watts * ((n as f64 - 0.5) * dx / (k * area) + 1.0 / (h * area));
    println!("  {n} cells: solved {rise:.6} K against the series' {closed:.6} K");
    assert!(
        (rise / closed - 1.0).abs() < 1e-5,
        "{rise:.6} K against {closed:.6}"
    );

    // And the gradient is the linear one a uniformly conducting bar has, cell by cell — a solve
    // that had stopped short would leave the far end lagging.
    let step =
        |i: usize| bar.temperature_at(i, 0, 0).to_si() - bar.temperature_at(i + 1, 0, 0).to_si();
    let (first, middle, last) = (step(0), step(n / 2), step(n - 2));
    println!("  steps along it: {first:.6}, {middle:.6}, {last:.6} K");
    assert!(
        (first / middle - 1.0).abs() < 1e-6 && (last / middle - 1.0).abs() < 1e-6,
        "the gradient is not uniform: {first:.6}, {middle:.6}, {last:.6}"
    );
}
