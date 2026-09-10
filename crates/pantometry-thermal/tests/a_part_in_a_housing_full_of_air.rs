//! A part inside a block's void, losing heat to the air around it.
//!
//! [`Solid3D::losing_from`] works on the block's six **outer** faces. A part rasterised inside a
//! block, or a bar with a clearance beside it, has surfaces interior to the grid — and nothing
//! could reach them, so such a part shed heat only by radiating to whatever faced it across the
//! gap. A housing full of air was a housing full of vacuum.
//!
//! It failed by producing nothing, which is the shape this workspace watches for: a scene could
//! state a film, run, audit, render and answer, with the film reaching no cell at all.
//! `a_bar_walled_in_by_void_cannot_shed_and_says_nothing` below is that state, kept as a test so
//! the fix has something to be measured against.
//!
//! # What is checked
//!
//! The lumped balance, which is exact here for the same reason it is in `a_block_that_generates`:
//! a copper bar 8 mm across has a Biot number of 8e-5 against a still-air film, so it is
//! isothermal and the algebra is the answer rather than an approximation to it. What is new is
//! that the conductance comes from the **grid** — `h · dx²` per face touching void — rather than
//! from a stated area, so the closed form counts faces.

use pantometry_core::{
    units::{Area, Length, Temperature, Time},
    Domain, Exchange, Substance,
};
use pantometry_thermal::{Environment, Face, Solid3D};

/// Copper with radiation switched off, so a balance is linear and a closed form is the problem's.
///
/// Radiation matters here for a reason worth keeping separate: the air term is convection only,
/// and what a surface exchanges with the surface across a gap is `find_gaps`. Turning emissivity
/// off leaves exactly the term under test.
fn grey() -> Substance {
    let mut s = Substance::copper();
    if let Some(t) = s.thermal.as_mut() {
        t.emissivity = 0.0;
    }
    s
}

/// A bar of `n` cells along x, one cell square, walled in by void on both sides along y.
///
/// The block is three cells deep in y and the bar is the middle one, so its two y faces are
/// **interior** — there is no outer face anywhere near them.
fn walled_in(n: usize, dx: f64, initial_c: f64) -> Solid3D {
    Solid3D::new(
        "bar",
        grey(),
        (n, 3, 1),
        Length::from_si(dx),
        Temperature::celsius(initial_c),
    )
    .empty(|_, j, _| j != 1)
}

/// **A bar walled in by void sheds nothing, and says nothing about it.**
///
/// The state this capability was added for, kept as a test rather than described. A film is stated
/// on `y-min`, which is a face made entirely of void: `cells_on` counts no solid cell there, so
/// nothing is charged and nothing is lost. The run completes.
///
/// This is not a defect to fix — a face of void really does cool nothing, and `pantometry-world`
/// refuses that spelling now. It is here because it is the thing `air_in` exists to make
/// expressible, and a capability whose absence is not written down is a capability nobody knows
/// they are missing.
#[test]
fn a_bar_walled_in_by_void_cannot_shed_and_says_nothing() {
    let mut bar = walled_in(8, 4e-3, 200.0).losing_from(
        Face::YMin,
        Environment {
            ambient: Temperature::celsius(20.0),
            convection_w_per_m2_k: 500.0,
            area: Area::from_si(1e-4),
        },
    );
    let start = bar.peak_temperature().to_si();
    let dt = Time::from_si(bar.max_stable_dt(Time::ZERO).to_si() * 0.5);
    for _ in 0..2000 {
        bar.step(Time::ZERO, dt, &mut Exchange::new())
            .expect("stable");
    }
    println!(
        "  a film on a face of void: {:.6} K after 2000 steps, {:.6e} J lost",
        bar.peak_temperature().to_si(),
        bar.lost_energy().to_si()
    );
    assert!(
        bar.peak_temperature().to_si() - start == 0.0,
        "something cooled it: {:.6} from {start:.6}",
        bar.peak_temperature().to_si()
    );
    assert!(bar.lost_energy().to_si() == 0.0, "it lost something");
    assert_eq!(bar.air(), None, "no air was asked for");
}

/// **Air in the void cools the surfaces facing it, at `h · dx²` a face.**
///
/// The closed form is the lumped balance with the conductance counted from the grid. The bar is
/// `n` cells long, one square, walled in along y — so every cell has two faces on void, and the
/// two end cells have a third each where the bar meets the block's x boundary... except that those
/// are the block's own outer faces and carry no void, so they do not. The count is asserted
/// separately, because a closed form built on the wrong number of faces agrees with a solver built
/// on the same wrong number.
#[test]
fn air_in_the_void_cools_what_faces_it() {
    let (n, dx, h) = (8usize, 4e-3, 500.0);
    let ambient = 20.0;
    let mut bar = walled_in(n, dx, 200.0).air_in(Temperature::celsius(ambient), h);

    assert_eq!(
        bar.air().map(|(t, h)| (t.to_si(), h)),
        Some((Temperature::celsius(ambient).to_si(), h))
    );
    // Two faces a cell, along y, and nothing else: the x ends and both z faces are the block's own
    // boundary, which is insulated because nothing was exposed.
    let faces = bar.faces_touching_void();
    println!("  {faces} faces of the bar touch void");
    assert_eq!(faces, 2 * n, "the bar has two void faces a cell");

    // The lumped balance. `Bi = h·(a/2)/k` is 8e-5 for 8 mm of copper against a 500 W/m²K film, so
    // the bar is isothermal and `ρcV dT/dt = −hA(T − T∞)` is exact rather than approximate.
    let g = h * (dx * dx) * faces as f64;
    let capacity = bar.heat_capacity().to_si();
    let tau = capacity / g;
    println!("  conductance {g:.6} W/K, capacity {capacity:.4} J/K, tau {tau:.3} s");

    let start = bar.mean_temperature().to_si();
    let air = Temperature::celsius(ambient).to_si();
    let mut elapsed = 0.0;
    let dt = Time::from_si(bar.max_stable_dt(Time::ZERO).to_si() * 0.5);
    while elapsed < tau {
        bar.step(Time::ZERO, dt, &mut Exchange::new())
            .expect("stable");
        elapsed += dt.to_si();
    }
    let want = air + (start - air) * (-elapsed / tau).exp();
    let got = bar.mean_temperature().to_si();
    println!("  after {elapsed:.3} s: {got:.4} K against the exponential's {want:.4} K");
    // The tolerance is the explicit scheme's own first order in time over one time constant, which
    // at `dt = τ/N` leaves about `dt/(2τ)` of the decay.
    let bound = dt.to_si() / (2.0 * tau) * (start - air);
    assert!(
        (got - want).abs() < bound.max(1e-9),
        "{got:.6} against {want:.6}, off {:.4} K where the step allows {bound:.4}",
        (got - want).abs()
    );

    // And the books balance: what left the cells arrived in the counter.
    let held = capacity * (start - got);
    println!(
        "  it gave up {:.6} J and the ledger says {:.6}",
        held,
        bar.lost_energy().to_si()
    );
    assert!(
        (bar.lost_energy().to_si() / held - 1.0).abs() < 1e-9,
        "the ledger says {:.9} J and the cells gave up {held:.9}",
        bar.lost_energy().to_si()
    );
}

/// **A cavity's air and an outer film add up**, and neither is charged twice.
///
/// The two paths are separate stores and separate terms, so the risk is arithmetic rather than
/// physics: a cell on the block's boundary that also touches void must carry both conductances,
/// once each. Checked at steady state, where the two in parallel have a closed form the single
/// terms do not.
#[test]
fn air_and_an_outer_film_are_two_conductances_in_parallel() {
    let (n, dx, h_air, h_film) = (8usize, 4e-3, 500.0, 300.0);
    let watts = 6.0;
    let ambient = 20.0;
    let air = Temperature::celsius(ambient);

    // The same bar, generating, with the outer film on `z-min` — which is the bar's own bottom, a
    // real solid face of the block — and air on the two y faces.
    let build = |with_air: bool, with_film: bool| {
        let mut b = walled_in(n, dx, ambient).dissipating(watts, |_, _, _| true);
        if with_film {
            b = b.losing_from(
                Face::ZMin,
                Environment {
                    ambient: air,
                    convection_w_per_m2_k: h_film,
                    area: Area::from_si(n as f64 * dx * dx),
                },
            );
        }
        if with_air {
            b = b.air_in(air, h_air);
        }
        b
    };

    let settle = |mut b: Solid3D| {
        let mut bus = Exchange::new();
        let mut previous = b.mean_temperature().to_si();
        for round in 0..8000 {
            let mut window = 0.0;
            while window < 1.0 {
                let dt = b.max_stable_dt(Time::ZERO);
                b.step(Time::ZERO, dt, &mut bus).expect("a stable step");
                window += dt.to_si();
            }
            let now = b.mean_temperature().to_si();
            if round > 0 && (now - previous).abs() / window < 1e-9 {
                return now - Temperature::celsius(ambient).to_si();
            }
            previous = now;
        }
        panic!("it never settled");
    };

    let both = settle(build(true, true));
    // The conductances, counted the way each is defined: the air from the grid, the film from the
    // area the caller stated. The half cell of copper behind the film is 2kA/dx, which at 8 mm of
    // copper is four orders above either, so it moves the answer by 0.02% and is left in.
    let g_air = h_air * dx * dx * (2 * n) as f64;
    let area = n as f64 * dx * dx;
    let half = 2.0
        * grey()
            .thermal
            .expect("copper conducts")
            .conductivity
            .to_si()
        * area
        / dx;
    let g_film = 1.0 / (1.0 / (h_film * area) + 1.0 / half);
    let closed = watts / (g_air + g_film);
    println!("  both paths: {both:.5} K against {closed:.5} K in parallel");
    assert!(
        (both / closed - 1.0).abs() < 2e-3,
        "{both:.5} K against a parallel {closed:.5} K"
    );

    // And each alone is the other's complement, which is what "in parallel" means and what a term
    // counted twice would break.
    let air_only = settle(build(true, false));
    let film_only = settle(build(false, true));
    println!("  air alone {air_only:.5} K, film alone {film_only:.5} K");
    assert!(
        (air_only / (watts / g_air) - 1.0).abs() < 2e-3,
        "air alone {air_only:.5} against {:.5}",
        watts / g_air
    );
    assert!(
        (film_only / (watts / g_film) - 1.0).abs() < 2e-3,
        "film alone {film_only:.5} against {:.5}",
        watts / g_film
    );
    assert!(
        both < air_only && both < film_only,
        "two paths must cool better than either: {both:.5} against {air_only:.5}/{film_only:.5}"
    );
}

/// **Air tightens the stability limit, and the limit it reports is one the march survives.**
///
/// The direction a sabotage does not find first. Leaving the air out of the limit makes the limit
/// *longer* than the flux needs — and an explicit boundary past its limit does not diverge loudly:
/// it oscillates about the air it is losing to while the conservation audit stays perfectly happy,
/// because the heat really did leave. `air_conductance` is read by both `film_flux` and
/// `loss_conductance_at`, and this is what says so.
#[test]
fn air_tightens_the_limit_and_the_march_survives_it() {
    let (n, dx) = (8usize, 4e-3);
    let bare = walled_in(n, dx, 400.0);
    let mut aired = walled_in(n, dx, 400.0).air_in(Temperature::celsius(20.0), 8000.0);
    let (loose, tight) = (
        bare.max_stable_dt(Time::ZERO).to_si(),
        aired.max_stable_dt(Time::ZERO).to_si(),
    );
    println!("  limit without air {loose:.6e} s, with it {tight:.6e} s");
    assert!(
        tight < loose,
        "air must tighten the limit: {tight:.3e} against {loose:.3e}"
    );

    // At its own limit, from uniform, every cell falls monotonically. An over-long limit goes
    // down, past the air, and back up.
    let mut previous = aired.peak_temperature().to_si();
    for step in 0..800 {
        let dt = aired.max_stable_dt(Time::ZERO);
        aired
            .step(Time::ZERO, dt, &mut Exchange::new())
            .expect("its own limit is a step it can take");
        let now = aired.peak_temperature().to_si();
        assert!(
            now <= previous + 1e-9,
            "it rose at step {step}: {previous:.6} then {now:.6} K"
        );
        previous = now;
    }
    println!("  800 steps at its own limit: down to {previous:.4} K");
    assert!(previous < 400.0 + 273.15, "nothing cooled: {previous:.4}");
}

/// **A core in a cavity sheds from cells the block's boundary never touches, and the books say so.**
///
/// The ledger walked only the block's own outer shell, on the reasoning that nothing else can lose
/// heat — true while a film is the only boundary flux, and false the moment a cavity has air in it.
/// A surface is 6% of the volume at 96³ and that shortcut is worth keeping, so it is kept *and*
/// switched off when there is air.
///
/// Every cell of this core is strictly interior: a 5×5×5 block whose only solid is the 3×3×3 at
/// its centre, so no core cell has an index of 0 or of `n − 1` on any axis. With the shortcut left
/// in, the cells shed and the counter never hears about it — the audit's books grow by exactly what
/// the air took, which reads as energy created from nothing.
#[test]
fn a_core_in_a_cavity_sheds_from_the_inside_and_the_ledger_hears_it() {
    let dx = 4e-3;
    let mut core = Solid3D::new(
        "core",
        grey(),
        (5, 5, 5),
        Length::from_si(dx),
        Temperature::celsius(300.0),
    )
    .empty(|i, j, k| ![i, j, k].iter().all(|a| (1..4).contains(a)))
    .air_in(Temperature::celsius(20.0), 400.0);

    // Every solid cell is off the block's boundary, so the old walk would have visited none of
    // them; and the core's own surface is 6 faces of 3x3 = 54.
    assert_eq!(core.void_cells(), 125 - 27, "only the 3x3x3 core is solid");
    assert_eq!(core.faces_touching_void(), 54, "the core has six 3x3 faces");

    let capacity = core.heat_capacity().to_si();
    let start = core.mean_temperature().to_si();
    let dt = Time::from_si(core.max_stable_dt(Time::ZERO).to_si() * 0.5);
    for _ in 0..600 {
        core.step(Time::ZERO, dt, &mut Exchange::new())
            .expect("stable");
    }
    let end = core.mean_temperature().to_si();
    let gave = capacity * (start - end);
    println!(
        "  the core fell {:.4} K, giving up {gave:.6} J; the ledger says {:.6}",
        start - end,
        core.lost_energy().to_si()
    );
    assert!(gave > 1.0, "it barely cooled: {gave:.6} J");
    assert!(
        (core.lost_energy().to_si() / gave - 1.0).abs() < 1e-9,
        "the ledger says {:.9} J and the cells gave up {gave:.9}",
        core.lost_energy().to_si()
    );
}
