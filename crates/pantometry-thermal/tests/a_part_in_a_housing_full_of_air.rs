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
//! The lumped balance, which is the scheme's own answer here rather than an approximation to it:
//! the bar is one cell across and every cell of it has the same faces on the air, so the block
//! holds it uniform whatever its Biot number — 2.5e-3 for copper 4 mm a cell under 500 W/m²K,
//! which is the half cell's share of the path and is in the closed form. What is new is that the
//! conductance comes from the **grid** — `h · dx²` per face touching void, in series with the half
//! cell behind it — rather than from a stated area, so the closed form counts faces.

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

/// One face's air, W/K: `h·dx²` in series with the half cell behind it, `2k·dx`, which is where
/// the air acts — at the surface, half a cell from the centre the block knows.
fn per_face(h: f64, k: f64, dx: f64) -> f64 {
    1.0 / (1.0 / (h * dx * dx) + 1.0 / (2.0 * k * dx))
}

/// Whether a residual is out: at or past `tolerance`, or not a number.
fn off_by(residual: f64, tolerance: f64) -> bool {
    residual.is_nan() || residual >= tolerance
}

/// Copper's conductivity, W/m·K.
fn copper_k() -> f64 {
    grey()
        .thermal
        .expect("copper conducts")
        .conductivity
        .to_si()
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

/// **Air in the void cools the surfaces facing it, at `h · dx²` a face behind its half cell.**
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

    // The lumped balance, `C dT/dt = −G(T − T∞)`, which the block holds the bar to exactly: it is
    // one cell across, so it stays uniform. Each face's air is behind the bar's half cell, `2k·dx`,
    // which is 3.2 W/K against the film's 8 mW/K and moves the conductance by 0.25% — the Biot
    // number, 2.5e-3 — small, and after a time constant 9.2e-4 of the drop, where the check on the
    // scheme's own answer below is met to 1.6e-15 of it.
    let g = per_face(h, copper_k(), dx) * faces as f64;
    let capacity = bar.heat_capacity().to_si();
    let tau = capacity / g;
    println!("  conductance {g:.6} W/K, capacity {capacity:.4} J/K, tau {tau:.3} s");

    let start = bar.mean_temperature().to_si();
    let air = Temperature::celsius(ambient).to_si();
    let mut elapsed = 0.0;
    let mut steps = 0usize;
    let dt = Time::from_si(bar.max_stable_dt(Time::ZERO).to_si() * 0.5);
    while elapsed < tau {
        bar.step(Time::ZERO, dt, &mut Exchange::new())
            .expect("stable");
        elapsed += dt.to_si();
        steps += 1;
    }
    // A time constant at half the limit is four hundred steps; none, and the check below would
    // compare the start with itself.
    assert!(
        steps >= 100,
        "only {steps} steps were taken over a time constant of {tau:.3} s"
    );
    let got = bar.mean_temperature().to_si();

    // **The scheme's own answer, exactly.** Every cell of the bar has the same two faces on the
    // air and its neighbours at its own temperature, so the bar stays uniform and a step
    // multiplies what it holds above the air by `1 − dt·G/C`: after `N` steps, that to the `N`th.
    // What is left is the rounding of each step's update, which a step makes in a different order
    // from this product's: four hundred steps at four ulps of 473 K each at most would be 5.1e-13
    // of the drop, and 1.6e-15 is measured, so 1e-12 is earned and is a billionth of what the air
    // charged against the centre would miss by.
    let factor = 1.0 - dt.to_si() * g / capacity;
    let mut excess = start - air;
    for _ in 0..steps {
        excess *= factor;
    }
    let discrete = air + excess;
    let rounding = (got - discrete).abs() / (start - air);
    println!(
        "  after {steps} steps: {got:.9} K against the scheme's own {discrete:.9} K, {rounding:.2e} \
         of the drop"
    );
    assert!(
        !off_by(rounding, 1e-12),
        "{got:.9} K after {steps} steps where each step's `1 − dt·G/C` leaves {discrete:.9} K"
    );

    let want = air + (start - air) * (-elapsed / tau).exp();
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
    // area the caller stated, each behind the half cell of copper between it and the cell centre.
    // That is 2kA/dx, four hundred times the air on a face at 500 W/m²K and six hundred and
    // sixty-eight times the film at 300, so it moves the air alone by 0.25%, both together by
    // 0.19% and the film alone by 0.15%. It is in all three.
    let g_air = per_face(h_air, copper_k(), dx) * (2 * n) as f64;
    let area = n as f64 * dx * dx;
    let half = 2.0 * copper_k() * area / dx;
    let g_film = 1.0 / (1.0 / (h_film * area) + 1.0 / half);
    let closed = watts / (g_air + g_film);

    // And each alone is the other's complement, which is what "in parallel" means and what a term
    // counted twice would break.
    let air_only = settle(build(true, false));
    let film_only = settle(build(false, true));

    // **The tolerance is the settling, and the bar's uniformity makes it nothing else.** Every
    // cell has the same faces, so the steady state is the closed form exactly; `settle` stops when
    // the mean moves by under a nanokelvin a second, which leaves at most that times the time
    // constant, 46 s at the longest, 4.6e-8 K — 1.3e-9 of the smallest rise; measured, 2.6e-10 to
    // 2.9e-10. 1e-7 is eighty times the bound, and the half cell left out of any of the three moves
    // it fifteen thousand times as far.
    let mut wrong = Vec::new();
    for (what, got, want) in [
        ("both paths", both, closed),
        ("air alone", air_only, watts / g_air),
        ("film alone", film_only, watts / g_film),
    ] {
        let miss = (got / want - 1.0).abs();
        println!("  {what}: {got:.7} K against {want:.7} K, {miss:.2e} apart");
        if off_by(miss, 1e-7) {
            wrong.push(format!(
                "{what}: {got:.7} K against {want:.7} K in parallel"
            ));
        }
    }
    if !(both < air_only && both < film_only) {
        wrong.push(format!(
            "two paths must cool better than either: {both:.5} against {air_only:.5}/{film_only:.5}"
        ));
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
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

/// The 5x5x5 block of `a_core_in_a_cavity_sheds_from_the_inside_and_the_ledger_hears_it`: a 3x3x3
/// copper core, every cell of it off the block's outer faces, in air at `h`.
fn cavity(h: f64) -> Solid3D {
    Solid3D::new(
        "core",
        grey(),
        (5, 5, 5),
        Length::from_si(4e-3),
        Temperature::celsius(300.0),
    )
    .empty(|i, j, k| ![i, j, k].iter().all(|a| (1..4).contains(a)))
    .air_in(Temperature::celsius(20.0), h)
}

/// **A core in a cavity takes the step its corners allow, and its corners are inside the block.**
///
/// The limit is the stiffest cell's, `1/max(dx·(Σk_f + G_air/dx)/C)`, and for a cell with `f`
/// faces on the void and the rest on copper that is `ρc dx²/((6 − f)k + f·a)`, where `a` is one
/// face's air behind its half cell, `1/(1/(h·dx) + 1/(2k))`, in W/m/K. Above `a = k`, which is
/// `h·dx = 2k`, the corners, with three faces on the air, are the stiffest: `ρc dx²/(3(k + a))`.
/// At `h = 4e5` that denominator is 2806 against an edge's 2672, a face's 2539 and the centre's
/// 2406. The half cell is also why that is all air can do: however stiff the film, a face of it
/// carries no more than the half cell, `2k`.
///
/// The limit with air was a pass over every cell. It is a pass over the cells facing the void
/// now, which are read from a list `resolve` fills, and a list that missed the corners would hand
/// the core the edge's limit, 1.050 times too long; one that missed them all, the centre's, 1.166
/// times. Every one of these cells is off the outer faces, so the faces' own pass cannot find them.
///
/// The tolerance is rounding, not physics: the closed form is the limit's own arithmetic in
/// another order: about twenty-three roundings on the limit's side, counting the capacity and the
/// faces' series means, and fourteen on the closed form's, the half cell seven of each — 18.5
/// epsilon at half an epsilon each at worst, against a bound of 32 epsilon. It measured none.
#[test]
fn a_core_in_a_cavity_takes_the_step_its_corners_allow() {
    let (h, dx) = (4e5, 4e-3);
    let copper = grey();
    let thermal = copper.thermal.expect("copper conducts");
    let (k, rho, c) = (
        thermal.conductivity.to_si(),
        copper.density.to_si(),
        thermal.specific_heat.to_si(),
    );
    let a = per_face(h, k, dx) / dx;
    assert!(a > k, "the corners set the limit only above a = k");
    let corners = rho * c * dx * dx / (3.0 * (k + a));
    let edges = rho * c * dx * dx / (4.0 * k + 2.0 * a);
    let centre = rho * c * dx * dx / (6.0 * k);

    let limit = cavity(h).max_stable_dt(Time::ZERO).to_si();
    println!(
        "  limit {limit:.9e} s, {:.1} epsilon from the corners' {corners:.9e}; the centre's {centre:.9e}",
        (limit / corners - 1.0).abs() / f64::EPSILON
    );
    assert!(
        (limit / corners - 1.0).abs() <= 32.0 * f64::EPSILON,
        "the core took {limit:.12e} s and its corners allow {corners:.12e}"
    );
    println!(
        "  the edges' limit is {:.4} times the corners', the centre's {:.4}",
        edges / corners,
        centre / corners
    );
    // Measured 1.0499: a list that missed the corners would be 5% long, which the check above sees
    // seven trillion times over. Asserted so that the case cannot drift to one in which it does not.
    assert!(
        edges > corners * 1.04,
        "the corners must set the limit, and by 4% over the edges: {:.4}",
        edges / corners
    );
}

/// **The cells facing the air are the void's as it is now, whatever order made it.**
///
/// `resolve` fills the list of them, and `empty`, `fill` and `air_in` each call it. Built three
/// ways — void then air, air then void, and air, then a void around one cell, then the core filled
/// back in — the same block reports the same limit to the bit. A list filled only where the air
/// was stated would be empty for the second and the third, and one filled again by `empty` but
/// not by `fill` would list one buried cell for the third; each hands the core its centre's limit
/// rather than its corners'.
#[test]
fn the_cells_facing_the_air_are_the_void_it_has_now() {
    let h = 4e5;
    let core = |i: usize, j: usize, k: usize| [i, j, k].iter().all(|a| (1..4).contains(a));
    let block = || {
        Solid3D::new(
            "core",
            grey(),
            (5, 5, 5),
            Length::from_si(4e-3),
            Temperature::celsius(300.0),
        )
    };
    let first = cavity(h).max_stable_dt(Time::ZERO).to_si();

    let air_first = block()
        .air_in(Temperature::celsius(20.0), h)
        .empty(move |i, j, k| !core(i, j, k))
        .max_stable_dt(Time::ZERO)
        .to_si();

    let mut refilled = block()
        .air_in(Temperature::celsius(20.0), h)
        .empty(|i, j, k| (i, j, k) != (2, 2, 2));
    refilled.fill(grey(), core);
    assert_eq!(
        refilled.void_cells(),
        125 - 27,
        "the core is back and nothing else"
    );
    let refilled = refilled.max_stable_dt(Time::ZERO).to_si();

    println!(
        "  void then air {first:.9e} s, air then void {air_first:.9e}, refilled {refilled:.9e}"
    );
    assert_eq!(
        air_first.to_bits(),
        first.to_bits(),
        "air stated before the void"
    );
    assert_eq!(
        refilled.to_bits(),
        first.to_bits(),
        "the core filled back into the void"
    );
}

/// **A film that sets the limit still sets it when the void has air.**
///
/// The limit with air visits two sets of cells, the exposed faces' and those facing the void, and
/// a block that has both can be set by either. A copper slab three cells square, four deep, with
/// its top layer void and a stiff film under it: the film's cells, which face no void, are the
/// stiffest — five copper faces and the film's `2k·h·dx/(2k + h·dx)`, 2673 W/m/K against the
/// centre's six faces, 2406, and the cells under the void's five and a still air's 0.032 behind
/// its half cell, 2005.
/// Air in the void moves none of that, so the limit is the same to the bit with it and without.
/// A limit that read only the cells facing the air, once there was air, would hand the slab the
/// centre's step, 1.11 times too long, and the one scene that ships with air takes its limit from
/// a cell that only conducts, so it would not see it.
#[test]
fn a_film_that_sets_the_limit_still_sets_it_with_air() {
    let dx = 4e-3;
    let slab = |film: bool, air: bool| {
        let mut b = Solid3D::new(
            "slab",
            grey(),
            (3, 3, 4),
            Length::from_si(dx),
            Temperature::celsius(300.0),
        )
        .empty(|_, _, k| k == 3);
        if film {
            b = b.losing_from(
                Face::ZMin,
                Environment {
                    ambient: Temperature::celsius(20.0),
                    convection_w_per_m2_k: 1e6,
                    area: Area::from_si(9.0 * dx * dx),
                },
            );
        }
        if air {
            b = b.air_in(Temperature::celsius(20.0), 8.0);
        }
        b.max_stable_dt(Time::ZERO).to_si()
    };
    let (bare, film, both) = (slab(false, false), slab(true, false), slab(true, true));
    println!("  conduction {bare:.9e} s, with the film {film:.9e}, and air as well {both:.9e}");
    assert!(
        film < bare / 1.1,
        "the film must set the limit: {film:e} against {bare:e}"
    );
    assert_eq!(
        both.to_bits(),
        film.to_bits(),
        "air moved a limit the film sets"
    );
}

/// Which outer layer of a block is void.
type Layer = fn(usize, usize, usize) -> bool;

/// **A cell that faces the air one way only is charged for it, whichever way that is.**
///
/// In the cavity above the stiffest cells are its corners, which face the void three ways at once
/// and between them all six, so a list that looked one way only would still find a corner and the
/// limit would not move. Here the void is one outer layer of a 3x3x3 copper block, and the block's
/// centre then has five copper faces and one on the air, on whichever axis and side the layer is:
/// `ρc dx²/(5k + a)`, with `a` the face's air behind its half cell as above. Every other cell has a
/// face fewer, so the centre sets the limit, and a list blind to that direction hands the block
/// `ρc dx²/(5k)`, 1.266 times too long. The six are one
/// arithmetic over exact terms in a different order, so they agree to the bit; the closed form is
/// held to the 32 epsilon above.
#[test]
fn a_cell_facing_the_air_one_way_sets_the_limit_whichever_way() {
    let (h, dx) = (4e5, 4e-3);
    let copper = grey();
    let thermal = copper.thermal.expect("copper conducts");
    let (k, rho, c) = (
        thermal.conductivity.to_si(),
        copper.density.to_si(),
        thermal.specific_heat.to_si(),
    );
    let closed = rho * c * dx * dx / (5.0 * k + per_face(h, k, dx) / dx);
    let layers: [(&str, Layer); 6] = [
        ("x = 0", |i, _, _| i == 0),
        ("x = 2", |i, _, _| i == 2),
        ("y = 0", |_, j, _| j == 0),
        ("y = 2", |_, j, _| j == 2),
        ("z = 0", |_, _, k| k == 0),
        ("z = 2", |_, _, k| k == 2),
    ];
    let mut first = None;
    for (layer, void) in layers {
        let limit = Solid3D::new(
            "block",
            grey(),
            (3, 3, 3),
            Length::from_si(dx),
            Temperature::celsius(300.0),
        )
        .empty(void)
        .air_in(Temperature::celsius(20.0), h)
        .max_stable_dt(Time::ZERO)
        .to_si();
        println!(
            "  void at {layer}: {limit:.9e} s, {:.1} epsilon from {closed:.9e}",
            (limit / closed - 1.0).abs() / f64::EPSILON
        );
        assert!(
            (limit / closed - 1.0).abs() <= 32.0 * f64::EPSILON,
            "void at {layer}: the block took {limit:.12e} s and its centre allows {closed:.12e}"
        );
        let first = *first.get_or_insert(limit);
        assert_eq!(
            limit.to_bits(),
            first.to_bits(),
            "void at {layer} against x = 0"
        );
    }
}
