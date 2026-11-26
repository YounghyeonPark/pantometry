//! `Solid3D::steady_state` of a sink: **answered by what it takes, and not by where the solve
//! starts.** A sink with a balance is answered from every start, at its closed form; one without is
//! refused from every start at one temperature below absolute zero — the closed form's root, with
//! each fourth power read as `T·|T|³` — and the refusal says the sinks take more than can reach
//! them.
//!
//! **The pairs and the films read a negative temperature's fourth power as even**, so the balance
//! had a mirror below zero, and a sink drawn through zero before its heat arrived settled on it.
//! Measured before this, from starts of 1, 100, 300, 600, 1000 and 1500 K, on the two blocks here —
//! a sink fed across a clearance by a body a film warms from a room at 1000 K, and a sink that film
//! feeds directly:
//!
//! - **Of ten sinks with a balance, five runs of sixty were refused below zero**: half a milliwatt
//!   from 1 K at −954.17 K, a milliwatt from 300 K at −898.56 K, one and a half from 100 K and from
//!   300 K at −827.16 K, two and a half from 100 K at −516.95 K. Each is the mirror of its sink's
//!   surface less the drop across its half cell: the milliwatt's surface is at 888.56 K and its
//!   cell 10 K below it, at 878.56 K, where the other five starts put it.
//! - **Of six without one, nineteen runs of thirty-six reported not converging**, and seventeen were
//!   refused below zero, anywhere from −0.32 K to −7009.6 K by start.
//!
//! # What is checked against what
//!
//! **Each answer is its chain's closed form**, every cell of it. The film's surface sits where
//! `εσA(T∞⁴ − s⁴)` is the sink's `P`; each cell is below the one before by `P` over what joins
//! them, half a cell `2kA/dx` from a surface and a whole one `kA/dx` from a cell; across the
//! clearance, the far surface `v` is where `C(u⁴ − v⁴) = P` with `C = σA/(1/ε_a + 1/ε_b − 1)`; and
//! the cell past the sink, which nothing crosses, sits at the sink's. **Each refusal is the same
//! chain read odd**, `g(v) = u⁴ − P/C` with `g(v) = v|v|³`, whose root is below zero where `u⁴` is
//! less than `P/C`, and the temperature refused is that chain's sink cell. Held to `SHAPE`, from
//! all six starts, so that an answer that moves with the start fails at the starts it moves at.
//!
//! **A sink the film feeds is held to the stop's own error at its own contraction**, which near the
//! most the film can bring, 51.03 mW, is more than `SHAPE`. The step on a body's level reads the
//! film's slope at the hotter of its cell and the room, here the room, `4εσA·T∞³`, where what the
//! film does is `4εσA·|s|³` at its surface, so each round leaves `ρ = 1 − |s|³/T∞³` of the level
//! still to go, and the stop at a microkelvin leaves `1e-6·ρ/(1 − ρ)`. Measured, the worst cell is
//! that to within 7%, and not over it, from 30 mW to 60 mW, on both sides of the limit: 1.75e-5 K
//! at 50 mW, answered, against 1.76e-5 predicted, and to the printed digits as before this;
//! 1.85e-5 K at 52 mW, refused, against 1.86e-5. Each is held to half again its prediction, or
//! `SHAPE` where that is more, which is everywhere but 50, 52 and 53 mW.
//! **A sink four times the most is not here**, 200 mW: its sweeps run out with its sink still
//! drawn down, anywhere from −1186.6 K to −1434.4 K by start against the chain's −1321.1 K, and
//! it is refused as not having settled, which says no more than was found.

use pantometry_core::substance::ThermalProps;
use pantometry_core::units::{
    Area, Density, Length, SpecificHeat, Temperature, ThermalConductivity, ThermalExpansion,
    STEFAN_BOLTZMANN,
};
use pantometry_core::Substance;
use pantometry_thermal::{Environment, Face, Solid3D};

/// How far a cell may be from its closed form, K, as `a_part_that_loses_nothing_keeps_the_heat_it_
/// holds` holds one. The solve stops when a round moves the field by less than a microkelvin,
/// which leaves `1e-6·ρ/(1 − ρ)` for a round that contracts the error by `ρ`: ten microkelvin is a
/// contraction of 0.91. Measured here, 1.0e-7 K at most across the clearance, and from the film
/// 5.7e-6 K where this is what it is held to, at 55 mW.
const SHAPE: f64 = 1e-5;

/// The starts: from near absolute zero to half again the room, through where each refusal above
/// was found.
const STARTS: [f64; 6] = [1.0, 100.0, 300.0, 600.0, 1000.0, 1500.0];

/// The room the film radiates from, K, with no convection, over a face of 1 mm².
const ROOM: f64 = 1000.0;
const AREA: f64 = 1e-6;
const DX: f64 = 1e-3;

fn stuff(k: f64, emissivity: f64) -> Substance {
    Substance::bulk("stuff", Density::from_si(1000.0)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(k),
        specific_heat: SpecificHeat::j_per_kg_k(4000.0),
        expansion: ThermalExpansion::ppm_per_k(50.0),
        emissivity,
    })
}

fn radiating() -> Environment {
    Environment {
        ambient: Temperature::from_si(ROOM),
        convection_w_per_m2_k: 0.0,
        area: Area::from_si(AREA),
    }
}

/// `y` with `y·|y|³ = x`: the fourth root, read odd. Square roots only, which round correctly.
fn odd_root(x: f64) -> f64 {
    if x < 0.0 {
        -(-x).sqrt().sqrt()
    } else {
        x.sqrt().sqrt()
    }
}

/// The film's surface when `p` W reaches the block through it from a face of emissivity 0.9:
/// `εσA(T∞⁴ − g(s)) = p`.
fn film_surface(p: f64) -> f64 {
    let room4 = (ROOM * ROOM) * (ROOM * ROOM);
    odd_root(room4 - p / (0.9 * STEFAN_BOLTZMANN.to_si() * AREA))
}

/// What a face of `k` between a cell and its neighbour conducts, W/K, and between a cell and its
/// surface.
fn whole(k: f64) -> f64 {
    k * AREA / DX
}
fn half(k: f64) -> f64 {
    2.0 * k * AREA / DX
}

/// **A sink fed across a clearance**, the block `a_part_that_loses_nothing_keeps_the_heat_it_holds`
/// solves: six cells of 1 mm in a row, the first three of `k = 5`, `ε = 0.9` warmed by the film,
/// the fourth empty, the fifth the sink, of `k = 0.05`, `ε = 0.05`, and the sixth of `k = 50`
/// beyond it, with nothing to lose to.
fn across(sink: f64, start: f64) -> Solid3D {
    let mut b = Solid3D::new(
        "block",
        stuff(5.0, 0.9),
        (6, 1, 1),
        Length::from_si(DX),
        Temperature::from_si(start),
    )
    .empty(|i, _, _| i == 3)
    .losing_from(Face::XMin, radiating());
    b.fill(stuff(0.05, 0.05), |i, _, _| i == 4);
    b.fill(stuff(50.0, 0.9), |i, _, _| i == 5);
    b.dissipating(-sink, |i, _, _| i == 4)
}

/// The cells of [`across`] by its chain, K, the empty one `NaN`.
fn across_closed(p: f64) -> [f64; 6] {
    let t0 = film_surface(p) - p / half(5.0);
    let t1 = t0 - p / whole(5.0);
    let t2 = t1 - p / whole(5.0);
    let u = t2 - p / half(5.0);
    let c = STEFAN_BOLTZMANN.to_si() * AREA / (1.0 / 0.9 + 1.0 / 0.05 - 1.0);
    let v = odd_root((u * u) * (u * u) - p / c);
    let t4 = v - p / half(0.05);
    [t0, t1, t2, f64::NAN, t4, t4]
}

/// **A sink the film feeds directly**: four cells of `k = 50`, `ε = 0.9` in a row, the film on the
/// first and the sink the last. At most `εσA·T∞⁴` = 51.03 mW can come in.
fn direct(sink: f64, start: f64) -> Solid3D {
    Solid3D::new(
        "block",
        stuff(50.0, 0.9),
        (4, 1, 1),
        Length::from_si(DX),
        Temperature::from_si(start),
    )
    .losing_from(Face::XMin, radiating())
    .dissipating(-sink, |i, _, _| i == 3)
}

/// How far a cell of [`direct`] taking `p` W may be from its chain: half again the stop's error at
/// its film's contraction, `1e-6·ρ/(1 − ρ)` with `ρ = 1 − |s|³/T∞³`, or `SHAPE` where that is
/// more. The half cell behind the film, 0.1 W/K, carries what it does unchanged beside the
/// film's 1e-5 W/K, so the film's slope is its surface's.
fn direct_allowed(p: f64) -> f64 {
    let s = film_surface(p).abs() / ROOM;
    let rho = 1.0 - s * s * s;
    SHAPE.max(1.5e-6 * rho / (1.0 - rho))
}

/// The cells of [`direct`] by its chain, K.
fn direct_closed(p: f64) -> [f64; 4] {
    let t0 = film_surface(p) - p / half(50.0);
    let t1 = t0 - p / whole(50.0);
    let t2 = t1 - p / whole(50.0);
    [t0, t1, t2, t2 - p / whole(50.0)]
}

/// `b` is answered, every cell of it at `expect` to `allowed`; the worst miss. A row, so cell `c` is
/// `(c, 0, 0)`. Only an empty cell is passed over, and `expect` has `NaN` there and nowhere else:
/// a closed form that came out short or not a number would check less without a word.
fn answered(why: &str, b: &Solid3D, expect: &[f64], allowed: f64) -> f64 {
    let t = b
        .steady_state()
        .unwrap_or_else(|e| panic!("{why}: there is a balance to strike: {e}"));
    assert_eq!(t.len(), expect.len(), "{why}: a closed form for every cell");
    let mut worst = 0.0_f64;
    for (c, (got, want)) in t.iter().zip(expect).enumerate() {
        assert_eq!(
            want.is_nan(),
            b.is_void(c, 0, 0),
            "{why}: cell {c}'s closed form is {want}"
        );
        if want.is_nan() {
            continue;
        }
        let miss = got - want;
        assert!(
            miss.abs() <= allowed,
            "{why}: cell {c} at {got:.6} K, {miss:+.3e} K from its closed form {want:.6} K, \
             against {allowed:.2e}"
        );
        worst = worst.max(miss.abs());
    }
    worst
}

/// `b` is refused as a balance below zero that the sinks set, at `expect` to `allowed`, naming
/// `cell` — and as one the sweeps settled on, not one they ran out before, said without the hedge
/// the mirror needed: "or more than reaches them in time for this solve" ended it, and a refusal
/// that only contains the words before that would pass with it. The miss.
fn refused_at(why: &str, b: &Solid3D, expect: f64, cell: &str, allowed: f64) -> f64 {
    let e = match b.steady_state() {
        Ok(t) => {
            panic!("{why}: there is no balance above absolute zero, and it was answered {t:?}")
        }
        Err(e) => e,
    };
    assert!(
        e.quantity.starts_with(&format!(
            "no balance above absolute zero was found: cell {cell} came to below it"
        )) && e
            .quantity
            .ends_with("drawn down by the sinks in its part, which take more than can reach them"),
        "{why}: refused for the wrong reason: {e}"
    );
    let miss = e.before - expect;
    assert!(
        miss.abs() <= allowed,
        "{why}: refused at {:.6} K, {miss:+.3e} K from the root below zero, {expect:.6} K, against \
         {allowed:.2e}",
        e.before
    );
    miss.abs()
}

/// **A sink fed across a clearance is answered from every start**, at its chain, where the even
/// reading refused five of the thirty-six runs here at the mirror of their balance.
#[test]
fn a_sink_fed_across_a_clearance_is_answered_from_every_start() {
    for sink in [0.1e-3, 0.5e-3, 1e-3, 1.5e-3, 2e-3, 2.5e-3] {
        let expect = across_closed(sink);
        assert!(expect[4] > 0.0, "{sink} W has a balance above zero");
        let mut worst = 0.0_f64;
        for start in STARTS {
            let why = format!("{} mW across a clearance from {start} K", sink * 1e3);
            worst = worst.max(answered(&why, &across(sink, start), &expect, SHAPE));
        }
        println!(
            "  {:4} mW across a clearance: the sink at {:.6} K from all six starts, the worst cell \
             {worst:.1e} K from its chain",
            sink * 1e3,
            expect[4]
        );
    }
}

/// **A sink fed across a clearance that takes more than can reach it is refused at one root from
/// every start**, its sink's cell where the chain read odd puts it — where the even reading
/// reported not converging, or ran out of sweeps below zero anywhere from −81.95 K to −174.2 K.
#[test]
fn a_sink_beyond_what_crosses_is_refused_at_one_root_from_every_start() {
    for sink in [2.7e-3, 3e-3, 5e-3] {
        let expect = across_closed(sink)[4];
        assert!(expect < 0.0, "{sink} W has no balance above zero");
        let mut worst = 0.0_f64;
        for start in STARTS {
            let why = format!("{} mW across a clearance from {start} K", sink * 1e3);
            // The cell past the sink sits at its temperature, and the sink's is named: the first
            // of the coldest, which every platform rounds alike.
            let miss = refused_at(&why, &across(sink, start), expect, "(4, 0, 0)", SHAPE);
            worst = worst.max(miss);
        }
        println!(
            "  {:4} mW across a clearance: refused at {expect:.6} K from all six starts, the worst \
             {worst:.1e} K from it",
            sink * 1e3
        );
    }
}

/// **A sink the film feeds directly is answered from every start, or refused at one root**, its
/// chain read odd at its film's surface: up to 50 mW answered, and from 52 mW, just past what the
/// film can bring, to 150 mW, whose surface is further below zero than the room is above it,
/// refused. The even reading refused sixty milliwatts by start anywhere from −15.42 K to
/// −7009.6 K, or as not converging, and fifty-two from −0.32 K to −781.55 K, or as not
/// converging.
#[test]
fn a_sink_the_film_feeds_is_answered_or_refused_by_what_it_takes() {
    for sink in [10e-3, 30e-3, 45e-3, 50e-3] {
        let (expect, allowed) = (direct_closed(sink), direct_allowed(sink));
        assert!(expect[3] > 0.0, "{sink} W has a balance above zero");
        let mut worst = 0.0_f64;
        for start in STARTS {
            let why = format!("{} mW from the film from {start} K", sink * 1e3);
            worst = worst.max(answered(&why, &direct(sink, start), &expect, allowed));
        }
        println!(
            "  {:4} mW from the film: the sink at {:.6} K from all six starts, the worst cell \
             {worst:.2e} K from its chain against {allowed:.2e}",
            sink * 1e3,
            expect[3]
        );
    }
    for sink in [52e-3, 53e-3, 55e-3, 60e-3, 70e-3, 100e-3, 150e-3] {
        let (expect, allowed) = (direct_closed(sink)[3], direct_allowed(sink));
        assert!(expect < 0.0, "{sink} W has no balance above zero");
        let mut worst = 0.0_f64;
        for start in STARTS {
            let why = format!("{} mW from the film from {start} K", sink * 1e3);
            let b = direct(sink, start);
            worst = worst.max(refused_at(&why, &b, expect, "(3, 0, 0)", allowed));
        }
        println!(
            "  {:4} mW from the film: refused at {expect:.6} K from all six starts, the worst \
             {worst:.2e} K from it against {allowed:.2e}",
            sink * 1e3
        );
    }
}
