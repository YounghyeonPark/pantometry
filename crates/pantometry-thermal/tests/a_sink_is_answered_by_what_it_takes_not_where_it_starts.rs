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
//! **The steady solve read the film's slope at the room**, the hotter of the sink's cell and its
//! room, `4εσA·T∞³`, where what the film does is `4εσA·|s|³` at its surface — in its sweep and in
//! the step on a body's level — so each round left about `ρ = 1 − |s|³/T∞³` of the level still to
//! go and the stop at a microkelvin left `1e-6·ρ/(1 − ρ)`. Measured, the worst cell was that to
//! within 7% from 30 mW to 60 mW: 1.75e-5 K at 50 mW, against 1.76e-5 predicted, and 1.85e-5 K at
//! 52 mW, refused, against 1.86e-5 — more than `SHAPE`, so that a first form of this test held
//! them to the prediction instead. From a little past `|s|³ = 2T∞³` each round overshot by more
//! than it corrected: 180 mW settled, and 185, 190 and 200 mW ran out of sweeps, two hundred
//! anywhere from −1186.6 K to −1434.4 K by start against its chain's −1321.1 K. A cell its film
//! alone holds fared worse: every load from 1 W to 100 W drawn from one reported that it did not
//! converge, from every start. Read at the surface by the step on the level, and in the sweep
//! wherever that is the larger, every sink here is held to `SHAPE` and every one without a
//! balance is refused at its root; a row of 32 cells, where the film is much of its cell's
//! diagonal, is not, and is recorded as open.
//!
//! **And one whose sweeps run out says so**, a long insulating bar's, which the slope does not
//! reach: its film does not radiate. **And one nothing can move says so**: a lone cell in a room at
//! absolute zero, its film with no convection and so no slope at zero, was answered where it
//! started, from a start at zero.

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
/// contraction of 0.91. Measured here, 4.5e-8 K at most across the clearance, 9.7e-7 K from the
/// film of the four cells, at 150 mW, and 6.19e-6 K where a film is much of its cell's diagonal,
/// on 16 cells.
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

/// The cells of [`direct`] by its chain, K.
fn direct_closed(p: f64) -> [f64; 4] {
    let t0 = film_surface(p) - p / half(50.0);
    let t1 = t0 - p / whole(50.0);
    let t2 = t1 - p / whole(50.0);
    [t0, t1, t2, t2 - p / whole(50.0)]
}

/// `b` is answered, every cell of it at `expect` to `SHAPE`; the worst miss. A row, so cell `c` is
/// `(c, 0, 0)`. Only an empty cell is passed over, and `expect` has `NaN` there and nowhere else:
/// a closed form that came out short or not a number would check less without a word.
fn answered(why: &str, b: &Solid3D, expect: &[f64]) -> f64 {
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
            miss.abs() <= SHAPE,
            "{why}: cell {c} at {got:.6} K, {miss:+.3e} K from its closed form {want:.6} K"
        );
        worst = worst.max(miss.abs());
    }
    worst
}

/// `b` is refused as a balance below zero that the sinks set, at `expect` to `SHAPE`, naming
/// `cell` — and as one the sweeps settled on, not one they ran out before, said without the hedge
/// the mirror needed: "or more than reaches them in time for this solve" ended it, and a refusal
/// that only contains the words before that would pass with it. The miss.
fn refused_at(why: &str, b: &Solid3D, expect: f64, cell: &str) -> f64 {
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
        miss.abs() <= SHAPE,
        "{why}: refused at {:.6} K, {miss:+.3e} K from the root below zero, {expect:.6} K",
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
            worst = worst.max(answered(&why, &across(sink, start), &expect));
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
            let miss = refused_at(&why, &across(sink, start), expect, "(4, 0, 0)");
            worst = worst.max(miss);
        }
        println!(
            "  {:4} mW across a clearance: refused at {expect:.6} K from all six starts, the worst \
             {worst:.1e} K from it",
            sink * 1e3
        );
    }
}

/// The film-fed sinks `sinks`, W, each answered from every start at its chain.
fn direct_answered(sinks: &[f64]) {
    for &sink in sinks {
        let expect = direct_closed(sink);
        assert!(expect[3] > 0.0, "{sink} W has a balance above zero");
        let mut worst = 0.0_f64;
        for start in STARTS {
            let why = format!("{} mW from the film from {start} K", sink * 1e3);
            worst = worst.max(answered(&why, &direct(sink, start), &expect));
        }
        println!(
            "  {:5} mW from the film: the sink at {:.6} K from all six starts, the worst cell \
             {worst:.2e} K from its chain",
            sink * 1e3,
            expect[3]
        );
    }
}

/// The film-fed sinks `sinks`, W, each refused from every start at its chain's root.
fn direct_refused(sinks: &[f64]) {
    for &sink in sinks {
        let expect = direct_closed(sink)[3];
        assert!(expect < 0.0, "{sink} W has no balance above zero");
        let mut worst = 0.0_f64;
        for start in STARTS {
            let why = format!("{} mW from the film from {start} K", sink * 1e3);
            worst = worst.max(refused_at(&why, &direct(sink, start), expect, "(3, 0, 0)"));
        }
        println!(
            "  {:5} mW from the film: refused at {expect:.6} K from all six starts, the worst \
             {worst:.2e} K from it",
            sink * 1e3
        );
    }
}

/// **A sink the film feeds directly is answered from every start**, its chain at the film's
/// surface, up to 50.98 mW, 53 µW short of the most the film can bring. Read at the room, the
/// film's slope left 50 mW 1.75e-5 K out.
#[test]
fn a_sink_the_film_feeds_is_answered_from_every_start() {
    direct_answered(&[10e-3, 30e-3, 45e-3, 50e-3, 50.98e-3]);
}

/// **Just past what the film can bring, a sink is refused at its root from every start**, its
/// chain read odd at the film's surface. The even reading refused sixty milliwatts by start
/// anywhere from −15.42 K to −7009.6 K, or as not converging, and fifty-two from −0.32 K to
/// −781.55 K, or as not converging; read at the room, the film's slope left fifty-two 1.85e-5 K
/// from its root.
#[test]
fn a_sink_just_past_what_the_film_can_bring_is_refused_at_its_root() {
    direct_refused(&[52e-3, 53e-3, 55e-3, 60e-3, 70e-3, 100e-3, 150e-3]);
}

/// **Far past it too**, from 185 mW to a watt, the surface more than twice as far below zero as
/// the room is above it at the last. Read at the room, the film's slope overshot here every
/// round: 185, 190 and 200 mW ran out of sweeps.
#[test]
fn a_sink_far_past_what_the_film_can_bring_is_refused_at_its_root() {
    direct_refused(&[185e-3, 200e-3, 300e-3, 1.0]);
}

/// **A cell its film alone holds is refused at its root, however far below zero**: one cell of
/// `k = 50` and a millimetre, its film from a room at 1000 K, from 1 W to 100 W drawn from it, its
/// root `s − P·dx/2kA` below the surface where `εσA(T∞⁴ − s·|s|³) = P`. Its film is all of its
/// diagonal. Read at the room by the sweep and the step on its level, every load reported that it
/// did not converge, from every start; by the sweep alone, two, five and a hundred watts did, and
/// thirty ran out of sweeps from 1500 K. The films' own slope settles each at its root.
#[test]
fn a_cell_its_film_alone_holds_is_refused_at_its_root() {
    let sigma = STEFAN_BOLTZMANN.to_si();
    for p in [1.0, 2.0, 5.0, 30.0, 100.0] {
        let room4 = (ROOM * ROOM) * (ROOM * ROOM);
        let expect = odd_root(room4 - p / (0.9 * sigma * AREA)) - p / half(50.0);
        let mut worst = 0.0_f64;
        for start in [1.0, 300.0, 1000.0, 1500.0] {
            let b = Solid3D::new(
                "block",
                stuff(50.0, 0.9),
                (1, 1, 1),
                Length::from_si(DX),
                Temperature::from_si(start),
            )
            .losing_from(Face::XMin, radiating())
            .dissipating(-p, |_, _, _| true);
            let why = format!("{p} W from a lone cell from {start} K");
            worst = worst.max(refused_at(&why, &b, expect, "(0, 0, 0)"));
        }
        println!(
            "  {p:5} W from a lone cell: refused at {expect:.6} K, the worst {worst:.2e} K from it"
        );
    }
}

/// **A film that is much of its cell's diagonal is answered to `SHAPE`**: a row of 4 and 16 cells
/// of `k = 1.76` and a centimetre, a film with no convection on the first from a room at 1200 K,
/// and 9.92 W drawn from that cell, so that every cell sits at `s − P·dx/2kA`, 318.38 K — the film
/// 18% of the first cell's diagonal, and half of it as the limit's tangent at the room reads it.
/// Read at the room by the sweep as well as by the step on the level, the stop left 1.19e-5 K at
/// 16 cells, and now 6.185e-6 K.
///
/// **At 32 cells it is not, and that is recorded as open.** The sweep reads the larger of the
/// film's slope and that tangent, which here is the tangent, and the stop leaves 1.23e-5 K, the
/// worst cell, measured; 2.43e-5 K with the tangent in the step on the level too, and 6.19e-6 K
/// with the film's slope in the sweep alone — which, from below a convex `T⁴`, overshoots:
/// `a_start_far_below_its_balance_is_answered` has where.
#[test]
fn a_film_that_is_much_of_its_cells_diagonal_is_answered() {
    let sigma = STEFAN_BOLTZMANN.to_si();
    let (k, dx, p, room) = (1.76, 1e-2, 9.92, 1200.0);
    let area = dx * dx;
    let room4 = (room * room) * (room * room);
    let t = odd_root(room4 - p / (0.9 * sigma * area)) - p / (2.0 * k * area / dx);
    for n in [4, 16] {
        let b = Solid3D::new(
            "row",
            stuff(k, 0.9),
            (n, 1, 1),
            Length::from_si(dx),
            Temperature::celsius(20.0),
        )
        .losing_from(
            Face::XMin,
            Environment {
                ambient: Temperature::from_si(room),
                convection_w_per_m2_k: 0.0,
                area: Area::from_si(area),
            },
        )
        .dissipating(-p, |i, _, _| i == 0);
        let worst = answered(&format!("{n} cells"), &b, &vec![t; n]);
        println!("  {n:2} cells at {t:.6} K: the worst cell {worst:.2e} K from it");
    }
}

/// **A cell nothing can move is refused, not answered**: a lone cell in a room at absolute zero,
/// its film with no convection, so that at a surface of zero it has no slope at all, drawing
/// 10 mW, which no film from a room at zero can bring. Measured before this: from zero, 1e-300 K
/// and 1e-120 K the solve passed it over, moved nothing else, and answered where it started,
/// 665.43 K above its root; from 1e-40 K to 300 K it was refused as settled below zero, but
/// wherever its steps had left it — −4.9e46 K from 1e-12 K, the first sweep's throw, and
/// −1589.45 K from 300 K — and not at its root. Refused now, every start, and never answered:
/// - **as a cell left out of balance with nothing to move it by**, carrying its 10 mW, from zero, from
///   1e-300 K and 1e-120 K, where its slope is below the least `f64`, and from 1e-12 K and 1e-6 K,
///   thrown as before and brought back to exactly zero by the next sweep, which reads the film's
///   slope there;
/// - **as drawn down by its sinks until its step overflowed** from 5e-99 K to 5e-23 K, measured,
///   where that first throw is past where `T⁴` overflows; refused there once as diverging, a poorer
///   word for it — from 3e-99 K down the slope is below the least `f64` and it is refused as
///   stuck, as from 1e-120 K here;
/// - **at its root**, settled, from 1e-3 K, 1 K and 300 K.
#[test]
fn a_cell_nothing_can_move_is_refused() {
    let sigma = STEFAN_BOLTZMANN.to_si();
    let root = odd_root(-0.01 / (0.9 * sigma * AREA)) - 0.01 / half(50.0);
    for start in [0.0, 1e-300, 1e-120, 1e-40, 1e-12, 1e-6, 1e-3, 1.0, 300.0] {
        let b = Solid3D::new(
            "block",
            stuff(50.0, 0.9),
            (1, 1, 1),
            Length::from_si(DX),
            Temperature::from_si(start),
        )
        .losing_from(
            Face::XMin,
            Environment {
                ambient: Temperature::from_si(0.0),
                convection_w_per_m2_k: 0.0,
                area: Area::from_si(AREA),
            },
        )
        .dissipating(-0.01, |_, _, _| true);
        let e = match b.steady_state() {
            Ok(t) => panic!("from {start} K: no balance, and answered {t:?}"),
            Err(e) => e,
        };
        println!("  from {start:e} K -> {e}");
        // Where it was left, as the refusal says: at zero, or where it started if nothing moved
        // it at all.
        let left_at = if start == 0.0 || start == 1e-300 || start == 1e-120 {
            format!("{start:e}")
        } else {
            "0e0".to_string()
        };
        let stuck = e.quantity
            == format!(
                "no balance was found: cell (0, 0, 0) is left -1e-2 W out of balance at {left_at} \
                 K, where nothing it exchanges has a slope to move it by"
            )
            && e.before == -0.01;
        let overflowed = e.quantity
            == "no balance above absolute zero was found: cell (0, 0, 0) was drawn down by the \
                sinks in its part until the step on its level overflowed";
        let settled = e.quantity.starts_with(
            "no balance above absolute zero was found: cell (0, 0, 0) came to below it",
        ) && e.quantity.ends_with("which take more than can reach them")
            && (e.before - root).abs() <= SHAPE;
        // What each start comes to, measured: the first sweep throws a start near zero by the
        // watts over a slope of nearly nothing.
        let right = if start >= 1e-3 {
            settled
        } else if start == 1e-40 {
            overflowed
        } else {
            stuck
        };
        assert!(
            right,
            "from {start} K: refused for the wrong reason, or not at its root {root} K: {e}"
        );
    }
}

/// **The cell named is the first left out of balance, where it is**: two rows of five cells along
/// z, the middle three conducting nothing, films with no convection on both ends from a room at
/// absolute zero, from a start at zero. With 10 mW drawn from the bottom and 20 mW from the top,
/// the bottom is named, carrying its 10 mW; with the top's 20 mW alone, the top, at (0, 0, 4) and
/// not at (4, 0, 0); and the three between, which have nothing to move and nothing out of
/// balance, are not.
#[test]
fn a_cell_nothing_can_move_is_named() {
    for (bottom, top, named, watts) in [
        (-0.01, -0.02, "(0, 0, 0)", "-1e-2"),
        (0.0, -0.02, "(0, 0, 4)", "-2e-2"),
    ] {
        let zero = || Environment {
            ambient: Temperature::from_si(0.0),
            convection_w_per_m2_k: 0.0,
            area: Area::from_si(AREA),
        };
        let mut b = Solid3D::new(
            "column",
            stuff(50.0, 0.9),
            (1, 1, 5),
            Length::from_si(DX),
            Temperature::from_si(0.0),
        )
        .losing_from(Face::ZMin, zero())
        .losing_from(Face::ZMax, zero());
        b.fill(stuff(0.0, 0.9), |_, _, k| (1..=3).contains(&k));
        let b = b
            .dissipating(bottom, |_, _, k| k == 0)
            .dissipating(top, |_, _, k| k == 4);
        let e = b
            .steady_state()
            .expect_err("no film from a room at zero can feed a sink");
        println!("  {bottom} W and {top} W -> {e}");
        assert_eq!(
            e.quantity,
            format!(
                "no balance was found: cell {named} is left {watts} W out of balance at 0e0 K, \
                 where nothing it exchanges has a slope to move it by"
            ),
            "refused for the wrong reason, or naming the wrong cell"
        );
        assert_eq!(e.before, if bottom != 0.0 { bottom } else { top });
    }
}

/// **A sink whose sweeps run out before it settles is refused as that**, and says no more than
/// was found: four hundred cells of an insulator, `k = 0.05`, a millimetre each, a film of
/// 10 W/m²·K that does not radiate on the first from a room at 300 K, and a milliwatt drawn from
/// the last, through `1/hA + dx/2kA + 399·dx/kA` = 8.09e6 K/W from the room, so that its root is
/// −7790.0 K. Measured: refused at −7734.9 K, the sweeps run out on the way, as they did before
/// the solve read the films' own slope — by construction, since this film is linear and its slope
/// is what the limit charges, to the bit. **It is held there by how slowly the solve settles a
/// long bar**, not by anything wrong: one that settled it would turn this red, and want a longer
/// bar.
#[test]
fn a_sink_whose_sweeps_run_out_says_so() {
    let b = Solid3D::new(
        "bar",
        stuff(0.05, 0.0),
        (400, 1, 1),
        Length::from_si(DX),
        Temperature::celsius(20.0),
    )
    .losing_from(
        Face::XMin,
        Environment {
            ambient: Temperature::from_si(300.0),
            convection_w_per_m2_k: 10.0,
            area: Area::from_si(AREA),
        },
    )
    .dissipating(-1e-3, |i, _, _| i == 399);
    let e = b.steady_state().unwrap_err();
    println!("  four hundred cells -> {e}");
    assert!(
        e.quantity
            .starts_with("no balance above absolute zero was found: cell (399, 0, 0) was below it")
            && e.quantity.ends_with(
                "drawn down by the sinks in its part, where 2000 sweeps had not settled"
            ),
        "refused for the wrong reason: {e}"
    );
    assert!(e.before < 0.0, "refused at {} K", e.before);
}
