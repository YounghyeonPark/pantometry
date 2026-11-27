//! `Solid3D::steady_state` from a start far below its balance: **answered there, at its closed
//! form**, every cell of it — and from a start at zero, where nothing has a slope, refused by name.
//!
//! **The sweep's diagonal is the larger of the films' own slope and the limit's tangent**, and the
//! slope alone was tried first. It is Newton on each cell, and from below a convex `T⁴` a Newton
//! step overshoots, and over-relaxed it overshoots by more; the tangent, read at the hotter of a
//! cell and its room, is at least the slope at and above zero and is what the sweep always read.
//! Measured with the slope alone, of the cases here that the tangent answered:
//!
//! - a strip of 64 cells radiating from four sides did not converge from 300 K or 150 K;
//! - a lone cell did not converge from 0 K or 1e-6 K to a room at 300 K, nor from 1e-6 K to one at
//!   20 K;
//! - a heated cell in a room at 0 K, from 1e-12 K to 3e-3 K, came to exactly zero, where nothing
//!   has a slope, and was refused 665.43 K from its balance.
//!
//! **And a step held short is not a balance.** The step on a body's level is a quarter of the
//! larger of its own level and the block's, so a block near absolute zero takes steps as small as
//! it is, and a round that moved it less than a microkelvin was taken as settled: two cells from
//! 0 K in a room at 20 K were answered at 4.59e-7 K, four at 2.3e-7 K to 3.7e-7 K, sixteen at
//! 1.3e-7 K to 5.5e-7 K, eight in a room at 10 K at 1.0e-8 K to 2.8e-8 K, two in a room at 1 K at
//! 2.9e-12 K, and two bodies across a clearance, a film on each, at 2.7e-7 K to 4.4e-7 K — each
//! where its balance is its room. That was before either change, and the slope alone answered the
//! rows the same. Two bodies with a film on the first alone were answered so for another reason,
//! and are in `a_step_judged_singular_is_not_a_balance`.

use pantometry_core::substance::ThermalProps;
use pantometry_core::units::{
    Area, Density, Length, SpecificHeat, Temperature, ThermalConductivity, ThermalExpansion,
    STEFAN_BOLTZMANN,
};
use pantometry_core::Substance;
use pantometry_thermal::{Environment, Face, Solid3D};

/// How far a cell may be from its closed form, K, as
/// `a_sink_is_answered_by_what_it_takes_not_where_it_starts` holds one: the solve stops when a
/// round moves the field by less than a microkelvin, which leaves `1e-6·ρ/(1 − ρ)` for a round that
/// contracts the error by `ρ`, and ten microkelvin is a contraction of 0.91. **Borrowed, and not
/// earned here**: measured, 6.4e-7 K at most, on the strip, and none of the cases here was seen to
/// contract that slowly. It is the stop's error the crate holds a steady cell to, and what this
/// file guards against does not lean on it: every failure recorded above is a refusal, a run out
/// of sweeps or a miss of a kelvin or more.
const SHAPE: f64 = 1e-5;

fn stuff(k: f64, emissivity: f64) -> Substance {
    Substance::bulk("stuff", Density::from_si(1000.0)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(k),
        specific_heat: SpecificHeat::j_per_kg_k(4000.0),
        expansion: ThermalExpansion::ppm_per_k(50.0),
        emissivity,
    })
}

/// A film with no convection from a room at `t` K, over `area` m².
fn room(t: f64, area: f64) -> Environment {
    Environment {
        ambient: Temperature::from_si(t),
        convection_w_per_m2_k: 0.0,
        area: Area::from_si(area),
    }
}

/// `b`, a row of `cells` cells of which `voids` are empty, is answered with every solid one at
/// `want` to `SHAPE`; the worst miss. The empty cells are counted, so that a block emptied by
/// mistake checks nothing without a word, and a miss that is not a number fails, as it should.
fn answered(why: &str, b: &Solid3D, cells: usize, voids: usize, want: f64) -> f64 {
    let t = b
        .steady_state()
        .unwrap_or_else(|e| panic!("{why}: there is a balance to strike: {e}"));
    assert_eq!(t.len(), cells, "{why}: a temperature for every cell");
    let (mut worst, mut skipped) = (0.0_f64, 0);
    for (c, got) in t.iter().enumerate() {
        if b.is_void(c, 0, 0) {
            skipped += 1;
            continue;
        }
        let miss = got - want;
        assert!(
            miss.abs() <= SHAPE,
            "{why}: cell {c} at {got:e} K, {miss:+.3e} K from its balance {want} K"
        );
        worst = worst.max(miss.abs());
    }
    assert_eq!(skipped, voids, "{why}: the empty cells");
    worst
}

/// **A strip radiating from four sides is answered from below its balance**: 64 cells of `k = 1`
/// and a centimetre, `ε = 0.9`, films with no convection on its four long faces from a room at
/// 300 K, its ends insulated, heated so that every surface is at 600 K — every cell at
/// `600 + εσ·dx·(600⁴ − 300⁴)/2k`, 631.0028 K, what crosses its half cell to each face. The films'
/// slope alone did not converge from 300 K or 150 K.
#[test]
fn a_strip_radiating_from_four_sides_is_answered_from_below() {
    let sigma = STEFAN_BOLTZMANN.to_si();
    let (n, dx, k) = (64, 1e-2, 1.0);
    let side = dx * dx * n as f64;
    let rise = (600.0_f64 * 600.0) * (600.0 * 600.0) - (300.0_f64 * 300.0) * (300.0 * 300.0);
    let watts = 0.9 * sigma * (4.0 * side) * rise;
    let want = 600.0 + 0.9 * sigma * dx * rise / (2.0 * k);
    for start in [150.0, 300.0, 600.0] {
        let b = Solid3D::new(
            "strip",
            stuff(k, 0.9),
            (n, 1, 1),
            Length::from_si(dx),
            Temperature::from_si(start),
        )
        .losing_from(Face::YMin, room(300.0, side))
        .losing_from(Face::YMax, room(300.0, side))
        .losing_from(Face::ZMin, room(300.0, side))
        .losing_from(Face::ZMax, room(300.0, side))
        .dissipating(watts, |_, _, _| true);
        let why = format!("the strip from {start} K");
        let worst = answered(&why, &b, n, 0, want);
        println!("  {why}: every cell at {want:.6} K, the worst {worst:.2e} K from it");
    }
}

/// **A lone cell below its room is answered from near zero, at the room**: a millimetre cube,
/// `ε = 0.9`, one film with no convection, nothing generated. The films' slope there is that of a
/// surface near zero, nearly nothing, and alone it threw the cell far past the room: it did not
/// converge from 0 K or 1e-6 K to a room at 300 K at `k = 400`, nor from 1e-6 K to one at 20 K at
/// `k = 50`.
#[test]
fn a_cell_below_its_room_is_answered_from_near_zero() {
    for (k, at, start) in [
        (400.0, 300.0, 0.0),
        (400.0, 300.0, 1e-6),
        (50.0, 20.0, 1e-6),
        (50.0, 300.0, 0.0),
    ] {
        let b = Solid3D::new(
            "cell",
            stuff(k, 0.9),
            (1, 1, 1),
            Length::from_si(1e-3),
            Temperature::from_si(start),
        )
        .losing_from(Face::XMin, room(at, 1e-6));
        let why = format!("k = {k} from {start:e} K to a room at {at} K");
        let worst = answered(&why, &b, 1, 0, at);
        println!("  {why}: {worst:.2e} K from it");
    }
}

/// **A heated cell in a room at absolute zero is answered from near zero**: a millimetre cube of
/// `k = 50`, `ε = 0.9`, one film with no convection from a room at 0 K, 10 mW generated in it — at
/// `s + P·dx/2kA`, 665.43 K, where `εσA·s⁴ = P`. The films' slope alone brought it to exactly zero
/// from 1e-12 K to 3e-3 K, and refused it there. **From exactly zero it is refused, by name**:
/// nothing it exchanges has a slope there, so nothing can move it, and the refusal says so,
/// carrying its 10 mW. It was answered at 0 K.
#[test]
fn a_heated_cell_in_a_room_at_zero_is_answered_from_near_zero() {
    let sigma = STEFAN_BOLTZMANN.to_si();
    let (p, area, dx, k) = (0.01, 1e-6, 1e-3, 50.0);
    let want = (p / (0.9 * sigma * area)).sqrt().sqrt() + p * dx / (2.0 * k * area);
    let cell = |start: f64| {
        Solid3D::new(
            "cell",
            stuff(k, 0.9),
            (1, 1, 1),
            Length::from_si(dx),
            Temperature::from_si(start),
        )
        .losing_from(Face::XMin, room(0.0, area))
        .dissipating(p, |_, _, _| true)
    };
    for start in [1e-12, 1e-6, 1e-3, 3e-3, 1.0] {
        let why = format!("10 mW from {start:e} K");
        let worst = answered(&why, &cell(start), 1, 0, want);
        println!("  {why}: at {want:.6} K, {worst:.2e} K from it");
    }
    let e = cell(0.0)
        .steady_state()
        .expect_err("nothing has a slope to move a cell at zero");
    println!("  10 mW from 0 K -> {e}");
    assert_eq!(
        e.quantity,
        "no balance was found: cell (0, 0, 0) is left 1e-2 W out of balance at 0e0 K, where \
         nothing it exchanges has a slope to move it by",
        "refused for the wrong reason"
    );
    assert_eq!(e.before, p);
}

/// **A block near absolute zero is not answered by a step held short**: rows of 1 mm cells of
/// `k = 50`, `ε = 0.9`, from 0 K, a film with no convection on the first from a cold room, nothing
/// generated — every cell at the room. Each was answered within a microkelvin of zero, the step on
/// its level held to a quarter of the little it had.
#[test]
fn a_block_near_zero_is_not_answered_by_a_step_held_short() {
    for (n, at) in [(2, 20.0), (4, 20.0), (16, 20.0), (8, 10.0), (2, 1.0)] {
        let b = Solid3D::new(
            "row",
            stuff(50.0, 0.9),
            (n, 1, 1),
            Length::from_si(1e-3),
            Temperature::from_si(0.0),
        )
        .losing_from(Face::XMin, room(at, 1e-6));
        let why = format!("{n} cells from 0 K in a room at {at} K");
        let worst = answered(&why, &b, n, 0, at);
        println!("  {why}: {worst:.2e} K from it");
    }
}

/// **Nor are the bodies a clearance joins**, whose step on their levels is one, held short as a
/// whole: two bodies of two such cells, a clearance between them, from 0 K, a film on each end
/// from a room at 20 K — every solid cell at the room. They were answered at 2.7e-7 K to
/// 4.4e-7 K.
///
/// **With a film on the first alone the step was not held short but not taken**: the second body,
/// which only the clearance reaches, has next to no slope to its pair near 0 K, and the step the
/// two take together was judged singular — answered at 3.27e-7 K and 6.60e-40 K. That is
/// `a_step_judged_singular_is_not_a_balance`.
#[test]
fn bodies_a_clearance_joins_near_zero_are_not_answered_by_a_step_held_short() {
    let b = Solid3D::new(
        "pair",
        stuff(50.0, 0.9),
        (5, 1, 1),
        Length::from_si(1e-3),
        Temperature::from_si(0.0),
    )
    .empty(|i, _, _| i == 2)
    .losing_from(Face::XMin, room(20.0, 1e-6))
    .losing_from(Face::XMax, room(20.0, 1e-6));
    assert!(
        (0..5).all(|i| b.is_void(i, 0, 0) == (i == 2)),
        "two bodies of two cells, the clearance between"
    );
    // The pair across it, without which the two are bodies alone, each at its room, and the step
    // they take together is never taken: the answer is the same, so this is what says it is that.
    let patches = b.gap_patches();
    assert!(
        patches.len() == 1 && patches[0].pairs == 1,
        "one pair across the clearance: {} patches",
        patches.len()
    );
    let why = "two bodies across a clearance from 0 K, each with a film from a room at 20 K";
    let worst = answered(why, &b, 5, 1, 20.0);
    println!("  {why}: {worst:.2e} K from it");
}
