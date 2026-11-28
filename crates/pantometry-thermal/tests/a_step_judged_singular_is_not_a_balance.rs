//! `Solid3D::steady_state`'s step on the levels of the bodies a clearance joins, **judged
//! singular**, and answered all the same: the part's own balance stands in for one of its bodies'
//! rows, and a body nothing moves steps alone.
//!
//! That step is Newton on the bodies' levels together, and the solve calls a pivot singular below
//! 1e-14 of the largest entry. Judged singular, the step was not taken and asked for nothing, and
//! a round the sweeps moved by less than a microkelvin was a balance. It is singular two ways.
//!
//! **A body nothing moves.** One whose every slope is a pair's or a radiating film's with no
//! convection, `4C·|T|³`, has next to none near absolute zero, so beside a warmer one its column
//! is nothing. Measured before the change, of the cases here:
//!
//! - two bodies of two cells across a clearance, a film on the first and the second reached only
//!   through the clearance, were answered within 2.1e-12 K of zero in a room at 1 K from 1e-40 K
//!   and below; in a room at 20 K from 1e-12 K and below, the first at 3.27e-7 K and the second
//!   within 6.6e-40 K of where it started; and in a room at 300 K they did not converge. A second
//!   body of one cell, the film on the second, and three bodies in a row, the same from 0 K;
//! - the first at its room of 20 K and the second at 0 K, the second answered at 2.97e-7 K, and the
//!   other way round, the first at 4.75e-7 K and 6.24e-7 K; three bodies, the last at 0 K beyond
//!   two at the room, the last at 2.97e-7 K; and the last of emissivity 1e-15, so that it is next
//!   to nothing to its pair until it is 2.6 times as hot, at 3.6e-22 K and 1.8e-17 K, and from
//!   starts with the two apart not converging, in either order;
//! - a nanowatt generated in the second body in a room at 0 K, answered with it at 4e-8 K at most,
//!   where it balances at 14.24 K; a milliwatt in a second body of one cell in a room at 20 K, at
//!   0 K, where it balances at 450.43 K; and in a room at 300 K, not converging. A picowatt under a
//!   film of 10 W/m²·K from a room at 0 K, at 4e-11 K where it balances at 2.15 K; 4e-25 W in a
//!   room at 0 K at 1.6e-23 K where it balances at 2.0e-3 K; 1e-20 W in a first body of
//!   emissivity 1e-15, the film on the last of three, answered 2.39 K to 303 K below its balance;
//!   and a milliwatt in a second body of one cell from 1e-12 K, beside a first at its room of
//!   1e-6 K under a film of 10 W/m²·K, at 1.6e-12 K where it balances at 383.61 K.
//!
//! **Rows that cancel.** A film so weak that the part's level has a slope below 1e-14 of its
//! pairs': no body's column is nothing, but each body's row reads its pair's slope, which a shift
//! of the whole part does not move, and the rows add to next to nothing. Two bodies under one were
//! answered where they started, 280 K and 300 K from their room; a microwatt in at the first cell
//! and out at the last, and a third body of low emissivity beside a second at 600 K, did not
//! converge; and two bodies from 0 K in a room at 1e-3 K or 1e-4 K were answered within 2.1e-24 K
//! of zero, their films' slope read as nothing below the rounding of the room.
//!
//! **Four forms of this change were measured wrong first**, and each case that showed it is here.
//! One took out every body whose column held nothing above 1e-14 and refused rows that cancel: it
//! refused the weak film at its own balance, which was answered, and a body whose column was just
//! above 1e-14 and whose pivot below — the second at 2.8e-5 K beside a first at its room of 1 K,
//! and six more. The next read the column the solve failed at, took a body out at 1e-10 of the
//! largest entry, and summed the part's balance from its cells' fluxes: the microwatt through the
//! weak film was answered 4.0e-3 K to 100 K wrong, its film rounded away in its cell's conduction,
//! and a third body of emissivity 1e-11 to 3e-14 was taken out though it was not nothing, and did
//! not converge or was answered 2.1e-5 K and 6.7e-5 K out. The third put the part's balance first,
//! read off what crosses its boundary, always in the first body's row, and took out a body where
//! the solve of that failed: a body held only by a weak pair, whose column only its own row sees,
//! was taken out when its row was the one replaced, and four layouts that were refused, or answered
//! a kelvin out, were answered up to 3624 K out; and a body taken out with nothing to move it by
//! asked for nothing. The fourth put the balance in the row of the body whose own slope is the
//! largest, and read a body taken out with no slope against an exact zero: a picowatt in a room at
//! 0 K did not converge, and a body at a balance of 0 K was thrown from it by what rounding left
//! of its total, a nanowatt started at its closed form refused as below zero — each answered
//! before. Now a body taken out is balanced within 1e-12 of the largest of the part's flows.
//!
//! **A part that loses nothing**, judged singular, was moved as a whole by the heat it already
//! held and asked for nothing: each balance of its step was over its own body's slope, next to
//! nothing near absolute zero, and two sealed bodies from 20 K and 0 K were answered at 20 K and
//! 2.6e-7 K where both balance at 10 K. Now the heat it is to put back stands in for the balance
//! of the body whose own slope is the largest, every row is brought to one, and a body nothing
//! moves is taken out as in a part that loses. A first form kept the slack and brought every
//! balance over the part's largest slope or coupling alike, and a weak first body made it take out
//! a body that moves; those cases are here too.

use pantometry_core::substance::ThermalProps;
use pantometry_core::units::{
    Area, Density, Length, SpecificHeat, Temperature, ThermalConductivity, ThermalExpansion,
    STEFAN_BOLTZMANN,
};
use pantometry_core::Substance;
use pantometry_thermal::{Environment, Face, Solid3D};

/// How far a cell may be from its closed form, K, as the crate's other steady-state tests hold
/// one: the stop at a microkelvin a round leaves `1e-6·ρ/(1 − ρ)` for a round contracting by `ρ`.
/// **Borrowed, and not earned here**: measured, 2.3e-6 K at most, four bodies under a film of
/// 10 W/m²·K from a room at 0 K, and 7.5e-9 K elsewhere; what this file guards against is a miss
/// of 1e-4 K or more, or a refusal.
const SHAPE: f64 = 1e-5;

const DX: f64 = 1e-3;
/// A cell's face, m², through which it conducts to its neighbour and a clearance's pair crosses.
const FACE: f64 = DX * DX;
/// The film's area, m², over the one cell it is on.
const AREA: f64 = 1e-6;
const K: f64 = 50.0;
const EMISSIVITY: f64 = 0.9;
/// A cell's face, its half cell, and a clearance's coefficient between two faces of `EMISSIVITY`.
const WHOLE: f64 = K * FACE / DX;
const HALF: f64 = 2.0 * K * FACE / DX;

fn stuff(emissivity: f64) -> Substance {
    Substance::bulk("stuff", Density::from_si(1000.0)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(K),
        specific_heat: SpecificHeat::j_per_kg_k(4000.0),
        expansion: ThermalExpansion::ppm_per_k(50.0),
        emissivity,
    })
}

/// A film of convection `h` from a room at `t` K, over `area` m².
fn film(t: f64, h: f64, area: f64) -> Environment {
    Environment {
        ambient: Temperature::from_si(t),
        convection_w_per_m2_k: h,
        area: Area::from_si(area),
    }
}

/// `y` with `y·|y|³ = x`. Square roots only, which round correctly.
fn odd_root(x: f64) -> f64 {
    if x < 0.0 {
        -(-x).sqrt().sqrt()
    } else {
        x.sqrt().sqrt()
    }
}

fn fourth(t: f64) -> f64 {
    (t * t) * (t * t)
}

/// A clearance's coefficient between faces of emissivities `a` and `b`, W/K⁴: `σA/(1/a + 1/b − 1)`.
fn pair(a: f64, b: f64) -> f64 {
    STEFAN_BOLTZMANN.to_si() * FACE / (1.0 / a + 1.0 / b - 1.0)
}

/// A row of 1 mm cells, `k = 50`, `ε = 0.9`, from `start`: the bodies of `cells` cells each, in
/// order, an empty cell between each two, so that a clearance joins each to the next; `on` on
/// `face`. A temperature for every cell, and which are empty.
fn row_on(bodies: &[usize], face: Face, on: Environment, start: f64) -> (Solid3D, Vec<bool>) {
    let mut void = Vec::new();
    for (n, cells) in bodies.iter().enumerate() {
        if n > 0 {
            void.push(true);
        }
        void.resize(void.len() + cells, false);
    }
    let empty = void.clone();
    let b = Solid3D::new(
        "row",
        stuff(EMISSIVITY),
        (void.len(), 1, 1),
        Length::from_si(DX),
        Temperature::from_si(start),
    )
    .empty(move |i, _, _| empty[i])
    .losing_from(face, on);
    (b, void)
}

/// [`row_on`] the first body's XMin face.
fn row_under(bodies: &[usize], on: Environment, start: f64) -> (Solid3D, Vec<bool>) {
    row_on(bodies, Face::XMin, on, start)
}

/// [`row_under`] a film with no convection from a room at `at`.
fn row(bodies: &[usize], at: f64, start: f64) -> (Solid3D, Vec<bool>) {
    row_under(bodies, film(at, 0.0, AREA), start)
}

/// Three bodies of two cells, `on` on `face`, the cells of body `n` from `starts[n]` and of
/// emissivity `emissivity[n]`.
fn three(face: Face, on: Environment, starts: [f64; 3], emissivity: [f64; 3]) -> Solid3D {
    let (mut b, _) = row_on(&[2, 2, 2], face, on, starts[0]);
    for n in 0..3 {
        if emissivity[n] != EMISSIVITY {
            b.fill(stuff(emissivity[n]), move |i, _, _| {
                i == 3 * n || i == 3 * n + 1
            });
        }
        for i in [3 * n, 3 * n + 1] {
            b.set_temperature(i, 0, 0, Temperature::from_si(starts[n]));
        }
    }
    b
}

const THREE: [bool; 8] = [false, false, true, false, false, true, false, false];

/// Bodies of `sizes` cells in a row, an empty cell between each two, body `n` of emissivity
/// `emissivity[n]` from `starts[n]`, `on` on `face`; and which cells are empty.
fn bodies(
    sizes: &[usize],
    emissivity: &[f64],
    face: Face,
    on: Environment,
    starts: &[f64],
) -> (Solid3D, Vec<bool>) {
    let (mut b, void) = row_on(sizes, face, on, starts[0]);
    let mut owner = Vec::new();
    for (n, cells) in sizes.iter().enumerate() {
        if n > 0 {
            owner.push(usize::MAX);
        }
        owner.resize(owner.len() + cells, n);
    }
    for n in 0..sizes.len() {
        if emissivity[n] != EMISSIVITY {
            let of = owner.clone();
            b.fill(stuff(emissivity[n]), move |i, _, _| of[i] == n);
        }
        for (i, o) in owner.iter().enumerate() {
            if *o == n {
                b.set_temperature(i, 0, 0, Temperature::from_si(starts[n]));
            }
        }
    }
    (b, void)
}

/// The cells of `row_under(&[2, nb], film(at, h, AREA), _)` with `p` W generated in its last cell,
/// K, the empty one `NaN`: all of it through the second body, across the clearance, through the
/// first and out through the film. Each face carries `p`: a cell's face `kA/dx`, half of one
/// `2kA/dx`, of the face's area or the film's; the clearance's coefficient `σA/(1/ε + 1/ε − 1)`;
/// and the film's surface `s` where `hA(s − T∞) + εσA(s⁴ − T∞⁴) = p` — by Newton from above
/// where there is convection, on a balance that is convex, and otherwise its fourth root.
fn chain(nb: usize, at: f64, h: f64, p: f64) -> Vec<f64> {
    let radiative = EMISSIVITY * STEFAN_BOLTZMANN.to_si() * AREA;
    let surface = if h == 0.0 {
        odd_root(fourth(at) + p / radiative)
    } else {
        let mut s = at + p / (h * AREA);
        for _ in 0..200 {
            let excess = h * AREA * (s - at) + radiative * (fourth(s) - fourth(at)) - p;
            let next = s - excess / (h * AREA + 4.0 * radiative * s * s * s);
            // `NaN` stops it too, which a bare `next < s` written as its negation would say.
            if next.is_nan() || next >= s {
                break;
            }
            s = next;
        }
        s
    };
    let mut t = vec![f64::NAN; 3 + nb];
    t[0] = surface + p / (2.0 * K * AREA / DX);
    t[1] = t[0] + p / WHOLE;
    let near = t[1] + p / HALF;
    let far = odd_root(fourth(near) + p / pair(EMISSIVITY, EMISSIVITY));
    t[3] = far + p / HALF;
    for i in 4..3 + nb {
        t[i] = t[i - 1] + p / WHOLE;
    }
    t
}

/// `b` is answered with every solid cell at `want` to `SHAPE`, `want` being `NaN` at the empty
/// cells `void` names and nowhere else; the worst miss. Counted against the block itself, so that
/// a closed form or a layout gone wrong checks less without a word.
fn answered(why: &str, b: &Solid3D, void: &[bool], want: &[f64]) -> f64 {
    let t = b
        .steady_state()
        .unwrap_or_else(|e| panic!("{why}: there is a balance to strike: {e}"));
    assert_eq!(t.len(), want.len(), "{why}: a closed form for every cell");
    let mut worst = 0.0_f64;
    for (c, (got, w)) in t.iter().zip(want).enumerate() {
        assert_eq!(b.is_void(c, 0, 0), void[c], "{why}: cell {c}'s layout");
        assert_eq!(w.is_nan(), void[c], "{why}: cell {c}'s closed form is {w}");
        if void[c] {
            continue;
        }
        let miss = got - w;
        assert!(
            miss.abs() <= SHAPE,
            "{why}: cell {c} at {got:e} K, {miss:+.3e} K from its balance {w} K"
        );
        worst = worst.max(miss.abs());
    }
    worst
}

/// What the room alone sets: every solid cell at `at`.
fn at_room(void: &[bool], at: f64) -> Vec<f64> {
    void.iter()
        .map(|v| if *v { f64::NAN } else { at })
        .collect()
}

/// A layout for [`sealed`]: each body's cells, its emissivity, how many times as dense as
/// [`stuff`] it is, and its start, K.
type Layout = (
    &'static [usize],
    &'static [f64],
    &'static [f64],
    &'static [f64],
);

/// A layout for [`drawn`]: each body's cells, its emissivity and how many times as dense it is,
/// the watts generated in the first cell and drawn from the last, and each body's start, K.
type Drawn = (
    &'static [usize],
    &'static [f64],
    &'static [f64],
    f64,
    &'static [f64],
);

/// Bodies of `sizes` cells in a row, an empty cell between each two and **no film anywhere**, so a
/// part that loses nothing: body `n` of emissivity `emissivity[n]`, `denser[n]` times as dense as
/// [`stuff`] and so of that many times its capacity, from `starts[n]`; and which cells are empty.
fn sealed(
    sizes: &[usize],
    emissivity: &[f64],
    denser: &[f64],
    starts: &[f64],
) -> (Solid3D, Vec<bool>) {
    let mut owner = Vec::new();
    for (n, cells) in sizes.iter().enumerate() {
        if n > 0 {
            owner.push(usize::MAX);
        }
        owner.resize(owner.len() + cells, n);
    }
    let void: Vec<bool> = owner.iter().map(|o| *o == usize::MAX).collect();
    let empty = void.clone();
    let mut b = Solid3D::new(
        "sealed",
        stuff(EMISSIVITY),
        (void.len(), 1, 1),
        Length::from_si(DX),
        Temperature::from_si(starts[0]),
    )
    .empty(move |i, _, _| empty[i]);
    for n in 0..sizes.len() {
        if emissivity[n] != EMISSIVITY || denser[n] != 1.0 {
            let dense = Substance::bulk("dense", Density::from_si(1000.0 * denser[n]))
                .with_thermal(stuff(emissivity[n]).thermal.expect("stuff is thermal"));
            let of = owner.clone();
            b.fill(dense, move |i, _, _| of[i] == n);
        }
        for (i, o) in owner.iter().enumerate() {
            if *o == n {
                b.set_temperature(i, 0, 0, Temperature::from_si(starts[n]));
            }
        }
    }
    (b, void)
}

/// **Bodies a clearance joins are answered from absolute zero**, every cell at its room: two of
/// two cells from 0 K and from starts near it, in rooms at 1 K, 20 K and 300 K; a second body of
/// one cell; three bodies in a row; and, from the other end, the film on the second alone.
#[test]
fn bodies_a_clearance_joins_are_answered_from_absolute_zero() {
    for at in [1.0, 20.0, 300.0] {
        for start in [0.0, 1e-300, 1e-40, 1e-12] {
            let (b, void) = row(&[2, 2], at, start);
            let why = format!("two bodies from {start:e} K, a room at {at} K");
            let worst = answered(&why, &b, &void, &at_room(&void, at));
            println!("  {why}: {worst:.2e} K from it");
        }
    }
    for at in [20.0, 300.0] {
        for (bodies, what) in [
            (&[2, 1][..], "a second body of one cell"),
            (&[2, 2, 2], "three"),
        ] {
            let (b, void) = row(bodies, at, 0.0);
            let why = format!("{what} from 0 K, a room at {at} K");
            let worst = answered(&why, &b, &void, &at_room(&void, at));
            println!("  {why}: {worst:.2e} K from it");
        }
        let (b, void) = row_on(&[2, 2], Face::XMax, film(at, 0.0, AREA), 0.0);
        let why = format!("the film on the second, from 0 K, a room at {at} K");
        let worst = answered(&why, &b, &void, &at_room(&void, at));
        println!("  {why}: {worst:.2e} K from it");
    }
}

/// **A body at absolute zero beside a warm one is answered**, every cell at its room: in a room at
/// 20 K, the first at the room and the second at 0 K, and the other way round; three bodies, the
/// first two at the room and the third at 0 K beyond them, so that two bodies with slopes are left
/// to step together; and the third of emissivity 1e-15, next to nothing to its pair until it is
/// 2.6 times as hot, from starts with the other two apart — and the same in mirror order, the film
/// on the last body and the first of emissivity 1e-15.
#[test]
fn a_body_at_zero_beside_a_warm_one_is_answered() {
    for (first, rest) in [(20.0, 0.0), (0.0, 20.0)] {
        let (mut b, void) = row(&[2, 2], 20.0, first);
        for i in 3..5 {
            b.set_temperature(i, 0, 0, Temperature::from_si(rest));
        }
        let why = format!("the first at {first} K and the second at {rest} K");
        let worst = answered(&why, &b, &void, &at_room(&void, 20.0));
        println!("  {why}: {worst:.2e} K from the room");
    }
    let b = three(
        Face::XMin,
        film(20.0, 0.0, AREA),
        [20.0, 20.0, 0.0],
        [EMISSIVITY; 3],
    );
    let why = "three bodies, the last at 0 K";
    let worst = answered(why, &b, &THREE, &at_room(&THREE, 20.0));
    println!("  {why}: {worst:.2e} K from the room");
    for at in [20.0, 300.0] {
        for (first, second) in [(300.0, 600.0), (600.0, 300.0), (at, at), (20.0, 1000.0)] {
            let on = film(at, 0.0, AREA);
            for (face, starts, emissivity, order) in [
                (
                    Face::XMin,
                    [first, second, 0.0],
                    [EMISSIVITY, EMISSIVITY, 1e-15],
                    "the last",
                ),
                (
                    Face::XMax,
                    [0.0, second, first],
                    [1e-15, EMISSIVITY, EMISSIVITY],
                    "the first",
                ),
            ] {
                let b = three(face, on, starts, emissivity);
                let why = format!(
                    "{order} of three of emissivity 1e-15 at 0 K, the others at {first} K and \
                     {second} K, a room at {at} K"
                );
                let worst = answered(&why, &b, &THREE, &at_room(&THREE, at));
                println!("  {why}: {worst:.2e} K from the room");
            }
        }
    }
}

/// **Heat generated beyond a clearance is answered from absolute zero**, at its chain: a nanowatt
/// and a microwatt in the second body in a room at 0 K, a nanowatt in a room at 300 K, a microwatt
/// and a milliwatt in a second body of one cell in rooms at 20 K and 300 K; a picowatt under a film
/// of 10 W/m²·K from a room at 0 K, where the first body's own step asks for little and the
/// second's for a great deal; 4e-25 W in a room at 0 K, where the first's slope is nothing at the
/// start; and a milliwatt in a second body of one cell from 1e-12 K, the first at its room of 1e-6 K
/// or 1e-7 K under a film of 10 W/m²·K, a block whose level is within a microkelvin of zero, so
/// that a step held to a quarter of it moves less than the stop reads. That was answered at
/// 1.6e-12 K, where it balances at 383.61 K.
#[test]
fn heat_generated_beyond_a_clearance_is_answered_from_absolute_zero() {
    let cases = [
        (2, 1e-9, 0.0, 0.0),
        (2, 1e-6, 0.0, 0.0),
        (2, 1e-9, 300.0, 0.0),
        (1, 1e-6, 20.0, 0.0),
        (1, 1e-3, 20.0, 0.0),
        (1, 1e-6, 300.0, 0.0),
        (1, 1e-3, 300.0, 0.0),
        (2, 1e-12, 0.0, 10.0),
        (2, 4e-25, 0.0, 0.0),
    ];
    for (nb, p, at, h) in cases {
        let (b, void) = row_under(&[2, nb], film(at, h, AREA), 0.0);
        let last = void.len() - 1;
        let b = b.dissipating(p, move |i, _, _| i == last);
        let want = chain(nb, at, h, p);
        let why = format!("{p:e} W in a second body of {nb}, a film of {h} W/m2K from {at} K");
        let worst = answered(&why, &b, &void, &want);
        println!(
            "  {why}: its far cell at {:e} K, {worst:.2e} K from it",
            want[last]
        );
    }
    // A block within a microkelvin of zero: the first at its room, the second from 1e-12 K.
    for at in [1e-6, 1e-7] {
        let (mut b, void) = row_under(&[2, 1], film(at, 10.0, AREA), at);
        b.set_temperature(3, 0, 0, Temperature::from_si(1e-12));
        let b = b.dissipating(1e-3, |i, _, _| i == 3);
        let want = chain(1, at, 10.0, 1e-3);
        let why =
            format!("1e-3 W in a second body of 1 from 1e-12 K, the first at its room of {at:e} K");
        let worst = answered(&why, &b, &void, &want);
        println!(
            "  {why}: its cell at {:e} K, {worst:.2e} K from it",
            want[3]
        );
    }
}

/// **Heat generated in a body nothing moves is answered**: three bodies, the first of emissivity
/// 1e-15 generating 1e-20 W in its first cell, a film of 1e-21 m² from a room at 300 K on the
/// last — so that the first's column is all but nothing beside the others' and their rows cancel.
/// Every cell at the chain from the film back: the surface `s` where `εσA(s⁴ − T∞⁴) = p`, each
/// half cell and face carrying `p`, and the two clearances' `C(u⁴ − v⁴) = p`, the first's of
/// emissivities 1e-15 and 0.9. Between the film's half cell and the first clearance the drops are
/// below what an `f64` resolves at 301.8 K, so what this sees is the film's rise of 1.80 K and the
/// first body's of 1.59 K across its pair. It was answered 2.39 K to 303 K below that.
#[test]
fn heat_generated_in_a_body_nothing_moves_is_answered() {
    let (p, area) = (1e-20, 1e-21);
    let s = odd_root(fourth(300.0) + p / (EMISSIVITY * STEFAN_BOLTZMANN.to_si() * area));
    let t7 = s + p / (2.0 * K * area / DX);
    let t6 = t7 + p / WHOLE;
    let t4 = odd_root(fourth(t6 + p / HALF) + p / pair(EMISSIVITY, EMISSIVITY)) + p / HALF;
    let t3 = t4 + p / WHOLE;
    let t1 = odd_root(fourth(t3 + p / HALF) + p / pair(1e-15, EMISSIVITY)) + p / HALF;
    let want = [t1 + p / WHOLE, t1, f64::NAN, t3, t4, f64::NAN, t6, t7];
    for start in [300.0, 301.0, 0.0] {
        let b = three(
            Face::XMax,
            film(300.0, 0.0, area),
            [start; 3],
            [1e-15, EMISSIVITY, EMISSIVITY],
        )
        .dissipating(p, |i, _, _| i == 0);
        let why = format!("1e-20 W in a first body of emissivity 1e-15, from {start} K");
        let worst = answered(&why, &b, &THREE, &want);
        println!("  {why}: at {:.6} K, {worst:.2e} K from it", want[0]);
    }
}

/// **A body whose pivot is nothing is answered**: two bodies of two cells, a film on the first
/// from a room of 1 K to 20 K, every cell at the room but the second's, a few tens of microkelvin
/// to half a millikelvin, so that the second's column is above 1e-14 of the largest entry and its
/// pivot, once the first's row is taken out of its, below. A first form of this change, reading
/// the columns alone, refused the first seven. Every cell at the room.
#[test]
fn a_body_whose_pivot_is_nothing_is_answered() {
    for (at, second) in [
        (1.0, 2.8e-5),
        (1.0, 3.0e-5),
        (1.0, 3.3e-5),
        (5.0, 1.4e-4),
        (5.0, 1.6e-4),
        (10.0, 2.8e-4),
        (10.0, 3.0e-4),
        (15.0, 4.2e-4),
        (20.0, 5.6e-4),
    ] {
        let (mut b, void) = row(&[2, 2], at, at);
        for i in 3..5 {
            b.set_temperature(i, 0, 0, Temperature::from_si(second));
        }
        let why = format!("the second at {second:e} K, a room at {at} K");
        let worst = answered(&why, &b, &void, &at_room(&void, at));
        println!("  {why}: {worst:.2e} K from the room");
    }
}

/// **A part whose rows cancel is answered**: two bodies of two cells across a clearance of 1 mm²,
/// under a film of 1e-21 m², 1e-22 m² or 1e-23 m², whose slope is 1.1e-15 to 1.1e-17 of the
/// pair's, from starts 280 K and 300 K from the room, and from the room itself and 1e-7 K above it.
/// Every cell at the room. The last two were answered before, the step being zero, and guard only
/// against a refusal — which a first form of this change gave them.
#[test]
fn a_part_whose_rows_cancel_is_answered() {
    for area in [1e-21, 1e-22, 1e-23] {
        for (at, start) in [
            (20.0, 300.0),
            (300.0, 20.0),
            (300.0, 600.0),
            (20.0, 20.0),
            (300.0, 300.0),
            (300.0, 300.0 + 1e-7),
        ] {
            let (b, void) = row_under(&[2, 2], film(at, 0.0, area), start);
            let why = format!("a film of {area:e} m2 from a room at {at} K, from {start} K");
            let worst = answered(&why, &b, &void, &at_room(&void, at));
            println!("  {why}: {worst:.2e} K from the room");
        }
    }
}

/// **Heat through a part whose rows cancel is answered**: a microwatt in at the first of two
/// bodies' cells and out at the last, under a film of 1e-21 m² to 1e-25 m² from a room at 300 K,
/// so that the film carries nothing and the first cell sits at the room; each face carries the
/// microwatt, and the clearance `C(u⁴ − v⁴)` of it. Where the part's balance was summed from its
/// cells' fluxes the film was rounded away in its cell's conduction.
#[test]
fn heat_through_a_part_whose_rows_cancel_is_answered() {
    let p = 1e-6;
    let t1 = 300.0 - p / WHOLE;
    let t3 = odd_root(fourth(t1 - p / HALF) - p / pair(EMISSIVITY, EMISSIVITY)) - p / HALF;
    let want = [300.0, t1, f64::NAN, t3, t3 - p / WHOLE];
    for area in [1e-21, 1e-23, 1e-25] {
        for start in [320.0, 280.0, 400.0] {
            let (b, void) = row_under(&[2, 2], film(300.0, 0.0, area), start);
            let b = b
                .dissipating(p, |i, _, _| i == 0)
                .dissipating(-p, |i, _, _| i == 4);
            let why = format!("a microwatt through, a film of {area:e} m2, from {start} K");
            let worst = answered(&why, &b, &void, &want);
            println!("  {why}: {worst:.2e} K from it");
        }
    }
}

/// **A warm body beside rows that cancel is answered**: three bodies under a film of 1e-21 m² from
/// a room at 300 K, the middle one at 600 K, one end body of emissivity 1e-11 to 1e-15 — the last,
/// the film on the first, or the first, the film on the last. Every cell at the room. A form of
/// this change that took a body out at 1e-10 of the largest entry took the end body out though it
/// was not nothing.
#[test]
fn a_warm_body_beside_rows_that_cancel_is_answered() {
    for (e, face) in [
        (1e-11, Face::XMin),
        (1e-13, Face::XMin),
        (3e-14, Face::XMin),
        (1e-11, Face::XMax),
        (1e-13, Face::XMax),
        (1e-15, Face::XMax),
    ] {
        let emissivity = if face == Face::XMin {
            [EMISSIVITY, EMISSIVITY, e]
        } else {
            [e, EMISSIVITY, EMISSIVITY]
        };
        let b = three(
            face,
            film(300.0, 0.0, 1e-21),
            [300.0, 600.0, 300.0],
            emissivity,
        );
        let why = format!("an end body of emissivity {e:e}, the film on {face:?}");
        let worst = answered(&why, &b, &THREE, &at_room(&THREE, 300.0));
        println!("  {why}: {worst:.2e} K from the room");
    }
}

/// **A body held only by a weak pair is answered**, every cell at its room: a body of emissivity
/// 1e-15 among others under a film too weak to show, in four layouts — four bodies, the second of
/// it and the last a kelvin warm under films of 1e-21 m² to 1e-25 m² from a room at 20 K; four,
/// the second a single cell, under a film of 1e-18 m² from 300 K; four, the last under a film of
/// 10 W/m²·K and 1e-23 m² from a room at 0 K; and three, the first of it, the middle at 600 K,
/// under films of 1e-21 m² to 1e-23 m² from 300 K. The part's balance stood in for the first
/// body's row, which for such a body is the only one its column is seen in, and they were answered
/// up to 3624 K out, or 1.6e-5 K and 7.4e-5 K, or did not converge.
#[test]
fn a_body_held_only_by_a_weak_pair_is_answered() {
    let weak = 1e-15;
    let mut cases = Vec::new();
    for area in [1e-21, 1e-23, 1e-25] {
        cases.push((
            vec![2, 2, 2, 2],
            vec![EMISSIVITY, weak, EMISSIVITY, EMISSIVITY],
            Face::XMax,
            film(20.0, 0.0, area),
            vec![20.0, 20.0, 20.0, 21.0],
            20.0,
        ));
    }
    cases.push((
        vec![1, 1, 2, 2],
        vec![EMISSIVITY, weak, EMISSIVITY, EMISSIVITY],
        Face::XMin,
        film(300.0, 0.0, 1e-18),
        vec![300.0, 601.0, 1e-6, 601.0],
        300.0,
    ));
    cases.push((
        vec![3, 2, 3, 2],
        vec![0.5, weak, 0.5, 0.5],
        Face::XMax,
        film(0.0, 10.0, 1e-23),
        vec![1.0, 300.0, 300.0, 0.1],
        0.0,
    ));
    for area in [1e-21, 1e-22, 3e-23, 1e-23] {
        cases.push((
            vec![2, 2, 2],
            vec![weak, EMISSIVITY, EMISSIVITY],
            Face::XMax,
            film(300.0, 0.0, area),
            vec![300.0, 600.0, 300.0],
            300.0,
        ));
    }
    for (sizes, emissivity, face, on, starts, at) in cases {
        let why = format!(
            "bodies {sizes:?} of emissivity {emissivity:?} from {starts:?} K, a film of {:e} m2 \
             on {face:?} from {at} K",
            on.area.to_si()
        );
        let (b, void) = bodies(&sizes, &emissivity, face, on, &starts);
        let worst = answered(&why, &b, &void, &at_room(&void, at));
        println!("  {why}: {worst:.2e} K from the room");
    }
}

/// **A body taken out with nothing to move it by has not settled**: two bodies from 0 K, the second
/// of emissivity 1e-15 under a film of 1e-18 m² from a room at 1e-3 K; and three from near it, of
/// emissivities 1e-15, 1e-11 and 1e-15, under one of 1e-21 m² from a room at 1 K. Every body was
/// taken out of the step with no slope to step on, and a round that moved nothing else stopped, at
/// 0 K and within a microkelvin of it. Every cell at the room.
#[test]
fn a_body_with_nothing_to_move_it_has_not_settled() {
    for (sizes, emissivity, on, starts, at) in [
        (
            vec![3, 3],
            vec![EMISSIVITY, 1e-15],
            film(1e-3, 0.0, 1e-18),
            vec![0.0, 0.0],
            1e-3,
        ),
        (
            vec![1, 3, 3],
            vec![1e-15, 1e-11, 1e-15],
            film(1.0, 0.0, 1e-21),
            vec![1e-6, 0.0, 0.0],
            1.0,
        ),
    ] {
        let why = format!("bodies {sizes:?} of emissivity {emissivity:?}, a room at {at} K");
        let (b, void) = bodies(&sizes, &emissivity, Face::XMax, on, &starts);
        let worst = answered(&why, &b, &void, &at_room(&void, at));
        println!("  {why}: {worst:.2e} K from the room");
    }
}

/// **A balance found is kept**, in a room at 0 K, where a body taken out of the step sits at a
/// balance of 0 K and has next to no slope: a nanowatt drawn from a body across a clearance of
/// emissivity 1e-11, started at its closed form; a microwatt across clearances of emissivity 1e-8
/// and 1e-13 into a sink under a film, from 1 K and again from what that answered; and a picowatt
/// shared by two cells of a body whose neighbour is an ulp from zero, from 0 K. What rounding left
/// of a balanced body's total, over its slope, threw it a quarter of the block's level, and the
/// first two were refused as below zero; read against nothing, the last did not converge. Each was
/// answered before. And the nanowatt again with the film on the other body, of 1e-6 m² to
/// 1e-14 m², from a room at its last cell's balance: what is left of the part once the cold body
/// is taken out balances only with what crosses to it, and without that was refused as below zero.
/// Every cell at its chain.
#[test]
fn a_balance_found_is_kept() {
    // A nanowatt from the second body's last cell into the first's last, across the clearance;
    // the film carries nothing, so the first sits at the room.
    let p = 1e-9;
    let t3 = odd_root(fourth(p / HALF) + p / pair(1e-11, EMISSIVITY)) + p / HALF;
    let want = [0.0, 0.0, f64::NAN, t3, t3 + p / WHOLE, t3 + 2.0 * p / WHOLE];
    let (mut b, void) = bodies(
        &[2, 3],
        &[1e-11, EMISSIVITY],
        Face::XMin,
        film(0.0, 0.0, 1e-18),
        &[0.0, t3],
    );
    for (i, w) in want.iter().enumerate() {
        if !void[i] {
            b.set_temperature(i, 0, 0, Temperature::from_si(*w));
        }
    }
    let b = b
        .dissipating(-p, |i, _, _| i == 1)
        .dissipating(p, |i, _, _| i == 5);
    let why = "a nanowatt across emissivity 1e-11, from its closed form";
    let worst = answered(why, &b, &void, &want);
    println!("  {why}: {worst:.2e} K from it");
    // A microwatt from the first body's first cell into the second's sink, the film on it
    // carrying nothing; from 1 K, and again from that answer.
    let p = 1e-6;
    let t2 = odd_root(fourth(p / HALF) + p / pair(1e-8, 1e-13)) + p / HALF;
    let want = [t2 + 2.0 * p / WHOLE, t2 + p / WHOLE, t2, f64::NAN, 0.0];
    let (b, void) = bodies(
        &[3, 1],
        &[1e-8, 1e-13],
        Face::XMax,
        film(0.0, 0.0, 1e-12),
        &[1.0, 1.0],
    );
    let mut b = b
        .dissipating(p, |i, _, _| i == 0)
        .dissipating(-p, |i, _, _| i == 4);
    let why = "a microwatt across emissivities 1e-8 and 1e-13, from 1 K";
    let worst = answered(why, &b, &void, &want);
    println!("  {why}: {worst:.2e} K from it");
    let first = b.steady_state().expect("answered just now");
    for (i, v) in first.iter().enumerate() {
        if !void[i] {
            b.set_temperature(i, 0, 0, Temperature::from_si(*v));
        }
    }
    let why = "the same, from what it was answered";
    let worst = answered(why, &b, &void, &want);
    println!("  {why}: {worst:.2e} K from it");
    // A picowatt and more shared by the first body's two cells, the film on the second.
    for p in [1e-12, 2e-10] {
        let (b, void) = bodies(
            &[2, 2],
            &[EMISSIVITY, EMISSIVITY],
            Face::XMax,
            film(0.0, 0.0, AREA),
            &[0.0, 0.0],
        );
        let b = b
            .dissipating(p, |i, _, _| i == 0)
            .dissipating(-p, |i, _, _| i == 1);
        let why = format!("{p:e} W from one cell to the next, from 0 K");
        let worst = answered(&why, &b, &void, &[p / WHOLE, 0.0, f64::NAN, 0.0, 0.0]);
        println!("  {why}: {worst:.2e} K from it");
    }
    // The nanowatt across emissivity 1e-11 again, the film on the second body's last cell instead,
    // from a room at that cell's balance, so it carries nothing there; from its closed form.
    let p = 1e-9;
    let t3 = odd_root(fourth(p / HALF) + p / pair(1e-11, EMISSIVITY)) + p / HALF;
    let want = [0.0, 0.0, f64::NAN, t3, t3 + p / WHOLE, t3 + 2.0 * p / WHOLE];
    for area in [1e-6, 1e-12, 1e-14] {
        let (mut b, void) = bodies(
            &[2, 3],
            &[1e-11, EMISSIVITY],
            Face::XMax,
            film(want[5], 0.0, area),
            &[0.0, t3],
        );
        for (i, w) in want.iter().enumerate() {
            if !void[i] {
                b.set_temperature(i, 0, 0, Temperature::from_si(*w));
            }
        }
        let b = b
            .dissipating(-p, |i, _, _| i == 1)
            .dissipating(p, |i, _, _| i == 5);
        let why =
            format!("a nanowatt across emissivity 1e-11, a film of {area:e} m2 on the second");
        let worst = answered(&why, &b, &void, &want);
        println!("  {why}: {worst:.2e} K from it");
    }
}

/// **A cold room is answered from absolute zero**: two bodies of two cells from 0 K, a film on the
/// first from a room at 1e-3 K, 3e-4 K or 1e-4 K. A cell below the rounding of its room has its
/// film's surface at exactly zero, `film_surface` starting from the room, so its film reads no
/// slope and the part's rows cancel to nothing, the step being singular with its balance in place;
/// a body is taken out of it, and the two warm past that rounding. They were answered within
/// 2.1e-24 K of zero, and refused by both earlier forms of this change. Every cell at the room.
#[test]
fn a_cold_room_is_answered_from_absolute_zero() {
    for at in [1e-3, 3e-4, 1e-4] {
        let (b, void) = row(&[2, 2], at, 0.0);
        let why = format!("two bodies from 0 K, a room at {at:e} K");
        let worst = answered(&why, &b, &void, &at_room(&void, at));
        println!("  {why}: {worst:.2e} K from it");
    }
}

/// **A part that loses nothing is answered at the level its heat sets**, judged singular: bodies a
/// clearance joins and nothing else reaches, every cell at the heat the starts hold over the
/// part's capacity, `Σ C·T₀ / Σ C`. Each balance of the step they take together was over its own
/// body's slope, and a body near absolute zero has next to none, so its row's entries came to as
/// much as 7e23 beside the heat row's halves, the step was judged singular, and the part was moved
/// as a whole by the heat it already held: nothing, and the step asked for nothing. Measured before
/// the change, of the cases here: two bodies of two cells from 20 K and 0 K answered at 20 K and
/// 2.6e-7 K where both balance at 10 K, and the other way round; from 1 K and 0 K at 1 K and
/// 1.6e-12 K, from 1e-3 K at 1e-3 K and 7.0e-25 K; with the second at 1e-12 K or 1e-6 K, at 20 K
/// and 2.6e-7 K or 1.26e-6 K; the second of emissivity 1e-11, from 300 K at 300 K and 1.7e-13 K
/// where they balance at 150 K; a first of three cells, at 20 K and 8.9e-8 K where they balance at
/// 15 K; three bodies, one or two of them from 20 K, the warm ones left there and the others within
/// 3.3e-7 K of zero; and one of three times the other's density, at 20 K and 3.9e-7 K where they
/// balance at 15 K, or at 5 K with the dense one cold. Fifteen of the eighteen were answered so.
/// Two others, a first body of one cell, its sweep moving it, were answered, and reach the step
/// all the same; and two bodies of one cell at absolute zero, with nothing to move them and
/// nothing to move, were answered at it, as they must be by a body taken out as balanced.
#[test]
fn a_part_that_loses_nothing_is_answered_at_the_level_its_heat_sets() {
    const E: f64 = EMISSIVITY;
    let cases: [Layout; 18] = [
        (&[1, 1], &[E, E], &[1.0, 1.0], &[0.0, 0.0]),
        (&[1, 1], &[E, E], &[1.0, 1.0], &[20.0, 0.0]),
        (&[1, 3], &[E, E], &[1.0, 1.0], &[20.0, 0.0]),
        (&[2, 2], &[E, E], &[1.0, 1.0], &[20.0, 0.0]),
        (&[2, 2], &[E, E], &[1.0, 1.0], &[0.0, 20.0]),
        (&[2, 2], &[E, E], &[1.0, 1.0], &[1.0, 0.0]),
        (&[2, 2], &[E, E], &[1.0, 1.0], &[1e-3, 0.0]),
        (&[2, 2], &[E, E], &[1.0, 1.0], &[20.0, 1e-12]),
        (&[2, 2], &[E, E], &[1.0, 1.0], &[20.0, 1e-6]),
        (&[2, 2], &[E, 1e-11], &[1.0, 1.0], &[300.0, 0.0]),
        (&[2, 2], &[E, 1e-11], &[1.0, 1.0], &[20.0, 0.0]),
        (&[3, 1], &[E, E], &[1.0, 1.0], &[20.0, 0.0]),
        (&[2, 2, 2], &[E, E, E], &[1.0, 1.0, 1.0], &[20.0, 0.0, 0.0]),
        (&[2, 2, 2], &[E, E, E], &[1.0, 1.0, 1.0], &[0.0, 20.0, 0.0]),
        (&[2, 2, 2], &[E, E, E], &[1.0, 1.0, 1.0], &[0.0, 0.0, 20.0]),
        (&[2, 2, 2], &[E, E, E], &[1.0, 1.0, 1.0], &[20.0, 20.0, 0.0]),
        (&[2, 2], &[E, E], &[3.0, 1.0], &[20.0, 0.0]),
        (&[2, 2], &[E, E], &[1.0, 3.0], &[20.0, 0.0]),
    ];
    for (sizes, emissivity, denser, starts) in cases {
        at_the_level_its_heat_sets(sizes, emissivity, denser, starts);
    }
}

/// [`sealed`] is answered with every solid cell at the level the heat its starts hold sets,
/// `Σ C·T₀ / Σ C`, each body's capacity its cells times how much denser it is.
fn at_the_level_its_heat_sets(sizes: &[usize], emissivity: &[f64], denser: &[f64], starts: &[f64]) {
    let (b, void) = sealed(sizes, emissivity, denser, starts);
    let (mut heat, mut capacity) = (0.0, 0.0);
    for n in 0..sizes.len() {
        let c = denser[n] * sizes[n] as f64;
        heat += c * starts[n];
        capacity += c;
    }
    let level = heat / capacity;
    let why = format!(
        "bodies of {sizes:?} cells, emissivities {emissivity:?}, {denser:?} as dense, sealed, \
         from {starts:?} K"
    );
    let worst = answered(&why, &b, &void, &at_room(&void, level));
    println!("  {why}: {worst:.2e} K from {level} K");
}

/// **A weak body in a part that loses nothing is answered**, judged singular: three or four
/// bodies, sealed, one of them of emissivity 1e-15 or 1e-11, its pair next to nothing beside the
/// others', and every cell at the level the heat the starts hold sets. A first form of this change
/// kept the step's slack and brought every balance over the part's largest slope or coupling
/// alike: with the weak body first, its column below 1e-14 of theirs, the solve failed at the last
/// of the strong ones, whose own shift their rows do not see, and it was taken out though it
/// moves. Measured on that form: the six of one cell from 20 K or 300 K, the weak one first, did
/// not converge, and two more from 0 K, 0 K and 20 K or 300 K; the three of two cells, the weak
/// one first, were answered up to 2812 K out; and four of one cell, the first a million times as
/// dense and the weak one last, left the last at 0 K where all balance at 4e-5 K. Measured on the
/// step as it was: those three of two cells refused as not converging, and the two from 0 K, 0 K
/// and the third warm, and two with the weak one last; three of two cells from 20 K, 300 K and
/// 300 K, the weak one first, answered where they started, where all balance at 206.67 K, as the
/// first form answered them too; and four of two cells, one a million times as dense, the weak one
/// of emissivity 1e-11 last, answered 20 K to 300 K out; and two of one, one, one and two cells,
/// the third of emissivity 1e-15 and the last a million times as dense, from a microkelvin, a
/// microkelvin, 0 K and 20 K or 300 K, answered with the first three left within a microkelvin of
/// zero, where all balance at 20.00 K and 300.00 K. Those two need the step a body taken out
/// takes: without it they do not converge.
#[test]
fn a_weak_body_in_a_part_that_loses_nothing_is_answered() {
    const E: f64 = EMISSIVITY;
    const W: f64 = 1e-15;
    // One layout to a line, as a table.
    #[rustfmt::skip]
    let cases: [Layout; 23] = [
        (&[1, 1, 1], &[W, E, E], &[1.0; 3], &[20.0, 0.0, 0.0]),
        (&[1, 1, 1], &[W, E, E], &[1.0; 3], &[300.0, 0.0, 0.0]),
        (&[1, 1, 1], &[W, E, E], &[1.0; 3], &[300.0, 20.0, 20.0]),
        (&[1, 1, 1, 1], &[W, E, E, E], &[1.0; 4], &[20.0, 0.0, 0.0, 0.0]),
        (&[1, 1, 1, 1], &[W, E, E, E], &[1.0; 4], &[300.0, 0.0, 0.0, 0.0]),
        (&[1, 1, 1, 1], &[W, E, E, E], &[1.0; 4], &[300.0, 20.0, 20.0, 20.0]),
        (&[1, 1, 1], &[W, E, E], &[1.0; 3], &[0.0, 0.0, 20.0]),
        (&[1, 1, 1], &[W, E, E], &[1.0; 3], &[0.0, 0.0, 300.0]),
        (&[2, 2, 2], &[W, E, E], &[1.0; 3], &[20.0, 0.0, 300.0]),
        (&[2, 2, 2, 2], &[W, E, E, E], &[1.0; 4], &[0.0, 0.0, 0.0, 300.0]),
        (&[2, 2, 2], &[W, E, E], &[1e6, 1.0, 1.0], &[20.0, 0.0, 300.0]),
        (&[1, 1, 1, 1], &[E, E, E, W], &[1e6, 1.0, 1.0, 1.0], &[0.0, 20.0, 20.0, 0.0]),
        (&[2, 2, 2], &[W, E, E], &[1.0; 3], &[20.0, 300.0, 300.0]),
        (&[1, 1, 1], &[E, E, W], &[1.0; 3], &[0.0, 20.0, 0.0]),
        (&[2, 2, 2], &[E, E, W], &[1.0; 3], &[300.0, 0.0, 20.0]),
        (&[1, 1, 1], &[E, E, W], &[1.0, 1.0, 1e6], &[0.0, 0.0, 20.0]),
        (&[1, 1, 1], &[E, E, W], &[1.0; 3], &[0.0, 0.0, 20.0]),
        (&[1, 1, 1, 1], &[E, E, E, W], &[1.0; 4], &[0.0, 0.0, 0.0, 300.0]),
        (&[2, 2, 2, 2], &[E, E, E, 1e-11], &[1.0, 1.0, 1.0, 1e6], &[0.0, 0.0, 0.0, 300.0]),
        (&[2, 2, 2, 2], &[E, E, E, 1e-11], &[1e6, 1.0, 1.0, 1.0], &[0.0, 1e-6, 0.0, 20.0]),
        (&[2, 2, 2, 2], &[E, E, E, 1e-11], &[1e6, 1.0, 1.0, 1.0], &[0.0, 1e-6, 0.0, 300.0]),
        (&[1, 1, 1, 2], &[E, E, W, E], &[1.0, 1.0, 1.0, 1e6], &[1e-6, 1e-6, 0.0, 20.0]),
        (&[1, 1, 1, 2], &[E, E, W, E], &[1.0, 1.0, 1.0, 1e6], &[1e-6, 1e-6, 0.0, 300.0]),
    ];
    for (sizes, emissivity, denser, starts) in cases {
        at_the_level_its_heat_sets(sizes, emissivity, denser, starts);
    }
}

/// [`sealed`] with `p` W generated in its first cell and drawn from its last.
fn drawn(
    sizes: &[usize],
    emissivity: &[f64],
    denser: &[f64],
    p: f64,
    starts: &[f64],
) -> (Solid3D, Vec<bool>) {
    let (b, void) = sealed(sizes, emissivity, denser, starts);
    let last = void.len() - 1;
    let b = b
        .dissipating(p, |i, _, _| i == 0)
        .dissipating(-p, move |i, _, _| i == last);
    (b, void)
}

/// The cells of [`drawn`] at the level the heat its starts hold sets, K, the empty ones `NaN`: the
/// chain `p` crosses — a cell's face, half of one, the clearance's fourth powers, read odd below
/// zero as the pairs read them, half a cell and a cell's face, to each next body — from the first
/// body's last cell, found by bisection on it, on which every cell rises; each cell's heat its
/// temperature times how much denser its body is.
fn through(
    sizes: &[usize],
    emissivity: &[f64],
    denser: &[f64],
    p: f64,
    starts: &[f64],
) -> Vec<f64> {
    let chain = |x: f64| {
        let first = sizes[0];
        let mut t: Vec<f64> = (0..first)
            .map(|i| x + (first - 1 - i) as f64 * p / WHOLE)
            .collect();
        let mut cell = x;
        for n in 1..sizes.len() {
            let near = cell - p / HALF;
            let odd = if near < 0.0 {
                -fourth(near)
            } else {
                fourth(near)
            };
            let far = odd_root(odd - p / pair(emissivity[n - 1], emissivity[n]));
            cell = far - p / HALF;
            t.extend([f64::NAN, cell]);
            for _ in 1..sizes[n] {
                cell -= p / WHOLE;
                t.push(cell);
            }
        }
        t
    };
    let mut weight = Vec::new();
    for (n, cells) in sizes.iter().enumerate() {
        if n > 0 {
            weight.push(0.0);
        }
        weight.resize(weight.len() + cells, denser[n]);
    }
    let held = |t: &[f64]| {
        t.iter()
            .zip(&weight)
            .filter(|(v, _)| !v.is_nan())
            .map(|(v, w)| v * w)
            .sum::<f64>()
    };
    let mut heat = 0.0;
    for n in 0..sizes.len() {
        heat += sizes[n] as f64 * denser[n] * starts[n];
    }
    let (mut low, mut high) = (-1e4_f64, 1e4_f64);
    assert!(held(&chain(low)) < heat && held(&chain(high)) > heat);
    for _ in 0..200 {
        let mid = 0.5 * (low + high);
        if held(&chain(mid)) < heat {
            low = mid;
        } else {
            high = mid;
        }
    }
    chain(0.5 * (low + high))
}

/// **Heat through a part that loses nothing is answered**, judged singular, and refused where it
/// has no balance above absolute zero: two to four bodies, sealed, a femtowatt to a microwatt
/// generated in the first cell and drawn from the last, every cell on the chain it crosses at the
/// level the heat the starts hold sets — of which `SHAPE` reads the level and the clearances, a
/// face's `p/WHOLE` and a half cell's `p/HALF` being below it. Where the chain at that level leaves
/// a cell below zero it is refused as that. Measured on the step as it was: two bodies from 40 K
/// and 0 K, and the other way round, answered as now; three from 40 K, 0 K and 0 K with a picowatt
/// not converging; three from 300 K with a nanowatt, and four with a microwatt, refused as below
/// zero, where they balance with their coldest cells at 99.99 K and 51.05 K; and of the six with
/// no balance above zero, four answered — from 1 K with a nanowatt through two bodies and through
/// three, a picowatt across emissivity 1e-11 from 40 K, and a picowatt between two bodies of one
/// cell, the first a million times as dense, from 0 K and 1 nK, with every cell within a nanokelvin
/// of zero — and two refused, as now. The last two are near zero because only there does the
/// step's ask decide the verdict: a step held to its cap moves by the cap, so one the stop can take
/// for a balance is one of a block within a few microkelvin of zero; and four bodies of one cell
/// from 0 K and 1e-12 K, a femtowatt through them, are refused only because a body taken out asks
/// for its step.
#[test]
fn heat_through_a_part_that_loses_nothing_is_answered() {
    const E: f64 = EMISSIVITY;
    #[rustfmt::skip]
    let balanced: [Drawn; 5] = [
        (&[2, 2], &[E, E], &[1.0; 2], 1e-9, &[40.0, 0.0]),
        (&[2, 2], &[E, E], &[1.0; 2], 1e-9, &[0.0, 40.0]),
        (&[2, 2, 2], &[E, E, E], &[1.0; 3], 1e-12, &[40.0, 0.0, 0.0]),
        (&[2, 2, 2], &[E, E, E], &[1.0; 3], 1e-9, &[300.0, 0.0, 0.0]),
        (&[2, 2, 2, 2], &[E, E, E, E], &[1.0; 4], 1e-6, &[300.0, 0.0, 0.0, 0.0]),
    ];
    for (sizes, emissivity, denser, p, starts) in balanced {
        let want = through(sizes, emissivity, denser, p, starts);
        let coldest = want
            .iter()
            .filter(|v| !v.is_nan())
            .fold(f64::INFINITY, |m, v| m.min(*v));
        assert!(coldest > 0.0, "a balance above absolute zero: {want:?}");
        let (b, void) = drawn(sizes, emissivity, denser, p, starts);
        let why = format!(
            "{p:e} W through bodies of {sizes:?} cells, emissivities {emissivity:?}, \
             {denser:?} as dense, sealed, from {starts:?} K"
        );
        let worst = answered(&why, &b, &void, &want);
        println!("  {why}: {worst:.2e} K from it, its coldest cell at {coldest} K");
    }
    #[rustfmt::skip]
    let unbalanced: [Drawn; 6] = [
        (&[2, 2], &[E, E], &[1.0; 2], 1e-9, &[0.0, 0.0]),
        (&[2, 2], &[E, E], &[1.0; 2], 1e-9, &[1.0, 0.0]),
        (&[2, 2], &[E, 1e-11], &[1.0; 2], 1e-12, &[40.0, 0.0]),
        (&[2, 2, 2], &[E, E, E], &[1.0; 3], 1e-9, &[1.0, 0.0, 0.0]),
        (&[1, 1], &[E, E], &[1e6, 1.0], 1e-12, &[0.0, 1e-9]),
        (&[1, 1, 1, 1], &[E, E, E, E], &[1.0; 4], 1e-15, &[0.0, 0.0, 0.0, 1e-12]),
    ];
    for (sizes, emissivity, denser, p, starts) in unbalanced {
        let want = through(sizes, emissivity, denser, p, starts);
        let coldest = want
            .iter()
            .filter(|v| !v.is_nan())
            .fold(f64::INFINITY, |m, v| m.min(*v));
        assert!(coldest < 0.0, "no balance above absolute zero: {want:?}");
        let why = format!(
            "{p:e} W through bodies of {sizes:?} cells, emissivities {emissivity:?}, \
             {denser:?} as dense, sealed, from {starts:?} K"
        );
        let (b, _) = drawn(sizes, emissivity, denser, p, starts);
        let e = b.steady_state().expect_err(&format!(
            "{why}: no balance above absolute zero, its coldest cell at {coldest} K"
        ));
        assert!(
            e.quantity.contains("no balance above absolute zero"),
            "{why}: refused as below absolute zero: {e}"
        );
        println!(
            "  {why}: refused, its coldest cell at {coldest} K: {}",
            e.quantity
        );
    }
}

/// Every solid cell of `b` started at `cells`, `void` naming the empty ones.
fn started(mut b: Solid3D, void: &[bool], cells: &[f64]) -> Solid3D {
    for (i, v) in cells.iter().enumerate() {
        if !void[i] {
            b.set_temperature(i, 0, 0, Temperature::from_si(*v));
        }
    }
    b
}

/// How far from 0 K the sink body of [`a_sink_at_a_balance_of_zero_is_not_refused`] may be read,
/// K: the band over which its own total — its microwatt in and out cancelling — does not change by
/// a rounding, ±0.037 K, computed at 80 digits by a review of this change. Within it the balance
/// cannot tell its level from zero, and a refusal as below zero is no more a reading of it.
const SINK_READ: f64 = 0.04;

/// **A sink at a balance of zero is not refused as below it**, judged singular: three bodies of
/// emissivities 0.001, 0.001 and 0.9, a microwatt generated in the second and drawn from the
/// third, whose film of 1e-21 m² with convection of 10 W/m²·K from a room at 0 K carries nothing
/// at the balance, its sink's cells at exactly 0 K — exact, the first two being held only through
/// pairs that carry nothing. The first two to `SHAPE`; the third at or above zero, within
/// `SINK_READ`. Measured before the change: from 300 K, 1 mK and 1000 K answered 5.9e-3 K out;
/// from the closed form with the third 0.1 K up, refused as below zero.
#[test]
fn a_sink_at_a_balance_of_zero_is_not_refused() {
    let p = 1e-6;
    let t4 = odd_root(fourth(p / HALF) + p / pair(0.001, EMISSIVITY)) + p / HALF;
    let want = [t4, t4, f64::NAN, t4, t4, f64::NAN, 0.0, 0.0];
    let mut froms = vec![[300.0, 300.0, 0.0, 0.001, 0.001, 0.0, 1000.0, 1000.0]];
    for up in [1e-3, 1e-2, 0.1] {
        froms.push([t4, t4, 0.0, t4, t4, 0.0, up, up]);
    }
    for from in froms {
        let (b, void) = bodies(
            &[2, 2, 2],
            &[0.001, 0.001, EMISSIVITY],
            Face::XMax,
            film(0.0, 10.0, 1e-21),
            &[0.0, 0.0, 0.0],
        );
        let b = started(b, &void, &from)
            .dissipating(p, |i, _, _| i == 4)
            .dissipating(-p, |i, _, _| i == 6);
        let why = format!("a microwatt into a sink at 0 K behind a film, from {from:?} K");
        let t = b
            .steady_state()
            .unwrap_or_else(|e| panic!("{why}: there is a balance to strike: {e}"));
        for c in [0, 1, 3, 4] {
            let miss = t[c] - want[c];
            assert!(
                miss.abs() <= SHAPE,
                "{why}: cell {c} {miss:+.3e} K from {} K",
                want[c]
            );
        }
        for c in [6, 7] {
            assert!(
                (-1e-4..=SINK_READ).contains(&t[c]),
                "{why}: the sink's cell {c} at {:e} K, beyond what its balance reads of zero",
                t[c]
            );
        }
        println!(
            "  {why}: the first two at {} K, the sink at {:e} K",
            t[0], t[6]
        );
    }
}

/// **A level its cells' conduction cancels in is read off its boundary**, judged singular: two
/// bodies of emissivity 1e-13, the first of two cells under a film of 1e-25 m² with convection of
/// 10 W/m²·K from a room at 0 K, the second of three, a milliwatt generated in its last cell and
/// drawn from its first. The second's only path out is a pair of 2.8e-27 W/K⁴ near absolute zero,
/// and its level is held by nothing else: what its cells' net fluxes sum to is the rounding of the
/// milliwatt they conduct, 5e-21 W, over a slope next to nothing. Measured before the change: from
/// 0 K and 1 K not converging, every round asking 1e10 K and more; from 1e-6 K refused as below
/// zero; and from the closed form with the second 1e-3 K to 0.1 K up, not converging or refused.
#[test]
fn a_level_its_cells_cancel_in_is_read_off_its_boundary() {
    let p = 1e-3;
    let want = [0.0, 0.0, f64::NAN, 0.0, p / WHOLE, 2.0 * p / WHOLE];
    let mut froms = Vec::new();
    for start in [0.0, 1e-6, 1.0] {
        froms.push([start, start, 0.0, start, start, start]);
    }
    for up in [1e-3, 1e-2, 5e-2, 0.1] {
        froms.push([0.0, 0.0, 0.0, up, want[4] + up, want[5] + up]);
    }
    for from in froms {
        let (b, void) = bodies(
            &[2, 3],
            &[1e-13, 1e-13],
            Face::XMin,
            film(0.0, 10.0, 1e-25),
            &[0.0, 0.0],
        );
        let b = started(b, &void, &from)
            .dissipating(-p, |i, _, _| i == 3)
            .dissipating(p, |i, _, _| i == 5);
        let why = format!("a milliwatt within a body near 0 K, from {from:?} K");
        let worst = answered(&why, &b, &void, &want);
        println!("  {why}: {worst:.2e} K from it");
    }
}

/// **A body taken out of a step judged singular is moved to its own balance**, the level at which
/// what crosses its boundary comes to nothing, not stepped on a slope next to nothing. Four bodies
/// of emissivities 0.5, 0.9, 0.001 and 1e-13, a film of 1e-25 m² with no convection from a room at
/// 0 K on the first, a picowatt generated in the third and drawn from the second, from 1 K, 0 K,
/// 0 K and 0 K: the third and fourth to `SHAPE` at 11.52 K; the first two, whose common level only
/// that film holds — εσA of 2.8e-33 W/K⁴ beside a picowatt, so that what f64 reads of it spans
/// kelvin — at or above zero. And three of emissivities 1e-11, 1e-15 and 1e-8, a film of 1e-21 m²
/// with convection of 10 W/m²·K on the last, a picowatt generated in the first and drawn from the
/// second, from where a review started them, the first at 11524 K: every cell to `SHAPE`. Measured
/// before the change: the four refused as below zero, the sink stepped down by a quarter of the
/// block's level a round and the source up by as little until the sink was drawn through zero; the
/// three answered 6.6e-3 K out. Started from the closed form with the second 1e-3 K to 0.1 K up,
/// the three are answered where they start, or 7.5e-3 K out, as before: its level a balance in f64
/// over that span, and no body taken out.
#[test]
fn a_body_taken_out_is_moved_to_its_own_balance() {
    let p = 1e-12;
    let t5 = odd_root(fourth(p / HALF) + p / pair(EMISSIVITY, 0.001)) + p / HALF;
    let t6 = t5 + p / WHOLE;
    let want = [0.0, f64::NAN, 0.0, 0.0, f64::NAN, t5, t6, f64::NAN, t6];
    let (b, _) = bodies(
        &[1, 2, 2, 1],
        &[0.5, EMISSIVITY, 0.001, 1e-13],
        Face::XMin,
        film(0.0, 0.0, 1e-25),
        &[1.0, 0.0, 0.0, 0.0],
    );
    let b = b
        .dissipating(-p, |i, _, _| i == 3)
        .dissipating(p, |i, _, _| i == 6);
    let why = "a picowatt across emissivities 0.9 and 0.001, from 1 K and 0 K";
    let t = b
        .steady_state()
        .unwrap_or_else(|e| panic!("{why}: there is a balance to strike: {e}"));
    for c in [5, 6, 8] {
        let miss = t[c] - want[c];
        assert!(
            miss.abs() <= SHAPE,
            "{why}: cell {c} {miss:+.3e} K from {} K",
            want[c]
        );
    }
    for c in [0, 2, 3] {
        assert!(t[c] >= -1e-4, "{why}: cell {c} at {:e} K", t[c]);
    }
    println!(
        "  {why}: the third at {} K, {t5} K; the first two at {:e} K",
        t[5], t[0]
    );
    let t2 = odd_root(fourth(p / HALF) + p / pair(1e-11, 1e-15)) + p / HALF;
    let t1 = t2 + p / WHOLE;
    let want = [t1, t1, t2, f64::NAN, 0.0, 0.0, f64::NAN, 0.0];
    let froms = [[
        11524.12400014308,
        11524.12401165568,
        11524.124000120011,
        0.0,
        0.0,
        0.0,
        0.0,
        3.8725919187219105e-121,
    ]];
    for from in froms {
        let (b, void) = bodies(
            &[3, 2, 1],
            &[1e-11, 1e-15, 1e-8],
            Face::XMax,
            film(0.0, 10.0, 1e-21),
            &[0.0, 0.0, 0.0],
        );
        let b = started(b, &void, &from)
            .dissipating(p, |i, _, _| i == 1)
            .dissipating(-p, |i, _, _| i == 4);
        let why = format!("a picowatt across emissivities 1e-11 and 1e-15, from {from:?} K");
        let worst = answered(&why, &b, &void, &want);
        println!("  {why}: {worst:.2e} K from it, the first at {t1} K");
    }
}

/// **A cell that generates nothing is not swept past what it exchanges with**: bodies of one or two
/// cells near absolute zero whose only losses are films of 1e-21 m² or 1e-25 m² with convection of
/// 10 W/m²·K from a room at 0 K, nothing generated, every cell at the room. A cell's sweep reads a
/// pair at its own end as `4C·|T|³`, nothing at 0 K, so a cell at 0 K fed across a clearance was
/// stepped by what it is fed over its film's slope alone: measured, one of emissivity 1e-11 fed
/// 9e-20 W from a neighbour at 20 K was stepped to 6.8e4 K in the first round, and the bodies beside
/// it thrown by thousands of kelvin after. The step as it was answered these, wandering there and
/// back; with a body taken out of a singular step moved to its own balance, a first form of this
/// change took the thrown neighbour for where to go, and none of the eleven converged. The middle
/// one of the first, moved from 0 K to its own balance at 427 K, had the last one's sweep throw it
/// to 1.9e6 K, and the field swung by 1e16 K.
#[test]
fn a_cell_that_generates_nothing_is_not_swept_past_what_it_exchanges_with() {
    const E: f64 = EMISSIVITY;
    /// Each body's cells and emissivity, the film's face and area, and each body's start, K.
    type Swept = (&'static [usize], &'static [f64], Face, f64, &'static [f64]);
    #[rustfmt::skip]
    let cases: [Swept; 6] = [
        (&[1, 1, 1], &[1e-11, 1e-11, 0.001], Face::XMax, 1e-21, &[1000.0, 0.0, 0.0]),
        (&[1, 1, 1], &[E, 1e-11, 0.001], Face::XMax, 1e-21, &[0.0, 20.0, 1000.0]),
        (&[2, 1, 1], &[E, 0.001, E], Face::XMax, 1e-21, &[0.0, 0.0, 1000.0]),
        (&[1, 1, 2], &[1e-11, E, 1e-11], Face::XMin, 1e-25, &[0.0, 20.0, 0.0]),
        (&[1, 1, 2], &[1e-11, 0.001, 1e-11], Face::XMin, 1e-25, &[1000.0, 0.0, 0.0]),
        (&[1, 1, 2], &[0.001, 0.001, 0.001], Face::XMin, 1e-25, &[1000.0, 0.0, 0.0]),
    ];
    for (sizes, emissivity, face, area, starts) in cases {
        let (b, void) = bodies(sizes, emissivity, face, film(0.0, 10.0, area), starts);
        let why = format!(
            "bodies of {sizes:?} cells, emissivities {emissivity:?}, a film of {area:e} m2 on \
             {face:?}, from {starts:?} K"
        );
        let worst = answered(&why, &b, &void, &at_room(&void, 0.0));
        println!("  {why}: {worst:.2e} K from the room");
    }
}

/// **A sink's part drawn down toward zero together is refused as having no balance above it**:
/// two bodies of one cell, a radiating film from a room at 0 K on the first, the second drawing a
/// microwatt or a picowatt across a clearance one of whose faces is of emissivity 1e-13 — no film
/// from a room at zero can bring it, so it has no balance above zero, its one balance far below.
/// With the first moved to its own balance every round, both fell toward zero together, the sink's
/// step held to a quarter of a level that went with them, until it asked 1.3e308 K and the field
/// was not a number: refused as having diverged, where the step as it was refused it as below zero.
/// Both are drawn down until that step overflows, and are refused as that, naming the sink's cell,
/// at `−∞`.
#[test]
fn a_sinks_part_drawn_down_together_is_refused_as_below_zero() {
    for (emissivity, p, starts) in [
        ([EMISSIVITY, 1e-13], 1e-6, [0.0, 1000.0]),
        ([1e-13, EMISSIVITY], 1e-12, [20.0, 1.0]),
    ] {
        let (b, _) = bodies(
            &[1, 1],
            &emissivity,
            Face::XMin,
            film(0.0, 0.0, AREA),
            &starts,
        );
        let b = b.dissipating(-p, |i, _, _| i == 2);
        let why = format!("{p:e} W drawn across emissivities {emissivity:?}, from {starts:?} K");
        let e = b
            .steady_state()
            .expect_err(&format!("{why}: no balance above absolute zero"));
        assert!(
            e.quantity
                == "no balance above absolute zero was found: cell (2, 0, 0) was drawn down by \
                    the sinks in its part until the step on its level overflowed"
                && e.before == f64::NEG_INFINITY,
            "{why}: refused as having no balance above absolute zero: {e}"
        );
        println!("  {why}: {}", e.quantity);
    }
}
