//! `Solid3D::steady_state` of a part that loses nothing, and of the bodies a clearance joins: the
//! level a part's own heat sets, a refusal where heat is generated with nowhere to go, and the
//! levels of bodies that see each other taken together.
//!
//! **A part that loses nothing balances at any level.** Every face term cancels under a uniform
//! shift, so the equations fix the shape of the field and not where it sits, and the solve
//! answered wherever its relaxation drifted. Measured before this, on blocks of 1 mm cells:
//!
//! - four cells insulated on every face, one at 100 °C and three at 20 °C, at **14.48 °C** —
//!   colder than the coldest of them — where they march to 40 °C, and at 14.48 °C again when the
//!   hot cell held four times the heat of the others, where they march to 65.71 °C;
//! - a body across a vacuum from one with a film, at **−8.05 °C** in a room at 0 °C, where it
//!   marches to the 60 °C it holds;
//! - two bodies that radiate only to each other, at **2653.69 °C**, where they march to 110 °C.
//!
//! `settle` put a block there, and `verify` measured what was still to come against it: a block
//! made like `19-a-coating-stops-the-heat` reported 18.69 K still to come where 10.58 K was.
//!
//! **Heat with nowhere to go was asked of the block, not of where the heat is.** A block with a
//! film, air or a clearance anywhere passed, so a block of `k = 0` given a film and a heated cell
//! returned the field it started from as its balance, and a film of zero, an area of zero, a
//! mounting or air of zero, a heated cell behind an insulator and a heated body across a vacuum
//! from the cooled one each swept two thousand times and reported that it had not converged.
//!
//! **Bodies a clearance joins each took their own step**, which closes an exchange from both ends
//! at once. Measured before this: under a film too weak to hold them, two radiating bodies came to
//! **2926.84 K** and called it a balance, in a room at 273.15 K, and under films of 1e-3 to 0.1
//! W/m²·K did not converge; a milliwatt in at one and out at the other answered 2923 K; three
//! bodies of one cell came to 2.7e17 K, and ten of two cells to 3055 K where they hold 307.15 K.
//!
//! # What is checked against what
//!
//! **The level is a conservation law.** The march conserves `Σ C·T` in a part that loses nothing,
//! so the steady state it goes to holds it, and the closed form for the level is that sum over the
//! part's capacity — with unequal capacities, so that a solve that averaged the cells would miss.
//! Held to the rounding of the shift that sets it. **The shape is each case's closed form**: one
//! level for a part with no source; Fourier's law across every face between a source and a sink,
//! `P/(k·dx)` a face; and across a clearance, Stefan–Boltzmann between the two surfaces, each half a
//! cell from its centre, carrying exactly what the source puts in.
//!
//! **Each refusal names a cell of the part and carries the part's watts**, and the part refused is
//! not the first the block lists, nor the first cell: a block whose cooled body is heated too names
//! the other, with that body's watts and not the block's; a block with two sealed parts names the
//! first. Along each axis, because a face that conducts nothing divides a body along all three.
//! **And the ways out still carry**: radiation with no convection, a mounting that is not zero, air
//! that is not, a clearance to the cooled body, a patch of film that reaches a body and one that
//! does not.

use pantometry_core::substance::ThermalProps;
use pantometry_core::units::{
    Area, Density, Length, SpecificHeat, Temperature, ThermalConductivity, ThermalExpansion,
    STEFAN_BOLTZMANN,
};
use pantometry_core::{Substance, Violation};
use pantometry_thermal::{Axis, Environment, Face, Solid3D};

/// How far a part's level may be from the heat it holds, K: sixteen times the relative rounding of
/// 400 K, `400ε`, which is 1.4e-12 K, about twenty-five ulps of 400 K. The last round adds one
/// shift to every cell, a rounding of half an ulp each, and the shift is a sum of `C·ΔT` over the
/// part divided by its capacity. Measured here, a third of an ulp at most.
const LEVEL: f64 = 16.0 * 400.0 * f64::EPSILON;

/// How far a cell may be from its closed form, K. The solve stops when a round moves the field by
/// less than a microkelvin, which leaves `1e-6·ρ/(1 − ρ)` for a round that contracts the error by
/// `ρ`: ten microkelvin is a contraction of 0.91. Measured here, 1.2e-6 K at most, in the block of
/// two materials in three dimensions; against levels that missed by 25 K and more.
const SHAPE: f64 = 1e-5;

/// How far what a clearance carries may be from what crosses it, as a share. A cell's `SHAPE`
/// through a pair's slope, `4C·T³`, is 1.1e-7 of a milliwatt at 390 K and 8e-8 of three at 504 K.
const CARRIED: f64 = 1e-6;

fn stuff(rho: f64, k: f64, emissivity: f64) -> Substance {
    Substance::bulk("stuff", Density::from_si(rho)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(k),
        specific_heat: SpecificHeat::j_per_kg_k(4000.0),
        expansion: ThermalExpansion::ppm_per_k(50.0),
        emissivity,
    })
}

/// A conductor that does not radiate.
fn grey() -> Substance {
    stuff(1000.0, 50.0, 0.0)
}

/// The same holding four times as much.
fn dense() -> Substance {
    stuff(4000.0, 50.0, 0.0)
}

/// A conductor that radiates, `k = 50` and `ε = 0.9` — what [`carries`] is written for.
fn black() -> Substance {
    stuff(1000.0, 50.0, 0.9)
}

/// A cell that conducts nothing.
fn insulator() -> Substance {
    stuff(1000.0, 0.0, 0.0)
}

fn film(h: f64, area: f64) -> Environment {
    Environment {
        ambient: Temperature::celsius(0.0),
        convection_w_per_m2_k: h,
        area: Area::from_si(area),
    }
}

const ROOM: f64 = 273.15;

/// A row of `n` cells of 1 mm at 20 °C.
fn row(n: usize, s: Substance) -> Solid3D {
    line(0, n, s)
}

/// A line of `n` cells of 1 mm at 20 °C along `axis`, 0, 1 or 2. Its `p`th cell is index `p`.
fn line(axis: usize, n: usize, s: Substance) -> Solid3D {
    let mut counts = [1, 1, 1];
    counts[axis] = n;
    Solid3D::new(
        "block",
        s,
        (counts[0], counts[1], counts[2]),
        Length::mm(1.0),
        Temperature::celsius(20.0),
    )
}

/// The `p`th cell along `axis`.
fn at(axis: usize, p: usize) -> (usize, usize, usize) {
    let mut c = [0, 0, 0];
    c[axis] = p;
    (c[0], c[1], c[2])
}

const LOW: [Face; 3] = [Face::XMin, Face::YMin, Face::ZMin];

/// A block shaped so that one cell is reached only going back along an axis: its name, its
/// counts, its void and that cell.
type Shape = (
    &'static str,
    (usize, usize, usize),
    fn(usize, usize, usize) -> bool,
    (usize, usize, usize),
);

/// A block with a cooled body and a sealed one: its name, the block, the cooled body's cells, the
/// sealed body's, and a cell that is a part of its own.
type Beside = (String, Solid3D, Vec<usize>, Vec<usize>, Option<usize>);
const AXES: [Axis; 3] = [Axis::X, Axis::Y, Axis::Z];
const NAMES: [&str; 3] = ["x", "y", "z"];

fn set(b: &mut Solid3D, i: usize, celsius: f64) {
    b.set_temperature(i, 0, 0, Temperature::celsius(celsius));
}

fn set_at(b: &mut Solid3D, cell: (usize, usize, usize), celsius: f64) {
    b.set_temperature(cell.0, cell.1, cell.2, Temperature::celsius(celsius));
}

/// Every cell's temperature, K, `NaN` for a void.
fn field(b: &Solid3D) -> Vec<f64> {
    let (nx, ny, nz) = b.counts();
    let mut out = Vec::new();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                out.push(b.temperature_at(i, j, k).to_si());
            }
        }
    }
    out
}

/// What the cells `which` hold over what they can, K — the level a part of them comes to.
fn level_of(b: &Solid3D, t: &[f64], which: &[usize]) -> f64 {
    let caps = b.cell_capacities();
    let held: f64 = which.iter().map(|c| caps[*c] * t[*c]).sum();
    held / which.iter().map(|c| caps[*c]).sum::<f64>()
}

/// The change in what the cells `which` hold, over their capacity, K: zero for a part whose heat
/// was conserved. Formed from the differences, as the solve forms its shift, so that it is not the
/// difference of two sums of four hundred kelvin.
fn drift(b: &Solid3D, t: &[f64], which: &[usize]) -> f64 {
    let caps = b.cell_capacities();
    let start = field(b);
    let moved: f64 = which.iter().map(|c| caps[*c] * (t[*c] - start[*c])).sum();
    moved / which.iter().map(|c| caps[*c]).sum::<f64>()
}

/// The largest of `misses` by size, and `NaN` if any is one.
fn worst(misses: impl Iterator<Item = f64>) -> f64 {
    misses.fold(0.0_f64, |w, d| {
        if d.is_nan() || w.is_nan() {
            f64::NAN
        } else if d.abs() > w.abs() {
            d
        } else {
            w
        }
    })
}

/// `t` holds `which` at `expect` K to `SHAPE`, and their heat to `LEVEL`.
fn holds(why: &str, b: &Solid3D, t: &[f64], which: &[usize], expect: f64) {
    let miss = worst(which.iter().map(|c| t[*c] - expect));
    let drifted = drift(b, t, which);
    println!(
        "  {why:56} {:.6} C, worst cell {miss:+.2e} K, heat {drifted:+.2e} K",
        expect - 273.15
    );
    assert!(
        miss.abs() <= SHAPE,
        "{why}: a cell is {miss:+.3e} K from {:.6} C",
        expect - 273.15
    );
    assert!(
        drifted.abs() <= LEVEL,
        "{why}: the part's heat moved by {drifted:+.3e} K of its capacity"
    );
}

/// `t` has every cell of `which` at `expect` K to `SHAPE`, with no claim about heat.
fn sits(why: &str, t: &[f64], which: &[usize], expect: f64) {
    let miss = worst(which.iter().map(|c| t[*c] - expect));
    println!("  {why:56} worst cell {miss:+.2e} K from {expect:.4} K");
    assert!(
        miss.abs() <= SHAPE,
        "{why}: a cell is {miss:+.3e} K from {expect:.4} K"
    );
}

/// What a clearance carries between a cell of conductivity `k_hot` and emissivity `e_hot` at
/// `hot` and one of `k_cold`, `e_cold` at `cold`, given that it carries `q`: Stefan–Boltzmann
/// between the two surfaces, `σA/(1/ε_a + 1/ε_b − 1)·(u⁴ − v⁴)`, each surface behind its half
/// cell, `2kA/dx`, by `q` over it.
fn crosses_between(hot: (f64, f64, f64), cold: (f64, f64, f64), q: f64) -> f64 {
    let ((t_hot, k_hot, e_hot), (t_cold, k_cold, e_cold)) = (hot, cold);
    let area = 1e-6;
    let coefficient = STEFAN_BOLTZMANN.to_si() * area / (1.0 / e_hot + 1.0 / e_cold - 1.0);
    let half = |k: f64| 2.0 * k * area / 1e-3;
    let (u, v) = (t_hot - q / half(k_hot), t_cold + q / half(k_cold));
    coefficient * ((u * u) * (u * u) - (v * v) * (v * v))
}

/// The clearance between the centres `hot` and `cold` of `black` carries `q`, to `CARRIED`.
fn carries(why: &str, hot: f64, cold: f64, q: f64) {
    carries_between(why, (hot, 50.0, 0.9), (cold, 50.0, 0.9), q);
}

/// The clearance between `hot` and `cold`, each `(T, k, ε)`, carries `q`, to `CARRIED`.
fn carries_between(why: &str, hot: (f64, f64, f64), cold: (f64, f64, f64), q: f64) {
    let crossed = crosses_between(hot, cold, q);
    let share = crossed / q - 1.0;
    println!("  {why:56} the clearance carries {share:+.2e} of {q:e} W too much");
    assert!(
        share.abs() <= CARRIED,
        "{why}: Stefan-Boltzmann across the surfaces is {crossed:e} W, not the {q:e} W put in"
    );
}

/// `refusal` says the heat warms without limit from the cells joined to `cell`, and carries `watts`.
fn warms(why: &str, refusal: &Violation, cell: (usize, usize, usize), watts: f64) {
    without_limit(why, refusal, "warms without limit", cell, watts);
}

fn without_limit(
    why: &str,
    refusal: &Violation,
    how: &str,
    cell: (usize, usize, usize),
    watts: f64,
) {
    println!("  {why:56} -> {refusal}");
    let (i, j, k) = cell;
    assert!(
        refusal.quantity.contains(how),
        "{why}: refused for the wrong reason: {refusal}"
    );
    assert!(
        refusal
            .quantity
            .contains(&format!("joined to ({i}, {j}, {k})")),
        "{why}: the refusal must name ({i}, {j}, {k}): {refusal}"
    );
    assert_eq!(
        refusal.before.to_bits(),
        watts.to_bits(),
        "{why}: the refusal must carry the part's {watts} W, carries {}",
        refusal.before
    );
    assert_eq!(refusal.site, "block", "{why}: refused at the wrong place");
}

fn solved(why: &str, b: &Solid3D) -> Vec<f64> {
    let t = b
        .steady_state()
        .unwrap_or_else(|e| panic!("{why}: there is a balance to strike: {e}"));
    assert!(
        t.iter().all(|v| v.is_finite()),
        "{why}: answered with a cell that is not a number: {t:?}"
    );
    t
}

fn refused(why: &str, b: &Solid3D) -> Violation {
    match b.steady_state() {
        Ok(t) => panic!("{why}: warms without limit and was answered {t:?}"),
        Err(e) => e,
    }
}

/// **A block that loses nothing comes to the heat it holds**, whatever its cells hold each.
#[test]
fn an_insulated_block_comes_to_the_heat_it_holds() {
    let mut b = row(4, grey());
    set(&mut b, 0, 100.0);
    let t = solved("four cells", &b);
    // Three at 20 °C and one at 100 °C, of equal capacity.
    holds(
        "four equal cells, one at 100 C",
        &b,
        &t,
        &[0, 1, 2, 3],
        313.15,
    );

    // The hot cell holds four times what each of the others does, so the level is not the mean of
    // the four temperatures — it is (4·373.15 + 4·293.15 + 2·293.15)/10. A second cell is the dense
    // one too, so the capacity is read for more than the first.
    let mut b = row(4, grey());
    b.fill(dense(), |i, _, _| i == 0 || i == 2);
    set(&mut b, 0, 100.0);
    let t = solved("two dense cells", &b);
    let expect = (4.0 * 373.15 + 4.0 * 293.15 + 2.0 * 293.15) / 10.0;
    holds(
        "two of four cells four times as dense",
        &b,
        &t,
        &[0, 1, 2, 3],
        expect,
    );
    assert!(
        (expect - level_of(&b, &field(&b), &[0, 1, 2, 3])).abs() < 1e-12,
        "the closed form is not the heat this block holds"
    );

    // In three dimensions, every cell somewhere else.
    let mut b = Solid3D::new(
        "block",
        grey(),
        (4, 3, 2),
        Length::mm(1.0),
        Temperature::celsius(20.0),
    );
    b.fill(dense(), |i, j, _| i + j == 3);
    for k in 0..2 {
        for j in 0..3 {
            for i in 0..4 {
                let c = (i + 4 * (j + 3 * k)) as f64;
                b.set_temperature(i, j, k, Temperature::celsius(10.0 + 7.0 * c % 31.0));
            }
        }
    }
    let all: Vec<usize> = (0..24).collect();
    let expect = level_of(&b, &field(&b), &all);
    let t = solved("in three dimensions", &b);
    holds(
        "4x3x2, two materials, every cell apart",
        &b,
        &t,
        &all,
        expect,
    );

    // A film that carries nothing is not a way out.
    let mut b = row(4, grey()).losing_from(Face::XMin, film(0.0, 1e-6));
    set(&mut b, 0, 100.0);
    let t = solved("a film of zero", &b);
    holds("a film of zero stated", &b, &t, &[0, 1, 2, 3], 313.15);
}

/// **A body is found the whole way round**: one per axis that can be walked whole only by going
/// back along that axis from where the search reached it, insulated and uneven, comes to the heat
/// all of it holds. A search that only went forward found two bodies in each, and each kept its
/// own heat.
#[test]
fn a_body_is_found_the_whole_way_round() {
    // The search starts at the lowest index and each far cell is reached last, back along its axis.
    let shapes: [Shape; 4] = [
        // (2,0) up to (2,1), then back along x to (1,1) and (0,1).
        (
            "back along x",
            (3, 2, 1),
            |i, j, _| j == 0 && i < 2,
            (0, 1, 0),
        ),
        // (0,0) up, along to (2,1), then back along y to (2,0).
        (
            "back along y",
            (3, 2, 1),
            |i, j, _| j == 0 && i == 1,
            (2, 0, 0),
        ),
        // (0,0,0) up, along to (2,0,1), then back along z to (2,0,0).
        (
            "back along z",
            (3, 1, 2),
            |i, _, k| k == 0 && i == 1,
            (2, 0, 0),
        ),
        // A column along y in the upper of two layers, joined by y faces alone: the faces of the
        // first layer are not the faces of the second, and a box joined every way round does not
        // say which a search read.
        (
            "along y, above the first layer",
            (1, 3, 2),
            |_, _, k| k == 0,
            (0, 2, 1),
        ),
    ];
    for (why, counts, void, far) in shapes {
        let mut b = Solid3D::new(
            "block",
            grey(),
            counts,
            Length::mm(1.0),
            Temperature::celsius(20.0),
        )
        .empty(void);
        b.set_temperature(far.0, far.1, far.2, Temperature::celsius(100.0));
        let (nx, ny, nz) = counts;
        let solid: Vec<usize> = (0..nx * ny * nz)
            .filter(|c| !b.is_void(c % nx, (c / nx) % ny, c / (nx * ny)))
            .collect();
        let expect = level_of(&b, &field(&b), &solid);
        let t = solved(why, &b);
        holds(why, &b, &t, &solid, expect);
    }
}

/// **A body that loses nothing keeps its heat beside one that does**, whichever comes first, and
/// whatever cuts it off: a vacuum, a cell that conducts nothing or a joint of no contact along each
/// axis, a patch of film that stops short of it.
#[test]
fn a_body_that_loses_nothing_keeps_its_heat_beside_one_that_does() {
    let mut cases: Vec<Beside> = vec![
        (
            "across a vacuum, the cooled body first".to_string(),
            row(5, grey())
                .empty(|i, _, _| i == 2)
                .losing_from(Face::XMin, film(10.0, 1e-6)),
            vec![0, 1],
            vec![3, 4],
            None,
        ),
        (
            "across a vacuum, the cooled body second".to_string(),
            row(5, grey())
                .empty(|i, _, _| i == 2)
                .losing_from(Face::XMax, film(10.0, 1e-6)),
            vec![3, 4],
            vec![0, 1],
            None,
        ),
    ];
    for a in 0..3 {
        let mut b = line(a, 5, grey()).losing_from(LOW[a], film(10.0, 1e-6));
        b.fill(insulator(), move |i, j, k| [i, j, k][a] == 2);
        cases.push((
            format!("behind a cell that conducts nothing, along {}", NAMES[a]),
            b,
            vec![0, 1],
            vec![3, 4],
            Some(2),
        ));
        let axis = AXES[a];
        cases.push((
            format!("behind a joint of no contact, along {}", NAMES[a]),
            line(a, 4, grey())
                .joined(move |face, i, j, k| (face == axis && [i, j, k][a] == 2).then_some(0.0))
                .losing_from(LOW[a], film(10.0, 1e-6)),
            vec![0, 1],
            vec![2, 3],
            None,
        ));
    }
    // A film on part of a face: the patch reaches the first body and the face reaches both.
    // Five cells by two, a void across the middle; the film on y = 0 over x in 0..2.
    cases.push((
        "beside a patch of film that stops short of it".to_string(),
        Solid3D::new(
            "block",
            grey(),
            (5, 2, 1),
            Length::mm(1.0),
            Temperature::celsius(20.0),
        )
        .empty(|i, _, _| i == 2)
        .losing_from_within(Face::YMin, [0, 0, 2, 1], film(10.0, 2e-6)),
        vec![0, 1, 5, 6],
        vec![3, 4, 8, 9],
        None,
    ));
    for (why, mut b, cooled, sealed, apart) in cases {
        // The sealed body uneven and denser in one cell, the cooled one warm.
        let (nx, ny, _) = b.counts();
        let cell = |c: usize| (c % nx, (c / nx) % ny, c / (nx * ny));
        set_at(&mut b, cell(sealed[0]), 100.0);
        let heavy = cell(sealed[1]);
        b.fill(dense(), move |i, j, k| (i, j, k) == heavy);
        set_at(&mut b, heavy, 20.0);
        set_at(&mut b, cell(cooled[0]), 60.0);
        let keeps = level_of(&b, &field(&b), &sealed);
        let apart_at = apart.map(|c| field(&b)[c]);
        let t = solved(&why, &b);
        holds(&format!("{why}: the sealed body"), &b, &t, &sealed, keeps);
        sits(&format!("{why}: the cooled body"), &t, &cooled, ROOM);
        // A cell that conducts nothing is a part of its own, and keeps its temperature exactly.
        if let (Some(c), Some(at)) = (apart, apart_at) {
            assert_eq!(
                t[c].to_bits(),
                at.to_bits(),
                "{why}: a cell that conducts nothing moved"
            );
        }
    }

    // A block of nothing but cells that conduct nothing is a part a cell, and every one stays.
    let mut b = row(4, insulator()).losing_from(Face::XMin, film(10.0, 1e-6));
    set(&mut b, 0, 100.0);
    let t = solved("nothing conducts", &b);
    let start = field(&b);
    assert!(
        t.iter()
            .zip(&start)
            .all(|(a, s)| a.to_bits() == s.to_bits()),
        "a block that conducts nothing moved: {t:?} from {start:?}"
    );
}

/// **Bodies that see only each other come to one level**, the heat the whole system holds: two,
/// one of them heavier; three in a row; three of one cell each, which their pairs alone hold; and
/// ten, in a block longer than any of them.
#[test]
fn bodies_that_see_only_each_other_come_to_one_level() {
    let heavy = || stuff(4000.0, 50.0, 0.9);

    let mut b = row(5, black()).empty(|i, _, _| i == 2);
    set(&mut b, 0, 100.0);
    set(&mut b, 3, 300.0);
    let all = [0, 1, 3, 4];
    let t = solved("two bodies", &b);
    holds("two bodies across a radiating vacuum", &b, &t, &all, 383.15);

    let mut b = row(5, black()).empty(|i, _, _| i == 2);
    b.fill(heavy(), |i, _, _| i == 4);
    set(&mut b, 0, 100.0);
    set(&mut b, 3, 300.0);
    let expect = level_of(&b, &field(&b), &all);
    let t = solved("two bodies, one heavier", &b);
    holds(
        "the same, the second body's far cell four times as dense",
        &b,
        &t,
        &all,
        expect,
    );

    let mut b = row(8, black()).empty(|i, _, _| i == 2 || i == 5);
    set(&mut b, 0, 100.0);
    set(&mut b, 4, 300.0);
    set(&mut b, 7, 200.0);
    let all = [0, 1, 3, 4, 6, 7];
    let expect = level_of(&b, &field(&b), &all);
    let t = solved("three bodies", &b);
    holds(
        "three bodies in a row, each seeing its neighbours",
        &b,
        &t,
        &all,
        expect,
    );

    // A cell with no face to a neighbour has only its pairs; its sweep took it 1250 K at once.
    let mut b = row(5, black()).empty(|i, _, _| i == 1 || i == 3);
    set(&mut b, 2, 300.0);
    b.fill(heavy(), |i, _, _| i == 4);
    let all = [0, 2, 4];
    let expect = level_of(&b, &field(&b), &all);
    let t = solved("three cells", &b);
    holds(
        "three bodies of one cell, one heavier, one at 300 C",
        &b,
        &t,
        &all,
        expect,
    );

    // A body at absolute zero beside one at 400 K: their heat sets 200 K. In a first form of the step
    // the bodies take together, a quarter of the cold one's own level held it still while the warm
    // one was swept down to it, to 0 K for both.
    let mut b = Solid3D::new(
        "block",
        black(),
        (3, 1, 1),
        Length::mm(1.0),
        Temperature::from_si(400.0),
    )
    .empty(|i, _, _| i == 1);
    b.set_temperature(2, 0, 0, Temperature::from_si(0.0));
    let t = solved("one at absolute zero", &b);
    holds(
        "a body at absolute zero and one at 400 K",
        &b,
        &t,
        &[0, 2],
        200.0,
    );

    // Two bodies holding 8e3 J/K each at 10 K and 20 K, whose pairs' slopes are 2e-11 W/K: the step
    // on their levels was judged singular and took nothing, and they answered where they began.
    let mut b = Solid3D::new(
        "block",
        stuff(1e9, 50.0, 0.05),
        (5, 1, 1),
        Length::mm(1.0),
        Temperature::from_si(10.0),
    )
    .empty(|i, _, _| i == 2);
    b.set_temperature(3, 0, 0, Temperature::from_si(20.0));
    b.set_temperature(4, 0, 0, Temperature::from_si(20.0));
    let t = solved("of great capacity", &b);
    holds(
        "two bodies of 8e3 J/K at 10 K and 20 K",
        &b,
        &t,
        &[0, 1, 3, 4],
        15.0,
    );

    // Ten bodies of two cells in a row of twenty-nine: each relaxed as two long, not twenty-nine.
    let mut b = row(29, black()).empty(|i, _, _| i % 3 == 2);
    set(&mut b, 0, 300.0);
    let all: Vec<usize> = (0..29).filter(|i| i % 3 != 2).collect();
    let expect = level_of(&b, &field(&b), &all);
    let t = solved("ten bodies", &b);
    holds(
        "ten bodies of two cells, one cell at 300 C",
        &b,
        &t,
        &all,
        expect,
    );
}

/// **What a clearance carries between bodies is Stefan–Boltzmann's**, whether the bodies lose
/// nothing — a source in one and as much a sink in the other — or the heat leaves through a film
/// on the far one. And under a film too weak to matter, bodies with no source sit at the room.
#[test]
fn a_clearance_carries_what_crosses_it_between_bodies() {
    // A source in the second body and as much a sink in the first, at two sizes: each body's own
    // face carries the watt, `P/(k·dx)` across it, and the clearance carries it between the two
    // surfaces, at levels the heat they hold sets.
    for p in [1e-3, 3e-3] {
        let b = row(5, black())
            .empty(|i, _, _| i == 2)
            .dissipating(-p, |i, _, _| i == 0)
            .dissipating(p, |i, _, _| i == 4);
        let why = format!("{p} W in at one body and out at the other");
        let t = solved(&why, &b);
        println!(
            "  {why:56} {:.4} K and {:.4} K",
            0.5 * (t[0] + t[1]),
            0.5 * (t[3] + t[4])
        );
        carries(&why, t[3], t[1], p);
        let drop = p / (50.0 * 1e-3);
        sits(
            &format!("{why}: the sink's face"),
            &[t[1] - t[0]],
            &[0],
            drop,
        );
        sits(
            &format!("{why}: the source's face"),
            &[t[4] - t[3]],
            &[0],
            drop,
        );
        let drifted = drift(&b, &t, &[0, 1, 3, 4]);
        assert!(
            drifted.abs() <= LEVEL,
            "{why}: the part's heat moved by {drifted:+.3e} K"
        );
    }

    // Heat in the second body, out through a film on the first: all of it crosses.
    let b = row(5, black())
        .empty(|i, _, _| i == 2)
        .losing_from(Face::XMin, film(10.0, 1e-6))
        .dissipating(0.1, |i, _, _| i == 4);
    let t = solved("out through the other body's film", &b);
    carries(
        "a tenth of a watt out through the other body",
        t[3],
        t[1],
        0.1,
    );

    // The first body's face bare under a film too weak to hold the pair: both at the room.
    for h in [1e-6, 1e-2, 0.1] {
        let mut b = row(5, black())
            .empty(|i, _, _| i == 2)
            .losing_from(Face::XMin, film(h, 1e-6));
        b.fill(grey(), |i, _, _| i == 0);
        set(&mut b, 0, 100.0);
        set(&mut b, 3, 300.0);
        let why = format!("two bodies under a bare film of {h} W/m2/K");
        let t = solved(&why, &b);
        sits(&why, &t, &[0, 1, 3, 4], ROOM);
    }

    // A sink fed only across a clearance, from a body a film warms by radiating to a room at
    // 1000 K, from a start at 300 K: its own cell runs down while the heat that would feed it is
    // still arriving. While the pairs read a negative temperature's fourth power as even, a step
    // on the bodies' levels held to a quarter of them was what kept it above zero long enough —
    // measured without the quarter, not converging — and a milliwatt here, from the same start,
    // was refused at the mirror of its balance. Read odd, it needs neither, and every start finds
    // the balance: `a_sink_is_answered_by_what_it_takes_not_where_it_starts`.
    let mut b = Solid3D::new(
        "block",
        stuff(1000.0, 5.0, 0.9),
        (6, 1, 1),
        Length::mm(1.0),
        Temperature::from_si(300.0),
    )
    .empty(|i, _, _| i == 3)
    .losing_from(
        Face::XMin,
        Environment {
            ambient: Temperature::from_si(1000.0),
            convection_w_per_m2_k: 0.0,
            area: Area::from_si(1e-6),
        },
    );
    b.fill(stuff(1000.0, 0.05, 0.05), |i, _, _| i == 4);
    b.fill(stuff(1000.0, 50.0, 0.9), |i, _, _| i == 5);
    let b = b.dissipating(-2e-3, |i, _, _| i == 4);
    let t = solved("a sink fed across a clearance", &b);
    carries_between(
        "two milliwatts to a sink across a clearance",
        (t[2], 5.0, 0.9),
        (t[4], 0.05, 0.05),
        2e-3,
    );

    // A body joined to the cooled one through another body's clearance, with a cell that conducts
    // nothing between them: three by three, the bodies in the first row, a void, the cooled one in
    // the third. Its heat crosses its own clearance; the other's carries nothing.
    let mut b = Solid3D::new(
        "block",
        black(),
        (3, 3, 1),
        Length::mm(1.0),
        Temperature::celsius(20.0),
    )
    .empty(|_, j, _| j == 1)
    .losing_from(Face::YMax, film(10.0, 3e-6));
    b.fill(stuff(1000.0, 0.0, 0.9), |i, j, _| i == 1 && j == 0);
    let b = b.dissipating(0.01, |i, j, _| i == 0 && j == 0);
    let t = solved("joined through another's clearance", &b);
    carries("joined through another body's clearance", t[0], t[6], 0.01);
    sits(
        "the other body, whose clearance carries nothing",
        &[t[2] - t[8]],
        &[0],
        0.0,
    );
    assert_eq!(
        t[1].to_bits(),
        293.15_f64.to_bits(),
        "the cell between them conducts nothing and moved"
    );
}

/// **A source and a sink in a part that loses nothing** have a balance: Fourier's law across every
/// face between them, at the level the part's heat sets. As much in as out to the rounding of
/// spreading it is a balance; a few times that rounding more is not.
#[test]
fn a_source_and_a_sink_in_a_part_that_loses_nothing() {
    // A bar of eight cells, a watt in at one end and out at the other. Every face carries the watt,
    // so each drops `P·dx/(k·dx²)` = 1/(50 × 1e-3) = 20 K, and the line is centred on the level.
    let p = 1e-3;
    let mut b = row(8, grey())
        .dissipating(p, |i, _, _| i == 0)
        .dissipating(-p, |i, _, _| i == 7);
    set(&mut b, 3, 100.0);
    let all: Vec<usize> = (0..8).collect();
    let level = level_of(&b, &field(&b), &all);
    let drop = p / (50.0 * 1e-3);
    let t = solved("a source and a sink", &b);
    let line: Vec<f64> = all
        .iter()
        .map(|c| t[*c] - drop * (3.5 - *c as f64))
        .collect();
    sits("a source and a sink, along the line", &line, &all, level);
    let drifted = drift(&b, &t, &all);
    assert!(
        drifted.abs() <= LEVEL,
        "the bar's heat moved by {drifted:+.3e} K"
    );

    // As much in as out, spread so that it does not add back to zero: still a balance, and each
    // face carries what the cells before it put in. The net is added as the block adds it, cell by
    // cell in index order; it was refused as warming without limit while the test was whether the
    // sum was zero.
    let b = row(10, grey())
        .dissipating(0.1, |i, _, _| i < 3)
        .dissipating(-0.1, |i, _, _| i >= 3);
    let each = |i: usize| if i < 3 { 0.1 / 3.0 } else { -0.1 / 7.0 };
    let net = (0..10).fold(0.0, |s, i| s + each(i));
    assert!(
        net != 0.0,
        "the case needs a net that rounding leaves, and this adds to {net}"
    );
    println!("  as much in as out, spread unevenly: {net:+.2e} W net");
    let t = solved("as much in as out, spread unevenly", &b);
    let mut through = 0.0;
    let mut faces = Vec::new();
    for i in 0..9 {
        through += each(i);
        faces.push((t[i] - t[i + 1]) - through / (50.0 * 1e-3));
    }
    sits(
        "as much in as out, each face against what came before",
        &faces,
        &(0..9).collect::<Vec<_>>(),
        0.0,
    );

    // A net of 3.6 times the rounding of one sum of its watts, under the 71 times a spread over 71
    // cells may leave: still a balance.
    let b = row(71, grey())
        .dissipating(0.7, |i, _, _| i < 38)
        .dissipating(-0.7, |i, _, _| i >= 38);
    let each = |i: usize| -> f64 {
        if i < 38 {
            0.7 / 38.0
        } else {
            -0.7 / 33.0
        }
    };
    let (net, gross) = (0..71).fold((0.0_f64, 0.0_f64), |(n, g), i| {
        (n + each(i), g + each(i).abs())
    });
    assert!(
        net.abs() > 3.0 * f64::EPSILON * gross && net.abs() <= 71.0 * f64::EPSILON * gross,
        "the case needs a net between a few and n roundings: {net:e} of {gross}"
    );
    solved("a balance that rounding leaves 3.6 roundings of", &b);

    // Three times what a spread over its cells may leave, more in than out: no balance, and the
    // part's net is the number.
    let p = 0.1;
    let less = p - 3.0 * 8.0 * f64::EPSILON * 2.0 * p;
    let b = row(8, grey())
        .dissipating(p, |i, _, _| i == 0)
        .dissipating(-less, |i, _, _| i == 7);
    let net = 0.0 + p + -less;
    assert!(
        net > 8.0 * f64::EPSILON * (p + less),
        "the case needs a net past what rounding leaves: {net:e}"
    );
    warms(
        "three roundings' worth more in than out",
        &refused("past the rounding", &b),
        (0, 0, 0),
        net,
    );

    // More in than out, by a part in a thousand.
    let b = row(4, grey())
        .dissipating(0.1, |i, _, _| i == 0)
        .dissipating(-0.0999, |i, _, _| i == 3);
    warms(
        "a tenth of a watt in and 0.0999 out",
        &refused("more in than out", &b),
        (0, 0, 0),
        0.1 + -0.0999,
    );

    // More out than in: it cools without limit, and says so.
    let b = row(4, grey()).dissipating(-0.1, |i, _, _| i == 2);
    without_limit(
        "a sink and nothing to feed it",
        &refused("a sink alone", &b),
        "cools without limit",
        (0, 0, 0),
        -0.1,
    );
}

/// **Heat generated where nothing carries it away is refused where it is**, before any sweep, and
/// what does carry it still does.
#[test]
fn heat_with_nowhere_to_go_is_refused_where_it_is() {
    let w = 0.1;
    let heated = |b: Solid3D, i: usize| b.dissipating(w, move |c, _, _| c == i);
    let cooled = || row(4, grey()).losing_from(Face::XMin, film(10.0, 1e-6));

    let mut cases: Vec<(String, Solid3D, (usize, usize, usize))> = vec![
        (
            "a block that conducts nothing, with a film".to_string(),
            heated(
                row(4, insulator()).losing_from(Face::XMin, film(10.0, 1e-6)),
                3,
            ),
            (3, 0, 0),
        ),
        (
            "a film of zero that does not radiate".to_string(),
            heated(row(4, grey()).losing_from(Face::XMin, film(0.0, 1e-6)), 3),
            (0, 0, 0),
        ),
        (
            "a film of no area".to_string(),
            heated(row(4, grey()).losing_from(Face::XMin, film(10.0, 0.0)), 3),
            (0, 0, 0),
        ),
        (
            "a film behind a mounting of zero".to_string(),
            heated(cooled().mounted_on(Face::XMin, 0.0), 3),
            (0, 0, 0),
        ),
        (
            "a film on a cell that conducts nothing".to_string(),
            heated(
                {
                    let mut b = cooled();
                    b.fill(insulator(), |i, _, _| i == 0);
                    b
                },
                3,
            ),
            (1, 0, 0),
        ),
        // The film on the heated cell itself, which conducts nothing: it reaches the film and
        // nothing behind the surface reaches it.
        (
            "a heated cell that conducts nothing, under the film".to_string(),
            heated(
                {
                    let mut b = cooled();
                    b.fill(insulator(), |i, _, _| i == 0);
                    b
                },
                0,
            ),
            (0, 0, 0),
        ),
        // Facing the cooled body across a radiating vacuum, from a cell that conducts nothing:
        // the pair has no half cell at that end, and carries nothing.
        (
            "a heated cell that conducts nothing, across a vacuum".to_string(),
            heated(
                {
                    let mut b = row(5, black())
                        .empty(|i, _, _| i == 2)
                        .losing_from(Face::XMin, film(10.0, 1e-6));
                    b.fill(stuff(1000.0, 0.0, 0.9), |i, _, _| i == 3);
                    b
                },
                3,
            ),
            (3, 0, 0),
        ),
        (
            "air of zero in the void".to_string(),
            heated(
                row(5, grey())
                    .empty(|i, _, _| i == 4)
                    .air_in(Temperature::celsius(0.0), 0.0),
                3,
            ),
            (0, 0, 0),
        ),
        (
            "across a vacuum from the cooled body".to_string(),
            heated(
                row(5, grey())
                    .empty(|i, _, _| i == 2)
                    .losing_from(Face::XMin, film(10.0, 1e-6)),
                4,
            ),
            (3, 0, 0),
        ),
        (
            "bodies that see only each other".to_string(),
            heated(row(5, black()).empty(|i, _, _| i == 2), 4),
            (0, 0, 0),
        ),
    ];
    // Along each axis: behind a cell that conducts nothing, and behind a joint of no contact.
    // The heat in the `p`th cell along the line, not the `p`th along x, which a line along y or z
    // does not have.
    let heated_along =
        |b: Solid3D, a: usize, p: usize| b.dissipating(w, move |i, j, k| [i, j, k][a] == p);
    for a in 0..3 {
        let mut b = line(a, 5, grey()).losing_from(LOW[a], film(10.0, 1e-6));
        b.fill(insulator(), move |i, j, k| [i, j, k][a] == 2);
        cases.push((
            format!("behind a cell that conducts nothing, along {}", NAMES[a]),
            heated_along(b, a, 4),
            at(a, 3),
        ));
        let axis = AXES[a];
        cases.push((
            format!("behind a joint of no contact, along {}", NAMES[a]),
            heated_along(
                line(a, 4, grey())
                    .joined(move |face, i, j, k| (face == axis && [i, j, k][a] == 2).then_some(0.0))
                    .losing_from(LOW[a], film(10.0, 1e-6)),
                a,
                3,
            ),
            at(a, 2),
        ));
    }
    for (why, b, cell) in &cases {
        warms(why, &refused(why, b), *cell, w);
    }

    // Both bodies heated, and only the second sealed: the second is named, with its own watts.
    let b = row(5, grey())
        .empty(|i, _, _| i == 2)
        .losing_from(Face::XMin, film(10.0, 1e-6))
        .dissipating(0.25, |i, _, _| i == 0)
        .dissipating(0.5, |i, _, _| i == 4);
    warms(
        "both heated, the second sealed",
        &refused("the second sealed", &b),
        (3, 0, 0),
        0.5,
    );
    // Two sealed parts heated, after one that is not: the first of them is named.
    let b = row(8, grey())
        .empty(|i, _, _| i == 2 || i == 5)
        .dissipating(0.2, |i, _, _| i == 4)
        .dissipating(0.3, |i, _, _| i == 7);
    warms(
        "two sealed parts heated, the second and third",
        &refused("two sealed", &b),
        (3, 0, 0),
        0.2,
    );

    // What carries the heat out still does, each the counterpart of a case above.
    let carry: Vec<(&str, Solid3D)> = vec![
        (
            "a film of zero that radiates",
            heated(
                row(4, stuff(1000.0, 50.0, 0.5)).losing_from(Face::XMin, film(0.0, 1e-6)),
                3,
            ),
        ),
        (
            "a film behind a mounting",
            heated(cooled().mounted_on(Face::XMin, 100.0), 3),
        ),
        (
            "air in the void",
            heated(
                row(5, grey())
                    .empty(|i, _, _| i == 4)
                    .air_in(Temperature::celsius(0.0), 10.0),
                3,
            ),
        ),
        (
            "across a radiating vacuum from the cooled body",
            heated(
                row(5, black())
                    .empty(|i, _, _| i == 2)
                    .losing_from(Face::XMin, film(10.0, 1e-6)),
                4,
            ),
        ),
        (
            "a film without limit",
            heated(
                row(4, grey()).losing_from(Face::XMin, film(f64::INFINITY, 1e-6)),
                3,
            ),
        ),
    ];
    for (why, b) in &carry {
        let t = solved(why, b);
        let (nx, _, _) = b.counts();
        let coldest = (0..nx)
            .filter(|i| !b.is_void(*i, 0, 0))
            .map(|i| t[i])
            .fold(f64::INFINITY, f64::min);
        println!("  {why:56} solved, coldest {:.4} C", coldest - 273.15);
        assert!(
            coldest > ROOM,
            "{why}: every cell of a heated block sits above its room, not at {coldest} K"
        );
    }
}

/// **A balance below absolute zero is refused**: a sink that takes more than can reach it has no
/// steady state, in a part that loses nothing and in one that does, and the cell named is the
/// sink's.
#[test]
fn a_balance_below_absolute_zero_is_refused() {
    let below = |why: &str, b: &Solid3D, cell: &str| -> f64 {
        let e = refused(why, b);
        println!("  {why:56} -> {e}");
        // A balance the sweeps settled on, said without the hedge the mirror needed: "or more than
        // reaches them in time for this solve" ended it, and "where 2000 sweeps had not settled"
        // ends one they ran out before.
        assert!(
            e.quantity
                .contains("no balance above absolute zero was found")
                && e.quantity.ends_with(
                    "drawn down by the sinks in its part, which take more than can reach them"
                )
                && e.quantity.contains(cell),
            "{why}: refused for the wrong reason: {e}"
        );
        e.before
    };

    // Two bodies that see only each other, a hundredth of a watt in one and out of the other.
    // Their heat sets a mean of 40 °C, and a pair carrying 0.01 W needs `u⁴ − v⁴` of 2.2e11 K⁴,
    // which no two temperatures above zero averaging 313 K have. **Refused at the balance's one
    // root**, each fourth power read as `T·|T|³`: four equal cells holding 4 × 313.15 K, each body
    // 0.2 K across from its source or sink to the cell beside the clearance and each surface 0.1 K
    // from its cell, so the two surfaces sum to 626.3 K and `C(u⁴ + v⁴) = 0.01` with `v` below
    // zero; the sink's cell is `v − 0.3`. Read even, `u⁴ − v⁴`, it was 0.0145 K colder.
    let mut b = row(5, black())
        .empty(|i, _, _| i == 2)
        .dissipating(-0.01, |i, _, _| i == 0)
        .dissipating(0.01, |i, _, _| i == 4);
    set(&mut b, 0, 100.0);
    let at = below("a sink that radiation cannot feed", &b, "(0, 0, 0)");
    let c = STEFAN_BOLTZMANN.to_si() * 1e-6 / (1.0 / 0.9 + 1.0 / 0.9 - 1.0);
    let sum = 2.0 * 313.15;
    let fourth = |t: f64| (t * t) * (t * t);
    // `C(u⁴ + v⁴)` with `u = sum − v` falls as `v` rises from `−sum/2` towards zero.
    let (mut lo, mut hi) = (-0.5 * sum, 0.0);
    for _ in 0..200 {
        let v = 0.5 * (lo + hi);
        if c * (fourth(sum - v) + fourth(v)) > 0.01 {
            lo = v;
        } else {
            hi = v;
        }
    }
    sits(
        "the sink's cell, at the root below zero",
        &[at],
        &[0],
        0.5 * (lo + hi) - 0.3,
    );

    // A cooled bar whose sink, at the far end from its film, takes a hundredth of a watt where the
    // film brings `hA·T∞` = 2.7e-3 W at most: the sink's cell at
    // `T∞ − P·(1/hA + dx/2kA + 3·dx/kA)` = 273.15 − 0.01·(1e5 + 10 + 60) = −727.55 K.
    let b = row(4, grey())
        .losing_from(Face::XMin, film(10.0, 1e-6))
        .dissipating(-0.01, |i, _, _| i == 3);
    let at = below("a sink a film cannot feed", &b, "(3, 0, 0)");
    let expect = ROOM - 0.01 * (1e5 + 10.0 + 60.0);
    sits("the sink's cell", &[at], &[0], expect);

    // A bar held at a room at absolute zero by films without limit at both ends is answered there,
    // to the solve's own error, and not refused: it stopped at −8.8e-7 K and was called diverged.
    let held = || Environment {
        ambient: Temperature::from_si(0.0),
        convection_w_per_m2_k: f64::INFINITY,
        area: Area::from_si(1e-6),
    };
    let b = row(8, grey())
        .losing_from(Face::XMin, held())
        .losing_from(Face::XMax, held());
    let t = solved("held at absolute zero", &b);
    sits(
        "held at absolute zero",
        &t,
        &(0..8).collect::<Vec<_>>(),
        0.0,
    );
    assert!(
        t.iter().all(|v| *v >= 0.0),
        "a field at absolute zero is answered at or above it: {t:?}"
    );

    // A sink the film can feed sits below the room and above zero.
    let b = row(4, grey())
        .losing_from(Face::XMax, film(10.0, 1e-6))
        .dissipating(-0.001, |i, _, _| i == 0);
    let t = solved("a sink the film can feed", &b);
    sits(
        "a sink the film can feed, at its own cell",
        &t,
        &[0],
        ROOM - 0.001 * (1e5 + 10.0 + 60.0),
    );
}

/// **More bodies in one part than the step on their levels takes is refused, by count**, rather than
/// left without the step and so without a level: 257 bodies of one cell in a row, each seeing its
/// neighbours across a cell of vacuum, are one part; 256 are solved.
#[test]
fn a_part_of_more_bodies_than_the_solve_couples_is_refused() {
    let bodies = |n: usize| {
        let mut b = row(2 * n - 1, black()).empty(|i, _, _| i % 2 == 1);
        set(&mut b, 0, 100.0);
        b
    };
    let e = refused("257 bodies", &bodies(257));
    println!("  257 bodies -> {e}");
    assert!(
        e.quantity.contains("257 bodies") && e.quantity.contains("256 at most"),
        "refused for the wrong reason: {e}"
    );
    assert_eq!(e.before, 257.0, "the refusal must carry the count");
    let b = bodies(256);
    let all: Vec<usize> = (0..511).filter(|i| i % 2 == 0).collect();
    let expect = level_of(&b, &field(&b), &all);
    let t = solved("256 bodies", &b);
    holds("256 bodies of one cell in a row", &b, &t, &all, expect);
}
