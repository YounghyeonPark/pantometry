//! A cell at the end of a clearance is given the step its whole row allows, not its largest term.
//!
//! An explicit step is stable, and keeps every cell within the temperatures the block started in,
//! when for every cell `θ = dt·ΣG/C ≤ 1`, the sum over everything that leaves it: its faces, its
//! film, its air and each clearance it radiates across. The limit added up the first three and
//! charged each clearance pair alone, so a cell's row was the largest of its terms rather than
//! their sum, and any cell that faced two clearances, or faced one and conducted, was handed a step
//! its row did not allow.
//!
//! Measured on the first case below, three 1 mm cells in a row at 900 K across two clearances: a
//! limit of 7.392 s where the middle cell's row allows 3.696 s, and forty steps at the limit it
//! reported took cells that started within 10 K of each other to 990.5 K and 708.9 K. The second,
//! two pairs of 10 mm cells across one clearance with each face conducting what the pair radiates,
//! went 179 K outside where it started in 40 steps. The books balanced throughout, because the
//! exchange is antisymmetric: it oscillated and grew, conserving energy.
//!
//! # What is checked against what
//!
//! The limit against the smallest row's closed form, `minᵢ Cᵢ/ΣGᵢ`, with each pair's `G` its
//! tangent `4kT³` at the hotter end, `k = σA/(1/ε₁ + 1/ε₂ − 1)`: for the three cells, whose middle
//! row is two pairs; the same with air in the clearances, two pairs and two faces of air; for the
//! conducting pairs, a face and a pair; for a cell with a film and a pair; for the centre of a cube
//! that faces a cell across void in each of the six directions; for a clearance between two
//! materials, whose row is the lighter end's; for two clearances one of which goes when an end
//! freezes, so the ends are grouped again from a list that changed; and for a block whose
//! clearances are filled in again, which has no pairs left to add.
//!
//! And then forty steps at the limit the block reports. A step at `θ ≤ 1` is a weighted mean of a
//! cell, what it sees and the room, so it keeps every cell within the range the block and its
//! surroundings started in — a theorem of the scheme, not a tolerance. A limit too long by enough
//! that the sharpest mode outgrows the relaxation shows as a cell outside that range: the old one
//! did, in the three cells, the conducting pairs, the cube and the air. The film's room sits 20 K
//! below the block rather than 600, because with 600 K of range the old limit's overshoot stayed
//! inside it. A limit too long by a little does not show here, and only the first check sees it:
//! in the three cells, 1.7% too long still passes the march.

use pantometry_core::substance::{FusionProps, ThermalProps};
use pantometry_core::units::{
    Area, Density, LatentHeat, Length, SpecificHeat, Temperature, ThermalConductivity,
    ThermalExpansion, Time,
};
use pantometry_core::{Domain, Exchange, Substance};
use pantometry_thermal::{Environment, Face, Solid3D};

/// The Stefan–Boltzmann constant, W/m²·K⁴, as CODATA 2018 states it exactly.
const SIGMA: f64 = 5.670_374_419e-8;
const EMISSIVITY: f64 = 0.9;
/// Density and specific heat, so a cell's capacity is `ρc·dx³`.
const RHO_C: (f64, f64) = (1000.0, 1000.0);

fn grey(k: f64) -> Substance {
    Substance::bulk("grey", Density::from_si(RHO_C.0)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(k),
        specific_heat: SpecificHeat::j_per_kg_k(RHO_C.1),
        expansion: ThermalExpansion::ppm_per_k(10.0),
        emissivity: EMISSIVITY,
    })
}

/// `grey` at 1 W/m·K, melting at 500 K, whose solid does not radiate at all and whose liquid
/// radiates as `grey` does: a clearance with a solid end carries nothing, and the block lists no
/// pair for it.
fn mirror() -> Substance {
    let props = |emissivity| ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(1.0),
        specific_heat: SpecificHeat::j_per_kg_k(RHO_C.1),
        expansion: ThermalExpansion::ppm_per_k(10.0),
        emissivity,
    };
    Substance::bulk("mirror", Density::from_si(RHO_C.0))
        .with_thermal(props(0.0))
        .with_fusion(
            FusionProps::new(Temperature::from_si(500.0), LatentHeat::from_si(3.0e5))
                .with_liquid(props(EMISSIVITY)),
        )
}

/// How much of `grey`'s capacity per cell `light` holds.
const LIGHT: f64 = 0.25;

/// `grey` at a quarter of its density, conducting what it does at 1 W/m·K.
fn light() -> Substance {
    Substance::bulk("light", Density::from_si(RHO_C.0 * LIGHT)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(1.0),
        specific_heat: SpecificHeat::j_per_kg_k(RHO_C.1),
        expansion: ThermalExpansion::ppm_per_k(10.0),
        emissivity: EMISSIVITY,
    })
}

fn cube(t: f64) -> f64 {
    t * t * t
}

/// One clearance pair's tangent conductance, W/K, with its hotter end at `hot`.
fn pair(dx: f64, hot: f64) -> f64 {
    let k = SIGMA * dx * dx / (2.0 / EMISSIVITY - 1.0);
    4.0 * k * cube(hot)
}

/// A film's conductance as the limit charges it: convection and the radiative tangent at the
/// hotter of the cell and the room, behind the half cell of conductivity `k`. W/K.
fn film(h: f64, dx: f64, k: f64, t: f64, ambient: f64) -> f64 {
    let a = dx * dx;
    let surface = h * a + 4.0 * EMISSIVITY * SIGMA * a * cube(t.max(ambient));
    1.0 / (1.0 / surface + dx / (2.0 * k * a))
}

/// Whether a residual is out: at or past `tolerance`, or not a number.
fn off(residual: f64, tolerance: f64) -> bool {
    residual.is_nan() || residual >= tolerance
}

struct Case {
    name: &'static str,
    block: Solid3D,
    /// The capacity of the cell whose row is the smallest, and that row's conductances, W/K.
    capacity: f64,
    row: f64,
    /// The temperatures every cell must stay within, K: the range the block starts in, and the
    /// room's if it has one.
    range: (f64, f64),
}

#[test]
fn a_clearance_end_is_given_the_step_its_whole_row_allows() {
    let mut cases = Vec::new();

    // Three cells in a row, 1 mm, across two clearances: the middle one radiates to both.
    let dx = 1e-3;
    let c = RHO_C.0 * RHO_C.1 * dx * dx * dx;
    let mut three = Solid3D::new(
        "three",
        grey(1.0),
        (5, 1, 1),
        Length::from_si(dx),
        Temperature::from_si(900.0),
    )
    .empty(|i, _, _| i == 1 || i == 3);
    three.set_temperature(2, 0, 0, Temperature::from_si(890.0));
    cases.push(Case {
        name: "a cell between two clearances",
        block: three.clone(),
        capacity: c,
        row: 2.0 * pair(dx, 900.0),
        range: (890.0, 900.0),
    });

    // The same with air in the clearances, at the coldest cell's temperature so the range is the
    // block's: every end of a clearance faces the air too, and the middle cell faces it on two
    // sides, each `h·dx²` behind its half cell, beside its two pairs.
    let h_air = 100.0;
    let air = 1.0 / (1.0 / (h_air * dx * dx) + dx / (2.0 * 1.0 * dx * dx));
    cases.push(Case {
        name: "a cell between two clearances full of air",
        block: three.air_in(Temperature::from_si(890.0), h_air),
        capacity: c,
        row: 2.0 * air + 2.0 * pair(dx, 900.0),
        range: (890.0, 900.0),
    });

    // A clearance between two materials, one holding a quarter of what the other does per cell:
    // the row that sets the limit is the lighter end's, its one pair over its own capacity.
    let mut two = Solid3D::new(
        "two materials",
        grey(1.0),
        (3, 1, 1),
        Length::from_si(dx),
        Temperature::from_si(900.0),
    )
    .empty(|i, _, _| i == 1);
    two.fill(light(), |i, _, _| i == 2);
    two.set_temperature(2, 0, 0, Temperature::from_si(890.0));
    cases.push(Case {
        name: "a clearance between two materials",
        block: two,
        capacity: c * LIGHT,
        row: pair(dx, 900.0),
        range: (890.0, 900.0),
    });

    // Two pairs of 10 mm cells across one clearance, each face conducting what the pair radiates
    // at 900 K, so the face and the pair are equal terms and their sum is twice either.
    let dx = 1e-2;
    let c = RHO_C.0 * RHO_C.1 * dx * dx * dx;
    let k = pair(dx, 900.0) / dx;
    let mut four = Solid3D::new(
        "four",
        grey(k),
        (5, 1, 1),
        Length::from_si(dx),
        Temperature::from_si(900.0),
    )
    .empty(|i, _, _| i == 2);
    for (i, t) in [(1, 895.0), (3, 890.0), (4, 885.0)] {
        four.set_temperature(i, 0, 0, Temperature::from_si(t));
    }
    cases.push(Case {
        name: "a clearance end that conducts",
        block: four,
        capacity: c,
        row: k * dx + pair(dx, 895.0),
        range: (885.0, 900.0),
    });

    // One cell with a film on its outer face and a clearance on its inner one, 10 mm: the row is
    // the film and the pair, and nothing conducts. The room 20 K below the block rather than at
    // 300 K, so that the march has a range an overshoot can leave: at 300 K the old limit, 1.77
    // times this row's, stayed inside 300 K to 900 K, and at 880 K it leaves 880 K to 900 K. The
    // film's tangent is taken at the hotter of the cell and the room, 900 K either way, so the row
    // is the same.
    let (h, ambient) = (20.0, 880.0);
    let mut filmed = Solid3D::new(
        "filmed",
        grey(k),
        (3, 1, 1),
        Length::from_si(dx),
        Temperature::from_si(900.0),
    )
    .empty(|i, _, _| i == 1)
    .losing_from(
        Face::XMin,
        Environment {
            ambient: Temperature::from_si(ambient),
            convection_w_per_m2_k: h,
            area: Area::from_si(dx * dx),
        },
    );
    filmed.set_temperature(2, 0, 0, Temperature::from_si(880.0));
    cases.push(Case {
        name: "a clearance end with a film",
        block: filmed,
        capacity: c,
        row: film(h, dx, k, 900.0, ambient) + pair(dx, 900.0),
        range: (ambient, 900.0),
    });

    // The centre of a 5×5×5 cube, 1 mm, whose six neighbours are void: it faces a cell across a
    // clearance in every direction, and conducts nowhere. A conductivity small enough that no
    // cell's faces come near six pairs, so the centre is the row that sets the limit.
    let dx = 1e-3;
    let c = RHO_C.0 * RHO_C.1 * dx * dx * dx;
    let mut cross = Solid3D::new(
        "cross",
        grey(1e-3),
        (5, 5, 5),
        Length::from_si(dx),
        Temperature::from_si(890.0),
    )
    .empty(|i, j, k| {
        let off: Vec<usize> = [i, j, k].into_iter().filter(|x| *x != 2).collect();
        off == [1] || off == [3]
    });
    cross.set_temperature(2, 2, 2, Temperature::from_si(900.0));
    cases.push(Case {
        name: "a cell facing six clearances",
        block: cross,
        capacity: c,
        row: 6.0 * pair(dx, 900.0),
        range: (890.0, 900.0),
    });

    // Two clearances, both of liquid ends, and then one end frozen between steps: its pair goes and
    // the other stays, so the list of pairs changes from two to one and the ends are grouped again
    // from it. The row that sets the limit is the remaining pair's, at its hotter end; the frozen
    // cell radiates nothing and conducts nowhere.
    let mut one_left = Solid3D::new(
        "one left",
        mirror(),
        (5, 1, 1),
        Length::from_si(dx),
        Temperature::from_si(900.0),
    )
    .empty(|i, _, _| i == 1 || i == 3);
    one_left.set_temperature(2, 0, 0, Temperature::from_si(600.0));
    one_left.set_temperature(0, 0, 0, Temperature::from_si(400.0));
    cases.push(Case {
        name: "a clearance that goes while another stays",
        block: one_left,
        capacity: c,
        row: pair(dx, 900.0),
        range: (400.0, 900.0),
    });

    // Two clearances filled in again, so the block has no void and no pairs: the limit is the faces'
    // alone, two to a middle cell, and the ends a block once had are not left behind with it. With
    // a film on one outer face, because without one the limit has nothing to look past its faces
    // for and never reads the ends at all. Radiation alone, which at 900 K and 1 mm adds 0.14 of a
    // face to the end cell's one, so a middle cell's two faces are still the row that sets it.
    let kf = 1.0;
    let mut refilled = Solid3D::new(
        "refilled",
        grey(kf),
        (5, 1, 1),
        Length::from_si(dx),
        Temperature::from_si(900.0),
    )
    .empty(|i, _, _| i == 1 || i == 3)
    .losing_from(
        Face::XMin,
        Environment {
            ambient: Temperature::from_si(ambient),
            convection_w_per_m2_k: 0.0,
            area: Area::from_si(dx * dx),
        },
    );
    refilled.fill(grey(kf), |i, _, _| i == 1 || i == 3);
    for (i, t) in [(1, 895.0), (2, 890.0), (3, 895.0)] {
        refilled.set_temperature(i, 0, 0, Temperature::from_si(t));
    }
    cases.push(Case {
        name: "clearances filled in again",
        block: refilled,
        capacity: c,
        row: 2.0 * kf * dx,
        range: (ambient, 900.0),
    });

    let mut wrong = Vec::new();
    for Case {
        name,
        block,
        capacity,
        row,
        range: (lo, hi),
    } in cases
    {
        // A dozen `+ − × ÷` a side in different orders, and `T³` by a different route; measured,
        // 2.2e-16 at most, and 1e-14 is what every closed-form limit in this crate is held to,
        // against a defect of a factor.
        let want = capacity / row;
        let limit = block.max_stable_dt(Time::from_si(0.0)).to_si();
        if off((limit / want - 1.0).abs(), 1e-14) {
            wrong.push(format!(
                "{name}: the limit is the smallest row's C/ΣG = {want:.6e} s, and the block reports \
                 {limit:.6e} s"
            ));
        }

        // Forty steps at whatever limit the block reports, asked again before each. A step at
        // `θ ≤ 1` is a weighted mean of the cell and what it sees, so nothing leaves the range the
        // block began in. Rounding can carry a cell an ulp past it in a step, 1.1e-13 K at 900 K,
        // so forty steps at most 4.5e-12 K; measured, never at all, and the defect this is for left
        // the range by 179 K.
        let mut marched = block.clone();
        let (nx, ny, nz) = marched.counts();
        let mut outside: f64 = 0.0;
        // Kept apart, because `max` passes over a `NaN` and the next cell would hide it.
        let mut nan = false;
        // And how far anything moved, because a march that moves nothing stays in any range.
        let mut moved: f64 = 0.0;
        for n in 0..40 {
            let dt = marched.max_stable_dt(Time::from_si(0.0));
            if let Err(v) = marched.step(Time::from_si(0.0), dt, &mut Exchange::new()) {
                wrong.push(format!("{name}: refused step {n} at its own limit: {v:?}"));
                break;
            }
            for k in 0..nz {
                for j in 0..ny {
                    for i in 0..nx {
                        if marched.is_void(i, j, k) {
                            continue;
                        }
                        let t = marched.temperature_at(i, j, k).to_si();
                        nan |= t.is_nan();
                        outside = outside.max(lo - t).max(t - hi);
                        moved = moved.max((t - block.temperature_at(i, j, k).to_si()).abs());
                    }
                }
            }
        }
        if nan || off(outside, 1e-11) {
            wrong.push(format!(
                "{name}: forty steps at its own limit took a cell {outside:.3e} K outside the \
                 {lo} K to {hi} K it started in"
            ));
        }
        // Every case starts at least 5 K from even, so a kelvin is a tenth of what it has to move.
        if moved.is_nan() || moved < 1.0 {
            wrong.push(format!(
                "{name}: forty steps moved no cell more than {moved:.3e} K, so staying in range \
                 says nothing"
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
