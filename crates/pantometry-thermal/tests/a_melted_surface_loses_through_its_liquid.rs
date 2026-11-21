//! A cell that has melted loses heat through the liquid it is, not through the solid it was.
//!
//! A film on an outer face acts at the surface, and a finite-volume cell knows only its centre, so
//! the half cell between them is part of the path: `G = 1/(1/(hA) + dx/(2kA))`. For a substance
//! with a liquid phase that `k` was the solid's whatever had melted — in the heat the film took out
//! and in the limit it set alike — while the faces between cells, the capacities and the rest of
//! the operator followed the phase.
//!
//! Measured on the cell of water below, under `h = 1e4`: the block reported a limit of 1.247 s where
//! its own is 3.590 s, and a step of a quarter of that took 13.80 mJ out through the film where the
//! film carries 4.79 mJ, leaving the cell at 1.40 °C rather than 3.75 °C. That is 2.88 times the
//! loss, the path through ice's half cell against the path through water's; half melted, it was
//! 1.40 times. The books balanced throughout, because the heat really did leave.
//!
//! # What is checked against what
//!
//! One exposed cell, one face and nothing radiating, so the film is that cell's only path out and
//! both of its uses have a closed form: the limit is the cell's capacity over everything leaving
//! it, `C/(K + G)` with `K` its faces to other cells, and a step of `dt` loses `dt·G·(T − T∞)` and
//! leaves the cell at `T − dt·G·(T − T∞)/C`, or a cell in the mush its fraction less what it lost
//! over its latent heat. For:
//!
//! - a cell of water;
//! - the same cell walled in by void at `(1, 2, 1)` of a larger block, where an index that is not
//!   the cell's reads a void;
//! - two cells half melted between steps, the film's `k` taken from the face between them rather
//!   than from the mixing rule, because the claim is that the film's half cell is the faces';
//! - a cell of ice, which this did not change, and which says the change did not simply swap one
//!   phase for the other;
//! - a half-melted cell taking heat off the bus in the step, which moves its phase before the
//!   sweep and after the limit: the film is the phase the step began in, as the faces and the
//!   limit are.

use pantometry_core::substance::ThermalProps;
use pantometry_core::units::{
    Area, Length, SpecificHeat, Temperature, ThermalConductivity, ThermalExpansion, Time,
};
use pantometry_core::{Domain, Exchange, Substance};
use pantometry_thermal::{Environment, Face, Solid3D, HEAT};

/// The film's coefficient, W/m²·K: high enough that the half cell is most of the path, which is
/// where the two phases' conductivities are told apart.
const H: f64 = 1.0e4;
/// The cell's side, m.
const DX: f64 = 1.0e-3;
/// Ice's conductivity and specific heat, as `Substance::ice` states them.
const ICE: (f64, f64) = (2.22, 2050.0);
/// The liquid this test names for it.
const WATER: (f64, f64) = (0.598, 4182.0);

/// Ice with water named as its liquid, neither of them radiating, so the film is convection
/// alone and its closed form has one term.
fn ice_and_water() -> Substance {
    let base = Substance::ice().with_emissivity(0.0);
    let fusion = base.fusion.expect("ice melts");
    base.with_fusion(fusion.with_liquid(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(WATER.0),
        specific_heat: SpecificHeat::j_per_kg_k(WATER.1),
        expansion: ThermalExpansion::ppm_per_k(207.0),
        emissivity: 0.0,
    }))
}

/// A block of it at `initial_c`.
fn block(counts: (usize, usize, usize), initial_c: f64) -> Solid3D {
    Solid3D::new(
        "block",
        ice_and_water(),
        counts,
        Length::mm(1.0),
        Temperature::celsius(initial_c),
    )
}

/// `block`, losing from `face` to `ambient_c` through `H`, with one cell's area: every case here
/// has exactly one solid cell on the face it exposes.
fn exposed(block: Solid3D, face: Face, ambient_c: f64) -> Solid3D {
    block.losing_from(
        face,
        Environment {
            ambient: Temperature::celsius(ambient_c),
            convection_w_per_m2_k: H,
            area: Area::from_si(DX * DX),
        },
    )
}

/// The film in series with a half cell of conductivity `k`, in W/K.
fn film(k: f64) -> f64 {
    let a = DX * DX;
    1.0 / (1.0 / (H * a) + DX / (2.0 * k * a))
}

/// `phi` of the way from the solid's value to the liquid's.
fn between(solid: f64, liquid: f64, phi: f64) -> f64 {
    solid + phi * (liquid - solid)
}

/// Whether a residual is out: at or past `tolerance`, or not a number. Written `>=` alone it let a
/// `NaN` through, because no comparison with one is true, so a step that lost a `NaN` and left a
/// mushy cell a `NaN` melted reported nothing.
fn off(residual: f64, tolerance: f64) -> bool {
    residual.is_nan() || residual >= tolerance
}

struct Case {
    name: &'static str,
    block: Solid3D,
    /// The exposed cell.
    at: (usize, usize, usize),
    /// What its half cell conducts, and what it holds, at the phase the step begins in.
    k: f64,
    capacity: f64,
    /// The conductance of its faces to other cells, W/K.
    faces: f64,
    /// Heat on the bus for the step, J.
    bus: f64,
    from_c: f64,
    ambient_c: f64,
    /// A conductivity this case exists to tell the right one from, and what it is.
    against: Option<(f64, &'static str)>,
}

#[test]
fn a_melted_cell_loses_through_its_own_half_cell() {
    let mass = Substance::ice().density.to_si() * DX * DX * DX;
    let latent = mass
        * Substance::ice()
            .fusion
            .expect("ice melts")
            .latent_heat
            .to_si();

    let walled = exposed(
        block((2, 3, 2), 5.0).empty(|i, j, k| (i, j, k) != (1, 2, 1)),
        Face::XMax,
        0.0,
    );
    let mut pair = exposed(block((2, 1, 1), 0.0), Face::XMin, -10.0);
    pair.set_melted_fraction(0, 0, 0, 0.5);
    pair.set_melted_fraction(1, 0, 0, 0.5);
    let pair_face = pair
        .face_conductance((0, 0, 0), (1, 0, 0))
        .expect("face neighbours")
        .to_si();
    let mut fed = exposed(block((1, 1, 1), 0.0), Face::XMin, -10.0);
    fed.set_melted_fraction(0, 0, 0, 0.5);

    let water = |name, block, at| Case {
        name,
        block,
        at,
        k: WATER.0,
        capacity: mass * WATER.1,
        faces: 0.0,
        bus: 0.0,
        from_c: 5.0,
        ambient_c: 0.0,
        against: Some((ICE.0, "the solid's")),
    };
    let cases = [
        water(
            "a cell of water",
            exposed(block((1, 1, 1), 5.0), Face::XMin, 0.0),
            (0, 0, 0),
        ),
        water("a cell of water walled in by void", walled, (1, 2, 1)),
        Case {
            name: "two cells half melted between steps",
            block: pair,
            at: (0, 0, 0),
            k: pair_face / DX,
            capacity: mass * between(ICE.1, WATER.1, 0.5),
            faces: pair_face,
            bus: 0.0,
            from_c: 0.0,
            ambient_c: -10.0,
            against: Some((ICE.0, "the solid's")),
        },
        Case {
            name: "a cell of ice",
            block: exposed(block((1, 1, 1), -5.0), Face::XMin, -10.0),
            at: (0, 0, 0),
            k: ICE.0,
            capacity: mass * ICE.1,
            faces: 0.0,
            bus: 0.0,
            from_c: -5.0,
            ambient_c: -10.0,
            against: None,
        },
        Case {
            name: "a half-melted cell taking heat off the bus",
            block: fed,
            at: (0, 0, 0),
            k: between(ICE.0, WATER.0, 0.5),
            capacity: mass * between(ICE.1, WATER.1, 0.5),
            faces: 0.0,
            bus: 0.4 * latent,
            from_c: 0.0,
            ambient_c: -10.0,
            against: Some((between(ICE.0, WATER.0, 0.9), "the phase after the bus heat")),
        },
    ];

    // Collected, so a cell wrong in every way says every way, and each with its numbers.
    let mut wrong = Vec::new();
    for case in cases {
        let Case {
            name,
            block,
            at: (i, j, k),
            k: conductivity,
            capacity,
            faces,
            bus,
            from_c,
            ambient_c,
            against,
        } = case;
        let phi = block.melted_fraction_at(i, j, k);
        let g = film(conductivity);
        // Required of each case that has one, so a case that cannot tell the right conductivity
        // from the wrong one does not pass for one that checked it.
        if let Some((other, what)) = against {
            let would = film(other);
            assert!(
                (would / g - 1.0).abs() > 0.25,
                "{name}: {what} half cell gives {would:e} W/K and its own {g:e}, too close for \
                 this case to tell them apart"
            );
        }

        // Each side is about ten `+ − × ÷`, each rounding by half an ulp, in different orders, so
        // they can disagree near 1e-15; measured, 4.4e-16 at most for the limit and the loss
        // alike, and 1e-14 is more than twenty times that, against a defect of a factor.
        let limit = block.max_stable_dt(Time::from_si(0.0)).to_si();
        let want = capacity / (faces + g);
        if off((limit / want - 1.0).abs(), 1e-14) {
            wrong.push(format!(
                "{name}: the limit is C/(K + G) = {want:.6} s, and the block reports {limit:.6} s"
            ));
        }

        // A quarter of the closed form's limit, so a block that reports the solid's limit still
        // takes the step and shows the heat it moved rather than only refusing.
        let dt = 0.25 * want;
        let (t0, ta) = (
            Temperature::celsius(from_c).to_si(),
            Temperature::celsius(ambient_c).to_si(),
        );
        let mut stepped = block.clone();
        let mut exchange = Exchange::new();
        if bus != 0.0 {
            exchange.publish(HEAT, bus);
        }
        if let Err(v) = stepped.step(Time::from_si(0.0), Time::from_si(dt), &mut exchange) {
            wrong.push(format!(
                "{name}: refused a quarter of its closed-form limit: {v:?}"
            ));
            continue;
        }
        let lost = stepped.lost_energy().to_si();
        let shed = dt * g * (t0 - ta);
        if off((lost / shed - 1.0).abs(), 1e-14) {
            wrong.push(format!(
                "{name}: a step of {dt:.4} s loses dt·G·(T − T∞) = {shed:e} J through its film, \
                 and the block lost {lost:e} J"
            ));
        }
        if phi == 0.0 || phi == 1.0 {
            // The enthalpy round trip rounds at 273 K, whose ulp is 5.7e-14 K; measured, exactly.
            let t = stepped.temperature_at(i, j, k).to_si();
            let expected = t0 - dt * g * (t0 - ta) / capacity;
            if off((t - expected).abs(), 1e-12) {
                wrong.push(format!(
                    "{name}: one step on it should be at {:.6} °C, and is at {:.6} °C",
                    expected - 273.15,
                    t - 273.15
                ));
            }
        } else {
            // A mushy cell changes phase instead, by what reached it over its latent heat. A
            // fraction of order one through an enthalpy rounds near 1e-16; measured, 1.1e-16.
            let melted = stepped.melted_fraction_at(i, j, k);
            let expected = phi + bus / latent - shed / latent;
            if off((melted - expected).abs(), 1e-12) {
                wrong.push(format!(
                    "{name}: one step on it should be {expected:.9} melted, and is {melted:.9}"
                ));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
