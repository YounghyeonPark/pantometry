//! A cell that has melted radiates as the liquid it is, not as the solid it was.
//!
//! A finish is a surface, and a surface that has melted is the liquid's: a black coat on a solid
//! does not survive it. A substance states its liquid's emissivity beside the liquid's conductivity
//! and capacity, and a cell's conductivity and capacity followed its phase — but the film on an
//! outer face and the radiation between two cells across a clearance both read the solid's
//! emissivity whatever had melted, in the heat they moved and in the limit they set alike. The two
//! agreed, so the block was stable, and the books balanced, because the heat really did move.
//!
//! Measured on the substance below, a black coat of 0.9 on a solid that melts bright, 0.1: one
//! cell of liquid at 900 K radiating to 300 K reported a limit of 30.21 s where its own is 245.2 s,
//! and shed 8.65 times what its surface sends; two cells of it across a void, at 900 K and 600 K,
//! exchanged 15.5 times their own pair's heat. The other way about, 0.1 under 0.9, it was 0.116
//! and 0.064 times.
//!
//! # What is checked against what
//!
//! Radiation alone, `h = 0`, and one path out of each cell, so every case has a closed form.
//!
//! **Through a film.** The limit is the cell's capacity over everything leaving it, `C/(K + G)`,
//! with `K` its faces to other cells and `G` the derivative `4εσT³A` in series with the half cell
//! behind the surface; and a step of `dt` loses `dt·G'·(T − T∞)`, with `G'` what the surface
//! radiates, `εσA(T_s⁴ − T∞⁴)`, over `T − T∞` — at the surface's temperature `T_s`, a half cell's
//! drop below the centre's, where the half cell brings it what it sends — leaving the cell at
//! `T − dt·G'·(T − T∞)/C`, or a cell in the mush its fraction less what it lost over its latent
//! heat. For a cell of liquid; the
//! same walled in by void at `(1, 2, 1)` of a larger block, where an index that is not the cell's
//! reads a void; two cells 30% melted between steps, radiating as the area-weighted mean of the
//! two surfaces; a cell of solid, which this did not change, and one whose liquid's emissivity is
//! not a number, which a solid cell has no business reading; a half-melted cell taking heat off
//! the bus, whose film is the phase the step began in — the heat goes in after the sweep, and once
//! went in before it and moved the phase under the film; and a cell melted between steps whose
//! liquid differs from its solid
//! in nothing but how it radiates, which alone made the block one that rebuilds as it melts.
//!
//! **Across a clearance.** Two cells facing each other along a grid line across void exchange
//! `k·(u⁴ − v⁴)`, `k = σA/(1/ε₁ + 1/ε₂ − 1)`, between their two facing surfaces, each behind its
//! own half cell of its own phase's conductivity — so a melted end's half cell is the liquid's —
//! and found here by bisection; the limit is the smaller capacity over the derivative at the
//! hotter end, `4k·T³`, which the half cells, in series with it, leave a bound; and a step moves
//! `dt` times what crosses from one to the other. For two cells of liquid; a liquid facing a solid, derived in full; an end melted between
//! steps, and one 30% melted, from each end and along each axis; both ends melted by a step, and
//! an end refilled with other matter, each of which the block derives in full; and a pair whose
//! ends radiate only once melted, which appears when both have and goes when one freezes again.
//!
//! `T⁴` is written `(T·T)·(T·T)` here; the block's film takes the factors `(T + T∞)(T² + T∞²)`
//! and its clearance `(u − v)(u + v)(u² + v²)`, and both find their surfaces by Newton where this
//! bisects. They can round apart, which is part of what the tolerances allow for.

use pantometry_core::substance::{FusionProps, ThermalProps};
use pantometry_core::units::{
    Area, Density, LatentHeat, Length, SpecificHeat, Temperature, ThermalConductivity,
    ThermalExpansion, Time,
};
use pantometry_core::{Domain, Exchange, Substance};
use pantometry_thermal::{Environment, Face, Solid3D, HEAT};

/// The Stefan–Boltzmann constant, W/m²·K⁴, as CODATA 2018 states it exactly.
const SIGMA: f64 = 5.670_374_419e-8;
/// The cell's side, m.
const DX: f64 = 1.0e-3;
/// Density, kg/m³, so a cell holds a milligram.
const RHO: f64 = 1000.0;
/// The solid's conductivity, specific heat and emissivity: a black coat.
const SOLID: (f64, f64, f64) = (2.2, 2000.0, 0.9);
/// The liquid's: it melts bright.
const LIQUID: (f64, f64, f64) = (0.6, 4000.0, 0.1);
/// Where it melts, K, and what that takes, J/kg.
const MELTS: (f64, f64) = (500.0, 3.0e5);

fn substance(solid: (f64, f64, f64), liquid: (f64, f64, f64)) -> Substance {
    let props = |(k, c, e): (f64, f64, f64)| ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(k),
        specific_heat: SpecificHeat::j_per_kg_k(c),
        expansion: ThermalExpansion::ppm_per_k(100.0),
        emissivity: e,
    };
    Substance::bulk("coated", Density::from_si(RHO))
        .with_thermal(props(solid))
        .with_fusion(
            FusionProps::new(Temperature::from_si(MELTS.0), LatentHeat::from_si(MELTS.1))
                .with_liquid(props(liquid)),
        )
}

/// A black coat on a solid that melts bright.
fn coated() -> Substance {
    substance(SOLID, LIQUID)
}

/// What the mirror's liquid radiates with.
const MIRROR_LIQUID: f64 = 0.9;

/// A solid that does not radiate at all, a perfect mirror, with `coated`'s conductivities and
/// capacities and a liquid that radiates: a clearance between two cells of it carries nothing
/// until both have melted.
fn mirror() -> Substance {
    substance((SOLID.0, SOLID.1, 0.0), (LIQUID.0, LIQUID.1, MIRROR_LIQUID))
}

/// What the painted matter radiates with.
const PAINT: f64 = 0.3;

/// Matter that conducts and holds what `coated`'s solid does, radiates at `PAINT`, and does not
/// melt.
fn painted() -> Substance {
    Substance::bulk("painted", Density::from_si(RHO)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(SOLID.0),
        specific_heat: SpecificHeat::j_per_kg_k(SOLID.1),
        expansion: ThermalExpansion::ppm_per_k(100.0),
        emissivity: PAINT,
    })
}

/// A block of `coated` at `kelvin`.
fn block(counts: (usize, usize, usize), kelvin: f64) -> Solid3D {
    block_of(coated(), counts, kelvin)
}

fn block_of(s: Substance, counts: (usize, usize, usize), kelvin: f64) -> Solid3D {
    Solid3D::new(
        "block",
        s,
        counts,
        Length::mm(1.0),
        Temperature::from_si(kelvin),
    )
}

/// `block`, radiating from `face` to `ambient` and nothing else, with one cell's area: every film
/// case here has exactly one solid cell on the face it exposes.
fn exposed(block: Solid3D, face: Face, ambient: f64) -> Solid3D {
    block.losing_from(
        face,
        Environment {
            ambient: Temperature::from_si(ambient),
            convection_w_per_m2_k: 0.0,
            area: Area::from_si(DX * DX),
        },
    )
}

fn fourth(t: f64) -> f64 {
    let s = t * t;
    s * s
}

/// A surface conductance in series with a half cell of conductivity `k`, W/K.
fn behind(g: f64, k: f64) -> f64 {
    1.0 / (1.0 / g + DX / (2.0 * k * DX * DX))
}

/// The film's derivative at `t` behind the half cell, and the conductance from a centre at `t` to
/// `ambient` through it: what the surface radiates over `t − ambient`, at the surface's own
/// temperature, where what the half cell brings it is what it sends. That temperature is found by
/// bisection, and it is below `t` by the half cell's drop — the surface and not the centre is what
/// radiates.
fn film(k: f64, emissivity: f64, t: f64, ambient: f64) -> (f64, f64) {
    let a = DX * DX;
    let hot = t.max(ambient);
    let half = 2.0 * k * a / DX;
    let sends = |s: f64| emissivity * SIGMA * a * (fourth(s) - fourth(ambient));
    let (mut lo, mut hi) = (t.min(ambient), t.max(ambient));
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        // What reaches the surface less what leaves it, which falls as the surface warms.
        if half * (t - mid) - sends(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (
        behind(4.0 * emissivity * SIGMA * a * hot * hot * hot, k),
        sends(0.5 * (lo + hi)) / (t - ambient),
    )
}

/// The parallel-plate exchange coefficient for one cell's face, W/K⁴.
fn pair(ea: f64, eb: f64) -> f64 {
    SIGMA * DX * DX / (1.0 / ea + 1.0 / eb - 1.0)
}

/// What crosses a clearance between two cells at `ta` and `tb`, W: the pair coefficient `k`
/// between their two facing surfaces, each behind its own half cell, `2k·A/dx`, of conductivity
/// `ka` and `kb`. What crosses is what each half cell carries between its centre and its surface,
/// and the surfaces send `k·(u⁴ − v⁴)` to each other; found by bisection between nothing and the
/// exchange at which the two surfaces meet, on which side of nothing the centres say.
fn crossed(k: f64, (ka, ta): (f64, f64), (kb, tb): (f64, f64)) -> f64 {
    let (ga, gb) = (2.0 * ka * DX, 2.0 * kb * DX);
    let meet = (ta - tb) / (1.0 / ga + 1.0 / gb);
    let (mut lo, mut hi) = if meet > 0.0 { (0.0, meet) } else { (meet, 0.0) };
    for _ in 0..200 {
        let q = 0.5 * (lo + hi);
        // What the surfaces send at this exchange, less the exchange, which falls as it grows.
        if k * (fourth(ta - q / ga) - fourth(tb + q / gb)) > q {
            lo = q;
        } else {
            hi = q;
        }
    }
    0.5 * (lo + hi)
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

fn mass() -> f64 {
    RHO * DX * DX * DX
}

fn latent() -> f64 {
    mass() * MELTS.1
}

struct Film {
    name: &'static str,
    block: Solid3D,
    /// The exposed cell.
    at: (usize, usize, usize),
    /// What its half cell conducts, what its surface radiates with and what it holds, at the phase
    /// the step begins in.
    k: f64,
    emissivity: f64,
    capacity: f64,
    /// The conductance of its faces to other cells, W/K.
    faces: f64,
    /// Heat on the bus for the step, J.
    bus: f64,
    from: f64,
    ambient: f64,
    /// An emissivity this case exists to tell the right one from, and what it is.
    against: Option<(f64, &'static str)>,
}

#[test]
fn a_melted_surface_radiates_through_its_film_as_its_liquid() {
    let walled = exposed(
        block((2, 3, 2), 900.0).empty(|i, j, k| (i, j, k) != (1, 2, 1)),
        Face::XMax,
        300.0,
    );
    // Not half, so a mixing law that is right at a half and nowhere else is not taken for this one.
    let mut mush = exposed(block((2, 1, 1), MELTS.0), Face::XMin, 300.0);
    mush.set_melted_fraction(0, 0, 0, 0.3);
    mush.set_melted_fraction(1, 0, 0, 0.3);
    let mut fed = exposed(block((1, 1, 1), MELTS.0), Face::XMin, 300.0);
    fed.set_melted_fraction(0, 0, 0, 0.5);
    // Its liquid conducts and holds what its solid does, and only radiates differently.
    let mut surface_only = exposed(
        block_of(
            substance(SOLID, (SOLID.0, SOLID.1, LIQUID.2)),
            (1, 1, 1),
            450.0,
        ),
        Face::XMin,
        300.0,
    );
    surface_only.set_temperature(0, 0, 0, Temperature::from_si(900.0));

    let liquid = |name, block, at| Film {
        name,
        block,
        at,
        k: LIQUID.0,
        emissivity: LIQUID.2,
        capacity: mass() * LIQUID.1,
        faces: 0.0,
        bus: 0.0,
        from: 900.0,
        ambient: 300.0,
        against: Some((SOLID.2, "the solid's")),
    };
    let half = between(SOLID.2, LIQUID.2, 0.5);
    let cases = [
        liquid(
            "a cell of liquid",
            exposed(block((1, 1, 1), 900.0), Face::XMin, 300.0),
            (0, 0, 0),
        ),
        liquid("a cell of liquid walled in by void", walled, (1, 2, 1)),
        Film {
            name: "two cells partly melted between steps",
            block: mush,
            at: (0, 0, 0),
            k: between(SOLID.0, LIQUID.0, 0.3),
            emissivity: between(SOLID.2, LIQUID.2, 0.3),
            capacity: mass() * between(SOLID.1, LIQUID.1, 0.3),
            // Two cells of one conductivity are joined by a face that conducts as either does, over
            // a cell's area and a cell's length: `k·DX`.
            faces: between(SOLID.0, LIQUID.0, 0.3) * DX,
            bus: 0.0,
            from: MELTS.0,
            ambient: 300.0,
            against: Some((SOLID.2, "the solid's")),
        },
        Film {
            name: "a cell of solid",
            block: exposed(block((1, 1, 1), 450.0), Face::XMin, 300.0),
            at: (0, 0, 0),
            k: SOLID.0,
            emissivity: SOLID.2,
            capacity: mass() * SOLID.1,
            faces: 0.0,
            bus: 0.0,
            from: 450.0,
            ambient: 300.0,
            against: None,
        },
        // A solid cell radiates as its solid whatever its liquid says, even a liquid whose emissivity
        // is not a number: it reads the liquid's nothing.
        Film {
            name: "a cell of solid whose liquid's emissivity is not a number",
            block: exposed(
                block_of(
                    substance(SOLID, (LIQUID.0, LIQUID.1, f64::NAN)),
                    (1, 1, 1),
                    450.0,
                ),
                Face::XMin,
                300.0,
            ),
            at: (0, 0, 0),
            k: SOLID.0,
            emissivity: SOLID.2,
            capacity: mass() * SOLID.1,
            faces: 0.0,
            bus: 0.0,
            from: 450.0,
            ambient: 300.0,
            against: None,
        },
        Film {
            name: "a half-melted cell taking heat off the bus",
            block: fed,
            at: (0, 0, 0),
            k: between(SOLID.0, LIQUID.0, 0.5),
            emissivity: half,
            capacity: mass() * between(SOLID.1, LIQUID.1, 0.5),
            faces: 0.0,
            bus: 0.4 * latent(),
            from: MELTS.0,
            ambient: 300.0,
            // Once the phase the bus heat left, which it reached before the film did; it goes in
            // after the sweep now, and `heat_off_the_bus_arrives_after_the_sweep.rs` holds the
            // phase it leaves to the next step.
            against: None,
        },
        Film {
            name: "a cell melted between steps, its liquid differing only in how it radiates",
            block: surface_only,
            at: (0, 0, 0),
            k: SOLID.0,
            emissivity: LIQUID.2,
            capacity: mass() * SOLID.1,
            faces: 0.0,
            bus: 0.0,
            from: 900.0,
            ambient: 300.0,
            against: Some((SOLID.2, "the solid's")),
        },
    ];

    // Collected, so a cell wrong in every way says every way, and each with its numbers.
    let mut wrong = Vec::new();
    for case in cases {
        let Film {
            name,
            block,
            at: (i, j, k),
            k: conductivity,
            emissivity,
            capacity,
            faces,
            bus,
            from,
            ambient,
            against,
        } = case;
        let phi = block.melted_fraction_at(i, j, k);
        let (tangent, secant) = film(conductivity, emissivity, from, ambient);
        // Required of each case that has one, so a case that cannot tell the right emissivity from
        // the wrong one does not pass for one that checked it.
        if let Some((other, what)) = against {
            let would = film(conductivity, other, from, ambient).1;
            assert!(
                (would / secant - 1.0).abs() > 0.25,
                "{name}: {what} surface gives {would:e} W/K and its own {secant:e}, too close for \
                 this case to tell them apart"
            );
        }

        // Each side is a dozen `+ − × ÷`, each rounding by half an ulp, in different orders, and
        // `T⁴` by a different route, and the surface found by Newton against bisection; measured,
        // 2.2e-16 at most for the limit and 6.7e-16 for the loss, and 1e-14 is fifteen times the
        // larger, against a defect of a factor.
        let limit = block.max_stable_dt(Time::from_si(0.0)).to_si();
        let want = capacity / (faces + tangent);
        if off((limit / want - 1.0).abs(), 1e-14) {
            wrong.push(format!(
                "{name}: the limit is C/(K + G) = {want:.6} s, and the block reports {limit:.6} s"
            ));
        }

        // A quarter of the closed form's limit or of the block's, whichever is shorter, so the heat a
        // step moves is checked whatever the limit says: the two readings' limits differ eight times
        // here, 245.2 s against 30.21 s, and a block reporting the shorter refuses a quarter of the
        // longer.
        // A limit that is not a positive number has already been reported, and gets the closed
        // form's, so a zero cannot make a step of nothing that loses nothing and agrees.
        let dt = 0.25 * if limit > 0.0 { limit.min(want) } else { want };
        let mut stepped = block.clone();
        let mut exchange = Exchange::new();
        if bus != 0.0 {
            exchange.publish(HEAT, bus);
        }
        if let Err(v) = stepped.step(Time::from_si(0.0), Time::from_si(dt), &mut exchange) {
            wrong.push(format!("{name}: refused a quarter of its own limit: {v:?}"));
            continue;
        }
        let lost = stepped.lost_energy().to_si();
        let shed = dt * secant * (from - ambient);
        if off((lost / shed - 1.0).abs(), 1e-14) {
            wrong.push(format!(
                "{name}: a step of {dt:.4} s loses dt·G'·(T − T∞) = {shed:e} J through its film, \
                 and the block lost {lost:e} J"
            ));
        }
        if phi == 0.0 || phi == 1.0 {
            // Through the enthalpy, which rounds at 900 K, where an ulp is 1.1e-13 K; measured,
            // exactly, and in a clearance one ulp.
            let t = stepped.temperature_at(i, j, k).to_si();
            let expected = from - shed / capacity;
            if off((t - expected).abs(), 1e-12) {
                wrong.push(format!(
                    "{name}: one step on it should be at {expected:.9} K, and is at {t:.9} K"
                ));
            }
        } else {
            // A fraction of order one through an enthalpy rounds near 1e-16; measured, exactly.
            let melted = stepped.melted_fraction_at(i, j, k);
            let expected = phi + bus / latent() - shed / latent();
            if off((melted - expected).abs(), 1e-14) {
                wrong.push(format!(
                    "{name}: one step on it should be {expected:.12} melted, and is {melted:.12}"
                ));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

struct Gap {
    name: &'static str,
    block: Solid3D,
    /// The two ends, low index first.
    ends: [(usize, usize, usize); 2],
    /// What each end radiates with, holds and conducts through its half cell, at the phase the
    /// step begins in.
    emissivity: [f64; 2],
    capacity: [f64; 2],
    conductivity: [f64; 2],
    /// An emissivity pair this case exists to tell the right one from, and what it is.
    against: ([f64; 2], &'static str),
}

#[test]
fn a_melted_surface_radiates_across_a_clearance_as_its_liquid() {
    let along = |counts: (usize, usize, usize), from: [f64; 2]| {
        let mut b = block(counts, from[0]).empty(|i, j, k| i + j + k == 1);
        let far = (counts.0 - 1, counts.1 - 1, counts.2 - 1);
        b.set_temperature(far.0, far.1, far.2, Temperature::from_si(from[1]));
        (b, far)
    };
    let (two_liquid, far) = along((3, 1, 1), [900.0, 600.0]);
    let liquids = Gap {
        name: "two cells of liquid across a void",
        block: two_liquid,
        ends: [(0, 0, 0), far],
        emissivity: [LIQUID.2; 2],
        capacity: [mass() * LIQUID.1; 2],
        conductivity: [LIQUID.0; 2],
        against: ([SOLID.2; 2], "the solid's at both ends"),
    };
    // Its far end frozen by a mutator, and then the whole block derived again by a fill of nothing,
    // so that the pair is the full derivation's: with unequal ends, which no other case gives it.
    let (mut mixed, far) = along((3, 1, 1), [900.0, 400.0]);
    mixed.fill(coated(), |_, _, _| false);
    let liquid_solid = Gap {
        name: "a liquid facing a solid across a void",
        block: mixed,
        ends: [(0, 0, 0), far],
        emissivity: [LIQUID.2, SOLID.2],
        capacity: [mass() * LIQUID.1, mass() * SOLID.1],
        conductivity: [LIQUID.0, SOLID.0],
        against: ([SOLID.2; 2], "the solid's at both ends"),
    };
    let solid_ends = [SOLID.2; 2];
    let mut cases = vec![liquids, liquid_solid];
    // An end melted between steps, so its pairs are put right by a mutator rather than by a full
    // derivation: from the low end of x, which walks up the line, and from the high end, which walks
    // down it; along y and along z; and partly melted, which leaves the end in the mush.
    for (name, counts, melt_far) in [
        (
            "an end melted between steps, from below in x",
            (3, 1, 1),
            false,
        ),
        (
            "an end melted between steps, from above in x",
            (3, 1, 1),
            true,
        ),
        (
            "an end melted between steps, from below in y",
            (1, 3, 1),
            false,
        ),
        (
            "an end melted between steps, from above in z",
            (1, 1, 3),
            true,
        ),
    ] {
        let (mut b, far) = along(counts, [400.0, 400.0]);
        let at = if melt_far { far } else { (0, 0, 0) };
        b.set_temperature(at.0, at.1, at.2, Temperature::from_si(900.0));
        let (e, c, conducts) = if melt_far {
            (
                [SOLID.2, LIQUID.2],
                [SOLID.1, LIQUID.1],
                [SOLID.0, LIQUID.0],
            )
        } else {
            (
                [LIQUID.2, SOLID.2],
                [LIQUID.1, SOLID.1],
                [LIQUID.0, SOLID.0],
            )
        };
        cases.push(Gap {
            name,
            block: b,
            ends: [(0, 0, 0), far],
            emissivity: e,
            capacity: [mass() * c[0], mass() * c[1]],
            conductivity: conducts,
            against: (solid_ends, "the pair as it was before the end melted"),
        });
    }
    let (mut b, far) = along((1, 3, 1), [400.0, 400.0]);
    b.set_melted_fraction(far.0, far.1, far.2, 0.3);
    cases.push(Gap {
        name: "an end partly melted between steps, from above in y",
        block: b,
        ends: [(0, 0, 0), far],
        emissivity: [SOLID.2, between(SOLID.2, LIQUID.2, 0.3)],
        capacity: [mass() * SOLID.1, mass() * between(SOLID.1, LIQUID.1, 0.3)],
        conductivity: [SOLID.0, between(SOLID.0, LIQUID.0, 0.3)],
        against: (solid_ends, "the pair as it was before the end melted"),
    });

    // Both ends melted by a step rather than by a mutator, so the pair the next step sees is the one
    // the step's own derivation lists, and that has to read the phase the step left. A mirror until
    // it melts, so before the step there is no pair at all; both ends a little below the melting
    // point, and heat off the bus enough to take each well into the liquid.
    let mut stepped = block_of(mirror(), (3, 1, 1), 499.9).empty(|i, _, _| i == 1);
    stepped.set_temperature(2, 0, 0, Temperature::from_si(499.0));
    let mut bus = Exchange::new();
    bus.publish(HEAT, 1.4);
    stepped
        .step(Time::from_si(0.0), Time::from_si(1.0), &mut bus)
        .expect("nothing conducts and nothing radiates yet, so any step is stable");
    assert!(
        [0, 2]
            .iter()
            .all(|i| stepped.melted_fraction_at(*i, 0, 0) == 1.0),
        "the heat off the bus did not melt both ends, so this case is not the one it says"
    );
    cases.push(Gap {
        name: "both ends melted by a step",
        block: stepped,
        ends: [(0, 0, 0), (2, 0, 0)],
        emissivity: [MIRROR_LIQUID; 2],
        capacity: [mass() * LIQUID.1; 2],
        conductivity: [LIQUID.0; 2],
        against: ([0.0; 2], "the mirror they were before the step"),
    });

    // An end given different matter by `fill`, which derives the block again in full: the pair has
    // to be the new matter's, so that derivation has to read its own emissivities.
    let (mut refilled, far) = along((3, 1, 1), [400.0, 400.0]);
    refilled.set_temperature(0, 0, 0, Temperature::from_si(450.0));
    refilled.fill(painted(), |i, _, _| i == 2);
    cases.push(Gap {
        name: "an end refilled with matter that radiates differently",
        block: refilled,
        ends: [(0, 0, 0), far],
        emissivity: [SOLID.2, PAINT],
        capacity: [mass() * SOLID.1; 2],
        conductivity: [SOLID.0; 2],
        against: (solid_ends, "the matter the end had before"),
    });

    let mut wrong = Vec::new();
    for case in cases {
        let Gap {
            name,
            block,
            ends: [a, b],
            emissivity: [ea, eb],
            capacity: [ca, cb],
            conductivity: [ka, kb],
            against: ([oa, ob], what),
        } = case;
        let k = pair(ea, eb);
        let would = pair(oa, ob);
        assert!(
            (would / k - 1.0).abs() > 0.25,
            "{name}: {what} give {would:e} W/K⁴ and its own ends {k:e}, too close for this case \
             to tell them apart"
        );
        let (ta, tb) = (
            block.temperature_at(a.0, a.1, a.2).to_si(),
            block.temperature_at(b.0, b.1, b.2).to_si(),
        );
        let (phi_a, phi_b) = (
            block.melted_fraction_at(a.0, a.1, a.2),
            block.melted_fraction_at(b.0, b.1, b.2),
        );

        let hot = ta.max(tb);
        let want = ca.min(cb) / (4.0 * k * hot * hot * hot);
        let limit = block.max_stable_dt(Time::from_si(0.0)).to_si();
        if off((limit / want - 1.0).abs(), 1e-14) {
            wrong.push(format!(
                "{name}: the limit is min(C)/(4kT³) = {want:.6} s, and the block reports {limit:.6} s"
            ));
        }

        // As for the film, a quarter of whichever limit is shorter.
        let dt = 0.25 * if limit > 0.0 { limit.min(want) } else { want };
        let mut stepped = block.clone();
        if let Err(v) = stepped.step(Time::from_si(0.0), Time::from_si(dt), &mut Exchange::new()) {
            wrong.push(format!("{name}: refused a quarter of its own limit: {v:?}"));
            continue;
        }
        let moved = dt * crossed(k, (ka, ta), (kb, tb));
        // Required of each case with an end that is not the solid, as the emissivities are above:
        // its half cells read as the solid's must move an end by a thousand times what an end is
        // held to, or the case cannot tell the phase's conductivity from the solid's.
        if ka != SOLID.0 || kb != SOLID.0 {
            let as_solid = dt * crossed(k, (SOLID.0, ta), (SOLID.0, tb));
            let apart = (moved - as_solid).abs() / ca.min(cb);
            if apart.is_nan() || apart <= 1e3 * 1e-12 {
                wrong.push(format!(
                    "{name}: the solid's half cells move an end by {apart:e} K from its own, too \
                     little for this case to tell them apart"
                ));
            }
        }
        for (end, (i, j, kk), t, phi, c, gain) in [
            ("near", a, ta, phi_a, ca, -moved),
            ("far", b, tb, phi_b, cb, moved),
        ] {
            // The tolerances, and what was measured against them, are the film's above: here 4.4e-16
            // for the limit, one ulp at 900 K, 1.1e-13 K, for a temperature, and a fraction exactly.
            if phi == 0.0 || phi == 1.0 {
                let got = stepped.temperature_at(i, j, kk).to_si();
                let expected = t + gain / c;
                if off((got - expected).abs(), 1e-12) {
                    wrong.push(format!(
                        "{name}: the {end} end should be at {expected:.9} K after a step that \
                         moves {moved:e} J, and is at {got:.9} K"
                    ));
                }
            } else {
                let got = stepped.melted_fraction_at(i, j, kk);
                let expected = phi + gain / latent();
                if off((got - expected).abs(), 1e-14) {
                    wrong.push(format!(
                        "{name}: the {end} end should be {expected:.12} melted after a step that \
                         moves {moved:e} J, and is {got:.12}"
                    ));
                }
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// **A pair that radiates only once melted appears when its ends have, and goes when one has not.**
///
/// A clearance between two surfaces either of which does not radiate at all carries nothing, so the
/// block lists no pair for it. With `mirror`, whose solid does not radiate and whose liquid does,
/// melting both ends between steps makes a pair where there was none, and freezing one end of a
/// pair takes it away — which changes the list, rather than one entry of it, and is the case a
/// mutator cannot put right in place.
#[test]
fn a_pair_that_radiates_only_once_melted_comes_and_goes_with_the_phase() {
    // The ends at different temperatures in every state below, so that a pair listed where none
    // belongs would move something.
    let mut b = block_of(mirror(), (3, 1, 1), 400.0).empty(|i, _, _| i == 1);
    b.set_temperature(0, 0, 0, Temperature::from_si(450.0));
    let mut wrong = Vec::new();
    let none = |b: &Solid3D, when: &str, wrong: &mut Vec<String>| {
        if !b.gap_patches().is_empty() {
            wrong.push(format!(
                "{when}: a pair is listed between surfaces that do not both radiate"
            ));
        }
        let limit = b.max_stable_dt(Time::from_si(0.0)).to_si();
        if limit != f64::INFINITY {
            wrong.push(format!(
                "{when}: nothing conducts and nothing radiates, and the limit is {limit} s"
            ));
        }
        // A second, through nothing: a pair listed at 0.9 and 0.9 would move an end by 0.36 K to
        // 15 K in it, across the states below. Not to the bit, because a melting cell's step goes
        // through its enthalpy even when nothing reaches it, and that round trip can round, by an
        // ulp of 1.1e-13 K at 900 K; measured, it did not.
        let mut s = b.clone();
        s.step(Time::from_si(0.0), Time::from_si(1.0), &mut Exchange::new())
            .expect("a block that moves nothing takes any step");
        for i in [0, 2] {
            let (t0, t1) = (
                b.temperature_at(i, 0, 0).to_si(),
                s.temperature_at(i, 0, 0).to_si(),
            );
            if off((t1 - t0).abs(), 1e-12) {
                wrong.push(format!(
                    "{when}: cell {i} moved from {t0} K to {t1} K with nothing to move it"
                ));
            }
        }
    };
    none(&b, "both ends solid", &mut wrong);
    b.set_temperature(0, 0, 0, Temperature::from_si(900.0));
    none(&b, "one end melted", &mut wrong);

    b.set_temperature(2, 0, 0, Temperature::from_si(600.0));
    if b.gap_patches().iter().map(|p| p.pairs).sum::<usize>() != 1 {
        wrong.push(format!(
            "both ends melted: one pair, and the block lists {:?}",
            b.gap_patches()
        ));
    }
    let k = pair(MIRROR_LIQUID, MIRROR_LIQUID);
    let c = mass() * LIQUID.1;
    let want = c / (4.0 * k * 900.0 * 900.0 * 900.0);
    let limit = b.max_stable_dt(Time::from_si(0.0)).to_si();
    if off((limit / want - 1.0).abs(), 1e-14) {
        wrong.push(format!(
            "both ends melted: the limit is C/(4kT³) = {want:.6} s, and the block reports {limit:.6} s"
        ));
    }
    let dt = 0.25 * if limit > 0.0 { limit.min(want) } else { want };
    let mut s = b.clone();
    match s.step(Time::from_si(0.0), Time::from_si(dt), &mut Exchange::new()) {
        Err(v) => wrong.push(format!("both ends melted: refused {dt} s: {v:?}")),
        Ok(()) => {
            let moved = dt * crossed(k, (LIQUID.0, 900.0), (LIQUID.0, 600.0));
            // As in the cases above, the liquid's half cells against the solid's.
            let as_solid = dt * crossed(k, (SOLID.0, 900.0), (SOLID.0, 600.0));
            if (moved - as_solid).abs() / c <= 1e3 * 1e-12 {
                wrong.push(format!(
                    "both ends melted: the solid's half cells move the hot end by {:e} K, too \
                     little to tell them apart",
                    (moved - as_solid).abs() / c
                ));
            }
            let got = s.temperature_at(0, 0, 0).to_si();
            if off((got - (900.0 - moved / c)).abs(), 1e-12) {
                wrong.push(format!(
                    "both ends melted: the hot end should be at {:.9} K, and is at {got:.9} K",
                    900.0 - moved / c
                ));
            }
        }
    }

    // And the other way, from a block built with the pair: both ends liquid, one frozen. Its own
    // block, so that the pair it starts with is the one a full derivation lists rather than one a
    // mutator made.
    let mut frozen = block_of(mirror(), (3, 1, 1), 900.0).empty(|i, _, _| i == 1);
    if frozen.gap_patches().iter().map(|p| p.pairs).sum::<usize>() != 1 {
        wrong.push(format!(
            "built liquid: one pair, and the block lists {:?}",
            frozen.gap_patches()
        ));
    }
    frozen.set_temperature(2, 0, 0, Temperature::from_si(400.0));
    none(&frozen, "one end frozen", &mut wrong);
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
