//! The air in a void cools a cell at its surface, through the half cell behind it.
//!
//! A finite-volume cell knows its centre, and the air touches its face, half a cell away. For an
//! outer face that distance has been a resistance in series since the cooled boundary was made
//! second order — `1/(1/(hA) + dx/(2kA))` — and the air in a void was charged `h·dx²` a face
//! against the centre, a half cell's drop from the surface's temperature. **First order**: in a
//! steady slab every cell facing the air sat too cold by exactly that drop, `q''·dx/2k`, and the
//! conservation audit was happy, because the heat really did leave.
//!
//! Measured on the slab below, 16 mm deep, `k = 0.6`, heated to shed 100 kW/m² into air at
//! `h = 200`: every cell 333.3 K too cold at 4 mm a cell and 41.67 K at 0.5 mm, halving with the
//! cell, where the same slab cooled through an outer face met the closed form to 1e-13 K at its
//! face and to the steady solve's stop everywhere.
//!
//! # What is checked against what
//!
//! **The slab, steady.** Uniform heating `q'''` over a depth `L` in a column of cells, insulated at
//! one end and facing a void full of air at the other. All of `q'''L` leaves through the air, so the
//! surface is at `T_air + q'''L/h` and the inside is the parabola `T_s + q'''(L² − x²)/2k`; a
//! cell-centred scheme whose boundary is exact sits **`q'''dx²/8k` above that parabola at every
//! cell**, which is the closed form, not a limit to approach. An outer face convecting at the same
//! `h` is the control, and the void is put on each of the six sides of the column, because the air
//! is counted from the grid's faces in each direction separately.
//!
//! **One step of one cell.** A cell whose only neighbours are void loses `dt·G·(T − T_air)` in a
//! step, with `G` the air's film in series with its half cell, a face at a time; and its limit is
//! its capacity over that `G`. Of a poor conductor; of a cell with air on two faces; of aluminium;
//! and of ice with water named as its liquid, three ways, because the half cell is the phase's as
//! the operator was derived: water, whose half cell is the liquid's and not the ice's; ice, whose
//! half cell is the solid's and not the water's; and a cell half melted taking heat off the bus,
//! whose half cell is the half-melted one the step began with — the heat goes in after the sweep.

use pantometry_core::substance::ThermalProps;
use pantometry_core::units::{
    Area, Density, Length, SpecificHeat, Temperature, ThermalConductivity, ThermalExpansion, Time,
};
use pantometry_core::{Domain, Exchange, Substance};
use pantometry_thermal::{Environment, Face, Solid3D, HEAT};

/// Density, kg/m³, and specific heat, J/kg·K.
const RHO: f64 = 1000.0;
const C: f64 = 4000.0;
/// Ice's conductivity and specific heat, as `Substance::ice` states them, and the liquid this test
/// names for it.
const ICE: (f64, f64) = (2.22, 2050.0);
const WATER: (f64, f64) = (0.598, 4182.0);

/// A substance that conducts `k` and does not radiate, so the air is the only thing at the
/// surface and a clearance's pairing has nothing to pair.
fn stuff(k: f64) -> Substance {
    Substance::bulk("stuff", Density::from_si(RHO)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(k),
        specific_heat: SpecificHeat::j_per_kg_k(C),
        expansion: ThermalExpansion::ppm_per_k(50.0),
        emissivity: 0.0,
    })
}

/// Ice with water named as its liquid, neither radiating.
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

/// `phi` of the way from the solid's value to the liquid's.
fn between(solid: f64, liquid: f64, phi: f64) -> f64 {
    solid + phi * (liquid - solid)
}

/// Whether a residual is out: at or past `tolerance`, or not a number.
fn off(residual: f64, tolerance: f64) -> bool {
    residual.is_nan() || residual >= tolerance
}

/// The air's film on one face of area `a`, in series with a half cell of `k`, W/K.
fn through_half_cell(h: f64, k: f64, dx: f64) -> f64 {
    let a = dx * dx;
    1.0 / (1.0 / (h * a) + dx / (2.0 * k * a))
}

/// Where the column's surface is: past the air, on which side of the column, or an outer face.
#[derive(Clone, Copy)]
enum Facing {
    /// The control: the column's own end, exposed to the room.
    Room,
    /// A void along `axis` (0, 1, 2 for x, y, z), past the column's high end or its low one.
    Air { axis: usize, high: bool },
}

/// For each cell of a column, its distance in cells from the insulated end, or `None` for void.
type Along = Box<dyn Fn(usize, usize, usize) -> Option<usize>>;

struct Slab {
    name: &'static str,
    k: f64,
    h: f64,
    facing: Facing,
    air: f64,
    /// Where the surface settles, K: the heating is what makes it.
    settles_at: f64,
}

/// How far a cell of the slab may be from the closed form, K. The steady solve stops when a
/// sweep moves the field by less than a microkelvin, and what that leaves is measured at 7.2e-6 K
/// at most here; against the air charged at the centre, which missed by 41 K and more on the
/// slabs of `k = 0.6`, 16.7 K on the one the air heats and 0.025 K on the conductor, at 32 cells.
const CELL: f64 = 1e-4;

/// **A slab facing air meets the cell-centred scheme's own closed form at every cell**, as one
/// facing the room through an outer face always has.
#[test]
fn a_slab_facing_air_sits_where_one_facing_a_room_does() {
    const DEPTH: f64 = 0.016;
    let poor = Slab {
        name: "",
        k: 0.6,
        h: 200.0,
        facing: Facing::Room,
        air: 300.0,
        settles_at: 800.0,
    };
    let past = |axis: usize, high: bool| Facing::Air { axis, high };
    let slabs = [
        Slab {
            name: "an outer face (the control)",
            ..poor
        },
        Slab {
            name: "air past x-max",
            facing: past(0, true),
            ..poor
        },
        Slab {
            name: "air past x-min",
            facing: past(0, false),
            ..poor
        },
        Slab {
            name: "air past y-max",
            facing: past(1, true),
            ..poor
        },
        Slab {
            name: "air past y-min",
            facing: past(1, false),
            ..poor
        },
        Slab {
            name: "air past z-max",
            facing: past(2, true),
            ..poor
        },
        Slab {
            name: "air past z-min",
            facing: past(2, false),
            ..poor
        },
        // A sink inside and air hotter than the slab: the heat comes in through the surface.
        Slab {
            name: "air hotter than the slab",
            facing: past(0, true),
            air: 600.0,
            settles_at: 400.0,
            ..poor
        },
        // A conductor, whose half cell is one part in 250 of its film's resistance at 4 mm.
        Slab {
            name: "a conductor",
            k: 200.0,
            h: 400.0,
            facing: past(0, true),
            settles_at: 350.0,
            ..poor
        },
    ];

    let mut wrong = Vec::new();
    for Slab {
        name,
        k,
        h,
        facing,
        air,
        settles_at,
    } in slabs
    {
        // What the surface sheds, W/m², and so the heating, W/m³, that sheds it. The block is told
        // the heating and not the temperature: it starts at the surface's, and a solve that moved
        // nothing would miss the parabola by the whole rise across the slab.
        let q = h * (settles_at - air);
        let heating = q / DEPTH;
        let mut previous = f64::NAN;
        let mut orders = 0;
        for n in [4usize, 8, 16, 32] {
            let dx = DEPTH / n as f64;
            // The column, and for each solid cell its distance from the insulated end, which is
            // what the parabola is written in. Twice as long where the void is part of it.
            let (block, along): (Solid3D, Along) = match facing {
                Facing::Room => {
                    let block = Solid3D::new(
                        "slab",
                        stuff(k),
                        (n, 1, 1),
                        Length::from_si(dx),
                        Temperature::from_si(settles_at),
                    )
                    .losing_from(
                        Face::XMax,
                        Environment {
                            ambient: Temperature::from_si(air),
                            convection_w_per_m2_k: h,
                            area: Area::from_si(dx * dx),
                        },
                    );
                    (block, Box::new(|i, _, _| Some(i)))
                }
                Facing::Air { axis, high } => {
                    let mut counts = [1usize; 3];
                    counts[axis] = 2 * n;
                    // Distance from the insulated end, in cells, or `None` for the void.
                    let along = move |i: usize, j: usize, k: usize| {
                        let p = [i, j, k][axis];
                        match (high, p < n) {
                            (true, true) => Some(p),
                            (false, false) => Some(2 * n - 1 - p),
                            _ => None,
                        }
                    };
                    let block = Solid3D::new(
                        "slab",
                        stuff(k),
                        (counts[0], counts[1], counts[2]),
                        Length::from_si(dx),
                        Temperature::from_si(settles_at),
                    )
                    .empty(move |i, j, k| along(i, j, k).is_none())
                    .air_in(Temperature::from_si(air), h);
                    (block, Box::new(along))
                }
            };
            let block = block.dissipating(q * dx * dx, |_, _, _| true);
            let (nx, ny, nz) = block.counts();
            // The block's own count of its void, not this test's map of it.
            if nx * ny * nz - block.void_cells() != n {
                wrong.push(format!(
                    "{name}, {n} cells: the block has {} solid cells",
                    nx * ny * nz - block.void_cells()
                ));
                continue;
            }
            let t = match block.steady_state() {
                Ok(t) => t,
                Err(v) => {
                    wrong.push(format!("{name}, {n} cells: no steady state: {v:?}"));
                    continue;
                }
            };
            let offset = heating * dx * dx / (8.0 * k);
            let parabola = |x: f64| settles_at + heating / (2.0 * k) * (DEPTH * DEPTH - x * x);
            let mut worst: f64 = 0.0;
            let mut last = f64::NAN;
            for kk in 0..nz {
                for j in 0..ny {
                    for i in 0..nx {
                        let Some(p) = along(i, j, kk) else { continue };
                        let x = (p as f64 + 0.5) * dx;
                        let got = t[i + nx * (j + ny * kk)];
                        let miss = got - (parabola(x) + offset);
                        if miss.is_nan() || miss.abs() > worst.abs() {
                            worst = miss;
                        }
                        if p == n - 1 {
                            last = got - parabola(x);
                        }
                    }
                }
            }
            if off(worst.abs(), CELL) {
                wrong.push(format!(
                    "{name}, {n} cells: a cell is {worst:+.3e} K from the parabola plus \
                     q'''dx²/8k = {offset:.3e} K"
                ));
            }
            // The order, at the cell beside the air: its distance from the continuous parabola is
            // `q'''dx²/8k` exactly when the surface is right, so it falls by four a halving, and
            // charged against the centre it fell by two. Held to what cells within `CELL` of the
            // closed form can say about a ratio and no closer — `5·CELL` over the error, to first
            // order — so it says nothing the cells have not, and says it as an order.
            let ratio = previous / last;
            println!(
                "  {name}, {n} cells: the face's cell {last:+.6e} K off, {ratio:.4} a halving"
            );
            if previous.is_finite() {
                orders += 1;
                let allowed = 5.0 * CELL / (last.abs() - CELL);
                if off((ratio - 4.0).abs(), allowed.abs()) {
                    wrong.push(format!(
                        "{name}, {n} cells: the error fell by {ratio:.4} a halving, not 4 within \
                         {allowed:.2e}"
                    ));
                }
            }
            previous = last;
        }
        if orders != 3 {
            wrong.push(format!("{name}: {orders} halvings were measured, not 3"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

struct Cell {
    name: &'static str,
    block: Solid3D,
    /// The cell, and how many of its faces touch the void.
    at: (usize, usize, usize),
    faces: usize,
    /// What its half cell conducts and what it holds, W/m·K and J/K, as the step begins.
    k: f64,
    capacity: f64,
    h: f64,
    dx: f64,
    from: f64,
    air: f64,
    /// Heat on the bus for the step, J, and the latent heat a mushy cell takes it into, J.
    bus: f64,
    latent: f64,
    /// A conductivity this case exists to tell the right one from, and whose it is.
    instead: Option<(f64, &'static str)>,
}

/// What a step may lose against `dt·G·(T − T_air)`, and the limit miss `C/G`, relative. Each side
/// is a dozen `+ − × ÷` in different orders: measured, 4.4e-16 at most for the limit and 2.2e-16
/// for the loss, so 1e-14 is twenty times that, against 6.0e-5 for the nearest wrong reading here,
/// aluminium's film charged against its centre.
const LOSS: f64 = 1e-14;

/// **A step loses what the air takes at the surface**, and the limit is the capacity over that.
#[test]
fn a_step_loses_what_the_air_takes_at_the_surface() {
    let alone = |k: f64, dx: f64, h: f64, from: f64, air: f64| -> Solid3D {
        Solid3D::new(
            "cell",
            stuff(k),
            (2, 1, 1),
            Length::from_si(dx),
            Temperature::from_si(from),
        )
        .empty(|i, _, _| i == 1)
        .air_in(Temperature::from_si(air), h)
    };
    // A millimetre cell of ice and water beside a millimetre of void with air at `air_c`.
    let frozen = |from_c: f64, air_c: f64| -> Solid3D {
        Solid3D::new(
            "cell",
            ice_and_water(),
            (2, 1, 1),
            Length::mm(1.0),
            Temperature::celsius(from_c),
        )
        .empty(|i, _, _| i == 1)
        .air_in(Temperature::celsius(air_c), 1e4)
    };
    let mass = Substance::ice().density.to_si() * 1e-9;
    let latent = mass
        * Substance::ice()
            .fusion
            .expect("ice melts")
            .latent_heat
            .to_si();
    let mut fed = frozen(0.0, -10.0);
    fed.set_melted_fraction(0, 0, 0, 0.5);
    let kelvin = |c: f64| Temperature::celsius(c).to_si();
    let simple = |name, block, faces, k: f64, dx: f64, h, from, air| Cell {
        name,
        block,
        at: if faces == 2 { (1, 0, 0) } else { (0, 0, 0) },
        faces,
        k,
        capacity: RHO * C * dx * dx * dx,
        h,
        dx,
        from,
        air,
        bus: 0.0,
        latent: 0.0,
        instead: None,
    };
    let cells = [
        simple(
            "a poor conductor",
            alone(0.6, 1e-3, 200.0, 900.0, 300.0),
            1,
            0.6,
            1e-3,
            200.0,
            900.0,
            300.0,
        ),
        simple(
            "a cell with air on two faces",
            Solid3D::new(
                "cell",
                stuff(0.6),
                (3, 1, 1),
                Length::from_si(1e-3),
                Temperature::from_si(900.0),
            )
            .empty(|i, _, _| i != 1)
            .air_in(Temperature::from_si(300.0), 200.0),
            2,
            0.6,
            1e-3,
            200.0,
            900.0,
            300.0,
        ),
        simple(
            "aluminium",
            alone(167.0, 2e-3, 10.0, 400.0, 300.0),
            1,
            167.0,
            2e-3,
            10.0,
            400.0,
            300.0,
        ),
        Cell {
            name: "a cell of water melted out of ice",
            block: frozen(5.0, 0.0),
            at: (0, 0, 0),
            faces: 1,
            k: WATER.0,
            capacity: mass * WATER.1,
            h: 1e4,
            dx: 1e-3,
            from: kelvin(5.0),
            air: kelvin(0.0),
            bus: 0.0,
            latent,
            instead: Some((ICE.0, "ice's half cell")),
        },
        Cell {
            name: "a cell of ice",
            block: frozen(-5.0, -10.0),
            at: (0, 0, 0),
            faces: 1,
            k: ICE.0,
            capacity: mass * ICE.1,
            h: 1e4,
            dx: 1e-3,
            from: kelvin(-5.0),
            air: kelvin(-10.0),
            bus: 0.0,
            latent,
            instead: Some((WATER.0, "water's half cell")),
        },
        // The bus heat melts it from a half to nine tenths in the step, after the air is asked, and
        // the air is asked about the cell the step began with — as its faces and its limit are.
        Cell {
            name: "a half-melted cell taking heat off the bus",
            block: fed,
            at: (0, 0, 0),
            faces: 1,
            k: between(ICE.0, WATER.0, 0.5),
            capacity: mass * between(ICE.1, WATER.1, 0.5),
            h: 1e4,
            dx: 1e-3,
            from: kelvin(0.0),
            air: kelvin(-10.0),
            bus: 0.4 * latent,
            latent,
            // Once the half cell of the phase the bus heat left, which it reached before the air
            // did; it goes in after the sweep now, and `heat_off_the_bus_arrives_after_the_sweep.rs`
            // holds the phase it leaves to the next step.
            instead: None,
        },
    ];

    let mut wrong = Vec::new();
    for Cell {
        name,
        block,
        at,
        faces,
        k,
        capacity,
        h,
        dx,
        from,
        air,
        bus,
        latent,
        instead,
    } in cells
    {
        let g = faces as f64 * through_half_cell(h, k, dx);
        let want = capacity / g;
        let limit = block.max_stable_dt(Time::ZERO).to_si();
        let missed = (limit / want - 1.0).abs();
        if off(missed, LOSS) {
            wrong.push(format!(
                "{name}: the limit is C/G = {want:.6e} s, and the block reports {limit:.6e} s"
            ));
        }
        // A quarter of the smaller of the two limits, so a block whose own limit is wrong still
        // takes the step and shows the heat it moved, rather than only refusing it.
        let dt = 0.25 * if limit > 0.0 { want.min(limit) } else { want };
        let shed = dt * g * (from - air);

        // **What this case exists to tell apart**, required to differ from the right one by a
        // hundred thousand times what the loss is held to, so that a case simplified until it
        // could not see its defect says so: the film charged against the centre, every face's film
        // behind one half cell, and the conductivity it is not.
        let mut others = vec![(
            "the centre's reading",
            dt * faces as f64 * h * dx * dx * (from - air),
        )];
        if faces > 1 {
            let one = 1.0 / (1.0 / (faces as f64 * h * dx * dx) + dx / (2.0 * k * dx * dx));
            others.push(("every face behind one half cell", dt * one * (from - air)));
        }
        if let Some((other, what)) = instead {
            others.push((
                what,
                dt * faces as f64 * through_half_cell(h, other, dx) * (from - air),
            ));
        }
        for (what, other) in others {
            let apart = (other / shed - 1.0).abs();
            if apart.is_nan() || apart <= 1e5 * LOSS {
                wrong.push(format!(
                    "{name}: {what} sheds {other:e} J and the surface {shed:e} J, too close for \
                     this case to tell them apart"
                ));
            }
        }

        let phi = block.melted_fraction_at(at.0, at.1, at.2);
        let mut stepped = block.clone();
        let mut exchange = Exchange::new();
        if bus != 0.0 {
            exchange.publish(HEAT, bus);
        }
        if let Err(v) = stepped.step(Time::ZERO, Time::from_si(dt), &mut exchange) {
            wrong.push(format!("{name}: refused a quarter of its own limit: {v:?}"));
            continue;
        }
        let lost = stepped.lost_energy().to_si();
        let loss = (lost / shed - 1.0).abs();
        if off(loss, LOSS) {
            wrong.push(format!(
                "{name}: a step of {dt:.4e} s loses dt·G·(T − T_air) = {shed:e} J, and the block \
                 lost {lost:e} J"
            ));
        }
        // An ulp at 900 K is 1.1e-13 K, and every case here measured none; a fraction of order one
        // through an enthalpy rounds near 1e-16, and measured 1.1e-16.
        if phi == 0.0 || phi == 1.0 {
            let t = stepped.temperature_at(at.0, at.1, at.2).to_si();
            let expected = from - shed / capacity;
            println!(
                "  {name}: limit {missed:.1e}, loss {loss:.1e}, temperature {:.1e} K",
                (t - expected).abs()
            );
            if off((t - expected).abs(), 1e-11) {
                wrong.push(format!(
                    "{name}: one step on it should be at {expected:.9} K, and is at {t:.9} K"
                ));
            }
        } else {
            // A mushy cell changes phase instead, by what reached it over its latent heat.
            let melted = stepped.melted_fraction_at(at.0, at.1, at.2);
            let expected = phi + bus / latent - shed / latent;
            println!(
                "  {name}: limit {missed:.1e}, loss {loss:.1e}, fraction {:.1e}",
                (melted - expected).abs()
            );
            if off((melted - expected).abs(), 1e-12) {
                wrong.push(format!(
                    "{name}: one step on it should be {expected:.9} melted, and is {melted:.9}"
                ));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
