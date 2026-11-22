//! A clearance radiates between its two surfaces, not between its two cell centres.
//!
//! Two solid cells facing each other across void exchange `σA(T₁⁴ − T₂⁴)/(1/ε₁ + 1/ε₂ − 1)`, and
//! the temperatures that radiate are the two facing **surfaces**', each half a cell from the centre
//! the block knows. The block read them at the centres. Every surface radiates at a temperature a
//! half cell's drop from its centre's, so the exchange came out wrong by an amount proportional to
//! `dx`: **first order**, and the conservation audit was happy, because what left one cell arrived
//! in the other.
//!
//! Measured on the pair of slabs below, 16 mm each across 4 mm, `k = 0.6`, `ε = 0.9`, 5 kW/m²
//! generated in one and shed from the other's outer face: the heated slab 24.8 K too cold at 4 mm a
//! cell and 3.13 K at 0.5 mm, falling by 1.99 a halving; with air in the clearance, 2.01, and the
//! other slab wrong too, 4.70 K at 4 mm and 0.615 K at 0.5 mm, because the air's share of each face
//! depends on where its surface is.
//!
//! # What is checked against what
//!
//! **Two slabs, steady.** Uniform heating `q'''` over a depth `L` of one slab, insulated at its far
//! end; a clearance; a second slab, which conducts what crosses to an outer face that convects and
//! radiates to a room. The second slab's inner surface `T_b` fixes everything: what its outer face
//! sheds is what it conducts, `k_b(T_b − T_o)/L`; what crosses is that plus what the air takes from
//! `T_b`; the first slab's surface is where `σ'(T_a⁴ − T_b⁴)` is what crosses; and the first slab's
//! heat is what crosses plus what its own air takes. One equation in `T_b`, solved by bisection.
//! A cell-centred scheme whose surfaces are exact sits **`q'''dx²/8k` above the first slab's
//! parabola at every cell**, and **on the second slab's line**: the closed form, not a limit to
//! approach. Across a vacuum and across air, the two slabs conducting alike and one ten times the
//! other.
//!
//! **One step of one pair.** Two cells across a cell of void: each surface is where what its half
//! cell brings it is what it sends to the other and to the air, found by bisection, and the step
//! moves each cell by what its half cell carried — so with air, and ends that conduct differently,
//! what one cell gives up is not what the other takes, and the books hear the difference as the
//! air's. Across a vacuum and across air; ends that conduct differently; air hotter than both,
//! which holds each surface above its own centre, and from a conductor to an insulator holds the
//! insulator's above both — which is what the limit has to allow for; an end that conducts
//! nothing, whose surface the air holds; and a black pair at 3000 K on insulators. Then at its own
//! limit to where it is going, every step between the coolest and the hottest of the cells and the
//! air — a theorem for a step at `θ ≤ 1`, not a tolerance — and `θ ≤ 1` itself checked row by row
//! before every step, with the secant the surfaces send at, which is what the limit has to bound
//! rather than a second copy of how it bounds it.
//!
//! **Air without limit.** It holds both surfaces at its own temperature, so a block in it has a
//! steady state with nothing crossing; `∞/∞` taken for a share made both surfaces not a number.

use pantometry_core::substance::ThermalProps;
use pantometry_core::units::{
    Area, Density, Length, SpecificHeat, Temperature, ThermalConductivity, ThermalExpansion, Time,
};
use pantometry_core::{Domain, Exchange, Substance};
use pantometry_thermal::{Environment, Face, Solid3D};

/// The Stefan–Boltzmann constant, W/m²·K⁴, as CODATA 2018 states it exactly.
const SIGMA: f64 = 5.670_374_419e-8;
/// Density, kg/m³, and specific heat, J/kg·K.
const RHO: f64 = 1000.0;
const C: f64 = 4000.0;

fn stuff(k: f64, emissivity: f64) -> Substance {
    Substance::bulk("stuff", Density::from_si(RHO)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(k),
        specific_heat: SpecificHeat::j_per_kg_k(C),
        expansion: ThermalExpansion::ppm_per_k(50.0),
        emissivity,
    })
}

fn fourth(t: f64) -> f64 {
    let s = t * t;
    s * s
}

/// Whether a residual is out: at or past `tolerance`, or not a number.
fn off(residual: f64, tolerance: f64) -> bool {
    residual.is_nan() || residual >= tolerance
}

/// The root of `f` between `lo` and `hi`, where `f` rises through it.
fn rising(lo: f64, hi: f64, f: impl Fn(f64) -> f64) -> f64 {
    let (mut lo, mut hi) = (lo, hi);
    for _ in 0..300 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

/// The parallel-plate coefficient per square metre, W/m²·K⁴, for two faces of `emissivity`.
fn plates(emissivity: f64) -> f64 {
    SIGMA / (2.0 / emissivity - 1.0)
}

/// One pair's surfaces and what crosses between them, W: each end's half cell `g`, W/K, from a
/// centre at `t`; the air on each face `h`, W/K, at `air`; and the pair's coefficient `pair`,
/// W/K⁴. Each surface is where what reaches it from its centre and from the air is what it sends
/// across: `u = (g_a t_a + h·air − q)/(g_a + h)`, `v = (g_b t_b + h·air + q)/(g_b + h)`, and
/// `q = pair·(u⁴ − v⁴)` falls as `q` grows, between nothing and the `q` at which the two meet.
fn crossed(pair: f64, (ga, ta): (f64, f64), (gb, tb): (f64, f64), h: f64, air: f64) -> [f64; 3] {
    let u = |q: f64| (ga * ta + h * air - q) / (ga + h);
    let v = |q: f64| (gb * tb + h * air + q) / (gb + h);
    let meet = (u(0.0) - v(0.0)) / (1.0 / (ga + h) + 1.0 / (gb + h));
    let q = rising(meet.min(0.0), meet.max(0.0), |q| {
        q - pair * (fourth(u(q)) - fourth(v(q)))
    });
    [q, u(q), v(q)]
}

struct Slabs {
    name: &'static str,
    /// What the heated slab and the other conduct, W/m·K.
    k: [f64; 2],
    /// The clearance's air, W/m²·K, none at zero.
    h: f64,
}

/// How far a cell may be from the closed form, K. The steady solve stops when a sweep moves the
/// field by under a microkelvin, and what that leaves is measured at 4.7e-6 K at most on these
/// slabs, so 2e-5 is four times it; against the centres' reading, which missed by 1.3 K and more at
/// 32 cells, and against a solve leaning on the wrong slopes, which stops further from its answer.
const CELL: f64 = 2e-5;

/// **Two slabs across a clearance meet the cell-centred scheme's own closed form at every cell.**
#[test]
fn two_slabs_across_a_clearance_sit_where_their_surfaces_say() {
    const DEPTH: f64 = 0.016;
    let (emissivity, room, outer_h) = (0.9, 300.0, 50.0);
    let shed = 5000.0;
    let heating = shed / DEPTH;
    let slabs = [
        Slabs {
            name: "across a vacuum",
            k: [0.6, 0.6],
            h: 0.0,
        },
        Slabs {
            name: "across air",
            k: [0.6, 0.6],
            h: 20.0,
        },
        Slabs {
            name: "across a vacuum, to a slab ten times the conductor",
            k: [0.6, 6.0],
            h: 0.0,
        },
        Slabs {
            name: "across air, to a slab ten times the conductor",
            k: [0.6, 6.0],
            h: 20.0,
        },
    ];

    let mut wrong = Vec::new();
    for Slabs {
        name,
        k: [ka, kb],
        h,
    } in slabs
    {
        // The closed form, from the second slab's inner surface `tb`.
        let sends = |t: f64| outer_h * (t - room) + emissivity * SIGMA * (fourth(t) - fourth(room));
        let through_b = |tb: f64| {
            let outer = rising(room, tb, |t| sends(t) - kb * (tb - t) / DEPTH);
            kb * (tb - outer) / DEPTH
        };
        let first_surface = |tb: f64| {
            let across = through_b(tb) + h * (tb - room);
            (fourth(tb) + across / plates(emissivity)).sqrt().sqrt()
        };
        let tb = rising(room, 2000.0, |tb| {
            h * (first_surface(tb) - room) + through_b(tb) + h * (tb - room) - shed
        });
        let (ta, carried) = (first_surface(tb), through_b(tb));
        println!("  {name}: the surfaces at {ta:.4} K and {tb:.4} K, {carried:.2} W/m² through");

        let mut previous = f64::NAN;
        let mut orders = 0;
        for n in [4usize, 8, 16, 32] {
            let dx = DEPTH / n as f64;
            let m = n / 4;
            let cells = 2 * n + m;
            let block = Solid3D::new(
                "slabs",
                stuff(ka, emissivity),
                (cells, 1, 1),
                Length::from_si(dx),
                Temperature::from_si(500.0),
            )
            .empty(|i, _, _| i >= n && i < n + m);
            let mut block = block.losing_from(
                Face::XMax,
                Environment {
                    ambient: Temperature::from_si(room),
                    convection_w_per_m2_k: outer_h,
                    area: Area::from_si(dx * dx),
                },
            );
            block.fill(stuff(kb, emissivity), |i, _, _| i >= n + m);
            let block = if h > 0.0 {
                block.air_in(Temperature::from_si(room), h)
            } else {
                block
            };
            let block = block.dissipating(shed * dx * dx, |i, _, _| i < n);
            if block.void_cells() != m {
                wrong.push(format!(
                    "{name}, {n} cells: {} void cells, not {m}",
                    block.void_cells()
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
            let offset = heating * dx * dx / (8.0 * ka);
            let parabola = |x: f64| ta + heating / (2.0 * ka) * (DEPTH * DEPTH - x * x);
            let mut worst: f64 = 0.0;
            let mut last = f64::NAN;
            for (i, got) in t.iter().enumerate() {
                let want = if i < n {
                    let x = (i as f64 + 0.5) * dx;
                    if i == n - 1 {
                        last = got - parabola(x);
                    }
                    parabola(x) + offset
                } else if i >= n + m {
                    tb - carried * ((i - n - m) as f64 + 0.5) * dx / kb
                } else {
                    continue;
                };
                let miss = got - want;
                if miss.is_nan() || miss.abs() > worst.abs() {
                    worst = miss;
                }
            }
            if off(worst.abs(), CELL) {
                wrong.push(format!(
                    "{name}, {n} cells: a cell is {worst:+.3e} K from the closed form"
                ));
            }
            // The order at the heated slab's cell beside the clearance, held to what cells within
            // `CELL` of the closed form can say about a ratio — `5·CELL` over the error, to first
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

struct Pair {
    name: &'static str,
    /// What the two ends conduct, W/m·K, and where they start, K.
    k: [f64; 2],
    from: [f64; 2],
    emissivity: f64,
    /// The clearance's air, W/m²·K, none at zero, and its temperature.
    h: f64,
    air: f64,
}

/// What a step may move a cell by, against what its half cell carried, K: an ulp at 3000 K is
/// 4.5e-13 K, and the most measured is 1.1e-13 K, an ulp at 600 K to 900 K.
const KELVIN: f64 = 1e-11;

/// **A step moves each end by what its half cell carried to the surface that radiates.**
#[test]
fn one_step_moves_what_the_surfaces_send() {
    let dx = 1e-3;
    let area = dx * dx;
    let capacity = RHO * C * dx * dx * dx;
    let poor = Pair {
        name: "",
        k: [0.6, 0.6],
        from: [900.0, 600.0],
        emissivity: 0.9,
        h: 0.0,
        air: 300.0,
    };
    let pairs = [
        Pair {
            name: "across a vacuum",
            ..poor
        },
        Pair {
            name: "across air",
            h: 50.0,
            ..poor
        },
        Pair {
            name: "across air, from an insulator to a conductor",
            k: [0.6, 50.0],
            h: 50.0,
            ..poor
        },
        // Each surface above its own centre, one of them above both. The centres' tangent still
        // bounds what each cell feels here; the limit is the sources' tangent all the same.
        Pair {
            name: "across air hotter than both",
            from: [600.0, 400.0],
            h: 50.0,
            air: 1500.0,
            ..poor
        },
        // The insulator's surface the air holds far above both centres, and its source is the
        // hotter: the limit has to read the far end's. The centres' tangent would charge the
        // conductor's row a third too little and step it at θ = 1.49.
        Pair {
            name: "across air hotter than both, from a conductor to an insulator",
            k: [50.0, 0.05],
            from: [400.0, 300.0],
            h: 200.0,
            air: 3000.0,
            ..poor
        },
        // Nothing reaches this end's surface from its cell, and the air holds it.
        Pair {
            name: "across air, from a cell that conducts nothing",
            k: [0.0, 0.6],
            h: 50.0,
            ..poor
        },
        Pair {
            name: "a black pair at 3000 K on insulators",
            k: [0.05, 0.05],
            from: [3000.0, 300.0],
            emissivity: 1.0,
            ..poor
        },
    ];

    let mut wrong = Vec::new();
    for Pair {
        name,
        k: [ka, kb],
        from: [ta, tb],
        emissivity,
        h,
        air,
    } in pairs
    {
        let mut block = Solid3D::new(
            "pair",
            stuff(ka, emissivity),
            (3, 1, 1),
            Length::from_si(dx),
            Temperature::from_si(ta),
        )
        .empty(|i, _, _| i == 1);
        block.fill(stuff(kb, emissivity), |i, _, _| i == 2);
        block.set_temperature(2, 0, 0, Temperature::from_si(tb));
        let block = if h > 0.0 {
            block.air_in(Temperature::from_si(air), h)
        } else {
            block
        };

        let pair = plates(emissivity) * area;
        let (ga, gb) = (2.0 * ka * dx, 2.0 * kb * dx);
        let face = h * area;
        let [q, u, v] = crossed(pair, (ga, ta), (gb, tb), face, air);

        // The limit: each end's row is the air on its face, behind its half cell, and the pair's
        // tangent at the hotter of the two centres and — with air hotter than both — the two
        // surfaces' sources, which the air holds above the centres.
        let behind = |g: f64| {
            if face > 0.0 && g > 0.0 {
                1.0 / (1.0 / face + 1.0 / g)
            } else {
                0.0
            }
        };
        let source = |g: f64, t: f64| (g * t + face * air) / (g + face);
        let centres = ta.max(tb);
        let hot = if face > 0.0 && air > centres {
            centres.max(source(ga, ta)).max(source(gb, tb))
        } else {
            centres
        };
        let row = |g: f64, hot: f64| capacity / (behind(g) + 4.0 * pair * hot * hot * hot);
        let want = row(ga, hot).min(row(gb, hot));
        let limit = block.max_stable_dt(Time::ZERO).to_si();
        let missed = (limit / want - 1.0).abs();
        // The same terms in another order, a dozen roundings a side: measured, 4.4e-16 at most.
        if off(missed, 1e-14) {
            wrong.push(format!(
                "{name}: the limit is {want:.9e} s, and the block reports {limit:.9e} s"
            ));
        }
        if hot != centres {
            let under = row(ga, centres).min(row(gb, centres));
            if (under / want - 1.0).abs() <= 1e-9 {
                wrong.push(format!(
                    "{name}: the centres' tangent gives {under:e} s, too close to tell apart"
                ));
            }
        }

        // **The theorem's own precondition, row by row.** A step of `dt` leaves each cell a convex
        // combination of the two and the air when `dt·(G_s + f²K)/C ≤ 1`: `G_s` its air behind its
        // half cell, `f` its half cell's share of its surface, and `K` what joins the two sources —
        // the two halves of the path and the secant the surfaces send at, in series. That is what
        // the limit has to bound, read off the surfaces themselves, and not the limit's formula.
        let theta = |ta: f64, tb: f64, dt: f64| -> [f64; 2] {
            let [_, u, v] = crossed(pair, (ga, ta), (gb, tb), face, air);
            let secant = pair * (u + v) * (u * u + v * v);
            let join = 1.0 / (1.0 / (ga + face) + 1.0 / (gb + face) + 1.0 / secant);
            let share = |g: f64| g / (g + face);
            [ga, gb].map(|g| dt * (behind(g) + share(g) * share(g) * join) / capacity)
        };
        // An ulp's worth of rounding in either reading, and nothing more.
        let most = 1.0 + 1e-12;
        for (end, at) in ["near", "far"].iter().zip(theta(ta, tb, limit)) {
            if off(at, most) {
                wrong.push(format!(
                    "{name}: at its own limit of {limit:.6e} s the {end} end's row is θ = {at:.6}"
                ));
            }
        }

        // A quarter of the smaller of the two limits, so a block whose own limit is wrong still
        // takes the step and shows the heat it moved, rather than only refusing it.
        let dt = 0.25 * if limit > 0.0 { want.min(limit) } else { want };
        // What each half cell carries, `g·(t − u)`, with `u` written out: a well-conducting end's
        // surface is within a millikelvin of its centre, and their difference would round.
        let near = ta - dt * ga * (face * (ta - air) + q) / (ga + face) / capacity;
        let far = tb + dt * gb * (face * (air - tb) + q) / (gb + face) / capacity;
        let to_air =
            dt * face * ((ga * (ta - air) - q) / (ga + face) + (gb * (tb - air) + q) / (gb + face));
        println!("  {name}: the surfaces at {u:.6} K and {v:.6} K");

        // **What this case exists to tell apart**, each wrong reading required to move the near
        // end by a hundred thousand times what a cell is held to: the centres' reading, the
        // surfaces solved without the air beside them, and with air, the whole exchange charged
        // to the cells rather than the share their half cells carry.
        let mut others = vec![{
            let q = pair * (fourth(ta) - fourth(tb));
            (
                "the centres' reading",
                ta - dt * (behind(ga) * (ta - air) + q) / capacity,
            )
        }];
        // Neither where nothing reaches the near surface from its cell: without the air that
        // surface has nothing at all, and the cell moves by nothing either way.
        if face > 0.0 && ga > 0.0 {
            let [bare, _, _] = crossed(pair, (ga, ta), (gb, tb), 0.0, air);
            others.push((
                "the surfaces without their air",
                ta - dt * (behind(ga) * (ta - air) + bare) / capacity,
            ));
            others.push((
                "the whole exchange charged to the cells",
                ta - dt * (behind(ga) * (ta - air) + q) / capacity,
            ));
        }
        for (what, other) in others {
            let apart = (other - near).abs();
            if apart.is_nan() || apart <= 1e5 * KELVIN {
                wrong.push(format!(
                    "{name}: {what} leaves the near end at {other} K and the surfaces at {near} K, \
                     too close for this case to tell them apart"
                ));
            }
        }

        let mut stepped = block.clone();
        if let Err(v) = stepped.step(Time::ZERO, Time::from_si(dt), &mut Exchange::new()) {
            wrong.push(format!("{name}: refused a quarter of its own limit: {v:?}"));
            continue;
        }
        let (got_near, got_far) = (
            stepped.temperature_at(0, 0, 0).to_si(),
            stepped.temperature_at(2, 0, 0).to_si(),
        );
        let lost = stepped.lost_energy().to_si();
        let moved = capacity * ((ta - got_near).abs() + (got_far - tb).abs());
        println!(
            "  {name}: {q:.6e} W across, ends {:.1e} K and {:.1e} K off, the air's {lost:.6e} J \
             against {to_air:.6e} J, {:.1e} of what moved, the limit {missed:.1e}",
            (got_near - near).abs(),
            (got_far - far).abs(),
            (lost - to_air).abs() / moved,
        );
        if off((got_near - near).abs(), KELVIN) {
            wrong.push(format!(
                "{name}: the near end should be at {near:.9} K, and is at {got_near:.9} K"
            ));
        }
        if off((got_far - far).abs(), KELVIN) {
            wrong.push(format!(
                "{name}: the far end should be at {far:.9} K, and is at {got_far:.9} K"
            ));
        }
        // Without air nothing leaves the pair, and the books say exactly nothing. With it, what
        // the air took from both surfaces, against the heat each cell gave up or took, which is
        // what the books are made of: measured, 2.0e-16 of it at most.
        if face == 0.0 {
            if lost != 0.0 {
                wrong.push(format!(
                    "{name}: nothing is beside the surfaces, and {lost:e} J left"
                ));
            }
        } else if off((lost - to_air).abs() / moved, 1e-14) {
            wrong.push(format!(
                "{name}: the air takes {to_air:e} J from the two surfaces, and the block lost \
                 {lost:e} J"
            ));
        }

        // At its own limit, read before each step, for forty steps: each cell lands between the
        // coolest and the hottest of where the two began and the air, to within the update's
        // rounding, four ulps of the largest.
        let (low, high) = if face > 0.0 {
            (ta.min(tb).min(air), ta.max(tb).max(air))
        } else {
            (ta.min(tb), ta.max(tb))
        };
        let slack = 4.0 * f64::EPSILON * high;
        let mut marched = block.clone();
        for step in 1..=40 {
            let dt = marched.max_stable_dt(Time::ZERO);
            let now = [
                marched.temperature_at(0, 0, 0).to_si(),
                marched.temperature_at(2, 0, 0).to_si(),
            ];
            let rows = theta(now[0], now[1], dt.to_si());
            if rows.iter().any(|at| off(*at, most)) {
                wrong.push(format!(
                    "{name}: step {step} at its own limit starts from {now:?} K with rows at θ = \
                     {rows:?}"
                ));
                break;
            }
            if let Err(v) = marched.step(Time::ZERO, dt, &mut Exchange::new()) {
                wrong.push(format!(
                    "{name}: step {step} at its own limit was refused: {v:?}"
                ));
                break;
            }
            let ends = [
                marched.temperature_at(0, 0, 0).to_si(),
                marched.temperature_at(2, 0, 0).to_si(),
            ];
            if ends
                .iter()
                .any(|t| t.is_nan() || *t < low - slack || *t > high + slack)
            {
                wrong.push(format!(
                    "{name}: step {step} at its own limit left the ends at {ends:?} K, outside \
                     {low} K to {high} K"
                ));
                break;
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// **Air without limit holds both surfaces at its own temperature, and a block in it settles.**
///
/// Two cells across a cell of void, one generating `P`, in air whose coefficient is infinite: each
/// surface is held at the air, so nothing crosses, and each cell loses only through its half cell,
/// `2k·dx`. The heated cell settles `P/(2k·dx)` above the air and the other at the air. Taking
/// `∞/∞` for the air's share of a surface made both surfaces not a number, the steady solve skip
/// both cells and, across two bodies, report that it diverged.
#[test]
fn a_block_in_air_without_limit_settles_with_nothing_crossing() {
    let (dx, k, watts, room) = (1e-3, 0.6, 0.05, 300.0);
    let block = Solid3D::new(
        "pair",
        stuff(k, 0.9),
        (3, 1, 1),
        Length::from_si(dx),
        Temperature::from_si(500.0),
    )
    .empty(|i, _, _| i == 1)
    .air_in(Temperature::from_si(room), f64::INFINITY)
    .dissipating(watts, |i, _, _| i == 0);
    let t = block
        .steady_state()
        .expect("air without limit still has a balance to strike");
    let heated = room + watts / (2.0 * k * dx);
    println!(
        "  the heated cell at {:.9} K against {heated:.9} K, the other at {:.9} K",
        t[0], t[2]
    );
    assert!(
        (t[0] - heated).abs() < CELL && (t[2] - room).abs() < CELL,
        "the cells settled at {} K and {} K, where the air holds them at {heated} K and {room} K",
        t[0],
        t[2]
    );
}
