//! A face that radiates sheds at its surface's temperature, not at its cell centre's.
//!
//! A finite-volume cell knows its centre, and a film acts at the surface, half a cell away. For a
//! convective film that distance is a resistance in series and nothing else, because `h` is the
//! same at either end of it, and putting the half cell in series made the cooled boundary second
//! order. A radiating film's conductance grows as `T³`, so *where* it is read matters: the block
//! read it at the centre, which the half cell's drop puts above the surface, and the face shed too
//! much by an amount proportional to `dx`. **First order**, and the conservation audit was happy,
//! because the heat really did leave.
//!
//! Measured on the slab below, 16 mm deep, generating heat throughout and radiating from one face
//! at 800 K with `ε = 0.9` and `k = 0.6`: every cell 10.06 K too cold at 1 mm, falling by 1.87,
//! 1.94, 1.97, 1.98 a halving. Behind a mounting it did not converge at all — the mounting's drop
//! does not shrink with the cell, so the surface was charged at least 10.2 K hot at every grid,
//! and the error went to −6.21 K instead of to zero.
//!
//! # What is checked against what
//!
//! **The slab, steady.** Uniform heating `q'''` over a depth `L`, every face insulated but one,
//! whose film sheds all of `q'''L`. The surface is at the root of
//! `h(T_s − T∞) + εσ(T_s⁴ − T∞⁴) = q'''L`, behind a mounting it is that much above the root again,
//! and the inside is the parabola `T_s + q'''(L² − x²)/2k`. A cell-centred scheme whose boundary is
//! exact sits **`q'''dx²/8k` above that parabola at every cell** — its half cell's straight line
//! against the parabola's curve, and the interior's second difference is exact on a quadratic — so
//! that is the closed form, not a limit to approach. A convecting face has always met it; this
//! holds the radiating ones to it, radiating alone, with convection, behind a mounting, at 1500 K,
//! and into a room hotter than the slab, which a sink cools.
//!
//! **One step of one cell.** The heat a step of `dt` loses is `dt` times what the surface sends,
//! `A·(h(T_s − T∞) + εσ(T_s⁴ − T∞⁴))`, at the temperature `T_s` at which what the half cell and
//! the mounting bring it is what it sends, found here by bisection; the cell then sits that many
//! joules lower. And the limit is what it was — the tangent at the hotter of the centre and the
//! room, which bounds the secant at the surface — so forty steps at it stay between the cell's
//! start and the room's, which an explicit step at `θ ≤ 1` does as a theorem rather than a
//! tolerance.

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

/// One face's film: its convection, its emissivity, the room, and the mounting it sits behind,
/// `None` for none.
#[derive(Clone, Copy)]
struct Surface {
    h: f64,
    emissivity: f64,
    ambient: f64,
    mount: Option<f64>,
}

impl Surface {
    /// What a square metre of it sends at `t`, W/m².
    fn sends(&self, t: f64) -> f64 {
        self.h * (t - self.ambient) + self.emissivity * SIGMA * (fourth(t) - fourth(self.ambient))
    }

    /// `block` with this film on `face`, one cell's area, and the mounting if there is one.
    fn on(&self, block: Solid3D, face: Face, dx: f64) -> Solid3D {
        let block = block.losing_from(
            face,
            Environment {
                ambient: Temperature::from_si(self.ambient),
                convection_w_per_m2_k: self.h,
                area: Area::from_si(dx * dx),
            },
        );
        match self.mount {
            Some(m) => block.mounted_on(face, m),
            None => block,
        }
    }
}

struct Slab {
    name: &'static str,
    k: f64,
    surface: Surface,
    /// Where the outer surface settles, K: the heating is what makes it.
    settles_at: f64,
}

/// **A radiating slab meets the cell-centred scheme's own closed form at every cell**, as a
/// convecting one always has.
#[test]
fn a_radiating_slab_sits_where_a_convecting_one_does() {
    const DEPTH: f64 = 0.016;
    let radiating = Surface {
        h: 0.0,
        emissivity: 0.9,
        ambient: 300.0,
        mount: None,
    };
    let slabs = [
        // The control: a convecting face, which the series form made second order long ago.
        Slab {
            name: "convecting only",
            k: 0.6,
            surface: Surface {
                h: 200.0,
                emissivity: 0.0,
                ..radiating
            },
            settles_at: 800.0,
        },
        Slab {
            name: "radiating",
            k: 0.6,
            surface: radiating,
            settles_at: 800.0,
        },
        Slab {
            name: "radiating and convecting",
            k: 0.6,
            surface: Surface {
                h: 50.0,
                ..radiating
            },
            settles_at: 800.0,
        },
        Slab {
            name: "radiating behind a mounting",
            k: 0.6,
            surface: Surface {
                mount: Some(2000.0),
                ..radiating
            },
            settles_at: 800.0,
        },
        Slab {
            name: "radiating at 1500 K",
            k: 0.6,
            surface: radiating,
            settles_at: 1500.0,
        },
        // A sink inside, and a room hotter than the slab: the surface takes heat in, and its own
        // temperature is between the centre's and the room's from the other side.
        Slab {
            name: "radiating in from a hotter room",
            k: 0.6,
            surface: Surface {
                ambient: 600.0,
                ..radiating
            },
            settles_at: 400.0,
        },
    ];

    let mut wrong = Vec::new();
    for Slab {
        name,
        k,
        surface,
        settles_at,
    } in slabs
    {
        // What the face sheds, W/m², and so the heating, W/m³, that sheds it.
        let q = surface.sends(settles_at);
        let heating = q / DEPTH;
        // The face the half cell reaches: the outer surface, and behind a mounting that plus the
        // mounting's drop. The outer surface is at `settles_at` by construction — the heating was
        // chosen as what it sends there — and that is not circular, because the block is told
        // the heating and not the temperature: it starts there, and a solve that moved nothing
        // would miss the parabola by 70 K to 3400 K.
        let inner = settles_at + surface.mount.map_or(0.0, |m| q / m);
        let mut previous = f64::NAN;
        for n in [4usize, 8, 16, 32] {
            let dx = DEPTH / n as f64;
            let block = Solid3D::new(
                "slab",
                stuff(k, surface.emissivity),
                (n, 1, 1),
                Length::from_si(dx),
                Temperature::from_si(settles_at),
            );
            let block = surface
                .on(block, Face::XMax, dx)
                .dissipating(heating * DEPTH * dx * dx, |_, _, _| true);
            let t = match block.steady_state() {
                Ok(t) => t,
                Err(v) => {
                    wrong.push(format!("{name}, {n} cells: no steady state: {v:?}"));
                    continue;
                }
            };
            let offset = heating * dx * dx / (8.0 * k);
            // The steady solve stops when a sweep moves the field by less than a microkelvin, and
            // what that leaves is the tolerance: measured, 9.6e-6 K at most over these cases, and
            // 1.6e-5 K at twice the cells, so 1e-4 K is ten times what is here, against an error
            // of 2.3 K and more when the surface was read at the centre.
            let mut worst: f64 = 0.0;
            for (i, got) in t.iter().enumerate() {
                let x = (i as f64 + 0.5) * dx;
                let parabola = inner + heating / (2.0 * k) * (DEPTH * DEPTH - x * x);
                let miss = got - (parabola + offset);
                if miss.is_nan() || miss.abs() > worst.abs() {
                    worst = miss;
                }
            }
            if off(worst.abs(), 1e-4) {
                wrong.push(format!(
                    "{name}, {n} cells: a cell is {worst:+.3e} K from the parabola plus \
                     q'''dx²/8k = {offset:.3e} K"
                ));
            }
            // The order itself, which the cells above bound only loosely — a cell 1e-4 K out moves
            // the ratio by up to 0.0075 at 800 K and 0.029 in from the hotter room, whose offset is
            // 0.017 K at 32 cells — so this is a check of its own, and the claim the file is named
            // for. Measured, every ratio here is within 2.2e-5 of four, and 0.01 is far
            // above that and far below the 1.9 to 2.6 the centre's reading gave.
            let centre = DEPTH - 0.5 * dx;
            let last = t[n - 1] - (inner + heating / (2.0 * k) * (DEPTH * DEPTH - centre * centre));
            let ratio = previous / last;
            println!("  {name}, {n} cells: error at the face {last:+.6e} K, {ratio:.4} a halving");
            if previous.is_finite() && off((ratio - 4.0).abs(), 0.01) {
                wrong.push(format!(
                    "{name}, {n} cells: the error fell by {ratio:.4} a halving, not 4"
                ));
            }
            previous = last;
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

struct Cell {
    name: &'static str,
    k: f64,
    dx: f64,
    surface: Surface,
    from: f64,
}

/// **No half cell at all, in the limit, and a block that refuses to step it.**
///
/// This was a case of the table below: a cell that conducts without limit, whose surface is its
/// centre, which is the reading the surface's replaced and the right one there. A block refuses
/// that conductivity now, by name and before a step moves anything — it reads as `NaN` once mixed
/// by phase, and between two cells a face of it carried nothing at all. What is left to hold is the
/// limit, asked of the operator rather than stepped: the capacity over the film's tangent alone,
/// with nothing in series. The tolerance is the table's.
#[test]
fn a_cell_that_conducts_without_limit_has_the_films_limit_and_is_refused() {
    let (dx, from) = (1e-3, 900.0);
    let surface = Surface {
        h: 0.0,
        emissivity: 0.9,
        ambient: 300.0,
        mount: None,
    };
    let a = dx * dx;
    let block = surface.on(
        Solid3D::new(
            "cell",
            stuff(f64::INFINITY, surface.emissivity),
            (1, 1, 1),
            Length::from_si(dx),
            Temperature::from_si(from),
        ),
        Face::XMin,
        dx,
    );
    let tangent = a * 4.0 * surface.emissivity * SIGMA * from * from * from;
    let want = RHO * C * dx * dx * dx / tangent;
    let limit = block.max_stable_dt(Time::ZERO).to_si();
    assert!(
        !off((limit / want - 1.0).abs(), 1e-14),
        "the limit is C/G = {want:.6} s, and the block reports {limit:.6} s"
    );
    let mut stepped = block.clone();
    let refusal = stepped
        .step(Time::ZERO, Time::from_si(0.25 * want), &mut Exchange::new())
        .expect_err("a conductivity without limit is refused");
    assert!(
        refusal.quantity.contains("stuff's conductivity") && refusal.before == f64::INFINITY,
        "the refusal names the conductivity and carries it: {refusal}"
    );
    assert_eq!(
        stepped.temperature_at(0, 0, 0).to_si(),
        from,
        "a refused step moved the cell"
    );
}

/// **A step loses what the surface sends**, at the surface's own temperature, and the limit that
/// charged the centre's tangent still bounds it.
#[test]
fn a_step_loses_what_the_surface_sends() {
    let cells = [
        Cell {
            name: "behind a mounting",
            k: 0.6,
            dx: 1e-3,
            surface: Surface {
                h: 0.0,
                emissivity: 0.9,
                ambient: 300.0,
                mount: Some(2000.0),
            },
            from: 900.0,
        },
        Cell {
            name: "radiating and convecting",
            k: 0.6,
            dx: 1e-3,
            surface: Surface {
                h: 50.0,
                emissivity: 0.9,
                ambient: 300.0,
                mount: None,
            },
            from: 900.0,
        },
        // The surface is between the centre and the room from below, so a search that started at
        // the centre would start on the wrong side of it.
        Cell {
            name: "in a hotter room",
            k: 0.6,
            dx: 1e-3,
            surface: Surface {
                h: 0.0,
                emissivity: 0.9,
                ambient: 1200.0,
                mount: None,
            },
            from: 300.0,
        },
        // The longest search there is: nine steps, the surface 1686 K under the centre.
        Cell {
            name: "a black face at 3000 K on an insulator",
            k: 0.05,
            dx: 1e-3,
            surface: Surface {
                h: 0.0,
                emissivity: 1.0,
                ambient: 300.0,
                mount: None,
            },
            from: 3000.0,
        },
        // A mounting that conducts without limit is no mounting.
        Cell {
            name: "behind a mounting that conducts without limit",
            k: 0.6,
            dx: 1e-3,
            surface: Surface {
                h: 0.0,
                emissivity: 0.9,
                ambient: 300.0,
                mount: Some(f64::INFINITY),
            },
            from: 900.0,
        },
        // The shortest: the half cell drops 6.6 mK.
        Cell {
            name: "aluminium a hundred kelvin above the room",
            k: 167.0,
            dx: 2e-3,
            surface: Surface {
                h: 10.0,
                emissivity: 0.1,
                ambient: 300.0,
                mount: None,
            },
            from: 400.0,
        },
    ];

    let mut wrong = Vec::new();
    for Cell {
        name,
        k,
        dx,
        surface,
        from,
    } in cells
    {
        let a = dx * dx;
        let capacity = RHO * C * dx * dx * dx;
        let block = surface.on(
            Solid3D::new(
                "cell",
                stuff(k, surface.emissivity),
                (1, 1, 1),
                Length::from_si(dx),
                Temperature::from_si(from),
            ),
            Face::XMin,
            dx,
        );
        // What lies between the centre and the outer surface, W/K.
        let half = 2.0 * k * a / dx;
        let behind = match surface.mount {
            Some(m) => 1.0 / (1.0 / half + 1.0 / (m * a)),
            None => half,
        };

        // The limit, unchanged: the capacity over the tangent at the hotter of the centre and the
        // room, in series with the same path. Measured, 4.4e-16 at most.
        let hot = from.max(surface.ambient);
        let tangent = a * (surface.h + 4.0 * surface.emissivity * SIGMA * hot * hot * hot);
        let want = capacity * (1.0 / behind + 1.0 / tangent);
        let limit = block.max_stable_dt(Time::ZERO).to_si();
        if off((limit / want - 1.0).abs(), 1e-14) {
            wrong.push(format!(
                "{name}: the limit is C/G = {want:.6} s, and the block reports {limit:.6} s"
            ));
        }

        // The outer surface, where what `behind` brings it is what it sends; bisection, from the
        // two temperatures it must lie between.
        let outer = |behind: f64, surface: &Surface| {
            let (mut lo, mut hi) = (from.min(surface.ambient), from.max(surface.ambient));
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                if behind * (from - mid) - a * surface.sends(mid) > 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            0.5 * (lo + hi)
        };
        let dt = 0.25 * if limit > 0.0 { limit.min(want) } else { want };
        let at = outer(behind, &surface);
        let shed = dt * a * surface.sends(at);

        // **What this case exists to tell apart**, each wrong reading in closed form and required
        // to differ from the right one by a hundred thousand times the tolerance below, so that a
        // case simplified until it could not see its defect says so: the surface read at the
        // centre, which is what this replaced and what a search started on the wrong side of the
        // root ends at; a mounting left out of the path; and convection left out of the surface.
        // Every case has something behind its surface: the one that had nothing, where the
        // centre's reading is the right one, is a test of its own, because a block refuses it.
        let mut others = Vec::new();
        let gap = from - surface.ambient;
        let centre = 1.0 / (1.0 / behind + 1.0 / (a * surface.sends(from) / gap));
        others.push(("the centre's reading", dt * centre * gap));
        if surface.mount.is_some_and(f64::is_finite) {
            let bare = Surface {
                mount: None,
                ..surface
            };
            others.push(("no mounting", dt * a * bare.sends(outer(half, &bare))));
        }
        if surface.h > 0.0 {
            let still = Surface { h: 0.0, ..surface };
            others.push(("no convection", dt * a * still.sends(outer(behind, &still))));
        }
        for (what, other) in others {
            let apart = (other / shed - 1.0).abs();
            if apart.is_nan() || apart <= 1e-9 {
                wrong.push(format!(
                    "{name}: {what} sheds {other:e} J and the surface {shed:e} J, too close for \
                     this case to tell them apart"
                ));
            }
        }

        let mut stepped = block.clone();
        let mut bus = Exchange::new();
        if let Err(v) = stepped.step(Time::ZERO, Time::from_si(dt), &mut bus) {
            wrong.push(format!("{name}: refused a quarter of its own limit: {v:?}"));
            continue;
        }
        // Newton in the block and bisection here meet the root within an ulp or so, and what is
        // sent moves as `4δT/T` with it; measured, 8.9e-16 at most, and 1e-14 is eleven times
        // that against the 0.86% the centre's reading cost a cell of liquid at 900 K.
        let lost = stepped.lost_energy().to_si();
        if off((lost / shed - 1.0).abs(), 1e-14) {
            wrong.push(format!(
                "{name}: a step of {dt:.4} s loses what the surface sends at {at:.6} K, \
                 {shed:e} J, and the block lost {lost:e} J"
            ));
        }
        let t = stepped.temperature_at(0, 0, 0).to_si();
        let expected = from - shed / capacity;
        // An ulp at 3000 K is 4.5e-13 K, and that is the most measured.
        if off((t - expected).abs(), 1e-11) {
            wrong.push(format!(
                "{name}: one step on it should be at {expected:.9} K, and is at {t:.9} K"
            ));
        }

        // Steps at the block's own limit, read before each, until the cell is at the room. A step
        // at `θ ≤ 1` is a convex combination of where it began and the room, so each lands between
        // the two — `0 ≤ (T′ − T∞)/(T − T∞) ≤ 1` — to within the update's rounding, four ulps of
        // the larger; and near the room, where the secant meets the tangent, `θ` reaches one and
        // the cell lands on the room. Measured, in four to ten steps.
        let mut marched = block.clone();
        let mut was = from;
        // Arrived at the room, failed on the way, or neither after forty steps.
        let mut ended = None;
        for step in 1..=40 {
            let dt = marched.max_stable_dt(Time::ZERO);
            if let Err(v) = marched.step(Time::ZERO, dt, &mut Exchange::new()) {
                wrong.push(format!(
                    "{name}: step {step} at its own limit was refused: {v:?}"
                ));
                ended = Some(Err(()));
                break;
            }
            let t = marched.temperature_at(0, 0, 0).to_si();
            let slack = 4.0 * f64::EPSILON * was.abs().max(surface.ambient.abs());
            let (low, high) = (was.min(surface.ambient), was.max(surface.ambient));
            if t.is_nan() || t < low - slack || t > high + slack {
                wrong.push(format!(
                    "{name}: step {step} at its own limit went from {was} K to {t} K, outside \
                     {low} K to {high} K"
                ));
                ended = Some(Err(()));
                break;
            }
            was = t;
            if (t - surface.ambient).abs() <= slack {
                ended = Some(Ok(step));
                break;
            }
        }
        match ended {
            Some(Ok(step)) => println!("  {name}: at the room after {step} steps at its own limit"),
            Some(Err(())) => {}
            None => wrong.push(format!(
                "{name}: forty steps at its own limit left it at {was} K, short of the room"
            )),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
