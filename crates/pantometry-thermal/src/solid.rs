//! Conduction through a block, in three dimensions.
//!
//! [`Bar1D`](crate::Bar1D) resolves a gradient along one axis, and that is the right model for a
//! rod, a wall or a beam landing on a strip. A heat sink is not any of those. Heat spreading
//! sideways out of a hot spot is exactly what a fin, a spreader plate and a mounting boss are
//! *for*, and a one-dimensional model cannot show it: it has nowhere for the heat to go but
//! along.
//!
//! # What the third dimension costs
//!
//! Twice, and the second time is the one that hurts.
//!
//! A block of `n` cells on a side is `n³` cells rather than `n`. And the explicit stability
//! limit tightens with each axis, because the limit is on the *sum* of what the three
//! directions do in one step:
//!
//! ```text
//!   1D      α·dt/dx² ≤ 1/2
//!   2D      α·dt/dx² ≤ 1/4
//!   3D      α·dt/dx² ≤ 1/6
//! ```
//!
//! So a 3D block at the same spacing takes three times as many steps as a bar, each of them n²
//! times more work. `Room` in `pantometry-acoustic` records the same trade for the wave equation and
//! reaches a factor of √3 rather than 3, because a wave's limit is on the wave speed and a
//! diffusion limit is on the diffusivity — one is linear in the sum and the other in its square
//! root.
//!
//! That is why [`LumpedMass`](crate::LumpedMass) and [`ThermalNetwork`](crate::ThermalNetwork)
//! are not going anywhere. In aluminium at a millimetre the step is **2.41 ms**, so a motor
//! housing over its two-thousand-second thermal time constant is **828,000 steps** — each of
//! them a sweep over however many cells the housing is, which at that spacing is around a
//! million. A graph of four nodes answers the same question immediately. Use the cheapest model
//! whose reduction still holds, and `LumpedMass::biot_number` is how you find out whether it
//! does.
//!
//! What a block is *for* is the case where the reduction does not hold: a hot spot, a spreader,
//! a gradient across a joint. Those are questions a lumped model cannot answer at any price.
//!
//! # What it is checked against
//!
//! A separable cosine mode is an **exact eigenvector of the discrete operator**, not merely an
//! approximate solution of the continuum one. On a cell-centred grid with insulated faces the
//! mode `cos(aπ(i+½)/nx)·cos(bπ(j+½)/ny)·cos(cπ(k+½)/nz)` decays by exactly the same factor
//! every step, and that factor is known in closed form. So the test is an equality at machine
//! precision rather than a tolerance, and it is sensitive to a swapped axis, a dropped term or a
//! wrong spacing in a way a smooth decaying blob would not be.
//!
//! The *continuum* rate is then a second test: the discrete eigenvalue approaches
//! `−α·π²·(1/Lx² + 1/Ly² + 1/Lz²)` at second order, and refining the grid quarters the error.
//! Rate rather than value, because a scheme that is first order where it claims to be second is
//! the defect this workspace has already shipped once.
//!
//! # A workaround that is gone
//!
//! This domain briefly offered `Domain::as_bodies` alongside `as_field`, so a viewer could get
//! its cells as a point cloud. That was not a design choice; it was cover for `pantometry-scene`'s
//! `Extent` being two-dimensional, which would have captured a block as its `z = 0` face.
//!
//! `Extent` is three-dimensional now, so the cover is unnecessary — and it was never free. A
//! domain that is two shapes at once makes the picture depend on whether somebody remembered to
//! set an extent, which is a mode nothing announces. It is a field, and only a field.

use glam::DVec3;
use pantometry_core::conserved::quantity;
use pantometry_core::Lattice;
use pantometry_core::Reading;
use pantometry_core::{
    units::{
        Area, Conductance, Energy, HeatCapacity, Length, LengthVec, Power, Temperature, Time,
        Volume,
    },
    Domain, Exchange, Ledger, ScalarField, Substance, Violation,
};

use crate::{Environment, HEAT};
use pantometry_units::STEFAN_BOLTZMANN;
use std::collections::BTreeMap;

/// The largest Fourier number an explicit three-dimensional sweep is stable at.
///
/// `1/(2d)` for `d` dimensions, from requiring the amplification factor of the worst-resolved
/// mode to stay inside the unit circle. Public because a caller sizing a grid needs it before
/// there is anything to ask.
pub const STABLE_FOURIER_3D: f64 = 1.0 / 6.0;

/// A rectangular block, conducting in three dimensions, of one material or of several.
///
/// Cells are **cubes** of a single spacing rather than boxes of three. That is a deliberate
/// restriction and the same one `Room::of_air` makes: an anisotropic cell makes the stability
/// limit anisotropic and the truncation error different along each axis, so the grid would be
/// resolving one direction better than another for a reason that had nothing to do with the
/// physics. A block that is longer than it is thick is more cells along, not longer cells.
///
/// Faces are **insulated**. Every boundary cell exchanges with the neighbours it has and no
/// others, which is what makes the total exactly conserved rather than conserved to a tolerance.
/// A block that should lose heat gets that from a domain on the other side of the bus.
///
/// [`fill`](Solid3D::fill) gives cells a different substance from the constructor's, which is what
/// makes a coating, a joint or a layered wall expressible. Read that method for the one number the
/// scheme turns on — the conductivity **on a face** is the harmonic mean of the two cells', not
/// A surface conductance put in series with the half cell of solid behind it.
///
/// `1/(1/G_film + dx/(2kA))`. The film is charged against the **surface**, and a finite-volume
/// cell knows only its centre, so the half cell between them is part of the path. Leaving it
/// out sheds too much heat and — measured — turns a second-order boundary into a first-order
/// one, which is the defect this workspace has twice found in an acoustic wall.
///
/// A conductivity that is not finite leaves the film alone. An infinite one is the right limit
/// rather than a special case, and a `NaN` one is what a substance with no thermal properties
/// reads as. A block steps neither — it refuses both by name, the first because between two
/// cells it carried nothing — so what reads this with one is the limit.
fn series_with_half_cell(film: f64, conductivity: f64, area: f64, dx: f64) -> f64 {
    if film <= 0.0 {
        return 0.0;
    }
    let half = 2.0 * conductivity * area / dx;
    if !half.is_finite() {
        return film;
    }
    if half <= 0.0 {
        return 0.0;
    }
    1.0 / (1.0 / film + 1.0 / half)
}

/// How many Newton steps `radiating_film` may take to find its surface. A bound that is never
/// the answer: the steps fall monotonically onto the root, so they stop on their own. Measured,
/// counting the step that finds nothing left to move: three for aluminium or copper near room
/// temperature, four or five for a cell of `k = 0.6` at 900 K to 1073 K, and nine and ten for a
/// black face at 3000 K on an insulator of `k = 0.05`, a millimetre and ten millimetres a cell.
/// Below zero they rise monotonically onto the root from the centre, in `film_surface_below_zero`:
/// over 10⁵ faces drawn as the unit test draws them, fourteen at most and 4.70 on average, and
/// that test holds its eight thousand to sixteen and 4.85. **Far below zero they grow with the
/// centre's distance**, radiation outweighing what is behind the surface, so each step covers a
/// share of the way: measured, 17 to 24 at −1e6 K and 53 to 59 at −1e12 K, and this bound is
/// reached, and is the answer, from −1e13 K to −1e14 K — where a solve has long been drawn down
/// without end, and is refused.
const SURFACE_STEPS: usize = 64;

/// The conductance, W/K, from a cell's centre to the room through a face that radiates.
///
/// `behind` is everything between the centre and the surface — the half cell and any mounting, in
/// series, W/K. The surface sheds `A·(h + εσ(T + T∞)(T² + T∞²))·(T − T∞)` at **its own**
/// temperature `T`, which is the one at which that equals what `behind` brings it:
///
/// ```text
/// behind · (T_c − T)  =  A·(h + εσ(T + T∞)(T² + T∞²)) · (T − T∞)
/// ```
///
/// and the answer is the series of `behind` and that secant, so that a flux computed from it is
/// exactly what reached the surface.
///
/// **It charged the surface's loss at the centre's temperature**, which is the same thing only
/// when the film is linear. A convective film's conductance does not depend on where it is read,
/// and that is why the series form made the convective boundary second order. A radiating one's
/// grows as `T³`, and the centre is hotter than the surface by what crosses the half cell, so
/// the face shed too much by an amount proportional to `dx`: **first order**. Measured on a slab
/// 16 mm deep generating `1.28 MW/m³` and radiating at 800 K from `ε = 0.9`, `k = 0.6`, the cells
/// sat 10.06 K too cold at 1 mm, falling by 1.87, 1.94, 1.97, 1.98 a halving. Behind a mounting
/// it did not converge at all: the mounting's drop does not shrink with the cell, so the surface
/// was charged at a temperature above its own by `q(1/h_mount + dx/2k)` — 27.3 K at 1 mm, 14.5 K
/// at 0.25 mm, and never less than the joint's 10.2 K however fine the grid — and the error went
/// to −6.21 K rather than to zero — the limit of the old balance as the cell vanishes, solved for,
/// with the grid at 0.25 mm still 8.77 K out. Now the error in each is `q'''dx²/8k`, which is what
/// a convective face has — the half cell's straight line against the parabola — and it falls by
/// 4.000 a halving.
///
/// **Newton from `max(T_c, T∞)`, which needs no safeguard.** `behind·(T − T_c)` plus the surface's
/// loss is increasing and convex in `T` — the loss's second derivative is `12εσT²A` — and its
/// root lies between `T_c` and `T∞`. Newton on such a function from a start at or beyond its root
/// never crosses it, so each step lands between the root and the last, and the first that does
/// not move down is where rounding has taken over. Only `+ − × ÷`, so the result is the same on
/// every platform.
///
/// **It cannot loosen the limit.** The secant at `T ≤ max(T_c, T∞)` is at most the tangent
/// `h + 4εσ·max(T_c, T∞)³` that `loss_conductance_at` charges, and putting the same `behind` in
/// series with each keeps the order, so the limit bounds this as it bounded the centre's secant.
/// Above the room that was the larger of the two; below it the surface's is, the surface being
/// nearer the room — a cell at 300 K in a room at 1200 K carries 1.09 times what it did — and the
/// bound holds on both sides because it is read at the hotter of the centre and the room.
///
/// An infinite `behind` puts the surface at the centre, which is the old form and the right one
/// there; a zero one carries nothing.
///
/// **A film without limit holds the surface at the room**, so what carries the heat is `behind`
/// alone, whatever the surface radiates. Solved for, it was `∞/∞`: the first Newton step was not a
/// number and so was the film, and a radiating face cooled without limit stepped its cell to `NaN`
/// while the limit — which reads the film through `series_with_half_cell`, and gets the half cell —
/// said the step was fine. Measured on four cells of `k = 50` radiating at `ε = 0.5`: `NaN` within
/// twenty steps at half the limit, and the steady solve reporting that it diverged.
///
/// **Below absolute zero a surface radiates `T·|T|³`**, so that its loss keeps rising with it, as
/// [`crossing`] reads a pair's. No body is there, but the steady solve passes through on the way
/// to an answer and stops there when there is none above it, and what it reads there decides
/// which: see [`Solid3D::steady_state`]. The balance is concave below zero and convex above, so
/// **where it is still positive at zero the surface is below zero**, and
/// [`film_surface_below_zero`] runs Newton up onto it from the centre, as the loop here runs down
/// onto a root above zero from the room. Where it is not, the surface is at or above zero, even
/// with the centre below it, to the rounding of the centre, and is found here exactly as it was;
/// a centre at or above zero has its surface there always.
#[inline(never)]
fn radiating_film(
    behind: f64,
    share: f64,
    convection: f64,
    radiative: f64,
    centre: f64,
    ambient: f64,
) -> f64 {
    if behind <= 0.0 {
        return 0.0;
    }
    if convection == f64::INFINITY {
        return behind;
    }
    if !behind.is_finite() {
        return film_secant(share, convection, radiative, centre, ambient);
    }
    let (surface, below) = film_surface(behind, share, convection, radiative, centre, ambient);
    let through = if below {
        film_secant(share, convection, radiative, surface, ambient)
    } else {
        film_secant_above_zero(share, convection, radiative, surface, ambient)
    };
    through * behind / (through + behind)
}

/// Where a radiating face's surface sits, for a finite `behind` above nothing and a finite film,
/// and whether it was found below absolute zero — see [`radiating_film`], which this is the
/// first half of, and [`radiating_film_slope`], which reads the same surface.
#[inline(always)]
fn film_surface(
    behind: f64,
    share: f64,
    convection: f64,
    radiative: f64,
    centre: f64,
    ambient: f64,
) -> (f64, bool) {
    // At or above zero, which is where every surface the loop below visits is.
    let secant = |t: f64| film_secant_above_zero(share, convection, radiative, t, ambient);
    // The balance at a surface at zero: what `behind` brings it from the centre, less what it
    // sheds there to the room.
    if centre < 0.0 && -behind * centre > secant(0.0) * ambient {
        let (surface, _) =
            film_surface_below_zero(behind, share, convection, radiative, centre, ambient);
        return (surface, true);
    }
    let mut surface = centre.max(ambient);
    for _ in 0..SURFACE_STEPS {
        let excess = secant(surface) * (surface - ambient) - behind * (centre - surface);
        let slope = share * (convection + 4.0 * radiative * surface * surface * surface) + behind;
        let next = surface - excess / slope;
        // `NaN` stops it too, which a bare `next < surface` written as its negation would say.
        if next.is_nan() || next >= surface {
            break;
        }
        surface = next;
    }
    (surface, false)
}

/// How fast what a radiating face carries, `radiating_film(…)·(T_c − T∞)`, rises with its cell's
/// centre, W/K: `behind` in series with the tangent of what the surface sheds, **at the surface**,
/// `A·(h + 4εσ·s²|s|)` — the derivative of `T·|T|³` on either side of zero.
///
/// From the surface's balance `behind·(T_c − s) = L(s)`: a rise of the centre moves the surface by
/// `behind/(behind + L′)` of it, and what crosses by `L′` times that. **The tangent at the hotter of
/// the centre and the room**, which `face_loss` charges and the stability limit needs as a bound,
/// is not this: for a cell below its room it reads `T∞³` where the surface is at `|s|³`, and
/// [`Solid3D::steady_state`] took it for the slope. A centre with nothing behind its surface has
/// the surface's own tangent, and a film without limit carries `behind` itself.
fn radiating_film_slope(
    behind: f64,
    share: f64,
    convection: f64,
    radiative: f64,
    centre: f64,
    ambient: f64,
) -> f64 {
    let tangent = |s: f64| share * (convection + 4.0 * radiative * s * s * s.abs());
    if behind <= 0.0 {
        return 0.0;
    }
    if convection == f64::INFINITY {
        return behind;
    }
    if !behind.is_finite() {
        return tangent(centre);
    }
    let (surface, _) = film_surface(behind, share, convection, radiative, centre, ambient);
    let at = tangent(surface);
    behind * at / (behind + at)
}

/// What a radiating face sheds per kelvin of its surface `t` above the room, W/K: the secant
/// [`film_secant_above_zero`] at or above absolute zero, and below it
/// `A·(h + εσ(t⁴ + T∞⁴)/(T∞ − t))`, the secant of `T·|T|³`, which meets it at zero. A surface
/// that is not a number takes the first.
fn film_secant(share: f64, convection: f64, radiative: f64, t: f64, ambient: f64) -> f64 {
    if t < 0.0 {
        let (t2, a2) = (t * t, ambient * ambient);
        share * (convection + radiative * (t2 * t2 + a2 * a2) / (ambient - t))
    } else {
        film_secant_above_zero(share, convection, radiative, t, ambient)
    }
}

/// The secant of a radiating face at a surface `t` at or above absolute zero, W/K:
/// `A·(h + εσ(t + T∞)(t² + T∞²))`.
#[inline(always)]
fn film_secant_above_zero(
    share: f64,
    convection: f64,
    radiative: f64,
    t: f64,
    ambient: f64,
) -> f64 {
    share * (convection + radiative * (t + ambient) * (t * t + ambient * ambient))
}

/// The surface of a radiating face below absolute zero, and how many steps it took: where what
/// `behind` brings it from a centre below zero, `behind·(T_c − T)`, is what it sheds,
/// `A·h·(T − T∞) + A·εσ(T·|T|³ − T∞⁴)`. For a face whose balance is still positive at zero, so
/// that the surface lies between the centre and zero.
///
/// **Newton up from the centre, which needs no safeguard**, as the one in `radiating_film` needs
/// none coming down. `behind·(T − T_c)` plus what the surface sheds rises with `T` and is concave
/// below zero, where what it sheds curves as `12εσA·T|T|`, so Newton from below its root never
/// crosses it, each step landing between the last and the root, and the first that does not move
/// up is where rounding has taken over. **A bracket from the room, halved where a step would
/// leave it, was tried first** and took fifty-four steps on a face whose root was an ulp above a
/// centre a few millikelvin below zero, in a room at zero: every step landed on the bracket's end,
/// and was halved. Only `+ − × ÷`.
fn film_surface_below_zero(
    behind: f64,
    share: f64,
    convection: f64,
    radiative: f64,
    centre: f64,
    ambient: f64,
) -> (f64, usize) {
    let mut surface = centre;
    let mut steps = 0;
    for _ in 0..SURFACE_STEPS {
        steps += 1;
        let sheds = film_secant(share, convection, radiative, surface, ambient);
        let excess = sheds * (surface - ambient) - behind * (centre - surface);
        let slope =
            share * (convection + 4.0 * radiative * surface * surface * surface.abs()) + behind;
        let next = surface - excess / slope;
        // `NaN` stops it too, which a bare `next > surface` written as its negation would say.
        if next.is_nan() || next <= surface {
            break;
        }
        surface = next;
    }
    (surface, steps)
}

/// How many steps `crossing` may take to find what a clearance carries. A bound that is never the
/// answer: over the eight thousand pairs the unit tests below draw — conductivities from 0.01 to
/// 500 W/m·K, cells from 0.1 to 30 mm, 200 K to 3000 K, with and without air — eleven at most and
/// 4.35 on average, and that test holds it to twelve. Below zero, in `crossing_below_zero`, over
/// 86 196 of 10⁵ pairs drawn there whose colder source is below zero: fourteen at most and 5.04
/// on average, and that test holds its eight thousand to sixteen and 5.1.
const CROSSING_STEPS: usize = 64;

/// What one pair of a clearance carries, from its end `a` to its end `b`, in a step's arithmetic.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Crossing {
    /// Watts from `a`'s surface to `b`'s; negative the other way.
    watts: f64,
    /// The share of `watts` that `a`'s cell gives up, and the share `b`'s takes: each the share of
    /// its surface's paths that runs through its half cell. One for both without air; with air,
    /// the air takes `watts·(from_a − to_b)` between them.
    from_a: f64,
    to_b: f64,
    /// The two surfaces, K, and the resistance from each one's source to it, K/W: what the slopes
    /// are made of, kept rather than the slopes because only the steady solve asks for them.
    surfaces: (f64, f64),
    behind: (f64, f64),
}

impl Crossing {
    /// How fast what `a`'s cell gives up through the pair moves with `a`'s own temperature, and
    /// what `b`'s takes with `b`'s, W/K — the slopes the steady solve leans on. From
    /// `dq·(1 + 4C(|u|³ r_a + |v|³ r_b)) = 4C(|u|³ dT*_a − |v|³ dT*_b)` and `dT* = f·dT_c`: the
    /// derivative of `T·|T|³`, which is `T⁴`'s at or above zero — see [`crossing`].
    fn slopes(&self, coefficient: f64) -> (f64, f64) {
        let ((u, v), (ra, rb)) = (self.surfaces, self.behind);
        let (ua, vb) = (
            4.0 * coefficient * u * u * u.abs(),
            4.0 * coefficient * v * v * v.abs(),
        );
        let stiff = 1.0 + ua * ra + vb * rb;
        (
            self.from_a * self.from_a * ua / stiff,
            self.to_b * self.to_b * vb / stiff,
        )
    }

    /// The other half of the same derivative: how fast what `a`'s cell **takes in** through the
    /// pair rises with `b`'s temperature, and what `b`'s takes in with `a`'s, W/K. From the same
    /// `dq`: what `a` gives up is `q·f_a`, so its rise with `T_b` is `f_a·f_b·4C|v|³/stiff`, and
    /// `b`'s with `T_a` is `f_a·f_b·4C|u|³/stiff`. Without air each is the other end's slope, so a
    /// shift of both ends together moves nothing across, as it must.
    fn couplings(&self, coefficient: f64) -> (f64, f64) {
        let ((u, v), (ra, rb)) = (self.surfaces, self.behind);
        let (ua, vb) = (
            4.0 * coefficient * u * u * u.abs(),
            4.0 * coefficient * v * v * v.abs(),
        );
        let stiff = 1.0 + ua * ra + vb * rb;
        let both = self.from_a * self.to_b;
        (both * vb / stiff, both * ua / stiff)
    }
}

/// What a clearance's pair exchanges, given the coefficient `σA/(1/ε_a + 1/ε_b − 1)` of its two
/// faces, each end's half cell, W/K, and its centre's temperature, and the air in the clearance
/// on each face, W/K, and its temperature, if it has any.
///
/// **The two surfaces radiate, not the two centres.** Each surface is half a cell from the
/// centre the block knows. Behind it is the half cell `G`, and beside it the air `H`; between
/// them they hold it as a source `T* = T_c − w·(T_c − T_air)`, with `w = H/(G + H)`, behind a
/// conductance `G + H`. So the surfaces are at `u = T*_a − q/(G_a + H)` and
/// `v = T*_b + q/(G_b + H)`, and `q = C(u⁴ − v⁴)`: one equation in `q`, decreasing, whose root
/// lies between nothing and the `q` at which the two surfaces meet. Of what reaches `a`'s surface,
/// the share `G_a/(G_a + H)` came through its half cell; the air's own share of the face is
/// `air_conductance`'s, which carries the rest.
///
/// **Below absolute zero a surface radiates `T·|T|³`**, and the equation keeps falling there.
/// No body is below zero, but the steady solve passes through on the way to its answer and stops
/// there when none is above, and it read `u⁴ − v⁴` as the two fourth powers, even in `T`: a
/// surface drawn through zero then took heat as though as warm as its mirror, and a sink fed
/// across a clearance settled at the mirror of a balance that exists — a milliwatt answered at
/// 878.56 K from five starts was refused at −898.56 K from the sixth, at 300 K. Read odd, `q` falls with
/// either surface's distance below the other everywhere, and the slopes below are its own. A pair
/// whose colder source is at or above zero has both surfaces there and reads exactly what it read;
/// one below is [`crossing_below_zero`].
///
/// **Newton inside that bracket, halving it where a step would leave it.** The residual is not
/// convex — its curvature is `12C(u²/G'_a² − v²/G'_b²)`, of either sign — so Newton alone is not
/// monotone here as it is in `radiating_film`, and the bracket is what makes it safe. It starts
/// at the secant between the two sources' answer. Written as `C·(u − v)(u + v)(u² + v²) − q` with
/// `u − v` formed directly, and the sources' difference from whichever end cancels least, so that
/// two surfaces close together are not the difference of two fourth powers. The unit tests below
/// hold it to its root within what a halving on `a`'s surface can resolve, and its slopes to their
/// derivatives. What is left against the exact root of the exact inputs is how those inputs round —
/// the air's shares, chiefly, where it holds the two sources close together with the centres far
/// apart. Only `+ − × ÷`.
fn crossing(coefficient: f64, a: (f64, f64), b: (f64, f64), air: Option<(f64, f64)>) -> Crossing {
    crossing_counted(coefficient, a, b, air).0
}

/// [`crossing`], and how many steps it took — which the unit tests below hold to its bound.
fn crossing_counted(
    coefficient: f64,
    (half_a, centre_a): (f64, f64),
    (half_b, centre_b): (f64, f64),
    air: Option<(f64, f64)>,
) -> (Crossing, usize) {
    let (h, ambient) = air.unwrap_or((0.0, 0.0));
    // A surface nothing reaches carries nothing and moves with nothing: no slope either.
    let nothing = Crossing {
        watts: 0.0,
        from_a: 1.0,
        to_b: 1.0,
        surfaces: (0.0, 0.0),
        behind: (0.0, 0.0),
    };
    // Each end: the resistance from its source to its surface, the source, the air's share of
    // it and the half cell's — or `None` for a surface nothing reaches, which carries nothing.
    let side = |half: f64, centre: f64| -> Option<(f64, f64, f64, f64)> {
        if h == 0.0 {
            if half <= 0.0 {
                return None;
            }
            let r = if half.is_finite() { 1.0 / half } else { 0.0 };
            return Some((r, centre, 0.0, 1.0));
        }
        // Air without limit holds the surface at its own temperature, whatever is behind it: the
        // half cell carries none of what crosses. Written out because `∞/∞` is not a share.
        if h == f64::INFINITY {
            return Some((0.0, ambient, 1.0, 0.0));
        }
        if !half.is_finite() {
            return Some((0.0, centre, 0.0, 1.0));
        }
        let through = half + h;
        let (w, f) = (h / through, half / through);
        Some((1.0 / through, centre - w * (centre - ambient), w, f))
    };
    let (Some((ra, sa, wa, fa)), Some((rb, sb, wb, fb))) =
        (side(half_a, centre_a), side(half_b, centre_b))
    else {
        return (nothing, 0);
    };
    // The sources' difference, from whichever form subtracts the smaller terms, since what each
    // rounds to is in proportion to them: the centres' difference less the air's shares of their
    // excesses, or the half cells' shares of those excesses. The first where the air holds little,
    // the second where it holds much — and where one end's air holds much and the other's little,
    // whichever the numbers say.
    let (ea, eb) = (centre_a - ambient, centre_b - ambient);
    let by_centres = (centre_a - centre_b).abs() + (wa * ea).abs() + (wb * eb).abs();
    let by_excess = (fa * ea).abs() + (fb * eb).abs();
    let apart = if by_centres <= by_excess {
        (centre_a - centre_b) - (wa * ea - wb * eb)
    } else {
        fa * ea - fb * eb
    };
    if apart == 0.0 {
        let still = Crossing {
            watts: 0.0,
            from_a: fa,
            to_b: fb,
            surfaces: (sa, sb),
            behind: (ra, rb),
        };
        return (still, 0);
    }
    // Solved from the hotter source, and the sign put back after.
    let (hot_r, hot, cold_r, cold, d0) = if apart > 0.0 {
        (ra, sa, rb, sb, apart)
    } else {
        (rb, sb, ra, sa, -apart)
    };
    let r = hot_r + cold_r;
    let mut steps = 0;
    let mut q = if cold < 0.0 {
        let (q, taken) = crossing_below_zero(coefficient, (hot_r, hot), (cold_r, cold), d0);
        steps = taken;
        q
    } else if r == 0.0 {
        // Nothing between either centre and its surface: the centres' own reading, exactly.
        coefficient * d0 * (hot + cold) * (hot * hot + cold * cold)
    } else {
        let (mut lo, mut hi) = (0.0, d0 / r);
        // From the secant between the two sources, which is the answer where nothing is behind
        // the surfaces and is off by the half cells' share of the path otherwise. It took fewer
        // steps on random pairs than the tangent at the hotter source, and cost a slit through a
        // block of 16³ 18.5% a step where the tangent cost 30%.
        let secant = coefficient * (hot + cold) * (hot * hot + cold * cold);
        let mut q = d0 / (r + 1.0 / secant);
        for _ in 0..CROSSING_STEPS {
            steps += 1;
            let (u, v) = (hot - q * hot_r, cold + q * cold_r);
            let gap = d0 - q * r;
            let spread = (u + v) * (u * u + v * v);
            let excess = coefficient * gap * spread - q;
            if excess > 0.0 {
                lo = q;
            } else if excess < 0.0 {
                hi = q;
            } else {
                break;
            }
            let bend =
                (cold_r - hot_r) * (u * u + v * v) + 2.0 * (u + v) * (v * cold_r - u * hot_r);
            let slope = coefficient * (gap * bend - r * spread) - 1.0;
            let mut next = q - excess / slope;
            if next == q {
                break;
            }
            // Written out rather than negated, so that a `NaN` step is one that leaves.
            if next.is_nan() || next <= lo || next >= hi {
                next = 0.5 * (lo + hi);
                if next <= lo || next >= hi {
                    break;
                }
            }
            q = next;
        }
        q
    };
    if apart < 0.0 {
        q = -q;
    }
    let crossed = Crossing {
        watts: q,
        from_a: fa,
        to_b: fb,
        surfaces: (sa - q * ra, sb + q * rb),
        behind: (ra, rb),
    };
    (crossed, steps)
}

/// What a pair carries from its hotter source `hot` to its colder `cold`, when `cold` is below
/// absolute zero, and how many steps it took: `crossing_counted`'s equation with each fourth power
/// read as `T·|T|³`. Each end is `(resistance from its source to its surface, source)`, and `d0`
/// the sources' difference, formed by the caller.
///
/// `C(g(u) − g(v))` is `C·(u − v)(u + v)(u² + v²)` while both surfaces are at or above zero, its
/// negative while both are below, and `C(u⁴ + v⁴)` while one is on each side, which has nothing
/// to cancel. Its derivative in `q` is `−4C(|u|³ r_hot + |v|³ r_cold) − 1` everywhere, so the
/// equation falls the whole way from nothing, where it is positive, to the `q` at which the
/// surfaces meet, where it is negative, and Newton inside that bracket, halving where a step would
/// leave it, finds its one root — as above zero, from the same secant. Only `+ − × ÷`.
fn crossing_below_zero(
    coefficient: f64,
    (hot_r, hot): (f64, f64),
    (cold_r, cold): (f64, f64),
    d0: f64,
) -> (f64, usize) {
    // `g(u) − g(v)` for `u` at or above `v`, given `u − v`.
    let apart = |u: f64, v: f64, gap: f64| {
        if v >= 0.0 {
            gap * (u + v) * (u * u + v * v)
        } else if u <= 0.0 {
            -gap * (u + v) * (u * u + v * v)
        } else {
            (u * u) * (u * u) + (v * v) * (v * v)
        }
    };
    let r = hot_r + cold_r;
    let reading = coefficient * apart(hot, cold, d0);
    if r == 0.0 {
        return (reading, 0);
    }
    let (mut lo, mut hi) = (0.0, d0 / r);
    // The secant between the two sources, `reading/d0`, in series with what is behind them.
    let mut q = reading / (1.0 + r * reading / d0);
    let mut steps = 0;
    for _ in 0..CROSSING_STEPS {
        steps += 1;
        let (u, v) = (hot - q * hot_r, cold + q * cold_r);
        let gap = d0 - q * r;
        let excess = coefficient * apart(u, v, gap) - q;
        if excess > 0.0 {
            lo = q;
        } else if excess < 0.0 {
            hi = q;
        } else {
            break;
        }
        let slope = -4.0 * coefficient * (u * u * u.abs() * hot_r + v * v * v.abs() * cold_r) - 1.0;
        let mut next = q - excess / slope;
        if next == q {
            break;
        }
        // Written out rather than negated, so that a `NaN` step is one that leaves.
        if next.is_nan() || next <= lo || next >= hi {
            next = 0.5 * (lo + hi);
            if next <= lo || next >= hi {
                break;
            }
        }
        q = next;
    }
    (q, steps)
}

/// How many Gauss-Seidel sweeps [`Solid3D::steady_state`] will take before reporting that it did
/// not converge.
///
/// A **bound, not a budget**: the linear cases here converge in tens and the radiative ones in
/// hundreds, and a block that has not settled in two thousand sweeps has something wrong with it
/// that a longer wait will not fix. Bounded so a pathological input reports rather than spins,
/// which is the same reason `ThermalNetwork::steady_state` bounds its Newton at eight.
const STEADY_SWEEPS: usize = 2000;

/// How many bodies one part of a block may hold for [`Solid3D::steady_state`] to solve it: the
/// size of the dense system the bodies a clearance joins take their step in together. Its cost
/// is the cube of this, `5.6e6` operations a round at the bound, against a scene's two bodies.
const STEADY_BODIES: usize = 256;

/// How far below absolute zero [`Solid3D::steady_state`] lets a cell stop and still call it zero,
/// K: ten times the 1e-5 K the tests of the steady state hold a cell to, which is itself ten times
/// the microkelvin a round has to move by less than. Below it, the field is refused.
const BELOW_ZERO: f64 = 1e-4;

/// One closed part of a block, as [`Solid3D::steady_state`] sees it — see `Solid3D::parts`.
struct Part {
    /// Whether any of its cells can lose heat.
    loses: bool,
    /// The watts generated in it, net.
    heat: f64,
    /// The watts generated in it with every source counted positive: the scale its net is
    /// rounded on.
    gross: f64,
    /// Whether any of its cells is a sink, generating negative watts.
    sinks: bool,
    /// How many cells it has.
    cells: usize,
    /// How many bodies it has.
    bodies: usize,
    /// Its lowest cell's index.
    first: usize,
}

impl Part {
    /// Whether it generates no heat net, to the rounding of spreading the watts over its cells
    /// and adding them back.
    ///
    /// `dissipating` stores `watts/n` in each cell, one rounding of at most `u·|each|`, and adding
    /// `n` of them back rounds by at most `(n − 1)·u·Σ|each|`, with `u = ε/2`: together
    /// `n·u·gross`, which `n·ε·gross` holds with a factor of two to spare. As much in as out,
    /// stated as two sources in one insulated part, is a balance that exists, and its net is
    /// rarely an exact zero once it is spread.
    fn balanced(&self) -> bool {
        self.heat.abs() <= self.cells as f64 * f64::EPSILON * self.gross
    }
}

/// One outer face of a block.
///
/// Named rather than indexed because a caller says which side of a part is exposed, and
/// `faces[3]` is a thing nobody can read back. The axis pairs are low and high along `x`, `y`
/// and `z` in the block's own coordinates — a [`Pose`](pantometry_core::Pose) is what puts those in
/// the world.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Face {
    /// The `x = 0` face.
    XMin,
    /// The far face along `x`.
    XMax,
    /// The `y = 0` face.
    YMin,
    /// The far face along `y`.
    YMax,
    /// The `z = 0` face.
    ZMin,
    /// The far face along `z`.
    ZMax,
}

/// Which of the three grid directions a face is normal to.
///
/// A [`Face`] names one of the six **outer** faces; this names an interior one, together with the
/// index of the cell on its high side. `Axis::X` with `i` is the face between `(i-1, j, k)` and
/// `(i, j, k)`, which is `kx`'s own indexing said in words rather than in arithmetic.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Axis {
    /// Faces normal to `x`.
    X,
    /// Faces normal to `y`.
    Y,
    /// Faces normal to `z`.
    Z,
}

/// The resolved seven-point operator, borrowed — see [`Solid3D::coefficients`].
///
/// Named rather than a tuple of five slices, because five slices of the same type in a row is a
/// signature where two can be swapped and nothing complains.
#[derive(Clone, Copy, Debug)]
pub struct Coefficients<'a> {
    /// Conductance of each face normal to x, `(nx+1)·ny·nz` of them.
    pub kx: &'a [f64],
    /// Conductance of each face normal to y, `nx·(ny+1)·nz`.
    pub ky: &'a [f64],
    /// Conductance of each face normal to z, `nx·ny·(nz+1)`.
    pub kz: &'a [f64],
    /// `dx/Cᵢ` for each cell, which turns a face flux into a rate of temperature. Zero in a void:
    /// nothing moves in a cell that holds nothing.
    pub mobility: &'a [f64],
    /// Watts a source puts in each cell. A cell emptied after its source was spread keeps its share
    /// here, until a later `dissipating` selects it and leaves it nothing, and generates none of
    /// it, which a kernel multiplying by the cell's zero `mobility` makes nothing of so long as the
    /// share is a number; the step refuses one that is not.
    pub source: &'a [f64],
}

impl Face {
    /// All six, in a fixed order — for a block exposed on every side.
    pub const ALL: [Face; 6] = [
        Face::XMin,
        Face::XMax,
        Face::YMin,
        Face::YMax,
        Face::ZMin,
        Face::ZMax,
    ];

    /// The two axis indices that lie *in* this face, in x-y-z order.
    ///
    /// A patch on `ZMin` is bounded in `[i, j]`, on `YMax` in `[i, k]`, on `XMin` in `[j, k]` —
    /// the same reading `ContactSpec` gives its own `from`/`to`, so a caller who has learned one
    /// has learned both.
    pub fn axes(&self) -> [usize; 2] {
        match self {
            Face::XMin | Face::XMax => [1, 2],
            Face::YMin | Face::YMax => [0, 2],
            Face::ZMin | Face::ZMax => [0, 1],
        }
    }

    /// Whether the cell at `(i, j, k)` of a `counts`-shaped block lies on this face.
    fn holds(&self, (i, j, k): (usize, usize, usize), counts: (usize, usize, usize)) -> bool {
        let (nx, ny, nz) = counts;
        match self {
            Face::XMin => i == 0,
            Face::XMax => i + 1 == nx,
            Face::YMin => j == 0,
            Face::YMax => j + 1 == ny,
            Face::ZMin => k == 0,
            Face::ZMax => k + 1 == nz,
        }
    }

    /// Every cell [`holds`](Face::holds) says lies on this face, once each, without visiting the
    /// ones it says do not — a plane of the block rather than all of it.
    fn cells(self, counts: (usize, usize, usize)) -> impl Iterator<Item = (usize, usize, usize)> {
        let (nx, ny, nz) = counts;
        let [a, b] = self.axes();
        let n = [nx, ny, nz];
        // An empty block has no face to hold anything. Without this a low face of a block with no
        // cells along x yields a cell at the origin, and a high face's `n - 1` underflows.
        let rows = if nx * ny * nz == 0 { 0 } else { n[b] };
        (0..rows).flat_map(move |q| {
            (0..n[a]).map(move |p| match self {
                Face::XMin => (0, p, q),
                Face::XMax => (nx - 1, p, q),
                Face::YMin => (p, 0, q),
                Face::YMax => (p, ny - 1, q),
                Face::ZMin => (p, q, 0),
                Face::ZMax => (p, q, nz - 1),
            })
        })
    }
}

/// the arithmetic one, and the difference is not a refinement away.
#[derive(Clone, Debug)]
pub struct Solid3D {
    name: String,
    /// Every substance in the block. Index 0 is the constructor's, and a block nobody has
    /// [`fill`](Solid3D::fill)ed holds only that one.
    materials: Vec<Substance>,
    /// Which of `materials` each cell is, indexed as `cells` is.
    which: Vec<u32>,
    /// Cell-centre temperatures, indexed `x + nx*(y + ny*z)`.
    cells: Vec<f64>,
    saved: Vec<f64>,
    counts: (usize, usize, usize),
    dx: Length,
    absorbed: f64,
    /// What [`stored_heat`](Solid3D::stored_heat) is measured from — see [`Bar1D`](crate::Bar1D)
    /// for why an enthalpy reference is chosen for precision rather than for physics.
    reference: f64,
    /// Face conductivities in W/m/K, rebuilt by [`resolve`](Solid3D::resolve): `kx` is indexed
    /// `i + (nx+1)*(j + ny*k)` and holds the face between cells `i-1` and `i`, so `kx[0]` and
    /// `kx[nx]` are the outer faces and are **zero** — which is the insulated boundary, stated in
    /// the conductance rather than in the stencil's index arithmetic.
    kx: Vec<f64>,
    ky: Vec<f64>,
    kz: Vec<f64>,
    /// Contact conductance on each face, W·m⁻²·K⁻¹, indexed exactly as `kx`/`ky`/`kz` are — see
    /// [`Solid3D::joined`]. `+∞` is a face with nothing stated, which is every face of a block
    /// nobody has joined.
    ///
    /// **Empty until a contact is placed**, and read only when it is not. A block with no joints
    /// carries no array and pays nothing, which matters because these are the same length as the
    /// face conductivities and those are the largest arrays here after the cells themselves.
    hx: Vec<f64>,
    hy: Vec<f64>,
    hz: Vec<f64>,
    /// The first contact conductance that cannot be used, if any was stated — negative, or not a
    /// number. Kept so the refusal can **name what happened**.
    ///
    /// Without it a bad joint takes the same road a substance with no conductivity takes, and the
    /// sweep says "substance has no diffusivity" to a caller whose substances are all fine. That
    /// message sends them to their materials, which is the one place the fault is not.
    bad_contact: Option<f64>,
    /// Per-cell `dx / C_i`, which is what multiplies `Σ k_f ΔT` to give a rate of temperature.
    mobility: Vec<f64>,
    /// Per-cell heat capacity `ρ_i c_i dx³`, in J/K. `NaN` where a substance does not say.
    capacity: Vec<f64>,
    /// Per-cell melting point in kelvin, or `+∞` for a substance that does not melt — so the phase
    /// branch is a comparison that is simply never true rather than an `Option` to unwrap per cell.
    melt_point: Vec<f64>,
    /// Per-cell latent heat in **joules**: `ρ L dx³`. Zero where there is no phase change, which is
    /// what the sweep tests.
    ///
    /// # It was kelvin, and two phases broke that
    ///
    /// `L/c_p` — 163 K for ice — kept the whole update in one unit and the numbers `O(1)` rather than
    /// `O(3e5)`, which is the reference-point argument `runtime/gpu` paid 1660× for. It is correct
    /// while a cell has **one** specific heat.
    ///
    /// With two it is not, and it fails quietly. The kelvin figure is normalised by the *solid's*
    /// `c_p` and then multiplied by the cell's *current* capacity, which for a liquid cell is the
    /// liquid's — so freezing a cell of water cost `4182/2050` = **2.04×** what it should, and the
    /// front came out **27% short** with nothing pointing at the latent heat. Joules do not have a
    /// normalisation to get wrong.
    latent: Vec<f64>,
    /// Per-cell capacity of the **solid** phase, `ρ c_s dx³`, in J/K.
    ///
    /// The enthalpy map needs both capacities by name: `C_s` below the melting point and `C_l` above
    /// it. `capacity` is the *mixed* one and is what the sweep divides a flux by, which is right —
    /// a fully solid or fully liquid cell has its own, and a mushy cell does not change temperature
    /// at all, so what it would have divided by never matters.
    cap_solid: Vec<f64>,
    /// Per-cell capacity of the **liquid** phase, `ρ c_l dx³`. Equal to `cap_solid` under the
    /// one-phase model.
    cap_liquid: Vec<f64>,
    /// How much of each cell has melted, `0` solid and `1` liquid. Derived from the temperature by
    /// [`resolve`](Solid3D::resolve) and then carried by the sweep, and by the mutators between
    /// steps, each of which re-derives what a fraction it moved decides — see
    /// [`rederive_cell`](Solid3D::rederive_cell).
    melted: Vec<f64>,
    /// The phase at the last [`checkpoint`](Domain::checkpoint), alongside `saved`.
    saved_melted: Vec<f64>,
    /// `max_i dx·Σ_f k_f / C_i`, the reciprocal of the largest stable step. See
    /// [`max_stable_dt`](Domain::max_stable_dt) for why it is a maximum over cells and not a
    /// function of the fastest material.
    ///
    /// Read only through [`conduction_rate`](Solid3D::conduction_rate), because between a
    /// mutator moving a phase and the next `resolve` it is not the answer.
    worst_rate: f64,
    /// Whether `worst_rate` is out of date: set when a mutator between steps moved one cell's
    /// phase and re-derived that cell, cleared by every `resolve`.
    ///
    /// The cell's faces and sums are put right on the spot, and the maximum is the one thing that
    /// is not, because finding it again is a pass over every cell, and a caller setting a region
    /// one cell at a time would pay a pass per cell. A full re-derivation per cell would have been
    /// worse: 6.2 ms on a two-phase block at 64³, so 27 minutes to declare every cell of one that
    /// way. So the pass is left to whoever next asks for the limit.
    worst_stale: bool,
    /// Per substance, the solid phase's and the liquid phase's conductivity, capacity per cell and
    /// emissivity, `((k_s, C_s, ε_s), (k_l, C_l, ε_l))` — the triples a cell's phase mixes between.
    /// Kept from `resolve` so one cell can be re-derived without rebuilding the rest.
    phase_props: Vec<(Phase, Phase)>,
    /// Each cell's conductivity, capacity and emissivity **as the operator was last derived**,
    /// `(k, C, ε)`: what `resolve` mixed them to, and `rederive_cell` since, and zero for void.
    ///
    /// The film reads its half cell's conductivity and its surface's emissivity here, and
    /// `find_gaps` the emissivities of a clearance's two ends, rather than mixing them again from
    /// `melted`: the faces, the capacities and the limit are built from these, and a reader that
    /// mixes for itself agrees with them only while nothing moves a phase between the limit and
    /// the sweep. Something did. The heat a step takes off the bus went in first, after the limit
    /// and before the sweep, and mixed again from `melted`, as the first version of the
    /// conductivity's fix did, a liquid conducting ten times its solid, melted by heat off the
    /// bus within a step at its limit, carried 7.2 times the conductance the limit had charged,
    /// and took out 0.272 J in that step, refreezing 83% of the cell. That heat goes in after the
    /// sweep now — see `step` — so for the film and the limit the two agree again. Not for every
    /// reader: the clearances' pass runs after the sweep, the source and the bus heat have all
    /// moved phases, and reads its half cells here, which is what keeps them the phase the step
    /// began in.
    derived: Vec<Phase>,
    /// Cells that are **not part of the block** — see [`Solid3D::empty`].
    void: Vec<bool>,
    /// How many solid cells lie on each face, in `Face` order. Filled by `resolve`; see
    /// [`Solid3D::cells_on`] for why it is a cache and not a count.
    solid_on: [usize; 6],
    /// The solid cells with at least one face on the void, in index order, when the block has
    /// [air](Solid3D::air_in) — the only solid cells off the exposed faces whose loss can be
    /// anything but zero. Empty without air. Filled by `resolve`, which every mutator of the void
    /// or the air calls, so the cells it lists are the ones the void has now — see
    /// `worst_rate_now`, which reads it.
    facing_air: Vec<(usize, usize, usize)>,
    /// Pairs of solid cells that face each other across a run of void, with the radiative
    /// exchange coefficient between them. Rebuilt whenever the void or the materials change, and,
    /// for a substance whose phases radiate differently, whenever a phase does — by `resolve`
    /// after a step, and by `regap` after a mutator.
    ///
    /// `(a, b, coefficient)`, `a < b`, where the coefficient is `σA/(1/ε₁ + 1/ε₂ − 1)` with each
    /// `ε` its cell's phase's, so the exchange is
    /// `coefficient · (Tₐ⁴ − T_b⁴)` — see [`Solid3D::empty`] for what this models and what it
    /// does not.
    gaps: Vec<(usize, usize, f64)>,
    /// Every cell that is an end of a clearance, in index order, with where its pairs are listed
    /// in `gap_end_pairs`. Built by `find_gaps` beside `gaps`, so that the limit can add up
    /// everything a cell radiates rather than charge each pair alone; `regap` changes coefficients
    /// and never which pairs there are, so it leaves this as it is.
    gap_ends: Vec<GapEnd>,
    /// Indices into `gaps`, each end's in pair order — see `gap_ends`.
    gap_end_pairs: Vec<usize>,
    /// Each cell's conduction sum, `Σ_f k_f`, kept so the stability limit can be rebuilt at
    /// the block's **current** temperature — see [`Solid3D::worst_rate_now`].
    face_sum: Vec<f64>,
    /// Which faces lose heat, and to what. Empty is the adiabatic block every scene had until
    /// now: six insulated faces and no steady state to reach.
    exposed: BTreeMap<Face, Environment>,
    /// A contact conductance between an outer face and whatever cools it, W·m⁻²·K⁻¹ — see
    /// [`Solid3D::mounted_on`]. Absent is a face bolted to its film with nothing in between.
    ///
    /// Its own map rather than a field on [`Environment`], because that struct is built by
    /// literal in nineteen places across this workspace and a joint is a property of the
    /// mounting rather than of the air.
    mounted: BTreeMap<Face, f64>,
    /// The box on a face that is actually exposed, `[lo0, lo1, hi0, hi1]` in that face's two axes
    /// — see [`Solid3D::losing_from_within`]. Absent is the whole face.
    ///
    /// Half-open, like every other pair of bounds in this workspace: a single cell is
    /// `hi = lo + 1`.
    patch: BTreeMap<Face, [usize; 4]>,
    /// The air filling this block's void, if it has any: `(ambient in kelvin, W·m⁻²·K⁻¹)`.
    ///
    /// Absent is a **vacuum**, which is what every void in this crate was until this existed — see
    /// [`Solid3D::air_in`] for what that cost.
    air: Option<(f64, f64)>,
    /// Joules given up to those environments over the run, and the reason the books still
    /// balance: the ledger is `stored + lost`, so what leaves the cells arrives in this
    /// counter and the total moves only by what crossed the bus. `LumpedMass` keeps the same
    /// pair for the same reason.
    lost: f64,
    /// Watts generated in each cell, held per cell rather than per region so a source and a
    /// material can be cut by different boxes without either knowing about the other.
    ///
    /// **The gap this closes is one the bus cannot**, by design: the plain channel carries an
    /// amount and no location, so heat arriving there spreads to a uniform rise over everything
    /// that can hold it. That is the only choice that adds no information, and it is the wrong
    /// answer for a die, a winding, a brake disc or a laser absorber — every real thing that
    /// dissipates does it *somewhere*. Symmetric with [`Solid3D::losing_from`], which takes energy
    /// out at a face; this puts it in at a region.
    source: Vec<f64>,
    /// The first cell whose source is not a number of watts, and that number, if any is — found
    /// by [`Solid3D::dissipating`], the one place a source is written, so a step pays nothing to
    /// refuse it. See `input_fault`.
    bad_source: Option<(usize, f64)>,
    /// The cells whose source takes heat out, in order: with the bus, the only way a step within
    /// its limit can take a cell below absolute zero. Found by [`Solid3D::dissipating`].
    sinks: Vec<usize>,
    /// Whether a cell's temperature was set since a step last read every one — by the
    /// constructor, `set_temperature`, `deposit`, `set_melted_fraction`, `release_mode` or a
    /// restore. A step reads the cells only when one was, so a block nobody touched between steps
    /// pays nothing to be asked.
    unchecked: bool,
    /// The fractions a step that takes heat out began with, to put back if it is refused.
    undo_melted: Vec<f64>,
    /// Joules generated over the run, and the reason the books still balance: the ledger is
    /// `stored + lost − supplied`, so a source that adds a joule to a cell adds one here too and
    /// the total moves only by what crossed the bus.
    supplied: f64,
    /// The saved counterpart of [`Solid3D::supplied`], for the same reason [`Solid3D::lost`] has
    /// one: an iterative sweep that rewinds the cells and not the counter would credit itself a
    /// sweep of generated heat per retry.
    saved_supplied: f64,
    /// The saved counterpart of [`Solid3D::lost`].
    ///
    /// Saved because `LumpedMass` learned this the expensive way: rewinding the cells and not
    /// the losses makes an iterative sweep grow its books by one sweep of shed heat per
    /// iteration, and the audit reports energy created from nothing. Measured there at 1567.6 J
    /// becoming 1600.9 J over forty advances.
    saved_lost: f64,
    /// The saved counterpart of what [`Solid3D::absorbed_energy`] reads. Not in the books, which
    /// are the cells, so nothing refused its absence: a block fed a hundredth of a joule a step
    /// read 0.03 J absorbed straight after being restored to a checkpoint at which it had taken
    /// none, and an iterated window counted every sweep it threw away.
    saved_absorbed: f64,
    /// Whether any material names a **liquid** phase whose properties differ from its solid: its
    /// conductivity, its specific heat or its emissivity.
    ///
    /// The one thing that decides whether the sweep has to rebuild its operator each step. A block of
    /// ice does; a block of aluminium, or of ice under the one-phase model, does not.
    two_phase: bool,
    /// Below this, something a cell exchanges with can freeze it in one step further than the
    /// capacity it has bounds: `T_m − L/(C_L − C_S)`, the highest of any material whose liquid holds
    /// more than its solid, and `−∞` where none does. See `charge_crossing`.
    freeze_below: f64,
    /// Above this, something a cell exchanges with can melt it in one step further than the
    /// capacity it has bounds: `T_m + L/(C_S − C_L)`, the lowest of any material whose solid holds
    /// more than its liquid, and `+∞` where none does.
    melt_above: f64,
    /// The coldest and the hottest solid cell, as `resolve` last found them, and `None` once
    /// anything has moved a cell since, or for a block whose bounds nothing can reach, which never
    /// reads it. Kept so that a two-phase block nowhere near either bound
    /// above does not pass over its cells twice a step to learn so: measured, that pass cost the
    /// limit of a 16³ block of water 6.8 µs of its 10.6, and its step 9%.
    cells_range: Option<(f64, f64)>,
    /// The previous state the sweep reads, kept between steps so the step does not allocate.
    ///
    /// `let old = self.cells.clone()` allocated and copied the whole grid every step. At 64³ that
    /// is 2 MB of memcpy against a 1.5 ms step, and at 128³ it is 16 MB — measurable, and free to
    /// remove once there is somewhere to put it.
    old: Vec<f64>,
    /// Whether any cell converts its rise through an enthalpy — see `resolve`.
    latent_anywhere: bool,
    /// The buffer the double-buffered sweep writes into, which becomes `cells` and hands back the
    /// array it replaced. Kept between steps so a step allocates nothing.
    ///
    /// Two earlier shapes cost more than the threading they enabled. Carrying the *rise* rather
    /// than the new value meant a third array written and read back; carrying the shed joules
    /// beside it meant sixteen bytes a cell instead of eight. The stencil is bandwidth-bound, so
    /// both showed up immediately: the first was slower than the sequential original at every grid
    /// below 96³.
    rise: Vec<f64>,
    /// What the film takes from each cell this step, in the sweep's units: written by the pass that
    /// totals it for the books and read by the sweep, so that a film is evaluated once a step and
    /// not twice. Only the cells `sheds` names are written or read.
    film: Vec<f64>,
}

/// One patch of a clearance: two facing surfaces and how far apart they are.
///
/// A gap is found pair by pair, along one axis, one cell face against the one directly across
/// from it. A *patch* is what those pairs add up to — a connected sheet of facing area at one
/// separation, which is the shape a view factor is a statement about.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GapPatch {
    /// How many facing cell pairs are in it.
    pub pairs: usize,
    /// The patch's extent across the gap, in metres — the bounding box of the facing cells.
    pub span: (f64, f64),
    /// How far the two surfaces are apart, in metres.
    pub distance: f64,
    /// Whether the facing cells fill their bounding box.
    ///
    /// [`view_factor`](GapPatch::view_factor) is exact for a rectangle and an **upper bound** for
    /// anything else, because a patch with a bite out of it sees less of itself than its bounding
    /// box does. Reported rather than corrected: an L-shaped clearance is not a shape the closed
    /// form covers, and a number invented for it would be worth less than knowing it is a bound.
    pub rectangular: bool,
}

impl GapPatch {
    /// The view factor from one of the facing surfaces to the other, `F₁₂`.
    ///
    /// The exact closed form for two equal, parallel, directly-opposed rectangles — which is the
    /// geometry the pairing produces by construction, since a pair is a cell and the cell across
    /// from it. With `X = a/c` and `Y = b/c`:
    ///
    /// ```text
    /// F = 2/(πXY) · [ ln √((1+X²)(1+Y²)/(1+X²+Y²))
    ///                 + X√(1+Y²)·atan(X/√(1+Y²)) + Y√(1+X²)·atan(Y/√(1+X²))
    ///                 − X·atan X − Y·atan Y ]
    /// ```
    ///
    /// **This is not what the exchange is charged**, and the difference is the point. What the
    /// block charges is `F̄ = 1`, which is exact when the sides of the gap are mirrors — and the
    /// block's own outer faces are exactly that, because an insulated boundary is implemented as
    /// a mirror and a mirror extends the two surfaces to infinity. `F₁₂` is what the same pair
    /// would exchange with **nothing** at the sides, open to space. The two agree when the gap is
    /// narrow compared with the surfaces and diverge without limit when it is not, so this is the
    /// number that says how much of the answer is resting on reading the boundary one way.
    pub fn view_factor(&self) -> f64 {
        let (x, y) = (self.span.0 / self.distance, self.span.1 / self.distance);
        // A patch with no extent sees nothing; two surfaces with nothing between them see each
        // other entirely. Both are limits of the form rather than sentinels, and both are here
        // because this is public: a caller can build a `GapPatch` this block would never produce.
        if !self.distance.is_finite() || self.span.0.is_nan() || self.span.1.is_nan() {
            return 0.0;
        }
        if self.distance <= 0.0 {
            return if self.span.0 > 0.0 && self.span.1 > 0.0 {
                1.0
            } else {
                0.0
            };
        }
        if x <= 0.0 || y <= 0.0 {
            return 0.0;
        }
        let (x2, y2) = (x * x, y * y);
        // **The closed form cancels itself away for a distant patch.** The bracket below is a
        // small difference of terms of order `X`, and the answer is of order `X²` — so at `X` of
        // 1e-3 four digits are already gone and the ratio to the exact `A/πc²` limit *turns
        // around*: measured 0.999993 at X = 1/300, 1.00009 at 1/1000 and 1.007 at 1/3000, growing
        // where it should be converging. A function that is quietly worse the further apart two
        // surfaces are is the shape of defect this repository looks for, so the small-argument
        // regime gets the series instead, `F = XY/π · (1 − (X²+Y²)/3 + O(X⁴))`, which is exact to
        // 1e-8 either side of the crossover.
        if x < 0.01 && y < 0.01 {
            return (x * y / std::f64::consts::PI * (1.0 - (x2 + y2) / 3.0)).clamp(0.0, 1.0);
        }
        let (rx, ry) = ((1.0 + x2).sqrt(), (1.0 + y2).sqrt());
        let t = ((1.0 + x2) * (1.0 + y2) / (1.0 + x2 + y2)).sqrt().ln()
            + x * ry * (x / ry).atan()
            + y * rx * (y / rx).atan()
            - x * x.atan()
            - y * y.atan();
        (2.0 / (std::f64::consts::PI * x * y) * t).clamp(0.0, 1.0)
    }
}

/// The phase a temperature names: liquid above the melting point, and at it or below, solid.
///
/// A temperature is not a state for a substance that melts — 0 °C is ice, water, or any mixture of
/// the two — so naming one has to name a phase too, and the choice is the unmelted one at the point
/// itself: a caller clamping a face below freezing means ice, and a caller clamping it above means
/// liquid. One rule for the constructor, `fill` and `set_temperature`.
fn phase_named_by(t: f64, point: f64) -> f64 {
    if t > point {
        1.0
    } else {
        0.0
    }
}

/// A running maximum of stability rates, raised by `rate`. `NaN` once any rate is.
///
/// NaN loses every comparison, so `max` alone would step straight past a substance that does not
/// say what it conducts, find the other cells' rates, and let the sweep write `NaN` into the cells
/// that have none. `a_block_that_cannot_be_stepped_still_cannot_after_a_phase_moves` holds the
/// search after a mutator to the refusal the derivation gives.
fn worse(worst: f64, rate: f64) -> f64 {
    if worst.is_nan() || rate.is_nan() {
        f64::NAN
    } else {
        worst.max(rate)
    }
}

/// `(lo, hi)` widened to take in `t`, by comparison: a `NaN` widens nothing, as `f64::min` and
/// `max` have it, without the instructions they spend on saying so.
#[inline]
fn widened((lo, hi): (f64, f64), t: f64) -> (f64, f64) {
    (if t < lo { t } else { lo }, if t > hi { t } else { hi })
}

/// What one phase of a substance conducts, holds per cell and radiates with, `(k, C, ε)`.
type Phase = (f64, f64, f64);

/// One cell at the end of a clearance, as the limit reads it: where it is, and where its pairs
/// are listed in `Solid3D::gap_end_pairs`.
///
/// Its coordinates are kept rather than divided out of its index each time, because the limit is
/// asked for every step and three integer divisions an end were most of what reading the ends
/// cost: measured against the limit's old loop, a step of a block with a clearance through its
/// middle 3% to 4% slower, 16³ to 64³, pinned, and within the noise of it with them kept.
#[derive(Clone, Debug)]
struct GapEnd {
    cell: usize,
    at: (usize, usize, usize),
    /// On the block's own boundary, the only place an exposed face can be.
    outer: bool,
    first: usize,
    count: usize,
}

/// A cell's conductivity, capacity and emissivity, `phi` of the way from its substance's solid
/// triple to its liquid triple — see [`Solid3D::mixed_at`].
///
/// **The emissivity is each phase's own where the cell is wholly that phase**, rather than the
/// same line evaluated at an end. The line at an end reads the other phase as well: a liquid whose
/// emissivity is not a number made `0.9 + 0·(NaN − 0.9)` of every solid cell, so a block of solid
/// that had never melted lost a `NaN` through its film, and before the liquid's emissivity was
/// read at all those cells were right. And it misses the end: `0.9 + (0.1 − 0.9)` is two ulps
/// under `0.1`. The conductivity and the capacity are left on the line, so a block that radiates
/// as it did is stepped to the bit as it was.
fn mix((solid, liquid): (Phase, Phase), phi: f64) -> Phase {
    (
        solid.0 + phi * (liquid.0 - solid.0),
        solid.1 + phi * (liquid.1 - solid.1),
        if phi <= 0.0 {
            solid.2
        } else if phi >= 1.0 {
            liquid.2
        } else {
            solid.2 + phi * (liquid.2 - solid.2)
        },
    )
}

/// The mean size of the solid cells' temperatures, K: what [`Solid3D::steady_state`] holds a body's
/// step to a quarter of, for a body at absolute zero that a quarter of its own level would hold
/// still. Added in index order, so for a block of one body whose cells are all above zero it is
/// that body's mean to the bit.
fn cell_level(t: &[f64], body: &[Option<usize>]) -> f64 {
    let (mut sum, mut cells) = (0.0, 0.0_f64);
    for (v, b) in t.iter().zip(body) {
        if b.is_some() {
            sum += v.abs();
            cells += 1.0;
        }
    }
    sum / cells.max(1.0)
}

/// The radiative exchange coefficient between two faces of area `area` facing each other, with
/// emissivities `a` and `b`: `σA/(1/ε₁ + 1/ε₂ − 1)`, the parallel-plate series — or `None` where
/// the pair carries nothing. A surface that does not radiate at all makes the pair carry nothing,
/// which is the right limit and also keeps the reciprocal finite.
///
/// Shared by `find_gaps`, which pairs a block's clearances, and `regap`, which puts one cell's
/// pairs right after its phase moved, so the two agree to the bit.
fn exchange(a: f64, b: f64, area: f64) -> Option<f64> {
    if a > 0.0 && b > 0.0 {
        let resistance = 1.0 / a + 1.0 / b - 1.0;
        if resistance > 0.0 {
            return Some(STEFAN_BOLTZMANN.to_si() * area / resistance);
        }
    }
    None
}

/// The series mean of the two half cells either side of a face. Written as `2ab/(a+b)` and guarded
/// at zero, because a perfect insulator is a legitimate fill and `0/0` is not a boundary condition.
fn series(a: f64, b: f64) -> f64 {
    let sum = a + b;
    if sum > 0.0 {
        2.0 * a * b / sum
    } else {
        0.0
    }
}

/// **A stated contact, in series with the two half cells.** `k` is the resistance
/// `dx/(2k_L) + dx/(2k_R)` written as a conductivity; a contact adds `1/h` to it, so
///
/// ```text
/// 1/k_face = 1/k_series + 1/(dx·h)
/// ```
///
/// `h = +∞` is the unstated face and returns `k` bit for bit, which is why an untouched block
/// cannot be perturbed by this existing. `h = 0` is a clearance and carries nothing. A `NaN`
/// conductance makes the face `NaN`, which reaches the limit and refuses the sweep — the same road
/// a substance that does not say what it conducts already takes.
fn contact(k: f64, h: f64, dx: f64) -> f64 {
    if h == f64::INFINITY {
        return k;
    }
    // A **negative** conductance is not a poor joint, it is a mistake: a face that would carry
    // heat up its own gradient. `NaN` reaches the limit and refuses the sweep, which is the road a
    // substance that does not say what it conducts already takes. Clamping it to an insulator
    // instead would model something nobody asked for and never say so — and an insulator is a
    // plausible enough answer to go unnoticed.
    if h.is_nan() || h < 0.0 {
        return f64::NAN;
    }
    // `h = 0` needs no case of its own: `1/(dx·0)` is `+∞` and `1/∞` is zero, which is the
    // clearance. Guarded on `k` because a void's conductivity is zero and `0/0` is not a boundary
    // condition.
    if k <= 0.0 {
        return 0.0;
    }
    1.0 / (1.0 / k + 1.0 / (dx * h))
}

impl Solid3D {
    /// A block of `counts` cubic cells of side `dx`, all starting at `initial`.
    ///
    /// Each count is forced to at least one. A block one cell thick in two directions is a
    /// legitimate thing to ask for and reduces exactly to a bar, which is how the closed-form
    /// tests check the three axes against each other.
    pub fn new(
        name: impl Into<String>,
        substance: Substance,
        counts: (usize, usize, usize),
        dx: Length,
        initial: Temperature,
    ) -> Solid3D {
        let counts = (counts.0.max(1), counts.1.max(1), counts.2.max(1));
        let cells = vec![initial.to_si(); counts.0 * counts.1 * counts.2];
        let mut block = Solid3D {
            name: name.into(),
            materials: vec![substance],
            which: vec![0; cells.len()],
            saved: cells.clone(),
            cells,
            counts,
            dx,
            absorbed: 0.0,
            reference: initial.to_si(),
            kx: Vec::new(),
            ky: Vec::new(),
            kz: Vec::new(),
            hx: Vec::new(),
            hy: Vec::new(),
            hz: Vec::new(),
            bad_contact: None,
            mobility: Vec::new(),
            capacity: Vec::new(),
            melt_point: Vec::new(),
            latent: Vec::new(),
            melted: Vec::new(),
            saved_melted: Vec::new(),
            cap_solid: Vec::new(),
            cap_liquid: Vec::new(),
            worst_rate: 0.0,
            worst_stale: false,
            phase_props: Vec::new(),
            derived: Vec::new(),
            void: vec![false; counts.0 * counts.1 * counts.2],
            solid_on: [0; 6],
            facing_air: Vec::new(),
            source: vec![0.0; counts.0 * counts.1 * counts.2],
            bad_source: None,
            sinks: Vec::new(),
            unchecked: true,
            undo_melted: Vec::new(),
            supplied: 0.0,
            saved_supplied: 0.0,
            gaps: Vec::new(),
            gap_ends: Vec::new(),
            gap_end_pairs: Vec::new(),
            face_sum: Vec::new(),
            exposed: BTreeMap::new(),
            mounted: BTreeMap::new(),
            patch: BTreeMap::new(),
            air: None,
            lost: 0.0,
            saved_lost: 0.0,
            saved_absorbed: 0.0,
            two_phase: false,
            freeze_below: f64::NEG_INFINITY,
            melt_above: f64::INFINITY,
            cells_range: None,
            old: Vec::new(),
            latent_anywhere: false,
            rise: Vec::new(),
            film: Vec::new(),
        };
        block.resolve();
        block
    }

    /// How much of one cell has melted: `0` entirely solid, `1` entirely liquid.
    ///
    /// Always zero for a substance with no [`FusionProps`](pantometry_core::substance::FusionProps), which is every
    /// entry in the catalogue but [`Substance::ice`].
    pub fn melted_fraction_at(&self, i: usize, j: usize, k: usize) -> f64 {
        let (nx, ny, nz) = self.counts;
        self.melted[i.min(nx - 1) + nx * (j.min(ny - 1) + ny * k.min(nz - 1))]
    }

    /// Declare how much of a cell has melted, for an initial condition a temperature cannot express.
    ///
    /// The counterpart to [`set_temperature`](Solid3D::set_temperature), and it exists because **a
    /// temperature is not a state** for a substance that melts: 0 °C is ice, water, or any mixture,
    /// and Stefan's problem starts from liquid at exactly the melting point. Without this the only
    /// way to say that is to start a hair above it, which puts sensible heat in the initial condition
    /// that the closed form does not have.
    ///
    /// Ignored for a substance that does not melt, and out of range is ignored, matching
    /// `set_temperature`.
    ///
    /// # The next step conducts through the phase declared here
    ///
    /// For a substance with a [liquid](pantometry_core::substance::FusionProps::liquid) phase, a
    /// cell's conductivity and capacity follow its fraction, so declaring one changes the
    /// operator, and this re-derives the cell and the faces it shares before it returns; the limit
    /// is found again by whatever next asks for it. It did not: the operator the constructor had
    /// derived from the temperature stayed
    /// until the end of the first step. Ice and water, built as water at 5 °C and declared ice at
    /// its melting point, reported a limit of 3.206 s for a block whose own was 0.423 s, and one
    /// step at the limit it reported took a cell set to −20 °C up to 0 °C and its neighbour down
    /// to −20.4 °C — past each other and out of the range the block had started in, with the
    /// books exact.
    ///
    /// # Supercooling is not representable, and this keeps it that way
    ///
    /// The state is one monotone number, so a fraction and a temperature cannot disagree: asking for
    /// liquid raises the temperature to at least the melting point and asking for solid lowers it to
    /// at most. Liquid below freezing and solid above it are real states of real matter and this
    /// model does not have them — the sharp-interface problem assumes the interface is *at* the
    /// melting point, which is what makes Neumann's solution its solution.
    pub fn set_melted_fraction(&mut self, i: usize, j: usize, k: usize, fraction: f64) {
        let Some(idx) = self.index(i, j, k) else {
            return;
        };
        if self.latent[idx] <= 0.0 {
            return;
        }
        self.cells_moved();
        let before = self.melted[idx];
        let phi = fraction.clamp(0.0, 1.0);
        self.melted[idx] = phi;
        let point = self.melt_point[idx];
        if phi >= 1.0 {
            self.cells[idx] = self.cells[idx].max(point);
        } else if phi <= 0.0 {
            self.cells[idx] = self.cells[idx].min(point);
        } else {
            self.cells[idx] = point;
        }
        self.rephased(idx, before);
    }

    /// The volume that has melted, summed over every cell.
    ///
    /// The integral quantity, and the one worth reading rather than a front *position*: a front is
    /// only a position if the problem is one-dimensional, and this is the same number in three. For
    /// a column of cells it gives the position anyway — `melted_volume / area` is how far in the
    /// front has reached, including the partial cell it is currently inside, which is what makes a
    /// front measurable to better than one cell.
    pub fn melted_volume(&self) -> Volume {
        Volume::from_si(self.melted.iter().sum::<f64>() * self.cell_volume())
    }

    /// Put a different substance in every cell the predicate accepts.
    ///
    /// Composable: call it once per layer, per coating, per inclusion. Cells nobody claims keep the
    /// constructor's substance, so a block is never partly undefined.
    ///
    /// # The harmonic mean, and why it is not a detail
    ///
    /// Heat crosses a face, and a face has two materials touching it. The conductance of the half
    /// cell either side is in **series**, so the conductivity that governs the face is
    ///
    /// ```text
    ///   k_face = 2 k_L k_R / (k_L + k_R)
    /// ```
    ///
    /// the harmonic mean, and the arithmetic mean `(k_L + k_R)/2` is not an alternative convention
    /// — it is wrong. With aluminium against borosilicate the two differ by **38×** (2.21 against
    /// 84.1 W/m/K), and the arithmetic one short-circuits the interface: a wall it models has 4.2%
    /// less resistance than its own layers add up to at 24 cells, and reaching the harmonic answer
    /// to 0.1% would take about a thousand.
    ///
    /// The harmonic mean is not merely better. With the material interface on a cell face it makes
    /// the discrete series resistance **exactly** `Σ Lᵢ/(kᵢA)` at every resolution, which is why
    /// `a_layered_wall.rs` is an equality and not a tolerance.
    ///
    /// "On a cell face" is the whole condition, and it is one this method cannot break: cells are
    /// whole, so an interface is always on a face. A scheme that placed a layer boundary partway
    /// through a cell would be first order there whichever mean it used, because the cell would be a
    /// mixture and no single conductivity describes one.
    ///
    /// # What it costs
    ///
    /// Five `f64` and a `u32` per cell beyond the temperatures: three face conductivities, a
    /// capacity, its reciprocal scaled by `dx`, and which substance the cell is. All rebuilt here
    /// and none of it in the sweep, which is the trade — the alternative is six harmonic means and
    /// six divisions per cell per step.
    ///
    /// A uniform block pays the same. That is deliberate: a fast path for one material would be a
    /// second implementation of the same physics, exercised by the tests that came before this
    /// method and by nothing after it.
    ///
    /// # A cell that changes substance takes the phase its temperature names
    ///
    /// As it would have in a block built of that substance, and as
    /// [`set_temperature`](Solid3D::set_temperature) names one. A cell that keeps its substance
    /// keeps its phase, so a block that has been melting can be filled again elsewhere without
    /// losing its front.
    ///
    /// It used to carry the old substance's fraction into the new one. Aluminium at 40 °C with a
    /// cell filled with a wax that melts at 28.15 °C held the wax as a solid at 40 °C, opened its
    /// books 0.199 J short — that cell's latent heat — and, with nothing flowing, came out of one
    /// step at 28.15 °C and 9.4% melted.
    pub fn fill(
        &mut self,
        substance: Substance,
        which: impl Fn(usize, usize, usize) -> bool,
    ) -> &mut Solid3D {
        let id = match self.materials.iter().position(|s| *s == substance) {
            Some(at) => at as u32,
            None => {
                self.materials.push(substance);
                (self.materials.len() - 1) as u32
            }
        };
        let (nx, ny, nz) = self.counts;
        let mut rephase = vec![false; self.cells.len()];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if which(i, j, k) {
                        let c = i + nx * (j + ny * k);
                        // New matter in the cell, so a fraction carried from what was there
                        // before is a statement about something else — or, from void, about
                        // nothing.
                        rephase[c] = self.which[c] != id || self.void[c];
                        self.which[c] = id;
                        // **Filling a cell puts something in it**, which includes a cell that
                        // held nothing. `empty` was one-way until this line, and `pantometry-
                        // world`'s own comment beside the `void` region says the opposite in as
                        // many words — "applied in the same order as any other region, later wins
                        // where they overlap, so a gap cut into a part reads the way a coating on
                        // a layer does". A scene written that way was refused for generating heat
                        // in cells it had just filled, which is the loud version; the quiet one is
                        // a part with a hole in it that nobody put there.
                        self.void[c] = false;
                    }
                }
            }
        }
        // A void cell filled is solid again, holding what it held as void, which nothing asked: a
        // cell set to −50 K while void, filled after a step, stepped on from −50 K.
        self.cells_moved();
        self.resolve_rephasing(&rephase);
        self
    }

    /// The substance in one cell. Out of range reads the nearest one in range.
    pub fn substance_at(&self, i: usize, j: usize, k: usize) -> &Substance {
        let (nx, ny, nz) = self.counts;
        let idx = i.min(nx - 1) + nx * (j.min(ny - 1) + ny * k.min(nz - 1));
        &self.materials[self.which[idx] as usize]
    }

    /// How many distinct substances are in the block. One until something has [`fill`]ed it.
    ///
    /// [`fill`]: Solid3D::fill
    pub fn substances(&self) -> usize {
        self.materials.len()
    }

    /// The resolved coefficients of the seven-point stencil, for an accelerator that wants to run
    /// **this** operator rather than one of its own.
    ///
    /// See [`Coefficients`] for what each array is and how long it is.
    ///
    /// # Why an accelerator should take these and not the material list
    ///
    /// Every difficult thing about a real block is **already in here**. Two materials meeting are
    /// a face conductance that is neither one's; a void is a face conductance of zero and a
    /// mobility of zero; a coating is a thin row of cells with its own `k`. A port that read the
    /// substances and rebuilt the operator would be a second implementation of `resolve`, and the
    /// first defect it had would be a physics difference nobody could see.
    ///
    /// It also fixes what a uniform-coefficient kernel cannot. `f·(Σ T_n − 6T)` has no per-cell
    /// `f`, so it conserves `Σ Cᵢ Tᵢ` only when every capacity is the same — which is why the GPU
    /// port could do one material and nothing else until it took these.
    pub fn coefficients(&self) -> Coefficients<'_> {
        Coefficients {
            kx: &self.kx,
            ky: &self.ky,
            kz: &self.kz,
            mobility: &self.mobility,
            source: &self.source,
        }
    }

    /// What an accelerator cannot run about this block, and why — empty when it can.
    ///
    /// A device is asked for **in the scene**, and a request that cannot be honoured is an error
    /// rather than a quiet fall back to the CPU: a run that silently changed where it ran is a run
    /// whose answer nobody asked for. This is the list of reasons, so the refusal can name them.
    pub fn unsupported_on_a_device(&self) -> Vec<&'static str> {
        let mut why = Vec::new();
        if !self.exposed.is_empty() {
            why.push("a surface film: `losing_from` is a boundary flux with no device pass yet");
        }
        if self.air.is_some() {
            why.push("air in the void: a per-cell boundary flux with no device pass yet");
        }
        if !self.gaps.is_empty() {
            why.push("radiation across a gap: a pair exchange, which is not a stencil");
        }
        if self.latent_anywhere {
            why.push("phase change: a rise is converted through an enthalpy, cell by cell");
        }
        why
    }

    /// The conductance of the face between two cells, or `None` if they are not face neighbours.
    ///
    /// `k_face · A / dx`, which for cubic cells is `k_face · dx`. This is the number the sweep
    /// actually uses, so a caller checking a joint against `Σ Lᵢ/(kᵢA)` by hand is checking the
    /// same arithmetic the march does rather than a restatement of it.
    pub fn face_conductance(
        &self,
        a: (usize, usize, usize),
        b: (usize, usize, usize),
    ) -> Option<Conductance> {
        let (nx, ny, _) = self.counts;
        self.index(a.0, a.1, a.2)?;
        self.index(b.0, b.1, b.2)?;
        let step = |p: usize, q: usize| (p as isize - q as isize).abs();
        let (di, dj, dk) = (step(a.0, b.0), step(a.1, b.1), step(a.2, b.2));
        let k = match (di, dj, dk) {
            (1, 0, 0) => self.kx[a.0.max(b.0) + (nx + 1) * (a.1 + ny * a.2)],
            (0, 1, 0) => self.ky[a.0 + nx * (a.1.max(b.1) + (ny + 1) * a.2)],
            (0, 0, 1) => self.kz[a.0 + nx * (a.1 + ny * a.2.max(b.2))],
            _ => return None,
        };
        Some(Conductance::from_si(k * self.dx.to_si()))
    }

    /// Rebuild everything derived from the materials: face conductivities, mobilities, the limit.
    ///
    /// Called by the constructor and by every [`fill`](Solid3D::fill), on the rule this workspace
    /// learned from `Puck::repack` — a mutator that leaves a cached solve stale reports a number
    /// that was true about the previous object, and it looks exactly like a number.
    fn resolve(&mut self) {
        self.resolve_rephasing(&[]);
    }

    /// [`resolve`](Solid3D::resolve), deriving the phase of each cell `rephase` marks from its
    /// temperature rather than carrying it. [`fill`](Solid3D::fill) marks the cells it put a
    /// different substance in; nothing else marks any.
    fn resolve_rephasing(&mut self, rephase: &[bool]) {
        self.count_faces();
        let (nx, ny, nz) = self.counts;
        let dx = self.dx.to_si();
        let volume = Volume::from_si(dx * dx * dx);

        // Per material, the **solid** triple and the **liquid** triple. A material whose
        // `FusionProps` names no liquid uses the solid triple for both, which is the one-phase
        // model. For its conductivity and capacity that is exact whenever the liquid sits at the
        // melting point, because a face with no temperature difference across it carries no heat
        // whatever its conductivity. For its emissivity it is not: a melted surface at the melting
        // point sends `εσ(T⁴ − T∞⁴)` whatever is behind it, so the one-phase model is a statement
        // that the finish survives melting, and a substance whose liquid radiates differently has
        // to say so.
        //
        // `NaN` conductivity and capacity for a substance that does not say, which `step` refuses
        // on rather than stepping with a plausible default; an emissivity of zero, which is what
        // every reader of a surface took one with no thermal properties to have.
        //
        // **The liquid's emissivity, not only its conductivity and capacity.** A finish is a
        // surface, and a surface that has melted is the liquid's: a black coat on a solid does not
        // survive it. This read the solid's whatever had melted, and so did the film and the
        // clearances — measured, a cell of liquid radiating at 900 K to 300 K, its solid 0.9 and its
        // liquid 0.1, lost 8.65 times what its own surface sends through its film, and two such
        // cells across a void exchanged 15.5 times their own pair's heat.
        let props: Vec<(Phase, Phase)> = self
            .materials
            .iter()
            .map(|s| {
                let solid = (
                    s.thermal.map_or(f64::NAN, |t| t.conductivity.to_si()),
                    s.heat_capacity(volume).map_or(f64::NAN, |c| c.to_si()),
                    s.thermal.map_or(0.0, |t| t.emissivity),
                );
                let liquid = match s.fusion.and_then(|f| f.liquid) {
                    Some(t) => (
                        t.conductivity.to_si(),
                        s.mass_of(volume).to_si() * t.specific_heat.to_si(),
                        t.emissivity,
                    ),
                    None => solid,
                };
                (solid, liquid)
            })
            .collect();

        // A liquid phase that differs from the solid is what makes the operator move with the front.
        // A material whose `liquid` is absent, or identical to its solid, needs no per-step rebuild.
        // An emissivity that differs is enough on its own: it moves the film's conductance and a
        // clearance's, which are as much the operator as the faces are.
        self.two_phase =
            self.materials
                .iter()
                .any(|s| match (s.fusion.and_then(|f| f.liquid), s.thermal) {
                    (Some(l), Some(t)) => {
                        l.conductivity != t.conductivity
                            || l.specific_heat != t.specific_heat
                            || l.emissivity != t.emissivity
                    }
                    _ => false,
                });
        // Where a step can carry a cell across its melting point further than the capacity it has
        // bounds — see `charge_crossing` — from the triples and the latent heat the cells are given.
        self.freeze_below = f64::NEG_INFINITY;
        self.melt_above = f64::INFINITY;
        for (s, ((_, solid, _), (_, liquid, _))) in self.materials.iter().zip(props.iter().copied())
        {
            let (Some(f), Some(_)) = (s.fusion, s.thermal) else {
                continue;
            };
            let latent = s.latent_energy(volume).map_or(0.0, |e| e.to_si());
            if latent > 0.0 {
                let point = f.melting_point.to_si();
                if liquid > solid {
                    self.freeze_below = self.freeze_below.max(point - latent / (liquid - solid));
                } else if solid > liquid {
                    self.melt_above = self.melt_above.min(point + latent / (solid - liquid));
                }
            }
        }

        // The phase state first, because everything below depends on it.
        //
        // Rebuilt from temperature only when there is not one yet, or for a cell `fill` has just
        // put new matter in. Once the sweep is carrying a fraction, re-deriving it would throw away
        // every partially melted cell — and this runs every step for a two-phase block, because a
        // cell's conductivity moves with its fraction.
        let fresh = self.melted.len() != self.cells.len();
        let carried = std::mem::take(&mut self.melted);
        self.melt_point = Vec::with_capacity(self.cells.len());
        self.latent = Vec::with_capacity(self.cells.len());
        self.melted = Vec::with_capacity(self.cells.len());
        for (c, prior) in carried
            .iter()
            .copied()
            .chain(std::iter::repeat(0.0))
            .take(self.cells.len())
            .enumerate()
        {
            let s = &self.materials[self.which[c] as usize];
            let (point, latent) = match (s.fusion, s.thermal) {
                (Some(f), Some(_)) => (
                    f.melting_point.to_si(),
                    s.latent_energy(volume).map_or(0.0, |e| e.to_si()),
                ),
                // A substance with fusion but no thermal properties cannot be stepped at all: its
                // capacity is `NaN` and the sweep refuses rather than guessing.
                _ => (f64::INFINITY, 0.0),
            };
            self.melt_point.push(point);
            self.latent.push(latent);
            self.melted.push(if fresh {
                phase_named_by(self.cells[c], point)
            } else {
                prior
            });
        }
        // After the loop and over only what is marked, rather than a test inside it: the loop runs
        // every step for a two-phase block, and only `fill` marks anything.
        for (c, _) in rephase.iter().enumerate().filter(|(_, marked)| **marked) {
            self.melted[c] = phase_named_by(self.cells[c], self.melt_point[c]);
        }

        // **Any cell with latent heat at all**, which is not the same question as `two_phase`.
        // `two_phase` asks whether the *operator* moves with the front; this asks whether
        // `add_kelvin` has to convert a rise through an enthalpy. A material that melts at
        // constant properties answers no to the first and yes to the second — and the fast sweep
        // path, gated on the wrong one, silently skipped the latent heat and let a freezing block
        // cool as if it were not freezing.
        //
        // **After the loop that fills `latent`, not before it.** The first attempt set this thirty
        // lines earlier, where the array is still the previous resolve's — or, on the first call,
        // empty. It read `false` every time and the same three tests failed the same way, which is
        // the useful part: the fix looked right and the placement was the bug.
        self.latent_anywhere = self.latent.iter().any(|l| *l > 0.0);
        // Only on a fresh build. Doing it every time would overwrite a checkpoint, and `resolve` is
        // called from the sweep now.
        if fresh {
            self.saved_melted.clone_from(&self.melted);
        }

        // **Mixed by melt fraction, arithmetically, and that choice is the one to argue about.**
        //
        // A half-melted cell is half of each, and across a cell the two phases sit side by side
        // rather than in series — heat crossing the mush passes through both at once, so a parallel
        // mean is what matches. The *face* between two cells stays the harmonic mean of whatever the
        // two cells came out as, because a face **is** series. Mixing one way and joining the other
        // is not an inconsistency; they are different geometries.
        //
        // It applies only inside the mush, one or two cells wide, and that is where the scheme's
        // remaining first-order error lives.
        //
        // The emissivity is mixed the same way and for the same reason: a surface half covered in
        // liquid is two surfaces side by side, and what it radiates is the sum of what each part
        // does, an area-weighted mean. Across a clearance that is exact rather than near: the two
        // patches' radiosity branches, `1/(ε_i A_i)` each, are in parallel and make one surface of
        // the mean, so the mean goes into the parallel-plate series rather than the two pairs'
        // coefficients being averaged — which for 0.9 and 0.1 half and half, facing 0.9, is 3.2%
        // low.
        //
        // **Where the front lies along the exposed face it is wrong, and by as much as the
        // conductivity is.** A surface melted or frozen by what it radiates is wholly one phase at
        // the face, and the mean says it is part of each: a liquid cell of 0.1 freezing under a
        // solid of 0.9 by its own radiation takes `ln(9)/0.8 · 0.9` = 2.47 times as long as a face
        // that turns solid at once. The time is one cell's, so the error falls with the cell, first
        // order — the same as the mixed conductivity's in that geometry, where it should be in
        // series.
        let mixed: Vec<Phase> = (0..self.cells.len())
            .map(|c| {
                // Nothing conducts nothing and holds nothing. Zeroing here rather than at every
                // reader is what makes void fall out of the rest of the arithmetic: the face
                // mean is already the harmonic one and guarded at zero, so a face touching void
                // carries zero without knowing what void is.
                if self.void[c] {
                    return (0.0, 0.0, 0.0);
                }
                mix(props[self.which[c] as usize], self.melted[c])
            })
            .collect();

        let k_of = |cell: usize| mixed[cell].0;

        // Read only where one was placed, so a block with no joints does no work and needs no
        // array. `joined` is what fills these, and it fills all three or none.
        let joined = !self.hx.is_empty();

        self.kx = vec![0.0; (nx + 1) * ny * nz];
        self.ky = vec![0.0; nx * (ny + 1) * nz];
        self.kz = vec![0.0; nx * ny * (nz + 1)];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let c = i + nx * (j + ny * k);
                    if i > 0 {
                        let f = i + (nx + 1) * (j + ny * k);
                        let k_f = series(k_of(c - 1), k_of(c));
                        self.kx[f] = if joined {
                            contact(k_f, self.hx[f], dx)
                        } else {
                            k_f
                        };
                    }
                    if j > 0 {
                        let f = i + nx * (j + (ny + 1) * k);
                        let k_f = series(k_of(c - nx), k_of(c));
                        self.ky[f] = if joined {
                            contact(k_f, self.hy[f], dx)
                        } else {
                            k_f
                        };
                    }
                    if k > 0 {
                        let f = i + nx * (j + ny * k);
                        let k_f = series(k_of(c - nx * ny), k_of(c));
                        self.kz[f] = if joined {
                            contact(k_f, self.hz[f], dx)
                        } else {
                            k_f
                        };
                    }
                }
            }
        }

        self.cap_solid = (0..self.cells.len())
            .map(|c| props[self.which[c] as usize].0 .1)
            .collect();
        self.cap_liquid = (0..self.cells.len())
            .map(|c| props[self.which[c] as usize].1 .1)
            .collect();
        self.capacity = vec![0.0; self.cells.len()];
        self.mobility = vec![0.0; self.cells.len()];
        self.face_sum = vec![0.0; self.cells.len()];
        self.worst_rate = 0.0;
        let mut nowhere = false;
        // The cells' range, for a block whose bounds anything can reach, and for no other: a block
        // of wax, whose bound is −630 K, stepped 8% slower with this taken on every cell.
        let (freezes, melts) = self.reachable();
        let ranged = freezes || melts;
        let mut range = (f64::INFINITY, f64::NEG_INFINITY);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let c = i + nx * (j + ny * k);
                    if ranged && !self.void[c] {
                        range = widened(range, self.cells[c]);
                    }
                    let cap = mixed[c].1;
                    self.capacity[c] = cap;
                    // A void cell has no capacity, and `dx/0` is an infinity that would poison
                    // the limit and the sweep alike. Zero mobility is the honest value: nothing
                    // moves in a cell that holds nothing.
                    self.mobility[c] = if self.void[c] { 0.0 } else { dx / cap };
                    let sum = self.kx[i + (nx + 1) * (j + ny * k)]
                        + self.kx[i + 1 + (nx + 1) * (j + ny * k)]
                        + self.ky[i + nx * (j + (ny + 1) * k)]
                        + self.ky[i + nx * (j + 1 + (ny + 1) * k)]
                        + self.kz[i + nx * (j + ny * k)]
                        + self.kz[i + nx * (j + ny * (k + 1))];
                    // NaN loses every comparison, so `max` would step straight past a substance that
                    // does not say what it conducts. Tested for, and the sweep is refused.
                    //
                    // An exposed cell also loses to air, and that conductance leaves the cell
                    // exactly as a face conductance does — so it belongs in the same sum. It is
                    // divided by `dx` because `mobility` is `dx/capacity`: the face terms are
                    // conductivities and this one is already a conductance.
                    self.face_sum[c] = sum;
                    let rate = self.mobility[c] * sum;
                    if rate.is_nan() {
                        nowhere = true;
                    } else {
                        self.worst_rate = self.worst_rate.max(rate);
                    }
                }
            }
        }
        if nowhere {
            self.worst_rate = f64::NAN;
        }
        self.cells_range = ranged.then_some(range);
        self.worst_stale = false;
        self.phase_props = props;
        // Kept rather than dropped, for the film and the clearances. Moving it costs no time —
        // measured, a step of a block that keeps it and reads nothing from it is the step of one
        // that drops it — and it holds twenty-four bytes a cell. Before `find_gaps`, which reads
        // its emissivities.
        self.derived = mixed;
        self.find_gaps();
        self.find_facing_air();
    }

    /// One cell's conductivity, capacity and emissivity, as `resolve` mixes them by its melted
    /// fraction — see the comment there — and zero for void.
    fn mixed_at(&self, c: usize) -> Phase {
        if self.void[c] {
            return (0.0, 0.0, 0.0);
        }
        mix(self.phase_props[self.which[c] as usize], self.melted[c])
    }

    /// The conductivity of the face between two cells of conductivities `low` and `high`, low
    /// index first, as `resolve` joins them: their series mean, with the contact stated on that
    /// face, if any, in series with it. `contacts` is the axis's `hx`, `hy` or `hz`, empty for a
    /// block nobody has joined, and `face` the face's index into it.
    fn face_conductivity(&self, contacts: &[f64], face: usize, low: f64, high: f64) -> f64 {
        let k_f = series(low, high);
        if contacts.is_empty() {
            k_f
        } else {
            contact(k_f, contacts[face], self.dx.to_si())
        }
    }

    /// Give a cell the capacity `cap`, and the mobility `resolve` gives it: `dx / cap`, and zero
    /// for void.
    fn hold(&mut self, c: usize, cap: f64) {
        self.capacity[c] = cap;
        self.mobility[c] = if self.void[c] {
            0.0
        } else {
            self.dx.to_si() / cap
        };
    }

    /// A cell's conduction sum `Σ_f k_f`, over its six faces in the order `resolve` adds them.
    fn face_sum_at(&self, i: usize, j: usize, k: usize) -> f64 {
        let (nx, ny, _) = self.counts;
        self.kx[i + (nx + 1) * (j + ny * k)]
            + self.kx[i + 1 + (nx + 1) * (j + ny * k)]
            + self.ky[i + nx * (j + (ny + 1) * k)]
            + self.ky[i + nx * (j + 1 + (ny + 1) * k)]
            + self.kz[i + nx * (j + ny * k)]
            + self.kz[i + nx * (j + ny * (k + 1))]
    }

    /// One cell's conduction rate, `mobility · Σ_f k_f`, in reciprocal seconds.
    fn rate_at(&self, c: usize) -> f64 {
        self.mobility[c] * self.face_sum[c]
    }

    /// `max_i mobility_i · Σ_f k_f` over every cell, or `NaN` if any cell's is not a number —
    /// what `resolve` finds, found again.
    fn scan_worst(&self) -> f64 {
        (0..self.cells.len()).fold(0.0, |worst, c| worse(worst, self.rate_at(c)))
    }

    /// The conduction part of the stability rate: what `resolve` found, or, after a mutator moved
    /// a phase since, the maximum found again from the per-cell terms it left right.
    fn conduction_rate(&self) -> f64 {
        if self.worst_stale {
            self.scan_worst()
        } else {
            self.worst_rate
        }
    }

    /// After a mutator between steps moved cell `c`'s melted fraction from `before`: re-derive
    /// what the fraction decides, which is something only when a substance here has a liquid
    /// phase that differs from its solid.
    fn rephased(&mut self, c: usize, before: f64) {
        if self.two_phase && self.melted[c] != before {
            self.rederive_cell(c);
        }
    }

    /// Re-derive what one cell's phase decides: its conductivity, capacity and emissivity, the
    /// faces it shares, the conduction sums of the cells across them, and the clearances it is an
    /// end of — the only entries one cell's phase reaches. The maximum over cells is left for
    /// [`conduction_rate`](Solid3D::conduction_rate) to find again; see `worst_stale` for why.
    ///
    /// **What `resolve` would make of the same state, to the bit**, and held to that by
    /// `a_phase_set_between_steps.rs`, which compares every face, capacity, mobility and limit
    /// with a full `resolve` after each of twelve scenarios of mutators — one substance and two,
    /// a void, contacts, a film, a clearance whose ends melt, and one melt beside the middle of a
    /// cube for each direction of the neighbours' sums — and the state a step on each leaves.
    /// `series`, `contact`, `mix` and `exchange` are shared; the
    /// rest is the same expressions written twice, because `resolve` runs every step of a
    /// two-phase block, and calling these helpers from its loops made that step 4% to 11% slower
    /// from 16³ to 64³, pinned to one core.
    fn rederive_cell(&mut self, c: usize) {
        let (nx, ny, nz) = self.counts;
        let (i, j, k) = (c % nx, (c / nx) % ny, c / (nx * ny));
        let here = self.mixed_at(c);
        let was = std::mem::replace(&mut self.derived[c], here).2;
        if here.2.to_bits() != was.to_bits() {
            self.regap(c);
        }
        if i > 0 {
            let f = i + (nx + 1) * (j + ny * k);
            self.kx[f] = self.face_conductivity(&self.hx, f, self.mixed_at(c - 1).0, here.0);
        }
        if i + 1 < nx {
            let f = i + 1 + (nx + 1) * (j + ny * k);
            self.kx[f] = self.face_conductivity(&self.hx, f, here.0, self.mixed_at(c + 1).0);
        }
        if j > 0 {
            let f = i + nx * (j + (ny + 1) * k);
            self.ky[f] = self.face_conductivity(&self.hy, f, self.mixed_at(c - nx).0, here.0);
        }
        if j + 1 < ny {
            let f = i + nx * (j + 1 + (ny + 1) * k);
            self.ky[f] = self.face_conductivity(&self.hy, f, here.0, self.mixed_at(c + nx).0);
        }
        if k > 0 {
            let f = i + nx * (j + ny * k);
            self.kz[f] = self.face_conductivity(&self.hz, f, self.mixed_at(c - nx * ny).0, here.0);
        }
        if k + 1 < nz {
            let f = i + nx * (j + ny * (k + 1));
            self.kz[f] = self.face_conductivity(&self.hz, f, here.0, self.mixed_at(c + nx * ny).0);
        }
        self.hold(c, here.1);
        self.face_sum[c] = self.face_sum_at(i, j, k);
        if i > 0 {
            self.face_sum[c - 1] = self.face_sum_at(i - 1, j, k);
        }
        if i + 1 < nx {
            self.face_sum[c + 1] = self.face_sum_at(i + 1, j, k);
        }
        if j > 0 {
            self.face_sum[c - nx] = self.face_sum_at(i, j - 1, k);
        }
        if j + 1 < ny {
            self.face_sum[c + nx] = self.face_sum_at(i, j + 1, k);
        }
        if k > 0 {
            self.face_sum[c - nx * ny] = self.face_sum_at(i, j, k - 1);
        }
        if k + 1 < nz {
            self.face_sum[c + nx * ny] = self.face_sum_at(i, j, k + 1);
        }
        self.worst_stale = true;
    }

    /// After cell `c`'s emissivity moved, put right the clearances it is one end of, as
    /// `find_gaps` would make them — the same coefficients, to the bit, in the same order.
    ///
    /// Its partners are found by walking: along each of the six directions, the first solid cell
    /// past one void cell or more, which is how `find_gaps` pairs a grid line. A cell with no void
    /// beside it is the end of no pair and costs nothing, and one that is costs a pass over the
    /// pairs for each partner — not over the block, so a region set one cell at a time still costs
    /// what its cells do. Where a pair appears or vanishes, because an end stops or starts
    /// radiating, its place in the list is what `find_gaps` decides, so that is called instead.
    fn regap(&mut self, c: usize) {
        let (nx, ny, nz) = self.counts;
        let (i, j, k) = (c % nx, (c / nx) % ny, c / (nx * ny));
        let now = self.derived[c].2;
        let area = self.dx.to_si() * self.dx.to_si();
        // For each axis, `c`'s position along it, the cells on it and its stride.
        let mut partners = Vec::new();
        for (at, n, stride) in [(i, nx, 1), (j, ny, nx), (k, nz, nx * ny)] {
            // Backwards, then forwards: the first non-void cell, and whether any void came first.
            let mut p = c;
            let mut m = at;
            while m > 0 {
                m -= 1;
                p -= stride;
                if !self.void[p] {
                    if p + stride != c {
                        partners.push((p, c));
                    }
                    break;
                }
            }
            let (mut p, mut m) = (c, at);
            while m + 1 < n {
                m += 1;
                p += stride;
                if !self.void[p] {
                    if p != c + stride {
                        partners.push((c, p));
                    }
                    break;
                }
            }
        }
        for (a, b) in partners {
            let other = self.derived[if a == c { b } else { a }].2;
            let listed = self.gaps.iter().position(|g| (g.0, g.1) == (a, b));
            match (listed, exchange(now, other, area)) {
                (Some(n), Some(coefficient)) => self.gaps[n].2 = coefficient,
                (None, None) => {}
                // Listed and no longer radiating, or radiating and not listed.
                _ => {
                    self.find_gaps();
                    return;
                }
            }
        }
    }

    /// Which solid cells have a face on the void, when there is air in it to lose to.
    ///
    /// A neighbour that is void, counted the way `air_conductance` counts it — six axis
    /// neighbours, none past the block's own outer faces — so every solid cell whose air
    /// conductance is not zero is listed. A listed one can still have none: a cell that conducts
    /// nothing has no half cell for the air to reach it through, and `charge` passes over a loss
    /// of zero.
    fn find_facing_air(&mut self) {
        self.facing_air.clear();
        if self.air.is_none() {
            return;
        }
        let (nx, ny, nz) = self.counts;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let c = i + nx * (j + ny * k);
                    if self.void[c] {
                        continue;
                    }
                    let faces = (i > 0 && self.void[c - 1])
                        || (i + 1 < nx && self.void[c + 1])
                        || (j > 0 && self.void[c - nx])
                        || (j + 1 < ny && self.void[c + nx])
                        || (k > 0 && self.void[c - nx * ny])
                        || (k + 1 < nz && self.void[c + nx * ny]);
                    if faces {
                        self.facing_air.push((i, j, k));
                    }
                }
            }
        }
    }

    /// Which solid cells face each other across void, and how strongly they radiate.
    ///
    /// Along each grid line, a run of void with solid at both ends is a gap, and the two cells
    /// at its ends see each other. Diagonals are not paired: a face is what radiates, and the
    /// grid's faces are axis-aligned.
    ///
    /// Each end radiates with the emissivity of its phase as the operator was derived, from
    /// `derived`, so a melted end is its liquid's. It read the substance's solid whatever had
    /// melted.
    fn find_gaps(&mut self) {
        if !self.void.iter().any(|v| *v) {
            self.gaps.clear();
            self.group_gap_ends();
            return;
        }
        let (nx, ny, nz) = self.counts;
        let dx = self.dx.to_si();
        let area = dx * dx;
        let emissivity = |c: usize| self.derived[c].2;

        // Written over the pairs found last time, in place, noting whether any is not the pair it
        // replaces — see the end, which is what this is for.
        let mut written = 0usize;
        let mut moved = false;
        // One pass per axis. `line` walks the indices of a single row, column or pillar.
        let mut walk = |line: Vec<usize>| {
            let mut last_solid: Option<usize> = None;
            let mut void_between = false;
            for c in line {
                if self.void[c] {
                    void_between = last_solid.is_some();
                    continue;
                }
                if let (Some(a), true) = (last_solid, void_between) {
                    if let Some(coefficient) = exchange(emissivity(a), emissivity(c), area) {
                        match self.gaps.get_mut(written) {
                            Some(old) => {
                                moved |= (old.0, old.1) != (a, c);
                                *old = (a, c, coefficient);
                            }
                            None => {
                                moved = true;
                                self.gaps.push((a, c, coefficient));
                            }
                        }
                        written += 1;
                    }
                }
                last_solid = Some(c);
                void_between = false;
            }
        };

        for k in 0..nz {
            for j in 0..ny {
                walk((0..nx).map(|i| i + nx * (j + ny * k)).collect());
            }
        }
        for k in 0..nz {
            for i in 0..nx {
                walk((0..ny).map(|j| i + nx * (j + ny * k)).collect());
            }
        }
        for j in 0..ny {
            for i in 0..nx {
                walk((0..nz).map(|k| i + nx * (j + ny * k)).collect());
            }
        }
        if written != self.gaps.len() {
            moved = true;
            self.gaps.truncate(written);
        }
        // **Grouped again only when these are not the pairs it was grouped from.** The grouping
        // is a statement about which pairs there are and in what order, not about how strongly
        // each radiates, and a two-phase block finds its clearances again after every step and
        // nearly always finds the same ones.
        if moved {
            self.group_gap_ends();
        }
    }

    /// Fill `gap_ends` and `gap_end_pairs` from `gaps`: each end cell, in index order, and the
    /// pairs it is in, in the order `gaps` lists them — so a row's sum is added in one order, the
    /// same on every run and every thread count.
    fn group_gap_ends(&mut self) {
        let mut ends: Vec<(usize, usize)> = self
            .gaps
            .iter()
            .enumerate()
            .flat_map(|(n, (a, b, _))| [(*a, n), (*b, n)])
            .collect();
        ends.sort_unstable();
        self.gap_ends.clear();
        self.gap_end_pairs.clear();
        let (nx, ny, nz) = self.counts;
        for (cell, n) in ends {
            match self.gap_ends.last_mut() {
                Some(end) if end.cell == cell => end.count += 1,
                _ => {
                    let (i, j, k) = (cell % nx, (cell / nx) % ny, cell / (nx * ny));
                    self.gap_ends.push(GapEnd {
                        cell,
                        at: (i, j, k),
                        outer: i == 0
                            || i + 1 == nx
                            || j == 0
                            || j + 1 == ny
                            || k == 0
                            || k + 1 == nz,
                        first: self.gap_end_pairs.len(),
                        count: 1,
                    });
                }
            }
            self.gap_end_pairs.push(n);
        }
    }

    /// The clearances in this block, as sheets of facing area rather than as cell pairs.
    ///
    /// A pair is what the exchange is computed on; a patch is what a **view factor** is a
    /// statement about, and the two are not the same object. Grouped by axis and separation, then
    /// by connectivity across the gap, so two unrelated clearances at the same width do not
    /// average into one patch that describes neither.
    ///
    /// Empty for a block with no void, which is every block that existed before there was one.
    pub fn gap_patches(&self) -> Vec<GapPatch> {
        let (nx, ny, _) = self.counts;
        let dx = self.dx.to_si();
        // A pair's axis and separation are recoverable from the two indices: the stride between
        // them says which way it ran, and the multiple of that stride says how far.
        let strides = [(0usize, 1usize), (1, nx), (2, nx * ny)];
        let mut sheets: BTreeMap<(usize, usize), Vec<(usize, usize)>> = BTreeMap::new();
        for (a, b, _) in &self.gaps {
            let delta = b - a;
            let Some(&(axis, stride)) = strides.iter().find(|(ax, st)| {
                delta % st == 0
                    && delta / st > 1
                    && match ax {
                        0 => delta < nx,
                        1 => delta < nx * ny,
                        _ => true,
                    }
            }) else {
                continue;
            };
            // The **void** between them, not the distance between their centres. A view factor is
            // a statement about two surfaces, and the surfaces are the faces bounding the gap: a
            // one-cell clearance puts them one cell apart while their centres are two.
            let cells = delta / stride - 1;
            let (i, j, k) = (a % nx, (a / nx) % ny, a / (nx * ny));
            // The two coordinates *across* the gap, which is what a patch is measured in.
            let across = match axis {
                0 => (j, k),
                1 => (i, k),
                _ => (i, j),
            };
            sheets.entry((axis, cells)).or_default().push(across);
        }

        let mut out = Vec::new();
        for ((_, cells), mut face) in sheets {
            face.sort_unstable();
            let mut seen = vec![false; face.len()];
            let index: BTreeMap<(usize, usize), usize> =
                face.iter().enumerate().map(|(n, p)| (*p, n)).collect();
            for start in 0..face.len() {
                if seen[start] {
                    continue;
                }
                // Four-connected flood fill over the facing cells, so a patch is a sheet somebody
                // could point at rather than every pair that happens to share a width.
                let mut stack = vec![start];
                seen[start] = true;
                let mut members = Vec::new();
                while let Some(n) = stack.pop() {
                    let (u, v) = face[n];
                    members.push((u, v));
                    for (du, dv) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
                        let (nu, nv) = (u as i64 + du, v as i64 + dv);
                        if nu < 0 || nv < 0 {
                            continue;
                        }
                        if let Some(&m) = index.get(&(nu as usize, nv as usize)) {
                            if !seen[m] {
                                seen[m] = true;
                                stack.push(m);
                            }
                        }
                    }
                }
                let (u0, u1) = (
                    members.iter().map(|m| m.0).min().unwrap_or(0),
                    members.iter().map(|m| m.0).max().unwrap_or(0),
                );
                let (v0, v1) = (
                    members.iter().map(|m| m.1).min().unwrap_or(0),
                    members.iter().map(|m| m.1).max().unwrap_or(0),
                );
                let (wu, wv) = (u1 + 1 - u0, v1 + 1 - v0);
                out.push(GapPatch {
                    pairs: members.len(),
                    span: (wu as f64 * dx, wv as f64 * dx),
                    distance: cells as f64 * dx,
                    rectangular: members.len() == wu * wv,
                });
            }
        }
        out
    }

    /// Mark cells as **nothing** — not a substance, not part of the block.
    ///
    /// A grid had no void until now, so the cells a part did not occupy were some other material,
    /// and an assembly of two parts in air was two parts buried in whatever the block was made
    /// of. `ARCHITECTURE.md` names it: "insulating it is a substance with a low conductivity,
    /// which is not the same thing". A low conductivity still conducts, still stores heat, and
    /// still sets a stability limit; nothing does none of those.
    ///
    /// A void cell holds no heat, conducts to nothing — every face it touches carries zero, which
    /// the harmonic mean already gives for a zero conductivity — takes no share of what arrives
    /// on the bus, and is left out of every average. Its temperature is **not a number**, because
    /// there is nothing there to have one, and `temperature_at` says so rather than returning a
    /// zero somebody would plot.
    ///
    /// A cell emptied after a source was spread over it **takes its share with it**: the step
    /// generates nothing there, [`generated_power`](Solid3D::generated_power) stops counting it,
    /// and the rest of the spread is not raised to make up the difference, because the block does
    /// not remember which box the watts were stated for. The share stays with the place, as the
    /// cell's temperature does, so [`fill`](Solid3D::fill)ing the cell again generates it again —
    /// unless a later [`dissipating`](Solid3D::dissipating) selected the cell, which replaces its
    /// share as it replaces a solid cell's, with nothing while the cell is void.
    ///
    /// # What crosses a gap, and what does not
    ///
    /// **Radiation does.** Two solid cells facing each other along a grid line across a run of
    /// void exchange `σA(T₁⁴ − T₂⁴)/(1/ε₁ + 1/ε₂ − 1)`, the parallel-plate series, with an
    /// exchange factor of **one**, at the temperatures of their two facing surfaces: each half a
    /// cell from its centre, behind its half cell, and beside the clearance's air if
    /// [`air_in`](Solid3D::air_in) put any there. Each `ε` is its cell's phase's: a cell that has
    /// melted radiates as its substance's
    /// [liquid](pantometry_core::substance::FusionProps::liquid) does, and one
    /// in the mush as the area-weighted mean of the two.
    ///
    /// That factor was written down here as a known approximation — a wide gap has a view factor
    /// well under one, so charging it as one was said to couple a wide gap too hard. Measuring it
    /// says otherwise, and the correction is worth more than the caveat was. The sides of a gap in
    /// this model are the **block's own outer faces**, and an insulated boundary is implemented
    /// as a mirror; a mirror puts an image of each surface beyond it and the images tile the
    /// plane, so the pair *is* two infinite parallel plates and one is exact at every width. What
    /// the old note described was a different geometry from the one the model has.
    ///
    /// The geometry it described is a real one, though — two parts floating in vacuum, open to
    /// space, where most of what leaves one surface does miss the other. Which of the two a scene
    /// means is a statement about its boundary that a grid cannot infer, so it is not chosen here.
    /// It is *reported*: [`gap_patches`](Solid3D::gap_patches) groups a clearance into the sheets
    /// a view factor is a statement about and [`GapPatch::view_factor`] gives the open-gap number,
    /// so the difference between the two readings is a factor somebody can see. For a 32 mm square
    /// part 16 mm under a lid it is **2.4x**.
    ///
    /// Still not here: a radiative-exchange solver. Side walls that are real material at a real
    /// temperature do not radiate into the gap at all, and that is a domain rather than a boundary
    /// condition.
    ///
    /// **Convection does not.** A gap full of air carries heat by moving that air, which needs a
    /// Rayleigh number and a correlation, and a correlation is not a closed form. So a gap in
    /// air is coupled *less* here than it really is, by however much the convection would have
    /// carried — and for a millimetre-scale gap at modest temperatures that is the same order as
    /// the radiation, so the answer is a lower bound rather than an estimate.
    ///
    /// **Conduction does not**, and that is right: there is nothing there to conduct through.
    ///
    /// It is still not a fluid. A part in air also loses heat to the room, and that is
    /// [`Solid3D::losing_from`], which is about the block's **outer** faces.
    pub fn empty(mut self, which: impl Fn(usize, usize, usize) -> bool) -> Solid3D {
        let (nx, ny, nz) = self.counts;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if which(i, j, k) {
                        self.void[i + nx * (j + ny * k)] = true;
                    }
                }
            }
        }
        self.resolve();
        self
    }

    /// Whether this cell is void — nothing rather than a substance.
    pub fn is_void(&self, i: usize, j: usize, k: usize) -> bool {
        self.index(i, j, k).map(|c| self.void[c]).unwrap_or(false)
    }

    /// How many cells hold nothing.
    pub fn void_cells(&self) -> usize {
        self.void.iter().filter(|v| **v).count()
    }

    /// A **contact conductance** on the interior faces `which` names, in W·m⁻²·K⁻¹.
    ///
    /// Called for every interior face as `which(axis, i, j, k)`, where the triple is the cell on
    /// the face's high side: `Axis::X` with `i` is the face between `(i-1, j, k)` and `(i, j, k)`.
    /// `Some(h)` places a contact there and `None` leaves the face as it was.
    ///
    /// # The resistance this adds, and why a face is the only place it can go
    ///
    /// A face already carries the two half cells in series, `dx/(2k_L) + dx/(2k_R)`, which is the
    /// harmonic mean written as a conductivity. A contact adds `1/h` to that resistance:
    ///
    /// ```text
    /// 1/k_face = 1/k_series + 1/(dx·h)
    /// ```
    ///
    /// so `h = ∞` returns the harmonic mean **exactly** — an unstated face is unchanged, not
    /// approximately unchanged — and `h = 0` carries nothing at all, which is a clearance.
    ///
    /// **This is the only way to state a layer thinner than a cell.** A bolted joint, a thermal
    /// interface material, a solder layer, an oxide: all of them are resistances of a few tens of
    /// microns in a part discretised at a millimetre. A layer of *material* is at least one cell
    /// thick, so the thinnest resistance a grid can state is `dx/k`, and a grid fine enough to
    /// resolve a real joint is a hundred times finer in every direction — a million times the work
    /// to carry a layer with no interesting field inside it. That is a **floor**, not a
    /// discretisation error: refining the mesh does not approach the right answer, it only makes
    /// the floor lower and the run longer.
    ///
    /// A contact has a resistance and no thickness, so there is no floor and no grid dependence.
    ///
    /// Measured on the power module that ships with this workspace: its 100 µm solder was written
    /// as a 1.5 mm region because 1.5 mm was the floor, carrying **fifteen times** its own
    /// resistance, and the interface material between the baseplate and the cold plate — the
    /// largest resistance in a real junction-to-ambient path — could not be written at all.
    ///
    /// # What it does not change
    ///
    /// **Stability.** A contact only ever lowers a face conductivity, so `Σ_f k_f` falls, so the
    /// explicit limit gets looser. Adding one can never make a step that was stable unstable.
    ///
    /// **Conservation.** The same conductance is read from both sides of the face, so what leaves
    /// one cell arrives in the other exactly, as it does without one.
    ///
    /// A face named more than once takes the **last** value, and a face on the outer boundary is
    /// ignored: there is no cell on the other side of it, and what a block loses at its boundary
    /// is [`losing_from`](Solid3D::losing_from).
    pub fn joined(mut self, which: impl Fn(Axis, usize, usize, usize) -> Option<f64>) -> Solid3D {
        let (nx, ny, nz) = self.counts;
        if self.hx.is_empty() {
            self.hx = vec![f64::INFINITY; (nx + 1) * ny * nz];
            self.hy = vec![f64::INFINITY; nx * (ny + 1) * nz];
            self.hz = vec![f64::INFINITY; nx * ny * (nz + 1)];
        }
        let mut bad = self.bad_contact;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let mut note = |h: f64| {
                        if (h.is_nan() || h < 0.0) && bad.is_none() {
                            bad = Some(h);
                        }
                        h
                    };
                    if i > 0 {
                        if let Some(h) = which(Axis::X, i, j, k) {
                            self.hx[i + (nx + 1) * (j + ny * k)] = note(h);
                        }
                    }
                    if j > 0 {
                        if let Some(h) = which(Axis::Y, i, j, k) {
                            self.hy[i + nx * (j + (ny + 1) * k)] = note(h);
                        }
                    }
                    if k > 0 {
                        if let Some(h) = which(Axis::Z, i, j, k) {
                            self.hz[i + nx * (j + ny * k)] = note(h);
                        }
                    }
                }
            }
        }
        self.bad_contact = bad;
        self.resolve();
        self
    }

    /// The contact conductance stated on one interior face, or `None` if none was.
    ///
    /// The counterpart of [`face_conductance`](Solid3D::face_conductance), which reports what the
    /// sweep uses — the two half cells **and** any contact, together. This reports the contact
    /// alone, so a caller can check the series by hand rather than infer it from the total.
    pub fn contact_at(&self, axis: Axis, i: usize, j: usize, k: usize) -> Option<f64> {
        let (nx, ny, nz) = self.counts;
        if self.hx.is_empty() || i >= nx || j >= ny || k >= nz {
            return None;
        }
        let h = match axis {
            Axis::X if i > 0 => self.hx[i + (nx + 1) * (j + ny * k)],
            Axis::Y if j > 0 => self.hy[i + nx * (j + (ny + 1) * k)],
            Axis::Z if k > 0 => self.hz[i + nx * (j + ny * k)],
            _ => return None,
        };
        h.is_finite().then_some(h)
    }

    /// Expose a face to an environment, so the block can lose heat through it.
    ///
    /// **Until this existed a `Solid3D` was adiabatic on all six faces**, which means no
    /// three-dimensional thermal scene could reach a steady state: every one of them warmed
    /// for as long as it ran. That is honest for a pulse and useless for the question a
    /// designer actually asks — *what temperature does this run at* — which is the question a
    /// chip package, a magnet busbar, a motor and a factory cell are all asking.
    ///
    /// The loss is [`Environment::loss_from`]'s, convective **and** radiative, applied per
    /// boundary cell at that cell's own temperature and the emissivity of its own phase. Per cell
    /// rather than to a mean, because that is what makes a gradient: the middle of a face runs
    /// hotter than its edge, and a lumped loss says it does not.
    ///
    /// The environment's `area` is the **whole face's** area; each boundary cell is charged its
    /// share, `area / cells_on_the_face`. Stating the face's area rather than a cell's keeps
    /// the number a caller writes independent of the grid they chose, which is what lets the
    /// same scene be refined without becoming a different problem — the property
    /// `pantometry-world`'s `verify` sweep depends on.
    ///
    /// Exposing a face **tightens the stability limit**, because a cell that can also lose to
    /// air has more conductance leaving it. [`Solid3D::max_stable_dt`] accounts for it; a
    /// caller stepping by hand past the returned limit is refused as ever.
    ///
    /// An environment no surface can have — a coefficient that is negative or not a number, an
    /// area that is negative or not finite, an ambient below absolute zero — is taken here and
    /// refused by [`step`](Domain::step), naming the face. A coefficient of `+∞` is a surface
    /// held at the ambient behind its half cell, and is not refused.
    pub fn losing_from(mut self, face: Face, environment: Environment) -> Solid3D {
        self.exposed.insert(face, environment);
        self.resolve();
        self
    }

    /// Expose a **box on a face** rather than the whole face, in that face's two axes.
    ///
    /// `within` is `[lo0, lo1, hi0, hi1]` over [`Face::axes`], half-open, so a single cell is
    /// `hi = lo + 1`. Everything outside it stays insulated.
    ///
    /// # Why a face entire was not enough
    ///
    /// **A part is bolted at pads.** A bracket carries heat from whatever is mounted on it down
    /// through two or three bolt bosses into a chassis, and that path — along the bracket to the
    /// pads — is the whole thermal question its shape answers. Cooled over its entire footprint
    /// there is no path: every cell sheds where it stands, the part falls together, and the field
    /// is a constant. Measured on `29-a-designed-bracket-becomes-cells`, which rasterises 4 100
    /// cells from an STL and holds a range of **0.067 K**.
    ///
    /// Stating a smaller `area` does not do it. The area is divided among the cells on the face,
    /// so a tenth of the area is a tenth of the conductance **spread over the whole face** — the
    /// right total and the wrong place, which is exactly the difference a shape is for.
    ///
    /// The environment's `area` is the **patch's** area under this, not the face's, for the same
    /// reason it is the face's under [`losing_from`](Solid3D::losing_from): a caller states what
    /// the part exposes and the grid divides it, so refining the mesh does not change the problem.
    pub fn losing_from_within(
        mut self,
        face: Face,
        within: [usize; 4],
        environment: Environment,
    ) -> Solid3D {
        self.patch.insert(face, within);
        self.exposed.insert(face, environment);
        self.resolve();
        self
    }

    /// The box on a face that is exposed, if it is not the whole face.
    pub fn patch_on(&self, face: Face) -> Option<[usize; 4]> {
        self.patch.get(&face).copied()
    }

    /// Fill this block's **void with air**, so the surfaces inside it convect.
    ///
    /// # A part in a housing could not lose heat to the air around it
    ///
    /// [`losing_from`](Solid3D::losing_from) works on the block's six **outer** faces. A part
    /// rasterised inside a block, or a bar with a clearance beside it, has surfaces that are
    /// interior to the grid — and until this existed nothing could reach them, so a part
    /// surrounded by void shed heat only by [radiating](Solid3D::empty) to whatever faced it
    /// across the gap. Measured: a copper bar voided in on both sides, with a film stated on the
    /// face its neighbours occupy, sat at its initial 200 °C for the whole run and the run
    /// completed and reported four figures.
    ///
    /// # Convection only, and the area is the grid's
    ///
    /// A transparent gas does not radiate. What a surface facing a clearance exchanges with the
    /// surface across it is the parallel-plate pairing `find_gaps` already computes, and that runs
    /// beside this rather than instead of it — both paths are real and they are in parallel.
    ///
    /// The area is `dx²` per touching face, counted from the grid. A caller states what an outer
    /// face exposes because a rasterised part covers less of it than the grid does; the faces
    /// touching an internal void are exactly the ones the grid has, so there is nothing to state
    /// and nothing to state wrongly. Each face's film is in series with the half cell behind it,
    /// as an outer face's is, so the air acts at the surface it touches rather than at the cell
    /// centre.
    ///
    /// # What it models and what it does not
    ///
    /// **Well-mixed air at a fixed temperature.** The gas has no state of its own: it does not
    /// warm, it does not move, and two cavities in one block share it. That is the same model
    /// `losing_from` uses for the air outside, and it is exact for a housing whose air is vented
    /// or whose walls are the heat sink; it is wrong for a sealed cavity small enough that the
    /// part heats its own air, which is a fluid problem and belongs to `pantometry-fluid`.
    ///
    /// **Absent is a vacuum**, which is what every void in this crate was, so no block changes
    /// unless it asks.
    ///
    /// Filling the void tightens the stability limit for every cell touching it, exactly as
    /// exposing a face does.
    ///
    /// A coefficient that is negative or not a number, or an ambient below absolute zero, is
    /// refused by [`step`](Domain::step). Both used to run: a negative or `NaN` coefficient lost
    /// nothing, and `steady_state` **panicked** on a `NaN` ambient. `+∞` holds the surfaces at the
    /// air's temperature and is not refused.
    pub fn air_in(mut self, ambient: Temperature, w_per_m2_k: f64) -> Solid3D {
        self.air = Some((ambient.to_si(), w_per_m2_k));
        self.resolve();
        self
    }

    /// The air in this block's void, as `(ambient, W·m⁻²·K⁻¹)`, if it has any.
    pub fn air(&self) -> Option<(Temperature, f64)> {
        self.air.map(|(t, h)| (Temperature::from_si(t), h))
    }

    /// The temperature field at which every cell balances, found without marching to it.
    ///
    /// # Why a march is not enough
    ///
    /// A design answer is a steady-state answer — a junction temperature, a margin, a rise above
    /// ambient — and the only way to get one here was to run for several time constants and check
    /// that it had stopped moving. That is expensive and it is easy to get wrong in the direction
    /// that looks like success: `29-a-designed-bracket-becomes-cells` needed **900 s** from a cold
    /// start, nine test binaries walk every shipped scene, and the app gate stopped fitting inside
    /// ten minutes. The run was shortened by starting it near the answer, which is a way of
    /// avoiding the question rather than answering it.
    ///
    /// # The same operator, not a second one
    ///
    /// The residual is `flux_at` — the function the sweep itself marches with — in watts, plus the
    /// source and the radiative pairs the sweep applies beside it. So this solves the balance the
    /// march converges to rather than a restatement of it, and the two agreeing is a statement
    /// about the arithmetic rather than about two people writing the same formula twice.
    ///
    /// **Successive over-relaxation with Newton on the diagonal**, matrix-free.
    /// `ThermalNetwork::steady_state` builds a dense Jacobian and factors it, which is right for a
    /// handful of nodes and impossible here: the bracket is 8 112 cells and its Jacobian would be
    /// 66 million entries. The discrete operator is diagonally dominant — a cell's diagonal is the
    /// sum of its face conductances plus whatever it loses to an environment, and its off-diagonals
    /// are those same face conductances — so sweeping in place converges, and the non-linear terms
    /// (a film's `T⁴`, a clearance's pair exchange) are taken at their slope, which is Newton.
    ///
    /// **Plain Gauss-Seidel is far too slow and it was measured being so.** Its spectral radius on
    /// an `n`-cell chain is `cos²(π/2n)`, so the sweeps needed grow as `n²`: a sixteen-cell bar
    /// took **over two thousand** and hit the bound. Over-relaxing by `ω = 2/(1 + sin(π/n))` — the
    /// optimal factor for a Poisson operator on that grid, with `n` the longest axis of the body
    /// being relaxed — makes it grow as `n` instead. The operator here is not pure Poisson, so
    /// that `ω` is an estimate rather than the optimum; it is clamped below 1.95 because
    /// over-relaxing past the true optimum diverges where under-relaxing only costs sweeps. It was
    /// the block's longest axis, and a body of two cells in a block twenty-nine long was taken past
    /// absolute zero by its first sweep.
    ///
    /// # Bodies, and the bodies a clearance joins
    ///
    /// **The slow mode is a body's level**, and each body that is a part on its own takes a step
    /// on it every round: what it loses over the slope of what it loses — **its films' own slope, at
    /// their surfaces**. It was the stability limit's tangent, at the hotter of each cell and its
    /// room, which reads a cell below its room as though its film were at the room: a sink fed by a
    /// radiating film stopped up to 1.75e-5 K from its balance near the most the film could bring,
    /// and one taking more than about three and a half times that ran out of sweeps, each round
    /// overshooting by more than it corrected. The sweep's diagonal reads the tangent at and above
    /// zero, where it bounds the slope, and the larger of the two below: the slope alone is Newton
    /// on each cell, which overshoots from below, and a strip radiating from four sides did not
    /// converge from 300 K to its 631.00 K. **The bodies a clearance
    /// joins take that step together**, Newton on all their levels at once, with what each pair
    /// carries moving with both its ends. Each body taking its own closes an exchange from both
    /// sides and moves it twice: two radiating bodies under a film too weak to hold them came to
    /// 2926.84 K in a room at 273.15 K and called it a balance, and under films a thousand to a
    /// hundred thousand times stronger did not converge. A cell whose pairs are most of what holds
    /// it takes a quarter of its own level a sweep at most, as a body does: three bodies of one
    /// cell each, one at 300 °C, settled at 2.7e17 K, the cold ones' first sweep asking for some
    /// 1250 K. A body's step is a quarter of the larger of its own level and the block's, so that a
    /// body near absolute zero is not held near it: of its own alone, a sink fed across a clearance
    /// found no balance above zero while the pairs read a negative temperature's fourth power as
    /// even, and in a first form of the step taken together a sealed pair at 400 K and 0 K
    /// answered 0 K for both. Read odd, as below, no block measured needs it. **A round whose step
    /// was held short has not settled**, however little it moved: two cells from 0 K in a room at
    /// 20 K, a step on their level a quarter of what little level they had, were answered at
    /// 4.59e-7 K. A round is a balance when what its steps asked for is below a microkelvin, as
    /// well as what they moved.
    ///
    /// # A part that loses nothing
    ///
    /// **Balances at any level, and the answer is the one holding the heat it has now.** A block
    /// insulated on every face, a body across a vacuum from the one with the film, an assembly
    /// whose bodies see each other and nothing else: every face term cancels under a uniform shift,
    /// so the equations fix the shape of the field and not where it sits. The march does fix it —
    /// it conserves `Σ C·T` — so that is what this holds: the level is the heat each such part
    /// holds now, over its capacity. **This answered wherever the relaxation drifted**, measured:
    /// four insulated cells, one at 100 °C and three at 20 °C, at 14.48 °C, colder than the coldest
    /// of them, where they march to 40 °C, and with a capacity four times the others' in the hot
    /// one, at 14.48 °C again rather than 65.71 °C; a body across a vacuum from a cooled one, at
    /// −8.05 °C in a room at 0 °C rather than its own 60 °C; two bodies radiating only to each
    /// other at 2653.69 °C rather than 110 °C. `settle` put a block there, and `verify` reported
    /// what was still to come against it.
    ///
    /// A part is the bodies joined by faces that conduct and by clearances that radiate, and it
    /// loses nothing when no film, air or mounting that carries reaches it through a cell that
    /// conducts — so a film of zero, a cell of zero conductivity and a joint of zero contact cut it
    /// off as a void does, and a cell that conducts nothing at all keeps its own temperature.
    ///
    /// # What it refuses
    ///
    /// **Heat generated where nothing carries it away.** There is no steady state; it warms without
    /// limit, and returning the last iterate would be a plausible temperature for a balance that
    /// was never struck. Asked of each part, naming a cell of it and carrying its watts. It was
    /// asked of the block — whether it had a film, air or a clearance anywhere — and a block of
    /// `k = 0` given a film returned the field it started from as its balance, while a film of
    /// zero, a heated body across a vacuum from the cooled one and a heated cell behind an
    /// insulator each swept two thousand times and reported not converging.
    ///
    /// **A balance below absolute zero.** A source of negative watts is a sink, and one that takes
    /// more than can reach it has no steady state; the field that balances it is not one a body
    /// can have. **The balance has one root, so one below zero says there is none above it.**
    /// Below zero a surface radiates `T·|T|³` — a pair's in `crossing`, a film's in
    /// `radiating_film` — so each term of a cell's balance falls with its own temperature and
    /// rises with its neighbours' there as above, and a change in one cell moves the others'
    /// balances by no more in all than it moves its own, the rest being what leaves its part;
    /// conduction, convection, air and mountings are linear and do so already. Between two fields
    /// that both balanced, the same holds of the slopes averaged along the way from one to the
    /// other, and summing their difference's signs against it leaves the difference of one sign
    /// over a whole part and none of it leaving: which a part that loses anything cannot have,
    /// and one that loses nothing cannot without holding other heat than it does. **Read even**,
    /// as the fourth power of a negative temperature,
    /// the pairs and films gave the balance a mirror below zero, and a sink fed across a clearance
    /// by a body a film warms from a room at 1000 K, drawn through zero before that heat arrived,
    /// settled on it. Measured on sinks fed that way and on sinks a film feeds directly, from six
    /// starts each: of ten with a balance, five runs of sixty were refused — a milliwatt answered
    /// at 878.56 K from five starts was refused at −898.56 K from the sixth, at 300 K — and of six
    /// without one, nineteen runs of thirty-six reported not converging and seventeen ran out of
    /// sweeps below zero anywhere from −0.32 K to −7009.6 K. Read odd, every run with a balance is
    /// answered, and
    /// every sink without one is refused at one temperature from all six starts — three milliwatts
    /// across the clearance at −630.06 K. **Held at zero below it**, as
    /// `ThermalNetwork::steady_state` holds a surface's radiation, the balance is flat there
    /// instead: measured, every sink without one was drawn down until the sweeps ran out, at
    /// −970 K to −87 120 K by start, and three runs with a balance were refused. A refusal where
    /// the sweeps ran out says so, and says only what was found. Within a tenth of a
    /// millikelvin below zero is the solve's own error around a room at absolute zero, and is
    /// answered as zero.
    ///
    /// **A block that melts.** The state of a mushy cell is its enthalpy and not its temperature,
    /// so "the field at which every cell balances" is not a statement about temperature alone.
    ///
    /// **A cell left out of balance with nothing to move it by.** A film with no convection to a
    /// room at absolute zero has no slope where its surface is at zero, and a cell held by nothing
    /// else is passed over there by every step the solve takes; a round that moves nothing else
    /// stops. A lone cell drawing 10 mW in a room at 0 K was answered at 0 K, 665.43 K above its
    /// root, from a start at zero. It is refused, naming the cell and carrying the watts.
    ///
    /// **Not converging**, which is a real answer and is not this one.
    ///
    /// **A number no block can have**, as [`step`](Domain::step) refuses it, by name. A face or
    /// air whose ambient was not a number made this **panic** rather than answer, inside the
    /// clamp that damps the uniform shift.
    pub fn steady_state(&self) -> Result<Vec<f64>, Violation> {
        let (nx, ny, nz) = self.counts;
        let dx = self.dx.to_si();
        if let Some(fault) = self.input_fault() {
            return Err(fault);
        }
        // The cells as the caller left them, as the step asks them: a part that loses nothing
        // comes to the heat it holds, and a cell below absolute zero holds heat no body can.
        if let Some(fault) = self.cell_fault() {
            return Err(fault);
        }
        if self.latent_anywhere {
            return Err(Violation::at(
                &self.name,
                "a block that melts has no steady temperature field: a mushy cell's state is its \
                 enthalpy, and this solves for temperature",
                f64::INFINITY,
            ));
        }
        // What the step generates, which a cell emptied after its source does not.
        let supplied = self.made_watts();
        // As the step refuses it, by what makes it unsteppable: this said "this block cannot be
        // stepped", for a substance with no thermal properties and for a bad joint alike.
        if self.conduction_rate().is_nan() {
            return Err(self.unsteppable());
        }

        // Which connected body each cell belongs to, six-connected across every face that
        // conducts. A block with a clearance in it holds two objects that conduct to nothing of
        // each other's, and so does one with an insulator through it, and each has a uniform mode
        // of its own — see the correction below.
        let body = self.bodies();
        let count = body.iter().filter_map(|b| *b).max().map_or(0, |m| m + 1);

        // **Which parts can lose heat**, and the heat generated in each — see `parts`. A part
        // that loses nothing and is given heat warms without limit, and is refused here rather
        // than swept two thousand times and reported as not converging; one given none keeps the
        // heat it holds, which the correction at the end of each round sees to.
        let (part, parts) = self.parts(&body, count);
        if let Some(p) = parts.iter().find(|p| !p.loses && !p.balanced()) {
            let (i, j, k) = (p.first % nx, (p.first / nx) % ny, p.first / (nx * ny));
            // Said the way it happens: a part whose sinks outweigh its sources cools without limit,
            // and was refused as warming.
            let why = if p.heat > 0.0 {
                format!(
                    "heat is generated in the cells joined to ({i}, {j}, {k}) and no face, void or \
                     clearance carries it away from them — nothing that loses heat is joined to \
                     them by anything that conducts — so there is no steady state: it warms \
                     without limit"
                )
            } else {
                format!(
                    "more heat is taken out of the cells joined to ({i}, {j}, {k}) than is \
                     generated in them and no face, void or clearance brings any in — nothing \
                     that loses heat is joined to them by anything that conducts — so there is no \
                     steady state: it cools without limit"
                )
            };
            return Err(Violation::at(&self.name, why, p.heat));
        }
        let closed = parts.iter().any(|p| !p.loses);

        // The bodies of each part, where each sits among them, and what each holds, J/K — the
        // unknowns of the step the bodies of a part take together, and its weights.
        let mut members: Vec<Vec<usize>> = vec![Vec::new(); parts.len()];
        let mut slot = vec![0; count];
        for b in 0..count {
            slot[b] = members[part[b]].len();
            members[part[b]].push(b);
        }
        let mut body_capacity = vec![0.0; count];
        let mut part_capacity = vec![0.0; parts.len()];
        for (c, b) in body.iter().enumerate() {
            if let Some(b) = b {
                body_capacity[*b] += self.capacity[c];
                part_capacity[part[*b]] += self.capacity[c];
            }
        }
        // **A bound on that step, refused rather than dropped.** It is a dense solve, so its
        // cost is the cube of the bodies a clearance joins into one part, and every shipped
        // scene has two at most. Past it a part would be left without its levels, which is
        // the drift this solve exists not to have.
        if let Some(bodies) = members.iter().find(|b| b.len() > STEADY_BODIES) {
            let first = body.iter().position(|b| *b == Some(bodies[0])).unwrap_or(0);
            let (i, j, k) = (first % nx, (first / nx) % ny, first / (nx * ny));
            return Err(Violation::at(
                &self.name,
                format!(
                    "the cells joined to ({i}, {j}, {k}) are {} bodies that see each other \
                     across clearances, and the steady solve couples {STEADY_BODIES} at most",
                    bodies.len()
                ),
                bodies.len() as f64,
            ));
        }

        // From where the block is. For the shape of the answer the start does not matter in a
        // linear problem, and in a radiative one it is a better guess than ambient whenever the
        // caller has already stepped, and no worse when they have not — the same reasoning
        // `ThermalNetwork::steady_state` gives. **For the level of a part that loses nothing it is
        // the answer**: such a part balances at any level, and the one its march reaches is the
        // one holding the heat it has now.
        let mut t = self.cells.clone();

        // **The over-relaxation factor, from each body's longest axis** rather than the block's:
        // what is relaxed is the body, and the factor is the one for a Poisson operator that long.
        // A body of two cells in a row of twenty-nine was over-relaxed by 1.81, as if it were
        // twenty-nine long, and its first sweep took a cell from 573 K to −114 K; past zero, its
        // pairs then read the fourth power of a negative temperature as a positive one, and
        // ten such bodies came to a balance with the first at −383.9 K and the rest at +383.9 K;
        // with fields below zero refused, they are refused as diverging. A body as long as its
        // block keeps the factor it had.
        let mut low = vec![(usize::MAX, usize::MAX, usize::MAX); count];
        let mut high = vec![(0, 0, 0); count];
        for (c, b) in body.iter().enumerate() {
            if let Some(b) = b {
                let (i, j, k) = (c % nx, (c / nx) % ny, c / (nx * ny));
                low[*b] = (low[*b].0.min(i), low[*b].1.min(j), low[*b].2.min(k));
                high[*b] = (high[*b].0.max(i), high[*b].1.max(j), high[*b].2.max(k));
            }
        }
        let omega: Vec<f64> = (0..count)
            .map(|b| {
                let (lo, hi) = (low[b], high[b]);
                let longest = (hi.0 - lo.0 + 1)
                    .max(hi.1 - lo.1 + 1)
                    .max(hi.2 - lo.2 + 1)
                    .max(2);
                (2.0 / (1.0 + (std::f64::consts::PI / longest as f64).sin())).clamp(1.0, 1.95)
            })
            .collect();

        // The net watts a clearance's pairs bring to each cell, and their slope. Held per sweep
        // rather than per cell: walking the list once fills both ends. The march's own crossing,
        // surfaces and all, so this is the balance it converges to.
        let mut gain = vec![0.0; t.len()];
        let mut slope = vec![0.0; t.len()];

        // What a part's sinks did to `c`. Where the sweeps settled, it is the balance's one root,
        // so the sinks take more than can reach them; this said "or more than reaches them in
        // time for this solve" while the pairs' even fourth powers gave it a mirror to settle on.
        // Where they had not, it is only what was found.
        let drawn_below = |c: usize, settled: bool| {
            let (i, j, k) = (c % nx, (c / nx) % ny, c / (nx * ny));
            if settled {
                format!(
                    "no balance above absolute zero was found: cell ({i}, {j}, {k}) came to \
                     below it, drawn down by the sinks in its part, which take more than can \
                     reach them"
                )
            } else {
                format!(
                    "no balance above absolute zero was found: cell ({i}, {j}, {k}) was below \
                     it, drawn down by the sinks in its part, where {STEADY_SWEEPS} sweeps had \
                     not settled"
                )
            }
        };

        for round in 0..STEADY_SWEEPS {
            for (a, b, coefficient) in &self.gaps {
                let pair = self.crossing_at(&t, *a, *b, *coefficient);
                let (at_a, at_b) = pair.slopes(*coefficient);
                gain[*a] -= pair.watts * pair.from_a;
                gain[*b] += pair.watts * pair.to_b;
                slope[*a] += at_a;
                slope[*b] += at_b;
            }

            let mut moved: f64 = 0.0;
            // The first cell this sweep left out of balance with nothing to move it by.
            let mut stuck: Option<(usize, f64)> = None;
            for k in 0..nz {
                for j in 0..ny {
                    for i in 0..nx {
                        let c = i + nx * (j + ny * k);
                        if self.void[c] {
                            continue;
                        }
                        // The residual in watts: the sweep's own operator, plus what the sweep
                        // applies beside it.
                        let r = self.flux_at(&t, (i, j, k), c) * dx + self.source[c] + gain[c];
                        // The diagonal: the face conductances this cell has, plus the slope of
                        // everything that leaves it. `face_sum` is conductivities, so it is the
                        // conductance divided by `dx`. **The limit's tangent, and below zero the
                        // films' slope where it is the larger** — see the step on a body's level,
                        // below. At or above zero the tangent is the larger, to the rounding, so the
                        // slope is not worked out there.
                        let tangent =
                            self.loss_conductance_at(Temperature::from_si(t[c]), (i, j, k));
                        let film = if t[c] < 0.0 {
                            tangent.max(self.loss_slope_at(t[c], (i, j, k), c))
                        } else {
                            tangent
                        };
                        let d = self.face_sum[c] * dx + film + slope[c];
                        // A cell with no conductance anywhere has no equation to solve, and a
                        // `NaN` one has no answer. Written out rather than as a negated
                        // comparison, which clippy asks for and is the clearer reading. **One
                        // left out of balance is kept**, for the exit to refuse: nothing here
                        // moves it, and a round that moves nothing else stops.
                        if d.is_nan() || d <= 0.0 {
                            if r != 0.0 && stuck.is_none() {
                                stuck = Some((c, r));
                            }
                            continue;
                        }
                        let mut step = body[c].map_or(1.0, |b| omega[b]) * r / d;
                        // **A cell held mostly by its pairs takes a quarter of its level at most**,
                        // as a body does below and for the same reason: a pair's slope is
                        // `4C·T³`, a poor guide far from the answer, and over-relaxed it is a
                        // worse one. Three bodies of one cell each across a radiating vacuum, one at
                        // 300 °C: the cold cells' first sweep asks for some 1250 K, past zero.
                        // Measured, they came to +12 772 K and −24 385 K and stopped moving, in a
                        // first form of this solve, and in this one they do not converge.
                        // A cell its faces or its film hold is not touched, so a block without
                        // clearances sweeps as it did. Its own level, unlike a body's below: a cell
                        // at absolute zero has no slope to its pairs, so it is not held here, and
                        // the body it is in is moved by the step below.
                        let cap = 0.25 * t[c].abs();
                        if slope[c] > 0.5 * d && cap >= 0.0 {
                            step = step.clamp(-cap, cap);
                        }
                        t[c] += step;
                        moved = moved.max(step.abs());
                    }
                }
            }

            for v in gain.iter_mut().chain(slope.iter_mut()) {
                *v = 0.0;
            }
            // Converged when a whole sweep moves the field by less than a microkelvin. Absolute
            // rather than relative to the span: a block settling *at* ambient has no span to be
            // relative to, and a microkelvin is four orders below anything any scene here asserts.
            // **The slow mode is a uniform shift, and it is one scalar equation.** Relaxation
            // converges at a rate set by the *smallest* eigenvalue, and for a block whose
            // conduction dwarfs what it loses that eigenvalue belongs to moving the whole field
            // together: a 4x4x4 aluminium block losing through one face has 0.0128 W/K to its air
            // against 0.668 W/K a face inside it, so the uniform mode decays by one part in a
            // thousand a sweep. Measured, the step fell from 0.104 K to 0.031 K over eighteen
            // hundred sweeps and the solve reported not converging.
            //
            // Face terms cancel under a uniform shift — every `T_j − T_i` is unchanged — so what
            // resists one is exactly what leaves the block. `Σr / Σ(loss slope)` is that balance,
            // and applying it kills the mode in a step. At the solution the residual is zero and
            // so is the shift, so this cannot move a converged field.
            //
            // **The loss slope is what the films do, at their surfaces** — `loss_slope_at`. It was
            // the stability limit's tangent, at the hotter of each cell and its room, which for a
            // body below its room reads `4εσA·T∞³` where the film moves by `4εσA·|s|³`: each round
            // then left about `1 − |s|³/T∞³` of the level still to go. Measured, on a sink a
            // radiating film from a room at 1000 K feeds: near the most the film can bring the
            // solve stopped 1.75e-5 K out at 50 mW, a round contracting by 0.946 and the stop at a
            // microkelvin leaving `ρ/(1 − ρ)` of one; from a little past `|s|³ = 2T∞³` each round
            // overshot by more than it corrected — 180 mW settled, and 185, 190 and 200 mW ran out
            // of sweeps — and a watt drawn from one cell reported that it did not converge. A
            // radiating sink with no balance, beside a film of 7 W/m²·K or none, ran out of sweeps
            // from 2 mW to 100 W, nine of nine. Now each settles, at 50 mW 5.1e-7 K out.
            //
            // **The sweep's diagonal is the larger of that slope and the tangent.** The tangent
            // bounds the films' slope wherever a surface is no further from zero than the hotter
            // of its cell and its room, so a cell at or above zero reads what the sweep always
            // read, to the bit, and one below zero the larger. Each alone fails somewhere. The
            // tangent alone reads a
            // cell far below zero at a fraction of its film: a watt drawn from a cell a film from
            // 1000 K alone holds was refused at its root, and two, five and a hundred reported not
            // converging. The slope alone is Newton on each cell, and from below a convex `T⁴` a
            // Newton step overshoots, over-relaxed by more: a strip of 64 cells radiating from four
            // sides did not converge from 300 K or 150 K to its 631.00 K, nor a lone cell from
            // 0 K or 1e-6 K to a room at 300 K, nor one from 1e-6 K to a room at 20 K; and a
            // heated cell in a room at 0 K, from 1e-12 K to 3e-3 K, came to exactly zero, where
            // nothing has a slope, and was refused 665.43 K from its balance. What the larger
            // costs is where the tangent reads a film as more of its cell's diagonal than it is:
            // on a row of 32 cells of 1.76 W/m·K and a centimetre drawing 9.92 W from the one a
            // film from 1200 K warms, the stop leaves 1.23e-5 K, as it did with the tangent in
            // the sweep alone, where the slope alone left 6.19e-6 K and the tangent in both
            // places 2.43e-5 K.
            let mut total = vec![0.0; count];
            let mut resist = vec![0.0; count];
            for k in 0..nz {
                for j in 0..ny {
                    for i in 0..nx {
                        let c = i + nx * (j + ny * k);
                        let Some(b) = body[c] else { continue };
                        total[b] += self.flux_at(&t, (i, j, k), c) * dx + self.source[c];
                        resist[b] += self.loss_slope_at(t[c], (i, j, k), c);
                    }
                }
            }
            // **A clearance's pair is a body's balance, not only its stiffness.** `flux_at` carries
            // conduction and the films; the pair exchange the sweep applies beside it has to be
            // added here too, or a body whose only path out is a clearance has a residual that
            // never falls. Measured: a bar radiating 0.02 W across a gap climbed to 1205 K and sat
            // there, its shift stuck at 2 K a sweep, because `total` said it was still holding all
            // 0.02 W however hot it got.
            //
            // A pair resists a shift only when its two ends are in *different* bodies: move one
            // object and the exchange changes, move both together and it does not. A pair inside
            // one body cancels, exactly as a conducting face does — without air. With it, what the
            // air takes between the two surfaces is the body's loss and stays in its total. Its
            // slope under a shift is left to `loss_slope_at`'s air, which carries all of it
            // where the two ends conduct alike and a few percent less where they do not: a slope
            // sets how fast the solve gets there, not where it stops.
            //
            // Between two bodies of one part the pair's couplings are kept as well, for the step
            // they take together below: `(a's body, b's body, a's rise with b, b's with a)`.
            let mut coupled: Vec<(usize, usize, f64, f64)> = Vec::new();
            for (a, b, coefficient) in &self.gaps {
                let (Some(ba), Some(bb)) = (body[*a], body[*b]) else {
                    continue;
                };
                let pair = self.crossing_at(&t, *a, *b, *coefficient);
                total[ba] -= pair.watts * pair.from_a;
                total[bb] += pair.watts * pair.to_b;
                if ba == bb {
                    continue;
                }
                let (at_a, at_b) = pair.slopes(*coefficient);
                resist[ba] += at_a;
                resist[bb] += at_b;
                if part[ba] == part[bb] && parts[part[ba]].bodies > 1 {
                    let (by_b, by_a) = pair.couplings(*coefficient);
                    coupled.push((ba, bb, by_b, by_a));
                }
            }

            // **Damped, because the slope of a `T⁴` term is a bad guide far from the answer.**
            // A body whose only path out is a clearance has a resistance of microwatts per
            // kelvin near ambient, so the first Newton step from a cold start asks for tens of
            // thousands of kelvin and the fourth power takes it to infinity. Measured: a bar
            // radiating 0.4 W through a pair asked for 93 000 K and the solve reported
            // diverging. Holding each step to a quarter of the body's own level makes the
            // approach geometric — a factor of 1.25 a sweep, so any starting point is tens of
            // sweeps away — and cannot pass through zero. The block's level, which `cap` takes
            // where it is the larger, can.
            let mut level = vec![0.0; count];
            let mut size = vec![0.0_f64; count];
            for (c, v) in t.iter().enumerate() {
                if let Some(b) = body[c] {
                    level[b] += v;
                    size[b] += 1.0;
                }
            }
            // **`clamp` panics on a bound that is not a number**, and a field that has stopped
            // being one has no level to hold a step to. It is the divergence the end of the
            // round reports, reached before the round ends. A cell given a temperature of
            // `NaN` reaches it with every input a block can have.
            //
            // **A quarter of the larger of its own level and the block's.** Of its own alone, a
            // body near absolute zero is held near it. Measured while the pairs read a negative
            // temperature's fourth power as even: a sink of two milliwatts fed across a clearance
            // by a body a film warms, from a start at 300 K, found no balance above zero, where
            // with this it settled near 685 K; and in a first form of the step the bodies of a
            // part take together, before its rows were brought to one, a sealed pair at 400 K and
            // 0 K answered 0 K for both, the cold one held still while the warm one was swept down
            // to it. **Read odd, nothing measured needs it**: that sink settles at 684.85 K from
            // every start with it and without it, and every test of this crate passes without it.
            // It stays until it is taken out on its own. A block of one body is its own reference,
            // so its bound is the one it had.
            let reference = cell_level(&t, &body);
            let cap = |b: usize| {
                let own = (level[b] / size[b].max(1.0)).abs();
                if own.is_nan() {
                    Err(Violation::at(
                        &self.name,
                        "the steady-state balance diverged",
                        f64::INFINITY,
                    ))
                } else {
                    Ok(0.25 * own.max(reference))
                }
            };
            // What each part that loses nothing has gained or lost since the start, as the
            // kelvin it would take over its capacity: `Σ C·(T₀ − T)`, put back below.
            let mut held = vec![0.0; parts.len()];
            if closed {
                for (c, v) in t.iter().enumerate() {
                    if let Some(b) = body[c] {
                        if !parts[part[b]].loses {
                            held[part[b]] += self.capacity[c] * (self.cells[c] - v);
                        }
                    }
                }
            }
            let mut shift = vec![0.0; count];
            // The most any level's step asked for this round, before it was held to `cap`: the
            // exit reads it beside what moved.
            let mut asked: f64 = 0.0;
            // A body that is a part on its own takes its own step: what it loses, over the slope
            // of what it loses.
            for b in 0..count {
                if parts[part[b]].bodies > 1 || resist[b] <= 0.0 {
                    continue;
                }
                let cap = cap(b)?;
                let wanted = total[b] / resist[b];
                asked = asked.max(wanted.abs());
                shift[b] = wanted.clamp(-cap, cap);
            }
            // **The bodies of a part of several take one step together**: Newton on their levels,
            // each body's balance against its own slope and the couplings of its pairs to the
            // others', and in a part that loses nothing the heat it holds as one more equation.
            // Each taking its own step, as a body alone does, closes an exchange from both ends
            // at once and moves it twice. Measured: two bodies across a radiating vacuum, one at
            // 100/20 °C and one at 300/20 °C, traded places every round and reported not
            // converging; a milliwatt in at one and out at the other, settling near 390 K and
            // 200 K, did not converge either; and under a film too weak to hold them, the same
            // two came to 2926.84 K and called it a balance, where the room is 273.15 K. The
            // step is held to a quarter of each body's level as one: every body's share of it
            // scaled together, so that what it says about their differences is kept. That bound
            // was written for the step a body took alone, and the 93 000 K it measured is not what
            // the step taken together asks for; what needed it was a sink fed across a clearance
            // from a body a film warms. Measured while the pairs read a negative temperature's
            // fourth power as even, two milliwatts drawn from a cell of `k = 0.05` that a
            // radiating film feeds through a clearance, from a start at 300 K: a balance near
            // 685 K with the bound, and not converging without it. Read odd, 684.85 K from every
            // start either way, and every test of this crate passes without it.
            for (p, bodies) in members.iter().enumerate() {
                if bodies.len() < 2 {
                    continue;
                }
                let closed_part = !parts[p].loses;
                let m = bodies.len();
                let n = m + usize::from(closed_part);
                let mut a = vec![0.0; n * n];
                let mut rhs = vec![0.0; n];
                for (s, b) in bodies.iter().enumerate() {
                    a[s * n + s] = resist[*b];
                    rhs[s] = total[*b];
                }
                for (ba, bb, by_b, by_a) in &coupled {
                    if part[*ba] == p {
                        a[slot[*ba] * n + slot[*bb]] -= by_b;
                        a[slot[*bb] * n + slot[*ba]] -= by_a;
                    }
                }
                // The heat the part holds, with a slack in every balance for what rounding
                // leaves of a net that is zero: the balances of a part that loses nothing add to
                // its net source, so without it they are one equation too many.
                //
                // **Each row in its own units, brought to one.** The balances are watts against
                // slopes in W/K, the heat row joules against capacities in J/K, and the solve
                // calls a pivot singular below 1e-14 of the largest entry of all: measured, two
                // sealed bodies of 8e3 J/K at 10 K and 20 K, whose slopes are 2e-11 W/K, were
                // judged singular, took no step, and answered their starting levels where their
                // heat sets 15 K. So each balance is divided by its body's slope, the slack by
                // their mean, and the heat row by the part's capacity.
                if closed_part {
                    let mut mean = 0.0;
                    for (s, b) in bodies.iter().enumerate() {
                        let r = if resist[*b] > 0.0 && resist[*b].is_finite() {
                            resist[*b]
                        } else {
                            1.0
                        };
                        for j in 0..m {
                            a[s * n + j] /= r;
                        }
                        a[s * n + m] = 1.0 / r;
                        rhs[s] /= r;
                        mean += r / m as f64;
                        a[m * n + s] = body_capacity[*b] / part_capacity[p];
                    }
                    for s in 0..m {
                        a[s * n + m] *= mean;
                    }
                    rhs[m] = held[p] / part_capacity[p];
                }
                let Some(step) = crate::network::solve(&mut a, &mut rhs, n) else {
                    // Judged singular all the same: a part that loses nothing still keeps its
                    // heat, moving as a whole, and the sweeps do the rest.
                    if closed_part {
                        for b in bodies {
                            shift[*b] = held[p] / part_capacity[p];
                        }
                    }
                    continue;
                };
                let mut scale: f64 = 1.0;
                for (s, b) in bodies.iter().enumerate() {
                    let cap = cap(*b)?;
                    asked = asked.max(step[s].abs());
                    if step[s].abs() > cap {
                        scale = scale.min(cap / step[s].abs());
                    }
                }
                for (s, b) in bodies.iter().enumerate() {
                    shift[*b] = scale * step[s];
                }
                // The heat row is scaled with the rest, so a round held short keeps only that
                // share of the heat it puts back; the next round asks for the rest, and at the
                // answer the scale is one. A quarter of a body's own level alone could make it zero
                // every round for a body at absolute zero — see `cap`, which is why that is a
                // quarter of the block's level too.
            }
            for (c, v) in t.iter_mut().enumerate() {
                if let Some(b) = body[c] {
                    *v += shift[b];
                }
            }
            moved = shift.iter().fold(moved, |m, s| m.max(s.abs()));

            // **A part that loses nothing keeps the heat it holds.** It balances at any level —
            // every face term cancels under a uniform shift, and nothing else acts on one — so the
            // equations do not say where it sits, and the relaxation left it wherever its sweeps
            // drifted. Measured before this: a block insulated on every face, one cell of four at
            // 100 °C and the rest at 20 °C, solved to 14.48 °C, colder than its coldest cell,
            // where its march goes to 40 °C; two radiating across a vacuum, at 2653.69 °C where
            // the march goes to 110 °C. What fixes the level is what the march conserves: the heat
            // in the part, `Σ C·T`, the cells' capacities and not their count. Each round puts
            // it back. At the answer it is already there, and this moves nothing. A part of
            // several bodies has put it back already, in the step they took together; this is
            // the part of one, which moves as a whole.
            if closed {
                for (c, v) in t.iter_mut().enumerate() {
                    if let Some(b) = body[c] {
                        let p = &parts[part[b]];
                        if !p.loses && p.bodies == 1 {
                            let shift = held[part[b]] / part_capacity[part[b]];
                            *v += shift;
                            moved = moved.max(shift.abs());
                        }
                    }
                }
            }

            // Before the exit rather than after it: `moved.max(NaN)` is `moved`, so a field that
            // had stopped being a number could stop moving and be returned.
            if t.iter().any(|v| !v.is_finite()) {
                return Err(Violation::at(
                    &self.name,
                    "the steady-state balance diverged",
                    f64::INFINITY,
                ));
            }
            // **And a level held short has not settled, however little it moved.** Its step is a
            // quarter of the larger of its own level and the block's, so a block near absolute
            // zero takes steps as small as it is. Measured: two cells from 0 K in a room at 20 K
            // were answered at 4.59e-7 K, where they balance at 20 K; four at 2.3e-7 K to
            // 3.7e-7 K, eight in a room at 10 K at 1.0e-8 K to 2.8e-8 K, and two bodies across a
            // clearance, a film on each, at 2.7e-7 K to 4.4e-7 K. What a step asked for has to be
            // below a microkelvin as well as what it moved. A cell's own quarter, in the sweep, is
            // not read: it binds only where the slope of a cell's pairs is more than its faces and
            // films carry, which near zero, where a pair's slope is all but nothing, is a cell with
            // neither — a body of one cell, whose level's step is read. **Nor is a step judged
            // singular**, which
            // asks for nothing and is not taken: a body only a clearance reaches has no slope to
            // its pair at exactly 0 K, and two such bodies from there, a film on the first from a
            // room at 20 K, are answered at 3.3e-7 K and 6.6e-40 K still — recorded as open.
            if round > 0 && moved < 1e-6 && asked < 1e-6 {
                // **A cell nothing moved is not balanced by that.** A cell whose every slope is
                // zero — no face that conducts, no pair, and a film with no convection to a room
                // at absolute zero, its surface at zero where `T·|T|³` is flat — is passed over
                // by the sweep and by the step on its level, and a round that moves nothing else
                // stops. Measured: a lone cell drawing 10 mW in a room at 0 K was answered at
                // 0 K from a start at zero, 665.43 K from its root below it; and from starts of
                // 1e-12 K and 1e-6 K the first sweep throws it — to −4.9e46 K from 1e-12 K — and
                // the next, reading the films' slope there, back to exactly zero. From 5e-99 K
                // to 5e-23 K, measured, that first throw is past where a fourth power overflows,
                // and the round reports the balance diverged: a refusal too, if a poorer word for
                // it. From 3e-99 K down the cell's slope is below the least `f64`, and it is
                // refused as left with nothing to move it by.
                if let Some((c, r)) = stuck {
                    let (i, j, k) = (c % nx, (c / nx) % ny, c / (nx * ny));
                    return Err(Violation::at(
                        &self.name,
                        format!(
                            "no balance was found: cell ({i}, {j}, {k}) is left {r:e} W out of \
                             balance at {:e} K, where nothing it exchanges has a slope to move it \
                             by",
                            t[c]
                        ),
                        r,
                    ));
                }
                // **A balance below absolute zero is not one.** A source of negative watts is a
                // sink, and a sink can take more out of a part than can reach it: in a part that
                // loses nothing, the level the heat it holds sets leaves the sink's cells below
                // zero. Measured: two bodies across a radiating vacuum, one given 0.01 W and the
                // other −0.01 W, balanced at −55.4 K — the field their march goes to, which
                // `step` does not refuse either. A shifted level converges there; before it, the
                // solve wandered and reported not converging.
                //
                // **And in a part with no sink it is the solve that went wrong**, which is said as
                // that. Nothing there can be colder than its coldest ambient, at or above zero, and
                // the balance's one root is not; a field below it is one the sweeps stopped at
                // without its being a root. It was a root of the pairs' even fourth powers. Asked of
                // the coldest solid cell, which in a part with a sink is the one the sinks draw
                // down, wherever the first cell happens to be.
                //
                // **What lies within the solve's own error of zero is zero.** A bar held at a room
                // at absolute zero by a film without limit at both ends stopped with a cell at
                // −8.8e-7 K and was refused as having diverged, where its answer is 0 K: the stop
                // leaves microkelvin either side, and on this side it is below a bound nothing
                // crosses. `BELOW_ZERO` is ten times what the tests here hold a cell to.
                let coldest = (0..t.len())
                    .filter(|c| body[*c].is_some())
                    .min_by(|a, b| t[*a].total_cmp(&t[*b]));
                if let Some(c) = coldest.filter(|c| t[*c] < -BELOW_ZERO) {
                    let (i, j, k) = (c % nx, (c / nx) % ny, c / (nx * ny));
                    let sinks = body[c].is_some_and(|b| parts[part[b]].sinks);
                    return Err(Violation::at(
                        &self.name,
                        if sinks {
                            drawn_below(c, true)
                        } else {
                            format!(
                                "the steady-state balance diverged, to cell ({i}, {j}, {k}) \
                                 below absolute zero, where no part without a sink can be"
                            )
                        },
                        t[c],
                    ));
                }
                for (c, v) in t.iter_mut().enumerate() {
                    if body[c].is_some() && *v < 0.0 {
                        *v = 0.0;
                    }
                }
                return Ok(t);
            }
        }

        // **Sinks that take more than can reach them draw their cells down without end**, and the
        // sweeps run out before anything settles: −50 mW out of a block whose film can bring it
        // 2 mW was refused as not converging. Said as what it is, as the converged case says it,
        // where the coldest cell of a part with a sink has gone below absolute zero.
        let coldest = (0..t.len())
            .filter(|c| body[*c].is_some())
            .min_by(|a, b| t[*a].total_cmp(&t[*b]));
        if let Some(c) = coldest
            .filter(|c| t[*c] < -BELOW_ZERO && body[*c].is_some_and(|b| parts[part[b]].sinks))
        {
            return Err(Violation::at(&self.name, drawn_below(c, false), t[c]));
        }
        Err(Violation::at(
            &self.name,
            "the steady-state balance did not converge in the sweeps allowed",
            supplied,
        ))
    }

    /// Which connected body each cell belongs to, or `None` for a void cell.
    ///
    /// Six-connected **across every face that conducts**. A block with a clearance in it holds two
    /// objects that conduct to nothing of each other's, and [`steady_state`](Solid3D::steady_state)
    /// needs to know which is which: a relaxation's slowest mode is moving a body as a whole, and
    /// there is one such mode per body.
    ///
    /// **A face of zero conductance divides a body as a void does.** An insulating fill and a joint
    /// of zero contact carry nothing, exactly as a clearance without radiation does, and each side
    /// has its own level. This walked the solid instead, so a block cut through by an insulator
    /// was one body to the correction — and a part of it with nothing to lose to was shifted by
    /// what the other part lost. Where no face is zero but at a void the bodies are the ones it
    /// found, in the same order.
    fn bodies(&self) -> Vec<Option<usize>> {
        let (nx, ny, nz) = self.counts;
        let mut label: Vec<Option<usize>> = vec![None; self.cells.len()];
        let mut next = 0;
        let at = |i: usize, j: usize, k: usize| i + nx * (j + ny * k);
        for k0 in 0..nz {
            for j0 in 0..ny {
                for i0 in 0..nx {
                    let start = at(i0, j0, k0);
                    if self.void[start] || label[start].is_some() {
                        continue;
                    }
                    label[start] = Some(next);
                    // Iterative, because a 128³ block is two million cells and a stack that deep
                    // is not a stack.
                    let mut stack = vec![(i0, j0, k0)];
                    while let Some((i, j, k)) = stack.pop() {
                        // `face` is the conductivity of the face crossed to reach the cell, which
                        // is zero at a void, so a void needs no test of its own.
                        let push = |i: usize,
                                    j: usize,
                                    k: usize,
                                    face: f64,
                                    s: &mut Vec<_>,
                                    l: &mut Vec<Option<usize>>| {
                            let c = at(i, j, k);
                            if face > 0.0 && l[c].is_none() {
                                l[c] = Some(next);
                                s.push((i, j, k));
                            }
                        };
                        let x = |i: usize| self.kx[i + (nx + 1) * (j + ny * k)];
                        let y = |j: usize| self.ky[i + nx * (j + (ny + 1) * k)];
                        let z = |k: usize| self.kz[i + nx * (j + ny * k)];
                        if i > 0 {
                            push(i - 1, j, k, x(i), &mut stack, &mut label);
                        }
                        if i + 1 < nx {
                            push(i + 1, j, k, x(i + 1), &mut stack, &mut label);
                        }
                        if j > 0 {
                            push(i, j - 1, k, y(j), &mut stack, &mut label);
                        }
                        if j + 1 < ny {
                            push(i, j + 1, k, y(j + 1), &mut stack, &mut label);
                        }
                        if k > 0 {
                            push(i, j, k - 1, z(k), &mut stack, &mut label);
                        }
                        if k + 1 < nz {
                            push(i, j, k + 1, z(k + 1), &mut stack, &mut label);
                        }
                    }
                    next += 1;
                }
            }
        }
        label
    }

    /// The closed parts of this block, as [`steady_state`](Solid3D::steady_state) needs them:
    /// which part each body is in, and for each part whether anything in it can lose heat, the
    /// heat generated in it, and the capacity that holds it.
    ///
    /// **A part is the bodies a clearance joins**, through a pair whose two ends both conduct — a
    /// pair behind a cell of zero conductivity carries nothing, as a face does. Radiation across a
    /// vacuum moves heat between two bodies and takes none away, so two bodies that see each other
    /// and nothing else are one closed system, and it is the system whose heat is conserved.
    ///
    /// **It loses heat** if any of its cells can — see `can_lose`. Which ones cannot is the
    /// question the old test asked of the whole block: whether it had a film, air or a clearance
    /// anywhere. A film that carries nothing passed it, and so did a heated body across a vacuum
    /// from the one with the film. Measured before this: a block of `k = 0` with a film and a
    /// heated cell returned the field it started from as its balance, and a film of zero, an
    /// area of zero, a mounting or air of zero, a heated cell behind an insulator and a heated body
    /// with no film beside one that had one each swept two thousand times and reported that the
    /// balance had not converged.
    fn parts(&self, body: &[Option<usize>], count: usize) -> (Vec<usize>, Vec<Part>) {
        let (nx, ny, _) = self.counts;
        // Union-find over the bodies, joined by the pairs that carry.
        let mut root: Vec<usize> = (0..count).collect();
        fn find(root: &mut [usize], mut b: usize) -> usize {
            while root[b] != b {
                root[b] = root[root[b]];
                b = root[b];
            }
            b
        }
        for (a, b, _) in &self.gaps {
            let (Some(ba), Some(bb)) = (body[*a], body[*b]) else {
                continue;
            };
            if self.derived[*a].0 > 0.0 && self.derived[*b].0 > 0.0 {
                let (ra, rb) = (find(&mut root, ba), find(&mut root, bb));
                root[ra.max(rb)] = ra.min(rb);
            }
        }
        // Numbered in the order their first cells come, so the first part refused is the one
        // holding the lowest cell.
        let mut number = vec![usize::MAX; count];
        let mut part = vec![0; count];
        let mut parts: Vec<Part> = Vec::new();
        for (c, b) in body.iter().enumerate() {
            let Some(b) = *b else { continue };
            let r = find(&mut root, b);
            if number[r] == usize::MAX {
                number[r] = parts.len();
                parts.push(Part {
                    loses: false,
                    heat: 0.0,
                    gross: 0.0,
                    sinks: false,
                    cells: 0,
                    bodies: 0,
                    first: c,
                });
            }
            part[b] = number[r];
            let p = &mut parts[number[r]];
            p.heat += self.source[c];
            p.gross += self.source[c].abs();
            p.sinks = p.sinks || self.source[c] < 0.0;
            p.cells += 1;
            p.loses = p.loses || self.can_lose((c % nx, (c / nx) % ny, c / (nx * ny)));
        }
        for b in 0..count {
            parts[part[b]].bodies += 1;
        }
        (part, parts)
    }

    /// Whether anything can leave this solid cell other than through its faces: the air in a void
    /// it faces, or a film on an outer face that convects or radiates, through whatever that face
    /// is mounted on and the half cell behind it.
    ///
    /// **Whether, not how much.** What `loss_conductance_at` charges is zero at every temperature
    /// exactly where this is false — a film, area, mounting or air of zero, a film with no
    /// convection on a face that does not radiate, a cell that conducts nothing — and asked
    /// separately because a radiating film's slope is zero at absolute zero, and a part that
    /// loses heat once it is warm is not one that cannot.
    fn can_lose(&self, cell: (usize, usize, usize)) -> bool {
        let (nx, ny, _) = self.counts;
        let (conductivity, _, emissivity) = self.derived[cell.0 + nx * (cell.1 + ny * cell.2)];
        let dx = self.dx.to_si();
        self.air_conductance(cell) > 0.0
            || self
                .exposed
                .iter()
                .filter(|(face, _)| self.exposed_at(**face, cell))
                .any(|(face, env)| {
                    let share = env.area.to_si() / self.cells_on(*face).max(1) as f64;
                    let acts = env.convection_w_per_m2_k > 0.0 || emissivity > 0.0;
                    // A film of one watt per square metre and kelvin in place of the one it has:
                    // what is asked here is whether the path behind it carries.
                    acts && series_with_half_cell(
                        self.outward(face, share, share),
                        conductivity,
                        share,
                        dx,
                    ) > 0.0
                })
    }

    /// Put this block **at** its steady state, and say how far it had to move.
    ///
    /// The answer in kelvin, as the largest change any cell took. A caller who wants to know
    /// whether a run has arrived can march it and then read this: a run that settled moves by
    /// nothing. A part of the block that loses nothing moves to where its own heat evens out, so
    /// its energy is the energy it had — see [`steady_state`](Solid3D::steady_state).
    pub fn settle(&mut self) -> Result<f64, Violation> {
        let field = self.steady_state()?;
        let mut moved: f64 = 0.0;
        for (c, v) in field.iter().enumerate() {
            moved = moved.max((v - self.cells[c]).abs());
            self.cells[c] = *v;
        }
        self.resolve();
        Ok(moved)
    }

    /// How many faces of the solid touch void, which is the area the air acts on.
    ///
    /// Zero for a block with no void, and for one whose void is all on the outside of it.
    pub fn faces_touching_void(&self) -> usize {
        let (nx, ny, nz) = self.counts;
        let mut n = 0;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if self.void[i + nx * (j + ny * k)] {
                        continue;
                    }
                    let mut touching = |i: usize, j: usize, k: usize| {
                        if self.void[i + nx * (j + ny * k)] {
                            n += 1;
                        }
                    };
                    if i > 0 {
                        touching(i - 1, j, k);
                    }
                    if i + 1 < nx {
                        touching(i + 1, j, k);
                    }
                    if j > 0 {
                        touching(i, j - 1, k);
                    }
                    if j + 1 < ny {
                        touching(i, j + 1, k);
                    }
                    if k > 0 {
                        touching(i, j, k - 1);
                    }
                    if k + 1 < nz {
                        touching(i, j, k + 1);
                    }
                }
            }
        }
        n
    }

    /// A **contact resistance** between an outer face and whatever cools it, W·m⁻²·K⁻¹.
    ///
    /// The boundary's half of [`joined`](Solid3D::joined), and the same physics on the one face
    /// that method cannot reach: `joined` puts a resistance between two cells, and a part bolted
    /// to a heatsink has one between its last cell and the film. Thermal grease is five to twenty
    /// thousand; a dry bolted joint a few hundred to a few thousand; a soldered baseplate
    /// hundreds of thousands.
    ///
    /// In series with the film **and** with the half cell of solid the film already sits behind,
    /// so the path from a boundary cell's centre outwards is
    ///
    /// ```text
    /// dx/(2kA)   +   1/(h_mount·A)   +   1/((h_conv + 4εσT³)·A)
    /// ```
    ///
    /// with the radiation read at the **outer face's** own temperature in the heat that moves,
    /// solved for so that what reaches that face is what it sheds, and as the tangent at the
    /// hotter of the centre and the room in the limit, which bounds it. It was read at the
    /// centre's in both, and the joint's drop does not shrink with the cell: a face radiating at
    /// 800 K behind `2000 W·m⁻²·K⁻¹` was charged at least 10.2 K hot at every grid and never
    /// converged.
    ///
    /// # Where the radiative term sits, and why it is behind the joint
    ///
    /// `Environment` derives radiation from the surface's emissivity, and that term stays in
    /// parallel with convection and in series with the joint — the coolant sees the *outer* face
    /// of the interface material. A surface buried under a filled joint does not radiate at all,
    /// so this is wrong in the limit of a hot, black, heavily greased face; it is worth stating
    /// because it is not worth measuring here. On the power module that ships with this
    /// workspace the cold face sits 5 K above a 40 °C plate and copper's emissivity is 0.03, so
    /// the radiative term is 0.03% of a 3000 W·m⁻²·K⁻¹ film.
    ///
    /// **Stability.** A mounting only lowers a boundary conductance, so the limit gets looser —
    /// the same argument [`joined`](Solid3D::joined) makes for an interior face.
    ///
    /// A face with no film is unaffected: there is nothing for the joint to be in series with,
    /// and a resistance to nowhere carries nothing either way.
    ///
    /// A negative or `NaN` conductance is refused by [`step`](Domain::step), naming the face. A
    /// negative one ran as a joint that carries nothing and a `NaN` one as no joint at all, and
    /// neither said so. `+∞` is the face bolted straight to its film.
    pub fn mounted_on(mut self, face: Face, w_per_m2_k: f64) -> Solid3D {
        self.mounted.insert(face, w_per_m2_k);
        self.resolve();
        self
    }

    /// The contact conductance between one outer face and its film, if one was stated.
    pub fn mounting_of(&self, face: Face) -> Option<f64> {
        self.mounted.get(&face).copied()
    }

    /// The first number this block was given that no block can have, as the refusal
    /// [`step`](Domain::step) and [`steady_state`](Solid3D::steady_state) give, or `None`.
    ///
    /// Asked of what was **stated**, every substance and every film, rather than of the operator
    /// `resolve` derived from them, because the operator is where these went quiet: a series mean
    /// guarded at zero reads a negative or `NaN` conductivity as an insulator, and `outward` reads
    /// a negative mounting as a joint that carries nothing. Each is a plausible block, so nothing
    /// downstream of it could tell. See `substance_fault` and `environment_fault` for the rest.
    ///
    /// **Every substance, whether or not a cell still holds it.** A block with a substance in it
    /// that cannot be one was built from a number somebody mistyped, and the cells it is in are
    /// not the place that mistake is. A substance with no thermal properties is the exception,
    /// as it always was: it is refused by the cells that hold it, as having no diffusivity, and a
    /// block filled over it entirely runs.
    ///
    /// A mounting or the air without limit is a face or a surface held at the film's ambient, and
    /// both are allowed, as an outer film without limit is. `NaN` is refused in each: a mounting
    /// of `NaN` was read as no mounting at all.
    fn input_fault(&self) -> Option<Violation> {
        let refuse = |what: String, value: f64| Some(Violation::at(&self.name, what, value));
        for s in &self.materials {
            if let Some((what, value)) = crate::substance_fault(s, true) {
                return refuse(what, value);
            }
        }
        // The cell size: zero, `NaN` and `+∞` were refused as a substance with no diffusivity, and
        // a cell a millimetre less than nothing ran — see `size_fault`.
        if let Some((what, value)) =
            crate::size_fault("the block's cell size", self.dx.to_si(), "metres")
        {
            return refuse(what, value);
        }
        // A source that is not a number of watts, by the cell it is in. `+∞` or `−∞` stepped its
        // cell there and answered `Ok`; `NaN` made the cell `NaN`, which the block's warmest and
        // coldest readings step over, so it read as a block at 292.92 K; and `steady_state` refused
        // each as a balance that diverged, or a `NaN` with nothing to lose it to as one that cools
        // without limit.
        if let Some((c, watts)) = self.bad_source {
            let (nx, ny) = (self.counts.0, self.counts.1);
            let (i, j, k) = (c % nx, (c / nx) % ny, c / (nx * ny));
            return refuse(
                format!("the source in cell ({i}, {j}, {k}) must be finite, in watts"),
                watts,
            );
        }
        for (face, environment) in &self.exposed {
            if let Some((what, value)) = crate::environment_fault(environment, true) {
                return refuse(format!("{face:?}: {what}"), value);
            }
        }
        for (face, h) in &self.mounted {
            if h.is_nan() || *h < 0.0 {
                return refuse(
                    format!(
                        "{face:?}: the mounting must be zero or more W/m²/K, or +inf for a face \
                         bolted straight to its film"
                    ),
                    *h,
                );
            }
        }
        if let Some((ambient, h)) = self.air {
            if h.is_nan() || h < 0.0 {
                return refuse(
                    "the air in the void must be zero or more W/m²/K, or +inf to hold its \
                     surfaces at its temperature"
                        .to_string(),
                    h,
                );
            }
            if !(ambient.is_finite() && ambient >= 0.0) {
                return refuse(
                    "the air in the void must be finite and at or above absolute zero, in kelvin"
                        .to_string(),
                    ambient,
                );
            }
        }
        None
    }

    /// Why a block whose conduction rate is not a number cannot be stepped, as
    /// [`step`](Domain::step) and [`steady_state`](Solid3D::steady_state) both refuse it.
    ///
    /// A bad joint and a substance with no conductivity both leave the limit unsayable, and they
    /// are not the same fault. Naming the first sends a caller to their materials, which is where a
    /// bad joint is not. The joint's refusal carries the joint, which it did not: it carried `+∞`,
    /// a value nobody stated.
    fn unsteppable(&self) -> Violation {
        match self.bad_contact {
            Some(h) if h.is_nan() => Violation::at(
                &self.name,
                "a joint's contact conductance is not a number",
                h,
            ),
            Some(h) => Violation::at(&self.name, "a joint's contact conductance is negative", h),
            None => Violation::at(&self.name, "substance has no diffusivity", f64::INFINITY),
        }
    }

    /// The first solid cell whose temperature nothing can be at — not a number, or below absolute
    /// zero — as the refusal [`step`](Domain::step) and [`steady_state`](Solid3D::steady_state)
    /// give, or `None`. See `state_fault`.
    ///
    /// **A void cell holds a placeholder, not a temperature, and has only to be a number.** One of
    /// `NaN` or `±∞` reached the solid beside it through a face that carries nothing, `0·(NaN − T)`
    /// being `NaN`, and the film's books by the same road: a void cell given `NaN`, or `deposit`ed
    /// into, whose capacity of nothing makes any heat infinite, stepped its neighbour to `NaN` and
    /// answered `Ok`, which the block's warmest and coldest readings step over.
    fn cell_fault(&self) -> Option<Violation> {
        let c = (0..self.cells.len()).find(|&c| {
            let t = self.cells[c];
            if self.void[c] {
                !t.is_finite()
            } else {
                crate::state_fault("", t).is_some()
            }
        })?;
        let (nx, ny) = (self.counts.0, self.counts.1);
        let (i, j, k) = (c % nx, (c / nx) % ny, c / (nx * ny));
        if self.void[c] {
            return Some(Violation::at(
                &self.name,
                format!(
                    "cell ({i}, {j}, {k}) is void, and what it holds must be a number: a face that \
                     carries nothing still multiplies it"
                ),
                self.cells[c],
            ));
        }
        crate::state_fault(
            &format!("cell ({i}, {j}, {k})'s temperature"),
            self.cells[c],
        )
        .map(|(what, value)| Violation::at(&self.name, what, value))
    }

    /// Generate `watts` **spread evenly over the cells `where_` selects**, replacing whatever
    /// those cells generated before.
    ///
    /// The total is the watts given, not watts per cell: a die dissipating 50 W dissipates 50 W
    /// whether the grid gives it eight cells or eight thousand, so a scene's answer does not move
    /// when its grid refines. That is what makes a source stated this way survive `verify`'s
    /// resolution sweep, and it is the opposite of the choice a per-cell figure would force.
    ///
    /// Void cells do not count toward the spread — nothing generates nothing — so a box drawn
    /// around a part and its clearance heats the part at the full rate rather than losing a share
    /// of it to the gap.
    ///
    /// That holds for the cells that are void when this is called. A cell emptied **after** it
    /// takes its share with it and the rest are not raised to cover it — see
    /// [`empty`](Solid3D::empty) — so empty first and generate second, which is the order a scene
    /// applies its keys in.
    ///
    /// **A void cell this selects is replaced too**, by nothing. It was skipped, and so kept the
    /// share an earlier call had given it before it was emptied, which filling it brought back:
    /// 5 W over four cells, the last emptied, then 5 W over all four again, made 6.25 W once the
    /// cell was filled, and 0 W over all four in place of the second call made 1.25 W.
    ///
    /// Selecting no solid cell is not an error here, because a caller building an assembly cell by
    /// cell passes through that state; the watts simply have nowhere to go and none are generated.
    /// A *scene* refuses it, because a scene saying 50 W and meaning none is a different mistake.
    pub fn dissipating(
        mut self,
        watts: f64,
        where_: impl Fn(usize, usize, usize) -> bool,
    ) -> Solid3D {
        let (nx, ny, nz) = self.counts;
        let mut chosen = Vec::new();
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let c = i + nx * (j + ny * k);
                    if !where_(i, j, k) {
                        continue;
                    }
                    if self.void[c] {
                        self.source[c] = 0.0;
                    } else {
                        chosen.push(c);
                    }
                }
            }
        }
        // Not returned from when nothing solid is chosen, as it was: a void cell replaced above
        // may have held the fault or a sink that are looked for again below.
        if !chosen.is_empty() {
            let each = watts / chosen.len() as f64;
            for c in chosen {
                self.source[c] = each;
            }
        }
        // Over every cell and not only these, because a later call that replaces a bad source
        // with a good one leaves no fault behind, and one that does not leaves the first.
        self.bad_source = self
            .source
            .iter()
            .position(|w| !w.is_finite())
            .map(|c| (c, self.source[c]));
        self.sinks = (0..self.source.len())
            .filter(|&c| self.source[c] < 0.0)
            .collect();
        self.resolve();
        self
    }

    /// How many **solid** cells a predicate selects.
    ///
    /// What a caller needs in order to refuse a source that would be generated nowhere: the block
    /// itself allows that state, because assembling cell by cell passes through it, and a scene
    /// does not.
    pub fn cells_on_where(&self, where_: &dyn Fn(usize, usize, usize) -> bool) -> usize {
        let (nx, ny, nz) = self.counts;
        let mut n = 0;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if !self.void[i + nx * (j + ny * k)] && where_(i, j, k) {
                        n += 1;
                    }
                }
            }
        }
        n
    }

    /// Whether anything in this block generates at all, so a block with no source pays nothing for
    /// the feature — which every scene written before it existed relies on.
    fn supplies(&self) -> bool {
        self.source.iter().any(|w| *w != 0.0)
    }

    /// Total watts generated, summed over the **solid** cells, which are the ones a step generates
    /// in.
    ///
    /// What a caller gets back is what it asked for, which is worth being able to check: a source
    /// spread over cells and then read back through a different route is exactly where a factor of
    /// the cell count hides.
    ///
    /// A cell emptied after a source was spread over it is not counted, because a step generates
    /// nothing there. This summed every cell: 5 W over four cells with the last emptied after read
    /// 5 W while each step made 3.75 W, so the one reading that could show watts going missing said
    /// none had — the books balance either way, since the joules are simply never made.
    pub fn generated_power(&self) -> Power {
        Power::from_si(self.made_watts())
    }

    /// The watts a step generates: every solid cell's source, and no void one's.
    ///
    /// Folded from `+0.0` rather than summed, because an empty `f64` sum is `−0.0` on stable and
    /// `+0.0` on 1.78, and a block that is all void would read a different zero by toolchain.
    fn made_watts(&self) -> f64 {
        self.source
            .iter()
            .zip(&self.void)
            .filter(|(_, void)| !**void)
            .fold(0.0, |total, (watts, _)| total + watts)
    }

    /// Joules generated over the run.
    pub fn generated_energy(&self) -> Energy {
        Energy::from_si(self.supplied)
    }

    /// Heat given up to the exposed faces' environments over the run.
    pub fn lost_energy(&self) -> Energy {
        Energy::from_si(self.lost)
    }

    /// Which faces are exposed, and to what.
    pub fn exposed_faces(&self) -> impl Iterator<Item = (Face, &Environment)> + '_ {
        self.exposed.iter().map(|(f, e)| (*f, e))
    }

    /// The stability rate at the block's **present** state: [`slope_rate_now`], each cell's row
    /// charged at the capacity it has, raised where a step can carry a cell across its melting point
    /// further than that capacity bounds — see `charge_crossing`.
    ///
    /// A block none of whose materials can be carried past its capacity — one phase, or a liquid
    /// holding more than its solid by too little for anything above absolute zero to reach, as
    /// octadecane's does — is not asked, and its limit is the slope's to the bit.
    ///
    /// [`slope_rate_now`]: Solid3D::slope_rate_now
    fn worst_rate_now(&self) -> f64 {
        let worst = self.slope_rate_now();
        let (freezes, melts) = self.reachable();
        if worst.is_nan() || (!freezes && !melts) {
            return worst;
        }
        self.charge_crossing(worst, freezes, melts)
    }

    /// Whether anything can reach the block's bound below, and its bound above: nothing a valid
    /// state holds is below `−ZERO_SLACK`, and no room or air below zero.
    fn reachable(&self) -> (bool, bool) {
        (
            self.freeze_below > -crate::ZERO_SLACK,
            self.melt_above < f64::INFINITY,
        )
    }

    /// The stability rate with every cell's row charged at the capacity it has now — the slope of
    /// its enthalpy: conduction, plus each exposed cell's loss conductance at that cell's own
    /// temperature, plus each clearance end's pairs at the hotter of each pair — all of a cell's
    /// added into its one row.
    ///
    /// State-dependent on purpose, which is what [`Domain::max_stable_dt`] means by "the
    /// largest step this domain can take **from `now`**". A block cooling from 1000 °C is
    /// stiffer at the start than at the end, and a limit that ignored that would be correct
    /// only at the end.
    ///
    /// Costs a pass over the cells of the exposed faces and, when it has air in its voids, over the
    /// solid cells facing them; and over the cells at the ends of its clearances, if it has any. With
    /// air it was a pass over every cell, 0.26 to 0.29 of a step. This said
    /// "the boundary cells" while the pass was over every cell, and the step asks for it before it
    /// moves anything, so on a small block that pass was 0.26 to 0.36 of what a step cost.
    fn slope_rate_now(&self) -> f64 {
        // Gaps as well as films: a pair of plates facing each other across void may have no
        // conducting face at all, so the radiative exchange is the *only* thing setting their
        // step. Returning early on `exposed` alone handed such a pair an infinite limit, and a
        // march at infinity is a NaN block reported as a substance with no diffusivity.
        let conducting = self.conduction_rate();
        if (self.exposed.is_empty() && self.gaps.is_empty() && self.air.is_none())
            || conducting.is_nan()
        {
            return conducting;
        }
        let mut worst = conducting;
        // A gap's radiative conductance, linearised at the **hotter** of the pair, which is the
        // conservative end: `dq/dT = 4σA T³/(1/ε₁+1/ε₂−1)` grows with temperature, so the hotter
        // cell's tangent bounds the pair's. Charged to both cells, because either could be the
        // one the step is too long for.
        //
        // **And of the two surfaces' sources, when air hotter than both centres is beside them.**
        // The pair radiates between its two surfaces, and both lie between the two sources `T*`
        // that hold them — each `T_c + w·(T_air − T_c)`, `w = H/(G + H)`, see `crossing` — so the
        // tangent at the hotter source bounds the secant between them. Air hotter than both
        // centres holds both sources above them, where the centres' tangent does not bound it;
        // air no hotter cannot, and then this is the centres' tangent to the bit. The sources and
        // not the air itself, which bounds them too and is far looser: air at 1500 K and
        // 50 W/m²·K over ends at 600 K and 400 K, in cells of 1 mm and `k = 0.6`, would charge the
        // pair 15.6 times the centres' tangent, where its hotter source charges it 1.19. The half
        // cells only add resistance to the path, so this charges more than the pair carries, never
        // less.
        //
        // **Added into each end's row, with everything else that leaves the cell** — its faces,
        // its film and its air — because the limit is a row's: `θᵢ = dt·ΣG/Cᵢ ≤ 1`. Each pair was
        // charged alone, as the largest of the cell's terms rather than their sum, so a cell that
        // also conducted, or faced two clearances, was given a step its row did not allow. Measured:
        // three 1 mm cells in a row at 900 K across two clearances reported 7.392 s where the middle
        // one's row allows 3.696 s, and forty steps at the limit it reported took cells that started
        // within 10 K of each other to 990.5 K and 708.9 K. Two pairs of 10 mm cells across one, each
        // face conducting what the pair radiates, went 179 K outside where they started in 40 steps.
        // The books balanced throughout: the exchange is antisymmetric.
        //
        // A cell's film and air are asked for only where it can have either: on the block's own
        // boundary, where an exposed face is, or with air in the void. Anywhere else the loss is a
        // sum over nothing, exactly zero, and asking for it at every end of a clearance through
        // the middle of a block cost a step 3.5% to 7.6%, 64³ to 16³, pinned.
        let dx = self.dx.to_si();
        let area = dx * dx;
        // The clearance's air as `crossing_at` reads it: W/K a face, and its temperature.
        let air = self
            .air
            .and_then(|(ambient, h)| (h > 0.0).then_some((h * area, ambient)));
        for end in &self.gap_ends {
            let c = end.cell;
            if self.capacity[c] > 0.0 {
                let radiated: f64 = self.gap_end_pairs[end.first..end.first + end.count]
                    .iter()
                    .map(|n| {
                        let (a, b, coefficient) = self.gaps[*n];
                        let centres = self.cells[a].max(self.cells[b]);
                        let mut hot = centres;
                        if let Some((h, ambient)) = air.filter(|(_, ambient)| *ambient > centres) {
                            for end in [a, b] {
                                let half = 2.0 * self.derived[end].0 * area / dx;
                                let w = if half.is_finite() {
                                    h / (half + h)
                                } else {
                                    0.0
                                };
                                let t = self.cells[end];
                                hot = hot.max(t + w * (ambient - t));
                            }
                        }
                        4.0 * coefficient * hot * hot * hot
                    })
                    .sum();
                let loss = if end.outer || self.air.is_some() {
                    self.loss_conductance_at(Temperature::from_si(self.cells[c]), end.at)
                } else {
                    0.0
                };
                worst = worst.max(self.mobility[c] * (self.face_sum[c] + (loss + radiated) / dx));
            }
        }
        // **Only the cells on an exposed face, and with air the solid cells facing the void**: any
        // other cell's loss is a sum over no face, exactly zero, and the pass over every cell spent
        // most of its time finding that out. A block's stability limit cost 0.26 of its step on
        // `24-a-power-module-junction-to-ambient`, whose one cooled face is 64 of its 512 cells,
        // and 0.36 on `29-a-designed-bracket-becomes-cells`.
        //
        // Exact, not near: the cells skipped are the ones the whole pass skipped by `loss == 0`
        // after computing their loss, and a maximum does not depend on the order its terms come
        // in. `worst` is not NaN here, and `max` returns the other side of a NaN, so it stays so.
        // A cell on two exposed faces, or on one and facing the air, is charged twice and counts
        // once.
        let planes = self
            .exposed
            .keys()
            .fold(worst, |w, face| self.charge(w, face.cells(self.counts)));
        // **The void's own cells are left out, and the whole pass charged them.** A void cell
        // conducts nothing, so there is no half cell for the air to reach it through and its air
        // conductance is exactly zero, as its mobility is: the whole pass skips it as `loss == 0`,
        // whatever the coefficient. A negative coefficient carries nothing on a solid cell
        // either, as a negative film on an outer face does, and an infinite one carries the half
        // cell, `2k·dx` a face. `air_in` takes both; `step` refuses the negative one before this
        // is asked, and the infinite one is a surface held at the air's temperature.
        let near = self.charge(planes, self.facing_air.iter().copied());
        // **And checked, on every debug run of every test that builds a block.** The two paths
        // rest on `loss_conductance_at` having no term anywhere else, and the shipped scenes barely
        // see the plane path at all: with its loop emptied, 33 of the 34 wrote the same bytes. A
        // term added there off the faces, the way the air was, would be left out of the limit and
        // make it too long, which does not diverge: it oscillates while the books balance.
        debug_assert_eq!(
            near.to_bits(),
            self.charge(worst, self.every_cell()).to_bits(),
            "{}: the exposed faces' cells and those facing the air set {near:e}, and every cell {:e}",
            self.name,
            self.charge(worst, self.every_cell())
        );
        near
    }

    /// `worst`, raised by every cell a step can carry across its melting point and past the coldest
    /// or the hottest thing it exchanges with: to the rate at which a step carries the cell's
    /// enthalpy toward that thing's, so that at the limit it lands on it. See `max_stable_dt`.
    ///
    /// `freezes` and `melts` say which of the block's two bounds anything can reach at all. One pass
    /// over the cells finds the coldest and hottest of them and of the rooms and the air, and a
    /// block that is nowhere past either bound pays that and nothing more.
    fn charge_crossing(&self, mut worst: f64, freezes: bool, melts: bool) -> f64 {
        let cells = self.cells_range.unwrap_or_else(|| self.range_of_cells());
        // **The kept range held to the cells, on every debug run**: a writer of a cell that does
        // not say so leaves it stale, and a stale range lets a cell past either bound unasked.
        if let (true, Some(kept)) = (cfg!(debug_assertions), self.cells_range) {
            let found = self.range_of_cells();
            assert!(
                kept.0.to_bits() == found.0.to_bits() && kept.1.to_bits() == found.1.to_bits(),
                "{}: a cell moved since `resolve` found the cells' range at {kept:?}, and nothing \
                 said so; they range over {found:?}",
                self.name
            );
        }
        let (lo, hi) = self
            .exposed
            .values()
            .map(|env| env.ambient.to_si())
            .chain(self.air.map(|(ambient, _)| ambient))
            .fold(cells, |(lo, hi), t| (lo.min(t), hi.max(t)));
        let freezes = freezes && lo < self.freeze_below;
        let melts = melts && hi > self.melt_above;
        if !freezes && !melts {
            return worst;
        }
        for c in 0..self.cells.len() {
            if !self.void[c] && self.latent[c] > 0.0 {
                worst = worst.max(self.crossing_rate(c, freezes, melts));
            }
        }
        worst
    }

    /// The coldest and the hottest solid cell, `(+∞, −∞)` for a block that is all void.
    fn range_of_cells(&self) -> (f64, f64) {
        self.cells
            .iter()
            .zip(&self.void)
            .filter(|(_, void)| !**void)
            .fold((f64::INFINITY, f64::NEG_INFINITY), |range, (t, _)| {
                widened(range, *t)
            })
    }

    /// A caller moved a cell: it is asked before a step reads it, and the range `resolve` found of
    /// the cells no longer holds.
    fn cells_moved(&mut self) {
        self.unchecked = true;
        self.cells_range = None;
    }

    /// What cell `c` asks of the limit beyond its row: where a step can freeze it past the coldest
    /// thing it exchanges with, `T_lo`, the rate `−F/(H(T) − H(T_lo))` at which the step carries its
    /// enthalpy `H` toward that thing's, `F` the watts the step moves into it; where a step can melt
    /// it past the hottest, `F/(H(T_hi) − H(T))`; and zero where neither.
    ///
    /// A step of `dt` leaves the cell at `H⁻¹(H(T) + dt·F)`, which is not below `T_lo` while
    /// `dt·(−F) ≤ H(T) − H(T_lo)`, and lands on it at equality. Within a step at the cell's row
    /// limit, `dt·ΣG ≤ C`, that holds in one phase, and across the melting point too unless `T_lo`
    /// is more than `L/(C_L − C_S)` below it — the same for a liquid cell at any temperature and a
    /// cell of any fraction, as the algebra of `H(T) + C·(T_lo − T) ≥ H(T_lo)` has it. So only a cell
    /// with liquid in it, of a liquid holding more than its solid, meeting something past that, is
    /// asked; and the same, the other way, of a cell with solid in it.
    ///
    /// `F` and not the row: a film or a clearance that radiates carries less than the tangent its
    /// row is charged at, from 900 K to 300 K 0.37 of it, and a bound charging the tangent shortened
    /// the limit of a cell that could not cross.
    fn crossing_rate(&self, c: usize, freezes: bool, melts: bool) -> f64 {
        let (point, latent) = (self.melt_point[c], self.latent[c]);
        let (solid, liquid) = (self.cap_solid[c], self.cap_liquid[c]);
        let colder = freezes && liquid > solid && self.melted[c] > 0.0;
        let warmer = melts && solid > liquid && self.melted[c] < 1.0;
        if !colder && !warmer {
            return 0.0;
        }
        // Past the cut-off, or the rate found is never the larger: asked to spare the work of what
        // the step carries, not for the answer. A cell with liquid in it is at or above its melting
        // point, so what lies past the cut-off is colder than it, and the same the other way.
        let (coldest, hottest) = self.input_range(c);
        let colder = colder && coldest < point - latent / (liquid - solid);
        let warmer = warmer && hottest > point + latent / (solid - liquid);
        if !colder && !warmer {
            return 0.0;
        }
        let watts = self.watts_into(c);
        let here = self.enthalpy_joules(c);
        if colder && watts < 0.0 {
            -watts / (here - solid * (coldest - point))
        } else if warmer && watts > 0.0 {
            watts / (latent + liquid * (hottest - point) - here)
        } else {
            0.0
        }
    }

    /// The watts a step moves into cell `c` at the present state: what its faces conduct in, less
    /// what its film and its air take, and what its clearances bring — the sweep's arithmetic and
    /// the clearances', read off `cells`. Neither a source nor the bus, which a step adds whatever
    /// its length and which `step` asks of on their own.
    fn watts_into(&self, c: usize) -> f64 {
        let (nx, ny, _) = self.counts;
        let cell = (c % nx, (c / nx) % ny, c / (nx * ny));
        let mut watts = self.dx.to_si() * self.flux_at(&self.cells, cell, c);
        if let Ok(n) = self.gap_ends.binary_search_by_key(&c, |end| end.cell) {
            let end = &self.gap_ends[n];
            for pair in &self.gap_end_pairs[end.first..end.first + end.count] {
                let (a, b, coefficient) = self.gaps[*pair];
                let crossing = self.crossing_at(&self.cells, a, b, coefficient);
                watts += if a == c {
                    -crossing.watts * crossing.from_a
                } else {
                    crossing.watts * crossing.to_b
                };
            }
        }
        watts
    }

    /// The coldest and the hottest of what cell `c` exchanges heat with in a step: each neighbour a
    /// face that carries anything joins it to, the room behind each of its exposed faces that
    /// carries anything, the air in the void where it touches it, and the far end of each of its
    /// clearances. `(+∞, −∞)` for a cell that exchanges with nothing.
    fn input_range(&self, c: usize) -> (f64, f64) {
        let (nx, ny, nz) = self.counts;
        let cell @ (i, j, k) = (c % nx, (c / nx) % ny, c / (nx * ny));
        let mut range = (f64::INFINITY, f64::NEG_INFINITY);
        let mut meet = |t: f64| range = (range.0.min(t), range.1.max(t));
        let (xr, yr) = ((nx + 1) * (j + ny * k), nx * (j + (ny + 1) * k));
        if i > 0 && self.kx[i + xr] > 0.0 {
            meet(self.cells[c - 1]);
        }
        if i + 1 < nx && self.kx[i + 1 + xr] > 0.0 {
            meet(self.cells[c + 1]);
        }
        if j > 0 && self.ky[i + yr] > 0.0 {
            meet(self.cells[c - nx]);
        }
        if j + 1 < ny && self.ky[i + nx + yr] > 0.0 {
            meet(self.cells[c + nx]);
        }
        if k > 0 && self.kz[c] > 0.0 {
            meet(self.cells[c - nx * ny]);
        }
        if k + 1 < nz && self.kz[c + nx * ny] > 0.0 {
            meet(self.cells[c + nx * ny]);
        }
        let (conductivity, _, emissivity) = self.derived[c];
        let here = Temperature::from_si(self.cells[c]);
        for (face, env) in &self.exposed {
            if self.exposed_at(*face, cell)
                && self.face_loss(face, env, here, conductivity, emissivity) > 0.0
            {
                meet(env.ambient.to_si());
            }
        }
        if let Some((ambient, _)) = self.air {
            if self.air_conductance(cell) > 0.0 {
                meet(ambient);
            }
        }
        // The far end of a clearance, where its centre reaches the pair: through a half cell that
        // conducts, on a pair that radiates, and not through air without limit, which holds each
        // surface at the air's own temperature — met above. A far end of a perfect insulator met
        // here would be a colder thing the cell does not meet, and let it past what it does.
        let held = self.air.is_some_and(|(_, h)| h == f64::INFINITY);
        if let Ok(n) = self.gap_ends.binary_search_by_key(&c, |end| end.cell) {
            let end = &self.gap_ends[n];
            for pair in &self.gap_end_pairs[end.first..end.first + end.count] {
                let (a, b, coefficient) = self.gaps[*pair];
                let far = if a == c { b } else { a };
                if coefficient > 0.0
                    && self.derived[c].0 > 0.0
                    && self.derived[far].0 > 0.0
                    && !held
                {
                    meet(self.cells[far]);
                }
            }
        }
        range
    }

    /// `worst`, raised by each of `cells` whose loss is not zero: its conduction and its loss at
    /// its own temperature, as a rate.
    fn charge(&self, mut worst: f64, cells: impl Iterator<Item = (usize, usize, usize)>) -> f64 {
        let (nx, ny, _) = self.counts;
        let dx = self.dx.to_si();
        for (i, j, k) in cells {
            let c = i + nx * (j + ny * k);
            let loss = self.loss_conductance_at(self.temperature_at(i, j, k), (i, j, k));
            if loss != 0.0 {
                worst = worst.max(self.mobility[c] * (self.face_sum[c] + loss / dx));
            }
        }
        worst
    }

    /// Every cell of the block, in index order.
    fn every_cell(&self) -> impl Iterator<Item = (usize, usize, usize)> {
        let (nx, ny, nz) = self.counts;
        (0..nz).flat_map(move |k| (0..ny).flat_map(move |j| (0..nx).map(move |i| (i, j, k))))
    }

    /// A boundary film's conductance with the mounting in series, both in W/K.
    ///
    /// **Two callers and one function, deliberately.** The boundary loss is computed twice — once
    /// by `loss_conductance_at` for the stability limit and once by `film_flux` for the heat that
    /// actually moves — and the first version of `mounted_on` changed only the first. The limit
    /// got looser, the flux did not, and a face bolted on through a perfect insulator shed 1004 J
    /// while reporting a step it had earned. Nothing about that looks wrong from either side.
    fn outward(&self, face: &Face, film: f64, share: f64) -> f64 {
        match self.mounted.get(face) {
            // A joint of zero carries nothing, and nothing is not a small number: a block that
            // cannot lose heat has no steady state, and one that loses a little has one far away.
            Some(m) if *m <= 0.0 => 0.0,
            Some(m) if m.is_finite() => 1.0 / (1.0 / film + 1.0 / (m * share)),
            // `+∞`, or no mounting at all: the film as it was, bit for bit rather than through a
            // reciprocal that need not round-trip.
            _ => film,
        }
    }

    /// What the clearance pair `(a, b)` carries at the temperatures `t`: each end behind its own
    /// half cell, `2k·dx` with the conductivity its faces are built from, and the clearance's air
    /// on both faces if the block has any — see `crossing`. The half cell is formed as
    /// `series_with_half_cell` forms it, so the air's share of the face here is the one
    /// `air_conductance` charges.
    fn crossing_at(&self, t: &[f64], a: usize, b: usize, coefficient: f64) -> Crossing {
        let dx = self.dx.to_si();
        let area = dx * dx;
        let half = |c: usize| 2.0 * self.derived[c].0 * area / dx;
        let air = self
            .air
            .and_then(|(ambient, h)| (h > 0.0).then_some((h * area, ambient)));
        crossing(coefficient, (half(a), t[a]), (half(b), t[b]), air)
    }

    /// The conductance from one solid cell to the air in the void touching it, in W/K.
    ///
    /// **One function, two callers**, for the reason `outward` above is: the boundary loss is
    /// computed once for the stability limit and once for the heat that actually moves, and a term
    /// that reaches one and not the other is a defect with no symptom on either side alone.
    ///
    /// Counted from the grid rather than from a stated area. A caller states what an *outer* face
    /// exposes because a rasterised part covers less of it than the grid does; the faces touching
    /// an internal void are exactly the ones the grid has, so there is nothing to state and
    /// nothing to state wrongly.
    ///
    /// **Each face's film in series with the half cell behind it**, `2k·dx`, as an outer face's
    /// is — with the conductivity the cell's faces are built from, from `derived`. It was
    /// `h·dx²` a face, charged against the cell centre, half a cell from the surface the air
    /// touches and a half cell's drop away from its temperature: first order. In a steady slab
    /// every cell facing the air sat too cold by exactly that drop, `q''·dx/2k`, or too hot
    /// where the heat came in.
    /// Measured on a slab 16 mm deep, `k = 0.6`, heated to shed 100 kW/m² into air at `h = 200`,
    /// 333.3 K too cold at 4 mm a cell and 41.67 K at 0.5 mm, halving with the cell, where the
    /// same slab cooled through an outer face met the scheme's own closed form to 1e-13 K at the
    /// face and to the steady solve's stop at every cell.
    /// `tests/the_air_in_a_void_cools_at_the_surface.rs` holds it there now.
    fn air_conductance(&self, cell: (usize, usize, usize)) -> f64 {
        let Some((_, h)) = self.air else {
            return 0.0;
        };
        let (nx, ny, nz) = self.counts;
        let (i, j, k) = cell;
        let dx = self.dx.to_si();
        let conductivity = self.derived[i + nx * (j + ny * k)].0;
        let mut faces = 0usize;
        let mut touching = |i: usize, j: usize, k: usize| {
            if self.void[i + nx * (j + ny * k)] {
                faces += 1;
            }
        };
        if i > 0 {
            touching(i - 1, j, k);
        }
        if i + 1 < nx {
            touching(i + 1, j, k);
        }
        if j > 0 {
            touching(i, j - 1, k);
        }
        if j + 1 < ny {
            touching(i, j + 1, k);
        }
        if k > 0 {
            touching(i, j, k - 1);
        }
        if k + 1 < nz {
            touching(i, j, k + 1);
        }
        if faces == 0 {
            return 0.0;
        }
        series_with_half_cell(h * dx * dx, conductivity, dx * dx, dx) * faces as f64
    }

    /// The film's flux out of a cell, in the sweep's own units — a conductivity times a
    /// temperature difference, so that `mobility · this` is a rate in kelvin per second.
    ///
    /// The **secant** conductance of the bare surface, which is the exact one for this flux,
    /// put in series with the half cell behind it. Positive means heat leaving. A face that
    /// radiates reads that secant at the surface's own temperature, which `radiating_film`
    /// solves for; one that only convects has a secant that is the same at every temperature,
    /// and keeps the arithmetic it always had.
    ///
    /// The half cell conducts as the cell's faces do, and the surface radiates as its phase does,
    /// both at the phase the operator was derived at — see `loss_conductance_at` for what reading
    /// the solid's cost, and `derived` for why not the phase of the moment.
    fn film_flux(&self, old: &[f64], cell: (usize, usize, usize), c: usize) -> f64 {
        let (conductivity, _, emissivity) = self.derived[c];
        let dx = self.dx.to_si();
        let here = Temperature::from_si(old[c]);
        let mut out = 0.0;
        for (face, env) in &self.exposed {
            if !self.exposed_at(*face, cell) {
                continue;
            }
            let share = env.area.to_si() / self.cells_on(*face).max(1) as f64;
            let gap = here.to_si() - env.ambient.to_si();
            if gap == 0.0 {
                continue;
            }
            let g = if emissivity == 0.0 {
                // Linear, so the film's conductance is one number whatever the surface's
                // temperature, and the series below is exact as it stands.
                let bare = Environment {
                    ambient: env.ambient,
                    convection_w_per_m2_k: env.convection_w_per_m2_k,
                    area: Area::from_si(share),
                }
                .loss_from(here, emissivity)
                .to_si();
                series_with_half_cell(
                    self.outward(face, bare / gap, share),
                    conductivity,
                    share,
                    dx,
                )
            } else {
                let half = 2.0 * conductivity * share / dx;
                let behind = if !half.is_finite() {
                    self.outward(face, f64::INFINITY, share)
                } else if half > 0.0 {
                    self.outward(face, half, share)
                } else {
                    0.0
                };
                radiating_film(
                    behind,
                    share,
                    env.convection_w_per_m2_k,
                    emissivity * STEFAN_BOLTZMANN.to_si(),
                    here.to_si(),
                    env.ambient.to_si(),
                )
            };
            out += g * gap / dx;
        }
        // The air in the void, if this block has any. **Convection only**: a transparent gas does
        // not radiate, and what a surface facing a clearance exchanges with the surface across it
        // is `find_gaps`' pairing, which runs beside this rather than instead of it.
        if let Some((ambient, _)) = self.air {
            let g = self.air_conductance(cell);
            if g > 0.0 {
                out += g * (old[c] - ambient) / dx;
            }
        }
        out
    }

    /// How many **solid** cells lie on a face.
    ///
    /// **It counted grid cells, and a designed part is mostly not there.** The stated `area_cm2`
    /// is divided among the cells on the face so that each carries its share; a void cell takes a
    /// share and loses nothing with it, so the area that actually cools is the stated one times
    /// the solid fraction of that face. Measured on `29-a-designed-bracket-becomes-cells`: 410 of
    /// the 676 cells on `z-min` are solid, so 16.5 cm² of stated cooling acted as **10.0**, and
    /// the block's 120 s drop came out 1.908 K against the 2.950 K its own lumped balance
    /// predicts — a ratio of 0.647 that tracks the area shortfall and nothing else.
    ///
    /// `cells_on_where` directly above has filtered voids since it was written, for the same
    /// reason and one line away.
    ///
    /// Read from a cache rather than counted, because this is called per cell per face per step.
    /// `void` is written in exactly one place — [`Solid3D::empty`] — and that calls `resolve`,
    /// which is where the cache is filled.
    fn cells_on(&self, face: Face) -> usize {
        self.solid_on[Solid3D::face_index(face)]
    }

    /// How many **solid** cells of a face are inside its patch, which is the number the stated
    /// area is divided among.
    ///
    /// Counted rather than cached, because a patch is stated once and this is read per cell per
    /// step only through `cells_on`, which is cached — see `count_faces`, where the patch is
    /// applied.
    fn count_patched(&self, face: Face) -> usize {
        let (nx, ny, nz) = self.counts;
        let mut n = 0;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if !self.void[i + nx * (j + ny * k)] && self.exposed_at(face, (i, j, k)) {
                        n += 1;
                    }
                }
            }
        }
        n
    }

    /// Whether this cell is on `face` **and** inside whatever patch that face states.
    ///
    /// The single place that answers it. `holds` alone was consulted from three: the flux, the
    /// stability limit and the count that divides the area — and a patch that reached two of the
    /// three would shed heat from cells it had excluded, or divide by a count it did not use.
    fn exposed_at(&self, face: Face, cell: (usize, usize, usize)) -> bool {
        if !face.holds(cell, self.counts) {
            return false;
        }
        match self.patch.get(&face) {
            None => true,
            Some([lo0, lo1, hi0, hi1]) => {
                let p = [cell.0, cell.1, cell.2];
                let [a, b] = face.axes();
                p[a] >= *lo0 && p[a] < *hi0 && p[b] >= *lo1 && p[b] < *hi1
            }
        }
    }

    /// Where a face's count sits in `solid_on`.
    fn face_index(face: Face) -> usize {
        match face {
            Face::XMin => 0,
            Face::XMax => 1,
            Face::YMin => 2,
            Face::YMax => 3,
            Face::ZMin => 4,
            Face::ZMax => 5,
        }
    }

    /// Count the solid cells on each of the six faces. Called by `resolve`.
    fn count_faces(&mut self) {
        let (nx, ny, nz) = self.counts;
        let solid = |s: &Solid3D, i: usize, j: usize, k: usize| !s.void[i + nx * (j + ny * k)];
        let mut out = [0usize; 6];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if !solid(self, i, j, k) {
                        continue;
                    }
                    if i == 0 {
                        out[0] += 1;
                    }
                    if i + 1 == nx {
                        out[1] += 1;
                    }
                    if j == 0 {
                        out[2] += 1;
                    }
                    if j + 1 == ny {
                        out[3] += 1;
                    }
                    if k == 0 {
                        out[4] += 1;
                    }
                    if k + 1 == nz {
                        out[5] += 1;
                    }
                }
            }
        }
        // **A patch changes the divisor, not only the set.** The stated area is spread over the
        // cells on a face, so a face with a patch on it must divide by the cells in the *patch* —
        // otherwise a bolt pad of nine cells on a face of six hundred carries a sixty-seventh of
        // the conductance the caller stated, and the run completes reporting four figures.
        for face in Face::ALL {
            if self.patch.contains_key(&face) {
                out[Solid3D::face_index(face)] = self.count_patched(face);
            }
        }
        self.solid_on = out;
    }

    /// The loss conductance an exposed cell carries, in W/K, summed over the faces it lies on.
    ///
    /// Evaluated at the temperature handed in, and the callers hand in **the cell's own current
    /// temperature** rather than the ambient. That is not a refinement, it is the difference
    /// between a limit that holds and one that does not: the radiative conductance is `4εσT³`
    /// and a part at 1000 °C carries **26 times** what the same surface carries at room
    /// temperature — measured, `ε = 0.9`, 136 against 5.14 W/m²·K. A limit linearised about
    /// ambient would hand a hot block a step an order too long, and an explicit boundary past
    /// its limit does not diverge loudly: it oscillates about the air it is losing to while the
    /// conservation audit stays perfectly happy, because the heat really did leave.
    ///
    /// Two things make it conservative rather than merely plausible, and each was measured
    /// wrong first.
    ///
    /// The **derivative** `4εσT³` rather than the secant `q/(T−T∞)`, because `T⁴` is convex and
    /// the derivative is the larger of the two — four times it in the hot limit.
    ///
    /// And at `max(T_cell, T_ambient)`, because the derivative only dominates the secant for
    /// `T ≥ T∞`. Below ambient it reverses without bound, and the first version of this used
    /// the cell alone: a 20 mm cell at 20 °C inside a 700 °C radiant enclosure was handed a
    /// step whose true ratio was **13**, and one accepted step took it to 8847 °C. `4T∞³`
    /// bounds `(T+T∞)(T²+T∞²)` for every `T ≤ T∞`, so the larger of the two points is safe on
    /// both sides.
    ///
    /// The film is put **in series with the half cell between the cell centre and the surface**,
    /// `2k·A/dx`. A finite-volume cell knows only its centre and the film acts at the surface, so
    /// the half cell between them is part of the path; charging the whole film against the centre
    /// sheds too much, by the cell Biot number `h·dx/2k`.
    ///
    /// **That `k` is the one the cell's faces are built from**, which for a substance with a
    /// liquid phase is mixed by what has melted — read from `derived`, so that inside a step it is
    /// the phase the faces and the limit have. This and `film_flux` used to read the solid's
    /// whatever had melted. Measured on one cell of water under `h = 1e4`: a limit of 1.247 s
    /// where its own is 3.590 s, and a step that shed 13.80 mJ where the film carries 4.79 mJ —
    /// 2.88 times, the path through ice's half cell against the path through water's. Half
    /// melted, 1.40 times. The books balanced, because the heat really did leave.
    /// `tests/a_melted_surface_loses_through_its_liquid.rs` holds both to the closed form.
    ///
    /// **And the emissivity is the phase's**, mixed by what has melted and read from `derived` for
    /// the same reason: a finish is a surface, and a surface that has melted is the liquid's. It was
    /// the solid's whatever had melted, here and in `film_flux` alike, so the two agreed and the
    /// block was stable and wrong. Measured on one cell of liquid at 900 K radiating to 300 K, its
    /// solid 0.9 and its liquid 0.1: a limit of 30.21 s where its own is 245.2 s, and a step that
    /// shed 8.65 times what its own surface sends. The other way about, 0.1 and 0.9, it shed 0.116
    /// times. `tests/a_melted_surface_radiates_as_its_liquid.rs` holds both to the closed form.
    ///
    /// **And it took two corrections to reach second order, not one.** The series form fixes
    /// the *space*: without it a review measured 2.08, 2.04, 2.02 per grid doubling against the
    /// interior's 4.0 — the signature this workspace has twice caught in an acoustic wall. With
    /// it the boundary was still first order, at 1.27, 1.71, 1.87, because the film was applied
    /// as a pass *after* the conduction sweep. That split's error carries a coefficient growing
    /// as `1/dx` while the step falls as `dx²`, so their product falls as `dx`. Applying the
    /// film inside the same explicit update — see `film_flux`, called from the sweep — makes
    /// them one operator, and `tests/the_cooled_boundary_order.rs` measures four per doubling at
    /// `Bi = 1.7` and at `Bi = 17`.
    ///
    /// **That is a convective face, and a radiating one needed a third.** After both corrections
    /// its film was still read at the centre's temperature, which a half cell's drop puts above
    /// the surface's, so it shed too much by an amount proportional to `dx` — 1.87 to 1.98 a
    /// halving. `film_flux` now reads it at the surface's, through `radiating_film`, and
    /// `tests/a_radiating_face_sheds_from_its_surface.rs` measures 4.000. This limit is unchanged:
    /// the tangent at the hotter of the centre and the room bounds the secant at the surface on
    /// either side of the room.
    ///
    /// **Nothing here may be non-zero for a solid cell off an exposed face that does not face the
    /// air.** `worst_rate_now` visits only those faces' cells and, with air, `facing_air`, and a
    /// term reaching any other cell must widen that pass, or the limit leaves it out. Its debug
    /// assertion is what says so.
    ///
    /// **A bound, and not the slope**: below the room it is more than a radiating film's, and it
    /// bounds it wherever the surface is no further from zero than the hotter of the cell and the
    /// room — every cell at or above zero. The steady solve reads
    /// [`loss_slope_at`](Solid3D::loss_slope_at) instead in the step on a body's level, and in its
    /// sweep the larger of the two for a cell below zero.
    fn loss_conductance_at(&self, at: Temperature, cell: (usize, usize, usize)) -> f64 {
        let (nx, ny, _) = self.counts;
        let (conductivity, _, emissivity) = self.derived[cell.0 + nx * (cell.1 + ny * cell.2)];
        self.air_conductance(cell)
            + self
                .exposed
                .iter()
                .filter(|(face, _)| self.exposed_at(**face, cell))
                .map(|(face, env)| self.face_loss(face, env, at, conductivity, emissivity))
                .sum::<f64>()
    }

    /// How fast what cell `c`'s films and air carry away rises with its own temperature `t`, W/K:
    /// the derivative of [`film_flux`](Solid3D::film_flux), read off the same faces, shares,
    /// half cells and mountings, for the step [`steady_state`](Solid3D::steady_state) takes on a
    /// body's level, and for its sweep below zero, where it can be the larger. A face that does
    /// not radiate is linear, and its slope is what [`face_loss`](Solid3D::face_loss) charges;
    /// summed as [`loss_conductance_at`](Solid3D::loss_conductance_at) sums, so that a cell none
    /// of whose faces radiates reads it to the bit. One that radiates is
    /// [`radiating_film_slope`]'s.
    /// **Not a bound**, as `loss_conductance_at` is for the stability limit: wherever a surface is
    /// no further from zero than the hotter of its cell and its room it is at most that, and less
    /// unless the surface is at the hotter one; further below zero, it is more.
    fn loss_slope_at(&self, t: f64, cell: (usize, usize, usize), c: usize) -> f64 {
        let (conductivity, _, emissivity) = self.derived[c];
        let dx = self.dx.to_si();
        let at = Temperature::from_si(t);
        self.air_conductance(cell)
            + self
                .exposed
                .iter()
                .filter(|(face, _)| self.exposed_at(**face, cell))
                .map(|(face, env)| {
                    if emissivity == 0.0 {
                        return self.face_loss(face, env, at, conductivity, emissivity);
                    }
                    let share = env.area.to_si() / self.cells_on(*face).max(1) as f64;
                    let half = 2.0 * conductivity * share / dx;
                    let behind = if !half.is_finite() {
                        self.outward(face, f64::INFINITY, share)
                    } else if half > 0.0 {
                        self.outward(face, half, share)
                    } else {
                        0.0
                    };
                    radiating_film_slope(
                        behind,
                        share,
                        env.convection_w_per_m2_k,
                        emissivity * STEFAN_BOLTZMANN.to_si(),
                        t,
                        env.ambient.to_si(),
                    )
                })
                .sum::<f64>()
    }

    /// One exposed face's part of [`loss_conductance_at`](Solid3D::loss_conductance_at), in W/K:
    /// its film linearised at `at` or at its room, whichever is hotter, with the mounting and the
    /// half cell of `conductivity` in series.
    fn face_loss(
        &self,
        face: &Face,
        env: &Environment,
        at: Temperature,
        conductivity: f64,
        emissivity: f64,
    ) -> f64 {
        let share = env.area.to_si() / self.cells_on(*face).max(1) as f64;
        let hot = at.to_si().max(env.ambient.to_si());
        let h =
            env.convection_w_per_m2_k + 4.0 * emissivity * STEFAN_BOLTZMANN.to_si() * hot.powi(3);
        // The mounting, between the film and the half cell — through `outward`, which
        // `film_flux` also calls, because the limit and the flux disagreeing about a
        // boundary is a defect with no symptom on either side alone.
        series_with_half_cell(
            self.outward(face, h * share, share),
            conductivity,
            share,
            self.dx.to_si(),
        )
    }

    /// How many cells along each axis.
    pub fn counts(&self) -> (usize, usize, usize) {
        self.counts
    }

    /// The cell side.
    pub fn spacing(&self) -> Length {
        self.dx
    }

    /// The block's extent, which is `counts × dx` — the outer faces, not the cell centres.
    pub fn size(&self) -> LengthVec {
        let (nx, ny, nz) = self.counts;
        LengthVec::from_si(DVec3::new(nx as f64, ny as f64, nz as f64) * self.dx.to_si())
    }

    /// The flat index of a cell, or `None` if any component is out of range.
    ///
    /// Returned rather than panicking because the natural way to write a stencil is to ask for a
    /// neighbour that may not exist, and a boundary is exactly where that happens.
    pub fn index(&self, i: usize, j: usize, k: usize) -> Option<usize> {
        let (nx, ny, nz) = self.counts;
        (i < nx && j < ny && k < nz).then(|| i + nx * (j + ny * k))
    }

    /// Where the centre of cell `(i, j, k)` is, in the block's own coordinates.
    pub fn centre_of(&self, i: usize, j: usize, k: usize) -> LengthVec {
        LengthVec::from_si(
            DVec3::new(i as f64 + 0.5, j as f64 + 0.5, k as f64 + 0.5) * self.dx.to_si(),
        )
    }

    /// The temperature of one cell. Out of range reads the nearest one in range.
    pub fn temperature_at(&self, i: usize, j: usize, k: usize) -> Temperature {
        let (nx, ny, nz) = self.counts;
        let idx = self
            .index(i.min(nx - 1), j.min(ny - 1), k.min(nz - 1))
            .expect("clamped indices are in range");
        // **Not a number where there is nothing.** A void cell has no temperature, and a zero
        // or an ambient here is a value somebody would plot, average or believe. Every reader in
        // this workspace that draws a field already skips a non-finite sample.
        if self.void[idx] {
            return Temperature::from_si(f64::NAN);
        }
        Temperature::from_si(self.cells[idx])
    }

    /// Set one cell, for an initial condition a constructor cannot express.
    ///
    /// **This does not change what the block has absorbed.** It is a statement about the initial
    /// state, not a delivery of heat, so `stored_heat` moves and `absorbed_energy` does not — and
    /// a simulation started this way and then audited will show the difference as its opening
    /// balance rather than as a leak. Use [`deposit`](Solid3D::deposit) for heat that arrived.
    ///
    /// Out of range is ignored rather than a panic: a caller writing a hot spot in a loop over a
    /// radius is the expected use, and clipping at the edge is what they mean.
    ///
    /// A temperature that names a different phase in a block with a liquid phase re-derives that
    /// cell's conductivity, capacity and faces, as
    /// [`set_melted_fraction`](Solid3D::set_melted_fraction) does and for the same reason, and at
    /// the cost of one cell: a region set one cell at a time costs what its cells do, not what the
    /// block does.
    pub fn set_temperature(&mut self, i: usize, j: usize, k: usize, t: Temperature) {
        if let Some(idx) = self.index(i, j, k) {
            self.cells[idx] = t.to_si();
            self.cells_moved();
            // A temperature is not a state for a substance that melts, so naming one names a
            // phase as well — see `phase_named_by`.
            if self.latent[idx] > 0.0 {
                let before = self.melted[idx];
                self.melted[idx] = phase_named_by(t.to_si(), self.melt_point[idx]);
                self.rephased(idx, before);
            }
        }
    }

    /// Put joules into one cell, as heat that arrived there.
    ///
    /// Counts toward [`absorbed_energy`](Solid3D::absorbed_energy), so the books balance. Out of
    /// range is ignored, which would silently lose energy — so it does not: the joules are
    /// refused, and nothing is added to either total.
    ///
    /// Heat that melts or freezes some of a cell in a block with a liquid phase re-derives the
    /// cell, as [`set_melted_fraction`](Solid3D::set_melted_fraction) does. It did not: in a
    /// column of eight 1 mm cells of ice at −5 °C, a cell given enough to melt it to water at 5 °C
    /// kept ice's 2.22 mW/K on each of its faces where ice against water is 0.94 mW/K, and a step
    /// of the column's limit, 0.423 s, took the cell to 0.10 °C where the same step through the
    /// right faces took it to 2.92 °C.
    pub fn deposit(&mut self, i: usize, j: usize, k: usize, joules: Energy) {
        let Some(idx) = self.index(i, j, k) else {
            return;
        };
        let before = self.melted[idx];
        let rise = joules.to_si() / self.capacity[idx];
        self.add_kelvin(idx, rise);
        self.absorbed += joules.to_si();
        self.cells_moved();
        self.rephased(idx, before);
    }

    /// Mean over every cell. Every cell has the same volume, so this is the volume average.
    pub fn mean_temperature(&self) -> Temperature {
        // Over what is there. Averaging in the void's placeholder would drag the number toward
        // a value nothing holds, and the more of the box is empty the further it would drag.
        let (sum, n) = (0..self.cells.len())
            .filter(|c| !self.void[*c])
            .fold((0.0, 0usize), |(s, n), c| (s + self.cells[c], n + 1));
        Temperature::from_si(if n == 0 { f64::NAN } else { sum / n as f64 })
    }

    /// The hottest cell — the number a hot spot exists to produce, and the one a lumped model
    /// reports as the mean.
    pub fn peak_temperature(&self) -> Temperature {
        Temperature::from_si(
            (0..self.cells.len())
                .filter(|c| !self.void[*c])
                .map(|c| self.cells[c])
                .fold(f64::MIN, f64::max),
        )
    }

    /// The coldest cell.
    pub fn coldest_temperature(&self) -> Temperature {
        Temperature::from_si(
            (0..self.cells.len())
                .filter(|c| !self.void[*c])
                .map(|c| self.cells[c])
                .fold(f64::MAX, f64::min),
        )
    }

    /// Heat taken from the bus over the run.
    pub fn absorbed_energy(&self) -> Energy {
        Energy::from_si(self.absorbed)
    }

    /// `α·dt/dx²` for the constructor's substance — the classic Fourier number.
    ///
    /// This is the quantity [`mode_amplification`](Solid3D::mode_amplification) is written in terms
    /// of, so it stays the textbook one and does not become shape- or fill-aware. For a block of
    /// one material with at least two cells on every axis it *is* the stability number, and
    /// `fourier_number(max_stable_dt) == 1/6` exactly.
    ///
    /// It is **not** the stability number for a thin block or a filled one, and those are the two
    /// cases where a caller sizing a step by hand needs the other one — see
    /// [`stability_ratio`](Solid3D::stability_ratio), which is what `step` actually enforces.
    pub fn fourier_number(&self, dt: Time) -> f64 {
        let Some(alpha) = self.materials[0].diffusivity() else {
            return f64::INFINITY;
        };
        alpha.to_si() * dt.to_si() / (self.dx.to_si() * self.dx.to_si())
    }

    /// `dt` as a fraction of the largest step this block is stable at — exactly one at the limit.
    ///
    /// The number `step` refuses on, and the honest one for a block that is thin or filled. It is
    /// `dt · maxᵢ (ΣGᵢ / Cᵢ)`, with `ΣGᵢ` everything that leaves cell `i`: a maximum over
    /// **cells**, because stability is a statement about a row of the update matrix and each cell
    /// has its own.
    pub fn stability_ratio(&self, dt: Time) -> f64 {
        dt.to_si() * self.worst_rate_now()
    }

    /// Each cell's own heat capacity, in J/K, in index order.
    ///
    /// For an accelerator that must divide joules by the capacity of the cell they landed in: a
    /// block of two materials has two, and a void has none. [`Solid3D::heat_capacity`] is their
    /// total and cannot answer that.
    pub fn cell_capacities(&self) -> &[f64] {
        &self.capacity
    }

    /// The heat capacity of the whole block, which for a filled one is a sum and not a product.
    pub fn heat_capacity(&self) -> HeatCapacity {
        HeatCapacity::from_si(self.capacity.iter().sum())
    }

    /// The exact per-step amplification of one separable cosine mode.
    ///
    /// `(a, b, c)` are half-wave counts along the three axes: `(1, 0, 0)` is the longest mode
    /// along x with the other two flat. Returns the factor the mode's amplitude is multiplied by
    /// in one step of `dt`, which is `1 + F·Σ(−4 sin²(mπ/2n))` — **exact**, not a linearisation,
    /// because that mode is an eigenvector of the discrete operator this domain steps with.
    ///
    /// Public because it is what makes this domain checkable without a reference implementation:
    /// a caller can predict an amplitude arbitrarily far ahead and compare. It is also what a
    /// grid designer wants, since a factor outside `(−1, 1]` is the instability itself.
    pub fn mode_amplification(&self, mode: (usize, usize, usize), dt: Time) -> f64 {
        let f = self.fourier_number(dt);
        1.0 + f * self.mode_eigenvalue_dx2(mode)
    }

    /// The discrete Laplacian eigenvalue of a mode, times `dx²` — dimensionless, in `[-12, 0]`.
    fn mode_eigenvalue_dx2(&self, mode: (usize, usize, usize)) -> f64 {
        let (nx, ny, nz) = self.counts;
        let term = |m: usize, n: usize| {
            if n <= 1 {
                // One cell across an axis is a mirror against itself: no gradient is
                // representable, so that direction contributes nothing at all.
                return 0.0;
            }
            let s = (m as f64 * std::f64::consts::PI / (2.0 * n as f64)).sin();
            -4.0 * s * s
        };
        term(mode.0, nx) + term(mode.1, ny) + term(mode.2, nz)
    }

    /// Fill the block with one separable cosine mode about a mean.
    ///
    /// The initial condition the closed-form tests use, and a genuinely useful one for anybody
    /// checking a grid: it is the only shape whose future this domain can state exactly.
    pub fn release_mode(&mut self, mode: (usize, usize, usize), mean: Temperature, amplitude: f64) {
        let (nx, ny, nz) = self.counts;
        let phase = |m: usize, n: usize, i: usize| {
            if n <= 1 {
                1.0
            } else {
                (m as f64 * std::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos()
            }
        };
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = i + nx * (j + ny * k);
                    self.cells[idx] = mean.to_si()
                        + amplitude
                            * phase(mode.0, nx, i)
                            * phase(mode.1, ny, j)
                            * phase(mode.2, nz, k);
                }
            }
        }
        // **Each cell's phase named by its temperature**, as `set_temperature` names one. It kept
        // the phase it had: water at 300 K released to 50 K stayed liquid there, an enthalpy below
        // its solid's at absolute zero, and one cell with nothing to meet came out of a step at
        // −19.37 K, its books unmoved — a step that moved no heat, which no limit bounds.
        let mut moved = false;
        for c in 0..self.cells.len() {
            if self.latent[c] > 0.0 {
                let named = phase_named_by(self.cells[c], self.melt_point[c]);
                moved |= named != self.melted[c];
                self.melted[c] = named;
            }
        }
        if moved && self.two_phase {
            self.resolve();
        }
        self.saved.clone_from(&self.cells);
        self.saved_melted.clone_from(&self.melted);
        self.cells_moved();
    }

    /// The amplitude of one mode currently present, by projection.
    ///
    /// The counterpart to [`release_mode`](Solid3D::release_mode), and what makes a decay
    /// measurable rather than merely visible. Cosine modes on this grid are orthogonal, so this
    /// is exact for a block holding one of them and is the correct coefficient for a block
    /// holding several.
    pub fn mode_amplitude(&self, mode: (usize, usize, usize)) -> f64 {
        let (nx, ny, nz) = self.counts;
        let phase = |m: usize, n: usize, i: usize| {
            if n <= 1 {
                1.0
            } else {
                (m as f64 * std::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos()
            }
        };
        // ⟨cos²⟩ is 1/2 per axis that actually varies, and 1 for an axis that cannot.
        let norm = [(mode.0, nx), (mode.1, ny), (mode.2, nz)]
            .iter()
            .map(|&(m, n)| if n <= 1 || m == 0 { 1.0 } else { 0.5 })
            .product::<f64>();
        let mut sum = 0.0;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = i + nx * (j + ny * k);
                    sum += self.cells[idx]
                        * phase(mode.0, nx, i)
                        * phase(mode.1, ny, j)
                        * phase(mode.2, nz, k);
                }
            }
        }
        sum / (self.cells.len() as f64 * norm)
    }

    /// The volume of the whole block.
    pub fn volume(&self) -> Volume {
        let (nx, ny, nz) = self.counts;
        Volume::from_si((nx * ny * nz) as f64 * self.cell_volume())
    }

    fn cell_volume(&self) -> f64 {
        let dx = self.dx.to_si();
        dx * dx * dx
    }

    /// Heat held, measured from the state the block started in.
    ///
    /// Weighted cell by cell, so a filled block's books balance for the same reason a uniform
    /// one's do: the sweep moves `G_f·ΔT` across a face and takes it off one side and puts it on
    /// the other, and this is the sum that is therefore constant.
    ///
    /// **In enthalpy, not in temperature**, so a melting front is on the books. A cell holding at its
    /// melting point while it absorbs 306 mJ per cubic millimetre has taken that heat in and its
    /// temperature says nothing about it; an audit reading temperature alone would call it a leak.
    fn stored_heat(&self) -> f64 {
        (0..self.cells.len())
            .map(|c| self.enthalpy_joules(c) - self.reference_enthalpy_joules(c))
            .sum()
    }

    /// The cell's state as one number, in kelvin above its melting point.
    ///
    /// Monotone and invertible, which is the whole reason the phase change needs no iteration and
    /// conserves exactly:
    ///
    /// ```text
    ///   solid    e = T − T_m ≤ 0
    ///   mush     e = φ·ℓ ∈ [0, ℓ]      at T = T_m
    ///   liquid   e = ℓ + T − T_m ≥ ℓ
    /// ```
    ///
    /// This is the enthalpy method, bookkept as a temperature and a fraction rather than as a single
    /// enthalpy field — identical arithmetic, and it keeps `cells` holding kelvin so everything that
    /// reads a temperature still can.
    ///
    /// It is **not** the apparent-heat-capacity method, and the difference is the failure that method
    /// has: smearing `L` over a narrow temperature interval lets a cell cross the whole interval in
    /// one step and skip the latent heat, which runs the front fast and conserves nothing. Here a
    /// step that overshoots the mush deposits the remainder as sensible heat on the far side, because
    /// the inverse map says where the energy goes rather than a branch guessing.
    fn enthalpy_joules(&self, c: usize) -> f64 {
        let (t, point, latent) = (self.cells[c], self.melt_point[c], self.latent[c]);
        if latent <= 0.0 {
            // No phase change. Measured from the block's own reference rather than from a melting
            // point at infinity, and that choice is precision: subtracting a nearby number keeps
            // digits that subtracting 273.15 from 293.15 does not.
            return self.capacity[c] * (t - self.reference);
        }
        if self.melted[c] >= 1.0 {
            latent + self.cap_liquid[c] * (t - point)
        } else if self.melted[c] <= 0.0 {
            self.cap_solid[c] * (t - point)
        } else {
            self.melted[c] * latent
        }
    }

    /// What [`enthalpy_kelvin`](Solid3D::enthalpy_kelvin) was when the block started, so the ledger
    /// reports what arrived rather than what is there.
    fn reference_enthalpy_joules(&self, c: usize) -> f64 {
        let (point, latent) = (self.melt_point[c], self.latent[c]);
        if latent <= 0.0 {
            return 0.0;
        }
        if self.reference > point {
            latent + self.cap_liquid[c] * (self.reference - point)
        } else {
            self.cap_solid[c] * (self.reference - point)
        }
    }

    /// Add `rise` kelvin of *enthalpy* to a cell, which is not the same as adding kelvin to it.
    ///
    /// The one place heat becomes state, so that a phase change cannot be forgotten at one of the
    /// several doors heat comes in through — the sweep, [`deposit`](Solid3D::deposit), and the plain
    /// channel all pass through here. A cell at its melting point takes the whole of it as melting
    /// and does not warm at all.
    /// The net flux into cell `c`, in kelvin per second of mobility — conduction across its six
    /// faces, less whatever the surface film takes.
    ///
    /// Factored out because the two apply paths differ and the arithmetic must not: a melting cell
    /// goes through the enthalpy and a plain one is a double-buffered write, and a second copy of
    /// a seven-point stencil is a second place for it to be wrong.
    ///
    /// **The film is part of this flux, not a pass after it.** Applying it separately is Lie
    /// splitting, and the split's error carries a coefficient that grows as `1/dx` while the step
    /// falls as `dx²` — so the product falls as `dx`, and the boundary came out **first order**
    /// while the interior was second. Measured before it moved in here: ratios 1.27, 1.71, 1.87
    /// per grid doubling, approaching two rather than four.
    ///
    /// Read off `old` like every other term, so the sweep stays a function of the state it began
    /// with.
    #[inline]
    fn flux_at(&self, old: &[f64], cell: (usize, usize, usize), c: usize) -> f64 {
        let flux = self.conducted(old, cell, c);
        if self.exposed.is_empty() && self.air.is_none() {
            flux
        } else {
            flux - self.film_flux(old, cell, c)
        }
    }

    /// `flux_at` inside a step, with the film read from what the books' pass found for it, `film`,
    /// rather than evaluated again. The same number to the bit: that pass read the same `old` and
    /// the same operator, and a cell `sheds` does not name has no film, so `flux − 0` is what
    /// `flux_at` would have returned and `flux` is that.
    ///
    /// **Every cell used to ask**, through `flux_at`, and asking walks the map of exposed faces to
    /// find that an interior cell is on none of them; the boundary was asked a second time by the
    /// books' pass. Measured, pinned to one core, a 32³ block of aluminium with a film on every
    /// face stepped in 1.35 ms and steps in 0.68 ms, every bit the same.
    #[inline]
    fn swept_at(&self, old: &[f64], cell: (usize, usize, usize), c: usize, film: &[f64]) -> f64 {
        let flux = self.conducted(old, cell, c);
        if (self.exposed.is_empty() && self.air.is_none()) || !self.sheds(cell) {
            flux
        } else {
            flux - film[c]
        }
    }

    /// Whether a cell can have a film: on the block's own boundary, or anywhere when the void holds
    /// air, whose faces are inside. The books' pass and the sweep both ask this, so they agree on
    /// which cells the film buffer holds.
    #[inline]
    fn sheds(&self, (i, j, k): (usize, usize, usize)) -> bool {
        let (nx, ny, nz) = self.counts;
        self.air.is_some()
            || i == 0
            || i + 1 == nx
            || j == 0
            || j + 1 == ny
            || k == 0
            || k + 1 == nz
    }

    /// What cell `c` conducts in across its six faces, in the sweep's units.
    #[inline]
    fn conducted(&self, old: &[f64], cell: (usize, usize, usize), c: usize) -> f64 {
        let (nx, ny, nz) = self.counts;
        let (i, j, k) = cell;
        let t = old[c];
        let (xr, yr) = ((nx + 1) * (j + ny * k), nx * (j + (ny + 1) * k));
        let mut flux = 0.0;
        if i > 0 {
            flux += self.kx[i + xr] * (old[c - 1] - t);
        }
        if i + 1 < nx {
            flux += self.kx[i + 1 + xr] * (old[c + 1] - t);
        }
        if j > 0 {
            flux += self.ky[i + yr] * (old[c - nx] - t);
        }
        if j + 1 < ny {
            flux += self.ky[i + nx + yr] * (old[c + nx] - t);
        }
        if k > 0 {
            flux += self.kz[c] * (old[c - nx * ny] - t);
        }
        if k + 1 < nz {
            flux += self.kz[c + nx * ny] * (old[c + nx * ny] - t);
        }
        flux
    }

    fn add_kelvin(&mut self, c: usize, rise: f64) {
        if self.latent[c] > 0.0 {
            // `rise` is what the cell *would* have warmed by; the joules it stands for are that times
            // whichever capacity the cell has now, which is the mixed one and is exactly right for a
            // cell that is wholly one phase. A mushy cell does not change temperature, so the flux
            // that reaches it has to be converted through the capacity it had when the flux was
            // computed — the same one.
            let joules = rise * self.capacity[c];
            let e = self.enthalpy_joules(c) + joules;
            let (t, phi) = self.state_from_enthalpy(c, e);
            self.cells[c] = t;
            self.melted[c] = phi;
        } else {
            self.cells[c] += rise;
        }
    }

    /// The inverse: a state in kelvin back to a temperature and a melted fraction.
    fn state_from_enthalpy(&self, c: usize, e: f64) -> (f64, f64) {
        let (point, latent) = (self.melt_point[c], self.latent[c]);
        if e <= 0.0 {
            (point + e / self.cap_solid[c], 0.0)
        } else if e >= latent {
            (point + (e - latent) / self.cap_liquid[c], 1.0)
        } else {
            (point, e / latent)
        }
    }
}

impl Domain for Solid3D {
    fn books_balance(&self) -> bool {
        true
    }

    fn name(&self) -> &str {
        &self.name
    }

    /// `minᵢ Cᵢ / ΣGᵢ`, `ΣGᵢ` everything that leaves cell `i` — the tightest row of the update
    /// matrix, and `dx²/(6α)` when the block is one material with at least two cells on every axis
    /// and nothing leaves it but through its faces.
    ///
    /// # Why a maximum over cells, and what summing the faces buys
    ///
    /// Stability is Gershgorin's condition on one row: the update is `I + dt·D⁻¹L`, row `i` has
    /// diagonal `−θᵢ` and off-diagonals summing to `+θᵢ` for `θᵢ = dt Σ_f G_f / Cᵢ`, so every
    /// eigenvalue is in `[1 − 2θ_max, 1]` and `θ_max ≤ 1` is the limit. The sum is over everything
    /// that leaves the cell — its faces, its film, its air and each clearance it radiates across,
    /// those linearised at the hotter end — and it was not: a clearance pair was charged alone, the
    /// largest of a cell's terms where the row is their sum, and three cells in a row across two
    /// clearances were given twice the step their middle row allows. See `worst_rate_now`.
    ///
    /// On a uniform grid that disc is **tight** rather than cautious — the sharpest representable
    /// mode saturates it — and this is where `1/2`, `1/4` and `1/6` come from. They are not three
    /// dimensions; they are however many axes have more than one cell, and the block is stepped at
    /// whichever applies. This used to report `dx²/(6α)` for every shape, which cost a bar-shaped
    /// block three times the steps for nothing.
    ///
    /// On a filled grid the value is the opposite one: the limit is usually far **looser** than
    /// `dx²/(6·α_max)`, because `k_f ≤ 2·min(k_L, k_R)` means a cell cannot be heated through a
    /// face faster than its worse side allows. Heat does not reach a fast material at that
    /// material's own rate; it reaches it at the rate the neighbour delivers. Measured on one
    /// aluminium cell embedded in borosilicate, the honest limit is **75×** the one aluminium's
    /// diffusivity would name, and that factor is wall-clock.
    ///
    /// It can in principle go the other way — up to a factor of two, when a cell's neighbours
    /// conduct better *and* store more — but that needs volumetric heat capacity to vary as widely
    /// as conductivity, and across solids it varies by one order of magnitude where conductivity
    /// varies by four. So the tightening is a bound this catalogue cannot reach, and the loosening
    /// is a saving any coating or inclusion gets.
    ///
    /// # A third of the bar's, for a cube
    ///
    /// A third of what [`Bar1D`](crate::Bar1D) reports for the same spacing and material. That
    /// factor is the reason `Schedule::Multirate` exists: a block and a lumped mass in one scene
    /// differ by five orders of magnitude in the step they can take, and a single global step
    /// would make the cheap domain pay the expensive one's bill.
    ///
    /// # Across a melting point, what the step carries
    ///
    /// A row charged at the capacity a cell has, the slope of its enthalpy, bounds a step within one
    /// phase. A step that carries the cell across its melting point moves its heat at that slope and
    /// gives it back at the other phase's, and where the phase it leaves holds more — water twice
    /// what ice does — it lands past everything it exchanges with. Measured at this limit before it
    /// said otherwise, each step `Ok`: one cell of water under films that only convect, to a room at
    /// 77 K, came to 35.71 K, 41 K below anything it met, and to a room at 0 K, −121.37 K; a solid
    /// holding twice what its liquid does, melting between cells of it at 1000 K, came to 1675 K.
    ///
    /// So a cell a step can carry past the coldest thing it meets, `T_lo`, is held to
    /// `dt·(−F) ≤ H(T) − H(T_lo)` as well, `F` the watts the step moves into it, and at the limit
    /// lands on `T_lo`; and the same to the hottest, the other way. A step within its row can do
    /// that only where `T_lo` is more than `L/(C_L − C_S)` below the melting point — for water,
    /// 116.70 K — or the hottest more than `L/(C_S − C_L)` above it, so anywhere else the limit is
    /// the rows' to the bit, and a block whose bounds nothing above absolute zero reaches is not
    /// asked: octadecane's is −630 K. `F`, and not the row, because a radiating film or clearance
    /// carries less than the tangent its row is charged at, from 900 K to 300 K 0.37 of it, and a
    /// cell that cannot cross is not held as one that can.
    ///
    /// # Stable is not accurate, and this is only the first
    ///
    /// At **exactly** this step the sharpest mode the grid can hold has an amplification factor of
    /// `−1`. It flips sign every step and never decays. That is what marginal stability means and
    /// it is not a defect — the scheme does not diverge there, which is the whole of what a
    /// stability limit claims.
    ///
    /// It does mean sharp initial data is carried badly. A point source excites that mode as hard
    /// as anything can, and the peak comes out **1.96×** the exact answer while the conservation
    /// audit stays exact to the last bit. At half this step it is 1.005×.
    ///
    /// So take a fraction of it when the initial condition is sharp. `Schedule::Multirate` divides
    /// by `ceil(dt / limit)`, or one more where that is an ulp over, and so usually lands
    /// comfortably inside, but a caller stepping by hand can sit exactly on it. `cargo run --example heat_in_three_dimensions` is the
    /// measurement.
    fn max_stable_dt(&self, _now: Time) -> Time {
        let rate = self.worst_rate_now();
        if rate.is_nan() {
            return Time::from_si(f64::INFINITY);
        }
        Time::from_si(1.0 / rate)
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        // Before the limit, because most of what this refuses leaves the limit looking fine: an
        // insulator's is infinite, and a negative capacity's is too.
        if let Some(fault) = self.input_fault() {
            return Err(fault);
        }
        // The step and the heat waiting for it, before anything moves: the heat goes in after the
        // sweep, and a refusal there would leave a block swept with its heat still on the bus.
        if let Some((what, value)) = crate::step_fault(dt).or_else(|| crate::heat_fault(bus)) {
            return Err(Violation::at(&self.name, what, value));
        }
        // The cells a caller set since a step last read them, before the limit reads them too: a
        // block started at −10 K, or given −50 K by `set_temperature`, stepped on from there.
        if self.unchecked {
            if let Some(fault) = self.cell_fault() {
                return Err(fault);
            }
            self.unchecked = false;
        }
        let ratio = self.stability_ratio(dt);
        if ratio.is_nan() {
            return Err(self.unsteppable());
        }
        // Reported in Fourier-number units, which is `dt/(6·max_stable_dt)` — the classic
        // `α·dt/dx²` for a block of one material with two cells on every axis, and the number the
        // limit is expressed in for every other block.
        if ratio > 1.0 + 6e-12 {
            return Err(Violation {
                quantity: "Fourier number".to_string(),
                site: format!("{} (explicit 3D conduction)", self.name),
                before: STABLE_FOURIER_3D,
                after: ratio * STABLE_FOURIER_3D,
                scale: STABLE_FOURIER_3D,
                tolerance: 1e-12,
            });
        }

        // The seven-point stencil in conductance form: `Cᵢ ΔTᵢ = dt Σ_f G_f (T_f − Tᵢ)`.
        //
        // It is the **face flux** that is computed, and each face is read once from each side with
        // opposite sign, so `Σ Cᵢ Tᵢ` is conserved to the last bit rather than to a tolerance —
        // and it stays that way when the capacities differ, which a stencil written as
        // `f·(Σ T_n − 6T)` cannot do because there is no per-cell `f` in it.
        let (nx, ny, nz) = self.counts;
        let cells = nx * ny * nz;
        let dts = dt.to_si();

        // This step's share of the bus, looked at, and taken once the step stands; and whether
        // anything can take a cell below absolute zero — heat taken out — for which the books and
        // the phases the step began with are kept, to put back if it does not stand. A block that
        // takes nothing out pays a compare.
        //
        // The cells move from here, so the range `resolve` found of them no longer holds; a
        // two-phase step finds it again as it ends, and one refused puts back the cells it had.
        self.cells_range = None;
        let share = bus.peek_share(HEAT, dt);
        let taking_out = share < 0.0 || !self.sinks.is_empty();
        let books = (self.lost, self.supplied, self.absorbed);
        if taking_out && self.latent_anywhere {
            self.undo_melted.clone_from(&self.melted);
        }

        // The state the sweep reads, into a buffer that outlives the step.
        let mut old = std::mem::take(&mut self.old);
        old.clear();
        old.extend_from_slice(&self.cells);

        // **What the film took out, totalled over the exposed faces.**
        //
        // A sum, so it is not in the parallel pass below: floating-point addition is not
        // associative and a total accumulated per thread would depend on the core count. Walked in
        // index order over the cells that can shed, which is a *surface* — at 96³ that is 6% of
        // the volume. Read off `old` and before anything is applied, and **kept**, so the sweep
        // subtracts the very number the books add rather than a second evaluation of it.
        let mut film = std::mem::take(&mut self.film);
        if !self.exposed.is_empty() || self.air.is_some() {
            // Sized only for a block that has a film, which is the only kind that reads it.
            if film.len() != cells {
                film.clear();
                film.resize(cells, 0.0);
            }
            for k in 0..nz {
                for j in 0..ny {
                    for i in 0..nx {
                        if !self.sheds((i, j, k)) {
                            continue;
                        }
                        let c = i + nx * (j + ny * k);
                        let shed = self.film_flux(&old, (i, j, k), c);
                        film[c] = shed;
                        if shed != 0.0 {
                            self.lost += dts * self.mobility[c] * shed * self.capacity[c];
                        }
                    }
                }
            }
        }

        if self.latent_anywhere {
            // A melting block applies its own rise cell by cell: `add_kelvin` converts through the
            // enthalpy and writes both a temperature and a melted fraction, and it reads the
            // cell's current state to do it. In place and sequential, which is what it was — and a
            // two-phase step already pays 2.6x to 5.2x for the `resolve` below, so this is not the
            // grid anyone is waiting on.
            for k in 0..nz {
                for j in 0..ny {
                    for i in 0..nx {
                        let c = i + nx * (j + ny * k);
                        let rise =
                            dts * self.mobility[c] * self.swept_at(&old, (i, j, k), c, &film);
                        self.add_kelvin(c, rise);
                    }
                }
            }
        } else {
            // **Double-buffered, and that is where the speed is.** The sweep reads `old` and writes
            // a fresh array which then *becomes* `cells`, so the step moves one read and one write
            // of the grid — where updating in place after a `clone` moved two of each.
            //
            // Measured: the stencil is bandwidth-bound at about 10 GB/s effective on this machine,
            // so halving the traffic is worth more than any number of cores. An earlier version
            // wrote the rise into a third array and read it back, and was **slower than the
            // sequential original** at every grid below 96³ for exactly that reason.
            //
            // **Split across cores, and bit-for-bit the same answer.** Each cell's new value is a
            // function of `old` and of read-only coefficients and of nothing another cell writes,
            // so contiguous chunks perform exactly the operations they performed sequentially, in
            // exactly the same order — see `pantometry_core::sweep`. Chunks are whole **planes**,
            // so the index arithmetic inside is the arithmetic it always was.
            let mut next = std::mem::take(&mut self.rise);
            next.clear();
            next.resize(cells, 0.0);
            {
                let this = &*self;
                let old = &old;
                let film = &film;
                pantometry_core::sweep::fill(&mut next, nx * ny, move |first, chunk| {
                    // The cell's coordinates carried forward rather than divided out: three
                    // integer divisions a cell cost 1.8x on the sequential path when this was
                    // first written that way.
                    let (mut i, mut j, mut k) = (first % nx, (first / nx) % ny, first / (nx * ny));
                    for (at, slot) in chunk.iter_mut().enumerate() {
                        let c = first + at;
                        *slot = old[c]
                            + dts * this.mobility[c] * this.swept_at(old, (i, j, k), c, film);
                        i += 1;
                        if i == nx {
                            i = 0;
                            j += 1;
                            if j == ny {
                                j = 0;
                                k += 1;
                            }
                        }
                    }
                });
            }
            // The new grid becomes the grid, and the one it replaced becomes next step's read
            // buffer. No copy either way.
            std::mem::swap(&mut self.cells, &mut next);
            self.rise = next;
        }

        // **Heat that does have a place**, added after the sweep rather than before it.
        //
        // The distinction is not cosmetic and this learned it by measuring. A forward step is
        // `T′ = T + (dt/C)(P + F(T))`: the source and the conduction are both evaluated at the
        // state the step *began* with. The first version applied the source to `cells` before
        // `old` was taken, so the stencil saw the generated joules and conducted a share of them
        // away inside the same step — a share equal to `G·dt/C`, which for an explicit sweep at
        // its own stability limit is about a half. **The steady state then depended on the
        // timestep**, which is the failure that matters: a bar generating at one end and cooled at
        // the other dropped `P·dx/(kA)` across every gap except the first, and exactly half of it
        // across that one.
        //
        // Splitting it costs nothing in order, unlike the film and the gap exchange above. Those
        // are fluxes proportional to `T`, so applying them separately is Lie splitting with a real
        // error; a constant source commutes with everything and `dT/dt = S` is solved exactly by
        // adding `S·dt`.
        if self.supplies() {
            let dts = dt.to_si();
            for c in 0..self.cells.len() {
                if self.source[c] == 0.0 || self.void[c] {
                    continue;
                }
                let joules = self.source[c] * dts;
                self.add_kelvin(c, joules / self.capacity[c]);
                self.supplied += joules;
            }
        }

        // **Heat off the plain channel**, which carries an amount and no location — added after
        // the sweep, as the source above is and for its reason.
        //
        // It went in first, before `old` was taken, so the film met the heated state and shed a
        // share of the heat in the step that brought it, `θ = dt·G/C` of it. Two things followed,
        // and the books balanced through both, because the heat really did leave:
        //
        // - **The steady state depended on the step.** Measured on one cell convecting at
        //   `h = 100` under constant bus power: it settled at 0.900, 0.500 and 0.100 of `P/G` above
        //   the air at a tenth, a half and nine tenths of its limit, and at its limit it stayed at
        //   the air.
        // - **The step outran its limit.** The limit is taken from the state the step begins in, and
        //   a radiating film met at a temperature the heat had raised is stiffer than that. A black
        //   cell of 1 mm and `k = 200` at 1000 K, stepped at its own limit while the bus raised it
        //   600 K, landed at −33.5 K, and raised 1000 K at −1982 K, every step accepted. This had
        //   been caught once and refused, by a re-check in the film's own pass after the sweep; the
        //   re-check went when the film moved into the sweep, and `a_cold_cell_and_a_bus_slug.rs`,
        //   which allowed a refusal or a step that left no cell below the air, passed on its second
        //   branch.
        //
        // After the sweep, the sweep reads the state the limit was taken from, and the heat is a
        // constant source over the step, which `dT/dt = S` adds exactly.
        //
        // **Spread evenly**, and this is the one place where the 3D domain must not copy the 1D
        // one. `Bar1D` puts placeless heat in its first cell, and that is defensible for a bar,
        // which has an end that a surface absorbing light would plausibly be. A block has six
        // faces and no distinguished cell, so choosing one would invent a location the bus never
        // carried — and a hot spot that came from a tie-break is worse than no hot spot, because
        // it looks like physics. Even spreading is the unique choice that adds no information.
        // Heat that *does* have a place arrives through `deposit` or over an `Interface`.
        // Spread to a uniform **rise**, which for a filled block means in proportion to each cell's
        // capacity rather than in equal joules. Equal joules would warm the low-capacity material
        // more and so would say where the heat landed, which the bus never carried.
        let gained = share;
        if gained != 0.0 {
            self.absorbed += gained;
            let rise = gained / self.capacity.iter().sum::<f64>();
            for c in 0..self.cells.len() {
                // Void takes no share: it has no capacity to contribute to the sum above and no
                // temperature to raise. Warming it would put joules where there is nothing to
                // hold them, and the ledger would then disagree with the block.
                if self.void[c] {
                    continue;
                }
                self.add_kelvin(c, rise);
            }
        }

        // **What the gaps carry.** A pair of solid cells facing each other across void exchange
        // radiation, and it is applied here — in the same explicit update as the conduction
        // sweep, read off the same `old` state — for the reason the surface film moved into the
        // sweep: a separate pass is Lie splitting, whose error carries a coefficient growing as
        // `1/dx` against a step falling as `dx²`, and the boundary comes out an order worse than
        // the interior.
        //
        // **Between the two surfaces, each behind its half cell** — see `crossing`. It was
        // `σA(T_a⁴ − T_b⁴)/(1/ε_a + 1/ε_b − 1)` read at the two centres, each a half cell's drop
        // from the surface that radiates: first order, 1.99 a halving measured on two slabs across
        // a vacuum and 2.01 with air between them.
        //
        // **Antisymmetric where nothing else is at the surfaces**, so the pair conserves exactly:
        // what leaves one arrives in the other, in the same statement, and `Σ Cᵢ Tᵢ` is unchanged
        // to the last bit as it is for a conduction face. With air in the clearance each surface
        // also loses to it, and what one cell gives up is not what the other takes: the air has
        // the difference, and the books hear it as they hear the film.
        if !self.gaps.is_empty() {
            let dts = dt.to_si();
            for (a, b, coefficient) in self.gaps.clone() {
                let pair = self.crossing_at(&old, a, b, coefficient);
                if pair.watts == 0.0 {
                    continue;
                }
                let leaving = pair.watts * pair.from_a * dts;
                let arriving = pair.watts * pair.to_b * dts;
                self.add_kelvin(a, -leaving / self.capacity[a]);
                self.add_kelvin(b, arriving / self.capacity[b]);
                if leaving != arriving {
                    self.lost += leaving - arriving;
                }
            }
        }

        // **Heat taken out that takes a cell below absolute zero is refused**, and the block, its
        // phases, its books and the bus are as the step found them. A block whose sink outran what
        // its film could bring in stepped twenty thousand times to −354.67 K and answered `Ok`
        // each time. Asked only of the cells heat was taken out of — the sinks, and every cell
        // when the bus takes heat out — because within its limit everything else in a step leaves
        // a cell between the coldest and the hottest of what it exchanges with, none below zero:
        // conduction, a film whose room is not, a clearance, whose rows the limit sums, and a
        // cell crossing its melting point, which the limit holds to what it meets — see
        // `max_stable_dt`.
        //
        // That last was not so, and a two-phase block was asked of every cell. A liquid that holds
        // more than its solid, as water does ice, came out of a step at its limit past what it
        // exchanged with, by up to the ratio: water at 1 °C between ice at 1 K came to −108.4 K,
        // and one cell of it under films to a room at 77 K to 36.6 K. Since the limit bounds it,
        // the pass found nothing, and it cost a two-phase step 1.1% to 3.3%.
        if taking_out {
            let below = |c: usize| !self.void[c] && self.cells[c] < 0.0;
            let found = if share < 0.0 {
                (0..self.cells.len()).find(|&c| below(c))
            } else {
                self.sinks.iter().copied().find(|&c| below(c))
            };
            if let Some(c) = found {
                let would = self.cells[c];
                self.cells.copy_from_slice(&old);
                if self.latent_anywhere {
                    self.melted.copy_from_slice(&self.undo_melted);
                }
                (self.lost, self.supplied, self.absorbed) = books;
                self.old = old;
                self.film = film;
                let (i, j, k) = (c % nx, (c / nx) % ny, c / (nx * ny));
                return Err(Violation::at(
                    &self.name,
                    format!(
                        "the heat taken out in this step would take cell ({i}, {j}, {k}) below \
                         absolute zero, to this many kelvin"
                    ),
                    would,
                ));
            }
        }
        // **And held to it, on every debug run of every test that steps a two-phase block**: no
        // cell comes out of a step the limit allowed colder than the coldest thing in or around
        // the block — a cell, a room or the air — unless heat was taken out of it, nor hotter than
        // the hottest unless heat was put in; one that does is a limit that does not bound a phase
        // change, which is what this was. Past either only by a microkelvin and what the step's
        // `1 + 6e-12` of its limit allows: that much more enthalpy moved, over the capacity it is
        // given back at, `6e-12·(C_L/C_S)` of the span at most, and `1e-10` of the largest thing
        // the step read covers a phase holding up to sixteen times the other's.
        if cfg!(debug_assertions) && self.two_phase {
            let around = self
                .exposed
                .values()
                .map(|env| env.ambient.to_si())
                .chain(self.air.map(|(ambient, _)| ambient));
            let (lo, hi, largest) = (0..old.len())
                .filter(|c| !self.void[*c])
                .map(|c| old[c])
                .chain(around)
                .fold(
                    (f64::INFINITY, f64::NEG_INFINITY, 0.0_f64),
                    |(lo, hi, l), t| (lo.min(t), hi.max(t), l.max(t.abs())),
                );
            let slack = crate::ZERO_SLACK + 1e-10 * largest;
            let past = |c: usize| {
                let t = self.cells[c];
                let colder = share >= 0.0 && self.source[c] >= 0.0 && t < lo - slack;
                let hotter = share <= 0.0 && self.source[c] <= 0.0 && t > hi + slack;
                !self.void[c] && (colder || hotter)
            };
            if let Some(c) = (0..self.cells.len()).find(|c| past(*c)) {
                panic!(
                    "{}: a step within the limit took cell {c} from {} K to {} K, outside {lo} K \
                     to {hi} K, the coldest and the hottest thing in or around the block",
                    self.name, old[c], self.cells[c]
                );
            }
        }
        let taken = bus.take_share(HEAT, dt);
        debug_assert_eq!(
            taken.to_bits(),
            share.to_bits(),
            "the share looked at was not taken"
        );

        // **A two-phase block re-derives its conductivities, because they moved.**
        //
        // A cell's `k` and `c` depend on how much of it has melted, so a front that advanced changed
        // the operator. `resolve` carries the fractions the sweep just produced rather than rebuilding
        // them from temperature, and recomputes the faces, the capacities and the limit from them.
        //
        // Only when there is a liquid phase to mix toward: `latent > 0` alone is the one-phase model,
        // whose properties do not depend on the fraction at all, and paying for a rebuild there would
        // slow every non-melting block for nothing.
        //
        // Measured, so the cost is on the record rather than assumed: a `resolve` is 1.6x a `step` at
        // 40 cells and 4.2x at 4096, so a two-phase sweep runs 2.6x to 5.2x a one-phase one. That is
        // the price of the simple version. An incremental update touching only the mush — one or two
        // cells wide — is the optimisation available if a problem ever needs it, and none does yet.
        // The film is applied inside the sweep above, in the same explicit update as
        // conduction — see . It used to run here, as a pass afterwards, and that
        // split cost the boundary an order.
        if self.two_phase {
            self.resolve();
        }
        // The read buffer goes back, so the next step reuses the allocation instead of making one,
        // and the film's with it.
        self.old = old;
        self.film = film;
        Ok(())
    }

    /// Heat gained since the start. The faces are insulated, so this is exactly what came in.
    fn ledger(&self) -> Ledger {
        // `stored + lost`, so the total moves only by what crossed the bus and the claim in
        // `books_balance` survives a block that sheds heat to air. `LumpedMass` carries the
        // same pair for the same reason; what differs is that this one keeps the claim,
        // because every joule it sheds is counted here in the same step it leaves a cell.
        // `stored + lost − supplied`. A source is energy entering from outside the domain, so it
        // is subtracted here for the same reason `lost` is added: the ledger is what the *bus*
        // moved, and a joule this block generated for itself never crossed it. Without the term a
        // dissipating block's books grow by its own output every step and the audit stops the run.
        // Three contributions rather than their sum, and the difference is the whole reason
        // `Ledger` records a scale. `add` raises that scale to the largest single entry, and the
        // audit judges a change against it — which is what makes a relative tolerance mean
        // anything when the net is near zero.
        //
        // Adding them here first threw that away. A block that starts at its own reference
        // temperature stores nothing, so a scene stating a source opened its books at **exactly
        // zero**, and the first 3.7e-11 J of rounding was judged a 100% change: the audit stopped
        // a correct run on its first step. The sum is the same; the scale is now the size of the
        // numbers the rounding actually lives on.
        Ledger::new()
            .with(quantity::ENERGY, self.stored_heat())
            .with(quantity::ENERGY, self.lost)
            .with(quantity::ENERGY, -self.supplied)
    }

    /// The temperatures **and** the phase, because a temperature alone is not a state.
    ///
    /// A cell at 0 °C is ice, water or any mixture, so saving `cells` and not `melted` would restore
    /// a block that was half melted as one that was entirely solid at the same temperature — losing
    /// 306 mJ per cubic millimetre with nothing to say it had gone. `Schedule::Iterative` and the
    /// audit's retry both restore, so this is a live path and not a precaution.
    fn checkpoint(&mut self) {
        self.saved.clone_from(&self.cells);
        self.saved_melted.clone_from(&self.melted);
        // The losses too, and `LumpedMass` records what it costs to forget: an iterative sweep
        // that rewinds the cells and not the counter grows its books by one sweep of shed heat
        // per iteration, and the audit reports energy created from nothing.
        self.saved_lost = self.lost;
        self.saved_supplied = self.supplied;
        self.saved_absorbed = self.absorbed;
    }

    fn restore(&mut self) {
        self.cells.clone_from(&self.saved);
        self.melted.clone_from(&self.saved_melted);
        // What was saved need not be what the last step read: a cell set below absolute zero,
        // saved, put right and stepped comes back from the checkpoint below it.
        self.cells_moved();
        self.lost = self.saved_lost;
        self.supplied = self.saved_supplied;
        self.absorbed = self.saved_absorbed;
        // A two-phase block's step ends by re-deriving the capacities, the conductances and the
        // limit from the fractions it leaves, so the ones a sweep thrown away left are re-derived
        // from the fractions just put back. Without it the next sweep began on the thrown-away
        // sweep's operator: ice at its melting point, fed and restored to a checkpoint before its
        // first step, reported a limit of 0.439 s where it had reported 0.423 s.
        if self.two_phase {
            self.resolve();
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    /// Peak, mean and coldest, in celsius, and what it has absorbed.
    ///
    /// All three ends of the distribution, because the whole reason to pay for a 3D grid is that
    /// they differ. A block reported by its mean alone is a `LumpedMass` that cost `n³` times as
    /// much, and the gap between peak and mean is the number that says whether the reduction
    /// would have been honest.
    fn readings(&self) -> Vec<Reading> {
        let mut out = vec![
            Reading::new(
                &self.name,
                "peak",
                self.peak_temperature().to_si() - 273.15,
                "C",
            ),
            Reading::new(
                &self.name,
                "mean",
                self.mean_temperature().to_si() - 273.15,
                "C",
            ),
            Reading::new(
                &self.name,
                "coldest",
                self.coldest_temperature().to_si() - 273.15,
                "C",
            ),
            Reading::new(&self.name, "absorbed", self.absorbed_energy().to_si(), "J"),
        ];
        // What the block made for itself, and **only** for a block that makes any. A source
        // nobody can see in the report is the shape this workspace calls a silent failure: the
        // run would be right and the reader would have no way to tell 45 W from 45 mW. Conditional
        // for the same reason `melted` is — a column of zeros in every other report costs a line
        // in all of them and tells nobody anything — and the condition is fixed at construction,
        // so it is not a mode that can surprise somebody mid-run.
        if self.supplies() {
            out.push(Reading::new(
                &self.name,
                "generated",
                self.generated_energy().to_si(),
                "J",
            ));
        }
        // The melted volume, and **only** for a block that can melt. A column of zeros in every
        // report tells a reader nothing and costs a line in all of them; the condition is a property
        // of the block's materials fixed at construction, so it is not a mode that can surprise
        // anybody mid-run. Without it a phase change is the one thing this domain does that a report
        // could not see, and the temperature would say a cell was holding still.
        if self.latent.iter().any(|l| *l > 0.0) {
            out.push(Reading::new(
                &self.name,
                "melted",
                self.melted_volume().to_si() * 1e9,
                "mm3",
            ));
        }
        out
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }

    /// A temperature field, so nothing above has to know this is a block.
    fn as_field(&self) -> Option<&dyn ScalarField> {
        Some(self)
    }
}

impl ScalarField for Solid3D {
    /// **Centred** — the values are cell averages and none sits on a boundary.
    ///
    /// Without this a caller sampling on this field's own grid gets its cell *boundaries*, and
    /// the peak of anything sharp comes out low. See [`Lattice`].
    fn lattice(&self) -> Lattice {
        Lattice::Centred
    }

    /// **Kelvin**, because that is what the cells hold. See [`Bar1D`](crate::Bar1D).
    fn unit(&self) -> &'static str {
        "K"
    }

    /// Trilinear between cell centres, clamped at the faces, and **masked over void**.
    ///
    /// Clamped rather than extrapolated: outside an insulated face the temperature is not
    /// defined, and continuing the gradient would draw a block hotter than any cell in it.
    ///
    /// Masked because `self.cells` still holds whatever an emptied cell held when it was emptied,
    /// and reading it raw is how a gap came out of every exporter as a piece of the block sitting
    /// at ambient forever — a plausible, unchanging number, which is worse than no number. The
    /// weights are renormalised over the solid corners, so a sample at a cell centre is that cell
    /// exactly, a sample inside a clearance is `NaN`, and a sample straddling the two is the
    /// material's own value rather than a blend with something that is not there.
    fn at(&self, p: LengthVec, t: Time) -> f64 {
        self.sample_with(p, t, &self.cells)
    }

    /// Central differences on the cell grid, mirrored at the faces **and at void**.
    ///
    /// A clearance is a boundary, and the block already knows what to do at one: mirror. Walking
    /// off the outer edge and walking into nothing are the same situation — there is no sample
    /// that way — so both return the asking cell and the difference becomes one-sided, which is
    /// the insulated condition rather than a slope towards a number that is not there.
    ///
    /// `NaN` inside a clearance, for the same reason [`temperature_at`](Solid3D::temperature_at)
    /// gives one: a gradient there is a value somebody would plot.
    fn gradient(&self, p: LengthVec, _t: Time, _h: Length) -> DVec3 {
        let (i, j, k) = self.nearest_cell(p);
        if self.absent(i, j, k) {
            return DVec3::splat(f64::NAN);
        }
        let d = 2.0 * self.dx.to_si();
        let m = |a: isize, b: isize, c: isize| self.mirrored_into(i, j, k, a, b, c);
        DVec3::new(
            (m(i + 1, j, k) - m(i - 1, j, k)) / d,
            (m(i, j + 1, k) - m(i, j - 1, k)) / d,
            (m(i, j, k + 1) - m(i, j, k - 1)) / d,
        )
    }

    /// `∇²T` on the seven-point stencil.
    ///
    /// The Laplacian of the *temperature*, which is what the trait asks for and is a statement
    /// about the field rather than about the material. It is the operator the sweep uses only when
    /// the block is one material; for a filled one the sweep uses `∇·(k∇T)`, and
    /// [`rate`](ScalarField::rate) is where that appears.
    fn laplacian(&self, p: LengthVec, _t: Time, _h: Length) -> f64 {
        let (i, j, k) = self.nearest_cell(p);
        if self.absent(i, j, k) {
            return f64::NAN;
        }
        let dx = self.dx.to_si();
        let centre = self.mirrored(i, j, k);
        // Each neighbour that is not there mirrors back to the centre and contributes nothing to
        // `sum - 6·centre`, which is exactly what an insulated face does — and exactly what a face
        // touching a clearance carries, since the harmonic face mean is already zero there.
        let m = |a: isize, b: isize, c: isize| self.mirrored_into(i, j, k, a, b, c);
        let sum = m(i - 1, j, k)
            + m(i + 1, j, k)
            + m(i, j - 1, k)
            + m(i, j + 1, k)
            + m(i, j, k - 1)
            + m(i, j, k + 1);
        (sum - 6.0 * centre) / (dx * dx)
    }

    /// `∂T/∂t = (1/Cᵢ)·Σ_f G_f (T_f − Tᵢ)`, which is `α∇²T` where the block is one material.
    ///
    /// The conductance form rather than `α·laplacian`, so that it is the sweep's own operator at a
    /// point in a filled block too — read off the same face conductivities, not reconstructed from
    /// a diffusivity the block may not have a single value of.
    ///
    /// Conduction only: heat arriving over the bus is a source this cannot see.
    fn rate(&self, p: LengthVec, _t: Time, _dt: Time) -> f64 {
        let (nx, ny, nz) = self.counts;
        let (i, j, k) = self.nearest_cell(p);
        let clamp = |v: isize, n: usize| v.clamp(0, n as isize - 1) as usize;
        let (i, j, k) = (clamp(i, nx), clamp(j, ny), clamp(k, nz));
        let c = i + nx * (j + ny * k);
        // Nothing does not warm at a rate, and the `is_finite` guard below would otherwise report
        // that it warms at **zero** — `mobility` is `1/C`, a cell with no capacity has none, and
        // `inf * 0.0` is the `NaN` that guard was written to catch coming from somewhere else.
        // A confident zero here is indistinguishable from a solid cell in equilibrium.
        if self.void[c] {
            return f64::NAN;
        }
        let t = self.cells[c];
        let (xr, yr) = ((nx + 1) * (j + ny * k), nx * (j + (ny + 1) * k));
        let mut flux = 0.0;
        if i > 0 {
            flux += self.kx[i + xr] * (self.cells[c - 1] - t);
        }
        if i + 1 < nx {
            flux += self.kx[i + 1 + xr] * (self.cells[c + 1] - t);
        }
        if j > 0 {
            flux += self.ky[i + yr] * (self.cells[c - nx] - t);
        }
        if j + 1 < ny {
            flux += self.ky[i + nx + yr] * (self.cells[c + nx] - t);
        }
        if k > 0 {
            flux += self.kz[c] * (self.cells[c - nx * ny] - t);
        }
        if k + 1 < nz {
            flux += self.kz[c + nx * ny] * (self.cells[c + nx * ny] - t);
        }
        let rate = self.mobility[c] * flux;
        if rate.is_finite() {
            rate
        } else {
            0.0
        }
    }
}

impl Solid3D {
    /// Whether the cell a stencil is asking about is not there — off the grid, or empty.
    ///
    /// The two cases are one case for a derivative: neither is a sample. Signed, so a stencil can
    /// ask about a neighbour it has already walked off the edge to reach.
    fn absent(&self, i: isize, j: isize, k: isize) -> bool {
        let (nx, ny, nz) = self.counts;
        let outside = |v: isize, n: usize| v < 0 || v >= n as isize;
        if outside(i, nx) || outside(j, ny) || outside(k, nz) {
            return true;
        }
        self.void[i as usize + nx * (j as usize + ny * k as usize)]
    }

    /// [`mirrored`](Solid3D::mirrored), but a neighbour that is **not there** mirrors back to the
    /// cell that asked rather than to the edge of the grid.
    ///
    /// The distinction only shows up inside a block: `mirrored` clamps, so a stencil at cell 3
    /// reaching into a clearance at cell 4 would be handed cell 4's frozen value. This hands back
    /// cell 3, which is the one-sided difference a boundary deserves.
    fn mirrored_into(&self, ci: isize, cj: isize, ck: isize, i: isize, j: isize, k: isize) -> f64 {
        if self.absent(i, j, k) {
            self.mirrored(ci, cj, ck)
        } else {
            self.mirrored(i, j, k)
        }
    }

    /// The value at a cell, with out-of-range indices mirrored back — an insulated face.
    ///
    /// A mirror rather than a zero, and the distinction is the whole boundary condition: a zero
    /// neighbour is a face held at absolute zero and would drain the block, where a mirror is a
    /// face with no gradient across it and so no flow through it.
    fn mirrored(&self, i: isize, j: isize, k: isize) -> f64 {
        let (nx, ny, nz) = self.counts;
        let clamp = |v: isize, n: usize| v.clamp(0, n as isize - 1) as usize;
        let idx = clamp(i, nx) + nx * (clamp(j, ny) + ny * clamp(k, nz));
        self.cells[idx]
    }

    /// The cell a point falls in, as signed indices so a stencil can walk off the edge.
    fn nearest_cell(&self, p: LengthVec) -> (isize, isize, isize) {
        let q = p.to_si() / self.dx.to_si();
        let one = |v: f64| {
            if v.is_nan() {
                0
            } else {
                v.floor().clamp(-1.0, 1e9) as isize
            }
        };
        (one(q.x), one(q.y), one(q.z))
    }
}

impl Solid3D {
    /// Sample this block's **geometry** with somebody else's values.
    ///
    /// The trilinear interpolation, the clamping at the faces and — the part that is easy to get
    /// wrong — the void handling, applied to a slice of kelvin the caller supplies instead of to
    /// `self.cells`.
    ///
    /// # Why this is public
    ///
    /// `pantometry-gpu` holds a block's values on a device and mirrors them back in `f32`. It has to
    /// answer [`ScalarField::at`] for the same grid, and the alternative is a second copy of the
    /// arithmetic above — where "a sample belongs to the cell it is in, and if that cell is empty
    /// there is nothing there to sample" is a rule that took a measurement to get right: 36 of 72
    /// empty cells in a two-layer gap came out solid before it existed. A second copy is a second
    /// chance at that.
    ///
    /// Read-only on purpose. A setter that let a caller overwrite a block's cells would let it put
    /// one in a state its own conservation bookkeeping never agreed to.
    ///
    /// `kelvin` shorter than the grid answers `NaN` beyond its end rather than panicking, because
    /// the caller is a mirror of a device and a short mirror is a bug worth *seeing* as a hole.
    pub fn sample_with(&self, p: LengthVec, _t: Time, kelvin: &[f64]) -> f64 {
        let (nx, ny, nz) = self.counts;
        let q = p.to_si() / self.dx.to_si() - DVec3::splat(0.5);
        // NaN spelled out rather than folded into a comparison: a visualiser can hand one over,
        // and it must not reach the cast below. Answered with a `NaN` rather than with cell zero —
        // a question about nowhere has no answer, and cell zero may itself be empty, in which case
        // the old fallback returned a frozen number for a point that was never asked about.
        if q.is_nan() {
            return f64::NAN;
        }
        let axis = |v: f64, n: usize| -> (usize, f64) {
            let last = n.saturating_sub(1);
            if v <= 0.0 {
                return (0, 0.0);
            }
            if v >= last as f64 {
                return (last, 0.0);
            }
            let i = v.floor();
            (i as usize, v - i)
        };
        let (i, fx) = axis(q.x, nx);
        let (j, fy) = axis(q.y, ny);
        let (k, fz) = axis(q.z, nz);

        // **A sample belongs to the cell it is in**, and if that cell is empty there is nothing
        // there to sample. Without this the masked weights below still answer — with the value of
        // whichever solid neighbour the sample leans towards — and a panel that samples across the
        // extent rather than on cell centres lands most of a clearance's first layer inside the
        // material's half. Measured on a two-layer gap: 36 of the 72 empty cells came out solid.
        let nearest = |v: f64, n: usize| (v.round().max(0.0) as usize).min(n.saturating_sub(1));
        if self.void[nearest(q.x, nx) + nx * (nearest(q.y, ny) + ny * nearest(q.z, nz))] {
            return f64::NAN;
        }
        let (i1, j1, k1) = (
            (i + 1).min(nx - 1),
            (j + 1).min(ny - 1),
            (k + 1).min(nz - 1),
        );

        let mut sum = 0.0;
        let mut weight = 0.0;
        for (a, wa) in [(i, 1.0 - fx), (i1, fx)] {
            for (b, wb) in [(j, 1.0 - fy), (j1, fy)] {
                for (c, wc) in [(k, 1.0 - fz), (k1, fz)] {
                    let w = wa * wb * wc;
                    let at = a + nx * (b + ny * c);
                    // A corner with no weight is not a corner, so a sample sitting on a cell
                    // centre never consults the neighbour it does not use — which is what makes
                    // an exactly-sampled grid exact rather than nearly so.
                    if w > 0.0 && !self.void[at] {
                        let Some(v) = kelvin.get(at) else {
                            return f64::NAN;
                        };
                        sum += w * v;
                        weight += w;
                    }
                }
            }
        }
        if weight > 0.0 {
            sum / weight
        } else {
            f64::NAN
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Face, Solid3D};
    use pantometry_core::{
        units::{Length, Temperature},
        Substance,
    };

    /// **A face visits the cells it holds**, each once, on blocks of every shape up to four on a
    /// side, flat and empty ones included. `cells` is what the stability limit visits, and a plane
    /// picked wrong or a face's two axes swapped need not show in a block whose face cells are all
    /// alike: every middle cell sets the same limit there, and a cell read past the end of a row is
    /// the start of the next.
    #[test]
    fn a_face_visits_the_cells_it_holds() {
        for nx in 0..=4 {
            for ny in 0..=4 {
                for nz in 0..=4 {
                    let counts = (nx, ny, nz);
                    for face in Face::ALL {
                        let mut visited: Vec<_> = face.cells(counts).collect();
                        visited.sort();
                        let mut held = Vec::new();
                        for k in 0..nz {
                            for j in 0..ny {
                                for i in 0..nx {
                                    if face.holds((i, j, k), counts) {
                                        held.push((i, j, k));
                                    }
                                }
                            }
                        }
                        held.sort();
                        assert_eq!(visited, held, "{face:?} on {counts:?}");
                    }
                }
            }
        }
    }

    /// **The cells listed as facing the air are the solid ones with an air conductance, each
    /// once, and asking again rebuilds the list rather than adding to it.**
    ///
    /// The stability limit reads this list in place of every cell, so a cell it leaves out is one
    /// whose loss the limit never sees. A direction checked wrong need not show in a limit: a
    /// cavity's stiffest cells are its corners, which face the void three ways at once and between
    /// them all six. So the list is checked cell by cell against `air_conductance`, which the flux
    /// and the limit both read, on a block whose void is scattered so that each of the six
    /// neighbours is the only void beside at least three cells. `resolve` runs every step for a
    /// two-phase block, and a list it did not clear would grow by itself each step.
    #[test]
    fn the_cells_facing_the_air_are_those_with_an_air_conductance() {
        let (nx, ny, nz) = (4, 5, 3);
        let mut block = Solid3D::new(
            "scattered",
            Substance::copper(),
            (nx, ny, nz),
            Length::from_si(1e-3),
            Temperature::celsius(20.0),
        )
        .empty(|i, j, k| (i + 2 * j + k) % 4 == 0)
        .air_in(Temperature::celsius(20.0), 1.0);

        let mut conducting = Vec::new();
        let mut alone = [0usize; 6];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    if block.is_void(i, j, k) {
                        continue;
                    }
                    if block.air_conductance((i, j, k)) > 0.0 {
                        conducting.push((i, j, k));
                    }
                    // Out of range reads as solid, so only the low side needs a guard.
                    let beside = [
                        i > 0 && block.is_void(i - 1, j, k),
                        block.is_void(i + 1, j, k),
                        j > 0 && block.is_void(i, j - 1, k),
                        block.is_void(i, j + 1, k),
                        k > 0 && block.is_void(i, j, k - 1),
                        block.is_void(i, j, k + 1),
                    ];
                    if beside.iter().filter(|b| **b).count() == 1 {
                        if let Some(way) = beside.iter().position(|b| *b) {
                            alone[way] += 1;
                        }
                    }
                }
            }
        }
        assert!(
            alone.iter().all(|n| *n >= 3),
            "each way the only void beside three cells: {alone:?}"
        );
        assert_eq!(
            block.facing_air, conducting,
            "the list against the air's own conductance"
        );

        // Asked again, twice over: a fill that changes no cell, and the air stated again.
        block.fill(Substance::copper(), |_, _, _| false);
        let block = block.air_in(Temperature::celsius(20.0), 1.0);
        assert_eq!(
            block.facing_air.len(),
            conducting.len(),
            "asked again, the list grew"
        );
        assert_eq!(
            block.facing_air, conducting,
            "asked again, the list changed"
        );
    }

    /// One pair, drawn from the range `crossing` is documented for: conductivities from 0.01 to
    /// 500 W/m·K and cells from 0.1 to 30 mm, a decade at a time; emissivities from 0.02 to one;
    /// centres from 200 K to 3000 K, a fifth of them within two kelvin of each other; and, on half
    /// of them, air at 200 K to 1500 K and 1 to 10⁴ W/m²·K. As `(C, (G_a, T_a), (G_b, T_b), air)`.
    type Pair = (f64, (f64, f64), (f64, f64), Option<(f64, f64)>);

    fn a_pair(index: u64) -> Pair {
        let mut rng = pantometry_core::Rng::for_index(0x7a1c_4a11, index);
        let mut decade = |lo: f64, hi: f64| 10f64.powf(rng.range(lo, hi));
        let dx = decade(-4.0, -1.5);
        let (ka, kb) = (decade(-2.0, 2.7), decade(-2.0, 2.7));
        let h = decade(0.0, 4.0);
        let area = dx * dx;
        let (ea, eb) = (rng.range(0.02, 1.0), rng.range(0.02, 1.0));
        let coefficient = super::STEFAN_BOLTZMANN.to_si() * area / (1.0 / ea + 1.0 / eb - 1.0);
        let ta = rng.range(200.0, 3000.0);
        let tb = if rng.unit() < 0.2 {
            ta + rng.range(-2.0, 2.0)
        } else {
            rng.range(200.0, 3000.0)
        };
        let air = (rng.unit() < 0.5).then(|| (h * area, rng.range(200.0, 1500.0)));
        (coefficient, (2.0 * ka * dx, ta), (2.0 * kb * dx, tb), air)
    }

    /// **`crossing` lands on the root of the balance it solves, and in a few steps.**
    ///
    /// The reference halves on `a`'s surface `u`, which `crossing` never solves for: what reaches
    /// it through its half cell, less what its air takes, is what crosses; `b`'s surface is where
    /// that, and what its own air brings, is what its half cell carries on; and the two surfaces
    /// send `C(u⁴ − v⁴)` across, which rises with `u`. Halving finds that to an ulp or so of `u`,
    /// which moves what crosses by `(G_a + H)` times an ulp, and the fourth powers round by an ulp
    /// of `C·u⁴`: `crossing` is held to eight of each, or to 1e-13 of what crosses where that is
    /// more. A Newton whose derivative is wrong still converges on most pairs, only slowly, so the
    /// steps are held as well.
    #[test]
    fn a_crossing_lands_on_its_root_in_a_few_steps() {
        let (mut most, mut taken, mut worst) = (0, 0, 0.0_f64);
        let pairs = 8000;
        for index in 0..pairs {
            let (c, (ga, ta), (gb, tb), air) = a_pair(index);
            let (h, room) = air.unwrap_or((0.0, 0.0));
            let (crossed, steps) = super::crossing_counted(c, (ga, ta), (gb, tb), air);
            most = most.max(steps);
            taken += steps;

            let across = |u: f64| ga * (ta - u) - h * (u - room);
            let far = |q: f64| (gb * tb + h * room + q) / (gb + h);
            let fourth = |t: f64| (t * t) * (t * t);
            let residual = |u: f64| c * (fourth(u) - fourth(far(across(u)))) - across(u);
            // Between where nothing crosses, `a`'s source, and where the two surfaces meet, their
            // sources' mean weighted by what holds each: outside that, `v` can pass through zero
            // and its fourth power turn the residual back.
            let source = |g: f64, t: f64| (g * t + h * room) / (g + h);
            let (sa, sb) = (source(ga, ta), source(gb, tb));
            let meet = ((ga + h) * sa + (gb + h) * sb) / (ga + gb + 2.0 * h);
            let (mut lo, mut hi) = (sa.min(meet), sa.max(meet));
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                if residual(mid) > 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            let u = 0.5 * (lo + hi);
            let reference = across(u);
            let v = far(reference);
            let resolution = 8.0 * f64::EPSILON * ((ga + h) * u + c * (fourth(u) + fourth(v)));
            let allowed = resolution.max(1e-13 * reference.abs());
            let miss = (crossed.watts - reference).abs();
            worst = worst.max(miss / allowed);
            assert!(
                miss <= allowed,
                "pair {index}: crossing gives {:e} W where its balance has its root at {reference:e} \
                 W, {miss:e} apart against {allowed:e}, in {steps} steps",
                crossed.watts
            );
        }
        println!(
            "  {pairs} pairs: {most} steps at most, {:.2} on average, the worst miss {worst:.2} of \
             what it is allowed",
            taken as f64 / pairs as f64
        );
        assert!(
            most <= 12,
            "{most} steps where the documentation says eleven at most"
        );
    }

    /// **The slopes the steady solve leans on are the derivatives of what crosses**: what `a`'s
    /// cell gives up against `a`'s centre, and what `b`'s takes against `b`'s, by central
    /// differences of `crossing` itself a ten-thousandth of the temperature apart. The difference
    /// is good to about `(δ/T)²` and to the rounding of what crosses over `δ`, so 1e-6 of the slope
    /// is a thousand times what either leaves, against slopes left out of a crossing, or charged
    /// by the share `f` rather than `f²`, or at the centres' tangent.
    #[test]
    fn a_crossings_slopes_are_its_derivatives() {
        for index in 0..2000 {
            let (c, (ga, ta), (gb, tb), air) = a_pair(index);
            let (at_a, at_b) = super::crossing(c, (ga, ta), (gb, tb), air).slopes(c);
            let gives = |t: f64| {
                let p = super::crossing(c, (ga, t), (gb, tb), air);
                p.watts * p.from_a
            };
            let takes = |t: f64| {
                let p = super::crossing(c, (ga, ta), (gb, t), air);
                p.watts * p.to_b
            };
            let (da, db) = (1e-4 * ta, 1e-4 * tb);
            let by_a = (gives(ta + da) - gives(ta - da)) / (2.0 * da);
            let by_b = -(takes(tb + db) - takes(tb - db)) / (2.0 * db);
            // A floor where a slope is nothing: a pair whose air takes almost all of a surface.
            let hot = ta.max(tb).max(air.map_or(0.0, |(_, t)| t));
            let floor = 1e-12 * 4.0 * c * hot * hot * hot;
            for (end, slope, by) in [("a", at_a, by_a), ("b", at_b, by_b)] {
                let apart = (slope - by).abs();
                assert!(
                    apart <= 1e-6 * by.abs() + floor,
                    "pair {index}: end {end}'s slope is {slope:e} W/K and what crosses moves by \
                     {by:e} W/K, {apart:e} apart"
                );
            }
        }
    }

    /// **The couplings are the cross derivatives of what crosses**: what `a`'s cell takes in moving
    /// with `b`'s temperature, and `b`'s with `a`'s — the off-diagonal the steady solve's step on
    /// the levels of several bodies is built from — against central differences, over the same
    /// eight thousand pairs and their air. And without air, each is the other end's slope.
    #[test]
    fn a_crossings_couplings_are_its_cross_derivatives() {
        for index in 0..2000 {
            let (c, (ga, ta), (gb, tb), air) = a_pair(index);
            let pair = super::crossing(c, (ga, ta), (gb, tb), air);
            let (by_b, by_a) = pair.couplings(c);
            let a_takes = |t: f64| {
                let p = super::crossing(c, (ga, ta), (gb, t), air);
                -p.watts * p.from_a
            };
            let b_takes = |t: f64| {
                let p = super::crossing(c, (ga, t), (gb, tb), air);
                p.watts * p.to_b
            };
            let (da, db) = (1e-4 * ta, 1e-4 * tb);
            let moved_by_b = (a_takes(tb + db) - a_takes(tb - db)) / (2.0 * db);
            let moved_by_a = (b_takes(ta + da) - b_takes(ta - da)) / (2.0 * da);
            let hot = ta.max(tb).max(air.map_or(0.0, |(_, t)| t));
            let floor = 1e-12 * 4.0 * c * hot * hot * hot;
            for (end, coupling, by) in [("a by b", by_b, moved_by_b), ("b by a", by_a, moved_by_a)]
            {
                let apart = (coupling - by).abs();
                assert!(
                    apart <= 1e-6 * by.abs() + floor,
                    "pair {index}: {end}'s coupling is {coupling:e} W/K and what crosses moves by \
                     {by:e} W/K, {apart:e} apart"
                );
            }
            if air.is_none() {
                let (at_a, at_b) = pair.slopes(c);
                assert_eq!(
                    (by_b.to_bits(), by_a.to_bits()),
                    (at_b.to_bits(), at_a.to_bits()),
                    "pair {index}: without air a coupling is the other end's slope"
                );
            }
        }
    }

    /// `T·|T|³`, what a surface below absolute zero radiates.
    fn odd_fourth(t: f64) -> f64 {
        (t * t) * (t * t) * t.signum()
    }

    /// The watts a pair's crossing can be resolved to, in ulps: an ulp of `a`'s surface `u` moves
    /// what crosses by `G_a + H` times it, the fourth powers round by an ulp of `C·u⁴` and `C·v⁴`,
    /// and **each source rounds by an ulp of its centre and of its air's share**, which is what
    /// `crossing` forms it from, `T − w(T − T_air)`: what crosses moves by at most `G + H` times
    /// the source, and the air's part of the source is `w·T_air`, so `(G + H)·|T| + H·T_air`. The
    /// first is an ulp of `u` where the centre is near its surface, as every pair `a_pair` draws
    /// has it; a centre at −2435 K beside air at 974 K has a source of 54 K, and the ulps of 2435 K
    /// are forty-five of 54 K. Counted as `(G + H)·T_air`, as first written, pairs with both
    /// sources within a kelvin of zero were allowed more than all they carry.
    fn resolved(
        c: f64,
        (ga, ta): (f64, f64),
        (gb, tb): (f64, f64),
        air: Option<(f64, f64)>,
        (u, v): (f64, f64),
    ) -> f64 {
        let (h, room) = air.unwrap_or((0.0, 0.0));
        let fourth = |t: f64| (t * t) * (t * t);
        (ga + h) * (u.abs() + ta.abs())
            + (gb + h) * tb.abs()
            + 2.0 * h * room
            + c * (fourth(u) + fourth(v))
    }

    /// One pair drawn as [`a_pair`] draws one, with its centres moved below absolute zero: `b`'s
    /// alone on a third of them, both on a third, `a`'s alone on the rest, each from 0 to 3000 K
    /// from zero and a fifth of them within two kelvin of it, and the other end at what `a_pair`
    /// gave it. The air, where there is any, is at or above zero, as a block's is, so a source it
    /// warms may not be below zero after all: [`a_crossing_below_zero_lands_on_its_root`] counts
    /// how many are.
    fn a_pair_below_zero(index: u64) -> Pair {
        let (c, (ga, ta), (gb, tb), air) = a_pair(index);
        let mut rng = pantometry_core::Rng::for_index(0x0b31_0e20, index);
        let mut below = || {
            if rng.unit() < 0.2 {
                -rng.range(0.0, 2.0)
            } else {
                -rng.range(0.0, 3000.0)
            }
        };
        let (ta, tb) = match index % 3 {
            0 => (ta, below()),
            1 => (below(), below()),
            _ => (below(), tb),
        };
        (c, (ga, ta), (gb, tb), air)
    }

    /// **Below absolute zero `crossing` lands on the root of the balance it solves, read odd, and in
    /// a few steps.** The same reference as [`a_crossing_lands_on_its_root_in_a_few_steps`] — a
    /// halving on `a`'s surface, which `crossing` never solves for — with each fourth power read as
    /// `T·|T|³`. Read so, the residual rises with `u` everywhere, so the halving needs no care about
    /// where `v` passes through zero. Held to eight ulps of what it is made of, as [`resolved`]
    /// counts them, which here includes each source's own rounding. **And the draws reach what they
    /// are for**: a pair whose colder source is below zero, with its surfaces at the root both
    /// below zero, one on each side, or both lifted to or above it by what crosses — the three
    /// readings `crossing_below_zero` has — counted and held to a share of the draws each, since a
    /// sabotage of one passes on draws that never reach it.
    #[test]
    fn a_crossing_below_zero_lands_on_its_root() {
        let (mut most, mut taken, mut worst) = (0, 0, 0.0_f64);
        let (mut below, mut astride, mut both, mut lifted) = (0, 0, 0, 0);
        let pairs = 8000;
        for index in 0..pairs {
            let (c, (ga, ta), (gb, tb), air) = a_pair_below_zero(index);
            let (h, room) = air.unwrap_or((0.0, 0.0));
            let source = |g: f64, t: f64| (g * t + h * room) / (g + h);
            let (sa, sb) = (source(ga, ta), source(gb, tb));
            let (crossed, steps) = super::crossing_counted(c, (ga, ta), (gb, tb), air);
            let (u, v) = crossed.surfaces;
            if sa.min(sb) < 0.0 {
                below += 1;
                most = most.max(steps);
                taken += steps;
                if u.max(v) < 0.0 {
                    both += 1;
                } else if u.min(v) < 0.0 {
                    astride += 1;
                } else {
                    lifted += 1;
                }
            }

            let across = |u: f64| ga * (ta - u) - h * (u - room);
            let far = |q: f64| (gb * tb + h * room + q) / (gb + h);
            let residual = |u: f64| c * (odd_fourth(u) - odd_fourth(far(across(u)))) - across(u);
            let meet = ((ga + h) * sa + (gb + h) * sb) / (ga + gb + 2.0 * h);
            let (mut lo, mut hi) = (sa.min(meet), sa.max(meet));
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                if residual(mid) > 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            let u = 0.5 * (lo + hi);
            let reference = across(u);
            let v = far(reference);
            // Eight ulps and no floor of a share of what crosses, which no effect here earns.
            let allowed = 8.0 * f64::EPSILON * resolved(c, (ga, ta), (gb, tb), air, (u, v));
            let miss = (crossed.watts - reference).abs();
            worst = worst.max(miss / allowed);
            assert!(
                miss <= allowed,
                "pair {index}: crossing gives {:e} W where its balance has its root at {reference:e} \
                 W, {miss:e} apart against {allowed:e}, in {steps} steps: {:?}, sources {sa} K and \
                 {sb} K, surfaces {u} K and {v} K",
                crossed.watts,
                (c, (ga, ta), (gb, tb), air)
            );
        }
        println!(
            "  {pairs} pairs, {below} with a source below zero: {most} steps at most there, {:.2} on \
             average; of those, {both} with both surfaces below zero, {astride} with one each \
             side and {lifted} with neither; the worst miss {worst:.2} of what it is allowed",
            taken as f64 / below.max(1) as f64
        );
        for (what, count, least) in [
            ("a source below zero", below, pairs / 2),
            ("both surfaces below zero", both, pairs / 10),
            ("a surface each side of zero", astride, pairs / 10),
            (
                "a source below zero and neither surface",
                lifted,
                pairs / 100,
            ),
        ] {
            assert!(count >= least, "{count} pairs of {pairs} with {what}");
        }
        // **The steps by their mean, and by a most with room**: twelve at most here is one draw's
        // tail, and over 10⁵ draws fourteen pairs took thirteen or fourteen; an ulp of an input
        // moves a step count, and these draws go through the platform's `pow`. The mean is what
        // moves with a defect that only slows the root: both below zero read as `v⁴ − u⁴`, its
        // two fourth powers cancelling, lands on it within the allowance and took 5.15 to 5.20
        // on average where this takes 5.02, and 5.04 over 10⁵.
        let mean = taken as f64 / below as f64;
        assert!(
            most <= 16 && mean <= 5.1,
            "{most} steps at most and {mean:.3} on average, where the documentation says fourteen \
             and 5.04 over 10⁵ pairs, held to sixteen and 5.1"
        );
    }

    /// **Below absolute zero the slopes and couplings are still the derivatives of what crosses**:
    /// [`a_crossings_slopes_are_its_derivatives`] and
    /// [`a_crossings_couplings_are_its_cross_derivatives`] on pairs below zero, the derivative of
    /// `T·|T|³` being `4|T|³`. A slope read as `4T³` there is negative, and the steady solve would
    /// step a cell the wrong way.
    #[test]
    fn below_zero_a_crossings_slopes_and_couplings_are_its_derivatives() {
        let (mut worst, mut floored, mut below, mut floor_share) = (0.0_f64, 0, 0, 0.0_f64);
        for index in 0..2000 {
            let (c, (ga, ta), (gb, tb), air) = a_pair_below_zero(index);
            let pair = super::crossing(c, (ga, ta), (gb, tb), air);
            if pair.surfaces.0.min(pair.surfaces.1) < 0.0 {
                below += 1;
            }
            let (at_a, at_b) = pair.slopes(c);
            let (by_b, by_a) = pair.couplings(c);
            let at = |t_a: f64, t_b: f64| super::crossing(c, (ga, t_a), (gb, t_b), air);
            // **A step that moves neither surface by more than a ten-thousandth of itself**, since
            // the fourth powers are of the surfaces and the difference is good to the square of
            // that. A ten-thousandth of the centre was not: a half cell of 1e-4 W/K holding a
            // surface at 49 K behind a centre at 2131 K was differenced across 0.43% of its
            // surface, and a step of `b`'s that moved `a`'s, at −30.87 K behind 7.6e-4 W/K, by
            // 3.4% of it missed by 2.3e-6. Found from a trial step, what crossing moves being
            // linear in it at this size, and relative to each surface however near zero: held to
            // a kelvin's ten-thousandth there, two surfaces at −0.40 K and −0.094 K missed by
            // 1.1e-6. A microkelvin is the floor, where a slope is `4C·1e-18` W/K.
            let (u, v) = pair.surfaces;
            let step = |moved: &dyn Fn(f64) -> (f64, f64), t: f64| {
                let trial = 1e-6 * t.abs().max(1.0);
                let (u1, v1) = moved(t + trial);
                let shift =
                    ((u1 - u).abs() / u.abs().max(1e-6)).max((v1 - v).abs() / v.abs().max(1e-6));
                if shift > 0.0 {
                    trial * 1e-4 / shift
                } else {
                    trial
                }
            };
            let da = step(&|t| at(t, tb).surfaces, ta);
            let db = step(&|t| at(ta, t).surfaces, tb);
            let central =
                |f: &dyn Fn(f64) -> f64, t: f64, d: f64| (f(t + d) - f(t - d)) / (2.0 * d);
            let gives = |t: f64| {
                let p = at(t, tb);
                p.watts * p.from_a
            };
            let takes = |t: f64| {
                let p = at(ta, t);
                p.watts * p.to_b
            };
            let a_takes = |t: f64| {
                let p = at(ta, t);
                -p.watts * p.from_a
            };
            let b_takes = |t: f64| {
                let p = at(t, tb);
                p.watts * p.to_b
            };
            // **The rounding of what crosses, over the interval it is differenced across**: of what
            // crosses and its fourth powers, and of the moved centre's source, which moves what
            // crosses by its slope; what rounds alike at both ends of the difference cancels. A
            // surface near zero is moved by a ten-thousandth of a kelvin, and what crosses is set
            // by the other end: a slope of 4.2e-9 W/K was differenced from 1.57 W, and an ulp of
            // that over the step is 1.1e-12 W/K.
            let fourth = |t: f64| (t * t) * (t * t);
            let room = air.map_or(0.0, |(_, t)| t);
            let rounding = |by: f64, t: f64, d: f64| {
                8.0 * f64::EPSILON
                    * (c * (fourth(u) + fourth(v)) + pair.watts.abs() + by.abs() * (t.abs() + room))
                    / d
            };
            let h = air.map_or(0.0, |(h, _)| h);
            for (what, derivative, by, t, d, g) in [
                ("a's slope", at_a, central(&gives, ta, da), ta, da, ga),
                ("b's slope", at_b, -central(&takes, tb, db), tb, db, gb),
                (
                    "a's coupling to b",
                    by_b,
                    central(&a_takes, tb, db),
                    tb,
                    db,
                    gb,
                ),
                (
                    "b's coupling to a",
                    by_a,
                    central(&b_takes, ta, da),
                    ta,
                    da,
                    ga,
                ),
            ] {
                let apart = (derivative - by).abs();
                let floor = rounding(by, t, d);
                let allowed = 1e-6 * by.abs() + floor;
                worst = worst.max(apart / allowed);
                if floor > 1e-6 * by.abs() {
                    floored += 1;
                }
                // **And the floor admits nothing the solve would notice**: a millionth of what
                // holds the moved end's surface, `G + H`, beside which the steady solve's diagonal
                // reads a pair's slope. Measured, 1.9e-8 of it at most.
                floor_share = floor_share.max(floor / (g + h));
                assert!(
                    floor <= 1e-6 * (g + h),
                    "pair {index}: {what}'s floor is {floor:e} W/K, against {:e} W/K holding its \
                     surface",
                    g + h
                );
                assert!(
                    derivative >= 0.0 && apart <= allowed,
                    "pair {index}: {what} is {derivative:e} W/K and what crosses moves by {by:e} \
                     W/K, {apart:e} apart: {:?}, surfaces {:?}",
                    (c, (ga, ta), (gb, tb), air),
                    pair.surfaces
                );
            }
        }
        // The rounding is what is allowed on a few derivatives only, or this holds nothing: 772
        // of them, measured — 483 at a pair with a surface within 10 K of zero, and the other 289
        // each a derivative through the end whose surface is the nearer zero, where `4C|T|³` is
        // next to nothing beside what the other end sends. Held to a tenth.
        println!(
            "  8000 derivatives at 2000 pairs, {below} with a surface below zero: the worst \
             {worst:.3} of what it is allowed, {floored} held to the rounding rather than to 1e-6, \
             the rounding {floor_share:.1e} at most of what holds the surface"
        );
        assert!(
            floored <= 800,
            "{floored} of 8000 derivatives held only to the rounding of what crosses, where 772 \
             were"
        );
        // And on pairs that are below zero: the slopes are those of the surfaces, which a pair
        // whose sources are below zero can lift above it.
        assert!(
            below >= 1000,
            "{below} pairs of 2000 with a surface below zero"
        );
    }

    /// **A film whose cell is below absolute zero lands on the root of its surface's balance, in a
    /// few steps.** A face drawn as a block can have one: a half cell of 0.01 to 500 W/m·K and 0.1
    /// to 30 mm, behind a mounting of 10² to 10⁵ W/m²·K on half of them; convection of none on a
    /// third and 1 to 10⁴ W/m²·K on the rest; an emissivity from 0.02 to one; a room at absolute
    /// zero on a fifth and up to 3000 K on the rest; and the centre 0 to 3000 K below zero, a fifth
    /// within two kelvin of it. The reference halves on the surface between the centre and the room
    /// where what `behind` brings it is what it sheds, `A·h·(T − T∞) + A·εσ(T·|T|³ − T∞⁴)`, which
    /// rises with `T`; the film is held to what that sheds to eight ulps of the terms it is made
    /// of, or 1e-13 of it, and where that surface is below zero, `film_surface_below_zero` to the
    /// twelve steps its bound documents. **And the draws reach what they are for**: a surface below
    /// zero, counted, as well as one lifted to or above it, which the loop above zero finds.
    #[test]
    fn a_film_below_zero_lands_on_its_root() {
        let sigma = super::STEFAN_BOLTZMANN.to_si();
        let (mut most, mut taken, mut worst, mut below) = (0, 0, 0.0_f64, 0);
        let faces = 8000;
        for index in 0..faces {
            let mut rng = pantometry_core::Rng::for_index(0xf11e_0b31, index);
            let mut decade = |lo: f64, hi: f64| 10f64.powf(rng.range(lo, hi));
            let dx = decade(-4.0, -1.5);
            let k = decade(-2.0, 2.7);
            let mount = decade(2.0, 5.0);
            let h = decade(0.0, 4.0);
            let share = dx * dx;
            let half = 2.0 * k * dx;
            let behind = if index % 2 == 0 {
                half
            } else {
                1.0 / (1.0 / half + 1.0 / (mount * share))
            };
            let convection = if index % 3 == 0 { 0.0 } else { h };
            let radiative = rng.range(0.02, 1.0) * sigma;
            // A fifth, and not a tenth, which would be even indices only and so never behind a
            // mounting.
            let ambient = if index % 5 == 0 {
                0.0
            } else {
                rng.range(0.0, 3000.0)
            };
            let centre = if rng.unit() < 0.2 {
                -rng.range(0.0, 2.0)
            } else {
                -rng.range(0.0, 3000.0)
            };
            let film = super::radiating_film(behind, share, convection, radiative, centre, ambient);

            let a4 = (ambient * ambient) * (ambient * ambient);
            let sheds =
                |t: f64| share * (convection * (t - ambient) + radiative * (odd_fourth(t) - a4));
            let residual = |t: f64| sheds(t) - behind * (centre - t);
            let (mut lo, mut hi) = (centre, ambient);
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                if residual(mid) > 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            let t = 0.5 * (lo + hi);
            let reference = sheds(t);
            let mut steps = 0;
            if t < 0.0 {
                below += 1;
                steps = super::film_surface_below_zero(
                    behind, share, convection, radiative, centre, ambient,
                )
                .1;
                most = most.max(steps);
                taken += steps;
            }
            let fourth = (t * t) * (t * t);
            // Eight ulps of the terms and no floor of a share of what crosses, which no effect
            // here earns: with one of 1e-13, it was the allowance on three faces in four.
            let allowed = 8.0
                * f64::EPSILON
                * (behind * (centre.abs() + t.abs())
                    + share * (convection * (t.abs() + ambient) + radiative * (fourth + a4)));
            let flux = film * (centre - ambient);
            let miss = (flux - reference).abs();
            worst = worst.max(miss / allowed);
            assert!(
                miss <= allowed,
                "face {index}: the film carries {flux:e} W where its surface's balance has its root \
                 at {reference:e} W, {miss:e} apart against {allowed:e}, in {steps} steps"
            );
        }
        println!(
            "  {faces} faces, {below} with the surface below zero: {most} steps at most there, \
             {:.2} on average; the worst miss {worst:.2} of what it is allowed",
            taken as f64 / below.max(1) as f64
        );
        // Both ways through: a centre below zero whose surface is lifted to or above it as well.
        for (what, count) in [("below zero", below), ("at or above zero", faces - below)] {
            assert!(
                count >= faces / 4,
                "{count} faces of {faces} with the surface {what}"
            );
        }
        // The steps by their mean, and by a most with room, as for the pair: over 10⁵ faces,
        // fourteen at most and 4.70 on average.
        let mean = taken as f64 / below as f64;
        assert!(
            most <= 16 && mean <= 4.85,
            "{most} steps at most and {mean:.3} on average, where the documentation says fourteen \
             and 4.70 over 10⁵ faces, held to sixteen and 4.85"
        );
    }

    /// **A radiating face's slope is the derivative of what it carries**, `radiating_film(…)·(T_c −
    /// T∞)` against its centre, by central differences, over faces drawn as
    /// [`a_film_below_zero_lands_on_its_root`] draws them but with the centre anywhere from 3000 K
    /// below zero to 3000 K above it, a fifth of them within two kelvin of zero, and the room from
    /// zero to 3000 K. The step moves the surface by a ten-thousandth of itself, since the fourth
    /// power is of the surface, and the difference is held to 1e-6 of the slope or the rounding of
    /// what crosses over the step. **And where it bounds it, the limit's tangent is no less**:
    /// where the surface is no further from zero than the hotter of centre and room, the slope is
    /// at most `behind` in series with `A·(h + 4εσ·max(T_c, T∞)³)`, which is what the steady solve's
    /// sweep reads and the stability limit charges. The draws reach every route: the surface below
    /// zero, at or above it with the centre below, and the centre above zero, each counted.
    #[test]
    fn a_films_slope_is_the_derivative_of_what_it_carries() {
        let sigma = super::STEFAN_BOLTZMANN.to_si();
        let (mut worst, mut floored) = (0.0_f64, 0);
        let (mut below, mut lifted, mut above) = (0, 0, 0);
        let faces = 4000;
        for index in 0..faces {
            let mut rng = pantometry_core::Rng::for_index(0x51_0be5, index);
            let mut decade = |lo: f64, hi: f64| 10f64.powf(rng.range(lo, hi));
            let dx = decade(-4.0, -1.5);
            let k = decade(-2.0, 2.7);
            let mount = decade(2.0, 5.0);
            let h = decade(0.0, 4.0);
            let share = dx * dx;
            let half = 2.0 * k * dx;
            let behind = if index % 2 == 0 {
                half
            } else {
                1.0 / (1.0 / half + 1.0 / (mount * share))
            };
            let convection = if index % 3 == 0 { 0.0 } else { h };
            let radiative = rng.range(0.02, 1.0) * sigma;
            let ambient = if index % 5 == 0 {
                0.0
            } else {
                rng.range(0.0, 3000.0)
            };
            let size = if rng.unit() < 0.2 {
                rng.range(0.0, 2.0)
            } else {
                rng.range(0.0, 3000.0)
            };
            // Drawn, not by index, which would tie the centre's side to the mounting's.
            let centre = if rng.unit() < 0.5 { size } else { -size };
            let (surface, _) =
                super::film_surface(behind, share, convection, radiative, centre, ambient);
            if centre >= 0.0 {
                above += 1;
            } else if surface < 0.0 {
                below += 1;
            } else {
                lifted += 1;
            }

            let carries = |t: f64| {
                super::radiating_film(behind, share, convection, radiative, t, ambient)
                    * (t - ambient)
            };
            let slope =
                super::radiating_film_slope(behind, share, convection, radiative, centre, ambient);
            // The surface moves by at most what the centre does, so a ten-thousandth of the
            // surface, or a microkelvin where it is nearer zero than that.
            let d = 1e-4 * surface.abs().max(1e-2);
            let by = (carries(centre + d) - carries(centre - d)) / (2.0 * d);
            let fourth = (surface * surface) * (surface * surface);
            let a4 = (ambient * ambient) * (ambient * ambient);
            let rounding = 8.0
                * f64::EPSILON
                * (carries(centre).abs()
                    + behind * (centre.abs() + surface.abs())
                    + share * (convection * (surface.abs() + ambient) + radiative * (fourth + a4)))
                / d;
            let apart = (slope - by).abs();
            let allowed = 1e-6 * by.abs() + rounding;
            worst = worst.max(apart / allowed);
            if rounding > 1e-6 * by.abs() {
                floored += 1;
            }
            assert!(
                slope >= 0.0 && apart <= allowed,
                "face {index}: the slope is {slope:e} W/K and what crosses moves by {by:e} W/K, \
                 {apart:e} apart against {allowed:e}: centre {centre} K, room {ambient} K, surface \
                 {surface} K"
            );
            let hot = centre.max(ambient);
            if surface.abs() <= hot {
                let tangent = share * (convection + 4.0 * radiative * hot * hot * hot);
                let bound = behind * tangent / (behind + tangent);
                assert!(
                    slope <= bound * (1.0 + 4.0 * f64::EPSILON),
                    "face {index}: the slope {slope:e} W/K is more than the limit's {bound:e} W/K: \
                     centre {centre} K, room {ambient} K, surface {surface} K"
                );
            }
        }
        println!(
            "  {faces} faces, {below} with the surface below zero, {lifted} lifted to or above it \
             and {above} with the centre above: the worst {worst:.3} of what it is allowed, \
             {floored} held to the rounding"
        );
        for (what, count) in [
            ("the surface below zero", below),
            ("the surface lifted to or above it", lifted),
            ("the centre above zero", above),
        ] {
            assert!(count >= faces / 10, "{count} faces of {faces} with {what}");
        }
        // 263 of them, measured: held to a tenth, or this holds nothing.
        assert!(
            floored <= faces / 10,
            "{floored} of {faces} slopes held only to the rounding, where 263 were"
        );
        // The three faces with no surface to find, which no draw above reaches: nothing behind it,
        // the surface's own tangent at the centre, read odd; a film without limit, the half cell;
        // nothing reaching it, nothing.
        let (share, h, radiative) = (1e-6, 10.0, 0.9 * sigma);
        for centre in [-300.0, 0.5, 900.0] {
            let own = share * (h + 4.0 * radiative * centre * centre * centre.abs());
            let at = super::radiating_film_slope(f64::INFINITY, share, h, radiative, centre, 600.0);
            assert_eq!(at.to_bits(), own.to_bits(), "nothing behind, {centre} K");
            let held =
                super::radiating_film_slope(0.1, share, f64::INFINITY, radiative, centre, 600.0);
            assert_eq!(held, 0.1, "a film without limit, {centre} K");
            let none = super::radiating_film_slope(0.0, share, h, radiative, centre, 600.0);
            assert_eq!(none, 0.0, "nothing behind the surface carries, {centre} K");
        }
        // And nothing reaching a surface at zero in a room at zero, where its own tangent is zero
        // too: nothing over nothing, but for the guard. Above, the tangent is not zero, and the
        // series gives zero with the guard or without it.
        let none = super::radiating_film_slope(0.0, share, 0.0, radiative, 0.0, 0.0);
        assert_eq!(
            none, 0.0,
            "nothing reaches the surface, at zero in a room at zero"
        );
    }

    /// **A cell none of whose faces radiates reads the slope the limit charges, to the bit**, so
    /// that the steady solve of a block with no radiating face takes the steps it took before the
    /// step on a body's level read its films' own slope: on a block of three cells a side, films on
    /// three faces — one behind a mounting — and air in a void on one of them, every solid cell at
    /// five temperatures either side of the room.
    #[test]
    fn a_cell_that_does_not_radiate_reads_the_limits_slope() {
        use pantometry_core::units::Area;
        let film = |h: f64| super::Environment {
            ambient: Temperature::celsius(20.0),
            convection_w_per_m2_k: h,
            area: Area::from_si(9e-6),
        };
        let block = Solid3D::new(
            "block",
            Substance::copper().with_emissivity(0.0),
            (3, 3, 3),
            Length::from_si(1e-3),
            Temperature::celsius(20.0),
        )
        // On a cooled face, so that (0, 1, 0) and (1, 2, 0) have air and two films: three terms,
        // whose order a sum rounds by.
        .empty(|i, j, k| (i, j, k) == (1, 1, 0))
        .air_in(Temperature::celsius(20.0), 25.0)
        .losing_from(Face::XMin, film(10.0))
        .losing_from(Face::YMax, film(1e4))
        .losing_from(Face::ZMin, film(300.0))
        .mounted_on(Face::ZMin, 2e3);
        let mut read = 0;
        for k in 0..3 {
            for j in 0..3 {
                for i in 0..3 {
                    if block.is_void(i, j, k) {
                        continue;
                    }
                    let c = i + 3 * (j + 3 * k);
                    for t in [-50.0, 1.0, 293.15, 400.0, 2000.0] {
                        let limit = block.loss_conductance_at(Temperature::from_si(t), (i, j, k));
                        let slope = block.loss_slope_at(t, (i, j, k), c);
                        assert_eq!(
                            slope.to_bits(),
                            limit.to_bits(),
                            "cell ({i}, {j}, {k}) at {t} K: {slope:e} W/K against {limit:e} W/K"
                        );
                        read += 1;
                    }
                }
            }
        }
        assert_eq!(read, 26 * 5, "every solid cell read");
        // **And it can see an order**: the same terms summed the other way, the air first and
        // then each face, round differently on at least one cell, or this compares nothing a
        // change of order would move. `face_loss` reads no temperature where nothing radiates, so
        // one temperature serves.
        let mut differs = 0;
        for k in 0..3 {
            for j in 0..3 {
                for i in 0..3 {
                    if block.is_void(i, j, k) {
                        continue;
                    }
                    let c = i + 3 * (j + 3 * k);
                    let (conductivity, _, emissivity) = block.derived[c];
                    let at = Temperature::from_si(400.0);
                    let other = block
                        .exposed
                        .iter()
                        .filter(|(face, _)| block.exposed_at(**face, (i, j, k)))
                        .fold(block.air_conductance((i, j, k)), |sum, (face, env)| {
                            sum + block.face_loss(face, env, at, conductivity, emissivity)
                        });
                    let limit = block.loss_conductance_at(at, (i, j, k));
                    if other.to_bits() != limit.to_bits() {
                        differs += 1;
                    }
                }
            }
        }
        assert!(
            differs >= 1,
            "no cell rounds the other order differently, so this cannot see one"
        );
    }

    /// **A cell's slope is the derivative of what its films and air carry**, `film_flux` itself
    /// differenced at the cell, on a block that reaches every route of `loss_slope_at` that a
    /// block which can be solved has, a conductivity without limit being refused by name: faces
    /// that radiate with no convection, with convection, and behind a mounting; a cell of
    /// `k = 0` on a radiating face, whose film carries nothing; air in a void on a cooled face;
    /// corner cells with three films, each from a room of its own; and cells below zero, below
    /// every room, between them and above them. **And
    /// where it bounds it, the limit's tangent is no less**: at or above zero the slope is at most
    /// what `loss_conductance_at` charges. Counted: the mounted radiating faces and the cell of
    /// `k = 0` are reached.
    #[test]
    fn a_cells_slope_is_the_derivative_of_what_its_films_carry() {
        use pantometry_core::substance::ThermalProps;
        use pantometry_core::units::{
            Area, Density, SpecificHeat, ThermalConductivity, ThermalExpansion,
        };
        let stuff = |k: f64| {
            Substance::bulk("stuff", Density::from_si(2000.0)).with_thermal(ThermalProps {
                conductivity: ThermalConductivity::w_per_m_k(k),
                specific_heat: SpecificHeat::j_per_kg_k(900.0),
                expansion: ThermalExpansion::ppm_per_k(10.0),
                emissivity: 0.9,
            })
        };
        // A room of its own for each, so that a face read at another's room is seen.
        let film = |h: f64, room: f64| super::Environment {
            ambient: Temperature::from_si(room),
            convection_w_per_m2_k: h,
            area: Area::from_si(9e-6),
        };
        let mut block = Solid3D::new(
            "block",
            stuff(0.2),
            (3, 3, 3),
            Length::from_si(1e-3),
            Temperature::from_si(300.0),
        )
        .empty(|i, j, k| (i, j, k) == (1, 1, 0))
        .air_in(Temperature::from_si(350.0), 25.0)
        .losing_from(Face::XMin, film(0.0, 250.0))
        .losing_from(Face::YMax, film(50.0, 300.0))
        .losing_from(Face::ZMin, film(10.0, 400.0))
        .mounted_on(Face::ZMin, 2e3);
        block.fill(stuff(0.0), |i, j, k| (i, j, k) == (0, 0, 1));
        let dx = 1e-3;
        let (mut mounted, mut insulated, mut worst) = (0, 0, 0.0_f64);
        for k in 0..3 {
            for j in 0..3 {
                for i in 0..3 {
                    if block.is_void(i, j, k) {
                        continue;
                    }
                    let c = i + 3 * (j + 3 * k);
                    let cell = (i, j, k);
                    if k == 0 && block.derived[c].0 > 0.0 {
                        mounted += 1;
                    }
                    if block.derived[c].0 == 0.0 {
                        insulated += 1;
                    }
                    for t in [-600.0, 150.0, 320.0, 900.0] {
                        let carries = |at: f64| {
                            let mut field = vec![300.0; 27];
                            field[c] = at;
                            block.film_flux(&field, cell, c) * dx
                        };
                        let d = 1e-4 * f64::abs(t);
                        let (above, below) = (carries(t + d), carries(t - d));
                        let by = (above - below) / (2.0 * d);
                        let slope = block.loss_slope_at(t, cell, c);
                        // An ulp or so of what crosses at each end, over the step.
                        let floor = 16.0 * f64::EPSILON * above.abs().max(below.abs()) / d;
                        let apart = (slope - by).abs();
                        let allowed = 1e-6 * by.abs() + floor;
                        worst = worst.max(apart / allowed);
                        assert!(
                            apart <= allowed,
                            "cell {cell:?} at {t} K: the slope is {slope:e} W/K and what its films \
                             carry moves by {by:e} W/K, {apart:e} apart"
                        );
                        if t >= 0.0 {
                            let limit = block.loss_conductance_at(Temperature::from_si(t), cell);
                            assert!(
                                slope <= limit * (1.0 + 4.0 * f64::EPSILON),
                                "cell {cell:?} at {t} K: the slope {slope:e} W/K is more than \
                                 the limit's {limit:e} W/K"
                            );
                        }
                    }
                }
            }
        }
        println!(
            "  26 cells at four temperatures, {mounted} behind the mounting, {insulated} of k = 0: \
             the worst {worst:.3} of what it is allowed"
        );
        assert!(
            mounted >= 5 && insulated == 1,
            "{mounted} mounted, {insulated} insulated"
        );
    }

    /// **Air without limit holds both surfaces at its own temperature**, so nothing crosses and
    /// nothing moves with either centre. `∞/∞` is not a share, and taken as one it made both
    /// surfaces not a number and the steady solve report that it diverged.
    #[test]
    fn air_without_limit_holds_both_surfaces_at_its_temperature() {
        let c = 4.6e-14;
        let pair = super::crossing(
            c,
            (1.2e-3, 900.0),
            (0.1, 600.0),
            Some((f64::INFINITY, 300.0)),
        );
        assert_eq!(pair.watts, 0.0, "{pair:?}");
        assert_eq!((pair.from_a, pair.to_b), (0.0, 0.0), "{pair:?}");
        assert_eq!(pair.slopes(c), (0.0, 0.0), "{pair:?}");
    }

    /// **Nothing behind either surface is the centres' reading**, exactly as it was: half cells
    /// that conduct without limit put each surface at its centre.
    #[test]
    fn nothing_behind_the_surfaces_is_the_centres_reading() {
        let c = 4.6e-14;
        let pair = super::crossing(c, (f64::INFINITY, 900.0), (f64::INFINITY, 600.0), None);
        let fourth = |t: f64| (t * t) * (t * t);
        let centres = c * (fourth(900.0) - fourth(600.0));
        assert!(
            (pair.watts / centres - 1.0).abs() < 1e-15,
            "{} W against the centres' {centres} W",
            pair.watts
        );
        assert_eq!(pair.surfaces, (900.0, 600.0));
    }

    /// **Nothing behind a surface below zero is its centre's reading, odd**: a pair whose half
    /// cells conduct without limit carries `C(g(T_a) − g(T_b))` between its centres, `g(T) =
    /// T·|T|³`, on each side of zero and with both below it; a film with nothing behind its face
    /// carries the odd secant at the centre, `A(h + εσ(T⁴ + T∞⁴)/(T∞ − T))`. No block reaches
    /// either — it refuses a conductivity without limit by name — so this is the only thing that
    /// does.
    #[test]
    fn nothing_behind_a_surface_below_zero_is_its_centres_odd_reading() {
        let c = 4.6e-14;
        for (ta, tb) in [
            (900.0, -600.0),
            (-600.0, 900.0),
            (-300.0, -800.0),
            (0.0, -50.0),
        ] {
            let pair = super::crossing(c, (f64::INFINITY, ta), (f64::INFINITY, tb), None);
            let centres = c * (odd_fourth(ta) - odd_fourth(tb));
            assert!(
                (pair.watts / centres - 1.0).abs() < 1e-15,
                "{ta} K and {tb} K: {} W against the centres' {centres} W",
                pair.watts
            );
            assert_eq!(pair.surfaces, (ta, tb));
        }
        let (share, h, radiative) = (1e-6, 10.0, 0.9 * super::STEFAN_BOLTZMANN.to_si());
        for (centre, room) in [(-300.0, 1000.0), (-1500.0, 300.0), (-0.5, 0.0)] {
            let film = super::radiating_film(f64::INFINITY, share, h, radiative, centre, room);
            let room4 = (room * room) * (room * room);
            let secant =
                share * (h + radiative * (odd_fourth(centre).abs() + room4) / (room - centre));
            assert!(
                (film / secant - 1.0).abs() < 1e-15,
                "a centre at {centre} K in a room at {room} K: {film} W/K against the secant's \
                 {secant} W/K"
            );
        }
    }
}
