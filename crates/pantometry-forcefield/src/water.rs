//! Rigid TIP3P water: the model's parameters, SETTLE for its constraints, and a box of it.
//!
//! This is step W2 of the explicit-water track: W1 built the periodic box and Ewald
//! ([`crate::periodic`], [`crate::ewald`]), W3 will compute benzene's hydration free energy in this
//! water and W4 the solvated complex.
//!
//! # The model
//!
//! **TIP3P as Jorgensen, Chandrasekhar, Madura, Impey and Klein published it**, *J. Chem. Phys.*
//! **79**, 926 (1983), Table I — **read secondarily**: the paper itself could not be opened (its
//! publisher and the one copy found both refused). The row below is as two readable sources
//! reproduce it, citing the paper, and they agree to every printed digit: Table 1 of Izadi,
//! Anandakrishnan and Onufriev, *J. Phys. Chem. Lett.* **5**, 3863 (2014) (read in arXiv:1408.1679),
//! and Wikipedia's "Water model" article. LAMMPS's "TIP3P water model" page gives the same
//! geometry and charges and the ε and σ the row implies:
//!
//! | r(OH) | ∠HOH | A | C | q(O) | q(H) |
//! | --- | --- | --- | --- | --- | --- |
//! | 0.9572 Å | 104.52° | 582.0 × 10³ kcal Å¹² mol⁻¹ | 595.0 kcal Å⁶ mol⁻¹ | −0.834 | +0.417 |
//!
//! The Lennard-Jones term is between oxygens only, `A/r¹² − C/r⁶`, so `ε = C²/(4A)` =
//! 0.152 07 kcal/mol and `σ = (A/C)^⅙` = 3.150 66 Å ([`epsilon`], [`sigma`]). **The hydrogens have
//! none.** CHARMM's modified TIP3P gives them a small one (LAMMPS lists it); that is another
//! model and is not this one. The geometry is rigid: three distances held, no bond or angle
//! term. Every pair inside one water is excluded, and the charges sum to zero exactly (doubling
//! 0.417 is exact in binary, and its double is the double nearest 0.834). The dipole is the closed
//! form `2 q_H r_OH cos(θ/2)` = 0.488 56 e Å = 2.347 D ([`dipole_moment`]), TIP3P's 2.35 D.
//!
//! The masses are the crate's ([`Element::mass`]): CIAAW's 15.999 and 1.008, a molecule of 18.015
//! g/mol.
//!
//! # Water with a UFF solute: a judgement
//!
//! TIP3P was fitted with OPLS-style solutes; here a solute is UFF ([`crate::energy`]), so the force
//! field is mixed and **the pair between a water oxygen and a UFF atom is an approximation** that
//! neither parameter set was fitted for. It is **UFF's own rule**, geometric in both the distance
//! and the depth (eqs 21b and 22): the oxygen enters as a UFF atom with `x = 2^⅙ σ` and `D = ε`,
//! and a hydrogen with `D = 0`, which makes every pair with it zero whatever it is combined with.
//! For: it is the rule the solute's own pairs use, so a water oxygen is one more atom type and
//! nothing else changes; and it keeps the long-range correction separable, `O(N)`
//! ([`crate::periodic`]). Against: AMBER and GAFF, and Mobley's hydration benchmark that W3 will
//! compare with, combine σ arithmetically (Lorentz–Berthelot). For benzene's aromatic carbon
//! against a water oxygen the two rules differ by 0.09% in the distance (3.690 against 3.694 Å) and
//! not at all in the depth, so the choice is small beside the mixing itself. Mobley's free
//! energies came from GAFF with TIP3P, a consistent pair; this crate's are UFF with TIP3P, which is
//! not.
//!
//! # SETTLE
//!
//! Miyamoto and Kollman, *J. Comput. Chem.* **13**, 952 (1992): the positions of a rigid
//! three-site molecule after an unconstrained step, solved analytically — **derived here from the
//! conditions the paper states, not transcribed**, and checked against SHAKE iterated to
//! convergence on the same input (`tests/a_rigid_water_against_closed_forms.rs`). The conditions:
//! the constraint impulses `m_i δ_i` are a combination of the three *old* bond vectors, as in
//! SHAKE, so they lie in the old molecular plane, sum to zero and exert no torque about any point
//! when applied at the old positions. In the frame of the old plane (`Z` its normal, `X` normal to
//! `Z` and to the new oxygen's offset from the new centre of mass, `Y = Z × X`):
//!
//! - every `Z` component is unchanged by the impulse, which fixes two of the three angles of the
//!   rigid triangle: `sin φ = z_O / r_a` and `sin ψ = (z_H2 − z_H1) / (2 r_c cos φ)`, with `r_a`
//!   the oxygen's distance from the centre of mass, `r_b` the H–H midpoint's and `r_c` half the H–H
//!   distance;
//! - the centre of mass does not move, which places the triangle;
//! - the torque about `Z` is zero, `Σ m_i (x⁰_i δy_i − y⁰_i δx_i) = 0`, which for the triangle
//!   turned by `θ` about `Z` is `α cos θ + β sin θ = γ` with `α = Σ m (x⁰ y′ − y⁰ x′)`,
//!   `β = Σ m (x⁰ x′ + y⁰ y′)` and `γ = Σ m (x⁰ y − y⁰ x)` (primes the triangle before the turn,
//!   unprimed the unconstrained positions). Of its two roots, the one with the larger `cos θ`.
//!
//! The torque about `X` and `Y` is then zero by itself, since every `δz` is. Those are six
//! conditions on six unknowns, which is why the answer is exact rather than iterated. A step
//! that turns a water past where `|sin φ| ≤ 1` has no solution and panics, naming the water.
//!
//! **Velocities** are made RATTLE-consistent separately ([`Settle::constrain_velocities`]): for
//! each water the three multipliers that remove the velocity along its three constraints, the
//! 3×3 system solved by Cramer's rule — the mass-weighted orthogonal projection onto the motions
//! that keep the triangle rigid, exactly.
//!
//! # A box of water
//!
//! [`WaterBox::lattice`] puts `k³` waters on a simple cubic lattice filling a cube at 0.997 g/cm³
//! ([`LIQUID_DENSITY`]), each centre of mass on a site and each turned by its own uniform random
//! rotation: a unit quaternion from four standard normals of `Rng::for_index(seed, water)`, which
//! is uniform on the rotations because the four-dimensional Gaussian is isotropic.
//! [`WaterBox::equilibrate`] then melts it by a stated protocol. Molecular dynamics never wraps, so
//! a water stays whole in the raw coordinates, which SETTLE requires.
//!
//! # The liquid, measured
//!
//! `tests/liquid_water_against_the_literature.rs`, in release: 216 and 512 waters at 298 K and
//! 0.997 g/cm³, `r_c` = 9 Å, Ewald at δ = 10⁻⁵, 2 fs. **Reported beside the literature, not
//! asserted**; the error bars are each series' own.
//!
//! - **U/N = −9.610 ± 0.014 and −9.612 ± 0.015 kcal/mol**, against Jorgensen's −9.86 (a cutoff
//!   Monte Carlo at the model's own density; quoted, not read) and the −9.67 that Izadi et al.'s
//!   ΔH_vap of 10.26 implies.
//! - **g_OO's first peak at 2.770 ± 0.004 and 2.782 ± 0.006 Å**, height 2.72 and 2.73 ± 0.02;
//!   Izadi et al. tabulate 2.77 Å.
//! - **P = +180 ± 37 and +139 ± 30 bar**: positive because 0.997 g/cm³ is 1.7% above TIP3P's own
//!   density, 0.980.
//! - **D**, 10⁻⁵ cm²/s, from NVE: 5.80 ± 0.27 at 304 K (216) and 5.04 ± 0.13 at 293 K (512) — each
//!   NVE run sits at the temperature its starting energy gives — and with Yeh and Hummer's
//!   correction `k_BTξ/(6πηL)` at each run's own temperature, ξ = 2.837297, η = 0.308 mPa s at
//!   298 K throughout (*J. Phys. Chem. B* **108**, 15873 (2004), read; η's own temperature
//!   dependence not included), 6.90 and 5.84, bracketing their 6.05–6.11 at 298 K. Mahoney and Jorgensen
//!   (*J. Chem. Phys.* **114**, 363 (2001), read) report 5.19 ± 0.08 uncorrected; experiment is
//!   2.30. TIP3P diffuses too fast, as it is known to.
//!
//! # The cost
//!
//! Release, one core, at 2 fs with Ewald at δ = 10⁻⁵ and `r_c` = 9 Å: 7.0 ms a step for 216
//! waters, 26.2 for 512, 55.2 for 1 000, of which SETTLE and one projection are 0.03, 0.07 and 0.14.
//! The force field is the cost.

use crate::ccd::{Element, ANGSTROM};
use crate::dynamics::{Bath, MolecularDynamics};
use crate::ewald::EwaldParameters;
use crate::periodic::{PeriodicBox, PeriodicForceField};
use crate::uff::KCAL_PER_MOL;
use pantometry_core::Rng;

/// The oxygen's charge, elementary charges (Table I).
pub const OXYGEN_CHARGE: f64 = -0.834;

/// Each hydrogen's charge, elementary charges (Table I).
pub const HYDROGEN_CHARGE: f64 = 0.417;

/// The O–H distance, metres (Table I: 0.9572 Å).
pub const OH_LENGTH: f64 = 0.9572e-10;

/// The H–O–H angle, degrees (Table I).
pub const HOH_ANGLE_DEGREES: f64 = 104.52;

/// Table I's `A`, kcal Å¹² mol⁻¹: the O–O repulsion `A/r¹²`.
pub const LJ_A: f64 = 582.0e3;

/// Table I's `C`, kcal Å⁶ mol⁻¹: the O–O dispersion `−C/r⁶`.
pub const LJ_C: f64 = 595.0;

/// Liquid water's density at 298 K and 1 atm, kg m⁻³: 0.997 g/cm³.
pub const LIQUID_DENSITY: f64 = 997.0;

/// The H–O–H angle, radians.
pub fn hoh_angle() -> f64 {
    HOH_ANGLE_DEGREES.to_radians()
}

/// The H–H distance, metres: `2 r_OH sin(θ/2)`, 1.5139 Å.
pub fn hh_length() -> f64 {
    2.0 * OH_LENGTH * (0.5 * hoh_angle()).sin()
}

/// The O–O Lennard-Jones well depth `ε = C²/(4A)`, joules per molecule: 0.152 07 kcal/mol.
pub fn epsilon() -> f64 {
    LJ_C * LJ_C / (4.0 * LJ_A) * KCAL_PER_MOL
}

/// The O–O Lennard-Jones diameter `σ = (A/C)^⅙`, metres: 3.150 66 Å.
pub fn sigma() -> f64 {
    cube_root(LJ_A / LJ_C).sqrt() * ANGSTROM
}

/// The oxygen as a UFF atom, `[x, D]` in metres and joules: `x = 2^⅙ σ`, the distance of the
/// minimum, and `D = ε`. See the module documentation for why water mixes with a solute by UFF's
/// rule.
pub fn oxygen_vdw() -> [f64; 2] {
    // x⁶ = 2σ⁶ = 2A/C.
    [cube_root(2.0 * LJ_A / LJ_C).sqrt() * ANGSTROM, epsilon()]
}

/// The dipole moment, the closed form `2 q_H r_OH cos(θ/2)`, elementary charges × metres.
pub fn dipole_moment() -> f64 {
    2.0 * HYDROGEN_CHARGE * OH_LENGTH * (0.5 * hoh_angle()).cos()
}

/// One water's charges, oxygen first.
pub const CHARGES: [f64; 3] = [OXYGEN_CHARGE, HYDROGEN_CHARGE, HYDROGEN_CHARGE];

/// One water's masses, kilograms, oxygen first.
pub fn masses() -> [f64; 3] {
    [Element::O.mass(), Element::H.mass(), Element::H.mass()]
}

/// The molecule's mass, kilograms.
pub fn molecular_mass() -> f64 {
    Element::O.mass() + 2.0 * Element::H.mass()
}

/// The rigid triangle about its centre of mass, metres, oxygen first: the oxygen on `+y`, the
/// hydrogens at `∓x`, in the `z = 0` plane — SETTLE's canonical frame.
pub fn body_frame() -> [[f64; 3]; 3] {
    let s = Settle::tip3p(Vec::new());
    [[0.0, s.ra, 0.0], [-s.rc, -s.rb, 0.0], [s.rc, -s.rb, 0.0]]
}

/// The cube root of `x > 0` by bisection on its bits' range, `+ − × ÷` only: the same on every
/// platform, where the platform's `cbrt` need not be. Within an ulp of the exact root.
fn cube_root(x: f64) -> f64 {
    let (mut lo, mut hi) = if x < 1.0 { (x, 1.0) } else { (1.0, x) };
    for _ in 0..2100 {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            break;
        }
        if mid * mid * mid < x {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    hi
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn unit(a: [f64; 3]) -> [f64; 3] {
    let l = dot(a, a).sqrt();
    [a[0] / l, a[1] / l, a[2] / l]
}

/// SETTLE for a set of identical rigid three-site molecules: see the module documentation.
#[derive(Clone, Debug, PartialEq)]
pub struct Settle {
    waters: Vec<[usize; 3]>,
    m_o: f64,
    m_h: f64,
    /// The oxygen's distance from the centre of mass.
    ra: f64,
    /// The H–H midpoint's distance from the centre of mass.
    rb: f64,
    /// Half the H–H distance.
    rc: f64,
}

impl Settle {
    /// SETTLE for molecules `waters` — atom indices, the apex first — of apex mass `m_o`, the two
    /// others `m_h` each (kilograms), apex distance `oh` and apex angle `angle` (metres, radians).
    ///
    /// # Panics
    ///
    /// If a mass or the distance is not positive and finite, the angle is not in `(0, π)`, or an
    /// index appears twice.
    pub fn new(waters: Vec<[usize; 3]>, oh: f64, angle: f64, m_o: f64, m_h: f64) -> Settle {
        assert!(
            [oh, m_o, m_h].iter().all(|x| x.is_finite() && *x > 0.0),
            "SETTLE's distance and masses must be positive and finite"
        );
        assert!(
            angle > 0.0 && angle < std::f64::consts::PI,
            "SETTLE's angle must be in (0, π)"
        );
        let mut seen: Vec<usize> = waters.iter().flatten().copied().collect();
        seen.sort_unstable();
        assert!(
            seen.windows(2).all(|w| w[0] != w[1]),
            "an atom is in two rigid molecules, or twice in one"
        );
        let height = oh * (0.5 * angle).cos();
        let ra = 2.0 * m_h * height / (m_o + 2.0 * m_h);
        Settle {
            waters,
            m_o,
            m_h,
            ra,
            rb: height - ra,
            rc: oh * (0.5 * angle).sin(),
        }
    }

    /// SETTLE for TIP3P waters: [`OH_LENGTH`], [`hoh_angle`] and the crate's oxygen and hydrogen
    /// masses.
    pub fn tip3p(waters: Vec<[usize; 3]>) -> Settle {
        Settle::new(
            waters,
            OH_LENGTH,
            hoh_angle(),
            Element::O.mass(),
            Element::H.mass(),
        )
    }

    /// The molecules, atom indices, apex first.
    pub fn waters(&self) -> &[[usize; 3]] {
        &self.waters
    }

    /// The apex mass and each other mass, kilograms.
    pub fn masses(&self) -> [f64; 2] {
        [self.m_o, self.m_h]
    }

    /// The three constrained distances of one molecule, metres: apex to each other atom, then
    /// between the two.
    pub fn lengths(&self) -> [f64; 3] {
        let h = self.ra + self.rb;
        let oh = (h * h + self.rc * self.rc).sqrt();
        [oh, oh, 2.0 * self.rc]
    }

    /// Moves every molecule of `after` — unconstrained positions one step on from `before`, where
    /// every molecule was rigid — to the rigid positions SETTLE gives. Atoms in no molecule are not
    /// touched.
    ///
    /// # Panics
    ///
    /// If `before` and `after` differ in length, or a molecule has turned too far in one step for
    /// a solution to exist.
    pub fn constrain_positions(&self, before: &[[f64; 3]], after: &mut [[f64; 3]]) {
        self.constrain_positions_where(before, after, |_| true);
    }

    pub(crate) fn constrain_positions_where(
        &self,
        before: &[[f64; 3]],
        after: &mut [[f64; 3]],
        active: impl Fn(usize) -> bool,
    ) {
        assert_eq!(before.len(), after.len(), "one old position per new one");
        for (w, &[o, h1, h2]) in self.waters.iter().enumerate() {
            if !active(w) {
                continue;
            }
            let moved = self.settle_one(
                [before[o], before[h1], before[h2]],
                [after[o], after[h1], after[h2]],
            );
            match moved {
                Some([a, b, c]) => {
                    after[o] = a;
                    after[h1] = b;
                    after[h2] = c;
                }
                None => panic!(
                    "SETTLE has no solution for the molecule of atoms {o}, {h1}, {h2}: it turned \
                     too far in one step"
                ),
            }
        }
    }

    /// One molecule: the rigid positions for unconstrained `new` one step from rigid `old`.
    fn settle_one(&self, old: [[f64; 3]; 3], new: [[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
        let m = [self.m_o, self.m_h, self.m_h];
        let total = self.m_o + 2.0 * self.m_h;
        let centre = |p: &[[f64; 3]; 3]| -> [f64; 3] {
            [0, 1, 2].map(|k| (m[0] * p[0][k] + m[1] * p[1][k] + m[2] * p[2][k]) / total)
        };
        let (d0, d1) = (centre(&old), centre(&new));
        let r0 = old.map(|p| sub(p, d0));
        let r1 = new.map(|p| sub(p, d1));
        // The old plane's normal, and the frame.
        let z = unit(cross(sub(old[1], old[0]), sub(old[2], old[0])));
        let x = unit(cross(r1[0], z));
        let y = cross(z, x);
        let frame = |p: [f64; 3]| [dot(p, x), dot(p, y), dot(p, z)];
        let p0 = r0.map(frame);
        let p1 = r1.map(frame);
        // Out of the old plane, nothing moves: φ and ψ.
        let sin_phi = p1[0][2] / self.ra;
        if sin_phi.is_nan() || sin_phi.abs() > 1.0 {
            return None;
        }
        let cos_phi = (1.0 - sin_phi * sin_phi).sqrt();
        let sin_psi = (p1[2][2] - p1[1][2]) / (2.0 * self.rc * cos_phi);
        if sin_psi.is_nan() || sin_psi.abs() > 1.0 {
            return None;
        }
        let cos_psi = (1.0 - sin_psi * sin_psi)
            .sqrt()
            .copysign(p1[2][0] - p1[1][0]);
        let a2 = [0.0, self.ra * cos_phi, self.ra * sin_phi];
        let mid = [0.0, -self.rb * cos_phi, -self.rb * sin_phi];
        let half = [
            self.rc * cos_psi,
            -self.rc * sin_psi * sin_phi,
            self.rc * sin_psi * cos_phi,
        ];
        let b2 = sub(mid, half);
        let c2 = [mid[0] + half[0], mid[1] + half[1], mid[2] + half[2]];
        let p2 = [a2, b2, c2];
        // In the plane, no torque about Z: α cos θ + β sin θ = γ.
        let (mut alpha, mut beta, mut gamma) = (0.0, 0.0, 0.0);
        for i in 0..3 {
            let ([x0, y0, _], [x2, y2, _], [xu, yu, _]) = (p0[i], p2[i], p1[i]);
            alpha += m[i] * (x0 * y2 - y0 * x2);
            beta += m[i] * (x0 * x2 + y0 * y2);
            gamma += m[i] * (x0 * yu - y0 * xu);
        }
        let r2 = alpha * alpha + beta * beta;
        let root = (r2 - gamma * gamma).max(0.0).sqrt();
        // The two roots; the one with the larger cos θ, the smaller turn.
        let (c_plus, s_plus) = (
            (alpha * gamma + beta * root) / r2,
            (beta * gamma - alpha * root) / r2,
        );
        let (c_minus, s_minus) = (
            (alpha * gamma - beta * root) / r2,
            (beta * gamma + alpha * root) / r2,
        );
        let (c, s) = if c_plus >= c_minus {
            (c_plus, s_plus)
        } else {
            (c_minus, s_minus)
        };
        Some(p2.map(|[px, py, pz]| {
            let (qx, qy) = (px * c - py * s, px * s + py * c);
            [0, 1, 2].map(|k| d1[k] + qx * x[k] + qy * y[k] + pz * z[k])
        }))
    }

    /// Removes from `velocities` every component along a constraint of a molecule at `at`, by the
    /// mass-weighted orthogonal projection: afterwards `(r_a − r_b)·(v_a − v_b) = 0` for each of
    /// the three pairs of every molecule, to rounding. Atoms in no molecule are not touched.
    ///
    /// # Panics
    ///
    /// If `at` and `velocities` differ in length.
    pub fn constrain_velocities(&self, at: &[[f64; 3]], velocities: &mut [[f64; 3]]) {
        self.constrain_velocities_where(at, velocities, |_| true);
    }

    pub(crate) fn constrain_velocities_where(
        &self,
        at: &[[f64; 3]],
        v: &mut [[f64; 3]],
        active: impl Fn(usize) -> bool,
    ) {
        assert_eq!(at.len(), v.len(), "one velocity per position");
        let (io, ih) = (1.0 / self.m_o, 1.0 / self.m_h);
        for (w, &[o, h1, h2]) in self.waters.iter().enumerate() {
            if !active(w) {
                continue;
            }
            // Constraint 1 is O–H1, 2 is O–H2 and 3 is H1–H2, each along r_a − r_b.
            let r = [sub(at[o], at[h1]), sub(at[o], at[h2]), sub(at[h1], at[h2])];
            let u = [sub(v[o], v[h1]), sub(v[o], v[h2]), sub(v[h1], v[h2])];
            let a = [
                [
                    (io + ih) * dot(r[0], r[0]),
                    io * dot(r[0], r[1]),
                    -ih * dot(r[0], r[2]),
                ],
                [
                    io * dot(r[0], r[1]),
                    (io + ih) * dot(r[1], r[1]),
                    ih * dot(r[1], r[2]),
                ],
                [
                    -ih * dot(r[0], r[2]),
                    ih * dot(r[1], r[2]),
                    2.0 * ih * dot(r[2], r[2]),
                ],
            ];
            let rhs = [0, 1, 2].map(|k| -dot(r[k], u[k]));
            let mu = solve3(a, rhs);
            for k in 0..3 {
                v[o][k] += io * (mu[0] * r[0][k] + mu[1] * r[1][k]);
                v[h1][k] += ih * (-mu[0] * r[0][k] + mu[2] * r[2][k]);
                v[h2][k] += ih * (-mu[1] * r[1][k] - mu[2] * r[2][k]);
            }
        }
    }
}

/// `A x = b` for a 3×3 `A` by Cramer's rule.
fn solve3(a: [[f64; 3]; 3], b: [f64; 3]) -> [f64; 3] {
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(a);
    [0, 1, 2].map(|c| {
        let mut m = a;
        for r in 0..3 {
            m[r][c] = b[r];
        }
        det(m) / d
    })
}

/// A cube of rigid TIP3P waters: its box, every atom's position (oxygen, then its two hydrogens,
/// water by water) and each water's atoms. See the module documentation.
#[derive(Clone, Debug, PartialEq)]
pub struct WaterBox {
    cell: PeriodicBox,
    positions: Vec<[f64; 3]>,
}

impl WaterBox {
    /// `per_side³` waters on a simple cubic lattice at [`LIQUID_DENSITY`], turned by rotations
    /// drawn from `seed`.
    ///
    /// # Panics
    ///
    /// If `per_side` is zero.
    pub fn lattice(per_side: usize, seed: u64) -> WaterBox {
        WaterBox::lattice_at(per_side, LIQUID_DENSITY, seed)
    }

    /// The same at `density` kg m⁻³.
    ///
    /// # Panics
    ///
    /// If `per_side` is zero or the density is not positive and finite.
    pub fn lattice_at(per_side: usize, density: f64, seed: u64) -> WaterBox {
        assert!(per_side > 0, "a box of no water");
        assert!(
            density.is_finite() && density > 0.0,
            "a density must be positive and finite"
        );
        let n = per_side * per_side * per_side;
        let side = cube_root(n as f64 * molecular_mass() / density);
        let spacing = side / per_side as f64;
        let body = body_frame();
        let mut positions = Vec::with_capacity(3 * n);
        for w in 0..n {
            let g = [
                w / (per_side * per_side),
                (w / per_side) % per_side,
                w % per_side,
            ];
            let site = g.map(|i| (i as f64 + 0.5) * spacing);
            let rotation = uniform_rotation(seed, w as u64);
            for p in body {
                let q = [0, 1, 2].map(|r| dot(rotation[r], p));
                positions.push([site[0] + q[0], site[1] + q[1], site[2] + q[2]]);
            }
        }
        WaterBox {
            cell: PeriodicBox::cubic(side),
            positions,
        }
    }

    /// The same box without every water that has an atom within `clearance` metres of an atom of
    /// `solute` (metres, any image; the distance by minimum image): the waters a solute placed in
    /// the box would overlap. The others keep their order and their bits. The solute's positions
    /// go first in a solvated system ([`PeriodicForceField::solvated`]), and these after them.
    ///
    /// # Panics
    ///
    /// If `clearance` is negative or not finite.
    pub fn without_overlaps(&self, solute: &[[f64; 3]], clearance: f64) -> WaterBox {
        assert!(
            clearance.is_finite() && clearance >= 0.0,
            "a clearance must be finite and not negative"
        );
        let c2 = clearance * clearance;
        let near = |p: [f64; 3]| {
            solute.iter().any(|s| {
                let d = self.cell.minimum_image(sub(p, *s));
                dot(d, d) < c2
            })
        };
        // By index: `chunks_exact(3)` draws clippy's `as_chunks`, stabilised in 1.88, and this
        // crate builds on 1.78.
        let positions = (0..self.count())
            .map(|w| &self.positions[3 * w..3 * w + 3])
            .filter(|w| !w.iter().any(|p| near(*p)))
            .flatten()
            .copied()
            .collect();
        WaterBox {
            cell: self.cell,
            positions,
        }
    }

    /// The box.
    pub fn cell(&self) -> PeriodicBox {
        self.cell
    }

    /// Every atom's position, metres.
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }

    /// How many waters.
    pub fn count(&self) -> usize {
        self.positions.len() / 3
    }

    /// Each water's atoms, oxygen first.
    pub fn waters(&self) -> Vec<[usize; 3]> {
        (0..self.count())
            .map(|w| [3 * w, 3 * w + 1, 3 * w + 2])
            .collect()
    }

    /// Every atom's mass, kilograms.
    pub fn masses(&self) -> Vec<f64> {
        (0..self.count()).flat_map(|_| masses()).collect()
    }

    /// The TIP3P force field of this box: Ewald with real-space cutoff `cutoff` (metres) at
    /// `accuracy` ([`EwaldParameters::for_accuracy`]).
    ///
    /// # Panics
    ///
    /// As [`EwaldParameters::for_accuracy`].
    pub fn force_field(&self, cutoff: f64, accuracy: f64) -> PeriodicForceField {
        let p = EwaldParameters::for_accuracy(&self.cell, cutoff, accuracy);
        PeriodicForceField::tip3p(self.cell, p, self.count())
    }

    /// Molecular dynamics of this box's atoms with every water held rigid by SETTLE, at rest.
    pub fn dynamics(&self) -> MolecularDynamics {
        MolecularDynamics::new(self.masses()).with_constraints(Settle::tip3p(self.waters()))
    }

    /// Melts the lattice at `temperature` kelvin and returns the dynamics to go on with, its
    /// velocities thermal: **0.5 ps at a 0.5 fs step in a 50 ps⁻¹ bath**, which takes out the
    /// lattice's strain without a step a hydrogen could cross a neighbour in, **then
    /// `picoseconds` at 2 fs in a 5 ps⁻¹ bath**. Velocities and kicks from `seed`. The positions
    /// are this box's, moved.
    pub fn equilibrate(
        &mut self,
        field: &PeriodicForceField,
        temperature: f64,
        picoseconds: f64,
        seed: u64,
    ) -> MolecularDynamics {
        let bath = |friction: f64, s: u64| Bath::Langevin {
            temperature,
            friction,
            seed: s,
        };
        let mut md = self.dynamics().with_bath(bath(50e12, seed)).thermalised(
            &self.positions,
            temperature,
            seed,
        );
        md.run(field, &mut self.positions, 0.5e-15, 1000);
        let mut md = md.with_bath(bath(5e12, seed ^ 0x5EED));
        let steps = (picoseconds * 1e-12 / 2e-15).round() as usize;
        md.run(field, &mut self.positions, 2e-15, steps);
        md
    }
}

/// A uniform random rotation, as the rows of its matrix: the unit quaternion of four standard
/// normals of `Rng::for_index(seed, index)`.
fn uniform_rotation(seed: u64, index: u64) -> [[f64; 3]; 3] {
    let mut rng = Rng::for_index(seed, index);
    let q = [
        rng.gaussian(),
        rng.gaussian(),
        rng.gaussian(),
        rng.gaussian(),
    ];
    let l = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    let [w, x, y, z] = q.map(|c| c / l);
    [
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - w * z),
            2.0 * (x * z + w * y),
        ],
        [
            2.0 * (x * y + w * z),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - w * x),
        ],
        [
            2.0 * (x * z - w * y),
            2.0 * (y * z + w * x),
            1.0 - 2.0 * (x * x + y * y),
        ],
    ]
}

/// The molecules' number density at `density` kg m⁻³, per cubic metre.
pub fn number_density(density: f64) -> f64 {
    density / molecular_mass()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::uff::AVOGADRO;

    #[test]
    fn avogadro_is_the_crates() {
        // The molecular mass is the masses' sum, 18.015 g/mol.
        assert!((molecular_mass() * AVOGADRO * 1e3 - 18.015).abs() < 1e-12);
    }
}
