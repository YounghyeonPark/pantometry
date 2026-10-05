//! A Boresch restraint: six harmonic terms that hold a ligand's position and orientation relative
//! to its receptor, and the analytic free energy of releasing them to the standard state.
//!
//! > S. Boresch, F. Tettinger, M. Leitgeb and M. Karplus, "Absolute binding free energies: a
//! > quantitative approach for their calculation", *J. Phys. Chem. B* **107**, 9535–9551 (2003),
//! > [doi:10.1021/jp0217839](https://doi.org/10.1021/jp0217839).
//!
//! **That paper was not opened here** — it is behind the publisher's paywall and no readable copy
//! was found. Everything below is **secondary**, taken from a paper whose full text was read:
//! J. Clark, M. Viner, …, J. Michel, "Comparison of receptor–ligand restraint schemes for
//! alchemical absolute binding free energy calculations", *J. Chem. Theory Comput.* **19**, 3686
//! (2023), [doi:10.1021/acs.jctc.3c00139](https://doi.org/10.1021/acs.jctc.3c00139), read in PMC
//! (PMC10308817) with its equations read off the published images: its Figure 3 for the six
//! coordinates, its eq 6 for the integral, and its eq 7, which it attributes to Boresch et al., for
//! the closed form. The task that asked for this module names that closed form "Boresch eq 32";
//! the equation number is the task's, not read.
//!
//! # The six coordinates (Clark et al., Figure 3)
//!
//! Three receptor atoms `a`, `b`, `c` and three ligand atoms `A`, `B`, `C`, chained
//! `c – b – a ⋯ A – B – C`:
//!
//! - the distance `r = |A − a|`;
//! - the angles `θ_A = ∠(b, a, A)` and `θ_B = ∠(a, A, B)`;
//! - the dihedrals `φ_A = (c, b, a, A)`, `φ_B = (b, a, A, B)` and `φ_C = (a, A, B, C)`, in the
//!   convention of [`dihedral`].
//!
//! Each is held by `½ K (ξ − ξ₀)²`, a dihedral's difference wrapped to (−π, π]. **The ½ is a
//! choice that the closed form fixes**: Clark et al.'s eq 8 writes their distance restraints as
//! `½ K (r − r₀)²`, and the `(2π k_BT)³` in eq 7 is six Gaussian integrals
//! `∫ e^(−K x²/2k_BT) dx = √(2π k_BT / K)` — the integral of `K x²` without the ½ would give
//! `(π k_BT)³`. (Mobley, Chodera and Dill, *J. Chem. Phys.* **125**, 084902 (2006), eq 13, write
//! `K₀ (ξ − ξ₀)²` with K₀ = 10 kcal mol⁻¹ Å⁻² or rad⁻², as read through a summary of its PMC
//! text; in this module's convention that is K = 20.)
//!
//! # The free energy of releasing it (Clark et al., eq 7; Boresch et al.)
//!
//! With the ligand decoupled, the restraint is the only thing holding it, and
//!
//! `ΔG°_release = −k_BT ln [ 8π² V° √(K_r K_θA K_θB K_φA K_φB K_φC) / (r₀² sin θ_A0 sin θ_B0 (2π k_BT)³) ]`
//!
//! is the free energy of letting it go to a standard volume `V°` per molecule and every
//! orientation (`8π²`). It is negative for any useful restraint: releasing gains entropy. The
//! free energy of **turning the restraint on** in the decoupled state is its negative.
//!
//! **Where it comes from, worked here so that the test can check each step.** The restrained
//! ligand's position is `A`'s, in spherical coordinates about `a` with the `a`–`b` axis as pole:
//! `r² sin θ_A dr dθ_A dφ_A`. Its orientation about `A` is `B`'s direction, polar angle `θ_B` from
//! `A → a` and azimuth `φ_B`, with measure `sin θ_B dθ_B dφ_B`, and `C`'s turn about `A–B`,
//! `dφ_C`; the three together measure `8π²`. So the restrained configurational integral is
//! Clark et al.'s eq 6,
//!
//! `Z_r = ∫ r² e^(−βu_r) dr ∫ sin θ_A e^(−βu_θA) dθ_A ∫ e^(−βu_φA) dφ_A · ∫ sin θ_B e^(−βu_θB) dθ_B ∫ e^(−βu_φB) dφ_B ∫ e^(−βu_φC) dφ_C`,
//!
//! and `ΔG°_release = −k_BT ln (8π² V° / Z_r)`. **Eq 7 is eq 6 with `r²` and `sin θ` taken at
//! their reference values and every Gaussian extended over the whole line.** Both approximations
//! have closed forms: over the line, `∫ r² e^(−K(r − r₀)²/2k_BT) dr = √(2π k_BT/K) (r₀² + k_BT/K)`
//! and `∫ sin θ e^(−K(θ − θ₀)²/2k_BT) dθ = √(2π k_BT/K) sin θ₀ e^(−k_BT/2K)`, exactly. So eq 7
//! differs from eq 6 by
//!
//! `+k_BT [ln(1 + k_BT/(K_r r₀²)) − k_BT/(2K_θA) − k_BT/(2K_θB)]`
//!
//! (eq 6 minus eq 7) plus the Gaussians' tails beyond `r < 0`, `θ ∉ [0, π]` and `|φ − φ₀| > π`.
//! [`Boresch::release_free_energy_extended`] is eq 7 with that correction, and
//! `tests/the_free_energy_against_closed_forms.rs` holds it to a quadrature of eq 6 and the
//! position and the orientation separately to quadratures **in Cartesian coordinates**, which
//! have no Jacobian to get wrong.
//!
//! # The standard state
//!
//! `V°` is the volume per molecule at 1 mol/L: `10⁻³ m³ / N_A` = 1660.539 Å³
//! ([`STANDARD_VOLUME`]).
//!
//! # Choosing the atoms
//!
//! [`Boresch::choose`] takes candidate triples for each side and keeps the pair whose four
//! angles — `θ_A`, `θ_B`, and the two that the dihedrals `φ_A` and `φ_C` hinge on, `∠(c, b, a)`
//! and `∠(A, B, C)` — are furthest from straight: the largest smallest sine, among pairs with `r₀`
//! in a stated range. A dihedral is undefined where one of its angles is straight, and the
//! restraint's force on it grows without bound there (Clark et al. §2.3 on "the collinearity of
//! any three contiguous anchor points"; they discard sets with `θ_A` or `θ_B` below 30° or above
//! 150°). The reference values are the coordinates at the positions given. This is a rule about
//! geometry only: it does not look for the stiffest receptor–ligand contacts, which Clark et al.'s
//! variance-based selection does, so a restraint chosen here may hold a ligand where it does not
//! sit for long.

use crate::minimise::{dihedral, DihedralRestraint};
use crate::uff::AVOGADRO;
use pantometry_units::BOLTZMANN;
use std::f64::consts::PI;

/// The volume per molecule at the standard concentration, 1 mol/L, in m³: `10⁻³ m³ / N_A`,
/// 1.660 539 × 10⁻²⁷ m³ = 1660.539 Å³.
pub const STANDARD_VOLUME: f64 = 1e-3 / AVOGADRO;

/// A Boresch restraint between receptor atoms `a`, `b`, `c` and ligand atoms `A`, `B`, `C`: see
/// the module documentation for the six coordinates and their order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Boresch {
    /// `[a, b, c]`, atom indices.
    pub receptor: [usize; 3],
    /// `[A, B, C]`, atom indices.
    pub ligand: [usize; 3],
    /// The reference values `[r₀, θ_A0, θ_B0, φ_A0, φ_B0, φ_C0]`: metres, then radians.
    pub reference: [f64; 6],
    /// `[K_r, K_θA, K_θB, K_φA, K_φB, K_φC]`, in the restraint `½ K (ξ − ξ₀)²`: J m⁻², then
    /// J rad⁻².
    pub force_constants: [f64; 6],
}

/// The angle at `v` between `p` and `q`, radians in [0, π], as `atan2(|u × w|, u · w)` with
/// `u = p − v` and `w = q − v`, and its gradient with respect to `p`, `v` and `q`.
fn angle_and_gradient(p: [f64; 3], v: [f64; 3], q: [f64; 3]) -> (f64, [[f64; 3]; 3]) {
    let u = sub(p, v);
    let w = sub(q, v);
    let (lu, lw) = (norm(u), norm(w));
    let theta = norm(cross(u, w)).atan2(dot(u, w));
    let (s, c) = theta.sin_cos();
    if s == 0.0 || lu == 0.0 || lw == 0.0 {
        return (theta, [[0.0; 3]; 3]);
    }
    // ∂cosθ/∂u = (ŵ − cosθ û)/|u|, and ∂θ = −∂cosθ / sinθ.
    let mut gp = [0.0; 3];
    let mut gq = [0.0; 3];
    for k in 0..3 {
        gp[k] = -(w[k] / lw - c * u[k] / lu) / (lu * s);
        gq[k] = -(u[k] / lu - c * w[k] / lw) / (lw * s);
    }
    let gv = [-(gp[0] + gq[0]), -(gp[1] + gq[1]), -(gp[2] + gq[2])];
    (theta, [gp, gv, gq])
}

/// The angle at `v` between `p` and `q`, radians.
fn angle(p: [f64; 3], v: [f64; 3], q: [f64; 3]) -> f64 {
    angle_and_gradient(p, v, q).0
}

impl Boresch {
    /// The restraint on `receptor` `[a, b, c]` and `ligand` `[A, B, C]` with its reference values
    /// taken at positions `at` (metres), and force constants `[K_r, K_θA, K_θB, K_φA, K_φB,
    /// K_φC]` (J m⁻², then J rad⁻²).
    ///
    /// # Panics
    ///
    /// If an index is past `at`, or a force constant is not positive and finite.
    pub fn at(
        at: &[[f64; 3]],
        receptor: [usize; 3],
        ligand: [usize; 3],
        force_constants: [f64; 6],
    ) -> Boresch {
        assert!(
            force_constants.iter().all(|k| k.is_finite() && *k > 0.0),
            "every force constant must be positive and finite"
        );
        let mut b = Boresch {
            receptor,
            ligand,
            reference: [0.0; 6],
            force_constants,
        };
        b.reference = b.coordinates(at);
        b
    }

    /// The six coordinates `[r, θ_A, θ_B, φ_A, φ_B, φ_C]` at `at`: metres, then radians.
    ///
    /// # Panics
    ///
    /// If an index is past `at`.
    pub fn coordinates(&self, at: &[[f64; 3]]) -> [f64; 6] {
        let [a, b, c] = self.receptor;
        let [la, lb, lc] = self.ligand;
        [
            norm(sub(at[la], at[a])),
            angle(at[b], at[a], at[la]),
            angle(at[a], at[la], at[lb]),
            dihedral(at, [c, b, a, la]),
            dihedral(at, [b, a, la, lb]),
            dihedral(at, [a, la, lb, lc]),
        ]
    }

    /// The four angles a straight line would make the dihedrals undefined at, radians:
    /// `[θ_A, θ_B, ∠(c, b, a), ∠(A, B, C)]` at `at`.
    pub fn hinge_angles(&self, at: &[[f64; 3]]) -> [f64; 4] {
        let [a, b, c] = self.receptor;
        let [la, lb, lc] = self.ligand;
        [
            angle(at[b], at[a], at[la]),
            angle(at[a], at[la], at[lb]),
            angle(at[c], at[b], at[a]),
            angle(at[la], at[lb], at[lc]),
        ]
    }

    fn dihedrals(&self) -> [DihedralRestraint; 3] {
        let [a, b, c] = self.receptor;
        let [la, lb, lc] = self.ligand;
        let atoms = [[c, b, a, la], [b, a, la, lb], [a, la, lb, lc]];
        [0, 1, 2].map(|k| DihedralRestraint {
            atoms: atoms[k],
            target: self.reference[3 + k],
            stiffness: self.force_constants[3 + k],
        })
    }

    /// The restraint's energy at `at`, joules per molecule: `Σ ½ K (ξ − ξ₀)²`, each dihedral's
    /// difference wrapped to (−π, π].
    pub fn energy(&self, at: &[[f64; 3]]) -> f64 {
        let x = self.coordinates(at);
        let mut e = 0.0;
        for ((x, x0), k) in x
            .iter()
            .zip(&self.reference)
            .zip(&self.force_constants)
            .take(3)
        {
            let d = x - x0;
            e += 0.5 * k * d * d;
        }
        for d in self.dihedrals() {
            e += d.energy_at(at);
        }
        e
    }

    /// Adds `scale` times the restraint's force to `forces` (newtons) and returns the restraint's
    /// **unscaled** energy — what `∂U/∂λ` is for a restraint scaled linearly by `λ`.
    ///
    /// The distance's and the angles' gradients are written here, the dihedrals' are
    /// [`DihedralRestraint::add_forces`]'s. Checked against central differences in the tests.
    ///
    /// # Panics
    ///
    /// If an index is past `at` or `forces`.
    pub fn add_forces(&self, at: &[[f64; 3]], forces: &mut [[f64; 3]], scale: f64) -> f64 {
        let [a, b, _] = self.receptor;
        let [la, lb, _] = self.ligand;
        let mut e = 0.0;
        // The distance.
        let d = sub(at[la], at[a]);
        let r = norm(d);
        let dr = r - self.reference[0];
        e += 0.5 * self.force_constants[0] * dr * dr;
        if r > 0.0 {
            let g = scale * self.force_constants[0] * dr / r;
            for k in 0..3 {
                forces[la][k] -= g * d[k];
                forces[a][k] += g * d[k];
            }
        }
        // The two angles, θ_A = ∠(b, a, A) and θ_B = ∠(a, A, B).
        for (slot, [p, v, q]) in [(1, [b, a, la]), (2, [a, la, lb])] {
            let (theta, grad) = angle_and_gradient(at[p], at[v], at[q]);
            let dt = theta - self.reference[slot];
            e += 0.5 * self.force_constants[slot] * dt * dt;
            let de = scale * self.force_constants[slot] * dt;
            for (atom, g) in [(p, grad[0]), (v, grad[1]), (q, grad[2])] {
                for k in 0..3 {
                    forces[atom][k] -= de * g[k];
                }
            }
        }
        // The three dihedrals: the unscaled energy, and the scaled force.
        for mut t in self.dihedrals() {
            e += t.energy_at(at);
            t.stiffness *= scale;
            t.add_forces(at, forces);
        }
        e
    }

    /// `ΔG°_release` at `temperature` kelvin, joules per molecule: the free energy of releasing
    /// the restraint on a decoupled ligand to the standard state, Clark et al.'s eq 7 (from
    /// Boresch et al.) — see the module documentation. Negative. The free energy of turning the
    /// restraint on in the decoupled state is its negative.
    ///
    /// # Panics
    ///
    /// If the temperature is not positive and finite.
    pub fn release_free_energy(&self, temperature: f64) -> f64 {
        assert!(
            temperature.is_finite() && temperature > 0.0,
            "a temperature must be positive"
        );
        let kt = BOLTZMANN.to_si() * temperature;
        let [r0, ta, tb, ..] = self.reference;
        let half_log_k: f64 = self.force_constants.iter().map(|k| 0.5 * k.ln()).sum();
        let log_argument = STANDARD_VOLUME.ln() + (8.0 * PI * PI).ln() + half_log_k
            - 2.0 * r0.ln()
            - ta.sin().ln()
            - tb.sin().ln()
            - 3.0 * (2.0 * PI * kt).ln();
        -kt * log_argument
    }

    /// [`Boresch::release_free_energy`] with the two approximations that separate eq 7 from eq 6
    /// undone in closed form: `r²` and `sin θ` averaged over their Gaussians instead of taken at
    /// the reference, `+k_BT [ln(1 + k_BT/(K_r r₀²)) − k_BT/(2K_θA) − k_BT/(2K_θB)]` more. Exact
    /// for the Gaussians extended over the whole line; what is left against eq 6 is their tails
    /// beyond the coordinates' ranges. See the module documentation.
    ///
    /// # Panics
    ///
    /// As [`Boresch::release_free_energy`].
    pub fn release_free_energy_extended(&self, temperature: f64) -> f64 {
        let kt = BOLTZMANN.to_si() * temperature;
        let [kr, ka, kb, ..] = self.force_constants;
        let r0 = self.reference[0];
        self.release_free_energy(temperature)
            + kt * ((kt / (kr * r0 * r0)).ln_1p() - kt / (2.0 * ka) - kt / (2.0 * kb))
    }

    /// The restraint, among every pair of a `receptor` candidate `[a, b, c]` and a `ligand`
    /// candidate `[A, B, C]` whose distance `r₀ = |A − a|` at `at` is within `distance` (metres,
    /// inclusive), whose four hinge angles ([`Boresch::hinge_angles`]) are furthest from straight:
    /// the largest smallest sine. Ties go to the first pair, receptor candidates outermost. Its
    /// reference values are the coordinates at `at` and its force constants `force_constants`, as
    /// for [`Boresch::at`]. `None` when no pair has `r₀` in range. See the module documentation for
    /// the rule and what it does not do.
    ///
    /// # Panics
    ///
    /// As [`Boresch::at`].
    pub fn choose(
        at: &[[f64; 3]],
        receptor: &[[usize; 3]],
        ligand: &[[usize; 3]],
        distance: (f64, f64),
        force_constants: [f64; 6],
    ) -> Option<Boresch> {
        let mut best: Option<(f64, Boresch)> = None;
        for &r in receptor {
            for &l in ligand {
                let b = Boresch::at(at, r, l, force_constants);
                let r0 = b.reference[0];
                if r0 < distance.0 || r0 > distance.1 {
                    continue;
                }
                let score = b
                    .hinge_angles(at)
                    .iter()
                    .map(|t| t.sin())
                    .fold(f64::INFINITY, f64::min);
                let better = match &best {
                    None => true,
                    Some((s, _)) => score > *s,
                };
                if better {
                    best = Some((score, b));
                }
            }
        }
        best.map(|(_, b)| b)
    }
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

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
