//! Bond-length constraints for any molecule: SHAKE for the positions and RATTLE for the
//! velocities, each step of the iteration solved exactly for its one bond.
//!
//! This is step W3 of the explicit-water track. W2 held water rigid with SETTLE
//! ([`crate::water::Settle`]), which solves one three-site triangle in closed form and nothing
//! else; a UFF solute's X–H bonds are free, and their period of about 11 fs is what held this
//! crate's step to 0.5 fs ([`crate::dynamics`], "The time step, measured"). Holding them lets
//! the solute run at the 2 fs the water does.
//!
//! # The method
//!
//! Ryckaert, Ciccotti and Berendsen, *J. Comput. Phys.* **23**, 327 (1977) for SHAKE, and
//! Andersen, *J. Comput. Phys.* **52**, 24 (1983) for RATTLE — **cited for the method, not
//! opened**; what follows is each condition's own requirement, derived here, and the tests check
//! the result rather than a transcription. A constraint holds atoms `i` and `j` at distance `d`.
//!
//! - **Positions.** After an unconstrained drift from positions `x⁰` that satisfied every
//!   constraint, the constraint impulse on a bond lies along that bond's **old** vector
//!   `r⁰ = x⁰_i − x⁰_j`, equal and opposite on its two atoms and inversely weighted by their
//!   masses: `x_i += g r⁰/m_i`, `x_j −= g r⁰/m_j`. With `s = x_i − x_j` the current vector,
//!   `w = 1/m_i + 1/m_j` and `μ = g w`, the bond is satisfied when
//!   `|r⁰|² μ² + 2 (s·r⁰) μ + (|s|² − d²) = 0`. **SHAKE linearises this and iterates; here each
//!   bond's quadratic is solved exactly**, by the root of smaller magnitude written without
//!   cancellation, `μ = −c / (b + sgn(b) √(b² − a c))` with `a = |r⁰|²`, `b = s·r⁰`, `c = |s|² − d²`.
//!   Bonds are visited in order and the sweep repeated (Gauss–Seidel, as SHAKE does) until every
//!   bond is within the tolerance. **A set of bonds that share no atom — every X–H bond of a
//!   molecule in which no heavy atom carries two hydrogens, benzene's six among them — is then
//!   solved in one sweep, to rounding**: each bond's update touches no other bond's atoms.
//! - **Velocities.** RATTLE's second half removes each bond's relative velocity along it:
//!   with `r = x_i − x_j` at the new positions and `u = v_i − v_j`, `κ = (r·u)/(w |r|²)`,
//!   `v_i −= κ r/m_i`, `v_j += κ r/m_j`, after which `r·(v_i − v_j) = 0`. Iterated over the bonds the
//!   same way, this is Gauss–Seidel on the mass-weighted normal equations, whose fixed point is
//!   the mass-weighted orthogonal projection onto the velocities that keep every bond's length —
//!   the same projection SETTLE's [`Settle::constrain_velocities`](crate::water::Settle::constrain_velocities)
//!   makes in closed form, and so it takes an isotropic Gaussian to Maxwell–Boltzmann on what is
//!   left. One sweep, exactly, for bonds that share no atom.
//! - **The tolerance**, [`Shake::DEFAULT_TOLERANCE`] = 10⁻¹²: a bond is converged when
//!   `|r² − d²| ≤ 2 tol d²`, its length within about `tol d` of `d`, and its velocity when
//!   `|r·u| ≤ tol |r| |u|`. Rounding alone leaves a bond in a box 25 Å across about `25 ε` ≈ 6 × 10⁻¹⁵
//!   of its length from `d`, so the tolerance is two hundred times above the floor it could not get
//!   below. A sweep corrects only the bonds outside it; a sweep that corrects none ends the
//!   iteration. [`Shake::MAX_SWEEPS`] without convergence panics, naming the worst bond.
//! - **Frozen atoms** are atoms of infinite mass: their inverse mass is taken as zero, so a
//!   constraint between a frozen and a free atom holds the free one at its distance from the
//!   frozen one, and moves nothing frozen. A constraint between two frozen atoms is skipped. SETTLE
//!   refuses a molecule frozen in part; SHAKE has no reason to.
//!
//! **The constrained ensemble needs no metric correction here.** For a set of bond-length
//! constraints that share no atom, the mass-weighted Gram matrix of their gradients is diagonal
//! with entries `1/m_i + 1/m_j`, a constant, so Fixman's determinant is a constant and the
//! constrained dynamics samples the rigid-bond Boltzmann distribution as it stands. For bonds that
//! share an atom it is not constant in general; it is not computed, and a free energy between two
//! states with the same constraints does not need it if it cancels, which it does when the
//! constrained atoms are the same in both and their environment does not change their geometry —
//! the case of an alchemical solute's own bonds.
//!
//! # Inside the integrator
//!
//! [`crate::MolecularDynamics::with_bond_constraints`] composes it with BAOAB exactly as SETTLE is
//! composed: after each drift, SHAKE and the velocity correction `Δq/h` it implies, then the
//! projection; after each kick and after `O`, the projection. **Each constraint removes one degree
//! of freedom** from [`crate::MolecularDynamics::degrees_of_freedom`].
//!
//! # What it was for: benzene's hydration free energy in TIP3P
//!
//! `tests/benzene_hydrated_in_tip3p.rs`, ignored, release, one core, 9.3 h. Benzene with its six
//! C–H held and 3c-3's QEq charges (±0.098 e), in 507 rigid TIP3P waters in a 24.94 Å cube — the
//! model's own one-atmosphere density, 33.00 nm⁻³ — `r_c` = 9 Å, Ewald at δ = 10⁻⁵; 2 fs in a
//! 1 ps⁻¹ bath at 298.15 K; thirteen windows of 20 ps discarded and 200 ps sampled every 100 fs,
//! each recording its energy at all 29 candidate states. Charges off at full van der Waals, then
//! van der Waals off through the soft core. **No window was inserted**: the lowest overlap was
//! 0.225, against 3c's threshold of 0.1. `ΔG_hyd = −ΔG_decouple`, with no vacuum leg, because the
//! decoupled state is benzene alone in vacuum with its intramolecular terms
//! ([`crate::PeriodicDecoupling`]).
//!
//! | | kcal/mol |
//! | --- | --- |
//! | charges off (BAR; TI +1.278) | +1.227 ± 0.038 |
//! | van der Waals off (BAR; TI +0.400), the long-range correction's +0.915 included | +0.343 ± 0.135 |
//! | decoupling (BAR; TI +1.678 ± 0.142) | +1.570 ± 0.140 |
//! | **ΔG_hydration** | **−1.57 ± 0.14** |
//! | GAFF/AM1-BCC in TIP3P, FreeSolv v0.52's calculated value | −0.81 ± 0.02 |
//! | OBC II with the same charges (3c-3): polar, and with its nonpolar term | −2.34, −1.13 |
//! | experiment, FreeSolv v0.52 | −0.90 ± 0.20 |
//!
//! **0.67 too negative, 2.7 of the combined errors.** Nothing is asserted against experiment.
//! The halves of every window give +1.38 ± 0.21 and +1.75 ± 0.18, 1.4σ apart, and TI is 0.11
//! above BAR, mostly the trapezoid's bias where `⟨∂U/∂λ_v⟩` turns over between 0.4 and 0.1. The
//! charges' −1.23 is about half of OBC II's −2.34 for the same charges: GB over-solvates them, as
//! 2c measured for small molecules. What limits the comparison: UFF's van der Waals mixed with
//! TIP3P's by UFF's geometric rule, which no one fitted (a UFF hydrogen's well is 0.044 kcal/mol);
//! QEq's charges, which are not the AM1-BCC charges TIP3P is usually paired with; a 200 ps window;
//! and a 507-water box, whose finite-size terms for a neutral solute are measured as small — its
//! own images −0.008 kcal/mol, the tinfoil dipole term 2 × 10⁻⁵ — while its density was not
//! re-equilibrated under a barostat.

use crate::ccd::{Component, Element};
use crate::energy::ForceField;

/// Bond-length constraints, held by SHAKE and RATTLE: see the module documentation.
#[derive(Clone, Debug, PartialEq)]
pub struct Shake {
    bonds: Vec<[usize; 2]>,
    lengths: Vec<f64>,
    tolerance: f64,
    independent: bool,
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

impl Shake {
    /// A bond is converged within this relative tolerance: see the module documentation.
    pub const DEFAULT_TOLERANCE: f64 = 1e-12;

    /// Sweeps over the bonds before the iteration gives up and panics.
    pub const MAX_SWEEPS: usize = 1000;

    /// Atoms `bonds[k]` held at distance `lengths[k]` metres, at [`Shake::DEFAULT_TOLERANCE`].
    ///
    /// # Panics
    ///
    /// If the two lists differ in length, a length is not positive and finite, a bond joins an
    /// atom to itself, or one pair of atoms is constrained twice.
    pub fn new(bonds: Vec<[usize; 2]>, lengths: Vec<f64>) -> Shake {
        assert_eq!(
            bonds.len(),
            lengths.len(),
            "one length per constrained bond"
        );
        assert!(
            lengths.iter().all(|d| d.is_finite() && *d > 0.0),
            "a constrained length must be positive and finite"
        );
        assert!(
            bonds.iter().all(|b| b[0] != b[1]),
            "a constraint joins two atoms"
        );
        let mut pairs: Vec<[usize; 2]> = bonds
            .iter()
            .map(|b| [b[0].min(b[1]), b[0].max(b[1])])
            .collect();
        pairs.sort_unstable();
        assert!(
            pairs.windows(2).all(|w| w[0] != w[1]),
            "a pair of atoms is constrained twice"
        );
        let mut atoms: Vec<usize> = bonds.iter().flatten().copied().collect();
        atoms.sort_unstable();
        let independent = atoms.windows(2).all(|w| w[0] != w[1]);
        Shake {
            bonds,
            lengths,
            tolerance: Shake::DEFAULT_TOLERANCE,
            independent,
        }
    }

    /// Every bond of `component` with a hydrogen at one end, held at its UFF natural length
    /// `r_IJ` — the rest length of `field`'s stretch for it, the one its bond term would relax
    /// to. `field` must be `component`'s ([`ForceField::new`], or a periodic field's
    /// [`bonded`](crate::PeriodicForceField::bonded) terms); its stretches are in the component's
    /// bond order.
    ///
    /// # Panics
    ///
    /// If `field` does not have one stretch per bond of `component`, on the same atoms.
    pub fn to_hydrogen(component: &Component, field: &ForceField) -> Shake {
        let stretches = field.stretches();
        assert_eq!(
            stretches.len(),
            component.bonds().len(),
            "one stretch per bond: the force field must be the component's"
        );
        let atoms = component.atoms();
        let (mut bonds, mut lengths) = (Vec::new(), Vec::new());
        for (b, s) in component.bonds().iter().zip(stretches) {
            assert_eq!(b.atoms, s.atoms, "the force field must be the component's");
            let [i, j] = b.atoms;
            if atoms[i].element == Element::H || atoms[j].element == Element::H {
                bonds.push(b.atoms);
                lengths.push(s.natural_length);
            }
        }
        Shake::new(bonds, lengths)
    }

    /// The same constraints converged to `tolerance` instead.
    ///
    /// # Panics
    ///
    /// If `tolerance` is not positive and finite.
    pub fn with_tolerance(mut self, tolerance: f64) -> Shake {
        assert!(
            tolerance.is_finite() && tolerance > 0.0,
            "a tolerance must be positive and finite"
        );
        self.tolerance = tolerance;
        self
    }

    /// The constrained pairs, atom indices.
    pub fn bonds(&self) -> &[[usize; 2]] {
        &self.bonds
    }

    /// Each pair's length, metres.
    pub fn lengths(&self) -> &[f64] {
        &self.lengths
    }

    /// The relative tolerance.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// How many constraints.
    pub fn len(&self) -> usize {
        self.bonds.len()
    }

    /// Whether there are none.
    pub fn is_empty(&self) -> bool {
        self.bonds.is_empty()
    }

    /// Whether no atom is in two constraints, so that one sweep solves them exactly.
    pub fn is_independent(&self) -> bool {
        self.independent
    }

    /// Each atom's inverse mass, zero where `frozen`.
    fn inverse_masses(masses: &[f64], frozen: &[bool]) -> Vec<f64> {
        masses
            .iter()
            .zip(frozen)
            .map(|(m, &f)| if f { 0.0 } else { 1.0 / m })
            .collect()
    }

    /// Moves `after` — unconstrained positions one drift on from `before`, where every
    /// constraint held — onto the constraints, by SHAKE with each bond's quadratic solved exactly.
    /// Atoms `frozen` marks are not moved. Returns the sweeps taken.
    ///
    /// # Panics
    ///
    /// If the lengths of the slices disagree or an atom is past them, if a bond has turned so far
    /// that its quadratic has no root, or after [`Shake::MAX_SWEEPS`] sweeps without convergence.
    pub fn constrain_positions(
        &self,
        masses: &[f64],
        frozen: &[bool],
        before: &[[f64; 3]],
        after: &mut [[f64; 3]],
    ) -> usize {
        let n = masses.len();
        assert!(
            frozen.len() == n && before.len() == n && after.len() == n,
            "one mass, mark and position per atom"
        );
        let inv = Shake::inverse_masses(masses, frozen);
        for sweep in 1..=Shake::MAX_SWEEPS {
            let mut moved = false;
            for (k, &[i, j]) in self.bonds.iter().enumerate() {
                let w = inv[i] + inv[j];
                if w == 0.0 {
                    continue;
                }
                let d2 = self.lengths[k] * self.lengths[k];
                let s = sub(after[i], after[j]);
                let c = dot(s, s) - d2;
                assert!(
                    c.is_finite(),
                    "SHAKE was given a position that is not finite, on the bond of atoms {i} and {j}"
                );
                if c.abs() <= 2.0 * self.tolerance * d2 {
                    continue;
                }
                let r0 = sub(before[i], before[j]);
                let (a, b) = (dot(r0, r0), dot(s, r0));
                let disc = b * b - a * c;
                let q = b + (disc.max(0.0)).sqrt().copysign(b);
                assert!(
                    disc >= 0.0 && q != 0.0,
                    "SHAKE has no solution for the bond of atoms {i} and {j}: it turned too far \
                     in one step"
                );
                let mu = -c / q;
                for x in 0..3 {
                    after[i][x] += mu * inv[i] / w * r0[x];
                    after[j][x] -= mu * inv[j] / w * r0[x];
                }
                moved = true;
            }
            if !moved {
                return sweep;
            }
        }
        let (k, worst) = self.worst_length(after);
        panic!(
            "SHAKE did not converge in {} sweeps: the bond of atoms {:?} is {worst:e} of its \
             length off",
            Shake::MAX_SWEEPS,
            self.bonds[k]
        );
    }

    /// The bond furthest from its length, and how far, relative.
    fn worst_length(&self, at: &[[f64; 3]]) -> (usize, f64) {
        self.bonds
            .iter()
            .enumerate()
            .map(|(k, &[i, j])| {
                let s = sub(at[i], at[j]);
                (
                    k,
                    (dot(s, s).sqrt() - self.lengths[k]).abs() / self.lengths[k],
                )
            })
            .fold((0, 0.0), |a, b| if b.1 > a.1 { b } else { a })
    }

    /// Removes from `velocities` every component along a constraint at positions `at`, by RATTLE's
    /// projection: afterwards `(r_i − r_j)·(v_i − v_j) = 0` for every constraint, to the tolerance.
    /// Atoms `frozen` marks are not touched. Returns the sweeps taken.
    ///
    /// # Panics
    ///
    /// If the lengths of the slices disagree or an atom is past them, or after
    /// [`Shake::MAX_SWEEPS`] sweeps without convergence.
    pub fn constrain_velocities(
        &self,
        masses: &[f64],
        frozen: &[bool],
        at: &[[f64; 3]],
        velocities: &mut [[f64; 3]],
    ) -> usize {
        let n = masses.len();
        assert!(
            frozen.len() == n && at.len() == n && velocities.len() == n,
            "one mass, mark, position and velocity per atom"
        );
        let inv = Shake::inverse_masses(masses, frozen);
        let v = velocities;
        for sweep in 1..=Shake::MAX_SWEEPS {
            let mut moved = false;
            for &[i, j] in &self.bonds {
                let w = inv[i] + inv[j];
                if w == 0.0 {
                    continue;
                }
                let r = sub(at[i], at[j]);
                let u = sub(v[i], v[j]);
                let along = dot(r, u);
                assert!(
                    along.is_finite(),
                    "RATTLE was given a position or velocity that is not finite, on the bond of atoms \
                     {i} and {j}"
                );
                let (r2, u2) = (dot(r, r), dot(u, u));
                if along * along <= self.tolerance * self.tolerance * r2 * u2 {
                    continue;
                }
                let kappa = along / (w * r2);
                for x in 0..3 {
                    v[i][x] -= kappa * inv[i] * r[x];
                    v[j][x] += kappa * inv[j] * r[x];
                }
                moved = true;
            }
            if !moved {
                return sweep;
            }
        }
        panic!(
            "RATTLE's velocities did not converge in {} sweeps",
            Shake::MAX_SWEEPS
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One bond, its quadratic solved exactly: the length to rounding in one sweep, then a sweep
    /// that moves nothing; the centre of mass where it was; and the velocity along the bond gone.
    #[test]
    fn one_bond_is_solved_in_one_sweep() {
        let masses = [12.0e-27, 1.0e-27];
        let frozen = [false, false];
        let before = [[0.0, 0.0, 0.0], [1.0e-10, 0.0, 0.0]];
        let mut after = [
            [0.01e-10, 0.02e-10, -0.01e-10],
            [1.05e-10, 0.1e-10, 0.03e-10],
        ];
        let s = Shake::new(vec![[0, 1]], vec![1.0e-10]);
        let com = |p: &[[f64; 3]; 2]| {
            [0, 1, 2].map(|k| (masses[0] * p[0][k] + masses[1] * p[1][k]) / (masses[0] + masses[1]))
        };
        let c0 = com(&after);
        assert_eq!(
            s.constrain_positions(&masses, &frozen, &before, &mut after),
            2
        );
        let d = sub(after[0], after[1]);
        assert!((dot(d, d).sqrt() - 1.0e-10).abs() <= 4.0 * f64::EPSILON * 1.0e-10);
        let c1 = com(&after);
        for k in 0..3 {
            assert!((c1[k] - c0[k]).abs() <= 4.0 * f64::EPSILON * 1.1e-10);
        }
        let mut v = [[3.0, -1.0, 2.0], [-40.0, 7.0, 1.0]];
        let u0 = sub(v[0], v[1]);
        s.constrain_velocities(&masses, &frozen, &after, &mut v);
        let u = sub(v[0], v[1]);
        // Rounding of the relative velocity before the projection, against the bond.
        assert!(dot(d, u).abs() <= 4.0 * f64::EPSILON * dot(d, d).sqrt() * dot(u0, u0).sqrt());
    }

    /// A velocity that is not finite is refused by name, rather than iterated until the sweeps run
    /// out and reported as a failure to converge.
    #[test]
    #[should_panic(expected = "RATTLE was given a position or velocity that is not finite")]
    fn a_velocity_that_is_not_finite_is_refused_by_name() {
        let s = Shake::new(vec![[0, 1]], vec![1.0e-10]);
        let at = [[0.0, 0.0, 0.0], [1.0e-10, 0.0, 0.0]];
        let mut v = [[f64::NAN, 0.0, 0.0], [0.0, 0.0, 0.0]];
        s.constrain_velocities(&[1e-26, 1e-27], &[false, false], &at, &mut v);
    }

    /// Three bonds that share atoms (a bent chain and its closing bond, a triangle) converge to the
    /// tolerance by iteration, and the counts say they are not independent.
    #[test]
    fn coupled_bonds_converge_to_the_tolerance() {
        let masses = [16.0e-27, 1.0e-27, 1.0e-27];
        let frozen = [false; 3];
        let before = [[0.0, 0.0, 0.0], [1.0e-10, 0.0, 0.0], [0.0, 1.0e-10, 0.0]];
        let lengths = vec![1.0e-10, 1.0e-10, 2f64.sqrt() * 1.0e-10];
        let s = Shake::new(vec![[0, 1], [0, 2], [1, 2]], lengths.clone());
        assert!(!s.is_independent());
        let mut after = [
            [0.02e-10, -0.01e-10, 0.0],
            [1.03e-10, 0.05e-10, 0.02e-10],
            [-0.04e-10, 1.02e-10, -0.03e-10],
        ];
        let sweeps = s.constrain_positions(&masses, &frozen, &before, &mut after);
        assert!(sweeps > 2, "{sweeps}");
        for (k, &[i, j]) in s.bonds().iter().enumerate() {
            let d = sub(after[i], after[j]);
            assert!((dot(d, d).sqrt() - lengths[k]).abs() <= 1.01e-12 * lengths[k]);
        }
    }
}
