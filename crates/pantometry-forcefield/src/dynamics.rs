//! Molecular dynamics: Newton's equations for the atoms of a force field, with or without a
//! Langevin heat bath.
//!
//! # The integrator, and why
//!
//! **BAOAB Langevin splitting** (Leimkuhler and Matthews, "Rational construction of stochastic
//! numerical methods for molecular sampling", *Appl. Math. Res. Express* **2013**, 34–56,
//! arXiv:1203.5428 — read there, its Appendix gives the method line by line). Langevin dynamics
//!
//! ```text
//! dq = M⁻¹ p dt,    dp = −∇U dt − γ p dt + √(2γ k_BT) M^½ dW
//! ```
//!
//! is split into a drift `A` (`dq = M⁻¹p dt`), a kick `B` (`dp = −∇U dt`) and an
//! Ornstein–Uhlenbeck part `O` (`dp = −γp dt + noise`), each solved exactly, and composed as
//!
//! ```text
//! B  p ← p + (h/2) F(q)
//! A  q ← q + (h/2) M⁻¹ p
//! O  p ← c₁ p + c₃ M^½ R,      c₁ = e^(−γh),  c₃ = √(k_BT (1 − c₁²))
//! A  q ← q + (h/2) M⁻¹ p
//! B  p ← p + (h/2) F(q)
//! ```
//!
//! with `R` a vector of independent standard normals — the paper's Appendix, symbols and all.
//! One force evaluation per step: the second `B`'s force is the next step's first.
//!
//! - **Why BAOAB rather than another order of the same three pieces.** Of the splittings the paper
//!   compares, BAOAB's invariant measure is the closest to the canonical one *in the positions*,
//!   which is what a binding question samples. Its leading error term (the paper's §3.3)
//!   is `δt²/8 (pᵀU″p − k_BT ΔU)`, and integrating out the momenta cancels it; for a harmonic
//!   potential the positions are sampled **exactly** at every stable step, `⟨q²⟩ = k_BT/k`. That
//!   is the closed form `tests/the_dynamics.rs` checks. The price is that the momenta at a whole
//!   step are not canonical: for a harmonic mode `⟨p²⟩ = m k_BT (1 − (ωh)²/4)` at the step
//!   boundary — exactly, as the same test checks — and `m k_BT` right after `O`. So
//!   [`MolecularDynamics::half_step_kinetic_energy`] is the kinetic energy to average for a temperature,
//!   and [`MolecularDynamics::kinetic_energy`] the one that belongs in the energy books.
//! - **Why not the molecular fluid's integrator.** `pantometry-molecular` applies the bath after a
//!   whole velocity-Verlet step, with its own discretisation of the OU part. That is a fine
//!   thermostat and a worse sampler; and this crate may not depend on it in any case.
//! - **Velocity Verlet is the same code with no bath**: [`Bath::Isolated`] skips `O` and drifts
//!   once by `h` instead of twice by `h/2`, which is BAOAB at γ = 0 without the extra rounding.
//!   It is symplectic and time-reversible, so its energy error stays inside a band of `O(h²)`
//!   instead of drifting: for a harmonic oscillator it conserves `½mv² + ½kx²(1 − (ωh)²/4)`
//!   exactly, so the energy's largest departure from its start at a turning point is
//!   `E₀ (ωh)²/4`. Also checked as a closed form.
//!
//! # Units
//!
//! SI, as the rest of the crate: metres, metres per second, kilograms, newtons, joules per
//! molecule, kelvin, seconds. The masses are the standard atomic weights ([`Element::atomic_weight`])
//! in grams per mole divided by the Avogadro constant.
//!
//! # The time step, measured
//!
//! Every X–H bond is free — no SHAKE or RATTLE — so the step is set by the fastest bond. UFF gives
//! aspirin's C–H `k` = 662 kcal mol⁻¹ Å⁻², a period of 11.5 fs with the C–H reduced mass. NVE
//! aspirin from its relaxed geometry, velocities at 300 K (seed 7, the rigid motion removed), one
//! picosecond at each step (`tests/the_dynamics.rs`, `the_step_size_measured`):
//!
//! | dt (fs) | RMS ΔE (kcal/mol) | max \|ΔE\| | drift (kcal/mol/ps) | ⟨T⟩ (K) |
//! | --- | --- | --- | --- | --- |
//! | 0.125 | 0.0040 | 0.0070 | −0.0003 | 237.9 |
//! | 0.25 | 0.0162 | 0.0281 | −0.0012 | 237.8 |
//! | 0.5 | 0.0654 | 0.1147 | −0.0059 | 237.7 |
//! | 1 | 0.2738 | 0.4468 | +0.0091 | 237.7 |
//! | 1.5 | 0.6637 | 1.0748 | −0.0143 | 237.5 |
//! | 2 | 1.3667 | 2.1629 | −0.0007 | 237.8 |
//! | 2.5 | 2.8027 | 4.5741 | −0.0448 | 237.4 |
//! | 3 | blew up at 27 fs | | | |
//!
//! The RMS error quarters per halving from 0.125 to 1 fs (ratios 4.05, 4.04, 4.19) and departs from
//! `h²` above 2 fs; the stability edge is between 2.5 and 3 fs, below velocity Verlet's `ωh = 2`
//! for that C–H alone (3.7 fs) because bends couple the stretches into faster modes. No step shows
//! a drift beyond its band's noise: what the table measures is the band. **0.5 fs is the default**
//! ([`Molecule::DEFAULT_TIME_STEP`](crate::Molecule::DEFAULT_TIME_STEP)): its worst departure,
//! 0.11 kcal/mol, is a fifth of `k_BT` at 300 K where 1 fs's is three quarters, and BAOAB's
//! whole-step kinetic bias on a C–H stretch, `(ωh)²/4`, is 1.9% against 7.5%. Constraints would
//! allow a step of about 2 fs and are not needed for that. (`⟨T⟩` is the same at every stable
//! step, which is the point of the column; it is not 300 K because the run starts at the minimum
//! with all its energy kinetic. A harmonic molecule would share it out to 150 K; aspirin's torsions
//! and methyl rotor are neither harmonic nor fast, and it reads 238 K. Not asserted.)
//!
//! # Randomness, and why a run is the same however it is cut up
//!
//! Every kick is drawn from [`Rng::for_index`]`(seed ^ KICK_STREAM, step · N + atom)`, with `N` the
//! number of atoms (frozen ones included) and `step` the count of steps this [`MolecularDynamics`] has
//! taken. So the noise atom 7 gets at step 4000 is a function of those three numbers and nothing
//! else: not of how many steps were taken per call, not of which other atoms are frozen, and not of
//! the order the atoms are visited in. A run of 300 steps in one call, three calls of 100 or 300
//! calls of one is the same run, to the bit — `tests/the_dynamics.rs` checks all three.
//! Initial velocities use `Rng::for_index(seed ^ VELOCITY_STREAM, atom)`, a different stream, so
//! a caller who passes one seed to both does not get its first kicks equal to its velocities.
//!
//! **Determinism, and where it stops**, is the crate's: bit-identical on one machine, not across
//! platforms. The Gaussian deviates call the platform's `ln` and `cos` (Box–Muller), `c₁` calls
//! `exp`, and the force field calls more.
//!
//! # Frozen atoms
//!
//! [`MolecularDynamics::with_frozen`] holds atoms where they are, as [`Minimiser::with_frozen`](crate::Minimiser::with_frozen)
//! does: a frozen atom's velocity is zero, it is never drifted, kicked or thermalised, and its
//! position is never written — so it stays put **to the bit**. The force on it is computed and
//! ignored. Its noise index is still reserved, so freezing atom 3 does not change the kicks atom 4
//! receives.
//!
//! # Initial velocities, and which momenta are removed
//!
//! [`MolecularDynamics::thermalised`] draws each free atom's velocity from the Maxwell–Boltzmann
//! distribution at `T`: each component an independent normal of variance `k_BT/m`. Then, **when no
//! atom is frozen**, it removes the centre-of-mass momentum and the angular momentum about the
//! centre of mass:
//!
//! - **Why both, for an isolated molecule.** In vacuum and without a bath, the force field's energy
//!   is invariant under translation and rotation, so velocity Verlet conserves the total linear and
//!   angular momentum. Whatever the draw put in them stays there for the whole run: a molecule that
//!   starts tumbling keeps tumbling, its rotational kinetic energy never reaches the internal modes,
//!   and the temperature read from all `3N` degrees of freedom is not the internal one. Removing
//!   them leaves `3N − 6` (`3N − 5` for a linear molecule, `0` for one atom) degrees of freedom that
//!   exchange energy, which is what [`MolecularDynamics::degrees_of_freedom`] then counts.
//! - **Why the removal does not bias the draw.** Both are orthogonal projections in the
//!   mass-weighted metric — the translations and, about the centre of mass, the rotations at fixed
//!   positions are linear subspaces orthogonal to each other, and `ω = I⁻¹L` is the least-squares
//!   rigid rotation — and the projection of an isotropic Gaussian onto a subspace is the isotropic
//!   Gaussian on it. So what is left is Maxwell–Boltzmann on the remaining degrees of freedom.
//! - **Why not when an atom is frozen.** A frozen atom is an external anchor: the forces it exerts
//!   are not balanced by a recoil, so neither momentum is conserved and there is nothing to remove.
//! - **Under a Langevin bath the removal is harmless and temporary.** The bath acts on every atom
//!   independently, so it does not conserve either momentum: within a few `1/γ` the centre of mass
//!   and the rotations are thermal again, the molecule diffuses and tumbles as a solute in a bath
//!   does, and the count is all `3 N_free` degrees of freedom.
//!
//! # The energy books
//!
//! The `O` step exchanges energy with the bath, and [`MolecularDynamics::thermostat_work`] adds up what it
//! put in: `Σ ½ m (|v′|² − |v|²)` over each `O`, atom by atom. Then `kinetic + potential −
//! thermostat work` changes only by what the deterministic `B` and `A` steps fail to conserve —
//! the integrator's error — and [`MolecularDynamics::ledger`] reports it as `energy` in the convention of
//! `pantometry_core::conserved`: three entries, so the audit judges a change against the largest
//! of them rather than against a sum that may be near zero. The potential is booked from its value
//! at the start ([`MolecularDynamics::potential_reference`]), because a force field's zero is
//! arbitrary and would otherwise be the largest entry, and the scale every change is judged on.

// Three-vectors as `[f64; 3]`, indexed by component as the formulas are written; the Jacobi
// rotations update two columns at once, which an iterator over one cannot.
#![allow(clippy::needless_range_loop)]

use crate::ccd::Element;
use crate::energy::ForceField;
use crate::uff::AVOGADRO;
use pantometry_core::conserved::quantity;
use pantometry_core::{Ledger, Rng};
use pantometry_units::BOLTZMANN;

/// The key mixed into a Langevin seed for the kicks. Any constant other than
/// [`VELOCITY_STREAM`] would do; this is the golden-ratio constant `pantometry-molecular` uses for
/// the same purpose.
pub const KICK_STREAM: u64 = 0x9E37_79B9_7F4A_7C15;

/// The key mixed into a seed for the initial velocities.
pub const VELOCITY_STREAM: u64 = 0xD1B5_4A32_D192_ED03;

impl Element {
    /// The standard atomic weight, the abridged value of the IUPAC Commission on Isotopic
    /// Abundances and Atomic Weights ("Abridged Standard Atomic Weights", ciaaw.org, the 2024
    /// table from the 2021 report: Prohaska et al., *Pure Appl. Chem.* **94**, 573 (2022)).
    /// Dimensionless; it is also the molar mass in grams per mole.
    ///
    /// The abridged values carry an uncertainty of about one in the last digit — 1.0080 ± 0.0002
    /// for hydrogen, 32.06 ± 0.02 for sulfur — which is the natural variation of isotopic
    /// composition, and far larger than anything the dynamics could resolve.
    pub fn atomic_weight(self) -> f64 {
        match self {
            Element::H => 1.0080,
            Element::C => 12.011,
            Element::N => 14.007,
            Element::O => 15.999,
            Element::F => 18.998,
            Element::P => 30.974,
            Element::S => 32.06,
            Element::Cl => 35.45,
            Element::Br => 79.904,
            Element::I => 126.90,
        }
    }

    /// The mass of one atom, kilograms: [`Element::atomic_weight`] grams per mole over the
    /// Avogadro constant. That takes the molar mass constant as exactly 1 g/mol; since 2019 it is
    /// measured, 0.999 999 999 65(30) g/mol (CODATA 2018), a difference of 3.5 × 10⁻¹⁰ against the
    /// weights' own 10⁻⁴.
    pub fn mass(self) -> f64 {
        self.atomic_weight() * 1e-3 / AVOGADRO
    }
}

/// Anything that gives an energy and the force on every atom: what [`MolecularDynamics`] integrates.
///
/// [`ForceField`] is one. A test can be another — an isotropic harmonic well whose dynamics is
/// known in closed form — which is how the integrator is checked apart from the force field.
pub trait Potential {
    /// The energy at positions `at` (metres), joules per molecule, with the force on each atom
    /// written into `forces` (newtons), which has one entry per atom.
    fn energy_and_forces(&self, at: &[[f64; 3]], forces: &mut [[f64; 3]]) -> f64;
}

impl Potential for ForceField {
    fn energy_and_forces(&self, at: &[[f64; 3]], forces: &mut [[f64; 3]]) -> f64 {
        let ev = self.evaluate(at);
        forces.copy_from_slice(&ev.forces);
        ev.energy.total
    }
}

/// What the atoms exchange energy with: nothing, or a heat bath.
///
/// Named `Bath` rather than `Thermostat` because `pantometry-molecular` has a `Thermostat` with the
/// same two ideas and different fields, and the prelude exports that one.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Bath {
    /// Nothing: velocity Verlet, the microcanonical ensemble, and the only setting in which the
    /// energy alone is conserved.
    Isolated,
    /// A Langevin bath, integrated by BAOAB. See [the module documentation](self).
    Langevin {
        /// The bath's temperature, kelvin.
        temperature: f64,
        /// The friction γ, per second. Its inverse is the time the bath takes to forget a
        /// velocity: 10¹² s⁻¹ (1 ps⁻¹) is a common choice for a solute in implicit water, gentle
        /// enough to leave the vibrations alone.
        friction: f64,
        /// The seed every kick is drawn from.
        seed: u64,
    },
}

/// The positions the cached forces were computed at, and those forces.
#[derive(Clone, Debug, PartialEq)]
struct Cache {
    at: Vec<[f64; 3]>,
    forces: Vec<[f64; 3]>,
    energy: f64,
}

/// Molecular dynamics on a [`Potential`], one step at a time: velocities, masses, a bath, and
/// the books. The positions belong to the caller and are passed to every step, as they are to
/// [`Minimiser::step`](crate::Minimiser::step).
///
/// Forces are kept between steps and reused when the positions passed are the ones the last step
/// left, bit for bit; otherwise they are recomputed. **The cache knows the positions, not the
/// potential**: after changing the potential — new charges, a solvent — call
/// [`MolecularDynamics::forget_forces`].
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularDynamics {
    masses: Vec<f64>,
    frozen: Vec<bool>,
    velocities: Vec<[f64; 3]>,
    bath: Bath,
    removed: usize,
    steps: u64,
    thermostat_work: f64,
    half_step_kinetic: f64,
    half_step_kinetic_per_atom: Vec<f64>,
    cache: Option<Cache>,
    reference: Option<f64>,
}

impl MolecularDynamics {
    /// Molecular dynamics of atoms with these masses (kilograms), at rest, in no bath and nothing
    /// frozen.
    ///
    /// # Panics
    ///
    /// If a mass is not positive and finite.
    pub fn new(masses: Vec<f64>) -> MolecularDynamics {
        assert!(
            masses.iter().all(|m| m.is_finite() && *m > 0.0),
            "every mass must be positive and finite"
        );
        let n = masses.len();
        MolecularDynamics {
            masses,
            frozen: vec![false; n],
            velocities: vec![[0.0; 3]; n],
            bath: Bath::Isolated,
            removed: 0,
            steps: 0,
            thermostat_work: 0.0,
            half_step_kinetic: f64::NAN,
            half_step_kinetic_per_atom: vec![f64::NAN; n],
            cache: None,
            reference: None,
        }
    }

    /// Molecular dynamics of atoms of these elements, each with its [`Element::mass`].
    pub fn for_elements(elements: impl IntoIterator<Item = Element>) -> MolecularDynamics {
        MolecularDynamics::new(elements.into_iter().map(Element::mass).collect())
    }

    /// The same dynamics in `bath`.
    pub fn with_bath(mut self, bath: Bath) -> MolecularDynamics {
        if let Bath::Langevin {
            temperature,
            friction,
            ..
        } = bath
        {
            assert!(
                temperature.is_finite() && temperature >= 0.0,
                "a bath's temperature must be finite and not negative"
            );
            assert!(
                friction.is_finite() && friction >= 0.0,
                "a bath's friction must be finite and not negative"
            );
        }
        self.bath = bath;
        self
    }

    /// The same dynamics holding every atom `i` with `frozen[i]` where it is, to the bit, with zero
    /// velocity. See [the module documentation](self).
    ///
    /// # Panics
    ///
    /// If the mask is not one entry per atom.
    pub fn with_frozen(mut self, frozen: Vec<bool>) -> MolecularDynamics {
        assert_eq!(frozen.len(), self.masses.len(), "one frozen flag per atom");
        for (v, &f) in self.velocities.iter_mut().zip(&frozen) {
            if f {
                *v = [0.0; 3];
            }
        }
        self.frozen = frozen;
        self.removed = 0;
        self
    }

    /// The same dynamics starting from these velocities (metres per second), one per atom; a
    /// frozen atom's is set to zero. No momentum is counted as removed.
    ///
    /// # Panics
    ///
    /// If `velocities` is not one per atom.
    pub fn with_velocities(mut self, velocities: Vec<[f64; 3]>) -> MolecularDynamics {
        assert_eq!(velocities.len(), self.masses.len(), "one velocity per atom");
        self.velocities = velocities;
        for (v, &f) in self.velocities.iter_mut().zip(&self.frozen) {
            if f {
                *v = [0.0; 3];
            }
        }
        self.removed = 0;
        self
    }

    /// The same dynamics with Maxwell–Boltzmann velocities at `temperature` kelvin for the atoms
    /// at `at`, drawn from `seed`, and — when no atom is frozen — the centre-of-mass and angular
    /// momentum removed. See [the module documentation](self) for why, and why that leaves the
    /// remaining degrees of freedom Maxwell–Boltzmann.
    ///
    /// # Panics
    ///
    /// If `at` is not one position per atom, or the temperature is negative or not finite.
    pub fn thermalised(
        mut self,
        at: &[[f64; 3]],
        temperature: f64,
        seed: u64,
    ) -> MolecularDynamics {
        assert_eq!(at.len(), self.masses.len(), "one position per atom");
        assert!(
            temperature.is_finite() && temperature >= 0.0,
            "a temperature must be finite and not negative"
        );
        let kt = BOLTZMANN.to_si() * temperature;
        for (i, v) in self.velocities.iter_mut().enumerate() {
            if self.frozen[i] {
                *v = [0.0; 3];
                continue;
            }
            let width = (kt / self.masses[i]).sqrt();
            let mut rng = Rng::for_index(seed ^ VELOCITY_STREAM, i as u64);
            *v = [
                rng.gaussian() * width,
                rng.gaussian() * width,
                rng.gaussian() * width,
            ];
        }
        self.removed = if self.frozen.iter().any(|&f| f) || self.masses.is_empty() {
            0
        } else {
            remove_rigid_motion(&self.masses, at, &mut self.velocities)
        };
        self
    }

    /// Recomputes the forces at the next step whatever the positions, and takes the books' zero of
    /// potential energy again there: for after the potential has changed. See
    /// [`MolecularDynamics::ledger`].
    pub fn forget_forces(&mut self) {
        self.cache = None;
        self.reference = None;
    }

    /// Makes sure the forces are those of `potential` at `at`, computing them if the positions are
    /// not the ones they were computed at. [`MolecularDynamics::step`] does this itself; calling it first
    /// gives [`MolecularDynamics::potential_energy`] and [`MolecularDynamics::ledger`] a value before any step.
    ///
    /// # Panics
    ///
    /// If `at` is not one position per atom.
    pub fn prepare(&mut self, potential: &(impl Potential + ?Sized), at: &[[f64; 3]]) {
        assert_eq!(at.len(), self.masses.len(), "one position per atom");
        if self.cache.as_ref().is_some_and(|c| c.at == at) {
            return;
        }
        let mut forces = vec![[0.0; 3]; at.len()];
        let energy = potential.energy_and_forces(at, &mut forces);
        self.reference.get_or_insert(energy);
        self.cache = Some(Cache {
            at: at.to_vec(),
            forces,
            energy,
        });
    }

    /// One step of `dt` seconds: BAOAB under a Langevin bath, velocity Verlet without one. Moves
    /// `at` in place; frozen atoms are not written.
    ///
    /// # Panics
    ///
    /// If `at` is not one position per atom, or `dt` is not positive and finite.
    pub fn step(&mut self, potential: &(impl Potential + ?Sized), at: &mut [[f64; 3]], dt: f64) {
        assert!(dt.is_finite() && dt > 0.0, "a time step must be positive");
        self.prepare(potential, at);
        let cache = self.cache.take().expect("prepared");
        let half = 0.5 * dt;

        // B: half a kick from the forces where the atoms are.
        self.kick(&cache.forces, half);
        match self.bath {
            Bath::Isolated => self.drift(at, dt),
            Bath::Langevin {
                temperature,
                friction,
                seed,
            } => {
                // A, O, A.
                self.drift(at, half);
                self.exchange_with_bath(temperature, friction, seed, dt);
                self.drift(at, half);
            }
        }
        // The forces where the atoms arrived, which the second B uses and the next step's first
        // B reuses.
        let mut forces = cache.forces;
        let energy = potential.energy_and_forces(at, &mut forces);
        // B: the other half kick.
        self.kick(&forces, half);
        self.cache = Some(Cache {
            at: at.to_vec(),
            forces,
            energy,
        });
        self.steps += 1;
    }

    /// `n` steps of `dt`. The same as calling [`MolecularDynamics::step`] `n` times, to the bit.
    pub fn run(
        &mut self,
        potential: &(impl Potential + ?Sized),
        at: &mut [[f64; 3]],
        dt: f64,
        n: usize,
    ) {
        for _ in 0..n {
            self.step(potential, at, dt);
        }
    }

    fn kick(&mut self, forces: &[[f64; 3]], h: f64) {
        for (i, v) in self.velocities.iter_mut().enumerate() {
            if self.frozen[i] {
                continue;
            }
            let a = h / self.masses[i];
            for c in 0..3 {
                v[c] += a * forces[i][c];
            }
        }
    }

    fn drift(&self, at: &mut [[f64; 3]], h: f64) {
        for (i, x) in at.iter_mut().enumerate() {
            if self.frozen[i] {
                continue;
            }
            for c in 0..3 {
                x[c] += h * self.velocities[i][c];
            }
        }
    }

    /// The exact Ornstein–Uhlenbeck solve over `dt`, adding what it does to the thermostat work.
    fn exchange_with_bath(&mut self, temperature: f64, friction: f64, seed: u64, dt: f64) {
        let c1 = (-friction * dt).exp();
        // 1 − c₁², without the cancellation of forming it from c₁ at small γh.
        let fresh = -(-2.0 * friction * dt).exp_m1();
        let kt = BOLTZMANN.to_si() * temperature;
        let n = self.masses.len() as u64;
        let mut work = 0.0;
        let mut kinetic = 0.0;
        for (i, v) in self.velocities.iter_mut().enumerate() {
            if self.frozen[i] {
                self.half_step_kinetic_per_atom[i] = 0.0;
                continue;
            }
            let m = self.masses[i];
            let mut rng = Rng::for_index(seed ^ KICK_STREAM, self.steps * n + i as u64);
            let width = (fresh * kt / m).sqrt();
            let before = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            for c in v.iter_mut() {
                *c = c1 * *c + width * rng.gaussian();
            }
            let after = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            work += 0.5 * m * (after - before);
            kinetic += 0.5 * m * after;
            self.half_step_kinetic_per_atom[i] = 0.5 * m * after;
        }
        self.thermostat_work += work;
        self.half_step_kinetic = kinetic;
    }

    /// The masses, kilograms.
    pub fn masses(&self) -> &[f64] {
        &self.masses
    }

    /// The frozen-atom mask, one entry per atom.
    pub fn frozen(&self) -> &[bool] {
        &self.frozen
    }

    /// The velocities, metres per second, at the last whole step.
    pub fn velocities(&self) -> &[[f64; 3]] {
        &self.velocities
    }

    /// The bath.
    pub fn bath(&self) -> Bath {
        self.bath
    }

    /// How many steps have been taken.
    pub fn steps(&self) -> u64 {
        self.steps
    }

    /// The kinetic energy at the last whole step, joules per molecule: the one that pairs with the
    /// potential energy at the same positions in the books.
    pub fn kinetic_energy(&self) -> f64 {
        self.velocities
            .iter()
            .zip(&self.masses)
            .map(|(v, m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum()
    }

    /// The kinetic energy right after the last step's `O`, joules per molecule; `NaN` without a
    /// Langevin bath or before the first step. BAOAB's momenta are canonical there — exactly, for a
    /// harmonic potential — and not at the whole step, so this is the one to average for a
    /// temperature. See [the module documentation](self).
    pub fn half_step_kinetic_energy(&self) -> f64 {
        self.half_step_kinetic
    }

    /// The same atom by atom, joules per molecule, in the atoms' order; zero for a frozen atom, and
    /// `NaN` without a Langevin bath or before the first step. What equipartition is checked on
    /// per element: every atom, whatever its mass, has `⟨½ m v²⟩ = (3/2) k_BT`.
    pub fn half_step_kinetic_energies(&self) -> &[f64] {
        &self.half_step_kinetic_per_atom
    }

    /// The potential energy at the positions of the last step or [`MolecularDynamics::prepare`], joules per
    /// molecule; `None` before either.
    pub fn potential_energy(&self) -> Option<f64> {
        self.cache.as_ref().map(|c| c.energy)
    }

    /// The energy the bath has put in since the start, joules per molecule; negative when it has
    /// taken energy out.
    pub fn thermostat_work(&self) -> f64 {
        self.thermostat_work
    }

    /// The degrees of freedom that exchange energy: three per free atom, less the momenta
    /// [`MolecularDynamics::thermalised`] removed — which only the microcanonical dynamics keeps removed. A
    /// Langevin bath re-thermalises them, so under one this is `3 N_free`.
    pub fn degrees_of_freedom(&self) -> usize {
        let free = 3 * self.frozen.iter().filter(|&&f| !f).count();
        match self.bath {
            Bath::Isolated => free - self.removed.min(free),
            Bath::Langevin { .. } => free,
        }
    }

    /// `2 KE / (n k_B)`, kelvin, from the kinetic energy at the whole step and
    /// [`MolecularDynamics::degrees_of_freedom`]; `NaN` with none.
    pub fn temperature(&self) -> f64 {
        self.temperature_of(self.kinetic_energy())
    }

    /// The same from [`MolecularDynamics::half_step_kinetic_energy`].
    pub fn half_step_temperature(&self) -> f64 {
        self.temperature_of(self.half_step_kinetic)
    }

    fn temperature_of(&self, kinetic: f64) -> f64 {
        let n = self.degrees_of_freedom();
        if n == 0 {
            return f64::NAN;
        }
        2.0 * kinetic / (n as f64 * BOLTZMANN.to_si())
    }

    /// The potential energy the books count from, joules per molecule: the potential at the first
    /// evaluation, or the first after [`MolecularDynamics::forget_forces`]; `None` before it.
    pub fn potential_reference(&self) -> Option<f64> {
        self.reference
    }

    /// `energy`: kinetic + (potential − [`MolecularDynamics::potential_reference`]) − thermostat
    /// work, as three entries; empty before the first step or [`MolecularDynamics::prepare`], since
    /// there is no potential energy to report.
    ///
    /// **Why the potential is counted from its start.** A force field's zero is arbitrary — UFF's
    /// aspirin sits at 29.6 kcal/mol at its minimum — and an audit judges a change against the
    /// largest entry. Booking the potential whole made that zero part of the scale. It let a bath
    /// whose work was half counted pass a 1e-2 audit, because the error was small against 30
    /// kcal/mol that never moves. Counted from the start, the entries are the energy that moves:
    /// the kinetic energy, the potential's change and the work.
    pub fn ledger(&self) -> Ledger {
        match (self.potential_energy(), self.reference) {
            (Some(potential), Some(reference)) => Ledger::new()
                .with(quantity::ENERGY, self.kinetic_energy())
                .with(quantity::ENERGY, potential - reference)
                .with(quantity::ENERGY, -self.thermostat_work),
            _ => Ledger::new(),
        }
    }

    /// The total linear momentum, kg m s⁻¹.
    pub fn momentum(&self) -> [f64; 3] {
        let mut p = [0.0; 3];
        for (v, m) in self.velocities.iter().zip(&self.masses) {
            for c in 0..3 {
                p[c] += m * v[c];
            }
        }
        p
    }

    /// The angular momentum about the centre of mass of the atoms at `at`, kg m² s⁻¹.
    ///
    /// # Panics
    ///
    /// If `at` is not one position per atom.
    pub fn angular_momentum(&self, at: &[[f64; 3]]) -> [f64; 3] {
        assert_eq!(at.len(), self.masses.len(), "one position per atom");
        let centre = centre_of_mass(&self.masses, at);
        let mut l = [0.0; 3];
        for ((x, v), m) in at.iter().zip(&self.velocities).zip(&self.masses) {
            let d = sub(*x, centre);
            let c = cross(d, *v);
            for k in 0..3 {
                l[k] += m * c[k];
            }
        }
        l
    }
}

fn centre_of_mass(masses: &[f64], at: &[[f64; 3]]) -> [f64; 3] {
    let total: f64 = masses.iter().sum();
    let mut r = [0.0; 3];
    for (x, m) in at.iter().zip(masses) {
        for c in 0..3 {
            r[c] += m * x[c];
        }
    }
    [r[0] / total, r[1] / total, r[2] / total]
}

/// Removes the centre-of-mass momentum and the angular momentum about the centre of mass, and
/// returns how many degrees of freedom that took: 3, plus the rank of the inertia tensor (3 for a
/// non-linear molecule, 2 for a linear one, 0 for one atom).
fn remove_rigid_motion(masses: &[f64], at: &[[f64; 3]], v: &mut [[f64; 3]]) -> usize {
    let total: f64 = masses.iter().sum();
    let mut p = [0.0; 3];
    for (vi, m) in v.iter().zip(masses) {
        for c in 0..3 {
            p[c] += m * vi[c];
        }
    }
    for vi in v.iter_mut() {
        for c in 0..3 {
            vi[c] -= p[c] / total;
        }
    }
    let centre = centre_of_mass(masses, at);
    let mut l = [0.0; 3];
    let mut inertia = [[0.0; 3]; 3];
    for ((x, vi), m) in at.iter().zip(v.iter()).zip(masses) {
        let d = sub(*x, centre);
        let c = cross(d, *vi);
        let dd = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        for a in 0..3 {
            l[a] += m * c[a];
            for b in 0..3 {
                let delta = if a == b { dd } else { 0.0 };
                inertia[a][b] += m * (delta - d[a] * d[b]);
            }
        }
    }
    let (values, vectors) = symmetric_eigen(inertia);
    let largest = values.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    // ω = I⁺ L through the eigen-decomposition, leaving out the axes with no moment: a linear
    // molecule's own axis, or every axis of a single atom. 1e-10 of the largest moment is far
    // above the rounding of an exactly zero one and far below any molecule's smallest.
    let mut omega = [0.0; 3];
    let mut rank = 0;
    for k in 0..3 {
        if values[k] > 1e-10 * largest && largest > 0.0 {
            rank += 1;
            let e = [vectors[0][k], vectors[1][k], vectors[2][k]];
            let along = (e[0] * l[0] + e[1] * l[1] + e[2] * l[2]) / values[k];
            for c in 0..3 {
                omega[c] += along * e[c];
            }
        }
    }
    for (x, vi) in at.iter().zip(v.iter_mut()) {
        let w = cross(omega, sub(*x, centre));
        for c in 0..3 {
            vi[c] -= w[c];
        }
    }
    3 + rank
}

/// The eigenvalues and eigenvectors (as columns) of a symmetric 3×3 matrix, by cyclic Jacobi
/// rotations. Deterministic: a fixed sweep order, stopping when the off-diagonal part is exactly
/// zero or after fifty sweeps.
fn symmetric_eigen(mut a: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..50 {
        let off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        if off == 0.0 {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if a[p][q] == 0.0 {
                continue;
            }
            let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            for k in 0..3 {
                let (akp, akq) = (a[k][p], a[k][q]);
                a[k][p] = c * akp - s * akq;
                a[k][q] = s * akp + c * akq;
            }
            for k in 0..3 {
                let (apk, aqk) = (a[p][k], a[q][k]);
                a[p][k] = c * apk - s * aqk;
                a[q][k] = s * apk + c * aqk;
            }
            for k in 0..3 {
                let (vkp, vkq) = (v[k][p], v[k][q]);
                v[k][p] = c * vkp - s * vkq;
                v[k][q] = s * vkp + c * vkq;
            }
        }
    }
    ([a[0][0], a[1][1], a[2][2]], v)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The Jacobi eigen-decomposition reconstructs its matrix, `A = V Λ Vᵀ`, with `V` orthonormal,
    /// on a full matrix and on a rank-two one (a linear molecule's inertia tensor).
    #[test]
    fn the_eigen_decomposition_reconstructs_the_matrix() {
        for a in [
            [[4.0, 1.0, -2.0], [1.0, 3.0, 0.5], [-2.0, 0.5, 6.0]],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            [[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]],
        ] {
            let (l, v) = symmetric_eigen(a);
            for i in 0..3 {
                for j in 0..3 {
                    let r: f64 = (0..3).map(|k| v[i][k] * l[k] * v[j][k]).sum();
                    assert!((r - a[i][j]).abs() < 1e-14, "{a:?}");
                    let o: f64 = (0..3).map(|k| v[k][i] * v[k][j]).sum();
                    let id = if i == j { 1.0 } else { 0.0 };
                    assert!((o - id).abs() < 1e-14, "{a:?}");
                }
            }
        }
    }
}
