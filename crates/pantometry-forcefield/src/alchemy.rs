//! A λ-dependent Hamiltonian that decouples one group of atoms — a ligand — from the rest:
//! electrostatics scaled linearly, van der Waals through a soft core, and a Boresch restraint
//! switched on linearly, each with its own λ. What [`crate::free_energy`] samples and integrates.
//!
//! # The states
//!
//! A [`Lambda`] has three components, each in [0, 1]: `restraint`, `electrostatics` and
//! `van_der_waals`. **`electrostatics = van_der_waals = 1` is the ligand fully coupled, and 0 is
//! decoupled; `restraint = 0` is off and 1 on.** Double decoupling (below) walks them one at a
//! time.
//!
//! # The Hamiltonian
//!
//! `U(x; λ) = U_rest(x) + λ_r U_B(x) + Σ_cross [λ_e C_ij(r) + U_sc(r; λ_v)]`
//!
//! - **`U_rest`** is the force field with the non-bonded pairs between the group and its
//!   environment taken out. **The group's own interactions are kept**, bonded and non-bonded,
//!   at every λ: this is **decoupling**, not annihilation. Decoupling is what a binding free
//!   energy by double decoupling needs, because the ligand's intramolecular energy is then the
//!   same in the decoupled state of both legs — in the complex and in solvent — and cancels from
//!   the cycle, so neither leg has to be converged against the ligand's internal conformations at
//!   λ = 0. It is what Mobley, Chodera and Dill (*J. Chem. Phys.* **125**, 084902 (2006)) did —
//!   "only the Lennard-Jones interactions of the ligand with its environment were eliminated", as
//!   read in its PMC text through a summary — and Mobley et al. 2007 on this very system
//!   (*J. Mol. Biol.* **371**, 1118, PMC2104542) uses the intermediate states of that work. In
//!   vacuum, which is all this module does yet, the choice matters less: an annihilated ligand's
//!   intramolecular Coulomb term is a constant the two legs would share too, but its van der Waals
//!   term is not, and annihilating it lets the ring's own atoms overlap.
//! - **`C_ij = 332.0637 q_i q_j / r`** ([`coulomb`]) for every pair across the partition, scaled
//!   linearly by `λ_e`, so `∂U/∂λ_e = Σ C_ij` exactly: the coupled cross Coulomb energy at every
//!   λ. Linear because a charge is turned off with the van der Waals still on — the protocol below
//!   — so the pair never comes near `r = 0` while it carries a charge. That is GROMACS's rule
//!   ("the Coulombic interaction is turned off linearly, rather than using soft-core
//!   interactions", when the two are switched sequentially; manual, "Free energy interactions").
//! - **`U_sc`, the soft-core Lennard-Jones**, for every pair across the partition, at `λ_v`:
//!
//!   `U_sc(r; λ) = λ · V_LJ(r_sc)`,   `r_sc⁶ = α σ⁶ (1 − λ)^p + r⁶`,
//!
//!   with `V_LJ` the pair's own UFF form `D[(x/r)¹² − 2(x/r)⁶]` and `σ⁶ = x⁶/2`, the zero of that
//!   form (`(C₁₂/C₆)^⅙` in GROMACS's words: `C₁₂ = D x¹²`, `C₆ = 2D x⁶`). **α = 0.5 and p = 1**
//!   ([`SoftCore::default`]).
//! - **`U_B`**, a [`Boresch`] restraint, scaled linearly by `λ_r`; none unless given.
//!
//! **The soft core, and where it was read.** The functional form is Beutler, Mark, van
//! Schaik, Gerber and van Gunsteren, *Chem. Phys. Lett.* **222**, 529 (1994), **which was not
//! opened here**: it is behind the publisher's paywall. The form as written above is the GROMACS
//! reference manual's free-energy section (secondary; read at manual.gromacs.org,
//! "Free energy interactions"), which attributes it to Beutler et al.:
//! `V_sc(r) = (1 − λ) V^A(r_A) + λ V^B(r_B)`, `r_B = (α σ_B⁶ (1 − λ)^p + r⁶)^⅙`, with state A the decoupled one here, so
//! `V^A = 0`. The values α = 0.5 and p = 1 are Mobley, Chodera and Dill 2006's — "the modified
//! functional form of Shirts and Pande with a soft core exponent of 1.0 and α = 0.5" (as read
//! through a summary of its PMC text) — which is this form with the λ outside raised to the first
//! power; GROMACS's own default power is 2 and its manual notes that `p = 1` gives smoother
//! derivatives.
//!
//! **What the soft core guarantees, exactly, and the tests check:** at `λ = 1`, `r_sc = r` and
//! the pair is UFF's own, **to the bit** — the code takes [`Pair`]'s own evaluation there,
//! because `x⁶/r⁶` and `(x/r)⁶` round differently; at `λ = 0` the energy and force are `0 × V`,
//! exactly zero; for `λ < 1`, `r_sc⁶ ≥ α σ⁶ (1 − λ)^p > 0`, so the energy and force are finite at
//! `r = 0`.
//!
//! **The derivatives**, with `q = r_sc⁶`, `s = x⁶/q` and `V = D(s² − 2s)`:
//!
//! - `dV/dq = 2D s(1 − s)/q`;
//! - `∂U_sc/∂r = λ dV/dq · 6r⁵`, used as `(∂U/∂r)/r = 6λ dV/dq r⁴` so that `r = 0` is no
//!   division;
//! - `∂U_sc/∂λ = V + λ dV/dq ∂q/∂λ`, `∂q/∂λ = −p α σ⁶ (1 − λ)^(p−1)`. The second term is the chain
//!   rule through `r_sc`, and leaving it out is one of the sabotages the tests catch.
//!
//! # Double decoupling
//!
//! The complex leg goes coupled and unrestrained → **restraint on** (`λ_r` 0 → 1) →
//! **charges off** (`λ_e` 1 → 0) → **van der Waals off** (`λ_v` 1 → 0), and the restraint's
//! release to the standard state is [`Boresch::release_free_energy`]. The solvent leg turns the
//! same two off in water, and
//!
//! `ΔG°_bind = ΔG_solvent − ΔG_complex − ΔG°_release`,
//!
//! with `ΔG_complex` and `ΔG_solvent` each the free energy of decoupling in its leg. Mobley et al.
//! 2007's order — "restrain the ligand harmonically … then annihilate the ligand's partial
//! charges, then decouple its Lennard-Jones interactions" — as read in PMC2104542.
//!
//! # What is not here
//!
//! **No solvent.** A [`ForceField`] with generalized Born is refused ([`AlchemyError::Solvated`]):
//! GB is quadratic in the charges and many-body through the Born radii, so decoupling the ligand
//! there means scaling its charges in the GB sum and deciding what its radii do as it vanishes —
//! step 3c-3's question, not this one's. **No group bonded to its environment**
//! ([`AlchemyError::BondedAcross`]): only non-bonded pairs are scaled.

use crate::boresch::Boresch;
use crate::dynamics::Potential;
use crate::energy::{coulomb, ForceField, Pair};
use std::fmt;

/// α of the soft core, [`SoftCore::default`]: Mobley, Chodera and Dill 2006's 0.5. See the module
/// documentation.
pub const SOFT_CORE_ALPHA: f64 = 0.5;

/// The power `p` of `(1 − λ)` in the soft core, [`SoftCore::default`]: Mobley, Chodera and Dill
/// 2006's 1. See the module documentation.
pub const SOFT_CORE_POWER: f64 = 1.0;

/// An alchemical state: three couplings, each in [0, 1]. See the module documentation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Lambda {
    /// The Boresch restraint's scale: 0 off, 1 on.
    pub restraint: f64,
    /// The group–environment electrostatics' scale: 1 coupled, 0 off.
    pub electrostatics: f64,
    /// The group–environment soft-core van der Waals' λ: 1 coupled, 0 off.
    pub van_der_waals: f64,
}

impl Lambda {
    /// Fully coupled, no restraint: the physical bound state.
    pub const COUPLED: Lambda = Lambda {
        restraint: 0.0,
        electrostatics: 1.0,
        van_der_waals: 1.0,
    };

    /// `(restraint, electrostatics, van_der_waals)`.
    pub const fn new(restraint: f64, electrostatics: f64, van_der_waals: f64) -> Lambda {
        Lambda {
            restraint,
            electrostatics,
            van_der_waals,
        }
    }

    /// The three components, `[restraint, electrostatics, van_der_waals]`.
    pub fn components(&self) -> [f64; 3] {
        [self.restraint, self.electrostatics, self.van_der_waals]
    }

    /// Whether every component is in [0, 1].
    pub fn is_valid(&self) -> bool {
        self.components().iter().all(|l| (0.0..=1.0).contains(l))
    }
}

/// The soft core's parameters: see the module documentation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SoftCore {
    /// α, dimensionless.
    pub alpha: f64,
    /// p, the power of `(1 − λ)`; at least 1, so that `∂U/∂λ` is finite at `λ = 1`.
    pub power: f64,
}

impl Default for SoftCore {
    /// α = [`SOFT_CORE_ALPHA`], p = [`SOFT_CORE_POWER`].
    fn default() -> SoftCore {
        SoftCore {
            alpha: SOFT_CORE_ALPHA,
            power: SOFT_CORE_POWER,
        }
    }
}

/// One soft-core pair's energy and derivatives at one `r` and `λ`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SoftCoreTerm {
    /// `U_sc`, joules per molecule.
    pub energy: f64,
    /// `(∂U_sc/∂r) / r`, J m⁻²: the force on the pair's first atom is `−` this times
    /// `r_first − r_second`. Finite at `r = 0` for `λ < 1`.
    pub de_dr_over_r: f64,
    /// `∂U_sc/∂λ`, joules per molecule.
    pub de_dlambda: f64,
}

impl SoftCore {
    /// The soft-core pair `pair` at separation `r` metres and coupling `lambda`: see the module
    /// documentation for the form and its derivatives. At `lambda == 1` the energy and the radial
    /// derivative are [`Pair`]'s own, bit for bit.
    ///
    /// # Panics
    ///
    /// If `lambda` is outside [0, 1] or the power is below 1.
    pub fn at(&self, pair: &Pair, r: f64, lambda: f64) -> SoftCoreTerm {
        assert!((0.0..=1.0).contains(&lambda), "λ must be in [0, 1]");
        assert!(
            self.power >= 1.0,
            "the soft core's power must be at least 1"
        );
        let x2 = pair.distance * pair.distance;
        let x6 = x2 * x2 * x2;
        let a = self.alpha * 0.5 * x6;
        let one_minus = 1.0 - lambda;
        let dq_dlambda = -self.power * a * one_minus.powf(self.power - 1.0);
        let r2 = r * r;
        if lambda == 1.0 {
            let (energy, de_dr) = pair.at(r);
            let q = r2 * r2 * r2;
            let s = x6 / q;
            let dv_dq = 2.0 * pair.well * s * (1.0 - s) / q;
            return SoftCoreTerm {
                energy,
                de_dr_over_r: de_dr / r,
                de_dlambda: energy + dv_dq * dq_dlambda,
            };
        }
        let q = r2 * r2 * r2 + a * one_minus.powf(self.power);
        let s = x6 / q;
        let v = pair.well * (s * s - 2.0 * s);
        let dv_dq = 2.0 * pair.well * s * (1.0 - s) / q;
        SoftCoreTerm {
            energy: lambda * v,
            de_dr_over_r: 6.0 * lambda * dv_dq * r2 * r2,
            de_dlambda: v + lambda * dv_dq * dq_dlambda,
        }
    }
}

/// One non-bonded pair between the group and its environment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrossPair {
    /// The pair's van der Waals parameters and its two atoms, `atoms[0] < atoms[1]`.
    pub pair: Pair,
    /// The two atoms' partial charges, elementary charges, in the order of `pair.atoms`.
    pub charges: [f64; 2],
}

/// The λ-dependent part of the energy at one state, and its gradient in λ.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Coupling {
    /// `U(x; λ) − U_rest(x)`, joules per molecule: everything λ touches, so that the difference
    /// between two states at one configuration has no `U_rest` in it to round.
    pub energy: f64,
    /// `[∂U/∂λ_r, ∂U/∂λ_e, ∂U/∂λ_v]`, joules per molecule.
    pub gradient: [f64; 3],
}

/// A Hamiltonian with a λ: what [`crate::free_energy`]'s windows sample. [`Decoupling`] is one; a
/// test's harmonic well another.
pub trait Alchemical {
    /// How many atoms it acts on.
    fn len(&self) -> usize;

    /// Whether it acts on no atom.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// `U(x; λ)`, joules per molecule, with the force on every atom written into `forces`
    /// (newtons).
    fn energy_and_forces(&self, at: &[[f64; 3]], lambda: Lambda, forces: &mut [[f64; 3]]) -> f64;

    /// The part of `U(x; λ)` that depends on λ, up to a constant the same at every state, and its
    /// gradient in λ: what TI averages and BAR differences. See [`Coupling`].
    fn coupling(&self, at: &[[f64; 3]], lambda: Lambda) -> Coupling;
}

/// An [`Alchemical`] Hamiltonian held at one state: a [`Potential`] that molecular dynamics can
/// integrate.
#[derive(Clone, Copy, Debug)]
pub struct AtLambda<'a, A: ?Sized> {
    /// The Hamiltonian.
    pub hamiltonian: &'a A,
    /// The state.
    pub lambda: Lambda,
}

impl<A: Alchemical + ?Sized> Potential for AtLambda<'_, A> {
    fn energy_and_forces(&self, at: &[[f64; 3]], forces: &mut [[f64; 3]]) -> f64 {
        self.hamiltonian.energy_and_forces(at, self.lambda, forces)
    }
}

/// Why a [`Decoupling`] could not be built.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AlchemyError {
    /// The force field has a solvent: decoupling under generalized Born is not here yet. See the
    /// module documentation.
    Solvated,
    /// A bond joins the group to its environment: these two atoms.
    BondedAcross {
        /// The bond's atoms.
        atoms: [usize; 2],
    },
    /// The group has no atom, or every atom.
    NothingToDecouple,
}

impl fmt::Display for AlchemyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AlchemyError::Solvated => f.write_str(
                "the force field has a solvent, and decoupling under generalized Born is not \
                 implemented",
            ),
            AlchemyError::BondedAcross { atoms } => write!(
                f,
                "atoms {} and {} are bonded across the group's boundary, and only non-bonded \
                 pairs are decoupled",
                atoms[0], atoms[1]
            ),
            AlchemyError::NothingToDecouple => {
                f.write_str("the group must have at least one atom and leave at least one out")
            }
        }
    }
}

impl std::error::Error for AlchemyError {}

/// A group decoupled from its environment by three λs: see the module documentation.
#[derive(Clone, Debug, PartialEq)]
pub struct Decoupling {
    atoms: usize,
    group: Vec<bool>,
    rest: Option<ForceField>,
    cross: Vec<CrossPair>,
    restraint: Option<Boresch>,
    soft_core: SoftCore,
}

impl Decoupling {
    /// `force_field` with the atoms `group` marks decoupled from the rest: its pairs across the
    /// partition taken out of it and scaled, every other term kept. No restraint, and
    /// [`SoftCore::default`].
    ///
    /// # Errors
    ///
    /// [`AlchemyError::Solvated`] for a force field with generalized Born,
    /// [`AlchemyError::BondedAcross`] when a bond crosses the partition, and
    /// [`AlchemyError::NothingToDecouple`] for an empty or whole group.
    ///
    /// # Panics
    ///
    /// If `group` is not one per atom.
    pub fn new(force_field: &ForceField, group: &[bool]) -> Result<Decoupling, AlchemyError> {
        if force_field.solvation().is_some() {
            return Err(AlchemyError::Solvated);
        }
        if !group.iter().any(|&g| g) || group.iter().all(|&g| g) {
            return Err(AlchemyError::NothingToDecouple);
        }
        let (rest, across) = force_field
            .split_across(group)
            .map_err(|atoms| AlchemyError::BondedAcross { atoms })?;
        let q = force_field.charges();
        let cross = across
            .into_iter()
            .map(|pair| CrossPair {
                pair,
                charges: [q[pair.atoms[0]], q[pair.atoms[1]]],
            })
            .collect();
        Ok(Decoupling {
            atoms: group.len(),
            group: group.to_vec(),
            rest: Some(rest),
            cross,
            restraint: None,
            soft_core: SoftCore::default(),
        })
    }

    /// A decoupling of `atoms` atoms that is nothing but the pairs `cross` between `group` and the
    /// rest: no other term. For a model system, or a test whose free energy is known.
    ///
    /// # Panics
    ///
    /// If `group` is not `atoms` long, or a pair does not join a group atom to one outside it.
    pub fn from_pairs(atoms: usize, group: Vec<bool>, cross: Vec<CrossPair>) -> Decoupling {
        assert_eq!(group.len(), atoms, "one mark per atom");
        for c in &cross {
            let [i, j] = c.pair.atoms;
            assert!(
                i < atoms && j < atoms && group[i] != group[j],
                "pair {i}–{j} does not cross the group's boundary"
            );
        }
        Decoupling {
            atoms,
            group,
            rest: None,
            cross,
            restraint: None,
            soft_core: SoftCore::default(),
        }
    }

    /// The same decoupling with a Boresch restraint, scaled by `λ_r`.
    ///
    /// # Panics
    ///
    /// If one of its atoms is past the atom count.
    pub fn with_restraint(mut self, restraint: Boresch) -> Decoupling {
        assert!(
            restraint
                .receptor
                .iter()
                .chain(&restraint.ligand)
                .all(|&a| a < self.atoms),
            "a restraint atom past the atom count"
        );
        self.restraint = Some(restraint);
        self
    }

    /// The same decoupling with another soft core.
    pub fn with_soft_core(mut self, soft_core: SoftCore) -> Decoupling {
        self.soft_core = soft_core;
        self
    }

    /// Which atoms are decoupled.
    pub fn group(&self) -> &[bool] {
        &self.group
    }

    /// The pairs that are scaled, in the order the force field held them.
    pub fn cross_pairs(&self) -> &[CrossPair] {
        &self.cross
    }

    /// Everything that is not scaled: the force field less the cross pairs; `None` for
    /// [`Decoupling::from_pairs`].
    pub fn rest(&self) -> Option<&ForceField> {
        self.rest.as_ref()
    }

    /// The restraint, if any.
    pub fn restraint(&self) -> Option<&Boresch> {
        self.restraint.as_ref()
    }

    /// The soft core.
    pub fn soft_core(&self) -> SoftCore {
        self.soft_core
    }

    /// The λ-dependent part at `at` and `lambda`, adding its forces to `forces` when given.
    fn coupling_into(
        &self,
        at: &[[f64; 3]],
        lambda: Lambda,
        mut forces: Option<&mut [[f64; 3]]>,
    ) -> Coupling {
        assert_eq!(at.len(), self.atoms, "one position per atom");
        assert!(lambda.is_valid(), "every λ must be in [0, 1]: {lambda:?}");
        let mut energy = 0.0;
        let mut gradient = [0.0; 3];
        for c in &self.cross {
            let [i, j] = c.pair.atoms;
            let d = [
                at[i][0] - at[j][0],
                at[i][1] - at[j][1],
                at[i][2] - at[j][2],
            ];
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            let sc = self.soft_core.at(&c.pair, r, lambda.van_der_waals);
            energy += sc.energy;
            gradient[2] += sc.de_dlambda;
            let mut radial = sc.de_dr_over_r;
            let [qi, qj] = c.charges;
            if qi != 0.0 && qj != 0.0 {
                let e = coulomb(qi, qj, r);
                gradient[1] += e;
                // Skipped at λ_e = 0, so that a pair at r = 0 is 0, not 0 × ∞.
                if lambda.electrostatics != 0.0 {
                    energy += lambda.electrostatics * e;
                    radial -= lambda.electrostatics * e / (r * r);
                }
            }
            if let Some(f) = forces.as_deref_mut() {
                for k in 0..3 {
                    f[i][k] -= radial * d[k];
                    f[j][k] += radial * d[k];
                }
            }
        }
        if let Some(b) = &self.restraint {
            let u = match forces {
                Some(f) => b.add_forces(at, f, lambda.restraint),
                None => b.energy(at),
            };
            gradient[0] = u;
            energy += lambda.restraint * u;
        }
        Coupling { energy, gradient }
    }
}

impl Alchemical for Decoupling {
    fn len(&self) -> usize {
        self.atoms
    }

    fn energy_and_forces(&self, at: &[[f64; 3]], lambda: Lambda, forces: &mut [[f64; 3]]) -> f64 {
        assert_eq!(forces.len(), self.atoms, "one force per atom");
        let rest = match &self.rest {
            Some(ff) => {
                let ev = ff.evaluate(at);
                forces.copy_from_slice(&ev.forces);
                ev.energy.total
            }
            None => {
                forces.iter_mut().for_each(|f| *f = [0.0; 3]);
                0.0
            }
        };
        rest + self.coupling_into(at, lambda, Some(forces)).energy
    }

    fn coupling(&self, at: &[[f64; 3]], lambda: Lambda) -> Coupling {
        self.coupling_into(at, lambda, None)
    }
}
