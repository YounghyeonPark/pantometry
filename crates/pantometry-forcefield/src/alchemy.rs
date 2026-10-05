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
//! `U(x; λ) = U_rest(x) + λ_r U_B(x) + Σ_cross [λ_e C_ij(r) + U_sc(r; λ_v)] + G(x; λ_e, λ_v)`,
//!
//! the last term the generalized Born solvation energy when the force field has one (below).
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
//!   vacuum the choice matters less: an annihilated ligand's
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
//! # Generalized Born
//!
//! A force field with OBC II ([`crate::solvation`]) has its solvation term taken out of `U_rest`
//! and evaluated as `G(x; λ_e, λ_v)` ([`GeneralizedBorn::decoupled`]):
//!
//! - **The group's charges are `λ_e q`, everywhere they enter eq 2**: its Born terms and its
//!   pairs among themselves go as `λ_e²`, its pairs with the environment as `λ_e`. A charge
//!   scaled is a charge scaled, so at fixed positions and `λ_v` the whole electrostatic part is
//!   `A + λ_e (B + Σ C_ij) + λ_e² D` exactly, and `∂U/∂λ_e = B + Σ C_ij + 2 λ_e D` is analytic. The
//!   group's intramolecular Coulomb term stays in `U_rest` at full charge: the decoupled ligand
//!   keeps its vacuum self-interaction, and loses only its solvation.
//! - **The descreening integral of eq 5 between atoms on opposite sides is scaled by `λ_v`**,
//!   that between atoms on the same side is not: `I_i = I_i^same + λ_v I_i^cross`. So the group's
//!   volume leaves the environment's Born radii together with its van der Waals, and the
//!   environment's leaves the group's. `∂G/∂λ_v = Σ_i ∂G/∂R_i · dR_i/dI_i · I_i^cross`, and the
//!   chain rule through the radii takes a cross pair's `dI/dr` times `λ_v`.
//! - **No nonpolar surface term.** [`crate::solvation::surface_area`] has no useful gradient,
//!   and is not in the force field ([`crate::solvation`], "The nonpolar part"), so it is in
//!   neither leg of a cycle; a caller can add `0.005 kcal mol⁻¹ Å⁻² × ΔA` at the end states as a
//!   separate estimate, and say so.
//!
//! **So the fully decoupled state, `λ_e = λ_v = 0`, is the environment alone in solvent plus the
//! group alone in vacuum**: the environment's radii are its own (the cross integrals are
//! multiplied by zero), and every GB term with a group charge in it is zero, so the group neither
//! screens nor is screened, nor is solvated. The tests hold that exactly. Both legs of a double
//! decoupling therefore end in the same state of the ligand — alone, in vacuum, with its
//! intramolecular terms — and the solvent leg, the group being the whole molecule, is its
//! electrostatic solvation free energy under the model.
//!
//! **What practice does, read.** openmmtools' `AbsoluteAlchemicalFactory`
//! (`openmmtools/alchemy/alchemy.py`, `_alchemically_modify_GBSAOBCForce`, read at commit
//! `f6ef22a8` of `github.com/choderalab/openmmtools`) builds the same end state: an alchemical
//! atom's charge is multiplied by `lambda_electrostatics` in the pair term, and its descreening of
//! any other atom (the computed value `I`) by the same factor, so a decoupled ligand neither
//! carries charge nor descreens. It differs on the path, in two ways chosen differently here: it
//! scales the Born self term `q²/B` linearly in λ rather than as the square of a scaled charge,
//! and it removes the descreening with the charges rather than with the van der Waals. Neither
//! changes ΔG, which depends on the end states only. Linking the descreening to `λ_v` keeps the
//! radii fixed while the charges are turned off, so that segment is exactly quadratic in `λ_e` —
//! which is what makes a Born ion's TI exact on a trapezoid — and takes the ligand's dielectric
//! cavity away with the volume that makes it.
//!
//! # What it measures: benzene and T4 lysozyme L99A under OBC II
//!
//! `tests/benzene_bound_in_generalized_born.rs`, ignored, release, one core. With the convention
//! above, `ΔG°_bind = ΔG_solvent − ΔG_complex − ΔG°_release`, each `ΔG` the free energy of
//! decoupling in its leg and `ΔG°_release` [`Boresch::release_free_energy`], negative:
//!
//! | | kcal/mol |
//! | --- | --- |
//! | solvent leg: benzene's charges off in OBC II (4 seeds × 11 windows × 100 ps) | +2.3395 ± 0.0003 |
//! | complex leg: restraint on | +0.718 ± 0.065 |
//! | complex leg: charges off | −5.172 ± 0.011 |
//! | complex leg: van der Waals off, the cavity's descreening with it | +8.501 ± 0.301 |
//! | complex leg (15 windows × 12 ps; TI +4.073 ± 0.313) | +4.047 ± 0.309 |
//! | restraint released to 1 M, analytic | −7.784 |
//! | **ΔG°_bind** | **+6.08 ± 0.31** |
//! | nonpolar, `0.005 kcal mol⁻¹ Å⁻² × ΔA`, ΔA = −293.3 Å² at the start | −1.47 |
//! | experiment (Mobley et al. 2007, Table 1) | −5.19 ± 0.16 |
//!
//! **The model does not bind benzene.** The electrostatics cost +7.51 (the solvent leg's +2.34
//! less the complex's −5.17), where in vacuum (3c-1) the complex's charge segment was +1.48, a
//! contribution of −1.48 to binding: OBC II with QEq charges
//! desolvates the ligand and the pocket by more than the cross terms give back, as 2c-2 and A-2
//! measured at fixed pose. The van der Waals segment is +8.50 against vacuum's +13.80. Removing
//! benzene's volume lets GB solvate the cavity's walls as water, which L99A's real cavity is not
//! (A-2: it is empty); A-2 measured that part of the polar term, the pocket desolvated by
//! benzene's volume, at +4.31 kcal/mol in the 10 Å binding at the crystal pose, most of the 5.3
//! difference. The rest — another pose, minimised under GB, and sampling — is not separated here. The solvent leg's hydration free energy, −2.34 polar and
//! −1.13 with the nonpolar term, against FreeSolv's −0.90 ± 0.20, is 0.23 too negative, about
//! the experiment's own uncertainty: for benzene the model's over-solvation of small molecules is
//! small. Nothing is asserted against experiment.
//!
//! # What is not here
//!
//! **No group bonded to its environment** ([`AlchemyError::BondedAcross`]): only non-bonded pairs
//! are scaled. **No nonpolar solvation term**, above.

use crate::boresch::Boresch;
use crate::dynamics::Potential;
use crate::energy::{coulomb, ForceField, Pair};
use crate::solvation::{DecoupledSolvation, GeneralizedBorn};
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

    /// [`Alchemical::coupling`] at each of `states`, at one configuration: **each entry that
    /// state's own [`Alchemical::coupling`], to the bit**. The default calls it once per state; an
    /// implementation may share what the states have in common, as [`Decoupling`] shares the
    /// generalized Born descreening.
    fn couplings(&self, at: &[[f64; 3]], states: &[Lambda]) -> Vec<Coupling> {
        states.iter().map(|&l| self.coupling(at, l)).collect()
    }
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
    /// A bond joins the group to its environment: these two atoms.
    BondedAcross {
        /// The bond's atoms.
        atoms: [usize; 2],
    },
    /// The group has no atom, or — in a force field without a solvent, where it would have
    /// nothing to be decoupled from — every atom.
    NothingToDecouple,
}

impl fmt::Display for AlchemyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AlchemyError::BondedAcross { atoms } => write!(
                f,
                "atoms {} and {} are bonded across the group's boundary, and only non-bonded \
                 pairs are decoupled",
                atoms[0], atoms[1]
            ),
            AlchemyError::NothingToDecouple => f.write_str(
                "the group must have at least one atom, and leave at least one out unless \
                     there is a solvent to decouple it from",
            ),
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
    solvent: Option<Solvent>,
}

/// The generalized Born term a [`Decoupling`] scales, and the charges it is evaluated with.
#[derive(Clone, Debug, PartialEq)]
struct Solvent {
    model: GeneralizedBorn,
    charges: Vec<f64>,
}

impl Decoupling {
    /// `force_field` with the atoms `group` marks decoupled from the rest: its pairs across the
    /// partition taken out of it and scaled, every other term kept. No restraint, and
    /// [`SoftCore::default`]. **A force field with generalized Born** has its solvation term taken
    /// out of the rest and scaled as the module documentation's "Generalized Born" says; the group
    /// may then be every atom, which decouples a molecule from the solvent alone (a solvent leg).
    ///
    /// # Errors
    ///
    /// [`AlchemyError::BondedAcross`] when a bond crosses the partition, and
    /// [`AlchemyError::NothingToDecouple`] for an empty group, or a whole one without a solvent.
    ///
    /// # Panics
    ///
    /// If `group` is not one per atom.
    pub fn new(force_field: &ForceField, group: &[bool]) -> Result<Decoupling, AlchemyError> {
        let solvated = force_field.solvation().is_some();
        if !group.iter().any(|&g| g) || (!solvated && group.iter().all(|&g| g)) {
            return Err(AlchemyError::NothingToDecouple);
        }
        let (rest, across) = force_field
            .split_across(group)
            .map_err(|atoms| AlchemyError::BondedAcross { atoms })?;
        let q = force_field.charges();
        let solvent = force_field.solvation().map(|model| Solvent {
            model: model.clone(),
            charges: q.to_vec(),
        });
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
            rest: Some(rest.without_solvation()),
            cross,
            restraint: None,
            soft_core: SoftCore::default(),
            solvent,
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
            solvent: None,
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

    /// The same decoupling with the generalized Born model `model`, evaluated with `charges`
    /// (elementary charges) and scaled as the module documentation's "Generalized Born" says, in
    /// place of any it had: for a model system of [`Decoupling::from_pairs`] — a Born ion, say,
    /// whose solvation free energy is known.
    ///
    /// # Panics
    ///
    /// If `model` or `charges` is not one per atom.
    pub fn with_solvation(mut self, model: GeneralizedBorn, charges: Vec<f64>) -> Decoupling {
        assert_eq!(model.len(), self.atoms, "one GB atom per atom");
        assert_eq!(charges.len(), self.atoms, "one charge per atom");
        self.solvent = Some(Solvent { model, charges });
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

    /// Everything that is not scaled: the force field less the cross pairs and less its solvation
    /// term, which is scaled; `None` for [`Decoupling::from_pairs`].
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

    /// The generalized Born model that is scaled, if the force field had one.
    pub fn solvation(&self) -> Option<&GeneralizedBorn> {
        self.solvent.as_ref().map(|s| &s.model)
    }

    /// The solvation term alone at each of `states`, from one descreening; with `forces`, there
    /// must be one state, and its force is added.
    fn solvation_at(
        &self,
        at: &[[f64; 3]],
        states: &[Lambda],
        forces: Option<&mut [[f64; 3]]>,
    ) -> Option<Vec<DecoupledSolvation>> {
        let s = self.solvent.as_ref()?;
        Some(match forces {
            Some(f) => {
                assert_eq!(states.len(), 1, "forces at one state");
                let l = states[0];
                vec![s.model.decoupled(
                    &s.charges,
                    &self.group,
                    at,
                    l.electrostatics,
                    l.van_der_waals,
                    f,
                )]
            }
            None => {
                let pairs: Vec<[f64; 2]> = states
                    .iter()
                    .map(|l| [l.electrostatics, l.van_der_waals])
                    .collect();
                s.model
                    .decoupled_energies(&s.charges, &self.group, at, &pairs)
            }
        })
    }

    /// The λ-dependent part at `at` and `lambda`, adding its forces to `forces` when given: the
    /// pairs and the restraint, then the solvation term.
    fn coupling_into(
        &self,
        at: &[[f64; 3]],
        lambda: Lambda,
        mut forces: Option<&mut [[f64; 3]]>,
    ) -> Coupling {
        let mut c = self.pairs_into(at, lambda, forces.as_deref_mut());
        if let Some(g) = self.solvation_at(at, &[lambda], forces) {
            add_solvation(&mut c, g[0]);
        }
        c
    }

    /// The cross pairs and the restraint at `at` and `lambda`, adding their forces to `forces`
    /// when given.
    fn pairs_into(
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

    /// Each state's pairs and restraint as [`Alchemical::coupling`] computes them, and the
    /// solvation term from one shared descreening ([`GeneralizedBorn::decoupled_energies`]):
    /// each entry is that state's [`Alchemical::coupling`] to the bit.
    fn couplings(&self, at: &[[f64; 3]], states: &[Lambda]) -> Vec<Coupling> {
        let solvation = self.solvation_at(at, states, None);
        states
            .iter()
            .enumerate()
            .map(|(k, &l)| {
                let mut c = self.pairs_into(at, l, None);
                if let Some(g) = &solvation {
                    add_solvation(&mut c, g[k]);
                }
                c
            })
            .collect()
    }
}

/// Adds a solvation term's energy and derivatives to a coupling.
fn add_solvation(c: &mut Coupling, g: DecoupledSolvation) {
    c.energy += g.energy;
    c.gradient[1] += g.d_electrostatics;
    c.gradient[2] += g.d_van_der_waals;
}
