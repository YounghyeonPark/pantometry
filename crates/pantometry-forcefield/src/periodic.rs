//! Periodic boundaries: an orthorhombic box, UFF's van der Waals cut off with its long-range
//! correction, Ewald electrostatics ([`crate::ewald`]), and the alchemical decoupling of a group
//! inside the box.
//!
//! This is step W1 of the explicit-water track, and W2 adds the water: [`PeriodicForceField::tip3p`]
//! is rigid TIP3P alone and [`PeriodicForceField::solvated`] a UFF solute in it, the model and its
//! constraints in [`crate::water`]. W3 is benzene's hydration and W4 the solvated complex.
//!
//! # The box
//!
//! **Orthorhombic**: three edges at right angles, of any lengths ([`PeriodicBox`]). It is enough
//! for a water box and a solvated protein, and it keeps the minimum image a rounding per axis,
//! `d − L round(d/L)`. A triclinic box — a truncated octahedron saves about 30% of the water
//! around a globular protein — is not here.
//!
//! **Minimum image** for every non-bonded pair: the real-space cutoff is at most half the shortest
//! edge, so a pair inside it has exactly one image there. **Wrapping** ([`PeriodicBox::wrap`])
//! puts a point in `[0, L)` per axis; the evaluation wraps for its own use, and nothing here
//! writes a wrapped position back. Molecular dynamics never wraps either, so a trajectory's
//! positions stay continuous — what a mean-square displacement or a restraint between molecules
//! wants — and the force field accepts positions in any image.
//!
//! # Molecules stay whole: unwrapped through their bonds
//!
//! **Each evaluation makes every molecule whole** ([`PeriodicForceField::whole`]) before the bonded
//! terms see it: walking each molecule's bond graph breadth-first from its lowest-numbered atom,
//! each atom is put at the image of itself nearest the atom it was reached from. The bonded
//! terms then run on those positions unchanged.
//!
//! The other choice was the minimum image inside every bonded term. Against it: the bend,
//! torsion and inversion of [`crate::angular`] read positions, not displacements, and would each
//! have to be rewritten; and the atomic virial `Σ r_i ⊗ F_i` of a bonded term is defined only on
//! whole molecules anyway. For the walk: it is `O(N)`, done once, leaves every bonded term as it
//! is checked, and **changes no bit** when no bond crosses a face, because an atom whose image is
//! already the nearest is taken as it is. It requires every bond to be shorter than half the
//! shortest edge, which a box twice a bond's length satisfies; a molecule as a whole may span the
//! box. The tests hold the energy to every atom wrapped separately, and to every atom moved by a
//! random whole number of box lengths.
//!
//! # Van der Waals: a plain cutoff and the long-range correction
//!
//! UFF's pair (eq 20, [`crate::energy`]) is `D [(x/r)¹² − 2 (x/r)⁶]`, which is exactly
//! Lennard-Jones's `4ε [(σ/r)¹² − (σ/r)⁶]` with `ε = D` and `σ = x / 2^⅙` — `4εσ¹² = D x¹²`,
//! `4εσ⁶ = 2 D x⁶` — so every Lennard-Jones result applies to it as it stands. **It is cut off at
//! the Ewald real-space cutoff, unshifted, and the rest is added as the analytic long-range
//! correction** of Allen and Tildesley (*Computer Simulation of Liquids*, 2nd ed., §2.8, eqs 2.144
//! and 2.145, which assume `g(r) = 1` beyond `r_c`):
//!
//! ```text
//! E_tail = (2π/V) Σ_i Σ_j ∫_{r_c}^∞ r² u_ij(r) dr = (2π/V) Σ_i Σ_j D_ij [x_ij¹²/(9 r_c⁹) − 2 x_ij⁶/(3 r_c³)]
//! P_tail = −(2π/3V²) Σ_i Σ_j ∫_{r_c}^∞ r³ u_ij′(r) dr = (2π/3V²) Σ_i Σ_j D_ij [4 x_ij¹²/(3 r_c⁹) − 4 x_ij⁶/r_c³]
//! ```
//!
//! over every ordered pair, `i = j` included, as the uniform-fluid correction counts them.
//! **UFF's geometric combination makes both sums separable**: `D_ij x_ij¹² = (√D_i x_i⁶)(√D_j x_j⁶)`
//! and `D_ij x_ij⁶ = (√D_i x_i³)(√D_j x_j³)`, so the double sums are squares of single ones, `O(N)`.
//! For one type, both are A&T's printed `(8/3)πNρεσ³[(σ/r_c)⁹/3 − (σ/r_c)³]` and
//! `(16/3)πρ²εσ³[(2/3)(σ/r_c)⁹ − (σ/r_c)³]`, which the tests check, with the integrals by
//! quadrature. Under Lorentz–Berthelot they do not factor, and are summed by atom class instead:
//! see the next section.
//!
//! Why plain, not shifted or switched: the correction is exact for a truncated potential with
//! `g = 1` beyond the cutoff, and needs no correction of its own as a shifted or switched one
//! does; and Ewald's real-space sum is cut the same way. **What it costs**: the energy jumps by
//! `u(r_c)` when a pair crosses the cutoff — for two `O_3` oxygens at 9 Å, `−4.1e-4` kcal/mol.
//!
//! **The pressure, for W2.** The virial ([`PeriodicEvaluation::virial`]) is the exact strain
//! derivative of this model's energy, `−∂U/∂ε`, and the correction's part of it is `E_tail` on
//! the diagonal (it goes as `1/V`). A&T's `P_tail` differs from that `E_tail/V` by
//! `(2π/3V²) Σ_i Σ_j r_c³ u_ij(r_c)` — the impulsive term of a potential that steps at `r_c`,
//! which no configuration's derivative contains and an average over the ensemble does. It has
//! the sign of `u(r_c)`: negative where the pair is attractive at the cutoff, as it is for every
//! UFF pair at any cutoff past `x_IJ`, so `P_tail` is the more negative of the two. So the
//! pressure of the untruncated fluid is `(Σ m v² + W)/3V − E_tail/V + P_tail`, with
//! [`PeriodicForceField::pressure_correction`] giving `P_tail`.
//!
//! # The force field
//!
//! [`PeriodicForceField`] is UFF in the box: the four bonded terms on whole molecules, every
//! non-bonded pair that is neither 1-2 nor 1-3 inside the cutoff by minimum image (a cell list
//! of cells at least half the cutoff wide, each pair visited once — or, on request, a Verlet
//! list with a skin that visits the same pairs in the same order, [`crate::neighbours`]), the
//! correction above, and [`Ewald`] with the same exclusions. It is a [`Potential`], so [`crate::MolecularDynamics`]
//! integrates it unchanged. It never builds the `N²` pair list [`ForceField::new`] keeps. The
//! atoms are the solute's, then the waters', then any ions' ([`PeriodicForceField::with_ions`]):
//! point charges with a van der Waals term, bonded to nothing, for a solvated complex
//! ([`crate::solvated`]). A protein is a solute like any other: its chain is made whole through its
//! peptide bonds as a water is through its two, so no bonded term needed a periodic version, and
//! `tests/benzene_bound_in_tip3p.rs` holds a pocket of 181L cut by every face to the vacuum force
//! field.
//!
//! # Another force field's non-bonded terms: overrides, Lorentz–Berthelot, scaled 1-4 pairs
//!
//! Three opt-in rules, so that a solute can carry another force field's van der Waals — GAFF's,
//! say — on UFF's bonded terms. **At their defaults every bit is as it was**: nothing is computed
//! differently, and the tests hold an evaluation with every rule set to its default to the plain
//! field's bits, and the plain field's default path to the bits HEAD `b6558df` computed before
//! these rules existed, on two systems whose every input is the same on every platform.
//!
//! - **[`PeriodicForceField::with_nonbonded`]** sets the `[x, D]` of any run of atoms, bonded
//!   molecules included, in UFF's form: `x` the minimum of the atom's pair with itself and `D` its
//!   depth, so AMBER's `R*/2` and `ε` are `x = 2 (R*/2)` and `D = ε`. Refused by name
//!   ([`NonbondedError`]) unless the lengths match, the run is inside the atoms, and every value is
//!   finite with `D ≥ 0` and `x > 0` — or `[0, 0]`, an atom with no van der Waals at all, as
//!   GAFF's `ho` is written and as this field holds TIP3P's hydrogens. `D = 0` makes every pair of
//!   the atom zero under either rule, whatever its `x`; `x = 0` with a well is refused, since it
//!   would make every pair of the atom zero and drop the well without a word.
//! - **[`PeriodicForceField::with_mixing`]**: [`Mixing::Geometric`], UFF's eqs 21b and 22,
//!   `x_ij = √(x_i x_j)`, the default; or [`Mixing::LorentzBerthelot`], `x_ij = (x_i + x_j)/2` and
//!   `D_ij = √(D_i D_j)`. `x` is `2^⅙ σ`, so the arithmetic mean of `x` is that of `σ`, and AMBER's
//!   `R*_ij = R*_i + R*_j`. Every pair goes through one function, so the cut-off pairs, the soft
//!   core of a decoupling, its own pairs in vacuum and its cross pairs all follow the rule.
//! - **[`PeriodicForceField::with_one_four`]`(s_v, s_c)`** scales the van der Waals and the Coulomb
//!   energy of every 1-4 pair — two atoms whose shortest path through the bonds is three bonds —
//!   by `s_v` and `s_c`: GAFF's and AMBER's are `½` and `1/1.2`. **Exclusion comes first**: a pair
//!   that is 1-2 or 1-3 by any path is excluded whatever other paths it has (so a four- or
//!   five-membered ring has no 1-4 pair across it), and a pair that is 1-4 by two paths (the para
//!   atoms of a six-membered ring) is scaled once. Like an excluded pair, a 1-4 pair is evaluated
//!   wherever it is, not only inside the cutoff, at its minimum image; the pair search leaves it
//!   out. The default `(1, 1)` is UFF's, every 1-4 pair a pair like any other.
//!
//! **A 1-4 pair's Coulomb term in Ewald is a partial exclusion.** The reciprocal sum and the self
//! term together contain `k_e q_i q_j erf(αr)/r` for every pair, whether or not it is excluded,
//! since `|S(k)|²` cannot tell pairs apart. Taking the pair out of the real-space loop and adding
//!
//! ```text
//! s_c · k_e q_i q_j erfc(αr)/r  −  (1 − s_c) · k_e q_i q_j erf(αr)/r
//! ```
//!
//! leaves it, with the reciprocal sum's `erf(αr)/r`, at `s_c (erfc + erf)/r = s_c/r` exactly. The
//! first part is booked as real space ([`EwaldEnergy::real`]) and the second as an exclusion
//! ([`EwaldEnergy::excluded`]); `s_c = 1` is the real-space pair and `s_c = 0` the excluded pair,
//! each uncut. The force is the derivative of the two along the pair, and the pair's virial is
//! `d ⊗ F` like any central pair's, `d` its minimum-image separation; the self term and the
//! background do not see it. The tests hold one molecule in growing boxes to its vacuum energy
//! with the 1-4 pairs scaled, less the tinfoil dipole term, the rest falling as `L⁻⁵`.
//!
//! **The long-range correction under Lorentz–Berthelot** does not factor: `D_ij x_ij¹²` is not a
//! product of one function of `i` and one of `j`. It is summed by class instead — atoms with the
//! same `[x, D]`, bit for bit, are one class `t` of `n_t` atoms — as
//! `Σ_i Σ_j D_ij x_ij^m = Σ_t Σ_u n_t n_u D_tu x_tu^m` for `m` = 12 and 6, which is exact, the
//! classes found in `O(N T)` and summed in `O(T²)` for `T` classes, not `O(N²)`. A solvent of
//! one molecule has a few classes, and a protein of ff14SB's types thirty. The decoupling's cross
//! part is the same sum over the group's classes and the rest's.
//!
//! **In a decoupling** ([`PeriodicDecoupling`]) the group's own 1-4 pairs are intramolecular and
//! stay in vacuum at every λ, scaled, with its other own pairs; no 1-4 pair crosses the partition,
//! since no bond does. So the decoupled state is still the rest in its box and the group alone in
//! vacuum under this field's non-bonded rules — which with any rule set is no longer the vacuum
//! [`ForceField`]'s, that being UFF's.
//!
//! # Decoupling inside the box
//!
//! [`PeriodicDecoupling`] is [`crate::Decoupling`] under periodic boundaries, an [`Alchemical`] for
//! [`crate::Windows`]: the group's van der Waals pairs with the rest through the same soft core,
//! its Coulomb pairs with the rest scaled linearly by `λ_e` — **in Ewald, as the cross terms of
//! the sum**: `E(all) − E(rest) − E(group)` is bilinear in the two sets of charges, and is its real
//! pairs across, `2 Re(S_rest* S_group)` in reciprocal space, and the cross part of the
//! background, the self terms cancelling exactly. So `∂U/∂λ_e` is that cross energy, at every λ.
//!
//! **The group's own non-bonded terms are not periodic**: its pairs among themselves, in vacuum,
//! all of them, on its whole molecules — the vacuum [`ForceField`]'s for that molecule, at the
//! default non-bonded rules — and no
//! interaction of the group with its own images, at any λ. A choice, and the reason: the decoupled
//! state is then the rest in its box plus the group alone in vacuum, **exactly the decoupled
//! state of [`crate::Decoupling`]**, so both legs of a cycle end in the same state of the ligand
//! whatever box each leg used. Kept periodic, the decoupled ligand would interact with its own
//! images, `−2π μ²/(3V)` for a dipole `μ` under tinfoil, different in each leg's box. The cost is
//! that the coupled end state leaves out that interaction too — which is a finite-size artefact,
//! and small: for a neutral ligand of dipole `μ` in a box `L` across, `2π μ²/(3L³)`. The correction
//! for the long-range van der Waals of the group with the rest is scaled linearly by `λ_v`; for
//! `r ≥ r_c` the soft core's `r_sc⁶ = r⁶ + α σ⁶ (1 − λ)` differs from `r⁶` by at most `α σ⁶/r_c⁶`,
//! `9e-4` for two `O_3` at 9 Å, so linear is the soft core's own tail to that.
//!
//! **For a neutral solute in water, W3, the choice costs nothing measurable.** Benzene's charges
//! sum to zero, so there is no Wigner term `ξ q²/2L` to correct, and its dipole would be zero by
//! symmetry — QEq's charges at 181L's crystal pose are not quite symmetric, so it is small rather
//! than zero, and with it the tinfoil term `2π μ²/3V`; what the coupled state leaves out is that,
//! benzene's quadrupole and its van der Waals with its own images, a box away.
//! `tests/benzene_hydrated_in_tip3p.rs` measures it at the start of the hydration run:
//! −0.008 kcal/mol, with a dipole of 0.021 e Å whose tinfoil term is 2 × 10⁻⁵. The coupled state is then benzene in water, and the decoupled
//! state the water in its box and benzene in vacuum — **the hydration free energy is minus the
//! free energy of decoupling, with no vacuum leg to run**, because benzene's intramolecular terms
//! are the same in both and in vacuum.
//!
//! **The coupling is evaluated from the cross terms alone** ([`Alchemical::couplings`] on a
//! [`PeriodicDecoupling`]): the group's pairs with the rest inside the cutoff — their distances by
//! the evaluation's own wrapped minimum image, so each is its bits — the reciprocal sum's cross
//! terms, the background's and the correction's, then each state's soft core over those pairs. No
//! pair of the rest is computed, so a sample at 29 states costs a sixth of one evaluation of the
//! 1 533-atom box. Its sums are the evaluation's λ-dependent terms in another order, which a unit
//! test holds to 64 ε of the evaluation's parts.
//!
//! # The cost, measured
//!
//! `tests/a_molecule_in_a_periodic_box.rs`, `the_cost`: release, one core,
//! `x86_64-pc-windows-gnu`, waters at 0.03343 Å⁻³ with their exclusions, `r_c` = 9 Å (7.8 Å in the
//! smallest box), one evaluation of the Ewald sum, milliseconds:
//!
//! | atoms | box (Å) | δ = 10⁻⁵: wave vectors, total (real, reciprocal) | δ = 10⁻⁶ |
//! | --- | --- | --- | --- |
//! | 375 | 15.5 | 257, 2.0 (1.5, 0.4) | 522, 2.2 (1.5, 0.8) |
//! | 1 029 | 21.7 | 423, 8.5 (6.8, 1.7) | 895, 10.1 (6.7, 3.4) |
//! | 3 000 | 31.0 | 1 051, 37.7 (25.8, 11.9) | 2 192, 51.1 (25.4, 25.7) |
//! | 6 591 | 40.4 | 1 955, 104 (56, 48) | 4 300, 160 (55, 105) |
//! | 12 288 | 49.7 | 3 309, 249 (97, 152) | 7 385, 434 (95, 338) |
//! | 24 000 | 62.1 | 5 533, 697 (198, 499) | 12 838, 1 347 (198, 1 149) |
//!
//! The whole force field on the 3 000-atom box is 40.2 ms at `10⁻⁵`. Real space is linear in `N`;
//! the reciprocal sum, at a fixed cutoff, goes as `N²` (here `N^1.8`). **So PME is needed at W4's
//! scale**, 20 000–30 000 atoms, where the reciprocal sum is most of a second an evaluation; at W3's
//! few thousand, real space is half the cost, and a neighbour list would save more than PME.
//!
//! **With the particle mesh** ([`PeriodicForceField::with_mesh`], [`crate::pme`]) the reciprocal
//! part of the 24 000-atom box is 33 ms at `10⁻⁵` and 23 ms at `10⁻⁶` with the order and grid
//! [`crate::PmeParameters::for_accuracy`] chooses, and 4.9 and 5.6 ms at 3 000 atoms, timed on a
//! machine under load, against 681 and 1 575 ms for the classical sum at 24 000; real space
//! is then nearly all of it. The decoupling's cross terms go through the mesh as well.
//!
//! # Not here
//!
//! **Smooth PME is opt-in** ([`PeriodicForceField::with_mesh`]); classical Ewald is the default.
//! **The neighbour list is opt-in** ([`PeriodicForceField::with_neighbour_list`],
//! [`crate::neighbours`]); without it the cell list is rebuilt at each evaluation, `O(N)`, and
//! with it the bits are the same. **No triclinic box, no barostat**; constraints are the dynamics' —
//! [`crate::water`]'s SETTLE and [`crate::shake`]'s bonds — not the force field's. **No implicit
//! solvent in the box**, and no QEq for a periodic system: the charges are given.

use crate::alchemy::{Alchemical, AlchemyError, Coupling, Lambda, SoftCore};
use crate::boresch::Boresch;
use crate::ccd::Component;
use crate::dynamics::Potential;
use crate::energy::{coulomb, ForceField, Pair, Unsupported};
use crate::ewald::{add_pair, norm, sub, Ewald, EwaldEnergy, EwaldParameters};
use crate::neighbours::NeighbourList;
use crate::uff::UffType;
use crate::water;
use std::f64::consts::PI;
use std::fmt;
use std::ops::Range;

/// How two atoms' van der Waals combine into their pair's: see the module documentation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Mixing {
    /// UFF's eqs 21b and 22: `x_ij = √(x_i x_j)`, `D_ij = √(D_i D_j)`. The default.
    #[default]
    Geometric,
    /// Lorentz–Berthelot, as AMBER and GAFF combine: `x_ij = (x_i + x_j)/2`, `D_ij = √(D_i D_j)`.
    LorentzBerthelot,
}

/// Why [`PeriodicForceField::with_nonbonded`] refused its values.
#[derive(Clone, Debug, PartialEq)]
pub enum NonbondedError {
    /// The run of atoms is not inside the force field's: `start..end` against `atoms` atoms.
    OutOfRange {
        /// The run's first atom.
        start: usize,
        /// One past its last.
        end: usize,
        /// How many atoms the force field has.
        atoms: usize,
    },
    /// The values are not one per atom of the run.
    Lengths {
        /// How many atoms the run has.
        atoms: usize,
        /// How many values were given.
        values: usize,
    },
    /// An atom's `x` or `D` is not finite.
    NotFinite {
        /// The atom.
        atom: usize,
    },
    /// An atom's `x` is negative, or zero with a well: only `[0, 0]` may have `x = 0`.
    Distance {
        /// The atom.
        atom: usize,
        /// Its `x`, metres.
        value: f64,
    },
    /// An atom's `D` is negative.
    Well {
        /// The atom.
        atom: usize,
        /// Its `D`, joules.
        value: f64,
    },
}

impl fmt::Display for NonbondedError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NonbondedError::OutOfRange { start, end, atoms } => write!(
                f,
                "atoms {start}..{end} are not inside the force field's {atoms} atoms"
            ),
            NonbondedError::Lengths { atoms, values } => write!(
                f,
                "{values} van der Waals values for {atoms} atoms: they must be one per atom"
            ),
            NonbondedError::NotFinite { atom } => {
                write!(f, "atom {atom}'s van der Waals is not finite")
            }
            NonbondedError::Distance { atom, value } => write!(
                f,
                "atom {atom}'s van der Waals distance x = {value:e} m is negative, or zero with a \
                 well: only an atom with no van der Waals at all may have x = 0"
            ),
            NonbondedError::Well { atom, value } => write!(
                f,
                "atom {atom}'s van der Waals well D = {value:e} J is negative"
            ),
        }
    }
}

impl std::error::Error for NonbondedError {}

/// Scaled 1-4 pairs ([`PeriodicForceField::with_one_four`]).
#[derive(Clone, Debug, PartialEq)]
struct OneFour {
    /// `[s_v, s_c]`, not both 1.
    scale: [f64; 2],
    /// Every 1-4 pair once, `[i, j]` with `i < j`, sorted.
    pairs: Vec<[usize; 2]>,
    /// The excluded pairs and the 1-4 pairs: what the pair search leaves out.
    skip: Exclusions,
}

/// The pairs three bonds apart that are neither 1-2 nor 1-3, once each, sorted.
fn one_four_pairs(n: usize, bonds: &[[usize; 2]], exclusions: &Exclusions) -> Vec<[usize; 2]> {
    let mut adjacent = vec![Vec::new(); n];
    for &[i, j] in bonds {
        adjacent[i].push(j);
        adjacent[j].push(i);
    }
    let mut out = Vec::new();
    for &[j, k] in bonds {
        for &i in &adjacent[j] {
            if i == k {
                continue;
            }
            for &l in &adjacent[k] {
                if l == j || l == i || exclusions.contains(i, l) {
                    continue;
                }
                out.push([i.min(l), i.max(l)]);
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// What the long-range correction needs of a set of atoms: under geometric mixing the two sums
/// `[Σ √D x⁶, Σ √D x³]`, and under Lorentz–Berthelot each class of identical `[x, D]` with its
/// count, in order of first appearance.
enum TailSet {
    Sums([f64; 2]),
    Classes(Vec<([f64; 2], f64)>),
}

/// `[Σ_{i∈A} Σ_{j∈B} D_ij x_ij¹², Σ_{i∈A} Σ_{j∈B} D_ij x_ij⁶]` for two sets summarised under one
/// mixing rule.
fn tail_products(a: &TailSet, b: &TailSet) -> [f64; 2] {
    match (a, b) {
        (TailSet::Sums(a), TailSet::Sums(b)) => [a[0] * b[0], a[1] * b[1]],
        (TailSet::Classes(a), TailSet::Classes(b)) => {
            let mut s = [0.0; 2];
            for &([xa, da], na) in a {
                for &([xb, db], nb) in b {
                    let x = 0.5 * (xa + xb);
                    let x2 = x * x;
                    let x6 = x2 * x2 * x2;
                    let w = na * nb * (da * db).sqrt();
                    s[0] += w * x6 * x6;
                    s[1] += w * x6;
                }
            }
            s
        }
        _ => unreachable!("both sets are summarised under the field's one rule"),
    }
}

/// An orthorhombic periodic box: three edges at right angles, metres. See the module documentation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PeriodicBox {
    lengths: [f64; 3],
}

impl PeriodicBox {
    /// The box with edges `lengths` (metres) along x, y and z.
    ///
    /// # Panics
    ///
    /// If an edge is not positive and finite.
    pub fn new(lengths: [f64; 3]) -> PeriodicBox {
        assert!(
            lengths.iter().all(|l| l.is_finite() && *l > 0.0),
            "a box's edges must be positive and finite"
        );
        PeriodicBox { lengths }
    }

    /// The cube of edge `side` metres.
    pub fn cubic(side: f64) -> PeriodicBox {
        PeriodicBox::new([side; 3])
    }

    /// The three edges, metres.
    pub fn lengths(&self) -> [f64; 3] {
        self.lengths
    }

    /// The volume, cubic metres.
    pub fn volume(&self) -> f64 {
        self.lengths[0] * self.lengths[1] * self.lengths[2]
    }

    /// Half the shortest edge: the longest cutoff a pair has one image inside.
    pub fn half_shortest_edge(&self) -> f64 {
        0.5 * self.lengths[0].min(self.lengths[1]).min(self.lengths[2])
    }

    /// The same box with each edge multiplied by its factor.
    pub fn scaled(&self, factors: [f64; 3]) -> PeriodicBox {
        PeriodicBox::new([0, 1, 2].map(|a| self.lengths[a] * factors[a]))
    }

    /// The image of displacement `d` (metres) nearest zero: `d − L round(d/L)` per axis, each
    /// component in `[−L/2, L/2]`.
    pub fn minimum_image(&self, d: [f64; 3]) -> [f64; 3] {
        [0, 1, 2].map(|a| {
            let l = self.lengths[a];
            d[a] - l * (d[a] / l).round()
        })
    }

    /// Point `p` (metres) moved by whole edges into `[0, L)` per axis.
    pub fn wrap(&self, p: [f64; 3]) -> [f64; 3] {
        [0, 1, 2].map(|a| {
            let l = self.lengths[a];
            let w = p[a] - l * (p[a] / l).floor();
            // A point a rounding below an edge can land on L itself.
            if w >= l {
                w - l
            } else {
                w
            }
        })
    }
}

/// The excluded pairs of a system, in a compressed table by atom: each atom's partners sorted.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Exclusions {
    start: Vec<usize>,
    partner: Vec<usize>,
    pairs: Vec<[usize; 2]>,
}

impl Exclusions {
    /// From a list of pairs in any order and orientation, repeats allowed.
    pub(crate) fn from_pairs(n: usize, pairs: &[[usize; 2]]) -> Exclusions {
        let mut sorted: Vec<[usize; 2]> = pairs
            .iter()
            .map(|&[i, j]| {
                assert!(
                    i < n && j < n && i != j,
                    "an excluded pair must join two atoms"
                );
                [i.min(j), i.max(j)]
            })
            .collect();
        sorted.sort_unstable();
        sorted.dedup();
        let mut lists = vec![Vec::new(); n];
        for &[i, j] in &sorted {
            lists[i].push(j);
            lists[j].push(i);
        }
        let mut start = Vec::with_capacity(n + 1);
        let mut partner = Vec::new();
        start.push(0);
        for mut l in lists {
            l.sort_unstable();
            partner.extend(l);
            start.push(partner.len());
        }
        Exclusions {
            start,
            partner,
            pairs: sorted,
        }
    }

    /// Whether `i` and `j` are excluded.
    pub(crate) fn contains(&self, i: usize, j: usize) -> bool {
        self.partner[self.start[i]..self.start[i + 1]]
            .binary_search(&j)
            .is_ok()
    }

    /// Every excluded pair once, `[i, j]` with `i < j`, sorted.
    pub(crate) fn pairs(&self) -> &[[usize; 2]] {
        &self.pairs
    }
}

/// The most cells a pair search makes: 2²² (4 194 304, 34 MB of offsets). Past it every axis has
/// at most 128 cells, which are then wider than half the cutoff and so still hold every pair, in
/// another order. Uncapped, a cutoff of 0.01 Å asked for 83.6 GB and the process aborted on the
/// allocation. No system this crate has built comes near the cap: W4's box, 62.35 × 65.47 ×
/// 74.82 Å at 9 Å, is 13 × 14 × 16 = 2 912 cells, so every existing search keeps its order.
pub(crate) const MOST_CELLS: f64 = 4_194_304.0;

/// How many cells the pairs nearer than one cell apart are searched over.
pub(crate) const REACH: usize = 2;

/// The cells along each axis of a pair search at `cutoff` in a box of edges `l`: as many as are at
/// least half the cutoff wide, capped at 128 an axis when there would be more than
/// [`MOST_CELLS`]; `None` when an axis has fewer than five, and every pair is visited instead.
pub(crate) fn cell_grid(l: [f64; 3], cutoff: f64) -> Option<[usize; 3]> {
    // Cells at least half the cutoff wide, and the 5³ about each. There is no 3³ branch of
    // cutoff-wide cells: three of those across is six half-cells, so it would never be taken.
    let across = l.map(|x| (2.0 * x / cutoff).floor());
    let mut cells = across.map(|c| c as usize);
    if across[0] * across[1] * across[2] > MOST_CELLS {
        cells = cells.map(|c| c.min(128));
    }
    if cells.iter().any(|&c| c < 2 * REACH + 1) {
        None
    } else {
        Some(cells)
    }
}

/// The cell of wrapped point `p` in a box of edges `l` cut into `cells`.
#[inline(always)]
pub(crate) fn cell_of(p: [f64; 3], l: [f64; 3], cells: [usize; 3]) -> [usize; 3] {
    [0, 1, 2].map(|a| ((p[a] / l[a] * cells[a] as f64) as usize).min(cells[a] - 1))
}

/// The minimum-image separation `w_i − w_j` of two wrapped points and its square: `d ∓ L` on an
/// axis where `|d| > L/2`.
#[inline(always)]
pub(crate) fn separation(
    wi: [f64; 3],
    wj: [f64; 3],
    l: [f64; 3],
    half: [f64; 3],
) -> ([f64; 3], f64) {
    let mut d = sub(wi, wj);
    for a in 0..3 {
        if d[a] > half[a] {
            d[a] -= l[a];
        } else if d[a] < -half[a] {
            d[a] += l[a];
        }
    }
    let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    (d, r2)
}

/// Calls `f(i, j, d, r²)` for every pair `i < j` whose minimum-image separation `d = r_i − r_j`
/// is shorter than `cutoff`, once each, in an order fixed by the positions' bits: by `i`, then by
/// the cells about `i`'s in a fixed order, then by `j` within a cell. A cell list: cells at least
/// half the cutoff wide and the 125 about each when the box is at least five such cells across
/// every axis, and every pair when it is fewer ([`cell_grid`]).
///
/// The minimum image of two wrapped points is `d ∓ L` when `|d| > L/2`: a comparison, not a
/// rounding, and exact. A pair exactly half a box apart keeps `d`, which is as near as `d − L`.
pub(crate) fn for_each_pair(
    cell: &PeriodicBox,
    at: &[[f64; 3]],
    cutoff: f64,
    mut f: impl FnMut(usize, usize, [f64; 3], f64),
) {
    let n = at.len();
    let l = cell.lengths();
    let half = l.map(|x| 0.5 * x);
    let w: Vec<[f64; 3]> = at.iter().map(|p| cell.wrap(*p)).collect();
    let rc2 = cutoff * cutoff;
    let mut visit = |i: usize, j: usize| {
        let (d, r2) = separation(w[i], w[j], l, half);
        if r2 < rc2 {
            f(i, j, d, r2);
        }
    };
    let reach = REACH;
    let Some(cells) = cell_grid(l, cutoff) else {
        for i in 0..n {
            for j in i + 1..n {
                visit(i, j);
            }
        }
        return;
    };
    let flat = |c: [usize; 3]| (c[0] * cells[1] + c[1]) * cells[2] + c[2];
    let total = cells[0] * cells[1] * cells[2];
    let home: Vec<[usize; 3]> = w.iter().map(|p| cell_of(*p, l, cells)).collect();
    let mut start = vec![0usize; total + 1];
    for c in &home {
        start[flat(*c) + 1] += 1;
    }
    for k in 0..total {
        start[k + 1] += start[k];
    }
    let mut fill = start.clone();
    let mut members = vec![0usize; n];
    for (i, c) in home.iter().enumerate() {
        let k = flat(*c);
        members[fill[k]] = i;
        fill[k] += 1;
    }
    let span = 2 * reach + 1;
    for (i, &c) in home.iter().enumerate() {
        for dx in 0..span {
            for dy in 0..span {
                for dz in 0..span {
                    let nb = [
                        (c[0] + cells[0] + dx - reach) % cells[0],
                        (c[1] + cells[1] + dy - reach) % cells[1],
                        (c[2] + cells[2] + dz - reach) % cells[2],
                    ];
                    let k = flat(nb);
                    for &j in &members[start[k]..start[k + 1]] {
                        if j > i {
                            visit(i, j);
                        }
                    }
                }
            }
        }
    }
}

/// `(2π/V) Σ_i Σ_j D_ij [x_ij¹²/(9 r_c⁹) − 2 x_ij⁶/(3 r_c³)]` for two sets of atoms given by
/// [`tail_products`]' `p = [Σ Σ D_ij x_ij¹², Σ Σ D_ij x_ij⁶]`: the long-range correction of every
/// ordered pair with one atom in each, joules. (Under geometric mixing `p` is a product of two sums
/// and `2 p₁` is `(2 a₁) b₁` exactly, a power of two scaling, so the bits are the factored form's.)
fn tail_energy(p: [f64; 2], cutoff: f64, volume: f64) -> f64 {
    let rc3 = cutoff * cutoff * cutoff;
    let rc9 = rc3 * rc3 * rc3;
    2.0 * PI / volume * (p[0] / (9.0 * rc9) - 2.0 * p[1] / (3.0 * rc3))
}

/// A&T's `P_tail` for the same products, pascals.
fn tail_pressure(p: [f64; 2], cutoff: f64, volume: f64) -> f64 {
    let rc3 = cutoff * cutoff * cutoff;
    let rc9 = rc3 * rc3 * rc3;
    2.0 * PI / (3.0 * volume * volume) * (4.0 * p[0] / (3.0 * rc9) - 4.0 * p[1] / rc3)
}

/// A periodic system's energy, by term, joules per molecule. `total` is accumulated separately.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PeriodicEnergy {
    /// Bond stretch.
    pub bond: f64,
    /// Angle bend.
    pub angle: f64,
    /// Torsion.
    pub torsion: f64,
    /// Inversion.
    pub inversion: f64,
    /// Van der Waals, the pairs inside the cutoff and any scaled 1-4 pairs (and, in a
    /// [`PeriodicDecoupling`], the group's own in vacuum and its soft-core pairs with the rest).
    pub van_der_waals: f64,
    /// The long-range van der Waals correction, `E_tail`.
    pub dispersion_correction: f64,
    /// Electrostatics: [`PeriodicEnergy::ewald`]'s total (and, in a [`PeriodicDecoupling`], the
    /// group's own Coulomb pairs in vacuum).
    pub electrostatic: f64,
    /// The Ewald sum by part; in a [`PeriodicDecoupling`], the rest's and the cross terms at λ.
    pub ewald: EwaldEnergy,
    /// The bonded terms, van der Waals, its correction and electrostatics: everything above but
    /// [`PeriodicEnergy::ewald`], whose total is in `electrostatic`.
    pub total: f64,
}

/// An energy, the force on every atom (newtons), and the virial `W_ab = −∂U/∂ε_ab` (joules), so
/// that `P_ab V = Σ m v_a v_b + W_ab`.
#[derive(Clone, Debug, PartialEq)]
pub struct PeriodicEvaluation {
    /// The energy.
    pub energy: PeriodicEnergy,
    /// The force on each atom, in the component's atom order.
    pub forces: Vec<[f64; 3]>,
    /// The virial: the strain derivative of exactly this energy, the cutoff's steps aside.
    pub virial: [[f64; 3]; 3],
}

/// UFF in an orthorhombic periodic box with Ewald electrostatics: see the module documentation.
#[derive(Clone, Debug, PartialEq)]
pub struct PeriodicForceField {
    /// The solute's bonded terms, on its atoms — the first `solute` — or none for water alone.
    bonded: Option<ForceField>,
    solute: usize,
    /// The rigid TIP3P waters after the solute, oxygen first.
    waters: Vec<[usize; 3]>,
    /// How many monatomic ions follow the waters ([`PeriodicForceField::with_ions`]).
    ions: usize,
    ewald: Ewald,
    charges: Vec<f64>,
    /// `[x_i, D_i]`, metres and joules.
    vdw: Vec<[f64; 2]>,
    /// How the pairs combine them.
    mixing: Mixing,
    bonds: Vec<[usize; 2]>,
    pub(crate) exclusions: Exclusions,
    /// The scaled 1-4 pairs, or `None` at the default `(1, 1)`.
    one_four: Option<OneFour>,
    /// `[atom, the atom it is reached from]`, breadth first, roots left out.
    tree: Vec<[usize; 2]>,
    /// The real-space pairs' Verlet list, when one is asked for ([`crate::neighbours`]).
    pub(crate) neighbours: Option<NeighbourList>,
}

/// What a [`PeriodicDecoupling`] gives the shared evaluation.
struct Split<'a> {
    group: &'a [bool],
    lambda: Lambda,
    soft_core: SoftCore,
    own: &'a [Pair],
    /// Each own pair's `[s_v, s_c]`: `[1, 1]` but for its 1-4 pairs.
    own_scale: &'a [[f64; 2]],
}

/// The shared evaluation's result.
struct Core {
    energy: PeriodicEnergy,
    forces: Vec<[f64; 3]>,
    virial: [[f64; 3]; 3],
    /// The λ-dependent part, the restraint aside, summed in the evaluation's own order: what the
    /// tests hold [`PeriodicDecoupling`]'s own `couplings` to.
    #[cfg_attr(not(test), allow(dead_code))]
    coupling: Coupling,
}

impl PeriodicForceField {
    /// UFF for `component` typed as `types`, in `cell`, with the cutoff and Ewald sum of
    /// `parameters` ([`EwaldParameters::for_accuracy`]) and every charge zero.
    ///
    /// # Errors
    ///
    /// As [`ForceField::new`].
    ///
    /// # Panics
    ///
    /// If `types` is not one per atom, or the cutoff is past half the box's shortest edge.
    pub fn new(
        component: &Component,
        types: &[UffType],
        cell: PeriodicBox,
        parameters: EwaldParameters,
    ) -> Result<PeriodicForceField, Unsupported> {
        PeriodicForceField::assemble(Some((component, types)), 0, cell, parameters)
    }

    /// Rigid TIP3P water alone ([`crate::water`]): `waters` molecules, atoms in the order oxygen,
    /// hydrogen, hydrogen, water by water, with Table I's charges and the O–O Lennard-Jones term
    /// as a UFF pair, and every pair inside a water excluded. No bonded term: the geometry is
    /// SETTLE's to hold ([`crate::water::Settle`]).
    ///
    /// # Panics
    ///
    /// If the cutoff is past half the box's shortest edge.
    pub fn tip3p(
        cell: PeriodicBox,
        parameters: EwaldParameters,
        waters: usize,
    ) -> PeriodicForceField {
        PeriodicForceField::assemble(None, waters, cell, parameters)
            .expect("water alone has nothing unsupported")
    }

    /// `component`, typed as `types` and with zero charges until
    /// [`PeriodicForceField::with_solute_charges`], followed by `waters` rigid TIP3P waters, as
    /// [`PeriodicForceField::tip3p`] builds them. A water oxygen and a UFF atom combine by UFF's
    /// geometric rule — a judgement, see [`crate::water`].
    ///
    /// # Errors
    ///
    /// As [`ForceField::new`].
    ///
    /// # Panics
    ///
    /// As [`PeriodicForceField::new`].
    pub fn solvated(
        component: &Component,
        types: &[UffType],
        waters: usize,
        cell: PeriodicBox,
        parameters: EwaldParameters,
    ) -> Result<PeriodicForceField, Unsupported> {
        PeriodicForceField::assemble(Some((component, types)), waters, cell, parameters)
    }

    fn assemble(
        solute: Option<(&Component, &[UffType])>,
        waters: usize,
        cell: PeriodicBox,
        parameters: EwaldParameters,
    ) -> Result<PeriodicForceField, Unsupported> {
        let (bonded, mut bonds, mut vdw, n_solute) = match solute {
            Some((component, types)) => {
                let bonded = ForceField::bonded_only(component, types)?;
                let bonds: Vec<[usize; 2]> = component.bonds().iter().map(|b| b.atoms).collect();
                let vdw: Vec<[f64; 2]> = types
                    .iter()
                    .map(|t| {
                        let p = t.parameters();
                        [p.vdw_distance, p.vdw_energy]
                    })
                    .collect();
                (Some(bonded), bonds, vdw, component.atoms().len())
            }
            None => (None, Vec::new(), Vec::new(), 0),
        };
        let mut charges = vec![0.0; n_solute];
        let mut rigid = Vec::with_capacity(waters);
        for w in 0..waters {
            let o = n_solute + 3 * w;
            bonds.push([o, o + 1]);
            bonds.push([o, o + 2]);
            vdw.push(water::oxygen_vdw());
            vdw.push([0.0, 0.0]);
            vdw.push([0.0, 0.0]);
            charges.extend(water::CHARGES);
            rigid.push([o, o + 1, o + 2]);
        }
        let n = n_solute + 3 * waters;
        let mut adjacent = vec![Vec::new(); n];
        for &[i, j] in &bonds {
            adjacent[i].push(j);
            adjacent[j].push(i);
        }
        let mut excluded = Vec::new();
        for (j, around) in adjacent.iter().enumerate() {
            for (p, &i) in around.iter().enumerate() {
                excluded.push([i, j]);
                for &k in &around[p + 1..] {
                    if k != i {
                        excluded.push([i, k]);
                    }
                }
            }
        }
        let mut seen = vec![false; n];
        let mut tree = Vec::new();
        for root in 0..n {
            if seen[root] {
                continue;
            }
            seen[root] = true;
            let mut queue = std::collections::VecDeque::from([root]);
            while let Some(a) = queue.pop_front() {
                for &b in &adjacent[a] {
                    if !seen[b] {
                        seen[b] = true;
                        tree.push([b, a]);
                        queue.push_back(b);
                    }
                }
            }
        }
        Ok(PeriodicForceField {
            bonded,
            solute: n_solute,
            waters: rigid,
            ions: 0,
            ewald: Ewald::new(cell, parameters),
            charges,
            vdw,
            mixing: Mixing::Geometric,
            bonds,
            exclusions: Exclusions::from_pairs(n, &excluded),
            one_four: None,
            tree,
            neighbours: None,
        })
    }

    /// The same force field with the van der Waals of atoms `atoms` set to `van_der_waals`, one
    /// `[x, D]` an atom in UFF's form — `x` the distance of the minimum of the atom's pair with
    /// itself, metres, and `D` its depth, joules; AMBER's `R*/2` and `ε` are `x = 2 (R*/2)` and
    /// `D = ε` — for any atoms, a bonded solute's, a water's or an ion's. See the module
    /// documentation. Everything else is unchanged; a [`PeriodicDecoupling`] takes the values its
    /// field has when it is built.
    ///
    /// # Errors
    ///
    /// [`NonbondedError`], naming what is wrong: the run is not inside the atoms, the values are not
    /// one per atom of it, or an atom's values are not finite, its `x` negative or zero with a
    /// well, or its `D` negative. `[0, 0]`, no van der Waals at all, is taken. Nothing is changed
    /// then.
    pub fn with_nonbonded(
        mut self,
        atoms: Range<usize>,
        van_der_waals: Vec<[f64; 2]>,
    ) -> Result<PeriodicForceField, NonbondedError> {
        let n = self.vdw.len();
        if atoms.start > atoms.end || atoms.end > n {
            return Err(NonbondedError::OutOfRange {
                start: atoms.start,
                end: atoms.end,
                atoms: n,
            });
        }
        if van_der_waals.len() != atoms.len() {
            return Err(NonbondedError::Lengths {
                atoms: atoms.len(),
                values: van_der_waals.len(),
            });
        }
        for (k, &[x, d]) in van_der_waals.iter().enumerate() {
            let atom = atoms.start + k;
            if !x.is_finite() || !d.is_finite() {
                return Err(NonbondedError::NotFinite { atom });
            }
            if x < 0.0 || (x == 0.0 && d != 0.0) {
                return Err(NonbondedError::Distance { atom, value: x });
            }
            if d < 0.0 {
                return Err(NonbondedError::Well { atom, value: d });
            }
        }
        self.vdw[atoms].copy_from_slice(&van_der_waals);
        Ok(self)
    }

    /// The same force field with its pairs' van der Waals combined by `mixing` ([`Mixing`]):
    /// every cut-off pair, every 1-4 pair, the long-range correction, and a decoupling's soft core,
    /// cross pairs and own pairs. [`Mixing::Geometric`], UFF's, is the default.
    pub fn with_mixing(mut self, mixing: Mixing) -> PeriodicForceField {
        self.mixing = mixing;
        self
    }

    /// How the pairs' van der Waals combine.
    pub fn mixing(&self) -> Mixing {
        self.mixing
    }

    /// The same force field with every 1-4 pair's van der Waals scaled by `van_der_waals` and its
    /// Coulomb energy by `coulomb` — GAFF's are `½` and `1/1.2` — the Coulomb term in Ewald as a
    /// partial exclusion. 1-2 and 1-3 pairs stay excluded. `(1, 1)`, UFF's, is the default and
    /// leaves every 1-4 pair a pair like any other, to the bit. See the module documentation. A
    /// neighbour list starts again empty, with the same skin.
    ///
    /// # Panics
    ///
    /// If a scale is negative or not finite.
    pub fn with_one_four(mut self, van_der_waals: f64, coulomb: f64) -> PeriodicForceField {
        assert!(
            [van_der_waals, coulomb]
                .iter()
                .all(|s| s.is_finite() && *s >= 0.0),
            "a 1-4 scale must be finite and not negative: {van_der_waals}, {coulomb}"
        );
        self.one_four = if van_der_waals == 1.0 && coulomb == 1.0 {
            None
        } else {
            let n = self.charges.len();
            let pairs = one_four_pairs(n, &self.bonds, &self.exclusions);
            let mut all = self.exclusions.pairs().to_vec();
            all.extend_from_slice(&pairs);
            Some(OneFour {
                scale: [van_der_waals, coulomb],
                pairs,
                skip: Exclusions::from_pairs(n, &all),
            })
        };
        self.neighbours = self.neighbours.as_ref().map(NeighbourList::emptied);
        self
    }

    /// The 1-4 scales `[van der Waals, Coulomb]`, `[1, 1]` by default.
    pub fn one_four(&self) -> [f64; 2] {
        self.one_four.as_ref().map_or([1.0, 1.0], |o| o.scale)
    }

    /// Every 1-4 pair that is scaled, once, `[i, j]` with `i < j`, sorted: none at the default
    /// scales, where a 1-4 pair is a pair like any other.
    pub fn one_four_pairs(&self) -> &[[usize; 2]] {
        self.one_four.as_ref().map_or(&[], |o| &o.pairs)
    }

    /// Each atom's van der Waals `[x, D]`, metres and joules, in UFF's form: its type's, or what
    /// [`PeriodicForceField::with_nonbonded`] set.
    pub fn van_der_waals(&self) -> &[[f64; 2]] {
        &self.vdw
    }

    /// The pairs the pair search leaves out: the excluded ones, and the scaled 1-4 pairs, which are
    /// evaluated on their own.
    pub(crate) fn skipped(&self) -> &Exclusions {
        self.one_four.as_ref().map_or(&self.exclusions, |o| &o.skip)
    }

    /// The same force field with partial charges `charges`, elementary charges.
    ///
    /// # Panics
    ///
    /// If `charges` is not one per atom.
    pub fn with_charges(mut self, charges: Vec<f64>) -> PeriodicForceField {
        assert_eq!(charges.len(), self.charges.len(), "one charge per atom");
        self.charges = charges;
        self
    }

    /// The same force field with the solute's partial charges `charges`, elementary charges, the
    /// waters' left as TIP3P's.
    ///
    /// # Panics
    ///
    /// If `charges` is not one per solute atom.
    pub fn with_solute_charges(mut self, charges: Vec<f64>) -> PeriodicForceField {
        assert_eq!(charges.len(), self.solute, "one charge per solute atom");
        self.charges[..self.solute].copy_from_slice(&charges);
        self
    }

    /// The same force field with `ions` appended after every atom it has, in order: each a point
    /// charge with a UFF-form van der Waals term `[x, D]` ([`crate::solvated::Ion`]), combined with
    /// every other atom by UFF's geometric rule as a water oxygen is, bonded to nothing, excluded
    /// from nothing. Counter-ions for a charged solute ([`crate::solvated`]). Every term is the
    /// force field's own — the cell list, Ewald and its background, the long-range correction —
    /// so an ion is one more atom and nothing else changes.
    pub fn with_ions(mut self, ions: &[crate::solvated::Ion]) -> PeriodicForceField {
        for ion in ions {
            self.charges.push(ion.charge);
            self.vdw.push([ion.vdw_distance, ion.vdw_energy]);
        }
        self.ions += ions.len();
        self.exclusions = Exclusions::from_pairs(self.charges.len(), self.exclusions.pairs());
        if let Some(o) = &mut self.one_four {
            o.skip = Exclusions::from_pairs(self.charges.len(), o.skip.pairs());
        }
        self.neighbours = self.neighbours.as_ref().map(NeighbourList::emptied);
        self
    }

    /// How many ions follow the waters: the last atoms.
    pub fn ion_atoms(&self) -> usize {
        self.ions
    }

    /// The rigid TIP3P waters, atom indices, oxygen first: what [`crate::water::Settle::tip3p`]
    /// is built from.
    pub fn rigid_waters(&self) -> &[[usize; 3]] {
        &self.waters
    }

    /// How many atoms the solute has: the first ones, before the waters.
    pub fn solute_atoms(&self) -> usize {
        self.solute
    }

    /// The same force field with the Ewald sum's reciprocal part by smooth particle-mesh Ewald on
    /// `mesh`'s order and grid ([`Ewald::with_mesh`], [`crate::pme`]): every other term, the
    /// decoupling's cross terms included, is unchanged.
    ///
    /// # Panics
    ///
    /// As [`crate::Pme::new`].
    pub fn with_mesh(mut self, mesh: crate::pme::PmeParameters) -> PeriodicForceField {
        self.ewald = self.ewald.with_mesh(mesh);
        self
    }

    /// The same force field in another box, with the same α, cutoff and integer wave vectors — or
    /// mesh — ([`Ewald::with_cell`]): for a derivative with respect to the box. A neighbour list
    /// starts again empty, with the same skin.
    ///
    /// # Panics
    ///
    /// If the cutoff is past half the new box's shortest edge.
    pub fn with_cell(&self, cell: PeriodicBox) -> PeriodicForceField {
        let mut f = self.clone();
        f.ewald = self.ewald.with_cell(cell);
        f.neighbours = self.neighbours.as_ref().map(NeighbourList::emptied);
        f
    }

    /// The box.
    pub fn cell(&self) -> PeriodicBox {
        self.ewald.cell()
    }

    /// The Ewald sum.
    pub fn ewald(&self) -> &Ewald {
        &self.ewald
    }

    /// The partial charges, elementary charges.
    pub fn charges(&self) -> &[f64] {
        &self.charges
    }

    /// Every excluded pair — 1-2 and 1-3 — once, `[i, j]` with `i < j`, sorted.
    pub fn excluded_pairs(&self) -> &[[usize; 2]] {
        self.exclusions.pairs()
    }

    /// The solute's bonded terms, a [`ForceField`] with no non-bonded pair on the first
    /// [`PeriodicForceField::solute_atoms`]; `None` for water alone.
    pub fn bonded(&self) -> Option<&ForceField> {
        self.bonded.as_ref()
    }

    /// The positions `at` with every molecule made whole: see the module documentation. An atom
    /// already at the image nearest the atom it is reached from keeps its bits.
    ///
    /// # Panics
    ///
    /// If `at` is not one position per atom.
    pub fn whole(&self, at: &[[f64; 3]]) -> Vec<[f64; 3]> {
        assert_eq!(at.len(), self.charges.len(), "one position per atom");
        let l = self.cell().lengths();
        let mut shift = vec![[0.0f64; 3]; at.len()];
        let mut out = at.to_vec();
        for &[b, a] in &self.tree {
            for k in 0..3 {
                let jump = (at[b][k] - at[a][k]) / l[k];
                shift[b][k] = shift[a][k] - l[k] * jump.round();
                if shift[b][k] != 0.0 {
                    out[b][k] = at[b][k] + shift[b][k];
                }
            }
        }
        out
    }

    /// `[Σ √D_i x_i⁶, Σ √D_i x_i³]` over the atoms `keep` accepts.
    fn tail_sums(&self, keep: impl Fn(usize) -> bool) -> [f64; 2] {
        let mut s = [0.0; 2];
        for (i, &[x, d]) in self.vdw.iter().enumerate() {
            if keep(i) {
                let x3 = x * x * x;
                s[0] += d.sqrt() * x3 * x3;
                s[1] += d.sqrt() * x3;
            }
        }
        s
    }

    /// The atoms `keep` accepts, summarised for the long-range correction under this field's
    /// mixing: their two sums, or their classes. See the module documentation.
    fn tail_set(&self, keep: impl Fn(usize) -> bool) -> TailSet {
        match self.mixing {
            Mixing::Geometric => TailSet::Sums(self.tail_sums(keep)),
            Mixing::LorentzBerthelot => {
                let mut classes: Vec<([f64; 2], f64)> = Vec::new();
                for (i, v) in self.vdw.iter().enumerate() {
                    if !keep(i) {
                        continue;
                    }
                    let bits = v.map(f64::to_bits);
                    match classes
                        .iter_mut()
                        .find(|(c, _)| c.map(f64::to_bits) == bits)
                    {
                        Some(class) => class.1 += 1.0,
                        None => classes.push((*v, 1.0)),
                    }
                }
                TailSet::Classes(classes)
            }
        }
    }

    /// The long-range van der Waals correction `E_tail`, joules: see the module documentation.
    pub fn dispersion_correction(&self) -> f64 {
        let s = self.tail_set(|_| true);
        tail_energy(
            tail_products(&s, &s),
            self.ewald.parameters().cutoff,
            self.cell().volume(),
        )
    }

    /// Allen and Tildesley's long-range correction to the pressure, `P_tail`, pascals: see the
    /// module documentation for how it differs from the virial's `E_tail/V`.
    pub fn pressure_correction(&self) -> f64 {
        let s = self.tail_set(|_| true);
        tail_pressure(
            tail_products(&s, &s),
            self.ewald.parameters().cutoff,
            self.cell().volume(),
        )
    }

    /// The van der Waals pair of atoms `i` and `j`, combined by the field's [`Mixing`]: under the
    /// default, as [`Pair::new`] combines two types.
    fn pair(&self, i: usize, j: usize) -> Pair {
        let ([xi, di], [xj, dj]) = (self.vdw[i], self.vdw[j]);
        Pair {
            atoms: [i, j],
            distance: match self.mixing {
                Mixing::Geometric => (xi * xj).sqrt(),
                Mixing::LorentzBerthelot => 0.5 * (xi + xj),
            },
            well: (di * dj).sqrt(),
        }
    }

    /// The energy at `at` (metres, any image), without the forces' use.
    ///
    /// # Panics
    ///
    /// If `at` is not one position per atom.
    pub fn energy(&self, at: &[[f64; 3]]) -> PeriodicEnergy {
        self.core(at, None, false).energy
    }

    /// The energy at `at` (metres, any image), the force on every atom and the virial.
    ///
    /// # Panics
    ///
    /// If `at` is not one position per atom.
    pub fn evaluate(&self, at: &[[f64; 3]]) -> PeriodicEvaluation {
        let c = self.core(at, None, true);
        PeriodicEvaluation {
            energy: c.energy,
            forces: c.forces,
            virial: c.virial,
        }
    }

    /// The whole evaluation, with or without a decoupled group.
    fn core(&self, at: &[[f64; 3]], split: Option<&Split>, forces_wanted: bool) -> Core {
        let n = self.charges.len();
        assert_eq!(at.len(), n, "one position per atom");
        let cell = self.cell();
        let parameters = self.ewald.parameters();
        let (cutoff, volume) = (parameters.cutoff, cell.volume());
        let marked = |i: usize| split.is_some_and(|s| s.group[i]);
        let (le, lv) = split.map_or((1.0, 1.0), |s| {
            (s.lambda.electrostatics, s.lambda.van_der_waals)
        });
        let mut e = PeriodicEnergy::default();
        let mut virial = [[0.0; 3]; 3];
        let mut coupling = Coupling {
            energy: 0.0,
            gradient: [0.0; 3],
        };

        // The bonded terms, on whole molecules.
        let whole = self.whole(at);
        let mut forces = vec![[0.0; 3]; n];
        if let Some(bonded) = &self.bonded {
            let b = bonded.evaluate(&whole[..self.solute]);
            e.bond = b.energy.bond;
            e.total += b.energy.bond;
            e.angle = b.energy.angle;
            e.total += b.energy.angle;
            e.torsion = b.energy.torsion;
            e.total += b.energy.torsion;
            e.inversion = b.energy.inversion;
            e.total += b.energy.inversion;
            forces[..self.solute].copy_from_slice(&b.forces);
        }
        for (r, f) in whole.iter().zip(&forces) {
            for a in 0..3 {
                for c in 0..3 {
                    virial[a][c] += r[a] * f[c];
                }
            }
        }

        // The pairs inside the cutoff: by the cell list, or by the neighbour list in the cell
        // list's order. The excluded pairs, and any scaled 1-4 pairs, are left out.
        let q = &self.charges;
        let skip = self.skipped();
        let pair_term = |i: usize, j: usize, d: [f64; 3], r2: f64| {
            if skip.contains(i, j) {
                return;
            }
            let (gi, gj) = (marked(i), marked(j));
            if gi && gj {
                return;
            }
            let across = gi != gj;
            let r = r2.sqrt();
            let pair = self.pair(i, j);
            let mut radial;
            if across {
                let s = split.expect("a pair across needs a group");
                let sc = s.soft_core.at(&pair, r, lv);
                e.van_der_waals += sc.energy;
                e.total += sc.energy;
                coupling.energy += sc.energy;
                coupling.gradient[2] += sc.de_dlambda;
                radial = sc.de_dr_over_r;
            } else {
                let (v, dv) = pair.at(r);
                e.van_der_waals += v;
                e.total += v;
                radial = dv / r;
            }
            let qq = q[i] * q[j];
            if qq != 0.0 {
                let (c, dc) = self.ewald.real_pair(qq, r);
                if across {
                    coupling.gradient[1] += c;
                    if le != 0.0 {
                        e.ewald.real += le * c;
                        e.ewald.total += le * c;
                        coupling.energy += le * c;
                        radial += le * dc / r;
                    }
                } else {
                    e.ewald.real += c;
                    e.ewald.total += c;
                    radial += dc / r;
                }
            }
            add_pair(&mut forces, &mut virial, i, j, d, radial);
        };
        match &self.neighbours {
            Some(list) => list.for_each_pair(&cell, at, cutoff, skip, pair_term),
            None => for_each_pair(&cell, at, cutoff, pair_term),
        }

        // The excluded pairs' correction, outside the group.
        for &[i, j] in self.exclusions.pairs() {
            let qq = q[i] * q[j];
            if qq == 0.0 || marked(i) || marked(j) {
                continue;
            }
            let d = cell.minimum_image(sub(at[i], at[j]));
            let r = norm(d);
            let (c, dc) = self.ewald.excluded_pair(qq, r);
            e.ewald.excluded += c;
            e.ewald.total += c;
            add_pair(&mut forces, &mut virial, i, j, d, dc / r);
        }

        // The scaled 1-4 pairs outside the group, wherever they are: `s_v` of the van der Waals,
        // and `s_c` of the screened Coulomb pair with `1 − s_c` of its smooth part taken back out
        // of the reciprocal sum, which leaves `s_c k_e q_i q_j / r` (see the module documentation).
        if let Some(o) = &self.one_four {
            let [sv, sc] = o.scale;
            for &[i, j] in &o.pairs {
                if marked(i) || marked(j) {
                    continue;
                }
                let d = cell.minimum_image(sub(at[i], at[j]));
                let r = norm(d);
                let (v, dv) = self.pair(i, j).at(r);
                e.van_der_waals += sv * v;
                e.total += sv * v;
                let mut radial = sv * dv / r;
                let qq = q[i] * q[j];
                if qq != 0.0 {
                    let (c, dc) = self.ewald.real_pair(qq, r);
                    let (x, dx) = self.ewald.excluded_pair(qq, r);
                    let (real, excluded) = (sc * c, (1.0 - sc) * x);
                    e.ewald.real += real;
                    e.ewald.excluded += excluded;
                    e.ewald.total += real + excluded;
                    radial += (sc * dc + (1.0 - sc) * dx) / r;
                }
                add_pair(&mut forces, &mut virial, i, j, d, radial);
            }
        }

        // Reciprocal space, self and background.
        if q.iter().any(|&x| x != 0.0) {
            let (rest, cross) = self.ewald.reciprocal(
                q,
                split.map(|s| (s.group, le)),
                at,
                forces_wanted.then_some(&mut forces[..]),
                &mut virial,
            );
            e.ewald.reciprocal = rest + le * cross;
            e.ewald.total += rest + le * cross;
            coupling.energy += le * cross;
            coupling.gradient[1] += cross;
        }
        let rest_charges: Vec<f64> = (0..n).map(|i| if marked(i) { 0.0 } else { q[i] }).collect();
        e.ewald.self_energy = self.ewald.self_energy(&rest_charges);
        e.ewald.total += e.ewald.self_energy;
        let q_rest: f64 = rest_charges.iter().sum();
        let q_group: f64 = (0..n).filter(|&i| marked(i)).map(|i| q[i]).sum();
        let unit = self.ewald.background(1.0);
        let (bg_rest, bg_cross) = (unit * q_rest * q_rest, 2.0 * unit * q_rest * q_group);
        e.ewald.background = bg_rest + le * bg_cross;
        e.ewald.total += bg_rest + le * bg_cross;
        coupling.energy += le * bg_cross;
        coupling.gradient[1] += bg_cross;

        // The long-range van der Waals correction.
        let s_rest = self.tail_set(|i| !marked(i));
        let s_group = self.tail_set(marked);
        let t_rest = tail_energy(tail_products(&s_rest, &s_rest), cutoff, volume);
        let t_cross = 2.0 * tail_energy(tail_products(&s_rest, &s_group), cutoff, volume);
        e.dispersion_correction = t_rest + lv * t_cross;
        e.total += t_rest + lv * t_cross;
        coupling.energy += lv * t_cross;
        coupling.gradient[2] += t_cross;

        for (a, row) in virial.iter_mut().enumerate() {
            row[a] += e.ewald.background + e.dispersion_correction;
        }
        e.electrostatic = e.ewald.total;

        // The group's own pairs, in vacuum: on the whole molecules, with no image taken, so
        // that a group of several molecules keeps their placement as given.
        // A 1-4 pair among them is scaled; every other has `[1, 1]`, which changes no bit.
        if let Some(s) = split {
            for (p, &[sv, sc]) in s.own.iter().zip(s.own_scale) {
                let [i, j] = p.atoms;
                let d = sub(whole[i], whole[j]);
                let r = norm(d);
                let (v, dv) = p.at(r);
                e.van_der_waals += sv * v;
                e.total += sv * v;
                let mut radial = sv * dv / r;
                if q[i] != 0.0 && q[j] != 0.0 {
                    let c = sc * coulomb(q[i], q[j], r);
                    e.electrostatic += c;
                    radial -= c / (r * r);
                }
                add_pair(&mut forces, &mut virial, i, j, d, radial);
            }
        }
        e.total += e.electrostatic;
        Core {
            energy: e,
            forces,
            virial,
            coupling,
        }
    }
}

impl Potential for PeriodicForceField {
    fn energy_and_forces(&self, at: &[[f64; 3]], forces: &mut [[f64; 3]]) -> f64 {
        let ev = self.evaluate(at);
        forces.copy_from_slice(&ev.forces);
        ev.energy.total
    }
}

/// A group decoupled from the rest of a [`PeriodicForceField`] by three λs: see the module
/// documentation and [`crate::alchemy`].
#[derive(Clone, Debug, PartialEq)]
pub struct PeriodicDecoupling {
    field: PeriodicForceField,
    group: Vec<bool>,
    own: Vec<Pair>,
    /// Each own pair's `[s_v, s_c]`, the field's 1-4 scales for a 1-4 pair and `[1, 1]` else.
    own_scale: Vec<[f64; 2]>,
    restraint: Option<Boresch>,
    soft_core: SoftCore,
}

impl PeriodicDecoupling {
    /// `field` with the atoms `group` marks decoupled from the rest: no restraint, and
    /// [`SoftCore::default`].
    ///
    /// # Errors
    ///
    /// [`AlchemyError::BondedAcross`] when a bond crosses the partition, and
    /// [`AlchemyError::NothingToDecouple`] when the group is empty or is every atom.
    ///
    /// # Panics
    ///
    /// If `group` is not one per atom.
    pub fn new(
        field: &PeriodicForceField,
        group: &[bool],
    ) -> Result<PeriodicDecoupling, AlchemyError> {
        let n = field.charges.len();
        assert_eq!(group.len(), n, "one mark per atom");
        if !group.iter().any(|&g| g) || group.iter().all(|&g| g) {
            return Err(AlchemyError::NothingToDecouple);
        }
        if let Some(&atoms) = field.bonds.iter().find(|b| group[b[0]] != group[b[1]]) {
            return Err(AlchemyError::BondedAcross { atoms });
        }
        let members: Vec<usize> = (0..n).filter(|&i| group[i]).collect();
        let mut own = Vec::new();
        let mut own_scale = Vec::new();
        let scale = field.one_four();
        let one_four = field.one_four_pairs();
        for (p, &i) in members.iter().enumerate() {
            for &j in &members[p + 1..] {
                if !field.exclusions.contains(i, j) {
                    own.push(field.pair(i, j));
                    let is_one_four = one_four.binary_search(&[i.min(j), i.max(j)]).is_ok();
                    own_scale.push(if is_one_four { scale } else { [1.0, 1.0] });
                }
            }
        }
        Ok(PeriodicDecoupling {
            field: field.clone(),
            group: group.to_vec(),
            own,
            own_scale,
            restraint: None,
            soft_core: SoftCore::default(),
        })
    }

    /// The same decoupling with a Boresch restraint, scaled by `λ_r`. **It reads the positions as
    /// they are given**, not their minimum image: molecular dynamics never wraps, so the
    /// receptor's and the ligand's atoms stay where a continuous trajectory put them.
    ///
    /// # Panics
    ///
    /// If one of its atoms is past the atom count.
    pub fn with_restraint(mut self, restraint: Boresch) -> PeriodicDecoupling {
        let n = self.group.len();
        assert!(
            restraint
                .receptor
                .iter()
                .chain(&restraint.ligand)
                .all(|&a| a < n),
            "a restraint atom past the atom count"
        );
        self.restraint = Some(restraint);
        self
    }

    /// The same decoupling with another soft core.
    pub fn with_soft_core(mut self, soft_core: SoftCore) -> PeriodicDecoupling {
        self.soft_core = soft_core;
        self
    }

    /// The force field.
    pub fn field(&self) -> &PeriodicForceField {
        &self.field
    }

    /// Which atoms are decoupled.
    pub fn group(&self) -> &[bool] {
        &self.group
    }

    /// The group's own non-bonded pairs, computed in vacuum at every λ, combined by the field's
    /// [`Mixing`]; a 1-4 pair among them is scaled by the field's 1-4 scales
    /// ([`PeriodicDecoupling::own_scales`]).
    pub fn own_pairs(&self) -> &[Pair] {
        &self.own
    }

    /// Each of [`PeriodicDecoupling::own_pairs`]' `[s_v, s_c]`: the field's 1-4 scales for a 1-4
    /// pair, and `[1, 1]` for every other.
    pub fn own_scales(&self) -> &[[f64; 2]] {
        &self.own_scale
    }

    /// The energy by term at `at` and `lambda`, the restraint aside.
    ///
    /// # Panics
    ///
    /// If `at` is not one position per atom, or a λ is outside [0, 1].
    pub fn energy(&self, at: &[[f64; 3]], lambda: Lambda) -> PeriodicEnergy {
        self.core(at, lambda, false).energy
    }

    fn core(&self, at: &[[f64; 3]], lambda: Lambda, forces: bool) -> Core {
        assert!(lambda.is_valid(), "every λ must be in [0, 1]: {lambda:?}");
        let split = Split {
            group: &self.group,
            lambda,
            soft_core: self.soft_core,
            own: &self.own,
            own_scale: &self.own_scale,
        };
        self.field.core(at, Some(&split), forces)
    }

    /// What the coupling at any state is made of, at `at`: the group's van der Waals pairs with the
    /// rest inside the cutoff and their distances, found by the same wrapped minimum image as the
    /// evaluation's pair search, so each distance is its bits; the cross Coulomb energy at full
    /// strength — those pairs' real-space terms, the reciprocal cross terms and the background's —
    /// and the long-range correction's cross part.
    fn cross_terms(&self, at: &[[f64; 3]]) -> CrossTerms {
        let f = &self.field;
        let n = f.charges.len();
        assert_eq!(at.len(), n, "one position per atom");
        let cell = f.cell();
        let parameters = f.ewald.parameters();
        let (cutoff, volume) = (parameters.cutoff, cell.volume());
        let rc2 = cutoff * cutoff;
        let l = cell.lengths();
        let half = l.map(|x| 0.5 * x);
        let w: Vec<[f64; 3]> = at.iter().map(|p| cell.wrap(*p)).collect();
        let q = &f.charges;
        let members: Vec<usize> = (0..n).filter(|&i| self.group[i]).collect();
        let mut pairs = Vec::new();
        let mut real = 0.0;
        for &g in &members {
            for j in 0..n {
                if self.group[j] {
                    continue;
                }
                let (a, b) = (g.min(j), g.max(j));
                let mut d = sub(w[a], w[b]);
                for x in 0..3 {
                    if d[x] > half[x] {
                        d[x] -= l[x];
                    } else if d[x] < -half[x] {
                        d[x] += l[x];
                    }
                }
                let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                if r2 >= rc2 || f.exclusions.contains(a, b) {
                    continue;
                }
                let r = r2.sqrt();
                pairs.push((f.pair(a, b), r));
                let qq = q[a] * q[b];
                if qq != 0.0 {
                    real += f.ewald.real_pair(qq, r).0;
                }
            }
        }
        let reciprocal = if q.iter().any(|&x| x != 0.0) {
            let mut unused = [[0.0; 3]; 3];
            f.ewald
                .reciprocal(q, Some((&self.group, 0.0)), at, None, &mut unused)
                .1
        } else {
            0.0
        };
        let rest_charges: Vec<f64> = (0..n)
            .map(|i| if self.group[i] { 0.0 } else { q[i] })
            .collect();
        let q_rest: f64 = rest_charges.iter().sum();
        let q_group: f64 = members.iter().map(|&i| q[i]).sum();
        let background = 2.0 * f.ewald.background(1.0) * q_rest * q_group;
        let s_rest = f.tail_set(|i| !self.group[i]);
        let s_group = f.tail_set(|i| self.group[i]);
        CrossTerms {
            pairs,
            coulomb: real + reciprocal + background,
            tail: 2.0 * tail_energy(tail_products(&s_rest, &s_group), cutoff, volume),
        }
    }

    /// The long-range van der Waals correction between the group and the rest at full coupling,
    /// joules: `2 (2π/V) Σ_{i∈G} Σ_{j∉G} D_ij [x_ij¹²/(9 r_c⁹) − 2 x_ij⁶/(3 r_c³)]`, combined by
    /// the field's [`Mixing`], the part [`PeriodicDecoupling`] scales by `λ_v`. It depends on
    /// nothing that moves, so it is exactly its own contribution to the free energy of decoupling.
    pub fn cross_dispersion_correction(&self) -> f64 {
        let f = &self.field;
        let s_rest = f.tail_set(|i| !self.group[i]);
        let s_group = f.tail_set(|i| self.group[i]);
        2.0 * tail_energy(
            tail_products(&s_rest, &s_group),
            f.ewald.parameters().cutoff,
            f.cell().volume(),
        )
    }
}

/// [`PeriodicDecoupling::cross_terms`]'s result.
struct CrossTerms {
    pairs: Vec<(Pair, f64)>,
    coulomb: f64,
    tail: f64,
}

impl Alchemical for PeriodicDecoupling {
    fn len(&self) -> usize {
        self.group.len()
    }

    fn energy_and_forces(&self, at: &[[f64; 3]], lambda: Lambda, forces: &mut [[f64; 3]]) -> f64 {
        assert_eq!(forces.len(), self.group.len(), "one force per atom");
        let c = self.core(at, lambda, true);
        forces.copy_from_slice(&c.forces);
        let mut energy = c.energy.total;
        if let Some(b) = &self.restraint {
            energy += lambda.restraint * b.add_forces(at, forces, lambda.restraint);
        }
        energy
    }

    /// [`Alchemical::couplings`] at the one state.
    fn coupling(&self, at: &[[f64; 3]], lambda: Lambda) -> Coupling {
        self.couplings(at, &[lambda])[0]
    }

    /// The λ-dependent part at each state from one pass over what they share: the cross pairs and
    /// their distances, the cross Coulomb energy and the correction's cross part
    /// ([`PeriodicDecoupling`]'s `cross_terms`), then each state's soft core over those pairs.
    /// **Not the whole evaluation**: nothing between two atoms of the rest is computed, so a state
    /// costs the group's pairs and one reciprocal sum is shared by all. The energy is
    /// `Σ U_sc(r; λ_v) + λ_e C + λ_v T` and the gradient `[U_B, C, Σ ∂U_sc/∂λ_v + T]`, with `C` the
    /// cross Coulomb energy and `T` the correction's cross part: the λ-dependent part of
    /// [`Alchemical::energy_and_forces`], whose evaluation sums the same terms in another order.
    fn couplings(&self, at: &[[f64; 3]], states: &[Lambda]) -> Vec<Coupling> {
        for l in states {
            assert!(l.is_valid(), "every λ must be in [0, 1]: {l:?}");
        }
        let cross = self.cross_terms(at);
        let restraint = self.restraint.as_ref().map(|b| b.energy(at));
        states
            .iter()
            .map(|l| {
                let mut energy = 0.0;
                let mut soft = 0.0;
                for (pair, r) in &cross.pairs {
                    let sc = self.soft_core.at(pair, *r, l.van_der_waals);
                    energy += sc.energy;
                    soft += sc.de_dlambda;
                }
                energy += l.electrostatics * cross.coulomb;
                energy += l.van_der_waals * cross.tail;
                let mut gradient = [0.0, cross.coulomb, soft + cross.tail];
                if let Some(u) = restraint {
                    gradient[0] = u;
                    energy += l.restraint * u;
                }
                Coupling { energy, gradient }
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ccd::ANGSTROM;
    use crate::water::WaterBox;

    /// The couplings from the cross terms alone are the evaluation's λ-dependent part, summed in
    /// another order: two waters of 27, carrying +0.1 e between them so that the background has a
    /// cross term, at five states, energy and both gradients within 64 ε of the evaluation's
    /// largest part — a few thousand terms each, none larger. And `coupling` is `couplings`' entry
    /// to the bit.
    #[test]
    fn the_cross_terms_are_the_evaluations_coupling() {
        cross_terms_against_the_evaluation(None);
    }

    /// The same through the particle mesh: the couplings take the cross term from one forward
    /// transform of the rest and the group as one complex grid, the evaluation from the same grid
    /// with the potentials' inverse transform as well, and they agree as the classical sum's do.
    /// On a grid of 16 × 8 × 32, so that each axis's own size is in the split of `C(m)` and `C(−m)`.
    #[test]
    fn the_cross_terms_are_the_evaluations_coupling_through_the_mesh() {
        cross_terms_against_the_evaluation(Some(crate::pme::PmeParameters {
            order: 5,
            grid: [16, 8, 32],
        }));
    }

    fn cross_terms_against_the_evaluation(mesh: Option<crate::pme::PmeParameters>) {
        let water = WaterBox::lattice(3, 0x3C0);
        let mut field = water.force_field(4.6 * ANGSTROM, 1e-8);
        if let Some(m) = mesh {
            field = field.with_mesh(m);
        }
        assert_eq!(field.ewald().mesh().is_some(), mesh.is_some());
        let mut q = field.charges().to_vec();
        q[0] += 0.1;
        let field = field.with_charges(q);
        let mut group = vec![false; water.positions().len()];
        group[..6].fill(true);
        let d = PeriodicDecoupling::new(&field, &group).unwrap();
        let at = water.positions();
        let e = field.energy(at);
        let scale = e.ewald.real.abs()
            + e.ewald.reciprocal.abs()
            + e.ewald.self_energy.abs()
            + e.van_der_waals.abs();
        let states = [
            Lambda::COUPLED,
            Lambda::new(0.0, 0.5, 1.0),
            Lambda::new(0.0, 0.0, 1.0),
            Lambda::new(0.0, 0.0, 0.4),
            Lambda::new(0.0, 0.0, 0.0),
        ];
        let fast = d.couplings(at, &states);
        for (l, c) in states.iter().zip(&fast) {
            let slow = d.core(at, *l, false).coupling;
            assert!(
                (c.energy - slow.energy).abs() <= 64.0 * f64::EPSILON * scale,
                "{l:?}"
            );
            for k in 1..3 {
                assert!(
                    (c.gradient[k] - slow.gradient[k]).abs() <= 64.0 * f64::EPSILON * scale,
                    "{l:?} {k}: {} against {}",
                    c.gradient[k],
                    slow.gradient[k]
                );
            }
            let one = d.coupling(at, *l);
            assert_eq!(one.energy.to_bits(), c.energy.to_bits());
            assert_eq!(one.gradient.map(f64::to_bits), c.gradient.map(f64::to_bits));
        }
        assert!(fast[0].gradient[1].abs() > 1e-3 * scale);
        // And the evaluation's reciprocal energy itself, which every comparison above, being a
        // difference or symmetric in the two sets, cannot tell from the rest and the group
        // swapped: the rest's own whole sum plus λ_e times the cross term of three whole sums.
        let reciprocal = |keep: &dyn Fn(usize) -> bool| {
            let q: Vec<f64> = (0..group.len())
                .map(|i| if keep(i) { field.charges()[i] } else { 0.0 })
                .collect();
            field.clone().with_charges(q).energy(at).ewald.reciprocal
        };
        let all = reciprocal(&|_| true);
        let rest = reciprocal(&|i| !group[i]);
        let own = reciprocal(&|i| group[i]);
        for l in [
            Lambda::COUPLED,
            Lambda::new(0.0, 0.5, 1.0),
            Lambda::new(0.0, 0.0, 1.0),
        ] {
            let e = d.core(at, l, false).energy.ewald.reciprocal;
            let want = rest + l.electrostatics * (all - rest - own);
            assert!(
                (e - want).abs() <= 64.0 * f64::EPSILON * scale,
                "{l:?}: {e} against {want}"
            );
        }
        assert!((all - rest).abs() > 1e-3 * scale);
    }
}
