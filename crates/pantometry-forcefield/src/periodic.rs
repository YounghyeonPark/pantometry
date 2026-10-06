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
//! quadrature.
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
//! of cells at least half the cutoff wide, each pair visited once), the correction above, and
//! [`Ewald`] with the same exclusions. It is a [`Potential`], so [`crate::MolecularDynamics`]
//! integrates it unchanged. It never builds the `N²` pair list [`ForceField::new`] keeps.
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
//! all of them, on its whole molecules — the vacuum [`ForceField`]'s for that molecule — and no
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
//! # Not here
//!
//! **No smooth PME**: classical Ewald only, the reference PME will be checked against; its cost
//! is measured above. **No neighbour list** with a skin: the cell list is rebuilt at each
//! evaluation, `O(N)`. **No triclinic box, no barostat**; constraints are the dynamics' —
//! [`crate::water`]'s SETTLE and [`crate::shake`]'s bonds — not the force field's. **No implicit
//! solvent in the box**, and no QEq for a periodic system: the charges are given.

use crate::alchemy::{Alchemical, AlchemyError, Coupling, Lambda, SoftCore};
use crate::boresch::Boresch;
use crate::ccd::Component;
use crate::dynamics::Potential;
use crate::energy::{coulomb, ForceField, Pair, Unsupported};
use crate::ewald::{add_pair, norm, sub, Ewald, EwaldEnergy, EwaldParameters};
use crate::uff::UffType;
use crate::water;
use std::f64::consts::PI;

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

/// Calls `f(i, j, d, r²)` for every pair `i < j` whose minimum-image separation `d = r_i − r_j`
/// is shorter than `cutoff`, once each, in an order fixed by the positions' bits: by `i`, then by
/// the cells about `i`'s in a fixed order, then by `j` within a cell. A cell list: cells at least
/// half the cutoff wide and the 125 about each when the box is at least five such cells across
/// every axis, and every pair when it is fewer.
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
        let mut d = sub(w[i], w[j]);
        for a in 0..3 {
            if d[a] > half[a] {
                d[a] -= l[a];
            } else if d[a] < -half[a] {
                d[a] += l[a];
            }
        }
        let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        if r2 < rc2 {
            f(i, j, d, r2);
        }
    };
    // Cells at least half the cutoff wide, and the 5³ about each. There is no 3³ branch of
    // cutoff-wide cells: three of those across is six half-cells, so it would never be taken.
    let reach = 2usize;
    let cells = l.map(|x| (2.0 * x / cutoff).floor() as usize);
    if cells.iter().any(|&c| c < 2 * reach + 1) {
        for i in 0..n {
            for j in i + 1..n {
                visit(i, j);
            }
        }
        return;
    }
    let index = |p: [f64; 3]| -> [usize; 3] {
        [0, 1, 2].map(|a| ((p[a] / l[a] * cells[a] as f64) as usize).min(cells[a] - 1))
    };
    let flat = |c: [usize; 3]| (c[0] * cells[1] + c[1]) * cells[2] + c[2];
    let total = cells[0] * cells[1] * cells[2];
    let home: Vec<[usize; 3]> = w.iter().map(|p| index(*p)).collect();
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

/// `(2π/V) Σ_i Σ_j D_ij [x_ij¹²/(9 r_c⁹) − 2 x_ij⁶/(3 r_c³)]` for two sets of atoms given by their
/// sums `[Σ √D x⁶, Σ √D x³]`: the long-range correction of every ordered pair with one atom in
/// each, joules.
fn tail_energy(a: [f64; 2], b: [f64; 2], cutoff: f64, volume: f64) -> f64 {
    let rc3 = cutoff * cutoff * cutoff;
    let rc9 = rc3 * rc3 * rc3;
    2.0 * PI / volume * (a[0] * b[0] / (9.0 * rc9) - 2.0 * a[1] * b[1] / (3.0 * rc3))
}

/// A&T's `P_tail` for the same sums, pascals.
fn tail_pressure(a: [f64; 2], b: [f64; 2], cutoff: f64, volume: f64) -> f64 {
    let rc3 = cutoff * cutoff * cutoff;
    let rc9 = rc3 * rc3 * rc3;
    2.0 * PI / (3.0 * volume * volume) * (4.0 * a[0] * b[0] / (3.0 * rc9) - 4.0 * a[1] * b[1] / rc3)
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
    /// Van der Waals, the pairs inside the cutoff (and, in a [`PeriodicDecoupling`], the group's
    /// own in vacuum and its soft-core pairs with the rest).
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
    ewald: Ewald,
    charges: Vec<f64>,
    /// `[x_i, D_i]`, metres and joules.
    vdw: Vec<[f64; 2]>,
    bonds: Vec<[usize; 2]>,
    exclusions: Exclusions,
    /// `[atom, the atom it is reached from]`, breadth first, roots left out.
    tree: Vec<[usize; 2]>,
}

/// What a [`PeriodicDecoupling`] gives the shared evaluation.
struct Split<'a> {
    group: &'a [bool],
    lambda: Lambda,
    soft_core: SoftCore,
    own: &'a [Pair],
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
            ewald: Ewald::new(cell, parameters),
            charges,
            vdw,
            bonds,
            exclusions: Exclusions::from_pairs(n, &excluded),
            tree,
        })
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

    /// The rigid TIP3P waters, atom indices, oxygen first: what [`crate::water::Settle::tip3p`]
    /// is built from.
    pub fn rigid_waters(&self) -> &[[usize; 3]] {
        &self.waters
    }

    /// How many atoms the solute has: the first ones, before the waters.
    pub fn solute_atoms(&self) -> usize {
        self.solute
    }

    /// The same force field in another box, with the same α, cutoff and integer wave vectors
    /// ([`Ewald::with_cell`]): for a derivative with respect to the box.
    ///
    /// # Panics
    ///
    /// If the cutoff is past half the new box's shortest edge.
    pub fn with_cell(&self, cell: PeriodicBox) -> PeriodicForceField {
        let mut f = self.clone();
        f.ewald = self.ewald.with_cell(cell);
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

    /// The long-range van der Waals correction `E_tail`, joules: see the module documentation.
    pub fn dispersion_correction(&self) -> f64 {
        let s = self.tail_sums(|_| true);
        tail_energy(s, s, self.ewald.parameters().cutoff, self.cell().volume())
    }

    /// Allen and Tildesley's long-range correction to the pressure, `P_tail`, pascals: see the
    /// module documentation for how it differs from the virial's `E_tail/V`.
    pub fn pressure_correction(&self) -> f64 {
        let s = self.tail_sums(|_| true);
        tail_pressure(s, s, self.ewald.parameters().cutoff, self.cell().volume())
    }

    /// The van der Waals pair of atoms `i` and `j`, combined as [`Pair::new`] combines them.
    fn pair(&self, i: usize, j: usize) -> Pair {
        let ([xi, di], [xj, dj]) = (self.vdw[i], self.vdw[j]);
        Pair {
            atoms: [i, j],
            distance: (xi * xj).sqrt(),
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

        // The pairs inside the cutoff.
        let q = &self.charges;
        for_each_pair(&cell, at, cutoff, |i, j, d, r2| {
            if self.exclusions.contains(i, j) {
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
        });

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
        let s_rest = self.tail_sums(|i| !marked(i));
        let s_group = self.tail_sums(marked);
        let t_rest = tail_energy(s_rest, s_rest, cutoff, volume);
        let t_cross = 2.0 * tail_energy(s_rest, s_group, cutoff, volume);
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
        if let Some(s) = split {
            for p in s.own {
                let [i, j] = p.atoms;
                let d = sub(whole[i], whole[j]);
                let r = norm(d);
                let (v, dv) = p.at(r);
                e.van_der_waals += v;
                e.total += v;
                let mut radial = dv / r;
                if q[i] != 0.0 && q[j] != 0.0 {
                    let c = coulomb(q[i], q[j], r);
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
        for (p, &i) in members.iter().enumerate() {
            for &j in &members[p + 1..] {
                if !field.exclusions.contains(i, j) {
                    own.push(field.pair(i, j));
                }
            }
        }
        Ok(PeriodicDecoupling {
            field: field.clone(),
            group: group.to_vec(),
            own,
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

    /// The group's own non-bonded pairs, computed in vacuum at every λ.
    pub fn own_pairs(&self) -> &[Pair] {
        &self.own
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
        let s_rest = f.tail_sums(|i| !self.group[i]);
        let s_group = f.tail_sums(|i| self.group[i]);
        CrossTerms {
            pairs,
            coulomb: real + reciprocal + background,
            tail: 2.0 * tail_energy(s_rest, s_group, cutoff, volume),
        }
    }

    /// The long-range van der Waals correction between the group and the rest at full coupling,
    /// joules: `2 (2π/V) Σ_{i∈G} Σ_{j∉G} D_ij [x_ij¹²/(9 r_c⁹) − 2 x_ij⁶/(3 r_c³)]`, the part
    /// [`PeriodicDecoupling`] scales by `λ_v`. It depends on nothing that moves, so it is exactly
    /// its own contribution to the free energy of decoupling.
    pub fn cross_dispersion_correction(&self) -> f64 {
        let f = &self.field;
        let s_rest = f.tail_sums(|i| !self.group[i]);
        let s_group = f.tail_sums(|i| self.group[i]);
        2.0 * tail_energy(
            s_rest,
            s_group,
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
        let water = WaterBox::lattice(3, 0x3C0);
        let field = water.force_field(4.6 * ANGSTROM, 1e-8);
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
    }
}
