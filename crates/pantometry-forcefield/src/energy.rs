//! The Universal Force Field's energy, all six of its terms, with their analytic forces.
//!
//! Bond stretch, van der Waals and electrostatics here, and angle bend, torsion and inversion in
//! [`crate::angular`], from Rappé et al., *J. Am. Chem. Soc.* **114**, 10024 (1992) — page and
//! equation numbers below are that paper's. [`ForceField`] holds every term of one molecule and
//! gives the energy by term and the force on every atom. Nothing here minimises; that is a later
//! step, and so is the charge model, without which the electrostatic term is zero by default.
//!
//! # Bond stretch (eq 1a, p. 10025)
//!
//! `E = ½ k_IJ (r − r_IJ)²`, harmonic. The natural length (eq 2, p. 10025) is
//!
//! `r_IJ = r_I + r_J + r_BO − r_EN`
//!
//! with the bond-order correction (eq 3, p. 10026) `r_BO = −λ (r_I + r_J) ln n`, λ = 0.1332
//! (p. 10026), and the electronegativity correction (eq 4, p. 10027)
//! `r_EN = r_I r_J (√χ_I − √χ_J)² / (χ_I r_I + χ_J r_J)`. The force constant (eq 6, p. 10027) is
//! `k_IJ = 664.12 Z_I* Z_J* / r_IJ³` kcal mol⁻¹ Å⁻², with Z* from Table I.
//!
//! **A decision: `r_EN` is subtracted.** The paper prints eq 2 as `+ r_EN`. Open Babel and RDKit
//! both subtract it, and the paper's own minimised structures support subtracting: dimethyl ether's
//! C–O natural length is 1.394 Å with the minus and 1.436 Å with the plus, and the paper's
//! minimised value is 1.410 Å (Fig. 6, p. 10032); trimethylamine's C–N is 1.451 Å against 1.463 Å,
//! minimised 1.471 Å (Fig. 5). With the minus both bonds stretch on minimisation (+0.016 and
//! +0.020 Å), which angle strain in a crowded centre would do; with the plus one stretches and one
//! compresses. **Evidence the other way, found while writing this:** p. 10026 says the amide C–N
//! bond order 1.41 was chosen to reproduce N-methylformamide's 1.366 Å "from the C_R and N_R
//! single bond radii". [`uff::assign`] types an amide's carbon `C_R` and its nitrogen `N_R` for
//! that reason, so this is the bond this crate computes: with `C_R` 0.729 and `N_R` 0.699 Å and
//! n = 1.41, the minus gives 1.35684 Å and the printed plus 1.36845 Å, against 1.366. The order
//! that would reproduce 1.366 exactly is 1.428 with the plus and 1.344 with the minus. **So the
//! paper's own sentence favours the printed sign.** The minus stays, as the two-implementation
//! consensus and the two minimised structures say; the question is open, and the minimisation step
//! is what settles it. See [`natural_length`].
//!
//! The bond order `n` is read from the dictionary's bond ([`bond_order`]): `SING` 1, `DOUB` 2,
//! `TRIP` 3; **an aromatic bond 1.5** whatever its Kekulé order (p. 10027: "Intra-ring bonds of
//! aromatic rings are assigned bond orders based on the number of π electrons, resulting in a bond
//! order of 1.5 for normal aromatic rings"); and **an amide C–N 1.41** (p. 10026) — the single bond
//! between a carbon carrying a C=O and a non-aromatic nitrogen that [`uff::assign`] typed `N_R` by
//! its amide rule.
//!
//! Two recorded choices at the edge of those rules, neither dictated by the paper:
//!
//! - **An N-acyl aromatic nitrogen is not an amide.** In N-acetylimidazole the ring nitrogen's lone
//!   pair is part of the ring's aromatic sextet, so it is not available to the carbonyl the way an
//!   amide nitrogen's is; the exocyclic N–C(=O) bond is written single and not aromatic, and gets
//!   order 1, not 1.41, and its carbon stays `C_2`.
//! - **A ring the dictionary flags aromatic gets 1.5 throughout, carbonyl or not.** If the CCD
//!   flags a pyridone or pyrimidinone ring aromatic, its ring N–C(=O) bonds get 1.5 rather than
//!   1.41. The paper's "normal aromatic rings" may not have been meant to cover these, and the
//!   dictionary's flag is taken as stated rather than second-guessed.
//!
//! # Van der Waals (eq 20, p. 10029)
//!
//! Lennard-Jones 12-6: `E = D_IJ { −2 (x_IJ/x)⁶ + (x_IJ/x)¹² }`, minimum `−D_IJ` at `x = x_IJ`,
//! with "standard geometric combination rules" (p. 10029): `x_IJ = √(x_I x_J)` (eq 21b) and
//! `D_IJ = √(D_I D_J)` (eq 22). **No cutoff**: every included pair is summed, which is exact and
//! costs O(N²), and N is a drug molecule's.
//!
//! # Electrostatics (eq 43, p. 10031)
//!
//! `E = 332.0637 Q_i Q_j / (ε R_ij)` kcal/mol, with Q in elementary charges, R in Å and ε = 1.
//! **The charges are an input, and default to zero.** The paper (p. 10031) obtained the valence
//! parameters without partial charges and takes charges, when it uses them, from charge
//! equilibration; that charge model is step 2. Until it arrives a molecule built here has no
//! electrostatic energy unless a caller supplies charges with [`ForceField::with_charges`], and
//! a dictionary's *formal* charges are deliberately not used, because they are not partial
//! charges.
//!
//! # Exclusions (p. 10031, "G. Nonbonded Exclusions")
//!
//! Van der Waals and electrostatics both exclude 1-2 and 1-3 pairs — atoms bonded, or bonded to a
//! common atom. Every other pair, 1-4 included, interacts at full strength. A pair joined by both a
//! short and a long path (in a small ring) is excluded if either path is short.
//!
//! # Units
//!
//! SI throughout: metres, joules per molecule, newtons. The paper's constants are converted once:
//! 1 Å = 10⁻¹⁰ m, and 1 kcal/mol = 4184 J / N_A ([`KCAL_PER_MOL`]). So `k_IJ` is in J m⁻² and the
//! Coulomb prefactor in J m per e². The two equation-level corrections, [`bond_order_correction`]
//! and [`electronegativity_correction`], are homogeneous of degree one in the radii and return
//! whatever length unit they are given, so they can be checked against the paper's ångström
//! directly.

use crate::angular::{inversion_parameters, torsion_parameters, Bend, Inversion, Torsion};
use crate::ccd::{BondOrder, Component, Element, ANGSTROM};
use crate::uff::{gmp_electronegativity, Hybridisation, UffType, KCAL_PER_MOL};
use std::fmt;

#[cfg(doc)]
use crate::uff;

/// λ in eq 3 (p. 10026), the bond-order correction's proportionality constant.
pub const LAMBDA: f64 = 0.1332;

/// The prefactor of eq 6 (p. 10027): `k_IJ = 664.12 Z_I* Z_J* / r_IJ³`, in kcal mol⁻¹ Å e⁻².
pub const BOND_FORCE_PREFACTOR: f64 = 664.12;

/// The Coulomb prefactor of eq 43 (p. 10031), in kcal mol⁻¹ Å e⁻²: `e² / (4π ε₀)` in those
/// units, which CODATA's ε₀ reproduces to the seven figures printed.
pub const COULOMB_KCAL: f64 = 332.0637;

/// The bond order of an aromatic ring bond (p. 10027).
pub const AROMATIC_BOND_ORDER: f64 = 1.5;

/// The bond order of an amide C–N bond (p. 10026).
pub const AMIDE_BOND_ORDER: f64 = 1.41;

/// A molecule this module will not compute an energy for, and why.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Unsupported {
    /// A sulfur whose bond orders sum past two, so not the divalent `S_3+2` (or the `S_2` and
    /// `S_R` with sum two). Sum four is a sulfoxide and six a sulfone, which UFF types `S_3+4`
    /// and `S_3+6` and this crate has no parameters for; three is a sulfonium or similar, which
    /// none of UFF's sulfur types as held here describes. [`uff::assign`] types each of them as if
    /// it were divalent, whose radius is wrong for it. Refused rather than computed wrong.
    HypervalentSulfur {
        /// The atom's name in the dictionary.
        atom: String,
        /// The sum of its bond orders, as written.
        valence: u32,
    },
}

impl fmt::Display for Unsupported {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Unsupported::HypervalentSulfur { atom, valence } => {
                let what = match valence {
                    4 => "a sulfoxide's, UFF type S_3+4",
                    6 => "a sulfone's, UFF type S_3+6",
                    _ => "neither divalent (S_3+2) nor a sulfoxide's or sulfone's",
                };
                write!(
                    f,
                    "sulfur {atom} has bond-order sum {valence} ({what}), which this force field \
                     has no parameters for"
                )
            }
        }
    }
}

impl std::error::Error for Unsupported {}

/// Eq 3: `r_BO = −λ (r_I + r_J) ln n`, in the unit of `r_i` and `r_j`. Exactly zero for a single
/// bond and negative — a shortening — above it.
pub fn bond_order_correction(r_i: f64, r_j: f64, n: f64) -> f64 {
    -LAMBDA * (r_i + r_j) * n.ln()
}

/// Eq 4: `r_EN = r_I r_J (√χ_I − √χ_J)² / (χ_I r_I + χ_J r_J)`, in the unit of `r_i` and `r_j`
/// (the χ unit cancels). Exactly zero between equal electronegativities, and positive otherwise.
pub fn electronegativity_correction(r_i: f64, chi_i: f64, r_j: f64, chi_j: f64) -> f64 {
    let d = chi_i.sqrt() - chi_j.sqrt();
    r_i * r_j * d * d / (chi_i * r_i + chi_j * r_j)
}

/// Eq 2, the natural length `r_IJ = r_I + r_J + r_BO − r_EN` of a bond of order `n` between types
/// `a` and `b`, in metres.
///
/// **The minus on `r_EN` departs from the paper as printed**; see the module documentation for
/// why, and for the one sentence in the paper that argues the other way.
pub fn natural_length(a: UffType, b: UffType, n: f64) -> f64 {
    let (ri, rj) = (a.table_i().r1, b.table_i().r1);
    let (ci, cj) = (
        gmp_electronegativity(a.element()),
        gmp_electronegativity(b.element()),
    );
    let angstrom =
        ri + rj + bond_order_correction(ri, rj, n) - electronegativity_correction(ri, ci, rj, cj);
    angstrom * ANGSTROM
}

/// Eq 6, the force constant `k_IJ = 664.12 Z_I* Z_J* / r_IJ³` for a bond of natural length
/// `natural` (metres) between types `a` and `b`, in J m⁻².
pub fn bond_force_constant(a: UffType, b: UffType, natural: f64) -> f64 {
    let r = natural / ANGSTROM;
    let kcal_per_angstrom2 = BOND_FORCE_PREFACTOR * a.table_i().z1 * b.table_i().z1 / (r * r * r);
    kcal_per_angstrom2 * KCAL_PER_MOL / (ANGSTROM * ANGSTROM)
}

/// Eq 43, the Coulomb energy of charges `q_i` and `q_j` (elementary charges) at `r` metres, ε = 1,
/// in joules per molecule.
pub fn coulomb(q_i: f64, q_j: f64, r: f64) -> f64 {
    COULOMB_KCAL * KCAL_PER_MOL * ANGSTROM * q_i * q_j / r
}

/// Whether atom `i` is aromatic: its own flag, or any of its bonds' — the reading
/// [`uff::assign`] makes.
fn aromatic(component: &Component, i: usize) -> bool {
    component.atoms()[i].aromatic || component.neighbours(i).any(|(_, b)| b.aromatic)
}

/// Whether carbon `c` carries a C=O.
fn carbonyl_carbon(component: &Component, c: usize) -> bool {
    let atoms = component.atoms();
    atoms[c].element == Element::C
        && component
            .neighbours(c)
            .any(|(o, b)| atoms[o].element == Element::O && b.order == BondOrder::Double)
}

/// The UFF bond order of `component`'s bond number `bond`, given the atoms' `types`: aromatic 1.5,
/// amide C–N 1.41, otherwise the dictionary's 1, 2 or 3. See the module documentation.
///
/// # Panics
///
/// If `bond` is not a bond of `component`, or `types` is shorter than its atom list.
pub fn bond_order(component: &Component, types: &[UffType], bond: usize) -> f64 {
    let b = &component.bonds()[bond];
    if b.aromatic {
        return AROMATIC_BOND_ORDER;
    }
    if b.order == BondOrder::Single {
        let [i, j] = b.atoms;
        for (n, c) in [(i, j), (j, i)] {
            if types[n] == UffType::NR && !aromatic(component, n) && carbonyl_carbon(component, c) {
                return AMIDE_BOND_ORDER;
            }
        }
    }
    f64::from(b.order.count())
}

/// One harmonic bond-stretch term.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Stretch {
    /// The two atoms, as indices into the component.
    pub atoms: [usize; 2],
    /// The UFF bond order `n`.
    pub order: f64,
    /// The natural length `r_IJ`, metres.
    pub natural_length: f64,
    /// The force constant `k_IJ`, J m⁻².
    pub force_constant: f64,
}

impl Stretch {
    /// The term between types `a` and `b` at bond order `n`.
    pub fn new(atoms: [usize; 2], a: UffType, b: UffType, n: f64) -> Stretch {
        let natural_length = natural_length(a, b, n);
        Stretch {
            atoms,
            order: n,
            natural_length,
            force_constant: bond_force_constant(a, b, natural_length),
        }
    }

    /// `½ k (r − r_IJ)²` at length `r` metres, in joules per molecule.
    pub fn energy(&self, r: f64) -> f64 {
        self.at(r).0
    }

    /// The energy and `dE/dr` at length `r`: the one place the formula is written.
    fn at(&self, r: f64) -> (f64, f64) {
        let d = r - self.natural_length;
        (0.5 * self.force_constant * d * d, self.force_constant * d)
    }
}

/// One van der Waals pair, with its geometrically combined parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pair {
    /// The two atoms, as indices into the component.
    pub atoms: [usize; 2],
    /// `x_IJ = √(x_I x_J)`, the distance of the minimum, metres.
    pub distance: f64,
    /// `D_IJ = √(D_I D_J)`, the depth of the minimum, joules per molecule.
    pub well: f64,
}

impl Pair {
    /// The pair between types `a` and `b`, combined by eqs 21b and 22.
    pub fn new(atoms: [usize; 2], a: UffType, b: UffType) -> Pair {
        let (pa, pb) = (a.parameters(), b.parameters());
        Pair {
            atoms,
            distance: (pa.vdw_distance * pb.vdw_distance).sqrt(),
            well: (pa.vdw_energy * pb.vdw_energy).sqrt(),
        }
    }

    /// `D_IJ { −2 (x_IJ/x)⁶ + (x_IJ/x)¹² }` at separation `x` metres, in joules per molecule.
    pub fn energy(&self, x: f64) -> f64 {
        self.at(x).0
    }

    /// The energy and `dE/dx` at separation `x`: the one place the formula is written.
    /// `dE/dx = 12 D (s⁶ − s¹²) / x` with `s = x_IJ / x`.
    fn at(&self, x: f64) -> (f64, f64) {
        let s = self.distance / x;
        let s6 = (s * s) * (s * s) * (s * s);
        let s12 = s6 * s6;
        (
            self.well * (s12 - 2.0 * s6),
            12.0 * self.well * (s6 - s12) / x,
        )
    }
}

/// An energy, by term and in total, in joules per molecule.
///
/// `total` is accumulated separately from the six terms, term by term as each is evaluated, so
/// that the parts summing to the whole is a check on the bookkeeping rather than a definition.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Energy {
    /// Bond stretch.
    pub bond: f64,
    /// Angle bend.
    pub angle: f64,
    /// Torsion.
    pub torsion: f64,
    /// Inversion.
    pub inversion: f64,
    /// Van der Waals.
    pub van_der_waals: f64,
    /// Electrostatics.
    pub electrostatic: f64,
    /// Everything above.
    pub total: f64,
}

impl Energy {
    /// Every field `NaN`: the energy of a molecule the force field refused.
    pub const NAN: Energy = Energy {
        bond: f64::NAN,
        angle: f64::NAN,
        torsion: f64::NAN,
        inversion: f64::NAN,
        van_der_waals: f64::NAN,
        electrostatic: f64::NAN,
        total: f64::NAN,
    };
}

/// An energy and the force on every atom, `−∂E/∂r`, in newtons.
#[derive(Clone, Debug, PartialEq)]
pub struct Evaluation {
    /// The energy.
    pub energy: Energy,
    /// The force on each atom, in the component's atom order.
    pub forces: Vec<[f64; 3]>,
}

/// All six of UFF's terms for one molecule, with every bend, torsion, inversion and pair resolved.
#[derive(Clone, Debug, PartialEq)]
pub struct ForceField {
    stretches: Vec<Stretch>,
    bends: Vec<Bend>,
    torsions: Vec<Torsion>,
    inversions: Vec<Inversion>,
    pairs: Vec<Pair>,
    charges: Vec<f64>,
}

impl ForceField {
    /// The terms of `component` typed as `types`: one [`Stretch`] per bond; one [`Bend`] per pair
    /// of bonds sharing an atom; a [`Torsion`] for every I–J–K–L about each bond the rules of
    /// [`torsion_parameters`] give a barrier, sharing it; three [`Inversion`]s for each
    /// three-coordinate atom [`inversion_parameters`] gives one; and one [`Pair`] for every pair
    /// of atoms that is neither 1-2 nor 1-3. Charges start at zero. Every list is in the
    /// component's atom and bond order, so the result is the same on every run.
    ///
    /// # Errors
    ///
    /// [`Unsupported::HypervalentSulfur`] for a sulfoxide or sulfone sulfur.
    ///
    /// # Panics
    ///
    /// If `types` is not one per atom.
    pub fn new(component: &Component, types: &[UffType]) -> Result<ForceField, Unsupported> {
        let atoms = component.atoms();
        let n = atoms.len();
        assert_eq!(types.len(), n, "one UFF type per atom");
        for (i, atom) in atoms.iter().enumerate() {
            if atom.element == Element::S {
                let valence: u32 = component.neighbours(i).map(|(_, b)| b.order.count()).sum();
                if valence > 2 {
                    return Err(Unsupported::HypervalentSulfur {
                        atom: atom.name.clone(),
                        valence,
                    });
                }
            }
        }
        let stretches: Vec<Stretch> = component
            .bonds()
            .iter()
            .enumerate()
            .map(|(k, b)| {
                let [i, j] = b.atoms;
                Stretch::new(b.atoms, types[i], types[j], bond_order(component, types, k))
            })
            .collect();
        // Each atom's neighbours with the stretch joining them, in bond order.
        let mut adjacent: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
        for (k, s) in stretches.iter().enumerate() {
            let [i, j] = s.atoms;
            adjacent[i].push((j, k));
            adjacent[j].push((i, k));
        }
        let length = |k: usize| stretches[k].natural_length;

        let mut bends = Vec::new();
        for (j, around) in adjacent.iter().enumerate() {
            for (p, &(i, ki)) in around.iter().enumerate() {
                for &(k, kk) in &around[p + 1..] {
                    bends.push(Bend::new(
                        [i, j, k],
                        [types[i], types[j], types[k]],
                        length(ki),
                        length(kk),
                    ));
                }
            }
        }

        let sp2 = |t: UffType| t.hybridisation() == Some(Hybridisation::Sp2);
        let has_sp2_neighbour =
            |a: usize, not: usize| adjacent[a].iter().any(|&(o, _)| o != not && sp2(types[o]));
        let mut torsions = Vec::new();
        for s in &stretches {
            let [j, k] = s.atoms;
            let Some(parameters) = torsion_parameters(
                types[j],
                types[k],
                s.order,
                has_sp2_neighbour(j, k),
                has_sp2_neighbour(k, j),
            ) else {
                continue;
            };
            let mut paths = Vec::new();
            for &(i, _) in &adjacent[j] {
                for &(l, _) in &adjacent[k] {
                    if i != k && l != j && i != l {
                        paths.push([i, j, k, l]);
                    }
                }
            }
            let m = paths.len();
            torsions.extend(paths.into_iter().map(|atoms| Torsion {
                atoms,
                parameters,
                torsions_about_bond: m,
            }));
        }

        let mut inversions = Vec::new();
        for (i, around) in adjacent.iter().enumerate() {
            if let [(a, _), (b, _), (c, _)] = around[..] {
                let to_o2 = around.iter().any(|&(o, _)| types[o] == UffType::O2);
                if let Some((force_constant, coefficients)) = inversion_parameters(types[i], to_o2)
                {
                    for [j, k, l] in [[b, c, a], [c, a, b], [a, b, c]] {
                        inversions.push(Inversion {
                            atoms: [i, j, k, l],
                            force_constant,
                            coefficients,
                        });
                    }
                }
            }
        }

        let mut excluded = vec![false; n * n];
        for i in 0..n {
            for (j, _) in component.neighbours(i) {
                excluded[i * n + j] = true;
                for (k, _) in component.neighbours(j) {
                    excluded[i * n + k] = true;
                }
            }
        }
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                if !excluded[i * n + j] && !excluded[j * n + i] {
                    pairs.push(Pair::new([i, j], types[i], types[j]));
                }
            }
        }
        Ok(ForceField {
            stretches,
            bends,
            torsions,
            inversions,
            pairs,
            charges: vec![0.0; n],
        })
    }

    /// The same terms with per-atom partial charges `charges`, in elementary charges.
    ///
    /// # Panics
    ///
    /// If `charges` is not one per atom.
    pub fn with_charges(mut self, charges: Vec<f64>) -> ForceField {
        assert_eq!(charges.len(), self.charges.len(), "one charge per atom");
        self.charges = charges;
        self
    }

    /// One bond-stretch term per bond, in the component's bond order.
    pub fn stretches(&self) -> &[Stretch] {
        &self.stretches
    }

    /// One angle-bend term per pair of bonds sharing an atom.
    pub fn bends(&self) -> &[Bend] {
        &self.bends
    }

    /// Every torsion term, grouped by central bond in the component's bond order.
    pub fn torsions(&self) -> &[Torsion] {
        &self.torsions
    }

    /// Every inversion term, three per centre that has one.
    pub fn inversions(&self) -> &[Inversion] {
        &self.inversions
    }

    /// Every non-bonded pair — neither 1-2 nor 1-3 — with `atoms[0] < atoms[1]`.
    pub fn pairs(&self) -> &[Pair] {
        &self.pairs
    }

    /// The partial charges, in elementary charges.
    pub fn charges(&self) -> &[f64] {
        &self.charges
    }

    /// The energy at positions `at` (metres), without the forces.
    ///
    /// # Panics
    ///
    /// If `at` is not one position per atom.
    pub fn energy(&self, at: &[[f64; 3]]) -> Energy {
        self.evaluate(at).energy
    }

    /// The energy at positions `at` (metres), and the analytic force on every atom.
    ///
    /// Each term's force is computed once and added to one atom and subtracted from the other, so
    /// the net force is zero to the rounding of the sums.
    ///
    /// # Panics
    ///
    /// If `at` is not one position per atom.
    pub fn evaluate(&self, at: &[[f64; 3]]) -> Evaluation {
        assert_eq!(at.len(), self.charges.len(), "one position per atom");
        let mut e = Energy::default();
        let mut forces = vec![[0.0; 3]; at.len()];
        // The force on atom i of a term E(r), r = |r_i − r_j|, is −(dE/dr) (r_i − r_j)/r, and
        // atom j's is its negative.
        fn push(forces: &mut [[f64; 3]], [i, j]: [usize; 2], d: [f64; 3], de_dr: f64, r: f64) {
            let g = -de_dr / r;
            for k in 0..3 {
                forces[i][k] += g * d[k];
                forces[j][k] -= g * d[k];
            }
        }
        for s in &self.stretches {
            let (d, r) = separation(at, s.atoms);
            let (energy, de_dr) = s.at(r);
            e.bond += energy;
            e.total += energy;
            push(&mut forces, s.atoms, d, de_dr, r);
        }
        for b in &self.bends {
            let energy = b.accumulate(at, &mut forces);
            e.angle += energy;
            e.total += energy;
        }
        for t in &self.torsions {
            let energy = t.accumulate(at, &mut forces);
            e.torsion += energy;
            e.total += energy;
        }
        for v in &self.inversions {
            let energy = v.accumulate(at, &mut forces);
            e.inversion += energy;
            e.total += energy;
        }
        for p in &self.pairs {
            let (d, r) = separation(at, p.atoms);
            let (vdw, mut de_dr) = p.at(r);
            e.van_der_waals += vdw;
            e.total += vdw;
            let [i, j] = p.atoms;
            let (qi, qj) = (self.charges[i], self.charges[j]);
            if qi != 0.0 && qj != 0.0 {
                let q = coulomb(qi, qj, r);
                de_dr -= q / r;
                e.electrostatic += q;
                e.total += q;
            }
            push(&mut forces, p.atoms, d, de_dr, r);
        }
        Evaluation { energy: e, forces }
    }
}

/// `r_i − r_j` and its length.
fn separation(at: &[[f64; 3]], [i, j]: [usize; 2]) -> ([f64; 3], f64) {
    let d = [
        at[i][0] - at[j][0],
        at[i][1] - at[j][1],
        at[i][2] - at[j][2],
    ];
    (d, (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
}
