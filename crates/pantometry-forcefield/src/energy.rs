//! The Universal Force Field's energy, all six of its terms, with their analytic forces.
//!
//! Bond stretch, van der Waals and electrostatics here, and angle bend, torsion and inversion in
//! [`crate::angular`], from Rappé et al., *J. Am. Chem. Soc.* **114**, 10024 (1992) — page and
//! equation numbers below are that paper's. [`ForceField`] holds every term of one molecule and
//! gives the energy by term and the force on every atom. Nothing here minimises — see
//! [`crate::minimise`] — and the charges are an input, zero by default; [`crate::qeq`] computes the
//! ones UFF prescribes.
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
//! consensus and the two minimised structures say. See [`natural_length`].
//!
//! **Measured once there was a minimiser** (`tests/the_papers_structures.rs`,
//! `tests/the_papers_barriers.rs`): the paper's own *minimised* numbers, which are what its
//! sentence about N-methylformamide's 1.366 Å is really about, agree with the minus. Of 24
//! heteronuclear bonds in its Figures 4–7, the minus is nearer the paper than the plus for all 24,
//! and nearer than no correction at all ([`ElectronegativitySign::Dropped`]) for all 24 — dimethyl
//! ether's C–O 1.4096 Å against the printed 1.410 (1.4497 with the plus), N-methylformamide's
//! C–N 1.3660 against 1.365 (1.3766), and a C=O 0.04 Å long with the plus. Of Table II's nine
//! simple rotors, the minus reproduces eight to the printed figures, and the plus misses five of
//! them. The natural length is not what the paper reports; the relaxed one is. **The code is not
//! changed on this evidence**: the sign is the maintainer's decision, and [`Variant`] exists so
//! that the comparison can be rerun.
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
//! equilibration. A molecule built here has no electrostatic energy unless a caller supplies
//! charges — its own with [`ForceField::with_charges`], or QEq's at a stated geometry with
//! [`ForceField::with_qeq_charges`], fixed there (see [`crate::qeq`]) — and a dictionary's *formal*
//! charges are deliberately not used as partial charges; QEq takes only their sum, as the total.
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
use crate::solvation::GeneralizedBorn;
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
    natural_length_with(a, b, n, ElectronegativitySign::Subtracted)
}

/// Which sign eq 2's `r_EN` takes. [`ElectronegativitySign::Subtracted`] is what the crate
/// computes; the other exists so that the question the module documentation leaves open can be
/// measured.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum ElectronegativitySign {
    /// `r_IJ = r_I + r_J + r_BO − r_EN`, as Open Babel and RDKit compute it.
    #[default]
    Subtracted,
    /// `r_IJ = r_I + r_J + r_BO + r_EN`, as eq 2 is printed. **Not the default and not a claim.**
    AddedAsPrinted,
    /// `r_IJ = r_I + r_J + r_BO`, no electronegativity correction at all: the baseline that says
    /// whether the comparison of the two signs measures the sign or only the correction's size.
    Dropped,
}

/// A reading of the paper other than this crate's, for measuring the evidence for this crate's.
/// `Variant::default()` is what [`ForceField::new`] builds. One question is left in it, the sign
/// of `r_EN`, kept as the evidence the module documentation cites: rerunning
/// `tests/the_papers_structures.rs` under it is how the sign was settled. (A second, the minimum
/// of the group-6 sp³–sp² torsion, was removed when the typing that made it matter was changed;
/// see [`crate::angular`].)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Variant {
    /// The sign of `r_EN` in eq 2.
    pub electronegativity: ElectronegativitySign,
}

/// [`natural_length`] with `r_EN` taken with either sign.
pub fn natural_length_with(a: UffType, b: UffType, n: f64, sign: ElectronegativitySign) -> f64 {
    let (ri, rj) = (a.table_i().r1, b.table_i().r1);
    let (ci, cj) = (
        gmp_electronegativity(a.element()),
        gmp_electronegativity(b.element()),
    );
    let r_en = electronegativity_correction(ri, ci, rj, cj);
    let angstrom = match sign {
        ElectronegativitySign::Subtracted => ri + rj + bond_order_correction(ri, rj, n) - r_en,
        ElectronegativitySign::AddedAsPrinted => ri + rj + bond_order_correction(ri, rj, n) + r_en,
        ElectronegativitySign::Dropped => ri + rj + bond_order_correction(ri, rj, n),
    };
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
        Stretch::with_sign(atoms, a, b, n, ElectronegativitySign::Subtracted)
    }

    /// [`Stretch::new`] with `r_EN` taken with either sign.
    pub fn with_sign(
        atoms: [usize; 2],
        a: UffType,
        b: UffType,
        n: f64,
        sign: ElectronegativitySign,
    ) -> Stretch {
        let natural_length = natural_length_with(a, b, n, sign);
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
    pub(crate) fn at(&self, x: f64) -> (f64, f64) {
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
/// `total` is accumulated separately from the terms, term by term as each is evaluated, so
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
    /// The generalized Born electrostatic solvation energy, ΔG_GB of [`crate::solvation`]: zero
    /// unless the force field was given a solvent ([`ForceField::with_solvation`]). Not one of
    /// UFF's six terms.
    pub solvation: f64,
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
        solvation: f64::NAN,
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
    elements: Vec<Element>,
    total_charge: i32,
    solvation: Option<GeneralizedBorn>,
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
        ForceField::with_variant(component, types, Variant::default())
    }

    /// [`ForceField::new`] under another reading of the two questions [`Variant`] names — for
    /// measuring them. `with_variant(c, t, Variant::default())` is `new(c, t)`.
    ///
    /// # Errors
    ///
    /// As [`ForceField::new`].
    ///
    /// # Panics
    ///
    /// As [`ForceField::new`].
    pub fn with_variant(
        component: &Component,
        types: &[UffType],
        variant: Variant,
    ) -> Result<ForceField, Unsupported> {
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
                Stretch::with_sign(
                    b.atoms,
                    types[i],
                    types[j],
                    bond_order(component, types, k),
                    variant.electronegativity,
                )
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
            elements: atoms.iter().map(|a| a.element).collect(),
            total_charge: atoms.iter().map(|a| a.charge).sum(),
            solvation: None,
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

    /// The same terms with **charge-equilibration (QEq) charges** computed at positions `at`
    /// (metres) — Rappé and Goddard's, the charges UFF prescribes — for the total charge the
    /// dictionary states atom by atom. See [`crate::qeq`].
    ///
    /// **The charges are fixed at `at`.** They are computed once, here, and the force field then
    /// treats them as constants: its electrostatic force is the derivative at fixed charges, which
    /// is how the QEq paper uses them ("solved once for a given structure"), and not the derivative
    /// of an energy whose charges move with the atoms. Call this again at a geometry that has moved
    /// far. [`ForceField::new`]'s zero charges stay the default, because UFF's valence parameters
    /// were obtained without charges (p. 10031).
    ///
    /// # Errors
    ///
    /// Any [`QeqError`](crate::qeq::QeqError).
    ///
    /// # Panics
    ///
    /// If `at` is not one position per atom, or a position is not finite.
    pub fn with_qeq_charges(self, at: &[[f64; 3]]) -> Result<ForceField, crate::qeq::QeqError> {
        assert_eq!(at.len(), self.charges.len(), "one position per atom");
        let charges = crate::qeq::Qeq::default().equilibrate(
            &self.elements,
            at,
            f64::from(self.total_charge),
        )?;
        Ok(self.with_charges(charges.charges))
    }

    /// The same terms **in implicit water**: [`GeneralizedBorn::new`] — OBC II, ε = 80, no salt —
    /// for its atoms, added to the energy as [`Energy::solvation`] with its analytic force. See
    /// [`crate::solvation`]. The solvation energy is computed from the force field's charges, so
    /// with the default zero charges it is zero: give it charges ([`ForceField::with_charges`] or
    /// [`ForceField::with_qeq_charges`]) as well.
    pub fn with_generalized_born(self) -> ForceField {
        let gb = GeneralizedBorn::new(&self.elements);
        self.with_solvation(gb)
    }

    /// The same terms with the solvation model `solvation` added to the energy and the force.
    /// The default, [`ForceField::new`], is vacuum: no solvation term.
    ///
    /// # Panics
    ///
    /// If `solvation` is not one atom per atom of the force field.
    pub fn with_solvation(mut self, solvation: GeneralizedBorn) -> ForceField {
        assert_eq!(solvation.len(), self.charges.len(), "one GB atom per atom");
        self.solvation = Some(solvation);
        self
    }

    /// The solvation model, if one was given; `None` is vacuum.
    pub fn solvation(&self) -> Option<&GeneralizedBorn> {
        self.solvation.as_ref()
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

    /// The same force field with only the terms that touch at least one atom `keep` marks, and no
    /// solvation: every term left out is among unmarked atoms alone. The terms kept are in the
    /// same order, so the force on a marked atom is the same sum, in the same order, as the whole
    /// force field's — bit for bit — and the energy differs by the left-out terms, which are a
    /// constant while the unmarked atoms do not move.
    ///
    /// # Panics
    ///
    /// If `keep` is not one per atom.
    pub(crate) fn touching(&self, keep: &[bool]) -> ForceField {
        assert_eq!(keep.len(), self.charges.len(), "one mark per atom");
        let any = |atoms: &[usize]| atoms.iter().any(|&a| keep[a]);
        ForceField {
            stretches: self
                .stretches
                .iter()
                .filter(|t| any(&t.atoms))
                .copied()
                .collect(),
            bends: self
                .bends
                .iter()
                .filter(|t| any(&t.atoms))
                .copied()
                .collect(),
            torsions: self
                .torsions
                .iter()
                .filter(|t| any(&t.atoms))
                .copied()
                .collect(),
            inversions: self
                .inversions
                .iter()
                .filter(|t| any(&t.atoms))
                .copied()
                .collect(),
            pairs: self
                .pairs
                .iter()
                .filter(|t| any(&t.atoms))
                .copied()
                .collect(),
            charges: self.charges.clone(),
            elements: self.elements.clone(),
            total_charge: self.total_charge,
            solvation: None,
        }
    }

    /// The same force field without the non-bonded pairs that join an atom `group` marks to one it
    /// does not, and those pairs, in the order the force field holds them: what
    /// [`crate::alchemy::Decoupling`] scales. Every other term is kept as it is, the group's own
    /// pairs included.
    ///
    /// # Errors
    ///
    /// The two atoms of the first bond that joins the group to the rest: a group bonded to its
    /// environment cannot be decoupled by scaling its non-bonded pairs alone.
    ///
    /// # Panics
    ///
    /// If `group` is not one per atom.
    pub(crate) fn split_across(
        &self,
        group: &[bool],
    ) -> Result<(ForceField, Vec<Pair>), [usize; 2]> {
        assert_eq!(group.len(), self.charges.len(), "one mark per atom");
        if let Some(s) = self
            .stretches
            .iter()
            .find(|s| group[s.atoms[0]] != group[s.atoms[1]])
        {
            return Err(s.atoms);
        }
        let (across, within): (Vec<Pair>, Vec<Pair>) = self
            .pairs
            .iter()
            .copied()
            .partition(|p| group[p.atoms[0]] != group[p.atoms[1]]);
        let mut rest = self.clone();
        rest.pairs = within;
        Ok((rest, across))
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
        if let Some(gb) = &self.solvation {
            let g = gb.accumulate(&self.charges, at, &mut forces);
            e.solvation += g;
            e.total += g;
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
