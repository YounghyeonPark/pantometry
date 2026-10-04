//! The binding energy of a ligand in a rigid protein pocket, and the ligand minimised there.
//!
//! [`Binding::new`] cuts a pocket out of a [`System`] around its ligand, gives the pocket and the
//! ligand their charges, and builds three force fields — the complex, the pocket alone and the
//! ligand alone — on the same atoms at the same positions. Then
//!
//! `ΔE_bind = E(complex) − E(pocket) − E(ligand)`
//!
//! at fixed geometry ([`Binding::interaction`]), the same difference of OBC II solvation energies
//! and of solvent-accessible areas ([`Binding::desolvation`]), and the ligand relaxed in the pocket
//! with every protein atom held ([`Binding::minimise_ligand`]). [`Binding::ligand_at`] moves the
//! ligand rigidly, which is what a path out of the pocket and a congener series need.
//! [`Binding::relaxing_hydrogens`] relaxes each of the three systems' own hydrogens on its frozen
//! heavy atoms, and every one of those quantities is then taken with each system at its own
//! positions.
//!
//! # The choices, and why
//!
//! **The pocket is whole residues.** A protein residue is in it when any of its heavy atoms is
//! within the cutoff of any ligand atom, hydrogens included on the ligand's side — the rule the
//! 2c-1 QEq test used, so its counts (6 Å: 18 residues and benzene, 322 atoms, +1) are this
//! module's too. Heavy atoms on the protein's side because they are what the crystal measured; the
//! hydrogens were placed. **Whole residues** so that no residue is split: its formal charge is an
//! integer, QEq gets a total it can be asked for, and no side chain is cut through a bond. Peptide
//! bonds to residues outside the pocket are cut; the terms that cut changes are among pocket atoms
//! only, so they are the same in the complex and the pocket alone and cancel exactly in ΔE_bind.
//! **The pocket is chosen once, around the ligand where it is when the binding is built** — the
//! crystal pose — and does not follow the ligand when [`Binding::ligand_at`] moves it: a ligand
//! moved out of the pocket meets no residue that was not in it. A path that leaves the pocket wants
//! a cutoff that covers the path.
//!
//! **The types are the whole system's.** [`uff::assign`] is run once on the whole [`System`] and
//! each atom keeps its type in every fragment, so a backbone nitrogen whose peptide partner was
//! cut away is still the amide `N_R` it is in the protein. Every term the types could change at
//! a cut is pocket-internal and cancels; what reaches ΔE_bind is each atom's van der Waals
//! parameters, which in UFF depend on the element alone for every type this system has (`O_2` and
//! `O_3` are both 3.500 Å and 0.060 kcal/mol, so the carboxylate typed one of each in 2c-1 binds as
//! two equivalent oxygens would).
//!
//! **Charges: QEq on the pocket and on the ligand, separately, each at its own total formal
//! charge** — fixed-charge force-field practice, in which a ligand's charges are its own and the
//! protein's are the protein's. **Polarisation of either partner by the other is therefore
//! excluded, and so is charge transfer between them.** Both are rigid-motion invariant, so a
//! ligand moved rigidly keeps exactly the charges it should have there, and nothing is
//! recomputed along a path. [`Binding::polarised_charges`] gives QEq on the whole complex, to
//! measure what the choice leaves out; the tests report it.
//!
//! **ΔE_bind in vacuum is exactly the protein–ligand cross terms.** The protein is rigid and the
//! ligand has no bond to it, so every bond, angle, torsion and inversion term and every
//! pocket–pocket and ligand–ligand pair is the same in the complex as in its fragment, and only
//! the pairs across the partition are left: van der Waals plus Coulomb, every pair, with no
//! exclusion (there is no bond to exclude by). [`Interaction`] carries both — the cross sum,
//! which has no cancellation in it, and the three-evaluation difference term by term, which
//! does — and the identity is tested.
//!
//! **Solvation does not decompose.** Generalized Born is many-body through the Born radii: the
//! ligand descreens the pocket's atoms and the pocket the ligand's, so the difference
//! `ΔG_GB(complex) − ΔG_GB(pocket) − ΔG_GB(ligand)` holds every pocket atom's self-energy change
//! as well as the screened cross terms. [`Desolvation`] reports the polar difference and the
//! nonpolar `0.005 kcal mol⁻¹ Å⁻² × ΔSASA` separately. **The pocket alone is solvated as a
//! fragment cut out of a protein**, so its absolute ΔG_GB is meaningless; only the difference is
//! read, and how far it has converged with the cutoff is measured, not assumed.
//!
//! **Minimisation in the pocket** holds every protein atom by [`Minimiser::with_frozen`], and
//! minimises the ligand's own energy plus the interaction, **in vacuum**, at the fixed charges.
//! The force field it minimises is the complex's with the terms among protein atoms alone left
//! out: they are a constant while the protein does not move, and leaving them out keeps a
//! thousand-kcal/mol constant out of the line search's arithmetic. The force on a ligand atom is
//! then bit for bit the whole complex's. Vacuum because the polar GB term on a truncated pocket is
//! the cutoff-sensitive one, and minimising on it would move the pose by an artefact of the cut.
//!
//! # Hydrogens relaxed in each system
//!
//! **Why.** Every hydrogen here was placed, not measured: the protein's by 2c-1's rules, whose rotor
//! step clears heavy atoms only, and the ligand's by superposing its dictionary template. Neither
//! placement sees the other partner's hydrogens, and at the crystal pose three of the nine T4
//! lysozyme L99A complexes of `tests/the_congener_series.rs` clash hydrogen to hydrogen: indene's
//! H12 is 1.41 Å from Val111's HG13 in 183L, and van der Waals is +224 kcal/mol. A crystal-pose
//! ΔE_bind there measures the placement.
//!
//! **What [`Binding::relaxing_hydrogens`] does.** In each of the three systems — the complex, the
//! pocket alone and the ligand alone — it minimises that system's own energy over its free
//! hydrogens with every other atom frozen ([`Minimiser::with_frozen`]), and then
//!
//! `ΔE_bind = E(complex, H relaxed) − E(pocket, its own H relaxed) − E(ligand, its own H relaxed)`.
//!
//! **Why this is the consistent definition.** ΔE_bind is the energy of the reaction pocket +
//! ligand → complex, and each side of a reaction is taken at its own minimum. Relaxing the complex
//! alone and taking the fragments at the complex's positions — what 2c-3a reported — leaves in the
//! fragments' energies whatever the complex's hydrogens bent to make room for the partner, and
//! calls that strain binding: the apo pocket would not have it. With each fragment at its own
//! minimum, the difference is the cross terms at the complex's positions plus a **reorganisation**,
//! `[E(pocket) at the complex's positions − at its own] + [the same for the ligand]`, the price of
//! the bound hydrogen arrangement, which [`Interaction::reorganisation`] reports and the tests
//! check is exactly that.
//!
//! **Each fragment starts from where the complex left it.** A minimiser finds the minimum downhill
//! from its start, and two different starts on the same flat landscape — a hydroxyl's or a
//! methyl's rotor far from the ligand — can stop at different minima for reasons that have nothing
//! to do with binding, which would put a rotamer's energy into ΔE_bind as noise. Started from the
//! complex, a fragment moves only where the partner's removal pushes it, so the reorganisation is
//! never negative (the minimiser does not raise the energy it minimises) and vanishes when the
//! partner exerts no force: far away, each fragment takes no step at all and its positions are the
//! complex's to the bit, which the tests assert. What this does not do is search: a fragment whose
//! own best arrangement is in another basin than the complex's is not found there.
//!
//! **Which hydrogens are free** ([`Binding::free_hydrogens`]): every ligand hydrogen, and each
//! pocket hydrogen whose heavy parent is within the cutoff of a ligand atom where the binding was
//! built — the pocket's own rule, atom by atom. The rest, on the outer side of residues the cutoff
//! reaches, are held where 2c-1 placed them in all three systems. **This is a truncation, and it
//! has a measured cost.** The held hydrogens do feel the ligand — at the relaxed complex the
//! partner's force on them is 0.001–0.14 kcal mol⁻¹ Å⁻¹ on four entries measured, nearly all of it
//! above the tolerance — so holding them leaves part of the pocket's response out. On benzene at
//! 6 Å, 106 of the pocket's 167 hydrogens are free; freeing all 167 moves the relaxed ΔE_bind from
//! −16.14 to −16.39 kcal/mol and the reorganisation from +2.62 to +4.03, and the minimised ΔE_bind
//! from −22.55 to −22.54. On the other eight at 6 Å the relaxed ΔE_bind moves by at most 0.10, and
//! at 8 Å by at most 0.03. The rule is kept because it is the pocket's own and because the outer
//! hydrogens relax against the cut as much as against the ligand — a backbone amide whose peptide
//! partner is gone, a rotor facing missing neighbours — not because they carry nothing.
//! [`Binding::relaxing_hydrogens_of`] takes any mask, which is how the cost was measured.
//!
//! **The tolerance is [`Binding::HYDROGEN_TOLERANCE`], 2e-3 kcal mol⁻¹ Å⁻¹, not the 1e-4 a
//! ligand alone is minimised to**, because at 1e-4 the minimiser stalled on the series, where the
//! rounding of a sum of 5·10⁴ terms or more can hide the last steps' decrease. That constant gives
//! the evidence, which is empirical, and the test that holds what the looser tolerance costs.
//!
//! **The charges are not recomputed**: QEq's, at the positions the binding was built with, as
//! everywhere in this module. **Minimised from a relaxed complex**, the pocket's free hydrogens move
//! with the ligand and each fragment is relaxed again afterwards ([`Binding::minimise_ligand`]), so
//! the binding stays relaxed; moving the ligand rigidly ([`Binding::ligand_at`]) gives an unrelaxed
//! binding, since the complex's hydrogens were relaxed against the ligand where it was.
//!
//! # What it measures, on benzene in T4 lysozyme L99A (PDB 181L)
//!
//! `tests/benzene_in_its_pocket.rs`, kcal/mol. At the crystal pose in the 6 Å pocket, ΔE_bind in
//! vacuum is −12.43: van der Waals −11.30, Coulomb −1.13. The polar desolvation is **+7.32, and
//! it is not small**. Generalized Born is quadratic in the charges, so the term splits exactly:
//! +2.27 is the pocket desolvated by benzene's volume (benzene's charges set to zero), +1.94 is
//! benzene's own (of the 2.27 it has alone), and +3.11 is the screened cross terms. GB counts the
//! empty apo cavity as solvent. The buried area is −294.6 Å², nonpolar −1.47. **The polar term
//! does not converge with the cutoff**: +7.32, +10.61, +12.21 at 6, 8 and 10 Å, and +14.96 with the
//! whole protein. Van der Waals goes −11.30, −13.05, −13.33 and −13.50. Benzene minimised in the
//! rigid pocket moves 0.81 Å heavy-atom RMSD (0.82 Å for the whole protein), and ΔE_bind goes from
//! −12.43 to −22.35. Freeing the protein's hydrogens as well moves it 0.83 Å, so the move is UFF's
//! and not an artefact of where 2c-1 placed the hydrogens — although the crystal pose does have a
//! 1.89 Å contact between a benzene H and Val111's HG13, which that rule, scoring heavy atoms
//! only, could not see.
//!
//! **With each system's hydrogens relaxed** the closest benzene–protein H–H is 2.07 Å, and ΔE_bind
//! in vacuum is −16.14: cross terms −18.76 (van der Waals −17.17, Coulomb −1.58) and a
//! reorganisation of +2.62 (pocket +1.06, benzene +1.56). The polar desolvation is +8.34 and the
//! buried area −291.4 Å². Minimised from there, benzene moves 0.90 Å RMSD, against 0.81 from the
//! crystal's hydrogens, and ΔE_bind is −22.55 against −22.35, of which +0.57 is reorganisation.
//!
//! # What it is not
//!
//! A binding *energy* at one geometry, not a free energy: no entropy, no protein flexibility
//! beyond its hydrogens, no ligand strain against its own solution minimum (the ligand alone keeps
//! its bound heavy atoms), no water in the pocket but GB's continuum.
//! With QEq charges, OBC II over-solvates small molecules (see [`crate::solvation`]), so the
//! solvated estimate is not quantitative.
//!
//! [`uff::assign`]: crate::uff::assign

use crate::ccd::{Atom, Bond, Component, Element, ANGSTROM};
use crate::energy::{coulomb, Energy, ForceField};
use crate::minimise::{Minimiser, Progress, Status};
use crate::pdb::{Part, System};
use crate::qeq::{Charges, Qeq, QeqError};
use crate::solvation::{
    intrinsic_radius, nonpolar_energy, surface_area, GeneralizedBorn, PROBE_RADIUS,
};
use crate::uff::{self, UffType};
use std::fmt;
use std::sync::Arc;

/// A rigid motion of a set of atoms about its own centroid `c`: `p ↦ c + R (p − c) + t`. About the
/// centroid, so that a rotation turns a ligand in place rather than swinging it across the pocket,
/// and a translation is where the centroid goes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidMotion {
    /// R, a proper rotation, row-major.
    pub rotation: [[f64; 3]; 3],
    /// t, metres.
    pub translation: [f64; 3],
}

impl RigidMotion {
    /// No motion.
    pub const IDENTITY: RigidMotion = RigidMotion {
        rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        translation: [0.0; 3],
    };

    /// A translation by `t` metres.
    pub fn translation(t: [f64; 3]) -> RigidMotion {
        RigidMotion {
            translation: t,
            ..RigidMotion::IDENTITY
        }
    }

    /// A rotation by `angle` radians about `axis` (right-handed), by Rodrigues' formula.
    ///
    /// # Panics
    ///
    /// If `axis` is zero or not finite.
    pub fn rotation(axis: [f64; 3], angle: f64) -> RigidMotion {
        let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        assert!(n > 0.0 && n.is_finite(), "a rotation needs a non-zero axis");
        let [x, y, z] = [axis[0] / n, axis[1] / n, axis[2] / n];
        let (s, c) = angle.sin_cos();
        let t = 1.0 - c;
        RigidMotion {
            rotation: [
                [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
                [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
                [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
            ],
            translation: [0.0; 3],
        }
    }

    /// The same rotation followed by a translation by `t` metres more.
    pub fn then_translate(mut self, t: [f64; 3]) -> RigidMotion {
        for (a, b) in self.translation.iter_mut().zip(t) {
            *a += b;
        }
        self
    }

    /// `c + R (p − c) + t`.
    pub fn apply(&self, centre: [f64; 3], p: [f64; 3]) -> [f64; 3] {
        let d = [p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]];
        let r = &self.rotation;
        let mut out = [0.0; 3];
        for (k, o) in out.iter_mut().enumerate() {
            *o = centre[k]
                + (r[k][0] * d[0] + r[k][1] * d[1] + r[k][2] * d[2])
                + self.translation[k];
        }
        out
    }
}

/// Why a binding could not be built.
#[derive(Clone, Debug, PartialEq)]
pub enum BindingError {
    /// The system has no ligand atom.
    NoLigand,
    /// No protein residue is within the cutoff of the ligand.
    EmptyPocket {
        /// The cutoff, metres.
        cutoff: f64,
    },
    /// A fragment's force field was refused — a hypervalent sulfur.
    Unsupported(crate::energy::Unsupported),
    /// QEq failed on the pocket or the ligand.
    Charges {
        /// `"pocket"` or `"ligand"`.
        fragment: &'static str,
        /// Why.
        error: QeqError,
    },
}

impl fmt::Display for BindingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BindingError::NoLigand => f.write_str("the system has no ligand atom"),
            BindingError::EmptyPocket { cutoff } => write!(
                f,
                "no protein residue has a heavy atom within {} Å of the ligand",
                cutoff / ANGSTROM
            ),
            BindingError::Unsupported(u) => write!(f, "{u}"),
            BindingError::Charges { fragment, error } => {
                write!(f, "QEq on the {fragment}: {error}")
            }
        }
    }
}

impl std::error::Error for BindingError {}

/// ΔE_bind in vacuum, joules per molecule: at one geometry shared by the three systems, or, for a
/// binding whose hydrogens are relaxed ([`Binding::relaxing_hydrogens`]), with each system at its
/// own.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Interaction {
    /// The van der Waals energy of every protein–ligand pair at the complex's positions, summed
    /// directly.
    pub van_der_waals: f64,
    /// The Coulomb energy (ε = 1) of every protein–ligand pair at the complex's positions, summed
    /// directly.
    pub electrostatic: f64,
    /// What each fragment's own energy loses when its hydrogens go from where they are in the
    /// complex to where they relax with the partner gone, term by term: `[E(pocket) at the
    /// complex's positions − E(pocket) at its own] + [the same for the ligand]`. **Exactly zero
    /// in every field, and not computed, when the three systems share positions** — the crystal
    /// pose, or any binding not relaxed. Never negative in its total for a relaxed binding, up to
    /// rounding: each fragment's relaxation starts from the complex's positions and the minimiser
    /// does not raise the energy it minimises.
    pub reorganisation: Energy,
    /// `E(complex) − E(pocket) − E(ligand)`, term by term, from three whole evaluations, each at
    /// its own system's positions: the definition. In exact arithmetic it is the cross terms in the
    /// van der Waals and Coulomb fields plus [`Interaction::reorganisation`] in every field. In
    /// floating point it carries the rounding of the pocket's whole energy, which the direct sums
    /// do not.
    pub difference: Energy,
}

impl Interaction {
    /// `van_der_waals + electrostatic + reorganisation.total`: ΔE_bind from the cross terms and the
    /// fragments' own hydrogens. For a binding at shared positions the last is exactly zero, and
    /// this is the cross terms alone, to the bit.
    pub fn total(&self) -> f64 {
        self.van_der_waals + self.electrostatic + self.reorganisation.total
    }
}

/// Where the three hydrogen relaxations of a [`Binding`] stopped: see
/// [`Binding::relaxing_hydrogens`].
#[derive(Clone, Debug, PartialEq)]
pub struct HydrogenRelaxation {
    /// Which atoms were free, complex order: hydrogens only. Every other atom was frozen in every
    /// system.
    pub free: Vec<bool>,
    /// The tolerance each relaxation was asked to converge to, newtons.
    pub tolerance: f64,
    /// The most steps each relaxation was allowed.
    pub max_steps: usize,
    /// The complex's: its free hydrogens relaxed — or, after [`Binding::minimise_ligand`], the
    /// ligand minimised with them.
    pub complex: Progress,
    /// The pocket alone's free hydrogens, relaxed from where the complex left them.
    pub pocket: Progress,
    /// The ligand alone's hydrogens, relaxed from where the complex left them.
    pub ligand: Progress,
}

/// The fragments' own positions, when they are not the complex's.
#[derive(Clone, Debug)]
struct Apart {
    relaxation: HydrogenRelaxation,
    pocket: Vec<[f64; 3]>,
    ligand: Vec<[f64; 3]>,
}

/// The solvation part of binding, each system at its own positions — one geometry, unless the
/// hydrogens are relaxed: OBC II and the nonpolar surface term, each as complex − pocket −
/// ligand. Energies in joules per molecule, areas in m².
#[derive(Clone, Debug, PartialEq)]
pub struct Desolvation {
    /// ΔG_GB of the complex.
    pub complex: f64,
    /// ΔG_GB of the pocket alone.
    pub pocket: f64,
    /// ΔG_GB of the ligand alone.
    pub ligand: f64,
    /// `complex − pocket − ligand`: the polar desolvation, positive when binding costs
    /// solvation.
    pub polar: f64,
    /// Each atom's solvent-accessible area in the complex, complex order.
    pub areas_bound: Vec<f64>,
    /// Each atom's area in its own fragment at the fragment's own positions, complex order.
    pub areas_apart: Vec<f64>,
    /// ΔSASA, `Σ (bound − apart)`: the area binding buries, negative.
    pub buried_area: f64,
    /// `0.005 kcal mol⁻¹ Å⁻² × ΔSASA` (OBC p. 392).
    pub nonpolar: f64,
}

/// What every [`Binding`] of one pocket shares: built once, behind an `Arc`, so that moving the
/// ligand copies positions and nothing else.
#[derive(Debug)]
struct Model {
    cutoff: f64,
    residues: Vec<String>,
    system_atoms: Vec<usize>,
    pocket_len: usize,
    elements: Vec<Element>,
    types: Vec<UffType>,
    charges: Vec<f64>,
    pocket_charge: i32,
    ligand_charge: i32,
    qeq_solves: [usize; 2],
    complex: ForceField,
    pocket: ForceField,
    ligand: ForceField,
    in_pocket: ForceField,
    crystal: Vec<[f64; 3]>,
    near_hydrogens: Vec<bool>,
}

/// A ligand in a rigid pocket: see the module documentation.
///
/// Atoms are in **complex order**: the pocket's atoms first, in the system's order, then the
/// ligand's. [`Binding::system_atoms`] maps each back to the [`System`].
#[derive(Clone, Debug)]
pub struct Binding {
    model: Arc<Model>,
    at: Vec<[f64; 3]>,
    apart: Option<Apart>,
}

impl Binding {
    /// The tolerance hydrogens are relaxed to here, newtons: 2e-3 kcal mol⁻¹ Å⁻¹, twenty times
    /// [`crate::Molecule::DEFAULT_TOLERANCE`]. **The justification is empirical.** At 1e-4 the
    /// minimiser stalled four times on the congener series — no step lowered the energy in floating
    /// point — at 1.6–4.3e-4 kcal mol⁻¹ Å⁻¹, while 1e-4 did converge on benzene's 6 Å pocket. At
    /// 2e-3 every one of the series' relaxations and minimisations from a relaxed complex converged
    /// (the nine entries at 6 and 8 Å), and on benzene at 6 Å the result agrees with 1e-4's to
    /// within half the last digit printed: a test asserts it (ΔE_bind and the reorganisation to
    /// 5e-3 kcal/mol, the minimised pose to 5e-4 Å).
    ///
    /// **Why a pocket stalls, as a scale and not a bound.** The line search accepts a step only
    /// when the energy falls, and a pocket's energy is a sum of 0.5–4·10⁵ terms whose rounding moves
    /// it under a 1e-9 Å nudge by an amount δE that varies from nudge to nudge: a median of 3.6e-12
    /// kcal/mol over 300 nudges of benzene's 6 Å pocket and a maximum of 2.1e-11 there, single draws
    /// of up to 1.8e-10 on the 8 Å pockets, and the largest of 300 about six times the first.
    /// Removing a residual force `g` along an N–H or O–H stretch (`k` ≈ 1190 kcal mol⁻¹ Å⁻²,
    /// measured) lowers the energy by `g²/2k`, so `√(4 k δE)` — from 1.3e-4 up to about 2e-3
    /// kcal mol⁻¹ Å⁻¹ over that spread of δE — is
    /// the order of the force below which the decrease can drown in the noise. That is where the
    /// stalls were, and it does not say where any one run will stop.
    pub const HYDROGEN_TOLERANCE: f64 = 2e-3 * crate::minimise::KCAL_PER_MOL_ANGSTROM;

    /// The pocket of `system` within `cutoff` metres of its ligand, at the system's coordinates,
    /// with QEq charges on the pocket and the ligand separately. See the module documentation for
    /// every choice.
    ///
    /// Costs two QEq solves, dense and O(N³) in the pocket's atoms: 1.4 s with `--release` at
    /// 6 Å (310 + 12 atoms), about 17 s at 10 Å (987) and 195 s for all of 181L (2616).
    ///
    /// # Errors
    ///
    /// [`BindingError`], naming what failed.
    pub fn new(system: &System, cutoff: f64) -> Result<Binding, BindingError> {
        let whole = system.component();
        let atoms = whole.atoms();
        let ligand = system.atoms_in(Part::Ligand);
        if ligand.is_empty() {
            return Err(BindingError::NoLigand);
        }
        let near = |i: usize| {
            ligand.iter().any(|&l| {
                let (a, b) = (atoms[i].at, atoms[l].at);
                let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() < cutoff
            })
        };
        let mut residues = Vec::new();
        let mut pocket = Vec::new();
        for r in system.residues() {
            if r.part == Part::Protein
                && r.atoms
                    .clone()
                    .any(|i| atoms[i].element != Element::H && near(i))
            {
                residues.push(r.label());
                pocket.extend(r.atoms.clone());
            }
        }
        if pocket.is_empty() {
            return Err(BindingError::EmptyPocket { cutoff });
        }
        let all_types = uff::assign(whole);
        let pocket_len = pocket.len();
        let system_atoms: Vec<usize> = pocket.iter().chain(&ligand).copied().collect();

        let fragment = |from: &[usize], name: &str| {
            let mut index = vec![usize::MAX; atoms.len()];
            for (k, &i) in from.iter().enumerate() {
                index[i] = k;
            }
            let kept: Vec<Atom> = from.iter().map(|&i| atoms[i].clone()).collect();
            let bonds: Vec<Bond> = whole
                .bonds()
                .iter()
                .filter(|b| index[b.atoms[0]] != usize::MAX && index[b.atoms[1]] != usize::MAX)
                .map(|b| Bond {
                    atoms: [index[b.atoms[0]], index[b.atoms[1]]],
                    ..b.clone()
                })
                .collect();
            let types: Vec<UffType> = from.iter().map(|&i| all_types[i]).collect();
            let component = Component::assembled(name.into(), kept, bonds);
            (component, types)
        };
        let (pc, pt) = fragment(&pocket, "pocket");
        let (lc, lt) = fragment(&ligand, "ligand");
        let (cc, ct) = fragment(&system_atoms, "complex");

        let elements: Vec<Element> = system_atoms.iter().map(|&i| atoms[i].element).collect();
        let crystal: Vec<[f64; 3]> = system_atoms.iter().map(|&i| atoms[i].at).collect();
        let charge_of = |c: &Component| c.atoms().iter().map(|a| a.charge).sum::<i32>();
        let (pocket_charge, ligand_charge) = (charge_of(&pc), charge_of(&lc));
        let solve = |fragment: &'static str, range: std::ops::Range<usize>, total: i32| {
            Qeq::default()
                .equilibrate(&elements[range.clone()], &crystal[range], f64::from(total))
                .map_err(|error| BindingError::Charges { fragment, error })
        };
        let qp = solve("pocket", 0..pocket_len, pocket_charge)?;
        let ql = solve("ligand", pocket_len..system_atoms.len(), ligand_charge)?;
        let charges: Vec<f64> = qp.charges.iter().chain(&ql.charges).copied().collect();

        let build = |c: &Component, t: &[UffType], q: Vec<f64>| {
            ForceField::new(c, t)
                .map(|ff| ff.with_charges(q))
                .map_err(BindingError::Unsupported)
        };
        let complex = build(&cc, &ct, charges.clone())?;
        let pocket_ff = build(&pc, &pt, qp.charges.clone())?;
        let ligand_ff = build(&lc, &lt, ql.charges.clone())?;
        let moving: Vec<bool> = (0..system_atoms.len()).map(|k| k >= pocket_len).collect();
        let in_pocket = complex.touching(&moving);
        // A hydrogen is near when it is the ligand's, or its heavy parent is within the cutoff of
        // the ligand: the pocket's own rule, atom by atom.
        let near_hydrogens: Vec<bool> = system_atoms
            .iter()
            .enumerate()
            .map(|(k, &i)| {
                atoms[i].element == Element::H
                    && (k >= pocket_len || whole.neighbours(i).any(|(j, _)| near(j)))
            })
            .collect();
        Ok(Binding {
            apart: None,
            at: crystal.clone(),
            model: Arc::new(Model {
                cutoff,
                residues,
                system_atoms,
                pocket_len,
                elements,
                types: ct,
                charges,
                pocket_charge,
                ligand_charge,
                qeq_solves: [qp.iterations, ql.iterations],
                complex,
                pocket: pocket_ff,
                ligand: ligand_ff,
                in_pocket,
                crystal,
                near_hydrogens,
            }),
        })
    }

    /// The cutoff the pocket was cut at, metres.
    pub fn cutoff(&self) -> f64 {
        self.model.cutoff
    }

    /// The pocket's residues, `A:LEU84`, in the system's order.
    pub fn residues(&self) -> &[String] {
        &self.model.residues
    }

    /// For each atom in complex order, its index in the [`System`].
    pub fn system_atoms(&self) -> &[usize] {
        &self.model.system_atoms
    }

    /// How many atoms the pocket has: complex indices `0..pocket_len()` are the pocket's.
    pub fn pocket_len(&self) -> usize {
        self.model.pocket_len
    }

    /// The ligand's atoms, as complex indices.
    pub fn ligand_range(&self) -> std::ops::Range<usize> {
        self.model.pocket_len..self.at.len()
    }

    /// Every atom's element, complex order.
    pub fn elements(&self) -> &[Element] {
        &self.model.elements
    }

    /// Every atom's UFF type, complex order — the whole system's typing.
    pub fn types(&self) -> &[UffType] {
        &self.model.types
    }

    /// Every atom's partial charge, e, complex order: QEq on the pocket, then QEq on the ligand.
    pub fn charges(&self) -> &[f64] {
        &self.model.charges
    }

    /// The pocket's and the ligand's total formal charges, e.
    pub fn formal_charges(&self) -> (i32, i32) {
        (self.model.pocket_charge, self.model.ligand_charge)
    }

    /// The QEq solves the pocket's and the ligand's charges took.
    pub fn qeq_solves(&self) -> [usize; 2] {
        self.model.qeq_solves
    }

    /// The complex's force field, with the charges.
    pub fn force_field(&self) -> &ForceField {
        &self.model.complex
    }

    /// The pocket's force field alone, its atoms the complex's first [`Binding::pocket_len`].
    pub fn pocket_force_field(&self) -> &ForceField {
        &self.model.pocket
    }

    /// The ligand's force field alone, its atoms the complex's [`Binding::ligand_range`] in order.
    pub fn ligand_force_field(&self) -> &ForceField {
        &self.model.ligand
    }

    /// Every atom's position now, metres, complex order.
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.at
    }

    /// The ligand's positions now, metres.
    pub fn ligand_positions(&self) -> &[[f64; 3]] {
        &self.at[self.model.pocket_len..]
    }

    /// The pocket alone's positions, metres, at which its energy, its solvation and its areas are
    /// taken: the complex's first [`Binding::pocket_len`], bit for bit, unless the hydrogens are
    /// relaxed, and then the pocket's own relaxed hydrogens on the same heavy atoms.
    pub fn pocket_alone_positions(&self) -> &[[f64; 3]] {
        match &self.apart {
            Some(a) => &a.pocket,
            None => &self.at[..self.model.pocket_len],
        }
    }

    /// The ligand alone's positions, metres: [`Binding::ligand_positions`], bit for bit, unless
    /// the hydrogens are relaxed, and then the ligand's own relaxed hydrogens on the same heavy
    /// atoms.
    pub fn ligand_alone_positions(&self) -> &[[f64; 3]] {
        match &self.apart {
            Some(a) => &a.ligand,
            None => self.ligand_positions(),
        }
    }

    /// The hydrogens [`Binding::relaxing_hydrogens`] frees, complex order: every hydrogen of the
    /// ligand, and every hydrogen of the pocket whose heavy parent is within the cutoff of a ligand
    /// atom where the binding was built — the pocket's own rule, applied to atoms instead of
    /// residues. See the module documentation for why the rest are held.
    pub fn free_hydrogens(&self) -> &[bool] {
        &self.model.near_hydrogens
    }

    /// Where the hydrogen relaxations stopped, if this binding's hydrogens are relaxed.
    pub fn hydrogen_relaxation(&self) -> Option<&HydrogenRelaxation> {
        self.apart.as_ref().map(|a| &a.relaxation)
    }

    /// The ligand's positions in the system it was built from — the crystal pose, for an entry.
    pub fn crystal_ligand(&self) -> &[[f64; 3]] {
        &self.model.crystal[self.model.pocket_len..]
    }

    /// The centroid of the ligand's atoms now, unweighted, metres.
    pub fn ligand_centroid(&self) -> [f64; 3] {
        centroid(self.ligand_positions())
    }

    /// The RMS distance of the ligand's heavy atoms from their crystal positions, metres, **without
    /// superposition**: the pocket is the frame, so a ligand that slid or turned in it has moved.
    pub fn ligand_rmsd(&self) -> f64 {
        let mut sum = 0.0;
        let mut count = 0usize;
        for k in self.ligand_range() {
            if self.model.elements[k] != Element::H {
                let (a, b) = (self.at[k], self.model.crystal[k]);
                sum += (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2);
                count += 1;
            }
        }
        (sum / count.max(1) as f64).sqrt()
    }

    /// The same binding with the ligand moved by `motion` about its current centroid; the pocket
    /// and every charge unchanged. Cheap: the force fields are shared.
    ///
    /// **Not relaxed**, even when `self` is: the complex's hydrogens were relaxed against the
    /// ligand where it was, so the moved binding's three systems share its positions again, as an
    /// unrelaxed binding's do, and [`Binding::relaxing_hydrogens`] relaxes them at the new pose.
    pub fn ligand_at(&self, motion: &RigidMotion) -> Binding {
        let c = self.ligand_centroid();
        let mut at = self.at.clone();
        for p in &mut at[self.model.pocket_len..] {
            *p = motion.apply(c, *p);
        }
        Binding {
            model: Arc::clone(&self.model),
            at,
            apart: None,
        }
    }

    /// The same binding with the ligand's atoms at `at` (metres, ligand order) — a conformer, or
    /// one atom displaced. The pocket and every charge unchanged. Not relaxed, as
    /// [`Binding::ligand_at`] is not.
    ///
    /// # Panics
    ///
    /// If `at` is not one position per ligand atom.
    pub fn with_ligand_positions(&self, at: &[[f64; 3]]) -> Binding {
        let n0 = self.model.pocket_len;
        assert_eq!(at.len(), self.at.len() - n0, "one position per ligand atom");
        let mut all = self.at.clone();
        all[n0..].copy_from_slice(at);
        Binding {
            model: Arc::clone(&self.model),
            at: all,
            apart: None,
        }
    }

    /// The same binding with **every** atom, protein and ligand, moved by `motion` about the
    /// ligand's centroid: the same complex in another frame, whose ΔE_bind is the same — which the
    /// tests check. A relaxed binding's fragments move with it and it stays relaxed.
    pub fn moved(&self, motion: &RigidMotion) -> Binding {
        let c = self.ligand_centroid();
        let turn = |at: &[[f64; 3]]| at.iter().map(|&p| motion.apply(c, p)).collect();
        Binding {
            model: Arc::clone(&self.model),
            at: turn(&self.at),
            apart: self.apart.as_ref().map(|a| Apart {
                relaxation: a.relaxation.clone(),
                pocket: turn(&a.pocket),
                ligand: turn(&a.ligand),
            }),
        }
    }

    /// The same binding with **each system's own hydrogens relaxed**, every heavy atom frozen in
    /// all three, the hydrogens [`Binding::free_hydrogens`] marks free: see
    /// [`Binding::relaxing_hydrogens_of`].
    pub fn relaxing_hydrogens(&self, max_steps: usize, tolerance: f64) -> Binding {
        self.relaxing_hydrogens_of(&self.model.near_hydrogens, max_steps, tolerance)
    }

    /// The same binding with each system's own hydrogens relaxed, the ones `free` marks (complex
    /// order) free and every other atom frozen ([`Minimiser::with_frozen`]), each relaxation for at
    /// most `max_steps` steps or until the largest force on a free atom is at most `tolerance`
    /// newtons. In order:
    ///
    /// 1. **the complex**, from its current positions;
    /// 2. **the pocket alone** and **the ligand alone**, each from where the complex left its
    ///    atoms, on its own force field.
    ///
    /// Then `ΔE_bind = E(complex) − E(pocket) − E(ligand)` is taken with each system at its own
    /// relaxed positions ([`Binding::interaction`]), and so are the solvation and the areas
    /// ([`Binding::desolvation`]). See the module documentation for why this is the consistent
    /// definition and why the fragments start from the complex. The charges are not recomputed:
    /// they are QEq's at the positions the binding was built with, as everywhere in this module.
    /// [`HydrogenRelaxation`] says where each stopped; nothing here asserts that they converged.
    /// **A tolerance below [`Binding::HYDROGEN_TOLERANCE`] may stall** rather than converge, for the
    /// reason given there.
    ///
    /// # Panics
    ///
    /// If `free` is not one per atom, or marks an atom that is not a hydrogen.
    pub fn relaxing_hydrogens_of(
        &self,
        free: &[bool],
        max_steps: usize,
        tolerance: f64,
    ) -> Binding {
        assert_eq!(free.len(), self.at.len(), "one mark per atom");
        for (k, &f) in free.iter().enumerate() {
            assert!(
                !f || self.model.elements[k] == Element::H,
                "atom {k} is marked free and is not a hydrogen"
            );
        }
        let mut at = self.at.clone();
        let complex = relax(&self.model.complex, free, &mut at, max_steps, tolerance);
        let mut b = Binding {
            model: Arc::clone(&self.model),
            at,
            apart: None,
        };
        b.relax_fragments(free.to_vec(), complex, max_steps, tolerance);
        b
    }

    /// Relaxes the pocket alone and the ligand alone from the complex's current positions, and
    /// records the three relaxations.
    fn relax_fragments(
        &mut self,
        free: Vec<bool>,
        complex: Progress,
        max_steps: usize,
        tolerance: f64,
    ) {
        let n0 = self.model.pocket_len;
        let mut pocket = self.at[..n0].to_vec();
        let mut ligand = self.at[n0..].to_vec();
        let p = relax(
            &self.model.pocket,
            &free[..n0],
            &mut pocket,
            max_steps,
            tolerance,
        );
        let l = relax(
            &self.model.ligand,
            &free[n0..],
            &mut ligand,
            max_steps,
            tolerance,
        );
        self.apart = Some(Apart {
            relaxation: HydrogenRelaxation {
                free,
                tolerance,
                max_steps,
                complex,
                pocket: p,
                ligand: l,
            },
            pocket,
            ligand,
        });
    }

    /// ΔE_bind in vacuum at the current positions — each system's own, for a relaxed binding: the
    /// cross terms summed directly, the fragments' reorganisation, and the three-evaluation
    /// difference. See [`Interaction`].
    pub fn interaction(&self) -> Interaction {
        let n0 = self.model.pocket_len;
        let (vdw, elec) = self.cross(&self.model.charges, None);
        let e_c = self.model.complex.energy(&self.at);
        let e_p = self.model.pocket.energy(self.pocket_alone_positions());
        let e_l = self.model.ligand.energy(self.ligand_alone_positions());
        let d = |c: f64, p: f64, l: f64| c - p - l;
        let reorganisation = match self.apart {
            None => Energy::default(),
            Some(_) => {
                let p0 = self.model.pocket.energy(&self.at[..n0]);
                let l0 = self.model.ligand.energy(&self.at[n0..]);
                let r = |p: f64, pr: f64, l: f64, lr: f64| (p - pr) + (l - lr);
                Energy {
                    bond: r(p0.bond, e_p.bond, l0.bond, e_l.bond),
                    angle: r(p0.angle, e_p.angle, l0.angle, e_l.angle),
                    torsion: r(p0.torsion, e_p.torsion, l0.torsion, e_l.torsion),
                    inversion: r(p0.inversion, e_p.inversion, l0.inversion, e_l.inversion),
                    van_der_waals: r(
                        p0.van_der_waals,
                        e_p.van_der_waals,
                        l0.van_der_waals,
                        e_l.van_der_waals,
                    ),
                    electrostatic: r(
                        p0.electrostatic,
                        e_p.electrostatic,
                        l0.electrostatic,
                        e_l.electrostatic,
                    ),
                    solvation: 0.0,
                    total: r(p0.total, e_p.total, l0.total, e_l.total),
                }
            }
        };
        Interaction {
            van_der_waals: vdw,
            electrostatic: elec,
            reorganisation,
            difference: Energy {
                bond: d(e_c.bond, e_p.bond, e_l.bond),
                angle: d(e_c.angle, e_p.angle, e_l.angle),
                torsion: d(e_c.torsion, e_p.torsion, e_l.torsion),
                inversion: d(e_c.inversion, e_p.inversion, e_l.inversion),
                van_der_waals: d(e_c.van_der_waals, e_p.van_der_waals, e_l.van_der_waals),
                electrostatic: d(e_c.electrostatic, e_p.electrostatic, e_l.electrostatic),
                solvation: 0.0,
                total: d(e_c.total, e_p.total, e_l.total),
            },
        }
    }

    /// The force the pocket puts on each ligand atom — `−∇` of the cross terms at the complex's
    /// positions with respect to that atom, which is `−∇ ΔE_bind` for an unrelaxed binding, since
    /// nothing else in ΔE_bind moves with the ligand — newtons, in ligand order.
    pub fn interaction_forces(&self) -> Vec<[f64; 3]> {
        let mut forces = vec![[0.0; 3]; self.at.len() - self.model.pocket_len];
        self.cross(&self.model.charges, Some(&mut forces));
        forces
    }

    /// The Coulomb energy (ε = 1) of every protein–ligand pair at the current positions with
    /// `charges` in place of the binding's own — QEq on the whole complex
    /// ([`Binding::polarised_charges`]), say. Joules per molecule.
    ///
    /// # Panics
    ///
    /// If `charges` is not one per atom.
    pub fn cross_electrostatic(&self, charges: &[f64]) -> f64 {
        assert_eq!(charges.len(), self.at.len(), "one charge per atom");
        self.cross(charges, None).1
    }

    /// The cross sums, van der Waals and Coulomb, and their force on the ligand if asked.
    fn cross(&self, charges: &[f64], mut forces: Option<&mut [[f64; 3]]>) -> (f64, f64) {
        let n0 = self.model.pocket_len;
        let (mut vdw, mut elec) = (0.0, 0.0);
        for p in self.model.in_pocket.pairs() {
            let [i, j] = p.atoms;
            // Pairs have i < j, and the ligand's indices are the largest.
            if i >= n0 || j < n0 {
                continue;
            }
            let (a, b) = (self.at[i], self.at[j]);
            let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            let (e, mut de_dr) = p.at(r);
            vdw += e;
            let q = coulomb(charges[i], charges[j], r);
            elec += q;
            de_dr -= q / r;
            if let Some(f) = forces.as_deref_mut() {
                for k in 0..3 {
                    f[j - n0][k] -= de_dr * d[k] / r;
                }
            }
        }
        (vdw, elec)
    }

    /// QEq on the whole complex at the current positions, at its total formal charge: the charges
    /// with each partner polarised by the other and charge transferred between them, which this
    /// module's charges leave out. For measuring that choice; not used by anything here.
    ///
    /// # Errors
    ///
    /// Any [`QeqError`].
    pub fn polarised_charges(&self) -> Result<Charges, QeqError> {
        Qeq::default().equilibrate(
            &self.model.elements,
            &self.at,
            f64::from(self.model.pocket_charge + self.model.ligand_charge),
        )
    }

    /// The ligand's own energy at its current positions.
    pub fn ligand_energy(&self) -> Energy {
        self.model.ligand.energy(self.ligand_positions())
    }

    /// The solvation part of binding at the current positions — each system's own, for a relaxed
    /// binding: OBC II on the complex, the pocket and the ligand with the binding's charges, and
    /// the solvent-accessible areas by
    /// [`surface_area`] at `points` points per atom, Bondi radii and a 1.4 Å probe. See
    /// [`Desolvation`].
    ///
    /// # Panics
    ///
    /// If `points` is zero.
    pub fn desolvation(&self, points: usize) -> Desolvation {
        self.desolvation_with(&self.model.charges, points)
    }

    /// [`Binding::desolvation`] with `charges` in place of the binding's own — the ligand's set to
    /// zero, say, to separate the pocket's desolvation by the ligand's volume from what the
    /// ligand's charges add.
    ///
    /// # Panics
    ///
    /// If `charges` is not one per atom, or `points` is zero.
    pub fn desolvation_with(&self, charges: &[f64], points: usize) -> Desolvation {
        assert_eq!(charges.len(), self.at.len(), "one charge per atom");
        let n0 = self.model.pocket_len;
        let el = &self.model.elements;
        let q = charges;
        let n = self.at.len();
        // Each system at its own positions: the complex's, and the pocket's and the ligand's alone.
        let systems: [(std::ops::Range<usize>, &[[f64; 3]]); 3] = [
            (0..n, &self.at),
            (0..n0, self.pocket_alone_positions()),
            (n0..n, self.ligand_alone_positions()),
        ];
        let gb = |(range, at): &(std::ops::Range<usize>, &[[f64; 3]])| {
            GeneralizedBorn::new(&el[range.clone()]).energy(&q[range.clone()], at)
        };
        let (complex, pocket, ligand) = (gb(&systems[0]), gb(&systems[1]), gb(&systems[2]));
        let radii: Vec<f64> = el.iter().map(|&e| intrinsic_radius(e)).collect();
        let area = |(range, at): &(std::ops::Range<usize>, &[[f64; 3]])| {
            surface_area(&radii[range.clone()], at, PROBE_RADIUS, points)
        };
        let areas_bound = area(&systems[0]);
        let mut areas_apart = area(&systems[1]);
        areas_apart.extend(area(&systems[2]));
        let buried_area: f64 = areas_bound
            .iter()
            .zip(&areas_apart)
            .map(|(b, a)| b - a)
            .sum();
        Desolvation {
            complex,
            pocket,
            ligand,
            polar: complex - pocket - ligand,
            areas_bound,
            areas_apart,
            buried_area,
            nonpolar: nonpolar_energy(buried_area),
        }
    }

    /// Minimises the ligand in the pocket from where it is now: every protein atom frozen
    /// ([`Minimiser::with_frozen`]), the ligand's own energy plus the interaction, in vacuum at
    /// the fixed charges, until the largest force on a ligand atom is at most `tolerance` newtons
    /// or `max_steps` steps have been taken. Says where it stopped.
    ///
    /// **For a relaxed binding** ([`Binding::relaxing_hydrogens`]) the free hydrogens of the pocket
    /// move with the ligand, every heavy atom of the pocket still frozen, so the complex stays at
    /// a minimum in everything its relaxation freed; then the pocket alone and the ligand alone
    /// are relaxed again from where the complex left them, with the relaxation's own steps and
    /// tolerance, and the binding is still relaxed. The progress returned is the complex's, and is
    /// also [`HydrogenRelaxation::complex`].
    pub fn minimise_ligand(&mut self, max_steps: usize, tolerance: f64) -> Progress {
        let n0 = self.model.pocket_len;
        if let Some(a) = self.apart.take() {
            let r = a.relaxation;
            let free: Vec<bool> = (0..self.at.len()).map(|k| k >= n0 || r.free[k]).collect();
            let p = relax(
                &self.model.complex,
                &free,
                &mut self.at,
                max_steps,
                tolerance,
            );
            self.relax_fragments(r.free, p, r.max_steps, r.tolerance);
            return p;
        }
        let frozen: Vec<bool> = (0..self.at.len()).map(|k| k < n0).collect();
        let mut m = Minimiser::new(tolerance).with_frozen(frozen);
        let ff = &self.model.in_pocket;
        let mut p = m.step(ff, &[], &mut self.at);
        while p.status == Status::Running && p.steps < max_steps {
            p = m.step(ff, &[], &mut self.at);
        }
        p
    }
}

/// Minimises `ff` over the atoms `free` marks, every other atom frozen, from `at` in place. The
/// force field minimised is `ff` less its terms among frozen atoms alone — a constant while they
/// do not move — so the force on a free atom is `ff`'s bit for bit, and [`Progress::energy`] is
/// the reduced force field's.
fn relax(
    ff: &ForceField,
    free: &[bool],
    at: &mut [[f64; 3]],
    max_steps: usize,
    tolerance: f64,
) -> Progress {
    let reduced = ff.touching(free);
    let frozen: Vec<bool> = free.iter().map(|f| !f).collect();
    let mut m = Minimiser::new(tolerance).with_frozen(frozen);
    let mut p = m.step(&reduced, &[], at);
    while p.status == Status::Running && p.steps < max_steps {
        p = m.step(&reduced, &[], at);
    }
    p
}

fn centroid(at: &[[f64; 3]]) -> [f64; 3] {
    let mut c = [0.0; 3];
    for p in at {
        for (ck, pk) in c.iter_mut().zip(p) {
            *ck += pk;
        }
    }
    let n = at.len().max(1) as f64;
    [c[0] / n, c[1] / n, c[2] / n]
}
