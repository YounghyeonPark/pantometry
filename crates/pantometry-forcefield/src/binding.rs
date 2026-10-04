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
//! # What it is not
//!
//! A binding *energy* at one geometry, not a free energy: no entropy, no protein flexibility, no
//! ligand strain against its own solution minimum, no water in the pocket but GB's continuum.
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

/// ΔE_bind in vacuum at one geometry, joules per molecule.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Interaction {
    /// The van der Waals energy of every protein–ligand pair, summed directly.
    pub van_der_waals: f64,
    /// The Coulomb energy (ε = 1) of every protein–ligand pair, summed directly.
    pub electrostatic: f64,
    /// `E(complex) − E(pocket) − E(ligand)`, term by term, from three whole evaluations: the
    /// definition, and in exact arithmetic `van_der_waals` and `electrostatic` in those two fields
    /// and zero in the other four. In floating point it carries the rounding of the pocket's whole
    /// energy, which the direct sums do not.
    pub difference: Energy,
}

impl Interaction {
    /// `van_der_waals + electrostatic`: ΔE_bind from the cross terms.
    pub fn total(&self) -> f64 {
        self.van_der_waals + self.electrostatic
    }
}

/// The solvation part of binding at one geometry: OBC II and the nonpolar surface term, each as
/// complex − pocket − ligand. Energies in joules per molecule, areas in m².
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
    /// Each atom's area in its own fragment, complex order.
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
}

/// A ligand in a rigid pocket: see the module documentation.
///
/// Atoms are in **complex order**: the pocket's atoms first, in the system's order, then the
/// ligand's. [`Binding::system_atoms`] maps each back to the [`System`].
#[derive(Clone, Debug)]
pub struct Binding {
    model: Arc<Model>,
    at: Vec<[f64; 3]>,
}

impl Binding {
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
        Ok(Binding {
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
    pub fn ligand_at(&self, motion: &RigidMotion) -> Binding {
        let c = self.ligand_centroid();
        let mut at = self.at.clone();
        for p in &mut at[self.model.pocket_len..] {
            *p = motion.apply(c, *p);
        }
        Binding {
            model: Arc::clone(&self.model),
            at,
        }
    }

    /// The same binding with the ligand's atoms at `at` (metres, ligand order) — a conformer, or
    /// one atom displaced. The pocket and every charge unchanged.
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
        }
    }

    /// The same binding with **every** atom, protein and ligand, moved by `motion` about the
    /// ligand's centroid: the same complex in another frame, whose ΔE_bind is the same — which the
    /// tests check.
    pub fn moved(&self, motion: &RigidMotion) -> Binding {
        let c = self.ligand_centroid();
        Binding {
            model: Arc::clone(&self.model),
            at: self.at.iter().map(|&p| motion.apply(c, p)).collect(),
        }
    }

    /// ΔE_bind in vacuum at the current positions: the cross terms summed directly, and the
    /// three-evaluation difference. See [`Interaction`].
    pub fn interaction(&self) -> Interaction {
        let n0 = self.model.pocket_len;
        let (vdw, elec) = self.cross(&self.model.charges, None);
        let e_c = self.model.complex.energy(&self.at);
        let e_p = self.model.pocket.energy(&self.at[..n0]);
        let e_l = self.model.ligand.energy(&self.at[n0..]);
        let d = |c: f64, p: f64, l: f64| c - p - l;
        Interaction {
            van_der_waals: vdw,
            electrostatic: elec,
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

    /// The force the pocket puts on each ligand atom — `−∇` of the cross terms with respect to
    /// that atom, which is `−∇ ΔE_bind` since nothing else in ΔE_bind moves with the ligand —
    /// newtons, in ligand order.
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

    /// The solvation part of binding at the current positions: OBC II on the complex, the pocket
    /// and the ligand with the binding's charges, and the solvent-accessible areas by
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
        let gb = |range: std::ops::Range<usize>| {
            GeneralizedBorn::new(&el[range.clone()]).energy(&q[range.clone()], &self.at[range])
        };
        let n = self.at.len();
        let (complex, pocket, ligand) = (gb(0..n), gb(0..n0), gb(n0..n));
        let radii: Vec<f64> = el.iter().map(|&e| intrinsic_radius(e)).collect();
        let area = |range: std::ops::Range<usize>| {
            surface_area(&radii[range.clone()], &self.at[range], PROBE_RADIUS, points)
        };
        let areas_bound = area(0..n);
        let mut areas_apart = area(0..n0);
        areas_apart.extend(area(n0..n));
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
    pub fn minimise_ligand(&mut self, max_steps: usize, tolerance: f64) -> Progress {
        let n0 = self.model.pocket_len;
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
