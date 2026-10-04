//! pantometry-forcefield: a small molecule, read from the chemical dictionary and typed for a
//! force field.
//!
//! The atomistic counterpart to `pantometry-protein`'s coarse network: one body per atom,
//! hydrogens included, with the bonds the chemistry has rather than the ones a cutoff implies.
//! The aim, over the steps that follow this one, is molecular mechanics on drug-sized molecules —
//! an energy, its minimum, and which conformations are stable — with every term checkable.
//!
//! **This step computes UFF's whole energy and minimises it.** What is here:
//!
//! - [`Component::from_ccd`] reads one entry of the wwPDB Chemical Component Dictionary, strictly:
//!   every malformed or unsupported input is a [`CcdError`] naming what it refused. See [`ccd`].
//! - [`uff::assign`] gives every atom its Universal Force Field type from its element, its bond
//!   orders and the dictionary's aromatic flags, and [`UffType::parameters`] gives each type's
//!   Table I row in SI. See [`uff`] for the citation, the table and the rules.
//! - [`ForceField`] holds all six of UFF's terms — harmonic bond stretch, angle bend, torsion,
//!   inversion, Lennard-Jones van der Waals, and Coulomb electrostatics with 1-2 and 1-3
//!   exclusions — and gives the energy, term by term, and the analytic force on every atom. See
//!   [`energy`] and [`angular`] for the equations, the page each is on, and the places this
//!   departs from the paper as printed: the sign of `r_EN`, the sign of the linear angle term,
//!   and a torsion switched off continuously rather than abruptly near a straight angle.
//! - [`ForceField::minimise`] relaxes a geometry: L-BFGS with an Armijo line search, so the
//!   energy never rises, to a stated force. [`DihedralRestraint`] holds a dihedral while the
//!   rest relaxes, which is how the paper's torsion barriers are compared. See [`minimise`].
//! - [`qeq`] computes charge-equilibration charges (Rappé and Goddard 1991), the ones UFF
//!   prescribes, at a stated geometry: the exact Slater Coulomb integrals, hydrogen's
//!   charge-dependent orbital iterated as the paper does, and each misprint and ambiguity in the
//!   paper settled by reproducing its tables. [`ForceField::with_qeq_charges`] puts them in the
//!   electrostatic term, fixed at that geometry; zero charges stay the default.
//! - [`solvation`] is implicit water: the generalized Born electrostatic solvation energy of
//!   Onufriev, Bashford and Case (2004), OBC II, with ε = 80 and optional salt, and its analytic
//!   force through the Born radii at fixed charges. [`ForceField::with_generalized_born`] adds it
//!   to the energy and the force; vacuum stays the default. A Shrake–Rupley surface gives the
//!   paper's nonpolar term beside it, as a number, not a force.
//! - [`System::from_pdb`] reads a protein chain and its ligand out of a PDB entry and makes them one
//!   [`Component`] the rest of this crate takes unchanged: hydrogens placed by superposing each
//!   residue's dictionary template on the crystal, bonds with orders (peptide bonds and disulfides
//!   included), formal charges at pH 7 by stated rules, and every atom marked [`Part::Protein`] or
//!   [`Part::Ligand`]. Strict, as the dictionary reader is: a missing heavy atom, a gap in a chain,
//!   a residue that is not the sequence's or an unnamed `HETATM` residue is a [`PdbError`] that
//!   names it. See [`pdb`].
//! - [`Molecule`] is the kernel [`Domain`]: the atoms as [`Bodies`] with their names and bonds,
//!   so the scene layer draws a ball-and-stick molecule without knowing what a molecule is, the
//!   energy terms and the force as readings, and **one minimiser iteration per step**, so a run
//!   is the molecule relaxing, frame by frame.
//!
//! # What it is checked against
//!
//! Facts about the molecule and the paper that this code did not produce. Aspirin (`AIN`) parses
//! to the atom count and composition its own `_chem_comp.formula` states; every hydrogen has one
//! bond; every heavy atom's bond orders sum to its valence, with an aromatic bond counted as one
//! and a half; the types of aspirin's thirteen heavy atoms are the ones a chemist reads off the
//! structure — its acid and ester oxygens resonant, as the paper types a conjugated ether oxygen
//! (see [`uff::assign`]); and the natural angles of `C_3`, `C_R` and `C_1` are the tetrahedral, trigonal and
//! linear angles of geometry. Every refusal has a test that feeds it the input it refuses.
//!
//! The energy terms against closed forms: a bond's energy is zero at its natural length and
//! `½ k Δr²` away from it; `C_3–C_3` is 1.514 Å; the paper's own worked Si–O correction is
//! reproduced; the van der Waals minimum is `−D_IJ` at `x_IJ` and crosses zero at `x_IJ / 2^(1/6)`;
//! two unit charges at 1 Å are 332.0637 kcal/mol; the exclusions leave exactly the 1-4 pair of a
//! four-atom chain. Every angle form is zero with zero slope at θ₀ and has curvature `K` there,
//! with `K` worked by hand for methane; ethane's nine torsions sum to a barrier of exactly
//! `V_C_3` = 2.119 kcal/mol; eq 17 by hand for ethylene and benzene; the inversion is zero when
//! planar and matches a hand value when not. The forces against central finite differences of
//! the energy, to a tolerance derived from the step; the energy against translation and rotation;
//! and the net force against zero.
//!
//! QEq against closed forms — every `ns` self-repulsion as an exact fraction, Roothaan's 1s–1s
//! `J(R)`, `J → 1/R`, very unequal exponents against independent references, the paper's
//! two-atom eq 18, both ends of the range clamping, the total charge, rigid-motion invariance —
//! and against its own paper: Table II's eighteen alkali halides, Table III's four hydrogens and
//! 36 Table IV charges at stated experimental geometries, each to its printed precision plus its
//! measured sensitivity to the geometry. See [`qeq`].
//!
//! Generalized Born against closed forms: eq 5's pair integral against its own quadrature, in
//! every branch; a lone ion's radius `ρ − 0.09 Å` and Born energy; two ions far apart as two Born
//! terms and a screened Coulomb, to a bound derived from the descreening; two at one point as one
//! of their summed charge; a charge inside a dielectric sphere, where HCT gives the Coulomb-field
//! radius exactly and Kirkwood's series shows it 5.73% long; eq 6 with eq 8's constants; the
//! paper's 30 Å bound; and the force against central differences. See [`solvation`].
//!
//! A protein against its own entry: T4 lysozyme L99A with benzene (PDB 181L) has every heavy atom
//! its templates name at the file's coordinates bit for bit, its sequence is `SEQRES`'s with the
//! two residues `REMARK 465` lists as unlocated, each residue's formula after protonation is the
//! textbook residue's adjusted by the rule, every heavy atom's bond orders meet its valence with its
//! formal charge, and the total, +9, is `Arg + Lys − Asp − Glu` counted from `SEQRES` and the
//! salt-bridged His31 the rule makes +1. Every
//! placed hydrogen keeps its template's bond length exactly and its angles within a bound its own
//! fit earns, no hydrogen is within 1.5 Å of a heavy atom, every superposition is no worse than one
//! built by hand from three atoms, and every refusal is fed the input it refuses.
//!
//! The minimiser against geometry whose minimum is known exactly: a diatomic relaxes to eq 2's
//! natural length, water to its two natural lengths and θ₀ (no non-bonded pair is left in
//! either), and methane to the regular tetrahedron — each to within the displacement its final
//! force allows, `|F| / λ_min` with λ_min the measured smallest Hessian eigenvalue. Then against
//! the paper's own minimised numbers: the relaxed torsion barriers of its Table II and the
//! structures of its Figures 4–7, compared row by row, and asserted only where the comparison is
//! earned.
//!
//! # Determinism, and where it stops
//!
//! No clock, no randomness, no hash order and no threads, so a run repeats bit for bit on one
//! machine, at any optimisation level. **Not across platforms.** IEEE-754 defines `+ − × ÷` and
//! `sqrt` exactly, but this crate also calls the platform's `cos`, `sin`, `ln`, `asin` and, under a
//! restraint, `atan2` — on every torsion evaluation as well as when a [`ForceField`] is built.
//! Those are not correctly rounded and are not required to be the same function on two
//! machines, so no digest of a minimum is pinned; see [`minimise`]. [`qeq`]'s integrals call
//! `exp`, with the same consequence for a charge, and [`solvation`] calls `exp`, `ln` and `tanh`.
//! [`pdb`]'s superposition is arithmetic and `sqrt`, but turning a rotor calls `sin` and `cos`.
//!
//! # What is deliberately not in it
//!
//! - **No dynamics.** A step of [`Molecule`] is a minimiser iteration, not a time step: there are
//!   no velocities, no temperature, and `dt` is not used. Molecular dynamics is a later step.
//! - **No conformer search.** The minimiser finds the minimum downhill from where it starts; the
//!   torsion scans hold a dihedral to find a barrier, and nothing looks for the global minimum.
//! - **No reading of the paper but this crate's in use**: [`Variant`] keeps one other — `r_EN`
//!   added as eq 2 prints it, or dropped — only as the evidence that settled the sign, which the
//!   paper's own minimised structures did (see [`energy`]). The other question 1d left open, the
//!   group-6 sp³–sp² torsion's minimum, was settled by typing a conjugated O or S resonant as the
//!   paper does, after which the row it asked about is reached only by an oxonium oxygen.
//! - **No ring perception, no bond-order guessing, no protonation of a dictionary entry.** The
//!   dictionary states aromaticity, bond orders and explicit hydrogens, and a second opinion
//!   computed here could only disagree with it. A file without bonds is refused, not guessed at.
//!   The one place protonation is decided is [`pdb`], for the twenty amino acids at pH 7, by a
//!   table of stated rules applied to the dictionary's templates — no pKa is computed, a histidine
//!   in a salt bridge is +1 by a stated geometric rule and any other is neutral, and no
//!   hydrogen-bond network
//!   is optimised.
//! - **No elements beyond H, C, N, O, F, P, S, Cl, Br and I**, and no hypervalent sulfur types —
//!   [`uff::assign`] documents the sulfone gap.
//! - **No charges unless asked for, no solvent unless asked for, no binding energies.**
//!   Electrostatics is computed
//!   from charges a caller supplies, and they default to zero, because UFF's valence parameters
//!   were fitted without them (p. 10031 of the UFF paper). QEq's are an option
//!   ([`ForceField::with_qeq_charges`]), **fixed at the geometry they were computed for**: the
//!   force treats them as constants, as the QEq paper uses them, so a minimisation with them does
//!   not re-equilibrate as it moves. Not here: charges that follow the geometry, and QEq's
//!   polarisation extensions. A hypervalent sulfur is refused by name ([`Unsupported`]) rather
//!   than given the wrong radius.
//! - **Implicit solvent is generalized Born only, its electrostatic part only in the force
//!   field, and opt-in** ([`ForceField::with_generalized_born`]). The nonpolar `0.005 × SASA` term
//!   is computed ([`solvation::surface_area`], [`solvation::nonpolar_energy`]) but not added to
//!   the energy a minimiser sees, because a point-counted surface has no gradient worth the name.
//!   Not here: Poisson–Boltzmann, explicit water, GB with a solute dielectric other than 1, a
//!   cutoff, and the binding energy itself (complex − protein − ligand), which is the next step.
//!   The radii and scale factors are AMBER's as OpenMM holds them, not read in primary, and Br and
//!   I are outside the set OBC was fitted with; see [`solvation`]. **With QEq charges, OBC II
//!   over-solvates** small molecules against experiment — mean −2.1, RMS 3.4 kcal/mol, ethers by
//!   about −4 to −6 — so absolute solvation and binding energies built on it are not
//!   quantitative.
//!
//! Nothing here opens a file: every crate in this workspace compiles to `wasm32`, so a caller
//! reads the text and passes it in.

#![deny(missing_docs)]

pub mod angular;
pub mod ccd;
pub mod energy;
pub mod minimise;
pub mod pdb;
pub mod qeq;
pub mod solvation;
pub mod uff;

pub use angular::{Bend, Inversion, Torsion};
pub use ccd::{Atom, Bond, BondOrder, CcdError, Component, Coordinates, Element};
pub use energy::{Energy, Evaluation, ForceField, Unsupported, Variant};
pub use minimise::{DihedralRestraint, Minimiser, Progress, Status};
pub use pdb::{Histidine, Part, PdbError, Placement, Residue, Selection, System};
pub use qeq::{Charges, Qeq, QeqError};
pub use solvation::{GeneralizedBorn, Rescaling};
pub use uff::{Parameters, TableI, UffType};

use pantometry_core::{Bodies, Domain, Exchange, Kind, Ledger, Reading, Violation};
use pantometry_units::{LengthVec, Qty, Time};

/// A molecule as a domain: atoms at places, typed, bonded, and drawable.
///
/// # How it moves
///
/// **Downhill, one minimiser iteration per step.** [`Domain::step`] takes one L-BFGS step of
/// [`Minimiser`] on the positions, so each frame of a run is the molecule further relaxed and
/// the `energy` reading never rises from one frame to the next. It stops moving when the largest
/// force on any atom is at most the tolerance ([`Molecule::DEFAULT_TOLERANCE`], or
/// [`Molecule::with_tolerance`]) or when no step lowers the energy, and says which. The step's
/// `dt` is not used: this is a minimisation, not dynamics, so [`Domain::max_stable_dt`] is
/// infinite. The ledger stays empty — the energy that leaves is not a flow to anywhere, and
/// claiming one would be inventing a sink.
///
/// [`Molecule::minimise`] runs the same iterations to convergence in one call.
///
/// # Readings
///
/// `atoms`, `heavy atoms`, `bonds` and `formal charge`, and the energy at the current positions in
/// **kcal/mol** — the paper's unit and the one every reader of a force field compares in — as
/// `energy`, `bond stretch`, `angle bend`, `torsion`, `inversion`, `van der Waals` and
/// `electrostatic`. `energy` is UFF's total, the sum of the six. `max force` and `rms force`
/// (kcal mol⁻¹ Å⁻¹, per-atom force vectors), `converged` (1 or 0) and `minimiser steps`. For a
/// molecule [`ForceField::new`] refuses, the energy and force readings are `NaN` and
/// [`Molecule::force_field`] says why. The last four are its [`Domain::diagnostics`]: they
/// describe the minimisation, not the molecule.
#[derive(Clone, Debug)]
pub struct Molecule {
    name: String,
    component: Component,
    types: Vec<UffType>,
    force_field: Result<ForceField, Unsupported>,
    at: Vec<[f64; 3]>,
    minimiser: Minimiser,
    saved: Option<(Vec<[f64; 3]>, Minimiser)>,
}

impl Molecule {
    /// The tolerance a new molecule minimises to: the largest force on any atom at most
    /// 1e-4 kcal mol⁻¹ Å⁻¹, in newtons. A bond of `k` ~ 700 kcal mol⁻¹ Å⁻² is then within 1.4e-7 Å
    /// of where the force would put it.
    pub const DEFAULT_TOLERANCE: f64 = 1e-4 * minimise::KCAL_PER_MOL_ANGSTROM;

    /// A molecule named `name`, at the component's coordinates, typed by [`uff::assign`], with
    /// its [`ForceField`] built and every partial charge zero.
    ///
    /// Never fails: a molecule the force field refuses is still read, typed and drawn, and
    /// [`Molecule::force_field`] carries the refusal.
    pub fn new(name: impl Into<String>, component: Component) -> Molecule {
        let types = uff::assign(&component);
        let force_field = ForceField::new(&component, &types);
        let at = component.atoms().iter().map(|a| a.at).collect();
        Molecule {
            name: name.into(),
            component,
            types,
            force_field,
            at,
            minimiser: Minimiser::new(Molecule::DEFAULT_TOLERANCE),
            saved: None,
        }
    }

    /// The same molecule minimising to `tolerance` newtons instead.
    pub fn with_tolerance(mut self, tolerance: f64) -> Molecule {
        self.minimiser = Minimiser::new(tolerance);
        self
    }

    /// Runs minimiser iterations until converged, stalled, or `max_steps` more have been taken,
    /// and says where it stopped; `None` for a molecule the force field refuses.
    pub fn minimise(&mut self, max_steps: usize) -> Option<Progress> {
        let ff = self.force_field.as_ref().ok()?;
        let mut last = None;
        for _ in 0..max_steps.max(1) {
            let p = self.minimiser.step(ff, &[], &mut self.at);
            last = Some(p);
            if p.status != Status::Running {
                break;
            }
        }
        last
    }

    /// The minimiser's status after the last step.
    pub fn status(&self) -> Status {
        self.minimiser.status()
    }

    /// The same molecule with partial charges `charges`, in elementary charges, one per atom.
    ///
    /// # Panics
    ///
    /// If `charges` is not one per atom.
    pub fn with_charges(mut self, charges: Vec<f64>) -> Molecule {
        assert_eq!(charges.len(), self.at.len(), "one charge per atom");
        if let Ok(ff) = self.force_field {
            self.force_field = Ok(ff.with_charges(charges));
        }
        // A different energy: the minimiser's history and status belong to the old one.
        self.minimiser = Minimiser::new(self.minimiser.tolerance());
        self
    }

    /// Its force field, or why there is none.
    pub fn force_field(&self) -> Result<&ForceField, &Unsupported> {
        self.force_field.as_ref()
    }

    /// The energy and forces at the current positions, or why there are none.
    pub fn evaluate(&self) -> Result<Evaluation, &Unsupported> {
        self.force_field().map(|ff| ff.evaluate(&self.at))
    }

    /// The component it was built from.
    pub fn component(&self) -> &Component {
        &self.component
    }

    /// Every atom's UFF type, in the component's atom order.
    pub fn types(&self) -> &[UffType] {
        &self.types
    }

    /// Where atom `i` is now, in metres.
    ///
    /// # Panics
    ///
    /// If `i` is at or past the atom count.
    pub fn at(&self, i: usize) -> [f64; 3] {
        self.at[i]
    }
}

impl Bodies for Molecule {
    fn count(&self) -> usize {
        self.at.len()
    }

    fn position(&self, i: usize) -> LengthVec {
        LengthVec::new(
            Qty::from_si(self.at[i][0]),
            Qty::from_si(self.at[i][1]),
            Qty::from_si(self.at[i][2]),
        )
    }

    /// The atomic number, so a picture is coloured by element — the one per-atom quantity every
    /// reader of a ball-and-stick model already knows how to read.
    fn value(&self, i: usize) -> f64 {
        f64::from(self.component.atoms()[i].element.atomic_number())
    }

    /// Dimensionless: an atomic number is a count.
    fn value_unit(&self) -> &'static str {
        ""
    }

    /// The atom's name in the dictionary, `C7`, `HO1` — what the entry and every bond in it call
    /// the atom, and so what a reader can match against the file.
    fn label(&self, i: usize) -> Option<String> {
        self.component.atoms().get(i).map(|a| a.name.clone())
    }

    /// Every bond the dictionary states, as index pairs. Real by construction: the reader refuses
    /// a bond that names an atom the entry does not have.
    fn bonds(&self) -> Vec<[u32; 2]> {
        self.component
            .bonds()
            .iter()
            .map(|b| [b.atoms[0] as u32, b.atoms[1] as u32])
            .collect()
    }
}

impl Domain for Molecule {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// Infinite: a step is a minimiser iteration, and its length in time means nothing. See
    /// [`Molecule`].
    fn max_stable_dt(&self, _now: Time) -> Time {
        Qty::from_si(f64::INFINITY)
    }

    /// One minimiser iteration; nothing, for a molecule the force field refuses. See
    /// [`Molecule`].
    fn step(&mut self, _t: Time, _dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        if let Ok(ff) = &self.force_field {
            self.minimiser.step(ff, &[], &mut self.at);
        }
        Ok(())
    }

    /// Empty: a minimisation lowers the energy without sending it anywhere, so there is no flow
    /// to record. See [`Molecule`].
    fn ledger(&self) -> Ledger {
        Ledger::new()
    }

    /// True and exact: it holds nothing and takes nothing from the bus.
    fn books_balance(&self) -> bool {
        true
    }

    /// The positions and the minimiser's whole state, so a restored run repeats bit for bit.
    fn checkpoint(&mut self) {
        self.saved = Some((self.at.clone(), self.minimiser.clone()));
    }

    fn restore(&mut self) {
        if let Some((at, minimiser)) = &self.saved {
            self.at.clone_from(at);
            self.minimiser = minimiser.clone();
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    fn readings(&self) -> Vec<Reading> {
        let atoms = self.component.atoms();
        let heavy = atoms.iter().filter(|a| a.element != Element::H).count();
        let charge: i32 = atoms.iter().map(|a| a.charge).sum();
        let ev = self.evaluate().ok();
        let e = ev.as_ref().map(|ev| ev.energy).unwrap_or(Energy::NAN);
        let kcal = |joules: f64| joules / uff::KCAL_PER_MOL;
        let (max_force, rms_force) = ev.as_ref().map_or((f64::NAN, f64::NAN), |ev| {
            let sq: Vec<f64> = ev
                .forces
                .iter()
                .map(|f| f[0] * f[0] + f[1] * f[1] + f[2] * f[2])
                .collect();
            let max = sq.iter().fold(0.0f64, |m, v| m.max(*v)).sqrt();
            let rms = (sq.iter().sum::<f64>() / sq.len().max(1) as f64).sqrt();
            (
                max / minimise::KCAL_PER_MOL_ANGSTROM,
                rms / minimise::KCAL_PER_MOL_ANGSTROM,
            )
        });
        let converged = if self.minimiser.status() == Status::Converged {
            1.0
        } else {
            0.0
        };
        vec![
            Reading::new(&self.name, "atoms", atoms.len() as f64, ""),
            Reading::new(&self.name, "heavy atoms", heavy as f64, ""),
            Reading::new(&self.name, "bonds", self.component.bonds().len() as f64, ""),
            Reading::new(&self.name, "formal charge", f64::from(charge), "e"),
            Reading::new(&self.name, "energy", kcal(e.total), "kcal/mol"),
            Reading::new(&self.name, "bond stretch", kcal(e.bond), "kcal/mol"),
            Reading::new(&self.name, "angle bend", kcal(e.angle), "kcal/mol"),
            Reading::new(&self.name, "torsion", kcal(e.torsion), "kcal/mol"),
            Reading::new(&self.name, "inversion", kcal(e.inversion), "kcal/mol"),
            Reading::new(
                &self.name,
                "van der Waals",
                kcal(e.van_der_waals),
                "kcal/mol",
            ),
            Reading::new(
                &self.name,
                "electrostatic",
                kcal(e.electrostatic),
                "kcal/mol",
            ),
            Reading::new(&self.name, "max force", max_force, "kcal/mol/Å"),
            Reading::new(&self.name, "rms force", rms_force, "kcal/mol/Å"),
            Reading::new(&self.name, "converged", converged, ""),
            Reading::new(
                &self.name,
                "minimiser steps",
                self.minimiser.steps() as f64,
                "",
            ),
        ]
    }

    /// The four that describe the minimisation rather than the molecule: `max force` and
    /// `rms force` are its gradient, the residual a minimiser drives towards zero, and
    /// `converged` and `minimiser steps` are its state. A sweep that compared them across runs as
    /// though they converged to something would be measuring the stopping rule.
    fn diagnostics(&self) -> &'static [&'static str] {
        &["max force", "rms force", "converged", "minimiser steps"]
    }

    fn as_bodies(&self) -> Option<&dyn Bodies> {
        Some(self)
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }
}
