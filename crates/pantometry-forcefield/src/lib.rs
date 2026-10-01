//! pantometry-forcefield: a small molecule, read from the chemical dictionary and typed for a
//! force field.
//!
//! The atomistic counterpart to `pantometry-protein`'s coarse network: one body per atom,
//! hydrogens included, with the bonds the chemistry has rather than the ones a cutoff implies.
//! The aim, over the steps that follow this one, is molecular mechanics on drug-sized molecules —
//! an energy, its minimum, and which conformations are stable — with every term checkable.
//!
//! **This step computes UFF's whole energy, but nothing moves yet.** What is here:
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
//! - [`Molecule`] is the kernel [`Domain`]: the atoms as [`Bodies`] with their names and bonds,
//!   so the scene layer draws a ball-and-stick molecule without knowing what a molecule is, and
//!   the energy terms as readings.
//!
//! # What it is checked against
//!
//! Facts about the molecule and the paper that this code did not produce. Aspirin (`AIN`) parses
//! to the atom count and composition its own `_chem_comp.formula` states; every hydrogen has one
//! bond; every heavy atom's bond orders sum to its valence, with an aromatic bond counted as one
//! and a half; the types of aspirin's thirteen heavy atoms are the ones a chemist reads off the
//! structure; and the natural angles of `C_3`, `C_R` and `C_1` are the tetrahedral, trigonal and
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
//! # What is deliberately not in it
//!
//! - **No motion and no minimisation**: [`Molecule::step`](Domain::step) leaves the atoms where
//!   the file put them, so nothing here has yet been compared with a structure or barrier the
//!   paper *minimised* — its Table II barriers are relaxed, and a rigid rotation is not.
//! - **The electronegativities are not from their primary source.** χ is transcribed from Open
//!   Babel, which copies RDKit; the paper it comes from has not been read. [`uff`] says so where
//!   the values are, and what the one partial check pins.
//! - **No ring perception, no bond-order guessing, no protonation.** The dictionary states
//!   aromaticity, bond orders and explicit hydrogens, and a second opinion computed here could
//!   only disagree with it. A file without bonds is refused, not guessed at.
//! - **No elements beyond H, C, N, O, F, P, S, Cl, Br and I**, and no hypervalent sulfur types —
//!   [`uff::assign`] documents the sulfone gap.
//! - **No charge model, no solvent, no binding energies.** Electrostatics is computed from
//!   charges a caller supplies, and they default to zero — the paper's charges come from charge
//!   equilibration, which is the next step. A hypervalent sulfur is refused by name
//!   ([`Unsupported`]) rather than given the wrong radius.
//!
//! Nothing here opens a file: every crate in this workspace compiles to `wasm32`, so a caller
//! reads the text and passes it in.

#![deny(missing_docs)]

pub mod angular;
pub mod ccd;
pub mod energy;
pub mod uff;

pub use angular::{Bend, Inversion, Torsion};
pub use ccd::{Atom, Bond, BondOrder, CcdError, Component, Coordinates, Element};
pub use energy::{Energy, Evaluation, ForceField, Unsupported};
pub use uff::{Parameters, TableI, UffType};

use pantometry_core::{Bodies, Domain, Exchange, Kind, Ledger, Reading, Violation};
use pantometry_units::{LengthVec, Qty, Time};

/// A molecule as a domain: atoms at places, typed, bonded, and drawable.
///
/// # How it moves
///
/// **It does not, yet.** [`Domain::step`] leaves every atom where the dictionary put it, so
/// [`Domain::max_stable_dt`] is infinite — honestly, because nothing moves. It has an energy now
/// — [`Molecule::evaluate`], and the readings — but no dynamics and no minimiser to use it, and
/// the ledger stays empty: an energy that never changes is not a flow anything balances against.
/// Until minimisation arrives a run of this domain is a picture of a parsed and typed molecule,
/// with its energy beside it.
///
/// # Readings
///
/// `atoms`, `heavy atoms`, `bonds` and `formal charge`, and the energy at the current positions in
/// **kcal/mol** — the paper's unit and the one every reader of a force field compares in — as
/// `energy`, `bond stretch`, `angle bend`, `torsion`, `inversion`, `van der Waals` and
/// `electrostatic`. `energy` is UFF's total, the sum of the six. For a molecule
/// [`ForceField::new`] refuses, the seven energy readings are `NaN` and [`Molecule::force_field`]
/// says why.
#[derive(Clone, Debug)]
pub struct Molecule {
    name: String,
    component: Component,
    types: Vec<UffType>,
    force_field: Result<ForceField, Unsupported>,
    at: Vec<[f64; 3]>,
    saved: Option<Vec<[f64; 3]>>,
}

impl Molecule {
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
            saved: None,
        }
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

    /// Infinite, and true for now: the atoms do not move, so no step is too long. This becomes a
    /// real limit — set by the stiffest bond's period — when the energy terms arrive.
    fn max_stable_dt(&self, _now: Time) -> Time {
        Qty::from_si(f64::INFINITY)
    }

    /// Nothing moves. See [`Molecule`].
    fn step(&mut self, _t: Time, _dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        Ok(())
    }

    /// Empty: nothing moves, so nothing flows, and an energy that cannot change has nothing to
    /// balance against. This becomes an entry when the molecule can move.
    fn ledger(&self) -> Ledger {
        Ledger::new()
    }

    /// True and exact: it holds nothing and takes nothing from the bus.
    fn books_balance(&self) -> bool {
        true
    }

    fn checkpoint(&mut self) {
        self.saved = Some(self.at.clone());
    }

    fn restore(&mut self) {
        if let Some(saved) = &self.saved {
            self.at.clone_from(saved);
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    fn readings(&self) -> Vec<Reading> {
        let atoms = self.component.atoms();
        let heavy = atoms.iter().filter(|a| a.element != Element::H).count();
        let charge: i32 = atoms.iter().map(|a| a.charge).sum();
        let e = self.evaluate().map(|ev| ev.energy).unwrap_or(Energy::NAN);
        let kcal = |joules: f64| joules / uff::KCAL_PER_MOL;
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
        ]
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
