//! pantometry-forcefield: a small molecule, read from the chemical dictionary and typed for a
//! force field.
//!
//! The atomistic counterpart to `pantometry-protein`'s coarse network: one body per atom,
//! hydrogens included, with the bonds the chemistry has rather than the ones a cutoff implies.
//! The aim, over the steps that follow this one, is molecular mechanics on drug-sized molecules —
//! an energy, its minimum, and which conformations are stable — with every term checkable.
//!
//! **This is the first step and it computes no energy.** What is here:
//!
//! - [`Component::from_ccd`] reads one entry of the wwPDB Chemical Component Dictionary, strictly:
//!   every malformed or unsupported input is a [`CcdError`] naming what it refused. See [`ccd`].
//! - [`uff::assign`] gives every atom its Universal Force Field type from its element, its bond
//!   orders and the dictionary's aromatic flags, and [`UffType::parameters`] gives each type's
//!   Table I row in SI. See [`uff`] for the citation, the table and the rules.
//! - [`Molecule`] is the kernel [`Domain`]: the atoms as [`Bodies`] with their names and bonds,
//!   so the scene layer draws a ball-and-stick molecule without knowing what a molecule is.
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
//! # What is deliberately not in it
//!
//! - **No energy, no forces, no motion — yet.** [`Molecule::step`](Domain::step) leaves the atoms
//!   where the file put them. Bond stretch, angle bend, torsion, inversion and van der Waals terms
//!   are the next step, and each will arrive with a closed form it is checked against.
//! - **No electronegativity.** UFF's natural bond length needs the GMP electronegativity χ, whose
//!   values come from a second paper not yet verified against its source. [`uff`] says so where
//!   the parameters are.
//! - **No ring perception, no bond-order guessing, no protonation.** The dictionary states
//!   aromaticity, bond orders and explicit hydrogens, and a second opinion computed here could
//!   only disagree with it. A file without bonds is refused, not guessed at.
//! - **No elements beyond H, C, N, O, F, P, S, Cl, Br and I**, and no hypervalent sulfur types —
//!   [`uff::assign`] documents the sulfone gap.
//! - **No solvent, no charges for electrostatics, no binding energies.** All of those are later,
//!   and a number that looked like any of them now would be invented.
//!
//! Nothing here opens a file: every crate in this workspace compiles to `wasm32`, so a caller
//! reads the text and passes it in.

#![deny(missing_docs)]

pub mod ccd;
pub mod uff;

pub use ccd::{Atom, Bond, BondOrder, CcdError, Component, Coordinates, Element};
pub use uff::{Parameters, TableI, UffType};

use pantometry_core::{Bodies, Domain, Exchange, Kind, Ledger, Reading, Violation};
use pantometry_units::{LengthVec, Qty, Time};

/// A molecule as a domain: atoms at places, typed, bonded, and drawable.
///
/// # How it moves
///
/// **It does not, yet.** [`Domain::step`] leaves every atom where the dictionary put it, so
/// [`Domain::max_stable_dt`] is infinite — honestly, because nothing moves — and the ledger is
/// empty, because the molecule holds no energy this crate can yet compute. Both change when the
/// force field's energy terms arrive; until then a run of this domain is a picture of a parsed and
/// typed molecule, and is meant to be.
#[derive(Clone, Debug)]
pub struct Molecule {
    name: String,
    component: Component,
    types: Vec<UffType>,
    at: Vec<[f64; 3]>,
    saved: Option<Vec<[f64; 3]>>,
}

impl Molecule {
    /// A molecule named `name`, at the component's coordinates, typed by [`uff::assign`].
    pub fn new(name: impl Into<String>, component: Component) -> Molecule {
        let types = uff::assign(&component);
        let at = component.atoms().iter().map(|a| a.at).collect();
        Molecule {
            name: name.into(),
            component,
            types,
            at,
            saved: None,
        }
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

    /// Empty: the molecule holds no quantity this crate can compute yet, and an entry it could
    /// not stand behind would be worse than none.
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
        vec![
            Reading::new(&self.name, "atoms", atoms.len() as f64, ""),
            Reading::new(&self.name, "heavy atoms", heavy as f64, ""),
            Reading::new(&self.name, "bonds", self.component.bonds().len() as f64, ""),
            Reading::new(&self.name, "formal charge", f64::from(charge), "e"),
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
