//! An AMBER-family non-bonded model on UFF's bonded terms: each atom's van der Waals by its type,
//! read out of an AMBER parameter file, its partial charge from a `.mol2`, and AMBER's two rules.
//!
//! **This crate ships no GAFF value.** GAFF's parameter file, `gaff.dat`, is AmberTools', under
//! GPL v3. This crate is MIT OR Apache-2.0, and that licence should not rest on reading a table of
//! numbers as facts rather than a work. So the caller fetches the file and hands its text to
//! [`VanDerWaals::from_parm`]. The same goes for the charges: AM1-BCC needs a semi-empirical
//! quantum calculation, which this crate does not do, so they come from a `.mol2` that a tool such
//! as antechamber, or a database such as FreeSolv, wrote ([`crate::mol2`]).
//!
//! # What a model is, and what it is not
//!
//! [`NonbondedModel::new`] puts on each atom of a [`Component`]:
//! - its type, the `.mol2`'s type column;
//! - that type's van der Waals, in UFF's form `[x, D]`, where `x = 2 R*/2` and `D = ε`;
//! - its charge, the `.mol2`'s.
//!
//! All three are matched onto the component by graph isomorphism ([`crate::mol2::Mol2::onto`]),
//! never by the file's order.
//!
//! [`NonbondedModel::apply`] sets them on a [`PeriodicForceField`] together with AMBER's rules,
//! Lorentz–Berthelot mixing and the 1-4 scales [`AMBER_ONE_FOUR`]. **The bonds, angles, torsions,
//! inversions and constraints stay UFF's**, so the model is a hybrid: GAFF's non-bonded terms on
//! UFF's bonded ones. That is the model `pantometry-forcefield`'s log P measurements were made
//! with, each from the two direct legs: benzene's 1.92 ± 0.23 against 2.13, and p-xylene's
//! 3.06 ± 0.22 against 3.15.
//!
//! # What is read from a parameter file, and what is refused
//!
//! Every `MOD4      RE` block of a parameter file, and every `NONBON` block of a `frcmod`: one line
//! a type, `type  R*/2  ε` in Å and kcal/mol, and the rest of the line a comment; the block ends
//! at a blank line. Every other section is skipped. Read from `gaff.dat` at ambermini `f7c421c4`
//! on 2026-11-05, it gives that file's 71 types, and the seven the excluded `tests/gaff/` writes
//! out agree with it to the bit. That is a measurement, not a committed check, because the file is
//! GPL. **Every refusal is a [`ParmError`] naming what it refused**:
//! - a file with neither block;
//! - a `MOD4` block of another kind (`SK`, `AC`, `SC`), whose numbers mean something else;
//! - a line short of its two numbers, or a number that is not finite;
//! - an `R*/2` or `ε` below zero, or `R*/2` zero with a well;
//! - one type given twice, in one block or two. A `frcmod` that overrides a parameter file's type
//!   is read on its own and merged by the caller ([`VanDerWaals::with_overrides`]), so that an
//!   override is a decision rather than whichever line came last.
//!
//! Equivalence lines, which give a type another's values (`parm10.dat`'s `C* CA CB …`), are not
//! read. A type that needs one is refused when a model asks for it ([`ModelError::NoSuchType`]),
//! not given a value it might not have.
//!
//! # From a `.mol2` and a parameter file to a force field
//!
//! Methanol, with a `.mol2` and a parameter block written for this example. The types are GAFF's
//! names, but the numbers are invented; with GAFF's own, the caller reads `gaff.dat`'s text instead
//! of `PARM`. Its atoms are written in a different order from the dictionary's, and the map
//! puts each charge on the right atom anyway.
//!
//! ```
//! use pantometry_forcefield::{
//!     uff, Component, EwaldParameters, Mixing, Mol2, NonbondedModel, PeriodicBox,
//!     PeriodicForceField, VanDerWaals, AMBER_ONE_FOUR,
//! };
//!
//! const MOL2: &str = "@<TRIPOS>MOLECULE\nmethanol\n 6 5 1 0 0\nSMALL\nUSER_CHARGES\n\n\
//!     @<TRIPOS>ATOM\n\
//!     1 O1 0 0 0 oh 1 MOL -0.60\n\
//!     2 H4 0 0 0 ho 1 MOL  0.40\n\
//!     3 C1 0 0 0 c3 1 MOL  0.11\n\
//!     4 H1 0 0 0 h1 1 MOL  0.03\n\
//!     5 H2 0 0 0 h1 1 MOL  0.03\n\
//!     6 H3 0 0 0 h1 1 MOL  0.03\n\
//!     @<TRIPOS>BOND\n1 1 2 1\n2 3 1 1\n3 3 4 1\n4 3 5 1\n5 3 6 1\n";
//! const PARM: &str = "INVENTED\n\nMOD4      RE\n\
//!     \x20 c3   1.90  0.10   invented\n\
//!     \x20 h1   1.40  0.02   invented\n\
//!     \x20 oh   1.70  0.20   invented\n\
//!     \x20 ho   0.00  0.00   no van der Waals\n\nEND\n";
//!
//! let methanol = Component::from_ccd(include_str!("../components/MOH.cif"))?;
//! let model = NonbondedModel::new(&Mol2::parse(MOL2)?, &methanol, &VanDerWaals::from_parm(PARM)?)?;
//! let oxygen = methanol.atoms().iter().position(|a| a.name == "O").unwrap();
//! assert_eq!(model.types[oxygen], "oh");
//! assert_eq!(model.charges[oxygen], -0.60);
//!
//! let cell = PeriodicBox::cubic(20e-10);
//! let ewald = EwaldParameters::for_accuracy(&cell, 9e-10, 1e-5);
//! let field = PeriodicForceField::solvated(&methanol, &uff::assign(&methanol), 2, cell, ewald)?;
//! let field = model.apply(field)?;
//! assert_eq!(field.mixing(), Mixing::LorentzBerthelot);
//! assert_eq!(field.one_four(), AMBER_ONE_FOUR);
//! assert_eq!(&field.charges()[..6], &model.charges[..]);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::fmt;

use crate::ccd::Component;
use crate::mol2::{Mol2, Mol2Error};
use crate::periodic::{Mixing, NonbondedError, PeriodicForceField};
use crate::uff::KCAL_PER_MOL;

/// AMBER's and GAFF's 1-4 scales, `[van der Waals, Coulomb]`: `SCNB` = 2 and `SCEE` = 1.2, so
/// `½` of the van der Waals and `1/1.2` of the Coulomb energy.
pub const AMBER_ONE_FOUR: [f64; 2] = [0.5, 1.0 / 1.2];

/// One ångström, metres.
const ANGSTROM: f64 = 1e-10;

/// Each type's van der Waals from an AMBER parameter file, `(type, R*/2 in Å, ε in kcal/mol)`, as
/// written: see the module documentation.
#[derive(Clone, Debug, PartialEq)]
pub struct VanDerWaals {
    entries: Vec<(String, f64, f64)>,
}

/// Why an AMBER parameter file's van der Waals were refused. Every variant names what it refused.
#[derive(Clone, Debug, PartialEq)]
pub enum ParmError {
    /// No `MOD4 RE` block.
    NoBlock,
    /// A `MOD4` block of a kind other than `RE`.
    Kind {
        /// The line, counting from one.
        line: usize,
        /// The kind it names.
        kind: String,
    },
    /// A line of the block short of its type and two numbers.
    Short {
        /// The line, counting from one.
        line: usize,
        /// What it said.
        text: String,
    },
    /// A number that is not finite.
    NotFinite {
        /// The line, counting from one.
        line: usize,
        /// What it said.
        text: String,
    },
    /// An `R*/2` or `ε` not finite or below zero, or an `R*/2` of zero with a well.
    Value {
        /// The type.
        kind: String,
        /// Its `R*/2`, Å.
        radius: f64,
        /// Its `ε`, kcal/mol.
        well: f64,
    },
    /// A type given twice.
    Repeated {
        /// The type.
        kind: String,
    },
}

impl fmt::Display for ParmError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParmError::NoBlock => write!(
                f,
                "no MOD4 RE block and no NONBON block: this file has no van der Waals"
            ),
            ParmError::Kind { line, kind } => write!(
                f,
                "line {line}: a MOD4 {kind} block, whose numbers are not R*/2 and ε; only RE is read"
            ),
            ParmError::Short { line, text } => {
                write!(f, "line {line}: {text:?} is not a type with R*/2 and ε")
            }
            ParmError::NotFinite { line, text } => {
                write!(f, "line {line}: {text:?} is not a finite number")
            }
            ParmError::Value { kind, radius, well } => write!(
                f,
                "type {kind}: R*/2 {radius} Å and ε {well} kcal/mol; each must be finite and not \
                 negative, and R*/2 may be zero only with no well"
            ),
            ParmError::Repeated { kind } => write!(f, "type {kind} is given twice"),
        }
    }
}

impl std::error::Error for ParmError {}

impl VanDerWaals {
    /// **Every `MOD4 RE` block of an AMBER parameter file, and every `NONBON` block of a
    /// `frcmod`**: see the module documentation.
    ///
    /// # Errors
    ///
    /// [`ParmError`], naming what it refused.
    pub fn from_parm(text: &str) -> Result<VanDerWaals, ParmError> {
        let mut entries: Vec<(String, f64, f64)> = Vec::new();
        let mut blocks = 0;
        let mut inside = false;
        for (k, line) in text.lines().enumerate() {
            let at = k + 1;
            let words: Vec<&str> = line.split_whitespace().collect();
            if words.len() == 1 && words[0] == "NONBON" {
                inside = true;
                blocks += 1;
                continue;
            }
            if words.first() == Some(&"MOD4") {
                match words.get(1) {
                    Some(&"RE") => {
                        inside = true;
                        blocks += 1;
                    }
                    other => {
                        return Err(ParmError::Kind {
                            line: at,
                            kind: other.map_or_else(String::new, |s| (*s).to_string()),
                        })
                    }
                }
                continue;
            }
            if !inside {
                continue;
            }
            if words.is_empty() {
                inside = false;
                continue;
            }
            if words.len() < 3 {
                return Err(ParmError::Short {
                    line: at,
                    text: line.into(),
                });
            }
            let number = |w: &str| -> Result<f64, ParmError> {
                match w.parse::<f64>() {
                    Ok(x) if x.is_finite() => Ok(x),
                    _ => Err(ParmError::NotFinite {
                        line: at,
                        text: w.into(),
                    }),
                }
            };
            let (radius, well) = (number(words[1])?, number(words[2])?);
            entries.push((words[0].to_string(), radius, well));
        }
        if blocks == 0 {
            return Err(ParmError::NoBlock);
        }
        VanDerWaals::from_entries(entries)
    }

    /// **Values given directly**, `(type, R*/2 in Å, ε in kcal/mol)`: for a caller that keeps its
    /// own table, held to the same checks as [`VanDerWaals::from_parm`]'s.
    ///
    /// # Errors
    ///
    /// [`ParmError::Value`] or [`ParmError::Repeated`].
    pub fn from_entries(entries: Vec<(String, f64, f64)>) -> Result<VanDerWaals, ParmError> {
        for (i, (kind, radius, well)) in entries.iter().enumerate() {
            if !radius.is_finite()
                || !well.is_finite()
                || *radius < 0.0
                || *well < 0.0
                || (*radius == 0.0 && *well != 0.0)
            {
                return Err(ParmError::Value {
                    kind: kind.clone(),
                    radius: *radius,
                    well: *well,
                });
            }
            if entries[..i].iter().any(|e| e.0 == *kind) {
                return Err(ParmError::Repeated { kind: kind.clone() });
            }
        }
        Ok(VanDerWaals { entries })
    }

    /// **These values with `overrides`' in place of any type both have**, and `overrides`' other
    /// types added after these: how a `frcmod` is laid over the file it modifies.
    pub fn with_overrides(&self, overrides: &VanDerWaals) -> VanDerWaals {
        let mut entries: Vec<(String, f64, f64)> = self
            .entries
            .iter()
            .map(|e| {
                overrides
                    .entries
                    .iter()
                    .find(|o| o.0 == e.0)
                    .unwrap_or(e)
                    .clone()
            })
            .collect();
        for o in &overrides.entries {
            if !self.entries.iter().any(|e| e.0 == o.0) {
                entries.push(o.clone());
            }
        }
        VanDerWaals { entries }
    }

    /// Every type, `(type, R*/2 in Å, ε in kcal/mol)`, in the order read.
    pub fn entries(&self) -> &[(String, f64, f64)] {
        &self.entries
    }

    /// Type `kind`'s `(R*/2 in Å, ε in kcal/mol)`, if the file has it.
    pub fn get(&self, kind: &str) -> Option<(f64, f64)> {
        self.entries
            .iter()
            .find(|e| e.0 == kind)
            .map(|e| (e.1, e.2))
    }

    /// Type `kind`'s van der Waals in UFF's form, `[x, D]` in metres and joules: `x = 2 R*/2`,
    /// `D = ε`. `None` if the file does not have it.
    pub fn uff_form(&self, kind: &str) -> Option<[f64; 2]> {
        self.get(kind)
            .map(|(r, e)| [2.0 * r * ANGSTROM, e * KCAL_PER_MOL])
    }
}

/// Why a [`NonbondedModel`] was not built, or not applied. Every variant names what it refused.
#[derive(Clone, Debug, PartialEq)]
pub enum ModelError {
    /// The `.mol2` is not the component's molecule, or did not map.
    Mol2(Mol2Error),
    /// An atom's type is not in the parameter file.
    NoSuchType {
        /// The component's atom, its index.
        atom: usize,
        /// Its name in the component.
        name: String,
        /// Its type.
        kind: String,
    },
    /// The `.mol2`'s charges do not add to the component's net charge, the sum of its atoms'
    /// formal charges, within `10⁻⁴ e` an atom: more than a file's rounding to four places.
    NetCharge {
        /// What the file's charges add to, elementary charges.
        sum: f64,
        /// The component's net charge.
        formal: i32,
    },
    /// The force field refused the values.
    Nonbonded(NonbondedError),
    /// The model's three columns are not the same length: its fields are public, so a model built
    /// by hand can be ragged.
    Ragged {
        /// How many types.
        types: usize,
        /// How many van der Waals.
        van_der_waals: usize,
        /// How many charges.
        charges: usize,
    },
    /// The model is not one entry per solute atom of the force field it was applied to.
    Atoms {
        /// The model's atoms.
        model: usize,
        /// The force field's solute atoms.
        solute: usize,
    },
}

impl fmt::Display for ModelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ModelError::Mol2(e) => write!(f, "{e}"),
            ModelError::NoSuchType { atom, name, kind } => {
                write!(
                    f,
                    "atom {name} (index {atom}) is type {kind}, which the parameter file has no \
                     van der Waals for"
                )?;
                if kind.contains('.') {
                    write!(
                        f,
                        "; {kind} is a SYBYL type, so the .mol2 was perhaps written without \
                         GAFF's types (antechamber's -at gaff)"
                    )?;
                }
                Ok(())
            }
            ModelError::NetCharge { sum, formal } => write!(
                f,
                "the .mol2's charges add to {sum:+.4} e, and the molecule's net charge is {formal:+}"
            ),
            ModelError::Nonbonded(e) => write!(f, "{e}"),
            ModelError::Ragged {
                types,
                van_der_waals,
                charges,
            } => write!(
                f,
                "a model of {types} types, {van_der_waals} van der Waals and {charges} charges"
            ),
            ModelError::Atoms { model, solute } => write!(
                f,
                "a model of {model} atoms on a force field whose solute has {solute}; for several \
                 copies of one molecule, apply NonbondedModel::repeated"
            ),
        }
    }
}

impl std::error::Error for ModelError {}

impl From<Mol2Error> for ModelError {
    fn from(e: Mol2Error) -> ModelError {
        ModelError::Mol2(e)
    }
}

/// **A molecule's AMBER-family non-bonded model**: each atom's type, van der Waals `[x, D]` in
/// UFF's form, and partial charge, in the component's order. See the module documentation.
#[derive(Clone, Debug, PartialEq)]
pub struct NonbondedModel {
    /// Each atom's type, from the `.mol2`.
    pub types: Vec<String>,
    /// Each atom's `[x, D]`, metres and joules, from the parameter file by its type.
    pub van_der_waals: Vec<[f64; 2]>,
    /// Each atom's charge, elementary charges, from the `.mol2`.
    pub charges: Vec<f64>,
}

impl NonbondedModel {
    /// **`mol2`'s types and charges on `c`, and each type's van der Waals from `parm`**, matched
    /// by graph isomorphism.
    ///
    /// # Errors
    ///
    /// [`ModelError::Mol2`] when `mol2` is not `c`'s molecule; [`ModelError::NetCharge`] when its
    /// charges do not add to `c`'s net charge; [`ModelError::NoSuchType`] when an atom's type is
    /// not in `parm`.
    pub fn new(
        mol2: &Mol2,
        c: &Component,
        parm: &VanDerWaals,
    ) -> Result<NonbondedModel, ModelError> {
        let mapped = mol2.onto(c)?;
        let sum: f64 = mapped.charges.iter().sum();
        let formal: i32 = c.atoms().iter().map(|a| a.charge).sum();
        if (sum - f64::from(formal)).abs() > 1e-4 * mapped.charges.len() as f64 {
            return Err(ModelError::NetCharge { sum, formal });
        }
        let van_der_waals = mapped
            .kinds
            .iter()
            .enumerate()
            .map(|(atom, kind)| {
                parm.uff_form(kind).ok_or_else(|| ModelError::NoSuchType {
                    atom,
                    name: c.atoms()[atom].name.clone(),
                    kind: kind.clone(),
                })
            })
            .collect::<Result<_, _>>()?;
        Ok(NonbondedModel {
            types: mapped.kinds,
            van_der_waals,
            charges: mapped.charges,
        })
    }

    /// How many atoms the model is for.
    pub fn atoms(&self) -> usize {
        self.types.len()
    }

    /// **`copies` of this model one after another**: for a solute that is several molecules of
    /// one kind written in turn, as a liquid's box is.
    pub fn repeated(&self, copies: usize) -> NonbondedModel {
        NonbondedModel {
            types: (0..copies).flat_map(|_| self.types.clone()).collect(),
            van_der_waals: self.van_der_waals.repeat(copies),
            charges: self.charges.repeat(copies),
        }
    }

    /// **The model on `field`'s solute, with AMBER's rules**: the solute's van der Waals and
    /// charges this model's, the waters' and any ions' left as they are, and the pairs mixed by
    /// Lorentz–Berthelot and every 1-4 pair scaled by [`AMBER_ONE_FOUR`]. **The rules are
    /// field-wide**: the solute's own pairs follow them, and so do its pairs with the waters and
    /// ions, which [`PeriodicForceField::solvated`] otherwise combines by UFF's geometric rule.
    /// `field`'s solute atoms come first, as [`PeriodicForceField::solvated`] puts them; for
    /// several copies of one molecule, apply [`NonbondedModel::repeated`].
    ///
    /// # Errors
    ///
    /// [`ModelError::Ragged`] when the model's columns are not one length;
    /// [`ModelError::Atoms`] when the model is not one entry per solute atom;
    /// [`ModelError::Nonbonded`] when the force field refuses a value. Nothing is changed then.
    pub fn apply(&self, field: PeriodicForceField) -> Result<PeriodicForceField, ModelError> {
        let n = self.atoms();
        if self.van_der_waals.len() != n || self.charges.len() != n {
            return Err(ModelError::Ragged {
                types: n,
                van_der_waals: self.van_der_waals.len(),
                charges: self.charges.len(),
            });
        }
        if n != field.solute_atoms() {
            return Err(ModelError::Atoms {
                model: n,
                solute: field.solute_atoms(),
            });
        }
        Ok(field
            .with_nonbonded(0..n, self.van_der_waals.clone())
            .map_err(ModelError::Nonbonded)?
            .with_solute_charges(self.charges.clone())
            .with_mixing(Mixing::LorentzBerthelot)
            .with_one_four(AMBER_ONE_FOUR[0], AMBER_ONE_FOUR[1]))
    }
}
