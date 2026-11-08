//! **GAFF's van der Waals for octanol's five types and benzene's two**, with where they were read,
//! GAFF's 1-4 scales, and each molecule's GAFF/AM1-BCC non-bonded model built from FreeSolv's files.
//! p-xylene's four types are among those seven (`c3`, `hc`, `ca`, `ha`), so its model needs no
//! value read beyond them, and so are phenol's (`ca`, `ha`, `oh`, `ho`).
//!
//! **Kept out of the published crate** (`Cargo.toml`'s `exclude`): the numbers are read from
//! `gaff.dat` in ambermini, which is AmberTools under GPL v3, as `tests/amber/`'s are. The values
//! may well be facts rather than a work; MIT OR Apache-2.0 should not rest on that. Used by
//! `liquid_octanol_with_gaff.rs` (octanol's), `benzene_partitioned_into_octanol.rs` (both) and
//! `benzene_bound_in_tip3p.rs`, which holds benzene's model as built here to the one variant 4's
//! staged run used.
//!
//! **GAFF's bonded terms are not used anywhere**: the bonds, angles, torsions, inversions, the
//! constraints' lengths and every conformer are UFF's. So each model is a hybrid: GAFF's
//! non-bonded terms on UFF's bonded ones.

#![allow(dead_code)]

use crate::freesolv;
use pantometry_forcefield::{Component, NonbondedModel, PeriodicForceField, VanDerWaals};

/// Where [`GAFF`] was read: ambermini at the commit W4's ff14SB tables cite.
pub const GAFF_SOURCE: &str = "github.com/choderalab/ambermini \
                               f7c421c46dca9f0e0f87b10a9c87985b446d8ff9 \
                               share/amber/dat/leap/parm/gaff.dat";

/// The SHA-256 of [`GAFF_SOURCE`]'s bytes, reproduced on 2026-10-21 (twice) and the one W4 cites.
pub const GAFF_SHA256: &str = "95e3c655c225064435742ab6414a87e6c30a559e8be60de0482eba55144d8433";

/// **GAFF 1.7's van der Waals for octanol's five types and benzene's two**,
/// `(type, R*/2 in Å, ε in kcal/mol)`, as `gaff.dat`'s `MOD4 RE` block writes them (its first
/// line: "AMBER General Force Field for organic molecules (Version 1.7, Nov 2013)"), lines 5892,
/// 5897, 5898, 5900, 5906, 5912 and 5913:
///
/// ```text
///   h1          1.3870  0.0157             Veenstra et al JCC,8,(1992),963
///   ha          1.4590  0.0150             Spellmeyer
///   hc          1.4870  0.0157             OPLS
///   ho          0.0000  0.0000             OPLS Jorgensen, JACS,110,(1988),1657
///   oh          1.7210  0.2104             OPLS
///   c3          1.9080  0.1094             OPLS
///   ca          1.9080  0.0860             OPLS
/// ```
///
/// `ca` and `ha` are the values `benzene_bound_in_tip3p.rs`'s `gaff_benzene` writes from the same
/// file, read again here rather than copied, and `the_direct_routes_benzene_is_variant_4s` there
/// holds the two to the bit.
pub const GAFF: [(&str, f64, f64); 7] = [
    ("c3", 1.9080, 0.1094),
    ("hc", 1.4870, 0.0157),
    ("h1", 1.3870, 0.0157),
    ("oh", 1.7210, 0.2104),
    ("ho", 0.0, 0.0),
    ("ca", 1.9080, 0.0860),
    ("ha", 1.4590, 0.0150),
];

/// GAFF's and AMBER's 1-4 scales: `SCNB` = 2 and `SCEE` = 1.2, so `½` of the van der Waals and
/// `1/1.2` of the Coulomb energy. The same values as the library's
/// [`pantometry_forcefield::AMBER_ONE_FOUR`], which `a_mol2_and_a_parameter_file_make_a_model.rs`
/// holds to them and which [`solute_in_water`] applies.
pub const ONE_FOUR: [f64; 2] = [0.5, 1.0 / 1.2];

/// [`GAFF`] as the library's table ([`VanDerWaals::from_entries`]), so that every model here is
/// built by the library's lines.
pub fn table() -> VanDerWaals {
    VanDerWaals::from_entries(
        GAFF.iter()
            .map(|&(k, r, e)| (k.to_string(), r, e))
            .collect(),
    )
    .expect("GAFF's values are valid")
}

/// GAFF's `[x, D]` in UFF's form, `x = 2 (R*/2)` and `D = ε`, for type `kind`.
pub fn gaff(kind: &str) -> [f64; 2] {
    table()
        .uff_form(kind)
        .unwrap_or_else(|| panic!("no GAFF type {kind} here"))
}

/// Each atom of `c`'s GAFF type, from FreeSolv's `mol2` of the same molecule by graph isomorphism.
pub fn types_onto(mol2: &freesolv::Mol2, c: &Component) -> Vec<String> {
    mol2.onto(c).expect("the same molecule").kinds
}

/// Each atom of the dictionary's octanol's GAFF type, from FreeSolv's `mobley_36119.mol2`.
pub fn gaff_types(c: &Component) -> Vec<String> {
    types_onto(&freesolv::octanol(), c)
}

/// **A molecule's GAFF/AM1-BCC non-bonded model**: each atom's GAFF `[x, D]` by its type in
/// FreeSolv's file, and FreeSolv's AM1-BCC charge, both mapped onto `c` by graph isomorphism. The
/// library's [`NonbondedModel`] with [`GAFF`]'s table, under the names this directory's tests use.
#[derive(Clone, Debug, PartialEq)]
pub struct Model {
    pub types: Vec<String>,
    pub vdw: Vec<[f64; 2]>,
    pub charges: Vec<f64>,
}

impl Model {
    fn from(mol2: &freesolv::Mol2, c: &Component) -> Model {
        let m = NonbondedModel::new(mol2, c, &table()).expect("the same molecule, GAFF's types");
        Model {
            types: m.types,
            vdw: m.van_der_waals,
            charges: m.charges,
        }
    }

    /// The library's model this is.
    pub fn library(&self) -> NonbondedModel {
        NonbondedModel {
            types: self.types.clone(),
            van_der_waals: self.vdw.clone(),
            charges: self.charges.clone(),
        }
    }

    /// Benzene's: `ca` and `ha` from `mobley_3053621.mol2`, and its AM1-BCC charges as FreeSolv ran
    /// them, summing to −0.0001 e.
    pub fn benzene(c: &Component) -> Model {
        Model::from(&freesolv::benzene(), c)
    }

    /// Octanol's: `c3`, `hc`, `h1`, `oh` and `ho` from `mobley_36119.mol2`, and its AM1-BCC
    /// charges, summing to zero.
    pub fn octanol(c: &Component) -> Model {
        Model::from(&freesolv::octanol(), c)
    }

    /// p-xylene's: `c3`, `hc`, `ca` and `ha` from `mobley_3234716.mol2`, and its AM1-BCC charges as
    /// FreeSolv ran them, summing to +0.0001 e.
    pub fn p_xylene(c: &Component) -> Model {
        Model::from(&freesolv::p_xylene(), c)
    }

    /// Phenol's: `ca`, `ha`, `oh` and `ho` from `mobley_20524.mol2`, and its AM1-BCC charges as
    /// FreeSolv ran them, summing to zero.
    pub fn phenol(c: &Component) -> Model {
        Model::from(&freesolv::phenol(), c)
    }
}

/// **GAFF/AM1-BCC benzene in TIP3P**: `field` a [`PeriodicForceField::solvated`] whose solute is
/// benzene `c`, its atoms first; benzene's [`Model`] on them, the waters TIP3P's — TIP3P's oxygen
/// is GAFF's `ow` to `gaff.dat`'s four places — and GAFF's two rules, which are field-wide:
/// Lorentz–Berthelot and [`ONE_FOUR`], so benzene's own 1-4 pairs are scaled in vacuum too. Variant 4's benzene, which the
/// staged run swapped in on W3's UFF leg; but here its own 21 1-4 pairs are scaled and its own
/// pairs are GAFF's and AM1-BCC's, where the staged run kept UFF's and QEq's.
pub fn benzene_in_water(field: PeriodicForceField, c: &Component) -> PeriodicForceField {
    solute_in_water(field, Model::benzene(c))
}

/// **Any solute's GAFF/AM1-BCC `model` in TIP3P**, by the same rules as [`benzene_in_water`],
/// which is this with benzene's model: `field` a [`PeriodicForceField::solvated`] whose solute's
/// atoms come first.
pub fn solute_in_water(field: PeriodicForceField, model: Model) -> PeriodicForceField {
    model
        .library()
        .apply(field)
        .expect("GAFF's values are valid, one per solute atom")
}
