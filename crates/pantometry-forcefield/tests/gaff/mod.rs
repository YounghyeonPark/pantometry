//! **GAFF's van der Waals for octanol's five types and benzene's two**, with where they were read,
//! GAFF's 1-4 scales, and each molecule's GAFF/AM1-BCC non-bonded model built from FreeSolv's files.
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
use pantometry_forcefield::uff::KCAL_PER_MOL;
use pantometry_forcefield::{Component, Mixing, PeriodicForceField};

const ANGSTROM: f64 = 1e-10;

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
/// `1/1.2` of the Coulomb energy.
pub const ONE_FOUR: [f64; 2] = [0.5, 1.0 / 1.2];

/// GAFF's `[x, D]` in UFF's form, `x = 2 (R*/2)` and `D = ε`, for type `kind`.
pub fn gaff(kind: &str) -> [f64; 2] {
    let &(_, r, e) = GAFF
        .iter()
        .find(|t| t.0 == kind)
        .unwrap_or_else(|| panic!("no GAFF type {kind} here"));
    [2.0 * r * ANGSTROM, e * KCAL_PER_MOL]
}

/// Each atom of `c`'s GAFF type, from FreeSolv's `mol2` of the same molecule by graph isomorphism.
pub fn types_onto(mol2: &freesolv::Mol2, c: &Component) -> Vec<String> {
    let (elements, bonds) = freesolv::component_graph(c);
    let mapped = freesolv::charges_onto(mol2, &elements, &bonds).expect("the same molecule");
    mapped
        .map
        .iter()
        .map(|&j| mol2.atoms[j].kind.clone())
        .collect()
}

/// Each atom of the dictionary's octanol's GAFF type, from FreeSolv's `mobley_36119.mol2`.
pub fn gaff_types(c: &Component) -> Vec<String> {
    types_onto(&freesolv::octanol(), c)
}

/// **A molecule's GAFF/AM1-BCC non-bonded model**: each atom's GAFF `[x, D]` by its type in
/// FreeSolv's file, and FreeSolv's AM1-BCC charge, both mapped onto `c` by graph isomorphism.
#[derive(Clone, Debug, PartialEq)]
pub struct Model {
    pub types: Vec<String>,
    pub vdw: Vec<[f64; 2]>,
    pub charges: Vec<f64>,
}

impl Model {
    fn from(mol2: &freesolv::Mol2, c: &Component) -> Model {
        let (elements, bonds) = freesolv::component_graph(c);
        let mapped = freesolv::charges_onto(mol2, &elements, &bonds).expect("the same molecule");
        let types: Vec<String> = mapped
            .map
            .iter()
            .map(|&j| mol2.atoms[j].kind.clone())
            .collect();
        Model {
            vdw: types.iter().map(|k| gaff(k)).collect(),
            charges: mapped.charges,
            types,
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
}

/// **GAFF/AM1-BCC benzene in TIP3P**: `field` a [`PeriodicForceField::solvated`] whose solute is
/// benzene `c`, its atoms first; benzene's [`Model`] on them, the waters TIP3P's — TIP3P's oxygen
/// is GAFF's `ow` to `gaff.dat`'s four places — and GAFF's two rules, which are field-wide:
/// Lorentz–Berthelot and [`ONE_FOUR`], so benzene's own 1-4 pairs are scaled in vacuum too. Variant 4's benzene, which the
/// staged run swapped in on W3's UFF leg; but here its own 21 1-4 pairs are scaled and its own
/// pairs are GAFF's and AM1-BCC's, where the staged run kept UFF's and QEq's.
pub fn benzene_in_water(field: PeriodicForceField, c: &Component) -> PeriodicForceField {
    let b = Model::benzene(c);
    let n = b.vdw.len();
    field
        .with_nonbonded(0..n, b.vdw)
        .expect("GAFF's values are valid")
        .with_solute_charges(b.charges)
        .with_mixing(Mixing::LorentzBerthelot)
        .with_one_four(ONE_FOUR[0], ONE_FOUR[1])
}
