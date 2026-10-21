//! **Liquid 1-octanol with GAFF's non-bonded terms on UFF's bonded ones**: model (b) of step 4c-1,
//! which needed [`PeriodicForceField::with_nonbonded`], [`PeriodicForceField::with_mixing`] and
//! [`PeriodicForceField::with_one_four`]. GAFF's Lennard-Jones for its five types (`c3`, `hc`,
//! `h1`, `oh`, `ho`), combined by Lorentz–Berthelot, GAFF's 1-4 scales (½ van der Waals, 1/1.2
//! Coulomb) and FreeSolv's AM1-BCC charges (`mobley_36119`); **GAFF's bonded terms are not used**:
//! the bonds, angles, torsions and inversions are UFF's, and so are the constraints' lengths and
//! the conformer the lattice is built from. So this is a hybrid, GAFF's non-bonded model on UFF's
//! geometry, and not GAFF.
//!
//! **Kept out of the published crate** (`Cargo.toml`'s `exclude`): [`GAFF`]'s numbers are read from
//! `gaff.dat` in ambermini, which is AmberTools under GPL v3, as `tests/amber/`'s are. The values
//! may well be facts rather than a work; MIT OR Apache-2.0 should not rest on that.
//!
//! Default, unoptimised: GAFF's types are the FreeSolv file's, mapped by graph isomorphism onto
//! the dictionary's `OC9`, and the ones a chemist reads off the molecule; and the gas phase's
//! potential, written out in the tests ([`nonbonded::Vacuum`]), has forces that are its energy's
//! derivatives. Ignored, release: [`the_liquid_with_gaff_measured`], 4c-1's protocol for this model.
//!
//! # What the liquid measured: within 6 kJ/mol of experiment, and hydrogen-bonded
//!
//! Run once on 2026-10-21, release, one core, 25 000 steps of 2 fs, 1 976 s; nothing asserted.
//! Beside 4c-1's two UFF models, from `liquid_octanol_in_a_box.rs`:
//!
//! | | (b) GAFF's LJ, LB, 1-4, AM1-BCC | (a) UFF, QEq | (c) UFF, AM1-BCC | experiment |
//! | --- | --- | --- | --- | --- |
//! | P at 0.8262 g/cm³ (bar) | −178 ± 46 | −930 ± 45 | −788 ± 56 | 1 |
//! | ΔH_vap (kJ/mol) | 77.3 ± 0.7 | 91.9 ± 0.9 | 97.3 ± 0.7 | 70.98 (CRC, read via PubChem) |
//! | U/N in the sample's halves (kcal/mol) | 1.369, 1.360 | 18.00, 17.89 | −1.23, −1.59 | |
//! | g_OO's first peak | 2.824 Å, 10.61; 1.79 O inside 3.71 Å | 3.39 Å, 0.82; 0.06 O | 3.37 Å, 2.47; 1.40 O | an alcohol's O–H···O near 2.7 Å (quoted) |
//!
//! - **The excess cohesion is mostly gone**: 6.3 kJ/mol (9%) over experiment against UFF's 21–26,
//!   and the pull at experiment's density a fifth of UFF's. `U_gas` is 19.238 ± 0.153 kcal/mol,
//!   eight runs of 200 ps in [`gas_phase`]'s potential.
//! - **The hydroxyls hydrogen-bond**: `ho` has no van der Waals, so its hydrogen reaches an oxygen,
//!   and the O–O peak is at 2.82 Å with 1.79 neighbours inside its minimum — two, as in a chain of
//!   hydrogen bonds — where UFF had none nearer than 3.3 Å.
//! - **Settled**: `U/N`'s halves agree to 0.009 kcal/mol and its τ is 2.4 samples of 0.1 ps; `P`'s
//!   halves, −120 ± 73 and −237 ± 23 bar, differ by 1.5 of their combined error.
//! - **What is not GAFF here**: the bonded terms and the geometry are UFF's, so the remaining
//!   6 kJ/mol is this hybrid's, and says nothing about GAFF's own octanol.

// The sums over atoms and axes are written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

mod freesolv;
mod nonbonded;
mod octanol;

use nonbonded::Vacuum;
use octanol::*;
use pantometry_forcefield::uff::KCAL_PER_MOL;
use pantometry_forcefield::{Component, Element, Mixing, PeriodicForceField, Potential};

/// Where [`GAFF`] was read: ambermini at the commit W4's ff14SB tables cite.
pub const GAFF_SOURCE: &str = "github.com/choderalab/ambermini \
                               f7c421c46dca9f0e0f87b10a9c87985b446d8ff9 \
                               share/amber/dat/leap/parm/gaff.dat";

/// The SHA-256 of [`GAFF_SOURCE`]'s bytes, reproduced on 2026-10-21 and the one W4 cites.
pub const GAFF_SHA256: &str = "95e3c655c225064435742ab6414a87e6c30a559e8be60de0482eba55144d8433";

/// **GAFF 1.7's van der Waals for octanol's five types**, `(type, R*/2 in Å, ε in kcal/mol)`, as
/// `gaff.dat`'s `MOD4 RE` block writes them (its first line: "AMBER General Force Field for
/// organic molecules (Version 1.7, Nov 2013)"), lines 5892–5912:
///
/// ```text
///   h1          1.3870  0.0157             Veenstra et al JCC,8,(1992),963
///   hc          1.4870  0.0157             OPLS
///   ho          0.0000  0.0000             OPLS Jorgensen, JACS,110,(1988),1657
///   oh          1.7210  0.2104             OPLS
///   c3          1.9080  0.1094             OPLS
/// ```
pub const GAFF: [(&str, f64, f64); 5] = [
    ("c3", 1.9080, 0.1094),
    ("hc", 1.4870, 0.0157),
    ("h1", 1.3870, 0.0157),
    ("oh", 1.7210, 0.2104),
    ("ho", 0.0, 0.0),
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

/// Each atom of the dictionary's octanol's GAFF type, from FreeSolv's `mobley_36119.mol2` by graph
/// isomorphism.
pub fn gaff_types(c: &Component) -> Vec<String> {
    let mol2 = freesolv::octanol();
    let (elements, bonds) = freesolv::component_graph(c);
    let mapped = freesolv::charges_onto(&mol2, &elements, &bonds).expect("the same molecule");
    mapped
        .map
        .iter()
        .map(|&j| mol2.atoms[j].kind.clone())
        .collect()
}

/// Octanol with FreeSolv's AM1-BCC charges.
fn am1bcc_octanol() -> Octanol {
    let qeq = Octanol::new(None);
    let (elements, bonds) = freesolv::component_graph(&qeq.component);
    let mapped =
        freesolv::charges_onto(&freesolv::octanol(), &elements, &bonds).expect("the same molecule");
    Octanol::new(Some(mapped.charges))
}

/// One octanol in vacuum with GAFF's non-bonded terms on UFF's bonded ones.
fn gas_phase(o: &Octanol) -> Vacuum {
    let values: Vec<[f64; 2]> = gaff_types(&o.component).iter().map(|k| gaff(k)).collect();
    Vacuum::new(
        &o.component,
        &o.types,
        values,
        o.charges.clone(),
        Mixing::LorentzBerthelot,
        ONE_FOUR,
    )
}

/// `field`, a box of whole octanols, with GAFF's van der Waals on every atom, Lorentz–Berthelot and
/// GAFF's 1-4 scales.
fn with_gaff(field: PeriodicForceField, o: &Octanol) -> PeriodicForceField {
    let one: Vec<[f64; 2]> = gaff_types(&o.component).iter().map(|k| gaff(k)).collect();
    let n = field.charges().len();
    assert_eq!(n % one.len(), 0, "whole octanols");
    let values: Vec<[f64; 2]> = (0..n).map(|i| one[i % one.len()]).collect();
    field
        .with_nonbonded(0..n, values)
        .expect("GAFF's values are valid")
        .with_mixing(Mixing::LorentzBerthelot)
        .with_one_four(ONE_FOUR[0], ONE_FOUR[1])
}

/// **GAFF's types are FreeSolv's file's, and the ones a chemist reads off octanol**: every carbon
/// `c3`; the oxygen `oh` and its hydrogen `ho`; the two hydrogens on the carbon bonded to the
/// oxygen `h1`, one electron-withdrawing neighbour; every other hydrogen `hc` — read off the
/// dictionary's bonds here, and equal to the file's mapped by graph isomorphism. And the values in
/// UFF's form are `x = 2 R*/2`, `D = ε`.
#[test]
fn gaffs_types_are_the_files() {
    let c = Component::from_ccd(OC9).expect("OC9");
    let kinds = gaff_types(&c);
    let mut neighbours = vec![Vec::new(); c.atoms().len()];
    for b in c.bonds() {
        let [i, j] = b.atoms;
        neighbours[i].push(j);
        neighbours[j].push(i);
    }
    let element = |i: usize| c.atoms()[i].element;
    for (i, a) in c.atoms().iter().enumerate() {
        let want = match a.element {
            Element::C => "c3",
            Element::O => "oh",
            Element::H => {
                let heavy = neighbours[i][0];
                if element(heavy) == Element::O {
                    "ho"
                } else if neighbours[heavy].iter().any(|&k| element(k) == Element::O) {
                    "h1"
                } else {
                    "hc"
                }
            }
            other => panic!("octanol has no {other:?}"),
        };
        assert_eq!(kinds[i], want, "{}", a.name);
    }
    let count = |k: &str| kinds.iter().filter(|x| *x == k).count();
    assert_eq!(
        [
            count("c3"),
            count("oh"),
            count("ho"),
            count("h1"),
            count("hc")
        ],
        [8, 1, 1, 2, 15]
    );
    assert_eq!(gaff("c3"), [2.0 * 1.9080 * ANGSTROM, 0.1094 * KCAL_PER_MOL]);
    assert_eq!(gaff("ho"), [0.0, 0.0]);
    println!("GAFF from {GAFF_SOURCE}, SHA-256 {GAFF_SHA256}");
}

/// **The gas phase's forces are its energy's derivatives**: [`nonbonded::Vacuum`] with GAFF's
/// terms on octanol's relaxed geometry, every atom along every axis by central differences, each
/// to the tolerance it earns — `64 ε` of the energy's terms over `h`, and the truncation by the
/// difference at `2h`. It is the potential [`the_liquid_with_gaff_measured`] samples the gas in.
#[test]
fn the_gas_phases_forces_are_its_energys_derivatives() {
    let o = am1bcc_octanol();
    let gas = gas_phase(&o);
    let n = o.at.len();
    let mut f = vec![[0.0; 3]; n];
    gas.energy_and_forces(&o.at, &mut f);
    let size = gas.nonbonded(&o.at).magnitude + gas.bonded.energy(&o.at).total.abs();
    let h = 1e-5 * ANGSTROM;
    let e = |at: &[[f64; 3]]| {
        let mut g = vec![[0.0; 3]; n];
        gas.energy_and_forces(at, &mut g)
    };
    let mut worst = 0.0f64;
    for i in 0..n {
        for a in 0..3 {
            let d = |s: f64| {
                let mut p = o.at.clone();
                p[i][a] += s;
                let mut m = o.at.clone();
                m[i][a] -= s;
                (e(&p) - e(&m)) / (2.0 * s)
            };
            let (d1, d2) = (d(h), d(2.0 * h));
            let tolerance = 2.0 * (d2 - d1).abs() / 3.0 + 64.0 * f64::EPSILON * size / h;
            assert!(
                (f[i][a] + d1).abs() <= tolerance,
                "atom {i} axis {a}: {:e} against {:e} ± {tolerance:e}",
                f[i][a],
                -d1
            );
            worst = worst.max((f[i][a] + d1).abs() / tolerance);
        }
    }
    println!("the worst |F + dE/dx| is {worst:.3} of its tolerance");
}

/// **The liquid, measured**, model (b): GAFF's van der Waals by Lorentz–Berthelot, GAFF's 1-4
/// scales and AM1-BCC's charges, on UFF's bonded terms — [`octanol::measure_liquid`], 4c-1's
/// protocol to the step, with the gas sampled in [`gas_phase`]'s potential. Nothing about the
/// model is asserted.
#[test]
#[ignore = "192 octanols melted and sampled, about an hour with --release: run with --release -- --ignored --nocapture"]
fn the_liquid_with_gaff_measured() {
    let o = am1bcc_octanol();
    let gas = gas_phase(&o);
    measure_liquid(
        "GAFF's van der Waals and 1-4 scales, Lorentz–Berthelot, AM1-BCC's charges, UFF's bonded terms",
        &o,
        &|field| with_gaff(field, &o),
        &gas,
    );
}

/// GAFF on a box of sixteen octanols at the liquid's density, with `ho` given `[x, 0]`.
fn small_gaff_box(o: &Octanol, ho_x: f64) -> (PeriodicForceField, Vec<[f64; 3]>) {
    let l = Liquid::lattice(o, [1, 4, 4], DENSITY, 0x6AFF, 6.5 * ANGSTROM, 1e-5);
    let kinds = gaff_types(&o.component);
    let n = l.at.len();
    let values: Vec<[f64; 2]> = (0..n)
        .map(|i| match kinds[i % kinds.len()].as_str() {
            "ho" => [ho_x, 0.0],
            k => gaff(k),
        })
        .collect();
    let field = l
        .field
        .clone()
        .with_nonbonded(0..n, values)
        .expect("GAFF's values are valid")
        .with_mixing(Mixing::LorentzBerthelot)
        .with_one_four(ONE_FOUR[0], ONE_FOUR[1]);
    (field, l.at)
}

/// **GAFF's tail is its double sum**: the sixteen-octanol box's five classes — `c3`, `hc`, `h1`,
/// `oh`, `ho`, with `hc` and `h1` sharing their `ε` — summed over every ordered pair of its 432
/// atoms, written out, energy and pressure, within `(N² + 64) ε` of the parts.
#[test]
fn gaffs_tail_is_its_double_sum() {
    let o = am1bcc_octanol();
    let (field, _) = small_gaff_box(&o, 0.0);
    let vdw = field.van_der_waals().to_vec();
    let mut classes: Vec<[u64; 2]> = vdw.iter().map(|x| x.map(f64::to_bits)).collect();
    classes.sort_unstable();
    classes.dedup();
    assert_eq!(classes.len(), 5);
    assert_eq!(gaff("hc")[1], gaff("h1")[1]);
    let rc = field.ewald().parameters().cutoff;
    let v = field.cell().volume();
    let ([e, p], [se, sp]) = nonbonded::tail_by_pairs(&vdw, &vdw, Mixing::LorentzBerthelot, rc, v);
    let n = vdw.len() as f64;
    println!(
        "GAFF: E_tail {:+.12e} J against {e:+.12e}; P_tail {:+.12e} Pa against {p:+.12e}",
        field.dispersion_correction(),
        field.pressure_correction()
    );
    assert!((field.dispersion_correction() - e).abs() <= (n * n + 64.0) * f64::EPSILON * se);
    assert!((field.pressure_correction() - p).abs() <= (n * n + 64.0) * f64::EPSILON * sp);
}

/// **`ho` as `[0, 0]` is `ho` as `[1 Å, 0]` to the bit**: the measurement was run with the second,
/// when the force field took only `x > 0`; with `D = 0` every pair of the atom is zero either way,
/// and the evaluation's every term, force and virial component, and both corrections, are the same
/// bits — so the trajectory, and what it measured, are too.
#[test]
fn ho_with_no_size_is_ho_with_one_angstrom() {
    let o = am1bcc_octanol();
    let (zero, at) = small_gaff_box(&o, 0.0);
    let (one, _) = small_gaff_box(&o, 1.0 * ANGSTROM);
    let bits = |f: &PeriodicForceField| {
        let ev = f.evaluate(&at);
        let e = &ev.energy;
        let mut out: Vec<u64> = [
            e.bond,
            e.angle,
            e.torsion,
            e.inversion,
            e.van_der_waals,
            e.dispersion_correction,
            e.electrostatic,
            e.total,
            f.dispersion_correction(),
            f.pressure_correction(),
        ]
        .map(f64::to_bits)
        .to_vec();
        out.extend(ev.forces.iter().flatten().map(|x| x.to_bits()));
        out.extend(ev.virial.iter().flatten().map(|x| x.to_bits()));
        out
    };
    assert_eq!(bits(&zero), bits(&one));
}
