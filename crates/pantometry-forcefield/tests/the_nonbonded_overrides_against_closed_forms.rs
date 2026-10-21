//! **Another force field's non-bonded terms in the periodic box**, against closed forms and exact
//! identities: [`PeriodicForceField::with_nonbonded`], [`PeriodicForceField::with_mixing`] and
//! [`PeriodicForceField::with_one_four`]. The values overridden with are invented
//! ([`nonbonded::invented`]), no force field's, so that this file carries nothing but arithmetic;
//! GAFF's are in `liquid_octanol_with_gaff.rs`, which is kept out of the published crate.
//!
//! - **Every rule at its default is the plain field's bits**: overrides equal to the UFF types' own
//!   values, written in two runs, with geometric mixing and `(1, 1)`, give the same energy by term,
//!   forces and virial, and the same decoupled energies and couplings, through classical Ewald, the
//!   mesh and the neighbour list.
//! - **One pair is Lorentz–Berthelot's formula** in σ and ε, energy and force, at four distances,
//!   and visibly not the geometric mean's.
//! - **The long-range correction under Lorentz–Berthelot is its double sum** over every ordered
//!   pair of atoms, written out, on a box of four classes; with one class it is Allen and
//!   Tildesley's printed form; and the decoupling's cross part is its own double sum.
//! - **The 1-4 pairs are the pairs three bonds apart by breadth-first search**: octanol's, benzene's
//!   (para carbons by two paths, once) and cyclopentane's, whose ring has no 1-4 pair across it.
//! - **The forces are `−∇E` and the virial is `−∂U/∂ε`** by central differences, with every rule on,
//!   through classical Ewald and the mesh.
//! - **A 1-4 pair's scaled Coulomb through Ewald**: one octanol in boxes of 40, 60 and 80 Å is its
//!   vacuum energy with the 1-4 pairs scaled, less the tinfoil term `2π k_e μ²/3V`, the rest
//!   falling as `L⁻⁵`, energy and forces; its van der Waals is the direct pair sum.
//! - **Decoupling with every rule on**: benzene in octanol, the end states the systems they say, the
//!   λ derivatives the whole energy's central differences, `couplings` `coupling`'s bits and the
//!   evaluation's λ-dependent part, and the decoupled reciprocal energy `rest + λ_e · cross`.
//! - **Determinism**: dynamics the same bits however it is cut, and the neighbour list the cell
//!   list's bits along it.

// The sums over atoms and axes are written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

mod common;
mod freesolv;
mod hydration;
mod nonbonded;
mod octanol;
mod pinned;

use nonbonded::{
    bonds, graph_distances, invented, one_four_by_search, tail_by_pairs, uff_values, Vacuum,
};
use octanol::{index, Liquid, Octanol, Solute, ANGSTROM, DENSITY, FS};
use pantometry_core::math::{exp, ln};
use pantometry_forcefield::ewald::COULOMB;
use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::water::{hoh_angle, OH_LENGTH};
use pantometry_forcefield::{
    Alchemical, Bath, Component, Element, EwaldParameters, Ion, Lambda, Mixing, NonbondedError,
    PeriodicBox, PeriodicDecoupling, PeriodicEnergy, PeriodicEvaluation, PeriodicForceField,
    PmeParameters, Potential,
};
use std::f64::consts::PI;

const EPS: f64 = f64::EPSILON;
/// GAFF's and AMBER's 1-4 scales, the van der Waals and the Coulomb: used here as two numbers that
/// are not 1.
const ONE_FOUR: [f64; 2] = [0.5, 1.0 / 1.2];
/// The small boxes: one octanol along its axis and four by four across, 432 atoms.
const SMALL: [usize; 3] = [1, 4, 4];
const SMALL_CUTOFF: f64 = 6.5 * ANGSTROM;
const KELVIN: f64 = 298.15;

/// `field` with every rule on: [`invented`] values for the atoms of `c` from atom `start`,
/// Lorentz–Berthelot, and [`ONE_FOUR`].
fn with_rules(field: PeriodicForceField, c: &Component, start: usize) -> PeriodicForceField {
    let n = c.atoms().len();
    field
        .with_nonbonded(start..start + n, invented(c))
        .expect("invented values are valid")
        .with_mixing(Mixing::LorentzBerthelot)
        .with_one_four(ONE_FOUR[0], ONE_FOUR[1])
}

/// `x` to the whole power `n` by repeated multiplication, `powi` being platform-dependent.
fn pow(x: f64, n: u32) -> f64 {
    (0..n).fold(1.0, |p, _| p * x)
}

/// `2^(−1/6)`: σ over the minimum's distance, by the kernel's `exp` and `ln`.
fn sigma_over_x() -> f64 {
    exp(-ln(2.0) / 6.0)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn len(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// The magnitudes the energy is a sum of, for its rounding.
fn magnitude(e: &PeriodicEnergy) -> f64 {
    let w = &e.ewald;
    e.bond.abs()
        + e.angle.abs()
        + e.torsion.abs()
        + e.inversion.abs()
        + e.van_der_waals.abs()
        + e.dispersion_correction.abs()
        + w.real.abs()
        + w.reciprocal.abs()
        + w.self_energy.abs()
        + w.excluded.abs()
        + w.background.abs()
}

fn largest_force(f: &[[f64; 3]]) -> f64 {
    f.iter().flatten().fold(0.0f64, |m, x| m.max(x.abs()))
}

/// A central difference and the tolerance it earns: three times the truncation term is
/// `|D(2h) − D(h)|`, and the energy's rounding, at most `64 ε` of the sum of its parts'
/// magnitudes, divided by `h`. W1's.
fn finite_difference(energy: impl Fn(f64) -> (f64, f64), h: f64) -> (f64, f64) {
    let d = |s: f64| {
        let (plus, scale) = energy(s);
        let (minus, _) = energy(-s);
        ((plus - minus) / (2.0 * s), scale)
    };
    let (d1, scale) = d(h);
    let (d2, _) = d(2.0 * h);
    (d1, 2.0 * (d2 - d1).abs() / 3.0 + 64.0 * EPS * scale / h)
}

/// The smallest distance from the cutoff, as a fraction of it, of any pair with an atom `only`
/// accepts, by minimum image: a difference that moves a pair across it sees the step.
fn nearest_to_cutoff(
    field: &PeriodicForceField,
    at: &[[f64; 3]],
    only: impl Fn(usize) -> bool,
) -> f64 {
    let cell = field.cell();
    let rc = field.ewald().parameters().cutoff;
    let mut nearest = f64::INFINITY;
    for i in 0..at.len() {
        for j in i + 1..at.len() {
            if only(i) || only(j) {
                let d = cell.minimum_image(sub(at[i], at[j]));
                nearest = nearest.min((len(d) - rc).abs() / rc);
            }
        }
    }
    nearest
}

type Bits = (Vec<u64>, Vec<[u64; 3]>, [[u64; 3]; 3]);

/// An evaluation as bits — every energy term, every force, the virial — so that `−0.0` is not
/// `0.0`, as `==` would have it.
fn evaluation_bits(e: &PeriodicEvaluation) -> Bits {
    let (t, w) = (&e.energy, &e.energy.ewald);
    let terms = [
        t.bond,
        t.angle,
        t.torsion,
        t.inversion,
        t.van_der_waals,
        t.dispersion_correction,
        t.electrostatic,
        t.total,
        w.real,
        w.reciprocal,
        w.self_energy,
        w.excluded,
        w.background,
        w.total,
    ];
    (
        terms.map(f64::to_bits).to_vec(),
        e.forces.iter().map(|f| f.map(f64::to_bits)).collect(),
        e.virial.map(|r| r.map(f64::to_bits)),
    )
}

fn small_box(o: &Octanol, seed: u64, accuracy: f64) -> Liquid {
    Liquid::lattice(o, SMALL, DENSITY, seed, SMALL_CUTOFF, accuracy)
}

/// The mesh `PmeParameters::for_accuracy` chooses for `field` at δ = 10⁻⁵.
fn mesh_for(field: &PeriodicForceField) -> PmeParameters {
    PmeParameters::for_accuracy(
        &field.cell(),
        &field.ewald().parameters(),
        1e-5,
        field.charges().len(),
    )
}

// ---------------------------------------------------------------------------------------------
// The defaults

/// **Every rule at its default is the plain field's bits.** On a box of sixteen octanols with
/// QEq's charges: the UFF types' own values written over every atom in two runs that meet at atom
/// 197 — mid-molecule, between a carbon and an oxygen, so a run off by one moves a value onto an atom of another type — geometric
/// mixing and 1-4 scales of `(1, 1)` set explicitly. The evaluation's every term, force and
/// virial component, the long-range corrections, and a decoupling of one octanol's energies,
/// forces and couplings at three states, all to the bit, through classical Ewald, the mesh, and
/// the neighbour list.
#[test]
fn every_rule_at_its_default_is_the_plain_fields_bits() {
    let o = Octanol::new(None);
    let l = small_box(&o, 0xB175, 1e-5);
    let n = l.at.len();
    let values = uff_values(&l.types);
    let plain = l.field.clone();
    let same = plain
        .clone()
        .with_nonbonded(0..197, values[..197].to_vec())
        .expect("valid")
        .with_nonbonded(197..n, values[197..].to_vec())
        .expect("valid")
        .with_mixing(Mixing::Geometric)
        .with_one_four(1.0, 1.0);
    assert_eq!(same.mixing(), Mixing::Geometric);
    assert_eq!(same.one_four(), [1.0, 1.0]);
    assert!(same.one_four_pairs().is_empty());
    let bits =
        |v: &[[f64; 2]]| -> Vec<[u64; 2]> { v.iter().map(|x| x.map(f64::to_bits)).collect() };
    assert_eq!(bits(same.van_der_waals()), bits(plain.van_der_waals()));
    // The two runs meet inside a molecule, between atoms of two types.
    assert_ne!(l.types[196], l.types[197]);
    let mesh = mesh_for(&plain);
    let group: Vec<bool> = (0..n).map(|i| i < 27).collect();
    let states = [
        Lambda::COUPLED,
        Lambda::new(0.0, 0.4, 0.7),
        Lambda::new(0.0, 0.0, 0.2),
    ];
    for (name, a, b) in [
        ("classical Ewald", plain.clone(), same.clone()),
        (
            "the mesh",
            plain.clone().with_mesh(mesh),
            same.clone().with_mesh(mesh),
        ),
        (
            "the neighbour list",
            plain
                .clone()
                .with_neighbour_list(PeriodicForceField::NEIGHBOUR_SKIN),
            same.clone()
                .with_neighbour_list(PeriodicForceField::NEIGHBOUR_SKIN),
        ),
    ] {
        assert_eq!(
            evaluation_bits(&a.evaluate(&l.at)),
            evaluation_bits(&b.evaluate(&l.at)),
            "{name}"
        );
        assert_eq!(
            a.dispersion_correction().to_bits(),
            b.dispersion_correction().to_bits()
        );
        assert_eq!(
            a.pressure_correction().to_bits(),
            b.pressure_correction().to_bits()
        );
        let (da, db) = (
            PeriodicDecoupling::new(&a, &group).expect("one octanol decouples"),
            PeriodicDecoupling::new(&b, &group).expect("one octanol decouples"),
        );
        assert_eq!(
            da.cross_dispersion_correction().to_bits(),
            db.cross_dispersion_correction().to_bits()
        );
        for l_ in states {
            let (mut fa, mut fb) = (vec![[0.0; 3]; n], vec![[0.0; 3]; n]);
            let ea = da.energy_and_forces(&l.at, l_, &mut fa);
            let eb = db.energy_and_forces(&l.at, l_, &mut fb);
            assert_eq!(ea.to_bits(), eb.to_bits(), "{name} {l_:?}");
            let fbits = |f: &[[f64; 3]]| -> Vec<[u64; 3]> {
                f.iter().map(|x| x.map(f64::to_bits)).collect()
            };
            assert_eq!(fbits(&fa), fbits(&fb), "{name} {l_:?}");
            let (ca, cb) = (da.coupling(&l.at, l_), db.coupling(&l.at, l_));
            assert_eq!(ca.energy.to_bits(), cb.energy.to_bits());
            assert_eq!(ca.gradient.map(f64::to_bits), cb.gradient.map(f64::to_bits));
        }
        assert!(da.own_scales().iter().all(|s| *s == [1.0, 1.0]));
        println!("{name}: the same bits");
    }
    assert_eq!(same, plain, "the same force field");
}

// ---------------------------------------------------------------------------------------------
// One pair, and the tail

/// Two TIP3P waters in a 30 Å cube with their charges off, the oxygens `r` apart along x with
/// `vdw[0]` and `vdw[1]`, every hydrogen `[1 Å, 0]`: so the oxygens' pair is the only van der
/// Waals there is.
fn two_waters(r: f64, vdw: [[f64; 2]; 2], mixing: Mixing) -> (PeriodicForceField, Vec<[f64; 3]>) {
    let cell = PeriodicBox::cubic(30.0 * ANGSTROM);
    let p = EwaldParameters::for_accuracy(&cell, 9.0 * ANGSTROM, 1e-5);
    let none = [1.0 * ANGSTROM, 0.0];
    let field = PeriodicForceField::tip3p(cell, p, 2)
        .with_charges(vec![0.0; 6])
        .with_nonbonded(0..6, vec![vdw[0], none, none, vdw[1], none, none])
        .expect("valid")
        .with_mixing(mixing);
    let half = 0.5 * hoh_angle();
    let water = |o: [f64; 3], away: f64| {
        [
            o,
            [
                o[0] + away * OH_LENGTH * half.cos(),
                o[1] + OH_LENGTH * half.sin(),
                o[2],
            ],
            [
                o[0] + away * OH_LENGTH * half.cos(),
                o[1] - OH_LENGTH * half.sin(),
                o[2],
            ],
        ]
    };
    let (a, y) = (10.0 * ANGSTROM, 15.0 * ANGSTROM);
    let mut at = water([a, y, y], -1.0).to_vec();
    at.extend(water([a + r, y, y], 1.0));
    (field, at)
}

/// **One pair is Lorentz–Berthelot's formula**, written in σ and ε: `σ_ij = (σ_i + σ_j)/2` with
/// `σ = x · 2^(−1/6)`, `ε_ij = √(ε_i ε_j)`, `E = 4ε [(σ/r)¹² − (σ/r)⁶]`, and the force on the first
/// oxygen `dE/dr` along +x, on the second the opposite, on every hydrogen none. Two oxygens of
/// 3.6 Å, 0.15 kcal/mol and 3.0 Å, 0.30 kcal/mol, at 3.0, 3.35, 4.0 and 6.0 Å. Each within
/// `64 ε` of the formula's two terms — `σ` by the kernel's `exp` and `ln` to about an ulp, raised to
/// the twelfth power — and **the geometric rule's field is the geometric formula** to the same, and
/// differs from it by more than a thousand times that bound at every distance.
#[test]
fn one_pair_is_lorentz_berthelots_formula() {
    let kcal = KCAL_PER_MOL;
    let vdw = [[3.6 * ANGSTROM, 0.15 * kcal], [3.0 * ANGSTROM, 0.30 * kcal]];
    let s = sigma_over_x();
    let (s1, s2) = (vdw[0][0] * s, vdw[1][0] * s);
    let eps = (vdw[0][1] * vdw[1][1]).sqrt();
    let formula = |sigma: f64, r: f64| {
        let q = sigma / r;
        let (q6, q12) = (pow(q, 6), pow(q, 12));
        let energy = 4.0 * eps * (q12 - q6);
        let de_dr = 4.0 * eps * (-12.0 * q12 + 6.0 * q6) / r;
        (energy, de_dr, 4.0 * eps * (q12 + q6))
    };
    for r in [3.0, 3.35, 4.0, 6.0] {
        let r = r * ANGSTROM;
        let mut values = Vec::new();
        for (mixing, sigma) in [
            (Mixing::LorentzBerthelot, 0.5 * (s1 + s2)),
            (Mixing::Geometric, (s1 * s2).sqrt()),
        ] {
            let (field, at) = two_waters(r, vdw, mixing);
            let ev = field.evaluate(&at);
            let (e, de_dr, size) = formula(sigma, r);
            let bound = 64.0 * EPS * size;
            println!(
                "{mixing:?} at {:.2} Å: {:+.12e} kcal/mol, the formula {:+.12e}; force {:+.6e} N \
                 against {:+.6e}",
                r / ANGSTROM,
                ev.energy.van_der_waals / kcal,
                e / kcal,
                ev.forces[0][0],
                de_dr
            );
            assert!((ev.energy.van_der_waals - e).abs() <= bound, "{mixing:?}");
            let force_bound = 64.0 * EPS * size * 12.0 / r;
            assert!((ev.forces[0][0] - de_dr).abs() <= force_bound);
            assert!((ev.forces[3][0] + de_dr).abs() <= force_bound);
            for i in [1, 2, 4, 5] {
                assert_eq!(ev.forces[i], [0.0; 3], "hydrogen {i} has no van der Waals");
            }
            for i in [0, 3] {
                assert!(
                    ev.forces[i][1].abs() <= force_bound && ev.forces[i][2].abs() <= force_bound
                );
            }
            values.push((e, bound));
        }
        let (lb, geometric) = (values[0], values[1]);
        assert!(
            (lb.0 - geometric.0).abs() > 1000.0 * lb.1,
            "the two rules part at {r:e}"
        );
    }
}

/// **The long-range correction under Lorentz–Berthelot is its double sum over every ordered pair**
/// of the box's 432 atoms, written out — six classes, so the class sum is not one class's square,
/// and among them one that shares a carbon's `D` and not its `x` and one its `x` and not its `D`,
/// so that classes keyed on either alone would merge two and be caught (GAFF's `hc` and `h1`
/// share their `ε`) — energy and pressure; and the geometric rule's, with the same overrides, is its own double sum.
/// Each within `(N² + 64) ε` of the sum of its two parts' magnitudes: every term of each part is
/// positive, so the rounding of a sum of `N²` of them in any order is below `N² ε` of it. The
/// decoupling of one octanol's cross part is its own double sum, both ways round. And the rules
/// part by more than a thousand times that bound.
#[test]
fn the_tail_under_lorentz_berthelot_is_its_double_sum() {
    let o = Octanol::new(None);
    let l = small_box(&o, 0x7A1, 1e-5);
    let n = l.at.len();
    let rc = l.field.ewald().parameters().cutoff;
    let v = l.cell().volume();
    let bound = |size: f64, terms: usize| (terms as f64 + 64.0) * EPS * size;
    let mut tails = Vec::new();
    for mixing in [Mixing::LorentzBerthelot, Mixing::Geometric] {
        let carbon = invented(&o.component)[0];
        assert_eq!(o.component.atoms()[0].element, Element::C);
        // The second octanol's eight carbons share the carbon's D, the third's its x.
        let field = with_rules(l.field.clone(), &l.component, 0)
            .with_mixing(mixing)
            .with_nonbonded(27..35, vec![[5.0 * ANGSTROM, carbon[1]]; 8])
            .expect("valid")
            .with_nonbonded(54..62, vec![[carbon[0], 0.05 * KCAL_PER_MOL]; 8])
            .expect("valid");
        let vdw = field.van_der_waals().to_vec();
        let mut classes: Vec<[u64; 2]> = vdw.iter().map(|x| x.map(f64::to_bits)).collect();
        classes.sort_unstable();
        classes.dedup();
        assert_eq!(
            classes.len(),
            6,
            "C, O, H on C, H on O, and the two that share half a carbon's"
        );
        let ([e, p], [se, sp]) = tail_by_pairs(&vdw, &vdw, mixing, rc, v);
        println!(
            "{mixing:?}: E_tail {:+.12e} J against the double sum {e:+.12e}; P_tail {:+.12e} Pa \
             against {p:+.12e}",
            field.dispersion_correction(),
            field.pressure_correction()
        );
        assert!((field.dispersion_correction() - e).abs() <= bound(se, n * n));
        assert!((field.pressure_correction() - p).abs() <= bound(sp, n * n));
        let group: Vec<bool> = (0..n).map(|i| i < 27).collect();
        let d = PeriodicDecoupling::new(&field, &group).expect("one octanol decouples");
        let ([cross, _], [size, _]) = tail_by_pairs(&vdw[..27], &vdw[27..], mixing, rc, v);
        let theirs = d.cross_dispersion_correction();
        println!(
            "  the cross part {theirs:+.12e} J against {:+.12e}",
            2.0 * cross
        );
        assert!((theirs - 2.0 * cross).abs() <= bound(2.0 * size, 27 * (n - 27)));
        tails.push((field.dispersion_correction(), bound(se, n * n)));
    }
    assert!((tails[0].0 - tails[1].0).abs() > 1000.0 * tails[0].1);
}

/// **With one class the tail is Allen and Tildesley's printed form**,
/// `(8/3)πNρεσ³[(σ/r_c)⁹/3 − (σ/r_c)³]` and `(16/3)πρ²εσ³[(2/3)(σ/r_c)⁹ − (σ/r_c)³]`, under
/// Lorentz–Berthelot: 81 atoms of TIP3P's box all given 3.3 Å and 0.12 kcal/mol. Within 128 ε:
/// `σ` to about three ulps through `exp` and `ln`, raised to the ninth power.
#[test]
fn one_class_is_allen_and_tildesleys_form() {
    let water = pantometry_forcefield::WaterBox::lattice(3, 0x0C1);
    let n = water.positions().len();
    let one = [3.3 * ANGSTROM, 0.12 * KCAL_PER_MOL];
    let field = water
        .force_field(4.6 * ANGSTROM, 1e-5)
        .with_nonbonded(0..n, vec![one; n])
        .expect("valid")
        .with_mixing(Mixing::LorentzBerthelot);
    let rc = field.ewald().parameters().cutoff;
    let (sigma, eps) = (one[0] * sigma_over_x(), one[1]);
    let rho = n as f64 / field.cell().volume();
    let s3 = pow(sigma / rc, 3);
    let e_at = 8.0 / 3.0 * PI * n as f64 * rho * eps * pow(sigma, 3) * (s3 * s3 * s3 / 3.0 - s3);
    let p_at = 16.0 / 3.0 * PI * rho * rho * eps * pow(sigma, 3) * (2.0 / 3.0 * s3 * s3 * s3 - s3);
    println!(
        "one class: E_tail {:.12e} against A&T's {e_at:.12e}; P_tail {:.12e} against {p_at:.12e}",
        field.dispersion_correction(),
        field.pressure_correction()
    );
    assert!((field.dispersion_correction() / e_at - 1.0).abs() <= 128.0 * EPS);
    assert!((field.pressure_correction() / p_at - 1.0).abs() <= 128.0 * EPS);
}

// ---------------------------------------------------------------------------------------------
// The 1-4 pairs

/// The cycloalkane `CₙH₂ₙ`: positions only to type it.
fn cycloalkane(n: usize) -> Component {
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let names: Vec<[String; 3]> = (0..n)
        .map(|k| [format!("C{k}"), format!("H{k}A"), format!("H{k}B")])
        .collect();
    for (k, nm) in names.iter().enumerate() {
        let t = 2.0 * PI * k as f64 / n as f64;
        let r = 0.75 * n as f64 / 3.0 + 0.5;
        let c = [r * t.cos(), r * t.sin(), 0.0];
        atoms.push((nm[0].as_str(), "C", c));
        let h = r + 0.95;
        atoms.push((nm[1].as_str(), "H", [h * t.cos(), h * t.sin(), 0.6]));
        atoms.push((nm[2].as_str(), "H", [h * t.cos(), h * t.sin(), -0.6]));
    }
    for k in 0..n {
        let next = &names[(k + 1) % n][0];
        bonds.push((names[k][0].as_str(), next.as_str(), "SING", false));
        bonds.push((names[k][0].as_str(), names[k][1].as_str(), "SING", false));
        bonds.push((names[k][0].as_str(), names[k][2].as_str(), "SING", false));
    }
    common::entry(&atoms, &bonds)
}

/// **The scaled pairs are the pairs three bonds apart by breadth-first search**, and nothing that
/// is excluded: octanol's — acyclic, so they are its `Σ_bonds (deg_j − 1)(deg_k − 1)` torsions'
/// ends, 66 — benzene's 21, whose para carbons are three bonds apart both ways round and appear
/// once; and the small rings, counted by hand. **Cyclopropane, 12**: every C–C is a bond and
/// every C–H 1-3, so only the hydrogens on two carbons, `3 × 2 × 2` — and a path round the
/// triangle comes back to its own start, which must be no pair. **Cyclobutane, 24**: each
/// hydrogen and the carbon across the ring, 8, and the hydrogens on neighbouring carbons,
/// `4 × 2 × 2`; the carbons across are 1-3 by two paths. **Cyclopentane, 40**, all with a
/// hydrogen: every two of its carbons are 1-2 or 1-3, so exclusion comes first and the ring has
/// no 1-4 pair across it.
#[test]
fn the_one_four_pairs_are_three_bonds_apart() {
    let octanol = Component::from_ccd(octanol::OC9).expect("OC9");
    let benzene = Component::from_ccd(hydration::BNZ).expect("BNZ");
    let (three, four, ring) = (cycloalkane(3), cycloalkane(4), cycloalkane(5));
    let cell = PeriodicBox::cubic(30.0 * ANGSTROM);
    let p = EwaldParameters::for_accuracy(&cell, 9.0 * ANGSTROM, 1e-5);
    for (name, c, want) in [
        ("octanol", &octanol, 66),
        ("benzene", &benzene, 21),
        ("cyclopropane", &three, 12),
        ("cyclobutane", &four, 24),
        ("cyclopentane", &ring, 40),
    ] {
        let n = c.atoms().len();
        let b = bonds(c);
        let field = PeriodicForceField::new(c, &uff::assign(c), cell, p)
            .expect("supported")
            .with_one_four(ONE_FOUR[0], ONE_FOUR[1]);
        let search = one_four_by_search(n, &b);
        assert_eq!(field.one_four_pairs(), &search[..], "{name}");
        assert_eq!(search.len(), want, "{name}");
        for pair in field.one_four_pairs() {
            assert!(field.excluded_pairs().binary_search(pair).is_err());
        }
        let d = graph_distances(n, &b);
        let excluded = (0..n)
            .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
            .filter(|&(i, j)| matches!(d[i][j], Some(1) | Some(2)))
            .count();
        assert_eq!(field.excluded_pairs().len(), excluded, "{name}");
        println!("{name}: {} 1-4 pairs, {} excluded", search.len(), excluded);
    }
    let b = bonds(&octanol);
    let mut degree = vec![0usize; octanol.atoms().len()];
    for &[i, j] in &b {
        degree[i] += 1;
        degree[j] += 1;
    }
    let torsions: usize = b
        .iter()
        .map(|&[j, k]| (degree[j] - 1) * (degree[k] - 1))
        .sum();
    assert_eq!(torsions, 66);
    let ring_pairs = PeriodicForceField::new(&ring, &uff::assign(&ring), cell, p)
        .unwrap()
        .with_one_four(0.5, 0.5);
    assert!(ring_pairs
        .one_four_pairs()
        .iter()
        .all(|&[i, j]| ring.atoms()[i].element == Element::H
            || ring.atoms()[j].element == Element::H));
    // Cyclopropane's are all H–H, cyclobutane's 8 C–H and 16 H–H.
    let kinds = |c: &Component| {
        let f = PeriodicForceField::new(c, &uff::assign(c), cell, p)
            .unwrap()
            .with_one_four(0.5, 0.5);
        let h = |i: usize| c.atoms()[i].element == Element::H;
        let both = f
            .one_four_pairs()
            .iter()
            .filter(|&&[i, j]| h(i) && h(j))
            .count();
        (f.one_four_pairs().len() - both, both)
    };
    assert_eq!(kinds(&three), (0, 12));
    assert_eq!(kinds(&four), (8, 16));
}

/// **A wrong override is refused by name**: a run past the atoms or backwards, values not one per
/// atom, a value that is not finite, an `x` that is negative or zero with a well, a `D` that is
/// negative. `[0, 0]`, no van der Waals at all, is taken, as GAFF's `ho` is written.
#[test]
fn a_wrong_override_is_refused_by_name() {
    let cell = PeriodicBox::cubic(30.0 * ANGSTROM);
    let p = EwaldParameters::for_accuracy(&cell, 9.0 * ANGSTROM, 1e-5);
    let field = PeriodicForceField::tip3p(cell, p, 2);
    let good = [3.0 * ANGSTROM, 1e-21];
    let refused = |atoms: std::ops::Range<usize>, values: Vec<[f64; 2]>| {
        field.clone().with_nonbonded(atoms, values).unwrap_err()
    };
    assert_eq!(
        refused(5..7, vec![good; 2]),
        NonbondedError::OutOfRange {
            start: 5,
            end: 7,
            atoms: 6
        }
    );
    #[allow(clippy::reversed_empty_ranges)]
    let backwards = refused(4..2, vec![]);
    assert!(matches!(
        backwards,
        NonbondedError::OutOfRange {
            start: 4,
            end: 2,
            ..
        }
    ));
    assert_eq!(
        refused(0..3, vec![good; 2]),
        NonbondedError::Lengths {
            atoms: 3,
            values: 2
        }
    );
    assert_eq!(
        refused(2..4, vec![good, [f64::NAN, 1e-21]]),
        NonbondedError::NotFinite { atom: 3 }
    );
    assert_eq!(
        refused(2..4, vec![good, [3e-10, f64::INFINITY]]),
        NonbondedError::NotFinite { atom: 3 }
    );
    assert_eq!(
        refused(1..2, vec![[0.0, 1e-21]]),
        NonbondedError::Distance {
            atom: 1,
            value: 0.0
        }
    );
    assert_eq!(
        refused(1..2, vec![[3e-10, -1e-22]]),
        NonbondedError::Well {
            atom: 1,
            value: -1e-22
        }
    );
    assert_eq!(
        refused(1..2, vec![[-1e-10, 0.0]]),
        NonbondedError::Distance {
            atom: 1,
            value: -1e-10
        }
    );
    let none = field
        .clone()
        .with_nonbonded(1..2, vec![[0.0, 0.0]])
        .expect("no van der Waals");
    assert_eq!(none.van_der_waals()[1], [0.0, 0.0]);
    for e in [
        refused(5..7, vec![good; 2]),
        refused(0..3, vec![good; 2]),
        refused(1..2, vec![[0.0, 1e-21]]),
    ] {
        println!("{e}");
        assert!(!e.to_string().is_empty());
    }
    // An empty run is no change; the last atom is in range.
    assert_eq!(field.clone().with_nonbonded(6..6, vec![]).unwrap(), field);
    assert!(field.clone().with_nonbonded(5..6, vec![good]).is_ok());
}

#[test]
#[should_panic(expected = "a 1-4 scale must be finite and not negative")]
fn a_negative_one_four_scale_is_refused() {
    let cell = PeriodicBox::cubic(30.0 * ANGSTROM);
    let p = EwaldParameters::for_accuracy(&cell, 9.0 * ANGSTROM, 1e-5);
    let _ = PeriodicForceField::tip3p(cell, p, 2).with_one_four(-0.5, 1.0);
}

// ---------------------------------------------------------------------------------------------
// Derivatives

/// **The forces are `−∇E` and the virial `−∂U/∂ε`**, with every rule on, on sixteen octanols with
/// QEq's charges: six atoms of two molecules — the hydroxyl's oxygen and hydrogen, two carbons and a
/// hydrogen of one, a carbon of the other — each along three axes, and the three diagonal strains,
/// by central differences to the tolerance each earns; the virial symmetric.
fn derivatives(mesh: bool) {
    let o = Octanol::new(None);
    let l = small_box(&o, 0xF0D, 1e-6);
    let mut field = with_rules(l.field.clone(), &l.component, 0);
    if mesh {
        let m = mesh_for(&field);
        field = field.with_mesh(m);
    }
    assert_eq!(field.ewald().mesh().is_some(), mesh);
    assert_eq!(field.one_four_pairs().len(), 66 * 16);
    let at = l.at.clone();
    let ev = field.evaluate(&at);
    let c = &o.component;
    let oxygen = index(c, "OAB");
    let atoms = [
        oxygen,
        index(c, "HAB"),
        index(c, "CAD"),
        index(c, "CAG"),
        index(c, "HAE1"),
        5 * 27 + index(c, "CAC"),
    ];
    let h = 1e-5 * ANGSTROM;
    let rc = field.ewald().parameters().cutoff;
    let mut worst = 0.0f64;
    for &i in &atoms {
        assert!(nearest_to_cutoff(&field, &at, |j| j == i) * rc > 4.0 * h);
        let fmax = ev.forces[i].iter().fold(0.0f64, |m, f| m.max(f.abs()));
        for a in 0..3 {
            let energy = |s: f64| {
                let mut moved = at.clone();
                moved[i][a] += s;
                let e = field.energy(&moved);
                (e.total, magnitude(&e))
            };
            let (d, tolerance) = finite_difference(energy, h);
            let f = ev.forces[i][a];
            assert!(
                (f + d).abs() <= tolerance,
                "atom {i} axis {a}: {f:e} against {:e} ± {tolerance:e}",
                -d
            );
            assert!(tolerance < 1e-4 * fmax, "{tolerance:e} against {fmax:e}");
            worst = worst.max((f + d).abs() / tolerance);
        }
    }
    println!("mesh {mesh}: the worst |F + dE/dx| is {worst:.3} of its tolerance");
    let hs = 1e-6;
    assert!(nearest_to_cutoff(&field, &at, |_| true) > 4.0 * hs);
    let w = ev.virial;
    for a in 0..3 {
        let energy = |s: f64| {
            let mut f = [1.0; 3];
            f[a] += s;
            let scaled: Vec<[f64; 3]> = at
                .iter()
                .map(|r| [r[0] * f[0], r[1] * f[1], r[2] * f[2]])
                .collect();
            let e = field.with_cell(field.cell().scaled(f)).energy(&scaled);
            (e.total, magnitude(&e))
        };
        let (d, tolerance) = finite_difference(energy, hs);
        println!(
            "mesh {mesh}: W_{a}{a} = {:.8e} J, −dU/dε = {:.8e} ± {tolerance:.1e}",
            w[a][a], -d
        );
        assert!((w[a][a] + d).abs() <= tolerance, "axis {a}");
        assert!(tolerance < 1e-5 * w[a][a].abs());
    }
    let wmax = w.iter().flatten().fold(0.0f64, |m, x| m.max(x.abs()));
    for a in 0..3 {
        for b in 0..3 {
            assert!((w[a][b] - w[b][a]).abs() <= 1e-10 * wmax);
        }
    }
}

#[test]
fn forces_and_virial_are_the_energys_derivatives_through_ewald() {
    derivatives(false);
}

#[test]
fn forces_and_virial_are_the_energys_derivatives_through_the_mesh() {
    derivatives(true);
}

// ---------------------------------------------------------------------------------------------
// One molecule in a large box

/// **A 1-4 pair's scaled Coulomb through Ewald, and its scaled van der Waals**: one octanol with
/// QEq's charges, every rule on, alone in cubes of 40, 60 and 80 Å with `r_c = L/2` — past the
/// molecule's 11 Å, so every pair of it is inside the cutoff and no image is — at δ = 10⁻¹².
/// Its van der Waals is [`Vacuum`]'s direct sum over every pair, the 1-4 pairs at `s_v`, within
/// the sum's rounding; its bonded terms are UFF's to the bit; and its electrostatics are
/// [`Vacuum`]'s Coulomb with the 1-4 pairs at `s_c`, less the tinfoil term `2π k_e μ²/3V`, with
/// the force `(4π k_e/3V) q_i μ` beside it — and what is left falls as `L⁻⁵`, energy and forces,
/// within a quarter from one box to the next. The 1-4 pairs' share of the vacuum Coulomb,
/// `(1 − s_c)` of it, is over a thousand times what is left in the smallest box, so a reciprocal
/// sum not corrected for them would be seen.
#[test]
fn one_octanol_in_a_large_box_is_its_vacuum_energy_with_the_one_four_pairs_scaled() {
    let o = Octanol::new(None);
    let c = &o.component;
    let n = c.atoms().len();
    let q = o.charges.clone();
    let at = o.at.clone();
    let vacuum = Vacuum::new(
        c,
        &o.types,
        invented(c),
        q.clone(),
        Mixing::LorentzBerthelot,
        ONE_FOUR,
    );
    let full = Vacuum::new(
        c,
        &o.types,
        invented(c),
        q.clone(),
        Mixing::LorentzBerthelot,
        [ONE_FOUR[0], 1.0],
    );
    let one_four_share = full.nonbonded(&at).coulomb - vacuum.nonbonded(&at).coulomb;
    let mu = [0, 1, 2].map(|a| (0..n).map(|i| q[i] * at[i][a]).sum::<f64>());
    let mu2 = mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2];
    let span = at
        .iter()
        .flat_map(|p| at.iter().map(move |r| len(sub(*p, *r))))
        .fold(0.0f64, f64::max);
    assert!(span < 12.0 * ANGSTROM, "span {span:e}");
    let mut rows = Vec::new();
    for side in [40.0, 60.0, 80.0] {
        let cell = PeriodicBox::cubic(side * ANGSTROM);
        let p = EwaldParameters::for_accuracy(&cell, 0.5 * side * ANGSTROM, 1e-12);
        let field = with_rules(
            PeriodicForceField::new(c, &o.types, cell, p)
                .expect("octanol")
                .with_charges(q.clone()),
            c,
            0,
        );
        // Centred, so that no atom is near a face.
        let centred: Vec<[f64; 3]> = at
            .iter()
            .map(|r| r.map(|x| x + 0.5 * side * ANGSTROM))
            .collect();
        let e = field.evaluate(&centred);
        let nb = vacuum.nonbonded(&centred);
        let bonded = vacuum.bonded.evaluate(&centred);
        assert_eq!(e.energy.bond, bonded.energy.bond);
        assert_eq!(e.energy.angle, bonded.energy.angle);
        assert_eq!(e.energy.torsion, bonded.energy.torsion);
        assert_eq!(e.energy.inversion, bonded.energy.inversion);
        let vdw_bound = nb.pairs as f64 * EPS * nb.magnitude;
        assert!(
            (e.energy.van_der_waals - nb.van_der_waals).abs() <= vdw_bound,
            "van der Waals {} against {}",
            e.energy.van_der_waals,
            nb.van_der_waals
        );
        let volume = cell.volume();
        let dipole = 2.0 * PI * COULOMB * mu2 / (3.0 * volume);
        let left = e.energy.electrostatic + dipole - nb.coulomb;
        let strength = 4.0 * PI * COULOMB / (3.0 * volume);
        let (mut off, mut tinfoil_max) = (0.0f64, 0.0f64);
        for i in 0..n {
            for a in 0..3 {
                let tinfoil = strength * q[i] * mu[a];
                let want = bonded.forces[i][a] + nb.forces[i][a] + tinfoil;
                off = off.max((e.forces[i][a] - want).abs());
                tinfoil_max = tinfoil_max.max(tinfoil.abs());
            }
        }
        println!(
            "octanol in {side} Å: Ewald − vacuum {:+.4e} kcal/mol, −2πμ²/3V {:+.4e}, left {:+.3e}; \
             forces: the dipole field's largest {:.3e} N, left {:.3e}; the 1-4 share {:+.4e}",
            (e.energy.electrostatic - nb.coulomb) / KCAL_PER_MOL,
            -dipole / KCAL_PER_MOL,
            left / KCAL_PER_MOL,
            tinfoil_max,
            off,
            one_four_share / KCAL_PER_MOL
        );
        rows.push((side, left, dipole, off, tinfoil_max));
    }
    assert!(one_four_share.abs() > 1000.0 * rows[0].1.abs());
    for w in rows.windows(2) {
        let (a, b) = (w[0], w[1]);
        let expected = pow(b.0 / a.0, 5);
        let (energy, force) = (a.1 / b.1, a.3 / b.3);
        println!(
            "{} → {} Å: left falls by {energy:.3}, forces' by {force:.3}, L⁻⁵ {expected:.3}",
            a.0, b.0
        );
        assert!((energy / expected - 1.0).abs() < 0.25, "energy {energy}");
        assert!((force / expected - 1.0).abs() < 0.25, "forces {force}");
        assert!(
            a.1.abs() < 0.25 * a.2 && a.3 < 0.25 * a.4,
            "the dipole term dominates"
        );
    }
}

// ---------------------------------------------------------------------------------------------
// Decoupling

/// Benzene (QEq, relaxed) at the centre of a [`SMALL`] octanol lattice less every octanol within
/// 2 Å of it, every rule on everywhere; and the systems its end states should be: benzene alone
/// in the box and in vacuum, the octanols alone in the box, each with the same rules.
struct Decoupled {
    field: PeriodicForceField,
    decoupling: PeriodicDecoupling,
    benzene_in_the_box: PeriodicForceField,
    benzene_in_vacuum: Vacuum,
    octanol_alone: PeriodicForceField,
    at: Vec<[f64; 3]>,
}

fn decoupled(mesh: bool) -> Decoupled {
    let b = hydration::Benzene::new(None);
    let o = Octanol::new(None);
    let (cell, _, _) = Liquid::cell_of(&o, SMALL, DENSITY);
    let centre = cell.lengths().map(|x| 0.5 * x);
    let l = Liquid::around(
        Solute {
            component: &b.component,
            types: &b.types,
            charges: &b.charges,
            at: b.centred_at(centre),
            clearance: 2.0 * ANGSTROM,
        },
        &o,
        SMALL,
        DENSITY,
        0xE5D,
        SMALL_CUTOFF,
        1e-5,
    );
    let p = l.field.ewald().parameters();
    let mut field = with_rules(l.field.clone(), &l.component, 0);
    let mut benzene_in_the_box = with_rules(
        PeriodicForceField::new(&b.component, &b.types, cell, p)
            .expect("benzene")
            .with_charges(b.charges.clone()),
        &b.component,
        0,
    );
    let parts: Vec<&Component> = (0..l.octanols).map(|_| &o.component).collect();
    let joined = octanol::joined(&parts);
    let mut octanol_alone = with_rules(
        PeriodicForceField::new(&joined, &l.types[12..], cell, p)
            .expect("octanol alone")
            .with_charges(l.charges[12..].to_vec()),
        &joined,
        0,
    );
    if mesh {
        let m = mesh_for(&field);
        field = field.with_mesh(m);
        benzene_in_the_box = benzene_in_the_box.with_mesh(m);
        octanol_alone = octanol_alone.with_mesh(m);
    }
    let group: Vec<bool> = (0..l.at.len()).map(|i| i < 12).collect();
    Decoupled {
        decoupling: PeriodicDecoupling::new(&field, &group).expect("benzene decouples"),
        field,
        benzene_in_the_box,
        benzene_in_vacuum: Vacuum::new(
            &b.component,
            &b.types,
            invented(&b.component),
            b.charges.clone(),
            Mixing::LorentzBerthelot,
            ONE_FOUR,
        ),
        octanol_alone,
        at: l.at.clone(),
    }
}

/// **The end states are the systems they say, with every rule on**: fully coupled, the whole box
/// less benzene alone in the box, plus benzene alone in vacuum by [`Vacuum`]'s direct sums — its
/// 21 1-4 pairs scaled, its pairs Lorentz–Berthelot's; fully decoupled, the octanols alone in the
/// box plus benzene in vacuum, energy and every force. Each within `10⁻¹²` of the parts'
/// magnitudes: the same terms summed in other orders. The decoupling's own pairs carry the 1-4
/// scales exactly on benzene's 21.
#[test]
fn the_end_states_are_the_systems_they_say_with_every_rule_on() {
    let d = decoupled(false);
    let n = d.at.len();
    let ones = d
        .decoupling
        .own_scales()
        .iter()
        .filter(|s| **s == ONE_FOUR)
        .count();
    assert_eq!(ones, 21);
    // Benzene's 66 pairs less its 12 bonds and 18 angles' ends.
    assert_eq!(d.decoupling.own_pairs().len(), 66 - 12 - 18);
    let mut f = vec![[0.0; 3]; n];
    let coupled = d
        .decoupling
        .energy_and_forces(&d.at, Lambda::COUPLED, &mut f);
    let whole = d.field.energy(&d.at);
    let in_box = d.benzene_in_the_box.energy(&d.at[..12]);
    let mut vf = vec![[0.0; 3]; 12];
    let vac = d.benzene_in_vacuum.energy_and_forces(&d.at[..12], &mut vf);
    let expected = whole.total - in_box.total + vac;
    let scale = magnitude(&whole);
    println!(
        "coupled {coupled:.10e} J against {expected:.10e}; benzene's own images {:+.4e} kcal/mol",
        (in_box.total - vac) / KCAL_PER_MOL
    );
    assert!((coupled - expected).abs() <= 1e-12 * scale);
    let decoupled = d
        .decoupling
        .energy_and_forces(&d.at, Lambda::new(0.0, 0.0, 0.0), &mut f);
    let solvent = d.octanol_alone.evaluate(&d.at[12..]);
    let expected = solvent.energy.total + vac;
    println!("decoupled {decoupled:.10e} J against {expected:.10e}");
    assert!((decoupled - expected).abs() <= 1e-12 * scale);
    let fscale = largest_force(&f);
    for i in 0..n {
        let want = if i < 12 {
            vf[i]
        } else {
            solvent.forces[i - 12]
        };
        for k in 0..3 {
            assert!(
                (f[i][k] - want[k]).abs() <= 1e-12 * fscale,
                "atom {i} axis {k}: {:e} against {:e}",
                f[i][k],
                want[k]
            );
        }
    }
}

/// **The λ derivatives are the whole energy's**, with every rule on, at four states — the charges
/// part-way, both part-way, the van der Waals part-way with the charges off, and near the
/// decoupled end where the soft core is softest — against central differences of
/// `energy_and_forces`, each to the tolerance it earns. **`couplings` is `coupling`'s bits** at
/// every state, and **the coupling is the evaluation's λ-dependent part**: `U(λ) − U(0)` is its
/// energy within twice the evaluations' rounding, `128 ε` of the parts.
#[test]
fn the_lambda_derivatives_are_the_whole_energys_with_every_rule_on() {
    let d = decoupled(false);
    let at = &d.at;
    let states = [
        Lambda::new(0.0, 0.5, 1.0),
        Lambda::new(0.0, 0.4, 0.8),
        Lambda::new(0.0, 0.0, 0.5),
        Lambda::new(0.0, 0.0, 0.1),
    ];
    let all = d.decoupling.couplings(at, &states);
    let scale = magnitude(&d.field.energy(at));
    let u = |l: Lambda| {
        let mut f = vec![[0.0; 3]; at.len()];
        d.decoupling.energy_and_forces(at, l, &mut f)
    };
    let u0 = u(Lambda::new(0.0, 0.0, 0.0));
    for (lambda, many) in states.iter().zip(&all) {
        let c = d.decoupling.coupling(at, *lambda);
        assert_eq!(c.energy.to_bits(), many.energy.to_bits());
        assert_eq!(
            c.gradient.map(f64::to_bits),
            many.gradient.map(f64::to_bits)
        );
        let du = u(*lambda) - u0;
        println!(
            "{lambda:?}: U(λ) − U(0) {du:+.10e} J, the coupling {:+.10e}",
            c.energy
        );
        assert!((du - c.energy).abs() <= 128.0 * EPS * scale);
        for (k, h) in [(1usize, 1e-3), (2, 1e-4)] {
            let at_k = lambda.components()[k];
            if at_k < 2.0 * h || at_k > 1.0 - 2.0 * h {
                continue;
            }
            let energy = |e: f64| {
                let mut x = *lambda;
                if k == 1 {
                    x.electrostatics += e;
                } else {
                    x.van_der_waals += e;
                }
                (u(x), scale)
            };
            let (dl, tolerance) = finite_difference(energy, h);
            println!(
                "  ∂U/∂λ[{k}] = {:+.10e} J, by differences {dl:+.10e} ± {tolerance:.1e}",
                c.gradient[k]
            );
            assert!((c.gradient[k] - dl).abs() <= tolerance, "{lambda:?} {k}");
            assert!(
                tolerance < 1e-3 * c.gradient[k].abs(),
                "the check can see a part in 10³"
            );
        }
    }
}

/// **The decoupled reciprocal energy is `rest + λ_e · cross`**, with every rule on, through
/// classical Ewald and the mesh: the field's reciprocal energy with only the rest's charges, with
/// only benzene's, and with all, and the decoupling's at three states, within `64 ε` of the parts.
/// A check that is not a difference, so it tells the rest from the group.
#[test]
fn the_decoupled_reciprocal_energy_is_the_rest_and_the_cross_terms() {
    for mesh in [false, true] {
        let d = decoupled(mesh);
        let at = &d.at;
        let q = d.field.charges().to_vec();
        let reciprocal = |keep: &dyn Fn(usize) -> bool| {
            let masked: Vec<f64> = (0..q.len())
                .map(|i| if keep(i) { q[i] } else { 0.0 })
                .collect();
            d.field
                .clone()
                .with_charges(masked)
                .energy(at)
                .ewald
                .reciprocal
        };
        let all = reciprocal(&|_| true);
        let rest = reciprocal(&|i| i >= 12);
        let own = reciprocal(&|i| i < 12);
        let scale = magnitude(&d.field.energy(at));
        for l in [
            Lambda::COUPLED,
            Lambda::new(0.0, 0.5, 1.0),
            Lambda::new(0.0, 0.0, 1.0),
        ] {
            let e = d.decoupling.energy(at, l).ewald.reciprocal;
            let want = rest + l.electrostatics * (all - rest - own);
            println!("mesh {mesh}, {l:?}: {e:+.12e} J against {want:+.12e}");
            assert!((e - want).abs() <= 64.0 * EPS * scale);
        }
        assert!((all - rest - own).abs() > 1e-6 * scale);
    }
}

// ---------------------------------------------------------------------------------------------
// Determinism

/// **The same bits however a run is cut, and with the neighbour list on or off**, with every rule
/// on: sixteen octanols, every X–H held, a 298 K bath, 12 steps of 1 fs in one call and cut as
/// 12 × 1 and 5 + 7, positions and velocities to the bit; and along that trajectory, at every
/// step, two fields with neighbour lists give the cell list's evaluation to the bit: one of skin
/// 0.05 Å, rebuilt at every step, and one of the default 2 Å, built once and reused — so both the
/// build's search, which leaves out the scaled 1-4 pairs, and the list's reuse are seen.
#[test]
fn the_same_bits_however_cut_and_with_the_list_on_or_off() {
    let o = Octanol::new(None);
    let l = small_box(&o, 0xD1CE, 1e-5);
    let field = with_rules(l.field.clone(), &l.component, 0);
    let listed = field.clone().with_neighbour_list(0.05 * ANGSTROM);
    let md = l
        .dynamics()
        .with_bath(Bath::Langevin {
            temperature: KELVIN,
            friction: 5e12,
            seed: 7,
        })
        .thermalised(&l.at, KELVIN, 7);
    let dt = FS;
    let mut whole = (md.clone(), l.at.clone());
    whole.0.run(&listed, &mut whole.1, dt, 12);
    for cuts in [vec![1usize; 12], vec![5, 7]] {
        let (mut m, mut at) = (md.clone(), l.at.clone());
        for k in cuts {
            m.run(&listed, &mut at, dt, k);
        }
        assert_eq!(at, whole.1);
        assert_eq!(m.velocities(), whole.0.velocities());
    }
    let moved = whole
        .1
        .iter()
        .zip(&l.at)
        .map(|(a, b)| len(sub(*a, *b)))
        .fold(0.0f64, f64::max);
    assert!(moved > 0.05 * ANGSTROM, "moved {moved:e}");
    // Along the trajectory, one step at a time, the listed field against the cell list.
    let rebuilt = field.clone().with_neighbour_list(0.05 * ANGSTROM);
    let reused = field
        .clone()
        .with_neighbour_list(PeriodicForceField::NEIGHBOUR_SKIN);
    let (mut m, mut at) = (md.clone(), l.at.clone());
    for step in 0..12 {
        let cells = evaluation_bits(&field.evaluate(&at));
        assert_eq!(
            evaluation_bits(&rebuilt.evaluate(&at)),
            cells,
            "step {step}"
        );
        assert_eq!(evaluation_bits(&reused.evaluate(&at)), cells, "step {step}");
        m.run(&field, &mut at, dt, 1);
    }
    let (a, b) = (
        rebuilt.neighbour_list().expect("a list"),
        reused.neighbour_list().expect("a list"),
    );
    println!(
        "skin 0.05 Å: {} builds in {} evaluations; 2 Å: {} in {}",
        a.builds, a.evaluations, b.builds, b.evaluations
    );
    assert!(a.builds >= 3, "{} builds", a.builds);
    assert_eq!((b.builds, b.evaluations), (1, 12));
    // And an evaluation is its own bits.
    let mut f1 = vec![[0.0; 3]; at.len()];
    let mut f2 = f1.clone();
    let e1 = field.energy_and_forces(&at, &mut f1);
    let e2 = field.energy_and_forces(&at, &mut f2);
    assert_eq!((e1.to_bits(), &f1), (e2.to_bits(), &f2));
}

/// The cross pairs' soft core follows the mixing rule: a decoupling whose field mixes geometrically
/// has other couplings, so the rule reaches the soft core and not only the evaluation's own pairs.
#[test]
fn the_soft_core_follows_the_mixing_rule() {
    let d = decoupled(false);
    let group: Vec<bool> = (0..d.at.len()).map(|i| i < 12).collect();
    let geometric =
        PeriodicDecoupling::new(&d.field.clone().with_mixing(Mixing::Geometric), &group)
            .expect("benzene");
    let l = Lambda::new(0.0, 0.0, 0.3);
    let (a, b) = (
        d.decoupling.coupling(&d.at, l),
        geometric.coupling(&d.at, l),
    );
    println!(
        "λ_v 0.3: ∂U/∂λ_v {:+.6e} J under LB, {:+.6e} geometric",
        a.gradient[2], b.gradient[2]
    );
    assert!((a.gradient[2] - b.gradient[2]).abs() > 1e-3 * a.gradient[2].abs());
}

// ---------------------------------------------------------------------------------------------
// Review: a face, ions after the rules, and HEAD's bits

/// **A scaled 1-4 pair across a face is taken at its minimum image.** Sixteen octanols with every
/// rule on, every molecule moved whole by 0.37 of the box so that some reach past a face; then a
/// copy with every coordinate past a face brought back by one edge — `p − L`, exact for `p` in
/// `[L, 2L]` by Sterbenz's lemma, so the two are the same points to the bit, one image apart. At
/// least one scaled 1-4 pair is then split by a face. The energy, forces and virial of the two
/// agree to the measured rounding the images leave, classical and by the mesh, held at `64 ε` of
/// the energy's parts and of the largest force and virial component; taken as the raw difference
/// instead of its image, a split 1-4 pair would be a box apart.
#[test]
fn a_one_four_pair_across_a_face_is_its_minimum_image() {
    let o = Octanol::new(None);
    let l = small_box(&o, 0x5DD1, 1e-5);
    let field = with_rules(l.field.clone(), &l.component, 0);
    let edge = l.cell().lengths();
    let whole: Vec<[f64; 3]> =
        l.at.iter()
            .map(|p| [0, 1, 2].map(|a| p[a] + 0.37 * edge[a]))
            .collect();
    let cut: Vec<[f64; 3]> = whole
        .iter()
        .map(|p| {
            [0, 1, 2].map(|a| {
                if p[a] >= edge[a] {
                    p[a] - edge[a]
                } else {
                    p[a]
                }
            })
        })
        .collect();
    for (p, q) in whole.iter().zip(&cut) {
        for a in 0..3 {
            assert!(p[a] == q[a] || p[a] - edge[a] == q[a]);
        }
    }
    let split = field
        .one_four_pairs()
        .iter()
        .filter(|&&[i, j]| (0..3).any(|a| (cut[i][a] - cut[j][a]).abs() > 0.5 * edge[a]))
        .count();
    println!(
        "{split} of {} scaled 1-4 pairs split by a face",
        field.one_four_pairs().len()
    );
    assert!(split > 0, "no 1-4 pair is split by a face");
    let m = mesh_for(&field);
    for (name, f) in [
        ("classical", field.clone()),
        ("mesh", field.clone().with_mesh(m)),
    ] {
        let (a, b) = (f.evaluate(&whole), f.evaluate(&cut));
        let scale = magnitude(&a.energy);
        let fmax = largest_force(&a.forces);
        let wmax = a
            .virial
            .iter()
            .flatten()
            .fold(0.0f64, |w, x| w.max(x.abs()));
        let de = (a.energy.total - b.energy.total).abs();
        let df = a
            .forces
            .iter()
            .zip(&b.forces)
            .flat_map(|(x, y)| (0..3).map(move |k| (x[k] - y[k]).abs()))
            .fold(0.0f64, f64::max);
        let mut dw = 0.0f64;
        for i in 0..3 {
            for j in 0..3 {
                dw = dw.max((a.virial[i][j] - b.virial[i][j]).abs());
            }
        }
        println!(
            "{name}: energy {:.2} ε of its parts, forces {:.2} ε of the largest, virial {:.2} ε of \
             its largest",
            de / (EPS * scale),
            df / (EPS * fmax),
            dw / (EPS * wmax)
        );
        assert!(de <= 64.0 * EPS * scale, "{name}: energy {de:e}");
        assert!(df <= 64.0 * EPS * fmax, "{name}: forces {df:e}");
        assert!(dw <= 64.0 * EPS * wmax, "{name}: virial {dw:e}");
    }
}

/// **Ions added after the rules are ions added before them**: a chloride after
/// `with_nonbonded`, `with_mixing` and `with_one_four`, and one before, are the same force field,
/// and give the same evaluation to the bit, with the cell list and with the neighbour list — the
/// scaled 1-4 pairs' table of what the pair search leaves out growing with the atoms.
#[test]
fn ions_after_the_rules_are_ions_before_them() {
    let o = Octanol::new(None);
    let l = small_box(&o, 0x10F, 1e-5);
    let after = with_rules(l.field.clone(), &l.component, 0).with_ions(&[Ion::chloride()]);
    let before = with_rules(
        l.field.clone().with_ions(&[Ion::chloride()]),
        &l.component,
        0,
    );
    assert_eq!(after, before);
    assert_eq!(after.one_four_pairs().len(), 66 * 16);
    let mut at = l.at.clone();
    at.push(l.cell().lengths().map(|x| 0.5 * x));
    assert_eq!(
        evaluation_bits(&after.evaluate(&at)),
        evaluation_bits(&before.evaluate(&at))
    );
    let skin = PeriodicForceField::NEIGHBOUR_SKIN;
    assert_eq!(
        evaluation_bits(&after.clone().with_neighbour_list(skin).evaluate(&at)),
        evaluation_bits(&before.evaluate(&at))
    );
}

/// The bits HEAD `b6558df` gave [`pinned::waters_and_an_ion`]: [`pinned::all_bits`].
const HEAD_WATER: [u64; 36] = [
    0x3c5e41e2b956f198,
    0xbc02c29144068c8f,
    0xbc556a09dfac89c6,
    0xbc4402d0f0084155,
    0x3c35364f640f0d8f,
    0xbc8d96509e75111b,
    0x3c8b807403d4b370,
    0xbbbd41aea4fc597c,
    0xbc556a09dfac89c6,
    0xbc02c29144068c8f,
    0xc1a1503af9e3a978,
    0x3c4083889f1466dc,
    0xbe213c3689b1730e,
    0xbe1aa73d2bd4149a,
    0xbe22c6df784c6d6a,
    0x3e116e61fc5c9e89,
    0x3e11a47ad390642c,
    0x3e061f9ae3b58969,
    0x3e12acb4ade6d93c,
    0x3e10da4fd337c49f,
    0x3e070d5eabaf1058,
    0x3e124b4f8d4c047b,
    0x3e020cd6eec93950,
    0xbde9b8e85679a70f,
    0xbe2c0b8820971a5a,
    0x3e5b3b9b5621f6bb,
    0x3e4d9ebf06b099fc,
    0x3c881b48dcbfc69c,
    0xbc3d944521c39450,
    0x3c49577c8088e5a7,
    0xbc3d944521c39454,
    0x3c8462950f0be2b6,
    0xbc4c4d391af73afd,
    0x3c49577c8088e5a7,
    0xbc4c4d391af73afd,
    0x3c60bbd577dce06d,
];

/// The bits HEAD `b6558df` gave [`pinned::ideal_octanols`]: [`pinned::nonbonded_bits`].
const HEAD_OCTANOL: [u64; 11] = [
    0xbc0d705c3fa8d07c,
    0xbc0886aa495ed52f,
    0xbc08d72364558300,
    0xbb5f79711767d548,
    0x3c0bbc9ea68e3f0c,
    0xbc50574c60e89950,
    0x3c4d659b9405254f,
    0xb54825cc5c3cd839,
    0xbc08d72364558300,
    0xbc0886aa495ed52f,
    0xc198429fb3410914,
];

/// [`pinned::digest`] of the 1 024 bits HEAD `b6558df` gave [`pinned::tails_over_cutoffs`]. The
/// tail rearranged as `p/9/r_c⁹` moves 6 of them, measured by that rearrangement of HEAD's own.
const HEAD_TAILS: u64 = 0x4d6c09f9ff9f11b9;

/// **The default path is HEAD's bits**, computed by HEAD `b6558df` before any of this existed, on
/// two systems whose every input is the same on every platform ([`pinned`]): 64 TIP3P waters and a
/// chloride — every energy term, both corrections, five atoms' forces and the virial — and nine
/// octanols at the dictionary's ideal coordinates with AM1-BCC's charges, their non-bonded terms;
/// and both systems' long-range corrections at 256 cutoffs each, since one cutoff's rounding can
/// hide a rearranged tail: two did, and sixteen caught it once in 64.
/// The plain field, the field with every rule set to its default explicitly, and the neighbour list.
/// A test that compares the new code with itself cannot see a refactor of what the two share; this
/// can.
#[test]
fn the_default_path_is_heads_bits() {
    let (water, at) = pinned::waters_and_an_ion();
    let n = at.len();
    let explicit = water
        .clone()
        .with_nonbonded(0..n, water.van_der_waals().to_vec())
        .expect("its own values")
        .with_mixing(Mixing::Geometric)
        .with_one_four(1.0, 1.0);
    let listed = water
        .clone()
        .with_neighbour_list(PeriodicForceField::NEIGHBOUR_SKIN);
    for (name, f) in [
        ("plain", &water),
        ("explicit", &explicit),
        ("listed", &listed),
    ] {
        assert_eq!(
            pinned::all_bits(f, &f.evaluate(&at)),
            HEAD_WATER,
            "water, {name}"
        );
    }
    let (octanols, at) = pinned::ideal_octanols();
    let listed = octanols
        .clone()
        .with_neighbour_list(PeriodicForceField::NEIGHBOUR_SKIN);
    for (name, f) in [("plain", &octanols), ("listed", &listed)] {
        assert_eq!(
            pinned::nonbonded_bits(f, &f.energy(&at)),
            HEAD_OCTANOL,
            "octanol, {name}"
        );
    }
    let tails = pinned::tails_over_cutoffs();
    assert_eq!(tails.len(), 1024);
    assert_eq!(pinned::digest(&tails), HEAD_TAILS);
}
