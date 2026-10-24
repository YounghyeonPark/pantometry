//! **Liquid n-hexadecane with GAFF's non-bonded terms on UFF's bonded ones**: step C-1, the oil of
//! benzene's oil:water partition. The molecule is the dictionary's `R16`; GAFF's `c3` and `hc` by
//! Lorentz–Berthelot, GAFF's 1-4 scales, and AM1-BCC charges built from FreeSolv's decane by
//! position (`tests/hexadecane/`). n-Decane (`D10`) is the fallback if hexadecane will not stay
//! liquid, measured by the same protocol.
//!
//! **Kept out of the published crate** (`Cargo.toml`'s `exclude`) with `tests/hexadecane/` and
//! `tests/gaff/`: GAFF's numbers are read from `gaff.dat`, GPL v3.
//!
//! Default, unoptimised, seconds each:
//!
//! - **the charge transfer reproduces octane, nonane and decane**: the construction on FreeSolv's
//!   own three files, every atom within [`HOMOLOGUE_BOUND`] of the file in octane and nonane, and
//!   decane exactly its own values plus the one shift;
//! - **hexadecane's charges are decane's and sum to zero**: every atom a decane value plus the shift,
//!   the positions summing to −0.0047 e before it, and to zero after;
//! - **the model is what it says**: GAFF's `c3` and `hc` read off the elements and equal to FreeSolv's
//!   files' types, on every atom of a box; Lorentz–Berthelot; GAFF's 1-4 scales on `9 (n − 1)` pairs a
//!   molecule, the closed form an n-alkane's 1-4 count is, and the breadth-first search's;
//! - **the box is at its density, with no overlap and every bond under half an edge**;
//! - **every held bond is at its length**, and **the degrees of freedom** are `3N − N_held − 3` in
//!   NVE and `3N − N_held` in a bath;
//! - **constrained NVE falls as `h²`** on four hexadecanes, with the excess over four shrinking;
//! - **what the liquid is judged by is right on closed cases**: the order parameter, the dihedral
//!   count and the mean-square displacement.
//!
//! Ignored, release: [`the_liquid_measured`] (and [`the_decane_liquid_measured`]), the melt, the
//! pressure at experiment's density, `ΔH_vap`, the carbons' `g(r)`, the self-diffusion coefficient
//! and the chains' order; [`the_cost_measured`]; [`the_h_squared_band_measured`].
//!
//! # What the hexadecane measured: not a liquid
//!
//! Run once on 2026-10-24, release, one core, 45 800 steps, 2 177 s (46.7 ms a step with the
//! records); nothing asserted. 98 molecules, 46.0 × 36.3 × 28.7 Å at 0.7701 g/cm³.
//!
//! | | measured | experiment / reference |
//! | --- | --- | --- |
//! | P at 0.7701 g/cm³ (bar) | −1179 ± 28; halves −1177, −1181 | 1 |
//! | ΔH_vap (kJ/mol) | 118.5 ± 2.6 | 81.35 (CRC via PubChem; NIST 81.3 ± 0.6) |
//! | U/N in the sample's halves (kcal/mol) | 15.563, 13.681 — still falling | |
//! | nematic order of the end-to-end vectors | 0.859; halves 0.884, 0.835 | isotropic set of 98: 0.081 ± 0.029; lattice 1 |
//! | trans fraction of C–C–C–C | 0.961; halves 0.948, 0.974 — rising | the start 1 |
//! | end-to-end (Å) | 18.27 | all-trans 18.96 |
//! | g(r) first peak, CH₃–CH₃ / interior / all C | 4.67 Å, 2.63 / 5.30 Å, 1.85 / 5.29 Å, 1.68 | |
//! | MSD of the centres at 5, 10, 20 ps | 2.08, 3.31, 5.77 Å²; log-log slope 0.67–0.80 | 1 for diffusion |
//! | D from 5–20 ps (m²/s) | 4.0 × 10⁻¹⁰ ± 0.35 (seven groups); Yeh–Hummer +0.49 | ≈ 3.9 × 10⁻¹⁰ (quoted, not read) |
//!
//! - **The verdict: an aligned, ordering phase, not a liquid.** The chains never lost the lattice's
//!   alignment — the order fell from 1 to 0.94 in 30 ps at 600 K and to 0.84 by the end, ten
//!   isotropic standard deviations away — and at 298.15 K they straighten, trans 0.85 → 0.98, as `U`
//!   falls 5 kcal/mol over 60 ps with no plateau. The MSD is not diffusive (slope 0.67–0.80), so the
//!   `D` that agrees with experiment's is a slope through sub-diffusion and proves nothing.
//! - **ΔH_vap is 37 kJ/mol over the liquid's** and between the liquid's 81 and a crystal's
//!   sublimation (the liquid's plus a fusion enthalpy of about 50, quoted): the box is partly
//!   crystalline in energy too. The tension at experiment's density says the model's own liquid,
//!   if it has one, is denser.
//! - **Two causes are not separated**: the melt (a parallel lattice of 20 Å rods at liquid density,
//!   which 30 ps at 600 K does not turn) and the model (UFF's torsions with GAFF's van der Waals,
//!   whose own melting point is not known and may be above 298 K: 6.8 K of margin in experiment).
//!   Either way this box is not the oil the partition needs.
//!
//! # And the decane: liquid in pressure and mobility, still aligned
//!
//! The same protocol, run once on 2026-10-24, 1 533 s (32.9 ms a step with the records); 112
//! molecules, 30.8 × 35.9 × 33.0 Å at 0.7255 g/cm³:
//!
//! | | measured | experiment / reference |
//! | --- | --- | --- |
//! | P at 0.7255 g/cm³ (bar) | −52 ± 31; halves −77, −26 | 1 |
//! | ΔH_vap (kJ/mol) | 58.3 ± 0.7 | 51.42 (CRC via PubChem) |
//! | U/N in the sample's halves (kcal/mol) | 12.937, 12.595 | |
//! | nematic order | 0.568; halves 0.586, 0.550 — falling | isotropic set of 112: 0.076 ± 0.027 |
//! | trans fraction | 0.825; halves 0.815, 0.836 | |
//! | g(r) first peak, CH₃–CH₃ / interior / all C | 4.15 Å, 1.67 / 5.38 Å, 1.31 / 5.27 Å, 1.30 | |
//! | MSD at 5, 10, 20 ps | 3.01, 4.96, 8.61 Å²; log-log slope 0.72–0.80 | |
//! | D from 5–20 ps (m²/s) | 6.1 × 10⁻¹⁰ ± 0.8; Yeh–Hummer +2.2 | not read |
//!
//! - **Liquid by its pressure, its g(r) and its cohesion**: experiment's density is the model's to
//!   within a few tens of bar, ΔH_vap is 6.9 kJ/mol (13%) high — octanol's GAFF hybrid was 6.3 —
//!   and the carbons' first peaks are 1.3, a liquid's.
//! - **But not yet isotropic, so not equilibrated**: the order is 0.57, eighteen isotropic standard
//!   deviations off and still falling, and the MSD has not reached slope one by 20 ps. The lattice's
//!   alignment is what is left, and 30 ps at 600 K took it only to 0.67.

// The sums over atoms and axes are written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

mod freesolv;
mod gaff;
mod hexadecane;
mod nonbonded;
mod octanol;

use gaff::{gaff, ONE_FOUR};
use hexadecane::*;
use octanol::{
    bathed, bonds_to_hydrogen, closest_contact, gas, kcal, pressure, production, Liquid, Stage,
    ANGSTROM, FS, KELVIN, KJ_PER_MOL, OVERLAP,
};
use pantometry_forcefield::uff::KCAL_PER_MOL;
use pantometry_forcefield::{
    Bath, Component, Element, Estimate, EwaldParameters, Mixing, MolecularDynamics, PeriodicBox,
    PeriodicForceField, PmeParameters, Shake,
};
use pantometry_units::BOLTZMANN;
use std::f64::consts::PI;

const EPS: f64 = f64::EPSILON;

/// The small boxes: one molecule along its axis, three by three across, 450 atoms.
const SMALL: [usize; 3] = [1, 3, 3];
/// Their cutoff: under half their shortest edge, 6.5 Å.
const SMALL_CUTOFF: f64 = 6.0 * ANGSTROM;

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn len(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// A Mol2's elements and bonds.
fn graph(m: &freesolv::Mol2) -> (Vec<Element>, Vec<[usize; 2]>) {
    (m.atoms.iter().map(|a| a.element).collect(), m.bonds.clone())
}

/// The small box of hexadecane at its density, GAFF's.
fn small_box(a: &Alkane, seed: u64) -> (Liquid, PeriodicForceField) {
    let l = Liquid::lattice(
        &a.molecule,
        SMALL,
        HEXADECANE.density,
        seed,
        SMALL_CUTOFF,
        1e-5,
    );
    let field = with_gaff(l.field.clone(), a);
    (l, field)
}

// ---------------------------------------------------------------------------------------------
// The charges

/// **The construction reproduces octane, nonane and decane**: [`transferred`] from decane's
/// positions applied to each of FreeSolv's three files' own graphs, against each file's charges.
/// Octane and nonane within [`HOMOLOGUE_BOUND`] on every atom — measured 0.001035 and 0.001193 e;
/// decane is its own values with the one shift on every atom, to the rounding of one addition, which
/// a construction that took three carbons from each end, or none of the shift, is not (0.0002 and
/// 0.0000156 e off).
#[test]
fn the_construction_reproduces_octane_nonane_and_decane() {
    let decane = decane_positions();
    println!("| file | n | Σq as written | Σq built, before the shift | shift | worst atom | where |\n| --- | --- | --- | --- | --- | --- | --- |");
    for (name, mol2) in [
        ("octane", freesolv::octane()),
        ("nonane", freesolv::nonane()),
        ("decane", freesolv::decane()),
    ] {
        let (elements, bonds) = graph(&mol2);
        let t = transferred(&decane, &elements, &bonds);
        let (mut worst, mut at) = (0.0f64, 0);
        for (i, a) in mol2.atoms.iter().enumerate() {
            let d = (t.charges[i] - a.charge).abs();
            if d > worst {
                worst = d;
                at = i;
            }
        }
        let written: f64 = mol2.atoms.iter().map(|a| a.charge).sum();
        println!(
            "| {name} | {} | {written:+.4} | {:+.4} | {:+.3e} | {worst:.6} | {} |",
            mol2.atoms.len(),
            t.net_before,
            t.shift,
            mol2.atoms[at].name
        );
        if name == "decane" {
            for (i, a) in mol2.atoms.iter().enumerate() {
                let off = t.charges[i] - a.charge - t.shift;
                assert!(
                    off.abs() <= 4.0 * EPS * (a.charge.abs() + t.shift.abs()),
                    "decane's {} is not its own value and the shift: {off:e}",
                    a.name
                );
            }
            assert!(t.shift != 0.0, "decane's shift is there");
        } else {
            assert!(
                worst <= HOMOLOGUE_BOUND,
                "{name}: {} is {worst} e off the file",
                mol2.atoms[at].name
            );
        }
    }
}

/// **Hexadecane's charges are decane's, and sum to zero**: every carbon and hydrogen of `R16` one of
/// decane's per-position values plus the shift, the outer four at each end their own end's and the
/// eight between at the interior's; before the shift they sum to −0.0047 e (decane's −0.0005, six
/// more interior CH₂ at −0.0007 each), the shift is +0.000094 e an atom, and after it the sum is zero
/// to the rounding of fifty terms.
#[test]
fn hexadecanes_charges_are_decanes_and_sum_to_zero() {
    let c = Component::from_ccd(R16).expect("R16");
    let elements: Vec<Element> = c.atoms().iter().map(|a| a.element).collect();
    let bonds: Vec<[usize; 2]> = c.bonds().iter().map(|b| b.atoms).collect();
    let d = decane_positions();
    let t = transferred(&d, &elements, &bonds);
    let (order, hydrogens) = chain(&elements, &bonds);
    assert_eq!(order.len(), 16);
    // Written out position by position, independently of `transferred`'s branches.
    let source = [0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 6, 7, 8, 9];
    println!("| position | atom | decane's | carbon | hydrogen |\n| --- | --- | --- | --- | --- |");
    for (k, &c_i) in order.iter().enumerate() {
        let (qc, qh) = d[source[k]];
        println!(
            "| {} | {} | C{} | {:+.6} | {:+.6} |",
            k + 1,
            c.atoms()[c_i].name,
            source[k] + 1,
            t.charges[c_i],
            t.charges[hydrogens[k][0]]
        );
        assert_eq!(
            t.charges[c_i].to_bits(),
            (qc + t.shift).to_bits(),
            "carbon {}",
            k + 1
        );
        for &h in &hydrogens[k] {
            assert_eq!(
                t.charges[h].to_bits(),
                (qh + t.shift).to_bits(),
                "hydrogen on {}",
                k + 1
            );
        }
    }
    let scale: f64 = t.charges.iter().map(|x| x.abs()).sum();
    let net: f64 = t.charges.iter().sum();
    println!(
        "before the shift {:+.6} e, the shift {:+.3e} e, after {net:+.2e} e (Σ|q| {scale:.4})",
        t.net_before, t.shift
    );
    assert!((t.net_before + 0.0047).abs() <= 1e-12, "{}", t.net_before);
    assert!((t.shift - 0.000094).abs() <= 1e-12, "{}", t.shift);
    assert!(
        net.abs() <= 64.0 * EPS * scale,
        "the molecule carries {net} e"
    );
}

// ---------------------------------------------------------------------------------------------
// The model

/// **The model is what it says**: GAFF's type of every atom of the small box, read off its element
/// here — `c3` for carbon, `hc` for hydrogen — is FreeSolv's: each of the three files' type column
/// is the same rule. Every atom of the box carries that type's `[x, D]` to the bit and its molecule's
/// transferred charge to the bit; the mixing is Lorentz–Berthelot and the 1-4 scales GAFF's, on
/// `9 (n − 1)` = 135 pairs a molecule — an n-alkane's 1-4 pairs are the dihedrals about its `n − 1`
/// C–C bonds, three by three — which the breadth-first search counts too.
#[test]
fn the_model_is_what_it_says() {
    for mol2 in [freesolv::octane(), freesolv::nonane(), freesolv::decane()] {
        for a in &mol2.atoms {
            let want = if a.element == Element::C { "c3" } else { "hc" };
            assert_eq!(a.kind, want, "{}: {}", mol2.name, a.name);
        }
    }
    let a = Alkane::new(R16);
    let (l, field) = small_box(&a, 0xC16);
    let n = a.molecule.at.len();
    assert_eq!(n, 50);
    assert_eq!(a.molecule.component.formula(), Some("C16 H34"));
    let vdw = field.van_der_waals();
    let q = field.charges();
    for (i, atom) in l.component.atoms().iter().enumerate() {
        let kind = match atom.element {
            Element::C => "c3",
            Element::H => "hc",
            other => panic!("{other:?}"),
        };
        assert_eq!(
            vdw[i].map(f64::to_bits),
            gaff(kind).map(f64::to_bits),
            "atom {i}"
        );
        assert_eq!(
            q[i].to_bits(),
            a.transfer.charges[i % n].to_bits(),
            "atom {i}"
        );
    }
    assert_eq!(field.mixing(), Mixing::LorentzBerthelot);
    assert_eq!(field.one_four(), ONE_FOUR);
    let per = nonbonded::one_four_by_search(n, &nonbonded::bonds(&a.molecule.component)).len();
    assert_eq!(per, 9 * 15);
    assert_eq!(field.one_four_pairs().len(), per * l.octanols);
    println!(
        "{} hexadecanes, {} atoms, {} 1-4 pairs scaled; c3 {:?}, hc {:?}",
        l.octanols,
        q.len(),
        field.one_four_pairs().len(),
        gaff("c3"),
        gaff("hc")
    );
}

// ---------------------------------------------------------------------------------------------
// The box

/// **The box is at its density to rounding, nothing overlaps, and every bond is under half an
/// edge**: the small box and both measurement boxes, each the molecules' mass over the volume within
/// `8 ε` of the density; the small box's own masses within the rounding of their sum; its closest
/// contact by brute force over the 27 images at least [`OVERLAP`] and the lattice's gap to within the
/// conformer's asymmetry; and every bond of every box shorter than half the box's shortest edge, so
/// that the field's minimum image makes each molecule whole.
#[test]
fn the_box_is_at_its_density_with_no_overlaps() {
    let hexadecane = Alkane::new(R16);
    let decane = Alkane::new(D10);
    for (a, oil, counts) in [
        (&hexadecane, &HEXADECANE, SMALL),
        (&hexadecane, &HEXADECANE, HEXADECANE.lattice),
        (&decane, &DECANE, DECANE.lattice),
    ] {
        let (cell, site, gap) = Liquid::cell_of(&a.molecule, counts, oil.density);
        let m = (counts[0] * counts[1] * counts[2]) as f64;
        let density = m * a.molecule.mass() / cell.volume();
        let longest = a
            .molecule
            .component
            .bonds()
            .iter()
            .map(|b| len(sub(a.molecule.at[b.atoms[0]], a.molecule.at[b.atoms[1]])))
            .fold(0.0f64, f64::max);
        println!(
            "{} {counts:?}: {} molecules, box {:.3} × {:.3} × {:.3} Å, site {:.3} × {:.3} × {:.3} Å, \
             gap {:.3} Å, density {:.10} g/cm³; longest bond {:.3} Å",
            oil.name,
            m,
            cell.lengths()[0] / ANGSTROM,
            cell.lengths()[1] / ANGSTROM,
            cell.lengths()[2] / ANGSTROM,
            site[0] / ANGSTROM,
            site[1] / ANGSTROM,
            site[2] / ANGSTROM,
            gap / ANGSTROM,
            density / 1e3,
            longest / ANGSTROM
        );
        assert!(
            (density / oil.density - 1.0).abs() <= 8.0 * EPS,
            "{counts:?}"
        );
        assert!(longest < cell.half_shortest_edge());
        assert!(gap >= OVERLAP, "a gap of {gap:e}");
    }
    let (l, _) = small_box(&hexadecane, 2);
    let bound = (l.masses.len() as f64 + 8.0) * EPS;
    println!(
        "the small box's own masses over its volume: {:+.2e} of the density off it, bound {bound:.1e}",
        l.density() / HEXADECANE.density - 1.0
    );
    assert!((l.density() / HEXADECANE.density - 1.0).abs() <= bound);
    let (r, i, j) = closest_contact(&l.at, &l.molecules, &l.cell());
    println!(
        "closest contact {:.4} Å, atoms {i} and {j} ({} and {}); gap {:.4} Å",
        r / ANGSTROM,
        l.component.atoms()[i].name,
        l.component.atoms()[j].name,
        l.gap / ANGSTROM
    );
    assert!(r >= OVERLAP, "an overlap: {:.3} Å", r / ANGSTROM);
    assert!(
        r >= 0.9 * l.gap && r <= 1.1 * l.gap,
        "{r:e} against {:e}",
        l.gap
    );
}

// ---------------------------------------------------------------------------------------------
// The constraints

/// **Every bond to a hydrogen is held, and at its length**: the X–H bonds counted from the elements
/// are the constraints' exactly, 34 a molecule, each at its UFF stretch's natural length; and after
/// 40 steps of constrained BAOAB at 2 fs from the lattice, under GAFF's terms, each is within twice
/// SHAKE's default tolerance of it, with no velocity along itself beyond twice RATTLE's.
#[test]
fn every_held_bond_is_at_its_length() {
    let a = Alkane::new(R16);
    let (l, field) = small_box(&a, 3);
    let mut from_elements = bonds_to_hydrogen(&l.component);
    assert_eq!(from_elements.len(), 34 * l.octanols);
    let shake = l.shake();
    let mut held: Vec<[usize; 2]> = shake.bonds().to_vec();
    let norm = |v: &mut Vec<[usize; 2]>| {
        for b in v.iter_mut() {
            *b = [b[0].min(b[1]), b[0].max(b[1])];
        }
        v.sort_unstable();
    };
    norm(&mut from_elements);
    norm(&mut held);
    assert_eq!(held, from_elements, "every X–H is held, and nothing else");
    let natural: std::collections::HashMap<[usize; 2], f64> = field
        .bonded()
        .expect("bonded")
        .stretches()
        .iter()
        .map(|s| {
            (
                [s.atoms[0].min(s.atoms[1]), s.atoms[0].max(s.atoms[1])],
                s.natural_length,
            )
        })
        .collect();
    for (b, d) in shake.bonds().iter().zip(shake.lengths()) {
        let key = [b[0].min(b[1]), b[0].max(b[1])];
        assert_eq!(natural[&key], *d, "{b:?} is held at its natural length");
    }
    let (md, at) = bathed(&l.dynamics(), &field, &l.at, 2.0 * FS, 40, 0x16C1);
    let v = md.velocities();
    let (mut worst, mut worst_v) = (0.0f64, 0.0f64);
    for &[i, j] in &from_elements {
        let d = natural[&[i, j]];
        let r = sub(at[i], at[j]);
        let u = sub(v[i], v[j]);
        worst = worst.max((len(r) - d).abs() / (Shake::DEFAULT_TOLERANCE * d));
        worst_v = worst_v.max(dot(r, u).abs() / (Shake::DEFAULT_TOLERANCE * len(r) * len(u)));
    }
    println!(
        "after 40 steps: the worst held bond {worst:.3} of SHAKE's tolerance off its length, the \
         worst relative velocity along a bond {worst_v:.3} of it"
    );
    assert!(worst <= 2.0, "a held bond is off its length: {worst}");
    assert!(
        worst_v <= 2.0,
        "a held bond has velocity along itself: {worst_v}"
    );
}

/// **The degrees of freedom are `3N − N_held − 3`** in NVE with the net momentum out, and
/// `3N − N_held` in a bath, `N_held` counted from the elements: 34 a molecule.
#[test]
fn the_degrees_of_freedom_are_three_n_less_the_held_less_three() {
    let a = Alkane::new(R16);
    let (l, _) = small_box(&a, 4);
    let n = l.at.len();
    let held = bonds_to_hydrogen(&l.component).len();
    assert_eq!(held, 34 * 9);
    let nve = l
        .dynamics()
        .thermalised(&l.at, KELVIN, 4)
        .without_net_momentum();
    assert_eq!(nve.degrees_of_freedom(), 3 * n - held - 3);
    let bathed = l.dynamics().with_bath(Bath::Langevin {
        temperature: KELVIN,
        friction: 1e12,
        seed: 4,
    });
    assert_eq!(bathed.degrees_of_freedom(), 3 * n - held);
    println!(
        "{n} atoms, {held} held: {} in NVE, {} in a bath",
        nve.degrees_of_freedom(),
        bathed.degrees_of_freedom()
    );
}

// ---------------------------------------------------------------------------------------------
// Constrained NVE

/// Four hexadecanes, a 1 × 2 × 2 lattice at 0.6 g/cm³ so that no contact is a clash a 2 fs step
/// cannot start from, turned by `seed`, moved into a box 56 Å across with `r_c` = 27 Å: the cluster
/// spans under 26 Å, so no pair reaches the cutoff and no image comes inside it, the potential is
/// smooth and NVE's energy error is the integrator's. GAFF's terms; every X–H held.
fn cluster(seed: u64) -> (PeriodicForceField, Vec<[f64; 3]>, MolecularDynamics) {
    let a = Alkane::new(R16);
    let l = Liquid::lattice(&a.molecule, [1, 2, 2], 600.0, seed, 4.0 * ANGSTROM, 1e-5);
    let side = 56.0 * ANGSTROM;
    let big = PeriodicBox::cubic(side);
    let p = EwaldParameters::for_accuracy(&big, 27.0 * ANGSTROM, 1e-5);
    let field = PeriodicForceField::new(&l.component, &l.types, big, p)
        .expect("the cluster")
        .with_charges(l.charges.clone());
    let field = with_gaff(field, &a);
    let small = l.cell().lengths();
    let shift = [0, 1, 2].map(|k| 0.5 * (side - small[k]));
    let at: Vec<[f64; 3]> =
        l.at.iter()
            .map(|p| [p[0] + shift[0], p[1] + shift[1], p[2] + shift[2]])
            .collect();
    let span = at
        .iter()
        .flat_map(|p| at.iter().map(move |q| len(sub(*p, *q))))
        .fold(0.0f64, f64::max);
    assert!(span < 26.0 * ANGSTROM, "span {span:e}");
    let shake = Shake::to_hydrogen(&l.component, field.bonded().expect("bonded"));
    let md = MolecularDynamics::new(l.masses.clone()).with_bond_constraints(shake);
    (field, at, md)
}

/// The RMS departure of `kinetic + potential` from its start over `steps` of `dt`, NVE, from
/// 298.15 K velocities with the net momentum out.
fn nve_rms(
    field: &PeriodicForceField,
    start: &[[f64; 3]],
    md: &MolecularDynamics,
    dt: f64,
    steps: usize,
    seed: u64,
) -> f64 {
    let mut at = start.to_vec();
    let mut md = md
        .clone()
        .thermalised(&at, KELVIN, seed)
        .without_net_momentum();
    md.prepare(field, &at);
    let e0 = md.kinetic_energy() + md.potential_energy().expect("prepared");
    let mut squares = 0.0;
    for _ in 0..steps {
        md.step(field, &mut at, dt);
        let e = md.kinetic_energy() + md.potential_energy().expect("stepped");
        squares += (e - e0) * (e - e0);
    }
    (squares / steps as f64).sqrt()
}

/// The three ratios of the RMS energy error per halving, 2 → 1, 1 → 0.5 and 0.5 → 0.25 fs, each over
/// 0.1 ps from the same start and velocities, and the RMS error at 2 fs.
fn ratios(seed: u64) -> ([f64; 3], f64) {
    let (field, at, md) = cluster(seed);
    let e = [(2.0, 50), (1.0, 100), (0.5, 200), (0.25, 400)]
        .map(|(h, n)| nve_rms(&field, &at, &md, h * FS, n, seed));
    ([e[0] / e[1], e[1] / e[2], e[2] / e[3]], e[0])
}

/// The bands, **measured** over sixteen starts by [`the_h_squared_band_measured`] (release, run
/// 2026-10-24; the debug run's start is the same bits): the ratio from 2 fs 4.2293–4.2504, from
/// 1 fs 4.0550–4.0603, from 0.5 fs 4.0144–4.0159, and the excess over four shrinking by
/// 0.2581–0.2646 per halving — octanol's to within its own bands. Each ratio's band is its measured
/// middle with three times its half-spread either side, and the shrink's 3.5 times, as octanol's
/// and W4's.
const H_SQUARED_FROM_TWO: std::ops::Range<f64> = 4.208..4.272;
const H_SQUARED_FROM_ONE: std::ops::Range<f64> = 4.0498..4.0655;
const H_SQUARED_FROM_HALF: std::ops::Range<f64> = 4.0129..4.0175;
const H_SQUARED_SHRINKS: std::ops::Range<f64> = 0.2499..0.2729;

/// **Constrained velocity Verlet on four hexadecanes**: every X–H held, NVE over 0.1 ps from 298 K
/// velocities. The RMS energy error falls by four per halving from 2 fs, and its excess over four
/// shrinks by four per halving, each held to its band.
#[test]
fn constrained_nve_on_hexadecane_falls_as_h_squared() {
    let (r, e2) = ratios(0xC16);
    let shrink = (r[2] - 4.0) / (r[1] - 4.0);
    println!(
        "RMS ΔE over 0.1 ps at 2 fs {:.3e} kcal/mol; ratios {:.4}, {:.4}, {:.4}; the excess \
         shrinks by {shrink:.4}",
        kcal(e2),
        r[0],
        r[1],
        r[2]
    );
    assert!(H_SQUARED_FROM_TWO.contains(&r[0]), "ratio {}", r[0]);
    assert!(H_SQUARED_FROM_ONE.contains(&r[1]), "ratio {}", r[1]);
    assert!(H_SQUARED_FROM_HALF.contains(&r[2]), "ratio {}", r[2]);
    assert!(
        H_SQUARED_SHRINKS.contains(&shrink),
        "the excess shrinks by {shrink}"
    );
}

/// The measurement the bands come from: sixteen starts, each a lattice turned by its own seed and
/// velocities from it.
#[test]
#[ignore = "sixteen starts of four NVE runs: run with --release -- --ignored --nocapture"]
fn the_h_squared_band_measured() {
    let mut lo = [f64::INFINITY; 4];
    let mut hi = [f64::NEG_INFINITY; 4];
    for start in 0..16u64 {
        let (r, e2) = ratios(0xC16 + start);
        let shrink = (r[2] - 4.0) / (r[1] - 4.0);
        println!(
            "start {start}: ratios {:.4} {:.4} {:.4}, (r₃ − 4)/(r₂ − 4) {shrink:.4}, RMS at 2 fs \
             {:.3e} kcal/mol",
            r[0],
            r[1],
            r[2],
            kcal(e2)
        );
        for (k, x) in [r[0], r[1], r[2], shrink].into_iter().enumerate() {
            lo[k] = lo[k].min(x);
            hi[k] = hi[k].max(x);
        }
    }
    println!(
        "from 2 fs {:.4}–{:.4}, from 1 fs {:.4}–{:.4}, from 0.5 fs {:.4}–{:.4}, excess shrink \
         {:.4}–{:.4}",
        lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], lo[3], hi[3]
    );
    for k in 0..4 {
        let (mid, half) = (0.5 * (lo[k] + hi[k]), 0.5 * (hi[k] - lo[k]));
        let f = if k == 3 { 3.5 } else { 3.0 };
        println!("band {k}: {:.4}..{:.4}", mid - f * half, mid + f * half);
    }
}

// ---------------------------------------------------------------------------------------------
// What the liquid is judged by, on closed cases

/// **The order parameter, the dihedrals and the mean-square displacement on cases with closed
/// answers**: a symmetric matrix with eigenvalues 3, 1 and −2, turned, gives 3; parallel vectors
/// give order one and three orthogonal axes zero; the relaxed all-trans conformer is all trans; and
/// one molecule of `n` moved `L` a frame along `x` while the rest stay has, relative to the box's
/// centre, the mean-square displacement `L² (n − 1)/n²` at a lag of one frame.
#[test]
fn what_the_liquid_is_judged_by_on_closed_cases() {
    // R = rotation by 0.3 about z then 0.7 about x; Q = R diag(3, 1, −2) Rᵀ.
    let (c1, s1, c2, s2) = (0.3f64.cos(), 0.3f64.sin(), 0.7f64.cos(), 0.7f64.sin());
    let rz = [[c1, -s1, 0.0], [s1, c1, 0.0], [0.0, 0.0, 1.0]];
    let rx = [[1.0, 0.0, 0.0], [0.0, c2, -s2], [0.0, s2, c2]];
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i][j] += rx[i][k] * rz[k][j];
            }
        }
    }
    let d = [3.0, 1.0, -2.0];
    let mut q = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                q[i][j] += r[i][k] * d[k] * r[j][k];
            }
        }
    }
    let top = largest_eigenvalue(q);
    println!("largest eigenvalue {top:.15} against 3");
    assert!((top - 3.0).abs() <= 64.0 * EPS * 3.0);
    let s = nematic_order(&[[0.0, 0.6, 0.8]; 7]);
    assert!((s - 1.0).abs() <= 64.0 * EPS, "{s}");
    let s = nematic_order(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    assert!(s.abs() <= 64.0 * EPS, "{s}");
    let (m, sd) = isotropic_order(98, 400, 0x150);
    println!("98 isotropic directions: S = {m:.4} ± {sd:.4}");

    let a = Alkane::new(R16);
    let (trans, ends) = conformations(&a.molecule.at, std::slice::from_ref(&(0..50)), &a.carbons);
    println!(
        "the relaxed conformer: {trans} of its dihedrals trans, ends {:.3} Å apart",
        ends / ANGSTROM
    );
    assert_eq!(trans, 1.0);

    let n = 5;
    let l = 1e-10;
    let frames: Vec<Vec<[f64; 3]>> = (0..4)
        .map(|t| {
            (0..n)
                .map(|i| [if i == 0 { t as f64 * l } else { 0.0 }, i as f64 * l, 0.0])
                .collect()
        })
        .collect();
    let got = msd(&frames, &[1, 2], None);
    for (k, lag) in [1.0, 2.0].iter().enumerate() {
        let want = lag * lag * l * l * (n - 1) as f64 / (n * n) as f64;
        assert!(
            (got[k] - want).abs() <= 1e-12 * want,
            "{:e} against {want:e}",
            got[k]
        );
    }
}

// ---------------------------------------------------------------------------------------------
// The liquid, measured

/// **The melt and the measurement**: 0.4 ps at 0.5 fs in a 50 ps⁻¹ bath at 298.15 K, which takes
/// the lattice's strain out; 30 ps at 600 K and 2 fs to lose the lattice's order — the chains start
/// parallel, and turning a 20 Å molecule takes longer than octanol's 10 ps hot stage gave it; 20 ps at
/// 298.15 K to settle; then 40 ps sampled in a 1 ps⁻¹ bath, long enough for a self-diffusion
/// coefficient from 5 to 20 ps of lag. 45 800 steps, about 1.7 h at the 128 ms a step
/// [`the_cost_measured`] read.
const PLAN: [Stage; 4] = [
    Stage {
        name: "strain out",
        picoseconds: 0.4,
        dt: 0.5 * FS,
        friction: 50e12,
        temperature: KELVIN,
    },
    Stage {
        name: "hot",
        picoseconds: 30.0,
        dt: 2.0 * FS,
        friction: 5e12,
        temperature: 600.0,
    },
    Stage {
        name: "settle",
        picoseconds: 20.0,
        dt: 2.0 * FS,
        friction: 5e12,
        temperature: KELVIN,
    },
    Stage {
        name: "sample",
        picoseconds: 40.0,
        dt: 2.0 * FS,
        friction: 1e12,
        temperature: KELVIN,
    },
];

/// One record of the melt: time (ps), `U/N` (kcal/mol), `P` (bar), `T` (K), the nematic order and
/// the trans fraction.
type Row = (f64, f64, f64, f64, f64, f64);

/// **The liquid, measured**, for `oil`: its lattice at its density, `r_c` = 9 Å, PME at δ = 10⁻⁵
/// with the neighbour list, every X–H held, GAFF's terms; melted by [`PLAN`]. Reported, nothing
/// about the model asserted: the melt in 2.5 ps blocks (`U/N`, `P`, `T`, the order and the trans
/// fraction); over the sample `U/N` and `P` with their errors and halves; `ΔH_vap = U_gas − U/N + RT`
/// with the gas sampled in [`gas_phase`]; `g(r)` of the terminal and the interior carbons; the
/// centres' mean-square displacement and the self-diffusion coefficient; the order parameter against
/// an isotropic set's and the trans fraction against the gas's all-trans start.
fn measure(oil: &Oil) {
    let started = std::time::Instant::now();
    let a = Alkane::new(oil.entry);
    let l = Liquid::lattice(
        &a.molecule,
        oil.lattice,
        oil.density,
        0xC1,
        9.0 * ANGSTROM,
        1e-5,
    );
    let field = production(with_gaff(l.field.clone(), &a), &l);
    let n = l.octanols as f64;
    println!(
        "\n{}: {} molecules, {} atoms, box {:.3} × {:.3} × {:.3} Å at {:.6} g/cm³",
        oil.name,
        l.octanols,
        l.at.len(),
        l.cell().lengths()[0] / ANGSTROM,
        l.cell().lengths()[1] / ANGSTROM,
        l.cell().lengths()[2] / ANGSTROM,
        l.density() / 1e3
    );
    let ends = a.ends();
    let mut owner = vec![0usize; l.at.len()];
    for (m, r) in l.molecules.iter().enumerate() {
        for i in r.clone() {
            owner[i] = m;
        }
    }
    let terminal: Vec<usize> = l
        .molecules
        .iter()
        .flat_map(|m| ends.map(|e| m.start + e))
        .collect();
    let interior: Vec<usize> = l
        .molecules
        .iter()
        .flat_map(|m| {
            a.carbons[1..a.carbons.len() - 1]
                .iter()
                .map(move |&c| m.start + c)
        })
        .collect();
    let carbons: Vec<usize> = l
        .molecules
        .iter()
        .flat_map(|m| a.carbons.iter().map(move |&c| m.start + c))
        .collect();
    let order_of = |at: &[[f64; 3]]| nematic_order(&directions(at, &l.molecules, ends));
    let (iso, iso_sd) = isotropic_order(l.octanols, 1000, 0x150);
    let mut at = l.at.clone();
    let (trans0, ends0) = conformations(&at, &l.molecules, &a.carbons);
    println!(
        "as built: order {:.4}, trans {trans0:.4}, ends {:.2} Å apart; an isotropic set of {} gives {iso:.4} ± {iso_sd:.4}",
        order_of(&at),
        ends0 / ANGSTROM,
        l.octanols
    );
    let mut md = l.dynamics();
    let mut series: Vec<Row> = Vec::new();
    let (mut sample_u, mut sample_p, mut sample_s, mut sample_trans, mut sample_ends) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut centres_of: Vec<Vec<[f64; 3]>> = Vec::new();
    let width = 0.05 * ANGSTROM;
    let bins = (l.cell().half_shortest_edge() / width) as usize;
    let mut g = [vec![0.0; bins], vec![0.0; bins], vec![0.0; bins]];
    let mut frames = 0usize;
    let mut time = 0.0;
    let mut sample_ms = 0.0;
    for (k, stage) in PLAN.iter().enumerate() {
        md = md.with_bath(Bath::Langevin {
            temperature: stage.temperature,
            friction: stage.friction,
            seed: 0xC10 + k as u64,
        });
        if k == 0 {
            md = md.thermalised(&at, stage.temperature, 0xC10);
        }
        md.forget_forces();
        let steps = (stage.picoseconds * 1e-12 / stage.dt).round() as usize;
        let every = ((0.1e-12 / stage.dt).round() as usize).max(1);
        let t0 = std::time::Instant::now();
        for s in 0..steps {
            md.step(&field, &mut at, stage.dt);
            if (s + 1) % every == 0 {
                time += every as f64 * stage.dt;
                let u = kcal(md.potential_energy().expect("stepped")) / n;
                let p = pressure(&field, &l, &at) / 1e5;
                let order = order_of(&at);
                let (trans, e2e) = conformations(&at, &l.molecules, &a.carbons);
                series.push((time * 1e12, u, p, md.half_step_temperature(), order, trans));
                if stage.name == "sample" {
                    sample_u.push(u);
                    sample_p.push(p);
                    sample_s.push(order);
                    sample_trans.push(trans);
                    sample_ends.push(e2e);
                    centres_of.push(centres(&at, &l.molecules, &l.masses));
                    if (s + 1) % (5 * every) == 0 {
                        let cell = l.cell();
                        let parts = [
                            g_of_r(&cell, &at, &owner, &terminal, &terminal, true, width, bins),
                            g_of_r(&cell, &at, &owner, &interior, &interior, true, width, bins),
                            g_of_r(&cell, &at, &owner, &carbons, &carbons, true, width, bins),
                        ];
                        for (acc, x) in g.iter_mut().zip(parts) {
                            for (y, z) in acc.iter_mut().zip(x) {
                                *y += z;
                            }
                        }
                        frames += 1;
                    }
                }
            }
        }
        let elapsed = t0.elapsed().as_secs_f64();
        if stage.name == "sample" {
            sample_ms = 1e3 * elapsed / steps as f64;
        }
        println!(
            "stage {k} ({}): {} ps at {} fs, {} K, γ {} /ps: {:.0} s, {:.1} ms a step with the records",
            stage.name,
            stage.picoseconds,
            stage.dt / FS,
            stage.temperature,
            stage.friction / 1e12,
            elapsed,
            1e3 * elapsed / steps as f64
        );
    }
    assert!(series.iter().all(|s| s.1.is_finite() && s.2.is_finite()));
    let shake = l.shake();
    let worst = shake
        .bonds()
        .iter()
        .zip(shake.lengths())
        .map(|(&[i, j], d)| (len(sub(at[i], at[j])) - d).abs() / d)
        .fold(0.0f64, f64::max);
    println!("the worst held bond is {worst:.2e} of its length off it");
    assert!(worst < 1e-10);

    println!("\nthe melt, in blocks of 2.5 ps:\n| t (ps) | U/N (kcal/mol) | P (bar) | T (K) | order | trans |\n| --- | --- | --- | --- | --- | --- |");
    for block in series.chunks(25) {
        let m = |f: &dyn Fn(&Row) -> f64| block.iter().map(f).sum::<f64>() / block.len() as f64;
        println!(
            "| {:.1}–{:.1} | {:.3} | {:.0} | {:.1} | {:.3} | {:.3} |",
            block[0].0 - 0.1,
            block[block.len() - 1].0,
            m(&|s| s.1),
            m(&|s| s.2),
            m(&|s| s.3),
            m(&|s| s.4),
            m(&|s| s.5)
        );
    }
    let half = sample_u.len() / 2;
    let report = |name: &str, xs: &[f64], digits: usize| {
        let (e, e1, e2) = (
            Estimate::of(xs),
            Estimate::of(&xs[..half]),
            Estimate::of(&xs[half..]),
        );
        println!(
            "{name} = {:.d$} ± {:.d$} (τ {:.1} samples of 0.1 ps); halves {:.d$} ± {:.d$}, {:.d$} ± {:.d$}",
            e.mean, e.error, e.tau, e1.mean, e1.error, e2.mean, e2.error, d = digits
        );
        e
    };
    println!();
    let u = report("U/N (kcal/mol)", &sample_u, 3);
    let p = report("P (bar)", &sample_p, 0);
    println!("  experiment at this density: 1 bar");
    let s = report("order", &sample_s, 4);
    println!(
        "  an isotropic set of {}: {iso:.4} ± {iso_sd:.4}; the lattice: 1",
        l.octanols
    );
    report("trans fraction", &sample_trans, 4);
    let e2e: Vec<f64> = sample_ends.iter().map(|x| x / ANGSTROM).collect();
    report("end-to-end (Å)", &e2e, 3);
    println!(
        "  the relaxed all-trans conformer: {:.3} Å",
        ends0 / ANGSTROM
    );

    let gas_potential = gas_phase(&a);
    let (ug, ug_err) = gas(&a.molecule, &gas_potential, 200.0);
    let rt = BOLTZMANN.to_si() * KELVIN;
    let dh = ug - u.mean * KCAL_PER_MOL + rt;
    let dh_err = (ug_err * ug_err + (u.error * KCAL_PER_MOL) * (u.error * KCAL_PER_MOL)).sqrt();
    println!(
        "U_gas = {:.3} ± {:.3} kcal/mol (eight runs of 200 ps); ΔH_vap = U_gas − U/N + RT = {:.2} ± \
         {:.2} kcal/mol = {:.1} ± {:.1} kJ/mol\n  experiment: {} kJ/mol at 25 °C ({})",
        kcal(ug),
        kcal(ug_err),
        kcal(dh),
        kcal(dh_err),
        dh / KJ_PER_MOL,
        dh_err / KJ_PER_MOL,
        oil.enthalpy_of_vaporization,
        oil.cited
    );

    for (name, gk) in ["terminal–terminal", "interior–interior", "carbon–carbon"]
        .iter()
        .zip(&g)
    {
        let mean: Vec<f64> = gk.iter().map(|x| x / frames as f64).collect();
        let (peak, min) = first_shell(&mean, width, [3.0 * ANGSTROM, 7.0 * ANGSTROM]);
        let tail: f64 = mean[bins - 40..].iter().sum::<f64>() / 40.0;
        println!(
            "g(r) {name}, intermolecular, {frames} frames: first peak at {:.2} Å, height {:.2}; \
             first minimum at {:.2} Å, {:.2}; mean over the last 2 Å {:.3}",
            peak[0] / ANGSTROM,
            peak[1],
            min[0] / ANGSTROM,
            min[1],
            tail
        );
    }

    // The self-diffusion coefficient from the centres' MSD: lags in frames of 0.1 ps.
    let lags: Vec<usize> = (1..=200).collect();
    let all = msd(&centres_of, &lags, None);
    for &t in &[10usize, 20, 50, 100, 150, 200] {
        let k = t - 1;
        let slope = (all[k].ln() - all[k / 2].ln()) / 2f64.ln();
        println!(
            "MSD at {:.0} ps: {:.3} Å², log-log slope from {:.1} ps {:.3}",
            t as f64 / 10.0,
            all[k] / (ANGSTROM * ANGSTROM),
            (k / 2 + 1) as f64 / 10.0,
            slope
        );
    }
    let fit = |m: &[f64]| {
        // Least squares over lags 5–20 ps.
        let pts: Vec<(f64, f64)> = (49..200)
            .map(|k| ((k + 1) as f64 * 0.1e-12, m[k]))
            .collect();
        let c = pts.len() as f64;
        let (sx, sy) = (
            pts.iter().map(|p| p.0).sum::<f64>(),
            pts.iter().map(|p| p.1).sum::<f64>(),
        );
        let (mx, my) = (sx / c, sy / c);
        let num: f64 = pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
        let den: f64 = pts.iter().map(|p| (p.0 - mx) * (p.0 - mx)).sum();
        num / den / 6.0
    };
    let d_all = fit(&all);
    let groups: Vec<f64> = (0..7)
        .map(|gi| {
            let members: Vec<usize> = (0..l.octanols).filter(|i| i % 7 == gi).collect();
            fit(&msd(&centres_of, &lags, Some(&members)))
        })
        .collect();
    let gm = groups.iter().sum::<f64>() / 7.0;
    let gse = (groups.iter().map(|x| (x - gm) * (x - gm)).sum::<f64>() / 6.0 / 7.0).sqrt();
    let v = l.cell().volume();
    let side = v.cbrt();
    let eta = oil.viscosity.0 * 1e-3;
    let yeh_hummer = 2.837297 * BOLTZMANN.to_si() * KELVIN / (6.0 * PI * eta * side);
    let bath = BOLTZMANN.to_si() * KELVIN / (a.molecule.mass() * 1e12);
    println!(
        "D from the MSD's slope over 5–20 ps: {:.3e} m²/s ({:.3} × 10⁻⁹), seven groups of molecules \
         {gm:.3e} ± {gse:.3e}; Yeh–Hummer's finite-box term at the experiment's {} mPa s ({} °C) and \
         L = V^(1/3) = {:.1} Å: +{:.2e}; the 1 ps⁻¹ bath's own k_BT/(Mγ) {:.2e}",
        d_all,
        d_all * 1e9,
        oil.viscosity.0,
        oil.viscosity.1,
        side / ANGSTROM,
        yeh_hummer,
        bath
    );
    println!(
        "  melting point {:.2} K: {:.1} K below 298.15 K; order over the sample {:.4} against an \
         isotropic {iso:.4} ± {iso_sd:.4}",
        oil.melting,
        KELVIN - oil.melting,
        s.mean
    );
    println!(
        "the sample at {sample_ms:.1} ms a step, records included; the whole run {:.0} s",
        started.elapsed().as_secs_f64()
    );
    let _ = p;
}

/// **Liquid n-hexadecane, measured**: [`measure`] on [`HEXADECANE`].
#[test]
#[ignore = "98 hexadecanes melted and sampled, about two hours with --release: run with --release -- --ignored --nocapture"]
fn the_liquid_measured() {
    measure(&HEXADECANE);
}

/// **Liquid n-decane, measured**, the fallback: [`measure`] on [`DECANE`].
#[test]
#[ignore = "112 decanes melted and sampled, about forty minutes with --release: run with --release -- --ignored --nocapture"]
fn the_decane_liquid_measured() {
    measure(&DECANE);
}

/// **Milliseconds a step** of the 98-hexadecane box at 2 fs, release, one core, after 0.2 ps at
/// 0.5 fs off the lattice: classical Ewald with the cell list, PME with the cell list, and PME with
/// the neighbour list; each the best of three runs of 20 steps, and one evaluation of each alone.
/// And what a 4c-2-style leg would cost at the last: 220 ps a window, thirteen windows.
#[test]
#[ignore = "a timing, release only: run with --release -- --ignored --nocapture"]
fn the_cost_measured() {
    let a = Alkane::new(R16);
    let l = Liquid::lattice(
        &a.molecule,
        HEXADECANE.lattice,
        HEXADECANE.density,
        0xC05,
        9.0 * ANGSTROM,
        1e-5,
    );
    let gaffed = with_gaff(l.field.clone(), &a);
    let p = gaffed.ewald().parameters();
    let mesh = PmeParameters::for_accuracy(&l.cell(), &p, 1e-5, l.at.len());
    println!(
        "{} molecules, {} atoms; α {:.4} /Å, {} wave vectors; PME {:?}",
        l.octanols,
        l.at.len(),
        p.alpha * ANGSTROM,
        gaffed.ewald().wave_vectors(),
        mesh
    );
    let pme = production(gaffed.clone(), &l);
    let (mut md, mut at) = bathed(&l.dynamics(), &pme, &l.at, 0.5 * FS, 400, 0xC05);
    let mut last = 0.0;
    for (name, field) in [
        ("classical Ewald, cell list", gaffed.clone()),
        ("PME, cell list", gaffed.clone().with_mesh(mesh)),
        ("PME, neighbour list", pme.clone()),
    ] {
        let mut eval = f64::INFINITY;
        for _ in 0..3 {
            let s = std::time::Instant::now();
            std::hint::black_box(field.evaluate(&at));
            eval = eval.min(s.elapsed().as_secs_f64());
        }
        md.forget_forces();
        let mut best = f64::INFINITY;
        for _ in 0..3 {
            let s = std::time::Instant::now();
            md.run(&field, &mut at, 2.0 * FS, 20);
            best = best.min(s.elapsed().as_secs_f64() / 20.0);
        }
        println!(
            "{name}: {:.1} ms a step at 2 fs; one evaluation {:.1} ms",
            1e3 * best,
            1e3 * eval
        );
        last = best;
    }
    let window = 110_000.0 * last;
    println!(
        "at the last: a 220 ps window {:.2} h, thirteen {:.1} h, before the couplings' cost",
        window / 3600.0,
        13.0 * window / 3600.0
    );
}
