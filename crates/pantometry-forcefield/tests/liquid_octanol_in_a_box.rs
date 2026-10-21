//! **Liquid 1-octanol in a periodic box**, step 4c-1 of the explicit-solvent track: the solvent in
//! which benzene's solvation free energy will be measured (4c-2), so that
//! `log P = (ΔG_hyd − ΔG_solv,oct) / (RT ln 10)` can be set beside benzene's 2.13. The molecule is
//! the dictionary's `OC9`, UFF throughout, its charges QEq's; dry octanol, no water in it.
//!
//! Default, unoptimised, on boxes of 16 molecules:
//!
//! - **the box is at its density to rounding, with no overlaps**: the mass over the volume against
//!   the density it was built for, and the closest contact of any two atoms of different molecules,
//!   or of a molecule and its own image, by brute force over the 27 images ([`octanol::closest_contact`]);
//! - **the typing is the molecule's**: `C_3`, `O_3` and `H_`, and the whole box typed as one
//!   component is each molecule typed alone;
//! - **every held bond is at its length**, and every X–H — counted from the elements, not from the
//!   constraints — is held, after constrained dynamics at 2 fs;
//! - **the degrees of freedom are `3N − N_held − 3`** in NVE with the net momentum out, and
//!   `3N − N_held` in a bath;
//! - **constrained NVE falls as `h²`** on a cluster of four octanols, with the excess over four
//!   shrinking by four per halving, each in a band [`the_h_squared_band_measured`] measured;
//! - **benzene decoupled from a small octanol box**: the end states are the systems they say, the
//!   long-range correction's cross part is its double sum over three octanol types, and `∂U/∂λ_e`
//!   and `∂U/∂λ_v` are the whole energy's central differences — W3's machinery in this solvent.
//!
//! Ignored, release: [`the_cost_measured`], and [`the_liquid_measured`] — 192 molecules at
//! 0.8262 g/cm³ melted and sampled, its pressure at that density, `U/N`, `ΔH_vap` against a gas of
//! single molecules sampled in vacuum, and `g_OO` — every number reported, none asserted; and
//! [`the_liquid_with_am1bcc_charges_measured`], the same with FreeSolv's charges.
//!
//! # What the liquid measured: UFF's octanol is too cohesive and has no hydrogen bonds
//!
//! Run once each on 2026-10-21, release, one core, 25 000 steps of 2 fs; nothing asserted:
//!
//! | | (a) UFF, QEq | (c) UFF, AM1-BCC | experiment |
//! | --- | --- | --- | --- |
//! | P at 0.8262 g/cm³ (bar) | −930 ± 45 | −788 ± 56 | 1 |
//! | ΔH_vap (kJ/mol) | 91.9 ± 0.9 | 97.3 ± 0.7, still rising | 70.98 (CRC, read via PubChem) |
//! | U/N in the sample's halves (kcal/mol) | 18.00, 17.89 | −1.23, −1.59 | |
//! | g_OO's first peak | 3.39 Å, 0.82; 0.06 O inside its minimum | 3.37 Å, 2.47; 1.40 O | an alcohol's O–H···O near 2.7 Å (quoted) |
//!
//! - **Both are 21–26 kJ/mol too cohesive and pull at −800 to −900 bar** at experiment's density,
//!   which a barostat would turn into a denser liquid. Swapping QEq's charges (±0.16 e on every
//!   C–H) for AM1-BCC's (±0.04) makes it *worse*, so it is not the charges.
//! - **Neither has a hydrogen bond.** No O–O contact closer than 3.3 Å in either; UFF gives the
//!   hydroxyl's hydrogen the van der Waals of every hydrogen, `x` = 2.886 Å, 0.044 kcal/mol, and its
//!   pair with an oxygen is `x` = 3.18 Å, repulsive at the 1.8 Å an H···O bond needs. GAFF's `ho` has
//!   none. So the cohesion is all van der Waals: summed over the molecule, UFF's `C₆` between two
//!   octanols is 1.28 times GAFF's, and its hydrogens' wells are 2.8 times GAFF's `hc`.
//! - **The temperature settles in 2.5 ps after the hot stage; U/N and P do not, in (c)**: its U/N
//!   falls 0.4 kcal/mol between the sample's halves, and τ is 2 ps. (a)'s halves agree to 0.11.
//!   Neither was melted long enough for an octanol whose viscosity is eight times water's.

// The sums over atoms and axes are written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

mod freesolv;
mod hydration;
mod octanol;

use octanol::*;
use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::{
    Alchemical, AtLambda, Bath, Estimate, EwaldParameters, Lambda, MolecularDynamics, PeriodicBox,
    PeriodicDecoupling, PeriodicEnergy, PeriodicForceField, PmeParameters, Shake, UffType,
};
use pantometry_units::BOLTZMANN;
use std::f64::consts::PI;

const EPS: f64 = f64::EPSILON;
const KELVIN: f64 = 298.15;
/// The small boxes' lattice: one molecule along its axis, four by four across.
const SMALL: [usize; 3] = [1, 4, 4];
/// The small boxes' cutoff: under half their shortest edge, 6.96 Å.
const SMALL_CUTOFF: f64 = 6.5 * ANGSTROM;
/// The measurement's lattice: 192 molecules, 41.6 × 36.4 × 33.6 Å.
const LIQUID: [usize; 3] = [3, 8, 8];
/// One kJ/mol in joules per molecule: the crate's `KCAL_PER_MOL` over the thermochemical 4.184.
const KJ_PER_MOL: f64 = KCAL_PER_MOL / 4.184;

fn kcal(j: f64) -> f64 {
    j / KCAL_PER_MOL
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn len(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// `x` to the whole power `n` by repeated multiplication: `powi` is not used for any input or
/// reference, its precision being unspecified and platform-dependent.
fn pow(x: f64, n: u32) -> f64 {
    (0..n).fold(1.0, |p, _| p * x)
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

/// The largest force component, for a scale.
fn largest_force(f: &[[f64; 3]]) -> f64 {
    f.iter().flatten().fold(0.0f64, |m, x| m.max(x.abs()))
}

/// A central difference and the tolerance it earns: three times the truncation term is
/// `|D(2h) − D(h)|`, and the energy's rounding, at most `64 ε` of the sum of its parts'
/// magnitudes, divided by `h`. W1's and W3's.
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

/// The small box: 16 octanols at the liquid's density, `r_c` = 6.5 Å, δ = 10⁻⁵.
fn small_box(o: &Octanol, seed: u64) -> Liquid {
    Liquid::lattice(o, SMALL, DENSITY, seed, SMALL_CUTOFF, 1e-5)
}

/// `steps` of constrained BAOAB at `dt` in a 5 ps⁻¹ bath at 298.15 K, from rest-thermalised
/// velocities, on `potential`; the dynamics and the positions.
fn bathed(
    template: &MolecularDynamics,
    potential: &impl pantometry_forcefield::Potential,
    start: &[[f64; 3]],
    dt: f64,
    steps: usize,
    seed: u64,
) -> (MolecularDynamics, Vec<[f64; 3]>) {
    let mut at = start.to_vec();
    let mut md = template
        .clone()
        .with_bath(Bath::Langevin {
            temperature: KELVIN,
            friction: 5e12,
            seed,
        })
        .thermalised(&at, KELVIN, seed);
    md.run(potential, &mut at, dt, steps);
    (md, at)
}

// ---------------------------------------------------------------------------------------------
// The molecule and the box

/// **The molecule is the one it says, typed as one would**: `OC9` is C₈H₁₈O with every hydrogen
/// on one bond, typed eight `C_3`, one `O_3` and eighteen `H_`; QEq's charges sum to zero; and the
/// whole small box typed as one component is each molecule's types in turn — typing is local, so
/// a liquid typed molecule by molecule is the liquid typed. FreeSolv's AM1-BCC charges map onto it
/// by element and bonding, equal under every one of its automorphisms.
#[test]
fn octanol_is_the_molecule_it_says() {
    let o = Octanol::new(None);
    let c = &o.component;
    assert_eq!(c.formula(), Some("C8 H18 O"));
    let count = |t: UffType| o.types.iter().filter(|&&x| x == t).count();
    assert_eq!(
        (count(UffType::C3), count(UffType::O3), count(UffType::H)),
        (8, 1, 18)
    );
    let q: f64 = o.charges.iter().sum();
    let scale: f64 = o.charges.iter().map(|x| x.abs()).sum();
    println!("QEq: Σq {q:+.2e}, Σ|q| {scale:.4}");
    assert!(q.abs() <= 64.0 * EPS * scale);

    let l = small_box(&o, 1);
    assert_eq!(
        uff::assign(&l.component),
        l.types,
        "typed whole, typed by molecule"
    );
    assert_eq!(l.octanols, 16);
    assert_eq!(l.component.atoms().len(), 16 * 27);

    let mol2 = freesolv::octanol();
    let (elements, bonds) = freesolv::component_graph(c);
    let mapped = freesolv::charges_onto(&mol2, &elements, &bonds).expect("the same molecule");
    // CH₃'s three hydrogens turned 3! ways and each CH₂'s two swapped, seven of them: 6 × 2⁷.
    assert_eq!(mapped.automorphisms, 768);
    assert_eq!(
        mapped.spread, 0.0,
        "AM1-BCC's charges are equal by symmetry"
    );
    println!("| atom | UFF | QEq | AM1-BCC |\n| --- | --- | --- | --- |");
    for i in 0..27 {
        println!(
            "| {} | {} | {:+.4} | {:+.4} |",
            c.atoms()[i].name,
            o.types[i].label(),
            o.charges[i],
            mapped.charges[i]
        );
    }
}

/// **The box is at its density to rounding, and nothing overlaps**: the small box and the
/// measurement's, each the molecules' mass over the volume within `8 ε` of [`DENSITY`], and the
/// small box's own atoms' masses within the rounding of their sum; and in the small box the
/// closest contact of two atoms of different
/// molecules, or of an atom and its own molecule's image, by brute force over the 27 images, is at
/// least [`OVERLAP`] and is the lattice's gap to within the conformer's asymmetry.
#[test]
fn the_box_is_at_its_density_with_no_overlaps() {
    let o = Octanol::new(None);
    for counts in [SMALL, LIQUID] {
        let (cell, site, gap) = Liquid::cell_of(&o, counts, DENSITY);
        let n = (counts[0] * counts[1] * counts[2]) as f64;
        let density = n * o.mass() / cell.volume();
        println!(
            "{counts:?}: box {:.3} × {:.3} × {:.3} Å, site {:.3} × {:.3} × {:.3} Å, gap {:.3} Å, \
             density {:.10} g/cm³",
            cell.lengths()[0] / ANGSTROM,
            cell.lengths()[1] / ANGSTROM,
            cell.lengths()[2] / ANGSTROM,
            site[0] / ANGSTROM,
            site[1] / ANGSTROM,
            site[2] / ANGSTROM,
            gap / ANGSTROM,
            density / 1e3
        );
        assert!((density / DENSITY - 1.0).abs() <= 8.0 * EPS, "{counts:?}");
    }
    let l = small_box(&o, 2);
    // Its own 432 masses summed in turn: a sum of `N` positive terms is within `(N − 1) ε/2` of
    // itself, so the bound is `N ε` and the 8 ε of the edges' product and the division.
    let bound = (l.masses.len() as f64 + 8.0) * EPS;
    println!(
        "the small box's own masses over its volume: {:+.2e} of the density off it, bound {bound:.1e}",
        l.density() / DENSITY - 1.0
    );
    assert!((l.density() / DENSITY - 1.0).abs() <= bound);
    let (r, i, j) = closest_contact(&l.at, &l.molecules, &l.cell());
    println!(
        "closest contact {:.4} Å, atoms {} and {} ({} and {}); gap {:.4} Å",
        r / ANGSTROM,
        i,
        j,
        l.component.atoms()[i].name,
        l.component.atoms()[j].name,
        l.gap / ANGSTROM
    );
    assert!(r >= OVERLAP, "an overlap: {:.3} Å", r / ANGSTROM);
    // The gap holds between neighbours along every axis, less what the turns move a hydrogen.
    assert!(
        r >= 0.9 * l.gap && r <= 1.1 * l.gap,
        "{r:e} against {:e}",
        l.gap
    );
}

/// **The brute force sees across every face, and a molecule's own image**: two one-atom molecules
/// 0.05 of an edge inside opposite faces are 0.1 of that edge apart, along each axis in turn, and
/// one molecule alone is its shortest edge from its own image. A contact check without its
/// images would report the first as 0.9 of the edge and the second as no contact at all.
#[test]
fn the_closest_contact_sees_across_every_face_and_its_own_image() {
    let cell = PeriodicBox::new([13.0 * ANGSTROM, 17.0 * ANGSTROM, 11.0 * ANGSTROM]);
    let l = cell.lengths();
    for axis in 0..3 {
        let mut a = [0.5 * l[0], 0.5 * l[1], 0.5 * l[2]];
        let mut b = a;
        a[axis] = 0.05 * l[axis];
        b[axis] = 0.95 * l[axis];
        let (r, i, j) = closest_contact(&[a, b], &[0..1, 1..2], &cell);
        assert_eq!((i, j), (0, 1));
        assert!(
            (r - 0.1 * l[axis]).abs() <= 1e-12 * l[axis],
            "axis {axis}: {r:e}"
        );
    }
    let (r, i, j) = closest_contact(
        &[[1e-10, 2e-10, 3e-10]],
        std::slice::from_ref(&(0..1)),
        &cell,
    );
    assert_eq!((i, j), (0, 0));
    assert!(
        (r - 11.0 * ANGSTROM).abs() <= 1e-12 * l[2],
        "its own image at {r:e}"
    );
}

// ---------------------------------------------------------------------------------------------
// The constraints

/// **Every bond to a hydrogen is held, and at its length**: the X–H bonds counted from the
/// elements are the constraints' bonds exactly, eighteen a molecule, each held at its UFF stretch's
/// natural length; and after 40 steps of constrained BAOAB at 2 fs from the lattice, each is
/// within twice SHAKE's default tolerance of that length and has no velocity along itself beyond
/// twice RATTLE's.
#[test]
fn every_held_bond_is_at_its_length() {
    let o = Octanol::new(None);
    let l = small_box(&o, 3);
    let mut from_elements = bonds_to_hydrogen(&l.component);
    assert_eq!(from_elements.len(), 18 * l.octanols);
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
    let natural: std::collections::HashMap<[usize; 2], f64> = l
        .field
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
    let (md, at) = bathed(&l.dynamics(), &l.field, &l.at, 2.0 * FS, 40, 0x4C1);
    let v = md.velocities();
    let (mut worst, mut worst_v) = (0.0f64, 0.0f64);
    for &[i, j] in &from_elements {
        let d = natural[&[i, j]];
        let r = sub(at[i], at[j]);
        let u = sub(v[i], v[j]);
        // In units of the default tolerance, not of `shake.tolerance()`: a looser SHAKE would
        // otherwise loosen its own check, which is how a first version of this test passed one.
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

/// **The degrees of freedom are `3N − N_held − 3`** in NVE with the net momentum taken out, and
/// `3N − N_held` in a bath, with `N_held` counted from the elements: eighteen a molecule.
#[test]
fn the_degrees_of_freedom_are_three_n_less_the_held_less_three() {
    let o = Octanol::new(None);
    let l = small_box(&o, 4);
    let n = l.at.len();
    let held = bonds_to_hydrogen(&l.component).len();
    assert_eq!(held, 18 * 16);
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

/// Four octanols, a 1 × 2 × 2 lattice at 0.6 g/cm³ — its gap 3.0 Å, so that no contact is a clash
/// a 2 fs step cannot start from — turned by `seed`, moved into a box 52 Å across with
/// `r_c` = 26 Å: no pair reaches the cutoff and no image comes inside it, so the potential is
/// smooth and NVE's energy error is the integrator's. Every X–H held, the methyl's three and each
/// methylene's two sharing their carbon.
fn cluster(seed: u64) -> (PeriodicForceField, Vec<[f64; 3]>, MolecularDynamics) {
    let o = Octanol::new(None);
    let l = Liquid::lattice(&o, [1, 2, 2], 600.0, seed, 4.0 * ANGSTROM, 1e-5);
    let side = 52.0 * ANGSTROM;
    let big = PeriodicBox::cubic(side);
    let p = EwaldParameters::for_accuracy(&big, 26.0 * ANGSTROM, 1e-5);
    let field = PeriodicForceField::new(&l.component, &l.types, big, p)
        .expect("the cluster")
        .with_charges(l.charges.clone());
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
    assert!(span < 22.0 * ANGSTROM, "span {span:e}");
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

/// The three ratios of the RMS energy error per halving, 2 → 1, 1 → 0.5 and 0.5 → 0.25 fs, each
/// over 0.1 ps from the same start and velocities, and the RMS error at 2 fs.
fn ratios(seed: u64) -> ([f64; 3], f64) {
    let (field, at, md) = cluster(seed);
    let e = [(2.0, 50), (1.0, 100), (0.5, 200), (0.25, 400)]
        .map(|(h, n)| nve_rms(&field, &at, &md, h * FS, n, seed));
    ([e[0] / e[1], e[1] / e[2], e[2] / e[3]], e[0])
}

/// The bands, **measured** over sixteen starts by [`the_h_squared_band_measured`] (release, run
/// 2026-10-21; the debug run's start is the same bits): the ratio from 2 fs 4.2445–4.2965, from
/// 1 fs 4.0588–4.0696, from 0.5 fs 4.0155–4.0181, and the excess over four shrinking by
/// 0.2587–0.2647 per halving, `(r₃ − 4)/(r₂ − 4)`, which an `h²` integrator with a positive `h⁴`
/// term puts at a quarter. **2 fs is nearer `h²` here than in W3 or W4** (4.229–4.508 and
/// 4.620–5.194 from 2 fs): every hydrogen is held, and what is left stiffest is a C–C stretch and
/// the bends with a hydrogen in them. Each ratio's band is its measured middle with three times
/// its half-spread either side, and the shrink's 3.5 times, as W4's.
const H_SQUARED_FROM_TWO: std::ops::Range<f64> = 4.19..4.35;
const H_SQUARED_FROM_ONE: std::ops::Range<f64> = 4.048..4.081;
const H_SQUARED_FROM_HALF: std::ops::Range<f64> = 4.0129..4.0207;
const H_SQUARED_SHRINKS: std::ops::Range<f64> = 0.2512..0.2722;

/// **Constrained velocity Verlet on four octanols**: every X–H held, NVE over 0.1 ps from 298 K
/// velocities. The RMS energy error falls by four per halving from 2 fs, and its excess over four
/// shrinks by four per halving, each held to its band.
#[test]
fn constrained_nve_on_octanol_falls_as_h_squared() {
    let (r, e2) = ratios(0x0C7);
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
/// velocities from it. `cargo test --release -p pantometry-forcefield --test
/// liquid_octanol_in_a_box the_h_squared_band_measured -- --ignored --nocapture`.
#[test]
#[ignore = "sixteen starts of four NVE runs: run with --release -- --ignored --nocapture"]
fn the_h_squared_band_measured() {
    let mut lo = [f64::INFINITY; 4];
    let mut hi = [f64::NEG_INFINITY; 4];
    for start in 0..16u64 {
        let (r, e2) = ratios(0x0C7 + start);
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
}

// ---------------------------------------------------------------------------------------------
// Benzene decoupled from octanol

/// Benzene (W3's: QEq charges, relaxed, its C–H held) at the centre of a [`SMALL`] lattice, less
/// every octanol with an atom within 2 Å of it, `r_c` = 6.5 Å, δ = 10⁻⁵.
fn benzene_in_octanol(seed: u64) -> (hydration::Benzene, Octanol, Liquid) {
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
        seed,
        SMALL_CUTOFF,
        1e-5,
    );
    (b, o, l)
}

/// The octanols of `l` alone, in the same box with the same Ewald parameters and charges.
fn octanol_alone(l: &Liquid, o: &Octanol) -> PeriodicForceField {
    let parts: Vec<&pantometry_forcefield::Component> =
        (0..l.octanols).map(|_| &o.component).collect();
    let c = joined(&parts);
    PeriodicForceField::new(
        &c,
        &l.types[l.solute..],
        l.cell(),
        l.field.ewald().parameters(),
    )
    .expect("octanol alone")
    .with_charges(l.charges[l.solute..].to_vec())
}

/// The small box after `steps` of constrained BAOAB at 2 fs at `lambda`.
fn moved(l: &Liquid, d: &PeriodicDecoupling, lambda: Lambda, steps: usize) -> Vec<[f64; 3]> {
    let potential = AtLambda {
        hamiltonian: d,
        lambda,
    };
    bathed(&l.dynamics(), &potential, &l.at, 2.0 * FS, steps, 17).1
}

/// **The end states are the systems they say**, at a configuration constrained dynamics left.
/// Fully coupled: the whole box less benzene alone in the same box, plus benzene alone in vacuum.
/// Fully decoupled: the octanols alone in their box plus benzene in vacuum, the energy and every
/// force. Each within `10⁻¹²` of the parts' magnitudes: the same terms, summed in other orders. And
/// the long-range correction's cross part is its double sum written out, every ordered pair of a
/// benzene atom and an octanol atom by UFF's geometric rule — three octanol types, so the sum the
/// force field factorises is not one type's square.
#[test]
fn benzene_in_octanol_end_states_are_the_systems_they_say() {
    let (b, o, l) = benzene_in_octanol(0xE5D);
    println!(
        "benzene and {} octanols in {:.2} × {:.2} × {:.2} Å",
        l.octanols,
        l.cell().lengths()[0] / ANGSTROM,
        l.cell().lengths()[1] / ANGSTROM,
        l.cell().lengths()[2] / ANGSTROM
    );
    assert!(l.octanols < 16 && l.octanols >= 10, "{} kept", l.octanols);
    let group: Vec<bool> = (0..l.at.len()).map(|i| i < l.solute).collect();
    let d = PeriodicDecoupling::new(&l.field, &group).expect("benzene decouples");
    let at = moved(&l, &d, Lambda::new(0.0, 0.5, 0.7), 20);
    let n = at.len();
    let mut f = vec![[0.0; 3]; n];
    let coupled = d.energy_and_forces(&at, Lambda::COUPLED, &mut f);
    let whole = l.field.energy(&at);
    let alone_in_the_box = PeriodicForceField::new(
        &b.component,
        &b.types,
        l.cell(),
        l.field.ewald().parameters(),
    )
    .unwrap()
    .with_charges(b.charges.clone())
    .energy(&at[..12]);
    let vacuum = b.vacuum.evaluate(&at[..12]);
    let expected = whole.total - alone_in_the_box.total + vacuum.energy.total;
    let scale = magnitude(&whole);
    println!(
        "coupled {coupled:.10e} J against {expected:.10e}; benzene's own images {:+.4e} kcal/mol",
        kcal(alone_in_the_box.total - vacuum.energy.total)
    );
    assert!((coupled - expected).abs() <= 1e-12 * scale);

    let rc = l.field.ewald().parameters().cutoff;
    let (rc3, v) = (rc * rc * rc, l.cell().volume());
    let mut tail = 0.0;
    for tb in &b.types {
        for to in &l.types[l.solute..] {
            let (pb, po) = (tb.parameters(), to.parameters());
            let x = (pb.vdw_distance * po.vdw_distance).sqrt();
            let dd = (pb.vdw_energy * po.vdw_energy).sqrt();
            let x3 = x * x * x;
            tail +=
                dd * (x3 * x3 * x3 * x3 / (9.0 * rc3 * rc3 * rc3) - 2.0 * x3 * x3 / (3.0 * rc3));
        }
    }
    let tail = 2.0 * 2.0 * PI / v * tail;
    let theirs = d.cross_dispersion_correction();
    println!(
        "the cross correction {:+.6e} kcal/mol, written out {:+.6e}",
        kcal(theirs),
        kcal(tail)
    );
    assert!((theirs - tail).abs() <= 1e-13 * tail.abs());

    let decoupled = d.energy_and_forces(&at, Lambda::new(0.0, 0.0, 0.0), &mut f);
    let solvent = octanol_alone(&l, &o).evaluate(&at[12..]);
    let expected = solvent.energy.total + vacuum.energy.total;
    println!("decoupled {decoupled:.10e} J against {expected:.10e}");
    assert!((decoupled - expected).abs() <= 1e-12 * scale);
    let fscale = largest_force(&f);
    for i in 0..n {
        let want = if i < 12 {
            vacuum.forces[i]
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

/// **`∂U/∂λ_e` and `∂U/∂λ_v` are the energy's derivatives** in octanol, at configurations
/// constrained dynamics produced at four states — the charges part-way, both part-way, the van der
/// Waals part-way with the charges off, and near the decoupled end, where the soft core is softest
/// — against central differences of the whole energy, each to the tolerance the difference earns.
#[test]
fn benzene_in_octanol_lambda_derivatives_are_the_whole_energys() {
    let (_, _, l) = benzene_in_octanol(0xDE1);
    let group: Vec<bool> = (0..l.at.len()).map(|i| i < l.solute).collect();
    let d = PeriodicDecoupling::new(&l.field, &group).expect("benzene decouples");
    for (lambda, steps) in [
        (Lambda::new(0.0, 0.5, 1.0), 20),
        (Lambda::new(0.0, 0.4, 0.8), 20),
        (Lambda::new(0.0, 0.0, 0.5), 20),
        (Lambda::new(0.0, 0.0, 0.1), 40),
    ] {
        let at = moved(&l, &d, lambda, steps);
        let c = d.coupling(&at, lambda);
        let scale = magnitude(&l.field.energy(&at));
        for (k, h) in [(1usize, 1e-3), (2, 1e-4)] {
            let at_k = lambda.components()[k];
            if at_k < 2.0 * h || at_k > 1.0 - 2.0 * h {
                continue;
            }
            let energy = |e: f64| {
                let mut x = lambda;
                if k == 1 {
                    x.electrostatics += e;
                } else {
                    x.van_der_waals += e;
                }
                let mut f = vec![[0.0; 3]; at.len()];
                (d.energy_and_forces(&at, x, &mut f), scale)
            };
            let (dl, tolerance) = finite_difference(energy, h);
            println!(
                "{lambda:?}: ∂U/∂λ[{k}] = {:+.10e} J, by differences {dl:+.10e} ± {tolerance:.1e}",
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

// ---------------------------------------------------------------------------------------------
// The liquid, measured

/// The pressure of the untruncated fluid at `at`, pascals, by the molecular virial: each
/// molecule's forces taken at its centre of mass, so that the constraint forces, internal to it,
/// drop out; the kinetic part the bath's `N_mol k_B T` exactly; and A&T's tail in place of the
/// virial's. W2's, with flexible molecules: the molecular virial is the atomic one less
/// `Σ_i (r_i − R_mol)·F_i`, on whole molecules, which molecular dynamics keeps them.
fn pressure(field: &PeriodicForceField, l: &Liquid, at: &[[f64; 3]]) -> f64 {
    let ev = field.evaluate(at);
    let mut w = ev.virial[0][0] + ev.virial[1][1] + ev.virial[2][2];
    for r in &l.molecules {
        let m: f64 = l.masses[r.clone()].iter().sum();
        let com = [0, 1, 2].map(|k| r.clone().map(|i| l.masses[i] * at[i][k]).sum::<f64>() / m);
        for i in r.clone() {
            w -= dot(sub(at[i], com), ev.forces[i]);
        }
    }
    let v = field.cell().volume();
    let n = l.molecules.len() as f64;
    (n * BOLTZMANN.to_si() * KELVIN + w / 3.0) / v - field.dispersion_correction() / v
        + field.pressure_correction()
}

/// The O–O histogram of one frame, normalised to `g(r)`, bins of `width` to `bins`.
fn g_oo(
    cell: &PeriodicBox,
    at: &[[f64; 3]],
    oxygens: &[usize],
    width: f64,
    bins: usize,
) -> Vec<f64> {
    let n = oxygens.len();
    let mut h = vec![0.0; bins];
    for a in 0..n {
        for b in a + 1..n {
            let d = cell.minimum_image(sub(at[oxygens[a]], at[oxygens[b]]));
            let k = (len(d) / width) as usize;
            if k < bins {
                h[k] += 1.0;
            }
        }
    }
    let v = cell.volume();
    for (k, x) in h.iter_mut().enumerate() {
        let (r0, r1) = (k as f64 * width, (k + 1) as f64 * width);
        let shell = 4.0 / 3.0 * PI * (pow(r1, 3) - pow(r0, 3));
        *x /= (n * (n - 1)) as f64 / 2.0 * shell / v;
    }
    h
}

/// The measurement's field: PME at δ = 10⁻⁵ by `PmeParameters::for_accuracy`, and the neighbour
/// list at its measured skin.
fn production_field(l: &Liquid) -> PeriodicForceField {
    let p = l.field.ewald().parameters();
    let mesh = PmeParameters::for_accuracy(&l.cell(), &p, 1e-5, l.at.len());
    l.field
        .clone()
        .with_mesh(mesh)
        .with_neighbour_list(PeriodicForceField::NEIGHBOUR_SKIN)
}

/// **Milliseconds a step** of the 192-molecule box at 2 fs, release, one core, after 0.2 ps at
/// 0.5 fs off the lattice: classical Ewald with the cell list, PME with the cell list, and PME with
/// the neighbour list; each the best of three runs of 20 steps, and one evaluation of each alone.
#[test]
#[ignore = "a timing, release only: run with --release -- --ignored --nocapture"]
fn the_cost_measured() {
    let o = Octanol::new(None);
    let l = Liquid::lattice(&o, LIQUID, DENSITY, 0x0C70, 9.0 * ANGSTROM, 1e-5);
    let p = l.field.ewald().parameters();
    let mesh = PmeParameters::for_accuracy(&l.cell(), &p, 1e-5, l.at.len());
    println!(
        "{} molecules, {} atoms; α {:.4} /Å, {} wave vectors; PME {:?}",
        l.octanols,
        l.at.len(),
        p.alpha * ANGSTROM,
        l.field.ewald().wave_vectors(),
        mesh
    );
    let pme = production_field(&l);
    let (mut md, mut at) = bathed(&l.dynamics(), &pme, &l.at, 0.5 * FS, 400, 0xC05);
    for (name, field) in [
        ("classical Ewald, cell list", l.field.clone()),
        ("PME, cell list", l.field.clone().with_mesh(mesh)),
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
    }
}

/// One record of the melt: the time (ps), `U/N` (kcal/mol), the pressure (bar), the temperature (K).
type Row = (f64, f64, f64, f64);

/// One stage of the melt: its length, step, bath and temperature.
struct Stage {
    name: &'static str,
    picoseconds: f64,
    dt: f64,
    friction: f64,
    temperature: f64,
}

/// **The melt and the measurement**: 0.4 ps at 0.5 fs in a 50 ps⁻¹ bath at 298.15 K, which takes
/// out the lattice's strain; then 10 ps at 500 K and 2 fs to lose the lattice's order; then 14 ps
/// at 298.15 K to settle, then 26 ps sampled in a 1 ps⁻¹ bath: 25 000 steps of 2 fs, chosen to fit
/// an hour at the 115 ms a step [`the_cost_measured`] first read. It read 46.5 ms the second time,
/// and the runs took 48–64 ms a step (26 minutes) and 48–116 (44 minutes): the same code, the
/// machine's speed changing under it. **They are short for this liquid**, whose viscosity is 7.3 mPa s at 25 °C (CRC Handbook p. 6-234, as PubChem quotes it),
/// eight times water's, and the melt's blocks are printed so that what has not settled shows.
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
        picoseconds: 10.0,
        dt: 2.0 * FS,
        friction: 5e12,
        temperature: 500.0,
    },
    Stage {
        name: "settle",
        picoseconds: 14.0,
        dt: 2.0 * FS,
        friction: 5e12,
        temperature: KELVIN,
    },
    Stage {
        name: "sample",
        picoseconds: 26.0,
        dt: 2.0 * FS,
        friction: 1e12,
        temperature: KELVIN,
    },
];

/// `U` of one molecule in vacuum at 298.15 K, joules: eight runs of `ps_each` after 10 ps each,
/// every X–H held, a 5 ps⁻¹ bath; the mean over the runs and its standard error from their spread.
fn gas(o: &Octanol, ps_each: f64) -> (f64, f64) {
    let means: Vec<f64> = (0..8u64)
        .map(|s| {
            let template =
                MolecularDynamics::new(o.masses.clone()).with_bond_constraints(o.shake());
            let (mut md, mut at) = bathed(&template, &o.vacuum, &o.at, 2.0 * FS, 5000, 0x6A5 + s);
            let steps = (ps_each * 1e-12 / (2.0 * FS)).round() as usize;
            let mut sum = 0.0;
            let mut count = 0usize;
            for k in 0..steps {
                md.step(&o.vacuum, &mut at, 2.0 * FS);
                if k % 50 == 0 {
                    sum += md.potential_energy().expect("stepped");
                    count += 1;
                }
            }
            sum / count as f64
        })
        .collect();
    let m = means.iter().sum::<f64>() / 8.0;
    let var = means.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / 7.0;
    (m, (var / 8.0).sqrt())
}

/// **The liquid, measured**, model (a): UFF with QEq's charges. [`LIQUID`]'s 192 molecules at
/// [`DENSITY`], `r_c` = 9 Å, PME at δ = 10⁻⁵ with the neighbour list, every X–H held, 2 fs; melted
/// by [`PLAN`]. Reported: every 0.1 ps of the melt, `U/N`, the pressure and the temperature in
/// 2.5 ps blocks, to show how long each takes to settle; over the sample, `U/N`, `P` with their
/// autocorrelation-aware errors and the sample's halves, `ΔH_vap = U_gas − U_liq/N + RT` with the
/// gas sampled in vacuum, and `g_OO`'s first peak and minimum. Nothing about the model is
/// asserted: only that the box has no overlap, every held bond stays held, and the run is finite.
#[test]
#[ignore = "192 octanols melted and sampled, about an hour with --release: run with --release -- --ignored --nocapture"]
fn the_liquid_measured() {
    liquid("UFF, QEq's charges", Octanol::new(None));
}

/// **The same with FreeSolv's AM1-BCC charges** on UFF's van der Waals and bonded terms
/// ([`freesolv::MOBLEY_36119`] by element and bonding): the half of model (b) that needs no
/// library change, which says how much of (a) is QEq's charges.
#[test]
#[ignore = "192 octanols melted and sampled, about an hour with --release: run with --release -- --ignored --nocapture"]
fn the_liquid_with_am1bcc_charges_measured() {
    let qeq = Octanol::new(None);
    let (elements, bonds) = freesolv::component_graph(&qeq.component);
    let mapped =
        freesolv::charges_onto(&freesolv::octanol(), &elements, &bonds).expect("the same molecule");
    liquid("UFF, AM1-BCC's charges", Octanol::new(Some(mapped.charges)));
}

/// [`the_liquid_measured`]'s protocol for octanol `o`.
fn liquid(name: &str, o: Octanol) {
    let started = std::time::Instant::now();
    let l = Liquid::lattice(&o, LIQUID, DENSITY, 0x0C71, 9.0 * ANGSTROM, 1e-5);
    let (r, _, _) = closest_contact(&l.at, &l.molecules, &l.cell());
    println!(
        "\n{name}: {} octanols, {} atoms, box {:.3} × {:.3} × {:.3} Å at {:.6} g/cm³; closest contact {:.3} Å",
        l.octanols,
        l.at.len(),
        l.cell().lengths()[0] / ANGSTROM,
        l.cell().lengths()[1] / ANGSTROM,
        l.cell().lengths()[2] / ANGSTROM,
        l.density() / 1e3,
        r / ANGSTROM
    );
    assert!(r >= OVERLAP);
    let field = production_field(&l);
    let n = l.octanols as f64;
    let oxygens: Vec<usize> = l
        .molecules
        .iter()
        .map(|m| m.start + index(&o.component, "OAB"))
        .collect();
    let mut at = l.at.clone();
    let mut md = l.dynamics();
    let mut series: Vec<Row> = Vec::new();
    let mut sample_u = Vec::new();
    let mut sample_p = Vec::new();
    let width = 0.02 * ANGSTROM;
    let bins = (l.cell().half_shortest_edge() / width) as usize;
    let mut g = vec![0.0; bins];
    let mut frames = 0usize;
    let mut time = 0.0;
    for (k, stage) in PLAN.iter().enumerate() {
        md = md.with_bath(Bath::Langevin {
            temperature: stage.temperature,
            friction: stage.friction,
            seed: 0x0C10 + k as u64,
        });
        if k == 0 {
            md = md.thermalised(&at, stage.temperature, 0x0C10);
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
                series.push((time * 1e12, u, p, md.half_step_temperature()));
                if stage.name == "sample" {
                    sample_u.push(u);
                    sample_p.push(p);
                    if (s + 1) % (5 * every) == 0 {
                        for (x, y) in g
                            .iter_mut()
                            .zip(g_oo(&l.cell(), &at, &oxygens, width, bins))
                        {
                            *x += y;
                        }
                        frames += 1;
                    }
                }
            }
        }
        let elapsed = t0.elapsed().as_secs_f64();
        println!(
            "stage {} ({}): {} ps at {} fs, {} K, γ {} /ps: {:.0} s, {:.1} ms a step",
            k,
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

    println!("\nthe melt, in blocks of 2.5 ps:\n| t (ps) | U/N (kcal/mol) | P (bar) | T (K) |\n| --- | --- | --- | --- |");
    for block in series.chunks(25) {
        let m = |f: &dyn Fn(&Row) -> f64| block.iter().map(f).sum::<f64>() / block.len() as f64;
        println!(
            "| {:.1}–{:.1} | {:.3} | {:.0} | {:.1} |",
            block[0].0 - 0.1,
            block[block.len() - 1].0,
            m(&|s| s.1),
            m(&|s| s.2),
            m(&|s| s.3)
        );
    }
    let u = Estimate::of(&sample_u);
    let p = Estimate::of(&sample_p);
    // The sample's halves, as a check on drift the autocorrelation cannot see.
    let half = sample_u.len() / 2;
    let (u1, u2) = (
        Estimate::of(&sample_u[..half]),
        Estimate::of(&sample_u[half..]),
    );
    let (p1, p2) = (
        Estimate::of(&sample_p[..half]),
        Estimate::of(&sample_p[half..]),
    );
    println!(
        "\nU/N = {:.3} ± {:.3} kcal/mol (τ {:.1} samples of 0.1 ps); halves {:.3} ± {:.3}, {:.3} ± {:.3}",
        u.mean, u.error, u.tau, u1.mean, u1.error, u2.mean, u2.error
    );
    println!(
        "P = {:.0} ± {:.0} bar (τ {:.1}); halves {:.0} ± {:.0}, {:.0} ± {:.0}; experiment at this density: 1 bar",
        p.mean, p.error, p.tau, p1.mean, p1.error, p2.mean, p2.error
    );
    let (ug, ug_err) = gas(&o, 200.0);
    let rt = BOLTZMANN.to_si() * KELVIN;
    let u_liq = u.mean * KCAL_PER_MOL;
    let dh = ug - u_liq + rt;
    let dh_err = (ug_err * ug_err + (u.error * KCAL_PER_MOL) * (u.error * KCAL_PER_MOL)).sqrt();
    println!(
        "U_gas = {:.3} ± {:.3} kcal/mol (eight runs of 200 ps); ΔH_vap = U_gas − U/N + RT = {:.2} \
         ± {:.2} kcal/mol = {:.1} ± {:.1} kJ/mol",
        kcal(ug),
        kcal(ug_err),
        kcal(dh),
        kcal(dh_err),
        dh / KJ_PER_MOL,
        dh_err / KJ_PER_MOL
    );
    println!("  experiment: 70.98 kJ/mol at 25 °C (CRC Handbook 95th ed. p. 6-142, as PubChem CID 957 quotes it, read)");
    let mean_g: Vec<f64> = g.iter().map(|x| x / frames as f64).collect();
    let lo = (2.3 * ANGSTROM / width) as usize;
    let hi = (3.4 * ANGSTROM / width) as usize;
    let b = (lo..hi)
        .max_by(|&i, &j| mean_g[i].total_cmp(&mean_g[j]))
        .expect("bins");
    let (y0, y1, y2) = (mean_g[b - 1], mean_g[b], mean_g[b + 1]);
    let shift = 0.5 * (y0 - y2) / (y0 - 2.0 * y1 + y2);
    let r_peak = (b as f64 + 0.5 + shift) * width;
    let min_hi = (4.5 * ANGSTROM / width) as usize;
    let m = (b..min_hi)
        .min_by(|&i, &j| mean_g[i].total_cmp(&mean_g[j]))
        .expect("bins");
    let second = (m..((6.5 * ANGSTROM / width) as usize).min(bins))
        .max_by(|&i, &j| mean_g[i].total_cmp(&mean_g[j]))
        .expect("bins");
    // The coordination number to the first minimum: the oxygens a hydroxyl's neighbours are.
    let rho = n / l.cell().volume();
    let coordination: f64 = (0..m)
        .map(|k| {
            let (r0, r1) = (k as f64 * width, (k + 1) as f64 * width);
            mean_g[k] * rho * 4.0 / 3.0 * PI * (pow(r1, 3) - pow(r0, 3))
        })
        .sum();
    println!(
        "g_OO over {frames} frames: first peak at {:.3} Å, height {:.2}; first minimum at {:.2} Å, \
         {:.3}; {:.2} oxygens inside it; the next maximum at {:.2} Å, {:.2}",
        r_peak / ANGSTROM,
        y1 - 0.25 * (y0 - y2) * shift,
        (m as f64 + 0.5) * width / ANGSTROM,
        mean_g[m],
        coordination,
        (second as f64 + 0.5) * width / ANGSTROM,
        mean_g[second]
    );
    println!("the whole run {:.0} s", started.elapsed().as_secs_f64());
}
