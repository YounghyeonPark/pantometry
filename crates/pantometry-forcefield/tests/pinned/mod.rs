//! **Two systems whose default evaluation is pinned to the bits HEAD `b6558df` computed**, before
//! the non-bonded overrides existed: so a refactor of the default path that the new code shares
//! with itself — the pair loop, the tail, the exclusions, the mesh — cannot pass by agreeing with
//! itself. This file compiles unchanged against that commit's API, and was run there to write the
//! pins (`the_nonbonded_overrides_against_closed_forms.rs`'s `HEAD_WATER` and `HEAD_OCTANOL`).
//!
//! **Every input is the same bits on every platform**, so CI's Linux, macOS and WebAssembly jobs
//! can hold the pins: positions from decimal literals, the dictionary's ideal coordinates and the
//! kernel's [`Rng`], never a minimiser or the platform's `sin`, `cos` or `cbrt`; charges TIP3P's
//! and FreeSolv's AM1-BCC literals; α, the cutoff and the mesh fixed, not chosen by
//! `for_accuracy`; and the reciprocal part by the mesh, whose twiddles are its own. Water and an
//! ion have no bonded term, so their forces and virial are pinned too; the octanols' bonded terms
//! read the platform's `cos`, so only their non-bonded energy terms are.

#![allow(dead_code)]

use crate::freesolv;
use crate::octanol::joined;
use pantometry_core::Rng;
use pantometry_forcefield::uff;
use pantometry_forcefield::{
    Component, EwaldParameters, Ion, PeriodicBox, PeriodicEnergy, PeriodicEvaluation,
    PeriodicForceField, PmeParameters,
};

pub const ANGSTROM: f64 = 1e-10;

/// α, cutoff and an unused `k_c` from a cutoff, by `+ − × ÷` alone.
fn parameters(cutoff: f64) -> EwaldParameters {
    EwaldParameters {
        alpha: 3.2 / cutoff,
        cutoff,
        k_cutoff: 1.0 / cutoff,
    }
}

/// **64 TIP3P waters and one chloride** in a 12 Å cube: the waters on a 4 × 4 × 4 grid of 3 Å, each
/// moved by up to ±0.3 Å an axis by `Rng::for_index`, all oriented alike, the hydrogens at decimal
/// offsets from the oxygen; the ion at a literal point between sites. `r_c` = 5.5 Å, the mesh of
/// order 5 on 16³. The plain field: no override, geometric mixing, every 1-4 pair (there are none)
/// a pair.
pub fn waters_and_an_ion() -> (PeriodicForceField, Vec<[f64; 3]>) {
    waters_and_an_ion_at(5.5 * ANGSTROM)
}

/// [`waters_and_an_ion`] with the cutoff `cutoff`.
pub fn waters_and_an_ion_at(cutoff: f64) -> (PeriodicForceField, Vec<[f64; 3]>) {
    let cell = PeriodicBox::cubic(12.0 * ANGSTROM);
    let field = PeriodicForceField::tip3p(cell, parameters(cutoff), 64)
        .with_ions(&[Ion::chloride()])
        .with_mesh(PmeParameters {
            order: 5,
            grid: [16, 16, 16],
        });
    let mut at = Vec::new();
    for k in 0..64u64 {
        let g = [k / 16, (k / 4) % 4, k % 4];
        let mut rng = Rng::for_index(0x9179, k);
        let o: [f64; 3] =
            [0, 1, 2].map(|a| (1.5 + 3.0 * g[a] as f64 + 0.6 * (rng.unit() - 0.5)) * ANGSTROM);
        at.push(o);
        at.push([o[0] + 0.9572 * ANGSTROM, o[1], o[2]]);
        at.push([
            o[0] - 0.2399872 * ANGSTROM,
            o[1] + 0.9266272 * ANGSTROM,
            o[2],
        ]);
    }
    at.push([3.0 * ANGSTROM, 3.0 * ANGSTROM, 4.5 * ANGSTROM]);
    (field, at)
}

/// **Nine octanols at the dictionary's ideal coordinates**, `OC9` as written and only translated,
/// on a 1 × 3 × 3 lattice of the conformer's extent plus 3 Å, with FreeSolv's AM1-BCC charges;
/// `r_c` = 0.45 of the shortest edge, the mesh of order 5 on 32 × 32 × 32. The plain field.
pub fn ideal_octanols() -> (PeriodicForceField, Vec<[f64; 3]>) {
    ideal_octanols_at(0.45)
}

/// [`ideal_octanols`] with the cutoff `fraction` of the shortest edge.
pub fn ideal_octanols_at(fraction: f64) -> (PeriodicForceField, Vec<[f64; 3]>) {
    let (c, at, charges, cell, shortest) = octanol_lattice();
    let field =
        PeriodicForceField::new(&c, &uff::assign(&c), cell, parameters(fraction * shortest))
            .expect("octanol")
            .with_charges(charges)
            .with_mesh(PmeParameters {
                order: 5,
                grid: [32, 32, 32],
            });
    (field, at)
}

/// [`ideal_octanols`]' component, positions, charges, box and shortest edge.
fn octanol_lattice() -> (Component, Vec<[f64; 3]>, Vec<f64>, PeriodicBox, f64) {
    let one = Component::from_ccd(crate::octanol::OC9).expect("OC9");
    let (elements, bonds) = freesolv::component_graph(&one);
    let q = freesolv::charges_onto(&freesolv::octanol(), &elements, &bonds)
        .expect("the same molecule")
        .charges;
    let ideal: Vec<[f64; 3]> = one.atoms().iter().map(|a| a.at).collect();
    let lo = [0, 1, 2].map(|k| ideal.iter().fold(f64::INFINITY, |m, p| m.min(p[k])));
    let hi = [0, 1, 2].map(|k| ideal.iter().fold(f64::NEG_INFINITY, |m, p| m.max(p[k])));
    let site = [0, 1, 2].map(|k| hi[k] - lo[k] + 3.0 * ANGSTROM);
    let counts = [1usize, 3, 3];
    let cell = PeriodicBox::new([0, 1, 2].map(|k| counts[k] as f64 * site[k]));
    let mut at = Vec::new();
    let mut charges = Vec::new();
    let mut parts = Vec::new();
    for i in 0..counts[0] {
        for j in 0..counts[1] {
            for k in 0..counts[2] {
                let shift = [i, j, k].map(|c| c as f64);
                for p in &ideal {
                    at.push([0, 1, 2].map(|a| p[a] - lo[a] + 1.5 * ANGSTROM + shift[a] * site[a]));
                }
                charges.extend_from_slice(&q);
                parts.push(&one);
            }
        }
    }
    let shortest = site[0]
        .min(counts[1] as f64 * site[1])
        .min(counts[2] as f64 * site[2]);
    (joined(&parts), at, charges, cell, shortest)
}

/// The non-bonded energy terms, and both long-range corrections, as bits.
pub fn nonbonded_bits(field: &PeriodicForceField, e: &PeriodicEnergy) -> Vec<u64> {
    let w = &e.ewald;
    [
        e.van_der_waals,
        e.dispersion_correction,
        e.electrostatic,
        w.real,
        w.reciprocal,
        w.self_energy,
        w.excluded,
        w.background,
        w.total,
        field.dispersion_correction(),
        field.pressure_correction(),
    ]
    .map(f64::to_bits)
    .to_vec()
}

/// **Both long-range corrections at 256 cutoffs**, from 0.30 to 0.45 of the shortest edge in
/// even steps, of both systems: the tail's arithmetic rounded 256 independent ways each. **Few
/// cutoffs do not suffice**: `p/(9 r_c⁹)` rearranged as `p/9/r_c⁹` rounded alike for both
/// pinned systems, and at sixteen cutoffs each in one value of 64. A rearrangement that rounds
/// alike at every one of these is one the pin cannot see, and no longer a likely one.
pub fn tails_over_cutoffs() -> Vec<u64> {
    let (c, _, _, cell, shortest) = octanol_lattice();
    let types = uff::assign(&c);
    let water_cell = PeriodicBox::cubic(12.0 * ANGSTROM);
    let mut out = Vec::new();
    for k in 0..256 {
        let fraction = 0.30 + 0.15 * k as f64 / 256.0;
        let water =
            PeriodicForceField::tip3p(water_cell, parameters(fraction * 12.0 * ANGSTROM), 64)
                .with_ions(&[Ion::chloride()]);
        let octanols = PeriodicForceField::new(&c, &types, cell, parameters(fraction * shortest))
            .expect("octanol");
        for f in [&water, &octanols] {
            out.push(f.dispersion_correction().to_bits());
            out.push(f.pressure_correction().to_bits());
        }
    }
    out
}

/// FNV-1a, 64-bit, over each value's eight bytes, little-endian: a digest of
/// [`tails_over_cutoffs`]' 1 024 values, short enough to pin.
pub fn digest(values: &[u64]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for v in values {
        for b in v.to_le_bytes() {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    h
}

/// [`nonbonded_bits`], then the total, the forces on atoms 0, 1, 2, 100 and the ion (the last),
/// and the virial: for a system with no bonded term.
pub fn all_bits(field: &PeriodicForceField, ev: &PeriodicEvaluation) -> Vec<u64> {
    let mut out = nonbonded_bits(field, &ev.energy);
    out.push(ev.energy.total.to_bits());
    let last = ev.forces.len() - 1;
    for i in [0, 1, 2, 100, last] {
        out.extend(ev.forces[i].map(f64::to_bits));
    }
    out.extend(ev.virial.iter().flatten().map(|x| x.to_bits()));
    out
}
