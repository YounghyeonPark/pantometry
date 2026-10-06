//! **UFF in a periodic box**: the force field's forces and virial against its energy, its
//! indifference to which image an atom is given in, its limit in a large box, the long-range van
//! der Waals correction against its integral, molecular dynamics in the box chunked any way, and
//! a group decoupled inside the box against the systems its end states are.
//!
//! The water here is UFF's — `O_3` and `H_`, with TIP3P's charges, −0.834 and +0.417 — flexible
//! and in no way a water model: W2 brings the model. What it gives this step is a dense, charged,
//! bonded periodic system whose every molecule can be cut by a face of the box.

// The sums over atoms and axes are written as the formulas' index loops, which read against
// them; iterator chains over two arrays at once would not.
#![allow(clippy::needless_range_loop)]

mod common;

use common::entry;
use pantometry_forcefield::alchemy::SoftCore;
use pantometry_forcefield::energy::Pair;
use pantometry_forcefield::ewald::COULOMB;
use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::{
    Alchemical, AlchemyError, Bath, Component, EwaldParameters, ForceField, Lambda,
    MolecularDynamics, PeriodicBox, PeriodicDecoupling, PeriodicEnergy, PeriodicForceField,
    Potential,
};

const ANGSTROM: f64 = 1e-10;
const EPS: f64 = f64::EPSILON;
const AIN: &str = include_str!("../components/AIN.cif");

fn splitmix(seed: u64, i: u64) -> u64 {
    let mut z = seed.wrapping_add((i + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn uniform(seed: u64, i: u64) -> f64 {
    (splitmix(seed, i) >> 11) as f64 / (1u64 << 53) as f64
}

/// `k³` waters on a grid filling a cube of side `side` Å, each turned by its own fixed angles
/// and moved off its site by up to 0.3 Å: positions in Å, oxygen first, then its two hydrogens.
fn water_sites(k: usize, side: f64, seed: u64) -> Vec<[[f64; 3]; 3]> {
    let spacing = side / k as f64;
    let (r, half) = (0.9572, 0.5 * 104.52f64.to_radians());
    let mut out = Vec::new();
    for i in 0..k * k * k {
        let g = [i / (k * k), (i / k) % k, i % k];
        let n = i as u64;
        let o = [0, 1, 2]
            .map(|a| (g[a] as f64 + 0.5) * spacing + 0.6 * (uniform(seed, 7 * n + a as u64) - 0.5));
        let (t, p, s) = (
            std::f64::consts::PI * uniform(seed, 7 * n + 3),
            2.0 * std::f64::consts::PI * uniform(seed, 7 * n + 4),
            2.0 * std::f64::consts::PI * uniform(seed, 7 * n + 5),
        );
        // An orthonormal frame from the three angles: u the bisector, v in the plane.
        let u = [t.sin() * p.cos(), t.sin() * p.sin(), t.cos()];
        let w0 = if u[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let c = common::cross(u, w0);
        let lc = common::len(c);
        let a = c.map(|x| x / lc);
        let b = common::cross(u, a);
        let v = [0, 1, 2].map(|m| s.cos() * a[m] + s.sin() * b[m]);
        let h1 = [0, 1, 2].map(|m| o[m] + r * (half.cos() * u[m] + half.sin() * v[m]));
        let h2 = [0, 1, 2].map(|m| o[m] + r * (half.cos() * u[m] - half.sin() * v[m]));
        out.push([o, h1, h2]);
    }
    out
}

/// A component of the waters `sites`, named `O1 H1A H1B O2 …`.
fn waters(sites: &[[[f64; 3]; 3]]) -> Component {
    let names: Vec<[String; 3]> = (0..sites.len())
        .map(|i| [format!("O{i}"), format!("H{i}A"), format!("H{i}B")])
        .collect();
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    for (n, s) in names.iter().zip(sites) {
        atoms.push((n[0].as_str(), "O", s[0]));
        atoms.push((n[1].as_str(), "H", s[1]));
        atoms.push((n[2].as_str(), "H", s[2]));
        bonds.push((n[0].as_str(), n[1].as_str(), "SING", false));
        bonds.push((n[0].as_str(), n[2].as_str(), "SING", false));
    }
    entry(&atoms, &bonds)
}

fn tip3p_charges(c: &Component) -> Vec<f64> {
    c.atoms()
        .iter()
        .map(|a| {
            if a.name.starts_with('O') {
                -0.834
            } else {
                0.417
            }
        })
        .collect()
}

fn positions(c: &Component) -> Vec<[f64; 3]> {
    c.atoms().iter().map(|a| a.at).collect()
}

/// A box of `k³` waters, its force field and its positions.
fn water_box(
    k: usize,
    side: f64,
    cutoff: f64,
    accuracy: f64,
    seed: u64,
) -> (PeriodicForceField, Vec<[f64; 3]>) {
    let c = waters(&water_sites(k, side, seed));
    let cell = PeriodicBox::cubic(side * ANGSTROM);
    let p = EwaldParameters::for_accuracy(&cell, cutoff * ANGSTROM, accuracy);
    let field = PeriodicForceField::new(&c, &uff::assign(&c), cell, p)
        .expect("water is supported")
        .with_charges(tip3p_charges(&c));
    (field, positions(&c))
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

/// A central difference and the tolerance it earns: three times the truncation term is
/// `|D(2h) − D(h)|`, and the energy's rounding, at most `64 ε` of the sum of its parts'
/// magnitudes, divided by `h`.
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
/// accepts: a finite difference that moves a pair across the cutoff sees the step, not the
/// derivative.
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
                let d = cell.minimum_image(common::sub(at[i], at[j]));
                nearest = nearest.min((common::len(d) - rc).abs() / rc);
            }
        }
    }
    nearest
}

#[test]
fn the_forces_are_the_negative_gradient_in_the_box() {
    // 27 waters, every pair by brute force (two cells across), and 64, by the cell list.
    for (k, side, cutoff) in [(3, 9.3, 4.6), (4, 12.4, 4.1)] {
        let (field, at) = water_box(k, side, cutoff, 1e-7, 0x3A7E);
        let h = 1e-5 * ANGSTROM;
        let ev = field.evaluate(&at);
        let mut worst = 0.0f64;
        for i in [0, 1, 2, 13, 26, 3 * k * k * k - 1] {
            // The moved atom's pairs stay on their side of the cutoff over ±2h.
            assert!(nearest_to_cutoff(&field, &at, |j| j == i) * cutoff * ANGSTROM > 4.0 * h);
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
                    "{k}³ waters, atom {i} axis {a}: {f:e} against {:e}, tolerance {tolerance:e}",
                    -d
                );
                // The check can see a part in 10⁴ of the force.
                assert!(tolerance < 1e-4 * f.abs(), "{tolerance:e} against {f:e}");
                worst = worst.max((f + d).abs() / tolerance);
            }
        }
        println!("{k}³ waters: worst |F + dE/dx| is {worst:.3} of its tolerance");
    }
}

#[test]
fn the_virial_is_the_strain_derivative_in_the_box() {
    let (field, at) = water_box(4, 12.4, 4.1, 1e-7, 0x7155);
    let h = 1e-6;
    // Every pair stays on its side of the cutoff over a strain of ±2h.
    assert!(nearest_to_cutoff(&field, &at, |_| true) > 4.0 * h);
    let ev = field.evaluate(&at);
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
        let (d, tolerance) = finite_difference(energy, h);
        println!(
            "W_{a}{a} = {:.6e} J, −dU/dε = {:.6e}, tolerance {tolerance:.1e}",
            w[a][a], -d
        );
        assert!((w[a][a] + d).abs() <= tolerance, "axis {a}");
        assert!(tolerance < 1e-6 * w[a][a].abs());
    }
    for a in 0..3 {
        for b in 0..3 {
            assert!(
                (w[a][b] - w[b][a]).abs() <= 1e-10 * w[0][0].abs(),
                "W is symmetric: {} against {}",
                w[a][b],
                w[b][a]
            );
        }
    }
    // The correction's part of the virial is E_tail on the diagonal, the model's own derivative;
    // A&T's P_tail differs from E_tail/V by the impulsive term, (2π/3V²) ΣΣ r_c³ u(r_c), written
    // here from the box's 64 oxygens and 128 hydrogens — negative, since u(r_c) < 0 for each.
    let v = field.cell().volume();
    let rc = field.ewald().parameters().cutoff;
    let (o, h) = (uff::UffType::O3, uff::UffType::H);
    let impulsive: f64 = [
        (o, o, 64.0 * 64.0),
        (o, h, 2.0 * 64.0 * 128.0),
        (h, h, 128.0 * 128.0),
    ]
    .iter()
    .map(|&(a, b, count)| count * rc.powi(3) * Pair::new([0, 1], a, b).energy(rc))
    .sum::<f64>()
        * 2.0
        * std::f64::consts::PI
        / (3.0 * v * v);
    let p_virial = ev.energy.dispersion_correction / v;
    println!(
        "P_tail {:.10e} Pa = E_tail/V {p_virial:.10e} + the impulsive term {impulsive:.10e}",
        field.pressure_correction()
    );
    assert!(impulsive < 0.0);
    assert!(
        (field.pressure_correction() - p_virial - impulsive).abs()
            <= 1e-12 * field.pressure_correction().abs()
    );
}

/// Every atom wrapped into the box on its own, so that molecules are cut by the faces, and every
/// atom moved by its own whole number of box lengths: the same energy and forces, to rounding.
#[test]
fn which_image_an_atom_is_given_in_changes_nothing_but_rounding() {
    let (field, at) = water_box(4, 12.4, 6.2, 1e-7, 0x1FA6);
    let cell = field.cell();
    let l = cell.lengths();
    let base = field.evaluate(&at);
    // Shifted by 0.37 of a box first, off the grid, so that wrapping cuts the molecules near a
    // face.
    let shifted: Vec<[f64; 3]> = at
        .iter()
        .map(|r| [0, 1, 2].map(|a| r[a] + 0.37 * l[a]))
        .collect();
    let wrapped: Vec<[f64; 3]> = shifted.iter().map(|r| cell.wrap(*r)).collect();
    let cut = (0..at.len() / 3)
        .filter(|&m| {
            let o = wrapped[3 * m];
            (1..3).any(|h| common::len(common::sub(wrapped[3 * m + h], o)) > 2.0 * ANGSTROM)
        })
        .count();
    assert!(cut > 5, "only {cut} molecules cut by a face");
    let scattered: Vec<[f64; 3]> = at
        .iter()
        .enumerate()
        .map(|(i, r)| {
            [0, 1, 2]
                .map(|a| r[a] + ((splitmix(0x5C, 3 * i as u64 + a as u64) % 9) as f64 - 4.0) * l[a])
        })
        .collect();
    let scale = magnitude(&base.energy);
    let fmax = base
        .forces
        .iter()
        .flatten()
        .fold(0.0f64, |m, f| m.max(f.abs()));
    for (what, other) in [
        ("translated by 0.37 of a box", shifted),
        ("wrapped atom by atom", wrapped),
        ("scattered over images", scattered),
    ] {
        let ev = field.evaluate(&other);
        let de = (ev.energy.total - base.energy.total).abs();
        assert!(
            de <= 1e-11 * scale,
            "{what}: energy moved by {de:e} of {scale:e}"
        );
        for (f0, f1) in base.forces.iter().zip(&ev.forces) {
            for a in 0..3 {
                assert!(
                    (f0[a] - f1[a]).abs() <= 1e-9 * fmax,
                    "{what}: a force moved"
                );
            }
        }
        // The virial too: its bonded part is Σ r ⊗ F on the molecules made whole, which a cut
        // molecule's raw positions would get wrong by a box length times its bonded forces.
        let wmax = base
            .virial
            .iter()
            .flatten()
            .fold(0.0f64, |m, w| m.max(w.abs()));
        let mut dw = 0.0f64;
        for a in 0..3 {
            for b in 0..3 {
                dw = dw.max((ev.virial[a][b] - base.virial[a][b]).abs());
            }
        }
        assert!(
            dw <= 1e-9 * wmax,
            "{what}: the virial moved by {dw:e} of {wmax:e}"
        );
        println!(
            "{what}: energy moved by {:.1e} of the parts, the virial by {:.1e} of its largest \
             ({cut} molecules cut)",
            de / scale,
            dw / wmax
        );
    }
    // A molecule cut by a face is made whole: its bonded energy is the uncut one's.
    let whole = field.whole(&field.whole(&at));
    assert_eq!(whole, at, "an uncut box is its own whole, bit for bit");
}

/// The force field's Ewald sum is [`pantometry_forcefield::Ewald::evaluate`]'s on the same charges,
/// positions and excluded pairs — which `excluded_pairs_are_left_out_exactly` holds to the exact
/// Coulomb sum of those pairs — with the cutoff at 1.4 Å, inside a water's 1.51 Å H–H, so that a
/// third of the excluded pairs are beyond it and are corrected all the same.
#[test]
fn the_field_corrects_excluded_pairs_beyond_the_cutoff() {
    let (field, at) = water_box(3, 9.3, 1.4, 1e-6, 0xE8C1);
    let rc = field.ewald().parameters().cutoff;
    let cell = field.cell();
    let beyond = field
        .excluded_pairs()
        .iter()
        .filter(|&&[i, j]| common::len(cell.minimum_image(common::sub(at[i], at[j]))) >= rc)
        .count();
    assert_eq!(beyond, 27, "one H–H pair per water beyond the cutoff");
    let mine = field.energy(&at).ewald;
    let reference = field
        .ewald()
        .evaluate(field.charges(), &at, field.excluded_pairs())
        .energy;
    let scale = reference.real.abs()
        + reference.reciprocal.abs()
        + reference.self_energy.abs()
        + reference.excluded.abs();
    println!(
        "excluded pairs beyond r_c: the field's Ewald {:.12e} J against {:.12e}",
        mine.total, reference.total
    );
    assert!((mine.total - reference.total).abs() <= 1e-13 * scale);
    assert!((mine.excluded - reference.excluded).abs() <= 1e-13 * scale);
}

/// Aspirin with neutral partial charges, alone in a box large enough that no atom sees another's
/// image inside the cutoff: its bonded terms are the vacuum force field's to the bit, its van der
/// Waals the vacuum's to rounding, and its electrostatics the vacuum's less the tinfoil dipole
/// term `2π k_e μ²/(3V)`, with the force `(4π k_e/3V) q_i μ` beside it — and what is left after
/// both falls as `L⁻⁵`.
#[test]
fn a_molecule_in_a_large_box_is_the_vacuum_force_field() {
    let c = Component::from_ccd(AIN).expect("AIN parses");
    let types = uff::assign(&c);
    let n = c.atoms().len();
    let raw: Vec<f64> = (0..n).map(|i| 0.3 * (1.7 * i as f64).cos()).collect();
    let mean = raw.iter().sum::<f64>() / n as f64;
    let q: Vec<f64> = raw.iter().map(|x| x - mean).collect();
    let vacuum = ForceField::new(&c, &types)
        .expect("aspirin is supported")
        .with_charges(q.clone());
    let at = positions(&c);
    let v = vacuum.evaluate(&at);
    let mu = [0, 1, 2].map(|a| (0..n).map(|i| q[i] * at[i][a]).sum::<f64>());
    let mu2 = mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2];
    let mut rows = Vec::new();
    for side in [40.0, 60.0, 80.0] {
        let cell = PeriodicBox::cubic(side * ANGSTROM);
        let p = EwaldParameters::for_accuracy(&cell, 0.5 * side * ANGSTROM, 1e-12);
        let field = PeriodicForceField::new(&c, &types, cell, p)
            .expect("aspirin is supported")
            .with_charges(q.clone());
        let e = field.evaluate(&at);
        assert_eq!(e.energy.bond, v.energy.bond);
        assert_eq!(e.energy.angle, v.energy.angle);
        assert_eq!(e.energy.torsion, v.energy.torsion);
        assert_eq!(e.energy.inversion, v.energy.inversion);
        assert!(
            (e.energy.van_der_waals - v.energy.van_der_waals).abs()
                <= 1e-13 * v.energy.van_der_waals.abs(),
            "van der Waals {} against {}",
            e.energy.van_der_waals,
            v.energy.van_der_waals
        );
        let volume = cell.volume();
        let dipole = 2.0 * std::f64::consts::PI * COULOMB * mu2 / (3.0 * volume);
        let left = e.energy.electrostatic + dipole - v.energy.electrostatic;
        // The forces, less the dipole term's own.
        let field_strength = 4.0 * std::f64::consts::PI * COULOMB / (3.0 * volume);
        let (mut off_max, mut tinfoil_max) = (0.0f64, 0.0f64);
        for i in 0..n {
            for a in 0..3 {
                let tinfoil = field_strength * q[i] * mu[a];
                off_max = off_max.max((e.forces[i][a] - v.forces[i][a] - tinfoil).abs());
                tinfoil_max = tinfoil_max.max(tinfoil.abs());
            }
        }
        println!(
            "aspirin in {side} Å: Ewald − vacuum = {:.4e} kcal/mol, −2πμ²/3V = {:.4e}, left \
             {:.3e}; forces: the dipole field's largest {:.3e} N, left {:.3e}; E_tail {:.3e}",
            (e.energy.electrostatic - v.energy.electrostatic) / KCAL_PER_MOL,
            -dipole / KCAL_PER_MOL,
            left / KCAL_PER_MOL,
            tinfoil_max,
            off_max,
            e.energy.dispersion_correction / KCAL_PER_MOL
        );
        rows.push((side, left, dipole, off_max, tinfoil_max));
    }
    // What is left after the dipole term falls as L⁻⁵, the energy's and the forces', each within
    // a quarter from one box to the next; the dipole term itself as L⁻³.
    for w in rows.windows(2) {
        let (a, b) = (w[0], w[1]);
        let expected = (b.0 / a.0).powi(5);
        let energy = a.1 / b.1;
        let force = a.3 / b.3;
        println!(
            "{} → {} Å: left falls by {energy:.3}, forces' by {force:.3}, L⁻⁵ {expected:.3}",
            a.0, b.0
        );
        assert!(
            (energy / expected - 1.0).abs() < 0.25,
            "energy: {energy} against {expected}"
        );
        assert!(
            (force / expected - 1.0).abs() < 0.25,
            "forces: {force} against {expected}"
        );
        assert!(
            a.1.abs() < 0.25 * a.2 && a.3 < 0.25 * a.4,
            "the dipole term dominates"
        );
    }
}

/// UFF's pair is Lennard-Jones's: `D[(x/r)¹² − 2(x/r)⁶] = 4ε[(σ/r)¹² − (σ/r)⁶]` with `ε = D` and
/// `σ = x/2^⅙`, to rounding. Its long-range correction is the integral of the pair from `r_c` to
/// infinity over the box's ordered pairs, by quadrature; and for one type, Allen and Tildesley's
/// printed forms for the energy and the pressure.
#[test]
fn the_long_range_correction_is_its_integral() {
    let types = [uff::UffType::O3, uff::UffType::H];
    for a in types {
        for b in types {
            let p = Pair::new([0, 1], a, b);
            let (eps, sigma) = (p.well, p.distance / 2f64.powf(1.0 / 6.0));
            for r in [0.8, 1.0, 1.3, 2.0, 3.5] {
                let r = r * p.distance;
                let lj = 4.0 * eps * ((sigma / r).powi(12) - (sigma / r).powi(6));
                let size = 4.0 * eps * ((sigma / r).powi(12) + (sigma / r).powi(6));
                assert!(
                    (p.energy(r) - lj).abs() <= 8.0 * EPS * size,
                    "{a:?}–{b:?} at {r:e}"
                );
            }
        }
    }
    // The integral ∫_{r_c}^∞ r² u(r) dr and ∫ r³ u′(r) dr by Simpson's rule in s = r_c/r.
    let integral = |p: &Pair, rc: f64, power: i32, derivative: bool| -> f64 {
        let n = 4000;
        let mut sum = 0.0;
        for k in 0..=n {
            let s = k as f64 / n as f64;
            let w = if k == 0 || k == n {
                1.0
            } else if k % 2 == 1 {
                4.0
            } else {
                2.0
            };
            if s == 0.0 {
                continue; // the integrand vanishes as s → 0
            }
            let r = rc / s;
            let x = p.distance / r;
            let x6 = x.powi(6);
            let value = if derivative {
                12.0 * p.well * (x6 - x6 * x6) / r
            } else {
                p.energy(r)
            };
            // dr = −r_c/s² ds
            sum += w * value * r.powi(power) * rc / (s * s);
        }
        sum / (3.0 * n as f64)
    };
    let (field, _) = water_box(4, 12.4, 6.0, 1e-6, 0x7A11);
    let rc = field.ewald().parameters().cutoff;
    let volume = field.cell().volume();
    let (n_o, n_h) = (64.0, 128.0);
    let pairs = [
        (types[0], types[0], n_o * n_o),
        (types[0], types[1], 2.0 * n_o * n_h),
        (types[1], types[1], n_h * n_h),
    ];
    let mut energy = 0.0;
    let mut pressure = 0.0;
    for (a, b, count) in pairs {
        let p = Pair::new([0, 1], a, b);
        energy += count * integral(&p, rc, 2, false);
        pressure += count * integral(&p, rc, 3, true);
    }
    energy *= 2.0 * std::f64::consts::PI / volume;
    pressure *= -2.0 * std::f64::consts::PI / (3.0 * volume * volume);
    println!(
        "E_tail {:.10e} J against the quadrature's {energy:.10e}; P_tail {:.10e} Pa against \
         {pressure:.10e}",
        field.dispersion_correction(),
        field.pressure_correction()
    );
    // Simpson's error at 4000 intervals on a smooth integrand: below 1e-12 relative.
    assert!((field.dispersion_correction() / energy - 1.0).abs() < 1e-10);
    assert!((field.pressure_correction() / pressure - 1.0).abs() < 1e-10);

    // One type: a box of oxygens alone, as O₂, against A&T's printed formulas.
    let sites: Vec<[[f64; 3]; 2]> = (0..27)
        .map(|i| {
            let o = [
                (i / 9) as f64 * 3.0 + 1.0,
                ((i / 3) % 3) as f64 * 3.0 + 1.0,
                (i % 3) as f64 * 3.0 + 1.0,
            ];
            [o, [o[0] + 1.21, o[1], o[2]]]
        })
        .collect();
    let names: Vec<[String; 2]> = (0..27)
        .map(|i| [format!("O{i}A"), format!("O{i}B")])
        .collect();
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    for (nm, s) in names.iter().zip(&sites) {
        atoms.push((nm[0].as_str(), "O", s[0]));
        atoms.push((nm[1].as_str(), "O", s[1]));
        bonds.push((nm[0].as_str(), nm[1].as_str(), "DOUB", false));
    }
    let c = entry(&atoms, &bonds);
    let t = uff::assign(&c);
    assert!(t.iter().all(|&x| x == t[0]), "one type");
    let cell = PeriodicBox::cubic(9.0 * ANGSTROM);
    let rc = 4.5 * ANGSTROM;
    let f = PeriodicForceField::new(&c, &t, cell, EwaldParameters::for_accuracy(&cell, rc, 1e-6))
        .expect("O₂ is supported");
    let p = Pair::new([0, 1], t[0], t[0]);
    let (eps, sigma) = (p.well, p.distance / 2f64.powf(1.0 / 6.0));
    let n = 54.0;
    let rho = n / cell.volume();
    let pi = std::f64::consts::PI;
    let s3 = (sigma / rc).powi(3);
    let e_at = 8.0 / 3.0 * pi * n * rho * eps * sigma.powi(3) * (s3 * s3 * s3 / 3.0 - s3);
    let p_at = 16.0 / 3.0 * pi * rho * rho * eps * sigma.powi(3) * (2.0 / 3.0 * s3 * s3 * s3 - s3);
    println!(
        "one type: E_tail {:.10e} against A&T's {e_at:.10e}; P_tail {:.10e} against {p_at:.10e}",
        f.dispersion_correction(),
        f.pressure_correction()
    );
    assert!((f.dispersion_correction() / e_at - 1.0).abs() < 1e-13);
    assert!((f.pressure_correction() / p_at - 1.0).abs() < 1e-13);
}

/// Langevin dynamics of the box, run in one call and cut into pieces: the same bits.
#[test]
fn dynamics_in_the_box_is_the_same_bits_however_it_is_cut() {
    let (field, at0) = water_box(3, 9.3, 4.6, 1e-5, 0xD1CE);
    let c = waters(&water_sites(3, 9.3, 0xD1CE));
    let md = MolecularDynamics::for_elements(c.atoms().iter().map(|a| a.element))
        .with_bath(Bath::Langevin {
            temperature: 300.0,
            friction: 1e12,
            seed: 17,
        })
        .thermalised(&at0, 300.0, 5);
    let dt = 0.5e-15;
    let mut whole = (md.clone(), at0.clone());
    whole.0.run(&field, &mut whole.1, dt, 24);
    for cuts in [vec![1usize; 24], vec![5, 19], vec![11, 1, 12]] {
        let (mut m, mut at) = (md.clone(), at0.clone());
        for k in cuts {
            m.run(&field, &mut at, dt, k);
        }
        assert_eq!(at, whole.1);
        assert_eq!(m.velocities(), whole.0.velocities());
    }
    // And the atoms moved, so the comparison was of something.
    let moved = whole
        .1
        .iter()
        .zip(&at0)
        .map(|(a, b)| common::len(common::sub(*a, *b)))
        .fold(0.0f64, f64::max);
    assert!(moved > 1e-3 * ANGSTROM, "moved {moved:e}");
    // An evaluation is its own bits.
    let mut f1 = vec![[0.0; 3]; at0.len()];
    let mut f2 = f1.clone();
    let e1 = field.energy_and_forces(&at0, &mut f1);
    let e2 = field.energy_and_forces(&at0, &mut f2);
    assert_eq!((e1.to_bits(), &f1), (e2.to_bits(), &f2));
}

/// The decoupling's systems: two waters decoupled from 25 others — two, so that the group has
/// non-bonded pairs of its own — the whole box's force field, the 25 alone in the same box, and
/// the two alone in vacuum.
struct Decoupled {
    decoupling: PeriodicDecoupling,
    field: PeriodicForceField,
    rest: PeriodicForceField,
    alone_in_the_box: PeriodicForceField,
    vacuum: ForceField,
    at: Vec<[f64; 3]>,
}

fn decoupled() -> Decoupled {
    let sites = water_sites(3, 9.3, 0xDEC0);
    let cell = PeriodicBox::cubic(9.3 * ANGSTROM);
    let p = EwaldParameters::for_accuracy(&cell, 4.6 * ANGSTROM, 1e-8);
    // The decoupled pair of waters carries +0.1 e and the next water −0.1 e, so that the group
    // and the rest are each charged and the background has a cross term; the box stays neutral.
    let build = |s: &[[[f64; 3]; 3]], tweak: &[(usize, f64)]| {
        let c = waters(s);
        let mut q = tip3p_charges(&c);
        for &(i, dq) in tweak {
            q[i] += dq;
        }
        let f = PeriodicForceField::new(&c, &uff::assign(&c), cell, p)
            .expect("water")
            .with_charges(q.clone());
        (c, f, q)
    };
    let (c, field, _) = build(&sites, &[(0, 0.1), (6, -0.1)]);
    let (_, rest, _) = build(&sites[2..], &[(0, -0.1)]);
    let (one, alone_in_the_box, q_one) = build(&sites[..2], &[(0, 0.1)]);
    let vacuum = ForceField::new(&one, &uff::assign(&one))
        .expect("water")
        .with_charges(q_one);
    let mut group = vec![false; c.atoms().len()];
    group[..6].fill(true);
    Decoupled {
        decoupling: PeriodicDecoupling::new(&field, &group).expect("two waters decouple"),
        field,
        rest,
        alone_in_the_box,
        vacuum,
        at: positions(&c),
    }
}

/// `∂U/∂λ_e` is the cross Coulomb energy of the Ewald sum — `E(all) − E(rest) − E(group)`, three
/// sums of the box with charges taken away — at every λ.
#[test]
fn the_electrostatic_gradient_is_the_ewald_sums_cross_terms() {
    let d = decoupled();
    let q = d.field.charges().to_vec();
    let mut q_rest = q.clone();
    q_rest[..6].fill(0.0);
    let mut q_group = vec![0.0; q.len()];
    q_group[..6].copy_from_slice(&q[..6]);
    let ewald = |charges: &[f64]| {
        d.field
            .clone()
            .with_charges(charges.to_vec())
            .energy(&d.at)
            .ewald
    };
    let (all, rest, group) = (ewald(&q), ewald(&q_rest), ewald(&q_group));
    let cross = all.total - rest.total - group.total;
    let scale = all.real.abs() + all.reciprocal.abs() + all.self_energy.abs() + all.excluded.abs();
    for lambda in [
        Lambda::COUPLED,
        Lambda::new(0.0, 0.4, 0.7),
        Lambda::new(0.0, 0.0, 0.2),
    ] {
        let c = d.decoupling.coupling(&d.at, lambda);
        println!(
            "λ = {lambda:?}: ∂U/∂λ_e = {:.10e} J, the sums' cross terms {cross:.10e}",
            c.gradient[1]
        );
        assert!(
            (c.gradient[1] - cross).abs() <= 1e-12 * scale,
            "{} against {cross}",
            c.gradient[1]
        );
    }
    assert!(
        cross.abs() > 1e-3 * scale,
        "the cross terms are not negligible"
    );
}

/// Fully coupled, the decoupling is the box's force field less what the water has with its own
/// images, and with its own pairs in vacuum; fully decoupled, it is the other 25 in the box and
/// the two in vacuum.
#[test]
fn the_end_states_are_the_systems_they_say() {
    let d = decoupled();
    let e = |l: Lambda| {
        let mut f = vec![[0.0; 3]; d.at.len()];
        d.decoupling.energy_and_forces(&d.at, l, &mut f)
    };
    let coupled = e(Lambda::COUPLED);
    let whole = d.field.energy(&d.at);
    let in_box = d.alone_in_the_box.energy(&d.at[..6]);
    let vac = d.vacuum.energy(&d.at[..6]);
    let expected = whole.total - in_box.total + vac.total;
    let scale = magnitude(&whole);
    println!(
        "coupled: {:.10e} J against {expected:.10e}; the waters' own images were {:.3e} kcal/mol",
        coupled,
        (in_box.total - vac.total) / KCAL_PER_MOL
    );
    assert!((coupled - expected).abs() <= 1e-12 * scale);
    let decoupled = e(Lambda::new(0.0, 0.0, 0.0));
    assert!(
        !d.decoupling.own_pairs().is_empty(),
        "the group has pairs of its own"
    );
    let expected = d.rest.energy(&d.at[6..]).total + vac.total;
    println!("decoupled: {decoupled:.10e} J against {expected:.10e}");
    assert!((decoupled - expected).abs() <= 1e-12 * scale);
}

/// The λ gradient against central differences in λ, and the forces at an intermediate state
/// against central differences in the positions.
#[test]
fn the_decoupled_derivatives_are_the_gradient() {
    let d = decoupled();
    let at = &d.at;
    let lambda = Lambda::new(0.0, 0.45, 0.55);
    let c = d.decoupling.coupling(at, lambda);
    let scale = magnitude(&d.field.energy(at));
    // The coupling is a sum over the group's pairs with the rest and the reciprocal cross terms,
    // a few thousand terms none larger than the box's own energy parts: its rounding is held at
    // 64 ε of those, `scale`, as the positions' is.
    for (k, h) in [(1usize, 1e-3), (2, 1e-4)] {
        let energy = |s: f64| {
            let mut l = lambda;
            if k == 1 {
                l.electrostatics += s;
            } else {
                l.van_der_waals += s;
            }
            (d.decoupling.coupling(at, l).energy, scale)
        };
        let (dl, tolerance) = finite_difference(energy, h);
        println!(
            "∂U/∂λ[{k}] = {:.10e}, by differences {dl:.10e} ± {tolerance:.1e}",
            c.gradient[k]
        );
        assert!((c.gradient[k] - dl).abs() <= tolerance);
        assert!(
            tolerance < 1e-3 * c.gradient[k].abs(),
            "the check can see a part in 10³"
        );
    }
    let mut forces = vec![[0.0; 3]; at.len()];
    d.decoupling.energy_and_forces(at, lambda, &mut forces);
    let h = 1e-5 * ANGSTROM;
    for i in [0, 2, 4, 6, 40] {
        for a in 0..3 {
            let energy = |s: f64| {
                let mut moved = at.clone();
                moved[i][a] += s;
                let mut f = vec![[0.0; 3]; at.len()];
                (
                    d.decoupling.energy_and_forces(&moved, lambda, &mut f),
                    scale,
                )
            };
            let (dx, tolerance) = finite_difference(energy, h);
            assert!(
                (forces[i][a] + dx).abs() <= tolerance,
                "atom {i} axis {a}: {:e} against {:e} ± {tolerance:e}",
                forces[i][a],
                -dx
            );
            assert!(tolerance < 1e-5 * forces[i][a].abs());
        }
    }
    // The soft core changes the pairs across, so it is in: at λ_v = 1 the cross pairs are UFF's.
    let soft = d
        .decoupling
        .clone()
        .with_soft_core(SoftCore {
            alpha: 0.9,
            power: 1.0,
        })
        .coupling(at, lambda);
    assert_ne!(soft.energy, c.energy);
}

#[test]
fn a_decoupling_refuses_what_it_cannot_decouple() {
    let d = decoupled();
    let n = d.at.len();
    let mut split_water = vec![false; n];
    split_water[0] = true;
    assert!(matches!(
        PeriodicDecoupling::new(&d.field, &split_water),
        Err(AlchemyError::BondedAcross { .. })
    ));
    assert_eq!(
        PeriodicDecoupling::new(&d.field, &vec![false; n]),
        Err(AlchemyError::NothingToDecouple)
    );
    assert_eq!(
        PeriodicDecoupling::new(&d.field, &vec![true; n]),
        Err(AlchemyError::NothingToDecouple)
    );
}

/// The cost of the Ewald sum against the number of atoms, at water's density, and of the whole
/// force field on a box of a thousand waters. Release, one core: `cargo test --release -p
/// pantometry-forcefield --test a_molecule_in_a_periodic_box the_cost -- --ignored --nocapture`.
#[test]
#[ignore = "a timing, release only"]
fn the_cost() {
    use pantometry_forcefield::Ewald;
    let best = |f: &mut dyn FnMut()| {
        let mut t = f64::INFINITY;
        for _ in 0..3 {
            let s = std::time::Instant::now();
            f();
            t = t.min(s.elapsed().as_secs_f64());
        }
        t
    };
    for accuracy in [1e-5, 1e-6] {
        for k in [5usize, 7, 10, 13, 16, 20] {
            let n_w = k * k * k;
            // 0.03343 waters per Å³: 997 kg m⁻³ at 18.015 g/mol.
            let side = (n_w as f64 / 0.033_43).cbrt();
            let c = waters(&water_sites(k, side, 0xC057));
            let q = tip3p_charges(&c);
            let at = positions(&c);
            let cell = PeriodicBox::cubic(side * ANGSTROM);
            let cutoff = (9.0f64).min(0.5 * side) * ANGSTROM;
            let p = EwaldParameters::for_accuracy(&cell, cutoff, accuracy);
            let excluded: Vec<[usize; 2]> = (0..n_w)
                .flat_map(|m| {
                    [
                        [3 * m, 3 * m + 1],
                        [3 * m, 3 * m + 2],
                        [3 * m + 1, 3 * m + 2],
                    ]
                })
                .collect();
            let ewald = Ewald::new(cell, p);
            let real_only = Ewald::new(
                cell,
                EwaldParameters {
                    k_cutoff: 1e-3 / ANGSTROM,
                    ..p
                },
            );
            assert_eq!(real_only.wave_vectors(), 0);
            let total = best(&mut || {
                std::hint::black_box(ewald.evaluate(&q, &at, &excluded));
            });
            let real = best(&mut || {
                std::hint::black_box(real_only.evaluate(&q, &at, &excluded));
            });
            println!(
                "δ = {accuracy:e}: {:5} atoms, {side:5.1} Å, r_c {:.1} Å, α {:.3} /Å, {:5} wave \
                 vectors: Ewald {:8.1} ms (real space and exclusions {:6.1}, reciprocal {:8.1})",
                3 * n_w,
                cutoff / ANGSTROM,
                p.alpha * ANGSTROM,
                ewald.wave_vectors(),
                1e3 * total,
                1e3 * real,
                1e3 * (total - real)
            );
            if k == 10 {
                let field = PeriodicForceField::new(&c, &uff::assign(&c), cell, p)
                    .expect("water")
                    .with_charges(q.clone());
                let whole = best(&mut || {
                    std::hint::black_box(field.evaluate(&at));
                });
                println!(
                    "    the whole force field on the same box: {:.1} ms",
                    1e3 * whole
                );
            }
        }
    }
}
