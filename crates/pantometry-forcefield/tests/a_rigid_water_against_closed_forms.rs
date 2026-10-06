//! **Rigid TIP3P water against closed forms**: the model is Jorgensen's Table I, SETTLE is the
//! constrained step SHAKE converges to and restores the geometry to rounding, the projected
//! velocities have nothing along a constraint, the energy error of constrained NVE falls as `h²`,
//! a bath shares its energy as `(3/2) k_BT` of translation and `(3/2) k_BT` of rotation per water,
//! the box's forces are its energy's gradient, and a run is the same bits however it is cut or
//! frozen. Every test here runs unoptimised; the liquid itself is
//! `tests/liquid_water_against_the_literature.rs`, in release.

// The sums over atoms and axes are written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

mod common;

use pantometry_forcefield::ewald::COULOMB;
use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::water::{self, Settle, WaterBox};
use pantometry_forcefield::{
    Bath, Estimate, Ewald, EwaldParameters, MolecularDynamics, PeriodicForceField, Potential,
};
use pantometry_units::BOLTZMANN;

const ANGSTROM: f64 = 1e-10;
const EPS: f64 = f64::EPSILON;
/// One debye in elementary charges × metres: 10⁻²¹/c C m over e.
const DEBYE: f64 = 1e-21 / 299_792_458.0 / 1.602_176_634e-19;

/// `x` to the whole power `n` by repeated multiplication: `powi`'s precision is unspecified and
/// differs between platforms, so no test input or reference is built with it.
fn pow(x: f64, n: u32) -> f64 {
    (0..n).fold(1.0, |p, _| p * x)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    common::sub(a, b)
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    common::dot(a, b)
}
fn len(a: [f64; 3]) -> f64 {
    common::len(a)
}

/// A uniform deviate in [0, 1) keyed by `(seed, i)`: a splitmix64 finaliser, as the W1 tests use.
fn uniform(seed: u64, i: u64) -> f64 {
    let mut z = seed.wrapping_add((i + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    (z >> 11) as f64 / (1u64 << 53) as f64
}

// ---------------------------------------------------------------------------------------------
// The model

/// Table I as it stands, and what follows from it: ε, σ, the charges' sum, the geometry of the
/// rigid triangle and the dipole.
#[test]
fn the_model_is_table_one() {
    assert_eq!(water::OXYGEN_CHARGE, -0.834);
    assert_eq!(water::HYDROGEN_CHARGE, 0.417);
    assert_eq!(water::OH_LENGTH, 0.9572e-10);
    assert_eq!(water::HOH_ANGLE_DEGREES, 104.52);
    assert_eq!((water::LJ_A, water::LJ_C), (582.0e3, 595.0));
    // The total charge is zero exactly, not to rounding.
    assert_eq!(water::CHARGES.iter().sum::<f64>(), 0.0);
    assert_eq!(
        water::OXYGEN_CHARGE + water::HYDROGEN_CHARGE + water::HYDROGEN_CHARGE,
        0.0
    );
    // ε and σ: A = 4εσ¹², C = 4εσ⁶.
    let (e, s) = (water::epsilon() / KCAL_PER_MOL, water::sigma() / ANGSTROM);
    println!("ε = {e:.6} kcal/mol, σ = {s:.6} Å");
    assert!((4.0 * e * pow(s, 12) / water::LJ_A - 1.0).abs() < 1e-14);
    assert!((4.0 * e * pow(s, 6) / water::LJ_C - 1.0).abs() < 1e-14);
    // Those two are the derivation, true of any A and C; this is the transcription, against a
    // third source: LAMMPS's TIP3P page prints ε = 0.1521 kcal/mol and σ = 3.1507 Å, which A and C
    // must give to the printed digits — half a unit in the last one.
    assert!((e - 0.1521).abs() <= 0.5e-4, "ε = {e}");
    assert!((s - 3.1507).abs() <= 0.5e-4, "σ = {s}");
    // The oxygen as a UFF atom: D [(x/r)¹² − 2 (x/r)⁶] is A/r¹² − C/r⁶ at every r, to 45 ulps
    // of the larger of its two terms — at 9 Å that is C/r⁶, not A/r¹².
    let [x, d] = water::oxygen_vdw();
    for r in [2.6, 2.8, 3.1, 3.5, 4.5, 9.0] {
        let u = pow(x / ANGSTROM / r, 6);
        let uff = d / KCAL_PER_MOL * (u * u - 2.0 * u);
        let table = water::LJ_A / pow(r, 12) - water::LJ_C / pow(r, 6);
        assert!(
            (uff - table).abs() < 1e-14 * (water::LJ_A / pow(r, 12) + water::LJ_C / pow(r, 6)),
            "r = {r}"
        );
    }
    // The rigid triangle: its three lengths, its angle and its centre of mass.
    let [o, h1, h2] = water::body_frame();
    let m = water::masses();
    let rel = |a: f64, b: f64| (a / b - 1.0).abs();
    assert!(rel(len(sub(o, h1)), water::OH_LENGTH) < 4.0 * EPS);
    assert!(rel(len(sub(o, h2)), water::OH_LENGTH) < 4.0 * EPS);
    assert!(rel(len(sub(h1, h2)), water::hh_length()) < 4.0 * EPS);
    let angle = (dot(sub(h1, o), sub(h2, o)) / pow(water::OH_LENGTH, 2)).acos();
    assert!((angle - water::hoh_angle()).abs() < 1e-14);
    for k in 0..3 {
        let c = m[0] * o[k] + m[1] * h1[k] + m[2] * h2[k];
        assert!(c.abs() < 4.0 * EPS * m[0] * water::OH_LENGTH);
    }
    // The dipole Σ q r against 2 q_H r_OH cos(θ/2), and TIP3P's 2.35 D.
    let mu: Vec<f64> = (0..3)
        .map(|k| water::CHARGES[0] * o[k] + water::CHARGES[1] * h1[k] + water::CHARGES[2] * h2[k])
        .collect();
    let mu = (mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2]).sqrt();
    println!(
        "μ = {:.6} e Å = {:.5} D; closed form {:.6} e Å",
        mu / ANGSTROM,
        mu / DEBYE,
        water::dipole_moment() / ANGSTROM
    );
    assert!(rel(mu, water::dipole_moment()) < 1e-14);
    assert!((water::dipole_moment() / DEBYE - 2.35).abs() < 0.005);
}

/// The box's force field is TIP3P written out: the O–O pair `A/r¹² − C/r⁶`, nothing on a
/// hydrogen, the Ewald sum of the waters' charges with each water's three pairs excluded — and a
/// UFF solute's pair with a water oxygen UFF's geometric rule.
#[test]
fn the_force_field_is_tip3p_written_out() {
    let wb = WaterBox::lattice(2, 0x11);
    let at = wb.positions().to_vec();
    let cell = wb.cell();
    let rc = 0.5 * cell.lengths()[0];
    let p = EwaldParameters::for_accuracy(&cell, rc, 1e-8);
    let field = PeriodicForceField::tip3p(cell, p, wb.count());
    assert_eq!(field.rigid_waters(), &wb.waters()[..]);
    assert_eq!(
        field.charges(),
        &(0..8).flat_map(|_| water::CHARGES).collect::<Vec<_>>()[..]
    );
    let e = field.energy(&at);
    // Van der Waals: every O–O pair inside the cutoff by minimum image.
    let mut lj = 0.0;
    for i in 0..8 {
        for j in i + 1..8 {
            let r = len(cell.minimum_image(sub(at[3 * i], at[3 * j]))) / ANGSTROM;
            if r < rc / ANGSTROM {
                lj += water::LJ_A / pow(r, 12) - water::LJ_C / pow(r, 6);
            }
        }
    }
    let got = e.van_der_waals / KCAL_PER_MOL;
    println!("van der Waals {got:.12} against {lj:.12} kcal/mol");
    assert!((got - lj).abs() < 1e-13 * lj.abs().max(1.0));
    // The long-range correction for one type, A&T's: (8/3)πNρεσ³[(σ/r_c)⁹/3 − (σ/r_c)³].
    let (eps, sig) = (water::epsilon(), water::sigma());
    let rho = 8.0 / cell.volume();
    let q = sig / rc;
    let tail = 8.0 / 3.0
        * std::f64::consts::PI
        * 8.0
        * rho
        * eps
        * pow(sig, 3)
        * (pow(q, 9) / 3.0 - pow(q, 3));
    assert!((e.dispersion_correction / tail - 1.0).abs() < 1e-13);
    // Electrostatics: the Ewald sum with the intramolecular pairs excluded.
    let excluded: Vec<[usize; 2]> = wb
        .waters()
        .iter()
        .flat_map(|&[o, a, b]| [[o, a], [o, b], [a, b]])
        .collect();
    let sum = Ewald::new(cell, p).evaluate(field.charges(), &at, &excluded);
    assert!(
        (e.electrostatic - sum.energy.total).abs()
            < 1e-13
                * (sum.energy.real.abs()
                    + sum.energy.reciprocal.abs()
                    + sum.energy.self_energy.abs())
    );
    // No bonded term.
    assert_eq!([e.bond, e.angle, e.torsion, e.inversion], [0.0; 4]);
    assert!(field.bonded().is_none());

    // A solute: methane among the waters, its pairs with each water oxygen by UFF's rule.
    let methane = common::entry(
        &[
            ("C", "C", [0.0, 0.0, 0.0]),
            ("H1", "H", [0.629, 0.629, 0.629]),
            ("H2", "H", [-0.629, -0.629, 0.629]),
            ("H3", "H", [-0.629, 0.629, -0.629]),
            ("H4", "H", [0.629, -0.629, -0.629]),
        ],
        &[
            ("C", "H1", "SING", false),
            ("C", "H2", "SING", false),
            ("C", "H3", "SING", false),
            ("C", "H4", "SING", false),
        ],
    );
    let types = uff::assign(&methane);
    let solvated = PeriodicForceField::solvated(&methane, &types, 8, cell, p).expect("methane");
    assert_eq!(solvated.rigid_waters()[0], [5, 6, 7]);
    let mut all: Vec<[f64; 3]> = methane
        .atoms()
        .iter()
        .map(|a| {
            let c = cell.lengths()[0] * 0.5;
            [a.at[0] + c, a.at[1] + c + 0.37 * ANGSTROM, a.at[2] + c]
        })
        .collect();
    all.extend(&at);
    let charges_off = solvated.clone().with_charges(vec![0.0; all.len()]);
    let alone = PeriodicForceField::tip3p(cell, p, 8).with_charges(vec![0.0; 24]);
    // Methane's own pairs are all 1-2 or 1-3, so what the solute adds is the pairs across.
    let across = charges_off.energy(&all).van_der_waals - alone.energy(&at).van_der_waals;
    let [xo, dox] = water::oxygen_vdw();
    let mut want = 0.0;
    for (i, t) in types.iter().enumerate() {
        let pi = t.parameters();
        for w in 0..8 {
            let r = len(cell.minimum_image(sub(all[i], at[3 * w])));
            if r < rc {
                let x = (pi.vdw_distance * xo).sqrt();
                let d = (pi.vdw_energy * dox).sqrt();
                let u = pow(x / r, 6);
                want += d * (u * u - 2.0 * u);
            }
        }
    }
    assert!(
        (across - want).abs() < 1e-12 * want.abs().max(1e-25),
        "{across:e} {want:e}"
    );
}

// ---------------------------------------------------------------------------------------------
// SETTLE

/// How far one water at `p` is from rigid, in units of the rounding of its coordinates: the
/// largest `|d − d₀|` of its three distances over `ε` times its largest coordinate. A position is
/// stored to half an ulp of itself, so a rigid water a few ångström from the origin cannot hold
/// its 0.96 Å bonds to an ulp of the bond, and a bound in ulps of the bond would be one in where
/// the origin is.
fn rigidity(p: [[f64; 3]; 3]) -> f64 {
    let l = Settle::tip3p(Vec::new()).lengths();
    let scale = p.iter().flatten().fold(0.0f64, |m, x| m.max(x.abs())) * EPS;
    let got = [
        len(sub(p[0], p[1])),
        len(sub(p[0], p[2])),
        len(sub(p[1], p[2])),
    ];
    (0..3)
        .map(|k| (got[k] - l[k]).abs() / scale)
        .fold(0.0, f64::max)
}

/// The bound on [`rigidity`]: SETTLE writes each coordinate as the centre of mass plus three
/// products, a few roundings of the coordinate's size each, and a distance takes the difference
/// of two such points in three components. 16 is that count with a margin; the tests print the
/// measured worst.
const RIGID: f64 = 16.0;

/// SHAKE, written here as Ryckaert, Ciccotti and Berendsen's iteration: each constraint in turn
/// corrected along its bond at the old positions, swept until no position changes. It solves the
/// equations SETTLE solves in closed form.
fn shake(old: &[[f64; 3]; 3], new: &[[f64; 3]; 3], m: [f64; 3], d: [f64; 3]) -> [[f64; 3]; 3] {
    let mut q = *new;
    let pairs = [(0usize, 1usize, d[0]), (0, 2, d[1]), (1, 2, d[2])];
    for _ in 0..1000 {
        let before = q;
        for &(a, b, l) in &pairs {
            let g = sub(q[a], q[b]);
            let r0 = sub(old[a], old[b]);
            let diff = dot(g, g) - l * l;
            let lambda = diff / (2.0 * (1.0 / m[a] + 1.0 / m[b]) * dot(g, r0));
            for k in 0..3 {
                q[a][k] -= lambda / m[a] * r0[k];
                q[b][k] += lambda / m[b] * r0[k];
            }
        }
        if q == before {
            break;
        }
    }
    q
}

/// Rigid waters, each given an unconstrained step of a thermal velocity at 2 fs plus a random
/// displacement up to 0.05 Å per atom: what one step of dynamics does to them.
fn perturbed(seed: u64, n: usize) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let wb = WaterBox::lattice(2, seed);
    let old: Vec<[f64; 3]> = wb.positions().iter().cycle().take(3 * n).copied().collect();
    let new: Vec<[f64; 3]> = old
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let kick = (0..3).map(|k| {
                let u = uniform(seed, 3 * i as u64 + k as u64) - 0.5;
                // ~1 km/s × 2 fs is 0.02 Å; the extra 0.05 Å is a force's.
                0.14 * u * ANGSTROM
            });
            let v: Vec<f64> = kick.collect();
            [p[0] + v[0], p[1] + v[1], p[2] + v[2]]
        })
        .collect();
    (old, new)
}

/// How near SHAKE's fixed point SETTLE's answer is, in ε × the largest coordinate: two different
/// computations of one solution, each a few dozen roundings, through a 3×3 system that is not
/// perfectly conditioned. A SETTLE that solved a different problem — a wrong root, a sign — would
/// differ by the size of the constraint correction itself, which the test also measures and
/// requires to be a million times larger.
const SHAKE_AGREES: f64 = 64.0;

/// After SETTLE every water has its three lengths to rounding ([`RIGID`]), its centre of mass
/// where the unconstrained step put it, and it is the position SHAKE converges to
/// ([`SHAKE_AGREES`]).
#[test]
fn settle_restores_the_geometry_and_is_shakes_answer() {
    let s = Settle::tip3p((0..64).map(|w| [3 * w, 3 * w + 1, 3 * w + 2]).collect());
    let lengths = s.lengths();
    let m = water::masses();
    let (mut worst_length, mut worst_shake) = (0.0f64, 0.0f64);
    let mut smallest_correction = f64::INFINITY;
    for seed in 0..8u64 {
        let (old, new) = perturbed(0x5E77 + seed, 64);
        let mut settled = new.clone();
        s.constrain_positions(&old, &mut settled);
        for w in 0..64 {
            let idx = [3 * w, 3 * w + 1, 3 * w + 2];
            let r = idx.map(|i| settled[i]);
            worst_length = worst_length.max(rigidity(r));
            for k in 0..3 {
                let c = |p: [[f64; 3]; 3]| m[0] * p[0][k] + m[1] * p[1][k] + m[2] * p[2][k];
                let (a, b) = (c(r), c(idx.map(|i| new[i])));
                // The rounding of a mass-weighted sum of positions: a few ulps of its terms.
                let terms: f64 = (0..3).map(|i| m[i] * r[i][k].abs()).sum();
                assert!((a - b).abs() < 16.0 * EPS * terms, "{:e}", (a - b) / terms);
            }
            let q = shake(&idx.map(|i| old[i]), &idx.map(|i| new[i]), m, lengths);
            let unit = r.iter().flatten().fold(0.0f64, |a, x| a.max(x.abs())) * EPS;
            for i in 0..3 {
                worst_shake = worst_shake.max(len(sub(q[i], r[i])) / unit);
                smallest_correction = smallest_correction.min(len(sub(new[idx[i]], r[i])) / unit);
            }
        }
    }
    println!(
        "SETTLE: worst length error {worst_length:.2}, worst distance from SHAKE \
         {worst_shake:.2}, smallest correction {smallest_correction:.3e}, each in ε × the \
         largest coordinate"
    );
    assert!(worst_length <= RIGID);
    assert!(worst_shake <= SHAKE_AGREES);
    // And the comparison was of something: every atom was moved by the constraint far more
    // than the agreement allows.
    assert!(smallest_correction > 1e6 * SHAKE_AGREES);
}

/// The projection leaves no velocity along a constraint, is idempotent, and does not touch a
/// rigid motion of a whole water.
#[test]
fn constrained_velocities_have_nothing_along_a_constraint() {
    let wb = WaterBox::lattice(3, 0xF00);
    let at = wb.positions().to_vec();
    let s = Settle::tip3p(wb.waters());
    let mut v: Vec<[f64; 3]> = (0..at.len())
        .map(|i| [0, 1, 2].map(|k| 1e3 * (uniform(9, 3 * i as u64 + k) - 0.5)))
        .collect();
    s.constrain_velocities(&at, &mut v);
    let mut worst = 0.0f64;
    for &[o, a, b] in &wb.waters() {
        for (i, j) in [(o, a), (o, b), (a, b)] {
            let r = sub(at[i], at[j]);
            let u = sub(v[i], v[j]);
            worst = worst.max(dot(r, u).abs() / (len(r) * (len(v[i]) + len(v[j]))));
        }
    }
    println!("worst |r·v| / (|r| |v|): {worst:.2e}");
    assert!(worst < 1e-14);
    let again = {
        let mut w = v.clone();
        s.constrain_velocities(&at, &mut w);
        w
    };
    for (x, y) in again.iter().zip(&v) {
        assert!(len(sub(*x, *y)) < 1e-13 * 1e3);
    }
    // A translation plus a rotation about each water's oxygen is rigid, and kept.
    let rigid: Vec<[f64; 3]> = at
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let o = at[3 * (i / 3)];
            let w = [3e12, -1e12, 2e12];
            let c = common::cross(w, sub(*p, o));
            [c[0] + 100.0, c[1] - 50.0, c[2] + 20.0]
        })
        .collect();
    let mut kept = rigid.clone();
    s.constrain_velocities(&at, &mut kept);
    for (x, y) in kept.iter().zip(&rigid) {
        assert!(len(sub(*x, *y)) < 1e-12 * len(*y));
    }
}

// ---------------------------------------------------------------------------------------------
// Dynamics

/// TIP3P waters in vacuum, every intermolecular pair, no cutoff and nothing periodic: a smooth
/// potential, so that the energy error of NVE is the integrator's and not a cutoff's steps.
/// Written out here from Table I, independently of the box's force field.
pub struct Cluster;

impl Potential for Cluster {
    fn energy_and_forces(&self, at: &[[f64; 3]], f: &mut [[f64; 3]]) -> f64 {
        let a = water::LJ_A * KCAL_PER_MOL * pow(ANGSTROM, 12);
        let c = water::LJ_C * KCAL_PER_MOL * pow(ANGSTROM, 6);
        let n = at.len() / 3;
        f.iter_mut().for_each(|x| *x = [0.0; 3]);
        let mut e = 0.0;
        for wi in 0..n {
            for wj in wi + 1..n {
                for p in 0..3 {
                    for q in 0..3 {
                        let (i, j) = (3 * wi + p, 3 * wj + q);
                        let d = sub(at[i], at[j]);
                        let r = len(d);
                        let qq = water::CHARGES[p] * water::CHARGES[q];
                        e += COULOMB * qq / r;
                        let mut de = -COULOMB * qq / (r * r);
                        if p == 0 && q == 0 {
                            e += a / pow(r, 12) - c / pow(r, 6);
                            de += -12.0 * a / pow(r, 13) + 6.0 * c / pow(r, 7);
                        }
                        let g = -de / r;
                        for k in 0..3 {
                            f[i][k] += g * d[k];
                            f[j][k] -= g * d[k];
                        }
                    }
                }
            }
        }
        e
    }
}

/// The largest and the RMS departure of `kinetic + potential` from its start over `steps` of
/// `dt` in NVE, from start number `start`.
pub fn nve_departure(dt: f64, steps: usize, start: u64) -> [f64; 2] {
    let wb = WaterBox::lattice(2, 0xE4E + start);
    let mut at = wb.positions().to_vec();
    let mut md = wb.dynamics().thermalised(&at, 300.0, 3 + start);
    // Six a water, less the rigid motion `thermalised` removed: NVE keeps it removed.
    assert_eq!(md.degrees_of_freedom(), 6 * wb.count() - 6);
    md.prepare(&Cluster, &at);
    let e0 = md.kinetic_energy() + md.potential_energy().expect("prepared");
    let (mut worst, mut squares) = (0.0f64, 0.0);
    for _ in 0..steps {
        md.step(&Cluster, &mut at, dt);
        let e = md.kinetic_energy() + md.potential_energy().expect("stepped");
        worst = worst.max((e - e0).abs());
        squares += (e - e0) * (e - e0);
    }
    [worst, (squares / steps as f64).sqrt()]
}

/// The band the RMS energy error's ratio per halving of the step is held to. Measured over sixteen
/// starts (eight waters, 300 K, 0.4 ps) at 1, 0.5 and 0.25 fs, the 32 ratios were 3.911 to 4.057;
/// the band is twice that half-width about 4. An error of order `h` or `h³` would read 2 or 8.
/// From 2 fs to 1 fs the ratios were 3.36 to 4.09: one start of the sixteen reads 3.36 there,
/// because 2 fs is past where `h²` alone describes it, which is why the test halves from 1 fs.
const H_SQUARED: std::ops::Range<f64> = 3.8..4.2;

/// The bound on the books' largest departure under a bath, joules: ten times constrained NVE's
/// RMS energy error from the same start at the same step, 5.06e-3 kcal/mol, as
/// `constrained_nve_energy_error_falls_as_h_squared` measures it — the integrator's own error,
/// which is all the books may show. **A constant, not that measurement taken again here**: a
/// sabotaged integrator has a larger NVE error too, and a bound that grew with it let a missing
/// `Δq/h` through at 0.2 of itself. Measured, the books moved 5.1e-3; the projection after `O` at
/// the positions before the drift moved them 0.67, the `Δq/h` correction left out 1.36 and
/// halved 0.68.
const BOOKS: f64 = 10.0 * 5.06e-3 * KCAL_PER_MOL;

/// Constrained velocity Verlet is symplectic: the energy's error stays in a band of `O(h²)`. Over
/// the same 0.4 ps the RMS departure falls by four per halving of the step, held to
/// [`H_SQUARED`], on one of the sixteen starts the band was measured on.
#[test]
fn constrained_nve_energy_error_falls_as_h_squared() {
    let fs = 1e-15;
    let e1 = nve_departure(1.0 * fs, 400, 0)[1];
    let e05 = nve_departure(0.5 * fs, 800, 0)[1];
    let e025 = nve_departure(0.25 * fs, 1600, 0)[1];
    let kcal = |x: f64| x / KCAL_PER_MOL;
    println!(
        "RMS ΔE over 0.4 ps: 1 fs {:.3e}, 0.5 fs {:.3e}, 0.25 fs {:.3e} kcal/mol; ratios {:.4}, {:.4}",
        kcal(e1),
        kcal(e05),
        kcal(e025),
        e1 / e05,
        e05 / e025
    );
    for ratio in [e1 / e05, e05 / e025] {
        assert!(H_SQUARED.contains(&ratio), "ratio {ratio}");
    }
}

/// The books balance with constraints under a bath: `kinetic + potential − thermostat work`
/// moves only by the integrator's error, which at 0.5 fs is the NVE band's size, while the work
/// itself — what `O` and its projection put in and take out — moves by far more. Booking the work
/// before the projection would leave the energy the projection removes from the bonds in the
/// books every step.
#[test]
fn the_books_balance_under_a_bath_with_constraints() {
    let wb = WaterBox::lattice(2, 0xE4E);
    let mut at = wb.positions().to_vec();
    let mut md = wb
        .dynamics()
        .with_bath(Bath::Langevin {
            temperature: 300.0,
            friction: 20e12,
            seed: 5,
        })
        .thermalised(&at, 300.0, 3);
    md.prepare(&Cluster, &at);
    let books = |m: &MolecularDynamics| {
        m.kinetic_energy() + m.potential_energy().expect("prepared") - m.thermostat_work()
    };
    let start = books(&md);
    let (mut worst, mut traffic, mut last) = (0.0f64, 0.0, 0.0);
    for _ in 0..800 {
        md.step(&Cluster, &mut at, 0.5e-15);
        worst = worst.max((books(&md) - start).abs());
        traffic += (md.thermostat_work() - last).abs();
        last = md.thermostat_work();
    }
    println!(
        "books moved at most {:.3e} kcal/mol against a bound of {:.3e}; the bath moved {:.3e} \
         through them",
        worst / KCAL_PER_MOL,
        BOOKS / KCAL_PER_MOL,
        traffic / KCAL_PER_MOL
    );
    assert!(worst < BOOKS, "{:e} kcal/mol", worst / KCAL_PER_MOL);
}

/// The 27-water box under a bath, its force field and a thermalised start, run `steps` of 0.5 fs.
fn bathed(
    friction: f64,
    seed: u64,
    steps: usize,
) -> (
    PeriodicForceField,
    MolecularDynamics,
    Vec<[f64; 3]>,
    WaterBox,
) {
    let wb = WaterBox::lattice(3, 0xBA7);
    let field = wb.force_field(4.6 * ANGSTROM, 1e-5);
    let mut at = wb.positions().to_vec();
    let mut md = wb
        .dynamics()
        .with_bath(Bath::Langevin {
            temperature: 300.0,
            friction,
            seed,
        })
        .thermalised(&at, 300.0, seed);
    md.run(&field, &mut at, 0.5e-15, steps);
    (field, md, at, wb)
}

/// Equipartition for rigid bodies: in a bath, the waters' translation and their rotation each
/// carry `(3/2) k_BT` a water on average, measured after `O` where BAOAB's momenta are canonical
/// — six degrees of freedom a water, which `degrees_of_freedom` counts. Each within four standard
/// errors from its series' own autocorrelation time.
#[test]
fn a_bath_gives_each_water_three_halves_kt_of_translation_and_of_rotation() {
    // One picosecond of bath first. With 0.1 ps the lattice was still giving up its strain: over
    // 20 ps from there, in release, translation read 1.0143 ± 0.0068 and rotation 1.0137 ± 0.0070;
    // from 1 ps, four starts read 0.985–1.002 and 0.991–1.011, each within 2σ.
    let (field, mut md, mut at, wb) = bathed(20e12, 41, 2000);
    assert_eq!(md.degrees_of_freedom(), 6 * 27);
    let kt = BOLTZMANN.to_si() * 300.0;
    let m = water::masses();
    let total_mass = water::molecular_mass();
    let (mut trans, mut rot) = (Vec::new(), Vec::new());
    for _ in 0..1500 {
        md.step(&field, &mut at, 1e-15);
        let v = md.half_step_velocities();
        let (mut t, mut r) = (0.0, 0.0);
        for &[o, a, b] in &wb.waters() {
            let com =
                [0, 1, 2].map(|k| (m[0] * v[o][k] + m[1] * v[a][k] + m[2] * v[b][k]) / total_mass);
            let ke_t = 0.5 * total_mass * dot(com, com);
            let ke: f64 = [o, a, b]
                .iter()
                .zip(m)
                .map(|(&i, mi)| 0.5 * mi * dot(v[i], v[i]))
                .sum();
            t += ke_t;
            r += ke - ke_t;
        }
        trans.push(t / (1.5 * 27.0 * kt));
        rot.push(r / (1.5 * 27.0 * kt));
    }
    for (name, series) in [("translation", &trans), ("rotation", &rot)] {
        let fifths: Vec<f64> = series
            .chunks(series.len() / 5)
            .map(|c| c.iter().sum::<f64>() / c.len() as f64)
            .collect();
        println!("{name} by fifths: {fifths:.3?}");
        let e = Estimate::of(series);
        let z = (e.mean - 1.0) / e.error;
        println!(
            "{name}: {:.4} ± {:.4} of (3/2) N k_BT, z = {z:+.2}, τ = {:.1}",
            e.mean, e.error, e.tau
        );
        assert!(z.abs() < 4.0, "{name}");
    }
}

/// No potential at all: free rigid waters, which overlap and pass through each other.
struct Nothing;

impl Potential for Nothing {
    fn energy_and_forces(&self, _: &[[f64; 3]], f: &mut [[f64; 3]]) -> f64 {
        f.iter_mut().for_each(|x| *x = [0.0; 3]);
        0.0
    }
}

/// How far, in standard errors, the bath may put each half of a free rigid water's kinetic
/// energy from `(3/2) k_BT`: four.
const SHARE_Z: f64 = 4.0;

/// **Equipartition to a part in a few hundred**: free rigid waters under no potential, in a bath
/// with `γh = 1`, so that each `O` replaces two thirds of the velocity and successive samples are
/// nearly independent — 40 000 of them where the interacting box above has about 75. Translation
/// and rotation each carry `(3/2) k_BT` a water; a noise width 2% wide, 4% in the energy, is many
/// standard errors out. With no potential nothing else can set the temperature, so this is the
/// bath and the projection alone.
#[test]
fn free_rigid_waters_take_the_baths_temperature_exactly() {
    let wb = WaterBox::lattice(3, 0xF8EE);
    let mut at = wb.positions().to_vec();
    let dt = 1e-15;
    let mut md = wb
        .dynamics()
        .with_bath(Bath::Langevin {
            temperature: 300.0,
            friction: 1.0 / dt,
            seed: 11,
        })
        .thermalised(&at, 300.0, 11);
    let kt = BOLTZMANN.to_si() * 300.0;
    let m = water::masses();
    let total_mass = water::molecular_mass();
    let n = wb.count() as f64;
    let (mut trans, mut rot) = (Vec::new(), Vec::new());
    for _ in 0..40_000 {
        md.step(&Nothing, &mut at, dt);
        let v = md.half_step_velocities();
        let (mut t, mut r) = (0.0, 0.0);
        for &[o, a, b] in &wb.waters() {
            let com =
                [0, 1, 2].map(|k| (m[0] * v[o][k] + m[1] * v[a][k] + m[2] * v[b][k]) / total_mass);
            let ke_t = 0.5 * total_mass * dot(com, com);
            let ke: f64 = [o, a, b]
                .iter()
                .zip(m)
                .map(|(&i, mi)| 0.5 * mi * dot(v[i], v[i]))
                .sum();
            t += ke_t;
            r += ke - ke_t;
        }
        trans.push(t / (1.5 * n * kt));
        rot.push(r / (1.5 * n * kt));
    }
    for (name, series) in [("translation", &trans), ("rotation", &rot)] {
        let e = Estimate::of(series);
        let z = (e.mean - 1.0) / e.error;
        println!(
            "free {name}: {:.5} ± {:.5} of (3/2) N k_BT, z = {z:+.2}, τ = {:.2}",
            e.mean, e.error, e.tau
        );
        assert!(z.abs() < SHARE_Z, "{name}");
    }
}

/// A rigid-water run under a bath is the same bits in one call or cut into several, the waters
/// stay rigid, and the velocities a step leaves — the ones the kinetic energy and the books are
/// read from — have nothing along a constraint. That last is its own check: a step whose final
/// kick is not projected follows the same trajectory, since the next step's first kick projects
/// at the same positions, and only the velocities it reports are wrong.
#[test]
fn a_rigid_run_is_the_same_bits_however_it_is_cut() {
    let (field, md, at0, wb) = bathed(5e12, 7, 200);
    let mut whole = (md.clone(), at0.clone());
    whole.0.run(&field, &mut whole.1, 2e-15, 24);
    let mut half_worst = 0.0f64;
    for cuts in [vec![1usize; 24], vec![5, 19], vec![11, 1, 12]] {
        let (mut m, mut at) = (md.clone(), at0.clone());
        for k in cuts {
            m.run(&field, &mut at, 2e-15, k);
            // After `O`, at the positions it was taken at, nothing along a constraint.
            let (p, u) = (m.half_step_positions(), m.half_step_velocities());
            for &[o, a, b] in &wb.waters() {
                for (i, j) in [(o, a), (o, b), (a, b)] {
                    let (d, w) = (sub(p[i], p[j]), sub(u[i], u[j]));
                    half_worst =
                        half_worst.max(dot(d, w).abs() / (len(d) * (len(u[i]) + len(u[j]))));
                }
            }
        }
        assert_eq!(at, whole.1);
        assert_eq!(m.velocities(), whole.0.velocities());
    }
    let (r, v) = (&whole.1, whole.0.velocities());
    let mut worst = 0.0f64;
    for &[o, a, b] in &wb.waters() {
        assert!(rigidity([r[o], r[a], r[b]]) <= RIGID);
        for (i, j) in [(o, a), (o, b), (a, b)] {
            let (d, u) = (sub(r[i], r[j]), sub(v[i], v[j]));
            worst = worst.max(dot(d, u).abs() / (len(d) * (len(v[i]) + len(v[j]))));
        }
    }
    println!("worst |r·v| / (|r| |v|): whole-step {worst:.2e}, after O {half_worst:.2e}");
    assert!(worst < 1e-14);
    assert!(half_worst < 1e-14);
}

/// Whole frozen waters keep their bits while the rest move; a water frozen in part is refused.
#[test]
fn frozen_waters_keep_their_bits() {
    let wb = WaterBox::lattice(3, 0xF2);
    let field = wb.force_field(4.6 * ANGSTROM, 1e-5);
    let mut at = wb.positions().to_vec();
    let frozen: Vec<bool> = (0..at.len()).map(|i| i < 9).collect();
    let mut md = wb
        .dynamics()
        .with_frozen(frozen)
        .with_bath(Bath::Langevin {
            temperature: 300.0,
            friction: 5e12,
            seed: 3,
        })
        .thermalised(&at, 300.0, 3);
    assert_eq!(md.degrees_of_freedom(), 6 * 24);
    let start = at.clone();
    md.run(&field, &mut at, 1e-15, 20);
    assert_eq!(&at[..9], &start[..9]);
    assert!(at[9..].iter().zip(&start[9..]).all(|(a, b)| a != b));
}

/// SETTLE moves a whole triangle; it has no answer for one held by a corner.
#[test]
#[should_panic(expected = "frozen whole or not at all")]
fn a_water_frozen_in_part_is_refused() {
    let wb = WaterBox::lattice(2, 1);
    let mut frozen = vec![false; 24];
    frozen[1] = true;
    let _ = wb.dynamics().with_frozen(frozen);
}

// ---------------------------------------------------------------------------------------------
// The box

/// The lattice is at the density it says, every water rigid, every centre of mass on its site,
/// the orientations a function of the seed, and uniform.
#[test]
fn the_lattice_is_at_its_density_and_rigid() {
    let wb = WaterBox::lattice(4, 0x1A7);
    let mass = 64.0 * water::molecular_mass();
    assert!((mass / wb.cell().volume() / water::LIQUID_DENSITY - 1.0).abs() < 8.0 * EPS);
    let side = wb.cell().lengths()[0];
    let spacing = side / 4.0;
    let m = water::masses();
    let p = wb.positions();
    for (w, &[o, a, b]) in wb.waters().iter().enumerate() {
        assert!(rigidity([p[o], p[a], p[b]]) <= RIGID);
        let g = [w / 16, (w / 4) % 4, w % 4];
        for k in 0..3 {
            let c = (m[0] * p[o][k] + m[1] * p[a][k] + m[2] * p[b][k]) / water::molecular_mass();
            assert!((c - (g[k] as f64 + 0.5) * spacing).abs() < 1e-14 * side);
        }
    }
    assert_eq!(WaterBox::lattice(4, 0x1A7), wb);
    assert_ne!(WaterBox::lattice(4, 0x1A8).positions(), wb.positions());
    // Uniform orientations: the bisectors' mean over 4096 waters is zero to its 1/√(3N).
    let big = WaterBox::lattice(16, 5);
    let q = big.positions();
    let (mut mean, mut fourth) = ([0.0; 3], [0.0; 3]);
    for &[o, a, b] in &big.waters() {
        let bis = sub([0, 1, 2].map(|k| 0.5 * (q[a][k] + q[b][k])), q[o]);
        let u = len(bis);
        for k in 0..3 {
            let c = bis[k] / u;
            mean[k] += c / 4096.0;
            fourth[k] += c * c * c * c / 4096.0;
        }
    }
    println!("bisector: mean {mean:.4?}, fourth moment {fourth:.4?} against 1/5");
    // A component of a uniform unit vector has mean 0 and variance 1/3, so 4σ of the mean is
    // 4/√(3·4096); its fourth power has mean 1/5 and variance 1/9 − 1/25, which tells a uniform
    // rotation from one biased toward the axes, as a quaternion uniform in a cube would be.
    for k in 0..3 {
        assert!(mean[k].abs() < 4.0 / (3.0f64 * 4096.0).sqrt(), "{mean:?}");
        let sd = ((1.0 / 9.0 - 1.0 / 25.0) / 4096.0f64).sqrt();
        assert!((fourth[k] - 0.2).abs() < 4.0 * sd, "{fourth:?}");
    }
}

/// W1's guarantees with TIP3P: the forces are the energy's gradient, by central differences to a
/// tolerance from the step's truncation and the energy's rounding; and the energy and forces do
/// not change beyond rounding when every atom is moved by its own whole number of box lengths,
/// which cuts every water.
#[test]
fn the_box_keeps_its_guarantees_with_tip3p() {
    let wb = WaterBox::lattice(3, 0x6A);
    let field = wb.force_field(4.6 * ANGSTROM, 1e-6);
    let at = wb.positions().to_vec();
    let ev = field.evaluate(&at);
    let w = &ev.energy.ewald;
    let scale = ev.energy.van_der_waals.abs()
        + ev.energy.dispersion_correction.abs()
        + w.real.abs()
        + w.reciprocal.abs()
        + w.self_energy.abs()
        + w.excluded.abs();
    let cell = field.cell();
    let rc = field.ewald().parameters().cutoff;
    let h = 1e-5 * ANGSTROM;
    let tested = [(0usize, 0usize), (1, 1), (2, 2), (40, 0), (41, 2), (80, 1)];
    for &(i, k) in &tested {
        // No pair of the moved atom within ten steps of the cutoff, so the largest move, 2h,
        // crosses none: a difference across the cutoff would see its step, not the derivative.
        let margin = (0..at.len())
            .filter(|&j| j != i)
            .map(|j| (len(cell.minimum_image(sub(at[i], at[j]))) - rc).abs())
            .fold(f64::INFINITY, f64::min);
        assert!(
            margin > 10.0 * h,
            "atom {i}: a pair {margin:e} m from the cutoff"
        );
        let e = |s: f64| {
            let mut p = at.clone();
            p[i][k] += s;
            field.energy(&p).total
        };
        let d1 = (e(h) - e(-h)) / (2.0 * h);
        let d2 = (e(2.0 * h) - e(-2.0 * h)) / (4.0 * h);
        let tol = 2.0 * (d2 - d1).abs() / 3.0 + 64.0 * EPS * scale / h;
        assert!(
            (-d1 - ev.forces[i][k]).abs() <= tol,
            "atom {i} axis {k}: {} against {}",
            -d1,
            ev.forces[i][k]
        );
    }
    let l = cell.lengths();
    let moved: Vec<[f64; 3]> = at
        .iter()
        .enumerate()
        .map(|(i, p)| {
            [0, 1, 2].map(|k| {
                let n = (uniform(0x1A, 3 * i as u64 + k as u64) * 7.0).floor() - 3.0;
                p[k] + n * l[k]
            })
        })
        .collect();
    let ev2 = field.evaluate(&moved);
    assert!((ev2.energy.total - ev.energy.total).abs() < 1e-14 * scale);
    let fmax = ev.forces.iter().map(|f| len(*f)).fold(0.0, f64::max);
    for (a, b) in ev2.forces.iter().zip(&ev.forces) {
        assert!(len(sub(*a, *b)) < 1e-12 * fmax);
    }
}
