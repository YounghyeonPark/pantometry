//! **The forces are minus the gradient of the energy, and the energy does not care where the
//! molecule is or which way it faces.**
//!
//! On aspirin, displaced from the dictionary's geometry by up to 0.05 Å per coordinate so that no
//! bond sits at a special length, and given small partial charges so that the Coulomb term's
//! gradient is exercised too. Every tolerance is computed from the geometry it is applied to, from
//! a per-term bound on rounding; the derivation is beside each.
//!
//! **The finite-difference test is the one that matters.** Translation and rotation invariance
//! of the energy hold by construction — every term reads only a separation — so those checks are
//! cheap regression checks on that construction. One exception, measured: the rotation test's
//! *force* comparison works at rounding resolution, far below the finite-difference tolerance, and
//! a sabotage scaling every x-component of the forces by 1 + 1e-9 was caught by it alone.

use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::{Component, Energy, ForceField};

const AIN: &str = include_str!("../components/AIN.cif");
const ANGSTROM: f64 = 1e-10;
const EPS: f64 = f64::EPSILON;

/// Aspirin's force field with deterministic, neutral partial charges of up to ±0.3 e, and the
/// dictionary's geometry moved by up to 0.05 Å per coordinate. No randomness: a fixed formula.
fn aspirin() -> (ForceField, Vec<[f64; 3]>) {
    let c = Component::from_ccd(AIN).expect("AIN parses");
    let n = c.atoms().len();
    let raw: Vec<f64> = (0..n).map(|i| 0.3 * (1.7 * i as f64).cos()).collect();
    let mean = raw.iter().sum::<f64>() / n as f64;
    let charges = raw.iter().map(|q| q - mean).collect();
    let ff = ForceField::new(&c, &uff::assign(&c))
        .expect("aspirin is supported")
        .with_charges(charges);
    let at = c
        .atoms()
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let mut p = a.at;
            for (k, x) in p.iter_mut().enumerate() {
                *x += 0.05 * ANGSTROM * (1.3 * (3 * i + k) as f64 + 0.7).sin();
            }
            p
        })
        .collect();
    (ff, at)
}

fn distance(at: &[[f64; 3]], [i, j]: [usize; 2]) -> f64 {
    (0..3)
        .map(|k| (at[i][k] - at[j][k]).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// One radial term's energy and its first three derivatives in `r`, in kcal/mol and Å, computed
/// here from the term's parameters rather than by the crate — and a bound on the rounding error
/// of the crate's evaluation of `e` and of `d1`, from the term's largest intermediate.
///
/// The rounding bounds, each `ε` times a multiple of a magnitude that occurs inside the term:
///
/// - `r` itself, from three differences, three squares, two sums and a square root, is good to
///   about `4 ε r`.
/// - **Bond**, `½ k d²` with `d = r − r₀`: `d` inherits `r`'s absolute error, so
///   `δE ≤ ε (4 k |d| r + |E|)` and `δf' ≤ ε (4 k r + |f'|)`. The `k |d| r` is the largest
///   intermediate, and it is what keeps the bound right at a relaxed geometry where `E` → 0.
/// - **Lennard-Jones**: `s = x/r` is good to `5 ε`, `s⁶` (three products) to `33 ε` and `s¹²` to
///   `67 ε`, relative. So `δE ≤ 70 ε D (s¹² + 2 s⁶)` and `δf' ≤ 75 ε · 12 D (s¹² + s⁶) / r` — the
///   two parts' magnitudes, not their difference, which vanishes at the minimum.
/// - **Coulomb**, `C qᵢ qⱼ / r`: `r`'s 4 ε and three roundings — `δE ≤ 8 ε |E|`, `δf' ≤ 9 ε |f'|`.
struct Term {
    atoms: [usize; 2],
    r: f64,
    e: f64,
    d1: f64,
    d2: f64,
    d3: f64,
    round_e: f64,
    round_d1: f64,
}

/// Every term of `ff` at `at`, in kcal/mol and Å.
fn terms(ff: &ForceField, at: &[[f64; 3]]) -> Vec<Term> {
    let mut out = Vec::new();
    for s in ff.stretches() {
        let r = distance(at, s.atoms) / ANGSTROM;
        let k = s.force_constant * ANGSTROM * ANGSTROM / KCAL_PER_MOL;
        let dr = r - s.natural_length / ANGSTROM;
        let (e, d1) = (0.5 * k * dr * dr, k * dr);
        out.push(Term {
            atoms: s.atoms,
            r,
            e,
            d1,
            d2: k,
            d3: 0.0,
            round_e: EPS * (4.0 * k * dr.abs() * r + e.abs()),
            round_d1: EPS * (4.0 * k * r + d1.abs()),
        });
    }
    let q = ff.charges();
    for p in ff.pairs() {
        let r = distance(at, p.atoms) / ANGSTROM;
        let (x, d) = (p.distance / ANGSTROM, p.well / KCAL_PER_MOL);
        let (a12, a6) = (d * x.powi(12), 2.0 * d * x.powi(6));
        out.push(Term {
            atoms: p.atoms,
            r,
            e: a12 / r.powi(12) - a6 / r.powi(6),
            d1: -12.0 * a12 / r.powi(13) + 6.0 * a6 / r.powi(7),
            d2: 156.0 * a12 / r.powi(14) - 42.0 * a6 / r.powi(8),
            d3: -2184.0 * a12 / r.powi(15) + 336.0 * a6 / r.powi(9),
            round_e: 70.0 * EPS * (a12 / r.powi(12) + a6 / r.powi(6)),
            round_d1: 75.0 * EPS * (12.0 * a12 / r.powi(13) + 6.0 * a6 / r.powi(7)),
        });
        let c = 332.0637 * q[p.atoms[0]] * q[p.atoms[1]];
        out.push(Term {
            atoms: p.atoms,
            r,
            e: c / r,
            d1: -c / (r * r),
            d2: 2.0 * c / r.powi(3),
            d3: -6.0 * c / r.powi(4),
            round_e: 8.0 * EPS * (c / r).abs(),
            round_d1: 9.0 * EPS * (c / (r * r)).abs(),
        });
    }
    out
}

/// A bound on the rounding error of one whole evaluation of the energy: every term's own, plus
/// the running sum, which adds at most one ε of a partial sum — itself at most `S = Σ|E|` — per
/// term.
fn evaluation_rounding(ts: &[Term]) -> f64 {
    let s: f64 = ts.iter().map(|t| t.e.abs()).sum();
    ts.iter().map(|t| t.round_e).sum::<f64>() + ts.len() as f64 * EPS * s
}

fn kcal(e: Energy) -> f64 {
    e.total / KCAL_PER_MOL
}

/// The largest distance of any atom from the origin, Å.
fn extent(at: &[[f64; 3]]) -> f64 {
    at.iter()
        .map(|p| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() / ANGSTROM)
        .fold(0.0, f64::max)
}

/// **Analytic force against a central difference of the energy, every coordinate of every atom.**
///
/// The central difference `(E(x+h) − E(x−h)) / 2h` differs from `dE/dx` by truncation
/// `h²/6 · |E'''|` and by rounding `δE / h`, where `δE` bounds the error of one energy evaluation
/// ([`evaluation_rounding`]).
///
/// - `|E'''|` along one coordinate of atom `i` is bounded by the sum, over the terms that touch
///   `i`, of `|f'''| + 3|f''|/r + 3|f'|/r²` — the third directional derivative of a radial
///   function `f(|r|)`, `f'''c³ + 3c(1−c²)(f''/r − f'/r²)`, with `|c| ≤ 1`. Evaluated at the base
///   point; within `h` = 1e-5 Å of it the derivatives change by a relative 1e-5 or less.
///
/// **And per atom, the check has to be able to see the atom's terms**: for every atom, at least
/// three quarters of the terms that touch it must have `|f'|` above 100 times that atom's
/// tolerance. So a term wrong by 1% of its own force is visible for most terms on every atom, not
/// only for the one clash that sets the largest force.
#[test]
fn the_force_is_minus_the_gradient_of_the_energy() {
    let (ff, at) = aspirin();
    let forces = ff.evaluate(&at).forces;
    let ts = terms(&ff, &at);
    let delta_e = evaluation_rounding(&ts);
    let h = 1e-5;
    let (mut worst_err, mut worst_tol) = (0.0f64, 0.0f64);
    let mut least_visible = 1.0f64;
    for i in 0..at.len() {
        let touching: Vec<&Term> = ts.iter().filter(|t| t.atoms.contains(&i)).collect();
        let m3: f64 = touching
            .iter()
            .map(|t| t.d3.abs() + 3.0 * t.d2.abs() / t.r + 3.0 * t.d1.abs() / (t.r * t.r))
            .sum();
        let tol = h * h / 6.0 * m3 + delta_e / h;
        let visible = touching.iter().filter(|t| t.d1.abs() > 100.0 * tol).count() as f64
            / touching.len() as f64;
        least_visible = least_visible.min(visible);
        assert!(
            visible >= 0.75,
            "atom {i}: only {:.0}% of its terms are 100x above tol {tol:e}",
            100.0 * visible
        );
        for k in 0..3 {
            let mut plus = at.clone();
            let mut minus = at.clone();
            plus[i][k] += h * ANGSTROM;
            minus[i][k] -= h * ANGSTROM;
            // The step actually taken, after rounding the displaced coordinate.
            let step = (plus[i][k] - minus[i][k]) / ANGSTROM;
            let numeric = -(kcal(ff.energy(&plus)) - kcal(ff.energy(&minus))) / step;
            let analytic = forces[i][k] * ANGSTROM / KCAL_PER_MOL;
            let err = (numeric - analytic).abs();
            assert!(
                err <= tol,
                "atom {i} axis {k}: analytic {analytic} numeric {numeric} err {err:e} tol {tol:e}"
            );
            worst_err = worst_err.max(err);
            worst_tol = worst_tol.max(tol);
        }
    }
    println!(
        "δE {delta_e:.2e} kcal/mol; worst error {worst_err:.2e}, largest tolerance \
         {worst_tol:.2e} kcal/mol/Å; least fraction of an atom's terms 100x above its tolerance \
         {least_visible:.2}"
    );
}

/// **Translation leaves the energy unchanged, to rounding.** Shifted by (1.7, −2.3, 0.9) Å. Each
/// shifted coordinate rounds by `ε (|p| + |t|)`, so each separation moves by at most
/// `2√3 ε (p_max + t_max)`; the energy by `Σ|f'|` times that, plus two evaluations' rounding.
/// A regression check: invariance holds by construction (see the file's documentation).
#[test]
fn translation_does_not_change_the_energy() {
    let (ff, at) = aspirin();
    let t = [1.7, -2.3, 0.9];
    let moved: Vec<[f64; 3]> = at
        .iter()
        .map(|p| {
            [
                p[0] + t[0] * ANGSTROM,
                p[1] + t[1] * ANGSTROM,
                p[2] + t[2] * ANGSTROM,
            ]
        })
        .collect();
    let (a, b) = (kcal(ff.energy(&at)), kcal(ff.energy(&moved)));
    let ts = terms(&ff, &at);
    let slope: f64 = ts.iter().map(|t| t.d1.abs()).sum();
    let t_max = t.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    let dr = 2.0 * 3f64.sqrt() * EPS * (extent(&at) + t_max);
    let tol = slope * dr + 2.0 * evaluation_rounding(&ts);
    println!(
        "translation: {a} vs {b}, |Δ| {:.2e}, tol {tol:.2e}",
        (a - b).abs()
    );
    assert!((a - b).abs() <= tol);
}

/// **Rotation leaves the energy unchanged, to rounding, and turns the forces with it.**
///
/// The rotation is three Euler angles multiplied out, so the matrix `R` is orthonormal only to
/// `δ = max|RᵀR − I|`, measured here. A separation `v` then comes out with length
/// `|v| (1 ± 1.5 δ)` (since `|vᵀ(RᵀR − I)v| ≤ 3δ|v|²`), and each rotated coordinate carries the
/// three products' rounding, `3√3 ε |p|`, so each separation moves by at most
/// `Δr = (3δ + 18ε) p_max`. The energy moves by `Σ|f'| Δr` plus two evaluations' rounding.
///
/// For the forces, `F(Rx)` against `R F(x)` per atom: each term's `f'` moves by `|f''| Δr` and its
/// direction by `2 Δr / r`, each evaluation has its own `δf'`, the per-atom sum adds one ε per
/// term, and applying `R` to `F(x)` costs `(1.5δ + 6ε)` of the atom's `Σ|f'|`. A regression check:
/// invariance holds by construction (see the file's documentation).
#[test]
fn rotation_does_not_change_the_energy() {
    let (ff, at) = aspirin();
    let (a, b, c) = (0.7f64, -1.1f64, 2.3f64);
    let rz = |t: f64| {
        [
            [t.cos(), -t.sin(), 0.0],
            [t.sin(), t.cos(), 0.0],
            [0.0, 0.0, 1.0],
        ]
    };
    let rx = |t: f64| {
        [
            [1.0, 0.0, 0.0],
            [0.0, t.cos(), -t.sin()],
            [0.0, t.sin(), t.cos()],
        ]
    };
    let mul = |m: [[f64; 3]; 3], n: [[f64; 3]; 3]| {
        let mut o = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                o[i][j] = (0..3).map(|k| m[i][k] * n[k][j]).sum();
            }
        }
        o
    };
    let r = mul(rz(a), mul(rx(b), rz(c)));
    let mut delta = 0.0f64;
    for i in 0..3 {
        for j in 0..3 {
            let rtr: f64 = (0..3).map(|k| r[k][i] * r[k][j]).sum();
            let id = if i == j { 1.0 } else { 0.0 };
            delta = delta.max((rtr - id).abs());
        }
    }
    let apply = |p: [f64; 3]| {
        let mut o = [0.0; 3];
        for (i, oi) in o.iter_mut().enumerate() {
            *oi = (0..3).map(|k| r[i][k] * p[k]).sum();
        }
        o
    };
    let turned: Vec<[f64; 3]> = at.iter().map(|&p| apply(p)).collect();
    let (e0, e1) = (ff.evaluate(&at), ff.evaluate(&turned));
    let (a, b) = (kcal(e0.energy), kcal(e1.energy));
    let ts = terms(&ff, &at);
    let slope: f64 = ts.iter().map(|t| t.d1.abs()).sum();
    let dr = (3.0 * delta + 18.0 * EPS) * extent(&at);
    let tol = slope * dr + 2.0 * evaluation_rounding(&ts);
    println!(
        "rotation: δ = {delta:.1e}; {a} vs {b}, |Δ| {:.2e}, tol {tol:.2e}",
        (a - b).abs()
    );
    assert!((a - b).abs() <= tol);

    let to_kcal = ANGSTROM / KCAL_PER_MOL;
    for (i, (f0, f1)) in e0.forces.iter().zip(&e1.forces).enumerate() {
        let touching: Vec<&Term> = ts.iter().filter(|t| t.atoms.contains(&i)).collect();
        let sum_f1: f64 = touching.iter().map(|t| t.d1.abs()).sum();
        let tol_i: f64 = touching
            .iter()
            .map(|t| t.d2.abs() * dr + 2.0 * t.d1.abs() * dr / t.r + 2.0 * t.round_d1)
            .sum::<f64>()
            + (1.5 * delta + 6.0 * EPS + 2.0 * touching.len() as f64 * EPS) * sum_f1;
        let rf = apply(*f0);
        for k in 0..3 {
            let diff = (rf[k] - f1[k]).abs() * to_kcal;
            assert!(diff <= tol_i, "atom {i} axis {k}: {diff:e} > {tol_i:e}");
        }
    }
}

/// **The forces sum to zero**: every term adds one vector to one atom and subtracts the same
/// vector from another, so what is left is the rounding of the per-atom sums — at most one ε of
/// each contribution per addition, bounded by `2 (N + n) ε Σ|f'|`.
#[test]
fn the_net_force_is_zero() {
    let (ff, at) = aspirin();
    let forces = ff.evaluate(&at).forces;
    let ts = terms(&ff, &at);
    let slope: f64 = ts.iter().map(|t| t.d1.abs()).sum::<f64>() * KCAL_PER_MOL / ANGSTROM;
    let tol = 2.0 * (ts.len() + at.len()) as f64 * EPS * slope;
    for k in 0..3 {
        let net: f64 = forces.iter().map(|f| f[k]).sum();
        assert!(net.abs() <= tol, "axis {k}: {net:e} > {tol:e}");
    }
}

/// **Aspirin at the dictionary's geometry: the decomposition sums to the total**, and the
/// numbers. `total` is accumulated separately from the three parts, so a term added to one and
/// not the other shows here; the two sums differ by their own summation rounding only, at most
/// `(N + 4) ε S`.
///
/// Measured: bond 2.139, van der Waals 173.276, electrostatic 0, total 175.415 kcal/mol, over 21
/// stretches and 157 pairs. **Most of the van der Waals is one clash in the dictionary's ideal
/// coordinates**: the acetyl carbonyl oxygen O4 sits 1.645 Å from the ring hydrogen H1 — a 1-6
/// pair whose `x_IJ` is 3.178 Å — and that pair alone is 133.9 kcal/mol. In the entry's model
/// (crystal) coordinates the two are about 4.2 Å apart. Not asserted: it is a fact about this
/// file's geometry, not about the force field, and minimisation is what should remove it.
#[test]
fn aspirins_energy_decomposes_into_its_terms() {
    let c = Component::from_ccd(AIN).expect("AIN parses");
    let ff = ForceField::new(&c, &uff::assign(&c)).expect("supported");
    let at: Vec<[f64; 3]> = c.atoms().iter().map(|a| a.at).collect();
    let e = ff.evaluate(&at).energy;
    let k = |x: f64| x / KCAL_PER_MOL;
    println!(
        "aspirin (AIN ideal), kcal/mol: bond {:.6}  van der Waals {:.6}  electrostatic {:.6}  \
         total {:.6}  ({} stretches, {} non-bonded pairs)",
        k(e.bond),
        k(e.van_der_waals),
        k(e.electrostatic),
        k(e.total),
        ff.stretches().len(),
        ff.pairs().len()
    );
    let ts = terms(&ff, &at);
    let s: f64 = ts.iter().map(|t| t.e.abs()).sum();
    let parts = k(e.bond) + k(e.van_der_waals) + k(e.electrostatic);
    assert!((parts - k(e.total)).abs() <= (ts.len() as f64 + 4.0) * EPS * s);
    // And the crate's sum agrees with the one this test computed from the parameters: two
    // evaluations, each within its rounding bound of the exact value.
    let mine: f64 = ts.iter().map(|t| t.e).sum();
    assert!((mine - k(e.total)).abs() <= 2.0 * evaluation_rounding(&ts));
    assert_eq!(e.electrostatic, 0.0, "no charges, no electrostatics");
    assert_eq!(ff.stretches().len(), 21);
}
