//! **The forces are minus the gradient of the energy, and the energy does not care where the
//! molecule is or which way it faces.**
//!
//! On aspirin, displaced from the dictionary's geometry by up to 0.05 Å per coordinate so that no
//! bond sits at a special length, and given small partial charges so that the Coulomb term's
//! gradient is exercised too; and on a small probe molecule built to reach what aspirin does not —
//! a torsion inside the collinear switch, a linear centre, and a phosphine's inversion. Every
//! tolerance is computed from the geometry it is applied to, from a per-term bound on rounding
//! and on the third derivative; the derivation is beside each.
//!
//! **The finite-difference test is the one that matters.** Translation and rotation invariance
//! of the energy hold by construction — every term reads only separations — so those checks are
//! cheap regression checks on that construction. One exception, measured in the bond step: the
//! rotation test's *force* comparison works at rounding resolution, far below the
//! finite-difference tolerance, and a sabotage scaling every x-component of the forces by
//! 1 + 1e-9 was caught by it alone.

use pantometry_forcefield::angular::BendForm;
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

/// A probe for what aspirin does not reach: H₂P–CH₂–CH(C≡N)–N=C=O, with the P–C1–C2 angle at 175°
/// so that every torsion about P–C1 and every P–C1–C2–X torsion is inside the collinear switch
/// (170°–180°); the C2–C3≡N and N=C=O angles at about 172° and 175°, off the linear form's
/// minimum; the C2–N2 bond an sp³–sp² bond whose sp² atom has no sp² neighbour, so the six-fold
/// row (b) is under the check; and the phosphorus three-coordinate and pyramidal, so its
/// inversion is on and away from its minimum.
fn probe() -> (ForceField, Vec<[f64; 3]>) {
    let atoms = [
        ("P1", "P", [-1.843, 0.161, 0.02]),
        ("HP1", "H", [-2.09, 1.27, 0.86]),
        ("HP2", "H", [-2.14, 1.06, -1.03]),
        ("C1", "C", [0.0, 0.0, 0.0]),
        ("H11", "H", [0.05, -0.55, 0.94]),
        ("H12", "H", [0.05, -0.6, -0.9]),
        ("C2", "C", [1.53, 0.0, 0.0]),
        ("H21", "H", [1.18, 0.55, 0.9]),
        ("N2", "N", [1.077, 0.776, -1.138]),
        ("C4", "C", [1.977, 1.276, -1.738]),
        ("O1", "O", [2.8896, 1.6972, -2.3347]),
        ("C3", "C", [2.28, -1.25, 0.05]),
        ("N1", "N", [3.009, -2.152, 0.09]),
    ];
    let bonds = [
        ("P1", "HP1", "SING"),
        ("P1", "HP2", "SING"),
        ("P1", "C1", "SING"),
        ("C1", "H11", "SING"),
        ("C1", "H12", "SING"),
        ("C1", "C2", "SING"),
        ("C2", "H21", "SING"),
        ("C2", "N2", "SING"),
        ("N2", "C4", "DOUB"),
        ("C4", "O1", "DOUB"),
        ("C2", "C3", "SING"),
        ("C3", "N1", "TRIP"),
    ];
    let mut s = String::from(
        "data_PRB\n_chem_comp.id PRB\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n\
         _chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.pdbx_aromatic_flag\n\
         _chem_comp_atom.pdbx_model_Cartn_x_ideal\n_chem_comp_atom.pdbx_model_Cartn_y_ideal\n\
         _chem_comp_atom.pdbx_model_Cartn_z_ideal\n",
    );
    for (name, element, [x, y, z]) in atoms {
        s += &format!("PRB {name} {element} 0 N {x} {y} {z}\n");
    }
    s += "loop_\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n\
          _chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n";
    for (a, b, order) in bonds {
        s += &format!("{a} {b} {order} N\n");
    }
    let c = Component::from_ccd(&s).expect("the probe parses");
    let n = c.atoms().len();
    let ff = ForceField::new(&c, &uff::assign(&c))
        .expect("supported")
        .with_charges((0..n).map(|i| 0.2 * (0.9 * i as f64).sin()).collect());
    let at = c.atoms().iter().map(|a| a.at).collect();
    (ff, at)
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

/// The sine of the angle at `j` between `i` and `k`.
fn sin_at(at: &[[f64; 3]], i: usize, j: usize, k: usize) -> f64 {
    let (u, v) = (sub(at[i], at[j]), sub(at[k], at[j]));
    let c = dot(u, v) / (len(u) * len(v));
    (1.0 - c * c).max(0.0).sqrt()
}

/// The angle at `j` between `i` and `k`, degrees.
fn angle_at(at: &[[f64; 3]], i: usize, j: usize, k: usize) -> f64 {
    let (u, v) = (sub(at[i], at[j]), sub(at[k], at[j]));
    (dot(u, v) / (len(u) * len(v)))
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees()
}

fn distance(at: &[[f64; 3]], [i, j]: [usize; 2]) -> f64 {
    len(sub(at[i], at[j]))
}

/// One term's part in the tolerance model, in kcal/mol and Å. Per atom of the term, in the order
/// of `atoms`:
///
/// - `grad`, the magnitude of the term's gradient on that atom;
/// - `m3`, a bound on the term's third derivative along each axis of that atom;
/// - `stiff`, a bound on how fast the term's force on that atom changes when the separations it
///   reads move — what the rotation test needs;
/// - `round_g`, a bound on the rounding of the term's force on that atom;
///
/// and `round_e`, a bound on the rounding of the term's energy.
///
/// **The radial terms** (bond, van der Waals, Coulomb) are worked analytically as in the bond step:
///
/// - `r` itself, from three differences, three squares, two sums and a square root, is good to
///   about `4 ε r`.
/// - **Bond**, `½ k d²` with `d = r − r₀`: `δE ≤ ε (4 k |d| r + |E|)` and
///   `δf' ≤ ε (4 k r + |f'|)` — `k |d| r` is the largest intermediate, and keeps the bound right
///   where `E` → 0.
/// - **Lennard-Jones**: `s = x/r` good to `5 ε`, `s⁶` to `33 ε`, `s¹²` to `67 ε`, relative, so
///   `δE ≤ 70 ε D (s¹² + 2 s⁶)` and `δf' ≤ 75 ε · 12 D (s¹² + s⁶) / r`.
/// - **Coulomb**: `δE ≤ 8 ε |E|`, `δf' ≤ 9 ε |f'|`.
/// - `m3 = |f‴| + 3|f″|/r + 3|f′|/r²` on every axis: the third directional derivative of a radial
///   `f(|r|)` is `f‴c³ + 3c(1−c²)(f″/r − f′/r²)`, `|c| ≤ 1`. `stiff = |f″| + 2|f′|/r`.
///
/// **The angular terms** (bend, torsion, inversion) are measured, through the term's own
/// `energy_at`, rather than differentiated by hand three times:
///
/// - `grad` by a central difference, h = 1e-5 Å.
/// - `m3` by the five-point stencil `(E(2H) − 2E(H) + 2E(−H) − E(−2H)) / 2H³`, H = 1e-3 Å, whose
///   truncation `H²/4 |E⁽⁵⁾|` is about `H²/ρ²` ≈ 1e-6 of it for a length scale ρ ~ 1 Å, and whose
///   rounding is `3 δE / H³`. **Taken twice over, plus that rounding**, the factor two covering the
///   stencil's truncation and the change of `E‴` within h of the base point with room to spare.
///   The truncation part of the finite-difference tolerance is a few percent of it (printed), so
///   this factor barely moves it.
/// - `stiff` as twice the sum over the term's atoms of the Frobenius norm of the Hessian block,
///   measured by the four-point mixed difference at 1e-4 Å, plus its rounding `4 δE / h²`.
/// - Rounding: each is a polynomial in one cosine-like quantity `c` (`cos θ`, `cos φ`, `sin ω`)
///   computed from difference vectors, which are exact to one rounding of themselves. `c` is then
///   good to `20 ε / sin` of the angle between the vectors whose cross product or dot product
///   makes it, so `δE ≤ |dE/dc|_max · 20 ε / sin + 10 ε M`, M the sum of the polynomial's parts'
///   sizes. The force is `dE/dc` times `∇c`, `|∇c| ≤ 2 / (r_min sin)`, each good to a few dozen ε:
///   `δF ≤ 100 ε |dE/dc|_max · 2 / (r_min sin)`.
struct Term {
    atoms: Vec<usize>,
    e: f64,
    round_e: f64,
    grad: Vec<f64>,
    m3: Vec<[f64; 3]>,
    stiff: Vec<f64>,
    round_g: Vec<f64>,
}

#[allow(clippy::too_many_arguments)]
fn radial(
    atoms: [usize; 2],
    r: f64,
    e: f64,
    d1: f64,
    d2: f64,
    d3: f64,
    round_e: f64,
    round_d1: f64,
) -> Term {
    let m3 = d3.abs() + 3.0 * d2.abs() / r + 3.0 * d1.abs() / (r * r);
    let stiff = d2.abs() + 2.0 * d1.abs() / r;
    Term {
        atoms: atoms.to_vec(),
        e,
        round_e,
        grad: vec![d1.abs(); 2],
        m3: vec![[m3; 3]; 2],
        stiff: vec![stiff; 2],
        round_g: vec![round_d1; 2],
    }
}

/// An angular term measured through `energy` (positions in metres → joules), with `de_dc` the
/// largest `|dE/dc|` and `magnitude` the size of its parts (both kcal/mol), and `sin` and `r_min`
/// the conditioning of its geometry.
fn angular(
    atoms: Vec<usize>,
    energy: &dyn Fn(&[[f64; 3]]) -> f64,
    at: &[[f64; 3]],
    de_dc: f64,
    magnitude: f64,
    sin: f64,
    r_min: f64,
) -> Term {
    let e_at = |p: &[[f64; 3]]| energy(p) / KCAL_PER_MOL;
    let e = e_at(at);
    let round_e = de_dc * 20.0 * EPS / sin + 10.0 * EPS * magnitude;
    let round_g = 100.0 * EPS * de_dc * 2.0 / (r_min * sin);
    let moved = |moves: &[(usize, usize, f64)]| {
        let mut p = at.to_vec();
        for &(i, k, d) in moves {
            p[i][k] += d * ANGSTROM;
        }
        e_at(&p)
    };
    let (h, big, h2) = (1e-5, 1e-3, 1e-4);
    let mut grad = Vec::new();
    let mut m3 = Vec::new();
    let mut stiff = Vec::new();
    for &i in &atoms {
        let mut g = [0.0; 3];
        let mut third = [0.0; 3];
        for k in 0..3 {
            g[k] = (moved(&[(i, k, h)]) - moved(&[(i, k, -h)])) / (2.0 * h);
            let s = moved(&[(i, k, 2.0 * big)]) - 2.0 * moved(&[(i, k, big)])
                + 2.0 * moved(&[(i, k, -big)])
                - moved(&[(i, k, -2.0 * big)]);
            third[k] = s.abs() / (2.0 * big * big * big);
        }
        let worst = third.iter().fold(0.0f64, |m, x| m.max(*x));
        m3.push([2.0 * worst + 3.0 * round_e / (big * big * big); 3]);
        grad.push(len(g));
        let mut sum = 0.0;
        for &j in &atoms {
            let mut frob = 0.0;
            for a in 0..3 {
                for b in 0..3 {
                    let hab = (moved(&[(i, a, h2), (j, b, h2)])
                        - moved(&[(i, a, h2), (j, b, -h2)])
                        - moved(&[(i, a, -h2), (j, b, h2)])
                        + moved(&[(i, a, -h2), (j, b, -h2)]))
                        / (4.0 * h2 * h2);
                    frob += hab * hab;
                }
            }
            sum += frob.sqrt();
        }
        stiff.push(2.0 * sum + 4.0 * round_e / (h2 * h2));
    }
    Term {
        round_g: vec![round_g; atoms.len()],
        atoms,
        e,
        round_e,
        grad,
        m3,
        stiff,
    }
}

/// Every term of `ff` at `at`, in kcal/mol and Å.
fn terms(ff: &ForceField, at: &[[f64; 3]]) -> Vec<Term> {
    let mut out = Vec::new();
    for s in ff.stretches() {
        let r = distance(at, s.atoms) / ANGSTROM;
        let k = s.force_constant * ANGSTROM * ANGSTROM / KCAL_PER_MOL;
        let dr = r - s.natural_length / ANGSTROM;
        let (e, d1) = (0.5 * k * dr * dr, k * dr);
        out.push(radial(
            s.atoms,
            r,
            e,
            d1,
            k,
            0.0,
            EPS * (4.0 * k * dr.abs() * r + e.abs()),
            EPS * (4.0 * k * r + d1.abs()),
        ));
    }
    let q = ff.charges();
    for p in ff.pairs() {
        let r = distance(at, p.atoms) / ANGSTROM;
        let (x, d) = (p.distance / ANGSTROM, p.well / KCAL_PER_MOL);
        let (a12, a6) = (d * x.powi(12), 2.0 * d * x.powi(6));
        out.push(radial(
            p.atoms,
            r,
            a12 / r.powi(12) - a6 / r.powi(6),
            -12.0 * a12 / r.powi(13) + 6.0 * a6 / r.powi(7),
            156.0 * a12 / r.powi(14) - 42.0 * a6 / r.powi(8),
            -2184.0 * a12 / r.powi(15) + 336.0 * a6 / r.powi(9),
            70.0 * EPS * (a12 / r.powi(12) + a6 / r.powi(6)),
            75.0 * EPS * (12.0 * a12 / r.powi(13) + 6.0 * a6 / r.powi(7)),
        ));
        let c = 332.0637 * q[p.atoms[0]] * q[p.atoms[1]];
        out.push(radial(
            p.atoms,
            r,
            c / r,
            -c / (r * r),
            2.0 * c / r.powi(3),
            -6.0 * c / r.powi(4),
            8.0 * EPS * (c / r).abs(),
            9.0 * EPS * (c / (r * r)).abs(),
        ));
    }
    let kcal = |j: f64| j / KCAL_PER_MOL;
    let r_of = |i: usize, j: usize| distance(at, [i, j]) / ANGSTROM;
    for b in ff.bends() {
        let [i, j, k] = b.atoms;
        let kk = kcal(b.force_constant);
        let (de_dc, magnitude) = match b.form {
            BendForm::General { c0, c1, c2 } => (
                kk * (c1.abs() + 4.0 * c2.abs()),
                kk * (c0.abs() + c1.abs() + 2.0 * c2.abs()),
            ),
            BendForm::Linear => (kk, 2.0 * kk),
            BendForm::Trigonal => (kk, kk),
        };
        out.push(angular(
            b.atoms.to_vec(),
            &|p| b.energy_at(p),
            at,
            de_dc,
            magnitude,
            sin_at(at, i, j, k).max(1e-3),
            r_of(i, j).min(r_of(j, k)),
        ));
    }
    for t in ff.torsions() {
        let [i, j, k, l] = t.atoms;
        let share = kcal(t.share());
        let n = f64::from(t.parameters.periodicity);
        out.push(angular(
            t.atoms.to_vec(),
            &|p| t.energy_at(p),
            at,
            0.5 * share * n * n,
            2.0 * share,
            sin_at(at, i, j, k).min(sin_at(at, j, k, l)),
            r_of(i, j).min(r_of(j, k)).min(r_of(k, l)),
        ));
    }
    for v in ff.inversions() {
        let [i, j, k, l] = v.atoms;
        let kk = kcal(v.force_constant) / 3.0;
        let [c0, c1, c2] = v.coefficients;
        // |dE/d sin ω| = K/3 |C₁ sin ω / cos ω + 4 C₂ sin ω| ≤ K/3 (|C₁| / cos ω + 4|C₂|), with
        // cos ω from this geometry, halved for the room the finite differences move it.
        let n = {
            let (u, w) = (sub(at[j], at[i]), sub(at[k], at[i]));
            [
                u[1] * w[2] - u[2] * w[1],
                u[2] * w[0] - u[0] * w[2],
                u[0] * w[1] - u[1] * w[0],
            ]
        };
        let axis = sub(at[l], at[i]);
        let sin_w = dot(n, axis) / (len(n) * len(axis));
        let cos_w = 0.5 * (1.0 - sin_w * sin_w).sqrt();
        out.push(angular(
            v.atoms.to_vec(),
            &|p| v.energy_at(p),
            at,
            kk * (c1.abs() / cos_w + 4.0 * c2.abs()),
            kk * (c0.abs() + c1.abs() + c2.abs()),
            sin_at(at, j, i, k),
            r_of(i, j).min(r_of(i, k)).min(r_of(i, l)),
        ));
    }
    out
}

/// A bound on the rounding error of one whole evaluation of the energy: every term's own, plus
/// the running sums, which add at most one ε of a partial sum — itself at most `S = Σ|E|` — per
/// term, twice (the term's own sum and the total).
fn evaluation_rounding(ts: &[Term]) -> f64 {
    let s: f64 = ts.iter().map(|t| t.e.abs()).sum();
    ts.iter().map(|t| t.round_e).sum::<f64>() + 2.0 * ts.len() as f64 * EPS * s
}

fn kcal(e: Energy) -> f64 {
    e.total / KCAL_PER_MOL
}

/// The largest distance of any atom from the origin, Å.
fn extent(at: &[[f64; 3]]) -> f64 {
    at.iter().map(|p| len(*p) / ANGSTROM).fold(0.0, f64::max)
}

/// Each touching term's position within `atoms`, for atom `i`.
fn touching(ts: &[Term], i: usize) -> Vec<(&Term, usize)> {
    ts.iter()
        .filter_map(|t| t.atoms.iter().position(|&a| a == i).map(|p| (t, p)))
        .collect()
}

/// **Analytic force against a central difference of the energy, every coordinate of every atom.**
///
/// The central difference `(E(x+h) − E(x−h)) / 2h` differs from `dE/dx` by truncation
/// `h²/6 · |E‴|` and by rounding `δE / h`, where `δE` bounds the error of one energy evaluation
/// ([`evaluation_rounding`]); `|E‴|` along one axis of atom `i` is bounded by the sum of the
/// touching terms' `m3`.
///
/// **And per atom, the check has to be able to see the atom's terms**: for every atom, at least
/// three quarters of the terms that touch it must have a gradient on it above 100 times that
/// atom's tolerance. So a term wrong by 1% of its own force is visible for most terms on every
/// atom, not only for the one clash that sets the largest force.
fn check_forces(name: &str, ff: &ForceField, at: &[[f64; 3]]) {
    let forces = ff.evaluate(at).forces;
    let ts = terms(ff, at);
    let delta_e = evaluation_rounding(&ts);
    let h = 1e-5;
    let (mut worst_err, mut worst_tol, mut worst_trunc) = (0.0f64, 0.0f64, 0.0f64);
    let mut least_visible = 1.0f64;
    for i in 0..at.len() {
        let near = touching(&ts, i);
        let mut tol = [0.0; 3];
        for (k, tk) in tol.iter_mut().enumerate() {
            let m3: f64 = near.iter().map(|(t, p)| t.m3[*p][k]).sum();
            *tk = h * h / 6.0 * m3 + delta_e / h;
            worst_trunc = worst_trunc.max(h * h / 6.0 * m3 / *tk);
        }
        let tol_i = tol.iter().fold(0.0f64, |m, x| m.max(*x));
        let visible = near
            .iter()
            .filter(|(t, p)| t.grad[*p] > 100.0 * tol_i)
            .count() as f64
            / near.len() as f64;
        least_visible = least_visible.min(visible);
        assert!(
            visible >= 0.75,
            "{name} atom {i}: only {:.0}% of its terms are 100x above tol {tol_i:e}",
            100.0 * visible
        );
        for k in 0..3 {
            let mut plus = at.to_vec();
            let mut minus = at.to_vec();
            plus[i][k] += h * ANGSTROM;
            minus[i][k] -= h * ANGSTROM;
            // The step actually taken, after rounding the displaced coordinate.
            let step = (plus[i][k] - minus[i][k]) / ANGSTROM;
            let numeric = -(kcal(ff.energy(&plus)) - kcal(ff.energy(&minus))) / step;
            let analytic = forces[i][k] * ANGSTROM / KCAL_PER_MOL;
            let err = (numeric - analytic).abs();
            assert!(
                err <= tol[k],
                "{name} atom {i} axis {k}: analytic {analytic} numeric {numeric} err {err:e} \
                 tol {:e}",
                tol[k]
            );
            worst_err = worst_err.max(err);
            worst_tol = worst_tol.max(tol[k]);
        }
    }
    println!(
        "{name}: {} terms; δE {delta_e:.2e} kcal/mol; worst error {worst_err:.2e}, largest \
         tolerance {worst_tol:.2e} kcal/mol/Å, truncation at most {:.0}% of a tolerance; least \
         fraction of an atom's terms 100x above its tolerance {least_visible:.2}",
        ts.len(),
        100.0 * worst_trunc
    );
}

#[test]
fn the_force_is_minus_the_gradient_of_the_energy() {
    let (ff, at) = aspirin();
    assert!(!ff.bends().is_empty() && !ff.torsions().is_empty() && !ff.inversions().is_empty());
    check_forces("aspirin", &ff, &at);
}

/// **The same on the probe, where the switch, the linear form and phosphorus's inversion are
/// live** — and the test first checks that they are, so it cannot pass by not reaching them.
#[test]
fn the_force_is_minus_the_gradient_where_aspirin_does_not_reach() {
    let (ff, at) = probe();
    let switched = ff
        .torsions()
        .iter()
        .filter(|t| {
            let [i, j, k, l] = t.atoms;
            let worst = angle_at(&at, i, j, k).max(angle_at(&at, j, k, l));
            worst > 170.0 && worst < 180.0
        })
        .count();
    assert!(switched >= 4, "{switched} torsions inside the switch");
    let linear = ff
        .bends()
        .iter()
        .filter(|b| b.form == BendForm::Linear)
        .count();
    assert_eq!(linear, 2, "the nitrile's and the isocyanate's");
    let six_fold = ff
        .torsions()
        .iter()
        .filter(|t| t.parameters.periodicity == 6)
        .count();
    assert_eq!(six_fold, 3, "C2-N2, row (b): three torsions");
    assert_eq!(ff.inversions().len(), 3, "phosphorus's");
    println!("probe: {switched} torsions inside the switch");
    check_forces("probe", &ff, &at);
}

/// **Phosphine near its own minimum, where the inversion's `cos ω` is small.** PH₃'s ω₀ is
/// 84.43°, `cos ω₀` ≈ 0.098, and there the inversion's `C₁ sin ω / cos ω` and `4 C₂ sin ω` nearly
/// cancel, so the derivative depends on `cos ω` being used as it is and not clamped. Bonds at
/// about the 93.8° of θ₀, nudged by up to 0.02 rad and 0.02 Å so no term sits exactly at its
/// minimum. The test first checks that every inversion here has `cos ω` below 0.15.
#[test]
fn the_force_is_minus_the_gradient_near_a_phosphines_minimum() {
    let text = [
        "data_PH3",
        "_chem_comp.id PH3",
        "loop_",
        "_chem_comp_atom.comp_id",
        "_chem_comp_atom.atom_id",
        "_chem_comp_atom.type_symbol",
        "_chem_comp_atom.charge",
        "_chem_comp_atom.pdbx_aromatic_flag",
        "_chem_comp_atom.pdbx_model_Cartn_x_ideal",
        "_chem_comp_atom.pdbx_model_Cartn_y_ideal",
        "_chem_comp_atom.pdbx_model_Cartn_z_ideal",
        "PH3 P P 0 N 0 0 0",
        "PH3 H1 H 0 N 1 0 0",
        "PH3 H2 H 0 N 0 1 0",
        "PH3 H3 H 0 N 0 0 1",
        "loop_",
        "_chem_comp_bond.atom_id_1",
        "_chem_comp_bond.atom_id_2",
        "_chem_comp_bond.value_order",
        "_chem_comp_bond.pdbx_aromatic_flag",
        "P H1 SING N",
        "P H2 SING N",
        "P H3 SING N",
        "",
    ]
    .join("\n");
    let c = Component::from_ccd(&text).expect("phosphine parses");
    let ff = ForceField::new(&c, &uff::assign(&c)).expect("supported");
    let theta = 93.8f64.to_radians();
    let alpha = ((theta.cos() + 0.5) / 1.5).sqrt().asin();
    let mut at = vec![[0.0; 3]];
    for (az, de, r) in [(0.0, 0.02, 1.40), (118.0, -0.01, 1.43), (243.0, 0.0, 1.42)] {
        let (a, e) = (f64::to_radians(az), alpha + de);
        at.push([
            r * e.cos() * a.cos() * ANGSTROM,
            r * e.cos() * a.sin() * ANGSTROM,
            r * e.sin() * ANGSTROM,
        ]);
    }
    assert_eq!(ff.inversions().len(), 3);
    for v in ff.inversions() {
        let [i, j, k, l] = v.atoms;
        let (u, w) = (sub(at[j], at[i]), sub(at[k], at[i]));
        let n = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        let axis = sub(at[l], at[i]);
        let s = dot(n, axis) / (len(n) * len(axis));
        let cos_w = (1.0 - s * s).sqrt();
        assert!(cos_w < 0.15, "cos ω = {cos_w}");
    }
    check_forces("phosphine", &ff, &at);
}

/// **Translation leaves the energy unchanged, to rounding.** Shifted by (1.7, −2.3, 0.9) Å. Each
/// shifted coordinate rounds by `ε (|p| + |t|)`, so each separation moves by at most
/// `dr = 2√3 ε (p_max + t_max)`; the energy by the sum over terms and atoms of the gradient times
/// that, plus two evaluations' rounding. A regression check: invariance holds by construction.
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
    let slope: f64 = ts.iter().flat_map(|t| t.grad.iter()).sum();
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
/// `|v| (1 ± 1.5 δ)`, and each rotated coordinate carries the three products' rounding,
/// `3√3 ε |p|`, so each separation moves by at most `Δr = (3δ + 18ε) p_max`. The energy moves by
/// the gradients times `Δr` plus two evaluations' rounding.
///
/// For the forces, `F(Rx)` against `R F(x)` per atom: each term's force on the atom moves by its
/// `stiff · Δr`, each evaluation has its own `round_g`, the per-atom sum adds one ε per term, and
/// applying `R` to `F(x)` costs `(1.5δ + 6ε)` of the atom's summed gradients.
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
    let slope: f64 = ts.iter().flat_map(|t| t.grad.iter()).sum();
    let dr = (3.0 * delta + 18.0 * EPS) * extent(&at);
    let tol = slope * dr + 2.0 * evaluation_rounding(&ts);
    println!(
        "rotation: δ = {delta:.1e}; {a} vs {b}, |Δ| {:.2e}, tol {tol:.2e}",
        (a - b).abs()
    );
    assert!((a - b).abs() <= tol);

    let to_kcal = ANGSTROM / KCAL_PER_MOL;
    for (i, (f0, f1)) in e0.forces.iter().zip(&e1.forces).enumerate() {
        let near = touching(&ts, i);
        let sum_g: f64 = near.iter().map(|(t, p)| t.grad[*p]).sum();
        let tol_i: f64 = near
            .iter()
            .map(|(t, p)| t.stiff[*p] * dr + 2.0 * t.round_g[*p])
            .sum::<f64>()
            + (1.5 * delta + 6.0 * EPS + 2.0 * near.len() as f64 * EPS) * sum_g;
        let rf = apply(*f0);
        for k in 0..3 {
            let diff = (rf[k] - f1[k]).abs() * to_kcal;
            assert!(diff <= tol_i, "atom {i} axis {k}: {diff:e} > {tol_i:e}");
        }
    }
}

/// **The forces sum to zero**: every term adds to its atoms vectors that sum to zero — a radial
/// term one vector and its negative, an angular term the negative of the others' sum on its
/// centre — so what is left is the rounding of the per-atom sums and of that one negation sum, at
/// most one ε of each contribution per addition: `2 (N + n) ε` times the summed gradients.
#[test]
fn the_net_force_is_zero() {
    for (name, (ff, at)) in [("aspirin", aspirin()), ("probe", probe())] {
        let forces = ff.evaluate(&at).forces;
        let ts = terms(&ff, &at);
        let slope: f64 =
            ts.iter().flat_map(|t| t.grad.iter()).sum::<f64>() * KCAL_PER_MOL / ANGSTROM;
        let tol = 2.0 * (ts.len() + at.len()) as f64 * EPS * slope;
        for k in 0..3 {
            let net: f64 = forces.iter().map(|f| f[k]).sum();
            assert!(net.abs() <= tol, "{name} axis {k}: {net:e} > {tol:e}");
        }
    }
}

/// **Aspirin at the dictionary's geometry: all six terms, summing to the total**, and the
/// numbers. `total` is accumulated separately from the six parts, so a term added to one and not
/// the other shows here; the two sums differ by their own summation rounding only, at most
/// `(N + 7) ε S`. And the crate's total agrees with the sum of this test's per-term energies —
/// the radial ones computed here from the parameters, the angular ones from each term's own
/// `energy_at` — to two evaluations' rounding.
///
/// Not asserted beyond that: the numbers are a fact about this file's geometry, not about the
/// force field, and minimisation is what should change them. Measured values are in the commit
/// and the changelog. **Most of the van der Waals is one clash in the dictionary's ideal
/// coordinates**: the acetyl carbonyl oxygen O4 sits 1.645 Å from the ring hydrogen H1, a 1-6
/// pair worth 133.9 kcal/mol by itself.
#[test]
fn aspirins_energy_decomposes_into_its_terms() {
    let c = Component::from_ccd(AIN).expect("AIN parses");
    let ff = ForceField::new(&c, &uff::assign(&c)).expect("supported");
    let at: Vec<[f64; 3]> = c.atoms().iter().map(|a| a.at).collect();
    let e = ff.evaluate(&at).energy;
    let k = |x: f64| x / KCAL_PER_MOL;
    println!(
        "aspirin (AIN ideal), kcal/mol: bond {:.6}  angle {:.6}  torsion {:.6}  inversion {:.6}  \
         van der Waals {:.6}  electrostatic {:.6}  total {:.6}  ({} stretches, {} bends, {} \
         torsions, {} inversions, {} non-bonded pairs)",
        k(e.bond),
        k(e.angle),
        k(e.torsion),
        k(e.inversion),
        k(e.van_der_waals),
        k(e.electrostatic),
        k(e.total),
        ff.stretches().len(),
        ff.bends().len(),
        ff.torsions().len(),
        ff.inversions().len(),
        ff.pairs().len()
    );
    let names = |a: &[usize]| {
        a.iter()
            .map(|&i| c.atoms()[i].name.as_str())
            .collect::<Vec<_>>()
            .join("-")
    };
    let mut seen = Vec::new();
    for t in ff.torsions() {
        let bond = [t.atoms[1], t.atoms[2]];
        if !seen.contains(&bond) {
            seen.push(bond);
            let here: f64 = ff
                .torsions()
                .iter()
                .filter(|u| [u.atoms[1], u.atoms[2]] == bond)
                .map(|u| k(u.energy_at(&at)))
                .sum();
            println!(
                "  torsion about {}: {:?}, V {:.4} kcal/mol, n {}, φ0 {}°, shared {} ways; \
                 {here:.4} kcal/mol here",
                names(&bond),
                t.parameters.case,
                k(t.parameters.barrier),
                t.parameters.periodicity,
                t.parameters.equilibrium.to_degrees(),
                t.torsions_about_bond
            );
        }
    }
    let ts = terms(&ff, &at);
    let s: f64 = ts.iter().map(|t| t.e.abs()).sum();
    let parts = k(e.bond)
        + k(e.angle)
        + k(e.torsion)
        + k(e.inversion)
        + k(e.van_der_waals)
        + k(e.electrostatic);
    assert!((parts - k(e.total)).abs() <= (ts.len() as f64 + 7.0) * EPS * s);
    let mine: f64 = ts.iter().map(|t| t.e).sum();
    assert!((mine - k(e.total)).abs() <= 2.0 * evaluation_rounding(&ts));
    assert_eq!(e.electrostatic, 0.0, "no charges, no electrostatics");
    assert_eq!(ff.stretches().len(), 21);
    assert_eq!(ff.bends().len(), 32, "Σ C(degree, 2) over aspirin's atoms");
    assert_eq!(ff.inversions().len(), 24);
}
