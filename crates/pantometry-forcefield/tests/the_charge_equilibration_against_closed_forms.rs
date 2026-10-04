//! **QEq's pieces against closed forms**: the Slater Coulomb integral against the exact results
//! that are known for it, the solver against the paper's two-atom solution (eq 18) and its
//! clamping rule (eqs 5 and 13), and the whole against the conservation and invariance it must
//! have. The comparison with the paper's own tables is `the_charge_equilibration_paper.rs`.

mod common;

use common::{index, positions};
use pantometry_forcefield::qeq::{
    self, charge_range, coulomb_integral, coulomb_matrix, principal_quantum_number,
    shielded_coulomb, slater_exponent, table_i, Qeq, BOHR, HARTREE_EV,
};
use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::{Component, Element, ForceField};

const ANGSTROM: f64 = 1e-10;
const EPS: f64 = f64::EPSILON;
const AIN: &str = include_str!("../components/AIN.cif");

fn aspirin() -> (Component, Vec<Element>, Vec<[f64; 3]>) {
    let c = Component::from_ccd(AIN).expect("AIN parses");
    let e = c.atoms().iter().map(|a| a.element).collect();
    let at = positions(&c);
    (c, e, at)
}

/// The Coulomb self-repulsion of a 1s Slater density is `5ζ/8` hartree, and of a 2s density
/// `93ζ/256` — the standard results (Roothaan, *J. Chem. Phys.* 19, 1445 (1951)). Two different
/// 1s exponents at one centre give `ζ_a ζ_b (ζ_a² + 3 ζ_a ζ_b + ζ_b²) / (ζ_a + ζ_b)³`, which is
/// `5ζ/8` when they are equal.
///
/// Tolerance: the `R = 0` branch is a sum of at most four positive-and-negative terms of size
/// ≤ ζ, so a few ulps of ζ; 8 ε ζ.
#[test]
fn at_one_centre_the_integral_is_the_textbook_self_repulsion() {
    for zeta in [0.5, 1.0, 1.0698, 2.0, 3.7] {
        let j1 = coulomb_integral(1, zeta, 1, zeta, 0.0);
        assert!(
            (j1 - 5.0 / 8.0 * zeta).abs() <= 8.0 * EPS * zeta,
            "1s {zeta}: {j1}"
        );
        let j2 = coulomb_integral(2, zeta, 2, zeta, 0.0);
        assert!(
            (j2 - 93.0 / 256.0 * zeta).abs() <= 8.0 * EPS * zeta,
            "2s {zeta}: {j2}"
        );
    }
    for (a, b) in [(1.0, 2.0), (0.7, 1.5), (1.0698, 0.8715)] {
        let exact = a * b * (a * a + 3.0 * a * b + b * b) / (a + b) * 1.0 / ((a + b) * (a + b));
        for j in [
            coulomb_integral(1, a, 1, b, 0.0),
            coulomb_integral(1, b, 1, a, 0.0),
        ] {
            assert!(
                (j - exact).abs() <= 8.0 * EPS * exact,
                "1s {a} {b}: {j} vs {exact}"
            );
        }
    }
}

/// `R → 0` is continuous with the one-centre value: the integral is even and smooth in `R`, so
/// `J(R) − J(0)` is `O(R²)` — measured, it is `−c R²` with `c` of order ζ³ — and the general form's
/// `ε/R` loss is below that down to `R = 10⁻³` bohr. Bound: `ζ³ R² + 64 ε / R`.
#[test]
fn the_integral_is_continuous_at_zero_distance() {
    for (n, zeta) in [(1, 1.0), (2, 0.8715), (3, 0.93165), (5, 1.0917)] {
        let j0 = coulomb_integral(n, zeta, n, zeta, 0.0);
        for r in [1e-1, 1e-2, 1e-3] {
            let j = coulomb_integral(n, zeta, n, zeta, r);
            let bound = zeta.powi(3) * r * r + 64.0 * EPS / r;
            assert!(
                (j - j0).abs() <= bound && j <= j0,
                "n = {n}, R = {r}: {j} against {j0} (bound {bound:e})"
            );
        }
    }
}

/// Two equal 1s exponents at any distance have the closed form
/// `J(R) = 1/R − e^(−2ζR) (1/R + 11ζ/8 + 3ζ²R/4 + ζ³R²/6)` (Roothaan 1951).
///
/// Tolerance: the closed form loses `ε/R` to its own cancellation and the general sum here about
/// as much; `64 ε (1 + 1/(ζR)) J`, which bounds both above `R = 0.5` bohr.
#[test]
fn two_equal_1s_densities_at_any_distance_are_roothaans_closed_form() {
    for zeta in [0.6f64, 1.0, 1.0698, 1.8] {
        for r in [0.5f64, 1.0, 1.7, 2.5, 4.0, 8.0, 15.0, 30.0] {
            let e = (-2.0 * zeta * r).exp();
            let exact = 1.0 / r
                - e * (1.0 / r
                    + 11.0 * zeta / 8.0
                    + 0.75 * zeta * zeta * r
                    + zeta.powi(3) * r * r / 6.0);
            let j = coulomb_integral(1, zeta, 1, zeta, r);
            let bound = 64.0 * EPS * (1.0 + 1.0 / (zeta * r)) * exact;
            assert!(
                (j - exact).abs() <= bound,
                "ζ {zeta}, R {r}: {j} vs {exact}"
            );
        }
    }
}

/// Far apart, the shielding vanishes and `J → 1/R` in atomic units: `14.39964/R` eV with R in Å,
/// the paper's eq 14 `14.4/R` to the three figures it prints.
///
/// The correction to `1/R` is exponentially small, bounded by the charge of either density
/// outside radius `R/2` (each density sees the other as a point charge where neither overlaps):
/// for every `n` up to 5 and ζ ≥ 0.8 at `R` = 120 bohr that is below `e^(−0.8·120)·poly ≪ ε`, so
/// `J·R = 1` to rounding: 8 ε.
#[test]
fn far_apart_the_integral_is_the_bare_coulomb_law() {
    for na in 1..=5 {
        for nb in 1..=5 {
            for (za, zb) in [(0.8, 0.8), (0.84034, 1.0917), (2.0, 0.85)] {
                let r = 120.0;
                let j = coulomb_integral(na, za, nb, zb, r);
                assert!(
                    (j * r - 1.0).abs() <= 8.0 * EPS,
                    "{na} {nb}: J·R = {}",
                    j * r
                );
            }
        }
    }
    // In the paper's units.
    let ev_angstrom = HARTREE_EV * BOHR / ANGSTROM;
    assert!((ev_angstrom - 14.399_645).abs() < 1e-6, "{ev_angstrom}");
    assert!((ev_angstrom - 14.4).abs() < 0.05, "14.4 to three figures");
    let r = 60.0 * ANGSTROM;
    let j = shielded_coulomb(
        Element::C,
        slater_exponent(Element::C),
        Element::O,
        slater_exponent(Element::O),
        r,
    );
    assert!(
        (j * 60.0 / ev_angstrom - 1.0).abs() <= 8.0 * EPS,
        "{}",
        j * 60.0
    );
}

/// The integral is symmetric — `∫ρ_A V_B = ∫ρ_B V_A` — though it is computed asymmetrically, A's
/// density against B's potential, so swapping the atoms runs different arithmetic; and it is
/// positive and below the bare `1/R` everywhere, since overlap only ever shields.
///
/// Tolerance for symmetry, derived: the shielding term `K ≤ 1/R` is an alternating binomial
/// expansion of `(u+v)^(2m−1) (u−v)^k`, whose terms can exceed their sum by up to `2^(2m+2n−2)`, so
/// each order is good to about `2^(2m+2n−2) ε / R` hartree; the bound is eight times that. The
/// worst measured is printed.
#[test]
fn the_integral_is_symmetric_positive_and_shielded() {
    let mut worst = 0.0f64;
    for na in 1..=5u32 {
        for nb in 1..=5u32 {
            for (za, zb) in [
                (1.0698, 0.8715),
                (0.84034, 1.0917),
                (0.95, 0.95),
                (2.07, 0.88),
            ] {
                for r in [0.3, 1.0, 1.8, 2.6, 4.0, 7.5, 20.0, 60.0] {
                    let ab = coulomb_integral(na, za, nb, zb, r);
                    let ba = coulomb_integral(nb, zb, na, za, r);
                    let bound = 8.0 * 2f64.powi(2 * (na + nb) as i32 - 2) * EPS / r;
                    worst = worst.max((ab - ba).abs() / bound);
                    assert!(
                        (ab - ba).abs() <= bound,
                        "{na}s({za}) {nb}s({zb}) at {r}: {ab} vs {ba}"
                    );
                    // Strictly below where the overlap is not lost in rounding; equal to 1/R,
                    // to rounding, at 60 bohr.
                    let below = if r < 30.0 {
                        ab < 1.0 / r
                    } else {
                        ab <= (1.0 + 4.0 * EPS) / r
                    };
                    assert!(ab > 0.0 && below, "{na} {nb} {r}: {ab}");
                }
            }
        }
    }
    println!("largest asymmetry, as a fraction of its bound: {worst:.3e}");
}

/// The exponents near each other — the case the textbook closed form, which divides by
/// `ζ_A − ζ_B`, loses everything in — agree with the equal-exponent value to first order:
/// `|J(ζ, ζ + h) − J(ζ, ζ)| ≤ |∂J/∂ζ_B| h`, and `|∂J/∂ζ| ≤ J/ζ·(1 + ζR)` for a density that
/// scales as ζ; the second-order term is below the first by a further `h`.
#[test]
fn nearly_equal_exponents_lose_nothing() {
    for n in 1..=5u32 {
        for r in [1.0, 2.5, 6.0] {
            let z = 0.9;
            let j = coulomb_integral(n, z, n, z, r);
            for h in [1e-4, 1e-8, 1e-12] {
                let jh = coulomb_integral(n, z, n, z + h, r);
                let bound = j / z * (1.0 + z * r) * h + 16.0 * EPS * j;
                assert!((jh - j).abs() <= bound, "n {n}, R {r}, h {h}: {jh} vs {j}");
            }
        }
    }
}

/// The matrix of eq 6′: symmetric, every element positive, the diagonal Table I's
/// idempotentials (hydrogen's by eq 21 at its charge), and every off-diagonal J the shielded
/// integral at that pair's distance with hydrogen's ζ by eq 20 — **including the J between a
/// hydrogen and another atom**, which a model with only `J_HH` charge-dependent would get wrong.
/// And positive definite for aspirin at its QEq charges, by Cholesky: eq 6′ then has one minimum.
#[test]
fn the_coulomb_matrix_is_the_papers() {
    let (c, elements, at) = aspirin();
    let q = qeq::charges(&c, &at).expect("aspirin equilibrates").charges;
    let n = elements.len();
    let j = coulomb_matrix(&elements, &at, &q);
    for a in 0..n {
        let row = table_i(elements[a]);
        let diagonal = if elements[a] == Element::H {
            (1.0 + q[a] / slater_exponent(Element::H)) * row.idempotential
        } else {
            row.idempotential
        };
        assert_eq!(j[a * n + a], diagonal, "J_{a}{a}");
        let zeta = |i: usize| {
            slater_exponent(elements[i]) + if elements[i] == Element::H { q[i] } else { 0.0 }
        };
        for b in 0..n {
            assert_eq!(j[a * n + b], j[b * n + a]);
            assert!(j[a * n + b] > 0.0);
            if a < b {
                let d: f64 = (0..3)
                    .map(|k| (at[a][k] - at[b][k]).powi(2))
                    .sum::<f64>()
                    .sqrt();
                let expected = shielded_coulomb(elements[a], zeta(a), elements[b], zeta(b), d);
                assert_eq!(j[a * n + b], expected, "J_{a}{b}");
                assert!(
                    j[a * n + b] < HARTREE_EV * BOHR / d,
                    "shielded below 14.4/R"
                );
            }
        }
    }
    // Cholesky.
    let mut l = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..=i {
            let s: f64 = (0..k).map(|p| l[i * n + p] * l[k * n + p]).sum();
            if i == k {
                let d = j[i * n + i] - s;
                assert!(d > 0.0, "not positive definite at {i}: {d}");
                l[i * n + i] = d.sqrt();
            } else {
                l[i * n + k] = (j[i * n + k] - s) / l[k * n + k];
            }
        }
    }
}

/// The principal quantum numbers, and eq 17′'s exponents against Table I's printed column: the
/// printed ζ is eq 17 with λ = 0.4913 to its four figures for eight of the ten elements,
/// hydrogen's is λ = ½, and oxygen's is neither — the discrepancy the module documents. "Fits"
/// allows the rounding of ζ, of λ = 0.4913 and of the paper's a₀; nitrogen needs λ's (its printed
/// 0.9089 is 1.3e-4 from what 0.4913 gives with CODATA's a₀).
#[test]
fn the_exponents_and_the_printed_column() {
    use Element::*;
    for (e, n) in [
        (H, 1),
        (C, 2),
        (N, 2),
        (O, 2),
        (F, 2),
        (P, 3),
        (S, 3),
        (Cl, 3),
        (Br, 4),
        (I, 5),
    ] {
        assert_eq!(principal_quantum_number(e), n, "{e}");
        let row = table_i(e);
        let r_bohr = row.radius * ANGSTROM / BOHR;
        let at = |lambda: f64| lambda * (2.0 * f64::from(n) + 1.0) / (2.0 * r_bohr);
        assert_eq!(slater_exponent(e), at(0.5), "{e}: eq 17′");
        // Half a unit in ζ's fourth decimal, plus λ's own rounding (0.4913 ± 0.00005, 1.02e-4
        // relative) and the paper's rounding of a₀ (1.4e-5 relative).
        let printed = 0.00005 + (1.02e-4 + 2e-5) * row.printed_zeta;
        let fits = |lambda: f64| (at(lambda) - row.printed_zeta).abs() <= printed;
        match e {
            H => assert!(fits(0.5) && !fits(0.4913), "H is λ = ½"),
            O => assert!(!fits(0.5) && !fits(0.4913), "O is neither"),
            _ => assert!(fits(0.4913) && !fits(0.5), "{e} is λ = 0.4913"),
        }
    }
}

/// The GMP electronegativities UFF's `r_EN` uses — transcribed from Open Babel before the QEq
/// paper was read — are Table I's χ, all ten, exactly: Table I prints them, citing the GMP paper
/// (its ref 9). One table, so eq 4 and QEq cannot disagree about an element.
#[test]
fn uffs_electronegativities_are_table_i() {
    for e in Element::ALL {
        assert_eq!(
            uff::gmp_electronegativity(e),
            table_i(e).electronegativity,
            "{e}"
        );
    }
}

/// Eq 5's printed ranges, and eq 5′'s.
#[test]
fn the_charge_ranges_are_eq_5s() {
    assert_eq!(charge_range(Element::C), (-4.0, 4.0));
    assert_eq!(charge_range(Element::O), (-2.0, 6.0));
    assert_eq!(charge_range(Element::H), (-1.0, 1.0));
}

/// The one-centre self-repulsion of an `ns` Slater density for every `n` the crate uses, against
/// an exact closed form derived a different way from the code's. With the radial density
/// `p(r) = α^(2n+1) r^(2n) e^(−αr)/(2n)!` (α = 2ζ) and its enclosed charge
/// `q(r) = 1 − e^(−αr) Σ_{k≤2n} (αr)^k/k!`, the self-repulsion is `2 ∫ p(r) q(r)/r dr` (each pair
/// counted once, the outer charge seeing the inner as a point), which is
///
/// `J/ζ = 4/(2n)! · [(2n−1)! − Σ_{k=0}^{2n} (2n−1+k)! / (k! 2^(2n+k))]`
///
/// — 5/8, 93/256, 793/3072, 26333/131072 and 43191/262144 for n = 1 to 5, as exact fractions
/// (worked in rational arithmetic; the first two are Roothaan's). The code computes the same
/// number as `∫ ρ_A V_B` with B's potential written out, a different sum.
///
/// Tolerance: the code's `R = 0` branch is a difference of two terms of size ≤ ζ; 16 ε ζ.
#[test]
fn every_ns_self_repulsion_is_the_exact_fraction() {
    let exact = [
        5.0 / 8.0,
        93.0 / 256.0,
        793.0 / 3072.0,
        26333.0 / 131072.0,
        43191.0 / 262144.0,
    ];
    let fact = |k: u32| (1..=k).fold(1.0f64, |f, i| f * f64::from(i));
    for n in 1..=5u32 {
        // The formula in floating point, as a check on the fractions' transcription.
        let sum: f64 = (0..=2 * n)
            .map(|k| fact(2 * n - 1 + k) / (fact(k) * 2f64.powi((2 * n + k) as i32)))
            .sum();
        let formula = 4.0 / fact(2 * n) * (fact(2 * n - 1) - sum);
        let c = exact[n as usize - 1];
        assert!(
            (formula - c).abs() <= 16.0 * EPS,
            "n = {n}: formula {formula} vs {c}"
        );
        for zeta in [0.4364, 0.8715, 1.0917, 2.5] {
            let j = coulomb_integral(n, zeta, n, zeta, 0.0);
            assert!(
                (j - c * zeta).abs() <= 16.0 * EPS * zeta,
                "{n}s, ζ = {zeta}: {j} against {}",
                c * zeta
            );
        }
    }
}

/// **The `|δ|R > 50` branch** — the antiderivative, used when the two exponents differ by a lot
/// at a distance — against references computed independently of this code:
///
/// - **1s–1s**, the closed form for unequal exponents, derived for this test with SymPy
///   (`∫ p_A(r) ⟨V_B⟩(r) dr` with the 1s potential, integrated symbolically, then evaluated at
///   40 digits with mpmath): 1s(10)–1s(0.5) at 6 bohr 0.16500795281709290, 1s(12)–1s(1) at 5 bohr
///   0.19994488423838100, 1s(9)–1s(0.6) at 7 bohr 0.14268889315024881.
/// - **Higher n**, a Fourier-space quadrature at 40 digits written by the numerics review
///   (`J = (2/π) ∫ F_A(k) F_B(k) sin(kR)/(kR) dk` with the densities' closed-form transforms),
///   which shares no code with this crate: 1s(10)–5s(0.5) at 6 bohr 0.098705349023525284,
///   2s(9)–3s(0.6) at 7 bohr 0.13612917504237427, 4s(8)–2s(0.7) at 8 bohr 0.12481847900814274.
///   The same quadrature gives the first 1s–1s point as 0.16500795281706081, 2e-13 from SymPy's,
///   which is the references' own agreement.
///
/// Tolerance: 1e-12 relative — five times the references' disagreement, and above the code's
/// rounding estimate `2^(2m+2n−2) ε / R` (3.8e-14 relative at the worst point here). Every point
/// is asserted to be in the branch. Both orders of the arguments are checked.
#[test]
fn the_far_unequal_branch_against_independent_references() {
    let points: [(u32, f64, u32, f64, f64, f64); 6] = [
        (1, 10.0, 1, 0.5, 6.0, 0.165_007_952_817_092_9),
        (1, 12.0, 1, 1.0, 5.0, 0.199_944_884_238_381),
        (1, 9.0, 1, 0.6, 7.0, 0.142_688_893_150_248_8),
        (1, 10.0, 5, 0.5, 6.0, 0.098_705_349_023_525_28),
        (2, 9.0, 3, 0.6, 7.0, 0.136_129_175_042_374_27),
        (4, 8.0, 2, 0.7, 8.0, 0.124_818_479_008_142_74),
    ];
    for (na, za, nb, zb, r, reference) in points {
        assert!((za - zb).abs() * r > 50.0, "not in the branch");
        for j in [
            coulomb_integral(na, za, nb, zb, r),
            coulomb_integral(nb, zb, na, za, r),
        ] {
            assert!(
                (j - reference).abs() <= 1e-12 * reference,
                "{na}s({za}) {nb}s({zb}) at {r}: {j} against {reference}"
            );
        }
    }
}

/// The upper end of eq 5's range, mirrored: ask C–F for a total of +5 and eq 18 with the total
/// puts carbon above +4, so it is fixed there and fluorine takes the rest, +1, exactly.
#[test]
fn a_charge_above_its_range_is_fixed_at_the_top() {
    let (c, f) = (Element::C, Element::F);
    let t = |e: Element| table_i(e);
    let r = 1.35 * ANGSTROM;
    let jcf = shielded_coulomb(c, slater_exponent(c), f, slater_exponent(f), r);
    let total = 5.0;
    let qc_free = (t(f).electronegativity - t(c).electronegativity
        + total * (t(f).idempotential - jcf))
        / (t(c).idempotential + t(f).idempotential - 2.0 * jcf);
    assert!(qc_free > 4.0, "the case needs C past +4: {qc_free}");
    let out = Qeq::default()
        .equilibrate(&[c, f], &[[0.0; 3], [r, 0.0, 0.0]], total)
        .expect("solves");
    assert_eq!(out.charges, vec![4.0, 1.0]);
    assert_eq!(out.at_bound, vec![0]);
}

/// Hydrogen's eq 5′ range, both ends, through the iteration: H–F with a total of +5 puts the
/// hydrogen at +1 (its orbital then ζ⁰ + 1, eq 20) and the fluorine at +4; H–N with −3.9 puts the
/// hydrogen at −1 and the nitrogen at −2.9. Each hydrogen is reported fixed, so each got there by
/// the clamp and not by the arithmetic.
#[test]
fn hydrogen_is_held_inside_minus_one_to_plus_one() {
    let (h, f, n) = (Element::H, Element::F, Element::N);
    let at = [[0.0; 3], [1.0 * ANGSTROM, 0.0, 0.0]];
    let up = Qeq::default()
        .equilibrate(&[h, f], &at, 5.0)
        .expect("solves");
    assert_eq!(up.at_bound, vec![0], "{:?}", up.charges);
    assert_eq!(up.charges[0], 1.0);
    assert!((up.charges[1] - 4.0).abs() <= 4.0 * EPS);
    let down = Qeq::default()
        .equilibrate(&[h, n], &at, -3.9)
        .expect("solves");
    assert_eq!(down.at_bound, vec![0], "{:?}", down.charges);
    assert_eq!(down.charges[0], -1.0);
    assert!((down.charges[1] + 2.9).abs() <= 4.0 * EPS);
}

/// Every error that can be reached is reached by an input that should reach it.
///
/// - **Infeasible**: aspirin asked for +100 e, beyond what its atoms' ranges sum to; and no atoms.
/// - **NonFinite**: a NaN total. Without the check it passed every comparison — `max` discards a
///   NaN and `<=` is false of one — and came back `Ok` with NaN charges.
/// - **Singular**: two hydrogens at the distance R* where their shielded J equals J⁰_HH, found by
///   bisection — 0.7915 Å, longer than H₂'s 0.741 Å bond. There eq 6′'s quadratic form for moving charge from one to the other,
///   `J_AA + J_BB − 2 J_AB`, is zero, and the 2 × 2 system has no unique solution. Closer than R*
///   the form is negative — QEq's energy has no minimum for that pair — which is a property of the
///   model (J⁰_HH = 13.89 eV is less than the 1s self-repulsion, 18.19 eV), not of this code.
#[test]
fn the_errors_are_reachable() {
    let (_, elements, at) = aspirin();
    let lowest: f64 = elements.iter().map(|&e| charge_range(e).0).sum();
    let highest: f64 = elements.iter().map(|&e| charge_range(e).1).sum();
    assert_eq!(
        Qeq::default().equilibrate(&elements, &at, 100.0),
        Err(qeq::QeqError::Infeasible {
            total_charge: 100.0,
            lowest,
            highest
        })
    );
    assert!(matches!(
        Qeq::default().equilibrate(&[], &[], 0.0),
        Err(qeq::QeqError::Infeasible { .. })
    ));

    let co = [Element::C, Element::O];
    let pair = [[0.0; 3], [1.128 * ANGSTROM, 0.0, 0.0]];
    assert_eq!(
        Qeq::default().equilibrate(&co, &pair, f64::NAN),
        Err(qeq::QeqError::NonFinite)
    );
    let ch = [Element::C, Element::H];
    assert_eq!(
        Qeq::default().equilibrate(&ch, &pair, f64::NAN),
        Err(qeq::QeqError::NonFinite)
    );

    let h = Element::H;
    let z = slater_exponent(h);
    let j0 = table_i(h).idempotential;
    let jhh = |r: f64| shielded_coulomb(h, z, h, z, r);
    let (mut lo, mut hi) = (0.0, 2.0 * ANGSTROM);
    assert!(jhh(lo) > j0 && jhh(hi) < j0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if jhh(mid) > j0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    println!(
        "R* = {:.6} Å, where J_HH = J⁰_HH = {j0} eV (J_HH(0) = {:.4} eV)",
        hi / ANGSTROM,
        jhh(0.0)
    );
    let out = Qeq::default().equilibrate(&[h, h], &[[0.0; 3], [hi, 0.0, 0.0]], 0.0);
    assert_eq!(out, Err(qeq::QeqError::Singular), "{out:?}");
    // Away from R*, the same pair solves: Q = 0 each, by symmetry.
    let fine = Qeq::default()
        .equilibrate(&[h, h], &[[0.0; 3], [1.5 * ANGSTROM, 0.0, 0.0]], 0.0)
        .expect("H₂ stretched");
    assert_eq!(fine.charges, vec![0.0, 0.0]);
}

/// Eq 18 (p. 3360): for two atoms with no hydrogen, one solve, and
/// `Q_A = (χ⁰_B − χ⁰_A) / (J_AA + J_BB − 2 J_AB)` exactly — here carbon monoxide's C and O at
/// 1.128 Å and C–F at 1.35 Å. Tolerance: the solve is a 2 × 2 elimination, a few roundings of
/// numbers of size ≤ 15 eV; 32 ε.
#[test]
fn two_atoms_without_hydrogen_are_eq_18() {
    for (b, r) in [(Element::O, 1.128), (Element::F, 1.35), (Element::Cl, 1.75)] {
        let elements = [Element::C, b];
        let at = [[0.0; 3], [r * ANGSTROM, 0.0, 0.0]];
        let out = Qeq::default()
            .equilibrate(&elements, &at, 0.0)
            .expect("solves");
        let (ta, tb) = (table_i(Element::C), table_i(b));
        let jab = shielded_coulomb(
            Element::C,
            slater_exponent(Element::C),
            b,
            slater_exponent(b),
            r * ANGSTROM,
        );
        let q = (tb.electronegativity - ta.electronegativity)
            / (ta.idempotential + tb.idempotential - 2.0 * jab);
        assert_eq!(out.iterations, 1);
        assert!(
            (out.charges[0] - q).abs() <= 32.0 * EPS,
            "{b}: {} vs {q}",
            out.charges[0]
        );
        assert!((out.charges[0] + out.charges[1]).abs() <= 2.0 * EPS, "{b}");
        assert!(out.at_bound.is_empty());
        // Eq 8: the common chemical potential is either atom's eq 7.
        let mu = ta.electronegativity + ta.idempotential * out.charges[0] + jab * out.charges[1];
        assert!((out.electronegativity - mu).abs() <= 64.0 * EPS * 15.0);
    }
}

/// Eqs 5 and 13: a charge pushed past its range is fixed at the boundary, and the rest solve with
/// its field folded into their electronegativity. Fluorine's range is −1 to +7; ask C–F for a
/// total of −3 and eq 18 with the total puts fluorine below −1, so it is fixed there and carbon
/// takes the rest, −2, exactly. With a third atom, the two free ones solve eq 18 with
/// `χ⁰F = χ⁰ − J_·F` (eq 13) and total −2.
#[test]
fn a_charge_outside_its_range_is_fixed_at_the_boundary() {
    let (c, f, o) = (Element::C, Element::F, Element::O);
    let x = |e: Element| table_i(e).electronegativity;
    let jj = |e: Element| table_i(e).idempotential;
    let z = slater_exponent;
    // Two atoms: the unclamped solution first, to show it is out of range.
    let r = 1.35 * ANGSTROM;
    let jcf = shielded_coulomb(c, z(c), f, z(f), r);
    let total = -3.0;
    // Q_C + Q_F = total, χ_C = χ_F.
    let qf_free = (x(c) - x(f) + total * (jj(c) - jcf)) / (jj(c) + jj(f) - 2.0 * jcf);
    assert!(qf_free < -1.0, "the case needs F past its range: {qf_free}");
    let out = Qeq::default()
        .equilibrate(&[c, f], &[[0.0; 3], [r, 0.0, 0.0]], total)
        .expect("solves");
    assert_eq!(out.charges, vec![-2.0, -1.0]);
    assert_eq!(out.at_bound, vec![1]);

    // Three atoms, F far enough out that it is the one fixed.
    let at = [
        [0.0; 3],
        [1.35 * ANGSTROM, 0.0, 0.0],
        [-1.2 * ANGSTROM, 0.0, 0.0],
    ];
    let out = Qeq::default()
        .equilibrate(&[c, f, o], &at, total)
        .expect("solves");
    assert_eq!(out.at_bound, vec![1], "{:?}", out.charges);
    assert_eq!(out.charges[1], -1.0);
    let d = |i: usize, k: usize| -> f64 {
        (0..3)
            .map(|a| (at[i][a] - at[k][a]).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    let (jcf, jof, jco) = (
        shielded_coulomb(c, z(c), f, z(f), d(0, 1)),
        shielded_coulomb(o, z(o), f, z(f), d(2, 1)),
        shielded_coulomb(c, z(c), o, z(o), d(0, 2)),
    );
    let (xc, xo) = (x(c) - jcf, x(o) - jof);
    let rest = total + 1.0;
    let qc = (xo - xc + rest * (jj(o) - jco)) / (jj(c) + jj(o) - 2.0 * jco);
    assert!(
        (out.charges[0] - qc).abs() <= 64.0 * EPS * 4.0,
        "{} vs {qc}",
        out.charges[0]
    );
    assert!((out.charges[0] + out.charges[2] - rest).abs() <= 8.0 * EPS * 4.0);
}

/// The total charge is the one asked for — the first row of eq 10 is that sum, so it holds to
/// the rounding of adding up the charges: `n ε Σ|Q|`. Aspirin neutral, as a cation and as an
/// anion; and an entry's total defaults to the sum of its formal charges (ammonium, +1).
#[test]
fn the_total_charge_is_conserved() {
    let (_, elements, at) = aspirin();
    for total in [0.0, 1.0, -1.0] {
        let q = Qeq::default()
            .equilibrate(&elements, &at, total)
            .expect("solves")
            .charges;
        let sum: f64 = q.iter().sum();
        let bound = q.len() as f64 * EPS * q.iter().map(|x| x.abs()).sum::<f64>();
        assert!((sum - total).abs() <= bound, "total {total}: {sum}");
    }
    let nh4 = "data_NH4\n_chem_comp.id NH4\n_chem_comp.pdbx_formal_charge 1\nloop_\n\
               _chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n\
               _chem_comp_atom.charge\n_chem_comp_atom.pdbx_aromatic_flag\n\
               _chem_comp_atom.pdbx_model_Cartn_x_ideal\n_chem_comp_atom.pdbx_model_Cartn_y_ideal\n\
               _chem_comp_atom.pdbx_model_Cartn_z_ideal\n\
               NH4 N N 1 N 0 0 0\nNH4 H1 H 0 N 0.593 0.593 0.593\nNH4 H2 H 0 N 0.593 -0.593 -0.593\n\
               NH4 H3 H 0 N -0.593 0.593 -0.593\nNH4 H4 H 0 N -0.593 -0.593 0.593\nloop_\n\
               _chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n\
               _chem_comp_bond.pdbx_aromatic_flag\nN H1 SING N\nN H2 SING N\nN H3 SING N\nN H4 SING N\n";
    let c = Component::from_ccd(nh4).expect("NH4 parses");
    let q = qeq::charges(&c, &positions(&c)).expect("solves").charges;
    let sum: f64 = q.iter().sum();
    assert!((sum - 1.0).abs() <= 5.0 * EPS * 2.0, "{sum}");
    // Tetrahedral: four equal hydrogens.
    for h in 2..5 {
        assert!((q[h] - q[1]).abs() <= 1e-14, "{q:?}");
    }
}

fn rotate(p: [f64; 3]) -> [f64; 3] {
    // Rotation by 0.7 rad about (1, 2, 2)/3, Rodrigues.
    let k = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0];
    let (c, s) = (0.7f64.cos(), 0.7f64.sin());
    let dot = k[0] * p[0] + k[1] * p[1] + k[2] * p[2];
    let cross = [
        k[1] * p[2] - k[2] * p[1],
        k[2] * p[0] - k[0] * p[2],
        k[0] * p[1] - k[1] * p[0],
    ];
    [0, 1, 2].map(|i| p[i] * c + cross[i] * s + k[i] * dot * (1.0 - c))
}

/// Charges depend only on distances, so translating and rotating the molecule leaves them where
/// they were — to what the distances' rounding moves them, and the iteration's tolerance: run at
/// 10⁻¹³ e so that each run is within `10⁻¹³ ρ/(1−ρ)` of the fixed point, ρ the contraction per
/// iteration (measured below 0.5), and a moved distance is off by ~ε·|r| ~ 1e-16 relative, worth
/// ~1e-15 e at `|∂Q/∂R|` ~ 1 e/Å. Bound: 1e-12 e.
#[test]
fn the_charges_do_not_see_where_the_molecule_is_or_which_way_it_faces() {
    let (_, elements, at) = aspirin();
    let qeq = Qeq {
        tolerance: 1e-13,
        ..Qeq::default()
    };
    let base = qeq.equilibrate(&elements, &at, 0.0).expect("solves");
    let moved: Vec<[f64; 3]> = at
        .iter()
        .map(|p| {
            let r = rotate(*p);
            [r[0] + 3e-10, r[1] - 2e-10, r[2] + 5e-10]
        })
        .collect();
    let other = qeq.equilibrate(&elements, &moved, 0.0).expect("solves");
    let worst = base
        .charges
        .iter()
        .zip(&other.charges)
        .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
    println!("largest change under a rigid motion: {worst:e} e");
    assert!(worst <= 1e-12, "{worst:e}");
}

/// Eq 8 at the solution: every atom's chemical potential (eq 7) is the same, computed with the
/// final charges in a fresh matrix. The matrix the last solve used was built from the previous
/// iterate, at most the tolerance away, and moving a hydrogen's charge by δ moves its J's by at
/// most `J⁰_HH δ/ζ⁰` on the diagonal and `|∂J/∂ζ| δ` ≲ 20 eV δ off it — so with charges below 1 e
/// and n atoms the spread is at most `2 · 20 eV · n · tolerance`, plus rounding.
#[test]
fn every_atom_has_the_same_electronegativity() {
    let (_, elements, at) = aspirin();
    let qeq = Qeq {
        tolerance: 1e-12,
        ..Qeq::default()
    };
    let out = qeq.equilibrate(&elements, &at, 0.0).expect("solves");
    let n = elements.len();
    let j = coulomb_matrix(&elements, &at, &out.charges);
    let chi: Vec<f64> = (0..n)
        .map(|a| {
            table_i(elements[a]).electronegativity
                + (0..n).map(|b| j[a * n + b] * out.charges[b]).sum::<f64>()
        })
        .collect();
    let spread = chi
        .iter()
        .fold(0.0f64, |m, x| m.max((x - out.electronegativity).abs()));
    let bound = 2.0 * 20.0 * n as f64 * qeq.tolerance + 1e-12;
    println!(
        "μ = {} eV, spread {spread:e} eV (bound {bound:e})",
        out.electronegativity
    );
    assert!(spread <= bound);
}

/// The hydrogen iteration from zero (p. 3361: "six to ten iterations for an initial guess of
/// zero"), reported at the tolerances that number could mean, and to this module's default.
#[test]
fn how_many_iterations_the_hydrogen_iteration_takes() {
    let (_, elements, at) = aspirin();
    for tolerance in [1e-2, 1e-3, 1e-4, 1e-10] {
        let qeq = Qeq {
            tolerance,
            ..Qeq::default()
        };
        let out = qeq.equilibrate(&elements, &at, 0.0).expect("solves");
        println!(
            "aspirin, tolerance {tolerance:e} e: {} solves",
            out.iterations
        );
    }
    let out = Qeq {
        max_iterations: 3,
        ..Qeq::default()
    }
    .equilibrate(&elements, &at, 0.0);
    assert!(
        matches!(out, Err(qeq::QeqError::NotConverged { iterations: 3, .. })),
        "{out:?}"
    );
}

/// **Aspirin's QEq charges at its relaxed geometry, and UFF's electrostatic energy with them
/// there** — reported, the charges fixed at that geometry. Nothing about aspirin's charges has a
/// closed form; what is asserted is what must hold: neutral, every charge inside its range, the
/// carbonyl and hydroxyl oxygens negative, the acid hydrogen the most positive hydrogen, and the
/// energy the sum of eq 43 over the included pairs.
#[test]
fn aspirin_at_its_minimum() {
    let (c, _, start) = aspirin();
    let types = uff::assign(&c);
    let ff = ForceField::new(&c, &types).expect("supported");
    let (_, at) = common::relaxed(&ff, &start);
    let out = qeq::charges(&c, &at).expect("aspirin equilibrates");
    println!("| atom | QEq charge (e) |");
    println!("| --- | --- |");
    for (a, q) in c.atoms().iter().zip(&out.charges) {
        println!("| {} | {q:+.4} |", a.name);
    }
    println!("{} solves", out.iterations);
    let sum: f64 = out.charges.iter().sum();
    assert!(sum.abs() < 1e-13);
    for (a, q) in c.atoms().iter().zip(&out.charges) {
        let (lo, hi) = charge_range(a.element);
        assert!(*q > lo && *q < hi, "{}", a.name);
    }
    for o in ["O1", "O2", "O3", "O4"] {
        assert!(out.charges[index(&c, o)] < 0.0, "{o}");
    }
    let ho1 = out.charges[index(&c, "HO1")];
    for (a, q) in c.atoms().iter().zip(&out.charges) {
        if a.element == Element::H && a.name != "HO1" {
            assert!(*q < ho1, "{} {q} vs HO1 {ho1}", a.name);
        }
    }
    let charged = ff.clone().with_qeq_charges(&at).expect("equilibrates");
    assert_eq!(charged.charges(), &out.charges[..]);
    let e = charged.energy(&at);
    let mut by_hand = 0.0;
    for p in charged.pairs() {
        let [i, j] = p.atoms;
        let d: f64 = (0..3)
            .map(|k| (at[i][k] - at[j][k]).powi(2))
            .sum::<f64>()
            .sqrt();
        by_hand += pantometry_forcefield::energy::coulomb(out.charges[i], out.charges[j], d);
    }
    println!(
        "electrostatic energy at the relaxed geometry with these charges: {:.4} kcal/mol \
         (total {:.4}, of which without charges {:.4})",
        e.electrostatic / KCAL_PER_MOL,
        e.total / KCAL_PER_MOL,
        ff.energy(&at).total / KCAL_PER_MOL
    );
    assert!((e.electrostatic - by_hand).abs() <= 1e-12 * by_hand.abs());
    // The other five terms do not see the charges.
    let zero = ff.energy(&at);
    assert_eq!(e.bond, zero.bond);
    assert_eq!(e.van_der_waals, zero.van_der_waals);
}
