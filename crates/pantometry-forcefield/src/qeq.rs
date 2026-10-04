//! Charge equilibration (QEq): partial charges from geometry, electronegativity and atomic size.
//!
//! > A. K. Rappé and W. A. Goddard III, "Charge equilibration for molecular dynamics
//! > simulations", *J. Phys. Chem.* **95**, 3358–3363 (1991).
//! > [doi:10.1021/j100161a070](https://doi.org/10.1021/j100161a070)
//!
//! UFF prescribes these charges for its electrostatic term (p. 10031 of the UFF paper). Page and
//! equation numbers below are the QEq paper's, read off the scanned pages.
//!
//! # The method
//!
//! The energy of atom A, expanded to second order in its charge (eq 1′, p. 3359), is
//! `E_A(Q) = E_A0 + χ⁰_A Q_A + ½ J⁰_AA Q_A²`, with χ⁰ the electronegativity and J⁰ the
//! idempotential, both from Table I. The molecule's electrostatic energy (eq 6′, p. 3359) is
//!
//! `E(Q₁…Q_N) = Σ_A (E_A0 + χ⁰_A Q_A) + ½ Σ_A Σ_B Q_A Q_B J_AB`,
//!
//! with `J_AA` the idempotential J⁰_AA (the sentence after eq 6′ says `J_AA(R) → J⁰_AA as R → 0`)
//! and `J_AB` the shielded Coulomb interaction below. Its derivative is each atom's chemical
//! potential (eq 7, p. 3359), `χ_A = χ⁰_A + Σ_B J_AB Q_B`, and equilibrium sets them all equal
//! (eq 8), `χ₁ = χ₂ = … = χ_N`, which with the total charge (eq 9), `Q_tot = Σ Q_i`, is N linear
//! equations in N charges (eq 10, `C Q = −D`).
//!
//! # A misprint in eq 12, corrected
//!
//! **Eq 12 prints `C_1j = Q_j`. It is `C_1j = 1`.** Row 1 of `C Q = −D` is the total-charge
//! condition, eq 9, and eq 11 gives it `D₁ = −Q_tot`, so it reads `Σ_j C_1j Q_j = Q_tot`; that is
//! eq 9 only if every `C_1j` is one. A row of the unknowns would make the system quadratic. The
//! other rows are as printed and follow from eqs 7 and 8: `χ_i = χ₁` for `i ≥ 2` is
//! `Σ_j (J_ij − J_1j) Q_j = −(χ⁰_i − χ⁰_1)`, which is eq 12's `C_ij = J_ij − J_1j` and eq 11's
//! `D_i = χ⁰_i − χ⁰_1`.
//!
//! # The allowed range of a charge (eqs 5, 5′ and 13)
//!
//! A charge is held inside the range the valence shell allows: from the shell full to the shell
//! empty. Eq 5 (p. 3359) gives three, `−7 < Q_Li < +1`, `−4 < Q_C < +4`, `−2 < Q_O < +6`, "etc.";
//! the rule they share — `−(8 − v)` to `+v` for `v` valence electrons in an s and p shell — is
//! extended here to the other seven elements ([`charge_range`]), and hydrogen's is eq 5′,
//! `−1 < Q_H < +1` (p. 3361). The procedure is the paper's (p. 3359): solve; fix every atom
//! outside its range at the boundary it crossed; solve the reduced system for the rest, with each
//! fixed atom's field folded into the others' electronegativity, `χ⁰F_A = χ⁰_A + Σ_{B fixed}
//! J_AB Q_B` (eq 13), and the fixed charge taken out of the total; repeat until no free atom is
//! out of range. An atom once fixed stays fixed, as the paper describes it.
//!
//! # The shielded Coulomb interaction (eqs 14–17, p. 3359–3360)
//!
//! `J_AB(R)` is the Coulomb integral between two normalised `ns` Slater densities (eq 15,
//! `φ = N r^(n−1) e^(−ζ r)`), one on each atom: `14.4/R` (eq 14) at large R, and finite as R → 0.
//! [`coulomb_integral`] computes it **exactly**, as a finite sum of closed forms, for n = 1 to 5:
//! `n` is 1 for H, 2 for C, N, O and F, 3 for P, S and Cl, 4 for Br and 5 for I — each element's
//! valence shell ([`principal_quantum_number`]). How it is done is in its documentation. The
//! units are the paper's: J in eV, R converted to bohr with CODATA's `a₀`
//! ([`BOHR`]) and hartree to eV with CODATA's [`HARTREE_EV`], so `J → 14.39964/R` (Å), which is
//! the paper's 14.4 to the three figures it prints.
//!
//! **The diagonal is Table I's idempotential J⁰_AA, not the Slater self-integral.** The paper
//! says why (p. 3359): J⁰ = IP − EA is measured, and a single Slater function's self-repulsion is
//! not what an atom's IP and EA see. They differ: carbon's 2s at eq 17′'s ζ = 0.87150 has
//! self-repulsion `(93/256) ζ` = 0.31660 hartree = 8.615 eV, against J⁰_C = 10.126 eV.
//!
//! # Which ζ: eq 17 with λ = ½, not Table I's printed column
//!
//! The orbital exponent comes from the covalent radius R_A of Table I (eq 17, p. 3360):
//! `ζ_A = λ (2n + 1) / (2 R_A)`, R_A in bohr. Section III fits λ to the dipole moments of twelve
//! alkali halides and finds 0.4913 (p. 3360), then says "rounding off to λ = ½ leads also to an
//! average error of 0.0018 e, and hence (17) becomes" eq 17′, `ζ_A = (2n + 1) / (4 R_A)`.
//! **Table I's printed ζ is λ = 0.4913**: for C, N, F, P, S, Cl, Br and I the λ each printed ζ
//! implies is 0.49123–0.49129. Hydrogen's printed 1.0698 is λ = ½ (0.50002), as the text below
//! eq 20 says it is. Oxygen's 0.9745 is neither: it implies λ = 0.49280, and is what λ = 0.4913
//! gives for R_O = 0.667 Å rather than the printed 0.669.
//!
//! **This module uses eq 17′, λ = ½, for every element**, because that is what reproduces the
//! paper's own charges and the printed column does not: hydrogen's charge in HF is 0.4623 with
//! λ = ½ against the paper's 0.462 (Table III, p. 3361), and 0.4568 with the printed ζ — which
//! misses by more than the printed figure allows; H₂O 0.3531 against 0.353 (0.3453), NH₃ 0.2412
//! against 0.243 (0.2295), CH₄ 0.1497 against 0.149 (0.1342). See
//! `tests/the_charge_equilibration_paper.rs` for every molecule compared. [`table_i`] carries the
//! printed ζ, as printed, so the comparison can be rerun.
//!
//! # Hydrogen (eqs 20–23, p. 3360–3361)
//!
//! Hydrogen's orbital is allowed to change size with its charge: `ζ_H(Q_H) = ζ⁰_H + Q_H`
//! (eq 20), with ζ⁰_H = 1.0698 bohr⁻¹ from eq 17′ and R_H = 0.371 Å, and so its idempotential
//! `J_HH(Q) = (1 + Q_H/ζ⁰_H) J⁰_HH` (eq 21). The changed ζ enters **every** J involving that
//! hydrogen, not only its own: eq 20 changes the orbital, and eq 14's integral is over the
//! orbital. The equations are then no longer linear, and the paper solves them by iteration
//! (p. 3361): "we use (20) with an estimated Q_H and iterate until all Q_H's are
//! self-consistent", from an initial guess of zero, in six to ten iterations.
//!
//! **That iteration is not the stationary point of eq 23.** Eq 23 writes hydrogen's energy as
//! `E_H0 + χ⁰_H Q_H + ½ J⁰_HH Q_H² (1 + Q_H/1.0698)`; its derivative is
//! `χ⁰_H + J⁰_HH Q_H (1 + 1.5 Q_H/ζ⁰_H)`, while the iteration's chemical potential is eq 7 with
//! eq 21's J, `χ⁰_H + J⁰_HH Q_H (1 + Q_H/ζ⁰_H) + …`. Only the iteration reproduces the paper's
//! numbers: with eq 23's derivative on the diagonal, hydrogen's charge in HF is 0.405, not 0.462.
//! **This module does what the paper did**, the iteration, and does not minimise eq 23.
//!
//! The iteration here: start every charge at zero; build every J at the current hydrogen
//! charges; solve (with the range clamping); repeat until no charge changes by more than
//! [`Qeq::tolerance`] (10⁻¹⁰ e by default) between two solves, or fail with
//! [`QeqError::NotConverged`] after [`Qeq::max_iterations`]. Undamped, as the paper describes.
//! Without hydrogen, the equations are linear and one solve is the answer.
//!
//! # What the paper's numbers cannot see
//!
//! Two decimals of a charge do not settle every choice. **The principal quantum number** is the
//! text's — "an atom whose outer valence orbital is ns, np, or nd" gets an ns Slater orbital
//! (p. 3359) — and Table IV cannot see it: chlorine with n = 2 gives HCl's hydrogen 0.3203, nearer
//! the printed 0.32 than n = 3's 0.3173. **Table II can**: its alkali halides are printed to three
//! decimals, and with n = 2 for chlorine NaCl, KCl and RbCl miss by 0.0011–0.0013 e (λ = ½
//! column) against a tolerance of 0.0005. Table II is also the only physical check of the 4s and
//! 5s integrals — K and Br are 4s, Rb and I 5s — and all eighteen of its charges for Na, K and Rb
//! with Cl, Br and I are reproduced to 0.00048 e (`tests/the_charge_equilibration_paper.rs`). **The ranges for N, F, P, S, Cl, Br and I** are the extension stated above; no
//! molecule compared reaches one. **Eq 21's ζ⁰_H** is eq 17′'s 1.06977 here, the printed 1.0698
//! to its four figures.
//!
//! # Fixed per geometry
//!
//! **The charges are fixed at the geometry they were computed for.** [`ForceField::with_qeq_charges`]
//! computes them once, at the positions given, and the force field then treats them as constants:
//! its electrostatic force is the derivative at fixed charges, which is the paper's usage ("solved
//! once for a given structure", p. 3359) and not the derivative of the equilibrated energy, whose
//! charges move with the atoms. Recompute them at a new geometry if the geometry has moved far.
//!
//! # Determinism
//!
//! The linear solve is Gaussian elimination with partial pivoting in a fixed order, with no
//! external crate. The integrals call the platform's `exp`, which is not correctly rounded, so a
//! charge repeats bit for bit on one machine and not necessarily across platforms — the same
//! promise as the rest of this crate.
//!
//! [`ForceField::with_qeq_charges`]: crate::energy::ForceField::with_qeq_charges

use crate::ccd::{Component, Element};
use std::fmt;

/// The Bohr radius a₀, in metres (CODATA 2018). The paper converts with 0.52917 Å (p. 3360),
/// 1.4 × 10⁻⁵ smaller in relative terms, which moves no ζ in its fourth figure.
pub const BOHR: f64 = 0.529_177_210_903e-10;

/// One hartree, in electronvolts (CODATA 2018).
pub const HARTREE_EV: f64 = 27.211_386_245_988;

/// λ of eq 17′ (p. 3360): the paper's fitted 0.4913, rounded to ½ as the paper rounds it.
pub const LAMBDA: f64 = 0.5;

/// Up to this `|ζ_A − ζ_B| R`, [`coulomb_integral`] expands its shielding term about the origin
/// of its integration rectangle; above it, about the corner where the weight is. Chosen by
/// measuring both in f64 against 60 digits; see [`coulomb_integral`].
pub const ORIGIN_EXPANSION_UP_TO: f64 = 2.0;

/// The largest principal quantum number [`coulomb_integral`] accepts: iodine's 5s.
pub const MAX_N: u32 = 5;

/// One row of the paper's Table I (p. 3359), as printed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Row {
    /// χ⁰, the electronegativity, eV.
    pub electronegativity: f64,
    /// J⁰, the idempotential, eV.
    pub idempotential: f64,
    /// R, the covalent radius, Å.
    pub radius: f64,
    /// ζ, the Slater exponent, bohr⁻¹, **as printed** — which is eq 17 with λ = 0.4913 and not
    /// what this module uses. See the module documentation and [`slater_exponent`].
    pub printed_zeta: f64,
}

/// `element`'s row of Table I (p. 3359), as printed.
///
/// Hydrogen's χ and J are the paper's fitted 4.5280 and 13.8904 eV (eq 22, p. 3361; Table I's
/// footnote b), "values for Q_H = 0".
pub fn table_i(element: Element) -> Row {
    let (electronegativity, idempotential, radius, printed_zeta) = match element {
        Element::H => (4.5280, 13.8904, 0.371, 1.0698),
        Element::C => (5.343, 10.126, 0.759, 0.8563),
        Element::N => (6.899, 11.760, 0.715, 0.9089),
        Element::O => (8.741, 13.364, 0.669, 0.9745),
        Element::F => (10.874, 14.948, 0.706, 0.9206),
        Element::P => (5.463, 8.000, 1.102, 0.8257),
        Element::S => (6.928, 8.972, 1.047, 0.8690),
        Element::Cl => (8.564, 9.892, 0.994, 0.9154),
        Element::Br => (7.790, 8.850, 1.141, 1.0253),
        Element::I => (6.822, 7.524, 1.333, 1.0726),
    };
    Row {
        electronegativity,
        idempotential,
        radius,
        printed_zeta,
    }
}

/// The principal quantum number of `element`'s valence shell, the `n` of its `ns` Slater
/// orbital (eq 15): 1 for H, 2 for C–F, 3 for P–Cl, 4 for Br, 5 for I.
pub fn principal_quantum_number(element: Element) -> u32 {
    element.period()
}

/// ζ⁰ of `element`, bohr⁻¹, by eq 17′ (p. 3360): `(2n + 1) / (4 R)`, R Table I's radius in
/// bohr. For hydrogen this is ζ⁰_H, 1.06977, the paper's 1.0698; eq 20 then adds Q_H.
pub fn slater_exponent(element: Element) -> f64 {
    let n = f64::from(principal_quantum_number(element));
    let r_bohr = table_i(element).radius * 1e-10 / BOHR;
    LAMBDA * (2.0 * n + 1.0) / (2.0 * r_bohr)
}

/// The range eq 5 allows `element`'s charge, `(lowest, highest)` in elementary charges: from the
/// valence shell full to empty. Hydrogen's is eq 5′, (−1, +1); the others are `(−(8 − v), +v)`
/// for `v` valence s and p electrons — the rule eq 5's three printed ranges (Li, C, O) share,
/// extended to N, F, P, S, Cl, Br and I.
pub fn charge_range(element: Element) -> (f64, f64) {
    let v = match element {
        Element::H => return (-1.0, 1.0),
        Element::C => 4.0,
        Element::N | Element::P => 5.0,
        Element::O | Element::S => 6.0,
        Element::F | Element::Cl | Element::Br | Element::I => 7.0,
    };
    (-(8.0 - v), v)
}

/// `k!` as a float, exact for every `k` used here (at most 20).
fn factorial(k: u32) -> f64 {
    (1..=k).fold(1.0, |f, i| f * f64::from(i))
}

/// `e^{−x} Σ_{k<s} x^k/k!`, the regularised upper incomplete gamma function Q(s, x) for an
/// integer `s ≥ 1`, as a finite sum of positive terms.
fn upper_sum(s: u32, x: f64) -> f64 {
    let (mut term, mut sum) = (1.0, 1.0);
    for k in 1..s {
        term *= x / f64::from(k);
        sum += term;
    }
    sum * (-x).exp()
}

/// The regularised lower incomplete gamma function P(s, x) by its series,
/// `x^s e^{−x}/s! Σ_k x^k / ((s+1)…(s+k))`: positive terms, each smaller than the last when
/// `x < s + 1`, which is where it is used.
fn lower_series(s: u32, x: f64) -> f64 {
    let mut lead = (-x).exp();
    for k in 1..=s {
        lead *= x / f64::from(k);
    }
    let (mut term, mut sum) = (1.0, 1.0);
    let mut k = 1.0;
    loop {
        term *= x / (f64::from(s) + k);
        sum += term;
        if term <= sum * 1e-17 {
            break;
        }
        k += 1.0;
    }
    lead * sum
}

/// `(P(s, x), Q(s, x))`, each computed from whichever form has no cancellation: the series for P
/// below `x = s`, where P is the smaller, and the finite sum for Q above, where Q is.
fn incomplete_gamma(s: u32, x: f64) -> (f64, f64) {
    if x < f64::from(s) {
        let p = lower_series(s, x);
        (p, 1.0 - p)
    } else {
        let q = upper_sum(s, x);
        (1.0 - q, q)
    }
}

/// The Coulomb integral between two normalised Slater densities, in **hartree**, at a distance
/// `r` **bohr** apart: `∫∫ ρ_A(r₁) ρ_B(r₂) / |r₁ − r₂|`, with `ρ = |φ|²` and
/// `φ = N r^(n−1) e^(−ζ r)` (eq 15, p. 3359) — principal quantum numbers `n_a`, `n_b` and
/// exponents `zeta_a`, `zeta_b` (bohr⁻¹). Exact, as a finite sum of closed forms, for every
/// `n` from 1 to [`MAX_N`] and any two exponents, equal or not.
///
/// # How
///
/// B's radial density is a gamma distribution, `β^(2n+1) s^(2n) e^(−βs)/(2n)!` with β = 2ζ_B, so
/// its potential is closed: `s V_B(s) = 1 − e^(−βs) Σ_{k<2n} b_k s^k` with
/// `b_k = β^k (2n − k) / (2n k!)`. Averaged over the sphere of radius `r` about A,
/// `⟨V_B⟩ = (1/(2rR)) ∫_{|r−R|}^{r+R} s V_B(s) ds`; the `1` gives `1/max(r, R)`, whose integral
/// against A's density is A's own potential at B, `V_A(R) = P(2m+1, αR)/R + (α/2m) Q(2m, αR)`
/// (α = 2ζ_A, m = n_a, P and Q the regularised incomplete gamma functions). The exponential
/// part is a double integral of polynomials times exponentials over the triangle
/// `|r − t| ≤ R ≤ r + t`, which `u = r + t`, `v = r − t` turn into a rectangle, `u ≥ R`,
/// `|v| ≤ R`, so it separates: `∫_R^∞ u^p e^(−σu) du` (σ = ζ_A + ζ_B, a finite sum) times
/// `∫_{−R}^{R} v^q e^(−δv) dv` (δ = ζ_A − ζ_B). The second is a power series in `δR` whose
/// terms all have one sign — so **equal or nearly equal exponents lose nothing**, where the
/// textbook closed form divides by `ζ_A − ζ_B`.
///
/// **Two expansions, and where each is used.** The polynomial `r^(2m−1) t^k` can be expanded
/// about the origin of `(u, v)` or about the corner of the rectangle where the weight
/// `e^(−σu − δv)` is largest — `u = R`, `v = ∓R`. About the origin, the binomial terms alternate,
/// and when one density is compact and the other diffuse the weight sits at the corner, where
/// the compact atom's coordinate is near zero and terms of size `R^(2m−1)` cancel down to it: an
/// **earlier version of this function, which used the origin expansion everywhere, was wrong by
/// 6.7 × 10⁻¹¹ hartree** (5.4 × 10⁻¹⁰ relative; 3.7 × 10⁻⁷ of the shielding term) for 4s(ζ 8)
/// against 2s(ζ 0.7) at 8 bohr, against two independent 40-digit quadratures. About the corner,
/// every binomial term is positive except one inner sum over `(2R − w)^j`, which falls by about
/// `1/(2|δ|R)` a term — good when the exponents differ and poor when they do not, where it
/// costs up to `2^j`. Measured in f64 against the same formula at 60 digits over nineteen cases
/// spanning both regimes, the two are comparable at `|δ|R` = 2, and **the origin expansion is used
/// up to [`ORIGIN_EXPANSION_UP_TO`] = 2 and the corner expansion above**; the worst error either
/// then makes over those cases is 4.0 × 10⁻¹⁴ relative. The exponentials are combined as
/// `e^(−2 min(ζ) R)` before they are multiplied, so nothing overflows at any distance.
///
/// At `r = 0` exactly it is `∫ ρ_A V_B`, a finite sum, and the general form is not used: it is a
/// difference of terms of order `1/R`, which loses `ε/R` near zero.
///
/// **Rounding.** Measured, not bounded: 4.0 × 10⁻¹⁴ relative at the worst over the nineteen
/// cases above, and below 10⁻¹² relative against independent references in the branch where the
/// earlier version failed (`tests/the_charge_equilibration_against_closed_forms.rs`). Computing
/// the integral both ways round, which runs different arithmetic, disagrees by at most 31% of the
/// estimate `2^(2m+2n−2) ε / R` hartree over that test's grid, which asserts eight times it.
///
/// # Panics
///
/// If `n_a` or `n_b` is not 1 to [`MAX_N`], an exponent is not positive and finite, or `r` is
/// negative or not finite.
pub fn coulomb_integral(n_a: u32, zeta_a: f64, n_b: u32, zeta_b: f64, r: f64) -> f64 {
    assert!(
        (1..=MAX_N).contains(&n_a) && (1..=MAX_N).contains(&n_b),
        "principal quantum numbers 1 to {MAX_N}, not {n_a} and {n_b}"
    );
    assert!(
        zeta_a > 0.0 && zeta_a.is_finite() && zeta_b > 0.0 && zeta_b.is_finite(),
        "Slater exponents must be positive and finite, not {zeta_a} and {zeta_b}"
    );
    assert!(
        r >= 0.0 && r.is_finite(),
        "a distance must be non-negative and finite, not {r}"
    );
    let (alpha, beta) = (2.0 * zeta_a, 2.0 * zeta_b);
    let (m, n) = (n_a, n_b);
    // B's potential: s V_B(s) = 1 − e^{−βs} Σ_{k<2n} b_k s^k.
    let b: Vec<f64> = (0..2 * n)
        .map(|k| beta.powi(k as i32) * f64::from(2 * n - k) / (f64::from(2 * n) * factorial(k)))
        .collect();
    // A's radial density is `prefactor · r^{2m} e^{−αr}`.
    let prefactor = alpha.powi(2 * m as i32 + 1) / factorial(2 * m);
    if r == 0.0 {
        // ∫ ρ_A V_B = ∫ p_A(r)/r dr − ∫ p_A(r) e^{−βr} Σ b_k r^{k−1} dr.
        let s = alpha + beta;
        let shielding: f64 = b
            .iter()
            .enumerate()
            .map(|(k, bk)| {
                let k = k as u32;
                bk * factorial(2 * m - 1 + k) / s.powi((2 * m + k) as i32)
            })
            .sum();
        return alpha / f64::from(2 * m) - prefactor * shielding;
    }
    let (p, _) = incomplete_gamma(2 * m + 1, alpha * r);
    let (_, q) = incomplete_gamma(2 * m, alpha * r);
    let potential_of_a = p / r + alpha / f64::from(2 * m) * q;

    let sigma = zeta_a + zeta_b;
    let delta = zeta_a - zeta_b;
    let a = (2 * m - 1) as usize;
    let top = a + (2 * n - 1) as usize;
    let mut shielding = 0.0;
    if delta.abs() * r <= ORIGIN_EXPANSION_UP_TO {
        // About the origin of (u, v): u ≥ R, |v| ≤ R.
        // e^{σR} ∫_R^∞ u^p e^{−σu} du = p!/σ^{p+1} Σ_{k≤p} (σR)^k/k!.
        let u: Vec<f64> = (0..=top)
            .map(|p| {
                let (mut term, mut sum) = (1.0, 1.0);
                for k in 1..=p {
                    term *= sigma * r / k as f64;
                    sum += term;
                }
                factorial(p as u32) / sigma.powi(p as i32 + 1) * sum
            })
            .collect();
        // e^{−|δ|R} ∫_{−R}^{R} v^q e^{−δv} dv = e^{−|δ|R} 2 R^{q+1} Σ_{j ≡ q mod 2} (−δR)^j /
        // (j! (q + j + 1)): one sign throughout.
        let x = delta.abs() * r;
        let v: Vec<f64> = (0..=top)
            .map(|q| {
                let mut sum = 0.0;
                let mut power = 1.0; // (−δR)^j / j!
                let mut j = 0usize;
                loop {
                    if j > 0 {
                        power *= -delta * r / j as f64;
                    }
                    if (q + j) & 1 == 0 {
                        let term = power / (q + j + 1) as f64;
                        sum += term;
                        if j as f64 > x && term.abs() <= 1e-17 * sum.abs() {
                            break;
                        }
                    }
                    j += 1;
                }
                2.0 * r.powi(q as i32 + 1) * sum * (-x).exp()
            })
            .collect();
        // r = (u+v)/2, t = (u−v)/2: r^{2m−1} t^k = 2^{−(2m−1+k)} Σ C(2m−1,i) C(k,j) (−1)^j
        // u^{2m−1−i+k−j} v^{i+j}.
        for (k, bk) in b.iter().enumerate() {
            let mut inner = 0.0;
            for i in 0..=a {
                for j in 0..=k {
                    let c = binomial(a, i) * binomial(k, j) * if j & 1 == 0 { 1.0 } else { -1.0 };
                    inner += c * u[a - i + k - j] * v[i + j];
                }
            }
            shielding += bk * inner / 2f64.powi((a + k) as i32);
        }
    } else {
        // About the corner where the weight e^{−σu − δv} is largest: u = R and v = −R for δ > 0
        // (A the more compact), v = +R for δ < 0. With u′ = u − R ≥ 0 and w = R ± v ∈ [0, 2R],
        // the compact atom's coordinate is (u′ + w)/2 and the other's (u′ + 2R − w)/2.
        let d = delta.abs();
        let x = 2.0 * d * r;
        // ∫_0^∞ u^p e^{−σu} du.
        let u: Vec<f64> = (0..=top)
            .map(|p| factorial(p as u32) / sigma.powi(p as i32 + 1))
            .collect();
        // ∫_0^{2R} w^q e^{−|δ|w} dw = q!/|δ|^{q+1} P(q+1, 2|δ|R): positive terms either way.
        let w: Vec<f64> = (0..=top)
            .map(|q| {
                let (p, _) = incomplete_gamma(q as u32 + 1, x);
                factorial(q as u32) / d.powi(q as i32 + 1) * p
            })
            .collect();
        for (k, bk) in b.iter().enumerate() {
            // Powers of (u′ + w) and (u′ + 2R − w): r^{2m−1} and t^k, compact one first.
            let (near, far) = if delta >= 0.0 { (a, k) } else { (k, a) };
            let mut inner = 0.0;
            for i in 0..=near {
                for j in 0..=far {
                    // (2R − w)^j = Σ_l C(j,l) (2R)^{j−l} (−w)^l: alternating, but w is
                    // concentrated within 1/|δ| ≪ 2R of zero here, so each term is smaller than
                    // the last by about 1/(2|δ|R) ≤ ¼.
                    let mut tail = 0.0;
                    for l in 0..=j {
                        let sign = if l & 1 == 0 { 1.0 } else { -1.0 };
                        tail += sign * binomial(j, l) * (2.0 * r).powi((j - l) as i32) * w[i + l];
                    }
                    inner += binomial(near, i) * binomial(far, j) * u[near - i + far - j] * tail;
                }
            }
            shielding += bk * inner / 2f64.powi((a + k) as i32);
        }
    }
    // dr dt = du dv / 2; the exponentials taken out of both forms combine to e^{−2 min(ζ) R}.
    let shielding = prefactor * shielding * 0.5 / (2.0 * r) * (-2.0 * zeta_a.min(zeta_b) * r).exp();
    potential_of_a - shielding
}

/// `C(n, k)` as a float, exact for the small arguments used here.
fn binomial(n: usize, k: usize) -> f64 {
    (0..k).fold(1.0, |c, i| c * (n - i) as f64 / (i + 1) as f64)
}

/// The shielded Coulomb interaction `J_AB` (eqs 14–17), in **eV**, between an atom of `a` with
/// Slater exponent `zeta_a` and one of `b` with `zeta_b` (bohr⁻¹), `r` **metres** apart:
/// [`coulomb_integral`] with each element's [`principal_quantum_number`], in the paper's units.
pub fn shielded_coulomb(a: Element, zeta_a: f64, b: Element, zeta_b: f64, r: f64) -> f64 {
    HARTREE_EV
        * coulomb_integral(
            principal_quantum_number(a),
            zeta_a,
            principal_quantum_number(b),
            zeta_b,
            r / BOHR,
        )
}

/// Each atom's Slater exponent at hydrogen charges `charges`: eq 17′ for every element, plus
/// `Q_H` for hydrogen (eq 20).
fn exponents(elements: &[Element], charges: &[f64]) -> Vec<f64> {
    elements
        .iter()
        .zip(charges)
        .map(|(&e, &q)| {
            let zeta = slater_exponent(e);
            if e == Element::H {
                zeta + q
            } else {
                zeta
            }
        })
        .collect()
}

/// The matrix `J_AB` of eq 6′, in eV, row-major, `n × n` for `n` atoms at positions `at`
/// (metres), with every hydrogen's orbital at the charge `charges` gives it.
///
/// Off the diagonal, [`shielded_coulomb`] at each pair's distance, with hydrogen's ζ by eq 20.
/// On it, Table I's idempotential J⁰ — **not** the Slater self-integral — and for hydrogen eq 21,
/// `(1 + Q_H/ζ⁰_H) J⁰_HH`. Symmetric by construction: each pair is computed once.
///
/// # Panics
///
/// If `at` or `charges` is not one per element, or a position is not finite.
pub fn coulomb_matrix(elements: &[Element], at: &[[f64; 3]], charges: &[f64]) -> Vec<f64> {
    let n = elements.len();
    assert_eq!(at.len(), n, "one position per atom");
    assert_eq!(charges.len(), n, "one charge per atom");
    let zeta = exponents(elements, charges);
    let mut j = vec![0.0; n * n];
    for a in 0..n {
        let row = table_i(elements[a]);
        j[a * n + a] = if elements[a] == Element::H {
            (1.0 + charges[a] / slater_exponent(Element::H)) * row.idempotential
        } else {
            row.idempotential
        };
        for b in a + 1..n {
            let d = [
                at[a][0] - at[b][0],
                at[a][1] - at[b][1],
                at[a][2] - at[b][2],
            ];
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            let jab = shielded_coulomb(elements[a], zeta[a], elements[b], zeta[b], r);
            j[a * n + b] = jab;
            j[b * n + a] = jab;
        }
    }
    j
}

/// Why charges could not be equilibrated.
#[derive(Clone, Debug, PartialEq)]
pub enum QeqError {
    /// The hydrogen iteration did not settle: after `iterations` solves the charges still moved by
    /// `change` elementary charges between the last two.
    NotConverged {
        /// Solves taken.
        iterations: usize,
        /// The largest change of any charge between the last two, e.
        change: f64,
    },
    /// The linear system has no unique solution. Reachable: two neutral hydrogens 0.7915 Å
    /// apart, where their shielded J equals the idempotential J⁰_HH and eq 6′'s quadratic form
    /// for moving charge between them, `J_AA + J_BB − 2 J_AB`, is zero. That is *longer* than
    /// H₂'s 0.741 Å bond, so at its bond length the form is negative at fixed orbitals; QEq's
    /// symmetric answer for H₂, zero, is then a stationary point and not a minimum.
    Singular,
    /// A charge came out NaN or infinite — from a total charge that is not a finite number.
    NonFinite,
    /// The total charge cannot be reached with every charge inside its eq 5 range.
    Infeasible {
        /// The total asked for, e.
        total_charge: f64,
        /// The lowest total the ranges allow, e.
        lowest: f64,
        /// The highest, e.
        highest: f64,
    },
}

impl fmt::Display for QeqError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            QeqError::NotConverged { iterations, change } => write!(
                f,
                "the hydrogen charges had not settled after {iterations} solves (last change \
                 {change:e} e)"
            ),
            QeqError::Singular => f.write_str("the charge-equilibration equations are singular"),
            QeqError::NonFinite => f.write_str("a charge came out NaN or infinite"),
            QeqError::Infeasible {
                total_charge,
                lowest,
                highest,
            } => write!(
                f,
                "a total charge of {total_charge} e is outside what the atoms' ranges allow \
                 ({lowest} to {highest} e)"
            ),
        }
    }
}

impl std::error::Error for QeqError {}

/// Equilibrated charges and what the solve says about them.
#[derive(Clone, Debug, PartialEq)]
pub struct Charges {
    /// One partial charge per atom, in elementary charges, summing to the total asked for.
    pub charges: Vec<f64>,
    /// The common chemical potential of the free atoms (eq 8), eV, from the last solve.
    pub electronegativity: f64,
    /// Linear solves taken: one without hydrogen, more with.
    pub iterations: usize,
    /// The atoms whose charge was fixed at an end of its eq 5 range, in index order.
    pub at_bound: Vec<usize>,
}

/// The QEq solver's convergence settings. [`Qeq::default`] is what [`charges`] and
/// [`ForceField::with_qeq_charges`](crate::energy::ForceField::with_qeq_charges) use.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Qeq {
    /// The hydrogen iteration has converged when no charge changes by more than this between two
    /// solves, in elementary charges. Default 10⁻¹⁰.
    pub tolerance: f64,
    /// Solves before giving up with [`QeqError::NotConverged`]. Default 100.
    pub max_iterations: usize,
}

impl Default for Qeq {
    fn default() -> Qeq {
        Qeq {
            tolerance: 1e-10,
            max_iterations: 100,
        }
    }
}

impl Qeq {
    /// The equilibrated charges of atoms `elements` at positions `at` (metres) with total charge
    /// `total_charge` (e). See the module documentation for the method.
    ///
    /// # Errors
    ///
    /// [`QeqError::Infeasible`] if the total is outside what the eq 5 ranges allow,
    /// [`QeqError::NonFinite`] if a charge comes out NaN or infinite (a NaN total),
    /// [`QeqError::Singular`] if the equations have no unique solution, and
    /// [`QeqError::NotConverged`] if the hydrogen iteration has not settled in
    /// [`Qeq::max_iterations`] solves.
    ///
    /// # Panics
    ///
    /// If `at` is not one per element, or a position is not finite.
    pub fn equilibrate(
        &self,
        elements: &[Element],
        at: &[[f64; 3]],
        total_charge: f64,
    ) -> Result<Charges, QeqError> {
        let n = elements.len();
        assert_eq!(at.len(), n, "one position per atom");
        assert!(
            at.iter().flatten().all(|x| x.is_finite()),
            "positions must be finite"
        );
        let ranges: Vec<(f64, f64)> = elements.iter().map(|&e| charge_range(e)).collect();
        let lowest: f64 = ranges.iter().map(|r| r.0).sum();
        let highest: f64 = ranges.iter().map(|r| r.1).sum();
        if n == 0 || total_charge < lowest || total_charge > highest {
            return Err(QeqError::Infeasible {
                total_charge,
                lowest,
                highest,
            });
        }
        let chi: Vec<f64> = elements
            .iter()
            .map(|&e| table_i(e).electronegativity)
            .collect();
        let has_hydrogen = elements.contains(&Element::H);
        let mut charges = vec![0.0; n];
        for iteration in 1..=self.max_iterations {
            let j = coulomb_matrix(elements, at, &charges);
            let (next, mu, at_bound) = solve(&j, &chi, &ranges, total_charge)?;
            // A NaN would pass every comparison below — `max` discards it and `<=` is false —
            // so it is refused by name before it can be called converged or reach an exponent.
            if !next.iter().all(|q| q.is_finite()) {
                return Err(QeqError::NonFinite);
            }
            let change = next
                .iter()
                .zip(&charges)
                .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
            charges = next;
            if !has_hydrogen || change <= self.tolerance {
                return Ok(Charges {
                    charges,
                    electronegativity: mu,
                    iterations: iteration,
                    at_bound,
                });
            }
            if iteration == self.max_iterations {
                return Err(QeqError::NotConverged {
                    iterations: iteration,
                    change,
                });
            }
        }
        Err(QeqError::NotConverged {
            iterations: 0,
            change: f64::NAN,
        })
    }
}

/// The QEq charges of `component` at positions `at` (metres), with [`Qeq::default`], and the
/// total charge the dictionary states atom by atom: the sum of its formal charges.
///
/// # Errors
///
/// As [`Qeq::equilibrate`].
///
/// # Panics
///
/// If `at` is not one position per atom.
pub fn charges(component: &Component, at: &[[f64; 3]]) -> Result<Charges, QeqError> {
    let elements: Vec<Element> = component.atoms().iter().map(|a| a.element).collect();
    let total: i32 = component.atoms().iter().map(|a| a.charge).sum();
    Qeq::default().equilibrate(&elements, at, f64::from(total))
}

/// One solve of eqs 10–13 at a fixed `J` (row-major, eV): the charges, the common
/// electronegativity, and the atoms fixed at a range boundary.
fn solve(
    j: &[f64],
    chi: &[f64],
    ranges: &[(f64, f64)],
    total: f64,
) -> Result<(Vec<f64>, f64, Vec<usize>), QeqError> {
    let n = chi.len();
    let mut fixed: Vec<Option<f64>> = vec![None; n];
    loop {
        let free: Vec<usize> = (0..n).filter(|&i| fixed[i].is_none()).collect();
        let fixed_total: f64 = fixed.iter().flatten().sum();
        if free.is_empty() {
            // Unreachable while the feasibility check holds, but stated rather than assumed.
            let (lowest, highest) = ranges
                .iter()
                .fold((0.0, 0.0), |(l, h), r| (l + r.0, h + r.1));
            return Err(QeqError::Infeasible {
                total_charge: total,
                lowest,
                highest,
            });
        }
        // Eq 13: each fixed atom's field folded into the free atoms' electronegativity.
        let chi_f: Vec<f64> = (0..n)
            .map(|a| {
                chi[a]
                    + (0..n)
                        .filter_map(|b| fixed[b].map(|q| j[a * n + b] * q))
                        .sum::<f64>()
            })
            .collect();
        let m = free.len();
        let first = free[0];
        // Eqs 10–12, with C_1j = 1 (eq 12 prints Q_j; see the module documentation).
        let mut c = vec![0.0; m * m];
        let mut rhs = vec![0.0; m];
        c[..m].fill(1.0);
        rhs[0] = total - fixed_total;
        for (k, &i) in free.iter().enumerate().skip(1) {
            for (l, &jj) in free.iter().enumerate() {
                c[k * m + l] = j[i * n + jj] - j[first * n + jj];
            }
            rhs[k] = -(chi_f[i] - chi_f[first]);
        }
        let q = gauss(c, rhs)?;
        let mut moved = false;
        for (k, &i) in free.iter().enumerate() {
            let (lo, hi) = ranges[i];
            if q[k] < lo {
                fixed[i] = Some(lo);
                moved = true;
            } else if q[k] > hi {
                fixed[i] = Some(hi);
                moved = true;
            }
        }
        if moved {
            continue;
        }
        let mut charges = vec![0.0; n];
        for (i, f) in fixed.iter().enumerate() {
            if let Some(v) = f {
                charges[i] = *v;
            }
        }
        for (k, &i) in free.iter().enumerate() {
            charges[i] = q[k];
        }
        let mu = chi[first] + (0..n).map(|b| j[first * n + b] * charges[b]).sum::<f64>();
        let at_bound = (0..n).filter(|&i| fixed[i].is_some()).collect();
        return Ok((charges, mu, at_bound));
    }
}

/// Solves `a x = b` (`a` row-major, square) by Gaussian elimination with partial pivoting.
/// [`QeqError::Singular`] when a pivot is at most 10⁻¹² of the largest entry.
fn gauss(mut a: Vec<f64>, mut b: Vec<f64>) -> Result<Vec<f64>, QeqError> {
    let m = b.len();
    let scale = a.iter().fold(0.0f64, |s, x| s.max(x.abs()));
    for col in 0..m {
        let mut pivot = col;
        for row in col + 1..m {
            if a[row * m + col].abs() > a[pivot * m + col].abs() {
                pivot = row;
            }
        }
        let size = a[pivot * m + col].abs();
        if size.is_nan() || size <= 1e-12 * scale {
            return Err(QeqError::Singular);
        }
        if pivot != col {
            for k in 0..m {
                a.swap(col * m + k, pivot * m + k);
            }
            b.swap(col, pivot);
        }
        for row in col + 1..m {
            let f = a[row * m + col] / a[col * m + col];
            if f != 0.0 {
                for k in col..m {
                    a[row * m + k] -= f * a[col * m + k];
                }
                b[row] -= f * b[col];
            }
        }
    }
    let mut x = vec![0.0; m];
    for row in (0..m).rev() {
        let mut s = b[row];
        for k in row + 1..m {
            s -= a[row * m + k] * x[k];
        }
        x[row] = s / a[row * m + row];
    }
    Ok(x)
}
