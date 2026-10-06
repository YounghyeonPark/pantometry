//! Ewald-summed electrostatics in an orthorhombic periodic box: real space, reciprocal space, the
//! self term, the excluded pairs and the neutralising background, with analytic forces and the
//! virial — and the complementary error function they need, written here.
//!
//! # The sum
//!
//! For charges `q_i` (elementary charges) at `r_i` in a box of volume `V`, repeated without end,
//! with the conducting ("tinfoil") boundary at infinity, the electrostatic energy is split by a
//! Gaussian of width `1/α` (Ewald 1921; the form below is de Leeuw, Perram and Smith, *Proc. R.
//! Soc. A* **373**, 27 (1980), as every textbook writes it — Allen and Tildesley, *Computer
//! Simulation of Liquids*, 2nd ed., §6.2, was the one at hand), in units of the Coulomb constant
//! `k_e` = 332.0637 kcal mol⁻¹ Å e⁻² ([`COULOMB`]):
//!
//! ```text
//! E = E_real + E_recip + E_self + E_excl + E_bg
//! E_real  =  Σ_{i<j, not excluded, |r_ij| < r_c}  q_i q_j erfc(α r_ij) / r_ij        (minimum image)
//! E_recip =  (2π/V) Σ_{k ≠ 0, |k| ≤ k_c}  exp(−k²/4α²)/k²  |S(k)|²,   S(k) = Σ_j q_j exp(i k·r_j)
//! E_self  = −(α/√π) Σ_i q_i²
//! E_excl  = −Σ_{excluded i<j}  q_i q_j erf(α r_ij) / r_ij                             (minimum image)
//! E_bg    = −π Q_net² / (2 V α²)
//! ```
//!
//! `k = 2π (h/L_x, k/L_y, l/L_z)` for integers. Each piece and why it is there:
//!
//! - **The excluded pairs.** UFF leaves out 1-2 and 1-3 pairs entirely (see [`crate::energy`]),
//!   and 1-4 pairs at full strength, so there is no scaled pair: a pair is in, or out. The
//!   reciprocal sum is over all pairs, excluded ones included, because `|S(k)|²` cannot tell them
//!   apart; `E_excl` takes their smooth part, `erf(αr)/r`, back out — exactly, at any `α`, so an
//!   excluded pair contributes nothing at all. Every excluded pair is corrected wherever it is,
//!   not only inside the cutoff.
//! - **The self term** removes each charge's interaction with its own Gaussian, which the
//!   reciprocal sum contains: `lim_{r→0} q² erf(αr)/r = 2α q²/√π`, halved because the pair sum is
//!   over `i < j`.
//! - **The background** is the energy of a uniform charge `−Q_net/V` neutralising a box whose
//!   charges do not sum to zero; the `k = 0` term, which is otherwise infinite, is left out by
//!   it. It is the term that makes a net charge's energy independent of `α`: a single charge in
//!   a cube of side `L` then has `E = k_e q² ξ / (2L)` with `ξ = −2.837297479480620`, the Wigner
//!   constant of the simple cubic lattice, which `tests/the_ewald_sum_against_closed_forms.rs`
//!   checks across `α`. It is a property of the model — every Ewald code with a net charge has
//!   it — and a finite-size artefact of a charged box, not a physical energy.
//! - **The tinfoil boundary.** Leaving out the `k = 0` term is the conducting boundary; a box
//!   carrying a dipole `M` in vacuum would add `2π M²/(3V)`. So the periodic energy of one neutral
//!   molecule of dipole `μ` in a large box is its non-periodic energy **less** `2π μ²/(3V)`, and
//!   the image terms after that fall as `L⁻⁵` (the cubic lattice sum of dipole–dipole terms is
//!   zero, and dipole–quadrupole is odd). That is the rate the tests check.
//!
//! # Forces and the virial
//!
//! Analytic, each piece differentiated as written: a pair's `dE/dr` for real space and the
//! excluded pairs, and for reciprocal space
//! `F_i = (4π/V) q_i Σ_k exp(−k²/4α²)/k² · k · Im(S(k)* exp(i k·r_i))`. The **virial** is
//! `W_ab = −∂E/∂ε_ab`, the derivative under a homogeneous strain of the box and every atom with
//! it, so that `P_ab V = Σ m v_a v_b + W_ab`: `Σ r_a f_b` over the pairs, and for reciprocal
//! space, from the same sum, `(2π/V) Σ_k exp(−k²/4α²)/k² |S|² [δ_ab − 2 k_a k_b (1/k² + 1/4α²)]`
//! (this is the strain derivative of `E_recip` at fixed `S`, since `k·r` does not change under the
//! strain, derived here and the form Nosé and Klein, *Mol. Phys.* **50**, 1055 (1983), give). The
//! background, `∝ 1/V`, contributes `E_bg δ_ab`; the self term does not depend on the box.
//!
//! # Choosing α, r_c and k_c: Kolafa and Perram's estimates
//!
//! Kolafa and Perram, *Mol. Simul.* **9**, 351 (1992), estimate the RMS error a truncation
//! leaves in the energy of a disordered system, with `Q = Σ q_i²` — **read secondarily**, as
//! Saffar Shamshirgar, Hess and Tornberg quote them (arXiv:1712.04718, eqs 2.1b and 2.2b); the
//! paper was not opened:
//!
//! ```text
//! δE_real  ≈ Q (r_c / 2V)^½ (α r_c)⁻² exp(−α² r_c²)
//! δE_recip ≈ Q α π⁻² K^(−3/2) exp(−(π K / α L)²)        K = k_c L / 2π
//! ```
//!
//! Both are proportional to `Q/r_c`, so **the accuracy parameter `δ` here is each estimate in
//! units of `k_e Q / r_c`**, dimensionless and independent of the charges:
//! [`EwaldParameters::for_accuracy`] takes `r_c` and `δ`, solves the first for `α` and the second
//! for `k_c` (bisection; `L = V^⅓` for a box that is not a cube, and a spherical cutoff in `k`),
//! and the expected RMS error of the energy is then about `√2 δ k_e Q / r_c`
//! ([`EwaldParameters::energy_error`]). On a box of a thousand TIP3P-like waters, 31 Å across,
//! `r_c` = 9 Å and `δ = 10⁻⁵` give `α` = 0.301 Å⁻¹ and 1 051 wave vectors.
//!
//! **What the estimates do not include, derived here**: the reciprocal truncation also leaves out
//! the diagonal of `|S(k)|²`, which for a disordered system tends to `Q` at large `k`, so
//! the omitted sum has a mean as well as a spread. Integrating `(2π/V) Q exp(−k²/4α²)/k²` over
//! `|k| > k_c` gives `k_e Q (α/√π) erfc(k_c/2α)` ([`EwaldParameters::reciprocal_bias`]): an
//! energy lowered by that, which is a constant at fixed charges and box — it moves no force and no
//! difference at fixed charges — and is larger than the spread by about `α L K^½`. **Measured, it
//! is the error**: on a disordered box of 96 charges the energy at each accuracy from `10⁻⁶` to
//! `10⁻¹¹` was below the tightest by 0.88–1.10 of the bias, ten times the estimated spread; on a lone
//! charge, where it is the whole truncation, it is 0.78–1.05 of the omitted shell. It has no cross
//! term between two sets of charges, so a decoupling's `∂U/∂λ_e` does not see it. The tests'
//! bounds include it. In a crystal it is not there as such: the off-diagonal terms reshape it, and
//! the estimates do not hold (see the Madelung test). Kolafa and Perram's real-space estimate
//! does hold: over sixteen disordered boxes the RMS error was 1.02–1.06 of it.
//!
//! Both sums are cut sharply, as the estimates assume: a pair at `r_c` and a wave vector at `k_c`
//! enter or leave whole. **The energy is therefore not continuous**: a pair of unit charges
//! crossing `r_c` moves it by `k_e erfc(α r_c)/r_c`, about `δ k_e (2V/r_c³)^½ / r_c` at the
//! chosen `α`.
//!
//! # erfc, written here
//!
//! [`erfc`] is original to this crate, written from published mathematics and fitted with
//! mpmath, with no library's code or table read: `tools/erfc-reference/generate.py` writes every
//! constant below and says how. For `|x| < 0.4375`, `erfc x = 1 − erf x` with erf by its Taylor
//! series (Abramowitz and Stegun 7.1.5), fourteen terms, where `erfc ≥ 0.535` so the subtraction
//! loses nothing. From 0.4375 to 6, `erfc x = exp(−x²) erfcx(x)`, the scaled complement by a
//! Chebyshev series on each of nine intervals (Clenshaw's recurrence), and beyond 6 the same with
//! `x erfcx(x)` a Chebyshev series in `u = 1/x²` on `[0, 1/36]`, where its value at `u = 0` is
//! `1/√π`, the asymptotic series' (A&S 7.1.23) first term. `exp(−x²)` is
//! [`pantometry_core::math::exp`] of `−x²` with `x²` carried exactly as a sum of two doubles
//! (Dekker's product, split by Veltkamp's constant `2²⁷ + 1`, no fused multiply–add), so that
//! the square's rounding, which would be `x²` ulps at large `x`, costs nothing. Negative
//! arguments use `erfc(−x) = 2 − erfc(x)`.
//!
//! **Its accuracy, against mpmath at 50 digits**: at most **3.21 ulp** over the 992 240 points
//! with a normal result of the million `generate.py --measure` was run on once — 2.74 from 0.4375
//! to 6, 3.21 above 6, 1.22 for `|x| < 0.4375` and 1.05 below −0.4375 — and the 438 reference points
//! `tests/the_ewald_sum_against_closed_forms.rs` carries, whose worst is 2.67, are held to
//! [`ERFC_ULP_BOUND`]. Above about 26.55 the result is subnormal, and it is not measured there;
//! above 27.3 it is zero. Not correctly rounded and not
//! claimed to be: the Ewald sum asks it for `10⁻⁶`, and gets a margin of nine orders. It is
//! **the same function on every platform**, built from `+ − × ÷` and the kernel's `exp`.
//!
//! # Determinism, and where it stops
//!
//! No threads, no hashing: every sum runs over its terms in an order fixed by the atom indices
//! and the integer wave vectors, so an evaluation is the same bits however often and from wherever
//! it is called. The reciprocal sum's phase factors start from **the platform's `sin` and `cos`**
//! of each atom's three fractional coordinates — 3N calls an evaluation — and are carried to
//! higher `k` by complex multiplication, so across platforms the reciprocal energy can differ in
//! its last bits, as UFF's torsion `cos` already makes every trajectory of this crate do.
//!
//! # Cost, and what is not here
//!
//! Classical Ewald: the reciprocal sum is `N` times the number of wave vectors, which at a fixed
//! accuracy grows as the box's volume, so at fixed `r_c` the cost goes as `N²`. Smooth particle-mesh
//! Ewald (Essmann et al., *J. Chem. Phys.* **103**, 8577 (1995)) is the production method and is
//! **not here**; this sum is the reference it will be checked against. See [`crate::periodic`]
//! for the measured cost.

use crate::ccd::ANGSTROM;
use crate::energy::COULOMB_KCAL;
use crate::periodic::{for_each_pair, Exclusions, PeriodicBox};
use crate::uff::KCAL_PER_MOL;
use pantometry_core::math::exp;
use std::f64::consts::{FRAC_2_SQRT_PI, PI};

/// The Coulomb constant `k_e = e²/(4π ε₀)` in joules × metres per elementary charge squared: eq 43's
/// 332.0637 kcal mol⁻¹ Å e⁻² in SI, the prefactor [`crate::energy::coulomb`] uses.
pub const COULOMB: f64 = COULOMB_KCAL * KCAL_PER_MOL * ANGSTROM;

/// The bound [`erfc`] is held to against mpmath's values in the tests, in ulps of the exact value:
/// the measured worst, 3.21, with a margin.
pub const ERFC_ULP_BOUND: f64 = 3.5;

/// Below this, `erfc = 1 − erf` by the series; at and above it, `exp(−x²) erfcx(x)`.
const SMALL: f64 = 0.4375;

/// At and above this, `erfcx` is a series in `1/x²`.
const FAR: f64 = 6.0;

/// Above this, `erfc x` is below half the smallest subnormal and is zero.
const UNDERFLOW: f64 = 27.3;

/// The intervals the Chebyshev series of `erfcx` cover, `[EDGES[i], EDGES[i + 1])`.
const EDGES: [f64; 10] = [0.4375, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0];

/// erf's Taylor coefficients, `(2/√π) (−1)ⁿ / (n! (2n + 1))`, `generate.py --tables`.
/// The first is `2/√π`, as the script computes it.
#[allow(clippy::approx_constant)]
const ERF_SERIES: [f64; 14] = [
    1.1283791670955126,
    -0.37612638903183754,
    0.11283791670955126,
    -0.026866170645131252,
    0.005223977625442188,
    -0.0008548327023450853,
    0.00012055332981789664,
    -1.492565035840625e-05,
    1.6462114365889248e-06,
    -1.6365844691234924e-07,
    1.4807192815879218e-08,
    -1.2290555301717928e-09,
    9.422759064650411e-11,
    -6.7113668551641105e-12,
];
/// How many of each row of [`ERFCX_TABLE`] are used.
const ERFCX_TERMS: [usize; 9] = [16, 15, 14, 14, 13, 13, 13, 15, 14];
/// Chebyshev coefficients of `erfcx` on each interval of [`EDGES`], `generate.py --tables`.
const ERFCX_TABLE: [[f64; 16]; 9] = [
    [
        0.5284696052263012,
        -0.1100111633368543,
        0.009860772787737379,
        -0.0007897379033177107,
        5.778492104114415e-05,
        -3.920685048681516e-06,
        2.4931202835501284e-07,
        -1.497702396459946e-08,
        8.552799194438779e-10,
        -4.6660918264601724e-11,
        2.4419185791636457e-12,
        -1.2300396821647068e-13,
        5.98092566855994e-15,
        -2.814185817661824e-16,
        1.2841147636410809e-17,
        -5.692954190414895e-19,
    ],
    [
        0.3711927943420716,
        -0.052799551940465685,
        0.00338076236159189,
        -0.00019896338415156294,
        1.0911979120987141e-05,
        -5.632372066195936e-07,
        2.7564656102315984e-08,
        -1.2864435796900215e-09,
        5.75189434369204e-11,
        -2.473180238194912e-12,
        1.0258794609509988e-13,
        -4.116181633727305e-15,
        1.6012087659273513e-16,
        -6.050937829047491e-18,
        2.225245070154316e-19,
        0.0,
    ],
    [
        0.28672723893640484,
        -0.03300680108481147,
        0.0017591477886840761,
        -8.788195718624664e-05,
        4.1516102756112e-06,
        -1.86697213290856e-07,
        8.033784461538863e-09,
        -3.321836144278999e-10,
        1.3243525744778221e-11,
        -5.105578892334744e-13,
        1.9079334307984007e-14,
        -6.925820759722678e-16,
        2.4466018310021198e-17,
        -8.424358002974198e-19,
        0.0,
        0.0,
    ],
    [
        0.23209235631383532,
        -0.02225125762219511,
        0.0010068780404189428,
        -4.3327905521453394e-05,
        1.783140732436e-06,
        -7.049777846240933e-08,
        2.687293297212188e-09,
        -9.906296971849274e-11,
        3.5405028248059043e-12,
        -1.2294649407955319e-13,
        4.1560224390309084e-15,
        -1.3698097771178357e-16,
        4.408482868322765e-18,
        -1.3871407881015812e-19,
        0.0,
        0.0,
    ],
    [
        0.1942820812745081,
        -0.01587924788333761,
        0.0006208296616878289,
        -2.3328859409806277e-05,
        8.456745835133116e-07,
        -2.9663701263111685e-08,
        1.0094101754255728e-09,
        -3.3394274912569374e-11,
        1.0760912211208994e-12,
        -3.383019820369886e-14,
        1.039102577686144e-15,
        -3.122195341440062e-17,
        9.187550752854276e-19,
        0.0,
        0.0,
        0.0,
    ],
    [
        0.16674094078984744,
        -0.011840246176534596,
        0.00040602706605389184,
        -1.3487936862426611e-05,
        4.351207693856568e-07,
        -1.3660344480880951e-08,
        4.1810593877867705e-10,
        -1.249586682885637e-11,
        3.651747963420464e-13,
        -1.0447676413633131e-14,
        2.929516020050166e-16,
        -8.058488951045153e-18,
        2.1765870687733285e-19,
        0.0,
        0.0,
        0.0,
    ],
    [
        0.1458678994299681,
        -0.00913882677497607,
        0.00027841748283800253,
        -8.265426896062553e-06,
        2.395159038287142e-07,
        -6.7849288502505444e-09,
        1.8813155805500703e-10,
        -5.111902689715732e-12,
        1.3625532798777336e-13,
        -3.565943642783521e-15,
        9.170794861729046e-17,
        -2.3194193951852586e-18,
        5.77285903440333e-20,
        0.0,
        0.0,
        0.0,
    ],
    [
        0.12316667236819653,
        -0.013112375841864296,
        0.0006836175123551547,
        -3.494745631765322e-05,
        1.7535988710219397e-06,
        -8.644752622984298e-08,
        4.190236831762333e-09,
        -1.9985303631241498e-10,
        9.385613797592563e-12,
        -4.342729824197335e-13,
        1.9808776494082428e-14,
        -8.911984215234841e-16,
        3.9566098188505085e-17,
        -1.7341995762262654e-18,
        7.507293529841913e-20,
        0.0,
    ],
    [
        0.10135070425173438,
        -0.008947312860317571,
        0.0003891909669137226,
        -1.6692555616261904e-05,
        7.063294204001365e-07,
        -2.9500674190544564e-08,
        1.2167353483368035e-09,
        -4.9577500529575074e-11,
        1.9964996401752625e-12,
        -7.948956254073044e-14,
        3.130090831054916e-15,
        -1.2194106991518773e-16,
        4.701332933672435e-18,
        -1.7942952656092543e-19,
        0.0,
        0.0,
    ],
];
/// Chebyshev coefficients of `x erfcx(x)` in `u = 1/x²` on `[0, 1/36]`, `generate.py --tables`.
const ERFCX_FAR: [f64; 13] = [
    0.5603874949378894,
    -0.0037644967870190307,
    3.698729804476604e-05,
    -5.912293529244797e-07,
    1.2927848568567902e-08,
    -3.554450382700824e-10,
    1.169123615029381e-11,
    -4.4516637926929027e-13,
    1.9171763060217432e-14,
    -9.178540101208144e-16,
    4.820238777304872e-17,
    -2.7476072489927693e-18,
    1.685409989600857e-19,
];

/// Clenshaw's recurrence for `Σ c_k T_k(t)`.
fn clenshaw(c: &[f64], t: f64) -> f64 {
    let mut b1 = 0.0;
    let mut b2 = 0.0;
    for &ck in c[1..].iter().rev() {
        let b0 = 2.0 * t * b1 - b2 + ck;
        b2 = b1;
        b1 = b0;
    }
    t * b1 - b2 + c[0]
}

/// `x²` as `hi + lo` exactly: Dekker's product, with Veltkamp's split of `x` into two halves of
/// 26 bits whose products are exact. Valid for `|x|` below about `2⁹⁹⁵`.
fn exact_square(x: f64) -> (f64, f64) {
    let p = x * x;
    let c = 134_217_729.0 * x;
    let hi = c - (c - x);
    let lo = x - hi;
    (p, ((hi * hi - p) + 2.0 * hi * lo) + lo * lo)
}

/// erf by its Taylor series, for `|x| < SMALL`.
fn erf_series(x: f64) -> f64 {
    let z = x * x;
    let mut s = 0.0;
    for &c in ERF_SERIES.iter().rev() {
        s = s * z + c;
    }
    x * s
}

/// `erfcx(x) = exp(x²) erfc(x)` for `SMALL ≤ x ≤ UNDERFLOW`.
fn erfcx(x: f64) -> f64 {
    if x >= FAR {
        let u = 1.0 / (x * x);
        // u in (0, 1/36]: t = 72 u − 1, in (−1, 1].
        return clenshaw(&ERFCX_FAR, (u - 1.0 / 72.0) * 72.0) / x;
    }
    let i = if x < 1.0 {
        0
    } else if x < 4.0 {
        (2.0 * x) as usize - 1
    } else if x < 5.0 {
        7
    } else {
        8
    };
    let (a, b) = (EDGES[i], EDGES[i + 1]);
    let t = (x - 0.5 * (a + b)) * (2.0 / (b - a));
    clenshaw(&ERFCX_TABLE[i][..ERFCX_TERMS[i]], t)
}

/// The three values one evaluation of the Ewald kernels needs at `x ≥ 0`: `erfc x`, `erf x` and
/// `exp(−x²)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct ErrorFunctions {
    pub erfc: f64,
    pub erf: f64,
    pub gaussian: f64,
}

/// [`ErrorFunctions`] at `x ≥ 0`, sharing the one `exp`.
pub(crate) fn error_functions(x: f64) -> ErrorFunctions {
    debug_assert!(x >= 0.0);
    let (hi, lo) = exact_square(x);
    let e = exp(-hi);
    // exp(−hi − lo) = e (1 − lo) to within lo², which is below 2⁻¹⁰⁰ of it here.
    let gaussian = e - e * lo;
    if x < SMALL {
        let erf = erf_series(x);
        return ErrorFunctions {
            erfc: 1.0 - erf,
            erf,
            gaussian,
        };
    }
    if x > UNDERFLOW {
        return ErrorFunctions {
            erfc: 0.0,
            erf: 1.0,
            gaussian,
        };
    }
    let r = erfcx(x);
    let erfc = e * (r - r * lo);
    ErrorFunctions {
        erfc,
        erf: 1.0 - erfc,
        gaussian,
    }
}

/// The complementary error function, `erfc x = (2/√π) ∫ₓ^∞ exp(−t²) dt`: see the module
/// documentation for the method and its measured accuracy (at most 3.21 ulp against mpmath).
/// `erfc(NaN)` is NaN, `erfc(+∞) = 0`, `erfc(−∞) = 2`.
pub fn erfc(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x >= 0.0 {
        if x < SMALL {
            1.0 - erf_series(x)
        } else if x > UNDERFLOW {
            0.0
        } else {
            let (hi, lo) = exact_square(x);
            let e = exp(-hi);
            let r = erfcx(x);
            e * (r - r * lo)
        }
    } else if x > -SMALL {
        1.0 + erf_series(-x)
    } else {
        2.0 - erfc(-x)
    }
}

/// The error function, `erf x = 1 − erfc x`: by its series for `|x| < 0.4375`, where that is more
/// accurate than the difference, and as `1 − erfc x` beyond.
pub fn erf(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x.abs() < SMALL {
        erf_series(x)
    } else if x > 0.0 {
        1.0 - erfc(x)
    } else {
        erfc(-x) - 1.0
    }
}

/// Bisection for the `y` in `[lo, hi]` where a decreasing `f` crosses zero: a fixed 200 halvings,
/// or until the bracket stops shrinking.
fn bisect(f: impl Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> f64 {
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            break;
        }
        if f(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    hi
}

/// The three numbers that fix an Ewald sum: the splitting parameter, the real-space cutoff and the
/// reciprocal-space cutoff. See the module documentation for how
/// [`EwaldParameters::for_accuracy`] chooses them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EwaldParameters {
    /// α, per metre: the inverse width of the screening Gaussian.
    pub alpha: f64,
    /// `r_c`, metres: real-space pairs at or beyond it are left out. At most half the box's
    /// shortest edge, so that a pair has one image inside it.
    pub cutoff: f64,
    /// `k_c`, per metre: wave vectors longer than it are left out.
    pub k_cutoff: f64,
}

impl EwaldParameters {
    /// α and `k_c` for real-space cutoff `cutoff` (metres) in `cell`, chosen so that each of Kolafa
    /// and Perram's two RMS energy-error estimates is `accuracy` in units of `k_e Q / r_c`. See the
    /// module documentation.
    ///
    /// # Panics
    ///
    /// If `accuracy` is not in `(0, 0.01]`, or `cutoff` is not positive and at most half the box's
    /// shortest edge.
    pub fn for_accuracy(cell: &PeriodicBox, cutoff: f64, accuracy: f64) -> EwaldParameters {
        assert!(
            accuracy > 0.0 && accuracy <= 0.01,
            "an Ewald accuracy must be in (0, 0.01]"
        );
        assert!(
            cutoff > 0.0 && cutoff <= cell.half_shortest_edge(),
            "the real-space cutoff must be positive and at most half the box's shortest edge"
        );
        let pre = (cutoff * cutoff * cutoff / (2.0 * cell.volume())).sqrt();
        let p = bisect(|p| pre / (p * p) * exp(-p * p) - accuracy, 0.5, 40.0);
        let alpha = p / cutoff;
        let side = cell.volume().cbrt();
        let reciprocal = |y: f64| {
            let k = y * alpha * side / PI;
            (p / (PI * PI)) / (k * k.sqrt()) * exp(-y * y) - accuracy
        };
        let y = bisect(reciprocal, 0.5, 40.0);
        EwaldParameters {
            alpha,
            cutoff,
            k_cutoff: 2.0 * alpha * y,
        }
    }

    /// Kolafa and Perram's RMS real-space energy error in units of `k_e Q / r_c`, in `cell`.
    pub fn real_space_error(&self, cell: &PeriodicBox) -> f64 {
        let p = self.alpha * self.cutoff;
        (self.cutoff.powi(3) / (2.0 * cell.volume())).sqrt() / (p * p) * exp(-p * p)
    }

    /// Kolafa and Perram's RMS reciprocal-space energy error in units of `k_e Q / r_c`, in `cell`,
    /// with `L = V^⅓`.
    pub fn reciprocal_error(&self, cell: &PeriodicBox) -> f64 {
        let y = self.k_cutoff / (2.0 * self.alpha);
        let k = self.k_cutoff * cell.volume().cbrt() / (2.0 * PI);
        (self.alpha * self.cutoff / (PI * PI)) / (k * k.sqrt()) * exp(-y * y)
    }

    /// The expected RMS error of the energy of charges whose squares sum to `sum_q2` (e²), joules:
    /// `k_e Q / r_c` times the two estimates added in quadrature.
    pub fn energy_error(&self, cell: &PeriodicBox, sum_q2: f64) -> f64 {
        let (r, k) = (self.real_space_error(cell), self.reciprocal_error(cell));
        COULOMB * sum_q2 / self.cutoff * (r * r + k * k).sqrt()
    }

    /// The mean of the omitted reciprocal terms for a disordered system, joules:
    /// `k_e Q (α/√π) erfc(k_c / 2α)`, by which the computed energy is low. See the module
    /// documentation; a constant at fixed charges.
    pub fn reciprocal_bias(&self, sum_q2: f64) -> f64 {
        COULOMB
            * sum_q2
            * self.alpha
            * 0.5
            * FRAC_2_SQRT_PI
            * erfc(self.k_cutoff / (2.0 * self.alpha))
    }
}

/// One wave vector of the half space, with what each evaluation needs of it.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Wave {
    /// The integers `(h, k, l)`.
    n: [i32; 3],
    /// `2π (h/L_x, k/L_y, l/L_z)`, per metre.
    k: [f64; 3],
    /// `(4π/V) exp(−k²/4α²)/k²`: the half space's factor 2 included.
    weight: f64,
    /// `2 (1/k² + 1/4α²)`, for the virial.
    stress: f64,
}

/// The electrostatic energy of a periodic system, by part, joules. `total` is accumulated
/// separately from the parts.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EwaldEnergy {
    /// `E_real`: the screened pairs inside the cutoff.
    pub real: f64,
    /// `E_recip`: the wave vectors.
    pub reciprocal: f64,
    /// `E_self`, negative.
    pub self_energy: f64,
    /// `E_excl`: the excluded pairs' smooth part taken back out, negative for like charges.
    pub excluded: f64,
    /// `E_bg`: the neutralising background, zero for a neutral system.
    pub background: f64,
    /// Everything above.
    pub total: f64,
}

/// An Ewald sum's energy, the force on every atom (newtons) and the virial `W_ab = −∂E/∂ε_ab`
/// (joules).
#[derive(Clone, Debug, PartialEq)]
pub struct EwaldEvaluation {
    /// The energy, by part.
    pub energy: EwaldEnergy,
    /// The force on each atom.
    pub forces: Vec<[f64; 3]>,
    /// The virial, symmetric.
    pub virial: [[f64; 3]; 3],
}

/// A classical Ewald sum for one box and one set of [`EwaldParameters`]: the wave vectors,
/// computed once. The charges and positions are given to each evaluation. See the module
/// documentation.
#[derive(Clone, Debug, PartialEq)]
pub struct Ewald {
    cell: PeriodicBox,
    parameters: EwaldParameters,
    waves: Vec<Wave>,
    extent: [usize; 3],
}

impl Ewald {
    /// The sum in `cell` with `parameters`: every wave vector of the half space
    /// (`h > 0`, or `h = 0` and `k > 0`, or `h = k = 0` and `l > 0`) with `|k| ≤ k_c`, in order of
    /// `h`, then `k`, then `l`.
    ///
    /// # Panics
    ///
    /// If the cutoff is past half the box's shortest edge, or α or `k_c` is not positive.
    pub fn new(cell: PeriodicBox, parameters: EwaldParameters) -> Ewald {
        assert!(
            parameters.cutoff > 0.0 && parameters.cutoff <= cell.half_shortest_edge(),
            "the real-space cutoff must be positive and at most half the box's shortest edge"
        );
        assert!(
            parameters.alpha > 0.0 && parameters.k_cutoff > 0.0,
            "α and k_c must be positive"
        );
        let l = cell.lengths();
        let extent = [0, 1, 2].map(|a| (parameters.k_cutoff * l[a] / (2.0 * PI)).floor() as usize);
        let mut integers = Vec::new();
        let [ex, ey, ez] = extent.map(|e| e as i32);
        for h in 0..=ex {
            for k in -ey..=ey {
                for m in -ez..=ez {
                    let upper = h > 0 || (h == 0 && (k > 0 || (k == 0 && m > 0)));
                    if !upper {
                        continue;
                    }
                    let kv = [0, 1, 2].map(|a| 2.0 * PI * f64::from([h, k, m][a]) / l[a]);
                    let k2 = kv[0] * kv[0] + kv[1] * kv[1] + kv[2] * kv[2];
                    if k2 <= parameters.k_cutoff * parameters.k_cutoff {
                        integers.push([h, k, m]);
                    }
                }
            }
        }
        let mut ewald = Ewald {
            cell,
            parameters,
            waves: Vec::new(),
            extent,
        };
        ewald.waves = ewald.waves_for(&integers);
        ewald
    }

    /// The wave vectors `integers` in this sum's box.
    fn waves_for(&self, integers: &[[i32; 3]]) -> Vec<Wave> {
        let l = self.cell.lengths();
        let a2 = 4.0 * self.parameters.alpha * self.parameters.alpha;
        let v = self.cell.volume();
        integers
            .iter()
            .map(|&n| {
                let k = [0, 1, 2].map(|a| 2.0 * PI * f64::from(n[a]) / l[a]);
                let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                Wave {
                    n,
                    k,
                    weight: 4.0 * PI / v * exp(-k2 / a2) / k2,
                    stress: 2.0 * (1.0 / k2 + 1.0 / a2),
                }
            })
            .collect()
    }

    /// The same sum, the same α, `r_c` and **the same integer wave vectors**, in another box: what
    /// a derivative with respect to the box needs, so that no wave vector enters or leaves as the
    /// box changes.
    ///
    /// # Panics
    ///
    /// If the cutoff is past half the new box's shortest edge.
    pub fn with_cell(&self, cell: PeriodicBox) -> Ewald {
        assert!(
            self.parameters.cutoff <= cell.half_shortest_edge(),
            "the real-space cutoff must be at most half the box's shortest edge"
        );
        let integers: Vec<[i32; 3]> = self.waves.iter().map(|w| w.n).collect();
        let mut ewald = Ewald {
            cell,
            parameters: self.parameters,
            waves: Vec::new(),
            extent: self.extent,
        };
        ewald.waves = ewald.waves_for(&integers);
        ewald
    }

    /// The box.
    pub fn cell(&self) -> PeriodicBox {
        self.cell
    }

    /// α, `r_c` and `k_c`.
    pub fn parameters(&self) -> EwaldParameters {
        self.parameters
    }

    /// How many wave vectors of the half space are summed.
    pub fn wave_vectors(&self) -> usize {
        self.waves.len()
    }

    /// `E_self = −k_e (α/√π) Σ q_i²`, joules.
    pub fn self_energy(&self, charges: &[f64]) -> f64 {
        let q2: f64 = charges.iter().map(|q| q * q).sum();
        -COULOMB * self.parameters.alpha * 0.5 * FRAC_2_SQRT_PI * q2
    }

    /// `E_bg = −k_e π Q_net² / (2 V α²)`, joules, for a net charge `net` (e).
    pub fn background(&self, net: f64) -> f64 {
        let a = self.parameters.alpha;
        -COULOMB * PI * net * net / (2.0 * self.cell.volume() * a * a)
    }

    /// A real-space pair of charge product `qq` (e²) at `r` metres: the energy
    /// `k_e qq erfc(αr)/r` and its `dE/dr`. Zero at and beyond the cutoff is the caller's to apply.
    pub fn real_pair(&self, qq: f64, r: f64) -> (f64, f64) {
        let a = self.parameters.alpha;
        let f = error_functions(a * r);
        let e = COULOMB * qq * f.erfc / r;
        (e, -(e + COULOMB * qq * a * FRAC_2_SQRT_PI * f.gaussian) / r)
    }

    /// An excluded pair's correction, `−k_e qq erf(αr)/r`, and its `dE/dr`.
    pub fn excluded_pair(&self, qq: f64, r: f64) -> (f64, f64) {
        let a = self.parameters.alpha;
        let f = error_functions(a * r);
        let e = -COULOMB * qq * f.erf / r;
        (e, -(e + COULOMB * qq * a * FRAC_2_SQRT_PI * f.gaussian) / r)
    }

    /// The reciprocal sum, with a partition. `group` marks the atoms whose charges are split off,
    /// and `lambda` scales their cross terms with the rest: returns `(E_rest, E_cross)` with
    /// `E_rest` the sum over the unmarked charges alone and `E_cross` the cross terms at full
    /// strength, and adds to `forces` and `virial` those of `E_rest + lambda E_cross`. The marked
    /// charges' own terms are not computed. With no group, `E_rest` is the whole sum and `E_cross`
    /// zero.
    pub(crate) fn reciprocal(
        &self,
        charges: &[f64],
        group: Option<(&[bool], f64)>,
        at: &[[f64; 3]],
        mut forces: Option<&mut [[f64; 3]]>,
        virial: &mut [[f64; 3]; 3],
    ) -> (f64, f64) {
        let n = at.len();
        let l = self.cell.lengths();
        // Phase factors exp(i 2π m x_a / L_a) for m = 0..=extent_a, per axis, laid out [m][atom].
        let mut tables: [(Vec<f64>, Vec<f64>); 3] = Default::default();
        for a in 0..3 {
            let m_max = self.extent[a];
            let (re, im) = &mut tables[a];
            re.resize((m_max + 1) * n, 0.0);
            im.resize((m_max + 1) * n, 0.0);
            for (i, p) in at.iter().enumerate() {
                let f = p[a] / l[a];
                let (s, c) = (2.0 * PI * (f - f.round())).sin_cos();
                re[i] = 1.0;
                im[i] = 0.0;
                for m in 1..=m_max {
                    let (pr, pi) = (re[(m - 1) * n + i], im[(m - 1) * n + i]);
                    re[m * n + i] = pr * c - pi * s;
                    im[m * n + i] = pr * s + pi * c;
                }
            }
        }
        let phase = |a: usize, m: i32, i: usize| -> (f64, f64) {
            let (re, im) = &tables[a];
            let k = m.unsigned_abs() as usize * n + i;
            if m >= 0 {
                (re[k], im[k])
            } else {
                (re[k], -im[k])
            }
        };
        let marked = |i: usize| group.is_some_and(|(g, _)| g[i]);
        let lambda = group.map_or(0.0, |(_, l)| l);
        let mut er = vec![0.0; n];
        let mut ei = vec![0.0; n];
        let (mut rest, mut cross) = (0.0, 0.0);
        for w in &self.waves {
            let (mut sr, mut si, mut gr, mut gi) = (0.0, 0.0, 0.0, 0.0);
            for i in 0..n {
                let (xr, xi) = phase(0, w.n[0], i);
                let (yr, yi) = phase(1, w.n[1], i);
                let (zr, zi) = phase(2, w.n[2], i);
                let (pr, pi) = (xr * yr - xi * yi, xr * yi + xi * yr);
                let (r, im) = (pr * zr - pi * zi, pr * zi + pi * zr);
                er[i] = r;
                ei[i] = im;
                let q = charges[i];
                if marked(i) {
                    gr += q * r;
                    gi += q * im;
                } else {
                    sr += q * r;
                    si += q * im;
                }
            }
            let own = w.weight * (sr * sr + si * si);
            let across = 2.0 * w.weight * (sr * gr + si * gi);
            rest += own;
            cross += across;
            let e = own + lambda * across;
            let k = w.k;
            for a in 0..3 {
                for b in 0..3 {
                    let delta = if a == b { 1.0 } else { 0.0 };
                    virial[a][b] += COULOMB * e * (delta - w.stress * k[a] * k[b]);
                }
            }
            if let Some(f) = forces.as_deref_mut() {
                for i in 0..n {
                    let q = charges[i];
                    if q == 0.0 {
                        continue;
                    }
                    let (tr, ti) = if marked(i) {
                        (lambda * sr, lambda * si)
                    } else {
                        (sr + lambda * gr, si + lambda * gi)
                    };
                    let g = 2.0 * COULOMB * q * w.weight * (tr * ei[i] - ti * er[i]);
                    for a in 0..3 {
                        f[i][a] += g * k[a];
                    }
                }
            }
        }
        (COULOMB * rest, COULOMB * cross)
    }

    /// The whole sum for `charges` (e) at `at` (metres, any image), with the pairs `excluded` left
    /// out: the energy by part, the force on every atom and the virial.
    ///
    /// # Panics
    ///
    /// If `charges` and `at` differ in length, or an excluded pair names an atom past the end.
    pub fn evaluate(
        &self,
        charges: &[f64],
        at: &[[f64; 3]],
        excluded: &[[usize; 2]],
    ) -> EwaldEvaluation {
        assert_eq!(charges.len(), at.len(), "one charge per position");
        let n = at.len();
        let exclusions = Exclusions::from_pairs(n, excluded);
        let mut forces = vec![[0.0; 3]; n];
        let mut virial = [[0.0; 3]; 3];
        let mut e = EwaldEnergy::default();
        let rc2 = self.parameters.cutoff * self.parameters.cutoff;
        for_each_pair(&self.cell, at, self.parameters.cutoff, |i, j, d, r2| {
            if r2 >= rc2 || exclusions.contains(i, j) {
                return;
            }
            let qq = charges[i] * charges[j];
            if qq == 0.0 {
                return;
            }
            let r = r2.sqrt();
            let (en, de_dr) = self.real_pair(qq, r);
            e.real += en;
            e.total += en;
            add_pair(&mut forces, &mut virial, i, j, d, de_dr / r);
        });
        for &[i, j] in exclusions.pairs() {
            let qq = charges[i] * charges[j];
            if qq == 0.0 {
                continue;
            }
            let d = self.cell.minimum_image(sub(at[i], at[j]));
            let r = norm(d);
            let (en, de_dr) = self.excluded_pair(qq, r);
            e.excluded += en;
            e.total += en;
            add_pair(&mut forces, &mut virial, i, j, d, de_dr / r);
        }
        if charges.iter().any(|&q| q != 0.0) {
            let (recip, _) = self.reciprocal(charges, None, at, Some(&mut forces), &mut virial);
            e.reciprocal = recip;
            e.total += recip;
        }
        e.self_energy = self.self_energy(charges);
        e.total += e.self_energy;
        let net: f64 = charges.iter().sum();
        e.background = self.background(net);
        e.total += e.background;
        for (a, row) in virial.iter_mut().enumerate() {
            row[a] += e.background;
        }
        EwaldEvaluation {
            energy: e,
            forces,
            virial,
        }
    }
}

/// Adds a pair term's force and virial: `de_dr_over_r` is `(dE/dr)/r`, `d` is `r_i − r_j`.
pub(crate) fn add_pair(
    forces: &mut [[f64; 3]],
    virial: &mut [[f64; 3]; 3],
    i: usize,
    j: usize,
    d: [f64; 3],
    de_dr_over_r: f64,
) {
    let g = -de_dr_over_r;
    for a in 0..3 {
        forces[i][a] += g * d[a];
        forces[j][a] -= g * d[a];
        for b in 0..3 {
            virial[a][b] += g * d[a] * d[b];
        }
    }
}

pub(crate) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub(crate) fn norm(d: [f64; 3]) -> f64 {
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}
