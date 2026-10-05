//! Implicit solvent: the generalized Born (GB) electrostatic solvation energy, with its analytic
//! force, and a solvent-accessible surface for the nonpolar part.
//!
//! > A. Onufriev, D. Bashford and D. A. Case, "Exploring protein native states and large-scale
//! > conformational changes with a modified generalized Born model", *Proteins* **55**, 383–394
//! > (2004). [doi:10.1002/prot.20033](https://doi.org/10.1002/prot.20033)
//!
//! Page and equation numbers below are that paper's ("OBC"), read off its pages. Everything not in
//! it — the closed form of the descreening integral, the scale factors, the radii — is marked
//! **secondary** where it is stated, with where it was taken from.
//!
//! # The energy (eqs 1–3, p. 384)
//!
//! The solvation free energy is split as `ΔG_solv = ΔG_el + ΔG_surf` (eq 1): an electrostatic part,
//! the work of charging the molecule in solvent less the work in vacuum, and a nonpolar part for
//! the uncharged molecule. GB approximates the first by (eq 2)
//!
//! `ΔG_GB = −½ Σ_ij q_i q_j / f_ij · (1 − e^(−κ f_ij) / ε_w)`,
//!
//! **over every ordered pair `i, j`, `i = j` included** — the `i = j` terms are each atom's Born
//! energy, the rest the change in the pair's Coulomb interaction — with Still's (eq 3)
//!
//! `f_ij = [r_ij² + R_i R_j exp(−r_ij² / 4 R_i R_j)]^½`,
//!
//! so `f_ii = R_i`, and `f → r` far apart. The paper's units are Gaussian; the energy here is that
//! times [`COULOMB_KCAL`] = 332.0637 kcal mol⁻¹ Å e⁻² (the crate's Coulomb constant), in SI.
//! **The solute's dielectric is 1** (p. 384: "the interior of the atom is assumed to be filled
//! uniformly with material of dielectric constant 1"), the same vacuum UFF's Coulomb term is
//! computed in, so the two add without double counting. **The solvent's is
//! [`WATER_DIELECTRIC`] = 80**, the paper's "80 for water at 300 K" (p. 384); 78.5 at 298 K would
//! change `1 − 1/ε` by 2.4 × 10⁻⁴ of itself. Salt enters as the Debye–Hückel `κ` of eq 2, with
//! `κ [Å⁻¹] ≈ 0.316 √[salt] [mol/L]` (p. 384, [`debye_kappa`]); the default is no salt, κ = 0.
//!
//! # The Born radii (eqs 4–8, p. 385)
//!
//! Each atom is a sphere of intrinsic radius ρ, used reduced by 0.09 Å, `ρ̃ = ρ − 0.09 Å`
//! (p. 385, [`RADIUS_OFFSET`]). The Hawkins–Cramer–Truhlar (HCT) radius is (eq 4)
//! `R⁻¹ = ρ̃⁻¹ − I`, with (eq 5) `I = (1/4π) ∫_VDW θ(|r| − ρ̃) r⁻⁴ d³r` over the other atoms'
//! spheres outside atom i's own. OBC rescales it (eq 6):
//!
//! `R⁻¹ = ρ̃⁻¹ − ρ⁻¹ tanh(αΨ − βΨ² + γΨ³)`, `Ψ = I ρ̃`,
//!
//! with OBC I `(α, β, γ) = (0.8, 0, 2.91)` (eq 7) and **OBC II `(1.0, 0.8, 4.85)` (eq 8), which
//! this module uses**. Since `tanh < 1`, `R < (ρ̃⁻¹ − ρ⁻¹)⁻¹` — "about 30 Å for a carbon atom"
//! (p. 386) — and R is always positive. HCT's eq 4 has no such bound: it diverges when `I`
//! reaches `ρ̃⁻¹` (p. 386), and [`Rescaling::Hct`] then gives a NaN radius rather than a
//! negative one.
//!
//! # The descreening integral (secondary)
//!
//! The paper says eq 5 is done pairwise over spheres "outside of the atom i" with a closed form in
//! Hawkins et al. and Schaefer–Froemmel (p. 385), and does not print it; nor does it print the
//! factors `S_j` by which each other atom's sphere is shrunk, to `S_j ρ̃_j`, to allow for overlap.
//! The form used here is OpenMM's `customgbforces.py` (its `I` computed value), and **it is
//! re-derived here and not taken on trust**: with `s = S_j ρ̃_j`, `r` the distance and atom i at the
//! origin, a shell of radius `u` about i lies inside sphere j over a fraction
//! `(s² − (r − u)²) / (4 r u)` of its area for `|r − s| < u < r + s`, all of it for `u < s − r`,
//! and none for `u < r − s`. So `I = ∫ frac(u) u⁻² du` from `max(ρ̃_i, |r − s|)` = L to `r + s` = U:
//!
//! `I = ½ (1/L − 1/U) + (s² − r²)/(8r) (1/L² − 1/U²) − ln(U/L)/(4r)`,
//!
//! zero when `r + s ≤ ρ̃_i`. **When atom i's sphere is wholly inside j's** — `ρ̃_i < s − r` — the
//! shells from `ρ̃_i` to `s − r` lie entirely in sphere j and add `1/ρ̃_i − 1/(s − r)`. OpenMM's
//! expression leaves that term out (its `L = max(ρ̃_i, |r − s|)` starts the integral at `s − r`);
//! it is in this one ([`descreening`]). Both forms were checked against a
//! quadrature of eq 5 — `tests/the_generalized_born_against_closed_forms.rs` does it again with its
//! own quadrature — and the inside term is what makes a probe inside a large sphere get the
//! Coulomb-field radius exactly (below).
//!
//! **Its derivative in `r` is only the explicit one**, at fixed L and U:
//! `dI/dr = −(s²/(8r²) + ⅛)(1/L² − 1/U²) + ln(U/L)/(4r²)`. The limits move with r, but the
//! integrand at `U = r + s` and at `L = |r − s|` is zero; and inside, the moving `L = s − r` and the
//! inside term's `−1/(s − r)` cancel exactly, since the integrand is one there.
//!
//! **The scale factors and radii are secondary**, read from OpenMM's `customgbforces.py` (the copy
//! the research saved; line numbers are that file's). `S`: H 0.85, C 0.72, N 0.79, O 0.85,
//! F 0.88, P 0.86, S 0.96 (`_SCREEN_PARAMETERS`' first column, lines 396–402), and its default 0.8
//! (line 403) for Cl, Br and I — the AMBER values from Hawkins, Cramer and Truhlar's fit. ρ: the
//! paper fitted α, β, γ "Based on the Bondi Radii Set" (Table I's caption, p. 387), and its one
//! OBC II trajectory is "parm99, GB^OBC (II), Bondi" (Table II, p. 388); OpenMM's `_bondi_radii`
//! (lines 238–251) is AMBER's version of that set — H 1.2, C 1.7, N 1.55, O 1.5, F 1.5, P 1.85,
//! S 1.8, Cl 1.7 Å, and 1.5 for anything else (line 240). **This is not what AMBER and OpenMM use
//! for OBC II by default**: OpenMM's `GBSAOBC2Force` ("equivalent to Amber igb=5", line 652)
//! inherits `GBSAOBC1Force`'s parameters, which are `_mbondi2_radii` (line 643) — Bondi's set with
//! a hydrogen on nitrogen at 1.3 Å, which the paper itself defines (Table II, p. 388) and pairs only
//! with OBC I there. The set the OBC II constants were fitted with is used here, so a hydrogen on
//! nitrogen is 1.2 Å; mbondi2 would make it 1.3. And **AMBER's "Bondi" is not Bondi's own
//! table**, which (J. Phys. Chem. **68**, 441 (1964), not read here; as tabulated secondarily) has
//! O 1.52, F 1.47, P 1.80 and Cl 1.75. The set the parameters were fitted with is used. **AMBER's
//! set has nothing for Br and I** — it gives them a default 1.5 Å, smaller than carbon — and
//! Bondi's 1.85 and 1.98 Å are used instead, secondary, and outside anything OBC was fitted to.
//! See [`intrinsic_radius`].
//!
//! # The force
//!
//! Analytic, with the charges fixed, by the chain rule through the Born radii: the explicit
//! derivative of eq 2 at fixed R, plus `Σ_k ∂G/∂R_k · dR_k/dI_k · ∂I_k/∂r`, where
//! `dR/dI = R² ρ̃/ρ · sech²(t) · (α − 2βΨ + 3γΨ²)` for OBC and `R²` for HCT, and
//! `∂f/∂r = r (1 − E/4)/f`, `∂f/∂R_i = R_j E (1 + r²/(4 R_i R_j)) / (2f)`, `E = e^(−r²/4R_iR_j)`.
//! Checked against central differences of the energy in
//! `tests/the_generalized_born_against_closed_forms.rs`.
//!
//! # The nonpolar part (p. 392), and why it is not in the force field
//!
//! `ΔG_surf = 0.005 kcal mol⁻¹ Å⁻² × A` (p. 392), A the solvent-accessible surface.
//! [`surface_area`] computes A by Shrake and Rupley's method (*J. Mol. Biol.* **79**, 351 (1973)):
//! points on each atom's sphere grown
//! by the probe, kept if no other grown sphere contains them, on a **golden-spiral point set**
//! — fixed, no randomness. **The probe radius is not in the paper**; 1.4 Å ([`PROBE_RADIUS`]) is
//! the conventional water probe, a choice. The count is a step function of the positions, so the
//! area has no useful gradient, and **it is reported, not put in [`ForceField`]**: an energy
//! whose force is not its gradient is what `tests/forces_are_the_gradient.rs` exists to refuse.
//! **Its convergence is measured**, against the exact area of two overlapping spheres over eight
//! pairs: the worst error is 1.86%, 0.49%, 0.29%, 0.096% and 0.025% at 100, 400, 1600, 6400 and
//! 25 600 points per atom. **No rate is claimed**: the ratios per fourfold step, 3.8, 1.7, 3.0 and
//! 3.8, are not a clean power of N, as a lattice meeting a cap's edge need not be.
//!
//! # What the numbers mean on real molecules
//!
//! With QEq charges at a UFF minimum, OBC II **over-solvates**: against FreeSolv's experimental
//! hydration free energies for eighteen molecules (`tests/the_generalized_born_on_small_molecules.rs`)
//! the polar plus nonpolar total is off by a mean of −2.1 and an RMS of 3.4 kcal/mol, the ethers
//! and anisole by about −4 to −6. QEq's charges are not the AMBER charges OBC was fitted with. So
//! **absolute solvation and binding energies built on this are not quantitative**; reported, not
//! asserted.
//!
//! # The cost, and why it is exact
//!
//! Evaluated directly, eq 2 with its force costs `N²` pair integrals for the radii, `N(N − 1)/2`
//! pair terms, and the `N²` pair integrals again for the chain rule, each integral with a
//! logarithm and each pair term with an exponential. [`GeneralizedBorn::accumulate`] computes the
//! same numbers with three changes, **none of which changes a bit**:
//!
//! - **each pair integral once**: [`descreening`] gives `I` and `dI/dr` together, and the chain
//!   rule reads the `dI/dr` the radii were computed with instead of computing the integral again;
//! - **the frozen atoms' descreening kept** ([`GeneralizedBorn::with_frozen`]): the integral
//!   between two atoms that do not move depends only on where they are, so it is computed once.
//!   What cannot be kept is said there: the radii, which every mobile atom changes, and so the
//!   pair terms between frozen atoms, which depend on them;
//! - **no salt, no `e^(−κf)`**: with κ = 0 it is `e^0`, which is exactly one, and the
//!   exponential is not called.
//!
//! Each kept or reused value is the very value the direct evaluation computes, and every sum takes
//! the same terms in the same order, so the energy, every radius and every force are the direct
//! evaluation's to the bit. The direct evaluation is kept as
//! [`GeneralizedBorn::accumulate_reference`] and [`GeneralizedBorn::born_radii_reference`],
//! and `tests/the_generalized_born_is_its_direct_sum.rs` holds the two equal, with no tolerance.
//! **There is no cutoff**: the descreening falls as `r⁻⁴` and the pair terms as `r⁻¹`, and a
//! truncated sum would be an approximation, which this is not.
//!
//! **Measured** (release, one core, `x86_64-pc-windows-gnu`), on step 3b's 987-atom complex —
//! benzene in 181L, a 6 Å zone of 322 mobile atoms in a 10 Å binding:
//!
//! | | direct | now |
//! | --- | --- | --- |
//! | one evaluation | 86.5 ms | 43.5–44.3 ms (57.6 without the kept terms) |
//! | the radii | 31.2 ms | 18.5 ms: 532 000 integrals, 35 ns each |
//! | the pair terms | 23.0 ms | 23.4 ms: 486 000 pairs, 48 ns each |
//! | the chain rule | 30.8 ms | 2.3 ms |
//! | a step of dynamics | 87.2 ms | 44.8–45.5 ms, 28 times vacuum's 1.6 ms |
//!
//! **What is left is the platform's `exp` and `ln`**: on this machine they cost 33 and 18.5 ns a
//! call, and one of each is in every pair term and pair integral. With both replaced by arithmetic
//! of no accuracy, for the timing alone, the step took 12.2 ms. A faster exponential and logarithm
//! would change the bits, and a cutoff would change the model; neither is done here.
//!
//! # Determinism
//!
//! `exp`, `ln` and `tanh` (and `sin`/`cos` for the surface's points) are the platform's, so an
//! energy repeats bit for bit on one machine and not necessarily across platforms — the crate's
//! promise.
//!
//! [`ForceField`]: crate::energy::ForceField

use crate::ccd::{Element, ANGSTROM};
use crate::energy::COULOMB_KCAL;
use crate::uff::KCAL_PER_MOL;

/// The reduction of an intrinsic radius, `ρ̃ = ρ − 0.09 Å` (p. 385), in metres.
pub const RADIUS_OFFSET: f64 = 0.09 * ANGSTROM;

/// The solvent's dielectric constant: "80 for water at 300 K" (p. 384).
pub const WATER_DIELECTRIC: f64 = 80.0;

/// The nonpolar surface tension of p. 392: 0.005 kcal mol⁻¹ Å⁻².
pub const SURFACE_TENSION_KCAL_PER_ANGSTROM2: f64 = 0.005;

/// The solvent probe radius for [`surface_area`], metres: 1.4 Å, the conventional water probe.
/// **A choice**: the OBC paper says "solvent-accessible surface" (p. 392) without a probe.
pub const PROBE_RADIUS: f64 = 1.4 * ANGSTROM;

/// The Debye–Hückel screening constant for a monovalent salt at `molar` mol/L, in m⁻¹:
/// `κ [Å⁻¹] ≈ 0.316 √[salt]` (p. 384).
pub fn debye_kappa(molar: f64) -> f64 {
    0.316 * molar.sqrt() / ANGSTROM
}

/// How the Born radius is obtained from the descreening integral `I`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Rescaling {
    /// HCT, eq 4: `R⁻¹ = ρ̃⁻¹ − I`. NaN once `I ≥ ρ̃⁻¹`, where eq 4 diverges (p. 386).
    Hct,
    /// OBC, eq 6: `R⁻¹ = ρ̃⁻¹ − ρ⁻¹ tanh(αΨ − βΨ² + γΨ³)`, `Ψ = I ρ̃`.
    Obc {
        /// α.
        alpha: f64,
        /// β.
        beta: f64,
        /// γ.
        gamma: f64,
    },
}

impl Rescaling {
    /// GB^OBC I, eq 7 (p. 385): α = 0.8, β = 0, γ = 2.91.
    pub const OBC_I: Rescaling = Rescaling::Obc {
        alpha: 0.8,
        beta: 0.0,
        gamma: 2.91,
    };

    /// GB^OBC II, eq 8 (p. 385): α = 1.0, β = 0.8, γ = 4.85. The default.
    pub const OBC_II: Rescaling = Rescaling::Obc {
        alpha: 1.0,
        beta: 0.8,
        gamma: 4.85,
    };

    /// The Born radius and `dR/dI` for intrinsic radius `rho`, reduced radius `rho_tilde` and
    /// integral `i`, in any one length unit (I in its inverse).
    fn radius(self, rho: f64, rho_tilde: f64, i: f64) -> (f64, f64) {
        match self {
            Rescaling::Hct => {
                let inverse = 1.0 / rho_tilde - i;
                if inverse > 0.0 {
                    let r = 1.0 / inverse;
                    (r, r * r)
                } else {
                    (f64::NAN, f64::NAN)
                }
            }
            Rescaling::Obc { alpha, beta, gamma } => {
                let psi = i * rho_tilde;
                let t = alpha * psi - beta * psi * psi + gamma * psi * psi * psi;
                let th = t.tanh();
                let r = 1.0 / (1.0 / rho_tilde - th / rho);
                let dt = alpha - 2.0 * beta * psi + 3.0 * gamma * psi * psi;
                (r, r * r * (1.0 - th * th) * dt * rho_tilde / rho)
            }
        }
    }
}

/// The intrinsic radius ρ of `element`, metres, before [`RADIUS_OFFSET`]. **Secondary**: AMBER's
/// "Bondi" set as OpenMM's `customgbforces.py` `_bondi_radii` holds it (lines 242–251) — the set
/// OBC's parameters were fitted with (Table I, p. 387) — H 1.2, C 1.7, N 1.55, O 1.5, F 1.5,
/// P 1.85, S 1.8, Cl 1.7 Å; and for **Br and I, which that set does not have** (it would give
/// them its default 1.5 Å, line 240), Bondi's 1964 values 1.85 and 1.98 Å, not read in primary
/// and not covered by OBC's fit. Not mbondi2: a hydrogen on nitrogen is 1.2 Å here. See the
/// module documentation.
pub fn intrinsic_radius(element: Element) -> f64 {
    let angstrom = match element {
        Element::H => 1.2,
        Element::C => 1.7,
        Element::N => 1.55,
        Element::O => 1.5,
        Element::F => 1.5,
        Element::P => 1.85,
        Element::S => 1.8,
        Element::Cl => 1.7,
        Element::Br => 1.85,
        Element::I => 1.98,
    };
    angstrom * ANGSTROM
}

/// HCT's descreening scale factor `S` of `element`: a neighbour's sphere enters the integral at
/// `S ρ̃`. **Secondary**, OpenMM's `customgbforces.py` `_SCREEN_PARAMETERS`, first column
/// (lines 396–402): H 0.85, C 0.72, N 0.79, O 0.85, F 0.88, P 0.86, S 0.96, and its default 0.8
/// (line 403) for Cl, Br and I.
pub fn scale_factor(element: Element) -> f64 {
    match element {
        Element::H => 0.85,
        Element::C => 0.72,
        Element::N => 0.79,
        Element::O => 0.85,
        Element::F => 0.88,
        Element::P => 0.86,
        Element::S => 0.96,
        Element::Cl | Element::Br | Element::I => 0.8,
    }
}

/// The pairwise descreening integral of eq 5 and its derivative in `r`: `(I, dI/dr)` for a sphere
/// of radius `s` whose centre is `r` from an atom of reduced radius `rho_tilde`, in any one length
/// unit. The closed form, the inside-sphere term and the derivative are derived in the module
/// documentation. Zero when the sphere lies inside `rho_tilde`; at `r = 0`, `1/ρ̃ − 1/s` (or zero)
/// with zero slope, where the general form would divide zero by zero.
pub fn descreening(r: f64, rho_tilde: f64, s: f64) -> (f64, f64) {
    let upper = r + s;
    if upper <= rho_tilde {
        return (0.0, 0.0);
    }
    if r == 0.0 {
        return (1.0 / rho_tilde - 1.0 / s, 0.0);
    }
    let lower = rho_tilde.max((r - s).abs());
    let (il, iu) = (1.0 / lower, 1.0 / upper);
    let squares = il * il - iu * iu;
    let log = (upper / lower).ln();
    let mut i = 0.5 * (il - iu) + (s * s - r * r) / (8.0 * r) * squares - log / (4.0 * r);
    if rho_tilde < s - r {
        i += 1.0 / rho_tilde - 1.0 / (s - r);
    }
    let di = -(s * s / (8.0 * r * r) + 0.125) * squares + log / (4.0 * r * r);
    (i, di)
}

/// Still's `f_GB` of eq 3, `[r² + R_i R_j exp(−r²/4R_iR_j)]^½`, in the unit of its arguments.
/// The same function [`GeneralizedBorn::energy`] evaluates.
pub fn still_distance(r: f64, r_i: f64, r_j: f64) -> f64 {
    still(r, r_i * r_j).0
}

/// Eq 3 for a pair at distance `r` with `R_i R_j = rr`: `(f, E)`, `E = exp(−r²/4R_iR_j)` being
/// what the derivatives also need. The one place eq 3 is written.
fn still(r: f64, rr: f64) -> (f64, f64) {
    let e = (-r * r / (4.0 * rr)).exp();
    ((r * r + rr * e).sqrt(), e)
}

/// A generalized Born model of one molecule: each atom's intrinsic radius and scale factor, the
/// rescaling, the solvent's dielectric and the salt. [`GeneralizedBorn::new`] is OBC II in water
/// with no salt.
///
/// Two models are equal when their radii, scales, rescaling, dielectric and κ are: the frozen
/// atoms' descreening that [`GeneralizedBorn::with_frozen`] keeps changes no result, and is not
/// compared.
#[derive(Clone, Debug)]
pub struct GeneralizedBorn {
    radii: Vec<f64>,
    scales: Vec<f64>,
    rescaling: Rescaling,
    solvent_dielectric: f64,
    kappa: f64,
    frozen: Option<FrozenPairs>,
}

impl PartialEq for GeneralizedBorn {
    fn eq(&self, other: &GeneralizedBorn) -> bool {
        self.radii == other.radii
            && self.scales == other.scales
            && self.rescaling == other.rescaling
            && self.solvent_dielectric == other.solvent_dielectric
            && self.kappa == other.kappa
    }
}

/// The descreening terms between two frozen atoms, computed once at the positions they are held
/// at: see [`GeneralizedBorn::with_frozen`].
#[derive(Clone)]
struct FrozenPairs {
    /// Each atom's place among the frozen atoms, or `usize::MAX` for a mobile one.
    rank: Vec<usize>,
    /// The frozen atoms, in atom order.
    atoms: Vec<usize>,
    /// Their positions' bits when the terms were computed, by rank.
    at: Vec<[u64; 3]>,
    /// `descreening(r_ij, ρ̃_i, s_j)`'s `I` for frozen `i` and `j`, at `rank_i · F + rank_j`;
    /// the diagonal is unused.
    integral: Vec<f64>,
    /// The same pairs' `dI/dr`.
    slope: Vec<f64>,
}

impl std::fmt::Debug for FrozenPairs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "FrozenPairs({} frozen atoms)", self.atoms.len())
    }
}

impl FrozenPairs {
    /// Whether every frozen atom of `at` is, to the bit, where the terms were computed.
    fn hold_at(&self, at: &[[f64; 3]]) -> bool {
        self.atoms
            .iter()
            .zip(&self.at)
            .all(|(&k, bits)| at[k].map(f64::to_bits) == *bits)
    }
}

impl GeneralizedBorn {
    /// OBC II ([`Rescaling::OBC_II`]) for atoms `elements`, with [`intrinsic_radius`] and
    /// [`scale_factor`], in water ([`WATER_DIELECTRIC`]) with no salt.
    pub fn new(elements: &[Element]) -> GeneralizedBorn {
        GeneralizedBorn::from_radii(
            elements.iter().map(|&e| intrinsic_radius(e)).collect(),
            elements.iter().map(|&e| scale_factor(e)).collect(),
        )
    }

    /// OBC II in water with no salt, for atoms of intrinsic radii `radii` (metres, before the
    /// 0.09 Å offset) and descreening scale factors `scales`.
    ///
    /// # Panics
    ///
    /// If the two are not the same length, a radius is not larger than the offset and finite, or a
    /// scale is negative or not finite.
    pub fn from_radii(radii: Vec<f64>, scales: Vec<f64>) -> GeneralizedBorn {
        assert_eq!(radii.len(), scales.len(), "one scale factor per radius");
        assert!(
            radii.iter().all(|&r| r > RADIUS_OFFSET && r.is_finite()),
            "every intrinsic radius must be finite and larger than the 0.09 Å offset"
        );
        assert!(
            scales.iter().all(|&s| s >= 0.0 && s.is_finite()),
            "every scale factor must be finite and non-negative"
        );
        GeneralizedBorn {
            radii,
            scales,
            rescaling: Rescaling::OBC_II,
            solvent_dielectric: WATER_DIELECTRIC,
            kappa: 0.0,
            frozen: None,
        }
    }

    /// The same model, told that the atoms `frozen` marks are held at their positions in `at`
    /// (metres): the descreening term of every ordered pair of two frozen atoms — `I` and `dI/dr`
    /// of eq 5's pair integral, which depend only on the two positions, ρ̃ and S — is computed
    /// here, once, and reused by every later evaluation instead of being computed again. **It is
    /// exact: every result is the same bits as without it**, because each reused term is the very
    /// value the evaluation would compute, and it enters every sum at the same place in the same
    /// order. The Born radii themselves are not kept: a frozen atom's radius changes whenever a
    /// mobile atom moves, so it is recomputed every time, from the kept terms and the others.
    ///
    /// **An evaluation at positions where a frozen atom is not, to the bit, where it was here
    /// falls back to computing every term**, so a stale term is never used; it is then only as
    /// slow as before. What it saves is the frozen–frozen share of the descreening, `F²` of the
    /// `N²` pair integrals for `F` frozen atoms of `N`.
    ///
    /// # Panics
    ///
    /// If `frozen` or `at` is not one per atom, or a frozen atom's position is not finite.
    pub fn with_frozen(mut self, frozen: &[bool], at: &[[f64; 3]]) -> GeneralizedBorn {
        let n = self.radii.len();
        assert_eq!(frozen.len(), n, "one frozen flag per atom");
        assert_eq!(at.len(), n, "one position per atom");
        let atoms: Vec<usize> = (0..n).filter(|&k| frozen[k]).collect();
        assert!(
            atoms.iter().all(|&k| at[k].iter().all(|x| x.is_finite())),
            "every frozen atom's position must be finite"
        );
        let mut rank = vec![usize::MAX; n];
        for (r, &k) in atoms.iter().enumerate() {
            rank[k] = r;
        }
        let f = atoms.len();
        let mut integral = vec![0.0; f * f];
        let mut slope = vec![0.0; f * f];
        let screened = self.screened_radii();
        for (ri, &i) in atoms.iter().enumerate() {
            let rho_tilde = self.radii[i] - RADIUS_OFFSET;
            for (rj, &j) in atoms.iter().enumerate() {
                if j != i {
                    let (t, d) = descreening(separation(at, i, j).1, rho_tilde, screened[j]);
                    integral[ri * f + rj] = t;
                    slope[ri * f + rj] = d;
                }
            }
        }
        self.frozen = Some(FrozenPairs {
            rank,
            at: atoms.iter().map(|&k| at[k].map(f64::to_bits)).collect(),
            atoms,
            integral,
            slope,
        });
        self
    }

    /// The same model without the frozen atoms' kept terms: every evaluation computes every term.
    pub fn without_frozen(mut self) -> GeneralizedBorn {
        self.frozen = None;
        self
    }

    /// How many atoms [`GeneralizedBorn::with_frozen`] was told are frozen; zero without it.
    pub fn frozen_count(&self) -> usize {
        self.frozen.as_ref().map_or(0, |f| f.atoms.len())
    }

    /// Whether an evaluation at `at` reuses the frozen atoms' kept terms: the model was given them
    /// ([`GeneralizedBorn::with_frozen`]), and every frozen atom of `at` is where it was then, to
    /// the bit. When it is not, every term is computed: the same result, more slowly.
    ///
    /// # Panics
    ///
    /// If `at` is not one position per atom.
    pub fn reuses_frozen_terms(&self, at: &[[f64; 3]]) -> bool {
        assert_eq!(at.len(), self.radii.len(), "one position per atom");
        self.frozen.as_ref().is_some_and(|f| f.hold_at(at))
    }

    /// Each atom's descreening sphere, `S ρ̃`: the expression the pair integral is called with.
    fn screened_radii(&self) -> Vec<f64> {
        self.scales
            .iter()
            .zip(&self.radii)
            .map(|(s, r)| s * (r - RADIUS_OFFSET))
            .collect()
    }

    /// The same model with the Born radii from `rescaling` instead.
    pub fn with_rescaling(mut self, rescaling: Rescaling) -> GeneralizedBorn {
        self.rescaling = rescaling;
        self
    }

    /// The same model in a solvent of dielectric constant `epsilon`.
    ///
    /// # Panics
    ///
    /// If `epsilon` is not at least 1 and finite.
    pub fn with_solvent_dielectric(mut self, epsilon: f64) -> GeneralizedBorn {
        assert!(
            epsilon >= 1.0 && epsilon.is_finite(),
            "a dielectric constant is at least 1"
        );
        self.solvent_dielectric = epsilon;
        self
    }

    /// The same model with Debye–Hückel screening `kappa`, m⁻¹ — [`debye_kappa`] gives it from a
    /// salt concentration.
    ///
    /// # Panics
    ///
    /// If `kappa` is negative or not finite.
    pub fn with_kappa(mut self, kappa: f64) -> GeneralizedBorn {
        assert!(
            kappa >= 0.0 && kappa.is_finite(),
            "κ must be non-negative and finite"
        );
        self.kappa = kappa;
        self
    }

    /// The number of atoms.
    pub fn len(&self) -> usize {
        self.radii.len()
    }

    /// Whether it has no atoms.
    pub fn is_empty(&self) -> bool {
        self.radii.is_empty()
    }

    /// The intrinsic radii ρ, metres.
    pub fn radii(&self) -> &[f64] {
        &self.radii
    }

    /// The scale factors S.
    pub fn scales(&self) -> &[f64] {
        &self.scales
    }

    /// The rescaling.
    pub fn rescaling(&self) -> Rescaling {
        self.rescaling
    }

    /// The solvent's dielectric constant.
    pub fn solvent_dielectric(&self) -> f64 {
        self.solvent_dielectric
    }

    /// κ, m⁻¹.
    pub fn kappa(&self) -> f64 {
        self.kappa
    }

    /// Every atom's effective Born radius at positions `at` (metres), metres.
    ///
    /// # Panics
    ///
    /// If `at` is not one position per atom.
    pub fn born_radii(&self, at: &[[f64; 3]]) -> Vec<f64> {
        self.descreen(at).radius
    }

    /// [`GeneralizedBorn::born_radii`] by the direct sum, every pair integral computed where it is
    /// used and nothing kept: **the reference the evaluation is held to, bit for bit**, kept as the
    /// code that was here before the evaluation was made faster. Ignores
    /// [`GeneralizedBorn::with_frozen`].
    ///
    /// # Panics
    ///
    /// If `at` is not one position per atom.
    pub fn born_radii_reference(&self, at: &[[f64; 3]]) -> Vec<f64> {
        self.radii_and_slopes_reference(at).0
    }

    /// The Born radii, each one's `dR/dI`, and every ordered pair's `dI/dr` that the kept frozen
    /// terms do not hold, for the chain rule to read instead of computing the pair integral again.
    // The loops over `j` are written as the direct sum's are, index by index, so that the two can
    // be read against each other term for term.
    #[allow(clippy::needless_range_loop)]
    fn descreen(&self, at: &[[f64; 3]]) -> Descreened<'_> {
        let n = self.radii.len();
        assert_eq!(at.len(), n, "one position per atom");
        let frozen = self.frozen.as_ref().filter(|f| f.hold_at(at));
        let kept = frozen.map_or(0, |f| f.atoms.len());
        let screened = self.screened_radii();
        let mut radius = Vec::with_capacity(n);
        let mut slope = Vec::with_capacity(n);
        let mut pairs = Vec::with_capacity(n * n.saturating_sub(1) - kept * kept.saturating_sub(1));
        let mut starts = Vec::with_capacity(n + 1);
        for i in 0..n {
            starts.push(pairs.len());
            let rho_tilde = self.radii[i] - RADIUS_OFFSET;
            let mut integral = 0.0;
            // The same terms in the same order as the direct sum: a kept term where both atoms
            // are frozen, the pair integral where either is not.
            match frozen.and_then(|f| f.row(i)) {
                Some((f, row, _)) => {
                    for j in 0..n {
                        if j == i {
                            continue;
                        }
                        match f.rank[j] {
                            usize::MAX => {
                                let (t, d) =
                                    descreening(separation(at, i, j).1, rho_tilde, screened[j]);
                                integral += t;
                                pairs.push(d);
                            }
                            rj => integral += row[rj],
                        }
                    }
                }
                None => {
                    for j in 0..n {
                        if j != i {
                            let (t, d) =
                                descreening(separation(at, i, j).1, rho_tilde, screened[j]);
                            integral += t;
                            pairs.push(d);
                        }
                    }
                }
            }
            let (r, d) = self.rescaling.radius(self.radii[i], rho_tilde, integral);
            radius.push(r);
            slope.push(d);
        }
        starts.push(pairs.len());
        Descreened {
            radius,
            slope,
            pairs,
            starts,
            frozen,
        }
    }

    /// The Born radii and each one's `dR/dI`, by the direct sum.
    fn radii_and_slopes_reference(&self, at: &[[f64; 3]]) -> (Vec<f64>, Vec<f64>) {
        let n = self.radii.len();
        assert_eq!(at.len(), n, "one position per atom");
        let mut radius = Vec::with_capacity(n);
        let mut slope = Vec::with_capacity(n);
        for i in 0..n {
            let rho_tilde = self.radii[i] - RADIUS_OFFSET;
            let mut integral = 0.0;
            for j in 0..n {
                if j != i {
                    let s = self.scales[j] * (self.radii[j] - RADIUS_OFFSET);
                    integral += descreening(separation(at, i, j).1, rho_tilde, s).0;
                }
            }
            let (r, d) = self.rescaling.radius(self.radii[i], rho_tilde, integral);
            radius.push(r);
            slope.push(d);
        }
        (radius, slope)
    }

    /// The electrostatic solvation energy ΔG_GB of eq 2 for partial charges `charges` (e) at
    /// positions `at` (metres), joules per molecule.
    ///
    /// # Panics
    ///
    /// If `charges` or `at` is not one per atom.
    pub fn energy(&self, charges: &[f64], at: &[[f64; 3]]) -> f64 {
        let mut forces = vec![[0.0; 3]; at.len()];
        self.accumulate(charges, at, &mut forces)
    }

    /// Adds the analytic force, `−∇ΔG_GB` at fixed charges, to `forces` (newtons) and returns the
    /// energy, joules per molecule. See the module documentation for the derivative.
    ///
    /// **The same bits as [`GeneralizedBorn::accumulate_reference`]**, the direct evaluation, in
    /// the energy and in every atom's force, with or without [`GeneralizedBorn::with_frozen`]: it
    /// computes each pair integral once and keeps its `dI/dr` for the chain rule, where the
    /// direct evaluation computes it twice; it reuses the kept frozen–frozen terms; and with no
    /// salt it does not evaluate `e^(−κf)`, which is `e^0`, exactly one. Every sum takes the same
    /// terms in the same order. See the module documentation, "The cost".
    ///
    /// # Panics
    ///
    /// If `charges`, `at` or `forces` is not one per atom.
    pub fn accumulate(&self, charges: &[f64], at: &[[f64; 3]], forces: &mut [[f64; 3]]) -> f64 {
        let n = self.radii.len();
        assert_eq!(charges.len(), n, "one charge per atom");
        assert_eq!(forces.len(), n, "one force per atom");
        let descreened = self.descreen(at);
        let (born, slope) = (&descreened.radius, &descreened.slope);
        let k = COULOMB_KCAL * KCAL_PER_MOL * ANGSTROM;
        let (eps, kappa) = (self.solvent_dielectric, self.kappa);
        // g(f) = (1 − e^{−κf}/ε)/f and g′(f). Without salt e^{−κf} = e^{−0} = 1 exactly, so the
        // exponential is not called: the same bits, at a twentieth of the cost of a pair.
        let g = |f: f64| {
            let screen = if kappa == 0.0 {
                1.0 / eps
            } else {
                (-kappa * f).exp() / eps
            };
            let g = (1.0 - screen) / f;
            (g, -g / f + kappa * screen / f)
        };
        let mut energy = 0.0;
        let mut d_born = vec![0.0; n];
        for i in 0..n {
            if charges[i] == 0.0 {
                continue;
            }
            // The i = j term of eq 2: f_ii = R_i.
            let (gi, dgi) = g(born[i]);
            let c = -0.5 * k * charges[i] * charges[i];
            energy += c * gi;
            d_born[i] += c * dgi;
            for j in i + 1..n {
                if charges[j] == 0.0 {
                    continue;
                }
                let (d, r) = separation(at, i, j);
                let rr = born[i] * born[j];
                let (f, e) = still(r, rr);
                let (gf, dgf) = g(f);
                // Both orderings, i j and j i.
                let c = -k * charges[i] * charges[j];
                energy += c * gf;
                let de_df = c * dgf;
                // ∂f/∂x_i = (1 − E/4) (x_i − x_j)/f: the force −∂E/∂x_i, and its opposite on j.
                let radial = -de_df * (1.0 - 0.25 * e) / f;
                for a in 0..3 {
                    forces[i][a] += radial * d[a];
                    forces[j][a] -= radial * d[a];
                }
                let common = de_df * e * (1.0 + r * r / (4.0 * rr)) / (2.0 * f);
                d_born[i] += common * born[j];
                d_born[j] += common * born[i];
            }
        }
        // The chain rule through every Born radius: ∂E/∂R_i · dR_i/dI_i · ∂I_i/∂r_ij, with each
        // ∂I_i/∂r_ij the one the radii were computed with.
        for i in 0..n {
            let c = d_born[i] * slope[i];
            if c == 0.0 {
                continue;
            }
            let mut own = descreened.pairs[descreened.starts[i]..descreened.starts[i + 1]].iter();
            let kept = descreened.frozen.and_then(|f| f.row(i));
            for j in 0..n {
                if j == i {
                    continue;
                }
                let di = match kept {
                    Some((f, _, row)) if f.rank[j] != usize::MAX => row[f.rank[j]],
                    _ => *own
                        .next()
                        .expect("one slope per pair the kept terms do not hold"),
                };
                let (d, r) = separation(at, i, j);
                if r == 0.0 {
                    continue;
                }
                let radial = -c * di / r;
                for a in 0..3 {
                    forces[i][a] += radial * d[a];
                    forces[j][a] -= radial * d[a];
                }
            }
        }
        energy
    }

    /// [`GeneralizedBorn::accumulate`] by the direct evaluation: every pair integral computed where
    /// it is used — twice, for the radii and for the chain rule — `e^(−κf)` called for every pair,
    /// and nothing kept. **The reference the evaluation is held to, bit for bit**, kept as the code
    /// that was here before it was made faster. Ignores [`GeneralizedBorn::with_frozen`].
    ///
    /// # Panics
    ///
    /// If `charges`, `at` or `forces` is not one per atom.
    pub fn accumulate_reference(
        &self,
        charges: &[f64],
        at: &[[f64; 3]],
        forces: &mut [[f64; 3]],
    ) -> f64 {
        let n = self.radii.len();
        assert_eq!(charges.len(), n, "one charge per atom");
        assert_eq!(forces.len(), n, "one force per atom");
        let (born, slope) = self.radii_and_slopes_reference(at);
        let k = COULOMB_KCAL * KCAL_PER_MOL * ANGSTROM;
        let (eps, kappa) = (self.solvent_dielectric, self.kappa);
        // g(f) = (1 − e^{−κf}/ε)/f and g′(f).
        let g = |f: f64| {
            let screen = (-kappa * f).exp() / eps;
            let g = (1.0 - screen) / f;
            (g, -g / f + kappa * screen / f)
        };
        let mut energy = 0.0;
        let mut d_born = vec![0.0; n];
        for i in 0..n {
            if charges[i] == 0.0 {
                continue;
            }
            // The i = j term of eq 2: f_ii = R_i.
            let (gi, dgi) = g(born[i]);
            let c = -0.5 * k * charges[i] * charges[i];
            energy += c * gi;
            d_born[i] += c * dgi;
            for j in i + 1..n {
                if charges[j] == 0.0 {
                    continue;
                }
                let (d, r) = separation(at, i, j);
                let rr = born[i] * born[j];
                let (f, e) = still(r, rr);
                let (gf, dgf) = g(f);
                // Both orderings, i j and j i.
                let c = -k * charges[i] * charges[j];
                energy += c * gf;
                let de_df = c * dgf;
                // ∂f/∂x_i = (1 − E/4) (x_i − x_j)/f: the force −∂E/∂x_i, and its opposite on j.
                let radial = -de_df * (1.0 - 0.25 * e) / f;
                for a in 0..3 {
                    forces[i][a] += radial * d[a];
                    forces[j][a] -= radial * d[a];
                }
                let common = de_df * e * (1.0 + r * r / (4.0 * rr)) / (2.0 * f);
                d_born[i] += common * born[j];
                d_born[j] += common * born[i];
            }
        }
        // The chain rule through every Born radius: ∂E/∂R_i · dR_i/dI_i · ∂I_i/∂r_ij.
        for i in 0..n {
            let c = d_born[i] * slope[i];
            if c == 0.0 {
                continue;
            }
            let rho_tilde = self.radii[i] - RADIUS_OFFSET;
            for j in 0..n {
                if j == i {
                    continue;
                }
                let (d, r) = separation(at, i, j);
                if r == 0.0 {
                    continue;
                }
                let s = self.scales[j] * (self.radii[j] - RADIUS_OFFSET);
                let (_, di) = descreening(r, rho_tilde, s);
                let radial = -c * di / r;
                for a in 0..3 {
                    forces[i][a] += radial * d[a];
                    forces[j][a] -= radial * d[a];
                }
            }
        }
        energy
    }
}

/// What [`GeneralizedBorn::descreen`] computed: the radii, their slopes `dR/dI`, and the `dI/dr`
/// of every ordered pair `(i, j)`, `j ≠ i`, in the order `i` then `j`, except the pairs of two
/// frozen atoms when the kept terms held — row `i` is `pairs[starts[i]..starts[i + 1]]`.
struct Descreened<'a> {
    radius: Vec<f64>,
    slope: Vec<f64>,
    pairs: Vec<f64>,
    starts: Vec<usize>,
    frozen: Option<&'a FrozenPairs>,
}

impl FrozenPairs {
    /// For a frozen atom `i`, its row of kept integrals and of kept slopes, by rank.
    fn row(&self, i: usize) -> Option<(&FrozenPairs, &[f64], &[f64])> {
        let f = self.atoms.len();
        match self.rank[i] {
            usize::MAX => None,
            r => Some((
                self,
                &self.integral[r * f..(r + 1) * f],
                &self.slope[r * f..(r + 1) * f],
            )),
        }
    }
}

/// `r_i − r_j` and its length.
fn separation(at: &[[f64; 3]], i: usize, j: usize) -> ([f64; 3], f64) {
    let d = [
        at[i][0] - at[j][0],
        at[i][1] - at[j][1],
        at[i][2] - at[j][2],
    ];
    (d, (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
}

/// The `n` points of the golden-spiral set on the unit sphere: `z = 1 − (2k + 1)/n`, azimuth
/// `k π (3 − √5)`. Each covers an equal area, `4π/n`, to the extent the spiral does; fixed, so a
/// surface repeats.
pub fn sphere_points(n: usize) -> Vec<[f64; 3]> {
    let turn = std::f64::consts::PI * (3.0 - 5f64.sqrt());
    (0..n)
        .map(|k| {
            let z = 1.0 - (2 * k + 1) as f64 / n as f64;
            let rho = (1.0 - z * z).max(0.0).sqrt();
            let phi = k as f64 * turn;
            [rho * phi.cos(), rho * phi.sin(), z]
        })
        .collect()
}

/// Each atom's solvent-accessible surface, m²: Shrake and Rupley's count over `points` golden-spiral
/// points ([`sphere_points`]) on each sphere of radius `radii[i] + probe`, keeping those inside no
/// other atom's grown sphere. Exact for an isolated atom; for overlapping ones it converges as the
/// count grows, which the tests measure against two spheres' exact area. The points are fixed in
/// the lab frame, so a rotated molecule's area differs by the discretisation.
///
/// # Panics
///
/// If `radii` and `at` differ in length, or `points` is zero.
pub fn surface_area(radii: &[f64], at: &[[f64; 3]], probe: f64, points: usize) -> Vec<f64> {
    assert_eq!(radii.len(), at.len(), "one radius per position");
    assert!(points > 0, "at least one point per sphere");
    let unit = sphere_points(points);
    let grown: Vec<f64> = radii.iter().map(|r| r + probe).collect();
    (0..at.len())
        .map(|i| {
            let neighbours: Vec<usize> = (0..at.len())
                .filter(|&j| j != i && separation(at, i, j).1 < grown[i] + grown[j])
                .collect();
            let free = unit
                .iter()
                .filter(|p| {
                    let x = [
                        at[i][0] + grown[i] * p[0],
                        at[i][1] + grown[i] * p[1],
                        at[i][2] + grown[i] * p[2],
                    ];
                    neighbours.iter().all(|&j| {
                        let d = [x[0] - at[j][0], x[1] - at[j][1], x[2] - at[j][2]];
                        d[0] * d[0] + d[1] * d[1] + d[2] * d[2] >= grown[j] * grown[j]
                    })
                })
                .count();
            4.0 * std::f64::consts::PI * grown[i] * grown[i] * free as f64 / points as f64
        })
        .collect()
}

/// The nonpolar solvation energy of p. 392, `0.005 kcal mol⁻¹ Å⁻² × A`, for a total
/// solvent-accessible surface `area` m², joules per molecule.
pub fn nonpolar_energy(area: f64) -> f64 {
    SURFACE_TENSION_KCAL_PER_ANGSTROM2 * KCAL_PER_MOL * area / (ANGSTROM * ANGSTROM)
}
