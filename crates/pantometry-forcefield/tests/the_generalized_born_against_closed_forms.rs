//! **Generalized Born against what is known about it without it**: eq 5 by quadrature, a lone ion's
//! Born energy, two ions far apart and at one point, a charge inside a dielectric sphere — where
//! HCT's radius is the Coulomb-field radius exactly and Kirkwood's series says how far that is from
//! the true reaction field — the paper's own constants and bound, rigid motions, and the force
//! against central differences of the energy.
//!
//! Every expected value is a closed form or a number printed in Onufriev, Bashford and Case,
//! *Proteins* **55**, 383 (2004), typed here from the paper, not called from the crate.

// The force checks walk atoms and axes as the formulas index them.
#![allow(clippy::needless_range_loop)]

use pantometry_forcefield::solvation::{
    self, descreening, sphere_points, still_distance, surface_area, GeneralizedBorn, Rescaling,
    PROBE_RADIUS, RADIUS_OFFSET, WATER_DIELECTRIC,
};
use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::{Component, Element, ForceField};

const AIN: &str = include_str!("../components/AIN.cif");
const ANGSTROM: f64 = 1e-10;
const EPS: f64 = f64::EPSILON;
/// e²/(4π ε₀) in kcal mol⁻¹ Å e⁻², the paper's Gaussian units made molar.
const K: f64 = 332.0637;

fn kcal(joules: f64) -> f64 {
    joules / KCAL_PER_MOL
}

/// The born energy `−(q²/2R)(1 − 1/ε) K`, kcal/mol, R in Å.
fn born(q: f64, r: f64, eps: f64) -> f64 {
    -0.5 * K * q * q * (1.0 - 1.0 / eps) / r
}

/// Eq 5 for one sphere by its own quadrature, with no closed form: `∫ frac(u) / u² du` from
/// `rho` to `r + s`, `frac(u)` the fraction of the shell of radius `u` about the atom inside a
/// sphere of radius `s` centred `r` away — one for `u < s − r`, zero for `u < r − s`, and
/// `(s² − (r − u)²) / (4 r u)` between. Composite Simpson on each smooth piece, 20 000 panels:
/// the integrand's fourth derivative is at most a few hundred over these lengths, so the
/// truncation is below 10⁻¹⁵.
fn eq5_by_quadrature(r: f64, rho: f64, s: f64) -> f64 {
    let frac = |u: f64| {
        if u <= (r - s).abs() {
            if s > r {
                1.0
            } else {
                0.0
            }
        } else if u >= r + s {
            0.0
        } else {
            (s * s - (r - u) * (r - u)) / (4.0 * r * u)
        }
    };
    let f = |u: f64| frac(u) / (u * u);
    let mut cuts = vec![rho, r + s];
    let mid = (r - s).abs();
    if mid > rho && mid < r + s {
        cuts.insert(1, mid);
    }
    let mut total = 0.0;
    for w in cuts.windows(2) {
        let (a, b) = (w[0], w[1]);
        if b <= a {
            continue;
        }
        let n = 20_000;
        let h = (b - a) / n as f64;
        let mut sum = f(a) + f(b);
        for k in 1..n {
            let u = a + k as f64 * h;
            sum += if k % 2 == 1 { 4.0 } else { 2.0 } * f(u);
        }
        total += sum * h / 3.0;
    }
    total
}

/// Every branch of the pair integral: apart (the lower limit `r − s`, outside ρ̃), overlapping
/// ρ̃ (the lower limit ρ̃), the atom's sphere wholly inside the other (`ρ̃ < s − r`), and the other
/// wholly inside the atom's (`r + s ≤ ρ̃`). Lengths in Å.
const CASES: [(f64, f64, f64, &str); 7] = [
    (3.0, 1.5, 1.2, "apart"),
    (1.5, 1.5, 1.2, "overlapping, from ρ̃"),
    (1.0, 1.5, 1.8, "overlapping ρ̃, s > r"),
    (2.0, 0.6, 1.0, "apart, a small ρ̃"),
    (0.7, 1.0, 3.0, "atom inside the sphere"),
    (1.2, 1.11, 2.91, "atom inside the sphere"),
    (0.3, 1.5, 1.0, "sphere inside the atom"),
];

#[test]
fn the_pair_integral_is_eq_5_by_quadrature() {
    for (r, rho, s, what) in CASES {
        let (closed, slope) = descreening(r, rho, s);
        let quad = eq5_by_quadrature(r, rho, s);
        println!("{what:28} r {r} ρ̃ {rho} s {s}: I = {closed:.15} (quadrature {quad:.15})");
        assert!(
            (closed - quad).abs() <= 1e-12 * quad.abs().max(1e-3),
            "{what}: {closed} against {quad}"
        );
        // The derivative against a central difference of the quadrature, the slope's error being
        // h²/6 I‴ ≤ 1e-8 of it at h = 1e-4 Å, and the quadrature's rounding 1e-15 / h.
        let h = 1e-4;
        let fd = (eq5_by_quadrature(r + h, rho, s) - eq5_by_quadrature(r - h, rho, s)) / (2.0 * h);
        assert!(
            (slope - fd).abs() <= 1e-7 * slope.abs().max(1e-3),
            "{what}: dI/dr {slope} against {fd}"
        );
    }
    // Wholly inside the atom's own excluded sphere: nothing.
    assert_eq!(descreening(0.3, 1.5, 1.0), (0.0, 0.0));
    // Concentric and larger: 1/ρ̃ − 1/s, flat.
    assert_eq!(descreening(0.0, 1.5, 3.0), (1.0 / 1.5 - 1.0 / 3.0, 0.0));
}

/// The inside case is the one OpenMM's expression drops: check that the term it drops is what
/// separates the two, so that the comparison above is about the term.
#[test]
fn the_inside_term_is_the_shells_wholly_in_the_sphere() {
    let (r, rho, s) = (0.7, 1.0, 3.0);
    let (with, _) = descreening(r, rho, s);
    let u = r + s;
    let l = s - r;
    let without = 0.5 * (1.0 / l - 1.0 / u)
        + (s * s - r * r) / (8.0 * r) * (1.0 / (l * l) - 1.0 / (u * u))
        - (u / l).ln() / (4.0 * r);
    let shells = 1.0 / rho - 1.0 / (s - r);
    assert!((with - without - shells).abs() <= 8.0 * EPS * with);
    assert!(
        shells > 0.4 * with,
        "the inside term is most of it here: {shells} of {with}"
    );
}

#[test]
fn a_lone_ion_has_its_reduced_radius_and_borns_energy() {
    // ρ typed here, Å: see `the_radii_and_scale_factors_are_the_tables` for where each is from.
    for (element, rho) in [
        (Element::H, 1.2),
        (Element::C, 1.7),
        (Element::N, 1.55),
        (Element::O, 1.5),
        (Element::S, 1.8),
        (Element::I, 1.98),
    ] {
        let rho = rho * ANGSTROM;
        for rescaling in [Rescaling::OBC_II, Rescaling::OBC_I, Rescaling::Hct] {
            let gb = GeneralizedBorn::new(&[element]).with_rescaling(rescaling);
            let at = [[0.3 * ANGSTROM, -ANGSTROM, 2.0 * ANGSTROM]];
            let r = gb.born_radii(&at)[0];
            let reduced = rho - 0.09 * ANGSTROM;
            assert!(
                (r - reduced).abs() <= 4.0 * EPS * reduced,
                "{element:?} {rescaling:?}: R {r} against ρ − 0.09 Å = {reduced}"
            );
            for q in [1.0, -1.0, 0.37] {
                let g = kcal(gb.energy(&[q], &at));
                let expected = born(q, reduced / ANGSTROM, WATER_DIELECTRIC);
                assert!(
                    (g - expected).abs() <= 1e-13 * expected.abs(),
                    "{element:?} q {q}: {g} against Born's {expected}"
                );
            }
        }
    }
    // Salt: eq 2's i = i term with Debye–Hückel screening, `−½ K q² (1 − e^(−κR)/ε) / R`, and
    // κ = 0.316 √[salt] Å⁻¹ (p. 384) at 0.15 M.
    let kappa = 0.316 * 0.15f64.sqrt();
    let salty = GeneralizedBorn::new(&[Element::O]).with_kappa(solvation::debye_kappa(0.15));
    let g = kcal(salty.energy(&[-0.8], &[[0.0; 3]]));
    let expected = -0.5 * K * 0.64 * (1.0 - (-kappa * 1.41).exp() / 80.0) / 1.41;
    println!("one −0.8 e oxygen in 0.15 M salt: {g:.6} kcal/mol (κ = {kappa:.5} Å⁻¹)");
    assert!((g - expected).abs() <= 1e-13 * expected.abs());
    // Sodium-sized, in kcal/mol: q = 1 at ρ̃ = 1.41 Å, −116.29.
    let gb = GeneralizedBorn::from_radii(vec![1.5 * ANGSTROM], vec![0.8]);
    let g = kcal(gb.energy(&[1.0], &[[0.0; 3]]));
    println!("one charge at ρ = 1.5 Å: {g:.4} kcal/mol");
    assert!((g - born(1.0, 1.41, 80.0)).abs() < 1e-12);
}

/// Far apart, each atom's radius tends to its ρ̃ and Still's f to r, so eq 2 tends to the two Born
/// energies and `−K (1 − 1/ε) q₁q₂ / r`, the change that turns the vacuum Coulomb energy into the
/// one screened by ε. The tolerance is the size of the two things that have not yet gone:
///
/// - **The descreening.** A sphere of radius `s` at distance `r` gives `I ≤ s³ / (3 (r − s)⁴)`
///   (its volume, `4π s³/3`, over 4π, at the nearest distance any of it is). Eq 6 then moves `1/R`
///   below `1/ρ̃` by `ρ⁻¹ tanh t ≤ ρ⁻¹ (α Ψ + γ Ψ³)` (β only lowers t), `Ψ = I ρ̃`; that times the
///   Born term's `½ K (1 − 1/ε) q²` bounds each self term's change, and the cross term's through
///   `R` is smaller and inside the next bound.
/// - **Still's exponential**: `1/r − 1/f ≤ R₁R₂ e^(−r²/4R₁R₂) / (2 r³)`, with `R ≤` the radius the
///   bound on I allows.
///
/// Plus rounding, `10⁻¹³` of the sum of the magnitudes.
#[test]
fn two_ions_far_apart_are_two_born_terms_and_a_screened_coulomb() {
    let elements = [Element::O, Element::C];
    let q = [1.0, -0.6];
    let gb = GeneralizedBorn::new(&elements);
    let rho = elements.map(|e| solvation::intrinsic_radius(e) / ANGSTROM);
    let tilde = rho.map(|r| r - 0.09);
    let s = [0.85 * tilde[0], 0.72 * tilde[1]];
    let (alpha, gamma) = (1.0, 4.85);
    for r in [10.0, 20.0, 40.0] {
        let at = [[0.0; 3], [r * ANGSTROM, 0.0, 0.0]];
        let g = kcal(gb.energy(&q, &at));
        let shield = 1.0 - 1.0 / WATER_DIELECTRIC;
        let limit =
            born(q[0], tilde[0], 80.0) + born(q[1], tilde[1], 80.0) - K * shield * q[0] * q[1] / r;
        let mut tolerance = 0.0;
        let mut r_max = [0.0; 2];
        for i in 0..2 {
            let other = s[1 - i];
            let bound_i = other.powi(3) / (3.0 * (r - other).powi(4));
            let psi = bound_i * tilde[i];
            let shift = (alpha * psi + gamma * psi.powi(3)) / rho[i];
            tolerance += 0.5 * K * shield * q[i] * q[i] * shift;
            r_max[i] = 1.0 / (1.0 / tilde[i] - shift);
        }
        let rr = r_max[0] * r_max[1];
        tolerance +=
            K * shield * (q[0] * q[1]).abs() * rr * (-r * r / (4.0 * rr)).exp() / (2.0 * r.powi(3));
        let magnitudes = born(q[0], tilde[0], 80.0).abs() + born(q[1], tilde[1], 80.0).abs();
        tolerance += 1e-13 * magnitudes;
        println!(
            "r {r:4} Å: ΔG {g:.10} against the limit {limit:.10}, off by {:.3e} within {tolerance:.3e} \
             (the screened Coulomb part is {:.4})",
            g - limit,
            -K * shield * q[0] * q[1] / r
        );
        assert!((g - limit).abs() <= tolerance, "r = {r} Å");
        // The bound is small beside what it bounds, so the check can see the cross term.
        assert!(tolerance < 1e-3 * (K * shield * q[0] * q[1] / r).abs());
    }
}

/// Two charges at one point with one radius are one charge, `q₁ + q₂`: the self terms and the two
/// cross terms of eq 2 make `(q₁ + q₂)²`, and that needs Still's `f(0) = R`, which is its
/// exponential's doing.
#[test]
fn two_coincident_ions_are_one_ion_of_their_combined_charge() {
    let rho = 1.6 * ANGSTROM;
    let at = [[1e-10, 2e-10, 3e-10], [1e-10, 2e-10, 3e-10]];
    let q = [0.7, 0.5];
    // No descreening (scale 0), so R = ρ̃ exactly.
    let alone = GeneralizedBorn::from_radii(vec![rho, rho], vec![0.0, 0.0]);
    let g = kcal(alone.energy(&q, &at));
    let expected = born(1.2, 1.6 - 0.09, 80.0);
    println!("coincident, no descreening: {g:.12} against {expected:.12}");
    assert!((g - expected).abs() <= 1e-13 * expected.abs());
    // Descreening each other: still one charge at whatever radius both have.
    let gb = GeneralizedBorn::from_radii(vec![rho, rho], vec![0.8, 0.8]);
    let radii = gb.born_radii(&at);
    assert_eq!(radii[0], radii[1]);
    let g = kcal(gb.energy(&q, &at));
    let expected = born(1.2, radii[0] / ANGSTROM, 80.0);
    assert!((g - expected).abs() <= 1e-13 * expected.abs());
    assert_eq!(still_distance(0.0, 1.3, 1.3), 1.3);
}

/// Kirkwood's reaction field of a charge at distance `d` from the centre of a sphere of radius `a`,
/// dielectric 1 inside and ε outside (J. G. Kirkwood, *J. Chem. Phys.* **2**, 351 (1934)), as the
/// effective radius GB would need: `−(q²/2) (1 − 1/ε) / R = −(q²/2a) Σ_n (n+1)(ε−1)/((n+1)ε + n)
/// (d/a)^(2n)`. Summed to convergence (the terms fall as `(d/a)²`).
fn kirkwood_radius(a: f64, d: f64, eps: f64) -> f64 {
    let x = (d / a) * (d / a);
    let mut sum = 0.0;
    let mut power = 1.0;
    for n in 0..400 {
        let n = n as f64;
        let term = (n + 1.0) * (eps - 1.0) / ((n + 1.0) * eps + n) * power;
        sum += term;
        if term < 1e-18 * sum {
            break;
        }
        power *= x;
    }
    (1.0 - 1.0 / eps) / (sum / a)
}

/// **A charge inside a sphere.** HCT's eq 4 computes `R⁻¹` as the Coulomb field's energy outside
/// the solute (p. 384: "the energy density of the electrostatic field … approximated as the energy
/// density of a Coulomb field"). For one probe atom inside one large sphere that is
/// `(1/4π) ∫_{|x|>a} |x − d|⁻⁴ d³x = ½ [a/(a² − d²) + ln((a+d)/(a−d)) / (2d)]`, whatever the probe's
/// own radius, and HCT must give it exactly — **with the inside-sphere term**. At a = 3 Å,
/// d = 1.2 Å that is R = 2.66714 Å. **The exact reaction field is not that**: Kirkwood's series at
/// ε = 80 gives 2.52265 Å, so the Coulomb-field approximation's radius is 5.73% long and its
/// energy 5.4% short. That is the model's error. **Only the first half of this test checks the
/// crate**, HCT against the Coulomb-field radius. The Kirkwood half uses only this file's functions:
/// it is a documented measurement of what the model gets wrong, held so the 5.73% stated in the
/// docs stays the number the two closed forms give, and it guards no crate code.
#[test]
fn a_charge_inside_a_sphere_gets_the_coulomb_field_radius_and_misses_kirkwood() {
    let (a, d): (f64, f64) = (3.0, 1.2);
    let cfa = 1.0 / (0.5 * (a / (a * a - d * d) + ((a + d) / (a - d)).ln() / (2.0 * d)));
    println!("Coulomb-field radius {cfa:.6} Å");
    assert!(
        (cfa - 2.6671).abs() < 5e-5,
        "the closed form, against the figure worked by hand"
    );
    for probe in [1.2, 1.5, 1.7] {
        // The large sphere as one dense atom: S = 1 and ρ̃ = a.
        let gb = GeneralizedBorn::from_radii(
            vec![probe * ANGSTROM, (a + 0.09) * ANGSTROM],
            vec![0.85, 1.0],
        );
        let at = [[0.0, d * ANGSTROM, 0.0], [0.0; 3]];
        let hct = gb.clone().with_rescaling(Rescaling::Hct).born_radii(&at)[0] / ANGSTROM;
        let obc = gb.born_radii(&at)[0] / ANGSTROM;
        println!("probe ρ {probe} Å: HCT {hct:.12} Å, OBC II {obc:.6} Å");
        assert!(
            (hct - cfa).abs() <= 1e-12 * cfa,
            "HCT {hct} against the Coulomb-field radius {cfa}"
        );
    }
    // The measurement, not a guard: nothing below calls the crate.
    let exact = kirkwood_radius(a, d, 80.0);
    let gap = cfa / exact - 1.0;
    println!(
        "Kirkwood at ε = 80: {exact:.6} Å; HCT's radius is {:.3}% long, its energy {:.3}% short",
        100.0 * gap,
        100.0 * (1.0 - exact / cfa)
    );
    // The series against its own two limits: the Born radius at the centre, and a = 3 at d = 0.
    assert!((kirkwood_radius(a, 0.0, 80.0) - a).abs() < 1e-14);
    assert!((exact - 2.5227).abs() < 5e-5);
    assert!(
        (gap - 0.0573).abs() < 5e-5,
        "the documented measurement no longer matches: {gap}"
    );
}

/// Eq 6 with eq 8's constants, typed from the paper, against the module's radius on a pair whose
/// integral is the one checked by quadrature above; OBC I and HCT beside it, so that the three are
/// distinguished on this input.
#[test]
fn obc_ii_is_eq_6_with_eq_8s_constants() {
    let gb = GeneralizedBorn::new(&[Element::C, Element::O]);
    assert_eq!(gb.rescaling(), Rescaling::OBC_II);
    assert_eq!(
        Rescaling::OBC_II,
        Rescaling::Obc {
            alpha: 1.0,
            beta: 0.8,
            gamma: 4.85
        }
    );
    assert_eq!(
        Rescaling::OBC_I,
        Rescaling::Obc {
            alpha: 0.8,
            beta: 0.0,
            gamma: 2.91
        }
    );
    assert_eq!(WATER_DIELECTRIC, 80.0);
    assert_eq!(RADIUS_OFFSET, 0.09 * ANGSTROM);
    let r = 1.25;
    let at = [[0.0; 3], [r * ANGSTROM, 0.0, 0.0]];
    let (rho, tilde) = (1.7, 1.61);
    let i = eq5_by_quadrature(r, tilde, 0.85 * 1.41);
    let psi = i * tilde;
    let eq6 = |a: f64, b: f64, g: f64| {
        1.0 / (1.0 / tilde - (a * psi - b * psi * psi + g * psi.powi(3)).tanh() / rho)
    };
    let expect = [
        (Rescaling::OBC_II, eq6(1.0, 0.8, 4.85)),
        (Rescaling::OBC_I, eq6(0.8, 0.0, 2.91)),
        (Rescaling::Hct, 1.0 / (1.0 / tilde - i)),
    ];
    for (rescaling, want) in expect {
        let got = gb.clone().with_rescaling(rescaling).born_radii(&at)[0] / ANGSTROM;
        println!("{rescaling:?}: carbon's R {got:.10} Å (eq 6: {want:.10}), Ψ = {psi:.4}");
        assert!((got - want).abs() <= 1e-11 * want);
    }
    assert!(
        (expect[0].1 - expect[1].1).abs() > 1e-3,
        "the two sets are distinguished here"
    );
}

/// p. 386: "the upper limit of an R_i computed from Eq. (6) is (ρ̃⁻¹ − ρ⁻¹)⁻¹, which is about 30 Å
/// for a carbon atom (ρ = 1.7 Å)". A carbon deep inside a dense sphere 50 Å across approaches it
/// from below. HCT's eq 4 has no such bound, and once overlapping neighbours push I past `1/ρ̃` it
/// has no radius at all: NaN, not a negative number.
#[test]
fn obc_radii_are_bounded_as_the_paper_says_and_hct_ones_are_not() {
    let bound: f64 = 1.0 / (1.0 / 1.61 - 1.0 / 1.7);
    println!("OBC's bound for carbon: {bound:.4} Å");
    assert!((bound - 30.0).abs() < 0.5, "'about 30 Å'");
    let gb = GeneralizedBorn::from_radii(vec![1.7 * ANGSTROM, 50.09 * ANGSTROM], vec![0.72, 1.0]);
    let at = [[0.0; 3], [0.0, 0.0, 0.5 * ANGSTROM]];
    let r = gb.born_radii(&at)[0] / ANGSTROM;
    println!("carbon inside a 50 Å sphere: OBC II R = {r:.4} Å");
    assert!(r < bound && r > 0.99 * bound);
    // Sixty spheres of ρ 1.5 Å, S 0.8, stacked 1 Å from a carbon: each gives I = 0.011476 Å⁻¹, and
    // sixty make 0.689, past 1/ρ̃ = 0.621 — the overlap overcounting eq 4 cannot survive.
    let mut radii = vec![1.7 * ANGSTROM];
    let mut at = vec![[0.0; 3]];
    for _ in 0..60 {
        radii.push(1.5 * ANGSTROM);
        at.push([1.0 * ANGSTROM, 0.0, 0.0]);
    }
    let scales = vec![0.8; radii.len()];
    let gb = GeneralizedBorn::from_radii(radii, scales);
    assert!(gb.clone().with_rescaling(Rescaling::Hct).born_radii(&at)[0].is_nan());
    let obc = gb.born_radii(&at)[0] / ANGSTROM;
    assert!(obc > 0.0 && obc < bound, "OBC stays bounded: {obc}");
}

/// Aspirin with deterministic neutral charges up to ±0.3 e, moved off the dictionary geometry by
/// up to 0.05 Å per coordinate.
fn aspirin() -> (Component, GeneralizedBorn, Vec<f64>, Vec<[f64; 3]>) {
    let c = Component::from_ccd(AIN).expect("AIN parses");
    let n = c.atoms().len();
    let raw: Vec<f64> = (0..n).map(|i| 0.3 * (1.7 * i as f64).cos()).collect();
    let mean = raw.iter().sum::<f64>() / n as f64;
    let charges = raw.iter().map(|q| q - mean).collect();
    let elements: Vec<Element> = c.atoms().iter().map(|a| a.element).collect();
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
    (c, GeneralizedBorn::new(&elements), charges, at)
}

/// A cluster built to reach every branch of the pair integral and both rescalings' slopes: a
/// large atom (ρ 3 Å, S 1) with two small ones inside its scaled sphere, one small atom wholly
/// inside another's, and overlapping pairs of both kinds. Positions in Å; charges deterministic.
fn cluster() -> (GeneralizedBorn, Vec<f64>, Vec<[f64; 3]>) {
    let atoms: [(f64, f64, [f64; 3]); 7] = [
        (3.0, 1.0, [0.0, 0.0, 0.0]),
        (1.2, 0.85, [0.4, 0.3, -0.2]),
        (1.5, 0.85, [-0.6, 0.5, 0.7]),
        (1.7, 0.72, [3.1, 0.4, 0.2]),
        (1.2, 0.85, [3.4, 1.0, 0.9]),
        (1.55, 0.79, [5.6, -0.8, 0.3]),
        (1.8, 0.96, [1.6, 3.9, -1.1]),
    ];
    let gb = GeneralizedBorn::from_radii(
        atoms.iter().map(|a| a.0 * ANGSTROM).collect(),
        atoms.iter().map(|a| a.1).collect(),
    );
    let charges = (0..atoms.len())
        .map(|i| 0.6 * (1.1 * i as f64 + 0.3).sin())
        .collect();
    let at = atoms.iter().map(|a| a.2.map(|x| x * ANGSTROM)).collect();
    (gb, charges, at)
}

/// How many ordered pairs of `gb` at `at` are in each branch of the pair integral: apart,
/// overlapping ρ̃, inside the other's sphere, the other inside its own.
fn branches(gb: &GeneralizedBorn, at: &[[f64; 3]]) -> [usize; 4] {
    let mut count = [0; 4];
    for i in 0..at.len() {
        let tilde = gb.radii()[i] - RADIUS_OFFSET;
        for j in 0..at.len() {
            if i == j {
                continue;
            }
            let s = gb.scales()[j] * (gb.radii()[j] - RADIUS_OFFSET);
            let d = [0, 1, 2].map(|k| at[i][k] - at[j][k]);
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            let b = if r + s <= tilde {
                3
            } else if tilde < s - r {
                2
            } else if r - s >= tilde {
                0
            } else {
                1
            };
            count[b] += 1;
        }
    }
    count
}

/// **The forces are minus the gradient.** Central differences of the energy, coordinate by
/// coordinate, at h = 10⁻⁴ Å. The tolerance on each, kcal mol⁻¹ Å⁻¹:
///
/// - truncation `h²/6 |E‴|`, `E‴` along that coordinate measured by the five-point stencil
///   `(E(2H) − 2E(H) + 2E(−H) − E(−2H)) / 2H³` at H = 10⁻³ Å and **taken twice**, plus the
///   stencil's rounding `3 δE / H³`;
/// - rounding `δE / h`, with `δE = 100 n ε M`: M bounds the sum of the magnitudes of eq 2's
///   terms, `½ K (Σ|q|)² / R_min` since `f ≥ √(R_i R_j)`, and `100 n ε` allows each of the n² terms
///   a few dozen roundings and the radii's sums of n terms their own.
///
/// And **the check has to be able to see something**: on every system most coordinates' force is
/// many tolerances large, and the system reaches every branch of the pair integral.
fn forces_are_the_gradient(name: &str, gb: &GeneralizedBorn, q: &[f64], at: &[[f64; 3]]) {
    let n = at.len();
    let mut forces = vec![[0.0; 3]; n];
    let e0 = gb.accumulate(q, at, &mut forces);
    assert!(e0.is_finite());
    let energy = |i: usize, k: usize, by: f64| {
        let mut moved = at.to_vec();
        moved[i][k] += by * ANGSTROM;
        kcal(gb.energy(q, &moved))
    };
    let r_min = gb
        .born_radii(at)
        .iter()
        .fold(f64::INFINITY, |m, r| m.min(*r))
        / ANGSTROM;
    let sum_q: f64 = q.iter().map(|x| x.abs()).sum();
    let delta_e = 100.0 * n as f64 * EPS * 0.5 * K * sum_q * sum_q / r_min;
    let (h, big_h) = (1e-4, 1e-3);
    let (mut worst, mut seen, mut total) = (0.0f64, 0, 0);
    for i in 0..n {
        for k in 0..3 {
            let fd = -(energy(i, k, h) - energy(i, k, -h)) / (2.0 * h);
            let third = (energy(i, k, 2.0 * big_h) - 2.0 * energy(i, k, big_h)
                + 2.0 * energy(i, k, -big_h)
                - energy(i, k, -2.0 * big_h))
                / (2.0 * big_h.powi(3));
            let m3 = 2.0 * third.abs() + 3.0 * delta_e / big_h.powi(3);
            let tolerance = h * h / 6.0 * m3 + delta_e / h;
            let analytic = forces[i][k] / (KCAL_PER_MOL / ANGSTROM);
            let miss = (analytic - fd).abs();
            assert!(
                miss <= tolerance,
                "{name}: atom {i} axis {k}: analytic {analytic} against {fd}, off {miss:e} > {tolerance:e}"
            );
            worst = worst.max(miss / tolerance);
            total += 1;
            if analytic.abs() > 100.0 * tolerance {
                seen += 1;
            }
        }
    }
    println!(
        "{name}: ΔG {:.6} kcal/mol; worst miss {worst:.3} of its tolerance; {seen} of {total} \
         coordinates have a force over 100 tolerances (δE = {delta_e:.2e} kcal/mol)",
        kcal(e0)
    );
    assert!(
        seen * 10 >= total * 8,
        "{name}: too few coordinates the check can see"
    );
}

#[test]
fn the_force_is_minus_the_gradient_on_aspirin() {
    let (_, gb, q, at) = aspirin();
    println!("aspirin's branches: {:?}", branches(&gb, &at));
    forces_are_the_gradient("aspirin, OBC II", &gb, &q, &at);
    let salty = gb.clone().with_kappa(solvation::debye_kappa(0.15));
    forces_are_the_gradient("aspirin, OBC II, 0.15 M salt", &salty, &q, &at);
    let hct = gb.with_rescaling(Rescaling::Hct);
    forces_are_the_gradient("aspirin, HCT", &hct, &q, &at);
}

#[test]
fn the_force_is_minus_the_gradient_on_every_branch() {
    let (gb, q, at) = cluster();
    let count = branches(&gb, &at);
    println!("cluster's branches (apart, overlapping ρ̃, inside, contains): {count:?}");
    assert!(
        count.iter().all(|&c| c > 0),
        "every branch reached: {count:?}"
    );
    forces_are_the_gradient("cluster, OBC II", &gb, &q, &at);
    forces_are_the_gradient(
        "cluster, OBC I",
        &gb.clone().with_rescaling(Rescaling::OBC_I),
        &q,
        &at,
    );
    forces_are_the_gradient("cluster, HCT", &gb.with_rescaling(Rescaling::Hct), &q, &at);
}

/// The energy reads only separations: a rigid motion moves it by rounding alone, and the forces
/// turn with the molecule.
#[test]
fn a_rigid_motion_changes_nothing() {
    let (_, gb, q, at) = aspirin();
    let mut f0 = vec![[0.0; 3]; at.len()];
    let e0 = gb.accumulate(&q, &at, &mut f0);
    let (c, s) = (0.6f64.cos(), 0.6f64.sin());
    let turned: Vec<[f64; 3]> = at
        .iter()
        .map(|p| {
            [
                c * p[0] - s * p[1] + 3.0 * ANGSTROM,
                s * p[0] + c * p[1] - 7.0 * ANGSTROM,
                p[2] + 11.0 * ANGSTROM,
            ]
        })
        .collect();
    let mut f1 = vec![[0.0; 3]; at.len()];
    let e1 = gb.accumulate(&q, &turned, &mut f1);
    println!(
        "aspirin: {:.12} and {:.12} kcal/mol moved",
        kcal(e0),
        kcal(e1)
    );
    assert!((e0 - e1).abs() <= 1e-12 * e0.abs());
    let scale = f0.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
    for (a, b) in f0.iter().zip(&f1) {
        let rotated = [c * a[0] - s * a[1], s * a[0] + c * a[1], a[2]];
        for k in 0..3 {
            assert!((rotated[k] - b[k]).abs() <= 1e-10 * scale);
        }
    }
    let net = f0
        .iter()
        .fold([0.0; 3], |m, f| [m[0] + f[0], m[1] + f[1], m[2] + f[2]]);
    for v in net {
        assert!(v.abs() <= 1e-12 * scale * at.len() as f64, "net force {v}");
    }
}

/// Opt-in: a force field is in vacuum unless given a solvent, and given one adds exactly the GB
/// model's energy and force to everything else.
#[test]
fn a_force_field_is_in_vacuum_until_given_water() {
    let (c, gb, q, at) = aspirin();
    let vacuum = ForceField::new(&c, &uff::assign(&c))
        .expect("aspirin")
        .with_charges(q.clone());
    assert!(vacuum.solvation().is_none());
    let e = vacuum.evaluate(&at);
    assert_eq!(e.energy.solvation, 0.0);
    let wet = vacuum.clone().with_generalized_born();
    assert_eq!(wet.solvation(), Some(&gb));
    let w = wet.evaluate(&at);
    let mut f = vec![[0.0; 3]; at.len()];
    let g = gb.accumulate(&q, &at, &mut f);
    assert_eq!(w.energy.solvation, g);
    assert!(
        (w.energy.total - e.energy.total - g).abs() <= 1e-13 * e.energy.total.abs().max(g.abs())
    );
    for i in 0..at.len() {
        for k in 0..3 {
            let d = w.forces[i][k] - e.forces[i][k] - f[i][k];
            assert!(d.abs() <= 1e-12 * f[i][k].abs().max(e.forces[i][k].abs()).max(1e-12));
        }
    }
    // Zero charges, zero solvation.
    let neutral = ForceField::new(&c, &uff::assign(&c))
        .expect("aspirin")
        .with_generalized_born();
    assert_eq!(neutral.evaluate(&at).energy.solvation, 0.0);
}

/// The area of sphere 1 (radius a) outside sphere 2 (radius b), centres d apart and overlapping:
/// `4πa² − 2πa h`, the cap's height `h = a − (d² + a² − b²) / (2d)`.
fn outside_cap(a: f64, b: f64, d: f64) -> f64 {
    let h = a - (d * d + a * a - b * b) / (2.0 * d);
    4.0 * std::f64::consts::PI * a * a - 2.0 * std::f64::consts::PI * a * h
}

/// **The surface converges to the exact area of two overlapping spheres**, and is exact for one.
/// Measured over eight pairs — a C–O, C–C, O–H and S–C bond's length and two contacts, each along
/// a different direction, so the lab-fixed point set meets each cap differently — grown by the
/// 1.4 Å probe, at point counts from 100 to 25 600. The worst error over the pairs at each count is
/// printed — measured 1.86%, 0.49%, 0.29%, 0.096% and 0.025%, no clean rate — and asserted is
/// that it falls at every step from 400 points and is below 0.5% at 1600, which the hydration table
/// uses: 0.005 kcal/mol on a nonpolar term of about 1 kcal/mol.
#[test]
fn the_surface_converges_to_two_spheres_exact_area() {
    let alone = surface_area(&[1.7 * ANGSTROM], &[[0.0; 3]], PROBE_RADIUS, 100)[0];
    let exact = 4.0 * std::f64::consts::PI * (3.1 * ANGSTROM).powi(2);
    assert!((alone - exact).abs() <= 4.0 * EPS * exact);
    // (radius a, radius b, distance, direction), Å.
    let pairs: [(f64, f64, f64, [f64; 3]); 8] = [
        (1.7, 1.5, 1.43, [0.0, -0.6, 0.8]),
        (1.7, 1.7, 1.53, [1.0, 0.0, 0.0]),
        (1.5, 1.2, 0.96, [0.36, 0.48, 0.8]),
        (1.8, 1.7, 1.82, [0.0, 0.0, 1.0]),
        (1.7, 1.2, 1.09, [-0.48, 0.6, 0.64]),
        (1.55, 1.2, 1.01, [0.6, 0.0, -0.8]),
        (1.7, 1.5, 3.4, [0.8, 0.6, 0.0]),
        (1.2, 1.2, 4.5, [0.0, 0.8, 0.6]),
    ];
    let mut previous = f64::INFINITY;
    for points in [100, 400, 1600, 6400, 25_600] {
        let mut worst = 0.0f64;
        for (a, b, d, u) in pairs {
            let want = outside_cap(a + 1.4, b + 1.4, d) + outside_cap(b + 1.4, a + 1.4, d);
            let at = [[0.0; 3], u.map(|x| x * d * ANGSTROM)];
            let got: f64 = surface_area(&[a * ANGSTROM, b * ANGSTROM], &at, PROBE_RADIUS, points)
                .iter()
                .sum::<f64>()
                / (ANGSTROM * ANGSTROM);
            worst = worst.max((got - want).abs() / want);
        }
        println!(
            "{points:6} points: worst {:.4}% off over the eight pairs",
            100.0 * worst
        );
        if points >= 400 {
            assert!(
                worst < previous,
                "{points} points: {worst} after {previous}"
            );
        }
        if points == 1600 {
            assert!(worst < 5e-3, "{points} points: {worst}");
        }
        previous = worst;
    }
}

/// **The radii and scale factors are the tables they are said to be**, typed here from the files
/// they come from, each opened: OpenMM's `customgbforces.py` as the research saved it (line
/// numbers that file's), and for Br and I, Bondi (1964) as tabulated secondarily.
#[test]
fn the_radii_and_scale_factors_are_the_tables() {
    // (element, ρ Å, its source, S, its source)
    let table: [(Element, f64, &str, f64, &str); 10] = [
        (
            Element::H,
            1.2,
            "_bondi_radii, line 243",
            0.85,
            "_SCREEN_PARAMETERS, line 396",
        ),
        (
            Element::C,
            1.7,
            "_bondi_radii, line 242",
            0.72,
            "_SCREEN_PARAMETERS, line 397",
        ),
        (
            Element::N,
            1.55,
            "_bondi_radii, line 245",
            0.79,
            "_SCREEN_PARAMETERS, line 398",
        ),
        (
            Element::O,
            1.5,
            "_bondi_radii, line 246",
            0.85,
            "_SCREEN_PARAMETERS, line 399",
        ),
        (
            Element::F,
            1.5,
            "_bondi_radii, line 247",
            0.88,
            "_SCREEN_PARAMETERS, line 400",
        ),
        (
            Element::P,
            1.85,
            "_bondi_radii, line 249",
            0.86,
            "_SCREEN_PARAMETERS, line 401",
        ),
        (
            Element::S,
            1.8,
            "_bondi_radii, line 250",
            0.96,
            "_SCREEN_PARAMETERS, line 402",
        ),
        (
            Element::Cl,
            1.7,
            "_bondi_radii, line 251",
            0.8,
            "_SCREEN_PARAMETERS default, line 403",
        ),
        (
            Element::Br,
            1.85,
            "Bondi 1964, secondary",
            0.8,
            "_SCREEN_PARAMETERS default, line 403",
        ),
        (
            Element::I,
            1.98,
            "Bondi 1964, secondary",
            0.8,
            "_SCREEN_PARAMETERS default, line 403",
        ),
    ];
    for (element, rho, from, s, s_from) in table {
        let got = solvation::intrinsic_radius(element) / ANGSTROM;
        assert!(
            (got - rho).abs() < 1e-12,
            "{element:?}: ρ {got} against {rho} ({from})"
        );
        assert_eq!(
            solvation::scale_factor(element),
            s,
            "{element:?}: S against {s} ({s_from})"
        );
    }
}

/// **Eq 2 with eq 3, typed here, for two bare ions at three separations, two dielectrics and two
/// salts.** Scale factors zero, so no descreening and R is exactly ρ̃: 1.41 and 1.61 Å. At 1, 2.3
/// and 4 Å Still's exponential is 0.896, 0.558 and 0.172, so its
/// exponent is seen; ε = 4 and 0.15 M salt reach the setter and the pair terms' screening.
#[test]
fn two_bare_ions_are_eq_2_with_eq_3_at_any_dielectric_and_salt() {
    let (q1, q2) = (0.7, -0.4);
    let (r1, r2): (f64, f64) = (1.41, 1.61);
    for eps in [80.0, 4.0] {
        for molar in [0.0, 0.15] {
            let gb =
                GeneralizedBorn::from_radii(vec![1.5 * ANGSTROM, 1.7 * ANGSTROM], vec![0.0, 0.0])
                    .with_solvent_dielectric(eps)
                    .with_kappa(solvation::debye_kappa(molar));
            assert_eq!(gb.solvent_dielectric(), eps);
            let kappa = 0.316 * molar.sqrt();
            for r in [1.0f64, 2.3, 4.0] {
                let f = (r * r + r1 * r2 * (-r * r / (4.0 * r1 * r2)).exp()).sqrt();
                let still = still_distance(r * ANGSTROM, r1 * ANGSTROM, r2 * ANGSTROM) / ANGSTROM;
                assert!(
                    (still - f).abs() <= 4.0 * EPS * f,
                    "eq 3 at {r} Å: {still} against {f}"
                );
                let term = |qq: f64, f: f64| -0.5 * K * qq * (1.0 - (-kappa * f).exp() / eps) / f;
                let want = term(q1 * q1, r1) + term(q2 * q2, r2) + 2.0 * term(q1 * q2, f);
                let got = kcal(gb.energy(&[q1, q2], &[[0.0; 3], [r * ANGSTROM, 0.0, 0.0]]));
                println!("ε {eps} salt {molar} M r {r} Å: {got:.12} against {want:.12}");
                assert!(
                    (got - want).abs() <= 1e-13 * want.abs(),
                    "ε {eps} salt {molar} r {r}: {got} against {want}"
                );
            }
        }
    }
}

/// p. 392's nonpolar term: 0.005 kcal mol⁻¹ Å⁻² over 100 Å² is 0.5 kcal/mol.
#[test]
fn the_nonpolar_term_is_five_thousandths_per_square_angstrom() {
    let a = solvation::nonpolar_energy(100.0 * ANGSTROM * ANGSTROM) / KCAL_PER_MOL;
    assert!((a - 0.5).abs() < 1e-14, "0.005 × 100 Å²: {a}");
}

/// The golden spiral's points are on the unit sphere, and at the midpoints of n equal-area bands,
/// `z = 1 − (2k + 1)/n`: so their z sum to zero, as an equal-area set's must (the sphere's
/// centroid), and each band `[z_k ± 1/n]` holds one. A set started at the pole sums to 1.
#[test]
fn the_point_set_is_n_equal_bands() {
    for n in [100, 1600, 1601] {
        let p = sphere_points(n);
        assert_eq!(p.len(), n);
        let sum: f64 = p.iter().map(|q| q[2]).sum();
        assert!(sum.abs() <= 4.0 * n as f64 * EPS, "n {n}: Σz = {sum}");
        for q in &p {
            let r = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
            assert!((r - 1.0).abs() <= 4.0 * EPS);
        }
        let mut z: Vec<f64> = p.iter().map(|q| q[2]).collect();
        z.sort_by(|a, b| b.partial_cmp(a).expect("finite"));
        for (k, zk) in z.iter().enumerate() {
            let top = 1.0 - 2.0 * k as f64 / n as f64;
            let bottom = 1.0 - 2.0 * (k + 1) as f64 / n as f64;
            assert!(
                *zk < top && *zk > bottom,
                "n {n}: point {k} at z {zk} outside its band"
            );
        }
    }
}
