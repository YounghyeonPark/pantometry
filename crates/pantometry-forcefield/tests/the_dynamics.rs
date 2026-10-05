//! **Molecular dynamics against closed forms**: velocity Verlet's shadow energy, BAOAB's exact
//! sampling of a harmonic well, the vibration of one bond, equipartition on aspirin, and the
//! determinism a run built from random numbers has to keep.
//!
//! The integrator is checked first on potentials whose dynamics is known exactly — an isotropic
//! harmonic well for one atom, written here as a [`Potential`], and a single C–H bond, which UFF
//! makes exactly harmonic in its length — and then on aspirin, where only averages and rates are
//! known: equipartition, and an energy error that falls as the square of the step.

mod common;

use common::{cross, dot, entry, positions, relaxed, sub, ANGSTROM};
use pantometry_core::conserved::{audit, quantity};
use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::{Bath, Component, Element, ForceField, MolecularDynamics, Potential};
use pantometry_units::BOLTZMANN;

const AIN: &str = include_str!("../components/AIN.cif");
const FS: f64 = 1e-15;

fn kb() -> f64 {
    BOLTZMANN.to_si()
}

/// Every atom tethered to its own origin by `½ k |x − x₀|²`: independent isotropic harmonic wells,
/// whose dynamics is linear and known in closed form.
struct Wells {
    k: f64,
    origins: Vec<[f64; 3]>,
}

impl Potential for Wells {
    fn energy_and_forces(&self, at: &[[f64; 3]], forces: &mut [[f64; 3]]) -> f64 {
        let mut e = 0.0;
        for ((x, o), f) in at.iter().zip(&self.origins).zip(forces.iter_mut()) {
            for c in 0..3 {
                let d = x[c] - o[c];
                e += 0.5 * self.k * d * d;
                f[c] = -self.k * d;
            }
        }
        e
    }
}

/// Mean, the standard error of the mean and the integrated autocorrelation time (in samples) of a
/// correlated series.
///
/// `σ² = 2 τ_int Var / N`, with `τ_int = ½ + Σ ρ(t)` summed over a self-consistent window — Sokal's
/// rule, stopping at the first `t ≥ 6 τ_int(t)` (A. D. Sokal, "Monte Carlo methods in statistical
/// mechanics", Cargèse lectures 1996, §3). Beyond that window the estimated `ρ(t)` is mostly noise,
/// and summing it would make the error bar itself noisy.
fn statistics(xs: &[f64]) -> (f64, f64, f64) {
    let n = xs.len();
    let mean = xs.iter().sum::<f64>() / n as f64;
    let var = xs.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / n as f64;
    let mut tau = 0.5;
    for t in 1..n / 2 {
        let c: f64 = (0..n - t)
            .map(|i| (xs[i] - mean) * (xs[i + t] - mean))
            .sum::<f64>()
            / ((n - t) as f64 * var);
        tau += c;
        if t as f64 >= 6.0 * tau {
            break;
        }
    }
    (mean, (2.0 * tau * var / n as f64).sqrt(), tau)
}

fn aspirin() -> (Component, ForceField) {
    let c = Component::from_ccd(AIN).expect("AIN parses");
    let ff = ForceField::new(&c, &uff::assign(&c)).expect("aspirin is supported");
    (c, ff)
}

/// Aspirin at the force field's minimum nearest the dictionary's coordinates.
fn relaxed_aspirin() -> (Component, ForceField, Vec<[f64; 3]>) {
    let (c, ff) = aspirin();
    let (_, at) = relaxed(&ff, &positions(&c));
    (c, ff, at)
}

fn elements(c: &Component) -> Vec<Element> {
    c.atoms().iter().map(|a| a.element).collect()
}

// ---------------------------------------------------------------------------------------------
// Masses
// ---------------------------------------------------------------------------------------------

/// **Every element's weight is CIAAW's**, typed here from the source rather than from the code:
/// "Abridged Standard Atomic Weights" (Commission on Isotopic Abundances and Atomic Weights,
/// ciaaw.org/abridged-atomic-weights.htm, the 2024 table, based on Prohaska et al., *Pure Appl.
/// Chem.* 94, 573 (2022)). Read 2026-10-05: H 1.0080(2), C 12.011(2), N 14.007(1), O 15.999(1),
/// F 18.998(1), P 30.974(1), S 32.06(2), Cl 35.45(1), Br 79.904(3), I 126.90(1). Exact equality,
/// because both are the same printed decimal; the formula-weight test below cannot see F, P, Cl,
/// Br or I in most entries and none of them closely, and this can.
#[test]
fn every_atomic_weight_is_the_ciaaw_abridged_value() {
    let table = [
        (Element::H, 1.0080),
        (Element::C, 12.011),
        (Element::N, 14.007),
        (Element::O, 15.999),
        (Element::F, 18.998),
        (Element::P, 30.974),
        (Element::S, 32.06),
        (Element::Cl, 35.45),
        (Element::Br, 79.904),
        (Element::I, 126.90),
    ];
    assert_eq!(table.len(), Element::ALL.len());
    for (e, weight) in table {
        assert_eq!(e.atomic_weight(), weight, "{}", e.symbol());
    }
}

/// **The masses add up to each dictionary entry's own formula weight**, for every entry in
/// `components/`.
///
/// `_chem_comp.formula_weight` is the dictionary's, computed by the wwPDB from its own table of
/// atomic weights, which is not this crate's: aspirin's 180.157 is the 2005-era 1.00794, 12.0107
/// and 15.9994, where [`Element::atomic_weight`] has CIAAW's abridged 1.0080, 12.011 and 15.999.
/// So the bound is earned from the uncertainty CIAAW states beside each abridged value (the same
/// table), summed over the atoms, plus half a unit in the dictionary's last printed digit.
#[test]
fn the_masses_add_up_to_every_entrys_formula_weight() {
    let files: [(&str, &str); 49] = [
        ("16R", include_str!("../components/16R.cif")),
        ("2F2", include_str!("../components/2F2.cif")),
        ("2ME", include_str!("../components/2ME.cif")),
        ("61G", include_str!("../components/61G.cif")),
        ("6AC", include_str!("../components/6AC.cif")),
        ("A1JFW", include_str!("../components/A1JFW.cif")),
        ("ACE", include_str!("../components/ACE.cif")),
        ("ACM", include_str!("../components/ACM.cif")),
        ("ACN", include_str!("../components/ACN.cif")),
        ("AIN", include_str!("../components/AIN.cif")),
        ("ALA", include_str!("../components/ALA.cif")),
        ("ARG", include_str!("../components/ARG.cif")),
        ("ASN", include_str!("../components/ASN.cif")),
        ("ASP", include_str!("../components/ASP.cif")),
        ("BNZ", include_str!("../components/BNZ.cif")),
        ("BZF", include_str!("../components/BZF.cif")),
        ("CCN", include_str!("../components/CCN.cif")),
        ("CYS", include_str!("../components/CYS.cif")),
        ("DEN", include_str!("../components/DEN.cif")),
        ("DMF", include_str!("../components/DMF.cif")),
        ("DMN", include_str!("../components/DMN.cif")),
        ("GLN", include_str!("../components/GLN.cif")),
        ("GLU", include_str!("../components/GLU.cif")),
        ("GLY", include_str!("../components/GLY.cif")),
        ("HIS", include_str!("../components/HIS.cif")),
        ("I4B", include_str!("../components/I4B.cif")),
        ("ILE", include_str!("../components/ILE.cif")),
        ("IND", include_str!("../components/IND.cif")),
        ("KEN", include_str!("../components/KEN.cif")),
        ("LEU", include_str!("../components/LEU.cif")),
        ("LYS", include_str!("../components/LYS.cif")),
        ("MEE", include_str!("../components/MEE.cif")),
        ("MET", include_str!("../components/MET.cif")),
        ("MOH", include_str!("../components/MOH.cif")),
        ("N4B", include_str!("../components/N4B.cif")),
        ("NME", include_str!("../components/NME.cif")),
        ("OXE", include_str!("../components/OXE.cif")),
        ("PEO", include_str!("../components/PEO.cif")),
        ("PHE", include_str!("../components/PHE.cif")),
        ("PRO", include_str!("../components/PRO.cif")),
        ("PXY", include_str!("../components/PXY.cif")),
        ("PYJ", include_str!("../components/PYJ.cif")),
        ("S2H", include_str!("../components/S2H.cif")),
        ("SER", include_str!("../components/SER.cif")),
        ("THR", include_str!("../components/THR.cif")),
        ("TME", include_str!("../components/TME.cif")),
        ("TRP", include_str!("../components/TRP.cif")),
        ("TYR", include_str!("../components/TYR.cif")),
        ("VAL", include_str!("../components/VAL.cif")),
    ];
    // CIAAW's stated uncertainty on each abridged value.
    let uncertainty = |e: Element| match e {
        Element::H => 0.0002,
        Element::C => 0.002,
        Element::N | Element::O | Element::F | Element::P => 0.001,
        Element::S => 0.02,
        Element::Cl | Element::I => 0.01,
        Element::Br => 0.003,
    };
    let mut worst = 0.0f64;
    for (code, text) in files {
        let c = Component::from_ccd(text).unwrap_or_else(|e| panic!("{code}: {e}"));
        let stated: f64 = text
            .lines()
            .find_map(|l| l.strip_prefix("_chem_comp.formula_weight"))
            .unwrap_or_else(|| panic!("{code} states no formula weight"))
            .trim()
            .parse()
            .unwrap_or_else(|e| panic!("{code}: {e}"));
        let sum: f64 = c.atoms().iter().map(|a| a.element.atomic_weight()).sum();
        let bound: f64 = c
            .atoms()
            .iter()
            .map(|a| uncertainty(a.element))
            .sum::<f64>()
            + 0.0005;
        assert!(
            (sum - stated).abs() <= bound,
            "{code}: the masses sum to {sum:.4}, the entry states {stated}, bound {bound:.4}"
        );
        // And the mass is the weight over N_A, to the rounding of one multiplication and one
        // division.
        for a in c.atoms() {
            let m = a.element.mass() * uff::AVOGADRO / 1e-3;
            assert!((m - a.element.atomic_weight()).abs() <= 4.0 * f64::EPSILON * m);
        }
        worst = worst.max((sum - stated).abs() / bound);
    }
    println!("worst |Σ weights − formula weight| / bound: {worst:.3}");
}

// ---------------------------------------------------------------------------------------------
// Velocity Verlet against closed forms
// ---------------------------------------------------------------------------------------------

/// **Velocity Verlet's energy error on a harmonic well is `E₀ (ωh)²/4`, and halving the step
/// quarters it.**
///
/// For `x″ = −ω²x` the scheme conserves `H̃ = ½mv² + ½kx²(1 − (ωh)²/4)` exactly, so
/// `E − H̃ = ½kx²(ωh)²/4 ≥ 0`. Released from rest at amplitude `A`, `E` starts at `E₀ = H̃/(1 −
/// (ωh)²/4)` — its largest value — and is lowest, at `H̃`, where the atom crosses the origin: the
/// largest departure from the start is `E₀(ωh)²/4`, a closed form rather than a rate. Over a
/// thousand periods the sampled crossing comes within a part in 10⁶ of the origin, so the measured
/// largest departure meets it to that, and can exceed it only by rounding.
#[test]
fn velocity_verlet_misses_the_energy_by_exactly_its_shadow() {
    let (m, k) = (Element::H.mass(), 500.0);
    let omega = (k / m).sqrt();
    let amplitude = 0.1 * ANGSTROM;
    let wells = Wells {
        k,
        origins: vec![[0.0; 3]],
    };
    let worst = |wh: f64| {
        let h = wh / omega;
        let mut at = vec![[amplitude, 0.0, 0.0]];
        let mut md = MolecularDynamics::new(vec![m]);
        md.prepare(&wells, &at);
        let e0 = md.kinetic_energy() + md.potential_energy().unwrap();
        let steps = (1000.0 * std::f64::consts::TAU / wh) as usize;
        let mut worst = 0.0f64;
        for _ in 0..steps {
            md.step(&wells, &mut at, h);
            let e = md.kinetic_energy() + md.potential_energy().unwrap();
            worst = worst.max((e0 - e) / e0);
        }
        worst
    };
    let mut previous = f64::NAN;
    for wh in [0.4, 0.2, 0.1] {
        let measured = worst(wh);
        let closed = wh * wh / 4.0;
        println!("ωh = {wh}: largest relative departure {measured:.9e}, closed form {closed:.9e}");
        assert!(
            measured <= closed * (1.0 + 1e-9) && measured >= closed * (1.0 - 1e-6),
            "ωh = {wh}: {measured:e} against {closed:e}"
        );
        if previous.is_finite() {
            let order = (previous / measured).log2();
            assert!((order - 2.0).abs() < 1e-5, "order {order}");
        }
        previous = measured;
    }
}

/// **One C–H bond vibrates at `√(k/μ)`, as velocity Verlet discretises it**, with `μ` from the two
/// elements' masses and `k` UFF's.
///
/// Placed on one axis and released from rest, the two atoms move only along it, and UFF's stretch
/// is exactly harmonic in the length: so the bond length obeys `r″ = −(k/μ)(r − r₀)` and velocity
/// Verlet's solution is `r − r₀ = A cos(ω̃ n h)` with `cos ω̃h = 1 − (ωh)²/2` — exactly, step by step,
/// not to an order in `h`. A mass taken from the wrong element moves `ω` by its square root.
#[test]
fn a_bond_vibrates_at_the_frequency_its_masses_and_stiffness_give() {
    let c = entry(
        &[("C1", "C", [0.0, 0.0, 0.0]), ("H1", "H", [1.2, 0.0, 0.0])],
        &[("C1", "H1", "SING", false)],
    );
    let ff = ForceField::new(&c, &uff::assign(&c)).expect("supported");
    assert_eq!(ff.stretches().len(), 1);
    assert!(ff.bends().is_empty() && ff.pairs().is_empty() && ff.torsions().is_empty());
    let s = ff.stretches()[0];
    let (mc, mh) = (Element::C.mass(), Element::H.mass());
    let mu = mc * mh / (mc + mh);
    let omega = (s.force_constant / mu).sqrt();
    let amplitude = 0.05 * ANGSTROM;
    let mut at = vec![[0.0; 3], [s.natural_length + amplitude, 0.0, 0.0]];
    let mut md = MolecularDynamics::for_elements(elements(&c));
    let h = 0.25 * FS;
    let discrete = (1.0 - (omega * h) * (omega * h) / 2.0).acos() / h;
    println!(
        "C–H: k = {:.1} kcal/mol/Å², μ = {:.4} u, period {:.3} fs",
        s.force_constant * ANGSTROM * ANGSTROM / KCAL_PER_MOL,
        mu / Element::H.mass() * Element::H.atomic_weight(),
        std::f64::consts::TAU / omega / FS
    );
    let mut worst = 0.0f64;
    for n in 1..=2000 {
        md.step(&ff, &mut at, h);
        let r = at[1][0] - at[0][0];
        let expected = amplitude * (discrete * n as f64 * h).cos();
        worst = worst.max((r - s.natural_length - expected).abs() / amplitude);
        assert_eq!([at[0][1], at[0][2], at[1][1], at[1][2]], [0.0; 4]);
    }
    // The bound is the rounding, summed linearly over the 2000 steps as a worst case:
    // - each step rounds the H position, 1.2e-10 m, to half an ulp, 1.3e-26 m, which is 2.6e-15
    //   of the 5e-12 m amplitude: 5.2e-12 over 2000 steps;
    // - the closed form's angle ω̃h = acos(1 − (ωh)²/2), with ωh = 0.136: rounding 1 − x by
    //   1.1e-16 moves acos by 1.1e-16 / √(2x) = 8.2e-16 rad, which over 2000 steps is a phase
    //   error of 1.6e-12 rad, the same of the amplitude.
    // That is 6.8e-12, so 1e-11; the measured worst is 4.8e-13. A drift scaled by 1 + 1e-6 moves ω̃
    // by 5e-7 of itself and the phase after 273 rad by 1.4e-4, seven orders above it.
    assert!(worst < 1e-11, "worst relative departure {worst:e}");
    println!("C–H against r₀ + A cos(ω̃nh): worst {worst:.3e} of the amplitude");
}

// ---------------------------------------------------------------------------------------------
// BAOAB against closed forms
// ---------------------------------------------------------------------------------------------

/// **BAOAB samples a harmonic well's positions exactly at any stable step**, and its momenta at
/// the closed forms its splitting gives.
///
/// For one atom in `½k|x|²` under BAOAB, the stationary distribution is Gaussian with, at a whole
/// step, `⟨x²⟩ = k_BT/k` per component **independent of h**, and `⟨m v²⟩ = k_BT (1 − (ωh)²/4)`;
/// right after `O`, `⟨m v²⟩ = k_BT` exactly. (The linear map's stationary covariance, solved in
/// closed form; Leimkuhler and Matthews note the harmonic oscillator's "special cancellations",
/// which is what makes it exact here.) Checked at ωh = 0.5, 1.0 and 1.9 — the last close to the
/// stability edge at 2, where the whole-step kinetic temperature is a tenth of the bath's and the
/// positions are still exact.
///
/// The bound on each mean is four standard errors, the standard error from the series' own
/// integrated autocorrelation time. ABOBA — the same pieces in another order — samples positions
/// exactly too, and its whole-step `⟨m v²⟩` is `k_BT/(1 − (ωh)²/4)`, 33% high at ωh = 1: the
/// momentum checks are what tell the orders apart.
#[test]
fn baoab_samples_a_harmonic_well_exactly_at_any_stable_step() {
    let (m, k, temperature) = (Element::C.mass(), 300.0, 300.0);
    let omega = (k / m).sqrt();
    let kt = kb() * temperature;
    let wells = Wells {
        k,
        origins: vec![[0.0; 3]],
    };
    for (wh, seed) in [(0.5, 11u64), (1.0, 12), (1.9, 13)] {
        let h = wh / omega;
        let mut at = vec![[0.0; 3]];
        let mut md = MolecularDynamics::new(vec![m])
            .with_bath(Bath::Langevin {
                temperature,
                friction: omega,
                seed,
            })
            .thermalised(&at, temperature, seed);
        md.run(&wells, &mut at, h, 2000);
        let steps = 150_000;
        let (mut x2, mut v2, mut half) = (Vec::new(), Vec::new(), Vec::new());
        for _ in 0..steps {
            md.step(&wells, &mut at, h);
            x2.push(k * dot(at[0], at[0]) / (3.0 * kt));
            let v = md.velocities()[0];
            v2.push(m * dot(v, v) / (3.0 * kt));
            half.push(2.0 * md.half_step_kinetic_energy() / (3.0 * kt));
        }
        for (what, series, exact) in [
            ("k⟨x²⟩/k_BT", &x2, 1.0),
            ("⟨mv²⟩/k_BT at the step", &v2, 1.0 - wh * wh / 4.0),
            ("⟨mv²⟩/k_BT after O", &half, 1.0),
        ] {
            let (mean, error, tau) = statistics(series);
            let z = (mean - exact) / error;
            println!(
                "ωh = {wh}: {what} = {mean:.5} ± {error:.5} (τ = {tau:.2} steps), exact \
                 {exact:.5}, z = {z:+.2}"
            );
            assert!(
                z.abs() < 4.0,
                "ωh = {wh}: {what} {mean} ± {error}, exact {exact}"
            );
        }
    }
}

/// **What a BAOAB step leaves in the books, on a harmonic well, in closed form.** For one atom in
/// `½k|x|²`, write `P = √m v`, `Q = √k x` and `a = ωh/2`, so the energy is `½(P² + Q²)`. Worked by
/// hand from the five lines:
///
/// - `BA` takes `(Q, P)` to `(Q₁, P₁)` with `P₁ = P − aQ`, `Q₁ = Q + aP₁`, and changes the energy
///   by exactly `(a²/2)(P₁² − Q²)`;
/// - `O` changes `P₁` to `P₂` at fixed `Q₁`, which is the bath's work `½(P₂² − P₁²)`, and is not in
///   the books' change;
/// - `AB` takes `(Q₁, P₂)` to `(Q₂, P₃)` and changes the energy by exactly `(a²/2)(Q₂² − P₂²)`.
///
/// So one step moves kinetic + potential − work by exactly
///
/// ```text
/// Δ = (a²/2) (P₁² − P₂² + Q₂² − Q²),
/// ```
///
/// whatever the noise was. Its leading part is `(ωh)²(P₁² − P₂²)/8`, the `h²` times the bath's
/// exchange that makes the books a random walk; in velocity Verlet `P₂ = P₁` and it is the shadow
/// Hamiltonian's band. Every quantity in it is readable from outside: `Q` and `Q₂` are the positions
/// before and after, `P₁² = m|v + (h/2)F/m|²` from the velocity before the step, and `P₂² = 2 KE`
/// after `O` ([`MolecularDynamics::half_step_kinetic_energies`]). So the books' change is checked
/// step by step, to rounding, rather than as a rate.
///
/// The bound is rounding: each side is a sum of a few terms of a few `k_BT` (the ledger's
/// totals also carry the bath's accumulated work, which over these 1000 steps is tens of `k_BT`),
/// each rounded to 1.1e-16 of itself, so 1e-13 of the larger of `k_BT` and the work is far from
/// both a correct step and any wrong one: a bath whose work is half counted is off by half the
/// exchange, about `√(γh) k_BT` = 0.7 `k_BT` here, and booking the kinetic energy after `O` is off by
/// a half kick, `O(ωh) k_BT`.
#[test]
fn a_step_leaves_exactly_its_closed_form_in_the_books() {
    let (m, k, temperature) = (Element::C.mass(), 300.0, 300.0);
    let omega = (k / m).sqrt();
    let kt = kb() * temperature;
    let wh = 0.5;
    let h = wh / omega;
    let a = wh / 2.0;
    let wells = Wells {
        k,
        origins: vec![[0.0; 3]],
    };
    let mut at = vec![[0.0; 3]];
    let mut md = MolecularDynamics::new(vec![m])
        .with_bath(Bath::Langevin {
            temperature,
            friction: omega,
            seed: 31,
        })
        .thermalised(&at, temperature, 31);
    md.run(&wells, &mut at, h, 200);
    let books = |md: &MolecularDynamics| md.ledger().get(quantity::ENERGY).unwrap();
    let mut worst = 0.0f64;
    let mut largest = 0.0f64;
    for _ in 0..1000 {
        let (x, v, before) = (at[0], md.velocities()[0], books(&md));
        md.step(&wells, &mut at, h);
        let p1 = [0, 1, 2].map(|c| v[c] - (h / 2.0) * (k / m) * x[c]);
        let p1_squared = m * dot(p1, p1);
        let p2_squared = 2.0 * md.half_step_kinetic_energies()[0];
        let q_squared = k * dot(x, x);
        let q2_squared = k * dot(at[0], at[0]);
        let closed = 0.5 * a * a * (p1_squared - p2_squared + q2_squared - q_squared);
        let measured = books(&md) - before;
        let scale = kt.max(md.thermostat_work().abs());
        worst = worst.max((measured - closed).abs() / scale);
        largest = largest.max(closed.abs() / kt);
    }
    println!(
        "books' change per step against its closed form: worst {worst:.2e} of the scale; the \
         change itself up to {largest:.3} k_BT"
    );
    assert!(
        largest > 1e-2,
        "the check has to be of something: {largest}"
    );
    assert!(worst < 1e-13, "{worst:e}");
}

/// **A bath at absolute zero is pure friction, and one BAOAB step is the five lines of the
/// paper's Appendix**, worked by hand. One atom released from rest at `A` in a well of stiffness
/// `k`, with `γh = ln 2` so that `c₁ = ½` and, at 0 K, `c₃ = 0`:
///
/// ```text
/// B  v₁ = −(h/2) kA/m         A  x₁ = A + (h/2) v₁
/// O  v₂ = v₁ / 2              A  x₂ = x₁ + (h/2) v₂
/// B  v₃ = v₂ − (h/2) k x₂/m
/// ```
///
/// The kinetic energy after `O` is `½ m v₂²`, the bath's work `½ m (v₂² − v₁²)`, and the position
/// and velocity at the step `x₂` and `v₃` — each to a few roundings. A kinetic energy read before
/// `O` instead of after it is the same in distribution for a harmonic well at any temperature, so
/// no statistical test can tell; this one can, by a factor of four.
#[test]
fn a_bath_at_absolute_zero_is_the_appendix_by_hand() {
    let (m, k, amplitude) = (Element::O.mass(), 400.0, 0.2 * ANGSTROM);
    let h = 0.7 * FS;
    let wells = Wells {
        k,
        origins: vec![[0.0; 3]],
    };
    let mut at = vec![[amplitude, 0.0, 0.0]];
    let mut md = MolecularDynamics::new(vec![m]).with_bath(Bath::Langevin {
        temperature: 0.0,
        friction: std::f64::consts::LN_2 / h,
        seed: 1,
    });
    md.step(&wells, &mut at, h);
    let v1 = -(h / 2.0) * k * amplitude / m;
    let x1 = amplitude + (h / 2.0) * v1;
    let v2 = v1 / 2.0;
    let x2 = x1 + (h / 2.0) * v2;
    let v3 = v2 - (h / 2.0) * k * x2 / m;
    let close = |a: f64, b: f64| (a - b).abs() <= 1e-14 * b.abs();
    assert!(close(md.half_step_kinetic_energy(), 0.5 * m * v2 * v2));
    assert!(close(md.thermostat_work(), 0.5 * m * (v2 * v2 - v1 * v1)));
    assert!(close(at[0][0], x2), "{} against {x2}", at[0][0]);
    assert!(close(md.velocities()[0][0], v3));
    assert_eq!([at[0][1], at[0][2]], [0.0, 0.0]);
}

/// **The bath's kicks are independent between atoms and between components.** Two atoms in
/// separate wells under one bath share nothing but the generator, so in the stationary state
/// `⟨x₀ · x₁⟩ = 0` and `⟨x₀ₓ x₀ᵧ⟩ = 0` exactly, where `⟨x²⟩ = k_BT/k`. Equipartition cannot see a
/// bath that kicks every atom the same way — each atom alone is still at the right temperature —
/// and this can: identical kicks make the two atoms' coordinates converge on each other, a
/// correlation of one. Four standard errors, from each product series' own autocorrelation.
#[test]
fn the_kicks_are_independent_between_atoms_and_components() {
    let (m, k, temperature) = (Element::C.mass(), 300.0, 300.0);
    let omega = (k / m).sqrt();
    let unit = kb() * temperature / k;
    let wells = Wells {
        k,
        origins: vec![[0.0; 3], [1e-9, 0.0, 0.0]],
    };
    let mut at = vec![[0.0; 3], [1e-9, 0.0, 0.0]];
    // Unequal masses, so a bath that used one mass for both would put them at different
    // temperatures: carbon and hydrogen, hydrogen at ωh = 1 and carbon at 0.29.
    let h = 1.0 / (k / Element::H.mass()).sqrt();
    let mut md = MolecularDynamics::new(vec![m, Element::H.mass()])
        .with_bath(Bath::Langevin {
            temperature,
            friction: omega,
            seed: 17,
        })
        .thermalised(&at, temperature, 17);
    md.run(&wells, &mut at, h, 2000);
    let (mut between, mut within) = (Vec::new(), Vec::new());
    let (mut x0, mut x1) = (Vec::new(), Vec::new());
    for _ in 0..100_000 {
        md.step(&wells, &mut at, h);
        let d1 = sub(at[1], wells.origins[1]);
        x0.push(dot(at[0], at[0]) / (3.0 * unit));
        x1.push(dot(d1, d1) / (3.0 * unit));
        between.push(dot(at[0], d1) / (3.0 * unit));
        within.push((at[0][0] * at[0][1] + at[1][1] * at[1][2] + d1[0] * at[0][2]) / (3.0 * unit));
    }
    for (what, series) in [("between atoms", &between), ("between components", &within)] {
        let (mean, error, tau) = statistics(series);
        println!("correlation {what}: {mean:+.5} ± {error:.5} (τ = {tau:.2} steps)");
        assert!((mean / error).abs() < 4.0, "{what}: {mean} ± {error}");
    }
    // Each atom samples k⟨x²⟩ = k_BT whatever its mass — BAOAB's positions are exact at any step.
    for (what, series) in [("carbon", &x0), ("hydrogen", &x1)] {
        let (mean, error, tau) = statistics(series);
        println!("{what}: k⟨x²⟩/k_BT = {mean:.5} ± {error:.5} (τ = {tau:.2} steps)");
        assert!(
            ((mean - 1.0) / error).abs() < 4.0,
            "{what}: {mean} ± {error}"
        );
    }
}

/// **Thermalised velocities are Maxwell–Boltzmann on the degrees of freedom left**: removing the
/// centre-of-mass and angular momentum is an orthogonal projection in the mass-weighted metric, so
/// `X = 2 KE / k_BT` is χ² with `k = 3N − 6` degrees of freedom — mean `k` and variance `2k`,
/// exactly. The draws are independent, so the bounds need no autocorrelation time:
///
/// - over 8000 seeds the mean of `X/k` has standard deviation `√(2 / (8000 k))` = 0.0021 for
///   aspirin, and the bound is four of those, 0.84%;
/// - the sample variance of `X`, over `2k`, has standard deviation `√((2 + 12/k) / 7999)` =
///   0.0166 — χ²'s excess kurtosis is `12/k` — and the bound is four of those, 6.6%.
///
/// A draw 5% cold is 24 standard deviations off the mean; a width at `2k_BT` is 100% off, and the
/// removed momenta left counted, `6/k` = 10.5%.
#[test]
fn thermalised_velocities_are_maxwell_boltzmann_on_what_is_left() {
    let (c, _, start) = relaxed_aspirin();
    let temperature = 300.0;
    let kt = kb() * temperature;
    let seeds = 8000u64;
    let k = (3 * c.atoms().len() - 6) as f64;
    let x: Vec<f64> = (0..seeds)
        .map(|seed| {
            let md = MolecularDynamics::for_elements(elements(&c)).thermalised(
                &start,
                temperature,
                seed,
            );
            2.0 * md.kinetic_energy() / kt
        })
        .collect();
    let n = seeds as f64;
    let mean = x.iter().sum::<f64>() / n;
    let variance = x.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (n - 1.0);
    let (sd_mean, sd_variance) = (
        (2.0 / (n * k)).sqrt(),
        ((2.0 + 12.0 / k) / (n - 1.0)).sqrt(),
    );
    let (m, v) = (mean / k, variance / (2.0 * k));
    println!(
        "⟨X⟩/k = {m:.5} (σ {sd_mean:.5}, z = {:+.2}); Var X/2k = {v:.4} (σ {sd_variance:.4}, \
         z = {:+.2})",
        (m - 1.0) / sd_mean,
        (v - 1.0) / sd_variance
    );
    assert!(((m - 1.0) / sd_mean).abs() < 4.0, "mean {m}");
    assert!(((v - 1.0) / sd_variance).abs() < 4.0, "variance {v}");
}

// ---------------------------------------------------------------------------------------------
// Aspirin
// ---------------------------------------------------------------------------------------------

/// **Equipartition on aspirin: `⟨KE⟩ = (3 N_free / 2) k_BT` under a Langevin bath**, read from the
/// kinetic energy right after `O` and counted over every free atom — the bath re-thermalises the
/// centre of mass and the rotations [`MolecularDynamics::thermalised`] removed, so all `3N` count.
/// With six ring carbons frozen, `3(N − 6)`.
///
/// **And per element**: every free atom, whatever its mass, has `⟨½ m v²⟩ = (3/2) k_BT`, so the
/// mean over aspirin's hydrogens, its carbons and its oxygens is each `(3/2) k_BT`. The total
/// alone cannot see a bath that kicks with the wrong mass: with the noise width built from the
/// mean free mass instead of each atom's own, hydrogens settle near 155 K and the heavy atoms
/// near 386 K, and the sum is exactly right.
///
/// A strong bath (20 ps⁻¹) so the kinetic energy decorrelates in tens of steps, the first
/// picosecond discarded, and four standard errors from the measured autocorrelation time, for
/// the total and for each element's mean.
#[test]
fn aspirin_shares_its_kinetic_energy_equally() {
    let (c, ff, start) = relaxed_aspirin();
    let temperature = 300.0;
    let kt = kb() * temperature;
    let n = c.atoms().len();
    let ring: Vec<usize> = ["C1", "C2", "C3", "C4", "C5", "C6"]
        .iter()
        .map(|name| common::index(&c, name))
        .collect();
    for frozen in [Vec::new(), ring] {
        let mut mask = vec![false; n];
        for &i in &frozen {
            mask[i] = true;
        }
        let mut at = start.clone();
        let mut md = MolecularDynamics::for_elements(elements(&c))
            .with_bath(Bath::Langevin {
                temperature,
                friction: 20e12,
                seed: 0xA591_0001,
            })
            .with_frozen(mask)
            .thermalised(&at, temperature, 0xA591_0002);
        let free = 3 * (n - frozen.len());
        assert_eq!(md.degrees_of_freedom(), free);
        let h = 0.5 * FS;
        md.run(&ff, &mut at, h, 2000);
        let mut series = Vec::new();
        let mut temperatures = Vec::new();
        let groups: Vec<(Element, Vec<usize>)> = [Element::H, Element::C, Element::O]
            .into_iter()
            .map(|e| {
                let members = (0..n)
                    .filter(|&i| c.atoms()[i].element == e && !md.frozen()[i])
                    .collect();
                (e, members)
            })
            .collect();
        let mut per_element: Vec<Vec<f64>> = vec![Vec::new(); groups.len()];
        for _ in 0..30_000 {
            md.step(&ff, &mut at, h);
            series.push(md.half_step_kinetic_energy() / (0.5 * free as f64 * kt));
            temperatures.push(md.half_step_temperature());
            let each = md.half_step_kinetic_energies();
            for ((_, members), out) in groups.iter().zip(per_element.iter_mut()) {
                let sum: f64 = members.iter().map(|&i| each[i]).sum();
                out.push(sum / (1.5 * kt * members.len() as f64));
            }
        }
        for ((e, members), out) in groups.iter().zip(&per_element) {
            assert!(!members.is_empty());
            let (mean, error, tau) = statistics(out);
            let z = (mean - 1.0) / error;
            println!(
                "{} frozen, {} ({}): ⟨½mv²⟩/(3/2 k_BT) = {mean:.4} ± {error:.4} (τ = {tau:.1}), \
                 z = {z:+.2}",
                frozen.len(),
                e.symbol(),
                members.len()
            );
            assert!(z.abs() < 4.0, "{}: {mean} ± {error}", e.symbol());
        }
        let (mean, error, tau) = statistics(&series);
        let z = (mean - 1.0) / error;
        println!(
            "{} frozen: ⟨KE⟩/(3N_free/2 k_BT) = {mean:.4} ± {error:.4} (τ = {tau:.1} steps), \
             z = {z:+.2}",
            frozen.len()
        );
        assert!(z.abs() < 4.0, "{mean} ± {error}");
        let t = temperatures.iter().sum::<f64>() / temperatures.len() as f64;
        assert!(
            (t / temperature - mean).abs() < 1e-9,
            "the temperature is the same average"
        );
    }
}

/// The step-size table: NVE aspirin at 300 K for a picosecond at each step, from the same start.
/// Printed, for the record in the CHANGELOG and the module documentation; the assertions on it are
/// in `the_energy_error_falls_as_the_power_of_the_step_it_should`.
#[test]
#[ignore = "a measurement that prints a table and asserts nothing: run with --release -- --ignored --nocapture"]
fn the_step_size_measured() {
    let (c, ff, start) = relaxed_aspirin();
    let base = MolecularDynamics::for_elements(elements(&c)).thermalised(&start, 300.0, 7);
    println!("| dt (fs) | RMS ΔE (kcal/mol) | max |ΔE| | drift (kcal/mol/ps) | ⟨T⟩ (K) |");
    for dt in [0.125, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0] {
        let h = dt * FS;
        let steps = (10_000.0 / dt) as usize / 10;
        let mut at = start.clone();
        let mut md = base.clone();
        md.prepare(&ff, &at);
        let e0 = md.kinetic_energy() + md.potential_energy().unwrap();
        let (mut sq, mut worst, mut tsum) = (0.0, 0.0f64, 0.0);
        let (mut st, mut se, mut stt, mut ste) = (0.0, 0.0, 0.0, 0.0);
        let mut blew = None;
        for s in 1..=steps {
            md.step(&ff, &mut at, h);
            let e = (md.kinetic_energy() + md.potential_energy().unwrap() - e0) / KCAL_PER_MOL;
            if !e.is_finite() || e.abs() > 1e3 {
                blew = Some(s);
                break;
            }
            let t = s as f64 * dt / 1000.0;
            sq += e * e;
            worst = worst.max(e.abs());
            tsum += md.temperature();
            st += t;
            se += e;
            stt += t * t;
            ste += t * e;
        }
        if let Some(s) = blew {
            println!(
                "| {dt} | blew up at step {s} ({:.0} fs) | | | |",
                s as f64 * dt
            );
            continue;
        }
        let nf = steps as f64;
        let slope = (nf * ste - st * se) / (nf * stt - st * st);
        println!(
            "| {dt} | {:.4} | {worst:.4} | {slope:+.4} | {:.1} |",
            (sq / nf).sqrt(),
            tsum / nf
        );
    }
}

/// RMS departure of aspirin's books from their start over `picoseconds` at step `dt` (fs), from
/// the relaxed geometry at 300 K, velocities from `seed`: in NVE the energy, under a 300 K bath
/// at 5 ps⁻¹ (kicks from `seed` too) kinetic + potential − thermostat work.
fn books_error(dt: f64, picoseconds: f64, seed: u64, bath: bool) -> [f64; 2] {
    let (c, ff, start) = relaxed_aspirin();
    let mut at = start.clone();
    let mut md = MolecularDynamics::for_elements(elements(&c));
    if bath {
        md = md.with_bath(Bath::Langevin {
            temperature: 300.0,
            friction: 5e12,
            seed,
        });
    }
    let mut md = md.thermalised(&start, 300.0, seed);
    md.prepare(&ff, &at);
    let e0 = md.ledger().get(quantity::ENERGY).unwrap();
    let steps = (picoseconds * 1000.0 / dt).round() as usize;
    let (mut sq, mut inc, mut last) = (0.0, 0.0, 0.0);
    for _ in 0..steps {
        md.step(&ff, &mut at, dt * FS);
        let e = md.ledger().get(quantity::ENERGY).unwrap() - e0;
        sq += e * e;
        inc += (e - last) * (e - last);
        last = e;
    }
    let n = steps as f64;
    [
        (sq / n).sqrt() / KCAL_PER_MOL,
        (inc / n).sqrt() / KCAL_PER_MOL,
    ]
}

/// The measured orders `log₂(error(1 fs)/error(0.5 fs))` and `log₂(error(0.5)/error(0.25))` over
/// 0.5 ps, from start `seed`: of the RMS departure from the start, then of the RMS change per
/// step.
fn orders(seed: u64, bath: bool) -> [f64; 4] {
    let e: Vec<[f64; 2]> = [1.0, 0.5, 0.25]
        .iter()
        .map(|&dt| books_error(dt, 0.5, seed, bath))
        .collect();
    let o = |k: usize, i: usize| (e[i][k] / e[i + 1][k]).log2();
    [o(0, 0), o(0, 1), o(1, 0), o(1, 1)]
}

/// How much each measured order varies from one start to another, over sixteen: the spread the
/// bounds in `the_energy_error_falls_as_the_power_of_the_step_it_should` are earned from.
#[test]
#[ignore = "a measurement that prints a table and asserts nothing: run with --release -- --ignored --nocapture"]
fn the_orders_measured_over_many_starts() {
    for bath in [false, true] {
        let runs: Vec<[f64; 4]> = (100..116).map(|seed| orders(seed, bath)).collect();
        for (k, what) in [(0, "RMS departure"), (2, "RMS change per step")] {
            let all: Vec<f64> = runs.iter().flat_map(|r| [r[k], r[k + 1]]).collect();
            let mean = all.iter().sum::<f64>() / all.len() as f64;
            let sd = (all.iter().map(|o| (o - mean) * (o - mean)).sum::<f64>()
                / (all.len() - 1) as f64)
                .sqrt();
            let (lo, hi) = all
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), o| {
                    (l.min(*o), h.max(*o))
                });
            println!(
                "{}, {what}: {} orders, mean {mean:.4}, sd {sd:.4}, range [{lo:.4}, {hi:.4}]",
                if bath { "Langevin books" } else { "NVE energy" },
                all.len()
            );
        }
    }
}

/// **The energy error falls as the power of the step it should — in NVE, and in the books under a
/// bath.** At 1, 0.5 and 0.25 fs over 0.5 ps, from one start.
///
/// - **NVE, the RMS departure of the energy from its start: `h²`.** Velocity Verlet conserves a
///   shadow Hamiltonian `H + O(h²)`, so the energy wanders in a band of that width.
/// - **NVE, the RMS change per step: `h³`.** The energy is the shadow's constant plus `h² g(q, p)`
///   for a smooth `g`, and one step moves `(q, p)` by `O(h)`.
/// - **Langevin, the RMS departure of the books — kinetic + potential − thermostat work — from
///   their start: `h²`.** Worked on a harmonic mode (`m = k = 1`): the half-step pair `BA` moves
///   the energy by `h²(p² − q²)/8 + O(h³)` and `AB` by minus that at the momentum it meets. In
///   velocity Verlet the two cancel; BAOAB puts `O` between them, which changes `p` to `p′` at
///   fixed `q`, so a step leaves `h²(p² − p′²)/8` in the books — a kick of either sign, of size
///   `h² √(γh)`. Over `T/h` steps that is a random walk of size `√(γT) h²`. **So the bath makes
///   the books a random walk** where NVE has a band, and one start's departure measures its order
///   with a spread seven times NVE's. The per-step change is not asserted under a bath: at these
///   steps its `h³` and `h^(5/2)` parts are comparable, and the measured order (2.9 ± 0.3) is
///   neither.
///
/// The bounds are earned from how much each measured order varies over sixteen starts
/// (`the_orders_measured_over_many_starts`; the numbers are in the CHANGELOG): about four standard
/// deviations, and above the worst of the 32. A bath whose work was not counted would leave the
/// books moving by the bath's work, which does not depend on the step: a departure order near 0.
/// **The Langevin bound, 0.8, cannot tell `h` from `h²` with any confidence**, and is not where the
/// bath's books are checked: `a_step_leaves_exactly_its_closed_form_in_the_books` checks each step's
/// change against its exact form on a harmonic well, to rounding.
#[test]
fn the_energy_error_falls_as_the_power_of_the_step_it_should() {
    let nve = orders(7, false);
    println!(
        "NVE: departure orders {:.3}, {:.3}; per-step orders {:.3}, {:.3}",
        nve[0], nve[1], nve[2], nve[3]
    );
    for order in &nve[..2] {
        assert!((order - 2.0).abs() < 0.15, "NVE departure order {order}");
    }
    for order in &nve[2..] {
        assert!((order - 3.0).abs() < 0.3, "NVE per-step order {order}");
    }
    let bath = orders(7, true);
    println!(
        "Langevin: departure orders {:.3}, {:.3}; per-step orders {:.3}, {:.3} (not asserted)",
        bath[0], bath[1], bath[2], bath[3]
    );
    for order in &bath[..2] {
        assert!(
            (order - 2.0).abs() < 0.8,
            "Langevin books departure order {order}"
        );
    }
}

/// **The books balance: kinetic + potential − thermostat work moves only by the integrator's
/// error**, while the bath moves energy in and out by far more. Started at 30 K in a 300 K bath, so
/// the bath has to put in about `(3N/2) k_B · 270 K` and its kinetic share — tens of kcal/mol —
/// while the books move by the NVE band at this step, a tenth of a kcal/mol.
#[test]
fn the_thermostat_work_balances_the_books() {
    let (c, ff, start) = relaxed_aspirin();
    let mut at = start.clone();
    let mut md = MolecularDynamics::for_elements(elements(&c))
        .with_bath(Bath::Langevin {
            temperature: 300.0,
            friction: 5e12,
            seed: 0xB00C,
        })
        .thermalised(&at, 30.0, 0xB00D);
    md.prepare(&ff, &at);
    let opening = md.ledger();
    let e0 = opening.get(quantity::ENERGY).unwrap();
    let mut worst = 0.0f64;
    for _ in 0..4000 {
        md.step(&ff, &mut at, 0.5 * FS);
        let e = md.ledger().get(quantity::ENERGY).unwrap();
        worst = worst.max((e - e0).abs());
    }
    let work = md.thermostat_work();
    println!(
        "bath work {:.3} kcal/mol, books moved at most {:.5} kcal/mol",
        work / KCAL_PER_MOL,
        worst / KCAL_PER_MOL
    );
    let kt = kb() * 300.0;
    assert!(
        work > 10.0 * kt,
        "warming from 30 K means putting energy in: {work}"
    );
    // The books move by the integrator's band: at 0.5 fs NVE aspirin's energy stays within
    // 0.115 kcal/mol of its start at 238 K (the step-size table), and the band scales with the
    // temperature, so 0.15 at 300 K. A kT (0.596 kcal/mol) is four times that, and the bath's work
    // is fifty kT.
    assert!(worst < kt, "the books moved {worst}, a kT is {kt}");
    audit("aspirin", &opening, &md.ledger(), 1e-2).expect("the books balance to the band");
}

// ---------------------------------------------------------------------------------------------
// Determinism and frozen atoms
// ---------------------------------------------------------------------------------------------

fn langevin_aspirin(seed: u64) -> (ForceField, Vec<[f64; 3]>, MolecularDynamics) {
    let (c, ff) = aspirin();
    let (_, at) = relaxed(&ff, &positions(&c));
    let md = MolecularDynamics::for_elements(elements(&c))
        .with_bath(Bath::Langevin {
            temperature: 300.0,
            friction: 1e12,
            seed,
        })
        .thermalised(&at, 300.0, seed);
    (ff, at, md)
}

/// **One seed, one run, to the bit — however it is cut up.** 300 steps in one call, as three calls
/// of 100, and as 300 calls of one, give the same positions, velocities and books bit for bit; a
/// different seed gives a different run.
#[test]
fn a_run_is_its_seed_however_it_is_chunked() {
    let (ff, start, md) = langevin_aspirin(42);
    let h = 0.5 * FS;
    let bits = |at: &[[f64; 3]], md: &MolecularDynamics| {
        let mut b: Vec<u64> = at.iter().flatten().map(|v| v.to_bits()).collect();
        b.extend(md.velocities().iter().flatten().map(|v| v.to_bits()));
        b.push(md.thermostat_work().to_bits());
        b.push(md.half_step_kinetic_energy().to_bits());
        b
    };
    let (mut a, mut ma) = (start.clone(), md.clone());
    ma.run(&ff, &mut a, h, 300);
    let (mut b, mut mb) = (start.clone(), md.clone());
    for _ in 0..3 {
        mb.run(&ff, &mut b, h, 100);
    }
    let (mut c, mut mc) = (start.clone(), md.clone());
    for _ in 0..300 {
        mc.step(&ff, &mut c, h);
    }
    assert_eq!(bits(&a, &ma), bits(&b, &mb), "three chunks of 100");
    assert_eq!(bits(&a, &ma), bits(&c, &mc), "300 single steps");
    // A fresh copy of the same start and seed: the same bits again.
    let (ff2, start2, md2) = langevin_aspirin(42);
    let (mut d, mut md2) = (start2, md2);
    md2.run(&ff2, &mut d, h, 300);
    assert_eq!(bits(&a, &ma), bits(&d, &md2), "rebuilt from scratch");

    let (_, start3, md3) = langevin_aspirin(43);
    let (mut e, mut md3) = (start3, md3);
    md3.run(&ff, &mut e, h, 300);
    assert_ne!(bits(&a, &ma), bits(&e, &md3), "another seed is another run");
}

/// **A frozen atom does not move, to the bit, and its velocity is zero** — and freezing one atom
/// does not change the kicks any other receives. Two atoms in independent wells: atom 0's
/// trajectory is the same bits whether atom 1 is frozen or free, because each kick is keyed by
/// `step · N + atom` with `N` counting frozen atoms too.
#[test]
fn frozen_atoms_stay_put_and_do_not_shift_the_noise() {
    let (c, ff, start) = relaxed_aspirin();
    let n = c.atoms().len();
    let mask: Vec<bool> = (0..n).map(|i| i % 3 == 0).collect();
    let mut at = start.clone();
    let mut md = MolecularDynamics::for_elements(elements(&c))
        .with_bath(Bath::Langevin {
            temperature: 300.0,
            friction: 1e12,
            seed: 5,
        })
        .with_frozen(mask.clone())
        .thermalised(&at, 300.0, 5);
    md.run(&ff, &mut at, 0.5 * FS, 500);
    let mut moved = 0;
    for i in 0..n {
        if mask[i] {
            assert_eq!(
                at[i].map(f64::to_bits),
                start[i].map(f64::to_bits),
                "frozen atom {i} moved"
            );
            assert_eq!(md.velocities()[i], [0.0; 3]);
        } else if at[i] != start[i] {
            moved += 1;
        }
    }
    assert_eq!(
        moved,
        mask.iter().filter(|f| !**f).count(),
        "every free atom moved"
    );

    let wells = Wells {
        k: 200.0,
        origins: vec![[0.0; 3], [1e-9, 0.0, 0.0]],
    };
    let run = |frozen: bool| {
        let mut at = vec![[1e-11, 0.0, 0.0], [1e-9, 2e-11, 0.0]];
        let mut md = MolecularDynamics::new(vec![Element::C.mass(); 2])
            .with_bath(Bath::Langevin {
                temperature: 300.0,
                friction: 1e13,
                seed: 9,
            })
            .with_frozen(vec![false, frozen])
            .with_velocities(vec![[0.0; 3]; 2]);
        md.run(&wells, &mut at, 1.0 * FS, 400);
        (at[0].map(f64::to_bits), at[1].map(f64::to_bits))
    };
    let (free, _) = run(false);
    let (beside_frozen, frozen_one) = run(true);
    assert_eq!(
        free, beside_frozen,
        "freezing atom 1 changed atom 0's noise"
    );
    assert_eq!(frozen_one, [1e-9, 2e-11, 0.0].map(f64::to_bits));
}

// ---------------------------------------------------------------------------------------------
// Initial velocities and rigid motion
// ---------------------------------------------------------------------------------------------

/// **A thermalised isolated molecule has no net momentum and no angular momentum**, and velocity
/// Verlet keeps both: its energy is invariant under translation and rotation, so the scheme
/// conserves them. It then counts `3N − 6` degrees of freedom — `3N − 5` for a linear molecule.
#[test]
fn thermalising_removes_the_rigid_motion_and_nve_keeps_it_removed() {
    let (c, ff, start) = relaxed_aspirin();
    let n = c.atoms().len();
    let mut at = start.clone();
    let mut md = MolecularDynamics::for_elements(elements(&c)).thermalised(&at, 300.0, 3);
    assert_eq!(md.degrees_of_freedom(), 3 * n - 6);
    // Each atom has its own element's mass, in the atoms' order: a permutation keeps the total.
    assert_eq!(md.masses().len(), n);
    for (i, a) in c.atoms().iter().enumerate() {
        assert_eq!(md.masses()[i], a.element.mass(), "{}", a.name);
    }
    // The scales the residues are judged on: one atom's thermal momentum, and its angular
    // momentum at the molecule's size.
    let p_scale = (Element::C.mass() * kb() * 300.0).sqrt();
    let l_scale = p_scale * 3.0 * ANGSTROM;
    let check = |md: &MolecularDynamics, at: &[[f64; 3]], what: &str| {
        let p = md.momentum();
        let l = md.angular_momentum(at);
        let (pn, ln) = (dot(p, p).sqrt() / p_scale, dot(l, l).sqrt() / l_scale);
        println!("{what}: |P| = {pn:.2e}, |L| = {ln:.2e} of one atom's");
        (pn, ln)
    };
    let (p0, l0) = check(&md, &at, "thermalised");
    assert!(p0 < 1e-13 && l0 < 1e-13);
    md.run(&ff, &mut at, 0.5 * FS, 2000);
    let (p1, l1) = check(&md, &at, "after 1 ps of NVE");
    assert!(p1 < 1e-10 && l1 < 1e-8, "{p1} {l1}");

    // Carbon dioxide on a line: two rotations, so 3N − 5 = 4, and no angular momentum about any
    // axis afterwards.
    let co2 = entry(
        &[
            ("O1", "O", [-1.16, 0.0, 0.0]),
            ("C1", "C", [0.0, 0.0, 0.0]),
            ("O2", "O", [1.16, 0.0, 0.0]),
        ],
        &[("O1", "C1", "DOUB", false), ("C1", "O2", "DOUB", false)],
    );
    let at = positions(&co2);
    let md = MolecularDynamics::for_elements(elements(&co2)).thermalised(&at, 300.0, 4);
    assert_eq!(md.degrees_of_freedom(), 4);
    let (p, l) = check(&md, &at, "CO₂ thermalised");
    assert!(p < 1e-13 && l < 1e-13);
    // One atom: nothing is left to move.
    let one = MolecularDynamics::new(vec![Element::C.mass()]).thermalised(&[[0.0; 3]], 300.0, 4);
    assert_eq!(one.degrees_of_freedom(), 0);
    assert_eq!(one.velocities()[0], [0.0; 3]);
}

/// **A translated and rotated start gives the translated and rotated trajectory**, to rounding —
/// in NVE only. Under a bath the kicks are vectors drawn per atom in the lab frame, so rotating the
/// start does not rotate the noise, and the trajectories are different samples of one ensemble.
#[test]
fn a_moved_start_gives_the_moved_trajectory() {
    let (c, ff, start) = relaxed_aspirin();
    let md = MolecularDynamics::for_elements(elements(&c)).thermalised(&start, 300.0, 21);
    // A rotation by 0.7 rad about (1, 2, 3)/√14, by Rodrigues' formula, and a 3 Å shift.
    let axis = {
        let a = [1.0, 2.0, 3.0];
        let l = dot(a, a).sqrt();
        [a[0] / l, a[1] / l, a[2] / l]
    };
    let (s, co) = 0.7f64.sin_cos();
    let rotate = |v: [f64; 3]| {
        let kxv = cross(axis, v);
        let kv = dot(axis, v);
        [0, 1, 2].map(|i| v[i] * co + kxv[i] * s + axis[i] * kv * (1.0 - co))
    };
    let shift = [3.0 * ANGSTROM, -2.0 * ANGSTROM, 1.0 * ANGSTROM];
    let moved = |p: [f64; 3]| {
        let r = rotate(p);
        [r[0] + shift[0], r[1] + shift[1], r[2] + shift[2]]
    };
    let mut a = start.clone();
    let mut ma = md.clone();
    let mut b: Vec<[f64; 3]> = start.iter().map(|p| moved(*p)).collect();
    let mut mb = md
        .clone()
        .with_velocities(md.velocities().iter().map(|v| rotate(*v)).collect());
    let mut worst = Vec::new();
    for _ in 0..200 {
        ma.step(&ff, &mut a, 0.5 * FS);
        mb.step(&ff, &mut b, 0.5 * FS);
        let w = a
            .iter()
            .zip(&b)
            .map(|(p, q)| {
                let d = sub(moved(*p), *q);
                dot(d, d).sqrt()
            })
            .fold(0.0f64, f64::max);
        worst.push(w / ANGSTROM);
    }
    println!(
        "moved-frame departure: {:.2e} Å after 1 step, {:.2e} after 50, {:.2e} after 200",
        worst[0], worst[49], worst[199]
    );
    assert!(worst[199] < 1e-9, "{:e}", worst[199]);
}
