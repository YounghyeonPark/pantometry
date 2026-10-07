//! **The alchemical free-energy engine against closed forms**: Bennett's acceptance ratio and
//! thermodynamic integration on two harmonic wells whose free-energy difference is
//! `(3/2) k_BT ln(k₁/k₀)` exactly; the Boresch correction against a quadrature of its own
//! integral and that integral against Cartesian quadratures that have no Jacobian to get wrong;
//! the soft core's exact limits and its derivatives against finite differences; one Lennard-Jones
//! particle decoupled from a fixed atom, against its configurational integral done by quadrature;
//! a schedule run backwards; and a campaign that is its seed however it is cut up.
//!
//! The long measurements — BAR's and TI's uncertainties calibrated over many MD seeds — are
//! ignored by default and print what they find:
//! `cargo test -p pantometry-forcefield --release --test the_free_energy_against_closed_forms -- --ignored --nocapture`.

use pantometry_core::Rng;
use pantometry_forcefield::alchemy::SoftCoreTerm;
use pantometry_forcefield::boresch::STANDARD_VOLUME;
use pantometry_forcefield::energy::Pair;
use pantometry_forcefield::free_energy::{bar, bar_correlated, bennett_chain, refine_schedule};
use pantometry_forcefield::uff::KCAL_PER_MOL;
use pantometry_forcefield::{
    Alchemical, Boresch, Coupling, CrossPair, Decoupling, Element, Estimate, Lambda, Protocol,
    Quadrature, SoftCore, UffType, Windows,
};
use pantometry_units::BOLTZMANN;
use std::f64::consts::PI;

const ANGSTROM: f64 = 1e-10;
const FS: f64 = 1e-15;
const KELVIN: f64 = 300.0;
const DEG: f64 = PI / 180.0;

fn kt() -> f64 {
    BOLTZMANN.to_si() * KELVIN
}

fn kcal(j: f64) -> f64 {
    j / KCAL_PER_MOL
}

/// Composite Simpson on `[a, b]` with `n` (even) intervals.
fn simpson(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    assert!(n & 1 == 0, "an even number of intervals");
    let h = (b - a) / n as f64;
    let mut s = f(a) + f(b);
    for i in 1..n {
        s += f(a + i as f64 * h) * if i % 2 == 1 { 4.0 } else { 2.0 };
    }
    s * h / 3.0
}

/// `x` wrapped to (−π, π].
fn wrap(x: f64) -> f64 {
    let mut d = x % (2.0 * PI);
    if d > PI {
        d -= 2.0 * PI;
    } else if d <= -PI {
        d += 2.0 * PI;
    }
    d
}

// ---------------------------------------------------------------------------------------------
// Two harmonic wells
// ---------------------------------------------------------------------------------------------

/// One atom in an isotropic well whose stiffness is `k₀ + λ (k₁ − k₀)`, the λ being the state's
/// `restraint` component: `U(λ) = (1 − λ) U₀ + λ U₁`, linear in λ, and
/// `F₁ − F₀ = (3/2) k_BT ln(k₁/k₀)` exactly.
struct Wells {
    k0: f64,
    k1: f64,
}

impl Wells {
    fn stiffness(&self, lambda: Lambda) -> f64 {
        self.k0 + lambda.restraint * (self.k1 - self.k0)
    }

    fn exact(&self, kt: f64) -> f64 {
        1.5 * kt * (self.k1 / self.k0).ln()
    }

    /// `⟨∂U/∂λ⟩` at `λ`, exactly: `½ (k₁ − k₀) ⟨|x|²⟩ = (3/2) k_BT (k₁ − k₀)/k(λ)`.
    fn mean_gradient(&self, kt: f64, lambda: f64) -> f64 {
        1.5 * kt * (self.k1 - self.k0) / (self.k0 + lambda * (self.k1 - self.k0))
    }
}

impl Alchemical for Wells {
    fn len(&self) -> usize {
        1
    }

    fn energy_and_forces(&self, at: &[[f64; 3]], lambda: Lambda, forces: &mut [[f64; 3]]) -> f64 {
        let k = self.stiffness(lambda);
        let x = at[0];
        forces[0] = [-k * x[0], -k * x[1], -k * x[2]];
        0.5 * k * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
    }

    fn coupling(&self, at: &[[f64; 3]], lambda: Lambda) -> Coupling {
        let x = at[0];
        let x2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
        Coupling {
            energy: 0.5 * self.stiffness(lambda) * x2,
            gradient: [0.5 * (self.k1 - self.k0) * x2, 0.0, 0.0],
        }
    }
}

/// The restraint component from 0 to 1 in `n` equal steps, the others at 1.
fn restraint_schedule(n: usize) -> Vec<Lambda> {
    (0..=n)
        .map(|k| Lambda::new(k as f64 / n as f64, 1.0, 1.0))
        .collect()
}

/// **BAR's estimate and its variance are calibrated on independent samples.** Two 3-D harmonic
/// wells with `k₁/k₀ = 4`, in reduced units, sampled exactly — Gaussian positions from
/// `Rng::for_index`, no dynamics — 400 from state 0 and 700 from state 1, so that `M ≠ 0`, for
/// 1600 seeds. `Δf = (3/2) ln 4` exactly. **The band's width is set by the seeds, not the samples
/// per seed**, so many small campaigns buy a band that sees a σ̂ wrong by 8%.
///
/// - The mean of the estimates is within four standard errors of the mean of `Δf` (BAR's
///   finite-sample bias is `O(1/N)`, and is inside that).
/// - The z-scores `(Δf̂ − Δf)/σ̂`, with Shirts's variance and with the delta method's, have mean 0
///   within `4/√1600` and variance 1 within `4 √(2/1599)` = 0.14, the spread of a sample variance
///   of 1600 unit normals. A variance formula without its `− (1/N_F + 1/N_R)` overstates σ̂ and
///   falls below the band.
/// - For two identical states the overlap is ½ and the estimate zero, exactly.
#[test]
fn bar_is_calibrated_on_independent_samples() {
    let (k0, k1) = (1.0f64, 4.0f64);
    let exact = 1.5 * (k1 / k0).ln();
    // Unequal, so that `M = ln(N_F/N_R)` is not zero and its sign is seen.
    let (nf, nr, seeds) = (400usize, 700usize, 1600u64);
    let draw = |rng: &mut Rng, k: f64| {
        let s = 1.0 / k.sqrt();
        let x = [rng.gaussian() * s, rng.gaussian() * s, rng.gaussian() * s];
        0.5 * (k1 - k0) * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
    };
    let (mut z, mut zc) = (Vec::new(), Vec::new());
    let mut estimates = Vec::new();
    for seed in 0..seeds {
        let mut f = Rng::for_index(0xBA2, 2 * seed);
        let mut r = Rng::for_index(0xBA2, 2 * seed + 1);
        let forward: Vec<f64> = (0..nf).map(|_| draw(&mut f, k0)).collect();
        let reverse: Vec<f64> = (0..nr).map(|_| draw(&mut r, k1)).collect();
        let b = bar(&forward, &reverse);
        let c = bar_correlated(&forward, &reverse);
        assert_eq!(b.delta.to_bits(), c.delta.to_bits());
        z.push((b.delta - exact) / b.error());
        zc.push((c.delta - exact) / c.error());
        estimates.push(b.delta);
    }
    let m = seeds as f64;
    let moments = |z: &[f64]| {
        let mean = z.iter().sum::<f64>() / m;
        (
            mean,
            z.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (m - 1.0),
        )
    };
    let (mean_z, var_z) = moments(&z);
    let (mean_zc, var_zc) = moments(&zc);
    let mean = estimates.iter().sum::<f64>() / m;
    let sd = (estimates.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (m - 1.0)).sqrt();
    println!(
        "BAR on exact samples: mean Δf {mean:.5} against {exact:.5} (sd {sd:.4}); Shirts's \
         variance: z mean {mean_z:+.3}, variance {var_z:.3}; the delta method's: z mean \
         {mean_zc:+.3}, variance {var_zc:.3}"
    );
    assert!(
        (mean - exact).abs() < 4.0 * sd / m.sqrt(),
        "mean {mean} against {exact}"
    );
    for (mz, vz) in [(mean_z, var_z), (mean_zc, var_zc)] {
        assert!(mz.abs() < 4.0 / m.sqrt(), "z mean {mz}");
        assert!(
            (vz - 1.0).abs() < 4.0 * (2.0 / (m - 1.0)).sqrt(),
            "z variance {vz}"
        );
    }

    // Bisection's last midpoint is a rounding from zero, not zero.
    let same: Vec<f64> = vec![0.0; 300];
    let b = bar(&same, &same);
    assert!(b.delta.abs() < 1e-15, "{}", b.delta);
    assert_eq!(b.overlap, 0.5);
}

/// One window of the wells below: candidate `g` and, per sample, the reduced energy at every
/// candidate.
struct Well {
    g: usize,
    energies: Vec<Vec<f64>>,
}

/// **Windows inserted where neighbours overlap too little recover a closed form**: seventeen
/// candidate states, isotropic harmonic wells in three dimensions with `k_g = k₀ 2^{5g/8}`, and a
/// schedule of the two ends alone. Each window is 2000 exact, independent samples of its own well,
/// with the reduced energy at every candidate, so that a window inserted later needs nothing of
/// the others run again. The overlap of two wells `c` apart in stiffness is 0.0010 for `c = 1024`
/// (the ends), 0.048 for 32 and 0.23 for `√32`, by quadrature of `p₀p₁/(p₀ + p₁)`: so the rule,
/// [`refine_schedule`] below 0.1, must insert candidate 8 in the first round, 4 and 12 in the
/// second, and stop — each decision several standard errors of an overlap from the threshold.
/// Asserted: exactly that schedule, those runs in that order and those findings; every interval
/// the rule left at 0.1 or above, and every one it split below; **BAR along the refined chain
/// within 4σ of `(3/2) ln(k₁₆/k₀) = 15 ln 2`**; and with one round, the schedule `0, 8, 16` with
/// both its intervals still below.
#[test]
fn windows_inserted_where_neighbours_overlap_little_recover_a_closed_form() {
    let k = |g: usize| (5.0 * g as f64 / 8.0 * 2f64.ln()).exp();
    let run = |g: usize| {
        let mut rng = Rng::for_index(0x1E5, g as u64);
        let width = 1.0 / k(g).sqrt();
        let energies = (0..2000)
            .map(|_| {
                let x = [0, 1, 2].map(|_| rng.gaussian() * width);
                let r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
                (0..17).map(|c| 0.5 * k(c) * r2).collect()
            })
            .collect();
        Well { g, energies }
    };
    let delta =
        |w: &Well, to: usize| -> Vec<f64> { w.energies.iter().map(|e| e[to] - e[w.g]).collect() };
    let overlap = |a: &Well, b: &Well| {
        let reverse: Vec<f64> = delta(b, a.g).iter().map(|x| -x).collect();
        bar(&delta(a, b.g), &reverse).overlap
    };
    let refine = |rounds: usize| {
        let mut schedule = vec![0, 16];
        let mut records = vec![run(0), run(16)];
        let mut ran = Vec::new();
        let found = refine_schedule(&mut schedule, &mut records, 0.1, rounds, overlap, |g| {
            ran.push(g);
            run(g)
        });
        (schedule, records, ran, found)
    };
    let (schedule, records, ran, found) = refine(5);
    for f in &found {
        println!(
            "round {}: {:?} at {:.4}, inserted {:?}",
            f.round, f.between, f.overlap, f.inserted
        );
    }
    assert_eq!(schedule, [0, 4, 8, 12, 16]);
    assert_eq!(ran, [8, 4, 12]);
    let rows: Vec<(usize, [usize; 2], Option<usize>)> = found
        .iter()
        .map(|f| (f.round, f.between, f.inserted))
        .collect();
    assert_eq!(
        rows,
        [
            (0, [0, 16], Some(8)),
            (1, [0, 8], Some(4)),
            (1, [8, 16], Some(12))
        ]
    );
    assert!(found.iter().all(|f| f.overlap < 0.1));
    for (r, g) in records.iter().zip(&schedule) {
        assert_eq!(r.g, *g);
    }
    for w in records.windows(2) {
        let o = overlap(&w[0], &w[1]);
        println!("{} → {}: overlap {o:.4}", w[0].g, w[1].g);
        assert!(o >= 0.1);
    }
    let chain: Vec<(Vec<f64>, Vec<f64>)> = (0..records.len())
        .map(|j| {
            let prev = if j > 0 {
                delta(&records[j], records[j - 1].g)
            } else {
                Vec::new()
            };
            let next = if j + 1 < records.len() {
                delta(&records[j], records[j + 1].g)
            } else {
                Vec::new()
            };
            (prev, next)
        })
        .collect();
    let f = bennett_chain(&chain);
    let exact = 15.0 * 2f64.ln();
    println!(
        "BAR along the refined chain {:.4} ± {:.4} against {exact:.4}",
        f.value, f.error
    );
    assert!((f.value - exact).abs() < 4.0 * f.error);
    let (schedule, records, ran, found) = refine(1);
    assert_eq!(schedule, [0, 8, 16]);
    assert_eq!(ran, [8]);
    assert_eq!(found.len(), 1);
    for w in records.windows(2) {
        assert!(
            overlap(&w[0], &w[1]) < 0.1,
            "a round fewer leaves them below"
        );
    }
}

/// **TI and BAR from MD windows recover two harmonic wells' `(3/2) k_BT ln(k₁/k₀)`.** One carbon
/// atom, `k₁ = 4k₀`, `ω₀h = 0.1` and `ω₁h = 0.2` at 1 fs, in a bath of `γ = 2ω₀`; eleven windows
/// in equal steps of λ. BAOAB samples a harmonic well's positions exactly at any stable step
/// (`tests/the_dynamics.rs`), so the windows' only errors are statistical and the quadrature's.
///
/// - **The quadrature's bias is computed, not hoped small**: each rule applied to the exact
///   `⟨∂U/∂λ⟩ = (3/2) k_BT (k₁ − k₀)/k(λ)`. Simpson's is asserted below a tenth of its σ̂, and
///   Simpson's TI then within 4σ̂ of the exact ΔF. The trapezoid's bias, +0.0105 `k_BT`, is 31
///   times Simpson's, and its TI is held within 4σ̂ of what the rule gives on the exact integrand.
/// - **BAR** within 4σ̂ of the exact ΔF, on every sample, its σ̂ the delta method's window by
///   window ([`bennett_chain`]).
#[test]
fn ti_and_bar_from_windows_recover_two_harmonic_wells() {
    let m = Element::C.mass();
    let h = 1.0 * FS;
    let omega0 = 0.1 / h;
    let k0 = m * omega0 * omega0;
    let wells = Wells { k0, k1: 4.0 * k0 };
    let kt = kt();
    let exact = wells.exact(kt);
    let schedule = restraint_schedule(10);
    let protocol = Protocol {
        time_step: h,
        temperature: KELVIN,
        friction: 2.0 * omega0,
        equilibration: 200,
        stride: 10,
        samples: 2000,
        seed: 0x7E115,
    };
    let mut w = Windows::new(
        &wells,
        schedule.clone(),
        &[[0.0; 3]],
        vec![m],
        vec![false],
        protocol,
    );
    assert!(w.run(&wells, u64::MAX));

    let quadrature_of_exact = |weights: &[f64]| -> f64 {
        schedule
            .iter()
            .zip(weights)
            .map(|(l, wk)| wk * wells.mean_gradient(kt, l.restraint))
            .sum()
    };
    let n = schedule.len();
    let trapezoid: Vec<f64> = (0..n)
        .map(|k| if k == 0 || k == n - 1 { 0.05 } else { 0.1 })
        .collect();
    let simpson_w: Vec<f64> = (0..n)
        .map(|k| {
            0.1 / 3.0
                * if k == 0 || k == n - 1 {
                    1.0
                } else if k % 2 == 1 {
                    4.0
                } else {
                    2.0
                }
        })
        .collect();
    let trap_expected = quadrature_of_exact(&trapezoid);
    let simp_expected = quadrature_of_exact(&simpson_w);

    let trap = w.thermodynamic_integration(Quadrature::Trapezoid).unwrap();
    let simp = w.thermodynamic_integration(Quadrature::Simpson).unwrap();
    let b = w.bennett_total();
    println!(
        "exact ΔF {:.5} kT; trapezoid {:.5} ± {:.5} (rule's bias {:+.5}); Simpson {:.5} ± {:.5} \
         (bias {:+.6}); BAR {:.5} ± {:.5}",
        exact / kt,
        trap.value / kt,
        trap.error / kt,
        (trap_expected - exact) / kt,
        simp.value / kt,
        simp.error / kt,
        (simp_expected - exact) / kt,
        b.value / kt,
        b.error / kt
    );
    assert!((simp_expected - exact).abs() < 0.1 * simp.error);
    assert!((simp.value - exact).abs() < 4.0 * simp.error);
    assert!((trap.value - trap_expected).abs() < 4.0 * trap.error);
    assert!((b.value - exact).abs() < 4.0 * b.error);
    // The intervals one by one are the same estimates the total sums, bit for bit.
    let by_interval = kt * w.bennett().iter().map(|x| x.delta).sum::<f64>();
    assert_eq!(by_interval.to_bits(), b.value.to_bits());
    // Each window's mean gradient against its own closed form, at 4σ̂.
    for (win, l) in w.windows().iter().zip(&schedule) {
        let g: Vec<f64> = win.samples().iter().map(|s| s.gradient[0]).collect();
        let e = Estimate::of(&g);
        let want = wells.mean_gradient(kt, l.restraint);
        assert!(
            (e.mean - want).abs() < 4.0 * e.error,
            "λ = {}: {} against {want}",
            l.restraint,
            e.mean
        );
    }
}

/// **The uncertainties over many MD seeds** — what calibrates the autocorrelation-corrected
/// errors of TI and of BAR, thinning included, and measures the covariance that summing the
/// adjacent intervals' variances leaves out. Printed, not asserted.
#[test]
#[ignore = "a measurement: 200 campaigns, run with --release -- --ignored --nocapture"]
fn the_uncertainties_over_many_md_seeds_measured() {
    let m = Element::C.mass();
    let h = 1.0 * FS;
    let omega0 = 0.1 / h;
    let k0 = m * omega0 * omega0;
    let wells = Wells { k0, k1: 4.0 * k0 };
    let kt = kt();
    let exact = wells.exact(kt);
    for (label, friction, stride) in [
        ("γ = 2ω₀, stride 10", 2.0 * omega0, 10u64),
        ("γ = 0.2ω₀, stride 2 (correlated)", 0.2 * omega0, 2),
    ] {
        let (mut zt, mut zb) = (Vec::new(), Vec::new());
        for seed in 0..200u64 {
            let protocol = Protocol {
                time_step: h,
                temperature: KELVIN,
                friction,
                equilibration: 500,
                stride,
                samples: 1000,
                seed,
            };
            let mut w = Windows::new(
                &wells,
                restraint_schedule(4),
                &[[0.0; 3]],
                vec![m],
                vec![false],
                protocol,
            );
            w.run(&wells, u64::MAX);
            let t = w.thermodynamic_integration(Quadrature::Simpson).unwrap();
            let b = w.bennett_total();
            // Simpson's bias on five states, from the exact integrand, removed.
            let s: Vec<f64> = (0..5)
                .map(|k| wells.mean_gradient(kt, k as f64 / 4.0))
                .collect();
            let simp = 0.25 / 3.0 * (s[0] + 4.0 * s[1] + 2.0 * s[2] + 4.0 * s[3] + s[4]);
            zt.push((t.value - simp) / t.error);
            zb.push((b.value - exact) / b.error);
        }
        let stats = |z: &[f64]| {
            let n = z.len() as f64;
            let mean = z.iter().sum::<f64>() / n;
            (
                mean,
                z.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0),
            )
        };
        let (mt, vt) = stats(&zt);
        let (mb, vb) = stats(&zb);
        println!(
            "{label}: TI z mean {mt:+.3} variance {vt:.3}; BAR z mean {mb:+.3} variance {vb:.3} \
             (a variance of 1 ± {:.2} is calibrated)",
            4.0 * (2.0f64 / 199.0).sqrt()
        );
    }
}

// ---------------------------------------------------------------------------------------------
// The Boresch restraint
// ---------------------------------------------------------------------------------------------

/// Six atoms, `c – b – a ⋯ A – B – C`, with every hinge angle well away from straight: indices
/// 0, 1, 2 are `a`, `b`, `c` and 3, 4, 5 are `A`, `B`, `C`.
fn six_atoms() -> Vec<[f64; 3]> {
    let p = |x: f64, y: f64, z: f64| [x * ANGSTROM, y * ANGSTROM, z * ANGSTROM];
    vec![
        p(0.0, 0.0, 0.0),
        p(-1.5, 0.2, 0.0),
        p(-2.1, 1.5, 0.3),
        p(2.4, 4.2, 0.7),
        p(3.1, 5.3, -0.4),
        p(4.5, 5.1, -0.9),
    ]
}

/// `[K_r, K_θ, K_θ, K_φ, K_φ, K_φ]` in SI from kcal mol⁻¹ Å⁻² and kcal mol⁻¹ rad⁻².
fn constants(kr: f64, ktheta: f64, kphi: f64) -> [f64; 6] {
    let r = kr * KCAL_PER_MOL / (ANGSTROM * ANGSTROM);
    let t = ktheta * KCAL_PER_MOL;
    let p = kphi * KCAL_PER_MOL;
    [r, t, t, p, p, p]
}

/// `∫ e^(−β K (x − x₀)²/2) w(x) dx` over `[a, b]`, Simpson: the restraint's Gaussian typed out,
/// used only where [`Boresch::energy`] cannot be evaluated — the tails beyond a coordinate's range
/// — and for the test's own radial tether. Zero for an empty interval.
fn gaussian(k: f64, x0: f64, a: f64, b: f64, w: impl Fn(f64) -> f64) -> f64 {
    if b <= a {
        return 0.0;
    }
    let beta = 1.0 / kt();
    simpson(
        |x| (-beta * 0.5 * k * (x - x0) * (x - x0)).exp() * w(x),
        a,
        b,
        4_000,
    )
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn unit(a: [f64; 3]) -> [f64; 3] {
    let n = norm(a);
    [a[0] / n, a[1] / n, a[2] / n]
}

fn angle_at(p: [f64; 3], v: [f64; 3], q: [f64; 3]) -> f64 {
    let (u, w) = (sub(p, v), sub(q, v));
    norm(cross(u, w)).atan2(u[0] * w[0] + u[1] * w[1] + u[2] * w[2])
}

/// The point `p₄` with `|p₃p₄| = length`, `∠(p₂, p₃, p₄) = angle` and the dihedral
/// `(p₁, p₂, p₃, p₄) = torsion` in [`pantometry_forcefield::minimise::dihedral`]'s convention:
/// the natural extension of the reference frame. Checked by [`eq6_factors`]'s round trip.
fn place(
    p1: [f64; 3],
    p2: [f64; 3],
    p3: [f64; 3],
    length: f64,
    angle: f64,
    torsion: f64,
) -> [f64; 3] {
    let bc = unit(sub(p3, p2));
    let n = unit(cross(sub(p2, p1), bc));
    let m = cross(n, bc);
    let d = [
        -length * angle.cos(),
        length * angle.sin() * torsion.cos(),
        TORSION_SIGN * length * angle.sin() * torsion.sin(),
    ];
    [0, 1, 2].map(|c| p3[c] + d[0] * bc[c] + d[1] * m[c] + d[2] * n[c])
}

/// Which way [`place`] turns for a positive dihedral, fixed by the round trip in [`eq6_factors`].
const TORSION_SIGN: f64 = 1.0;

/// The six atoms with the receptor's where `at0` has them, the ligand's `|AB|`, `|BC|` and
/// `∠(A, B, C)` as `at0` has them, and the six Boresch coordinates `xi`.
fn configuration(at0: &[[f64; 3]], xi: [f64; 6]) -> Vec<[f64; 3]> {
    let (a, b, c) = (at0[0], at0[1], at0[2]);
    let (ab, bc) = (norm(sub(at0[4], at0[3])), norm(sub(at0[5], at0[4])));
    let abc = angle_at(at0[3], at0[4], at0[5]);
    let la = place(c, b, a, xi[0], xi[1], xi[3]);
    let lb = place(b, a, la, ab, xi[2], xi[4]);
    let lc = place(a, la, lb, bc, abc, xi[5]);
    vec![a, b, c, la, lb, lc]
}

/// Eq 6's six one-dimensional integrals (Clark et al. 2023) over each coordinate's own range,
/// **of `e^(−β U)` with `U` [`Boresch::energy`] itself**, at a configuration with that one
/// coordinate moved and the other five at their reference: `[I_r, I_θA, I_θB, I_φA, I_φB, I_φC]`,
/// with the `r²` and `sin θ` of the measure. `r` from 20σ below `r₀` (below which nothing is left
/// to integrate, and `A` would meet `a`), `θ` over [0, π], `φ` over the period about its
/// reference. The configurations are checked to give back their coordinates first.
fn eq6_factors(b: &Boresch, at0: &[[f64; 3]]) -> [f64; 6] {
    let x0 = b.reference;
    let mut moved = x0;
    moved[0] *= 1.07;
    moved[1] -= 0.2;
    moved[2] += 0.15;
    moved[3] += 2.5;
    moved[4] -= 0.4;
    moved[5] += 0.9;
    for xi in [x0, moved] {
        let back = b.coordinates(&configuration(at0, xi));
        assert!((back[0] / xi[0] - 1.0).abs() < 1e-13, "{back:?} {xi:?}");
        for k in 1..6 {
            assert!(wrap(back[k] - xi[k]).abs() < 1e-12, "{k}: {back:?} {xi:?}");
        }
    }
    let beta = 1.0 / kt();
    let s = (kt() / b.force_constants[0]).sqrt();
    let range = |k: usize| match k {
        0 => (x0[0] - 20.0 * s, x0[0] + 40.0 * s),
        1 | 2 => (0.0, PI),
        _ => (x0[k] - PI, x0[k] + PI),
    };
    let weight = |k: usize, x: f64| match k {
        0 => x * x,
        1 | 2 => x.sin(),
        _ => 1.0,
    };
    [0, 1, 2, 3, 4, 5].map(|k| {
        let (lo, hi) = range(k);
        simpson(
            |x| {
                // At θ = 0 or π the next frame is undefined, and the measure is zero.
                let w = weight(k, x);
                if w == 0.0 || (k == 1 || k == 2) && (x <= 0.0 || x >= PI) {
                    return 0.0;
                }
                let mut xi = x0;
                xi[k] = x;
                (-beta * b.energy(&configuration(at0, xi))).exp() * w
            },
            lo,
            hi,
            4_000,
        )
    })
}

/// What eq 6's integrals leave out of the Gaussians over the whole line, by quadrature of the
/// typed Gaussian over the excluded ranges: `r` below the range of [`eq6_factors`], `θ` outside
/// [0, π] (where `sin θ` is negative), and `φ` more than π from its reference.
fn eq6_tails(b: &Boresch) -> [f64; 6] {
    let k = b.force_constants;
    let x0 = b.reference;
    let sigma = |i: usize| (kt() / k[i]).sqrt();
    [0, 1, 2, 3, 4, 5].map(|i| match i {
        0 => gaussian(
            k[0],
            x0[0],
            x0[0] - 40.0 * sigma(0),
            x0[0] - 20.0 * sigma(0),
            |r| r * r,
        ),
        1 | 2 => {
            gaussian(k[i], x0[i], x0[i] - 40.0 * sigma(i), 0.0, f64::sin)
                + gaussian(k[i], x0[i], PI, x0[i] + 40.0 * sigma(i), f64::sin)
        }
        _ => {
            gaussian(k[i], x0[i], x0[i] - 40.0 * sigma(i), x0[i] - PI, |_| 1.0)
                + gaussian(k[i], x0[i], x0[i] + PI, x0[i] + 40.0 * sigma(i), |_| 1.0)
        }
    })
}

/// `8π²`, the orientations' measure, by the same quadrature: eq 6's orientational integral with
/// every force constant zero, `∫₀^π sin θ dθ ∫ dφ ∫ dφ`.
fn orientations() -> f64 {
    let one = simpson(|_| 1.0, -PI, PI, 4_000);
    simpson(f64::sin, 0.0, PI, 4_000) * one * one
}

/// **The standard volume is a litre per mole**: `10⁻³ m³ / N_A`, N_A = 6.022 140 76 × 10²³ exactly
/// (SI 2019), is 1660.539 067 Å³.
#[test]
fn the_standard_volume_is_a_litre_per_mole() {
    let v = STANDARD_VOLUME / (ANGSTROM * ANGSTROM * ANGSTROM);
    assert!((v - 1_660.539_067_173_847).abs() < 1e-9, "{v}");
}

/// **Boresch's closed form is eq 6 with two approximations, and each is checked exactly.** At
/// `K_r` = 10 kcal mol⁻¹ Å⁻² and every angle's 10 kcal mol⁻¹ rad⁻² — the size Mobley et al. use —
/// eq 6 is integrated by quadrature of `e^(−β U)` with `U` the restraint's own
/// [`Boresch::energy`], factor by factor ([`eq6_factors`]), with `8π²` from the same quadrature
/// ([`orientations`]), and set against:
///
/// - **eq 7 with its approximations undone in closed form**
///   ([`Boresch::release_free_energy_extended`]). That is exact for the Gaussians over the whole
///   line, so it is held to eq 6's integrals **plus the excluded tails, by quadrature**
///   ([`eq6_tails`]), to 1e-9 `k_BT`; measured, 4e-13. The tails are 2.3e-5 of the integral at
///   these angles and the `r²` correction `ln(1 + k_BT/K_r r₀²)` is 2.5e-3 `k_BT`, four and six
///   orders above the bound. (The bound was once `e^(−min(θ₀, π − θ₀)²/2σ²)`, 2.5e-3, a hundred
///   times the tails it bounded, and caught a missing `r²` correction by 2%.)
/// - **eq 7 itself** ([`Boresch::release_free_energy`]), printed.
///
/// Because the integrand is [`Boresch::energy`], a restraint whose `½` went missing — in the
/// energy, the forces and the dihedrals alike, so that it is still its own gradient — is a
/// restraint of `2K` integrated against a closed form of `K`, and fails here.
#[test]
fn the_boresch_release_is_its_own_integral_in_closed_form() {
    let at = six_atoms();
    let b = Boresch::at(&at, [0, 1, 2], [3, 4, 5], constants(10.0, 10.0, 10.0));
    let kt = kt();
    let eight_pi2 = orientations();
    assert!(
        (eight_pi2 / (8.0 * PI * PI) - 1.0).abs() < 1e-12,
        "{eight_pi2}"
    );
    let i = eq6_factors(&b, &at);
    let tails = eq6_tails(&b);
    let z: f64 = i.iter().product();
    let line: f64 = i.iter().zip(&tails).map(|(a, t)| a + t).product();
    let numeric = -kt * (eight_pi2 * STANDARD_VOLUME / z).ln();
    let over_the_line = -kt * (eight_pi2 * STANDARD_VOLUME / line).ln();
    let extended = b.release_free_energy_extended(KELVIN);
    let eq7 = b.release_free_energy(KELVIN);
    let theta = b.reference.map(|x| x / DEG);
    println!(
        "θ_A0 {:.1}°, θ_B0 {:.1}°, r₀ {:.3} Å: eq 6 by quadrature {:.6} kcal/mol, with the tails \
         {:.9}; eq 7 extended {:.9}; eq 7 {:.6}; the tails {:.3e} of the integral; extended less \
         quadrature {:+.2e} kT",
        theta[1],
        theta[2],
        b.reference[0] / ANGSTROM,
        kcal(numeric),
        kcal(over_the_line),
        kcal(extended),
        kcal(eq7),
        line / z - 1.0,
        (extended - over_the_line) / kt
    );
    assert!(
        (extended - over_the_line).abs() <= 1e-9 * kt,
        "eq 7 extended {extended} against eq 6 and its tails {over_the_line}"
    );
    assert!(numeric < 0.0, "releasing a restraint gains entropy");
}

/// **Eq 6's Jacobian, from Cartesian coordinates**: the integral over `A`'s position of
/// `e^(−β U)`, with `U` [`Boresch::energy`] and `B` and `C` carried on their reference `θ_B`,
/// `φ_B` and `φ_C`, so that only the distance, `θ_A` and `φ_A` terms are not zero, done on a grid
/// in `A`'s Cartesian coordinates, against `I_r I_θA I_φA`; and the integral over `B`'s position
/// (with `C` carried on its reference `φ_C`) of `e^(−β U)` times a stiff radial tether on
/// `|B − A|`, against `I_θB I_φB` times the tether's own radial integral. A Cartesian grid has
/// unit Jacobian, so the `r² sin θ_A` and `sin θ_B` of eq 6 are measured, not assumed.
///
/// The constants are stiff (40 kcal mol⁻¹ Å⁻², 200 kcal mol⁻¹ rad⁻²) so that a grid of a few
/// tens of thousands of points covers ±8σ at a spacing of σ/2 in a frame aligned with the
/// restraint. The trapezoid rule on a smooth integrand decaying faster than any power converges
/// faster than any power; the grid is repeated at 0.8 of the spacing, and both are asserted
/// within 1e-9 of the one-dimensional product — eight orders below a missing `sin θ` at these
/// angles.
#[test]
fn the_boresch_integral_has_its_jacobian() {
    let at0 = six_atoms();
    let b = Boresch::at(&at0, [0, 1, 2], [3, 4, 5], constants(40.0, 200.0, 200.0));
    let kt = kt();
    let beta = 1.0 / kt;
    let i = eq6_factors(&b, &at0);
    let k = b.force_constants;
    let x0 = b.reference;
    // A frame aligned with `centre − origin`, and a grid in it of ±8σ at spacing `scale` σ.
    let grid = |origin: [f64; 3],
                centre: [f64; 3],
                along: f64,
                across: f64,
                scale: f64,
                f: &dyn Fn([f64; 3]) -> f64| {
        let e0 = unit(sub(centre, origin));
        let e1 = unit(cross(e0, [0.3, 0.1, 1.0]));
        let e2 = cross(e0, e1);
        let (hu, hv) = (scale * along, scale * across);
        let (nu, nv) = ((8.0 / scale).ceil() as i64, (8.0 / scale).ceil() as i64);
        let mut sum = 0.0;
        for iu in -nu..=nu {
            for iv in -nv..=nv {
                for iw in -nv..=nv {
                    let (u, v, w) = (iu as f64 * hu, iv as f64 * hv, iw as f64 * hv);
                    let p = [0, 1, 2].map(|c| centre[c] + u * e0[c] + v * e1[c] + w * e2[c]);
                    sum += f(p);
                }
            }
        }
        sum * hu * hv * hv
    };
    let sr = (kt / k[0]).sqrt();
    let st = (kt / k[1]).sqrt();
    let (ab, bc) = (norm(sub(at0[4], at0[3])), norm(sub(at0[5], at0[4])));
    let abc = angle_at(at0[3], at0[4], at0[5]);
    // The position: A over space, B and C placed on their reference θ_B, φ_B, φ_C.
    let position = |p: [f64; 3]| {
        let mut at = at0.clone();
        at[3] = p;
        at[4] = place(at[1], at[0], at[3], ab, x0[2], x0[4]);
        at[5] = place(at[0], at[3], at[4], bc, abc, x0[5]);
        (-beta * b.energy(&at)).exp()
    };
    let want_position = i[0] * i[1] * i[3];
    // The orientation: B over space about A, C on its reference φ_C, and a radial tether
    // |B − A| about its reference length, which is the test's and not the restraint's.
    let kd = k[0];
    let orientation = |p: [f64; 3]| {
        let mut at = at0.clone();
        at[4] = p;
        at[5] = place(at[0], at[3], at[4], bc, abc, x0[5]);
        let dd = norm(sub(p, at0[3])) - ab;
        (-beta * (b.energy(&at) + 0.5 * kd * dd * dd)).exp()
    };
    let radial = gaussian(kd, ab, 0.0, ab + 40.0 * sr, |r| r * r);
    let want_orientation = i[2] * i[4] * radial;
    for scale in [0.5, 0.4] {
        let gp = grid(at0[0], at0[3], sr, st * x0[0], scale, &position);
        let go = grid(at0[3], at0[4], sr, st * ab, scale, &orientation);
        let (ep, eo) = (gp / want_position - 1.0, go / want_orientation - 1.0);
        println!("spacing {scale} σ: position {ep:+.2e}, orientation {eo:+.2e}");
        assert!(ep.abs() < 1e-9, "position: {ep:e}");
        assert!(eo.abs() < 1e-9, "orientation: {eo:e}");
    }
}

/// Central differences of `f` at step `h` and at `h/2` against `analytic`, both within `bound`.
/// A central difference's error is `h² f‴/6` plus a rounding of about `ε |f| / h`; each call
/// states a step at which both are orders below its bound, and an analytic derivative that is
/// wrong — a missing chain-rule factor, a sign — is off by a part in a few, not in 10⁷.
fn converges(f: impl Fn(f64) -> f64, x: f64, h: f64, analytic: f64, bound: f64) {
    let fd = |h: f64| (f(x + h) - f(x - h)) / (2.0 * h);
    let (e1, e2) = ((fd(h) - analytic).abs(), (fd(h / 2.0) - analytic).abs());
    assert!(
        e1 <= bound && e2 <= bound,
        "at {x:e}: analytic {analytic:e}, errors {e1:e} and {e2:e}, bound {bound:e}"
    );
}

/// **The restraint's force is minus its gradient**, scaled: on the six atoms moved off the
/// reference, `add_forces` at scale 0.7 against central differences of 0.7 × `energy` in every
/// coordinate of every atom, the bound 1e-7 of the largest force component; the energy it
/// returns is `energy`'s, unscaled, bit for bit; and `energy` is `Σ ½ K (ξ − ξ₀)²` of its own
/// [`Boresch::coordinates`], each dihedral's difference wrapped, written out here, to 4 ε — the
/// `½` that [`Boresch::release_free_energy`]'s `(2π k_BT)³` assumes.
#[test]
fn the_restraint_force_is_its_gradient() {
    let at0 = six_atoms();
    let b = Boresch::at(&at0, [0, 1, 2], [3, 4, 5], constants(10.0, 10.0, 10.0));
    let mut at = at0.clone();
    at[3][1] += 0.4 * ANGSTROM;
    at[4][2] -= 0.3 * ANGSTROM;
    at[5][0] += 0.5 * ANGSTROM;
    at[1][2] += 0.2 * ANGSTROM;
    let scale = 0.7;
    let mut forces = vec![[0.0; 3]; 6];
    let e = b.add_forces(&at, &mut forces, scale);
    assert_eq!(e.to_bits(), b.energy(&at).to_bits());
    let x = b.coordinates(&at);
    let by_hand: f64 = (0..6)
        .map(|k| {
            let d = if k < 3 {
                x[k] - b.reference[k]
            } else {
                wrap(x[k] - b.reference[k])
            };
            0.5 * b.force_constants[k] * d * d
        })
        .sum();
    assert!(
        (e - by_hand).abs() <= 4.0 * f64::EPSILON * by_hand,
        "{e:e} against ½ K Δξ² {by_hand:e}"
    );
    let largest = forces.iter().flatten().fold(0.0f64, |m, f| m.max(f.abs()));
    for atom in 0..6 {
        for c in 0..3 {
            let f = |x: f64| {
                let mut p = at.clone();
                p[atom][c] = x;
                scale * b.energy(&p)
            };
            converges(
                f,
                at[atom][c],
                1e-4 * ANGSTROM,
                -forces[atom][c],
                1e-7 * largest,
            );
        }
    }
    // Net force zero: the restraint is internal to the six atoms.
    for c in 0..3 {
        let net: f64 = forces.iter().map(|f| f[c]).sum();
        assert!(net.abs() < 1e-12 * largest, "net {net:e}");
    }
}

/// **The rule picks the pair furthest from straight**: of two receptor triples, one with `θ_A`
/// near 175°, the other's is chosen, and nothing is chosen when no `r₀` is in range.
#[test]
fn the_anchor_rule_avoids_a_straight_hinge() {
    let mut at = six_atoms();
    // A seventh atom on the far side of `a` from `A`, nearly in line: θ_A ≈ 175° for (a, 6, c).
    let a = at[0];
    let la = at[3];
    let u = [a[0] - la[0], a[1] - la[1], a[2] - la[2]];
    at.push([
        a[0] + 0.3 * u[0] + 0.05 * ANGSTROM,
        a[1] + 0.3 * u[1],
        a[2] + 0.3 * u[2],
    ]);
    let k = constants(10.0, 10.0, 10.0);
    let chosen = Boresch::choose(
        &at,
        &[[0, 6, 2], [0, 1, 2]],
        &[[3, 4, 5]],
        (2.0 * ANGSTROM, 8.0 * ANGSTROM),
        k,
    )
    .unwrap();
    assert_eq!(chosen.receptor, [0, 1, 2]);
    let straight = Boresch::at(&at, [0, 6, 2], [3, 4, 5], k).hinge_angles(&at)[0];
    assert!(straight > 170.0 * DEG, "{}", straight / DEG);
    assert!(Boresch::choose(
        &at,
        &[[0, 1, 2]],
        &[[3, 4, 5]],
        (6.0 * ANGSTROM, 8.0 * ANGSTROM),
        k
    )
    .is_none());
}

// ---------------------------------------------------------------------------------------------
// The soft core
// ---------------------------------------------------------------------------------------------

fn carbon_pair() -> Pair {
    Pair::new([0, 1], UffType::C3, UffType::C3)
}

/// **At λ = 1 the soft core is UFF's pair to the bit, at λ = 0 it is nothing, exactly, and for
/// λ < 1 it is finite at r = 0 with the value its formula gives there.** The pair's radial
/// derivative at λ = 1 is the documented `12 D (s⁶ − s¹²)/x` over `r`, written here in the same
/// order. At `r = 0`, `r_sc⁶ = α σ⁶ (1 − λ)` and `σ⁶ = x⁶/2`, so `s = 2/(α(1 − λ))` and
/// `U = λ D (s² − 2s)`, a closed form held to 4 ε.
#[test]
fn the_soft_core_is_the_pair_at_one_and_nothing_at_zero() {
    let p = carbon_pair();
    let sc = SoftCore::default();
    assert_eq!((sc.alpha, sc.power), (0.5, 1.0));
    for r in [2.5, 3.0, 3.851, 4.5, 8.0, 15.0] {
        let r = r * ANGSTROM;
        let t: SoftCoreTerm = sc.at(&p, r, 1.0);
        assert_eq!(t.energy.to_bits(), p.energy(r).to_bits(), "r {r:e}");
        let s = p.distance / r;
        let s6 = (s * s) * (s * s) * (s * s);
        let s12 = s6 * s6;
        let de_dr = 12.0 * p.well * (s6 - s12) / r;
        assert_eq!(t.de_dr_over_r.to_bits(), (de_dr / r).to_bits());
    }
    for r in [0.0, 0.5, 2.0, 3.851, 6.0] {
        let t = sc.at(&p, r * ANGSTROM, 0.0);
        assert_eq!(t.energy, 0.0);
        assert_eq!(t.de_dr_over_r, 0.0);
        assert!(t.de_dlambda.is_finite());
    }
    for lambda in [0.05, 0.3, 0.5, 0.9, 0.999] {
        let t = sc.at(&p, 0.0, lambda);
        let s = 2.0 / (sc.alpha * (1.0 - lambda));
        let want = lambda * p.well * (s * s - 2.0 * s);
        assert!(t.energy.is_finite() && t.de_dlambda.is_finite());
        assert_eq!(t.de_dr_over_r, 0.0, "no force at the centre, by symmetry");
        assert!(
            (t.energy - want).abs() <= 4.0 * f64::EPSILON * want.abs(),
            "λ {lambda}: {} against {want}",
            t.energy
        );
    }
}

/// **The soft core's derivatives are its energy's**, in λ and in r, by central differences within
/// 1e-7 of the term's scale (see [`converges`]), on a grid of λ and r that covers the core, the well and the tail;
/// and at λ = 1 the λ-derivative against a one-sided second-order difference from below.
/// Dropping the chain rule's `λ dV/dq ∂q/∂λ` from `∂U/∂λ` is off by more than the energy itself
/// at small r.
#[test]
fn the_soft_cores_derivatives_are_its_energys() {
    let p = carbon_pair();
    let sc = SoftCore::default();
    for lambda in [0.02, 0.2, 0.5, 0.8, 0.98] {
        for r in [0.3, 1.5, 2.8, 3.4, 4.2, 7.0] {
            let r = r * ANGSTROM;
            let t = sc.at(&p, r, lambda);
            let scale = t.energy.abs().max(t.de_dlambda.abs()).max(1e-3 * p.well);
            converges(
                |l| sc.at(&p, r, l).energy,
                lambda,
                // `r_sc⁶ = r⁶ + α σ⁶ (1 − λ)`: in the core the energy's scale in λ is `1 − λ`.
                1e-5 * (1.0 - lambda),
                t.de_dlambda,
                1e-7 * scale,
            );
            converges(
                |x| sc.at(&p, x, lambda).energy,
                r,
                1e-5 * ANGSTROM,
                t.de_dr_over_r * r,
                1e-7 * scale / ANGSTROM,
            );
        }
    }
    for r in [2.8, 3.851, 5.0] {
        let r = r * ANGSTROM;
        let t = sc.at(&p, r, 1.0);
        let f = |l: f64| sc.at(&p, r, l).energy;
        let one_sided = |h: f64| (3.0 * f(1.0) - 4.0 * f(1.0 - h) + f(1.0 - 2.0 * h)) / (2.0 * h);
        let (e1, e2) = (
            (one_sided(1e-5) - t.de_dlambda).abs(),
            (one_sided(5e-6) - t.de_dlambda).abs(),
        );
        let scale = t.energy.abs().max(t.de_dlambda.abs());
        assert!(e1 < 1e-6 * scale && e2 < 1e-6 * scale, "{e1:e} {e2:e}");
    }
}

// ---------------------------------------------------------------------------------------------
// One Lennard-Jones particle, decoupled from a fixed atom
// ---------------------------------------------------------------------------------------------

/// A carbon atom held at the origin (atom 0, frozen) and a second (atom 1) tethered by
/// `½ k_t |x − c|²` to a point `d₀` = 3.2 Å from it, inside the pair's repulsive wall: the
/// soft-core pair decoupled, the tether λ-independent.
struct Tethered {
    decoupling: Decoupling,
    centre: [f64; 3],
    stiffness: f64,
}

const TETHER_DISTANCE: f64 = 3.2 * ANGSTROM;
const TETHER_WIDTH: f64 = 0.3 * ANGSTROM;

fn tethered() -> Tethered {
    let pair = carbon_pair();
    Tethered {
        decoupling: Decoupling::from_pairs(
            2,
            vec![false, true],
            vec![CrossPair {
                pair,
                charges: [0.0, 0.0],
            }],
        ),
        centre: [TETHER_DISTANCE, 0.0, 0.0],
        stiffness: kt() / (TETHER_WIDTH * TETHER_WIDTH),
    }
}

impl Alchemical for Tethered {
    fn len(&self) -> usize {
        2
    }

    fn energy_and_forces(&self, at: &[[f64; 3]], lambda: Lambda, forces: &mut [[f64; 3]]) -> f64 {
        let mut e = self.decoupling.energy_and_forces(at, lambda, forces);
        for k in 0..3 {
            let d = at[1][k] - self.centre[k];
            e += 0.5 * self.stiffness * d * d;
            forces[1][k] -= self.stiffness * d;
        }
        e
    }

    fn coupling(&self, at: &[[f64; 3]], lambda: Lambda) -> Coupling {
        self.decoupling.coupling(at, lambda)
    }
}

/// `F(decoupled) − F(coupled)` for [`tethered`], by quadrature of its configurational integral:
/// about the fixed atom, with the tether's centre on the polar axis, the angular integral is
/// `∫ sin θ e^(β k_t r d₀ cos θ) dθ = 2 sinh(a)/a`, `a = β k_t r d₀`, so
/// `Z₁ = 2π ∫ r² e^(−βV(r)) [e^(−βk_t(r − d₀)²/2) − e^(−βk_t(r + d₀)²/2)] / (β k_t r d₀) dr`,
/// and `Z₀ = (2π k_BT/k_t)^(3/2)`. Simpson from 1.5 Å, where `e^(−βV)` is below 10⁻³⁰⁰, to
/// `d₀ + 14σ`. The same quadrature with `V = 0` is held to `Z₀` at 1e-12, which checks it.
fn decoupling_free_energy_by_quadrature(t: &Tethered) -> f64 {
    let beta = 1.0 / kt();
    let (kk, d0) = (t.stiffness, TETHER_DISTANCE);
    let p = carbon_pair();
    let z = |with_pair: bool, from: f64| {
        simpson(
            |r| {
                if r == 0.0 {
                    return 0.0;
                }
                let v = if with_pair { p.energy(r) } else { 0.0 };
                let shell = ((-beta * kk * 0.5 * (r - d0).powi(2)).exp()
                    - (-beta * kk * 0.5 * (r + d0).powi(2)).exp())
                    / (beta * kk * r * d0);
                2.0 * PI * r * r * (-beta * v).exp() * shell
            },
            from,
            d0 + 14.0 * TETHER_WIDTH,
            40_000,
        )
    };
    let z0 = (2.0 * PI / (beta * kk)).powf(1.5);
    let z0_quadrature = z(false, 0.0);
    assert!(
        (z0_quadrature / z0 - 1.0).abs() < 1e-12,
        "{z0_quadrature} {z0}"
    );
    let z1 = z(true, 1.5 * ANGSTROM);
    -kt() * (z0 / z1).ln()
}

/// λ_v from `from` to `to` in `n` equal steps, no restraint, charges off.
fn vdw_schedule(from: f64, to: f64, n: usize) -> Vec<Lambda> {
    (0..=n)
        .map(|k| Lambda::new(0.0, 0.0, from + (to - from) * k as f64 / n as f64))
        .collect()
}

fn lj_protocol(seed: u64) -> Protocol {
    Protocol {
        time_step: 1.0 * FS,
        temperature: KELVIN,
        friction: 5e13,
        equilibration: 500,
        stride: 20,
        samples: 800,
        seed,
    }
}

fn lj_windows(schedule: Vec<Lambda>, seed: u64) -> (Tethered, Windows) {
    let t = tethered();
    let start = [[0.0; 3], t.centre];
    let mut w = Windows::new(
        &t,
        schedule,
        &start,
        vec![Element::C.mass(); 2],
        vec![true, false],
        lj_protocol(seed),
    );
    assert!(w.run(&t, u64::MAX));
    (t, w)
}

/// **Decoupling one Lennard-Jones particle from a fixed atom, against its configurational
/// integral.** The particle is tethered 3.2 Å from the fixed atom, inside the pair's repulsive
/// wall (its zero is at 3.43 Å), so the coupled state is pushed out and decoupling releases it:
/// `−0.506 k_BT` by quadrature, about fifteen of this run's σ̂ from zero. Twenty-one windows of
/// λ_v from 1 to 0, 800 samples every 20 fs in each.
///
/// BAR and Simpson's TI are each asserted within 4σ̂ of the quadrature; the trapezoid is printed.
/// Neither the soft-core path nor the integrator is in the reference: both endpoints are exactly
/// the plain pair and nothing ([`the_soft_core_is_the_pair_at_one_and_nothing_at_zero`]). The
/// step, 1 fs at `ωh` ≈ 0.06 on the wall's curvature, could leave a configurational bias; over
/// 48 seeds ([`the_particle_over_many_seeds_measured`]) BAR's mean is 0.6 standard errors of the
/// mean from the quadrature at 1 fs and 0.15 at 0.5 fs, so none is resolved.
#[test]
fn decoupling_a_particle_matches_its_configurational_integral() {
    let (t, w) = lj_windows(vdw_schedule(1.0, 0.0, 20), 0x11A);
    let exact = decoupling_free_energy_by_quadrature(&t);
    let kt = kt();
    let b = w.bennett_total();
    let simp = w.thermodynamic_integration(Quadrature::Simpson).unwrap();
    let trap = w.thermodynamic_integration(Quadrature::Trapezoid).unwrap();
    let worst_overlap = w
        .bennett()
        .iter()
        .map(|x| x.overlap)
        .fold(f64::INFINITY, f64::min);
    println!(
        "quadrature {:.4} kT; BAR {:.4} ± {:.4}; TI Simpson {:.4} ± {:.4}, trapezoid {:.4} ± \
         {:.4}; smallest overlap {worst_overlap:.3}",
        exact / kt,
        b.value / kt,
        b.error / kt,
        simp.value / kt,
        simp.error / kt,
        trap.value / kt,
        trap.error / kt
    );
    assert!(
        exact < -0.5 * kt,
        "the reference is not small: {}",
        exact / kt
    );
    assert!((b.value - exact).abs() < 4.0 * b.error);
    assert!((simp.value - exact).abs() < 4.0 * simp.error);
}

/// **The particle over many seeds and two steps**: the mean of BAR's and TI's estimates against
/// the quadrature, their spread against the σ̂ each reports, at 1 fs and 0.5 fs. Printed.
#[test]
#[ignore = "a measurement: 2 × 48 campaigns, run with --release -- --ignored --nocapture"]
fn the_particle_over_many_seeds_measured() {
    let exact = decoupling_free_energy_by_quadrature(&tethered());
    let kt = kt();
    for (h, stride, samples) in [(1.0, 20u64, 800usize), (0.5, 40, 800), (1.0, 100, 800)] {
        let t = tethered();
        let start = [[0.0; 3], t.centre];
        let (mut b, mut ti, mut eb, mut et) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for seed in 0..48u64 {
            let protocol = Protocol {
                time_step: h * FS,
                stride,
                samples,
                equilibration: (500.0 / h) as u64,
                ..lj_protocol(seed)
            };
            let mut w = Windows::new(
                &t,
                vdw_schedule(1.0, 0.0, 20),
                &start,
                vec![Element::C.mass(); 2],
                vec![true, false],
                protocol,
            );
            w.run(&t, u64::MAX);
            let x = w.bennett_total();
            let y = w.thermodynamic_integration(Quadrature::Simpson).unwrap();
            b.push(x.value / kt);
            eb.push(x.error / kt);
            ti.push(y.value / kt);
            et.push(y.error / kt);
        }
        let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        let sd = |v: &[f64]| {
            let m = mean(v);
            (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (v.len() as f64 - 1.0)).sqrt()
        };
        println!(
            "h {h} fs, stride {stride}: quadrature {:.4}; BAR mean {:.4} (± {:.4} of the mean), \
             spread {:.4} against σ̂ {:.4}; TI mean {:.4} (± {:.4}), spread {:.4} against σ̂ {:.4}",
            exact / kt,
            mean(&b),
            sd(&b) / (b.len() as f64).sqrt(),
            sd(&b),
            mean(&eb),
            mean(&ti),
            sd(&ti) / (ti.len() as f64).sqrt(),
            sd(&ti),
            mean(&et)
        );
    }
}

/// **The cycle closes**: the same particle coupled again, λ_v from 0 to 1, from another seed,
/// is minus the decoupling, within four of the two runs' combined standard errors, by BAR and by
/// Simpson's TI.
#[test]
fn a_schedule_run_backwards_closes_the_cycle() {
    let (_, forward) = lj_windows(vdw_schedule(1.0, 0.0, 20), 0xF0);
    let (_, backward) = lj_windows(vdw_schedule(0.0, 1.0, 20), 0xBAC);
    for (label, f, b) in [
        ("BAR", forward.bennett_total(), backward.bennett_total()),
        (
            "TI",
            forward
                .thermodynamic_integration(Quadrature::Simpson)
                .unwrap(),
            backward
                .thermodynamic_integration(Quadrature::Simpson)
                .unwrap(),
        ),
    ] {
        let gap = f.value + b.value;
        let sigma = (f.error * f.error + b.error * b.error).sqrt();
        println!(
            "{label}: forward {:.4}, backward {:.4}, sum {:.4} ± {:.4} kT",
            f.value / kt(),
            b.value / kt(),
            gap / kt(),
            sigma / kt()
        );
        assert!(gap.abs() < 4.0 * sigma, "{label}");
    }
}

// ---------------------------------------------------------------------------------------------
// Determinism
// ---------------------------------------------------------------------------------------------

fn sample_bits(w: &Windows) -> Vec<u64> {
    w.windows()
        .iter()
        .flat_map(|win| {
            win.samples().iter().flat_map(|s| {
                [
                    s.step,
                    s.gradient[0].to_bits(),
                    s.gradient[1].to_bits(),
                    s.gradient[2].to_bits(),
                    s.to_previous.to_bits(),
                    s.to_next.to_bits(),
                ]
            })
        })
        .chain(
            w.windows()
                .iter()
                .flat_map(|win| win.positions().iter().flatten().map(|x| x.to_bits())),
        )
        .collect()
}

/// **A campaign is its seed, however it is cut up, and resumes where it stopped.** Five windows
/// of the tethered particle run in one call, in calls of seven steps (which divide neither the
/// equilibration nor the stride), and stopped half-way, cloned, and finished: the same bits in
/// every sample and every final position. Each window is the same alone as among the others, and
/// two windows at the same state, from the same start, differ, which only their seeds can make
/// them do.
#[test]
fn a_campaign_is_its_seed_however_it_is_cut_up() {
    let t = tethered();
    let start = [[0.0; 3], t.centre];
    let protocol = Protocol {
        samples: 40,
        ..lj_protocol(0xD1CE)
    };
    let build = || {
        Windows::new(
            &t,
            vdw_schedule(1.0, 0.0, 4),
            &start,
            vec![Element::C.mass(); 2],
            vec![true, false],
            protocol,
        )
    };
    let mut whole = build();
    assert!(whole.run(&t, u64::MAX));
    // Bounded, so that a campaign that never completes in pieces fails instead of hanging.
    let mut pieces = build();
    let calls = protocol.steps_per_window() / 7 + 1;
    for _ in 0..calls {
        pieces.run(&t, 7);
    }
    assert!(
        pieces.is_complete(),
        "not complete in {calls} calls of seven steps"
    );
    assert_eq!(sample_bits(&whole), sample_bits(&pieces));

    let mut half = build();
    half.run(&t, protocol.steps_per_window() / 2);
    assert!(!half.is_complete());
    let mut resumed = half.clone();
    drop(half);
    assert!(resumed.run(&t, u64::MAX));
    assert_eq!(sample_bits(&whole), sample_bits(&resumed));

    let mut alone = pantometry_forcefield::Window::new(
        &t,
        &vdw_schedule(1.0, 0.0, 4),
        3,
        &start,
        vec![Element::C.mass(); 2],
        vec![true, false],
        protocol,
    );
    alone.advance(&t, u64::MAX);
    assert_eq!(alone.samples(), whole.windows()[3].samples());
    // Two windows at the same state, from the same start: only their seeds tell them apart.
    let same = Lambda::new(0.0, 0.0, 0.5);
    let twins = Windows::new(
        &t,
        vec![same, same],
        &start,
        vec![Element::C.mass(); 2],
        vec![true, false],
        Protocol {
            samples: 5,
            ..protocol
        },
    );
    let mut twins = twins;
    assert!(twins.run(&t, u64::MAX));
    assert_ne!(
        twins.windows()[0].positions(),
        twins.windows()[1].positions(),
        "two windows of one campaign drew the same noise"
    );
    // The frozen atom kept its bits in every window.
    for win in whole.windows() {
        assert_eq!(win.positions()[0], [0.0; 3]);
    }
}

/// **On correlated samples the delta method's variance is calibrated, and Shirts's is not.** The
/// two wells of [`bar_is_calibrated_on_independent_samples`], each set drawn as a stationary AR(1)
/// chain — `x_t = φ x_{t−1} + √(1 − φ²) σ ξ_t` per component, so every sample is exactly from its
/// state and successive ones are correlated, `φ = 0.8` — 800 per state, 1600 seeds. `|x|²`'s
/// correlation is `φ²` per step, a statistical inefficiency of `(1 + φ²)/(1 − φ²)` = 4.6.
///
/// [`bar_correlated`]'s z-scores have mean 0 within `4/√1600` and variance 1 within
/// `4√(2/1599)` = 0.14, so a σ̂ 8% wrong either way is outside it; [`bar`]'s, which assumes
/// independence, have a variance above 2.5.
#[test]
fn on_correlated_samples_the_delta_method_is_calibrated() {
    let (k0, k1) = (1.0f64, 4.0f64);
    let exact = 1.5 * (k1 / k0).ln();
    let (n, seeds, phi) = (800usize, 1600u64, 0.8f64);
    let chain = |rng: &mut Rng, k: f64| -> Vec<f64> {
        let s = 1.0 / k.sqrt();
        let fresh = (1.0 - phi * phi).sqrt();
        let mut x = [rng.gaussian() * s, rng.gaussian() * s, rng.gaussian() * s];
        (0..n)
            .map(|_| {
                for c in &mut x {
                    *c = phi * *c + fresh * s * rng.gaussian();
                }
                0.5 * (k1 - k0) * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
            })
            .collect()
    };
    let (mut zs, mut zc) = (Vec::new(), Vec::new());
    for seed in 0..seeds {
        let forward = chain(&mut Rng::for_index(0xC0, 2 * seed), k0);
        let reverse = chain(&mut Rng::for_index(0xC0, 2 * seed + 1), k1);
        let c = bar_correlated(&forward, &reverse);
        // Shirts's variance from the same solve: `1/S − (1/N_F + 1/N_R)`, `S = O₀₁ N_F`, which
        // `bar_is_calibrated_on_independent_samples` holds `bar` to.
        let shirts = 1.0 / (c.overlap * n as f64) - 2.0 / n as f64;
        zs.push((c.delta - exact) / shirts.sqrt());
        zc.push((c.delta - exact) / c.error());
    }
    let m = seeds as f64;
    let moments = |z: &[f64]| {
        let mean = z.iter().sum::<f64>() / m;
        (
            mean,
            z.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (m - 1.0),
        )
    };
    let (ms, vs) = moments(&zs);
    let (mc, vc) = moments(&zc);
    println!(
        "AR(1) samples, φ = {phi}: Shirts's z mean {ms:+.3} variance {vs:.3}; the delta \
         method's z mean {mc:+.3} variance {vc:.3}"
    );
    assert!(mc.abs() < 4.0 / m.sqrt(), "z mean {mc}");
    assert!(
        (vc - 1.0).abs() < 4.0 * (2.0 / (m - 1.0)).sqrt(),
        "z variance {vc}"
    );
    assert!(vs > 2.5, "independence assumed: z variance {vs}");
}

/// **Along a chain of states the total's error includes the covariance between intervals that
/// share a window, and is calibrated.** Five states of the two wells, `k` linear in λ from `k₀`
/// to `4k₀`, each window 150 exact, independent samples carrying its energy differences to both
/// neighbours, as a window's do; 1200 seeds. [`bennett_chain`]'s z-scores against
/// `(3/2) ln 4` have mean 0 within `4/√1200` and variance 1 within `4√(2/1199)` = 0.16. Summing the
/// four intervals' own variances instead — what the window-by-window delta method replaced —
/// gives a variance of z above that band, which is asserted too, so the test can tell the two
/// apart.
#[test]
fn along_a_chain_the_error_includes_the_shared_windows() {
    let (k0, k1) = (1.0f64, 4.0f64);
    let exact = 1.5 * (k1 / k0).ln();
    let (states, n, seeds) = (5usize, 150usize, 1200u64);
    let k = |j: usize| k0 + (k1 - k0) * j as f64 / (states - 1) as f64;
    let (mut zc, mut zs) = (Vec::new(), Vec::new());
    for seed in 0..seeds {
        let chain: Vec<(Vec<f64>, Vec<f64>)> = (0..states)
            .map(|j| {
                let mut rng = Rng::for_index(0xC4A1 + seed, j as u64);
                let s = 1.0 / k(j).sqrt();
                let mut prev = Vec::with_capacity(n);
                let mut next = Vec::with_capacity(n);
                for _ in 0..n {
                    let x = [rng.gaussian() * s, rng.gaussian() * s, rng.gaussian() * s];
                    let x2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
                    prev.push(if j > 0 {
                        0.5 * (k(j - 1) - k(j)) * x2
                    } else {
                        f64::NAN
                    });
                    next.push(if j + 1 < states {
                        0.5 * (k(j + 1) - k(j)) * x2
                    } else {
                        f64::NAN
                    });
                }
                (prev, next)
            })
            .collect();
        let total = bennett_chain(&chain);
        let summed: f64 = (0..states - 1)
            .map(|j| {
                let reverse: Vec<f64> = chain[j + 1].0.iter().map(|x| -x).collect();
                bar_correlated(&chain[j].1, &reverse).variance
            })
            .sum();
        zc.push((total.value - exact) / total.error);
        zs.push((total.value - exact) / summed.sqrt());
    }
    let m = seeds as f64;
    let moments = |z: &[f64]| {
        let mean = z.iter().sum::<f64>() / m;
        (
            mean,
            z.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (m - 1.0),
        )
    };
    let (mc, vc) = moments(&zc);
    let (ms, vs) = moments(&zs);
    let band = 4.0 * (2.0 / (m - 1.0)).sqrt();
    println!(
        "chain of {states}: window by window z mean {mc:+.3} variance {vc:.3}; intervals' \
         variances summed z mean {ms:+.3} variance {vs:.3}; band ±{band:.3}"
    );
    assert!(mc.abs() < 4.0 / m.sqrt(), "z mean {mc}");
    assert!((vc - 1.0).abs() < band, "z variance {vc}");
    assert!(
        vs > 1.0 + band,
        "the summed variances are not distinguished: {vs}"
    );
}

// ---------------------------------------------------------------------------------------------
// The windows' own uncertainties, over many campaigns, on a schedule with a corner
// ---------------------------------------------------------------------------------------------

/// One atom in a well whose stiffness is driven first by `λ_r` and then by `λ_e`:
/// `k = k₀ + λ_r (k_m − k₀) + (1 − λ_e)(k₁ − k_m)`. From `(λ_r, λ_e) = (0, 1)` through the corner
/// `(1, 1)` to `(1, 0)` the stiffness goes `k₀ → k_m → k₁`, and `ΔF = (3/2) k_BT ln(k₁/k₀)`
/// whatever the path.
struct CornerWells {
    k0: f64,
    km: f64,
    k1: f64,
}

impl CornerWells {
    fn stiffness(&self, l: Lambda) -> f64 {
        self.k0 + l.restraint * (self.km - self.k0) + (1.0 - l.electrostatics) * (self.k1 - self.km)
    }

    /// `⟨∇_λ U⟩` at `l`, exactly: `⟨|x|²⟩ = 3 k_BT/k`.
    fn mean_gradient(&self, kt: f64, l: Lambda) -> [f64; 3] {
        let x2 = 3.0 * kt / self.stiffness(l);
        [
            0.5 * (self.km - self.k0) * x2,
            -0.5 * (self.k1 - self.km) * x2,
            0.0,
        ]
    }
}

impl Alchemical for CornerWells {
    fn len(&self) -> usize {
        1
    }

    fn energy_and_forces(&self, at: &[[f64; 3]], lambda: Lambda, forces: &mut [[f64; 3]]) -> f64 {
        let k = self.stiffness(lambda);
        let x = at[0];
        forces[0] = [-k * x[0], -k * x[1], -k * x[2]];
        0.5 * k * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
    }

    fn coupling(&self, at: &[[f64; 3]], lambda: Lambda) -> Coupling {
        let x = at[0];
        let x2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
        Coupling {
            energy: 0.5 * self.stiffness(lambda) * x2,
            gradient: [
                0.5 * (self.km - self.k0) * x2,
                -0.5 * (self.k1 - self.km) * x2,
                0.0,
            ],
        }
    }
}

/// **The windows' own uncertainties are calibrated, TI's and BAR's, on a schedule with a
/// corner.** [`CornerWells`] with `k₀ : k_m : k₁ = 1 : 2 : 4`, three windows —
/// `(λ_r, λ_e) = (0, 1)`, the corner `(1, 1)`, and `(1, 0)` — so the middle window's trapezoid
/// weight spans two λ components and both BAR intervals meet at the corner. 600 campaigns
/// through [`Windows`], each small: 200 samples per window, every 10 steps after 500, at
/// `ω₀h = 0.1` in a bath of `γ = 2ω₀`.
///
/// - **Unbiased**: TI's deviation from the trapezoid applied to the exact `⟨∇_λU⟩` (the rule's
///   bias, +0.17 `k_BT` at a corner this coarse, is the schedule's and is taken out), and BAR's
///   from `(3/2) k_BT ln 4`, each average within four standard errors of zero over the campaigns.
/// - **Calibrated**: the z-scores `(estimate − reference)/σ̂` have variance 1 within
///   `4√(2/599)` = 0.23. An error ×3 or ×½ is a variance of 1/9 or 4, and one 30% wrong either way
///   leaves the band. Their mean is not asserted: σ̂ grows with the estimate here — a campaign
///   that sampled large `|x|²` has both — so z's mean sits at −0.11 for an unbiased estimate, which
///   is what the deviations above check instead. With 100 steps of equilibration from rest the
///   variance read 1.14 for both: the first samples remembered the start.
/// - Simpson's rule refuses the schedule: it is not one straight line.
#[test]
fn the_windows_uncertainties_are_calibrated_around_a_corner() {
    let m = Element::C.mass();
    let h = 1.0 * FS;
    let omega0 = 0.1 / h;
    let k0 = m * omega0 * omega0;
    let wells = CornerWells {
        k0,
        km: 2.0 * k0,
        k1: 4.0 * k0,
    };
    let kt = kt();
    let exact = 1.5 * kt * 4.0f64.ln();
    let schedule = vec![
        Lambda::new(0.0, 1.0, 1.0),
        Lambda::new(1.0, 1.0, 1.0),
        Lambda::new(1.0, 0.0, 1.0),
    ];
    // The trapezoid on the exact integrand: w_k = (λ_{k+1} − λ_{k−1})/2, one-sided at the ends.
    let n = schedule.len();
    let trapezoid: f64 = (0..n)
        .map(|k| {
            let lo = schedule[k.saturating_sub(1)].components();
            let hi = schedule[(k + 1).min(n - 1)].components();
            let g = wells.mean_gradient(kt, schedule[k]);
            (0..3).map(|c| 0.5 * (hi[c] - lo[c]) * g[c]).sum::<f64>()
        })
        .sum();
    let seeds = 600u64;
    let (mut zt, mut zb) = (Vec::new(), Vec::new());
    let (mut dt, mut db) = (Vec::new(), Vec::new());
    for seed in 0..seeds {
        let protocol = Protocol {
            time_step: h,
            temperature: KELVIN,
            friction: 2.0 * omega0,
            equilibration: 500,
            stride: 10,
            samples: 200,
            seed: 0xC0E2 + seed,
        };
        let mut w = Windows::new(
            &wells,
            schedule.clone(),
            &[[0.0; 3]],
            vec![m],
            vec![false],
            protocol,
        );
        assert!(w.run(&wells, u64::MAX));
        if seed == 0 {
            assert!(w.thermodynamic_integration(Quadrature::Simpson).is_none());
        }
        let t = w.thermodynamic_integration(Quadrature::Trapezoid).unwrap();
        let b = w.bennett_total();
        zt.push((t.value - trapezoid) / t.error);
        zb.push((b.value - exact) / b.error);
        dt.push((t.value - trapezoid) / kt);
        db.push((b.value - exact) / kt);
    }
    let count = seeds as f64;
    let moments = |z: &[f64]| {
        let mean = z.iter().sum::<f64>() / count;
        (
            mean,
            z.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (count - 1.0),
        )
    };
    let (mt, vt) = moments(&zt);
    let (mb, vb) = moments(&zb);
    let ((at, st), (ab, sb)) = (moments(&dt), moments(&db));
    let (st, sb) = ((st / count).sqrt(), (sb / count).sqrt());
    let band = 4.0 * (2.0 / (count - 1.0)).sqrt();
    println!(
        "around a corner, {seeds} campaigns: TI z mean {mt:+.3} variance {vt:.3}; BAR z mean \
         {mb:+.3} variance {vb:.3}; band ±{band:.3}; mean deviations TI {at:+.5} ± {st:.5}, BAR \
         {ab:+.5} ± {sb:.5} kT; the trapezoid's own bias {:+.4} kT",
        (trapezoid - exact) / kt
    );
    for (label, dev, se, var) in [("TI", at, st, vt), ("BAR", ab, sb, vb)] {
        assert!(dev.abs() < 4.0 * se, "{label} biased: {dev} ± {se}");
        assert!((var - 1.0).abs() < band, "{label} z variance {var}");
    }
}
