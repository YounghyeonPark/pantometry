//! **Liquid TIP3P at 298 K and 0.997 g/cm³, set beside the literature.** Every number here is
//! reported, not asserted: the potential energy per molecule, the O–O radial distribution
//! function's first peak, the self-diffusion coefficient with Yeh and Hummer's finite-size
//! correction, and the pressure — each with a standard error from its own autocorrelation (or from
//! independent blocks) — and the cost of a step. Release only, and long:
//!
//! ```text
//! cargo test --release -p pantometry-forcefield --test liquid_water_against_the_literature -- --ignored --nocapture
//! ```
//!
//! What is asserted is only what must hold whatever the model does: the waters stay rigid, the
//! run is finite, and the bath's temperature is reached within its error.

#![allow(clippy::needless_range_loop)]

use pantometry_forcefield::uff::KCAL_PER_MOL;
use pantometry_forcefield::water::{self, Settle, WaterBox};
use pantometry_forcefield::{Bath, Estimate, MolecularDynamics, PeriodicForceField};
use pantometry_units::BOLTZMANN;
use std::f64::consts::PI;

const ANGSTROM: f64 = 1e-10;
const FS: f64 = 1e-15;
const TEMPERATURE: f64 = 298.0;
/// The Wigner constant of the simple cubic lattice, Yeh and Hummer's ξ.
const XI: f64 = 2.837_297;
/// TIP3P's shear viscosity at 298 K, Pa s: Yeh and Hummer 2004, (3.08 ± 0.1) × 10⁻⁴ kg m⁻¹ s⁻¹,
/// at 33.00 nm⁻³ with PME.
const ETA_TIP3P: f64 = 3.08e-4;

/// `x` to the whole power `n` by repeated multiplication: `powi`'s precision is unspecified and
/// differs between platforms, so no test input or reference is built with it.
fn pow(x: f64, n: u32) -> f64 {
    (0..n).fold(1.0, |p, _| p * x)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The pressure of the untruncated fluid at `at`, pascals, by the molecular virial: each water's
/// forces taken at its centre of mass, so that the constraint forces, internal to it, drop out;
/// the kinetic part the bath's `N k_B T` exactly; and A&T's tail in place of the virial's.
fn pressure(field: &PeriodicForceField, at: &[[f64; 3]], waters: &[[usize; 3]]) -> f64 {
    let ev = field.evaluate(at);
    let m = water::masses();
    let mut w = ev.virial[0][0] + ev.virial[1][1] + ev.virial[2][2];
    for idx in waters {
        let com = [0, 1, 2]
            .map(|k| (0..3).map(|a| m[a] * at[idx[a]][k]).sum::<f64>() / water::molecular_mass());
        for &i in idx {
            w -= dot(sub(at[i], com), ev.forces[i]);
        }
    }
    let v = field.cell().volume();
    let n = waters.len() as f64;
    (n * BOLTZMANN.to_si() * TEMPERATURE + w / 3.0) / v - field.dispersion_correction() / v
        + field.pressure_correction()
}

/// `D` from oxygen positions saved every `every` seconds: the mean-square displacement over every
/// origin 0.5 ps apart, fitted by least squares from 2 to 10 ps, over 6.
fn diffusion(frames: &[Vec<[f64; 3]>], every: f64) -> f64 {
    let stride = (0.5e-12 / every).round() as usize;
    let (lo, hi) = (
        (2e-12 / every).round() as usize,
        (10e-12 / every).round() as usize,
    );
    let mut msd = vec![0.0; hi + 1];
    let mut count = vec![0usize; hi + 1];
    let mut origin = 0;
    while origin + hi < frames.len() {
        for lag in lo..=hi {
            let (a, b) = (&frames[origin], &frames[origin + lag]);
            let s: f64 = a
                .iter()
                .zip(b)
                .map(|(p, q)| dot(sub(*q, *p), sub(*q, *p)))
                .sum();
            msd[lag] += s / a.len() as f64;
            count[lag] += 1;
        }
        origin += stride;
    }
    let pts: Vec<(f64, f64)> = (lo..=hi)
        .map(|l| (l as f64 * every, msd[l] / count[l] as f64))
        .collect();
    let n = pts.len() as f64;
    let (mx, my) = (
        pts.iter().map(|p| p.0).sum::<f64>() / n,
        pts.iter().map(|p| p.1).sum::<f64>() / n,
    );
    let slope = pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum::<f64>()
        / pts.iter().map(|p| (p.0 - mx) * (p.0 - mx)).sum::<f64>();
    slope / 6.0
}

/// The O–O histogram of one frame, bins of `width` to half the box.
fn histogram(
    field: &PeriodicForceField,
    at: &[[f64; 3]],
    n: usize,
    width: f64,
    bins: usize,
) -> Vec<f64> {
    let cell = field.cell();
    let mut h = vec![0.0; bins];
    for i in 0..n {
        for j in i + 1..n {
            let d = cell.minimum_image(sub(at[3 * i], at[3 * j]));
            let b = (dot(d, d).sqrt() / width) as usize;
            if b < bins {
                h[b] += 1.0;
            }
        }
    }
    // g(r) = pairs in the shell / (N(N−1)/2 × shell volume / V).
    let v = cell.volume();
    for (b, x) in h.iter_mut().enumerate() {
        let (r0, r1) = (b as f64 * width, (b + 1) as f64 * width);
        let shell = 4.0 / 3.0 * PI * (pow(r1, 3) - pow(r0, 3));
        *x /= (n * (n - 1)) as f64 / 2.0 * shell / v;
    }
    h
}

/// The largest departure of any water from its three lengths, in units of `ε` times the water's
/// largest coordinate — the rounding a position is stored to, which grows as a water diffuses
/// away from the origin. `tests/a_rigid_water_against_closed_forms.rs` earns the bound, 16.
fn rigidity(at: &[[f64; 3]], waters: &[[usize; 3]]) -> f64 {
    let l = Settle::tip3p(Vec::new()).lengths();
    let mut worst = 0.0f64;
    for &[o, a, b] in waters {
        let scale = [o, a, b]
            .iter()
            .flat_map(|&i| at[i])
            .fold(0.0f64, |m, x| m.max(x.abs()))
            * f64::EPSILON;
        for (k, (i, j)) in [(o, a), (o, b), (a, b)].into_iter().enumerate() {
            let d = sub(at[i], at[j]);
            worst = worst.max((dot(d, d).sqrt() - l[k]).abs() / scale);
        }
    }
    worst
}

fn liquid(per_side: usize, seed: u64, nvt_ps: f64, nve_ps: f64) {
    let n = per_side.pow(3);
    let mut wb = WaterBox::lattice(per_side, seed);
    let side = wb.cell().lengths()[0];
    let cutoff = (9.0 * ANGSTROM).min(0.5 * side);
    let field = wb.force_field(cutoff, 1e-5);
    let p = field.ewald().parameters();
    let q2: f64 = field.charges().iter().map(|q| q * q).sum();
    println!(
        "\n{n} waters, box {:.3} Å, r_c {:.2} Å, α {:.4} /Å, {} wave vectors, reciprocal bias {:.2e} kcal/mol a water",
        side / ANGSTROM,
        cutoff / ANGSTROM,
        p.alpha * ANGSTROM,
        field.ewald().wave_vectors(),
        p.reciprocal_bias(q2) / KCAL_PER_MOL / n as f64
    );
    let start = std::time::Instant::now();
    let md = wb.equilibrate(&field, TEMPERATURE, 20.0, seed);
    println!(
        "equilibrated 0.5 + 20 ps in {:.0} s",
        start.elapsed().as_secs_f64()
    );
    let mut at = wb.positions().to_vec();
    let waters = wb.waters();
    let dt = 2.0 * FS;

    // NVT, a gentle bath: the energy, the structure and the pressure, and D under the bath.
    let mut md = md.with_bath(Bath::Langevin {
        temperature: TEMPERATURE,
        friction: 1e12,
        seed: seed ^ 0xA7,
    });
    let steps = (nvt_ps * 1e-12 / dt).round() as usize;
    let (width, bins) = (0.02 * ANGSTROM, (0.5 * side / (0.02 * ANGSTROM)) as usize);
    let mut energy = Vec::new();
    let mut temperature = Vec::new();
    let mut press = Vec::new();
    let mut g_frames: Vec<Vec<f64>> = Vec::new();
    let mut o_nvt: Vec<Vec<[f64; 3]>> = Vec::new();
    let t0 = std::time::Instant::now();
    for s in 0..steps {
        md.step(&field, &mut at, dt);
        if s % 10 == 0 {
            energy.push(md.potential_energy().expect("stepped") / n as f64 / KCAL_PER_MOL);
            temperature.push(md.half_step_temperature());
        }
        if s % 50 == 0 {
            g_frames.push(histogram(&field, &at, n, width, bins));
            press.push(pressure(&field, &at, &waters) / 1e5);
            o_nvt.push(waters.iter().map(|w| at[w[0]]).collect());
        }
    }
    let per_step = t0.elapsed().as_secs_f64() / steps as f64;
    assert!(energy.iter().all(|e| e.is_finite()));
    let rigid_nvt = rigidity(&at, &waters);

    // NVE from the last NVT state with the centre-of-mass drift taken out: D without a bath. The
    // temperature is the library's count, 6N − 3 here, not one written in this file.
    let v = md.velocities().to_vec();
    let mut nve: MolecularDynamics = md
        .with_bath(Bath::Isolated)
        .with_velocities(v)
        .without_net_momentum();
    assert_eq!(nve.degrees_of_freedom(), 6 * n - 3);
    let nve_steps = (nve_ps * 1e-12 / dt).round() as usize;
    let mut o_nve: Vec<Vec<[f64; 3]>> = Vec::new();
    let mut nve_t = Vec::new();
    let mut nve_e = Vec::new();
    nve.prepare(&field, &at);
    for s in 0..nve_steps {
        nve.step(&field, &mut at, dt);
        if s % 50 == 0 {
            o_nve.push(waters.iter().map(|w| at[w[0]]).collect());
            nve_t.push(nve.temperature());
            nve_e.push(
                (nve.kinetic_energy() + nve.potential_energy().expect("stepped")) / KCAL_PER_MOL,
            );
        }
    }
    let rigid_nve = rigidity(&at, &waters);

    // Report.
    let e = Estimate::of(&energy);
    let t = Estimate::of(&temperature);
    let pr = Estimate::of(&press);
    println!(
        "NVT {nvt_ps} ps at 2 fs, γ = 1/ps: {:.2} ms a step; ⟨T⟩ after O {:.2} ± {:.2} K",
        1e3 * per_step,
        t.mean,
        t.error
    );
    println!("  rigidity after NVT {rigid_nvt:.2}, after NVE {rigid_nve:.2} (ε × the largest coordinate)");
    println!(
        "  U/N = {:.4} ± {:.4} kcal/mol (τ {:.1} samples); with the reciprocal bias put back {:.4}",
        e.mean,
        e.error,
        e.tau,
        e.mean + p.reciprocal_bias(q2) / KCAL_PER_MOL / n as f64
    );
    println!("    Jorgensen et al. 1983: −9.86, Monte Carlo NPT with a spherical cutoff (quoted, not read); Izadi et al. 2014: ΔH_vap 10.26, so about −9.67");
    println!(
        "  P = {:.0} ± {:.0} bar (τ {:.1})",
        pr.mean, pr.error, pr.tau
    );
    println!("    Yeh and Hummer 2004, PME at 33.00 nm⁻³ (0.987 g/cm³), N = 256: −2.9 ± 2.5 bar; here 33.33 nm⁻³");
    // g(r): the mean, its peak by a parabola through the three highest bins, and errors from four blocks.
    let peak = |g: &[f64]| -> (f64, f64) {
        let lo = (2.4 * ANGSTROM / width) as usize;
        let hi = (3.2 * ANGSTROM / width) as usize;
        let b = (lo..hi)
            .max_by(|&i, &j| g[i].total_cmp(&g[j]))
            .expect("bins");
        let (y0, y1, y2) = (g[b - 1], g[b], g[b + 1]);
        let shift = 0.5 * (y0 - y2) / (y0 - 2.0 * y1 + y2);
        let r = (b as f64 + 0.5 + shift) * width;
        (r, y1 - 0.25 * (y0 - y2) * shift)
    };
    let mean_g = |fr: &[Vec<f64>]| -> Vec<f64> {
        (0..bins)
            .map(|b| fr.iter().map(|g| g[b]).sum::<f64>() / fr.len() as f64)
            .collect()
    };
    let (r_peak, g_peak) = peak(&mean_g(&g_frames));
    let blocks: Vec<(f64, f64)> = g_frames
        .chunks(g_frames.len() / 4)
        .take(4)
        .map(|c| peak(&mean_g(c)))
        .collect();
    let spread = |xs: Vec<f64>| {
        let m = xs.iter().sum::<f64>() / xs.len() as f64;
        (xs.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (xs.len() * (xs.len() - 1)) as f64)
            .sqrt()
    };
    println!(
        "  g_OO first peak at {:.3} ± {:.3} Å, height {:.3} ± {:.3} (bins 0.02 Å, four blocks)",
        r_peak / ANGSTROM,
        spread(blocks.iter().map(|b| b.0).collect()) / ANGSTROM,
        g_peak,
        spread(blocks.iter().map(|b| b.1).collect())
    );
    println!("    Izadi, Anandakrishnan and Onufriev 2014, Table 3: first peak at 2.77 Å");
    let first_min = {
        let g = mean_g(&g_frames);
        let lo = (r_peak / width) as usize;
        let hi = (4.2 * ANGSTROM / width) as usize;
        let b = (lo..hi)
            .min_by(|&i, &j| g[i].total_cmp(&g[j]))
            .expect("bins");
        ((b as f64 + 0.5) * width, g[b])
    };
    println!(
        "  g_OO first minimum near {:.2} Å, {:.3}",
        first_min.0 / ANGSTROM,
        first_min.1
    );
    // Diffusion: blocks of 25 ps, each its own fit.
    let every = 50.0 * dt;
    let block = (25e-12 / every).round() as usize;
    let d_of = |frames: &[Vec<[f64; 3]>]| -> Estimate {
        let ds: Vec<f64> = frames
            .chunks(block)
            .filter(|c| c.len() == block)
            .map(|c| diffusion(c, every) * 1e4 / 1e-5)
            .collect();
        println!("    blocks: {ds:.3?}");
        Estimate::of(&ds)
    };
    // Yeh and Hummer's k_BT ξ/(6πηL) at each run's own temperature: the bath's for the run under
    // it and the NVE run's mean for that one. η is their 298 K value at both: its own temperature
    // dependence, which would raise the correction of a warmer run further, is not included.
    let correction_at =
        |t: f64| BOLTZMANN.to_si() * t * XI / (6.0 * PI * ETA_TIP3P * side) * 1e4 / 1e-5;
    let d_nvt = d_of(&o_nvt);
    let d_nve = d_of(&o_nve);
    let te = Estimate::of(&nve_t);
    let correction = correction_at(te.mean);
    let correction_nvt = correction_at(t.mean);
    let drift = (nve_e.last().expect("samples") - nve_e[0]) / nve_ps;
    println!(
        "NVE {nve_ps} ps: ⟨T⟩ {:.2} ± {:.2} K, total energy drift {:+.3} kcal/mol/ps over {} molecules",
        te.mean, te.error, drift, n
    );
    println!(
        "  D (10⁻⁵ cm²/s): NVE {:.3} ± {:.3}, under the 1/ps bath {:.3} ± {:.3}; Yeh–Hummer correction +{:.3} (η = 0.308 mPa s, L = {:.2} Å)",
        d_nve.mean,
        d_nve.error,
        d_nvt.mean,
        d_nvt.error,
        correction,
        side / ANGSTROM
    );
    println!(
        "  D∞ from NVE = {:.3} ± {:.3} at {:.1} K; under the bath {:.3} ± {:.3} at {:.1} K (correction +{:.3})",
        d_nve.mean + correction,
        d_nve.error,
        te.mean,
        d_nvt.mean + correction_nvt,
        d_nvt.error,
        t.mean,
        correction_nvt
    );
    println!("    Yeh and Hummer 2004 (PME, 33.00 nm⁻³): D_PBC 5.123(27) at N = 256, D₀ 6.05 by fit, 6.11 corrected");
    println!("    Mahoney and Jorgensen 2001 (267 waters, 9 Å cutoff, uncorrected): 5.19 ± 0.08 NPT, 5.06 ± 0.09 NVT");
    println!("    experiment 2.30 (as Mahoney and Jorgensen cite it)");
    assert!(rigid_nvt <= 16.0 && rigid_nve <= 16.0);
    assert!((t.mean - TEMPERATURE).abs() < 4.0 * t.error);
}

/// The liquid, 216 waters: 0.5 + 20 ps of melting, 100 ps NVT and 100 ps NVE.
#[test]
#[ignore = "liquid water in release, about fifteen minutes"]
fn liquid_tip3p_216() {
    liquid(6, 0x7193, 100.0, 100.0);
}

/// The same with 512 waters, a second box size for the finite-size correction: 50 ps NVT.
#[test]
#[ignore = "liquid water in release, about three quarters of an hour"]
fn liquid_tip3p_512() {
    liquid(8, 0x7194, 50.0, 100.0);
}

/// Milliseconds a step of rigid TIP3P at 2 fs, classical Ewald at δ = 10⁻⁵ and `r_c` = 9 Å, for
/// 216 and 1000 waters, and the parts: one force evaluation, and the constraints.
#[test]
#[ignore = "a timing, release only"]
fn the_cost_of_a_step() {
    for per_side in [6usize, 8, 10] {
        let mut wb = WaterBox::lattice(per_side, 3);
        let n = wb.count();
        let field = wb.force_field(9.0 * ANGSTROM, 1e-5);
        let mut md = wb.equilibrate(&field, TEMPERATURE, 0.2, 1);
        let mut at = wb.positions().to_vec();
        let mut best = f64::INFINITY;
        for _ in 0..3 {
            let s = std::time::Instant::now();
            md.run(&field, &mut at, 2.0 * FS, 20);
            best = best.min(s.elapsed().as_secs_f64() / 20.0);
        }
        let mut eval = f64::INFINITY;
        for _ in 0..3 {
            let s = std::time::Instant::now();
            std::hint::black_box(field.evaluate(&at));
            eval = eval.min(s.elapsed().as_secs_f64());
        }
        let settle = Settle::tip3p(wb.waters());
        let mut v = md.velocities().to_vec();
        let mut moved = at.clone();
        let s = std::time::Instant::now();
        for _ in 0..100 {
            for (p, q) in moved.iter_mut().zip(&v) {
                for k in 0..3 {
                    p[k] = std::hint::black_box(p[k] + 0.0 * q[k]);
                }
            }
            settle.constrain_positions(&at, &mut moved);
            settle.constrain_velocities(&at, &mut v);
        }
        let constraints = s.elapsed().as_secs_f64() / 100.0;
        println!(
            "{n:5} waters ({} atoms), {} wave vectors: {:.2} ms a step; one evaluation {:.2} ms; SETTLE and one projection {:.3} ms",
            3 * n,
            field.ewald().wave_vectors(),
            1e3 * best,
            1e3 * eval,
            1e3 * constraints
        );
    }
}
