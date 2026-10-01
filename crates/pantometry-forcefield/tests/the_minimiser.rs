//! **The minimiser finds minima that are known exactly, never raises the energy, and does not care
//! where the molecule is.**
//!
//! Every geometric tolerance here comes from the force the minimisation stopped at and a
//! stiffness: a displacement δx from a minimum leaves a force `H δx`, so `|δx| ≤ |F| / λ_min`, and
//! where a term's force has a direction of its own (a bond's along the bond, an angle's across it)
//! the bound is that term's alone. None of them is a number chosen to pass.

mod common;

use common::{
    angle, assert_minimised, distance, entry, hessian_extremes, index, len, max_force, positions,
    ANGSTROM, DEG,
};
use pantometry_core::Domain;
use pantometry_forcefield::energy::natural_length;
use pantometry_forcefield::minimise::{dihedral, DihedralRestraint, KCAL_PER_MOL_ANGSTROM};
use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::{Component, ForceField, Minimiser, Molecule, Status, UffType};

const AIN: &str = include_str!("../components/AIN.cif");
const MOH: &str = include_str!("../components/MOH.cif");
const EPS: f64 = f64::EPSILON;

fn force_field(c: &Component) -> ForceField {
    ForceField::new(c, &uff::assign(c)).expect("supported")
}

/// The tolerance the closed-form tests minimise to: 1e-6 kcal mol⁻¹ Å⁻¹, in newtons. That is at
/// the rounding floor for a molecule as stiff as methanol, so a stall is accepted when the
/// independently evaluated force is below that floor — see [`common::stall_floor`] — and the
/// bounds below use the force actually reached either way.
const TIGHT: f64 = 1e-6 * KCAL_PER_MOL_ANGSTROM;

/// [`TIGHT`] in kcal mol⁻¹ Å⁻¹.
const TIGHT_KCAL: f64 = 1e-6;

/// **Two carbons relax to eq 2's 1.514 Å.** A lone C–C single bond has only its stretch — the
/// pair is 1-2, so excluded from van der Waals — so the force on each atom is `k |r − r₀|` along
/// the bond, and the final force bounds the distance from 1.514 Å exactly: `|r − r₀| ≤ F / k`.
/// Plus the rounding of `r` itself, 4 ε r. **A check of where the minimum is, not of the
/// minimiser's convergence**: a one-dimensional quadratic is solved by the second step at any
/// tolerance, so this test cannot tell a minimiser that stops early from one that does not —
/// water and methane can.
#[test]
fn a_diatomic_relaxes_to_its_natural_length() {
    let c = entry(
        &[("C1", "C", [0.0, 0.0, 0.0]), ("C2", "C", [1.9, 0.3, -0.2])],
        &[("C1", "C2", "SING", false)],
    );
    let ff = force_field(&c);
    let mut at = positions(&c);
    let p = ff.minimise(&mut at, 1000, TIGHT);
    let k = ff.stretches()[0].force_constant * ANGSTROM * ANGSTROM / KCAL_PER_MOL;
    let f = max_force(&ff, &at);
    let r = distance(&at, 0, 1);
    println!(
        "diatomic: {:?} after {} steps, r = {r:.15} Å, F = {f:.2e}",
        p.status, p.steps
    );
    assert_eq!(p.status, Status::Converged);
    assert!((r - 1.514).abs() <= f / k + 4.0 * EPS * r, "{r}");
}

/// **Water relaxes to its natural O–H lengths and `O_3`'s θ₀ = 104.51°.** Nothing in it is
/// non-bonded (both H are 1-3), so the energy is two stretches and one bend. On each hydrogen the
/// stretch pushes along its bond and the bend pushes across it, at right angles, so with `F_H`
/// the force on that hydrogen: `|r − r₀| ≤ |F_H| / k`, and `|dE/dθ| = |F_H⊥| r ≤ |F_H| r`. Near
/// θ₀ the bend is `K Δθ (1 + O(Δθ))`, so `|Δθ| ≤ 1.01 |F_H| r / K` — the 1.01 is far more than
/// the `O(Δθ)` correction at Δθ ~ 1e-8 rad.
#[test]
fn water_relaxes_to_its_natural_lengths_and_angle() {
    let c = entry(
        &[
            ("O", "O", [0.0, 0.0, 0.0]),
            ("H1", "H", [0.9, 0.2, 0.1]),
            ("H2", "H", [-0.3, 1.1, -0.2]),
        ],
        &[("O", "H1", "SING", false), ("O", "H2", "SING", false)],
    );
    let ff = force_field(&c);
    assert!(ff.pairs().is_empty() && ff.torsions().is_empty() && ff.inversions().is_empty());
    let mut at = positions(&c);
    let p = ff.minimise(&mut at, 1000, TIGHT);
    assert_minimised(&ff, &at, &p, TIGHT_KCAL);
    let r0 = natural_length(UffType::O3, UffType::H, 1.0) / ANGSTROM;
    let k = ff.stretches()[0].force_constant * ANGSTROM * ANGSTROM / KCAL_PER_MOL;
    let kb = ff.bends()[0].force_constant / KCAL_PER_MOL;
    let forces = ff.evaluate(&at).forces;
    let to = ANGSTROM / KCAL_PER_MOL;
    let theta = angle(&at, 1, 0, 2);
    println!(
        "water: {} steps, r = {:.12} {:.12} Å (r₀ {r0:.12}), θ = {theta:.10}°",
        p.steps,
        distance(&at, 0, 1),
        distance(&at, 0, 2)
    );
    let mut theta_bound = 0.0f64;
    for h in [1, 2] {
        let fh = len(forces[h]) * to;
        let r = distance(&at, 0, h);
        assert!((r - r0).abs() <= fh / k + 4.0 * EPS * r, "H{h}: {r}");
        theta_bound = theta_bound.max(1.01 * fh * r / kb);
    }
    assert!(
        ((theta - 104.51) * DEG).abs() <= theta_bound + 1e-14,
        "θ = {theta}, bound {theta_bound:e} rad"
    );
}

/// **Methane relaxes to the regular tetrahedron with its natural C–H length.** `C_3`'s θ₀ is
/// 109.47°, not quite the tetrahedral 109.4712°, and six angles about one centre cannot all be
/// 109.47°; by symmetry the minimum is the regular tetrahedron, every angle `acos(−1/3)`, and
/// the angle terms push each H across its bond, so the bonds sit at r₀. Started from a distorted
/// methane. The bound: the final force `F` (all atoms, as one vector) is `H δx`, so
/// `|δx| ≤ |F| / λ_min`, λ_min the smallest Hessian eigenvalue with rigid motions removed,
/// measured at the end; then each bond moves by at most `√2 |δx|` and each angle by at most
/// `√6 |δx| / r_min` (the angle's gradient over its three atoms). The measured λ_min is a
/// central difference good to far better than the factor 0.9 it is taken with.
#[test]
fn methane_relaxes_to_the_regular_tetrahedron() {
    let c = entry(
        &[
            ("C", "C", [0.0, 0.0, 0.0]),
            ("H1", "H", [0.7, 0.6, 0.65]),
            ("H2", "H", [-0.7, -0.5, 0.7]),
            ("H3", "H", [-0.6, 0.65, -0.7]),
            ("H4", "H", [0.55, -0.7, -0.6]),
        ],
        &[
            ("C", "H1", "SING", false),
            ("C", "H2", "SING", false),
            ("C", "H3", "SING", false),
            ("C", "H4", "SING", false),
        ],
    );
    let ff = force_field(&c);
    let mut at = positions(&c);
    let p = ff.minimise(&mut at, 1000, TIGHT);
    assert_minimised(&ff, &at, &p, TIGHT_KCAL);
    let to = ANGSTROM / KCAL_PER_MOL;
    let f_all = ff
        .evaluate(&at)
        .forces
        .iter()
        .map(|f| len(*f) * len(*f))
        .sum::<f64>()
        .sqrt()
        * to;
    let (lambda, _) = hessian_extremes(&ff, &at);
    let dx = f_all / (0.9 * lambda);
    let r0 = natural_length(UffType::C3, UffType::H, 1.0) / ANGSTROM;
    let tetrahedral = (-1.0f64 / 3.0).acos() / DEG;
    let r_min = (1..5)
        .map(|h| distance(&at, 0, h))
        .fold(f64::INFINITY, f64::min);
    println!(
        "methane: {} steps, |F| {f_all:.2e}, λ_min {lambda:.3} kcal/mol/Å², |δx| ≤ {dx:.2e} Å",
        p.steps
    );
    for h in 1..5 {
        let r = distance(&at, 0, h);
        assert!(
            (r - r0).abs() <= 2f64.sqrt() * dx + 4.0 * EPS * r,
            "H{h}: {r}"
        );
    }
    for a in 1..5 {
        for b in a + 1..5 {
            let t = angle(&at, a, 0, b);
            assert!(
                ((t - tetrahedral) * DEG).abs() <= 6f64.sqrt() * dx / r_min + 1e-14,
                "H{a}-C-H{b}: {t}"
            );
        }
    }
}

/// **The energy never rises, step by step, and a step that is taken lowers it** — the Armijo
/// condition's guarantee, read off aspirin from the dictionary's ideal coordinates, which start
/// with a 134 kcal/mol clash. And when it stops, it says why: converged with the largest force at
/// most the tolerance (checked by evaluating the force independently), or stalled.
#[test]
fn the_energy_never_rises() {
    let c = Component::from_ccd(AIN).expect("AIN parses");
    let ff = force_field(&c);
    let mut at = positions(&c);
    let tol = 1e-4 * KCAL_PER_MOL_ANGSTROM;
    let mut m = Minimiser::new(tol);
    let mut last = ff.energy(&at).total;
    let mut moved = 0;
    for _ in 0..5000 {
        let before = at.clone();
        let p = m.step(&ff, &[], &mut at);
        assert!(
            p.energy <= last,
            "the energy rose: {} after {}",
            p.energy,
            last
        );
        if at != before {
            assert!(
                p.energy < last,
                "a step was taken without lowering the energy"
            );
            moved += 1;
        }
        last = p.energy;
        if p.status != Status::Running {
            break;
        }
    }
    let f = max_force(&ff, &at);
    println!(
        "aspirin: {:?} after {moved} steps, E {:.6} kcal/mol, max force {f:.2e} kcal/mol/Å",
        m.status(),
        last / KCAL_PER_MOL
    );
    assert_eq!(m.status(), Status::Converged);
    assert!(f <= 1e-4, "{f}");
}

/// **A converged minimiser does not move, and one given positions it did not leave starts its
/// history again rather than mixing two problems.**
#[test]
fn a_converged_minimiser_stays_put() {
    let c = Component::from_ccd(MOH).expect("MOH parses");
    let ff = force_field(&c);
    let mut at = positions(&c);
    let p = ff.minimise(&mut at, 2000, 1e-4 * KCAL_PER_MOL_ANGSTROM);
    assert_eq!(p.status, Status::Converged);
    let mut m = Minimiser::new(1e-4 * KCAL_PER_MOL_ANGSTROM);
    let held = at.clone();
    for _ in 0..3 {
        let q = m.step(&ff, &[], &mut at);
        assert_eq!(q.status, Status::Converged);
        assert_eq!(at, held);
    }
}

/// **The minimiser keeps its curvature history between steps.** Each step hands the positions
/// back as metres, and the history is kept only if the next step is given exactly those; a
/// comparison that went through a unit conversion would never match, so every step would start
/// again as steepest descent — still monotone, still converging, and slow. That happened while
/// this was written: aspirin took more than 6000 steps and stopped at 18.806 kcal/mol instead of
/// 18.573. After eight steps of aspirin from the dictionary, the history is full.
#[test]
fn the_curvature_history_survives_between_steps() {
    let c = Component::from_ccd(AIN).expect("AIN parses");
    let ff = force_field(&c);
    let mut at = positions(&c);
    let mut m = Minimiser::new(1e-4 * KCAL_PER_MOL_ANGSTROM);
    let mut most = 0;
    for _ in 0..12 {
        m.step(&ff, &[], &mut at);
        most = most.max(m.history_len());
    }
    assert!(most >= 8, "history reached only {most}");
}

/// The pairwise distances, Å.
fn distances(at: &[[f64; 3]]) -> Vec<f64> {
    let mut d = Vec::new();
    for i in 0..at.len() {
        for j in i + 1..at.len() {
            d.push(distance(at, i, j));
        }
    }
    d
}

/// **Moving and turning the start does not change the minimum.** Methanol from the dictionary,
/// and the same molecule shifted by (3.1, −1.7, 0.4) Å and turned by three Euler angles: both
/// minimised to 1e-6 kcal mol⁻¹ Å⁻¹. The paths differ by rounding, so the two minima are not
/// bit-equal; each is within `|F| / λ_min` of the true minimum, so every interatomic distance
/// agrees to `√2 (|δx_a| + |δx_b|)`, and the energies to `Σ |F|² / (2 λ_min)`, plus the
/// energy's rounding (1e-12 kcal/mol is a thousand ε of methanol's few kcal/mol).
#[test]
fn moving_or_turning_the_start_does_not_change_the_minimum() {
    let c = Component::from_ccd(MOH).expect("MOH parses");
    let ff = force_field(&c);
    let mut a = positions(&c);
    let (s1, c1) = 0.7f64.sin_cos();
    let (s2, c2) = (-1.1f64).sin_cos();
    let (s3, c3) = 2.3f64.sin_cos();
    let turn = |p: [f64; 3]| {
        let p = [c1 * p[0] - s1 * p[1], s1 * p[0] + c1 * p[1], p[2]];
        let p = [p[0], c2 * p[1] - s2 * p[2], s2 * p[1] + c2 * p[2]];
        let p = [c3 * p[0] - s3 * p[1], s3 * p[0] + c3 * p[1], p[2]];
        [
            p[0] + 3.1 * ANGSTROM,
            p[1] - 1.7 * ANGSTROM,
            p[2] + 0.4 * ANGSTROM,
        ]
    };
    let mut b: Vec<[f64; 3]> = a.iter().map(|p| turn(*p)).collect();
    let pa = ff.minimise(&mut a, 5000, TIGHT);
    let pb = ff.minimise(&mut b, 5000, TIGHT);
    assert_minimised(&ff, &a, &pa, TIGHT_KCAL);
    assert_minimised(&ff, &b, &pb, TIGHT_KCAL);
    let to = ANGSTROM / KCAL_PER_MOL;
    let total = |at: &[[f64; 3]]| {
        ff.evaluate(at)
            .forces
            .iter()
            .map(|f| len(*f) * len(*f))
            .sum::<f64>()
            .sqrt()
            * to
    };
    let (fa, fb) = (total(&a), total(&b));
    let (lambda, _) = hessian_extremes(&ff, &a);
    let lambda = 0.9 * lambda;
    let bound = 2f64.sqrt() * (fa + fb) / lambda;
    let worst = distances(&a)
        .iter()
        .zip(distances(&b))
        .map(|(x, y)| (x - y).abs())
        .fold(0.0f64, f64::max);
    let de = (pa.energy - pb.energy).abs() / KCAL_PER_MOL;
    println!(
        "methanol: λ_min {lambda:.3}, worst distance difference {worst:.2e} Å (bound {bound:.2e}), \
         energies differ by {de:.2e} kcal/mol"
    );
    assert!(worst <= bound + 1e-12, "{worst} > {bound}");
    assert!(de <= (fa * fa + fb * fb) / (2.0 * lambda) + 1e-12);
}

/// **Determinism: the same start gives the same minimum, bit for bit**, and a digest of it is
/// printed so that a debug and a release run can be compared — which is how "the same across
/// optimisation levels" was checked when this was written (the two digests matched; see the
/// changelog). Not pinned, because the force field's constructor calls the platform's `ln`,
/// `cos` and `asin`, and a pinned digest would assert those are the same function everywhere.
#[test]
fn the_same_start_gives_the_same_bits() {
    let c = Component::from_ccd(AIN).expect("AIN parses");
    let ff = force_field(&c);
    let run = || {
        let mut at = positions(&c);
        let p = ff.minimise(&mut at, 3000, 1e-4 * KCAL_PER_MOL_ANGSTROM);
        (at, p)
    };
    let (a, pa) = run();
    let (b, pb) = run();
    assert_eq!(a, b);
    assert_eq!(pa, pb);
    let mut digest: u64 = 0xcbf2_9ce4_8422_2325;
    for p in &a {
        for v in p {
            for byte in v.to_bits().to_le_bytes() {
                digest ^= u64::from(byte);
                digest = digest.wrapping_mul(0x0100_0000_01b3);
            }
        }
    }
    println!(
        "aspirin minimum after {} steps: FNV-1a digest {digest:016x}",
        pa.steps
    );
}

/// **The restraint's force is minus the gradient of its energy.** A central difference at h and
/// at 2h: their difference is three times the truncation of the one at h (the error is `c h²`
/// to leading order), so `|D(h) − D(2h)|` bounds it; plus the rounding of the two energies, each
/// good to 10 ε of the restraint energy plus that of the angle (`k · 20 ε · |Δφ|`), over h.
#[test]
fn the_restraints_force_is_its_gradient() {
    let c = Component::from_ccd(MOH).expect("MOH parses");
    let at = positions(&c);
    let atoms = [
        index(&c, "HO"),
        index(&c, "O"),
        index(&c, "C"),
        index(&c, "H1"),
    ];
    let r = DihedralRestraint {
        atoms,
        target: dihedral(&at, atoms) + 0.4,
        stiffness: 100.0 * KCAL_PER_MOL,
    };
    let with = vec![[0.0; 3]; at.len()];
    let mut forces = with.clone();
    let e0 = r.add_forces(&at, &mut forces);
    assert!((e0 - r.energy_at(&at)).abs() <= 4.0 * EPS * e0);
    let to = ANGSTROM / KCAL_PER_MOL;
    let k = r.stiffness / KCAL_PER_MOL;
    let round = 10.0 * EPS * e0 / KCAL_PER_MOL + k * 20.0 * EPS * 0.4;
    for &a in &atoms {
        for axis in 0..3 {
            let d = |h: f64| {
                let mut p = at.clone();
                let mut q = at.clone();
                p[a][axis] += h * ANGSTROM;
                q[a][axis] -= h * ANGSTROM;
                (r.energy_at(&p) - r.energy_at(&q)) / KCAL_PER_MOL / (2.0 * h)
            };
            let (d1, d2) = (d(1e-5), d(2e-5));
            let analytic = -(forces[a][axis] - with[a][axis]) * to;
            let tol = (d1 - d2).abs() + 2.0 * round / 1e-5;
            assert!(
                (analytic - d1).abs() <= tol,
                "atom {a} axis {axis}: {analytic} vs {d1}, tol {tol:e}"
            );
        }
    }
}

/// **A molecule relaxes one minimiser step per kernel step, and says so in its readings**: the
/// energy reading never rises, the positions change until it converges, `converged` turns to 1,
/// and `max force` is then at most the tolerance. A checkpoint taken mid-run and restored repeats
/// the next step bit for bit.
#[test]
fn a_molecule_relaxes_one_step_at_a_time() {
    use pantometry_core::{Exchange, Schedule, Simulation};
    use pantometry_units::Time;
    let c = Component::from_ccd(MOH).expect("MOH parses");
    let mut m = Molecule::new("methanol", c);
    let reading = |m: &Molecule, label: &str| {
        m.readings()
            .into_iter()
            .find(|r| r.label == label)
            .unwrap_or_else(|| panic!("no {label}"))
            .value
    };
    let start: Vec<[f64; 3]> = (0..6).map(|i| m.at(i)).collect();
    let mut bus = Exchange::default();
    m.step(Time::s(0.0), Time::s(1.0), &mut bus).expect("steps");
    let after: Vec<[f64; 3]> = (0..6).map(|i| m.at(i)).collect();
    assert_ne!(start, after, "the first step moves the atoms");

    m.checkpoint();
    m.step(Time::s(0.0), Time::s(1.0), &mut bus).expect("steps");
    let once: Vec<[f64; 3]> = (0..6).map(|i| m.at(i)).collect();
    m.restore();
    m.step(Time::s(0.0), Time::s(1.0), &mut bus).expect("steps");
    let again: Vec<[f64; 3]> = (0..6).map(|i| m.at(i)).collect();
    assert_eq!(once, again, "a restored minimiser repeats its step");

    let mut sim = Simulation::new(Schedule::OneWay).with(m);
    let mut last = f64::INFINITY;
    for _ in 0..2000 {
        sim.advance(Time::s(1.0)).expect("runs");
        let m = sim.domain_as::<Molecule>("methanol").expect("itself");
        let e = reading(m, "energy");
        assert!(e <= last, "{e} after {last}");
        last = e;
        if reading(m, "converged") == 1.0 {
            break;
        }
    }
    let m = sim.domain_as::<Molecule>("methanol").expect("itself");
    assert_eq!(reading(m, "converged"), 1.0);
    assert!(reading(m, "max force") <= 1e-4);
    assert!(reading(m, "rms force") <= reading(m, "max force"));
    println!(
        "methanol in a simulation: converged after {} steps, E {:.6} kcal/mol",
        reading(m, "minimiser steps"),
        last
    );
}

/// **The two-loop recursion is the BFGS update, written out densely.** For three curvature pairs
/// in six dimensions, `two_loop` must equal `−H g` with `H` built by applying
/// `H ← (I − ρ s yᵀ) H (I − ρ y sᵀ) + ρ s sᵀ`, `ρ = 1/(y·s)`, oldest pair first, to
/// `H₀ = (s·y / y·y) I` from the newest pair — Nocedal's definition, which is a different
/// computation from the recursion and shares none of its loops. A wrong pairing of the second
/// loop's α with its pair, a different `H₀`, or a direction that ignores the pairs all give a
/// different vector. The pairs are `y = A s` for a fixed symmetric matrix `A` whose eigenvalues lie
/// between 1 and 3 (diagonal 2, off-diagonals of 0.15), so `s·y > 0` and nothing is
/// ill-conditioned; each side is a few hundred roundings on quantities of order one, so they
/// agree to 1e-12 of the result's size.
#[test]
fn the_two_loop_recursion_is_the_bfgs_update() {
    use pantometry_forcefield::minimise::two_loop;
    let n = 6;
    let a: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    if i == j {
                        2.0
                    } else {
                        0.15 * (((i + j) % 3) as f64 - 1.0)
                    }
                })
                .collect()
        })
        .collect();
    let mul = |m: &Vec<Vec<f64>>, v: &[f64]| -> Vec<f64> {
        m.iter()
            .map(|row| row.iter().zip(v).map(|(x, y)| x * y).sum())
            .collect()
    };
    let dot = |u: &[f64], v: &[f64]| -> f64 { u.iter().zip(v).map(|(x, y)| x * y).sum() };
    let pairs: Vec<(Vec<f64>, Vec<f64>)> = (0..3)
        .map(|k| {
            let s: Vec<f64> = (0..n)
                .map(|i| (1.3 * (i + 7 * k) as f64 + 0.4).sin())
                .collect();
            let y = mul(&a, &s);
            (s, y)
        })
        .collect();
    let g: Vec<f64> = (0..n).map(|i| (0.7 * i as f64 - 1.1).cos()).collect();

    let (s_new, y_new) = pairs.last().expect("pairs");
    let gamma = dot(s_new, y_new) / dot(y_new, y_new);
    let mut h: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { gamma } else { 0.0 }).collect())
        .collect();
    for (s, y) in &pairs {
        let rho = 1.0 / dot(y, s);
        // V = I − ρ y sᵀ; H ← Vᵀ H V + ρ s sᵀ.
        let v: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| (if i == j { 1.0 } else { 0.0 }) - rho * y[i] * s[j])
                    .collect()
            })
            .collect();
        let mut next = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                let mut t = 0.0;
                for k in 0..n {
                    for l in 0..n {
                        t += v[k][i] * h[k][l] * v[l][j];
                    }
                }
                next[i][j] = t + rho * s[i] * s[j];
            }
        }
        h = next;
    }
    let dense: Vec<f64> = mul(&h, &g).iter().map(|x| -x).collect();
    let recursion = two_loop(&pairs, &g);
    let size = dot(&dense, &dense).sqrt();
    let worst = dense
        .iter()
        .zip(&recursion)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0f64, f64::max);
    println!("two-loop against dense BFGS: worst difference {worst:.2e} of {size:.3}");
    assert!(worst <= 1e-12 * size, "{worst:e}");
}

/// **The history buys steps: L-BFGS reaches aspirin's minimum in fewer steps than the same
/// minimiser without memory**, which is steepest descent with the same Armijo search. If the
/// history were built and never used, the two runs would be the same arithmetic and take exactly
/// the same number of steps, so `fewer` needs no threshold to mean something. Both to
/// 1e-4 kcal mol⁻¹ Å⁻¹; steepest descent is given up to 20000 steps and counted as 20000 if it
/// has not converged by then.
#[test]
fn the_history_buys_steps() {
    let c = Component::from_ccd(AIN).expect("AIN parses");
    let ff = force_field(&c);
    let tol = 1e-4 * KCAL_PER_MOL_ANGSTROM;
    let run = |memory: usize, budget: usize| {
        let mut at = positions(&c);
        let mut m = Minimiser::with_memory(tol, memory);
        let mut p = m.step(&ff, &[], &mut at);
        while p.status == Status::Running && p.steps < budget {
            p = m.step(&ff, &[], &mut at);
        }
        p
    };
    let lbfgs = run(pantometry_forcefield::minimise::MEMORY, 20_000);
    let steepest = run(0, 20_000);
    println!(
        "aspirin to 1e-4: L-BFGS {:?} in {} steps at {:.6} kcal/mol; steepest descent {:?} in {} \
         steps at {:.6}",
        lbfgs.status,
        lbfgs.steps,
        lbfgs.energy / KCAL_PER_MOL,
        steepest.status,
        steepest.steps,
        steepest.energy / KCAL_PER_MOL
    );
    assert_eq!(lbfgs.status, Status::Converged);
    assert!(lbfgs.steps < steepest.steps);
}
