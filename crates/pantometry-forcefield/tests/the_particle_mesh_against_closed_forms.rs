//! **Smooth particle-mesh Ewald against what it must equal**: the B-spline against its
//! truncated-power closed form, the Euler factor `|b(m)|²` against its cotangent series and the
//! tangent numbers, the mesh's energy against the interpolation's own Fourier series summed with
//! no grid, the convergence rates the spline order predicts, crystals and a lone charge on the
//! grid exactly, the derived bias and spread against measurement, the forces and the virial as
//! the derivatives of this energy, the net force as its translation derivative, a decoupling's
//! cross terms, and the bits.
//!
//! Point charges, with exclusions by hand, except for the decoupling and the water box, which are
//! TIP3P. The classical sum converged — `k_c` at twelve times α, where `exp(−k²/4α²)` is `e⁻³⁶` —
//! is the reference only where the claim is convergence: a different discretisation of the same
//! lattice sum, converged to a stated precision.

// The sums over atoms and axes are written as the formulas' index loops, which read against
// them; iterator chains over two arrays at once would not.
#![allow(clippy::needless_range_loop)]

use pantometry_forcefield::ewald::{erfc, COULOMB};
use pantometry_forcefield::pme::{b_spline, euler_denominator, MAX_ORDER, MIN_ORDER};
use pantometry_forcefield::{
    Alchemical, Ewald, EwaldEnergy, EwaldParameters, Lambda, PeriodicBox, PeriodicDecoupling,
    PeriodicEnergy, PeriodicForceField, Pme, PmeParameters, WaterBox,
};
use std::f64::consts::PI;

const ANGSTROM: f64 = 1e-10;
const EPS: f64 = f64::EPSILON;

fn splitmix(seed: u64, i: u64) -> u64 {
    let mut z = seed.wrapping_add((i + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Uniform in [0, 1), from a seed and an index.
fn uniform(seed: u64, i: u64) -> f64 {
    (splitmix(seed, i) >> 11) as f64 / (1u64 << 53) as f64
}

/// A disordered neutral box of `n` charges of ±0.8 e at least 1.2 Å apart, edges `lengths`.
fn disordered_in(n: usize, lengths: [f64; 3], seed: u64) -> (Vec<f64>, Vec<[f64; 3]>, PeriodicBox) {
    let cell = PeriodicBox::new(lengths);
    let mut at: Vec<[f64; 3]> = Vec::new();
    let mut i = 0u64;
    while at.len() < n {
        let p = [0, 1, 2].map(|a| lengths[a as usize] * uniform(seed, 3 * i + a));
        i += 1;
        let close = at.iter().any(|o| {
            let d = cell.minimum_image([p[0] - o[0], p[1] - o[1], p[2] - o[2]]);
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() < 1.2 * ANGSTROM
        });
        if !close {
            at.push(p);
        }
    }
    let q = (0..n)
        .map(|k| if k % 2 == 0 { 0.8 } else { -0.8 })
        .collect();
    (q, at, cell)
}

/// The same in a cube.
fn disordered(n: usize, side: f64, seed: u64) -> (Vec<f64>, Vec<[f64; 3]>, PeriodicBox) {
    disordered_in(n, [side; 3], seed)
}

/// `α`, a real-space cutoff, and `k_c` at twelve times α: the classical sum converged to `e⁻³⁶`
/// in its reciprocal part, the reference of a convergence.
fn converged(alpha: f64, cutoff: f64) -> EwaldParameters {
    EwaldParameters {
        alpha,
        cutoff,
        k_cutoff: 12.0 * alpha,
    }
}

fn mesh(order: usize, grid: [usize; 3]) -> PmeParameters {
    PmeParameters { order, grid }
}

/// `x^k` by multiplication: never `powi`, whose rounding Rust leaves to the platform.
fn power(x: f64, k: usize) -> f64 {
    let mut v = 1.0;
    for _ in 0..k {
        v *= x;
    }
    v
}

fn binomial(n: usize, k: usize) -> f64 {
    let mut b = 1.0;
    for i in 0..k {
        b = b * (n - i) as f64 / (i + 1) as f64;
    }
    b
}

fn factorial(n: usize) -> f64 {
    (2..=n).fold(1.0, |f, i| f * i as f64)
}

/// The cardinal B-spline by its closed form, `M_n(u) = Σ_k (−1)^k C(n, k) (u − k)₊^{n−1} / (n−1)!`,
/// and its derivative, the same sum with `n − 2`: `(value, derivative, the sum of the terms'
/// magnitudes)`, the last the scale its cancellation is measured against.
fn truncated_power(n: usize, u: f64) -> (f64, f64, f64) {
    let (mut v, mut d, mut scale) = (0.0, 0.0, 0.0);
    for k in 0..=n {
        let x = u - k as f64;
        if x <= 0.0 {
            continue;
        }
        let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
        let c = binomial(n, k);
        let t = c * power(x, n - 1) / factorial(n - 1);
        let s = c * power(x, n - 2) / factorial(n - 2);
        v += sign * t;
        d += sign * s;
        scale += t + s;
    }
    (v, d, scale)
}

/// The B-spline the mesh spreads with is the closed form, at every order: each of the `n`
/// weights `M_n(t + j)` and derivatives `M_n′(t + j)` against the truncated-power sum, which is
/// written from the definition and shares nothing with the recursion; the weights sum to one and
/// the derivatives to zero — the partition of unity and its derivative — and the derivative is
/// `M_{n−1}(u) − M_{n−1}(u − 1)` from the order below, to the bit, since that is how it is built.
#[test]
fn the_b_spline_is_its_truncated_power_sum() {
    let mut worst = 0.0f64;
    for n in MIN_ORDER..=MAX_ORDER {
        for t in [0.0, 1e-9, 0.125, 0.3, 0.5, 0.61803, 0.875, 0.999_999] {
            let (w, dw) = b_spline(n, t);
            assert_eq!(w.len(), n);
            let (mut sum, mut slope, mut magnitude) = (0.0, 0.0, 0.0);
            for j in 0..n {
                let (v, d, scale) = truncated_power(n, t + j as f64);
                // The closed form cancels: its own rounding is a few ε of its terms' magnitudes.
                let tolerance = 8.0 * EPS * scale + 4.0 * EPS;
                assert!(
                    (w[j] - v).abs() <= tolerance,
                    "M_{n}({}) = {} against {v}",
                    t + j as f64,
                    w[j]
                );
                assert!(
                    (dw[j] - d).abs() <= tolerance,
                    "M_{n}′({}) = {} against {d}",
                    t + j as f64,
                    dw[j]
                );
                worst = worst.max((w[j] - v).abs().max((dw[j] - d).abs()) / tolerance);
                sum += w[j];
                slope += dw[j];
                magnitude += dw[j].abs();
            }
            if n > MIN_ORDER {
                // M_n′(u) = M_{n−1}(u) − M_{n−1}(u − 1), M_{n−1} zero off its n − 1 points.
                let (below, _) = b_spline(n - 1, t);
                for j in 0..n {
                    let upper = if j < n - 1 { below[j] } else { 0.0 };
                    let lower = if j > 0 { below[j - 1] } else { 0.0 };
                    assert_eq!(dw[j].to_bits(), (upper - lower).to_bits(), "n {n} j {j}");
                }
            }
            assert!(
                (sum - 1.0).abs() <= 2.0 * n as f64 * EPS,
                "n {n} t {t}: Σ M = {sum}"
            );
            assert!(
                slope.abs() <= 2.0 * EPS * magnitude,
                "n {n} t {t}: Σ M′ = {slope}"
            );
        }
    }
    println!("B-splines: worst |M − closed form| is {worst:.3} of its tolerance");
}

/// `|Σ_{s=1}^{n−1} M_n(s) e^{iθs}|`, `θ = 2πm/K`, in closed form: Poisson's summation makes it
/// `sinⁿ(θ/2) |Σ_p (θ/2 − πp)⁻ⁿ|`, and the series of `cot`'s derivatives sum that: with
/// `s = sin(θ/2)` and `c = cos(θ/2)`, `c` for `n = 3`, `1 − 2s²/3` for 4, `c(1 − s²/3)` for 5 and
/// `1 − s² + 2s⁴/15` for 6.
fn cotangent_series(n: usize, theta: f64) -> f64 {
    let (s, c) = (0.5 * theta).sin_cos();
    let s2 = s * s;
    match n {
        3 => c.abs(),
        4 => 1.0 - 2.0 * s2 / 3.0,
        5 => (c * (1.0 - s2 / 3.0)).abs(),
        6 => 1.0 - s2 + 2.0 * s2 * s2 / 15.0,
        _ => unreachable!(),
    }
}

/// The Euler exponential spline's factor against its closed form, and its zero.
///
/// - `|b(m)|²` as the mesh holds it, for orders 3 to 6 on grids of 8, 16 and 64 at every `m` but
///   `K/2`, is one over the square of [`cotangent_series`] to `1e-13` of itself.
/// - **At `m = K/2`** the denominator is `|Σ_s M_n(s)(−1)^s|²`: for odd `n` exactly zero, since
///   `M_n(s) = M_n(n − s)` and `(−1)^s = −(−1)^{n−s}`; for even `n` the square of the coefficient
///   of `x^{n−1}` in `tan x` — `1/3`, `2/15`, `17/315`, `62/2835`, `1382/155925` for 4 to 12. The
///   mesh leaves that plane out at every order, so its `|b|²` there is zero, not infinite.
#[test]
fn the_euler_factor_is_its_cotangent_series_and_its_zero_is_left_out() {
    let mut worst = 0.0f64;
    let cell = PeriodicBox::cubic(10.0 * ANGSTROM);
    for n in 3..=6 {
        for k in [8usize, 16, 64] {
            let pme = Pme::new(cell, 0.3 / ANGSTROM, mesh(n, [k; 3]));
            for j in 0..k {
                let b2 = pme.b_squared(0, j);
                if 2 * j == k {
                    assert_eq!(b2, 0.0, "the Nyquist plane is left out");
                    continue;
                }
                let d = cotangent_series(n, 2.0 * PI * j as f64 / k as f64);
                let want = 1.0 / (d * d);
                assert!(
                    (b2 / want - 1.0).abs() <= 1e-13,
                    "n {n}, K {k}, m {j}: |b|² {b2} against {want}"
                );
                worst = worst.max((b2 / want - 1.0).abs());
                for axis in 1..3 {
                    assert_eq!(pme.b_squared(axis, j).to_bits(), b2.to_bits());
                }
            }
        }
    }
    let tangent = [
        (4usize, 1.0 / 3.0),
        (6, 2.0 / 15.0),
        (8, 17.0 / 315.0),
        (10, 62.0 / 2835.0),
        (12, 1382.0 / 155_925.0),
    ];
    for (n, t) in tangent {
        for k in [16usize, 32] {
            let d = euler_denominator(n, k / 2, k);
            assert!(
                (d / (t * t) - 1.0).abs() <= 1e-13,
                "n {n}: {d} against {}",
                t * t
            );
        }
    }
    for n in [3usize, 5, 7, 9, 11] {
        let d = euler_denominator(n, 8, 16);
        println!("odd order {n}: the denominator at K/2 is {d:e}");
        assert!(d <= 1e-30, "n {n}: {d:e}");
        // And a neighbour is not: the zero is the Nyquist plane's alone.
        assert!(euler_denominator(n, 7, 16) > 1e-6);
    }
    println!("|b(m)|²: worst relative difference from the cotangent series {worst:.2e}");
}

/// The interpolated structure factor of one charge on one axis: `b(m) Σ_k M_n(u − k) e^{2πimk/K}`
/// is exactly `e^{2πimu/K} (1 + Σ_{p≠0} r_p e^{−2πipu}) / (1 + Σ_{p≠0} r_p)`, `r_p = (x/(x − p))ⁿ`,
/// `x = m/K`, by Poisson's summation — the module documentation's derivation. Here to `|p| =
/// 1000`, where the omitted terms are below `x^n 1000^{1−n}`. `frac` is the fractional
/// coordinate in `[0, 1)` and `u = K frac`.
fn interpolated(n: usize, m: f64, k: usize, frac: f64) -> (f64, f64) {
    let x = m / k as f64;
    let u = k as f64 * frac;
    let (s1, c1) = (2.0 * PI * u).sin_cos();
    // e^{−2πipu} for p = 1, 2, … by rotation, and its conjugate for −p.
    let (mut er, mut ei) = (1.0, 0.0);
    let (mut sr, mut si, mut r) = (1.0, 0.0, 0.0);
    for p in 1..=1000 {
        let (nr, ni) = (er * c1 + ei * s1, ei * c1 - er * s1);
        er = nr;
        ei = ni;
        let plus = power(x / (x - p as f64), n);
        let minus = power(x / (x + p as f64), n);
        sr += plus * er + minus * er;
        si += plus * ei - minus * ei;
        r += plus + minus;
    }
    let (sp, cp) = (2.0 * PI * m * frac).sin_cos();
    let (fr, fi) = (sr / (1.0 + r), si / (1.0 + r));
    (cp * fr - sp * fi, cp * fi + sp * fr)
}

/// The signed vector a grid index stands for, `None` on the Nyquist plane.
fn signed(j: usize, k: usize) -> Option<f64> {
    if 2 * j == k {
        None
    } else if 2 * j < k {
        Some(j as f64)
    } else {
        Some(j as f64 - k as f64)
    }
}

/// The mesh's reciprocal energy written without the mesh: `(2π/V) Σ_k exp(−k²/4α²)/k² |S̃(m)|²`
/// over the grid's vectors, `m_a` in `(−K_a/2, K_a/2)`, with `S̃(m) = Σ_j q_j Π_a` [`interpolated`]
/// — no spreading, no transform, no `b(m)` table. With `shear = Some((a, b, ε))` every wave
/// vector is the one of the box sheared by `x_a → x_a + ε x_b`, `k_b → k_b − ε k_a`, at fixed
/// fractional coordinates: the strain a virial differentiates. Joules.
fn series_energy(
    q: &[f64],
    at: &[[f64; 3]],
    cell: &PeriodicBox,
    alpha: f64,
    p: PmeParameters,
    shear: Option<(usize, usize, f64)>,
) -> f64 {
    let l = cell.lengths();
    let grid = p.grid;
    // Per atom, per axis, per grid index: the interpolated phase.
    let factors: Vec<[Vec<(f64, f64)>; 3]> = at
        .iter()
        .map(|r| {
            [0, 1, 2].map(|a| {
                let f = r[a] / l[a];
                let frac = f - f.floor();
                (0..grid[a])
                    .map(|j| {
                        signed(j, grid[a])
                            .map_or((0.0, 0.0), |m| interpolated(p.order, m, grid[a], frac))
                    })
                    .collect()
            })
        })
        .collect();
    let mut sum = 0.0;
    for j1 in 0..grid[0] {
        let Some(m1) = signed(j1, grid[0]) else {
            continue;
        };
        for j2 in 0..grid[1] {
            let Some(m2) = signed(j2, grid[1]) else {
                continue;
            };
            for j3 in 0..grid[2] {
                let Some(m3) = signed(j3, grid[2]) else {
                    continue;
                };
                if j1 == 0 && j2 == 0 && j3 == 0 {
                    continue;
                }
                let k0 = [0, 1, 2].map(|a| 2.0 * PI * [m1, m2, m3][a] / l[a]);
                let mut k = k0;
                if let Some((a, b, eps)) = shear {
                    k[b] -= eps * k0[a];
                }
                let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                let (mut sr, mut si) = (0.0, 0.0);
                for (i, f) in factors.iter().enumerate() {
                    let (ar, ai) = f[0][j1];
                    let (br, bi) = f[1][j2];
                    let (cr, ci) = f[2][j3];
                    let (xr, xi) = (ar * br - ai * bi, ar * bi + ai * br);
                    sr += q[i] * (xr * cr - xi * ci);
                    si += q[i] * (xr * ci + xi * cr);
                }
                sum += (-k2 / (4.0 * alpha * alpha)).exp() / k2 * (sr * sr + si * si);
            }
        }
    }
    COULOMB * 2.0 * PI / cell.volume() * sum
}

/// **The mesh is the interpolation it claims to be**: its reciprocal energy is [`series_energy`],
/// summed from the interpolation's Fourier series with no grid, no transform and no `b(m)` table,
/// to `1e-12` of itself — orders 5, 6 and 8, in a box of 13 × 9.5 × 17 Å on a grid of
/// 16 × 8 × 32, so that every axis has its own length and its own grid, with a net charge.
#[test]
fn the_mesh_is_its_interpolations_fourier_series() {
    let (mut q, at, cell) = disordered_in(24, [13.0, 9.5, 17.0].map(|x| x * ANGSTROM), 0x5E71);
    q[3] += 0.3;
    let alpha = 0.3 / ANGSTROM;
    let p = EwaldParameters {
        alpha,
        cutoff: 4.0 * ANGSTROM,
        k_cutoff: alpha,
    };
    for order in [5usize, 6, 8] {
        let grid = mesh(order, [16, 8, 32]);
        let e = Ewald::new(cell, p)
            .with_mesh(grid)
            .evaluate(&q, &at, &[])
            .energy;
        let s = series_energy(&q, &at, &cell, alpha, grid, None);
        println!(
            "order {order}: the mesh {:.15e} J, the series {s:.15e}, off {:.1e}",
            e.reciprocal,
            e.reciprocal / s - 1.0
        );
        assert!((e.reciprocal / s - 1.0).abs() <= 1e-12, "order {order}");
    }
}

/// Crystals: positions of `basis` (fractions of the cubic cell `a`, with charges) repeated `m[k]`
/// times along axis `k`, every position moved by `shift`, and the box.
fn crystal(
    basis: &[([f64; 3], f64)],
    a: f64,
    m: [usize; 3],
    shift: [f64; 3],
) -> (Vec<f64>, Vec<[f64; 3]>, PeriodicBox) {
    let mut q = Vec::new();
    let mut at = Vec::new();
    for i in 0..m[0] {
        for j in 0..m[1] {
            for k in 0..m[2] {
                for &(f, c) in basis {
                    q.push(c);
                    let g = [i, j, k];
                    at.push([0, 1, 2].map(|x| (g[x] as f64 + f[x]) * a + shift[x]));
                }
            }
        }
    }
    (q, at, PeriodicBox::new(m.map(|n| n as f64 * a)))
}

fn rock_salt() -> Vec<([f64; 3], f64)> {
    let fcc = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];
    let mut b: Vec<([f64; 3], f64)> = fcc.iter().map(|&f| (f, 1.0)).collect();
    b.extend(fcc.iter().map(|&f| ([f[0] + 0.5, f[1], f[2]], -1.0)));
    b
}

fn caesium_chloride() -> Vec<([f64; 3], f64)> {
    vec![([0.0, 0.0, 0.0], 1.0), ([0.5, 0.5, 0.5], -1.0)]
}

/// Every real-space term the Ewald sum with `p` leaves out, by brute force, joules: `q_i q_j
/// erfc(αr)/r` over every image pair at `r_c ≤ r < r_c + 7/α`, half the double sum.
fn omitted_real(q: &[f64], at: &[[f64; 3]], cell: &PeriodicBox, p: &EwaldParameters) -> f64 {
    let l = cell.lengths();
    let reach = p.cutoff + 7.0 / p.alpha;
    let images = l.map(|x| (reach / x).ceil() as i32 + 1);
    let mut real = 0.0;
    for i in 0..q.len() {
        for j in 0..q.len() {
            for nx in -images[0]..=images[0] {
                for ny in -images[1]..=images[1] {
                    for nz in -images[2]..=images[2] {
                        let n = [nx, ny, nz];
                        let d = [0, 1, 2].map(|a| at[i][a] - at[j][a] + f64::from(n[a]) * l[a]);
                        let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                        if r >= p.cutoff && r < reach {
                            real += 0.5 * q[i] * q[j] * erfc(p.alpha * r) / r;
                        }
                    }
                }
            }
        }
    }
    COULOMB * real
}

/// **Madelung's constants through the mesh, exactly, with every ion on a grid point**: the Euler
/// spline is exact at integer `u`, so a crystal whose ions sit on grid points has the exact
/// structure factor at every vector of the grid, and its energy is the lattice sum's but for the
/// vectors beyond the grid, `e⁻¹⁶⁰` here, and the real-space pairs beyond `r_c`, summed back by
/// brute force. NaCl, 64 ions, `1.747564594633182`, and CsCl, 16, `1.762674773070988`, on the
/// nearest-neighbour distance: to `1e-12` at orders 3, 4, 6 and 8 — whatever the order, which is
/// the point. Moved off the grid, the same crystals are not exact, and their error falls with the
/// grid at the rate the convergence test holds.
#[test]
fn madelungs_constants_on_the_grid_are_exact_at_every_order() {
    let r0 = 2.8 * ANGSTROM;
    let alpha = 0.35 / ANGSTROM;
    let cases = [
        (
            "NaCl",
            rock_salt(),
            2.0 * r0,
            [2usize; 3],
            1.747_564_594_633_182,
            4.0,
        ),
        (
            "CsCl",
            caesium_chloride(),
            2.0 * r0 / 3f64.sqrt(),
            [2; 3],
            1.762_674_773_070_988,
            1.0,
        ),
    ];
    for (name, basis, a, m, exact, per_cell) in cases {
        let pairs = per_cell * (m[0] * m[1] * m[2]) as f64;
        for (shifted, shift) in [
            (false, [0.0; 3]),
            (true, [0.37, 0.11, 0.23].map(|x| x * ANGSTROM)),
        ] {
            let (q, at, cell) = crystal(&basis, a, m, shift);
            let p = EwaldParameters {
                alpha,
                cutoff: cell.half_shortest_edge(),
                k_cutoff: alpha,
            };
            let back = omitted_real(&q, &at, &cell, &p);
            for order in [3usize, 4, 6, 8] {
                let e = Ewald::new(cell, p)
                    .with_mesh(mesh(order, [32; 3]))
                    .evaluate(&q, &at, &[])
                    .energy;
                let madelung = -(e.total + back) * r0 / (pairs * COULOMB);
                println!(
                    "{name}{}, order {order}: M = {madelung:.15}, off {:.1e}",
                    if shifted { " off the grid" } else { "" },
                    madelung - exact
                );
                if shifted {
                    // Off the grid the spline is not exact: the error is far above the 6.7e-16
                    // on it, 4.9e-12 at order 8 and more below.
                    assert!((madelung - exact).abs() > 1e-12, "the shift is seen");
                } else {
                    assert!(
                        (madelung - exact).abs() <= 1e-12,
                        "{name} order {order}: {madelung} against {exact}"
                    );
                }
            }
        }
    }
}

/// The lattice sum of a lone unit charge over the vectors a grid leaves out — those with some
/// `|m_a| ≥ K_a/2`, the Nyquist planes among them — out to `3K`, where `exp(−k²/4α²)` is below
/// `e⁻²⁰⁰` here: `k_e (2π/V) Σ exp(−k²/4α²)/k²`, joules.
fn omitted_vectors(cell: &PeriodicBox, alpha: f64, grid: [usize; 3]) -> f64 {
    let l = cell.lengths();
    let reach = grid.map(|k| 3 * k as i64);
    let mut sum = 0.0;
    for m1 in -reach[0]..=reach[0] {
        for m2 in -reach[1]..=reach[1] {
            for m3 in -reach[2]..=reach[2] {
                let m = [m1, m2, m3];
                if (0..3).all(|a| 2 * m[a].unsigned_abs() < grid[a] as u64) {
                    continue;
                }
                let k = [0, 1, 2].map(|a| 2.0 * PI * m[a] as f64 / l[a]);
                let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                sum += (-k2 / (4.0 * alpha * alpha)).exp() / k2;
            }
        }
    }
    COULOMB * 2.0 * PI / cell.volume() * sum
}

/// **A lone charge in a cube is the Wigner lattice, `ξ = −2.837297479480620`, exactly on a grid
/// point**, for the same reason, at every order: on 32³ to `1e-12`, where the vectors beyond the
/// grid weigh `e⁻⁴⁹`; and on 16³, where they do not, it is the lattice less exactly those vectors,
/// summed here, Nyquist planes included. Off the grid the error is the mesh's self term alone,
/// whose mean over positions is [`Pme::self_bias`] in closed form — checked in a box of
/// 15 × 12 × 18 Å, so that each axis's length is its own. Averaged over 4 × 4 × 4 offsets in a grid
/// cell, the error is that bias to within the harmonics the offsets do not cancel, those at
/// multiples of four per cell, whose share is about the alias ratio at `p = 4` over the one at
/// `p = 1`, below `4⁻ⁿ`: held to twice that, at α = 0.25 Å⁻¹, where the vectors beyond the grid are
/// nothing, and at 0.5 Å⁻¹, where at order 6 they are a part of the bias many times the bound.
#[test]
fn a_lone_charge_is_the_wigner_lattice_on_the_grid_and_the_bias_off_it() {
    const XI: f64 = -2.837_297_479_480_62;
    let side = 15.0 * ANGSTROM;
    let cell = PeriodicBox::cubic(side);
    let exact = COULOMB * XI / (2.0 * side);
    let p = EwaldParameters::for_accuracy(&cell, 7.5 * ANGSTROM, 1e-9);
    for order in [3usize, 4, 5, 6, 8] {
        let ewald = Ewald::new(cell, p).with_mesh(mesh(order, [32; 3]));
        let h = side / 32.0;
        let e = ewald
            .evaluate(&[1.0], &[[3.0 * h, 5.0 * h, 0.0]], &[])
            .energy
            .total;
        println!(
            "lone charge on a grid point, order {order}: E/exact − 1 = {:.2e}",
            e / exact - 1.0
        );
        // The nearest image is 2 r_c away, erfc(α L) below 1e-30, and the vectors beyond the grid
        // weigh e⁻⁴⁹: what is left is the rounding of the parts, the self term 2.3 times the total.
        assert!(
            (e / exact - 1.0).abs() <= 1e-12,
            "order {order}: {e} against {exact}"
        );
        let coarse = Ewald::new(cell, p).with_mesh(mesh(order, [16; 3]));
        let h = side / 16.0;
        let e = coarse
            .evaluate(&[1.0], &[[3.0 * h, 5.0 * h, 0.0]], &[])
            .energy
            .total;
        let left_out = omitted_vectors(&cell, p.alpha, [16; 3]);
        println!(
            "    and on 16³: E/exact − 1 = {:.2e}, with the vectors beyond the grid back {:.1e}",
            e / exact - 1.0,
            (e + left_out) / exact - 1.0
        );
        assert!(
            left_out > 1e-6 * exact.abs(),
            "the vectors beyond the grid count"
        );
        assert!(
            ((e + left_out) / exact - 1.0).abs() <= 1e-12,
            "order {order} on 16³"
        );
    }
    let cell = PeriodicBox::new([15.0, 12.0, 18.0].map(|x| x * ANGSTROM));
    let l = cell.lengths();
    for (alpha, order) in [3usize, 4, 5, 6, 8]
        .map(|n| (0.25 / ANGSTROM, n))
        .into_iter()
        .chain([(0.5 / ANGSTROM, 4), (0.5 / ANGSTROM, 6)])
    {
        let p = converged(alpha, 6.0 * ANGSTROM);
        let reference = Ewald::new(cell, p)
            .evaluate(&[1.0], &[[0.0; 3]], &[])
            .energy
            .reciprocal;
        let ewald = Ewald::new(cell, p).with_mesh(mesh(order, [16; 3]));
        let mut sum = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    let t = [i, j, k];
                    let x = [0, 1, 2].map(|a| (3.0 + t[a] as f64 / 4.0) * l[a] / 16.0);
                    sum += ewald.evaluate(&[1.0], &[x], &[]).energy.reciprocal - reference;
                }
            }
        }
        let mean = sum / 64.0;
        let bias = ewald.mesh().expect("a mesh").self_bias(1.0);
        let bound = 2.0 / power(4.0, order);
        println!(
            "lone charge off the grid, α {:.2} /Å, order {order}: mean error / bias − 1 = {:.2e} \
             (bound {bound:.1e}); the vectors beyond the grid {:.2} of the bias",
            alpha * ANGSTROM,
            mean / bias - 1.0,
            -omitted_vectors(&cell, alpha, [16; 3]) / bias
        );
        assert!((mean / bias - 1.0).abs() <= bound, "order {order}");
    }
}

/// `log₂` of an error's fall when the grid doubles, from `(coarse, fine)`.
fn rate(coarse: f64, fine: f64) -> f64 {
    (coarse / fine).log2()
}

/// **The convergence rates the order predicts, measured.** The energy against the converged
/// classical sum, from 32³ to 64³ at α = 0.35 Å⁻¹, on crystals moved off the grid. The alias ratios
/// are `O(xⁿ)`, and for odd `n` the `±p` pair cancels in their sum `R`, leaving `η ≈ 2i xⁿ sin 2πu`
/// imaginary to leading order. **Where every charge has the same offset in its cell** — each
/// charge's own term, and NaCl here, whose ions are all a whole number of grid spacings apart — the
/// error is a common factor `|Π(1 + η_a)|² − 1`, whose leading term is `|η|²` or `Re η`, and the
/// energy falls as `K^{−2⌈n/2⌉}`: 4, 4, 6, 6, 8 for orders 3, 4, 5, 6, 8, measured 3.86, 3.89, 5.88,
/// 5.94, 8.00, held to `[r − 0.4, r + 0.6]`. **With several offsets the pairing is not
/// guaranteed**: CsCl's three offsets measure 4.03, 4.08, 6.07, 6.19, 8.31 here, but a crystal of four
/// ions at generic fractions measures 3.6–4.5 for order 4 and 4.5–6.2 for order 5 over doublings to
/// 128³ (an independent numpy implementation measured about `n`, 3.24, 4.12, 5.54, 6.17 for 3 to 6),
/// so for CsCl only the even orders' `n` is asserted, and the odd orders are printed. A crystal's
/// error sits at its Bragg vectors, whose `x = m/K` does not depend on α: from 16³ to 32³ the rates
/// are 4.4–6.7 at α = 0.35 and at 0.2 Å⁻¹ alike, not yet asymptotic.
///
/// The forces, on eight disordered boxes of 40 charges, RMS over every component of every box,
/// fall as `K^{−(n−1)}` at every order: differentiating the alias terms brings down `p/x`. Measured
/// 1.97, 3.10, 4.03, 5.16 and 7.36 for orders 3, 4, 5, 6 and 8 from 16³ to 32³ in a 12 Å cube at
/// α = 0.2 Å⁻¹, and 1.92–2.10, 3.00–3.11, 3.97–4.14, 5.02–5.16 and 6.0–6.3 for order 7, on every
/// doubling to 128³ — asymptotic; held to `[n − 1.3, n − 0.2]`, which excludes both neighbouring
/// orders.
#[test]
fn the_errors_fall_at_the_rates_the_order_predicts() {
    let r0 = 2.8 * ANGSTROM;
    let alpha = 0.2 / ANGSTROM;
    let shift = [0.37, 0.11, 0.23].map(|x| x * ANGSTROM);
    for (name, basis, a, m) in [
        ("NaCl", rock_salt(), 2.0 * r0, [2usize; 3]),
        ("CsCl", caesium_chloride(), 2.0 * r0 / 3f64.sqrt(), [3; 3]),
    ] {
        let (q, at, cell) = crystal(&basis, a, m, shift);
        let p = converged(0.35 / ANGSTROM, cell.half_shortest_edge());
        let reference = Ewald::new(cell, p).evaluate(&q, &at, &[]).energy.reciprocal;
        for order in [3usize, 4, 5, 6, 8] {
            let error = |k: usize| {
                let e = Ewald::new(cell, p)
                    .with_mesh(mesh(order, [k; 3]))
                    .evaluate(&q, &at, &[])
                    .energy
                    .reciprocal;
                (e - reference).abs()
            };
            let (coarse, fine) = (error(32), error(64));
            let r = rate(coarse, fine);
            let predicted = (order.div_ceil(2) * 2) as f64;
            let asserted = name == "NaCl" || order % 2 == 0;
            println!(
                "{name}, order {order}: energy error {coarse:.2e} → {fine:.2e} J, rate {r:.2} \
                 ({}{predicted})",
                if asserted {
                    "predicted "
                } else {
                    "not asserted; the pairing would give "
                }
            );
            if asserted {
                assert!(
                    (predicted - 0.4..=predicted + 0.6).contains(&r),
                    "{name} order {order}: rate {r}"
                );
            }
        }
    }
    for order in [3usize, 4, 5, 6, 8] {
        let (mut coarse, mut fine) = (0.0, 0.0);
        for seed in 0..8u64 {
            let (q, at, cell) = disordered(40, 12.0 * ANGSTROM, 0x5A7E + seed);
            let p = converged(alpha, 6.0 * ANGSTROM);
            let reference = Ewald::new(cell, p).evaluate(&q, &at, &[]).forces;
            for (k, total) in [(16usize, &mut coarse), (32, &mut fine)] {
                let f = Ewald::new(cell, p)
                    .with_mesh(mesh(order, [k; 3]))
                    .evaluate(&q, &at, &[])
                    .forces;
                for (x, y) in f.iter().zip(&reference) {
                    for c in 0..3 {
                        *total += (x[c] - y[c]) * (x[c] - y[c]);
                    }
                }
            }
        }
        let r = 0.5 * rate(coarse, fine);
        println!(
            "forces, order {order}: rate {r:.2} (predicted {})",
            order - 1
        );
        let n = order as f64;
        assert!((n - 1.3..=n - 0.2).contains(&r), "order {order}: rate {r}");
    }
}

/// The disordered boxes the estimates are measured on: 96 charges in 18 × 15.5 × 21 Å, so that each
/// axis's length is its own, and 48 of them.
const ESTIMATE_BOXES: u64 = 48;

/// **The estimates against measurement**, over 48 disordered boxes of 96 charges in an
/// 18 × 15.5 × 21 Å box, `r_c` = 7.5 Å at `δ = 10⁻⁶`: the mean error of the reciprocal energy is
/// [`Pme::self_bias`], and the RMS about it is [`Pme::energy_error`], each derived in the module
/// documentation and neither fitted.
///
/// - **The mean** is held to 0.85–1.15 of the bias; measured 0.986–1.023 here. Its standard error is
///   `σ/(√N |b|)`, under 2% at 48 boxes: the 0.95–0.98 sixteen boxes in a cube gave was 1.3 of
///   their standard errors low, not a systematic shortfall.
/// - **The RMS about it**: `N` boxes measure an RMS to a relative standard error of `1/√(2N)`,
///   10.2% at 48, so each case is held to `1 ± 3/√(2N)`, 0.69–1.31; and the geometric mean of the
///   six, whose error is `1/√(2N)/√6` if they were independent — they are not, sharing their boxes,
///   which makes the band generous rather than tight — to `1 ± 3/√(12N)`, 0.875–1.125. A σ wrong by
///   a factor 1.5 either way is outside that.
/// - **The bias is the larger part of the error at even orders**, where it and σ both fall as
///   `K⁻ⁿ` and their ratio is fixed (measured 2.2–2.4 on 40 charges, 10–14 on 96): held to three
///   times σ there. At odd orders the bias falls as `K^{−(n+1)}` and σ more slowly, so the ratio
///   falls with the grid and is not asserted.
#[test]
fn the_bias_and_the_spread_are_the_error_on_disordered_boxes() {
    let cases = [
        (3usize, 32usize),
        (4, 16),
        (4, 32),
        (5, 16),
        (6, 16),
        (8, 16),
    ];
    let mut sum = vec![0.0; cases.len()];
    let mut sum_sq = vec![0.0; cases.len()];
    let mut estimates = vec![(0.0, 0.0); cases.len()];
    for seed in 0..ESTIMATE_BOXES {
        let (q, at, cell) =
            disordered_in(96, [18.0, 15.5, 21.0].map(|x| x * ANGSTROM), 0xC0DE + seed);
        let sum_q2: f64 = q.iter().map(|x| x * x).sum();
        let p = EwaldParameters::for_accuracy(&cell, 7.5 * ANGSTROM, 1e-6);
        let reference = Ewald::new(cell, converged(p.alpha, p.cutoff))
            .evaluate(&q, &at, &[])
            .energy
            .reciprocal;
        for (c, &(order, k)) in cases.iter().enumerate() {
            let ewald = Ewald::new(cell, p).with_mesh(mesh(order, [k; 3]));
            let e = ewald.evaluate(&q, &at, &[]).energy.reciprocal;
            let m = ewald.mesh().expect("a mesh");
            let (bias, sigma) = (m.self_bias(sum_q2), m.energy_error(sum_q2));
            estimates[c] = (bias, sigma);
            sum[c] += e - reference;
            sum_sq[c] += (e - reference - bias) * (e - reference - bias);
        }
    }
    let n = ESTIMATE_BOXES as f64;
    let band = 3.0 / (2.0 * n).sqrt();
    let mut log_sum = 0.0;
    for (c, &(order, k)) in cases.iter().enumerate() {
        let (bias, sigma) = estimates[c];
        let mean = sum[c] / n / bias;
        let rms = (sum_sq[c] / n).sqrt() / sigma;
        println!(
            "order {order}, {k}³: bias {bias:.3e} J, σ {sigma:.3e} J, {:.1} σ; mean error {mean:.3} \
             of the bias, RMS about it {rms:.3} of σ",
            bias.abs() / sigma
        );
        assert!(
            (0.85..1.15).contains(&mean),
            "order {order} {k}³: mean {mean}"
        );
        assert!(
            (1.0 - band..1.0 + band).contains(&rms),
            "order {order} {k}³: RMS {rms} outside 1 ± {band}"
        );
        if order % 2 == 0 {
            assert!(bias.abs() > 3.0 * sigma, "order {order} {k}³");
        }
        log_sum += rms.ln();
    }
    let pooled = (log_sum / cases.len() as f64).exp();
    let pooled_band = 3.0 / (2.0 * n * cases.len() as f64).sqrt();
    println!("the six cases' geometric mean: {pooled:.3} of σ (band 1 ± {pooled_band:.3})");
    assert!((1.0 - pooled_band..1.0 + pooled_band).contains(&pooled));
}

/// The estimated error [`PmeParameters::for_accuracy`] holds to its target, `(σ² + b²)^½`, in units
/// of `k_e Q / r_c`.
fn targeted(cell: &PeriodicBox, p: &EwaldParameters, m: PmeParameters) -> f64 {
    let pme = Pme::new(*cell, p.alpha, m);
    let (s, b) = (pme.energy_error(1.0), pme.self_bias(1.0));
    (s * s + b * b).sqrt() * p.cutoff / COULOMB
}

/// **[`PmeParameters::for_accuracy`] delivers its accuracy**: at `δ` = 10⁻⁴ and 10⁻⁵, `r_c` = 7.5 Å,
/// on sixteen disordered boxes of 96 charges in 18 × 15.5 × 21 Å, the RMS error of the reciprocal
/// energy at the order and grid it chooses — the error itself, about zero, its mean included —
/// is at most `δ k_e Q / r_c`, with the `1 + 3/√(2N)` an RMS over `N` boxes needs. The grid is a
/// power of two, so the estimate is under `δ` by whatever the last doubling left.
#[test]
fn the_chosen_mesh_delivers_its_accuracy() {
    for accuracy in [1e-4, 1e-5] {
        // Every box has the same edges and the same charges, so one choice serves them all.
        let cell = PeriodicBox::new([18.0, 15.5, 21.0].map(|x| x * ANGSTROM));
        let p = EwaldParameters::for_accuracy(&cell, 7.5 * ANGSTROM, accuracy);
        let m = PmeParameters::for_accuracy(&cell, &p, accuracy, 96);
        let mut sum_sq = 0.0;
        let mut unit = 0.0;
        for seed in 0..16u64 {
            let (q, at, cell) =
                disordered_in(96, [18.0, 15.5, 21.0].map(|x| x * ANGSTROM), 0xDE11 + seed);
            let sum_q2: f64 = q.iter().map(|x| x * x).sum();
            let reference = Ewald::new(cell, converged(p.alpha, p.cutoff))
                .evaluate(&q, &at, &[])
                .energy
                .reciprocal;
            let e = Ewald::new(cell, p)
                .with_mesh(m)
                .evaluate(&q, &at, &[])
                .energy
                .reciprocal;
            sum_sq += (e - reference) * (e - reference);
            unit = COULOMB * sum_q2 / p.cutoff;
        }
        let rms = (sum_sq / 16.0).sqrt() / unit;
        let estimate = targeted(&cell, &p, m);
        println!(
            "δ = {accuracy:e}: chose {m:?}, estimate {estimate:.2e}; the error's RMS {rms:.2e}, \
             {:.2} of δ",
            rms / accuracy
        );
        assert!(
            rms <= accuracy * (1.0 + 3.0 / 32f64.sqrt()),
            "δ = {accuracy:e}"
        );
    }
}

/// [`PmeParameters::for_accuracy`] reaches what it is asked for and no more: the grid it chooses
/// for each order has its estimate `(σ² + b²)^½` at or under the target, its spacings within a
/// factor of two of each other, and the same grid with the axis it doubled last halved again does
/// not reach it; and of the five orders it takes the one [`PmeParameters::cost`] puts cheapest. On
/// a box of 20 × 18 × 23 Å, `r_c` = 8 Å, at `δ` = 10⁻⁴ and 10⁻⁵.
#[test]
fn the_chosen_grid_is_the_smallest_that_reaches_the_target() {
    let cell = PeriodicBox::new([20.0, 18.0, 23.0].map(|x| x * ANGSTROM));
    let l = cell.lengths();
    for accuracy in [1e-4, 1e-5] {
        let p = EwaldParameters::for_accuracy(&cell, 8.0 * ANGSTROM, accuracy);
        let estimate = |m: PmeParameters| targeted(&cell, &p, m);
        for order in 4..=8 {
            let grid = PmeParameters::grid_for(&cell, &p, accuracy, order);
            let e = estimate(mesh(order, grid));
            assert!(e <= accuracy, "order {order}: {e:e} against {accuracy:e}");
            // Doubling the coarsest axis keeps every spacing within a factor of two of the others.
            let spacing = [0, 1, 2].map(|a| l[a] / grid[a] as f64);
            let (lo, hi) = spacing
                .iter()
                .fold((f64::INFINITY, 0.0f64), |(lo, hi), &x| {
                    (lo.min(x), hi.max(x))
                });
            assert!(hi <= 2.0 * lo, "order {order}: spacings {spacing:?}");
            // The axis doubled last is the one with the finest spacing, the first such on a tie.
            let mut finest = 0;
            for a in 1..3 {
                if l[a] / (grid[a] as f64) < l[finest] / (grid[finest] as f64) {
                    finest = a;
                }
            }
            let mut smaller = grid;
            smaller[finest] /= 2;
            if smaller[finest] >= order.next_power_of_two() {
                let s = estimate(mesh(order, smaller));
                assert!(
                    s > accuracy,
                    "order {order}: {smaller:?} reaches {s:e} already"
                );
            }
            println!("δ = {accuracy:e}, order {order}: {grid:?}, estimate {e:.2e}");
        }
        // The choice minimises the cost model over the five orders' grids, at both sizes.
        for atoms in [3000, 24_000] {
            let chosen = PmeParameters::for_accuracy(&cell, &p, accuracy, atoms);
            let costs: Vec<(usize, f64)> = (4..=8)
                .map(|n| {
                    (
                        n,
                        mesh(n, PmeParameters::grid_for(&cell, &p, accuracy, n)).cost(atoms),
                    )
                })
                .collect();
            println!("δ = {accuracy:e}, {atoms} atoms: chose {chosen:?}; costs {costs:?}");
            assert!(estimate(chosen) <= accuracy);
            assert_eq!(
                chosen.grid,
                PmeParameters::grid_for(&cell, &p, accuracy, chosen.order)
            );
            for &(n, c) in &costs {
                assert!(
                    chosen.cost(atoms) <= c,
                    "order {n} is cheaper: {c} against {}",
                    chosen.cost(atoms)
                );
            }
        }
    }
}

/// A central difference and the tolerance it earns: three times the truncation term is
/// `|D(2h) − D(h)|`, and rounding adds `ε_E / h` with `ε_E` the energy's rounding, bounded by `64 ε`
/// times the sum of the parts' magnitudes.
fn finite_difference(energy: impl Fn(f64) -> (f64, f64), h: f64) -> (f64, f64) {
    let d = |s: f64| {
        let (plus, scale) = energy(s);
        let (minus, _) = energy(-s);
        ((plus - minus) / (2.0 * s), scale)
    };
    let (d1, scale) = d(h);
    let (d2, _) = d(2.0 * h);
    (d1, 2.0 * (d2 - d1).abs() / 3.0 + 64.0 * EPS * scale / h)
}

fn magnitude(e: &EwaldEnergy) -> f64 {
    e.real.abs() + e.reciprocal.abs() + e.self_energy.abs() + e.excluded.abs() + e.background.abs()
}

/// **The forces are the gradient of the mesh's energy** — not of the classical sum's, which they
/// differ from by the mesh's error — against central differences of it, at orders 4 and 6 in a
/// cube and at order 5 in a box of 14 × 12 × 16 Å on a grid of 16 × 8 × 32, where each axis's
/// `K_a/L_a` is its own, with a net charge (the background) and exclusions (their correction), each
/// to the tolerance the differences earn and that tolerance under a part in 10⁵ of the force.
#[test]
fn the_forces_are_the_negative_gradient_of_the_mesh_energy() {
    for (order, grid, lengths, cutoff) in [
        (4usize, [16usize; 3], [14.0; 3], 7.0),
        (6, [32; 3], [14.0; 3], 7.0),
        (5, [16, 8, 32], [14.0, 12.0, 16.0], 6.0),
    ] {
        let (mut q, at, cell) = disordered_in(40, lengths.map(|x| x * ANGSTROM), 0xF0CE);
        q[0] += 0.3;
        let excluded = [[0, 1], [2, 3], [3, 4], [2, 4], [10, 30]];
        let p = EwaldParameters::for_accuracy(&cell, cutoff * ANGSTROM, 1e-6);
        let ewald = Ewald::new(cell, p).with_mesh(mesh(order, grid));
        let ev = ewald.evaluate(&q, &at, &excluded);
        let h = 1e-5 * ANGSTROM;
        let mut worst = 0.0f64;
        for i in [0, 1, 2, 7, 10, 23, 39] {
            for a in 0..3 {
                let energy = |s: f64| {
                    let mut moved = at.clone();
                    moved[i][a] += s;
                    let e = ewald.evaluate(&q, &moved, &excluded).energy;
                    (e.total, magnitude(&e))
                };
                let (d, tolerance) = finite_difference(energy, h);
                let f = ev.forces[i][a];
                assert!(
                    (f + d).abs() <= tolerance,
                    "order {order}: atom {i} axis {a}: force {f:e}, −dE/dx {:e}, tolerance \
                     {tolerance:e}",
                    -d
                );
                assert!(tolerance < 1e-5 * f.abs(), "{tolerance:e} against {f:e}");
                worst = worst.max((f + d).abs() / tolerance);
            }
        }
        println!("order {order}, {grid:?}: worst |F + dE/dx| is {worst:.3} of its tolerance");
    }
}

/// **The net force is not zero, and it is exactly the energy's derivative under a rigid
/// translation**: moving every charge by the same vector changes where each sits in its grid
/// cell, so the mesh's energy depends on the system's place and `Σ_i F_i = −dE/dd`. Against a
/// central difference along each axis, to the tolerance it earns; measured against the RMS force,
/// it falls as the grid is refined, and the classical sum's is rounding.
#[test]
fn the_net_force_is_the_derivative_under_translation() {
    let (q, at, cell) = disordered(40, 12.0 * ANGSTROM, 0x7E57);
    let p = converged(0.3 / ANGSTROM, 6.0 * ANGSTROM);
    let classical = Ewald::new(cell, p).evaluate(&q, &at, &[]);
    let rms =
        |f: &[[f64; 3]]| (f.iter().flatten().map(|x| x * x).sum::<f64>() / f.len() as f64).sqrt();
    let net = |f: &[[f64; 3]]| [0, 1, 2].map(|a| f.iter().map(|x| x[a]).sum::<f64>());
    let largest = |v: [f64; 3]| v.iter().map(|x| x.abs()).fold(0.0, f64::max);
    let scale = rms(&classical.forces);
    let n0 = largest(net(&classical.forces));
    println!("classical: net force {:.1e} of the RMS force", n0 / scale);
    let mut last = f64::INFINITY;
    for k in [16usize, 32] {
        let ewald = Ewald::new(cell, p).with_mesh(mesh(4, [k; 3]));
        let ev = ewald.evaluate(&q, &at, &[]);
        let total = net(&ev.forces);
        for a in 0..3 {
            let energy = |s: f64| {
                let moved: Vec<[f64; 3]> = at
                    .iter()
                    .map(|r| {
                        let mut r = *r;
                        r[a] += s;
                        r
                    })
                    .collect();
                let e = ewald.evaluate(&q, &moved, &[]).energy;
                (e.total, magnitude(&e))
            };
            let (d, tolerance) = finite_difference(energy, 1e-4 * ANGSTROM);
            println!(
                "{k}³, axis {a}: Σ F = {:.6e} N, −dE/dd = {:.6e} ± {tolerance:.1e}; {:.2e} of the \
                 RMS force",
                total[a],
                -d,
                total[a].abs() / scale
            );
            assert!((total[a] + d).abs() <= tolerance, "{k}³ axis {a}");
            assert!(tolerance < 1e-3 * total[a].abs(), "the check can see it");
        }
        let size = largest(total);
        assert!(size < last, "{k}³: the net force did not fall");
        assert!(size > 1e3 * n0, "the mesh's net force is not rounding");
        last = size;
    }
}

/// **Translation moves the energy by the mesh's error, at its rate**: shifting every charge by a
/// fraction of a cell, the RMS change over eight boxes falls as the energy's spread does. **At even
/// orders that is `K⁻ⁿ`, asymptotic**: measured 4.03, 4.05, 4.01 for order 4 and 6.01, 6.12, 6.01 for
/// order 6 over the doublings from 16³ to 128³ at α = 0.2 Å⁻¹, 7.83 for order 8 from 16³ to 32³;
/// held to `[n − 0.4, n + 0.6]` from 16³ to 32³. **At odd orders it is not settled and is not
/// asserted**: 4.20 and 6.29 for orders 3 and 5 from 16³ to 32³, where the bias's `K^{−(n+1)}` is
/// still the larger part, then 2.0–2.8, 3.9–4.6 and 5.9–6.6 for orders 3, 5 and 7 on the doublings to
/// 128³ at α = 0.2 and 0.1 Å⁻¹ — between `n − 1` and `n`, and printed. The classical sum moves by
/// rounding.
#[test]
fn translation_moves_the_energy_by_the_meshs_error_at_its_rate() {
    let alpha = 0.2 / ANGSTROM;
    for order in [3usize, 4, 5, 6, 8] {
        let mut moved_by = [0.0; 2];
        for seed in 0..8u64 {
            let (q, at, cell) = disordered(40, 12.0 * ANGSTROM, 0x5A7E + seed);
            let p = converged(alpha, 6.0 * ANGSTROM);
            for (slot, k) in [16usize, 32].into_iter().enumerate() {
                let shift = [0.31, 0.77, 0.52].map(|x| x * 12.0 / k as f64 * ANGSTROM);
                let moved: Vec<[f64; 3]> = at
                    .iter()
                    .map(|r| [0, 1, 2].map(|a| r[a] + shift[a]))
                    .collect();
                let ewald = Ewald::new(cell, p).with_mesh(mesh(order, [k; 3]));
                let a = ewald.evaluate(&q, &at, &[]).energy.reciprocal;
                let b = ewald.evaluate(&q, &moved, &[]).energy.reciprocal;
                moved_by[slot] += (b - a) * (b - a);
                if order == 4 && seed == 0 && k == 16 {
                    let c = Ewald::new(cell, p);
                    let d = c.evaluate(&q, &moved, &[]).energy.reciprocal
                        - c.evaluate(&q, &at, &[]).energy.reciprocal;
                    println!("classical: moved by {:.1e} of itself", d.abs() / a.abs());
                    assert!(d.abs() <= 1e-12 * a.abs());
                }
            }
        }
        let r = 0.5 * rate(moved_by[0], moved_by[1]);
        let n = order as f64;
        if order % 2 == 0 {
            println!("order {order}: translation moves the energy at rate {r:.2} (predicted {n})");
            assert!((n - 0.4..=n + 0.6).contains(&r), "order {order}: {r}");
        } else {
            println!("order {order}: translation moves the energy at rate {r:.2} (not asserted)");
        }
    }
}

/// **The virial is the strain derivative of the mesh's energy**, with the real space switched
/// off (`r_c` = 0.5 Å, where no pair is: they are at least 1.2 Å apart) so that what is compared is
/// the reciprocal part alone. Not smaller: the cell list allocates `(2L/r_c)³` cells, 83.6 GB at
/// 0.01 Å.
/// The diagonal against central differences of the mesh's own energy in a stretched box, the same
/// grid, every fractional coordinate fixed — exactly the strain, since `b(m)` depends on integers.
/// The off-diagonal against central differences of [`series_energy`] under a simple shear, which
/// the orthorhombic mesh cannot represent, and which is the mesh's energy (the series test). On a
/// 16³ grid and on 16 × 8 × 32, where each axis's own grid is in each vector.
#[test]
fn the_virial_is_the_strain_and_shear_derivative_of_the_mesh_energy() {
    let (q, at, cell) = disordered_in(24, [11.0, 9.5, 12.5].map(|x| x * ANGSTROM), 0x51A1);
    let alpha = 0.3 / ANGSTROM;
    let p = EwaldParameters {
        alpha,
        cutoff: 0.5 * ANGSTROM,
        k_cutoff: alpha,
    };
    for grid in [mesh(6, [16, 16, 16]), mesh(6, [16, 8, 32])] {
        let ewald = Ewald::new(cell, p).with_mesh(grid);
        let ev = ewald.evaluate(&q, &at, &[]);
        assert_eq!(ev.energy.real, 0.0, "no pair inside the cutoff");
        let w = ev.virial;
        let h = 1e-6;
        for a in 0..3 {
            let energy = |s: f64| {
                let mut f = [1.0; 3];
                f[a] += s;
                let scaled: Vec<[f64; 3]> =
                    at.iter().map(|r| [0, 1, 2].map(|b| r[b] * f[b])).collect();
                let e = ewald
                    .with_cell(cell.scaled(f))
                    .evaluate(&q, &scaled, &[])
                    .energy;
                (e.reciprocal, e.reciprocal.abs())
            };
            let (d, tolerance) = finite_difference(energy, h);
            println!(
                "W_{a}{a} = {:.8e} J, −dE/dε = {:.8e} ± {tolerance:.1e}",
                w[a][a], -d
            );
            assert!((w[a][a] + d).abs() <= tolerance, "axis {a}");
            assert!(tolerance < 1e-6 * w[a][a].abs());
        }
        for (a, b) in [(0, 1), (0, 2), (1, 2)] {
            let energy = |s: f64| {
                let e = series_energy(&q, &at, &cell, alpha, grid, Some((a, b, s)));
                (e, e.abs())
            };
            let (d, tolerance) = finite_difference(energy, h);
            // dU/dε = Σ ∂U/∂r_a r_b = −W_ba.
            let wab = w[b][a];
            println!(
                "W_{a}{b} = {wab:.8e} J, −dE/dε = {:.8e} ± {tolerance:.1e}",
                -d
            );
            assert!((wab + d).abs() <= tolerance, "W_{a}{b}");
            assert!(tolerance < 1e-4 * wab.abs(), "the check can see W_{a}{b}");
            assert!(
                (w[a][b] - w[b][a]).abs() <= 1e-12 * w[0][0].abs(),
                "symmetric"
            );
        }
    }
}

fn periodic_magnitude(e: &PeriodicEnergy) -> f64 {
    let w = &e.ewald;
    e.bond.abs()
        + e.angle.abs()
        + e.torsion.abs()
        + e.inversion.abs()
        + e.van_der_waals.abs()
        + e.dispersion_correction.abs()
        + w.real.abs()
        + w.reciprocal.abs()
        + w.self_energy.abs()
        + w.excluded.abs()
        + w.background.abs()
}

/// Two TIP3P waters decoupled from 25, through the mesh: the box of [`WaterBox::lattice`] with
/// +0.1 e on the first water and −0.1 e on the third, so that the group and the rest are each
/// charged and the background has a cross term, on a grid of 16 × 8 × 32, so that each axis's own
/// size is in the split of `C(m)` and `C(−m)`.
fn decoupled() -> (PeriodicDecoupling, PeriodicForceField, Vec<[f64; 3]>) {
    let water = WaterBox::lattice(3, 0xDEC0);
    let field = water.force_field(4.6 * ANGSTROM, 1e-8);
    let mut q = field.charges().to_vec();
    q[0] += 0.1;
    q[6] -= 0.1;
    let field = field.with_charges(q).with_mesh(mesh(6, [16, 8, 32]));
    let mut group = vec![false; water.positions().len()];
    group[..6].fill(true);
    let d = PeriodicDecoupling::new(&field, &group).expect("two waters decouple");
    (d, field, water.positions().to_vec())
}

/// **A decoupling through the mesh.** The reciprocal cross term is one complex grid carrying the
/// rest in its real part and the group in its imaginary — `E(all) − E(rest) − E(group)` exactly,
/// because the energy is a quadratic form in the grid's charges. Checked as:
///
/// - `∂U/∂λ_e`, the cross Coulomb energy, is three whole Ewald sums of the box with charges taken
///   away, `E(all) − E(rest) − E(group)`, to `1e-12` of their parts, at three states;
/// - **the decoupled reciprocal energy itself** is the rest's own whole sum plus `λ_e` times that
///   cross term, at `λ_e` = 1, 0.4 and 0, to `1e-12` of the parts. The checks above are all
///   differences or symmetric in the two sets, so they cannot tell the rest from the group: with
///   the two swapped in the spread and the gather, the energy at `λ_e = 1` is 94% wrong, and every
///   one of them passes;
/// - the couplings at five states differ as the whole energy does, `U(λ) − U(λ₀)` from five
///   evaluations of the box, to `1e-12` of the parts;
/// - `∂U/∂λ_e` and `∂U/∂λ_v` against central differences of the coupling, and the forces at an
///   intermediate state, on atoms of the group and of the rest, against central differences of
///   the energy, each to the tolerance it earns.
#[test]
fn a_decoupling_through_the_mesh_is_its_whole_energy() {
    let (d, field, at) = decoupled();
    assert!(field.ewald().mesh().is_some(), "the mesh is in use");
    let q = field.charges().to_vec();
    let mut q_rest = q.clone();
    q_rest[..6].fill(0.0);
    let mut q_group = vec![0.0; q.len()];
    q_group[..6].copy_from_slice(&q[..6]);
    let ewald = |charges: &[f64]| {
        field
            .clone()
            .with_charges(charges.to_vec())
            .energy(&at)
            .ewald
    };
    let (all, rest, group) = (ewald(&q), ewald(&q_rest), ewald(&q_group));
    let cross = all.total - rest.total - group.total;
    let scale = magnitude(&all);
    for lambda in [
        Lambda::COUPLED,
        Lambda::new(0.0, 0.4, 0.7),
        Lambda::new(0.0, 0.0, 0.2),
    ] {
        let c = d.coupling(&at, lambda);
        println!(
            "λ = {lambda:?}: ∂U/∂λ_e = {:.12e} J, the sums' cross terms {cross:.12e}",
            c.gradient[1]
        );
        assert!((c.gradient[1] - cross).abs() <= 1e-12 * scale);
    }
    println!(
        "the cross terms are {:.2e} of the sums' parts",
        cross.abs() / scale
    );
    assert!(
        cross.abs() > 1e-5 * scale,
        "the cross terms are not negligible"
    );
    let cross_reciprocal = all.reciprocal - rest.reciprocal - group.reciprocal;
    for lambda in [
        Lambda::COUPLED,
        Lambda::new(0.0, 0.4, 1.0),
        Lambda::new(0.0, 0.0, 1.0),
    ] {
        let e = d.energy(&at, lambda).ewald.reciprocal;
        let want = rest.reciprocal + lambda.electrostatics * cross_reciprocal;
        println!(
            "λ_e = {}: the decoupled reciprocal energy {e:.12e} J, the rest's own sum plus λ_e \
             times the cross term {want:.12e}, off {:.1e} of the parts",
            lambda.electrostatics,
            (e - want).abs() / scale
        );
        assert!(
            (e - want).abs() <= 1e-12 * scale,
            "λ_e = {}",
            lambda.electrostatics
        );
    }
    assert!(
        (rest.reciprocal - all.reciprocal).abs() > 1e-3 * scale,
        "the rest's own sum is not the whole box's"
    );

    let whole_scale = periodic_magnitude(&field.energy(&at));
    let states = [
        Lambda::COUPLED,
        Lambda::new(0.0, 0.5, 1.0),
        Lambda::new(0.0, 0.0, 1.0),
        Lambda::new(0.0, 0.0, 0.4),
        Lambda::new(0.0, 0.0, 0.0),
    ];
    let couplings = d.couplings(&at, &states);
    let energy = |l: Lambda| {
        let mut f = vec![[0.0; 3]; at.len()];
        d.energy_and_forces(&at, l, &mut f)
    };
    let u0 = energy(states[0]);
    for (l, c) in states.iter().zip(&couplings).skip(1) {
        let du = energy(*l) - u0;
        let dc = c.energy - couplings[0].energy;
        println!("{l:?}: U − U(coupled) = {du:.12e} J, from the couplings {dc:.12e}");
        assert!((du - dc).abs() <= 1e-12 * whole_scale, "{l:?}");
    }

    let lambda = Lambda::new(0.0, 0.45, 0.55);
    let c = d.coupling(&at, lambda);
    for (k, h) in [(1usize, 1e-3), (2, 1e-4)] {
        let energy = |s: f64| {
            let mut l = lambda;
            if k == 1 {
                l.electrostatics += s;
            } else {
                l.van_der_waals += s;
            }
            (d.coupling(&at, l).energy, whole_scale)
        };
        let (dl, tolerance) = finite_difference(energy, h);
        println!(
            "∂U/∂λ[{k}] = {:.10e}, by differences {dl:.10e} ± {tolerance:.1e}",
            c.gradient[k]
        );
        assert!((c.gradient[k] - dl).abs() <= tolerance);
        assert!(tolerance < 1e-3 * c.gradient[k].abs());
    }
    let mut forces = vec![[0.0; 3]; at.len()];
    d.energy_and_forces(&at, lambda, &mut forces);
    let h = 1e-5 * ANGSTROM;
    for i in [0, 2, 4, 6, 40] {
        for a in 0..3 {
            let energy = |s: f64| {
                let mut moved = at.clone();
                moved[i][a] += s;
                let mut f = vec![[0.0; 3]; at.len()];
                (d.energy_and_forces(&moved, lambda, &mut f), whole_scale)
            };
            let (dx, tolerance) = finite_difference(energy, h);
            assert!(
                (forces[i][a] + dx).abs() <= tolerance,
                "atom {i} axis {a}: {:e} against {:e} ± {tolerance:e}",
                forces[i][a],
                -dx
            );
            // Against the atom's force, not one component, which can be near zero.
            let size = forces[i].iter().map(|x| x * x).sum::<f64>().sqrt();
            assert!(
                tolerance < 1e-5 * size,
                "atom {i}: {tolerance:e} against {size:e}"
            );
        }
    }
}

/// **A water box through the mesh is the water box**: 64 lattice waters, the force field with
/// the mesh [`PmeParameters::grid_for`] gives at `δ = 10⁻⁵`, at every order
/// [`PmeParameters::for_accuracy`] considers. Every term but the reciprocal one is the classical
/// field's to the bit, and the reciprocal energy's error against the converged classical sum is
/// under `δ k_e Q / r_c` at every one of those grids — what the target promises.
///
/// **Neither estimate describes water's error, and neither is used as its tolerance.** Each water's
/// charges sit a bond apart, so most of each charge's own term is cancelled by its neighbours': the
/// error is a fraction of [`Pme::self_bias`] of either sign (measured −0.33 to 0.36 here, and −0.33
/// to 0.40 at `δ = 10⁻⁶`), and, being correlated, up to 3.8 of [`Pme::energy_error`]'s σ (3.9 at
/// `10⁻⁶`) — so a `4σ` tolerance would not be earned. Both ratios are printed. Measured, the error
/// is 0.04–0.24 of `δ`.
#[test]
fn a_water_box_through_the_mesh_is_the_classical_box() {
    let water = WaterBox::lattice(4, 0x3A7E);
    let at = water.positions().to_vec();
    let field = water.force_field(6.0 * ANGSTROM, 1e-5);
    let cell = field.cell();
    let p = field.ewald().parameters();
    let q = field.charges().to_vec();
    let sum_q2: f64 = q.iter().map(|x| x * x).sum();
    let reference = Ewald::new(cell, converged(p.alpha, p.cutoff))
        .evaluate(&q, &at, field.excluded_pairs())
        .energy
        .reciprocal;
    let a = field.energy(&at);
    println!(
        "classical at δ = 1e-5: off {:.2e} J, its reciprocal bias {:.2e}",
        a.ewald.reciprocal - reference,
        p.reciprocal_bias(sum_q2)
    );
    let chosen = PmeParameters::for_accuracy(&cell, &p, 1e-5, at.len());
    for order in 4..=8 {
        let grid = mesh(order, PmeParameters::grid_for(&cell, &p, 1e-5, order));
        let meshed = field.clone().with_mesh(grid);
        let b = meshed.energy(&at);
        let m = meshed.ewald().mesh().expect("a mesh");
        let (bias, sigma) = (m.self_bias(sum_q2), m.energy_error(sum_q2));
        let off = b.ewald.reciprocal - reference;
        println!(
            "{grid:?}{}: off {off:.2e} J, {:.2} of the bias {bias:.2e}, {:.2} of σ {sigma:.2e}",
            if grid == chosen { " (chosen)" } else { "" },
            off / bias,
            off / sigma
        );
        let unit = COULOMB * sum_q2 / p.cutoff;
        println!("    {:.3} of δ", off.abs() / unit / 1e-5);
        assert!(
            off.abs() <= 1e-5 * unit,
            "order {order}: the mesh delivers δ"
        );
        assert_eq!(a.van_der_waals.to_bits(), b.van_der_waals.to_bits());
        assert_eq!(a.ewald.real.to_bits(), b.ewald.real.to_bits());
        assert_eq!(a.ewald.self_energy.to_bits(), b.ewald.self_energy.to_bits());
        assert_eq!(a.ewald.excluded.to_bits(), b.ewald.excluded.to_bits());
        assert_eq!(a.ewald.background.to_bits(), b.ewald.background.to_bits());
        assert_eq!(
            a.dispersion_correction.to_bits(),
            b.dispersion_correction.to_bits()
        );
    }
}

/// The bits [`an_evaluation_is_its_own_bits_on_every_platform`] pins: the reciprocal energy, the
/// force `F[7][2]` and the virial `W[0][1]`.
const PIN: [u64; 3] = [0x3c6deaee1184bec7, 0x3e36a21d54f62ed3, 0xbc3624200079d4b3];

/// **The same bits, here and on every platform.** The mesh is arithmetic, `sqrt`, `floor` and
/// the kernel's `exp` — its twiddles are its own `sin` and `cos` — so an evaluation repeats to
/// the bit, a second mesh built the same way gives the same bits, and the reciprocal energy, a
/// force and a virial component of a fixed configuration are pinned to the bit as the platform
/// that wrote them computed them, which CI's Linux, macOS and WebAssembly jobs then hold. The
/// classical sum stays the default.
#[test]
fn an_evaluation_is_its_own_bits_on_every_platform() {
    let (q, at, cell) = disordered_in(30, [12.0, 11.0, 13.0].map(|x| x * ANGSTROM), 0xB175);
    let p = EwaldParameters::for_accuracy(&cell, 5.5 * ANGSTROM, 1e-6);
    assert!(
        Ewald::new(cell, p).mesh().is_none(),
        "the classical sum is the default"
    );
    let grid = mesh(5, [16, 16, 32]);
    let a = Ewald::new(cell, p)
        .with_mesh(grid)
        .evaluate(&q, &at, &[[0, 1]]);
    let b = Ewald::new(cell, p)
        .with_mesh(grid)
        .evaluate(&q, &at, &[[1, 0]]);
    assert_eq!(a, b);
    let bits = [
        a.energy.reciprocal.to_bits(),
        a.forces[7][2].to_bits(),
        a.virial[0][1].to_bits(),
    ];
    println!("pin: {bits:#018x?}");
    assert_eq!(bits, PIN);
}

/// The cost against classical Ewald, and what each is accurate to, on the W1/W2 water boxes of
/// 3 000 and 24 000 atoms at `r_c` = 9 Å, against the classical sum at `k_c = 9α`. The
/// reciprocal times are the whole sum's less real space's, each the best of seven. The energy
/// errors are as computed: neither sum's bias estimate is taken out, since for water — charges a
/// bond apart — neither describes it. Release, one
/// core: `cargo test --release -p pantometry-forcefield --test
/// the_particle_mesh_against_closed_forms the_cost -- --ignored --nocapture`. Reported, not
/// asserted.
#[test]
#[ignore = "a timing, release only"]
fn the_cost() {
    let best = |f: &mut dyn FnMut()| {
        let mut t = f64::INFINITY;
        for _ in 0..7 {
            let s = std::time::Instant::now();
            f();
            t = t.min(s.elapsed().as_secs_f64());
        }
        t
    };
    let rms = |a: &[[f64; 3]], b: &[[f64; 3]]| {
        let (mut d, mut f) = (0.0, 0.0);
        for (x, y) in a.iter().zip(b) {
            for c in 0..3 {
                d += (x[c] - y[c]) * (x[c] - y[c]);
                f += y[c] * y[c];
            }
        }
        (d / f).sqrt()
    };
    for per_side in [10usize, 20] {
        let water = WaterBox::lattice(per_side, 0xC057);
        let at = water.positions().to_vec();
        let cell = water.cell();
        let n = at.len();
        for accuracy in [1e-5, 1e-6] {
            let field = water.force_field(9.0 * ANGSTROM, accuracy);
            let q = field.charges().to_vec();
            let ex = field.excluded_pairs().to_vec();
            let p = field.ewald().parameters();
            let reference = Ewald::new(
                cell,
                EwaldParameters {
                    k_cutoff: 9.0 * p.alpha,
                    ..p
                },
            )
            .evaluate(&q, &at, &ex);
            let r = reference.energy.reciprocal;
            let classical = field.ewald().clone();
            let real_only = Ewald::new(
                cell,
                EwaldParameters {
                    k_cutoff: 1e-3 / ANGSTROM,
                    ..p
                },
            );
            let t_real = best(&mut || {
                std::hint::black_box(real_only.evaluate(&q, &at, &ex));
            });
            let t_classical = best(&mut || {
                std::hint::black_box(classical.evaluate(&q, &at, &ex));
            });
            let c = classical.evaluate(&q, &at, &ex);
            println!(
                "{n} atoms, {:.1} Å, δ = {accuracy:e}, α {:.3} /Å: real space and exclusions \
                 {:.1} ms; classical, {} waves: {:.1} ms ({:.1} reciprocal), energy {:.1e} of \
                 E_recip, forces {:.1e} RMS",
                cell.lengths()[0] / ANGSTROM,
                p.alpha * ANGSTROM,
                1e3 * t_real,
                classical.wave_vectors(),
                1e3 * t_classical,
                1e3 * (t_classical - t_real),
                (c.energy.reciprocal - r).abs() / r.abs(),
                rms(&c.forces, &reference.forces)
            );
            let chosen = PmeParameters::for_accuracy(&cell, &p, accuracy, n);
            for order in 4..=8 {
                let grid = PmeParameters::grid_for(&cell, &p, accuracy, order);
                let meshed = classical.clone().with_mesh(mesh(order, grid));
                let t = best(&mut || {
                    std::hint::black_box(meshed.evaluate(&q, &at, &ex));
                });
                let e = meshed.evaluate(&q, &at, &ex);
                println!(
                    "    order {order}, {grid:?}{}: {:.1} ms ({:.2} reciprocal), energy {:.1e} \
                     of E_recip, forces {:.1e} RMS",
                    if chosen == mesh(order, grid) {
                        " (chosen)"
                    } else {
                        ""
                    },
                    1e3 * t,
                    1e3 * (t - t_real),
                    (e.energy.reciprocal - r).abs() / r.abs(),
                    rms(&e.forces, &reference.forces)
                );
            }
        }
    }
}
