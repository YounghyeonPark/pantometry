//! Smooth particle-mesh Ewald: the reciprocal part of [`crate::ewald`]'s sum by cardinal B-splines
//! on a grid and a fast Fourier transform, in `O(N + K log K)` where the classical sum is
//! `O(N × waves)`.
//!
//! The method is Essmann, Perera, Berkowitz, Darden, Lee and Pedersen, *J. Chem. Phys.* **103**,
//! 8577 (1995), which smooths the Lagrange interpolation of Darden, York and Pedersen, *J. Chem.
//! Phys.* **98**, 10089 (1993). Every formula below is derived here from those papers' idea — the
//! structure factor's phases interpolated by B-splines, the result a convolution — and no
//! library's code was read.
//!
//! # What it computes
//!
//! [`crate::ewald`]'s reciprocal energy, written over all integer vectors `m` with
//! `m̃ = (m₁/L₁, m₂/L₂, m₃/L₃)` (so `k = 2π m̃`):
//!
//! ```text
//! E_recip = (1/2πV) Σ_{m ≠ 0} w(m) |S(m)|²,   w(m) = exp(−π² m̃²/α²) / m̃²,   S(m) = Σ_j q_j exp(2πi m·u_j/K)
//! ```
//!
//! with `u_{a j} = K_a r_{a j} / L_a` the scaled fractional coordinate on a grid of `K_a` points
//! along axis `a`. It is the classical `(2π/V) Σ_k exp(−k²/4α²)/k² |S(k)|²` with `k = 2πm̃`.
//!
//! # The interpolation
//!
//! **The cardinal B-spline** `M_n` of order `n` is the box `M_1` on `[0, 1)` convolved with itself
//! `n − 1` times: `M_2(u) = 1 − |u − 1|` on `[0, 2]`, and in closed form
//! `M_n(u) = Σ_{k=0}^{n} (−1)^k C(n, k) (u − k)₊^{n−1} / (n − 1)!`, the `n`-th difference of a
//! truncated power. It obeys
//!
//! ```text
//! M_n(u) = [u M_{n−1}(u) + (n − u) M_{n−1}(u − 1)] / (n − 1),      M_n′(u) = M_{n−1}(u) − M_{n−1}(u − 1)
//! ```
//!
//! — the second because `d/du (f ∗ M_1)(u) = f(u) − f(u − 1)`; the first from the closed form,
//! writing `(u − k)₊^{n−1} = (u − k)(u − k)₊^{n−2}` and splitting `C(n, k)` by Pascal's rule into the
//! two sums that are `M_{n−1}(u)` and `M_{n−1}(u − 1)`. `M_n` is supported on `[0, n]`, is piecewise
//! a polynomial of degree `n − 1`, is `n − 2` times continuously differentiable, and its integer
//! translates sum to one. [`b_spline`] evaluates the `n` translates that are non-zero at a point by
//! the recursion; the tests hold it to the closed form.
//!
//! **The Euler exponential spline.** Ask for `exp(2πi m u/K) ≈ b(m) Σ_k M_n(u − k) exp(2πi m k/K)`
//! to be exact at every integer `u`. Putting `k = u − s` makes the sum
//! `exp(2πi m u/K) Σ_s M_n(s) exp(−2πi m s/K)`, so
//!
//! ```text
//! 1 / b(m) = Σ_{s=1}^{n−1} M_n(s) exp(−2πi m s/K),          |b(m)|² = 1 / |Σ_{s=1}^{n−1} M_n(s) exp(2πi m s/K)|²
//! ```
//!
//! ([`euler_denominator`] is the modulus squared of that sum). It is Essmann et al.'s
//! `exp(2πi(n−1)m/K) / Σ_{k=0}^{n−2} M_n(k+1) exp(2πi m k/K)`, the two related by `M_n(s) = M_n(n − s)`.
//! Only `|b|²` enters the energy.
//!
//! **How good it is, in closed form.** Poisson's summation turns the spline sum into
//! `exp(2πi m u/K) Σ_p e^{−2πi p u} M̂_n(θ − 2πp)` with `θ = 2πm/K` and `M̂_n(ω) = e^{−iωn/2}
//! (sin(ω/2)/(ω/2))ⁿ`, and the ratio of each alias to the main term is real:
//! `r_p = (x/(x − p))ⁿ` with `x = m/K`. So the interpolated phase of one charge is exactly
//!
//! ```text
//! b(m) Σ_k M_n(u − k) e^{2πimk/K} = e^{2πimu/K} (1 + η(m, u)),    1 + η = (1 + Σ_{p≠0} r_p e^{−2πipu}) / (1 + Σ_{p≠0} r_p)
//! ```
//!
//! which is zero error at every grid point and `O((m/K)ⁿ)` between them. In three dimensions the
//! factors multiply. `tests/the_particle_mesh_against_closed_forms.rs` sums the energy that way,
//! with no grid, no spline and no transform, and holds the mesh to it.
//!
//! **For odd `n` the leading error is imaginary.** `r₁ = −xⁿ(1 + nx + …)` and `r₋₁ = xⁿ(1 − nx + …)`
//! when `n` is odd, so `R = Σ_{p≠0} r_p` is `O(x^{n+1})` and `η ≈ 2i xⁿ sin 2πu`: an error of phase
//! rather than of amplitude. **Where every charge has the same offset in its grid cell** — each
//! charge's own term, whose mean is [`Pme::self_bias`], and a crystal whose ions are all a whole
//! number of grid spacings apart — the error is the common factor `|Π(1 + η_a)|² − 1`, whose
//! leading terms are `|η|²` and `Re η`, and it falls as `K^{−2⌈n/2⌉}`: the same rate for orders 3
//! and 4, and for 5 and 6, which the tests measure on NaCl. **That is all the pairing covers.**
//! Charges at several offsets have cross terms `η_i* + η_j`, which a phase error enters at first
//! order: a crystal of four ions at generic fractions measures 3.6–4.5 for order 4 and 4.5–6.2 for
//! order 5 over doublings to 128³ (an independent implementation, about `n`), and the spread of
//! disordered charges, about `n + 1` for odd orders from 16³ to 32³ where the bias still leads, falls
//! to 2.2–2.7, 4.5–4.7 and 6.2–6.8 for orders 3, 5 and 7 by 128³ — between `n − 1` and `n`, not
//! settled. For even orders the rate is `n` throughout, and the forces' is `n − 1` at every order.
//!
//! # The mesh
//!
//! Spreading every charge onto the grid, `Q(k) = Σ_j q_j Π_a M_n(u_{aj} − k_a)` (periodically in
//! each `k_a`), makes the interpolated structure factor `b₁b₂b₃ F(Q)(−m)`, with `F` the discrete
//! transform `F(Q)(m) = Σ_k Q(k) e^{−2πi m·k/K}`, and so
//!
//! ```text
//! E_recip = ½ Σ_m G(m) |F(Q)(m)|²,      G(m) = exp(−π² m̃²/α²) / (π V m̃²) × |b₁(m₁)|² |b₂(m₂)|² |b₃(m₃)|²
//! ```
//!
//! over the grid's `K₁K₂K₃` vectors, `m_a` taken in `(−K_a/2, K_a/2)`. `G` is the influence
//! function, built once per box. Then:
//!
//! - **The potential on the grid** is `φ(k) = ∂E/∂Q(k) = Σ_m G(m) Re[F(Q)(m) e^{2πi m·k/K}]`, the
//!   unnormalised inverse transform of `G F(Q)`, real because `G(m) = G(−m)`.
//! - **The force** on atom `i` is `−Σ_k ∂Q(k)/∂r_i φ(k)`, with `∂/∂r_a M_n(u_a − k_a) =
//!   (K_a/L_a) M_n′(u_a − k_a)`: the analytic derivative of the spline, so the forces are the exact
//!   gradient of this energy. They are not the interpolated classical force, and they do not sum to
//!   zero: see below.
//! - **The virial** `W_ab = −∂E/∂ε_ab`. Under a homogeneous strain every `u` is fixed, so `Q` and
//!   `F(Q)` are, and only `V` and `m̃` move: each `m`'s energy `E_m = ½ G|F|²` contributes
//!   `E_m [δ_ab − 2 m̃_a m̃_b (1/m̃² + π²/α²)]`, the classical stress at `k = 2πm̃`. Exact for this
//!   energy, since `b(m)` depends on integers alone.
//!
//! # The planes at `m_a = K_a/2`, and the zero of `b(m)` there
//!
//! **For odd `n`, `1/b(K/2)` is exactly zero**: `Σ_s M_n(s)(−1)^s` pairs `s` with `n − s`, whose
//! signs are opposite when `n` is odd and whose values are equal. So `|b(K/2)|²` is infinite. For
//! even `n` it is finite, `1/T²` with `T` the coefficient of `x^{n−1}` in `tan x` — `1/3`, `2/15`,
//! `17/315` for `n` = 4, 6, 8 — which the tests check. Even then the Nyquist vector has no sign:
//! `K/2` and `−K/2` are the same grid point and different vectors under a shear.
//!
//! **Every vector with some `m_a = K_a/2` is left out**, for every order: `G` is zero there. Its
//! weight is at most `exp(−π² (K/2L)²/α²)`, and there the interpolation is no interpolation at all
//! — at `x = 1/2` the first alias ratio `r₁` is `±1`, as large as the term — so what is lost is of
//! the order of the error the grid already has there, and it is counted: [`Pme::self_bias`] has it
//! among the vectors left out. A crystal with every ion on a grid point, where the spline is exact
//! and the Nyquist planes are the only loss on the grid, gives Madelung's constants to 6.7e-16
//! (NaCl) and 1.8e-14 (CsCl). A
//! neighbour's average in its place, the other choice, would interpolate a pole for odd `n`.
//!
//! # The transform, and why it is here rather than the kernel's
//!
//! A three-dimensional complex transform, radix two, written in this module with **its twiddle
//! factors from this module's own `sin` and `cos`**: Taylor's series on `[0, π/4]` and the exact
//! symmetries of the circle, `+ − × ÷` only. [`pantometry_core::transform`] computes its twiddles
//! from the platform's `cos` and `sin`, at every butterfly, so a mesh built on it would differ
//! across platforms in its last bits; and it is one- and two-dimensional, and costs a `sin_cos` per
//! butterfly, 2.4 million per transform of a 64³ grid. So **the mesh is the same bits on every
//! platform**, as the classical sum is not (its phases start from the platform's `sin_cos`): the
//! splines are arithmetic, `G` uses the kernel's [`exp`], and the transform is arithmetic. Lengths
//! are powers of two.
//!
//! **Decoupling.** The energy is a quadratic form in `Q`, which is linear in the charges, so with
//! the charges split into a rest and a group, `E(all) − E(rest) − E(group) = Σ_m G Re(F_rest* F_group)`
//! exactly. The two grids are transformed as one complex grid — the rest in the real part, the
//! group in the imaginary — and separated by `F_rest(m) = [C(m) + C(−m)*]/2`,
//! `F_group(m) = [C(m) − C(−m)*]/2i`: the cross term, the rest's own energy, and with one inverse
//! transform of `G(F_rest + λF_group) + i λ G F_rest` the potentials both sets of atoms feel, cost
//! what one evaluation costs. `E(group)` is never computed: [`crate::PeriodicDecoupling`] puts the
//! group's own pairs in vacuum.
//!
//! # Choosing the grid and the order
//!
//! [`PmeParameters::for_accuracy`] takes the classical sum's α, its cutoff and the same accuracy
//! `δ` in units of `k_e Q / r_c`, and for each order from 4 to 8 doubles the axis with the coarsest
//! spacing until the estimated error of uncorrelated charges, **the spread and the bias together**,
//! `(σ² + b²)^½` from [`Pme::energy_error`] and [`Pme::self_bias`], is at most `δ`; of those, it
//! takes the one [`PmeParameters::cost`] puts cheapest for the number of atoms. **It delivers
//! `δ`**: on sixteen disordered boxes the RMS error at its choice is 0.99 and 0.96 of `δ` at 10⁻⁴
//! and 10⁻⁵, and on lattice water 0.04–0.24 of it. A first version held σ alone to `δ`, which, the
//! bias being ten times σ on disordered charges, did not. **The estimates are derived here**: for
//! uncorrelated charges the interpolation error of a pair's term at `m` is `η_i* + η_j`, whose
//! cross-pair terms have random phases, so
//!
//! ```text
//! σ_E ≈ k_e (Q / π V) [Σ_m w(m)² (⟨|η(m)|²⟩ + ⟨η(m)⟩²)]^½
//! ```
//!
//! with the averages over a charge's position in closed form per axis: `⟨1 + η_a⟩ = 1/(1 + R_a)` and
//! `⟨|1 + η_a|²⟩ = (1 + Σ_p r_p²)/(1 + R_a)²`, `R_a = Σ_{p≠0} r_p`, the alias sums truncated at
//! `|p| = 64` with the rest of `R` by its integral, and the products over the axes expanded so that
//! nothing cancels. Each charge's own term, `|1 + η|²`, has a mean as well, [`Pme::self_bias`]:
//! for uncorrelated charges the energy is off by that per `Σq²`, which, like the classical sum's
//! [`EwaldParameters::reciprocal_bias`], has no cross term between two sets of charges. **It is a
//! constant only on average**: each charge's own term varies with its place in its grid cell, and
//! that variation is what makes the energy depend on where the system is and the net force not
//! zero. The estimate leaves the variation out.
//!
//! **Both are for uncorrelated charges, and molecules are not.** In water each charge has two of
//! the opposite sign a bond away, and their pair terms cancel most of each one's own: measured on
//! 64 lattice waters over orders 4 to 8 at the grids `grid_for` gives, the mesh's error is −0.33 to
//! 0.40 of [`Pme::self_bias`] and, being correlated, up to 3.9 σ. The target, which counts the
//! whole bias, is conservative there. The classical sum's bias estimate fails on water the same
//! way, and the spread is an estimate of the same kind as Kolafa and Perram's.
//!
//! # Measured
//!
//! `tests/the_particle_mesh_against_closed_forms.rs`:
//!
//! - **The rates.** NaCl's energy, every ion at one offset, against the converged classical sum from
//!   32³ to 64³: 3.86, 3.89, 5.88, 5.94 and 8.00 for orders 3, 4, 5, 6 and 8 — `2⌈n/2⌉`; CsCl's at even
//!   orders, 4.08, 6.19, 8.31 — `n`. The forces on disordered charges, RMS: 1.97, 3.10, 4.03, 5.16 and
//!   7.36 for orders 3, 4, 5, 6 and 8 from 16³ to 32³, and `n − 1` on every doubling to 128³ — since
//!   differentiating the alias terms brings down `p/x`.
//! - **The estimates**, over 48 disordered boxes of 96 charges in 18 × 15.5 × 21 Å: the mean error
//!   0.986–1.023 of [`Pme::self_bias`], the RMS about it 0.886–1.107 of [`Pme::energy_error`], their
//!   geometric mean over six cases 0.973.
//! - **The cost**, release, one core, `x86_64-pc-windows-gnu`, TIP3P lattice boxes at `r_c` = 9 Å,
//!   the reciprocal part alone, at the order and grid [`PmeParameters::for_accuracy`] chose, with
//!   each sum's error against the classical sum at `k_c = 9α`, from one run on a loaded machine
//!   (real space 31 ms at 3 000 atoms and 234–257 ms at 24 000, against 26 and 195 unloaded):
//!
//!   | atoms | δ | classical: waves, ms, energy, RMS force | mesh: order, grid, ms, energy, RMS force |
//!   | --- | --- | --- | --- |
//!   | 3 000 | 10⁻⁵ | 1 051, 17.8, 2.6e-3, 5.2e-4 | 6, 32³, 4.9, 3.3e-5, 1.0e-5 |
//!   | 3 000 | 10⁻⁶ | 2 192, 35.0, 3.1e-4, 9.2e-5 | 8, 32³, 5.6, 1.3e-5, 2.7e-6 |
//!   | 24 000 | 10⁻⁵ | 5 537, 681, 5.2e-3, 1.3e-3 | 6, 64³, 33.5, 1.6e-5, 6.1e-6 |
//!   | 24 000 | 10⁻⁶ | 12 838, 1 575, 1.6e-3, 2.3e-4 | 8, 64³, 22.9, 5.9e-6, 1.3e-6 |
//!
//!   The energies relative to the reciprocal energy, the forces to the RMS force. Real space and the
//!   exclusions are the same for both and now nearly all of an evaluation. At every setting the
//!   chosen mesh is more accurate than the classical sum at the same `δ`, in energy by 24 to 330
//!   times and in forces by 34 to 210, and at 24 000 atoms 20 to 70 times cheaper.
//!
//! //! # Not here
//!
//! A triclinic box; a real-to-complex transform (the mesh transforms a complex grid with a zero
//! imaginary part, twice the arithmetic a real transform needs, except when decoupling, where the
//! imaginary part carries the group); threads; any grid length but a power of two; and a force
//! target for the grid, which is chosen from the energy's estimates. The transform is this
//! module's own, private: moving it to [`pantometry_core::transform`] would change the kernel's
//! public API, and is a separate decision.

// The splines, grids and axes are written as the formulas' index loops, which read against them.
#![allow(clippy::needless_range_loop)]

use crate::ewald::{EwaldParameters, COULOMB};
use crate::periodic::PeriodicBox;
use pantometry_core::math::exp;
use std::f64::consts::PI;

/// The lowest spline order a mesh accepts: `M_3`'s derivative is continuous, `M_2`'s is not.
pub const MIN_ORDER: usize = 3;

/// The highest spline order a mesh accepts.
pub const MAX_ORDER: usize = 12;

/// The longest grid edge, points: 4 096.
pub const MAX_GRID: usize = 1 << 12;

/// The spline order and the grid of a smooth particle-mesh Ewald sum.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PmeParameters {
    /// `n`, the order of the cardinal B-spline: a piecewise polynomial of degree `n − 1` over `n`
    /// grid intervals. From [`MIN_ORDER`] to [`MAX_ORDER`].
    pub order: usize,
    /// `K_a`, the grid points along each axis, each a power of two, at least the order and at most
    /// [`MAX_GRID`].
    pub grid: [usize; 3],
}

impl PmeParameters {
    /// The order and grid for an Ewald sum with `ewald`'s α and cutoff in `cell`, so that the
    /// estimated error of uncorrelated charges' reciprocal energy, **the spread and the bias
    /// together**, `(σ² + b²)^½` from [`Pme::energy_error`] and [`Pme::self_bias`], is at most
    /// `accuracy` in units of `k_e Q / r_c` — the classical [`EwaldParameters::for_accuracy`]'s
    /// units — and, among the orders 4 to 8 that reach it, the one [`PmeParameters::cost`] puts
    /// cheapest for `atoms` atoms. See the module documentation.
    ///
    /// **Unlike the classical sum's, the target includes the bias.** The classical sum's bias is
    /// a constant at fixed charges; the mesh's is a constant only on average, each charge's own
    /// term varying with its place in its cell, and it is the larger part — ten times the spread on
    /// disordered charges. For a molecule's charges, which cancel most of it, the target is
    /// conservative.
    ///
    /// # Panics
    ///
    /// If `accuracy` is not in `(0, 0.01]`, or no grid of at most [`MAX_GRID`] points an edge
    /// reaches it.
    pub fn for_accuracy(
        cell: &PeriodicBox,
        ewald: &EwaldParameters,
        accuracy: f64,
        atoms: usize,
    ) -> PmeParameters {
        assert!(
            accuracy > 0.0 && accuracy <= 0.01,
            "a mesh accuracy must be in (0, 0.01]"
        );
        let mut best: Option<(f64, PmeParameters)> = None;
        for order in 4..=8 {
            let candidate = PmeParameters {
                order,
                grid: PmeParameters::grid_for(cell, ewald, accuracy, order),
            };
            let cost = candidate.cost(atoms);
            if !matches!(best, Some((c, _)) if c <= cost) {
                best = Some((cost, candidate));
            }
        }
        best.expect("five orders were tried").1
    }

    /// The cost model [`PmeParameters::for_accuracy`] minimises, for `atoms` atoms:
    /// `atoms × n³ + 2 K log₂ K`, `K = K₁K₂K₃` — the spline's `n³` points spread and gathered per
    /// atom, and two complex transforms, weighted as they were measured on one core, release,
    /// `x86_64-pc-windows-gnu`: about 1.3 ns per atom and `n³`, and 2.7 ns per `K log₂ K`. A first
    /// model, `4 n³` against `3 K log K`, put the spline three times too dear. It is a model, and
    /// `tests/the_particle_mesh_against_closed_forms.rs`'s `the_cost` times what it chose against
    /// the other orders.
    pub fn cost(&self, atoms: usize) -> f64 {
        let [k1, k2, k3] = self.grid;
        let points = (k1 * k2 * k3) as f64;
        let n3 = (self.order * self.order * self.order) as f64;
        let log2 = (k1.trailing_zeros() + k2.trailing_zeros() + k3.trailing_zeros()) as f64;
        atoms as f64 * n3 + 2.0 * points * log2
    }

    /// The smallest grid of order `order` whose estimate `(σ² + b²)^½` reaches `accuracy` (see
    /// [`PmeParameters::for_accuracy`]): from the order's next power of two on every axis, the axis
    /// with the coarsest spacing doubled until it does.
    ///
    /// # Panics
    ///
    /// If `order` is outside [`MIN_ORDER`]..=[`MAX_ORDER`], or no grid of at most [`MAX_GRID`]
    /// points an edge reaches `accuracy`.
    pub fn grid_for(
        cell: &PeriodicBox,
        ewald: &EwaldParameters,
        accuracy: f64,
        order: usize,
    ) -> [usize; 3] {
        assert!(
            (MIN_ORDER..=MAX_ORDER).contains(&order),
            "a spline order must be from {MIN_ORDER} to {MAX_ORDER}"
        );
        let l = cell.lengths();
        let mut grid = [order.next_power_of_two(); 3];
        loop {
            let p = PmeParameters { order, grid };
            let spread = spread_coefficient(cell, ewald.alpha, &p);
            let bias = bias_coefficient(cell, ewald.alpha, &p);
            let e = (spread * spread + bias * bias).sqrt() * ewald.cutoff;
            if e <= accuracy {
                return grid;
            }
            let mut coarsest = 0;
            for a in 1..3 {
                if l[a] / grid[a] as f64 > l[coarsest] / grid[coarsest] as f64 {
                    coarsest = a;
                }
            }
            assert!(
                grid[coarsest] < MAX_GRID,
                "no grid of at most {MAX_GRID} points an edge reaches {accuracy:e} at order {order}"
            );
            grid[coarsest] *= 2;
        }
    }
}

/// `M_n(t + j)` and `M_n′(t + j)` for `j = 0, …, n − 1`, the cardinal B-spline of order `order` and
/// its derivative at the `n` points of the form `t + j` inside its support: the weights a charge at
/// fractional offset `t ∈ [0, 1)` past grid point `k₀` gives the points `k₀ − j`. By the recursion
/// of the module documentation, the derivative from order `n − 1`.
///
/// # Panics
///
/// If `order` is outside [`MIN_ORDER`]..=[`MAX_ORDER`].
pub fn b_spline(order: usize, t: f64) -> (Vec<f64>, Vec<f64>) {
    assert!(
        (MIN_ORDER..=MAX_ORDER).contains(&order),
        "a spline order must be from {MIN_ORDER} to {MAX_ORDER}"
    );
    let mut w = [0.0; MAX_ORDER];
    let mut dw = [0.0; MAX_ORDER];
    weights(order, t, &mut w, &mut dw);
    (w[..order].to_vec(), dw[..order].to_vec())
}

/// [`b_spline`] into fixed arrays.
fn weights(order: usize, t: f64, w: &mut [f64; MAX_ORDER], dw: &mut [f64; MAX_ORDER]) {
    // Order 2: M_2(t) = t, M_2(t + 1) = 1 − t.
    let mut m = [0.0; MAX_ORDER];
    m[0] = t;
    m[1] = 1.0 - t;
    for p in 3..=order {
        if p == order {
            // M_p′(t + j) = M_{p−1}(t + j) − M_{p−1}(t + j − 1), from the order below.
            dw[0] = m[0];
            for j in 1..p - 1 {
                dw[j] = m[j] - m[j - 1];
            }
            dw[p - 1] = -m[p - 2];
        }
        // M_p(t + j) = [(t + j) M_{p−1}(t + j) + (p − t − j) M_{p−1}(t + j − 1)] / (p − 1), from
        // the top down so that each M_{p−1} is read before it is overwritten.
        let d = (p - 1) as f64;
        m[p - 1] = (1.0 - t) * m[p - 2] / d;
        for j in (1..p - 1).rev() {
            m[j] = ((t + j as f64) * m[j] + ((p - j) as f64 - t) * m[j - 1]) / d;
        }
        m[0] = t * m[0] / d;
    }
    w[..order].copy_from_slice(&m[..order]);
}

/// `|Σ_{s=1}^{n−1} M_n(s) exp(2πi m s/K)|²` for order `order` on a grid of `k` points: one over the
/// Euler exponential spline's `|b(m)|²`. Zero at `m = K/2` for odd orders; see the module
/// documentation.
///
/// # Panics
///
/// If `order` is outside [`MIN_ORDER`]..=[`MAX_ORDER`] or `k` is not a power of two.
pub fn euler_denominator(order: usize, m: usize, k: usize) -> f64 {
    assert!(k.is_power_of_two(), "a grid edge must be a power of two");
    let (w, _) = b_spline(order, 0.0);
    let (mut re, mut im) = (0.0, 0.0);
    for (s, &ms) in w.iter().enumerate().skip(1) {
        let (c, si) = unit_root((m % k) * s % k, k);
        re += ms * c;
        im += ms * si;
    }
    re * re + im * im
}

/// `cos x` and `sin x` for `x ∈ [0, π/4]` by Taylor's series, eleven terms each: the last omitted
/// term is below `(π/4)²³/23!`, `4e-25`.
fn octant(x: f64) -> (f64, f64) {
    let z = x * x;
    let (mut c, mut s) = (0.0, 0.0);
    // Horner in z from the top: cos = Σ (−1)^k z^k/(2k)!, sin = x Σ (−1)^k z^k/(2k+1)!, each
    // factorial exact (to 22!, each one's odd part is below 2⁵³), each coefficient one correctly
    // rounded `÷`.
    for k in (0..=10u32).rev() {
        let even = factorial(2 * k);
        let odd = even * f64::from(2 * k + 1);
        let sign = if k & 1 == 0 { 1.0 } else { -1.0 };
        c = c * z + sign / even;
        s = s * z + sign / odd;
    }
    (c, x * s)
}

/// `n!` as a double, exact to `22!`.
fn factorial(n: u32) -> f64 {
    let mut f = 1.0;
    for i in 2..=n {
        f *= f64::from(i);
    }
    f
}

/// `cos(2πj/n)` and `sin(2πj/n)` for `n` a power of two: the angle reduced exactly to `[0, π/4]`
/// by the quadrant and the octant's reflection, and [`octant`] there. The same bits on every
/// platform; `cos²+ sin² = 1` to two ulps.
fn unit_root(j: usize, n: usize) -> (f64, f64) {
    debug_assert!(n.is_power_of_two());
    let j = j % n;
    match n {
        1 => return (1.0, 0.0),
        2 => return if j == 0 { (1.0, 0.0) } else { (-1.0, 0.0) },
        _ => {}
    }
    let q = n / 4;
    let (quadrant, r) = (j / q, j % q);
    // (c, s) of 2πr/n with 0 ≤ r < n/4: past the octant, cos(π/2 − y) = sin y.
    let (c, s) = if 2 * r <= q {
        octant(std::f64::consts::TAU / n as f64 * r as f64)
    } else {
        let (c, s) = octant(std::f64::consts::TAU / n as f64 * (q - r) as f64);
        (s, c)
    };
    match quadrant {
        0 => (c, s),
        1 => (-s, c),
        2 => (-c, -s),
        _ => (s, -c),
    }
}

/// One axis of the transform: its length and `exp(−2πi j/n)` for `j < n/2`.
#[derive(Clone, Debug, PartialEq)]
struct Line {
    n: usize,
    re: Vec<f64>,
    im: Vec<f64>,
}

impl Line {
    fn new(n: usize) -> Line {
        let (re, im) = (0..n / 2)
            .map(|j| {
                let (c, s) = unit_root(j, n);
                (c, -s)
            })
            .unzip();
        Line { n, re, im }
    }

    /// The transform of one line in place, `Σ_k x_k exp(∓2πi jk/n)`, unnormalised either way:
    /// `inverse` takes the plus sign.
    fn transform(&self, re: &mut [f64], im: &mut [f64], inverse: bool) {
        let n = self.n;
        if n < 2 {
            return;
        }
        let mut j = 0usize;
        for i in 1..n {
            let mut bit = n >> 1;
            while j & bit != 0 {
                j ^= bit;
                bit >>= 1;
            }
            j |= bit;
            if i < j {
                re.swap(i, j);
                im.swap(i, j);
            }
        }
        let sign = if inverse { -1.0 } else { 1.0 };
        let mut len = 2;
        while len <= n {
            let half = len / 2;
            let step = n / len;
            let mut base = 0;
            while base < n {
                for k in 0..half {
                    let (wr, wi) = (self.re[k * step], sign * self.im[k * step]);
                    let (a, b) = (base + k, base + k + half);
                    let vr = re[b] * wr - im[b] * wi;
                    let vi = re[b] * wi + im[b] * wr;
                    let (ur, ui) = (re[a], im[a]);
                    re[a] = ur + vr;
                    im[a] = ui + vi;
                    re[b] = ur - vr;
                    im[b] = ui - vi;
                }
                base += len;
            }
            len <<= 1;
        }
    }
}

/// The three-dimensional transform of a grid stored `[k₁][k₂][k₃]`, `k₃` fastest: the third axis,
/// then the second, then the first, in a fixed order.
fn transform_3d(lines: &[Line; 3], re: &mut [f64], im: &mut [f64], inverse: bool) {
    let [n1, n2, n3] = [lines[0].n, lines[1].n, lines[2].n];
    for row in 0..n1 * n2 {
        let s = row * n3..(row + 1) * n3;
        lines[2].transform(&mut re[s.clone()], &mut im[s], inverse);
    }
    let longest = n1.max(n2);
    let (mut br, mut bi) = (vec![0.0; longest], vec![0.0; longest]);
    for i1 in 0..n1 {
        for i3 in 0..n3 {
            for i2 in 0..n2 {
                let t = (i1 * n2 + i2) * n3 + i3;
                br[i2] = re[t];
                bi[i2] = im[t];
            }
            lines[1].transform(&mut br[..n2], &mut bi[..n2], inverse);
            for i2 in 0..n2 {
                let t = (i1 * n2 + i2) * n3 + i3;
                re[t] = br[i2];
                im[t] = bi[i2];
            }
        }
    }
    let plane = n2 * n3;
    for p in 0..plane {
        for i1 in 0..n1 {
            br[i1] = re[i1 * plane + p];
            bi[i1] = im[i1 * plane + p];
        }
        lines[0].transform(&mut br[..n1], &mut bi[..n1], inverse);
        for i1 in 0..n1 {
            re[i1 * plane + p] = br[i1];
            im[i1 * plane + p] = bi[i1];
        }
    }
}

/// The signed vector `m` a grid index `j` on an axis of `k` points stands for, `(−K/2, K/2)`, or
/// `None` on the Nyquist plane `j = K/2`.
fn signed(j: usize, k: usize) -> Option<f64> {
    if 2 * j == k {
        None
    } else if 2 * j < k {
        Some(j as f64)
    } else {
        Some(j as f64 - k as f64)
    }
}

/// A smooth particle-mesh Ewald sum's reciprocal part for one box, one α and one
/// [`PmeParameters`]: the influence function, built once. See the module documentation.
#[derive(Clone, Debug, PartialEq)]
pub struct Pme {
    cell: PeriodicBox,
    alpha: f64,
    parameters: PmeParameters,
    lines: [Line; 3],
    /// `|b_a(m)|²` per axis by grid index, zero on the Nyquist plane.
    b_squared: [Vec<f64>; 3],
    /// `G(m)` by grid index, zero at `m = 0` and wherever an `m_a = K_a/2`.
    influence: Vec<f64>,
}

impl Pme {
    /// The mesh for `cell` and splitting parameter `alpha` (per metre).
    ///
    /// # Panics
    ///
    /// If the order is outside [`MIN_ORDER`]..=[`MAX_ORDER`], a grid edge is not a power of two
    /// from the order to [`MAX_GRID`], or α is not positive.
    pub fn new(cell: PeriodicBox, alpha: f64, parameters: PmeParameters) -> Pme {
        let n = parameters.order;
        assert!(
            (MIN_ORDER..=MAX_ORDER).contains(&n),
            "a spline order must be from {MIN_ORDER} to {MAX_ORDER}"
        );
        assert!(
            parameters
                .grid
                .iter()
                .all(|&k| k.is_power_of_two() && k >= n && k <= MAX_GRID),
            "each grid edge must be a power of two from the order to {MAX_GRID}"
        );
        assert!(alpha > 0.0, "α must be positive");
        let grid = parameters.grid;
        let lines = grid.map(Line::new);
        let b_squared = grid.map(|k| {
            (0..k)
                .map(|j| {
                    if 2 * j == k {
                        0.0
                    } else {
                        1.0 / euler_denominator(n, j, k)
                    }
                })
                .collect::<Vec<f64>>()
        });
        let l = cell.lengths();
        let v = cell.volume();
        let a2 = alpha * alpha;
        let mut influence = vec![0.0; grid[0] * grid[1] * grid[2]];
        for j1 in 0..grid[0] {
            let Some(m1) = signed(j1, grid[0]) else {
                continue;
            };
            let x1 = m1 / l[0];
            for j2 in 0..grid[1] {
                let Some(m2) = signed(j2, grid[1]) else {
                    continue;
                };
                let x2 = m2 / l[1];
                let b12 = b_squared[0][j1] * b_squared[1][j2];
                for j3 in 0..grid[2] {
                    let Some(m3) = signed(j3, grid[2]) else {
                        continue;
                    };
                    if j1 == 0 && j2 == 0 && j3 == 0 {
                        continue;
                    }
                    let x3 = m3 / l[2];
                    let mm = x1 * x1 + x2 * x2 + x3 * x3;
                    influence[(j1 * grid[1] + j2) * grid[2] + j3] =
                        exp(-PI * PI * mm / a2) / (PI * v * mm) * b12 * b_squared[2][j3];
                }
            }
        }
        Pme {
            cell,
            alpha,
            parameters,
            lines,
            b_squared,
            influence,
        }
    }

    /// The same mesh — the same order, grid and α — in another box: what a derivative with
    /// respect to the box needs, so that no grid point moves relative to the charges.
    pub fn with_cell(&self, cell: PeriodicBox) -> Pme {
        Pme::new(cell, self.alpha, self.parameters)
    }

    /// The order and grid.
    pub fn parameters(&self) -> PmeParameters {
        self.parameters
    }

    /// The box.
    pub fn cell(&self) -> PeriodicBox {
        self.cell
    }

    /// α, per metre.
    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    /// `|b(m)|²` on axis `axis` at grid index `j` (`m = j`, or `j − K` past `K/2`), as the
    /// influence function uses it: zero at `j = K/2`, whose plane is left out.
    ///
    /// # Panics
    ///
    /// If `axis` is past 2 or `j` past the grid.
    pub fn b_squared(&self, axis: usize, j: usize) -> f64 {
        self.b_squared[axis][j]
    }

    /// The estimated RMS error of the reciprocal energy, joules, for charges whose squares sum to
    /// `sum_q2` (e²), from the interpolation's alias terms over uncorrelated charges: see the
    /// module documentation. It leaves out [`Pme::self_bias`], the mean.
    pub fn energy_error(&self, sum_q2: f64) -> f64 {
        COULOMB * sum_q2 * spread_coefficient(&self.cell, self.alpha, &self.parameters)
    }

    /// The mean error of the reciprocal energy, joules, for charges whose squares sum to `sum_q2`
    /// at uniformly distributed positions: each charge's own term `q² Σ_m w(m) (⟨|1 + η|²⟩ − 1) /
    /// 2πV` over the grid's vectors, less its terms at every vector the grid leaves out (beyond
    /// it, and its Nyquist planes). A constant at fixed charges on average only — each charge's own
    /// term varies with its place in its cell — and for uncorrelated charges only: see the module
    /// documentation.
    pub fn self_bias(&self, sum_q2: f64) -> f64 {
        COULOMB * sum_q2 * bias_coefficient(&self.cell, self.alpha, &self.parameters)
    }

    /// The reciprocal sum, with a partition, as [`crate::Ewald`]'s classical one: `group` marks the
    /// atoms whose charges are split off and `lambda` scales their cross terms with the rest.
    /// Returns `(E_rest, E_cross)` and adds to `forces` and `virial` those of `E_rest + λ E_cross`.
    pub(crate) fn reciprocal(
        &self,
        charges: &[f64],
        group: Option<(&[bool], f64)>,
        at: &[[f64; 3]],
        forces: Option<&mut [[f64; 3]]>,
        virial: &mut [[f64; 3]; 3],
    ) -> (f64, f64) {
        let n = self.parameters.order;
        let grid = self.parameters.grid;
        let [k1, k2, k3] = grid;
        let total = k1 * k2 * k3;
        let l = self.cell.lengths();
        let marked = |i: usize| group.is_some_and(|(g, _)| g[i]);
        let lambda = group.map_or(0.0, |(_, x)| x);

        // Each charge's spline, per axis: the grid points k₀ − j and the weights there.
        let spline = |p: &[f64; 3]| {
            let mut index = [[0usize; MAX_ORDER]; 3];
            let mut value = [[0.0; MAX_ORDER]; 3];
            let mut slope = [[0.0; MAX_ORDER]; 3];
            for a in 0..3 {
                let f = p[a] / l[a];
                let u = grid[a] as f64 * (f - f.floor());
                let base = u.floor();
                weights(n, u - base, &mut value[a], &mut slope[a]);
                let k0 = base as i64;
                for j in 0..n {
                    index[a][j] = (k0 - j as i64).rem_euclid(grid[a] as i64) as usize;
                }
            }
            (index, value, slope)
        };

        // Spread: the rest in the real part, the group in the imaginary.
        let mut re = vec![0.0; total];
        let mut im = vec![0.0; total];
        for (i, &q) in charges.iter().enumerate() {
            if q == 0.0 {
                continue;
            }
            let (index, value, _) = spline(&at[i]);
            let target = if marked(i) { &mut im } else { &mut re };
            for j1 in 0..n {
                let row1 = index[0][j1] * k2;
                let q1 = q * value[0][j1];
                for j2 in 0..n {
                    let row2 = (row1 + index[1][j2]) * k3;
                    let q12 = q1 * value[1][j2];
                    for j3 in 0..n {
                        target[row2 + index[2][j3]] += q12 * value[2][j3];
                    }
                }
            }
        }
        transform_3d(&self.lines, &mut re, &mut im, false);

        // The energy by vector, the virial, and the spectrum of the potentials.
        let want = forces.is_some();
        let (mut pr, mut pi) = if want {
            (vec![0.0; total], vec![0.0; total])
        } else {
            (Vec::new(), Vec::new())
        };
        let split = group.is_some();
        let (mut rest, mut cross) = (0.0, 0.0);
        let a2 = self.alpha * self.alpha;
        for j1 in 0..k1 {
            let n1 = (k1 - j1) % k1;
            for j2 in 0..k2 {
                let n2 = (k2 - j2) % k2;
                for j3 in 0..k3 {
                    let t = (j1 * k2 + j2) * k3 + j3;
                    let g = self.influence[t];
                    if g == 0.0 {
                        continue;
                    }
                    let (ar, ai, br, bi) = if split {
                        let s = (n1 * k2 + n2) * k3 + (k3 - j3) % k3;
                        let (cr, ci, nr, ni) = (re[t], im[t], re[s], im[s]);
                        (
                            0.5 * (cr + nr),
                            0.5 * (ci - ni),
                            0.5 * (ci + ni),
                            -0.5 * (cr - nr),
                        )
                    } else {
                        (re[t], im[t], 0.0, 0.0)
                    };
                    let own = 0.5 * g * (ar * ar + ai * ai);
                    let across = g * (ar * br + ai * bi);
                    rest += own;
                    cross += across;
                    let e = own + lambda * across;
                    let m = [
                        signed(j1, k1).expect("G is zero on the Nyquist plane") / l[0],
                        signed(j2, k2).expect("G is zero on the Nyquist plane") / l[1],
                        signed(j3, k3).expect("G is zero on the Nyquist plane") / l[2],
                    ];
                    let mm = m[0] * m[0] + m[1] * m[1] + m[2] * m[2];
                    let stress = 2.0 * (1.0 / mm + PI * PI / a2);
                    for a in 0..3 {
                        for b in 0..3 {
                            let delta = if a == b { 1.0 } else { 0.0 };
                            virial[a][b] += COULOMB * e * (delta - stress * m[a] * m[b]);
                        }
                    }
                    if want {
                        // G(A + λB) for the rest, and i λ G A for the group.
                        pr[t] = g * (ar + lambda * br) - lambda * g * ai;
                        pi[t] = g * (ai + lambda * bi) + lambda * g * ar;
                    }
                }
            }
        }

        if let Some(f) = forces {
            transform_3d(&self.lines, &mut pr, &mut pi, true);
            for (i, &q) in charges.iter().enumerate() {
                if q == 0.0 {
                    continue;
                }
                let phi = if marked(i) { &pi } else { &pr };
                let (index, value, slope) = spline(&at[i]);
                let mut grad = [0.0; 3];
                for j1 in 0..n {
                    let row1 = index[0][j1] * k2;
                    let (w1, d1) = (value[0][j1], slope[0][j1]);
                    for j2 in 0..n {
                        let row2 = (row1 + index[1][j2]) * k3;
                        let (w2, d2) = (value[1][j2], slope[1][j2]);
                        for j3 in 0..n {
                            let v = phi[row2 + index[2][j3]];
                            let (w3, d3) = (value[2][j3], slope[2][j3]);
                            grad[0] += d1 * w2 * w3 * v;
                            grad[1] += w1 * d2 * w3 * v;
                            grad[2] += w1 * w2 * d3 * v;
                        }
                    }
                }
                for a in 0..3 {
                    f[i][a] -= COULOMB * q * grad[a] * grid[a] as f64 / l[a];
                }
            }
        }
        (COULOMB * rest, COULOMB * cross)
    }
}

/// Per grid index on an axis of `k` points, `(⟨η⟩, ⟨|η − ⟨η⟩|²⟩)` over a charge's position:
/// `−R/(1 + R)` and `Σ_{p≠0} r_p² / (1 + R)²`, from the alias ratios `r_p = (x/(x − p))ⁿ`,
/// `x = m/K`, summed to `|p| = 64` and the rest of `R` by its integral. Both are returned
/// directly rather than as `⟨1 + η⟩` and `⟨|1 + η|²⟩`, whose differences from one cancel. Index
/// `K/2` is not used.
fn axis_moments(order: usize, k: usize) -> Vec<(f64, f64)> {
    const P: i64 = 64;
    (0..k)
        .map(|j| {
            let Some(m) = signed(j, k) else {
                return (0.0, 0.0);
            };
            if m == 0.0 {
                return (0.0, 0.0);
            }
            let x = m / k as f64;
            let (mut r, mut r2) = (0.0, 0.0);
            for p in (-P..=P).filter(|&p| p != 0) {
                let ratio = x / (x - p as f64);
                let mut v = 1.0;
                for _ in 0..order {
                    v *= ratio;
                }
                r += v;
                r2 += v * v;
            }
            // Σ_{|p| > P} (x/(x − p))ⁿ: for even n, 2 xⁿ / ((n − 1) Pⁿ⁻¹) to leading order; for odd
            // n the two sides cancel to that order.
            if order & 1 == 0 {
                let mut xp = 1.0;
                for _ in 0..order {
                    xp *= x / P as f64;
                }
                r += 2.0 / (order as f64 - 1.0) * P as f64 * xp;
            }
            let d = 1.0 + r;
            (-r / d, r2 / (d * d))
        })
        .collect()
}

/// For one vector `m` from its three axes' [`axis_moments`]: `(⟨η⟩, ⟨|η|²⟩)` for
/// `1 + η = Π_a (1 + η_a)`, the axes independent. Expanded so that nothing cancels:
/// `⟨η⟩ = Π(1 + δ_a) − 1` term by term, and `⟨|η|²⟩ = ⟨η⟩² + Π(μ_a² + v_a) − Π μ_a²`, the second
/// difference as its sum over the non-empty sets of axes of `Π v_a` times the others' `μ²`.
fn vector_moments(axes: [(f64, f64); 3]) -> (f64, f64) {
    let [(d1, v1), (d2, v2), (d3, v3)] = axes;
    let mean = d1 + d2 + d3 + d1 * d2 + d1 * d3 + d2 * d3 + d1 * d2 * d3;
    let (s1, s2, s3) = (
        (1.0 + d1) * (1.0 + d1),
        (1.0 + d2) * (1.0 + d2),
        (1.0 + d3) * (1.0 + d3),
    );
    let spread = v1 * s2 * s3
        + s1 * v2 * s3
        + s1 * s2 * v3
        + v1 * v2 * s3
        + v1 * s2 * v3
        + s1 * v2 * v3
        + v1 * v2 * v3;
    (mean, mean * mean + spread)
}

/// `b / (k_e Σq²)` per metre: [`Pme::self_bias`] without its charges. Over every vector out to
/// where `w` is below `e⁻⁴⁰` of its largest or to the grid's edge, whichever is further.
fn bias_coefficient(cell: &PeriodicBox, alpha: f64, parameters: &PmeParameters) -> f64 {
    let n = parameters.order;
    let grid = parameters.grid;
    let l = cell.lengths();
    let a2 = alpha * alpha;
    let moments = grid.map(|k| axis_moments(n, k));
    let reach = [0, 1, 2]
        .map(|a| ((40f64.sqrt() * alpha * l[a] / PI).ceil() as i64 + 1).max(grid[a] as i64 / 2));
    let mut inside = 0.0;
    let mut outside = 0.0;
    for m1 in -reach[0]..=reach[0] {
        for m2 in -reach[1]..=reach[1] {
            for m3 in -reach[2]..=reach[2] {
                if m1 == 0 && m2 == 0 && m3 == 0 {
                    continue;
                }
                let m = [m1, m2, m3];
                let x = [0, 1, 2].map(|a| m[a] as f64 / l[a]);
                let mm = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
                let w = exp(-PI * PI * mm / a2) / mm;
                let zone = (0..3).all(|a| 2 * m[a].unsigned_abs() < grid[a] as u64);
                if zone {
                    let j = [0, 1, 2].map(|a| m[a].rem_euclid(grid[a] as i64) as usize);
                    // ⟨|1 + η|²⟩ − 1 = 2⟨η⟩ + ⟨|η|²⟩.
                    let (eta, eta2) = vector_moments([0, 1, 2].map(|a| moments[a][j[a]]));
                    inside += w * (2.0 * eta + eta2);
                } else {
                    outside += w;
                }
            }
        }
    }
    (inside - outside) / (2.0 * PI * cell.volume())
}

/// `σ_E / (k_e Q)` per metre: the module documentation's spread estimate without its charges.
fn spread_coefficient(cell: &PeriodicBox, alpha: f64, parameters: &PmeParameters) -> f64 {
    let n = parameters.order;
    let grid = parameters.grid;
    let l = cell.lengths();
    let a2 = alpha * alpha;
    let moments = grid.map(|k| axis_moments(n, k));
    let mut sum = 0.0;
    for j1 in 0..grid[0] {
        let Some(m1) = signed(j1, grid[0]) else {
            continue;
        };
        let x1 = m1 / l[0];
        for j2 in 0..grid[1] {
            let Some(m2) = signed(j2, grid[1]) else {
                continue;
            };
            let x2 = m2 / l[1];
            for j3 in 0..grid[2] {
                let Some(m3) = signed(j3, grid[2]) else {
                    continue;
                };
                if j1 == 0 && j2 == 0 && j3 == 0 {
                    continue;
                }
                let x3 = m3 / l[2];
                let mm = x1 * x1 + x2 * x2 + x3 * x3;
                let w = exp(-PI * PI * mm / a2) / mm;
                let (eta, eta2) = vector_moments([moments[0][j1], moments[1][j2], moments[2][j3]]);
                sum += w * w * (eta2 + eta * eta);
            }
        }
    }
    sum.sqrt() / (PI * cell.volume())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The twiddles: exact at the quarter turns, `√½` at the eighth, on the circle to two ulps,
    /// and the double angle of each the one at twice its index to a few ulps.
    #[test]
    fn the_roots_of_unity_are_on_the_circle() {
        for n in [4usize, 8, 64, 1024] {
            assert_eq!(unit_root(0, n), (1.0, 0.0));
            assert_eq!(unit_root(n / 4, n), (0.0, 1.0));
            assert_eq!(unit_root(n / 2, n), (-1.0, 0.0));
            assert_eq!(unit_root(3 * n / 4, n), (0.0, -1.0));
            if n >= 8 {
                let (c, s) = unit_root(n / 8, n);
                assert!((c - std::f64::consts::FRAC_1_SQRT_2).abs() <= f64::EPSILON);
                assert!((s - std::f64::consts::FRAC_1_SQRT_2).abs() <= f64::EPSILON);
            }
            for j in 0..n {
                let (c, s) = unit_root(j, n);
                assert!((c * c + s * s - 1.0).abs() <= 4.0 * f64::EPSILON, "{j}/{n}");
                if 2 * j < n {
                    let (c2, s2) = unit_root(2 * j, n);
                    assert!((c * c - s * s - c2).abs() <= 8.0 * f64::EPSILON, "{j}/{n}");
                    assert!((2.0 * c * s - s2).abs() <= 8.0 * f64::EPSILON, "{j}/{n}");
                }
            }
        }
    }

    /// One complex exponential lands in one bin, in each axis's transform of a 4 × 8 × 16 grid,
    /// forward with weight `K` and back to the input times `K`.
    #[test]
    fn one_frequency_lands_in_one_bin() {
        let grid = [4usize, 8, 16];
        let lines = grid.map(Line::new);
        let total = grid[0] * grid[1] * grid[2];
        for f in [[1usize, 3, 5], [0, 0, 1], [3, 7, 15], [2, 4, 8]] {
            let mut re = vec![0.0; total];
            let mut im = vec![0.0; total];
            for j1 in 0..grid[0] {
                for j2 in 0..grid[1] {
                    for j3 in 0..grid[2] {
                        let t = (j1 * grid[1] + j2) * grid[2] + j3;
                        // exp(2πi f·j/K), the product of three roots.
                        let (a, b) = unit_root(f[0] * j1, grid[0]);
                        let (c, d) = unit_root(f[1] * j2, grid[1]);
                        let (e, g) = unit_root(f[2] * j3, grid[2]);
                        let (x, y) = (a * c - b * d, a * d + b * c);
                        re[t] = x * e - y * g;
                        im[t] = x * g + y * e;
                    }
                }
            }
            let (re0, im0) = (re.clone(), im.clone());
            transform_3d(&lines, &mut re, &mut im, false);
            let peak = (f[0] * grid[1] + f[1]) * grid[2] + f[2];
            for t in 0..total {
                let want = if t == peak { total as f64 } else { 0.0 };
                assert!(
                    (re[t] - want).abs() < 1e-12 * total as f64
                        && im[t].abs() < 1e-12 * total as f64,
                    "{f:?}: bin {t} is ({}, {})",
                    re[t],
                    im[t]
                );
            }
            transform_3d(&lines, &mut re, &mut im, true);
            for t in 0..total {
                assert!((re[t] / total as f64 - re0[t]).abs() < 1e-13);
                assert!((im[t] / total as f64 - im0[t]).abs() < 1e-13);
            }
        }
    }
}
