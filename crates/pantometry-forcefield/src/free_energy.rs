//! Free energies from λ windows: molecular dynamics at each state of a schedule, and the two
//! estimators that turn what it sampled into a free-energy difference with an uncertainty —
//! thermodynamic integration and Bennett's acceptance ratio.
//!
//! # The windows
//!
//! [`Windows::new`] puts one [`Window`] at each state of a schedule of [`Lambda`]s. Each is
//! [`MolecularDynamics`] — BAOAB under a Langevin bath, as everywhere in this crate — on the
//! [`Alchemical`] Hamiltonian held at its state ([`AtLambda`]), from the same start, with its own
//! velocities and kicks. After `equilibration` steps it takes a [`Sample`] every `stride` steps
//! until it has `samples`: the gradient `∂U/∂λ` at its own state, and the energy differences to
//! the states either side, `U(λ_{k±1}) − U(λ_k)` at the same configuration, from the λ-dependent
//! part alone ([`Alchemical::coupling`]) so that nothing that cancels is rounded.
//!
//! **Deterministic per window, independent of chunking, and resumable.** Window `k`'s seed is
//! [`window_seed`]`(seed, k)`, which is all its velocities and kicks are drawn from, so a window is
//! the same whatever the others do and in whatever order they run. A sample is taken on the
//! dynamics' own step count, so a window advanced in one call or in many records the same
//! samples, to the bit ([`MolecularDynamics`]'s own guarantee carried one level up). Everything a
//! window is lives in the [`Window`] — positions, velocities, the step count and the samples —
//! so a [`Windows`] that has stopped part-way resumes by being advanced again, and a clone of it
//! is a checkpoint.
//!
//! **Constrained windows.** [`Windows::with_dynamics`] builds every window from a template
//! [`MolecularDynamics`] — its masses, frozen atoms, rigid waters and held bonds — given the
//! protocol's bath and thermalised from the window's seed; [`Windows::new`] is the same with an
//! unconstrained template, to the bit. A constraint is a property of the dynamics, not of the
//! Hamiltonian, so the same [`Alchemical`] serves both. Both end states carry the same
//! constraints, which is what a free energy between them needs: for bonds that share no atom, as a
//! solute's bonds to hydrogen do not, the constrained distribution needs no metric correction at
//! all ([`crate::shake`]).
//!
//! # Thermodynamic integration
//!
//! `ΔF = F(λ_last) − F(λ_first) = ∫ ⟨∂U/∂λ⟩_λ · dλ` along the schedule, a path of straight
//! segments in the three λs. [`Quadrature::Trapezoid`] is the trapezoid rule on that path: window
//! `k` contributes `⟨∇_λU · w_k⟩` with `w_k = (λ_{k+1} − λ_{k−1})/2` (one-sided at the ends),
//! and works on any schedule. [`Quadrature::Simpson`] is Simpson's rule, weights `(1, 4, 2, …,
//! 4, 1)/3`, for a schedule that is one straight line in equal steps with an even number of them.
//! Each window's term is the mean of one scalar series, `∇_λU · w_k`, and its standard error is
//! [`Estimate::of`]'s, from that series' own autocorrelation time; the windows are independent, so
//! the variances add. **The quadrature's own error is not in the uncertainty**: it is a bias, set
//! by the schedule, and a test of TI has to bound it separately.
//!
//! # Bennett's acceptance ratio
//!
//! > C. H. Bennett, *J. Comput. Phys.* **22**, 245 (1976) — **not opened here** (publisher's
//! > paywall); M. R. Shirts, E. Bair, G. Hooker and V. S. Pande, *Phys. Rev. Lett.* **91**, 140601
//! > (2003) — **not opened here** either: no readable copy was found, and it is not on arXiv.
//!
//! **Both are secondary, read from** M. R. Shirts and J. D. Chodera, "Statistically optimal
//! analysis of samples from multiple equilibrium states", *J. Chem. Phys.* **129**, 124105 (2008),
//! arXiv:0801.1426, which restates BAR as its two-state case (its eq 11 with K = 2, and Appendix
//! E for the equivalence and Shirts et al.'s variance). In reduced units (`u = U/k_BT`), with
//! `N_F` samples from state 0, `N_R` from state 1, `Δu(x) = u₁(x) − u₀(x)` **for both**, and
//! `M = ln(N_F/N_R)`, the estimate `Δf = f₁ − f₀` solves
//!
//! `Σ_F 1/(1 + e^(M + Δu − Δf)) = Σ_R 1/(1 + e^(−M − Δu + Δf))`,
//!
//! whose left side rises and right side falls with `Δf`, so the root is unique. [`bar`] finds it
//! **to the last bit**, from a bracket where each side's terms are all below `e^(−50)` of the
//! other's: Newton's step wherever it lands inside the bracket and bisection wherever it does not,
//! each evaluation moving one end of the bracket, until no float is left strictly inside it. No
//! tolerance, no iteration count to choose, the same answer on every run; about a tenth of the
//! evaluations bisection alone took. The asymptotic variance (Shirts et al. 2003, as restated) is
//!
//! `σ²(Δf) = [Σ_n 1/(2 + 2 cosh(M + Δu_n − Δf))]⁻¹ − (1/N_F + 1/N_R)`,
//!
//! the sum over every sample of both sets. **It holds for independent samples**, and [`bar`]
//! reports it; the overlap [`Bar::overlap`] is the two-state MBAR overlap-matrix element
//! `O₀₁ = (1/N_F) Σ_n 1/(2 + 2 cosh(…))` — the same sum — which is ½ for two identical states and
//! falls to 0 as they part.
//!
//! **Windows give neither independent samples nor independent intervals**, and the variance the
//! windows report is not that formula. Window `k` is the forward set of interval `k` and the
//! reverse set of interval `k − 1`, with the same configurations in both, so the two intervals'
//! estimates are correlated, and summing their variances leaves out a covariance. **Measured**, on
//! the particle of `tests/the_free_energy_against_closed_forms.rs` over 48 seeds: the summed
//! variances, from samples thinned to every `⌈2τ⌉`-th, put σ̂ at 0.0198 `k_BT` against a spread of
//! 0.0345 between seeds — 1.7 times too small — and 0.0102 against 0.0166 at a stride five times
//! longer, where the samples are close to independent. The thinning was not the fault; the
//! covariance was. With the delta method below, the same 48 seeds give σ̂ 0.0282 against a spread
//! of 0.0304, and 0.0131 against 0.0135 at the longer stride. On two harmonic wells over 200 MD
//! campaigns the z-scores' variance is 0.99 at a friction that decorrelates in a few samples and
//! 1.29 where successive samples are strongly correlated (γ = 0.2ω, stride 2): there σ̂ is about
//! 12% small, which is the autocorrelation-time estimator's, not BAR's (TI's is 1.35 there).
//! On exact, independent samples along a chain of five states, over 1200 seeds, the summed
//! variances understate σ̂ by 1.26 (z variance 1.555 against 0.982 window by window): the
//! covariance alone, with nothing thinned. The particle's 1.7 is that and the thinning together.
//!
//! **So the variance is the delta method's, window by window** ([`bar_correlated`],
//! [`bennett_chain`], [`Windows::bennett_total`]). Perturbing BAR's equation, `δΔf = −(δa − δb)/S` with `a` and `b`
//! its two sums and `S = Σ_n 1/(2 + 2 cosh(M + Δu_n − Δf))` their slope. `a` is a sum over one
//! window's samples and `b` over the next's, so the total `Σ_k Δf_k` moves by a sum, over
//! windows, of one series per window,
//!
//! `h_j(x) = −F(M_j + Δu_{j→j+1}(x) − Δf_j)/S_j + F(−M_{j−1} − Δu_{j−1→j}(x) + Δf_{j−1})/S_{j−1}`,
//!
//! `F(y) = 1/(1 + eʸ)`, the first term absent for the last window and the second for the first.
//! The windows are independent, so `σ²(ΣΔf) = Σ_j N_j² σ²(h̄_j)`, with `σ(h̄_j)` [`Estimate::of`]'s
//! standard error of the series' mean from its own autocorrelation time. Each interval's variance
//! alone is the same with only its own two terms. **Every sample is used**: nothing is thinned. For
//! independent samples of two states this is asymptotically Shirts et al.'s variance — it is the
//! robust ("sandwich") form of the same estimator's variance — and the tests hold both to the
//! spread of exact samples. **On independent samples it is about 4% conservative**: over 1600
//! seeds of 400 and 700 samples its z variance is 0.916 where Shirts's is 0.963, because
//! [`Estimate::of`]'s `τ`, a windowed sum of noisy correlations, is clipped at ½ from below and so
//! errs upward only. On correlated samples, where `τ` is well above ½, it reads 1.001.
//!
//! **Not MBAR.** MBAR would use every window's samples at every state, which needs each sample's
//! energy at every state and a pseudo-inverse of a `K × K` matrix for its covariance; on a
//! schedule chosen so that neighbours overlap and others barely do, what it adds over BAR between
//! neighbours is small, and BAR is what Mobley et al. 2007 used on this system. Not here, and
//! nothing here prevents adding it: a window's samples would need the energies at every state.

use crate::alchemy::{Alchemical, AtLambda, Lambda};
use crate::complex::Estimate;
use crate::dynamics::{Bath, MolecularDynamics};
use pantometry_units::BOLTZMANN;

/// The key mixed with a window's index into its seed: see [`window_seed`].
pub const WINDOW_STREAM: u64 = 0xA076_1D64_78BD_642F;

/// Window `index`'s seed for a campaign seeded `seed`: `seed ^ (index + 1) · WINDOW_STREAM`,
/// wrapping. It seeds both the window's velocities and its kicks, which [`MolecularDynamics`]
/// draws from separate streams of it; distinct for every window of one campaign.
pub fn window_seed(seed: u64, index: usize) -> u64 {
    seed ^ (index as u64 + 1).wrapping_mul(WINDOW_STREAM)
}

/// How each window is run.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Protocol {
    /// The dynamics step, seconds.
    pub time_step: f64,
    /// The bath's temperature, kelvin; also the `k_BT` of the estimators.
    pub temperature: f64,
    /// The bath's friction, per second.
    pub friction: f64,
    /// Steps discarded at the start of each window.
    pub equilibration: u64,
    /// Steps between samples.
    pub stride: u64,
    /// Samples per window.
    pub samples: usize,
    /// The campaign's seed; window `k` uses [`window_seed`]`(seed, k)`.
    pub seed: u64,
}

impl Protocol {
    /// `k_BT`, joules.
    pub fn kt(&self) -> f64 {
        BOLTZMANN.to_si() * self.temperature
    }

    /// The steps a window takes in all: `equilibration + stride × samples`.
    pub fn steps_per_window(&self) -> u64 {
        self.equilibration + self.stride * self.samples as u64
    }
}

/// What a window records at one step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sample {
    /// The window's dynamics' step count.
    pub step: u64,
    /// `[∂U/∂λ_r, ∂U/∂λ_e, ∂U/∂λ_v]` at the window's state, joules per molecule.
    pub gradient: [f64; 3],
    /// `U(λ_{k−1}) − U(λ_k)` at this configuration, joules per molecule; `NaN` for the first
    /// window.
    pub to_previous: f64,
    /// `U(λ_{k+1}) − U(λ_k)`, joules per molecule; `NaN` for the last window.
    pub to_next: f64,
}

/// One state of a schedule and the dynamics run at it: see the module documentation.
#[derive(Clone, Debug, PartialEq)]
pub struct Window {
    index: usize,
    lambda: Lambda,
    previous: Option<Lambda>,
    next: Option<Lambda>,
    protocol: Protocol,
    md: MolecularDynamics,
    at: Vec<[f64; 3]>,
    samples: Vec<Sample>,
}

impl Window {
    /// Window `index` of `schedule` for `hamiltonian`: atoms of `masses` (kilograms) at `start`
    /// (metres), those `frozen` marks held, thermalised at the protocol's temperature from
    /// [`window_seed`]`(protocol.seed, index)`, in a Langevin bath of its friction.
    ///
    /// # Panics
    ///
    /// If `index` is past the schedule, a state is not valid, or the lengths do not agree.
    pub fn new(
        hamiltonian: &(impl Alchemical + ?Sized),
        schedule: &[Lambda],
        index: usize,
        start: &[[f64; 3]],
        masses: Vec<f64>,
        frozen: Vec<bool>,
        protocol: Protocol,
    ) -> Window {
        Window::with_dynamics(
            hamiltonian,
            schedule,
            index,
            start,
            &MolecularDynamics::new(masses).with_frozen(frozen),
            protocol,
        )
    }

    /// The same window with its dynamics built from `dynamics` — its masses, frozen atoms and
    /// constraints ([`MolecularDynamics::with_constraints`],
    /// [`MolecularDynamics::with_bond_constraints`]) — given the protocol's bath and thermalised
    /// from [`window_seed`]`(protocol.seed, index)`. [`Window::new`] is this with
    /// `MolecularDynamics::new(masses).with_frozen(frozen)`, to the bit. `start` must satisfy the
    /// constraints, and every rigid molecule be whole.
    ///
    /// # Panics
    ///
    /// As [`Window::new`], and if `dynamics` has taken a step.
    pub fn with_dynamics(
        hamiltonian: &(impl Alchemical + ?Sized),
        schedule: &[Lambda],
        index: usize,
        start: &[[f64; 3]],
        dynamics: &MolecularDynamics,
        protocol: Protocol,
    ) -> Window {
        assert!(
            index < schedule.len(),
            "window {index} of {}",
            schedule.len()
        );
        assert!(
            schedule.iter().all(Lambda::is_valid),
            "every λ must be in [0, 1]"
        );
        assert_eq!(start.len(), hamiltonian.len(), "one position per atom");
        assert_eq!(
            dynamics.masses().len(),
            start.len(),
            "one mass per position"
        );
        assert_eq!(dynamics.steps(), 0, "a window's dynamics starts unstepped");
        assert!(protocol.stride > 0, "a sample every zero steps");
        let seed = window_seed(protocol.seed, index);
        let mut md = dynamics
            .clone()
            .with_bath(Bath::Langevin {
                temperature: protocol.temperature,
                friction: protocol.friction,
                seed,
            })
            .thermalised(start, protocol.temperature, seed);
        let lambda = schedule[index];
        md.prepare(
            &AtLambda {
                hamiltonian,
                lambda,
            },
            start,
        );
        Window {
            index,
            lambda,
            previous: index.checked_sub(1).map(|i| schedule[i]),
            next: schedule.get(index + 1).copied(),
            protocol,
            md,
            at: start.to_vec(),
            samples: Vec::new(),
        }
    }

    /// Up to `max_steps` more steps, stopping once the window has all its samples; returns how
    /// many were taken. The same samples, to the bit, however the steps are cut into calls.
    /// `hamiltonian` must be the one the window was built for.
    pub fn advance(&mut self, hamiltonian: &(impl Alchemical + ?Sized), max_steps: u64) -> u64 {
        let potential = AtLambda {
            hamiltonian,
            lambda: self.lambda,
        };
        let (eq, stride) = (self.protocol.equilibration, self.protocol.stride);
        let mut taken = 0;
        // Clippy on current stable suggests `u64::is_multiple_of`, stabilised in 1.87; this crate
        // builds on 1.78.
        #[allow(clippy::manual_is_multiple_of)]
        while taken < max_steps && !self.is_complete() {
            self.md
                .step(&potential, &mut self.at, self.protocol.time_step);
            taken += 1;
            let s = self.md.steps();
            if s > eq && (s - eq) % stride == 0 {
                // The window's own state and its neighbours at one configuration, in one call, so
                // that a Hamiltonian can share what they have in common.
                let mut states = vec![self.lambda];
                states.extend(self.previous);
                states.extend(self.next);
                let c = hamiltonian.couplings(&self.at, &states);
                let here = c[0];
                let mut others = c[1..].iter();
                let mut to = |other: Option<Lambda>| {
                    other.map_or(f64::NAN, |_| {
                        others.next().expect("one coupling per neighbour").energy - here.energy
                    })
                };
                let to_previous = to(self.previous);
                let to_next = to(self.next);
                self.samples.push(Sample {
                    step: s,
                    gradient: here.gradient,
                    to_previous,
                    to_next,
                });
            }
        }
        taken
    }

    /// Whether the window has all its samples.
    pub fn is_complete(&self) -> bool {
        self.samples.len() >= self.protocol.samples
    }

    /// Its index in the schedule.
    pub fn index(&self) -> usize {
        self.index
    }

    /// Its state.
    pub fn lambda(&self) -> Lambda {
        self.lambda
    }

    /// The samples so far, in order.
    pub fn samples(&self) -> &[Sample] {
        &self.samples
    }

    /// The positions now, metres.
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.at
    }

    /// The dynamics.
    pub fn dynamics(&self) -> &MolecularDynamics {
        &self.md
    }
}

/// How [`Windows::thermodynamic_integration`] integrates: see the module documentation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Quadrature {
    /// The trapezoid rule along the schedule's path: any schedule.
    Trapezoid,
    /// Simpson's rule: a schedule that is one straight line in equal steps, an even number of
    /// them.
    Simpson,
}

/// A free-energy difference and its standard error: joules per molecule, except from
/// [`bennett_chain`], which is in reduced units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeEnergy {
    /// The estimate.
    pub value: f64,
    /// Its standard error.
    pub error: f64,
}

/// One BAR estimate, in reduced units: see [`bar`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bar {
    /// `Δf = (F₁ − F₀)/k_BT`.
    pub delta: f64,
    /// Its asymptotic variance (Shirts et al. 2003), `(k_BT)⁻²`-less: reduced.
    pub variance: f64,
    /// The two-state overlap `O₀₁`, ½ for identical states.
    pub overlap: f64,
    /// Samples in the forward set.
    pub forward: usize,
    /// Samples in the reverse set.
    pub reverse: usize,
}

impl Bar {
    /// `√variance`.
    pub fn error(&self) -> f64 {
        self.variance.max(0.0).sqrt()
    }
}

/// `1/(1 + eˣ)`, without overflow for large `|x|`.
fn fermi(x: f64) -> f64 {
    if x > 0.0 {
        let e = (-x).exp();
        e / (1.0 + e)
    } else {
        1.0 / (1.0 + x.exp())
    }
}

/// Bennett's acceptance ratio between two states, reduced units: `forward` holds
/// `Δu = u₁(x) − u₀(x)` on samples from state 0 and `reverse` **the same difference,
/// `u₁(x) − u₀(x)`**, on samples from state 1. The variance is Shirts et al.'s and holds for
/// independent samples; [`bar_correlated`] is the one for a correlated series. See the module documentation for the equation, the bisection and the variance.
///
/// # Panics
///
/// If either set is empty or holds a value that is not finite.
pub fn bar(forward: &[f64], reverse: &[f64]) -> Bar {
    assert!(
        !forward.is_empty() && !reverse.is_empty(),
        "BAR needs samples from both states"
    );
    assert!(
        forward.iter().chain(reverse).all(|w| w.is_finite()),
        "every energy difference must be finite"
    );
    let (nf, nr) = (forward.len() as f64, reverse.len() as f64);
    let m = (nf / nr).ln();
    // g rises with Δf: the forward sum rises and the reverse sum falls. Its slope is
    // S = Σ_n F(x_n) F(−x_n) over both sets, x = M + Δu − Δf.
    let g = |df: f64| -> (f64, f64) {
        let (mut f, mut slope) = (0.0, 0.0);
        for w in forward {
            let x = m + w - df;
            f += fermi(x);
            slope += fermi(x) * fermi(-x);
        }
        let mut r = 0.0;
        for w in reverse {
            let x = m + w - df;
            r += fermi(-x);
            slope += fermi(x) * fermi(-x);
        }
        (f - r, slope)
    };
    let lowest = forward
        .iter()
        .chain(reverse)
        .copied()
        .fold(f64::INFINITY, f64::min);
    let highest = forward
        .iter()
        .chain(reverse)
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let margin = 50.0 + m.abs();
    let (mut lo, mut hi) = (lowest - margin, highest + margin);
    // Newton's step where it lands inside the bracket, bisection where it does not; every
    // evaluation moves one end of the bracket to it, so the bracket shrinks at every step, and the
    // search stops when no float strictly inside it is left to try. Deterministic, and no
    // tolerance to choose.
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
    let mut x = (0.5 * (mean(forward) + mean(reverse))).clamp(lo, hi);
    let delta = loop {
        let (value, slope) = g(x);
        if value == 0.0 {
            break x;
        }
        if value < 0.0 {
            lo = x;
        } else {
            hi = x;
        }
        let newton = x - value / slope;
        let next = if newton > lo && newton < hi && newton.is_finite() {
            newton
        } else {
            0.5 * (lo + hi)
        };
        if next <= lo || next >= hi {
            // Nothing strictly inside the bracket: lo and hi are adjacent floats, or equal.
            break if g(lo).0.abs() <= g(hi).0.abs() {
                lo
            } else {
                hi
            };
        }
        x = next;
    };
    let s: f64 = forward
        .iter()
        .chain(reverse)
        .map(|w| {
            let x = m + w - delta;
            // 1/(2 + 2 cosh x) = σ(x)(1 − σ(x)), without overflow.
            fermi(x) * fermi(-x)
        })
        .sum();
    Bar {
        delta,
        variance: 1.0 / s - (1.0 / nf + 1.0 / nr),
        overlap: s / nf,
        forward: forward.len(),
        reverse: reverse.len(),
    }
}

/// `[F(M + Δu − Δf)]` over the forward set and `[F(−M − Δu + Δf)]` over the reverse set, and
/// `S`: the terms of BAR's equation at its root, from which [`bar_correlated`] and
/// [`Windows::bennett_total`] build the delta method's variance.
fn influence(forward: &[f64], reverse: &[f64], b: &Bar) -> (Vec<f64>, Vec<f64>, f64) {
    let m = (forward.len() as f64 / reverse.len() as f64).ln();
    let f: Vec<f64> = forward.iter().map(|w| fermi(m + w - b.delta)).collect();
    let r: Vec<f64> = reverse.iter().map(|w| fermi(-m - w + b.delta)).collect();
    (f, r, b.overlap * forward.len() as f64)
}

/// The variance of the sum of a series, from [`Estimate::of`]'s standard error of its mean: `(N σ)²`.
fn sum_variance(xs: &[f64]) -> f64 {
    let e = Estimate::of(xs);
    let n = xs.len() as f64;
    (n * e.error) * (n * e.error)
}

/// [`bar`] on two **correlated** series — every sample of one window each, say — with the
/// variance from the delta method on each side's series, each from its own autocorrelation time
/// ([`Estimate::of`]): `σ² = [N_F² σ²(F̄) + N_R² σ²(R̄)]/S²` with `F` and `R` the two sides' terms at
/// the root and `S` their slope. See the module documentation. The estimate and the overlap are
/// [`bar`]'s.
///
/// # Panics
///
/// As [`bar`].
pub fn bar_correlated(forward: &[f64], reverse: &[f64]) -> Bar {
    let b = bar(forward, reverse);
    let (f, r, s) = influence(forward, reverse, &b);
    Bar {
        variance: (sum_variance(&f) + sum_variance(&r)) / (s * s),
        ..b
    }
}

/// BAR along a chain of states, reduced units: `windows[j]` holds window `j`'s samples as
/// `(u_{j−1} − u_j, u_{j+1} − u_j)` per sample, in reduced units (the first window's first series
/// and the last's second are not read). The value is the sum of each interval's [`bar`]
/// estimate, from every sample; the error is the delta method's, window by window — one series
/// `h_j` per window, each from its own autocorrelation time ([`Estimate::of`]), so that the two
/// intervals a window belongs to are correlated as they are. See the module documentation for
/// `h_j`.
///
/// # Panics
///
/// If there are fewer than two windows, or a window has no sample, or as [`bar`].
pub fn bennett_chain(windows: &[(Vec<f64>, Vec<f64>)]) -> FreeEnergy {
    let n = windows.len();
    assert!(n >= 2, "a chain needs two states");
    let terms: Vec<(Vec<f64>, Vec<f64>, f64, f64)> = (0..n - 1)
        .map(|k| {
            let forward = &windows[k].1;
            // The reverse set's u₁ − u₀ is minus what window k + 1 has towards k.
            let reverse: Vec<f64> = windows[k + 1].0.iter().map(|x| -x).collect();
            let b = bar(forward, &reverse);
            let (tf, tr, s) = influence(forward, &reverse, &b);
            (tf, tr, s, b.delta)
        })
        .collect();
    let mut variance = 0.0;
    for j in 0..n {
        let len = windows[j].0.len().max(windows[j].1.len());
        let h: Vec<f64> = (0..len)
            .map(|i| {
                let mut x = 0.0;
                if j + 1 < n {
                    let (tf, _, s, _) = &terms[j];
                    x -= tf[i] / s;
                }
                if j > 0 {
                    let (_, tr, s, _) = &terms[j - 1];
                    x += tr[i] / s;
                }
                x
            })
            .collect();
        variance += sum_variance(&h);
    }
    FreeEnergy {
        value: terms.iter().map(|t| t.3).sum(),
        error: variance.sqrt(),
    }
}

/// A schedule's windows: see the module documentation.
#[derive(Clone, Debug, PartialEq)]
pub struct Windows {
    schedule: Vec<Lambda>,
    protocol: Protocol,
    windows: Vec<Window>,
}

impl Windows {
    /// One [`Window`] at each state of `schedule`, every one from `start`, as [`Window::new`].
    ///
    /// # Panics
    ///
    /// If the schedule has fewer than two states, or as [`Window::new`].
    pub fn new(
        hamiltonian: &(impl Alchemical + ?Sized),
        schedule: Vec<Lambda>,
        start: &[[f64; 3]],
        masses: Vec<f64>,
        frozen: Vec<bool>,
        protocol: Protocol,
    ) -> Windows {
        Windows::with_dynamics(
            hamiltonian,
            schedule,
            start,
            &MolecularDynamics::new(masses).with_frozen(frozen),
            protocol,
        )
    }

    /// One [`Window`] at each state of `schedule`, every one from `start`, as
    /// [`Window::with_dynamics`]: for constrained dynamics — rigid water, bonds to hydrogen —
    /// and anything else `dynamics` carries. [`Windows::new`] is this with
    /// `MolecularDynamics::new(masses).with_frozen(frozen)`, to the bit.
    ///
    /// # Panics
    ///
    /// If the schedule has fewer than two states, or as [`Window::with_dynamics`].
    pub fn with_dynamics(
        hamiltonian: &(impl Alchemical + ?Sized),
        schedule: Vec<Lambda>,
        start: &[[f64; 3]],
        dynamics: &MolecularDynamics,
        protocol: Protocol,
    ) -> Windows {
        assert!(schedule.len() >= 2, "a schedule needs two states");
        let windows = (0..schedule.len())
            .map(|k| Window::with_dynamics(hamiltonian, &schedule, k, start, dynamics, protocol))
            .collect();
        Windows {
            schedule,
            protocol,
            windows,
        }
    }

    /// Advances every window by up to `max_steps` steps, in order, and says whether all are
    /// complete. Call again to continue.
    pub fn run(&mut self, hamiltonian: &(impl Alchemical + ?Sized), max_steps: u64) -> bool {
        for w in &mut self.windows {
            w.advance(hamiltonian, max_steps);
        }
        self.is_complete()
    }

    /// Whether every window has all its samples.
    pub fn is_complete(&self) -> bool {
        self.windows.iter().all(Window::is_complete)
    }

    /// The schedule.
    pub fn schedule(&self) -> &[Lambda] {
        &self.schedule
    }

    /// The protocol.
    pub fn protocol(&self) -> &Protocol {
        &self.protocol
    }

    /// The windows, in the schedule's order.
    pub fn windows(&self) -> &[Window] {
        &self.windows
    }

    /// Window `k`'s mutable state, for advancing one window alone.
    ///
    /// # Panics
    ///
    /// If `k` is past the schedule.
    pub fn window_mut(&mut self, k: usize) -> &mut Window {
        &mut self.windows[k]
    }

    /// The trapezoid weights `w_k` of the module documentation, one per window.
    fn trapezoid_weights(&self) -> Vec<[f64; 3]> {
        let s = &self.schedule;
        let n = s.len();
        (0..n)
            .map(|k| {
                let lo = s[k.saturating_sub(1)].components();
                let hi = s[(k + 1).min(n - 1)].components();
                [0, 1, 2].map(|c| 0.5 * (hi[c] - lo[c]))
            })
            .collect()
    }

    /// Simpson's weights, or `None` when the schedule is not one straight line in an even number
    /// of equal steps.
    // Clippy on current stable suggests `usize::is_multiple_of`, stabilised in 1.87; this crate
    // builds on 1.78.
    #[allow(clippy::manual_is_multiple_of)]
    fn simpson_weights(&self) -> Option<Vec<[f64; 3]>> {
        let s = &self.schedule;
        let n = s.len();
        if n < 3 || (n - 1) % 2 != 0 {
            return None;
        }
        let step = {
            let (a, b) = (s[0].components(), s[1].components());
            [b[0] - a[0], b[1] - a[1], b[2] - a[2]]
        };
        for k in 1..n {
            let (a, b) = (s[k - 1].components(), s[k].components());
            let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            // Equal to rounding of the schedule's own arithmetic: 1e-12 of a unit range.
            if (0..3).any(|c| (d[c] - step[c]).abs() > 1e-12) {
                return None;
            }
        }
        Some(
            (0..n)
                .map(|k| {
                    let c = if k == 0 || k == n - 1 {
                        1.0
                    } else if k % 2 == 1 {
                        4.0
                    } else {
                        2.0
                    };
                    step.map(|x| x * c / 3.0)
                })
                .collect(),
        )
    }

    /// Each window's TI term, `⟨∇_λU · w_k⟩` with its [`Estimate`], joules: the series the
    /// estimate is built on, so that a caller can see what each window contributed. `None` for
    /// [`Quadrature::Simpson`] on a schedule it does not fit, or a window with no sample.
    pub fn ti_terms(&self, rule: Quadrature) -> Option<Vec<Estimate>> {
        let weights = match rule {
            Quadrature::Trapezoid => self.trapezoid_weights(),
            Quadrature::Simpson => self.simpson_weights()?,
        };
        self.windows
            .iter()
            .zip(&weights)
            .map(|(w, wk)| {
                let ys: Vec<f64> = w
                    .samples
                    .iter()
                    .map(|s| (0..3).map(|c| s.gradient[c] * wk[c]).sum())
                    .collect();
                (!ys.is_empty()).then(|| Estimate::of(&ys))
            })
            .collect()
    }

    /// `F(λ_last) − F(λ_first)` by thermodynamic integration with `rule`, and its standard error
    /// from each window's autocorrelation time: see the module documentation. `None` as
    /// [`Windows::ti_terms`].
    pub fn thermodynamic_integration(&self, rule: Quadrature) -> Option<FreeEnergy> {
        let terms = self.ti_terms(rule)?;
        Some(FreeEnergy {
            value: terms.iter().map(|t| t.mean).sum(),
            error: terms.iter().map(|t| t.error * t.error).sum::<f64>().sqrt(),
        })
    }

    /// Interval `k`'s two series, reduced: window `k`'s `Δu_{k→k+1}` and window `k + 1`'s
    /// `u_{k+1} − u_k`, which is minus what it recorded towards `k`.
    fn interval(&self, k: usize) -> (Vec<f64>, Vec<f64>) {
        let kt = self.protocol.kt();
        let forward = self.windows[k]
            .samples
            .iter()
            .map(|s| s.to_next / kt)
            .collect();
        let reverse = self.windows[k + 1]
            .samples
            .iter()
            .map(|s| -s.to_previous / kt)
            .collect();
        (forward, reverse)
    }

    /// BAR between each pair of neighbouring windows, reduced units, on every sample, each with
    /// its own delta-method variance ([`bar_correlated`]). One per interval.
    ///
    /// # Panics
    ///
    /// If a window has no sample.
    pub fn bennett(&self) -> Vec<Bar> {
        (0..self.windows.len() - 1)
            .map(|k| {
                let (f, r) = self.interval(k);
                bar_correlated(&f, &r)
            })
            .collect()
    }

    /// `F(λ_last) − F(λ_first)` by BAR, joules: [`bennett_chain`] on the windows' samples, so that
    /// the covariance between intervals that share a window is in the error. See the module
    /// documentation.
    ///
    /// # Panics
    ///
    /// As [`Windows::bennett`].
    pub fn bennett_total(&self) -> FreeEnergy {
        let kt = self.protocol.kt();
        let chain: Vec<(Vec<f64>, Vec<f64>)> = self
            .windows
            .iter()
            .map(|w| {
                (
                    w.samples.iter().map(|s| s.to_previous / kt).collect(),
                    w.samples.iter().map(|s| s.to_next / kt).collect(),
                )
            })
            .collect();
        let f = bennett_chain(&chain);
        FreeEnergy {
            value: kt * f.value,
            error: kt * f.error,
        }
    }
}
