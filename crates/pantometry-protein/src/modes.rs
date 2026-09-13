//! The motions a structure has, and how far each residue moves.
//!
//! # The six that are not motions
//!
//! A protein floating in nothing can be moved and turned without stretching a single spring, so
//! six directions in the `3n`-dimensional space cost no energy: three translations and three
//! rotations. They are **exact** null vectors of the Hessian, not small eigenvalues, because the
//! potential is a function of interatomic distances and a rigid motion changes none of them.
//!
//! Everything downstream divides by an eigenvalue, so those six have to be taken out first, and
//! *counted* rather than assumed:
//!
//! - a structure whose nodes lie on a straight line has **five**, because turning it about its
//!   own axis does nothing,
//! - a network that has fallen into two pieces has **twelve**, and its lowest non-rigid mode is
//!   one piece drifting away from the other.
//!
//! [`Modes::rigid_body_modes`] is what the count is, and [`Modes::separation`] is the evidence
//! that the count is right: the ratio between the first real eigenvalue and the largest of the
//! ones called zero. For a protein it is ten orders of magnitude or more, so nothing here rests
//! on where a threshold was put.
//!
//! # What a fluctuation is
//!
//! At equilibrium the covariance of the displacements is `k_BT H⁺`, so
//!
//! ```text
//! ⟨Δr_i²⟩ = k_BT Σ_{k not rigid} (1/λ_k) (u²_k,3i + u²_k,3i+1 + u²_k,3i+2)
//! ```
//!
//! The stiffness `γ` divides every one of these equally, so it sets the scale and changes no
//! shape — which is what makes [`correlation`] against a measured B-factor a prediction rather
//! than a fit.

use crate::network::Network;
use crate::spectrum::{Spectrum, Symmetric};
use pantometry_units::{Qty, Stiffness, Temperature, BOLTZMANN};

/// The motions of one [`Network`].
#[derive(Clone, Debug, PartialEq)]
pub struct Modes {
    spectrum: Spectrum,
    nodes: usize,
    rigid: usize,
}

impl Modes {
    /// Diagonalise the network's Hessian: the **whole complex**, ligand included.
    ///
    /// These are the motions to animate, and for a network with no ligand they are also the
    /// motions to take statistics from. With a ligand they are not —
    /// [`Modes::of_the_protein`] is, and [`Network::protein_hessian`] says why in one measured
    /// number.
    ///
    /// # Panics
    ///
    /// If the eigensolver did not converge, which for a Hessian means something is wrong with the
    /// structure rather than with the budget — see [`Spectrum::converged`].
    pub fn of(network: &Network) -> Modes {
        Modes::from_hessian(network.hessian(), network.nodes())
    }

    /// The motions of the **protein's own coordinates**, with the ligand integrated out.
    ///
    /// This is what a residue's fluctuation is a property of, and the two agree exactly when
    /// there is no ligand. See [`Network::protein_hessian`].
    ///
    /// # Panics
    ///
    /// As [`Modes::of`].
    pub fn of_the_protein(network: &Network) -> Modes {
        Modes::from_hessian(network.protein_hessian(), network.residues())
    }

    fn from_hessian(h: Symmetric, nodes: usize) -> Modes {
        let spectrum = h.eigen();
        assert!(
            spectrum.converged(),
            "the Hessian of {nodes} nodes did not diagonalise in the sweep budget, leaving an \
             off-diagonal norm of {:e}",
            spectrum.off_diagonal_norm()
        );
        let largest = spectrum.values().last().copied().unwrap_or(0.0);
        // Where to cut is not a judgement call as long as the gap is wide, and `separation`
        // reports how wide it was. The rounding floor sits at about `10⁻¹²` of the largest
        // eigenvalue and a protein's first real mode at about `10⁻⁴`; `√ε` is between them, four
        // orders clear of each.
        let floor = f64::EPSILON.sqrt() * largest.abs();
        let rigid = spectrum
            .values()
            .iter()
            .filter(|v| v.abs() <= floor)
            .count();
        Modes {
            spectrum,
            nodes,
            rigid,
        }
    }

    /// How many nodes the network had.
    pub fn nodes(&self) -> usize {
        self.nodes
    }

    /// How many eigenvalues came back at zero: the rigid-body motions.
    ///
    /// Six for a protein. Anything else is a statement about the structure — five for a
    /// collinear one, a multiple of six for a network in pieces — and never about this code.
    pub fn rigid_body_modes(&self) -> usize {
        self.rigid
    }

    /// How many modes are left after the rigid ones.
    pub fn len(&self) -> usize {
        self.spectrum.order() - self.rigid
    }

    /// Whether there are none, which happens only for a structure of one or two nodes.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The ratio of the first real eigenvalue to the largest of the ones called rigid.
    ///
    /// The number that says whether [`Modes::rigid_body_modes`] is a measurement or a
    /// coincidence. Infinite when every rigid eigenvalue came back exactly zero.
    pub fn separation(&self) -> f64 {
        if self.rigid == 0 || self.rigid >= self.spectrum.order() {
            return f64::INFINITY;
        }
        let noise = self.spectrum.values()[..self.rigid]
            .iter()
            .fold(0.0f64, |m, v| m.max(v.abs()));
        if noise == 0.0 {
            f64::INFINITY
        } else {
            self.spectrum.values()[self.rigid] / noise
        }
    }

    /// The eigenvalue of mode `k`, counting from the softest non-rigid one.
    ///
    /// # Panics
    ///
    /// If `k` is at or past [`Modes::len`].
    pub fn eigenvalue(&self, k: usize) -> Stiffness {
        assert!(k < self.len(), "mode {k} of {}", self.len());
        Qty::from_si(self.spectrum.values()[self.rigid + k])
    }

    /// The shape of mode `k`: `3n` components, three per node, unit length.
    ///
    /// **Its sign means nothing.** `v` and `−v` are the same motion, and the rule
    /// [`Spectrum::vector`] uses to pick between them is a tie whenever the structure is
    /// symmetric. Anything read off a mode has to be invariant under `v → −v`, which is why
    /// [`Modes::overlap`] is an absolute value.
    ///
    /// # Panics
    ///
    /// If `k` is at or past [`Modes::len`].
    pub fn shape(&self, k: usize) -> &[f64] {
        assert!(k < self.len(), "mode {k} of {}", self.len());
        self.spectrum.vector(self.rigid + k)
    }

    /// `|cos|` between this mode `k` and `other`'s mode `j`.
    ///
    /// One for the same motion, zero for unrelated ones. Absolute because of the sign: a
    /// comparison that used the signed cosine would report `−1` — as unrelated as it gets — for
    /// two computations of the identical motion.
    ///
    /// # Panics
    ///
    /// If either index is out of range, or if the two mode sets are of different sizes.
    pub fn overlap(&self, k: usize, other: &Modes, j: usize) -> f64 {
        let (a, b) = (self.shape(k), other.shape(j));
        assert_eq!(a.len(), b.len(), "modes of structures of different sizes");
        a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>().abs()
    }

    /// `|cos|` between mode `k` and a direction in `3n` space.
    ///
    /// The comparison [`Structure::direction_to`](crate::structure::Structure::direction_to) exists for: how much of what a protein actually
    /// did is along one mode the model says it has. One is the whole motion, and the scale to
    /// read it against is what a *random* direction would give, which is `1/√3n` — `0.039` for a
    /// 214-residue protein.
    ///
    /// Absolute, like [`Modes::overlap`], and for the same reason: a mode's sign is arbitrary, so
    /// a signed comparison would call a perfect match `−1`.
    ///
    /// # Panics
    ///
    /// If `k` is out of range or `direction` is not `3n` long.
    pub fn overlap_with(&self, k: usize, direction: &[f64]) -> f64 {
        let shape = self.shape(k);
        assert_eq!(shape.len(), direction.len(), "a direction of the same size");
        shape
            .iter()
            .zip(direction)
            .map(|(a, b)| a * b)
            .sum::<f64>()
            .abs()
    }

    /// `⟨Δr²⟩` for every node, in `m²`.
    ///
    /// # Panics
    ///
    /// If `temperature` is not positive.
    pub fn fluctuations(&self, temperature: Temperature) -> Vec<f64> {
        let kt = BOLTZMANN.to_si() * temperature.to_si();
        assert!(kt > 0.0, "a temperature of {} K", temperature.to_si());
        let mut out = vec![0.0; self.nodes];
        for k in self.rigid..self.spectrum.order() {
            let lambda = self.spectrum.values()[k];
            let u = self.spectrum.vector(k);
            for (i, slot) in out.iter_mut().enumerate() {
                let s =
                    u[3 * i] * u[3 * i] + u[3 * i + 1] * u[3 * i + 1] + u[3 * i + 2] * u[3 * i + 2];
                *slot += kt * s / lambda;
            }
        }
        out
    }

    /// The temperature factor each node would have in a PDB file, in `m²`.
    ///
    /// `B = 8π²⟨Δr²⟩/3`, the same relation [`Structure::experimental_fluctuations`](crate::structure::Structure::experimental_fluctuations) inverts, so
    /// the two are directly comparable. A PDB file writes these in `Å²`; this is SI like
    /// everything else here.
    ///
    /// # Panics
    ///
    /// If `temperature` is not positive.
    pub fn b_factors(&self, temperature: Temperature) -> Vec<f64> {
        let factor = 8.0 * std::f64::consts::PI * std::f64::consts::PI / 3.0;
        self.fluctuations(temperature)
            .into_iter()
            .map(|x| factor * x)
            .collect()
    }

    /// The spectrum underneath, rigid modes and all.
    pub fn spectrum(&self) -> &Spectrum {
        &self.spectrum
    }
}

/// Pearson's correlation between two series.
///
/// This is how an elastic network model is judged against a crystal: not by whether the
/// fluctuations come out the right *size* — the stiffness sets that and nothing predicts it —
/// but by whether the same residues are the mobile ones.
///
/// **Read the result beside a baseline and not beside a remembered range.** What this crate
/// measures on its own four deposited structures is `0.659`, `0.571`, `0.395` and `−0.404`; what
/// the distance from the centroid measures on the same four is `0.699`, `0.804`, `0.402` and
/// `−0.301`. A correlation here is a weak instrument, and
/// `tests/a_protein_a_crystallographer_measured.rs` is where that is written down in numbers.
///
/// `None` if the two are different lengths, empty, or if either is constant — a constant series
/// has no correlation with anything, and the division that produces one would return `NaN`, which
/// compares false against every threshold and so passes a test written with `<`.
pub fn correlation(a: &[f64], b: &[f64]) -> Option<f64> {
    if a.len() != b.len() || a.is_empty() {
        return None;
    }
    let n = a.len() as f64;
    let (mean_a, mean_b) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let mut cov = 0.0;
    let mut var_a = 0.0;
    let mut var_b = 0.0;
    for (x, y) in a.iter().zip(b) {
        cov += (x - mean_a) * (y - mean_b);
        var_a += (x - mean_a) * (x - mean_a);
        var_b += (y - mean_b) * (y - mean_b);
    }
    if var_a == 0.0 || var_b == 0.0 {
        return None;
    }
    Some(cov / (var_a * var_b).sqrt())
}
