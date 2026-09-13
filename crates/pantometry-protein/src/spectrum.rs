//! Every eigenvalue and eigenvector of a real symmetric matrix.
//!
//! A normal-mode analysis is an eigenproblem and nothing else: the Hessian of a harmonic
//! potential is real and symmetric, its eigenvalues are the squared frequencies, and its
//! eigenvectors are the directions the structure moves in. Everything this crate says about a
//! protein is a statement about that spectrum, so the spectrum is the first thing that has to be
//! right.
//!
//! # Why it is written here rather than pulled in
//!
//! This workspace resolves twelve external crates and `deny.toml` gates every one; the same
//! crates go to `wasm32` and to Rust 1.78. A dense symmetric eigensolver is two hundred lines of
//! arithmetic with no dependencies and a closed form to check every part of it against, which is
//! a much smaller thing to carry than a linear-algebra stack.
//!
//! # The method: cyclic Jacobi
//!
//! Repeatedly pick an off-diagonal pair `(p, q)` and apply the plane rotation that makes it
//! exactly zero. Each rotation is orthogonal, so the eigenvalues never move; each one reduces the
//! sum of the squared off-diagonal entries by `2 a_pq²`, so the matrix walks to diagonal and the
//! accumulated rotations are the eigenvectors.
//!
//! It is chosen over the faster tridiagonalise-then-QL route for two reasons that matter here and
//! would not in general:
//!
//! - **It is short enough to read.** Convergence is monotone and needs no shift strategy, no
//!   deflation and no special case for a repeated eigenvalue — and a protein's Hessian is
//!   degenerate by construction, with six exact zeros.
//! - **Orthogonality is structural.** The eigenvectors are a product of orthogonal rotations, so
//!   `VᵀV = I` to rounding whatever the spectrum looks like. That matters because the six
//!   rigid-body modes have to be *projected out*, and projecting with a basis that is not quite
//!   orthonormal leaves a residue that looks like physics.
//!
//! The cost is `O(n³)` per sweep and about ten sweeps; [`Symmetric::eigen`] says what that is in
//! seconds.
//!
//! # No tolerance to earn
//!
//! The off-diagonal norm falls monotonically in exact arithmetic and plateaus in floating point at
//! the rounding level of the rotations themselves. So the stopping rule is not a number somebody
//! picked: a sweep that fails to reduce the off-diagonal norm is the last useful sweep, and
//! [`Spectrum::off_diagonal_norm`] reports where it stopped rather than asserting that this was
//! small enough. [`Spectrum::converged`] is false only when the sweep budget ran out while
//! progress was still being made.

/// How many sweeps [`Symmetric::eigen`] allows before it reports failure.
///
/// Jacobi converges quadratically once the off-diagonal entries are small, and the matrices here
/// reach the floating-point plateau in about ten sweeps. Sixty is not a tuned number — it is far
/// enough above that for the budget never to be what stops a well-posed problem, so that
/// [`Spectrum::converged`] coming back false means something went wrong rather than that the
/// caller was impatient.
const SWEEP_BUDGET: usize = 60;

/// A real symmetric matrix, stored whole.
///
/// Both triangles are kept and kept equal: [`Symmetric::set`] and [`Symmetric::add`] write both,
/// and the only constructor that could produce an asymmetric one, [`Symmetric::from_rows`],
/// refuses.
#[derive(Clone, Debug, PartialEq)]
pub struct Symmetric {
    n: usize,
    a: Vec<f64>,
}

impl Symmetric {
    /// An `n × n` matrix of zeros.
    pub fn zeros(n: usize) -> Symmetric {
        Symmetric {
            n,
            a: vec![0.0; n * n],
        }
    }

    /// The matrix whose rows are `rows`, laid out row by row.
    ///
    /// `None` unless `rows` is exactly `n²` long and **exactly** symmetric. Exactly, rather than
    /// within a tolerance: a caller that built the matrix from a symmetric expression has
    /// bit-identical pairs, and one that did not has a bug, and this is the only place it can be
    /// caught.
    pub fn from_rows(n: usize, rows: Vec<f64>) -> Option<Symmetric> {
        if rows.len() != n * n {
            return None;
        }
        for i in 0..n {
            for j in (i + 1)..n {
                if rows[i * n + j] != rows[j * n + i] {
                    return None;
                }
            }
        }
        Some(Symmetric { n, a: rows })
    }

    /// The order `n`.
    pub fn order(&self) -> usize {
        self.n
    }

    /// The entry at `(i, j)`.
    ///
    /// # Panics
    ///
    /// If either index is at or past [`Symmetric::order`].
    pub fn get(&self, i: usize, j: usize) -> f64 {
        assert!(i < self.n && j < self.n, "({i}, {j}) is outside {}", self.n);
        self.a[i * self.n + j]
    }

    /// Writes `value` at `(i, j)` **and** at `(j, i)`.
    ///
    /// # Panics
    ///
    /// If either index is at or past [`Symmetric::order`].
    pub fn set(&mut self, i: usize, j: usize, value: f64) {
        assert!(i < self.n && j < self.n, "({i}, {j}) is outside {}", self.n);
        self.a[i * self.n + j] = value;
        self.a[j * self.n + i] = value;
    }

    /// Adds `value` at `(i, j)` and at `(j, i)` — once in total on the diagonal, where those are
    /// the same entry.
    ///
    /// This is how a network Hessian is assembled: every spring contributes to four blocks, and a
    /// version of this that doubled the diagonal would produce a matrix that is still symmetric,
    /// still positive semi-definite, and wrong everywhere.
    ///
    /// # Panics
    ///
    /// If either index is at or past [`Symmetric::order`].
    pub fn add(&mut self, i: usize, j: usize, value: f64) {
        assert!(i < self.n && j < self.n, "({i}, {j}) is outside {}", self.n);
        self.a[i * self.n + j] += value;
        if i != j {
            self.a[j * self.n + i] += value;
        }
    }

    /// `Σ a_ii`, which is also `Σ λ_k` for any symmetric matrix.
    pub fn trace(&self) -> f64 {
        (0..self.n).map(|i| self.a[i * self.n + i]).sum()
    }

    /// `Σ a_ij²` over every entry, which is also `Σ λ_k²`.
    pub fn frobenius_squared(&self) -> f64 {
        self.a.iter().map(|x| x * x).sum()
    }

    /// `A x`.
    ///
    /// # Panics
    ///
    /// If `x` is not [`Symmetric::order`] long.
    pub fn multiply(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.n, "a vector of the matrix's own order");
        (0..self.n)
            .map(|i| (0..self.n).map(|j| self.a[i * self.n + j] * x[j]).sum())
            .collect()
    }

    /// Every eigenvalue and eigenvector, ascending by eigenvalue.
    ///
    /// # What it costs
    ///
    /// Measured on this machine, release, on a matrix with no exploitable structure:
    ///
    /// | `n` | | sweeps |
    /// | --- | --- | --- |
    /// | 150 | 0.05 s | 13 |
    /// | 300 | 0.31 s | 13 |
    /// | 600 | 3.7 s | 14 |
    /// | 900 | 13.1 s | 15 |
    ///
    /// The sweep count barely moves, so this is the `n³` the method promises, and the jump from
    /// 300 to 600 is steeper than cubic because the matrix stops fitting in cache.
    ///
    /// A protein of `N` residues gives `n = 3N`, which puts the useful range at a few hundred
    /// residues: ubiquitin's 76 is a tenth of a second and lysozyme's 129 is under one. Doubling
    /// the residue count costs eight times as much, so this is not the routine to reach for at a
    /// thousand residues without splitting the problem up first.
    pub fn eigen(&self) -> Spectrum {
        self.eigen_within(SWEEP_BUDGET)
    }

    /// [`Symmetric::eigen`] with the sweep budget named, for a caller that wants to stop it short.
    pub fn eigen_within(&self, max_sweeps: usize) -> Spectrum {
        let n = self.n;
        let mut a = self.a.clone();
        // Column-major, so that eigenvector `k` is the contiguous run at `k * n`.
        let mut v = vec![0.0; n * n];
        for k in 0..n {
            v[k * n + k] = 1.0;
        }

        let mut sweeps = 0;
        let converged;
        let mut previous = f64::INFINITY;
        let mut off = off_diagonal_norm(&a, n);
        loop {
            // A sweep that does not reduce the off-diagonal norm is the last useful sweep: the
            // fall is monotone in exact arithmetic, so a plateau is the rotations' own rounding.
            if off == 0.0 || off >= previous {
                converged = true;
                break;
            }
            if sweeps == max_sweeps {
                converged = false;
                break;
            }
            previous = off;
            sweeps += 1;

            for p in 0..n {
                for q in (p + 1)..n {
                    let apq = a[p * n + q];
                    if apq == 0.0 {
                        continue;
                    }
                    // The root of `t² + 2τt − 1 = 0` with the smaller magnitude, written so that
                    // neither branch subtracts two nearly equal numbers.
                    let tau = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
                    let root = (1.0 + tau * tau).sqrt();
                    let t = if tau >= 0.0 {
                        1.0 / (tau + root)
                    } else {
                        -1.0 / (root - tau)
                    };
                    let c = 1.0 / (1.0 + t * t).sqrt();
                    let s = t * c;

                    // A ← JᵀAJ, as the columns and then the rows.
                    for k in 0..n {
                        let (akp, akq) = (a[k * n + p], a[k * n + q]);
                        a[k * n + p] = c * akp - s * akq;
                        a[k * n + q] = s * akp + c * akq;
                    }
                    for k in 0..n {
                        let (apk, aqk) = (a[p * n + k], a[q * n + k]);
                        a[p * n + k] = c * apk - s * aqk;
                        a[q * n + k] = s * apk + c * aqk;
                    }
                    // `t` was chosen to annihilate this pair, so whatever the arithmetic left
                    // there is rounding. Keeping it would make the stopping rule above a
                    // measurement of noise rather than of progress.
                    a[p * n + q] = 0.0;
                    a[q * n + p] = 0.0;

                    // V ← VJ.
                    let (cp, cq) = (p * n, q * n);
                    for i in 0..n {
                        let (vip, viq) = (v[cp + i], v[cq + i]);
                        v[cp + i] = c * vip - s * viq;
                        v[cq + i] = s * vip + c * viq;
                    }
                }
            }
            off = off_diagonal_norm(&a, n);
        }

        let diagonal: Vec<f64> = (0..n).map(|i| a[i * n + i]).collect();
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&i, &j| diagonal[i].total_cmp(&diagonal[j]));
        let values: Vec<f64> = order.iter().map(|&i| diagonal[i]).collect();
        let mut vectors = vec![0.0; n * n];
        for (k, &from) in order.iter().enumerate() {
            vectors[k * n..(k + 1) * n].copy_from_slice(&v[from * n..(from + 1) * n]);
            fix_sign(&mut vectors[k * n..(k + 1) * n]);
        }

        Spectrum {
            n,
            values,
            vectors,
            sweeps,
            converged,
            off_diagonal: off,
        }
    }
}

/// `√(Σ_{i≠j} a_ij²)`, the quantity a Jacobi sweep drives to zero.
fn off_diagonal_norm(a: &[f64], n: usize) -> f64 {
    let mut sum = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            sum += a[i * n + j] * a[i * n + j];
        }
    }
    (2.0 * sum).sqrt()
}

/// Makes the largest-magnitude component positive, a tie going to the lower index.
///
/// An eigenvector's sign is arbitrary, and "arbitrary" in a library that promises bit-identical
/// results across platforms has to mean "fixed by a rule" rather than "whatever the arithmetic
/// happened to produce". This does not make a *degenerate* subspace's basis unique — nothing can,
/// because there is no preferred basis in one — so anything read off such a subspace has to be
/// invariant under a rotation inside it.
fn fix_sign(vector: &mut [f64]) {
    let mut best = 0;
    for i in 1..vector.len() {
        if vector[i].abs() > vector[best].abs() {
            best = i;
        }
    }
    if vector.get(best).is_some_and(|&x| x < 0.0) {
        for x in vector.iter_mut() {
            *x = -*x;
        }
    }
}

/// Every eigenvalue and eigenvector of a [`Symmetric`], ascending by eigenvalue.
#[derive(Clone, Debug, PartialEq)]
pub struct Spectrum {
    n: usize,
    values: Vec<f64>,
    vectors: Vec<f64>,
    sweeps: usize,
    converged: bool,
    off_diagonal: f64,
}

impl Spectrum {
    /// The order of the matrix this came from, which is also how many eigenpairs there are.
    pub fn order(&self) -> usize {
        self.n
    }

    /// The eigenvalues, ascending.
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// Eigenvector `k`, unit length, belonging to `values()[k]`.
    ///
    /// # Panics
    ///
    /// If `k` is at or past [`Spectrum::order`].
    pub fn vector(&self, k: usize) -> &[f64] {
        assert!(k < self.n, "mode {k} of {}", self.n);
        &self.vectors[k * self.n..(k + 1) * self.n]
    }

    /// How many sweeps it took.
    pub fn sweeps(&self) -> usize {
        self.sweeps
    }

    /// Whether the sweeps ran out before the off-diagonal norm stopped falling.
    ///
    /// False is not "the answer is slightly off" — it is "this is a partly diagonalised matrix and
    /// its diagonal is not the spectrum". A caller that ignores it gets numbers for any input.
    pub fn converged(&self) -> bool {
        self.converged
    }

    /// `√(Σ_{i≠j} a_ij²)` of what was left when it stopped.
    ///
    /// Reported rather than asserted against a tolerance, because what counts as small depends on
    /// the matrix: compare it against `√(`[`Symmetric::frobenius_squared`]`)`.
    pub fn off_diagonal_norm(&self) -> f64 {
        self.off_diagonal
    }

    /// The largest `‖A v_k − λ_k v_k‖` over every mode.
    ///
    /// The definition of an eigenpair, measured. This is the check that notices a correct
    /// eigenvalue carrying the wrong vector, which every invariant on the eigenvalues alone — the
    /// trace, the Frobenius norm, a closed-form spectrum — is blind to.
    ///
    /// # Panics
    ///
    /// If `matrix` is not the order this spectrum came from.
    pub fn residual(&self, matrix: &Symmetric) -> f64 {
        assert_eq!(matrix.order(), self.n, "the matrix this came from");
        let mut worst: f64 = 0.0;
        for k in 0..self.n {
            let v = self.vector(k);
            let av = matrix.multiply(v);
            let norm: f64 = av
                .iter()
                .zip(v)
                .map(|(a, x)| (a - self.values[k] * x).powi(2))
                .sum::<f64>()
                .sqrt();
            worst = worst.max(norm);
        }
        worst
    }

    /// The largest `|v_jᵀv_k − δ_jk|` over every pair.
    pub fn orthogonality_error(&self) -> f64 {
        let mut worst: f64 = 0.0;
        for j in 0..self.n {
            for k in j..self.n {
                let dot: f64 = self
                    .vector(j)
                    .iter()
                    .zip(self.vector(k))
                    .map(|(a, b)| a * b)
                    .sum();
                let expected = if j == k { 1.0 } else { 0.0 };
                worst = worst.max((dot - expected).abs());
            }
        }
        worst
    }
}
