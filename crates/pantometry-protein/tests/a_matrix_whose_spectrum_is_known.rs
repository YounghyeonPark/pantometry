//! The eigensolver against matrices whose eigenvalues are known in closed form.
//!
//! Convention 1 of this workspace is to check against a closed form and never against another
//! implementation, and a dense symmetric eigensolver is unusually well served by it: three
//! families of matrix have their whole spectrum written down, and two invariants tie the
//! eigenvalues to the input without reference to any of them.
//!
//! | | what it pins |
//! | --- | --- |
//! | The free chain, `4 sin²(jπ/2N)` | The spectrum, **including the one exact zero** — the 1D shadow of the six a protein's Hessian has |
//! | The Dirichlet chain, `4 sin²(jπ/2(N+1))` | The spectrum *and* the eigenvectors, `sin(ijπ/(N+1))`, which no invariant below can see |
//! | A constructed `QΛQᵀ` | A spectrum chosen to be nasty: a triple root, a zero, and two values a part in `10⁸` apart |
//! | `Σλ = tr A`, `Σλ² = ‖A‖²_F` | Every eigenvalue of an arbitrary matrix, with no closed form needed |
//! | `‖Av − λv‖` | The **pairing**. Every line above is satisfied by a solver that shuffles the eigenvectors |
//!
//! # The one that is not a closed form, and why it is here anyway
//!
//! [`the_residual_and_the_orthogonality_are_at_the_rounding_level`] compares against `ε‖A‖`. That
//! is not a tolerance somebody chose: an orthogonal similarity in floating point has backward
//! error `O(ε‖A‖)` per rotation, so a residual at that level is the arithmetic's floor and a
//! residual above it is a bug. The constant in front is the thing to argue about, and
//! [`rotation_floor`] derives it from how many rotations actually ran rather than picking a
//! number — which it has to, because the number picked first was wrong twice.
//!
//! # What the eigenvectors are compared up to
//!
//! Up to sign, and that is not a convenience. `v` and `−v` are the same mode, and there is no
//! rule that can pick between them continuously: [`Spectrum::vector`] makes the
//! largest-magnitude component positive, which is a *tie* whenever the structure is symmetric,
//! and a symmetric structure is the normal case here rather than the exotic one. Measured on the
//! twenty-bead chain below: every eigenvector matches the closed form to `1e-15`, and **five of
//! the twenty come back with the opposite sign** because the two largest components are equal in
//! magnitude and rounding decides which one the rule sees.
//!
//! So the first version of that test failed at `n = 20` with an error of `0.615`, which is not a
//! rounding level and is not a defect either. The consequence for everything built on this
//! crate is the line worth keeping: **anything read off a mode has to be invariant under
//! `v → −v`** — an overlap between two modes is `|cos|`, never `cos`.

use pantometry_protein::spectrum::Symmetric;
use std::f64::consts::PI;

/// The rounding floor for a quantity accumulated over every rotation the solver performed.
///
/// A sweep applies a rotation to each of the `n(n−1)/2` pairs, so each *column* is touched
/// `n−1` times per sweep and each touch is one rounding. `sweeps · (n−1) · ε` is that count,
/// and it is the bound every measurement in this file sits one to three orders of magnitude
/// below — the numbers are printed so that a change which quietly eats into that margin is
/// visible rather than merely still passing.
///
/// The first version of this file wrote `n · ε`, which traced to nothing and failed at `n = 6`
/// and `n = 25` by factors of 1.3 and 1.03. A bound that is nearly right by luck is the shape a
/// tolerance takes just before somebody multiplies it by ten.
fn rotation_floor(n: usize, sweeps: usize) -> f64 {
    rotations(n, sweeps) as f64 * f64::EPSILON
}

/// How many rotations touch one column: `n−1` per sweep, and at least one, so that a matrix
/// that was already diagonal still gets a floor rather than zero.
fn rotations(n: usize, sweeps: usize) -> usize {
    sweeps.max(1) * n.max(2).saturating_sub(1)
}

/// The floor for a comparison against a closed form that has to be **evaluated**, which
/// [`rotation_floor`] alone is too tight for.
///
/// `4 sin²(jπ/2N)` is exact as mathematics and is not exact as arithmetic. `sin` is not
/// correctly rounded and is not the same function on two platforms: a libm is normally within one
/// ulp and normally no better, squaring doubles that relative error, and rounding the square adds
/// a half. Three ulps of the value, and the value is at most the largest eigenvalue, which is at
/// most `‖A‖_F`.
///
/// **This is the only place in the crate where a platform difference can appear.** Everything the
/// library itself computes uses `+ - * / sqrt`, every one of which IEEE-754 defines exactly, so
/// its answers are bit-identical everywhere and convention 3 holds. What varies is the reference.
///
/// Measured, and the reason this exists: at `n = 2` this machine got an error of exactly `0` and
/// CI's Windows runner got `8.88e-16` — two ulps of 2.0 — against a rotations-only floor of
/// `4.44e-16`. A bound that holds on the machine that wrote it is not a bound.
fn evaluated_floor(n: usize, sweeps: usize) -> f64 {
    (rotations(n, sweeps) + 3) as f64 * f64::EPSILON
}

/// The floor for an **invariant that is a sum**, which [`rotation_floor`] alone is too tight for.
///
/// `Σλ²` and `‖A‖²_F` are the same exact number computed two ways, and the two ways add up
/// different numbers of terms: `n` against `n²`. Even with an exact eigendecomposition the two
/// results differ by the summations' own rounding — one for each product and one for each
/// addition — so the floor is that count as well as the rotations.
///
/// Measured, and the reason this exists separately: at `n = 2` and one sweep the Frobenius
/// comparison came out at `3.47e-16` against a rotations-only floor of `2.22e-16`. One rotation
/// cannot be responsible for a `1.6 ε` disagreement; four products and four additions can.
fn sum_floor(terms: usize, n: usize, sweeps: usize) -> f64 {
    (2 * terms + rotations(n, sweeps)) as f64 * f64::EPSILON
}

/// `Σ λ` and `Σ λ²` against the trace and the Frobenius norm, as relative errors, with the
/// sweep count the floor needs.
fn invariants(m: &Symmetric) -> (f64, f64, usize) {
    let s = m.eigen();
    let sum: f64 = s.values().iter().sum();
    let squares: f64 = s.values().iter().map(|x| x * x).sum();
    let scale = m.frobenius_squared().sqrt();
    (
        (sum - m.trace()).abs() / scale,
        (squares - m.frobenius_squared()).abs() / m.frobenius_squared(),
        s.sweeps(),
    )
}

/// A matrix that is reproducible without an RNG and has no structure the solver could exploit.
///
/// `sin` of a pair of indices: bounded, irrational, and symmetric by construction.
fn a_matrix_with_no_pattern(n: usize) -> Symmetric {
    let mut m = Symmetric::zeros(n);
    for i in 0..n {
        for j in i..n {
            let (x, y) = (i as f64, j as f64);
            m.set(i, j, (1.0 + x * 7.3 + y * 2.9 + x * y * 0.11).sin());
        }
    }
    m
}

/// **A chain of beads with free ends: the eigenvalues are `4 sin²(jπ/2N)` exactly.**
///
/// This is the elastic network model in one dimension — `N` nodes, a unit spring between
/// neighbours, nothing holding the ends — and it is the closed form that matters most here,
/// because of what the `j = 0` term is. `4 sin²(0) = 0`: the chain can slide, and that mode is an
/// **exact** null vector of the matrix, not a small eigenvalue. A protein's Hessian has six of
/// those for the same reason, and telling them from the lowest real mode is the whole of a normal
/// mode analysis. A solver that returned `1e-9` where the closed form says `0` would make the six
/// rigid-body modes indistinguishable from a floppy hinge.
#[test]
fn a_free_chain_has_the_spectrum_of_a_sine_squared_and_one_exact_zero() {
    for n in [2usize, 3, 5, 8, 13, 40] {
        let mut m = Symmetric::zeros(n);
        for i in 0..(n - 1) {
            m.add(i, i, 1.0);
            m.add(i + 1, i + 1, 1.0);
            m.add(i, i + 1, -1.0);
        }
        let s = m.eigen();
        assert!(s.converged(), "n = {n} did not converge");

        let scale = m.frobenius_squared().sqrt();
        // Weyl for the solver's half, plus the reference's own evaluation. The sliding mode
        // below keeps the tighter of the two, because zero does not go through `sin`.
        let exact_zero = rotation_floor(n, s.sweeps()) * scale;
        let floor = evaluated_floor(n, s.sweeps()) * scale;
        let mut worst: f64 = 0.0;
        for j in 0..n {
            let exact = 4.0 * (j as f64 * PI / (2.0 * n as f64)).sin().powi(2);
            worst = worst.max((s.values()[j] - exact).abs());
        }
        assert!(
            worst <= floor,
            "n = {n}: worst eigenvalue error {worst:e} against 4 sin^2(j pi / 2N), and the \
             arithmetic's floor is {floor:e}"
        );

        // The sliding mode, separately, because the line above measures it against a *scale* and
        // this one says it is zero at the level a sum of `n` unit terms can be zero at.
        assert!(
            s.values()[0].abs() <= exact_zero,
            "n = {n}: the sliding mode came back at {:e}, not zero, over {exact_zero:e}",
            s.values()[0]
        );
        println!(
            "  free chain n = {n:3}: worst {worst:.3e}/{floor:.3e}, sliding mode {:.3e}/{exact_zero:.3e}",
            s.values()[0]
        );
    }
}

/// **A chain held at both ends: `4 sin²(jπ/2(N+1))`, and the eigenvectors are sampled sines.**
///
/// The free chain above pins the eigenvalues and says nothing about which vector goes with which.
/// This one does, because the eigenvectors of the Dirichlet tridiagonal are known too —
/// `v_j[i] ∝ sin((i+1)(j+1)π/(N+1))` — so a solver that produced the right spectrum attached to
/// the wrong vectors fails here and passes everything else in this file except the residual.
#[test]
fn a_held_chain_matches_its_eigenvectors_as_well_as_its_eigenvalues() {
    for n in [2usize, 3, 7, 20] {
        let mut m = Symmetric::zeros(n);
        for i in 0..n {
            m.set(i, i, 2.0);
            if i + 1 < n {
                m.set(i, i + 1, -1.0);
            }
        }
        let s = m.eigen();
        let scale = m.frobenius_squared().sqrt();
        let floor = evaluated_floor(n, s.sweeps()) * scale;

        let mut worst_value: f64 = 0.0;
        let mut worst_vector: f64 = 0.0;
        let mut flipped = 0usize;
        for j in 0..n {
            let exact = 4.0 * ((j + 1) as f64 * PI / (2.0 * (n + 1) as f64)).sin().powi(2);
            worst_value = worst_value.max((s.values()[j] - exact).abs());

            // The closed-form vector, normalised. Its sign is not fixed and cannot be — see the
            // module documentation — so the comparison takes the better of the two, and counts
            // how often that is the flipped one so the number is on the record rather than
            // hidden by the `min`.
            let mut exact_v: Vec<f64> = (0..n)
                .map(|i| ((i + 1) as f64 * (j + 1) as f64 * PI / (n + 1) as f64).sin())
                .collect();
            let norm: f64 = exact_v.iter().map(|x| x * x).sum::<f64>().sqrt();
            for x in exact_v.iter_mut() {
                *x /= norm;
            }
            let agreeing = s
                .vector(j)
                .iter()
                .zip(&exact_v)
                .map(|(g, w)| (g - w).abs())
                .fold(0.0f64, f64::max);
            let opposed = s
                .vector(j)
                .iter()
                .zip(&exact_v)
                .map(|(g, w)| (g + w).abs())
                .fold(0.0f64, f64::max);
            if opposed < agreeing {
                flipped += 1;
            }
            worst_vector = worst_vector.max(agreeing.min(opposed));
        }
        assert!(
            worst_value <= floor,
            "n = {n}: eigenvalue error {worst_value:e} over {floor:e}"
        );
        // Not an assertion about *how many* flip — that is rounding's business and would change
        // with the compiler. It is an assertion that the comparison above is doing work: if no
        // vector ever disagreed in sign, taking the minimum of the two would be dead code and a
        // solver that negated every vector would pass.
        if n == 20 {
            assert!(
                flipped > 0,
                "no sign disagreed, so the two-sided comparison is untested"
            );
        }
        // An eigenvector's own error is larger than an eigenvalue's whenever two eigenvalues are
        // close, by `ε‖A‖ / gap` — and this chain's gap narrows as `1/n²`, which is the `n²` here.
        let vector_floor = (n * n) as f64 * f64::EPSILON * scale;
        assert!(
            worst_vector <= vector_floor,
            "n = {n}: eigenvector error {worst_vector:e} over the gap-limited floor {vector_floor:e}"
        );
        println!(
            "  held chain n = {n:3}: value {worst_value:.3e}, vector {worst_vector:.3e}, {flipped}/{n} opposed"
        );
    }
}

/// **A spectrum chosen to be difficult, planted in a matrix and read back out.**
///
/// `Q Λ Qᵀ` with `Q` a product of rotations at angles nothing here can simplify. This is a closed
/// form in the sense that matters — `Λ` is known exactly because it was *chosen*, not because a
/// second routine computed it — and it is the only check here that can carry a spectrum the
/// physical families above cannot produce:
///
/// - **a triple root**, which is what a degenerate subspace does to a solver that assumes
///   distinct eigenvalues,
/// - **an exact zero** sitting next to nonzero values rather than at the end,
/// - **two values a part in `10⁸` apart**, which is what a nearly-degenerate pair does, and
/// - **a negative eigenvalue**, which a Hessian at a real minimum never has and a Hessian at a
///   saddle does.
#[test]
fn a_planted_spectrum_comes_back_including_a_triple_root_and_a_zero() {
    let planted = [-3.5, 0.0, 1.0, 1.0, 1.0, 2.0, 2.000_000_02, 9.75];
    let n = planted.len();

    let mut q = vec![0.0; n * n];
    for i in 0..n {
        q[i * n + i] = 1.0;
    }
    // Rotate every pair by an angle that is not a fraction of pi, so no entry of `Q` is a number
    // the arithmetic can get exactly right by luck.
    for p in 0..n {
        for r in (p + 1)..n {
            let theta = 0.37 + 0.19 * p as f64 + 0.11 * r as f64;
            let (c, s) = (theta.cos(), theta.sin());
            for k in 0..n {
                let (a, b) = (q[k * n + p], q[k * n + r]);
                q[k * n + p] = c * a - s * b;
                q[k * n + r] = s * a + c * b;
            }
        }
    }
    // A = Q diag(planted) Q^T, symmetrised by construction: the expression for (i, j) and (j, i)
    // is the same sum in the same order.
    let mut m = Symmetric::zeros(n);
    for i in 0..n {
        for j in i..n {
            let entry: f64 = (0..n)
                .map(|k| q[i * n + k] * planted[k] * q[j * n + k])
                .sum();
            m.set(i, j, entry);
        }
    }

    let s = m.eigen();
    assert!(s.converged());
    // Rotations only. `A` was *built* in floating point -- each entry is a sum of `n` triple
    // products, so it carries about `2n ε` of relative error and is not quite the matrix whose
    // eigenvalues are exactly `planted`; Weyl turns that into another `2n ε ‖A‖` of eigenvalue
    // error. At `n = 8` that is `16 ε ‖A‖` against a floor of about `49 ε ‖A‖`, and the
    // measured error is `1.42e-14` against a floor of `1.18e-13`, so the term is inside the
    // margin and the floor stays where it is. Written down because the sine tests above had the
    // same omission and it took a Windows runner to find it.
    let floor = rotation_floor(n, s.sweeps()) * m.frobenius_squared().sqrt();
    for (k, &want) in planted.iter().enumerate() {
        let got = s.values()[k];
        assert!(
            (got - want).abs() <= floor,
            "planted {want} came back as {got} ({:e} over a floor of {floor:e})",
            (got - want).abs()
        );
    }
    // The triple root's *basis* is arbitrary, so nothing may be asserted about the individual
    // vectors — only that they are an orthonormal basis of the right subspace, which the
    // residual below covers and the orthogonality test covers.
    assert!(
        s.orthogonality_error() <= floor,
        "{:e}",
        s.orthogonality_error()
    );
    assert!(s.residual(&m) <= floor, "{:e}", s.residual(&m));
    println!(
        "  planted: worst {:.3e}, residual {:.3e}, orthogonality {:.3e}",
        planted
            .iter()
            .enumerate()
            .map(|(k, w)| (s.values()[k] - w).abs())
            .fold(0.0, f64::max),
        s.residual(&m),
        s.orthogonality_error()
    );
}

/// **`Σλ = tr A` and `Σλ² = ‖A‖²_F`, on a matrix with no closed form at all.**
///
/// Both are exact for any real symmetric matrix, and they hold for matrices the families above
/// cannot reach. They are weak on their own — a solver that returned the diagonal of a matrix it
/// never rotated would satisfy the first — which is why they are here beside the residual rather
/// than instead of it.
#[test]
fn the_eigenvalues_sum_to_the_trace_and_square_to_the_frobenius_norm() {
    for n in [1usize, 2, 6, 25, 60] {
        let m = a_matrix_with_no_pattern(n);
        let (trace_error, frobenius_error, sweeps) = invariants(&m);
        // The trace adds `n` terms on each side; the Frobenius norm adds `n²` on one of them.
        let trace_floor = sum_floor(n, n, sweeps);
        let frobenius_floor = sum_floor(n * n, n, sweeps);
        assert!(
            trace_error <= trace_floor,
            "n = {n}: trace off by {trace_error:e} over {trace_floor:e}"
        );
        assert!(
            frobenius_error <= frobenius_floor,
            "n = {n}: Frobenius off by {frobenius_error:e} over {frobenius_floor:e}"
        );
        println!(
            "  n = {n:3}: trace {trace_error:.3e}/{trace_floor:.3e}, Frobenius {frobenius_error:.3e}/{frobenius_floor:.3e}"
        );
    }
}

/// **Every eigenpair satisfies `Av = λv`, and the vectors are orthonormal.**
///
/// This is the one check in this file that sees the *pairing*. A spectrum whose eigenvalues are
/// all correct and whose vectors have been permuted passes the closed forms on the values, passes
/// the trace, passes the Frobenius norm, and fails here.
#[test]
fn the_residual_and_the_orthogonality_are_at_the_rounding_level() {
    for n in [1usize, 2, 6, 25, 60] {
        let m = a_matrix_with_no_pattern(n);
        let s = m.eigen();
        assert!(s.converged(), "n = {n}");
        let relative = rotation_floor(n, s.sweeps());
        let floor = relative * m.frobenius_squared().sqrt();
        let residual = s.residual(&m);
        let orthogonality = s.orthogonality_error();
        assert!(
            residual <= floor,
            "n = {n}: residual {residual:e} over {floor:e}"
        );
        assert!(
            orthogonality <= relative,
            "n = {n}: orthogonality {orthogonality:e} over {relative:e}"
        );
        println!(
            "  n = {n:3}: residual {residual:.3e} (floor {floor:.3e}), orthogonality \
             {orthogonality:.3e} (floor {relative:.3e}), {} sweeps",
            s.sweeps()
        );
    }
}

/// **A diagonal matrix is already solved, and a zero matrix is not a special case.**
///
/// Both are edge cases a sweep loop can get wrong in the same way — by doing a rotation with a
/// zero denominator — and both have their answer written on them.
#[test]
fn a_matrix_that_is_already_diagonal_is_returned_untouched() {
    let mut m = Symmetric::zeros(4);
    for (i, v) in [3.0, -1.0, 0.0, 7.5].into_iter().enumerate() {
        m.set(i, i, v);
    }
    let s = m.eigen();
    assert_eq!(s.sweeps(), 0, "a diagonal matrix needs no sweep");
    assert_eq!(s.values(), &[-1.0, 0.0, 3.0, 7.5]);
    assert_eq!(s.vector(0), &[0.0, 1.0, 0.0, 0.0]);

    let zero = Symmetric::zeros(3);
    let s = zero.eigen();
    assert!(s.converged());
    assert_eq!(s.values(), &[0.0, 0.0, 0.0]);
    assert_eq!(s.orthogonality_error(), 0.0);
}

/// **The sign rule is a rule, and it is the only thing in this file that tests it.**
///
/// Every closed-form comparison above is up to sign, because a mode's sign is arbitrary; that is
/// correct physics and it leaves [`Spectrum::vector`]'s stated rule — the largest-magnitude
/// component is positive — checked by nothing. Deleting the rule passed the whole of the rest of
/// this file, measured.
///
/// What the rule buys is small and worth being exact about: the same matrix gives the same signs
/// on every platform and under every optimisation level. What it does not buy is any stability
/// against a change to the matrix, and the sign-flip counts printed above are that in numbers.
#[test]
fn the_largest_component_of_every_eigenvector_is_positive() {
    let mut checked = 0;
    for m in [
        a_matrix_with_no_pattern(3),
        a_matrix_with_no_pattern(12),
        a_matrix_with_no_pattern(31),
    ] {
        let s = m.eigen();
        for k in 0..s.order() {
            let v = s.vector(k);
            let mut big = 0;
            for i in 1..v.len() {
                if v[i].abs() > v[big].abs() {
                    big = i;
                }
            }
            assert!(
                v[big] >= 0.0,
                "mode {k} of {} has its largest component, {}, negative",
                s.order(),
                v[big]
            );
            checked += 1;
        }
    }
    // A loop that measured nothing would satisfy the assertion above for free.
    assert_eq!(checked, 3 + 12 + 31);
    println!("  {checked} eigenvectors, all with the largest component positive");
}

/// **A budget that runs out says so.**
///
/// The failure this guards is the one that reads as success: a partly diagonalised matrix has a
/// diagonal, the diagonal is a list of numbers, and nothing about the list says it is not the
/// spectrum. One sweep on a matrix that needs eight has to come back `converged() == false` with
/// an off-diagonal norm that is visibly not small.
#[test]
fn a_spectrum_that_ran_out_of_sweeps_reports_it() {
    let m = a_matrix_with_no_pattern(30);
    let full = m.eigen();
    assert!(
        full.converged() && full.sweeps() >= 3,
        "{} sweeps",
        full.sweeps()
    );

    let stopped = m.eigen_within(1);
    assert!(
        !stopped.converged(),
        "one sweep on a matrix that took {}",
        full.sweeps()
    );
    assert_eq!(stopped.sweeps(), 1);
    assert!(
        stopped.off_diagonal_norm() > 1e-3 * m.frobenius_squared().sqrt(),
        "an unconverged run left {:e}, which is small enough to be mistaken for an answer",
        stopped.off_diagonal_norm()
    );
    println!(
        "  {} sweeps to converge; after 1, off-diagonal {:.3e} of {:.3e}",
        full.sweeps(),
        stopped.off_diagonal_norm(),
        m.frobenius_squared().sqrt()
    );
}

/// **The same matrix gives the same bits, and `add` on the diagonal adds once.**
///
/// Convention 3 is that nothing is random and results are bit-for-bit reproducible; this solver
/// consults no clock and no generator, so the claim reduces to the arithmetic being in a fixed
/// order, which is what the first half checks. The second half is the assembly mistake
/// [`Symmetric::add`]'s documentation names: a version that wrote both triangles unconditionally
/// would double every diagonal entry, and the matrix would still be symmetric.
#[test]
fn the_same_matrix_gives_the_same_answer_and_the_diagonal_is_added_once() {
    let m = a_matrix_with_no_pattern(12);
    assert_eq!(m.eigen(), m.eigen());

    let mut d = Symmetric::zeros(2);
    d.add(0, 0, 1.0);
    d.add(0, 1, 2.0);
    assert_eq!(d.get(0, 0), 1.0, "the diagonal took the value twice");
    assert_eq!(d.get(0, 1), 2.0);
    assert_eq!(d.get(1, 0), 2.0);

    assert!(Symmetric::from_rows(2, vec![1.0, 2.0, 2.0, 3.0]).is_some());
    assert!(
        Symmetric::from_rows(2, vec![1.0, 2.0, 2.000_000_1, 3.0]).is_none(),
        "an asymmetric matrix was accepted"
    );
    assert!(Symmetric::from_rows(2, vec![1.0, 2.0, 2.0]).is_none());
}
