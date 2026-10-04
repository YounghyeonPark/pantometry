//! **Benzene bound in T4 lysozyme L99A's cavity (PDB 181L): the binding energy in a rigid pocket,
//! and the ligand minimised there.** Step 2c-2.
//!
//! Every assertion is an identity or a bound derived beside it — the cross terms, the far field,
//! a change of frame, a frozen atom, a gradient, and what burying can do to an area and a Born
//! radius — and none is a comparison with experiment. The numbers that nothing bounds in
//! advance (the binding energy itself, the pose, the solvation terms, the cutoff's effect) are
//! printed, and the two measurements that take minutes are ignored by default:
//! `cargo test -p pantometry-forcefield --release --test benzene_in_its_pocket -- --ignored`.
//!
//! **Rounding, in every tolerance here**: each term — a bond, a pair — is the same value, bit for
//! bit, in the complex as in its fragment (same positions, same parameters, same charges, same
//! arithmetic), so the only rounding left in `E(complex) − E(pocket) − E(ligand)` is that of the
//! running sums and the two subtractions: a sum of `n` terms is within `n ε Σ|t|` of exact (first
//! order in ε; `n ε` is below 10⁻¹⁰ for every sum here), and each subtraction adds one ε of its
//! operands.

mod protein;

use pantometry_forcefield::energy::{coulomb, Pair};
use pantometry_forcefield::minimise::{Progress, KCAL_PER_MOL_ANGSTROM};
use pantometry_forcefield::solvation::GeneralizedBorn;
use pantometry_forcefield::uff::KCAL_PER_MOL;
use pantometry_forcefield::{Binding, Element, ForceField, Minimiser, Qeq, RigidMotion, Status};
use protein::*;
use std::sync::OnceLock;

const ANGSTROM: f64 = 1e-10;
const EPS: f64 = f64::EPSILON;

/// The pocket's cutoff, Å: 2c-1's, so its counts (18 residues, 322 atoms with benzene, +1) are
/// checked again here.
const POCKET: f64 = 6.0;

/// The minimiser's tolerance: the largest force on a ligand atom at most 1e-4 kcal mol⁻¹ Å⁻¹,
/// [`pantometry_forcefield::Molecule::DEFAULT_TOLERANCE`]'s.
const TOLERANCE: f64 = 1e-4 * KCAL_PER_MOL_ANGSTROM;

/// Points per atom for the surface: 1600, where the two-sphere check in
/// `the_generalized_born_against_closed_forms.rs` measured a worst error of 0.29%.
const POINTS: usize = 1600;

fn kcal(joules: f64) -> f64 {
    joules / KCAL_PER_MOL
}

/// The binding at the crystal pose, built once for every test in this file: two QEq solves, 7 s
/// unoptimised.
fn crystal() -> &'static Binding {
    static B: OnceLock<Binding> = OnceLock::new();
    B.get_or_init(|| {
        Binding::new(&system(), POCKET * ANGSTROM).unwrap_or_else(|e| panic!("181L: {e}"))
    })
}

/// The crystal binding with benzene minimised in the pocket, once.
fn minimised() -> &'static (Binding, Progress) {
    static M: OnceLock<(Binding, Progress)> = OnceLock::new();
    M.get_or_init(|| {
        let mut b = crystal().clone();
        let p = b.minimise_ligand(5000, TOLERANCE);
        (b, p)
    })
}

/// `|a − b|`.
fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    len(sub(a, b))
}

/// For one force field at `at`: each energy field's term count and `Σ|t|`, joules — bond, angle,
/// torsion, inversion, van der Waals, electrostatic — computed term by term through the public
/// per-term energies, independently of [`ForceField::evaluate`].
fn magnitudes(ff: &ForceField, at: &[[f64; 3]]) -> [(usize, f64); 6] {
    let mut m = [(0usize, 0.0f64); 6];
    let mut add = |k: usize, e: f64| {
        m[k].0 += 1;
        m[k].1 += e.abs();
    };
    for s in ff.stretches() {
        add(0, s.energy(distance(at[s.atoms[0]], at[s.atoms[1]])));
    }
    for b in ff.bends() {
        add(1, b.energy_at(at));
    }
    for t in ff.torsions() {
        add(2, t.energy_at(at));
    }
    for v in ff.inversions() {
        add(3, v.energy_at(at));
    }
    let q = ff.charges();
    for p in ff.pairs() {
        let [i, j] = p.atoms;
        let r = distance(at[i], at[j]);
        add(4, p.energy(r));
        if q[i] != 0.0 && q[j] != 0.0 {
            add(5, coulomb(q[i], q[j], r));
        }
    }
    m
}

/// `n ε Σ|t|` for one sum.
fn sum_rounding((n, s): (usize, f64)) -> f64 {
    n as f64 * EPS * s
}

/// Every protein–ligand pair, built in the test from the types and the charges with no
/// exclusion list: `(pocket atom, ligand atom, its pair, its Coulomb energy, r)`.
fn cross_pairs(b: &Binding) -> Vec<(usize, usize, Pair, f64, f64)> {
    let (at, q, t) = (b.positions(), b.charges(), b.types());
    let mut out = Vec::new();
    for p in 0..b.pocket_len() {
        for l in b.ligand_range() {
            let r = distance(at[p], at[l]);
            out.push((
                p,
                l,
                Pair::new([p, l], t[p], t[l]),
                coulomb(q[p], q[l], r),
                r,
            ));
        }
    }
    out
}

/// The allowance for `E(complex) − E(pocket) − E(ligand)` against the cross sum, field by field
/// in [`magnitudes`]' order and then the total: every running sum's `n ε Σ|t|`, the cross sum's,
/// and one ε of each operand of the two subtractions.
fn difference_allowance(b: &Binding) -> [f64; 7] {
    let n0 = b.pocket_len();
    let at = b.positions();
    let c = magnitudes(b.force_field(), at);
    let p = magnitudes(b.pocket_force_field(), &at[..n0]);
    let l = magnitudes(b.ligand_force_field(), &at[n0..]);
    let (ec, ep, el) = (
        b.force_field().energy(at),
        b.pocket_force_field().energy(&at[..n0]),
        b.ligand_force_field().energy(&at[n0..]),
    );
    let x = cross_pairs(b);
    let cross = [
        (
            x.len(),
            x.iter().map(|t| t.2.energy(t.4).abs()).sum::<f64>(),
        ),
        (x.len(), x.iter().map(|t| t.3.abs()).sum::<f64>()),
    ];
    let fields = |e: &pantometry_forcefield::Energy| {
        [
            e.bond,
            e.angle,
            e.torsion,
            e.inversion,
            e.van_der_waals,
            e.electrostatic,
        ]
    };
    let (fc, fp, fl) = (fields(&ec), fields(&ep), fields(&el));
    let mut out = [0.0; 7];
    for k in 0..6 {
        let cross_part = if k >= 4 {
            sum_rounding(cross[k - 4])
        } else {
            0.0
        };
        out[k] = sum_rounding(c[k])
            + sum_rounding(p[k])
            + sum_rounding(l[k])
            + cross_part
            + 2.0 * EPS * (fc[k].abs() + fp[k].abs() + fl[k].abs());
    }
    // The total runs over every term at once.
    let all = |m: &[(usize, f64); 6]| {
        m.iter()
            .fold((0usize, 0.0f64), |(n, s), &(a, b)| (n + a, s + b))
    };
    out[6] = sum_rounding(all(&c))
        + sum_rounding(all(&p))
        + sum_rounding(all(&l))
        + sum_rounding((2 * x.len(), cross[0].1 + cross[1].1))
        + 2.0 * EPS * (ec.total.abs() + ep.total.abs() + el.total.abs());
    out
}

/// **The pocket is 2c-1's**: 18 residues with a heavy atom within 6 Å of benzene, 310 protein
/// atoms and benzene's 12, total formal charge +1, benzene 0 — the counts
/// `qeq_converges_on_the_binding_pocket` asserts for the same rule.
#[test]
fn the_pocket_is_whole_residues_within_the_cutoff() {
    let b = crystal();
    assert_eq!(b.residues().len(), 18, "{:?}", b.residues());
    assert_eq!((b.pocket_len(), b.positions().len()), (310, 322));
    assert_eq!(b.formal_charges(), (1, 0));
    // Whole residues: every atom of every pocket residue is in, and nothing else of the protein.
    let s = system();
    let mut want: Vec<usize> = s
        .residues()
        .iter()
        .filter(|r| b.residues().contains(&r.label()))
        .flat_map(|r| r.atoms.clone())
        .collect();
    want.extend(s.atoms_in(pantometry_forcefield::Part::Ligand));
    assert_eq!(b.system_atoms(), &want[..]);
    // And each one in has a heavy atom within the cutoff, each one out has none.
    let c = s.component();
    let ligand = s.atoms_in(pantometry_forcefield::Part::Ligand);
    for r in s
        .residues()
        .iter()
        .filter(|r| r.part == pantometry_forcefield::Part::Protein)
    {
        let nearest = r
            .atoms
            .clone()
            .filter(|&i| c.atoms()[i].element != Element::H)
            .flat_map(|i| ligand.iter().map(move |&l| (i, l)))
            .map(|(i, l)| distance(c.atoms()[i].at, c.atoms()[l].at))
            .fold(f64::INFINITY, f64::min);
        assert_eq!(
            b.residues().contains(&r.label()),
            nearest < POCKET * ANGSTROM,
            "{} nearest heavy atom {} Å",
            r.label(),
            nearest / ANGSTROM
        );
    }
    // Each fragment's charges sum to its formal charge. The allowance is the test's own sum,
    // `n ε Σ|q|`, plus the solve's residual on eq 9's row, `Σ Q = Q_tot`: QEq's Gaussian
    // elimination with partial pivoting is backward stable, `|r| ≤ γ_{3n} ρ ‖C‖ ‖Q‖` (Higham,
    // *Accuracy and Stability of Numerical Algorithms*, 2nd ed., Theorem 9.4), and with every
    // entry of eq 9's row one the row's part is about `3 n ε ρ Σ|Q|`; ρ, the growth factor, is
    // taken as 8, and the residuals measured are printed beside it.
    let n0 = b.pocket_len();
    for (range, want) in [(0..n0, 1.0), (b.ligand_range(), 0.0)] {
        let q = &b.charges()[range];
        let n = q.len() as f64;
        let abs: f64 = q.iter().map(|x| x.abs()).sum();
        let allowance = n * EPS * abs + 3.0 * n * EPS * 8.0 * abs;
        let sum: f64 = q.iter().sum();
        assert!(
            (sum - want).abs() <= allowance,
            "charges sum to {sum}, not {want}: off {:e} > {allowance:e}",
            (sum - want).abs()
        );
        eprintln!(
            "{} charges sum to {want} within {:.1e} e, allowance {allowance:.1e}",
            q.len(),
            (sum - want).abs()
        );
    }
    eprintln!(
        "pocket: {:?}; QEq solves {:?}",
        b.residues(),
        b.qeq_solves()
    );
}

/// **The charges are QEq on each fragment alone**: the pocket's are what [`Qeq::equilibrate`]
/// gives for the pocket's own atoms, read here from the [`System`], at its own formal charge, and
/// the ligand's likewise — bit for bit, since it is the same arithmetic on the same numbers. A
/// binding that took one fragment's charges from the wrong atoms, or at the wrong total, would
/// put the wrong charges on every term downstream and still satisfy every identity below.
///
/// [`System`]: pantometry_forcefield::System
#[test]
fn the_charges_are_qeq_on_each_fragment_alone() {
    let b = crystal();
    let s = system();
    let atoms = s.component().atoms();
    let n0 = b.pocket_len();
    for range in [0..n0, b.ligand_range()] {
        let ids = &b.system_atoms()[range.clone()];
        let elements: Vec<Element> = ids.iter().map(|&i| atoms[i].element).collect();
        let at: Vec<[f64; 3]> = ids.iter().map(|&i| atoms[i].at).collect();
        let total: i32 = ids.iter().map(|&i| atoms[i].charge).sum();
        let q = Qeq::default()
            .equilibrate(&elements, &at, f64::from(total))
            .expect("QEq on the fragment");
        let got: Vec<u64> = b.charges()[range].iter().map(|x| x.to_bits()).collect();
        let want: Vec<u64> = q.charges.iter().map(|x| x.to_bits()).collect();
        assert_eq!(got, want, "{} atoms at total {total}", ids.len());
    }
}

/// **ΔE_bind in vacuum is exactly the protein–ligand cross terms.** The protein is rigid and has
/// no bond to the ligand, so `E(complex) − E(pocket) − E(ligand)` is zero in its bond, angle,
/// torsion and inversion fields and the cross van der Waals and Coulomb in the other two — to
/// the rounding of the sums ([`difference_allowance`]). The cross sum is checked twice: the
/// crate's, against one built here from the types and charges over **every** pocket × ligand
/// pair, with no exclusion list, so a cross pair the force field wrongly excluded is caught.
///
/// **What this does and does not check.** It holds the bookkeeping: that the three fragments are
/// built from the same terms, that their difference is the cross pairs and nothing else, and that
/// every pocket × ligand pair is covered. It does **not** check the pair formula: the test's own
/// sum uses the crate's [`Pair::energy`] and [`coulomb`], so a wrong Lennard-Jones or Coulomb
/// would pass here, and is held instead by `the_energy_terms_against_closed_forms.rs` (the
/// minimum `−D` at `x`, the zero at `x/2^(1/6)`, 332.0637 kcal/mol for two unit charges at 1 Å).
/// What it does catch is a cross sum that disagrees with the force fields, which a 1% error in
/// one of them alone is: of the far-field and identity tests, only this one sees it.
#[test]
fn the_vacuum_binding_energy_is_exactly_the_cross_terms() {
    let b = crystal();
    let i = b.interaction();
    let allowance = difference_allowance(b);
    let d = i.difference;
    let checks = [
        ("bond", d.bond, 0.0),
        ("angle", d.angle, 0.0),
        ("torsion", d.torsion, 0.0),
        ("inversion", d.inversion, 0.0),
        ("van der Waals", d.van_der_waals, i.van_der_waals),
        ("electrostatic", d.electrostatic, i.electrostatic),
        ("total", d.total, i.total()),
    ];
    let mut worst = 0.0f64;
    for ((name, got, want), tol) in checks.iter().zip(allowance) {
        let err = (got - want).abs();
        assert!(
            err <= tol,
            "{name}: difference {got:e} J against {want:e} J, off {err:e} > {tol:e}"
        );
        worst = worst.max(err / tol);
    }
    // The crate's cross sum against every pair, built here.
    let x = cross_pairs(b);
    assert_eq!(x.len(), b.pocket_len() * b.ligand_range().len());
    let vdw: f64 = x.iter().map(|t| t.2.energy(t.4)).sum();
    let elec: f64 = x.iter().map(|t| t.3).sum();
    let s_vdw: f64 = x.iter().map(|t| t.2.energy(t.4).abs()).sum();
    let s_elec: f64 = x.iter().map(|t| t.3.abs()).sum();
    let n = x.len() as f64;
    assert!(
        (vdw - i.van_der_waals).abs() <= 2.0 * n * EPS * s_vdw,
        "{vdw} against {}",
        i.van_der_waals
    );
    assert!(
        (elec - i.electrostatic).abs() <= 2.0 * n * EPS * s_elec,
        "{elec} against {}",
        i.electrostatic
    );
    eprintln!(
        "benzene in the 6 Å pocket, crystal pose: ΔE_bind {:.4} kcal/mol = van der Waals {:.4} + \
         electrostatic {:.4}; the difference of three evaluations agrees to {:.2e} kcal/mol, \
         {:.3} of its allowance at worst",
        kcal(i.total()),
        kcal(i.van_der_waals),
        kcal(i.electrostatic),
        kcal((d.total - i.total()).abs()),
        worst
    );
}

/// A bound on the polar desolvation of a binding whose partners are far apart, joules, given the
/// bound `coulomb_bound` on the vacuum cross Coulomb energy.
///
/// With `k′ = k (1 − 1/ε)`, `ΔG_GB = −½ k′ Σ_ij q_i q_j / f_ij` over ordered pairs, `i = j` included,
/// and the polar desolvation is three parts:
///
/// - **The cross pairs**, `−k′ Σ_{i∈pocket, j∈ligand} q_i q_j / f_ij`. Since `f ≥ r` and
///   `f − r = R_iR_j E / (f + r) ≤ R_iR_j E / 2r` (`E = e^(−r²/4R_iR_j)`), this is within
///   `k′ Σ |q_i q_j| R_iR_j E / 2r³` of `(1 − 1/ε)` times the vacuum cross Coulomb, which is within
///   `coulomb_bound`.
/// - **Each fragment's own terms, through its Born radii only.** By the mean value theorem on the
///   segment from the fragment's radii to the complex's, the change is at most
///   `Σ_i |ΔR_i| sup |∂G/∂R_i|`, and over the box of radii between the two sets
///   `|∂G/∂R_i| ≤ ½ k′ [q_i² / R_i,min² + 2 Σ_j |q_i q_j| R_j,max (1 + x) e^(−x) / 2r³]` with
///   `x = r² / 4 R_i,max R_j,max` — from `∂f/∂R_i = R_j E (1 + r²/4R_iR_j) / 2f`, `f ≥ r`, and
///   `(1 + x) e^(−x)` falling in `x`. The radii are the ones computed here, which are what the
///   energies use.
/// - **Rounding**: each of the three evaluations sums `n` terms, each within about ten ε of
///   itself (an `exp`, a `sqrt`, a few products and quotients), so within `(n + 10) ε Σ|t|`, with
///   `Σ|t| ≤ k′ Σ |q_i q_j| / r` (and `q_i² / R_i` on the diagonal); and two subtractions.
fn far_polar_bound(b: &Binding, coulomb_bound: f64) -> f64 {
    let (n0, n) = (b.pocket_len(), b.positions().len());
    let (at, q, el) = (b.positions(), b.charges(), b.elements());
    let radii = |r: std::ops::Range<usize>| GeneralizedBorn::new(&el[r.clone()]).born_radii(&at[r]);
    let bound = radii(0..n);
    let mut apart = radii(0..n0);
    apart.extend(radii(n0..n));
    let lo: Vec<f64> = bound.iter().zip(&apart).map(|(a, b)| a.min(*b)).collect();
    let hi: Vec<f64> = bound.iter().zip(&apart).map(|(a, b)| a.max(*b)).collect();
    let k = coulomb(1.0, 1.0, 1.0) * (1.0 - 1.0 / 80.0);
    let same = |i: usize, j: usize| (i < n0) == (j < n0);
    let (mut cross, mut fragments, mut sum_abs, mut terms) = (0.0, 0.0, 0.0, 0usize);
    for i in 0..n {
        let mut slope = q[i] * q[i] / (lo[i] * lo[i]);
        sum_abs += q[i] * q[i] / lo[i];
        terms += 1;
        for j in 0..n {
            if j == i {
                continue;
            }
            let r = distance(at[i], at[j]);
            let qq = (q[i] * q[j]).abs();
            if same(i, j) {
                let x = r * r / (4.0 * hi[i] * hi[j]);
                slope += 2.0 * qq * hi[j] * (1.0 + x) * (-x).exp() / (2.0 * r.powi(3));
            } else if i < j {
                let rr = hi[i] * hi[j];
                cross += qq * rr * (-r * r / (4.0 * rr)).exp() / (2.0 * r.powi(3));
            }
            sum_abs += qq / r;
            terms += 1;
        }
        fragments += (bound[i] - apart[i]).abs() * 0.5 * slope;
    }
    let rounding = 2.0 * (terms as f64 + 10.0) * EPS * 0.5 * k * sum_abs;
    (1.0 - 1.0 / 80.0) * coulomb_bound + k * (cross + fragments) + rounding
}

/// **Far away, ΔE_bind → 0**, inside a bound derived from the two tails. Benzene's QEq charges
/// sum to zero to rounding, so the Coulomb tail has no monopole: with `c` the ligand's centroid,
/// `a` its largest atom distance from `c`, `R` the nearest pocket atom's distance from `c`, `Q`
/// and `μ = Σ q_j (x_j − c)` the ligand's charge and dipole, the first-order Taylor expansion of
/// `1/|x − c − d|` in `d` with its remainder bounded by the Hessian of `1/r` (spectral norm
/// `2/r³`) gives the ligand's potential at a pocket atom as at most
/// `|Q|/R + |μ|/R² + a² Σ|q_j| / (R − a)³`, and the Coulomb energy is at most `k Σ_i |q_i|` times
/// that. Each van der Waals pair is at `r ≥ R − a`, past its minimum, where
/// `|D (s¹² − 2 s⁶)| ≤ 3 D s⁶`. The three-evaluation difference gets its rounding allowance on top.
///
/// **The Coulomb bound is loose, measured: about 470 times the value at 1000 Å.** It bounds the
/// pocket's side by `Σ_i |q_i|`, with no cancellation among the pocket's 310 charges, where the
/// true tail is the pocket's net +1 against benzene's small crystal dipole. A tighter bound would
/// expand the pocket's side too, and its remainder grows with the pocket's own extent (about
/// 12 Å here), which at these distances gives a looser bound, not a tighter one; none is claimed.
/// So this test sees a far field that does not decay, such as stale positions or an unmatched
/// term, and **does not see a Coulomb energy wrong by 1%**. That is caught only by
/// `the_vacuum_binding_energy_is_exactly_the_cross_terms`, when the error is in one place and not
/// the other.
#[test]
fn far_away_the_binding_energy_vanishes_inside_its_tails() {
    let b = crystal();
    let u = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0];
    for far in [1e3, 1e4] {
        let moved = b.ligand_at(&RigidMotion::translation([
            far * ANGSTROM * u[0],
            far * ANGSTROM * u[1],
            far * ANGSTROM * u[2],
        ]));
        let at = moved.positions();
        let q = moved.charges();
        let c = moved.ligand_centroid();
        let lig = moved.ligand_range();
        let a = lig.clone().map(|j| distance(at[j], c)).fold(0.0, f64::max);
        let r_min = (0..moved.pocket_len())
            .map(|i| distance(at[i], c))
            .fold(f64::INFINITY, f64::min);
        let total_q: f64 = lig.clone().map(|j| q[j]).sum();
        let mut mu = [0.0; 3];
        for j in lig.clone() {
            for k in 0..3 {
                mu[k] += q[j] * (at[j][k] - c[k]);
            }
        }
        let mu = len(mu);
        let abs_l: f64 = lig.clone().map(|j| q[j].abs()).sum();
        let abs_p: f64 = (0..moved.pocket_len()).map(|i| q[i].abs()).sum();
        let k = coulomb(1.0, 1.0, 1.0);
        let coulomb_bound = k
            * abs_p
            * (total_q.abs() / r_min + mu / (r_min * r_min) + a * a * abs_l / (r_min - a).powi(3));
        let vdw_bound: f64 = cross_pairs(&moved)
            .iter()
            .map(|t| 3.0 * t.2.well * (t.2.distance / (r_min - a)).powi(6))
            .sum();
        let i = moved.interaction();
        let allowance = difference_allowance(&moved)[6];
        assert!(
            i.van_der_waals.abs() <= vdw_bound,
            "{far} Å: van der Waals {:e} past {vdw_bound:e}",
            i.van_der_waals
        );
        assert!(
            i.electrostatic.abs() <= coulomb_bound,
            "{far} Å: Coulomb {:e} past {coulomb_bound:e}",
            i.electrostatic
        );
        assert!(
            i.difference.total.abs() <= coulomb_bound + vdw_bound + allowance,
            "{far} Å: difference {:e}",
            i.difference.total
        );
        // Solvation too: no grown sphere of one partner reaches the other's, so every atom's
        // surface points and neighbours are its fragment's and its area the same bits; and the
        // polar desolvation inside `far_polar_bound`.
        let d = moved.desolvation(POINTS);
        assert_eq!(d.areas_bound, d.areas_apart, "{far} Å");
        assert_eq!(d.buried_area, 0.0);
        let polar_bound = far_polar_bound(&moved, coulomb_bound);
        assert!(
            d.polar.abs() <= polar_bound,
            "{far} Å: polar desolvation {:e} past {polar_bound:e}",
            d.polar
        );
        eprintln!(
            "{far} Å away: polar desolvation {:.3e} kcal/mol (bound {:.3e})",
            kcal(d.polar),
            kcal(polar_bound)
        );
        eprintln!(
            "{far} Å away: van der Waals {:.3e} (bound {:.3e}), Coulomb {:.3e} (bound {:.3e}, \
             of which the charge {:.1e} e and dipole {:.3e} e Å), difference {:.3e} kcal/mol",
            kcal(i.van_der_waals),
            kcal(vdw_bound),
            kcal(i.electrostatic),
            kcal(coulomb_bound),
            total_q,
            mu / ANGSTROM,
            kcal(i.difference.total)
        );
    }
}

/// **A rigid motion is the one it is named for**, against closed forms. A quarter turn about z
/// centred at (1, 2, 3) takes (2, 2, 3) to (1, 3, 3); and about an oblique axis, every point goes
/// where Rodrigues' vector formula `v cos θ + (k × v) sin θ + k (k · v)(1 − cos θ)`, written here,
/// sends it — which a transposed matrix (the inverse rotation) or a rotation about the origin
/// fails. Each coordinate is a few operations on numbers no larger than `X`, so the two agree to
/// `16 ε X`. And `ligand_at` turns the ligand about its own centroid, which therefore stays put:
/// the centroid of the moved atoms is the old one plus `R` times the mean of `p − c`, which is
/// zero but for the rounding of the first centroid, `(n − 1) ε X` per coordinate for `n` atoms;
/// the moved atoms carry `10 ε X` each and the second centroid another `(n − 1) ε X`, so the
/// centroid moves by at most `(2n + 8) ε X`.
#[test]
fn a_rigid_motion_is_the_one_its_name_says() {
    let quarter = RigidMotion::rotation([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
    let got = quarter.apply([1.0, 2.0, 3.0], [2.0, 2.0, 3.0]);
    for (g, w) in got.iter().zip([1.0, 3.0, 3.0]) {
        assert!((g - w).abs() <= 16.0 * EPS * 3.0, "{got:?}");
    }
    let axis = [0.3, -0.5, 0.8];
    let n = len(axis);
    let k = [axis[0] / n, axis[1] / n, axis[2] / n];
    let theta = 2.1f64;
    let (sin, cos) = theta.sin_cos();
    let centre = [4.0, -1.5, 2.5];
    let motion = RigidMotion::rotation(axis, theta).then_translate([0.7, 0.2, -0.4]);
    for p in [[1.0, 2.0, 3.0], [-3.5, 0.25, 7.0], [10.0, -8.0, 0.5]] {
        let v = sub(p, centre);
        let kv = cross(k, v);
        let along = dot(k, v) * (1.0 - cos);
        let want: Vec<f64> = (0..3)
            .map(|i| centre[i] + v[i] * cos + kv[i] * sin + k[i] * along + motion.translation[i])
            .collect();
        let got = motion.apply(centre, p);
        let x = p
            .iter()
            .chain(&centre)
            .chain(&got)
            .fold(0.0f64, |m, a| m.max(a.abs()));
        for i in 0..3 {
            assert!(
                (got[i] - want[i]).abs() <= 16.0 * EPS * x,
                "{p:?}: {got:?} against {want:?}"
            );
        }
    }
    let b = crystal();
    let turned = b.ligand_at(&RigidMotion::rotation(axis, theta));
    let (c0, c1) = (b.ligand_centroid(), turned.ligand_centroid());
    let n = b.ligand_range().len() as f64;
    let x = b
        .ligand_positions()
        .iter()
        .chain(turned.ligand_positions())
        .flatten()
        .fold(0.0f64, |m, a| m.max(a.abs()));
    for i in 0..3 {
        assert!(
            (c1[i] - c0[i]).abs() <= (2.0 * n + 8.0) * EPS * x,
            "the centroid moved {:e} m on axis {i}",
            c1[i] - c0[i]
        );
    }
    // And it did turn: benzene's atoms moved.
    assert!(
        turned.ligand_rmsd() > 0.5 * ANGSTROM,
        "{}",
        turned.ligand_rmsd()
    );
}

/// **A change of frame changes nothing**: the whole complex turned 2.1 rad about an oblique axis
/// and moved 13.6 Å. Every cross distance moves by rounding only: each moved coordinate is a few
/// operations on numbers no larger than `X`, the largest coordinate or centroid in either frame,
/// so good to `10 ε X` per coordinate (one subtraction, three products and two sums, two more
/// sums, and the rotation's own entries a few ε each), and a distance between two such points
/// to `2√3 · 10 ε X` plus `4 ε r` for computing it in each frame. Each pair's energy then moves
/// by at most `|dE/dr|` times that, and the two cross sums round by `n ε Σ|t|` each.
#[test]
fn a_rigid_motion_of_the_whole_complex_leaves_the_binding_energy() {
    let b = crystal();
    let m = b.moved(
        &RigidMotion::rotation([0.3, -0.5, 0.8], 2.1).then_translate([
            7.0 * ANGSTROM,
            -3.0 * ANGSTROM,
            11.0 * ANGSTROM,
        ]),
    );
    let biggest = |b: &Binding| {
        b.positions()
            .iter()
            .flat_map(|p| p.iter().map(|x| x.abs()))
            .fold(0.0f64, f64::max)
    };
    let x = biggest(b).max(biggest(&m)).max(len(b.ligand_centroid()));
    let pairs = cross_pairs(b);
    let mut allowance = 0.0;
    let mut s = 0.0;
    for (_, _, p, e, r) in &pairs {
        let s6 = (p.distance / r).powi(6);
        let slope = (12.0 * p.well * (s6 - s6 * s6) / r).abs() + e.abs() / r;
        allowance += slope * (2.0 * 3f64.sqrt() * 10.0 * EPS * x + 8.0 * EPS * r);
        s += p.energy(*r).abs() + e.abs();
    }
    allowance += 2.0 * 2.0 * pairs.len() as f64 * EPS * s;
    let (before, after) = (b.interaction(), m.interaction());
    let err = (after.total() - before.total()).abs();
    assert!(
        err <= allowance,
        "ΔE_bind {:e} J moved to {:e} J, off {err:e} > {allowance:e}",
        before.total(),
        after.total()
    );
    // The difference of three evaluations too, with each frame's own allowance.
    let both = difference_allowance(b)[6] + difference_allowance(&m)[6];
    let err_d = (after.difference.total - before.difference.total).abs();
    assert!(
        err_d <= allowance + both,
        "{err_d:e} > {:e}",
        allowance + both
    );
    // Generalized Born: each of its three energies moves, to first order, by at most
    // `Σ_i |F_i| δp` with `δp = √3 · 10 ε X` each atom's position error, F the analytic GB force
    // (checked against central differences in `the_generalized_born_against_closed_forms.rs`),
    // and each evaluation in each frame rounds by `(n + 10) ε Σ|t|` (see `gb_rebuilt`).
    let (gb0, gb1) = (b.desolvation(1).polar, m.desolvation(1).polar);
    let n0 = b.pocket_len();
    let n = b.positions().len();
    let delta_p = 3f64.sqrt() * 10.0 * EPS * x;
    let mut gb_allowance = 0.0;
    for range in [0..n, 0..n0, n0..n] {
        let el = &b.elements()[range.clone()];
        let q = &b.charges()[range.clone()];
        for frame in [b, &m] {
            let at = &frame.positions()[range.clone()];
            let mut forces = vec![[0.0; 3]; at.len()];
            GeneralizedBorn::new(el).accumulate(q, at, &mut forces);
            let (_, sum_abs, terms) = gb_rebuilt(el, q, at);
            gb_allowance += forces.iter().map(|f| len(*f)).sum::<f64>() * delta_p
                + (terms as f64 + 10.0) * EPS * sum_abs;
        }
    }
    let gb_err = (gb1 - gb0).abs();
    assert!(
        gb_err <= gb_allowance,
        "polar desolvation moved {gb_err:e} > {gb_allowance:e}"
    );
    eprintln!(
        "frame change: ΔE_bind moved by {:.2e} kcal/mol ({:.2e} of its allowance); the polar \
         desolvation by {:.2e} kcal/mol ({:.2e} of its allowance)",
        kcal(err),
        err / allowance,
        kcal(gb_err),
        gb_err / gb_allowance
    );
}

/// ΔG_GB of atoms `el` with charges `q` at `at`, rebuilt here — eq 2 over every ordered pair,
/// `i = j` included, with Still's eq 3 typed in the test — from the crate's Born radii
/// ([`GeneralizedBorn::born_radii`], whose closed forms are held in
/// `the_generalized_born_against_closed_forms.rs`): `(energy, Σ|t|, terms)`, joules. Water, no
/// salt: `−½ k (1 − 1/80) q_i q_j / f_ij`.
fn gb_rebuilt(el: &[Element], q: &[f64], at: &[[f64; 3]]) -> (f64, f64, usize) {
    let r = GeneralizedBorn::new(el).born_radii(at);
    let k = coulomb(1.0, 1.0, 1.0) * (1.0 - 1.0 / 80.0);
    let (mut e, mut s, mut n) = (0.0, 0.0, 0usize);
    for i in 0..el.len() {
        for j in 0..el.len() {
            let f = if i == j {
                r[i]
            } else {
                let d = distance(at[i], at[j]);
                let rr = r[i] * r[j];
                (d * d + rr * (-d * d / (4.0 * rr)).exp()).sqrt()
            };
            let t = -0.5 * k * q[i] * q[j] / f;
            e += t;
            s += t.abs();
            n += 1;
        }
    }
    (e, s, n)
}

/// **The desolvation is rebuilt from its parts**, at the crystal pose and at the minimised one:
/// ΔG_GB of the complex, the pocket and the ligand against [`gb_rebuilt`] on each fragment's own
/// atoms at the binding's current positions, to `2 (n + 10) ε Σ|t|` — both sums of `n` terms each
/// within about ten ε of itself — and the polar term as their difference, with the three
/// allowances added; each atom's area against [`surface_area`] run here on the complex and on
/// each fragment, bit for bit; and ΔSASA as the sum over every atom of bound minus apart, bit for
/// bit, in the same order.
///
/// [`surface_area`]: pantometry_forcefield::solvation::surface_area
#[test]
fn the_desolvation_is_rebuilt_from_its_parts() {
    use pantometry_forcefield::solvation::{intrinsic_radius, surface_area, PROBE_RADIUS};
    for (name, b) in [("crystal", crystal()), ("minimised", &minimised().0)] {
        let d = b.desolvation(POINTS);
        let (n0, n) = (b.pocket_len(), b.positions().len());
        let mut allowances = 0.0;
        let mut rebuilt = [0.0; 3];
        for (k, (range, got)) in [(0..n, d.complex), (0..n0, d.pocket), (n0..n, d.ligand)]
            .into_iter()
            .enumerate()
        {
            let (e, sum_abs, terms) = gb_rebuilt(
                &b.elements()[range.clone()],
                &b.charges()[range.clone()],
                &b.positions()[range],
            );
            let allowance = 2.0 * (terms as f64 + 10.0) * EPS * sum_abs;
            assert!(
                (got - e).abs() <= allowance,
                "{name} fragment {k}: {got:e} against {e:e}, off {:e} > {allowance:e}",
                (got - e).abs()
            );
            eprintln!(
                "{name} fragment {k}: ΔG_GB {:.4} kcal/mol, rebuilt to {:.1e} kcal/mol, {:.2e} of \
                 its allowance",
                kcal(got),
                kcal((got - e).abs()),
                (got - e).abs() / allowance
            );
            allowances += allowance;
            rebuilt[k] = e;
        }
        let polar = rebuilt[0] - rebuilt[1] - rebuilt[2];
        assert!(
            (d.polar - polar).abs() <= allowances + 4.0 * EPS * rebuilt[0].abs(),
            "{name}: polar {:e} against {polar:e}",
            d.polar
        );
        let radii: Vec<f64> = b.elements().iter().map(|&e| intrinsic_radius(e)).collect();
        let area = |r: std::ops::Range<usize>| {
            surface_area(&radii[r.clone()], &b.positions()[r], PROBE_RADIUS, POINTS)
        };
        let bound = area(0..n);
        let mut apart = area(0..n0);
        apart.extend(area(n0..n));
        assert_eq!(d.areas_bound, bound, "{name}");
        assert_eq!(d.areas_apart, apart, "{name}");
        let buried: f64 = bound.iter().zip(&apart).map(|(b, a)| b - a).sum();
        assert_eq!(d.buried_area, buried, "{name}");
        eprintln!(
            "{name}: polar {:+.4} kcal/mol, ΔSASA {:.2} Å²",
            kcal(d.polar),
            d.buried_area / (ANGSTROM * ANGSTROM)
        );
    }
}

/// **A frozen atom does not move, to the bit**: after benzene is minimised in the pocket, every
/// protein atom is where it started, bit for bit — while benzene moved, and the minimiser
/// converged, so the test cannot pass by nothing moving.
#[test]
fn frozen_protein_atoms_do_not_move_to_the_bit() {
    let (b, p) = minimised();
    let c = crystal();
    let n0 = c.pocket_len();
    assert_eq!(p.status, Status::Converged, "{p:?}");
    for i in 0..n0 {
        assert_eq!(
            b.positions()[i].map(f64::to_bits),
            c.positions()[i].map(f64::to_bits),
            "protein atom {i} moved"
        );
    }
    assert!(b.ligand_rmsd() > 0.1 * ANGSTROM, "{}", b.ligand_rmsd());
}

/// **The same with the mask alone, on aspirin**: three atoms frozen while the rest relaxes; the
/// three stay bit for bit, every other atom's force goes below the tolerance, and the three
/// carry the force the frozen ones are held against — so the mask, not a stationary point, is
/// what kept them.
#[test]
fn a_frozen_mask_holds_its_atoms_on_a_small_molecule() {
    let c = pantometry_forcefield::Component::from_ccd(include_str!("../components/AIN.cif"))
        .expect("AIN parses");
    let ff = ForceField::new(&c, &pantometry_forcefield::uff::assign(&c)).expect("aspirin");
    let mut at: Vec<[f64; 3]> = c
        .atoms()
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let mut p = a.at;
            p[0] += 0.1 * ANGSTROM * (1.3 * i as f64).sin();
            p
        })
        .collect();
    // Each frozen coordinate moved by the first of a sequence of offsets, up to 0.5 Å, after which
    // `(x / Å) · Å` is not `x` — the round trip the minimiser's ångström vector makes — so a mask
    // that rebuilt frozen positions from that vector instead of copying them would move them, and
    // this test would see it. **Measured**: the round trip fails only where `x / Å` lies in a band
    // near the bottom of its binade (about [1, 1.2) · 2^k), for 6–11% of uniform values overall
    // and for none within 0.04 Å of 1.69 Å; so the offsets reach 0.5 Å, and each frozen atom is
    // required to have at least one coordinate that fails.
    let round_trips = |x: f64| (x / ANGSTROM) * ANGSTROM == x;
    for (i, p) in at.iter_mut().enumerate().take(3) {
        for (k, x) in p.iter_mut().enumerate() {
            let base = *x;
            for j in 0..1000 {
                let phase = 2.1 * (3 * i + k) as f64 + 0.4 + 0.77 * j as f64;
                let trial = base + 0.5 * ANGSTROM * phase.cos();
                if !round_trips(trial) {
                    *x = trial;
                    break;
                }
            }
        }
    }
    for p in &at[..3] {
        assert!(p.iter().any(|&x| !round_trips(x)), "{p:?}");
    }
    let start = at.clone();
    let frozen: Vec<bool> = (0..at.len()).map(|i| i < 3).collect();
    let mut m = Minimiser::new(TOLERANCE).with_frozen(frozen.clone());
    let mut p = m.step(&ff, &[], &mut at);
    while p.status == Status::Running && p.steps < 5000 {
        p = m.step(&ff, &[], &mut at);
    }
    assert_eq!(p.status, Status::Converged, "{p:?}");
    let forces = ff.evaluate(&at).forces;
    for i in 0..at.len() {
        if frozen[i] {
            assert_eq!(at[i].map(f64::to_bits), start[i].map(f64::to_bits));
        } else {
            assert!(
                len(forces[i]) <= TOLERANCE,
                "atom {i}: {:e}",
                len(forces[i])
            );
        }
    }
    let held = (0..3).map(|i| len(forces[i])).fold(0.0, f64::max);
    assert!(held > 100.0 * TOLERANCE, "{held:e}");
}

/// **Minimised in the pocket, the force on every benzene atom — the protein's on it included — is
/// below the tolerance**, evaluated by the complex's whole force field rather than the reduced
/// one the minimiser used. Reported: the pose and the interaction before and after.
#[test]
fn the_minimised_ligand_feels_no_force_from_anything() {
    let (b, p) = minimised();
    let c = crystal();
    let forces = b.force_field().evaluate(b.positions()).forces;
    let worst = b
        .ligand_range()
        .map(|l| len(forces[l]))
        .fold(0.0f64, f64::max);
    assert!(worst <= TOLERANCE, "{:e} N", worst);
    // What it minimised went down: the ligand's own energy plus the interaction.
    let before = c.ligand_energy().total + c.interaction().total();
    let after = b.ligand_energy().total + b.interaction().total();
    assert!(after < before, "{after} against {before}");
    let (i0, i1) = (c.interaction(), b.interaction());
    let mut gas = c.ligand_positions().to_vec();
    let g = c.ligand_force_field().minimise(&mut gas, 5000, TOLERANCE);
    eprintln!(
        "minimised in the pocket, {} steps: largest ligand force {:.2e} kcal/mol/Å; heavy-atom \
         RMSD from the crystal {:.3} Å, centroid moved {:.3} Å; ΔE_bind {:.3} → {:.3} kcal/mol \
         (van der Waals {:.3} → {:.3}, electrostatic {:.3} → {:.3}); benzene's own energy {:.3} → \
         {:.3}, against {:.3} at its vacuum minimum ({:?})",
        p.steps,
        worst / KCAL_PER_MOL_ANGSTROM,
        b.ligand_rmsd() / ANGSTROM,
        distance(b.ligand_centroid(), c.ligand_centroid()) / ANGSTROM,
        kcal(i0.total()),
        kcal(i1.total()),
        kcal(i0.van_der_waals),
        kcal(i1.van_der_waals),
        kcal(i0.electrostatic),
        kcal(i1.electrostatic),
        kcal(c.ligand_energy().total),
        kcal(b.ligand_energy().total),
        kcal(g.energy),
        g.status
    );
}

/// The derivatives of one cross pair's energy, `f = D (s¹² − 2 s⁶) + C / r`: `(f′, f″, f‴)`.
fn pair_derivatives(p: &Pair, coulomb_energy: f64, r: f64) -> (f64, f64, f64) {
    let s6 = (p.distance / r).powi(6);
    let s12 = s6 * s6;
    let d = p.well;
    let c = coulomb_energy * r;
    (
        12.0 * d * (s6 - s12) / r - c / (r * r),
        d * (156.0 * s12 - 84.0 * s6) / (r * r) + 2.0 * c / r.powi(3),
        d * (672.0 * s6 - 2184.0 * s12) / r.powi(3) - 6.0 * c / r.powi(4),
    )
}

/// **The interaction force on each benzene atom is −∇ of ΔE_bind with respect to it**, against a
/// central difference of the cross sum (which is ΔE_bind, by
/// [`the_vacuum_binding_energy_is_exactly_the_cross_terms`], and has none of the pocket's
/// rounding in it). Tolerance as in `forces_are_the_gradient.rs`: truncation `h²/6 |E‴|`, with
/// `|E‴|` along an axis bounded by `|f‴| + 3|f″|/r + 3|f′|/r²` summed over the atom's pairs, plus
/// rounding `δE / h`, with each pair's energy good to `70 ε D (s¹² + 2 s⁶) + 8 ε |C/r|` and the sum
/// to `2 n ε Σ|t|`. And the check must see the force: every atom's force is over a hundred
/// tolerances.
#[test]
fn the_interaction_force_is_minus_the_gradient_of_the_binding_energy() {
    let b = crystal();
    let pairs = cross_pairs(b);
    let mut delta_e = 0.0;
    let mut s = 0.0;
    for (_, _, p, e, r) in &pairs {
        let s6 = (p.distance / r).powi(6);
        delta_e += 70.0 * EPS * p.well * (s6 * s6 + 2.0 * s6) + 8.0 * EPS * e.abs();
        s += p.energy(*r).abs() + e.abs();
    }
    delta_e += 2.0 * pairs.len() as f64 * EPS * s;
    let delta_e = kcal(delta_e);
    let analytic = b.interaction_forces();
    let h = 1e-5;
    let (mut worst, mut least_visible) = (0.0f64, f64::INFINITY);
    for (k, l) in b.ligand_range().enumerate() {
        let m3: f64 = pairs
            .iter()
            .filter(|t| t.1 == l)
            .map(|(_, _, p, e, r)| {
                let (d1, d2, d3) = pair_derivatives(p, *e, *r);
                let ra = r / ANGSTROM;
                // Per Å: the derivatives are in J m⁻ⁿ.
                kcal(d3.abs() * ANGSTROM.powi(3))
                    + 3.0 * kcal(d2.abs() * ANGSTROM.powi(2)) / ra
                    + 3.0 * kcal(d1.abs() * ANGSTROM) / (ra * ra)
            })
            .sum();
        let tol = h * h / 6.0 * m3 + delta_e / h;
        let force = len(analytic[k]) * ANGSTROM / KCAL_PER_MOL;
        least_visible = least_visible.min(force / tol);
        assert!(
            force > 100.0 * tol,
            "atom {l}: force {force} against tol {tol}"
        );
        for axis in 0..3 {
            let mut plus = b.ligand_positions().to_vec();
            let mut minus = plus.clone();
            plus[k][axis] += h * ANGSTROM;
            minus[k][axis] -= h * ANGSTROM;
            let step = (plus[k][axis] - minus[k][axis]) / ANGSTROM;
            let e = |at: &[[f64; 3]]| kcal(b.with_ligand_positions(at).interaction().total());
            let numeric = -(e(&plus) - e(&minus)) / step;
            let exact = analytic[k][axis] * ANGSTROM / KCAL_PER_MOL;
            let err = (numeric - exact).abs();
            assert!(
                err <= tol,
                "atom {l} axis {axis}: analytic {exact} numeric {numeric}, off {err:e} > {tol:e}"
            );
            worst = worst.max(err / tol);
        }
    }
    eprintln!(
        "interaction force against central differences: worst {worst:.2e} of a tolerance; the \
         weakest atom's force is {least_visible:.0} tolerances"
    );
}

/// **What burying benzene does to areas and Born radii, exactly, and the polar desolvation it
/// costs.**
///
/// Exact, by construction: every atom's surface points are the same in the complex and in its
/// fragment, and the complex only adds spheres that can cover them, so **no atom's area grows on
/// binding** and ΔSASA < 0. Every descreening term is the integral of a non-negative integrand,
/// the complex only adds terms to each atom's sum, and OBC's `R(I)` is increasing — `dt/dΨ =
/// α − 2βΨ + 3γΨ²` has a negative discriminant, `4β² − 12αγ` = 2.56 − 58.2 — so **no Born radius
/// shrinks on binding**. That check calls [`GeneralizedBorn::born_radii`] on the fragments
/// directly: it holds the solvation module to the property on this system, and says nothing about
/// [`Binding::desolvation`], which `the_desolvation_is_rebuilt_from_its_parts` holds.
///
/// **The polar term is not small, and the test says where it comes from.** Generalized Born is
/// quadratic in the charges at fixed radii, so the polar desolvation splits exactly into the
/// pocket's own part (the ligand uncharged: the pocket desolvated by benzene's volume), the
/// ligand's own part (the pocket uncharged) and a cross part. Asserted: the pocket's and the
/// ligand's own parts are positive — each partner solvates less bound, which follows term by term
/// for each atom's Born term from the radii growing, and is not a theorem for the pair terms
/// between charges of opposite sign — and the ligand's part is no more than the solvation it has
/// alone, which is weak and is said to be. Reported: the total is +7.3 kcal/mol at 6 Å, not
/// small, and most of it is not benzene's.
#[test]
fn burying_benzene_buries_area_and_costs_polar_solvation() {
    let b = crystal();
    let d = b.desolvation(POINTS);
    for (k, (bound, apart)) in d.areas_bound.iter().zip(&d.areas_apart).enumerate() {
        assert!(bound <= apart, "atom {k}: {bound} > {apart}");
    }
    assert!(d.buried_area < 0.0);
    // Born radii, computed here.
    let (n0, n) = (b.pocket_len(), b.positions().len());
    let radii = |r: std::ops::Range<usize>| {
        GeneralizedBorn::new(&b.elements()[r.clone()]).born_radii(&b.positions()[r])
    };
    let bound = radii(0..n);
    let mut apart = radii(0..n0);
    apart.extend(radii(n0..n));
    for (k, (rb, ra)) in bound.iter().zip(&apart).enumerate() {
        assert!(
            rb >= ra,
            "atom {k}: Born radius {rb} bound against {ra} apart"
        );
    }
    let zeroed = |range: std::ops::Range<usize>| {
        let mut q = b.charges().to_vec();
        q[range].fill(0.0);
        b.desolvation_with(&q, 1).polar
    };
    let pocket_part = zeroed(n0..n);
    let ligand_part = zeroed(0..n0);
    let cross_part = d.polar - pocket_part - ligand_part;
    assert!(pocket_part > 0.0, "{pocket_part}");
    assert!(ligand_part > 0.0, "{ligand_part}");
    // Below the whole solvation it has alone: this says only that the complex's GB energy of
    // benzene's charges alone is negative, which a positive-definite GB form makes it.
    assert!(
        ligand_part <= -d.ligand,
        "ligand's own desolvation {ligand_part} against its solvation {}",
        d.ligand
    );
    let ligand_area: f64 = d.areas_apart[n0..].iter().sum();
    eprintln!(
        "desolvation: polar {:+.3} kcal/mol = pocket by benzene's volume {:+.3} + benzene's own \
         {:+.3} (of the {:.3} it has alone) + cross {:+.3}; ΔSASA {:.1} Å² of benzene's {:.1} \
         alone, nonpolar {:+.3}; largest Born radius growth {:.3} Å",
        kcal(d.polar),
        kcal(pocket_part),
        kcal(ligand_part),
        kcal(d.ligand),
        kcal(cross_part),
        d.buried_area / (ANGSTROM * ANGSTROM),
        ligand_area / (ANGSTROM * ANGSTROM),
        kcal(d.nonpolar),
        bound
            .iter()
            .zip(&apart)
            .map(|(a, b)| a - b)
            .fold(0.0f64, f64::max)
            / ANGSTROM
    );
}

/// One row of the cutoff table.
fn row(b: &Binding) -> String {
    let i = b.interaction();
    let d = b.desolvation(POINTS);
    let mut m = b.clone();
    let p = m.minimise_ligand(5000, TOLERANCE);
    let j = m.interaction();
    let q = b.polarised_charges().expect("QEq on the complex");
    let transfer: f64 = q.charges[b.ligand_range()].iter().sum();
    let moved = b
        .ligand_range()
        .map(|k| (q.charges[k] - b.charges()[k]).abs())
        .fold(0.0f64, f64::max);
    let polarised = b.cross_electrostatic(&q.charges);
    // The identity at this cutoff too.
    let allowance = difference_allowance(b)[6];
    assert!((i.difference.total - i.total()).abs() <= allowance);
    assert_eq!(p.status, Status::Converged);
    format!(
        "| {:.0} Å | {} | {} | {:+} | {:.3} | {:.3} | {:.3} | {:+.3} | {:.1} | {:+.3} | {:.3} | {:.3} | \
         {:.3} | {:+.4} | {:.4} | {:.3} |",
        b.cutoff() / ANGSTROM,
        b.residues().len(),
        b.positions().len(),
        b.formal_charges().0,
        kcal(i.van_der_waals),
        kcal(i.electrostatic),
        kcal(i.total()),
        kcal(d.polar),
        d.buried_area / (ANGSTROM * ANGSTROM),
        kcal(d.nonpolar),
        kcal(i.total() + d.polar + d.nonpolar),
        m.ligand_rmsd() / ANGSTROM,
        kcal(j.total()),
        transfer,
        moved,
        kcal(polarised)
    )
}

const HEADER: &str = "| cutoff | residues | atoms | pocket charge | vdW | elec | ΔE vacuum | \
     polar | ΔSASA Å² | nonpolar | ΔE + ΔΔG_solv | minimised RMSD Å | ΔE minimised | \
     QEq on complex: charge to benzene e | its largest benzene charge change e | its cross \
     elec |";

/// **The cutoff's effect, measured**: 6, 8 and 10 Å, each with the counts asserted and the
/// vacuum identity checked again. Energies in kcal/mol. Ignored: QEq on the 10 Å pocket and on
/// its complex is 987 atoms twice, about a minute with `--release`.
#[test]
#[ignore = "QEq on pockets of up to 987 atoms; about a minute with --release -- --ignored"]
fn the_cutoffs_effect() {
    let s = system();
    eprintln!("{HEADER}");
    for (cutoff, residues, atoms, charge) in
        [(6.0, 18, 322, 1), (8.0, 43, 719, 3), (10.0, 60, 987, 6)]
    {
        let b = Binding::new(&s, cutoff * ANGSTROM).expect("pocket");
        assert_eq!(
            (
                b.residues().len(),
                b.positions().len(),
                b.formal_charges().0
            ),
            (residues, atoms, charge)
        );
        eprintln!("{}", row(&b));
    }
}

/// **The whole protein as the pocket**: every residue (a cutoff past the protein's size), the
/// reference the table converges towards — or does not. Ignored: QEq on 2616 atoms, twice, is
/// minutes even with `--release`.
#[test]
#[ignore = "QEq on all 2616 atoms, twice; minutes with --release -- --ignored"]
fn the_whole_protein_as_the_pocket() {
    let s = system();
    let b = Binding::new(&s, 1000.0 * ANGSTROM).expect("whole protein");
    assert_eq!(b.residues().len(), 162);
    eprintln!("{HEADER}");
    eprintln!("{}", row(&b));
}

/// **How much of the pose's move is the placed hydrogens'**: the same minimisation with every
/// protein hydrogen free as well, only the protein's heavy atoms frozen. Reported. Ignored: it
/// evaluates the whole complex's force field thousands of times.
#[test]
#[ignore = "thousands of whole-complex evaluations; seconds with --release -- --ignored"]
fn the_placed_hydrogens_share_of_the_move() {
    let b = crystal();
    let el = b.elements();
    let n0 = b.pocket_len();
    let frozen: Vec<bool> = (0..el.len())
        .map(|i| i < n0 && el[i] != Element::H)
        .collect();
    let mut at = b.positions().to_vec();
    let mut m = Minimiser::new(TOLERANCE).with_frozen(frozen);
    let ff = b.force_field();
    let mut p = m.step(ff, &[], &mut at);
    while p.status == Status::Running && p.steps < 20000 {
        p = m.step(ff, &[], &mut at);
    }
    let relaxed = b.with_ligand_positions(&at[n0..]);
    let mut contacts: Vec<(f64, usize, usize)> = Vec::new();
    for l in b.ligand_range() {
        for i in 0..n0 {
            contacts.push((distance(b.positions()[l], b.positions()[i]), l, i));
        }
    }
    contacts.sort_by(|a, b| a.0.total_cmp(&b.0));
    let (d, l, i) = contacts[0];
    let s = system();
    eprintln!(
        "closest crystal contact: benzene {:?} to {} at {:.3} Å; with the protein's hydrogens \
         free too ({:?}, {} steps) benzene moves {:.3} Å heavy-atom RMSD, against {:.3} Å with \
         them frozen",
        el[l],
        s.component().atoms()[b.system_atoms()[i]].name,
        d / ANGSTROM,
        p.status,
        p.steps,
        relaxed.ligand_rmsd() / ANGSTROM,
        minimised().0.ligand_rmsd() / ANGSTROM
    );
}
