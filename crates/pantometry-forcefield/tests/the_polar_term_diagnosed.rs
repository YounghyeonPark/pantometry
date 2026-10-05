//! **Why OBC II's polar desolvation does not converge with the pocket's cutoff, measured.** Step
//! A-2.
//!
//! Benzene (181L) and n-butylbenzene (186L), the series' largest ligand, at the crystal pose. GB
//! is a quadratic form in the charges at fixed Born radii, and a zero charge changes no radius, so
//! the polar term splits exactly into the pocket desolvated by the ligand's volume (the ligand's
//! charges zero), the ligand's own (the pocket's zero) and the screened cross terms. Each part is
//! then attributed to the pocket's atoms — an atom's own Born term and half of each pair term it
//! is in — and binned by distance from the ligand, which says where the term comes from and
//! whether a protein atom far away is desolvated by it. Then three questions:
//!
//! - **does the descreening integral decay as it should** with the ligand far from an atom: each
//!   atom past 10 Å against the asymptote `Σ s³/3r⁴` of a sphere of radius `s` at distance `r`
//!   (eq 5 with `r⁻⁴` constant over the sphere);
//! - **does OBC II's rescaling amplify a buried atom's response**: `d(1/R)/dI` is
//!   `(ρ̃/ρ) sech² t · t′(Ψ)` against HCT's 1, and its linearisation against the actual change;
//! - **is the growth the cut or the charges**: the whole protein's QEq charges on each cut
//!   pocket, and a cut pocket's charges inside the whole protein's dielectric (every other
//!   protein atom present and uncharged).
//!
//! Asserted: only that the attribution sums to the part it attributes, within the rounding of
//! its sums. **That allowance is the a-priori worst case and it is loose**: 3.2e-5 kcal/mol at
//! 6 Å, where the shares agree to 4.4e-13, and 9.5e-2 for the whole protein, where they agree to
//! 4.0e-11 — so it is the cut pockets' rows that would see an attribution that drops a term.
//! Everything else is printed and recorded in the `binding` module's documentation and the
//! changelog. Ignored: QEq on all 2616 protein atoms for each ligand, 751 s for both tests with
//! `cargo test -p pantometry-forcefield --release --test the_polar_term_diagnosed -- --ignored`.

mod protein;

use pantometry_forcefield::energy::coulomb;
use pantometry_forcefield::solvation::{
    descreening, intrinsic_radius, scale_factor, still_distance, GeneralizedBorn, RADIUS_OFFSET,
};
use pantometry_forcefield::uff::KCAL_PER_MOL;
use pantometry_forcefield::{Binding, Component, Element, Selection, System};
use protein::*;
use std::collections::HashSet;

const ANGSTROM: f64 = 1e-10;

fn kcal(joules: f64) -> f64 {
    joules / KCAL_PER_MOL
}

/// 186L, n-butylbenzene, built as `the_congener_series.rs` builds it.
fn butylbenzene() -> System {
    let mut t = templates();
    t.push(Component::from_ccd(include_str!("../components/N4B.cif")).expect("N4B"));
    System::from_pdb(
        include_str!("../components/186L.pdb"),
        &t,
        &Selection::new("N4B").dropping(&["HOH", "CL", "HED"]),
    )
    .expect("186L")
}

/// The polar term and its three parts, kcal/mol.
#[derive(Clone, Copy)]
struct Parts {
    polar: f64,
    pocket: f64,
    ligand: f64,
    cross: f64,
}

/// `b`'s polar term with `charges`, split by zeroing each partner's charges in turn.
fn split(b: &Binding, charges: &[f64]) -> Parts {
    let (n0, n) = (b.pocket_len(), b.positions().len());
    let polar = |q: &[f64]| kcal(b.desolvation_with(q, 1).polar);
    let zeroed = |r: std::ops::Range<usize>| {
        let mut q = charges.to_vec();
        q[r].fill(0.0);
        q
    };
    let (all, pocket, ligand) = (polar(charges), polar(&zeroed(n0..n)), polar(&zeroed(0..n0)));
    Parts {
        polar: all,
        pocket,
        ligand,
        cross: all - pocket - ligand,
    }
}

/// One pocket atom's share: its distance from the nearest ligand atom (Å), its Born term's and its
/// half of the pair terms' change between the pocket alone and the complex (the pocket part), its
/// half of the cross terms, the ligand's descreening of it `ΔI` (Å⁻¹), OBC II's amplification in
/// the apo pocket, `q²`, its Born radius apart and bound (Å), and the `Σ|t|` and term count of what
/// went into it (kcal/mol).
struct Share {
    distance: f64,
    born: f64,
    pairs: f64,
    cross: f64,
    delta_i: f64,
    amplification: f64,
    q2: f64,
    apart: f64,
    bound: f64,
    magnitude: f64,
    terms: usize,
}

/// `(ρ̃/ρ) sech² t · t′(Ψ)` for OBC II, eq 8: `d(1/R)/dI` against HCT's 1.
fn amplification(e: Element, i: f64) -> f64 {
    let rho = intrinsic_radius(e);
    let rho_tilde = rho - RADIUS_OFFSET;
    let psi = i * rho_tilde;
    let t = psi - 0.8 * psi * psi + 4.85 * psi.powi(3);
    let dt = 1.0 - 1.6 * psi + 3.0 * 4.85 * psi * psi;
    let th = t.tanh();
    rho_tilde / rho * (1.0 - th * th) * dt
}

/// Eq 5's integral for atom `i` over the atoms `from`, metres⁻¹.
fn integral(el: &[Element], at: &[[f64; 3]], i: usize, from: std::ops::Range<usize>) -> f64 {
    let rho_tilde = intrinsic_radius(el[i]) - RADIUS_OFFSET;
    from.filter(|&j| j != i)
        .map(|j| {
            let s = scale_factor(el[j]) * (intrinsic_radius(el[j]) - RADIUS_OFFSET);
            descreening(len(sub(at[i], at[j])), rho_tilde, s).0
        })
        .sum()
}

/// Every pocket atom's [`Share`] with `q` (complex order), at the crystal pose.
fn shares(b: &Binding, q: &[f64]) -> Vec<Share> {
    let (n0, n) = (b.pocket_len(), b.positions().len());
    let (at, el) = (b.positions(), b.elements());
    let bound = GeneralizedBorn::new(el).born_radii(at);
    let apart = GeneralizedBorn::new(&el[..n0]).born_radii(&at[..n0]);
    let k = coulomb(1.0, 1.0, 1.0) * (1.0 - 1.0 / 80.0);
    (0..n0)
        .map(|i| {
            let mut magnitude = 0.0;
            let mut add = |t: f64| {
                magnitude += t.abs();
                t
            };
            let born =
                add(-0.5 * k * q[i] * q[i] / bound[i]) - add(-0.5 * k * q[i] * q[i] / apart[i]);
            let mut pairs = 0.0;
            for j in (0..n0).filter(|&j| j != i) {
                let r = len(sub(at[i], at[j]));
                let c = -0.5 * k * q[i] * q[j];
                pairs += add(c / still_distance(r, bound[i], bound[j]))
                    - add(c / still_distance(r, apart[i], apart[j]));
            }
            let mut cross = 0.0;
            for l in n0..n {
                let r = len(sub(at[i], at[l]));
                cross += add(-k * q[i] * q[l] / still_distance(r, bound[i], bound[l]));
            }
            Share {
                distance: (n0..n)
                    .map(|l| len(sub(at[i], at[l])))
                    .fold(f64::INFINITY, f64::min)
                    / ANGSTROM,
                born: kcal(born),
                pairs: kcal(pairs),
                cross: kcal(cross),
                delta_i: integral(el, at, i, n0..n) * ANGSTROM,
                amplification: amplification(el[i], integral(el, at, i, 0..n0)),
                q2: q[i] * q[i],
                apart: apart[i] / ANGSTROM,
                bound: bound[i] / ANGSTROM,
                magnitude: kcal(magnitude),
                terms: 2 * n0 + n - n0,
            }
        })
        .collect()
}

/// `2 (n² + 10) ε Σ|t|` for one ΔG_GB evaluation of atoms `el` with charges `q` at `at`: its `n²`
/// terms `½ k (1 − 1/80) |q_i q_j| / f_ij`, each within about ten ε of itself, summed — the
/// allowance `benzene_in_its_pocket.rs` derives for the same sum. Kcal/mol.
fn evaluation_rounding(el: &[Element], q: &[f64], at: &[[f64; 3]]) -> f64 {
    let r = GeneralizedBorn::new(el).born_radii(at);
    let k = coulomb(1.0, 1.0, 1.0) * (1.0 - 1.0 / 80.0);
    let n = el.len();
    let mut sum = 0.0;
    for i in 0..n {
        sum += 0.5 * k * q[i] * q[i] / r[i];
        for j in (0..n).filter(|&j| j != i) {
            let d = len(sub(at[i], at[j]));
            sum += 0.5 * k * (q[i] * q[j]).abs() / still_distance(d, r[i], r[j]);
        }
    }
    2.0 * ((n * n) as f64 + 10.0) * f64::EPSILON * kcal(sum)
}

/// Asserts that the shares sum to the pocket part and the cross terms of `parts`, within the
/// rounding of the shares' own sums (`n ε Σ|t|`) and of the evaluations each part is a difference
/// of ([`evaluation_rounding`]: the complex three times, the pocket alone and the ligand alone
/// once each, with the binding's charges bounding every partial charge set's `Σ|t|`), and prints
/// the shells.
fn report(label: &str, b: &Binding, q: &[f64]) -> Vec<Share> {
    let parts = split(b, q);
    let rows = shares(b, q);
    let (n0, at, el) = (b.pocket_len(), b.positions(), b.elements());
    let magnitude: f64 = rows.iter().map(|r| r.magnitude).sum();
    let terms: usize = rows.iter().map(|r| r.terms).sum();
    let allowance = terms as f64 * f64::EPSILON * magnitude
        + 3.0 * evaluation_rounding(el, q, at)
        + evaluation_rounding(&el[..n0], &q[..n0], &at[..n0])
        + evaluation_rounding(&el[n0..], &q[n0..], &at[n0..]);
    let pocket: f64 = rows.iter().map(|r| r.born + r.pairs).sum();
    let cross: f64 = rows.iter().map(|r| r.cross).sum();
    assert!(
        (pocket - parts.pocket).abs() <= allowance,
        "{label}: shares {pocket} against the pocket part {}",
        parts.pocket
    );
    assert!(
        (cross - parts.cross).abs() <= allowance,
        "{label}: shares {cross} against the cross terms {}",
        parts.cross
    );
    eprintln!(
        "\n{label}: polar {:+.3} = pocket by volume {:+.3} + ligand's own {:+.3} + cross {:+.3} \
         kcal/mol; the shares sum to the pocket part within {:.1e} and the cross terms within \
         {:.1e}, against {allowance:.1e}",
        parts.polar,
        parts.pocket,
        parts.ligand,
        parts.cross,
        (pocket - parts.pocket).abs(),
        (cross - parts.cross).abs()
    );
    eprintln!(
        "| shell Å | atoms | Σq² e² | ΣΔI Å⁻¹ | Born terms | pair terms | pocket part | cross | \
         mean OBC amplification |"
    );
    for w in [0.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, f64::INFINITY].windows(2) {
        let shell: Vec<&Share> = rows
            .iter()
            .filter(|r| r.distance >= w[0] && r.distance < w[1])
            .collect();
        let sum = |f: &dyn Fn(&Share) -> f64| shell.iter().map(|r| f(r)).sum::<f64>();
        eprintln!(
            "| {}–{} | {} | {:.2} | {:.4} | {:+.3} | {:+.3} | {:+.3} | {:+.3} | {:.2} |",
            w[0],
            w[1],
            shell.len(),
            sum(&|r| r.q2),
            sum(&|r| r.delta_i),
            sum(&|r| r.born),
            sum(&|r| r.pairs),
            sum(&|r| r.born + r.pairs),
            sum(&|r| r.cross),
            sum(&|r| r.amplification) / shell.len().max(1) as f64
        );
    }
    rows
}

/// The far field, the amplification, and where in the whole protein the term comes from.
fn whole_protein(name: &str, b: &Binding, rows: &[Share]) {
    let (n0, n) = (b.pocket_len(), b.positions().len());
    let (at, el) = (b.positions(), b.elements());
    let mut worst: f64 = 0.0;
    let mut count = 0;
    for (i, r) in rows.iter().enumerate() {
        if r.distance < 10.0 {
            continue;
        }
        let asymptote: f64 = (n0..n)
            .map(|l| {
                let s = scale_factor(el[l]) * (intrinsic_radius(el[l]) - RADIUS_OFFSET);
                s.powi(3) / (3.0 * len(sub(at[i], at[l])).powi(4))
            })
            .sum();
        worst = worst.max((integral(el, at, i, n0..n) / asymptote - 1.0).abs());
        count += 1;
    }
    eprintln!(
        "{name}: {count} atoms past 10 Å, the ligand's descreening of each within {:.2}% of Σ s³/3r⁴",
        100.0 * worst
    );
    let mut a: Vec<f64> = rows.iter().map(|r| r.amplification).collect();
    a.sort_by(f64::total_cmp);
    let m = a.len();
    eprintln!(
        "{name}: OBC II's d(1/R)/dI against HCT's in the apo protein: {:.3} to {:.3}, median \
         {:.3}, deciles {:.3} and {:.3}",
        a[0],
        a[m - 1],
        a[m / 2],
        a[m / 10],
        a[9 * m / 10]
    );
    let far: Vec<&Share> = rows.iter().filter(|r| r.distance > 10.0).collect();
    let actual: f64 = far
        .iter()
        .map(|r| r.q2 * (1.0 / r.bound - 1.0 / r.apart))
        .sum();
    let linear: f64 = far
        .iter()
        .map(|r| -r.q2 * r.amplification * r.delta_i)
        .sum();
    let hct: f64 = far.iter().map(|r| -r.q2 * r.delta_i).sum();
    eprintln!(
        "{name}: past 10 Å, Σ q² Δ(1/R) is {actual:.4e} e² Å⁻¹; OBC II linearised {linear:.4e}, \
         HCT {hct:.4e}"
    );
    for within in [6.0, 8.0, 10.0, 15.0, 20.0, f64::INFINITY] {
        let inner = rows.iter().filter(|r| r.distance < within);
        let (pocket, cross, born) = inner.fold((0.0, 0.0, 0.0), |(p, c, s), r| {
            (p + r.born + r.pairs, c + r.cross, s + r.born)
        });
        eprintln!(
            "{name}: protein atoms within {within} Å of the ligand carry pocket part {pocket:+.3} \
             (Born terms {born:+.3}) and cross {cross:+.3}"
        );
    }
}

/// The whole protein's QEq charges on `b`'s atoms, by system index.
fn charges_from(whole: &Binding, b: &Binding) -> Vec<f64> {
    let by_system: std::collections::HashMap<usize, f64> = whole
        .system_atoms()
        .iter()
        .copied()
        .zip(whole.charges().iter().copied())
        .collect();
    b.system_atoms().iter().map(|i| by_system[i]).collect()
}

/// The whole protein's charges with every atom not in `b` uncharged, complex order of `whole`.
fn only_inside(whole: &Binding, b: &Binding) -> Vec<f64> {
    let inside: HashSet<usize> = b.system_atoms().iter().copied().collect();
    whole
        .system_atoms()
        .iter()
        .zip(whole.charges())
        .map(|(i, &q)| if inside.contains(i) { q } else { 0.0 })
        .collect()
}

fn diagnose(name: &str, s: &System) {
    let whole = Binding::new(s, 1000.0 * ANGSTROM).expect("the whole protein");
    let rows = report(
        &format!("{name}, the whole protein"),
        &whole,
        whole.charges(),
    );
    whole_protein(name, &whole, &rows);
    for cutoff in [6.0, 8.0, 10.0] {
        let b = Binding::new(s, cutoff * ANGSTROM).expect("a pocket");
        report(&format!("{name} {cutoff} Å, its own QEq"), &b, b.charges());
        let q = charges_from(&whole, &b);
        let p = split(&b, &q);
        let moved = q
            .iter()
            .zip(b.charges())
            .map(|(a, c)| (a - c).abs())
            .fold(0.0f64, f64::max);
        eprintln!(
            "{name} {cutoff} Å with the whole protein's QEq charges: polar {:+.3} = {:+.3} + \
             {:+.3} + {:+.3}; the largest charge differs by {moved:.3} e",
            p.polar, p.pocket, p.ligand, p.cross
        );
    }
}

/// **Benzene and n-butylbenzene, split at 6, 8 and 10 Å and with the whole protein**, the whole
/// protein's shells, far field and amplification, and each cut pocket with the whole protein's
/// charges. Ignored: QEq on 2616 atoms for each.
#[test]
#[ignore = "QEq on all 2616 protein atoms for each of two ligands; minutes with --release -- --ignored"]
fn the_polar_term_split_at_each_cutoff() {
    diagnose("181L", &system());
    diagnose("186L", &butylbenzene());
}

/// **A cut pocket's charges in the whole protein's dielectric**: benzene's pocket at 6 to 20 Å,
/// its own fragment against the whole protein with every atom outside the pocket present and
/// uncharged. Separates what the cut does to the Born radii from the charges it leaves out.
/// Ignored: QEq on up to 2616 atoms, six times.
#[test]
#[ignore = "QEq on pockets of up to 2616 atoms, six times; minutes with --release -- --ignored"]
fn a_cut_pockets_charges_in_the_whole_dielectric() {
    let s = system();
    let whole = Binding::new(&s, 1000.0 * ANGSTROM).expect("the whole protein");
    eprintln!(
        "| cutoff | atoms | own fragment: polar | pocket | ligand's own | cross | in the whole \
         dielectric: polar | pocket | ligand's own | cross |"
    );
    for cutoff in [6.0, 8.0, 10.0, 15.0, 20.0] {
        let b = Binding::new(&s, cutoff * ANGSTROM).expect("a pocket");
        let own = split(&b, b.charges());
        let inside = split(&whole, &only_inside(&whole, &b));
        eprintln!(
            "| {cutoff} Å | {} | {:+.3} | {:+.3} | {:+.3} | {:+.3} | {:+.3} | {:+.3} | {:+.3} | \
             {:+.3} |",
            b.positions().len(),
            own.polar,
            own.pocket,
            own.ligand,
            own.cross,
            inside.polar,
            inside.pocket,
            inside.ligand,
            inside.cross
        );
    }
    let all = split(&whole, whole.charges());
    eprintln!(
        "| whole | {} | {:+.3} | {:+.3} | {:+.3} | {:+.3} | | | | |",
        whole.positions().len(),
        all.polar,
        all.pocket,
        all.ligand,
        all.cross
    );
}
