//! **Decoupling under OBC II generalized Born, and benzene's absolute binding free energy to T4
//! lysozyme L99A by double decoupling** — step 3c-3.
//!
//! The default tests hold what must be exact or analytic in the coupling of a decoupled group to
//! generalized Born ([`pantometry_forcefield::alchemy`], "Generalized Born"), on 181L's pocket at
//! 3 Å (benzene and four residues, 82 atoms) with QEq charges, and on one and two ions whose GB
//! energy is written out here:
//!
//! - fully coupled, the solvation term is the force field's to the bit, and the whole energy is
//!   the solvated force field's to a traced allowance;
//! - with the charges off, the group's charges are in no term, to the bit;
//! - fully decoupled, the energy is the pocket alone in solvent plus benzene alone in vacuum — the
//!   solvation term and its forces to the bit, the total to its rounding;
//! - `∂U/∂λ` and the forces are the energy's derivatives;
//! - the frozen buffer's kept descreening is used, changes no bit, and is dropped when a frozen
//!   atom moves; and the states of one configuration evaluated together are each evaluated alone;
//! - a Born ion decoupled through [`Windows`] gives its Born energy `−(q²/2R)(1 − 1/ε)` by TI and
//!   by BAR, and two ions follow the closed form of eqs 2–8 at every `(λ_e, λ_v)`.
//!
//! The measurements are ignored and run once in release; they print, and assert only what must
//! hold of any run — the buffer's bits, the kept descreening, nothing to decouple in the solvent
//! leg's van der Waals. Each window is written to a file in `PANTOMETRY_3C3_DIR` (default: the
//! system's temporary directory) as it runs, and a window already there is read instead, so a
//! stopped run resumes; a second run in the same directory refuses. The solvent leg first, since
//! the complex leg closes the cycle with its result:
//!
//! `cargo test -p pantometry-forcefield --release --test benzene_bound_in_generalized_born the_solvent_leg_measured -- --ignored --nocapture`
//!
//! `cargo test -p pantometry-forcefield --release --test benzene_bound_in_generalized_born the_complex_leg_measured -- --ignored --nocapture`
//!
//! **What they measured** (one core of the machine this was written on): the solvent leg
//! +2.3395 ± 0.0003 kcal/mol in 178 s; the complex leg +4.047 ± 0.309 in 5.11 h; with the
//! restraint's release, −7.784, **ΔG°_bind = +6.08 ± 0.31 kcal/mol against experiment's −5.19 ±
//! 0.16** — +4.61 with a nonpolar estimate at the end states. See `pantometry_forcefield::alchemy`
//! for the table and what limits it.

mod protein;

use pantometry_forcefield::energy::COULOMB_KCAL;
use pantometry_forcefield::solvation::{descreening, RADIUS_OFFSET, WATER_DIELECTRIC};
use pantometry_forcefield::uff::KCAL_PER_MOL;
use pantometry_forcefield::{
    Alchemical, Binding, Complex, Decoupling, Element, GeneralizedBorn, Lambda, Protocol,
    Quadrature, Solvent, System, Windows,
};
use pantometry_units::BOLTZMANN;
use protein::*;
use std::sync::OnceLock;

const ANGSTROM: f64 = 1e-10;
const FS: f64 = 1e-15;
const KELVIN: f64 = 300.0;
const EPS: f64 = f64::EPSILON;

fn kcal(j: f64) -> f64 {
    j / KCAL_PER_MOL
}

/// The Coulomb constant of eq 2 in SI, J m e⁻²: 332.0637 kcal mol⁻¹ Å e⁻².
fn coulomb_si() -> f64 {
    COULOMB_KCAL * KCAL_PER_MOL * ANGSTROM
}

fn the_system() -> &'static System {
    static S: OnceLock<System> = OnceLock::new();
    S.get_or_init(system)
}

/// 181L's pocket at 3 Å, once.
fn small_binding() -> &'static Binding {
    static B: OnceLock<Binding> = OnceLock::new();
    B.get_or_init(|| Binding::new(the_system(), 3.0 * ANGSTROM).expect("181L at 3 Å"))
}

fn ligand_mask(b: &Binding) -> Vec<bool> {
    (0..b.positions().len())
        .map(|k| k >= b.pocket_len())
        .collect()
}

fn same_forces(a: &[[f64; 3]], b: &[[f64; 3]], what: &str) {
    for (k, (f, g)) in a.iter().zip(b).enumerate() {
        for c in 0..3 {
            // `==`, not the bits: a term that is `±0` may leave a zero of either sign.
            assert!(
                f[c] == g[c],
                "{what}: atom {k} {c}: {:e} against {:e}",
                f[c],
                g[c]
            );
        }
    }
}

/// **Fully coupled, the solvation term is the force field's own, to the bit**:
/// `GeneralizedBorn::decoupled` at `(1, 1)` is `accumulate`, energy and every force, on 181L's
/// 3 Å pocket and benzene with QEq charges. **And the whole energy is the solvated force field's**
/// to a traced allowance, the same terms summed in another order: `n ε Σ|t|` with `n` the number
/// of terms (each GB pair counted once more for its two orderings) and `Σ|t|` bounded by every
/// field's magnitude plus the cross terms' twice.
#[test]
fn coupled_the_solvation_term_is_the_force_fields_to_the_bit() {
    let b = small_binding();
    let ff = b.force_field().clone().with_generalized_born();
    let gb = ff.solvation().unwrap();
    let at = b.positions();
    let n = at.len();
    let mask = ligand_mask(b);
    let mut whole = vec![[0.0; 3]; n];
    let e = gb.accumulate(b.charges(), at, &mut whole);
    let mut part = vec![[0.0; 3]; n];
    let d = gb.decoupled(b.charges(), &mask, at, 1.0, 1.0, &mut part);
    assert_eq!(d.energy.to_bits(), e.to_bits());
    for (f, g) in part.iter().zip(&whole) {
        for c in 0..3 {
            assert_eq!(f[c].to_bits(), g[c].to_bits());
        }
    }
    assert!(
        e < -10.0 * KCAL_PER_MOL,
        "a real solvation energy: {}",
        kcal(e)
    );

    let dec = Decoupling::new(&ff, &mask).unwrap();
    assert!(
        dec.rest().unwrap().solvation().is_none(),
        "the rest is in vacuum"
    );
    let mut forces = vec![[0.0; 3]; n];
    let total = dec.energy_and_forces(at, Lambda::COUPLED, &mut forces);
    let reference = ff.evaluate(at);
    let w = reference.energy;
    let (vdw, elec) = b.cross_terms_at(at);
    let magnitude = w.bond.abs()
        + w.angle.abs()
        + w.torsion.abs()
        + w.inversion.abs()
        + w.van_der_waals.abs()
        + w.electrostatic.abs()
        + w.solvation.abs()
        + 2.0 * (vdw.abs() + elec.abs());
    let terms = (ff.stretches().len()
        + ff.bends().len()
        + ff.torsions().len()
        + ff.inversions().len()
        + 2 * ff.pairs().len()
        + n * n) as f64;
    let allowance = terms * EPS * magnitude;
    println!(
        "coupled: {:.9} against the solvated force field's {:.9} kcal/mol, gap {:.1e} against \
         {:.1e}; ΔG_GB {:.4} kcal/mol",
        kcal(total),
        kcal(w.total),
        kcal((total - w.total).abs()),
        kcal(allowance),
        kcal(e)
    );
    assert!((total - w.total).abs() <= allowance);
    let largest = reference
        .forces
        .iter()
        .flatten()
        .fold(0.0f64, |m, f| m.max(f.abs()));
    for (f, g) in forces.iter().zip(&reference.forces) {
        for c in 0..3 {
            assert!((f[c] - g[c]).abs() <= terms * EPS * 10.0 * largest);
        }
    }
}

/// **With the charges off, the group's charges are in no term**: at `λ_e = 0` (and `λ_v = 1`)
/// the solvation term is `accumulate` with benzene's charges set to zero, energy and forces, and
/// the coupling energy is the cross van der Waals plus that, with no cross Coulomb in it — each
/// to the bit. Its `∂/∂λ_e` there is the cross Coulomb plus the GB terms linear in the group's
/// charges, `−k Σ_{i∈L, j∉L} q_i q_j g(f_ij)` over both orderings, rebuilt here from the radii.
#[test]
fn with_the_charges_off_the_group_carries_no_charge() {
    let b = small_binding();
    let ff = b.force_field().clone().with_generalized_born();
    let gb = ff.solvation().unwrap();
    let at = b.positions();
    let n = at.len();
    let mask = ligand_mask(b);
    let zeroed: Vec<f64> = b
        .charges()
        .iter()
        .zip(&mask)
        .map(|(&q, &l)| if l { 0.0 } else { q })
        .collect();
    let mut want = vec![[0.0; 3]; n];
    let e = gb.accumulate(&zeroed, at, &mut want);
    let mut got = vec![[0.0; 3]; n];
    let d = gb.decoupled(b.charges(), &mask, at, 0.0, 1.0, &mut got);
    assert_eq!(d.energy.to_bits(), e.to_bits());
    same_forces(&got, &want, "λ_e = 0");

    let dec = Decoupling::new(&ff, &mask).unwrap();
    let (vdw, elec) = b.cross_terms_at(at);
    let c = dec.coupling(at, Lambda::new(0.0, 0.0, 1.0));
    assert_eq!(c.energy.to_bits(), (vdw + e).to_bits());

    // The linear GB terms, from the radii the model reports and eq 3.
    let born = gb.born_radii(at);
    let k = coulomb_si();
    let q = b.charges();
    let mut linear = 0.0;
    let mut absolute = 0.0;
    for i in 0..n {
        for j in 0..n {
            if mask[i] && !mask[j] {
                let r = distance(at[i], at[j]);
                let rr = born[i] * born[j];
                let f = (r * r + rr * (-r * r / (4.0 * rr)).exp()).sqrt();
                let t = -k * q[i] * q[j] * (1.0 - 1.0 / WATER_DIELECTRIC) / f;
                linear += t;
                absolute += t.abs();
            }
        }
    }
    let want = elec + linear;
    let pairs = (b.pocket_len() * (n - b.pocket_len())) as f64;
    assert!(
        (c.gradient[1] - want).abs() <= 4.0 * pairs * EPS * (absolute + elec.abs()),
        "∂U/∂λ_e at 0: {:e} against {:e}",
        c.gradient[1],
        want
    );
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

/// **Fully decoupled, it is the pocket alone in solvent and benzene alone in vacuum.** At
/// `λ_e = λ_v = 0` the solvation term is the pocket's own GB — its own atoms, charges and
/// positions — to the bit, its force on every pocket atom the pocket's own to the bit and on every
/// benzene atom zero; and the whole energy is the pocket's solvated force field plus benzene's
/// vacuum force field, to their rounding (`n ε Σ|t|`), with every force the same to `n ε` times
/// the largest. The restraint is off and its λ irrelevant here.
#[test]
fn fully_decoupled_it_is_the_pocket_in_solvent_and_benzene_in_vacuum() {
    let b = small_binding();
    let ff = b.force_field().clone().with_generalized_born();
    let gb = ff.solvation().unwrap();
    let at = b.positions();
    let n = at.len();
    let n0 = b.pocket_len();
    let mask = ligand_mask(b);
    let pocket = GeneralizedBorn::new(&b.elements()[..n0]);
    let mut want = vec![[0.0; 3]; n0];
    let e = pocket.accumulate(&b.charges()[..n0], &at[..n0], &mut want);
    let mut got = vec![[0.0; 3]; n];
    let d = gb.decoupled(b.charges(), &mask, at, 0.0, 0.0, &mut got);
    assert_eq!(d.energy.to_bits(), e.to_bits());
    same_forces(&got[..n0], &want, "the pocket");
    assert!(got[n0..].iter().flatten().all(|&f| f == 0.0));
    // Not zero: the pocket's radii still answer to benzene's volume coming back.
    assert!(d.d_van_der_waals != 0.0);

    let dec = Decoupling::new(&ff, &mask).unwrap();
    let mut forces = vec![[0.0; 3]; n];
    let total = dec.energy_and_forces(at, Lambda::new(0.0, 0.0, 0.0), &mut forces);
    let pocket_ff = b.pocket_force_field().clone().with_generalized_born();
    let p = pocket_ff.evaluate(&at[..n0]);
    let l = b.ligand_force_field().evaluate(&at[n0..]);
    let fields = |x: &pantometry_forcefield::Energy| {
        x.bond.abs()
            + x.angle.abs()
            + x.torsion.abs()
            + x.inversion.abs()
            + x.van_der_waals.abs()
            + x.electrostatic.abs()
            + x.solvation.abs()
    };
    let terms = (ff.stretches().len()
        + ff.bends().len()
        + ff.torsions().len()
        + ff.inversions().len()
        + ff.pairs().len()
        + n * n) as f64;
    let allowance = terms * EPS * (fields(&p.energy) + fields(&l.energy));
    let want = p.energy.total + l.energy.total;
    println!(
        "decoupled: {:.9} against pocket in GB + benzene in vacuum {:.9} kcal/mol, gap {:.1e} \
         against {:.1e}",
        kcal(total),
        kcal(want),
        kcal((total - want).abs()),
        kcal(allowance)
    );
    assert!((total - want).abs() <= allowance);
    let reference: Vec<[f64; 3]> = p.forces.iter().chain(&l.forces).copied().collect();
    let largest = reference
        .iter()
        .flatten()
        .fold(0.0f64, |m, f| m.max(f.abs()));
    for (f, g) in forces.iter().zip(&reference) {
        for c in 0..3 {
            assert!((f[c] - g[c]).abs() <= terms * EPS * 10.0 * largest);
        }
    }
}

/// **The forces are the energy's gradient and `∂U/∂λ` its derivative in each λ**, at three states
/// that put each GB derivative to work — charges part-way with the volume present (where `∂/∂λ_v`
/// is taken one-sided, second order), both part-way, and the volume part-way with the charges off — on the 3 Å complex with benzene moved 0.3 Å off
/// its pose so that no term sits at a minimum. Central differences at 1e-5 Å and its half, each
/// within 1e-6 of the largest force on the atoms checked (benzene's and the six pocket atoms
/// nearest it); in λ at 1e-5, within 1e-6 of the derivative. The truncation, `h² U‴/6`, and the
/// rounding, `ε|U|/h`, are each orders below that.
#[test]
fn the_forces_and_the_lambda_derivatives_are_the_energys() {
    let b = small_binding();
    let ff = b.force_field().clone().with_generalized_born();
    let mask = ligand_mask(b);
    let dec = Decoupling::new(&ff, &mask).unwrap();
    let n0 = b.pocket_len();
    let mut at = b.positions().to_vec();
    for p in &mut at[n0..] {
        p[0] += 0.3 * ANGSTROM;
    }
    let centre = at[n0];
    let mut near: Vec<usize> = (0..n0).collect();
    near.sort_by(|&i, &j| {
        distance(at[i], centre)
            .partial_cmp(&distance(at[j], centre))
            .unwrap()
    });
    let mut atoms: Vec<usize> = (n0..at.len()).collect();
    atoms.extend(&near[..6]);
    for l in [
        Lambda::new(0.0, 0.4, 1.0),
        Lambda::new(0.0, 0.6, 0.7),
        Lambda::new(0.0, 0.0, 0.35),
    ] {
        let mut forces = vec![[0.0; 3]; at.len()];
        dec.energy_and_forces(&at, l, &mut forces);
        let largest = atoms
            .iter()
            .flat_map(|&k| forces[k])
            .fold(0.0f64, |m, f| m.max(f.abs()));
        let energy = |p: &[[f64; 3]]| {
            let mut scratch = vec![[0.0; 3]; p.len()];
            dec.energy_and_forces(p, l, &mut scratch)
        };
        for &k in &atoms {
            for c in 0..3 {
                for h in [1e-5 * ANGSTROM, 0.5e-5 * ANGSTROM] {
                    let (mut plus, mut minus) = (at.clone(), at.clone());
                    plus[k][c] += h;
                    minus[k][c] -= h;
                    let fd = (energy(&plus) - energy(&minus)) / (2.0 * h);
                    let err = (fd + forces[k][c]).abs();
                    assert!(
                        err <= 1e-6 * largest,
                        "{l:?} atom {k} {c}: {err:e} of {largest:e}"
                    );
                }
            }
        }
        let g = dec.coupling(&at, l).gradient;
        for c in 1..3 {
            let shifted = |x: f64| {
                let mut m = l.components();
                m[c] = x;
                dec.coupling(&at, Lambda::new(m[0], m[1], m[2])).energy
            };
            let x = l.components()[c];
            let h = 1e-5;
            // Central inside, and the second-order one-sided difference at an end of [0, 1].
            let fd = if x + h > 1.0 {
                (3.0 * shifted(x) - 4.0 * shifted(x - h) + shifted(x - 2.0 * h)) / (2.0 * h)
            } else if x - h < 0.0 {
                (-3.0 * shifted(x) + 4.0 * shifted(x + h) - shifted(x + 2.0 * h)) / (2.0 * h)
            } else {
                (shifted(x + h) - shifted(x - h)) / (2.0 * h)
            };
            assert!(
                (fd - g[c]).abs() <= 1e-6 * g[c].abs().max(1e-3 * KCAL_PER_MOL),
                "{l:?} ∂U/∂λ[{c}]: {fd:e} against {:e}",
                g[c]
            );
        }
    }
}

/// **The frozen buffer's kept descreening is still used, changes no bit, and is dropped when a
/// frozen atom moves.** A complex with only benzene mobile, under OBC II: its decoupling's model
/// holds the pocket as frozen; at the start and with benzene moved its kept terms are used, and
/// with one pocket atom moved by 0.1 Å they are not; in each case the decoupled energy, both
/// derivatives and every force at four states are the same model's without the kept terms, to the
/// bit. **And the states of one configuration evaluated together** (`Alchemical::couplings`,
/// which shares the descreening) are each the state evaluated alone, to the bit.
#[test]
fn the_frozen_cache_is_used_and_changes_no_bit() {
    let b = small_binding();
    let mask = ligand_mask(b);
    let c = Complex::with_mobile(b, mask.clone(), Solvent::GeneralizedBorn).expect("the zone");
    let dec = Decoupling::new(c.potential(), &mask).unwrap();
    let kept = dec.solvation().expect("solvated");
    assert_eq!(kept.frozen_count(), b.pocket_len());
    let bare = kept.clone().without_frozen();
    let start = c.positions().to_vec();
    let mut moved = start.clone();
    for p in &mut moved[b.pocket_len()..] {
        p[1] += 0.2 * ANGSTROM;
    }
    let mut stale = moved.clone();
    stale[3][2] += 0.1 * ANGSTROM;
    let states = [[1.0, 1.0], [0.5, 1.0], [0.0, 0.6], [0.0, 0.0]];
    for (at, reuses) in [(&start, true), (&moved, true), (&stale, false)] {
        assert_eq!(kept.reuses_frozen_terms(at), reuses);
        for [le, lv] in states {
            let n = at.len();
            let (mut f1, mut f2) = (vec![[0.0; 3]; n], vec![[0.0; 3]; n]);
            let x = kept.decoupled(b.charges(), &mask, at, le, lv, &mut f1);
            let y = bare.decoupled(b.charges(), &mask, at, le, lv, &mut f2);
            assert_eq!(x.energy.to_bits(), y.energy.to_bits());
            assert_eq!(x.d_electrostatics.to_bits(), y.d_electrostatics.to_bits());
            assert_eq!(x.d_van_der_waals.to_bits(), y.d_van_der_waals.to_bits());
            same_forces(&f1, &f2, "kept against bare");
        }
        let lambdas = [
            Lambda::new(0.0, 1.0, 1.0),
            Lambda::new(1.0, 1.0, 1.0),
            Lambda::new(1.0, 0.25, 1.0),
            Lambda::new(1.0, 0.0, 0.8),
            Lambda::new(1.0, 0.0, 0.8),
            Lambda::new(1.0, 0.0, 0.0),
        ];
        let together = dec.couplings(at, &lambdas);
        for (l, t) in lambdas.iter().zip(&together) {
            let alone = dec.coupling(at, *l);
            assert_eq!(t.energy.to_bits(), alone.energy.to_bits(), "{l:?}");
            for k in 0..3 {
                assert_eq!(t.gradient[k].to_bits(), alone.gradient[k].to_bits());
            }
        }
    }
}

/// A Born ion: one atom of intrinsic radius `rho` and scale `scale`, charge `q`, decoupled from the
/// solvent alone.
fn born_ion(rho: f64, scale: f64, q: f64) -> Decoupling {
    Decoupling::from_pairs(1, vec![true], Vec::new())
        .with_solvation(GeneralizedBorn::from_radii(vec![rho], vec![scale]), vec![q])
}

/// **A Born ion decouples to its Born energy, by TI and by BAR, through `Windows`.** A lone atom's
/// descreening integral is zero, so its Born radius is `ρ̃ = ρ − 0.09 Å` exactly (OBC's `tanh 0 = 0`)
/// and its GB energy `G = −(q²/2ρ̃)(1 − 1/ε) × 332.0637 kcal/mol` at every position. With its charge
/// scaled, `U(λ_e) = λ_e² G`, so the integrand `∂U/∂λ_e = 2λ_e G` is linear in λ and **the
/// trapezoid is exact**, as Simpson is; every sample is the same number, so BAR's `Δu` is a
/// constant and its root is that constant. Five states of `λ_e` from 1 to 0: the decoupling free
/// energy is `−G` to the rounding of a few sums — five weighted terms for TI, four roots for BAR,
/// each a few roundings of `G`'s size — so to 8 ε of `G` (measured: at most 1.69 ε). A
/// chloride-sized ion (ρ 1.7 Å) of charge −1 and one of charge +0.5 and ρ 1.2 Å. And `∂U/∂λ_v` is
/// zero: there is nothing to descreen.
///
/// **The totals cannot see the path**, only its ends: a Born self term linear in `λ_e` with its
/// derivative consistent, `U = λ_e G`, integrates to the same `−G`. **What tells the path is the
/// per-sample check** that each window's `∂U/∂λ_e` is `2 λ_e G`, the square of a scaled charge's
/// derivative, to 4 ε — and that sabotage is caught there, and only there.
#[test]
fn a_born_ion_decouples_to_its_born_energy() {
    for (rho, q) in [(1.7 * ANGSTROM, -1.0), (1.2 * ANGSTROM, 0.5)] {
        let ion = born_ion(rho, 0.8, q);
        let schedule: Vec<Lambda> = [1.0, 0.75, 0.5, 0.25, 0.0]
            .iter()
            .map(|&e| Lambda::new(0.0, e, 1.0))
            .collect();
        let protocol = Protocol {
            time_step: 1.0 * FS,
            temperature: KELVIN,
            friction: 1e12,
            equilibration: 10,
            stride: 5,
            samples: 20,
            seed: 0x3C03,
        };
        let mut w = Windows::new(
            &ion,
            schedule,
            &[[0.0; 3]],
            vec![Element::Cl.mass()],
            vec![false],
            protocol,
        );
        assert!(w.run(&ion, u64::MAX));
        let born =
            -0.5 * coulomb_si() * q * q * (1.0 - 1.0 / WATER_DIELECTRIC) / (rho - RADIUS_OFFSET);
        let exact = -born;
        let ti = w.thermodynamic_integration(Quadrature::Trapezoid).unwrap();
        let simpson = w.thermodynamic_integration(Quadrature::Simpson).unwrap();
        let bar = w.bennett_total();
        println!(
            "Born ion ρ {:.1} Å, q {q:+}: −G = {:.9} kcal/mol; TI {:.9}, Simpson {:.9}, BAR {:.9}",
            rho / ANGSTROM,
            kcal(exact),
            kcal(ti.value),
            kcal(simpson.value),
            kcal(bar.value)
        );
        for (what, v) in [
            ("trapezoid", ti.value),
            ("Simpson", simpson.value),
            ("BAR", bar.value),
        ] {
            assert!(
                (v - exact).abs() <= 8.0 * EPS * exact.abs(),
                "{what}: {:e} against {:e}",
                v,
                exact
            );
        }
        for win in w.windows() {
            for s in win.samples() {
                assert_eq!(s.gradient[2], 0.0);
                let want = 2.0 * win.lambda().electrostatics * born;
                assert!((s.gradient[1] - want).abs() <= 4.0 * EPS * born.abs());
            }
        }
    }
}

/// **The measurements' interval BAR is `Windows::bennett`'s, to the bit**: the Born ion's windows,
/// read into the records the measurements analyse, give each interval's `Bar` — estimate, variance
/// and overlap — with the same bits. The ion's `Δu` is a nonzero constant, so a reverse set of the
/// wrong sign moves the root.
#[test]
fn the_interval_bar_is_the_windows_own() {
    let ion = born_ion(1.7 * ANGSTROM, 0.8, -1.0);
    let schedule: Vec<Lambda> = [1.0, 0.5, 0.0]
        .iter()
        .map(|&e| Lambda::new(0.0, e, 1.0))
        .collect();
    let protocol = Protocol {
        time_step: 1.0 * FS,
        temperature: KELVIN,
        friction: 1e12,
        equilibration: 10,
        stride: 5,
        samples: 20,
        seed: 0x3C03,
    };
    let mut w = Windows::new(
        &ion,
        schedule.clone(),
        &[[0.0; 3]],
        vec![Element::Cl.mass()],
        vec![false],
        protocol,
    );
    assert!(w.run(&ion, u64::MAX));
    let n = schedule.len();
    // Candidate k is window k; each record holds its energy relative to its own state, 0, at
    // itself and at its neighbours.
    let records: Vec<Recorded> = w
        .windows()
        .iter()
        .enumerate()
        .map(|(k, win)| {
            let states: Vec<usize> = (k.saturating_sub(1)..(k + 2).min(n)).collect();
            let energies = win
                .samples()
                .iter()
                .map(|s| {
                    states
                        .iter()
                        .map(|&g| match g.cmp(&k) {
                            std::cmp::Ordering::Less => s.to_previous,
                            std::cmp::Ordering::Equal => 0.0,
                            std::cmp::Ordering::Greater => s.to_next,
                        })
                        .collect()
                })
                .collect();
            Recorded {
                grid: k,
                lambda: win.lambda(),
                states,
                gradients: win.samples().iter().map(|s| s.gradient).collect(),
                energies,
                seconds: 0.0,
            }
        })
        .collect();
    let want = w.bennett();
    let kt = protocol.kt();
    for k in 0..n - 1 {
        let got = interval_bar(&records[k], &records[k + 1], kt);
        assert!(
            want[k].delta.abs() > 1.0,
            "a Δu far from zero: {}",
            want[k].delta
        );
        assert_eq!(got.delta.to_bits(), want[k].delta.to_bits(), "interval {k}");
        assert_eq!(got.variance.to_bits(), want[k].variance.to_bits());
        assert_eq!(got.overlap.to_bits(), want[k].overlap.to_bits());
    }
}

/// **The cycle has its signs**: with decoupling free energies of +2 in solvent and +10 in the
/// complex, and a release of −7 (negative, as every useful restraint's is), the binding free energy
/// is `2 − 10 + 7 = −1`; a release counted with the other sign would give −15, and the legs
/// swapped +15.
#[test]
fn the_cycle_has_its_signs() {
    assert_eq!(binding_free_energy(2.0, 10.0, -7.0), -1.0);
    assert_eq!(binding_free_energy(0.0, 0.0, -7.0), 7.0);
    assert_eq!(binding_free_energy(3.0, 0.0, 0.0), 3.0);
    assert_eq!(binding_free_energy(0.0, 4.0, 0.0), -4.0);
}

/// OBC II's radius (eq 6, 8) for intrinsic radius `rho` and integral `i`, written out here.
fn obc_ii(rho: f64, i: f64) -> f64 {
    let rt = rho - RADIUS_OFFSET;
    let psi = i * rt;
    let t = psi - 0.8 * psi * psi + 4.85 * psi * psi * psi;
    1.0 / (1.0 / rt - t.tanh() / rho)
}

/// **Two ions follow eqs 2–8 written out, at every `(λ_e, λ_v)`.** Ion A (the group: ρ 1.7 Å,
/// S 0.8, q +0.5) and ion B (ρ 1.5 Å, S 0.85, q −0.7) 2.6 Å apart, so that each one's sphere
/// reaches into the other's. Each radius is OBC II's with `I = λ_v ×` the pair integral of the
/// other (the crate's `descreening`, itself checked against quadrature elsewhere), and
/// `G = −½k[(λ_e q_A)² g(R_A) + q_B² g(R_B)] − k λ_e q_A q_B g(f)`, `g(x) = (1 − 1/ε)/x`. Energy
/// and `∂/∂λ_e` to 1e-13 of themselves, a few roundings of a short sum. `∂/∂λ_v` against a
/// difference of the written-out energy at `h` = 1e-4 — central inside [0, 1], second-order
/// one-sided at its ends, as the benzene test takes them — within its own error: the truncation
/// `h² M₃/6` (central) or `h² M₃/3` (one-sided), with `M₃` the largest `|∂³G/∂λ_v³|` on a 0.01
/// grid by a third difference, doubled for what the grid misses, plus the rounding `4 ε max|G|/h`.
///
/// **And with both ions frozen**, so that the one kept descreening term crosses the partition, the
/// decoupled energy, both derivatives and both forces are the uncached model's to the bit at every
/// state.
#[test]
fn two_ions_follow_the_closed_form_at_every_state() {
    let (ra, sa, qa) = (1.7 * ANGSTROM, 0.8, 0.5);
    let (rb, sb, qb) = (1.5 * ANGSTROM, 0.85, -0.7);
    let r = 2.6 * ANGSTROM;
    let at = [[0.0; 3], [r, 0.0, 0.0]];
    let gb = GeneralizedBorn::from_radii(vec![ra, rb], vec![sa, sb]);
    let group = [true, false];
    let k = coulomb_si();
    let g = |x: f64| (1.0 - 1.0 / WATER_DIELECTRIC) / x;
    let closed = |le: f64, lv: f64| -> (f64, f64) {
        let ia = lv * descreening(r, ra - RADIUS_OFFSET, sb * (rb - RADIUS_OFFSET)).0;
        let ib = lv * descreening(r, rb - RADIUS_OFFSET, sa * (ra - RADIUS_OFFSET)).0;
        let (big_a, big_b) = (obc_ii(ra, ia), obc_ii(rb, ib));
        let rr = big_a * big_b;
        let f = (r * r + rr * (-r * r / (4.0 * rr)).exp()).sqrt();
        let e = -0.5 * k * ((le * qa) * (le * qa) * g(big_a) + qb * qb * g(big_b))
            - k * le * qa * qb * g(f);
        let de = -k * le * qa * qa * g(big_a) - k * qa * qb * g(f);
        (e, de)
    };
    assert!(descreening(r, rb - RADIUS_OFFSET, sa * (ra - RADIUS_OFFSET)).0 > 0.0);
    for le in [0.0, 0.3, 1.0] {
        for lv in [0.0, 0.45, 1.0] {
            let mut f = [[0.0; 3]; 2];
            let d = gb.decoupled(&[qa, qb], &group, &at, le, lv, &mut f);
            let (e, de) = closed(le, lv);
            assert!(
                (d.energy - e).abs() <= 1e-13 * e.abs(),
                "({le}, {lv}): {e:e}"
            );
            assert!(
                (d.d_electrostatics - de).abs() <= 1e-13 * de.abs(),
                "({le}, {lv}): {de:e}"
            );
            // The third derivative's bound on [0, 1], and the energy's size, from the closed form.
            let u = |x: f64| closed(le, x).0;
            let grid = 0.01;
            let m3 = 2.0
                * (0..97)
                    .map(|i| {
                        let x = i as f64 * grid;
                        ((u(x + 3.0 * grid) - 3.0 * u(x + 2.0 * grid) + 3.0 * u(x + grid) - u(x))
                            / (grid * grid * grid))
                            .abs()
                    })
                    .fold(0.0f64, f64::max);
            let size = (0..=100)
                .map(|i| u(i as f64 * grid).abs())
                .fold(0.0f64, f64::max);
            let h = 1e-4;
            let (dv, truncation) = if lv + h > 1.0 {
                (
                    (3.0 * u(lv) - 4.0 * u(lv - h) + u(lv - 2.0 * h)) / (2.0 * h),
                    h * h * m3 / 3.0,
                )
            } else if lv - h < 0.0 {
                (
                    (-3.0 * u(lv) + 4.0 * u(lv + h) - u(lv + 2.0 * h)) / (2.0 * h),
                    h * h * m3 / 3.0,
                )
            } else {
                ((u(lv + h) - u(lv - h)) / (2.0 * h), h * h * m3 / 6.0)
            };
            let allowance = truncation + 4.0 * EPS * size / h;
            assert!(
                (d.d_van_der_waals - dv).abs() <= allowance,
                "({le}, {lv}): ∂/∂λ_v {:e} against {dv:e}, allowance {allowance:e}",
                d.d_van_der_waals
            );

            // Both frozen: the one kept term crosses the partition.
            let kept = gb.clone().with_frozen(&[true, true], &at);
            assert!(kept.reuses_frozen_terms(&at));
            let mut fk = [[0.0; 3]; 2];
            let dk = kept.decoupled(&[qa, qb], &group, &at, le, lv, &mut fk);
            assert_eq!(dk.energy.to_bits(), d.energy.to_bits(), "({le}, {lv})");
            assert_eq!(dk.d_electrostatics.to_bits(), d.d_electrostatics.to_bits());
            assert_eq!(dk.d_van_der_waals.to_bits(), d.d_van_der_waals.to_bits());
            same_forces(&fk, &f, "both ions frozen");
        }
    }
}

/// **The solvent leg's molecule, decoupled from the solvent alone**: benzene's own force field
/// under OBC II with the group every atom. Fully coupled the energy is the solvated force field's
/// to its rounding; at `λ_e = 0` the solvation term is exactly zero, so the energy is the vacuum
/// force field's to the bit; and `∂U/∂λ_v` is exactly zero at every state. **That last is a
/// consequence of the construction, not a finding**: with the group the whole molecule there is no
/// cross pair and no cross descreening, so nothing depends on `λ_v`. What it checks is that the
/// code adds nothing where nothing is.
#[test]
fn the_solvent_legs_molecule_decouples_from_the_solvent_alone() {
    let b = small_binding();
    let n0 = b.pocket_len();
    let ligand = b.ligand_force_field().clone();
    let solvated = ligand.clone().with_generalized_born();
    let at = &b.positions()[n0..];
    let all = vec![true; at.len()];
    let dec = Decoupling::new(&solvated, &all).unwrap();
    assert!(dec.cross_pairs().is_empty());
    let mut f = vec![[0.0; 3]; at.len()];
    let coupled = dec.energy_and_forces(at, Lambda::COUPLED, &mut f);
    let want = solvated.evaluate(at).energy;
    let size = want.bond.abs()
        + want.angle.abs()
        + want.torsion.abs()
        + want.inversion.abs()
        + want.van_der_waals.abs()
        + want.electrostatic.abs()
        + want.solvation.abs();
    assert!((coupled - want.total).abs() <= 200.0 * EPS * size);
    for l in [
        Lambda::new(0.0, 0.0, 1.0),
        Lambda::new(0.0, 0.0, 0.3),
        Lambda::new(0.0, 0.0, 0.0),
    ] {
        let c = dec.coupling(at, l);
        assert_eq!(c.energy, 0.0, "{l:?}");
        let e = dec.energy_and_forces(at, l, &mut f);
        assert_eq!(e.to_bits(), ligand.evaluate(at).energy.total.to_bits());
    }
    for le in [1.0, 0.5, 0.0] {
        assert_eq!(dec.coupling(at, Lambda::new(0.0, le, 0.6)).gradient[2], 0.0);
    }
    println!(
        "benzene's ΔG_GB at the crystal pose with its QEq charges: {:.4} kcal/mol",
        kcal(want.solvation)
    );
}

// ---------------------------------------------------------------------------------------------
// The measurements: the solvent leg, the complex leg, and the cycle
// ---------------------------------------------------------------------------------------------

/// Where the measurements write their windows: `PANTOMETRY_3C3_DIR`, or `pantometry-3c3` in the
/// system's temporary directory. Never the repository.
fn results_dir() -> std::path::PathBuf {
    let dir = std::env::var_os("PANTOMETRY_3C3_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("pantometry-3c3"));
    std::fs::create_dir_all(&dir).expect("the results directory");
    dir
}

/// A lock that makes a second run in the same directory refuse, rather than write into the files
/// the first is writing. Removed when dropped; a run that was killed leaves it, and says so.
struct Lock(std::path::PathBuf);

impl Lock {
    fn take(dir: &std::path::Path, leg: &str) -> Lock {
        let path = dir.join(format!("{leg}.lock"));
        std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap_or_else(|e| {
                panic!(
                    "{} exists ({e}): another run is writing here, or one was killed — remove it \
                     by hand once no run is",
                    path.display()
                )
            });
        Lock(path)
    }
}

impl Drop for Lock {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// Appends a line to the leg's log and prints it.
fn log(dir: &std::path::Path, leg: &str, line: &str) {
    use std::io::Write;
    println!("{line}");
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(dir.join(format!("{leg}.log")))
        .expect("the log");
    writeln!(f, "{line}").expect("the log");
}

/// The candidate states, in path order: the restraint on, the charges off, the van der Waals off.
/// Every window of the complex leg is one of them, and records its energy at every candidate
/// between its two neighbours, so a window inserted at a candidate between two that ran needs
/// neither of them run again.
fn grid() -> Vec<Lambda> {
    let mut g = Vec::new();
    for r in [0.0, 0.05, 0.1, 0.175, 0.25, 0.375, 0.5, 0.75, 1.0] {
        g.push(Lambda::new(r, 1.0, 1.0));
    }
    for e in [0.75, 0.5, 0.25, 0.0] {
        g.push(Lambda::new(1.0, e, 1.0));
    }
    for v in [
        0.95, 0.9, 0.825, 0.75, 0.675, 0.6, 0.525, 0.45, 0.375, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05,
        0.0,
    ] {
        g.push(Lambda::new(1.0, 0.0, v));
    }
    g
}

/// The complex leg's first schedule, as indices into [`grid`]: the restraint on in five states
/// (3c-1 had six), the charges off in two more (3c-1: four), the van der Waals off in eight (3c-1:
/// twelve) — every other one of 3c-1's spacing, roughly, since every interval there overlapped by
/// 0.38–0.50. Fifteen windows.
const COMPLEX_SCHEDULE: [usize; 15] = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28];

/// An interval whose overlap is below this gets a window at the candidate between its two ends.
/// Klimovich, Shirts and Mobley (*J. Comput.-Aided Mol. Des.* **29**, 397 (2015), PMC4420631)
/// find 0.03 tolerable "with enough samples"; a window here has 1000 correlated samples, a few
/// hundred independent ones, so the threshold is set at three times that and more.
const OVERLAP_THRESHOLD: f64 = 0.1;

/// The solvent leg's schedule: the complex leg's charges and van der Waals states, without the
/// restraint, which the solvent leg has none of.
fn solvent_grid() -> Vec<Lambda> {
    grid()
        .into_iter()
        .filter(|l| l.restraint == 1.0)
        .map(|l| Lambda::new(0.0, l.electrostatics, l.van_der_waals))
        .collect()
}

/// What one window recorded, read back from its file.
#[derive(Clone, Debug)]
struct Recorded {
    /// Its candidate's index.
    grid: usize,
    lambda: Lambda,
    /// The candidates whose energies each sample holds, in this order.
    states: Vec<usize>,
    gradients: Vec<[f64; 3]>,
    /// Per sample, the coupling energy at each of `states`, joules.
    energies: Vec<Vec<f64>>,
    seconds: f64,
}

impl Recorded {
    /// Reduced `u(to) − u(own)` per sample.
    fn delta(&self, to: usize, kt: f64) -> Vec<f64> {
        let own = self.column(self.grid);
        let other = self.column(to);
        self.energies
            .iter()
            .map(|e| (e[other] - e[own]) / kt)
            .collect()
    }

    fn column(&self, g: usize) -> usize {
        self.states
            .iter()
            .position(|&s| s == g)
            .unwrap_or_else(|| panic!("window {} did not record candidate {g}", self.grid))
    }

    /// The same window with only the samples in `range`.
    fn slice(&self, range: std::ops::Range<usize>) -> Recorded {
        Recorded {
            gradients: self.gradients[range.clone()].to_vec(),
            energies: self.energies[range].to_vec(),
            ..self.clone()
        }
    }
}

/// BAR between neighbouring windows `a` and `b`, reduced, with the delta method's variance: the
/// forward set is `u_b − u_a` on `a`'s samples, and **the reverse set the same difference,
/// `u_b − u_a`, on `b`'s — minus what `b` recorded towards `a`**. The one place the overlap that
/// decides an insertion, and every interval the analysis prints, are computed.
/// `the_interval_bar_is_the_windows_own` holds it to [`Windows::bennett`] to the bit.
fn interval_bar(a: &Recorded, b: &Recorded, kt: f64) -> pantometry_forcefield::Bar {
    let reverse: Vec<f64> = b.delta(a.grid, kt).iter().map(|x| -x).collect();
    pantometry_forcefield::free_energy::bar_correlated(&a.delta(b.grid, kt), &reverse)
}

/// **The cycle**: `ΔG°_bind = ΔG_solvent − ΔG_complex − ΔG°_release`, any one unit. `solvent` and
/// `complex` are the free energies of decoupling the ligand in each leg — from solvent to vacuum,
/// and from the unrestrained bound complex to the restrained decoupled one — and `release` is
/// `Boresch::release_free_energy`, the free energy of releasing the decoupled ligand's restraint to
/// the standard state, which is negative. Equivalently `+ ΔG°_restrain` with
/// `ΔG°_restrain = −release`. Held by `the_cycle_has_its_signs`.
fn binding_free_energy(solvent: f64, complex: f64, release: f64) -> f64 {
    solvent - complex - release
}

/// A leg's description: the Hamiltonian, where it starts, and how its windows run.
struct Leg<'a> {
    name: &'a str,
    hamiltonian: &'a Decoupling,
    start: Vec<[f64; 3]>,
    masses: Vec<f64>,
    frozen: Vec<bool>,
    protocol: Protocol,
    grid: Vec<Lambda>,
    /// Whether each window records every candidate (cheap legs), or only those between its
    /// neighbours.
    record_all: bool,
}

impl Leg<'_> {
    fn header(&self, g: usize, states: &[usize]) -> String {
        let p = &self.protocol;
        let l = self.grid[g];
        format!(
            "# leg {} candidate {g} lambda {:e} {:e} {:e} protocol {:e} {:e} {:e} {} {} {} {} \
             states {}",
            self.name,
            l.restraint,
            l.electrostatics,
            l.van_der_waals,
            p.time_step,
            p.temperature,
            p.friction,
            p.equilibration,
            p.stride,
            p.samples,
            p.seed,
            states
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<_>>()
                .join(",")
        )
    }

    /// The candidates a window at schedule position `pos` must record.
    fn needed(&self, schedule: &[usize], pos: usize) -> Vec<usize> {
        if self.record_all {
            return (0..self.grid.len()).collect();
        }
        let lo = schedule[pos.saturating_sub(1)];
        let hi = schedule[(pos + 1).min(schedule.len() - 1)];
        (lo..=hi).collect()
    }

    fn path(&self, dir: &std::path::Path, g: usize) -> std::path::PathBuf {
        dir.join(format!("{}_{g:02}.txt", self.name))
    }

    /// Window `g`'s record: read from its file if it is there and was run as this one would be,
    /// with at least the candidates `needed`; otherwise run, written sample by sample to a
    /// `.partial` file, and renamed when complete.
    fn window(&self, dir: &std::path::Path, g: usize, needed: &[usize]) -> Recorded {
        let path = self.path(dir, g);
        if let Ok(text) = std::fs::read_to_string(&path) {
            let r = parse(&text);
            let first = text.lines().next().unwrap_or("");
            let want = self.header(g, &r.states);
            assert_eq!(first, want, "{} was run differently", path.display());
            assert!(
                needed.iter().all(|n| r.states.contains(n)),
                "{} lacks a candidate it now needs",
                path.display()
            );
            assert_eq!(
                r.energies.len(),
                self.protocol.samples,
                "{}",
                path.display()
            );
            return r;
        }
        self.run(dir, g, needed)
    }

    fn run(&self, dir: &std::path::Path, g: usize, states: &[usize]) -> Recorded {
        use pantometry_forcefield::free_energy::window_seed;
        use pantometry_forcefield::{AtLambda, Bath, MolecularDynamics};
        use std::io::Write;
        let p = self.protocol;
        let lambda = self.grid[g];
        let seed = window_seed(p.seed, g);
        let t = std::time::Instant::now();
        let mut md = MolecularDynamics::new(self.masses.clone())
            .with_bath(Bath::Langevin {
                temperature: p.temperature,
                friction: p.friction,
                seed,
            })
            .with_frozen(self.frozen.clone())
            .thermalised(&self.start, p.temperature, seed);
        let potential = AtLambda {
            hamiltonian: self.hamiltonian,
            lambda,
        };
        let mut at = self.start.clone();
        md.prepare(&potential, &at);
        let partial = self.path(dir, g).with_extension("partial");
        let mut file = std::io::BufWriter::new(std::fs::File::create(&partial).expect("partial"));
        writeln!(file, "{}", self.header(g, states)).unwrap();
        let lambdas: Vec<Lambda> = std::iter::once(lambda)
            .chain(states.iter().map(|&s| self.grid[s]))
            .collect();
        let total = p.steps_per_window();
        let mut taken = 0usize;
        log(
            dir,
            self.name,
            &format!(
                "{}: window at candidate {g} λ = ({}, {}, {}) begins, {total} steps, recording {} \
                 states",
                self.name,
                lambda.restraint,
                lambda.electrostatics,
                lambda.van_der_waals,
                states.len()
            ),
        );
        // Clippy on current stable suggests `u64::is_multiple_of`, stabilised in 1.87; this crate
        // builds on 1.78.
        #[allow(clippy::manual_is_multiple_of)]
        while taken < p.samples {
            md.step(&potential, &mut at, p.time_step);
            let s = md.steps();
            if s > p.equilibration && (s - p.equilibration) % p.stride == 0 {
                let c = self.hamiltonian.couplings(&at, &lambdas);
                let mut line = format!(
                    "{s} {:e} {:e} {:e}",
                    c[0].gradient[0], c[0].gradient[1], c[0].gradient[2]
                );
                for x in &c[1..] {
                    line.push_str(&format!(" {:e}", x.energy));
                }
                writeln!(file, "{line}").unwrap();
                taken += 1;
                if taken % 100 == 0 {
                    file.flush().unwrap();
                    log(
                        dir,
                        self.name,
                        &format!(
                            "{}: candidate {g}: {taken} of {} samples, step {s} of {total}, {:.0} s",
                            self.name,
                            p.samples,
                            t.elapsed().as_secs_f64()
                        ),
                    );
                }
            }
        }
        // The buffer kept its bits.
        for (k, &f) in self.frozen.iter().enumerate() {
            if f {
                assert_eq!(at[k], self.start[k], "frozen atom {k}");
            }
        }
        let seconds = t.elapsed().as_secs_f64();
        writeln!(file, "# done seconds {seconds:e}").unwrap();
        drop(file);
        std::fs::rename(&partial, self.path(dir, g)).expect("rename");
        log(
            dir,
            self.name,
            &format!(
                "{}: window at candidate {g} done in {seconds:.0} s ({:.2} ms per step)",
                self.name,
                seconds * 1e3 / total as f64
            ),
        );
        parse(&std::fs::read_to_string(self.path(dir, g)).unwrap())
    }
}

/// A window's file, read back.
fn parse(text: &str) -> Recorded {
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().expect("a header").split_whitespace().collect();
    let field = |name: &str| header.iter().position(|&w| w == name).expect(name) + 1;
    let grid: usize = header[field("candidate")].parse().unwrap();
    let l = field("lambda");
    let lambda = Lambda::new(
        header[l].parse().unwrap(),
        header[l + 1].parse().unwrap(),
        header[l + 2].parse().unwrap(),
    );
    let states: Vec<usize> = header[field("states")]
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect();
    let (mut gradients, mut energies, mut seconds) = (Vec::new(), Vec::new(), f64::NAN);
    for line in lines {
        if let Some(rest) = line.strip_prefix("# done seconds ") {
            seconds = rest.trim().parse().unwrap();
            continue;
        }
        let v: Vec<f64> = line
            .split_whitespace()
            .skip(1)
            .map(|x| x.parse().unwrap())
            .collect();
        assert_eq!(v.len(), 3 + states.len(), "a sample line: {line}");
        gradients.push([v[0], v[1], v[2]]);
        energies.push(v[3..].to_vec());
    }
    assert!(seconds.is_finite(), "a complete window says so");
    Recorded {
        grid,
        lambda,
        states,
        gradients,
        energies,
        seconds,
    }
}

/// Every window of `schedule` (candidate indices, in path order), run or read, then a window at
/// the candidate between each pair whose BAR overlap is below [`OVERLAP_THRESHOLD`], once. Returns
/// the final schedule and its records.
fn run_leg(leg: &Leg<'_>, dir: &std::path::Path, mut schedule: Vec<usize>) -> Vec<Recorded> {
    let _lock = Lock::take(dir, leg.name);
    let kt = leg.protocol.kt();
    let mut records: Vec<Recorded> = (0..schedule.len())
        .map(|pos| leg.window(dir, schedule[pos], &leg.needed(&schedule, pos)))
        .collect();
    let mut inserted = Vec::new();
    for k in 0..schedule.len() - 1 {
        let (a, b) = (&records[k], &records[k + 1]);
        let bar = interval_bar(a, b, kt);
        if bar.overlap < OVERLAP_THRESHOLD {
            let between = (schedule[k] + schedule[k + 1]) / 2;
            log(
                dir,
                leg.name,
                &format!(
                    "{}: overlap {:.3} between candidates {} and {}, below {OVERLAP_THRESHOLD}: {}",
                    leg.name,
                    bar.overlap,
                    schedule[k],
                    schedule[k + 1],
                    if between > schedule[k] {
                        format!("inserting candidate {between}")
                    } else {
                        "no candidate between them".to_string()
                    }
                ),
            );
            if between > schedule[k] {
                inserted.push((k + 1, between));
            }
        }
    }
    for (offset, (pos, g)) in inserted.into_iter().enumerate() {
        schedule.insert(pos + offset, g);
        let needed = leg.needed(&schedule, pos + offset);
        records.insert(pos + offset, leg.window(dir, g, &needed));
    }
    records
}

/// A leg's numbers, in kcal/mol: BAR's and TI's totals with their errors, and BAR by segment.
struct Analysis {
    bar: (f64, f64),
    ti: (f64, f64),
    segments: Vec<(f64, f64)>,
}

/// BAR along `records`, reduced, by the delta method window by window.
fn chain(records: &[Recorded], kt: f64) -> pantometry_forcefield::FreeEnergy {
    let n = records.len();
    let windows: Vec<(Vec<f64>, Vec<f64>)> = (0..n)
        .map(|j| {
            let to_previous = if j > 0 {
                records[j].delta(records[j - 1].grid, kt)
            } else {
                Vec::new()
            };
            let to_next = if j + 1 < n {
                records[j].delta(records[j + 1].grid, kt)
            } else {
                Vec::new()
            };
            (to_previous, to_next)
        })
        .collect();
    pantometry_forcefield::free_energy::bennett_chain(&windows)
}

/// `ln mean(exp(x))` without overflow.
fn log_mean_exp(xs: &[f64]) -> f64 {
    let m = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    m + (xs.iter().map(|x| (x - m).exp()).sum::<f64>() / xs.len() as f64).ln()
}

/// Prints the leg's table and returns its numbers. `segments` names runs of windows by their
/// first and last schedule positions.
fn analyse(
    dir: &std::path::Path,
    leg: &str,
    records: &[Recorded],
    kt: f64,
    segments: &[(&str, usize, usize)],
) -> Analysis {
    use pantometry_forcefield::Estimate;
    let n = records.len();
    let to_kcal = |x: f64| kcal(kt * x);
    // TI on the trapezoid along the path.
    let weights: Vec<[f64; 3]> = (0..n)
        .map(|k| {
            let lo = records[k.saturating_sub(1)].lambda.components();
            let hi = records[(k + 1).min(n - 1)].lambda.components();
            [0, 1, 2].map(|c| 0.5 * (hi[c] - lo[c]))
        })
        .collect();
    let terms: Vec<Estimate> = records
        .iter()
        .zip(&weights)
        .map(|(r, w)| {
            let ys: Vec<f64> = r
                .gradients
                .iter()
                .map(|g| (0..3).map(|c| g[c] * w[c]).sum())
                .collect();
            Estimate::of(&ys)
        })
        .collect();
    let ti = (
        kcal(terms.iter().map(|t| t.mean).sum()),
        kcal(terms.iter().map(|t| t.error * t.error).sum::<f64>().sqrt()),
    );
    log(dir, leg, &format!("{leg}: the windows"));
    log(
        dir,
        leg,
        "| window | candidate | λ_r | λ_e | λ_v | TI term (kcal/mol) | τ (samples) | wall (s) | → next: BAR (kcal/mol) | overlap | EXP fwd | EXP rev |",
    );
    log(
        dir,
        leg,
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    );
    let (mut fwd, mut rev) = (0.0, 0.0);
    for k in 0..n {
        let r = &records[k];
        let next = if k + 1 < n {
            let s = &records[k + 1];
            let f = r.delta(s.grid, kt);
            let reverse: Vec<f64> = s.delta(r.grid, kt).iter().map(|x| -x).collect();
            let bar = interval_bar(r, s, kt);
            let neg: Vec<f64> = f.iter().map(|x| -x).collect();
            let ef = -log_mean_exp(&neg);
            let er = log_mean_exp(&reverse);
            fwd += ef;
            rev += er;
            format!(
                "{:+.3} ± {:.3} | {:.3} | {:+.3} | {:+.3}",
                to_kcal(bar.delta),
                to_kcal(bar.error()),
                bar.overlap,
                to_kcal(ef),
                to_kcal(er)
            )
        } else {
            "| | |".to_string()
        };
        log(
            dir,
            leg,
            &format!(
                "| {k} | {} | {} | {} | {} | {:+.3} ± {:.3} | {:.1} | {:.0} | {next} |",
                r.grid,
                r.lambda.restraint,
                r.lambda.electrostatics,
                r.lambda.van_der_waals,
                kcal(terms[k].mean),
                kcal(terms[k].error),
                terms[k].tau,
                r.seconds
            ),
        );
    }
    let total = chain(records, kt);
    let bar = (to_kcal(total.value), to_kcal(total.error));
    let mut seg = Vec::new();
    for &(name, a, b) in segments {
        let c = chain(&records[a..=b], kt);
        let ti_seg: f64 = {
            // The segment's TI with its own end weights.
            let sub = &records[a..=b];
            let m = sub.len();
            (0..m)
                .map(|k| {
                    let lo = sub[k.saturating_sub(1)].lambda.components();
                    let hi = sub[(k + 1).min(m - 1)].lambda.components();
                    let w = [0, 1, 2].map(|c| 0.5 * (hi[c] - lo[c]));
                    let ys: Vec<f64> = sub[k]
                        .gradients
                        .iter()
                        .map(|g| (0..3).map(|c| g[c] * w[c]).sum())
                        .collect();
                    Estimate::of(&ys).mean
                })
                .sum()
        };
        log(
            dir,
            leg,
            &format!(
                "{leg}: {name}: BAR {:+.3} ± {:.3}, TI {:+.3} kcal/mol",
                to_kcal(c.value),
                to_kcal(c.error),
                kcal(ti_seg)
            ),
        );
        seg.push((to_kcal(c.value), to_kcal(c.error)));
    }
    let half = records[0].energies.len() / 2;
    let first: Vec<Recorded> = records.iter().map(|r| r.slice(0..half)).collect();
    let second: Vec<Recorded> = records
        .iter()
        .map(|r| r.slice(half..r.energies.len()))
        .collect();
    let (h1, h2) = (chain(&first, kt), chain(&second, kt));
    let halves = [
        (to_kcal(h1.value), to_kcal(h1.error)),
        (to_kcal(h2.value), to_kcal(h2.error)),
    ];
    let exp = (to_kcal(fwd), to_kcal(rev));
    log(
        dir,
        leg,
        &format!(
            "{leg}: BAR {:+.3} ± {:.3}, TI (trapezoid) {:+.3} ± {:.3} kcal/mol; first half of \
             the samples {:+.3} ± {:.3}, second half {:+.3} ± {:.3}; EXP forward {:+.3}, reverse \
             {:+.3}",
            bar.0,
            bar.1,
            ti.0,
            ti.1,
            halves[0].0,
            halves[0].1,
            halves[1].0,
            halves[1].1,
            exp.0,
            exp.1
        ),
    );
    Analysis {
        bar,
        ti,
        segments: seg,
    }
}

/// The 10 Å binding with its hydrogens relaxed, as 3b and 3c-1 start.
fn relaxed_binding() -> Binding {
    Binding::new(the_system(), 10.0 * ANGSTROM)
        .expect("the binding")
        .relaxing_hydrogens(20_000, Binding::HYDROGEN_TOLERANCE)
}

/// The Boresch force constants of 3c-1: 20 kcal mol⁻¹ Å⁻² and rad⁻² in `½ K (ξ − ξ₀)²`, Mobley,
/// Chodera and Dill's K₀ = 10 in `K₀ (ξ − ξ₀)²`.
fn force_constants() -> [f64; 6] {
    let r = 20.0 * KCAL_PER_MOL / (ANGSTROM * ANGSTROM);
    let a = 20.0 * KCAL_PER_MOL;
    [r, a, a, a, a, a]
}

/// 3c-1's anchor rule, the same code: receptor candidates `(N, CA, C)` and `(C, CA, N)` of every
/// pocket residue, ligand candidates every path of three bonded heavy atoms, `r₀` in [3, 8] Å.
fn anchors(b: &Binding, at: &[[f64; 3]]) -> pantometry_forcefield::Boresch {
    let s = the_system();
    let name = |k: usize| s.atom_name(b.system_atoms()[k]);
    let mut receptor = Vec::new();
    for span in b.residue_spans() {
        let find = |n: &str| span.clone().find(|&k| name(k) == n);
        if let (Some(n), Some(ca), Some(c)) = (find("N"), find("CA"), find("C")) {
            receptor.push([n, ca, c]);
            receptor.push([c, ca, n]);
        }
    }
    let n0 = b.pocket_len();
    let heavy = |k: usize| b.elements()[k] != Element::H;
    let bonds: Vec<[usize; 2]> = b
        .ligand_force_field()
        .stretches()
        .iter()
        .map(|t| [t.atoms[0] + n0, t.atoms[1] + n0])
        .filter(|&[i, j]| heavy(i) && heavy(j))
        .collect();
    let neighbours = |k: usize| -> Vec<usize> {
        bonds
            .iter()
            .filter_map(|&[i, j]| {
                if i == k {
                    Some(j)
                } else if j == k {
                    Some(i)
                } else {
                    None
                }
            })
            .collect()
    };
    let mut ligand = Vec::new();
    for a in b.ligand_range().filter(|&k| heavy(k)) {
        for bb in neighbours(a) {
            for c in neighbours(bb) {
                if c != a {
                    ligand.push([a, bb, c]);
                }
            }
        }
    }
    pantometry_forcefield::Boresch::choose(
        at,
        &receptor,
        &ligand,
        (3.0 * ANGSTROM, 8.0 * ANGSTROM),
        force_constants(),
    )
    .expect("an anchor pair within 3–8 Å")
}

/// The solvent leg's protocol: four seeds, each 20 ps discarded and 100 ps sampled every 20 fs
/// per window — affordable because benzene alone costs microseconds a step.
fn solvent_protocol(seed: u64) -> Protocol {
    Protocol {
        time_step: 0.5 * FS,
        temperature: KELVIN,
        friction: 1e12,
        equilibration: 40_000,
        stride: 40,
        samples: 5000,
        seed,
    }
}

const SOLVENT_SEEDS: [u64; 4] = [0x3C3_501, 0x3C3_502, 0x3C3_503, 0x3C3_504];

/// **Benzene decoupled from OBC II water alone**, the solvent leg: its charges off in two steps
/// and its van der Waals in eight — the complex leg's states, without the restraint. **The eight
/// van der Waals windows sample a segment that is identically zero** (the group is the whole
/// molecule, so nothing depends on `λ_v`); they are kept so that the two legs run one schedule and
/// the zero is seen on sampled configurations, at about 130 of the leg's 178 s. Four seeds of
/// 20 + 100 ps per window, every window recording every state. Printed: each seed's table, the
/// mean over seeds with the seeds' spread, the electrostatic part against ⟨ΔG_GB⟩ at both ends and
/// its exponential average, the nonpolar term, and the hydration free energy against FreeSolv.
/// Asserted: only that the van der Waals segment is exactly zero, as it must be with nothing to
/// decouple from. Written to `solvent_*.txt` and `solvent.log` in [`results_dir`], and the result
/// to `solvent_summary.txt` for the complex leg.
#[test]
#[ignore = "benzene alone in OBC II, 4 seeds × 11 windows × 120 ps, a few minutes with --release: run with --release -- --ignored --nocapture"]
fn the_solvent_leg_measured() {
    use pantometry_forcefield::solvation::{
        intrinsic_radius, nonpolar_energy, surface_area, PROBE_RADIUS,
    };
    let dir = results_dir();
    let t = std::time::Instant::now();
    let b = relaxed_binding();
    let ligand = b.ligand_force_field().clone().with_generalized_born();
    let start = b.ligand_alone_positions().to_vec();
    let all = vec![true; start.len()];
    let d = Decoupling::new(&ligand, &all).unwrap();
    let masses: Vec<f64> = b.elements()[b.pocket_len()..]
        .iter()
        .map(|e| e.mass())
        .collect();
    let grid = solvent_grid();
    // The complex leg's charges and van der Waals states: λ_e 1, 0.5, 0, then λ_v 0.9 … 0.
    let pick = |e: f64, v: f64| {
        grid.iter()
            .position(|l| l.electrostatics == e && l.van_der_waals == v)
            .unwrap()
    };
    let mut schedule = vec![pick(1.0, 1.0), pick(0.5, 1.0), pick(0.0, 1.0)];
    for v in [0.9, 0.75, 0.6, 0.45, 0.3, 0.2, 0.1, 0.0] {
        schedule.push(pick(0.0, v));
    }
    let kt = BOLTZMANN.to_si() * KELVIN;
    let mut results = Vec::new();
    for (i, &seed) in SOLVENT_SEEDS.iter().enumerate() {
        let name = format!("solvent{i}");
        let leg = Leg {
            name: &name,
            hamiltonian: &d,
            start: start.clone(),
            masses: masses.clone(),
            frozen: vec![false; start.len()],
            protocol: solvent_protocol(seed),
            grid: grid.clone(),
            record_all: true,
        };
        let records = run_leg(&leg, &dir, schedule.clone());
        // The last window with the charges on its way off; after it, only the van der Waals moves.
        let off = records
            .iter()
            .rposition(|r| r.lambda.van_der_waals == 1.0)
            .unwrap();
        for r in &records {
            assert!(
                r.gradients.iter().all(|g| g[2] == 0.0),
                "nothing to decouple"
            );
        }
        let a = analyse(
            &dir,
            &name,
            &records,
            kt,
            &[
                ("charges off", 0, off),
                ("van der Waals off", off, records.len() - 1),
            ],
        );
        // ΔG_GB at full charge on every sample of the two charge ends: E(λ_e = 1) − E(λ_e = 0),
        // which is the solvation term alone.
        let (one, zero) = (pick(1.0, 1.0), pick(0.0, 1.0));
        let solvation = |r: &Recorded| -> Vec<f64> {
            r.energies
                .iter()
                .map(|e| e[r.column(one)] - e[r.column(zero)])
                .collect()
        };
        let coupled = solvation(&records[0]);
        let decoupled = solvation(&records[2]);
        let mean = |x: &[f64]| x.iter().sum::<f64>() / x.len() as f64;
        // Zwanzig from the coupled end: ΔG_dec = −kT ln ⟨e^{βG}⟩₁.
        let zw: Vec<f64> = coupled.iter().map(|g| g / kt).collect();
        let zwanzig = -kt * log_mean_exp(&zw);
        log(
            &dir,
            &name,
            &format!(
                "{name}: ⟨ΔG_GB⟩ at λ_e = 1 {:.4}, at λ_e = 0 (full charges on its samples) {:.4} \
                 kcal/mol; linear response −(⟨G⟩₁ + ⟨G⟩₀)/2 = {:+.4}; Zwanzig from the coupled \
                 end {:+.4}; against the charges' BAR {:+.4} kcal/mol",
                kcal(mean(&coupled)),
                kcal(mean(&decoupled)),
                kcal(-(mean(&coupled) + mean(&decoupled)) / 2.0),
                kcal(zwanzig),
                a.segments[0].0
            ),
        );
        results.push(a);
    }
    let m = results.len() as f64;
    let mean = |f: &dyn Fn(&Analysis) -> f64| results.iter().map(f).sum::<f64>() / m;
    let spread = |f: &dyn Fn(&Analysis) -> f64| {
        let mu = mean(f);
        (results
            .iter()
            .map(|a| (f(a) - mu) * (f(a) - mu))
            .sum::<f64>()
            / (m - 1.0))
            .sqrt()
    };
    let bar = mean(&|a| a.bar.0);
    let bar_sd = spread(&|a| a.bar.0);
    let bar_se = mean(&|a| a.bar.1);
    let ti = mean(&|a| a.ti.0);
    // The nonpolar term at the start, Bondi radii + 1.4 Å, 1600 points.
    let radii: Vec<f64> = b.elements()[b.pocket_len()..]
        .iter()
        .map(|&e| intrinsic_radius(e))
        .collect();
    let area: f64 = surface_area(&radii, &start, PROBE_RADIUS, 1600)
        .iter()
        .sum();
    let np = kcal(nonpolar_energy(area));
    log(
        &dir,
        "solvent",
        &format!(
            "solvent leg over {} seeds: BAR {bar:+.4} kcal/mol, seeds' spread {bar_sd:.4} \
             (standard error of the mean {:.4}), each seed's own σ̂ {bar_se:.4}; TI {ti:+.4}. \
             Hydration free energy, polar only: {:+.3}; with the nonpolar 0.005 × {:.1} Å² = \
             {np:+.3}: {:+.3} kcal/mol, against FreeSolv's −0.90 ± 0.20. {:.0} s",
            results.len(),
            bar_sd / m.sqrt(),
            -bar,
            area / (ANGSTROM * ANGSTROM),
            -bar + np,
            t.elapsed().as_secs_f64()
        ),
    );
    std::fs::write(
        dir.join("solvent_summary.txt"),
        format!(
            "{bar:e} {:e} {np:e}\n",
            (bar_sd / m.sqrt()).max(bar_se / m.sqrt())
        ),
    )
    .unwrap();
}

/// The complex leg's protocol: 3c-1's — 0.5 fs, 1 ps⁻¹, 2 ps discarded, 10 ps sampled every 10 fs.
fn complex_protocol() -> Protocol {
    Protocol {
        time_step: 0.5 * FS,
        temperature: KELVIN,
        friction: 1e12,
        equilibration: 4000,
        stride: 20,
        samples: 1000,
        seed: 0x3C03,
    }
}

/// **Benzene decoupled from T4 lysozyme L99A under OBC II**, the complex leg, and the binding free
/// energy it closes with the solvent leg. The 3b setup — a 6 Å mobile zone in a 10 Å binding,
/// hydrogens relaxed, the zone minimised under GB (2000 steps, as 3b's GB run), the buffer frozen
/// and its descreening kept — and 3c-1's Boresch rule and protocol, on [`COMPLEX_SCHEDULE`] with a
/// window added at the candidate between any two whose overlap is below [`OVERLAP_THRESHOLD`].
/// Each window is written to `complex_NN.txt` in [`results_dir`] as it runs, and a window already
/// there is read rather than run, so a stopped run resumes; progress is in `complex.log`.
/// Printed: each window and interval, the three segments, the totals by BAR and TI, the halves and
/// EXP both ways, the restraint's release, the nonpolar estimate, and the cycle. Asserted: only
/// that the buffer kept its bits and the kept descreening held.
#[test]
#[ignore = "15 or more windows of 12 ps on 987 atoms under OBC II, about six hours with --release: run with --release -- --ignored --nocapture"]
fn the_complex_leg_measured() {
    use pantometry_forcefield::minimise::KCAL_PER_MOL_ANGSTROM;
    use pantometry_forcefield::solvation::{
        intrinsic_radius, nonpolar_energy, surface_area, PROBE_RADIUS,
    };
    let dir = results_dir();
    let t = std::time::Instant::now();
    let b = relaxed_binding();
    let mut c = Complex::new(&b, 6.0 * ANGSTROM, Solvent::GeneralizedBorn).expect("the zone");
    let p = c.minimise(2000, Binding::HYDROGEN_TOLERANCE);
    let start = c.positions().to_vec();
    log(
        &dir,
        "complex",
        &format!(
            "complex: {} atoms, {} mobile; minimised under OBC II {:?} in {} steps to {:.2e} \
             kcal/mol/Å ({:.0} s)",
            start.len(),
            c.mobile_count(),
            p.status,
            p.steps,
            p.max_force / KCAL_PER_MOL_ANGSTROM,
            t.elapsed().as_secs_f64()
        ),
    );
    let restraint = anchors(&b, &start);
    log(
        &dir,
        "complex",
        &format!(
            "Boresch: receptor {:?}, ligand {:?}; r₀ {:.3} Å, θ_A {:.1}°, θ_B {:.1}°, φ {:.1}° \
             {:.1}° {:.1}°; release to 1 M {:+.3} kcal/mol (extended {:+.3})",
            restraint.receptor,
            restraint.ligand,
            restraint.reference[0] / ANGSTROM,
            restraint.reference[1].to_degrees(),
            restraint.reference[2].to_degrees(),
            restraint.reference[3].to_degrees(),
            restraint.reference[4].to_degrees(),
            restraint.reference[5].to_degrees(),
            kcal(restraint.release_free_energy(KELVIN)),
            kcal(restraint.release_free_energy_extended(KELVIN))
        ),
    );
    let mask = ligand_mask(&b);
    let d = Decoupling::new(c.potential(), &mask)
        .expect("benzene decouples")
        .with_restraint(restraint);
    assert!(
        d.solvation().unwrap().reuses_frozen_terms(&start),
        "the buffer's descreening is kept at the start"
    );
    let leg = Leg {
        name: "complex",
        hamiltonian: &d,
        start: start.clone(),
        masses: b.elements().iter().map(|e| e.mass()).collect(),
        frozen: c.frozen(),
        protocol: complex_protocol(),
        grid: grid(),
        record_all: false,
    };
    let records = run_leg(&leg, &dir, COMPLEX_SCHEDULE.to_vec());
    let kt = BOLTZMANN.to_si() * KELVIN;
    let last =
        |pred: &dyn Fn(&Lambda) -> bool| records.iter().rposition(|r| pred(&r.lambda)).unwrap();
    let on = last(&|l| l.electrostatics == 1.0 && l.van_der_waals == 1.0);
    let off = last(&|l| l.van_der_waals == 1.0);
    let a = analyse(
        &dir,
        "complex",
        &records,
        kt,
        &[
            ("restraint on", 0, on),
            ("charges off", on, off),
            ("van der Waals off", off, records.len() - 1),
        ],
    );
    let release = kcal(restraint.release_free_energy(KELVIN));
    // The nonpolar estimate at the start: γ (A_complex − A_pocket − A_benzene), Bondi + 1.4 Å.
    let radii: Vec<f64> = b.elements().iter().map(|&e| intrinsic_radius(e)).collect();
    let n0 = b.pocket_len();
    let area =
        |r: &[f64], x: &[[f64; 3]]| -> f64 { surface_area(r, x, PROBE_RADIUS, 1600).iter().sum() };
    let buried =
        area(&radii, &start) - area(&radii[..n0], &start[..n0]) - area(&radii[n0..], &start[n0..]);
    let np = kcal(nonpolar_energy(buried));
    let wall: f64 = records.iter().map(|r| r.seconds).sum();
    log(
        &dir,
        "complex",
        &format!(
            "complex leg: BAR {:+.3} ± {:.3}, TI {:+.3} ± {:.3} kcal/mol over {} windows, {:.0} s \
             of windows ({:.2} h); ΔA on binding {:.1} Å², nonpolar {np:+.3} kcal/mol",
            a.bar.0,
            a.bar.1,
            a.ti.0,
            a.ti.1,
            records.len(),
            wall,
            wall / 3600.0,
            buried / (ANGSTROM * ANGSTROM)
        ),
    );
    match std::fs::read_to_string(dir.join("solvent_summary.txt")) {
        Ok(text) => {
            let v: Vec<f64> = text
                .split_whitespace()
                .map(|x| x.parse().unwrap())
                .collect();
            let (solvent, solvent_err) = (v[0], v[1]);
            let bind = binding_free_energy(solvent, a.bar.0, release);
            let err = (solvent_err * solvent_err + a.bar.1 * a.bar.1).sqrt();
            log(
                &dir,
                "complex",
                &format!(
                    "ΔG°_bind = ΔG_solvent − ΔG_complex − ΔG°_release = {solvent:+.3} − ({:+.3}) − \
                     ({release:+.3}) = {bind:+.3} ± {err:.3} kcal/mol (polar GB only); with the \
                     nonpolar estimate {:+.3}; experiment −5.19 ± 0.16",
                    a.bar.0,
                    bind + np
                ),
            );
        }
        Err(e) => log(
            &dir,
            "complex",
            &format!(
                "no solvent leg in {} ({e}): run the_solvent_leg_measured",
                dir.display()
            ),
        ),
    }
}
