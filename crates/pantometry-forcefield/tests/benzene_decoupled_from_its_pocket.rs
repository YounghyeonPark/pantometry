//! **The decoupling Hamiltonian on T4 lysozyme L99A and benzene (PDB 181L)**, step 3c-1: what has
//! to hold exactly when the group is a real ligand in a real pocket with QEq charges, and — ignored,
//! run once in release — benzene decoupled from the pocket in vacuum with a Boresch restraint, the
//! engine exercised end to end.
//!
//! The default tests run on 181L's pocket at 3 Å (benzene and four residues, 82 atoms):
//! `∂U/∂λ_e` is the coupled cross Coulomb energy to the bit at every state, and is
//! `Binding::cross_terms_at`'s; the soft core at λ_v = 1 is the binding's cross van der Waals to
//! the bit; full coupling is the force field to its rounding; the decoupled ligand feels nothing,
//! exactly; the forces at an intermediate state are the energy's gradient; and what cannot be
//! decoupled is refused by name.
//!
//! The demonstration: `cargo test -p pantometry-forcefield --release --test
//! benzene_decoupled_from_its_pocket -- --ignored --nocapture`. **A vacuum leg has no solvent
//! counterpart here and its number cannot be compared with experiment.** It is reported, not
//! asserted.

mod protein;

use pantometry_forcefield::energy::coulomb;
use pantometry_forcefield::minimise::KCAL_PER_MOL_ANGSTROM;
use pantometry_forcefield::uff::KCAL_PER_MOL;
use pantometry_forcefield::{
    Alchemical, AlchemyError, Binding, Boresch, Complex, Decoupling, Lambda, Protocol, Quadrature,
    Solvent, System, Windows,
};
use pantometry_units::BOLTZMANN;
use protein::*;
use std::sync::OnceLock;
use std::time::Instant;

const ANGSTROM: f64 = 1e-10;
const FS: f64 = 1e-15;
const KELVIN: f64 = 300.0;
const EPS: f64 = f64::EPSILON;

fn kcal(j: f64) -> f64 {
    j / KCAL_PER_MOL
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

/// The Boresch force constants of the demonstration and the tests, in this crate's `½ K (ξ − ξ₀)²`:
/// 20 kcal mol⁻¹ Å⁻² and 20 kcal mol⁻¹ rad⁻², which is Mobley, Chodera and Dill's K₀ = 10 in their
/// `K₀ (ξ − ξ₀)²` (see [`pantometry_forcefield::boresch`]).
fn force_constants() -> [f64; 6] {
    let r = 20.0 * KCAL_PER_MOL / (ANGSTROM * ANGSTROM);
    let a = 20.0 * KCAL_PER_MOL;
    [r, a, a, a, a, a]
}

/// The anchor rule applied to a binding at `at`: receptor candidates `(N, CA, C)` and `(C, CA,
/// N)` of every pocket residue, ligand candidates every path `A–B–C` of three bonded ligand heavy
/// atoms, `r₀` in [3, 8] Å; [`Boresch::choose`] keeps the pair whose hinges are furthest from
/// straight.
fn anchors(b: &Binding, at: &[[f64; 3]]) -> Boresch {
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
    let heavy = |k: usize| b.elements()[k] != pantometry_forcefield::Element::H;
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
    Boresch::choose(
        at,
        &receptor,
        &ligand,
        (3.0 * ANGSTROM, 8.0 * ANGSTROM),
        force_constants(),
    )
    .expect("an anchor pair within 3–8 Å")
}

/// **`∂U/∂λ_e` is the coupled cross Coulomb energy, to the bit, at every state**, with QEq's
/// charges on 181L's 3 Å pocket and benzene: the same bits at seven states, equal to the test's
/// own sum of `coulomb(q_i, q_j, r)` over every ligand–protein pair and to
/// `Binding::cross_terms_at`'s, and **the soft core at `λ_v = 1` is the binding's cross van der
/// Waals to the bit** (the coupling energy with the charges off). Linear in `λ_e`: the coupling
/// energy at `λ_e` less that at 0 is `λ_e` times it, to the rounding of two sums over the pairs.
#[test]
fn the_electrostatic_gradient_is_the_coupled_cross_coulomb_exactly() {
    let b = small_binding();
    let ff = b.force_field();
    let d = Decoupling::new(ff, &ligand_mask(b)).unwrap();
    let at = b.positions();
    let n0 = b.pocket_len();
    let q = b.charges();
    let (mut own, mut absolute, mut pairs) = (0.0, 0.0, 0usize);
    for p in ff.pairs() {
        let [i, j] = p.atoms;
        if (i < n0) != (j < n0) {
            let r = {
                let v = [
                    at[i][0] - at[j][0],
                    at[i][1] - at[j][1],
                    at[i][2] - at[j][2],
                ];
                (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
            };
            let e = coulomb(q[i], q[j], r);
            own += e;
            absolute += e.abs() + p.energy(r).abs();
            pairs += 1;
        }
    }
    let (vdw, elec) = b.cross_terms_at(at);
    assert!(elec.abs() > 1e-3 * KCAL_PER_MOL, "a real Coulomb term");
    for l in [
        Lambda::new(0.0, 1.0, 1.0),
        Lambda::new(0.0, 0.0, 1.0),
        Lambda::new(0.0, 0.3, 1.0),
        Lambda::new(1.0, 0.7, 1.0),
        Lambda::new(0.0, 0.0, 0.4),
        Lambda::new(0.0, 0.0, 0.0),
        Lambda::new(0.5, 0.5, 0.05),
    ] {
        let g = d.coupling(at, l).gradient[1];
        assert_eq!(g.to_bits(), own.to_bits(), "{l:?}");
        assert_eq!(g.to_bits(), elec.to_bits(), "{l:?}");
    }
    let off = d.coupling(at, Lambda::new(0.0, 0.0, 1.0)).energy;
    assert_eq!(off.to_bits(), vdw.to_bits());
    // Two sums of `pairs` terms each, interleaved pair by pair: each rounds by at most
    // `pairs ε Σ|t|`.
    for le in [0.25, 0.5, 0.9] {
        let on = d.coupling(at, Lambda::new(0.0, le, 1.0)).energy;
        let want = le * elec;
        assert!(
            ((on - off) - want).abs() <= 2.0 * pairs as f64 * EPS * absolute,
            "λ_e {le}: {:e}",
            (on - off) - want
        );
    }
    println!(
        "cross terms at 3 Å: van der Waals {:.4}, Coulomb {:.4} kcal/mol; {} cross pairs",
        kcal(vdw),
        kcal(elec),
        d.cross_pairs().len()
    );
}

/// **Fully coupled, the decoupling is the force field; decoupled, the ligand feels nothing.** At
/// `Lambda::COUPLED` the energy and every force are the force field's, to a traced allowance:
/// the two sum the same terms in another order, so they differ by at most `n ε Σ|t|` with `n` the
/// number of terms and `Σ|t|` their absolute sum, bounded here by `n ε` times the sum of each
/// energy field's magnitude plus the cross terms'. At `λ_e = λ_v = 0` with no restraint the
/// coupling energy is exactly zero and every force is the rest's to the bit.
#[test]
fn coupled_it_is_the_force_field_and_decoupled_it_is_nothing() {
    let b = small_binding();
    let ff = b.force_field();
    let d = Decoupling::new(ff, &ligand_mask(b)).unwrap();
    let at = b.positions();
    let n = at.len();
    let whole = ff.evaluate(at);
    let mut forces = vec![[0.0; 3]; n];
    let e = d.energy_and_forces(at, Lambda::COUPLED, &mut forces);
    let w = whole.energy;
    let (vdw, elec) = b.cross_terms_at(at);
    let magnitude = w.bond.abs()
        + w.angle.abs()
        + w.torsion.abs()
        + w.inversion.abs()
        + w.van_der_waals.abs()
        + w.electrostatic.abs()
        + 2.0 * (vdw.abs() + elec.abs());
    let terms = (ff.stretches().len()
        + ff.bends().len()
        + ff.torsions().len()
        + ff.inversions().len()
        + 2 * ff.pairs().len()) as f64;
    let allowance = terms * EPS * magnitude;
    println!(
        "coupled: decoupling {:.9} against force field {:.9} kcal/mol, gap {:.1e} against an \
         allowance of {:.1e}",
        kcal(e),
        kcal(w.total),
        kcal((e - w.total).abs()),
        kcal(allowance)
    );
    assert!((e - w.total).abs() <= allowance);
    let largest = whole
        .forces
        .iter()
        .flatten()
        .fold(0.0f64, |m, f| m.max(f.abs()));
    for (f, g) in forces.iter().zip(&whole.forces) {
        for c in 0..3 {
            assert!((f[c] - g[c]).abs() <= terms * EPS * 10.0 * largest);
        }
    }

    let none = Lambda::new(0.0, 0.0, 0.0);
    assert_eq!(d.coupling(at, none).energy, 0.0);
    let mut decoupled = vec![[0.0; 3]; n];
    let e0 = d.energy_and_forces(at, none, &mut decoupled);
    let rest = d.rest().unwrap().evaluate(at);
    assert_eq!(e0.to_bits(), rest.energy.total.to_bits());
    for (f, g) in decoupled.iter().zip(&rest.forces) {
        for c in 0..3 {
            assert_eq!(f[c].to_bits(), g[c].to_bits());
        }
    }
}

/// **At an intermediate state the force is the energy's gradient**, Boresch restraint included:
/// λ = (0.6, 0.4, 0.7) on the 3 Å complex with benzene moved 0.3 Å off its crystal pose, every
/// coordinate of every ligand atom and of the three receptor anchors, against central differences
/// at 1e-5 Å and its half, each within 1e-6 of the largest force component on those atoms. A
/// difference's truncation there is `h² U‴/6`, about 1e-30 m² × 10¹³ J m⁻³, and its rounding
/// `ε |U|/h` about 10⁻¹⁹ N, both far below 10⁻¹⁵ N. And `∂U/∂λ` in each component against a
/// central difference in that λ.
#[test]
fn at_an_intermediate_state_the_force_is_the_gradient() {
    let b = small_binding();
    let mut at = b.positions().to_vec();
    for p in &mut at[b.pocket_len()..] {
        p[0] += 0.3 * ANGSTROM;
    }
    let restraint = anchors(b, b.positions());
    let d = Decoupling::new(b.force_field(), &ligand_mask(b))
        .unwrap()
        .with_restraint(restraint);
    let l = Lambda::new(0.6, 0.4, 0.7);
    let mut forces = vec![[0.0; 3]; at.len()];
    d.energy_and_forces(&at, l, &mut forces);
    let mut atoms: Vec<usize> = b.ligand_range().collect();
    atoms.extend(restraint.receptor);
    let largest = atoms
        .iter()
        .flat_map(|&k| forces[k])
        .fold(0.0f64, |m, f| m.max(f.abs()));
    let energy = |p: &[[f64; 3]]| {
        let mut scratch = vec![[0.0; 3]; p.len()];
        d.energy_and_forces(p, l, &mut scratch)
    };
    for &k in &atoms {
        for c in 0..3 {
            let fd = |h: f64| {
                let (mut plus, mut minus) = (at.clone(), at.clone());
                plus[k][c] += h;
                minus[k][c] -= h;
                (energy(&plus) - energy(&minus)) / (2.0 * h)
            };
            for h in [1e-5 * ANGSTROM, 0.5e-5 * ANGSTROM] {
                let err = (fd(h) + forces[k][c]).abs();
                assert!(
                    err <= 1e-6 * largest,
                    "atom {k} {c}: {err:e} of {largest:e}"
                );
            }
        }
    }
    let g = d.coupling(&at, l).gradient;
    for c in 0..3 {
        let shifted = |x: f64| {
            let mut m = l.components();
            m[c] = x;
            d.coupling(&at, Lambda::new(m[0], m[1], m[2])).energy
        };
        let x = l.components()[c];
        let h = 1e-5;
        let fd = (shifted(x + h) - shifted(x - h)) / (2.0 * h);
        assert!(
            (fd - g[c]).abs() <= 1e-6 * g[c].abs().max(1e-3 * KCAL_PER_MOL),
            "∂U/∂λ[{c}]: {fd:e} against {:e}",
            g[c]
        );
    }
    println!(
        "anchors: receptor {:?}, ligand {:?}; r₀ {:.3} Å, θ {:.1}° {:.1}°",
        restraint.receptor,
        restraint.ligand,
        restraint.reference[0] / ANGSTROM,
        restraint.reference[1].to_degrees(),
        restraint.reference[2].to_degrees()
    );
}

/// **What cannot be decoupled is refused by name**: a solvated force field, a group that a bond
/// crosses (one benzene carbon alone), and a group that is empty or everything.
#[test]
fn what_cannot_be_decoupled_is_refused() {
    let b = small_binding();
    let ff = b.force_field();
    let mask = ligand_mask(b);
    assert_eq!(
        Decoupling::new(&ff.clone().with_generalized_born(), &mask),
        Err(AlchemyError::Solvated)
    );
    let mut one = vec![false; mask.len()];
    let carbon = b
        .ligand_range()
        .find(|&k| b.elements()[k] == pantometry_forcefield::Element::C)
        .unwrap();
    one[carbon] = true;
    match Decoupling::new(ff, &one) {
        Err(AlchemyError::BondedAcross { atoms }) => assert!(atoms.contains(&carbon)),
        other => panic!("{other:?}"),
    }
    assert_eq!(
        Decoupling::new(ff, &vec![false; mask.len()]),
        Err(AlchemyError::NothingToDecouple)
    );
    assert_eq!(
        Decoupling::new(ff, &vec![true; mask.len()]),
        Err(AlchemyError::NothingToDecouple)
    );
}

// ---------------------------------------------------------------------------------------------
// The demonstration
// ---------------------------------------------------------------------------------------------

/// The complex leg's schedule: the restraint on in six states, the charges off in four more, the
/// van der Waals off in twelve, denser at both ends.
fn schedule() -> Vec<Lambda> {
    let mut s = Vec::new();
    for r in [0.0, 0.1, 0.25, 0.5, 0.75, 1.0] {
        s.push(Lambda::new(r, 1.0, 1.0));
    }
    for e in [0.75, 0.5, 0.25, 0.0] {
        s.push(Lambda::new(1.0, e, 1.0));
    }
    for v in [0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.0] {
        s.push(Lambda::new(1.0, 0.0, v));
    }
    s
}

/// **Benzene decoupled from T4 lysozyme L99A in vacuum**, with the 3b setup — a 6 Å mobile zone
/// in a 10 Å binding, hydrogens relaxed and the zone minimised, the buffer frozen — and a Boresch
/// restraint chosen by the rule. 22 windows at 0.5 fs in a 1 ps⁻¹ bath, 2 ps discarded and
/// 10 ps sampled every 10 fs in each. Printed: each window's TI term and each interval's BAR
/// with its overlap, the three segments and their sum, the restraint's analytic release, and the
/// cost. Asserted: only that the buffer kept its bits and every window finished.
#[test]
#[ignore = "22 windows of 12 ps on 987 atoms, about fifteen minutes with --release: run with --release -- --ignored --nocapture"]
fn benzene_decoupled_in_vacuum_measured() {
    let s = the_system();
    let t = Instant::now();
    let b = Binding::new(s, 10.0 * ANGSTROM)
        .expect("the binding")
        .relaxing_hydrogens(20_000, Binding::HYDROGEN_TOLERANCE);
    let mut c = Complex::new(&b, 6.0 * ANGSTROM, Solvent::Vacuum).expect("the zone");
    let p = c.minimise(20_000, Binding::HYDROGEN_TOLERANCE);
    println!(
        "binding and zone: {} atoms, {} mobile; minimised {:?} in {} steps to {:.2e} kcal/mol/Å \
         ({:.0} s)",
        c.positions().len(),
        c.mobile_count(),
        p.status,
        p.steps,
        p.max_force / KCAL_PER_MOL_ANGSTROM,
        t.elapsed().as_secs_f64()
    );
    let start = c.positions().to_vec();
    let restraint = anchors(&b, &start);
    let name = |k: usize| {
        let i = b.system_atoms()[k];
        s.atom_name(i).to_string()
    };
    println!(
        "Boresch: receptor {} {} {} (complex {:?}), ligand {} {} {}; r₀ {:.3} Å, θ_A {:.1}°, \
         θ_B {:.1}°, φ {:.1}° {:.1}° {:.1}°; hinges {:?}°",
        name(restraint.receptor[0]),
        name(restraint.receptor[1]),
        name(restraint.receptor[2]),
        restraint.receptor,
        name(restraint.ligand[0]),
        name(restraint.ligand[1]),
        name(restraint.ligand[2]),
        restraint.reference[0] / ANGSTROM,
        restraint.reference[1].to_degrees(),
        restraint.reference[2].to_degrees(),
        restraint.reference[3].to_degrees(),
        restraint.reference[4].to_degrees(),
        restraint.reference[5].to_degrees(),
        restraint
            .hinge_angles(&start)
            .map(|x| (x.to_degrees() * 10.0).round() / 10.0)
    );
    let d = Decoupling::new(c.potential(), &ligand_mask(&b))
        .expect("benzene decouples")
        .with_restraint(restraint);
    let masses: Vec<f64> = b.elements().iter().map(|e| e.mass()).collect();
    let protocol = Protocol {
        time_step: 0.5 * FS,
        temperature: KELVIN,
        friction: 1e12,
        equilibration: 4000,
        stride: 20,
        samples: 1000,
        seed: 0x3C01,
    };
    let schedule = schedule();
    let mut w = Windows::new(&d, schedule.clone(), &start, masses, c.frozen(), protocol);
    let mut seconds = Vec::new();
    for k in 0..schedule.len() {
        let t = Instant::now();
        w.window_mut(k).advance(&d, u64::MAX);
        seconds.push(t.elapsed().as_secs_f64());
        assert!(w.windows()[k].is_complete());
    }
    for win in w.windows() {
        for (k, &m) in c.mobile().iter().enumerate() {
            if !m {
                assert_eq!(win.positions()[k], start[k], "frozen atom {k}");
            }
        }
    }
    let kt = BOLTZMANN.to_si() * KELVIN;
    let terms = w.ti_terms(Quadrature::Trapezoid).unwrap();
    let bars = w.bennett();
    println!(
        "| window | λ_r | λ_e | λ_v | TI term (kcal/mol) | τ (samples) | s | → BAR Δ (kcal/mol) | overlap |"
    );
    println!("| --- | --- | --- | --- | --- | --- | --- | --- | --- |");
    for (k, l) in schedule.iter().enumerate() {
        let bar = bars.get(k).map_or(String::from("|"), |x| {
            format!(
                "{:+.3} ± {:.3} | {:.3}",
                kcal(kt * x.delta),
                kcal(kt * x.error()),
                x.overlap
            )
        });
        println!(
            "| {k} | {} | {} | {} | {:+.3} ± {:.3} | {:.1} | {:.0} | {bar} |",
            l.restraint,
            l.electrostatics,
            l.van_der_waals,
            kcal(terms[k].mean),
            kcal(terms[k].error),
            terms[k].tau,
            seconds[k]
        );
    }
    // The three segments by BAR: intervals 0–4 the restraint, 5–8 the charges, 9–20 the van der
    // Waals.
    let segment = |range: std::ops::Range<usize>| {
        let v: f64 = bars[range.clone()].iter().map(|x| x.delta).sum();
        let e: f64 = bars[range].iter().map(|x| x.variance).sum::<f64>().sqrt();
        (kcal(kt * v), kcal(kt * e))
    };
    let (r, re) = segment(0..5);
    let (q, qe) = segment(5..9);
    let (v, ve) = segment(9..21);
    let total = w.bennett_total();
    let ti = w.thermodynamic_integration(Quadrature::Trapezoid).unwrap();
    let release = restraint.release_free_energy(KELVIN);
    println!(
        "BAR by segment (each segment's intervals' variances summed): restraint on {r:+.3} ± \
         {re:.3}, charges off {q:+.3} ± {qe:.3}, van der Waals off {v:+.3} ± {ve:.3} kcal/mol"
    );
    println!(
        "complex leg (restrained decoupling, coupled → decoupled): BAR {:+.3} ± {:.3}, TI \
         (trapezoid) {:+.3} ± {:.3} kcal/mol; the restraint's release to 1 M, analytic: {:+.3} \
         kcal/mol",
        kcal(total.value),
        kcal(total.error),
        kcal(ti.value),
        kcal(ti.error),
        kcal(release)
    );
    let total_s: f64 = seconds.iter().sum();
    let steps = protocol.steps_per_window() as f64;
    println!(
        "cost: {:.0} s for {} windows, {:.1} s per window of {} steps, {:.2} ms per step",
        total_s,
        schedule.len(),
        total_s / schedule.len() as f64,
        protocol.steps_per_window(),
        total_s / (schedule.len() as f64 * steps) * 1e3
    );
}
