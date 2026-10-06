//! **Generalized Born's evaluation is its direct sum, to the bit.** Step 3c-2 made OBC II faster
//! without approximating it: [`GeneralizedBorn::accumulate`] computes each pair integral once and
//! keeps its `dI/dr` for the chain rule, reuses the descreening between two frozen atoms that
//! [`GeneralizedBorn::with_frozen`] computed once, and skips `e^(−κf)` without salt, where it is
//! one. None of that may change a bit. Every test here holds the evaluation to
//! [`GeneralizedBorn::accumulate_reference`] and [`GeneralizedBorn::born_radii_reference`] — the
//! code that was there before, kept — **bit for bit**: the energy, every Born radius and every
//! component of every force, with no tolerance, on aspirin, on a cluster that reaches every branch
//! of the pair integral (two of its atoms coincident), and on T4 lysozyme L99A's pocket around
//! benzene (PDB 181L), each under several frozen masks, rescalings, dielectrics and salt.
//!
//! The kept terms are only valid while the frozen atoms are where they were, so the tests also
//! move them — by 0.1 Å and by one unit in the last place — and require the evaluation to notice,
//! to compute every term, and to be the reference's still.
//!
//! In release and ignored, the same on the 987-atom complex of step 3b (a 6 Å zone in a 10 Å
//! binding) and the cost of a step, which is what the step was for:
//! `cargo test -p pantometry-forcefield --release --test the_generalized_born_is_its_direct_sum
//! -- --ignored --nocapture`.

mod protein;

use pantometry_forcefield::solvation::{self, GeneralizedBorn, Rescaling};
use pantometry_forcefield::{
    Bath, Binding, Complex, Component, Element, ForceField, MolecularDynamics, Potential, Record,
    Solvent, System,
};
use protein::*;
use std::sync::OnceLock;
use std::time::Instant;

const AIN: &str = include_str!("../components/AIN.cif");
const ANGSTROM: f64 = 1e-10;
const FS: f64 = 1e-15;

fn bits(v: &[[f64; 3]]) -> Vec<u64> {
    v.iter().flatten().map(|x| x.to_bits()).collect()
}

fn scalar_bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

/// A force array that is not zero, so that "adds to" is checked as well as "computes": the same
/// starting values go into both evaluations.
fn start(n: usize) -> Vec<[f64; 3]> {
    (0..n)
        .map(|i| [0, 1, 2].map(|k| 1e-11 * (0.37 * (3 * i + k) as f64).sin()))
        .collect()
}

/// **`gb` at `at` is its reference, bit for bit**: the Born radii, the energy alone, and the
/// energy with every force added to [`start`]. Returns whether the kept frozen terms were used.
fn is_the_reference(name: &str, gb: &GeneralizedBorn, q: &[f64], at: &[[f64; 3]]) -> bool {
    let n = at.len();
    assert_eq!(
        scalar_bits(&gb.born_radii(at)),
        scalar_bits(&gb.born_radii_reference(at)),
        "{name}: the Born radii"
    );
    let (mut fast, mut slow) = (start(n), start(n));
    let e = gb.accumulate(q, at, &mut fast);
    let e_ref = gb.accumulate_reference(q, at, &mut slow);
    assert!(e.is_finite(), "{name}: {e}");
    assert_eq!(e.to_bits(), e_ref.to_bits(), "{name}: the energy");
    assert_eq!(
        gb.energy(q, at).to_bits(),
        e_ref.to_bits(),
        "{name}: the energy alone"
    );
    for k in 0..n {
        assert_eq!(
            bits(&[fast[k]]),
            bits(&[slow[k]]),
            "{name}: atom {k}'s force"
        );
    }
    // The forces are not what they started as: the check compared something.
    assert_ne!(bits(&fast), bits(&start(n)), "{name}: no force was added");
    gb.reuses_frozen_terms(at)
}

/// The masks every system is checked under, by name: nothing kept, nothing frozen, everything
/// frozen, every third atom, the first half, and every atom but the last.
fn masks(n: usize) -> Vec<(&'static str, Option<Vec<bool>>)> {
    vec![
        ("no kept terms", None),
        ("none frozen", Some(vec![false; n])),
        ("all frozen", Some(vec![true; n])),
        (
            "every third frozen",
            Some((0..n).map(|k| k % 3 == 0).collect()),
        ),
        (
            "first half frozen",
            Some((0..n).map(|k| k < n / 2).collect()),
        ),
        (
            "all but the last frozen",
            Some((0..n).map(|k| k + 1 < n).collect()),
        ),
    ]
}

/// Every mask of [`masks`] on `gb`, at `at` and at `at` with the mobile atoms moved, each
/// required to be the reference; and with a mask, at `at` with one frozen atom moved by 0.1 Å
/// and by one unit in the last place, which must not reuse the kept terms and must still be the
/// reference.
fn under_every_mask(name: &str, gb: &GeneralizedBorn, q: &[f64], at: &[[f64; 3]]) {
    let n = at.len();
    for (label, mask) in masks(n) {
        let model = match &mask {
            None => gb.clone(),
            Some(m) => gb.clone().with_frozen(m, at),
        };
        let frozen = mask.clone().unwrap_or_else(|| vec![false; n]);
        let kept = frozen.iter().any(|&f| f);
        assert_eq!(model.frozen_count(), frozen.iter().filter(|&&f| f).count());
        let used = is_the_reference(&format!("{name}, {label}"), &model, q, at);
        assert_eq!(
            used,
            mask.is_some(),
            "{name}, {label}: the kept terms at their own positions"
        );
        // The mobile atoms moved: the kept terms still hold, and are still the reference.
        let mut moved = at.to_vec();
        for (k, p) in moved.iter_mut().enumerate() {
            if !frozen[k] {
                p[0] += 0.03 * ANGSTROM * (k as f64 + 1.0).sin();
                p[2] -= 0.02 * ANGSTROM;
            }
        }
        let used = is_the_reference(&format!("{name}, {label}, mobile moved"), &model, q, &moved);
        assert_eq!(used, mask.is_some(), "{name}, {label}: mobile atoms moved");
        if !kept {
            continue;
        }
        // A frozen atom moved — by 0.1 Å, and by one ulp of one coordinate — is noticed: the kept
        // terms are not used, and the result is the reference's.
        let k = frozen.iter().position(|&f| f).expect("a frozen atom");
        for (how, by) in [("0.1 Å", 0.1 * ANGSTROM), ("one ulp", 0.0)] {
            let mut shifted = moved.clone();
            shifted[k][1] = if by == 0.0 {
                f64::from_bits(shifted[k][1].to_bits() + 1)
            } else {
                shifted[k][1] + by
            };
            let label = format!("{name}, {label}, frozen atom {k} moved by {how}");
            assert!(
                !is_the_reference(&label, &model, q, &shifted),
                "{label}: kept terms reused"
            );
        }
    }
}

/// Aspirin with deterministic neutral charges up to ±0.3 e, moved off the dictionary geometry by
/// up to 0.05 Å per coordinate: the closed-form tests' molecule.
fn aspirin() -> (Vec<Element>, Vec<f64>, Vec<[f64; 3]>) {
    let c = Component::from_ccd(AIN).expect("AIN parses");
    let n = c.atoms().len();
    let raw: Vec<f64> = (0..n).map(|i| 0.3 * (1.7 * i as f64).cos()).collect();
    let mean = raw.iter().sum::<f64>() / n as f64;
    let at = c
        .atoms()
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let mut p = a.at;
            for (k, x) in p.iter_mut().enumerate() {
                *x += 0.05 * ANGSTROM * (1.3 * (3 * i + k) as f64 + 0.7).sin();
            }
            p
        })
        .collect();
    (
        c.atoms().iter().map(|a| a.element).collect(),
        raw.iter().map(|q| q - mean).collect(),
        at,
    )
}

/// **Aspirin, under every mask, in OBC II, OBC I and HCT, at ε = 80 and 4, and in 0.15 M salt.**
#[test]
fn aspirin_is_its_direct_sum() {
    let (elements, q, at) = aspirin();
    let base = GeneralizedBorn::new(&elements);
    let models = [
        ("OBC II", base.clone()),
        ("OBC I", base.clone().with_rescaling(Rescaling::OBC_I)),
        ("HCT", base.clone().with_rescaling(Rescaling::Hct)),
        ("OBC II at ε = 4", base.clone().with_solvent_dielectric(4.0)),
        (
            "OBC II in 0.15 M salt",
            base.clone().with_kappa(solvation::debye_kappa(0.15)),
        ),
    ];
    for (name, gb) in &models {
        under_every_mask(&format!("aspirin, {name}"), gb, &q, &at);
    }
}

/// **A cluster that reaches every branch of the pair integral**: the seven atoms with which the
/// closed-form tests reach every branch — a large sphere with two small ones inside it, one small
/// atom inside another, pairs overlapping both ways — and an eighth **at the same point as the fourth**, so that a pair at
/// `r = 0` is in the sums: its integral takes the `r = 0` branch and the chain rule skips it,
/// which is where a kept slope could be read against the wrong pair.
#[test]
fn a_cluster_with_every_branch_and_a_coincident_pair_is_its_direct_sum() {
    let atoms: [(f64, f64, [f64; 3]); 8] = [
        (3.0, 1.0, [0.0, 0.0, 0.0]),
        (1.2, 0.85, [0.4, 0.3, -0.2]),
        (1.5, 0.85, [-0.6, 0.5, 0.7]),
        (1.7, 0.72, [3.1, 0.4, 0.2]),
        (1.2, 0.85, [3.4, 1.0, 0.9]),
        (1.55, 0.79, [5.6, -0.8, 0.3]),
        (1.8, 0.96, [1.6, 3.9, -1.1]),
        (1.5, 0.85, [3.1, 0.4, 0.2]),
    ];
    let gb = GeneralizedBorn::from_radii(
        atoms.iter().map(|a| a.0 * ANGSTROM).collect(),
        atoms.iter().map(|a| a.1).collect(),
    );
    let q: Vec<f64> = (0..atoms.len())
        .map(|i| 0.6 * (1.1 * i as f64 + 0.3).sin())
        .collect();
    let at: Vec<[f64; 3]> = atoms.iter().map(|a| a.2.map(|x| x * ANGSTROM)).collect();
    assert_eq!(bits(&[at[3]]), bits(&[at[7]]), "two atoms at one point");
    for (name, gb) in [
        ("OBC II", gb.clone()),
        ("HCT", gb.clone().with_rescaling(Rescaling::Hct)),
        ("salt", gb.clone().with_kappa(solvation::debye_kappa(0.5))),
    ] {
        under_every_mask(&format!("cluster, {name}"), &gb, &q, &at);
    }
    // The coincident pair split across the mask, each way round: one of the two frozen alone, so
    // that its `r = 0` slope is in a mobile atom's row and in a frozen atom's.
    for k in [3, 7] {
        let mask: Vec<bool> = (0..atoms.len()).map(|j| j == k).collect();
        let model = gb.clone().with_frozen(&mask, &at);
        assert!(is_the_reference(
            &format!("cluster, atom {k} frozen alone"),
            &model,
            &q,
            &at
        ));
    }
}

/// 181L as a system, once.
fn the_system() -> &'static System {
    static S: OnceLock<System> = OnceLock::new();
    S.get_or_init(system)
}

/// 181L's pocket at 3 Å: benzene and four side chains, 82 atoms, QEq charges.
fn small_binding() -> &'static Binding {
    static B: OnceLock<Binding> = OnceLock::new();
    B.get_or_init(|| Binding::new(the_system(), 3.0 * ANGSTROM).expect("181L at 3 Å"))
}

/// The backbone's atom names: what stays frozen in the small complex.
const BACKBONE: [&str; 6] = ["N", "CA", "C", "O", "H", "HA"];

/// The small complex's mobile mask of `tests/the_complex_in_motion.rs`: benzene and the four side
/// chains free, every backbone atom frozen.
fn small_mask(b: &Binding) -> Vec<bool> {
    let s = the_system();
    (0..b.positions().len())
        .map(|k| k >= b.pocket_len() || !BACKBONE.contains(&s.atom_name(b.system_atoms()[k])))
        .collect()
}

/// **T4 lysozyme's pocket around benzene, with its QEq charges**, under every mask and under the
/// small complex's own — the backbones frozen — in OBC II and in salt.
#[test]
fn the_pocket_is_its_direct_sum() {
    let b = small_binding();
    let (q, at) = (b.charges(), b.positions());
    let gb = GeneralizedBorn::new(b.elements());
    under_every_mask("181L at 3 Å", &gb, q, at);
    let frozen: Vec<bool> = small_mask(b).iter().map(|m| !m).collect();
    assert!(frozen.iter().any(|&f| f) && frozen.iter().any(|&f| !f));
    for (name, gb) in [
        ("OBC II", gb.clone()),
        ("salt", gb.clone().with_kappa(solvation::debye_kappa(0.15))),
    ] {
        let model = gb.with_frozen(&frozen, at);
        assert!(is_the_reference(
            &format!("181L at 3 Å, backbones frozen, {name}"),
            &model,
            q,
            at
        ));
    }
}

/// **A complex in generalized Born keeps its frozen atoms' terms, and uses them all the way
/// through a run.** [`Complex::with_mobile`] gives its potential's model the buffer as frozen, at
/// the binding's positions; after minimising the mobile zone and 40 steps of dynamics the buffer
/// has not moved, so the kept terms still hold — the fast path is the one taken — and the
/// evaluation there is the reference's, bit for bit. In vacuum there is no model.
#[test]
fn a_complex_keeps_its_buffers_terms_through_a_run() {
    let b = small_binding();
    let mut c = Complex::with_mobile(b, small_mask(b), Solvent::GeneralizedBorn).unwrap();
    let gb = c.potential().solvation().expect("a solvent").clone();
    assert_eq!(gb.frozen_count(), b.positions().len() - c.mobile_count());
    assert!(
        gb.reuses_frozen_terms(c.positions()),
        "at the binding's positions"
    );
    c.minimise(200, Binding::HYDROGEN_TOLERANCE);
    assert!(gb.reuses_frozen_terms(c.positions()), "after minimising");
    let mut md = c.dynamics(
        Bath::Langevin {
            temperature: 300.0,
            friction: 1e12,
            seed: 0x3C2,
        },
        300.0,
        0x3C3,
    );
    let mut record = Record::new(&c, 10);
    c.run(&mut md, 0.5 * FS, 40, &mut record);
    assert!(is_the_reference(
        "181L at 3 Å after 40 steps",
        &gb,
        c.potential().charges(),
        c.positions()
    ));
    let vacuum = Complex::with_mobile(b, small_mask(b), Solvent::Vacuum).unwrap();
    assert!(vacuum.potential().solvation().is_none());
    the_run_is_the_direct_run("181L at 3 Å", b, &c, 60);
}

/// A complex's potential with its generalized Born taken by the direct evaluation: the vacuum
/// potential's terms, then [`GeneralizedBorn::accumulate_reference`] added to its forces and its
/// total — the order [`ForceField::evaluate`] adds the solvation in.
struct Direct<'a> {
    vacuum: &'a ForceField,
    gb: &'a GeneralizedBorn,
}

impl Potential for Direct<'_> {
    fn energy_and_forces(&self, at: &[[f64; 3]], forces: &mut [[f64; 3]]) -> f64 {
        let ev = self.vacuum.evaluate(at);
        forces.copy_from_slice(&ev.forces);
        let mut total = ev.energy.total;
        total += self
            .gb
            .accumulate_reference(self.vacuum.charges(), at, forces);
        total
    }
}

/// `steps` steps of BAOAB from `start` on `potential`, the buffer `frozen`: the positions,
/// velocities and books, as bits.
fn trajectory(
    potential: &dyn Potential,
    md: &MolecularDynamics,
    start: &[[f64; 3]],
    steps: usize,
) -> Vec<u64> {
    let (mut md, mut at) = (md.clone(), start.to_vec());
    for _ in 0..steps {
        md.step(potential, &mut at, 0.5 * FS);
    }
    let mut b = bits(&at);
    b.extend(bits(md.velocities()));
    b.push(
        md.potential_energy()
            .expect("stepped or prepared")
            .to_bits(),
    );
    b
}

/// **A complex's run in generalized Born is the run the direct evaluation gives, to the bit**:
/// `steps` steps of `c`'s dynamics on its potential, and on [`Direct`] — the same terms, its
/// solvation by the reference — from the same start, seed and bath.
fn the_run_is_the_direct_run(name: &str, b: &Binding, c: &Complex, steps: usize) {
    let vacuum = Complex::with_mobile(b, c.mobile().to_vec(), Solvent::Vacuum).unwrap();
    let gb = c.potential().solvation().expect("a solvent");
    let direct = Direct {
        vacuum: vacuum.potential(),
        gb,
    };
    let md = c.dynamics(
        Bath::Langevin {
            temperature: 300.0,
            friction: 1e12,
            seed: 0x3C4,
        },
        300.0,
        0x3C5,
    );
    let fast = trajectory(c.potential(), &md, c.positions(), steps);
    let slow = trajectory(&direct, &md, c.positions(), steps);
    assert_eq!(fast, slow, "{name}: {steps} steps");
    let n = 3 * c.positions().len();
    assert_ne!(
        fast[..n],
        bits(c.positions())[..],
        "{name}: the run went somewhere"
    );
}

/// **Two models that differ only in the kept terms are equal**, and the kept terms come off.
#[test]
fn the_kept_terms_change_no_model() {
    let (elements, _, at) = aspirin();
    let gb = GeneralizedBorn::new(&elements);
    let mask: Vec<bool> = (0..at.len()).map(|k| k % 2 == 0).collect();
    let kept = gb.clone().with_frozen(&mask, &at);
    assert_eq!(kept, gb);
    assert_eq!(kept.frozen_count(), mask.iter().filter(|&&f| f).count());
    let off = kept.without_frozen();
    assert_eq!(off.frozen_count(), 0);
    assert!(!off.reuses_frozen_terms(&at));
}

/// **The cost of a step under OBC II, and the evaluation the 3b run took, bit for bit.** The
/// 987-atom complex of step 3b — benzene in 181L, a 6 Å mobile zone (322 atoms) in a 10 Å binding,
/// the binding's hydrogens unrelaxed — in OBC II: the evaluation is the reference's at the
/// binding's positions and after 20 steps, and the reuse holds through them. Then the cost:
/// milliseconds per step of dynamics in vacuum and in OBC II, and one evaluation of the reference.
/// **Asserted about the cost only that each saving is taken, by a margin**: with the kept terms the
/// evaluation is under 0.9 of its cost without them (measured 0.75 with the platform's `exp` and
/// `ln`, 0.77–0.79 with the kernel's), and without them under 0.85 of the reference's (0.68, and
/// 0.64–0.67). A fast path that silently stopped being taken gives the same bits, so the default
/// tests cannot see it, and this is what does. The numbers are printed.
#[test]
#[ignore = "QEq on 987 atoms and timed dynamics, about a minute with --release: run with --release -- --ignored --nocapture"]
fn the_3b_complex_is_its_direct_sum_and_costs_measured() {
    let s = the_system();
    let b = Binding::new(s, 10.0 * ANGSTROM).expect("181L at 10 Å");
    let mut gb_complex = Complex::new(&b, 6.0 * ANGSTROM, Solvent::GeneralizedBorn).unwrap();
    let mut vacuum = Complex::new(&b, 6.0 * ANGSTROM, Solvent::Vacuum).unwrap();
    assert_eq!(b.positions().len(), 987);
    assert_eq!(gb_complex.mobile_count(), 322);
    let gb = gb_complex.potential().solvation().unwrap().clone();
    let q = gb_complex.potential().charges().to_vec();
    assert!(is_the_reference(
        "the 3b complex",
        &gb,
        &q,
        gb_complex.positions()
    ));
    the_run_is_the_direct_run("the 3b complex", &b, &gb_complex, 20);
    let bath = Bath::Langevin {
        temperature: 300.0,
        friction: 1e12,
        seed: 0x3C2,
    };
    let dt = 0.5 * FS;
    let mut md = gb_complex.dynamics(bath, 300.0, 0x3C3);
    let mut record = Record::new(&gb_complex, 1 << 40);
    gb_complex.run(&mut md, dt, 20, &mut record);
    assert!(is_the_reference(
        "the 3b complex after 20 steps",
        &gb,
        &q,
        gb_complex.positions()
    ));
    let steps = 200;
    let t = Instant::now();
    gb_complex.run(&mut md, dt, steps, &mut record);
    let wet = t.elapsed().as_secs_f64() * 1e3 / steps as f64;
    assert!(gb.reuses_frozen_terms(gb_complex.positions()));
    let mut md = vacuum.dynamics(bath, 300.0, 0x3C3);
    let mut record = Record::new(&vacuum, 1 << 40);
    vacuum.run(&mut md, dt, 20, &mut record);
    let t = Instant::now();
    vacuum.run(&mut md, dt, 1000, &mut record);
    let dry = t.elapsed().as_secs_f64() * 1e3 / 1000.0;
    let at = gb_complex.positions().to_vec();
    let time = |f: &dyn Fn() -> f64, reps: usize| {
        let t = Instant::now();
        let mut sink = 0.0;
        for _ in 0..reps {
            sink += f();
        }
        assert!(sink.is_finite());
        t.elapsed().as_secs_f64() * 1e3 / reps as f64
    };
    let plain = gb.clone().without_frozen();
    let kept_ms = time(
        &|| gb.accumulate(&q, &at, &mut vec![[0.0; 3]; at.len()]),
        20,
    );
    let plain_ms = time(
        &|| plain.accumulate(&q, &at, &mut vec![[0.0; 3]; at.len()]),
        20,
    );
    let reference_ms = time(
        &|| gb.accumulate_reference(&q, &at, &mut vec![[0.0; 3]; at.len()]),
        5,
    );
    let radii_ms = time(&|| gb.born_radii(&at)[0], 20);
    let radii_reference_ms = time(&|| gb.born_radii_reference(&at)[0], 5);
    println!(
        "one OBC II evaluation of 987 atoms: {kept_ms:.2} ms with the buffer's terms kept, \
         {plain_ms:.2} ms without, {reference_ms:.2} ms the reference; the radii alone \
         {radii_ms:.2} ms against {radii_reference_ms:.2}"
    );
    println!(
        "dynamics: {dry:.3} ms/step in vacuum, {wet:.2} ms/step in OBC II, {:.1} times; \
         {:.2} ns/day",
        wet / dry,
        86_400.0 / (wet * 1e-3) * dt * 1e9
    );
    assert!(
        kept_ms < 0.9 * plain_ms,
        "the kept terms save nothing: {kept_ms:.2} ms with them, {plain_ms:.2} ms without"
    );
    assert!(
        plain_ms < 0.85 * reference_ms,
        "one integral per pair saves nothing: {plain_ms:.2} ms against {reference_ms:.2} ms"
    );
}
