//! **T4 lysozyme L99A and benzene in motion (PDB 181L): a mobile zone around the ligand, a frozen
//! buffer around that, and whether the ligand stays bound.** Step 3b.
//!
//! The default tests check what must hold exactly or as a closed form, on a small complex built
//! from 181L in well under a second: the buffer keeps its bits, the degrees of freedom are
//! `3 N_mobile`, a run is its seed however it is cut into calls, the NVE energy error falls as the
//! step squared, every observable's definition on a hand-built input, the B-factor conversion on
//! a typed value, and the autocorrelation estimator on a series whose correlation time is known.
//!
//! The trajectories take minutes and are ignored by default; each prints its measurements and
//! asserts only what must hold (the buffer's bits, the mobile zone's rule):
//! `cargo test -p pantometry-forcefield --release --test the_complex_in_motion -- --ignored --nocapture`.
//! **Nothing about how the ligand moves, or about the B-factors, is asserted**: those are
//! measurements, reported.

mod protein;

use pantometry_core::conserved::quantity;
use pantometry_core::Rng;
use pantometry_forcefield::complex::{
    self, contacts, mean_square_displacement, rmsd, CONTACT_DISTANCE, LINING_DISTANCE,
};
use pantometry_forcefield::energy::{coulomb, Pair};
use pantometry_forcefield::minimise::KCAL_PER_MOL_ANGSTROM;
use pantometry_forcefield::uff::KCAL_PER_MOL;
use pantometry_forcefield::{
    Bath, Binding, Complex, ComplexError, Component, Element, Estimate, ForceField,
    MolecularDynamics, Record, Selection, Solvent, System,
};
use pantometry_units::BOLTZMANN;
use protein::*;
use std::collections::HashMap;
use std::sync::OnceLock;
use std::time::Instant;

const ANGSTROM: f64 = 1e-10;
const FS: f64 = 1e-15;
const EPS: f64 = f64::EPSILON;

/// The step every run here takes: [`pantometry_forcefield::Molecule::DEFAULT_TIME_STEP`]'s 0.5 fs.
const DT: f64 = 0.5 * FS;

/// The temperature, kelvin.
const KELVIN: f64 = 300.0;

/// The bath's friction for the trajectories: 1 ps⁻¹, gentle enough to leave the vibrations and the
/// ligand's rattling their own time scales.
const FRICTION: f64 = 1e12;

/// A centroid displacement above this, from the crystal position, is the ligand leaving: 3 Å,
/// three times what benzene's rattling in its own cavity is expected to reach and about the
/// cavity's own radius, so a ligand past it is no longer in the site the crystal shows.
const LEAVES: f64 = 3.0 * ANGSTROM;

fn kcal(joules: f64) -> f64 {
    joules / KCAL_PER_MOL
}

/// 181L, built once.
fn the_system() -> &'static System {
    static S: OnceLock<System> = OnceLock::new();
    S.get_or_init(system)
}

// ---------------------------------------------------------------------------------------------
// The small complex the default tests run on
// ---------------------------------------------------------------------------------------------

/// The backbone's atom names: what stays frozen in the small complex.
const BACKBONE: [&str; 6] = ["N", "CA", "C", "O", "H", "HA"];

/// 181L's pocket at 3 Å — benzene and Leu84, Val87, Val111 and Leu121, 82 atoms, two QEq solves in
/// a quarter of a second unoptimised — once.
fn small_binding() -> &'static Binding {
    static B: OnceLock<Binding> = OnceLock::new();
    B.get_or_init(|| Binding::new(the_system(), 3.0 * ANGSTROM).expect("181L at 3 Å"))
}

/// The small complex's mobile mask: benzene and the four side chains, every backbone atom frozen.
/// At 3 Å every residue's backbone is bonded across the cut, so a cutoff's zone of whole residues
/// would be refused; the side chains are not, and they are a real mobile zone of 55 atoms — five
/// elements' worth of protein free beside the ligand.
fn small_mask(b: &Binding) -> Vec<bool> {
    let s = the_system();
    (0..b.positions().len())
        .map(|k| k >= b.pocket_len() || !BACKBONE.contains(&s.atom_name(b.system_atoms()[k])))
        .collect()
}

fn small(solvent: Solvent) -> Complex {
    let b = small_binding();
    Complex::with_mobile(b, small_mask(b), solvent).expect("no side-chain atom is at the cut")
}

/// The small complex, its mobile zone minimised once (vacuum).
fn small_minimised() -> &'static Complex {
    static C: OnceLock<Complex> = OnceLock::new();
    C.get_or_init(|| {
        let mut c = small(Solvent::Vacuum);
        let p = c.minimise(5000, Binding::HYDROGEN_TOLERANCE);
        assert_eq!(p.status, pantometry_forcefield::Status::Converged, "{p:?}");
        c
    })
}

fn langevin(seed: u64, friction: f64) -> Bath {
    Bath::Langevin {
        temperature: KELVIN,
        friction,
        seed,
    }
}

fn bits(at: &[[f64; 3]]) -> Vec<u64> {
    at.iter().flatten().map(|v| v.to_bits()).collect()
}

fn frame_bits(r: &Record) -> Vec<u64> {
    r.frames()
        .iter()
        .flat_map(|f| {
            [
                f.step,
                f.rmsd.to_bits(),
                f.site_rmsd.to_bits(),
                f.turn.to_bits(),
                f.centroid_displacement.to_bits(),
                f.cavity_distance.to_bits(),
                f.contacts as u64,
                f.van_der_waals.to_bits(),
                f.electrostatic.to_bits(),
                f.temperature.to_bits(),
                f.books.to_bits(),
            ]
        })
        .collect()
}

// ---------------------------------------------------------------------------------------------
// The zone
// ---------------------------------------------------------------------------------------------

/// **A zone that would leave no buffer, or free an atom bonded across the cut, is refused by
/// name.** At 3 Å every residue of 181L's pocket is bonded to residues the cut left out — checked
/// here by building each residue's bonds from the system, not the binding — so a cutoff's zone of
/// whole residues is refused at its first residue. `Binding::bonded_outside` is held atom by atom
/// to that scan, and then **every** atom it marks — each residue's backbone `N` and `C`, bonded to
/// the neighbours the cut left out — is freed alone beside the side chains, and each is refused,
/// naming its own residue. A cutoff below every residue's heavy atoms gives the ligand alone.
#[test]
fn a_zone_without_a_buffer_or_across_the_cut_is_refused() {
    let b = small_binding();
    let s = the_system();
    for mobile in [0.0, -1.0, 3.0, 4.0] {
        assert!(
            matches!(
                Complex::new(b, mobile * ANGSTROM, Solvent::Vacuum),
                Err(ComplexError::NoBuffer { .. })
            ),
            "{mobile} Å"
        );
    }
    // Every residue at 3 Å has a backbone bond leaving the binding, read from the whole system,
    // and the binding's own mask is that scan, atom by atom.
    let kept: Vec<usize> = b.system_atoms().to_vec();
    let leaves = |k: usize| {
        let i = b.system_atoms()[k];
        s.component().neighbours(i).any(|(j, _)| !kept.contains(&j))
    };
    for k in 0..b.positions().len() {
        assert_eq!(b.bonded_outside()[k], leaves(k), "atom {k}");
    }
    for (span, label) in b.residue_spans().iter().zip(b.residues()) {
        assert!(
            span.clone().any(leaves),
            "{label} is whole inside the binding"
        );
    }
    let refused = Complex::new(b, 2.99 * ANGSTROM, Solvent::Vacuum).unwrap_err();
    println!("{refused}");
    assert!(
        matches!(&refused, ComplexError::MobileAtTheCut { residue } if b.residues().contains(residue)),
        "{refused}"
    );
    let mut names = Vec::new();
    for (span, label) in b.residue_spans().iter().zip(b.residues()) {
        for k in span.clone().filter(|&k| leaves(k)) {
            names.push(s.atom_name(b.system_atoms()[k]).to_string());
            let mut mask = small_mask(b);
            mask[k] = true;
            let refused = Complex::with_mobile(b, mask, Solvent::Vacuum).unwrap_err();
            assert!(
                matches!(&refused, ComplexError::MobileAtTheCut { residue } if residue == label),
                "freeing {label} {}: {refused}",
                s.atom_name(b.system_atoms()[k])
            );
        }
    }
    names.sort();
    names.dedup();
    assert_eq!(names, ["C", "N"], "the cut goes through peptide bonds only");
    let alone = Complex::new(b, 0.5 * ANGSTROM, Solvent::Vacuum).unwrap();
    assert_eq!(alone.mobile_count(), 12, "benzene alone");
    assert!(alone.mobile_residues().is_empty());
    assert_eq!(
        small(Solvent::Vacuum).mobile_count(),
        12 + b.pocket_len() - 4 * BACKBONE.len()
    );
}

/// **A cutoff's zone is the pocket's rule, at debug speed.** The mask [`Complex::zone`] frees in
/// the 3 Å binding is held, residue by residue, to the pocket `Binding::new` itself cuts at the
/// same radius: at 2.7 Å that is Val87 alone (its nearest heavy atom is 2.60 Å from benzene's
/// atoms, the other three's 2.89–2.97), and at 2 Å there is none, which `Binding::new` refuses as
/// an empty pocket. The radii are chosen so that the rule's two parts each matter: a radius half
/// as large again reaches all four residues, and a rule that let a residue's hydrogens count
/// would take in Val111, one of whose hydrogens is 1.89 Å from a benzene atom.
#[test]
fn a_cutoffs_zone_is_the_pockets_rule_at_debug_speed() {
    let b = small_binding();
    for radius in [2.0, 2.7] {
        let zone = Complex::zone(b, radius * ANGSTROM);
        let freed: Vec<String> = b
            .residue_spans()
            .iter()
            .zip(b.residues())
            .filter(|(span, _)| (*span).clone().any(|k| zone[k]))
            .map(|(_, r)| r.clone())
            .collect();
        for (span, _) in b.residue_spans().iter().zip(b.residues()) {
            assert!(
                span.clone().all(|k| zone[k]) || span.clone().all(|k| !zone[k]),
                "whole residues"
            );
        }
        assert!(
            b.ligand_range().all(|k| zone[k]),
            "the ligand is always free"
        );
        let pocket: Vec<String> = match Binding::new(the_system(), radius * ANGSTROM) {
            Ok(p) => p.residues().to_vec(),
            Err(pantometry_forcefield::BindingError::EmptyPocket { .. }) => Vec::new(),
            Err(e) => panic!("{e}"),
        };
        println!("{radius} Å: {freed:?}");
        assert_eq!(freed, pocket, "{radius} Å");
    }
    assert_eq!(
        Complex::zone(b, 2.7 * ANGSTROM)
            .iter()
            .filter(|&&m| m)
            .count(),
        12 + 16,
        "benzene and Val87's sixteen atoms"
    );
}

/// **The buffer keeps its bits, and the degrees of freedom are `3 N_mobile`.** 300 steps under a
/// bath: every frozen atom's position is the same bits as before, its velocity is exactly zero,
/// and every mobile atom has moved. With atoms frozen nothing is removed from the draw, so the
/// count is three per mobile atom, under the bath and without it.
#[test]
fn the_buffer_keeps_its_bits_and_the_mobile_atoms_have_the_degrees_of_freedom() {
    let mut c = small_minimised().clone();
    let start = c.positions().to_vec();
    let mut md = c.dynamics(langevin(11, FRICTION), KELVIN, 12);
    assert_eq!(md.degrees_of_freedom(), 3 * c.mobile_count());
    for (m, e) in md.masses().iter().zip(small_binding().elements()) {
        assert_eq!(*m, e.mass(), "each atom its own element's mass");
    }
    assert_eq!(
        c.dynamics(Bath::Isolated, KELVIN, 12).degrees_of_freedom(),
        3 * c.mobile_count()
    );
    let mut record = Record::new(&c, 50);
    c.run(&mut md, DT, 300, &mut record);
    assert_eq!(record.frames().len(), 6);
    for (k, &m) in c.mobile().iter().enumerate() {
        let (now, then) = (c.positions()[k], start[k]);
        if m {
            assert_ne!(bits(&[now]), bits(&[then]), "mobile atom {k} did not move");
        } else {
            assert_eq!(bits(&[now]), bits(&[then]), "frozen atom {k} moved");
            assert_eq!(md.velocities()[k], [0.0; 3]);
            assert_eq!(record.mean_square_fluctuation(k), 0.0);
        }
    }
    assert!(c.mobile().iter().any(|&m| !m), "the buffer is not empty");
}

/// **One seed, one run, to the bit — however it is cut into calls.** 120 steps in one call, in
/// four of 30, in 120 of one, and from a complex rebuilt from the file: the same positions,
/// velocities, frames and fluctuations, bit for bit. A frame interval that does not divide a
/// chunk's length (every 7 steps here) still records the same frames, because a frame is keyed
/// by the dynamics' own step count. Another seed is another run.
#[test]
fn a_run_is_its_seed_however_it_is_cut_into_calls() {
    let base = small_minimised().clone();
    let go = |c: &Complex, seed: u64, chunks: &[usize]| {
        let mut c = c.clone();
        let mut md = c.dynamics(langevin(seed, FRICTION), KELVIN, seed + 1);
        let mut r = Record::new(&c, 7);
        for &n in chunks {
            c.run(&mut md, DT, n, &mut r);
        }
        let mut b = bits(c.positions());
        b.extend(bits(md.velocities()));
        b.extend(frame_bits(&r));
        b.extend((0..c.positions().len()).map(|k| r.mean_square_fluctuation(k).to_bits()));
        b
    };
    let one = go(&base, 3, &[120]);
    assert_eq!(one, go(&base, 3, &[30; 4]), "four calls of 30");
    assert_eq!(one, go(&base, 3, &[1; 120]), "120 calls of one");
    let mut rebuilt = Complex::with_mobile(
        &Binding::new(&system(), 3.0 * ANGSTROM).unwrap(),
        small_mask(small_binding()),
        Solvent::Vacuum,
    )
    .unwrap();
    rebuilt.minimise(5000, Binding::HYDROGEN_TOLERANCE);
    assert_eq!(one, go(&rebuilt, 3, &[120]), "rebuilt from the file");
    assert_ne!(one, go(&base, 4, &[120]), "another seed");
}

/// The sum, the sum of magnitudes and the count of the terms of `ff` at `at` whose atoms `keep`
/// accepts, each through its public per-term energy.
fn magnitude(ff: &ForceField, at: &[[f64; 3]], keep: impl Fn(&[usize]) -> bool) -> (f64, f64, f64) {
    let (mut sum, mut size, mut n) = (0.0, 0.0, 0.0);
    let mut add = |e: f64| {
        sum += e;
        size += e.abs();
        n += 1.0;
    };
    let r = |i: usize, j: usize| {
        ((at[i][0] - at[j][0]).powi(2)
            + (at[i][1] - at[j][1]).powi(2)
            + (at[i][2] - at[j][2]).powi(2))
        .sqrt()
    };
    for t in ff.stretches().iter().filter(|t| keep(&t.atoms)) {
        add(t.energy(r(t.atoms[0], t.atoms[1])));
    }
    for t in ff.bends().iter().filter(|t| keep(&t.atoms)) {
        add(t.energy_at(at));
    }
    for t in ff.torsions().iter().filter(|t| keep(&t.atoms)) {
        add(t.energy_at(at));
    }
    for t in ff.inversions().iter().filter(|t| keep(&t.atoms)) {
        add(t.energy_at(at));
    }
    let q = ff.charges();
    for p in ff.pairs().iter().filter(|t| keep(&t.atoms)) {
        let [i, j] = p.atoms;
        add(p.energy(r(i, j)));
        if q[i] != 0.0 && q[j] != 0.0 {
            add(coulomb(q[i], q[j], r(i, j)));
        }
    }
    (sum, size, n)
}

/// **The potential is the whole complex's on every mobile atom, bit for bit, in vacuum and in
/// generalized Born.** It is the complex's force field with the terms among frozen atoms left out,
/// in the same order, so a mobile atom's force is the same sum in the same order; under OBC II the
/// solvation of every atom is added to both.
///
/// **And the frozen-only terms are gone.** No term of the potential has all of its atoms frozen;
/// the whole complex has such terms (bonds, angles and pairs inside the frozen backbones); and at
/// both geometries the whole complex's energy less the potential's is those terms' sum, computed
/// here term by term through the public per-term energies. The allowance is traced: each of the
/// three sums carries at most `n ε Σ|t|` of rounding (first order in ε), and the subtraction one ε
/// of its larger operand. OBC II is the same value in both evaluations, the same function of the
/// same inputs, so it enters only as one more term of `Σ|t|`.
#[test]
fn the_potential_is_the_complex_on_the_mobile_atoms() {
    let b = small_binding();
    let mut moved = small_minimised().positions().to_vec();
    for (k, p) in moved.iter_mut().enumerate() {
        if small_minimised().mobile()[k] {
            p[0] += 0.05 * ANGSTROM;
        }
    }
    for solvent in [Solvent::Vacuum, Solvent::GeneralizedBorn] {
        let c = small(solvent);
        let whole = match solvent {
            Solvent::Vacuum => b.force_field().clone(),
            Solvent::GeneralizedBorn => b.force_field().clone().with_generalized_born(),
        };
        assert_eq!(
            c.potential().solvation().is_some(),
            solvent == Solvent::GeneralizedBorn
        );
        let frozen = |atoms: &[usize]| atoms.iter().all(|&a| !c.mobile()[a]);
        let z = c.potential();
        let any_frozen_only = z.stretches().iter().any(|t| frozen(&t.atoms))
            || z.bends().iter().any(|t| frozen(&t.atoms))
            || z.torsions().iter().any(|t| frozen(&t.atoms))
            || z.inversions().iter().any(|t| frozen(&t.atoms))
            || z.pairs().iter().any(|t| frozen(&t.atoms));
        assert!(
            !any_frozen_only,
            "{solvent:?}: a term among frozen atoms alone was kept"
        );
        for at in [c.positions(), &moved[..]] {
            let (zone, all) = (c.potential().evaluate(at), whole.evaluate(at));
            for (k, &m) in c.mobile().iter().enumerate() {
                if m {
                    assert_eq!(
                        bits(&[zone.forces[k]]),
                        bits(&[all.forces[k]]),
                        "{solvent:?}, atom {k}"
                    );
                }
            }
            let (left_out, lo_size, lo_n) = magnitude(&whole, at, frozen);
            let (_, all_size, all_n) = magnitude(&whole, at, |_| true);
            let (_, zone_size, zone_n) = magnitude(z, at, |_| true);
            let g = whole
                .solvation()
                .map_or(0.0, |gb| gb.energy(whole.charges(), at).abs());
            assert!(
                lo_n > 0.0 && left_out != 0.0,
                "the whole complex has frozen-only terms"
            );
            let gap = all.energy.total - zone.energy.total;
            let allowance = EPS
                * ((all_n + 1.0) * (all_size + g)
                    + (zone_n + 1.0) * (zone_size + g)
                    + lo_n * lo_size
                    + all.energy.total.abs().max(zone.energy.total.abs()));
            println!(
                "{solvent:?}: left out {lo_n} terms, {:.6} kcal/mol; gap − sum {:.2e} against {:.2e}",
                kcal(left_out),
                gap - left_out,
                allowance
            );
            assert!(
                (gap - left_out).abs() <= allowance,
                "{solvent:?}: the gap {gap} is not the frozen-only terms' {left_out}"
            );
        }
    }
}

/// **Every frame is its positions.** Twenty steps with a frame after each: every observable
/// recomputed here from the positions the step left — the RMSD, the site RMSD and the centroid
/// against the binding's crystal positions, the cavity centre from the lining atoms where they
/// are now, the contacts, the cross terms summed pair by pair here, the temperature after `O`,
/// the books, and the turn through the same function, so for the turn only the wiring (its
/// definition is held by the closed forms in `the_site_rmsd_is_nearest_and_the_turn_is_the_angle`)
/// — and every atom's mean-square fluctuation from the twenty positions, to the
/// rounding of a sum of twenty.
#[test]
fn every_frame_is_its_positions() {
    let mut c = small_minimised().clone();
    let b = small_binding();
    let crystal = b.crystal_positions().to_vec();
    let mut md = c.dynamics(langevin(31, FRICTION), KELVIN, 32);
    let mut r = Record::new(&c, 1);
    let mut seen: Vec<Vec<[f64; 3]>> = Vec::new();
    let protein: Vec<usize> = (0..b.pocket_len())
        .filter(|&k| b.elements()[k] != Element::H)
        .collect();
    for _ in 0..20 {
        c.run(&mut md, DT, 1, &mut r);
        let at = c.positions();
        seen.push(at.to_vec());
        let f = r.frames().last().unwrap();
        let pick =
            |ks: &[usize], from: &[[f64; 3]]| ks.iter().map(|&k| from[k]).collect::<Vec<_>>();
        let lig = pick(c.ligand_heavy_atoms(), at);
        let lig0 = pick(c.ligand_heavy_atoms(), &crystal);
        let near = |a: [f64; 3], b: [f64; 3]| {
            ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
        };
        let rms = (lig
            .iter()
            .zip(&lig0)
            .map(|(p, q)| near(*p, *q).powi(2))
            .sum::<f64>()
            / lig.len() as f64)
            .sqrt();
        let scale = 64.0 * EPS * 30.0 * ANGSTROM;
        assert!((f.rmsd - rms).abs() <= scale);
        let site = (lig
            .iter()
            .map(|p| {
                lig0.iter()
                    .map(|q| near(*p, *q).powi(2))
                    .fold(f64::INFINITY, f64::min)
            })
            .sum::<f64>()
            / lig.len() as f64)
            .sqrt();
        assert!((f.site_rmsd - site).abs() <= scale);
        assert_eq!(f.turn, complex::in_plane_turn(&lig, &lig0));
        let (cl, c0) = (complex::centroid(&lig), complex::centroid(&lig0));
        assert!((f.centroid_displacement - near(cl, c0)).abs() <= scale);
        let lining = complex::centroid(&pick(c.lining(), at));
        assert!((f.cavity_distance - near(cl, lining)).abs() <= scale);
        let count = lig
            .iter()
            .map(|l| {
                protein
                    .iter()
                    .filter(|&&k| near(*l, at[k]) < CONTACT_DISTANCE)
                    .count()
            })
            .sum::<usize>();
        assert_eq!(f.contacts, count);
        // The cross terms summed here, every protein–ligand pair from the types and the charges
        // with no exclusion list, as `benzene_in_its_pocket.rs` builds them: the same terms in
        // another order, so to `n ε Σ|t|` of each sum.
        let (t, q) = (b.types(), b.charges());
        let (mut vdw, mut elec, mut size, mut n) = (0.0, 0.0, 0.0, 0.0);
        for p in 0..b.pocket_len() {
            for l in b.ligand_range() {
                let r = near(at[p], at[l]);
                let (e, c) = (
                    Pair::new([p, l], t[p], t[l]).energy(r),
                    coulomb(q[p], q[l], r),
                );
                vdw += e;
                elec += c;
                size += e.abs() + c.abs();
                n += 1.0;
            }
        }
        assert!(
            (f.van_der_waals - vdw).abs() <= 2.0 * n * EPS * size,
            "{} {vdw}",
            f.van_der_waals
        );
        assert!(
            (f.electrostatic - elec).abs() <= 2.0 * n * EPS * size,
            "{} {elec}",
            f.electrostatic
        );
        assert_eq!(f.temperature, md.half_step_temperature());
        assert_eq!(f.books, md.ledger().get(quantity::ENERGY).unwrap());
        assert_eq!(f.step, seen.len() as u64);
    }
    for k in 0..crystal.len() {
        let n = seen.len() as f64;
        let mean: Vec<f64> = (0..3)
            .map(|d| seen.iter().map(|p| p[k][d]).sum::<f64>() / n)
            .collect();
        let msf = seen
            .iter()
            .map(|p| (0..3).map(|d| (p[k][d] - mean[d]).powi(2)).sum::<f64>())
            .sum::<f64>()
            / n;
        // The record sums displacements from the crystal and subtracts; its rounding is a few ε
        // of the squared displacement, which is at most an ångström's square here.
        let allowance = 64.0 * EPS * ANGSTROM * ANGSTROM;
        assert!(
            (r.mean_square_fluctuation(k) - msf).abs() <= allowance,
            "atom {k}: {} against {msf}",
            r.mean_square_fluctuation(k)
        );
        if c.mobile()[k] {
            assert!(msf > 0.0);
        }
    }
}

/// The RMS departure of the books from their start over `femtoseconds` of NVE at step `dt` from
/// the minimised small complex, velocities at 300 K from `seed`.
fn nve_error(dt: f64, femtoseconds: f64, seed: u64) -> f64 {
    let mut c = small_minimised().clone();
    let mut md = c.dynamics(Bath::Isolated, KELVIN, seed);
    // The books at the start: the kinetic energy of the draw, the potential counted from itself.
    md.prepare(c.potential(), c.positions());
    let e0 = md.ledger().get(quantity::ENERGY).unwrap();
    let steps = (femtoseconds * FS / dt).round() as usize;
    let mut r = Record::new(&c, 1);
    c.run(&mut md, dt, steps, &mut r);
    let sq: f64 = r.frames().iter().map(|f| (f.books - e0).powi(2)).sum();
    (sq / r.frames().len() as f64).sqrt()
}

/// `log₂(error(1 fs)/error(0.5 fs))` and `log₂(error(0.5)/error(0.25))` over 100 fs from `seed`.
fn nve_orders(seed: u64) -> [f64; 2] {
    let e: Vec<f64> = [1.0, 0.5, 0.25]
        .iter()
        .map(|&h| nve_error(h * FS, 100.0, seed))
        .collect();
    [(e[0] / e[1]).log2(), (e[1] / e[2]).log2()]
}

/// **NVE on the mobile zone: the energy error falls as the step squared.** The frozen buffer is a
/// fixed external potential, so the mobile atoms' energy is conserved by the exact dynamics and
/// velocity Verlet conserves a shadow of it, `H + O(h²)`: the RMS departure from the start falls
/// fourfold per halving. Over 100 fs at 1, 0.5 and 0.25 fs from one start.
///
/// **The bound, 0.15, is earned** from sixteen starts (`the_nve_orders_measured_over_many_starts`):
/// their 32 orders have mean 2.041 and standard deviation 0.025, range [2.015, 2.070], so the
/// bound is the mean plus four of them (2.14), above the worst.
///
/// **What it cannot see is a potential that is wrong and conservative**: a zone that left the
/// buffer's pull out would conserve its own energy just as well. That is
/// `the_potential_is_the_complex_on_the_mobile_atoms`' job.
#[test]
fn the_mobile_zones_energy_error_falls_as_the_step_squared() {
    let o = nve_orders(21);
    println!("NVE orders {:.4}, {:.4}", o[0], o[1]);
    for order in o {
        assert!((order - 2.0).abs() < 0.15, "order {order}");
    }
}

/// The spread of the NVE orders over sixteen starts: what the bound above is earned from.
#[test]
#[ignore = "a measurement that prints and asserts nothing: run with --release -- --ignored --nocapture"]
fn the_nve_orders_measured_over_many_starts() {
    let all: Vec<f64> = (100..116).flat_map(nve_orders).collect();
    let mean = all.iter().sum::<f64>() / all.len() as f64;
    let sd = (all.iter().map(|o| (o - mean).powi(2)).sum::<f64>() / (all.len() - 1) as f64).sqrt();
    let (lo, hi) = all
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), o| {
            (l.min(*o), h.max(*o))
        });
    println!(
        "{} orders: mean {mean:.4}, sd {sd:.4}, range [{lo:.4}, {hi:.4}]",
        all.len()
    );
}

// ---------------------------------------------------------------------------------------------
// The observables, on inputs whose answers are known
// ---------------------------------------------------------------------------------------------

/// **The RMSD without superposition, on known motions.** A translation by `t` gives `|t|`, to the
/// rounding of a difference of coordinates: each component of `(p + t) − p` is within one ε of
/// `|p|` of exact, so the bound is a few ε of the largest coordinate. A ring of radius `a` turned
/// by `θ` about its own axis moves every atom by `2a sin(θ/2)`, so its RMSD is that — for a
/// six-membered ring turned by 60°, which lays it on itself with its atoms permuted, the RMSD is
/// `a` exactly as a closed form, though the ring looks the same. That is how a benzene spinning
/// in its cavity registers.
#[test]
fn the_rmsd_without_superposition_is_the_motion() {
    let a = 1.39 * ANGSTROM;
    let centre = [25.0 * ANGSTROM, 5.0 * ANGSTROM, 4.0 * ANGSTROM];
    let ring = |turn: f64| -> Vec<[f64; 3]> {
        (0..6)
            .map(|k| {
                let phi = k as f64 * std::f64::consts::PI / 3.0 + turn;
                [
                    centre[0] + a * phi.cos(),
                    centre[1] + a * phi.sin(),
                    centre[2],
                ]
            })
            .collect()
    };
    let start = ring(0.0);
    let t = [3.0e-11, 4.0e-11, 12.0e-11];
    let moved: Vec<[f64; 3]> = start
        .iter()
        .map(|p| [p[0] + t[0], p[1] + t[1], p[2] + t[2]])
        .collect();
    let largest = 27.0 * ANGSTROM;
    let r = rmsd(&moved, &start);
    assert!((r - 13.0e-11).abs() <= 8.0 * EPS * largest, "{r}");
    assert_eq!(rmsd(&start, &start), 0.0);
    for theta in [0.3, std::f64::consts::PI / 3.0, std::f64::consts::PI] {
        let r = rmsd(&ring(theta), &start);
        let exact = 2.0 * a * (theta / 2.0).sin();
        // The ring's positions carry the rounding of `cos`, `sin` and the centre's addition: a
        // few ε of the largest coordinate per component.
        assert!(
            (r - exact).abs() <= 16.0 * EPS * largest,
            "θ = {theta}: {r} against {exact}"
        );
    }
    // Blind to the permutation: a turn by 60° lays the ring on its own sites, and a turn by θ < 30°
    // leaves every atom nearest its own site, at `2a sin(θ/2)` from it; by 30° it is halfway.
    let carbons = [Element::C; 6];
    let site = |theta: f64| complex::site_rmsd(&ring(theta), &start, &carbons);
    assert!(site(std::f64::consts::PI / 3.0) <= 16.0 * EPS * largest);
    for theta in [0.1, 0.4, std::f64::consts::PI / 6.0] {
        let exact = 2.0 * a * (theta / 2.0).sin();
        assert!(
            (site(theta) - exact).abs() <= 16.0 * EPS * largest,
            "θ = {theta}"
        );
    }
    // A site of another element is not a site: with one of the six called nitrogen and the ring
    // turned by 60°, the nitrogen's only site is its own, an edge `a` away, and the carbon that
    // lands on the nitrogen's site is an edge from the nearest carbon one. Two atoms at `a`, four
    // at zero: `a √(2/6)`.
    let mut kinds = carbons;
    kinds[0] = Element::N;
    let one = complex::site_rmsd(&ring(std::f64::consts::PI / 3.0), &start, &kinds);
    assert!(
        (one - a / 3f64.sqrt()).abs() <= 16.0 * EPS * largest,
        "{one}"
    );
    let c = complex::centroid(&start);
    for k in 0..3 {
        assert!((c[k] - centre[k]).abs() <= 4.0 * EPS * largest);
    }
}

/// **The site RMSD is nearest, not one-to-one, and so it has a ceiling.** A ring of radius `a`
/// turned by `θ` is at `2a sin(ψ/2)` from its nearest sites, `ψ = |θ|` folded to [0°, 30°], so a
/// ring turning uniformly reads, on average over the turn, `⟨site RMSD²⟩ = 2a²⟨1 − cos ψ⟩ =
/// 2a²(1 − 3/π)` — 0.174 Å² at 1.39 Å — however freely it turns. Averaged here over 6000 equal
/// turns at the cells' midpoints, which put the kinks at 30° on cell edges: the midpoint rule's
/// error is then at most `h² max|f″| / 24 = h² · 2a²/24` for a mean, plus the sum's rounding.
/// And two atoms on one site both read zero, where a one-to-one assignment would put one of them
/// on the other site, `d/√2` RMS.
///
/// The in-plane turn of the same ring is `θ` itself, translated or not — a turn of 60°, which the
/// site RMSD cannot see, is 60° here.
#[test]
fn the_site_rmsd_is_nearest_and_the_turn_is_the_angle() {
    let a = 1.39 * ANGSTROM;
    let centre = [25.0 * ANGSTROM, 5.0 * ANGSTROM, 4.0 * ANGSTROM];
    let ring = |turn: f64, shift: f64| -> Vec<[f64; 3]> {
        (0..6)
            .map(|k| {
                let phi = k as f64 * std::f64::consts::PI / 3.0 + turn;
                [
                    centre[0] + a * phi.cos() + shift,
                    centre[1] + a * phi.sin() - shift,
                    centre[2] + 0.5 * shift,
                ]
            })
            .collect()
    };
    let start = ring(0.0, 0.0);
    let carbons = [Element::C; 6];
    let n = 6000;
    let h = 2.0 * std::f64::consts::PI / n as f64;
    let mean = (0..n)
        .map(|k| complex::site_rmsd(&ring((k as f64 + 0.5) * h, 0.0), &start, &carbons).powi(2))
        .sum::<f64>()
        / n as f64;
    let exact = 2.0 * a * a * (1.0 - 3.0 / std::f64::consts::PI);
    let allowance = h * h * 2.0 * a * a / 24.0 + 1e-11 * a * a;
    println!(
        "uniform turn: ⟨site RMSD²⟩ = {:.6} Å², 2a²(1 − 3/π) = {:.6} Å²",
        mean / (ANGSTROM * ANGSTROM),
        exact / (ANGSTROM * ANGSTROM)
    );
    assert!((mean - exact).abs() <= allowance, "{mean} against {exact}");
    // Two atoms on one site.
    let sites = [[0.0; 3], [1.4 * ANGSTROM, 0.0, 0.0]];
    let both = [[0.0; 3], [0.0; 3]];
    assert_eq!(complex::site_rmsd(&both, &sites, &[Element::C; 2]), 0.0);
    assert!((rmsd(&both, &sites) - 1.4 * ANGSTROM / 2f64.sqrt()).abs() <= 4.0 * EPS * ANGSTROM);
    // The turn.
    let scale = 64.0 * EPS * 30.0 * ANGSTROM / a;
    for theta in [0.0, 0.2, -0.7, std::f64::consts::PI / 3.0, 2.5] {
        for shift in [0.0, 0.8 * ANGSTROM] {
            let t = complex::in_plane_turn(&ring(theta, shift), &start);
            // The ring's normal is +z or −z as the eigenvector falls; the angle's sign with it.
            assert!(
                (t.abs() - f64::abs(theta)).abs() <= scale,
                "θ = {theta}, shift {shift}: {t}"
            );
        }
    }
    // Not rigid: one atom of the six turned by α about the centre, the rest where they were. With
    // the centroid held, the circular mean of one α and five zeros would be
    // `atan2(sin α, 5 + cos α)`, which a turn read off one atom, or averaged as plain angles,
    // does not give.
    let alpha: f64 = 0.6;
    let mut one = start.clone();
    one[3] = [
        centre[0] + a * (std::f64::consts::PI + alpha).cos(),
        centre[1] + a * (std::f64::consts::PI + alpha).sin(),
        centre[2],
    ];
    // The centroid moves with the atom; the closed form is for angles about the moved centroid,
    // so it is computed here from the same projected offsets, by hand.
    let c1 = complex::centroid(&one);
    let (mut sum_sin, mut sum_cos) = (0.0, 0.0);
    for (p, q) in one.iter().zip(&start) {
        let u = [p[0] - c1[0], p[1] - c1[1]];
        let u0 = [q[0] - centre[0], q[1] - centre[1]];
        let angle = (u0[0] * u[1] - u0[1] * u[0]).atan2(u0[0] * u[0] + u0[1] * u[1]);
        sum_sin += angle.sin();
        sum_cos += angle.cos();
    }
    let expected = sum_sin.atan2(sum_cos);
    let got = complex::in_plane_turn(&one, &start);
    assert!(
        (got.abs() - expected.abs()).abs() <= scale,
        "one atom turned: {got} against {expected}"
    );
    // And it is neither the one atom's angle nor the plain mean of the angles.
    assert!((got.abs() - alpha).abs() > 0.1 && (got.abs() - alpha / 6.0).abs() > 1e-3);
    let t1 = complex::in_plane_turn(&ring(0.2, 0.0), &start);
    let t2 = complex::in_plane_turn(&ring(-0.7, 0.0), &start);
    assert!(t1 * t2 < 0.0, "opposite turns have opposite signs");
}

/// **Contacts on a constructed geometry: strictly below the cutoff, every pair counted.** One
/// atom at the origin and others on the axes at 2, 3.9, 4.0 and 4.1 Å: two contacts, the one at
/// exactly 4 Å not one. Then a 2 × 3 grid at 1 Å spacing against itself at 1.5 Å: each point is
/// within 1.5 Å of itself and its edge neighbours (1 Å) and not its diagonal ones (√2 Å is below
/// 1.5, so those count too) — counted here by hand as 6 + 2·7 + 2·4 = 28.
#[test]
fn contacts_are_counted_on_a_constructed_geometry() {
    let origin = [[0.0; 3]];
    let others = [
        [2.0 * ANGSTROM, 0.0, 0.0],
        [0.0, 3.9 * ANGSTROM, 0.0],
        [0.0, 0.0, 4.0 * ANGSTROM],
        [-4.1 * ANGSTROM, 0.0, 0.0],
    ];
    assert_eq!(contacts(&origin, &others, CONTACT_DISTANCE), 2);
    assert_eq!(contacts(&others, &origin, CONTACT_DISTANCE), 2);
    let grid: Vec<[f64; 3]> = (0..2)
        .flat_map(|i| (0..3).map(move |j| [i as f64 * ANGSTROM, j as f64 * ANGSTROM, 0.0]))
        .collect();
    // Ordered pairs (p, q), p = q included: 6 with themselves; edge neighbours, 7 unordered
    // pairs (3 across, 4 along), twice; diagonals at √2 Å, 4 unordered pairs, twice. The pairs
    // two apart along the long side are 2 Å and do not count.
    assert_eq!(contacts(&grid, &grid, 1.5 * ANGSTROM), 6 + 2 * 7 + 2 * 4);
}

/// The atoms of 181L as the file has them: `(chain, residue number, atom name) → (residue name,
/// position in metres, B in Å²)`, read with string operations and nothing from the crate.
type FileAtoms = HashMap<(char, i32, String), (String, [f64; 3], f64)>;

fn atoms_of_the_file(pdb: &str) -> FileAtoms {
    atom_lines(pdb)
        .map(|l| {
            let num = |a, b| columns(l, a, b).trim().parse::<f64>().unwrap();
            let chain = columns(l, 22, 22).chars().next().unwrap();
            let number: i32 = columns(l, 23, 26).trim().parse().unwrap();
            let name = columns(l, 13, 16).trim().to_string();
            let at = [
                num(31, 38) * ANGSTROM,
                num(39, 46) * ANGSTROM,
                num(47, 54) * ANGSTROM,
            ];
            (
                (chain, number, name),
                (columns(l, 18, 20).trim().to_string(), at, num(61, 66)),
            )
        })
        .collect()
}

/// **At the crystal pose the observables are the file's.** The small complex, not minimised,
/// sits at the crystal coordinates: the RMSD and the centroid displacement are exactly zero, the
/// interaction is the binding's own cross sum bit for bit, and the contacts and the cavity distance
/// are what the file gives when they are recomputed here from its `ATOM` and `HETATM` records —
/// over the residues the binding holds, since a 3 Å pocket does not hold every atom within 5 Å.
#[test]
fn at_the_crystal_pose_the_observables_are_the_files() {
    let c = small(Solvent::Vacuum);
    let b = small_binding();
    let md = c.dynamics(Bath::Isolated, KELVIN, 1);
    let f = c.observe(&md);
    assert_eq!(f.rmsd, 0.0);
    assert_eq!(f.site_rmsd, 0.0);
    assert_eq!(f.turn, 0.0);
    assert_eq!(f.centroid_displacement, 0.0);
    let i = b.interaction();
    assert_eq!(f.van_der_waals.to_bits(), i.van_der_waals.to_bits());
    assert_eq!(f.electrostatic.to_bits(), i.electrostatic.to_bits());
    assert!(f.temperature.is_nan(), "no bath, no step: no temperature");

    let file = atoms_of_the_file(PDB_181L);
    let held: Vec<i32> = b
        .residues()
        .iter()
        .map(|r| r[5..].parse().unwrap())
        .collect();
    let benzene: Vec<[f64; 3]> = file
        .values()
        .filter(|(res, _, _)| res == "BNZ")
        .map(|(_, at, _)| *at)
        .collect();
    assert_eq!(benzene.len(), 6);
    let protein: Vec<[f64; 3]> = file
        .iter()
        .filter(|((_, n, _), (res, _, _))| {
            // 181L has no hydrogens, so every protein record is a heavy atom.
            res != "BNZ" && res != "HOH" && res != "CL" && res != "HED" && held.contains(n)
        })
        .map(|(_, (_, at, _))| *at)
        .collect();
    let mut expected = 0;
    for l in &benzene {
        for p in &protein {
            let d = ((l[0] - p[0]).powi(2) + (l[1] - p[1]).powi(2) + (l[2] - p[2]).powi(2)).sqrt();
            if d < CONTACT_DISTANCE {
                expected += 1;
            }
        }
    }
    println!(
        "contacts at the crystal pose: {} (file: {expected})",
        f.contacts
    );
    assert_eq!(f.contacts, expected);
    assert!(expected > 0);
    let lining: Vec<[f64; 3]> = protein
        .iter()
        .copied()
        .filter(|p| {
            benzene.iter().any(|l| {
                ((l[0] - p[0]).powi(2) + (l[1] - p[1]).powi(2) + (l[2] - p[2]).powi(2)).sqrt()
                    < LINING_DISTANCE
            })
        })
        .collect();
    assert_eq!(lining.len(), c.lining().len());
    let (cl, cb) = (complex::centroid(&lining), complex::centroid(&benzene));
    let d = ((cl[0] - cb[0]).powi(2) + (cl[1] - cb[1]).powi(2) + (cl[2] - cb[2]).powi(2)).sqrt();
    println!(
        "cavity distance at the crystal pose: {:.4} Å over {} lining atoms",
        f.cavity_distance / ANGSTROM,
        lining.len()
    );
    // The same atoms summed in another order: a few ε of the coordinates.
    assert!((f.cavity_distance - d).abs() <= 64.0 * EPS * 30.0 * ANGSTROM);
}

/// **`⟨|u|²⟩ = 3B/(8π²)`, on a typed value.** Benzene's C1 in 181L has B = 20.05 Å² (column 61 of
/// its `HETATM` record); `3 × 20.05 / (8π²)` = 0.761 808 649 510 827 1 Å², worked to forty digits
/// with `π` typed to forty. And `B = 8π²/3` Å² is one Å² exactly, up to the rounding of `π²`.
#[test]
fn a_b_factor_is_three_over_eight_pi_squared_of_a_mean_square_displacement() {
    let a2 = ANGSTROM * ANGSTROM;
    let u2 = mean_square_displacement(20.05 * a2) / a2;
    // The forty-digit value to f64's sixteen: the seventeenth, 1, is far inside the bound.
    assert!((u2 - 0.761_808_649_510_827).abs() <= 4.0 * EPS, "{u2}");
    let pi = std::f64::consts::PI;
    assert!((mean_square_displacement(8.0 * pi * pi / 3.0) - 1.0).abs() <= 4.0 * EPS);
    let file = atoms_of_the_file(PDB_181L);
    assert_eq!(file[&('A', 400, "C1".to_string())].2, 20.05);
}

/// **The autocorrelation estimator, on a series whose correlation time is known.** An AR(1)
/// series `x_{t+1} = φ x_t + √(1 − φ²) ξ_t` has `ρ(t) = φ^t`, so `τ_int = ½ + Σ φ^t =
/// (1 + φ)/(2(1 − φ))`: 4.5 at φ = 0.8, and ½ for independent draws. Over N = 200 000 samples the
/// estimate's relative standard deviation is `√(2(2M + 1)/N)` with `M` the window, about 6τ
/// (Sokal, §3): 2.4% at φ = 0.8. The bound is four of those. The window's own bias, `e^(−M/τ)`,
/// is 0.25% and inside it. The mean, zero, is within four of the standard errors it reports.
#[test]
fn the_estimator_finds_a_known_correlation_time() {
    let n = 200_000;
    for (phi, seed) in [(0.0, 1u64), (0.8, 2)] {
        let mut rng = Rng::for_index(seed, 0);
        let mut x = 0.0;
        let xs: Vec<f64> = (0..n)
            .map(|_| {
                x = phi * x + (1.0f64 - phi * phi).sqrt() * rng.gaussian();
                x
            })
            .collect();
        let e = Estimate::of(&xs);
        let tau = (1.0 + phi) / (2.0 * (1.0 - phi));
        let window = 6.0 * tau;
        let spread = (2.0 * (2.0 * window + 1.0) / n as f64).sqrt();
        println!(
            "φ = {phi}: τ = {:.4} against {tau}, ±{:.4}; mean {:.5} ± {:.5}",
            e.tau,
            spread * tau,
            e.mean,
            e.error
        );
        assert!(
            (e.tau / tau - 1.0).abs() < 4.0 * spread,
            "τ {} against {tau}",
            e.tau
        );
        // The error itself against its closed form `√(2τ Var/N)`, Var = 1 by construction. It goes
        // as `√τ`, so its spread is half τ's, and the sample variance's own, `√(2/N)`, is far
        // smaller: inside the same allowance.
        let expected = (2.0 * tau / n as f64).sqrt();
        assert!(
            (e.error / expected - 1.0).abs() < 4.0 * spread,
            "error {} against {expected}",
            e.error
        );
        assert!(e.mean.abs() < 4.0 * e.error);
        assert_eq!(e.samples, n);
    }
    // **Over many seeds, the window's constant shows.** One series' τ is uncertain by 2.4%, so a
    // window stopped at `3τ` instead of `6τ`, biased by −4.6%, passes one series. The truncation
    // bias is `φ^M / ((1 − φ) τ)` for a window `M`: 0.27% at `M ≈ 6τ = 27` and 4.9% at `M ≈ 13.5`.
    // Over 32 seeds of 100 000 samples the mean of `τ̂/τ` has a standard error of 0.51%,
    // measured from the seeds' own spread, and the bound is four of them: 2.1%.
    let tau = (1.0 + 0.8) / (2.0 * (1.0 - 0.8));
    let ratios: Vec<f64> = (0..32u64)
        .map(|seed| {
            let mut rng = Rng::for_index(0xA1, seed);
            let mut x = 0.0;
            let xs: Vec<f64> = (0..100_000)
                .map(|_| {
                    x = 0.8 * x + 0.6 * rng.gaussian();
                    x
                })
                .collect();
            Estimate::of(&xs).tau / tau
        })
        .collect();
    let m = ratios.iter().sum::<f64>() / ratios.len() as f64;
    let se = (ratios.iter().map(|r| (r - m).powi(2)).sum::<f64>()
        / ((ratios.len() - 1) * ratios.len()) as f64)
        .sqrt();
    println!("32 seeds: ⟨τ̂/τ⟩ = {m:.5} ± {se:.5}");
    assert!((m - 1.0).abs() < 4.0 * se, "⟨τ̂/τ⟩ = {m} ± {se}");
    let flat = Estimate::of(&[2.5; 10]);
    assert_eq!((flat.mean, flat.error, flat.tau), (2.5, 0.0, 0.5));
}

// ---------------------------------------------------------------------------------------------
// The trajectories: measured, printed, not asserted
// ---------------------------------------------------------------------------------------------

/// A binding at `buffer` Å with each system's hydrogens relaxed (A-1), as a trajectory starts.
fn relaxed(system: &System, buffer: f64) -> Binding {
    let t = Instant::now();
    let b = Binding::new(system, buffer * ANGSTROM)
        .expect("the binding")
        .relaxing_hydrogens(20_000, Binding::HYDROGEN_TOLERANCE);
    let r = b.hydrogen_relaxation().unwrap();
    println!(
        "binding at {buffer} Å: {} atoms; hydrogens relaxed: complex {:?} in {} steps ({:.1} s)",
        b.positions().len(),
        r.complex.status,
        r.complex.steps,
        t.elapsed().as_secs_f64()
    );
    b
}

/// What one trajectory measured.
struct Run {
    complex: Complex,
    record: Record,
    ms_per_step: f64,
}

/// The mobile zone of `binding` at `mobile` Å minimised, thermalised at 300 K, `equilibrate` ps
/// discarded and `produce` ps recorded every 10 fs.
#[allow(clippy::too_many_arguments)]
fn trajectory(
    label: &str,
    binding: &Binding,
    mobile: f64,
    solvent: Solvent,
    equilibrate: f64,
    produce: f64,
    max_minimiser_steps: usize,
    seed: u64,
) -> Run {
    let mut c = Complex::new(binding, mobile * ANGSTROM, solvent).expect("the zone");
    let t = Instant::now();
    let p = c.minimise(max_minimiser_steps, Binding::HYDROGEN_TOLERANCE);
    let md0 = c.dynamics(Bath::Isolated, KELVIN, 0);
    let after = c.observe(&md0);
    println!(
        "{label}: {} mobile atoms in {} residues of {}; minimised {:?} in {} steps to {:.2e} \
         kcal/mol/Å ({:.0} s); ligand RMSD {:.3} Å, interaction {:.2} kcal/mol",
        c.mobile_count(),
        c.mobile_residues().len(),
        c.positions().len(),
        p.status,
        p.steps,
        p.max_force / KCAL_PER_MOL_ANGSTROM,
        t.elapsed().as_secs_f64(),
        after.rmsd / ANGSTROM,
        kcal(after.interaction())
    );
    let start = c.positions().to_vec();
    let mut md = c.dynamics(langevin(seed, FRICTION), KELVIN, seed ^ 0x5EED);
    let mut warm = Record::new(&c, 1 << 40);
    let per_ps = (1e-12 / DT).round() as usize;
    c.run(
        &mut md,
        DT,
        (equilibrate * per_ps as f64) as usize,
        &mut warm,
    );
    let mut record = Record::new(&c, 20);
    let steps = (produce * per_ps as f64) as usize;
    let t = Instant::now();
    c.run(&mut md, DT, steps, &mut record);
    let ms_per_step = t.elapsed().as_secs_f64() * 1e3 / steps as f64;
    // What must hold: the buffer kept its bits over the whole run.
    for (k, &m) in c.mobile().iter().enumerate() {
        if !m {
            assert_eq!(
                bits(&[c.positions()[k]]),
                bits(&[start[k]]),
                "frozen atom {k}"
            );
        }
    }
    Run {
        complex: c,
        record,
        ms_per_step,
    }
}

/// The statistics of a run, printed as a table row per observable.
fn report(label: &str, run: &Run) {
    let r = &run.record;
    let frames = r.frames();
    let ns_per_day = 86_400.0 / (run.ms_per_step * 1e-3) * DT * 1e9;
    println!(
        "{label}: {} frames over {:.0} ps, {:.2} ms/step, {:.2} ns/day",
        frames.len(),
        frames.len() as f64 * r.every() as f64 * DT * 1e12,
        run.ms_per_step,
        ns_per_day
    );
    println!("| observable | mean | ± | τ (frames) | min | max |");
    let row = |name: &str, xs: Vec<f64>| {
        let e = Estimate::of(&xs);
        let (lo, hi) = xs
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), x| {
                (l.min(*x), h.max(*x))
            });
        println!(
            "| {name} | {:.4} | {:.4} | {:.1} | {lo:.4} | {hi:.4} |",
            e.mean, e.error, e.tau
        );
    };
    row("ligand RMSD (Å)", r.series(|f| f.rmsd / ANGSTROM));
    row("ligand site RMSD (Å)", r.series(|f| f.site_rmsd / ANGSTROM));
    row(
        "centroid displacement (Å)",
        r.series(|f| f.centroid_displacement / ANGSTROM),
    );
    row(
        "distance to cavity centre (Å)",
        r.series(|f| f.cavity_distance / ANGSTROM),
    );
    row("contacts < 4 Å", r.series(|f| f.contacts as f64));
    row(
        "van der Waals (kcal/mol)",
        r.series(|f| kcal(f.van_der_waals)),
    );
    row("Coulomb (kcal/mol)", r.series(|f| kcal(f.electrostatic)));
    row(
        "interaction (kcal/mol)",
        r.series(|f| kcal(f.interaction())),
    );
    row("mobile temperature (K)", r.series(|f| f.temperature));
    let books = r.series(|f| kcal(f.books));
    let drift: Vec<f64> = books.iter().map(|b| b - books[0]).collect();
    row("books − first frame (kcal/mol)", drift.clone());
    // The in-plane orientation: which 60° sector the ring is in, how far inside it, and how often
    // it changes sector. A ring that jumps between its own sites sits near the sectors' centres;
    // one that turns freely fills them evenly, an RMS of 60°/√12 = 17.3° and a sixth of the frames
    // in each 10° bin.
    let sixty = std::f64::consts::PI / 3.0;
    let turns = r.series(|f| f.turn);
    let sector = |phi: f64| (phi / sixty).round();
    let inside: Vec<f64> = turns
        .iter()
        .map(|&p| (p - sector(p) * sixty).to_degrees())
        .collect();
    let changes = turns
        .windows(2)
        .filter(|w| (sector(w[0]) - sector(w[1])).rem_euclid(6.0) != 0.0)
        .count();
    let (mut net, mut path) = (0.0, 0.0);
    for w in turns.windows(2) {
        let d = (w[1] - w[0] + std::f64::consts::PI).rem_euclid(2.0 * std::f64::consts::PI)
            - std::f64::consts::PI;
        net += d;
        path += d.abs();
    }
    let mut bins = [0usize; 6];
    for d in &inside {
        bins[(((d + 30.0) / 10.0).floor() as usize).min(5)] += 1;
    }
    let rms = (inside.iter().map(|d| d * d).sum::<f64>() / inside.len() as f64).sqrt();
    let mut occupied = [0usize; 6];
    for &p in &turns {
        occupied[sector(p).rem_euclid(6.0) as usize] += 1;
    }
    println!(
        "in-plane turn: {changes} sector changes in {} frames; net {:+.0}°, path {:.0}°; within its \
         sector RMS {rms:.1}° (uniform 17.3°); 10° bins from −30°: {:?} of {}; frames per sector \
         (0°, 60°, …): {:?}",
        turns.len(),
        net.to_degrees(),
        path.to_degrees(),
        bins,
        inside.len(),
        occupied
    );
    let left = frames.iter().position(|f| f.centroid_displacement > LEAVES);
    println!(
        "leaves (centroid > {} Å): {}; books moved {:+.3} kcal/mol first frame to last, kT = {:.3}",
        LEAVES / ANGSTROM,
        match left {
            Some(i) => format!("at frame {i}"),
            None => "never".into(),
        },
        drift[drift.len() - 1],
        kcal(BOLTZMANN.to_si() * KELVIN)
    );
}

/// Pearson's r and its 95% interval by Fisher's transform.
fn pearson(x: &[f64], y: &[f64]) -> (f64, f64, f64) {
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = x.iter().map(|a| (a - mx).powi(2)).sum();
    let syy: f64 = y.iter().map(|b| (b - my).powi(2)).sum();
    let r = sxy / (sxx * syy).sqrt();
    let z = r.atanh();
    let w = 1.96 / (n - 3.0).sqrt();
    (r, (z - w).tanh(), (z + w).tanh())
}

/// The B-factor comparison: the ligand's fluctuation atom by atom against `3B/(8π²)`, and the
/// mobile protein's heavy atoms' as a correlation, by atom and by residue.
fn compare_b_factors(label: &str, system: &System, pdb: &str, run: &Run) {
    let file = atoms_of_the_file(pdb);
    let c = &run.complex;
    let b = c.binding();
    let a2 = ANGSTROM * ANGSTROM;
    let key = |k: usize| {
        let i = b.system_atoms()[k];
        let r = system.residue_of(i);
        (r.chain, r.number, system.atom_name(i).to_string())
    };
    println!("{label}: the ligand against its B-factors (Å²)");
    println!("| atom | B | 3B/(8π²) | simulated ⟨|r − ⟨r⟩|²⟩ | ratio |");
    let (mut sum_b, mut sum_md) = (0.0, 0.0);
    for &k in c.ligand_heavy_atoms() {
        let (chain, number, name) = key(k);
        let bf = file[&(chain, number, name.clone())].2;
        let u2 = mean_square_displacement(bf * a2) / a2;
        let msf = run.record.mean_square_fluctuation(k) / a2;
        sum_b += u2;
        sum_md += msf;
        println!(
            "| {name} | {bf:.2} | {u2:.3} | {msf:.3} | {:.2} |",
            msf / u2
        );
    }
    let n = c.ligand_heavy_atoms().len() as f64;
    let site = Estimate::of(&run.record.series(|f| (f.site_rmsd / ANGSTROM).powi(2)));
    // Beside it, the two numbers that say what the site RMSD can report at all, for this ligand's
    // crystal geometry: a planar ring turning freely reads `2a²(1 − 3/π)` whatever it does, with
    // `a` the crystal ring's mean radius; and an isotropic Gaussian of the crystal's own
    // `3B/(8π²)` about each site reads back only part of itself, because a large displacement is
    // nearer a neighbour's site. The second by Monte Carlo over 20 000 draws.
    let crystal: Vec<[f64; 3]> = c
        .ligand_heavy_atoms()
        .iter()
        .map(|&k| b.crystal_positions()[k])
        .collect();
    let kinds: Vec<Element> = c
        .ligand_heavy_atoms()
        .iter()
        .map(|&k| b.elements()[k])
        .collect();
    let centre = complex::centroid(&crystal);
    let radius = crystal
        .iter()
        .map(|p| {
            ((p[0] - centre[0]).powi(2) + (p[1] - centre[1]).powi(2) + (p[2] - centre[2]).powi(2))
                .sqrt()
        })
        .sum::<f64>()
        / crystal.len() as f64
        / ANGSTROM;
    let ceiling = 2.0 * radius * radius * (1.0 - 3.0 / std::f64::consts::PI);
    let u2 = sum_b / n;
    let width = (u2 / 3.0).sqrt() * ANGSTROM;
    let mut rng = Rng::for_index(0xB0, 0);
    let draws = 20_000;
    let mut read = 0.0;
    for _ in 0..draws {
        let moved: Vec<[f64; 3]> = crystal
            .iter()
            .map(|p| {
                [
                    p[0] + width * rng.gaussian(),
                    p[1] + width * rng.gaussian(),
                    p[2] + width * rng.gaussian(),
                ]
            })
            .collect();
        read += (complex::site_rmsd(&moved, &crystal, &kinds) / ANGSTROM).powi(2);
    }
    println!(
        "site RMSD² can report: a ring of radius {radius:.3} Å turning freely reads {ceiling:.3} \
         Å²; a Gaussian of ⟨|u|²⟩ = {u2:.3} Å² reads back {:.3} of itself",
        read / draws as f64 / u2
    );
    println!(
        "ligand ⟨site RMSD²⟩ (blind to which atom is on which site): {:.3} ± {:.3} Å² (τ {:.0} frames), ratio to 3B/(8π²) {:.2}",
        site.mean,
        site.error,
        site.tau,
        site.mean / (sum_b / n)
    );
    println!(
        "ligand mean: 3B/(8π²) {:.3} Å² (RMS {:.3} Å), simulated {:.3} Å² (RMS {:.3} Å), ratio {:.2}",
        sum_b / n,
        (sum_b / n).sqrt(),
        sum_md / n,
        (sum_md / n).sqrt(),
        sum_md / sum_b
    );
    let mut by_atom = (Vec::new(), Vec::new());
    let mut by_residue: Vec<(String, f64, f64, usize)> = Vec::new();
    for (span, residue) in b.residue_spans().iter().zip(b.residues()) {
        if !c.mobile()[span.start] {
            continue;
        }
        let (mut sb, mut sm, mut count) = (0.0, 0.0, 0);
        for k in span.clone() {
            if b.elements()[k] == Element::H {
                continue;
            }
            let bf = file[&key(k)].2;
            let u2 = mean_square_displacement(bf * a2) / a2;
            let msf = run.record.mean_square_fluctuation(k) / a2;
            by_atom.0.push(u2);
            by_atom.1.push(msf);
            sb += u2;
            sm += msf;
            count += 1;
        }
        by_residue.push((residue.clone(), sb / count as f64, sm / count as f64, count));
    }
    let (r, lo, hi) = pearson(&by_atom.0, &by_atom.1);
    let mean_b = by_atom.0.iter().sum::<f64>() / by_atom.0.len() as f64;
    let mean_md = by_atom.1.iter().sum::<f64>() / by_atom.1.len() as f64;
    println!(
        "mobile protein heavy atoms: {} atoms, mean 3B/(8π²) {mean_b:.3} Å², simulated {mean_md:.3} \
         Å²; r = {r:+.3}, 95% [{lo:+.3}, {hi:+.3}]",
        by_atom.0.len()
    );
    let x: Vec<f64> = by_residue.iter().map(|r| r.1).collect();
    let y: Vec<f64> = by_residue.iter().map(|r| r.2).collect();
    let (r, lo, hi) = pearson(&x, &y);
    println!(
        "by residue: {} residues, r = {r:+.3}, 95% [{lo:+.3}, {hi:+.3}]",
        by_residue.len()
    );
    for (name, u2, msf, count) in &by_residue {
        println!("  {name}: {count} heavy atoms, 3B/(8π²) {u2:.3}, simulated {msf:.3} Å²");
    }
}

/// **A cutoff's zone is the pocket's rule**: in a 10 Å binding, the 6 Å zone is exactly the
/// residues `Binding::new` cuts at 6 Å — 2c-1's eighteen, 322 atoms with benzene's — and at 3 Å
/// it is the four the 3 Å pocket has. Two bindings at different cutoffs, compared by their labels.
#[test]
#[ignore = "QEq on 987 atoms, about fifteen seconds with --release: run with --release -- --ignored"]
fn a_cutoffs_zone_is_the_pockets_rule() {
    let s = the_system();
    let b = Binding::new(s, 10.0 * ANGSTROM).unwrap();
    for (mobile, residues, atoms) in [(6.0, 18, 322), (3.0, 4, 82)] {
        let pocket = Binding::new(s, mobile * ANGSTROM).unwrap();
        let c = Complex::new(&b, mobile * ANGSTROM, Solvent::Vacuum).unwrap();
        assert_eq!(c.mobile_residues(), pocket.residues());
        assert_eq!(c.mobile_residues().len(), residues);
        assert_eq!(c.mobile_count(), atoms);
        assert_eq!(pocket.positions().len(), atoms);
    }
}

/// **Benzene in 181L at 300 K, in vacuum: 100 ps with a 6 Å mobile zone in a 10 Å binding.** The
/// main measurement of step 3b. Asserted: the zone is the 6 Å pocket of 2c-1 exactly (its 18
/// residues and 322 atoms, benzene's included), and the buffer keeps its bits. Everything else is
/// printed.
#[test]
#[ignore = "100 ps of dynamics, about six minutes with --release: run with --release -- --ignored --nocapture"]
fn benzene_stays_bound_measured() {
    let s = the_system();
    let b = relaxed(s, 10.0);
    let six = Binding::new(s, 6.0 * ANGSTROM).unwrap();
    let c = Complex::new(&b, 6.0 * ANGSTROM, Solvent::Vacuum).unwrap();
    assert_eq!(c.mobile_residues(), six.residues());
    assert_eq!(c.mobile_count(), six.positions().len());
    assert_eq!(c.mobile_count(), 322);
    let run = trajectory(
        "6 Å in 10 Å, vacuum",
        &b,
        6.0,
        Solvent::Vacuum,
        10.0,
        100.0,
        20_000,
        0x3B01,
    );
    report("6 Å in 10 Å, vacuum", &run);
    compare_b_factors("6 Å in 10 Å, vacuum", s, PDB_181L, &run);
}

/// **How much the zone matters**: a mobile radius of 8 Å (which needs a 14 Å binding: Trp126's
/// side chain reaches past 10 Å and its backbone would be cut), and the same 6 Å zone in that 14 Å
/// binding, so the buffer's effect is apart from the zone's. 50 ps each.
#[test]
#[ignore = "two runs of 50 ps, about twenty minutes with --release: run with --release -- --ignored --nocapture"]
fn the_zone_measured() {
    let s = the_system();
    let b = relaxed(s, 14.0);
    for mobile in [6.0, 8.0] {
        let label = format!("{mobile} Å in 14 Å, vacuum");
        let run = trajectory(
            &label,
            &b,
            mobile,
            Solvent::Vacuum,
            10.0,
            50.0,
            20_000,
            0x3B02,
        );
        report(&label, &run);
        compare_b_factors(&label, s, PDB_181L, &run);
    }
}

/// **Benzene in generalized Born**: the 6 Å zone in the 10 Å binding with OBC II over all 987
/// atoms. About sixty times the vacuum step's cost, so 2 ps of equilibration and 8 ps recorded:
/// what is affordable, and a tenth of the vacuum run's statistics.
#[test]
#[ignore = "10 ps of solvated dynamics, about forty minutes with --release: run with --release -- --ignored --nocapture"]
fn benzene_in_generalized_born_measured() {
    let s = the_system();
    let b = relaxed(s, 10.0);
    let label = "6 Å in 10 Å, OBC II";
    let run = trajectory(
        label,
        &b,
        6.0,
        Solvent::GeneralizedBorn,
        2.0,
        8.0,
        2000,
        0x3B03,
    );
    report(label, &run);
    compare_b_factors(label, s, PDB_181L, &run);
}

/// **A second ligand: n-butylbenzene in 186L**, the series' largest and best binder, in a 6 Å
/// zone as benzene's main run, 50 ps. Its 10 Å binding cuts Phe153's backbone, so the buffer is
/// the smallest of 12 and 14 Å that the cut allows.
#[test]
#[ignore = "50 ps of dynamics, about four minutes with --release: run with --release -- --ignored --nocapture"]
fn n_butylbenzene_measured() {
    const PDB_186L: &str = include_str!("../components/186L.pdb");
    let mut templates = templates();
    templates.push(Component::from_ccd(include_str!("../components/N4B.cif")).expect("N4B parses"));
    let s = System::from_pdb(
        PDB_186L,
        &templates,
        &Selection::new("N4B").dropping(&["HOH", "CL", "HED"]),
    )
    .expect("186L");
    // In a 10 Å binding the 6 Å zone is refused: Phe153 is bonded across the cut. The smallest buffer that holds it, of 12 and 14 Å, is the one run.
    let unrelaxed = Binding::new(&s, 10.0 * ANGSTROM).unwrap();
    println!(
        "10 Å: {:?}",
        Complex::new(&unrelaxed, 6.0 * ANGSTROM, Solvent::Vacuum).err()
    );
    let buffer = [12.0, 14.0]
        .into_iter()
        .find(|&r| {
            let b = Binding::new(&s, r * ANGSTROM).unwrap();
            let zone = Complex::new(&b, 6.0 * ANGSTROM, Solvent::Vacuum);
            println!("{r} Å: {:?}", zone.as_ref().err());
            zone.is_ok()
        })
        .expect("a buffer that holds the zone");
    let b = relaxed(&s, buffer);
    let label = format!("n-butylbenzene, 6 Å in {buffer} Å, vacuum");
    let label = label.as_str();
    let run = trajectory(label, &b, 6.0, Solvent::Vacuum, 10.0, 50.0, 20_000, 0x3B04);
    report(label, &run);
    compare_b_factors(label, &s, PDB_186L, &run);
}

/// **Equipartition over the mobile atoms, per element.** Benzene and the four residues within
/// 3 Å of it, free in an 8 Å binding — the smallest cutoff zone the cut allows — under a 20 ps⁻¹
/// bath so the kinetic energy decorrelates in tens of steps: every element's `⟨½mv²⟩` is
/// `(3/2) k_BT`, and the whole zone's kinetic energy is `(3 N_mobile / 2) k_BT`, each within four
/// standard errors from its own autocorrelation time. The degrees of freedom are `3 N_mobile`.
/// The frozen buffer contributes nothing: its atoms' kinetic energies are exactly zero.
///
/// **Ignored, and run by no gate.** The property it would catch first — a noise width that is not
/// each atom's own mass's — is held in the default run by `tests/the_dynamics.rs`'
/// `aspirin_shares_its_kinetic_energy_equally` (3a), not by this test. This one checks the same
/// property on the complex's zone, when it is run by hand.
#[test]
#[ignore = "30 000 steps on an 82-atom zone, about thirty seconds with --release: run with --release -- --ignored --nocapture"]
fn equipartition_holds_per_element_over_the_mobile_atoms() {
    let b = Binding::new(the_system(), 8.0 * ANGSTROM).unwrap();
    let mut c = Complex::new(&b, 3.0 * ANGSTROM, Solvent::Vacuum).unwrap();
    c.minimise(20_000, Binding::HYDROGEN_TOLERANCE);
    let kt = BOLTZMANN.to_si() * KELVIN;
    let mut md: MolecularDynamics = c.dynamics(langevin(0xE901, 20e12), KELVIN, 0xE902);
    let n_mobile = c.mobile_count();
    assert_eq!(md.degrees_of_freedom(), 3 * n_mobile);
    let mut warm = Record::new(&c, 1 << 40);
    c.run(&mut md, DT, 2000, &mut warm);
    let groups: Vec<(Element, Vec<usize>)> =
        [Element::H, Element::C, Element::N, Element::O, Element::S]
            .into_iter()
            .map(|e| {
                (
                    e,
                    (0..c.positions().len())
                        .filter(|&k| c.mobile()[k] && b.elements()[k] == e)
                        .collect::<Vec<_>>(),
                )
            })
            .filter(|(_, m)| !m.is_empty())
            .collect();
    let mut per: Vec<Vec<f64>> = vec![Vec::new(); groups.len()];
    let mut total = Vec::new();
    for _ in 0..30_000 {
        c.run(&mut md, DT, 1, &mut warm);
        let each = md.half_step_kinetic_energies();
        for ((_, members), out) in groups.iter().zip(per.iter_mut()) {
            let s: f64 = members.iter().map(|&k| each[k]).sum();
            out.push(s / (1.5 * kt * members.len() as f64));
        }
        let frozen: f64 = (0..each.len())
            .filter(|&k| !c.mobile()[k])
            .map(|k| each[k])
            .sum();
        assert_eq!(frozen, 0.0);
        total.push(md.half_step_kinetic_energy() / (1.5 * n_mobile as f64 * kt));
    }
    for ((e, members), out) in groups.iter().zip(&per) {
        let est = Estimate::of(out);
        let z = (est.mean - 1.0) / est.error;
        println!(
            "{} ({}): ⟨½mv²⟩/(3/2 k_BT) = {:.4} ± {:.4} (τ = {:.1}), z = {z:+.2}",
            e.symbol(),
            members.len(),
            est.mean,
            est.error,
            est.tau
        );
        assert!(z.abs() < 4.0, "{}: {est:?}", e.symbol());
    }
    let est = Estimate::of(&total);
    let z = (est.mean - 1.0) / est.error;
    println!(
        "{n_mobile} mobile atoms: ⟨KE⟩/(3N/2 k_BT) = {:.4} ± {:.4} (τ = {:.1}), z = {z:+.2}",
        est.mean, est.error, est.tau
    );
    assert!(z.abs() < 4.0, "{est:?}");
}
