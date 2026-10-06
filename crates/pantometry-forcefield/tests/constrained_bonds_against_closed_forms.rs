//! **Bonds held by SHAKE and RATTLE, against exact identities**: benzene's six C–H bonds stay at
//! UFF's natural length to rounding and every velocity is tangent to them, beside rigid TIP3P
//! water held by SETTLE; bonds that share atoms converge to the stated tolerance; the degrees of
//! freedom are counted for both kinds of constraint and for frozen atoms, and a bath fills exactly
//! that many with `k_BT/2` each; constrained NVE's energy error falls as `h²` at the step the
//! constraints were added for; a frozen atom holds its partner at the bond's length; and a run is
//! the same bits however it is cut. Every test here runs unoptimised.

// The sums over atoms and axes are written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

mod hydration;

use hydration::*;
use pantometry_forcefield::uff::KCAL_PER_MOL;
use pantometry_forcefield::water::Settle;
use pantometry_forcefield::{Bath, Estimate, MolecularDynamics, Potential, Shake};
use pantometry_units::BOLTZMANN;

const EPS: f64 = f64::EPSILON;

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn len(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn largest_coordinate(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().chain(&b).fold(0.0f64, |m, x| m.max(x.abs()))
}

/// A held length may be off by this many ε of the larger coordinate of its two atoms: each
/// position is stored to half an ulp of itself, so a bond far from the origin cannot be closer to
/// its length than that, whatever solved it. W2's bound for SETTLE, the same reason.
const ROUNDING: f64 = 16.0;

/// SHAKE's tolerance as the module states it, written here rather than read from
/// [`Shake::tolerance`], so that a looser default fails the bounds built on it.
const TOL: f64 = 1e-12;

/// The bound on a held bond's relative velocity along itself, `|r·(v_i − v_j)|`: rounding, or
/// RATTLE's own tolerance, whichever is larger. A bond already within `TOL` of tangent is not
/// corrected, so tolerance is all the library promises; bonds that share no atom are in practice
/// solved in one exact sweep, and measured at rounding (the worst is printed).
fn tangent_bound(r: [f64; 3], vi: [f64; 3], vj: [f64; 3]) -> f64 {
    let rounding = ROUNDING * EPS * len(r) * (len(vi) + len(vj));
    rounding.max(1.01 * TOL * len(r) * len(sub(vi, vj)))
}

/// Benzene in 27 lattice waters less those it overlaps, 9.3 Å across, `r_c` = 4.6 Å, δ = 10⁻⁵.
fn small_box(seed: u64) -> (Benzene, Solvated) {
    let b = Benzene::new(None);
    let s = Solvated::new(&b, 3, seed, 2.3 * ANGSTROM, 4.6 * ANGSTROM, 1e-5);
    (b, s)
}

/// The UFF natural lengths of benzene's C–H bonds are what `to_hydrogen` holds, the six are the
/// bonds to a hydrogen and share no atom — one sweep solves them — and the relaxed start is on
/// them.
#[test]
fn benzenes_bonds_to_hydrogen_are_six_at_ufs_natural_length() {
    let b = Benzene::new(None);
    let s = b.shake();
    assert_eq!(s.len(), 6);
    assert!(s.is_independent());
    for (k, &[i, j]) in s.bonds().iter().enumerate() {
        let elements = [i, j].map(|a| b.component.atoms()[a].element);
        assert!(elements.contains(&pantometry_forcefield::Element::H));
        let stretch = b
            .vacuum
            .stretches()
            .iter()
            .find(|t| t.atoms == [i, j])
            .expect("a stretch for the bond");
        assert_eq!(s.lengths()[k], stretch.natural_length);
        let r = len(sub(b.at[i], b.at[j]));
        assert!(
            (r - s.lengths()[k]).abs() <= ROUNDING * EPS * largest_coordinate(b.at[i], b.at[j])
        );
    }
    println!("C–H held at {:.6} Å", s.lengths()[0] / ANGSTROM);
}

/// Under a bath at 2 fs, with SETTLE on every water and SHAKE on benzene's C–H: after every step
/// each C–H is at its length and each water rigid, to rounding; the whole-step velocities have
/// nothing along a C–H bond; and the velocities after `O` have nothing along one at the positions
/// they were projected at.
#[test]
fn held_bonds_stay_at_their_length_to_rounding_and_their_velocities_are_tangent() {
    let (b, s) = small_box(0x5A1);
    let field = &s.field;
    let mut at = s.at.clone();
    let mut md = s
        .dynamics(&b)
        .with_bath(Bath::Langevin {
            temperature: 298.15,
            friction: 5e12,
            seed: 3,
        })
        .thermalised(&at, 298.15, 3);
    let shake = b.shake();
    let settle = Settle::tip3p(field.rigid_waters().to_vec());
    let lengths = settle.lengths();
    let tangent = |at: &[[f64; 3]], v: &[[f64; 3]], [i, j]: [usize; 2]| {
        let r = sub(at[i], at[j]);
        let u = sub(v[i], v[j]);
        (dot(r, u).abs(), tangent_bound(r, v[i], v[j]))
    };
    let (mut worst_length, mut worst_velocity) = (0.0f64, 0.0f64);
    for _ in 0..60 {
        md.step(field, &mut at, 2.0 * FS);
        for (k, &[i, j]) in shake.bonds().iter().enumerate() {
            let off = (len(sub(at[i], at[j])) - shake.lengths()[k]).abs();
            let bound = ROUNDING * EPS * largest_coordinate(at[i], at[j]);
            worst_length = worst_length.max(off / (EPS * largest_coordinate(at[i], at[j])));
            assert!(off <= bound, "C–H {k}: {off:e} m off");
            let (along, bound) = tangent(&at, md.velocities(), [i, j]);
            let v = md.velocities();
            let r = len(sub(at[i], at[j]));
            worst_velocity = worst_velocity.max(along / (EPS * r * (len(v[i]) + len(v[j]))));
            assert!(
                along <= bound,
                "C–H {k}: whole-step velocity along the bond"
            );
            let (along, bound) =
                tangent(md.half_step_positions(), md.half_step_velocities(), [i, j]);
            assert!(along <= bound, "C–H {k}: velocity after O along the bond");
        }
        for &[o, h1, h2] in field.rigid_waters() {
            for (pair, d) in [[o, h1], [o, h2], [h1, h2]].iter().zip(lengths) {
                let off = (len(sub(at[pair[0]], at[pair[1]])) - d).abs();
                assert!(off <= ROUNDING * EPS * largest_coordinate(at[pair[0]], at[pair[1]]));
            }
        }
    }
    println!(
        "worst C–H length {worst_length:.2} ε × its coordinate, worst velocity along a bond \
         {worst_velocity:.2} ε (bound {ROUNDING})"
    );
}

/// No potential: free particles, so that each C–H pair is a free rigid dumbbell.
struct Nothing;

impl Potential for Nothing {
    fn energy_and_forces(&self, _: &[[f64; 3]], f: &mut [[f64; 3]]) -> f64 {
        f.iter_mut().for_each(|x| *x = [0.0; 3]);
        0.0
    }
}

/// Every C–C bond of the ring held as well as the C–H: twelve constraints, each carbon in three,
/// so the sweeps iterate. Every length within the tolerance of its own after every step of NVE in
/// vacuum, and every velocity within it of tangent; and the count is `3N − 12`, less the six
/// rigid motions NVE keeps removed.
#[test]
fn bonds_that_share_atoms_converge_to_their_tolerance() {
    let b = Benzene::new(None);
    let mut bonds = Vec::new();
    let mut lengths = Vec::new();
    for s in b.vacuum.stretches() {
        bonds.push(s.atoms);
        lengths.push(len(sub(b.at[s.atoms[0]], b.at[s.atoms[1]])));
    }
    let shake = Shake::new(bonds, lengths);
    assert!(!shake.is_independent());
    assert_eq!(shake.len(), 12);
    assert_eq!(shake.tolerance(), TOL, "the tolerance the module states");
    let tol = TOL;
    let mut at = b.at.clone();
    let mut md = MolecularDynamics::new(b.masses.clone())
        .with_bond_constraints(shake.clone())
        .thermalised(&at, 298.15, 9);
    assert_eq!(md.degrees_of_freedom(), 3 * 12 - 12 - 6);
    for _ in 0..100 {
        md.step(&b.vacuum, &mut at, 2.0 * FS);
        for (k, &[i, j]) in shake.bonds().iter().enumerate() {
            let d = shake.lengths()[k];
            let r = len(sub(at[i], at[j]));
            // |r² − d²| ≤ 2 tol d² puts r within tol d of d to first order; rounding on top.
            let bound = 1.01 * tol * d + ROUNDING * EPS * largest_coordinate(at[i], at[j]);
            assert!(
                (r - d).abs() <= bound,
                "bond {k}: {:e} of its length",
                (r - d) / d
            );
            let u = sub(md.velocities()[i], md.velocities()[j]);
            let along = dot(sub(at[i], at[j]), u).abs();
            let v = md.velocities();
            assert!(
                along <= 1.01 * tol * r * len(u) + ROUNDING * EPS * r * (len(v[i]) + len(v[j])),
                "bond {k}: velocity along it"
            );
        }
    }
}

/// **Shared atoms, and the energy**: the twelve-bond benzene of
/// [`bonds_that_share_atoms_converge_to_their_tolerance`] in NVE in vacuum, 0.2 ps from the same
/// velocities at 1 fs and 0.5 fs. The RMS energy error falls by four, held to
/// [`H_SQUARED_FROM_ONE`]: the velocity correction `Δq/h` of an atom in three bonds is its one
/// displacement, once, which nothing about the lengths or the tangency would show.
#[test]
fn bonds_that_share_atoms_conserve_energy_as_h_squared() {
    let b = Benzene::new(None);
    let (mut bonds, mut lengths) = (Vec::new(), Vec::new());
    for s in b.vacuum.stretches() {
        bonds.push(s.atoms);
        lengths.push(len(sub(b.at[s.atoms[0]], b.at[s.atoms[1]])));
    }
    let shake = Shake::new(bonds, lengths);
    let rms = |dt: f64, steps: usize| {
        let mut at = b.at.clone();
        let mut md = MolecularDynamics::new(b.masses.clone())
            .with_bond_constraints(shake.clone())
            .thermalised(&at, 298.15, 4);
        md.prepare(&b.vacuum, &at);
        let e0 = md.kinetic_energy() + md.potential_energy().unwrap();
        let mut squares = 0.0;
        for _ in 0..steps {
            md.step(&b.vacuum, &mut at, dt);
            let e = md.kinetic_energy() + md.potential_energy().unwrap() - e0;
            squares += e * e;
        }
        (squares / steps as f64).sqrt()
    };
    let (e1, e05) = (rms(1.0 * FS, 200), rms(0.5 * FS, 400));
    println!(
        "twelve held bonds: RMS ΔE {:.3e} and {:.3e} kcal/mol at 1 and 0.5 fs, ratio {:.4}",
        e1 / KCAL_PER_MOL,
        e05 / KCAL_PER_MOL,
        e1 / e05
    );
    assert!(
        H_SQUARED_FROM_ONE.contains(&(e1 / e05)),
        "ratio {}",
        e1 / e05
    );
}

/// **Held bonds alone under a bath**, no rigid water: benzene as six free C–H dumbbells under no
/// potential, `γh = 1`. The velocities after `O` have nothing along a bond at the positions they
/// were projected at, and the kinetic energy there is `15 k_BT`, the thirty counted degrees of
/// freedom's, within four standard errors — an `O` left unprojected would put `18 k_BT` there.
#[test]
fn held_bonds_alone_take_the_baths_temperature() {
    let b = Benzene::new(None);
    let dt = 1.0 * FS;
    let mut at = b.at.clone();
    let mut md = MolecularDynamics::new(b.masses.clone())
        .with_bond_constraints(b.shake())
        .with_bath(Bath::Langevin {
            temperature: 298.15,
            friction: 1.0 / dt,
            seed: 31,
        })
        .thermalised(&at, 298.15, 31);
    assert_eq!(md.degrees_of_freedom(), 30);
    let kt = BOLTZMANN.to_si() * 298.15;
    let mut share = Vec::new();
    for _ in 0..8000 {
        md.step(&Nothing, &mut at, dt);
        share.push(md.half_step_kinetic_energy() / (15.0 * kt));
        let (x, v) = (md.half_step_positions(), md.half_step_velocities());
        for &[i, j] in b.shake().bonds() {
            let r = sub(x[i], x[j]);
            let along = dot(r, sub(v[i], v[j])).abs();
            assert!(along <= tangent_bound(r, v[i], v[j]));
        }
    }
    let e = Estimate::of(&share);
    let z = (e.mean - 1.0) / e.error;
    println!(
        "held bonds alone: {:.5} ± {:.5} of 15 k_BT, z = {z:+.2}",
        e.mean, e.error
    );
    assert!(z.abs() < 4.0);
}

/// The degrees of freedom: three a free atom, less three a rigid water, less one a held bond with
/// an atom free to move, less — in NVE — the rigid motion `thermalised` removed. A bond to a
/// frozen atom still holds its partner, so it still costs one; a bond frozen at both ends costs
/// none.
#[test]
fn the_degrees_of_freedom_count_both_kinds_of_constraint() {
    let (b, s) = small_box(0xD0F);
    let n = s.at.len();
    let w = s.waters;
    let bath = Bath::Langevin {
        temperature: 298.15,
        friction: 1e12,
        seed: 1,
    };
    let md = s.dynamics(&b).with_bath(bath);
    assert_eq!(md.degrees_of_freedom(), 3 * n - 3 * w - 6);
    let nve = s.dynamics(&b).thermalised(&s.at, 298.15, 1);
    assert_eq!(nve.degrees_of_freedom(), 3 * n - 3 * w - 6 - 6);
    let [c, h] = b.shake().bonds()[0];
    let mut frozen = vec![false; n];
    frozen[h] = true;
    let one = s.dynamics(&b).with_bath(bath).with_frozen(frozen.clone());
    assert_eq!(one.degrees_of_freedom(), 3 * (n - 1) - 3 * w - 6);
    frozen[c] = true;
    let both = s.dynamics(&b).with_bath(bath).with_frozen(frozen);
    assert_eq!(both.degrees_of_freedom(), 3 * (n - 2) - 3 * w - 5);
}

/// **Equipartition over exactly the counted degrees of freedom**: no potential, so benzene is six
/// free C–H dumbbells and the waters free rigid bodies, in a bath with `γh = 1` so that samples
/// are nearly independent. Each dumbbell's centre of mass carries `(3/2) k_BT` and its rotation
/// `k_BT` — two degrees of freedom, the third taken by the bond — each water `3 k_BT`, and the
/// temperature read from all of it through `degrees_of_freedom` is the bath's. Each within four
/// standard errors.
#[test]
fn a_bath_fills_exactly_the_counted_degrees_of_freedom() {
    let (b, s) = small_box(0xE9);
    let dt = 1.0 * FS;
    let mut at = s.at.clone();
    let mut md = s
        .dynamics(&b)
        .with_bath(Bath::Langevin {
            temperature: 298.15,
            friction: 1.0 / dt,
            seed: 13,
        })
        .thermalised(&at, 298.15, 13);
    let kt = BOLTZMANN.to_si() * 298.15;
    let m = &s.masses;
    let bonds = b.shake().bonds().to_vec();
    let (mut trans, mut rot, mut water, mut temperature) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for _ in 0..12_000 {
        md.step(&Nothing, &mut at, dt);
        let v = md.half_step_velocities();
        let ke = |i: usize| 0.5 * m[i] * dot(v[i], v[i]);
        let (mut t, mut r) = (0.0, 0.0);
        for &[i, j] in &bonds {
            let mass = m[i] + m[j];
            let com = [0, 1, 2].map(|k| (m[i] * v[i][k] + m[j] * v[j][k]) / mass);
            let k_t = 0.5 * mass * dot(com, com);
            t += k_t;
            r += ke(i) + ke(j) - k_t;
        }
        trans.push(t / (6.0 * 1.5 * kt));
        rot.push(r / (6.0 * kt));
        let w: f64 = (12..at.len()).map(ke).sum();
        water.push(w / (s.waters as f64 * 3.0 * kt));
        temperature.push(md.half_step_temperature() / 298.15);
    }
    for (name, series) in [
        ("dumbbell translation", &trans),
        ("dumbbell rotation", &rot),
        ("water", &water),
        ("temperature", &temperature),
    ] {
        let e = Estimate::of(series);
        let z = (e.mean - 1.0) / e.error;
        println!(
            "{name}: {:.5} ± {:.5} of its share, z = {z:+.2}, τ = {:.2}",
            e.mean, e.error, e.tau
        );
        assert!(z.abs() < 4.0, "{name}");
    }
}

/// A bond with one atom frozen holds the other at its length from it: the free hydrogen moves, the
/// frozen carbon keeps its bits, and the distance is the bond's to rounding after every step.
#[test]
fn a_frozen_atom_holds_its_partner_at_the_bonds_length() {
    let b = Benzene::new(None);
    let n = b.at.len();
    let shake = b.shake();
    let [c, h] = shake.bonds()[0];
    let d = shake.lengths()[0];
    let mut frozen = vec![false; n];
    frozen[c] = true;
    let mut at = b.at.clone();
    let start = at.clone();
    let mut md = MolecularDynamics::new(b.masses.clone())
        .with_bond_constraints(shake)
        .with_frozen(frozen)
        .with_bath(Bath::Langevin {
            temperature: 298.15,
            friction: 5e12,
            seed: 2,
        })
        .thermalised(&at, 298.15, 2);
    for _ in 0..100 {
        md.step(&b.vacuum, &mut at, 2.0 * FS);
        assert_eq!(at[c], start[c]);
        let r = len(sub(at[c], at[h]));
        assert!((r - d).abs() <= ROUNDING * EPS * largest_coordinate(at[c], at[h]));
    }
    assert!(
        len(sub(at[h], start[h])) > 1e-3 * ANGSTROM,
        "the hydrogen moved"
    );
}

/// An atom in a rigid water and in a held bond is refused, in either order of building.
#[test]
#[should_panic(expected = "in a rigid molecule and in a constrained bond")]
fn an_atom_held_twice_is_refused() {
    let (b, s) = small_box(0x2);
    let o = s.field.rigid_waters()[0][0];
    let _ = s
        .dynamics(&b)
        .with_bond_constraints(Shake::new(vec![[0, o]], vec![3.0 * ANGSTROM]));
}

/// A run with SETTLE and SHAKE together under a bath is the same bits in one call or cut into
/// pieces.
#[test]
fn a_constrained_run_is_the_same_bits_however_it_is_cut() {
    let (b, s) = small_box(0xC07);
    let field = &s.field;
    let md = s
        .dynamics(&b)
        .with_bath(Bath::Langevin {
            temperature: 298.15,
            friction: 5e12,
            seed: 21,
        })
        .thermalised(&s.at, 298.15, 21);
    let mut whole = (md.clone(), s.at.clone());
    whole.0.run(field, &mut whole.1, 2.0 * FS, 12);
    for cuts in [vec![1usize; 12], vec![5, 7], vec![3, 1, 8]] {
        let mut part = (md.clone(), s.at.clone());
        for c in cuts {
            part.0.run(field, &mut part.1, 2.0 * FS, c);
        }
        assert_eq!(part.1, whole.1);
        assert_eq!(part.0.velocities(), whole.0.velocities());
        assert_eq!(part.0.steps(), 12);
    }
}

/// Benzene and the waters nearest it, alone in a box 30 Å across with `r_c` = 14 Å: no pair of the
/// cluster reaches the cutoff and no image comes inside it, so the potential is smooth and the
/// energy error of NVE is the integrator's.
fn cluster(seed: u64) -> (Benzene, Solvated) {
    let b = Benzene::new(None);
    let near = Solvated::new(&b, 3, seed, 2.3 * ANGSTROM, 4.6 * ANGSTROM, 1e-5);
    // The waters whose oxygen is within 5 Å of benzene's centroid, moved into the big box.
    let side = 30.0 * ANGSTROM;
    let centre = near.cell().lengths()[0] * 0.5;
    let shift = 0.5 * side - centre;
    let mut kept = Vec::new();
    for &[o, h1, h2] in near.field.rigid_waters() {
        let d = sub(near.at[o], [centre; 3]);
        if len(d) < 5.0 * ANGSTROM {
            kept.push([o, h1, h2]);
        }
    }
    let big = pantometry_forcefield::PeriodicBox::cubic(side);
    let p = pantometry_forcefield::EwaldParameters::for_accuracy(&big, 14.0 * ANGSTROM, 1e-5);
    let field = pantometry_forcefield::PeriodicForceField::solvated(
        &b.component,
        &b.types,
        kept.len(),
        big,
        p,
    )
    .expect("the cluster")
    .with_solute_charges(b.charges.clone());
    let moved = |p: [f64; 3]| p.map(|x| x + shift);
    let mut at: Vec<[f64; 3]> = near.at[..12].iter().map(|&p| moved(p)).collect();
    let mut masses = b.masses.clone();
    for w in &kept {
        for &i in w {
            at.push(moved(near.at[i]));
            masses.push(near.masses[i]);
        }
    }
    let mut group = vec![false; at.len()];
    group[..12].fill(true);
    let waters = kept.len();
    (
        b,
        Solvated {
            field,
            at,
            masses,
            waters,
            group,
        },
    )
}

/// The RMS departure of `kinetic + potential` from its start over `steps` of `dt`, NVE.
fn nve_rms(b: &Benzene, s: &Solvated, dt: f64, steps: usize, seed: u64) -> f64 {
    let mut at = s.at.clone();
    let mut md = s.dynamics(b).thermalised(&at, 298.15, seed);
    md.prepare(&s.field, &at);
    let e0 = md.kinetic_energy() + md.potential_energy().expect("prepared");
    let mut squares = 0.0;
    for _ in 0..steps {
        md.step(&s.field, &mut at, dt);
        let e = md.kinetic_energy() + md.potential_energy().expect("stepped");
        squares += (e - e0) * (e - e0);
    }
    (squares / steps as f64).sqrt()
}

/// The three ratios of the RMS energy error per halving, 2 → 1, 1 → 0.5 and 0.5 → 0.25 fs, each
/// run over 0.2 ps from the same start and velocities.
fn ratios(b: &Benzene, s: &Solvated, seed: u64) -> [f64; 3] {
    let e = [(2.0, 100), (1.0, 200), (0.5, 400), (0.25, 800)]
        .map(|(h, n)| nve_rms(b, s, h * FS, n, seed));
    [e[0] / e[1], e[1] / e[2], e[2] / e[3]]
}

/// The bands the RMS energy error's ratio per halving is held to, from 2 fs to 1 fs and from 1 fs
/// to 0.5 fs. **Measured** over sixteen starts ([`the_h_squared_band_measured`], release): 4.229 to
/// 4.508 from 2 fs, and 4.050 to 4.109 from 1 fs. Above 4, and more so from 2 fs, because the next
/// term is `h⁴` and positive: `RMS ∝ h² (1 + c h²)` with the second ratio's `c` ≈ 0.02 fs⁻² predicts
/// 4.24 for the first. The bands are 3.2 and 6.8 times the measured spreads. **They catch an error
/// of order `h` or `h³` with a large coefficient, which would read 2 or 8, and not one with a small
/// coefficient**: `Δq/h` scaled by 0.99 read 4.375 and 4.177, inside both. [`H_SQUARED_SHRINKS`]
/// is the check that sees that.
const H_SQUARED_FROM_TWO: std::ops::Range<f64> = 3.9..4.8;
const H_SQUARED_FROM_ONE: std::ops::Range<f64> = 3.8..4.2;

/// **The excess over 4 shrinks by four per halving**: `(r₃ − 4)/(r₂ − 4)`, with `r₂` the ratio from
/// 1 fs and `r₃` the ratio from 0.5 fs to 0.25 fs. For `RMS ∝ h² (1 + c h²)` the excess is about
/// `3c h²`, a quarter of itself one halving on. Measured over the sixteen starts, 0.2365 to 0.2593.
/// A defect of lower order in `h` with a small coefficient goes the other way — negligible at 2 fs,
/// dominant at 0.25 — so the excess grows: `Δq/h` scaled by 0.99 read 0.491 to 1.322 over the same
/// sixteen. The band is 0.25 ± 0.10: the prediction, nine times the measured half-spread, and short
/// of the sabotage's least by 0.14. [`H_SQUARED_FROM_HALF`] holds `r₃` itself: measured 4.012 to
/// 4.027, and 4.029 to 4.279 under the sabotage.
const H_SQUARED_SHRINKS: std::ops::Range<f64> = 0.15..0.35;
const H_SQUARED_FROM_HALF: std::ops::Range<f64> = 3.97..4.05;

/// **Constrained velocity Verlet at the step the constraints were added for**: benzene with its
/// C–H held and the ten waters nearest it, rigid, NVE over 0.2 ps from 298 K velocities. The RMS
/// energy error falls by four per halving from 2 fs, held to [`H_SQUARED_FROM_TWO`],
/// [`H_SQUARED_FROM_ONE`] and [`H_SQUARED_FROM_HALF`], and its excess over four shrinks by four
/// per halving, held to [`H_SQUARED_SHRINKS`]; at 2 fs it is printed in kcal/mol.
///
/// **This test is the whole guard of SHAKE's velocity correction `Δq/h`** for bonds that share no
/// atom. Removing it entirely passes every Langevin and equipartition test here and everything in
/// `benzene_hydrated_in_tip3p.rs`: the projection after the next kick takes out what it would have,
/// and only an energy that is meant to be conserved shows the difference.
#[test]
fn constrained_nve_at_two_femtoseconds_falls_as_h_squared() {
    let (b, s) = cluster(0xA2);
    let e2 = nve_rms(&b, &s, 2.0 * FS, 100, 5);
    let r = ratios(&b, &s, 5);
    let shrink = (r[2] - 4.0) / (r[1] - 4.0);
    println!(
        "benzene and {} waters, RMS ΔE over 0.2 ps at 2 fs {:.3e} kcal/mol; ratios {:.4}, {:.4}, \
         {:.4}; the excess shrinks by {shrink:.4}",
        s.waters,
        e2 / KCAL_PER_MOL,
        r[0],
        r[1],
        r[2]
    );
    assert!(H_SQUARED_FROM_TWO.contains(&r[0]), "ratio {}", r[0]);
    assert!(H_SQUARED_FROM_ONE.contains(&r[1]), "ratio {}", r[1]);
    assert!(H_SQUARED_FROM_HALF.contains(&r[2]), "ratio {}", r[2]);
    assert!(
        H_SQUARED_SHRINKS.contains(&shrink),
        "the excess shrinks by {shrink}"
    );
}

/// The measurement the `H_SQUARED` bands come from: sixteen starts, each a lattice turned by its
/// own seed and velocities from their own, the three ratios of each and how the excess shrinks. Release:
/// `cargo test --release -p pantometry-forcefield --test constrained_bonds_against_closed_forms
/// the_h_squared_band_measured -- --ignored --nocapture`.
#[test]
#[ignore = "sixteen starts of four NVE runs: run with --release -- --ignored --nocapture"]
fn the_h_squared_band_measured() {
    let mut lo = [f64::INFINITY; 4];
    let mut hi = [f64::NEG_INFINITY; 4];
    for start in 0..16u64 {
        let (b, s) = cluster(0xA2 + start);
        let r = ratios(&b, &s, 5 + start);
        let shrink = (r[2] - 4.0) / (r[1] - 4.0);
        println!(
            "start {start}: ratios {:.4} {:.4} {:.4}, (r₃ − 4)/(r₂ − 4) {shrink:.4}",
            r[0], r[1], r[2]
        );
        for (k, x) in [r[0], r[1], r[2], shrink].into_iter().enumerate() {
            lo[k] = lo[k].min(x);
            hi[k] = hi[k].max(x);
        }
    }
    println!("from 2 fs {:.4}–{:.4}, from 1 fs {:.4}–{:.4}, from 0.5 fs {:.4}–{:.4}, excess shrink {:.4}–{:.4}",
        lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], lo[3], hi[3]);
}
