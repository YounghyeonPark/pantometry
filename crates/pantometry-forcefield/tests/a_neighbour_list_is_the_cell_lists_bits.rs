//! **A Verlet neighbour list changes no bit**: [`PeriodicForceField::with_neighbour_list`] against
//! the cell list it replaces, along molecular-dynamics trajectories that cross many builds, under
//! every kind of constraint the dynamics has.
//!
//! The identity is legitimate here, not a comparison of two implementations of one idea: both
//! enumerate the same set of pairs — every pair inside the cutoff — and the list reproduces the
//! cell list's order, so the sums are the same additions. What could go wrong is the set, and the
//! unit tests in `src/neighbours.rs` hold it to a brute-force search over 27 images and to the
//! constructed cases at the rule's edge.
//!
//! Default, unoptimised, about 18 s — 40 to 60 steps each, the skins small enough for a build
//! every three to six:
//!
//! - **W1**: UFF's flexible water, 64 molecules, no constraint, every molecule cut by the faces;
//! - **W2**: rigid TIP3P, 125 waters, under SETTLE;
//! - **W3**: benzene in TIP3P under SETTLE and SHAKE, through [`PeriodicDecoupling`] at a state
//!   between the ends, in a box too small for cells and one large enough, the couplings too;
//! - **W4**: 181L's 3 Å pocket and benzene, through the mesh, the protein frozen.
//!
//! In each, every step's positions, velocities and potential energy are the same bits with and
//! without the list; every tenth step's whole evaluation is — every energy term, every force and
//! the virial — and every twentieth the decoupling's, at five states, with its couplings; and the
//! number of builds is the one an independent replay of the rule predicts from the trajectory,
//! frozen atoms included. Ignored: the same at hundreds of steps and larger skins, with 216
//! waters and the 5 Å pocket and its chloride (271 s unoptimised); and, in release, the
//! cost on W4's box.

// The sums over atoms and axes are written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

mod common;
mod hydration;
mod protein;

use hydration::{Benzene, Solvated};
use pantometry_forcefield::uff;
use pantometry_forcefield::water::WaterBox;
use pantometry_forcefield::{
    Alchemical, AtLambda, Bath, Binding, Component, EwaldParameters, Flexibility, Lambda,
    MolecularDynamics, PeriodicBox, PeriodicDecoupling, PeriodicEnergy, PeriodicEvaluation,
    PeriodicForceField, Potential, SolvatedComplex, Solvation,
};

const ANGSTROM: f64 = 1e-10;
const FS: f64 = 1e-15;

/// Every energy term as bits.
fn energy_bits(e: &PeriodicEnergy) -> [u64; 14] {
    let w = &e.ewald;
    [
        e.bond,
        e.angle,
        e.torsion,
        e.inversion,
        e.van_der_waals,
        e.dispersion_correction,
        e.electrostatic,
        e.total,
        w.real,
        w.reciprocal,
        w.self_energy,
        w.excluded,
        w.background,
        w.total,
    ]
    .map(f64::to_bits)
}

/// A trajectory's end, positions and velocities, as bits.
type EndBits = (Vec<[u64; 3]>, Vec<[u64; 3]>);

fn bits(v: &[[f64; 3]]) -> Vec<[u64; 3]> {
    v.iter().map(|p| p.map(f64::to_bits)).collect()
}

/// Two evaluations are the same bits: every term, every force, the virial.
fn same_evaluation(a: &PeriodicEvaluation, b: &PeriodicEvaluation, what: &str) {
    assert_eq!(energy_bits(&a.energy), energy_bits(&b.energy), "{what}");
    assert_eq!(bits(&a.forces), bits(&b.forces), "{what}");
    assert_eq!(
        a.virial.map(|r| r.map(f64::to_bits)),
        b.virial.map(|r| r.map(f64::to_bits)),
        "{what}"
    );
}

/// The rule replayed from the trajectory, apart from the force field: a build whenever an atom's
/// displacement since the last build, `d − L round(d/L)` per axis, is longer than half the skin.
struct Replay {
    skin: f64,
    lengths: [f64; 3],
    reference: Option<Vec<[f64; 3]>>,
    builds: u64,
}

impl Replay {
    fn new(skin: f64, cell: PeriodicBox) -> Replay {
        Replay {
            skin,
            lengths: cell.lengths(),
            reference: None,
            builds: 0,
        }
    }

    fn see(&mut self, at: &[[f64; 3]]) {
        let l = self.lengths;
        let moved = match &self.reference {
            None => true,
            Some(r) => at.iter().zip(r).any(|(p, q)| {
                let mut u2 = 0.0;
                for a in 0..3 {
                    let d = p[a] - q[a];
                    let u = d - l[a] * (d / l[a]).round();
                    u2 += u * u;
                }
                u2.sqrt() > 0.5 * self.skin
            }),
        };
        if moved {
            self.builds += 1;
            self.reference = Some(at.to_vec());
        }
    }
}

/// `steps` of `md` with `plain` and with `listed` side by side from `at`, `dt` a step: positions,
/// velocities and potential energy the same bits after every step, `check` at every `every`th,
/// and the listed field's builds — read by `status` — exactly the replay's at `skin`, the skin
/// the test asked for (metres), not the one the list reports. Returns the builds and the largest
/// displacement of any atom.
#[allow(clippy::too_many_arguments)]
fn side_by_side(
    plain: &dyn Potential,
    listed: &dyn Potential,
    status: &dyn Fn() -> pantometry_forcefield::NeighbourListStatus,
    skin: f64,
    cell: PeriodicBox,
    md: MolecularDynamics,
    at: &[[f64; 3]],
    dt: f64,
    steps: usize,
    every: usize,
    check: &dyn Fn(&[[f64; 3]]),
) -> (u64, f64) {
    let (mut ma, mut mb) = (md.clone(), md);
    let (mut a, mut b) = (at.to_vec(), at.to_vec());
    let mut replay = Replay::new(skin, cell);
    println!("the list's skin {:e} m, asked for {skin:e}", status().skin);
    assert_eq!(
        status().skin.to_bits(),
        skin.to_bits(),
        "the skin asked for"
    );
    ma.prepare(plain, &a);
    mb.prepare(listed, &b);
    replay.see(&b);
    check(&b);
    for k in 1..=steps {
        ma.step(plain, &mut a, dt);
        mb.step(listed, &mut b, dt);
        replay.see(&b);
        assert_eq!(bits(&a), bits(&b), "positions at step {k}");
        assert_eq!(bits(ma.velocities()), bits(mb.velocities()), "step {k}");
        assert_eq!(
            ma.potential_energy().map(f64::to_bits),
            mb.potential_energy().map(f64::to_bits),
            "step {k}"
        );
        assert_eq!(status().builds, replay.builds, "builds at step {k}");
        if k % every == 0 {
            check(&b);
        }
    }
    let moved = b
        .iter()
        .zip(at)
        .map(|(p, q)| common::len(common::sub(*p, *q)))
        .fold(0.0f64, f64::max);
    (status().builds, moved)
}

fn bath(seed: u64) -> Bath {
    Bath::Langevin {
        temperature: 300.0,
        friction: 5e12,
        seed,
    }
}

// W1: UFF's flexible water, from `tests/a_molecule_in_a_periodic_box.rs`.

fn splitmix(seed: u64, i: u64) -> u64 {
    let mut z = seed.wrapping_add((i + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn uniform(seed: u64, i: u64) -> f64 {
    (splitmix(seed, i) >> 11) as f64 / (1u64 << 53) as f64
}

/// `k³` flexible UFF waters on a grid filling a cube of side `side` Å, each turned by its own
/// angles and moved off its site by up to 0.3 Å, with TIP3P's charges: W1's box.
fn uff_water(k: usize, side: f64, cutoff: f64) -> (PeriodicForceField, Vec<[f64; 3]>, Component) {
    let spacing = side / k as f64;
    let (r, half) = (0.9572, 0.5 * 104.52f64.to_radians());
    let seed = 0xD1CE;
    let mut atoms = Vec::new();
    let mut names = Vec::new();
    for i in 0..k * k * k {
        let g = [i / (k * k), (i / k) % k, i % k];
        let n = i as u64;
        let o = [0, 1, 2]
            .map(|a| (g[a] as f64 + 0.5) * spacing + 0.6 * (uniform(seed, 7 * n + a as u64) - 0.5));
        let (t, p, s) = (
            std::f64::consts::PI * uniform(seed, 7 * n + 3),
            2.0 * std::f64::consts::PI * uniform(seed, 7 * n + 4),
            2.0 * std::f64::consts::PI * uniform(seed, 7 * n + 5),
        );
        let u = [t.sin() * p.cos(), t.sin() * p.sin(), t.cos()];
        let w0 = if u[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let c = common::cross(u, w0);
        let lc = common::len(c);
        let a = c.map(|x| x / lc);
        let b = common::cross(u, a);
        let v = [0, 1, 2].map(|m| s.cos() * a[m] + s.sin() * b[m]);
        let h1 = [0, 1, 2].map(|m| o[m] + r * (half.cos() * u[m] + half.sin() * v[m]));
        let h2 = [0, 1, 2].map(|m| o[m] + r * (half.cos() * u[m] - half.sin() * v[m]));
        names.push([format!("O{i}"), format!("H{i}A"), format!("H{i}B")]);
        atoms.push([o, h1, h2]);
    }
    let mut list = Vec::new();
    let mut bonds = Vec::new();
    for (n, s) in names.iter().zip(&atoms) {
        list.push((n[0].as_str(), "O", s[0]));
        list.push((n[1].as_str(), "H", s[1]));
        list.push((n[2].as_str(), "H", s[2]));
        bonds.push((n[0].as_str(), n[1].as_str(), "SING", false));
        bonds.push((n[0].as_str(), n[2].as_str(), "SING", false));
    }
    let c = common::entry(&list, &bonds);
    let charges: Vec<f64> = c
        .atoms()
        .iter()
        .map(|a| {
            if a.name.starts_with('O') {
                -0.834
            } else {
                0.417
            }
        })
        .collect();
    let cell = PeriodicBox::cubic(side * ANGSTROM);
    let p = EwaldParameters::for_accuracy(&cell, cutoff * ANGSTROM, 1e-5);
    let field = PeriodicForceField::new(&c, &uff::assign(&c), cell, p)
        .expect("water is supported")
        .with_charges(charges);
    let at = common::positions(&c);
    (field, at, c)
}

/// Steepest descent, no atom moved more than 0.02 Å an iteration: W1's box is built with
/// contacts its UFF hydrogens repel at thousands of kcal/mol/Å, which a 300 K start would turn
/// into a temperature of tens of thousands of kelvin.
fn relax(field: &PeriodicForceField, at: &mut [[f64; 3]], iterations: usize) {
    for _ in 0..iterations {
        let f = field.evaluate(at).forces;
        let largest = f.iter().map(|x| common::len(*x)).fold(0.0f64, f64::max);
        let scale = 0.02 * ANGSTROM / largest;
        for (p, x) in at.iter_mut().zip(&f) {
            for a in 0..3 {
                p[a] += scale * x[a];
            }
        }
    }
}

/// A field and its listed copy evaluate to the same bits at `at`. The listed one is the
/// trajectory's own, at positions it has just been evaluated at, so the check builds nothing.
fn fields_agree<'a>(
    plain: &'a PeriodicForceField,
    listed: &'a PeriodicForceField,
) -> impl Fn(&[[f64; 3]]) + 'a {
    move |at| same_evaluation(&plain.evaluate(at), &listed.evaluate(at), "evaluate")
}

/// **W1: flexible UFF water, no constraint**: 64 molecules in 12.4 Å at `r_c` = 4.6 Å, five
/// cells of half the cutoff an axis, so that the 125 about each cell are the whole box; relaxed,
/// then `steps` of 0.5 fs in a 300 K bath. Returns the builds.
fn flexible_water(relaxation: usize, steps: usize, skin: f64, every: usize) -> u64 {
    let (plain, mut at, c) = uff_water(4, 12.4, 4.6);
    relax(&plain, &mut at, relaxation);
    let listed = plain.clone().with_neighbour_list(skin * ANGSTROM);
    let md = MolecularDynamics::for_elements(c.atoms().iter().map(|a| a.element))
        .with_bath(bath(11))
        .thermalised(&at, 300.0, 11);
    let check = fields_agree(&plain, &listed);
    let (builds, moved) = side_by_side(
        &plain,
        &listed,
        &|| listed.neighbour_list().unwrap(),
        skin * ANGSTROM,
        plain.cell(),
        md,
        &at,
        0.5 * FS,
        steps,
        every,
        &check,
    );
    println!(
        "W1, {skin} Å: {builds} builds in {} evaluations, an atom moved {:.2} Å",
        steps + 1,
        moved / ANGSTROM
    );
    builds
}

/// **W2: rigid TIP3P under SETTLE**: `per_side³` waters at `cutoff` Å, `steps` of 1 fs.
fn rigid_water(per_side: usize, cutoff: f64, steps: usize, skin: f64, every: usize) -> u64 {
    let water = WaterBox::lattice(per_side, 0x5E77);
    let plain = water.force_field(cutoff * ANGSTROM, 1e-5);
    let listed = plain.clone().with_neighbour_list(skin * ANGSTROM);
    let at = water.positions().to_vec();
    let md = water
        .dynamics()
        .with_bath(bath(12))
        .thermalised(&at, 300.0, 12);
    let check = fields_agree(&plain, &listed);
    let (builds, moved) = side_by_side(
        &plain,
        &listed,
        &|| listed.neighbour_list().unwrap(),
        skin * ANGSTROM,
        plain.cell(),
        md,
        &at,
        FS,
        steps,
        every,
        &check,
    );
    println!(
        "W2, {} waters, {skin} Å: {builds} builds in {} evaluations, an atom moved {:.2} Å",
        water.count(),
        steps + 1,
        moved / ANGSTROM
    );
    builds
}

/// The decoupling's energy by term at every state, its forces at the first — the trajectory's —
/// and its couplings at every state, the same bits with and without the list.
fn decouplings_agree<'a>(
    plain: &'a PeriodicDecoupling,
    listed: &'a PeriodicDecoupling,
    states: &'a [Lambda],
) -> impl Fn(&[[f64; 3]]) + 'a {
    move |at| {
        for &l in states {
            assert_eq!(
                energy_bits(&plain.energy(at, l)),
                energy_bits(&listed.energy(at, l)),
                "{l:?}"
            );
        }
        let (mut f1, mut f2) = (vec![[0.0; 3]; at.len()], vec![[0.0; 3]; at.len()]);
        let e1 = plain.energy_and_forces(at, states[0], &mut f1);
        let e2 = listed.energy_and_forces(at, states[0], &mut f2);
        assert_eq!(e1.to_bits(), e2.to_bits());
        assert_eq!(bits(&f1), bits(&f2));
        let (c1, c2) = (plain.couplings(at, states), listed.couplings(at, states));
        for (a, b) in c1.iter().zip(&c2) {
            assert_eq!(a.energy.to_bits(), b.energy.to_bits());
            assert_eq!(a.gradient.map(f64::to_bits), b.gradient.map(f64::to_bits));
        }
    }
}

/// The trajectory's state first, whose forces are checked; the energy alone at the others.
const STATES: [Lambda; 5] = [
    Lambda::new(0.0, 0.5, 0.7),
    Lambda::COUPLED,
    Lambda {
        restraint: 1.0,
        electrostatics: 0.5,
        van_der_waals: 1.0,
    },
    Lambda {
        restraint: 1.0,
        electrostatics: 0.0,
        van_der_waals: 0.4,
    },
    Lambda {
        restraint: 1.0,
        electrostatics: 0.0,
        van_der_waals: 0.0,
    },
];

/// **W3: benzene in TIP3P under SETTLE and SHAKE, decoupled half-way**, in `per_side³` lattice
/// waters less those it overlaps, at `cutoff` Å: `steps` of 1 fs at `λ = (0, 0.5, 0.7)`, and
/// every check through the decoupling at five states and its couplings.
fn benzene_in_water(per_side: usize, cutoff: f64, steps: usize, skin: f64, every: usize) -> u64 {
    let b = Benzene::new(None);
    let s = Solvated::new(&b, per_side, 0x3B, 2.3 * ANGSTROM, cutoff * ANGSTROM, 1e-5);
    let plain = s.decoupling();
    let field = s.field.clone().with_neighbour_list(skin * ANGSTROM);
    let listed = PeriodicDecoupling::new(&field, &s.group).unwrap();
    let lambda = Lambda::new(0.0, 0.5, 0.7);
    let md = s
        .dynamics(&b)
        .with_bath(bath(13))
        .thermalised(&s.at, 300.0, 13);
    let check = decouplings_agree(&plain, &listed, &STATES);
    let (builds, _) = side_by_side(
        &AtLambda {
            hamiltonian: &plain,
            lambda,
        },
        &AtLambda {
            hamiltonian: &listed,
            lambda,
        },
        &|| listed.field().neighbour_list().unwrap(),
        skin * ANGSTROM,
        s.cell(),
        md,
        &s.at,
        FS,
        steps,
        every,
        &check,
    );
    println!(
        "W3, {} atoms at {cutoff} Å, {skin} Å: {builds} builds in {} steps",
        s.at.len(),
        steps + 1
    );
    builds
}

/// **W4: a pocket of 181L at `radius` Å and benzene, the protein frozen**, through the mesh, at
/// `cutoff` Å: W4's dynamics — SETTLE, SHAKE on every bond to a hydrogen of the complex, the
/// protein's atoms held — `steps` of 1 fs at `λ = (0, 0.5, 0.7)`. **A frozen atom never moves**:
/// its bits are the start's at the last check, which is the last step, and the builds are the
/// replay's, whose displacements are the frozen atoms' zeros and the rest's.
fn pocket_frozen(radius: f64, cutoff: f64, steps: usize, skin: f64, every: usize) -> u64 {
    assert_eq!(steps % every, 0, "the last check is at the last step");
    let system = protein::system();
    let pocket = Binding::new(&system, radius * ANGSTROM).expect("a pocket of 181L");
    let solvation = Solvation {
        margin: 3.5 * ANGSTROM,
        cutoff: cutoff * ANGSTROM,
        seed: 0xE5D,
        ..Solvation::w4()
    };
    let s = SolvatedComplex::new(&pocket, &solvation).expect("the pocket in water");
    assert!(s.field().ewald().mesh().is_some());
    let plain = s.decoupling();
    let field = s.field().clone().with_neighbour_list(skin * ANGSTROM);
    let listed = PeriodicDecoupling::new(&field, &s.ligand_mask()).unwrap();
    let lambda = Lambda::new(0.0, 0.5, 0.7);
    let md = s
        .dynamics(Flexibility::Frozen)
        .with_bath(bath(14))
        .thermalised(s.positions(), 300.0, 14);
    let frozen: Vec<bool> = md.frozen().to_vec();
    let held = frozen.iter().filter(|&&f| f).count();
    assert!(held > 50, "{held} frozen");
    let check = decouplings_agree(&plain, &listed, &STATES);
    let at = s.positions().to_vec();
    let end = std::cell::RefCell::new(Vec::new());
    let (builds, moved) = side_by_side(
        &AtLambda {
            hamiltonian: &plain,
            lambda,
        },
        &AtLambda {
            hamiltonian: &listed,
            lambda,
        },
        &|| listed.field().neighbour_list().unwrap(),
        skin * ANGSTROM,
        s.cell(),
        md,
        &at,
        FS,
        steps,
        every,
        &|x: &[[f64; 3]]| {
            check(x);
            *end.borrow_mut() = x.to_vec();
        },
    );
    let end = end.into_inner();
    for (i, f) in frozen.iter().enumerate() {
        if *f {
            assert_eq!(
                end[i].map(f64::to_bits),
                at[i].map(f64::to_bits),
                "atom {i}"
            );
        }
    }
    let mobile_moved = (0..at.len()).any(|i| !frozen[i] && end[i] != at[i]);
    assert!(mobile_moved, "the mobile atoms moved");
    println!(
        "W4, the {radius} Å pocket, {} atoms, {held} frozen, {skin} Å: {builds} builds in {} steps, an atom moved {:.2} Å",
        at.len(),
        steps + 1,
        moved / ANGSTROM
    );
    builds
}

/// **W1, flexible water**: 60 steps of 0.5 fs with a 0.3 Å skin, so that the list is built
/// several times in them.
#[test]
fn flexible_water_without_constraints_is_the_same_bits() {
    let builds = flexible_water(100, 60, 0.3, 10);
    assert!(builds >= 3, "{builds} builds");
}

/// **W2, rigid water under SETTLE**: 125 waters in 15.5 Å at 4.6 Å, six cells an axis; 60 steps
/// of 1 fs with a 0.6 Å skin.
#[test]
fn rigid_water_under_settle_is_the_same_bits() {
    let builds = rigid_water(5, 4.6, 60, 0.6, 10);
    assert!(builds >= 3, "{builds} builds");
}

/// **W3, benzene under SETTLE and SHAKE through the decoupling**: the 27-water box of
/// `tests/benzene_hydrated_in_tip3p.rs` at 4.6 Å, too small for cells, and 64 waters at 4.6 Å,
/// five cells an axis; 60 steps of 1 fs with a 0.6 Å skin, the decoupling checked every 20.
#[test]
fn benzene_in_water_under_settle_and_shake_is_the_same_bits() {
    for per_side in [3, 4] {
        let builds = benzene_in_water(per_side, 4.6, 60, 0.6, 20);
        assert!(builds >= 3, "{per_side}: {builds} builds");
    }
}

/// **W4, 181L's 3 Å pocket and benzene with the protein frozen** (the pocket and benzene 82
/// atoms, neutral; 832 with the water; `r_c` = 6 Å): 40 steps of 1 fs with a 0.8 Å skin, the
/// decoupling checked every 20.
#[test]
fn the_pocket_with_its_protein_frozen_is_the_same_bits() {
    let builds = pocket_frozen(3.0, 6.0, 40, 0.8, 20);
    assert!(builds >= 3, "{builds} builds");
}

/// **The same, longer**: what the default tests were first written as — hundreds of steps with
/// skins of 0.3–1 Å, a build every six to ten steps — with 216 waters at 6 Å and 181L's 5 Å pocket
/// and its chloride at 5.5 Å. Nothing the short ones do not check; more of it.
#[test]
#[ignore = "about 4.5 minutes unoptimised: run with -- --ignored --nocapture"]
fn the_long_trajectories_are_the_same_bits() {
    assert!(flexible_water(300, 300, 0.3, 10) >= 4);
    assert!(rigid_water(6, 6.0, 300, 1.0, 10) >= 4);
    for (per_side, cutoff) in [(3, 4.6), (5, 6.0)] {
        assert!(benzene_in_water(per_side, cutoff, 300, 1.0, 50) >= 4);
    }
    assert!(pocket_frozen(5.0, 5.5, 200, 1.0, 50) >= 3);
}

/// **The cost on W4's box**: the whole of 181L and benzene in 9 160 TIP3P waters with nine
/// chlorides, 30 105 atoms, through the mesh at δ = 10⁻⁶ and `r_c` = 9 Å — W4's potential — in
/// release on one core. Melted with the protein frozen, 0.15 ps at 0.5 fs in a 50 ps⁻¹ bath, with
/// a list (which changes no bit); then from that one state, for the cell list and each skin:
/// milliseconds an evaluation with no build, an evaluation that builds, and 50 steps of 2 fs in a
/// 1 ps⁻¹ bath with the builds they took. Each list's end — positions and velocities — is
/// asserted to be the cell list's to the bit before its row is printed, so a printed row is one
/// whose trajectory was. The timings are reported, not asserted: a timing is the machine's.
#[test]
#[ignore = "W4's 30 105-atom box, about seven minutes with --release: run with --release -- --ignored --nocapture"]
fn the_cost_on_w4s_box_measured() {
    let clock = std::time::Instant::now();
    let b = Binding::new(&protein::system(), 100.0 * ANGSTROM).expect("the whole of 181L");
    let s = SolvatedComplex::new(&b, &Solvation::w4()).expect("181L in water");
    let n = s.positions().len();
    println!(
        "{n} atoms in {:?} Å, built in {:.0} s",
        s.cell()
            .lengths()
            .map(|x| (x / ANGSTROM * 100.0).round() / 100.0),
        clock.elapsed().as_secs_f64()
    );
    let plain = s.field().clone();
    let template = s.dynamics(Flexibility::Frozen);
    let mut at = s.positions().to_vec();
    let melt_field = plain.clone().with_neighbour_list(1.5 * ANGSTROM);
    let mut md = template
        .clone()
        .with_bath(Bath::Langevin {
            temperature: 298.15,
            friction: 50e12,
            seed: 0xE01,
        })
        .thermalised(&at, 298.15, 0xE01);
    let t = std::time::Instant::now();
    md.run(&melt_field, &mut at, 0.5 * FS, 300);
    println!(
        "melted 0.15 ps at 0.5 fs: {:.0} s, T {:.1} K, {} builds",
        t.elapsed().as_secs_f64(),
        md.temperature(),
        melt_field.neighbour_list().unwrap().builds
    );
    let velocities = md.velocities().to_vec();
    let md = template
        .with_velocities(velocities)
        .with_bath(Bath::Langevin {
            temperature: 298.15,
            friction: 1e12,
            seed: 0xE02,
        });
    let ms = |f: &dyn Fn(), k: usize| {
        let t = std::time::Instant::now();
        for _ in 0..k {
            f();
        }
        t.elapsed().as_secs_f64() / k as f64 * 1e3
    };
    let steps = 50;
    // The cell list's end, positions and velocities, which every list's must be before its row
    // is printed.
    let mut reference: Option<EndBits> = None;
    println!("| | pairs kept | ms an evaluation | ms one that builds | ms a step at 2 fs | builds in {steps} steps |");
    println!("| --- | --- | --- | --- | --- | --- |");
    for skin in [None, Some(1.0), Some(1.5), Some(2.0), Some(2.5)] {
        let field = match skin {
            None => plain.clone(),
            Some(x) => plain.clone().with_neighbour_list(x * ANGSTROM),
        };
        field.evaluate(&at);
        let evaluation = ms(&|| drop(field.evaluate(&at)), 3);
        let building = match skin {
            None => evaluation,
            Some(x) => ms(
                &|| {
                    drop(
                        plain
                            .clone()
                            .with_neighbour_list(x * ANGSTROM)
                            .evaluate(&at),
                    )
                },
                3,
            ),
        };
        let mut m = md.clone();
        let mut x = at.clone();
        m.prepare(&field, &x);
        let before = field.neighbour_list().map_or(0, |l| l.builds);
        let t = std::time::Instant::now();
        m.run(&field, &mut x, 2.0 * FS, steps);
        let step = t.elapsed().as_secs_f64() / steps as f64 * 1e3;
        let end = (bits(&x), bits(m.velocities()));
        match &reference {
            None => reference = Some(end),
            Some(r) => assert!(
                *r == end,
                "{skin:?}: the cell list's trajectory, to the bit"
            ),
        }
        let status = field.neighbour_list();
        println!(
            "| {} | {} | {evaluation:.1} | {building:.1} | {step:.1} | {} |",
            skin.map_or("cell list".to_string(), |x| format!("skin {x} Å")),
            status.map_or("—".to_string(), |l| l.pairs.to_string()),
            status.map_or("—".to_string(), |l| (l.builds - before).to_string()),
        );
    }
    println!("in all {:.0} s", clock.elapsed().as_secs_f64());
}
