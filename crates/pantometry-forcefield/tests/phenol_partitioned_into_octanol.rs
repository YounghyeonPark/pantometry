//! **Phenol decoupled from 1-octanol and from water, for its log P**: the third solute through
//! 4c-2's legs, and the first polar one. Benzene's log P came out at 1.92 ± 0.23 against
//! experiment's 2.13 (`benzene_partitioned_into_octanol.rs`) and p-xylene's at 3.06 ± 0.22 against
//! 3.15 (`p_xylene_partitioned_into_octanol.rs`). Both are hydrocarbons, so neither has asked the
//! charges for much: what each leg measured is mostly the van der Waals segment. Phenol is benzene
//! with a hydroxyl, a hydrogen to give and an oxygen to take one, in water and in octanol both, and
//! its hydration free energy is −6.60 kcal/mol where benzene's is −0.90. This asks whether the
//! method carries to a solute whose answer is in its charges.
//!
//! It brings no GAFF type the two legs do not already have: `ca` and `ha` are benzene's, and `oh`
//! and `ho` are octanol's own hydroxyl.
//!
//! **Kept out of the published crate** (`Cargo.toml`'s `exclude`), as benzene's and p-xylene's
//! are: GAFF's values are read from `gaff.dat`, which is GPL v3.
//!
//! # The model and the boxes
//!
//! Everything is benzene's and p-xylene's, built by the same lines (`tests/partition/`'s `_of`
//! builders):
//!
//! - **The solute**: `IPH` of the Chemical Component Dictionary, PHENOL, `C6 H6 O`. It is typed for
//!   UFF, relaxed in UFF's vacuum with its AM1-BCC charges, and every X–H is held, the O–H among
//!   them, as octanol's is.
//! - **Its non-bonded terms**: GAFF's `ca`, `ha`, `oh` and `ho` and FreeSolv's AM1-BCC charges
//!   (`mobley_20524.mol2`, `tests/freesolv/`), mapped by graph isomorphism.
//! - **The rules**: Lorentz–Berthelot and GAFF's 1-4 scales.
//! - **The octanol box**: 4c-1's 192-site lattice, compensated for the solute's volume ([`VOLUME`]).
//! - **The water box**: W3's.
//! - **The melts, the 29 candidates, the 13 windows, the insertion, the protocol and the seeds**:
//!   benzene's. The results go under their own `PANTOMETRY_OCTANOL_DIR`. Every window and start
//!   carries a digest of the model, so nothing of benzene's or p-xylene's can be read back as
//!   phenol's.
//!
//! # What it is compared with, registered before the run
//!
//! Experiment's `log K_ow` is **1.46**: C. Hansch, A. Leo and D. Hoekman, *Exploring QSAR:
//! Hydrophobic, Electronic and Steric Constants*, ACS 1995, p. 20. It was read on 2026-11-08 as
//! PubChem CID 996's HSDB record quotes it, the source benzene's 2.13 and p-xylene's 3.15 are from.
//! The record's ICSC and HMDB entries are 1.46 too, and its Sangster entry 1.6. FreeSolv's hydration
//! free energy is −6.60 ± 0.20 by experiment and −5.71 ± 0.03 for GAFF/AM1-BCC on GAFF's bonded
//! terms: the model is known to hydrate phenol 0.89 kcal/mol too weakly, which by itself would put
//! log P 0.65 high. Whether octanol's leg carries the same error is what the run says.
//!
//! **The gate, written before any window ran.** Let `σ` be log P's error with **both legs'** charge
//! segments' seed spread in it: benzene's and p-xylene's reseeded the octanol's alone, but phenol's
//! hydration free energy is mostly its water leg's charge segment, so that one is reseeded too.
//! Then:
//! - the verdict is **uninformative** if `2σ ≥ 0.46`, whatever log P is. A method blind to the
//!   hydroxyl gives this method's benzene, 1.92, not experiment's 2.13, so a band reaching 1.92 from
//!   1.46 would pass it. Experiment's benzene less phenol is 0.67, and was the ceiling this gate
//!   first had: between 0.46 and 0.67 it would have called a result of 1.92 a method that carries.
//!   σ takes the larger of the seeds' half-difference and BAR's, so it widens as the seeds
//!   disagree, and this is its ceiling;
//! - otherwise the method **carries to phenol** if `|log P − 1.46| ≤ 2σ`;
//! - it **does not** if `|log P − 1.46| > 3σ`;
//! - between the two, it is not settled.
//!
//! σ leaves out three things, and they are stated beside the verdict rather than added to it. The
//! van der Waals segments' own seed spread is not measured. A difference of two seeds is one degree
//! of freedom. And experiment's own sources differ by 0.14 (1.46 against Sangster's 1.6).
//!
//! The difference from benzene, `Δ log P = log P(phenol) − log P(benzene)`, is reported beside it
//! against experiment's −0.67, the same reading on the same two legs, with whether it is below
//! benzene's by more than twice its error: whether the hydroxyl is seen at all. It is reported, not
//! gated. Nothing in this file asserts either number.
//!
//! **Much of this is the water leg's known error measured again.** FreeSolv's 0.89 kcal/mol alone
//! would put log P at about 2.11 if octanol's leg were exact; a result near there says the octanol
//! leg is right and the water's is FreeSolv's, not that the method has failed in a new way.
//!
//! # What it measured
//!
//! Not run yet.

// The sums over atoms and axes are written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

mod freesolv;
mod gaff;
mod hydration;
mod nonbonded;
mod octanol;
mod partition;

use gaff::{gaff, Model, ONE_FOUR};
use hydration::{Benzene, Solvated};
use nonbonded::Vacuum;
use octanol::{joined, Liquid, ANGSTROM, DENSITY, FS, KELVIN};
use pantometry_forcefield::free_energy::refine_schedule;
use pantometry_forcefield::uff::AVOGADRO;
use pantometry_forcefield::{
    Alchemical, AtLambda, Bath, Component, Element, Estimate, Lambda, Mixing, MolecularDynamics,
    PeriodicDecoupling, PeriodicForceField, Potential, Protocol,
};
use pantometry_units::BOLTZMANN;
use partition::*;

/// Phenol, `IPH` of the Chemical Component Dictionary, byte for byte.
const IPH: &str = include_str!("../components/IPH.cif");

/// Its atoms.
const N: usize = 13;

/// **Phenol's volume in the solvent, m³ a molecule**: 94.11 g/mol over 1.071 g/cm³ (O'Neil, *The
/// Merck Index*, 13th ed., Merck 2001, p. 1299, as PubChem CID 996's HSDB record quotes it, read
/// 2026-11-08), 87.9 cm³/mol. It stands in for the partial molar volume in octanol, as benzene's 89
/// and p-xylene's 123.3 do. Phenol melts at 41 °C and the record gives no temperature, so this is
/// most likely the crystal's; the liquid's 1.0545 at 45 °C (the record's PAC entry) would be 89.2,
/// 1.6% more. [`the_octanol_is_at_its_density_around_phenol`] holds the octanol to its density
/// around whatever this is; nothing here checks the number itself.
const VOLUME: f64 = 94.11 / 1.071 * 1e-6 / AVOGADRO;

/// Experiment's `log K_ow` for phenol: see the module documentation.
const LOG_P_PHENOL: f64 = 1.46;

/// **The gate's ceiling on 2σ**, 0.46: the smaller of experiment's benzene less phenol, 0.67, and
/// this method's benzene less experiment's phenol, 1.92 − 1.46. A method blind to the hydroxyl
/// would give this method's benzene, not experiment's, so a band reaching 1.92 from 1.46 passes it;
/// for p-xylene the first was the smaller, and for phenol it is the second.
fn ceiling() -> f64 {
    (LOG_P - LOG_P_PHENOL).min(BENZENE_LOG_P.0 - LOG_P_PHENOL)
}

/// Benzene's log P from its two direct legs over two seeds of the octanol charge segment, and its
/// error (`the_octanol_charges_reseeded_measured`, run 2026-10-24): what [`Δ log
/// P`](partition_measured) is taken from.
const BENZENE_LOG_P: (f64, f64) = (1.92, 0.23);

/// Phenol with FreeSolv's AM1-BCC charges, relaxed in UFF's vacuum with them and put on its
/// constraints; and its model.
fn phenol() -> (Benzene, Model) {
    let c = Component::from_ccd(IPH).expect("IPH");
    let model = Model::phenol(&c);
    (Benzene::from_ccd(IPH, Some(model.charges.clone())), model)
}

/// Each atom's neighbours, from the bonds.
fn neighbours(c: &Component) -> Vec<Vec<usize>> {
    let mut next = vec![Vec::new(); c.atoms().len()];
    for b in c.bonds() {
        next[b.atoms[0]].push(b.atoms[1]);
        next[b.atoms[1]].push(b.atoms[0]);
    }
    next
}

/// Every atom's bond distance from atom `from`, by breadth-first search on the bonds.
fn distances_from(c: &Component, from: usize) -> Vec<usize> {
    let next = neighbours(c);
    let mut d = vec![usize::MAX; next.len()];
    d[from] = 0;
    let mut queue = std::collections::VecDeque::from([from]);
    while let Some(a) = queue.pop_front() {
        for &b in &next[a] {
            if d[b] == usize::MAX {
                d[b] = d[a] + 1;
                queue.push_back(b);
            }
        }
    }
    d
}

/// Each atom's GAFF type by its place in the molecule, not by FreeSolv's file: the oxygen is `oh`
/// and the hydrogen on it `ho`; a carbon is `ca` and a hydrogen on one `ha`.
fn types_by_structure(c: &Component) -> Vec<&'static str> {
    let next = neighbours(c);
    (0..next.len())
        .map(|i| match c.atoms()[i].element {
            Element::O => "oh",
            Element::C => "ca",
            Element::H if c.atoms()[next[i][0]].element == Element::O => "ho",
            Element::H => "ha",
            other => panic!("phenol has no {other:?}"),
        })
        .collect()
}

/// **Each atom's AM1-BCC charge by its place**, as FreeSolv's file gives it: a ring carbon by its
/// bond distance from the carbon that carries the oxygen (`ipso` 0, `ortho` 1, `meta` 2, `para` 3),
/// a ring hydrogen by its carbon's, and the hydroxyl's two atoms by name. The two ortho and the two
/// meta positions are each one charge, as the file writes them.
fn charge_by_place(c: &Component) -> Vec<f64> {
    let next = neighbours(c);
    let o = (0..next.len())
        .find(|&i| c.atoms()[i].element == Element::O)
        .expect("phenol's oxygen");
    let ipso = *next[o]
        .iter()
        .find(|&&j| c.atoms()[j].element == Element::C)
        .expect("the carbon the oxygen is on");
    let d = distances_from(c, ipso);
    (0..next.len())
        .map(|i| match c.atoms()[i].element {
            Element::O => -0.4989,
            Element::H if next[i][0] == o => 0.4186,
            Element::C => [0.1230, -0.1852, -0.0940, -0.1659][d[i]],
            Element::H => [f64::NAN, 0.1409, 0.1337, 0.1324][d[next[i][0]]],
            other => panic!("phenol has no {other:?}"),
        })
        .collect()
}

/// Every pair of `c`'s atoms by its bond distance: `(pairs three bonds apart, pairs more than
/// three apart)`.
fn pairs_by_distance(c: &Component) -> (usize, usize) {
    let n = c.atoms().len();
    let (mut three, mut beyond) = (0, 0);
    for i in 0..n {
        let d = distances_from(c, i);
        for j in i + 1..n {
            match d[j] {
                3 => three += 1,
                x if x > 3 => beyond += 1,
                _ => {}
            }
        }
    }
    (three, beyond)
}

// ---------------------------------------------------------------------------------------------
// Default checks

/// The small boxes, as benzene's file builds them, around phenol: the field, the decoupling, the
/// solute alone in the box and in vacuum, the solvent alone, and a dynamics template.
struct Small {
    field: PeriodicForceField,
    decoupling: PeriodicDecoupling,
    solute_in_the_box: PeriodicForceField,
    solute_in_vacuum: Vacuum,
    solvent_alone: PeriodicForceField,
    at: Vec<[f64; 3]>,
    solvent: usize,
    template: MolecularDynamics,
}

/// The rules written out, for the systems the end states are held to.
fn rules(field: PeriodicForceField, vdw: Vec<[f64; 2]>, q: Vec<f64>) -> PeriodicForceField {
    let n = vdw.len();
    field
        .with_nonbonded(0..n, vdw)
        .expect("valid")
        .with_charges(q)
        .with_mixing(Mixing::LorentzBerthelot)
        .with_one_four(0.5, 1.0 / 1.2)
}

fn vacuum(b: &Benzene, m: &Model) -> Vacuum {
    Vacuum::new(
        &b.component,
        &b.types,
        m.vdw.clone(),
        m.charges.clone(),
        Mixing::LorentzBerthelot,
        [0.5, 1.0 / 1.2],
    )
}

fn small_octanol(seed: u64) -> (Small, Liquid) {
    let (b, m) = phenol();
    let o = partition::octanol();
    let (l, field, _) = octanol_box_of(&b, &m, VOLUME, &o, SMALL, SMALL_CUTOFF, 1e-5, seed, false);
    let cell = l.cell();
    let p = l.field.ewald().parameters();
    let om = Model::octanol(&o.component);
    let parts: Vec<&Component> = (0..l.octanols).map(|_| &o.component).collect();
    let (mut vdw, mut q) = (Vec::new(), Vec::new());
    for _ in 0..l.octanols {
        vdw.extend_from_slice(&om.vdw);
        q.extend_from_slice(&om.charges);
    }
    let group: Vec<bool> = (0..l.at.len()).map(|i| i < N).collect();
    let small = Small {
        decoupling: PeriodicDecoupling::new(&field, &group).expect("phenol decouples"),
        solute_in_the_box: rules(
            PeriodicForceField::new(&b.component, &b.types, cell, p).expect("phenol"),
            m.vdw.clone(),
            m.charges.clone(),
        ),
        solute_in_vacuum: vacuum(&b, &m),
        solvent_alone: rules(
            PeriodicForceField::new(&joined(&parts), &l.types[N..], cell, p).expect("octanol"),
            vdw,
            q,
        ),
        field,
        at: l.at.clone(),
        solvent: l.octanols,
        template: l.dynamics(),
    };
    (small, l)
}

/// Phenol in 64 lattice waters less those within 2.3 Å of it, 12.4 Å across, `r_c` = 6.0 Å, its
/// model swapped in: p-xylene's small water box, which phenol, the shorter of the two, fits.
fn small_water(seed: u64) -> Small {
    let (b, m) = phenol();
    let s = with_model(
        Solvated::new(&b, 4, seed, 2.3 * ANGSTROM, 6.0 * ANGSTROM, 1e-5),
        &m,
    );
    let cell = s.cell();
    let p = s.field.ewald().parameters();
    Small {
        decoupling: s.decoupling(),
        field: s.field.clone(),
        solute_in_the_box: rules(
            PeriodicForceField::new(&b.component, &b.types, cell, p).expect("phenol"),
            m.vdw.clone(),
            m.charges.clone(),
        ),
        solute_in_vacuum: vacuum(&b, &m),
        solvent_alone: PeriodicForceField::tip3p(cell, p, s.waters)
            .with_mixing(Mixing::LorentzBerthelot)
            .with_one_four(0.5, 1.0 / 1.2),
        at: s.at.clone(),
        solvent: s.waters,
        template: s.dynamics(&b),
    }
}

/// `s` after `steps` of constrained BAOAB at 2 fs at `lambda`.
fn moved(s: &Small, lambda: Lambda, steps: usize) -> Vec<[f64; 3]> {
    let mut at = s.at.clone();
    let mut md = s
        .template
        .clone()
        .with_bath(Bath::Langevin {
            temperature: KELVIN,
            friction: 5e12,
            seed: 17,
        })
        .thermalised(&at, KELVIN, 17);
    md.run(
        &AtLambda {
            hamiltonian: &s.decoupling,
            lambda,
        },
        &mut at,
        2.0 * FS,
        steps,
    );
    at
}

/// **The model is what it says**, in both small boxes.
/// - Every atom's GAFF `[x, D]` is that of its type **by its place in the molecule**
///   ([`types_by_structure`]), not by FreeSolv's type column, and FreeSolv's types agree. `ho`'s is
///   `[0, 0]`: the hydroxyl hydrogen is a bare charge, as octanol's is.
/// - The charges are FreeSolv's to the bit, **each by the atom's place** ([`charge_by_place`]):
///   the ipso carbon +0.1230, the ortho −0.1852, the meta −0.0940 and the para −0.1659; their
///   hydrogens +0.1409, +0.1337 and +0.1324; the oxygen −0.4989 and its hydrogen +0.4186. They sum
///   to zero, and none is QEq's. **The geometry was relaxed with them too**: relaxing with QEq's
///   instead passed every other check here, since the field's charges are swapped in afterwards.
/// - The mixing is Lorentz–Berthelot, with GAFF's 1-4 scales.
/// - The solute's own 1-4 pairs are the pairs three bonds apart, counted on the graph, and its
///   own pairs those and the ones further.
/// - **The same model is in both boxes.**
#[test]
fn the_model_is_what_it_says() {
    let (oct, _) = small_octanol(0xE5D);
    let wat = small_water(0xE5E);
    let (b, m) = phenol();
    assert_eq!(b.component.atoms().len(), N);
    let kinds = types_by_structure(&b.component);
    assert_eq!(
        m.types,
        kinds.iter().map(|k| k.to_string()).collect::<Vec<_>>(),
        "FreeSolv's types are the structure's"
    );
    let want = charge_by_place(&b.component);
    // The geometry was relaxed with these charges, not QEq's: nothing below sees which, since the
    // field's charges are swapped in after, and a relaxation with QEq's passed every check here.
    // Read off the force field the minimiser ran on, not off `b.charges`, which is a second copy
    // of the same input and would not see QEq's put on that field alone.
    assert_eq!(
        b.vacuum
            .charges()
            .iter()
            .map(|x| x.to_bits())
            .collect::<Vec<_>>(),
        want.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        "phenol was relaxed with other charges than its model's"
    );
    let qeq = Benzene::from_ccd(IPH, None).charges;
    let (three, beyond) = pairs_by_distance(&b.component);
    println!("phenol: {three} pairs three bonds apart, {beyond} further");
    for (name, s) in [("octanol", &oct), ("water", &wat)] {
        let f = &s.field;
        let (vdw, q) = (f.van_der_waals(), f.charges());
        for i in 0..N {
            assert_eq!(
                vdw[i].map(f64::to_bits),
                gaff(kinds[i]).map(f64::to_bits),
                "{name} {i}"
            );
            assert_eq!(
                q[i].to_bits(),
                want[i].to_bits(),
                "{name}: atom {i} ({}) carries {} e and its place {} e",
                kinds[i],
                q[i],
                want[i]
            );
            assert_ne!(q[i].to_bits(), qeq[i].to_bits(), "{name}: QEq's charge");
        }
        let net: f64 = q[..N].iter().sum();
        assert!(net.abs() < 1e-12, "{name}: phenol's {net} e");
        assert_eq!(f.mixing(), Mixing::LorentzBerthelot, "{name}");
        assert_eq!(f.one_four(), ONE_FOUR, "{name}");
        assert_eq!(s.decoupling.own_pairs().len(), three + beyond, "{name}");
        assert_eq!(
            s.decoupling
                .own_scales()
                .iter()
                .filter(|x| **x == ONE_FOUR)
                .count(),
            three,
            "{name}"
        );
        println!(
            "{name}: phenol and {} solvent molecules, {} atoms",
            s.solvent,
            q.len()
        );
    }
    let bits =
        |v: &[[f64; 2]]| -> Vec<[u64; 2]> { v.iter().map(|x| x.map(f64::to_bits)).collect() };
    assert_eq!(
        bits(&oct.field.van_der_waals()[..N]),
        bits(&wat.field.van_der_waals()[..N])
    );
    let charge_bits = |f: &PeriodicForceField| -> Vec<u64> {
        f.charges()[..N].iter().map(|x| x.to_bits()).collect()
    };
    assert_eq!(charge_bits(&oct.field), charge_bits(&wat.field));
    assert_eq!(oct.decoupling.own_scales(), wat.decoupling.own_scales());
}

/// **The octanol is at its density in the volume phenol leaves it**, in the production's box, the
/// fixed point reached; and the compensation does something.
#[test]
fn the_octanol_is_at_its_density_around_phenol() {
    let (b, m) = phenol();
    let o = partition::octanol();
    let build = |compensate| {
        octanol_box_of(
            &b,
            &m,
            VOLUME,
            &o,
            octanol::LIQUID,
            9.0 * ANGSTROM,
            1e-5,
            0x0C71,
            compensate,
        )
    };
    let (l, _, lattice) = build(true);
    let rho = l.octanols as f64 * o.mass() / (l.cell().volume() - VOLUME);
    println!(
        "{} of 192 octanols kept; the lattice at {:.6} g/cm³; the octanol at {:.10} g/cm³",
        l.octanols,
        lattice / 1e3,
        rho / 1e3
    );
    assert!(l.octanols < 192, "phenol took somebody's place");
    assert!((rho - DENSITY).abs() <= 1e3 * EPS * DENSITY, "{rho}");
    let (plain, _, _) = build(false);
    let plain_rho = plain.octanols as f64 * o.mass() / (plain.cell().volume() - VOLUME);
    println!(
        "uncompensated: {} kept, {:+.2}%",
        plain.octanols,
        100.0 * (plain_rho / DENSITY - 1.0)
    );
    assert!((plain_rho - DENSITY).abs() > 1e-3 * DENSITY);
}

/// **The end states are the systems they say**, in both small boxes, at a configuration
/// constrained dynamics left: fully coupled, the whole field less phenol alone in the box plus
/// phenol alone in vacuum with its own terms; fully decoupled, the solvent alone in its box plus
/// phenol in vacuum, energy and every force. Each within `10⁻¹²` of the parts' magnitudes: the same
/// terms summed in other orders. Benzene's check, on the third solute.
#[test]
fn the_end_states_are_the_systems_they_say() {
    for (name, s) in [
        ("octanol", small_octanol(0xE5D).0),
        ("water", small_water(0xE5E)),
    ] {
        let at = moved(&s, Lambda::new(0.0, 0.5, 0.7), 20);
        let n = at.len();
        let mut f = vec![[0.0; 3]; n];
        let coupled = s.decoupling.energy_and_forces(&at, Lambda::COUPLED, &mut f);
        let whole = s.field.energy(&at);
        let in_box = s.solute_in_the_box.energy(&at[..N]);
        let mut vf = vec![[0.0; 3]; N];
        let vac = s.solute_in_vacuum.energy_and_forces(&at[..N], &mut vf);
        let expected = whole.total - in_box.total + vac;
        let scale = magnitude(&whole);
        println!("{name}: coupled {coupled:.10e} J against {expected:.10e}");
        assert!(
            (coupled - expected).abs() <= 1e-12 * scale,
            "{name} coupled"
        );
        let decoupled = s
            .decoupling
            .energy_and_forces(&at, Lambda::new(0.0, 0.0, 0.0), &mut f);
        let solvent = s.solvent_alone.evaluate(&at[N..]);
        let expected = solvent.energy.total + vac;
        println!("{name}: decoupled {decoupled:.10e} J against {expected:.10e}");
        assert!(
            (decoupled - expected).abs() <= 1e-12 * scale,
            "{name} decoupled"
        );
        let fscale = largest_force(&f);
        for i in 0..n {
            let want = if i < N { vf[i] } else { solvent.forces[i - N] };
            for k in 0..3 {
                assert!(
                    (f[i][k] - want[k]).abs() <= 1e-12 * fscale,
                    "{name}: atom {i} axis {k}: {:e} against {:e}",
                    f[i][k],
                    want[k]
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// The measurements: written and costed, not run here

/// The production's legs around phenol.
fn legs() -> [System; 2] {
    let (b, m) = phenol();
    [water_leg_of(&b, &m), octanol_leg_of(&b, &m, VOLUME)]
}

/// One leg at protocol `base` in `ROOT/NAME/TAG/`: benzene's `run_leg`, with phenol named.
fn run_leg(root: &std::path::Path, sys: &System, tag: &str, base: Protocol) -> LegResult {
    let leg_dir = root.join(sys.name);
    std::fs::create_dir_all(&leg_dir).expect("the leg's directory");
    let _lock = Lock::take(&leg_dir, tag);
    let start = melted(sys, &leg_dir, MELT_SEED);
    let dir = leg_dir.join(tag);
    std::fs::create_dir_all(&dir).expect("the windows' directory");
    let p = seeded(base, sys.name);
    let kt = p.kt();
    let leg = Leg::new(sys.name, &sys.decoupling, &sys.template, start, p, grid());
    log(
        &dir,
        sys.name,
        &format!(
            "{}: phenol and {} solvent molecules, {} atoms, box {:?} Å, protocol {p:?}, model \
             {:016x}",
            sys.name,
            sys.solvent,
            sys.built.len(),
            sys.field.cell().lengths().map(|x| x / ANGSTROM),
            leg.model
        ),
    );
    let mut schedule: Vec<usize> = SCHEDULE.to_vec();
    let mut records: Vec<Recorded> = schedule.iter().map(|&g| window(&leg, &dir, g)).collect();
    let found = refine_schedule(
        &mut schedule,
        &mut records,
        OVERLAP_THRESHOLD,
        INSERTION_ROUNDS,
        |a, b| interval_bar(a, b, kt).overlap,
        |g| window(&leg, &dir, g),
    );
    for f in &found {
        log(
            &dir,
            sys.name,
            &format!(
                "round {}: overlap {:.3} between candidates {} and {}: {}",
                f.round,
                f.overlap,
                f.between[0],
                f.between[1],
                f.inserted
                    .map_or("nothing between them".to_string(), |g| format!(
                        "inserted {g}"
                    ))
            ),
        );
    }
    report(&dir, sys.name, &records, kt)
}

/// `RT ln 10`, kcal/mol, at 298.15 K.
fn rt_ln10() -> f64 {
    kcal(BOLTZMANN.to_si() * KELVIN) * std::f64::consts::LN_10
}

/// **Both legs at `base` and phenol's log P**: the water first, then the octanol, every window
/// resumable; then log P from the two direct legs against 1.46, and `Δ log P` from benzene's
/// against −0.67. One seed: the gate reads [`the_charges_reseeded_measured`]'s.
fn partition_measured(tag: &str, base: Protocol) {
    let root = results_dir();
    let t = std::time::Instant::now();
    let [water, oct] = legs();
    let water = run_leg(&root, &water, tag, base);
    let oct = run_leg(&root, &oct, tag, base);
    let (g_hyd, e_hyd) = (-water.bar.0, water.bar.1);
    let (g_oct, e_oct) = (-oct.bar.0, oct.bar.1);
    let log_p = -(g_oct - g_hyd) / rt_ln10();
    let err = (e_oct * e_oct + e_hyd * e_hyd).sqrt() / rt_ln10();
    let lines = [
        format!(
            "ΔG_hyd {g_hyd:+.3} ± {e_hyd:.3} kcal/mol (TI {:+.3}, halves {:+.3}, {:+.3}, lowest \
             overlap {:.3}); FreeSolv's GAFF/AM1-BCC −5.71 ± 0.03, experiment −6.60 ± 0.20",
            -water.ti.0, -water.halves[0], -water.halves[1], water.overlap
        ),
        format!(
            "ΔG_solv,oct {g_oct:+.3} ± {e_oct:.3} kcal/mol (TI {:+.3}, halves {:+.3}, {:+.3}, \
             lowest overlap {:.3})",
            -oct.ti.0, -oct.halves[0], -oct.halves[1], oct.overlap
        ),
        format!(
            "log P {log_p:.2} ± {err:.2} from the two direct legs, one seed; experiment \
             {LOG_P_PHENOL} (Hansch, Leo and Hoekman 1995, via PubChem CID 996); Δ log P from \
             benzene's {:.2} ± {:.2}: {:+.2} ± {:.2} against experiment's {:+.2}; dry octanol; {:.0} s",
            BENZENE_LOG_P.0,
            BENZENE_LOG_P.1,
            log_p - BENZENE_LOG_P.0,
            (err * err + BENZENE_LOG_P.1 * BENZENE_LOG_P.1).sqrt(),
            LOG_P_PHENOL - LOG_P,
            t.elapsed().as_secs_f64()
        ),
    ];
    for line in lines {
        log(&root, tag, &line);
    }
}

/// **Phenol's log P at W3's protocol** (release, one core). Both legs, each window 20 ps discarded
/// and 200 ps sampled ([`LONG`]), into `water/long/` and `octanol/long/`.
/// `PANTOMETRY_OCTANOL_DIR=C:/pantometry-w4/phenol cargo test --release --locked -p
/// pantometry-forcefield --test phenol_partitioned_into_octanol the_partition_measured --
/// --ignored --nocapture`.
#[test]
#[ignore = "both legs at W3's protocol, thirteen windows each: about a day with --release"]
fn the_partition_measured() {
    partition_measured("long", LONG);
}

/// Refuses unless `schedule` holds every window of [`SCHEDULE`] and ends at the decoupled end.
fn assert_whole(schedule: &[usize], leg: &str) {
    let last = grid().len() - 1;
    assert!(
        SCHEDULE.iter().all(|g| schedule.contains(g)) && schedule.last() == Some(&last),
        "{leg} is not a finished leg: it has {schedule:?}, and a finished one has every one of \
         {SCHEDULE:?} and ends at {last}"
    );
}

/// **A leg killed partway is not whole**: [`assert_whole`] passes the schedule a finished leg
/// leaves, with an inserted window too, and refuses by name one missing any of [`SCHEDULE`]'s
/// windows, one stopped before the decoupled end, and one that is only the charge segment.
#[test]
#[cfg(not(target_family = "wasm"))]
fn a_leg_killed_partway_is_not_whole() {
    let full: Vec<usize> = SCHEDULE.to_vec();
    assert_whole(&full, "full");
    let mut inserted = full.clone();
    inserted.insert(3, 9);
    assert_whole(&inserted, "inserted");
    let mut missing = full.clone();
    missing.remove(6);
    let short: Vec<usize> = full[..full.len() - 1].to_vec();
    for (name, s) in [
        ("missing", missing),
        ("short", short),
        ("charges", charge_segment(&SCHEDULE)),
    ] {
        let said = refusal(|| assert_whole(&s, name));
        assert!(said.contains("is not a finished leg"), "{name}: {said}");
    }
}

/// The charge segment's last candidate: the electrostatics off and the van der Waals whole.
const LAST_CHARGE: usize = 8;

/// **A leg's charge segment, read from its schedule**: the candidates up to and including
/// [`LAST_CHARGE`]. Benzene's and p-xylene's files write it as `[0, 4, 8]`, which an insertion
/// between 0 and 8 makes wrong, and the reseed would then refuse a day into the run, over a
/// schedule the leg chose for itself. Phenol's charge segment is where an insertion is likeliest.
fn charge_segment(schedule: &[usize]) -> Vec<usize> {
    schedule
        .iter()
        .copied()
        .take_while(|&g| g <= LAST_CHARGE)
        .collect()
}

/// **The charge segment is read from the schedule**: every candidate up to [`LAST_CHARGE`] has the
/// van der Waals whole and every one after it the electrostatics off and the van der Waals less
/// than whole, so the segment is all charge and the rest none; [`SCHEDULE`]'s segment is 0, 4 and
/// 8, and the same schedule with windows inserted at 2, 6 and 9 has 0, 2, 4, 6 and 8.
#[test]
fn the_charge_segment_is_read_from_the_schedule() {
    let g = grid();
    for (i, l) in g.iter().enumerate() {
        if i <= LAST_CHARGE {
            assert_eq!(l.van_der_waals, 1.0, "candidate {i}");
        } else {
            assert_eq!(l.electrostatics, 0.0, "candidate {i}");
            assert!(l.van_der_waals < 1.0, "candidate {i}");
        }
    }
    assert_eq!(g[0].electrostatics, 1.0);
    assert_eq!(g[LAST_CHARGE].electrostatics, 0.0);
    assert_eq!(charge_segment(&SCHEDULE), [0, 4, 8]);
    let mut inserted = SCHEDULE.to_vec();
    for c in [2, 6, 9] {
        let at = inserted
            .iter()
            .position(|&x| x > c)
            .expect("a later candidate");
        inserted.insert(at, c);
    }
    assert_eq!(charge_segment(&inserted), [0, 2, 4, 6, 8]);
}

/// **The gate's ceiling is this method's benzene, not experiment's**: 1.92 − 1.46 = 0.46, below
/// experiment's 2.13 − 1.46 = 0.67.
#[test]
fn the_gate_ceiling_is_the_methods_benzene() {
    assert!((ceiling() - 0.46).abs() < 1e-12, "{}", ceiling());
    assert!(ceiling() < LOG_P - LOG_P_PHENOL);
}

/// The octanol leg's reseeded windows' seed: benzene's reseed's.
const RESEED: u64 = 0x4C2_5EED;

/// The water leg's: neither leg's `long/` seed and not the octanol's reseed.
const WATER_RESEED: u64 = 0x4C2_5EEE;

/// A leg's `long/`, checked whole and read back, with its directory locked for the reseed.
struct Finished<'a> {
    sys: &'a System,
    leg_dir: std::path::PathBuf,
    start: Vec<[f64; 3]>,
    long_p: Protocol,
    records: Vec<Recorded>,
    /// How many of `records` are the charge segment ([`charge_segment`]).
    charges: usize,
    _lock: Lock,
}

/// **`sys`'s `long/`, refused unless finished**: no `long.lock` left by a run still inserting
/// windows or killed while it did, every one of [`SCHEDULE`]'s windows, the decoupled end, and every
/// file read back against its header.
fn finished<'a>(root: &std::path::Path, sys: &'a System) -> Finished<'a> {
    let leg_dir = root.join(sys.name);
    // A leg in its insertion rounds has all thirteen windows and passes `assert_whole`, and the
    // reseed takes its own lock, so the two would not exclude each other. `run_leg` holds
    // `long.lock` until it returns, and a run killed partway leaves it: either way, not finished.
    let lock = leg_dir.join("long.lock");
    assert!(
        !lock.exists(),
        "{} exists: that leg is still running or was killed, and is not a finished leg",
        lock.display()
    );
    let _lock = Lock::take(&leg_dir, "reseed");
    let start = read_back_start(sys, &leg_dir);
    let long_p = seeded(LONG, sys.name);
    let long = Leg::new(
        sys.name,
        &sys.decoupling,
        &sys.template,
        start.clone(),
        long_p,
        grid(),
    );
    let long_dir = leg_dir.join("long");
    let schedule = finished_schedule(&long_dir);
    assert_whole(&schedule, &format!("{}/long/", sys.name));
    let charges = charge_segment(&schedule).len();
    assert_eq!(
        schedule[charges - 1],
        LAST_CHARGE,
        "{}/long/'s schedule {schedule:?} has no window at the charge segment's end",
        sys.name
    );
    let records = schedule
        .iter()
        .map(|&g| read_back(&long, &long_dir, g))
        .collect();
    Finished {
        sys,
        leg_dir,
        start,
        long_p,
        records,
        charges,
        _lock,
    }
}

/// `f`'s charge segment again at `seed`, into `NAME/reseed/`, each window resumable: the same
/// leg, model, start, protocol and windows as `long/`, only the seed changed.
fn reseeded(f: &Finished, seed: u64) -> Vec<Recorded> {
    assert_ne!(f.long_p.seed, seed, "{}: a second seed", f.sys.name);
    let dir = f.leg_dir.join("reseed");
    std::fs::create_dir_all(&dir).expect("the reseed's directory");
    let p = Protocol { seed, ..f.long_p };
    let leg = Leg::new(
        f.sys.name,
        &f.sys.decoupling,
        &f.sys.template,
        f.start.clone(),
        p,
        grid(),
    );
    let windows: Vec<usize> = f.records[..f.charges].iter().map(|r| r.grid).collect();
    log(
        &dir,
        "reseed",
        &format!(
            "{}: the charge segment {windows:?} reseeded, {:#x} for long/'s {:#x}; protocol \
             {p:?}, model {:016x}",
            f.sys.name, p.seed, f.long_p.seed, leg.model
        ),
    );
    windows.iter().map(|&g| window(&leg, &dir, g)).collect()
}

/// **`f`'s ΔG of decoupling with both seeds of its charge segment**, kcal/mol, and its error,
/// logged in `NAME/reseed/`: the segment `a` from `long/` and `b` reseeded, taken as their mean
/// with the error `max(|d|/2, √(e_a² + e_b²)/2)` for their difference `d`, which is how benzene's
/// 1.92 ± 0.23 and p-xylene's 3.06 ± 0.22 were taken; and the van der Waals segment from `long/`.
fn with_both_seeds(f: &Finished, again: &[Recorded], seed: u64) -> (f64, f64) {
    let kt = f.long_p.kt();
    let k = f.charges;
    let (long_charges, long_vdw) = (&f.records[..k], &f.records[k - 1..]);
    let (whole, charges, vdw) = (
        chain(&f.records, kt),
        chain(long_charges, kt),
        chain(long_vdw, kt),
    );
    assert!(
        (whole.value - (charges.value + vdw.value)).abs()
            <= 64.0 * EPS * (charges.value.abs() + vdw.value.abs()),
        "{}: the segments add to the leg",
        f.sys.name
    );
    let (a, ea) = segment(long_charges, kt);
    let (b, eb) = segment(again, kt);
    let (v, ev) = segment(long_vdw, kt);
    let d = b - a;
    let combined = (ea * ea + eb * eb).sqrt();
    let mean = 0.5 * (a + b);
    let e_mean = (0.5 * d.abs()).max(0.5 * combined);
    let dir = f.leg_dir.join("reseed");
    for (name, rs, seed) in [
        ("long/", long_charges, f.long_p.seed),
        ("reseed/", again, seed),
    ] {
        let (g, e) = segment(rs, kt);
        let h = halves_of(rs, kt);
        let taus: Vec<String> = rs
            .iter()
            .map(|r| {
                let tau = Estimate::of(&r.gradients.iter().map(|x| x[1]).collect::<Vec<_>>()).tau;
                format!("{} τ {tau:.1}", r.grid)
            })
            .collect();
        log(
            &dir,
            "reseed",
            &format!(
                "{} {name} (seed {seed:#x}): the charge segment's ΔG {g:+.3} ± {e:.3} kcal/mol, \
                 halves {:+.3}, {:+.3}; {}",
                f.sys.name,
                h[0],
                h[1],
                taus.join(", ")
            ),
        );
    }
    log(
        &dir,
        "reseed",
        &format!(
            "{}: the seeds differ by {d:+.3} kcal/mol, {:.2} of their combined quoted error \
             {combined:.3}; with both, the charge segment {mean:+.3} ± {e_mean:.3} and the van der \
             Waals segment {v:+.3} ± {ev:.3}",
            f.sys.name,
            d.abs() / combined
        ),
    );
    (mean + v, (e_mean * e_mean + ev * ev).sqrt())
}

/// **Both legs' charge segments at a second seed, and log P with both, against the gate**
/// (release, one core). p-xylene's reseed, on phenol, and on the water leg too: phenol's hydration
/// free energy is mostly its charge segment, where benzene's and p-xylene's were mostly not, so
/// that segment's seed spread is in `σ` rather than beside it. Before any window runs, **both legs
/// are checked finished** ([`finished`]), so neither's hours are spent on a run the other would
/// refuse. Then each leg's charge segment, read from its own schedule ([`charge_segment`]), is run
/// again at its own seed into `NAME/reseed/`; each leg is taken with both seeds
/// ([`with_both_seeds`]); and log P and its `σ` are read against the gate in the module
/// documentation, with whether the result is distinguishable from benzene's.
///
/// Nothing is asserted but that the segments add to each leg.
#[test]
#[ignore = "three octanol windows and three water windows at W3's protocol: about 7 h with --release"]
fn the_charges_reseeded_measured() {
    let root = results_dir();
    let t = std::time::Instant::now();
    let [water, oct] = legs();
    let (fw, fo) = (finished(&root, &water), finished(&root, &oct));
    let (water_again, oct_again) = (reseeded(&fw, WATER_RESEED), reseeded(&fo, RESEED));
    let (gw, ew) = with_both_seeds(&fw, &water_again, WATER_RESEED);
    let (go, eo) = with_both_seeds(&fo, &oct_again, RESEED);
    let (g_hyd, g_oct) = (-gw, -go);
    let log_p = -(g_oct - g_hyd) / rt_ln10();
    let sigma = (eo * eo + ew * ew).sqrt() / rt_ln10();
    let off = (log_p - LOG_P_PHENOL).abs();
    let verdict = if 2.0 * sigma >= ceiling() {
        "uninformative: 2σ reaches this method's own benzene from experiment's phenol"
    } else if off <= 2.0 * sigma {
        "the method carries to phenol"
    } else if off > 3.0 * sigma {
        "the method does not carry to phenol"
    } else {
        "not settled"
    };
    let delta = log_p - BENZENE_LOG_P.0;
    let e_delta = (sigma * sigma + BENZENE_LOG_P.1 * BENZENE_LOG_P.1).sqrt();
    let apart = if -delta > 2.0 * e_delta {
        "below benzene by more than 2σ: the hydroxyl is seen"
    } else if delta.abs() <= 2.0 * e_delta {
        "within 2σ of benzene: not distinguished from it"
    } else {
        "above benzene by more than 2σ"
    };
    let lines = [
        format!(
            "ΔG_solv,oct {g_oct:+.3} ± {eo:.3}, ΔG_hyd {g_hyd:+.3} ± {ew:.3} kcal/mol, each with \
             both seeds of its charge segment; FreeSolv's ΔG_hyd for GAFF/AM1-BCC −5.71 ± 0.03, \
             experiment −6.60 ± 0.20"
        ),
        format!(
            "log P {log_p:.2} ± {sigma:.2} against {LOG_P_PHENOL}, {:.2} σ, 2σ {:.2} against the \
             ceiling {:.2}: {verdict}. Not in σ: the van der Waals segments' seed spread, \
             unmeasured; a two-seed difference's one degree of freedom; experiment's 1.46 against \
             Sangster's 1.6",
            off / sigma,
            2.0 * sigma,
            ceiling()
        ),
        format!(
            "Δ log P from benzene's {:.2} ± {:.2}: {delta:+.2} ± {e_delta:.2} against experiment's \
             {:+.2}, {apart}; reported and not gated; {:.0} s",
            BENZENE_LOG_P.0,
            BENZENE_LOG_P.1,
            LOG_P_PHENOL - LOG_P,
            t.elapsed().as_secs_f64()
        ),
    ];
    for line in lines {
        log(&root, "reseed", &line);
    }
}

/// **The melts and the pilot** (release, one core, in the background): benzene's
/// `the_pilot_measured` on phenol. Each leg is melted into `NAME/start.txt`, the production's own
/// start, and three windows of [`PILOT`] run at candidates 0, 8 and 22 in `NAME/pilot/`, each
/// one's τ in samples of 0.1 ps.
#[test]
#[ignore = "both boxes melted and three 32 ps windows of each: a few hours with --release"]
fn the_pilot_measured() {
    let root = results_dir();
    for sys in legs() {
        let leg_dir = root.join(sys.name);
        std::fs::create_dir_all(&leg_dir).expect("the leg's directory");
        let _lock = Lock::take(&leg_dir, "pilot");
        let start = melted(&sys, &leg_dir, MELT_SEED);
        let dir = leg_dir.join("pilot");
        std::fs::create_dir_all(&dir).expect("the pilot's directory");
        let p = seeded(PILOT, sys.name);
        let leg = Leg::new(sys.name, &sys.decoupling, &sys.template, start, p, grid());
        log(
            &dir,
            "pilot",
            &format!(
                "\n{}: phenol and {} solvent molecules, {} atoms; τ in samples of 0.1 ps\n| \
                 candidate | λ_e | λ_v | ⟨∂U/∂λ⟩ (kcal/mol) | its τ | ms a step |\n| --- | --- | \
                 --- | --- | --- | --- |",
                sys.name,
                sys.solvent,
                sys.built.len()
            ),
        );
        for g in [0usize, 8, 22] {
            let r = window(&leg, &dir, g);
            let c = if g == 0 { 1 } else { 2 };
            let grad = Estimate::of(&r.gradients.iter().map(|x| x[c]).collect::<Vec<_>>());
            log(
                &dir,
                "pilot",
                &format!(
                    "| {g} | {} | {} | {:+.3} ± {:.3} (∂/∂λ_{}) | {:.1} | {:.1} |",
                    r.lambda.electrostatics,
                    r.lambda.van_der_waals,
                    kcal(grad.mean),
                    kcal(grad.error),
                    if c == 1 { "e" } else { "v" },
                    grad.tau,
                    1e3 * r.seconds / p.steps_per_window() as f64
                ),
            );
        }
    }
}
