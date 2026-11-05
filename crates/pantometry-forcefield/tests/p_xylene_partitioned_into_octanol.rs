//! **p-xylene decoupled from 1-octanol and from water, for its log P**: the second solute through
//! 4c-2's legs. Benzene's log P came out at 1.92 ± 0.23 from the two direct legs against
//! experiment's 2.13 (`benzene_partitioned_into_octanol.rs`), which is one molecule; this asks
//! whether the method, not the molecule, is what agreed. p-xylene is benzene with two methyl
//! groups across the ring, so it carries the two types benzene does not (`c3`, `hc`) and none that
//! octanol does not already have here.
//!
//! **Kept out of the published crate** (`Cargo.toml`'s `exclude`), as benzene's is: GAFF's values
//! are read from `gaff.dat`, which is GPL v3.
//!
//! # The model and the boxes
//!
//! Everything is benzene's, built by the same lines (`tests/partition/`'s `_of` builders):
//!
//! - **The solute**: `PXY` of the Chemical Component Dictionary, which `the_congener_series.rs`
//!   already holds to the molecule its entry names. It is typed for UFF, relaxed in UFF's vacuum
//!   with its AM1-BCC charges, and every X–H is held.
//! - **Its non-bonded terms**: GAFF's `c3`, `hc`, `ca` and `ha` and FreeSolv's AM1-BCC charges
//!   (`mobley_3234716.mol2`, `tests/freesolv/`), mapped by graph isomorphism.
//! - **The rules**: Lorentz–Berthelot and GAFF's 1-4 scales.
//! - **The octanol box**: 4c-1's 192-site lattice, compensated for the solute's volume. That
//!   volume is liquid p-xylene's molar volume, 106.16 / 0.86104 = 123.3 cm³/mol ([`VOLUME`]).
//! - **The water box**: W3's.
//! - **The melts, the 29 candidates, the 13 windows, the insertion, the protocol and the seeds**:
//!   benzene's. The results go under their own `PANTOMETRY_OCTANOL_DIR`. Every window and start
//!   carries a digest of the model, so nothing of benzene's can be read back as p-xylene's.
//!
//! [`benzenes_legs_are_the_bits_they_were`] holds benzene's two production boxes and models, and
//! the box benzene's `octanol_box` builds for its other tests, to the digests benzene's own runs
//! wrote. So the generalisation moved none of the bits those builders give benzene. Benzene's small
//! boxes, built by the same wrapper at other sizes, are not pinned; the wrapper forwards every
//! argument in order, and that is what pinning its production-size box checks.
//!
//! # What it is compared with, registered before the run
//!
//! Experiment's `log K_ow` is **3.15**: C. Hansch, A. Leo and D. Hoekman, *Exploring QSAR:
//! Hydrophobic, Electronic and Steric Constants*, ACS 1995, p. 43. It was read on 2026-11-02 as
//! PubChem CID 7809's HSDB record quotes it, the same source benzene's 2.13 is from. The record's
//! ICSC entry is 3.15 too, and its Sangster entry 3.42. FreeSolv's hydration free energy is −0.80
//! by experiment and −0.66 ± 0.03 for GAFF/AM1-BCC on GAFF's bonded terms.
//!
//! **The gate, written before any window ran.** Let `σ` be log P's error with the octanol charge
//! segment's seed spread in it, as benzene's two-seed number is. Then:
//! - the verdict is **uninformative** if `2σ ≥ 1.02`, whatever log P is: that is experiment's
//!   p-xylene less benzene, so a band that wide passes a method blind to the methyls. σ takes the
//!   larger of the seeds' half-difference and BAR's, so it widens as the seeds disagree, and this is
//!   its ceiling;
//! - otherwise the method **carries to p-xylene** if `|log P − 3.15| ≤ 2σ`;
//! - it **does not** if `|log P − 3.15| > 3σ`;
//! - between the two, it is not settled.
//!
//! σ leaves out three things, and they are stated beside the verdict rather than added to it. The
//! van der Waals segment's own seed spread is not measured. A difference of two seeds is one degree
//! of freedom. And experiment's own sources differ by 0.27 (3.15 against Sangster's 3.42).
//!
//! The difference from benzene, `Δ log P = log P(p-xylene) − log P(benzene)`, is reported beside
//! it against experiment's 1.02, the same reading on the same two legs. Dry octanol's and the hybrid
//! model's errors should partly cancel from it, so it is the sharper of the two, and it is reported,
//! not gated. Nothing in this file asserts either number.
//!
//! # What it measured
//!
//! Run on `x86_64-pc-windows-gnu` from the binary built at `6a4f1db`. The two legs took 107 542 s
//! and finished on 2026-11-03. The second seed took 15 499 s and finished on 2026-11-04. Its first
//! attempt died in its first window on a failed allocation: the machine had reached its commit
//! limit, with 11 GB in the kernel's paged pool. It was rerun from that window after a reboot, and
//! nothing it had written was kept. No window was inserted on either leg. In kcal/mol:
//!
//! | | BAR | TI | halves | lowest overlap |
//! | --- | --- | --- | --- | --- |
//! | water, decoupled | +1.099 ± 0.188 | +1.273 ± 0.183 | +1.226, +0.956 | 0.198 |
//! | octanol, decoupled, one seed | +5.364 ± 0.284 | +5.583 ± 0.293 | +5.046, +5.757 | 0.243 |
//!
//! - **Hydration:** ΔG_hyd is −1.099 ± 0.188, against FreeSolv's −0.66 ± 0.03 for GAFF/AM1-BCC on
//!   GAFF's bonded terms and −0.80 ± 0.60 by experiment.
//! - **One seed:** log P is 3.13 ± 0.25.
//! - **The charge segment:** +0.726 ± 0.190 at `long/`'s seed, whose candidate 0 has τ = 145
//!   samples, so about 14 independent ones. It is +0.541 ± 0.138 reseeded. The seeds differ by
//!   −0.185, 0.79 of their combined quoted error. So the segment's error is BAR's, ±0.117, not
//!   the half-difference.
//! - **Both seeds:** the charge segment is +0.633 ± 0.117 and the van der Waals segment
//!   +4.638 ± 0.210. ΔG_solv,oct is −5.271 ± 0.241.
//! - **log P is 3.06 ± 0.22 against 3.15**, a miss of 0.41σ. 2σ = 0.44 is under the 1.02 that would
//!   make it uninformative, so **the method carries to p-xylene**.
//! - **Δ log P** from benzene's 1.92 ± 0.23 is +1.14 ± 0.32 against experiment's +1.02. It is
//!   reported, not gated.
//!
//! σ still leaves out the three things above. Sangster's 3.42 would be a miss of 1.6σ, which the
//! gate still calls carried.

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

/// p-xylene, `PXY` of the Chemical Component Dictionary, byte for byte.
const PXY: &str = include_str!("../components/PXY.cif");

/// Its atoms.
const N: usize = 18;

/// **p-xylene's volume in the solvent, m³ a molecule**: liquid p-xylene's molar volume, 106.16
/// g/mol over 0.86104 g/cm³ at 20 °C (O'Neil, *The Merck Index*, RSC 2013, p. 1876, as PubChem CID
/// 7809's HSDB record quotes it, read 2026-11-02), 123.3 cm³/mol. It stands in for the partial molar
/// volume in octanol, as benzene's 89 does. The same computation gives benzene 89.2 from its own
/// record's 0.8756.
const VOLUME: f64 = 106.16 / 0.86104 * 1e-6 / AVOGADRO;

/// Experiment's `log K_ow` for p-xylene: see the module documentation.
const LOG_P_XYLENE: f64 = 3.15;

/// Benzene's log P from its two direct legs over two seeds of the octanol charge segment, and its
/// error (`the_octanol_charges_reseeded_measured`, run 2026-10-24): what [`Δ log
/// P`](partition_measured) is taken from.
const BENZENE_LOG_P: (f64, f64) = (1.92, 0.23);

/// p-xylene with FreeSolv's AM1-BCC charges, relaxed in UFF's vacuum with them and put on its
/// constraints; and its model.
fn p_xylene() -> (Benzene, Model) {
    let c = Component::from_ccd(PXY).expect("PXY");
    let model = Model::p_xylene(&c);
    (Benzene::from_ccd(PXY, Some(model.charges.clone())), model)
}

/// Each atom's GAFF type by its place in the molecule, not by FreeSolv's file: a carbon with four
/// neighbours is a methyl's `c3` and its hydrogens `hc`; a ring carbon is `ca` and its hydrogen
/// `ha`.
fn types_by_structure(c: &Component) -> Vec<&'static str> {
    let n = c.atoms().len();
    let mut degree = vec![0usize; n];
    let mut partner = vec![usize::MAX; n];
    for b in c.bonds() {
        for (i, j) in [(b.atoms[0], b.atoms[1]), (b.atoms[1], b.atoms[0])] {
            degree[i] += 1;
            partner[i] = j;
        }
    }
    (0..n)
        .map(|i| match c.atoms()[i].element {
            Element::C if degree[i] == 4 => "c3",
            Element::C => "ca",
            Element::H if degree[partner[i]] == 4 => "hc",
            Element::H => "ha",
            other => panic!("p-xylene has no {other:?}"),
        })
        .collect()
}

/// Every pair of `c`'s atoms by its bond distance: `(pairs three bonds apart, pairs more than
/// three apart)`, by breadth-first search on the bonds.
fn pairs_by_distance(c: &Component) -> (usize, usize) {
    let n = c.atoms().len();
    let mut next = vec![Vec::new(); n];
    for b in c.bonds() {
        next[b.atoms[0]].push(b.atoms[1]);
        next[b.atoms[1]].push(b.atoms[0]);
    }
    let (mut three, mut beyond) = (0, 0);
    for i in 0..n {
        let mut d = vec![usize::MAX; n];
        d[i] = 0;
        let mut queue = std::collections::VecDeque::from([i]);
        while let Some(a) = queue.pop_front() {
            for &b in &next[a] {
                if d[b] == usize::MAX {
                    d[b] = d[a] + 1;
                    queue.push_back(b);
                }
            }
        }
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

/// **Benzene's legs are the bits they were**: the production's octanol and water boxes, as built,
/// and their models, against the digests benzene's own melts wrote into `start.txt` on 2026-10-21
/// (`C:\pantometry-w4\octanol\{octanol,water}\start.txt`). Benzene's builders are the `_of`
/// builders with benzene's model and volume now, and this is what says that changed nothing.
///
/// **On the platform those runs were made on only**, `x86_64-pc-windows-gnu`. The boxes are not
/// the same bits elsewhere, and `tests/pinned/` says why: a minimiser and the platform's `cos` are
/// in them. W3's lattice turns each water by `Rng::gaussian`, which is the platform's `ln` and
/// `cos`, about SETTLE's frame, which is its `cos`. CI measured this at `6a4f1db`:
/// - the water box as built was `ba16a354…` on Windows (MSVC) and `47588d36…` on macOS, where it is
///   `a16ebf27…` here and on Linux;
/// - the octanol box and both models were the same bits on all three;
/// - WebAssembly failed one of the comparisons, and its output did not say which.
///
/// A pin that holds only where the runs were made is what the claim needs, since the claim is
/// about those runs' bits.
#[test]
#[cfg(all(target_os = "windows", target_env = "gnu"))]
fn benzenes_legs_are_the_bits_they_were() {
    for (sys, atoms, built, model) in [
        (
            octanol_leg(),
            5088,
            0x6a00_3b69_8247_a78a_u64,
            0xbfd0_ccdb_52ef_40bc_u64,
        ),
        (
            water_leg(),
            1533,
            0xa16e_bf27_64c1_48cf,
            0xb80d_6cc9_c032_7aff,
        ),
    ] {
        println!(
            "{}: {} atoms, built {:016x}, model {:016x}",
            sys.name,
            sys.built.len(),
            positions_digest(&sys.built),
            model_digest(&sys.field)
        );
        assert_eq!(sys.built.len(), atoms, "{}", sys.name);
        assert_eq!(positions_digest(&sys.built), built, "{} as built", sys.name);
        assert_eq!(model_digest(&sys.field), model, "{}'s model", sys.name);
    }
    // The wrapper benzene's other tests and the decane leg build through, at the production's
    // size: the same lattice as the octanol leg's.
    let (l, _, _) = octanol_box(
        &benzene(),
        &partition::octanol(),
        octanol::LIQUID,
        9.0 * ANGSTROM,
        1e-5,
        0x0C71,
        true,
    );
    assert_eq!(
        positions_digest(&l.at),
        0x6a00_3b69_8247_a78a,
        "octanol_box"
    );
}

/// The small boxes, as benzene's file builds them, around p-xylene: the field, the decoupling, the
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
    let (b, m) = p_xylene();
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
        decoupling: PeriodicDecoupling::new(&field, &group).expect("p-xylene decouples"),
        solute_in_the_box: rules(
            PeriodicForceField::new(&b.component, &b.types, cell, p).expect("p-xylene"),
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

/// p-xylene in 64 lattice waters less those within 2.3 Å of it, 12.4 Å across, `r_c` = 6.0 Å, its
/// model swapped in. Benzene's small water box is 27 waters at 4.6 Å; p-xylene is 2 Å longer.
fn small_water(seed: u64) -> Small {
    let (b, m) = p_xylene();
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
            PeriodicForceField::new(&b.component, &b.types, cell, p).expect("p-xylene"),
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
///   ([`types_by_structure`]), not by FreeSolv's type column, and FreeSolv's types agree.
/// - The charges are FreeSolv's to the bit, **each by the atom's place**. A ring carbon bonded to a
///   methyl carries −0.0817 and any other −0.1273. One methyl carbon carries −0.0529 and the other
///   −0.0528. Methyl hydrogens carry +0.0432 and ring hydrogens +0.1298. They sum to +0.0001 e, and
///   none is QEq's.
/// - The mixing is Lorentz–Berthelot, with GAFF's 1-4 scales.
/// - The solute's own 1-4 pairs are the pairs three bonds apart, counted on the graph, and its
///   own pairs those and the ones further.
/// - **The same model is in both boxes.**
#[test]
fn the_model_is_what_it_says() {
    let (oct, _) = small_octanol(0xE5D);
    let wat = small_water(0xE5E);
    let (b, m) = p_xylene();
    assert_eq!(b.component.atoms().len(), N);
    let kinds = types_by_structure(&b.component);
    assert_eq!(
        m.types,
        kinds.iter().map(|k| k.to_string()).collect::<Vec<_>>(),
        "FreeSolv's types are the structure's"
    );
    let qeq = Benzene::from_ccd(PXY, None).charges;
    let (three, beyond) = pairs_by_distance(&b.component);
    println!("p-xylene: {three} pairs three bonds apart, {beyond} further");
    for (name, s) in [("octanol", &oct), ("water", &wat)] {
        let f = &s.field;
        let (vdw, q) = (f.van_der_waals(), f.charges());
        for i in 0..N {
            assert_eq!(
                vdw[i].map(f64::to_bits),
                gaff(kinds[i]).map(f64::to_bits),
                "{name} {i}"
            );
            let ipso =
                b.component.bonds().iter().any(|bd| {
                    bd.atoms.contains(&i) && kinds[bd.atoms[0] + bd.atoms[1] - i] == "c3"
                });
            let allowed: &[f64] = match kinds[i] {
                "ca" if ipso => &[-0.0817],
                "ca" => &[-0.1273],
                "c3" => &[-0.0529, -0.0528],
                "hc" => &[0.0432],
                _ => &[0.1298],
            };
            assert!(
                allowed.contains(&q[i]),
                "{name}: atom {i} ({}{}) carries {} e",
                kinds[i],
                if ipso { ", bonded to a methyl" } else { "" },
                q[i]
            );
            assert_ne!(q[i].to_bits(), qeq[i].to_bits(), "{name}: QEq's charge");
        }
        let mut methyls: Vec<u64> = (0..N)
            .filter(|&i| kinds[i] == "c3")
            .map(|i| q[i].to_bits())
            .collect();
        methyls.sort_unstable();
        methyls.dedup();
        assert_eq!(
            methyls.len(),
            2,
            "{name}: one methyl at each of FreeSolv's two charges"
        );
        let net: f64 = q[..N].iter().sum();
        assert!((net - 0.0001).abs() < 1e-12, "{name}: p-xylene's {net} e");
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
            "{name}: p-xylene and {} solvent molecules, {} atoms",
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

/// **The octanol is at its density in the volume p-xylene leaves it**, in the production's box,
/// the fixed point reached; and the compensation does something.
#[test]
fn the_octanol_is_at_its_density_around_p_xylene() {
    let (b, m) = p_xylene();
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
    assert!(l.octanols < 192, "p-xylene took somebody's place");
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
/// constrained dynamics left: fully coupled, the whole field less p-xylene alone in the box plus
/// p-xylene alone in vacuum with its own terms; fully decoupled, the solvent alone in its box
/// plus p-xylene in vacuum, energy and every force. Each within `10⁻¹²` of the parts' magnitudes:
/// the same terms summed in other orders. Benzene's check, on the second solute.
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
// The measurements: written, costed, the first two run once

/// The production's legs around p-xylene.
fn legs() -> [System; 2] {
    let (b, m) = p_xylene();
    [water_leg_of(&b, &m), octanol_leg_of(&b, &m, VOLUME)]
}

/// One leg at protocol `base` in `ROOT/NAME/TAG/`: benzene's `run_leg`, with p-xylene named.
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
            "{}: p-xylene and {} solvent molecules, {} atoms, box {:?} Å, protocol {p:?}, model \
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

/// **Both legs at `base` and p-xylene's log P**: the water first, then the octanol, every window
/// resumable; then log P from the two direct legs against 3.15, and `Δ log P` from benzene's
/// against 1.02. One seed: the gate reads [`the_octanol_charges_reseeded_measured`]'s.
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
             overlap {:.3}); FreeSolv's GAFF/AM1-BCC −0.66 ± 0.03, experiment −0.80 ± 0.60",
            -water.ti.0, -water.halves[0], -water.halves[1], water.overlap
        ),
        format!(
            "ΔG_solv,oct {g_oct:+.3} ± {e_oct:.3} kcal/mol (TI {:+.3}, halves {:+.3}, {:+.3}, \
             lowest overlap {:.3})",
            -oct.ti.0, -oct.halves[0], -oct.halves[1], oct.overlap
        ),
        format!(
            "log P {log_p:.2} ± {err:.2} from the two direct legs, one seed; experiment \
             {LOG_P_XYLENE} (Hansch, Leo and Hoekman 1995, via PubChem CID 7809); Δ log P from \
             benzene's {:.2} ± {:.2}: {:+.2} ± {:.2} against experiment's {:+.2}; dry octanol; {:.0} s",
            BENZENE_LOG_P.0,
            BENZENE_LOG_P.1,
            log_p - BENZENE_LOG_P.0,
            (err * err + BENZENE_LOG_P.1 * BENZENE_LOG_P.1).sqrt(),
            LOG_P_XYLENE - LOG_P,
            t.elapsed().as_secs_f64()
        ),
    ];
    for line in lines {
        log(&root, tag, &line);
    }
}

/// **p-xylene's log P at W3's protocol** (release, one core; run once, in the module documentation). Both legs, each
/// window 20 ps discarded and 200 ps sampled ([`LONG`]), into `water/long/` and `octanol/long/`.
/// `PANTOMETRY_OCTANOL_DIR=C:/pantometry-w4/p-xylene cargo test --release --locked -p
/// pantometry-forcefield --test p_xylene_partitioned_into_octanol the_partition_measured --
/// --ignored --nocapture`.
#[test]
#[ignore = "both legs at W3's protocol, thirteen windows each: about a day with --release, run once (2026-11-03)"]
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
        ("charges", CHARGE_SEGMENT.to_vec()),
    ] {
        let said = refusal(|| assert_whole(&s, name));
        assert!(said.contains("is not a finished leg"), "{name}: {said}");
    }
}

/// The octanol leg's charge segment, candidates 0, 4 and 8: benzene's reseeded windows.
const CHARGE_SEGMENT: [usize; 3] = [0, 4, 8];

/// The reseeded windows' seed: benzene's reseed's.
const RESEED: u64 = 0x4C2_5EED;

/// **The octanol leg's charge segment at a second seed, and log P with both, against the gate**
/// (release, one core; run once, the verdict in the module documentation). Benzene's reseed, on
/// p-xylene: the same leg, model, start
/// and protocol as `octanol/long/`, only the seed changed, the windows at [`CHARGE_SEGMENT`] into
/// `octanol/reseed/`. Before any window runs, every file it reads back is checked against its
/// header, and **both legs are checked whole**: every one of [`SCHEDULE`]'s windows there, and the
/// decoupled end, candidate 28. A leg killed partway leaves whole windows that read back cleanly,
/// and without that a verdict would be printed for a decoupling that stopped partway. Then:
/// - the charge segment `a` from `long/` and `b` reseeded, and their difference `d`;
/// - the segment taken as their mean, its error `max(|d|/2, √(e_a² + e_b²)/2)`, which is how
///   benzene's 1.92 ± 0.23 was taken;
/// - log P and its `σ` with that, against the gate in the module documentation, and the verdict.
///
/// Nothing is asserted but that the segments add to the leg.
#[test]
#[ignore = "three octanol windows at W3's protocol: about 5 h with --release, run once (2026-11-04)"]
fn the_octanol_charges_reseeded_measured() {
    let root = results_dir();
    let t = std::time::Instant::now();
    let [water, sys] = legs();
    let leg_dir = root.join(sys.name);
    let (long_dir, dir) = (leg_dir.join("long"), leg_dir.join("reseed"));
    let _lock = Lock::take(&leg_dir, "reseed");
    let start = read_back_start(&sys, &leg_dir);
    let long_p = seeded(LONG, sys.name);
    assert_ne!(long_p.seed, RESEED, "a second seed");
    let kt = long_p.kt();
    let long = Leg::new(
        sys.name,
        &sys.decoupling,
        &sys.template,
        start.clone(),
        long_p,
        grid(),
    );
    let schedule = finished_schedule(&long_dir);
    assert_whole(&schedule, "octanol/long/");
    assert!(
        schedule.get(..3) == Some(&CHARGE_SEGMENT[..]),
        "long/'s schedule {schedule:?} does not begin with the charge segment"
    );
    let long_records: Vec<Recorded> = schedule
        .iter()
        .map(|&g| read_back(&long, &long_dir, g))
        .collect();
    let water_dir = root.join(water.name);
    let wleg = Leg::new(
        water.name,
        &water.decoupling,
        &water.template,
        read_back_start(&water, &water_dir),
        seeded(LONG, water.name),
        grid(),
    );
    let water_long = water_dir.join("long");
    let water_records: Vec<Recorded> = finished_schedule(&water_long)
        .iter()
        .map(|&g| read_back(&wleg, &water_long, g))
        .collect();
    assert_whole(&finished_schedule(&water_long), "water/long/");

    std::fs::create_dir_all(&dir).expect("the reseed's directory");
    let p = Protocol {
        seed: RESEED,
        ..long_p
    };
    let leg = Leg::new(sys.name, &sys.decoupling, &sys.template, start, p, grid());
    log(
        &dir,
        "reseed",
        &format!(
            "{}: the charge segment {CHARGE_SEGMENT:?} reseeded, {:#x} for long/'s {:#x}; \
             protocol {p:?}, model {:016x}",
            sys.name, p.seed, long_p.seed, leg.model
        ),
    );
    let records: Vec<Recorded> = CHARGE_SEGMENT
        .iter()
        .map(|&g| window(&leg, &dir, g))
        .collect();

    let (whole, charges, vdw) = (
        chain(&long_records, kt),
        chain(&long_records[..3], kt),
        chain(&long_records[2..], kt),
    );
    assert!(
        (whole.value - (charges.value + vdw.value)).abs()
            <= 64.0 * EPS * (charges.value.abs() + vdw.value.abs()),
        "the segments add to the leg"
    );
    let (a, ea) = segment(&long_records[..3], kt);
    let (b, eb) = segment(&records, kt);
    let (v, ev) = segment(&long_records[2..], kt);
    let d = b - a;
    let combined = (ea * ea + eb * eb).sqrt();
    let mean = 0.5 * (a + b);
    let e_mean = (0.5 * d.abs()).max(0.5 * combined);
    let g_oct = -(mean + v);
    let e_oct = (e_mean * e_mean + ev * ev).sqrt();
    let (w, ew) = segment(&water_records, kt);
    let g_hyd = -w;
    let log_p = -(g_oct - g_hyd) / rt_ln10();
    let sigma = (e_oct * e_oct + ew * ew).sqrt() / rt_ln10();
    let off = (log_p - LOG_P_XYLENE).abs();
    let verdict = if 2.0 * sigma >= LOG_P_XYLENE - LOG_P {
        "uninformative: 2σ is experiment's whole difference from benzene"
    } else if off <= 2.0 * sigma {
        "the method carries to p-xylene"
    } else if off > 3.0 * sigma {
        "the method does not carry to p-xylene"
    } else {
        "not settled"
    };
    let delta = log_p - BENZENE_LOG_P.0;
    let e_delta = (sigma * sigma + BENZENE_LOG_P.1 * BENZENE_LOG_P.1).sqrt();
    let mut lines = Vec::new();
    for (name, rs, seed) in [
        ("long/", &long_records[..3], long_p.seed),
        ("reseed/", &records[..], p.seed),
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
        lines.push(format!(
            "{name} (seed {seed:#x}): the charge segment's ΔG {g:+.3} ± {e:.3} kcal/mol, halves \
             {:+.3}, {:+.3}; {}",
            h[0],
            h[1],
            taus.join(", ")
        ));
    }
    lines.push(format!(
        "the seeds differ by {d:+.3} kcal/mol, {:.2} of their combined quoted error {combined:.3}",
        d.abs() / combined
    ));
    lines.push(format!(
        "with both seeds: the charge segment {mean:+.3} ± {e_mean:.3}, the van der Waals segment \
         {v:+.3} ± {ev:.3}; ΔG_solv,oct {g_oct:+.3} ± {e_oct:.3}, ΔG_hyd {g_hyd:+.3} ± {ew:.3}; \
         log P {log_p:.2} ± {sigma:.2} against {LOG_P_XYLENE}, {:.2} σ: {verdict}. Not in σ: the \
         van der Waals segment's seed spread, unmeasured; a two-seed difference's one degree of \
         freedom; experiment's 3.15 against Sangster's 3.42",
        off / sigma
    ));
    lines.push(format!(
        "Δ log P from benzene's {:.2} ± {:.2}: {delta:+.2} ± {e_delta:.2} against experiment's \
         {:+.2}, reported and not gated; {:.0} s",
        BENZENE_LOG_P.0,
        BENZENE_LOG_P.1,
        LOG_P_XYLENE - LOG_P,
        t.elapsed().as_secs_f64()
    ));
    for line in lines {
        log(&dir, "reseed", &line);
    }
}

/// **The melts and the pilot** (release, one core, in the background): benzene's
/// `the_pilot_measured` on p-xylene. Each leg is melted into `NAME/start.txt`, the production's own
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
                "\n{}: p-xylene and {} solvent molecules, {} atoms; τ in samples of 0.1 ps\n| \
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
