//! **Benzene decoupled from 1-octanol and from water with one model, for its log P**: step 4c-2 of
//! the explicit-solvent track. Benzene is GAFF/AM1-BCC on UFF's bonded terms — `ca` and `ha` by
//! Lorentz–Berthelot, GAFF's 1-4 scales, FreeSolv's charges — the model variant 4's staged run
//! swapped in on W3's water leg, where it hydrated at −0.82 ± 0.15 kcal/mol. The octanol is 4c-1's
//! with GAFF's non-bonded terms (`liquid_octanol_with_gaff.rs`). Dry octanol only.
//!
//! **Kept out of the published crate** (`Cargo.toml`'s `exclude`), with `tests/gaff/`, which holds
//! GAFF's values: they are read from `gaff.dat`, which is GPL v3.
//!
//! **The legs' machinery is in `tests/partition/`**: the model, the boxes, the melts, the digests,
//! the windows, the lock and the report, moved there unchanged for step C-2's decane leg
//! (`benzene_partitioned_into_decane.rs`), which reads this file's water leg back through them. The
//! one addition, an observation written beside each sample, is absent from every window here.
//!
//! # The thermodynamics
//!
//! `log P = −(ΔG_solv,oct − ΔG_hyd) / (RT ln 10)`, each `ΔG_solv = −ΔG_decouple`, benzene decoupled
//! W3's way in each box — charges off, then the van der Waals through the soft core — so that the
//! decoupled state is the solvent alone in its box and benzene alone in vacuum **with its own
//! terms**. Those are the same Hamiltonian at both ends of both legs: UFF's bonded terms, and its
//! own non-bonded pairs GAFF's, Lorentz–Berthelot's, 1-4 scaled and AM1-BCC's, because the field's
//! mixing and 1-4 scales are field-wide and octanol needs them. So no vacuum leg is needed and
//! benzene's own terms cancel from both, except through how the solvent changes benzene's
//! conformations, which for a ring with its C–H held is small. **The staged route kept benzene's own
//! pairs UFF's and QEq's**, unscaled: so the direct water leg is the same cross model as the staged
//! one and differs from it only in that; [`the_partition_measured`] reports both water numbers.
//! Both decouplings are at fixed volume and the same concentration at each end, so the standard
//! state cancels from `log P`, a ratio of concentrations.
//!
//! Experiment: `log K_ow = 2.13` (C. Hansch, A. Leo and D. Hoekman, *Exploring QSAR: Hydrophobic,
//! Electronic and Steric Constants*, ACS 1995, p. 18, as PubChem CID 241's HSDB record quotes it,
//! read 2026-10-21; its HMDB and ICSC records give 2.13 too, and the record's Sangster entry 2.31).
//! So `ΔG_oct − ΔG_hyd = −2.13 RT ln 10 = −2.906 kcal/mol` at 298.15 K. That is **wet** octanol,
//! water-saturated, about 0.27 water by mole fraction (quoted, not read): this file builds dry,
//! which leaves out waters clustered on the hydroxyls around a solute with no hydrogen bond to
//! make. `PeriodicForceField::solvated` takes TIP3P after a many-molecule solute already, so wet
//! octanol needs no library change. Not built.
//!
//! # The boxes
//!
//! - **Octanol**: 4c-1's lattice, 192 sites at 3 × 8 × 8, benzene at its centre and every octanol
//!   with an atom within 2 Å of it removed; **the lattice built denser so that the octanol is at
//!   0.8262 g/cm³ in the volume benzene leaves it**, `N_oct m / (V − v_B)`, `v_B` liquid benzene's
//!   89 cm³/mol (quoted), solved to a fixed point in the number removed ([`octanol_box`]); `r_c`
//!   9 Å, PME at δ = 10⁻⁵ and the neighbour list, every X–H held. Uncompensated, the 188 octanols
//!   kept would be 1.8% light: at a compressibility of order 10⁻⁴ bar⁻¹ (quoted), some 200 bar
//!   more tension, and `P v_B` about 0.4 kcal/mol.
//! - **Water**: W3's box, 512 lattice waters at TIP3P's own density less benzene's five, `r_c` 9
//!   Å, classical Ewald at δ = 10⁻⁵, with benzene's model swapped in ([`gaff::benzene_in_water`]).
//!
//! Each is melted with benzene grown in through the soft core, as W3 grows it ([`OCTANOL_MELT`],
//! [`WATER_MELT`]), and the melt kept in `start.txt` under a header with digests of the box as
//! built and of the model.
//!
//! # The windows
//!
//! W3's: 29 candidates, the charges off in eight steps and the van der Waals in twenty; thirteen
//! windows first, the charges off in two intervals and the van der Waals by 0.1; a window inserted
//! at the candidate between any pair overlapping below 0.1, twice. Each window 20 ps discarded and
//! 200 ps sampled every 100 fs ([`LONG`]), or 4 + 16 ps every 50 fs ([`SHORT`]), as variant 4.
//! Written W4's way under `PANTOMETRY_OCTANOL_DIR`: every window `window_NN.partial` until whole,
//! then renamed; its header the candidate, λ, the protocol, the atom count, the number of
//! candidates, and digests of the start and **of the model** — every atom's `[x, D]` and charge,
//! the mixing, the 1-4 scales and the Ewald parameters — so a window run under another model or
//! from another start is refused by name; and a lock.
//!
//! **The charge segment reseeded** ([`the_octanol_charges_reseeded_measured`], written, not run):
//! the octanol leg's candidates 0, 4 and 8 again, from the same start under the same model and
//! protocol, with only the campaign seed changed ([`RESEED`]), into `octanol/reseed/`. Its windows
//! with charges on had τ of 53 and 35 samples in `long/`, and the leg's halves differ by 0.67
//! against a quoted 0.235: the second seed measures what one seed's error leaves out.
//!
//! # What the melts and the pilot measured
//!
//! [`the_pilot_measured`], run once on 2026-10-21, release, one core, exit 0 after 12 748 s; nothing
//! asserted. Its `start.txt` files are the production's starts.
//!
//! - **The octanol needed its 60 ps at 298.15 K**: `U` per octanol in 2.5 ps blocks fell from 2.31
//!   kcal/mol at the settle's start to 1.48 at 20 ps and was flat from about 40 ps, 1.18–1.39 over
//!   the last 20 (4c-1's GAFF liquid: 1.36); `P` −110 to +60 bar over the same 20 ps. The settle's
//!   halves, 1.88 ± 0.22 and 1.38 ± 0.04, show the drift the first half still had.
//! - **τ**, 300 samples of 0.1 ps after 2 ps, in samples: with the charges on, `∂U/∂λ_e` 7.2 in
//!   water and **22.2 in octanol**; at `λ_v` 1 uncharged 0.6 and 1.1; at `λ_v` 0.3, 10.2 and 3.6. So
//!   200 ps a window gives about 45 independent samples where octanol is slowest, and 16 ps about
//!   four: [`SHORT`] is a test of the machinery there, not a measurement.
//!
//! Default, unoptimised, on boxes of sixteen octanol sites and 27 waters: the model is what it
//! says; the end states are the systems they say; the λ derivatives are the whole energy's; every
//! state of `couplings` is its own `coupling`; the decoupled reciprocal energy is the rest's plus
//! `λ_e` times the cross term; and a killed leg resumes to the same bits. `benzene_bound_in_tip3p.rs`
//! holds benzene's model as built here to variant 4's to the bit.

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
use pantometry_forcefield::free_energy::{refine_schedule, window_seed};
use pantometry_forcefield::{
    Alchemical, AtLambda, Bath, Component, Element, Estimate, Lambda, Mixing, MolecularDynamics,
    PeriodicDecoupling, PeriodicForceField, PmeParameters, Potential, Protocol,
};
use pantometry_units::BOLTZMANN;
use partition::*;

// ---------------------------------------------------------------------------------------------
// Default checks, on small boxes

/// Benzene in a [`SMALL`] octanol box, GAFF's everywhere, through classical Ewald or the mesh; and
/// the systems its end states should be, each built here from the rules written out: benzene alone
/// in the box and in vacuum ([`Vacuum`]'s direct sums), the octanols alone in the box.
struct Small {
    field: PeriodicForceField,
    decoupling: PeriodicDecoupling,
    benzene_in_the_box: PeriodicForceField,
    benzene_in_vacuum: Vacuum,
    solvent_alone: PeriodicForceField,
    at: Vec<[f64; 3]>,
    solvent: usize,
    liquid: Option<Liquid>,
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

fn small_octanol(seed: u64, mesh: bool) -> Small {
    let (b, o) = (benzene(), octanol());
    let (l, mut field, _) = octanol_box(&b, &o, SMALL, SMALL_CUTOFF, 1e-5, seed, false);
    let cell = l.cell();
    let p = l.field.ewald().parameters();
    let bm = Model::benzene(&b.component);
    let om = Model::octanol(&o.component);
    let mut benzene_in_the_box = rules(
        PeriodicForceField::new(&b.component, &b.types, cell, p).expect("benzene"),
        bm.vdw.clone(),
        bm.charges.clone(),
    );
    let parts: Vec<&Component> = (0..l.octanols).map(|_| &o.component).collect();
    let joined = joined(&parts);
    let (mut vdw, mut q) = (Vec::new(), Vec::new());
    for _ in 0..l.octanols {
        vdw.extend_from_slice(&om.vdw);
        q.extend_from_slice(&om.charges);
    }
    let mut solvent_alone = rules(
        PeriodicForceField::new(&joined, &l.types[12..], cell, p).expect("octanol alone"),
        vdw,
        q,
    );
    if mesh {
        let m = PmeParameters::for_accuracy(&cell, &p, 1e-6, l.at.len());
        field = field.with_mesh(m);
        benzene_in_the_box = benzene_in_the_box.with_mesh(m);
        solvent_alone = solvent_alone.with_mesh(m);
    }
    let group: Vec<bool> = (0..l.at.len()).map(|i| i < 12).collect();
    Small {
        decoupling: PeriodicDecoupling::new(&field, &group).expect("benzene decouples"),
        field,
        benzene_in_the_box,
        benzene_in_vacuum: Vacuum::new(
            &b.component,
            &b.types,
            bm.vdw,
            bm.charges,
            Mixing::LorentzBerthelot,
            [0.5, 1.0 / 1.2],
        ),
        solvent_alone,
        at: l.at.clone(),
        solvent: l.octanols,
        template: l.dynamics(),
        liquid: Some(l),
    }
}

/// Benzene in 27 lattice waters less those within 2.3 Å of it, 9.3 Å across, `r_c` = 4.6 Å, its
/// model swapped in; and its end states' systems, the waters alone with the rules written out.
fn small_water(seed: u64) -> Small {
    let b = benzene();
    let s = with_gaff_benzene(
        Solvated::new(&b, 3, seed, 2.3 * ANGSTROM, 4.6 * ANGSTROM, 1e-5),
        &b,
    );
    let cell = s.cell();
    let p = s.field.ewald().parameters();
    let bm = Model::benzene(&b.component);
    Small {
        decoupling: s.decoupling(),
        field: s.field.clone(),
        benzene_in_the_box: rules(
            PeriodicForceField::new(&b.component, &b.types, cell, p).expect("benzene"),
            bm.vdw.clone(),
            bm.charges.clone(),
        ),
        benzene_in_vacuum: Vacuum::new(
            &b.component,
            &b.types,
            bm.vdw,
            bm.charges,
            Mixing::LorentzBerthelot,
            [0.5, 1.0 / 1.2],
        ),
        solvent_alone: PeriodicForceField::tip3p(cell, p, s.waters)
            .with_mixing(Mixing::LorentzBerthelot)
            .with_one_four(0.5, 1.0 / 1.2),
        at: s.at.clone(),
        solvent: s.waters,
        liquid: None,
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

/// **The model is what it says**, in both small boxes: benzene's every carbon `ca` and hydrogen
/// `ha` ([`gaff`], by element, not by FreeSolv's types), every octanol atom its GAFF type, the
/// waters TIP3P's; benzene's charges FreeSolv's — every one ±0.1301 e but C1's −0.1302, summing to
/// −0.0001 e, and **not QEq's** — octanol's summing to zero with the oxygen at −0.6023 and its
/// hydrogen +0.3984; Lorentz–Berthelot and GAFF's 1-4 scales, the 1-4 pairs benzene's 21 and
/// octanol's 66 each; and **benzene the same model in both boxes**: its `[x, D]`, charges and own
/// pairs' scales the same bits.
#[test]
fn the_model_is_what_it_says() {
    let oct = small_octanol(0xE5D, false);
    let wat = small_water(0xE5E);
    let b = benzene();
    let qeq = Benzene::new(None).charges;
    for (name, s) in [("octanol", &oct), ("water", &wat)] {
        let f = &s.field;
        let vdw = f.van_der_waals();
        let q = f.charges();
        for i in 0..12 {
            let kind = match b.component.atoms()[i].element {
                Element::C => "ca",
                Element::H => "ha",
                other => panic!("benzene has no {other:?}"),
            };
            assert_eq!(
                vdw[i].map(f64::to_bits),
                gaff(kind).map(f64::to_bits),
                "{name} {i}"
            );
            assert!(
                q[i] == -0.1302 || q[i] == -0.1301 || q[i] == 0.1301,
                "{name}: benzene's atom {i} carries {} e",
                q[i]
            );
            assert_ne!(q[i].to_bits(), qeq[i].to_bits(), "{name}: QEq's charge");
        }
        let net: f64 = q[..12].iter().sum();
        assert!((net + 0.0001).abs() < 1e-12, "{name}: benzene's {net} e");
        assert_eq!(f.mixing(), Mixing::LorentzBerthelot, "{name}");
        assert_eq!(f.one_four(), ONE_FOUR, "{name}");
        let per = if name == "octanol" { 66 } else { 0 };
        assert_eq!(f.one_four_pairs().len(), 21 + per * s.solvent, "{name}");
        println!(
            "{name}: benzene and {} solvent molecules, {} atoms, {} 1-4 pairs scaled",
            s.solvent,
            q.len(),
            f.one_four_pairs().len()
        );
    }
    // The octanols, molecule by molecule, and the waters.
    let o = octanol();
    let kinds = gaff::gaff_types(&o.component);
    let l = oct.liquid.as_ref().expect("a liquid");
    for m in &l.molecules[1..] {
        let q: Vec<f64> = oct.field.charges()[m.clone()].to_vec();
        for (k, i) in m.clone().enumerate() {
            assert_eq!(
                oct.field.van_der_waals()[i].map(f64::to_bits),
                gaff(&kinds[k]).map(f64::to_bits)
            );
        }
        assert!(q.iter().sum::<f64>().abs() < 1e-12);
        let oxygen = index_of(&o.component, Element::O);
        assert_eq!(q[oxygen], -0.6023);
        let ho = o
            .component
            .bonds()
            .iter()
            .filter(|bd| bd.atoms.contains(&oxygen))
            .map(|bd| bd.atoms[0] + bd.atoms[1] - oxygen)
            .find(|&h| o.component.atoms()[h].element == Element::H)
            .expect("the hydroxyl's hydrogen");
        assert_eq!(q[ho], 0.3984);
    }
    let [xo, dox] = pantometry_forcefield::water::oxygen_vdw();
    for (k, v) in wat.field.van_der_waals()[12..].iter().enumerate() {
        let want = if k % 3 == 0 { [xo, dox] } else { [0.0, 0.0] };
        assert_eq!(v.map(f64::to_bits), want.map(f64::to_bits));
    }
    // One benzene model in both boxes.
    let bits =
        |v: &[[f64; 2]]| -> Vec<[u64; 2]> { v.iter().map(|x| x.map(f64::to_bits)).collect() };
    assert_eq!(
        bits(&oct.field.van_der_waals()[..12]),
        bits(&wat.field.van_der_waals()[..12])
    );
    assert_eq!(
        oct.field.charges()[..12]
            .iter()
            .map(|x| x.to_bits())
            .collect::<Vec<_>>(),
        wat.field.charges()[..12]
            .iter()
            .map(|x| x.to_bits())
            .collect::<Vec<_>>()
    );
    assert_eq!(oct.decoupling.own_scales(), wat.decoupling.own_scales());
    assert_eq!(
        oct.decoupling
            .own_scales()
            .iter()
            .filter(|s| **s == ONE_FOUR)
            .count(),
        21
    );
}

fn index_of(c: &Component, e: Element) -> usize {
    c.atoms()
        .iter()
        .position(|a| a.element == e)
        .expect("the element")
}

/// **The octanol is at its density in the volume benzene leaves it**, in the production's box: 192
/// sites, benzene at the centre, every octanol within 2 Å of it removed, the lattice built denser
/// to make up — `N_oct m / (V − v_B)` is 0.8262 g/cm³ within its rounding — and the fixed point
/// reached, the number removed the same at the density it gives.
#[test]
fn the_octanol_is_at_its_density_around_benzene() {
    let (b, o) = (benzene(), octanol());
    let (l, _, lattice) = octanol_box(&b, &o, octanol::LIQUID, 9.0 * ANGSTROM, 1e-5, 0x0C71, true);
    let v = l.cell().volume();
    let rho = l.octanols as f64 * o.mass() / (v - BENZENE_VOLUME);
    println!(
        "{} of 192 octanols kept; the lattice at {:.6} g/cm³; the octanol at {:.10} g/cm³ in \
         the volume benzene leaves; box {:.3} × {:.3} × {:.3} Å",
        l.octanols,
        lattice / 1e3,
        rho / 1e3,
        l.cell().lengths()[0] / ANGSTROM,
        l.cell().lengths()[1] / ANGSTROM,
        l.cell().lengths()[2] / ANGSTROM
    );
    assert!(l.octanols < 192, "benzene took somebody's place");
    assert!((rho - DENSITY).abs() <= 1e3 * EPS * DENSITY, "{rho}");
    // Plain: what it would have been with no compensation.
    let (plain, _, _) = octanol_box(&b, &o, octanol::LIQUID, 9.0 * ANGSTROM, 1e-5, 0x0C71, false);
    let plain_rho = plain.octanols as f64 * o.mass() / (plain.cell().volume() - BENZENE_VOLUME);
    println!(
        "uncompensated: {} kept, {:.5} g/cm³, {:+.2}%",
        plain.octanols,
        plain_rho / 1e3,
        100.0 * (plain_rho / DENSITY - 1.0)
    );
    assert!(
        (plain_rho - DENSITY).abs() > 1e-3 * DENSITY,
        "the compensation does something"
    );
}

/// **The end states are the systems they say**, in both small boxes, at a configuration
/// constrained dynamics left. Fully coupled: the whole field less benzene alone in the box, plus
/// benzene alone in vacuum with its own terms — UFF's bonded, its 36 own pairs Lorentz–Berthelot's
/// and AM1-BCC's with its 21 1-4 pairs scaled, by [`Vacuum`]'s direct sums; fully decoupled, the
/// solvent alone in its box plus benzene in vacuum, energy and every force. Each within `10⁻¹²` of
/// the parts' magnitudes: the same terms summed in other orders.
#[test]
fn the_end_states_are_the_systems_they_say() {
    for (name, s) in [
        ("octanol", small_octanol(0xE5D, false)),
        ("water", small_water(0xE5E)),
    ] {
        let at = moved(&s, Lambda::new(0.0, 0.5, 0.7), 20);
        let n = at.len();
        assert_eq!(s.decoupling.own_pairs().len(), 36);
        let mut f = vec![[0.0; 3]; n];
        let coupled = s.decoupling.energy_and_forces(&at, Lambda::COUPLED, &mut f);
        let whole = s.field.energy(&at);
        let in_box = s.benzene_in_the_box.energy(&at[..12]);
        let mut vf = vec![[0.0; 3]; 12];
        let vac = s.benzene_in_vacuum.energy_and_forces(&at[..12], &mut vf);
        let expected = whole.total - in_box.total + vac;
        let scale = magnitude(&whole);
        println!(
            "{name}: coupled {coupled:.10e} J against {expected:.10e}; benzene's own images {:+.4e} \
             kcal/mol",
            kcal(in_box.total - vac)
        );
        assert!(
            (coupled - expected).abs() <= 1e-12 * scale,
            "{name} coupled"
        );
        let decoupled = s
            .decoupling
            .energy_and_forces(&at, Lambda::new(0.0, 0.0, 0.0), &mut f);
        let solvent = s.solvent_alone.evaluate(&at[12..]);
        let expected = solvent.energy.total + vac;
        println!("{name}: decoupled {decoupled:.10e} J against {expected:.10e}");
        assert!(
            (decoupled - expected).abs() <= 1e-12 * scale,
            "{name} decoupled"
        );
        let fscale = largest_force(&f);
        for i in 0..n {
            let want = if i < 12 {
                vf[i]
            } else {
                solvent.forces[i - 12]
            };
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

/// **`∂U/∂λ_e` and `∂U/∂λ_v` are the whole energy's derivatives**, in both small boxes, at four
/// states — the charges part-way, both part-way, the van der Waals part-way with no charge, and
/// near the decoupled end where the soft core is softest — each at a configuration dynamics at
/// that state left, against central differences of `energy_and_forces`, each to the tolerance it
/// earns and that tolerance a part in 10³ of the derivative.
#[test]
fn the_lambda_derivatives_are_the_whole_energys() {
    for (name, s) in [
        ("octanol", small_octanol(0xDE1, false)),
        ("water", small_water(0xDE2)),
    ] {
        for (lambda, steps) in [
            (Lambda::new(0.0, 0.5, 1.0), 20),
            (Lambda::new(0.0, 0.4, 0.8), 20),
            (Lambda::new(0.0, 0.0, 0.5), 20),
            (Lambda::new(0.0, 0.0, 0.1), 40),
        ] {
            let at = moved(&s, lambda, steps);
            let c = s.decoupling.coupling(&at, lambda);
            let scale = magnitude(&s.field.energy(&at));
            for (k, h) in [(1usize, 1e-3), (2, 1e-4)] {
                let at_k = lambda.components()[k];
                if at_k < 2.0 * h || at_k > 1.0 - 2.0 * h {
                    continue;
                }
                let energy = |e: f64| {
                    let mut l = lambda;
                    if k == 1 {
                        l.electrostatics += e;
                    } else {
                        l.van_der_waals += e;
                    }
                    let mut f = vec![[0.0; 3]; at.len()];
                    (s.decoupling.energy_and_forces(&at, l, &mut f), scale)
                };
                let (dl, tolerance) = finite_difference(energy, h);
                println!(
                    "{name} {lambda:?}: ∂U/∂λ[{k}] = {:+.10e} J, by differences {dl:+.10e} ± \
                     {tolerance:.1e}",
                    c.gradient[k]
                );
                assert!(
                    (c.gradient[k] - dl).abs() <= tolerance,
                    "{name} {lambda:?} {k}"
                );
                assert!(
                    tolerance < 1e-3 * c.gradient[k].abs(),
                    "the check can see a part in 10³"
                );
            }
        }
    }
}

/// **Every state of `couplings` is its own `coupling`, to the bit**, at all 29 candidates and in a
/// shuffled call; and `U(λ) − U(0)` is the coupling's energy within twice the evaluations'
/// rounding, `128 ε` of the parts, in both small boxes.
#[test]
fn every_state_of_couplings_is_its_own_coupling() {
    for (name, s) in [
        ("octanol", small_octanol(0xC0B, false)),
        ("water", small_water(0xC0C)),
    ] {
        let at = moved(&s, Lambda::new(0.0, 0.0, 0.6), 20);
        let states = grid();
        let all = s.decoupling.couplings(&at, &states);
        let some = s
            .decoupling
            .couplings(&at, &[states[20], states[3], states[20]]);
        for (k, l) in states.iter().enumerate() {
            let one = s.decoupling.coupling(&at, *l);
            assert_eq!(
                one.energy.to_bits(),
                all[k].energy.to_bits(),
                "{name} {l:?}"
            );
            assert_eq!(
                one.gradient.map(f64::to_bits),
                all[k].gradient.map(f64::to_bits)
            );
        }
        assert_eq!(some[0].energy.to_bits(), all[20].energy.to_bits());
        assert_eq!(some[1].energy.to_bits(), all[3].energy.to_bits());
        assert_eq!(
            some[2].gradient.map(f64::to_bits),
            all[20].gradient.map(f64::to_bits)
        );
        let scale = magnitude(&s.field.energy(&at));
        let mut f = vec![[0.0; 3]; at.len()];
        let u0 = s
            .decoupling
            .energy_and_forces(&at, Lambda::new(0.0, 0.0, 0.0), &mut f);
        let mut worst = 0.0f64;
        for (k, l) in states.iter().enumerate() {
            let du = s.decoupling.energy_and_forces(&at, *l, &mut f) - u0;
            let r = (du - all[k].energy).abs() / (128.0 * EPS * scale);
            worst = worst.max(r);
            assert!(r <= 1.0, "{name} {l:?}: {du:e} against {:e}", all[k].energy);
        }
        println!("{name}: U(λ) − U(0) against the coupling, the worst {worst:.3} of 128 ε");
    }
}

/// **The decoupled reciprocal energy is `rest + λ_e · cross`**, with GAFF's model and both rules
/// on: the field's reciprocal energy with only the rest's charges, with only benzene's, and with
/// all, against the decoupling's at three states, within `64 ε` of the parts — the octanol box
/// through classical Ewald and through the mesh the production uses, the water box classical. Not
/// a difference, so it tells the rest from the group; and the cross term is required to be seen.
#[test]
fn the_decoupled_reciprocal_energy_is_the_rest_and_the_cross_terms() {
    for (name, s) in [
        ("octanol, classical", small_octanol(0x2EC, false)),
        ("octanol, mesh", small_octanol(0x2EC, true)),
        ("water, classical", small_water(0x2ED)),
    ] {
        let at = &s.at;
        let q = s.field.charges().to_vec();
        let reciprocal = |keep: &dyn Fn(usize) -> bool| {
            let masked: Vec<f64> = (0..q.len())
                .map(|i| if keep(i) { q[i] } else { 0.0 })
                .collect();
            s.field
                .clone()
                .with_charges(masked)
                .energy(at)
                .ewald
                .reciprocal
        };
        let all = reciprocal(&|_| true);
        let rest = reciprocal(&|i| i >= 12);
        let own = reciprocal(&|i| i < 12);
        let scale = magnitude(&s.field.energy(at));
        for l in [
            Lambda::COUPLED,
            Lambda::new(0.0, 0.5, 1.0),
            Lambda::new(0.0, 0.0, 1.0),
        ] {
            let e = s.decoupling.energy(at, l).ewald.reciprocal;
            let want = rest + l.electrostatics * (all - rest - own);
            println!(
                "{name}, {l:?}: {e:+.12e} J against {want:+.12e}, {:.3} of 64 ε",
                (e - want).abs() / (64.0 * EPS * scale)
            );
            assert!((e - want).abs() <= 64.0 * EPS * scale, "{name} {l:?}");
        }
        assert!(
            (all - rest - own).abs() > 1e-6 * scale,
            "{name}: the cross term"
        );
    }
}

/// A fresh directory under the system's temporary one, named by the process and `tag`.
#[cfg(not(target_family = "wasm"))]
fn fresh_dir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("pantometry-4c2-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("a temporary directory");
    d
}

/// The resume test's melt: a few steps growing benzene in and a few coupled.
#[cfg(not(target_family = "wasm"))]
const TEST_MELT: [MeltStage; 2] = [
    stage("grow", GROW_HALF, 0.0015, 0.5 * FS, 50e12, KELVIN),
    stage("coupled", Lambda::COUPLED, 0.001, 0.5 * FS, 50e12, KELVIN),
];

/// **A killed leg resumes to the same bits**, through this file's own [`melted`], [`Leg`] and
/// [`window`] on the small octanol box, a melt of five steps and windows of eight 0.5 fs steps and
/// three samples at candidates 0 and 8, in temporary directories:
///
/// - **uninterrupted**: the melt and both windows;
/// - **killed and resumed**: the melt and window 0, then what a run killed during window 8 leaves —
///   its `.partial`, a header and one sample — then the leg again. Window 0 and `start.txt` are
///   read, their bytes unchanged; the partial is begun afresh and renamed; and **every gradient and
///   energy of both records is the uninterrupted run's to the bit**.
///
/// And each guard refuses by name: the melt asked for with another seed; a window asked for from a
/// start one ulp away; **a window asked for from the same start under another model** — benzene's
/// charges QEq's — which only the header's model digest can tell; and a `start.txt` cut short.
#[cfg(not(target_family = "wasm"))]
#[test]
fn a_killed_leg_resumes_to_the_same_bits() {
    let s = small_octanol(0x2E5, false);
    let l = s.liquid.as_ref().expect("a liquid");
    let sys = octanol_system(l, s.field.clone(), &TEST_MELT);
    let protocol = Protocol {
        time_step: 0.5e-15,
        temperature: KELVIN,
        friction: 50e12,
        equilibration: 2,
        stride: 2,
        samples: 3,
        seed: 0x2E6,
    };
    let run = |dir: &std::path::Path, schedule: &[usize]| -> Vec<Recorded> {
        let start = melted(&sys, dir, MELT_SEED);
        let leg = Leg::new(
            "octanol",
            &sys.decoupling,
            &sys.template,
            start,
            protocol,
            grid(),
        );
        schedule.iter().map(|&g| window(&leg, dir, g)).collect()
    };
    let whole = fresh_dir("whole");
    let full = run(&whole, &[0, 8]);
    let killed = fresh_dir("killed");
    run(&killed, &[0]);
    let bytes = |p: std::path::PathBuf| std::fs::read(p).expect("a file");
    let (finished, start_bytes) = (
        bytes(killed.join("window_00.txt")),
        bytes(killed.join("start.txt")),
    );
    let text = std::fs::read_to_string(whole.join("window_08.txt")).unwrap();
    let partial: String = text.lines().take(2).map(|l| format!("{l}\n")).collect();
    std::fs::write(killed.join("window_08.partial"), partial).unwrap();
    let resumed = run(&killed, &[0, 8]);
    assert!(
        bytes(killed.join("window_00.txt")) == finished,
        "window 0 read, not run"
    );
    assert!(bytes(killed.join("start.txt")) == start_bytes);
    assert!(!killed.join("window_08.partial").exists());
    let f64s = |v: &[f64]| -> Vec<u64> { v.iter().map(|x| x.to_bits()).collect() };
    assert_eq!(full.len(), resumed.len());
    for (a, b) in full.iter().zip(&resumed) {
        assert_eq!(a.grid, b.grid);
        assert_eq!(a.gradients.len(), protocol.samples);
        let flat = |r: &Recorded| -> Vec<u64> {
            let mut v = f64s(&r.gradients.iter().flatten().copied().collect::<Vec<_>>());
            for e in &r.energies {
                v.extend(f64s(e));
            }
            v
        };
        assert_eq!(flat(a), flat(b), "candidate {}", a.grid);
    }
    println!(
        "resumed: {} records, every number the uninterrupted run's",
        resumed.len()
    );

    let why = refusal(|| drop(melted(&sys, &killed, MELT_SEED + 1)));
    assert!(why.contains("melted differently"), "{why}");
    let start = melted(&sys, &killed, MELT_SEED);
    let mut moved = start.clone();
    moved[0][0] = f64::from_bits(moved[0][0].to_bits() + 1);
    let leg = Leg::new(
        "octanol",
        &sys.decoupling,
        &sys.template,
        moved,
        protocol,
        grid(),
    );
    let why = refusal(|| drop(window(&leg, &killed, 0)));
    assert!(why.contains("was run differently"), "{why}");
    // The same start under another model: benzene's charges QEq's.
    let mut q = sys.field.charges().to_vec();
    q[..12].copy_from_slice(&Benzene::new(None).charges);
    let other = PeriodicDecoupling::new(&sys.field.clone().with_charges(q), sys.decoupling.group())
        .expect("benzene decouples");
    let leg = Leg::new("octanol", &other, &sys.template, start, protocol, grid());
    let why = refusal(|| drop(window(&leg, &killed, 0)));
    assert!(why.contains("was run differently"), "another model: {why}");
    let cut = fresh_dir("cut");
    std::fs::write(
        cut.join("start.txt"),
        &start_bytes[..start_bytes.len() - 30],
    )
    .unwrap();
    let why = refusal(|| drop(melted(&sys, &cut, MELT_SEED)));
    assert!(why.contains("not whole"), "{why}");
    for d in [whole, killed, cut] {
        let _ = std::fs::remove_dir_all(d);
    }
}

/// **The reseeded windows' headers differ from `long/`'s only in the seed**, at each of the charge
/// segment's candidates: the same leg and start, [`LONG`] seeded as `long/` was and [`reseeded`],
/// one word of the header apart, the seed's, `long/`'s `0x4c20c71` against [`RESEED`]'s; every
/// other field of the protocol the same. And **every window seed the reseed gives is apart from
/// every one the two `long/` legs give** at all 29 candidates, from the melts' and from the three
/// campaign seeds themselves.
#[test]
fn the_reseed_differs_from_long_only_in_its_seed() {
    let s = small_octanol(0x2E5, false);
    let long_p = seeded(LONG, "octanol");
    let p = reseeded(long_p);
    assert_eq!(Protocol { seed: 0, ..p }, Protocol { seed: 0, ..long_p });
    assert_ne!(p.seed, long_p.seed);
    let leg = |p: Protocol| {
        Leg::new(
            "octanol",
            &s.decoupling,
            &s.template,
            s.at.clone(),
            p,
            grid(),
        )
    };
    let (long, re) = (leg(long_p), leg(p));
    for g in CHARGE_SEGMENT {
        let (a, b) = (header(&long, g), header(&re, g));
        let (wa, wb): (Vec<&str>, Vec<&str>) = (
            a.split_whitespace().collect(),
            b.split_whitespace().collect(),
        );
        assert_eq!(wa.len(), wb.len());
        let differ: Vec<usize> = (0..wa.len()).filter(|&i| wa[i] != wb[i]).collect();
        assert_eq!(differ.len(), 1, "candidate {g}:\n{a}\n{b}");
        let i = differ[0];
        assert_eq!(
            wa[i - 7],
            "protocol",
            "the seed is the protocol's seventh word"
        );
        assert_eq!(wa[i], "0x4c20c71");
        assert_eq!(wb[i], format!("{RESEED:#x}"));
        println!(
            "candidate {g}: the headers differ only at word {i}, {} against {}",
            wa[i], wb[i]
        );
    }
    let mut others: Vec<u64> = vec![MELT_SEED];
    for name in ["octanol", "water"] {
        let base = seeded(LONG, name).seed;
        others.push(base);
        others.extend((0..grid().len()).map(|g| window_seed(base, g)));
    }
    others.extend((0..OCTANOL_MELT.len()).map(|k| window_seed(MELT_SEED, k)));
    assert!(!others.contains(&RESEED));
    for g in CHARGE_SEGMENT {
        assert!(!others.contains(&window_seed(RESEED, g)), "candidate {g}");
    }
    assert_ne!(RESEED_TAG, "long");
}

/// **The reseed refuses `long/`'s files, and resumes**, through the measurement's own
/// [`read_back_start`], [`read_back`], [`finished_schedule`], [`reseeded`] and [`window`] on the
/// small octanol box, a melt of five steps and windows of eight 0.5 fs steps and three samples at
/// the charge segment's candidates, `long/` seeded as the production's:
///
/// - a start or a window that is not there is **refused by name and nothing is written** — no
///   melt, no window;
/// - `long/`'s window copied into `reseed/` is refused, and the reseed's read back as `long/`'s;
/// - **the seed reaches the dynamics**: every reseeded record differs from `long/`'s;
/// - a reseed killed during its last window — its `.partial` a header and one sample — resumes to
///   the uninterrupted reseed's bits, the finished windows' bytes unchanged.
#[cfg(not(target_family = "wasm"))]
#[test]
fn the_reseed_refuses_long_and_resumes() {
    let s = small_octanol(0x2E5, false);
    let l = s.liquid.as_ref().expect("a liquid");
    let sys = octanol_system(l, s.field.clone(), &TEST_MELT);
    let tiny = Protocol {
        time_step: 0.5e-15,
        temperature: KELVIN,
        friction: 50e12,
        equilibration: 2,
        stride: 2,
        samples: 3,
        seed: seeded(LONG, "octanol").seed,
    };
    let root = fresh_dir("reseed");
    let why = refusal(|| drop(read_back_start(&sys, &root)));
    assert!(why.contains("is not there"), "{why}");
    assert!(is_empty(&root), "no melt was run");
    let start = melted(&sys, &root, MELT_SEED);
    let bits = |at: &[[f64; 3]]| -> Vec<u64> { at.iter().flatten().map(|x| x.to_bits()).collect() };
    assert_eq!(bits(&read_back_start(&sys, &root)), bits(&start));
    let leg = |p: Protocol| {
        Leg::new(
            "octanol",
            &sys.decoupling,
            &sys.template,
            start.clone(),
            p,
            grid(),
        )
    };
    let (long, re) = (leg(tiny), leg(reseeded(tiny)));
    let dirs = ["long", RESEED_TAG, "killed"].map(|d| root.join(d));
    for d in &dirs {
        std::fs::create_dir_all(d).unwrap();
    }
    let [long_dir, dir, killed] = &dirs;
    let long_records: Vec<Recorded> = CHARGE_SEGMENT
        .iter()
        .map(|&g| window(&long, long_dir, g))
        .collect();
    assert_eq!(finished_schedule(long_dir), CHARGE_SEGMENT);

    let why = refusal(|| drop(read_back(&re, dir, 0)));
    assert!(why.contains("is not there"), "{why}");
    assert!(is_empty(dir), "no window was run");
    std::fs::copy(long_dir.join("window_00.txt"), dir.join("window_00.txt")).unwrap();
    let why = refusal(|| drop(window(&re, dir, 0)));
    assert!(why.contains("was run differently"), "long/'s file: {why}");
    std::fs::remove_file(dir.join("window_00.txt")).unwrap();
    let whole: Vec<Recorded> = CHARGE_SEGMENT
        .iter()
        .map(|&g| window(&re, dir, g))
        .collect();
    let why = refusal(|| drop(read_back(&long, dir, 0)));
    assert!(
        why.contains("was run differently"),
        "the reseed's file: {why}"
    );
    for (a, b) in long_records.iter().zip(&whole) {
        assert_eq!(a.grid, b.grid);
        assert_ne!(record_bits(a), record_bits(b), "candidate {}", a.grid);
    }

    for &g in &CHARGE_SEGMENT[..2] {
        window(&re, killed, g);
    }
    let bytes = |p: std::path::PathBuf| std::fs::read(p).expect("a file");
    let finished = [0, 4].map(|g| bytes(killed.join(format!("window_{g:02}.txt"))));
    let text = std::fs::read_to_string(dir.join("window_08.txt")).unwrap();
    let partial: String = text.lines().take(2).map(|l| format!("{l}\n")).collect();
    std::fs::write(killed.join("window_08.partial"), partial).unwrap();
    let resumed: Vec<Recorded> = CHARGE_SEGMENT
        .iter()
        .map(|&g| window(&re, killed, g))
        .collect();
    for (k, g) in [0, 4].into_iter().enumerate() {
        assert!(
            bytes(killed.join(format!("window_{g:02}.txt"))) == finished[k],
            "window {g} read, not run"
        );
    }
    assert!(!killed.join("window_08.partial").exists());
    for (a, b) in whole.iter().zip(&resumed) {
        assert_eq!(record_bits(a), record_bits(b), "candidate {}", a.grid);
    }
    println!(
        "reseeded: {} records, each apart from long/'s and resumed to the bit",
        resumed.len()
    );
    let _ = std::fs::remove_dir_all(&root);
}

// ---------------------------------------------------------------------------------------------
// The measurements: written, costed, not run

/// One leg at protocol `base` (its seed [`seeded`]), in `ROOT/NAME/TAG/`: the melt from
/// `ROOT/NAME/start.txt` ([`melted`]), W3's thirteen windows, the insertion, and the report, under
/// the lock `ROOT/NAME/TAG.lock`.
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
            "{}: benzene and {} solvent molecules, {} atoms, box {:?} Å, protocol {p:?}, model \
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

/// **Both legs at `base` and benzene's log P**: the water first, then the octanol, each through
/// [`run_leg`], every window resumable; then `log P` from the two direct legs, and from the octanol
/// leg with the staged water route's −0.82 ± 0.15, against experiment's 2.13. Nothing about the
/// model is asserted.
fn partition_measured(tag: &str, base: Protocol) {
    let root = results_dir();
    let t = std::time::Instant::now();
    let water = run_leg(&root, &water_leg(), tag, base);
    let oct = run_leg(&root, &octanol_leg(), tag, base);
    let rt_ln10 = kcal(BOLTZMANN.to_si() * KELVIN) * std::f64::consts::LN_10;
    let (g_hyd, e_hyd) = (-water.bar.0, water.bar.1);
    let (g_oct, e_oct) = (-oct.bar.0, oct.bar.1);
    let log_p = |hyd: f64| -(g_oct - hyd) / rt_ln10;
    let err = |e: f64| (e_oct * e_oct + e * e).sqrt() / rt_ln10;
    let lines = [
        format!(
            "ΔG_hyd (direct) {g_hyd:+.3} ± {e_hyd:.3} kcal/mol (TI {:+.3}, halves {:+.3}, {:+.3}, \
             lowest overlap {:.3}); the staged route's {:+.2} ± {:.2}; FreeSolv's GAFF/AM1-BCC \
             −0.81 ± 0.02, experiment −0.90 ± 0.20",
            -water.ti.0,
            -water.halves[0],
            -water.halves[1],
            water.overlap,
            HYDRATION_STAGED.0,
            HYDRATION_STAGED.1
        ),
        format!(
            "ΔG_solv,oct {g_oct:+.3} ± {e_oct:.3} kcal/mol (TI {:+.3}, halves {:+.3}, {:+.3}, \
             lowest overlap {:.3})",
            -oct.ti.0, -oct.halves[0], -oct.halves[1], oct.overlap
        ),
        format!(
            "ΔG_oct − ΔG_hyd {:+.3} (direct), {:+.3} (staged water); experiment {:+.3} = −2.13 RT \
             ln 10",
            g_oct - g_hyd,
            g_oct - HYDRATION_STAGED.0,
            -LOG_P * rt_ln10
        ),
        format!(
            "log P {:.2} ± {:.2} from the two direct legs; {:.2} ± {:.2} with the staged water; \
             experiment {LOG_P} (Hansch, Leo and Hoekman 1995, via PubChem CID 241); dry octanol; \
             {:.0} s",
            log_p(g_hyd),
            err(e_hyd),
            log_p(HYDRATION_STAGED.0),
            err(HYDRATION_STAGED.1),
            t.elapsed().as_secs_f64()
        ),
    ];
    for line in lines {
        log(&root, tag, &line);
    }
}

/// **Benzene's log P at W3's protocol** (release, one core; run once, at `3cb1427`, exit 0 on
/// 2026-10-23 after 70 562 s, the water leg's windows read back from `c51c91e`'s run): log P
/// 2.06 ± 0.20 from the two direct legs and 2.28 ± 0.20 with the staged water, against 2.13.
/// Nothing here asserts it. Both legs, each window 20 ps
/// discarded and 200 ps sampled ([`LONG`]), under `PANTOMETRY_OCTANOL_DIR` in `water/long/` and
/// `octanol/long/`. `PANTOMETRY_OCTANOL_DIR=C:/pantometry-w4/octanol cargo test --release --locked
/// -p pantometry-forcefield --test benzene_partitioned_into_octanol the_partition_measured --
/// --ignored --nocapture`.
#[test]
#[ignore = "both legs at W3's protocol, thirteen windows each: about a day with --release, run once (2026-10-23)"]
fn the_partition_measured() {
    partition_measured("long", LONG);
}

/// **The same at variant 4's short protocol**, 4 + 16 ps a window ([`SHORT`]), in `*/short/`.
#[test]
#[ignore = "both legs at 4 + 16 ps a window: hours with --release, not run"]
fn the_partition_short_measured() {
    partition_measured("short", SHORT);
}

// ---------------------------------------------------------------------------------------------
// The octanol leg's charge segment, reseeded

/// **The charge segment's windows**: candidates 0, 4 and 8, `λ_e` 1, 0.5 and 0 at full van der
/// Waals. The first two are the octanol leg's slowest windows (τ 53 and 35 samples in `long/`), and
/// 8 closes the segment, where the van der Waals segment begins.
const CHARGE_SEGMENT: [usize; 3] = [0, 4, 8];

/// **The reseeded campaign's seed**: window `g` is seeded `window_seed(RESEED, g)`. Not the octanol
/// leg's `0x4C2_0C71` ([`seeded`]), nor the water leg's `0x4C2_3A7E`, nor the melts' `0x4C2_E01`;
/// [`the_reseed_differs_from_long_only_in_its_seed`] holds every window seed it gives apart from
/// every seed those three give.
const RESEED: u64 = 0x4C2_5EED;

/// The reseeded windows' directory, `ROOT/octanol/reseed/`, beside `long/` and never in it.
const RESEED_TAG: &str = "reseed";

/// `p` with only its seed changed, to [`RESEED`].
fn reseeded(p: Protocol) -> Protocol {
    Protocol { seed: RESEED, ..p }
}

/// **The octanol leg's charge segment at a second seed, and log P with it** (release, one core;
/// written, not run). The same leg, model, start (`octanol/start.txt`, header-checked) and
/// protocol as `long/`, [`LONG`] at 20 + 200 ps, with only the seed changed, to [`RESEED`]: the
/// windows at [`CHARGE_SEGMENT`] through [`window`], resumable and headed, into
/// `PANTOMETRY_OCTANOL_DIR/octanol/reseed/` under the lock `octanol/reseed.lock`. Before any of
/// them runs, every file it reads is read and checked: the octanol leg's `long/` windows and the
/// water leg's, through [`read_back`], so a missing one fails in seconds and none is run.
///
/// It prints the charge segment's BAR from both seeds and their difference; a single run's spread
/// that difference implies, `|d|/√2`, and the whole octanol leg's with the van der Waals segment's
/// quoted error added, against the halves' difference and the BAR's quoted error; and log P with
/// the charge segment the two seeds' mean, its error `max(|d|/2, √(e₁² + e₂²)/2)` — the seed
/// spread's standard error of the mean, or the BAR's, the larger — the van der Waals segment and
/// the water leg read from `long/`. **One difference is one degree of freedom**: the spread it
/// implies is known to a factor of a few. Nothing is asserted but that the segments add to the
/// leg.
///
/// About 4.5 h: `long/`'s three windows took 5 328, 5 342 and 5 566 s.
/// `PANTOMETRY_OCTANOL_DIR=C:/pantometry-w4/octanol cargo test --release --locked -p
/// pantometry-forcefield --test benzene_partitioned_into_octanol
/// the_octanol_charges_reseeded_measured -- --ignored --nocapture`.
///
/// **Run once**, at `86b5eb3`, exit 0 on 2026-10-24 after 15 763 s. The charge segment is +0.793 ±
/// 0.111 from long/ and +0.410 ± 0.069 reseeded: 0.384 apart, 2.9 of their quoted errors. With the
/// two seeds' mean, log P is 1.92 ± 0.23 from the direct legs and 2.14 ± 0.23 with the staged
/// water, against 2.13. Nothing here asserts it.
#[test]
#[ignore = "three octanol windows at W3's protocol: about 4.5 h with --release, run once (2026-10-24)"]
fn the_octanol_charges_reseeded_measured() {
    let root = results_dir();
    let t = std::time::Instant::now();
    let sys = octanol_leg();
    let leg_dir = root.join(sys.name);
    let (long_dir, dir) = (leg_dir.join("long"), leg_dir.join(RESEED_TAG));
    assert_ne!(long_dir, dir, "the reseed never writes into long/");
    let _lock = Lock::take(&leg_dir, RESEED_TAG);
    let start = read_back_start(&sys, &leg_dir);
    let long_p = seeded(LONG, sys.name);
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
    assert!(
        schedule.get(..3) == Some(&CHARGE_SEGMENT[..]),
        "long/'s schedule {schedule:?} does not begin with the charge segment {CHARGE_SEGMENT:?}"
    );
    let long_records: Vec<Recorded> = schedule
        .iter()
        .map(|&g| read_back(&long, &long_dir, g))
        .collect();
    let water = water_leg();
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
    assert!(water_records.len() >= 2, "water/long/ is a finished leg");

    std::fs::create_dir_all(&dir).expect("the reseed's directory");
    let p = reseeded(long_p);
    let leg = Leg::new(sys.name, &sys.decoupling, &sys.template, start, p, grid());
    log(
        &dir,
        RESEED_TAG,
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

    // The segments add to the leg: the same intervals' BARs summed in another order.
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
    let to_kcal = |x: f64| kcal(kt * x);
    let mut lines = Vec::new();
    for (name, rs, seed) in [
        ("long/", &long_records[..3], long_p.seed),
        ("reseed/", &records[..], p.seed),
    ] {
        let (g, e) = segment(rs, kt);
        let h = halves_of(rs, kt);
        let mut line = format!(
            "{name} (seed {seed:#x}): the charge segment's ΔG {g:+.3} ± {e:.3} kcal/mol, halves \
             {:+.3}, {:+.3};",
            h[0], h[1]
        );
        for k in 0..rs.len() {
            let tau = Estimate::of(&rs[k].gradients.iter().map(|x| x[1]).collect::<Vec<_>>()).tau;
            line.push_str(&format!(" candidate {} τ {tau:.1}", rs[k].grid));
            if k + 1 < rs.len() {
                let b = interval_bar(&rs[k], &rs[k + 1], kt);
                line.push_str(&format!(
                    ", → {} {:+.3} ± {:.3} overlap {:.3};",
                    rs[k + 1].grid,
                    to_kcal(b.delta),
                    to_kcal(b.error()),
                    b.overlap
                ));
            }
        }
        lines.push(line);
    }
    let (a, ea) = segment(&long_records[..3], kt);
    let (b, eb) = segment(&records, kt);
    let (v, ev) = segment(&long_records[2..], kt);
    let d = b - a;
    let combined = (ea * ea + eb * eb).sqrt();
    let single = d.abs() / std::f64::consts::SQRT_2;
    let leg_spread = (single * single + ev * ev).sqrt();
    let leg_halves = halves_of(&long_records, kt);
    let quoted = to_kcal(whole.error);
    lines.push(format!(
        "the charge segment's seeds differ by {d:+.3} kcal/mol, {:.2} of their combined quoted \
         error {combined:.3}; a single run's spread from it |d|/√2 = {single:.3}, one degree of \
         freedom",
        d.abs() / combined
    ));
    lines.push(format!(
        "the whole octanol leg: long/'s {:+.3} ± {quoted:.3} (BAR, quoted), halves {:+.3}, {:+.3}, \
         differing by {:.3}, which implies a single run's spread of half that, {:.3}; the van der \
         Waals segment {v:+.3} ± {ev:.3} (one seed, its own spread not measured); the seeds' \
         spread with it, √((|d|/√2)² + e_v²) = {leg_spread:.3}",
        to_kcal(whole.value),
        leg_halves[0],
        leg_halves[1],
        (leg_halves[0] - leg_halves[1]).abs(),
        0.5 * (leg_halves[0] - leg_halves[1]).abs()
    ));
    // log P with the charge segment the two seeds' mean.
    let mean = 0.5 * (a + b);
    let e_mean = (0.5 * d.abs()).max(0.5 * combined);
    let g_oct = -(mean + v);
    let e_oct = (e_mean * e_mean + ev * ev).sqrt();
    let (w, ew) = segment(&water_records, kt);
    let g_hyd = -w;
    let rt_ln10 = kcal(BOLTZMANN.to_si() * KELVIN) * std::f64::consts::LN_10;
    let log_p = |g_oct: f64, hyd: f64| -(g_oct - hyd) / rt_ln10;
    let err = |e_oct: f64, e_hyd: f64| (e_oct * e_oct + e_hyd * e_hyd).sqrt() / rt_ln10;
    let e_long = quoted;
    lines.push(format!(
        "long/ alone: ΔG_solv,oct {:+.3} ± {e_long:.3}, ΔG_hyd {g_hyd:+.3} ± {ew:.3}, log P {:.2} \
         ± {:.2}",
        -to_kcal(whole.value),
        log_p(-to_kcal(whole.value), g_hyd),
        err(e_long, ew)
    ));
    lines.push(format!(
        "with the two seeds: the charge segment {mean:+.3} ± {e_mean:.3} (the larger of |d|/2 = \
         {:.3} and the BAR's {:.3}); ΔG_solv,oct {g_oct:+.3} ± {e_oct:.3}; log P {:.2} ± {:.2} \
         from the direct legs, {:.2} ± {:.2} with the staged water's {:+.2} ± {:.2}; experiment \
         {LOG_P}; {:.0} s",
        0.5 * d.abs(),
        0.5 * combined,
        log_p(g_oct, g_hyd),
        err(e_oct, ew),
        log_p(g_oct, HYDRATION_STAGED.0),
        err(e_oct, HYDRATION_STAGED.1),
        HYDRATION_STAGED.0,
        HYDRATION_STAGED.1,
        t.elapsed().as_secs_f64()
    ));
    for line in lines {
        log(&dir, RESEED_TAG, &line);
    }
}

/// **The melts and the pilot** (release, one core, in the background): each leg melted into
/// `NAME/start.txt` under `PANTOMETRY_OCTANOL_DIR` — the production's own start, read back by the
/// legs — with `U`, `P` and `T` logged in 2.5 ps blocks (`NAME/melt.log`); then three windows of
/// [`PILOT`] at candidates 0, 8 and 22 (the charges on, the van der Waals on with no charge, and
/// `λ_v` 0.3 in the soft core) in `NAME/pilot/`, and each one's τ, in samples of 0.1 ps, of its
/// `∂U/∂λ` and of its `ΔU` to the candidates either side.
#[test]
#[ignore = "both boxes melted and three 32 ps windows of each: about two hours with --release"]
fn the_pilot_measured() {
    let root = results_dir();
    for sys in [water_leg(), octanol_leg()] {
        let leg_dir = root.join(sys.name);
        std::fs::create_dir_all(&leg_dir).expect("the leg's directory");
        let _lock = Lock::take(&leg_dir, "pilot");
        let start = melted(&sys, &leg_dir, MELT_SEED);
        let dir = leg_dir.join("pilot");
        std::fs::create_dir_all(&dir).expect("the pilot's directory");
        let p = seeded(PILOT, sys.name);
        let kt = p.kt();
        let leg = Leg::new(sys.name, &sys.decoupling, &sys.template, start, p, grid());
        log(
            &dir,
            "pilot",
            &format!(
                "\n{}: τ in samples of 0.1 ps\n| candidate | λ_e | λ_v | ⟨∂U/∂λ⟩ (kcal/mol) | its τ | τ of βΔU to g−1 | to g+1 | ms a step |\n| --- | --- | --- | --- | --- | --- | --- | --- |",
                sys.name
            ),
        );
        for g in [0usize, 8, 22] {
            let r = window(&leg, &dir, g);
            let c = if g == 0 { 1 } else { 2 };
            let grad = Estimate::of(&r.gradients.iter().map(|x| x[c]).collect::<Vec<_>>());
            let tau_to = |h: Option<usize>| {
                h.filter(|&h| h < leg.grid.len())
                    .map_or("".to_string(), |h| {
                        format!("{:.1}", Estimate::of(&r.delta(h, kt)).tau)
                    })
            };
            log(
                &dir,
                "pilot",
                &format!(
                    "| {g} | {} | {} | {:+.3} ± {:.3} (∂/∂λ_{}) | {:.1} | {} | {} | {:.1} |",
                    r.lambda.electrostatics,
                    r.lambda.van_der_waals,
                    kcal(grad.mean),
                    kcal(grad.error),
                    if c == 1 { "e" } else { "v" },
                    grad.tau,
                    tau_to(g.checked_sub(1)),
                    tau_to(Some(g + 1)),
                    1e3 * r.seconds / p.steps_per_window() as f64
                ),
            );
        }
    }
}

/// **The cost** (release, one core, a few minutes; writes nothing): each production box as built,
/// 200 steps at 0.5 fs at `λ_v` 0.5 with no charge off the lattice, then 300 steps at 2 fs timed
/// there, and the couplings at all 29 candidates timed twenty times; the production at both
/// protocols projected from them, thirteen windows and seventeen, with the melt.
#[test]
#[ignore = "both production boxes built and stepped a few hundred times: minutes with --release"]
fn the_cost_measured() {
    let mut hours = Vec::new();
    for sys in [water_leg(), octanol_leg()] {
        let lambda = Lambda::new(0.0, 0.0, 0.5);
        let potential = AtLambda {
            hamiltonian: &sys.decoupling,
            lambda,
        };
        let mut at = sys.built.clone();
        let mut md = sys
            .template
            .clone()
            .with_bath(Bath::Langevin {
                temperature: KELVIN,
                friction: 50e12,
                seed: 1,
            })
            .thermalised(&at, KELVIN, 1);
        md.run(&potential, &mut at, 0.5 * FS, 200);
        let mut md = md.with_bath(Bath::Langevin {
            temperature: KELVIN,
            friction: 1e12,
            seed: 2,
        });
        let steps = 300;
        let t = std::time::Instant::now();
        md.run(&potential, &mut at, 2.0 * FS, steps);
        let step = t.elapsed().as_secs_f64() / steps as f64;
        let states = grid();
        let t = std::time::Instant::now();
        for _ in 0..20 {
            std::hint::black_box(sys.decoupling.couplings(&at, &states));
        }
        let couplings = t.elapsed().as_secs_f64() / 20.0;
        let melt_steps: f64 = sys
            .melt
            .iter()
            .map(|s| (s.picoseconds * 1e-12 / s.dt).round())
            .sum();
        println!(
            "{}: {} atoms, box {:.2?} Å; a step {:.1} ms, the couplings at all {} candidates {:.1} ms; \
             the melt {melt_steps} steps, {:.2} h",
            sys.name,
            at.len(),
            sys.field.cell().lengths().map(|x| x / ANGSTROM),
            1e3 * step,
            states.len(),
            1e3 * couplings,
            melt_steps * step / 3600.0
        );
        for (name, p) in [("LONG", LONG), ("SHORT", SHORT)] {
            let window = p.steps_per_window() as f64 * step + p.samples as f64 * couplings;
            println!(
                "  {name}: a window {:.2} h; 13 windows {:.1} h, 17 {:.1} h",
                window / 3600.0,
                13.0 * window / 3600.0,
                17.0 * window / 3600.0
            );
            hours.push((sys.name, name, window, melt_steps * step));
        }
    }
    for p in ["LONG", "SHORT"] {
        let total: f64 = hours
            .iter()
            .filter(|h| h.1 == p)
            .map(|h| 13.0 * h.2 + h.3)
            .sum();
        let most: f64 = hours
            .iter()
            .filter(|h| h.1 == p)
            .map(|h| 17.0 * h.2 + h.3)
            .sum();
        println!(
            "{p}, both legs with their melts: {:.1} h at thirteen windows each, {:.1} h at seventeen",
            total / 3600.0,
            most / 3600.0
        );
    }
}

/// The digest is of the model: every number in it moves it. Each of benzene's `x`, an octanol's
/// charge, the mixing and the 1-4 scale changed alone changes [`model_digest`].
#[test]
fn the_model_digest_sees_the_model() {
    let s = small_water(0xD16);
    let f = s.field.clone();
    let d = model_digest(&f);
    let mut v = f.van_der_waals().to_vec();
    v[0][0] = f64::from_bits(v[0][0].to_bits() + 1);
    let n = v.len();
    let mut q = f.charges().to_vec();
    q[20] = f64::from_bits(q[20].to_bits() + 1);
    for (what, g) in [
        ("an x", f.clone().with_nonbonded(0..n, v).expect("valid")),
        ("a charge", f.clone().with_charges(q)),
        ("the mixing", f.clone().with_mixing(Mixing::Geometric)),
        ("the 1-4 scale", f.clone().with_one_four(0.5, 1.0)),
    ] {
        assert_ne!(model_digest(&g), d, "{what}");
    }
}
