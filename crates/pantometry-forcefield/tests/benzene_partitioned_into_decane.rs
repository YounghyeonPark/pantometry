//! **Benzene decoupled from n-decane, for its alkane:water log P**: step C-2 of the explicit-solvent
//! track. Benzene is the octanol and water legs' model to the bit — GAFF's `ca` and `ha` by
//! Lorentz–Berthelot, GAFF's 1-4 scales `(0.5, 1/1.2)`, FreeSolv's AM1-BCC charges, on UFF's bonded
//! terms (`tests/partition/`). The decane is C-1's: GAFF's `c3` and `hc`, the AM1-BCC charges
//! transferred by position (`tests/hexadecane/`), in C-1b's isotropic start.
//!
//! **Kept out of the published crate** (`Cargo.toml`'s `exclude`), with `tests/gaff/`,
//! `tests/hexadecane/` and `tests/partition/`: GAFF's numbers are read from `gaff.dat`, GPL v3.
//!
//! # The thermodynamics
//!
//! `log P(decane/water) = −(ΔG_solv,dec − ΔG_hyd) / (RT ln 10)`, each `ΔG_solv = −ΔG_decouple`,
//! benzene decoupled W3's way in each box, so that the decoupled state is the solvent alone in its
//! box and benzene alone in vacuum with its own terms. Benzene's own terms are the same Hamiltonian at
//! both ends of both legs and cancel. Both legs are at fixed volume and the same concentration at
//! each end, so the standard state cancels.
//!
//! **The water leg is read, not run**: 4c-2's, in `PANTOMETRY_OCTANOL_DIR/water/` — its `start.txt`
//! and its `long/` windows through their own headers, so a window run under another model, from
//! another start or at another seed is refused by name, and a missing one is refused rather than
//! run. Its ΔG_hyd is −1.130 ± 0.144 kcal/mol direct; the staged route's −0.82 ± 0.15 is printed
//! beside it, as the octanol leg printed it.
//!
//! **What it is compared against**, both read in the design study and neither decane:
//!
//! - **n-hexadecane:water, 2.14**, at 25 °C: benzene's `log L¹⁶ = 2.786` (its hexadecane:air
//!   partition, Abraham's descriptor as Stolov et al. 2014 give it) plus `log K_aw`, the air:water
//!   partition from its Henry's constant (PubChem CID 241's HSDB record), `K_aw = 0.2273`:
//!   `2.786 + log₁₀ 0.2273 = 2.143` ([`HEXADECANE_WATER`]).
//! - **olive oil:saline, 2.23**, at 37 °C (M. L. Gargas et al., *Toxicol. Appl. Pharmacol.* 98, 87
//!   (1989)).
//!
//! **No measured value for decane itself was read.** Linear alkanes from hexane to hexadecane
//! partition benzene alike to within a few tenths, which is why hexadecane's is the comparison; it is
//! not decane's.
//!
//! # The box
//!
//! C-1b's isotropic start — 112 decane sites, each a conformer with random gauche defects in a
//! uniformly random orientation, packed apart without turning ([`hexadecane::isotropic`]) — with
//! benzene at the box's centre, and every decane with an atom nearer to a benzene atom than their
//! packing radii' sum removed, so that what is left is overlap-free at the packer's own criterion.
//! **The box then shrunk so that the decanes are at 0.7255 g/cm³ in the volume benzene leaves
//! them**, `N_dec m / (V − v_B)`, `v_B` 89 cm³/mol as 4c-2's: every edge and every molecule's
//! centre scaled by `(V'/V)^(1/3)`, each molecule rigid and unturned, and the packer run again with
//! benzene held ([`decane_box`]). 4c-2's octanol solved a fixed point in the number removed instead;
//! here rebuilding the start at another density repacks it from scratch, and the number removed
//! flipped between rounds and never settled. Uncompensated, one decane removed of 112 for benzene's
//! 89 cm³/mol leaves the 111 0.5% light. `r_c` 9 Å, PME at δ = 10⁻⁵ with the neighbour
//! list, every X–H held.
//!
//! # The melt, and its gate
//!
//! [`DECANE_MELT`]: 1 ps at 0.5 fs in a 50 ps⁻¹ bath, which takes the packing's strain out (C-1b's);
//! 5 ps at 600 K, which C-1b measured to be more than enough (its order was in the isotropic band
//! from the first 2.5 ps block at 600 K, and `U` flat); and **200 ps at 298.15 K**, where C-1b's `U`
//! settled in about 20 ps and its trans fraction had not settled in 200. `U`, `P`, the nematic order
//! and the trans fraction every 0.1 ps, logged in 2.5 ps blocks and kept in `melt.txt` beside
//! `start.txt`.
//!
//! **The windows run only if the settle passes a measured criterion** ([`settle_verdict`]), read
//! from `melt.txt` every time the start is used, so a rerun cannot step past a failed one. Over the
//! settle less its first 20 ps:
//!
//! - **the halves of `U` and of `P` agree within three of their combined errors** (each half's error
//!   from its own autocorrelation time). Three is the threshold C-1b's reading drew: its `U/N` halves
//!   1.8 σ apart were called settled and its trans fraction's 3.2 σ not.
//! - **the order's halves agree the same way**, so it is not rising; and a rising order is refused by
//!   that name;
//! - **the order's second half is at most the isotropic band plus [`ORDER_ALLOWANCE`] = 0.04**. The
//!   band is an ideal random set's mean plus three standard deviations, about 0.157 for 111 decanes;
//!   a liquid of rods has short-range orientational correlation that a random set has not, and
//!   C-1b's liquid measured 0.157 ± 0.007. So about 0.20 passes a liquid like C-1b's with six of its
//!   standard errors to spare, and refuses any residual alignment: Maier–Saupe's order at the
//!   isotropic–nematic transition is 0.43, and C-1's lattice-aligned decane was 0.57.
//!
//! The trans fraction is reported with its halves and not gated: C-1b found it still rising after
//! 200 ps, and the brief gates the three above.
//!
//! # The windows
//!
//! 4c-2's: the 29 candidates, the thirteen windows of [`SCHEDULE`], a window inserted between any
//! pair overlapping below 0.1, twice; [`LONG`] at 20 + 200 ps, with this leg's own seed
//! ([`DECANE_SEED`]). Resumable and headed through `tests/partition/`'s `window`: each window
//! `window_NN.partial` until whole, its header the candidate, λ, the protocol, the atom count, the
//! candidates and digests of the start and of the model, and a lock. **The nematic order every
//! sample**, as a diagnostic: the end-to-end vectors' order is written beside each sample and
//! reported per window, and nothing asserts it.
//!
//! Written under `PANTOMETRY_DECANE_DIR` (`start.txt`, `melt.txt`, `melt.log`, `long/`); nothing is
//! written under `PANTOMETRY_OCTANOL_DIR`.
//!
//! # Checks
//!
//! Default, unoptimised, on a small decane box of twelve sites at `r_c` 6.5 Å unless named: both
//! models are what they say, and benzene's is the octanol and water legs' to the bit; the decane is
//! at its density with benzene in, and nothing overlaps, in the production's box; the end states are
//! the systems they say; the λ derivatives are the whole energy's; every state of `couplings` is its
//! own `coupling`; the decoupled reciprocal energy is the rest's plus `λ_e` times the cross term; the
//! settle's gate on closed cases; and a killed leg resumes to the same bits and refuses another
//! model, seed or start.

// The sums over atoms and axes are written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

mod freesolv;
mod gaff;
mod hexadecane;
mod hydration;
mod nonbonded;
mod octanol;
mod partition;

use gaff::{gaff, Model, ONE_FOUR};
use hexadecane::{
    conformations, directions, isotropic, isotropic_order, nematic_order, pack, packing_radius,
    worst_contact, Alkane, D10, DECANE, START_TRANS,
};
use hydration::{Benzene, Solvated};
use nonbonded::Vacuum;
use octanol::{joined, production, Liquid, ANGSTROM, FS, KELVIN};
use pantometry_forcefield::free_energy::refine_schedule;
use pantometry_forcefield::{
    qeq, Alchemical, AtLambda, Bath, Component, Element, Estimate, EwaldParameters, Lambda, Mixing,
    PeriodicBox, PeriodicDecoupling, PeriodicForceField, PmeParameters, Potential, Protocol,
};
use pantometry_units::BOLTZMANN;
use partition::*;
use std::ops::Range;

/// The small decane box: one site along the chain, three by four across, 15.4 × 15.4 × 16.5 Å.
const SMALL: [usize; 3] = [1, 3, 4];
/// Its cutoff, under half its shortest edge.
const SMALL_CUTOFF: f64 = 6.5 * ANGSTROM;

/// C-1b's start seed: the production's box is built from it.
const START_SEED: u64 = 0xC1B;

/// The packer's budget, C-1b's.
const PACKING_ITERATIONS: usize = 20_000;

/// **This leg's campaign seed**: window `g` is seeded `window_seed(DECANE_SEED, g)`. Not 4c-2's
/// water `0x4C2_3A7E`, octanol `0x4C2_0C71`, reseed `0x4C2_5EED` or melt `0x4C2_E01`.
const DECANE_SEED: u64 = 0xC2_DEC;

/// The decane melt's seed, stage `k` seeded `window_seed(DECANE_MELT_SEED, k)`.
const DECANE_MELT_SEED: u64 = 0xC2_E01;

/// How often the melt records, seconds: every 0.1 ps, as 4c-2's melts.
const RECORD_EVERY: f64 = 0.1e-12;

/// **The decane leg's melt**: see the module documentation.
const DECANE_MELT: [MeltStage; 3] = [
    stage("strain out", Lambda::COUPLED, 1.0, 0.5 * FS, 50e12, KELVIN),
    stage("hot", Lambda::COUPLED, 5.0, 2.0 * FS, 5e12, 600.0),
    stage("settle", Lambda::COUPLED, 200.0, 2.0 * FS, 5e12, KELVIN),
];

/// How much of the settle the gate leaves out, picoseconds: C-1b's `U` settled in about 20 ps.
const SETTLE_DISCARD: f64 = 20.0;

/// How far apart two halves may be, in their combined errors, and still agree: see the module
/// documentation.
const HALVES_AGREE: f64 = 3.0;

/// **How far above the isotropic band the order may sit**: see the module documentation.
const ORDER_ALLOWANCE: f64 = 0.04;

/// **n-Hexadecane:water log P for benzene, 25 °C**: `log L¹⁶ + log₁₀ K_aw`, 2.786 and 0.2273 (see the
/// module documentation), 2.143.
fn hexadecane_water() -> f64 {
    2.786 + 0.2273f64.log10()
}

/// The figure [`hexadecane_water`] is quoted as.
const HEXADECANE_WATER: f64 = 2.14;

/// **Olive oil:saline log P for benzene, 37 °C**: Gargas et al. (1989).
const OLIVE_OIL_SALINE: f64 = 2.23;

fn len(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

// ---------------------------------------------------------------------------------------------
// The box

/// What [`decane_box`] built: the box with benzene first, its GAFF field (classical Ewald), the
/// density all its sites would have in the box it was scaled to, the linear scale, how many decanes
/// benzene displaced, every atom's packing radius, and the packer's iterations after the scaling.
struct DecaneBox {
    liquid: Liquid,
    field: PeriodicForceField,
    lattice: f64,
    scale: f64,
    removed: usize,
    radii: Vec<f64>,
    repacked: usize,
}

/// **Benzene in C-1b's isotropic decane start**: [`isotropic`]'s start of `counts` sites from
/// `seed` at [`DECANE`]'s density, exactly C-1b's; benzene at the box's centre, and every decane with
/// an atom nearer to a benzene atom than their packing radii' sum by minimum image removed. GAFF
/// everywhere — benzene's [`Model`], the decanes' `c3`/`hc` and their transferred charges — by
/// Lorentz–Berthelot with GAFF's 1-4 scales.
///
/// With `compensate`, **the box then shrunk so that the decanes kept are at [`DECANE`]'s density in
/// the volume benzene leaves them**: `V' = N_dec m / ρ + v_B`, every edge scaled by `(V'/V)^(1/3)`
/// and every molecule's centre with it, each molecule rigid and unturned, then [`pack`]ed again with
/// benzene held, so that no contact is nearer than the packer's radii. The number removed is decided
/// once, on C-1b's start, and the volume follows from it: a fixed point in the number removed, as
/// 4c-2's octanol lattice has, does not exist here — rebuilding the start at another density
/// repacks it from scratch, and measured on the production's box the number removed then flipped
/// between rounds and never settled in ten.
fn decane_box(
    b: &Benzene,
    a: &Alkane,
    counts: [usize; 3],
    cutoff: f64,
    accuracy: f64,
    seed: u64,
    compensate: bool,
) -> DecaneBox {
    let sites = counts[0] * counts[1] * counts[2];
    let o = &a.molecule;
    let brad: Vec<f64> = b
        .component
        .atoms()
        .iter()
        .map(|x| packing_radius(x.element))
        .collect();
    let s = isotropic(
        a,
        counts,
        DECANE.density,
        seed,
        START_TRANS,
        cutoff,
        accuracy,
        PACKING_ITERATIONS,
    );
    let l = &s.liquid;
    let cell = l.cell();
    let bat = b.centred_at(cell.lengths().map(|x| 0.5 * x));
    let keep: Vec<usize> = (0..l.octanols)
        .filter(|&m| {
            !l.molecules[m].clone().any(|i| {
                bat.iter()
                    .zip(&brad)
                    .any(|(q, rq)| len(cell.minimum_image(sub(l.at[i], *q))) < s.radii[i] + rq)
            })
        })
        .collect();
    let volume = cell.volume();
    let wanted = keep.len() as f64 * o.mass() / DECANE.density + BENZENE_VOLUME;
    let scale = if compensate {
        (wanted / volume).cbrt()
    } else {
        1.0
    };
    let edges = cell.lengths().map(|x| scale * x);
    let small = PeriodicBox::new(edges);
    // Benzene and the kept decanes as rigid bodies about their centres, the decanes' centres
    // scaled with the box and benzene's at its new centre, held.
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut bodies = vec![b.centred_at([0.0; 3])];
    let mut radii_of = vec![brad.clone()];
    let mut centres = vec![edges.map(|x| 0.5 * x)];
    let mut rotations = vec![identity];
    for &m in &keep {
        let r = l.molecules[m].clone();
        let com = [0, 1, 2].map(|k| {
            r.clone()
                .map(|i| o.masses[i - r.start] * l.at[i][k])
                .sum::<f64>()
                / o.mass()
        });
        bodies.push(s.bodies[m].clone());
        radii_of.push(r.clone().map(|i| s.radii[i]).collect());
        centres.push(com.map(|x| scale * x));
        rotations.push(s.rotations[m]);
    }
    let mut fixed = vec![false; bodies.len()];
    fixed[0] = true;
    let packed = pack(
        &bodies,
        &radii_of,
        centres,
        &rotations,
        &fixed,
        &small,
        seed ^ 0xC2,
        PACKING_ITERATIONS,
    );
    assert!(
        packed.worst >= 1.0,
        "the scaled box left an overlap at {} of the radii after {} iterations",
        packed.worst,
        packed.iterations
    );
    // Each decane moved by whole edges so that its centre is in the box.
    let mut at: Vec<[f64; 3]> = Vec::new();
    let mut molecules: Vec<Range<usize>> = Vec::new();
    let mut offset = 0;
    for (k, body) in bodies.iter().enumerate() {
        let c = packed.centres[k];
        let shift = sub(small.wrap(c), c);
        molecules.push(offset..offset + body.len());
        for i in 0..body.len() {
            let p = packed.at[offset + i];
            at.push([p[0] + shift[0], p[1] + shift[1], p[2] + shift[2]]);
        }
        offset += body.len();
    }
    let bm = Model::benzene(&b.component);
    let decane_charges = a.transfer.charges.clone();
    let mut parts: Vec<&Component> = vec![&b.component];
    let mut types = b.types.clone();
    let mut charges = bm.charges.clone();
    let mut masses = b.masses.clone();
    let mut vdw = bm.vdw.clone();
    for _ in &keep {
        parts.push(&o.component);
        types.extend(o.types.iter().cloned());
        charges.extend_from_slice(&decane_charges);
        masses.extend_from_slice(&o.masses);
        vdw.extend(a.vdw());
    }
    let component = joined(&parts);
    let p = EwaldParameters::for_accuracy(&small, cutoff, accuracy);
    let n = vdw.len();
    let field = PeriodicForceField::new(&component, &types, small, p)
        .expect("benzene in decane")
        .with_charges(charges.clone())
        .with_nonbonded(0..n, vdw)
        .expect("GAFF's values are valid")
        .with_mixing(Mixing::LorentzBerthelot)
        .with_one_four(ONE_FOUR[0], ONE_FOUR[1]);
    let liquid = Liquid {
        component,
        types,
        charges,
        field: field.clone(),
        at,
        masses,
        molecules,
        solute: b.at.len(),
        octanols: keep.len(),
        site: l.site.map(|x| scale * x),
        gap: l.gap,
    };
    DecaneBox {
        removed: sites - keep.len(),
        liquid,
        field,
        lattice: sites as f64 * o.mass() / small.volume(),
        scale,
        radii: radii_of.concat(),
        repacked: packed.iterations,
    }
}

/// [`System`] for benzene in decane: `d` under `field`, the decoupling of benzene's twelve atoms.
fn decane_system(d: &DecaneBox, field: PeriodicForceField, melt: &'static [MeltStage]) -> System {
    let l = &d.liquid;
    let group: Vec<bool> = (0..l.at.len()).map(|i| i < l.solute).collect();
    System {
        name: "decane",
        decoupling: PeriodicDecoupling::new(&field, &group).expect("benzene decouples"),
        field,
        template: l.dynamics(),
        built: l.at.clone(),
        molecules: l.molecules.clone(),
        masses: l.masses.clone(),
        solvent: l.octanols,
        melt,
    }
}

/// The production's decane box: [`DECANE`]'s 112 sites from C-1b's seed, `r_c` 9 Å, compensated.
fn production_box(b: &Benzene, a: &Alkane) -> DecaneBox {
    decane_box(b, a, DECANE.lattice, 9.0 * ANGSTROM, 1e-5, START_SEED, true)
}

/// The production's decane leg: [`production_box`] with PME at δ = 10⁻⁵ and the neighbour list.
fn decane_leg(b: &Benzene, a: &Alkane) -> (System, DecaneBox) {
    let d = production_box(b, a);
    let field = production(d.field.clone(), &d.liquid);
    (decane_system(&d, field, &DECANE_MELT), d)
}

/// What the decane is judged by, on a box whose molecule 0 is benzene: the decanes' end-to-end
/// order and their trans fraction.
struct Watch {
    solvent: Vec<Range<usize>>,
    ends: [usize; 2],
    carbons: Vec<usize>,
}

impl Watch {
    fn new(sys: &System, a: &Alkane) -> Watch {
        Watch {
            solvent: sys.molecules[1..].to_vec(),
            ends: a.ends(),
            carbons: a.carbons.clone(),
        }
    }

    /// The nematic order of the decanes' end-to-end vectors.
    fn order(&self, at: &[[f64; 3]]) -> f64 {
        nematic_order(&directions(at, &self.solvent, self.ends))
    }

    /// The fraction of their C–C–C–C dihedrals that are trans.
    fn trans(&self, at: &[[f64; 3]]) -> f64 {
        conformations(at, &self.solvent, &self.carbons).0
    }

    /// The isotropic band for this many decanes: an ideal random set's mean and standard deviation
    /// over a thousand sets, C-1b's.
    fn band(&self) -> (f64, f64) {
        isotropic_order(self.solvent.len(), 1000, 0x150)
    }
}

// ---------------------------------------------------------------------------------------------
// The melt

/// One record of the melt: its stage, and the time (ps), `U` per decane (kcal/mol), `P` (bar), `T`
/// (K), the nematic order and the trans fraction.
type MeltRow = (usize, [f64; 6]);

/// The decane melt's header: [`melt_header`]'s, and how often it records.
fn decane_melt_header(sys: &System, seed: u64, every: f64) -> String {
    format!("{} record {}", melt_header(sys, seed), exact(every))
}

/// The melt's record as text: the header, a row a record, and a trailer `# end N DIGEST`.
fn record_text(header: &str, rows: &[MeltRow]) -> String {
    let mut text = format!("{header}\n");
    for (k, r) in rows {
        text.push_str(&format!("{k}"));
        for x in r {
            text.push(' ');
            text.push_str(&exact(*x));
        }
        text.push('\n');
    }
    text.push_str(&format!(
        "# end {} {:016x}\n",
        rows.len(),
        digest(rows.iter().flat_map(|r| r.1))
    ));
    text
}

/// A melt record read back, or why it is not whole or not this melt's.
fn parse_record(text: &str, header: &str) -> Result<Vec<MeltRow>, String> {
    let mut lines = text.lines();
    let head = lines.next().ok_or("an empty record")?;
    if head != header {
        return Err("it was melted differently".into());
    }
    let mut rows = Vec::new();
    for line in lines {
        if let Some(end) = line.strip_prefix("# end ") {
            let w: Vec<&str> = end.split_whitespace().collect();
            let n: Option<usize> = w.first().and_then(|x| x.parse().ok());
            let d = w.get(1).and_then(|x| u64::from_str_radix(x, 16).ok());
            return match (n, d) {
                (Some(n), Some(d))
                    if n == rows.len() && d == digest(rows.iter().flat_map(|r: &MeltRow| r.1)) =>
                {
                    Ok(rows)
                }
                _ => Err(format!("its trailer is not its rows': {line:?}")),
            };
        }
        let w: Vec<&str> = line.split_whitespace().collect();
        let k: Option<usize> = w.first().and_then(|x| x.parse().ok());
        let v: Result<Vec<f64>, _> = w.iter().skip(1).map(|x| x.parse::<f64>()).collect();
        match (k, v) {
            (Some(k), Ok(v)) if v.len() == 6 => {
                rows.push((k, [v[0], v[1], v[2], v[3], v[4], v[5]]))
            }
            _ => return Err(format!("not a record: {line:?}")),
        }
    }
    Err("no trailer: the record is not whole".into())
}

/// **`sys` melted from its positions as built**, with its record: read back from `start.txt` and
/// `melt.txt` in `dir` when a run wrote them under the same header ([`decane_melt_header`]), a file
/// melted differently or not whole refused by name; otherwise run now under a lock `melt` — each
/// stage from the last, the bath thermalised once at the start, every `every` seconds `U` per
/// decane, `P`, `T`, the order and the trans fraction recorded — and written atomically, the record
/// first. The record is logged to `melt.log` in 2.5 ps blocks.
fn melted_decane(
    sys: &System,
    dir: &std::path::Path,
    seed: u64,
    watch: &Watch,
    every: f64,
) -> (Vec<[f64; 3]>, Vec<MeltRow>) {
    let start_path = dir.join("start.txt");
    let record_path = dir.join("melt.txt");
    let want = decane_melt_header(sys, seed, every);
    if let Ok(text) = std::fs::read_to_string(&start_path) {
        let (head, at) = parse_positions(&text).unwrap_or_else(|why| {
            panic!(
                "{} is not whole ({why}): remove it and melt.txt to melt again",
                start_path.display()
            )
        });
        assert_eq!(
            head,
            want,
            "{} was melted differently: remove it, melt.txt and every window run from it, to melt \
             again",
            start_path.display()
        );
        assert_eq!(at.len(), sys.built.len());
        let rows = std::fs::read_to_string(&record_path)
            .map_err(|e| e.to_string())
            .and_then(|t| parse_record(&t, &want))
            .unwrap_or_else(|why| {
                panic!(
                    "{} is not whole ({why}): remove start.txt to melt again",
                    record_path.display()
                )
            });
        return (at, rows);
    }
    let _lock = Lock::take(dir, "melt");
    let t = std::time::Instant::now();
    let mut at = sys.built.clone();
    let mut md = sys.template.clone();
    let n = sys.solvent as f64;
    let mut rows: Vec<MeltRow> = Vec::new();
    let mut time = 0.0;
    for (k, st) in sys.melt.iter().enumerate() {
        let s = pantometry_forcefield::free_energy::window_seed(seed, k);
        md = md.with_bath(Bath::Langevin {
            temperature: st.temperature,
            friction: st.friction,
            seed: s,
        });
        if k == 0 {
            md = md.thermalised(&at, st.temperature, s);
        }
        md.forget_forces();
        let potential = AtLambda {
            hamiltonian: &sys.decoupling,
            lambda: st.lambda,
        };
        let steps = (st.picoseconds * 1e-12 / st.dt).round() as usize;
        let stride = ((every / st.dt).round() as usize).max(1);
        let t0 = std::time::Instant::now();
        for i in 0..steps {
            md.step(&potential, &mut at, st.dt);
            if (i + 1) % stride == 0 {
                time += stride as f64 * st.dt;
                rows.push((
                    k,
                    [
                        time * 1e12,
                        kcal(md.potential_energy().expect("stepped")) / n,
                        pressure(sys, &at) / 1e5,
                        md.half_step_temperature(),
                        watch.order(&at),
                        watch.trans(&at),
                    ],
                ));
            }
        }
        let secs = t0.elapsed().as_secs_f64();
        log(
            dir,
            "melt",
            &format!(
                "{} stage {k} ({}): {} ps at {} fs, {} K, γ {} /ps: {secs:.0} s, {:.1} ms a step",
                sys.name,
                st.name,
                st.picoseconds,
                st.dt / FS,
                st.temperature,
                st.friction / 1e12,
                1e3 * secs / steps.max(1) as f64
            ),
        );
    }
    assert!(rows.iter().all(|r| r.1.iter().all(|x| x.is_finite())));
    log(
        dir,
        "melt",
        &format!(
            "{}: the melt in blocks of 2.5 ps\n| stage | t (ps) | U/N (kcal/mol) | P (bar) | T (K) | order | trans |\n| --- | --- | --- | --- | --- | --- | --- |",
            sys.name
        ),
    );
    let per_block = ((2.5e-12 / every).round() as usize).max(1);
    for block in rows.chunks(per_block) {
        let m = |c: usize| block.iter().map(|r| r.1[c]).sum::<f64>() / block.len() as f64;
        log(
            dir,
            "melt",
            &format!(
                "| {} | {:.1}–{:.1} | {:.4} | {:.0} | {:.1} | {:.4} | {:.4} |",
                sys.melt[block[block.len() - 1].0].name,
                block[0].1[0] - every * 1e12,
                block[block.len() - 1].1[0],
                m(1),
                m(2),
                m(3),
                m(4),
                m(5)
            ),
        );
    }
    let partial = record_path.with_extension("partial");
    std::fs::write(&partial, record_text(&want, &rows)).expect("the record");
    std::fs::rename(&partial, &record_path).expect("rename");
    write_positions(&start_path, &want, &at);
    log(
        dir,
        "melt",
        &format!("{} melted in {:.0} s", sys.name, t.elapsed().as_secs_f64()),
    );
    (at, rows)
}

/// Two halves of `xs`: each one's estimate, and how far apart their means are in their combined
/// errors.
fn halves(xs: &[f64]) -> (Estimate, Estimate, f64) {
    let h = xs.len() / 2;
    let (a, b) = (Estimate::of(&xs[..h]), Estimate::of(&xs[h..]));
    let apart = (a.mean - b.mean).abs() / (a.error * a.error + b.error * b.error).sqrt();
    (a, b, apart)
}

/// **The settle's gate**: over stage `stage`'s rows less the first `discard`, the halves of `U` and
/// of `P` within [`HALVES_AGREE`] of their combined errors; the order's halves the same, refused as
/// *rising* when its second half is the higher; and the order's second half at most `band` plus
/// [`ORDER_ALLOWANCE`]. The trans fraction reported, not gated. `Ok` with the lines that say why, or
/// `Err` with them and what failed.
fn settle_verdict(
    rows: &[MeltRow],
    stage: usize,
    discard: usize,
    band: f64,
) -> Result<Vec<String>, Vec<String>> {
    let kept: Vec<&[f64; 6]> = rows
        .iter()
        .filter(|r| r.0 == stage)
        .skip(discard)
        .map(|r| &r.1)
        .collect();
    if kept.len() < 8 {
        return Err(vec![format!(
            "the settle has {} records after the {discard} discarded: too few to judge",
            kept.len()
        )]);
    }
    let column = |c: usize| -> Vec<f64> { kept.iter().map(|r| r[c]).collect() };
    let mut lines = vec![format!(
        "the settle's gate over {} records, {:.1}–{:.1} ps: halves within {HALVES_AGREE} of their \
         combined errors, the order's second half at most the band {band:.4} + {ORDER_ALLOWANCE} = \
         {:.4}",
        kept.len(),
        kept[0][0],
        kept[kept.len() - 1][0],
        band + ORDER_ALLOWANCE
    )];
    let mut failed = Vec::new();
    for (c, what, gated) in [
        (1, "U/N (kcal/mol)", true),
        (2, "P (bar)", true),
        (4, "order", true),
        (5, "trans fraction", false),
    ] {
        let (a, b, apart) = halves(&column(c));
        lines.push(format!(
            "  {what}: halves {:.4} ± {:.4} (τ {:.1}) and {:.4} ± {:.4} (τ {:.1}), {apart:.2} σ \
             apart{}",
            a.mean,
            a.error,
            a.tau,
            b.mean,
            b.error,
            b.tau,
            if gated { "" } else { "; reported, not gated" }
        ));
        let agree = apart <= HALVES_AGREE;
        if gated && !agree {
            if c == 4 && b.mean > a.mean {
                failed.push(format!(
                    "the order is rising: {:.4} → {:.4}, {apart:.2} σ",
                    a.mean, b.mean
                ));
            } else {
                failed.push(format!("{what}'s halves disagree, {apart:.2} σ apart"));
            }
        }
        let under = b.mean <= band + ORDER_ALLOWANCE;
        if c == 4 && !under {
            failed.push(format!(
                "the order's second half {:.4} is above the band and its allowance, {:.4}",
                b.mean,
                band + ORDER_ALLOWANCE
            ));
        }
    }
    if failed.is_empty() {
        lines.push("  the settle passes".into());
        Ok(lines)
    } else {
        lines.push(format!("  the settle FAILS: {}", failed.join("; ")));
        Err(lines)
    }
}

// ---------------------------------------------------------------------------------------------
// The windows

/// This leg at protocol `base`: [`LONG`] or another, with [`DECANE_SEED`].
fn decane_protocol(base: Protocol) -> Protocol {
    Protocol {
        seed: DECANE_SEED,
        ..base
    }
}

/// One window of the decane leg, its order observed every sample.
fn decane_window(leg: &Leg, dir: &std::path::Path, g: usize, watch: &Watch) -> Recorded {
    let order = |at: &[[f64; 3]]| watch.order(at);
    window_observed(leg, dir, g, Some(&order))
}

/// One leg's result without its table: BAR along the chain, TI, the halves and the lowest overlap.
fn summarise(records: &[Recorded], kt: f64) -> LegResult {
    let total = chain(records, kt);
    let (ti, ti_error, _) = trapezoid(records);
    let overlap = (0..records.len() - 1)
        .map(|k| interval_bar(&records[k], &records[k + 1], kt).overlap)
        .fold(f64::INFINITY, f64::min);
    LegResult {
        bar: (kcal(kt * total.value), kcal(kt * total.error)),
        ti: (ti, ti_error),
        halves: halves_of(records, kt),
        overlap,
    }
}

/// Where the decane leg writes: `PANTOMETRY_DECANE_DIR`, or `pantometry-decane` in the system's
/// temporary directory. Never the repository.
fn decane_dir() -> std::path::PathBuf {
    let dir = std::env::var_os("PANTOMETRY_DECANE_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("pantometry-decane"));
    std::fs::create_dir_all(&dir).expect("the results directory");
    dir
}

/// **4c-2's water leg, read back and never run**: `PANTOMETRY_OCTANOL_DIR/water/start.txt` and the
/// `long/` windows, each through its own header ([`read_back_start`], [`read_back`]), so that one
/// run under another model, from another start or at another seed is refused by name and a
/// missing one is refused rather than run. Nothing is written there. Its schedule must hold
/// [`SCHEDULE`].
fn water_read_back() -> (Vec<usize>, LegResult) {
    let root = std::path::PathBuf::from(std::env::var_os("PANTOMETRY_OCTANOL_DIR").expect(
        "PANTOMETRY_OCTANOL_DIR names the octanol campaign's directory, whose water/ is read back",
    ));
    let sys = water_leg();
    let dir = root.join(sys.name);
    let start = read_back_start(&sys, &dir);
    let p = seeded(LONG, sys.name);
    let leg = Leg::new(sys.name, &sys.decoupling, &sys.template, start, p, grid());
    let long = dir.join("long");
    let schedule = finished_schedule(&long);
    assert!(
        SCHEDULE.iter().all(|g| schedule.contains(g)),
        "{}'s schedule {schedule:?} does not hold {SCHEDULE:?}",
        long.display()
    );
    let records: Vec<Recorded> = schedule
        .iter()
        .map(|&g| read_back(&leg, &long, g))
        .collect();
    (schedule, summarise(&records, p.kt()))
}

/// **Benzene's log P(decane/water) at `base`**: the water leg read back first ([`water_read_back`]),
/// then the decane leg under `PANTOMETRY_DECANE_DIR` — the melt and its record ([`melted_decane`]),
/// the settle's gate ([`settle_verdict`]), refused and reported if it fails; W3's thirteen windows
/// and the insertion, each window's order observed; the table, the order per window, and log P from
/// the direct legs and with the staged water, against hexadecane:water's 2.14 and olive oil:saline's
/// 2.23. Nothing about the model is asserted.
fn decane_partition_measured(tag: &str, base: Protocol) {
    let t = std::time::Instant::now();
    let (water_schedule, water) = water_read_back();
    let root = decane_dir();
    let _lock = Lock::take(&root, tag);
    log(
        &root,
        tag,
        &format!(
            "the water leg read back from PANTOMETRY_OCTANOL_DIR/water/long/, windows \
             {water_schedule:?}: ΔG_decouple BAR {:+.3} ± {:.3}, TI {:+.3} ± {:.3}, halves {:+.3}, \
             {:+.3}, lowest overlap {:.3}",
            water.bar.0,
            water.bar.1,
            water.ti.0,
            water.ti.1,
            water.halves[0],
            water.halves[1],
            water.overlap
        ),
    );
    let b = benzene();
    let a = Alkane::new(D10);
    let (sys, d) = decane_leg(&b, &a);
    let watch = Watch::new(&sys, &a);
    let (iso, sd) = watch.band();
    let band = iso + 3.0 * sd;
    log(
        &root,
        tag,
        &format!(
            "decane: benzene and {} decanes ({} removed), {} atoms, box {:?} Å, C-1b's start scaled by \
             {:.6} and repacked in {} iterations; an isotropic set of {} gives {iso:.4} ± {sd:.4}, band ≤ {band:.4}; model \
             {:016x}",
            sys.solvent,
            d.removed,
            sys.built.len(),
            sys.field.cell().lengths().map(|x| x / ANGSTROM),
            d.scale,
            d.repacked,
            sys.solvent,
            model_digest(&sys.field)
        ),
    );
    let (start, rows) = melted_decane(&sys, &root, DECANE_MELT_SEED, &watch, RECORD_EVERY);
    let discard = (SETTLE_DISCARD * 1e-12 / RECORD_EVERY).round() as usize;
    let verdict = settle_verdict(&rows, sys.melt.len() - 1, discard, band);
    let (lines, passed) = match verdict {
        Ok(l) => (l, true),
        Err(l) => (l, false),
    };
    for line in &lines {
        log(&root, tag, line);
    }
    assert!(
        passed,
        "the settle did not pass its gate, so no window runs: {}",
        lines.join("\n")
    );
    let dir = root.join(tag);
    std::fs::create_dir_all(&dir).expect("the windows' directory");
    let p = decane_protocol(base);
    let kt = p.kt();
    let leg = Leg::new(sys.name, &sys.decoupling, &sys.template, start, p, grid());
    log(
        &dir,
        sys.name,
        &format!("{}: protocol {p:?}, model {:016x}", sys.name, leg.model),
    );
    let mut schedule: Vec<usize> = SCHEDULE.to_vec();
    let mut records: Vec<Recorded> = schedule
        .iter()
        .map(|&g| decane_window(&leg, &dir, g, &watch))
        .collect();
    let found = refine_schedule(
        &mut schedule,
        &mut records,
        OVERLAP_THRESHOLD,
        INSERTION_ROUNDS,
        |x, y| interval_bar(x, y, kt).overlap,
        |g| decane_window(&leg, &dir, g, &watch),
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
    let dec = report(&dir, sys.name, &records, kt);
    log(
        &dir,
        sys.name,
        &format!(
            "the decanes' nematic order in each window, a diagnostic (band ≤ {band:.4})\n| \
             candidate | λ_e | λ_v | order | τ (samples) | halves |\n| --- | --- | --- | --- | --- \
             | --- |"
        ),
    );
    for r in &records {
        let (h1, h2, apart) = halves(&r.observed);
        let e = Estimate::of(&r.observed);
        log(
            &dir,
            sys.name,
            &format!(
                "| {} | {} | {} | {:.4} ± {:.4} | {:.1} | {:.4}, {:.4} ({apart:.1} σ) |",
                r.grid,
                r.lambda.electrostatics,
                r.lambda.van_der_waals,
                e.mean,
                e.error,
                e.tau,
                h1.mean,
                h2.mean
            ),
        );
    }
    let rt_ln10 = kcal(BOLTZMANN.to_si() * KELVIN) * std::f64::consts::LN_10;
    let (g_hyd, e_hyd) = (-water.bar.0, water.bar.1);
    let (g_dec, e_dec) = (-dec.bar.0, dec.bar.1);
    let log_p = |hyd: f64| -(g_dec - hyd) / rt_ln10;
    let err = |e: f64| (e_dec * e_dec + e * e).sqrt() / rt_ln10;
    let lines = [
        format!(
            "ΔG_hyd (direct, read back) {g_hyd:+.3} ± {e_hyd:.3} kcal/mol (TI {:+.3}, halves {:+.3}, \
             {:+.3}, lowest overlap {:.3}); the staged route's {:+.2} ± {:.2}",
            -water.ti.0,
            -water.halves[0],
            -water.halves[1],
            water.overlap,
            HYDRATION_STAGED.0,
            HYDRATION_STAGED.1
        ),
        format!(
            "ΔG_solv,dec {g_dec:+.3} ± {e_dec:.3} kcal/mol (TI {:+.3}, halves {:+.3}, {:+.3}, \
             lowest overlap {:.3})",
            -dec.ti.0, -dec.halves[0], -dec.halves[1], dec.overlap
        ),
        format!(
            "ΔG_dec − ΔG_hyd {:+.3} (direct), {:+.3} (staged water)",
            g_dec - g_hyd,
            g_dec - HYDRATION_STAGED.0
        ),
        format!(
            "log P(decane/water) {:.2} ± {:.2} from the two direct legs; {:.2} ± {:.2} with the \
             staged water; against n-hexadecane:water {HEXADECANE_WATER} (25 °C, log L¹⁶ 2.786 + \
             log K_aw, {:.3}) and olive oil:saline {OLIVE_OIL_SALINE} (37 °C, Gargas 1989); no \
             measured decane:water value was read; {:.0} s",
            log_p(g_hyd),
            err(e_hyd),
            log_p(HYDRATION_STAGED.0),
            err(HYDRATION_STAGED.1),
            hexadecane_water(),
            t.elapsed().as_secs_f64()
        ),
    ];
    for line in lines {
        log(&root, tag, &line);
    }
}

/// **Benzene's log P(decane/water) at W3's protocol** (release, one core; written, costed, not
/// run). The water leg read back from 4c-2's run; the decane leg melted, gated, and thirteen
/// windows of 20 + 200 ps ([`LONG`]), under `PANTOMETRY_DECANE_DIR` in `start.txt`, `melt.txt`
/// and `long/`. About 14 h: 13.9 h at the 32.5 ms a step
/// [`the_decane_leg_cost_measured`] read. `PANTOMETRY_DECANE_DIR=C:/pantometry-w4/decane
/// PANTOMETRY_OCTANOL_DIR=C:/pantometry-w4/octanol cargo test --release --locked -p
/// pantometry-forcefield --test benzene_partitioned_into_decane the_decane_partition_measured --
/// --ignored --nocapture`.
#[test]
#[ignore = "the decane leg melted and thirteen windows at W3's protocol: about 14 h with --release, not run"]
fn the_decane_partition_measured() {
    decane_partition_measured("long", LONG);
}

/// **The water leg read back**, alone (release, seconds; writes nothing): [`water_read_back`]'s
/// table, so that the production's first step can be seen to read before anything runs.
/// `PANTOMETRY_OCTANOL_DIR=C:/pantometry-w4/octanol cargo test --release --locked -p
/// pantometry-forcefield --test benzene_partitioned_into_decane the_water_leg_read_back --
/// --ignored --nocapture`.
#[test]
#[ignore = "reads 4c-2's water leg back from PANTOMETRY_OCTANOL_DIR: seconds with --release"]
fn the_water_leg_read_back() {
    let (schedule, w) = water_read_back();
    println!(
        "water/long/ windows {schedule:?}: ΔG_decouple BAR {:+.3} ± {:.3}, TI {:+.3} ± {:.3}, \
         halves {:+.3}, {:+.3}, lowest overlap {:.3}; ΔG_hyd {:+.3}",
        w.bar.0, w.bar.1, w.ti.0, w.ti.1, w.halves[0], w.halves[1], w.overlap, -w.bar.0
    );
}

/// **The cost** (release, one core, a few minutes; writes nothing): the production's box as built,
/// 200 steps at 0.5 fs at `λ_v` 0.5 with no charge, then 300 at 2 fs timed there; the couplings at
/// all 29 candidates, the order and the pressure each timed twenty times; and the leg projected
/// from them — the melt with its records, a window, thirteen and seventeen.
#[test]
#[ignore = "the production box built and stepped a few hundred times: minutes with --release"]
fn the_decane_leg_cost_measured() {
    let t = std::time::Instant::now();
    let b = benzene();
    let a = Alkane::new(D10);
    let (sys, d) = decane_leg(&b, &a);
    let watch = Watch::new(&sys, &a);
    println!(
        "benzene and {} decanes ({} removed), {} atoms, box {:.2?} Å, scaled by {:.6} and repacked \
         in {} iterations, built in {:.1} s",
        sys.solvent,
        d.removed,
        sys.built.len(),
        sys.field.cell().lengths().map(|x| x / ANGSTROM),
        d.scale,
        d.repacked,
        t.elapsed().as_secs_f64()
    );
    let potential = AtLambda {
        hamiltonian: &sys.decoupling,
        lambda: Lambda::new(0.0, 0.0, 0.5),
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
    let s = std::time::Instant::now();
    md.run(&potential, &mut at, 2.0 * FS, steps);
    let step = s.elapsed().as_secs_f64() / steps as f64;
    let states = grid();
    let time = |f: &dyn Fn()| {
        let s = std::time::Instant::now();
        for _ in 0..20 {
            f();
        }
        s.elapsed().as_secs_f64() / 20.0
    };
    let couplings = time(&|| {
        std::hint::black_box(sys.decoupling.couplings(&at, &states));
    });
    let order = time(&|| {
        std::hint::black_box(watch.order(&at));
    });
    let record = time(&|| {
        std::hint::black_box((pressure(&sys, &at), watch.trans(&at)));
    }) + order;
    let melt: f64 = sys
        .melt
        .iter()
        .map(|st| {
            let n = (st.picoseconds * 1e-12 / st.dt).round();
            n * step + (st.picoseconds * 1e-12 / RECORD_EVERY).round() * record
        })
        .sum();
    let p = decane_protocol(LONG);
    let window = p.steps_per_window() as f64 * step + p.samples as f64 * (couplings + order);
    println!(
        "a step at 2 fs {:.1} ms; the couplings at all {} candidates {:.1} ms; the order {:.3} ms; a \
         melt record (P, order, trans) {:.1} ms",
        1e3 * step,
        states.len(),
        1e3 * couplings,
        1e3 * order,
        1e3 * record
    );
    println!(
        "the melt {:.2} h; a window {:.2} h; 13 windows {:.1} h, 17 {:.1} h; the leg with its melt \
         {:.1} h at thirteen, {:.1} h at seventeen",
        melt / 3600.0,
        window / 3600.0,
        13.0 * window / 3600.0,
        17.0 * window / 3600.0,
        (melt + 13.0 * window) / 3600.0,
        (melt + 17.0 * window) / 3600.0
    );
}

// ---------------------------------------------------------------------------------------------
// Default checks

/// Benzene in a [`SMALL`] decane box, GAFF's everywhere, through classical Ewald or the mesh; and the
/// systems its end states should be, each built here from the rules written out.
struct Small {
    field: PeriodicForceField,
    decoupling: PeriodicDecoupling,
    benzene_in_the_box: PeriodicForceField,
    benzene_in_vacuum: Vacuum,
    solvent_alone: PeriodicForceField,
    at: Vec<[f64; 3]>,
    d: DecaneBox,
    alkane: Alkane,
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

fn small_decane(seed: u64, mesh: bool) -> Small {
    let b = benzene();
    let a = Alkane::new(D10);
    let d = decane_box(&b, &a, SMALL, SMALL_CUTOFF, 1e-5, seed, false);
    let l = &d.liquid;
    let cell = l.cell();
    let p = d.field.ewald().parameters();
    let bm = Model::benzene(&b.component);
    let mut benzene_in_the_box = rules(
        PeriodicForceField::new(&b.component, &b.types, cell, p).expect("benzene"),
        bm.vdw.clone(),
        bm.charges.clone(),
    );
    let parts: Vec<&Component> = (0..l.octanols).map(|_| &a.molecule.component).collect();
    let (mut vdw, mut q) = (Vec::new(), Vec::new());
    for _ in 0..l.octanols {
        vdw.extend(a.vdw());
        q.extend_from_slice(&a.transfer.charges);
    }
    let mut solvent_alone = rules(
        PeriodicForceField::new(&joined(&parts), &l.types[12..], cell, p).expect("decane alone"),
        vdw,
        q,
    );
    let mut field = d.field.clone();
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
        d,
        alkane: a,
    }
}

/// `s` after `steps` of constrained BAOAB at 0.5 fs in a 50 ps⁻¹ bath at `lambda`: the start is a
/// packing's, at 0.8 of the contact radii, whose strain C-1b's first stage takes out at that step.
fn moved(s: &Small, lambda: Lambda, steps: usize) -> Vec<[f64; 3]> {
    let mut at = s.at.clone();
    let mut md =
        s.d.liquid
            .dynamics()
            .with_bath(Bath::Langevin {
                temperature: KELVIN,
                friction: 50e12,
                seed: 17,
            })
            .thermalised(&at, KELVIN, 17);
    md.run(
        &AtLambda {
            hamiltonian: &s.decoupling,
            lambda,
        },
        &mut at,
        0.5 * FS,
        steps,
    );
    at
}

/// Benzene's model alone, as a digest: its twelve atoms' `[x, D]` and charges, the mixing and the
/// 1-4 scales — what the octanol and water legs' benzene must be to the bit.
fn benzene_digest(f: &PeriodicForceField) -> u64 {
    digest(
        f.van_der_waals()[..12]
            .iter()
            .flatten()
            .copied()
            .chain(f.charges()[..12].iter().copied())
            .chain([f64::from(u8::from(f.mixing() == Mixing::LorentzBerthelot))])
            .chain(f.one_four()),
    )
}

/// **Both models are what they say**, on the small decane box. Benzene: every carbon `ca` and
/// hydrogen `ha` by element, FreeSolv's AM1-BCC charges ([`Model::benzene`], which
/// `benzene_bound_in_tip3p.rs` holds to variant 4's to the bit) and not QEq's; and **the octanol and
/// water legs' benzene to the bit** — its `[x, D]`, charges, mixing and 1-4 scales by
/// [`benzene_digest`], and its own pairs' scales — against 4c-2's two small boxes built by
/// `tests/partition/`'s own code. The decane: every atom's type read off its element, `c3` or `hc`,
/// FreeSolv's decane file's own types; each atom's `[x, D]` GAFF's and its charge the transferred
/// one, to the bit, and **not QEq's**; each molecule neutral; Lorentz–Berthelot and GAFF's 1-4
/// scales, on benzene's 21 1-4 pairs and `9 (n − 1)` = 81 a decane.
#[test]
fn both_models_are_what_they_say() {
    let s = small_decane(0xC2A, false);
    let f = &s.field;
    let b = benzene();
    let a = &s.alkane;
    let qeq_benzene = Benzene::new(None).charges;
    let bm = Model::benzene(&b.component);
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
            "{i}"
        );
        assert_eq!(
            q[i].to_bits(),
            bm.charges[i].to_bits(),
            "benzene's atom {i}"
        );
        assert_ne!(
            q[i].to_bits(),
            qeq_benzene[i].to_bits(),
            "QEq's charge on {i}"
        );
    }
    // The decanes, molecule by molecule.
    let o = &a.molecule;
    let ideal: Vec<[f64; 3]> = o.component.atoms().iter().map(|x| x.at).collect();
    let qeq_decane = qeq::charges(&o.component, &ideal).expect("QEq").charges;
    let file = freesolv::decane();
    for x in &file.atoms {
        let want = if x.element == Element::C { "c3" } else { "hc" };
        assert_eq!(x.kind, want, "FreeSolv's decane, {}", x.name);
    }
    let l = &s.d.liquid;
    assert_eq!(l.molecules.len(), 1 + l.octanols);
    for m in &l.molecules[1..] {
        assert_eq!(m.len(), 32);
        let mut net = 0.0;
        let mut scale = 0.0;
        for (k, i) in m.clone().enumerate() {
            let kind = match o.component.atoms()[k].element {
                Element::C => "c3",
                Element::H => "hc",
                other => panic!("decane has no {other:?}"),
            };
            assert_eq!(a.kinds[k], kind);
            assert_eq!(
                vdw[i].map(f64::to_bits),
                gaff(kind).map(f64::to_bits),
                "atom {i}"
            );
            assert_eq!(q[i].to_bits(), a.transfer.charges[k].to_bits(), "atom {i}");
            assert_ne!(
                q[i].to_bits(),
                qeq_decane[k].to_bits(),
                "QEq's charge on atom {i}"
            );
            net += q[i];
            scale += q[i].abs();
        }
        assert!(net.abs() <= 64.0 * EPS * scale, "a decane carries {net} e");
    }
    assert_eq!(f.mixing(), Mixing::LorentzBerthelot);
    assert_eq!(f.one_four(), ONE_FOUR);
    assert_eq!(f.one_four_pairs().len(), 21 + 81 * l.octanols);
    // One benzene model in all three legs.
    let (ob, oo) = (benzene(), octanol());
    let (_, octanol_field, _) = octanol_box(
        &ob,
        &oo,
        partition::SMALL,
        partition::SMALL_CUTOFF,
        1e-5,
        0xC2B,
        false,
    );
    let water = with_gaff_benzene(
        Solvated::new(&b, 3, 0xC2C, 2.3 * ANGSTROM, 4.6 * ANGSTROM, 1e-5),
        &b,
    );
    let mine = benzene_digest(f);
    for (name, other) in [("octanol", &octanol_field), ("water", &water.field)] {
        assert_eq!(benzene_digest(other), mine, "benzene in {name}");
        let group: Vec<bool> = (0..other.charges().len()).map(|i| i < 12).collect();
        let dc = PeriodicDecoupling::new(other, &group).expect("benzene decouples");
        assert_eq!(dc.own_scales(), s.decoupling.own_scales(), "{name}");
    }
    assert_eq!(
        s.decoupling
            .own_scales()
            .iter()
            .filter(|x| **x == ONE_FOUR)
            .count(),
        21
    );
    println!(
        "benzene and {} decanes ({} removed), {} atoms, {} 1-4 pairs; benzene's model {mine:016x} in \
         decane, octanol and water",
        l.octanols,
        s.d.removed,
        q.len(),
        f.one_four_pairs().len()
    );
}

/// **The decane is at its density around benzene, and nothing overlaps**, in the production's box:
/// 112 sites, benzene at the centre, the decanes it overlaps removed, the box shrunk to make up —
/// `N_dec m / (V − v_B)` is 0.7255 g/cm³ within its rounding — and not so without the
/// compensation, which is C-1b's start unscaled and needs no packing. No two atoms nearer than their packing radii' sum by brute
/// force over the 27 images, benzene included; **every molecule whole**: each decane's bonds its
/// relaxed conformer's and benzene's its own to rounding, every bond under half the shortest edge,
/// and the molecules the atoms in order, benzene's twelve first.
#[test]
fn the_decane_is_at_its_density_around_benzene_and_nothing_overlaps() {
    let b = benzene();
    let a = Alkane::new(D10);
    let d = production_box(&b, &a);
    let l = &d.liquid;
    let cell = l.cell();
    let v = cell.volume();
    let m = a.molecule.mass();
    let rho = l.octanols as f64 * m / (v - BENZENE_VOLUME);
    println!(
        "{} of 112 decanes kept; C-1b's start scaled by {:.6}, repacked in {} iterations, all 112 \
         sites at {:.6} g/cm³; the decane at {:.10} g/cm³ in the volume \
         benzene leaves; box {:.3} × {:.3} × {:.3} Å",
        l.octanols,
        d.scale,
        d.repacked,
        d.lattice / 1e3,
        rho / 1e3,
        cell.lengths()[0] / ANGSTROM,
        cell.lengths()[1] / ANGSTROM,
        cell.lengths()[2] / ANGSTROM
    );
    assert!(l.octanols < 112, "benzene took somebody's place");
    assert_eq!(d.removed, 112 - l.octanols);
    assert!(
        (rho - DECANE.density).abs() <= 1e3 * EPS * DECANE.density,
        "{rho}"
    );
    assert!(d.scale < 1.0, "the box shrinks");
    let plain = decane_box(
        &b,
        &a,
        DECANE.lattice,
        9.0 * ANGSTROM,
        1e-5,
        START_SEED,
        false,
    );
    let pl = &plain.liquid;
    assert_eq!(plain.scale, 1.0);
    assert_eq!(
        plain.repacked, 0,
        "C-1b's start needs no packing with benzene in"
    );
    let plain_rho = pl.octanols as f64 * m / (pl.cell().volume() - BENZENE_VOLUME);
    println!(
        "uncompensated: {} kept, {:.5} g/cm³, {:+.2}%",
        pl.octanols,
        plain_rho / 1e3,
        100.0 * (plain_rho / DECANE.density - 1.0)
    );
    assert!(
        (plain_rho - DECANE.density).abs() > 1e-3 * DECANE.density,
        "the compensation does something"
    );

    let (w, i, j) = worst_contact(&l.at, &l.molecules, &d.radii, &cell);
    println!(
        "the worst pair over 27 images {w:.5} of its radii' sum, atoms {i} and {j}, {:.3} Å apart",
        w * (d.radii[i] + d.radii[j]) / ANGSTROM
    );
    assert!(w >= 1.0, "an overlap at {w} of the radii");
    let mut next = 0;
    for r in &l.molecules {
        assert_eq!(r.start, next, "the molecules are the atoms in order");
        next = r.end;
    }
    assert_eq!(next, l.at.len());
    assert_eq!(l.molecules[0], 0..12);
    let mut worst = 0.0f64;
    let mut longest = 0.0f64;
    let mut bond = |here: f64, own: f64| {
        worst = worst.max((here / own - 1.0).abs());
        longest = longest.max(here);
    };
    for x in b.component.bonds() {
        let [p, q] = x.atoms;
        bond(len(sub(l.at[p], l.at[q])), len(sub(b.at[p], b.at[q])));
    }
    for r in &l.molecules[1..] {
        for x in a.molecule.component.bonds() {
            let [p, q] = x.atoms;
            bond(
                len(sub(l.at[r.start + p], l.at[r.start + q])),
                len(sub(a.molecule.at[p], a.molecule.at[q])),
            );
        }
    }
    println!(
        "every bond its molecule's own to {worst:.1e}; the longest {:.3} Å, half the shortest edge \
         {:.3} Å",
        longest / ANGSTROM,
        cell.half_shortest_edge() / ANGSTROM
    );
    assert!(worst <= 64.0 * EPS);
    assert!(longest < cell.half_shortest_edge());
}

/// **The end states are the systems they say**, on the small box at a configuration constrained
/// dynamics left. Fully coupled: the whole field less benzene alone in the box, plus benzene alone in
/// vacuum with its own terms by [`Vacuum`]'s direct sums; fully decoupled, the decanes alone in their
/// box plus benzene in vacuum, energy and every force. Each within `10⁻¹²` of the parts' magnitudes:
/// the same terms summed in other orders.
#[test]
fn the_end_states_are_the_systems_they_say() {
    let s = small_decane(0xC2D, false);
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
    println!("coupled {coupled:.10e} J against {expected:.10e}");
    assert!((coupled - expected).abs() <= 1e-12 * scale, "coupled");
    let decoupled = s
        .decoupling
        .energy_and_forces(&at, Lambda::new(0.0, 0.0, 0.0), &mut f);
    let solvent = s.solvent_alone.evaluate(&at[12..]);
    let expected = solvent.energy.total + vac;
    println!("decoupled {decoupled:.10e} J against {expected:.10e}");
    assert!((decoupled - expected).abs() <= 1e-12 * scale, "decoupled");
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
                "atom {i} axis {k}: {:e} against {:e}",
                f[i][k],
                want[k]
            );
        }
    }
}

/// **`∂U/∂λ_e` and `∂U/∂λ_v` are the whole energy's derivatives**, on the small box, at four states
/// — the charges part-way, both part-way, the van der Waals part-way with no charge, and near the
/// decoupled end — each at a configuration dynamics at that state left, against central differences
/// of `energy_and_forces`, each to the tolerance it earns and that tolerance a part in 10³ of the
/// derivative.
#[test]
fn the_lambda_derivatives_are_the_whole_energys() {
    let s = small_decane(0xC2E, false);
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
                "{lambda:?}: ∂U/∂λ[{k}] = {:+.10e} J, by differences {dl:+.10e} ± {tolerance:.1e}",
                c.gradient[k]
            );
            assert!((c.gradient[k] - dl).abs() <= tolerance, "{lambda:?} {k}");
            assert!(
                tolerance < 1e-3 * c.gradient[k].abs(),
                "the check can see a part in 10³"
            );
        }
    }
}

/// **Every state of `couplings` is its own `coupling`, to the bit**, at all 29 candidates and in a
/// shuffled call; and `U(λ) − U(0)` is the coupling's energy within `128 ε` of the parts, on the
/// small box.
#[test]
fn every_state_of_couplings_is_its_own_coupling() {
    let s = small_decane(0xC2F, false);
    let at = moved(&s, Lambda::new(0.0, 0.0, 0.6), 20);
    let states = grid();
    let all = s.decoupling.couplings(&at, &states);
    let some = s
        .decoupling
        .couplings(&at, &[states[20], states[3], states[20]]);
    for (k, l) in states.iter().enumerate() {
        let one = s.decoupling.coupling(&at, *l);
        assert_eq!(one.energy.to_bits(), all[k].energy.to_bits(), "{l:?}");
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
        assert!(r <= 1.0, "{l:?}: {du:e} against {:e}", all[k].energy);
    }
    println!("U(λ) − U(0) against the coupling, the worst {worst:.3} of 128 ε");
}

/// **The decoupled reciprocal energy is `rest + λ_e · cross`**: the field's reciprocal energy with
/// only the decanes' charges, with only benzene's, and with all, against the decoupling's at three
/// states, within `64 ε` of the parts — through classical Ewald and through the mesh the production
/// uses. Not a difference, so it tells the rest from the group; and the cross term must be seen.
#[test]
fn the_decoupled_reciprocal_energy_is_the_rest_and_the_cross_terms() {
    for (name, s) in [
        ("classical", small_decane(0xC30, false)),
        ("mesh", small_decane(0xC30, true)),
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

/// Rows of a settle for the gate's closed cases: `n` records, each column a constant plus a drift
/// `slope` a record plus noise uniform in `±spread`, from `Rng::for_index(seed, record)`.
fn synthetic(n: usize, seed: u64, columns: [(f64, f64, f64); 6]) -> Vec<MeltRow> {
    (0..n)
        .map(|i| {
            let mut rng = pantometry_core::Rng::for_index(seed, i as u64);
            let mut r = [0.0; 6];
            for (c, &(level, slope, spread)) in columns.iter().enumerate() {
                r[c] = level + slope * i as f64 + spread * (2.0 * rng.unit() - 1.0);
            }
            (2, r)
        })
        .collect()
}

/// **The settle's gate on closed cases**: 1 800 records of white noise about constant `U`, `P` and
/// an order at 0.15 pass; and each is refused by the name of what failed — `U` drifting by 0.2
/// kcal/mol over the settle, `P` by 200 bar, the order rising from 0.12 to 0.18, an order flat at
/// 0.22 above a band of 0.157 and its 0.04, and a settle too short to judge. The drifts are each a few
/// tens of the halves' errors, and the trans fraction drifting is reported and not refused.
#[test]
fn the_settle_gate_on_closed_cases() {
    let band = 0.157;
    let flat = [
        (0.0, 0.01, 0.0),
        (12.6, 0.0, 0.3),
        (-30.0, 0.0, 300.0),
        (298.15, 0.0, 10.0),
        (0.15, 0.0, 0.01),
        (0.85, 0.0, 0.01),
    ];
    let n = 2000;
    let run =
        |cols: [(f64, f64, f64); 6]| settle_verdict(&synthetic(n, 0x6A7E, cols), 2, 200, band);
    let ok = run(flat).unwrap_or_else(|l| panic!("a flat settle refused:\n{}", l.join("\n")));
    println!("{}", ok.join("\n"));
    let mut trans = flat;
    trans[5].1 = 0.05 / n as f64;
    assert!(run(trans).is_ok(), "the trans fraction is not gated");
    for (what, column, slope, level, says) in [
        (
            "U drifting",
            1,
            0.2 / n as f64,
            None,
            "U/N (kcal/mol)'s halves disagree",
        ),
        (
            "P drifting",
            2,
            200.0 / n as f64,
            None,
            "P (bar)'s halves disagree",
        ),
        (
            "the order rising",
            4,
            0.06 / n as f64,
            Some(0.12),
            "the order is rising",
        ),
        ("the order high", 4, 0.0, Some(0.22), "above the band"),
    ] {
        let mut cols = flat;
        cols[column].1 = slope;
        if let Some(x) = level {
            cols[column].0 = x;
        }
        let why = run(cols).expect_err(what).join("\n");
        assert!(why.contains(says), "{what}: {why}");
        println!("{what}: {}", why.lines().last().unwrap_or(""));
    }
    let short = settle_verdict(&synthetic(205, 0x6A7E, flat), 2, 200, band).expect_err("short");
    assert!(short.join(" ").contains("too few"), "{short:?}");
}

/// A fresh directory under the system's temporary one, named by the process and `tag`.
#[cfg(not(target_family = "wasm"))]
fn fresh_dir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("pantometry-c2-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("a temporary directory");
    d
}

/// The resume test's melt: a few steps of strain out and a few settling.
#[cfg(not(target_family = "wasm"))]
const TEST_MELT: [MeltStage; 2] = [
    stage(
        "strain out",
        Lambda::COUPLED,
        0.002,
        0.5 * FS,
        50e12,
        KELVIN,
    ),
    stage("settle", Lambda::COUPLED, 0.0015, 0.5 * FS, 5e12, KELVIN),
];

/// **A killed leg resumes to the same bits**, through this file's own [`melted_decane`],
/// [`decane_window`] and `tests/partition/`'s [`Leg`] on the small box, a melt of seven steps
/// recording every step, and windows of eight 0.5 fs steps and three samples at candidates 0 and 8,
/// in temporary directories:
///
/// - **uninterrupted**: the melt and both windows;
/// - **killed and resumed**: the melt and window 0, then what a run killed during window 8 leaves —
///   its `.partial`, a header, one sample and its order — then the leg again. Window 0, `start.txt`
///   and `melt.txt` are read, their bytes unchanged; the partial is begun afresh and renamed; and
///   **every gradient, energy and order of both records, and every row of the melt's record, is the
///   uninterrupted run's to the bit**.
///
/// And each guard refuses by name: the melt asked for with another seed; a window asked for at
/// another seed, from a start one ulp away, and **from the same start under another model** — the
/// decanes' charges QEq's — which only the header's model digest can tell; a window run without its
/// order; a `start.txt` cut short, and a `melt.txt` cut short.
#[cfg(not(target_family = "wasm"))]
#[test]
fn a_killed_leg_resumes_to_the_same_bits() {
    let s = small_decane(0xC31, false);
    let sys = decane_system(&s.d, s.field.clone(), &TEST_MELT);
    let watch = Watch::new(&sys, &s.alkane);
    let every = 0.5 * FS;
    let protocol = Protocol {
        time_step: 0.5e-15,
        temperature: KELVIN,
        friction: 50e12,
        equilibration: 2,
        stride: 2,
        samples: 3,
        seed: DECANE_SEED,
    };
    let run = |dir: &std::path::Path, schedule: &[usize]| -> (Vec<Recorded>, Vec<MeltRow>) {
        let (start, rows) = melted_decane(&sys, dir, DECANE_MELT_SEED, &watch, every);
        let leg = Leg::new(
            "decane",
            &sys.decoupling,
            &sys.template,
            start,
            protocol,
            grid(),
        );
        let records = schedule
            .iter()
            .map(|&g| decane_window(&leg, dir, g, &watch))
            .collect();
        (records, rows)
    };
    let whole = fresh_dir("whole");
    let (full, full_rows) = run(&whole, &[0, 8]);
    assert_eq!(full_rows.len(), 7, "a record every step");
    let killed = fresh_dir("killed");
    run(&killed, &[0]);
    let bytes = |p: std::path::PathBuf| std::fs::read(p).expect("a file");
    let finished = bytes(killed.join("window_00.txt"));
    let start_bytes = bytes(killed.join("start.txt"));
    let record_bytes = bytes(killed.join("melt.txt"));
    let text = std::fs::read_to_string(whole.join("window_08.txt")).unwrap();
    let partial: String = text.lines().take(3).map(|l| format!("{l}\n")).collect();
    assert!(partial
        .lines()
        .nth(2)
        .unwrap_or("")
        .starts_with("# observed "));
    std::fs::write(killed.join("window_08.partial"), partial).unwrap();
    let (resumed, resumed_rows) = run(&killed, &[0, 8]);
    assert!(
        bytes(killed.join("window_00.txt")) == finished,
        "window 0 read, not run"
    );
    assert!(bytes(killed.join("start.txt")) == start_bytes);
    assert!(bytes(killed.join("melt.txt")) == record_bytes);
    assert!(!killed.join("window_08.partial").exists());
    let row_bits = |rows: &[MeltRow]| -> Vec<u64> {
        rows.iter()
            .flat_map(|r| r.1.map(f64::to_bits))
            .chain(rows.iter().map(|r| r.0 as u64))
            .collect()
    };
    assert_eq!(row_bits(&full_rows), row_bits(&resumed_rows));
    assert_eq!(full.len(), resumed.len());
    for (a, b) in full.iter().zip(&resumed) {
        assert_eq!(a.grid, b.grid);
        assert_eq!(a.gradients.len(), protocol.samples);
        assert_eq!(a.observed.len(), protocol.samples);
        assert_eq!(record_bits(a), record_bits(b), "candidate {}", a.grid);
        assert_eq!(
            a.observed.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            b.observed.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            "candidate {}'s order",
            a.grid
        );
    }
    println!(
        "resumed: {} records and {} melt rows, every number the uninterrupted run's; the order {:.4} \
         at the first sample",
        resumed.len(),
        resumed_rows.len(),
        resumed[0].observed[0]
    );

    let why = refusal(|| {
        drop(melted_decane(
            &sys,
            &killed,
            DECANE_MELT_SEED + 1,
            &watch,
            every,
        ))
    });
    assert!(why.contains("melted differently"), "{why}");
    let (start, _) = melted_decane(&sys, &killed, DECANE_MELT_SEED, &watch, every);
    let refused = |leg: &Leg, what: &str, says: &str| {
        let why = refusal(|| drop(decane_window(leg, &killed, 0, &watch)));
        assert!(why.contains(says), "{what}: {why}");
    };
    let other_seed = Protocol {
        seed: DECANE_SEED + 1,
        ..protocol
    };
    refused(
        &Leg::new(
            "decane",
            &sys.decoupling,
            &sys.template,
            start.clone(),
            other_seed,
            grid(),
        ),
        "another seed",
        "was run differently",
    );
    let mut moved = start.clone();
    moved[20][0] = f64::from_bits(moved[20][0].to_bits() + 1);
    refused(
        &Leg::new(
            "decane",
            &sys.decoupling,
            &sys.template,
            moved,
            protocol,
            grid(),
        ),
        "a start one ulp away",
        "was run differently",
    );
    // The same start under another model: the decanes' charges QEq's.
    let o = &s.alkane.molecule;
    let ideal: Vec<[f64; 3]> = o.component.atoms().iter().map(|x| x.at).collect();
    let qeq_decane = qeq::charges(&o.component, &ideal).expect("QEq").charges;
    let mut q = sys.field.charges().to_vec();
    for k in 12..q.len() {
        q[k] = qeq_decane[(k - 12) % 32];
    }
    let other = PeriodicDecoupling::new(&sys.field.clone().with_charges(q), sys.decoupling.group())
        .expect("benzene decouples");
    refused(
        &Leg::new(
            "decane",
            &other,
            &sys.template,
            start.clone(),
            protocol,
            grid(),
        ),
        "another model",
        "was run differently",
    );
    // A window whose file has no order.
    let leg = Leg::new(
        "decane",
        &sys.decoupling,
        &sys.template,
        start,
        protocol,
        grid(),
    );
    let plain = fresh_dir("plain");
    window(&leg, &plain, 0);
    std::fs::copy(plain.join("window_00.txt"), killed.join("window_04.txt")).unwrap();
    let text = std::fs::read_to_string(killed.join("window_04.txt")).unwrap();
    std::fs::write(
        killed.join("window_04.txt"),
        text.replacen(
            " candidate 0 lambda 0e0 1e0 1e0 ",
            " candidate 4 lambda 0e0 5e-1 1e0 ",
            1,
        ),
    )
    .unwrap();
    let why = refusal(|| drop(decane_window(&leg, &killed, 4, &watch)));
    assert!(why.contains("without its observation"), "{why}");
    // Files cut short.
    let cut = fresh_dir("cut");
    std::fs::write(
        cut.join("start.txt"),
        &start_bytes[..start_bytes.len() - 30],
    )
    .unwrap();
    let why = refusal(|| drop(melted_decane(&sys, &cut, DECANE_MELT_SEED, &watch, every)));
    assert!(why.contains("not whole"), "{why}");
    std::fs::write(cut.join("start.txt"), &start_bytes).unwrap();
    std::fs::write(
        cut.join("melt.txt"),
        &record_bytes[..record_bytes.len() - 30],
    )
    .unwrap();
    let why = refusal(|| drop(melted_decane(&sys, &cut, DECANE_MELT_SEED, &watch, every)));
    assert!(why.contains("melt.txt is not whole"), "{why}");
    for d in [whole, killed, cut, plain] {
        let _ = std::fs::remove_dir_all(d);
    }
}

/// **The seeds are apart**: this leg's campaign seed and every window seed it gives at all 29
/// candidates, and its melt's stage seeds, from every seed 4c-2's legs, reseed and melts give.
#[test]
fn the_seeds_are_apart_from_4c2s() {
    use pantometry_forcefield::free_energy::window_seed;
    let mut theirs: Vec<u64> = vec![MELT_SEED, 0x4C2_5EED];
    for name in ["octanol", "water"] {
        let base = seeded(LONG, name).seed;
        theirs.push(base);
        theirs.extend((0..grid().len()).map(|g| window_seed(base, g)));
    }
    theirs.extend((0..grid().len()).map(|g| window_seed(0x4C2_5EED, g)));
    theirs.extend((0..OCTANOL_MELT.len()).map(|k| window_seed(MELT_SEED, k)));
    let mut mine: Vec<u64> = vec![DECANE_SEED, DECANE_MELT_SEED];
    mine.extend((0..grid().len()).map(|g| window_seed(DECANE_SEED, g)));
    mine.extend((0..DECANE_MELT.len()).map(|k| window_seed(DECANE_MELT_SEED, k)));
    for s in &mine {
        assert!(!theirs.contains(s), "{s:#x}");
    }
    assert!((hexadecane_water() - HEXADECANE_WATER).abs() < 0.005);
    println!(
        "{} seeds of this leg apart from {} of 4c-2's; hexadecane:water {:.4}",
        mine.len(),
        theirs.len(),
        hexadecane_water()
    );
}
