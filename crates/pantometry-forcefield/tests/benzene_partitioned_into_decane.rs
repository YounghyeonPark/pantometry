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
//! # The settle failed, and is extended
//!
//! **Run 2026-10-25, the settle failed its gate**: over 1 800 records, 26.1–206.0 ps, `U/N`'s halves
//! 0.32 σ apart and `P`'s 0.29 σ, but the order's 0.0737 ± 0.0047 and 0.1314 ± 0.0069 (τ 35 and 40
//! samples), **6.96 σ apart, "rising"**; its second half was under the band's 0.1969. One reading,
//! *a hypothesis and not a finding*: the order left 600 K fully random and relaxed up to the
//! liquid's own orientational correlation, C-1b's 0.157.
//!
//! **What follows was registered before any of it was run, and is not to be changed after**
//! ([`DECANE_EXTENSION`], [`extension_verdict`]):
//!
//! - **the settle continued 300 ps** at 298.15 K, 2 fs, its own 5 ps⁻¹ bath, from the melt's end
//!   in `start.txt`. The melt kept its positions and not its velocities, so they are thermalised
//!   afresh at 298.15 K from `window_seed(`[`EXTENSION_SEED`]`, 0)`, which also seeds the bath. The
//!   same five columns every 0.1 ps, into `extension.txt`, the end into `extended.txt`, both headed
//!   ([`extension_header`]: the stage, the seed, the start's digest and the melt's whole header)
//!   and digested, and read back every time the start is used. `melt.txt` and `start.txt` are read
//!   and never written;
//! - **the gate on the extension's last 200 ps alone**, records 1 000 to 3 000, by
//!   [`settle_verdict`]'s three criteria unchanged; a record of any other length is refused. The
//!   gate over the whole extension is printed beside it and decides nothing;
//! - if it passes, **the windows start from the extension's end**, whose digest is in every
//!   window's header, so that no window from the melt's end can be read as one of these; if it
//!   fails, the leg stops and nothing is retried with another window.
//!
//! **Run 2026-10-25, the extension failed its gate too**: over its last 200 ps `U/N`'s halves
//! were 0.14 σ apart and `P`'s 0.15 σ, the order's 0.2277 ± 0.0053 and 0.2103 ± 0.0044, 2.55 σ
//! apart and so not rising — but **its second half, 0.2103, was over the band's 0.1969**. The leg
//! stopped there, as registered, and no window ran.
//!
//! # The order at four times the size: registered, then run
//!
//! Is about 0.21 a finite-size effect — a correlated isotropic liquid's order, which falls as
//! `1/√N` — or the model's own tendency to order? **Registered before any of it was run, and not to
//! be changed after** ([`the_decane_order_at_four_times_the_size_measured`]):
//!
//! - **the box**: pure decane, **444 molecules**, four times the partition box's 111, by C-1b's
//!   start ([`hexadecane::isotropic_sites`]): C-1b's 2 × 7 × 8 lattice replicated 2 × 2 × 1, one site
//!   left empty in each replica ([`four_times_empty`]), every site shrunk by `(444/448)^(1/3)` so
//!   that the 444 are at 0.7255 g/cm³, C-1b's conformers, rotations, packer and radii, seed
//!   [`FOUR_TIMES_START_SEED`]; 14 208 atoms, 61.4 × 71.7 × 32.9 Å;
//! - **the run** ([`FOUR_TIMES_MELT`]): 1 ps of strain out, 5 ps at 600 K, then **300 ps at
//!   298.15 K**, `U`, `P`, `T`, the order and the trans fraction every 0.1 ps, seed
//!   [`FOUR_TIMES_SEED`];
//! - **the reading** ([`finite_size_verdict`]): S₄₄₄ over the last 200 ps with its τ-corrected
//!   error, and the halves of `U` and `P`;
//! - **the rule** ([`finite_size_outcome`]), S₁₁₁ = 0.2103 ([`S_111`]): **finite-size (isotropic
//!   liquid)** if S₄₄₄ ≤ 0.137 = 1.3 × S₁₁₁/√4; **genuine order** if S₄₄₄ ≥ 0.158 = 0.75 × S₁₁₁;
//!   **inconclusive** between; and **not settled**, whatever S₄₄₄, unless the halves of `U` and of
//!   `P` each agree within 3 of their combined errors;
//! - **if finite-size**, the partition leg's settle gate becomes [`finite_size_gate`]: the order's
//!   halves within 3 σ as before, `U` and `P` as before, and the order's second half at most
//!   [`FINITE_SIZE_CEILING`] = 0.30, which still refuses a nematic (the transition's 0.43, C-1's
//!   0.57). It is to be applied to the 111-decane extension's existing record without rerunning it.
//!
//! **Run 2026-10-26, the verdict is finite-size**: S₄₄₄ = 0.0668 ± 0.0062 over the last 200 ps,
//! under half the finite-size threshold, with `U`'s halves 2.12 σ apart and `P`'s 2.65 σ. The
//! order's own halves rose from 0.0565 to 0.0770, 3.58 σ, which the rule does not read; both are
//! inside the 444 band. **So [`partition_gate`] is [`finite_size_gate`]**, and the gate as first
//! registered is printed beside it and decides nothing.
//!
//! An uncorrelated random set of 444 directions shows 0.0386 ± 0.0140 (of 111, 0.0765 ± 0.0268),
//! and is reported beside the reading. Written under `PANTOMETRY_DECANE_FINITE_SIZE_DIR`, which must
//! be named: `built.txt`, `record.txt`, `end.txt`, headed ([`four_times_header`]) and digested,
//! and `run.log` and `finite_size.log`.
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
//! Written under `PANTOMETRY_DECANE_DIR` (`start.txt`, `melt.txt`, `melt.log`, `extended.txt`,
//! `extension.txt`, `extend.log`, `long/`); nothing is written under `PANTOMETRY_OCTANOL_DIR`.
//!
//! # Checks
//!
//! Default, unoptimised, on a small decane box of twelve sites at `r_c` 6.5 Å unless named: both
//! models are what they say, and benzene's is the octanol and water legs' to the bit; the decane is
//! at its density with benzene in, and nothing overlaps, in the production's box; the end states are
//! the systems they say; the λ derivatives are the whole energy's; every state of `couplings` is its
//! own `coupling`; the decoupled reciprocal energy is the rest's plus `λ_e` times the cross term; the
//! settle's gate and the extension's on closed cases; the 444-decane start at its density, whole and
//! isotropic, and deterministic; the finite-size rule's four outcomes at their edges, and the
//! finite-size gate on closed cases; the extension starts from the melt's end and
//! the windows from the extension's; and a killed leg — in its extension or its windows — resumes to
//! the same bits and refuses another model, seed or start.

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
    conformations, directions, isotropic, isotropic_order, isotropic_sites, nematic_order, pack,
    packing_radius, with_gaff, worst_contact, Alkane, D10, DECANE, START_TRANS,
};
use hydration::{Benzene, Solvated};
use nonbonded::Vacuum;
use octanol::{joined, pressure as liquid_pressure, production, Liquid, ANGSTROM, FS, KELVIN};
use pantometry_forcefield::free_energy::refine_schedule;
use pantometry_forcefield::{
    qeq, Alchemical, AtLambda, Bath, Component, Element, Estimate, EwaldParameters, Lambda, Mixing,
    MolecularDynamics, PeriodicBox, PeriodicDecoupling, PeriodicForceField, PmeParameters,
    Potential, Protocol,
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

/// **The settle continued, pre-registered after the first settle failed and before any of this was
/// run**: 300 ps more at 298.15 K, 2 fs, in the settle's own 5 ps⁻¹ bath, from the melt's end.
const DECANE_EXTENSION: [MeltStage; 1] = [stage(
    "extend",
    Lambda::COUPLED,
    300.0,
    2.0 * FS,
    5e12,
    KELVIN,
)];

/// The extension's seed: its bath and its fresh velocities are `window_seed(EXTENSION_SEED, 0)`.
const EXTENSION_SEED: u64 = 0xC2_E02;

/// **What of the extension the gate reads, picoseconds: its last 200**, records 1 000 to 3 000 at
/// 0.1 ps. Pre-registered with the extension.
const EXTENSION_GATED: f64 = 200.0;

/// How far apart two halves may be, in their combined errors, and still agree: see the module
/// documentation.
const HALVES_AGREE: f64 = 3.0;

/// **How far above the isotropic band the order may sit**: see the module documentation.
const ORDER_ALLOWANCE: f64 = 0.04;

// The finite-size test, registered after the extension's gate failed and before any of it ran:
// see the module documentation.

/// **The finite-size box's lattice**: C-1b's [`DECANE`] lattice, 2 × 7 × 8, replicated 2 × 2 × 1.
const FOUR_TIMES: [usize; 3] = [4, 14, 8];

/// How many decanes it holds: four times the partition box's 111.
const FOUR_TIMES_MOLECULES: usize = 444;

/// The 111 the order is compared with: the partition box's decanes.
const ONCE_MOLECULES: usize = 111;

/// **The four sites of [`FOUR_TIMES`] left empty**, one in each 2 × 7 × 8 replica: its site
/// `(1, 3, 4)` in the replica's own indices, one of the four nearest the replica's centre, as the
/// partition box's benzene at its centre took one decane's place. 448 sites less four is 444.
fn four_times_empty() -> Vec<usize> {
    let [_, nj, nk] = FOUR_TIMES;
    let [ri, rj, _] = DECANE.lattice;
    let mut empty = Vec::new();
    for a in 0..2 {
        for b in 0..2 {
            let (i, j, k) = (a * ri + 1, b * rj + 3, 4);
            empty.push((i * nj + j) * nk + k);
        }
    }
    empty
}

/// The finite-size box's start seed.
const FOUR_TIMES_START_SEED: u64 = 0xC1B4;

/// The finite-size run's seed: stage `k` bathed with `window_seed(FOUR_TIMES_SEED, k)`.
const FOUR_TIMES_SEED: u64 = 0xC2_F5;

/// **The finite-size run**, registered: 1 ps of strain out, 5 ps at 600 K, then **300 ps at
/// 298.15 K**: [`DECANE_MELT`]'s first two stages, then one stage at 298.15 K as long as the
/// 111-decane extension, in the same 5 ps⁻¹ bath at 2 fs.
const FOUR_TIMES_MELT: [MeltStage; 3] = [
    stage("strain out", Lambda::COUPLED, 1.0, 0.5 * FS, 50e12, KELVIN),
    stage("hot", Lambda::COUPLED, 5.0, 2.0 * FS, 5e12, 600.0),
    stage("settle", Lambda::COUPLED, 300.0, 2.0 * FS, 5e12, KELVIN),
];

/// What of the 300 ps at 298.15 K is read, picoseconds: the last 200, as the extension's gate.
const FOUR_TIMES_READ: f64 = 200.0;

/// **S₁₁₁, the order of the 111 decanes**: the extension's last 100 ps, its gated window's second
/// half, 0.2103 as `the_decane_partition_measured` printed it on 2026-10-25 (0.210 to three places).
/// The registered thresholds are 0.137 and 0.158: this value gives both to three places, where 0.210
/// itself would sit on the rounding edge of each (0.1365 and 0.1575).
const S_111: f64 = 0.2103;

/// **Finite-size**: S₄₄₄ at most this times `S₁₁₁/√(444/111)`, which an isotropic liquid's order,
/// falling as `1/√N`, would give with 30% to spare.
const FINITE_SIZE_FACTOR: f64 = 1.3;

/// **Genuine order**: S₄₄₄ at least this fraction of `S₁₁₁`, which an order that is the model's and
/// not the box's would keep.
const ORDERED_FRACTION: f64 = 0.75;

/// **The partition leg's order ceiling if the verdict is finite-size** ([`finite_size_gate`]).
const FINITE_SIZE_CEILING: f64 = 0.30;

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
    let (at, rows) = integrate(sys, sys.melt, &sys.built, seed, watch, every, dir, "melt");
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

/// **`stages` run from `from`**, each from the last: stage `k` bathed with `window_seed(seed, k)`
/// and the velocities thermalised once, before the first, from `window_seed(seed, 0)`; every
/// `every` seconds `U` per decane, `P`, `T`, the order and the trans fraction recorded, the time
/// counted from zero; each stage's cost and the record in 2.5 ps blocks logged to `name.log`. The
/// end positions and the record. The melt's and the extension's one loop.
#[allow(clippy::too_many_arguments)]
fn integrate(
    sys: &System,
    stages: &[MeltStage],
    from: &[[f64; 3]],
    seed: u64,
    watch: &Watch,
    every: f64,
    dir: &std::path::Path,
    name: &str,
) -> (Vec<[f64; 3]>, Vec<MeltRow>) {
    integrate_with(
        sys.name,
        &sys.template,
        sys.solvent,
        |st: &MeltStage| AtLambda {
            hamiltonian: &sys.decoupling,
            lambda: st.lambda,
        },
        |at: &[[f64; 3]]| pressure(sys, at),
        stages,
        from,
        seed,
        watch,
        every,
        dir,
        name,
    )
}

/// **[`integrate`]'s loop on any potential**: `label`'s `stages` from `from` under `template`'s
/// dynamics, each stage's potential `potential_of` it, `U` per `molecules`, and `pressure_of` the
/// pressure in pascals. The decane leg's melt and extension come here through [`integrate`], and
/// the finite-size box of pure decane ([`four_times_measured`]) directly, so both are one loop.
#[allow(clippy::too_many_arguments)]
fn integrate_with<P: Potential>(
    label: &str,
    template: &MolecularDynamics,
    molecules: usize,
    potential_of: impl Fn(&MeltStage) -> P,
    pressure_of: impl Fn(&[[f64; 3]]) -> f64,
    stages: &[MeltStage],
    from: &[[f64; 3]],
    seed: u64,
    watch: &Watch,
    every: f64,
    dir: &std::path::Path,
    name: &str,
) -> (Vec<[f64; 3]>, Vec<MeltRow>) {
    let mut at = from.to_vec();
    let mut md = template.clone();
    let n = molecules as f64;
    let mut rows: Vec<MeltRow> = Vec::new();
    let mut time = 0.0;
    for (k, st) in stages.iter().enumerate() {
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
        let potential = potential_of(st);
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
                        pressure_of(&at) / 1e5,
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
            name,
            &format!(
                "{} stage {k} ({}): {} ps at {} fs, {} K, γ {} /ps: {secs:.0} s, {:.1} ms a step",
                label,
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
        name,
        &format!(
            "{}: the {name} in blocks of 2.5 ps\n| stage | t (ps) | U/N (kcal/mol) | P (bar) | T (K) | order | trans |\n| --- | --- | --- | --- | --- | --- | --- |",
            label
        ),
    );
    let per_block = ((2.5e-12 / every).round() as usize).max(1);
    for block in rows.chunks(per_block) {
        let m = |c: usize| block.iter().map(|r| r.1[c]).sum::<f64>() / block.len() as f64;
        log(
            dir,
            name,
            &format!(
                "| {} | {:.1}–{:.1} | {:.4} | {:.0} | {:.1} | {:.4} | {:.4} |",
                stages[block[block.len() - 1].0].name,
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
    (at, rows)
}

/// The stages as a header gives them: each one's `λ_e/λ_v/ps/dt/γ/T/name`, as [`melt_header`]'s.
fn stages_text(stages: &[MeltStage]) -> String {
    let st: Vec<String> = stages
        .iter()
        .map(|s| {
            format!(
                "{}/{}/{}/{}/{}/{}/{}",
                exact(s.lambda.electrostatics),
                exact(s.lambda.van_der_waals),
                exact(s.picoseconds),
                exact(s.dt),
                exact(s.friction),
                exact(s.temperature),
                s.name.replace(' ', "_")
            )
        })
        .collect();
    st.join(",")
}

/// **The extension's header**: its stages, its seed, how often it records, the digest of the
/// positions it starts from, and the whole header of the melt that left them. An extension of
/// another melt, from other positions, at another seed or for another length is another header.
fn extension_header(
    sys: &System,
    melt_header: &str,
    from: &[[f64; 3]],
    stages: &[MeltStage],
    seed: u64,
    every: f64,
) -> String {
    format!(
        "# extended {} stages {} seed {seed:#x} record {} from {:016x} after {melt_header}",
        sys.name,
        stages_text(stages),
        exact(every),
        positions_digest(from)
    )
}

/// **The settle continued**: `stages` from `from` — the melt's end, read back from `start.txt` —
/// with velocities thermalised afresh from `window_seed(seed, 0)`, because the melt kept its
/// positions and not its velocities. Read back from `extended.txt` and `extension.txt` in `dir`
/// when a run wrote them under the same header ([`extension_header`]), a file extended differently
/// or not whole refused by name; otherwise run now under a lock `extend`, logged to `extend.log`,
/// and written atomically, the record first. Nothing the melt wrote is touched.
#[allow(clippy::too_many_arguments)]
fn extended_decane(
    sys: &System,
    dir: &std::path::Path,
    melt_header: &str,
    from: &[[f64; 3]],
    stages: &[MeltStage],
    seed: u64,
    watch: &Watch,
    every: f64,
) -> (Vec<[f64; 3]>, Vec<MeltRow>) {
    let start_path = dir.join("extended.txt");
    let record_path = dir.join("extension.txt");
    let want = extension_header(sys, melt_header, from, stages, seed, every);
    if let Ok(text) = std::fs::read_to_string(&start_path) {
        let (head, at) = parse_positions(&text).unwrap_or_else(|why| {
            panic!(
                "{} is not whole ({why}): remove it and extension.txt to extend again",
                start_path.display()
            )
        });
        assert_eq!(
            head,
            want,
            "{} was extended differently: remove it, extension.txt and every window run from it, \
             to extend again",
            start_path.display()
        );
        assert_eq!(at.len(), sys.built.len());
        let rows = std::fs::read_to_string(&record_path)
            .map_err(|e| e.to_string())
            .and_then(|t| parse_record(&t, &want))
            .unwrap_or_else(|why| {
                panic!(
                    "{} is not whole ({why}): remove extended.txt to extend again",
                    record_path.display()
                )
            });
        return (at, rows);
    }
    let _lock = Lock::take(dir, "extend");
    let t = std::time::Instant::now();
    let (at, rows) = integrate(sys, stages, from, seed, watch, every, dir, "extend");
    let partial = record_path.with_extension("partial");
    std::fs::write(&partial, record_text(&want, &rows)).expect("the record");
    std::fs::rename(&partial, &record_path).expect("rename");
    write_positions(&start_path, &want, &at);
    log(
        dir,
        "extend",
        &format!(
            "{} extended in {:.0} s",
            sys.name,
            t.elapsed().as_secs_f64()
        ),
    );
    (at, rows)
}

/// Where the windows start from, and the records that put them there.
struct Started {
    /// The melt's end, `start.txt`: the failed settle's last positions.
    melted: Vec<[f64; 3]>,
    /// The melt's record, `melt.txt`.
    melt: Vec<MeltRow>,
    /// The extension's record, `extension.txt`.
    extension: Vec<MeltRow>,
    /// **The windows' start**: the extension's end, `extended.txt`.
    start: Vec<[f64; 3]>,
}

/// **The leg's start**: the melt ([`melted_decane`]) and then its extension ([`extended_decane`])
/// from the melt's end, under `dir`; the windows start from the extension's end. The production and
/// the resume check both come here.
#[allow(clippy::too_many_arguments)]
fn decane_start(
    sys: &System,
    dir: &std::path::Path,
    melt_seed: u64,
    extension: &[MeltStage],
    extension_seed: u64,
    watch: &Watch,
    every: f64,
) -> Started {
    let (melted, melt) = melted_decane(sys, dir, melt_seed, watch, every);
    let head = decane_melt_header(sys, melt_seed, every);
    let (start, extension) = extended_decane(
        sys,
        dir,
        &head,
        &melted,
        extension,
        extension_seed,
        watch,
        every,
    );
    Started {
        melted,
        melt,
        extension,
        start,
    }
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
    settle_verdict_under(rows, stage, discard, OrderCeiling::Band(band))
}

/// What the order's second half may not exceed in the settle's gate.
#[derive(Clone, Copy, Debug)]
enum OrderCeiling {
    /// The isotropic band plus [`ORDER_ALLOWANCE`]: the gate as first registered.
    Band(f64),
    /// A fixed ceiling: [`FINITE_SIZE_CEILING`], registered for the case that the 444-decane box
    /// finds 0.21 a finite-size effect ([`finite_size_gate`]).
    Fixed(f64),
}

impl OrderCeiling {
    fn value(self) -> f64 {
        match self {
            OrderCeiling::Band(band) => band + ORDER_ALLOWANCE,
            OrderCeiling::Fixed(c) => c,
        }
    }

    /// How the gate's first line names it.
    fn named(self) -> String {
        match self {
            OrderCeiling::Band(band) => format!(
                "the band {band:.4} + {ORDER_ALLOWANCE} = {:.4}",
                band + ORDER_ALLOWANCE
            ),
            OrderCeiling::Fixed(c) => format!("the finite-size ceiling {c:.2}"),
        }
    }

    /// What a failure says.
    fn above(self, second: f64) -> String {
        match self {
            OrderCeiling::Band(band) => format!(
                "the order's second half {second:.4} is above the band and its allowance, {:.4}",
                band + ORDER_ALLOWANCE
            ),
            OrderCeiling::Fixed(c) => format!(
                "the order's second half {second:.4} is above the finite-size ceiling, {c:.2}"
            ),
        }
    }
}

/// [`settle_verdict`] with the order's second half held under `ceiling`; every other criterion,
/// and every line under [`OrderCeiling::Band`], as [`settle_verdict`]'s.
fn settle_verdict_under(
    rows: &[MeltRow],
    stage: usize,
    discard: usize,
    ceiling: OrderCeiling,
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
         combined errors, the order's second half at most {}",
        kept.len(),
        kept[0][0],
        kept[kept.len() - 1][0],
        ceiling.named()
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
        let under = b.mean <= ceiling.value();
        if c == 4 && !under {
            failed.push(ceiling.above(b.mean));
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

/// **The extension's pre-registered gate**: the record must be [`DECANE_EXTENSION`]'s whole length
/// at `every` — refused by name when it is shorter, longer or not one stage — and then
/// [`settle_verdict`]'s criteria, unchanged, over its last [`EXTENSION_GATED`] ps alone.
fn extension_verdict(rows: &[MeltRow], every: f64, band: f64) -> Result<Vec<String>, Vec<String>> {
    extension_verdict_under(rows, every, OrderCeiling::Band(band))
}

/// **The partition leg's settle gate if the 444-decane box finds 0.21 a finite-size effect**,
/// registered with that test and before any of it ran: [`extension_verdict`]'s record and window
/// unchanged, the `U`, `P` and order halves within [`HALVES_AGREE`] unchanged, and the order's second
/// half at most [`FINITE_SIZE_CEILING`] = 0.30 in place of the band plus 0.04. 0.30 still refuses a
/// nematic: Maier–Saupe's order at the transition is 0.43, and C-1's aligned decane was 0.57.
///
/// It was to be applied to the 111-decane extension's existing record, without rerunning, only if
/// [`finite_size_verdict`] said finite-size. It did, on 2026-10-26, so it is
/// [`partition_gate`].
fn finite_size_gate(rows: &[MeltRow], every: f64) -> Result<Vec<String>, Vec<String>> {
    extension_verdict_under(rows, every, OrderCeiling::Fixed(FINITE_SIZE_CEILING))
}

/// **The gate [`decane_partition_measured`] holds the extension to before any window runs**:
/// [`finite_size_gate`], since the 444-decane run's registered verdict was finite-size (S₄₄₄ 0.0668
/// ± 0.0062, 2026-10-26). [`extension_verdict`], the gate as first registered, is still printed
/// beside it and decides nothing.
fn partition_gate(rows: &[MeltRow], every: f64) -> Result<Vec<String>, Vec<String>> {
    finite_size_gate(rows, every)
}

/// [`extension_verdict`] with the order held under `ceiling`.
fn extension_verdict_under(
    rows: &[MeltRow],
    every: f64,
    ceiling: OrderCeiling,
) -> Result<Vec<String>, Vec<String>> {
    let st = &DECANE_EXTENSION[0];
    let whole = (st.picoseconds * 1e-12 / every).round() as usize;
    let gated = (EXTENSION_GATED * 1e-12 / every).round() as usize;
    if rows.len() != whole || rows.iter().any(|r| r.0 != 0) {
        return Err(vec![format!(
            "the extension's record has {} records where its {} ps at {} ps make {whole}, one \
             stage: not the record the gate was registered for",
            rows.len(),
            st.picoseconds,
            every * 1e12
        )]);
    }
    settle_verdict_under(rows, 0, whole - gated, ceiling)
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
    let started = decane_start(
        &sys,
        &root,
        DECANE_MELT_SEED,
        &DECANE_EXTENSION,
        EXTENSION_SEED,
        &watch,
        RECORD_EVERY,
    );
    let discard = (SETTLE_DISCARD * 1e-12 / RECORD_EVERY).round() as usize;
    let lines = match settle_verdict(&started.melt, sys.melt.len() - 1, discard, band) {
        Ok(l) | Err(l) => l,
    };
    log(
        &root,
        tag,
        "the first settle's gate, read back from melt.txt and reported: it failed, and the \
         extension's gate decides",
    );
    for line in &lines {
        log(&root, tag, line);
    }
    let lines = match settle_verdict(&started.extension, 0, 0, band) {
        Ok(l) | Err(l) => l,
    };
    log(
        &root,
        tag,
        "the extension's whole 300 ps under the same criteria, for information and not gated",
    );
    for line in &lines {
        log(&root, tag, line);
    }
    let lines = match extension_verdict(&started.extension, RECORD_EVERY, band) {
        Ok(l) | Err(l) => l,
    };
    log(
        &root,
        tag,
        "the extension's gate as first registered, band plus 0.04, read back from extension.txt: \
         it failed on 2026-10-25, and is reported and decides nothing",
    );
    for line in &lines {
        log(&root, tag, line);
    }
    let (lines, passed) = match partition_gate(&started.extension, RECORD_EVERY) {
        Ok(l) => (l, true),
        Err(l) => (l, false),
    };
    log(
        &root,
        tag,
        "the extension's gate since the 444-decane verdict, finite-size: the same last 200 ps, \
         the order at most 0.30",
    );
    for line in &lines {
        log(&root, tag, line);
    }
    assert!(
        passed,
        "the extension did not pass the finite-size gate, so no window runs: {}",
        lines.join("\n")
    );
    let dir = root.join(tag);
    std::fs::create_dir_all(&dir).expect("the windows' directory");
    let p = decane_protocol(base);
    let kt = p.kt();
    let leg = Leg::new(
        sys.name,
        &sys.decoupling,
        &sys.template,
        started.start,
        p,
        grid(),
    );
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

/// **Benzene's log P(decane/water) at W3's protocol** (release, one core). The water leg read
/// back from 4c-2's run; the decane leg melted (`start.txt`, `melt.txt`: run 2026-10-25, and its
/// settle failed its gate), the settle extended by 300 ps (`extended.txt`, `extension.txt`) and
/// gated on the extension's last 200 ps by [`partition_gate`], and thirteen windows of 20 + 200 ps ([`LONG`]) from the
/// extension's end, in `long/`, all under `PANTOMETRY_DECANE_DIR`. About 15.3 h: the extension
/// 1.4 h at the 33.7 ms a step the settle took, and the windows 13.9 h at the 32.5 ms a step
/// [`the_decane_leg_cost_measured`] read. `PANTOMETRY_DECANE_DIR=C:/pantometry-w4/decane
/// PANTOMETRY_OCTANOL_DIR=C:/pantometry-w4/octanol cargo test --release --locked -p
/// pantometry-forcefield --test benzene_partitioned_into_decane the_decane_partition_measured --
/// --ignored --nocapture`.
#[test]
#[ignore = "the decane leg melted, extended and thirteen windows at W3's protocol: about 15 h with --release"]
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
    let cost = |st: &MeltStage| {
        let n = (st.picoseconds * 1e-12 / st.dt).round();
        n * step + (st.picoseconds * 1e-12 / RECORD_EVERY).round() * record
    };
    let melt: f64 = sys.melt.iter().map(cost).sum();
    let extension: f64 = DECANE_EXTENSION.iter().map(cost).sum();
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
        "the melt {:.2} h; its extension {:.2} h; a window {:.2} h; 13 windows {:.1} h, 17 {:.1} h; \
         the leg with its melt {:.1} h at thirteen, {:.1} h at seventeen",
        melt / 3600.0,
        extension / 3600.0,
        window / 3600.0,
        13.0 * window / 3600.0,
        17.0 * window / 3600.0,
        (melt + 13.0 * window) / 3600.0,
        (melt + 17.0 * window) / 3600.0
    );
}

// ---------------------------------------------------------------------------------------------
// The order at four times the size

/// What the finite-size reading says, by the registered rule ([`finite_size_outcome`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FiniteSize {
    /// S₄₄₄ at most the finite-size threshold, `U` and `P` stationary.
    Isotropic,
    /// S₄₄₄ at least the ordered threshold, `U` and `P` stationary.
    Ordered,
    /// Between the two, `U` and `P` stationary.
    Inconclusive,
    /// `U`'s or `P`'s halves further apart than [`HALVES_AGREE`]: no verdict.
    NotSettled,
}

impl FiniteSize {
    fn named(self) -> &'static str {
        match self {
            FiniteSize::Isotropic => "finite-size (isotropic liquid)",
            FiniteSize::Ordered => "genuine order",
            FiniteSize::Inconclusive => "inconclusive",
            FiniteSize::NotSettled => "not settled",
        }
    }
}

/// **The registered thresholds**, finite-size and ordered: `1.3 × S₁₁₁/√(444/111)` and
/// `0.75 × S₁₁₁`, 0.1367 and 0.1577 — 0.137 and 0.158 to three places, as registered.
fn finite_size_thresholds() -> (f64, f64) {
    let ratio = FOUR_TIMES_MOLECULES as f64 / ONCE_MOLECULES as f64;
    (
        FINITE_SIZE_FACTOR * S_111 / ratio.sqrt(),
        ORDERED_FRACTION * S_111,
    )
}

/// **The registered rule**: *not settled* unless the halves of `U` and of `P` are each within
/// [`HALVES_AGREE`] of their combined errors (a distance that is not a number is not within it), or
/// if the order is not a number; then *finite-size* at or under the first of
/// [`finite_size_thresholds`], *genuine order* at or over the second, and *inconclusive* between.
fn finite_size_outcome(order: f64, u_apart: f64, p_apart: f64) -> FiniteSize {
    let (isotropic, ordered) = finite_size_thresholds();
    if !(u_apart <= HALVES_AGREE && p_apart <= HALVES_AGREE) || !order.is_finite() {
        FiniteSize::NotSettled
    } else if order <= isotropic {
        FiniteSize::Isotropic
    } else if order >= ordered {
        FiniteSize::Ordered
    } else {
        FiniteSize::Inconclusive
    }
}

/// **The finite-size reading**: the record must hold [`FOUR_TIMES_MELT`]'s 300 ps at 298.15 K whole
/// at `every`, refused by name otherwise; over its last [`FOUR_TIMES_READ`] ps, S₄₄₄'s mean with its
/// τ-corrected error and its halves, and the halves of `U`, `P` and the trans fraction; then
/// [`finite_size_outcome`] on S₄₄₄'s mean and the `U` and `P` halves. The outcome, and the lines that
/// say why, the last of them the verdict.
fn finite_size_verdict(rows: &[MeltRow], every: f64) -> Result<(FiniteSize, Vec<String>), String> {
    let last = FOUR_TIMES_MELT.len() - 1;
    let st = &FOUR_TIMES_MELT[last];
    let whole = (st.picoseconds * 1e-12 / every).round() as usize;
    let read = (FOUR_TIMES_READ * 1e-12 / every).round() as usize;
    let settle: Vec<&[f64; 6]> = rows.iter().filter(|r| r.0 == last).map(|r| &r.1).collect();
    if settle.len() != whole {
        return Err(format!(
            "the record's stage {last} has {} records where its {} ps at {} ps make {whole}: not \
             the record the rule was registered for",
            settle.len(),
            st.picoseconds,
            every * 1e12
        ));
    }
    let kept = &settle[whole - read..];
    let column = |c: usize| -> Vec<f64> { kept.iter().map(|r| r[c]).collect() };
    let order = Estimate::of(&column(4));
    let (o1, o2, o_apart) = halves(&column(4));
    let (u1, u2, u_apart) = halves(&column(1));
    let (p1, p2, p_apart) = halves(&column(2));
    let (t1, t2, t_apart) = halves(&column(5));
    let outcome = finite_size_outcome(order.mean, u_apart, p_apart);
    let (isotropic, ordered) = finite_size_thresholds();
    let half = |what: &str, a: &Estimate, b: &Estimate, apart: f64, gated: &str| {
        format!(
            "  {what}: halves {:.4} ± {:.4} (τ {:.1}) and {:.4} ± {:.4} (τ {:.1}), {apart:.2} σ \
             apart{gated}",
            a.mean, a.error, a.tau, b.mean, b.error, b.tau
        )
    };
    let settled = if outcome == FiniteSize::NotSettled {
        "not both within"
    } else {
        "both within"
    };
    let lines = vec![
        format!(
            "the rule, registered before the run, over the last {FOUR_TIMES_READ} ps at {} K ({} \
             records, {:.1}–{:.1} ps): finite-size (isotropic liquid) if S₄₄₄ ≤ {isotropic:.3} = \
             {FINITE_SIZE_FACTOR} × S₁₁₁/√4; genuine order if S₄₄₄ ≥ {ordered:.3} = \
             {ORDERED_FRACTION} × S₁₁₁; inconclusive between; S₁₁₁ = {S_111} from the 111-decane \
             extension; not settled unless the halves of U and of P each agree within \
             {HALVES_AGREE} of their combined errors",
            st.temperature,
            kept.len(),
            kept[0][0],
            kept[kept.len() - 1][0]
        ),
        format!(
            "  S₄₄₄, the nematic order of the 444 end-to-end vectors: {:.4} ± {:.4} (τ {:.1} \
             samples)",
            order.mean, order.error, order.tau
        ),
        half("order", &o1, &o2, o_apart, "; reported, not in the rule"),
        half("U/N (kcal/mol)", &u1, &u2, u_apart, ""),
        half("P (bar)", &p1, &p2, p_apart, ""),
        half("trans fraction", &t1, &t2, t_apart, "; reported, not in the rule"),
        format!(
            "the verdict: {}: S₄₄₄ {:.4} ± {:.4} against finite-size ≤ {isotropic:.3} and genuine \
             order ≥ {ordered:.3}; U's halves {u_apart:.2} σ and P's {p_apart:.2} σ apart, {settled} \
             {HALVES_AGREE}",
            outcome.named(),
            order.mean,
            order.error
        ),
    ];
    Ok((outcome, lines))
}

/// A field as [`integrate_with`]'s potential, for a box with nothing decoupled.
struct Plain<'a>(&'a PeriodicForceField);

impl Potential for Plain<'_> {
    fn energy_and_forces(&self, at: &[[f64; 3]], forces: &mut [[f64; 3]]) -> f64 {
        Potential::energy_and_forces(self.0, at, forces)
    }
}

/// **The finite-size box**: [`isotropic_sites`] on [`FOUR_TIMES`] less [`four_times_empty`], at
/// [`DECANE`]'s density — C-1b's start procedure, its conformers at [`START_TRANS`], its uniform
/// rotations, its packer and radii and budget — from [`FOUR_TIMES_START_SEED`]; `r_c` 9 Å, Ewald at
/// δ = 10⁻⁵ before the production's mesh. Pure decane: no benzene.
fn four_times_start(a: &Alkane) -> hexadecane::Start {
    isotropic_sites(
        a,
        FOUR_TIMES,
        &four_times_empty(),
        DECANE.density,
        FOUR_TIMES_START_SEED,
        START_TRANS,
        9.0 * ANGSTROM,
        1e-5,
        PACKING_ITERATIONS,
    )
}

/// The finite-size box's field: the decanes' GAFF model ([`with_gaff`]), with the production's mesh
/// and neighbour list.
fn four_times_field(a: &Alkane, l: &Liquid) -> PeriodicForceField {
    production(with_gaff(l.field.clone(), a), l)
}

/// What the finite-size box is judged by: every molecule is a decane.
fn four_times_watch(l: &Liquid, a: &Alkane) -> Watch {
    Watch {
        solvent: l.molecules.clone(),
        ends: a.ends(),
        carbons: a.carbons.clone(),
    }
}

/// **The finite-size run's header**: the lattice, the sites left empty, the start's seed, trans
/// probability and density, the stages, the run's seed, how often it records, the atom count and
/// digests of the box as built and of the model.
fn four_times_header(l: &Liquid, field: &PeriodicForceField, every: f64) -> String {
    format!(
        "# finite-size decane molecules {} sites {:?} empty {:?} start-seed {FOUR_TIMES_START_SEED:#x} \
         trans {} density {} stages {} seed {FOUR_TIMES_SEED:#x} record {} atoms {} built {:016x} \
         model {:016x}",
        l.octanols,
        FOUR_TIMES,
        four_times_empty(),
        exact(START_TRANS),
        exact(DECANE.density),
        stages_text(&FOUR_TIMES_MELT),
        exact(every),
        l.at.len(),
        positions_digest(&l.at),
        model_digest(field)
    )
}

/// Where the finite-size run writes: `PANTOMETRY_DECANE_FINITE_SIZE_DIR`, which must be named, so
/// that it never writes anywhere it was not sent.
fn four_times_dir() -> std::path::PathBuf {
    let dir = std::path::PathBuf::from(std::env::var_os("PANTOMETRY_DECANE_FINITE_SIZE_DIR").expect(
        "PANTOMETRY_DECANE_FINITE_SIZE_DIR names the directory the finite-size run writes, and only it",
    ));
    std::fs::create_dir_all(&dir).expect("the finite-size directory");
    dir
}

/// **444 decanes, run and read by the registered rule**, under [`four_times_dir`]: read back from
/// `record.txt` and `end.txt` when a run wrote them under the same header ([`four_times_header`]),
/// either not whole or written differently refused by name; otherwise run now under a lock `run`,
/// the start first written to `built.txt`, [`FOUR_TIMES_MELT`] through [`integrate_with`] — `U` per
/// decane, `P`, `T`, the order and the trans fraction every 0.1 ps, logged to `run.log` in 2.5 ps
/// blocks — and the record and the end written atomically, the record first. Then
/// [`finite_size_verdict`] and the uncorrelated bands of 444 and 111, logged to `finite_size.log`.
/// Nothing about the model is asserted.
fn four_times_measured() {
    let t = std::time::Instant::now();
    let dir = four_times_dir();
    let a = Alkane::new(D10);
    let s = four_times_start(&a);
    let l = &s.liquid;
    assert_eq!(l.octanols, FOUR_TIMES_MOLECULES);
    let field = four_times_field(&a, l);
    let watch = four_times_watch(l, &a);
    let want = four_times_header(l, &field, RECORD_EVERY);
    let (record_path, end_path) = (dir.join("record.txt"), dir.join("end.txt"));
    let cell = l.cell();
    log(
        &dir,
        "finite_size",
        &format!(
            "{} decanes, {} atoms, box {:?} Å at {:.6} g/cm³; packed in {} iterations, {} kicks, worst \
             pair {:.5} of the radii; built in {:.1} s",
            l.octanols,
            l.at.len(),
            cell.lengths().map(|x| x / ANGSTROM),
            l.octanols as f64 * a.molecule.mass() / cell.volume() / 1e3,
            s.packed.iterations,
            s.packed.kicks,
            s.packed.worst,
            t.elapsed().as_secs_f64()
        ),
    );
    let rows = if let Ok(text) = std::fs::read_to_string(&record_path) {
        let rows = parse_record(&text, &want).unwrap_or_else(|why| {
            panic!(
                "{} is not this run's or not whole ({why}): remove it and end.txt to run again",
                record_path.display()
            )
        });
        let (head, end) = std::fs::read_to_string(&end_path)
            .map_err(|e| e.to_string())
            .and_then(|t| parse_positions(&t))
            .unwrap_or_else(|why| panic!("{} is not whole ({why})", end_path.display()));
        assert_eq!(head, want, "{} was run differently", end_path.display());
        assert_eq!(end.len(), l.at.len());
        log(&dir, "finite_size", "read back from record.txt and end.txt");
        rows
    } else {
        let _lock = Lock::take(&dir, "run");
        write_positions(&dir.join("built.txt"), &want, &l.at);
        let (end, rows) = integrate_with(
            "decane ×4",
            &l.dynamics(),
            l.octanols,
            |_: &MeltStage| Plain(&field),
            |at: &[[f64; 3]]| liquid_pressure(&field, l, at),
            &FOUR_TIMES_MELT,
            &l.at,
            FOUR_TIMES_SEED,
            &watch,
            RECORD_EVERY,
            &dir,
            "run",
        );
        let partial = record_path.with_extension("partial");
        std::fs::write(&partial, record_text(&want, &rows)).expect("the record");
        std::fs::rename(&partial, &record_path).expect("rename");
        write_positions(&end_path, &want, &end);
        rows
    };
    let (iso, sd) = isotropic_order(FOUR_TIMES_MOLECULES, 1000, 0x150);
    let (iso_once, sd_once) = isotropic_order(ONCE_MOLECULES, 1000, 0x150);
    log(
        &dir,
        "finite_size",
        &format!(
            "for information: an uncorrelated random set of {FOUR_TIMES_MOLECULES} directions gives \
             {iso:.4} ± {sd:.4} (band ≤ {:.4}), of {ONCE_MOLECULES} {iso_once:.4} ± {sd_once:.4} (band \
             ≤ {:.4}), a ratio of {:.3} against √4's 2",
            iso + 3.0 * sd,
            iso_once + 3.0 * sd_once,
            iso_once / iso
        ),
    );
    let (_, lines) = finite_size_verdict(&rows, RECORD_EVERY).unwrap_or_else(|why| panic!("{why}"));
    for line in &lines {
        log(&dir, "finite_size", line);
    }
    log(
        &dir,
        "finite_size",
        &format!("the whole in {:.0} s", t.elapsed().as_secs_f64()),
    );
}

/// **The order at four times the size** (release, one core, about 5.8 h; writes only
/// `PANTOMETRY_DECANE_FINITE_SIZE_DIR`): [`four_times_measured`]. 154 500 steps and 3 060 records
/// at the 132.1 ms a step and 114.4 ms a record that
/// [`the_decane_order_at_four_times_the_size_cost_measured`] read with throttling off and
/// AboveNormal priority; 238.0 and 258.3 ms, 10.4 h, with it on. Registered before it ran, with
/// [`finite_size_verdict`]'s rule and [`finite_size_gate`]. `PANTOMETRY_DECANE_FINITE_SIZE_DIR=
/// C:/pantometry-w4/decane_finite_size cargo test --release --locked -p pantometry-forcefield
/// --test benzene_partitioned_into_decane the_decane_order_at_four_times_the_size_measured --
/// --ignored --nocapture`.
#[test]
#[ignore = "444 decanes, 154 500 steps: about six hours with --release"]
fn the_decane_order_at_four_times_the_size_measured() {
    four_times_measured();
}

/// **The finite-size run's cost** (release, one core, minutes; writes nothing): the 444-decane box
/// built, 200 steps at 0.5 fs in a 50 ps⁻¹ bath off the packing, then 300 at 2 fs timed; a record's
/// pressure, order and trans fraction timed ten times; and [`FOUR_TIMES_MELT`] projected from them.
#[test]
#[ignore = "the 444-decane box built and stepped a few hundred times: minutes with --release"]
fn the_decane_order_at_four_times_the_size_cost_measured() {
    let t = std::time::Instant::now();
    let a = Alkane::new(D10);
    let s = four_times_start(&a);
    let l = &s.liquid;
    let field = four_times_field(&a, l);
    let watch = four_times_watch(l, &a);
    println!(
        "{} decanes, {} atoms, box {:.2?} Å; packed in {} iterations and {} kicks; built in {:.1} s",
        l.octanols,
        l.at.len(),
        l.cell().lengths().map(|x| x / ANGSTROM),
        s.packed.iterations,
        s.packed.kicks,
        t.elapsed().as_secs_f64()
    );
    let mut at = l.at.clone();
    let mut md = l
        .dynamics()
        .with_bath(Bath::Langevin {
            temperature: KELVIN,
            friction: 50e12,
            seed: 1,
        })
        .thermalised(&at, KELVIN, 1);
    md.run(&field, &mut at, 0.5 * FS, 200);
    let mut md = md.with_bath(Bath::Langevin {
        temperature: KELVIN,
        friction: 5e12,
        seed: 2,
    });
    let steps = 300;
    let s0 = std::time::Instant::now();
    md.run(&field, &mut at, 2.0 * FS, steps);
    let step = s0.elapsed().as_secs_f64() / steps as f64;
    let s1 = std::time::Instant::now();
    for _ in 0..10 {
        std::hint::black_box((
            liquid_pressure(&field, l, &at),
            watch.order(&at),
            watch.trans(&at),
        ));
    }
    let record = s1.elapsed().as_secs_f64() / 10.0;
    let (mut total_steps, mut records) = (0.0, 0.0);
    for st in &FOUR_TIMES_MELT {
        total_steps += (st.picoseconds * 1e-12 / st.dt).round();
        records += (st.picoseconds * 1e-12 / RECORD_EVERY).round();
    }
    // The strain-out stage's steps are taken at the 2 fs step's cost: the same force evaluation.
    let whole = total_steps * step + records * record;
    println!(
        "a step at 2 fs {:.1} ms; a record (P, order, trans) {:.1} ms; {total_steps} steps and \
         {records} records: {:.2} h",
        1e3 * step,
        1e3 * record,
        whole / 3600.0
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

/// **The extension's gate reads exactly its registered window**, records 1 000 to 3 000 of 3 000
/// at 0.1 ps. White noise about constant `U` and `P`, with the order rising over the first 1 000
/// from 0.07 to 0.15 — the first settle's shape — and flat at 0.15 after: **passes**, its line
/// naming 2 000 records over 100.1–300.0 ps. The order flat for 1 000 and then rising from 0.12 to
/// 0.18 is refused as rising. **The window's edges, exactly**: a `U` that is not a number at record
/// 999 passes, and at record 1 000 or 2 999 fails, since a half's mean that holds it is not a
/// number either. (A trend inside a half is no such probe: it lengthens that half's τ until its
/// error covers it, and the ramp above gated whole is 1.8 σ.) A record of 2 999 or 3 001 records,
/// a 200 ps one, or one whose first record is of another stage is refused as not the registered
/// one.
#[test]
fn the_extension_gate_reads_its_registered_window() {
    let band = 0.157;
    let n = 3000;
    let flat = [
        (0.1, 0.1, 0.0),
        (12.9, 0.0, 0.3),
        (60.0, 0.0, 300.0),
        (298.15, 0.0, 10.0),
        (0.15, 0.0, 0.01),
        (0.83, 0.0, 0.01),
    ];
    let rows = |order: &dyn Fn(usize) -> f64| -> Vec<MeltRow> {
        synthetic(n, 0xE47, flat)
            .into_iter()
            .enumerate()
            .map(|(i, (_, mut r))| {
                r[4] += order(i);
                (0, r)
            })
            .collect()
    };
    let early = rows(&|i| {
        if i < 1000 {
            -0.08 * (1000 - i) as f64 / 1000.0
        } else {
            0.0
        }
    });
    let ok = extension_verdict(&early, RECORD_EVERY, band)
        .unwrap_or_else(|l| panic!("a settled last 200 ps refused:\n{}", l.join("\n")));
    println!("{}", ok.join("\n"));
    assert!(
        ok[0].contains("over 2000 records, 100.1–300.0 ps"),
        "{}",
        ok[0]
    );
    let whole = match settle_verdict(&early, 0, 0, band) {
        Ok(l) | Err(l) => l,
    };
    println!("gated whole, for information: {}", whole[3]);
    for (record, passes) in [(999, true), (1000, false), (2999, false)] {
        let mut poisoned = rows(&|_| 0.0);
        poisoned[record].1[1] = f64::NAN;
        let v = extension_verdict(&poisoned, RECORD_EVERY, band);
        assert_eq!(v.is_ok(), passes, "a NaN at record {record}: {v:?}");
    }
    let late = rows(&|i| {
        if i < 1000 {
            -0.03
        } else {
            -0.03 + 0.06 * (i - 1000) as f64 / 2000.0
        }
    });
    let why = extension_verdict(&late, RECORD_EVERY, band).expect_err("rising late");
    assert!(why.join("\n").contains("the order is rising"), "{why:?}");
    println!("rising late: {}", why.last().unwrap());
    let flat_rows = rows(&|_| 0.0);
    assert!(extension_verdict(&flat_rows, RECORD_EVERY, band).is_ok());
    let mut other_stage = flat_rows.clone();
    other_stage[0].0 = 1;
    for (what, r) in [
        ("2 999", flat_rows[..2999].to_vec()),
        ("3 001", [&flat_rows[..], &flat_rows[..1]].concat()),
        ("200 ps", flat_rows[..2000].to_vec()),
        ("two stages", other_stage),
    ] {
        let why = extension_verdict(&r, RECORD_EVERY, band).expect_err(what);
        assert!(
            why.join(" ")
                .contains("not the record the gate was registered for"),
            "{what}: {why:?}"
        );
    }
}

/// **The 444-decane start is what it says, and deterministic**: [`four_times_start`] built twice.
///
/// - **444 decanes, 32 atoms each, in order**, the four [`four_times_empty`] sites one per
///   2 × 7 × 8 replica and distinct;
/// - **at [`DECANE`]'s density**: `444 m / V` within 16 ε, and every edge C-1b's box's replicated
///   2 × 2 × 1 times `(444/448)^(1/3)`, within 16 ε, the closed form for 448 sites shrunk to hold
///   444 at one density;
/// - **nothing overlaps, every molecule whole**: no two atoms nearer than their packing radii' sum by
///   brute force over the 27 images with each molecule taken whole ([`worst_contact`]); every bond
///   its relaxed conformer's to 64 ε and under half the shortest edge;
/// - **isotropic**: the end-to-end vectors' order inside the band an uncorrelated random set of 444
///   shows, its mean plus three standard deviations; and the trans fraction the generator's count
///   to the bit;
/// - **deterministic**: the second build is the first to the bit, every coordinate.
#[test]
fn the_four_times_start_is_isotropic_at_its_density_and_whole() {
    let a = Alkane::new(D10);
    let t = std::time::Instant::now();
    let s = four_times_start(&a);
    let built = t.elapsed().as_secs_f64();
    let l = &s.liquid;
    let o = &a.molecule;
    let empty = four_times_empty();
    let [_, nj, nk] = FOUR_TIMES;
    let replica = |e: usize| {
        let (i, j) = (e / (nj * nk), (e / nk) % nj);
        (i / DECANE.lattice[0], j / DECANE.lattice[1])
    };
    let mut replicas: Vec<(usize, usize)> = empty.iter().map(|&e| replica(e)).collect();
    replicas.sort_unstable();
    replicas.dedup();
    assert_eq!(
        replicas.len(),
        4,
        "one empty site in each replica: {empty:?}"
    );
    assert_eq!(
        FOUR_TIMES.iter().product::<usize>() - empty.len(),
        FOUR_TIMES_MOLECULES
    );
    assert_eq!(l.octanols, FOUR_TIMES_MOLECULES);
    assert_eq!(l.molecules.len(), FOUR_TIMES_MOLECULES);
    let mut next = 0;
    for r in &l.molecules {
        assert_eq!(r.start, next, "the molecules are the atoms in order");
        assert_eq!(r.len(), 32);
        next = r.end;
    }
    assert_eq!(next, l.at.len());

    let cell = l.cell();
    let density = l.octanols as f64 * o.mass() / cell.volume();
    let (once, _, _) = Liquid::cell_of(o, DECANE.lattice, DECANE.density);
    let shrink = (FOUR_TIMES_MOLECULES as f64 / 448.0).cbrt();
    let edges: Vec<f64> = (0..3)
        .map(|k| {
            let want = (FOUR_TIMES[k] / DECANE.lattice[k]) as f64 * once.lengths()[k] * shrink;
            cell.lengths()[k] / want - 1.0
        })
        .collect();
    println!(
        "built in {built:.1} s: {} iterations, {} kicks; {} decanes in {:.3?} Å at {:.10} g/cm³; the \
         edges C-1b's × (2, 2, 1) × (444/448)^(1/3) to {:.1e}",
        s.packed.iterations,
        s.packed.kicks,
        l.octanols,
        cell.lengths().map(|x| x / ANGSTROM),
        density / 1e3,
        edges.iter().fold(0.0f64, |m, x| m.max(x.abs()))
    );
    assert!(
        (density / DECANE.density - 1.0).abs() <= 16.0 * EPS,
        "{density}"
    );
    for (k, e) in edges.iter().enumerate() {
        assert!(e.abs() <= 16.0 * EPS, "edge {k}: {e:e}");
    }

    let (w, i, j) = worst_contact(&l.at, &l.molecules, &s.radii, &cell);
    println!(
        "the worst pair over 27 images, molecules whole, {w:.5} of its radii' sum, atoms {i} and {j}"
    );
    assert!(w >= 1.0, "an overlap at {w} of the radii");
    // A bond's two ends are each placed by a few roundings at the size of their coordinates — the
    // centre plus the arm, the conversion to metres, the shift into the box — so its length is its
    // conformer's to about `ε |x| / b` and not to `ε`: a 1.09 Å bond 70 Å from the origin carries
    // 60 times the relative rounding of one at the origin. Each bond is held to eight of that, and
    // the conformer's own rotation to 64 ε.
    let mut worst = 0.0f64;
    let mut longest = 0.0f64;
    for r in &l.molecules {
        for b in o.component.bonds() {
            let [x, y] = b.atoms;
            let (p, q) = (l.at[r.start + x], l.at[r.start + y]);
            let here = len(sub(p, q));
            let own = len(sub(o.at[x], o.at[y]));
            let reach = p.iter().chain(&q).fold(0.0f64, |m, c| m.max(c.abs()));
            let rounding = 64.0 * EPS + 8.0 * EPS * reach / own;
            worst = worst.max((here / own - 1.0).abs() / rounding);
            longest = longest.max(here);
        }
    }
    println!("every bond its conformer's to {worst:.3} of its rounding bound");
    assert!(
        worst <= 1.0,
        "a bond {worst} of its rounding bound off its conformer's"
    );
    assert!(longest < cell.half_shortest_edge());

    let ends = a.ends();
    let order = nematic_order(&directions(&l.at, &l.molecules, ends));
    let (iso, sd) = isotropic_order(FOUR_TIMES_MOLECULES, 1000, 0x150);
    let gauche: usize = s.gauche.iter().sum();
    let dihedrals = 7 * l.octanols;
    let (trans, _) = conformations(&l.at, &l.molecules, &a.carbons);
    println!(
        "order {order:.4} against an uncorrelated set of 444's {iso:.4} ± {sd:.4}, band ≤ {:.4}; \
         {gauche} gauche of {dihedrals}, trans {trans:.4}",
        iso + 3.0 * sd
    );
    assert!(order <= iso + 3.0 * sd, "the start's order {order}");
    assert_eq!(
        trans.to_bits(),
        ((dihedrals - gauche) as f64 / dihedrals as f64).to_bits()
    );

    let again = four_times_start(&a);
    let bits = |at: &[[f64; 3]]| -> Vec<u64> { at.iter().flatten().map(|x| x.to_bits()).collect() };
    assert_eq!(
        bits(&again.liquid.at),
        bits(&l.at),
        "the same seed, the same bits"
    );
}

/// Rows of the finite-size run for the rule's closed cases: [`FOUR_TIMES_MELT`]'s 10, 50 and
/// `settle` records at 0.1 ps, each column as [`synthetic`]'s.
fn four_times_rows(settle: usize, seed: u64, columns: [(f64, f64, f64); 6]) -> Vec<MeltRow> {
    synthetic(60 + settle, seed, columns)
        .into_iter()
        .enumerate()
        .map(|(i, (_, r))| (usize::from(i >= 10) + usize::from(i >= 60), r))
        .collect()
}

/// **The finite-size rule on closed cases**, each of its four outcomes and its edges:
///
/// - the thresholds are the registered ones: 0.137 and 0.158 to three places, from 1.3 × S₁₁₁/√4
///   and 0.75 × S₁₁₁;
/// - **the rule itself** ([`finite_size_outcome`]): at the finite-size threshold exactly,
///   finite-size, and one ulp above, inconclusive; at the ordered threshold exactly, genuine order,
///   and one ulp below, inconclusive; `U`'s or `P`'s halves 3 σ apart exactly settled, one ulp more
///   not settled, and not a number not settled; an order that is not a number, not settled;
/// - **on a record** ([`finite_size_verdict`]): white noise about an order of 0.08, finite-size; of
///   0.147, inconclusive; of 0.20, genuine order; `U` drifting 0.2 kcal/mol or `P` 200 bar over
///   the window, not settled even with the order at 0.08. **The window's edges**: a `U` that is not a
///   number at the settle's record 999 changes nothing, and at 1 000 or 2 999 makes it not settled.
///   A settle of 2 999 or 3 001 records is refused. The verdict line names the outcome and both
///   thresholds.
#[test]
fn the_finite_size_rule_on_closed_cases() {
    let (lo, hi) = finite_size_thresholds();
    println!("the thresholds: finite-size ≤ {lo:.6}, genuine order ≥ {hi:.6}");
    assert_eq!(
        (format!("{lo:.3}"), format!("{hi:.3}")),
        ("0.137".to_string(), "0.158".to_string()),
        "the registered thresholds"
    );
    let up = |x: f64| f64::from_bits(x.to_bits() + 1);
    let down = |x: f64| f64::from_bits(x.to_bits() - 1);
    for (order, u, p, want) in [
        (0.0, 0.0, 0.0, FiniteSize::Isotropic),
        (lo, 0.0, 0.0, FiniteSize::Isotropic),
        (up(lo), 0.0, 0.0, FiniteSize::Inconclusive),
        (down(hi), 0.0, 0.0, FiniteSize::Inconclusive),
        (hi, 0.0, 0.0, FiniteSize::Ordered),
        (0.5, 0.0, 0.0, FiniteSize::Ordered),
        (lo, HALVES_AGREE, HALVES_AGREE, FiniteSize::Isotropic),
        (lo, up(HALVES_AGREE), 0.0, FiniteSize::NotSettled),
        (hi, 0.0, up(HALVES_AGREE), FiniteSize::NotSettled),
        (0.147, f64::NAN, 0.0, FiniteSize::NotSettled),
        (0.147, 0.0, f64::NAN, FiniteSize::NotSettled),
        (f64::NAN, 0.0, 0.0, FiniteSize::NotSettled),
    ] {
        assert_eq!(
            finite_size_outcome(order, u, p),
            want,
            "order {order:e}, U {u:e} σ, P {p:e} σ"
        );
    }

    let settle = 3000;
    let flat = [
        (0.1, 0.1, 0.0),
        (12.5, 0.0, 0.3),
        (-45.0, 0.0, 300.0),
        (298.15, 0.0, 10.0),
        (0.08, 0.0, 0.01),
        (0.855, 0.0, 0.01),
    ];
    let judged = |cols: [(f64, f64, f64); 6]| {
        finite_size_verdict(&four_times_rows(settle, 0xF5, cols), RECORD_EVERY)
            .unwrap_or_else(|why| panic!("{why}"))
    };
    let (outcome, lines) = judged(flat);
    println!("{}", lines.join("\n"));
    assert_eq!(outcome, FiniteSize::Isotropic);
    assert!(
        lines[0].contains("2000 records, 106.1–306.0 ps"),
        "{}",
        lines[0]
    );
    let last = lines.last().unwrap();
    assert!(
        last.starts_with("the verdict: finite-size (isotropic liquid): ")
            && last.contains("finite-size ≤ 0.137 and genuine order ≥ 0.158"),
        "{last}"
    );
    for (level, want) in [
        (0.147, FiniteSize::Inconclusive),
        (0.20, FiniteSize::Ordered),
    ] {
        let mut cols = flat;
        cols[4].0 = level;
        let (outcome, lines) = judged(cols);
        assert_eq!(outcome, want, "{level}");
        assert!(lines.last().unwrap().contains(want.named()));
    }
    for (what, column, slope) in [("U", 1, 0.2 / 2000.0), ("P", 2, 200.0 / 2000.0)] {
        let mut cols = flat;
        cols[column].1 = slope;
        let (outcome, lines) = judged(cols);
        assert_eq!(outcome, FiniteSize::NotSettled, "{what} drifting");
        println!("{what} drifting: {}", lines.last().unwrap());
    }
    for (record, want) in [
        (999, FiniteSize::Isotropic),
        (1000, FiniteSize::NotSettled),
        (2999, FiniteSize::NotSettled),
    ] {
        let mut rows = four_times_rows(settle, 0xF5, flat);
        rows[60 + record].1[1] = f64::NAN;
        let (outcome, _) = finite_size_verdict(&rows, RECORD_EVERY).expect("judged");
        assert_eq!(outcome, want, "a NaN at the settle's record {record}");
    }
    for n in [2999, 3001] {
        let why = finite_size_verdict(&four_times_rows(n, 0xF5, flat), RECORD_EVERY)
            .expect_err("another length");
        assert!(
            why.contains("not the record the rule was registered for"),
            "{why}"
        );
    }
}

/// **The partition leg's finite-size gate on closed cases** ([`finite_size_gate`]), on records of
/// the extension's 3 000: white noise about stationary `U` and `P` and an order of 0.21 — the
/// extension's — **passes**, which the gate as first registered refuses as above the band; 0.29
/// passes; **0.31 is refused** as above the finite-size ceiling; and **an order rising** from 0.18
/// to 0.24 over the gated 200 ps, all of it under 0.30, is refused as rising. And
/// [`partition_gate`], the production's, passes 0.21 and refuses 0.31 as this gate does.
#[test]
fn the_finite_size_gate_for_the_partition_leg() {
    let band = 0.1569;
    let flat = [
        (0.1, 0.1, 0.0),
        (12.5, 0.0, 0.3),
        (-45.0, 0.0, 300.0),
        (298.15, 0.0, 10.0),
        (0.21, 0.0, 0.01),
        (0.855, 0.0, 0.01),
    ];
    let rows = |level: f64, rise: f64| -> Vec<MeltRow> {
        let mut cols = flat;
        cols[4].0 = level;
        synthetic(3000, 0x9A7E, cols)
            .into_iter()
            .enumerate()
            .map(|(i, (_, mut r))| {
                if i >= 1000 {
                    r[4] += rise * (i - 1000) as f64 / 2000.0;
                }
                (0, r)
            })
            .collect()
    };
    let ok = finite_size_gate(&rows(0.21, 0.0), RECORD_EVERY)
        .unwrap_or_else(|l| panic!("0.21 stationary refused:\n{}", l.join("\n")));
    println!("{}", ok.join("\n"));
    assert!(
        ok[0].contains("at most the finite-size ceiling 0.30"),
        "{}",
        ok[0]
    );
    let first = extension_verdict(&rows(0.21, 0.0), RECORD_EVERY, band).expect_err("the band");
    assert!(first.join("\n").contains("above the band"), "{first:?}");
    assert!(finite_size_gate(&rows(0.29, 0.0), RECORD_EVERY).is_ok());
    let high = finite_size_gate(&rows(0.31, 0.0), RECORD_EVERY).expect_err("0.31");
    assert!(
        high.join("\n")
            .contains("above the finite-size ceiling, 0.30"),
        "{high:?}"
    );
    let rising = finite_size_gate(&rows(0.18, 0.06), RECORD_EVERY).expect_err("rising");
    println!("rising: {}", rising.last().unwrap());
    assert!(
        rising.join("\n").contains("the order is rising"),
        "{rising:?}"
    );
    assert!(!rising.join("\n").contains("above the finite-size ceiling"));
    // The production's gate is this one: 0.21 stationary, which the band refuses, passes, and 0.31
    // does not.
    assert!(partition_gate(&rows(0.21, 0.0), RECORD_EVERY).is_ok());
    assert!(partition_gate(&rows(0.31, 0.0), RECORD_EVERY).is_err());
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

/// The resume test's extension: four steps at the settle's bath.
#[cfg(not(target_family = "wasm"))]
const TEST_EXTENSION: [MeltStage; 1] = [stage(
    "extend",
    Lambda::COUPLED,
    0.002,
    0.5 * FS,
    5e12,
    KELVIN,
)];

/// **A killed leg resumes to the same bits**, through this file's own [`decane_start`],
/// [`decane_window`] and `tests/partition/`'s [`Leg`] on the small box, a melt of seven steps and an
/// extension of four, each recording every step, and windows of eight 0.5 fs steps and three
/// samples at candidates 0 and 8, in temporary directories:
///
/// - **uninterrupted**: the melt, its extension and both windows. **The extension starts from
///   `start.txt` to the bit**: its first record's `U` and order are one step from there with
///   velocities from `window_seed(EXTENSION_SEED, 0)`, rebuilt here, its clock starts at zero, and
///   its header names `start.txt`'s digest. **The windows start from `extended.txt`**, not from
///   `start.txt`; the two windows' headers differ in the start's digest and nowhere else, and a
///   window run from `start.txt` is refused by the extension's leg;
/// - **killed during the extension**: the melt whole and a partial record left, then the leg again
///   — the melt's files unchanged, and every row and record the uninterrupted run's to the bit;
/// - **killed during a window**: the melt, the extension and window 0, then what a run killed
///   during window 8 leaves — its `.partial`, a header, one sample and its order — then the leg
///   again. Window 0 and the melt's and the extension's four files are read, their bytes unchanged; the partial is begun afresh and renamed; and
///   **every gradient, energy and order of both records, and every row of the melt's and the extension's records, is the
///   uninterrupted run's to the bit**.
///
/// And each guard refuses by name: the melt asked for with another seed; a window asked for at
/// another seed, from a start one ulp away, and **from the same start under another model** — the
/// decanes' charges QEq's — which only the header's model digest can tell; a window run without its
/// order; a `start.txt` cut short, and a `melt.txt` cut short; the extension at another seed, and
/// an `extension.txt` cut short.
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
    let begin = |dir: &std::path::Path| {
        decane_start(
            &sys,
            dir,
            DECANE_MELT_SEED,
            &TEST_EXTENSION,
            EXTENSION_SEED,
            &watch,
            every,
        )
    };
    let run = |dir: &std::path::Path, schedule: &[usize]| -> (Vec<Recorded>, Vec<MeltRow>) {
        let s = begin(dir);
        let leg = Leg::new(
            "decane",
            &sys.decoupling,
            &sys.template,
            s.start,
            protocol,
            grid(),
        );
        let records = schedule
            .iter()
            .map(|&g| decane_window(&leg, dir, g, &watch))
            .collect();
        (records, [s.melt, s.extension].concat())
    };
    let whole = fresh_dir("whole");
    let (full, full_rows) = run(&whole, &[0, 8]);
    assert_eq!(
        full_rows.len(),
        7 + 4,
        "a record every step, the melt's and the extension's"
    );
    let bytes = |p: std::path::PathBuf| std::fs::read(p).expect("a file");
    let positions = |p: std::path::PathBuf| {
        parse_positions(&std::fs::read_to_string(p).expect("a file"))
            .expect("whole")
            .1
    };
    let bits = |at: &[[f64; 3]]| -> Vec<u64> { at.iter().flatten().map(|x| x.to_bits()).collect() };

    // The extension starts from the melt's end, start.txt, to the bit, with the stated
    // re-thermalisation: its first record is one step from there, rebuilt here.
    let started = begin(&whole);
    let melt_end = positions(whole.join("start.txt"));
    assert_eq!(bits(&started.melted), bits(&melt_end));
    let ext = TEST_EXTENSION[0];
    let s0 = pantometry_forcefield::free_energy::window_seed(EXTENSION_SEED, 0);
    let mut md = sys
        .template
        .clone()
        .with_bath(Bath::Langevin {
            temperature: ext.temperature,
            friction: ext.friction,
            seed: s0,
        })
        .thermalised(&melt_end, ext.temperature, s0);
    let mut at = melt_end.clone();
    md.step(
        &AtLambda {
            hamiltonian: &sys.decoupling,
            lambda: ext.lambda,
        },
        &mut at,
        ext.dt,
    );
    let first = started.extension[0].1;
    let u = kcal(md.potential_energy().expect("stepped")) / sys.solvent as f64;
    assert_eq!(
        first[0].to_bits(),
        (every * 1e12).to_bits(),
        "the extension's clock from zero"
    );
    assert_eq!(first[1].to_bits(), u.to_bits(), "the extension's first U");
    assert_eq!(
        first[4].to_bits(),
        watch.order(&at).to_bits(),
        "its first order"
    );
    let ext_head = std::fs::read_to_string(whole.join("extension.txt")).unwrap();
    let ext_head = ext_head.lines().next().unwrap();
    assert!(
        ext_head.contains(&format!(
            "from {:016x} after # melted",
            positions_digest(&melt_end)
        )),
        "{ext_head}"
    );
    // The windows start from the extension's end, extended.txt, and their header says so.
    let ext_end = positions(whole.join("extended.txt"));
    assert_eq!(bits(&started.start), bits(&ext_end));
    assert_ne!(bits(&ext_end), bits(&melt_end));
    let leg_at = |start: &[[f64; 3]]| {
        Leg::new(
            "decane",
            &sys.decoupling,
            &sys.template,
            start.to_vec(),
            protocol,
            grid(),
        )
    };
    let (new, old) = (leg_at(&ext_end), leg_at(&melt_end));
    let (dn, d_old) = (positions_digest(&ext_end), positions_digest(&melt_end));
    assert_ne!(header(&new, 0), header(&old, 0));
    assert_eq!(
        header(&old, 0).replacen(
            &format!("start {d_old:016x}"),
            &format!("start {dn:016x}"),
            1
        ),
        header(&new, 0),
        "the start's digest is the whole difference"
    );
    let old_dir = fresh_dir("old-start");
    decane_window(&old, &old_dir, 0, &watch);
    let why = refusal(|| drop(decane_window(&new, &old_dir, 0, &watch)));
    assert!(
        why.contains("was run differently"),
        "a window from the melt's end: {why}"
    );
    println!(
        "the extension from start.txt {d_old:016x} to extended.txt {dn:016x}; its first U \
         {u:.6} kcal/mol rebuilt to the bit"
    );

    // Killed during the extension: the melt is whole, and the extension left a partial record.
    let mid = fresh_dir("mid-extension");
    melted_decane(&sys, &mid, DECANE_MELT_SEED, &watch, every);
    let melt_bytes = (bytes(mid.join("start.txt")), bytes(mid.join("melt.txt")));
    std::fs::write(mid.join("extension.partial"), "# extended, cut short\n").unwrap();
    let (resumed, resumed_rows) = run(&mid, &[0, 8]);
    assert!(
        (bytes(mid.join("start.txt")), bytes(mid.join("melt.txt"))) == melt_bytes,
        "the extension writes nothing the melt wrote"
    );
    let row_bits = |rows: &[MeltRow]| -> Vec<u64> {
        rows.iter()
            .flat_map(|r| r.1.map(f64::to_bits))
            .chain(rows.iter().map(|r| r.0 as u64))
            .collect()
    };
    assert_eq!(
        row_bits(&full_rows),
        row_bits(&resumed_rows),
        "killed in the extension"
    );
    for (a, b) in full.iter().zip(&resumed) {
        assert_eq!(
            record_bits(a),
            record_bits(b),
            "killed in the extension: {}",
            a.grid
        );
    }
    assert_eq!(full.len(), resumed.len());

    let killed = fresh_dir("killed");
    run(&killed, &[0]);
    let finished = bytes(killed.join("window_00.txt"));
    let start_bytes = bytes(killed.join("start.txt"));
    let record_bytes = bytes(killed.join("melt.txt"));
    let extended_bytes = bytes(killed.join("extended.txt"));
    let extension_bytes = bytes(killed.join("extension.txt"));
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
    assert!(bytes(killed.join("extended.txt")) == extended_bytes);
    assert!(bytes(killed.join("extension.txt")) == extension_bytes);
    assert!(!killed.join("window_08.partial").exists());
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
    // The extension asked for at another seed, and its record cut short.
    let why = refusal(|| {
        drop(decane_start(
            &sys,
            &killed,
            DECANE_MELT_SEED,
            &TEST_EXTENSION,
            EXTENSION_SEED + 1,
            &watch,
            every,
        ))
    });
    assert!(why.contains("extended differently"), "{why}");
    std::fs::write(cut.join("melt.txt"), &record_bytes).unwrap();
    std::fs::write(cut.join("extended.txt"), &extended_bytes).unwrap();
    std::fs::write(
        cut.join("extension.txt"),
        &extension_bytes[..extension_bytes.len() - 30],
    )
    .unwrap();
    let why = refusal(|| drop(begin(&cut)));
    assert!(why.contains("extension.txt is not whole"), "{why}");
    for d in [whole, killed, cut, plain, old_dir, mid] {
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
    let mut mine: Vec<u64> = vec![DECANE_SEED, DECANE_MELT_SEED, EXTENSION_SEED];
    mine.extend((0..grid().len()).map(|g| window_seed(DECANE_SEED, g)));
    mine.extend((0..DECANE_MELT.len()).map(|k| window_seed(DECANE_MELT_SEED, k)));
    mine.extend((0..DECANE_EXTENSION.len()).map(|k| window_seed(EXTENSION_SEED, k)));
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
